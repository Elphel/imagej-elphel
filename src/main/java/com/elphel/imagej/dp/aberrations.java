package com.elphel.imagej.dp;
/**
** -----------------------------------------------------------------------------**
** aberrations.java
**
** Programs for experimenting with aberrations correction
** 
** 
**
** Copyright (C) 2010 Elphel, Inc.
**
** -----------------------------------------------------------------------------**
**  
**  focus_tuning.java is free software: you can redistribute it and/or modify
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
import ij.io.*;
import ij.plugin.filter.GaussianBlur;
import ij.process.*;
import ij.gui.*;

import java.awt.*;
import java.awt.event.*;
import java.io.*;

import ij.plugin.frame.*;

import java.util.List;
import java.util.ArrayList;



import java.lang.Integer;

import javax.swing.*;

import com.elphel.imagej.common.DoubleFHT;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.jp4.JP46_Reader_camera;


public class aberrations extends PlugInFrame implements ActionListener {
  /**
	 * 
	 */
	private static final long serialVersionUID = -7799747215741901892L;
Panel panel1,panel2;
//  Panel panel;
  JP46_Reader_camera jp4_instance;
  ShowDoubleFloatArrays SDFA_instance;
  DoubleFHT fht_instance;
  static Frame instance;
 public static int DEBUG_LEVEL = 2;
 public static int MASTER_DEBUG_LEVEL = 2;
 public static boolean UPDATE_STATUS= true; // update ImageJ status info
 public static int NUMBER_OF_MEASUREMENTS = 100;
 public ImagePlus  imp_camera=null;
 public ImagePlus  imp_sel=null;
 public String title_src;
 public static int FFTSize=  256;
 public static int mapFFTSize=64; // used to find where grid covers the image
 public double GAUSS_WIDTH=0.4; //0 - use Hamming window
 public static int FFTOverlap=32;
 public static int INTERPOLATE_SUBDIV=   1; // interploate kernels between defined by FFTOverlap

 public static int subDivFreq=1; // increase FFT size to have higher resolution in frequency domain
 private static double [] Hamming;
 public static double [][] input_bayer; 
 public static double [][][] quarter_bayers; 

 public static int    simul_patternType=1;  // 0 - vanilla linear, 1 - curved,
 public static double simul_patternModifier=2.0;
 public static double simul_freqX1=0.0241; //0.028; //0.056; 
 public static double simul_freqY1=0.0055; //0.003; //0.006; 
 public static double simul_phase1=0.0; // radians
 public static double simul_freqX2=-0.0001; //-0.003; //-0.006; 
 public static double simul_freqY2=0.0231; //0.028; //0.056; 
 public static double simul_phase2=0.0; // radians
 public static double [] simul_corr={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; // Second order correction to the mesh (perspective+distortions)
 public static int    simul_subdiv=8; // 16; // subdivide pixels in each direction
 public static double simul_fill=0.5; // part of the (center) pixel area being "phptosensitive"
 public static boolean sim_hiRes=true;  
 public static double[][] simul_pixels=null;
 public static double corr_gamma=0.2; // pattern detection: gamma applied to PS
 public static double corr_sigma=1.5; // pattern detection: high-pass filter (0.0 - none) gamma(PS)
 public static double model_highpass=1.5; // model correlation high-pass filter (relative to pattern fundamental frequency - average of 2)
 public static double corrRingWidth=0.4;   // ring (around r=0.5 dist to opposite corr) width , center circle r=0.5*corrRingWidth
 public static double minCorrContrast=5.0; // Discrimination threshold between good and bad pattern correleation
 
 public static int    diffSpectrCorr=2; // maximal distance between maximum on spectrum and predicted maximum on autocorrelation of gamma(|spectrum|)
 public static double shrinkClusters=0.0; //  Shrink clusters by this ratio (remove lowest) after initial separation
 public static int    multiplesToTry=4; // try this number of m0.300aximums proportionally farther from 0,0 than the two closest (increase precision)

 public static double  deviation=1.0;     // when looking for maximums - maximal distance from predicted from the lower order one
 public static int     deviation_steps=6; // maximal iterations when looking for local maximum

 public static double  deconvInvert =  0.008; // 0.015; //0.03; // when FFT component is less than this fraction of the maximal value, replace 1/z with Z
 public static boolean centerForG2 = true;  // Align pattern to phases for the diagonal (both greens) sub-array
 public static boolean correctBayerRed=     true;
 public static boolean correctBayerBlue=    true;
 public static boolean correctBayerGreen1=  false;
 public static boolean correctBayerGreen2=  false;
 public static boolean correctBayerDiagonal=false;
 public static boolean correctBayerChecker= true;

 public static int     referenceComponent=  5; // component to calculate lateral chromatic from (0 - G1, 1 - R, 2 - B, 3 - G2,4 - diagonal greens, 5 - checker greens)
 public static boolean equalizeGreens=      true;   // equalize 2 greens in Bayer mosaic

 public static int     kernelHalSizeSingleColor=6; /// crop result convolution kernels to (kernelHalSizeSingleColor+1) * (kernelHalSizeSingleColor+1) for individual Bayer components
 public static int     kernelHalSizeDualColors= 9; /// crop result convolution kernels to (kernelHalSizeDualColors+1)  * (kernelHalSizeDualColors+1) for combined green Bayer components


 public static boolean forwardOTF=          true;

 public static double [][] kernels=         new double[6][];  // valid only for a small region of the frame

 public static int     PSF_subpixel=4;         // sub-pixel decimation 
// public static int     OTF_multiple=2;        // 0 - use each pixel once, 1 - add first negatives (4), 2 - second positives()4)
// public static double  PSF_minContrast=0.4;   // minimal instance contrast to use in binning (compared to the one at [0,0]
// public static boolean PSF_useHamming=true;   // multiply separated OTF instance by Hamming window
 public static double  PSF_minContrast=0.85;     // minimal instance contrast to use in binning (compared to the one at [0,0]
 public static double  PSF_windowFrac=0.5; //0.75;      // reduce the PSF cell size to this part of the area connecting first negative clones
 public static boolean PSF_useHamming=true;     // multiply separated OTF instance by Hamming window

 public static boolean PSF_symm180=false;      // make OTF center-symmetrical (around centroid that is defined by lateral chromatic aberration)

 public static int     OTF_FFT_size=32;        // initially will be increased by PSF_subpixel
 public static double  OTF_cutoff_energy=0.9; // 0.5;  // use frequency points that have OTF_cutoff_energy of the total to determine ellipse for limiting frequency responce
 public static double  OTF_ellipse_scale=0.9; // 2.5;  // size of elliptical window relative to the cluster defined by OTF_cutoff_energy
 public static boolean OTF_ellipse_gauss=true;  // size of elliptical window relative to the cluster defined by OTF_cutoff_energy

 public static double  OTF_deconvInvert =  0.007; // 0.015; // 0.01; // when FFT component is less than this fraction of the maximal value, replace 1/z with Z
/* TODO: check why combined greens result in center at x=-0.532/2 , y= -.0272/2 */
 public static boolean PSF_ignoreChromatic= false; // ignore lateral chromatic aberration (center OTF to 0,0)
 
 
 public static boolean OTF_fold =          false; // fold high frequency to lower when downsampling pixels (before inverse FFT)

 public static double  PSF_cutoff_energy=0.9;  // Limit result kernel to proportional of the PSF, calculate initial cluster shape by this cutoff energy
 public static double  PSF_ellipse_scale=0.5; // 0.3; //1.0;  // size of elliptical window to limuit reverse PSF as proportional to direct one
 public static double  RPSF_min_mask_threshold=0.01; // completely zero reversed kernel elements where elliptical mask is below this threshold

 public static double  RPSF_sigma_to_radius=0.2; // variable blurring - sigma will be proportional distance from the center
 public static double  RPSF_var_sigma_scale=0.8; // reduce variable sigma in the center from uniuform one

 public static double  OTF_zerofreq_size=  2.0; // used for filtering oversampling artifacts - size of zero freq maximum (if absent on simulated model PS)
 public static double  OTF_smoothPS=       2.5; // smooth model PS for rejecting aliases (0 - no smouth, >0 additional Gauss before FFT smaller than normal by this ratio)
 public static double  OTF_threshold_high= 0.02;// used for filtering oversampling artifacts - relative to max PS value to make filter completele rejecting
 public static double  OTF_threshold_low= 0.002;// used for filtering oversampling artifacts - relative to max PS to make filter completely transmissive

 

 public static double  XMASK_threshold=    0.01; // if 0 use normalize amplitude, if 0..1 - make binary: 1.0 if > threshold, 0.0 - otherwise
 public static double  XMASK_radius=        0.6; // low-pass result with low pass filter (fraction of average mesh period)
 public static int     XMASK_hsize=         128; // 2d histogram size (size/2 probably a good guess),
 public static double  XMASK_percentile=    0.1; // use this percentile (0.0..1.0)) value for given radius as a target
 public static double  XMASK_maxGain=       5.0; // maximal gain for low components
 public static double  XMASK_exaggerate=    1.0; // Exaggerate correction mask (Math,pow)

 public static double  PSF_smoothSeparate=  0.2; //125;  // low pass filter width when separation individual PSF
 public static double  PSF_thresholdSeparate=  0.2; // do not try to compensate for adjacent PSF clones if model sigma is less than this fraction of the one used for smoothing
 public static double  PSF_topCenter=        0.75; // consider only points above this fraction of the peak to find the centroid
 public static boolean PSF_removeNegtative=  false; // remove negative when separating composite PSF (will need low-pass filtering)
 public static double  PSF_sigmaToRadius=    0; //0.1; // variable-sigma blurring to reduce high frequencies more for the pixels farther from the PSF center
 public static double  PSF_wings_energy=     0.8; // fraction of energy in the pixels to be used
 public static double  PSF_wings_ellipse_scale=2.0; // increase wings cutoff ellipse by this from one defined by the  cutoff energy
 public static int     PSF_kernel_size=      32; // kernel (to be stored) size (per color component) 


// public static double  PSF_wings_min_mask_threshold=0.003; // zero output element if elliptical Gauss mask is below this threshold
// public static int     RPSF_kernel_size=     32; //  64; //32;  // size of deconvolution kernel
 public static int     RPSF_kernel_size=     64; //  64; //32;  // size of deconvolution kernel
 public static boolean SHOW_PSF=             false; // show combined PSF kernels (per-color and/or composite, as defined in SHOW_INDIVIDUAL, SHOW_COMPOSITE)
 public static boolean SHOW_MTF=             false; // calculate/show MTF (see notes to SHOW_PSF)
 public static boolean SHOW_INVERTED=        false; // show inverted kernels (unfiltered), same notes
 public static boolean SHOW_FILTERED=        true; // filter and show inverted kernels
 public static boolean SHOW_REDUCED_ALIASES= false; // calculate kernels with suppressed sampling aliases patterns
 public static boolean SHOW_GAUSSIANS=       true; // create gaussian kernels with the same centers as inverted ones (low noise, use in low details areas)
 public static boolean SHOW_REDUCED_ALIASES_GAUSSIAN= false; // calculate gaussian kernels with suppressed sampling aliases patterns
// public static boolean SHOW_INDIVIDUAL=      false; // for each of the kernels above - show per-color images
// public static boolean SHOW_COMPOSITE=       true; // for each of the kernels above - show single composite image
 
 
 public static double  DECONCV_BLUR_INDIVIDUAL=2.0;// 1.8;
 public static double  DECONCV_BLUR_DIAGONAL=  2.0;
 public static double  DECONCV_BLUR_CHECKER=   1.6;// 1.4;

 public static double  DECONCV_ALIASREJ_INDIVIDUAL=2.0;//
 public static double  DECONCV_ALIASREJ_DIAGONAL=  2.0;
 public static double  DECONCV_ALIASREJ_CHECKER=   2.0;//


 public static double  NONLIN_SIGMA=           5.0;  // sigma for the nonlinear filtering (higher the sigma, father from the edges extends the PSF correection)
 public static double  NONLIN_MIN=             0.01; //0.008; // 0.02; // minimal low-pass filtered squared difference between the corrected and original pixels to trigger sharpness enhancement
 public static double  NONLIN_MAX=             0.05;  //0.08; //0.3; // squared low-pass filtered difference between the corrected and original pixels, so abopve that level 100% corrected image is used
 public static double  NONLIN_THRESHOLD=       0.01; //when blurred intensity is below this value, use it as a denominator



 public static double  DECONCV_GAMMA=          0.5;
// public static double  DECONCV_SATURATION=     2.0;
 public static double  DECONCV_HEADROOM=       0.1;
 public static double  DECONCV_LOWPERC=        0.05;// 0.005;
 public static double  DECONCV_LOWLEVEL=       0.1;// 0.05;
 public static double  DECONCV_HIGHPERC=       0.05; //0.005;
 public static double  DECONCV_HIGHLEVEL=      0.9; //0.95;
 public static boolean PSF_enableModelSubtract=false; // not needed with Escher pattern, period >2* PSF width

 public static boolean MASK_BAYER_ALIASES= true;
 public static boolean CALC_BAYER_WEIGHTS= true;


 public static double BALANCE_RED =     1.1; // manual color balance, gain 1.0 matches 0.0.255.0 range of the unput Bayer data
 public static double BALANCE_BLUE =    1.4;
 public static double GAIN_GREEN =      1.0;
 public static double WEIGHT_SCALE_R =  1.0; // additional correction for the weights of colors for different sub-pixels in a Bayer cell
 public static double WEIGHT_SCALE_B =  1.0; // WEIGHT_SCALE_G=1.0
 public static double COLOR_SIGMA =     2.0; // Gaussian sigma to low-pass color components when calculating "smooth" color
 public static double YCbCr_Gamma=      0.50;
 public static double YCbCr_minLin=     0.003;
 public static double YCbCr_Kr=         0.299;
 public static double YCbCr_Kb=         0.114;
 public static double COLOR_SATURATION= 2.0;
 public static boolean USE_FIRST_Y=     true;

 public static double  DEBAYER_THRESHOLD=0.2; // Measuring maximal spectral component amplitude in mid-frequencies. If below this - use default de-bayer mask
 public static double  DEBAYER_WIDTH_GREEN=  1.5;
 public static double  DEBAYER_WIDTH_REDBLUE=  1.5;


 public static double  DEBAYER_GAMMA=    0.3;
 public static double  DEBAYER_RZ=       0.1; //0.15;  // later chnage them to space pixels
 public static double  DEBAYER_RA=       0.1 ; //0.15; // 0.25;   // fraction of the distance to diagonal (closest) alias
 public static double  DEBAYER_SIGMA=    0.5; //0.6; //0.50;  // relative to distance to nearest alias
 public static double  DEBAYER_DECAY=   -0.5 ; //0.5 //0.8 ; //0.4;  // relative to distance to nearest alias
 public static double  DEBAYER_RADIUS_POWER=   0.3; // Divide ray values by the radius to this power
 public static double  DEBAYER_MAINTOALIAS=  0.6; // 1.0; //0.7; // 0.6; //0.5; // relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)

 public static double  DEBAYER_FARTHEST_MAX=   0.5; // fartherst absolute maximum on a ray to count
 public static double  DEBAYER_MASK_BLUR=    2.0; // sigma for gaussian blur of the green and red/blue masks

 public static boolean DEBAYER_LO_GREEN=    true; // combine alias-reject "scissors" with lopass filter for greens
 public static boolean DEBAYER_LO_POSTGREEN=true;// combine alias-reject "scissors" with lopass filter for greens applied to red/blue
 public static boolean DEBAYER_LO_REDBLUE=  true;// combine alias-reject "scissors" with lopass filter for greens applied to red/blue


 public static boolean DEBAYER_TEST_FROMIMAGE=true;// false;
 public static boolean DEBAYER_TEST_MASKSPLIT=true;//  false;
 public static int     DEBAYER_FFT_SIZE=           64; //128;

 public static int     INTERPOLATE_INSIZE=       32;
 public static int     INTERPOLATE_STEP=          4;
 public static int     INTERPOLATE_ADDTOP=       17;
 public static int     INTERPOLATE_ADDLEFT=      17;
 public static int     INTERPOLATE_ADDRIGHT=     18;
 public static int     INTERPOLATE_ADDBOTTOM=    17;
 public static double  INTERPOLATE_EXTRAPOLATE= 0.0;

 public static int     INVERSE_DIRECT_SIZE=    32; // size (side of square) of direct PSF kernel
 public static int     INVERSE_REVERSE_SIZE=   64; // size (side of square) of reverse PSF kernel
 public static boolean INVERSE_FILTER=       true; // apply variable-sigma filtering to the inverted PSF


 public static int     SPLIT_OVERSAMPLE=    2;
 public static int     SPLIT_ADDTOP=       32;
 public static int     SPLIT_ADDLEFT=      32;
 public static int     SPLIT_ADDRIGHT=     32;
 public static int     SPLIT_ADDBOTTOM=    32;


 public static int     CONVOLVE_FFT_SIZE=  128;



 public static double [][] averagedKernels= new double[6][];  // Averaged over rectangular area (FFTSize/2 overlap

 public static boolean [][] simulation_barray=null;  // global over-sampled bit image array (2-d, square) used to extract particular Bayer components

 public static int        patternSize=512; // size of the side of the square pattern bitmap 
 public static boolean [] bPattern=null; // pattern bitmap
 public static boolean [] colorsToCorrect=new boolean[6];

 public static double [][] deconvKernels={null,null,null,null,null,null,null};
 public static double []   encodedDirectKernels=null;
 public static double []   encodedInvertKernels=null;
 public static String [] componentColorNames={"green1","red","blue","green2", "greens (diagonal)", "greens (checker)"};
 public static String [] stackColorNames={"red","green","blue"};

       static File dir;
 public static String [] filenames;
 public static double [][][][][] filePatternMaps;
 public static boolean [][][]    fileTilesMaps; // map of 2*FFTSize x 2*FFTSize squares with 2*FFTOverlap step, true means that that tile can be used for PSF
 

 public static double [][][][] PSFKernelMap=null;
 public static double [][][][] rPSFKernelMap=null;
 public static double [][][][] gaussianKernelMap=null;

 public static ImagePlus imp_convolved=null;
 public static ImagePlus imp_gaussian=null;

 public static int CHANNEL_RED=    1;
 public static int CHANNEL_BLUE=   2;
 public static int CHANNEL_GREEN=  5;
 public static int CHANNEL_CHECKER=5;

 public static double DOUBLE_DEBUG_RESULT; // just for temporary passing results from inside the functions

 public static ImageStack convolutionKernelStack=null; // select to use for image convolution
  public aberrations() {
    super("Aberrations");
    if (IJ.versionLessThan("1.39t")) return;
    if (instance!=null) {
      instance.toFront();
      return;
    }
    instance = this;
    addKeyListener(IJ.getInstance());
    setLayout(new FlowLayout());
    panel1 = new Panel();
    panel1.setLayout(new GridLayout(4, 1, 50, 5));

    addButton("Configure",panel1);
//    addButton("Average");
//    addButton("Select");
    addButton("Find pattern",panel1);
    addButton("Map image",panel1);
    addButton("Process files",panel1);
    add(panel1);
    panel2 = new Panel();
    panel2.setLayout(new GridLayout(4, 3, 5, 5));
    addButton("Find distorted",panel2);
    addButton("Simulate",panel2);
    addButton("Simulate Quarters",panel2);
    addButton("Split Bayer",panel2);
    addButton("Get PSF",panel2);
    addButton("Test Convolve",panel2);
    addButton("Map PSF",panel2);
    addButton("Process Kernels",panel2);
    addButton("Interpolate Kernels",panel2);
    addButton("Inverse Stack",panel2);
    addButton("Gaussian Stack",panel2);
    addButton("Split Image",panel2);
    addButton("Debayer Image",panel2);
    addButton("Select kernel stack",panel2);
    addButton("Convolve with stack",panel2);
    addButton("Read Deconvs",panel2);
    addButton("Convolve Image",panel2);
    addButton("Read Gaussian",panel2);
    addButton("Combine pair",panel2);
    addButton("Test",panel2);
    addButton("Colors",panel2);
    addButton("Test Debayer",panel2);
//    addButton("Test Debayer0");
    add(panel2);
    pack();
    GUI.center(this);
    setVisible(true);
    jp4_instance=       new JP46_Reader_camera();
    fht_instance=       new DoubleFHT();
    SDFA_instance=      new ShowDoubleFloatArrays();


    Hamming=initHamming(FFTSize);

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
    int i,j,l,iq;
    double [][] pattern;
    double [] patternCorr=new double[6]; // second order non-linear pattern correction (perspective+distortion)
    double [][][] quarter_patterns;
    double [] quarterHamming;
    double [][] inverted=new double[6][];
    String label = e.getActionCommand();
//    boolean [] colorsToCorrect={correctBayerGreen1,correctBayerRed,correctBayerBlue,correctBayerGreen2,correctBayerDiagonal,correctBayerChecker};
    colorsToCorrect[0]=correctBayerGreen1;
    colorsToCorrect[1]=correctBayerRed;
    colorsToCorrect[2]=correctBayerBlue;
    colorsToCorrect[3]=correctBayerGreen2;
    colorsToCorrect[4]=correctBayerDiagonal;
    colorsToCorrect[5]=correctBayerChecker;
    for (referenceComponent=5;(referenceComponent>=0) && (!colorsToCorrect[referenceComponent]); referenceComponent--);

    for (i=0;i<inverted.length;i++) inverted[i]=null;
    for (i=0;i<inverted.length;i++) kernels[i]=null;

    if (label==null) return;
/* ======================================================================== */
    if (label.equals("Configure")) {
      if (showConfigureDialog()) {
      Hamming=initHamming(FFTSize);

      }
      return;
    }
    DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
/*
    if (label.equals("Average")) {
      if (showAverageDialog()) {
      }
      return;

    } else  */
/* ======================================================================== */
    if (label.equals("Select")) {
      imp_sel= selectImage(2*FFTSize, false);
      return;

/* ======================================================================== */
    } else if  (label.equals("Split Bayer")) {
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
  
//    imp_sel= selectImage(2*FFTSize, false);
      imp_sel= selectImage(0, true);
/* second pass - no windowes */
      Roi roi_sel= imp_sel.getRoi();
      if (roi_sel==null){
        imp_sel.setRoi(0, 0, imp_sel.getWidth(), imp_sel.getHeight());
        roi_sel= imp_sel.getRoi();
      }
      Rectangle rroi=roi_sel.getBounds();
      input_bayer=splitBayer (imp_sel,rroi,equalizeGreens);
      if (!colorsToCorrect[4] && (input_bayer.length>4)) input_bayer[4]=null; // do not use memory fro diagonal greens here
      int sel_width= rroi.width/2;
      int sel_height=rroi.height/2;

//      for (i=0;i<4;i++) if (!colorsToCorrect[i]) input_bayer[i]=null; // leave composite greens even if disabled
//      input_bayer= normalizeAndWindow (input_bayer, Hamming);
      if (subDivFreq>1) {
         input_bayer= extendFFTInput (input_bayer,sel_width, subDivFreq);
         sel_width*= subDivFreq;
         sel_height*=subDivFreq;
      }
      if ((PSF_subpixel>1) && (OTF_zerofreq_size>=0)){
        input_bayer= oversampleFFTInput (input_bayer,sel_width,PSF_subpixel);
        sel_width*=PSF_subpixel;
        sel_height*=PSF_subpixel;
        if (colorsToCorrect[5]) input_bayer=combineCheckerGreens (input_bayer,   // pixel arrays after oversampleFFTInput() or extractSimulPatterns())
                                                                  sel_width, 
                                                                  PSF_subpixel); // same as used in oversampleFFTInput() - oversampling ratio


      }
      for (i=0;i<4;i++) if (!colorsToCorrect[i]) input_bayer[i]=null; // leave composite greens even if disabled
  

      if (DEBUG_LEVEL>2) System.out.println ("Colors to process ("+input_bayer.length+")"+
                                             " " +componentColorNames[0]+"-"+colorsToCorrect[0]+
                                             ", "+componentColorNames[1]+"-"+colorsToCorrect[1]+
                                             ", "+componentColorNames[2]+"-"+colorsToCorrect[2]+
                                             ", "+componentColorNames[3]+"-"+colorsToCorrect[3]+
                                             ", "+componentColorNames[4]+"-"+colorsToCorrect[4]+
                                             ", "+componentColorNames[5]+"-"+colorsToCorrect[5]);
      if (DEBUG_LEVEL>2) System.out.println ("sel_width="+sel_width+" sel_height="+sel_height);
      if (colorsToCorrect[1] && colorsToCorrect[2] && colorsToCorrect[5]) {
        ImageStack bayerStack= combineRGBCorrectBayer ( input_bayer[1],
                                                        input_bayer[5],
                                                        input_bayer[2],
                                                             sel_width,   // image width
                                                        PSF_subpixel/2);  // half Bayer period (GR/BG)
        ImagePlus imp_bayerStack = new ImagePlus(imp_sel.getTitle()+"bayer-stack", bayerStack);
        imp_bayerStack.getProcessor().resetMinAndMax();
        imp_bayerStack.show();
        if (DEBUG_LEVEL>2) SDFA_instance.showArrays(input_bayer, sel_width,sel_height, imp_sel.getTitle());
      } else {
        SDFA_instance.showArrays(input_bayer, sel_width,sel_height, imp_sel.getTitle());
//      SDFA_instance.showArrays(input_bayer, imp_sel.getTitle());
      }
      return;
/* ======================================================================== */
    } else if  (label.equals("Get PSF")) {


      double averagePeriod=1.0/Math.sqrt((simul_freqX1*simul_freqX1+simul_freqY1*simul_freqY1+simul_freqX2*simul_freqX2+simul_freqY2*simul_freqY2)/2.0);

      if (!showPSFDialog()) return;
      if (!showRPSFDialog()) return; // second dialog for filtering the inverted kernel

      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;

      imp_sel= selectImage(2*FFTSize, false);
      input_bayer=splitBayer (imp_sel,equalizeGreens);

/* adding functionality of the other commands */

      double[][] distortedPattern= findPatternDistorted(   input_bayer, // pixel array to process (no widowing!)
                                                            corr_gamma,
                                                            corr_sigma,  // pattern detection: high-pass filter (0.0 - none) gamma(PS)
                                                        diffSpectrCorr,
                                                        shrinkClusters,
                                                        multiplesToTry,
                                                             deviation,
                                                       deviation_steps,
                                                                  true, // boolean greens, // this is a pattern for combined greens (diagonal), adjust results accordingly
                                                   imp_sel.getTitle()); // title prefix to use for debug  images
      if (distortedPattern==null) return;

      int patternSize=512; /* TODO - calculate? */
      boolean [] bPattern= patternGenerator(patternSize,simul_patternType, simul_patternModifier);



       simulation_barray=  simulatePatternFullPattern(bPattern,
                                              distortedPattern[0][0],
                                              distortedPattern[0][1],
                                              distortedPattern[0][2],
                                              distortedPattern[1][0],
                                              distortedPattern[1][1],
                                              distortedPattern[1][2],
                                              distortedPattern[2], //
                                              simul_subdiv,
                                              FFTSize,
                                              centerForG2);


//      for (i=0;i<4;i++) if (!colorsToCorrect[i]) input_bayer[i]=null; // leave composite greens even if disabled

      input_bayer= normalizeAndWindow (input_bayer, Hamming);
      if (subDivFreq>1)  input_bayer= extendFFTInput (input_bayer,subDivFreq);
      if ((PSF_subpixel>1) && (OTF_zerofreq_size>=0)) {
        input_bayer= oversampleFFTInput (input_bayer,PSF_subpixel);
        if (colorsToCorrect[5])  input_bayer=combineCheckerGreens (input_bayer,   // pixel arrays after oversampleFFTInput() or extractSimulPatterns())
                                                                   PSF_subpixel); // same as used in oversampleFFTInput() - oversampling ratio
      }
      for (i=0;i<4;i++) if (!colorsToCorrect[i]) input_bayer[i]=null; // leave composite greens even if disabled

      if (DEBUG_LEVEL>3) {
        if (OTF_zerofreq_size>=0) SDFA_instance.showArrays(input_bayer, FFTSize*subDivFreq*PSF_subpixel, FFTSize*subDivFreq*PSF_subpixel, imp_sel.getTitle());
        else                      SDFA_instance.showArrays(input_bayer, FFTSize*subDivFreq, FFTSize*subDivFreq, imp_sel.getTitle());
      }

      boolean goodSimul= (simulation_barray != null) && (simulation_barray.length>=simul_subdiv*2*FFTSize);
 
      if (goodSimul) {
        simul_pixels= extractSimulPatterns (simulation_barray, // high resolution boolean pattern array
                                                   simul_fill, // part of the (center) pixel area being "phptosensitive"
                                                 simul_subdiv, // boolean pixels to real pixels resolution
                                                 PSF_subpixel, // subdivide pixels
                                         FFTSize*PSF_subpixel, // number of Bayer cells in width of the square selection (half number of pixels)
                                                          0.0,    // selection center, X (in pixels)
                                                          0.0);   // selection center, y (in pixels)
        double [] fullHamming=initHamming( FFTSize*PSF_subpixel);
        simul_pixels= normalizeAndWindow (simul_pixels, fullHamming);
        if (subDivFreq>1) {
               simul_pixels= extendFFTInput (simul_pixels,subDivFreq);

        }
        if ((PSF_subpixel>1) && (OTF_zerofreq_size>=0)) {
                      if (colorsToCorrect[5])  simul_pixels=combineCheckerGreens (simul_pixels,   // pixel arrays after oversampleFFTInput() or extractSimulPatterns())
                                                                          PSF_subpixel); // same as used in oversampleFFTInput() - oversampling ratio
        }

//        if (DEBUG_LEVEL>3) SDFA_instance.showArrays(simul_pixels, FFTSize*subDivFreq*PSF_subpixel, FFTSize*subDivFreq*PSF_subpixel, imp_sel.getTitle()+"-SIM");
        if (DEBUG_LEVEL>2) SDFA_instance.showArrays(simul_pixels, FFTSize*subDivFreq*PSF_subpixel, FFTSize*subDivFreq*PSF_subpixel, imp_sel.getTitle()+"-SIM");
        if (DEBUG_LEVEL>2) System.out.println ( " input_bayer.length="+input_bayer.length+" simul_pixels.length="+simul_pixels.length);

        for (i=0;(i<input_bayer.length) && (i<simul_pixels.length);i++) if ((colorsToCorrect[i]) && (input_bayer[i]!=null)){
          if (DEBUG_LEVEL>2) System.out.println ( "Color "+componentColorNames[i]+" XMAX_radius ("+XMASK_radius+") is re-calculated into bayer pixels as "+ ((i==4)?Math.sqrt(2):1.0)*0.5*averagePeriod*XMASK_radius);
          if (OTF_zerofreq_size<0) {

            inverted[i]=limitedInverseOfFHTDiffSize(input_bayer[i],
                                         simul_pixels[i],
                                            deconvInvert,
                                              forwardOTF,
                                imp_sel.getTitle()+"-"+i);
          } else {
            inverted[i]=limitedInverseOfFHT(input_bayer[i],
                                         simul_pixels[i],
                         FFTSize*subDivFreq*PSF_subpixel,
                                                  (i==5),     //    boolean checker // checkerboard pattern in the source file (use when filtering)
                                            deconvInvert,
                                              forwardOTF,
//new below
                                            PSF_subpixel,
                                       OTF_zerofreq_size,
                                            OTF_smoothPS,
                                      OTF_threshold_high,
                                       OTF_threshold_low,
                                         XMASK_threshold,
((i==4)?Math.sqrt(2):1.0)*0.5*averagePeriod*XMASK_radius,
                                             XMASK_hsize,
                                        XMASK_percentile,
                                           XMASK_maxGain,
                                        XMASK_exaggerate,
                                imp_sel.getTitle()+"-"+i);

          }
        }
        if (DEBUG_LEVEL>1) SDFA_instance.showArrays(inverted, FFTSize*subDivFreq*PSF_subpixel, FFTSize*subDivFreq*PSF_subpixel, imp_sel.getTitle()+"_Combined-PSF");
      }

/* correct composite greens */
/* Here we divide wave vectors by PSF_subpixel as the pixels are already added */
      double [][]wVectors=  {{2.0*distortedPattern[0][0]/PSF_subpixel,                     2.0*distortedPattern[0][1]/PSF_subpixel}, // pattern was in sensor pixels, not bayer sub-pixel
                             {2.0*distortedPattern[1][0]/PSF_subpixel,                     2.0*distortedPattern[1][1]/PSF_subpixel}};


      double [][] wVrotMatrix= {{0.5,0.5},{-0.5,0.5}};

      double [][]wVectors4= new double [2][2];
      for (i=0;i<2;i++) for (j=0;j<2;j++) {
         wVectors4[i][j]=0.0;
         for (l=0;l<2;l++) wVectors4[i][j]+=wVectors[i][l]*wVrotMatrix[l][j];
      }



      double [][] PSF_shifts=         new double [input_bayer.length][];    // X/Y shift of the PSF array, in Bayer component pioxel coordinates (same as PSF arrays)
      double [][] PSF_centroids=      new double [input_bayer.length][];    // X/Y coordinates of the centroids of PSF in Bayer component pioxel coordinates (same as PSF arrays) (after they were optionally shifted)
      double [][] lateralChromatic=   new double [input_bayer.length][]; // X/Y coordinates of the centroids of Bayer component PSF in sensor pixel coordinates
      double [][]kernelsForFFT=       new double [input_bayer.length][];
      double [][] psf_inverted=       new double [input_bayer.length][];
      double [][] psf_inverted_masked=new double [input_bayer.length][];
//      double [][] restoredInput=      new double [input_bayer.length][];

      double [] lateralChromaticAbs=new double [input_bayer.length];
      double [] zeroVector={0.0,0.0};
      for (i=input_bayer.length-1;i>=0;i--) {
        if (colorsToCorrect[i]) {
          PSF_shifts[i]=       zeroVector.clone();
          PSF_centroids[i]=    zeroVector.clone();
          lateralChromatic[i]= zeroVector.clone();
        } else {
          PSF_shifts[i]=       null;
          PSF_centroids[i]=    null;
          lateralChromatic[i]= null;
        }
        lateralChromaticAbs[i]=0.0;
        kernelsForFFT[i]=null;
        psf_inverted[i]=null;
        psf_inverted_masked[i]=null;
 //       restoredInput[i]=null;
      }
//      int [][]  clusterMask;
//      int referenceComponent=4; /*  use 5 late, made global r */
/* Start with referenceComponent */
      i= referenceComponent;
      if (DEBUG_LEVEL>2) {
                 System.out.println("Before: color Component "+i+" PSF_shifts["+i+"][0]="+IJ.d2s(PSF_shifts[i][0],3)+
                                                " PSF_shifts["+i+"][1]="+IJ.d2s(PSF_shifts[i][1],3));
      }  

      kernels[i]=combinePSF (inverted[i], // Square array of pixels with multiple repeated PSF (alternating sign)
               (i==4)?wVectors4:wVectors, // two wave vectors, lengths in cycles/pixel (pixels match pixel array)
                                       1, // already applied PSF_subpixel,         // sub-pixel decimation 
                         PSF_minContrast,   // minimal instance contrast to use in binning
                          PSF_windowFrac,    // reduce the PSF cell size to this part of the area connecting first negative clones
                          PSF_useHamming,
                             PSF_symm180, //PSF_symm180,
                                    true, //PSF_ignoreChromatic
                           PSF_shifts[i],  // centerXY[] - will be modified inside combinePSF() if PSF_ignoreChromatic is true
                         PSF_centroids[i], // will return array of XY coordinates of the result centroid
                  PSF_enableModelSubtract,  // generate/subtract gaussian models (not needed if no overlap between pos/neg clones)
                       PSF_smoothSeparate,  // low-pass filter multiple opposite-sign PSF instaces for separation (width relative to distance to the opposite sign PSF)
                    PSF_thresholdSeparate,  // threshold for locating zero-crossing
                            PSF_topCenter,  // consider only points above this fraction of the peak to find the centroid
                      PSF_removeNegtative,  // remove PSF negative values  when separating composite PSF (will need low-pass filtering)
                        PSF_sigmaToRadius,   // 0.4;  variable-sigma blurring to reduce high frequencies more for the pixels farther from the PSF center
                          PSF_wings_energy, // fraction of energy in the pixels to be used
                   PSF_wings_ellipse_scale,
//              PSF_wings_min_mask_threshold, // zero output element if elliptical Gauss mask is below this threshold
                 imp_sel.getTitle()+"_"+i,
                          (DEBUG_LEVEL>4));
      if (DEBUG_LEVEL>2)     System.out.println("After-1: color Component "+i+"    PSF_shifts["+i+"][0]="+IJ.d2s(PSF_shifts   [i][0],3)+"    PSF_shifts["+i+"][1]="+IJ.d2s(   PSF_shifts[i][1],3));
      if (DEBUG_LEVEL>2)     System.out.println("After-1: color Component "+i+" PSF_centroids["+i+"][0]="+IJ.d2s(PSF_centroids[i][0],3)+" PSF_centroids["+i+"][1]="+IJ.d2s(PSF_centroids[i][1],3));

      if ((referenceComponent==4) && !PSF_ignoreChromatic) { /* Recalculate center to pixels from greens (diagonal)) and supply it to other colors (lateral chromatic aberration correction) */
        for (j=0;j<input_bayer.length;j++) if ((colorsToCorrect[j]) && (j!=referenceComponent)) {
          PSF_shifts[j]=shiftSensorToBayer (shiftBayerToSensor(PSF_shifts[i],4,PSF_subpixel),j,PSF_subpixel);
          if (DEBUG_LEVEL>2)       System.out.println("After-2 (recalc): color Component "+j+" PSF_shifts["+j+"][0]="+IJ.d2s(PSF_shifts[j][0],3)+" PSF_shifts["+j+"][1]="+IJ.d2s(PSF_shifts[j][1],3));
        }
      }
      lateralChromatic[i]=shiftBayerToSensor ( PSF_shifts[i][0]+PSF_centroids[i][0],
                                               PSF_shifts[i][1]+PSF_centroids[i][1],
                                               i,
                                               PSF_subpixel);
      lateralChromaticAbs[i]=Math.sqrt(lateralChromatic[i][0]*lateralChromatic[i][0]+lateralChromatic[i][1]*lateralChromatic[i][1]);

/* Now process all the other components */

      for (i=0; i<input_bayer.length;i++) if ((i!=referenceComponent) && (colorsToCorrect[i])) {
        if (DEBUG_LEVEL>2) {
                   System.out.println("Before: color Component "+i+" PSF_shifts["+i+"][0]="+IJ.d2s(PSF_shifts[i][0],3)+
                                                  " PSF_shifts["+i+"][1]="+IJ.d2s(PSF_shifts[i][1],3));
        }  
        kernels[i]=combinePSF (inverted[i], // Square array of pixels with multiple repeated PSF (alternating sign)
                 (i==4)?wVectors4:wVectors, // two wave vectors, lengths in cycles/pixel (pixels match pixel array)
                                         1, // already applied PSF_subpixel,         // sub-pixel decimation 
                           PSF_minContrast,   // minimal instance contrast to use in binning
                            PSF_windowFrac,    // reduce the PSF cell size to this part of the area connecting first negative clones
                            PSF_useHamming,
                               PSF_symm180, //PSF_symm180,
                       PSF_ignoreChromatic, //PSF_ignoreChromatic
                             PSF_shifts[i],  // centerXY[] - will be modified inside combinePSF() if PSF_ignoreChromatic is true
                          PSF_centroids[i], // will return array of XY coordinates of the result centroid
                   PSF_enableModelSubtract,  // generate/subtract gaussian models (not needed if no overlap between pos/neg clones)
                        PSF_smoothSeparate,  // low-pass filter multiple opposite-sign PSF instaces for separation (width relative to distance to the opposite sign PSF)
                     PSF_thresholdSeparate,  // threshold for locating zero-crossing
                             PSF_topCenter,  // consider only points above this fraction of the peak to find the centroid
                       PSF_removeNegtative,  // remove PSF negative values  when separating composite PSF (will need low-pass filtering)
                         PSF_sigmaToRadius,   // 0.4;  variable-sigma blurring to reduce high frequencies more for the pixels farther from the PSF center
                          PSF_wings_energy, // fraction of energy in the pixels to be used
                   PSF_wings_ellipse_scale,
//              PSF_wings_min_mask_threshold, // zero output element if elliptical Gauss mask is below this threshold
                  imp_sel.getTitle()+"_"+i,
                           (DEBUG_LEVEL>4));
        if (DEBUG_LEVEL>2)     System.out.println("After-1: color Component "+i+"    PSF_shifts["+i+"][0]="+IJ.d2s(PSF_shifts   [i][0],3)+"    PSF_shifts["+i+"][1]="+IJ.d2s(   PSF_shifts[i][1],3));
        if (DEBUG_LEVEL>2)     System.out.println("After-1: color Component "+i+" PSF_centroids["+i+"][0]="+IJ.d2s(PSF_centroids[i][0],3)+" PSF_centroids["+i+"][1]="+IJ.d2s(PSF_centroids[i][1],3));
        lateralChromatic[i]=shiftBayerToSensor ( PSF_shifts[i][0]+PSF_centroids[i][0],
                                                 PSF_shifts[i][1]+PSF_centroids[i][1],
                                                                                    i,
                                                                          PSF_subpixel);
        lateralChromaticAbs[i]=Math.sqrt((lateralChromatic[i][0]-lateralChromatic[referenceComponent][0])*(lateralChromatic[i][0]-lateralChromatic[referenceComponent][0])+
                                         (lateralChromatic[i][1]-lateralChromatic[referenceComponent][1])*(lateralChromatic[i][1]-lateralChromatic[referenceComponent][1]));
      }
      if (DEBUG_LEVEL>1) {
        for (i=0;i<PSF_shifts.length;i++) if (colorsToCorrect[i]){
          if (DEBUG_LEVEL>2) {
            System.out.println("Color Component "+i+" PSF_subpixel="+PSF_subpixel+
                                                    " PSF_ignoreChromatic="+PSF_ignoreChromatic+
                                                    " PSF_symm180="+PSF_symm180);
            System.out.println(                     " PSF_shifts["+i+"][0]="+IJ.d2s(PSF_shifts[i][0],3)+
                                                    " PSF_shifts["+i+"][1]="+IJ.d2s(PSF_shifts[i][1],3)+
                                                    " PSF_centroids["+i+"][0]="+IJ.d2s(PSF_centroids[i][0],3)+
                                                    " PSF_centroids["+i+"][1]="+IJ.d2s(PSF_centroids[i][1],3));
            System.out.println("  lateralChromatic["+i+"][0]="+IJ.d2s(lateralChromatic[i][0],3)+
                               "  lateralChromatic["+i+"][1]="+IJ.d2s(lateralChromatic[i][1],3));
          }
        }
        if (colorsToCorrect[referenceComponent]) for (i=0;i<colorsToCorrect.length;i++) if ((colorsToCorrect[i])&& (i!=referenceComponent)){

         System.out.println(componentColorNames[i]+" lateral chromatic (from green) "+IJ.d2s(lateralChromaticAbs[i],3)+"pix:  ["+i+"][0]="+IJ.d2s(lateralChromatic[i][0]-lateralChromatic[referenceComponent][0],3)+
                            "  ["+i+"][1]="+IJ.d2s(lateralChromatic[i][1]-lateralChromatic[referenceComponent][1],3));
        }
        System.out.println("Lateral shift green from simulation "+IJ.d2s(lateralChromaticAbs[referenceComponent],3)+"pix:  ["+referenceComponent+"][0]="+IJ.d2s(lateralChromatic[referenceComponent][0],3)+
                                                                                                                    "  ["+referenceComponent+"][1]="+IJ.d2s(lateralChromatic[referenceComponent][1],3));
      }

/* testing encode/decodeKernels */
/*
      encodedDirectKernels=encodeKernels (kernels, // arrays of color component kernels (or nulls), each kernel as a 1-d array
                               2*RPSF_kernel_size,  // size of the side of the square of the result
                                            true);//    boolean normalize);
      if (DEBUG_LEVEL>1) {
         SDFA_instance.showArrays(encodedDirectKernels, 2*RPSF_kernel_size, 2*RPSF_kernel_size, imp_sel.getTitle()+"encoded");
     
      }
      int greensToExtract= colorsToCorrect[5]?2:(colorsToCorrect[4]?1:0);
      double [][] decodedDirectKernels=decodeKernels (encodedDirectKernels, // combined kernels
                                                        2*RPSF_kernel_size, // twice bigger - just for testing // int ksize, // output size (side of square) of the extracted individual kernels
                                                        greensToExtract); // extract combined greens (component 4), when false - components 0 and 3

      if (DEBUG_LEVEL>1) {
         SDFA_instance.showArrays(decodedDirectKernels, 2*RPSF_kernel_size, 2*RPSF_kernel_size, imp_sel.getTitle()+"decoded");
     
      }
*/
      for (i=0;i<input_bayer.length;i++) if (colorsToCorrect[i]){
        kernelsForFFT[i]=resizeForFFT(kernels[i],OTF_FFT_size*PSF_subpixel);
      }
      if (DEBUG_LEVEL>1) {
         SDFA_instance.showArrays(kernelsForFFT,OTF_FFT_size*PSF_subpixel, OTF_FFT_size*PSF_subpixel, imp_sel.getTitle()+"K-"+deconvInvert);
        for (i=0;i<kernelsForFFT.length;i++) if (kernelsForFFT[i]!=null){
          System.out.println("Direct kernel "+i+" shot noise factor= "+kernelShotNoiseFactor (kernelsForFFT[i]));
        }
      }
//-------------------------

//      double [][] psf_inverted= {null,null,null,null,null};
      for (i=0;i<input_bayer.length;i++) if (colorsToCorrect[i]){
        psf_inverted[i]=cleanupAndReversePSF (kernelsForFFT[i],  // input pixels
                                              OTF_deconvInvert,
                                             OTF_cutoff_energy,
                                             OTF_ellipse_scale,
                                             OTF_ellipse_gauss,
                                                             1, //          PSF_subpixel,
                                                      OTF_fold,
                                                        "E-"+i);          // just for the plot names
      }
/* Find direct kernel approximation ellipse, increase it, mirror center around 0,0 and use it as a mask for the reversed kernel */
//      double [][] psf_inverted_masked= {null,null,null,null,null};
      for (i=0;i<input_bayer.length;i++) if (colorsToCorrect[i]){
        psf_inverted_masked[i]=maskReversePSFKernel
                                          (kernelsForFFT[i], // direct PSF function, square array, may be proportionally larger than reversed
                                            psf_inverted[i], // reversed psf, square array
                                          PSF_cutoff_energy, // fraction of energy in the pixels to be used
                                          PSF_ellipse_scale,
                                    RPSF_min_mask_threshold, // zero output element if elliptical Gauss mask is below this threshold
                                                     "ME-"+i);
      }
      for (i=0;i<input_bayer.length;i++){
        if (colorsToCorrect[i]){
          deconvKernels[i]=resizeForFFT(psf_inverted_masked[i],RPSF_kernel_size);
        } else {
          deconvKernels[i]=null;
        }
      }
      if (DEBUG_LEVEL>2) {
        SDFA_instance.showArrays(deconvKernels,RPSF_kernel_size, RPSF_kernel_size, "IK-"+deconvInvert+"-"+PSF_subpixel);
      }
      if (DEBUG_LEVEL>1) {
        for (i=0;i<deconvKernels.length;i++) if (deconvKernels[i]!=null){
          System.out.println("Inverted kernel "+i+" shot noise factor= "+kernelShotNoiseFactor (deconvKernels[i]));
        }
      }

/* testing encode/decodeKernels */
      encodedInvertKernels=encodeKernels (deconvKernels, // arrays of color component kernels (or nulls), each kernel as a 1-d array
                                     2*RPSF_kernel_size,  // size of the side of the square of the result
                                                   true);//    boolean normalize);
      if (DEBUG_LEVEL>2) {
         SDFA_instance.showArrays(encodedInvertKernels, 2*RPSF_kernel_size, 2*RPSF_kernel_size, imp_sel.getTitle()+"encoded");
     
      }
      int greensToExtract= colorsToCorrect[5]?2:(colorsToCorrect[4]?1:0);
      int decodedInvertKernelsSize=((greensToExtract==1)?2:1)*RPSF_kernel_size;
      double [][] decodedInvertKernels=decodeKernels (encodedInvertKernels, // combined kernels
                                                  decodedInvertKernelsSize, // twice bigger if diagonal greens // int ksize, // output size (side of square) of the extracted individual kernels
                                                           greensToExtract); // extract combined greens (component 4), when false - components 0 and 3

      if (DEBUG_LEVEL>1) {
         SDFA_instance.showArrays(decodedInvertKernels, decodedInvertKernelsSize, decodedInvertKernelsSize, imp_sel.getTitle()+"decoded");
      }
/* Add filtering of the iverted kernels */
      if (RPSF_sigma_to_radius<=0) return;
      double sigma=DECONCV_BLUR_INDIVIDUAL;
      double [] smoothKernel;
      double [] centroidXY;
      double [] variableSigmas;
      int size;
      for (i=0;i<decodedInvertKernels.length;i++) if (decodedInvertKernels[i]!=null) {
        switch (i) {
          case 0:
          case 1:
          case 2:
          case 3:
             sigma=DECONCV_BLUR_INDIVIDUAL;
             break;
          case 4:
              sigma=DECONCV_BLUR_DIAGONAL;
             break;
          case 5:
             sigma=DECONCV_BLUR_CHECKER;
             break;
        }
        smoothKernel=lowPassGauss(kernelsForFFT[i], 2*sigma,true);
        size= (int) Math.sqrt(decodedInvertKernels[i].length);
        centroidXY=extractCentroidFromReverseKernel(smoothKernel,   // square array of direct PSF
                                                             0.5); // fraction of the maximal value to use as a bottom of the part used for centroid calculation
        variableSigmas= createSigmasFromCenter(size, // side of square
                               RPSF_sigma_to_radius, // variable blurring - sigma will be proportional distance from the center
                         sigma*RPSF_var_sigma_scale, //blurring in the center sigma(r)=sqrt((sigma_to_radius*r)^2)+center_sigma^2)
                                       -centroidXY[0], // coordinates of the center (0:0 - size/2: size/2)
                                       -centroidXY[1]);
        decodedInvertKernels[i]=variableGaussBlurr(decodedInvertKernels[i], // input square pixel array, preferrably having many exact zeros (they will be skipped)
                                                            variableSigmas, // array of sigmas to be used for each pixel, matches pixels[]
                                                                       3.5, // drop calculatin if farther then nSigma
                                                                         0, // int WOICenterX, // window of interest in pixels[] array - do not generate data outside it
                                                                         0, // int WOICenterY, // 
                                                                      size, //int WOIWidth, reduce later
                                                                      size); //int WOIHeight)
/* normalize deconvolve kernels */
        normalizeKernel(decodedInvertKernels[i]); // in-place
      }
      if (DEBUG_LEVEL>1) {
        SDFA_instance.showArrays(decodedInvertKernels, decodedInvertKernelsSize, decodedInvertKernelsSize, imp_sel.getTitle()+"filtered");
        for (i=0;i<decodedInvertKernels.length;i++) if (decodedInvertKernels[i]!=null){
          System.out.println("Filtered inverted kernel "+i+" shot noise factor= "+kernelShotNoiseFactor (decodedInvertKernels[i]));
        }
      }
      return;
/* ======================================================================== */
    } else if (label.equals("Test Convolve")) {
      if (!showDeconvDialog()) return;

      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      if (encodedInvertKernels==null) {
         IJ.showMessage("Error","No inverted kernels available, please run 'Get PSF' first");
         return;
      }
      imp_sel= selectImage(2*FFTSize, false);
      input_bayer=splitBayer (imp_sel,equalizeGreens);
//      input_bayer= normalizeAndWindow (input_bayer, Hamming);
      if (subDivFreq>1)  input_bayer= extendFFTInput (input_bayer,subDivFreq);
      if ((PSF_subpixel>1) && (OTF_zerofreq_size>=0)) {
        input_bayer= oversampleFFTInput (input_bayer,PSF_subpixel);
        if (colorsToCorrect[5])  input_bayer=combineCheckerGreens (input_bayer,   // pixel arrays after oversampleFFTInput() or extractSimulPatterns())
                                                                   PSF_subpixel); // same as used in oversampleFFTInput() - oversampling ratio
      }

      int greensToExtract= colorsToCorrect[5]?2:(colorsToCorrect[4]?1:0);
      int decodedInvertKernelsSize=((greensToExtract==1)?2:1)*RPSF_kernel_size;
      double [][] decodedInvertKernels=decodeKernels (encodedInvertKernels, // combined kernels
                                                  decodedInvertKernelsSize, // twice bigger if diagonal greens // int ksize, // output size (side of square) of the extracted individual kernels
                                                           greensToExtract); // extract combined greens (component 4), when false - components 0 and 3
      if (DEBUG_LEVEL>2) {
//         SDFA_instance.showArrays(encodedInvertKernels, 2*RPSF_kernel_size, 2*RPSF_kernel_size, imp_sel.getTitle()+"encoded");
         SDFA_instance.showArrays(encodedInvertKernels, 2*RPSF_kernel_size, 2*RPSF_kernel_size, "encoded-kernels");
     
      }
      if (DEBUG_LEVEL>2) {
         SDFA_instance.showArrays(decodedInvertKernels, decodedInvertKernelsSize, decodedInvertKernelsSize, imp_sel.getTitle()+"decoded");
      }




      double [][] restoredInput=            new double [input_bayer.length][];
      double [][] smoothDeconvKernels=      new double [input_bayer.length][];
      double [][] smoothInputBayers=        new double [input_bayer.length][];
      double [][] gaussianKernels=          new double [input_bayer.length][];
      double [][] gaussianInput=            new double [input_bayer.length][];
      double [][] aliasRejectMask=          new double [input_bayer.length][];
      int      [] aliasMaskSize=            new int    [input_bayer.length];
      boolean  [] isChecker=                new boolean [input_bayer.length];
      double [][] smoothDeconvKernelsFHT=   new double [input_bayer.length][];
      double [][] gaussianKernelsFHT=       new double [input_bayer.length][];
      double [][] variableSigmas=           new double [input_bayer.length][];

      double      gaussianBlurr= DECONCV_BLUR_INDIVIDUAL;
      double      aliasReject=   DECONCV_ALIASREJ_INDIVIDUAL;
      double []   centroidXY;
      for (i=input_bayer.length-1;i>=0;i--) {
        restoredInput[i]=      null;
        smoothDeconvKernels[i]=null;
        gaussianKernels[i]=    null;
        gaussianInput[i]=      null;
        aliasRejectMask[i]=    null;
        aliasMaskSize[i]=         0;
        smoothDeconvKernelsFHT[i]=null;
        gaussianKernelsFHT[i]=    null;
        variableSigmas[i]=    null;
        if (deconvKernels[i]!=null) {
          switch (i) {
            case 0:
            case 1:
            case 2:
            case 3:
               gaussianBlurr=DECONCV_BLUR_INDIVIDUAL;
               aliasReject=  DECONCV_ALIASREJ_INDIVIDUAL;
               isChecker[i]=false;
               break;
            case 4:
               gaussianBlurr=DECONCV_BLUR_DIAGONAL;
               aliasReject=  DECONCV_ALIASREJ_DIAGONAL;
               isChecker[i]=false;
               break;
            case 5:
               gaussianBlurr=DECONCV_BLUR_CHECKER;
               aliasReject=  DECONCV_ALIASREJ_CHECKER;
               isChecker[i]=true;
               break;
          }
          smoothDeconvKernels[i]=lowPassGauss(decodedInvertKernels[i], 2*gaussianBlurr,true);
          smoothInputBayers[i]=  lowPassGauss(input_bayer[i],   2*gaussianBlurr,true);
          gaussianKernels[i]=extractLateralChromaticFromReverseKernel(smoothDeconvKernels[i],gaussianBlurr,0.5);


          if (RPSF_sigma_to_radius>0.0) { // recalculate smoothDeconvKernels[i]
            centroidXY=extractCentroidFromReverseKernel(smoothDeconvKernels[i],   // square array of reversed PSF
                                                                           0.5); // fraction of the maximal value to use as a bottom of the part used for centroid calculation

            variableSigmas[i]= createSigmasFromCenter((int) Math.sqrt(smoothDeconvKernels[i].length), // side of square
                                                                                RPSF_sigma_to_radius, // variable blurring - sigma will be proportional distance from the center
                                                                  RPSF_var_sigma_scale*gaussianBlurr, //blurring in the center sigma(r)=sqrt((sigma_to_radius*r)^2)+center_sigma^2)
                                                                                       centroidXY[0], // coordinates of the center (0:0 - size/2: size/2)
                                                                                       centroidXY[1]);
            if (DEBUG_LEVEL>1) {
              SDFA_instance.showArrays(smoothDeconvKernels[i], decodedInvertKernelsSize, decodedInvertKernelsSize, imp_sel.getTitle()+"pre-smooth"+i);
              SDFA_instance.showArrays(variableSigmas[i],      decodedInvertKernelsSize, decodedInvertKernelsSize, imp_sel.getTitle()+"sigmas"+i);
            }

            smoothDeconvKernels[i]=variableGaussBlurr(decodedInvertKernels[i], // input square pixel array, preferrably having many exact zeros (they will be skipped)
                                                            variableSigmas[i], // array of sigmas to be used for each pixel, matches pixels[]
                                                                            4, // drop calculatin if farther then nSigma
                                                                            0, // int WOICenterX, // window of interest in pixels[] array - do not generate data outside it
                                                                            0, // int WOICenterY, // 
                                                     decodedInvertKernelsSize, //int WOIWidth, reduce later
                                                     decodedInvertKernelsSize); //int WOIHeight)
          }

          if (DEBUG_LEVEL>1) {
            SDFA_instance.showArrays(smoothDeconvKernels[i], decodedInvertKernelsSize, decodedInvertKernelsSize, imp_sel.getTitle()+"smooth"+i);
          }
/* normalize deconvolved kernels? */
          normalizeKernel(smoothDeconvKernels[i]); // in-place
          if (DEBUG_LEVEL>1) {
            SDFA_instance.showArrays(smoothDeconvKernels[i], decodedInvertKernelsSize, decodedInvertKernelsSize, imp_sel.getTitle()+"normalized"+i);
          }




/* This test does not reject aliases in originals, just in kernels */

          if (aliasReject>0.0) {
             aliasMaskSize[i]= (int) Math.sqrt(smoothDeconvKernels[i].length);
             aliasRejectMask[i]= createAliasReject (aliasMaskSize[i],  // size of the mask
                                                        isChecker[i],  // checkerboard pattern in the source file (use when filtering)
                                                        PSF_subpixel,  // measured array is sampled at 1/oversample frequency than model (will add more parameters later)
                                                         aliasReject); // width of rejection areas on the spectrum (the smaller, the sharper rejection)
             smoothDeconvKernels[i]= rejectByMask (smoothDeconvKernels[i],  // square input data
                                                       aliasRejectMask[i], // mask to multiply FHT
                                                                    true); // image is centered around the center of the square (use swapQuadrants)
             gaussianKernels[i]=     rejectByMask (    gaussianKernels[i],  // square input data
                                                       aliasRejectMask[i], // mask to multiply FHT
                                                                    true); // image is sentered around the center of the square (use swapQuadrants)
          }
        }
      }
/* Now test the results by convolving with the measured mesh (with sub-pixels) */
//    smoothRejectedClonesPixels= lowPassGauss(rejectedClonesPixels, smoothSigma, true); /* not yet used */
      if (DEBUG_LEVEL>2) {
         if (aliasReject>0.0) SDFA_instance.showArrays(aliasRejectMask, decodedInvertKernelsSize, decodedInvertKernelsSize, imp_sel.getTitle()+"alias mask");
         SDFA_instance.showArrays(smoothDeconvKernels, decodedInvertKernelsSize, decodedInvertKernelsSize, imp_sel.getTitle()+"smooth-decoded");
         SDFA_instance.showArrays(gaussianKernels,     decodedInvertKernelsSize, decodedInvertKernelsSize, imp_sel.getTitle()+"gaussian");
      }


      double [] slidingMask=null;
      int inputSize;

      for (i=0;i<input_bayer.length;i++)  if (colorsToCorrect[i]){
        inputSize= (int) Math.sqrt(input_bayer[i].length); 
        smoothDeconvKernelsFHT[i]=extendFFTInput (smoothDeconvKernels[i], 2);
        fht_instance.swapQuadrants(smoothDeconvKernelsFHT[i]);
        fht_instance.transform(smoothDeconvKernelsFHT[i]);
        gaussianKernelsFHT[i]=    extendFFTInput (gaussianKernels[i], 2);
        fht_instance.swapQuadrants(gaussianKernelsFHT[i]);
        fht_instance.transform(gaussianKernelsFHT[i]);
        if (slidingMask==null) slidingMask=getSlidingMask(RPSF_kernel_size); // all colors - the same size
        restoredInput[i]= convolveByFHT(input_bayer[i], // rectangular image array to be convolved
                                             inputSize, // image width,
                             smoothDeconvKernelsFHT[i], // FHT transform of the DOUBLE-SIZE convolution kernel)
                                          slidingMask); // sliding mask (made of 0.25(cos(X)+1)*(cos(y)+1), or null (will be calculated)

          gaussianInput[i]=convolveByFHT(input_bayer[i], // rectangular image array to be convolved
                                              inputSize, // image width,
                                  gaussianKernelsFHT[i], // FHT transform of the DOUBLE-SIZE convolution kernel)
                                            slidingMask); // sliding mask (made of 0.25(cos(X)+1)*(cos(y)+1), or null (will be calculated)
      }



      int maskColor;
      for (maskColor=restoredInput.length-1;maskColor<=0;maskColor--) if ((restoredInput[maskColor]!=null) && (gaussianInput[maskColor]!=null)) break;
      double [] nonlinMask=null;
      if (maskColor>=0) {
        nonlinMask=createNonlinMask(restoredInput[maskColor],  // image with PSF deconvolution applied
                                    gaussianInput[maskColor],  // image convolved with shiftted gaussian (reject subsampling modulation, compensate lateral chromatic aberration)
                                                NONLIN_SIGMA,  // gaussian sigma to create mask for selecting between restored/noisy and original/smooth
                                                  NONLIN_MIN,  // minimal value of the low-pass filtered squared difference where mask >0.0
                                                  NONLIN_MAX, // maximal value of the low-pass filtered squared difference where mask <1.0                                              );
                                            NONLIN_THRESHOLD);
      }


      if (DEBUG_LEVEL>2) {
//        for (i=0;i<input_bayer.length;i++) if (!colorsToCorrect[i]) input_bayer[i]=null;
        SDFA_instance.showArrays(smoothInputBayers, FFTSize*subDivFreq*PSF_subpixel, FFTSize*subDivFreq*PSF_subpixel, imp_sel.getTitle()+"-orig");
        SDFA_instance.showArrays(gaussianInput,     FFTSize*subDivFreq*PSF_subpixel, FFTSize*subDivFreq*PSF_subpixel, imp_sel.getTitle()+"-gauss");
      }
      if (DEBUG_LEVEL>1) {
//        for (i=0;i<input_bayer.length;i++) if (!colorsToCorrect[i]) input_bayer[i]=null;
        SDFA_instance.showArrays(restoredInput,     FFTSize*subDivFreq*PSF_subpixel, FFTSize*subDivFreq*PSF_subpixel, imp_sel.getTitle()+"-proc");
      }


     double [][] combinedMaskedImages= combineMaskImages (restoredInput,
                                                          gaussianInput,
                                                             nonlinMask);
     if (DEBUG_LEVEL>1) {


/* converts array of Bayer planes (including combined greens) into array of RGB triplets */
  double [][] doubleColorProc=    mergeBayer(restoredInput, // array af color components arrays, some may be null
                           FFTSize*subDivFreq*PSF_subpixel, // component width (not the result image)
                                              PSF_subpixel);

/* simple white and brightness/contrast balance */
  double [][] balColorProc= balanceDoubleRGB (doubleColorProc, // pixel array of triplets
                                              DECONCV_LOWPERC, // fraction of pixels to become min (lowLevel). Use negative value to preserve original zeros
                                             DECONCV_LOWLEVEL, // level, so lowPerc of all pixels will be lower than that
                                             DECONCV_HIGHPERC, // fraction of pixels to become max
                                            DECONCV_HIGHLEVEL);// level, so highPerc of all pixels will be higher than that
/* Gamma correction for color images */
  double [][] gamColorProc= gammaDoubleRGB (balColorProc, // pixel array of triplets
                                           DECONCV_GAMMA,
                                        DECONCV_HEADROOM);


 showDoubleColor(gamColorProc, // array of color triplets
                                              null, //ImagePlus imp, // ImagePlus to reuse or null for the new one 
                   FFTSize*subDivFreq*PSF_subpixel,
                    FFTSize*subDivFreq*PSF_subpixel,
                       imp_sel.getTitle()+"-color");

   }

   if (DEBUG_LEVEL>1) {

/* converts array of Bayer planes (including combined greens) into array of RGB triplets */
  double [][] doubleColorGaussian=    mergeBayer(gaussianInput, // array af color components arrays, some may be null
                           FFTSize*subDivFreq*PSF_subpixel, // component width (not the result image)
                                              PSF_subpixel);

/* simple white and brightness/contrast balance */
  double [][] balColorGaussian= balanceDoubleRGB (doubleColorGaussian, // pixel array of triplets
                                              DECONCV_LOWPERC, // fraction of pixels to become min (lowLevel). Use negative value to preserve original zeros
                                             DECONCV_LOWLEVEL, // level, so lowPerc of all pixels will be lower than that
                                             DECONCV_HIGHPERC, // fraction of pixels to become max
                                            DECONCV_HIGHLEVEL);// level, so highPerc of all pixels will be higher than that
/* Gamma correction for color images */
  double [][] gamColorGaussian= gammaDoubleRGB (balColorGaussian, // pixel array of triplets
                                           DECONCV_GAMMA,
                                        DECONCV_HEADROOM);


 showDoubleColor(gamColorGaussian, // array of color triplets
                                              null, //ImagePlus imp, // ImagePlus to reuse or null for the new one 
                   FFTSize*subDivFreq*PSF_subpixel,
                    FFTSize*subDivFreq*PSF_subpixel,
                       imp_sel.getTitle()+"-color-gauss");
  }
  if (DEBUG_LEVEL>1) {

  double [][] doubleColorOrig=    mergeBayer(smoothInputBayers, // array af color components arrays, some may be null
                               FFTSize*subDivFreq*PSF_subpixel, // component width (not the result image)
                                                 PSF_subpixel);
/* simple white and brightness/contrast balance */
  double [][] balColorOrig= balanceDoubleRGB (doubleColorOrig, // pixel array of triplets
                                              DECONCV_LOWPERC, // fraction of pixels to become min (lowLevel). Use negative value to preserve original zeros
                                             DECONCV_LOWLEVEL, // level, so lowPerc of all pixels will be lower than that
                                             DECONCV_HIGHPERC, // fraction of pixels to become max
                                            DECONCV_HIGHLEVEL);// level, so highPerc of all pixels will be higher than that
/* Gamma correction for color images */
  double [][] gamColorOrig= gammaDoubleRGB (balColorOrig, // pixel array of triplets
                                           DECONCV_GAMMA,
                                        DECONCV_HEADROOM);


 showDoubleColor(gamColorOrig, // array of color triplets
                                                   null, //ImagePlus imp, // ImagePlus to reuse or null for the new one 
                        FFTSize*subDivFreq*PSF_subpixel,
                        FFTSize*subDivFreq*PSF_subpixel,
                       imp_sel.getTitle()+"-Orig-color");

  }
  if (DEBUG_LEVEL>0) {

/* converts array of Bayer planes (including combined greens) into array of RGB triplets */
  double [][] doubleColorComb=    mergeBayer(combinedMaskedImages, // array af color components arrays, some may be null
                                  FFTSize*subDivFreq*PSF_subpixel, // component width (not the result image)
                                                    PSF_subpixel);

/* simple white and brightness/contrast balance */
  double [][] balColorComb= balanceDoubleRGB (doubleColorComb, // pixel array of triplets
                                              DECONCV_LOWPERC, // fraction of pixels to become min (lowLevel). Use negative value to preserve original zeros
                                             DECONCV_LOWLEVEL, // level, so lowPerc of all pixels will be lower than that
                                             DECONCV_HIGHPERC, // fraction of pixels to become max
                                            DECONCV_HIGHLEVEL);// level, so highPerc of all pixels will be higher than that
/* Gamma correction for color images */
  double [][] gamColorComb= gammaDoubleRGB (balColorComb, // pixel array of triplets
                                           DECONCV_GAMMA,
                                        DECONCV_HEADROOM);


 showDoubleColor(gamColorComb, // array of color triplets
                                                  null, //ImagePlus imp, // ImagePlus to reuse or null for the new one 
                       FFTSize*subDivFreq*PSF_subpixel,
                       FFTSize*subDivFreq*PSF_subpixel,
                           imp_sel.getTitle()+"-Combined");
    }
    return;
    
/* process combined greens to find the pattern, generate pattern, apply to colors*/
/*
    } else if  (label.equals("Find and Apply")) { 
      return;
    } else if  (label.equals("Process Rectangle")) { 
      return;
*/

/* ======================================================================== */
    } else if (label.equals("Simulate")) {
      if (!showSimulDialog()) return;
      String thisTitle="Sim";
      if (imp_sel!=null) thisTitle=imp_sel.getTitle();
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
   if (DEBUG_LEVEL>2) {
          System.out.println("Corr   "+
                             " Ax="+     IJ.d2s(simul_corr[0],5)+
                             " Bx="+     IJ.d2s(simul_corr[1],5)+
                             " Cx="+     IJ.d2s(simul_corr[2],5)+
                             " Ay="+     IJ.d2s(simul_corr[3],5)+
                             " By="+     IJ.d2s(simul_corr[4],5)+
                             " Cy="+     IJ.d2s(simul_corr[5],5)+
                             " Dx="+     IJ.d2s(simul_corr[6],5)+
                             " Ex="+     IJ.d2s(simul_corr[7],5)+
                             " Dy="+     IJ.d2s(simul_corr[8],5)+
                             " Ey="+     IJ.d2s(simul_corr[9],5));
    }
//      quarterHamming=initHamming(FFTSize/2);
      double [] fullHamming=initHamming( FFTSize*PSF_subpixel);
//      for (i=0;i<fullHamming.length;i++) fullHamming[i]=1.0; //
/*  Create high-res boolean simulation array */
//      int patternSize=512; /* TODO - calculate? */
      bPattern= patternGenerator(patternSize,simul_patternType, simul_patternModifier);
      simulation_barray=  simulatePatternFullPattern(bPattern,
                                              simul_freqX1,
                                              simul_freqY1,
                                              simul_phase1,
                                              simul_freqX2,
                                              simul_freqY2,
                                              simul_phase2,
                                              simul_corr, /* apply mesh distortion here */
                                              simul_subdiv,
                                              FFTSize,
                                              centerForG2);
      simul_pixels= extractSimulPatterns (simulation_barray,   // high resolution boolean pattern array
                                                simul_fill, // part of the (center) pixel area being "phptosensitive"
                                               simul_subdiv,    // boolean pixels to real pixels resolution
                                               PSF_subpixel,   // subdivide pixels
                                       FFTSize*PSF_subpixel,    // number of Bayer cells in width of the square selection (half number of pixels)
                                                        0.0,    // selection center, X (in pixels)
                                                        0.0);   // selection center, y (in pixels)

      if ((PSF_subpixel>1) && (OTF_zerofreq_size>=0)) {
        if (colorsToCorrect[5])  simul_pixels=combineCheckerGreens (simul_pixels,   // pixel arrays after oversampleFFTInput() or extractSimulPatterns())
                                                                   PSF_subpixel); // same as used in oversampleFFTInput() - oversampling ratio
      }


      if (DEBUG_LEVEL>1) {
          SDFA_instance.showArrays(simul_pixels[1], FFTSize*subDivFreq*PSF_subpixel, FFTSize*subDivFreq*PSF_subpixel, thisTitle+"-S-red"+simul_patternType+"-"+simul_patternModifier);
          SDFA_instance.showArrays(simul_pixels[4], FFTSize*subDivFreq*PSF_subpixel, FFTSize*subDivFreq*PSF_subpixel, thisTitle+"-S-diag"+simul_patternType+"-"+simul_patternModifier);
          SDFA_instance.showArrays(simul_pixels[5], FFTSize*subDivFreq*PSF_subpixel, FFTSize*subDivFreq*PSF_subpixel, thisTitle+"-S-checker"+simul_patternType+"-"+simul_patternModifier);
      }
      simul_pixels= normalizeAndWindow (simul_pixels, fullHamming);
      if (subDivFreq>1)  simul_pixels= extendFFTInput (simul_pixels,subDivFreq);
      if (DEBUG_LEVEL>1) {
          SDFA_instance.showArrays(simul_pixels[1], FFTSize*subDivFreq*PSF_subpixel, FFTSize*subDivFreq*PSF_subpixel, thisTitle+"-WS-red"+simul_patternType+"-"+simul_patternModifier);
          SDFA_instance.showArrays(simul_pixels[4], FFTSize*subDivFreq*PSF_subpixel, FFTSize*subDivFreq*PSF_subpixel, thisTitle+"-WS-diag"+simul_patternType+"-"+simul_patternModifier);
          SDFA_instance.showArrays(simul_pixels[5], FFTSize*subDivFreq*PSF_subpixel, FFTSize*subDivFreq*PSF_subpixel, thisTitle+"-WS-checker"+simul_patternType+"-"+simul_patternModifier);
      }
      return;
/* ======================================================================== */
    } else if (label.equals("Simulate Quarters")) {
      if (!showSimulDialog()) return;
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      String thisTitle="Sim";
      if (imp_sel!=null) thisTitle=imp_sel.getTitle();

//   boolean hiRes=true;
   int thisSize=FFTSize*(sim_hiRes?PSF_subpixel:1);
   int thisSubdiv=sim_hiRes?PSF_subpixel:1;
   if (DEBUG_LEVEL>2) {
          System.out.println("Corr (1000x)  "+
                             " Ax="+     IJ.d2s(1000*simul_corr[0],5)+
                             " Bx="+     IJ.d2s(1000*simul_corr[1],5)+
                             " Cx="+     IJ.d2s(1000*simul_corr[2],5)+
                             " Ay="+     IJ.d2s(1000*simul_corr[3],5)+
                             " By="+     IJ.d2s(1000*simul_corr[4],5)+
                             " Cy="+     IJ.d2s(1000*simul_corr[5],5)+
                             " Dx="+     IJ.d2s(1000*simul_corr[6],5)+
                             " Ex="+     IJ.d2s(1000*simul_corr[7],5)+
                             " Dy="+     IJ.d2s(1000*simul_corr[8],5)+
                             " Ey="+     IJ.d2s(1000*simul_corr[9],5));
    }

      quarterHamming=initHamming(thisSize/2);
      double [] fullHamming=initHamming(thisSize);

/*  Create high-res boolean simulation array */

      int patternSize=512; /* TODO - calculate? */
      boolean [] bPattern= patternGenerator(patternSize,simul_patternType, simul_patternModifier);
      simulation_barray=  simulatePatternFullPattern(bPattern,
                                              simul_freqX1,
                                              simul_freqY1,
                                              simul_phase1,
                                              simul_freqX2,
                                              simul_freqY2,
                                              simul_phase2,
                                              simul_corr, /* apply mesh distortion here */
                                              simul_subdiv,
                                              FFTSize,
                                              centerForG2);
//      imp_sel= selectImage(2*FFTSize, false); 
      quarter_bayers =new double[9][][];
      quarter_patterns =new double[9][][];
/* Still do for the full area */
      input_bayer=  extractSimulPatterns (simulation_barray,   // high resolution boolean pattern array
                                                 simul_fill, // part of the (center) pixel area being "phptosensitive"
                                               simul_subdiv,   // boolean pixels to real pixels resolution
                                                 thisSubdiv, // subdivide pixels
                                                   thisSize, // number of Bayer cells in width of the square selection (half number of pixels)
//                                                          1, // subdivide pixels
//                                                    FFTSize,    // number of Bayer cells in width of the square selection (half number of pixels)
                                                        0.0,    // selection center, X (in pixels)
                                                        0.0);   // selection center, y (in pixels)
      if (DEBUG_LEVEL>1) SDFA_instance.showArrays(input_bayer,thisSize*subDivFreq, thisSize*subDivFreq,"SIM_FULL");
      input_bayer= normalizeAndWindow (input_bayer,fullHamming);
      if (subDivFreq>1)  input_bayer= extendFFTInput (input_bayer,subDivFreq);
      if (DEBUG_LEVEL>1) SDFA_instance.showArrays(input_bayer,thisSize*subDivFreq, thisSize*subDivFreq,"SIM_FULL_WINDOW");
      pattern =  findPattern(input_bayer[4],
                        thisSize*subDivFreq,
                                 corr_gamma,
                                 corr_sigma,  // pattern detection: high-pass filter (0.0 - none) gamma(PS)
                             diffSpectrCorr,
                             shrinkClusters,
                             multiplesToTry,
                                  deviation,
                            deviation_steps,
                                       true,
                            thisTitle+"SPF");
      int x0=0;
      int y0=0;
      for (iq=0; iq<9;iq++) {
        switch (iq) {
/*
          case 0: x0=-thisSize/4; y0=-thisSize/4;  break; 
          case 1: x0=+thisSize/4; y0=-thisSize/4;  break; 
          case 2: x0=-thisSize/4; y0= thisSize/4;  break; 
          case 3: x0=+thisSize/4; y0= thisSize/4;  break; 
          case 4: x0= 0;         y0= 0;            break; 
          case 5: x0= 0;         y0=-thisSize/4;  break; 
          case 6: x0=-thisSize/4; y0= 0;          break; 
          case 7: x0=+thisSize/4; y0= 0;          break; 
          case 8: x0= 0;         y0= thisSize/4;  break; 
*/
          case 0: x0=-FFTSize/4; y0=-FFTSize/4;  break; 
          case 1: x0=+FFTSize/4; y0=-FFTSize/4;  break; 
          case 2: x0=-FFTSize/4; y0= FFTSize/4;  break; 
          case 3: x0=+FFTSize/4; y0= FFTSize/4;  break; 
          case 4: x0= 0;         y0= 0;          break; 
          case 5: x0= 0;         y0=-FFTSize/4;  break; 
          case 6: x0=-FFTSize/4; y0= 0;          break; 
          case 7: x0=+FFTSize/4; y0= 0;          break; 
          case 8: x0= 0;         y0= FFTSize/4;  break; 


        }
        quarter_bayers[iq]= extractSimulPatterns (simulation_barray,   // high resolution boolean pattern array
                                                         simul_fill, // part of the (center) pixel area being "phptosensitive"
                                                       simul_subdiv,   // boolean pixels to real pixels resolution
                                                         thisSubdiv,
                                                          thisSize/2,    // number of Bayer cells in width of the square selection (half number of pixels)
                                                             2.0*x0,    // selection center, X (in pixels)
                                                             2.0*y0);   // selection center, Y (in pixels)

        if (DEBUG_LEVEL>1) SDFA_instance.showArrays(quarter_bayers[iq],thisSize*subDivFreq/2, thisSize*subDivFreq/2, "Sim"+iq+"X"+x0+"Y"+y0);
        quarter_bayers[iq]= normalizeAndWindow (quarter_bayers[iq], quarterHamming);
        if (subDivFreq>1)  quarter_bayers[iq]= extendFFTInput (quarter_bayers[iq],subDivFreq);
        if (DEBUG_LEVEL>2) SDFA_instance.showArrays(quarter_bayers[iq],thisSize*subDivFreq/2, thisSize*subDivFreq/2, "SimW"+iq+"X"+x0+"Y"+y0);
        quarter_patterns[iq] =  findPattern(quarter_bayers[iq][4],
                                            thisSize*subDivFreq/2,
                                                      corr_gamma,
                                                      corr_sigma,  // pattern detection: high-pass filter (0.0 - none) gamma(PS)
                                                  diffSpectrCorr,
                                                  shrinkClusters,
                                                  multiplesToTry,
                                                       deviation,
                                                 deviation_steps,
                                                            true,
                                                      "SPQ_"+iq);

      }
      if (DEBUG_LEVEL>1) {
          System.out.println("Sim Full area"+
                             " W0x="+     IJ.d2s(pattern[0][0],4)+
                             " W0y="+     IJ.d2s(pattern[0][1],4)+
                             " W0_phase="+IJ.d2s(pattern[0][2],2)+
                             " W1x="+     IJ.d2s(pattern[1][0],4)+
                             " W1y="+     IJ.d2s(pattern[1][1],4)+
                             " W1_phase="+IJ.d2s(pattern[1][2],2));
        for (iq=0; iq<9;iq++) {
          System.out.println("Sim Quarter="+iq+
                             " W0x="+     IJ.d2s(quarter_patterns[iq][0][0],4)+
                             " W0y="+     IJ.d2s(quarter_patterns[iq][0][1],4)+
                             " W0_phase="+IJ.d2s(quarter_patterns[iq][0][2],2)+
                             " W1x="+     IJ.d2s(quarter_patterns[iq][1][0],4)+
                             " W1y="+     IJ.d2s(quarter_patterns[iq][1][1],4)+
                             " W1_phase="+IJ.d2s(quarter_patterns[iq][1][2],2));
        }
      }
/* Filter pattern coefficients to make sure they all match between quadrtants (match to the center one?)*/
     boolean patternsMatchedInitially=matchPatterns(quarter_patterns,quarter_patterns[4]);    // uses pattern in the center quadrant
     if (DEBUG_LEVEL>1) {
          System.out.println(patternsMatchedInitially?"Sim: All quadrant wave vectors matched initially, no correction needed":"Some quadrant wave vectors were adjusted to match");
     }
     
      patternCorr=calcPatternNonLinear(quarter_patterns); // divide results by ,(FFTSize/2)^2
      if (DEBUG_LEVEL>1) { /* increase LEVEL later */
          System.out.println("Pre- (1000x)   "+
                             " Ax="+     IJ.d2s(1000*patternCorr[0]/(FFTSize/2),5)+
                             " Bx="+     IJ.d2s(1000*patternCorr[1]/(FFTSize/2),5)+
                             " Cx="+     IJ.d2s(1000*patternCorr[2]/(FFTSize/2),5)+
                             " Ay="+     IJ.d2s(1000*patternCorr[3]/(FFTSize/2),5)+
                             " By="+     IJ.d2s(1000*patternCorr[4]/(FFTSize/2),5)+
                             " Cy="+     IJ.d2s(1000*patternCorr[5]/(FFTSize/2),5)+
                             " Dx="+     IJ.d2s(1000*patternCorr[6]/(FFTSize/2),5)+
                             " Ex="+     IJ.d2s(1000*patternCorr[7]/(FFTSize/2),5)+
                             " Dy="+     IJ.d2s(1000*patternCorr[8]/(FFTSize/2),5)+
                             " Ey="+     IJ.d2s(1000*patternCorr[9]/(FFTSize/2),5));
      }

      patternCorr=refinePatternNonLinear(quarter_patterns, // [tl,tr,bl,br, center][wv0, wv1][x,y,phase]
                                              patternCorr, //[ax,bx,cx,ay,by,cy]
                                               thisSize/2 ); // distance to quadrats center in sensor pixels ==FFTSize/2

      for (i=0;i<patternCorr.length;i++)patternCorr[i]/= (thisSize/2);

     if (DEBUG_LEVEL>1) { /* increase LEVEL later */
          System.out.println("Corr (1000x)   "+
                             " Ax="+     IJ.d2s(1000*patternCorr[0],5)+
                             " Bx="+     IJ.d2s(1000*patternCorr[1],5)+
                             " Cx="+     IJ.d2s(1000*patternCorr[2],5)+
                             " Ay="+     IJ.d2s(1000*patternCorr[3],5)+
                             " By="+     IJ.d2s(1000*patternCorr[4],5)+
                             " Cy="+     IJ.d2s(1000*patternCorr[5],5)+
                             " Dx="+     IJ.d2s(1000*patternCorr[6],5)+
                             " Ex="+     IJ.d2s(1000*patternCorr[7],5)+
                             " Dy="+     IJ.d2s(1000*patternCorr[8],5)+
                             " Ey="+     IJ.d2s(1000*patternCorr[9],5));
      }

//      simul_corr=patternCorr.clone();

      return;

/* ======================================================================== */
    } else if (label.equals("Find pattern")) {
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      imp_sel= selectImage(2*FFTSize, false);
      input_bayer=splitBayer (imp_sel,equalizeGreens);

      for (i=0;i<4;i++) if (!colorsToCorrect[i]) input_bayer[i]=null; // leave composite greens even if disabled
      input_bayer= normalizeAndWindow (input_bayer, Hamming);
      if (subDivFreq>1)  input_bayer= extendFFTInput (input_bayer,subDivFreq);
      if (DEBUG_LEVEL>2) SDFA_instance.showArrays(input_bayer,FFTSize*subDivFreq, FFTSize*subDivFreq, imp_sel.getTitle());
      boolean pattern4detected=false;
      for (i=0;i<input_bayer.length;i++) if (colorsToCorrect[i] || (i==4)) {
        pattern=null;
        if (i==4) {
          pattern =findPattern(input_bayer[i], FFTSize*subDivFreq, corr_gamma, corr_sigma, diffSpectrCorr, shrinkClusters, multiplesToTry, deviation, deviation_steps, (i==4), imp_sel.getTitle()+"P_"+i);
/* copy patern parameters to globals so they can be viewed/edited */
          if (pattern!=null) {
            simul_freqX1=pattern[0][0];
            simul_freqY1=pattern[0][1];
            simul_phase1=pattern[0][2];
            simul_freqX2=pattern[1][0];
            simul_freqY2=pattern[1][1];
            simul_phase2=pattern[1][2];
            pattern4detected=true;
          }
        } else if (colorsToCorrect[i])  {
/* Just for testing, not needed in final - will use only combined greens ? */
          pattern =findPattern(input_bayer[i], FFTSize*subDivFreq, corr_gamma, corr_sigma, diffSpectrCorr, shrinkClusters, multiplesToTry, deviation, deviation_steps, (i==4), imp_sel.getTitle()+"P_"+i);
        }
        if (pattern!=null) {
           System.out.println("Component: "+componentColorNames[i]+" pattern wave vectors (in original sensor pixels): ");
           System.out.println(" W0x="+     IJ.d2s(pattern[0][0],4)+
                              " W0y="+     IJ.d2s(pattern[0][1],4)+
                              " W0_phase="+IJ.d2s(pattern[0][2],2)+" radians");
           System.out.println(" W1x="+     IJ.d2s(pattern[1][0],4)+
                              " W1y="+     IJ.d2s(pattern[1][1],4)+
                              " W1_phase="+IJ.d2s(pattern[1][2],2)+"radians");
        } else {
           System.out.println("Component: "+componentColorNames[i]+" - FAILURE to find pattern");
        }
      }
      if (!pattern4detected) return;
/* Create simulation pattern (no distortion correction here ). Will use parameters from dual greens - old ones if current failed */
/* generate pattern bitmap image if it does not exist yet */
      if (bPattern==null) bPattern= patternGenerator(patternSize,simul_patternType, simul_patternModifier);
      simulation_barray=  simulatePatternFullPattern(bPattern,
                                              simul_freqX1,
                                              simul_freqY1,
                                              simul_phase1,
                                              simul_freqX2,
                                              simul_freqY2,
                                              simul_phase2,
                                              null, // no mesh distortion here
                                              simul_subdiv,
                                              FFTSize,
                                              true);
      double [][]  sim_pix= extractSimulPatterns (simulation_barray,   // high resolution boolean pattern array
                                                         simul_fill, // part of the (center) pixel area being "phptosensitive"
                                                       simul_subdiv,   // boolean pixels to real pixels resolution
                                                                  1,
                                                            FFTSize,    // number of Bayer cells in width of the square selection (half number of pixels)
                                                                0.0,    // selection center, X (in pixels)
                                                                0.0);   // selection center, y (in pixels)
      sim_pix[4]= normalizeAndWindow (sim_pix[4], Hamming);
      if (subDivFreq>1)  sim_pix[4]= extendFFTInput (sim_pix[4],subDivFreq);
      if (DEBUG_LEVEL>1) {
          SDFA_instance.showArrays(sim_pix[4], FFTSize*subDivFreq, FFTSize*subDivFreq, imp_sel.getTitle()+"SIM"+simul_patternType+"-"+simul_patternModifier);
      }
      double [] model_corr=correlateWithModel (input_bayer[4],  // measured pixel array
                                                   sim_pix[4],  // simulated (model) pixel array)
                                                          0.0,  // double sigma,   // Sigma for high pass filtering
                                           imp_sel.getTitle());
      double [][] convMatrix= {{1.0,-1.0},{1.0,1.0}}; // from pixel WV to greens2 
      double [][] WVpixel={{simul_freqX1,simul_freqY1},{simul_freqX2,simul_freqY2}};
      double [][] WVgreens=matrix2x2_mul(matrix2x2_scale(WVpixel,2.0),matrix2x2_invert(convMatrix)); //matrix2x2_scale(WVpixel,2.0) - bayer decimation
      double [] corrMaxXY={0.0,0.0}; /* TODO - find actual centroid? */

      if (DEBUG_LEVEL>2) {
        System.out.println("  WVpixel[0][0]="+IJ.d2s( WVpixel[0][0],4)+"  WVpixel[0][1]="+IJ.d2s( WVpixel[0][1],4));
        System.out.println("  WVpixel[1][0]="+IJ.d2s( WVpixel[1][0],4)+"  WVpixel[1][1]="+IJ.d2s( WVpixel[1][1],4));
        System.out.println(" WVgreens[0][0]="+IJ.d2s(WVgreens[0][0],4)+" WVgreens[0][1]="+IJ.d2s(WVgreens[0][1],4));
        System.out.println(" WVgreens[1][0]="+IJ.d2s(WVgreens[1][0],4)+" WVgreens[1][1]="+IJ.d2s(WVgreens[1][1],4));
      }



      double contrast= correlationContrast (model_corr,    // square pixel array
                                              WVgreens,    // wave vectors (same units as the pixels array)
                                         corrRingWidth,   // ring (around r=0.5 dist to opposite corr) width
                                          corrMaxXY[0],    //  x0,              // center coordinates
                                          corrMaxXY[1],    //y0,
                                    imp_sel.getTitle());   // title base for optional plots names
      System.out.println("Pattern correlation contrast= "+IJ.d2s(contrast,3)+ ", threshold is "+minCorrContrast);
      return;
/* ======================================================================== */
    } else if (label.equals("Find distorted")) {
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      imp_sel= selectImage(2*FFTSize, false); 
      input_bayer=splitBayer (imp_sel,equalizeGreens);

 /* returns array of 3 arrays: first two are 3-element wave vectors (x,y,phase), last - 3-rd order correction coefficients */
     double[][] distortedPattern= findPatternDistorted(   input_bayer, // pixel array to process (no widowing!)
                                                           corr_gamma,
                                                           corr_sigma,  // pattern detection: high-pass filter (0.0 - none) gamma(PS)
                                                       diffSpectrCorr,
                                                       shrinkClusters,
                                                       multiplesToTry,
                                                            deviation,
                                                      deviation_steps,
                                                                 true, // boolean greens, // this is a pattern for combined greens (diagonal), adjust results accordingly
                                                  imp_sel.getTitle()); // title prefix to use for debug  images

      if (DEBUG_LEVEL>1) { /* increase LEVEL later */
          System.out.println("Pattern: "+
                             " W0x="+     IJ.d2s(distortedPattern[0][0],4)+
                             " W0y="+     IJ.d2s(distortedPattern[0][1],4)+
                             " W0_phase="+IJ.d2s(distortedPattern[0][2],2)+
                             " W1x="+     IJ.d2s(distortedPattern[1][0],4)+
                             " W1y="+     IJ.d2s(distortedPattern[1][1],4)+
                             " W1_phase="+IJ.d2s(distortedPattern[1][2],2));

          System.out.println("Corr (1000x)   "+
                             " Ax="+     IJ.d2s(1000*distortedPattern[2][0],5)+
                             " Bx="+     IJ.d2s(1000*distortedPattern[2][1],5)+
                             " Cx="+     IJ.d2s(1000*distortedPattern[2][2],5)+
                             " Ay="+     IJ.d2s(1000*distortedPattern[2][3],5)+
                             " By="+     IJ.d2s(1000*distortedPattern[2][4],5)+
                             " Cy="+     IJ.d2s(1000*distortedPattern[2][5],5)+
                             " Dx="+     IJ.d2s(1000*distortedPattern[2][6],5)+
                             " Ex="+     IJ.d2s(1000*distortedPattern[2][7],5)+
                             " Dy="+     IJ.d2s(1000*distortedPattern[2][8],5)+
                             " Ey="+     IJ.d2s(1000*distortedPattern[2][9],5));
      }
/* copy patern parameters to globals so they can be viewed/edited - use center pattern, not the full one - center is used when calculating corrections*/

      simul_freqX1=distortedPattern[0][0];
      simul_freqY1=distortedPattern[0][1];
      simul_phase1=distortedPattern[0][2];
      simul_freqX2=distortedPattern[1][0];
      simul_freqY2=distortedPattern[1][1];
      simul_phase2=distortedPattern[1][2];
      simul_corr=distortedPattern[2].clone();
      return;
/* ======================================================================== */
    } else if (label.equals("Process files")) {
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      JFileChooser fc= new JFileChooser();
      fc.setMultiSelectionEnabled(true);

      if (dir==null) { // global
	  String sdir = OpenDialog.getDefaultDirectory();
	  if (sdir!=null)
	      dir = new File(sdir);
      }
      if (dir!=null)
	  fc.setCurrentDirectory(dir);
      int returnVal = fc.showOpenDialog(IJ.getInstance());
      if (returnVal!=JFileChooser.APPROVE_OPTION)
      return;
      File[] files = fc.getSelectedFiles();
      if (files.length==0) { // getSelectedFiles does not work on some JVMs
	  files = new File[1];
	  files[0] = fc.getSelectedFile();
      }
      String path = fc.getCurrentDirectory().getPath()+Prefs.getFileSeparator();
      dir = fc.getCurrentDirectory();
      System.out.println("path= "+path+", files:");

      filenames=new String[files.length]; // global
      filePatternMaps=new double[filenames.length][][][][];
      fileTilesMaps=  new boolean[filenames.length][][]; // map of 2*FFTSize x 2*FFTSize squares with 2*FFTOverlap step, true means that that tile can be used for PSF
      int nFile;
      for (nFile=0;nFile<files.length;nFile++) {
        filenames[nFile]= files[nFile]+"";
        filePatternMaps[nFile]=null;
        fileTilesMaps[nFile]=null;
      }
      for (nFile=0;nFile<files.length;nFile++) {
        System.out.println(nFile+": "+filenames[nFile]);
        imp_sel=jp4_instance.open(
          "", // path,
            filenames[nFile],
          "",  //arg - not used in JP46 reader
          true, // un-apply camera color gains
         imp_sel); // reuse the same image window
         if (bPattern==null) bPattern= patternGenerator(patternSize,simul_patternType, simul_patternModifier);
         filePatternMaps[nFile]=scanImageForPatterns(imp_sel,
                                               mapFFTSize,
                                               corr_gamma,
                                               corr_sigma,
                                           diffSpectrCorr,
                                           shrinkClusters,
                                           multiplesToTry,
                                                deviation,
                                          deviation_steps,
                                                 bPattern,
                                           model_highpass,
                                            corrRingWidth,
                                          minCorrContrast,
                                                        1); // debug level to be used while scanning cells
         fileTilesMaps[nFile]= remapSquares(filePatternMaps[nFile], // [][]map of either null or 2 wave vectors
                                                      mapFFTSize/2, // step of initial map
                                                        mapFFTSize, // size of square used in scanning of initial map (should be multiple of map step)
                                                        FFTOverlap, // step of the new map (should be multiple of map step)
                                                           FFTSize); // size of square used in sthe new map (should be multiple of map step)


      }


/* Show number of files that cover each squqre*/
      int nTileX,nTileY;
      ImageProcessor ip_sel=imp_sel.getProcessor();
      float [] src_pixels;
      src_pixels=(float[]) ip_sel.getPixels();
//      float [] masked_pixels=new float[src_pixels.length];
//      double src_pixels_average=0.0;
//      int npX,npY;
      int imp_sel_width= imp_sel.getWidth();
      int imp_sel_height=imp_sel.getHeight();
//      int nimg;
//      int mTileY=filePatternMaps[0].length;
//      int mTileX=filePatternMaps[0][0].length;
      int ix,iy;

/* TODO: use correlation contrast as a measure of the "quality" of the cell in a file, then use the "best" file for each square

/* Show map other way, make same size as the image */
      float [] floatOrigMap=new float [src_pixels.length];
      int [][] intOriginalMap= mapToNumFiles (filePatternMaps,2);
      int tileSize=mapFFTSize;
      for (i=0;i<floatOrigMap.length;i++) {
        iy=i / imp_sel_width;
        ix=i % imp_sel_width;
        nTileY=iy/tileSize;
        nTileX=ix/tileSize;
        if ((nTileY<intOriginalMap.length) && (nTileX<intOriginalMap[0].length)) floatOrigMap[i]= (float) intOriginalMap[nTileY][nTileX];
        else floatOrigMap[i]= 0.0F;
      }
      ImageProcessor ip_masked1=new FloatProcessor(imp_sel_width,imp_sel_height);
      ip_masked1.setPixels(floatOrigMap);
      ip_masked1.resetMinAndMax();
      ImagePlus imp_masked1=  new ImagePlus("num_img1", ip_masked1);
      imp_masked1.show();
/* Show map for 2*FFTSize squares with 2*FFTOverlap step */
      float [] floatRemap=new float [src_pixels.length];
      int [][] intRemap= mapToNumFiles (fileTilesMaps,FFTSize/FFTOverlap);
      tileSize=2*FFTOverlap;
      for (i=0;i<floatRemap.length;i++) {
        iy=i / imp_sel_width;
        ix=i % imp_sel_width;
        nTileY=iy/tileSize;
        nTileX=ix/tileSize;
        if ((nTileY<intRemap.length) && (nTileX<intRemap[0].length)) floatRemap[i]= (float) intRemap[nTileY][nTileX];
        else floatRemap[i]= 0.0F;
      }
      ImageProcessor ip_masked2=new FloatProcessor(imp_sel_width,imp_sel_height);
      ip_masked2.setPixels(floatRemap);
      ip_masked2.resetMinAndMax();
      ImagePlus imp_masked2=  new ImagePlus("num_img2", ip_masked2);
      imp_masked2.show();

      return;

/* ======================================================================== */

    } else if (label.equals("Map image")) {
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      imp_sel = WindowManager.getCurrentImage();
      if (imp_sel==null){
        IJ.showMessage("Error","There are no images open\nProcess canceled");
        return;
      }
      if (bPattern==null) bPattern= patternGenerator(patternSize,simul_patternType, simul_patternModifier);
      double [][][][] patternMap=scanImageForPatterns(imp_sel,
                                                   mapFFTSize,
                                                   corr_gamma,
                                                   corr_sigma,
                                               diffSpectrCorr,
                                               shrinkClusters,
                                               multiplesToTry,
                                                    deviation,
                                              deviation_steps,
                                                     bPattern,
                                               model_highpass,
                                                corrRingWidth,
                                              minCorrContrast,
                                                            1); // debug level to be used while scanning cells
      int nTileX,nTileY;
      if (DEBUG_LEVEL>1) {
         int numPatternCells=0;
         for (nTileY=0;nTileY<patternMap.length;nTileY++) for (nTileX=0;nTileX<patternMap[0].length;nTileX++) if (patternMap[nTileY][nTileX]!=null) numPatternCells++;
         System.out.println("Finished mapping, created array["+patternMap.length+"]["+patternMap[0].length+"][][], "+
                             numPatternCells+" cells (of "+(patternMap.length*patternMap[0].length)+") with pattern detected");
      }
/* Filter results based on correlation with the actual pattern */

/* Show masked image - where it failed to find pattern*/
      double maskingContrast=0.1;
      ImageProcessor ip_sel=imp_sel.getProcessor();
      float [] src_pixels;
      src_pixels=(float[]) ip_sel.getPixels();
      float [] masked_pixels=new float[src_pixels.length];
      double src_pixels_average=0.0;
      for (i=0;i<src_pixels.length; i++) src_pixels_average+=src_pixels[i];
      src_pixels_average/=src_pixels.length;
      double aa=src_pixels_average*(1.0-maskingContrast);
      for (i=0;i<src_pixels.length; i++) masked_pixels[i]= (float) (src_pixels[i]*maskingContrast+aa);
      int npX,npY;
      int imp_sel_width= imp_sel.getWidth();
      int imp_sel_height=imp_sel.getHeight();
      for (nTileY=0;nTileY<patternMap.length;nTileY++) for (nTileX=0;nTileX<patternMap[0].length;nTileX++) if (patternMap[nTileY][nTileX]!=null) {
        for (npY=0;npY<2*mapFFTSize;npY++) for (npX=0;npX<2*mapFFTSize;npX++) {
          i=(nTileY*mapFFTSize+npY)*imp_sel_width + (nTileX*mapFFTSize+npX);
          masked_pixels[i]= src_pixels[i];
        }
      }
      ImageProcessor ip_masked=new FloatProcessor(imp_sel_width,imp_sel_height);
      ip_masked.setPixels(masked_pixels);
      ip_masked.resetMinAndMax();
      ImagePlus imp_masked=  new ImagePlus(imp_sel.getTitle()+"_masked", ip_masked);
      imp_masked.show();
      return;

/* ======================================================================== */

    } else if (label.equals("Map PSF")) {
      int loop_debug_level=1;
      if (!showPSFDialog()) return;
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      imp_sel = WindowManager.getCurrentImage();
      if (imp_sel==null){
        IJ.showMessage("Error","There are no images open\nProcess canceled");
        return;
      }
      if (bPattern==null) bPattern= patternGenerator(patternSize,simul_patternType, simul_patternModifier);
      double [][][][] patternMap=scanImageForPatterns(imp_sel,
                                                   mapFFTSize,
                                                   corr_gamma,
                                                   corr_sigma,
                                               diffSpectrCorr,
                                               shrinkClusters,
                                               multiplesToTry,
                                                    deviation,
                                              deviation_steps,
                                                     bPattern,
                                               model_highpass,
                                                corrRingWidth,
                                              minCorrContrast,
                                                            1); // debug level to be used while scanning cells
      int nTileX,nTileY;
      int numPatternCells=0;
      for (nTileY=0;nTileY<patternMap.length;nTileY++) for (nTileX=0;nTileX<patternMap[0].length;nTileX++) if (patternMap[nTileY][nTileX]!=null) numPatternCells++;
      if (DEBUG_LEVEL>1) {
         System.out.println("Finished mapping, created array["+patternMap.length+"]["+patternMap[0].length+"][][], "+
                             numPatternCells+" cells (of "+(patternMap.length*patternMap[0].length)+") with pattern detected");
      }
/* Filter results based on correlation with the actual pattern */
      boolean [][]   PSFBooleanMap; // map of 2*FFTSize x 2*FFTSize squares with 2*FFTOverlap step, true means that that tile can be used for PSF
      PSFBooleanMap= remapSquares(patternMap, // [][]map of either null or 2 wave vectors
                         mapFFTSize/2, // step of initial map
                           mapFFTSize, // size of square used in scanning of initial map (should be multiple of map step)
                           FFTOverlap, // step of the new map (should be multiple of map step)
                              FFTSize); // size of square used in sthe new map (should be multiple of map step)
//      nTileY=PSFBooleanMap.length;
//      nTileX=PSFBooleanMap[0].length;
      numPatternCells=0;
      for (nTileY=0;nTileY<PSFBooleanMap.length;nTileY++) for (nTileX=0;nTileX<PSFBooleanMap[0].length;nTileX++) if (PSFBooleanMap[nTileY][nTileX]) numPatternCells++;
      if (DEBUG_LEVEL>1) {
         System.out.println("Remapped for PSF measurment- converted to an array["+PSFBooleanMap.length+"]["+PSFBooleanMap[0].length+"], "+
                             numPatternCells+" cells (of "+(PSFBooleanMap.length*PSFBooleanMap[0].length)+") with pattern detected");
      }
//      double [][][][] PSFKernelMap=new double[PSFBooleanMap.length][PSFBooleanMap[0].length][][];
      PSFKernelMap=new double[PSFBooleanMap.length][PSFBooleanMap[0].length][][]; //PSFKernelMap - global
      for (nTileY=0;nTileY<PSFKernelMap.length;nTileY++) for (nTileX=0;nTileX<PSFKernelMap[0].length;nTileX++) PSFKernelMap[nTileY][nTileX]=null;
      DEBUG_LEVEL=loop_debug_level;
      int ncell=1;
      int x0,y0;
//    double [] test01;
//    double [][] test02;
      for (nTileY=0;nTileY<PSFBooleanMap.length;nTileY++) for (nTileX=0;nTileX<PSFBooleanMap[0].length;nTileX++) if (PSFBooleanMap[nTileY][nTileX]) {
        if (MASTER_DEBUG_LEVEL>0) IJ.showStatus("Processing tile["+nTileY+"]["+nTileY+"] ("+ncell+" of "+numPatternCells+")");
        y0=FFTOverlap*2*nTileY;
        x0=FFTOverlap*2*nTileX;
        PSFKernelMap[nTileY][nTileX]=getPSFKernels (          imp_sel,
                                                            2*FFTSize,   // size in pixels (twice FFTSize)
                                                                   x0,      // top left corner X (pixels)
                                                                   y0,      // top left corner Y (pixels)
                                                       equalizeGreens,
                                                           corr_gamma, // gamma to use for power spectrum for correlation
                                                           corr_sigma, // high-pass gaussian filter sigma when correlating power spectrum
                                                       diffSpectrCorr, // maximal distance between maximum on spectrum and predicted maximum on autocorrelation of gamma(|spectrum|)
                                                       shrinkClusters, // Shrink clusters by this ratio (remove lowest) after initial separation
                                                       multiplesToTry, // try this number of maximums proportionally farther from 0,0 than the two closest (increase precision)
                                                            deviation, // when looking for maximums - maximal distance from predicted from the lower order one
                                                      deviation_steps, // maximal iterations when looking for local maximum
                                                             bPattern, // pattern bitmask (calculate once)
                                                         simul_subdiv, // simulatin boolean pixels per sensor pixel
                                                          centerForG2,
                                                           subDivFreq, // add zeros arond (increase frequency resolution, probably will not be used here)
                                                 initHamming(FFTSize), //=initHamming( fft_size) calculate once
                                    initHamming(FFTSize*PSF_subpixel), //=initHamming( fft_size*subpixel);
                                                         PSF_subpixel, // use finer grid than actual pixels 
                                                    OTF_zerofreq_size,
                                                           simul_fill, // part of the (center) pixel area being "photosensitive"
                                                      colorsToCorrect, // color channnels to process (should be at least 6 long)
                                                         deconvInvert, //  fraction of the maximal value to be used to limit zeros
// next 3 used for filtering aliases 
                                                         OTF_smoothPS, // 0 - none, otherwise Gauss width = FFT size/2/smoothPS
                                                   OTF_threshold_high,  // reject completely if energy is above this part of maximal
                                                    OTF_threshold_low,  // leave intact if energy is below this part of maximal
// Next 6 were used to filter out X-shaped artifacts on the PSF (when used plain checkerboard pattern), probbaly not needed anymore
                                                      XMASK_threshold, // if 0 use normalize amplitude, if 0..1 - make binary: 1.0 if > threshold, 0.0 - otherwise -1 - disable mask
                                                          XMASK_radius, // low-pass result with low pass filter (should be later defined automatically)
                                                           XMASK_hsize,  // 2d histogram size (size/2 probably a good guess),
                                                      XMASK_percentile, // use this percentile (0.0..1.0)) value for given radius as a target
                                                         XMASK_maxGain, // maximal gain for low components
                                                      XMASK_exaggerate, // exaggerate correction mask with Math.pow()) 
                                                                     5, // int referenceComp, // number of color component to reference lateral chromatic aberration to (now 4 - checkerboard greens)
                                                        PSF_useHamming,  // if false - use gausian
                                                       PSF_minContrast,    // minimal instance contrast to use in binning
                                                        PSF_windowFrac,     // reduce the PSF cell size to this part of the area connecting first negative clones
                                                           PSF_symm180,
                                                   PSF_ignoreChromatic,
                                               PSF_enableModelSubtract,  // generate/subtract gaussian models (not needed if no overlap between pos/neg clones)
                                                    PSF_smoothSeparate, // low-pass filter multiple opposite-sign PSF instaces for separation (width relative to distance to the opposite sign PSF)
                                                 PSF_thresholdSeparate,  // threshold for locating zero-crossing
                                                         PSF_topCenter,     // consider only points above that fraction of the local max to find centroid
                                                   PSF_removeNegtative,  // remove PSF negative values  when separating composite PSF (will need low-pass filtering)
                                                     PSF_sigmaToRadius,   // 0.4;  variable-sigma blurring to reduce high frequencies more for the pixels farther from the PSF center
                                                      PSF_wings_energy, // fraction of energy in the pixels to be used
                                               PSF_wings_ellipse_scale);
/* TODO: Add normalization of kernels */
// 256x256 here
/*
        if (MASTER_DEBUG_LEVEL>1) {
          SDFA_instance.showArrays(PSFKernelMap[nTileY][nTileX],
          (int) Math.sqrt(kernelLength(PSFKernelMap[nTileY][nTileX])),
          (int) Math.sqrt(kernelLength(PSFKernelMap[nTileY][nTileX])),
          imp_sel.getTitle()+"-full-X"+x0+"Y"+y0);
        }
*/

        if (kernelLength(PSFKernelMap[nTileY][nTileX])>(PSF_kernel_size*PSF_kernel_size)) PSFKernelMap[nTileY][nTileX]=resizeForFFT(PSFKernelMap[nTileY][nTileX],PSF_kernel_size); // shrink before normalizing
/*
        if (MASTER_DEBUG_LEVEL>1) SDFA_instance.showArrays(PSFKernelMap[nTileY][nTileX],
          PSF_kernel_size,
          PSF_kernel_size,
          imp_sel.getTitle()+"X"+x0+"Y"+y0);
*/
        normalizeKernel(PSFKernelMap[nTileY][nTileX]); // in-place
/*
        if (MASTER_DEBUG_LEVEL>1) SDFA_instance.showArrays(PSFKernelMap[nTileY][nTileX],
          PSF_kernel_size,
          PSF_kernel_size,
          imp_sel.getTitle()+"-norm-X"+x0+"Y"+y0);
*/

        if (kernelLength(PSFKernelMap[nTileY][nTileX])<(PSF_kernel_size*PSF_kernel_size)) PSFKernelMap[nTileY][nTileX]=resizeForFFT(PSFKernelMap[nTileY][nTileX],PSF_kernel_size); // expand after normalizing
        ncell++;
      }
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
 /* merge 2-d array of kernels into large 2d-array of pixels (i.e. to be shown with showBayer()) */
//      int mergedWidth= PSFKernelMap[0].length*PSF_kernel_size;
//      int mergedHeight=PSFKernelMap.length*PSF_kernel_size;

      SDFA_instance.showImageStack( mergeKernelsToStack(PSFKernelMap), imp_sel.getTitle()+"-PSF_KERNEL");
      return;
    } else if (label.equals("Read Deconvs")) {
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      ImagePlus imp_kernels = WindowManager.getCurrentImage();
      if (imp_kernels==null){
        IJ.showMessage("Error","There is no image to process");
        return;
      }
      ImageProcessor ip_kernels=imp_kernels.getProcessor();
      double [][][][] kMap;
      if (imp_kernels.getStackSize()>1) {
        kMap=splitSquareKernelsFromStack(imp_kernels.getStack(), // Image stack, each slice consists of square kernels of one channel
                                               RPSF_kernel_size, // size of each kernel (should be square)
                                 channelNumbers(colorsToCorrect));  // four-element array of color channel numbers (usually {1,2,5,-1})
      } else {
        kMap= splitSquareKernelsFromOne((float[]) ip_kernels.getPixels(), // flattened kernels, same as produced by mergeSquareKernelsToOne()
                                                        RPSF_kernel_size, // size of each kernel (should be square)
                                                  imp_kernels.getWidth(), // width of the flattened kernels, in pixels (should be multiple of 2*size)
                                         channelNumbers(colorsToCorrect));  // four-element array of color channel numbers (usually {1,2,5,-1})
      }
      if (DEBUG_LEVEL>2) SDFA_instance.showImageStack( mergeKernelsToStack(kMap), imp_kernels.getTitle()+"-test");
      rPSFKernelMap=kMap;
      if (DEBUG_LEVEL>1) System.out.println ( "Read image "+imp_kernels.getTitle()+" as reversed PSF kernels (rPSFKernelMap), "+kMap[0].length+"x"+kMap.length);
      return;

    } else if (label.equals("Read Gaussian")) {
      imp_gaussian=WindowManager.getCurrentImage();
      if ((imp_gaussian==null) || (imp_gaussian.getStackSize()<3)) {
        IJ.showMessage("Error","Please select image stack of 3 (float) slices (r,g,b), resulted from convolving input bayer image with Gaussian kernels");
        return;
      }
      if (DEBUG_LEVEL>1) System.out.println ( "Read image "+imp_gaussian.getTitle()+" as the image convolved with Gaussian kernels");

      return;
    } else if (label.equals("Combine pair")) {
      if ((imp_gaussian==null) || (imp_gaussian.getStackSize()<3)) {
        IJ.showMessage("Error","Wrong/empty Gaussian convolved image, please use 'Read Gaussian' command first");
        return;
      }
      imp_convolved=WindowManager.getCurrentImage();
      if ((imp_convolved==null) || (imp_convolved.getStackSize()<3)) {
        IJ.showMessage("Error","Please select image stack of 3 (float) slices (r,g,b), resulted from convolving input bayer image with the reversed PSF kernels");
        return;
      }
      if (DEBUG_LEVEL>1) System.out.println ( "Read image "+imp_convolved.getTitle()+" as the image convolved with the reversed PSF kernels");
      ImageStack stack_convolved=imp_convolved.getStack();
      ImageStack stack_gaussian=imp_gaussian.getStack();
      if (!showCombinePairDialog()) return;
      int imgWidth= imp_convolved.getWidth();
      int imgHeight=imp_convolved.getHeight();
      float [] diffGreens=new float[imgWidth*imgHeight];

/* find number of the green channel - should be called "green", if none - use last */
     int greenChn=2;
     for (i=0;i<3;i++) if (stack_convolved.getSliceLabel(i+1).equals("green")){
       greenChn=i;
       break;
     }

//      int thersholdColor=2; // green
      float d;
      float max=0.0f;
      float [] hipassPixels=(float[]) stack_convolved.getPixels(greenChn+1);
      float [] lopassPixels=(float[]) stack_gaussian.getPixels(greenChn+1);
      for (i=0;i<lopassPixels.length;i++) {
        d=hipassPixels[i]-lopassPixels[i];
        diffGreens[i]=d*d;
        if (max<lopassPixels[i]) max=lopassPixels[i];
      }
      max*=(float) NONLIN_THRESHOLD;
      for (i=0;i<lopassPixels.length;i++) {
        diffGreens[i]/=(float) Math.max(max,lopassPixels[i]);
      }
      FloatProcessor fp_diffGreens = new FloatProcessor(imgWidth, imgHeight, diffGreens, null);
      GaussianBlur gb = new GaussianBlur();
      gb.blurFloat(fp_diffGreens, NONLIN_SIGMA, NONLIN_SIGMA, 0.01);

//      ImagePlus imp_mask=new ImagePlus("mask", fp_diffGreens);
//      fp_diffGreens.resetMinAndMax();
//      imp_mask.show();

       d= (float) ( 1.0/(NONLIN_MAX-NONLIN_MIN));
       if (NONLIN_MAX>NONLIN_MIN) {
         for (i=0;i<diffGreens.length;i++) {
           if (diffGreens[i]<NONLIN_MIN) diffGreens[i]=0.0f;
           else if (diffGreens[i]>NONLIN_MAX) diffGreens[i]=1.0f;
           else diffGreens[i]=d*(diffGreens[i]- (float) NONLIN_MIN);
         }
       }

      ImagePlus imp_mask=new ImagePlus("mask", fp_diffGreens);
      fp_diffGreens.resetMinAndMax();
      imp_mask.show();

/* Combine 2 stacks and a mask */
      ImageStack stack_combo= combineStacksWithMask (stack_gaussian,
                                                    stack_convolved, 
                                                         diffGreens);
      ImagePlus imp_stack_combo = new ImagePlus(imp_convolved.getTitle()+"combo-rgb", stack_combo);
      imp_stack_combo.getProcessor().resetMinAndMax();
      imp_stack_combo.show();
      return;


    } else if (label.equals("Test")) {
      ImagePlus imp_test = WindowManager.getCurrentImage();
      if (imp_test==null){
        IJ.showMessage("Error","There is no image to process");
        return;
      }
      ImageProcessor ip_test=imp_test.getProcessor();
      float [] test_fpixels=(float[]) ip_test.getPixels();
      double [] test_pixels=new double [test_fpixels.length];
      int test_size=(int) Math.sqrt(test_fpixels.length);
      for (i=0;i<test_pixels.length;i++) test_pixels[i]=test_fpixels[i];
      test_pixels= normalizeAndWindow (test_pixels, initHamming(test_size),false);

      double [][][] otf=psf2otf (test_pixels,true);
      double [][] aphase=new double [2][test_fpixels.length];
      int otf_size=otf[0].length;
      int ix,iy,is;
      for (i=0;i<aphase[0].length;i++) {
        ix= (i + otf_size/2) % otf_size;
        iy= (i/otf_size + otf_size/2) % otf_size;
        if (iy > otf_size/2) {
          iy= otf_size-iy;
          ix=(otf_size-ix) % otf_size;
          is=-1;
        } else is=1;
        aphase[0][i]=   otf[iy][ix][0];
        aphase[1][i]=is*otf[iy][ix][1];
      }
      SDFA_instance.showArrays(aphase, otf_size, otf_size,imp_test.getTitle()+"-OTF");

      double [] restoredPixels= ampPhaseToSpace(otf,true);
      SDFA_instance.showArrays(restoredPixels, otf_size, otf_size,imp_test.getTitle()+"-BACK");
/* Create Bayer pattern */
      float [][] bayer_patterns=new float [4][test_size*test_size];
      for (iy=0;iy<test_size;iy++) for (ix=0;ix<test_size;ix++) {
        bayer_patterns[0][iy*test_size+ix]=(((iy % PSF_subpixel)==0) && ((ix % PSF_subpixel)==(PSF_subpixel/2) ))?1.0f:0.0f;
        bayer_patterns[1][iy*test_size+ix]=((((iy % PSF_subpixel)==0) && ((ix % PSF_subpixel)==0 )) || (((iy % PSF_subpixel)==(PSF_subpixel/2)) && ((ix % PSF_subpixel)==(PSF_subpixel/2) )))?1.0f:0.0f;
        bayer_patterns[2][iy*test_size+ix]=(((iy % PSF_subpixel)==(PSF_subpixel/2)) && ((ix % PSF_subpixel)==0 ))?1.0f:0.0f;
        bayer_patterns[3][iy*test_size+ix]=(((iy % PSF_subpixel)==0) && ((ix % PSF_subpixel)==0 ))?1.0f:0.0f;
      }
      ImageStack bayer_stack=new ImageStack(test_size,test_size);
      bayer_stack.addSlice("Bayer_red",    bayer_patterns[0]);
      bayer_stack.addSlice("Bayer_greens", bayer_patterns[1]);
      bayer_stack.addSlice("Bayer_blue",   bayer_patterns[2]);
      bayer_stack.addSlice("Bayer_zero",   bayer_patterns[3]);
      SDFA_instance.showImageStack(bayer_stack, "Bayer_stack_RGB_zero");
      return;
/* ======================================================================== */
    } else if (label.equals("Interpolate Kernels")) {

//      int loop_debug_level=1;
      if (!showInterpolateKernelsDialog()) return;
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
      ImageStack interpolatedStack= interpolateKernelStack(imp_kernels.getStack(), // Image stack, each slice consists of square kernels of one channel
                                                               INTERPOLATE_INSIZE, // size of each kernel (should be square)
                                                                 INTERPOLATE_STEP, // number of subdivisions form input to output
                                                              INTERPOLATE_ADDLEFT, // add this number of kernel columns to the output on the left of the existent/interpolated
                                                               INTERPOLATE_ADDTOP, // add this number of kernel rows to the output above the existent/interpolated
                                                             INTERPOLATE_ADDRIGHT, // add this number of kernel columns to the output on the right of the existent/interpolated
                                                            INTERPOLATE_ADDBOTTOM, // add this number of kernel rows to the output below the existent/interpolated
                                                          INTERPOLATE_EXTRAPOLATE,  // 0 - duplicate, 1.0 - extrapolate outside of the known kernels
                                                                    UPDATE_STATUS); // update status info

      ImagePlus imp_interpolated_stack = new ImagePlus(imp_kernels.getTitle()+"-"+INTERPOLATE_STEP+ "X-interpolated", interpolatedStack);
      imp_interpolated_stack.getProcessor().resetMinAndMax();
      imp_interpolated_stack.show();
      return;
/* ======================================================================== */
    } else if (label.equals("Inverse Stack")) {

//      int loop_debug_level=1;
      if (!showInverseStackDialog()) return;
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
      double [] sigmas=null;
      if (INVERSE_FILTER) {
        sigmas=new double[3];
        sigmas[0]=DECONCV_BLUR_INDIVIDUAL;
        sigmas[1]=DECONCV_BLUR_INDIVIDUAL;
        sigmas[2]=DECONCV_BLUR_CHECKER;
      }
      ImageStack invertedStack= reversePSFKernelStack(imp_kernels.getStack(), //  stack of 3 32-bit (float) images, made of square kernels
                                                         INVERSE_DIRECT_SIZE, // size (side of square) of direct PSF kernel
                                                        INVERSE_REVERSE_SIZE, // size (side of square) of reverse PSF kernel
                                                            OTF_deconvInvert,    // deconvInvert
                                                           OTF_cutoff_energy,  // OTF_cutoff_energy
                                                           OTF_ellipse_scale,  // ellipse mask size relative to the cluster
                                                           OTF_ellipse_gauss,  // use Gauss instead of Hamming for ellipse mask
                                                           PSF_cutoff_energy,  // PSF_cutoff_energy
                                                           PSF_ellipse_scale,  // ellipse mask size relative to the cluster
                                                     RPSF_min_mask_threshold,  // zero output element if elliptical Gauss mask is below this threshold
// Optional variable-sigma blurring parameters
                                                                      sigmas, // array of sigmas in the center, matching stacks sequence. Null if no blurring is needed
                                                        RPSF_var_sigma_scale,  // scale sigma in the center when using variable sigma
                                                        RPSF_sigma_to_radius,  // sigma-to-radius ratio (0.0 to disable variable blur)
                                                               UPDATE_STATUS);  // update status info

      ImagePlus imp_invertedStack = new ImagePlus(imp_kernels.getTitle()+"-rPSF", invertedStack);
      imp_invertedStack.getProcessor().resetMinAndMax();
      imp_invertedStack.show();
      return;

/* ======================================================================== */
    } else if (label.equals("Gaussian Stack")) {
//      int loop_debug_level=1;
//      if (!showInverseStackDialog()) return;
      if (!showGaussianStackDialog()) return;
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
      double [] sigmas=null;
      sigmas=new double[3];
      sigmas[0]=DECONCV_BLUR_INDIVIDUAL;
      sigmas[1]=DECONCV_BLUR_INDIVIDUAL;
      sigmas[2]=DECONCV_BLUR_CHECKER;
      ImageStack gaussianStack=  generateGaussianStackFromDirect(imp_kernels.getStack(), // stack of 3 32-bit (float) images, made of square kernels
                                                                    INVERSE_DIRECT_SIZE, // size (side of square) of direct PSF kernel
                                                                   INVERSE_REVERSE_SIZE, // size (side of square) of the Gaussian kernel
                                                                      PSF_cutoff_energy,  // OTF_cutoff_energy
                                                                                 sigmas, // array of sigmas in the center, matching stacks sequence. Null if no blurring is needed
                                                                          UPDATE_STATUS);  // update status info

      ImagePlus imp_gaussianStack = new ImagePlus(imp_kernels.getTitle()+"-Gaussian", gaussianStack);
      imp_gaussianStack.getProcessor().resetMinAndMax();
      imp_gaussianStack.show();
      return;
/* ======================================================================== */
    } else if (label.equals("Split Image")) {
//      int loop_debug_level=1;
      if (!showSplitBayerToStackDialog()) return;
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      ImagePlus imp_src = WindowManager.getCurrentImage();
      if (imp_src==null){
        IJ.showMessage("Error","Bayer Image required");
        return;
      }
      ImageStack sourceStack= bayerToStack(imp_src, // source bayer image, linearized, 32-bit (float))
                                  SPLIT_OVERSAMPLE, // multiple samples per pixel in each direction (2)
                                     SPLIT_ADDLEFT, // add this number of scan lines above the image (reducing border effects)
                                      SPLIT_ADDTOP, // 
                                    SPLIT_ADDRIGHT, // 
                                   SPLIT_ADDBOTTOM);
      ImagePlus imp_srcStack = new ImagePlus(imp_src.getTitle()+"-EXP", sourceStack);
      imp_srcStack.getProcessor().resetMinAndMax();
      imp_srcStack.show();
      return;

/* ======================================================================== */
    } else if (label.equals("Debayer Image")) {
//      int loop_debug_level=1;
      if (!aliasScissorsStackDialog()) return;
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      ImagePlus imp_src = WindowManager.getCurrentImage();
      if ((imp_src==null) || (imp_src.getStackSize()<3)){
        IJ.showMessage("Error","Bayer image stack required");
        return;
      }
    
      ImageStack imageStack= aliasScissorsStack(imp_src.getStack(),  // stack with 3 colors/slices with the image
                                                  DEBAYER_FFT_SIZE, // 64 - fft size
                                                  SPLIT_OVERSAMPLE, // number of image pixels/ sensor pixels (each direction) == 2
                                                 DEBAYER_THRESHOLD, // no high frequencies - use default uniform filter
                                               DEBAYER_WIDTH_GREEN, // realtive width of the low-pass filter for greens (1.0 goes to zero at exactrly the closest alias )
                                             DEBAYER_WIDTH_REDBLUE, // same for individual (re and blue) color components
                                                     DEBAYER_GAMMA, // power function applied to the amplitudes before generating spectral masks
                                                        DEBAYER_RZ, // for green mask - rays start radius from center, relative to distance to the first alias
                                                        DEBAYER_RA, // for green mask - rays radius around aliases, relative to distance to the first alias
                                                     DEBAYER_SIGMA, // for green mask - reduce value of the far amplitudes, relative to distance to the first alias
                                                     DEBAYER_DECAY, // for green mask - exponential decay of the ray value, relative to distance to the first alias
                                              DEBAYER_FARTHEST_MAX, // fartherst absolute maximum on a ray to count
                                              DEBAYER_RADIUS_POWER, // divide ray values by the radius to this power
                                               DEBAYER_MAINTOALIAS, // relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
                                                 DEBAYER_MASK_BLUR, // for both masks  sigma for gaussian blur of the binary masks (<0 -do not use "scissors")
                                                  DEBAYER_LO_GREEN, //  combine alias-reject "scissors" with lopass filter for greens
                                              DEBAYER_LO_POSTGREEN, // combine alias-reject "scissors" with lopass filter for greens applied to red/blue
                                                DEBAYER_LO_REDBLUE, // combine alias-reject "scissors" with lopass filter for greens applied to red/blue
                                                     UPDATE_STATUS);// update status info
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
      ImageStack convolvedStack= convolveStackWithKarnelStack (imp_src.getStack(),  // stack with 3 colors/slices with the image
                                                           convolutionKernelStack, // stack with 3 colors/slices convolution kernels
                                                                CONVOLVE_FFT_SIZE, // 128 - fft size, kernel size should be size/2 
                                                                    UPDATE_STATUS); // update status info

      ImagePlus imp_convolvedStack = new ImagePlus(imp_src.getTitle()+"-convolved", convolvedStack);
      imp_convolvedStack.getProcessor().resetMinAndMax();
      imp_convolvedStack.show();
      return;

/* ======================================================================== */
    } else if (label.equals("Process Kernels")) {
      int loop_debug_level=1;
      if (!showRPSFDialog()) return;
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      ImagePlus imp_kernels = WindowManager.getCurrentImage();
      if (imp_kernels==null){
        IJ.showMessage("Error","There is no image to process");
        return;
      }
      ImageProcessor ip_kernels=imp_kernels.getProcessor();
      if (imp_kernels.getStackSize()>1) {
        PSFKernelMap=splitSquareKernelsFromStack(imp_kernels.getStack(), // Image stack, each slice consists of square kernels of one channel
                                                        PSF_kernel_size, // size of each kernel (should be square)
                                         channelNumbers(colorsToCorrect));  // four-element array of color channel numbers (usually {1,2,5,-1})
      } else {
        PSFKernelMap= splitSquareKernelsFromOne((float[]) ip_kernels.getPixels(), // flattened kernels, same as produced by mergeSquareKernelsToOne()
                                                                 PSF_kernel_size, // size of each kernel (should be square)
                                                          imp_kernels.getWidth(), // width of the flattened kernels, in pixels (should be multiple of 2*size)
                                                 channelNumbers(colorsToCorrect));  // four-element array of color channel numbers (usually {1,2,5,-1})
      }

/*    float [] kernels_pixels;
      kernels_pixels=(float[])ip_kernels.getPixels();    
      PSFKernelMap= splitSquareKernelsFromOne(kernels_pixels, // flattened kernels, same as produced by mergeSquareKernelsToOne()
                                             PSF_kernel_size, // size of each kernel (should be square)
                                      imp_kernels.getWidth(), // width of the flattened kernels, in pixels (should be multiple of 2*size)
                             channelNumbers(colorsToCorrect));  // four-element array of color channel numbers (usually {1,2,5,-1}) */

/* testing kernel interpolation  */
      if (INTERPOLATE_SUBDIV>1) { // interpolate in-place
        DEBUG_LEVEL=loop_debug_level;
        PSFKernelMap= interpolateKernels (PSFKernelMap,  // 2d array of component kernels (some components may be null)
                                    INTERPOLATE_SUBDIV, // number of subdivisions (i.e. subdiv=2 will insert interpolated kernel between each 2)
                                         UPDATE_STATUS);  // update status info
        DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      }

/* just testing - convert kernels to 3 color merged ones */
      if (SHOW_PSF ) {
        SDFA_instance.showImageStack( mergeKernelsToStack(PSFKernelMap), imp_kernels.getTitle()+"-PSF"+((INTERPOLATE_SUBDIV>1)?"-INTERPOLATED":""));
      }
/* Calulate MTF for each cell (may be moved to other command) */
      if (SHOW_MTF) {
        double [][][][] MTFs=calculateMTFfromPSF(PSFKernelMap, // 2-d array of direct psf kernels
                                             RPSF_kernel_size); // size (side of square) of reverse PSF kernel
        SDFA_instance.showImageStack( mergeKernelsToStack(MTFs), imp_kernels.getTitle()+"-MTF");
     }

     if (!SHOW_INVERTED &&!SHOW_FILTERED &&!SHOW_REDUCED_ALIASES &&!SHOW_GAUSSIANS &&!SHOW_REDUCED_ALIASES_GAUSSIAN) return; // nothing to do


/* now reverse each kernel in the array */
     DEBUG_LEVEL=loop_debug_level;
     rPSFKernelMap = reversePSFKernels(PSFKernelMap, // 2-d array of direct psf kernels
                                   RPSF_kernel_size, // size (side of square) of reverse PSF kernel
                                   OTF_deconvInvert,    // deconvInvert
                                  OTF_cutoff_energy,  // OTF_cutoff_energy
                                  OTF_ellipse_scale,  // ellipse mask size relative to the cluster
                                  OTF_ellipse_gauss,  // use Gauss instead of Hamming for ellipse mask
                                  PSF_cutoff_energy,  // OTF_cutoff_energy
                                  PSF_ellipse_scale,  // ellipse mask size relative to the cluster
                            RPSF_min_mask_threshold,  // zero output element if elliptical Gauss mask is below this threshold
                                      UPDATE_STATUS);   // update status info
     DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
 /* merge 2-d array of kernels into large 2d-array of pixels (i.e. to be shown with showBayer()) */

//     int mergedRWidth= rPSFKernelMap[0].length*RPSF_kernel_size;
//     int mergedRHeight=rPSFKernelMap.length*RPSF_kernel_size;
//     double [][] mergedRKernelsThree;
//     double [] mergedRKernelsOne;

     if (SHOW_INVERTED) {
//       if (SHOW_INDIVIDUAL) SDFA_instance.showArrays(mergeSquareKernelsToThree(rPSFKernelMap), mergedRWidth, mergedRHeight, imp_kernels.getTitle()+"-rPSF_KERNEL");
//       if (SHOW_COMPOSITE)  SDFA_instance.showArrays(mergeSquareKernelsToOne(rPSFKernelMap), 2*mergedRWidth, 2*mergedRHeight, imp_kernels.getTitle()+"-rPSF_KERNEL");
        SDFA_instance.showImageStack( mergeKernelsToStack(rPSFKernelMap), imp_kernels.getTitle()+"-rPSF_KERNEL");
     }
     if (!SHOW_FILTERED &&!SHOW_REDUCED_ALIASES &&!SHOW_GAUSSIANS &&!SHOW_REDUCED_ALIASES_GAUSSIAN) return; // nothing more to do

/* Blur deconvolution kernels with either just Gaussian or variable sigma (farther from center - more the blurring) */
     DEBUG_LEVEL=loop_debug_level;
/* TODO: use direct kernels (mirrored) for the centers*/
     variableBlurDeconvolutionKernels (rPSFKernelMap,  // 2-d array of reversed psf kernels
                                       DECONCV_BLUR_INDIVIDUAL,  // sigma for individual bayer components (1/4 pixels))
                                       DECONCV_BLUR_DIAGONAL,  // sigma for diagonal greens (square array rotated 45 degrees) (1/2 pixels))
                                       DECONCV_BLUR_CHECKER,  // sigma for checkerboard  greens (normal array, half pixels 0)
                                       RPSF_var_sigma_scale,  // scale sigma in the center when using variable sigma
                                       RPSF_sigma_to_radius,  // sigma-to-radius ratio (0.0 to disable variable blur)
                                       UPDATE_STATUS);  // update status info
     DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
     if (SHOW_FILTERED) {
//       if (SHOW_INDIVIDUAL) SDFA_instance.showArrays(mergeSquareKernelsToThree(rPSFKernelMap),   mergedRWidth,   mergedRHeight, imp_kernels.getTitle()+"-rPSF_KERNEL-BLURRED");
//       if (SHOW_COMPOSITE)  SDFA_instance.showArrays(mergeSquareKernelsToOne(  rPSFKernelMap), 2*mergedRWidth, 2*mergedRHeight, imp_kernels.getTitle()+"-rPSF_KERNEL-BLURRED");
       SDFA_instance.showImageStack( mergeKernelsToStack(rPSFKernelMap), imp_kernels.getTitle()+"-rPSF_KERNEL-BLURRED");
     }
     if (!SHOW_REDUCED_ALIASES &&!SHOW_GAUSSIANS &&!SHOW_REDUCED_ALIASES_GAUSSIAN) return; // nothing more to do

/* Create gaussian kerels with the same centers as the deconvolution ones (to re-create original image with just lateral chromatic aberration corrected use in smooth areas to reduce noise */
     DEBUG_LEVEL=loop_debug_level;
//     gaussianKernelMap= generateGaussianKernels (rPSFKernelMap,  // 2-d array of reversed psf kernels
//                                       DECONCV_BLUR_INDIVIDUAL,  // sigma for individual bayer components (1/4 pixels))
//                                         DECONCV_BLUR_DIAGONAL,  // sigma for diagonal greens (square array rotated 45 degrees) (1/2 pixels))
//                                          DECONCV_BLUR_CHECKER,  // sigma for checkerboard  greens (normal array, half pixels 0)
//                                                 UPDATE_STATUS);  // update status info


     gaussianKernelMap= generateGaussianKernelsFromDirect (PSFKernelMap,  // 2-d array of direct psf kernels
                                                       RPSF_kernel_size, // size (side of square) of reverse PSF kernel
                                                DECONCV_BLUR_INDIVIDUAL,  // sigma for individual bayer components (1/4 pixels))
                                                  DECONCV_BLUR_DIAGONAL,  // sigma for diagonal greens (square array rotated 45 degrees) (1/2 pixels))
                                                   DECONCV_BLUR_CHECKER,  // sigma for checkerboard  greens (normal array, half pixels 0)
                                                          UPDATE_STATUS);  // update status info

     DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
     if (SHOW_GAUSSIANS) {
//        if (SHOW_INDIVIDUAL) SDFA_instance.showArrays(mergeSquareKernelsToThree(gaussianKernelMap),   mergedRWidth,   mergedRHeight, imp_kernels.getTitle()+"-Gaussian");
//        if (SHOW_COMPOSITE)  SDFA_instance.showArrays(mergeSquareKernelsToOne(  gaussianKernelMap), 2*mergedRWidth, 2*mergedRHeight, imp_kernels.getTitle()+"-Gaussian");
        SDFA_instance.showImageStack( mergeKernelsToStack(gaussianKernelMap), imp_kernels.getTitle()+"-Gaussian");
     }
/* gaussian can also be with/without alias rejection */
     if (!SHOW_REDUCED_ALIASES &&!SHOW_REDUCED_ALIASES_GAUSSIAN) return; // nothing more to do

     if (SHOW_REDUCED_ALIASES) {
       DEBUG_LEVEL=loop_debug_level;
       aliasesRejectInKernels (rPSFKernelMap,  // 2-d array of reversed psf kernels
                 DECONCV_ALIASREJ_INDIVIDUAL,  // sigma for individual bayer components (1/4 pixels))
                   DECONCV_ALIASREJ_DIAGONAL,  // sigma for diagonal greens (square array rotated 45 degrees) (1/2 pixels))
                    DECONCV_ALIASREJ_CHECKER,  // sigma for checkerboard  greens (normal array, half pixels 0)
                                PSF_subpixel,  // measured array is sampled at 1/oversample frequency than model
                               UPDATE_STATUS);  // update status info
       DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
//       if (SHOW_INDIVIDUAL) SDFA_instance.showArrays(mergeSquareKernelsToThree(rPSFKernelMap),   mergedRWidth,   mergedRHeight, imp_kernels.getTitle()+"-rPSF_KERNEL-REJ");
//       if (SHOW_COMPOSITE)  SDFA_instance.showArrays(mergeSquareKernelsToOne(  rPSFKernelMap), 2*mergedRWidth, 2*mergedRHeight, imp_kernels.getTitle()+"-rPSF_KERNEL-REJ");
       SDFA_instance.showImageStack( mergeKernelsToStack(rPSFKernelMap), imp_kernels.getTitle()+"-rPSF_KERNEL-REJ");
     }

     if (SHOW_REDUCED_ALIASES_GAUSSIAN) {
       DEBUG_LEVEL=loop_debug_level;
       aliasesRejectInKernels(gaussianKernelMap,  // 2-d array of reversed psf kernels
                    DECONCV_ALIASREJ_INDIVIDUAL,  // sigma for individual bayer components (1/4 pixels))
                      DECONCV_ALIASREJ_DIAGONAL,  // sigma for diagonal greens (square array rotated 45 degrees) (1/2 pixels))
                       DECONCV_ALIASREJ_CHECKER,  // sigma for checkerboard  greens (normal array, half pixels 0)
                                   PSF_subpixel,  // measured array is sampled at 1/oversample frequency than model
                                  UPDATE_STATUS);  // update status info
       DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
//       if (SHOW_INDIVIDUAL) SDFA_instance.showArrays(mergeSquareKernelsToThree(gaussianKernelMap),   mergedRWidth,   mergedRHeight, imp_kernels.getTitle()+"-Gaussian-REJ");
//       if (SHOW_COMPOSITE)  SDFA_instance.showArrays(mergeSquareKernelsToOne(  gaussianKernelMap), 2*mergedRWidth, 2*mergedRHeight, imp_kernels.getTitle()+"-Gaussian-REJ");
       SDFA_instance.showImageStack( mergeKernelsToStack(gaussianKernelMap), imp_kernels.getTitle()+"-Gaussian-REJ");
     }
     return;

/*
 public static boolean SHOW_PSF=             true; // show combined PSF kernels (per-color and/or composite, as defined in SHOW_INDIVIDUAL, SHOW_COMPOSITE)
 public static boolean SHOW_MTF=             true; // calculate/show MTF (see notes to SHOW_PSF)
 public static boolean SHOW_INVERTED=        true; // show inverted kernels (unfiltered), same notes
 public static boolean SHOW_FILTERED=        true; // filter and show inverted kernels
 public static boolean SHOW_REDUCED_ALIASES= true; // calculate kernels with suppressed sampling aliases patterns
 public static boolean SHOW_GAUSSIANS=       true; // create gaussian kernels with the same centers as inverted ones (low noise, use in low details areas)
 public static boolean SHOW_INDIVIDUAL=      true; // for each of the kernels above - show per-color images
 public static boolean SHOW_COMPOSITE=       true; // for each of the kernels above - show single composite image

*/
    } else if (label.equals("Convolve Image")) {
      if (rPSFKernelMap==null) {
        IJ.showMessage("Error","No convolution kernel array defined.\nProcess canceled");
        return;
      }
      if (!showDeconvDialog()) return;
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;

//      imp_sel= selectImage(2*FFTSize, false);
      imp_sel= selectImage(0, true);
      Roi roi_sel= imp_sel.getRoi();
      if (roi_sel==null){
        imp_sel.setRoi(0, 0, imp_sel.getWidth(), imp_sel.getHeight());
        roi_sel= imp_sel.getRoi();
      }
      int kernelStep=FFTOverlap*PSF_subpixel/2; // 64
      int kernelMargins=FFTSize/FFTOverlap/2; // 4
      int slidingFFTSize=128;
      int outWidth= imp_sel.getWidth()* PSF_subpixel/2;
      int compMask=0;
      float [][]bayerWeights=null;
      if (CALC_BAYER_WEIGHTS) bayerWeights=new float[2][];
      for (i=0;i<colorsToCorrect.length;i++) if (colorsToCorrect[i]) compMask|=(1<<i);
      float [][] convolved= convolveImageWithKernels(imp_sel, // 32-bit (float) linearized image with Bayer mosaic
                                         roi_sel.getBounds(), // area to process (if not all the image) in pixels
                                              equalizeGreens, // equalize two green gains (currently - independently for each tile)
                                               rPSFKernelMap, // 2d array of convolution kernels 
                                                PSF_subpixel, // measured kernels are sampled at subpixel higher than Bayer (su2.000bpixel/2 - than sensor) //4
                                                  kernelStep, // distance between kernel samples, in sensor pixels // 64 // 32
                                               kernelMargins, // number of kernel steps from the image 0,0 points to the center of the (0,0) kernel // 4
                                              slidingFFTSize, // size of sliding FFT (twice the kernel size) here:128
                                                    compMask, // bitmask of color components to process
                                          MASK_BAYER_ALIASES, // filter Bayer
/* Parameters related to reducing mosaic alias artifacts, can be removed if that filtering will be preformed before aberration correction.
     That just requires that lateral chtomatic aberration is small relative to FFT size, so spectrums amplitudes are similar */
                                               DEBAYER_GAMMA, // power function applied to the amplitudes before generating spectral masks
                                                  DEBAYER_RZ, // for green mask - rays start radius from center, relative to distance to the first alias
                                                  DEBAYER_RA, // for green mask - rays radius around aliases, relative to distance to the first alias
                                               DEBAYER_SIGMA, // for green mask - reduce value of the far amplitudes, relative to distance to the first alias
                                               DEBAYER_DECAY, // for green mask - exponential decay of the ray value, relative to distance to the first alias
                                        DEBAYER_FARTHEST_MAX, // for green mask - rays will be extyended from maf back to the center, max up to this distance (relative)
                                        DEBAYER_RADIUS_POWER, // for green mask -  divide ray values by the radius to this power
                                         DEBAYER_MAINTOALIAS,// relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
                                           DEBAYER_MASK_BLUR, // for both masks  sigma for gaussian blur of the binary masks
                                            DEBAYER_LO_GREEN, //  combine alias-reject "scissors" with lopass filter for greens
                                        DEBAYER_LO_POSTGREEN, // combine alias-reject "scissors" with lopass filter for greens applied to red/blue
                                          DEBAYER_LO_REDBLUE, // combine alias-reject "scissors" with lopass filter for greens applied to red/blue
/* end of demosaic parameters */
                                                bayerWeights,
                                              UPDATE_STATUS);  // update status info
//      SDFA_instance.showArrays(convolved,outWidth,outHeight,imp_sel.getTitle()+"-CONVOLVED");

//      ImageStack convolvedStack=new ImageStack(outWidth,outHeight,3);
      ImageStack convolvedStack= combineRGBCorrectBayer (convolved[1],
                                                        convolved[5],
                                                        convolved[2],
                                                            outWidth,   // image width
                                                                   0); // corrected already in convolveImageWithKernels()  PSF_subpixel/2);  // half Bayer period (GR/BG)
      if (bayerWeights!=null) {
        convolvedStack.addSlice("R-weight",bayerWeights[0]);
        convolvedStack.addSlice("B-weight",bayerWeights[1]);
      }
      ImagePlus imp_stack_convolved = new ImagePlus(imp_sel.getTitle()+"convolved-rgb", convolvedStack);
      imp_stack_convolved.getProcessor().resetMinAndMax();
      imp_stack_convolved.show();
     return;


/* ======================================================================== */
    } else if (label.equals("Colors")) {

      imp_convolved=WindowManager.getCurrentImage();
      if ((imp_convolved==null) || (imp_convolved.getStackSize()<3)) {
        IJ.showMessage("Error","Please select image stack of 5 (float) slices (r,g,b) and 2 weight slices");
        return;
      }
      if (!showColorProcessDialog()) return;



      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      if (DEBUG_LEVEL>1) System.out.println ( "Read image "+imp_convolved.getTitle()+" as the image convolved with the reversed PSF kernels");
      ImageStack stack_convolved=imp_convolved.getStack();

/* Set stack sequence to r-g-b */
// public static String [] stackColorNames={"red","green","blue"};

/* find number of the green channel - should be called "green", if none - use last */
     int [] rgbNumbers= {0,0,0};
     for (j=0;j<3;j++) {
       for (i=1;i<=3;i++) if (stack_convolved.getSliceLabel(i).equals(stackColorNames[j])){
         rgbNumbers[j]=i;
         break;
       }
     }
     if (DEBUG_LEVEL>1) {
       System.out.println ( "Input file color slice numbers:");
       System.out.println ( "  red -   slice "+((rgbNumbers[0]>0)?rgbNumbers[0]:"missing"));
       System.out.println ( "  green - slice "+((rgbNumbers[1]>0)?rgbNumbers[1]:"missing"));
       System.out.println ( "  blue -  slice "+((rgbNumbers[2]>0)?rgbNumbers[2]:"missing"));
     }
     for (i=0;i<3;i++) if (rgbNumbers[i]<=0) {
       System.out.println ( stackColorNames[i]+ "  slice is missing in the input file. Please check slice names");
       return;
     }
     while ((rgbNumbers[0]!=1) || (rgbNumbers[1]!=2) ||(rgbNumbers[2]!=3)) {
       if      (rgbNumbers[0]==1) swapStackSlices(stack_convolved,2,3);
       else if (rgbNumbers[2]==3) swapStackSlices(stack_convolved,1,2);
       else                       swapStackSlices(stack_convolved,1,3);
       for (j=0;j<3;j++) {
         for (i=1;i<=3;i++) if (stack_convolved.getSliceLabel(i).equals(stackColorNames[j])){
           rgbNumbers[j]=i;
           break;
         }
       }
     }




      processColorsWeights(stack_convolved,
                               255.0/PSF_subpixel/PSF_subpixel,
                               BALANCE_RED,     // 1.0; // manual color balance, gain 1.0 matches 0.0.255.0 range of the unput Bayer data
                               BALANCE_BLUE,    // 1.0;
                               GAIN_GREEN,      // 1.0;
                               WEIGHT_SCALE_R,  // 1.0; // additional correction for the weights of colors for different sub-pixels in a Bayer cell
                               WEIGHT_SCALE_B,  // 1.0; // WEIGHT_SCALE_G=1.0
                               COLOR_SIGMA,     // 2.0; // Gaussian sigma to low-pass color components when calculating "smooth" color
                               YCbCr_Gamma,     // 0.53;
                               YCbCr_minLin,    //    0.003;
                               YCbCr_Kr,        // 0.299;
                               YCbCr_Kb,        // 0.114;
                               COLOR_SATURATION, // Pr2.0
                               COLOR_SATURATION // Pb2.0
                           );

//USE_FIRST_Y
     YPrPbToRGB(stack_convolved,
                       YCbCr_Kr,        // 0.299;
                       YCbCr_Kb,        // 0.114;
                USE_FIRST_Y?9:8,        //  int sliceY,
                              6, // int slicePr,
                              7// int slicePb
                           );

      imp_convolved.getProcessor().resetMinAndMax();
      imp_convolved.updateAndDraw();

      SDFA_instance.showImageStackThree(stack_convolved, imp_convolved.getTitle()+"restored_rgb");



//BALANCE_RED
     return;

/* ======================================================================== */
    } else if (label.equals("Test Debayer")) {

      ImagePlus imp_debayer=WindowManager.getCurrentImage();
      if (imp_debayer==null) {
        IJ.showMessage("Error","Please select image stack of 5 (float) slices (r,g,b) and 2 weight slices");
        return;
      }
      if (!showDeBayerDialog()) return;
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      float [] fpixels_debayer= (float[]) imp_debayer.getProcessor().getPixels();
      double [][] pixels_debayer=new double [3][fpixels_debayer.length];
      for (i=0;i<fpixels_debayer.length;i++) pixels_debayer[1][i]=fpixels_debayer[i];
      pixels_debayer[0]=pixels_debayer[1].clone();
      pixels_debayer[2]=pixels_debayer[1].clone();
      if (DEBAYER_TEST_MASKSPLIT) {
        pixels_debayer= normalizeAndWindow (pixels_debayer, initHamming((int)Math.sqrt(pixels_debayer[1].length)),false);
        pixels_debayer=extendFFTInputTo (pixels_debayer, DEBAYER_FFT_SIZE);
        if (DEBUG_LEVEL>2)   SDFA_instance.showArrays(pixels_debayer,  "DT-all");
        for (i=0;i<DEBAYER_FFT_SIZE;i++) for (j=0;j<DEBAYER_FFT_SIZE;j++) {
          pixels_debayer[1][i*DEBAYER_FFT_SIZE+j]*= ((((i % PSF_subpixel)==0) && ((j % PSF_subpixel)==0 )) || (((i % PSF_subpixel)==(PSF_subpixel/2)) && ((j % PSF_subpixel)==(PSF_subpixel/2) )))?1.0:0.0;
          pixels_debayer[0][i*DEBAYER_FFT_SIZE+j]*=  (((i % PSF_subpixel)==0) && ((j % PSF_subpixel)==(PSF_subpixel/2) ))?1.0:0.0;
          pixels_debayer[2][i*DEBAYER_FFT_SIZE+j]*=  (((i % PSF_subpixel)==(PSF_subpixel/2)) && ((j % PSF_subpixel)==0 ))?1.0:0.0;
        }
        if (DEBUG_LEVEL>1)   SDFA_instance.showArrays(pixels_debayer,  "DT-pixels");
      }
/* swap quadreants and perform FHT */
      double [][] amps=null;
      if (DEBUG_LEVEL>1)  amps=new double[3][];

      for (i=0;i<3;i++) {
        fht_instance.swapQuadrants(pixels_debayer[i]);
        fht_instance.transform(pixels_debayer[i]);
        if (amps!=null) amps[i]=fht_instance.calculateAmplitude(pixels_debayer[i]);

      }
   PolarSpectrums pol_instace=new PolarSpectrums(
                 DEBAYER_FFT_SIZE, // size of the square array, centar is at size/2, size/2, only top half+line will be used
                 Math.PI, //2*Math.PI, // i.e. Math.PI, 2*Math.PI
                 DEBAYER_FFT_SIZE/2-2, // width of the polar array - should be <= size/2-2
                 0.5, //0.75, //2.0, //0.5, // maximal step in pixels on the maxRadius for 1 angular step (i.e. 0.5)
                 4 );// angular symmetry - 0- none,1 - pi corresponds to integer, 2 - pi/2 corresponds to integer, n - pi/n corresponds to integer angular step 
      double [][][]  lopass=  createAliasFilters (DEBAYER_WIDTH_GREEN, // realtive width of the low-pass filter for greens (1.0 goes to zero at exactrly the closest alias )
                                                DEBAYER_WIDTH_REDBLUE, // same for individual (re and blue) color components
                                                     DEBAYER_FFT_SIZE, // side of the square
                                                         PSF_subpixel); // should be 4 now

      double [][] both_masks= aliasScissors(pixels_debayer[1], // fht array for green, will be masked in-place
                                                 PSF_subpixel, // measured kernels are sampled at subpixel higher than Bayer (subpixel/2 - than sensor) //4
                                            DEBAYER_THRESHOLD, // no high frequencies - use default uniform filter
                                                DEBAYER_GAMMA, // power function applied to the amplitudes before generating spectral masks
                                                   DEBAYER_RZ, // for green mask - rays start radius from center, relative to distance to the first alias
                                                   DEBAYER_RA, // for green mask - rays radius around aliases, relative to distance to the first alias
                                                DEBAYER_SIGMA, // for green mask - reduce value of the far amplitudes, relative to distance to the first alias
                                                DEBAYER_DECAY, // for green mask - exponential decay of the ray value, relative to distance to the first alias
                                         DEBAYER_FARTHEST_MAX, // fartherst absolute maximum on a ray to count
                                         DEBAYER_RADIUS_POWER, // divide ray values by the radius to this power
                                          DEBAYER_MAINTOALIAS, // relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
                                            DEBAYER_MASK_BLUR, // for both masks  sigma for gaussian blur of the binary masks (<0 -do not use "scissors")
                                             DEBAYER_LO_GREEN, //  combine alias-reject "scissors" with lopass filter for greens
                                         DEBAYER_LO_POSTGREEN, // combine alias-reject "scissors" with lopass filter for greens applied to red/blue
                                           DEBAYER_LO_REDBLUE, // combine alias-reject "scissors" with lopass filter for greens applied to red/blue
                           (DEBAYER_DECAY<0)?pol_instace:null,
                                                       lopass,
                                                  DEBUG_LEVEL); // internal debug level
      double [] green_mask=   both_masks[0];
      double [] red_blue_mask=both_masks[1];
//      if (DEBAYER_TEST_FROMIMAGE) {  // - DEBAYER_TEST_FROMIMAGE - OBSOLETE
        
        pixels_debayer[1]=fht_instance.multiply(pixels_debayer[1],green_mask,false);
        pixels_debayer[0]=fht_instance.multiply(pixels_debayer[0],red_blue_mask,false);
        pixels_debayer[2]=fht_instance.multiply(pixels_debayer[2],red_blue_mask,false);

        fht_instance.inverseTransform(pixels_debayer[1]);
        fht_instance.swapQuadrants(pixels_debayer[1]);

        fht_instance.inverseTransform(pixels_debayer[0]);
        fht_instance.swapQuadrants(pixels_debayer[0]);

        fht_instance.inverseTransform(pixels_debayer[2]);
        fht_instance.swapQuadrants(pixels_debayer[2]);

       SDFA_instance.showArrays(pixels_debayer,  "filtered");

/* Swap quadrants in masks and display them */
      if (DEBUG_LEVEL>1){
        fht_instance.swapQuadrants(green_mask);
        fht_instance.swapQuadrants(red_blue_mask);
        SDFA_instance.showArrays(green_mask,     "G-mask");
        SDFA_instance.showArrays(red_blue_mask,  "RB-mask");
        if (amps!=null) {
/**normailze amplitudes, apply gamma */
          double dmax=0.0;
          for (i=0;i<amps.length;i++) {
            for (j=0;j<amps[i].length;j++) if (amps[i][j]>dmax) dmax=amps[i][j];
            dmax=1.0/dmax;
            for (j=0;j<amps[i].length;j++) amps[i][j]= Math.pow(amps[i][j]*dmax,DEBAYER_GAMMA);
          }
          SDFA_instance.showArrays(amps,  "A");
          for(i=0;i<amps[1].length;i++){
             amps[1][i]*=green_mask[i];
             amps[0][i]*=red_blue_mask[i];
             amps[2][i]*=red_blue_mask[i];
          }
          SDFA_instance.showArrays(amps,  "A-masked");
       }
      }
      return;

/* ======================================================================== */

    }
  }

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
 




  public double kernelShotNoiseFactor (double[] kernel) {
    double s=0,s2=0;
    for (int i=0; i<kernel.length;i++) {
      s+=kernel[i];
      s2+=kernel[i]*kernel[i];
    }
    return Math.sqrt(s2)/Math.abs(s);
  }
/*

*/
// temporary using float implementation in ImageJ - re-write to directly use double [] arrays
  public void  blurDouble(double[] pixels,
                                int width,
                            double sigmaX,
                            double sigmaY,
                         double precision) {
//  public void  blurFloat(red_blue_mask, DEBAYER_MASK_BLUR, DEBAYER_MASK_BLUR, 0.01);
      int i;
      int height = pixels.length/width;
      float [] fpixels=new float [pixels.length];
      for (i=0;i<pixels.length;i++) fpixels[i]= (float) pixels[i];
      FloatProcessor fp = new FloatProcessor(width, height, fpixels, null);
      GaussianBlur gb = new GaussianBlur();
      gb.blurFloat(fp, sigmaX, sigmaY, precision);
      for (i=0;i<pixels.length;i++) pixels[i]=fpixels[i];
 }

  public double [] calcRedBlueAliasMask (double [] green_amp,
                                        double [] green_mask, // may be null if amp_pixels is already masked
                                                int subpixel,
                                          double mainToAlias){// relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)

    int length=green_amp.length;
    int size = (int) Math.sqrt(length);
    int hsize=size/2;
    int i,j,index,index_back,x,y;
    int aliasX=size/subpixel;
    double [] amp=green_amp.clone();
    if (green_mask!=null) for (i=0;i<amp.length;i++) amp[i]*=green_mask[i];
    double [] mask= new double [length];
    for (i=0;i<length;i++) mask[i]=0.0;
/* Combine into mask by comparing pixels[] from the zero and 7 aliases */
    double d;
    int nAlias;
    int [][] aliasMap={{-1,-1},{-1,0},{-1,1},
                       { 0,-1},       { 0,1},
                       { 1,-1},{ 1,0},{ 1,1}};
    for (i=0;i<=hsize;i++) for (j=0;j<size;j++) {
      index=i*size+j;
      index_back=((size-i) % size) * size + ((size-j) % size);
      d=amp[index]*mainToAlias;
      if (d>0.0) {
        mask[index]=1.0;
        mask[index_back]=1.0;
//        isGreater=true;
        for(nAlias=0;nAlias<aliasMap.length; nAlias++) {
          y=i-aliasX*aliasMap[nAlias][0];
          x=j-aliasX*aliasMap[nAlias][1];
          if (y>hsize) {
            y=size-y;
            x=size-x;
          }
          if ((x>=0)&& (y>=0) && (amp[y*size+x]>d)) {
            mask[index]=0.0;
            mask[index_back]=0.0;
            break;
          }
        }
      }
    }
    return mask;
  }


  public double [] calcGreensAliasMask (double [] amp_pixels, // normalized amplitude spectrum, (0,0) in the center
                                                 int subpixel,
                                       double radius_zero_rel,
                                      double radius_alias_rel,
                                             double sigma_rel,
                                             double decay_rel,
                                      double farterst_max_rel,
                                          double radius_power, //divide by radius to this power
                                          double mainToAlias){// relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)

    int length=amp_pixels.length;
    int size = (int) Math.sqrt(length);
    int hsize=size/2;
    double [] mask=  new double [length];
    double [] pixels=new double [length];
    double [] mamp_pixels= amp_pixels.clone();
    int []    pixels_num=new int [length];
    int i,j;
    for (i=0;i<length;i++) pixels[i]=0.0;
    for (i=0;i<length;i++) mask[i]=0.0;
    for (i=0;i<length;i++) pixels_num[i]=0;
    int aliasX=size/subpixel;
    double aliasR=(Math.sqrt(2)*aliasX);
    int maxR=(int) aliasR+1;
    double radius_zero = aliasR*radius_zero_rel;
    double radius_alias= aliasR*radius_alias_rel;
    double sigma=        aliasR*sigma_rel;
    double decay=        aliasR*decay_rel;
    double farterst_max= aliasR*farterst_max_rel;
    double [] gaussian=new double[hsize+1];
    int x,y,absX,xy0,xy1;
    double r2=maxR*maxR;
    double [] dr     = new double[maxR+1];
    int    [] indices= new int [maxR+1];
    double [] dr_scale=new double[maxR+1];
    if (radius_zero> aliasX) radius_zero= aliasX;
    if (radius_alias>aliasX) radius_alias=aliasX;
    double radius_alias2=radius_alias*radius_alias;
    double x2,y2;
    double k=-0.5/(sigma*sigma);
    dr_scale[0]=1.0; 
    for (i=1;i<dr_scale.length;i++)  dr_scale[i]=1.0/((radius_power!=0)?Math.pow(i,radius_power):1.0);
    for (i=0;i<gaussian.length;i++) gaussian[i]=Math.exp(k*i*i);
    for (i=0;i<=hsize;i++) {
      y=hsize-i;
      y2=y*y;
      for (j=0;j<size;j++) {
        x=j-hsize;
        absX=Math.abs(x);
        x2=x*x;
        if ((x2+y2)>r2) mamp_pixels[i*size+j]=-1.0; // mask out
        else {
          mamp_pixels[i*size+j]*=gaussian[y]*gaussian[absX];
        }
      }
    }
    xy0=(int) (hsize-aliasX-radius_alias-0.5);
    if (xy0<1) xy0=1;
    xy1=(int) (hsize-aliasX+radius_alias+0.5);
    if (xy1>hsize) xy1=hsize;

    if (DEBUG_LEVEL>3) {
      System.out.println ( " calcGreensAliasMask() xy0="+xy0+" xy1="+xy1);
    }

    for (i=xy0;i<=xy1;i++) {
      y=i-(hsize-aliasX);
      y2=y*y;
      for (j=xy0;j<=xy1;j++) {
        x=j-(hsize-aliasX);
        x2=x*x;
        if ((x2+y2)<=radius_alias2) {
          mamp_pixels[    i*size+j]=-1.0; // mask out
          mamp_pixels[(i+1)*size-j]=-1.0; // mask out
        }
      }
    }
    if (DEBUG_LEVEL>2)   SDFA_instance.showArrays(mamp_pixels,  "mamp"); // only top half+1 will be used
    int ia,r,index,index_back;

//    int maxA=(int) Math.round(Math.PI*maxR);
    int maxA=(int) Math.round(2.0*Math.PI*maxR);
    double cos,sin;
    double d,p;
    double dk=Math.exp(-1/decay);// fraction when radius increments by 1
    double max;
    int    rmax, rfar;
    double a;
    double ka=Math.PI/maxA;
    double dx,dy;
    int rMaxInterpolate= (int) aliasR; // above that no interpolation. Maybe just set to some small number?
    double p_tl, p_tr, p_br; // points fro interpolation near p
    int index01,index10,index11;
    if (DEBUG_LEVEL>3) {
          System.out.println ( " calcGreensAliasMask() maxA="+maxA+" maxR="+maxR+" subpixel="+ subpixel+" radius_zero="+radius_zero+" radius_alias="+radius_alias+" sigma="+sigma+" decay="+ decay);
    }
    
    for (ia=0;ia<=maxA;ia++) {
      a=ka*ia;
      cos=Math.cos(a);
      sin=Math.sin(a);
      d=0.0;
      max=0.0;
      rmax=0;
      rfar=0;
      for (r=0;r<=maxR;r++ ) {
        dy=r*sin;
        if (dy<0) break;
        y=(int) Math.round(dy);
//        if (y<0) break;
        i=(hsize-y);
        if (i<0) break;
        dx=r*cos;
        x=(int) Math.round(dx);
        dy-=y;
        dx-=x;
        j=hsize+x;
        if (j<0) break;
        if (j>=size) break;
        index=i*size+j;
//        if ((i==hsize) && (j<hsize) && (j>0)) index=i*size+(size-j); // Not needed - this line is filled all
        p=mamp_pixels[index];
/* interpolate? otherwise high-angular frequency artifacts*/
        if (p<0) break;

        if (r<rMaxInterpolate) {
          index01=index;
          index10=index+1;
          index11=index+1;
          if (dx<0) {
            index10=index-1;
            dx=-dx;
          } else {
            index10=index+1;
          }
          if (dy<0) {
            index01=index+size;
            index11=index+size;
            dy=-dy;
            if (y==0) { // special case, crossing to the bottom (uninitialized) half - mirror around the center
              index01= size*size-index01;
              index11= size*size-index11;
            }
          } else {
            index01=index-size;
            index11=index-size;
          }

          p_tl=mamp_pixels[index01];
          if (p_tl<0) break;
          p_tr=mamp_pixels[index11];
          if (p_tr<0) break;
          p_br=mamp_pixels[index10];
          if (p_br<0) break;
//          p=(1-dy)*((1-dx)*p_tl+dx*p_tr) + dy* ((1-dx)*p+dx*p_br);
          p=dy*((1-dx)*p_tl+dx*p_tr) + (1-dy)* ((1-dx)*p+dx*p_br);
        }

        rfar=r;
        if (r<radius_zero) d=0;
        else d= (d*dk)+p;
        dr[r]=d*d;
        indices[r]=index;
        if ((r<=farterst_max) && (d>max)) {
          max=d;
          rmax=r;
        }

//        pixels[index]=d; // will write several times to the same element
//        if (max<d) max=d;
      }
      for (r=0;r<=rfar;r++ ) {
        d=dr_scale[r]*((r<rmax)?(max*max):dr[r]);
        index=indices[r];
        pixels[index]+=d; // will write several times to the same element
        pixels_num[index]++;
      }
    }
// farterst_max
//    if (DEBUG_LEVEL>3)   SDFA_instance.showArrays(pixels_debayer,  "field"); // only top half+1 will be used
    for (i=0;i<length;i++) if (pixels_num[i]>0) {
    }//pixels[i] /=pixels_num[i];
    for (i=0;i<=hsize;i++) for (j=(i==hsize)?hsize:0;j<size;j++) {
      index=i*size+j;
      index_back=((size-i) % size) * size + ((size-j) % size);
      pixels[index] = (pixels_num[index]>0)?Math.sqrt(pixels[index]/pixels_num[index]):0.0; // is it needed if we just compare values?
      pixels[index_back]=pixels[index];
    }
    if (DEBUG_LEVEL>2)   SDFA_instance.showArrays(pixels,  "field"); // only top half+1 will be used

//    System.out.println("polar_green_mask_pixels[-1]="+polar_green_mask_pixels[-1]); // make an error

// Combine into mask by comparing pixels[] from the zero and 7 aliases 
    int nAlias;
    int [][] aliasMap={{-2,-2},{-2,0},{-2,2},
                           {-1,-1},{-1,1},
                       { 0,-2},       { 0,2},
                           { 1,-1},{ 1,1},
                       { 2,-2},{ 2,0},{ 2,2}};
//    boolean isGreater;
    for (i=0;i<=hsize;i++) for (j=0;j<size;j++) {
      index=i*size+j;
      index_back=((size-i) % size) * size + ((size-j) % size);
      d=pixels[index]*mainToAlias;
      if (d>0.0) {
        mask[index]=1.0;
        mask[index_back]=1.0;
//        isGreater=true;
        for(nAlias=0;nAlias<aliasMap.length; nAlias++) {
          y=i-aliasX*aliasMap[nAlias][0];
          x=j-aliasX*aliasMap[nAlias][1];
          if (y>hsize) {
            y=size-y;
            x=size-x;
          }
          if ((x>=0)&& (y>=0) && (pixels[y*size+x]>d)) { // -3968
            mask[index]=0.0;
            mask[index_back]=0.0;
            break;
          }
        }
      }
    }
//    if (DEBUG_LEVEL>1)   SDFA_instance.showArrays(mask,size,size,   "old_mask");
    return mask;

  }

  public double [] calcRedBlueAliasMaskRays (double [] green_amp, // both halves are needed ??
                                            double [] green_mask, // may be null if amp_pixels is already masked
                                      PolarSpectrums pol_instace, // initialized instance (if null - skip rays processing)
                                               double mainToAlias,// relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
                                                     double bonus, // scale far pixels as (1.0+bonus*r/rmax)
                                                   int this_debug){// relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)

    int length=green_amp.length;
    int size = (int) Math.sqrt(length);
    int hsize=size/2;
    int subpixel=4; // hardwired - when changing it will need to change alias maps
    int aliasX=size/subpixel;
    int i,j,index,index_back,x,y;
    double [] amp=green_amp.clone();
    if (green_mask!=null) for (i=0;i<amp.length;i++) amp[i]*=green_mask[i];
    double [] mask= new double [length];
    for (i=0;i<length;i++) mask[i]=0.0;
/* Combine into mask by comparing pixels[] from the zero and 7 aliases */
    double d;
    int nAlias;
    int [][] aliasMapRedBlue={{-2,-2},{-2,-1},{-2,0},{-2,1},
                              {-1,-2},{-1,-1},{-1,0},{-1,1},
                              { 0,-2},{ 0,-1},       { 0,1},
                              { 1,-2},{ 1,-1},{ 1,0},{ 1,1}};

/*    int [][] aliasMap={{-1,-1},{-1,0},{-1,1},
                       { 0,-1},       { 0,1},
                       { 1,-1},{ 1,0},{ 1,1}};*/

/* First step - mask out all the pixels where at least one of the alias amplitude is above the main one */

    for (i=0;i<=hsize;i++) for (j=0;j<size;j++) {
      index=i*size+j;
      index_back=((size-i) % size) * size + ((size-j) % size);
      d=amp[index]*mainToAlias;
      if (d>0.0) {
        mask[index]=1.0;
        mask[index_back]=1.0;
//        isGreater=true;
        for(nAlias=0;nAlias<aliasMapRedBlue.length; nAlias++) {
          y=(i-aliasX*aliasMapRedBlue[nAlias][0]+size) % size;
          x=(j-aliasX*aliasMapRedBlue[nAlias][1]+size) % size;
          if (y>hsize) {
            y=size-y;
            x=(size-x)%size;
          }
          if (amp[y*size+x]>d) {
            mask[index]=0.0;
            mask[index_back]=0.0;
            break;
          }
        }
      }
    }
    if (pol_instace==null) return mask;
/* Now apply mask to amplitudes and use ray processing (same as with greens)*/
    for (i=0;i<amp.length;i++) amp[i]*=mask[i];
    double [] polar_amp=pol_instace.cartesianToPolar(amp);
    if (this_debug>2)   SDFA_instance.showArrays(polar_amp.clone(),pol_instace.getWidth(),pol_instace.getHeight(),  "RB-polar-amp");
    double k= bonus/pol_instace.getWidth();
    for (i=0;i<pol_instace.getHeight();i++) for (j=0;j<pol_instace.getWidth();j++) polar_amp[i*pol_instace.getWidth()+j]*=1.0+k*j;
    double [] polar_mask_pixels=pol_instace.genPolarRedBlueMask(polar_amp,0); // 0 - just 1.0/0.0, 1 - "analog"
    double [] cart_mask_pixels= pol_instace.polarToCartesian (polar_mask_pixels,size,0.0);
    if (this_debug>2) {
       SDFA_instance.showArrays(polar_amp,  pol_instace.getWidth(),pol_instace.getHeight(),     "RB-amp-bonus");
       SDFA_instance.showArrays(polar_mask_pixels,pol_instace.getWidth(),pol_instace.getHeight(), "pRBm");
       SDFA_instance.showArrays(cart_mask_pixels,size,size,   "cRBm");
    }
    if (this_debug>2) {
       double [] polar_mask_pixels1=pol_instace.genPolarRedBlueMask(polar_amp,1);
       double [] cart_mask_pixels1= pol_instace.polarToCartesian (polar_mask_pixels1,size,0.0);
       SDFA_instance.showArrays(polar_mask_pixels1,pol_instace.getWidth(),pol_instace.getHeight(), "pRBm1");
       SDFA_instance.showArrays(cart_mask_pixels1,size,size,   "cRBm1");
    }
    return cart_mask_pixels;

  }





 public double [] calcGreensAliasMaskRays (double [] amp_pixels, // normalized amplitude spectrum, (0,0) in the center
                                     PolarSpectrums pol_instace, // initialized instance
                                                   double bonus, // scale far pixels as (1.0+bonus*r/rmax)
                                                 int this_debug){// relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)

    int length=amp_pixels.length;
    int size = (int) Math.sqrt(length);
    double [] polar_amp_pixels=pol_instace.cartesianToPolar(amp_pixels);
    if (this_debug>2)   SDFA_instance.showArrays(polar_amp_pixels.clone(),pol_instace.getWidth(),pol_instace.getHeight(),  "polar-amp");

    double k= bonus/pol_instace.getWidth();
    for (int i=0;i<pol_instace.getHeight();i++) for (int j=0;j<pol_instace.getWidth();j++) polar_amp_pixels[i*pol_instace.getWidth()+j]*=1.0+k*j;
    double [] polar_green_mask_pixels=pol_instace.genPolarGreenMask(polar_amp_pixels,0); // 0 - just 1.0/0.0, 1 - "analog"
    double [] cart_green_mask_pixels= pol_instace.polarToCartesian (polar_green_mask_pixels,size,0.0);
    if (this_debug>2) {
       SDFA_instance.showArrays(polar_amp_pixels,  pol_instace.getWidth(),pol_instace.getHeight(),     "amp-bonus");
       SDFA_instance.showArrays(polar_green_mask_pixels,pol_instace.getWidth(),pol_instace.getHeight(), "pgm");
       SDFA_instance.showArrays(cart_green_mask_pixels,size,size,   "cgm");
    }

    if (this_debug>2) {
       double [] polar_green_mask_pixels1=pol_instace.genPolarGreenMask(polar_amp_pixels,1);
       double [] cart_green_mask_pixels1= pol_instace.polarToCartesian (polar_green_mask_pixels1,size,0.0);
       SDFA_instance.showArrays(polar_green_mask_pixels1,pol_instace.getWidth(),pol_instace.getHeight(), "PGM1");
       SDFA_instance.showArrays(cart_green_mask_pixels1,size,size,   "CGM1");
    }
    return cart_green_mask_pixels;

  }


 public double [] calcGreensAliasMaskRaysDebug (double [] amp_pixels, // normalized amplitude spectrum, (0,0) in the center
                                     PolarSpectrums pol_instace, // initialized instance
                                                   double bonus, // scale far pixels as (1.0+bonus*r/rmax)
                                                 int this_debug){// relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)

    int length=amp_pixels.length;
    int size = (int) Math.sqrt(length);
    double [] polar_amp_pixels=pol_instace.cartesianToPolar(amp_pixels);

    if (this_debug>2)   SDFA_instance.showArrays(polar_amp_pixels,pol_instace.getWidth(),pol_instace.getHeight(),  "polar-amp"); // only top half+1 will be used
    if (this_debug>3) {
      double [] cart_amp_pixels= pol_instace.polarToCartesian (polar_amp_pixels );
      SDFA_instance.showArrays(cart_amp_pixels,size,size,  "cart-amp"); // only top half+1 will be used
    }


/*
    if (this_debug>2) {
      double [] polar_mamp_pixels=pol_instace.cartesianToPolar(mamp_pixels);
      SDFA_instance.showArrays(polar_mamp_pixels,pol_instace.getWidth(),pol_instace.getHeight(),  "polar-mamp"); // only top half+1 will be used
    }
    if (this_debug>2) {
      double [] polar_field_pixels=pol_instace.cartesianToPolar(pixels);
      SDFA_instance.showArrays(polar_field_pixels,pol_instace.getWidth(),pol_instace.getHeight(),  "polar-field"); // only top half+1 will be used
    }
*/
// generate image with pixels equal to number of the polar pixels mapping to the same cartesian cell
    if (this_debug>2) {
      double [] same_cart_pixels= pol_instace.testMapsLengths(0);
      SDFA_instance.showArrays(same_cart_pixels,pol_instace.getWidth(),pol_instace.getHeight(),  "same_cart");

// generate image with pixels equal to number of green alias cells in the same polar grid
      double [] green_alias_pixels= pol_instace.testMapsLengths(1);
      SDFA_instance.showArrays(green_alias_pixels,pol_instace.getWidth(),pol_instace.getHeight(),  "green_alias");
      double [] cart_same_cart_pixels= pol_instace.polarToCartesian (same_cart_pixels);
      SDFA_instance.showArrays(cart_same_cart_pixels,size,size,  "cart-same_cart");
      double [] cart_green_alias_pixels= pol_instace.polarToCartesian (green_alias_pixels);
      SDFA_instance.showArrays(cart_green_alias_pixels,size,size,  "cart-green_alias");
    }
    if (this_debug>2) {
      int iaa=0;
      int irr=0;
      double [] test_green_map= pol_instace.testGreenMap(iaa,irr); // (int ia, int ir);
      if (this_debug>1)   SDFA_instance.showArrays(test_green_map,pol_instace.getWidth(),pol_instace.getHeight(), "tgm_"+iaa+"_"+irr);
      double [] cart_test_green_map= pol_instace.polarToCartesian (test_green_map);
      if (this_debug>1)   SDFA_instance.showArrays(cart_test_green_map,size,size,   "ctgm_"+iaa+"_"+irr);
    }
    if (this_debug>2) {
      double [] polar_green_mask_pixels=pol_instace.genPolarGreenMask(polar_amp_pixels,1);
      SDFA_instance.showArrays(polar_green_mask_pixels,pol_instace.getWidth(),pol_instace.getHeight(), "PGM");
      double [] cart_green_mask_pixels= pol_instace.polarToCartesian (polar_green_mask_pixels,size,0.0);
      SDFA_instance.showArrays(cart_green_mask_pixels,size,size,   "CGM");
    }


    double [] polar_green_mask_pixels0=pol_instace.genPolarGreenMask(polar_amp_pixels,0);
    if (this_debug>1)   SDFA_instance.showArrays(polar_green_mask_pixels0,pol_instace.getWidth(),pol_instace.getHeight(), "PGM0");
    double [] cart_green_mask_pixels0= pol_instace.polarToCartesian (polar_green_mask_pixels0,size,0.0);
    if (this_debug>1)   SDFA_instance.showArrays(cart_green_mask_pixels0,size,size,   "CGM0");


   
 //   double sigma=pol_instace.getWidth()*0.5;
/*
    double [] polar_amp_windowed=polar_amp_pixels.clone();
    oneDPolarFFTWindow(polar_amp_windowed, // half circle  polar array
                   pol_instace.getWidth(), // width (number of radial pixels in the polarPixels array)
                                    sigma); // multiply each radial ray by the gaussian window (if>0.0)
  
    double [] polar_amp_spectrum=oneDPolarFFTamp(polar_amp_windowed, // half circle  polar array
                                          pol_instace.getWidth()); // width (number of radial pixels in the polarPixels array)
    int polar_amp_spectrum_width=polar_amp_spectrum.length/pol_instace.getHeight();
*/
//int hwndr=5; // sliding average filterwindow half-width
//int hwnda=5; // sliding average filterwindow half-width
//double threshold=0.025;
//int ignoreNear=6;
//int ignoreFar=50;
//int noSkipFar=35; //.no minimums should be skipped closer than that (>2*ignoreNear)
//double tolerance=1.5; // period should not pol_instace.getHeight()deviate from the average more than that


//double fracPer=0.25;   // patch +/- this fraction of the periodic
//double patchCoeff=0.75;  // patch this part of the full 

//double bonus=0.02;
double abonus=bonus/pol_instace.getWidth();

    if (this_debug>1) {
       SDFA_instance.showArrays(polar_amp_pixels,  pol_instace.getWidth(),pol_instace.getHeight(),     "polar-amp"); // only top half+1 will be used
 //      SDFA_instance.showArrays(polar_amp_windowed,pol_instace.getWidth(),pol_instace.getHeight(),     "amp-wind"); // only top half+1 will be used
 //      SDFA_instance.showArrays(polar_amp_spectrum,polar_amp_spectrum_width,pol_instace.getHeight(),   "SPECTRUM");
/*
      double [] amp_radialHighpass=polar_amp_pixels.clone();
      radialHighpassSliding(amp_radialHighpass,
                        pol_instace.getWidth(),
                                         hwndr); // half averaging window width
//       SDFA_instance.showArrays(amp_radialHighpass,pol_instace.getWidth(),pol_instace.getHeight(),   "R-HIPASS"+hwndr);

      double [] amp_angularLowpass=amp_radialHighpass.clone();
      angularLowpassSliding(amp_angularLowpass,
                        pol_instace.getWidth(),
                                         hwnda); // half averaging window width
       SDFA_instance.showArrays(amp_angularLowpass,pol_instace.getWidth(),pol_instace.getHeight(),   "A-LOPASS"+hwnda);


      double [] radial_minmax=amp_angularLowpass.clone();
      polarDiscriminateMinMax(radial_minmax,
                     pol_instace.getWidth(),
                                  threshold, // half averaging window width
                                      true);
      SDFA_instance.showArrays(radial_minmax,pol_instace.getWidth(),pol_instace.getHeight(),   "MINMAX"+threshold);
      double [][] minData=polarDetectPeriodic (radial_minmax,
                                      pol_instace.getWidth(),
                                                  ignoreNear, // ignore wrong minimums closer to the center (not more that 1 period there)
                                                   ignoreFar, // ignore wrong minimums farter from the center_for_g2
                                                   noSkipFar, //.no minimums should be skipped closer than that (>2*ignoreNear)
                                                   tolerance); // period should not deviate from the average more than that
      double [] radial_minmax_filtered=radial_minmax.clone();
      for (int i=0;i<pol_instace.getHeight();i++) for (int j=0;j<pol_instace.getWidth();j++) radial_minmax_filtered[i*pol_instace.getWidth()+j]*=(int) minData[i][0];

      SDFA_instance.showArrays(radial_minmax_filtered,pol_instace.getWidth(),pol_instace.getHeight(),   "FILTERED-MINMAX");
*/
      double [] amp_radialPatched=polar_amp_pixels.clone();
/*
      polarPatchMins(amp_radialPatched,
                pol_instace.getWidth(),
                               minData,   // pairs of number of mins/period for each angle
                               fracPer,   // patch +/- this fraction of the periodic
                            patchCoeff);  // patch this part of the full 
       SDFA_instance.showArrays(amp_radialPatched,  pol_instace.getWidth(),pol_instace.getHeight(),     "amp-patched"); // only top half+1 will be used
*/
      for (int i=0;i<pol_instace.getHeight();i++) for (int j=0;j<pol_instace.getWidth();j++) amp_radialPatched[i*pol_instace.getWidth()+j]*=1.0+abonus*j;
       SDFA_instance.showArrays(amp_radialPatched,  pol_instace.getWidth(),pol_instace.getHeight(),     "amp-mod"); // only top half+1 will be used

       double [] polar_green_mask_pixels1=pol_instace.genPolarGreenMask(amp_radialPatched,0);
       if (this_debug>1)   SDFA_instance.showArrays(polar_green_mask_pixels1,pol_instace.getWidth(),pol_instace.getHeight(), "PGM1");
       double [] cart_green_mask_pixels1= pol_instace.polarToCartesian (polar_green_mask_pixels1,size,0.0);
       if (this_debug>1)   SDFA_instance.showArrays(cart_green_mask_pixels1,size,size,   "CGM1");

       double [] polar_green_mask_pixels2=pol_instace.genPolarGreenMask(amp_radialPatched,1);
       if (this_debug>1)   SDFA_instance.showArrays(polar_green_mask_pixels2,pol_instace.getWidth(),pol_instace.getHeight(), "PGM1");
       double [] cart_green_mask_pixels2= pol_instace.polarToCartesian (polar_green_mask_pixels2,size,0.0);
       if (this_debug>1)   SDFA_instance.showArrays(cart_green_mask_pixels2,size,size,   "CGM1");



    return cart_green_mask_pixels1;


    }





    return cart_green_mask_pixels0;

  }


/*
  void oneDPolarFFTWindow(double [] polarPixels, // half circle  polar array
                                      int width, // width (number of radial pixels in the polarPixels array)
                                   double sigma){ // multiply each radial ray by the gaussian window (if>0.0)
    int height=polarPixels.length/width;
    double k;
    if (sigma<=0.0) return;
    k=0.5/sigma/sigma;
    int i,ia;
    double[] gaussian=new double [width];
    for (i=0;i< width;i++) gaussian [i]=Math.exp (-(i*i*k));
    int index=0;
    for (ia=0;ia<height;ia++) for (i=0;i<width;i++) polarPixels[index++]*=gaussian [i];
 }

  double[] oneDPolarFFTamp(double [] polarPixels, // half circle  polar array
                                     int width){ // width (number of radial pixels in the polarPixels array)
    int fftSize=1;
    while (fftSize<(2*width)) fftSize<<=1;
    int height=polarPixels.length/width;
    int i,j,ia;
    int hSize=fftSize/2;
    double [] result = new double [height*(hSize+1)];
    Complex [] crow= new Complex [fftSize];
    for (ia=0;ia<height;ia++) {
      crow[0]=    new Complex(polarPixels[ia*width],0.0);
      crow[hSize]=new Complex(0.0,0.0);
      for (i=1;i<hSize;i++) {
        if (i>=width) {
          crow[i]=    new Complex(0.0,0.0);
        } else {
          crow[i]=    new Complex(polarPixels[ia*width+i],0.0);
        }
        crow[fftSize-i]=crow[i];
      }
      crow=fft(crow);
      for (i=0;i<hSize;i++) {
        result[ia*(hSize+1)+i]=crow[i].abs();
      }
    }
    return result;
 }
*/
 void radialLowpassSliding(double []polarArray,
                                      int width,
                                       int hwnd // half averaging window width
                                          ) {
   int height=polarArray.length/width;
   double []line=new double [width];
   int j,j0,j1,ia;
   int base;
   double k=1/(hwnd*2+1);
   for (ia=0;ia<height;ia++) {
      base=ia*width;
      line[0]=k*polarArray[base];
      for (j=1;j<=hwnd;j++) line[0]+=2*k*polarArray[base+j];
      for (j=1;j<width;j++) {
        j0=j-hwnd-1;
        j1=j+hwnd;
        if (j0<0) j0=-j0;
        if (j1>=width) j1=width-1;
        line[j]=line[j-1]+k*(polarArray[base+j1]-polarArray[base+j0]);
      }
      for (j=0;j<width;j++) polarArray[base++]=line[j];
   }
 }

 void radialHighpassSliding(double []polarArray,
                                      int width,
                                       int hwnd // half averaging window width
                                          ) {
   int height=polarArray.length/width;
   double []line=new double [width];
   int j,j0,j1,ia;
   int base;
   double k=1.0/(hwnd*2+1);
   for (ia=0;ia<height;ia++) {
      base=ia*width;
      line[0]=k*polarArray[base];
      for (j=1;j<=hwnd;j++) line[0]+=2*k*polarArray[base+j];
//      System.out.println ( "radialHighpassSliding() ia="+ia+" line[0]="+line[0]+" width="+width+" hwnd="+hwnd);
      for (j=1;j<width;j++) {
        j0=j-hwnd-1;
        j1=j+hwnd;
        if (j0<0) j0=-j0;
        if (j1>=width) j1=width-1;
        line[j]=line[j-1]+k*(polarArray[base+j1]-polarArray[base+j0]);
      }
      for (j=0;j<width;j++) polarArray[base++]-=line[j];
   }
 }


/* assuming angular-periodic */
 void angularLowpassSliding(double []polarArray,
                                      int width,
                                       int hwnd // half averaging window width
                                          ) {
   int height=polarArray.length/width;
   double []line=new double [height];
   int i,ir;
//   int base;
   double k=1.0/(hwnd*2+1);
   for (ir=0;ir<width;ir++) {
//      base=ia*width;
      line[0]=0.0;
      for (i=-hwnd;i<=hwnd;i++) line[0]+=k*polarArray[width*((i+height)%height)+ir];
//      System.out.println ( "angularLowpassSliding() ir="+ir+" line[0]="+line[0]+" width="+width+" hwnd="+hwnd);
      for (i=1;i<height;i++) {
        line[i]=line[i-1]+k*(polarArray[width*((i+hwnd)%height)+ir]-polarArray[width*((i-hwnd+height-1)%height)+ir]);
      }
      for (i=0;i<height;i++) polarArray[width*i+ir]=line[i];
   }
 }

 void polarDiscriminateMinMax(double []polarArray,
                                        int width,
                                 double threshold, // half averaging window width
                             boolean startWithMax
                                          ) {
   int height=polarArray.length/width;
   int []line=new int [width];
   int ia,i;
   double min,max;
   int base;
   boolean lookingForMax=startWithMax; //looking for maximum
   int lastMin,lastMax;
   double d;
   for (ia=0;ia<height;ia++) {
     base=width*ia;
     lastMin=0;
     lastMax=0;
     max=polarArray[base+lastMax];
     min=max;
     for (i=0;i<width;i++) line[i]=0;

     for (i=1;i<width;i++) {
       d=polarArray[base+i];
       if (lookingForMax) {
          if (d>max) {
            max=d;
            lastMax=i;
          } else if (d<(max-threshold)) {
            line[lastMax]=1;
            min=d;
            lastMin=i;
            lookingForMax=false;
          }
       } else { // looking for min
          if (d<min) {
            min=d;
            lastMin=i;
          } else if (d>(min+threshold)) {
            line[lastMin]=-1;
            max=d;
            lastMax=i;
            lookingForMax=true;
          }
       }
     }
     if (lookingForMax) line[lastMax]=1;
     else               line[lastMin]=-1;
     for (i=1;i<width;i++) {
       polarArray[base+i]=line[i]; 
     }
   }
 }


 void polarPatchMins(double []polarArray,
                               int width,
                        double [][] mins,   // pairs of number of mins/period for each angle
                          double fracPer,   // patch +/- this fraction of the periodic
                       double patchCoeff){  // patch this part of the full 
   int height=polarArray.length/width;
   int ia,n,i,i0,i1;
   int base;
   double d0,k;
   for (ia=0;ia<height;ia++) {
     for (n=1;n<=mins[ia][0];n++) {
       i0=(int) Math.round (mins[ia][1]*(n-fracPer));
       i1=(int) Math.round (mins[ia][1]*(n+fracPer));
       base=ia*width;
       d0=polarArray[base+i0];
       k=(polarArray[base+i1]-d0)/(i1-i0);
       for (i=i0+1;i<i1;i++) {
         polarArray[base+i]=polarArray[base+i]*(1.0-patchCoeff)+patchCoeff*(d0+k*(i-i0));
       }
     }
   }
 }


 // pairs - number of qualified minimums, period (number==0 - no minimums)
 double [][] polarDetectPeriodic (double []minMaxArray,
                                        int width,
                                   int ignoreNear, // ignore wrong minimums closer to the center (not more that 1 period there)
                                    int ignoreFar, // ignore wrong minimums farter from the center_for_g2
                                    int noSkipFar, //.no minimums should be skipped closer than that (>2*ignoreNear)
                                 double tolerance // period should not deviate from the average more than that
                        ) {
   int height=minMaxArray.length/width;
   int []allMins= new int [width];
   double [][]result= new double [height][2];
   int ia,i;
   int base;
   int numInThis;
   int ilow,ihigh;
   double period;
   int nCycles;
   double err,d;
   for (ia=0;ia<height;ia++) {
     result[ia][0]=0.0;
     result[ia][1]=0.0;
     base=width*ia;
     numInThis=0;
     for (i=0;i<width;i++) if (minMaxArray[base++]<0) allMins[numInThis++]=i;
     ilow=0;
     while ((ilow<numInThis) && (allMins[ilow]<ignoreNear)) ilow++;
     ihigh=numInThis-1;
     while ((ihigh>=0) && (allMins[ihigh]>ignoreFar)) ihigh--;
     if ((ilow<numInThis) && (ihigh>=0)) { // there are sum minimums
       if (ihigh==ilow) {
         period=allMins[ihigh];
         nCycles=1;
       } else {
         nCycles=(int)Math.round(1.0*(ihigh-ilow)*allMins[ihigh]/(allMins[ihigh]-allMins[ilow]));
         period=1.0*allMins[ihigh]/nCycles;
         err=0;
         for (i=ilow;i<=ihigh;i++) {
           d=period*(nCycles+i-ihigh)-allMins[i];
           err+=d*d;
//      System.out.println ( "polarDetectPeriodic() ia"+ia+" i="+i+" d="+d);

         }
         err=Math.sqrt(err/(ihigh-ilow+1));
         if (err>tolerance) continue; 
       }
/* extend limits */
       while ((ihigh< (numInThis-1)) && (Math.abs(allMins[ihigh+1]- period*nCycles)<1.5*tolerance) ) {
         ihigh++;
         nCycles++;
       }
/* See if required far minimums are present */
       if ((period* (nCycles+1)+tolerance)< noSkipFar) continue; // failed to meet no-skip up to noSkipFar
/* near minimums only extra allower, but not skipping the needed, and no extras after first needed period */
       while ((ilow>0) && ((nCycles-ihigh+ilow)>1) && (Math.abs(allMins[ilow]- period*(nCycles-ihigh+ilow))<1.5*tolerance)) ilow--;
       if (((nCycles-ihigh+ilow)>1) || (Math.abs(allMins[ilow]- period)>1.5*tolerance)) continue;
/* Verify each actual period */
       result[ia][0]=nCycles;
       for (i=ilow+1;i<ihigh;i++) if (Math.abs(allMins[i]- allMins[i-1]-period)>1.5*tolerance) {
         result[ia][0]=0.0;
         break;
       }
       if (result[ia][0]==0.0) continue;
       result[ia][1]=1.0*allMins[ihigh]/nCycles;
//      System.out.println ( "polarDetectPeriodic():"+ia+":"+result[ia][0]+":"+result[ia][1]);
     }
   }   
   return result;
 }



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


  public void  processColorsWeights(ImageStack stack,
                                        double scale,     // initila maximal pixel value (16))
                                  double balance_red,     // 1.0; // manual color balance, gain 1.0 matches 0.0.255.0 range of the unput Bayer data
                                 double balance_blue,    // 1.0;
                                         double gain,      // 1.0;
                               double weight_scale_r,  // 1.0; // additional correction for the weights of colors for different sub-pixels in a Bayer cell
                               double weight_scale_b,  // 1.0; // WEIGHT_SCALE_G=1.0
                                        double sigma,     // 2.0; // Gaussian sigma to low-pass color components when calculating "smooth" color
                                        double gamma,     // 0.53;
                                       double minLin,
                                           double Kr,        // 0.299;
                                           double Kb,        // 0.114;
                                   double saturation_red,
                                   double saturation_blue // 2.0
                           ) {
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
      double gain_red= balance_red* gain/scale;
      double gain_blue=balance_blue*gain/scale;
      double gain_green=gain/scale;
      double gamma_a=Math.pow(minLin,gamma)*(1.0-gamma);
      gamma_a=gamma_a/(1.0-gamma_a);
      double gamma_linK=(1.0+gamma_a)*gamma*Math.pow(minLin,gamma)/minLin;



      for (i=0;i<length;i++) {
        fpixels_r[i]=(float) linGamma(gamma, gamma_a, gamma_linK, minLin, fpixels_r[i]*gain_red);
        fpixels_g[i]=(float) linGamma(gamma, gamma_a, gamma_linK, minLin, fpixels_g[i]*gain_green);
        fpixels_b[i]=(float) linGamma(gamma, gamma_a, gamma_linK, minLin, fpixels_b[i]*gain_blue);
      }
/* Convert to YPbPr */
      double Y,Pb,Pr;
      double Kg=1.0-Kr-Kb;
      double Sb=0.5/(1.0-Kb)*saturation_blue;
      double Sr=0.5/(1.0-Kr)*saturation_red;
      double Yr,Yg,Yb,Wr,Wg,Wb,S;
/* coefficients to find Y from Pb, Pr and a color (R,G or B)
 Yr = R- Pr*KPrR
 Yb = B- Pb*KPbB
 Yg = G+ Pr*KPrG  + Pb*KPbG
 */
      double KPrR= -(2.0*(1-Kr))/saturation_red;
      double KPbB= -(2.0*(1-Kb))/saturation_blue;
      double KPrG=  2.0*Kr*(1-Kr)/Kg/saturation_red;
      double KPbG=  2.0*Kb*(1-Kb)/Kg/saturation_blue;
      if (DEBUG_LEVEL>1) {
          System.out.println ( " processColorsWeights() gain_red="+gain_red+" gain_green="+gain_green+" gain_blue="+gain_blue);
          System.out.println ( " processColorsWeights() gamma="+gamma+      " minLin="+minLin+" gamma_a="+gamma_a+" gamma_linK="+gamma_linK);
          System.out.println ( " processColorsWeights() Kr="+Kr+" Kg="+Kg+" Kb="+Kb+" Sr="+Sr+" Sb="+Sb);
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
        Y=Kr*fpixels_r[i]+Kg*fpixels_g[i]+Kb*fpixels_b[i];
        fpixels_pb[i] = (float) (Sb*(fpixels_b[i]-Y));
        fpixels_pr[i] = (float) (Sr*(fpixels_r[i]-Y));
        fpixels_y0[i]=(float) Y;
      }
/* Low-pass filter Pb and Pr */
      FloatProcessor fp_pb = new FloatProcessor(width, height, fpixels_pb, null);
      FloatProcessor fp_pr = new FloatProcessor(width, height, fpixels_pr, null);
      GaussianBlur gb = new GaussianBlur();
      gb.blurFloat(fp_pb, sigma, sigma, 0.01);
      gb.blurFloat(fp_pr, sigma, sigma, 0.01);
      stack.addSlice("Pr", fpixels_pr);
      stack.addSlice("Pb", fpixels_pb);
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
          S=1.0/(Wr*(weight_scale_r-1.0)+Wb*(weight_scale_b-1.0)+1.0);
          Wr*=S*weight_scale_r;
          Wb*=S*weight_scale_b;
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
      stack.addSlice("Y",  fpixels_y);
      stack.addSlice("Y0", fpixels_y0); // not filtered by low-pass, preliminary (for comaprison only)

      if (DEBUG_LEVEL>2) {
        stack.addSlice("Yr",fpixels_yR);
        stack.addSlice("Yg",fpixels_yG);
        stack.addSlice("Yb",fpixels_yB);
      }

  }
// fpixels_r[i]=(float) linGamma(gamma, gamma_a, gamma_linK, minLin, fpixels_r[i]*gain_red);
  public double linGamma(double gamma, double a, double k, double x0, double x) {
    if (x<0) return 0.0;

    if (x<=x0) return k*x;
    return (1.0+a)*Math.pow(x,gamma)-a;
//  return x;
  }


/* Combine 2 stacks and a mask */
  public ImageStack combineStacksWithMask (ImageStack stack_bg,
                                           ImageStack stack_fg, 
                                                 float [] mask ) {
    
    ImageStack stack=new ImageStack(stack_bg.getWidth(),stack_bg.getHeight());
    int slice,i;
    float [] fpixels;
    float [] fpixels_bg;
    float [] fpixels_fg;
    for (slice=1; slice <=stack_bg.getSize(); slice++) {
       fpixels_bg= (float[])stack_bg.getPixels(slice);
       fpixels_fg= (float[])stack_fg.getPixels(slice);
       fpixels=new float [fpixels_bg.length];
       for (i=0;i<fpixels_bg.length;i++) fpixels[i]= mask[i]*fpixels_fg[i]+(1.0f-mask[i])*fpixels_bg[i];
       stack.addSlice(stack_fg.getSliceLabel(slice), fpixels);
    }
    return stack;
  }


  public ImageStack combineRGBCorrectBayer ( double [] red_pixels,
                                           double [] green_pixels,
                                            double [] blue_pixels,
                                                         int width,   // image width
                                                      int rb_shift){  // half Bayer period (GR/BG)
    if (red_pixels==null) {
       System.out.println("combineRGBCorrectBayer(): red_pixels==null");
       return null;
    }
    if (green_pixels==null) {
       System.out.println("combineRGBCorrectBayer(): green_pixels==null");
       return null;
    }
    if (blue_pixels==null) {
       System.out.println("combineRGBCorrectBayer(): blue_pixels==null");
       return null;
    }
    if ((red_pixels.length!=green_pixels.length) || (blue_pixels.length!=green_pixels.length)) {
       System.out.println("Different lengths: red_pixels.length="+red_pixels.length+
                                        " green_pixels.length="+green_pixels.length+
                                          " blue_pixels.length="+blue_pixels.length);
       return null;
    }
    int height=green_pixels.length/width;
    ImageStack stack=new ImageStack(width,height);
    float [] fpixels= new float [green_pixels.length];
    int iy,ix,index;
    for (index=0;index<green_pixels.length;index++) {
      ix=index%width;
      if (ix<rb_shift) fpixels[index]=0.0f;
      else fpixels[index]=(float) red_pixels[index-rb_shift];
    }
    stack.addSlice("red", fpixels);
    fpixels= new float [green_pixels.length];
    for (index=0;index<green_pixels.length;index++) {
      fpixels[index]=(float) green_pixels[index];
    }
    stack.addSlice("green", fpixels);
    fpixels= new float [green_pixels.length];
    int shift=rb_shift*width;
    for (index=0;index<green_pixels.length;index++) {
      iy=index/width;
      if (iy<rb_shift) fpixels[index]=0.0f;
      else fpixels[index]=(float) blue_pixels[index-shift];
    }
    stack.addSlice("blue", fpixels);
    return stack;
  }




  public ImageStack combineRGBCorrectBayer ( float [] red_pixels,
                                           float [] green_pixels,
                                            float [] blue_pixels,
                                                        int width,   // image width
                                                     int rb_shift){  // half Bayer period (GR/BG)
    if (red_pixels==null) {
       System.out.println("combineRGBCorrectBayer(): red_pixels==null");
       return null;
    }
    if (green_pixels==null) {
       System.out.println("combineRGBCorrectBayer(): green_pixels==null");
       return null;
    }
    if (blue_pixels==null) {
       System.out.println("combineRGBCorrectBayer(): blue_pixels==null");
       return null;
    }
    if ((red_pixels.length!=green_pixels.length) || (blue_pixels.length!=green_pixels.length)) {
       System.out.println("Different lengths: red_pixels.length="+red_pixels.length+
                                        " green_pixels.length="+green_pixels.length+
                                          " blue_pixels.length="+blue_pixels.length);
       return null;
    }
    int height=green_pixels.length/width;
    ImageStack stack=new ImageStack(width,height);
    float []fpixels= new float [green_pixels.length];
    int iy,ix,index;
    for (index=0;index<green_pixels.length;index++) {
      ix=index%width;
      if (ix<rb_shift) fpixels[index]=0.0f;
      else fpixels[index]=red_pixels[index-rb_shift];
    }
    stack.addSlice("red", fpixels);
    stack.addSlice("green", green_pixels.clone());
    fpixels= new float [green_pixels.length];
    int shift=rb_shift*width;
    for (index=0;index<green_pixels.length;index++) {
      iy=index/width;
      if (iy<rb_shift) fpixels[index]=0.0f;
      else fpixels[index]=blue_pixels[index-shift];
    }
    stack.addSlice("blue", fpixels);
    return stack;
  }

/* returns 2 masks (0:0 in the top left corner, match fht) [0] - for greens, [1] - for red/blue */
/* Possible improvements: - 1 make the initial green mask (or actually "fan"-like image) to have sharper ends.
                             2. detect periodic (line of spots) on the spectrum aplitudes (corresponds to thin lines) and use this
                                info to confirm this area to belong to the main spectrum */

  public double [][] aliasScissors(double [] green_fht, // fht array for green, will be masked in-place
                                           int subpixel, // measured kernels are sampled at subpixel higher than Bayer (subpixel/2 - than sensor) //4
                               double debayer_threshold, // no high frequencies - use default uniform filter
//                          double debayer_relative_width, // Debayer lopass filter width (relative to distance to the nearest alias)
                                   double debayer_gamma, // power function applied to the amplitudes before generating spectral masks
                                      double debayer_rz, // for green mask - rays start radius from center, relative to distance to the first alias
                                      double debayer_ra, // for green mask - rays radius around aliases, relative to distance to the first alias
                                   double debayer_sigma, // for green mask - reduce value of the far amplitudes, relative to distance to the first alias
                                   double debayer_decay, // for green mask - exponential decay of the ray value, relative to distance to the first alias
                            double debayer_farthest_max, // for green mask - rays will be extyended from maf back to the center, max up to this distance (relative)
                            double debayer_radius_power, // for green mask -  divide ray values by the radius to this power
                                     double mainToAlias,// relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
                               double debayer_mask_blur, // for both masks  sigma for gaussian blur of the binary masks (<0 -do not use "scissors")
                               boolean debayer_lo_green, //  combine alias-reject "scissors" with lopass filter for greens
                           boolean debayer_lo_postgreen, // combine alias-reject "scissors" with lopass filter for greens applied to red/blue
                             boolean debayer_lo_redblue, // combine alias-reject "scissors" with lopass filter for greens applied to red/blue
                             PolarSpectrums pol_instace, // initialized instance or null
                                   double [][][] lopass, // [1.0,scaled][green,redBlue][size*size] - externally prepared arrays, centered in the center
                                         int this_debug){ // internal debug level
    int length=green_fht.length;
    int size=(int) Math.sqrt(length);
    double [] green_mask=null;
    double [] red_blue_mask=null;
    double [] green_amp=fht_instance.calculateAmplitude(green_fht);
    int i,j;
/**normailze amplitudes, apply gamma */
    double dmax=0.0;
    for (i=0;i<green_amp.length;i++) if (green_amp[i]>dmax) dmax=green_amp[i];
    dmax=1.0/dmax;
    for (i=0;i<green_amp.length;i++) green_amp[i]= Math.pow(green_amp[i]*dmax,debayer_gamma);
    if (this_debug>2)   SDFA_instance.showArrays(green_amp,  "DT-gam"); // only top half+1 will be used
    double midRangeSpectral=pol_instace.maxAmpInRing (green_amp);
    boolean useFancyDebayer=(midRangeSpectral>=debayer_threshold);

    DOUBLE_DEBUG_RESULT= midRangeSpectral; 

    if (useFancyDebayer && (debayer_mask_blur>=0)) { /* calculate and apply "scissors" masks */
      if (pol_instace!=null) {
        green_mask= calcGreensAliasMaskRays (green_amp, // normalized amplitude spectrum, (0,0) in the center
                                           pol_instace, // initialized instance
                                        -debayer_decay, // hack - here it is "bonus"
                                            this_debug);//
      }
      /*else {
        green_mask= calcGreensAliasMask (green_amp,
                                          subpixel,
                                        debayer_rz,
                                        debayer_ra,
                                     debayer_sigma,
                                     debayer_decay,
                              debayer_farthest_max, // fartherst absolute maximum on a ray to count
                              debayer_radius_power,// divide ray values by the radius to this power
//                                     mainToAlias);// relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
                                              1.0);// relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
      }
*/
      if (debayer_mask_blur>0) {
        if (this_debug>3) SDFA_instance.showArrays(green_mask,  "G-raw");
        blurDouble(green_mask,   size, debayer_mask_blur, debayer_mask_blur, 0.01);
        if (this_debug>3) SDFA_instance.showArrays(green_mask,  "G-blurred");
      }
      double [] green_mask_post;
      if (debayer_lo_postgreen && !debayer_lo_green) {
        green_mask_post=green_mask.clone();
        for (i=0;i<green_mask_post.length;i++) green_mask_post[i]*=lopass[1][0][i]; //  scaled, green - was green_lopass[i];
      } else if (!debayer_lo_postgreen && debayer_lo_green) {
        green_mask_post=green_mask.clone();
        for (i=0;i<green_mask.length;i++) green_mask[i]*=lopass[1][0][i]; //  scaled, green - was green_lopass[i];
      } else  { // both the same
        if (debayer_lo_green) for (i=0;i<green_mask.length;i++) green_mask[i]*=lopass[1][0][i]; //  scaled, green - was green_lopass[i];
         green_mask_post=green_mask;
      }

      if (pol_instace!=null) {
/* Maybe here we need to unmasked (wide bandwidth) green_amp? */
        red_blue_mask= calcRedBlueAliasMaskRays (green_amp, // both halves are needed ??
                                           green_mask_post, // may be null if amp_pixels is already masked
                                               pol_instace, // initialized instance (if null - skip rays processing)
                                               mainToAlias,// relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
                                            -debayer_decay, // scale far pixels as (1.0+bonus*r/rmax)
                                               this_debug);// relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
      }/* else {
        red_blue_mask= calcRedBlueAliasMask (green_amp,
                                       green_mask_post,
                                              subpixel,
                                           mainToAlias);// relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
      }
*/



/* add    double mainToAlias){// relative main/alias amplitudes to enable pixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out) */

      if (debayer_mask_blur>0) {
        if (this_debug>3) SDFA_instance.showArrays(red_blue_mask,  "RB-raw");
        blurDouble(red_blue_mask, size,debayer_mask_blur, debayer_mask_blur, 0.01);
        if (this_debug>3) SDFA_instance.showArrays(red_blue_mask,  "RB-blurred");
      }
      if (debayer_lo_redblue) for (i=0;i<red_blue_mask.length;i++) red_blue_mask[i]*=lopass[1][1][i]; //  scaled, red-blue - was red_blue_lopass[i];
    } else { // debayer_mask_blur<0 : use default masks
       green_mask=lopass[1][0].clone(); //green_lopass.clone(); variable (wide) filter here)
       red_blue_mask=lopass[1][1].clone(); //red_blue_lopass.clone();
       if (!useFancyDebayer) for (i=0;i<green_mask.length;i++) { // no high-frequency componnets detected - reduce noise by extra (narrow) filtering
         green_mask[i]*=   lopass[0][0][i]; // *=   green_lopass[i];
         red_blue_mask[i]*=lopass[0][1][i]; // *=red_blue_lopass[i];
       }
    }
/* Swap quadrants in the masks to match FHT arrays (0:0 in the top left corner) */
    fht_instance.swapQuadrants(green_mask);
    fht_instance.swapQuadrants(red_blue_mask);
/* return both masks */
    double [][] result =new double [2][];
    result[0]= green_mask;
    result[1]= red_blue_mask;
//    if (this_debug>3) SDFA_instance.showArrays(result,  "before_norm_masks");


/* normalize masks to have exactly 1.0 at 0:0 - it can be reduced by blurring */
    for (i=0;i<result.length;i++) {
      dmax=1.0/result[i][0];
      for (j=0;j<result[i].length;j++) result[i][j]*=dmax;
    }
//    if (this_debug>3) SDFA_instance.showArrays(result,  "masks");
    return result;
  }

  public float [][] convolveImageWithKernels(ImagePlus imp, // 32-bit (float) linearized image with Bayer mosaic
                                               Rectangle r, // area to process (if not all the image) in pixels
                                   boolean equalize_greens, // equalize two green gains (currently - independently for each tile)
                                   double [][][][] kernels, // 2d array of convolution kernels 
//                            double [][][][][][] amplPhases, // 2d array of per component ampliteudes/full phases (for interpolation)
                                              int subpixel, // measured kernels are sampled at subpixel higher than Bayer (subpixel/2 - than sensor) //4
                                            int kernelStep, // distance between kernel samples, in sensor pixels // 64 // 32
                                         int kernelMargins, // number of kernel steps from the image 0,0 points to the center of the (0,0) kernel // 4
                                               int fftSize, // size of sliding FFT (twice the kernel size) here:128
                                              int compMask, // bitmask of color components to process
                                       boolean filterBayer, // compensate for Bayer aliases
/* Parameters related to reducing mosaic alias artifacts, can be removed if that filtering will be preformed before aberration correction.
     That just requires that lateral chtomatic aberration is small relative to FFT size, so spectrums amplitudes are similar */
                                      double debayer_gamma, // power function applied to the amplitudes before generating spectral masks
                                         double debayer_rz, // for green mask - rays start radius from center, relative to distance to the first alias
                                         double debayer_ra, // for green mask - rays radius around aliases, relative to distance to the first alias
                                      double debayer_sigma, // for green mask - reduce value of the far amplitudes, relative to distance to the first alias
                                      double debayer_decay, // for green mask - exponential decay of the ray value, relative to distance to the first alias
                               double debayer_farthest_max, // for green mask - rays will be extyended from maf back to the center, max up to this distance (relative)
                               double debayer_radius_power, // for green mask -  divide ray values by the radius to this power
                                        double mainToAlias, // relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
                                  double debayer_mask_blur, // for both masks  sigma for gaussian blur of the binary masks
                                  boolean debayer_lo_green, //  combine alias-reject "scissors" with lopass filter for greens
                              boolean debayer_lo_postgreen, // combine alias-reject "scissors" with lopass filter for greens applied to red/blue
                                boolean debayer_lo_redblue, // combine alias-reject "scissors" with lopass filter for greens applied to red/blue
/* end of demosaic parameters */
                                 float [][] bayerWeightsRB, // if non-null, should be a three element array, each will be assigned weight-red and weight_blue
                                      boolean updateStatus){  // update status info

    ImageProcessor ip=imp.getProcessor();
    Rectangle rroi=new Rectangle (r);
    if (rroi.x<0) rroi.x=0;
    if (rroi.y<0) rroi.y=0;
    rroi.x&=~1;
    rroi.x&=~1;
    if ((rroi.x+rroi.width )>ip.getWidth ()) rroi.width= ip.getWidth() -rroi.x;
    if ((rroi.y+rroi.height)>ip.getHeight()) rroi.height=ip.getHeight()-rroi.y;
    rroi.width &= ~1;
    rroi.height &= ~1;
    int i,j,k,l,chn;
    int nChn=0;
    for (i=compMask; i!=0;i>>=1) nChn++;
    float [][] result = new float [nChn][];
    int outWidth= ip.getWidth()* subpixel/2;
    int outHeight=ip.getHeight()*subpixel/2;

    for (i=0;i<nChn;i++) if ((compMask & (1<<i))!=0) {
      result[i]=new float [outWidth*outHeight];
      for (j=0;j<result[i].length;j++) result[i][j]=(float) 0.0;
    }

    int imgStep=fftSize/4/(subpixel/2); // Tile overlap in sensor pixels (16)cosMask
    int xTile0=rroi.x/kernelStep;
    int yTile0=rroi.y/kernelStep;
    int xTile1=(rroi.x+rroi.width-1)/kernelStep+1;
    int yTile1=(rroi.y+rroi.height-1)/kernelStep+1;
    int xTile,yTile,xSubTile,ySubTile,iy,ix;

//    double [] slidingMask=getSlidingMask(kernelStep*subpixel); // 32*4=128
    int slidingMaskSize=imgStep*subpixel; // 16*4=64
    double [] slidingMask=getSlidingMask(slidingMaskSize); //(imgStep*subpixel); // 16*4=64

    double [][] bayer_patterns=null;
/* Later computetion could be greately decresed by at least using half-size fft, or just convolving directly (only 4x4=16 points are needed) */
    boolean hasRGB=((compMask & (1<< CHANNEL_RED))!=0) && ((compMask & (1<< CHANNEL_BLUE))!=0) && ((compMask & (1<< CHANNEL_GREEN))!=0); // has channels 1,2, and 5
    boolean createBayerWeights=(bayerWeightsRB!=null) && (bayerWeightsRB.length>=2) && hasRGB;
    float [][] bayerWeights=null; // same first dimension as result[][], in the end will be re-calculated to a two-=element bayerWeightsRB[2][]




    if (createBayerWeights) {
      bayerWeights = new float [nChn][];
      for (i=0;i<nChn;i++) if ((compMask & (1<<i))!=0) {
        bayerWeights[i]=new float [outWidth*outHeight];
        for (j=0;j<bayerWeights[i].length;j++) bayerWeights[i][j]=(float) 0.0;
      }
      for (i=0;i<2;i++); for (j=0;j<bayerWeights[i].length;j++) bayerWeights[i][j]=(float) 0.0;
      bayer_patterns=new double [2][slidingMaskSize*slidingMaskSize];
      for (iy=0;iy<slidingMaskSize;iy++) for (ix=0;ix<slidingMaskSize;ix++) {
        bayer_patterns[0][iy*slidingMaskSize+ix]=((( iy % subpixel)==0) && ((ix % subpixel)==0 ))?1.0:0.0; // all colors but checker
        bayer_patterns[1][iy*slidingMaskSize+ix]=((((iy % subpixel)==0) && ((ix % subpixel)==0 )) || (((iy % subpixel)==(subpixel/2)) && ((ix % subpixel)==(subpixel/2) )))?1.0:0.0; // checker greens
      }
      bayer_patterns=extendFFTInputTo (bayer_patterns, fftSize);
/* prepare bayer paterns for convolution - perform FHT transform */
      for (i=0;i<bayer_patterns.length;i++) {
                fht_instance.swapQuadrants(bayer_patterns[i]);
                fht_instance.transform(bayer_patterns[i]);
      }
    }
    int rb_shift=subpixel/2; //half Bayer period
    double [][] input_bayer;
    Rectangle rTile=new Rectangle (2*imgStep,2*imgStep); // 64x64
    int subTiles=kernelStep/imgStep; // 4 
//    double [][][][][][] amplPhasesCache=new double[subTiles+1][subTiles][][][][]; // first dimension HORIZONTAL here
    double [][][][]     fhtCache=       new double[subTiles+1][subTiles][nChn][];      // first dimension HORIZONTAL here [X][Y]
    int [][][] cachedTile=  new int [2][2][2]; // [x][y][x,y]
    int [][][] tileCorners= new int [2][2][2]; // [x][y][x,y] which kernels are needed for the 4 corners of the current tile
    double [][][][] fht_corners=new double [2][2][nChn][]; // kernels, extended and converted to amplitudes/full phases fro the 4 corners of the tile
                                                           // first index horizontal
    double [] points = new double [subTiles+1];       // interpolation coefficients, i.e. {0.0, 0.25, 0.5, 0.75, 1.0}
    double [][] fht_line = new double [subTiles+1][]; // result of interpolation of FHT arrays
    for (i=0;i<=subTiles;i++) points[i]=(1.0*i)/subTiles;
    double [] kernel;
    double [] product;
  //  double [] test;
  //  double [][] test_arr;
    if (DEBUG_LEVEL>1) System.out.println("convolveImageWithKernels(): yTile0="+yTile0+" xTile0="+xTile0+" yTile1="+yTile1+" xTile1="+xTile1+" compMask="+compMask);
if (DEBUG_LEVEL>2){
  SDFA_instance.showArrays(slidingMask,  "slidingMask");
}
    double debayer_relative_width_green=1.0; // realtive width of the low-pass filter for greens (1.0 goes to zero at exactrly the closest alias )
    double debayer_relative_width_redblue=1.0;// same for individual (re and blue) color components
    double [][][] lopass= createAliasFilters (debayer_relative_width_green, // realtive width of the low-pass filter for greens (1.0 goes to zero at exactrly the closest alias )
                                            debayer_relative_width_redblue, // same for individual (re and blue) color components
                                                                   fftSize, // side of the square
                                                                  subpixel); // should be 4 now
    double [][] both_masks=null; /* green and r/b "scissors" masks to reduce sampling aliases */

    PolarSpectrums pol_instace=new PolarSpectrums(
                fftSize, // size of the square array, centar is at size/2, size/2, only top half+line will be used
                 Math.PI, //2*Math.PI, // i.e. Math.PI, 2*Math.PI
                fftSize/2-2, // width of the polar array - should be <= size/2-2
                 0.5, //0.75, //2.0, //0.5, // maximal step in pixels on the maxRadius for 1 angular step (i.e. 0.5)
                 4 );// angular symmetry - 0- none,1 - pi corresponds to integer, 2 - pi/2 corresponds to integer, n - pi/n corresponds to integer angular step 






    for (i=0;i<2;i++) for (j=0;j<2;j++)  for (k=0;k<2;k++) cachedTile[i][j][k]=-1;
    for (yTile=yTile0;yTile<=yTile1;yTile++) {
      for (xTile=xTile0;xTile<=xTile1;xTile++) {
        if (updateStatus) IJ.showStatus("Convolving image, tile "+((yTile-yTile0)*(xTile1-xTile0+1)+(xTile-xTile0)+1)+" of "+((yTile1-yTile0+1)*(xTile1-xTile0+1)));
    if (DEBUG_LEVEL>1) System.out.println("convolveImageWithKernels(): yTile="+yTile+" xTile="+xTile);
        for (i=0;i<2;i++) for (j=0;j<2;j++) {
          tileCorners[j][i][0]=xTile-kernelMargins+j;
          tileCorners[j][i][1]=yTile-kernelMargins+i;
/* fix margins */
          if (tileCorners[j][i][0]<0) tileCorners[j][i][0]=0;
          if (tileCorners[j][i][1]<0) tileCorners[j][i][1]=0;
          if (tileCorners[j][i][0]>=kernels[0].length) tileCorners[j][i][0]=kernels[0].length-1;
          if (tileCorners[j][i][1]>=kernels.length)    tileCorners[j][i][1]=kernels.length-1;
        }
if (DEBUG_LEVEL>1) System.out.println(" tileCorners[0][0]="+tileCorners[0][0][0]+"/"+tileCorners[0][0][1]+
                                      " tileCorners[0][1]="+tileCorners[0][1][0]+"/"+tileCorners[0][1][1]+
                                      " tileCorners[1][0]="+tileCorners[1][0][0]+"/"+tileCorners[1][0][1]+
                                      " tileCorners[1][1]="+tileCorners[1][1][0]+"/"+tileCorners[1][1][1]);
if (DEBUG_LEVEL>1) System.out.println(" cachedTile[0][0]=" +cachedTile[0][0][0]+"/"+ cachedTile[0][0][1]+
                                      " cachedTile[0][1]=" +cachedTile[0][1][0]+"/"+ cachedTile[0][1][1]+
                                      " cachedTile[1][0]=" +cachedTile[1][0][0]+"/"+ cachedTile[1][0][1]+
                                      " cachedTile[1][1]=" +cachedTile[1][1][0]+"/"+ cachedTile[1][1][1]);

/* not dealing with missing kernels here - that can be done before */
/* compare new corners with the cache */
        if ((tileCorners[0][0][0]!=cachedTile[0][0][0]) ||
            (tileCorners[0][0][1]!=cachedTile[0][0][1]) ||
            (tileCorners[0][1][0]!=cachedTile[0][1][0]) ||
            (tileCorners[0][1][1]!=cachedTile[0][1][1]) ||
            (tileCorners[1][0][0]!=cachedTile[1][0][0]) ||
            (tileCorners[1][0][1]!=cachedTile[1][0][1]) ||
            (tileCorners[1][1][0]!=cachedTile[1][1][0]) ||
            (tileCorners[1][1][1]!=cachedTile[1][1][1])) { // at least some of the cache needs to be re-calculated
if (DEBUG_LEVEL>1) System.out.println("new tile");
/* Is it next tile in a row? */
          if ((tileCorners[0][0][0]==cachedTile[1][0][0]) &&
              (tileCorners[0][0][1]==cachedTile[1][0][1]) &&
              (tileCorners[0][1][0]==cachedTile[1][1][0]) &&
              (tileCorners[0][1][1]==cachedTile[1][1][1])) { /**yes, just copy last interpolated column to the first */
              fhtCache[0]=fhtCache[subTiles];
if (DEBUG_LEVEL>1) System.out.println("same row");

          } else { /* no, need to recalculate first column too */
if (DEBUG_LEVEL>1) System.out.println("new row");
//            amplPhasesCache[0]=new double[subTiles][nChn][][][];
            fhtCache[0]=new double[subTiles][nChn][];
            for (chn=0;chn<nChn;chn++) {
              if ((compMask & (1<<chn))!=0)  {
    if (DEBUG_LEVEL>1) System.out.println("tileCorners[0][0][1]="+tileCorners[0][0][1]+" tileCorners[0][0][0]="+tileCorners[0][0][0]+ " chn="+chn);

                kernel=extendFFTInputTo (kernels[tileCorners[0][0][1]][tileCorners[0][0][0]][chn], fftSize);
                fht_instance.swapQuadrants(kernel);
                fht_instance.transform(kernel);
                fht_corners[0][0][chn]=kernel.clone();
                kernel=extendFFTInputTo (kernels[tileCorners[0][1][1]][tileCorners[0][1][0]][chn], fftSize);
                fht_instance.swapQuadrants(kernel);
                fht_instance.transform(kernel);
                fht_corners[0][1][chn]=kernel.clone();
                fht_line= interpolateFHT (fht_corners[0][0][chn],    // first FHT array
                                          fht_corners[0][1][chn],    // second FHT array
                                                      points);   // array of interpolation points - 0.0 - fht0, 1.0 - fht1
                for (k=0;k<subTiles;k++) fhtCache[0][k][chn]= fht_line[k];
              } else for (k=0;k<subTiles;k++) fhtCache[0][k][chn]= null;
            }
          }
/* Update what is cached */
          cachedTile[0][0][0]=tileCorners[0][0][0];
          cachedTile[0][0][1]=tileCorners[0][0][1];
          cachedTile[0][1][0]=tileCorners[0][1][0];
          cachedTile[0][1][1]=tileCorners[0][1][1];
/* interpolate the last column */
          fhtCache[subTiles]=new double[subTiles][nChn][];

          for (chn=0;chn<nChn;chn++) {
            if ((compMask & (1<<chn))!=0)  {
              kernel=extendFFTInputTo (kernels[tileCorners[1][0][1]][tileCorners[1][0][0]][chn], fftSize);
              fht_instance.swapQuadrants(kernel);
              fht_instance.transform(kernel);
              fht_corners[1][0][chn]=kernel.clone();
              kernel=extendFFTInputTo (kernels[tileCorners[1][1][1]][tileCorners[1][1][0]][chn], fftSize);
              fht_instance.swapQuadrants(kernel);
              fht_instance.transform(kernel);
              fht_corners[1][1][chn]=kernel.clone();
              fht_line= interpolateFHT (fht_corners[1][0][chn],    // first FHT array
                                        fht_corners[1][1][chn],    // second FHT array
                                                    points);   // array of interpolation points - 0.0 - fht0, 1.0 - fht1
              for (k=0;k<subTiles;k++) fhtCache[subTiles][k][chn]= fht_line[k];
            } else for (k=0;k<subTiles;k++) fhtCache[subTiles][k][chn]= null;
          }
/* Update what is cached */
          cachedTile[1][0][0]=tileCorners[1][0][0];
          cachedTile[1][0][1]=tileCorners[1][0][1];
          cachedTile[1][1][0]=tileCorners[1][1][0];
          cachedTile[1][1][1]=tileCorners[1][1][1];
/* Interpolate horizontally */
          for (chn=0;chn<nChn;chn++) {
            if ((compMask & (1<<chn))!=0)  {
              for (k=0;k<subTiles;k++) {// Interpolate rows
                fht_line= interpolateFHT (fhtCache[0][k][chn],    // first FHT array
                                   fhtCache[subTiles][k][chn],    // second FHT array
                                                       points);   // array of interpolation points - 0.0 - fht0, 1.0 - fht1
                for (l=1;l<subTiles;l++) fhtCache[l][k][chn]= fht_line[l];
              }
            } else  for (k=0;k<subTiles;k++) for (l=1;l<subTiles;l++) fhtCache[l][k][chn]=null;
          }

/* Convert cached kernels to FHT for convolution */
        } // end of at least some of the cache needs to be re-calculated
          else System.out.println("same tile");


//    double [][] test_arr;


if ((DEBUG_LEVEL>1) && (yTile==yTile0) && (xTile==xTile0)){
//  int test_index=0;
//  test_arr=new double[fhtCache.length*fhtCache[0].length][];
  for (k=0;k<fhtCache.length;k++) for (l=0;l<fhtCache[0].length;l++) {
     kernel=fhtCache[k][l][5].clone();
     fht_instance.inverseTransform(kernel);
     fht_instance.swapQuadrants(kernel);
     SDFA_instance.showArrays(kernel,  "t5-k"+k+"-l"+l);
  }
}



/* fhtCache array for this tile is current now, can be used for the convolution with the input image data [x][y]*/
        for (ySubTile=0;ySubTile<subTiles;ySubTile++) for (xSubTile=0;xSubTile<subTiles;xSubTile++) {
/* if the subtile intersects with the roi */
          rTile.x=xTile*kernelStep+xSubTile*imgStep-imgStep;
          rTile.y=yTile*kernelStep+ySubTile*imgStep-imgStep;
          if ((rTile.x>=(rroi.x-2*imgStep)) &&
              (rTile.y>=(rroi.y-2*imgStep)) &&
              (rTile.x< (rroi.x+rroi.width)) &&
              (rTile.y< (rroi.y+rroi.height))) {
            input_bayer= splitBayer (imp, rTile, equalize_greens);
            input_bayer= oversampleFFTInput (input_bayer,subpixel);
            if ((compMask & (1<< CHANNEL_CHECKER))!=0)  input_bayer=combineCheckerGreens (input_bayer,   // pixel arrays after oversampleFFTInput() or extractSimulPatterns())
                                                                              subpixel); // same as used in oversampleFFTInput() - oversampling ratio
//if ((DEBUG_LEVEL>1) && (yTile==yTile0) && (xTile==xTile0)){
//  SDFA_instance.showArrays(input_bayer,  "x"+xSubTile+ "y"+ySubTile);
//}

//            for (chn=0;chn<nChn;chn++) {
            for (chn=nChn-1;chn>=0;chn--) { /* start with green - before red and blue */
              if ((compMask & (1<<chn))==0) {
                input_bayer[chn]=null;
              } else {
/* Multiply expanded component tile by a sliding mask */
                for (l=0;l<input_bayer[chn].length;l++) input_bayer[chn][l]*=slidingMask[l];
/* extend input with zeros */
                input_bayer[chn]=extendFFTInputTo (input_bayer[chn], fftSize);
/* make direct FHT of the input image tile */
if ((DEBUG_LEVEL>2) && (yTile==yTile0) && (xTile==xTile0)){
  SDFA_instance.showArrays(input_bayer[chn],  "x"+xSubTile+ "y"+ySubTile+"-"+chn);
}
                fht_instance.swapQuadrants(input_bayer[chn]);
                fht_instance.transform(input_bayer[chn]);
                if (filterBayer) {
//  if ((DEBUG_LEVEL>2) && (yTile==yTile0) && (xTile==xTile0) && (ySubTile==0) && (xSubTile==0) ){
  if ((DEBUG_LEVEL>2) && (yTile==yTile0)  && (xTile==xTile0)){
//                 SDFA_instance.showArrays(input_bayer[chn],  "pre-"+chn);
                 SDFA_instance.showArrays(input_bayer[chn],  "pre-"+chn+"-y"+ySubTile+"-x"+xSubTile);
  }
                  if (chn== CHANNEL_CHECKER) { // should be before red/blue
    double debayer_threshold=0.0; // old implementation - will always use debayer

                    both_masks= aliasScissors(input_bayer[chn], // fht array for green, will be masked in-place
                                                      subpixel, // measured kernels are sampled at subpixel higher than Bayer (subpixel/2 - than sensor) //4
                                             debayer_threshold, // no high frequencies - use default uniform filter
//                                        debayer_relative_width, // Debayer lopass filter width (relative to distance to the nearest alias)
                                                 debayer_gamma, // power function applied to the amplitudes before generating spectral masks
                                                    debayer_rz, // for green mask - rays start radius from center, relative to distance to the first alias
                                                    debayer_ra, // for green mask - rays radius around aliases, relative to distance to the first alias
                                                 debayer_sigma, // for green mask - reduce value of the far amplitudes, relative to distance to the first alias
                                                 debayer_decay, // for green mask - exponential decay of the ray value, relative to distance to the first alias
                                          debayer_farthest_max, // for green mask - rays will be extyended from maf back to the center, max up to this distance (relative)
                                          debayer_radius_power, // for green mask -  divide ray values by the radius to this power
                                                   mainToAlias,// relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
                                             debayer_mask_blur, // for both masks  sigma for gaussian blur of the binary masks (<0 -do not use "scissors")
                                              debayer_lo_green, //  combine alias-reject "scissors" with lopass filter for greens
                                          debayer_lo_postgreen, // combine alias-reject "scissors" with lopass filter for greens applied to red/blue
                                            debayer_lo_redblue, // combine alias-reject "scissors" with lopass filter for greens applied to red/blue
                            (debayer_decay<0)?pol_instace:null,
                                                        lopass, // [1.0,scaled][green,redBlue][size*size] - externally prepared arrays, centered in the center
                                                            1); // internal debug level ((DEBUG_LEVEL>2) && (yTile==yTile0) && (xTile==xTile0))?3:1;
                  }


                  input_bayer[chn]=fht_instance.multiply(input_bayer[chn],both_masks[(chn==CHANNEL_CHECKER)?0:1],false);
//  if ((DEBUG_LEVEL>2) && (yTile==yTile0) && (xTile==xTile0) && (ySubTile==0) && (xSubTile==0) ){
  if ((DEBUG_LEVEL>2) && (yTile==yTile0)  && (xTile==xTile0)){
                 SDFA_instance.showArrays(input_bayer[chn],  "post-"+chn+"-y"+ySubTile+"-x"+xSubTile);
  }
                }

//test=input_bayer[chn].clone();
                product= fht_instance.multiply(input_bayer[chn], fhtCache[xSubTile][ySubTile][chn], false);
                fht_instance.inverseTransform(product);
                fht_instance.swapQuadrants(product);

//                if ((DEBUG_LEVEL>3)&&(tileY==0)&&(tileX==0)) SDFA_instance.showArrays(product, size, size, "product");
/* Add result to the output (float) array, compensate for red and blue Bayer shift */
                for (i=0;i<fftSize;i++) {
                  iy=(rTile.y-imgStep)*subpixel/2+i+((chn==2)?rb_shift:0); //product is twice larger (each direction) than rTile, the centers match
                  if ((iy>=0) && (iy<outHeight)) {
                    for (j=0;j<fftSize;j++) {
                      ix=(rTile.x-imgStep)*subpixel/2+j+((chn==1)?rb_shift:0);
                      if ((ix>=0) && (ix<outWidth)) {
                        result [chn][iy*outWidth+ix]+=(float) product[i*fftSize+j];
                      }
                    }
                  }
                }
/* Create "weights" by convoving same kernels with bayer pattern (not shifted, will shift when adding to the result ) */
                if (createBayerWeights) {
                  product= fht_instance.multiply(bayer_patterns[(chn==CHANNEL_CHECKER)?1:0], fhtCache[xSubTile][ySubTile][chn], false);
                  fht_instance.inverseTransform(product);
                  fht_instance.swapQuadrants(product);
/* Add result to the output array, compensate for Bayer shift for red and blue colors */
                  for (i=0;i<fftSize;i++) {
                    iy=(rTile.y-imgStep)*subpixel/2+i+((chn==2)?rb_shift:0); //product is twice larger (each direction) than rTile, the centers match
                    if ((iy>=0) && (iy<outHeight)) {
                      for (j=0;j<fftSize;j++) {
                        ix=(rTile.x-imgStep)*subpixel/2+j+((chn==1)?rb_shift:0);
                        if ((ix>=0) && (ix<outWidth)) {
                          bayerWeights [chn][iy*outWidth+ix]+=(float) product[i*fftSize+j];
                        }
                      }
                    }
                  }
                } // if (createBayerWeights)
              }
            }
//if ((DEBUG_LEVEL>1) && (yTile==yTile0) && (xTile==xTile0)){
//  SDFA_instance.showArrays(input_bayer,  "x"+xSubTile+ "y"+ySubTile);
//}
          }
        } // for (ySubTile=0;ySubTile<subTiles;ySubTile++) for (xSubTile=0;xSubTile<subTiles;xSubTile++)
      } //for (xTile=xTile0;xTile<=xTile1;xTile++)
    } // for (yTile=yTile0;yTile<=yTile1;yTile++) {
    double dR,dG,dB;
    if (createBayerWeights) { // calculate bayerWeightsRB
     for (i=0;i<bayerWeights[CHANNEL_GREEN].length;i++) {
       dR=bayerWeights[CHANNEL_RED][i];
       dG=bayerWeights[CHANNEL_GREEN][i];
       dB=bayerWeights[CHANNEL_BLUE][i];
       if (dR<0) dR=0;
       if (dG<0) dG=0;
       if (dB<0) dB=0;
       dG+=dR+dB;
       bayerWeights[CHANNEL_RED][i]= (float) (dR/dG);
       bayerWeights[CHANNEL_BLUE][i]=(float) (dB/dG);
     }
     bayerWeightsRB[0]=bayerWeights[CHANNEL_RED];
     bayerWeightsRB[1]=bayerWeights[CHANNEL_BLUE];
    }
    return result;
  }



/* old way - calculating amp/phases before interpolation */
  public float [][] convolveImageWithKernelsAmpPhases(ImagePlus imp, // 32-bit (float) linearized image with Bayer mosaic
                                               Rectangle r, // area to process (if not all the image) in pixels
                                   boolean equalize_greens, // equalize two green gains (currently - independently for each tile)
                                   double [][][][] kernels, // 2d array of convolution kernels 
//                            double [][][][][][] amplPhases, // 2d array of per component ampliteudes/full phases (for interpolation)
                                              int subpixel, // measured kernels are sampled at subpixel higher than Bayer (subpixel/2 - than sensor) //4
                                            int kernelStep, // distance between kernel samples, in sensor pixels // 64 // 32
                                         int kernelMargins, // number of kernel steps from the image 0,0 points to the center of the (0,0) kernel // 4
                                               int fftSize, // size of sliding FFT (twice the kernel size) here:128
                                              int compMask, // bitmask of color components to process
                                       boolean filterBayer, // compensate for Bayer aliases
/* Parameters related to reducing mosaic alias artifacts, can be removed if that filtering will be preformed before aberration correction.
     That just requires that lateral chtomatic aberration is small relative to FFT size, so spectrums amplitudes are similar */
                                      double debayer_gamma, // power function applied to the amplitudes before generating spectral masks
                                         double debayer_rz, // for green mask - rays start radius from center, relative to distance to the first alias
                                         double debayer_ra, // for green mask - rays radius around aliases, relative to distance to the first alias
                                      double debayer_sigma, // for green mask - reduce value of the far amplitudes, relative to distance to the first alias
                                      double debayer_decay, // for green mask - exponential decay of the ray value, relative to distance to the first alias
                               double debayer_farthest_max, // for green mask - rays will be extyended from maf back to the center, max up to this distance (relative)
                               double debayer_radius_power, // for green mask -  divide ray values by the radius to this power
                                        double mainToAlias, // relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
                                  double debayer_mask_blur, // for both masks  sigma for gaussian blur of the binary masks
                                  boolean debayer_lo_green, //  combine alias-reject "scissors" with lopass filter for greens
                              boolean debayer_lo_postgreen, // combine alias-reject "scissors" with lopass filter for greens applied to red/blue
                                boolean debayer_lo_redblue, // combine alias-reject "scissors" with lopass filter for greens applied to red/blue
/* end of demosaic parameters */
                                 float [][] bayerWeightsRB, // if non-null, should be a three element array, each will be assigned weight-red and weight_blue
                                      boolean updateStatus){  // update status info
    ImageProcessor ip=imp.getProcessor();
    Rectangle rroi=new Rectangle (r);
    if (rroi.x<0) rroi.x=0;
    if (rroi.y<0) rroi.y=0;
    rroi.x&=~1;
    rroi.x&=~1;
    if ((rroi.x+rroi.width )>ip.getWidth ()) rroi.width= ip.getWidth() -rroi.x;
    if ((rroi.y+rroi.height)>ip.getHeight()) rroi.height=ip.getHeight()-rroi.y;
    rroi.width &= ~1;
    rroi.height &= ~1;
    int i,j,k,l,chn;
    int nChn=0;
    for (i=compMask; i!=0;i>>=1) nChn++;
    float [][] result = new float [nChn][];
    int outWidth= ip.getWidth()* subpixel/2;
    int outHeight=ip.getHeight()*subpixel/2;

    for (i=0;i<nChn;i++) if ((compMask & (1<<i))!=0) {
      result[i]=new float [outWidth*outHeight];
      for (j=0;j<result[i].length;j++) result[i][j]=(float) 0.0;
    }

    int imgStep=fftSize/4/(subpixel/2); // Tile overlap in sensor pixels (16)cosMask
    int xTile0=rroi.x/kernelStep;
    int yTile0=rroi.y/kernelStep;
    int xTile1=(rroi.x+rroi.width-1)/kernelStep+1;
    int yTile1=(rroi.y+rroi.height-1)/kernelStep+1;
    int xTile,yTile,xSubTile,ySubTile,iy,ix;

//    double [] slidingMask=getSlidingMask(kernelStep*subpixel); // 32*4=128
    int slidingMaskSize=imgStep*subpixel; // 16*4=64
    double [] slidingMask=getSlidingMask(slidingMaskSize); //(imgStep*subpixel); // 16*4=64

    double [][] bayer_patterns=null;
/* Later computetion could be greately decresed by at least using half-size fft, or just convolving directly (only 4x4=16 points are needed) */
    boolean hasRGB=((compMask & (1<< CHANNEL_RED))!=0) && ((compMask & (1<< CHANNEL_BLUE))!=0) && ((compMask & (1<< CHANNEL_GREEN))!=0); // has channels 1,2, and 5
    boolean createBayerWeights=(bayerWeightsRB!=null) && (bayerWeightsRB.length>=2) && hasRGB;
    float [][] bayerWeights=null; // same first dimension as result[][], in the end will be re-calculated to a two-=element bayerWeightsRB[2][]
    if (createBayerWeights) {
      bayerWeights = new float [nChn][];
      for (i=0;i<nChn;i++) if ((compMask & (1<<i))!=0) {
        bayerWeights[i]=new float [outWidth*outHeight];
        for (j=0;j<bayerWeights[i].length;j++) bayerWeights[i][j]=(float) 0.0;
      }
      for (i=0;i<2;i++); for (j=0;j<bayerWeights[i].length;j++) bayerWeights[i][j]=(float) 0.0;
      bayer_patterns=new double [2][slidingMaskSize*slidingMaskSize];
      for (iy=0;iy<slidingMaskSize;iy++) for (ix=0;ix<slidingMaskSize;ix++) {
        bayer_patterns[0][iy*slidingMaskSize+ix]=((( iy % subpixel)==0) && ((ix % subpixel)==0 ))?1.0:0.0; // all colors but checker
        bayer_patterns[1][iy*slidingMaskSize+ix]=((((iy % subpixel)==0) && ((ix % subpixel)==0 )) || (((iy % subpixel)==(subpixel/2)) && ((ix % subpixel)==(subpixel/2) )))?1.0:0.0; // checker greens
      }
      bayer_patterns=extendFFTInputTo (bayer_patterns, fftSize);
/* prepare bayer paterns for convolution - perform FHT transform */
      for (i=0;i<bayer_patterns.length;i++) {
                fht_instance.swapQuadrants(bayer_patterns[i]);
                fht_instance.transform(bayer_patterns[i]);
      }
    }
    int rb_shift=subpixel/2; //half Bayer period
    double [][] input_bayer;
    Rectangle rTile=new Rectangle (2*imgStep,2*imgStep); // 64x64
    int subTiles=kernelStep/imgStep; // 4 
    double [][][][][][] amplPhasesCache=new double[subTiles+1][subTiles][][][][]; // first dimension HORIZONTAL here
    double [][][][]     fhtCache=       new double[subTiles][subTiles][nChn][];      // first dimension HORIZONTAL here
//    boolean interpolated=false; // outside of the margins use closest kernel, do not interpolate
    int [][][] cachedTile=  new int [2][2][2];
    int [][][] tileCorners= new int [2][2][2]; // which kernels are needed for the 4 corners of the current tile
    double [][][][][][] amplPhases=new double [2][2][nChn][][][]; // kernels, extended and converted to amplitudes/full phases fro the 4 corners of the tile
                                                                  // first index horizontal
    double [] kernel;
    double [] product;
//    double [] test;
    if (DEBUG_LEVEL>1) System.out.println("convolveImageWithKernels(): yTile0="+yTile0+" xTile0="+xTile0+" yTile1="+yTile1+" xTile1="+xTile1+" compMask="+compMask);
if (DEBUG_LEVEL>2){
  SDFA_instance.showArrays(slidingMask,  "slidingMask");
}

    double debayer_relative_width_green=1.0; // realtive width of the low-pass filter for greens (1.0 goes to zero at exactrly the closest alias )
    double debayer_relative_width_redblue=1.0;// same for individual (re and blue) color components
    double [][][] lopass= createAliasFilters (debayer_relative_width_green, // realtive width of the low-pass filter for greens (1.0 goes to zero at exactrly the closest alias )
                                            debayer_relative_width_redblue, // same for individual (re and blue) color components
                                                                   fftSize, // side of the square
                                                                  subpixel); // should be 4 now


    double [][] both_masks=null; /* green and r/b "scissors" masks to reduce sampling aliases */
    PolarSpectrums pol_instace=new PolarSpectrums(
                fftSize, // size of the square array, centar is at size/2, size/2, only top half+line will be used
                 Math.PI, //2*Math.PI, // i.e. Math.PI, 2*Math.PI
                 fftSize/2-2, // width of the polar array - should be <= size/2-2
                 0.5, //0.75, //2.0, //0.5, // maximal step in pixels on the maxRadius for 1 angular step (i.e. 0.5)
                 4 );// angular symmetry - 0- none,1 - pi corresponds to integer, 2 - pi/2 corresponds to integer, n - pi/n corresponds to integer angular step 

    for (i=0;i<2;i++) for (j=0;j<2;j++)  for (k=0;k<2;k++) cachedTile[i][j][k]=-1;
    for (yTile=yTile0;yTile<=yTile1;yTile++) {
//      if (updateStatus) IJ.showStatus("Convolving image, row"+(yTile-yTile0+1)+" of "+(yTile1-yTile0+1));
//      rTile.y=(yTile-1)*kernelStep;
      for (xTile=xTile0;xTile<=xTile1;xTile++) {
        if (updateStatus) IJ.showStatus("Convolving image, tile "+((yTile-yTile0)*(xTile1-xTile0+1)+(xTile-xTile0)+1)+" of "+((yTile1-yTile0+1)*(xTile1-xTile0+1)));
    if (DEBUG_LEVEL>1) System.out.println("convolveImageWithKernels(): yTile="+yTile+" xTile="+xTile);
        for (i=0;i<2;i++) for (j=0;j<2;j++) {
          tileCorners[i][j][0]=xTile-kernelMargins+j;
          tileCorners[i][j][1]=yTile-kernelMargins+i;
/* fix margins */
          if (tileCorners[i][j][0]<0) tileCorners[i][j][0]=0;
          if (tileCorners[i][j][1]<0) tileCorners[i][j][1]=0;
          if (tileCorners[i][j][0]>=kernels[0].length) tileCorners[i][j][0]=kernels[0].length-1;
          if (tileCorners[i][j][1]>=kernels.length)    tileCorners[i][j][1]=kernels.length-1;
        }
if (DEBUG_LEVEL>1) System.out.println(" tileCorners[0][0]="+tileCorners[0][0][0]+"/"+tileCorners[0][0][1]+
                                      " tileCorners[0][1]="+tileCorners[0][1][0]+"/"+tileCorners[0][1][1]+
                                      " tileCorners[1][0]="+tileCorners[1][0][0]+"/"+tileCorners[1][0][1]+
                                      " tileCorners[1][1]="+tileCorners[1][1][0]+"/"+tileCorners[1][1][1]);
if (DEBUG_LEVEL>1) System.out.println(" cachedTile[0][0]=" +cachedTile[0][0][0]+"/"+ cachedTile[0][0][1]+
                                      " cachedTile[0][1]=" +cachedTile[0][1][0]+"/"+ cachedTile[0][1][1]+
                                      " cachedTile[1][0]=" +cachedTile[1][0][0]+"/"+ cachedTile[1][0][1]+
                                      " cachedTile[1][1]=" +cachedTile[1][1][0]+"/"+ cachedTile[1][1][1]);

/* not dealing with missing kernels here - that can be done before */
/* compare new corners with the cache */
        if ((tileCorners[0][0][0]!=cachedTile[0][0][0]) ||
            (tileCorners[0][0][1]!=cachedTile[0][0][1]) ||
            (tileCorners[0][1][0]!=cachedTile[0][1][0]) ||
            (tileCorners[0][1][1]!=cachedTile[0][1][1]) ||
            (tileCorners[1][0][0]!=cachedTile[1][0][0]) ||
            (tileCorners[1][0][1]!=cachedTile[1][0][1]) ||
            (tileCorners[1][1][0]!=cachedTile[1][1][0]) ||
            (tileCorners[1][1][1]!=cachedTile[1][1][1])) { // at least some of the cache needs to be re-calculated
/* Is it next tile in a row? */
          if ((tileCorners[0][0][0]==cachedTile[1][0][0]) &&
              (tileCorners[0][0][1]==cachedTile[1][0][1]) &&
              (tileCorners[0][1][0]==cachedTile[1][1][0]) &&
              (tileCorners[0][1][1]==cachedTile[1][1][1])) { /**yes, just copy last interpolated column to the first */
            amplPhasesCache[0]=amplPhasesCache[subTiles];
          } else { /* no, need to recalculate first column too */
            amplPhasesCache[0]=new double[subTiles][nChn][][][];
            for (chn=0;chn<nChn;chn++) {
              if ((compMask & (1<<chn))!=0)  {
    if (DEBUG_LEVEL>1) System.out.println("tileCorners[0][0][1]="+tileCorners[0][0][1]+" tileCorners[0][0][0]="+tileCorners[0][0][0]+ " chn="+chn);

                kernel=extendFFTInputTo (kernels[tileCorners[0][0][1]][tileCorners[0][0][0]][chn], fftSize);
                fht_instance.swapQuadrants(kernel);
                fht_instance.transform(kernel);
/*
if ((DEBUG_LEVEL>2) && (yTile==yTile0) && (xTile==xTile0)){
  test=kernel.clone();
  fht_instance.inverseTransform(test);
  fht_instance.swapQuadrants(test);
  SDFA_instance.showArrays(test,  "test-"+chn);
}
*/

                amplPhases[0][0][chn]=amplPhase(FHT2FFTHalf (kernel, fftSize)); 
                kernel=extendFFTInputTo (kernels[tileCorners[0][1][1]][tileCorners[0][1][0]][chn], fftSize);
                fht_instance.swapQuadrants(kernel);
                fht_instance.transform(kernel);
                amplPhases[0][1][chn]=amplPhase(FHT2FFTHalf (kernel, fftSize)); 
                amplPhasesCache[0][0][chn]=amplPhases[0][0][chn]; // copy top left corner
                for (k=1;k<subTiles;k++) {// Interpolate the rest of the column
                  amplPhasesCache[0][k][chn]= interpolateAmplPhase(amplPhases[0][0][chn], // first set of amplitudes/phases
                                                                   amplPhases[0][1][chn], // second set of amplitudes/phases
                                                                       (1.0*k)/subTiles); // interpolation coefficient
                }
              } else for (k=0;k<subTiles;k++) amplPhasesCache[0][k][chn]=null;
            }
          }
/* Update what is cached */
          cachedTile[0][0][0]=tileCorners[0][0][0];
          cachedTile[0][0][1]=tileCorners[0][0][1];
          cachedTile[0][1][0]=tileCorners[0][1][0];
          cachedTile[0][1][1]=tileCorners[0][1][1];
/* interpolate the last column */
          amplPhasesCache[subTiles]=new double[subTiles][nChn][][][];
          for (chn=0;chn<nChn;chn++) {
            if ((compMask & (1<<chn))!=0)  {
              kernel=extendFFTInputTo (kernels[tileCorners[1][0][1]][tileCorners[1][0][0]][chn], fftSize);
              fht_instance.swapQuadrants(kernel);
              fht_instance.transform(kernel);
              amplPhases[1][0][chn]=amplPhase(FHT2FFTHalf (kernel, fftSize)); 
              kernel=extendFFTInputTo (kernels[tileCorners[1][1][1]][tileCorners[1][1][0]][chn], fftSize);
              fht_instance.swapQuadrants(kernel);
              fht_instance.transform(kernel);
              amplPhases[1][1][chn]=amplPhase(FHT2FFTHalf (kernel, fftSize)); 
              amplPhasesCache[subTiles][0][chn]=amplPhases[1][0][chn]; // copy top right corner
              for (k=1;k<subTiles;k++) {// Interpolate the rest of the column
                amplPhasesCache[subTiles][k][chn]= interpolateAmplPhase(amplPhases[1][0][chn], // first set of amplitudes/phases
                                                                        amplPhases[1][1][chn], // second set of amplitudes/phases
                                                                            (1.0*k)/subTiles); // interpolation coefficient
              }
            } else  for (k=0;k<subTiles;k++) amplPhasesCache[subTiles][k][chn]=null;
          }
/* Update what is cached */
          cachedTile[1][0][0]=tileCorners[1][0][0];
          cachedTile[1][0][1]=tileCorners[1][0][1];
          cachedTile[1][1][0]=tileCorners[1][1][0];
          cachedTile[1][1][1]=tileCorners[1][1][1];
/* Interpolate horizontally */
          for (l=1;l<subTiles;l++) {
            amplPhasesCache[l]=new double[subTiles][nChn][][][];
            for (chn=0;chn<nChn;chn++) {
              if ((compMask & (1<<chn))!=0)  {
                for (k=0;k<subTiles;k++) {// Interpolate rows
                  amplPhasesCache[l][k][chn]= interpolateAmplPhase(amplPhasesCache[0]       [k][chn], // first set of amplitudes/phases
                                                                   amplPhasesCache[subTiles][k][chn], // second set of amplitudes/phases
                                                                                   (1.0*l)/subTiles); // interpolation coefficient
                }
              } else  for (k=0;k<subTiles;k++) amplPhasesCache[l][k][chn]=null; //null pointer
            }
          }
/*
if ((DEBUG_LEVEL>2) && (yTile==yTile0) && (xTile==xTile0) ){
                   System.out.println("* cachedTile[0][0]=" +cachedTile[0][0][0]+"/"+ cachedTile[0][0][1]+
                                      "* cachedTile[0][1]=" +cachedTile[0][1][0]+"/"+ cachedTile[0][1][1]+
                                      "* cachedTile[1][0]=" +cachedTile[1][0][0]+"/"+ cachedTile[1][0][1]+
                                      "* cachedTile[1][1]=" +cachedTile[1][1][0]+"/"+ cachedTile[1][1][1]);
  int testX,testY;
  double [][][] ap;
  for (testY=0;testY<2;testY++) for (testX=0;testX<2;testX++) {
              kernel=extendFFTInputTo (kernels[cachedTile[testX][testY][1]][cachedTile[testX][testY][0]][1], fftSize);
              fht_instance.swapQuadrants(kernel);
              fht_instance.transform(kernel);
              ap=amplPhase(FHT2FFTHalf (kernel, fftSize));
              test=ampPhaseToFHT(ap);
              fht_instance.inverseTransform(test);
              fht_instance.swapQuadrants(test);
              SDFA_instance.showArrays(test,  "cornerY"+testY+"X"+testX);
  }
}
*/


/* Convert cached kernels to FHT for convolution */
          for (k=0;k<subTiles;k++) for (l=0;l<subTiles;l++) for (chn=0;chn<nChn;chn++){
            if ((compMask & (1<<chn))!=0)  {
              fhtCache[k][l][chn]= ampPhaseToFHT(amplPhasesCache[l][k][chn]);
/*
if ((DEBUG_LEVEL>2) && (yTile==yTile0) && (xTile==xTile0) && (chn==1)){
  test=fhtCache[k][l][chn].clone();
  fht_instance.inverseTransform(test);
  fht_instance.swapQuadrants(test);
  SDFA_instance.showArrays(test,  "ik"+k+ ":"+l+"-"+chn);
}
*/
            } else fhtCache[k][l][chn]=null;
          }
        } // end of at least some of the cache needs to be re-calculated
/* fhtCache array for this tile is current now, can be used for the convolution with the input image data */
        for (ySubTile=0;ySubTile<subTiles;ySubTile++) for (xSubTile=0;xSubTile<subTiles;xSubTile++) {
/* if the subtile intersects with the roi */
          rTile.x=xTile*kernelStep+xSubTile*imgStep-imgStep;
          rTile.y=yTile*kernelStep+ySubTile*imgStep-imgStep;
          if ((rTile.x>=(rroi.x-2*imgStep)) &&
              (rTile.y>=(rroi.y-2*imgStep)) &&
              (rTile.x< (rroi.x+rroi.width)) &&
              (rTile.y< (rroi.y+rroi.height))) {
            input_bayer= splitBayer (imp, rTile, equalize_greens);
            input_bayer= oversampleFFTInput (input_bayer,subpixel);
            if ((compMask & (1<< CHANNEL_CHECKER))!=0)  input_bayer=combineCheckerGreens (input_bayer,   // pixel arrays after oversampleFFTInput() or extractSimulPatterns())
                                                                              subpixel); // same as used in oversampleFFTInput() - oversampling ratio
//if ((DEBUG_LEVEL>1) && (yTile==yTile0) && (xTile==xTile0)){
//  SDFA_instance.showArrays(input_bayer,  "x"+xSubTile+ "y"+ySubTile);
//}

//            for (chn=0;chn<nChn;chn++) {
            for (chn=nChn-1;chn>=0;chn--) { /* start with green - before red and blue */
              if ((compMask & (1<<chn))==0) {
                input_bayer[chn]=null;
              } else {
/* Multiply expanded component tile by a sliding mask */
                for (l=0;l<input_bayer[chn].length;l++) input_bayer[chn][l]*=slidingMask[l];
/* extend input with zeros */
                input_bayer[chn]=extendFFTInputTo (input_bayer[chn], fftSize);
/* make direct FHT of the input image tile */
if ((DEBUG_LEVEL>2) && (yTile==yTile0) && (xTile==xTile0)){
  SDFA_instance.showArrays(input_bayer[chn],  "x"+xSubTile+ "y"+ySubTile+"-"+chn);
}
                fht_instance.swapQuadrants(input_bayer[chn]);
                fht_instance.transform(input_bayer[chn]);
                if (filterBayer) {
//  if ((DEBUG_LEVEL>2) && (yTile==yTile0) && (xTile==xTile0) && (ySubTile==0) && (xSubTile==0) ){
  if ((DEBUG_LEVEL>2) && (yTile==yTile0)  && (xTile==xTile0)){
//                 SDFA_instance.showArrays(input_bayer[chn],  "pre-"+chn);
                 SDFA_instance.showArrays(input_bayer[chn],  "pre-"+chn+"-y"+ySubTile+"-x"+xSubTile);
  }
      double debayer_threshold=0.0; // old way
                  if (chn== CHANNEL_CHECKER) { // shouold be before red/blue
                    both_masks= aliasScissors(input_bayer[chn], // fht array for green, will be masked in-place
                                                      subpixel, // measured kernels are sampled at subpixel higher than Bayer (subpixel/2 - than sensor) //4
                                             debayer_threshold, // no high frequencies - use default uniform filter
                                                 debayer_gamma, // power function applied to the amplitudes before generating spectral masks
                                                    debayer_rz, // for green mask - rays start radius from center, relative to distance to the first alias
                                                    debayer_ra, // for green mask - rays radius around aliases, relative to distance to the first alias
                                                 debayer_sigma, // for green mask - reduce value of the far amplitudes, relative to distance to the first alias
                                                 debayer_decay, // for green mask - exponential decay of the ray value, relative to distance to the first alias
                                          debayer_farthest_max, // for green mask - rays will be extyended from maf back to the center, max up to this distance (relative)
                                          debayer_radius_power, // for green mask -  divide ray values by the radius to this power
                                                   mainToAlias, // relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
                                             debayer_mask_blur, // for both masks  sigma for gaussian blur of the binary masks (<0 -do not use "scissors")
                                              debayer_lo_green, //  combine alias-reject "scissors" with lopass filter for greens
                                          debayer_lo_postgreen, // combine alias-reject "scissors" with lopass filter for greens applied to red/blue
                                            debayer_lo_redblue, // combine alias-reject "scissors" with lopass filter for greens applied to red/blue
                            (debayer_decay<0)?pol_instace:null,
                                                        lopass,
                                                            1); // internal debug level ((DEBUG_LEVEL>2) && (yTile==yTile0) && (xTile==xTile0))?3:1;
                  }
                  input_bayer[chn]=fht_instance.multiply(input_bayer[chn],both_masks[(chn==CHANNEL_CHECKER)?0:1],false);

//  if ((DEBUG_LEVEL>2) && (yTile==yTile0) && (xTile==xTile0) && (ySubTile==0) && (xSubTile==0) ){
  if ((DEBUG_LEVEL>2) && (yTile==yTile0)  && (xTile==xTile0)){
                 SDFA_instance.showArrays(input_bayer[chn],  "post-"+chn+"-y"+ySubTile+"-x"+xSubTile);
  }
                }

//test=input_bayer[chn].clone();
                product= fht_instance.multiply(input_bayer[chn], fhtCache[xSubTile][ySubTile][chn], false);
                fht_instance.inverseTransform(product);
                fht_instance.swapQuadrants(product);
/*
if ((DEBUG_LEVEL>2) && (yTile==yTile0) && (xTile==xTile0)){
  SDFA_instance.showArrays(product,  "rslt-x"+xSubTile+ "y"+ySubTile+"-"+chn);
  test=fhtCache[xSubTile][ySubTile][chn].clone();
  fht_instance.inverseTransform(test);
  fht_instance.swapQuadrants(test);
  SDFA_instance.showArrays(test,  "kernel-x"+xSubTile+ "y"+ySubTile+"-"+chn);
}
*/

//                if ((DEBUG_LEVEL>3)&&(tileY==0)&&(tileX==0)) SDFA_instance.showArrays(product, size, size, "product");
/* Add result to the output (float) array, compensate for red and blue Bayer shift */
                for (i=0;i<fftSize;i++) {
                  iy=(rTile.y-imgStep)*subpixel/2+i+((chn==2)?rb_shift:0); //product is twice larger (each direction) than rTile, the centers match
                  if ((iy>=0) && (iy<outHeight)) {
                    for (j=0;j<fftSize;j++) {
                      ix=(rTile.x-imgStep)*subpixel/2+j+((chn==1)?rb_shift:0);
                      if ((ix>=0) && (ix<outWidth)) {
                        result [chn][iy*outWidth+ix]+=(float) product[i*fftSize+j];
                      }
                    }
                  }
                }
/* Create "weights" by convoving same kernels with bayer pattern (not shifted, will shift when adding to the result ) */
                if (createBayerWeights) {
                  product= fht_instance.multiply(bayer_patterns[(chn==CHANNEL_CHECKER)?1:0], fhtCache[xSubTile][ySubTile][chn], false);
                  fht_instance.inverseTransform(product);
                  fht_instance.swapQuadrants(product);
/* Add result to the output array, compensate for Bayer shift for red and blue colors */
                  for (i=0;i<fftSize;i++) {
                    iy=(rTile.y-imgStep)*subpixel/2+i+((chn==2)?rb_shift:0); //product is twice larger (each direction) than rTile, the centers match
                    if ((iy>=0) && (iy<outHeight)) {
                      for (j=0;j<fftSize;j++) {
                        ix=(rTile.x-imgStep)*subpixel/2+j+((chn==1)?rb_shift:0);
                        if ((ix>=0) && (ix<outWidth)) {
                          bayerWeights [chn][iy*outWidth+ix]+=(float) product[i*fftSize+j];
                        }
                      }
                    }
                  }
                } // if (createBayerWeights)
              }
            }
//if ((DEBUG_LEVEL>1) && (yTile==yTile0) && (xTile==xTile0)){
//  SDFA_instance.showArrays(input_bayer,  "x"+xSubTile+ "y"+ySubTile);
//}
          }
        } // for (ySubTile=0;ySubTile<subTiles;ySubTile++) for (xSubTile=0;xSubTile<subTiles;xSubTile++)
      } //for (xTile=xTile0;xTile<=xTile1;xTile++)
    } // for (yTile=yTile0;yTile<=yTile1;yTile++) {
    double dR,dG,dB;
    if (createBayerWeights) { // calculate bayerWeightsRB
     for (i=0;i<bayerWeights[CHANNEL_GREEN].length;i++) {
       dR=bayerWeights[CHANNEL_RED][i];
       dG=bayerWeights[CHANNEL_GREEN][i];
       dB=bayerWeights[CHANNEL_BLUE][i];
       if (dR<0) dR=0;
       if (dG<0) dG=0;
       if (dB<0) dB=0;
       dG+=dR+dB;
       bayerWeights[CHANNEL_RED][i]= (float) (dR/dG);
       bayerWeights[CHANNEL_BLUE][i]=(float) (dB/dG);
     }
     bayerWeightsRB[0]=bayerWeights[CHANNEL_RED];
     bayerWeightsRB[1]=bayerWeights[CHANNEL_BLUE];
    }
    return result;
  }

/* subtract bayer aliases, FHT array should be made from the Bayer starting from pixel at (0,0) - either single color or checker pattern */ 
  double [][][] createAliasFilters (double debayer_relative_width_green, // realtive width of the low-pass filter for greens (1.0 goes to zero at exactrly the closest alias )
                                  double debayer_relative_width_redblue, // same for individual (re and blue) color components
                                                               int size, // side of the square
                                                           int subpixel){ // should be 4 now
    int i;
    double [] cosMask= createCosMask (size,  subpixel);   // oversampling
    double [][] [] lopass =new double [2][2][];
//System.out.println("size="+size+", subpixel="+subpixel);
//    for (i=0;i<cosMask.length;i++) System.out.println("cosMask["+i+"]="+cosMask[i]);
    lopass[0][0]=new double [size*size];
    for (i=0;i<lopass[0][0].length;i++) lopass[0][0][i]=1.0;
    lopass[0][1]=lopass[0][0].clone();
    lopass[1][0]=lopass[0][0].clone();
    lopass[1][1]=lopass[0][0].clone();
    maskBayerAliases (lopass[0][0],   // FHT array to be filtered
                           cosMask,   // cosine mask array
                              true); // this fht array is for the checkerboard greens
    maskBayerAliases (lopass[0][1],   // FHT array to be filtered
                           cosMask,   // cosine mask array
                             false); // this fht array is for the checkerboard greens
   
 // [1.0,scaled][green,redBlue][size*size] - externally prepared arrays, centered in the center
    cosMask= createCosMask ((int) Math.round(size*debayer_relative_width_green),  subpixel);   // oversampling
    maskBayerAliases (lopass[1][0],   // FHT array to be filtered
                           cosMask,   // cosine mask array
                              true); // this fht array is for the checkerboard greens
    cosMask= createCosMask ((int) Math.round(size*debayer_relative_width_redblue),  subpixel);   // oversampling
    maskBayerAliases (lopass[1][1],   // FHT array to be filtered
                           cosMask,   // cosine mask array
                             false); // this fht array is for the checkerboard greens
    fht_instance.swapQuadrants(lopass[0][0]);
    fht_instance.swapQuadrants(lopass[0][1]);
    fht_instance.swapQuadrants(lopass[1][0]);
    fht_instance.swapQuadrants(lopass[1][1]);

    if (DEBUG_LEVEL>1) {
      SDFA_instance.showArrays(lopass[0], size,size, "narrow");
      SDFA_instance.showArrays(lopass[1], size,size, "wide");
    }
    return lopass;
  }




  void maskBayerAliases (double [] fht,   // FHT array to be filtered
//                                int subdiv,   // oversampling
                         double [] cosMask,   // cosine mask array
                         boolean isChecker) { // this fht array is for the checkerboard greens
     int size= (int) Math.sqrt(fht.length);
     int iy,ix, ix1,iy1;
//     int tsize= size*(isChecker?2:1)/subdiv;
     int tsize= (cosMask.length-1)/(isChecker?1:2);
     int index=0;
     int hsizeM1=(size/2)-1;
//    System.out.println("askBayerAliases(),tsize="+tsize+" size="+size);
//    for (ix=0;ix<cosMask.length;ix++) System.out.println("cosMask["+ix+"]="+cosMask[ix]);
//    for (ix=0;ix<fht.length;ix++) System.out.println("fht["+ix+"]="+fht[ix]);
     if (isChecker) {

       for (iy=0;iy<size;iy++) {
         iy1=(iy+hsizeM1)%size -hsizeM1;
         for (ix=0;ix<size;ix++) {
           ix1=(ix+hsizeM1)%size -hsizeM1;
           if (((ix1+iy1)>-tsize) &&
               ((ix1-iy1)>-tsize) &&
               ((ix1+iy1)<=tsize) &&
               ((ix1-iy1)<=tsize)) fht[index++]*=cosMask[Math.abs(ix1+iy1)]*cosMask[Math.abs(ix1-iy1)];
           else fht[index++]=0.0;
         }
       }

     } else {
       for (iy=0;iy<size;iy++) {
         iy1=(iy+hsizeM1)%size -hsizeM1;
         for (ix=0;ix<size;ix++) {
           ix1=(ix+hsizeM1)%size -hsizeM1;
           if ((iy1>-tsize) && (iy1<=tsize) && (ix1>-tsize) && (ix1<=tsize)) fht[index++]*=cosMask[2*Math.abs(iy1)]*cosMask[2*Math.abs(ix1)];
           else fht[index++]=0.0;
         }
       }
     }
  }

  double [] createCosMask (int fftsize,   // FHT array to be filtered - just length is used
                            int subdiv   // oversampling
                                      ) { // this fht array is for the checkerboard greens
     int size= 2*fftsize/subdiv;
     double [] cosMask=new double [size+1];
     for (int i=0;i<=size;i++) cosMask[i]=0.5*(1.0+Math.cos(i*Math.PI/size));
     return cosMask; 
  }

  public void aliasesRejectInKernels (double [][][][] rkernels,  // 2-d array of reversed psf kernels
                                       double sigma_individual,  // sigma for individual bayer components (1/4 pixels))
                                         double sigma_diagonal,  // sigma for diagonal greens (square array rotated 45 degrees) (1/2 pixels))
                                          double sigma_checker,  // sigma for checkerboard  greens (normal array, half pixels 0)
                                                  int subpixel,  // measured array is sampled at subpixel higher than Bayer (subpixel/2 - than sensor)
                                         boolean updateStatus){  // update status info
   int i,j,k;
   int tilesY=rkernels.length;
   int tilesX=rkernels[0].length;
   int nChn;
   double sigma=sigma_individual;
//   double [] smoothKernel;
//   double [] centroidXY;
//   double [] variableSigmas;
   double [][]aliasRejectMask={null,null,null};
   int size;
   int aliasMaskType=0;
   for (i=0;i<tilesY;i++) for (j=0;j<tilesX;j++) {
     if (rkernels[i][j]!=null) {
       nChn=rkernels[i][j].length;
       if (updateStatus) IJ.showStatus("Removing sampling aliases, tile "+(i*tilesY+j+1)+" of "+(tilesY*tilesX));
       for (k=0;k<nChn;k++) if ( rkernels[i][j][k]!=null) {
         switch (k) {
            case 0:
            case 1:
            case 2:
            case 3:
               sigma=sigma_individual;
               aliasMaskType=0;
               break;
            case 4:
               sigma=sigma_diagonal;
               aliasMaskType=1;
               break;
            case 5:
               aliasMaskType=2;
               sigma=sigma_checker;
               break;
          }
          size= (int) Math.sqrt(rkernels[i][j][k].length);
          if (aliasRejectMask[aliasMaskType]==null) { /* Create this alias mask */
             aliasRejectMask[aliasMaskType]= createAliasReject (size,  // size of the mask
                                                  (aliasMaskType==2),  // checkerboard pattern in the source file (use when filtering)
                                                            subpixel,  // measured array is sampled at 1/oversample frequency than model (will add more parameters later)
                                                               sigma); // width of rejection areas on the spectrum (the smaller, the sharper rejection)
          }
          rkernels[i][j][k]=rejectByMask (rkernels[i][j][k],  // square input data
                             aliasRejectMask[aliasMaskType], // mask to multiply FHT
                                                       true); // image is centered around the center of the square (use swapQuadrants)
/* normalize result kernels */
          normalizeKernel(rkernels[i][j][k]); // in-place
       }
     } 
   }
  }

/**TODO: Modify - divide second by first and use the result for phases. When calculating full phase - use product of the amplitudes */

  public double [][][][] interpolateKernels (double [][][][] kernels,  // 2d array of component kernels (some components may be null)
                                                          int subdiv, // number of subdivisions (i.e. subdiv=2 will insert interpolated kernel between each 2)
                                                boolean updateStatus){  // update status info
     int origWidth= kernels[0].length;
     int origHeight=kernels.length;
     int width=(origWidth-1)*subdiv+1;
     int height=(origHeight-1)*subdiv+1;
     double [][][][][][] amplPhases=new double[height][width][][][][];
     double [] kernel;
     double [][][]fft;
     int i,j,k,size,l;
     double subdivFrac=1.0/subdiv;
     for (i=0;i<height;i++) for (j=0;j<width;j++) amplPhases[i][j]=null;
/* calculate amplitudes/phases for input kernels */
     for (i=0;i<origHeight;i++) {
       if (updateStatus) IJ.showStatus("Calculating amplitudes/phases, row "+(i+1)+" of "+origHeight);
       for (j=0;j<origWidth;j++) {
         if (kernels[i][j]!=null) {
           amplPhases[i*subdiv][j*subdiv]= new double [kernels[i][j].length][][][];
           for (k=0;k<kernels[i][j].length;k++) {
             if (kernels[i][j][k]==null) amplPhases[i*subdiv][j*subdiv][k]=null;
             else {
               kernel=kernels[i][j][k].clone();
               size=(int)Math.sqrt(kernel.length);
               fht_instance.swapQuadrants(kernel);
               fht_instance.transform(kernel);
               fft=  FHT2FFTHalf (kernel, size);
               amplPhases[i*subdiv][j*subdiv][k]=amplPhase(fft); 
             }
           }
         }
       }
     }
/* create bi-linear interpolation, first fill columns by interpoalting vertically*/
     for (i=0;i<(origHeight-1);i++) {
       if (updateStatus) IJ.showStatus("Interpolating vertically, row "+(i+1)+" of "+(origHeight-1));
       for (j=0;j<origWidth;j++) if ((amplPhases[i*subdiv][j*subdiv]!=null) && (amplPhases[(i+1)*subdiv][j*subdiv]!=null)){
         for (l=1;l<subdiv;l++) {
           amplPhases[i*subdiv+l][j*subdiv]=new double [amplPhases[i*subdiv][j*subdiv].length][][][];
           for (k=0;k<amplPhases[i*subdiv][j*subdiv].length;k++) {
              if ((amplPhases[i*subdiv][j*subdiv][k]==null) || (amplPhases[(i+1)*subdiv][j*subdiv][k]==null)) {
                amplPhases[i*subdiv+l][j*subdiv][k]=null;
              } else {
//      if (DEBUG_LEVEL>1) {
//                 System.out.println("interpolateKernels(): i="+i+" j="+j+" l="+l+" k="+k+" i*subdiv="+(i*subdiv)+" j*subdiv="+(j*subdiv));
//      }  

                amplPhases[i*subdiv+l][j*subdiv][k]=interpolateAmplPhase(amplPhases[i*    subdiv][j*subdiv][k], // first set of amplitudes/phases
                                                                         amplPhases[(i+1)*subdiv][j*subdiv][k], // second set of amplitudes/phases
                                                                                                 l*subdivFrac); // interpolation coefficient
              }
           }
         }
       }
     }
/* now interpolate all the remaining kernels */
     for (i=0;i<height;i++) {
       if (updateStatus) IJ.showStatus("Interpolating horizontally, row "+(i+1)+" of "+height);
       for (j=0;j<(origWidth-1);j++) if ((amplPhases[i][j*subdiv]!=null) && (amplPhases[i][(j+1)*subdiv]!=null)){
         for (l=1;l<subdiv;l++) {
           amplPhases[i][j*subdiv+l]=new double [amplPhases[i][j*subdiv].length][][][];
           for (k=0;k<amplPhases[i][j*subdiv].length;k++) {
             if ((amplPhases[i][j*subdiv][k]==null) || (amplPhases[i][(j+1)*subdiv][k]==null)) {
               amplPhases[i][j*subdiv+l][k]=null;
             } else {
               amplPhases[i][j*subdiv+l][k]=interpolateAmplPhase(amplPhases[i][j*    subdiv][k], // first set of amplitudes/phases
                                                                 amplPhases[i][(j+1)*subdiv][k], // second set of amplitudes/phases
                                                                                 l*subdivFrac); // interpolation coefficient
             }
           }
         }
       }
     }
/* Convert each amplitude/phase to space */
     double [][][][] allKernels=new double[height][width][][];
     for (i=0;i<height;i++) {
       if (updateStatus) IJ.showStatus("converting to space, row "+(i+1)+" of "+height);
       for (j=0;j<width;j++) {
         if (amplPhases[i][j]==null) allKernels[i][j]=null;
         else {
           allKernels[i][j]=new double [amplPhases[i][j].length][];
           for (k=0;k<allKernels[i][j].length;k++) {
             if (amplPhases[i][j][k]==null) allKernels[i][j][k]=null;
             else allKernels[i][j][k]=ampPhaseToSpace(amplPhases[i][j][k], true);
           }
         }
       }
     }
     return allKernels;
  }
/* TODO: use direct kernels (mirrored) for the centers*/
  public void variableBlurDeconvolutionKernels (double [][][][] rkernels,  // 2-d array of reversed psf kernels
                                                 double sigma_individual,  // sigma for individual bayer components (1/4 pixels))
                                                   double sigma_diagonal,  // sigma for diagonal greens (square array rotated 45 degrees) (1/2 pixels))
                                                    double sigma_checker,  // sigma for checkerboard  greens (normal array, half pixels 0)
                                                      double sigma_scale,  // scale sigma in the center when using variable sigma
                                                    double sigmaToRadius,  // sigma-to-radius ratio (0.0 to disable variable blur)
                                                    boolean updateStatus){  // update status info
   int i,j,k;
   int tilesY=rkernels.length;
   int tilesX=rkernels[0].length;
   int nChn;
   double sigma=sigma_individual;
   double [] smoothKernel;
   double [] centroidXY;
   double [] variableSigmas;
   int size;
   for (i=0;i<tilesY;i++){
     if (updateStatus) IJ.showStatus("Filtering deconvolution kernels, row "+(i+1)+" of "+tilesY);
     for (j=0;j<tilesX;j++) {
       if (rkernels[i][j]!=null) {
         nChn=rkernels[i][j].length;
//         if (updateStatus) IJ.showStatus("Filtering deconvolution kernels, tile "+(i*tilesY+j+1)+" of "+(tilesY*tilesX));
         for (k=0;k<nChn;k++) if ( rkernels[i][j][k]!=null) {
            switch (k) {
              case 0:
              case 1:
              case 2:
              case 3:
                 sigma=sigma_individual;
                 break;
              case 4:
                 sigma=sigma_diagonal;
                 break;
              case 5:
                 sigma=sigma_checker;
                 break;
            }
            smoothKernel=lowPassGauss(rkernels[i][j][k], 2*sigma,true);
            size= (int) Math.sqrt(rkernels[i][j][k].length);
            if (sigmaToRadius>0.0) { // recalculate smoothDeconvKernels[i]
              centroidXY=extractCentroidFromReverseKernel(smoothKernel,   // square array of reversed PSF
                                                                   0.5); // fraction of the maximal value to use as a bottom of the part used for centroid calculation
              variableSigmas= createSigmasFromCenter((int) Math.sqrt(rkernels[i][j][k].length), // side of square
                                                                                 sigmaToRadius, // variable blurring - sigma will be proportional distance from the center
                                                                             sigma*sigma_scale, //blurring in the center sigma(r)=sqrt((sigma_to_radius*r)^2)+center_sigma^2)
                                                                                 centroidXY[0], // coordinates of the center (0:0 - size/2: size/2)
                                                                                 centroidXY[1]);
              smoothKernel=variableGaussBlurr(rkernels[i][j][k], // input square pixel array, preferrably having many exact zeros (they will be skipped)
                                                 variableSigmas, // array of sigmas to be used for each pixel, matches pixels[]
                                                            3.5, // drop calculatin if farther then nSigma
                                                              0, // int WOICenterX, // window of interest in pixels[] array - do not generate data outside it
                                                              0, // int WOICenterY, // 
                                                           size, //int WOIWidth, reduce later
                                                           size); //int WOIHeight)
            }
            rkernels[i][j][k]=smoothKernel;
/* normalize deconvolve kernels */
            normalizeKernel(rkernels[i][j][k]); // in-place
         }
       } 
     }
   }
}




/* Generates gaussian kernels from the deconvolution ones, preserving center (lateral chromatic) */
/* seems not working - need to blur more before looking for a maximum - try direct kernels and mirroring center*/

  public double [][][][] generateGaussianKernels (double [][][][] rkernels,  // 2-d array of reversed psf kernels
                                                   double sigma_individual,  // sigma for individual bayer components (1/4 pixels))
                                                     double sigma_diagonal,  // sigma for diagonal greens (square array rotated 45 degrees) (1/2 pixels))
                                                      double sigma_checker,  // sigma for checkerboard  greens (normal array, half pixels 0)
                                                      boolean updateStatus){ // update status info
   int i,j,k;
   int tilesY=rkernels.length;
   int tilesX=rkernels[0].length;
   int nChn;
   double [][][][] gaussians= new double [tilesY][tilesX][][];
   double sigma=sigma_individual;
   double [] smoothKernel;
//   double [] centroidXY;
//   double [] variableSigmas;
   for (i=0;i<tilesY;i++) for (j=0;j<tilesX;j++) {
     if (rkernels[i][j]==null) gaussians[i][j]=null;
     else {
       nChn=rkernels[i][j].length;
       gaussians[i][j]=new double[nChn][];
       if (updateStatus) IJ.showStatus("Creating gaussian kernels, tile "+(i*tilesY+j+1)+" of "+(tilesY*tilesX));
       for (k=0;k<nChn;k++) if ( rkernels[i][j][k]!=null) {
         switch (k) {
            case 0:
            case 1:
            case 2:
            case 3:
               sigma=sigma_individual;
               break;
            case 4:
               sigma=sigma_diagonal;
               break;
            case 5:
               sigma=sigma_checker;
               break;
          }
          smoothKernel=lowPassGauss(rkernels[i][j][k], 2*sigma,true);
//          centroidXY=extractCentroidFromReverseKernel(smoothKernel,   // square array of reversed PSF
//                                                               0.5); // fraction of the maximal value to use as a bottom of the part used for centroid calculation
          gaussians[i][j][k]=extractLateralChromaticFromReverseKernel(smoothKernel,sigma,0.5);
       }
     } 
   }
   return gaussians;
}

  public double [][][][] generateGaussianKernelsFromDirect (double [][][][] kernels,  // 2-d array of direct psf kernels
                                                                   int gaussianSize,
                                                            double sigma_individual,  // sigma for individual bayer components (1/4 pixels))
                                                              double sigma_diagonal,  // sigma for diagonal greens (square array rotated 45 degrees) (1/2 pixels))
                                                               double sigma_checker,  // sigma for checkerboard  greens (normal array, half pixels 0)
                                                               boolean updateStatus){ // update status info
   int i,j,k;
   int tilesY=kernels.length;
   int tilesX=kernels[0].length;
   int nChn;
   double [][][][] gaussians= new double [tilesY][tilesX][][];
   double sigma=sigma_individual;
   double [] smoothKernel;
//   double [] centroidXY;
//   double [] variableSigmas;
   for (i=0;i<tilesY;i++) for (j=0;j<tilesX;j++) {
     if (kernels[i][j]==null) gaussians[i][j]=null;
     else {
       nChn=kernels[i][j].length;
       gaussians[i][j]=new double[nChn][];
       if (updateStatus) IJ.showStatus("Creating gaussian kernels, tile "+(i*tilesY+j+1)+" of "+(tilesY*tilesX));
       for (k=0;k<nChn;k++) if ( kernels[i][j][k]!=null) {
         switch (k) {
            case 0:
            case 1:
            case 2:
            case 3:
               sigma=sigma_individual/2;
               break;
            case 4:
               sigma=sigma_diagonal/2;
               break;
            case 5:
               sigma=sigma_checker/2;
               break;
          }
          smoothKernel=lowPassGauss(kernels[i][j][k], 2*sigma,true);

/* ????? */
/*          centroidXY=extractCentroidFromReverseKernel(smoothKernel,   // square array of reversed PSF
                                                               0.5); // fraction of the maximal value to use as a bottom of the part used for centroid calculation */
          gaussians[i][j][k]=extractLateralChromaticFromDirectKernel(smoothKernel,gaussianSize,sigma,0.5);
       }
     } 
   }
   return gaussians;
}






  public double [] createSigmasFromCenter(int size, // side of square
                            double sigma_to_radius, // variable blurring - sigma will be proportional distance from the center
                               double center_sigma, //blurring in the center sigma(r)=sqrt((sigma_to_radius*r)^2)+center_sigma^2)
                                    double centerX, // coordinates of the center (0:0 - size/2: size/2)
                                    double centerY) {
    double [] sigmas = new double [size*size];
    int i,j;
    double x,y;
    double center_sigma2=center_sigma*center_sigma;
    double sigma_to_radius2=sigma_to_radius*sigma_to_radius;
    for (i=0;i<size;i++) for (j=0;j<size;j++) {
      y=i-size/2-centerY;
      x=j-size/2-centerX;
      sigmas[i*size+j]=Math.sqrt((x*x+y*y)*sigma_to_radius2+ center_sigma2);
    }
    return sigmas;
  }




/*
double [] psf2mtf (double [] psf, // square array of psf pixels
                boolean centered, // true if PSF center is in the center of the array, false if it is at point 0
               boolean normalize ) { // normalize so mtf(0)=1.0
*/

public double [][][][] calculateMTFfromPSF(double [][][][] PSFKernels, // 2-d array of direct psf kernels
                                                             int size // size (side of square) of reverse PSF kernel
                                                                   ){
   int i,j,k;
   int tilesY=PSFKernels.length;
   int tilesX=PSFKernels[0].length;
   int nChn;
   double [][][][] MTFs= new double [tilesY][tilesX][][];
   for (i=0;i<tilesY;i++) for (j=0;j<tilesX;j++) {
     if (PSFKernels[i][j]==null) MTFs[i][j]=null;
     else {
       nChn=PSFKernels[i][j].length;
       MTFs[i][j]=resizeForFFT(PSFKernelMap[i][j],size);
       for (k=0;k<nChn;k++) if (MTFs[i][j][k]!=null) {
          MTFs[i][j][k]= psf2mtf (MTFs[i][j][k], // square array of psf pixels
                                           true, // true if PSF center is in the center of the array, false if it is at point 0
                                           true);  // normalize so mtf(0)=1.0

       }
     } 

   }
   return MTFs;
}


public double [][][][] reversePSFKernels(double [][][][] PSFKernels, // 2-d array of direct psf kernels
                                                           int size, // size (side of square) of reverse PSF kernel
                                                 double invertRange,    // deconvInvert
                                           double otf_cutoff_energy,  // OTF_cutoff_energy
                                           double otf_ellipse_scale,  // ellipse mask size relative to the cluster
                                          boolean otf_ellipse_gauss,  // use Gauss instead of Hamming for ellipse mask
                                           double psf_cutoff_energy,  // OTF_cutoff_energy
                                           double psf_ellipse_scale,  // ellipse mask size relative to the cluster
                                     double rpsf_min_mask_threshold,  // zero output element if elliptical Gauss mask is below this threshold
                                               boolean updateStatus   // update status info
                                                                   ){
   int i,j,k;
   int tilesY=PSFKernels.length;
   int tilesX=PSFKernels[0].length;
   int nChn;
   double [][][][] rPSFKernels= new double [tilesY][tilesX][][];
   for (i=0;i<tilesY;i++){
     if (updateStatus) IJ.showStatus("Reversing PSF, row "+(i+1)+" of "+tilesY);
     for (j=0;j<tilesX;j++) {
       if (PSFKernels[i][j]==null) rPSFKernels[i][j]=null;
       else {
//         if (updateStatus) IJ.showStatus("Reversing PSF, tile "+(i*tilesY+j+1)+" of "+(tilesY*tilesX));
         nChn=PSFKernels[i][j].length;
         rPSFKernels[i][j]=resizeForFFT(PSFKernelMap[i][j],size);
         for (k=0;k<nChn;k++) if (rPSFKernels[i][j][k]!=null) {
/* reverse PSF kernel */
/* TODO: convert cleanupAndReversePSF() to double FHT*/
          rPSFKernels[i][j][k]= cleanupAndReversePSF (rPSFKernels[i][j][k],  // input pixels
                                                               invertRange,    // deconvInvert
                                                         otf_cutoff_energy,  // OTF_cutoff_energy
                                                         otf_ellipse_scale,  // ellipse mask size relative to the cluster
                                                         otf_ellipse_gauss,  // use Gauss instead of Hamming for ellipse mask
                                                                         1,  // decimate frequency to match Bayer component pixels pitch
                                                                     false,  // fold high frequency into low, when false - use Hamming to cut off high frequencies
                                                                        ""); // just for the plot names
/* Find direct kernel approximation ellipse, increase it, mirror center around 0,0 and use it as a mask for the reversed kernel */
          rPSFKernels[i][j][k]=maskReversePSFKernel( PSFKernels[i][j][k], // direct PSF function, square array, may be proportionally larger than reversed
                                                    rPSFKernels[i][j][k], // reversed psf, square array
                                                       psf_cutoff_energy, // fraction of energy in the pixels to be used
                                                       psf_ellipse_scale,
                                                 rpsf_min_mask_threshold, // zero output element if elliptical Gauss mask is below this threshold
                                                                     "");
          normalizeKernel(rPSFKernels[i][j]); // in-place
/* Add filtering of high-frequeny far from 0 */



         }
       } 
     }
   }
   return rPSFKernels;
}




  public int [] channelNumbers (boolean [] activeChannels) {
    int nChn=0;
    int i;
    for (i=0;i<activeChannels.length;i++) if (activeChannels[i]) nChn++;
    int [] channels = new int [nChn];
    nChn=0;
    for (i=0;i<activeChannels.length;i++) if (activeChannels[i]) channels[nChn++]=i;
    return channels;
  }


/*************************************************************************************/

  public double [][] getPSFKernels ( ImagePlus imp,
                                          int size,   // size in pixels (twice FFTSize)
                                            int x0,      // top left corner X (pixels)
                                            int y0,      // top left corner Y (pixels)
                            boolean equalizeGreens,
                                      double gamma, // gamma to use for power spectrum for correlation
                                      double sigma, // high-pass gaussian filter sigma when correlating power spectrum
                              int diff_spectr_corr, // maximal distance between maximum on spectrum and predicted maximum on autocorrelation of gamma(|spectrum|)
                            double shrink_clusters, // Shrink clusters by this ratio (remove lowest) after initial separation
                              int multiples_to_try, // try this number of maximums proportionally farther from 0,0 than the two closest (increase precision)
                                  double deviation, // when looking for maximums - maximal distance from predicted from the lower order one
                               int deviation_steps, // maximal iterations when looking for local maximum
                               boolean [] bPattern, // pattern bitmask (calculate once)
                                        int subdiv, // simulatin boolean pixels per sensor pixel
                               boolean centerForG2,
                                    int subDivFreq, // add zeros arond (increase frequency resolution, probably will not be used here)
                                 double [] Hamming, //=initHamming( fft_size) calculate once
                             double [] fullHamming, //=initHamming( fft_size*subpixel);
                                      int subpixel, // use finer grid than actual pixels 
                              double zerofreq_size,
                                  double simulfill, // part of the (center) pixel area being "photosensitive"
                        boolean [] colorsToProcess, // color channnels to process (should be at least 6 long)
                               double deconvInvert, //  fraction of the maximal value to be used to limit zeros
// next 3 used for filtering aliases 
                                   double smoothPS, // 0 - none, otherwise Gauss width = FFT size/2/smoothPS
                             double threshold_high,  // reject completely if energy is above this part of maximal
                              double threshold_low,  // leave intact if energy is below this part of maximal
// Next 6 were used to filter out X-shaped artifacts on the PSF (when used plain checkerboard pattern), probbaly not needed anymore
                                 double X_threshold, // if 0 use normalize amplitude, if 0..1 - make binary: 1.0 if > threshold, 0.0 - otherwise -1 - disable mask
                                    double X_radius, // low-pass result with low pass filter (should be later defined automatically)
                                       int  X_hsize,  // 2d histogram size (size/2 probably a good guess),
                                double X_percentile, // use this percentile (0.0..1.0)) value for given radius as a target
                                   double X_maxGain, // maximal gain for low components
                                double X_exaggerate, // exaggerate correction mask with Math.pow()) 
                                  int referenceComp, // number of color component to reference lateral chromatic aberration to (now 4 - checkerboard greens)
                                boolean useHamming,  // if false - use gausian
                               double   minContrast,    // minimal instance contrast to use in binning
                                double   windowFrac,     // reduce the PSF cell size to this part of the area connecting first negative clones
                                   boolean  symm180,
                           boolean  ignoreChromatic,
                        boolean  enableModelSubtract,  // generate/subtract gaussian models (not needed if no overlap between pos/neg clones)
                             double   smoothSeparate, // low-pass filter multiple opposite-sign PSF instaces for separation (width relative to distance to the opposite sign PSF)
                          double   thresholdSeparate,  // threshold for locating zero-crossing
                                  double   topCenter,     // consider only points above that fraction of the local max to find centroid
                            boolean  removeNegtative,  // remove PSF negative values  when separating composite PSF (will need low-pass filtering)
                              double   sigmaToRadius,   // 0.4;  variable-sigma blurring to reduce high frequencies more for the pixels farther from the PSF center
                               double   wings_energy, // fraction of energy in the pixels to be used
                        double   wings_ellipse_scale
                                   ){
      if (imp==null) return null; // Maybe convert to double pixel array once to make it faster?
      String title=imp.getTitle()+"X"+x0+"Y"+y0;
      Rectangle PSFCell=new Rectangle (x0,y0,size,size);
      int fft_size=size/2;
      double [][] input_bayer=splitBayer (imp,PSFCell,equalizeGreens);
//      int greensToProcess=4;
      int i,j,l;
/* Calculate pattern parameters, including distortion */
      double[][] distortedPattern= findPatternDistorted(input_bayer, // pixel array to process (no widowing!)
                                                              gamma,
                                                              sigma,  // pattern detection: high-pass filter (0.0 - none) gamma(PS)
                                                   diff_spectr_corr,
                                                    shrink_clusters,
                                                   multiples_to_try,
                                                          deviation,
                                                    deviation_steps,
                                                              true, //(greensToProcess==4), // boolean greens, // this is a pattern for combined greens (diagonal), adjust results accordingly
                                                              title); // title prefix to use for debug  images
      if (distortedPattern==null) return null;

      boolean[][] simulation_barray=  simulatePatternFullPattern(bPattern,
                                                   distortedPattern[0][0],
                                                   distortedPattern[0][1],
                                                   distortedPattern[0][2],
                                                   distortedPattern[1][0],
                                                   distortedPattern[1][1],
                                                   distortedPattern[1][2],
                                                      distortedPattern[2], //
                                                                   subdiv,
                                                                 fft_size,
                                                              centerForG2);
//      for (i=0;i<4;i++) if (!colorsToProcess[i]) input_bayer[i]=null; // leave composite greens even if disabled
      input_bayer= normalizeAndWindow (input_bayer, Hamming);
      if (subDivFreq>1)  input_bayer= extendFFTInput (input_bayer,subDivFreq);
      if ((subpixel>1) && (zerofreq_size>=0)) {
        input_bayer= oversampleFFTInput (input_bayer,subpixel);
        if (colorsToProcess[5])  input_bayer=combineCheckerGreens (input_bayer,   // pixel arrays after oversampleFFTInput() or extractSimulPatterns())
                                                                   subpixel); // same as used in oversampleFFTInput() - oversampling ratio
      }
      for (i=0;i<4;i++) if (!colorsToProcess[i]) input_bayer[i]=null; // leave composite greens even if disabled
      if (DEBUG_LEVEL>3) {
        if (zerofreq_size>=0) SDFA_instance.showArrays(input_bayer, fft_size*subDivFreq*subpixel, fft_size*subDivFreq*subpixel, title);
        else                      SDFA_instance.showArrays(input_bayer, fft_size*subDivFreq, fft_size*subDivFreq, title);
      }
      double [][] simul_pixels= extractSimulPatterns (simulation_barray, // high resolution boolean pattern array
                                                              simulfill, // part of the (center) pixel area being "photosensitive"
                                                                 subdiv, // boolean pixels to real pixels resolution
                                                               subpixel, // subdivide pixels
                                                      fft_size*subpixel, // number of Bayer cells in width of the square selection (half number of pixels)
                                                                    0.0,    // selection center, X (in pixels)
                                                                    0.0);   // selection center, y (in pixels)
      simul_pixels= normalizeAndWindow (simul_pixels, fullHamming);
      if (subDivFreq>1) {
          simul_pixels= extendFFTInput (simul_pixels,subDivFreq);
      }
      if ((subpixel>1) && (zerofreq_size>=0)) {
                      if (colorsToProcess[5])  simul_pixels=combineCheckerGreens (simul_pixels,   // pixel arrays after oversampleFFTInput() or extractSimulPatterns())
                                                                          subpixel); // same as used in oversampleFFTInput() - oversampling ratio
      }
      if (DEBUG_LEVEL>2) SDFA_instance.showArrays(simul_pixels, fft_size*subDivFreq*subpixel, fft_size*subDivFreq*subpixel, title+"-SIM");
      if (DEBUG_LEVEL>2) System.out.println ( " input_bayer.length="+input_bayer.length+" simul_pixels.length="+simul_pixels.length);

      double averagePeriod=1.0/Math.sqrt((distortedPattern[0][0]*distortedPattern[0][0]+distortedPattern[0][1]*distortedPattern[0][1]+
                                          distortedPattern[1][0]*distortedPattern[1][0]+distortedPattern[1][1]*distortedPattern[1][1])/2.0);

      double [][] inverted=new double[colorsToProcess.length][];

      for (i=0;(i<input_bayer.length) && (i<simul_pixels.length);i++) if ((colorsToProcess[i]) && (input_bayer[i]!=null)){
        if (DEBUG_LEVEL>2) System.out.println ( "Color "+componentColorNames[i]+" XMAX_radius ("+XMASK_radius+") is re-calculated into bayer pixels as "+ ((i==4)?Math.sqrt(2):1.0)*0.5*averagePeriod*XMASK_radius);
        if (zerofreq_size<0) {
            inverted[i]=limitedInverseOfFHTDiffSize(input_bayer[i],
                                         simul_pixels[i],
                                            deconvInvert,
                                                    true, //      forwardOTF,
                                            title+"-"+i);
        } else {
          inverted[i]=limitedInverseOfFHT(input_bayer[i],
                                         simul_pixels[i],
                         fft_size*subDivFreq*subpixel,
                                                  (i==5),     //    boolean checker // checkerboard pattern in the source file (use when filtering)
                                            deconvInvert,
                                                    true, //      forwardOTF,
                                                subpixel,
                                           zerofreq_size,
                                                smoothPS,
                                          threshold_high,
                                           threshold_low,
                                             X_threshold,
    ((i==4)?Math.sqrt(2):1.0)*0.5*averagePeriod*X_radius,
                                                 X_hsize,
                                            X_percentile,
                                               X_maxGain,
                                            X_exaggerate,
                                             title+"-"+i);

        }
      }
      if (DEBUG_LEVEL>1) SDFA_instance.showArrays(inverted, fft_size*subDivFreq*subpixel, fft_size*subDivFreq*subpixel, title+"_Combined-PSF");
/* correct composite greens */
/* Here we divide wave vectors by subpixel as the pixels are already added */
      double [][]wVectors=  {{2.0*distortedPattern[0][0]/subpixel,                     2.0*distortedPattern[0][1]/subpixel}, // pattern was in sensor pixels, not bayer sub-pixel
                             {2.0*distortedPattern[1][0]/subpixel,                     2.0*distortedPattern[1][1]/subpixel}};
      double [][] wVrotMatrix= {{0.5,0.5},{-0.5,0.5}};
      double [][]wVectors4= new double [2][2];
      for (i=0;i<2;i++) for (j=0;j<2;j++) {
         wVectors4[i][j]=0.0;
         for (l=0;l<2;l++) wVectors4[i][j]+=wVectors[i][l]*wVrotMatrix[l][j];
      }
      double [][] PSF_shifts=         new double [input_bayer.length][];    // X/Y shift of the PSF array, in Bayer component pioxel coordinates (same as PSF arrays)
      double [][] PSF_centroids=      new double [input_bayer.length][];    // X/Y coordinates of the centroids of PSF in Bayer component pioxel coordinates (same as PSF arrays) (after they were optionally shifted)
      double [][] lateralChromatic=   new double [input_bayer.length][]; // X/Y coordinates of the centroids of Bayer component PSF in sensor pixel coordinates
      double [][] kernelsForFFT=      new double [input_bayer.length][];
      double [][] psf_inverted=       new double [input_bayer.length][];
      double [][] psf_inverted_masked=new double [input_bayer.length][];
      double [] lateralChromaticAbs=new double [input_bayer.length];
      double [] zeroVector={0.0,0.0};
      for (i=input_bayer.length-1;i>=0;i--) {
        if (colorsToProcess[i]) {
          PSF_shifts[i]=       zeroVector.clone();
          PSF_centroids[i]=    zeroVector.clone();
          lateralChromatic[i]= zeroVector.clone();
        } else {
          PSF_shifts[i]=       null;
          PSF_centroids[i]=    null;
          lateralChromatic[i]= null;
        }
        lateralChromaticAbs[i]=0.0;
        kernelsForFFT[i]=null;
        psf_inverted[i]=null;
        psf_inverted_masked[i]=null;
      }
//      int [][]  clusterMask;
/* Start with referenceComp */
      i= referenceComp;
      if (DEBUG_LEVEL>2) {
                 System.out.println("Before: color Component "+i+" PSF_shifts["+i+"][0]="+IJ.d2s(PSF_shifts[i][0],3)+
                                                " PSF_shifts["+i+"][1]="+IJ.d2s(PSF_shifts[i][1],3));
      }  

      kernels[i]=combinePSF (inverted[i], // Square array of pixels with multiple repeated PSF (alternating sign)
               (i==4)?wVectors4:wVectors, // two wave vectors, lengths in cycles/pixel (pixels match pixel array)
                                       1, // already applied subpixel,         // sub-pixel decimation 
                             minContrast,   // minimal instance contrast to use in binning
                              windowFrac,    // reduce the PSF cell size to this part of the area connecting first negative clones
                              useHamming,
                                 symm180, //PSF_symm180,
                                    true, //PSF_ignoreChromatic
                           PSF_shifts[i],  // centerXY[] - will be modified inside combinePSF() if PSF_ignoreChromatic is true
                         PSF_centroids[i], // will return array of XY coordinates of the result centroid
                     enableModelSubtract,  // generate/subtract gaussian models (not needed if no overlap between pos/neg clones)
                          smoothSeparate,  // low-pass filter multiple opposite-sign PSF instaces for separation (width relative to distance to the opposite sign PSF)
                       thresholdSeparate,  // threshold for locating zero-crossing
                               topCenter,  // consider only points above this fraction of the peak to find the centroid
                         removeNegtative,  // remove PSF negative values  when separating composite PSF (will need low-pass filtering)
                           sigmaToRadius,   // 0.4;  variable-sigma blurring to reduce high frequencies more for the pixels farther from the PSF center
                            wings_energy, // fraction of energy in the pixels to be used
                     wings_ellipse_scale,
                             title+"_"+i,
                         (DEBUG_LEVEL>4));
      if (DEBUG_LEVEL>2)     System.out.println("After-1: color Component "+i+"    PSF_shifts["+i+"][0]="+IJ.d2s(PSF_shifts   [i][0],3)+"    PSF_shifts["+i+"][1]="+IJ.d2s(   PSF_shifts[i][1],3));
      if (DEBUG_LEVEL>2)     System.out.println("After-1: color Component "+i+" PSF_centroids["+i+"][0]="+IJ.d2s(PSF_centroids[i][0],3)+" PSF_centroids["+i+"][1]="+IJ.d2s(PSF_centroids[i][1],3));

      if ((referenceComp==4) && !PSF_ignoreChromatic) { /* Recalculate center to pixels from greens (diagonal)) and supply it to other colors (lateral chromatic aberration correction) */
        for (j=0;j<input_bayer.length;j++) if ((colorsToProcess[j]) && (j!=referenceComp)) {
          PSF_shifts[j]=shiftSensorToBayer (shiftBayerToSensor(PSF_shifts[i],4,subpixel),j,subpixel);
          if (DEBUG_LEVEL>2)       System.out.println("After-2 (recalc): color Component "+j+" PSF_shifts["+j+"][0]="+IJ.d2s(PSF_shifts[j][0],3)+" PSF_shifts["+j+"][1]="+IJ.d2s(PSF_shifts[j][1],3));
        }
      }
      lateralChromatic[i]=shiftBayerToSensor ( PSF_shifts[i][0]+PSF_centroids[i][0],
                                               PSF_shifts[i][1]+PSF_centroids[i][1],
                                               i,
                                               subpixel);
      lateralChromaticAbs[i]=Math.sqrt(lateralChromatic[i][0]*lateralChromatic[i][0]+lateralChromatic[i][1]*lateralChromatic[i][1]);
/* Now process all the other components */
      for (i=0; i<input_bayer.length;i++) if ((i!=referenceComp) && (colorsToProcess[i])) {
        if (DEBUG_LEVEL>2) {
                   System.out.println("Before: color Component "+i+" PSF_shifts["+i+"][0]="+IJ.d2s(PSF_shifts[i][0],3)+
                                                  " PSF_shifts["+i+"][1]="+IJ.d2s(PSF_shifts[i][1],3));
        }  
        kernels[i]=combinePSF (inverted[i], // Square array of pixels with multiple repeated PSF (alternating sign)
                 (i==4)?wVectors4:wVectors, // two wave vectors, lengths in cycles/pixel (pixels match pixel array)
                                         1, // already applied subpixel,         // sub-pixel decimation 
                               minContrast,   // minimal instance contrast to use in binning
                                windowFrac,    // reduce the PSF cell size to this part of the area connecting first negative clones
                                useHamming,
                                   symm180, //PSF_symm180,
                           ignoreChromatic, //PSF_ignoreChromatic
                             PSF_shifts[i],  // centerXY[] - will be modified inside combinePSF() if PSF_ignoreChromatic is true
                           PSF_centroids[i], // will return array of XY coordinates of the result centroid
                       enableModelSubtract,  // generate/subtract gaussian models (not needed if no overlap between pos/neg clones)
                            smoothSeparate,  // low-pass filter multiple opposite-sign PSF instaces for separation (width relative to distance to the opposite sign PSF)
                         thresholdSeparate,  // threshold for locating zero-crossing
                                 topCenter,  // consider only points above this fraction of the peak to find the centroid
                           removeNegtative,  // remove PSF negative values  when separating composite PSF (will need low-pass filtering)
                             sigmaToRadius,   // 0.4;  variable-sigma blurring to reduce high frequencies more for the pixels farther from the PSF center
                              wings_energy, // fraction of energy in the pixels to be used
                       wings_ellipse_scale,
                               title+"_"+i,
                            (DEBUG_LEVEL>4));
        if (DEBUG_LEVEL>2)     System.out.println("After-1: color Component "+i+"    PSF_shifts["+i+"][0]="+IJ.d2s(PSF_shifts   [i][0],3)+"    PSF_shifts["+i+"][1]="+IJ.d2s(   PSF_shifts[i][1],3));
        if (DEBUG_LEVEL>2)     System.out.println("After-1: color Component "+i+" PSF_centroids["+i+"][0]="+IJ.d2s(PSF_centroids[i][0],3)+" PSF_centroids["+i+"][1]="+IJ.d2s(PSF_centroids[i][1],3));
        lateralChromatic[i]=shiftBayerToSensor ( PSF_shifts[i][0]+PSF_centroids[i][0],
                                                 PSF_shifts[i][1]+PSF_centroids[i][1],
                                                                                    i,
                                                                          subpixel);
        lateralChromaticAbs[i]=Math.sqrt((lateralChromatic[i][0]-lateralChromatic[referenceComp][0])*(lateralChromatic[i][0]-lateralChromatic[referenceComp][0])+
                                         (lateralChromatic[i][1]-lateralChromatic[referenceComp][1])*(lateralChromatic[i][1]-lateralChromatic[referenceComp][1]));
      }
      if (DEBUG_LEVEL>1) {
        for (i=0;i<PSF_shifts.length;i++) if (colorsToProcess[i]){
          if (DEBUG_LEVEL>2) {
            System.out.println("Color Component "+i+" subpixel="+subpixel+
                                                    " PSF_ignoreChromatic="+PSF_ignoreChromatic+
                                                    " PSF_symm180="+PSF_symm180);
            System.out.println(                     " PSF_shifts["+i+"][0]="+IJ.d2s(PSF_shifts[i][0],3)+
                                                    " PSF_shifts["+i+"][1]="+IJ.d2s(PSF_shifts[i][1],3)+
                                                    " PSF_centroids["+i+"][0]="+IJ.d2s(PSF_centroids[i][0],3)+
                                                    " PSF_centroids["+i+"][1]="+IJ.d2s(PSF_centroids[i][1],3));
            System.out.println("  lateralChromatic["+i+"][0]="+IJ.d2s(lateralChromatic[i][0],3)+
                               "  lateralChromatic["+i+"][1]="+IJ.d2s(lateralChromatic[i][1],3));
          }
        }
        if (colorsToProcess[referenceComp]) for (i=0;i<colorsToProcess.length;i++) if ((colorsToProcess[i])&& (i!=referenceComp)){
          System.out.println(componentColorNames[i]+" lateral chromatic (from green) "+IJ.d2s(lateralChromaticAbs[i],3)+"pix:  ["+i+"][0]="+IJ.d2s(lateralChromatic[i][0]-lateralChromatic[referenceComp][0],3)+
                            "  ["+i+"][1]="+IJ.d2s(lateralChromatic[i][1]-lateralChromatic[referenceComp][1],3));
        }
        System.out.println("Lateral shift green from simulation "+IJ.d2s(lateralChromaticAbs[referenceComp],3)+"pix:  ["+referenceComp+"][0]="+IJ.d2s(lateralChromatic[referenceComp][0],3)+
                                                                                                                    "  ["+referenceComp+"][1]="+IJ.d2s(lateralChromatic[referenceComp][1],3));
      }
      return kernels;
  }








/*************************************************************************************/


public double [] createNonlinMask(double [] hipassPixels,  // image with PSF deconvolution applied
                                  double [] lopassPixels,  // image convolved with shiftted gaussian (reject subsampling modulation, compensate lateral chromatic aberration)
                                        double maskSigma,  // gaussian sigma to create mask for selecting between restored/noisy and original/smooth
                                           double maskLo,  // minimal value of the low-pass filtered squared difference where mask >0.0
                                           double maskHi, // maximal value of the low-pass filtered squared difference where mask <1.0
                                        double threshold){
  if ((hipassPixels==null) || (lopassPixels==null)) return null;
  int length=hipassPixels.length;
  if (lopassPixels.length!=length) return null;
  int size=(int) Math.sqrt(length);
  double [] mask=new double [length];
  int i;
  double d;
  double max=0.0;
  for (i=0;i<length;i++) {
    d=hipassPixels[i]-lopassPixels[i];
    mask[i]=d*d;
    if (max<lopassPixels[i]) max=lopassPixels[i];
  }
  max*=threshold;
  for (i=0;i<length;i++) {
    mask[i]/=Math.max(max,lopassPixels[i]);
  }
  if (DEBUG_LEVEL>1) SDFA_instance.showArrays(mask, size,size, "diff_(createNonlinMask)");
  mask= lowPassGauss(mask,   2*maskSigma, true);
  if (DEBUG_LEVEL>1) SDFA_instance.showArrays(mask, size,size, "diff-smooth_(createNonlinMask)");
  double k=1.0/(maskHi-maskLo);
   if (maskHi>maskLo) {
    for (i=0;i<length;i++) {
       if (mask[i]<maskLo) mask[i]=0.0;
       else if (mask[i]>maskHi) mask[i]=1.0;
       else mask[i]=k*(mask[i]-maskLo);
    }
  }
  if (DEBUG_LEVEL>1) SDFA_instance.showArrays(mask, size,size, "mask_(createNonlinMask)");
  return mask;
}

/* Combine 2 images and a mask - when mask is greater or equal 1.0 - use frist image pixels,
    when mask is less or equal 0 - use second image pixels, between - linear interpolate */

public double [][] combineMaskImages (double [][] firstPixels,
                                       double [][] secondPixels,
                                         double [] mask) {
  if ((firstPixels==null) || (secondPixels==null) || (mask==null)) return null;
  int length=firstPixels.length;
  if (secondPixels.length!=length) return null;
  int i;
  double [][] result = new double [length][];
  for (i=0;i<length;i++) {
    result[i]=combineMaskImages (firstPixels[i], secondPixels[i], mask); // OK to pass null
  }
  return result;
}



public double [] combineMaskImages (double [] firstPixels,
                                     double [] secondPixels,
                                     double [] mask) {
  if ((firstPixels==null) || (secondPixels==null) || (mask==null)) return null;
  int length=firstPixels.length;
  if ((secondPixels.length!=length) || (mask.length!=length)) return null;
  double [] result=new double [length];
  int i;
  double m;
  for (i=0;i<length;i++) {
    m=mask[i];
    if (m<0.0) m=0.0;
    else if (m>1.0) m=1.0;
    result[i]=m*firstPixels[i]+(1.0-m)*secondPixels[i];
  }
  return result;
}



/* Create gaussian kernel with specified sigma, same center (lateral chromatic) as in the input kernel */
public double [] extractLateralChromaticFromReverseKernel(double [] smoothReverseKernel, double sigma) {
   return extractLateralChromaticFromReverseKernel(smoothReverseKernel, sigma, 0.5);
}

public double [] extractLateralChromaticFromReverseKernel(double [] smoothReverseKernel,   // square array of reversed PSF
                                                                           double sigma,   // gauss sigma to be used in the result array
                                                                        double fraction) { // fraction of the maximal value to use as a bottom of the part used for centroid calculation
  if (smoothReverseKernel==null) return null;
  int length=smoothReverseKernel.length;
  int size=(int) Math.sqrt(length);
  int i,j;
  double [] centroidXY=  extractCentroidFromReverseKernel(smoothReverseKernel,   // square array of reversed PSF
                                                                     fraction); // fraction of the maximal value to use as a bottom of the part used for centroid calculation


  if (DEBUG_LEVEL>1)  System.out.println("extractLateralChromaticFromReverseKernel(): centroidXY[0]="+IJ.d2s(centroidXY[0],3)+", centroidXY[1]="+IJ.d2s(centroidXY[1],3));


  double nSigma2=(4*sigma)*(4*sigma);
  double x,y,x2,y2;
  double [] gaussian=new double[length];
  double [] gaussianX=new double[size];
  double [] gaussianY=new double[size];
  double minsigma=0.1;
  double k=(sigma<minsigma)?(0.5/minsigma/minsigma):(0.5/sigma/sigma);
  for ( i=0; i<size;i++) {
    x=i-size/2-centroidXY[0];
    x2=x*x;
    if (x2>nSigma2) gaussianX[i]=0.0;
    else gaussianX[i]=Math.exp(-x2*k);
    y=i-size/2-centroidXY[1];
    y2=y*y;
    if (y2>nSigma2) gaussianY[i]=0.0;
    else gaussianY[i]=Math.exp(-y2*k);
  }
/*
  if (DEBUG_LEVEL>1) {
    for ( i=0; i<size;i++) {
        System.out.println("extractLateralChromaticFromReverseKernel(), I="+I+": centroidXY[0]="+IJ.d2s(centroidXY[0],3)+", centroidXY[1]="+IJ.d2s(centroidXY[1],3));
    }
  }
*/
  double sum=0.0;
  double d;
  for ( i=0; i<size;i++) for (j=0;j<size;j++) {
    d=gaussianX[j]*gaussianY[i];
    sum+=d;
    gaussian[i*size+j]=d;
  }
  k=1.0/sum;
  for (i=0;i<length;i++) gaussian[i]*=k;
  return gaussian;
}

/* Create gaussian kernel with specified sigma, same center (lateral chromatic) as in the input kernel */
public double [] extractLateralChromaticFromDirectKernel(double [] smoothDirectKernel, int gaussianSize, double sigma) {
   return extractLateralChromaticFromDirectKernel(smoothDirectKernel, gaussianSize, sigma, 0.5);
}

public double [] extractLateralChromaticFromDirectKernel(double [] smoothDirectKernel,   // square array of reversed PSF
                                                                             int size,
                                                                         double sigma,   // gauss sigma to be used in the result array
                                                                      double fraction) { // fraction of the maximal value to use as a bottom of the part used for centroid calculation
  if (smoothDirectKernel==null) return null;
//  int length=smoothDirectKernel.length;
//  int size=(int) Math.sqrt(length);
  int length=size*size;
  int i,j;
  double [] centroidXY=  extractCentroidFromReverseKernel(smoothDirectKernel,   // square array of direct PSF
                                                                   fraction); // fraction of the maximal value to use as a bottom of the part used for centroid calculation
  centroidXY[0]=-centroidXY[0];
  centroidXY[1]=-centroidXY[1];

  if (DEBUG_LEVEL>1)  System.out.println("extractLateralChromaticFromDirectKernel(): centroidXY[0]="+IJ.d2s(centroidXY[0],3)+", centroidXY[1]="+IJ.d2s(centroidXY[1],3));


  double nSigma2=(4*sigma)*(4*sigma);
  double x,y,x2,y2;
  double [] gaussian=new double[length];
  double [] gaussianX=new double[size];
  double [] gaussianY=new double[size];
  double minsigma=0.1;
  double k=(sigma<minsigma)?(0.5/minsigma/minsigma):(0.5/sigma/sigma);
  for ( i=0; i<size;i++) {
    x=i-size/2-centroidXY[0];
    x2=x*x;
    if (x2>nSigma2) gaussianX[i]=0.0;
    else gaussianX[i]=Math.exp(-x2*k);
    y=i-size/2-centroidXY[1];
    y2=y*y;
    if (y2>nSigma2) gaussianY[i]=0.0;
    else gaussianY[i]=Math.exp(-y2*k);
  }
/*
  if (DEBUG_LEVEL>1) {
    for ( i=0; i<size;i++) {
        System.out.println("extractLateralChromaticFromReverseKernel(), I="+I+": centroidXY[0]="+IJ.d2s(centroidXY[0],3)+", centroidXY[1]="+IJ.d2s(centroidXY[1],3));
    }
  }
*/
  double sum=0.0;
  double d;
  for ( i=0; i<size;i++) for (j=0;j<size;j++) {
    d=gaussianX[j]*gaussianY[i];
    sum+=d;
    gaussian[i*size+j]=d;
  }
  k=1.0/sum;
  for (i=0;i<length;i++) gaussian[i]*=k;
  return gaussian;
}


public double [] extractCentroidFromReverseKernel(double [] smoothReverseKernel,   // square array of reversed PSF
                                                                double fraction) { // fraction of the maximal value to use as a bottom of the part used for centroid calculation
  if (smoothReverseKernel==null) return null;
  int length=smoothReverseKernel.length;
  int size=(int) Math.sqrt(length);
/* find absolute maximum  */
  int i;
  double max=smoothReverseKernel[0];
  int imax=0;
  for (i=0; i<length;i++) if (max<smoothReverseKernel[i]) {
     max=smoothReverseKernel[i];
     imax=i;
  }
  int [][]  clusterMask =    findClusterOnPSF(smoothReverseKernel, // PSF function, square array (use smooth array)
                                                        -fraction, // fraction of energy in the pixels to be used (or minimal level if it is negative)
                                                     (imax % size),  // location of a start point, x-coordinate
                                                     (imax / size),  // location of a start point, y-coordinate
                                                           "");
  double [] centroidXY= calcCentroidFromCenter(smoothReverseKernel, // use original array (mask from the smoothed one)
                                                       clusterMask, // integer mask -0 - don't use this pixel, 1 - use it
                                                          fraction);// subtract level below fraction*max

  return centroidXY;

}










 /* returns array of 3 arrays: first two are 3-element wave vectors (x,y,phase), last - 3-rd order correction coefficients */
  public double[][] findPatternDistorted(double [][] bayer_pixels, // pixel array to process (no widowing!), two greens will be used
                                                   double gamma, // gamma to use for power spectrum for correlation
                                                   double sigma, // high-pass gaussian filter sigma when correlating power spectrum
                                           int diff_spectr_corr, // maximal distance between maximum on spectrum and predicted maximum on autocorrelation of gamma(|spectrum|)
                                         double shrink_clusters, // Shrink clusters by this ratio (remove lowest) after initial separation
                                           int multiples_to_try, // try this number of maximums proportionally farther from 0,0 than the two closest (increase precision)
                                               double deviation, // when looking for maximums - maximal distance from predicted from the lower order one
                                            int deviation_steps, // maximal iterations when looking for local maximum
                                                 boolean greens, // this is a pattern for combined greens (diagonal), adjust results accordingly
                                                 String title) { // title prefix to use for debug  images
      if (bayer_pixels==null) return null;
      if (bayer_pixels.length<4) return null;
      if (bayer_pixels[0]==null) return null;
      if (bayer_pixels[3]==null) return null;
      int size2=bayer_pixels[0].length;
      int size=(int) Math.sqrt(size2);
      int hsize=size/2;
      int hsize2=hsize*hsize;
      double [] quarterHamming=initHamming(hsize);
      double [] patternCorr=new double[6]; // second order non-linear pattern correction (perspective+distortion)
      double [] green0= new double[hsize2];
      double [] green3= new double[hsize2];

      double [][]   quarter_pixels =new double[9][];
      double [][][] quarter_patterns =new double[9][][];
      int [] quarterIndex={0, // top left
                      size/2, // top right
                 size*size/2, // bottom left
             (size+1)*size/2, // bottom right
             (size+1)*size/4, // center
                      size/4, // top
                 size*size/4,// left
             (size+2)*size/4,   // right
           (2*size+1)*size/4};   // bottom

      int i,j,iq;
      int index,qindex;
      if (DEBUG_LEVEL>2) SDFA_instance.showArrays(bayer_pixels, size, size, title+"-bayer");
      for (iq=0; iq<9;iq++) {
//        quarter_pixels[iq]=new double[hsize2];
        index=quarterIndex[iq];
        qindex=0;
        for (i=0;i<hsize;i++) {
          for (j=0;j<hsize;j++) { //quarter_pixels[iq][qindex++]=input_pixels[index++];
            green0[qindex]=  bayer_pixels[0][index];
            green3[qindex++]=bayer_pixels[3][index++];
          }
          quarter_pixels[iq]=combineDiagonalGreens (green0, green3,  hsize, hsize);
          index+=hsize; // jump to the nexty line
        }
        quarter_pixels[iq]= normalizeAndWindow (quarter_pixels[iq], quarterHamming);
        if (DEBUG_LEVEL>2) SDFA_instance.showArrays(quarter_pixels[iq],hsize, hsize, title+"-new"+iq);
        quarter_patterns[iq]   =findPattern(quarter_pixels[iq],
                                                         hsize,
                                                         gamma,
                                                         sigma,  // pattern detection: high-pass filter (0.0 - none) gamma(PS)
                                              diff_spectr_corr,
                                               shrink_clusters,
                                              multiples_to_try,
                                                     deviation,
                                               deviation_steps,
                                                        greens,
                                                 title+"Q_"+iq);
       
        if (quarter_patterns[iq]==null) return null;
      }
      if (DEBUG_LEVEL>2) {
        for (iq=0; iq<9;iq++) {
          System.out.println("Quarter="+iq+
                             " W0x="+     IJ.d2s(quarter_patterns[iq][0][0],4)+
                             " W0y="+     IJ.d2s(quarter_patterns[iq][0][1],4)+
                             " W0_phase="+IJ.d2s(quarter_patterns[iq][0][2],2)+
                             " W1x="+     IJ.d2s(quarter_patterns[iq][1][0],4)+
                             " W1y="+     IJ.d2s(quarter_patterns[iq][1][1],4)+
                             " W1_phase="+IJ.d2s(quarter_patterns[iq][1][2],2));
        }
      }
/* Filter pattern coefficients to make sure they all match between quadrants (match to the center one)*/
     boolean patternsMatchedInitially=matchPatterns(quarter_patterns,quarter_patterns[4]);    // use center pattern
     if (DEBUG_LEVEL>2) {
        System.out.println(patternsMatchedInitially?"All quadrant wave vectors matched initially, no correction needed":"Some quadrant wave vectors were adjusted to match");
     }
     
      patternCorr=calcPatternNonLinear(quarter_patterns); // divide results by ,(FFTSize/2)^2 - only first 5 patterns are used
      if (DEBUG_LEVEL>2) { /* increase LEVEL later */
          System.out.println("Pre- (1000x)   "+
                             " Ax="+     IJ.d2s(1000*patternCorr[0]/(FFTSize/2),5)+
                             " Bx="+     IJ.d2s(1000*patternCorr[1]/(FFTSize/2),5)+
                             " Cx="+     IJ.d2s(1000*patternCorr[2]/(FFTSize/2),5)+
                             " Ay="+     IJ.d2s(1000*patternCorr[3]/(FFTSize/2),5)+
                             " By="+     IJ.d2s(1000*patternCorr[4]/(FFTSize/2),5)+
                             " Cy="+     IJ.d2s(1000*patternCorr[5]/(FFTSize/2),5)+
                             " Dx="+     IJ.d2s(1000*patternCorr[6],5)+
                             " Ex="+     IJ.d2s(1000*patternCorr[7],5)+
                             " Dy="+     IJ.d2s(1000*patternCorr[8],5)+
                             " Ey="+     IJ.d2s(1000*patternCorr[9],5));
      }


      patternCorr=refinePatternNonLinear(quarter_patterns, // [tl,tr,bl,br, center][wv0, wv1][x,y,phase]
                                              patternCorr, //[ax,bx,cx,ay,by,cy]
                                                   hsize ); // distance to quadrats center in sensor pixels ==FFTSize/2



//      for (i=0;i<patternCorr.length;i++)patternCorr[i]/= hsize;
      for (i=0;i<6;i++)patternCorr[i]/= hsize; /* Not linear Dx,Ex, Dy,Ey! */
 
      if (DEBUG_LEVEL>2) { /* increase LEVEL later */
          System.out.println("Corr (1000x)   "+
                             " Ax="+     IJ.d2s(1000*patternCorr[0],5)+
                             " Bx="+     IJ.d2s(1000*patternCorr[1],5)+
                             " Cx="+     IJ.d2s(1000*patternCorr[2],5)+
                             " Ay="+     IJ.d2s(1000*patternCorr[3],5)+
                             " By="+     IJ.d2s(1000*patternCorr[4],5)+
                             " Cy="+     IJ.d2s(1000*patternCorr[5],5)+
                             " Dx="+     IJ.d2s(1000*patternCorr[6],5)+
                             " Ex="+     IJ.d2s(1000*patternCorr[7],5)+
                             " Dy="+     IJ.d2s(1000*patternCorr[8],5)+
                             " Ey="+     IJ.d2s(1000*patternCorr[9],5));
      }

      double [][]result=new double [3][];
      result[0]=quarter_patterns[4][0].clone();
      result[1]=quarter_patterns[4][1].clone();
      result[2]=patternCorr.clone();
      return result;
}


/*************************************************************************************/



/* converting coordinates between sensor pixels and Bayer sub-arrays.
   color 0 - green, top left (normally not used)
   color 1 - red, top right
   color 2 - blue, bottom left
   color 3 - green, bottom right (normally not used)
   color 4 - both greens, diagonal array
   color 5 - both greens, checkerboard array
   subpixel - number of sub-pixels in each direction (vert/hor) of each Bayer sub-array (currently using 4) */
 public double [] bayerToSensor ( double [] xy,
                                  int color,
                                  int subPixel
                                 ) {
    return bayerToSensor (xy[0], xy[1], color, subPixel);
 }

 public double [] bayerToSensor ( double x,
                                  double y,
                                  int color,
                                  int subPixel
                                 ) {
    double [] xy=new double[2];
    switch (color) {
      case 0:xy[0]=2.0*x/subPixel;     xy[1]=2.0*y/subPixel;     break;
      case 1:xy[0]=2.0*x/subPixel+1.0; xy[1]=2.0*y/subPixel;     break;
      case 2:xy[0]=2.0*x/subPixel;     xy[1]=2.0*y/subPixel+1.0; break;
      case 3:xy[0]=2.0*x/subPixel+1.0; xy[1]=2.0*y/subPixel+1.0; break;
      case 4:xy[0]=(x+y)/subPixel;     xy[1]= (y-x)/subPixel;    break;
    }
    return xy;
}
 public double [] sensorToBayer ( double [] xy,
                                  int color,
                                  int subPixel) {
    return sensorToBayer (xy[0], xy[1], color, subPixel);
 }


 public double [] sensorToBayer ( double x,
                                  double y,
                                  int color,
                                  int subPixel) {
    double [] xy=new double[2];
    switch (color) {
      case 0:xy[0]=0.5*x*    subPixel; xy[1]=0.5*y*    subPixel; break;
      case 1:xy[0]=0.5*(x-1)*subPixel; xy[1]=0.5*y*    subPixel; break;
      case 2:xy[0]=0.5*x    *subPixel; xy[1]=0.5*(y-1)*subPixel; break;
      case 3:xy[0]=0.5*(x-1)*subPixel; xy[1]=0.5*(y-1)*subPixel; break;
      case 4:xy[0]=0.5*(x-y)*subPixel; xy[1]=0.5*(x+y)*subPixel; break;
    }
    return xy;
}

/* shift (like lateral chromatic aberration) in Bayer component to sensor pixels */

 public double [] shiftBayerToSensor ( double [] dxy,
                                           int color,
                                        int subPixel) {
    return shiftBayerToSensor (dxy[0], dxy[1], color, subPixel);
 }

 public double [] shiftBayerToSensor ( double dx,
                                       double dy,
                                       int color,
                                      int subPixel) {
    double [] dxy=new double[2];
    switch (color) {
      case 5:
      case 0:
      case 1:
      case 2:
      case 3:dxy[0]=2.0*dx/subPixel;  dxy[1]= 2.0*dy/subPixel;  break;
      case 4:dxy[0]=(dx+dy)/subPixel; dxy[1]= (dy-dx)/subPixel; break;
    }
    if (DEBUG_LEVEL>2)  System.out.println("shiftBayerToSensor(), color="+color+" subPixel="+subPixel+" ["+IJ.d2s(dx,3)+"/"+IJ.d2s(dy,3)+"] ->["+IJ.d2s(dxy[0],3)+"/"+IJ.d2s(dxy[1],3)+"]");
    return dxy;
}

 public double [] shiftSensorToBayer ( double [] dxy,
                                  int color,
                                  int subPixel) {
    return shiftSensorToBayer (dxy[0], dxy[1], color, subPixel);
 }
 public double [] shiftSensorToBayer ( double dx,
                                       double dy,
                                       int color,
                                     int subPixel) {
    double [] dxy=new double[2];
    switch (color) {
      case 5:
      case 0:
      case 1:
      case 2:
      case 3:dxy[0]=0.5*dx*subPixel;      dxy[1]=0.5*dy*subPixel; break;
      case 4:dxy[0]=0.5*(dx-dy)*subPixel; dxy[1]=0.5*(dx+dy)*subPixel; break;
    }
    if (DEBUG_LEVEL>2)  System.out.println("shiftSensorToBayer(), color="+color+" subPixel="+subPixel+" ["+IJ.d2s(dx,3)+"/"+IJ.d2s(dy,3)+"] ->["+IJ.d2s(dxy[0],3)+"/"+IJ.d2s(dxy[1],3)+"]");

    return dxy;
}




 /* convolves pixels array wi th the kernel, resulting in the same size array as the pixels (no roll-over or extending)
     kernel is considered cenered at floor(size/2) */

 public double [] convolveWithKernel (double [] pixels,
                                       int pixelsWidth,
                                      double [] kernel,
                                       int kernelWidth) {
   int pixelsHeight=pixels.length/pixelsWidth;
   int kernelHeight=kernel.length/kernelWidth;
   double [] result = new double [pixelsHeight*pixelsWidth];
   int i,xp,yp,xk,yk,xk0,xk1,yk0,yk1,xr,yr;
   int ykc=kernelHeight/2;
   int xkc=kernelWidth/2;
   double a;
   for (i=0;i<result.length;i++) result[i]=0.0;
   for (yp=0;yp<pixelsHeight;yp++) {
     yk0=0; if (yk0<(ykc-yp)) yk0=(ykc-yp);
     yk1=kernelHeight; if (yk1>(pixelsHeight-yp+ykc)) yk1=(pixelsHeight-yp+ykc);
     for (xp=0;xp<pixelsWidth;xp++) {
       xk0=0; if (xk0<(xkc-xp)) xk0=(xkc-xp);
       xk1=kernelWidth;  if (xk1>(pixelsWidth -xp+xkc)) xk1=(pixelsWidth -xp+xkc);
       for (yk=yk0;yk<yk1;yk++) {
         yr=yp+yk-ykc;
         for (xk=xk0;xk<xk1; xk++) {
           xr=xp+xk-xkc;
//           result[yr*pixelsWidth+xr]+=pixels[yp*pixelsWidth+xp]*kernel[yk*kernelWidth+xk];
           a=pixels[yp*pixelsWidth+xp];
           a*=kernel[yk*kernelWidth+xk];
   if (((yr*pixelsWidth+xr)<0) || ((yr*pixelsWidth+xr)>=result.length)) {
        System.out.println(" yp="+yp+" xp="+xp+" ykc="+ykc+" xkc="+xkc+" yk0="+yk0+" yk1="+yk1+" xk0="+xk0+" xk1="+xk1);
        System.out.println(" yk="+yk+" xk="+xk+" yr="+yr+" xr="+xr);

   }
           result[yr*pixelsWidth+xr]+=a;
         }
       }
     }
   }
   return result;
 }

/*

yr=yp+yk-ykc
yrmin=yp+yk0-ykc
yp+yk0-ykc>=0
yk0>=ykc-yp
yr>=0
yr<height
yp+yk1-ykc<=height

yk1<=height-yp+ykc


*/



/* Makes sure all pairs of wavevectors match, swap and/or negate if needed.
    Returns true if no correction was required, false - if some vector(s) were adjusted.*/
  boolean matchPatterns(double [][][] patterns) {
    return matchPatterns(patterns, patterns[0].clone());
  }

  boolean matchPatterns(double [][][] patterns, double [][] targetPattern) {
    int n,i,j;
    double [][] sp=new double[2][2];
    double [] target_lengths= new double[2];
    double [] current_lengths=new double[2];
    double [] swap_wv;
    double    swap;
    boolean noCorrectionWasNeeded=true;
    if (patterns==null) return false;
    if (targetPattern==null) return false;
    for (n=0;n<patterns.length; n++) {
      for (i=0;i<2;i++) {
        target_lengths[i]= Math.sqrt(targetPattern[i][0]*targetPattern[i][0]+targetPattern[i][1]*targetPattern[i][1]);
        current_lengths[i]=Math.sqrt(patterns[n][i][0]*patterns[n][i][0]+patterns[n][i][1]*patterns[n][i][1]);

      }
      for (i=0;i<2;i++) for (j=0;j<2;j++){
        sp[i][j]=    (patterns[n][i][0]*targetPattern[j][0]+patterns[n][i][1]*targetPattern[j][1])/current_lengths[i]/ target_lengths[j];
      }
/* Swap vectors regardless of sign (parallel/anti-parallel)*/
      if ((Math.abs(sp[0][0])<Math.abs(sp[0][1])) || (Math.abs(sp[0][0])<Math.abs(sp[1][0]))) {
        noCorrectionWasNeeded=false;
        if (DEBUG_LEVEL>0) System.out.println("Swapped wave vectors in quadrant "+n);
        swap_wv=patterns[n][0];   patterns[n][0]=patterns[n][1]; patterns[n][1]=swap_wv;
        swap=sp[0][0];            sp[0][0]=sp[1][0];             sp[1][0]=swap;
        swap=sp[0][1];            sp[0][1]=sp[1][1];             sp[1][1]=swap;
      }
/* Now correct vector signs if needed */
      for (i=0;i<2;i++) {
        if (sp[i][i] <0) {
          noCorrectionWasNeeded=false;
          if (DEBUG_LEVEL>0) System.out.println("Changing wave vector "+(i+1)+" direction in quadrant "+n);
          for (j=0;j<patterns[n][i].length;j++) patterns[n][i][j]=-patterns[n][i][j]; /// Will negate phase if available
        }
      }
    }
    return noCorrectionWasNeeded;
  }

/* converts 2 wave vectors (WVx,WVy,phase) into two checker pattern vectors (VPx,VPy, phase)
    phase is in the point x=0,y=0*/
  double [][] waveVectorsToPatternVectors(double [] wv0, double [] wv1) {
    double [][] v=new double [2][3];
    double vect_wv0_x_wv1=wv0[0]*wv1[1]-wv0[1]*wv1[0];
//    v[0][0]= wv0[1]/vect_wv0_x_wv1;
//    v[0][1]=-wv0[0]/vect_wv0_x_wv1;
//    v[1][0]=-wv1[1]/vect_wv0_x_wv1;
//    v[1][1]= wv1[0]/vect_wv0_x_wv1;

    v[0][0]= wv1[1]/vect_wv0_x_wv1;
    v[0][1]=-wv1[0]/vect_wv0_x_wv1;
    v[1][0]=-wv0[1]/vect_wv0_x_wv1;
    v[1][1]= wv0[0]/vect_wv0_x_wv1;

    v[0][2]=wv0[2];
    v[1][2]=wv1[2];
/* "white square" center had coordinates
   -wv0[0]
*/


   return v;
  }


/**
  Refining non-linear mesh matching by comparing phases in the centers of 4 quadrants and the very center.
  Can only compensate to a fraction of mesh period (TBD - total range, probably +/-half period),
  so non-linear coefficients should be already known to that precision
  9 measurements are used here - top-left, top-right,bottom-left, bottom-right, center, top, left, right,bottom

*/

  double [] refinePatternNonLinear(double [][][] qp, // [tl,tr,bl,br, center][wv0, wv1][x,y,phase]
                                   double [] nonlin, //[ax,bx,cx,ay,by,cy]
                                   int size ) { // distance to quadrats center in sensor pixels ==FFTSize/2

    int iq,i,j;
    double [][] xy=new double [qp.length][2];
    double [] uv= new double [2];
    double [] duv=new double [2];
    
//    double [][][] wl=new double [5][2][2]; // Wave length vectors - same direction as wavevectors, length=distance between wavefronts
    double [][][] wp=new double [9][2][3]; // pattern vectors (with phase)
    double x1,y1;
    double xq=0;
    double yq=0;
// probably only wl[4] is needed
//    for (iq=0;iq<5;iq++) for (i=0;i<2;i++) for (j=0;j<2;j++) wl[iq][i][j]=qp[iq][i][j]/(qp[iq][i][0]*qp[iq][i][0]+qp[iq][i][1]*qp[iq][i][1]);
    for (iq=0;iq<9;iq++) for (i=0;i<2;i++) for (j=0;j<2;j++) wp[iq]= waveVectorsToPatternVectors(qp[iq][0], qp[iq][1]);

      if (DEBUG_LEVEL>2) { /* increase LEVEL later */
         for (iq=0;iq<5;iq++){
          System.out.println(" wp["+iq+"][0][0]="+     IJ.d2s(wp[iq][0][0],4)+
                             " wp["+iq+"][0][1]="+     IJ.d2s(wp[iq][0][1],4)+
                             " wp["+iq+"][0][2]="+     IJ.d2s(wp[iq][0][2],4)+"("+ IJ.d2s(wp[iq][0][2]/2/Math.PI,4)+")"+
                             " wp["+iq+"][1][0]="+     IJ.d2s(wp[iq][1][0],4)+
                             " wp["+iq+"][1][1]="+     IJ.d2s(wp[iq][1][1],4)+
                             " wp["+iq+"][1][2]="+     IJ.d2s(wp[iq][1][2],4)+"("+ IJ.d2s(wp[iq][1][2]/2/Math.PI,4)+")");
         }
      }
    for (iq=0;iq<9;iq++) {
      if (iq==4) continue; // nothing to calculate in the center
      switch (iq){
        case 0: xq=-1.0;yq=-1.0; break;
        case 1: xq= 1.0;yq=-1.0; break;
        case 2: xq=-1.0;yq= 1.0; break;
        case 3: xq= 1.0;yq= 1.0; break;
        case 4: xq= 0.0;yq= 0.0; break;
        case 5: xq= 0.0;yq=-1.0; break;
        case 6: xq=-1.0;yq= 0.0; break;
        case 7: xq= 1.0;yq= 0.0; break;
        case 8: xq= 0.0;yq= 1.0; break;
      }
//      x1=size*(xq + nonlin[0]*xq*xq+ nonlin[1]*yq*yq+ nonlin[2]*xq*yq +nonlin[6]*xq+ nonlin[7]*yq); // in pixels
//      y1=size*(yq + nonlin[3]*xq*xq+ nonlin[4]*yq*yq+ nonlin[5]*xq*yq +nonlin[8]*xq+ nonlin[9]*yq); // in pixels

      x1=size*(xq + nonlin[0]*xq*xq+ nonlin[1]*yq*yq+ 2* nonlin[2]*xq*yq +nonlin[6]*xq+ nonlin[7]*yq); // in pixels
      y1=size*(yq + nonlin[3]*xq*xq+ nonlin[4]*yq*yq+ 2* nonlin[5]*xq*yq +nonlin[8]*xq+ nonlin[9]*yq); // in pixels


//      for (i=0;i<2;i++) for (j=0;j<2;j++) wl[iq][i][j]=qp[iq][i][j]/(qp[iq][i][0]*qp[iq][i][0]+qp[iq][i][1]*qp[iq][i][1]);
/* convert x1,y1 into wp vector coordiantes */
      uv[0]=(wp[4][1][1]*x1-wp[4][1][0]*y1)/(wp[4][0][0]*wp[4][1][1]-wp[4][0][1]*wp[4][1][0]); // wl in center vectors, not local !
      uv[1]=(wp[4][0][0]*y1-wp[4][0][1]*x1)/(wp[4][0][0]*wp[4][1][1]-wp[4][0][1]*wp[4][1][0]);
//      duv[0]=uv[0]-Math.round(uv[0])- (qp[iq][0][2]-qp[4][0][2])/(Math.PI*2);
//      duv[1]=uv[1]-Math.round(uv[1])- (qp[iq][1][2]-qp[4][1][2])/(Math.PI*2);
/* Actually phases seem to be the same? */
      duv[0]=uv[0]-Math.round(uv[0])- (wp[iq][0][2]-wp[4][0][2])/(Math.PI*2);
      duv[1]=uv[1]-Math.round(uv[1])- (wp[iq][1][2]-wp[4][1][2])/(Math.PI*2);


      if (DEBUG_LEVEL>2) { /* increase LEVEL later */
          System.out.println("iq=  "+ iq+
                             " x1="+     IJ.d2s(x1,4)+
                             " y1="+     IJ.d2s(y1,4)+
                             " uv[0]"+   IJ.d2s(uv[0],4)+
                             " uv[1]"+   IJ.d2s(uv[1],4)+
                             " duv[0]"+   IJ.d2s(duv[0],4)+
                             " duv[1]"+   IJ.d2s(duv[1],4));
      }
/* re-normalize phases to be +/- 0.5 (+/-PI) range */
      duv[0]=duv[0]-Math.round(duv[0]);
      duv[1]=duv[1]-Math.round(duv[1]);
      if (DEBUG_LEVEL>2) { /* increase LEVEL later */
          System.out.println("iq=  "+ iq+" ------- "+
                             " duv[0]"+   IJ.d2s(duv[0],4)+
                             " duv[1]"+   IJ.d2s(duv[1],4));
      }
/* Fix half period vertical/half period horizontal shift - is that needed?*/
      if (Math.abs(duv[0])>0.25) {
        duv[0]+=0.5;
        duv[1]+=0.5;
        duv[0]=duv[0]-Math.round(duv[0]);
        duv[1]=duv[1]-Math.round(duv[1]);
        if (DEBUG_LEVEL>2) {
          System.out.println("Correct phase shift >0.25 in quadrant "+ iq+ ", now"+
                             " duv[0]"+   IJ.d2s(duv[0],4)+
                             " duv[1]"+   IJ.d2s(duv[1],4));
        }
      }

/* Verify here that phase adjustment is within range, fail otherwise */
      if ((Math.abs(duv[0])>0.5) || (Math.abs(duv[1])>0.5)) {
        if (DEBUG_LEVEL>0) {
          System.out.println("Error: in quadrant "+ iq+" - attempted to adjust phase too much (>+/- pi/2), keeping initial parematers");
        }
        return nonlin;
      }
/* convert duv to x,y */
//      xy[iq][0]=x1+wl[iq][0][0]*duv[0]+wl[iq][1][0]*duv[1];
//      xy[iq][1]=y1+wl[iq][0][1]*duv[0]+wl[iq][1][1]*duv[1];
//      xy[iq][0]=x1+wl[4][0][0]*duv[0]+wl[4][1][0]*duv[1];
//      xy[iq][1]=y1+wl[4][0][1]*duv[0]+wl[4][1][1]*duv[1];

/*      xy[iq][0]=x1+wp[4][0][0]*duv[0]+wp[4][1][0]*duv[1];
      xy[iq][1]=y1+wp[4][0][1]*duv[0]+wp[4][1][1]*duv[1]; */

      xy[iq][0]=x1-wp[4][0][0]*duv[0]-wp[4][1][0]*duv[1];
      xy[iq][1]=y1-wp[4][0][1]*duv[0]-wp[4][1][1]*duv[1];


      if (DEBUG_LEVEL>2) { /* increase LEVEL later */
          System.out.println(" xy["+iq+"][0]="+   IJ.d2s(xy[iq][0],4)+
                             " xy["+iq+"][1]="+   IJ.d2s(xy[iq][1],4));
      }
/* convert xy to non-linear differences, remove pixels dimensions - quadrats +/- 1.0 */
      xy[iq][0]=(xy[iq][0]-size*xq)/size;
      xy[iq][1]=(xy[iq][1]-size*yq)/size;

      if (DEBUG_LEVEL>2) { /* increase LEVEL later */
          System.out.println("Quadrant "+iq+": Original non-linear difference was"+
                             " x="+   IJ.d2s(nonlin[0]*xq*xq+ nonlin[1]*yq*yq+ nonlin[2]*xq*yq,4)+
                             " y="+   IJ.d2s(nonlin[3]*xq*xq+ nonlin[4]*yq*yq+ nonlin[5]*xq*yq,4));
          System.out.println("Quadrant "+iq+": Refined  non-linear difference is"+
                             " x="+   IJ.d2s(xy[iq][0],4)+
                             " y="+   IJ.d2s(xy[iq][1],4));
      }

    }

/* Do the refinement itself - recalculate coefficients minimizing errors in 5 points */
/**

      x1=size*(xq + nonlin[0]*xq*xq+ nonlin[1]*yq*yq+ 2* nonlin[2]*xq*yq +nonlin[6]*xq+ nonlin[7]*yq); // in pixels
      y1=size*(yq + nonlin[3]*xq*xq+ nonlin[4]*yq*yq+ 2* nonlin[5]*xq*yq +nonlin[8]*xq+ nonlin[9]*yq); // in pixels
    dx=AX*x^2 +BX*y^2 +2*CX*x*y + DX*x +EX *y
    dy=AY*x^2 +BY*y^2 +2*CY*x*y + DY*x +EY *y
      
    dx0= AX +BX +2*CX -DX -EX
    dy0= AY +BY +2*CY -DY -EY

    dx1= AX +BX -2*CX +DX -EX
    dy1= AY +BY -2*CY +DY -EY

    dx2= AX +BX -2*CX -DX +EX
    dy2= AY +BY -2*CY -DY +EY

    dx3= AX +BX +2*CX +DX +EX
    dy3= AY +BY +2*CY +DY +EY

    dx5=    +BX           -EX
    dy5=    +BY           -EY

    dx6= AX           -DX
    dy6= AY           -DY

    dx7= AX           +DX
    dy7= AY           +DY

    dx8=    +BX           +EX
    dy8=    +BY           +EY

//---------
-(1)   dx0= AX +BX +2*CX -DX -EX
-(2)   dx1= AX +BX -2*CX +DX -EX
+(3)   dx2= AX +BX -2*CX -DX +EX
+(4)   dx3= AX +BX +2*CX +DX +EX
-(5)   dx5=    +BX           -EX
 (6)   dx6= AX           -DX
 (7)   dx7= AX           +DX
+(8)   dx8=    +BX           +EX
-(1)+(2)-(3)+(4)-(6)+(7) DX=()-dx0+dx1-dx2+dx3-dx6+dx7)/6
-(1)-(2)+(3)+(4)-(5)+(8) EX=()-dx0-dx1+dx2+dx3-dx5+dx8)/6
-(1)-(2)+(3)+(4)-(5)+(8) EY=()-dy0-dy1+dy2+dy3-dy5+dy8)/6
-(1)+(2)-(3)+(4)-(6)+(7) DY=()-dy0+dy1-dy2+dy3-dy6+dy7)/6

    AX+BX+2*CX = p1x= (dx0+dx3)/2
    AY+BY+2*CY = p1y= (dy0+dy3)/2
    AX+BX-2*CX = p2x= (dx1+dx2)/2
    AY+BY-2*CY = p2y= (dy1+dy2)/2
       BX      = p3x= (dx5+dx8)/2
       BY      = p3y= (dy5+dy8)/2
    AX         = p4x= (dx6+dx7)/2
    AY         = p4y= (dy6+dy7)/2
// minimizing sum of squares of errors:
    CX= (p1x-p2x)/4
    CY= (p1y-p2y)/4
    AX= (3*p4x-2*p3x+p1x+p2x)/5
    AY= (3*p4y-2*p3y+p1y+p2y)/5
    BX= (3*p3x-2*p4x+p1x+p2x)/5
    BY= (3*p3y-2*p4y+p1y+p2y)/5
*/
    double p1x =(xy[0][0]+xy[3][0])/2;
    double p1y =(xy[0][1]+xy[3][1])/2;

    double p2x =(xy[1][0]+xy[2][0])/2;
    double p2y =(xy[1][1]+xy[2][1])/2;

    double p3x =(xy[5][0]+xy[8][0])/2;
    double p3y =(xy[5][1]+xy[8][1])/2;

    double p4x =(xy[6][0]+xy[7][0])/2;
    double p4y =(xy[6][1]+xy[7][1])/2;
   
    double [] rslt=new double[10];
    rslt[2]=(p1x-p2x)/4; //CX
    rslt[5]=(p1y-p2y)/4; //CY
    rslt[0]= (3*p4x-2*p3x+p1x+p2x)/5 ; //AX
    rslt[3]= (3*p4y-2*p3y+p1y+p2y)/5; // AY
    rslt[1]= (3*p3x-2*p4x+p1x+p2x)/5; // BX
    rslt[4]= (3*p3y-2*p4y+p1y+p2y)/5; // BY
/*
-(1)+(2)-(3)+(4)-(6)+(7) DX=(-dx0+dx1-dx2+dx3-dx6+dx7)/6
-(1)-(2)+(3)+(4)-(5)+(8) EX=(-dx0-dx1+dx2+dx3-dx5+dx8)/6
-(1)-(2)+(3)+(4)-(5)+(8) EY=(-dy0-dy1+dy2+dy3-dy5+dy8)/6
-(1)+(2)-(3)+(4)-(6)+(7) DY=(-dy0+dy1-dy2+dy3-dy6+dy7)/6
*/
    rslt[6]=(-xy[0][0]+xy[1][0]-xy[2][0]+xy[3][0]-xy[6][0]+xy[7][0])/6; // DX=(-dx0+dx1-dx2+dx3-dx6+dx7)/6
    rslt[7]=(-xy[0][0]-xy[1][0]+xy[2][0]+xy[3][0]-xy[5][0]+xy[8][0])/6; // EX=(-dx0-dx1+dx2+dx3-dx5+dx8)/6
    rslt[8]=(-xy[0][1]+xy[1][1]-xy[2][1]+xy[3][1]-xy[6][1]+xy[7][1])/6; // EY=(-dy0-dy1+dy2+dy3-dy5+dy8)/6
    rslt[9]=(-xy[0][1]-xy[1][1]+xy[2][1]+xy[3][1]-xy[5][1]+xy[8][1])/6; // DY=(-dy0+dy1-dy2+dy3-dy6+dy7)/6

    return rslt; // temporary
  }


/**
    Matching non-linear mesh using wave vectors in the centers of four quadrants, phases are ignored here
    processes wave vectors for the four quadrants (top-left,tr,bl,br) and calculates second degree polynominal correction
    X1= Ax*X^2 + Bx*Y^2 + 2Cx*X*Y
    Y1= Ay*X^2 + By*Y^2 + 2Cy*X*Y
    where X and Y are normalized so top left corner is (-1,-1), top right - (+1,-1) , bottom left - (-1,+1), bottom right - (+1,+1)
    returns array of 6 elements {Ax,Bx,Cx,Ay,By,Cy}
*/
  double [] calcPatternNonLinear(double [][][] qp) {
    int iq;
/*
    double [][][] qp = new double [4][2][2];
   for (iq=0;iq<4;iq++) for (i=0;i<2;i++)for (j=0;j<2;j++){
      qp[iq][i][j]=qp_in[iq][i][j]/(qp_in[iq][i][0]*qp_in[iq][i][0]+qp_in[iq][i][1]*qp_in[iq][i][1]);
    }
*/
/* Calculate center WV */
/*    double [][] centerWV={{0.25*(qp[0][0][0]+qp[1][0][0]+qp[2][0][0]+qp[3][0][0]),
                           0.25*(qp[0][0][1]+qp[1][0][1]+qp[2][0][1]+qp[3][0][1])},
                          {0.25*(qp[0][1][0]+qp[1][1][0]+qp[2][1][0]+qp[3][1][0]),
                           0.25*(qp[0][1][1]+qp[1][1][1]+qp[2][1][1]+qp[3][1][1])}};*/


    double [][][] mCorners=    new double [4][][] ; // only two 2x2 are needed, other two - just for debugging
//    double [][][] mCornAverage=new double [2][2][2] ; // average of matrix pairs (0/3, 1/2)
/*
    mCorners[0]= matrixToConvertTwoPairs(qp[0][0], qp[0][1], centerWV[0],centerWV[1]);
    mCorners[1]= matrixToConvertTwoPairs(qp[1][0], qp[1][1], centerWV[0],centerWV[1]);
    mCorners[2]= matrixToConvertTwoPairs(centerWV[0],centerWV[1],qp[2][0], qp[2][1]);
    mCorners[3]= matrixToConvertTwoPairs(centerWV[0],centerWV[1],qp[3][0], qp[3][1]);
    for (iq=0;iq<2;iq++)for (i=0;i<2;i++) for (j=0;j<2;j++) {
     mCornAverage[iq][i][j]= 0.5* (mCorners[iq][i][j]+mCorners[3-iq][i][j]);
    }
*/
    mCorners[0]= matrixToConvertTwoPairs(qp[0][0], qp[0][1], qp[4][0], qp[4][1]);
    mCorners[1]= matrixToConvertTwoPairs(qp[1][0], qp[1][1], qp[4][0], qp[4][1]);
    mCorners[2]= matrixToConvertTwoPairs(qp[2][0], qp[2][1], qp[4][0], qp[4][1]);
    mCorners[3]= matrixToConvertTwoPairs(qp[3][0], qp[3][1], qp[4][0], qp[4][1]);

/**
x,y - image coordinates (distorted), with no center shift (later use effective radius)
x1,y1 - linear mesh

x1=x+Ax*x^2+Bx*y^2+2*Cx*x*y+Dx*x+Ex*y
y1=y+Ay*x^2+By*y^2+2*Cy*x*y+Dy*x+Ey*y

dx1/dx=1+2*Ax*x+2*Cx*y+Dx
dx1/dy=0+2*Cx*x+2*Bx*y+Ex

dy1/dx=0+2*Ay*x+2*Cy*y+Dy
dy1/dy=1+2*Cy*x+2*By*y+Ey

| dx1 | = |  1+2*Ax*x+2*Cx*y+Dx     2*Cx*x+2*Bx*y+Ex | * | dx |
| dy1 |   |    2*Ay*x+2*Cy*y+Dy   1+2*Cy*x+2*By*y+Ey |   | dy |

for x=-1,y=-1
M[0]=|  1-2*Ax-2*Cx+Dx    -2*Cx-2*Bx+Ex |
     |   -2*Ay-2*Cy+Dy   1-2*Cy-2*By+Ey |

for x=+1,y=-1
M[1]=|  1+2*Ax-2*Cx+Dx    +2*Cx-2*Bx+Ex |
     |   +2*Ay-2*Cy+Dy   1+2*Cy-2*By+Ey |

for x=-1,y=+1
M[2]=|  1-2*Ax+2*Cx+Dx    -2*Cx+2*Bx+Ex |
     |   -2*Ay+2*Cy+Dy   1-2*Cy+2*By+Ey |

for x=+1,y=+1
M[2]=|  1+2*Ax+2*Cx+Dx    +2*Cx+2*Bx+Ex |
     |   +2*Ay+2*Cy+Dy   1+2*Cy+2*By+Ey |


*/
   double [] det=  {mCorners[0][0][0]*mCorners[0][1][1]-mCorners[0][0][1]*mCorners[0][1][0],
                    mCorners[1][0][0]*mCorners[1][1][1]-mCorners[1][0][1]*mCorners[1][1][0],
                    mCorners[2][0][0]*mCorners[2][1][1]-mCorners[2][0][1]*mCorners[2][1][0],
                    mCorners[3][0][0]*mCorners[3][1][1]-mCorners[3][0][1]*mCorners[3][1][0]};

   double [][][] M={{{ mCorners[0][1][1]/det[0], -mCorners[0][1][0]/det[0]},
                     {-mCorners[0][0][1]/det[0],  mCorners[0][0][0]/det[0]}},
                    {{ mCorners[1][1][1]/det[1], -mCorners[1][1][0]/det[1]},
                     {-mCorners[1][0][1]/det[1],  mCorners[1][0][0]/det[1]}},
                    {{ mCorners[2][1][1]/det[2], -mCorners[2][1][0]/det[2]},
                     {-mCorners[2][0][1]/det[2],  mCorners[2][0][0]/det[2]}},
                    {{ mCorners[3][1][1]/det[3], -mCorners[3][1][0]/det[3]},
                     {-mCorners[3][0][1]/det[3],  mCorners[3][0][0]/det[3]}}};
/**
   Overdefined - 16 equations for 10 unknowns
 1   M[0][0][0]=1-2*Ax-2*Cx+Dx
 2   M[0][0][1]= -2*Cx-2*Bx+Ex
 3   M[0][1][0]= -2*Ay-2*Cy+Dy
 4   M[0][1][1]=1-2*Cy-2*By+Ey
 5   M[1][0][0]=1+2*Ax-2*Cx+Dx
 6   M[1][0][1]=  2*Cx-2*Bx+Ex
 7   M[1][1][0]=  2*Ay-2*Cy+Dy
 8   M[1][1][1]=1+2*Cy-2*By+Ey

 9   M[2][0][0]=1-2*Ax+2*Cx+Dx
10   M[2][0][1]= -2*Cx+2*Bx+Ex
11   M[2][1][0]= -2*Ay+2*Cy+Dy
12   M[2][1][1]=1-2*Cy+2*By+Ey
13   M[3][0][0]=1+2*Ax+2*Cx+Dx
14   M[3][0][1]=  2*Cx+2*Bx+Ex
15   M[3][1][0]=  2*Ay+2*Cy+Dy
16   M[3][1][1]=1+2*Cy+2*By+Ey

( 1)   M[0][0][0]=1-2*Ax     -2*Cx                +Dx
( 2)   M[0][0][1]=      -2*Bx-2*Cx                      +Ex
( 3)   M[0][1][0]=                -2*Ay     -2*Cy    +Dy
( 4)   M[0][1][1]=1                    -2*By-2*Cy          +Ey
( 5)   M[1][0][0]=1+2*Ax     -2*Cx                +Dx
( 6)   M[1][0][1]=      -2*Bx+2*Cx                      +Ex
( 7)   M[1][1][0]=                +2*Ay     -2*Cy    +Dy
( 8)   M[1][1][1]=1                    -2*By+2*Cy          +Ey

( 9)   M[0][0][0]=1-2*Ax     +2*Cx                +Dx
(10)   M[0][0][1]=      +2*Bx-2*Cx                      +Ex
(11)   M[0][1][0]=                -2*Ay     +2*Cy    +Dy
(12)   M[0][1][1]=1                    +2*By-2*Cy          +Ey
(13)   M[1][0][0]=1+2*Ax     +2*Cx                +Dx
(14)   M[1][0][1]=      +2*Bx+2*Cx                      +Ex
(15)   M[1][1][0]=                +2*Ay     +2*Cy    +Dy
(16)   M[1][1][1]=1                    +2*By+2*Cy          +Ey

 ( 1)+( 5)+( 9)+(13)                     Dx=( M[0][0][0]+M[1][0][0]+M[2][0][0]+M[3][0][0])/4-1
 ( 2)+( 6)+(10)+(14)                     Ex=( M[0][0][1]+M[1][0][1]+M[2][0][1]+M[3][0][1])/4
 ( 3)+( 7)+(11)+(15)                     Dy=( M[0][1][0]+M[1][1][0]+M[2][1][0]+M[3][1][0])/4
 ( 4)+( 8)+(12)+(16)                     Ey=( M[0][1][1]+M[1][1][1]+M[2][1][1]+M[3][1][1])/4-1

-( 1)+( 5)-( 9)+(13)                     Ax=(-M[0][0][0]+M[1][0][0]-M[2][0][0]+M[3][0][0])/8
-( 2)-( 6)+(10)+(14)                     Bx=(-M[0][0][1]-M[1][0][1]+M[2][0][1]+M[3][0][1])/8
-( 2)+( 6)-(10)+(14)-( 1)-( 5)+( 9)+(13) Cx=(-M[0][0][1]+M[1][0][1]-M[2][0][1]+M[3][0][1]-M[0][0][0]-M[1][0][0]+M[2][0][0]+M[3][0][0])/16
-( 3)+( 7)-(11)+(15)                     Ay=(-M[0][1][0]+M[1][1][0]-M[2][1][0]+M[3][1][0])/8
-( 4)-( 8)+(12)+(16)                     By=(-M[0][1][1]-M[1][1][1]+M[2][1][1]+M[3][1][1])/8
-( 4)+( 8)-(12)+(16)-( 3)-( 7)+(11)+(15) Cy=(-M[0][1][1]+M[1][1][1]-M[2][1][1]+M[3][1][1]-M[0][1][0]-M[1][1][0]+M[2][1][0]+M[3][1][0])/16


*/
    double [] rslt = {(-M[0][0][0]+M[1][0][0]-M[2][0][0]+M[3][0][0])/8, // Ax
                      (-M[0][0][1]-M[1][0][1]+M[2][0][1]+M[3][0][1])/8, // Bx
                      (-M[0][0][1]+M[1][0][1]-M[2][0][1]+M[3][0][1]-M[0][0][0]-M[1][0][0]+M[2][0][0]+M[3][0][0])/16,  // Cx
                      (-M[0][1][0]+M[1][1][0]-M[2][1][0]+M[3][1][0])/8,  // Ay
                      (-M[0][1][1]-M[1][1][1]+M[2][1][1]+M[3][1][1])/8,  // By
                      (-M[0][1][1]+M[1][1][1]-M[2][1][1]+M[3][1][1]-M[0][1][0]-M[1][1][0]+M[2][1][0]+M[3][1][0])/16,   // Cy
                      ( M[0][0][0]+M[1][0][0]+M[2][0][0]+M[3][0][0])/4-1, // Dx
                      ( M[0][0][1]+M[1][0][1]+M[2][0][1]+M[3][0][1])/4,   // Ex
                      ( M[0][1][0]+M[1][1][0]+M[2][1][0]+M[3][1][0])/4,   // Dy
                      ( M[0][1][1]+M[1][1][1]+M[2][1][1]+M[3][1][1])/4-1  // Ey
                       };

   if (DEBUG_LEVEL>2) { /* increase LEVEL later */
/*          System.out.println("Center   "+
                             " W0x="+     IJ.d2s(centerWV[0][0],4)+
                             " W0y="+     IJ.d2s(centerWV[0][1],4)+
                             " W1x="+     IJ.d2s(centerWV[1][0],4)+
                             " W1y="+     IJ.d2s(centerWV[1][1],4));*/
        for (iq=0; iq<4;iq++) {
          System.out.println("Matrix= "+iq+
                             " M00="+     IJ.d2s(mCorners[iq][0][0],4)+
                             " M01="+     IJ.d2s(mCorners[iq][0][1],4)+
                             " M10="+     IJ.d2s(mCorners[iq][1][0],4)+
                             " M11="+     IJ.d2s(mCorners[iq][1][1],4));
        }
/*
        for (iq=0; iq<2;iq++) {
          System.out.println("Average "+iq+
                             " M00="+     IJ.d2s(mCornAverage[iq][0][0],4)+
                             " M01="+     IJ.d2s(mCornAverage[iq][0][1],4)+
                             " M10="+     IJ.d2s(mCornAverage[iq][1][0],4)+
                             " M11="+     IJ.d2s(mCornAverage[iq][1][1],4));
        }
*/
        for (iq=0; iq<4;iq++) {
          System.out.println(" M["+iq+"][0][0]="+  IJ.d2s(M[iq][0][0],4)+
                             " M["+iq+"][0][1]="+  IJ.d2s(M[iq][0][1],4)+
                             " M["+iq+"][1][0]="+  IJ.d2s(M[iq][1][0],4)+
                             " M["+iq+"][1][1]="+  IJ.d2s(M[iq][1][1],4));
        }
    }
   if (DEBUG_LEVEL>2) { /* increase LEVEL later */
          System.out.println("Corr   "+
                             " Ax="+     IJ.d2s(rslt[0],5)+
                             " Bx="+     IJ.d2s(rslt[1],5)+
                             " Cx="+     IJ.d2s(rslt[2],5)+
                             " Ay="+     IJ.d2s(rslt[3],5)+
                             " By="+     IJ.d2s(rslt[4],5)+
                             " Cy="+     IJ.d2s(rslt[5],5)+
                             " Dx="+     IJ.d2s(rslt[6],5)+
                             " Ex="+     IJ.d2s(rslt[7],5)+
                             " Dy="+     IJ.d2s(rslt[8],5)+
                             " Ey="+     IJ.d2s(rslt[9],5));
    }
    return rslt;
  };



/* calculates 2x2 matrix that converts two pairs of vectors: u2=M*u1, v2=M*v1*/
  double [][] matrixToConvertTwoPairs(double [] u1, double [] v1, double [] u2, double [] v2) {
    double [][] rslt= {{(u2[0]*v1[1]-v2[0]*u1[1])/(u1[0]*v1[1]-v1[0]*u1[1]),
                        (v2[0]*u1[0]-u2[0]*v1[0])/(u1[0]*v1[1]-v1[0]*u1[1])},
                       {(u2[1]*v1[1]-v2[1]*u1[1])/(u1[0]*v1[1]-v1[0]*u1[1]),
                        (v2[1]*u1[0]-u2[1]*v1[0])/(u1[0]*v1[1]-v1[0]*u1[1])}};
    return rslt;
  }

  public double [] convert2d_1d(double [][] pixels){
    int i,j;
    int width=pixels[0].length;
    double [] rslt=new double[pixels.length*pixels[0].length];
    for (i=0;i<pixels.length;i++) for (j=0;j<width;j++) rslt[i*width+j]=pixels[i][j];
    return rslt;
  }

  public int [] convert2d_1d(int [][] pixels){
    int i,j;
    int width=pixels[0].length;
    int [] rslt=new int[pixels.length*pixels[0].length];
    for (i=0;i<pixels.length;i++) for (j=0;j<width;j++) rslt[i*width+j]=pixels[i][j];
    return rslt;
  }


/* pixels should be a square array, zero is in the center (/center+0.5 for even dimensions) */
  public double [] calcCentroidFromCenter(double [] pixels) {return calcCentroidFromCenter(pixels, (int[]) null, 0.0);}
  public double [] calcCentroidFromCenter(double [] pixels, // square pixel array
                                              int[][] mask, // integer mask -0 - don't use this pixel, 1 - use it
                                           double refLevel) { // subtract this fraction of maximal level from all pixels
    return calcCentroidFromCenter(pixels, convert2d_1d(mask), refLevel);
  }
  public double [] calcCentroidFromCenter(double [] pixels, // square pixel array
                                                int[] mask, // integer mask -0 - don't use this pixel, 1 - use it
                                           double refLevel) { // subtract this fraction of maximal leve from all pixels
    int size = (int) Math.sqrt ( pixels.length);
    int c= size/2;
    double S0=0.0;
    double SX=0.0;
    double SY=0.0;
    double x,y,p;
    int i,j,indx;
    double maxValue = 0.0;
    if (refLevel>0.0) for (i=0;i<pixels.length;i++) if (((mask==null) || (mask[i]>0)) && (pixels[i] > maxValue)) maxValue=pixels[i];

    double minValue=refLevel*maxValue;

    for (i=0;i<size;i++) {
      y=i-c;
      for (j=0;j<size;j++) {
        indx=i*size+j;
        if ((mask==null) || (mask[indx]>0)) {
          x=j-c;
          p=pixels[indx]-minValue;
          if (p>0.0) { // with mask mis-match there could be negative total mask
            S0+=p;
            SX+=p*x;
            SY+=p*y;
          }
        }
      }
    }
    double [] result={SX/S0,SY/S0};
/*
    if ((DEBUG_LEVEL>1) && (mask!=null) && (minValue>0)) {
      float [] floatPixels=new float[pixels.length];
      indx=0;
      for (i=0;i<size;i++) for (j=0;j<size;j++){
        if (mask[indx]==0) floatPixels[indx]=0;
        else               floatPixels[indx]=(float) (pixels[indx]-minValue);
        indx++;
      }
      ImageProcessor ip=new FloatProcessor(size,size);
      ip.setPixels(floatPixels);
      ip.resetMinAndMax();
      ImagePlus imp=  new ImagePlus("centroid", ip);
      imp.show();
    }
*/
    return result;

  }
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
  public double [] cleanupAndReversePSF (double [] psf_pixels,  // input pixels
                                         double invertRange,    // deconvInvert
                                         double cutoff_energy,  // OTF_cutoff_energy
                                         double ellipse_scale,  // ellipse mask size relative to the cluster
                                        boolean ellipse_gauss,  // use Gauss instead of Hamming for ellipse mask
                                         int    decimate,       // decimate frequency to match Bayer component pixels pitch
                                         boolean fold,          // fold high frequency into low, when false - use Hamming to cut off high frequencies
                                         String title           // just for the plot names
                                        ) {
    int size=(int) Math.sqrt(psf_pixels.length);
    double[][][] fft_complex;
    int i,j,ix,iy;
    double a,k,r,r2,k2;

//    float []floatPixels=new float[psf_pixels.length];
    double [] cpixels=psf_pixels.clone();
/* Swapping quadrants, so the center will be 0,0 */
    fht_instance.swapQuadrants(cpixels);
/* get to frequency domain */
    fht_instance.transform(cpixels);
/* Convert from FHT to complex FFT - avoid that in the future, process FHT directly*/
    fft_complex= FHT2FFTHalf (cpixels,size);
    double [][]fft_energy=new double[(size/2)+1][size];
    for (i=0;i<(size/2+1);i++) for (j=0;j<size;j++) {
      fft_energy[i][j]=fft_complex[i][j][0]*fft_complex[i][j][0]+fft_complex[i][j][1]*fft_complex[i][j][1];
    }

    int  [][] clusterPS = findClusterOnPS(fft_energy, cutoff_energy,title);

    double [] ellipse_coeff = findEllipseOnPS(fft_energy, clusterPS, title);
/* create ellipse window using Hamming */
/* TODO: scale radius */
    double [][] ellipseMask=new double [size/2+1][size];
    k2=1/ellipse_scale/ellipse_scale;
    for (i=0;i<(size/2+1);i++) for (j=0;j<size;j++) {
       iy=(i==size/2)?-i:i;
       ix=(j>=(size/2))?(j-size):j;
       if (iy<0) ix=-ix;
       r2=ellipse_coeff[0]*ix*ix+ellipse_coeff[1]*iy*iy+ellipse_coeff[2]*ix*iy;
       if (ellipse_gauss){
         ellipseMask[i][j]=Math.exp(-k2*r2);
       } else {
         r=Math.sqrt(r2)/ellipse_scale;
         ellipseMask[i][j]=(r>1.0)?0.0:(0.54+0.46*Math.cos(r*Math.PI));
       }
    }

/* optionally display selection */
/*
    if (DEBUG_LEVEL>2) {
      ImageProcessor ip_ellipse0 = new FloatProcessor(size,size/2+1);
      float [] ellipsePixels0 = new float [size*(size/2+1)];
      for (i=0;i<ellipsePixels0.length;i++) {
        iy=i/size;
        ix=i%size;
        ellipsePixels0[i]= (float) ellipseMask[iy][ix];
      }
      ip_ellipse0.setPixels(ellipsePixels0);
      ip_ellipse0.resetMinAndMax();
      ImagePlus imp_ellipse0= new ImagePlus(title+"_EL0-MASK_"+cutoff_energy+"-"+ellipse_scale, ip_ellipse0);
      imp_ellipse0.show();
    }
*/

    if (DEBUG_LEVEL>2) {
      ImageProcessor ip_ellipse = new FloatProcessor(size,size);
      float [] ellipsePixels = new float [size*size];
      for (i=0;i<ellipsePixels.length;i++) {
        iy=i/size-size/2;
        ix=i%size-size/2;
        if (iy<0) {
          ix=-ix;
          iy=-iy;
        }
        ix= (ix+size) % size;
        ellipsePixels[i]= (float) ellipseMask[iy][ix];
      }
      ip_ellipse.setPixels(ellipsePixels);
      ip_ellipse.resetMinAndMax();
      ImagePlus imp_ellipse= new ImagePlus(title+"_EL-MASK_"+cutoff_energy+"-"+ellipse_scale, ip_ellipse);
      imp_ellipse.show();
    }
/*
    if (DEBUG_LEVEL>2) {
      ImageProcessor ip_ellipse1 = new FloatProcessor(size,size);
      float [] ellipsePixels1 = new float [size*size];
      for (i=0;i<ellipsePixels1.length;i++) {
        iy=i/size;
        ix=i%size;
        if (iy>size/2) {
          ix=(size-ix)%size;
          iy= size-iy;
        }
        ellipsePixels1[i]= (float) ellipseMask[iy][ix];
      }
      ip_ellipse1.setPixels(ellipsePixels1);
      ip_ellipse1.resetMinAndMax();
      ImagePlus imp_ellipse1= new ImagePlus(title+"_EL1-MASK_"+cutoff_energy+"-"+ellipse_scale, ip_ellipse1);
      imp_ellipse1.show();
    }
*/

/* inverse fft_complex */
    if (invertRange>0.0) {

/// Invert Z for large values, but make them Z - for small ones. So it will be a mixture of correlation and deconvolution
// here the targets are round, but what will th\be the colrrect way fo assymmetrical ones?

/// First - find maximal value

   
      double fft_max=0;
      for (i=0;i<fft_complex.length; i++) for (j=0;j<fft_complex[0].length;j++) {
        r2=fft_complex[i][j][0]*fft_complex[i][j][0]+fft_complex[i][j][1]*fft_complex[i][j][1];
        if (r2>fft_max) fft_max=r2;
      }
      k=Math.sqrt(fft_max)*invertRange;
      k2=k*k;
      for (i=0;i<fft_complex.length; i++) for (j=0;j<fft_complex[0].length;j++) {
        r=Math.sqrt(fft_complex[i][j][0]*fft_complex[i][j][0]+fft_complex[i][j][1]*fft_complex[i][j][1]);
        a=-Math.atan2(fft_complex[i][j][1],fft_complex[i][j][0]); /// was zero for circular targets)
        r=r/(r*r+k2);
        fft_complex[i][j][0]=r*Math.cos(a);
        fft_complex[i][j][1]=r*Math.sin(a);
      }
/* multiply by ellipse window */
      for (i=0;i<fft_complex.length; i++) for (j=0;j<fft_complex[0].length;j++) {
        fft_complex[i][j][0]*=ellipseMask[i][j];
        fft_complex[i][j][1]*=ellipseMask[i][j];
      }
    } else { // Do just the division (low power frequencies will be masked out by ellipse window)
      for (i=0;i<fft_complex.length; i++) for (j=0;j<fft_complex[0].length;j++) if (ellipseMask[i][j]>=0.0){
        r2=fft_complex[i][j][0]*fft_complex[i][j][0]+fft_complex[i][j][1]*fft_complex[i][j][1];
        fft_complex[i][j][0]*= ellipseMask[i][j]/r2;
        fft_complex[i][j][1]*=-ellipseMask[i][j]/r2;
      } else {
        fft_complex[i][j][0]=0.0;
        fft_complex[i][j][1]=0.0;
      }
    }

    double [] pixels=null;
/* convert back original dimension array if there was no decimation or debug is set (in that case both sizes arrays will be converted) */
    if ((decimate==1) || (DEBUG_LEVEL>1)){ // invert pixel-decimated OTF (as input to this function)
/* Convert fft array back to fht array and 
      set fht pixels with new values */
      pixels=FFTHalf2FHT (fft_complex,size);
/* optionally show the result FHT*/
/* transform to space */
      fht_instance.inverseTransform(pixels);
      fht_instance.swapQuadrants(pixels);
/*   return inverted psf pixels */
   }
   if (decimate>1) {
/* create folded array so the transformed back inverted otf will have the same pixel resolution as original arrays */
     int outSize=size/decimate;
     double [] foldedHamming=new double [outSize];
     for (i=0;i<outSize;i++) foldedHamming[i]=0.54+0.46*Math.cos(2*Math.PI*i/outSize);
     double [][][] fft_complex_folded = new double [outSize/2+1][outSize][2];
     for (i=0;i<fft_complex_folded.length; i++) for (j=0;j<fft_complex_folded[0].length;j++) {
       fft_complex_folded[i][j][0]=0.0;
       fft_complex_folded[i][j][1]=0.0;
     }
     if (fold || (decimate==1)) {
       for (i=0;i<size;i++) for (j=0;j<size;j++) {
         iy=i%(size/decimate);
         ix=j%(size/decimate);
         if (iy<=(size/decimate/2)){
           if (i<=size/2) {
             fft_complex_folded[iy][ix][0]+=fft_complex[i][j][0];
             fft_complex_folded[iy][ix][1]+=fft_complex[i][j][1];
           } else  {
             fft_complex_folded[iy][ix][0]+=fft_complex[size-i][(size-j)%size][0];
             fft_complex_folded[iy][ix][1]-=fft_complex[size-i][(size-j)%size][1];
           }
         }
       }
     } else {
       for (iy=0;iy<(outSize/2+1);iy++) {
         i=iy;
         for (ix=0;ix<outSize;ix++) {
           j=ix+((ix > (outSize/2))?(size-outSize):0);
           fft_complex_folded[iy][ix][0]+= foldedHamming[ix]* foldedHamming[iy]*fft_complex[i][j][0];
           fft_complex_folded[iy][ix][1]+= foldedHamming[ix]* foldedHamming[iy]*fft_complex[i][j][1];
         }
/* not sure about Im() sign */
         fft_complex_folded[iy][outSize/2][0]= 0.5*(fft_complex_folded[iy][outSize/2][0]+foldedHamming[outSize/2]* foldedHamming[iy]*fft_complex[i][size-outSize/2][0]);
         fft_complex_folded[iy][outSize/2][1]= 0.5*(fft_complex_folded[iy][outSize/2][1]-foldedHamming[outSize/2]* foldedHamming[iy]*fft_complex[i][size-outSize/2][1]);
       }
       for (ix=1;ix<outSize;ix++) if (ix!=outSize/2){
         j=ix+((ix > (outSize/2))?(size-outSize):0);
         fft_complex_folded[outSize/2][ix][0]= 0.5*(fft_complex_folded[outSize/2][ix][0]+foldedHamming[outSize/2]* foldedHamming[ix]*fft_complex[outSize/2][size-j][0]);
         fft_complex_folded[outSize/2][ix][1]= 0.5*(fft_complex_folded[outSize/2][ix][1]-foldedHamming[outSize/2]* foldedHamming[ix]*fft_complex[outSize/2][size-j][1]);
       }
       fft_complex_folded[outSize/2][outSize/2][0]= 0.5*foldedHamming[outSize/2]*foldedHamming[outSize/2]*( 
                                                          fft_complex[outSize/2][outSize/2][0]
                                                         +fft_complex[outSize/2][size-outSize/2][0]);
       fft_complex_folded[outSize/2][outSize/2][1]= 0.0;
     }
/* Convert folded fft array back to fht array and  set fht pixels with new values */
     pixels=FFTHalf2FHT (fft_complex_folded,size/decimate);
/* transform to space */
     fht_instance.inverseTransform(pixels);
     fht_instance.swapQuadrants(pixels);
/*   return inverted psf pixels */
   }  
   return pixels;
  }




/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
 /* Seems it was an error - we do not need to scale ellipse here because rPSF kernel is not decimated from the PSF - it just inverted zero-padded one */

 public double [] maskReversePSFKernel (double [] psf_pixels, // direct PSF function, square array, may be proportionally larger than reversed
                                        double []rpsf_pixels, // reversed psf, square array
                                        double cutoff_energy, // fraction of energy in the pixels to be used
                                        double ellipse_scale,
                                   double min_mask_threshold, // zero output element if elliptical Gauss mask is below this threshold
                                                String title)
 {
//    int psf_size=(int)Math.sqrt(psf_pixels.length);
    int rpsf_size=(int)Math.sqrt(rpsf_pixels.length);
    double [] masked_rpsf=new double[rpsf_size*rpsf_size];
//    double decimation= ((double) psf_size)/rpsf_size;
    int  [][]selection=   findClusterOnPSF(psf_pixels, cutoff_energy, title);
    double [] ellipse_coeff=findEllipseOnPSF(psf_pixels,  selection,    title);
    int ix,iy;
    double x,y,r2;
    int indx=0;
    double k2=1/ellipse_scale/ellipse_scale;
    double m;

    for (iy=0;iy<rpsf_size;iy++) {
//      y=(iy-rpsf_size/2)*decimation+ellipse_coeff[1];  // scale to the original psf (and ellipse_coeff), move center opposite to that of direct kernel (psf)
      y=iy-rpsf_size/2+ellipse_coeff[1];  // move center opposite to that of direct kernel (psf)
      for (ix=0;ix<rpsf_size;ix++) {
//        x=(ix-rpsf_size/2)*decimation+ellipse_coeff[0]; // scale to the original psf (and ellipse_coeff), move center opposite to that of direct kernel (psf)
        x=ix -rpsf_size/2 +ellipse_coeff[0]; //  move center opposite to that of direct kernel (psf)
        r2=ellipse_coeff[2]*x*x+ellipse_coeff[3]*y*y+ellipse_coeff[4]*x*y;
        m=Math.exp(-k2*r2);
        masked_rpsf[indx]=(m>=min_mask_threshold)?(rpsf_pixels[indx]*Math.exp(-k2*r2)):0.0;
        indx++;
      }
    }

    if (DEBUG_LEVEL>2) {
      ImageProcessor ip_ellipse = new FloatProcessor(rpsf_size,rpsf_size);
      float [] ellipsePixels = new float [rpsf_size*rpsf_size];
      indx=0;
      for (iy=0;iy<rpsf_size;iy++) {
//        y=(iy-rpsf_size/2)*decimation+ellipse_coeff[1];  // scale to the original psf (and ellipse_coeff), move center opposite to that of direct kernel (psf)
        y=iy-rpsf_size/2 +ellipse_coeff[1];  //  move center opposite to that of direct kernel (psf)
        for (ix=0;ix<rpsf_size;ix++) {
//          x=(ix-rpsf_size/2)*decimation+ellipse_coeff[0]; // scale to the original psf (and ellipse_coeff), move center opposite to that of direct kernel (psf)
          x=ix-rpsf_size/2 +ellipse_coeff[0]; // move center opposite to that of direct kernel (psf)
          r2=ellipse_coeff[2]*x*x+ellipse_coeff[3]*y*y+ellipse_coeff[4]*x*y;
          m=Math.exp(-k2*r2);
          ellipsePixels[indx++]=(float)((m>=min_mask_threshold)?(Math.exp(-k2*r2)):0.0);
        }
      }
      ip_ellipse.setPixels(ellipsePixels);
      ip_ellipse.resetMinAndMax();
      ImagePlus imp_ellipse= new ImagePlus(title+"_RPSF-MASK_"+cutoff_energy+"-"+ellipse_scale, ip_ellipse);
      imp_ellipse.show();
    }
    return masked_rpsf;
 }

 public double [] maskReversePSFKernel( double []rpsf_pixels, // reversed psf, square array
                                     double [] ellipse_coeff, // ellipse coefficients from _direct_ kernel
                                        double ellipse_scale,
                                   double min_mask_threshold) // zero output element if elliptical Gauss mask is below this threshold
 {
    int rpsf_size=(int)Math.sqrt(rpsf_pixels.length);
    double [] masked_rpsf=new double[rpsf_size*rpsf_size];
    int ix,iy;
    double x,y,r2;
    int indx=0;
    double k2=1/ellipse_scale/ellipse_scale;
    double m;
    for (iy=0;iy<rpsf_size;iy++) {
      y=iy-rpsf_size/2+ellipse_coeff[1];  // move center opposite to that of direct kernel (psf)
      for (ix=0;ix<rpsf_size;ix++) {
        x=ix -rpsf_size/2 +ellipse_coeff[0]; //  move center opposite to that of direct kernel (psf)
        r2=ellipse_coeff[2]*x*x+ellipse_coeff[3]*y*y+ellipse_coeff[4]*x*y;
        m=Math.exp(-k2*r2);
        masked_rpsf[indx]=(m>=min_mask_threshold)?(rpsf_pixels[indx]*Math.exp(-k2*r2)):0.0;
        indx++;
      }
    }
    return masked_rpsf;
 }

 public void showMmaskReversePSFKernel( int rpsf_size, // size of kernel to show
                                     double [] ellipse_coeff, // ellipse coefficients from direct kernel
                                        double ellipse_scale,
                                   double min_mask_threshold, // zero output element if elliptical Gauss mask is below this threshold
                                                String title)
 {
    int ix,iy;
    double x,y,r2;
    int indx=0;
    double k2=1/ellipse_scale/ellipse_scale;
    double m;
    ImageProcessor ip_ellipse = new FloatProcessor(rpsf_size,rpsf_size);
    float [] ellipsePixels = new float [rpsf_size*rpsf_size];
    indx=0;
    for (iy=0;iy<rpsf_size;iy++) {
      y=iy-rpsf_size/2 +ellipse_coeff[1];  //  move center opposite to that of direct kernel (psf)
      for (ix=0;ix<rpsf_size;ix++) {
        x=ix-rpsf_size/2 +ellipse_coeff[0]; // move center opposite to that of direct kernel (psf)
        r2=ellipse_coeff[2]*x*x+ellipse_coeff[3]*y*y+ellipse_coeff[4]*x*y;
        m=Math.exp(-k2*r2);
        ellipsePixels[indx++]=(float)((m>=min_mask_threshold)?(Math.exp(-k2*r2)):0.0);
      }
    }
    ip_ellipse.setPixels(ellipsePixels);
    ip_ellipse.resetMinAndMax();
    ImagePlus imp_ellipse= new ImagePlus(title+"_RPSF-MASK_"+ellipse_scale, ip_ellipse);
    imp_ellipse.show();
 }





/* find ellipse approximating section of the PSF, scale ellipse and use it as a mask to remove PSF far wings */
 public double [] cutPSFWings (double [] psf_pixels, // direct PSF function, square array, may be proportionally larger than reversed
                            double cutoff_energy, // fraction of energy in the pixels to be used
                            double ellipse_scale,
                       double min_mask_threshold, // zero output element if elliptical Gauss mask is below this threshold
                                    String title)
 {
    int psf_size=(int)Math.sqrt(psf_pixels.length);
    double [] masked_psf=new double[psf_size*psf_size];
    int  [][]selection=   findClusterOnPSF(psf_pixels, cutoff_energy, title);
    double [] ellipse_coeff=findEllipseOnPSF(psf_pixels,  selection,    title);
    int ix,iy;
    double x,y,r2;
    int indx=0;
    double k2=1/ellipse_scale/ellipse_scale;
    double m;

    for (iy=0;iy<psf_size;iy++) {
      y=(iy-psf_size/2)-ellipse_coeff[1];  // scale to the original psf (and ellipse_coeff)
      for (ix=0;ix<psf_size;ix++) {
        x=(ix-psf_size/2)-ellipse_coeff[0]; // scale to the original psf (and ellipse_coeff)
        r2=ellipse_coeff[2]*x*x+ellipse_coeff[3]*y*y+ellipse_coeff[4]*x*y;
        m=Math.exp(-k2*r2);
        masked_psf[indx]=(m>=min_mask_threshold)?(psf_pixels[indx]*Math.exp(-k2*r2)):0.0;
        indx++;
      }
    }

    if (DEBUG_LEVEL>2) {
      ImageProcessor ip_ellipse = new FloatProcessor(psf_size,psf_size);
      float [] ellipsePixels = new float [psf_size*psf_size];
      indx=0;
      for (iy=0;iy<psf_size;iy++) {
        y=(iy-psf_size/2)+ellipse_coeff[1];  // scale to the original psf (and ellipse_coeff), move center opposite to that of direct kernel (psf)
        for (ix=0;ix<psf_size;ix++) {
          x=(ix-psf_size/2)+ellipse_coeff[0]; // scale to the original psf (and ellipse_coeff), move center opposite to that of direct kernel (psf)
          r2=ellipse_coeff[2]*x*x+ellipse_coeff[3]*y*y+ellipse_coeff[4]*x*y;
          m=Math.exp(-k2*r2);
          ellipsePixels[indx++]=(float)((m>=min_mask_threshold)?(Math.exp(-k2*r2)):0.0);
        }
      }
      ip_ellipse.setPixels(ellipsePixels);
      ip_ellipse.resetMinAndMax();
      ImagePlus imp_ellipse= new ImagePlus(title+"_PSFWINGS-MASK_"+cutoff_energy+"-"+ellipse_scale, ip_ellipse);
      imp_ellipse.show();
    }
    return masked_psf;
 }







 /* finds cluster on the PSF (with the center at specidfied point)  by flooding from the specified center, so total energy is cutoff_energy fraction
     returns integer array (same dimensions as input) with 1 - selected, 0 - not selected
     cutoff_energy: if positive - specifies fraction of total energy, if negative -cutoff_energy is the minimal value of the pixel to be included 
     UPDATE: follows gradient from the start point to a local maximum if "cutoff_energy" is negative" */
 private int [][] findClusterOnPSF(double [] psf, // PSF function, square array
                                   double cutoff_energy, // fraction of energy in the pixels to be used
                                   String title) {
    int size=(int) Math.sqrt(psf.length);
    return findClusterOnPSF(psf, // PSF function, square array
                            cutoff_energy, // fraction of energy in the pixels to be used
                            size/2, // X0
                            size/2, // Y0
                            title);
 }



 private int [][] findClusterOnPSF(double [] psf, // PSF function, square array
                                   double cutoff_energy, // fraction of energy in the pixels to be used (or minimal level if it is negative)
                                   int startX,  // location of a start point, x-coordinate
                                   int startY,  // location of a start point, y-coordinate
                                   String title) {
    int i,j;
    int ix,iy,ix1,iy1,maxX, maxY;
    List <Integer> pixelList=new ArrayList<Integer>(100);
    Integer Index;
    int size=(int) Math.sqrt(psf.length);
    int [][]clusterMap=new int[size][size];
    double full_energy=0.0;
    int [][] dirs={{-1,0},{-1,-1},{0,-1},{1,-1},{1,0},{1,1},{0,1},{-1,1}};
    ix=startX;
    iy=startY;
    Index=iy*size + ix;
    double maxValue=psf[Index];
/* Make ix,iy to start from the maximal value on PSF */
    Index=0;
    for (i=0;i<size;i++) for (j=0;j<size;j++) {
      full_energy+=psf[Index];
      clusterMap[i][j]=0;
      if (psf[Index]>maxValue){
        maxValue=psf[Index];
        ix=j;
        iy=i;
      }
      Index++;
    }
    double threshold=full_energy*((cutoff_energy>0)?cutoff_energy:1.0); // no limit for negative values of cutoff_energy
    double minValue=0.0; // no limit if total energy is controlled
    double cluster_energy=0.0;
    int clusterSize=0;
    boolean noNew=true;
    if (cutoff_energy<0) { // find nearest local maximum following gradient
      ix=startX;
      iy=startY;
      maxValue=psf[iy*size + ix];
      for (noNew=false;noNew==false;){
        noNew=true;
        for (j=0;j<dirs.length;j++) if (((iy > 0 )        || (dirs[j][1]>=0)) &&
                                        ((iy < (size-1) ) || (dirs[j][1]<=0)) &&
                                        ((ix > 0 )        || (dirs[j][0]>=0)) &&
                                        ((ix < (size-1) ) || (dirs[j][0]<=0))){
          ix1= ix+dirs[j][0];
          iy1= iy+dirs[j][1];
          if (psf[iy1*size+ix1]>maxValue) {
            noNew=false;
            maxValue= psf[iy1*size+ix1];
            ix=ix1;
            iy=iy1;
           break;
          }
        }
      }
      minValue=maxValue*(-cutoff_energy);
    }

    maxX=0;
    maxY=0;
    int listIndex;
    Index=iy*size + ix;
    pixelList.clear();
    pixelList.add (Index);
    clusterSize++;
    clusterMap[iy][ix]=1;
    cluster_energy+=psf[Index];
    noNew=true;
    while ((pixelList.size()>0) &&  (cluster_energy<threshold)) { // will break from the loop if  (psf[Index] <minValue)
/* Find maximal new neighbor */
      maxValue=0.0;
      listIndex=0;
      while (listIndex<pixelList.size()) {
        Index=pixelList.get(listIndex);
        iy=Index/size;
        ix=Index%size;
        noNew=true;
        for (j=0;j<8;j++) if (((iy > 0 ) || (dirs[j][1]>=0)) && ((iy < (size-1) ) || (dirs[j][1]<=0))){
          ix1=(ix+dirs[j][0]+size) % size;
          iy1= iy+dirs[j][1];
          if (clusterMap[iy1][ix1]==0) {
            noNew=false;
            if (psf[iy1*size+ix1]>maxValue) {
              maxValue= psf[iy1*size+ix1];
              maxX=ix1;
              maxY=iy1;
            }
          }
        }
        if (noNew) pixelList.remove(listIndex);  //  remove current list element
        else       listIndex++;     // increase list index
      }
      if (maxValue==0.0) { // Should
        System.out.println("findClusterOnPSF: - should not get here - no points around >0, and threshold is not reached yet.");
        break;
      }
/* Add this new point to the list */
      if (psf[Index]<minValue) break; // break if the condition was value, not total energy
      Index=maxY*size + maxX;
      pixelList.add (Index);
      clusterSize++;
      clusterMap[maxY][maxX]=1;
      cluster_energy+=psf[Index];
      
    } // end of while ((pixelList.size()>0) &&  (cluster_energy<threshold))
    if (DEBUG_LEVEL>3)   System.out.println("findClusterOnPSF: cluster size is "+clusterSize);
    if (DEBUG_LEVEL>6) {
      ImageProcessor ip2 = new FloatProcessor(size,size);
      float [] floatPixels = new float [size*size];
      for (i=0;i<floatPixels.length;i++) {
        floatPixels[i]=(float) psf[i];
      }
      ip2.setPixels(floatPixels);
      ip2.resetMinAndMax();
      ImagePlus imp2= new ImagePlus(title+"_PSF1_"+cutoff_energy, ip2);
      imp2.show();
    }
    if (DEBUG_LEVEL>5) {
      ImageProcessor ip = new FloatProcessor(size,size);
      float [] floatPixels = new float [size*size];
      for (i=0;i<floatPixels.length;i++) {
        floatPixels[i]=(float) clusterMap[i/size][i%size];
      }
      ip.setPixels(floatPixels);
      ip.resetMinAndMax();
      ImagePlus imp= new ImagePlus(title+"_PSF-SEL_"+cutoff_energy, ip);
      imp.show();
    }
    return clusterMap;
}

 /* calculates ellipse (with the center at DC) that interpolates area of the points defined by flooding from the initial centert,
     so total energy is cutoff_energy fraction
     returns {x0,y0,a,b,c} , where a*x^2+b*y^2 + c*x*y=r^2 , so r^2 can be used for a window that removes high far pixels*/

 private double [] findEllipseOnPSF(double [] psf,   // Point Spread Function (may be off-center)
                                    int    [][] selection, // 0/1 - selected/not selected
                                    String title) {
    int i,j;
    double x,y;
    int size=(int) Math.sqrt(psf.length);
    double SX=0.0;
    double SY=0.0;
    double SX2=0.0;
    double SY2=0.0;
    double SXY=0.0;
    double S0=0.0;
    double d,k;
    double area=0; // selection area
/* find centyer */

    for (i=0;i<size;i++) {
      y=i-size/2;
      for (j=0;j<size;j++) if (selection[i][j]>0){
        x=j-size/2;
        d=psf[i*size+j];
        S0+=d;
        SX+=x*d;
        SY+=y*d;
        area+=1.0;
      }
    }
    double centerX=SX/S0;
    double centerY=SY/S0;
    if (DEBUG_LEVEL>5) {
       System.out.println("findEllipseOnPSF: title="+title+" area="+area+" S0="+S0+" SX="+SX+" SY="+SY+" centerX="+centerX+" centerY="+centerY);
    }

/* second puss (could all be done in a single) */
    SX2=0.0;
    SY2=0.0;
    SXY=0.0;
    for (i=0;i<size;i++) {
      y=i-size/2-centerY;
      for (j=0;j<size;j++) if (selection[i][j]>0){
        x=j-size/2-centerX;
        d=psf[i*size+j];
        SX2+=x*x*d;
        SY2+=y*y*d;
        SXY+=x*y*d;
      }
    }
    if (DEBUG_LEVEL>5) {
       System.out.println("findEllipseOnPXF: title="+title+" SX2="+SX2+" SY2="+SY2+" SXY="+SXY);
    }
    k=Math.PI*Math.PI/(2.0*S0*area*area);
    double [] result = {centerX,centerY,k*SY2,k*SX2,-2*k*SXY};
    if (DEBUG_LEVEL>3) {
       System.out.println("findEllipseOnPS: title="+title+" x0="+result[0]+" y0="+result[1]+" a="+result[2]+" b="+result[3]+" c="+result[4]);
    }
    return result;
 }


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */




 /* finds cluster (with the center at DC)  by flooding from DC, so total energy is cutoff_energy fraction
     returns integer array (same dimensions as input) with 1 - selected, 0 - not selected */
 private int [][] findClusterOnPS(double [][] ps, // half power spectrum, starting from 0.0 (DC)
                                   double cutoff_energy, // fraction of energy in the pixels to be used
                                   String title) {
    int i,j;
    List <Integer> pixelList=new ArrayList<Integer>(100);
    Integer Index;
    int size=ps[0].length;
    int [][]clusterMap=new int[size/2+1][size];
    double full_energy=0.0;
    int [][] dirs={{-1,0},{-1,-1},{0,-1},{1,-1},{1,0},{1,1},{0,1},{-1,1}};
    for (i=0;i<(size/2+1);i++) for (j=0;j<size;j++) {
      full_energy+=((i%(size/2))==0)?ps[i][j]:(2*ps[i][j]); /* first and last line are counted once, others - twice */
      clusterMap[i][j]=0;
    }
    double threshold=full_energy*cutoff_energy;
    double cluster_energy=0.0;
    double maxValue;
    int ix,iy,ix1,iy1,maxX, maxY;
    int clusterSize=0;
    ix=0;
    iy=0;
    maxX=0;
    maxY=0;
    int listIndex;
    Index=iy*size + ix;
    pixelList.clear();
    pixelList.add (Index);
    clusterSize++;
    clusterMap[iy][ix]=1;
    cluster_energy+=ps[iy][ix];
    boolean noNew=true;
    while ((pixelList.size()>0) &&  (cluster_energy<threshold)) {
/* Find maximal new neighbor */
      maxValue=0.0;
      listIndex=0;
      while (listIndex<pixelList.size()) {
        Index=pixelList.get(listIndex);
        iy=Index/size;
        ix=Index%size;
        noNew=true;
        for (j=0;j<8;j++) if (((iy > 0 ) || (dirs[j][1]>=0)) && ((iy < (size/2) ) || (dirs[j][1]<=0))){
          ix1=(ix+dirs[j][0]+size) % size;
          iy1= iy+dirs[j][1];
          if (clusterMap[iy1][ix1]==0) {
            noNew=false;
            if (ps[iy1][ix1]>maxValue) {
              maxValue= ps[iy1][ix1];
              maxX=ix1;
              maxY=iy1;
            }
          }
        }
        if (noNew) pixelList.remove(listIndex);  //  remove current list element
        else       listIndex++;     // increase list index
      }
      if (maxValue==0.0) { // Should
        System.out.println("findClusterOnPS: - should not get here - no points around >0, and threshold is not reached yet.");
        break;
      }
/* Add this new point to the list */
      Index=maxY*size + maxX;
      pixelList.add (Index);
      clusterSize++;
      clusterMap[maxY][maxX]=1;
      cluster_energy+=((maxY%(size/2))==0)?ps[maxY][maxX]:(2*ps[maxY][maxX]);
    } // end of while ((pixelList.size()>0) &&  (cluster_energy<threshold))
    if (DEBUG_LEVEL>3)   System.out.println("findClusterOnPS: cluster size is "+clusterSize);
    if (DEBUG_LEVEL>6) {
      ImageProcessor ip2 = new FloatProcessor(size,size/2+1);
      float [] floatPixels = new float [size*(size/2+1)];
      for (i=0;i<floatPixels.length;i++) {
        floatPixels[i]=(float) ps[i/size][i%size];
      }
      ip2.setPixels(floatPixels);
      ip2.resetMinAndMax();
      ImagePlus imp2= new ImagePlus(title+"_PS1_"+cutoff_energy, ip2);
      imp2.show();
    }
    if (DEBUG_LEVEL>6) {
      ImageProcessor ip1 = new FloatProcessor(size,size);
      float [] floatPixels = new float [size*size];
      for (i=0;i<floatPixels.length;i++) {
        iy=i/size-size/2;
        ix=i%size-size/2;
        if (iy<0) {
          ix=-ix;
          iy=-iy;
        }
        ix= (ix+size) % size;
        floatPixels[i]=(float) ps[iy][ix];
      }
      ip1.setPixels(floatPixels);
      ip1.resetMinAndMax();
      ImagePlus imp1= new ImagePlus(title+"_PS_"+cutoff_energy, ip1);
      imp1.show();
    }

    if (DEBUG_LEVEL>5) {
      ImageProcessor ip = new FloatProcessor(size,size);
      float [] floatPixels = new float [size*size];
      for (i=0;i<floatPixels.length;i++) {
        iy=i/size-size/2;
        ix=i%size-size/2;
        if (iy<0) {
          ix=-ix;
          iy=-iy;
        }
        ix= (ix+size) % size;
        floatPixels[i]=(float) clusterMap[iy][ix];
      }
      ip.setPixels(floatPixels);
      ip.resetMinAndMax();
      ImagePlus imp= new ImagePlus(title+"_SEL_"+cutoff_energy, ip);
      imp.show();
    }
    return clusterMap;
}

 /* calculates ellipse (with the center at DC) that interpolates area of the points defined by flooding from DC, so total energy is cutoff_energy fraction
     returns {a,b,c} , where a*x^2+b*y^2 + c*x*y=r^2 , so r^2 can be used for a window that removes high frequancy components that are too low to be useful*/

 private double [] findEllipseOnPS(double [][] ps,   // half power spectrum, starting from 0.0 (DC)
                                   int    [][] selection, // 0/1 - selected/not selected
                                   String title) {
    int i,j;
    double x,y;
    int size=ps[0].length;
    double SX2=0.0;
    double SY2=0.0;
    double SXY=0.0;
    double S0=0.0;
    double k=2.0;
    double d;
    double area=0; // selection area
    for (i=0;i<(size/2+1);i++) {
      k=((i%(size/2))==0)?1.0:2.0;
      y=i;
      for (j=0;j<size;j++) if (selection[i][j]>0){
        x=(j>(size/2))?(j-size):j;
        d=k*ps[i][j];
        S0+=d;
        SX2+=x*x*d;
        SY2+=y*y*d;
        SXY+=x*y*d;
        area+=1.0;
      }
    }
    if (DEBUG_LEVEL>5) {
       System.out.println("findEllipseOnPS: title="+title+" area="+area+" S0="+S0+" SX2="+SX2+" SY2="+SY2+" SXY="+SXY);
    }
//    k=Math.PI*Math.PI/(2.0*S0*S0*area*area);
//    double [] result = {k*SY2,k*SX2,2*k*SXY};
    k=Math.PI*Math.PI/(2.0*S0*area*area);
    double [] result = {k*SY2,k*SX2,-2*k*SXY};
    if (DEBUG_LEVEL>3) {
       System.out.println("findEllipseOnPS: title="+title+" a="+result[0]+" b="+result[1]+" c="+result[2]);
    }
    return result;
 }

/* Trying to remove aliasing artifacts when the decimated (pixel resolution) image is deconvolved with full resolution (sub-pixel resolution)
    model pattern. This effect is also easily visible if the decimated model is deconvolved with the same one art full resolution.
    Solution is to clone the power spectrum of the full resolution model with the shifts to match oversampling (15 clones for the 4x oversampling),
    And add them together (adding also zero frequerncy point - it might be absent o0n the model) but not include the original (true one) and
    use the result to create a rejectiobn mask - if the energy was high, (multiplicative) mask should be zero at those points. */

  public double [][] maskAliases (double [][][] cmplx, // complex spectrum, [size/2+1][size]
                                      boolean checker, // checkerboard pattern in the source file (use when filtering)
                                       int oversample,  // measured array is sampled at 1/oversample frequency than model (will add more parameters later)
                                 double zerofreq_size,  // add rejection of zero frequency (~2-3pix)
                                double threshold_high,  // reject completely if energy is above this part of maximal
                                 double threshold_low){  // leave intact if energy is below this part of maximal
//    double th=Math.sqrt(threshold_high);
//    double tl=Math.sqrt(threshold_low);
    double th=threshold_high*threshold_high;
    double tl=threshold_low*threshold_low;
    int size=cmplx[0].length;
    double [][] ps=new double [size/2+1][size];
    int i,j, zs ,ix,iy, cloneNx, cloneNy, cloneX, cloneY, n;
    int cloneStep=size/oversample;
    double psMax=0.0;

    if (DEBUG_LEVEL>3) System.out.println("maskAliases(cmplx,"+checker+","+oversample+","+zerofreq_size+","+threshold_high+","+threshold_low);


/* generating power spectrum for the high-res complex spectrum, find maximum value and normalize */
    for (i=0;i<=size/2; i++) for (j=0;j<size; j++) {
      ps[i][j]=cmplx[i][j][0]*cmplx[i][j][0]+cmplx[i][j][1]*cmplx[i][j][1];
      if (ps[i][j]>psMax) psMax=ps[i][j];
    }
    double k=1/psMax;
    for (i=0;i<=size/2; i++) for (j=0;j<size; j++) ps[i][j]*=k;

/* Add maximum at (0,0) */
    if (zerofreq_size>0.0) {
      zs=(int) (4*zerofreq_size);
      k=0.5/(zerofreq_size*zerofreq_size);
      if (zs>ps.length) zs =ps.length;
      for (iy=0;iy<=zs;iy++) for (ix=-zs; ix <= zs; ix++) {
        j=(ix+size)%size;
        ps[iy][j]+=Math.exp(-k*(iy*iy+ix*ix));
      }
    }
    double [][] mask=new double [size/2+1][size];
    for (i=0;i<=size/2; i++) for (j=0;j<size; j++) mask[i][j]=0.0;
/* clone spectrums */
    int dup=1; // duplicate (mirror) points on all lines but first and last
    for (i=0;i<=size/2;i++) {
       dup=((i>0) && (i!=(size/2)))?2:1;
       for (j=0;j<size;j++) for (n=0;n<dup;n++){
         iy=(n>0)?(-i):i;
         ix=(n>0)?(-(j>(size/2)?(j-size):j))  :  (j>(size/2)?(j-size):j);
         for (cloneNy=0;cloneNy<oversample;cloneNy++) for (cloneNx=0;cloneNx<oversample;cloneNx++)
           if (((cloneNy!=0) || (cloneNx!=0)) && // not a zero point
               (!checker ||                      // use all if it is not a checkerboard pattren
                (((cloneNx ^ cloneNy) & 1)==0) )) { // remove clones in a checker pattern
           cloneY=iy+cloneNy*cloneStep;
           cloneY=(cloneY+size)%size;
           if (cloneY<=(size/2)) {// only calculate half
              cloneX=ix+cloneNx*cloneStep;
              cloneX=(cloneX+size)%size;
              mask[cloneY][cloneX]+=ps[i][j];
           }
         }
       }
    }

/* debug show the mask */
      if (DEBUG_LEVEL>2) {
        float [] fpixelsMask = new float [(size/2+1)*size];
        for (i=0;i<=size/2;i++) for (j=0; j<size;j++) fpixelsMask[i*size+j]=(float) mask[i][j];
        ImageProcessor ip_mask = new FloatProcessor(size,size/2+1);
        ip_mask.setPixels(fpixelsMask);
        ip_mask.resetMinAndMax();
        ImagePlus imp_mask= new ImagePlus("MASK", ip_mask);
        imp_mask.show();
      }
/*
      if (DEBUG_LEVEL>2) {
        float [] fpixels = new float [(size/2+1)*size];
        for (i=0;i<=size/2;i++) for (j=0; j<size;j++) fpixels[i*size+j]=(float) ps[i][j];
        ImageProcessor ip1 = new FloatProcessor(size,size/2+1);
        ip1.setPixels(fpixels);
        ip1.resetMinAndMax();
        ImagePlus imp1= new ImagePlus("PS0-"+zerofreq_size+"-"+threshold_low+"-"+threshold_high, ip1);
        imp1.show();
      }
      if (DEBUG_LEVEL>10) { 
        float [] maskPixels = new float [size*size];
        for (i=0;i<=size/2;i++) {
          iy=(i+size/2)%size;
          for (j=0; j<size;j++) {
            ix=(j+size/2)%size;
            maskPixels[iy*size+ix]=(float) ps[i][j];
            if ((i>0) && (iy>0)) {
              maskPixels[((size-iy)%size)*size+((size-ix)%size)]=(float) ps[i][j];
            }

          }
        }
        ImageProcessor ip_fhtMask = new FloatProcessor(size,size);
        ip_fhtMask.setPixels(maskPixels);
        ip_fhtMask.resetMinAndMax();
        ImagePlus imp_fhtMask= new ImagePlus("PS-"+zerofreq_size+"-"+threshold_low+"-"+threshold_high, ip_fhtMask);
        imp_fhtMask.show();
      }
*/


/* make mask of cloned power spectrums */
    double a;
    for (i=0;i<=size/2;i++) for (j=0;j<size;j++) {
      if      (mask[i][j]<tl)  mask[i][j]=1.0;
      else if (mask[i][j]>th) mask[i][j]=0.0;
      else { // make smooth transition
        a=(2.0 * mask[i][j] - th - tl)/(th - tl);
        mask[i][j]=0.5*(1.0-a*a*a);
      }
    }
    return mask;
  }

  public double [][] calcMeshFFTCorrectionMap(double [][][] cmplx, // complex spectrum, [size/2+1][size]
                                                 double threshold, // if 0 use normalize amplitude, if 0..1 - make binary: 1.0 if > threshold, 0.0 - otherwise
                                                    double radius, // low-pass result with low pass filter (should be later defined automatically)
                                                       int  hsize,  // 2d histogram size (size/2 probably a good guess),
                                                double percentile, // use this percentile (0.0..1.0)) value for given radius as a target
                                                   double maxGain, // maximal gain for low components
                                                 double exaggerate, // exaggerate correction mask with Math.pow()) 
                                                    String title ){
    int i,j,indx,ir,ir1,iv,ix,iy;
    double max,d,r,d0;
    int size=cmplx[0].length;
    float [] fmask= new float [size*size];
    indx=0;
    max=0;
    for (i=0;i<(size/2+1);i++) for (j=0;j<size;j++) {
      d=Math.sqrt(cmplx[i][j][0]*cmplx[i][j][0]+cmplx[i][j][1]*cmplx[i][j][1]);
      fmask[indx++] = (float) d;
      if (max <d ) max=d;
    }
    if (threshold ==0)  for (i=0;i<indx;i++) fmask[i]/= max;
    else {
      max*=threshold;
      for (i=0;i<indx;i++) fmask[i]= (float) ((fmask[i]>max)?1.0:0.0);
    }
    for (i=1;i<fmask.length/2;i++) fmask[fmask.length-i]=fmask[i];
/* low-pass filter by FFT+window */
    ImageProcessor ip_mask = new FloatProcessor(size,size);
    ip_mask.setPixels(fmask);
    if (DEBUG_LEVEL>7) {
      ip_mask.resetMinAndMax();
      ImagePlus imp_mask=  new ImagePlus(title+"-mask_before_filter-"+threshold, ip_mask);
      imp_mask.show();
    }
    FHT fht_mask =  new FHT(ip_mask);
// No need to swap quadrants
    fht_mask.transform();
    float [] fpixels = (float []) fht_mask.getPixels();
/* multiply by mask */
    double [] gauss= new double [size];
    for (i=0;i<=size/2;i++) {
      d=i/radius;
      gauss[i]=Math.exp(-d*d);
      if (i>0)  gauss[size-i]=gauss[i];
    }
    indx=0;
    for (i=0;i<size;i++) for (j=0;j<size;j++) {
      d=gauss[i]*gauss[j];
      fpixels[indx]*=(float)d;
      indx++;
    }
    fht_mask.setPixels (fpixels);
/* transform back */
    fht_mask.inverseTransform();
    fmask=(float []) fht_mask.getPixels();
    if (DEBUG_LEVEL>7) {
      fht_mask.resetMinAndMax();
      ImagePlus imp_mask_after=  new ImagePlus(title+"-mask_afer_filter-"+threshold, fht_mask);
      imp_mask_after.show();
    }

/* extract filtered smooth mask - shold have "X"-shape */
    double [][] dmask = new double [size/2+1][size];
    indx=0;
    d=1.0/fmask[0];
    for (i=0;i<(size/2+1);i++) for (j=0;j<size;j++) dmask[i][j]=fmask[indx++]*d;


    if ((DEBUG_LEVEL>3) && (title!="")) { /* Increase debug level later */
        float [] maskPixels0 = new float [size*size];
//        int ix,iy;
        for (i=0;i<=size/2;i++) {
          iy=(i+size/2)%size;
          for (j=0; j<size;j++) {
            ix=(j+size/2)%size;
            maskPixels0[iy*size+ix]=(float) dmask[i][j];
            if ((i>0) && (iy>0)) {
              maskPixels0[((size-iy)%size)*size+((size-ix)%size)]=(float) dmask[i][j];
            }
          }
        }
        ImageProcessor ip_rsltMask0 = new FloatProcessor(size,size);
        ip_rsltMask0.setPixels(maskPixels0);
        ip_rsltMask0.resetMinAndMax();
        ImagePlus imp_rsltMask0= new ImagePlus(title+"_RM0-"+threshold+"-"+radius+"-"+hsize+"-"+percentile+"-"+maxGain, ip_rsltMask0);
        imp_rsltMask0.show();
    }




/* Build a 2-d histogram for value vs. radius. Use radius step of 1 and values from 0 to max in (size/2 steps) */
    double [][] hist2d=new double [hsize][hsize];
    for (i=0;i<hsize;i++) for (j=0;j<hsize;j++) hist2d[i][j]=0.0;
    double kr=2.0*hsize/size;
    for (i=0;i<=(size/2);i++) for (j=0;j<size;j++){
      ix=j-((j>size/2)?size:0);
      ir=(int) Math.round(kr*Math.sqrt(i*i+ix*ix));
      if (ir<0) ir=0;
      else if (ir>=hsize) ir=hsize-1;
      iv= (int) Math.round (dmask[i][j]*(hsize-1));
      if (iv<0) iv=0;
      else if (iv>=hsize) iv=hsize-1;
      hist2d[ir][iv]+=((i>0)&&(i<=(size/2)))?2.0:1.0;
//     if (ir==1)  System.out.println("kr="+kr+" ir="+ir+" iv="+iv+ " i="+i+" j="+j+" ix="+ix);

    }
/* Smooth histogram ???*/    
/// optionally show the histogram
    if (DEBUG_LEVEL>3) {
      float [] fhist=new float [hsize*hsize];
      indx=0;
      for (i=0;i<hsize;i++) for (j=0;j<hsize;j++) fhist[indx++]= (float) Math.pow(hist2d[i][j],0.2);
      ImageProcessor ip_hist = new FloatProcessor(hsize,hsize);
      ip_hist.setPixels(fhist);
      ip_hist.resetMinAndMax();
      ImagePlus imp_hist= new ImagePlus(title+"_histogram", ip_hist);
      imp_hist.show();
    }



/* calculate target value vs. radius function using percentile value and a 2-d histogram. First radius with o bin having more than
    a percentile limits radius - all components farther will be masked out */
    double [] v_vs_r=new double [hsize];
    for (i=0;i<hsize;i++) v_vs_r[i]=0.0;
    for (i=0;i<hsize;i++) {
      d=0.0;
      for (j=0;j<hsize;j++) d+=hist2d[i][j];
      max=d*percentile;
      d=0.0;
      d0=0.0;
      for (j=0;(j<hsize) && (d<max);j++ ) {
         d0=d;
         d+=hist2d[i][j];
      }
      j--;
      if (j<1) v_vs_r[i]=0.0;
      else if ((d0==0.0) && (j==(hsize-1))) v_vs_r[i]=1.0;
      else                                  v_vs_r[i]=1.0*(j-1)/(hsize-1) + ((d>d0)?((max-d0)/(d-d0)/(hsize-1)):0.0);
//      System.out.println("(v_vs_r["+i+"]="+v_vs_r[i]+ " i="+i+" j="+j+" d="+d+" d0="+d0+" max="+max);
      if (j==0) break; // or just ==0.0? - will leave values for higher r equal to 0.0 (as initialized))
//      System.out.println("(v_vs_r["+i+"]="+v_vs_r[i]);
    }
//    System.out.println("(v_vs_r["+(i-1)+"] was the last one, broke from the cycle");

/* calculate mask array itself */
    double rmaxGain2=1/maxGain/maxGain;
    for (i=0;i<=(size/2);i++) for (j=0;j<size;j++){
      ix=j-((j>size/2)?size:0);
      r=kr*Math.sqrt(i*i+ix*ix);
      ir=(int) Math.floor(r); /// low for interpolation
      ir1=ir+1; // high for the interpolation
      if (ir<0) ir=0;
      else if (ir>=hsize) ir=hsize-1;
/* linear interpolate value for this r using table */
      d=v_vs_r[ir];
      if ((ir1>r) && (ir1<hsize)) {
        d+=(v_vs_r[ir1]-d)*(r-ir);
      }
/* calculate mask value (and limit by maxGain)*/
      if (d>0){
        if (dmask[i][j] >0.0) {
/* variant with possible sharp edges */
//        dmask[i][j]=d/dmask[i][j];
//        if (dmask[i][j]> maxGain)  dmask[i][j]=maxGain;
/* variant with smooth edges */
        d=dmask[i][j]/d;
        dmask[i][j]=d/(rmaxGain2+d*d);
        } else dmask[i][j]=maxGain;
        if (exaggerate!=1.0) dmask[i][j]=Math.pow(dmask[i][j],exaggerate);
      } else dmask[i][j]=0.0;
      
    }
/* optionally show the result mask*/
    if ((DEBUG_LEVEL>2) && (title!="")) { /* Increase debug level later */
        float [] maskPixels = new float [size*size];
//        int ix,iy;
        for (i=0;i<=size/2;i++) {
          iy=(i+size/2)%size;
          for (j=0; j<size;j++) {
            ix=(j+size/2)%size;
            maskPixels[iy*size+ix]=(float) dmask[i][j];
            if ((i>0) && (iy>0)) {
              maskPixels[((size-iy)%size)*size+((size-ix)%size)]=(float) dmask[i][j];
            }
          }
        }
        ImageProcessor ip_rsltMask = new FloatProcessor(size,size);
        ip_rsltMask.setPixels(maskPixels);
        ip_rsltMask.resetMinAndMax();
        ImagePlus imp_rsltMask= new ImagePlus(title+"_RMSK-"+threshold+"-"+radius+"-"+hsize+"-"+percentile+"-"+maxGain, ip_rsltMask);
        imp_rsltMask.show();
    }
    return dmask;
  }

  public double[] limitedInverseOfFHTDiffSize(double [] measuredPixels,
                                                 double [] modelPixels,
                                                   double deconvInvert,
                                                   boolean forward_OTF,
                                                          String title) {
     int size=(int) Math.sqrt(measuredPixels.length);
     int modelSize=(int) Math.sqrt(modelPixels.length);
     int oversample= modelSize/size;
     double [] rslt =   new double [modelSize*modelSize];
     double [] subModel=new double [size*size];
     double [] subRslt= new double [size*size];
     int i,j,k,l,indx;
     for (k=0;k<oversample;k++) for (l=0;l<oversample;l++) {
        indx=0;
        for (i=0;i<size;i++) for (j=0;j<size;j++) {
          subModel[indx++]=modelPixels[modelSize*(oversample*i+k)+ (oversample*j+l)];
        }
        subRslt=limitedInverseOfFHT(measuredPixels,
                                          subModel,
                                              size,
                                             false,  // checkerboard pattern in the source file (use when filtering)
                                      deconvInvert,
                                       forward_OTF,
                                                 1, // int oversample
                                               0.0, // double zerofreq_size
                                               0.0, // OTF_smoothPS,
                                               2.0, // double threshold_high,  // reject completely if energy is above this part of maximal
                                               2.0, // double threshold_low,  // leave intact if energy is below this part of maximal
                                                "");
        indx=0;
        for (i=0;i<size;i++) for (j=0;j<size;j++) {
//          rslt[modelSize*(oversample*i+(oversample-k-1))+ (oversample*j+(oversample-l-1))]=subRslt[indx++];
          rslt[modelSize*(oversample*i+k)+ (oversample*j+l)]=subRslt[indx++];
        }
     }
     if (DEBUG_LEVEL>1) {
       float [] floatPixels=new float[modelSize*modelSize];
       for (i=0;i<floatPixels.length; i++)   floatPixels[i]=(float) rslt[i];
        ImageProcessor ip = new FloatProcessor(modelSize,modelSize);
       ip.setPixels(floatPixels);
//    float [] direct_target=(float[])ip_target.getPixels();
       ip.resetMinAndMax();
       ImagePlus imp=  new ImagePlus(title+"_R_"+deconvInvert, ip);
       imp.show();
    }
    return rslt;
  }


/* TODO: remove normalization on max when dividing */
  public double[] limitedInverseOfFHT(double [] measuredPixels,
                                         double [] modelPixels,
                                                      int size,
                                               boolean checker,  // checkerboard pattern in the source file (use when filtering)
                                           double deconvInvert,
                                           boolean forward_OTF,
                                                 String title) {
      return limitedInverseOfFHT(measuredPixels,
                                    modelPixels,
                                           size,
                                        checker,
                                   deconvInvert,
                                    forward_OTF,
                                              1, // int oversample
                                            0.0, // double zerofreq_size
                                            0.0, // OTF_smoothPS,
                                            2.0, // double threshold_high,  // reject completely if energy is above this part of maximal
                                            2.0, // double threshold_low,  // leave intact if energy is below this part of maximal
                                          title);
  }
// pixels - windowed (around center at [size/2, size/2]) measured pixels (one of the bayer components

  public double[] limitedInverseOfFHT(double [] measuredPixels,  // measured pixel array
                                         double [] modelPixels,  // simulated (model) pixel array)
                                                      int size,  // FFT size
                                               boolean checker,  // checkerboard pattern in the source file (use when filtering)
                                           double deconvInvert,  //  fraction of the maximal value to be used to limit zeros
                                           boolean forward_OTF,  // divide measured by simulated when true, simulated by measured - when false
                                                int oversample,  // measured array is sampled at 1/oversample frequency than model (will add more parameters later)
                                          double zerofreq_size,  // add rejection of zero frequency (~2-3pix)
                                               double smoothPS,       // 0 - none, otherwise Gauss width = FFT size/2/smoothPS
                                         double threshold_high,  // reject completely if energy is above this part of maximal
                                          double threshold_low,  // leave intact if energy is below this part of maximal
                                                 String title) { // title base for optional plots names
      return limitedInverseOfFHT(measuredPixels,
                                    modelPixels,
                                           size,
                                        checker,
                                   deconvInvert,
                                    forward_OTF,
                                     oversample,  // measured array is sampled at 1/oversample frequency than model (will add more parameters later)
                                  zerofreq_size,  // add rejection of zero frequency (~2-3pix)
                                       smoothPS,       // 0 - none, otherwise Gauss width
                                 threshold_high,  // reject completely if energy is above this part of maximal
                                  threshold_low,  // leave intact if energy is below this part of maximal
                                           -1.0, // if 0 use normalize amplitude, if 0..1 - make binary: 1.0 if > threshold, 0.0 - otherwise -1 - disable mask
                                            0.0, // low-pass result with low pass filter (should be later defined automatically)
                                              0,  // 2d histogram size (size/2 probably a good guess),
                                            0.0, // use this percentile (0.0..1.0)) value for given radius as a target
                                              0, // maximal gain for low components
                                            1.0, // exaggerate correction mask with Math.pow()) 
                                          title);
  }

  public double[] limitedInverseOfFHT(double [] measuredPixels,  // measured pixel array
                                         double [] modelPixels,  // simulated (model) pixel array)
                                                      int size,  // FFT size
                                               boolean checker,  // checkerboard pattern in the source file (use when filtering)
                                           double deconvInvert,  //  fraction of the maximal value to be used to limit zeros
                                           boolean forward_OTF,  // divide measured by simulated when true, simulated by measured - when false
                                                int oversample,  // measured array is sampled at 1/oversample frequency than model (will add more parameters later)
                                          double zerofreq_size,  // add rejection of zero frequency (~2-3pix)
                                               double smoothPS,       // 0 - none, otherwise Gauss width = FFT size/2/smoothPS
                                         double threshold_high,  // reject completely if energy is above this part of maximal
                                          double threshold_low,  // leave intact if energy is below this part of maximal
                                              double threshold, // if 0 use normalize amplitude, if 0..1 - make binary: 1.0 if > threshold, 0.0 - otherwise -1 - disable mask
                                                 double radius, // low-pass result with low pass filter (should be later defined automatically)
                                                    int  hsize,  // 2d histogram size (size/2 probably a good guess),
                                             double percentile, // use this percentile (0.0..1.0)) value for given radius as a target
                                                double maxGain, // maximal gain for low components
                                             double exaggerate, // exaggerate correction mask with Math.pow()) 
                                                 String title) { // title base for optional plots names

    double [] denominatorPixels= forward_OTF? modelPixels.clone():    measuredPixels.clone();
    double [] nominatorPixels=   forward_OTF? measuredPixels.clone(): modelPixels.clone();

    double[][][] fft_complex,fft_nominator;
    int i,j,indx;
    double DCLevel=0.0;
    double a,k,r,r2,k2;

    for (i=0;i<denominatorPixels.length; i++)  DCLevel+=denominatorPixels[i];
    DCLevel/=(size*size);
    for (i=0;i<denominatorPixels.length; i++)  denominatorPixels[i]-=DCLevel;
    double [] denominatorPixelsSmooth=denominatorPixels.clone();
    fht_instance.swapQuadrants(denominatorPixels);
    fht_instance.transform(denominatorPixels);
// Convert from FHT to complex FFT
    fft_complex= FHT2FFTHalf (denominatorPixels,size);
    double [][][] fft_complex_smooth=fft_complex;
    if (smoothPS>0.0) {
/* Make a new FFT of the windowed denominator - window would widen PS used for aliases rejection mask maskAliases() */
      double [] gaussWindow = new double[size];
      double gaussK=4.0*smoothPS*smoothPS/size/size;
      for (i=0;i<=size/2; i++) {
         gaussWindow[size/2-i]=Math.exp(-gaussK*i*i);
         if (i<size/2) gaussWindow[size/2+i]=gaussWindow[size/2-i];
      }
      for (indx=0;indx<denominatorPixels.length;indx++) {
        i=indx/size;
        j=indx%size;
        denominatorPixelsSmooth[indx]*=gaussWindow[i]*gaussWindow[j];
      }
      fht_instance.swapQuadrants(denominatorPixelsSmooth);
      fht_instance.transform(denominatorPixelsSmooth);
// Convert from FHT to complex FFT
      fft_complex_smooth= FHT2FFTHalf (denominatorPixelsSmooth,size);
    }
    double [][] mask= null;
    if ((oversample>1) && (threshold_low<1.0)) {
/* create mask */
       mask= maskAliases (fft_complex_smooth,   // complex spectrum, [size/2+1][size]
                                     checker, // checkerboard pattern in the source file (use when filtering)
                                  oversample,   // measured array is sampled at 1/oversample frequency than model (will add more parameters later)
                               zerofreq_size,   // add rejection of zero frequency (~2-3pix)
                              threshold_high,   // reject completely if energy is above this part of maximal
                               threshold_low);  // leave intact if energy is below this part of maximal
/* debug show the mask */
      if ((DEBUG_LEVEL>3) && (title!="")) { /* Increase debug level later */
        float [] maskPixels = new float [size*size];
        int ix,iy;
        for (i=0;i<=size/2;i++) {
          iy=(i+size/2)%size;
          for (j=0; j<size;j++) {
            ix=(j+size/2)%size;
            maskPixels[iy*size+ix]=(float) mask[i][j];
            if ((i>0) && (iy>0)) {
              maskPixels[((size-iy)%size)*size+((size-ix)%size)]=(float) mask[i][j];
            }
          }
        }
        ImageProcessor ip_fhtMask = new FloatProcessor(size,size);
        ip_fhtMask.setPixels(maskPixels);
        ip_fhtMask.resetMinAndMax();
        ImagePlus imp_fhtMask= new ImagePlus(title+"_MASK_"+zerofreq_size+"-"+threshold_low+"-"+threshold_high+"-", ip_fhtMask);
        imp_fhtMask.show();
      }
    }
  double [][] xMask=null;
/* Calculate "X"-like mask to reduce artifacts caused by different representation of the FFT components in the model mesh"*/  
  if (threshold>=0.0) {
    xMask=calcMeshFFTCorrectionMap(fft_complex, // complex spectrum, [size/2+1][size]
                                     threshold, // if 0 use normalize amplitude, if 0..1 - make binary: 1.0 if > threshold, 0.0 - otherwise
                                        radius, // low-pass result with low pass filter (should be later defined automatically)
                                         hsize,  // 2d histogram size (size/2 probably a good guess),
                                    percentile, // use this percentile (0.0..1.0)) value for given radius as a target
                                       maxGain, // maximal gain for low components
                                    exaggerate, // exaggerate correction mask with Math.pow()) 
                                         title );
  }
/// deconvInvert
/// Now tricky thing. Invert Z for large values, but make them Z - for small ones. So it will be a mizture of correlation and deconvolution
// here the targets are round, but what will th\be the colrrect way fo assymmetrical ones?
/// First - find maximal value
    double fft_max=0;
    for (i=0;i<fft_complex.length; i++) for (j=0;j<fft_complex[0].length;j++) {
      r2=fft_complex[i][j][0]*fft_complex[i][j][0]+fft_complex[i][j][1]*fft_complex[i][j][1];
      if (r2>fft_max) fft_max=r2;
    }
    k=Math.sqrt(fft_max)*deconvInvert;
    k2=k*k;
    for (i=0;i<fft_complex.length; i++) for (j=0;j<fft_complex[0].length;j++) {
      r=Math.sqrt(fft_complex[i][j][0]*fft_complex[i][j][0]+fft_complex[i][j][1]*fft_complex[i][j][1]);
      a=-Math.atan2(fft_complex[i][j][1],fft_complex[i][j][0]); /// was zero for circular targets)
      r=r/(r*r+k2); /* still seems to be the best */
//    r=(r<k)?0.0:1/r;   // works OK
//      r=(r<k)?1.0:1/r; // won't work
      fft_complex[i][j][0]=r*Math.cos(a);
      fft_complex[i][j][1]=r*Math.sin(a);
    }
    if ((DEBUG_LEVEL>5) && (title!="")) {
        ImageProcessor ip_fht1 = new FloatProcessor(size,size);
        ip_fht1.setPixels(floatFFTHalf2FHT (fft_complex,size));
        ip_fht1.resetMinAndMax();
        ImagePlus imp_fht1= new ImagePlus(title+"_1/denom_FHT_nomask"+deconvInvert, ip_fht1);
        imp_fht1.show();
    }

/* debug show before masking */
      if ((oversample>1.0) &&(DEBUG_LEVEL>3) && (mask!=null)  && (title!="")) { /* Increase debug level later */
        float [] maskPixels1 = new float [size*size];
        int ix,iy;
        for (i=0;i<=size/2;i++) {
          iy=(i+size/2)%size;
          for (j=0; j<size;j++) {
            ix=(j+size/2)%size;
            a=Math.sqrt(fft_complex[i][j][0]*fft_complex[i][j][0]+fft_complex[i][j][1]*fft_complex[i][j][1]);
            maskPixels1[iy*size+ix]=(float) a;
            if ((i>0) && (iy>0)) {
              maskPixels1[((size-iy)%size)*size+((size-ix)%size)]=(float) a;
            }
          }
        }
        ImageProcessor ip_fhtMask1 = new FloatProcessor(size,size);
        ip_fhtMask1.setPixels(maskPixels1);
        ip_fhtMask1.resetMinAndMax();
        ImagePlus imp_fhtMask1= new ImagePlus(title+"_UNMASKED_"+zerofreq_size+"-"+threshold_low+"-"+threshold_high+"-", ip_fhtMask1);
        imp_fhtMask1.show();
      }

/* multiply 1/denom by the filter mask if defined */
    if (mask!=null) {
      if (DEBUG_LEVEL>2) {
         System.out.println("Multiplying denominator by a mask");
      }

      for (i=0;i<fft_complex.length; i++) for (j=0;j<fft_complex[0].length;j++) {
        fft_complex[i][j][0]*=mask[i][j];
        fft_complex[i][j][1]*=mask[i][j];
      }
    }
/* multiply 1/denom by the xMask mask if defined */
    if (xMask!=null) {
      if (DEBUG_LEVEL>2) {
         System.out.println("Multiplying denominator by xMask");
      }

      for (i=0;i<fft_complex.length; i++) for (j=0;j<fft_complex[0].length;j++) {
        fft_complex[i][j][0]*=xMask[i][j];
        fft_complex[i][j][1]*=xMask[i][j];
      }
    }


/* debug show after masking */
      if ((oversample>1.0) &&(DEBUG_LEVEL>3) && (mask!=null) && (title!="")) { /* Increase debug level later */
        float [] maskPixels2 = new float [size*size];
        int ix,iy;
        for (i=0;i<=size/2;i++) {
          iy=(i+size/2)%size;
          for (j=0; j<size;j++) {
            ix=(j+size/2)%size;
            a=Math.sqrt(fft_complex[i][j][0]*fft_complex[i][j][0]+fft_complex[i][j][1]*fft_complex[i][j][1]);
            maskPixels2[iy*size+ix]=(float) a;
            if ((i>0) && (iy>0)) {
              maskPixels2[((size-iy)%size)*size+((size-ix)%size)]=(float) a;
            }
          }
        }
        ImageProcessor ip_fhtMask2 = new FloatProcessor(size,size);
        ip_fhtMask2.setPixels(maskPixels2);
        ip_fhtMask2.resetMinAndMax();
        ImagePlus imp_fhtMask2= new ImagePlus(title+"_MASKED_"+zerofreq_size+"-"+threshold_low+"-"+threshold_high+"-", ip_fhtMask2);
        imp_fhtMask2.show();
      }
/// if nominatorPixels are defined, convert them and mutiply the result by it (convolve)
    if (nominatorPixels!=null) {
       if (nominatorPixels.length!=denominatorPixels.length) {
         IJ.showMessage("Error","Arrays have different sizes - nominatorPixels.length="+nominatorPixels.length+", denominatorPixels.length="+denominatorPixels.length);
         return null;
       }
/// optionally show the result
       if ((DEBUG_LEVEL>5) && (title!="")) {
        ImageProcessor ip_fht1 = new FloatProcessor(size,size);
        ip_fht1.setPixels(floatFFTHalf2FHT (fft_complex,size));
        ip_fht1.resetMinAndMax();
        ImagePlus imp_fht1= new ImagePlus(title+"_1/denom_FHT_masked"+deconvInvert, ip_fht1);
        imp_fht1.show();
      }

//       floatPixels=new float[nominatorPixels.length];

       for (i=0;i<nominatorPixels.length; i++)  DCLevel+=nominatorPixels[i];
       DCLevel/=(size*size);
       for (i=0;i<nominatorPixels.length; i++)  nominatorPixels[i]-=DCLevel;
// Swapping quadrants, so the center will be 0,0
       fht_instance.swapQuadrants(nominatorPixels);
// get to frequency domain
       fht_instance.transform(nominatorPixels);

// Convert from FHT to complex FFT
       fft_nominator= FHT2FFTHalf (nominatorPixels,size);
// multiply fft_complex by fft_nominator
       for (i=0;i<fft_complex.length; i++) for (j=0;j<fft_complex[0].length;j++) {
         a=                   fft_complex[i][j][0]*fft_nominator[i][j][0]-fft_complex[i][j][1]*fft_nominator[i][j][1];
         fft_complex[i][j][1]=fft_complex[i][j][0]*fft_nominator[i][j][1]+fft_complex[i][j][1]*fft_nominator[i][j][0];
         fft_complex[i][j][0]=a;
       }   
    }

// set fht_target pixels with new values 
    double [] pixels=FFTHalf2FHT (fft_complex,size);
/// transform to space
    fht_instance.inverseTransform(pixels);
    fht_instance.swapQuadrants(pixels);
   return pixels;
  }



  public double correlationContrast ( double [] pixels,       // square pixel array
                                      double [][] wVectors,   // wave vectors (same units as the pixels array)
                                      double ringWidth,       // ring (around r=0.5 dist to opposite corr) width
                                      double x0,              // center coordinates
                                      double y0,
                                      String title) { // title base for optional plots names
    int size=(int) Math.sqrt(pixels.length);
    double [] xy= new double [2];
    double [] uv; 
    double r2,d;
    int i,j;
/* opposite sign correlation points in uv are at uv=(0,-0.5),(0,0.5), (-0.5,0) and (0.5,0), with radius of (1/2)
    selecting center circle and a ring from 0.25 to 0.75 of the distance to opposite sign correlations */
    
    double r2WingsOuter= 0.0625*(1.0+ringWidth)*(1.0+ringWidth);
    double r2WingsInner= 0.0625*(1.0-ringWidth)*(1.0-ringWidth);
    double r2Center=0.0625*(ringWidth)*(ringWidth);

    double valCenter=0.0;
    double valWings=0.0;
    double numCenter=0.0;
    double numWings=0.0;
    for (i=0;i<size;i++) {
      xy[1]=i-size/2-y0;
      for (j=0;j<size;j++) {
        xy[0]=j-size/2-x0;
        uv=matrix2x2_mul(wVectors,xy);
        r2=uv[0]*uv[0]+uv[1]*uv[1];
        if (r2<=r2WingsOuter) {
          d=pixels[i*size+j];
          if (r2<=r2Center){
            valCenter+=d*d;
            numCenter+=1.0;
          } else if (r2>r2WingsInner){
            valWings+=d*d;
            numWings+=1.0;
          }
        }
      }
    }
    if ((numWings==0.0) || (numCenter==0.0)) {
      System.out.println("Not enough data for correlation contrast: numCenter="+numCenter+" numWings="+numWings+
                         " valCenter="+IJ.d2s(valCenter,2)+" valWings="+IJ.d2s(valWings,2));
      return -1.0;
    }
    double contrast=Math.sqrt((valCenter/numCenter)/(valWings/numWings));
    if (DEBUG_LEVEL>2) {
      System.out.println("Correlation contrast is "+contrast);
      float [] floatPixels=new float[size*size];
      int index;
      for (i=0;i<size;i++) {
        xy[1]=i-size/2-y0;
        for (j=0;j<size;j++) {
          xy[0]=j-size/2-x0;
          uv=matrix2x2_mul(wVectors,xy);
          r2=uv[0]*uv[0]+uv[1]*uv[1];
          index=i*size+j;
/*          r=Math.sqrt(r2);
          r-=Math.floor(r);
          floatPixels[index]=(float) r;*/

          if (((r2<=r2WingsOuter) && (r2>r2WingsInner)) || (r2<=r2Center)){
            floatPixels[index]=(float) pixels[index];
          } else {
            floatPixels[index]=0.0F;
          }
        }
      }
      ImageProcessor ip = new FloatProcessor(size,size);
      ip.setPixels(floatPixels);
      ip.resetMinAndMax();
      ImagePlus imp= new ImagePlus(title+"_CORR_MASK", ip);
      imp.show();
    }
    return contrast;
  }

  public double[] correlateWithModel (double [] imagePixels,  // measured pixel array
                                      double [] modelPixels,  // simulated (model) pixel array)
                                              double sigma,   // Sigma for high pass filtering
                                              String title) { // title base for optional plots names

    if (imagePixels.length!=modelPixels.length) {
      IJ.showMessage("Error","Arrays have different sizes - imagePixels.length="+imagePixels.length+", modelPixels.length="+ modelPixels.length);
      return null;
    }
    int size = (int) Math.sqrt(imagePixels.length);
    ImageProcessor ip,ip_model;
    FHT fht,fht_model;
    double[][][] fft_complex,fft_model;
    int i,j;
    double a;

    float [] floatImagePixels=new float[size*size];
/* convert to float for image processor; */
    for (i=0;i<(size*size); i++) floatImagePixels[i]=(float) imagePixels[i];
    ip = new FloatProcessor(size,size);
    ip.setPixels(floatImagePixels);
    fht =  new FHT(ip);
// Swapping quadrants, so the center will be 0,0
    fht.swapQuadrants();
// get to frequency domain
    fht.transform();
    floatImagePixels=(float []) fht.getPixels();
  
    if((DEBUG_LEVEL>5) && (title!="")) {
      ImageProcessor ip_fht = new FloatProcessor(size,size);
      ip_fht.setPixels(floatImagePixels);
      ip_fht.resetMinAndMax();
      ImagePlus imp_fht= new ImagePlus(title+"_FHT_image", ip_fht);
      imp_fht.show();
    }
// Convert from FHT to complex FFT
    fft_complex= FHT2FFTHalf (fht,size);
    float [] floatModelPixels=new float[size*size];
// convert to float for image processor;
    for (i=0;i<(size*size); i++) floatModelPixels[i]=(float) modelPixels[i];
    ip_model = new FloatProcessor(size,size);
    ip_model.setPixels(floatModelPixels);
    fht_model =  new FHT(ip_model);
// Swapping quadrants, so the center will be 0,0
    fht_model.swapQuadrants();
// get to frequency domain
    fht_model.transform();
    floatModelPixels=(float []) fht_model.getPixels();
    if ((DEBUG_LEVEL>5) && (title!="")) {
      ImageProcessor ip_fht_model = new FloatProcessor(size,size);
      ip_fht_model.setPixels(floatModelPixels);
      ip_fht_model.resetMinAndMax();
      ImagePlus imp_fht_model= new ImagePlus(title+"_FHT_model", ip_fht_model);
      imp_fht_model.show();
    }
// Convert from FHT to complex FFT
    fft_model= FHT2FFTHalf (fht_model,size);

// multiply fft_complex by fft_nominator
      for (i=0;i<fft_complex.length; i++) for (j=0;j<fft_complex[0].length;j++) {
        a=                    fft_complex[i][j][0]*fft_model[i][j][0]+fft_complex[i][j][1]*fft_model[i][j][1]; // already changed Im() sign
        fft_complex[i][j][1]=-fft_complex[i][j][0]*fft_model[i][j][1]+fft_complex[i][j][1]*fft_model[i][j][0]; // already changed Im() sign
        fft_complex[i][j][0]=a;
      }   
/* Add sigma high=pass filtering here */
// Convert fft array back to fht array and
// set fht_target pixels with new values 
    fht.setPixels (floatFFTHalf2FHT (fft_complex,size));
/// optionally show the result
    if ((DEBUG_LEVEL>5) && (title!="")) {
      ImageProcessor ip_fht2 = new FloatProcessor(size,size);
      ip_fht2.setPixels(floatFFTHalf2FHT (fft_complex,size));
      ip_fht2.resetMinAndMax();
      ImagePlus imp_fht2= new ImagePlus(title+"-corr-sigma"+sigma, ip_fht2);
      imp_fht2.show();
    }

/// transform to space

   fht.inverseTransform();
   fht.swapQuadrants();
   fht.resetMinAndMax();

//   ImagePlus imp= new ImagePlus(title, ip_fht);
   if ((DEBUG_LEVEL>1) && (title!="")) {
     ImagePlus imp_corr= new ImagePlus(title+"_Correlated_filt-"+sigma, fht);
     imp_corr.show();
   }
//   return direct_target;
   floatImagePixels =(float[])fht.getPixels();
   double [] pixels=new double[floatImagePixels.length];
   for (i=0;i<floatImagePixels.length;i++) pixels[i]=floatImagePixels[i];
   return pixels;
  }


/* convert multi-file maps (both null/wave vectors and just boolean) into 2-d array of "number of files that includes this tile" */

  public int [][] mapToNumFiles (double[][][][][] maps,
                                      int square) {
    boolean [][][] bmaps=new boolean[maps.length][maps[0].length][maps[0][0].length];
    int i,j,n;
    for (n=0;n<maps.length;n++) for (i=0;i<maps[0].length;i++) for (j=0;j<maps[0][0].length;j++) bmaps[n][i][j]=(maps[n][i][j]!=null);
    return mapToNumFiles (bmaps,square);
  }

  public int [][] mapToNumFiles (boolean[][][] maps,
                                         int square) {//  size of the result tile (in units of maps)
     int mapSizeY=maps[0].length+   (square-1);
     int mapSizeX=maps[0][0].length+(square-1);
     boolean [][]bmap= new boolean [mapSizeY][mapSizeX];
     int [][] nMap=    new int     [mapSizeY][mapSizeX];
     int nTileY, nTileX, y, x,n;
     for (nTileY=0;nTileY<mapSizeY;nTileY++) for (nTileX=0;nTileX<mapSizeX;nTileX++) nMap[nTileY][nTileX]=0;
     for (n=0;n<maps.length;n++){
/* create a map for one file */     
       for (y=0;y<mapSizeY;y++) for (x=0;x<mapSizeX;x++) bmap[y][x]=false;
       for (nTileY=0;nTileY<maps[0].length;nTileY++) for (nTileX=0;nTileX<maps[0][0].length;nTileX++) if (maps[n][nTileY][nTileX]){
         for (y=0;y<square;y++) for (x=0;x<square;x++) bmap[nTileY+y][nTileX+x]=true;
       }
       for (y=0;y<mapSizeY;y++) for (x=0;x<mapSizeX;x++) if (bmap[y][x]) nMap[y][x]++;
     }
    return nMap;
  }
  /* build a map of available overlapping squares of diffrent size than initially mapped. New size and step should be multiple of initial step */
  public boolean [][] remapSquares (double [][][][] map, // [][]map of either null or 2 wave vectors
                                            int mapStep, // step of initial map
                                          int mapSquare, // size of square used in scanning of initial map (should be multiple of map step)
                                            int newStep, // step of the new map (should be multiple of map step)
                                          int newSquare){ // size of square used in sthe new map (should be multiple of map step)
     int nSteps=mapSquare/mapStep;
     int mapSizeY=map.length+   (nSteps-1);
     int mapSizeX=map[0].length+(nSteps-1);
     boolean [][]bmap= new boolean [mapSizeY][mapSizeX];
     int nTileY,nTileX,y,x;
     if (DEBUG_LEVEL>1) {
         System.out.println("remapSquares, map.length= "+map.length+" map[0].length="+map[0].length+" nSteps="+nSteps+" mapSizeY="+mapSizeY+" mapSizeX="+mapSizeX);
     }

/* create full map from initial */     
     for (y=0;y<mapSizeY;y++) for (x=0;x<mapSizeX;x++) bmap[y][x]=false;
     for (nTileY=0;nTileY<map.length;nTileY++) for (nTileX=0;nTileX<map[0].length;nTileX++) if (map[nTileY][nTileX]!=null){
       for (y=0;y<nSteps;y++) for (x=0;x<nSteps;x++) bmap[nTileY+y][nTileX+x]=true;
     }
/* build output map */
     int rsltStep=newStep/mapStep;
     int rsltSquare=newSquare/mapStep;

     int rsltSizeY=(mapSizeY-rsltSquare)/rsltStep+1; /* got zero when overlap was set to 16 - less than initial scan */
     int rsltSizeX=(mapSizeX-rsltSquare)/rsltStep+1;

     if (DEBUG_LEVEL>1) {
         System.out.println("remapSquares, rsltStep= "+rsltStep+" rsltSquare="+rsltSquare+" rsltSizeY="+rsltSizeY+" rsltSizeX="+rsltSizeX);
     }


     boolean [][]rslt= new boolean [rsltSizeY][rsltSizeX];
     for (nTileY=0;nTileY<rsltSizeY;nTileY++) for (nTileX=0;nTileX<rsltSizeX;nTileX++) {
       rslt[nTileY][nTileX]=true;
       for (y=0;y<rsltSquare;y++) for (x=0;x<rsltSquare;x++) if (!bmap[rsltStep*nTileY+y][rsltStep*nTileX+x]) rslt[nTileY][nTileX]=false;
     }
     return rslt;
  }


/********************************************************************/
  public double [][][][] scanImageForPatterns(ImagePlus imp,
                                                   int size, // FFT size
                                               double gamma, // gamma to use for power spectrum for correlation
                                               double sigma, // high-pass gaussian filter sigma when correlating power spectrum
                                       int diff_spectr_corr, // maximal distance between maximum on spectrum and predicted maximum on autocorrelation of gamma(|spectrum|)
                                     double shrink_clusters, // Shrink clusters by this ratio (remove lowest) after initial separation
                                       int multiples_to_try, // try this number of maximums proportionally farther from 0,0 than the two closest (increase precision)
                                           double deviation, // when looking for maximums - maximal distance from predicted from the lower order one
                                        int deviation_steps, // maximal iterations when looking for local maximum
                                        boolean [] bPattern,
                                            double highpass,
                                          double  ringWidth,
                                        double  minContrast,
                                            int debug_level) {  // debug level to use while iterating through steps
      if (imp==null){
        IJ.showMessage("Error","No image specified\nProcess canceled");
        return null;
      }
      Roi roi= imp.getRoi();
      Rectangle selection;
      if (imp.getType() !=ImagePlus.GRAY32 ) {
        if ((imp.getType() ==ImagePlus.GRAY8 ) ||
            (imp.getType() ==ImagePlus.GRAY16) ) {
           IJ.showStatus("Converting source image to gray 32 bits (float)");
           new ImageConverter(imp).convertToGray32();
        } else {
          IJ.showMessage("Error","Image should be Bayer array as a grayscale (8,16 or 32 bits)");
          return null;
        }
      }
      if (roi==null){
        selection=new Rectangle(0, 0, imp.getWidth(), imp.getHeight());
      } else {
        selection=roi.getBounds();
      }
      int mapWidth=  (imp.getWidth()/2-size)/(size/2)+1; // 50% overlap
      int mapHeight= (imp.getHeight()/2-size)/(size/2)+1; // 50% overlap
      String title=imp.getTitle();
      if (DEBUG_LEVEL>1) {
           System.out.println("Mapping image "+title+" with squares of "+(size*2)+"x"+(size*2)+"pixels with 50% overlap, covering total area of "+
                              ((mapWidth+1)*size)+"x"+((mapHeight+1)*size));
      }
      double [][][][] patternMap=new double [mapHeight][mapWidth][][];
      Rectangle mapCell=new Rectangle (0,0,size*2,size*2);
      int nTileX,nTileY,yc,xc;
      int wasDebug=DEBUG_LEVEL;
      double [][] pixels;
      double [] hamming=initHamming(size);
/* Prepare for filtering - common part outside of the iteration */
      double [][]  sim_pix;
      int this_simul_subdiv=2;
      double [] model_corr;
      double [][] convMatrix= {{1.0,-1.0},{1.0,1.0}}; // from greens2 to pixel WV
      double [][] invConvMatrix= matrix2x2_scale(matrix2x2_invert(convMatrix),2.0);
      double [][] WVgreens;
      double contrast;
      DEBUG_LEVEL=debug_level; // modify DEBUG_LEVEL to mute it while scanning many cells
      for (nTileY=0;nTileY<mapHeight;nTileY++) {
        if (DEBUG_LEVEL>0) IJ.showStatus("Mapping row "+(nTileY+1)+" (of "+mapHeight+")");
        for (nTileX=0;nTileX<mapHeight;nTileX++) patternMap[nTileY][nTileX]=null;
        yc=(nTileY+1)*size;
        if ((yc<selection.y) || (yc>=selection.y+selection.height)) continue;
        mapCell.y=nTileY*size;
        for (nTileX=0;nTileX<mapWidth;nTileX++) {
          xc=(nTileX+1)*size;
          if ((xc<selection.x) || (xc>=selection.x+selection.width)) continue;
          mapCell.x=nTileX*size;
          pixels=splitBayer(imp, mapCell,equalizeGreens);
          pixels[4]= normalizeAndWindow (pixels[4], hamming);
          patternMap[nTileY][nTileX]=findPattern(pixels[4],
                                                      size,
                                                     gamma,
                                                     sigma,
                                          diff_spectr_corr,
                                           shrink_clusters,
                                          multiples_to_try,
                                                 deviation,
                                           deviation_steps,
                                                      true,
                                                     title); // title - will not be used
/* Now verify by correlating with the actual pattern */
          if ((patternMap[nTileY][nTileX]!=null) && (minContrast>0.0)) {
            simulation_barray=  simulatePatternFullPattern(bPattern,
                                   patternMap[nTileY][nTileX][0][0],
                                   patternMap[nTileY][nTileX][0][1],
                                   patternMap[nTileY][nTileX][0][2],
                                   patternMap[nTileY][nTileX][1][0],
                                   patternMap[nTileY][nTileX][1][1],
                                   patternMap[nTileY][nTileX][1][2],
                                                               null, // no mesh distortion here
                                                  this_simul_subdiv, // simul_subdiv, - do not need high quality here
                                                               size,
                                                               true); // center fro greens
            sim_pix= extractSimulPatterns (simulation_barray,   // high resolution boolean pattern array
                                                  simul_fill, // part of the (center) pixel area being "phptosensitive"
                                           this_simul_subdiv,   // boolean pixels to real pixels resolution
                                                           1,
                                                        size,    // number of Bayer cells in width of the square selection (half number of pixels)
                                                         0.0,    // selection center, X (in pixels)
                                                         0.0);   // selection center, y (in pixels)
            sim_pix[4]= normalizeAndWindow (sim_pix[4], hamming);
            model_corr=correlateWithModel (pixels[4],  // measured pixel array
                                          sim_pix[4],  // simulated (model) pixel array)
                                                 0.0,  // double sigma,   // Sigma for high pass filtering
                                       imp.getTitle());
            WVgreens=matrix2x2_mul(patternMap[nTileY][nTileX],invConvMatrix);
            contrast= correlationContrast (model_corr,    // square pixel array
                                             WVgreens,    // wave vectors (same units as the pixels array)
                                            ringWidth,   // ring (around r=0.5 dist to opposite corr) width
                                                  0.0,    //  x0,              // center coordinates
                                                  0.0,    //y0,
                                               title);   // title base for optional plots names
//      System.out.println("Pattern correlation contrast= "+IJ.d2s(contrast,3)+ ", threshold is "+minCorrContrast);
            if (!(contrast >= minContrast)) patternMap[nTileY][nTileX]=null; // still getting NaN sometimes
          }
        }
      }
      DEBUG_LEVEL=wasDebug; // restore original debug level
//componentColorNames
//  public double[][] splitBayer (ImagePlus imp, Rectangle r) {
//imp_src.getWidth(), imp_src.getHeight(),imp_src.getTitle()
      return patternMap;
//
  }


// pixels - windowed (around center at [size/2, size/2]) measured pixels (one of the bayer components
// returns [2][3] array - 2 vectors of (fx,fy, phase)
// bug - changes pixels[]
// returns null on failure
  public double[][] findPattern(double [] input_pixels, // pixel array to process
                                              int size, // FFT size
                                          double gamma, // gamma to use for power spectrum for correlation
                                          double sigma, // high-pass gaussian filter sigma when correlating power spectrum
                                  int diff_spectr_corr, // maximal distance between maximum on spectrum and predicted maximum on autocorrelation of gamma(|spectrum|)
                                double shrink_clusters, // Shrink clusters by this ratio (remove lowest) after initial separation
                                  int multiples_to_try, // try this number of maximums proportionally farther from 0,0 than the two closest (increase precision)
                                      double deviation, // when looking for maximums - maximal distance from predicted from the lower order one
                                   int deviation_steps, // maximal iterations when looking for local maximum
                                        boolean greens, // this is a pattern for combined greens (diagonal), adjust results accordingly
                                          String title) { // title prefix to use for debug  images
    double [] pixels=input_pixels.clone();
    double [][]result=new double [2][3];
//System.out.println("pixels.length="+pixels.length); //4096

    ImageProcessor ip, ip1;
    FHT fht, fht1;
    double[][][] fft_complex,fft_corr;
    double[][] fft_gamma;
    int i,j;
    double DCLevel=0.0;
    double a;
    float []floatPixels=new float[pixels.length];
    for (i=0;i<pixels.length; i++)  DCLevel+=pixels[i];
    DCLevel/=(size*size);
    for (i=0;i<pixels.length; i++)  pixels[i]-=DCLevel;
// convert to float for image processor;
    for (i=0;i<pixels.length; i++) floatPixels[i]=(float) pixels[i];
    ip = new FloatProcessor(size,size);
    ip.setPixels(floatPixels);
    if (DEBUG_LEVEL>8) {
      ip.resetMinAndMax();
      ImagePlus imp_direct=  new ImagePlus(title+"_Direct_"+gamma, ip);
      imp_direct.show();
    }
    fht =  new FHT(ip);
// Swapping quadrants, so the center will be 0,0
    fht.swapQuadrants();
// get to frequency domain
    fht.transform();
    if (DEBUG_LEVEL>5) {
      floatPixels=(float []) fht.getPixels();
      ImageProcessor ip_fht = new FloatProcessor(size,size);
      ip_fht.setPixels(floatPixels);
      ip_fht.resetMinAndMax();
      ImagePlus imp_fht= new ImagePlus(title+"_FHT_"+deconvInvert, ip_fht);
      imp_fht.show();
    }

// Convert from FHT to complex FFT
    fft_complex= FHT2FFTHalf (fht,size);
// will need fft_complex  again later for later phase pattern measurements, calculate fft_gamma for correlation (pattern 2 frequencies measurement)
    fft_gamma=new double [size][size];
    floatPixels=new float[pixels.length];
    DCLevel=0.0;
    for (i=0;i<fft_complex.length; i++) for (j=0;j<fft_complex[0].length;j++) {
      fft_gamma[i][j]=Math.pow(fft_complex[i][j][0]*fft_complex[i][j][0]+fft_complex[i][j][1]*fft_complex[i][j][1],gamma);
      DCLevel+=fft_gamma[i][j];
      floatPixels[i*size+j]=(float) fft_gamma[i][j];
    }
    DCLevel/=(fft_complex.length*fft_complex[0].length);
    for (i=0;i<fft_complex.length; i++) for (j=0;j<fft_complex[0].length;j++) {
      floatPixels[i*size+j]-=DCLevel;
      if ((i>0)&& (i<(size/2))){
        floatPixels[(size-i)*size+((size-j)%size)]=floatPixels[i*size+j];
      }
    }

/* TODO:  maybe it is better to find the pattern frequencies just here, without converting back.
    After rejecting low frequencies, there seem to be just 2 nice maximums - easy to extract*/
// now perform direct FFT of gamma(power spectrum)  
    ip1 = new FloatProcessor(size,size);
    ip1.setPixels(floatPixels);
    if (DEBUG_LEVEL>7) {
      ip1.resetMinAndMax();
      ImagePlus imp1=  new ImagePlus(title+"_gamma(ps)_"+gamma, ip1);
      imp1.show();
    }
    fht1 =  new FHT(ip1);
// Swapping quadrants, so the center will be 0,0
    fht1.swapQuadrants();
    fht1.transform();
    fft_corr= FHT2FFTHalf (fht1,size);
    double[] highPassFilter=new double[fft_complex[0].length];
    double expK=(sigma>0)?(1.0/(2*sigma*sigma)):0.0;
    for (j=0;j<=fft_complex[0].length/2;j++) {
     highPassFilter[j]=(expK>0.0)?(1.0-Math.exp(-(expK*j*j))):1.0;
     if (j>0) highPassFilter[highPassFilter.length-j]=highPassFilter[j];
    }

    for (i=0;i<fft_complex.length; i++) for (j=0;j<fft_complex[0].length;j++) {
      fft_corr[i][j][0]=highPassFilter[i]*highPassFilter[j]*(fft_corr[i][j][0]*fft_corr[i][j][0]+fft_corr[i][j][1]*fft_corr[i][j][1]);
      fft_corr[i][j][1]=0.0;
    }

// Convert fft array back to fht array and
// set fht_target pixels with new values 
    fht1.setPixels (floatFFTHalf2FHT (fft_corr,size));   /* FIXME: - done, there is no difference as Im()==0 */
/// optionally show the result
    if (DEBUG_LEVEL>7) {
      ImageProcessor ip_fht2 = new FloatProcessor(size,size);
      ip_fht2.setPixels(floatFFTHalf2FHT (fft_corr,size));
      ip_fht2.resetMinAndMax();
      ImagePlus imp_fht2= new ImagePlus(title+"_fht_corr_"+gamma, ip_fht2);
      imp_fht2.show();
    }

/// transform to space
   fht1.inverseTransform();
   floatPixels=(float []) fht1.getPixels();
   a=1/floatPixels[0];
   for (i=0; i<floatPixels.length; i++){
      floatPixels[i]*=a;
   }
   fht1.setPixels(floatPixels);

//System.out.println("2:y="+y+" x="+x+" base_b="+base_b+" base="+base);

   fht1.swapQuadrants();

   if (DEBUG_LEVEL>2) {
     fht1.resetMinAndMax();
     ImagePlus imp_corr= new ImagePlus(title+"_corr_"+gamma, fht1);
     imp_corr.show();
   }
//   return direct_target;
   floatPixels =(float[])fht1.getPixels();
   for (i=0;i<floatPixels.length;i++) pixels[i]=floatPixels[i];


   int [][] max2OnSpectrum=  findFirst2MaxOnSpectrum (fft_complex, // complex, top half, starting from 0,0
                                                      1,   // skip +- from (0,0) and previous max - add parameter to dialog?
                                                      0.5); // 0.5 - 30deg. orthogonality of 2 vectors - 1.0 - perpendicular, 0.0 - parallel - add parameter to dialog?
/* Trying to filter out unreasonable maximums (if there is no pattern at all) */
   double maxFrequency=0.25*fft_complex.length;
   if ((Math.abs(max2OnSpectrum[0][0])>maxFrequency) ||
       (Math.abs(max2OnSpectrum[0][1])>maxFrequency) ||
       (Math.abs(max2OnSpectrum[1][0])>maxFrequency) ||
       (Math.abs(max2OnSpectrum[1][1])>maxFrequency)) {
       if (DEBUG_LEVEL>2) {
         System.out.println("Failed to detect pattern, as frequecy is above limit="+IJ.d2s(maxFrequency,2));
         System.out.println("Maximum 1 on spectrum:  x="+IJ.d2s(max2OnSpectrum[0][0],4)+" y="+IJ.d2s(max2OnSpectrum[0][1],4));
         System.out.println("Maximum 2 on spectrum:  x="+IJ.d2s(max2OnSpectrum[1][0],4)+" y="+IJ.d2s(max2OnSpectrum[1][1],4));
       }
       return null;
   }
   if (DEBUG_LEVEL>6) {
     System.out.println("Maximum 1 on spectrum:  x="+IJ.d2s(max2OnSpectrum[0][0],4)+" y="+IJ.d2s(max2OnSpectrum[0][1],4));
     System.out.println("Maximum 2 on spectrum:  x="+IJ.d2s(max2OnSpectrum[1][0],4)+" y="+IJ.d2s(max2OnSpectrum[1][1],4));
   }

   int [][] startPoints={{max2OnSpectrum[0][0]+max2OnSpectrum[1][0], max2OnSpectrum[0][1]+max2OnSpectrum[1][1]},
                         {max2OnSpectrum[0][0]-max2OnSpectrum[1][0], max2OnSpectrum[0][1]-max2OnSpectrum[1][1]}};
   if (startPoints[1][1] <0) { /* startPoints[1][1] > 0 anyway */
     startPoints[1][0]= -startPoints[1][0];
     startPoints[1][1]= -startPoints[1][1];
   }
   if (DEBUG_LEVEL>2) {
     System.out.println("Predicted correlation maximum 1 from spectrum:  x="+IJ.d2s(startPoints[0][0],4)+" y="+IJ.d2s(startPoints[0][1],4));
     System.out.println("Predicted correlation maximum 2 from spectrum:  x="+IJ.d2s(startPoints[1][0],4)+" y="+IJ.d2s(startPoints[1][1],4));
   }

   double[][] max2=  findFirst2MaxOnCorrelation(pixels, startPoints, diff_spectr_corr, shrink_clusters,  multiples_to_try, deviation, deviation_steps);

/**TODO:  get out on failure */
   if (max2==null) {
     if (DEBUG_LEVEL>1)  System.out.println("Failed to find a pattern");
     return null;
   }
/* these are combined greens, convert vectors to original pixel space) */
    if (greens) { 
      double [][] rotMatrix= {{1.0,-1.0},{1.0,1.0}};
      double [][] max2orig= matrix2x2_mul(max2,rotMatrix);
      for (i=0;i<2;i++) for (j=0;j<2;j++) result[i][j]=max2orig[i][j]; // result is [2][3], max2orig is [2][2]
      if (DEBUG_LEVEL>2) {
        System.out.println("Corrected to original pixels[0]  x="+IJ.d2s(result[0][0],4)+" y="+IJ.d2s(result[0][1],4));
        System.out.println("Corrected to original pixels[1]  x="+IJ.d2s(result[1][0],4)+" y="+IJ.d2s(result[1][1],4));
      }
    } else {
      for (i=0;i<2;i++) for (j=0;j<2;j++) result[i][j]=max2[i][j]; // result is [2][3], max2 is [2][2]
    }
/* Calculate locations of the maximums on FFT (corresponding to the diagonals of the checkerboard pattern) */
    double [][] maxOnFFT = {{2*size*(max2[0][0]-max2[1][0]),2*size*(max2[0][1]-max2[1][1])},
                            {2*size*(max2[0][0]+max2[1][0]),2*size*(max2[0][1]+max2[1][1])}};
/*  We have only one half of the FFT data  so rotate 180-degrees around the center if the point is in the bottom half*/
    double [] maxPhases=new double[2];
    boolean [] invertPhaseSign={false,false};
    int maxIndex; // iterate through the two maximums on FFT for phase measurement
    int [][]      interpolateXY= new int [2][2]; 
    double [][]   interpolateKxy=new double [2][2];
    double [][][] interpolatePhases=new double [2][2][2];
    double [][][] interpolateAmplitudes=new double [2][2][2];// Maybe use it? if the amplitudes are very different?
    int ix,iy;
    boolean phaseCorr; // phase shift before averaging, to prevent rollover
    for (maxIndex=0;maxIndex<2;maxIndex++) {
      if (maxOnFFT[maxIndex][1]<0) {
        invertPhaseSign[maxIndex]=true;
        maxOnFFT[maxIndex][0]=-maxOnFFT[maxIndex][0];
        maxOnFFT[maxIndex][1]=-maxOnFFT[maxIndex][1];
      }
      interpolateXY [maxIndex][1] = (int) maxOnFFT[maxIndex][1];
      interpolateKxy[maxIndex][1] = maxOnFFT[maxIndex][1] - interpolateXY [maxIndex][1];
      if (maxOnFFT[maxIndex][0]<0) {
        interpolateXY  [maxIndex][0] = (int) (size+ maxOnFFT[maxIndex][0]);
        interpolateKxy [maxIndex][0] = (size+ maxOnFFT[maxIndex][0]) - interpolateXY[maxIndex][0];
      } else {
        interpolateXY  [maxIndex][0] = (int) maxOnFFT[maxIndex][0];
        interpolateKxy [maxIndex][0] = maxOnFFT[maxIndex][0] - interpolateXY[maxIndex][0];
      }
      for (j=0;j<2;j++) {
        ix=(interpolateXY[maxIndex][0]+j) % size;
        for (i=0;i<2;i++) {
          iy=interpolateXY[maxIndex][1]+i;
          interpolateAmplitudes[maxIndex][i][j]= Math.sqrt(fft_complex[iy][ix][0]*fft_complex[iy][ix][0]+fft_complex[iy][ix][1]*fft_complex[iy][ix][1]);
          interpolatePhases[maxIndex][i][j]= Math.atan2((invertPhaseSign[maxIndex]?-1.0:1.0)*fft_complex[iy][ix][1], fft_complex[iy][ix][0]);
          if (DEBUG_LEVEL>5) {
            System.out.println("maxIndex="+maxIndex+" ix="+ix+" iy="+iy+" phase="+IJ.d2s(interpolatePhases[maxIndex][i][j],4)+" amplitude="+interpolateAmplitudes[maxIndex][i][j]);
          }
        }
      }
      phaseCorr=false;
      if ((interpolatePhases[maxIndex][0][0]> Math.PI/2) ||
          (interpolatePhases[maxIndex][0][1]> Math.PI/2) ||
          (interpolatePhases[maxIndex][1][0]> Math.PI/2) ||
          (interpolatePhases[maxIndex][1][1]> Math.PI/2) ||
          (interpolatePhases[maxIndex][0][0]<-Math.PI/2) ||
          (interpolatePhases[maxIndex][0][1]<-Math.PI/2) ||
          (interpolatePhases[maxIndex][1][0]<-Math.PI/2) ||
          (interpolatePhases[maxIndex][1][1]<-Math.PI/2)) {
          phaseCorr=true;
          interpolatePhases[maxIndex][0][0]+= (interpolatePhases[maxIndex][0][0]<0)?Math.PI:-Math.PI;
          interpolatePhases[maxIndex][0][1]+= (interpolatePhases[maxIndex][0][1]<0)?Math.PI:-Math.PI;
          interpolatePhases[maxIndex][1][0]+= (interpolatePhases[maxIndex][1][0]<0)?Math.PI:-Math.PI;
          interpolatePhases[maxIndex][1][1]+= (interpolatePhases[maxIndex][1][1]<0)?Math.PI:-Math.PI;
          if (DEBUG_LEVEL>5) {
            System.out.println("Shifting phases by PI/2 before averaging (to avoid rollover)");
          }
      }

      maxPhases[maxIndex]=       interpolateKxy[maxIndex][1] *(interpolateKxy[maxIndex][0] * interpolatePhases[maxIndex][1][1] + (1.0-interpolateKxy[maxIndex][0])* interpolatePhases[maxIndex][1][0])+
                          (1.0 - interpolateKxy[maxIndex][1])*(interpolateKxy[maxIndex][0] * interpolatePhases[maxIndex][0][1] + (1.0-interpolateKxy[maxIndex][0])* interpolatePhases[maxIndex][0][0]);
          if (phaseCorr) maxPhases[maxIndex]+=(maxPhases[maxIndex]<0)?Math.PI:-Math.PI;
      if (DEBUG_LEVEL>5) {
         System.out.println("kx="+IJ.d2s(interpolateKxy [maxIndex][0],4)+ " ky="+IJ.d2s(interpolateKxy [maxIndex][1],4));
      }
      if (DEBUG_LEVEL>2) {
        System.out.println("maxIndex="+maxIndex+" phase="+IJ.d2s(maxPhases[maxIndex],4));
      }
    }
    double [] checkerPhases= findCheckerPhases(max2, maxPhases); /* may be different for greens==true . No, the same */
    for (i=0;i<2;i++) result[i][2]=checkerPhases[i];
    if (DEBUG_LEVEL>2)  System.out.println();
   return result;
  }
/// returns -pixels[0].length/2<x<=pixels[0].length/2
///                          0<=y<=pixels[0].length/2
int [][] findFirst2MaxOnSpectrum (double [][][] pixels, // complex, top half, starting from 0,0
/* May need to reduce the skip_around to be able to handle smaller number of pattern periods? Or re-try if failed? Guess somehow?*/
                                  int skip_around, // skip +- from (0,0) and previous max 
                                  double minOrtho) { // 0.5 - 30deg. ortopgonality of 2 vectors - 1.0 - perpendicular, 0.0 - parallel
   int [][] max2={{0,0},{0,0}};
   double thisMax=0.0;
   int x,y,sx1,sx2;
   double p,a;
/* find first point */
   for (y=0;y<pixels.length;y++) for (x=0;x<pixels[0].length; x++ ) {
     p=pixels[y][x][0]*pixels[y][x][0]+pixels[y][x][1]*pixels[y][x][1];
     if (p>thisMax) {
       if ((y<=skip_around) && ((x<=skip_around) || (x>=pixels[0].length-skip_around))) continue; /* too close to [0,0] */
       max2[0][0]=x;
       max2[0][1]=y;
       thisMax=p;
     }
   }
   thisMax=0.0;
   sx1=(max2[0][0]>(pixels[0].length/2))?(max2[0][0]-pixels[0].length):max2[0][0]; /* y is always positive here */
/* find second point */
// Maybe also check if it is a local maximum (not on the border with protected area (0, first point, y=0, ...) ?
   for (y=0;y<pixels.length;y++) for (x=0;x<pixels[0].length; x++ ) {
     p=pixels[y][x][0]*pixels[y][x][0]+pixels[y][x][1]*pixels[y][x][1];
     if (p>thisMax) {
/* Is this a local maximum? */
       if (y>0) {
         if (                             p< (pixels[y-1][x  ][0]*pixels[y-1][x  ][0]+pixels[y-1][x  ][1]*pixels[y-1][x  ][1]))  continue;
         if ((x>0) &&                    (p< (pixels[y-1][x-1][0]*pixels[y-1][x-1][0]+pixels[y-1][x-1][1]*pixels[y-1][x-1][1]))) continue;
         if ((x<(pixels[0].length-1)) && (p< (pixels[y-1][x+1][0]*pixels[y-1][x+1][0]+pixels[y-1][x+1][1]*pixels[y-1][x+1][1]))) continue;
       }
       if (y< (pixels.length-1)) {
         if (                             p< (pixels[y+1][x  ][0]*pixels[y+1][x  ][0]+pixels[y+1][x  ][1]*pixels[y+1][x  ][1]))  continue;
         if ((x>0) &&                    (p< (pixels[y+1][x-1][0]*pixels[y+1][x-1][0]+pixels[y+1][x-1][1]*pixels[y+1][x-1][1]))) continue;
         if ((x<(pixels[0].length-1)) && (p< (pixels[y+1][x+1][0]*pixels[y+1][x+1][0]+pixels[y+1][x+1][1]*pixels[y+1][x+1][1]))) continue;
       }
       if ((x>0) &&                      (p< (pixels[y  ][x-1][0]*pixels[y  ][x-1][0]+pixels[y  ][x-1][1]*pixels[y  ][x-1][1]))) continue;
       if ((x<(pixels[0].length-1)) &&   (p< (pixels[y  ][x+1][0]*pixels[y  ][x+1][0]+pixels[y  ][x+1][1]*pixels[y  ][x+1][1]))) continue;
       if ((y<=skip_around) && ((x<=skip_around) || (x>=pixels[0].length-skip_around))) {
         if (DEBUG_LEVEL>5) {
           System.out.println("rejecting point ["+x+","+y+"] it is too close to [0,0]");
         }
          continue; /* too close to [0,0] */
       }

       if ((y<=skip_around) && ((x<=skip_around) || (x>=pixels[0].length-skip_around))) {
         if (DEBUG_LEVEL>5) {
           System.out.println("rejecting point ["+x+","+y+"] it is too close to [0,0]");
         }
          continue; /* too close to [0,0] */
       }
       if (((y<=(max2[0][1]+skip_around)) && (y>=(max2[0][1]-skip_around))) &&
           (((x<=(max2[0][0]+skip_around)) && (x>=(max2[0][0]-skip_around))) ||
             (x>=(max2[0][0]+pixels[0].length-skip_around)) || 
             (x<=(max2[0][0]-pixels[0].length+skip_around)))) {
         if (DEBUG_LEVEL>5) {
           System.out.println("rejecting point ["+x+","+y+"] as it is too close to the first one - ["+max2[0][0]+"("+sx1+"),"+max2[0][1]+"]");
         }
            continue; /* too close to first maximum */
       }
       sx2=(x>(pixels[0].length/2))?(x-pixels[0].length):x;
       a=(sx1*y -max2[0][1]*sx2);
       a=a*a/(sx1*sx1+max2[0][1]*max2[0][1])/(sx2*sx2+y*y);
       if (a < (minOrtho*minOrtho)) { /* vectors are too close to parallel */
         if (DEBUG_LEVEL>5) {
           System.out.println("rejecting point ["+x+"("+sx2+"),"+y+"] as the vector is too close to the first one - ["+max2[0][0]+"("+sx1+"),"+max2[0][1]+"]");
           System.out.println("pixels.length="+pixels.length+" pixels[0].length="+pixels[0].length);
         }
         continue;
       }
       max2[1][0]=x;
       max2[1][1]=y;
       thisMax=p;
     }
   }
   if (  thisMax==0.0) {
      System.out.println("Failed to find a second maximum");
      return null;
   }
   if (max2[0][0]>(pixels[0].length/2)) max2[0][0]-=pixels[0].length;
   if (max2[1][0]>(pixels[0].length/2)) max2[1][0]-=pixels[0].length;
   return max2;
}



/* Can it handle negative y if the refined maximum goes there? (maximal value on positive Y) */
  private double[][] findFirst2MaxOnCorrelation(double [] pixels,
                                            int [][] startPoints,
                                            int diff_spectr_corr,
                                          double shrink_clusters,
                                            int multiples_to_try,
                                                double deviation,
                                             int deviation_steps) {
    double reasonbleFrequency=2.0; // reject frequencies below that
    int size =(int) Math.sqrt (pixels.length);
    int [][] imax =startPoints.clone();
    int [][] imax2 =new int [2*multiples_to_try][2];
    boolean  []maxDefined=new boolean [2*multiples_to_try];

    double  [] maxValues =new double [startPoints.length];
    double  [] max2Values =new double [2];
    double [][] max2 =new double [2*multiples_to_try][2];
    int lim=size/2-2; // safety measure
    int nmax=0;
    int x=1;
    int y=0;
    int indx;
    int [] dirs = {-1, -size-1, -size, -size+1, 1, size+1, size, size-1};

    boolean isMax;
    int i,j,k, xmn,xmx,ymn,ymx;
    int halfSize=size/2;
    double [] vlengths={Math.sqrt(startPoints[0][0]*startPoints[0][0]+startPoints[0][1]*startPoints[0][1]),
                        Math.sqrt(startPoints[1][0]*startPoints[1][0]+startPoints[1][1]*startPoints[1][1])};
    boolean tooFar=false;
    int numVect;
    
/* Look for the maximal values around startPoints (+/-diff_spectr_corr )*/
    for (nmax=0; nmax<startPoints.length; nmax++) {
       ymn=imax[nmax][1]-diff_spectr_corr;
       ymx=imax[nmax][1]+diff_spectr_corr; 
       if (ymx>lim) ymx=lim;
       xmn=imax[nmax][0]-diff_spectr_corr;
       if (xmn<-lim) xmx=-lim;
       xmx=imax[nmax][0]+diff_spectr_corr;
       if (xmx>lim) xmx=lim;
       indx=(size+1)*size/2 + imax[nmax][1] *size+imax[nmax][0];
       if ((Math.abs(imax[nmax][0])>lim) || (Math.abs(imax[nmax][1])>lim)) {
         if (DEBUG_LEVEL>2) {
           System.out.println("Bad start point  imax[0][0]="+imax[0][0]+" imax[0][1]="+imax[0][1]+" imax[1][0]="+imax[1][0]+" imax[1][1]="+imax[1][1]+ " lim="+lim);
         }
         return null;
       }
//       if ((indx<0) || (indx>=pixels.length)) {
//          System.out.println(" imax[0][0]="+imax[0][0]+" imax[0][1]="+imax[0][1]+" imax[1][0]="+imax[1][0]+" imax[1][1]="+imax[1][1]+ " lim="+lim);
//       }
       maxValues[nmax]=pixels[indx];
       for (y=ymn;y<=ymx;y++) for (x=xmn;x<=xmx; x++) {
         indx=(size+1)*size/2 + y *size+x;
         if (pixels[indx]>maxValues[nmax]) {
/* Make sure it is closer to this point than to any other or [0,0] */
           tooFar=false;
           for (numVect=0;numVect<startPoints.length;numVect++) {
             if (Math.abs((x-imax[nmax][0])*imax[numVect][0]+(y-imax[nmax][1])*imax[numVect][1])>0.5*vlengths[nmax]*vlengths[numVect]) {
               tooFar=true;
               break;
             }
           }
           if (tooFar) {
             if (DEBUG_LEVEL>5) {
               System.out.println("rejecting point ["+x+","+y+"] as the vector is closer to other max than ["+imax[nmax][0]+","+imax[nmax][0]+"]"+
                                  " in the (+/-) direction: ["+imax[numVect][0]+","+imax[numVect][0]+"]");
             }
             continue;
           }
           maxValues[nmax]=pixels[indx];
           imax[nmax][0]=x;
           imax[nmax][1]=y;
         }
       }
/* Make sure the maximum in the scanned area is also a local maximum */
       isMax=true;
       indx=(size+1)*size/2 + imax[nmax][1] *size+imax[nmax][0];
       for (j=0;j<7;j++) if (pixels[indx]<pixels[indx+dirs[j]]) {
         isMax=false;
         break;
       }
       if (!isMax) {
         if (DEBUG_LEVEL>1) {
           System.out.println("This should not happen:");
           System.out.println("Maximum is not a local maximum - BUG or consider changing diff_spectr_corr="+diff_spectr_corr);
           System.out.println("point #"+(nmax+1)+" (of 2), x0="+startPoints[nmax][0]+" y0="+startPoints[nmax][1]+ " x="+imax[nmax][0]+" y="+imax[nmax][1]);
         }
/* Maybe return from here with null?*/
         while (!isMax && (indx>0) && (indx>pixels.length)) {
           isMax=true;
           for (j=0;j<7;j++) if (pixels[indx]<pixels[indx+dirs[j]]) {
             isMax=false;
             indx+=dirs[j];
             imax[nmax][1]=(indx / size) - halfSize;
             imax[nmax][1]=(indx % size) - halfSize;
            break;
           }
         }
         if (!isMax) {
           if (DEBUG_LEVEL>1) {
             System.out.println("Maximum still not reached, bailing out");
             System.out.println("point #"+(nmax+1)+" (of 2), x0="+startPoints[nmax][0]+" y0="+startPoints[nmax][1]+ " x="+imax[nmax][0]+" y="+imax[nmax][1]);
           }
           return null;
         } else {
           if (DEBUG_LEVEL>2) {
             System.out.println("point #"+(nmax+1)+" (of 2), corrected local maximum is  x="+imax[nmax][0]+" y="+imax[nmax][1]);
           }
         }

       }
    }
/* Sort maximums so first vector to second vector will be clockwise (positive y is downwards) */
    j=0; k=1;
    if ((imax[j][0]*imax[k][1]-imax[k][0]*imax[j][1])<0) {
      j=1;k=0;
    }
    imax2[0][0]=imax[j][0];
    imax2[0][1]=imax[j][1];
    imax2[1][0]=imax[k][0];
    imax2[1][1]=imax[k][1];


/* Now define maximal radius of cluster (~0.7 of the average distance from 0,0 to the 2 start points */
    int maxX2Y2=0;
    for (i=0;i<2;i++) for (j=0;j<2;j++) maxX2Y2+=imax2[i][j]*imax2[i][j];
    maxX2Y2/=4;
//System.out.println("maxX2Y2="+maxX2Y2);

    for (i=0;i<2*multiples_to_try;i++)  maxDefined[i]=(i<2);

    nmax=2*multiples_to_try; /* but only the first two are known by now */
    for (i=0;i<2;i++) {
if (DEBUG_LEVEL>5) System.out.println("i="+i+" x="+imax2[i][0]+" y="+imax2[i][1]+" value="+max2Values[i]);
    }

    int clusterNumber;
    List <Integer> pixelList=new ArrayList<Integer>(100);
    Integer Index, NewIndex, HillIndex;
    double cx,cy,cm,minInCluster,f;
    
    int []clusterMap=new int[pixels.length];
    for (i=0;i<clusterMap.length;i++) clusterMap[i]=0; /// 0 - unused, -1 - "do not use"
    int listIndex;
    boolean noHills;
    int clusterSize;
    boolean isLocalMax;
    int pair;
    for (clusterNumber=0;clusterNumber<nmax;clusterNumber++) {
      pair=clusterNumber/2+1;
      if (!maxDefined[clusterNumber] && maxDefined[clusterNumber-2]) {
/* We do not know the seed for this maximum, but the previous (of the same direction) may be known */
        x= (int) (max2[clusterNumber-2][0]*pair)/(pair-1);
        y= (int) (max2[clusterNumber-2][1]*pair)/(pair-1);
        if ((x>(-lim+1)) && (x<(lim-1)) && (y>(-lim+1)) && (y<(lim-1))) {
          Index=(y+halfSize)*size + x+halfSize;
/* there should be local maximum not more than "deviation_steps" steps from the x,y */
          isLocalMax=false;
//          i=(int) (2*deviation+1);
          i=deviation_steps;
          while ((i>0) && !isLocalMax) {
            isLocalMax=true;
            for (j=0;j<dirs.length;j++) {
              NewIndex=Index+dirs[j];
              if ((NewIndex>=0) && (NewIndex<clusterMap.length) && (pixels[NewIndex]>pixels[Index])) {
                isLocalMax=false;
                Index=NewIndex;
                i--;
if (DEBUG_LEVEL>5) System.out.println("i="+i+" x="+((Index % size) - halfSize)+" y="+((Index / size) - halfSize)+" value="+pixels[Index]);
                break;
              }
            }
          }
          if (isLocalMax && (clusterMap[Index]==0)) { // not yet used
             imax2[clusterNumber][0]= (Index % size) - halfSize;
             imax2[clusterNumber][1]= (Index / size) - halfSize;
             maxDefined[clusterNumber]=true;   
          }

        }
      }
      if (maxDefined[clusterNumber]) { // skip if seed for the cluster is not defined
/* Grow cluster around maximum, find centroid */
        Index=(imax2[clusterNumber][1]+halfSize)*size + imax2[clusterNumber][0]+halfSize;
        pixelList.clear();
        pixelList.add (Index);
        clusterMap[Index]=clusterNumber+1;
        listIndex=0;
        while (listIndex<pixelList.size() ) {
          Index=pixelList.get(listIndex++);
          for (j=0;j<dirs.length;j++) {
            NewIndex=Index+dirs[j];
            if ((NewIndex>=0) && (NewIndex<clusterMap.length) && (clusterMap[NewIndex]==0) && (pixels[NewIndex]<pixels[Index])) {
/* did we get too far?*/           
             y=(NewIndex/size) - halfSize - imax2[clusterNumber][1];
             x=(NewIndex % size) - halfSize - imax2[clusterNumber][0];
//System.out.println(" dy="+y+" dx="+x+" dx*dx+dy*dy="+(x*x+y*y));
             if ((x*x+y*y) <= maxX2Y2) {
/* See if there is any neighbor of the new pixel that is higher and not yet marked (prevent rivers flowing between hills) */
                noHills=true;
                for (k=0;k<dirs.length;k++) {
                  HillIndex=NewIndex+dirs[k];
                  if ((HillIndex>=0) && (HillIndex<clusterMap.length) && (clusterMap[HillIndex]!=(clusterNumber+1)) && (pixels[HillIndex]>pixels[NewIndex])) {
                    noHills=false;
                    break;
                  }
                }
                if (noHills) {
                  pixelList.add (NewIndex);
                  clusterMap[NewIndex]=clusterNumber+1;
//System.out.println("NewIndex="+NewIndex+" y="+(NewIndex/size - halfSize)+" x="+((NewIndex % size) - halfSize)+" new pixel="+pixels[NewIndex]+" old pixel="+pixels[Index]);
                }
              }
            }
          }
        }
/* Shrink clusters to a fraction of initial size */
// TODO: shring to a value (between min and max) if there is a sharp maximum??
//, double shrink_clusters
        if (shrink_clusters==0.0) { // use "smart" size
           clusterSize=(int) Math.sqrt(5* pixelList.size()); // use proportional size
        } else if (shrink_clusters<0) {
           clusterSize=(int)(- shrink_clusters ); // use specified size
        } else {
          clusterSize=(int) (pixelList.size()*shrink_clusters); // use proportional size
        }
        if (clusterSize<5) clusterSize=5;
        while (pixelList.size()>clusterSize) {
          i=0;
          f=pixels[pixelList.get(i)];
          for (j=1;j<pixelList.size();j++) if (pixels[pixelList.get(j)]<f){
           i=j;
           f=pixels[pixelList.get(j)];
          }
          clusterMap[pixelList.get(i)]=-1; // Do not use looking for the next cluster
          pixelList.remove(i);
        }



/* now find centroid of the cluster */
        minInCluster=pixels[pixelList.get(0)];
        for (i=1;i<pixelList.size();i++) if (minInCluster>pixels[pixelList.get(i)]) minInCluster = pixels[pixelList.get(i)];
        cx=0.0; cy=0.0; cm=0.0;
        for (i=0;i<pixelList.size();i++) {
          j=pixelList.get(i);
          y=j / size - halfSize;
          x=j % size - halfSize;
          f=pixels[j]-minInCluster;
          cm+=f;
          cx+=f*x;
          cy+=f*y;
        }
        cx/=cm;
        cy/=cm;
        max2[clusterNumber][0]=cx;
        max2[clusterNumber][1]=cy;
        f=0.0;
        if (pair>1) {
          cx=max2[clusterNumber-2][0]*pair/(pair-1)-max2[clusterNumber][0];
          cy=max2[clusterNumber-2][1]*pair/(pair-1)-max2[clusterNumber][1];
          f=Math.sqrt(cx*cx+cy*cy);
/* Verify deviation here */
          if (f>deviation) maxDefined[clusterNumber]=false;
        }
if (DEBUG_LEVEL>6) System.out.println("pixelList.size()="+pixelList.size()+" centroid sum="+cm);
if (DEBUG_LEVEL>5) System.out.println("clusterNumber="+clusterNumber+" x="+max2[clusterNumber][0]+" y="+max2[clusterNumber][1] + " x0="+(max2[clusterNumber][0]/pair)+" y0="+(max2[clusterNumber][1]/pair)+" deviat="+f);
        if ((cm==0.0) || (pixelList.size()<3)) maxDefined[clusterNumber]=false;
/* Filter out unreasonably low frequencies*/
        if ((max2[clusterNumber][0]*max2[clusterNumber][0]+max2[clusterNumber][1]*max2[clusterNumber][1])<(reasonbleFrequency*reasonbleFrequency)) {
          if (DEBUG_LEVEL>2) System.out.println("Frequency too low:clusterNumber="+clusterNumber+" x="+max2[clusterNumber][0]+" y="+max2[clusterNumber][1]+ ", minimal allowed frequency is "+reasonbleFrequency);
          maxDefined[clusterNumber]=false;
        }
      }
    }
/* Average (or just use farthest?) multiple maximums */
    if (DEBUG_LEVEL>2){
      float [] dbg_pixels=new float[clusterMap.length];
      for (j=0;j<dbg_pixels.length;j++) dbg_pixels[j]=(float)clusterMap[j];
      dbg_pixels[(size+1)*size/2]=-1; // mark center

      ImageProcessor ip=new FloatProcessor(size,size);
      ip.setPixels(dbg_pixels);
      ip.resetMinAndMax();
      ImagePlus imp=  new ImagePlus("clusters", ip);
      imp.show();
    }

//    for (i=0;i<2;i++) for (j=0;j<2;j++) max2[i][j]=imax2[i][j];
    double [][]maxFinal=new double[2][2];
    boolean [] definedFinal={false,false};
    for (i=0;i<2;i++) {
      maxFinal[i][0]=0.0;
      maxFinal[i][1]=0.0;
      j=0;
      for (pair=0;pair<nmax/2;pair++) {
        if (maxDefined[i+2*pair]) {
          j++;
          maxFinal[i][0]+=max2[i+2*pair][0]/(pair+1);
          maxFinal[i][1]+=max2[i+2*pair][1]/(pair+1);
          definedFinal[i]=true;
//          System.out.println("i="+i+" pair="+pair+" j="+j+" maxFinal["+i+"][0]="+maxFinal[i][0]+" maxFinal["+i+"][1]="+maxFinal[i][1]);
        }
      }
      maxFinal[i][0]/=j;
      maxFinal[i][1]/=j;
      if (j==0) definedFinal[i]=false;
/* These two vectors correspond to checker diagonals and they are calculated for the correlation space.
 Actual frequencies for the checker board will be 1/4 of these, also divide by FFT size so the result will be in cycles per pixel */
      maxFinal[i][0]/=size*4;
      maxFinal[i][1]/=size*4;
    }
    if (DEBUG_LEVEL>2) {
      System.out.println("Checkerboard frequency[0]  x="+IJ.d2s(maxFinal[0][0],4)+" y="+IJ.d2s(maxFinal[0][1],4));
      System.out.println("Checkerboard frequency[1]  x="+IJ.d2s(maxFinal[1][0],4)+" y="+IJ.d2s(maxFinal[1][1],4));
//      System.out.println();
    }
    if (!definedFinal[0] || !definedFinal[1]) {
      if (DEBUG_LEVEL>2) {
        System.out.println("Undefined frequency(ies)");
      }
      return null;
    }
    return maxFinal;
  }


/**

f1,f2 - vectors of the checker board
g1,g2 - wave vectors going diagonally through the centers of whites
g1=f1-f2, g2=f1+f2
p1,p2 - their phases
calculating location of the white center as intesection of the two "waves"
l1, l2 - vectors connecting two wave fronts perpendicular
W = location of the white center
l1=  g1/abs(g1)^2
l2=  g2/abs(g2)^2

R90 = | 0 1|
      |-1 0| - rotation matrix

(1)W= -l1 *p1/(2*pi) + l1*R90*u = -l2 *p2/(2*pi) + l2*R90*v

(3) -l1X*p1/(2*pi) +l1Y*u = -l2X*p2/(2*pi) +l2Y*v
(4) -l1Y*p1/(2*pi) -l1X*u = -l2Y*p2/(2*pi) -l2X*v

  l1Y*u  - l2Y*v + (l2X*p2-l1X+p1)/(2*pi) =0
 -l1X*u  + l2X*v + (l2Y*p2-l1Y*p1)/(2*pi) =0

l1X*  l1Y*u  - l1X*l2Y*v + l1X*(l2X*p2-l1X*p1)/(2*pi) =0
l1Y* -l1X*u  + l1Y*l2X*v + l1Y*(l2Y*p2-l1Y*p1)/(2*pi) =0

l1X*  l1Y*u  - l1X*l2Y*v + l1X*(l2X*p2-l1X*p1)/(2*pi) + l1Y* -l1X*u  + l1Y*l2X*v + l1Y*(l2Y*p2-l1Y*p1)/(2*pi) =0
             - l1X*l2Y*v + l1X*(l2X*p2-l1X*p1)/(2*pi) +                l1Y*l2X*v + l1Y*(l2Y*p2-l1Y*p1)/(2*pi) =0
            (l1Y*l2X - l1X*l2Y)*v    + l1X*(l2X*p2-l1X*p1)/(2*pi) +  l1Y*(l2Y*p2-l1Y*p1)/(2*pi) =0
            v= (l1X*(l2X*p2-l1X+p1)/(2*pi) +  l1Y*(l2Y*p2-l1Y*p1)/(2*pi))/(l1X*l2Y-l1Y*l2X)

            v= (l1X*(l2X*p2-l1X*p1) +  l1Y*(l2Y*p2-l1Y*p1))/(l1X*l2Y-l1Y*l2X)/(2*pi)

W=  -l2 *p2/(2*pi) + l2*R90*v
Wx= -l2X *p2/(2*pi) + l2Y*v
Wy= -l2Y *p2/(2*pi) - l2X*v

Wx= -l2X *p2/(2*pi) + l2Y* (l1X*(l2X*p2-l1X*p1) +  l1Y*(l2Y*p2-l1Y*p1))/(l1X*l2Y-l1Y*l2X)/(2*pi)
Wy= -l2Y *p2/(2*pi) - l2X* (l1X*(l2X*p2-l1X*p1) +  l1Y*(l2Y*p2-l1Y*p1))/(l1X*l2Y-l1Y*l2X)/(2*pi)

Wx= (-l2X *p2 + l2Y* (l1X*(l2X*p2-l1X*p1) +  l1Y*(l2Y*p2-l1Y*p1))/(l1X*l2Y-l1Y*l2X))/(2*pi)
Wy= (-l2Y *p2 - l2X* (l1X*(l2X*p2-l1X*p1) +  l1Y*(l2Y*p2-l1Y*p1))/(l1X*l2Y-l1Y*l2X))/(2*pi)

Next step - find phases of f1, f2

distance from W to the wave front f1 going through (0,0) is a scalar product (*)
 (-W) * f1/abs(f1)
linear period of f1 is 1/abs(f1)
phase of wave f1 will be 
 2*pi* ((-W) * f1/abs(f1))/ (1/abs(f1))
=
 2*pi* ((-W)* f1)

so
 phase(f1)= -2*pi*(Wx*f1x+Wy*f1y)
 phase(f2)= -2*pi*(Wx*f2x+Wy*f2y)

*/

  private double [] findCheckerPhases(double [][] WVectors, double [] P) {
    double [][] DWVectors = {{WVectors[0][0]-WVectors[1][0],WVectors[0][1]-WVectors[1][1]},
                             {WVectors[0][0]+WVectors[1][0],WVectors[0][1]+WVectors[1][1]}};
    if (DEBUG_LEVEL>3)  System.out.println("      DWVectors[0][0]="+IJ.d2s(DWVectors[0][0],4)+"  DWVectors[0][1]="+IJ.d2s(DWVectors[0][1],4));
    if (DEBUG_LEVEL>3)  System.out.println("      DWVectors[1][0]="+IJ.d2s(DWVectors[1][0],4)+"  DWVectors[1][1]="+IJ.d2s(DWVectors[1][1],4));
    double [] DWVectorsAbs2={DWVectors[0][0]*DWVectors[0][0]+DWVectors[0][1]*DWVectors[0][1],
                             DWVectors[1][0]*DWVectors[1][0]+DWVectors[1][1]*DWVectors[1][1]};
    if (DEBUG_LEVEL>3)  System.out.println("      sqrt(DWVectorsAbs2[0])="+IJ.d2s(Math.sqrt(DWVectorsAbs2[0]),4)+"  sqrt(DWVectorsAbs2[1])="+IJ.d2s(Math.sqrt(DWVectorsAbs2[1]),4));
    double [][] DL= {{DWVectors[0][0]/DWVectorsAbs2[0],DWVectors[0][1]/DWVectorsAbs2[0]},
                     {DWVectors[1][0]/DWVectorsAbs2[1],DWVectors[1][1]/DWVectorsAbs2[1]}};
    if (DEBUG_LEVEL>3)  System.out.println("      DL[0][0]="+IJ.d2s(DL[0][0],4)+"  DL[0][1]="+IJ.d2s(DL[0][1],4));
    if (DEBUG_LEVEL>3)  System.out.println("      DL[1][0]="+IJ.d2s(DL[1][0],4)+"  DL[1][1]="+IJ.d2s(DL[1][1],4));

//         v= (l1X     *(l2X     *p2  -l1X     *p1  ) +  l1Y     *(l2Y     *p2  -l1Y     *p1  ))/(l1X     *l2Y     -l1Y     *l2X     )/(2*pi)
    double v= (DL[0][0]*(DL[1][0]*P[1]-DL[0][0]*P[0]) +  DL[0][1]*(DL[1][1]*P[1]-DL[0][1]*P[0]))/(DL[0][0]*DL[1][1]-DL[0][1]*DL[1][0])/(2*Math.PI);
    if (DEBUG_LEVEL>2)  System.out.println("v="+IJ.d2s(v,4));


//            Wx= -l2X      *p2    / (2*pi)      + l2Y      * v
//            Wy= -l2Y      *p2    / (2*pi)      - l2X      * v
    double [] WC={-DL[1][0] * P[1] / (2*Math.PI) + DL[1][1] * v,
                  -DL[1][1] * P[1] / (2*Math.PI) - DL[1][0] * v};
    if (DEBUG_LEVEL>2)  System.out.println("WC[0]="+IJ.d2s(WC[0],4)+"  WC[1]="+IJ.d2s(WC[1],4));

// phase(f1)=         -2*pi      *(Wx    * f1x            + Wy    * f1y           )
// phase(f2)=         -2*pi      *(Wx    * f2x            + Wy    * f2y           )

    double [] phases={-2*Math.PI *(WC[0] * WVectors[0][0] + WC[1] * WVectors[0][1]),
                      -2*Math.PI *(WC[0] * WVectors[1][0] + WC[1] * WVectors[1][1])};
    if (DEBUG_LEVEL>2)  System.out.println("phases[0]="+IJ.d2s(phases[0],4)+"  phases[1]="+IJ.d2s(phases[1],4));
    return phases;
  }


  public boolean showStackConvolutionDialog() {
    int i;
    GenericDialog gd = new GenericDialog("Stack convolution parameters");
    gd.addNumericField("Convolution FFT size (twice the kernel size)",                             CONVOLVE_FFT_SIZE,     0); // 128
    gd.addCheckbox    ("Update ImageJ status",                                                         UPDATE_STATUS);

    gd.addNumericField("Debug Level:",                          MASTER_DEBUG_LEVEL,      0);
    gd.showDialog();
    if (gd.wasCanceled()) return false;
    CONVOLVE_FFT_SIZE=1;
    for (i=(int) gd.getNextNumber(); i >1; i>>=1) CONVOLVE_FFT_SIZE <<=1; /* make it to be power of 2 */
    UPDATE_STATUS=               gd.getNextBoolean();
    MASTER_DEBUG_LEVEL=    (int) gd.getNextNumber();
    return true;
 }  



  public boolean showInterpolateKernelsDialog() {
    GenericDialog gd = new GenericDialog("Interpolate kernels parameters");
    gd.addNumericField("Input kernel size",                                                            INTERPOLATE_INSIZE,     0); // 64
    gd.addNumericField("Interpolation step between original kernels",                                  INTERPOLATE_STEP,       0); //4
    gd.addNumericField("Interpolation add on the top (in output, subdivided steps)",                   INTERPOLATE_ADDTOP,     0); //16
    gd.addNumericField("Interpolation add on the left (in output, subdivided steps)",                  INTERPOLATE_ADDLEFT,    0); //16
    gd.addNumericField("Interpolation add on the right (in output, subdivided steps)",                 INTERPOLATE_ADDRIGHT,   0); //16
    gd.addNumericField("Interpolation add on the bottom (in output, subdivided steps)",                INTERPOLATE_ADDBOTTOM,  0); //16
    gd.addNumericField("Interpolation: extrapolate margins - 0.0 - duplicate, 1.0 - full extrapolate", INTERPOLATE_EXTRAPOLATE,3); //1.0
    gd.addCheckbox    ("Update ImageJ status",                                                         UPDATE_STATUS);

    gd.addNumericField("Debug Level:",                          MASTER_DEBUG_LEVEL,      0);
    gd.showDialog();
    if (gd.wasCanceled()) return false;
    INTERPOLATE_INSIZE=    (int) gd.getNextNumber();
    INTERPOLATE_STEP=      (int) gd.getNextNumber();
    INTERPOLATE_ADDTOP=    (int) gd.getNextNumber();
    INTERPOLATE_ADDLEFT=   (int) gd.getNextNumber();
    INTERPOLATE_ADDRIGHT=  (int) gd.getNextNumber();
    INTERPOLATE_ADDBOTTOM= (int) gd.getNextNumber();
    INTERPOLATE_EXTRAPOLATE=     gd.getNextNumber();
    UPDATE_STATUS=               gd.getNextBoolean();
    MASTER_DEBUG_LEVEL=    (int) gd.getNextNumber();
    return true;
 }  

  public boolean showSplitBayerToStackDialog() {
    GenericDialog gd = new GenericDialog("Interpolate kernels parameters");
    gd.addNumericField("Interpolation step between original kernels",                                  SPLIT_OVERSAMPLE, 0); //2
    gd.addNumericField("Interpolation add on the top (in output, subdivided steps)",                   SPLIT_ADDTOP,     0); //32
    gd.addNumericField("Interpolation add on the left (in output, subdivided steps)",                  SPLIT_ADDLEFT,    0); //32
    gd.addNumericField("Interpolation add on the right (in output, subdivided steps)",                 SPLIT_ADDRIGHT,   0); //32
    gd.addNumericField("Interpolation add on the bottom (in output, subdivided steps)",                SPLIT_ADDBOTTOM,  0); //32

    gd.addNumericField("Debug Level:",                          MASTER_DEBUG_LEVEL,      0);
    gd.showDialog();
    if (gd.wasCanceled()) return false;
    SPLIT_OVERSAMPLE=   (int) gd.getNextNumber();
    SPLIT_ADDTOP=       (int) gd.getNextNumber();
    SPLIT_ADDLEFT=      (int) gd.getNextNumber();
    SPLIT_ADDRIGHT=     (int) gd.getNextNumber();
    SPLIT_ADDBOTTOM=    (int) gd.getNextNumber();
    MASTER_DEBUG_LEVEL= (int) gd.getNextNumber();
    return true;
 }  

  public boolean showGaussianStackDialog() {
    int i;
    GenericDialog gd = new GenericDialog("Gaussian stack generation parameters");
    gd.addNumericField("Direct kernel size:",                                            INVERSE_DIRECT_SIZE, 0); // 32
    gd.addNumericField("Direct kernel size:",                                            INVERSE_REVERSE_SIZE, 0); //64

    gd.addNumericField("PSF cutoff energy (used to determine bounding ellipse)",          PSF_cutoff_energy, 3); //0.9;  // Limit result kernel to proportional of the PSF, calculate initial cluster shape by this cutoff energy

    gd.addNumericField("Gaussian blur for individual colors",                             DECONCV_BLUR_INDIVIDUAL, 3); //1.8
    gd.addNumericField("Gaussian blur for checkerboard greens",                           DECONCV_BLUR_CHECKER,    3); //1.4

    gd.addCheckbox    ("Update ImageJ status",                                            UPDATE_STATUS);
    gd.addNumericField("Debug Level:",                                                    MASTER_DEBUG_LEVEL, 0);
    gd.showDialog();
    if (gd.wasCanceled()) return false;

    INVERSE_DIRECT_SIZE=1;
    for (i=(int) gd.getNextNumber(); i >1; i>>=1) INVERSE_DIRECT_SIZE <<=1; /* make it to be power of 2 */
    INVERSE_REVERSE_SIZE=1;
    for (i=(int) gd.getNextNumber(); i >1; i>>=1) INVERSE_REVERSE_SIZE <<=1; /* make it to be power of 2 */
    PSF_cutoff_energy=             gd.getNextNumber();
    DECONCV_BLUR_INDIVIDUAL=       gd.getNextNumber();
    DECONCV_BLUR_CHECKER=          gd.getNextNumber();
    UPDATE_STATUS=                 gd.getNextBoolean();
    MASTER_DEBUG_LEVEL=      (int) gd.getNextNumber();
    return true;
  }




  public boolean showInverseStackDialog() {
    int i;
    GenericDialog gd = new GenericDialog("PSF stack inversion parameters");
    gd.addNumericField("Direct kernel size:",                                            INVERSE_DIRECT_SIZE, 0); // 32
    gd.addNumericField("Direct kernel size:",                                            INVERSE_REVERSE_SIZE, 0); //64

    gd.addNumericField("OTF cutoff energy (used to determine bounding ellipse)",          OTF_cutoff_energy, 3); //0.6; use frequency points that have OTF_cutoff_energy of the total to determine ellipse for limiting frequency responce
    gd.addNumericField("OTF size of elliptical window relative to cluster size",          OTF_ellipse_scale, 3); //1.5;  // size of elliptical window relative to the cluster defined by OTF_cutoff_energy
    gd.addCheckbox    ("Use Gauss window (instead of Hamming) for Ellipse mask",          OTF_ellipse_gauss); //   // size of elliptical window relative to the cluster defined by OTF_cutoff_energy
    gd.addNumericField("OTF deconvolution parameter ",                                    OTF_deconvInvert, 3); //0.01; // when FFT component is less than this fraction of the maximal value, replace 1/z with Z
    gd.addNumericField("PSF cutoff energy (used to determine bounding ellipse)",          PSF_cutoff_energy, 3); //0.9;  // Limit result kernel to proportional of the PSF, calculate initial cluster shape by this cutoff energy
    gd.addNumericField("Reversed PSF elliptical mask relative to the direct PSF ellipse", PSF_ellipse_scale, 3); //1.5;  // size of elliptical window to limuit reverse PSF as proportional to direct one
    gd.addNumericField("Reversed PSF mask thershold (completely ignore far pixels)",      RPSF_min_mask_threshold, 3); //0.01; // completely zero reversed kernel elements where elliptical mask is below this threshold

    gd.addCheckbox    ("Filter inverted PSF",                                             INVERSE_FILTER);

    gd.addNumericField("Reversed PSF sigma to radius ratio (0 to disable)",               RPSF_sigma_to_radius, 3); //0.2- variable blurring - sigma will be proportional distance from the center
    gd.addNumericField("Reversed PSF: scale variable sigma (in the center) from the uniform one",  RPSF_var_sigma_scale, 3); //=0.8; // reduce variable sigma in the center from uniuform one

    gd.addNumericField("Gaussian blur for individual colors",                             DECONCV_BLUR_INDIVIDUAL, 3); //1.8
    gd.addNumericField("Gaussian blur for checkerboard greens",                           DECONCV_BLUR_CHECKER,    3); //1.4


    gd.addCheckbox    ("Update ImageJ status",                                            UPDATE_STATUS);
    gd.addNumericField("Debug Level:",                                                    MASTER_DEBUG_LEVEL, 0);
    gd.showDialog();
    if (gd.wasCanceled()) return false;

    INVERSE_DIRECT_SIZE=1;
    for (i=(int) gd.getNextNumber(); i >1; i>>=1) INVERSE_DIRECT_SIZE <<=1; /* make it to be power of 2 */
    INVERSE_REVERSE_SIZE=1;
    for (i=(int) gd.getNextNumber(); i >1; i>>=1) INVERSE_REVERSE_SIZE <<=1; /* make it to be power of 2 */


    OTF_cutoff_energy=             gd.getNextNumber();
    OTF_ellipse_scale=             gd.getNextNumber();
    OTF_ellipse_gauss=             gd.getNextBoolean();
    OTF_deconvInvert=              gd.getNextNumber();
    PSF_cutoff_energy=             gd.getNextNumber();
    PSF_ellipse_scale=             gd.getNextNumber();
    RPSF_min_mask_threshold=       gd.getNextNumber();
    INVERSE_FILTER=                gd.getNextBoolean();
    RPSF_sigma_to_radius=          gd.getNextNumber();
    RPSF_var_sigma_scale=          gd.getNextNumber();

    DECONCV_BLUR_INDIVIDUAL=       gd.getNextNumber();
    DECONCV_BLUR_CHECKER=          gd.getNextNumber();
    UPDATE_STATUS=                 gd.getNextBoolean();
    MASTER_DEBUG_LEVEL=      (int) gd.getNextNumber();
    return true;
  }

  public boolean aliasScissorsStackDialog() {
    int i;
    GenericDialog gd = new GenericDialog("De-bayer parameters");
    gd.addNumericField("Interpolation step between original kernels",      SPLIT_OVERSAMPLE,      0); //2
    gd.addNumericField("Debayer threshold (lower use default filtering)",  DEBAYER_THRESHOLD,     3) ; //=0.2; Measuring maximal spectral component amplitude in mid-frequencies. If below this - use default de-bayer mask
    gd.addNumericField("Debayer lo-pass relative width for green",         DEBAYER_WIDTH_GREEN,   3); //1.5
    gd.addNumericField("Debayer lo-pass relative width for red/blue",      DEBAYER_WIDTH_REDBLUE, 3); //1.5

    gd.addNumericField("Debayer mask - power for the ampliutude",      DEBAYER_GAMMA,   3); //0.3  
    gd.addNumericField("Min frequency radius from zero (fraction)",    DEBAYER_RZ,      3); //0.0 - recalculate to space domain (pixels)
    gd.addNumericField("Min frequency radius from aliases (fraction)", DEBAYER_RA,      3); //0.25
    gd.addNumericField("Frequency sigma (reduce far pixels, fraction)",DEBAYER_SIGMA,   3); //0.5
    gd.addNumericField("Wings decay (frequency pixels, fraction)",     DEBAYER_DECAY,   3); //0.5
    gd.addNumericField("Fartherst absolute maximum on a ray to count, fraction)",  DEBAYER_FARTHEST_MAX,   3); //0.5 fartherst absolute maximum on a ray to count
    gd.addNumericField("Divide data by radius to this power",          DEBAYER_RADIUS_POWER,   3); //0.5  // divide ray values by the radius to this power
    gd.addNumericField("Relative alais strength to mask out point",    DEBAYER_MAINTOALIAS,3); //0.5; // relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
    gd.addNumericField("Gaussian blur sigma for the alias-rejecting masks",  DEBAYER_MASK_BLUR,   3); //2.0

    gd.addCheckbox    ("Combine alias-reject 'scissors' with lo-pass filter for greens",DEBAYER_LO_GREEN);    // combine alias-reject "scissors" with lopass filter for greens
    gd.addCheckbox    ("Same, but only for generating red/blue mask ",DEBAYER_LO_POSTGREEN);                  // combine alias-reject "scissors" with lopass filter for greens applied to red/blue
    gd.addCheckbox    ("Combine alias-reject 'scissors' with lo-pass filter for red/blue",DEBAYER_LO_REDBLUE);// combine alias-reject "scissors" with lopass filter for greens applied to red/blue
    gd.addNumericField("Debayer FFT Size (128)",                       DEBAYER_FFT_SIZE,0); // 128
    gd.addCheckbox    ("Update ImageJ status",                                            UPDATE_STATUS);
    gd.addNumericField("Debug Level:",                            MASTER_DEBUG_LEVEL,      0);
    gd.showDialog();
    if (gd.wasCanceled()) return false;
    SPLIT_OVERSAMPLE=  (int) gd.getNextNumber();
    DEBAYER_THRESHOLD=       gd.getNextNumber();
    DEBAYER_WIDTH_GREEN=     gd.getNextNumber();
    DEBAYER_WIDTH_REDBLUE=   gd.getNextNumber();
    DEBAYER_GAMMA=           gd.getNextNumber();
    DEBAYER_RZ=              gd.getNextNumber();
    DEBAYER_RA=              gd.getNextNumber();
    DEBAYER_SIGMA=           gd.getNextNumber();
    DEBAYER_DECAY=           gd.getNextNumber();
    DEBAYER_FARTHEST_MAX=    gd.getNextNumber();
    DEBAYER_RADIUS_POWER=    gd.getNextNumber();
    DEBAYER_MAINTOALIAS=     gd.getNextNumber();
    DEBAYER_MASK_BLUR=       gd.getNextNumber();
    DEBAYER_LO_GREEN=        gd.getNextBoolean();
    DEBAYER_LO_POSTGREEN=    gd.getNextBoolean();
    DEBAYER_LO_REDBLUE=      gd.getNextBoolean();

    DEBAYER_FFT_SIZE=1;
    for (i=(int) gd.getNextNumber(); i >1; i>>=1) DEBAYER_FFT_SIZE <<=1; /* make it to be power of 2 */
    UPDATE_STATUS=           gd.getNextBoolean();
    MASTER_DEBUG_LEVEL=(int) gd.getNextNumber();
    return true;
  }


  public boolean showCombinePairDialog() {
    GenericDialog gd = new GenericDialog("Combining high-re and low-res images");

    gd.addNumericField("Nonlinear filter mask sigma" ,          NONLIN_SIGMA,            3); // 5.0-  sigma for the nonlinear filtering (higher the sigma, father from the edges extends the PSF correection)
    gd.addNumericField("Nonlinear filter mask min. level",      NONLIN_MIN,              3); // 0.01  minimal low-pass filtered squared difference between the corrected and original pixels to trigger sharpness enhancement
    gd.addNumericField("Nonlinear filter mask max. level",      NONLIN_MAX,              3); // 0.15 squared low-pass filtered difference between the corrected and original pixels, so abopve that level 100% corrected image is used
    gd.addNumericField("Nonlinear filter threshold",            NONLIN_THRESHOLD,        3); // 0.01 when blurred intensity is below this value, use it as a denominator


    gd.addNumericField("Debug Level:",                          MASTER_DEBUG_LEVEL,      0);
    gd.showDialog();
    if (gd.wasCanceled()) return false;

    NONLIN_SIGMA=            gd.getNextNumber();
    NONLIN_MIN=              gd.getNextNumber();
    NONLIN_MAX=              gd.getNextNumber();
    NONLIN_THRESHOLD=        gd.getNextNumber();

    MASTER_DEBUG_LEVEL=      (int) gd.getNextNumber();
    return true;
  }


  public boolean showDeconvDialog() {
    GenericDialog gd = new GenericDialog("Deconv parameters");
    gd.addNumericField("Gaussian blur for individual colors",   DECONCV_BLUR_INDIVIDUAL, 3); //1.8
    gd.addNumericField("Gaussian blur for diagonal greens",     DECONCV_BLUR_DIAGONAL,   3); //2.0
    gd.addNumericField("Gaussian blur for checkerboard greens", DECONCV_BLUR_CHECKER,    3); //1.4

    gd.addNumericField("Reversed PSF sigma to radius ratio (0 to disable)",              RPSF_sigma_to_radius, 3); //0.2- variable blurring - sigma will be proportional distance from the center
    gd.addNumericField("Reversed PSF: scale variable sigma (in the center) from the uniform one",  RPSF_var_sigma_scale, 3); //=0.8; // reduce variable sigma in the center from uniuform one
    gd.addNumericField("Reject sampliang aliases for individual colors",   DECONCV_ALIASREJ_INDIVIDUAL, 3); //1.0
    gd.addNumericField("Reject sampliang aliases for diagonal greens",     DECONCV_ALIASREJ_DIAGONAL,   3); //1.0
    gd.addNumericField("Reject sampliang aliases for checkerboard greens", DECONCV_ALIASREJ_CHECKER,    3); //1.5

    gd.addNumericField("Nonlinear filter mask sigma" ,          NONLIN_SIGMA,            3); // 5.0-  sigma for the nonlinear filtering (higher the sigma, father from the edges extends the PSF correection)
    gd.addNumericField("Nonlinear filter mask min. level",      NONLIN_MIN,              3); // 0.01  minimal low-pass filtered squared difference between the corrected and original pixels to trigger sharpness enhancement
    gd.addNumericField("Nonlinear filter mask max. level",      NONLIN_MAX,              3); // 0.15 squared low-pass filtered difference between the corrected and original pixels, so abopve that level 100% corrected image is used
    gd.addNumericField("Nonlinear filter threshold",            NONLIN_THRESHOLD,        3); // 0.01 when blurred intensity is below this value, use it as a denominator

    gd.addNumericField("Results gamma",                         DECONCV_GAMMA,           3); //0.5
//    gd.addNumericField("Color saturation",                      DECONCV_SATURATION,      3); //2.0
    gd.addNumericField("Color saturation headroom",             DECONCV_HEADROOM,        3); //0.1
    gd.addNumericField("Color balance lowPerc",                 DECONCV_LOWPERC,         3); //0.005
    gd.addNumericField("Color balance lowLevel",                DECONCV_LOWLEVEL,        3); //0.05
    gd.addNumericField("Color balance highPerc",                DECONCV_HIGHPERC,        3); //0.005
    gd.addNumericField("Color balance highLevel",               DECONCV_HIGHLEVEL,       3); //0.95
    gd.addCheckbox    ("Mask out Bayer aliases",                MASK_BAYER_ALIASES);
    gd.addCheckbox    ("Calculate Bayer weights",               CALC_BAYER_WEIGHTS);

    gd.addNumericField("Debayer threshold (lower use default filtering)",  DEBAYER_THRESHOLD,3) ; //=0.2; Measuring maximal spectral component amplitude in mid-frequencies. If below this - use default de-bayer mask
    gd.addNumericField("Debayer lo-pass relative width for green",         DEBAYER_WIDTH_GREEN,   3); //1.5
    gd.addNumericField("Debayer lo-pass relative width for red/blue",      DEBAYER_WIDTH_REDBLUE, 3); //1.5
    gd.addNumericField("Debayer mask - power for the ampliutude",      DEBAYER_GAMMA,   3); //0.3  
    gd.addNumericField("Min frequency radius from zero (fraction)",    DEBAYER_RZ,      3); //0.0 - recalculate to space domain (pixels)
    gd.addNumericField("Min frequency radius from aliases (fraction)", DEBAYER_RA,      3); //0.25
    gd.addNumericField("Frequency sigma (reduce far pixels, fraction)",DEBAYER_SIGMA,   3); //0.5
    gd.addNumericField("Wings decay (frequency pixels, fraction)",     DEBAYER_DECAY,   3); //0.5
    gd.addNumericField("Fartherst absolute maximum on a ray to count, fraction)",  DEBAYER_FARTHEST_MAX,   3); //0.5 fartherst absolute maximum on a ray to count
    gd.addNumericField("Didide data by radius to this power",          DEBAYER_RADIUS_POWER,   3); //0.5  // divide ray values by the radius to this power
    gd.addNumericField("Gaussian blur sigma for the alias-rejecting masks",  DEBAYER_MASK_BLUR,   3); //2.0
    gd.addCheckbox    ("Combine alias-reject 'scissors' with lo-pass filter for greens",DEBAYER_LO_GREEN);    // combine alias-reject "scissors" with lopass filter for greens
    gd.addCheckbox    ("Same, but only for generating red/blue mask ",DEBAYER_LO_POSTGREEN);                  // combine alias-reject "scissors" with lopass filter for greens applied to red/blue
    gd.addCheckbox    ("Combine alias-reject 'scissors' with lo-pass filter for red/blue",DEBAYER_LO_REDBLUE);// combine alias-reject "scissors" with lopass filter for greens applied to red/blue

  
    gd.addNumericField("Debug Level:",                          MASTER_DEBUG_LEVEL,      0);
    gd.showDialog();
    if (gd.wasCanceled()) return false;
    DECONCV_BLUR_INDIVIDUAL= gd.getNextNumber();
    DECONCV_BLUR_DIAGONAL=   gd.getNextNumber();
    DECONCV_BLUR_CHECKER=    gd.getNextNumber();

    RPSF_sigma_to_radius=    gd.getNextNumber();
    RPSF_var_sigma_scale=    gd.getNextNumber();


    DECONCV_ALIASREJ_INDIVIDUAL= gd.getNextNumber();
    DECONCV_ALIASREJ_DIAGONAL=   gd.getNextNumber();
    DECONCV_ALIASREJ_CHECKER=    gd.getNextNumber();

    NONLIN_SIGMA=            gd.getNextNumber();
    NONLIN_MIN=              gd.getNextNumber();
    NONLIN_MAX=              gd.getNextNumber();
    NONLIN_THRESHOLD=        gd.getNextNumber();

    DEBAYER_THRESHOLD=       gd.getNextNumber();
    DEBAYER_WIDTH_GREEN=     gd.getNextNumber();
    DEBAYER_WIDTH_REDBLUE=   gd.getNextNumber();

    DECONCV_GAMMA=           gd.getNextNumber();
//    DECONCV_SATURATION=      gd.getNextNumber();
    DECONCV_HEADROOM=        gd.getNextNumber();
    DECONCV_LOWPERC=         gd.getNextNumber();
    DECONCV_LOWLEVEL=        gd.getNextNumber();
    DECONCV_HIGHPERC=        gd.getNextNumber();
    DECONCV_HIGHLEVEL=       gd.getNextNumber();
    MASK_BAYER_ALIASES=      gd.getNextBoolean();
    CALC_BAYER_WEIGHTS=      gd.getNextBoolean();

    DEBAYER_GAMMA=          gd.getNextNumber();
    DEBAYER_RZ=             gd.getNextNumber();
    DEBAYER_RA=             gd.getNextNumber();
    DEBAYER_SIGMA=          gd.getNextNumber();
    DEBAYER_DECAY=          gd.getNextNumber();
    DEBAYER_FARTHEST_MAX=   gd.getNextNumber();
    DEBAYER_RADIUS_POWER=   gd.getNextNumber();
    DEBAYER_MASK_BLUR=      gd.getNextNumber();
    DEBAYER_LO_GREEN= gd.getNextBoolean();
    DEBAYER_LO_POSTGREEN= gd.getNextBoolean();
    DEBAYER_LO_REDBLUE= gd.getNextBoolean();

    MASTER_DEBUG_LEVEL=      (int) gd.getNextNumber();
    return true;
  }

  public boolean aliasScissorsDialog() {
    int i;
    GenericDialog gd = new GenericDialog("De-bayer parameters");
    gd.addNumericField("Debayer threshold (lower use default filtering)",  DEBAYER_THRESHOLD,3) ; //=0.2; Measuring maximal spectral component amplitude in mid-frequencies. If below this - use default de-bayer mask
    gd.addNumericField("Debayer lo-pass relative width for green",         DEBAYER_WIDTH_GREEN,   3); //1.5
    gd.addNumericField("Debayer lo-pass relative width for red/blue",      DEBAYER_WIDTH_REDBLUE, 3); //1.5
    gd.addNumericField("Debayer mask - power for the ampliutude",      DEBAYER_GAMMA,   3); //0.3  
    gd.addNumericField("Min frequency radius from zero (fraction)",    DEBAYER_RZ,      3); //0.0 - recalculate to space domain (pixels)
    gd.addNumericField("Min frequency radius from aliases (fraction)", DEBAYER_RA,      3); //0.25
    gd.addNumericField("Frequency sigma (reduce far pixels, fraction)",DEBAYER_SIGMA,   3); //0.5
    gd.addNumericField("Wings decay (frequency pixels, fraction)",     DEBAYER_DECAY,   3); //0.5
    gd.addNumericField("Fartherst absolute maximum on a ray to count, fraction)",  DEBAYER_FARTHEST_MAX,   3); //0.5 fartherst absolute maximum on a ray to count
    gd.addNumericField("Divide data by radius to this power",          DEBAYER_RADIUS_POWER,   3); //0.5  // divide ray values by the radius to this power
    gd.addNumericField("Relative alais strength to mask out point",    DEBAYER_MAINTOALIAS,3); //0.5; // relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
    gd.addNumericField("Gaussian blur sigma for the alias-rejecting masks",  DEBAYER_MASK_BLUR,   3); //2.0

    gd.addCheckbox    ("Combine alias-reject 'scissors' with lo-pass filter for greens",DEBAYER_LO_GREEN);    // combine alias-reject "scissors" with lopass filter for greens
    gd.addCheckbox    ("Same, but only for generating red/blue mask ",DEBAYER_LO_POSTGREEN);                  // combine alias-reject "scissors" with lopass filter for greens applied to red/blue
    gd.addCheckbox    ("Combine alias-reject 'scissors' with lo-pass filter for red/blue",DEBAYER_LO_REDBLUE);// combine alias-reject "scissors" with lopass filter for greens applied to red/blue
    gd.addNumericField("Debayer FFT Size (64)",                       DEBAYER_FFT_SIZE,0); // 128

    gd.addNumericField("Debug Level:",                            MASTER_DEBUG_LEVEL,      0);
    gd.showDialog();
    if (gd.wasCanceled()) return false;
    DEBAYER_THRESHOLD=       gd.getNextNumber();
    DEBAYER_WIDTH_GREEN=     gd.getNextNumber();
    DEBAYER_WIDTH_REDBLUE=   gd.getNextNumber();
    DEBAYER_GAMMA=          gd.getNextNumber();
    DEBAYER_RZ=             gd.getNextNumber();
    DEBAYER_RA=             gd.getNextNumber();
    DEBAYER_SIGMA=          gd.getNextNumber();
    DEBAYER_DECAY=          gd.getNextNumber();
    DEBAYER_FARTHEST_MAX=   gd.getNextNumber();
    DEBAYER_RADIUS_POWER=   gd.getNextNumber();
    DEBAYER_MAINTOALIAS=   gd.getNextNumber();
    DEBAYER_MASK_BLUR=      gd.getNextNumber();

    DEBAYER_LO_GREEN= gd.getNextBoolean();
    DEBAYER_LO_POSTGREEN= gd.getNextBoolean();
    DEBAYER_LO_REDBLUE= gd.getNextBoolean();

    DEBAYER_TEST_FROMIMAGE= gd.getNextBoolean();
    DEBAYER_TEST_MASKSPLIT= gd.getNextBoolean();
    DEBAYER_FFT_SIZE=1;
    for (i=(int) gd.getNextNumber(); i >1; i>>=1) DEBAYER_FFT_SIZE <<=1; /* make it to be power of 2 */

    MASTER_DEBUG_LEVEL=      (int) gd.getNextNumber();
    return true;
  }

  public boolean showDeBayerDialog() {
    int i;
    GenericDialog gd = new GenericDialog("De-bayer parameters");
    gd.addNumericField("Debayer threshold (lower use default filtering)",  DEBAYER_THRESHOLD,3) ; //=0.2; Measuring maximal spectral component amplitude in mid-frequencies. If below this - use default de-bayer mask
    gd.addNumericField("Debayer lo-pass relative width for green",         DEBAYER_WIDTH_GREEN,   3); //1.5
    gd.addNumericField("Debayer lo-pass relative width for red/blue",      DEBAYER_WIDTH_REDBLUE, 3); //1.5
    gd.addNumericField("Debayer mask - power for the ampliutude",      DEBAYER_GAMMA,   3); //0.3  
    gd.addNumericField("Min frequency radius from zero (fraction)",    DEBAYER_RZ,      3); //0.0 - recalculate to space domain (pixels)
    gd.addNumericField("Min frequency radius from aliases (fraction)", DEBAYER_RA,      3); //0.25
    gd.addNumericField("Frequency sigma (reduce far pixels, fraction)",DEBAYER_SIGMA,   3); //0.5
    gd.addNumericField("Wings decay (frequency pixels, fraction)",     DEBAYER_DECAY,   3); //0.5
    gd.addNumericField("Fartherst absolute maximum on a ray to count, fraction)",  DEBAYER_FARTHEST_MAX,   3); //0.5 fartherst absolute maximum on a ray to count
    gd.addNumericField("Divide data by radius to this power",          DEBAYER_RADIUS_POWER,   3); //0.5  // divide ray values by the radius to this power
    gd.addNumericField("Relative alais strength to mask out point",    DEBAYER_MAINTOALIAS,3); //0.5; // relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
    gd.addNumericField("Gaussian blur sigma for the alias-rejecting masks",  DEBAYER_MASK_BLUR,   3); //2.0

    gd.addCheckbox    ("Combine alias-reject 'scissors' with lo-pass filter for greens",DEBAYER_LO_GREEN);    // combine alias-reject "scissors" with lopass filter for greens
    gd.addCheckbox    ("Same, but only for generating red/blue mask ",DEBAYER_LO_POSTGREEN);                  // combine alias-reject "scissors" with lopass filter for greens applied to red/blue
    gd.addCheckbox    ("Combine alias-reject 'scissors' with lo-pass filter for red/blue",DEBAYER_LO_REDBLUE);// combine alias-reject "scissors" with lopass filter for greens applied to red/blue

    gd.addCheckbox    ("Start from image (false - amplitude)",         DEBAYER_TEST_FROMIMAGE);
    gd.addCheckbox    ("Mask source with window and green pattern",    DEBAYER_TEST_MASKSPLIT);
    gd.addNumericField("Debayer FFT Size (128)",                       DEBAYER_FFT_SIZE,0); // 128

    gd.addNumericField("Debug Level:",                            MASTER_DEBUG_LEVEL,      0);
    gd.showDialog();
    if (gd.wasCanceled()) return false;
    DEBAYER_THRESHOLD=       gd.getNextNumber();
    DEBAYER_WIDTH_GREEN=     gd.getNextNumber();
    DEBAYER_WIDTH_REDBLUE=   gd.getNextNumber();
    DEBAYER_GAMMA=          gd.getNextNumber();
    DEBAYER_RZ=             gd.getNextNumber();
    DEBAYER_RA=             gd.getNextNumber();
    DEBAYER_SIGMA=          gd.getNextNumber();
    DEBAYER_DECAY=          gd.getNextNumber();
    DEBAYER_FARTHEST_MAX=   gd.getNextNumber();
    DEBAYER_RADIUS_POWER=   gd.getNextNumber();
    DEBAYER_MAINTOALIAS=   gd.getNextNumber();
    DEBAYER_MASK_BLUR=      gd.getNextNumber();

    DEBAYER_LO_GREEN= gd.getNextBoolean();
    DEBAYER_LO_POSTGREEN= gd.getNextBoolean();
    DEBAYER_LO_REDBLUE= gd.getNextBoolean();

    DEBAYER_TEST_FROMIMAGE= gd.getNextBoolean();
    DEBAYER_TEST_MASKSPLIT= gd.getNextBoolean();
    DEBAYER_FFT_SIZE=1;
    for (i=(int) gd.getNextNumber(); i >1; i>>=1) DEBAYER_FFT_SIZE <<=1; /* make it to be power of 2 */

    MASTER_DEBUG_LEVEL=      (int) gd.getNextNumber();
    return true;
  }





  public boolean showColorProcessDialog() {
    GenericDialog gd = new GenericDialog("Color processing parameters");
    gd.addNumericField("Color balance, red-to-green",             BALANCE_RED,     3); //1.8
    gd.addNumericField("Color balance, blue-to-green",            BALANCE_BLUE,    3); //1.8
    gd.addNumericField("Gain green",                              GAIN_GREEN,      3); //1.8
    gd.addNumericField("Weight scale RED  (which color to use)",  WEIGHT_SCALE_R,  3); //1.8
    gd.addNumericField("Weight scale BLUE (which color to use)",  WEIGHT_SCALE_B,  3); //1.8
    gd.addNumericField("Color low-pass sigma (pixels)",           COLOR_SIGMA,     3); //1.8
    gd.addNumericField("YCbCr gamma",                             YCbCr_Gamma,     3); //0.53
    gd.addNumericField("Minimal linear value to apply gammaY",    YCbCr_minLin,    3); //0.53
    gd.addNumericField("YCbCr Kb",                                YCbCr_Kb,        3); //0.114
    gd.addNumericField("YCbCr Kr",                                YCbCr_Kr,        3); //0.299
    gd.addNumericField("Color saturation",                        COLOR_SATURATION,3); //2.0
    gd.addCheckbox    ("Use first approximation for Y",           USE_FIRST_Y);

    gd.addNumericField("Debug Level:",                          MASTER_DEBUG_LEVEL,      0);
    gd.showDialog();
    if (gd.wasCanceled()) return false;

    BALANCE_RED=                gd.getNextNumber();
    BALANCE_BLUE=               gd.getNextNumber();
    GAIN_GREEN=                 gd.getNextNumber();
    WEIGHT_SCALE_R=             gd.getNextNumber();
    WEIGHT_SCALE_B=             gd.getNextNumber();
    COLOR_SIGMA=                gd.getNextNumber();
    YCbCr_Gamma=                gd.getNextNumber();
    YCbCr_minLin=               gd.getNextNumber();
    YCbCr_Kb=                   gd.getNextNumber();
    YCbCr_Kr=                   gd.getNextNumber();
    COLOR_SATURATION=           gd.getNextNumber();
    USE_FIRST_Y=                gd.getNextBoolean();
    MASTER_DEBUG_LEVEL=      (int) gd.getNextNumber();
    return true;
  }

  public boolean showPSFDialog() {
    int i;
    GenericDialog gd = new GenericDialog("PSF measurement parameters");
//    gd.addNumericField("FFT_Size:",                                                   FFTSize, 0);
//    gd.addNumericField("Gauss to window ratio (0 - use Hamming:",                     GAUSS_WIDTH, 3);//0.4; //0 - use Hamming window
//    gd.addNumericField("Frequncy subdivide:",                                         subDivFreq, 0);
    gd.addNumericField("Invert deconvolution if less than",                           deconvInvert, 3);
//    gd.addNumericField("Gamma value for pattern frequency measurement:",              corr_gamma, 3);
//    gd.addNumericField("Location diff between spectrum and correlation:",             diffSpectrCorr, 0); //2
//    gd.addNumericField("Shrink clusters after initial separation (0.0 for sqrt(5*s), negative - absolute size):",   shrinkClusters, 3); // 0.0; //  Shrink clusters by this ratio (remove lowest) after initial separation
//    gd.addNumericField("Number secondary maximums to try:",                           multiplesToTry, 0); //4 -  try this number of maximums proportionally farther from 0,0 than the two closest (increase precision)
//    gd.addNumericField("Deviation",                                                   deviation, 3); // 1.0 - when looking for maximums - maximal distance from predicted from the lower order one
//    gd.addNumericField("Deviation steps",                                             deviation_steps, 0); // 6; maximal iterations when looking for local maximum
    gd.addNumericField("Decimate PSF before binning:",                                PSF_subpixel,    0); // OTF sub-pixel decimation
//    gd.addNumericField("Show multiple OTF instances:",                                OTF_multiple,    0); // (debug feature, normally 0) // 0 - use each pixel once, 1 - add first negatives (4), 2 - second positives()4)
    gd.addNumericField("Minimal PSF contrast to use",                                 PSF_minContrast, 3); // 0.1 - minimal instance contrast to use in binning (compared to the one at [0,0]
    gd.addNumericField("PSF cell reduction from centers of negative clones",          PSF_windowFrac, 3); // 0.9 reduce the PSF cell size to this part of the area connecting first negative clones
    gd.addCheckbox    ("Multiply PSF cell by Hamming window",                         PSF_useHamming); //  true;
    gd.addCheckbox    ("Force PSF center- symmetrical (around centroid)",              PSF_symm180); //  true; // make OTF center-symmetrical (around centroid that is defined by lateral chromatic aberration)
    gd.addNumericField("OTF FFT size (initilaly will be increased by decimation factor)",  OTF_FFT_size, 0); // 16  initially will be increased by PSF_subpixel
    gd.addNumericField("OTF cutoff energy (used to determine bounding ellipse)",      OTF_cutoff_energy, 3); //0.6; use frequency points that have OTF_cutoff_energy of the total to determine ellipse for limiting frequency responce
    gd.addNumericField("OTF size of elliptical window relative to cluster size",      OTF_ellipse_scale, 3); //1.5;  // size of elliptical window relative to the cluster defined by OTF_cutoff_energy
    gd.addCheckbox    ("Use Gauss window (instead of Hamming) for Ellipse mask",      OTF_ellipse_gauss); //   // size of elliptical window relative to the cluster defined by OTF_cutoff_energy
    gd.addNumericField("OTF deconvolution parameter ",                                OTF_deconvInvert, 3); //0.01; // when FFT component is less than this fraction of the maximal value, replace 1/z with Z

    gd.addNumericField("OTF zero frequency size on power spectrum ",                  OTF_zerofreq_size, 3); //2.0;
    gd.addNumericField("OTF smouth PS to generate alias rejection mask (0 - none)",   OTF_smoothPS,      3); //2.5 - smooth model PS for rejecting aliases (0 - no smouth, >0 additional Gauss )

    gd.addNumericField("OTF relative high value of PS for rejection mask ",           OTF_threshold_high, 3); //0.1
    gd.addNumericField("OTF relative low  value of PS for rejection mask ",           OTF_threshold_low,  3); //0.01; // when FFT component is less than this fraction of the maximal value, replace 1/z with Z


    gd.addNumericField("XMASK threshold: <0 - disabled, 0 - use amplitude, 0..1 - make binary ",           XMASK_threshold,  3); // 0.01; // if 0 use normalize amplitude, if 0..1 - make binary: 1.0 if > threshold, 0.0 - otherwise
    gd.addNumericField("XMASK: PS filtering radius (derive from model later) ",        XMASK_radius,     3); // 40.0; // low-pass result with low pass filter (should be later defined automatically)
    gd.addNumericField("XMASK:histogram size (use FFTSize/2?) ",                       XMASK_hsize,      0); // 128; // 2d histogram size (size/2 probably a good guess),
    gd.addNumericField("XMASK: percentile to determine weak directions ",              XMASK_percentile, 3); // 0.1; // use this percentile (0.0..1.0)) value for given radius as a target
    gd.addNumericField("XMASK: boost low components up to this gain (too high can go to 0) ",           XMASK_maxGain,    3); // 5.0; // maximal gain for low components
    gd.addNumericField("XMASK: exaggerate correction mask (Math.pow()) ",              XMASK_exaggerate,  3); // 1.0 exaggerate correction mask

    gd.addNumericField("PSF separation: low-pass filter width (to PSF half-period) ",  PSF_smoothSeparate,  3); // 0.125 low pass filter width (relative to PSF pitch) when separation individual PSF
    gd.addNumericField("PSF: Do not try to compensate clones if they are too sharp (sigma as part of smoothing one)", PSF_thresholdSeparate,  3); //0.1   do not try to compensate for adjacent PSF clones if model sigma is less than this fraction of the one used for smoothing

    gd.addNumericField("PSF separation: threshold to find the PSF maximum",            PSF_topCenter,  3); // 0.75 consider only points above this fraction of the peak to find the centroid
    gd.addCheckbox    ("Remove negative PSF areas while removing clones",              PSF_removeNegtative); //true; // remove negative when separating composite PSF (will need low-pass filtering)

    gd.addNumericField("PSF variable Gauss blurring (farther from center, higher the sigma",    PSF_sigmaToRadius,3); // 0.4 variable-sigma blurring to reduce high frequencies more for the pixels farther from the PSF center
    gd.addNumericField("PSF wings energy (searching for ellipse approximation)",       PSF_wings_energy, 3); //  0.8 fraction of energy in the pixels to be used
    gd.addNumericField("PSF wings ellipse scale (multiply PSF by elliptical gaussian)",PSF_wings_ellipse_scale, 3);// 2.0 increase wings cutoff ellipse by this from one defined by the  cutoff energy
    gd.addNumericField("PSF kernel size (side of square array to be stored) ",         PSF_kernel_size,  0); // 64 -kernel (to be stored) size (per color component)

//    gd.addNumericField("PSF wings threshold ",    PSF_wings_min_mask_threshold, 3); //=0.003 zero output element if elliptical Gauss mask is below this threshold

    gd.addCheckbox     ("Ignore lateral chromatic aberrations, center PSF",              PSF_ignoreChromatic); // true; // ignore lateral chromatic aberration (center OTF to 0,0)
    gd.addCheckbox     ("Fold high frequencies to low before inverse FFT of inverted OTF",OTF_fold); // false; // fold high frequency to lower when downsampling pixels (before inverse FFT)
    gd.addNumericField("PSF cutoff energy (used to determine bounding ellipse)",          PSF_cutoff_energy, 3); //0.9;  // Limit result kernel to proportional of the PSF, calculate initial cluster shape by this cutoff energy
    gd.addNumericField("Reversed PSF elliptical mask relative to the direct PSF ellipse", PSF_ellipse_scale, 3); //1.5;  // size of elliptical window to limuit reverse PSF as proportional to direct one
    gd.addNumericField("Reversed PSF mask thershold (completely ignore far pixels)",      RPSF_min_mask_threshold, 3); //0.01; // completely zero reversed kernel elements where elliptical mask is below this threshold
//    gd.addNumericField("Reversed PSF sigma to radius ratio (0 to disable)",              RPSF_sigma_to_radius, 3); //0.2- variable blurring - sigma will be proportional distance from the center
//    gd.addNumericField("Reversed PSF minimal sigma (in the center)",                      RPSF_min_sigma, 3); //1.0;  blurring in the center sigma(r)=sqrt((RPSF_sigma_to_radius*r)^2)+RPSF_min_sigma^2)
    gd.addCheckbox     ("Calculate kernel for direct (not inverted) OTF",                 forwardOTF); //  true;
    gd.addNumericField("Debug Level:",                                                    MASTER_DEBUG_LEVEL, 0);
    gd.showDialog();
    if (gd.wasCanceled()) return false;
//    FFTSize=4;
//    for (i=(int) gd.getNextNumber(); i >4; i>>=1) FFTSize <<=1; /* make FFTSize to be power of 2*/
//    GAUSS_WIDTH=                   gd.getNextNumber();
//    subDivFreq=              (int) gd.getNextNumber();
    deconvInvert=                  gd.getNextNumber();  //0.05; // when FFT component is lass than this fraction of the maximal value, replace 1/z with Z
//    corr_gamma=                    gd.getNextNumber();  //0.2; lower the value - higher harmonics will participate in pattern frequency measurements
//    diffSpectrCorr=          (int) gd.getNextNumber();
//    shrinkClusters=                gd.getNextNumber(); // 0.5; //  Shrink clusters by this ratio (remove lowest) after initial separation
//    multiplesToTry=          (int) gd.getNextNumber();
//    deviation=                     gd.getNextNumber();
//    deviation_steps=         (int) gd.getNextNumber();
    PSF_subpixel=            (int) gd.getNextNumber();
//    OTF_multiple=            (int) gd.getNextNumber();
    PSF_minContrast=               gd.getNextNumber();
    PSF_windowFrac=                gd.getNextNumber();
    PSF_useHamming=                gd.getNextBoolean();
    PSF_symm180=                   gd.getNextBoolean();
    OTF_FFT_size=1;
    for (i=(int) gd.getNextNumber(); i >1; i>>=1) OTF_FFT_size <<=1; /* make it to be power of 2 */
    OTF_cutoff_energy=             gd.getNextNumber();
    OTF_ellipse_scale=             gd.getNextNumber();
    OTF_ellipse_gauss=             gd.getNextBoolean();
    OTF_deconvInvert=              gd.getNextNumber();

    OTF_zerofreq_size=             gd.getNextNumber();
    OTF_smoothPS=                  gd.getNextNumber();
    OTF_threshold_high=            gd.getNextNumber();
    OTF_threshold_low=             gd.getNextNumber();

    XMASK_threshold=               gd.getNextNumber();
    XMASK_radius=                  gd.getNextNumber();
    XMASK_hsize=             (int) gd.getNextNumber();
    XMASK_percentile=              gd.getNextNumber();
    XMASK_maxGain=                 gd.getNextNumber();
    XMASK_exaggerate=              gd.getNextNumber();

    PSF_smoothSeparate=            gd.getNextNumber();
    PSF_thresholdSeparate=         gd.getNextNumber();
    PSF_topCenter=                 gd.getNextNumber();
    PSF_removeNegtative=           gd.getNextBoolean();
    PSF_sigmaToRadius=             gd.getNextNumber();
    PSF_wings_energy=              gd.getNextNumber();
    PSF_wings_ellipse_scale=       gd.getNextNumber();
    PSF_kernel_size=         (int) gd.getNextNumber();

    PSF_ignoreChromatic=           gd.getNextBoolean();
    OTF_fold=                      gd.getNextBoolean();
    PSF_cutoff_energy=             gd.getNextNumber();
    PSF_ellipse_scale=             gd.getNextNumber();
    RPSF_min_mask_threshold=       gd.getNextNumber();
//    RPSF_sigma_to_radius=          gd.getNextNumber();
//    RPSF_min_sigma=                gd.getNextNumber();

    forwardOTF=                    gd.getNextBoolean();
    MASTER_DEBUG_LEVEL=      (int) gd.getNextNumber();
//    if (DEBUG_LEVEL>1)   System.out.println("showPSFDialog: DEBUG_LEVEL="+DEBUG_LEVEL+" MASTER_DEBUG_LEVEL="+MASTER_DEBUG_LEVEL+ " RPSF_min_mask_threshold="+RPSF_min_mask_threshold);
    return true;
  }

  public boolean showConfigureDialog() {
    int i;
    GenericDialog gd = new GenericDialog("Parameters");
    gd.addStringField ("Filename prefix:                   ", jp4_instance.getTitle(), 20);
    gd.addStringField ("Camera address:                    ", jp4_instance.getURL(), 20);
    gd.addNumericField("FFT_Size:",                           FFTSize, 0);
    gd.addNumericField("FFT_Size for mapping image (areas covered):", mapFFTSize, 0); //64 used to find where grid covers the image

    gd.addNumericField("Gauss to window ratio (0 - use Hamming:",GAUSS_WIDTH, 3);//0.4; //0 - use Hamming window
    gd.addNumericField("FFT Overlap (=FFT_Size recommended):", FFTOverlap, 0);

    gd.addNumericField("Frequncy subdivide:",                 subDivFreq, 0);
    gd.addNumericField("Invert deconvolution if less than",   deconvInvert, 3);
    gd.addNumericField("Gamma value for pattern frequency measurement:",   corr_gamma, 3);
    gd.addNumericField("Sigma value for gauss high-pass pattern filtering:",   corr_sigma, 3);// 1.5; // pattern detection: high-pass filter (0.0 - none) gamma(PS)

    gd.addNumericField("Location diff between spectrum and correlation:",  diffSpectrCorr, 0); //2
    gd.addNumericField("Shrink clusters after initial separation (0.0 for sqrt(5*s), negative - absolute size):",   shrinkClusters, 3); // 0.0; //  Shrink clusters by this ratio (remove lowest) after initial separation
    gd.addNumericField("Number secondary maximums to try:",   multiplesToTry, 0); //4 -  try this number of maximums proportionally farther from 0,0 than the two closest (increase precision)

    gd.addNumericField("Deviation",                           deviation, 3); // 1.0 - when looking for maximums - maximal distance from predicted from the lower order one
    gd.addNumericField("Deviation steps",                     deviation_steps, 0); // 6; maximal iterations when looking for local maximum


    gd.addNumericField("High-pass filter for correlation with the model",model_highpass, 3); //1.5 model correlation high-pass filter (relative to pattern fundamental frequency - average of 2)
    gd.addNumericField("Model correlation relative ring width",corrRingWidth, 3);   //0.4, ring (around r=0.5 dist to opposite corr) width , center circle r=0.5*corrRingWidth
    gd.addNumericField("Minimal pattern correlation contrast" ,minCorrContrast, 3);   //5.0; // Discrimination threshold between good and bad pattern correleation

    gd.addCheckbox     ("Calculate correction for RED component",               correctBayerRed); //  true;
    gd.addCheckbox     ("Calculate correction for BLUE component",              correctBayerBlue); //  true;
    gd.addCheckbox     ("Calculate correction for GREEN components (diagonal)", correctBayerDiagonal); // false;
    gd.addCheckbox     ("Calculate correction for GREEN components (checker)",  correctBayerChecker); //  true;
    gd.addCheckbox     ("Calculate correction for GREEN1 component (individual)",correctBayerGreen1); //  false;
    gd.addCheckbox     ("Calculate correction for GREEN2 component (individual)",correctBayerGreen2); //  false;
    gd.addNumericField ("Color component used as base for lateral chromatic (4,5)" ,referenceComponent, 0);   //4; // (change to 5 later) component to calculate lateral chromatic from (0 - G1, 1 - R, 2 - B, 3 - G2,4 - diagonal greens, 5 - checker greens)
    gd.addCheckbox     ("Equalize average values of the two greens in Bayer mosaic",equalizeGreens); //  true; Equalize average values of the two greens in Bayer mosaic

    gd.addNumericField("Half (size-1) of the result convolution kernel, individual colors:",  kernelHalSizeSingleColor, 0); //6 crop result coinvolution kernels to (kernelHalSizeSingleColor+1) * (kernelHalSizeSingleColor+1) for individual Bayer components
    gd.addNumericField("Half (size-1) of the result convolution kernel, combined greens:",    kernelHalSizeDualColors, 0); //9; /// crop result coinvolution kernels to (kernelHalSizeDualColors+1)  * (kernelHalSizeDualColors+1) for combined green Bayer components


    gd.addNumericField("Decimate PSF before binning:",                           PSF_subpixel,    0); // OTF sub-pixel decimation
    gd.addNumericField("Reversed PSF  kernel size",                              RPSF_kernel_size, 0); //   32 - size of deconvolution kernel

    gd.addNumericField("Debug Level:",      MASTER_DEBUG_LEVEL, 0);
    gd.showDialog();
    if (gd.wasCanceled()) return false;
    jp4_instance.setTitle(gd.getNextString());
    jp4_instance.setURL(gd.getNextString());
    FFTSize=4;
    for (i=(int) gd.getNextNumber(); i >4; i>>=1) FFTSize <<=1; /* make FFTSize to be power of 2*/
    mapFFTSize=4;
    for (i=(int) gd.getNextNumber(); i >4; i>>=1) mapFFTSize <<=1; /* make FFTSize to be power of 2*/
    GAUSS_WIDTH=              gd.getNextNumber();
    FFTOverlap=         (int) gd.getNextNumber();
    subDivFreq=         (int) gd.getNextNumber();
    deconvInvert=             gd.getNextNumber();  //0.05; // when FFT component is lass than this fraction of the maximal value, replace 1/z with Z
    corr_gamma=               gd.getNextNumber();  //0.2; lower the value - higher harmonics will participate in pattern frequency measurements
    corr_sigma=               gd.getNextNumber();  // 1.5; // pattern detection: high-pass filter (0.0 - none) gamma(PS)
    diffSpectrCorr=     (int) gd.getNextNumber();
    shrinkClusters=           gd.getNextNumber(); // 0.5; //  Shrink clusters by this ratio (remove lowest) after initial separation
    multiplesToTry=     (int) gd.getNextNumber();
    deviation=                gd.getNextNumber();
    deviation_steps=    (int) gd.getNextNumber();
    model_highpass=           gd.getNextNumber();
    corrRingWidth=            gd.getNextNumber();
    minCorrContrast=          gd.getNextNumber();
    correctBayerRed=          gd.getNextBoolean();
    correctBayerBlue=         gd.getNextBoolean();
    correctBayerDiagonal=     gd.getNextBoolean();
    correctBayerChecker=      gd.getNextBoolean();
    correctBayerGreen1=       gd.getNextBoolean();
    correctBayerGreen2=       gd.getNextBoolean();
    referenceComponent=(int)  gd.getNextNumber();
    equalizeGreens=       gd.getNextBoolean();

    kernelHalSizeSingleColor=(int) gd.getNextNumber();
    kernelHalSizeDualColors= (int) gd.getNextNumber();

    PSF_subpixel=            (int) gd.getNextNumber();
    RPSF_kernel_size=        (int) gd.getNextNumber();
    colorsToCorrect[0]=correctBayerGreen1;
    colorsToCorrect[1]=correctBayerRed;
    colorsToCorrect[2]=correctBayerBlue;
    colorsToCorrect[3]=correctBayerGreen2;
    colorsToCorrect[4]=correctBayerDiagonal;
    colorsToCorrect[5]=correctBayerChecker;
    for (referenceComponent=5;(referenceComponent>=0) && (!colorsToCorrect[referenceComponent]); referenceComponent--);

    MASTER_DEBUG_LEVEL=      (int) gd.getNextNumber();
    return true;
  }

  public boolean showAverageDialog() {
    GenericDialog gd = new GenericDialog("FPN parameters");
    gd.addNumericField ("Number of measurements:            ", NUMBER_OF_MEASUREMENTS, 0);
    gd.showDialog();
    if (gd.wasCanceled()) return false;
    NUMBER_OF_MEASUREMENTS= (int) gd.getNextNumber();
//    imp_camera=acquireAndDisplay(imp_camera);
    imp_camera=acquireAndAverage(imp_camera,NUMBER_OF_MEASUREMENTS);

    return true;
  }

  public boolean showSimulDialog() {
    int i;
    double scale=1000.0;
    GenericDialog gd = new GenericDialog("Simulated pattern parameters");
    gd.addNumericField ("Pattern type (0-linear, 1 - curved, ...)", simul_patternType, 0);
    gd.addNumericField ("Pattern modifier (i.e. scale)      ", simul_patternModifier, 3); // 1.0
    gd.addNumericField ("Frequency 1, X component:          ", simul_freqX1, 4);
    gd.addNumericField ("Frequency 1, Y component:          ", simul_freqY1, 4);
    gd.addNumericField ("Phase 1, (radians):                ", simul_phase1, 4);
    gd.addNumericField ("Frequency 2, X component:          ", simul_freqX2, 4);
    gd.addNumericField ("Frequency 2, Y component:          ", simul_freqY2, 4);
    gd.addNumericField ("Phase 2, (radians):                ", simul_phase2, 4);
    gd.addNumericField ("Distorttion corr, Ax (1/"+((int) scale)+" pix):", scale*simul_corr[0], 4);
    gd.addNumericField ("                  Bx (1/"+((int) scale)+" pix):", scale*simul_corr[1], 4);
    gd.addNumericField ("                  Cx (1/"+((int) scale)+" pix):", scale*simul_corr[2], 4);
    gd.addNumericField ("                  Ay (1/"+((int) scale)+" pix):", scale*simul_corr[3], 4);
    gd.addNumericField ("                  By (1/"+((int) scale)+" pix):", scale*simul_corr[4], 4);
    gd.addNumericField ("                  Cy (1/"+((int) scale)+" pix):", scale*simul_corr[5], 4);
    gd.addNumericField ("                  Dx (1/"+((int) scale)+" pix):", scale*simul_corr[6], 4);
    gd.addNumericField ("                  Ex (1/"+((int) scale)+" pix):", scale*simul_corr[7], 4);
    gd.addNumericField ("                  Dy (1/"+((int) scale)+" pix):", scale*simul_corr[8], 4);
    gd.addNumericField ("                  Ey (1/"+((int) scale)+" pix):", scale*simul_corr[9], 4);

    gd.addCheckbox     ("Center pattern for combined greens", centerForG2); // true;  // Align pattern to phases for the diagonal (both greens) sub-array
    gd.addNumericField ("Subdivide pixels by:            ", simul_subdiv, 0);
    gd.addNumericField ("Photosensitive center part of each simulated pixel", simul_fill, 4); //0.5;  part of the (center) pixel area being "phptosensitive"
    gd.addCheckbox     ("Simulate hi-res patterns", sim_hiRes); // true;  // Align pattern to phases for the diagonal (both greens) sub-array
    gd.addNumericField("Debug Level:",      MASTER_DEBUG_LEVEL, 0);
    gd.showDialog();
    if (gd.wasCanceled()) return false;
    simul_patternType=(int) gd.getNextNumber();
    simul_patternModifier=  gd.getNextNumber();
    simul_freqX1=           gd.getNextNumber();
    simul_freqY1=           gd.getNextNumber();
    simul_phase1=           gd.getNextNumber();
    simul_freqX2=           gd.getNextNumber();
    simul_freqY2=           gd.getNextNumber();
    simul_phase2=           gd.getNextNumber();
    for (i=0;i<simul_corr.length;i++) simul_corr[i]= gd.getNextNumber()/scale;
    centerForG2=gd.getNextBoolean();
    simul_subdiv= (int) gd.getNextNumber();
    simul_fill=             gd.getNextNumber();
    sim_hiRes=              gd.getNextBoolean();
    MASTER_DEBUG_LEVEL=      (int) gd.getNextNumber();
    return true;
  }

  public boolean showRPSFDialog() {
    int i;
    GenericDialog gd = new GenericDialog("PSF reversal parameters");
    gd.addNumericField("Interpoalete kernels between measured ones (defined by FFTOverlap)", INTERPOLATE_SUBDIV, 0); //2 - interploate kernels between defined by FFTOverlap
    gd.addNumericField("OTF cutoff energy (used to determine bounding ellipse)",          OTF_cutoff_energy, 3); //0.6; use frequency points that have OTF_cutoff_energy of the total to determine ellipse for limiting frequency responce
    gd.addNumericField("OTF size of elliptical window relative to cluster size",          OTF_ellipse_scale, 3); //1.5;  // size of elliptical window relative to the cluster defined by OTF_cutoff_energy
    gd.addCheckbox    ("Use Gauss window (instead of Hamming) for Ellipse mask",          OTF_ellipse_gauss); //   // size of elliptical window relative to the cluster defined by OTF_cutoff_energy
    gd.addNumericField("OTF deconvolution parameter ",                                    OTF_deconvInvert, 3); //0.01; // when FFT component is less than this fraction of the maximal value, replace 1/z with Z
    gd.addNumericField("Direct PSF kernel size (side of square array to be stored) ",     PSF_kernel_size,  0); // 64 -kernel (to be stored) size (per color component)
    gd.addNumericField("PSF cutoff energy (used to determine bounding ellipse)",          PSF_cutoff_energy, 3); //0.9;  // Limit result kernel to proportional of the PSF, calculate initial cluster shape by this cutoff energy
    gd.addNumericField("Reversed PSF elliptical mask relative to the direct PSF ellipse", PSF_ellipse_scale, 3); //1.5;  // size of elliptical window to limuit reverse PSF as proportional to direct one
    gd.addNumericField("Reversed PSF mask thershold (completely ignore far pixels)",      RPSF_min_mask_threshold, 3); //0.01; // completely zero reversed kernel elements where elliptical mask is below this threshold

    gd.addNumericField("Reversed PSF sigma to radius ratio (0 to disable)",               RPSF_sigma_to_radius, 3); //0.2- variable blurring - sigma will be proportional distance from the center
    gd.addNumericField("Reversed PSF: scale variable sigma (in the center) from the uniform one",  RPSF_var_sigma_scale, 3); //=0.8; // reduce variable sigma in the center from uniuform one

    gd.addNumericField("Gaussian blur for individual colors",   DECONCV_BLUR_INDIVIDUAL, 3); //1.8
    gd.addNumericField("Gaussian blur for diagonal greens",     DECONCV_BLUR_DIAGONAL,   3); //2.0
    gd.addNumericField("Gaussian blur for checkerboard greens", DECONCV_BLUR_CHECKER,    3); //1.4

    gd.addNumericField("Reject sampliang aliases for individual colors",   DECONCV_ALIASREJ_INDIVIDUAL, 3); //1.0
    gd.addNumericField("Reject sampliang aliases for diagonal greens",     DECONCV_ALIASREJ_DIAGONAL,   3); //1.0
    gd.addNumericField("Reject sampliang aliases for checkerboard greens", DECONCV_ALIASREJ_CHECKER,    3); //1.5

    gd.addNumericField("Reversed PSF  kernel size",                                       RPSF_kernel_size, 0); //   32 - size of deconvolution kernel
    gd.addCheckbox    ("Update ImageJ status",                                            UPDATE_STATUS);
    gd.addCheckbox    ("Show PSF",                                                        SHOW_PSF); // show combined PSF kernels (per-color and/or composite, as defined in SHOW_INDIVIDUAL, SHOW_COMPOSITE)
    gd.addCheckbox    ("Show MTF",                                                        SHOW_MTF);  // calculate/show MTF (see notes to SHOW_PSF)
    gd.addCheckbox    ("Show inverted PSF kernels",                                       SHOW_INVERTED);        // show inverted kernels (unfiltered), same notes
    gd.addCheckbox    ("Show filtered inverted PSF kernels",                              SHOW_FILTERED);        // filter and show inverted kernels
    gd.addCheckbox    ("Show inverted PSF with reduced sampling aliases pattern",         SHOW_REDUCED_ALIASES); // calculate kernels with suppressed sampling aliases patterns
    gd.addCheckbox    ("Show gaussian kernels (same centers as in filtered",              SHOW_GAUSSIANS);       // create gaussian kernels with the same centers as inverted ones (low noise, use in low details areas)
    gd.addCheckbox    ("Show gaussian kernels with reduced sampling aliases pattern",     SHOW_REDUCED_ALIASES_GAUSSIAN); // calculate kernels with suppressed sampling aliases patterns
//    gd.addCheckbox    ("Show per-component images for the PSF/MTF/...data",               SHOW_INDIVIDUAL);      // for each of the kernels above - show per-color images
//    gd.addCheckbox    ("Show a single composite images for the PSF/MTF/...data",          SHOW_COMPOSITE);       // for each of the kernels above - show single composite image


    gd.addNumericField("Debug Level:",                                                    MASTER_DEBUG_LEVEL, 0);
    gd.showDialog();
    if (gd.wasCanceled()) return false;
    INTERPOLATE_SUBDIV=      (int) gd.getNextNumber();
    OTF_cutoff_energy=             gd.getNextNumber();
    OTF_ellipse_scale=             gd.getNextNumber();
    OTF_ellipse_gauss=             gd.getNextBoolean();
    OTF_deconvInvert=              gd.getNextNumber();

    PSF_kernel_size=         (int) gd.getNextNumber();
    PSF_cutoff_energy=             gd.getNextNumber();
    PSF_ellipse_scale=             gd.getNextNumber();
    RPSF_min_mask_threshold=       gd.getNextNumber();
    RPSF_sigma_to_radius=          gd.getNextNumber();
    RPSF_var_sigma_scale=          gd.getNextNumber();

    DECONCV_BLUR_INDIVIDUAL= gd.getNextNumber();
    DECONCV_BLUR_DIAGONAL=   gd.getNextNumber();
    DECONCV_BLUR_CHECKER=    gd.getNextNumber();


    DECONCV_ALIASREJ_INDIVIDUAL= gd.getNextNumber();
    DECONCV_ALIASREJ_DIAGONAL=   gd.getNextNumber();
    DECONCV_ALIASREJ_CHECKER=    gd.getNextNumber();

    RPSF_kernel_size=1;
    for (i=(int) gd.getNextNumber(); i >1; i>>=1) RPSF_kernel_size <<=1; /* make it to be power of 2 */
    UPDATE_STATUS=                gd.getNextBoolean();
    SHOW_PSF=                     gd.getNextBoolean();
    SHOW_MTF=                     gd.getNextBoolean();
    SHOW_INVERTED=                gd.getNextBoolean();
    SHOW_FILTERED=                gd.getNextBoolean();
    SHOW_REDUCED_ALIASES=         gd.getNextBoolean();
    SHOW_GAUSSIANS=               gd.getNextBoolean();
    SHOW_REDUCED_ALIASES_GAUSSIAN=gd.getNextBoolean();
//    SHOW_INDIVIDUAL=              gd.getNextBoolean();
//    SHOW_COMPOSITE=               gd.getNextBoolean();

    MASTER_DEBUG_LEVEL=      (int) gd.getNextNumber();
    return true;
  }


  public ImagePlus acquireAndDisplay(ImagePlus imp_src) {
     ImagePlus imp=jp4_instance.openURL(imp_src);
//     if (imp!=null) measureAndDisplay(imp);
     return imp;
 }

  public ImagePlus acquireAndAverage (ImagePlus imp_src, int number) {
    int i,j;
    ImagePlus imp=imp_src;
    ImageProcessor ip=null;
    float [] averagePixels=null;
    float [] pixels=null;
    for (i=0;i<number;i++) {
      imp=jp4_instance.openURL(imp);
      if (imp==null) return null;
      if (i==0) ip = imp.getProcessor();
      pixels = (float[])ip.getPixels();
      if (i==0) averagePixels= pixels.clone();
      else for (j=0; (j < pixels.length) && (j < averagePixels.length); j++ ) averagePixels[j]+=pixels[j];
      IJ.showStatus("Averaging image: "+i);
    }
    for (j=0; j < pixels.length; j++ ) averagePixels[j]/=number; 
    ip.setPixels(averagePixels);
    return imp;
 }


  public ImagePlus selectImage(int square, boolean ALLOW_FULL_IMAGE){
    ImagePlus imp_src = WindowManager.getCurrentImage();
    if (imp_src==null){
      IJ.showMessage("Error","There are no images open\nProcess canceled");
      return null;
    }
    title_src=imp_src.getTitle();
    Roi roi_src= imp_src.getRoi();    
    if (imp_src.getType() !=ImagePlus.GRAY32 ) {
      if ((imp_src.getType() ==ImagePlus.GRAY8 ) ||
          (imp_src.getType() ==ImagePlus.GRAY16) ) {
         IJ.showStatus("Converting source image to gray 32 bits (float)");
         new ImageConverter(imp_src).convertToGray32();
      } else {
        IJ.showMessage("Error","Image should be Bayer array as a grayscale (8,16 or 32 bits)");
        return null;
      }
    }
    if (roi_src==null){
      if (ALLOW_FULL_IMAGE) {
        imp_src.setRoi(0, 0, imp_src.getWidth(), imp_src.getHeight());
        roi_src= imp_src.getRoi();    
      } else {
        IJ.showMessage("Error","No selection (line or rectangle) in the source image.\n"+
                               "You may allow processing of the full image in \"Configure\"");
        return null; /* Full image selected */
      }
    }
    Rectangle r=roi_src.getBounds();
/* align ROI to Bayer */
    r.width &= ~1;
    r.height &= ~1;
    r.x &= ~1;
    r.y &= ~1;

/* make selection square for FFT */
    if (square>0) {
      r.x+=((r.width -  square) / 2) & ~1;
      r.y+=((r.height - square) / 2) & ~1;
      r.width= square;
      r.height=square;
    }
//    selecWidth=r.width;
//    selecHeight=r.height;
    imp_src.setRoi(r);
    return imp_src;
  }
/* converts array of Bayer planes (including combined greens) into array of RGB triplets */
public double [][] mergeBayer(double [][] bayer_pixels, // array af color components arrays, some may be null
                                    int componentWidth, // component width (not the result image)
                                         int subSample) {
   if (bayer_pixels==null) return null ;
   if (subSample<2) {
     if ( (bayer_pixels.length<4) ||
          (bayer_pixels[0]==null)||
          (bayer_pixels[1]==null)||
          (bayer_pixels[2]==null)||
          (bayer_pixels[3]==null) ) return null ;
   } else {
     if ( (bayer_pixels.length<6) ||
          (bayer_pixels[1]==null)||
          (bayer_pixels[2]==null)||
          (bayer_pixels[5]==null) ) return null ; // modify to use alternative greens (diagonal, "4"?)
   }
   int componentHeight=bayer_pixels[5].length/componentWidth;
   int width=componentWidth*((subSample>1)?1:2);
   int height=componentHeight*((subSample>1)?1:2);
   double [][] doubleRGB=new double [width*height][3];
   int i,j,k,index;
   int dIndexR=subSample/2;
   int dIndexB=width*subSample/2;
  
   if (subSample<2) {
     for (i=0;i<height;i++) for (j=0;j<width; j++) for (k=0;k<3;k++) doubleRGB[i*width+j][k]=0.0;
     for (i=0;i<componentHeight;i++) {
       for (j=0;j<componentWidth;j++) {
         index=i*width+j;
         doubleRGB[(2*i  )*width+(2*j  )][1]=bayer_pixels[0][index]; // green1
         doubleRGB[(2*i+1)*width+(2*j+1)][1]=bayer_pixels[3][index]; // green1
         doubleRGB[(2*i  )*width+(2*j+1)][0]=bayer_pixels[1][index]; // red
         doubleRGB[(2*i+1)*width+(2*j  )][2]=bayer_pixels[2][index]; // blue
       }
     }
/* Add simple bi-linear demosaic here? */
   } else {
     for (i=0;i<subSample/2;i++)  for (j=0;j<width;      j++) doubleRGB[i*width+j][2]=0.0;//blue
     for (i=0;i<height;     i++)  for (j=0;j<subSample/2;j++) doubleRGB[i*width+j][0]=0.0;//red
     for (i=0;i<height;i++) {
       for (j=0;j<width;j++) {
         index=i*width+j;
         doubleRGB[index][1]=bayer_pixels[5][index]; // greens
         if ((j+dIndexR)<width ) doubleRGB[index+dIndexR][0]=bayer_pixels[1][index]; // red
         if ((i+dIndexR)<height) doubleRGB[index+dIndexB][2]=bayer_pixels[2][index]; // blue
       }
     }
   }
   return doubleRGB; 
 }

public ImagePlus showDoubleColor(double [][] doubleRGB, // array of color triplets
                                         ImagePlus imp, // ImagePlus to reuse or null for the new one 
                                             int width,
                                            int height,
                                          String title) {
  int [] pixels=doubleRGB2color(doubleRGB);
  if (pixels==null) return null;
  ImageProcessor ip;
  if (imp!=null) {
    ip=imp.getProcessor();
    ip.setPixels(pixels);
    ip.resetMinAndMax();
  } else {
    ip= new ColorProcessor(width,height);
    ip.setPixels(pixels);
    ip.resetMinAndMax();
    imp= new ImagePlus(title, ip);
    imp.show();
  }
  return imp;
}


/* convert array of double RGB triplets (in the range of 0.0..1.0) to an integer array to be set to IJ ColorProcessor */
 public int [] doubleRGB2color(double [][] doubleRGB) {
   if ((doubleRGB==null) || (doubleRGB[0]==null) || (doubleRGB[0].length<3)) return null;
   int i;
   int r,g,b;
   int size=doubleRGB.length;
   int [] intColor= new int [size];
   for (i=0; i<size; i++) {
     r= (int) Math.round(255.0*Math.max(Math.min(doubleRGB[i][0],1.0),0.0));
     g= (int) Math.round(255.0*Math.max(Math.min(doubleRGB[i][1],1.0),0.0));
     b= (int) Math.round(255.0*Math.max(Math.min(doubleRGB[i][2],1.0),0.0));
     intColor[i]=(r<<16)+(g<<8)+b;
   }
   return intColor;
 }
/* simple white and brightness/contrast balance */
public double [][] balanceDoubleRGB (double [][] doubleRGB, // pixel array of triplets
                                      double lowPerc,       // fraction of pixels to become min (lowLevel). Use negative value to preserve original zeros
                                      double lowLevel,      // level, so lowPerc of all pixels will be lower than that
                                      double highPerc,      // fraction of pixels to become max
                                      double highLevel){    // level, so highPerc of all pixels will be higher than that
    int histSize=1000;
    if ((doubleRGB==null) || (doubleRGB[0]==null) || (doubleRGB[0].length<3)) return null;
    int size=doubleRGB.length;
    double minRGB, maxRGB, k,d, di;
    int n,i;
    int []      histogram=   new int [histSize+1];
    double [][] balancedRGB= new double [size][3];
    double      highPercOriginal,lowPercOriginal;
    for (n=0;n<3;n++) {
      minRGB=doubleRGB[0][n];
      maxRGB=minRGB;
      for (i=1;i<size;i++) {
        if      (maxRGB<doubleRGB[i][n]) maxRGB=doubleRGB[i][n];
        else if (minRGB>doubleRGB[i][n]) minRGB=doubleRGB[i][n];
      }
      k=(1.0*(histSize-0.01))/(maxRGB-minRGB);
/* create histogram */
      for (i=0;i<=histSize;i++) histogram[i]=0;
      for (i=0;i<size;i++) histogram[1+ (int)(k*(doubleRGB[i][n]-minRGB))]++;
/* make it cumulative */
      for (i=1;i<=histSize;i++) histogram[i]+=histogram[i-1];
/* Find the upper level that matches highPerc */
      d=histogram[histSize]*(1.0-highPerc);
      if (d<0.0) d=0.0;
      for (i=histSize;(i>0) && (histogram[i]>d);i--);
      highPercOriginal=maxRGB;
      if (i<histSize) {
        di=i;
        if (histogram[i+1]>histogram[i]) di+= (d-histogram[i])/(histogram[i+1]-histogram[i]);
        highPercOriginal=minRGB+ ((di*(maxRGB-minRGB))/histSize);
      } // else keep highPercOriginal=maxRGB
/* Find the lower level that matches lowPerc (or use 0) if negative */
      lowPercOriginal=0.0;
      if (lowPerc>=0) {
        d=histogram[histSize]*lowPerc;
        for (i=1;(i<=histSize) && (histogram[i]<=d);i++); // histogram[0]=0;
        if (i>histSize) lowPercOriginal=maxRGB;
        else {
          di=i;
          if (histogram[i]>histogram[i-1]) di+= (d-histogram[i-1])/(histogram[i]-histogram[i-1]);
          lowPercOriginal=minRGB+ ((di*(maxRGB-minRGB))/histSize);
        }
      } // else keep lowPercOriginal=0
/* scale data into the new range */
      k=(highLevel-lowLevel)/(highPercOriginal-lowPercOriginal);
      for (i=0;i<size;i++)  balancedRGB[i][n]=lowLevel+k*(doubleRGB[i][n]-lowPercOriginal);
    }
    return balancedRGB;
}

/* Gamma correction for color images */
public double [][] gammaDoubleRGB (double [][] doubleRGB, // pixel array of triplets
                                              double gamma,
                                              double headroom) { // part of the full range reserved for saturation boost after gamma
    if ((doubleRGB==null) || (doubleRGB[0]==null) || (doubleRGB[0].length<3)) return null;
    int size=doubleRGB.length;
    int i;
    double Kb=0.114,Kr=0.299;
    double [][] gammaRGB= new double [size][3];
    double Y, YG;
    double k=1.0-headroom;
/*  Preserving same ratio of the color components after the gamma. It is probably not correct */
    for (i=0;i<size;i++) {
       Y=Kr*doubleRGB[i][0]+(1-Kr-Kb)*doubleRGB[i][1]+Kb*doubleRGB[i][2];
       YG=(Y>0)?Math.pow(Y,gamma):0.0;
       gammaRGB[i][0]=k*YG*doubleRGB[i][0]/Y;
       gammaRGB[i][1]=k*YG*doubleRGB[i][1]/Y;
       gammaRGB[i][2]=k*YG*doubleRGB[i][2]/Y;
    }
    return gammaRGB;
}





/**==========================================*/

 /* Use ROI */
  public double[][] splitBayer (ImagePlus imp, boolean equalize_greens) {
    Roi roi_src= imp.getRoi();    
    if (roi_src==null){
        imp.setRoi(0, 0, imp.getWidth(), imp.getHeight());
        roi_src= imp.getRoi();
    }
    return splitBayer (imp, roi_src.getBounds(),equalize_greens);
  }

/* Suppiy rectangle */
  public double[][] splitBayer (ImagePlus imp, Rectangle r, boolean equalize_greens) {
    ImageProcessor ip=imp.getProcessor();
    float [] pixels;
    pixels=(float[])ip.getPixels();    
    if (DEBUG_LEVEL>10) IJ.showMessage("splitBayer","r.width="+r.width+
                              "\nr.height="+r.height+
                              "\nr.x="+r.x+
                              "\nr.y="+r.y+
                              "\nlength="+pixels.length);
    int x,y,base,base_b,bv,i;
    int half_height=r.height>>1;
    int half_width=r.width>>1;
    int full_width= imp.getWidth();  // full image width
    int full_height=imp.getHeight(); // full image height
    int numColors=(half_height==half_width)?5:4;
    int pixX,pixY;
    double [][] bayer_pixels=new double[numColors][half_height * half_width];
//      base=r.width*((y<<1)+bv);
    for (y=0; y<half_height; y++) for (bv=0;bv<2;bv++){
      pixY=(y*2)+bv+r.y;
      base_b=half_width*y;
      if ((pixY<0) || (pixY>=full_height)) {
        if (bv==0) for (x=0; x<half_width; x++) {
          bayer_pixels[0][base_b]= 0.0;
          bayer_pixels[1][base_b]= 0.0;
          base_b++;
        } else  for (x=0; x<half_width; x++) {
          bayer_pixels[2][base_b]= 0.0;
          bayer_pixels[3][base_b]= 0.0;
          base_b++;
        }
      } else {
        base=full_width*((y*2)+bv+r.y)+r.x;
        pixX=r.x;
        if (bv==0) for (x=0; x<half_width; x++) {
          if ((pixX<0) || (pixX>=(full_width-1))) {
            bayer_pixels[0][base_b]= 0.0;
            bayer_pixels[1][base_b]= 0.0;
            base+=2;
          } else {
            bayer_pixels[0][base_b]= pixels[base++];
            bayer_pixels[1][base_b]= pixels[base++];
          }
          base_b++;
          pixX+=2;
        } else  for (x=0; x<half_width; x++) {
          if ((pixX<0) || (pixX>=(full_width-1))) {
            bayer_pixels[2][base_b]= 0.0;
            bayer_pixels[3][base_b]= 0.0;
            base+=2;
          } else {
            bayer_pixels[2][base_b]= pixels[base++];
            bayer_pixels[3][base_b]= pixels[base++];
          }
          base_b++;
          pixX+=2;
        }
      }
    }
    double g0=0,g3=0;
    if (equalize_greens) {
      for (i=0;i<bayer_pixels[0].length;i++) {
         g0+=bayer_pixels[0][i];
         g3+=bayer_pixels[3][i];
      }
      if ((g0>0.0) && (g3>0.0)) {
         g0=Math.sqrt(g3/g0);
         g3=1.0/g0;
         for (i=0;i<bayer_pixels[0].length;i++) {
           bayer_pixels[0][i]*=g0;
           bayer_pixels[3][i]*=g3;
         }
      }
      if (DEBUG_LEVEL>3) System.out.println("greens equalized, g0*="+g0+", g3*="+g3);
    }

    if (numColors>4) bayer_pixels[4]=combineDiagonalGreens (bayer_pixels[0], bayer_pixels[3],  half_width, half_height);
    return bayer_pixels;
  }
  
  public double [] combineDiagonalGreens (double [] green0, double []green3, int half_width, int half_height) {
    int y,x,base;
    int base_b=0;
    double [] result= new double [green0.length];
    for (y=0;y<half_height/2; y++){
      base=half_height*half_width/2+ y* (half_width+1);
      for (x=0; x<half_width/2; x++) {
        result[base_b++]=green0[base];
        base-=half_width;
        result[base_b++]=green3[base++];
      }
      base=half_height*half_width/2+ y* (half_width+1);
      for (x=0; x<half_width/2; x++) {
//System.out.println("2:y="+y+" x="+x+" base_b="+base_b+" base="+base);
        result[base_b++]=green3[base++];
        result[base_b++]=green0[base];
        base-=half_width;
      }
    }
    return result;
  }
  public void normalizeKernel(double [][] kernel) {
    int i;
    for (i=0;i<kernel.length;i++) if (kernel[i]!=null) normalizeKernel(kernel[i]);
  }

  public void normalizeKernel(double [] kernel) {
//    if (kernel==null) return null;
    int i;
    double s=0;
    for (i=0;i<kernel.length;i++) s+= kernel[i];
    s=1.0/s;
    for (i=0;i<kernel.length;i++) kernel[i]*=s;
  }

 /* merge 2-d array of kernels into large 2d-array of pixels (i.e. to be shown with showBayer()) */
  public double [][] mergeSquareKernelsToThree(double [][][][] kernels) { 
    if (kernels==null) return null;
    int tilesY=kernels.length;
    int tilesX=kernels[0].length;
    int i,j,k;
    double [][]kernel=null;
    for (i=0;(i<tilesY) && (kernel==null);i++)  for (j=0;(j<tilesX) && (kernel==null);j++)  kernel=kernels[i][j];
    if (kernel==null) return null;
    double [][]merged = new double [kernel.length][];
    int length=0;
    for (i=0;i<kernel.length;i++) if (kernel[i]!=null){
      length=kernel[i].length;
      break;
    }
    int size=(int) Math.sqrt(length);
    int outWidth= size*tilesX;
    int outHeight=size*tilesY;
    for (i=0;i<kernel.length;i++) {
      if (kernel[i]!=null) merged[i]=new double[outWidth*outHeight];
      else                 merged[i]=null;
    }
    int x,y, index;
    if (DEBUG_LEVEL>1) {
       System.out.println("mergeSquareKernelstoThree(): tilesY="+tilesY+" tilesX="+tilesX);
    }

    for (i=0; i<tilesY;i++)  for (j=0;j<tilesX;j++) for (k=0;k<merged.length;k++) if (merged[k]!=null) {

//     if (DEBUG_LEVEL>1) { System.out.println("mergeSquareKernels(): i="+i+" j="+j+" k="+k); }
      for (y=0;y<size;y++) for (x=0;x<size;x++) {
        index=((i*size+y)*outWidth)+(j*size+x);
//         if (DEBUG_LEVEL>1) {
//          System.out.println("mergeSquareKernels(): i="+i+" j="+j+" k="+k+" y="+y+" x="+x+" index="+index);
//        }
        if (kernels[i][j]==null) merged[k][index]=0.0;
        else {
               merged[k][index]=kernels[i][j][k][y*size+x];
//          if (DEBUG_LEVEL>1) System.out.println("mergeSquareKernels(): merged["+k+"]["+index+"]=kernels["+i+"]["+j+"]["+k+"]["+(y*size+x)+"]="+kernels[i][j][k][y*size+x]);
        }
      }
    }
    return merged;
  }
/*
  public void SDFA_instance.showImageStack(ImageStack stack, String title) {
      ImagePlus imp_stack = new ImagePlus(title, stack);
      imp_stack.getProcessor().resetMinAndMax();
      imp_stack.show();
  }
*/
  public ImageStack aliasScissorsStack (ImageStack imageStack,  // stack with 3 colors/slices with the image
                                                     int size, // 64 - fft size
                                                 int subpixel, // number of image pixels/ sensor pixels (each direction) == 2
                                     double debayer_threshold, // no high frequencies - use default uniform filter
                          double debayer_relative_width_green, // realtive width of the low-pass filter for greens (1.0 goes to zero at exactrly the closest alias )
                        double debayer_relative_width_redblue, // same for individual (re and blue) color components
                                         double debayer_gamma, // power function applied to the amplitudes before generating spectral masks
                                            double debayer_rz, // for green mask - rays start radius from center, relative to distance to the first alias
                                            double debayer_ra, // for green mask - rays radius around aliases, relative to distance to the first alias
                                         double debayer_sigma, // for green mask - reduce value of the far amplitudes, relative to distance to the first alias
                                         double debayer_decay, // for green mask - exponential decay of the ray value, relative to distance to the first alias
                                  double debayer_farthest_max, // for green mask - rays will be extyended from maf back to the center, max up to this distance (relative)
                                  double debayer_radius_power, // for green mask -  divide ray values by the radius to this power
                                           double mainToAlias,// relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
                                     double debayer_mask_blur, // for both masks  sigma for gaussian blur of the binary masks (<0 -do not use "scissors")
                                     boolean debayer_lo_green, //  combine alias-reject "scissors" with lopass filter for greens
                                 boolean debayer_lo_postgreen, // combine alias-reject "scissors" with lopass filter for greens applied to red/blue
                                   boolean debayer_lo_redblue, // combine alias-reject "scissors" with lopass filter for greens applied to red/blue
                                         boolean updateStatus) // update status info
  {
int yTileDebug=6;// 3; //4;
int xTileDebug=6 ; //4; //23;

    if (imageStack==null) return null;
//    String [] chnNames={"red","blue","green"};
    int imgWidth=imageStack.getWidth();
    int imgHeight=imageStack.getHeight();
    int length=imgWidth*imgHeight;
    int step=size/2;
    int tilesX=imgWidth/step-1; // horizontal number of overlapping tiles in the source image (should be expanded from the registerd one by "step" in each direction)
    int tilesY=imgHeight/step-1; // vertical number of overlapping tiles in the source image (should be expanded from the registerd one by "step" in each direction)
    int nChn=imageStack.getSize();
    int i,chn,tileX,tileY;
/* find number of the green channel - should be called "green", if none - use last */
    int greenChn=nChn-1;
    for (chn=0;chn<nChn;chn++) if (imageStack.getSliceLabel(chn+1).equals("green")){
      greenChn=chn;
      break;
    }
    float [][] outPixels=new float[nChn][length]; // same as input
    for (chn=0;chn<nChn;chn++) for (i=0;i<length;i++) outPixels[chn][i]=0.0f;
//    float [][] pixels=(float[][]) imageStack.getImageArray();
    float [][] pixels= new float[nChn][];
    for (chn=0;chn<nChn;chn++) pixels[chn]= (float[]) imageStack.getPixels(chn+1);
    
    double [][] tile=        new double[nChn][size * size ];
    double [] slidingWindow= getSlidingMask(size); // 64x64
    double [][] both_masks;
    PolarSpectrums pol_instace=new PolarSpectrums(
                 size, // size of the square array, centar is at size/2, size/2, only top half+line will be used
                 Math.PI, //2*Math.PI, // i.e. Math.PI, 2*Math.PI
                 size/2-2, // width of the polar array - should be <= size/2-2
                 0.5, //0.75, //2.0, //0.5, // maximal step in pixels on the maxRadius for 1 angular step (i.e. 0.5)
                 4 );// angular symmetry - 0- none,1 - pi corresponds to integer, 2 - pi/2 corresponds to integer, n - pi/n corresponds to integer angular step 
/* TODO: Separate those widths upstream*/
    double [][][] lopass= createAliasFilters (debayer_relative_width_green, // realtive width of the low-pass filter for greens (1.0 goes to zero at exactrly the closest alias )
                                            debayer_relative_width_redblue, // same for individual (re and blue) color components
                                                                      size, // side of the square
                                                                2* subpixel); // should be 4 now
double [] debayer_energy=null;
if (DEBUG_LEVEL>1) {
  debayer_energy=new double[tilesY*tilesX];
}

    for (tileY=0;tileY<tilesY;tileY++) {
      if (updateStatus) IJ.showStatus("Reducing sampling aliases, row "+(tileY+1)+" of "+tilesY);
      for (tileX=0;tileX<tilesX;tileX++) {
/* Read source image tiles (for each channel), swap corners and calculate FHT */
if ((tileY==yTileDebug) && (tileX==xTileDebug)) DEBUG_LEVEL=4;
else DEBUG_LEVEL=2;
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
if ((tileY==yTileDebug) && (tileX==xTileDebug)) SDFA_instance.showArrays (tile.clone(),size,size, "x"+(tileX*step)+"_y"+(tileY*step));
        for (chn=0;chn<nChn;chn++){
           fht_instance.swapQuadrants(tile[chn]);
           fht_instance.transform(tile[chn]);
        }
if ((tileY==yTileDebug) && (tileX==xTileDebug)) SDFA_instance.showArrays (tile.clone(),size,size, "tile-fht");
        both_masks= aliasScissors(tile[greenChn], // fht array for green, will be masked in-place
                                      2*subpixel, // measured kernels are sampled at subpixel higher than Bayer (subpixel/2 - than sensor) //4
                               debayer_threshold, // no high frequencies - use default uniform filter
//                          debayer_relative_width, // Debayer lopass filter width (relative to distance to the nearest alias)
                                   debayer_gamma, // power function applied to the amplitudes before generating spectral masks
                                      debayer_rz, // for green mask - rays start radius from center, relative to distance to the first alias
                                      debayer_ra, // for green mask - rays radius around aliases, relative to distance to the first alias
                                   debayer_sigma, // for green mask - reduce value of the far amplitudes, relative to distance to the first alias
                                   debayer_decay, // for green mask - exponential decay of the ray value, relative to distance to the first alias
                            debayer_farthest_max, // for green mask - rays will be extyended from maf back to the center, max up to this distance (relative)
                            debayer_radius_power, // for green mask -  divide ray values by the radius to this power
                                     mainToAlias,// relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
                               debayer_mask_blur, // for both masks  sigma for gaussian blur of the binary masks (<0 -do not use "scissors")
                                debayer_lo_green, //  combine alias-reject "scissors" with lopass filter for greens
                            debayer_lo_postgreen, // combine alias-reject "scissors" with lopass filter for greens applied to red/blue
                              debayer_lo_redblue, // combine alias-reject "scissors" with lopass filter for greens applied to red/blue
              (debayer_decay<0)?pol_instace:null,
                                          lopass,
                 ((tileY==yTileDebug) && (tileX==xTileDebug))?4:1);
//                                               1); // internal debug level ((DEBUG_LEVEL>2) && (yTile==yTile0) && (xTile==xTile0))?3:1;
if ((tileY==yTileDebug) && (tileX==xTileDebug)) {
   SDFA_instance.showArrays (tile.clone(),size,size, "A00");
   SDFA_instance.showArrays (both_masks.clone(),size,size, "masks");
   
}
if (DEBUG_LEVEL>1) {
  debayer_energy[tileY*tilesX+tileX]=DOUBLE_DEBUG_RESULT;
//      System.out.println("debayer_energy["+tileY+"]["+tileX+"]="+DOUBLE_DEBUG_RESULT);

//  DOUBLE_DEBUG_RESULT= maxAmpInRing (green_amp); 
}


        for (chn=0;chn<nChn;chn++) {
//          if (chn!=greenChn) fht_instance.multiply(tile[chn],both_masks[1],false); // green is already multiplied
          tile[chn]=fht_instance.multiply(tile[chn],both_masks[(chn==greenChn)?0:1],false);
          fht_instance.inverseTransform(tile[chn]);
          fht_instance.swapQuadrants(tile[chn]);
/* accumulate result */
          accumulateSquareTile(outPixels[chn], //  float pixels array to accumulate tile
                                    tile[chn], // data to accumulate to the pixels array
                                     imgWidth, // width of pixels array
                                   tileX*step, // left corner X
                                   tileY*step); // top corner Y
        }
if ((tileY==yTileDebug) && (tileX==xTileDebug)) SDFA_instance.showArrays (tile.clone(),size,size, "B00");

      }
    }
/* prepare result stack to return */
    ImageStack outStack=new ImageStack(imgWidth,imgHeight);
    for (chn=0;chn<nChn;chn++) {
      outStack.addSlice(imageStack.getSliceLabel(chn+1), outPixels[chn]);
    }
/*
if (DEBUG_LEVEL>1) {
   SDFA_instance.showArrays (debayer_energy,tilesX,tilesY, "Debayer-Energy");
}
*/


    return outStack;
  }




  /* convolve image stack with the kernel stack using FHT. kernels ahould be (size/2)*(size/2) - currently 64x64, then image will be split into same 
      (size/2)*(size/2) overlappig by step=size/4 segments. Both are zero-padded to size x size, so after convolution the result will not roll over, and
      processed 128x128 result arrays are accumulated in the output stack.
      The input image should be properly extended by size/4 in each direction (and so the kernel arrays should match it) - that would minimize border effects.*/
   public ImageStack convolveStackWithKarnelStack (ImageStack imageStack,  // stack with 3 colors/slices with the image
                                                 ImageStack kernelStack, // stack with 3 colors/slices convolution kernels
                                                               int size, // 128 - fft size, kernel size should be size/2 
                                                   boolean updateStatus) // update status info
  {
    if ((imageStack==null) || (kernelStack==null)) return null;
    int imgWidth=imageStack.getWidth();
    int imgHeight=imageStack.getHeight();
    int length=imgWidth*imgHeight;
    int step=size/4;
    int kernelSize=size/2;
    int tilesX=imgWidth/step-1; // horizontal number of overlapping tiles in the source image (should be expanded from the registerd one by "step" in each direction)
    int tilesY=imgHeight/step-1; // vertical number of overlapping tiles in the source image (should be expanded from the registerd one by "step" in each direction)
    int kernelWidth=kernelStack.getWidth();
    int kernelNumHor=kernelWidth/(size/2);

    int nChn=imageStack.getSize();
    float [][] outPixels=new float[nChn][length]; // same as input
    int i,chn,tileX,tileY;
    for (chn=0;chn<nChn;chn++) for (i=0;i<length;i++) outPixels[chn][i]=0.0f;
    float [] pixels;
    float [] kernelPixels;
    double [] kernel=       new double[kernelSize*kernelSize];
    double [] inTile=       new double[kernelSize*kernelSize];
    double [] outTile=      new double[size    * size   ];
    double [] doubleKernel= new double[size    * size   ];
    double [] slidingWindow=getSlidingMask(kernelSize); // 64x64

    for (chn=0;chn<nChn;chn++) {
      pixels=      (float[]) imageStack.getPixels(chn+1);
      kernelPixels=(float[]) kernelStack.getPixels(chn+1);
      for (tileY=0;tileY<tilesY;tileY++) {
        if (updateStatus) IJ.showStatus("Convolving image with kernels, channel "+(chn+1)+" of "+nChn+", row "+(tileY+1)+" of "+tilesY);
        for (tileX=0;tileX<tilesX;tileX++) {
/* Read source image tile */
           extractSquareTile( pixels, // source pixel array,
                              inTile, // will be filled, should have correct size before call
                       slidingWindow, // window (same size as the kernel)
                            imgWidth, // width of pixels array
                          tileX*step, // left corner X
                          tileY*step); // top corner Y
/* zero pad twice the original size*/
          outTile=extendFFTInputTo (inTile, size);
//if ((tileY==tilesY/2) && (tileX==tilesX/2))  SDFA_instance.showArrays(outTile,size,size, "input-"+chn);
/* FHT transform of the source image data*/
          fht_instance.swapQuadrants(outTile);
          fht_instance.transform(    outTile);
/* read convolution kernel */
          extractOneKernel(kernelPixels, //  array of combined square kernels, each 
                                  kernel, // will be filled, should have correct size before call
                            kernelNumHor, // number of kernels in a row
//                        tileX*kernelSize, // horizontal number of kernel to extract
//                        tileY*kernelSize); // vertical number of kernel to extract
                                    tileX, // horizontal number of kernel to extract
                                    tileY); // vertical number of kernel to extract
/* zero pad twice the original size*/
          doubleKernel=extendFFTInputTo (kernel, size);
//if ((tileY==tilesY/2) && (tileX==tilesX/2))  SDFA_instance.showArrays(doubleKernel,size,size, "doubleKernel-"+chn);
/* FHT transform of the kernel */
          fht_instance.swapQuadrants(doubleKernel);
          fht_instance.transform(    doubleKernel);
/* multiply in frequnecy domain */
          outTile=     fht_instance.multiply(outTile, doubleKernel, false);
/* FHT inverse transform of the product - back to space domain */
          fht_instance.inverseTransform(outTile);
          fht_instance.swapQuadrants(outTile);
/* accumulate result */
//if ((tileY==tilesY/2) && (tileX==tilesX/2))  SDFA_instance.showArrays(outTile,size,size, "out-"+chn);
          accumulateSquareTile(outPixels[chn], //  float pixels array to accumulate tile
                                      outTile, // data to accumulate to the pixels array
                                     imgWidth, // width of pixels array
                               (tileX-1)*step, // left corner X
                               (tileY-1)*step); // top corner Y
        }
      }
    }
/* prepare result stack to return */
    ImageStack outStack=new ImageStack(imgWidth,imgHeight);
    for (chn=0;chn<nChn;chn++) {
      outStack.addSlice(imageStack.getSliceLabel(chn+1), outPixels[chn]);
    }
    return outStack;
  }

/* Convert source Bayer pattern (GR/BG) image to higher resolution, add margins by duplicating pattern around  */
  public ImageStack  bayerToStack(ImagePlus imp, // source bayer image, linearized, 32-bit (float))
                                 int oversample, // multiple samples per pixel in each direction (2)
                                    int addLeft, // add this number of scan lines above the image (reducing border effects)
                                     int addTop, // 
                                   int addRight, // 
                                  int addBottom){

    if (imp==null) return null;
    String [] chnNames={"red","blue","green"};
    int nChn=chnNames.length;
    ImageProcessor ip=imp.getProcessor();
    int inWidth=imp.getWidth();
    int inHeight=imp.getHeight();
    int outHeight=inHeight*oversample+addTop+addBottom;
    int outWidth=inWidth*oversample+addLeft+addRight;
    int outLength=outWidth*outHeight;

    float [][] outPixels=new float[nChn][outLength];
    float [] pixels = (float[]) ip.getPixels();
    int chn,y,x,i,index;
    int bayerPeriod=2*oversample;
    int ovrWidth= inWidth*oversample;
    int ovrHeight=inHeight*oversample;
    for (chn=0;chn<nChn;chn++) for (i=0;i<outPixels[chn].length;i++) outPixels[chn][i]=0.0f;
/* Can be optimized - now it calculate input address fro all those 0-es */
    for (index=0; index<outLength; index++) {
      y=(index / outWidth)-addTop;
      x=(index % outWidth)-addLeft;
      if (y<0) y= (bayerPeriod-((-y) % bayerPeriod))%bayerPeriod;
      else if (y>=ovrHeight) y= ovrHeight-bayerPeriod +((y-ovrHeight) % bayerPeriod);
      if (x<0) x= (bayerPeriod-((-x) % bayerPeriod))%bayerPeriod;
      else  if (x>=ovrWidth) x= ovrWidth-bayerPeriod +((x-ovrWidth) % bayerPeriod);
      if (((y% oversample)==0) && ((x% oversample)==0)) {
        x/=oversample;
        y/=oversample;
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


/* interpolate kernels minimizing memory image - use directly the image stack (32-bit, float) with kernels.
    Add kernels around by either replication or extrapolation to compensate for "margins" in the original; kernels */
  public ImageStack interpolateKernelStack(ImageStack kernelStack, // Image stack, each slice consists of square kernels of one channel
                                                         int size, // size of each kernel (should be square)
                                                       int subdiv, // number of subdivisions form input to output
                                                      int addLeft, // add this number of kernel columns to the output on the left of the existent/interpolated
                                                       int addTop, // add this number of kernel rows to the output above the existent/interpolated
                                                     int addRight, // add this number of kernel columns to the output on the right of the existent/interpolated
                                                    int addBottom, // add this number of kernel rows to the output below the existent/interpolated
                                               double extrapolate,  // 0 - duplicate, 1.0 - extrapolate outside of the known kernels
                                             boolean updateStatus) // update status info
  {
    if (kernelStack==null) return null;
    int inWidth=kernelStack.getWidth();
    int inHeight=kernelStack.getHeight();
    int inTilesX=inWidth/size;
    int inTilesY=inHeight/size;
    int outTilesX= (inTilesX-1)*subdiv +1 + addLeft + addRight;
    int outTilesY= (inTilesY-1)*subdiv +1 + addTop + addBottom;
    int nChn=kernelStack.getSize();
    float [][] outPixels=new float[nChn][outTilesX*size*outTilesY*size];
    float [] pixels;
    int i,j,chn;
    int xTile0=(addLeft>0)?-1:0;
    int xTile1=inTilesX+((addRight>0)?0:-1);
    int yTile0=(addTop>0)?-1:0;
    int yTile1=inTilesY+((addBottom>0)?0:-1);
    int tileY,tileX; //,subTileY,subTileX;

    int tileWidth, tileHeight; // for inner cells (subdiv+1)*(subdiv+1), for outer includes exte row/column fro extrapolation
//    int maxTileWidth= Math.max(subdiv,1+Math.max(addRight,addLeft));
//    int maxTileHeight=Math.max(subdiv,1+Math.max(addBottom,addTop));
    boolean lastColumn=false;  //last column - inverse convert and copy the last column of rectangleFHT to the result array
    boolean lastRow=false;     //last row - interpolate, inverse convert and copy the last row of rectangleFHT to the result array

    double [] pointsVert;
    double [] pointsHor;
    double [][] fhtLine;
    double extraScale=extrapolate/subdiv;
    int [] outTopLeft=new int [2]; // top left kernel in the output array
    int [] inTopLeft=new int [2]; // top left kernel in the input array
    double [][] firstFHTColumn=null;
    double [][] secondFHTColumn=null;
    double [][][] cornerFHT=new double[2][2][size*size]; //[y][x][pixel] 
    double [] swapArray=null;

    for (chn=0;chn<nChn;chn++) {
      pixels=(float[]) kernelStack.getPixels(chn+1);
      for (tileY=yTile0;tileY<yTile1;tileY++) {
        if (updateStatus) IJ.showStatus("Interpolating kernels, channel "+kernelStack.getSliceLabel(chn+1)+", row "+(tileY-yTile0+1)+" of "+(yTile1-yTile0));
        lastRow=(tileY==(yTile1-1));
        if (tileY<0) {
          inTopLeft[1]=0;
          tileHeight=addTop;
          outTopLeft[1]=0;
          pointsVert=new double[tileHeight];
          for (i=0;i<tileHeight;i++)  pointsVert[i]=(i-tileHeight)*extraScale; // negative values
        } else if (tileY>=(inTilesY-1)){
          inTopLeft[1]=tileY-1;
          tileHeight=addBottom+1; // always last row, if got here at all (addBottom>0)
          outTopLeft[1]=addTop+subdiv*tileY;
          pointsVert=new double[tileHeight];
          for (i=0;i<tileHeight;i++)  pointsVert[i]=1.0+i*extraScale;
        } else {
          inTopLeft[1]=tileY;
          tileHeight=subdiv+ (lastRow?1:0); // last tile row includes bottom outpout kernel row
          outTopLeft[1]=addTop+subdiv*tileY;
          pointsVert=new double[tileHeight];
          for (i=0;i<tileHeight;i++) pointsVert[i]=(1.0*i)/tileHeight;
        }
        firstFHTColumn=null; // invalidate
        secondFHTColumn=null; // invalidate
        for (tileX=xTile0;tileX<xTile1;tileX++) {
    if (DEBUG_LEVEL>2)  System.out.println(" interpolateKernelStack(): chn="+chn+" tileY="+tileY+" tileX="+tileX);

          lastColumn=(tileX==(xTile1-1));
          if (tileX<0) {
            inTopLeft[0]=0;
            tileWidth=addLeft;
            outTopLeft[0]=0;
            pointsHor=new double[tileWidth];
            for (i=0;i<tileWidth;i++)  pointsHor[i]=(i-tileWidth)*extraScale; // negative values
          } else if (tileX>=(inTilesX-1)){
            inTopLeft[0]=tileX-1;
            tileWidth=addRight+1; // always last columnw, if got here at all (addRight>0)
            outTopLeft[0]=addLeft+subdiv*tileX;
            pointsHor=new double[tileWidth];
            for (i=0;i<tileWidth;i++)  pointsHor[i]=1.0+ i*extraScale;
            // else keep both firstFHTColumn and secondFHTColumn
    if (DEBUG_LEVEL>2)  System.out.println("last column: tileX="+tileX);
          } else {
            inTopLeft[0]=tileX;
            tileWidth=subdiv+ (lastColumn?1:0); // last tile column includes rightmost outpout kernel column
            outTopLeft[0]=addLeft+subdiv*tileX;
            pointsHor=new double[tileWidth];
            for (i=1;i<tileWidth;i++)  pointsHor[i]=(1.0*i)/tileWidth;
//    if (DEBUG_LEVEL>2)  System.out.println("else: tileX="+tileX);
            if (tileX!=0) {
              firstFHTColumn=secondFHTColumn;
              secondFHTColumn=null; // invalidate
//    if (DEBUG_LEVEL>2)  System.out.println(" secondFHTColumn==null");
/* swap columns, the new second one will be just reused */
              swapArray=cornerFHT[0][0];
              cornerFHT[0][0]=cornerFHT[0][1];
              cornerFHT[0][1]=swapArray;
              swapArray=cornerFHT[1][0];
              cornerFHT[1][0]=cornerFHT[1][1];
              cornerFHT[1][1]=swapArray;

            } // else keep both firstFHTColumn and secondFHTColumn
          }
    if (DEBUG_LEVEL>2)  System.out.println(" interpolateKernelStack(): tileHeight="+tileHeight+" tileWidth="+tileWidth+" inTopLeft[0]="+inTopLeft[0]+" inTopLeft[1]="+inTopLeft[1]+
                                                                                             " outTopLeft[0]="+outTopLeft[0]+" outTopLeft[1]="+outTopLeft[1]);

          if (firstFHTColumn==null) { /* First colum needs to be input and calculated*/
             extractOneKernel(          pixels, //  array of combined square kernels, each 
                               cornerFHT[0][0], // will be filled, should have correct size before call
                                      inTilesX, // number of kernels in a row
                                  inTopLeft[0], // horizontal number of kernel to extract
                                  inTopLeft[1]); // vertical number of kernel to extract
             extractOneKernel(          pixels, //  array of combined square kernels, each 
                               cornerFHT[1][0], // will be filled, should have correct size before call
                                      inTilesX, // number of kernels in a row
                                  inTopLeft[0], // horizontal number of kernel to extract
                                inTopLeft[1]+1); // vertical number of kernel to extract
/* convert to frequency domain */
             fht_instance.swapQuadrants(cornerFHT[0][0]);
             fht_instance.transform(    cornerFHT[0][0]);
             fht_instance.swapQuadrants(cornerFHT[1][0]);
             fht_instance.transform(    cornerFHT[1][0]);
/* inter/extrapolate the column */
             firstFHTColumn=interpolateFHT (cornerFHT[0][0],    // first FHT array
                                            cornerFHT[1][0],    // second FHT array
                                                 pointsVert,    // array of interpolation points - 0.0 - fht0, 1.0 - fht1
                                                     false);   // OK not to clone, so corners will be referenced?
    if (DEBUG_LEVEL>2)  System.out.println(" firstFHTColumn.length="+firstFHTColumn.length);
          }
          if (secondFHTColumn==null) { /* Last colum needs to be input and calculated*/
             extractOneKernel(          pixels, //  array of combined square kernels, each 
                               cornerFHT[0][1], // will be filled, should have correct size before call
                                      inTilesX, // number of kernels in a row
                                inTopLeft[0]+1, // horizontal number of kernel to extract
                                  inTopLeft[1]); // vertical number of kernel to extract
             extractOneKernel(          pixels, //  array of combined square kernels, each 
                               cornerFHT[1][1], // will be filled, should have correct size before call
                                      inTilesX, // number of kernels in a row
                                inTopLeft[0]+1, // horizontal number of kernel to extract
                                inTopLeft[1]+1); // vertical number of kernel to extract
/* convert to frequency domain */
             fht_instance.swapQuadrants(cornerFHT[0][1]);
             fht_instance.transform(    cornerFHT[0][1]);
             fht_instance.swapQuadrants(cornerFHT[1][1]);
             fht_instance.transform(    cornerFHT[1][1]);
/* inter/extrapolate the column */
             secondFHTColumn=interpolateFHT (cornerFHT[0][1],    // first FHT array
                                             cornerFHT[1][1],    // second FHT array
                                                  pointsVert,    // array of interpolation points - 0.0 - fht0, 1.0 - fht1
                                                       false);   // OK not to clone, so corners will be referenced?

    if (DEBUG_LEVEL>2)  {
        System.out.println(" secondFHTColumn.length="+secondFHTColumn.length);
       for (i=0;i<pointsVert.length;i++) System.out.println(""+pointsVert[i]);
       System.out.println("");
    }
          }
/* interpolate horizontally */
/* TODO: calculate top-left corner in output array */
/*
     if ((DEBUG_LEVEL>1) &&(tileY==0)) {
        SDFA_instance.showArrays(firstFHTColumn,size,size, "firstFHTColumn");
        SDFA_instance.showArrays(secondFHTColumn,size,size, "secondFHTColumn");
        DEBUG_LEVEL=4;
        return null;
     }
*/
          for (i=0;i<tileHeight;i++) {
    if (DEBUG_LEVEL>2)  System.out.print("i="+i);

             fhtLine=interpolateFHT ( firstFHTColumn[i],    // first FHT array
                                     secondFHTColumn[i],    // second FHT array
                                              pointsHor,    // array of interpolation points - 0.0 - fht0, 1.0 - fht1
                                                   true); //clone ends
    if (DEBUG_LEVEL>2)  System.out.print(": ");
             for (j=0;j<tileWidth;j++) {
    if (DEBUG_LEVEL>2)  System.out.print(j);
               fht_instance.inverseTransform(fhtLine[j]);
               fht_instance.swapQuadrants   (fhtLine[j]);
               storeOneKernel(           outPixels[chn], // float [] array of combined square kernels - will be filled
                                             fhtLine[j], // square kernel to store
                                              outTilesX, // number of kernels in a row
                                        outTopLeft[0]+j, // horizontal number of kernel to store
                                       outTopLeft[1]+i); // vertical number of kernel to store
             }
    if (DEBUG_LEVEL>2)  System.out.println("");

          }
        }
      }
    }    
/* prepare result stack to return */
    ImageStack outStack=new ImageStack(outTilesX*size,outTilesY*size);
    for (chn=0;chn<nChn;chn++) {
      outStack.addSlice(kernelStack.getSliceLabel(chn+1), outPixels[chn]);
    }
    return outStack;
  }

   public ImageStack  reversePSFKernelStack(ImageStack PSFStack, // stack of 3 32-bit (float) images, made of square kernels
                                                     int dSize, // size (side of square) of direct PSF kernel
                                                     int rSize, // size (side of square) of reverse PSF kernel
                                            double invertRange,    // deconvInvert
                                      double otf_cutoff_energy,  // OTF_cutoff_energy
                                      double otf_ellipse_scale,  // ellipse mask size relative to the cluster
                                     boolean otf_ellipse_gauss,  // use Gauss instead of Hamming for ellipse mask
                                      double psf_cutoff_energy,  // OTF_cutoff_energy
                                      double psf_ellipse_scale,  // ellipse mask size relative to the cluster
                                double rpsf_min_mask_threshold,  // zero output element if elliptical Gauss mask is below this threshold
// Optional variable-sigma blurring parameters
                                              double [] sigmas, // array of sigmas in the center, matching stacks sequence. Null if no blurring is needed
                                            double sigma_scale,  // scale sigma in the center when using variable sigma
                                          double sigmaToRadius,  // sigma-to-radius ratio (0.0 to disable variable blur)
                                          boolean updateStatus){  // update status info
    if (PSFStack==null) return null;
    int inWidth=PSFStack.getWidth();
    int inHeight=PSFStack.getHeight();
    int tilesX=inWidth/dSize;
    int tilesY=inHeight/dSize;
    int nChn=PSFStack.getSize();
    float [][] outPixels=new float[nChn][tilesX*rSize*tilesY*rSize];
    float [] pixels;
    double [] kernel= new double[dSize*dSize];
    double [] rKernel=new double[rSize*rSize];

    int  [][]selection;
    double [] ellipse_coeff;
    double [] variableSigmas;

    int chn,tileY,tileX;
    for (chn=0;chn<nChn;chn++) {
      pixels=(float[]) PSFStack.getPixels(chn+1);
      for (tileY=0;tileY<tilesY;tileY++) {
        if (updateStatus) IJ.showStatus("Reversing PSF, channel "+(chn+1)+" of "+nChn+", row "+(tileY+1)+" of "+tilesY);
        for (tileX=0;tileX<tilesX;tileX++) {
          extractOneKernel( pixels, //  array of combined square kernels, each 
                            kernel, // will be filled, should have correct size before call
                            tilesX, // number of kernels in a row
                             tileX, // horizontal number of kernel to extract
                             tileY); // vertical number of kernel to extract
          rKernel=resizeForFFT(kernel,rSize);
/* reverse PSF kernel */
/* TODO: convert cleanupAndReversePSF() to double FHT*/
          rKernel= cleanupAndReversePSF (rKernel,  // input pixels
                                     invertRange,    // deconvInvert
                               otf_cutoff_energy,  // OTF_cutoff_energy
                               otf_ellipse_scale,  // ellipse mask size relative to the cluster
                               otf_ellipse_gauss,  // use Gauss instead of Hamming for ellipse mask
                                               1,  // decimate frequency to match Bayer component pixels pitch
                                           false,  // fold high frequency into low, when false - use Hamming to cut off high frequencies
                                              ""); // just for the plot names
/* Find direct kernel approximation ellipse, increase it, mirror center around 0,0 and use it as a mask for the reversed kernel */
          selection=    findClusterOnPSF(kernel,  psf_cutoff_energy, "");
          ellipse_coeff=findEllipseOnPSF(kernel,  selection, ""); // coefficients for direct PSF, for rPSF [0] and [1] need to be opposite size
          rKernel= maskReversePSFKernel(rKernel, // reversed psf, square array
                                  ellipse_coeff, // ellipse coefficients from _direct_ kernel
                              psf_ellipse_scale,
                        rpsf_min_mask_threshold); // zero output element if elliptical Gauss mask is below this threshold

          normalizeKernel(rKernel); // in-place
/* Apply varuable blur to inversed kernel, using (and reversing sign) the center X,Y from the direct kernel */
          if (sigmas!=null) {
              variableSigmas= createSigmasFromCenter(rSize, // side of square
                                             sigmaToRadius, // variable blurring - sigma will be proportional distance from the center
                                   sigmas[chn]*sigma_scale, //blurring in the center sigma(r)=sqrt((sigma_to_radius*r)^2)+center_sigma^2)
                                         -ellipse_coeff[0], // coordinates of the center (0:0 - size/2: size/2)
                                         -ellipse_coeff[1]);
              rKernel=variableGaussBlurr(          rKernel, // input square pixel array, preferrably having many exact zeros (they will be skipped)
                                            variableSigmas, // array of sigmas to be used for each pixel, matches pixels[]
                                                       3.5, // drop calculatin if farther then nSigma
                                                         0, // int WOICenterX, // window of interest in pixels[] array - do not generate data outside it
                                                         0, // int WOICenterY, // 
                                                     rSize, //int WOIWidth, reduce later
                                                     rSize); //int WOIHeight)

          }
          storeOneKernel( outPixels[chn], // float [] array of combined square kernels - will be filled
                                 rKernel, // square kernel to store
                                  tilesX, // number of kernels in a row
                                  tileX, // horizontal number of kernel to store
                                  tileY); // vertical number of kernel to store
        }
      } 
    }
/* prepare result stack to return */
    ImageStack outStack=new ImageStack(tilesX*rSize,tilesY*rSize);
    for (chn=0;chn<nChn;chn++) {
      outStack.addSlice(PSFStack.getSliceLabel(chn+1), outPixels[chn]);
    }
    return outStack;
}

 public ImageStack  generateGaussianStackFromDirect(ImageStack PSFStack, // stack of 3 32-bit (float) images, made of square kernels
                                                              int dSize, // size (side of square) of direct PSF kernel
                                                              int size, // size (side of square) of the Gaussian kernel
                                               double psf_cutoff_energy,  // OTF_cutoff_energy
                                                       double [] sigmas, // array of sigmas in the center, matching stacks sequence. Null if no blurring is needed
                                          boolean updateStatus){  // update status info
    if (PSFStack==null) return null;
    int inWidth=PSFStack.getWidth();
    int inHeight=PSFStack.getHeight();
    int tilesX=inWidth/dSize;
    int tilesY=inHeight/dSize;
    int nChn=PSFStack.getSize();
    float [][] outPixels=new float[nChn][tilesX*size*tilesY*size];
    float [] pixels;
    double [] kernel= new double[dSize*dSize];
    int  [][]selection;
    double [] ellipse_coeff;
    int chn,tileY,tileX;
    double nSigma2;
    double x,y,x2,y2;
    int       length=   size*size;
    double [] gaussian= new double[length];
    double [] gaussianX=new double[size];
    double [] gaussianY=new double[size];
    double minsigma=0.1;
    double k,d,sum;
    int i,j;
    for (chn=0;chn<nChn;chn++) {
      pixels=(float[]) PSFStack.getPixels(chn+1);
      for (tileY=0;tileY<tilesY;tileY++) {
        if (updateStatus) IJ.showStatus("Generating Gaussians, channel "+(chn+1)+" of "+nChn+", row "+(tileY+1)+" of "+tilesY);
        for (tileX=0;tileX<tilesX;tileX++) {
          extractOneKernel( pixels, //  array of combined square kernels, each 
                            kernel, // will be filled, should have correct size before call
                            tilesX, // number of kernels in a row
                             tileX, // horizontal number of kernel to extract
                             tileY); // vertical number of kernel to extract
/* Find direct kernel approximation ellipse, mirror center around 0,0 */
          selection=    findClusterOnPSF(kernel,  psf_cutoff_energy, "");
          ellipse_coeff=findEllipseOnPSF(kernel,  selection, ""); // coefficients for direct PSF, for rPSF [0] and [1] need to be opposite size
          nSigma2=(4*sigmas[chn])*(4*sigmas[chn]);
          k=(sigmas[chn]<minsigma)?(0.5/(minsigma*minsigma)):(0.5/(sigmas[chn]*sigmas[chn]));
          for ( i=0; i<size;i++) {
            x=i-size/2+ellipse_coeff[0];
            x2=x*x;
            if (x2>nSigma2) gaussianX[i]=0.0;
            else gaussianX[i]=Math.exp(-x2*k);
            y=i-size/2+ellipse_coeff[1];
            y2=y*y;
            if (y2>nSigma2) gaussianY[i]=0.0;
            else gaussianY[i]=Math.exp(-y2*k);
          }
          sum=0.0;
          for ( i=0; i<size;i++) for (j=0;j<size;j++) {
            d=gaussianX[j]*gaussianY[i];
            sum+=d;
            gaussian[i*size+j]=d;
          }
          k=1.0/sum;
          for (i=0;i<length;i++) gaussian[i]*=k;
          storeOneKernel( outPixels[chn], // float [] array of combined square kernels - will be filled
                                gaussian, // square kernel to store
                                  tilesX, // number of kernels in a row
                                   tileX, // horizontal number of kernel to store
                                   tileY); // vertical number of kernel to store
        }
      } 
    }
/* prepare result stack to return */
    ImageStack outStack=new ImageStack(tilesX*size,tilesY*size);
    for (chn=0;chn<nChn;chn++) {
      outStack.addSlice(PSFStack.getSliceLabel(chn+1), outPixels[chn]);
    }
    return outStack;
}



/* Used in interpolateKernelStack() */  
  void storeOneKernel(float [] pixels, // float [] array of combined square kernels - will be filled
                        double [] kernel, // square kernel to store
                              int numHor, // number of kernels in a row
                               int xTile, // horizontal number of kernel to store
                               int yTile) { // vertical number of kernel to store
    int length=kernel.length;
    int size=(int) Math.sqrt(length);
    int i,j;
    int pixelsWidth=numHor*size;
    int base=(yTile*pixelsWidth+xTile)*size;
    for (i=0;i<size;i++) for (j=0;j<size;j++) pixels[base+i*pixelsWidth+j]= (float) kernel[i*size+j];
  }

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
/* accumulate square tile to the pixel array (tile may extend beyond the array, will be cropped) */
  void accumulateSquareTile(float [] pixels, //  float pixels array to accumulate tile
                             double [] tile, // data to accumulate to the pixels array
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
         if ((x>=0) && (x<width)) pixels[y*width+x]+=tile [index];
         index++;
        }
      }
    }
  }






  public ImageStack mergeKernelsToStack(double [][][][] kernels) { 
    if (kernels==null) return null;
    int tilesY=kernels.length;
    int tilesX=kernels[0].length;
    int i,j,k,nChn, chn,x,y,index;
    double [][]kernel=null;
    for (i=0;(i<tilesY) && (kernel==null);i++)  for (j=0;(j<tilesX) && (kernel==null);j++)  kernel=kernels[i][j];
    if (kernel==null) return null;
    int length=0;
    for (i=0;i<kernel.length;i++) if (kernel[i]!=null){
      length=kernel[i].length;
      break;
    }
    int [] channels = new int [kernel.length];
    nChn=0;
    for (i=0;i<kernel.length;i++) if (kernel[i]!=null)  if (nChn<channels.length) channels[nChn++]=i;
    int size=(int) Math.sqrt(length);
    int outWidth= size*tilesX;
    int outHeight=size*tilesY;

    ImageStack stack=new ImageStack(outWidth,outHeight);
    float [] fpixels;
    for (chn=0;chn<nChn;chn++) {
      fpixels= new float [outWidth*outHeight];
      k=channels[chn];
      for (i=0; i<tilesY;i++)  for (j=0;j<tilesX;j++) {
        for (y=0;y<size;y++) for (x=0;x<size;x++) {
          index=((i*size+y)*outWidth)+(j*size+x);
          if (kernels[i][j]==null) fpixels[index]=0.0f;
          else {
            fpixels[index]= (float) kernels[i][j][k][y*size+x];
          }
        }
      }
      stack.addSlice("channel"+k, fpixels);
    }
    return stack;
  }



  public double [][][][] splitSquareKernelsFromStack(ImageStack stack, // Image stack, each slice consists of square kernels of one channel
                                                             int size, // size of each kernel (should be square)
                                                       int [] channels  // four-element array of color channel numbers (usually {1,2,5,-1})
                                                                  ) { 
    if (stack==null) return null;
    int width=stack.getWidth();
    int height=stack.getHeight();
    int tilesX=width/size;
    int tilesY=height/size;


    int nChn=-1;
    int i,j,k, chn;
    for (chn=0;chn<channels.length;chn++) if (nChn<channels[chn]) nChn = channels[chn];
    nChn++;
    if (nChn==0) return null; // no channels defined
    double [][][][] kernels=new double [tilesY][tilesX][][];
    double [][]kernel=new double [nChn][];
    for (chn=0;chn<nChn;chn++) kernel[chn]=null;
    for (i=0;i<channels.length;i++) if (channels[i]>=0) kernel[channels[i]]=new double [size*size];
//    for (i=0; i<tilesY;i++)  for (j=0;j<tilesX;j++) kernels[i][j]=null;

    boolean nonzero=false;
    int x,y, index, l;
    double d;
    float [][] fpixels= new float[nChn][];

    if (DEBUG_LEVEL>1)  System.out.println("splitSquareKernelsFromStack(): nChn="+nChn);
    for  (chn=0;chn<channels.length;chn++) if (channels[chn]>=0){
//      if (DEBUG_LEVEL>1)  System.out.println("splitSquareKernelsFromStack(): chn="+chn);
      fpixels[chn]=(float []) stack.getPixels(chn+1);
    }
    for (i=0; i<tilesY;i++)  for (j=0;j<tilesX;j++) {
      kernels[i][j]=null;
      for  (chn=0;chn<channels.length;chn++) if (channels[chn]>=0){
        k=channels[chn];
        nonzero=false;
//         if (DEBUG_LEVEL>1) { System.out.println("mergeSquareKernels(): i="+i+" j="+j+" k="+k); }
/* extract one set of kernels */
        for (y=0;y<size;y++) for (x=0;x<size;x++) {
          index=((i*size+y)*width)+(j*size+x);
//             if (DEBUG_LEVEL>1) {
//              System.out.println("mergeSquareKernels(): i="+i+" j="+j+" k="+k+" y="+y+" x="+x+" index="+index);
//            }
          d=fpixels[chn][index];
          kernel[k][y*size+x]=d;
          if (d!=0.0) nonzero=true;
        }
        if (nonzero) {
          if (kernels[i][j]==null) {
            kernels[i][j]=new double [nChn][];
            for (l=0;l<nChn;l++) kernels[i][j][l]=null;
          }
          kernels[i][j][k]=kernel[k].clone();
        }
      }
    }
    return kernels;
  }

/* flattens array of kernels to be saved as a tiff file */
  public double [] mergeSquareKernelsToOne(double [][][][] kernels) { 
    if (kernels==null) return null;
    int tilesY=kernels.length;
    int tilesX=kernels[0].length;
    int i,j,k;
    double [][]kernel=null;
    for (i=0;(i<tilesY) && (kernel==null);i++)  for (j=0;(j<tilesX) && (kernel==null);j++)  kernel=kernels[i][j];
    if (kernel==null) return null;
    int [] channels = {-1,-1,-1,-1};
    int length=0;
    for (i=0;i<kernel.length;i++) if (kernel[i]!=null){
      length=kernel[i].length;
      break;
    }
    j=0;
    for (i=0;i<kernel.length;i++) if (kernel[i]!=null)  if (j<channels.length) channels[j++]=i;
    int size=(int) Math.sqrt(length);
    int outWidth= size*tilesX;
    int outHeight=size*tilesY;
    double [] merged = new double [4*outWidth*outHeight]; // heap space
    for (i=0;i<merged.length;i++)merged[i]=0.0;

    int x,y, index, cx,cy;
    if (DEBUG_LEVEL>1) {
       System.out.println("mergeSquareKernelsToOne(): tilesY="+tilesY+" tilesX="+tilesX);
       System.out.println("mergeSquareKernelsToOne(): channels[0]="+channels[0]+" channels[1]="+channels[1]+" channels[2]="+channels[2]+" channels[3]="+channels[3]);
    }

//    for (i=0; i<tilesY;i++)  for (j=0;j<tilesX;j++) for (k=0;k<channels.length;k++) if (channels[k]>-0) {
    for (i=0; i<tilesY;i++)  for (j=0;j<tilesX;j++) for (cy=0;cy<2;cy++) for (cx=0;cx<2;cx++) {
      k=channels[2*cy+cx];
      if (k >= 0) {
//       if (DEBUG_LEVEL>1) { System.out.println("mergeSquareKernels(): i="+i+" j="+j+" k="+k); }
        for (y=0;y<size;y++) for (x=0;x<size;x++) {
          index=(((2*i+cy)*size+y)*2*outWidth)+((2*j+cx)*size+x);
//           if (DEBUG_LEVEL>1) {
//            System.out.println("mergeSquareKernels(): i="+i+" j="+j+" k="+k+" y="+y+" x="+x+" index="+index);
//          }
          if (kernels[i][j]==null) merged[index]=0.0;
          else {
               merged[index]=kernels[i][j][k][y*size+x];
//            if (DEBUG_LEVEL>1) System.out.println("mergeSquareKernels(): merged["+k+"]["+index+"]=kernels["+i+"]["+j+"]["+k+"]["+(y*size+x)+"]="+kernels[i][j][k][y*size+x]);
          }
        }
      }
    }
    return merged;
  }

  public double [][][][] splitSquareKernelsFromOne(double [] flattened, // flattened kernels, same as produced by mergeSquareKernelsToOne()
                                                              int size, // size of each kernel (should be square)
                                                             int width, // width of the flattened kernels, in pixels (should be multiple of 2*size)
                                                       int [] channels  // four-element array of color channel numbers (usually {1,2,5,-1})
                                                                  ) { 
    if (flattened==null) return null;
    int height=flattened.length/width;
    int tilesX=width/2/size;
    int tilesY=height/2/size;
    if (flattened.length!=(tilesY*tilesX*4*size*size)) return null; // not multiple ofkernel cells
    int nChn=-1;
    int i,j,k;
    for (i=0;i<channels.length;i++) if (nChn<channels[i]) nChn = channels[i];
    nChn++;
    if (nChn==0) return null; // no channels defined
    double [][][][] kernels=new double [tilesY][tilesX][][];
    double [][]kernel=new double [nChn][];
    for (i=0;i<nChn;i++) kernel[i]=null;
    for (i=0;i<channels.length;i++) if (channels[i]>=0) kernel[channels[i]]=new double [size*size];
    boolean nonzero=false;
    int x,y, index, cx, cy, l;
    double d;
    for (i=0; i<tilesY;i++)  for (j=0;j<tilesX;j++) {
//      nonzero=false;
      kernels[i][j]=null;
      for (cy=0;cy<2;cy++) for (cx=0;cx<2;cx++) {
//        k=channels[2*cy+cx];
        k=2*cy+cx;
        k= (k>=channels.length)?-1:channels[k];
        if (k >= 0) {
          nonzero=false;
//         if (DEBUG_LEVEL>1) { System.out.println("mergeSquareKernels(): i="+i+" j="+j+" k="+k); }
/* extrac one set of kernels */
          for (y=0;y<size;y++) for (x=0;x<size;x++) {
//            index=(((2*i+cy)*size+y)*2*outWidth)+((2*j+cx)*size+x);
            index=(((2*i+cy)*size+y)*width)+((2*j+cx)*size+x);
//             if (DEBUG_LEVEL>1) {
//              System.out.println("mergeSquareKernels(): i="+i+" j="+j+" k="+k+" y="+y+" x="+x+" index="+index);
//            }
            d=flattened[index];
            kernel[k][y*size+x]=d;
            if (d!=0.0) nonzero=true;
          }
          if (nonzero) {
            if (kernels[i][j]==null) {
              kernels[i][j]=new double [nChn][];
              for (l=0;l<nChn;l++) kernels[i][j][l]=null;
            }
            kernels[i][j][k]=kernel[k].clone();
          }
        }
      }
//      if (nonzero) kernels[i][j]=kernel.clone();
//      else         kernels[i][j]=null;
    }
    return kernels;
  }


// Do not know how to combine methods for double[] and float []
  public double [][][][] splitSquareKernelsFromOne(float  [] flattened, // flattened kernels, same as produced by mergeSquareKernelsToOne()
                                                              int size, // size of each kernel (should be square)
                                                             int width, // width of the flattened kernels, in pixels (should be multiple of 2*size)
                                                       int [] channels  // four-element array of color channel numbers (usually {1,2,5,-1})
                                                                  ) { 
    if (flattened==null) return null;
    int height=flattened.length/width;
    int tilesX=width/2/size;
    int tilesY=height/2/size;
    if (flattened.length!=(tilesY*tilesX*4*size*size)) return null; // not multiple of the kernel cells
    int nChn=-1;
    int i,j,k;
    for (i=0;i<channels.length;i++) if (nChn<channels[i]) nChn = channels[i];
    nChn++;
    if (nChn==0) return null; // no channels defined
    double [][][][] kernels=new double [tilesY][tilesX][][];
    double [][]kernel=new double [nChn][];
    for (i=0;i<nChn;i++) kernel[i]=null;
    for (i=0;i<channels.length;i++) if (channels[i]>=0) kernel[channels[i]]=new double [size*size];
    boolean nonzero=false;
    int x,y, index, cx,cy, l;
    double d;
    for (i=0; i<tilesY;i++)  for (j=0;j<tilesX;j++) {
//      nonzero=false;
      kernels[i][j]=null;
      for (cy=0;cy<2;cy++) for (cx=0;cx<2;cx++) {
//        k=channels[2*cy+cx];
        k=2*cy+cx;
        k= (k>=channels.length)?-1:channels[k];
        if (k >= 0) {
          nonzero=false;
//         if (DEBUG_LEVEL>1) { System.out.println("mergeSquareKernels(): i="+i+" j="+j+" k="+k); }
/* extract one set of kernels */
          for (y=0;y<size;y++) for (x=0;x<size;x++) {
//            index=(((2*i+cy)*size+y)*2*outWidth)+((2*j+cx)*size+x);
            index=(((2*i+cy)*size+y)*width)+((2*j+cx)*size+x);
//             if (DEBUG_LEVEL>1) {
//              System.out.println("mergeSquareKernels(): i="+i+" j="+j+" k="+k+" y="+y+" x="+x+" index="+index);
//            }
            d=flattened[index];
            kernel[k][y*size+x]=d;
            if (d!=0.0) nonzero=true;
          }
          if (nonzero) {
            if (kernels[i][j]==null) {
              kernels[i][j]=new double [nChn][];
              for (l=0;l<nChn;l++) kernels[i][j][l]=null;
            }
            kernels[i][j][k]=kernel[k].clone();
          }
        }
      }
//      if (nonzero) kernels[i][j]=kernel.clone();
//      else         kernels[i][j]=null;
    }
    return kernels;
  }




/* Decode combined kernels into individual bayer component ones */
  public double[][] decodeKernels (double[] encodedKernels, // combined kernels
                                                 int ksize, // output size (side of square) of the extracted individual kernels
                                         int combineGreens) { // extract combined greens: 0 - individual, 1 - diagonal, 2 - checker
    boolean [] enComponent={(combineGreens==0),true,true,(combineGreens==0),(combineGreens==1),(combineGreens==2)};
    int size= (int) Math.sqrt(encodedKernels.length);
    int n,i;
    int outX,outY, kernelX,kernelY;
    int dx=0;
    int dy=0;
    double [][] kernels=new double [enComponent.length][];
    for (n=0;n<kernels.length;n++) {
      if (enComponent[n]) {
        kernels[n]=new double [ksize*ksize];
       for (i=0;i<kernels[n].length;i++) kernels[n][i]=0.0;
      } else {
        kernels[n]=null;
        continue;
      }
      if (n!=4) {
        switch (n) {
         case 5: 
         case 0: dy=0; dx=0; break;
         case 1: dy=0; dx=1; break;
         case 2: dy=1; dx=0; break;
         case 3: dy=1; dx=1; break;
        }
        for (outY=dy;outY<size;outY+=2) {
          kernelY=(outY+ksize-size/2)/2;
          for (outX=dx;outX<size;outX+=2) {
            kernelX=(outX+ksize-size/2)/2;
            kernels[n][kernelY*ksize+kernelX]=encodedKernels[outY*size+outX];
          }
        }
      } else { //n>=4
        for (outY=0;outY<size;outY++) for (outX=(outY & 1); outX<size;outX+=2) {
          kernelX= (outX-outY+ksize)/2;       //           ((outX-size/2) - (outY-size/2))/2 + ksize/2;
          kernelY= (outX+outY+ksize-size)/2;  //    ((outX-size/2) + (outY-size/2))/2 + ksize/2;
          if ((kernelY>=0) && (kernelY<ksize) && (kernelX>=0) && (kernelX<ksize)) {
            kernels[n][kernelY*ksize+kernelX]=encodedKernels[outY*size+outX];
          }
        }
      }
    }
    return kernels;
  }

/* Combine kernels (direct PSF or reversed one) for R,B and combined G into a single array, where each color
    component kernel uses corresponding pixels places  */
  public double[] encodeKernels (double[][]kernels, // arrays of color component kernels (or nulls), each kernel as a 1-d array
                                          int size,  // size of the side of the square of the result
                                 boolean normalize) {
     int n,i;
     int ksize=0;
     boolean [] enComponent=new boolean[kernels.length];
     for (i=0;i<enComponent.length;i++) enComponent[i]=(kernels[i]!=null);
/* disable individual greens if any of the composite ones is defined */
     if (enComponent[4] || enComponent[5]) {
       enComponent[0]=false;
       enComponent[3]=false;
     }
/* disable composite diagonal if composite "checker" is defined */
     if (enComponent[5]) {
       enComponent[4]=false;
     }
     for (i=0;i<kernels.length;i++) if (kernels[i]!=null) {
        ksize=(int) Math.sqrt(kernels[i].length);
        break;
     }
     if (ksize==0) {
       if (DEBUG_LEVEL>1) System.out.println("encodeKernels - nothing to encode, all kernels are null");
       return null;
     }
     double [] encodedKernels = new double[size*size];
     int x,y,outX,outY, kernelX,kernelY;
     int dx=0;
     int dy=0;
     double d,sum;
//     int   num;
     for (i=0;i<encodedKernels.length;i++) encodedKernels[i]=0.0;
     for (n=0;n<kernels.length;n++) if (enComponent[n]) {
       sum=0.0;
//       num=0;
       if (n!=4) {
         switch (n) {
          case 5: 
          case 0: dy=0; dx=0; break;
          case 1: dy=0; dx=1; break;
          case 2: dy=1; dx=0; break;
          case 3: dy=1; dx=1; break;
         }
         
         for (y=0;y<size/2;y++) {
           outY=2*y+dy;
           kernelY=y+((2*ksize-size)/4);
           if ((kernelY>=0) && (kernelY<ksize)) for (x=0;x<size/2;x++) {
             outX=2*x+dx;
             kernelX=x+((2*ksize-size)/4);
             if ((kernelX>=0) && (kernelX<ksize)) {
               d=kernels[n][kernelY*ksize+kernelX];
               encodedKernels[outY*size+outX]=d;
               sum+=d;
//               num++;
             }
           }
         }
         if (normalize) {
           d=1.0/sum;
           for (outY=dy;outY<size;outY+=2) for (outX=dx;outX<size;outX+=2) encodedKernels[outY*size+outX]*=d;
         }
       } else { // n==4
         for (outY=0;outY<size;outY++) for (outX=(outY & 1); outX<size;outX+=2) {
           kernelX= (outX-outY+ksize)/2;       //           ((outX-size/2) - (outY-size/2))/2 + ksize/2;
           kernelY= (outX+outY+ksize-size)/2;  //    ((outX-size/2) + (outY-size/2))/2 + ksize/2;
           if ((kernelY>=0) && (kernelY<ksize) && (kernelX>=0) && (kernelX<ksize)) {
             d=kernels[n][kernelY*ksize+kernelX];
             encodedKernels[outY*size+outX]=d;
             sum+=d;
//             num++;
           }
         }
         if (normalize) {
           d=1.0/sum;
           for (outY=0;outY<size;outY++) for (outX=(outY & 1); outX<size; outX+=2) encodedKernels[outY*size+outX]*=d;
         }
       }
     }
     return encodedKernels;
  }






//=============================================

/* create higher resolution boolean array, including margins  - 4 pixels each side. Size - number of Bayer cells*/

  public boolean [][]   simulatePatternFull(double freqX1,
                                     double freqY1,
                                     double phase1,
                                     double freqX2,
                                     double freqY2,
                                     double phase2,
                                     int subdiv,
                                     int size,
                                     boolean center_for_g2) {
    return simulatePatternFull(freqX1, freqY1, phase1, freqX2, freqY2, phase2, null, subdiv, size, center_for_g2);
  }
  public boolean [][]   simulatePatternFull(double freqX1,
                                     double freqY1,
                                     double phase1,
                                     double freqX2,
                                     double freqY2,
                                     double phase2,
                                     double [] corr,
                                     int subdiv,
                                     int size,
                                     boolean center_for_g2) {
    int i,j;
    int fullSize=subdiv*(size+4)*2;
    boolean [][] barray=new boolean [fullSize][fullSize];
    double xl,yl,x,y,p1,p2;
    if (DEBUG_LEVEL>1) {
      System.out.println("simulatePatternFull:");
      System.out.println(" Ax="+IJ.d2s(corr[0],5)+" Bx="+IJ.d2s(corr[1],5)+" Cx="+IJ.d2s(corr[2],5)+" Dx="+IJ.d2s(corr[6],5)+" Ex="+IJ.d2s(corr[7],5));
      System.out.println(" Ay="+IJ.d2s(corr[3],5)+" By="+IJ.d2s(corr[4],5)+" Cy="+IJ.d2s(corr[5],5)+" Dy="+IJ.d2s(corr[8],5)+" Ey="+IJ.d2s(corr[9],5));
    }
    for (i=0;i<fullSize;i++) {
      yl=(i-0.5*fullSize)/subdiv-(center_for_g2?0.5:1.0); // center in the middle of Bayer
      for (j=0;j<fullSize;j++) {
        xl=(j-0.5*fullSize)/subdiv-(center_for_g2?0.5:1.0); // center in the middle of Bayer

/* apply second order poliniminal correction to x,y
    x=xl+Ax*xl^2+Bx*yl^2+2*Cx*xl*yl;
    y=xl+Ay*xl^2+By*yl^2+2*Cy*xl*yl; */
        if (corr==null) {
          x=xl;
          y=yl;
        } else {
          x=xl + corr[0]*xl*xl + corr[1]*yl*yl + 2* corr[2]*xl*yl + corr[6]*xl + corr[7]*yl;
          y=yl + corr[3]*xl*xl + corr[4]*yl*yl + 2* corr[5]*xl*yl + corr[8]*xl + corr[9]*yl;
        }
        p1=y*freqY1+x*freqX1+(phase1/(Math.PI*2));
        p1-=Math.floor(p1);
        p2=y*freqY2+x*freqX2+(phase2/(Math.PI*2));
        p2-=Math.floor(p2);
        barray[i][j]=!(((p1<0.25) || (p1>=0.75) )^ ((p2<0.25) || (p2>=0.75) ));
      }
    }
    return barray;
  }

  public boolean [][]   simulatePatternFullPattern(
                                     boolean []pattern, // square pattern
                                     double freqX1,
                                     double freqY1,
                                     double phase1,
                                     double freqX2,
                                     double freqY2,
                                     double phase2,
                                     double [] corr,
                                     int subdiv,
                                     int size,
                                     boolean center_for_g2) {
    int patternSize= (pattern!=null)?((int) Math.sqrt(pattern.length)):0;
    double twicePatternSize=2*patternSize;
    int i,j;
    int fullSize=subdiv*(size+4)*2;
    boolean [][] barray=new boolean [fullSize][fullSize];
    double xl,yl; //,x,y;//,p1,p2;


    double [][] xy2uv= {{freqX1,freqY1},
                        {freqX2,freqY2}};

    if (DEBUG_LEVEL>2) {
      System.out.println("simulatePatternFullPattern:");
      System.out.println(" Ax="+IJ.d2s(corr[0],5)+" Bx="+IJ.d2s(corr[1],5)+" Cx="+IJ.d2s(corr[2],5)+" Dx="+IJ.d2s(corr[6],5)+" Ex="+IJ.d2s(corr[7],5));
      System.out.println(" Ay="+IJ.d2s(corr[3],5)+" By="+IJ.d2s(corr[4],5)+" Cy="+IJ.d2s(corr[5],5)+" Dy="+IJ.d2s(corr[8],5)+" Ey="+IJ.d2s(corr[9],5));
    }

    if (DEBUG_LEVEL>2) {
      System.out.println("simulatePatternFullPattern:  xy2uv[0][0]="+IJ.d2s(xy2uv[0][0],4)+" xy2uv[0][1]="+IJ.d2s(xy2uv[0][1],4));
      System.out.println("simulatePatternFullPattern:  xy2uv[1][0]="+IJ.d2s(xy2uv[1][0],4)+" xy2uv[1][1]="+IJ.d2s(xy2uv[1][1],4));
    }

    double []uv, xy;
    xy=new double [2];
    
    double [] phases={phase1/(Math.PI*2)+0.25,phase2/(Math.PI*2)+0.25}; // period=1.0;
    int iu,iv;
    boolean invert;
    for (i=0;i<fullSize;i++) {
      yl=(i-0.5*fullSize)/subdiv-(center_for_g2?0.5:1.0); // center in the middle of Bayer
      for (j=0;j<fullSize;j++) {
        xl=(j-0.5*fullSize)/subdiv-(center_for_g2?0.5:1.0); // center in the middle of Bayer

/* apply second order poliniminal correction to x,y
    x=xl+Ax*xl^2+Bx*yl^2+2*Cx*xl*yl;
    y=xl+Ay*xl^2+By*yl^2+2*Cy*xl*yl; */
        if (corr==null) {
          xy[0]=xl;
          xy[1]=yl;
        } else {
          xy[0]=xl + corr[0]*xl*xl + corr[1]*yl*yl + 2* corr[2]*xl*yl + corr[6]*xl + corr[7]*yl;
          xy[1]=yl + corr[3]*xl*xl + corr[4]*yl*yl + 2* corr[5]*xl*yl + corr[8]*xl + corr[9]*yl;
        }
        uv= matrix2x2_mul(xy2uv, xy);
        uv= vector_add(uv,phases);
        uv[0]-=Math.floor(uv[0]);
        uv[1]-=Math.floor(uv[1]);
        invert=false;
        if (uv[0]>=0.5){
          invert=!invert;
          uv[0]-=0.5;
        }
        if (uv[1]>=0.5){
          invert=!invert;
          uv[1]-=0.5;
        }
        if (pattern==null) {
          barray[i][j]=!invert;
        } else {
          iu= (int) Math.round(uv[0]*twicePatternSize);
          iv= (int) Math.round(uv[1]*twicePatternSize);
          if ((iu<0) || (iu>=patternSize)) {
            invert=!invert;
            iu=(iu+patternSize)% patternSize;
          }
          if ((iv<0) || (iv>=patternSize)) {
            invert=!invert;
            iv=(iv+patternSize)% patternSize;
          }
          barray[i][j]=invert ^ pattern[iv*patternSize + iu];
        }
      }
    }
    return barray;
  }


  public boolean [] patternGenerator(int size,
                            int patternNumber,
                       double patternModifier) {
    boolean [] pattern=new boolean [size*size];
    int i,j,index,k;
    boolean p;
    double a,r,r2,h;
    double qSize=size/4;
    switch (patternNumber) {
      case 1:
              a=patternModifier*(Math.sqrt(2)-1.0);
              r=(a*a+1)/(2*a)*qSize;
              r2=r*r;
              h=Math.sqrt(r2-qSize*qSize);
              if (a>1.0) h=-h;
              double [][] pattern1Centers={{qSize,         -h},
                                           { size+h,       qSize},
                                           { size-qSize,   size+h},
                                           {-h,            size-qSize}};
              index=0;
              for (i=0;i<size;i++) for (j=0;j<size;j++) {
                p=true;
                for (k=0;k<pattern1Centers.length;k++) if ((((i-pattern1Centers[k][1])*(i-pattern1Centers[k][1])+(j-pattern1Centers[k][0])*(j-pattern1Centers[k][0])))<r2) p=false;
                pattern[index++]=p;
              }
              break;
      case 2:
              index=0;
              for (i=0;i<size;i++) for (j=0;j<size;j++) {
                p= (i>=0.3*size) && (i<0.7*size) && (j>=0.3*size) && (j<0.7*size);
                pattern[index++]=p;
              }
              break;
      case 3:
              index=0;
              for (i=0;i<size;i++) for (j=0;j<size;j++) {
                p= (i>=0.1*size) && (i<0.9*size) && (j>=0.1*size) && (j<0.9*size);
                pattern[index++]=p;
              }
              break;
      default: for (index=0;index<pattern.length;index++) pattern[index]=true;
    }
    return pattern;
  }


/* make it faster when outSubdiv =2*n (usually so) */
/* TODO: cleanup shifts - they seem now to work correctly */
  public double [][] extractSimulPatterns (boolean [][] barray,   // high resolution boolean pattern array
                                            double  simul_fill, // part of the (center) pixel area being "phptosensitive"
                                                    int subdiv,   // boolean pixels to real pixels resolution
                                                 int outSubdiv,  // subdivide output pixels
                                                      int size,    // number of Bayer cells in width of the square selection (half number of pixels)
                                                     double x0,    // selection center, X (in pixels)
                                                     double y0) {
    int sampleWidth=(int) (Math.sqrt(simul_fill)*subdiv);
    int sampleN=sampleWidth*sampleWidth;
    if      (sampleWidth<1)     sampleWidth=1;
    else if (sampleWidth>subdiv)sampleWidth=subdiv;
    double sampleAverage=0.5*sampleN;

    int n,i,j;
    int fullSize=barray.length;
    double [][] simul_pixels=new double [5][size*size];
    int ix,iy, iy0,ix0,px,py;
    double bx,by;
    double s;
    double span=((double) size)/outSubdiv;
    int sampLow=-sampleWidth/2;
    int sampHigh=sampLow+sampleWidth;
    for (n=0;n<4;n++) {
      bx=(n&1)-0.5+0.5;        // last 0.5 to make same center as for dual greens
      by=((n>>1) & 1)-0.5-0.5;// last 0.5 to make same center as for dual greens
      for (iy=0;iy<size;iy++) {
        iy0=(fullSize/2) + (int) ((-span+y0+by  +1.5  +2.0*iy/outSubdiv)*subdiv);
        for (ix=0;ix<size;ix++) {
          ix0=(fullSize/2) + (int) ((-span+x0+bx+0.5   +2.0*ix/outSubdiv)*subdiv);
          s=0.0;
          for (py=iy0+sampLow;py<iy0+sampHigh;py++) for (px=ix0+sampLow;px<ix0+sampHigh;px++) {
            if (barray[py][px]) s+=1.0;
          }
          simul_pixels[n][iy*size+ix]= (s-sampleAverage)/sampleAverage;
        }
      }

    }
    if (outSubdiv>1) {
if (DEBUG_LEVEL>2)System.out.println("Generating combined greens pattern greens from scratch");
      n=4;
      bx=0.0;
      by=0.0;
      for (iy=0;iy<size;iy++) {
        for (ix=0;ix<size;ix++) {
          iy0=(fullSize/2) + (int) ((-span+y0+by-1+1.5  +1.0*(size+iy-ix)/outSubdiv)*subdiv);
          ix0=(fullSize/2) + (int) ((-span+x0+bx  +0.5  +1.0*(iy+ix)/outSubdiv)*subdiv);
          s=0.0;
          for (py=iy0+sampLow;py<iy0+sampHigh;py++) for (px=ix0+sampLow;px<ix0+sampHigh;px++) {
            if (barray[py][px]) s+=1.0;
          }
          simul_pixels[n][iy*size+ix]= (s-sampleAverage)/sampleAverage;
        }
      }
    } else { // just reuse available greens
if (DEBUG_LEVEL>2)System.out.println("Generating combined greens pattern from individual greens");
/* now combine greens - same as in splitBayer() */

      int base, base_b;
      base_b=0;
      for (i=0;i<size/2; i++){
        base=size*size/2+ i* (size+1);
        for (j=0; j<size/2; j++) {
          simul_pixels[4][base_b++]=simul_pixels[0][base];
          base-=size;
          simul_pixels[4][base_b++]=simul_pixels[3][base++];
        }
        base=size*size/2+ i* (size+1);
        for (j=0; j<size/2; j++) {
//System.out.println("2:y="+y+" x="+x+" base_b="+base_b+" base="+base);
          simul_pixels[4][base_b++]=simul_pixels[3][base++];
          simul_pixels[4][base_b++]=simul_pixels[0][base];
          base-=size;
        }
      }
    }
 if (DEBUG_LEVEL>2) {
    System.out.println("extractSimulPatterns, x0="+x0+" y0="+y0+" fullSize="+fullSize+" size="+size+" subdiv="+subdiv+" outSubdiv="+outSubdiv);
    System.out.println(" sampLow="+sampLow+" sampHigh="+sampHigh+" span="+span+" size="+size);
    for (n=0;n<simul_pixels.length;n++) {
      s=0.0;
      for (i=0;i<simul_pixels[n].length;i++) s+=simul_pixels[n][i];
      System.out.println(" component="+i+" sum of pixels="+s);
    }
 }

    if (DEBUG_LEVEL>2) SDFA_instance.showArrays(simul_pixels,size,size, "SIMUL");

    return simul_pixels;
  }

  public double [][] old_extractSimulPatterns (boolean [][] barray,   // high resolution boolean pattern array
                                                   int subdiv,   // boolean pixels to real pixels resolution
                                                    int size,    // number of Bayer cells in width of the square selection (half of the number of the pixels)
                                                      int x0,    // selection center, X (in bayer cells)
                                                      int y0) {  // selection center, Y (in bayer cells)
    return extractSimulPatterns (barray, subdiv, size, x0, y0, 0, 0);
  }

  public double [][] extractSimulPatterns (boolean [][] barray,   // high resolution boolean pattern array
                                                   int subdiv,   // boolean pixels to real pixels resolution
                                                    int size,    // number of Bayer cells in width of the square selection (half number of pixels)
                                                      int x0,    // selection center, X (in bayer cells)
                                                      int y0,    // selection center, Y (in bayer cells)
                                                      int dx,    // sub-pixel shift in x-direction (+/- 2*subdiv)
                                                      int dy) {
    int n,i,j,k,l;
    int []indx=new int[4];
    int fullSize=barray.length;
    double [][] simul_pixels=new double [5][size*size];
//  double xl,yl,x,y,p1,p2;
    int i0,j0,i1,j1;
    double d= 1.0/subdiv/subdiv;
    i0=(fullSize/2) + (2*y0*subdiv) -(size*subdiv);
    j0=(fullSize/2) + (2*x0*subdiv) -(size*subdiv);
    i1=i0+size*subdiv*2;
    j1=j0+size*subdiv*2;
 if (DEBUG_LEVEL>2) {
    System.out.println("extractSimulPatterns, x0="+x0+" y0="+y0+" i0="+i0+" j0="+j0+" i1="+i1+" j1="+j1+" fullSize="+fullSize+" size="+size);
 }
    for (n=0;n<4;n++) indx[n]=0;
    for (i=i0;i<i1;i+=subdiv) for (j=j0;j<j1;j+=subdiv) {
      n=(((i/subdiv) & 1)<<1) | ((j/subdiv) & 1);
      simul_pixels[n][indx[n]]=0.0;
      for (k=i;k<(i+subdiv);k++) for (l=j;l<(j+subdiv);l++) if (barray[k+dy][l+dx])  simul_pixels[n][indx[n]]+=d;
      indx[n]++;
    }
/* now combine greens - same as in splitBayer() */
    int base, base_b;
    base_b=0;
    for (i=0;i<size/2; i++){
      base=size*size/2+ i* (size+1);
      for (j=0; j<size/2; j++) {
        simul_pixels[4][base_b++]=simul_pixels[0][base];
        base-=size;
        simul_pixels[4][base_b++]=simul_pixels[3][base++];
      }
      base=size*size/2+ i* (size+1);
      for (j=0; j<size/2; j++) {
//System.out.println("2:y="+y+" x="+x+" base_b="+base_b+" base="+base);
        simul_pixels[4][base_b++]=simul_pixels[3][base++];
        simul_pixels[4][base_b++]=simul_pixels[0][base];
        base-=size;
      }
    }
    return simul_pixels;
  }
//=============================================



//=============================================



  public double [][]   simulatePattern(double freqX1,
                                     double freqY1,
                                     double phase1,
                                     double freqX2,
                                     double freqY2,
                                     double phase2,
                                     int subdiv,
                                     int size,
                                     boolean center_for_g2) {
    return simulatePattern(freqX1, freqY1, phase1, freqX2, freqY2, phase2, null, subdiv, size, center_for_g2);
  }
  public double [][]   simulatePattern(double freqX1,
                                     double freqY1,
                                     double phase1,
                                     double freqX2,
                                     double freqY2,
                                     double phase2,
                                     double [] corr,
                                     int subdiv,
                                     int size,
                                     boolean center_for_g2) {
  int n,i,j,k,l;
  int []indx=new int[4];
  int fullSize=subdiv*size*2;
  boolean [][] barray=new boolean [fullSize][fullSize];
  double [][] simul_pixels=new double [5][size*size];
  double xl,yl,x,y,p1,p2;
  for (i=0;i<fullSize;i++) {
    yl=(i-0.5*fullSize)/subdiv-(center_for_g2?0.5:1.0); // center in the middle of Bayer
    for (j=0;j<fullSize;j++) {
      xl=(j-0.5*fullSize)/subdiv-(center_for_g2?0.5:1.0); // center in the middle of Bayer

/* apply second order poliniminal correction to x,y
    x=xl+Ax*xl^2+Bx*yl^2+Cx*xl*yl;
    y=xl+Ay*xl^2+By*yl^2+Cy*xl*yl; */
      if (corr==null) {
        x=xl;
        y=yl;
      } else {
//        x=xl + corr[0]*xl*xl + corr[1]*yl*yl + corr[2]*xl*yl;
//        y=yl + corr[3]*xl*xl + corr[4]*yl*yl + corr[5]*xl*yl;
          x=xl + corr[0]*xl*xl + corr[1]*yl*yl + 2* corr[2]*xl*yl + corr[6]*xl + corr[7]*yl;
          y=yl + corr[3]*xl*xl + corr[4]*yl*yl + 2* corr[5]*xl*yl + corr[8]*xl + corr[9]*yl;

      }
      p1=y*freqY1+x*freqX1+(phase1/(Math.PI*2));
      p1-=Math.floor(p1);
      p2=y*freqY2+x*freqX2+(phase2/(Math.PI*2));
      p2-=Math.floor(p2);
      barray[i][j]=!(((p1<0.25) || (p1>=0.75) )^ ((p2<0.25) || (p2>=0.75) ));
    }
  }
  double d= 1.0/subdiv/subdiv;
  for (n=0;n<4;n++) indx[n]=0;
  for (i=0;i<fullSize;i+=subdiv) for (j=0;j<fullSize;j+=subdiv) {
    n=(((i/subdiv) & 1)<<1) | ((j/subdiv) & 1);
    simul_pixels[n][indx[n]]=0.0;
    for (k=i;k<(i+subdiv);k++) for (l=j;l<(j+subdiv);l++) if (barray[k][l])  simul_pixels[n][indx[n]]+=d;
    indx[n]++;
  }
/// now combine greens - same as in splitBayer()
  int base, base_b;
  base_b=0;
  for (i=0;i<size/2; i++){
    base=size*size/2+ i* (size+1);
    for (j=0; j<size/2; j++) {
      simul_pixels[4][base_b++]=simul_pixels[0][base];
      base-=size;
      simul_pixels[4][base_b++]=simul_pixels[3][base++];
    }
    base=size*size/2+ i* (size+1);
    for (j=0; j<size/2; j++) {
//System.out.println("2:y="+y+" x="+x+" base_b="+base_b+" base="+base);
      simul_pixels[4][base_b++]=simul_pixels[3][base++];
      simul_pixels[4][base_b++]=simul_pixels[0][base];
      base-=size;
    }
  }
  return simul_pixels;
 }



  public double[][] normalizeAndWindow (double [][] pixels, double [] hamming) {
    return normalizeAndWindow (pixels, hamming, true);
  }

  public double[] normalizeAndWindow (double [] pixels, double [] hamming) {
    return normalizeAndWindow (pixels, hamming, true);
  }


  public double[][] normalizeAndWindow (double [][] pixels, double [] hamming, boolean removeDC) {
    int i;
    for (i=0;i<pixels.length;i++)  if (pixels[i]!=null) pixels[i]=normalizeAndWindow (pixels[i],  hamming, removeDC);
    return pixels;
  }


  public double[] normalizeAndWindow (double [] pixels, double [] hamming, boolean removeDC) {
    int j;
    double s=0.0;
    if (pixels==null) return null;
    if (removeDC) {
      for (j=0;j<pixels.length;j++) s+=pixels[j];
      s/=pixels.length;
    }
    for (j=0;j<pixels.length;j++) pixels[j]=(pixels[j]-s)*hamming[j];
    return pixels;
  }


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

  public double [][] extendFFTInput (double[][] input_pixels,
                                                   int width,   // width of the image
                                              int subDivFreq) {
    double [][] pixels=new double[input_pixels.length][];
    int i;
    for (i=0;i<pixels.length;i++) pixels[i]= extendFFTInput (input_pixels[i], width, subDivFreq);
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





/* inserts zeros between pixels */ 
  public double [][] oversampleFFTInput (double[][] input_pixels,
                                                       int ratio) {
    double [][] pixels=new double[input_pixels.length][];
    int i;
    for (i=0;i<pixels.length;i++) pixels[i]= oversampleFFTInput (input_pixels[i], ratio);
    return pixels;
  }

  public double [][] oversampleFFTInput (double[][] input_pixels,
                                                       int width,   // width of the image
                                                       int ratio) {
    double [][] pixels=new double[input_pixels.length][];
    int i;
    for (i=0;i<pixels.length;i++) pixels[i]= oversampleFFTInput (input_pixels[i], width, ratio);
    return pixels;
  }

  public double [] oversampleFFTInput (double[] input_pixels, int ratio) {
    if (input_pixels==null) return null;
    int width=(int) Math.sqrt(input_pixels.length);
    return oversampleFFTInput (input_pixels,
                                      width,   // width of the image
                                      ratio);
  }

  public double [] oversampleFFTInput (double[] input_pixels,
                                                   int width,   // width of the image
                                                   int ratio) {
    if (input_pixels==null) return null;
 if (DEBUG_LEVEL>2) System.out.println ("oversampleFFTInput(), width="+width+" ratio="+ratio+" input_pixels.length="+input_pixels.length);
    double [] pixels=new double[input_pixels.length*ratio*ratio];
    int i,j,x,y;
//    int size=(int) Math.sqrt(input_pixels.length);
    int height=input_pixels.length/width;
    for (i=0;i<pixels.length;i++) pixels[i]=0.0;
    j=0;
    for (y=0;y<height;y++) {
      i=width*ratio*ratio*y;
      for (x=0;x<width;x++) {
        pixels[i]=input_pixels[j++];
        i+=ratio;
      }
    }
    if (DEBUG_LEVEL>2) System.out.println ("oversampleFFTInput(), pixels.length="+pixels.length);
    return pixels;
  }




/* Combine both greens as a checkerboard pattern (after oversampleFFTInput()) */

  public double [][] combineCheckerGreens (double[][] input_pixels,   // pixel arrays after oversampleFFTInput() or extractSimulPatterns())
                                                         int ratio) { // same as used in oversampleFFTInput() - oversampling ratio
    int width=(int) Math.sqrt(input_pixels[0].length);
    return combineCheckerGreens (input_pixels,   // pixel arrays after oversampleFFTInput() or extractSimulPatterns())
                                        width,   // width of the image
                                        ratio);
  }

  public double [][] combineCheckerGreens (double[][] input_pixels,   // pixel arrays after oversampleFFTInput() or extractSimulPatterns())
                                                         int width,   // width of the image
                                                         int ratio) { // same as used in oversampleFFTInput() - oversampling ratio
    if (DEBUG_LEVEL>5) System.out.println ("combineCheckerGreens(), ratio="+ratio+" input_pixels.length="+input_pixels.length);

    if ((ratio<2) ||
        (input_pixels==null) ||
        ((input_pixels.length>5) && (input_pixels[5]!=null)) ||
        (input_pixels.length<4) ||
        (input_pixels[0]==null) ||
        (input_pixels[3]==null)) return input_pixels;
//    int size=(int) Math.sqrt(input_pixels[0].length);
    int height=input_pixels[0].length/width;
//    double [][] pixels={input_pixels[0],input_pixels[1],input_pixels[2],input_pixels[3],input_pixels[4],null};
    int i,j;
    double [][] pixels={null,null,null,null,null,null};
    for (i=0;i<input_pixels.length;i++) pixels[i]=input_pixels[i];
    pixels[5]= new double[input_pixels[0].length];
    int index=0;
    int index_diff=(width+1)*ratio/2;
    double d;
    for (i=0;i<height;i++) for (j=0;j<width;j++) {
      d=input_pixels[0][index];
      if ((i>=ratio) && (j>=ratio)) d=0.5*(d+input_pixels[3][index-index_diff]);
      pixels[5][index++]=d;
    }

    if (DEBUG_LEVEL>5) {
      for (j=0;j<pixels.length;j++) if (pixels[j]!=null) {
        d=0.0;
        for (i=0;i<pixels[j].length;i++) d+=pixels[j][i];
        System.out.println ("combineCheckerGreens(),  sum of pixels["+j+"]="+d);
      }
    }

    return pixels;
  }




  public double[][] swapQuadrants (double[][] pixels) {
    int i;
    for (i=0;i<pixels.length;i++) pixels[i]=swapQuadrants (pixels[i]);
    return pixels;
  }

  public double[] swapQuadrants (double[] pixels) {
    int n,m, i,j;
    double d;
    n=0;
    for (i=1; i<pixels.length; i=i<<2) n++;
    m=(1<<(n-1)) | (1<< (2*n-1)) ;
    for (j=0;j< (pixels.length>>1);j++) {
      d=pixels[j];
      pixels[j]= pixels[j^ m];
      pixels[j^ m]=d;
    }
    return pixels;
  }

double [][] matrix2x2_invert(double [][] m ){
    double det=m[0][0]*m[1][1]-m[0][1]*m[1][0];
    double [][] rslt= {{ m[1][1]/det,  -m[0][1]/det},
                       {-m[1][0]/det,   m[0][0]/det}};
    return rslt;
}

double [][] matrix2x2_mul(double [][] a, double [][] b ){
    double [][] rslt={{a[0][0]*b[0][0]+a[0][1]*b[1][0], a[0][0]*b[0][1]+a[0][1]*b[1][1]},
                      {a[1][0]*b[0][0]+a[1][1]*b[1][0], a[1][0]*b[0][1]+a[1][1]*b[1][1]}};
    return rslt;
}

double [] matrix2x2_mul(double [][] a, double [] b ){
    double [] rslt={a[0][0]*b[0]+a[0][1]*b[1],
                    a[1][0]*b[0]+a[1][1]*b[1]};
    return rslt;
}

double [][] matrix2x2_scale(double [][] a, double  b ){
    double [][] rslt={{a[0][0]*b, a[0][1]*b},
                      {a[1][0]*b, a[1][1]*b}};
    return rslt;
}

double [][] matrix_add(double [][] a, double [][]  b ){
    double [][] rslt= new double [a.length][a[0].length];
    int i,j;
    for (i=0;i<rslt.length;i++) for (j=0;j<rslt[0].length;j++) rslt[i][j]=a[i][j]+b[i][j];
    return rslt;
}

double []    vector_add(double [] a, double [] b ){
    double [] rslt= new double [a.length];
    int i;
    for (i=0;i<rslt.length;i++) rslt[i]=a[i]+b[i];
    return rslt;
}



//===============
/* Low pass filter using convolution with gauss (source should have zero in the center of the array)  */
/* Replaced by using double FHT */
double [] lowPassGaussFloat(double [] pixels, double gaussSigma, boolean centered) {
   if (pixels==null) return null;
   int size = (int) Math.sqrt(pixels.length);
   int i,j,indx;
   double [] gauss=new double [size];
   double k=4.0*gaussSigma*gaussSigma/size/size;
   for (i=0;i<=size/2;i++) {
     gauss[i]=Math.exp(-(k*i*i));
     if (i>0) gauss[size-i]=gauss[i];
   }
   float [] floatPixels = new float [pixels.length];
   for (i=0; i<pixels.length;i++) floatPixels[i]= (float) pixels[i];
   ImageProcessor ip = new FloatProcessor(size,size);
   ip.setPixels(floatPixels);
   FHT fht =  new FHT(ip);
// Swapping quadrants, so the center will be 0,0
   if (centered) fht.swapQuadrants();
// get to frequency domain
   fht.transform();
   floatPixels=(float []) fht.getPixels();
   for (i=0;i<size;i++) for (j=0;j<size;j++) {
      indx=i*size+j;
      floatPixels[indx]= (float) (floatPixels[indx]*gauss[i]*gauss[j]);
   }
   ip.setPixels(floatPixels);
   fht.inverseTransform();
   if (centered) fht.swapQuadrants();
   floatPixels=(float []) fht.getPixels();
   double [] result=new double [floatPixels.length];
   for (i=0; i<pixels.length;i++) result[i]=floatPixels[i];
   return result;
}

double [] lowPassGauss(double [] pixels, double gaussSigma, boolean centered) {
   double [] clonedPixels=pixels.clone();
   if (lowPassGaussInplace(clonedPixels, gaussSigma, centered)) return clonedPixels;
   else return null;
}
boolean lowPassGaussInplace(double [] pixels, double gaussSigma, boolean centered) {
   if (pixels==null) return false;
   int size = (int) Math.sqrt(pixels.length);
   int i,j,indx;
   double [] gauss=new double [size];
   double k=4.0*gaussSigma*gaussSigma/size/size;
   for (i=0;i<=size/2;i++) {
     gauss[i]=Math.exp(-(k*i*i));
     if (i>0) gauss[size-i]=gauss[i];
   }
// Swapping quadrants, so the center will be 0,0
   if (centered) fht_instance.swapQuadrants(pixels);
// get to frequency domain
   fht_instance.transform(pixels);
   for (i=0;i<size;i++) for (j=0;j<size;j++) {
      indx=i*size+j;
      pixels[indx]*=gauss[i]*gauss[j];
   }
   fht_instance.inverseTransform(pixels);
   if (centered) fht_instance.swapQuadrants(pixels);
   return true;
}

double [] psf2mtf (double [] psf, // square array of psf pixels
                boolean centered, // true if PSF center is in the center of the array, false if it is at point 0
               boolean normalize ) { // normalize so mtf(0)=1.0
   int size=(int) Math.sqrt(psf.length);
   double [] mtf=psf.clone();
   if (centered) fht_instance.swapQuadrants(mtf);
// get to frequency domain
   fht_instance.transform(mtf);
   mtf=fht_instance.calculateAmplitude(mtf);
   if (normalize && (mtf[(size+1)*size/2]!=0.0)) {
     double k=1.0/mtf[(size+1)*size/2];
     for (int i=0;i<mtf.length;i++) mtf[i]*=k;
   }
   return mtf;
}
/* calculate amplitudes and full phases for 1/2 of the frequencies (as the input is real) */

double [][][] psf2otf (double [] psf, // square array of psf pixels
                boolean centered ) { // true if PSF center is in the center of the array, false if it is at point 0
   int size=(int) Math.sqrt(psf.length);
   double [] otf=psf.clone();
   if (centered) fht_instance.swapQuadrants(otf);
// get to frequency domain
   fht_instance.transform(otf);
   double [][][]cmplx=  FHT2FFTHalf (otf, size);
   cmplx= amplPhase(cmplx); 
   return cmplx;
}





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

double [] convolveByFHT(double [] pixels, // rectangular image array to be convolved
                               int width, // image width,
                      double []fhtKernel, // FHT transform of the DOUBLE-SIZE convolution kernel)
                    double []slidingMask) { // sliding mask (made of 0.25(cos(X)+1)*(cos(y)+1), or null (will be calculated)
    if ((pixels==null) || (fhtKernel==null)) return null;
    int length=fhtKernel.length;
    int size=(int) Math.sqrt(length);
    int size_q=size/4;
    int size_h=size/2;
    int size_3q=3*size/4;
    int height=pixels.length/width;
    double [] result = new double [width*height];
    double [] buffer = new double [length];
    double [] product=null;
    int i,j,i1,j1;
    for (i=0;i<result.length;i++) result[i]=0.0;
    if ((slidingMask==null) || (slidingMask.length!=(length/4))) slidingMask=getSlidingMask(size_h);
    int tileY,tileX;
    int nTileX=width/size_q+1;
    if (width>((nTileX-1)*size_q)) nTileX++;
    if (nTileX<2) nTileX=2;
    int nTileY=height/size_q+1;
    if (height>((nTileY-1)*size_q))nTileY++;
    if (nTileY<2) nTileY=2;
    if (DEBUG_LEVEL>3) SDFA_instance.showArrays(fhtKernel, size, size, "fhtKernel");
    if (DEBUG_LEVEL>3) SDFA_instance.showArrays(slidingMask, size_h, size_h, "slidingMask");
    for (tileY=0;tileY<nTileY;tileY++) for (tileX=0;tileX<nTileX;tileX++) {
/* initilaize buffer */
      for (i=0;i<length;i++) buffer[i]=0.0;
      for (i=size_q;i<size_3q;i++) {
         i1=i+size_q*(tileY-2);
         if ((i1<0) || (i1>=height)) continue;
         for (j=size_q;j<size_3q;j++) {
           j1=j+size_q*(tileX-2);
           if ((j1<0) || (j1>=width)) continue;
           buffer[i*size+j]=pixels[i1*width+j1]*slidingMask[(i-size_q)*size_h+(j-size_q)];
         }
      }
      if ((DEBUG_LEVEL>3)&&(tileY==0)&&(tileX==0)) SDFA_instance.showArrays(buffer, size, size, "buffer");
/* make direct FHT of buffer */
// Swapping quadrants, so the center will be 0,0
      fht_instance.swapQuadrants(buffer);
// get to frequency domain
      fht_instance.transform(buffer);
      product=     fht_instance.multiply(buffer, fhtKernel, false);
      fht_instance.inverseTransform(product);
      fht_instance.swapQuadrants(product);
      if ((DEBUG_LEVEL>3)&&(tileY==0)&&(tileX==0)) SDFA_instance.showArrays(product, size, size, "product");

/* Add buffer data to the result array */
      for (i=0;i<size;i++) {
         i1=i+size_q*(tileY-2);
         if ((i1<0) || (i1>=height)) continue;
         for (j=0;j<size;j++) {
           j1=j+size_q*(tileX-2);
           if ((j1<0) || (j1>=width)) continue;
           result[i1*width+j1]+=product[i*size+j];
         }
      }
      if ((DEBUG_LEVEL>3)&&(tileY==0)&&(tileX==0)) SDFA_instance.showArrays(result, width, height, "result00");

    }
    return result;
}

/*
      if (DEBUG_LEVEL>1) {
        for (i=0;i<input_bayer.length;i++) if (!colorsToCorrect[i]) input_bayer[i]=null;
        SDFA_instance.showArrays(smoothInputBayers, FFTSize*subDivFreq*PSF_subpixel, FFTSize*subDivFreq*PSF_subpixel, imp_sel.getTitle()+"-orig");

*/

/* reject Bayer and oversampling aliases by masking out certain frequency components. mask has 1.0 at frequency 0,0 */


public double [] rejectByMask (double [] pixels,  // square input data
                               double [] mask, // mask to multiply FHT
                               boolean centered){ // image is sentered around the center of the square (use swapQuadrants)
   if ((pixels==null) || (mask==null)) return null;
   int length = pixels.length;
   if (length != mask.length) return null;
   int size = (int) Math.sqrt(length);
   double [] result=new double [length];
   int i;
   float [] floatPixels = new float [length];
   for (i=0; i<length;i++) floatPixels[i]= (float) pixels[i];
   ImageProcessor ip = new FloatProcessor(size,size);
   ip.setPixels(floatPixels);
   FHT fht =  new FHT(ip);
// Swapping quadrants, so the center will be 0,0
   if (centered) fht.swapQuadrants();
// get to frequency domain
   fht.transform();
   floatPixels=(float []) fht.getPixels();
   for (i=0;i<length;i++) {
      floatPixels[i]= (float) (floatPixels[i]*mask[i]);
   }
   ip.setPixels(floatPixels);
   fht.inverseTransform();
   if (centered) fht.swapQuadrants();
   floatPixels=(float []) fht.getPixels();
   for (i=0; i<length;i++) result[i]=floatPixels[i];
   return result;
 }

public double [] createAliasReject (int    size,  // size of the mask
                                boolean checker,  // checkerboard pattern in the source file (use when filtering)
                                 int oversample,  // measured array is sampled at 1/oversample frequency than model (will add more parameters later)
                             double gaussSigma){ // width of rejection areas on the spectrum (the smaller, the sharper rejection)
   int length = size*size;
   double [] gauss=new double [size];
   double [] mask=new double [length];
   double k=0.5/gaussSigma/gaussSigma;
   int maxI=(int) gaussSigma+1;
   int cloneStep=size/oversample;
   int i,j;
   for (i=0; i<=size/2;i++) {
     if (i>maxI) gauss[i]=0.0;
     else gauss[i]=Math.exp(-(k*i*i));
     if (i>0) gauss[size-i]=gauss[i];
   }
   
   
   for (i=0;i<length;i++) mask[i]=0.0;
/* Clone gaussians */
   double d=0.0;
   int cloneNy,cloneNx,cloneY,cloneX;
   for (i=0;i<size;i++) if (gauss[i]!=0.0) for (j=0;j<size;j++) if (gauss[j]!=0.0){
     d=gauss[i]*gauss[j];
     for (cloneNy=0;cloneNy<oversample;cloneNy++) for (cloneNx=0;cloneNx<oversample;cloneNx++)
       if (((cloneNy!=0) || (cloneNx!=0)) && // not a zero point
           (!checker ||                      // use all if it is not a checkerboard pattren
          (((cloneNx ^ cloneNy) & 1)==0) )) { // remove clones in a checker pattern
             cloneY=(i+cloneNy*cloneStep) % size;
             cloneX=(j+cloneNx*cloneStep) % size;
       mask[cloneY*size+cloneX]+=d;
     }
   }
/* reverse and make DC coefficient ==1.0 */
  double max=mask[0];
  double min=mask[0]; // mask[0] _is_ min 
  for (i=1;i<length;i++)  if (max<mask[i]) max = mask[i];
  k=1.0/(max-min);
/*
   if (DEBUG_LEVEL>1)  {
      System.out.println("createAliasReject("+size+", "+checker+", "+oversample+", "+gaussSigma+")");
      System.out.println("mask[0]="+mask[0]+", k="+k);
    }
*/
  for (i=0;i<length;i++) {
       mask[i]=1.0-k*(mask[i]-min);
  }
  return mask;
 }


//-------------------------------------------------

double [] binPSF(double [] pixels,
                 double [][] g,
                 int outSize,
                 int      decimate,     // sub-pixel decimation 
                 double minContrast,
                 double [] centerXY,    // coordinates (x,y) of the center point (will be alway subtracted)
                 double[] symmXY,       // coordinates (x,y) of the center of symmetry (to combine with 180 if enabled by symm180)
                 int pass,              // mostly for debug purposes
                 String title,
                 boolean debug  ) {
    int multiple=2;         // 0 - use each pixel once, 1 - add first negatives (4), 2 - second positives()4)
    int pixelSize=(int) Math.sqrt(pixels.length);
    int halfOutSize=outSize/2;
    int indx,i,j,outIndex,ix,iy;
    double x,y,xc,yc,uc,vc,u,v,p,q,d, du, dv, dp,dq, xr,yr, overThreshold;
    int np,nq;
    int PSF_sign=1;
    double [] contrastCache=new double[pixelSize*pixelSize];
    double [] debugPixels=null;
    if (debug)  debugPixels=new double[pixelSize*pixelSize];

    double det_g=g[0][0]*g[1][1]-g[0][1]*g[1][0];
    double [][] xy2uv= {{-2.0*g[0][1]/det_g,  2.0*g[0][0]/det_g},
                        {-2.0*g[1][1]/det_g,  2.0*g[1][0]/det_g}};
    double [][] uv2xy= matrix2x2_scale(matrix2x2_invert(xy2uv),2); // real pixels are twice
    double [] pixelsPSF       =new double [outSize*outSize];  
    int    [] pixelsPSFCount  =new int    [outSize*outSize];
    double [] pixelsPSFWeight =new double [outSize*outSize];  
    double [] center=centerXY;
    for (i=0;i<contrastCache.length;i++) {
      contrastCache[i]=-1.0;
    }

//    double overThreshold;
//    double threshold=minContrast*contrastAtXYDebug(1, pixels, pixelSize, 0.0, 0.0,  g, contrastCache);
    double threshold=minContrast*contrastAtXY(1, pixels, pixelSize, 0.0, 0.0,  g, contrastCache);



   if (debug)  {
      System.out.println("binPSF title="+title+" g[0][0]="+IJ.d2s(g[0][0],4)+" g[0][1]="+IJ.d2s(g[0][1],4));
      System.out.println("binPSF title="+title+" g[1][0]="+IJ.d2s(g[1][0],4)+" g[1][1]="+IJ.d2s(g[1][1],4));
      System.out.println("  center[0]="+center[0]+"  center[1]="+center[1]);
      System.out.println("  decimate="+decimate+"  threshold="+threshold);
    }




    if (center==null) {
      center = new double[2];
      center[0]=0.0;
      center[1]=0.0;
    }
    for (i=0;i<pixelsPSF.length;i++) {
      pixelsPSF[i]=0.0;
      pixelsPSFCount[i]=0;
      pixelsPSFWeight[i]=0.0;
    }

    for (indx=0;indx<pixels.length;indx++) {
      y= indx / pixelSize- pixelSize/2;
      x= indx % pixelSize- pixelSize/2;
      u= xy2uv[0][0]*x + xy2uv[0][1]*y;
      v= xy2uv[1][0]*x + xy2uv[1][1]*y;
      p=u+v;
      q=u-v;
      np=(int)Math.floor((1+p)/2);
      nq=(int)Math.floor((1+q)/2);
//      if (debug)  debugPixels[indx]=(int)Math.floor((1+q)/2);
/* see if the point is in the cell of positive or negative OTF instance */
      PSF_sign= (((np + nq) & 1)==0)?1:-1;
/* find x,y coordinates of the center of the cell */
      uc=0.5*(np+nq);
      vc=0.5*(np-nq);
//    xc=g[0][0]*uc + g[1][0]*vc;
//    yc=g[0][1]*uc + g[1][1]*vc;

      yc=-g[0][0]*uc - g[1][0]*vc;
      xc= g[0][1]*uc + g[1][1]*vc;


//if (debug) debugPixels[indx]=p/2-Math.round(p/2);

/* See if this cell has enough contrast */
      overThreshold=contrastAtXY(PSF_sign,pixels, pixelSize, xc,yc,  g, contrastCache);
//if (debug) debugPixels[indx]=overThreshold;
      if (overThreshold<threshold) {
         if (debug) debugPixels[indx]=0.0;
//         if (debug) debugPixels[indx]=yc;
      } else {
//     if (debug) debugPixels[indx]=yc;

/* Do binning itself here */
        d=PSF_sign*PSFAtXY(pixels, pixelSize, x,y);

/* map to the segment around 0,0 */        
        dp=p/2-Math.round(p/2);
        dq=q/2-Math.round(q/2);
/* dp, dq are between +/- 0.5 - use them for Hamming windowing -NOT HERE, moved later*/
/*
        if (useWindow) {
          d*= (0.54+0.46*Math.cos(dp*2.0*Math.PI));
          d*= (0.54+0.46*Math.cos(dq*2.0*Math.PI));
        }
*/
        du=(dp+dq)/2;
        dv=(dp-dq)/2;
//        if (debug) debugPixels[indx]=du;
//if (debug) debugPixels[indx]=uv2xy[0][0]*(0+du) + uv2xy[0][1]*(0+dv);

/* bin this point to the center and some (positive) duplicates if enabled */
        for (i=-(multiple/2); i<=(multiple/2); i++) for (j=-(multiple/2); j<=(multiple/2); j++) {
          xr= uv2xy[0][0]*(j+du) + uv2xy[0][1]*(i+dv);
          yr= uv2xy[1][0]*(j+du) + uv2xy[1][1]*(i+dv);
          xr= Math.round(decimate*xr-center[0]);
          yr= Math.round(decimate*yr-center[1]);

//          xr=Math.round((ignoreChromatic?(-symmXY[0]):0.0)+decimate*(g[0][0]*(j+du) + g[1][0]*(i+dv)));
//          yr=Math.round((ignoreChromatic?(-symmXY[1]):0.0)+decimate*(g[0][1]*(j+du) + g[1][1]*(i+dv)));
/* does it fit into output array ? */
          if ((yr>=-halfOutSize) && (yr<halfOutSize) && (xr>=-halfOutSize) && (xr<halfOutSize)) {
            outIndex=outSize*(outSize/2+ ((int) yr))+(outSize/2)+((int) xr);
            pixelsPSFCount[outIndex]++;
            pixelsPSF[outIndex]+=d*overThreshold;
//            pixelsPSF[outIndex]+=xr*overThreshold;
//            pixelsPSF[outIndex]+=overThreshold;
            pixelsPSFWeight[outIndex]+=overThreshold;
          }
        }
/* bin this to center-symmetrical point if enabled */
        if (symmXY!=null) {
          for (i=-(multiple/2); i<=(multiple/2); i++) for (j=-(multiple/2); j<=(multiple/2); j++) {
            xr= uv2xy[0][0]*(j+du) + uv2xy[0][1]*(i+dv);
            yr= uv2xy[1][0]*(j+du) + uv2xy[1][1]*(i+dv);
            xr= Math.round(symmXY[0]*2.0-decimate*xr-center[0]);
            yr= Math.round(symmXY[1]*2.0-decimate*yr-center[1]);
//            xr=Math.round((ignoreChromatic?1.0:2.0)*symmXY[0]-(decimate*(g[0][0]*(j+du) + g[1][0]*(i+dv))));
//            yr=Math.round((ignoreChromatic?1.0:2.0)*symmXY[1]-(decimate*(g[0][1]*(j+du) + g[1][1]*(i+dv))));
// does it fit into output array ?
            if ((yr>=-halfOutSize) && (yr<halfOutSize) && (xr>=-halfOutSize) && (xr<halfOutSize)) {
              outIndex=outSize*(outSize/2+ ((int) yr))+(outSize/2)+((int) xr);
              pixelsPSFCount[outIndex]++;
              pixelsPSF[outIndex]+=d*overThreshold;
              pixelsPSFWeight[outIndex]+=overThreshold;
            }
          }
        }
/* Now bin this point to the negative duplicates if enabled (debug feature). Normally it will be skipped */
        if (multiple>0) for (i=-((multiple+1)/2); i<((multiple+1)/2); i++) for (j=-((multiple+1)/2); j<((multiple+1)/2); j++) {
          xr= uv2xy[0][0]*(j+du+0.5) + uv2xy[0][1]*(i+dv+0.5);
          yr= uv2xy[1][0]*(j+du+0.5) + uv2xy[1][1]*(i+dv+0.5);
          xr= Math.round(decimate*xr-center[0]);
          yr= Math.round(decimate*yr-center[1]);
//          xc=Math.round((ignoreChromatic?(-symmXY[0]):0.0)+decimate*(g[0][0]*(j+du+0.5) + g[1][0]*(i+dv+0.5)));
//          yc=Math.round((ignoreChromatic?(-symmXY[1]):0.0)+decimate*(g[0][1]*(j+du+0.5) + g[1][1]*(i+dv+0.5)));
// does it fit into output array ?
          if ((yr>=-halfOutSize) && (yr<halfOutSize) && (xr>=-halfOutSize) && (xr<halfOutSize)) {
            outIndex=outSize*(outSize/2+ ((int) yr))+(outSize/2)+((int) xr);
            pixelsPSFCount[outIndex]++;
            pixelsPSF[outIndex]-=d*overThreshold;
            pixelsPSFWeight[outIndex]+=overThreshold;
          }
        }
/* bin this to center-symmetrical point if enabled */
/* Now bin this point to the negative duplicates if enabled (debug feature). Normally it will be skipped */
        if (symmXY!=null) {
          if (multiple>0) for (i=-((multiple+1)/2); i<((multiple+1)/2); i++) for (j=-((multiple+1)/2); j<((multiple+1)/2); j++) {
            xr= uv2xy[0][0]*(j+du+0.5) + uv2xy[0][1]*(i+dv+0.5);
            yr= uv2xy[1][0]*(j+du+0.5) + uv2xy[1][1]*(i+dv+0.5);
            xr= Math.round(symmXY[0]*2.0-decimate*xr-center[0]);
            yr= Math.round(symmXY[1]*2.0-decimate*yr-center[1]);

//            xc=Math.round((ignoreChromatic?1.0:2.0)*symmXY[0]-(decimate*(g[0][0]*(j+du+0.5) + g[1][0]*(i+dv+0.5))));
//            yc=Math.round((ignoreChromatic?1.0:2.0)*symmXY[1]-(decimate*(g[0][1]*(j+du+0.5) + g[1][1]*(i+dv+0.5))));
// does it fit into output array ?
            if ((yr>=-halfOutSize) && (yr<halfOutSize) && (xr>=-halfOutSize) && (xr<halfOutSize)) {
              outIndex=outSize*(outSize/2+ ((int) yr))+(outSize/2)+((int) xr);
              pixelsPSFCount[outIndex]++;
              pixelsPSF[outIndex]+=d*overThreshold;
              pixelsPSFWeight[outIndex]+=overThreshold;
            }
          }
        }
      }
    }


    for (i=0;i<pixelsPSF.length;i++) {
      if (pixelsPSFWeight[i]>0.0) pixelsPSF[i]/=pixelsPSFWeight[i];
    }
/* Interpolate  missing points (pixelsPSFCount[i]==0) */

    for (i=0;i<pixelsPSF.length;i++) if (pixelsPSFWeight[i]==0.0){
      iy=i/outSize;
      ix=i%outSize;
      if ((ix>0)&&(ix<(outSize-1))&&(iy>0)&&(iy<(outSize-1))) {
        if ((pixelsPSFWeight[(iy-1)*outSize+ix  ]>0.0) &&
            (pixelsPSFWeight[(iy+1)*outSize+ix  ]>0.0) &&
            (pixelsPSFWeight[(iy  )*outSize+ix-1]>0.0) &&
            (pixelsPSFWeight[(iy  )*outSize+ix+1]>0.0)) {
           if (DEBUG_LEVEL>5) System.out.println("Interpolating missing OTF point at x="+ix+" y="+iy);
           pixelsPSF[i]=
            0.25*(pixelsPSF[(iy-1)*outSize+ix  ]+
                  pixelsPSF[(iy+1)*outSize+ix  ]+
                  pixelsPSF[(iy  )*outSize+ix-1]+
                  pixelsPSF[(iy  )*outSize+ix+1]);
        }
      }
    }
/* optionally show original array with masked out low-contrast cells */
    if ((DEBUG_LEVEL>1) && (pass==1)) {
      float [] floatPixelsPSFDbg=new float[pixelsPSF.length];
      for (j=0;j<pixelsPSF.length;j++) floatPixelsPSFDbg[j]=(float)pixelsPSF[j];
      ImageProcessor ip_PSFdbg=new FloatProcessor(outSize,outSize);
      ip_PSFdbg.setPixels(floatPixelsPSFDbg);
      ip_PSFdbg.resetMinAndMax();
      ImagePlus imp_PSFdbg=  new ImagePlus(title+"_Used-PSF", ip_PSFdbg);
      imp_PSFdbg.show();
    }

    if (debug) {
      float [] floatPixelsDbg=new float[debugPixels.length];
      for (j=0;j<debugPixels.length;j++) floatPixelsDbg[j]=(float)debugPixels[j];
      ImageProcessor ip_dbg=new FloatProcessor(pixelSize,pixelSize);
      ip_dbg.setPixels(floatPixelsDbg);
      ip_dbg.resetMinAndMax();
      ImagePlus imp_dbg=  new ImagePlus(title+"_mask_PSF", ip_dbg);
      imp_dbg.show();

      float [] floatPixelsPSFCount=new float [pixelsPSF.length];
      for (j=0;j<floatPixelsPSFCount.length;j++) floatPixelsPSFCount[j]=(float)pixelsPSFCount[j];
      ImageProcessor ip_count=new FloatProcessor(outSize,outSize);
      ip_count.setPixels(floatPixelsPSFCount);
      ip_count.resetMinAndMax();
      ImagePlus imp_count=  new ImagePlus(title+"_PSF_bin_count", ip_count);
      imp_count.show();

      float [] floatPixelsPSFWeight=new float [pixelsPSF.length];
      for (j=0;j<floatPixelsPSFWeight.length;j++) floatPixelsPSFWeight[j]=(float)pixelsPSFWeight[j];
      ImageProcessor ip_weight=new FloatProcessor(outSize,outSize);
      ip_weight.setPixels(floatPixelsPSFWeight);
      ip_weight.resetMinAndMax();
      ImagePlus imp_weight=  new ImagePlus(title+"_PSF_bin_weight", ip_weight);
      imp_weight.show();

      float [] floatContrastCache=new float [contrastCache.length];
      for (j=0;j<floatContrastCache.length;j++) floatContrastCache[j]=(float)((contrastCache[j]>=0.0)?contrastCache[j]:-0.00001);
      ImageProcessor ip_ContrastCache=new FloatProcessor(pixelSize,pixelSize);
      ip_ContrastCache.setPixels(floatContrastCache);
      ip_ContrastCache.resetMinAndMax();
      ImagePlus imp_ContrastCache=  new ImagePlus(title+"_ContrastCache", ip_ContrastCache);
      imp_ContrastCache.show();

// contrastCache
    }
    return pixelsPSF;
  }

//====================






double [] combinePSF (double []pixels,         // Square array of pixels with multiple repeated PSF (alternating sign)
                      double [][] wVectors,    // two wave vectors, lengths in cycles/pixel (pixels match pixel array)
                      int      decimate,       // sub-pixel decimation 
                      double   minContrast,    // minimal instance contrast to use in binning
                      double   windowFrac,     // reduce the PSF cell size to this part of the area connecting first negative clones
                      boolean  useWindow,
                      boolean  symm180,
                      boolean  ignoreChromatic,
                      double[] centerXY,         // coordinates (x,y) of the center point (will update if ignoreChromatic is true)
                      double [] centroid_xy,    // RETURNS centroid of the result array (should be small) if ignoreChromatic is true
                      boolean enableModelSubtract,  // generate/subtract gaussian models (not needed if no overlap between pos/neg clones)
                      double   smoothSeparate, // low-pass filter multiple opposite-sign PSF instaces for separation (width relative to distance to the opposite sign PSF)
                      double   thresholdSeparate,  // threshold for locating zero-crossing
                      double   topCenter,     // consider only points above that fraction of the local max to find centroid
                      boolean  removeNegtative,  // remove PSF negative values  when separating composite PSF (will need low-pass filtering)
                      double   sigmaToRadius,   // 0.4;  variable-sigma blurring to reduce high frequencies more for the pixels farther from the PSF center
                      double   wings_energy, // fraction of energy in the pixels to be used
                      double   wings_ellipse_scale,
                      String title,
                      boolean debug)
  {
    if (pixels==null) return null;
    //    double [] contrastCache=new double[pixelSize*pixelSize];
    int i,j;

    if (DEBUG_LEVEL>2) {
      System.out.println("combinePSF title="+title+" wV[0][0]="+IJ.d2s(wVectors[0][0],4)+" wV[0][1]="+IJ.d2s(wVectors[0][1],4));
      System.out.println("combinePSF title="+title+" wV[1][0]="+IJ.d2s(wVectors[1][0],4)+" wV[1][1]="+IJ.d2s(wVectors[1][1],4));
    }

/* vectors perpendicular to the checkerboard edges, lengths equal to the periods */
    double [][] f= {{wVectors[0][0]/(wVectors[0][0]*wVectors[0][0]+wVectors[0][1]*wVectors[0][1]),
                     wVectors[0][1]/(wVectors[0][0]*wVectors[0][0]+wVectors[0][1]*wVectors[0][1])},
                    {wVectors[1][0]/(wVectors[1][0]*wVectors[1][0]+wVectors[1][1]*wVectors[1][1]),
                     wVectors[1][1]/(wVectors[1][0]*wVectors[1][0]+wVectors[1][1]*wVectors[1][1])}};
    if (DEBUG_LEVEL>2) {
      System.out.println("combinePSF title="+title+" f[0][0]="+IJ.d2s(f[0][0],4)+" f[0][1]="+IJ.d2s(f[0][1],4));
      System.out.println("combinePSF title="+title+" f[1][0]="+IJ.d2s(f[1][0],4)+" f[1][1]="+IJ.d2s(f[1][1],4));
    }

/* vectors parallel to checkerboard edges, lenghs equal to the period along those lines */
    double l2f1=   f[0][0]*f[0][0]+f[0][1]*f[0][1];
    double l2f2=   f[1][0]*f[1][0]+f[1][1]*f[1][1];
    double pf1f2  =f[0][1]*f[1][0]-f[1][1]*f[0][0];
    double [][]g0= {{f[0][1]*l2f2/pf1f2,  -f[0][0]*l2f2/pf1f2},
                   {f[1][1]*l2f1/pf1f2,  -f[1][0]*l2f1/pf1f2}};
    if (DEBUG_LEVEL>2) {
      System.out.println("combinePSF title="+title+" g0[0][0]="+IJ.d2s(g0[0][0],4)+" g[0][1]="+IJ.d2s(g0[0][1],4));
      System.out.println("combinePSF title="+title+" g0[1][0]="+IJ.d2s(g0[1][0],4)+" g[1][1]="+IJ.d2s(g0[1][1],4));
    }
/* calculate vectors connecting centers of the "positive" PSF copies */

    double [][] g= {{0.5*(g0[0][0]+g0[1][0]), 0.5*(g0[0][1]+g0[1][1])},
                    {0.5*(g0[0][0]-g0[1][0]), 0.5*(g0[0][1]-g0[1][1])}};

    if (DEBUG_LEVEL>2) {
      System.out.println("combinePSF title="+title+" g[0][0]="+IJ.d2s(g[0][0],4)+" g[0][1]="+IJ.d2s(g[0][1],4));
      System.out.println("combinePSF title="+title+" g[1][0]="+IJ.d2s(g[1][0],4)+" g[1][1]="+IJ.d2s(g[1][1],4));
    }

/*   double det_g=g[0][0]*g[1][1]-g[0][1]*g[1][0];
   double [][] xy2uv= {{-2.0*g[0][1]/det_g,  2.0*g[0][0]/det_g},
                        {-2.0*g[1][1]/det_g,  2.0*g[1][0]/det_g}};*/

/// =================

/* calculate outSize to be able to use FFT here */
    double sizeNegatives= Math.max(Math.max(Math.abs(g[0][0]+ g[1][0]),Math.abs(g[0][1]+ g[1][1])),
                                   Math.max(Math.abs(g[0][0]- g[1][0]),Math.abs(g[0][1]- g[1][1])));
    double scaleSize=2.5; /// Will include next positive centers and overlap
    int outSize;
    for (outSize=8;outSize<scaleSize*sizeNegatives; outSize<<=1);
    int halfOutSize=outSize/2;
    if (DEBUG_LEVEL>2) {
      System.out.println("sizeNegatives="+sizeNegatives+ " scaled="+ (scaleSize*sizeNegatives)+" outSize="+outSize+" halfOutSize="+halfOutSize);
    }

    double [] pixelsPSF= binPSF(pixels,
                                     g,
                               outSize,
                              decimate,  // sub-pixel decimation (now not used as the input pixels array already has subpixel resolution)
                           minContrast,  // minimal contrast of PSF clones
                              centerXY,  //  coordinates (x,y) of the center point
                                  null,  // coordinates of the center of symmetry - not applicable
                                     1, // pass 1
                                 title,
                                 debug);
//                                 true);

    double distToNegativeClones=0.5*Math.sqrt(((g[0][0]+g[1][0])*(g[0][0]+g[1][0])+
                                               (g[0][1]+g[1][1])*(g[0][1]+g[1][1])+
                                               (g[0][0]-g[1][0])*(g[0][0]-g[1][0])+
                                               (g[0][1]-g[1][1])*(g[0][1]-g[1][1]))/2.0);
    if (DEBUG_LEVEL>2) {
      System.out.println("distToNegativeClones="+distToNegativeClones+ " gaussWidth="+ distToNegativeClones*smoothSeparate);
    }
    double smoothSigma=distToNegativeClones*smoothSeparate;

    double [] smoothPixelsPSF= lowPassGauss(pixelsPSF, smoothSigma, true);
/* find amplitude of smoothed pixel array */
    double smoothMin=0.0;
    double smoothMax=0.0;
    for (i=0;i<smoothPixelsPSF.length;i++) {
      if      (smoothPixelsPSF[i] > smoothMax) smoothMax=smoothPixelsPSF[i];
      else if (smoothPixelsPSF[i] < smoothMin) smoothMin=smoothPixelsPSF[i];
    }
    double smoothAmplitude=(smoothMax-smoothMin)/2;


    int [][]  clusterMask = findClusterOnPSF(smoothPixelsPSF, // PSF function, square array (use smooth array)
                                                  -topCenter, // fraction of energy in the pixels to be used (or minimal level if it is negative)
                                                   outSize/2,  // location of a start point, x-coordinate
                                                   outSize/2,  // location of a start point, y-coordinate
                                                       title);
    double [] centroidXY=       calcCentroidFromCenter(pixelsPSF, // use original array (mask from the smoothed one)
// --centroidXY is in function call arguments
//    centroidXY=            calcCentroidFromCenter(pixelsPSF, // use original array (mask from the smoothed one)
                                                clusterMask, // integer mask -0 - don't use this pixel, 1 - use it
                                                  topCenter);// subtract level below topCenter*max
    double [] centroidXY_smooth=calcCentroidFromCenter(smoothPixelsPSF, // use smooth - not final, just for clones rejection
                                                clusterMask, // integer mask -0 - don't use this pixel, 1 - use it
                                                  topCenter);// subtract level below topCenter*max

    if (DEBUG_LEVEL>2) System.out.println("Centroid after first binPSF: x="+IJ.d2s(centroidXY[0],3)+" y="+IJ.d2s(centroidXY[1],3)+" center was at x="+IJ.d2s(centerXY[0],3)+" y="+IJ.d2s(centerXY[1],3));

/*
    pixelsPSF= binPSF(   pixels,
                              g,
                        outSize,
                       decimate,  // sub-pixel decimation (now not used as the input pixels array already has subpixel resolution)
                    minContrast,  // minimal contrast of PSF clones
     ignoreChromatic? centroidXY:null,
      symm180?centroidXY:null,
                              2,
                          title,
                          debug);
*/
/* Re-bin results with the new center if ignoreChromatic is true, update centerXY[](shift of the result PSF array) and centroidXY[] (center of the optionally shifted PDF array) */
    if (ignoreChromatic) {
      if (centerXY!=null) {
        centerXY[0]+=centroidXY[0];
        centerXY[1]+=centroidXY[1];
      }
      pixelsPSF= binPSF(   pixels,
                                g,
                          outSize,
                         decimate,  // sub-pixel decimation (now not used as the input pixels array already has subpixel resolution)
                      minContrast,  // minimal contrast of PSF clones
                         centerXY,  // now includes centroid from the pass 1
//       ignoreChromatic? centroidXY:null,
       symm180?centroidXY:null,
                                2, // pass2
                            title,
                            debug);

/*  recalculate centroids  */
      smoothPixelsPSF= lowPassGauss(pixelsPSF, smoothSigma, true);
      smoothMin=0.0;
      smoothMax=0.0;
      for (i=0;i<smoothPixelsPSF.length;i++) {
        if      (smoothPixelsPSF[i] > smoothMax) smoothMax=smoothPixelsPSF[i];
        else if (smoothPixelsPSF[i] < smoothMin) smoothMin=smoothPixelsPSF[i];
      }
      smoothAmplitude=(smoothMax-smoothMin)/2;
      clusterMask = findClusterOnPSF(smoothPixelsPSF, // PSF function, square array (use smooth array)
                                          -topCenter, // fraction of energy in the pixels to be used (or minimal level if it is negative)
                                           outSize/2,  // location of a start point, x-coordinate
                                           outSize/2,  // location of a start point, y-coordinate
                                               title);
      centroidXY= calcCentroidFromCenter(pixelsPSF, // use original array (mask from the smoothed one)
                                         clusterMask, // integer mask -0 - don't use this pixel, 1 - use it
                                           topCenter);// subtract level below topCenter*max
      centroidXY_smooth=calcCentroidFromCenter(smoothPixelsPSF, // use smooth - not final, just for clones rejection
                                                     clusterMask, // integer mask -0 - don't use this pixel, 1 - use it
                                                       topCenter);// subtract level below topCenter*max
      if (DEBUG_LEVEL>2) System.out.println("Centroid after second binPSF: x="+IJ.d2s(centroidXY[0],3)+" y="+IJ.d2s(centroidXY[1],3)+" center was at x="+IJ.d2s(centerXY[0],3)+" y="+IJ.d2s(centerXY[1],3));

    }
/* compensate center point and/or add centre-symmetrical points if enabled */
    double [] rejectedClonesPixels=null;
    double [] smoothRejectedClonesPixels=null;
    double [][] modelPSFVectors={{0.5*(g[0][0]+g[1][0]),0.5*(g[0][1]+g[1][1])},
                                 {0.5*(g[0][0]-g[1][0]),0.5*(g[0][1]-g[1][1])}};
    double [] modelPixels=null;

    if (enableModelSubtract) {
/* This method only can work between opposite sign clones, so when PSF is elliptical with the big axis in the middle  between negative clones,
    it will remain unnoticed in the model, model will have two approximately the same sigmas.
    Such example - tilt09.tif (low right corner)  */

      double [][]section = new double [2][];
      int sectionLength=50;
      section[0]=  createSection (smoothPixelsPSF, // square pixel array
                                  sectionLength, // length of connecting section line, pixels
                                centroidXY_smooth[0], // x of the center maximum (positive), relative to the center of array
                                centroidXY_smooth[1],
                          0.5*(g[0][0]+g[1][0]),
                          0.5*(g[0][1]+g[1][1]));
      section[1]=  createSection (smoothPixelsPSF, // square pixel array
                                  sectionLength, // length of connecting section line, pixels
                                centroidXY_smooth[0], // x of the center maximum (positive), relative to the center of array
                                centroidXY_smooth[1],
                          0.5*(g[0][0]-g[1][0]),
                          0.5*(g[0][1]-g[1][1]));

      if (DEBUG_LEVEL>2) for (i=0;i<section.length;i++) {
        System.out.println("Section="+i+ " section length="+ section[0].length+" =========================================");
        for (j=0;j<section[0].length;j++)System.out.println(section[i][j]);
      }
      double [][] gaussApprox=new double [2][];
      gaussApprox[0]= approximateOppositeGauss (section[0], // Array of pixes along a section, connecting positive and negative maximums
             0.5*Math.sqrt((g[0][0]+g[1][0])* (g[0][0]+g[1][0])+(g[0][1]+g[1][1])*(g[0][1]+g[1][1])),
                         smoothSigma  // width of the "smoothing" filter applied to the data before sestion (in pixels)
                                                           );
      gaussApprox[1]= approximateOppositeGauss (section[1], // Array of pixes along a section, connecting positive and negative maximums
             0.5*Math.sqrt((g[0][0]-g[1][0])* (g[0][0]-g[1][0])+(g[0][1]-g[1][1])*(g[0][1]-g[1][1])),
                         smoothSigma  // width of the "smoothing" filter applied to the data before sestion (in pixels)
                                                         );
/* Now we have Gauss approximation along vectors connecting positive PSF with negative clones (4 values)),
    we can create model PSF for all clones ( 4 negatives and 4 positives) around the center one and
    subtract them from (not smoothed!) composite pixelsPSF[]). During Gauss approximation we paid most
    attention to wide part of the model - it is what influences most the subtraction result, overlapping
    with the center (main) PSF instance. We'll smooth with variable filter the result anyway, later */
//    double [][] modelPSFVectors={{0.5*(g[0][0]+g[1][0]),0.5*(g[0][1]+g[1][1])},
//                                 {0.5*(g[0][0]-g[1][0]),0.5*(g[0][1]-g[1][1])}};
    modelPixels=new double [pixelsPSF.length];
    double scaleModelPSF = 1.0; /* calculate later */
    int woi_width=outSize; /* reduce later */
    int    signModelPSF;
    double [] averageSigmas={Math.sqrt(0.5*(gaussApprox[0][0]*gaussApprox[0][0]+gaussApprox[0][1]*gaussApprox[0][1])),
                             Math.sqrt(0.5*(gaussApprox[1][0]*gaussApprox[1][0]+gaussApprox[1][1]*gaussApprox[1][1]))};
    for (i=0;i<modelPixels.length; i++) modelPixels[i]=0.0;
/* See if we got reasonable sigmas for 1ll 4 directions, otherwise we have very sharp PSF and will skip the next step */
    double sigmaThreshold=thresholdSeparate*smoothSigma;
    if ((gaussApprox[0][0]<sigmaThreshold) || (gaussApprox[0][1]<sigmaThreshold) || (gaussApprox[1][0]<sigmaThreshold) || (gaussApprox[0][1]<sigmaThreshold)) {
      if (DEBUG_LEVEL>2) {
        System.out.println("Approximated PSF model is too sharp, compared to the smoothing value, skipping clones compensation");
        System.out.println("sigmaThreshold="+sigmaThreshold+
                         " gaussApprox[0][0]="+gaussApprox[0][0]+
                         " gaussApprox[0][1]="+gaussApprox[0][1]+
                         " gaussApprox[1][0]="+gaussApprox[1][0]+
                         " gaussApprox[1][1]="+gaussApprox[1][1]);
      }
    } else {
      scaleModelPSF=smoothAmplitude*Math.sqrt(((averageSigmas[0]*averageSigmas[0]+smoothSigma*smoothSigma)/(averageSigmas[0]*averageSigmas[0])) *
                                              ((averageSigmas[1]*averageSigmas[1]+smoothSigma*smoothSigma)/(averageSigmas[1]*averageSigmas[1])));

      if (DEBUG_LEVEL>2) {
        System.out.println("sigmaThreshold="+sigmaThreshold+
                         " gaussApprox[0][0]="+gaussApprox[0][0]+
                         " gaussApprox[0][1]="+gaussApprox[0][1]+
                         " gaussApprox[1][0]="+gaussApprox[1][0]+
                         " gaussApprox[1][1]="+gaussApprox[1][1]+
                         " scaleModelPSF="+scaleModelPSF);
      }


/* Seems it is too much - with that the value between clones is exactly zero, trying scale - either 0.5 or sqrt(0.5) - maybe there is an error? */
//    scaleModelPSF*=Math.sqrt(0.5);
      scaleModelPSF*=0.5;
      if (DEBUG_LEVEL>2) {
        System.out.println("smoothAmplitude="+smoothAmplitude+" scaleModelPSF="+scaleModelPSF);
      }
      for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) if ((i!=0) || (j!=0)) {
        signModelPSF= (((i^j) & 1)!=0)?-1:1;
        addModelPSF(modelPixels, // square pixel array where the model PSF is added
            scaleModelPSF*signModelPSF, // scale model PSF (it has 1.0 in the center) before additing to pixels [] array
                centroidXY_smooth[0]+modelPSFVectors[0][0]*i+modelPSFVectors[1][0]*j, // model PSF center X-coordinate (in pixels[] units, from the center of the array )
                centroidXY_smooth[1]+modelPSFVectors[0][1]*i+modelPSFVectors[1][1]*j, // same for Y
                                     0, // window of interest in pixels[] array - do not generate data outside it
                                     0, // 
                             woi_width, //
                             woi_width, //
                       modelPSFVectors, // vectors that connect center of PSF with two oppositre sign clones
                        gaussApprox[0], // Gauss widths along vectors[0][] in p[ositive [0] and negative [1] directions
                        gaussApprox[1]); // Gauss widths along vectors[1][] in p[ositive [0] and negative [1] directions
      }
    }
/* Subtract clone models from original (later can be done in a single pass using addModelPSF() with original and reversed sign) */
    rejectedClonesPixels=new double [pixelsPSF.length];
    for (j=0;j<modelPixels.length;j++) rejectedClonesPixels[j]=pixelsPSF[j]-modelPixels[j];
/* Limit PSF to positive values only */
    if (removeNegtative) for (j=0;j<modelPixels.length;j++) if (rejectedClonesPixels[j]<0.0) rejectedClonesPixels[j]=0.0;

    smoothRejectedClonesPixels= lowPassGauss(rejectedClonesPixels, smoothSigma, true); /* not yet used */

/* Recalulate center */

    clusterMask = findClusterOnPSF(rejectedClonesPixels, // PSF function, square array (use smooth array)
                                             -topCenter, // fraction of energy in the pixels to be used (or minimal level if it is negative)
                                              outSize/2,  // location of a start point, x-coordinate
                                              outSize/2,  // location of a start point, y-coordinate
                                                 title);
    centroidXY=       calcCentroidFromCenter(rejectedClonesPixels, // use original array (mask from the smoothed one)
                                                      clusterMask, // integer mask -0 - don't use this pixel, 1 - use it
                                                        topCenter);// subtract level below topCenter*max
  } else {
     rejectedClonesPixels=pixelsPSF; // Maybe fo the opposite?
  }
  maskClonesPSF(rejectedClonesPixels, // square pixel array where the model PSF is added
                          windowFrac, // multiply window by this value
                     centroidXY[0], // Center of the remaining single PSF
                     centroidXY[1], // same for Y
                    modelPSFVectors, // vectors that connect center of PSF with two oppositre sign clones
                          useWindow);  // use Hamming window, if false - just cut sharp

  if (wings_energy>0.0) {
     rejectedClonesPixels=cutPSFWings (rejectedClonesPixels, // direct PSF function, square array, may be proportionally larger than reversed
                            wings_energy, // fraction of energy in the pixels to be used
                            wings_ellipse_scale,
                            0.003, // wings_min_mask_threshold, // zero output element if elliptical Gauss mask is below this threshold
                            title+"-w");
  }



  double [] sigmas=createSigmasRadius(rejectedClonesPixels, // input square pixel array, preferrably having many exact zeros (they will be skipped)
                                             sigmaToRadius, // sigma is proportional to the distance from the center
                                           centroidXY[0], // model PSF center X-coordinate (in pixels[] units, from the center of the array )
                                           centroidXY[1], // same for Y
                                                         0, // int WOICenterX, // window of interest in pixels[] array - do not generate data outside it
                                                         0, // int WOICenterY, // 
                                                   outSize, //int WOIWidth, reduce later
                                                   outSize); //int WOIHeight)


  double max1=0;
//  for (i=0;i<rejectedClonesPixels.length;i++) if (rejectedClonesPixels[i]>max1) max1=rejectedClonesPixels[i];
  for (i=0;i<smoothPixelsPSF.length;i++) if (smoothPixelsPSF[i]>max1) max1=smoothPixelsPSF[i];
  double minSigma=0.5;
  double varSigmaTop=1.0 ; //0.7;
//  double t1=max1*varSigmaTop;
  double kk;

  for (i=0;i<sigmas.length;i++) {
         kk=smoothPixelsPSF[i]/max1;
         if (kk>varSigmaTop) sigmas[i]=minSigma;  
         else                sigmas[i] = minSigma+ sigmas[i]*((varSigmaTop-kk)*(varSigmaTop-kk)/varSigmaTop/varSigmaTop);
  }
  double [] varFilteredPSF=variableGaussBlurr(rejectedClonesPixels, // input square pixel array, preferrably having many exact zeros (they will be skipped)
                                                            sigmas, // array of sigmas to be used for each pixel, matches pixels[]
                                                               3.5, // drop calculatin if farther then nSigma
                                                                 0, // int WOICenterX, // window of interest in pixels[] array - do not generate data outside it
                                                                 0, // int WOICenterY, // 
                                                           outSize, //int WOIWidth, reduce later
                                                           outSize); //int WOIHeight)


    if (DEBUG_LEVEL>2) {
/* Sigmas are 0 here ??? */
      if (sigmaToRadius>0.0) {
        float [] floatPixelsSigmas=new float[sigmas.length];
        for (j=0;j<sigmas.length;j++) floatPixelsSigmas[j]=(float) sigmas[j];
        ImageProcessor ip_Sigmas=new FloatProcessor(outSize,outSize);
        ip_Sigmas.setPixels(floatPixelsSigmas);
        ip_Sigmas.resetMinAndMax();
        ImagePlus imp_Sigmas=  new ImagePlus(title+"_Sigmas", ip_Sigmas);
        imp_Sigmas.show();
      }
      if (enableModelSubtract) {

        float [] floatPixelsModel=new float[modelPixels.length];
        for (j=0;j<modelPixels.length;j++) floatPixelsModel[j]=(float) modelPixels[j];
        ImageProcessor ip_Model=new FloatProcessor(outSize,outSize);
        ip_Model.setPixels(floatPixelsModel);
        ip_Model.resetMinAndMax();
        ImagePlus imp_Model=  new ImagePlus(title+"_Model", ip_Model);
        imp_Model.show();

        float [] floatPixelsDiff=new float[modelPixels.length];
        for (j=0;j<modelPixels.length;j++) floatPixelsDiff[j]=(float) rejectedClonesPixels[j];
        ImageProcessor ip_Diff=new FloatProcessor(outSize,outSize);
        ip_Diff.setPixels(floatPixelsDiff);
        ip_Diff.resetMinAndMax();
        ImagePlus imp_Diff=  new ImagePlus(title+"_Diff", ip_Diff);
        imp_Diff.show();

        float [] floatPixelsDiffSmooth=new float[modelPixels.length];
        for (j=0;j<modelPixels.length;j++) floatPixelsDiffSmooth[j]=(float) smoothRejectedClonesPixels[j];
        ImageProcessor ip_DiffSmooth=new FloatProcessor(outSize,outSize);
        ip_DiffSmooth.setPixels(floatPixelsDiffSmooth);
        ip_DiffSmooth.resetMinAndMax();
        ImagePlus imp_DiffSmooth=  new ImagePlus(title+"_DiffSmooth", ip_DiffSmooth);
        imp_DiffSmooth.show();
      }
/*
      float [] floatPixelsDiff1=new float[modelPixels.length];
      for (j=0;j<modelPixels.length;j++) floatPixelsDiff1[j]=(float) (pixelsPSF[j]-1.41421356*modelPixels[j]);
      ImageProcessor ip_Diff1=new FloatProcessor(outSize,outSize);
      ip_Diff1.setPixels(floatPixelsDiff1);
      ip_Diff1.resetMinAndMax();
      ImagePlus imp_Diff1=  new ImagePlus(title+"_Diff1", ip_Diff1);
      imp_Diff1.show();
*/


      System.out.println("title="+title+" center X(pix)="+centroidXY_smooth[0]+"(smooth) center Y(pix)="+centroidXY_smooth[1]+"(smooth)");
      System.out.println("title="+title+" center X(pix)="+centroidXY[0]+"          center Y(pix)="+centroidXY[1]);
/*
      float [] floatPSFClusterMask=new float[pixelsPSF.length];
//      for (j=0;j<pixelsPSF.length;j++) floatPSFClusterMask[j]=(float)clusterMask[j/outSize][j%outSize];
      int [] cm1d=convert2d_1d(clusterMask);
      for (j=0;j<pixelsPSF.length;j++) floatPSFClusterMask[j]=(float)cm1d[j];
      ImageProcessor ip_PSFClusterMask=new FloatProcessor(outSize,outSize);
      ip_PSFClusterMask.setPixels(floatPSFClusterMask);
      ip_PSFClusterMask.resetMinAndMax();
      ImagePlus imp_PSFClusterMask=  new ImagePlus(title+"_PSFcluster-"+topCenter, ip_PSFClusterMask);
      imp_PSFClusterMask.show();


      float [] floatSmoothPixelsPSF=new float[pixelsPSF.length];
      for (j=0;j<pixelsPSF.length;j++) floatSmoothPixelsPSF[j]=(float)smoothPixelsPSF[j];
      ImageProcessor ip_smoothPSFdbg=new FloatProcessor(outSize,outSize);
      ip_smoothPSFdbg.setPixels(floatSmoothPixelsPSF);
      ip_smoothPSFdbg.resetMinAndMax();
      ImagePlus imp_smoothPSFdbg=  new ImagePlus(title+"_PSFsmooth-"+smoothSeparate, ip_smoothPSFdbg);
      imp_smoothPSFdbg.show();
*/
    }
//                      double   sigmaToRadius,   // 0.4;  variable-sigma blurring to reduce high frequencies more for the pixels farther from the PSF center


//    return  pixelsPSF;
//    return  rejectedClonesPixels;
    centroid_xy[0]=centroidXY[0];
    centroid_xy[1]=centroidXY[1];
    return  varFilteredPSF;
  }

/* Trying variable-kernel filtering so the high frequencies are rejected more when farther from the center
    Can not use two 1-d passes, as radius will change. But it still can work with small enough sigmaToRadius values if performance will be too low with brute force 2-d filtering  */

  double [] variableGaussBlurr0 (double []pixels, // input square pixel array, preferrably having many exact zeros (they will be skipped)
                         double sigmaToRadius, // sigma is proportional to the distance from the center
//                         double radiusSigma, // sigma will grow proportionally to radius^2, at this radius sigma will be equal to radius
                                double nSigma, // drop calculatin if farther then nSigma
                                    double xc, // model PSF center X-coordinate (in pixels[] units, from the center of the array )
                                    double yc, // same for Y
                               int WOICenterX, // window of interest in pixels[] array - do not generate data outside it
                               int WOICenterY, // 
                                 int WOIWidth, //
                                int WOIHeight){ //
    int size = (int) Math.sqrt(pixels.length);
    double [] result =new double [size*size];
    double [] gauss= new double [2*size];
    int x0= (size-WOIWidth)/2 +WOICenterX;
    int y0= (size-WOIHeight)/2+WOICenterY;
    int x1=x0+WOIWidth;
    int y1=x0+WOIHeight;
    int i,ix,iy,max_i;
    double r, sum,k,x,y,sigma,d,dky,scale;
    int xk0,xk1,yk0,yk1, ikx,iky;
    for (i=0;i<result.length;i++) result[i]=0.0;

    if (x0<0) x0=0; if (x1>size) x1=size; if (y0<0) y0=0; if (y1>size) y1=size;
    for (iy=0;iy<size;iy++) {
      y=(iy-size/2)-yc;
      for (ix=0;ix<size;ix++) {
        d=pixels[iy*size+ix];
        if (d!=0.0) {
          x=(ix-size/2)-xc;
          r=Math.sqrt(x*x+y*y);
//          sigma=r*sigmaToRadius;
//          sigma=r*r/radiusSigma;
          sigma=(r*sigmaToRadius)+1;
          if (sigma==0.0) {
            result[iy*size+ix]+=d; // just copy input data, no convolving
          } else {
            max_i= (int) (sigma*nSigma+1);
            k=1.0/(2.0*sigma*sigma);
            if (max_i>=gauss.length) max_i=gauss.length-1;
            sum=-0.5; // 0 is counted twice
            for (i=0; i<=max_i; i++) {
              gauss[i]=Math.exp(-k*i*i);
              sum+= gauss[i]; // could use - more errors for small values of gamma 1/Math.sqrt(2*Math.PI*sigma*sigma)
            }
            scale=0.5/sum;
            for (i=0; i<=max_i; i++) gauss[i]*=scale;
            yk0=-max_i; if (yk0<(y0-iy)) yk0=y0-iy;
            yk1= max_i; if (yk1>=(y1-iy)) yk1=y0-iy-1;
            xk0=-max_i; if (xk0<(x0-ix)) xk0=x0-ix;
            xk1= max_i; if (xk1>=(x1-ix)) xk1=x0-ix-1;
/*if ((DEBUG_LEVEL>1)&& (r<2)) {
      System.out.println(" iy="+iy+" ix="+ix+" x="+x+" y="+y+" r="+r+" sum="+sum+" scale="+scale+" k="+k+" max_i="+max_i+" sigma="+sigma);
      System.out.println(" yk0="+yk0+" yk1="+yk1+" xk0="+xk0+" xk1="+xk1);
      for (ii=0;ii<=max_i;ii++) System.out.println(" gauss["+ii+"]="+gauss[ii]);
      sum=0;
      for (iky=yk0;iky<=yk1;iky++) for (ikx=xk0;ikx<=xk1;ikx++) sum+=gauss[Math.abs(iky)]*gauss[Math.abs(ikx)];
      System.out.println(" sum of gauss="+sum);

}*/
            for (iky=yk0;iky<=yk1;iky++) {
              dky=d*gauss[Math.abs(iky)];
              for (ikx=xk0;ikx<=xk1;ikx++) {
                result[(iy+iky)*size+ix+ikx]+=dky*gauss[Math.abs(ikx)];
/*if ((DEBUG_LEVEL>1)&& (r<2)) {
      System.out.println("     iky="+iky+" ikx="+ikx+" dky="+dky+" iy+iky="+(iy+iky)+" ix+ikx="+(ix+ikx)+" dky*gauss[Math.abs(ikx)]="+(dky*gauss[Math.abs(ikx)])+" result[(iy+iky)*size+ix+ikx]="+result[(iy+iky)*size+ix+ikx]);
}*/

              }
            }
          }
        }
      }
    }
    return result;
  }

/* create aray (to be used with variableGaussBlurr() ) of per-pixel sigma values for gauss blur, proportional to distance from the specified center */
  double [] createSigmasRadius (double []pixels, // input square pixel array, preferrably having many exact zeros (they will be skipped)
                           double sigmaToRadius, // sigma is proportional to the distance from the center
                                      double xc, // model PSF center X-coordinate (in pixels[] units, from the center of the array )
                                      double yc, // same for Y
                                 int WOICenterX, // window of interest in pixels[] array - do not generate data outside it
                                 int WOICenterY, // 
                                   int WOIWidth, //
                                  int WOIHeight) {
    int size = (int) Math.sqrt(pixels.length);
    double [] sigmas =new double [size*size];
    int x0= (size-WOIWidth)/2 +WOICenterX;
    int y0= (size-WOIHeight)/2+WOICenterY;
    int x1=x0+WOIWidth;
    int y1=x0+WOIHeight;
    int i,ix,iy;
    double r,x,y;
    for (i=0;i<sigmas.length;i++) sigmas[i]=0.0;
    if (x0<0) x0=0; if (x1>size) x1=size; if (y0<0) y0=0; if (y1>size) y1=size;
    for (iy=0;iy<size;iy++) {
      y=(iy-size/2)-yc;
      for (ix=0;ix<size;ix++) {
        x=(ix-size/2)-xc;
        r=Math.sqrt(x*x+y*y);
//          sigma=r*sigmaToRadius;
//          sigma=r*r/radiusSigma;
//          sigmas[iy*size+ix]=(r*sigmaToRadius)+1;
        sigmas[iy*size+ix]=(r*sigmaToRadius);
      }
    }


    return sigmas;
  }


  double [] variableGaussBlurr (double []pixels, // input square pixel array, preferrably having many exact zeros (they will be skipped)
                                double []sigmas, // array of sigmas to be used for each pixel, matches pixels[]
                                double nSigma, // drop calculatin if farther then nSigma
                               int WOICenterX, // window of interest in pixels[] array - do not generate data outside it
                               int WOICenterY, // 
                                 int WOIWidth, //
                                int WOIHeight){ //
    int size = (int) Math.sqrt(pixels.length);
    double [] result =new double [size*size];
    double [] gauss= new double [2*size];
    int x0= (size-WOIWidth)/2 +WOICenterX;
    int y0= (size-WOIHeight)/2+WOICenterY;
    int x1=x0+WOIWidth;
    int y1=x0+WOIHeight;
    int i,ix,iy,max_i;
    double sum,k,sigma,d,gy,scale,g;
    int xk0,xk1,yk0,yk1, ikx,iky, index;
    for (i=0;i<result.length;i++) result[i]=0.0;
//debug  
/*
    double [] dbg_result1 =new double [size*size];
    double [] dbg_result2 =new double [size*size];
    for (i=0;i<result.length;i++) dbg_result1[i]=0.0;
    for (i=0;i<result.length;i++) dbg_result2[i]=0.0;
*/
if (DEBUG_LEVEL>2) {
      System.out.println(" variableGaussBlurr(), x0="+x0+" y0="+y0+" x1="+x1+" y1="+y1);
}
    if (x0<0) x0=0; if (x1>size) x1=size; if (y0<0) y0=0; if (y1>size) y1=size;
    for (iy=0;iy<size;iy++) {
      for (ix=0;ix<size;ix++) {
        d=pixels[iy*size+ix];
        if (d!=0.0) {
          sigma=sigmas[iy*size+ix];
          if (sigma==0.0) {
            result[iy*size+ix]+=d; // just copy input data, no convolving
          } else {
/* opposite to "normal" convolution we have diffrent kernel for each point, so we need to make sure that two points with the same values but
    diffrent sigma values will not move "energy" from one to another. For this we can do accumulation both ways - from the source point to all
     points "reachable" by the kernel (proportional to the pixel value) and also in opposite direction - from those other points to the current
     pointer (where kernel is centered) with the value proportional to that othre point  */

            max_i= (int) (sigma*nSigma+1);
            k=1.0/(2.0*sigma*sigma);
            if (max_i>=gauss.length) max_i=gauss.length-1;
            sum=-0.5; // 0 is counted twice
            for (i=0; i<=max_i; i++) {
              gauss[i]=Math.exp(-k*i*i);
              sum+= gauss[i]; // could use - more errors for small values of gamma 1/Math.sqrt(2*Math.PI*sigma*sigma)
            }
            scale=0.5/sum;
            for (i=0; i<=max_i; i++) gauss[i]*=scale;
            yk0=-max_i; if (yk0<(y0-iy)) yk0=y0-iy;
            yk1= max_i; if (yk1>=(y1-iy)) yk1=y1-iy-1;
            xk0=-max_i; if (xk0<(x0-ix)) xk0=x0-ix;
            xk1= max_i; if (xk1>=(x1-ix)) xk1=x1-ix-1;

//if ((DEBUG_LEVEL>1)&& (Math.abs(iy-128)<3)&& (Math.abs(ix-128)<3)) {
/*if ((DEBUG_LEVEL>1)&& (Math.abs(iy-16)<3)&& (Math.abs(ix-16)<3)) {
      System.out.println(" iy="+iy+" ix="+ix+" sum="+sum+" scale="+scale+" k="+k+" max_i="+max_i+" sigma="+sigma);
      System.out.println(" yk0="+yk0+" yk1="+yk1+" xk0="+xk0+" xk1="+xk1);
      for (ii=0;ii<=max_i;ii++) System.out.println(" gauss["+ii+"]="+gauss[ii]);
      sum=0;
      for (iky=yk0;iky<=yk1;iky++) for (ikx=xk0;ikx<=xk1;ikx++) sum+=gauss[Math.abs(iky)]*gauss[Math.abs(ikx)];
      System.out.println(" sum of gauss="+sum);
}*/

            for (iky=yk0;iky<=yk1;iky++) {
              gy=gauss[Math.abs(iky)]/2; // Extra /2 because we'll calculate the convolution twice from the [ix,iy] and to [ix,iy]
              for (ikx=xk0;ikx<=xk1;ikx++) {
                index=(iy+iky)*size+ix+ikx;
                g=gy*gauss[Math.abs(ikx)];
                result[index]+=d*g;
                result[iy*size+ix]+=pixels[index]*g;
/*
                dbg_result1[index]+=d*g;
                dbg_result2[iy*size+ix]+=pixels[index]*g;
*/
/*if ((DEBUG_LEVEL>1)&& (r<2)) {
      System.out.println("     iky="+iky+" ikx="+ikx+" dky="+dky+" iy+iky="+(iy+iky)+" ix+ikx="+(ix+ikx)+" dky*gauss[Math.abs(ikx)]="+(dky*gauss[Math.abs(ikx)])+" result[(iy+iky)*size+ix+ikx]="+result[(iy+iky)*size+ix+ikx]);
}*/

              }
            }
          }
        }
      }
    }
/*
    SDFA_instance.showArrays(dbg_result1, size,size, "dbg_result1");
    SDFA_instance.showArrays(dbg_result2, size,size, "dbg_result2");
    SDFA_instance.showArrays(result,      size,size, "result");
*/
    return result;
  }
















/* generate a single instance of the model PSF (4 diffrent Gauss widths along 4 directions) based on 2 intre-center vectors */
  double [] addModelPSF(double [] pixels, // square pixel array where the model PSF is added
                            double scale, // scale model PSF (it has 1.0 in the center) before additing to pixels [] array
                               double xc, // model PSF center X-coordinate (in pixels[] units, from the center of the array )
                               double yc, // same for Y
                          int WOICenterX, // window of interest in pixels[] array - do not generate data outside it
                          int WOICenterY, // 
                            int WOIWidth, //
                           int WOIHeight, //
                      double[][] vectors, // vectors that connect center of PSF with two oppositre sign clones
                       double [] widths0, // Gauss widths along vectors[0][] in p[ositive [0] and negative [1] directions
                       double [] widths1 // Gauss widths along vectors[1][] in p[ositive [0] and negative [1] directions
                                        ) {
    int ix,iy;
    int size = (int) Math.sqrt (pixels.length);
    double [] xy= new double[2];
    double [] uv;
    double d;
/* matrix that converts u,v (lengths along the) 2 input vectors connecting opposite sign PSFs into x,y coordinates */
    double [] vectorsLengths={Math.sqrt(vectors[0][0]*vectors[0][0]+vectors[0][1]*vectors[0][1]),
                              Math.sqrt(vectors[1][0]*vectors[1][0]+vectors[1][1]*vectors[1][1])};
    double [][] uv2xy= {{vectors[0][0]/vectorsLengths[0],vectors[1][0]/vectorsLengths[1]},
                        {vectors[0][1]/vectorsLengths[0],vectors[1][1]/vectorsLengths[1]}};
    double [][] xy2uv=  matrix2x2_invert(uv2xy);
    int x0= (size-WOIWidth)/2 +WOICenterX;
    int y0= (size-WOIHeight)/2+WOICenterY;
    int x1=x0+WOIWidth;
    int y1=x0+WOIHeight;
    if (x0<0) x0=0; if (x1>size) x1=size; if (y0<0) y0=0; if (y1>size) y1=size;
    double [][]widths2={{widths0[0]*widths0[0],widths0[1]*widths0[1]},
                        {widths1[0]*widths1[0],widths1[1]*widths1[1]}};
/* See if all widths are non-zero - some could be after approximation if the PSF was too far/narrow for that method used */
    if ((widths0[0]==0.0) || (widths0[1]==0.0) || (widths1[0]==0.0) || (widths1[1]==0.0)) {
      return pixels; // return with no changes
    }
    for (iy=y0;iy<y1;iy++) {
      xy[1]=(iy-size/2)-yc;
      for (ix=x0;ix<x1;ix++) {
        xy[0]=(ix-size/2)-xc;
        uv=matrix2x2_mul(xy2uv, xy);
        d=scale;
        if (uv[0]>0) d*=Math.exp(-((uv[0]*uv[0]/widths2[0][0])));
        else         d*=Math.exp(-((uv[0]*uv[0]/widths2[0][1])));
        if (uv[1]>0) d*=Math.exp(-((uv[1]*uv[1]/widths2[1][0])));
        else         d*=Math.exp(-((uv[1]*uv[1]/widths2[1][1])));
        pixels[iy*size+ix]+=d;
      }
    }
    return pixels;
  }

/* zeroes out area outside of the area bound by 4 negative clones (or a fraction of it), either sharp or with Hamming
 */
  double [] maskClonesPSF(double [] pixels, // square pixel array where the model PSF is added
                         double windowPart, // multiply window by this value
                                 double xc, // Center of the remaining single PSF
                                 double yc, // same for Y
                        double[][] vectors, // vectors that connect center of PSF with two oppositre sign clones
                       boolean  useHamming  // use Hamming window, if false - just cut sharp
                                        ) {
    int ix,iy;
    int size = (int) Math.sqrt (pixels.length);
    double [] xy= new double[2];
    double [] uv;
    /* matrix that converts u,v (lengths along the) 2 input vectors connecting opposite sign PSFs into x,y coordinates */
    double [][] uv2xy= {{vectors[0][0]*windowPart,vectors[1][0]*windowPart},
                        {vectors[0][1]*windowPart,vectors[1][1]*windowPart}};
    double [][] xy2uv=  matrix2x2_invert(uv2xy);
    for (iy=0;iy<size;iy++) {
      xy[1]=(iy-size/2)-yc;
      for (ix=0;ix<size;ix++) {
        xy[0]=(ix-size/2)-xc;
        uv=matrix2x2_mul(xy2uv, xy);
        if ((Math.abs(uv[0])>1.0) || (Math.abs(uv[1])>1.0)) pixels[iy*size+ix]=0.0;
        else if (useHamming) {
          pixels[iy*size+ix]*=(0.54+0.46*Math.cos(uv[0]*Math.PI))*(0.54+0.46*Math.cos(uv[1]*Math.PI));
        }
      }
    }
    return pixels;
  }





/*  builds profile between positive PSF maximum (center[]) and negative (along vector[]) and finds
     approximation by two pieces of gauss functions, measuring zero point and derivative in that point.
     threshold is used to remove center uncertainty range (if any): two crossings of positive and
     negative thresholds (nearest to the PSF opposite maximums) are connected straight - that defines
     zero crossing and the derivative  */
/*
  double [] approximateOppositeGauss (double [] pixels, // square pixel array
                                      int         size, // length of connecting section line, pixels
                                      double threshold, // uncertainty range - fraction of max
                                      double x0, // x of the center maximum (positive), relative to the center of array
                                      double y0,
                                      double vx, // vector from the positive center to the negative one
                                      double vy) {
    return null;
  }
*/
  double [] approximateOppositeGauss (double [] section, // Array of pixes along a section, connecting positive and negative maximums
                                      double sectionLength, // length of section in pixelsPSF
                                      double initialGauss   // width of the "smoothing" filter applied to the data before sestion (in pixels)
                                         ) {
    double [] resultGauss={0.0,0.0}; // first is from the side of the positive PSF, second - from negative
    int size=section.length-1;
    int i,n;
/* First we need to find the maximum on the derivative (it is negative, so maximum is where derivative has the smallest absolute value)
    That would help to find (and compensate for ) the assymmetry of the PSF in this section. Using zero crossing is more prone to errors,
    artifact low-frequency "interference" pattern with the preriod of the grid shifts the values up/down, and with relatively lower
    derivative value it moves zero crossing point significantly. So watching extrenums on derivative seems to work better (may be some other errors)*/
    double normalizeScale=2.0/(section[0]-section[size]);
    double [] sectionDerivative=new double[size+1];
    for (i=1; i<size; i++) sectionDerivative[i]=0.5*(section[i+1]-section[i-1])*normalizeScale;
    sectionDerivative[0]= sectionDerivative[1];
    sectionDerivative[size]= sectionDerivative[size-1];
    double [] minsOnDerivative    = new double[2];
    int    [] minsOnDerivativePos = new int[2];
    double maxOnDerivative;
    int    maxOnDerivativePos;
    for (n=1;n<2;n++) {
      minsOnDerivative[n]=0.0;
      minsOnDerivativePos[n]=0;
    }
    for (i=1;i<size;i++) {
       n=(i>=size/2)?1:0;
       if (sectionDerivative[i]<minsOnDerivative[n]) {
         minsOnDerivative[n]=sectionDerivative[i];
         minsOnDerivativePos[n]=i;
       }
    }
    if ( (minsOnDerivativePos[0]==0) ||  (minsOnDerivativePos[1]==0)) {
      System.out.println("ERROR - failed to find two minimums on digital the derivative of the section between PSFs");
      if (DEBUG_LEVEL>1) for (i=1; i<size; i++)System.out.println("sectionDerivative["+i+"]="+sectionDerivative[i]);
      return resultGauss; // {0.0,0.0};
    }
    maxOnDerivative= minsOnDerivative[0];
    maxOnDerivativePos=0;
    for (i=minsOnDerivativePos[0]+1;i<minsOnDerivativePos[1];i++) {
       if (sectionDerivative[i]>maxOnDerivative) {
         maxOnDerivative=sectionDerivative[i];
         maxOnDerivativePos=i;
       }
    }

    if (DEBUG_LEVEL>2) {
       System.out.println(" minsOnDerivative[0]="+minsOnDerivative[0]+" minsOnDerivative[1]="+minsOnDerivative[1]+
                          " minsOnDerivativePos[0]="+minsOnDerivativePos[0]+" minsOnDerivativePos[1]="+minsOnDerivativePos[1]+
                          " maxOnDerivative="+maxOnDerivative+" maxOnDerivativePos="+maxOnDerivativePos);
//      if (DEBUG_LEVEL>1) for (i=1; i<size; i++)System.out.println("sectionDerivative["+i+"]="+sectionDerivative[i]);
    }


    if ( maxOnDerivativePos==0) {
      System.out.println("ERROR - failed to find maximum (lowest negative) on the digital derivative of the section between PSFs");
/* there is a stricter limit on the value, not just 0 - we know input is convolved with the Gauss*/
      return resultGauss; // {0.0,0.0};
    }
/* Now validate that the derivative value makes sense - if it is positive - need to smooth more, nothing to do here?
    Or just filter it right away?
    It is possible - just connect the two pieces of section (negative fro the second one), they should have the matching ends. Then - 1d-convolve
    with Gauss, adding to the initial one as sqrt(w1^2 + w2^2), It can probably make sense to use the same width (steps) as originally was applied.
    That we may want to improve later.
    */
    if ( maxOnDerivative>=0.0) {
      System.out.println("ERROR - maximum (lowest negative) on the digital derivative is positive - increase initial smoothing or read comments in the code above.");
      return resultGauss; // {0.0,0.0};
    }
/* We can improve location of this maximum by second degree polinominal approximation in the vicinity of this point - later */
    double derivMax=1.0*maxOnDerivativePos;
/* At this point both positive and negative should have the same derivative, equal to 1/2 of what we see - negative maxOnDerivative
   exp(-((x/width_positive)^2))/width_positive^2 should be equal to maxOnDerivative * k^2/(2*((k^2-1)*derivMax+size) when width is
   measured in units of the section.
   Don't know of any analytical solution, do it numerically*/

/* convert everything to 0..1 interval */
    double x0=derivMax/size;
//    double k= (size-derivMax)/derivMax; // ratio of the Gauss width from the negative to the gauss width from the positive
    double k= (1-x0)/x0; // ratio of the Gauss width from the negative to the gauss width from the positive
    double w0=initialGauss/sectionLength; // Gauss width of the smoothing filter already applied.
    double derMax=-maxOnDerivative*size; // make it positive
    double targetVal=derMax/(4*x0);
    double wMax=1.0/x0/x0/Math.E;
    double wMin=Math.exp (-x0*x0/w0/w0)/w0/w0;
    double w,v,wl,wh,vl,vh,step;
    double minStep;

    if (DEBUG_LEVEL>2) {
       System.out.println("Starting Gauss approximation, x0="+x0+" w0="+w0+" derMax="+derMax+" targetVal="+targetVal);
    }
    if (targetVal>wMax) {
       System.out.println("Derivative is too high to use, requested "+targetVal+", possible <="+wMax);
       w=x0;
    } else if (targetVal<wMin){
       System.out.println("Derivative is too low to use, requested "+targetVal+", possible >="+wMin);
       w=w0;
    } else {
//      w=Math.sqrt(x0*w0); //start of approximation - reasonable solution should be between these points (x0 and w0)
/* exp (-(x/w)^2)/w^2  is a monotonic function of w on this interval, so we can just divide interval - no need to increase performance by using derivative */
      wl=w0;
      wh=x0;
      w=wl;
      minStep=0.0001*(wh-wl); // final precision determined by this
      vl=Math.exp (-x0*x0/wl/wl)/wl/wl;
      vh=Math.exp (-x0*x0/wh/wh)/wh/wh;
      step=vh-vl;
      while (step>minStep) {
        if (DEBUG_LEVEL>2) {
          System.out.println("Approximating: wl="+wl+" wh="+wh+" vl="+vl+" vh="+vh+" targetVal="+targetVal+" minStep="+minStep+" last step="+step);
        }
        w=wl+(wh-wl)*(targetVal-vl)/(vh-vl);
        v=Math.exp (-x0*x0/w/w)/w/w;
        if (v>targetVal) {
          step=wh-w;
          wh=w;
          vh=v;
        } else {
          step=w-wl;
          wl=w;
          vl=v;
        }
      }
    }
/* removing "smoothing" Gauss from the result:*/
    w=Math.sqrt(w*w-w0*w0);
    resultGauss[0]=  w*sectionLength;
    resultGauss[1]=k*w*sectionLength;
    if (DEBUG_LEVEL>1) {
      System.out.println("Corrected for initial smoothing: w="+w+" resultGauss[0]="+resultGauss[0]+"pix, resultGauss[1]="+resultGauss[1]+"pix");
    }


    return resultGauss;
  }



/*
 
public static void showTrace(String msg)
{
  if (msg.length() > 0) System.out.println(msg);
  System.out.println("Trace: " + 
           "file " + new Throwable().getStackTrace()[1].getFileName() +
           " class " + new Throwable().getStackTrace()[1].getClassName() +
           " method " + new Throwable().getStackTrace()[1].getMethodName() +
           " line " + new Throwable().getStackTrace()[1].getLineNumber());
}

*/
  double [] createSection (    double [] pixels, // square pixel array
                               int         size, // length of connecting section line, pixels
                                      double x0, // x of the center maximum (positive), relative to the center of array
                                      double y0,
                                      double vx, // vector from the positive center to the negative one
                                      double vy) {

    double [] section=new double [size+1];
    int pixelSize=(int)Math.sqrt(pixels.length);
    int i,ix,iy;
    double x,y,dx,dy;
    double sx=vx/size;
    double sy=vy/size;
    double xc=x0+pixelSize/2;
    double yc=y0+pixelSize/2;
    double [][] corners=new double [2][2];
    double []   sides=new   double [2];
    
    for (i=0;i<=size;i++) {
      x=xc+i*sx;
      y=yc+i*sy;
      if (x<0)x=0.0;
      else if (x>(pixelSize-2.0)) x=pixelSize-2.0;
      if (y<0)y=0.0;
      else if (y>(pixelSize-2.0)) y=pixelSize-2.0;
      ix=(int) x;
      dx=x-ix;
      iy=(int) y;
      dy=y-iy;
/* bilinear interpolation */
      corners[0][0]=pixels[pixelSize*(iy)   + (ix  )];
      corners[0][1]=pixels[pixelSize*(iy)   + (ix+1)];
      corners[1][0]=pixels[pixelSize*(iy+1) + (ix  )];
      corners[1][1]=pixels[pixelSize*(iy+1) + (ix+1)];
      sides[0]=corners[0][0]*(1-dx)+corners[0][1]*dx;
      sides[1]=corners[1][0]*(1-dx)+corners[1][1]*dx;
      section[i]=sides[0]*(1-dy)+sides[1]*dy;
    }
    return section;

  }



  private double PSFAtXY(double [] pixels, int size, double x, double y) {
      int ix=(int) Math.round(x);
      int iy=(int) Math.round(y);
      if      (ix <  -size/2) ix=-size/2;
      else if (ix >=  size/2) ix= size/2-1;
      if      (iy <  -size/2) iy=-size/2;
      else if (iy >=  size/2) iy= size/2-1;
//      if ((ix<= -size/2) || (ix>= size/2) || (iy<= -size/2) || (iy>= size/2)) return 0.0;
      int index=size* (size/2 + iy)+ size/2 + ix;
      if ((index<0) || (index > pixels.length)) {
System.out.println("PSFAtXY error, x="+IJ.d2s(x,0)+" y="+IJ.d2s(y,0)+ " index="+(size*(size/2 + (int) Math.round(y))+ size/2 + (int) Math.round(x))+ " pixels.length="+pixels.length);
      }
//      return pixels[size* (size/2 + (int) Math.round(y))+ size/2 + (int) Math.round(x)]; // out of bounds
      return pixels[index];
  }
  private double contrastAtXY(int sign, double [] pixels, int size, double x, double y, double [][] g, double [] cache) {
//      int ir= (int) Math.round(0.125*Math.min(Math.max(g[0][0],g[1][0]),Math.max(g[0][1],g[1][1]))); // sample at square 1 1/2x1/2 of the grid "square"
//      int ir= (int) Math.round(0.25*Math.min(Math.max(g[0][0],g[1][0]),Math.max(g[0][1],g[1][1]))); // sample at square 1 1/2x1/2 of the grid "square"
      int ir= (int) Math.round(0.2*Math.min(Math.max(Math.abs(g[0][0]),Math.abs(g[1][0])),Math.max(Math.abs(g[0][1]),Math.abs(g[1][1])))); // sample at square 1 1/2x1/2 of the grid "square"

      int ix=(int) Math.round(x);
      int iy=(int) Math.round(y);
      if      (ix <  -size/2) ix=-size/2;
      else if (ix >=  size/2) ix= size/2-1;
      if      (iy <  -size/2) iy=-size/2;
      else if (iy >=  size/2) iy= size/2-1;
      int index= size* (size/2 + iy)+ size/2 + ix;
//      if ((cache!=null) && (cache[index]>=0)) return sign*cache[index];
      if ((cache!=null) && (cache[index]>=0)) return cache[index];
      double rslt=0.0;
      int i,j;
      for (i=-ir;i<=ir;i++) for (j=-ir;j<=ir;j++) {
        rslt+=     PSFAtXY(pixels,size,j+ix,i+iy) -
            0.25* (PSFAtXY(pixels,size,j+ix+(g[0][0]+ g[1][0])/2  ,i+iy+(g[0][1]+ g[1][1])/2)+
                   PSFAtXY(pixels,size,j+ix+(g[0][0]- g[1][0])/2  ,i+iy+(g[0][1]- g[1][1])/2)+
                   PSFAtXY(pixels,size,j+ix-(g[0][0]+ g[1][0])/2  ,i+iy-(g[0][1]+ g[1][1])/2)+
                   PSFAtXY(pixels,size,j+ix-(g[0][0]- g[1][0])/2  ,i+iy-(g[0][1]- g[1][1])/2));

      }
      rslt=rslt*sign;
      cache[index] = (rslt>0.0)?rslt:0.0;
 //     System.out.println("contrastAtXY("+(ix+size/2)+","+(iy+size/2)+ ")="+rslt);
      return rslt/ir/ir;
  }

  public double[] cropAndHammiongKenel (double[]full_kernel, int halfsize) {
     int full_size=(int) Math.sqrt(full_kernel.length);
     int size=(2*halfsize+1);
     int i,y,x,index,full_index;
     double [] kernel = new double[size*size];
     double [] hamming_line=new double [size];
     double a,k;
     if (GAUSS_WIDTH<=0) {
  //     for (i=0; i<size; i++) hamming_line[i]= (0.54-0.46*Math.cos((i*2.0*Math.PI)/size));
       for (i=0; i<size; i++) hamming_line[i]= (0.54-0.46*Math.cos((i*2.0*Math.PI)/(size-1)));
     } else {
       k=2.0/(size*GAUSS_WIDTH);
       for (i=0; i<size; i++) {
         a=(i-size/2)*k;
         hamming_line[i]= Math.exp( - a*a);
       }
     }
     index=0;
     for (y=0; y<size; y++) {
       full_index= full_size* (full_size/2-halfsize +y) + (full_size/2 -halfsize);
       for (x=0; x<size; x++) kernel[index++]=hamming_line[y]*hamming_line[x]*full_kernel[full_index++];
     }
     return kernel;
  }

  public int kernelLength(double[][]kernels) {
    if (kernels==null) return 0;
    for (int i=0; i<kernels.length;i++) if (kernels[i]!=null) return kernels[i].length;
    return 0;
  }

/* extends/shrinks image to make it square for FFT */
  public double[][] resizeForFFT (double[][]kernels, int size) {
    if (kernels==null) return null;
    double [][]result=new double [kernels.length][];
    for (int i=0;i<kernels.length;i++) {
      if (kernels[i]!=null) result[i]=resizeForFFT(kernels[i],size);
      else result[i]=null;
    }
    return result;
  }


  public double[] resizeForFFT (double[]kernel, int size) {
     int ksize=(int) Math.sqrt(kernel.length);
     double [] kernelForFFT = new double[size*size];
     int i,j,index, full_index;
     if (DEBUG_LEVEL>10) System.out.println("resizeForFFT: new size="+size+" old size="+ksize);
     index=0;
     if (size==ksize) {
        return kernel.clone();
     } else if (size>ksize) {
       for (full_index=0;full_index<kernelForFFT.length; full_index++) kernelForFFT [full_index]=0.0;
       for (i=0;i<ksize; i++) {
         full_index=size* (size/2- ksize/2 + i) +size/2-ksize/2;
         for (j=0;j<ksize; j++) kernelForFFT[full_index++]=kernel[index++];
       }
     } else {
       for (i=0; i<size; i++) {
         full_index= ksize* (ksize/2-(size/2) +i) + (ksize/2-(size/2));
         for (j=0; j<size; j++) kernelForFFT[index++]=kernel[full_index++];
       }
     }
     return kernelForFFT;
  }



  public double[] initHamming(int size) {
    double [] hamming =new double [size*size];
    double [] hamming_line=new double [size];
    double a,k;
    int i,j;
    if (GAUSS_WIDTH<=0) {
      for (i=0; i<size; i++) hamming_line[i]= (0.54-0.46*Math.cos((i*2.0*Math.PI)/size));
    } else {
      k=2.0/(size*GAUSS_WIDTH);
      for (i=0; i<size; i++) {
         a=(i-size/2)*k;
         hamming_line[i]= Math.exp( - a*a);
      }
    }
    for (i=0; i<size; i++) for (j=0; j<size; j++){
       hamming[size*i+j]=hamming_line[i]*hamming_line[j];
    }
    return hamming;
  }

//GAUSS_WIDTH

/* converts FHT results (frequency space) to complex numbers of [fftsize/2+1][fftsize] */
public double[][][] FHT2FFTHalf (FHT fht, int fftsize) {
   float[] fht_pixels=(float[])fht.getPixels();
   double[][][] fftHalf=new double[(fftsize>>1)+1][fftsize][2];
   int row1,row2,col1,col2;

   for (row1=0;row1<=(fftsize>>1);row1++) {
     row2=(fftsize-row1) %fftsize;
     for (col1=0;col1<fftsize;col1++) {
       col2=(fftsize-col1) %fftsize;
//       fftHalf[row1][col1]=   complex( 0.5*(fht_pixels[row1*fftsize+col1] + fht_pixels[row2*fftsize+col2]),
//                                       0.5*(fht_pixels[row2*fftsize+col2] - fht_pixels[row1*fftsize+col1]));
       fftHalf[row1][col1][0]=   0.5*(fht_pixels[row1*fftsize+col1] + fht_pixels[row2*fftsize+col2]);
       fftHalf[row1][col1][1]=   0.5*(fht_pixels[row2*fftsize+col2] - fht_pixels[row1*fftsize+col1]);
     }
   }
   return fftHalf;
}

public double[][][] FHT2FFTHalf (double [] fht_pixels, int fftsize) {
   double[][][] fftHalf=new double[(fftsize>>1)+1][fftsize][2];
   int row1,row2,col1,col2;

   for (row1=0;row1<=(fftsize>>1);row1++) {
     row2=(fftsize-row1) %fftsize;
     for (col1=0;col1<fftsize;col1++) {
       col2=(fftsize-col1) %fftsize;
       fftHalf[row1][col1][0]=   0.5*(fht_pixels[row1*fftsize+col1] + fht_pixels[row2*fftsize+col2]);
       fftHalf[row1][col1][1]=   0.5*(fht_pixels[row2*fftsize+col2] - fht_pixels[row1*fftsize+col1]);
     }
   }
   return fftHalf;
}




/* converts FFT arrays of complex numbers of [fftsize/2+1][fftsize] to FHT arrays */
public float[] floatFFTHalf2FHT (double [][][] fft, int fftsize) {
   float[] fht_pixels=new float [fftsize*fftsize];
   int row1,row2,col1,col2;
   for (row1=0;row1<=(fftsize>>1);row1++) {
     row2=(fftsize-row1) %fftsize;
     for (col1=0;col1 < fftsize;col1++) {
       col2=(fftsize-col1) %fftsize;
       fht_pixels[row1*fftsize+col1]=(float)(fft[row1][col1][0]-fft[row1][col1][1]);
       fht_pixels[row2*fftsize+col2]=(float)(fft[row1][col1][0]+fft[row1][col1][1]);
     }
   }
   return fht_pixels;
}

public double[] FFTHalf2FHT (double [][][] fft, int fftsize) {
   double[] fht_pixels=new double [fftsize*fftsize];
   int row1,row2,col1,col2;
   for (row1=0;row1<=(fftsize>>1);row1++) {
     row2=(fftsize-row1) %fftsize;
     for (col1=0;col1 < fftsize;col1++) {
       col2=(fftsize-col1) %fftsize;
       fht_pixels[row1*fftsize+col1]=(double) (fft[row1][col1][0]-fft[row1][col1][1]);
       fht_pixels[row2*fftsize+col2]=(double) (fft[row1][col1][0]+fft[row1][col1][1]);
     }
   }
   return fht_pixels;
}

/* Interpolate between the two space-domain pixel arrays (i.e. PSFs) by interpolating amplitudes and phases for each frequency component*/
public double [] interpolatePSF(double [] pixels0,  // first square pixel array
                                double [] pixels1,  // second square pixel array
                                         double k,  // interpolation coefficient (0.0 - use ampPhase0, 1.0 - ampPhase1)
                                 boolean centered){ // true if result should have 0:0 in the center
  double [][][] ampPhase0=spaceToAmplPhase(pixels0, centered );
  double [][][] ampPhase1=spaceToAmplPhase(pixels1, centered );
  double [][][] ampPhase= interpolateAmplPhase(ampPhase0, // first set of amplitudes/phases
                                               ampPhase1, // second set of amplitudes/phases
                                                       k); // interpolation coefficient (0.0 - use ampPhase0, 1.0 - ampPhase1)
  return ampPhaseToSpace(ampPhase, centered);
}
/* Interpolate amplitudes/Phases*/
public double [][][] interpolateAmplPhase(double [][][] ampPhase0, // first set of amplitudes/phases
                                          double [][][] ampPhase1, // second set of amplitudes/phases
                                                         double k) { // interpolation coefficient (0.0 - use ampPhase0, 1.0 - ampPhase1)
  int i,j;
  if      (k==0.0) return ampPhase0.clone();
  else if (k==1.0) return ampPhase1.clone();
  double [][][] result=new double[ampPhase0.length][ampPhase0[0].length][2];
  for (i=0;i<ampPhase0.length;i++) for (j=0;j<ampPhase0[0].length;j++) {
     result[i][j][0]=Math.pow(ampPhase0[i][j][0],(1.0-k))*Math.pow(ampPhase1[i][j][0],k);
     result[i][j][1]= ampPhase0[i][j][1]*(1.0-k) + ampPhase1[i][j][1]*k;
  }
  return result;  
}


/* convert amplitude/phase into complex data and then inverse FFT */
public double [] ampPhaseToSpace(double [][][] ampPhase,
                                       boolean centered ) { // true if result should have 0:0 in the center
  double [][][] fft=new double[ampPhase.length][ampPhase[0].length][2];
  int i,j;
  for (i=0;i<ampPhase.length;i++) for (j=0;j<ampPhase[0].length;j++) {
    fft[i][j][0]=ampPhase[i][j][0]*Math.cos(ampPhase[i][j][1]);
    fft[i][j][1]=ampPhase[i][j][0]*Math.sin(ampPhase[i][j][1]);
  }
  double [] result=FFTHalf2FHT(fft,fft[0].length);
  fht_instance.inverseTransform(result);
  if (centered) fht_instance.swapQuadrants(result);
  return result;
}

/* convert amplitude/phase into FHT */
public double [] ampPhaseToFHT(double [][][] ampPhase) {
  double [][][] fft=new double[ampPhase.length][ampPhase[0].length][2];
  int i,j;
  for (i=0;i<ampPhase.length;i++) for (j=0;j<ampPhase[0].length;j++) {
    fft[i][j][0]=ampPhase[i][j][0]*Math.cos(ampPhase[i][j][1]);
    fft[i][j][1]=ampPhase[i][j][0]*Math.sin(ampPhase[i][j][1]);
  }
  double [] result=FFTHalf2FHT(fft,fft[0].length);
  return result;
}
 public double [] interpolateFHT (double [] fht0,    // first FHT array
                                    double [] fht1,  // second FHT array
                                    double  k        // interpolation coefficient ( 0.0 - fht0, 1.0 - fht1 )
                                    ){
    double [] points={k};
    double [][]result=interpolateFHT (fht0,    // first FHT array
                                      fht1,    // second FHT array
                                    points,   // array of interpolation points - 0.0 - fht0, 1.0 - fht1
                                      true);   // array of interpolation points - 0.0 - fht0, 1.0 - fht1
    return result[0];
 }


/* returns array of interpolated FHTs between fht0 and fht1, endpoints if present (0.0, 1.0) are referenced, not cloned */
 public double [][] interpolateFHT (double [] fht0,    // first FHT array
                                    double [] fht1,    // second FHT array
                                  double [] points){   // array of interpolation points - 0.0 - fht0, 1.0 - fht1
     return interpolateFHT (fht0,    // first FHT array
                            fht1,    // second FHT array
                          points,   // array of interpolation points - 0.0 - fht0, 1.0 - fht1
                           false);   // do not clone 0.0 and 1.0 (ends)

 }
 public double [][] interpolateFHT (double [] fht0,    // first FHT array
                                    double [] fht1,    // second FHT array
                                  double [] points,   // array of interpolation points - 0.0 - fht0, 1.0 - fht1
                               boolean cloneTrivial   ){
    double [] fht_div=fht_instance.divide(fht1,fht0);
    int size=(int) Math.sqrt(fht0.length);
    int hsize=size/2;
    double [][][] aphase= new double[hsize+1][size][2];
    double [][][] amp01=  new double[hsize+1][size][2]; /* squared amplituides of fht0 and fht1 */
    double [][]   phase=  new double[hsize+1][size]; /* +/-pi phase of the first array */
    double[][][]fft0=    FHT2FFTHalf (fht0, size);
    double[][][]fft1=    FHT2FFTHalf (fht1, size);
    double[][][]fft_div= FHT2FFTHalf (fht_div, size);
    int i,j,k;
    double a,c,p;
/* use mul for amplitudes, div - for phases */
    for (i=0;i<=hsize;i++) for (j=0;j<size;j++) {
      amp01[i][j][0]= fft0[i][j][0]*fft0[i][j][0]+fft0[i][j][1]*fft0[i][j][1];
      amp01[i][j][1]= fft1[i][j][0]*fft1[i][j][0]+fft1[i][j][1]*fft1[i][j][1];
      if (amp01[i][j][0]>0.0) phase[i][j]=Math.atan2(fft0[i][j][1],fft0[i][j][0]);
      else phase[i][j]=0.0;
      aphase[i][j][0]=amp01[i][j][0]*amp01[i][j][1]; // product of squared amplitudes is OK for phase restorations, we just need to know where amplitudes are higher
      a=              fft_div[i][j][0]*fft_div[i][j][0]+fft_div[i][j][1]*fft_div[i][j][1];
      if (a>0.0) aphase[i][j][1]=Math.atan2(fft_div[i][j][1],fft_div[i][j][0]);
      else aphase[i][j][1]=0.0;
    }
    aphase[0][0][1]=0.0;
/* calculate full phases */    
    fullPhase(aphase);
    double [][]result=new double[points.length][];
    for (k=0;k<result.length;k++) {
       if      (points[k]==0.0) result[k]=cloneTrivial?fht0.clone():fht0;
       else if (points[k]==1.0) result[k]=cloneTrivial?fht1.clone():fht1;
       else { /* interpolate */
         c=points[k];
         for (i=0;i<=hsize;i++) for (j=0;j<size;j++) {
           if ((amp01[i][j][0]==0.0) || (amp01[i][j][1]==0.0)) a=0.0;
/* Extrapolation is defined here only in the direction of decreasing of the spectral amplitudes (outside, to the wider PSF), so additional limit to prevent division of small values */
/* Seems to wor, possible improvements: 1-filter spectr in high-freq areas. 2 - use farther inner points fro farther approximation */

           else if ((c<0.0) && (amp01[i][j][0]>amp01[i][j][1])) a=Math.sqrt(amp01[i][j][0]);
           else if ((c>1.0) && (amp01[i][j][0]<amp01[i][j][1])) a=Math.sqrt(amp01[i][j][1]);
           else a= Math.pow(amp01[i][j][0],0.5*(1.0-c))*Math.pow(amp01[i][j][1],0.5*c);
           p= phase[i][j]+c*aphase[i][j][1];
           fft0[i][j][0]=a*Math.cos(p);
           fft0[i][j][1]=a*Math.sin(p);
         }
         result[k]=FFTHalf2FHT(fft0,size);
       }
    }
    return result;  
 }



/* perform direct FFT on the square pixel array and then convert result to amplitude/full phase */
public double [][][] spaceToAmplPhase(double [] pixels,
                                      boolean centered ) { // true if result should have 0:0 in the center
  int size = (int) Math.sqrt(pixels.length);
  double [] fht=pixels.clone();
  if (centered) fht_instance.swapQuadrants(fht);
  fht_instance.transform(fht);
  double[][][] cmplx= FHT2FFTHalf (fht, size);
  return amplPhase(cmplx);
}


/* Convert complex FFT (half) into amplitude+full phase by flooding from the 0,0 and searching for the smallest phase increment. Are loops possible? */
public double [][][] amplPhase(double [][][] fft) {
  int size = fft[0].length;
  int hsize=fft.length-1;
  double [][][] aphase= new double[hsize+1][size][2];
  boolean [][]  map=    new boolean[hsize+1][size];
  int i,j;
/* calculate amplitudes and -pi..+pi phases */
  for (i=0;i<=hsize;i++) for (j=0;j<size;j++) {
    aphase[i][j][0]=Math.sqrt(fft[i][j][0]*fft[i][j][0]+fft[i][j][1]*fft[i][j][1]);
    if (aphase[i][j][0]>0.0) aphase[i][j][1]=Math.atan2(fft[i][j][1],fft[i][j][0]);
    else aphase[i][j][1]=0.0;
    map[i][j]=false;
  }
  aphase[0][0][1]=0.0;
  fullPhase(aphase);
  return aphase;
}

/* replace +/-pi  phase with the full phase, using amplitude to guide grouth of the covered area, so amplitude and phase does not need to be a pair from the same FFT array */
public void fullPhase(double [][][] aphase) {
  int size = aphase[0].length;
  int hsize=aphase.length-1;
  boolean [][]  map=    new boolean[hsize+1][size];
  int i,j;
  aphase[0][0][1]=0.0;
  int ix,iy,ix1,iy1,ix1n,iy1n;
  List <Integer> pixelList=new ArrayList<Integer>(100);
  Integer Index;
  int [][] dirs={{-1,0},{-1,-1},{0,-1},{1,-1},{1,0},{1,1},{0,1},{-1,1}};
  ix=0;
  iy=0;
//  int clusterSize=0;
  boolean noNew=true;
  int maxX=0;
  int maxY=0;
  int oldX=0;
  int oldY=0;
  boolean oldConj=false;
  int listIndex;
  Index=iy*size + ix;
  pixelList.clear();
  pixelList.add (Index);
//  clusterSize++;
  map[iy][ix]=true;
  noNew=true;
  double phase, oldPhase, fullCyclesPhase;
  double maxValue=-1.0;

  while (pixelList.size()>0) {
/* Find maximal new neighbor */
    maxValue=-1.0;
    listIndex=0;
    while (listIndex<pixelList.size()) {
      Index=pixelList.get(listIndex);
      iy=Index/size;
      ix=Index%size;
      noNew=true;
      for (j=0;j<8;j++) {
        ix1=(ix+dirs[j][0]+size) % size;
        iy1=(iy+dirs[j][1]+size) % size;
        if ((iy1>hsize) || (((iy1==0) || (iy1==hsize)) && (ix1> hsize))) {
          ix1n=(size-ix1)%size;
          iy1n=(size-iy1)%size;
        } else { /* But phase will be opposite sign */
          ix1n=ix1;
          iy1n=iy1;
        }
        if (!map[iy1n][ix1n]) {
          noNew=false;
          if (aphase[iy1n][ix1n][0]>maxValue) {
            maxValue= aphase[iy1n][ix1n][0];
            maxX=ix1n;
            maxY=iy1n;


//    if (DEBUG_LEVEL>4)      System.out.println(" amplPhase(): iy="+iy+ " ix="+ix+" maxY="+maxY+" maxX="+maxX);

          }
        }
      }
      if (noNew) pixelList.remove(listIndex);  //  remove current list element
      else       listIndex++;     // increase list index
    }


    if (pixelList.size()==0) break;
/* To calculate the phase - find already processed neighbor with the highest amplitude */
    maxValue=-1.0;
    for (j=0;j<8;j++) {
      ix1=(maxX+dirs[j][0]+size) % size;
      iy1=(maxY+dirs[j][1]+size) % size;
      if ((iy1>hsize) || (((iy1==0) || (iy1==hsize)) && (ix1> hsize))) {
        ix1n=(size-ix1)%size;
        iy1n=(size-iy1)%size;
      } else { /* But phase will be opposite sign */
        ix1n=ix1;
        iy1n=iy1;
      }
//      if (DEBUG_LEVEL>5) System.out.println(" amplPhase()-- : iy1="+iy1+ " ix1="+ix1+" map["+iy1n+"]["+ix1n+"]="+map[iy1n][ix1n]);
      if (map[iy1n][ix1n]) {
        if (aphase[iy1n][ix1n][0]>maxValue) {
          maxValue= aphase[iy1n][ix1n][0];
          oldX=ix1n;
          oldY=iy1n;
          oldConj=(iy1!=iy1n) || (ix1!=ix1n); // point on the other half (conjugate)
//    if (DEBUG_LEVEL>3)      System.out.println(" amplPhase(): iy1="+iy1+ " ix1="+ix1+" oldY="+oldY+" oldX="+oldX);

        }
      }
    }
//    if (DEBUG_LEVEL>4)   System.out.println(" amplPhase():Old:["+oldConj+"] ("+pixelList.size()+")  "+oldX+":"+oldY+" "+IJ.d2s(aphase[oldY][oldX][0],2)+":"+ IJ.d2s(aphase[oldY][oldX][1],2)+
//                                     " New:"+maxX+":"+maxY+" "+IJ.d2s(aphase[maxY][maxX][0],2)+":"+ IJ.d2s(aphase[maxY][maxX][1],2)+
//                                     " Diff="+IJ.d2s((aphase[maxY][maxX][1]-(oldConj?-1:1)*aphase[oldY][oldX][1]),2));

/* Calculate the phase from the closest neighbor */
    oldPhase=(oldConj?-1:1)*aphase[oldY][oldX][1];
    fullCyclesPhase=2*Math.PI*Math.floor(oldPhase/(2*Math.PI)+0.5);
    oldPhase-=fullCyclesPhase; // +/- pi
    phase=aphase[maxY][maxX][1];
    if      ((phase - oldPhase) > Math.PI) fullCyclesPhase-=2*Math.PI;
    else if ((oldPhase - phase) > Math.PI) fullCyclesPhase+=2*Math.PI;
    aphase[maxY][maxX][1]+=fullCyclesPhase;


    if (DEBUG_LEVEL>3) {
      System.out.println(" amplPhase():Old:["+oldConj+"] ("+pixelList.size()+")  "+oldX+":"+oldY+" "+IJ.d2s(aphase[oldY][oldX][0],2)+":"+ IJ.d2s(aphase[oldY][oldX][1],2)+
                                     " New:"+maxX+":"+maxY+" "+IJ.d2s(aphase[maxY][maxX][0],2)+":"+ IJ.d2s(aphase[maxY][maxX][1],2)+
                                     " Diff="+IJ.d2s((aphase[maxY][maxX][1]-(oldConj?-1:1)*aphase[oldY][oldX][1]),2));
    }


//IJ.d2s(1000*patternCorr[7]/(FFTSize/2),5)+

/* Add this new point to the list */
    Index=maxY*size + maxX;
    pixelList.add (Index);
//    clusterSize++;
    map[maxY][maxX]=true;
  } // end of while (pixelList.size()>0)
/* Fix remaining phases for y=0 and y=hsize */
  for (i=1;i<hsize;i++) {
    aphase[0][size-i][1]=-aphase[0][i][1];
    aphase[hsize][size-i][1]=-aphase[hsize][i][1];
  }
//  return aphase;
}
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

// -------------------------------------
}
