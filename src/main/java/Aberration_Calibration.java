/**
 ** -----------------------------------------------------------------------------**
 ** Aberration_Calibration.java
 **
 ** Measurement of the aberrations (array of PSF), preparation of the calibration files
 ** 
 **
 ** Copyright (C) 2010-2011 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  Aberration_Calibration.java is free software: you can redistribute it and/or modify
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
import ij.process.*;
import ij.gui.*;
import ij.plugin.frame.*;
import ij.text.TextWindow;

import java.awt.*;
import java.awt.event.*;
import java.io.*; // FIle
import java.util.Properties;

import javax.swing.*; // TODO: modify methods that depend on it, use class CalibrationFileManagement

import java.util.*; 
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;








//import FocusingField.FocusingFieldMeasurement;
//import FocusingField.MeasuredSample;
import Jama.Matrix;  // Download here: http://math.nist.gov/javanumerics/jama/

public class Aberration_Calibration extends PlugInFrame implements ActionListener {
	private static final long serialVersionUID = 1040236897357482595L;
	private Panel panel1,panel2,panel3,panelDirs,panelLens;
//	private Panel panelFlat;
	private Panel panelConf1,panelConf2,panelRun, panelDistortions, panelFitDistortions, panelProcessDistortions,panelAberrations ;
	private Panel panelCorrectGrid;
	private Panel panelFocusing,panelFocusing1;
	private Panel panelCurvature;
	private Panel panelGoniometer;
	private Panel panelPixelMapping, panelStereo,panelStereo1;
	
	private showDoubleFloatArrays SDFA_INSTANCE; // just for debugging?
	JP46_Reader_camera JP4_INSTANCE;
	static Frame instance;
	static Properties PROPERTIES=new Properties();
	public static int     DEBUG_LEVEL =     1;
	public static int     MASTER_DEBUG_LEVEL = 1;
	public static boolean SHOW_AS_STACKS= true; //  Show debug images as stacks (false - individual)

	public static int     FFT_SIZE=       256;
	public static int     MAP_FFT_SIZE=    64; // used to find where grid covers the image
	public static double  GAUSS_WIDTH=    0.4; //0 - use Hamming window - initWindowFunction()
	public static int     FFT_OVERLAP=     32; // createPSFMap()
	public static int     PSF_SUBPIXEL=     4; // sub-pixel decimation 
	public static boolean PSF_SAVE_FILE= true; // save PSF array to a multi-slice TIFF file
	public static int     THREADS_MAX=    100; // testing multi-threading, limit maximal number of threads
	public static boolean UPDATE_STATUS= true; // update ImageJ status info
    public static MatchSimulatedPattern matchSimulatedPattern=null; //=new MatchSimulatedPattern();
    public static PixelMapping PIXEL_MAPPING=null;
    public static PolynomialApproximation polynomialApproximation=new PolynomialApproximation(2); // just to force recompile
    public static int     LAST_FRAME_NUMBER=-1;
	public static SimulationPattern.SimulParameters SIMUL = new SimulationPattern.SimulParameters (
			1024, //public static int SIMUL.patternSize=512; // size of the side of the square pattern bitmap 
			1,      // pattern_type:  0 - vanilla linear, 1 - curved,
			2.0,    // pattern_modifier
			0.241,  // freq_x1 (these particular parameters - sample only, not used in real calculations)
			0.0055, // freq_y1
			0.0,    // phase1 (radians)
			-0.0001,// freq_x2;
			0.0231, // freq_y2;
			0.0,    // phase2;
			4,      // subdiv - subdivide pixels in each direction
			0.75,   // fill - part of the (center) pixel area being "photosensitive"
			true,   // center_for_g2: align pattern to phases for the diagonal (both greens) sub-array
			0.005,    // 2.0,    // bPatternSigma, // blur bPattern with this sigma
			0.1,     // barraySigma // blur barray with this sigma, multiplied by subdiv
			0.1,     // smallestSubPix, // subdivide pixels down to that fraction (linear) when simulating
			0.05,    // bitmapNonuniforityThreshold // subdivide pixels until difference between the corners is below this value
			0.5,     // offsetX, // debug - add to X during simulation, in pixels
			0.5      // offsetY // debug - add to Y during simulation, in pixels
	);

	public static EyesisAberrations.InterpolateParameters INTERPOLATE= new EyesisAberrations.InterpolateParameters (
			32, // insize -      size of each kernel (should be square)
			4,  // step   -      number of subdivisions from input to output
			17, // add_top     - add this number of kernel rows to the output above the existent/interpolated
			17, // add_left    - add this number of kernel columns to the output on the left of the existent/interpolated
			18, // add_right   - add this number of kernel columns to the output on the right of the existent/interpolated
			17, // add_bottom  - add this number of kernel rows to the output below the existent/interpolated
			0.0 // extrapolate - 0.0 - duplicate, 1.0 - extrapolate outside of the known kernels
	);

	public static EyesisAberrations.InverseParameters INVERSE = new EyesisAberrations.InverseParameters  (
			32,    // dSize - size (side of square) of direct PSF kernel
			64,    // rSize - size (side of square) of reverse PSF kernel
			0.007, // invertRange - when FFT component is less than this fraction of the maximal value, replace 1/z with Z
			0.9,   // otfCutoffEnergy - use frequency points that have INVERSE.otfCutoffEnergy of the total to determine ellipse for limiting frequency responce
			3.0,   // otfEllipseScale - size of elliptical window relative to the cluster defined by INVERSE.otfCutoffEnergy
			true,  // otfEllipseGauss
			0.9,   // psfCutoffEnergy - Limit result kernel to proportional of the PSF, calculate initial cluster shape by this cutoff energy
			2.5,   // psfEllipseScale - size of elliptical window to limuit reverse PSF as proportional to direct one
			0.01,  // rpsfMinMaskThreshold -completely zero reversed kernel elements where elliptical mask is below this threshold
			true,  // filter - apply variable-sigma filtering to the inverted PSF
			2.0,   // blurIndividual;
			2.0,   // blurDiagonal;
			1.6,   // blurChecker;
			2.0,   // gaussianSigmaIndividual;
			2.0,   // gaussianSigmaDiagonal;
			1.6,   // gaussianSigmaChecker;
			0.8,   // sigmaScale - reduce variable sigma in the center from uniform one
			0.2,   // sigmaToRadius - variable blurring - sigma will be proportional distance from the center
			true,  // filterDirect - apply variable blurring of the direct kernels (before inversion)
			0.4,   // sigmaScaleDirect - reduce variable sigma in the center from uniform one (defined for inverted kernels above);
			0.2    // public double sigmaToRadiusDirect increase blurring farther form the center;

	);
	public static EyesisAberrations.PSFParameters PSF_PARS = new EyesisAberrations.PSFParameters (
			0.85,  // minContrast - minimal instance contrast to use in binning (compared to the one at [0,0]
			0.5,   // windowFrac  - reduce the PSF cell size to this part of the area connecting first negative clones
			true,  // useWindow   - multiply separated OTF instance by window function (Hamming or Gaussian)
			false, // symm180     - force PSF center-symmetrical (around centroid that is defined by lateral chromatic aberration
			false, // ignoreChromatic - ignore lateral chromatic aberration (center PSF to 0,0)
			0.2,   // smoothSeparate  - low pass filter width when separating individual PSF instances
			0.75,  // topCenter - consider only points above this fraction of the peak to find the centroid
			0.0,   // sigmaToRadius - variable-sigma blurring to reduce high frequencies more for the pixels farther from the PSF center
			0.8,   // wingsEnergy - fraction of energy in the pixels to be used (where is it used?)
			4.0,   // wingsEllipseScale - increase wings cutoff ellipse by this from one defined by the  cutoff energy
			.75,   // minDefinedArea // minimal (weighted) fraction of the defined patter pixels in the FFT area
			false, //boolean approximateGrid; // approximate grid with polynomial
			true, // was FALSE before centerPSF
			1.0,   // psfParameters.mask1_sigma,
			0.25,   // psfParameters.mask1_threshold,
			1.0,   // psfParameters.gaps_sigma,
			0.25    // mask_denoise - somewhat depends on psfParameters.mask1_threshold

	);


	public static EyesisAberrations.OTFFilterParameters OTF_FILTER = new  EyesisAberrations.OTFFilterParameters(
			0.008, // deconvInvert - when FFT component is less than this fraction of the maximal value, replace 1/z with Z
			2.0,   // zerofreqSize - used for filtering oversampling artifacts - size of zero freq maximum (if absent on simulated model PS)
			2.5,   // smoothPS - smooth model PS for rejecting aliases (0 - no smooth, >0 additional Gauss before FFT smaller than normal by this ratio)
			0.02,  // thresholdHigh -used for filtering oversampling artifacts - relative to max PS value to make filter completely rejecting
			0.002  // thresholdLow - used for filtering oversampling artifacts - relative to max PS to make filter completely transmissive
	);

	public static MatchSimulatedPattern.PatternDetectParameters PATTERN_DETECT = new MatchSimulatedPattern.PatternDetectParameters (
			0.4, // GAUSS_WIDTH=    0.4; //0 - use Hamming window - initWindowFunction()
			0.2, // corrGamma - pattern detection: gamma applied to PS
			1.5, // corrSigma - pattern detection: high-pass filter (0.0 - none) gamma(PS)
			2,   // diffSpectrCorr - maximal distance between maximum on spectrum and predicted maximum on autocorrelation of gamma(|spectrum|)
			0.0, // shrinkClusters - Shrink clusters by this ratio (remove lowest) after initial separation
			4,   // multiplesToTry - try this number of m0.300aximums proportionally farther from 0,0 than the two closest (increase precision)
			1.0, // deviation - when looking for maximums - maximal distance from predicted from the lower order one
			6,   // deviationSteps -maximal iterations when looking for local maximum
			1.5, // highpass - model correlation high-pass filter (relative to pattern fundamental frequency - average of 2)
			0.4, // corrRingWidth - ring (around r=0.5 dist to opposite corr) width , center circle r=0.5*PATTERN_DETECT.corrRingWidth
			5.0, // minCorrContrast -  discrimination threshold between good and bad pattern correlation
			0.0, // minGridPeriod
			0.0  // maxGridPeriod
	);
	
	
	public static EyesisAberrations.ColorComponents COMPONENTS = new EyesisAberrations.ColorComponents(
			false, // green1 (colorsToCorrect[0])
			true,  // red
			true,  // blue
			false, // green2
			false, // diagonal
			true,  // checker
			5,     // referenceComponent (checker)
			true   // equalizeGreens
	);
	public static ShowResults SHOW_RESULTS =new ShowResults (
			false, //  showPSF -  show combined PSF kernels (per-color and/or composite, as defined in SHOW_INDIVIDUAL, SHOW_COMPOSITE)
			false, //  showMTF -  calculate/show MTF (see notes to SHOW_RESULTS.showPSF)
			false, //  showInverted - show inverted kernels (unfiltered), same notes
			true,  //  showFiltered -  filter and show inverted kernels
			true   //  showGaussians -  create Gaussian kernels with the same centers as inverted ones (low noise, use in low details areas)
	);
	public static EyesisAberrations.MultiFilePSF MULTIFILE_PSF = new EyesisAberrations.MultiFilePSF(
			0.025, // overexposedMaxFraction -  allowed fraction of the overexposed pixels in the PSF kernel measurement area 
			0.5,  // public double  weightOnBorder=0.01;  
			0.2, //0.05,  // radiusDiffLow; // do not remove partial kernel cell if radius differs from average less than by this fraction
			0.25,   // radiusDiffHigh;  // remove this cell even if it is the only one
			1.0,    // shiftToRadiusContrib; // Center shift (in pixels) addition to the difference relative to radius difference (in pixels)
			2.0,    // sharpBonusPower; // increase weight of the "sharp" kernels by dividing weight by radius to this power
			0.1,    // maxFracDiscardWorse=0.1; // discard up to this fraction of samples that have larger radius (i.e. falling on the target seam that may only make PSF larger)
			0.5,    // maxFracDiscardAll=0.5; // continue removing outlayers (combined radius and shift), removing not more that this fraction (including maxFracDiscardWorse)
			0.8, //,   // internalBonus // cell having 8 around will "seem" twice better than having none (radiusDiff* twice higher)
			0.75,  // validateThreshold -       fraction of full PSF "energy"
			false, // validateShowEllipse -     show ellipse parameters of partial PSF arrays
			true,  // showWeights -             show image indicating frame coverage
			false  // fill missing PSF kernels from nearest existent ones

	);
	
	public static ProcessCalibrationFilesParameters PROCESS_PARAMETERS = new ProcessCalibrationFilesParameters(
		"jp46",              // sourceFileExtension, 
		"tiff",              // kernelFileExtension, 
		"PSF_",              // kernelFilePrefix
		"PSF-RAW_",          // psfRawPrefix,
		"PSF-INTERPOLATED_", // psfInterpoaltedPrefix,
		"RPSF_",             // rpsfPrefix,
		"GAUSSIAN_",         // gaussianPrefix,
		"", //"source",            // sourceSuperDirectory; // having 1-1, 1-2,.. 3-3 subdirs
		"", //"partial_kernels",   // partialKernelsSuperDirectory; // having 1-1, 1-2,.. 3-3 subdirs
		"", //"kernels",           //  kernelsDirectory,              // results with Gaussian and deconvolution kernels for all channels

		true,                // processAllChannels;           // if true - process all channels, otherwise only enabled in processChannels[]
		false,               // processChannels[0] 1-1
		false,               // processChannels[1] 1-2
		false,               // processChannels[2] 1-3
		false,               // processChannels[3] 2-1
		false,               // processChannels[4] 2-2
		false,               // processChannels[5] 2-3
		false,               // processChannels[6] 3-1
		false,               // processChannels[7] 3-2
		false,               // processChannels[8] 3-3
		true,                // keepOld;                   // do not re-calculate existent partial kernels, only the new ones
		true,  //false       // selectFiles                // select individual files to process
		
		true,                // processSourceImages,           // process source calibration files		
		true,                // combinePSFfiles,               // combine partial PSF kernels 
		true,                // interpolatePSFkernel,          // interpolate PSF kernels (fail if missing??)
		true,                // invertKernels,                 // invert interpolated kernels
		true,                // gaussianKernels                // create Gaussian kernels
		true,                // useXML;                        // save/restore settings as xml file
		true                 // saveSettings;                  // save current settings in results directory
		
    );

    public static  FlatFieldParameters FLATFIELD_PARAMETERS = new FlatFieldParameters (
    		"jp4",
		    true,   // normalize,
		 // to filter good/bad images in each sub-band    	
	    	0.99,   // overExpValue= 0.99;
	    	0.025,  // overExpFrac=  0.025;
	    	0.25,   // underExpValue=0.25;
	    	0.5,    // underExpFrac= 0.5;
	    	true,   // noTiltEdges = do not apply tilt to edge sections - to get closer to the corners
	    	0,      // functionType: 0 polynomial, 1 - power
	    	1,      // functionModifier - additional function modifier
	       	0.75,   // section34=0.5; // location of 4-th and 5-th section (ratio from o to 1 and from 0 to 2 (1-3-0-4-2) 

            -1.0, //0.5,    // centerWeight, //("-1" all same weight) weight for the error function will be proportional to r^2 (r - half smallest dimension) plus weight in the center
	    	true,   // LM_auto=true;   // automatically iterate (false open - dialogs)
	    	0.001,  // LM_lambdaInitial=0.001;
	    	8.0,    // LM_lambdaStepUp=   8.0; // multiply lambda by this if result is worse
	    	0.5,    // LM_lambdaStepDown= 0.5; // multiply lambda by this if result is better
	    	0.0001, // LM_thresholdFinish=0.0001; // stop iterations if 2 last steps had less improvement (but not worsening ) 
	    	100,    // LM_numIterations=  100; // maximal number of iterations 

	    	3.0,    // fatZero
		    64,     // margin_left,
		    64,     // margin_right,
		    64,     // margin_top,
		    64,     // margin_bottom,
		    2,      // decimate,
		    16,     // sampleWidth,
		    256.0,  // highPassSigma,
	    	0.25,   //  maxTilt,       // real life 0.16

			"",     // flatFieldDirectory,           // results with flat field calibration files
			true,   // eyesisMode
			true,   // processAllChannels,           // if true - process all channels, otherwise only enabled in processChannels[]       
			false,  // processChannels11,
			false,  // processChannels12,
			false,  // processChannels13,
			false,  // processChannels21,
			false,  // processChannels22,
			false,  // processChannels23,
			false,  // processChannels31,
			false,  // processChannels32,
			false,  // processChannels33,
			true,   // useXML,
			true    // saveSettings
	);

public static MatchSimulatedPattern.DistortionParameters DISTORTION =new MatchSimulatedPattern.DistortionParameters(
		  64, //32, // use 64 for less artifacts, // correlationSize
		  0.75,// reduce to 0.5 when correlationSize==64 // correlationGaussWidth
		  false, // boolean absoluteCorrelationGaussWidth=false; // do not scale correlationGaussWidth when the FFT size is increased  
		  0, //zeros - // leave this number of zeros on the margins of the window (toatal from both sides). If correlationGaussWidth>0 will 
	        // additionally multiply by Hamming
		  128, // FFT size
		  0.5, //fftGaussWidth
		  0.0, //phaseCorrelationFraction
		  1.5, // 2.5, //6.0, // 2.0, // 0.0, // correlationHighPassSigma, - pixels in frequency domain
		  0.6, //2.0, //0.5, //0.0, //  correlationLowPassSigma, - fraction of the frequency range
		  0.4,  // correlationRingWidth- ring (around r=0.5 dist to opposite corr) width , center circle r=0.5*PATTERN_DETECT.corrRingWidth
		  3.0, //  correlationMaxOffset,     // maximal distance between predicted and actual pattern node
		  2.0, // increase back to .5? was needed with fisheye. 5.0, //	double correlationMinContrast,   // minimal contrast for the pattern to pass
		  2.5, // correlationMinInitialContrast,   // minimal contrast for the pattern of the center (initial point)
		  150,  // minimalPatternCluster minimal pattern cluster size (0 - disable retries)
		  2.0, // scaleMinimalInitialContrast increase/decrease minimal contrast if initial cluster is >0 but less than minimalPatternCluster
		  0.5, //  when searching for grid, step this amount of the FFTSize
          4, // public int    patternSubdiv;
		  0.0, // double correlationDx; // not saved
		  0.0,  // double correlationDy; // not saved
		  300, // gridSize
		  0, //2,    // loop_debug_level
		  true, // refineCorrelations
		  true, // fastCorrelationOnFirstPass, // use fast (less precise) correlation on first pass
		  false, // fastCorrelationOnFinalPass, // use fast (less precise) correlation on refine pass
		  0.02, // bPatternSigma; // overwrites SimulationParameters.bPatternSigma
		  0.5,  // barraySigma // blur barray with this sigma, multiplied by subdiv
		  2.5,  // 0.8,	// correlationWeightSigma, // sigma (in pixels) for maximum approximation 
		  2.0,  // 2.0	// correlationRadiusScale // maximal radius to consider, in sigmas (if 0 - use sigma as radius)
		  2,    // 6,   //public int    correlationRadius;    // radius (green pixel) of the correlation maximum to use for x/y measurement
		  0.8,   // double correlationThreshold; // fraction of the value of the maximum for the point to be included in centroid calculation
		  16, //	public int    correlationSubdiv;    // Total subdivision of the correlation maximum (linear and FFT)
		  4, //1 // 4 // correlationFFTSubdiv
	      true,// correlationAverageOnRefine, // average position between neighbor samples
	      false, // boolean refineInPlace;       // Update coordinates of the grid points as they are recalculated (false - then update all at once)
	      0.5, // averageOrthoDist,     // distance to up/down/right left neighbors (0.5)
	      0.4, //averageOrthoWeight,   // weight of 4 ortho neighbors (combined) - 0.4), weight of center -s 1.0-averageOrthoWeight-averageDiagWeight
		  0.5, // averageDiagDist,     // distance to diagonal neighbors (projection on x/y) (0.5)
	      0.4, //averageDiagWeight   // weight of 4 diagonal neighbors (combined) - 0.4)
	      true, // useQuadratic
	      true,//  boolean removeLast,         // remove outer (unreliable) row of nodes
	      3, //  int    numberExtrapolated  // add this number of extrapolated nodes
	      4.0, //  public double extrapolationSigma;  // use instead of the correlationWeightSigma during final extrapolation
	      1.2,//  double minUVSpan           // Minimal u/v span in correlation window that triggers increase of the correlation FFT size
	      true, // boolean flatFieldCorrection,  // compensate grid uneven intensity (vignetting, illumination)
	      1.0,  // double flatFieldExtarpolate, // extrapolate flat field intensity map (relative to the average grid period)
	      1.0,  // double flatFieldBlur,        // blur the intensity map (relative to the average grid period)
	      0.1,  // double flatFieldMin;    // do not try to compensate if intensity less than this part of maximal
	      1.0,  // double flatFieldShrink=1.0;     // Shrink before extrapolating intensity map (relative to the average grid period) 
	      4.0,  // double flatFieldExpand=3.0;     // Expand during extrapolation (relative to the average grid period)
	      1.0,  // double flatFieldSigmaRadius=1.0;// Extrapolation weight effective radius (relative to the average grid period)
	      1.5,  // double flatFieldExtraRadius=1.5;// Consider pixels in a square with the side twice this (relative to flatFieldSigmaRadius)
	      2.0,  // multiply the average grid period to determine the area for averaging the grig brightness
	      false //legacyMode
		);
    static int [] viewMap={0,0,0,0, 0,0,0,0,
    	                   0,0,0,0, 0,0,0,0,
    	         	       0,0,0,0, 0,0,0,0,
    	         	       1,1};
	public static Distortions LENS_DISTORTIONS;
  
    public static Distortions.PatternParameters PATTERN_PARAMETERS=new Distortions.PatternParameters(
    		viewMap,
    		1, // initial number of stations
    		3022.6, // double patternWidth;  // pattern full width in mm
    		2667.0, // double patternHeight; // pattern full height in mm
    		41.6667, // patternHalfPeriod;    // distance between opposite sign nodes
    		5.0     // double patternTilt;   // pattern tilt (degrees) - U clockwise from X-right (V clockwise from Y-down)
    		);
//    public static Distortions.LensDistortionParameters LENS_DISTORTION_PARAMETERS=new Distortions.LensDistortionParameters(

//    public static Distortions.LensDistortionParameters LENS_DISTORTION_PARAMETERS=LENS_DISTORTIONS.new LensDistortionParameters(
  	  public static Distortions.LensDistortionParameters LENS_DISTORTION_PARAMETERS=(new Distortions()).new LensDistortionParameters(
    		4.5, // double focalLength
    		2.2, // double pixelSize (um)
    		2.8512, //double distortionRadius mm - half width of the sensor
    		0.0, //public double distortionA8=0.0; // r^8 (normalized to focal length or to sensor half width?)
    		0.0, //public double distortionA7=0.0; // r^7 (normalized to focal length or to sensor half width?)
    		0.0, //public double distortionA6=0.0; // r^6 (normalized to focal length or to sensor half width?)
    		0.0, //public double distortionA5=0.0; // r^5 (normalized to focal length or to sensor half width?)
    		0.0, // double distortionA // r^4 (normalized to focal length or to sensor half width?)
    		0.0, // double distortionB // r^3
    		0.0, // double distortionC // r^2
    		// orientation/position parameters
    		0.0, // double yaw;    // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - clockwise from top
    		0.0, // double pitch   // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - up
    		0.0, // double roll    // angle in degrees rotation around camera optical axis (perpendicular to pattern if yaw==0, pitch==0), positive - clockwise
    		0.0, // double x0      // lens axis from pattern center, mm (to the right)
    		0.0, // double y0      // lens axis from pattern center, mm (down)
    		0.0, // double z0      // lens axis from pattern center, mm (away from the camera perpendicular to the patter plane)
    		2360,// double distance// distance from the lens input pupil to the pattern plane along the camera axis, mm 
    		1296.0, // double px0     // lens axis from sensor, horizontal, from left (pixels)
    		968.0, // double py0     // lens axis from sensor, vertical, from top (pixels)
    		true//  boolean flipVertical // acquired image is mirrored vertically (mirror used)
    		);
//    public static double [] defaultGoniometerPosition={0.0, 0.0, 2360};
    public static Distortions.EyesisCameraParameters EYESIS_CAMERA_PARAMETERS=new Distortions.EyesisCameraParameters(
    		1,    //int numStations,
    		true, //false, // boolean isTripod=false; // when true - make goniometerHorizontal rotation around "vertical" axis and "goniometerAxial" - around 
	    	0.0, // double goniometerHorizontal, // goniometer rotation around "horizontal" axis (tilting from the target - positive)
	    	0.0, // double goniometerAxial, // goniometer rotation around Eyesis axis (clockwise in plan - positive 
			1, // 26, // 1,   // int numSubCameras,
	    	0.0, // double interAxisDistance, // distance in mm between two goniometer axes
	    	0.0, //double interAxisAngle,    // angle in degrees between two goniometer axes minus 90. negative if "vertical" axis is rotated
	    	                            // clockwise when eyesis is in 'normal' position, looking to the target
	    	0.0, //double horAxisErrPhi,   // angle in degrees "horizontal" goniometer axis is rotated around target Y axis from target X axis (CW)
	    	0.0, //double horAxisErrPsi,   // angle in degrees "horizontal" goniometer axis is rotated around moving X axis (up)
	    	0.0,
	    	0.0,
	    	0.0, 0.0, 2360, //double [] GXYZ // coordinates (in mm) of the goniometer horizontal axis closest to the moving one in target system
	       	2592, // int sensorWidth=      2592;
	    	1936, //int sensorHeight=     1936;
	    	2,    //int    shrinkGridForMask=2; //shrink detected grids by one point for/vert this number of times before calculating masks
	    	-2.0, // double maskBlurSigma=    2.0;   // blur sensor masks (>0 - pixels, <0 - in grid units)
	    	4,    // int    decimateMasks
	    	0.1,  // double badNodeThreshold=0.1; // filter out grid nodes with difference from quadratically predicted from 8 neighbors in pixels
    		1,    // int  maxBadNeighb; // maximal number of bad nodes around the corrected one to fix
    		50,   // int minimalValidNodes
        	1,    // int   weightMultiImageMode=1; // increase weight for multi-image sets (0 - do not increase, 1 - multiply by number of images in a set)
        	1.0,  // public double  weightMultiExponent= 1.0; // if( >0) use grid diameter to scale weights of this image
        	1.0,  // public double  weightDiameterExponent=1.0;
        	1.0,  // public double weightYtoX=1.0; // relative Y-to-X errors weight (to somewhat compensate for rectabular shape of the sensor)
    		0.4,  //minimalGridContrast
	    	4.0,  // public double shrinkBlurSigma = 4.0;
	    	0.5,  // public double shrinkBlurLevel = 0.5;
	       -1.0, // double balanceChannelWeightsMode
	    	2.0, // public double removeOverRMS=2.0;            // error is multiplied by weight function before comparison (more permissive on the borders
	    	4.0 //public double removeOverRMSNonweighted=4.0; // error is not multiplied (no more permissions on tyhe borders
	    	);

    //
    static double [][] LASER_UV={
//    	{-30.5,-11.5}, // top left
//    	{26.5,-16.5},  // top right
//    	{-27.5,35.5},  // bottom left
//    	{32.5,29.5}};  // bottom right
// settings for the office wall    	
	{-30.5,-20.5}, // top left
	{26.5,-25.5},  // top right
	{-27.5,26.5},  // bottom left
	{32.5,20.5}};  // bottom right
    public static MatchSimulatedPattern.LaserPointer LASER_POINTERS= new MatchSimulatedPattern.LaserPointer (
	    1.06,     //	public double headLasersTilt=  1.06; // degrees, right laser lower than left laser
    	0.05,      // minimalIntensity
    	1.5,      // maximalIntensity
    	30,       //overexposedRadius
    	1.0,      //public double lowpassSigma;    // low pass sigma, in pixels
    	20.0,    //public double highpassSigma;   // high pass sigma, in pixels
    	0.8,      // public double headLowpassSigma;    // low pass sigma, in pixels for optical head lasers
    	1.0,      // quadraticScaleSigma; // find local maximum by quadratic intrepolating pixels around maximal value (relative to lkow pass sigma) 
		4,        // algorithmNumber
    	3,        // closestOffender
    	200,       // fartherstOffender
    	0.05,     // fatZero
    	0.6,      // double greenFloor;      // when dividing by green, add this fraction of maximal value (decrease green accordingly)
    	true,    // boolean useOther=false; // when true - use red and other color, when false - only red
    	true,     // boolean otherGreen=true; // other color is green (false - blue)
    	0.1,      // public double threshold;
    	false,    // public boolean swapUV; // first
    	false,    // public boolean flipU;
    	false,    // public boolean flipV;
    	true,     // public boolean whiteOnly; // verify laser is on the white pattern cell
    	0.6,      // public double  maxOffsetFromCenter; // maximal offset of the laser spot from the center, relative to cell radius
    	LASER_UV, // double [][] laserUVMap; // first index - number of pointer points
    	1.5,      // public double laserSignalToNoise=4.0; // Minimal signal-to-noise ratio for laser pointers
    	10,       // public double localMaxRadius=10; // sensor pix. currently uses just square (2*localMaxRadius+1)**2 
    	true,     // public boolean usePatternFilter=true; // Filter laser positions by likely pattern white cells
    	2,        // public int decimatePatternFilter=1; // reduce resolution for pattern filter
    	40.0,     // public double localContrastSigma=40; // use to calculate local level and contrast
    	0.8,      // public double localToGlobalContrast=0.8; // 0 - same contarst normalization fro the whole image, 1.0 - pure local
    	4.0,      // public double patternLowPassSigma=4.0; // filter normalized patetrn before thresholding 
    	0.2,      // public double patternThreshold=0.2; // fraction of dispersion (same positive for white cells, negative for black ones) 
    	30.0,     // public double maximalCellSize=30.0; // White cells should have black pixels in all 4 quadrants not farther than this
    	3,        // public int numPasses=3; // number of black/white alternations of the surrounding cells to use in quadrant filtering
    	false,    // public boolean bordersOK=false; // frame border as good cell for quadrant filter
    	0.1,      // public double blurredMaskThreshold=0.1; // select only areas with multiple pattern white cells
    	2.0,      // public double maskGrow=4.0; // grow final mask (pixels)              
    	1         // public int    debugLevel=1;
    );
    
    public static Distortions.RefineParameters REFINE_PARAMETERS = new Distortions.RefineParameters();
    public static Distortions.DistortionCalibrationData DISTORTION_CALIBRATION_DATA=null; 
//    public static Distortions.FittingStrategy FITTING_STRATEGY=null; 
    
//	public static boolean ADVANCED_MODE=false;
	public static boolean ADVANCED_MODE=true;
	public static boolean MORE_BUTTONS=false; //true;
	public static boolean MORE_BUTTONS1=false; //true;

	static File DIR;
	public static String DEFAULT_DIRECTORY=null;

	public ImagePlus  imp_sel=null;
//	public ImagePlus imp_distortions=null;
	//   public static boolean [] bPattern=null; // pattern bitmap
	public static double [][][][] PSF_KERNEL_MAP=null; // remove?
	public static String [] COMPONENT_COLOR_NAMES={"green1","red","blue","green2", "greens (diagonal)", "greens (checker)"};
	public static String [] STACK_COLOR_NAMES={"red","green","blue"};
/**
 * DIST_ARRAY:
 *  [v][u][0][0] - pixel x of the grid node u,v
 *  [v][u][0][1] - pixel y of the grid node u,v
 *  [v][u][1][0] - Wave vector 1 x component of the grid node u,v
 *  [v][u][1][1] - Wave vector 1 y component of the grid node u,v
 *  [v][u][2][0] - Wave vector 2 x component of the grid node u,v
 *  [v][u][2][1] - Wave vector 2 y component of the grid node u,v
 *  
*/
//	public static Distortions LENS_DISTORTIONS;
// Moved to MatchSimulatedPattern	
//	public static double [][][][] DIST_ARRAY=null;
//	public static Rectangle DIST_SELECTION=null;
//	public int [] UV_INDEX=null; // array containing index of the pattern UV (scanline order, U first), or -1 for the areas with no pattern
// End of Moved to MatchSimulatedPattern	
	SyncCommand SYNC_COMMAND=new SyncCommand();
	public float [][] SIM_ARRAY=null; // first index - 0 - main, 1 - shifted by 0.5 pixel diagonally (to extract checker greens)

	public static CalibrationHardwareInterface.LaserPointers LASERS=new CalibrationHardwareInterface.LaserPointers(LASER_POINTERS);
	public static CalibrationHardwareInterface.CamerasInterface CAMERAS=new CalibrationHardwareInterface.CamerasInterface(26,LASERS);
	public static CalibrationHardwareInterface.FocusingMotors MOTORS=new CalibrationHardwareInterface.FocusingMotors();
	
	public static Distortions.DistortionProcessConfiguration DISTORTION_PROCESS_CONFIGURATION=new Distortions.DistortionProcessConfiguration();
	
	public static LensAdjustment.FocusMeasurementParameters FOCUS_MEASUREMENT_PARAMETERS= new LensAdjustment.FocusMeasurementParameters(MOTORS.curpos);
	public static CalibrationHardwareInterface.GoniometerMotors GONIOMETER_MOTORS= new CalibrationHardwareInterface.GoniometerMotors();
	public static FocusingField FOCUSING_FIELD=null;
	//GoniometerParameters
	public static Goniometer.GoniometerParameters GONIOMETER_PARAMETERS= new Goniometer.GoniometerParameters(GONIOMETER_MOTORS);
	
	

	public static CalibrationHardwareInterface.UVLEDandLasers UV_LED_LASERS=new CalibrationHardwareInterface.UVLEDandLasers(FOCUS_MEASUREMENT_PARAMETERS);
	
	public static LensAdjustment LENS_ADJUSTMENT = new LensAdjustment();
	
	public static EyesisAberrations.AberrationParameters ABERRATIONS_PARAMETERS=new EyesisAberrations.AberrationParameters();
	public static EyesisAberrations EYESIS_ABERRATIONS; // need Distortions to be set up
	
	public static Goniometer GONIOMETER=null; 
//	new CalibrationHardwareInterface.LaserPointers();
	public class SyncCommand{
	    public boolean isRunning=      false;
	    public AtomicInteger stopRequested=  new AtomicInteger(0); // 0 - not requested, 1 - ASAP, 2 - gracefully
	    public String  buttonLabel="";
	}
	public Aberration_Calibration() {
		super("Aberration_Calibration");
		if (IJ.versionLessThan("1.43q")) return;
		if (instance!=null) {
			instance.toFront();
			return;
		}
		EYESIS_ABERRATIONS=new EyesisAberrations(SYNC_COMMAND.stopRequested,
				ABERRATIONS_PARAMETERS); // need Distortions to be set up
//		System.out.println("SwingUtilities.isEventDispatchThread()="+SwingUtilities.isEventDispatchThread());
		instance = this;
		addKeyListener(IJ.getInstance());
//		setLayout(new GridLayout(ADVANCED_MODE?8:5, 1));
//		setLayout(new GridLayout(ADVANCED_MODE?9:6, 1));
		setLayout(new GridLayout(ADVANCED_MODE?20:20, 1));
		Color color_configure=     new Color(200, 200,160);
		Color color_process=       new Color(180, 180, 240);
		Color color_conf_process=  new Color(180, 240, 240);
		Color color_restore=       new Color(180, 240, 180);
		Color color_stop=          new Color(255, 160, 160);
		Color color_goniometer=    new Color(160, 250, 160);
		Color color_lenses=        new Color(250, 250, 160);
		Color color_bundle=        new Color(250, 160, 250);
		Color color_report=        new Color(160, 250, 250);
		Color color_debug=         new Color(250, 250,  80);
		Color color_aberration=    new Color(60,  250, 250);


		panelRun = new Panel();
		panelRun.setLayout(new GridLayout(1, 0, 5, 5));
		addButton("Process Calibration Files",panelRun);
		addButton("Save",panelRun);
		addButton("Save Selected",panelRun);
		addButton("Restore",panelRun,color_restore);
		addButton("Restore no autoload",panelRun);
		addButton("Restore SFE Latest",panelRun,color_restore);
		addButton("List SFE",panelRun,color_report);
		addButton("Stop",panelRun,color_stop);
		addButton("Abort",panelRun,color_stop);
		
		add(panelRun);
		
		panelDirs= new Panel();
		panelDirs.setLayout(new GridLayout(1, 0, 5, 5));
		addButton("Select source directory",      panelDirs);
		addButton("Select intermediate directory",panelDirs);
		addButton("Select results directory",     panelDirs);
		add(panelDirs);
		
		

		panelConf1 = new Panel();
		panelConf1.setLayout(new GridLayout(1, 0, 5, 5));
		addButton("Configure Globals",panelConf1);
		addButton("Conf. Components",panelConf1);
		addButton("Conf. Multifile",panelConf1,color_configure);
		addButton("Conf. Simulation",panelConf1);
		addButton("Conf. Pattern Detection",panelConf1);
		add(panelConf1);

		panelConf2 = new Panel();
		panelConf2.setLayout(new GridLayout(1, 0, 5, 5));
		addButton("Conf. PSF",panelConf2);
		addButton("Conf. OTF Filter",panelConf2);
		addButton("Conf. Interpolation",panelConf2,color_configure);
		addButton("Conf. Inversion",panelConf2,color_configure);
		addButton("Conf. Results",panelConf2);
		add(panelConf2);
//		/panelAberrations

		if (ADVANCED_MODE){
			panel1 = new Panel();
			panel1.setLayout(new GridLayout(1, 0, 5, 5)); // rows, columns, vgap, hgap
			addButton("Map PSF",panel1);
			addButton("Scan and Map",panel1);
			addButton("Combine PSF files",panel1);
			add(panel1);

			panel2 = new Panel();
			panel2.setLayout(new GridLayout(1, 0, 5, 5));
			addButton("Interpolate Stack",panel2);
			addButton("Invert Stack",panel2);
			addButton("Gaussian Stack",panel2);
			add(panel2);

			panel3 = new Panel();
			panel3.setLayout(new GridLayout(1, 0, 5, 5));
			addButton("Configure",panel3);
			addButton("Configure Simul.",panel3);
			add(panel3);

		}

		panelDistortions = new Panel();
		panelDistortions.setLayout(new GridLayout(1, 0, 5, 5));
		if (ADVANCED_MODE) addButton("Correlation tests",panelDistortions);
		addButton("Configure Distortion",panelDistortions,color_configure);
		addButton("Distortion",panelDistortions);
		addButton("Grid Brightness",panelDistortions);
		addButton("Configure Pointers",panelDistortions,color_configure);
		addButton("List Pointers",panelDistortions,color_report);
		addButton("Simulate Full",panelDistortions);
		add(panelDistortions);
		
		
		panelFitDistortions = new Panel();
		panelFitDistortions.setLayout(new GridLayout(1, 0, 5, 5));
//		addButton("Configure Fit",panelFitDistortions);
		addButton("Select Grid Files",panelFitDistortions,color_configure);
		addButton("Edit Calibration",panelFitDistortions,color_configure);
		addButton("Save Calibration",panelFitDistortions,color_bundle);
		addButton("Restore Calibration",panelFitDistortions,color_restore);
		addButton("New Strategy",panelFitDistortions,color_bundle);
		addButton("Edit Strategy",panelFitDistortions,color_configure);
		addButton("Save Strategy",panelFitDistortions,color_bundle);
		addButton("Restore Strategy",panelFitDistortions,color_restore);
		addButton("Run LMA",panelFitDistortions,color_bundle);
		addButton("Bad Nodes",panelFitDistortions,color_bundle);
		addButton("Debug deriv",panelFitDistortions,color_debug);
		add(panelFitDistortions);
/*		
		panelFlat = new Panel();
		panelFlat.setLayout(new GridLayout(1, 0, 5, 5));
		addButton("Configure Flat Field",panelFlat);
		addButton("Remove FF Source Dirs",panelFlat);
		addButton("Add FF Source Dirs",panelFlat);
		addButton("Accumulate Flat Fields",panelFlat);
		addButton("Process Flat Field",panelFlat);
		add(panelFlat);
*/		
//panelLens		
		panelLens = new Panel();
		panelLens.setLayout(new GridLayout(1, 0, 5, 5));
		addButton("Configure Pattern Dimensions",panelLens,color_configure);
		addButton("Configure Distortion/Location",panelLens,color_configure);
		addButton("Configure Eyesis4pi",panelLens,color_configure);
		addButton("List Eyesis4pi",panelLens,color_report);
		addButton("Process Lens Distortion",panelLens);
		addButton("Configure Lasers",panelLens,color_configure);
		addButton("Manual laser pointers",panelLens,color_debug);
		addButton("Configure Cameras",panelLens,color_configure);
		addButton("Cameras Settings",panelLens,color_configure);
		addButton("Test Cameras",panelLens,color_debug);
		addButton("Test No Lasers",panelLens,color_debug);
		add(panelLens);
//panelProcessDistortions
		panelProcessDistortions = new Panel();
		panelProcessDistortions.setLayout(new GridLayout(1, 0, 5, 5));
		addButton("Configure Process Distortions",panelProcessDistortions,color_configure);
		addButton("Quick get&show",panelProcessDistortions,color_debug);
		addButton("Acquire",panelProcessDistortions);
		addButton("Calculate grids",panelProcessDistortions,color_goniometer);
		addButton("Convert X/Y slices to color",panelProcessDistortions); // move elsewhere
		addButton("Grid candidate",panelProcessDistortions); // just for testing
		add(panelProcessDistortions);
//panelCorrectGrid
		panelCorrectGrid = new Panel();
		panelCorrectGrid.setLayout(new GridLayout(1, 0, 5, 5));
		addButton("Calculate Sensor Masks",panelCorrectGrid,color_bundle);
		addButton("Save Sensor Masks",panelCorrectGrid,color_bundle);
		addButton("Restore Sensor Masks",panelCorrectGrid,color_restore);
		
		addButton("Reset Grid",panelCorrectGrid,color_bundle);
		addButton("Reset Margins",panelCorrectGrid,color_bundle);
		addButton("Restore Grid",panelCorrectGrid,color_restore);
		addButton("Reset Variations",panelCorrectGrid,color_restore);
		if (MORE_BUTTONS1) addButton("Correct Grid0",panelCorrectGrid,color_process);
		addButton("Correct Grid",panelCorrectGrid,color_process);
		addButton("Reset Variations",panelCorrectGrid,color_bundle);
		addButton("Grid Diffs",panelCorrectGrid,color_process);
		addButton("Save Grid",panelCorrectGrid,color_bundle);
		
		addButton("Reset Sensor",panelCorrectGrid,color_bundle);
		addButton("Restore Sensor",panelCorrectGrid,color_restore);
		if (MORE_BUTTONS1) addButton("Correct Sensor Old",panelCorrectGrid,color_process);
		addButton("Correct Sensor",panelCorrectGrid,color_process);
		addButton("Save Sensor",panelCorrectGrid,color_bundle);
//		addButton("TestIpl",panelCorrectGrid);
		
		addButton("List Calibration",panelCorrectGrid,color_report);
		addButton("Reproject",panelCorrectGrid,color_debug);
		
		add(panelCorrectGrid);
///panelFocusing
		
		panelFocusing1 = new Panel();
		panelFocusing1.setLayout(new GridLayout(1, 0, 5, 5));
		addButton("Reset Focusing",panelFocusing1);
if (MORE_BUTTONS){		
		addButton("Get Focusing Grid",panelFocusing1);
		addButton("Update Focusing Grid",panelFocusing1);
		addButton("Focusing PSF",panelFocusing1);
		addButton("Focusing New PSF",panelFocusing1);
		addButton("Focusing Acquire PSF",panelFocusing1);
}
if (MORE_BUTTONS) {		
		addButton("Probe around",panelFocusing1);
}		
//        addButton("Reset Histories",panelFocusing1);
		addButton("Manual Pre-focus",panelFocusing1,color_lenses);
		addButton("List Pre-focus",panelFocusing1);
		addButton("Manual Focus/Tilt",panelFocusing1,color_lenses);
		addButton("Fine Focus",panelFocusing1);
		addButton("Calibrate Distance",panelFocusing1);
		addButton("Show Grid",panelFocusing1);
		addButton("Lasers Toggle",panelFocusing1,color_debug);
		addButton("UV on", panelFocusing1);
		addButton("UV off",panelFocusing1,color_lenses);
		addButton("Camera Power Cycled",panelFocusing1,color_lenses);
		addButton("Acquire&Save",panelFocusing1);
		addButton("No-move measure",panelFocusing1,color_lenses);
		
		add(panelFocusing1);

		
		panelFocusing = new Panel();
		panelFocusing.setLayout(new GridLayout(1, 0, 5, 5));
		addButton("Configure Focusing",panelFocusing,color_configure);
		if (MORE_BUTTONS) {		
			addButton("Head Orientation",panelFocusing);
		}
		addButton("Lens Center",panelFocusing,color_process);
//		if (MORE_BUTTONS) {		
			addButton("Find Grid",panelFocusing,color_lenses);
//		}
        addButton("Select WOI",panelFocusing,color_lenses);
        addButton("Reset Histories",panelFocusing,color_lenses);
        addButton("Motors Home",panelFocusing,color_lenses);
		addButton("Auto Pre-focus",panelFocusing,color_process);
		addButton("Scan Calib",panelFocusing,color_process);
		
		addButton("Auto Focus/Tilt",panelFocusing,color_process);
//		addButton("List Pre-focus",panelFocusing);
		addButton("Focus Average",panelFocusing,color_report);
		addButton("Temp. Scan",panelFocusing,color_process);
		//
		addButton("List History",panelFocusing,color_report);
		addButton("Show PSF",panelFocusing,color_report);
		add(panelFocusing);
		
// panelCurvature		
		panelCurvature=new Panel();
		panelCurvature.setLayout(new GridLayout(1, 0, 5, 5));
		addButton("Scan Calib LMA", panelCurvature,color_process);
		addButton("Save History",   panelCurvature,color_debug);
		addButton("Restore History",panelCurvature,color_restore);
		addButton("History RMS",    panelCurvature,color_report);
		addButton("Modify LMA",     panelCurvature,color_configure);
		addButton("Load strategies", panelCurvature,color_restore);
		addButton("Organize strategies", panelCurvature,color_configure);
		addButton("Save strategies", panelCurvature,color_bundle);
		addButton("LMA History",    panelCurvature,color_process);
		addButton("List curv pars", panelCurvature,color_debug);
		addButton("List curv data", panelCurvature,color_debug);
		addButton("List qualB",     panelCurvature,color_report);
		addButton("List curv",      panelCurvature,color_report);
		addButton("Show curv corr", panelCurvature,color_report);
		addButton("Test measurement", panelCurvature,color_debug);
		addButton("Optimize qualB", panelCurvature,color_debug);
		addButton("Focus/Tilt LMA", panelCurvature,color_process);
		add(panelCurvature);
		
	//panelGoniometer
		
		panelGoniometer = new Panel();
		panelGoniometer.setLayout(new GridLayout(1, 0, 5, 5));
		addButton("Configure Goniometer",panelGoniometer,color_configure); 
		addButton("Goniometer Scan",panelGoniometer,color_conf_process);
		addButton("Filter Grids",panelGoniometer,color_bundle);
		addButton("Update Image Set",panelGoniometer);
		
		addButton("Remove Outlayers",panelGoniometer,color_bundle);
		addButton("Remove Sets",panelGoniometer,color_bundle);
		
		addButton("Update Sets Orientation",panelGoniometer);
		addButton("List Image Sets",panelGoniometer,color_report);
		addButton("Re-calibrate Grids",panelGoniometer,color_bundle);
		addButton("Re-calibrate Set",panelGoniometer,color_bundle);
		addButton("Get Orientation",panelGoniometer);
		addButton("Test Hinted Grid",panelGoniometer);
		addButton("Test Hinted Grid Cameras",panelGoniometer);
		addButton("Simulate Grid View",panelGoniometer);
		addButton("Show grid/hint",    panelGoniometer,color_debug);
//		addButton("Test Progress",    panelGoniometer);
		
		add(panelGoniometer);

		
		panelPixelMapping = new Panel();
		panelPixelMapping.setLayout(new GridLayout(1, 0, 5, 5));
		addButton("Load Pixel Mapping",panelPixelMapping); 
		addButton("List Mapping Parameters",panelPixelMapping,color_report);
		addButton("Test Direct Mapping",panelPixelMapping); 
		addButton("Test Equirectangular Mapping",panelPixelMapping); 
		addButton("Crop Equirectangular Mapping",panelPixelMapping);
		addButton("Generate & Save Equirectangular",panelPixelMapping);
		addButton("Load Equirectangular Maps",panelPixelMapping);
		addButton("Show Maps Overlap",panelPixelMapping);
		addButton("Test Lanczos",panelPixelMapping);
		addButton("Warp Image",panelPixelMapping);
		addButton("Warp Files",panelPixelMapping);
		addButton("Pattern Flat-Field",panelPixelMapping,color_conf_process);
		addButton("Remove Specular",panelPixelMapping,color_conf_process);
		
		addButton("Flat-Field",panelPixelMapping,color_conf_process);
		
		add(panelPixelMapping);
		
		panelAberrations = new Panel();
		panelAberrations.setLayout(new GridLayout(1, 0, 5, 5));
		addButton("Configure aberrations",panelAberrations,color_configure);
		addButton("Select Channels",panelAberrations,color_configure);
		addButton("Partial Kernels",panelAberrations,color_aberration);
		addButton("Combine Kernels",panelAberrations,color_aberration);
		addButton("Interpolate Kernels",panelAberrations,color_aberration);
		addButton("Invert Kernels",panelAberrations,color_aberration);
		addButton("Create Plane Map",panelAberrations);
		addButton("Test Plane Map",panelAberrations);
		add(panelAberrations);

		panelStereo= new Panel();
		panelStereo.setLayout(new GridLayout(1, 0, 5, 5));
		addButton("Create Intermaps",panelStereo);
		addButton("Test Intermaps",panelStereo);
		addButton("Show Sobel",panelStereo);
		addButton("Edge Thinning",panelStereo);
//		addButton("Segmenting",panelStereo);
		addButton("Distance from Edges",panelStereo);
		addButton("Edge Areas",panelStereo);
		addButton("Vacuum Edges",panelStereo);
		addButton("Test Edges",panelStereo);
		add(panelStereo);
		
		
		panelStereo1= new Panel();
		panelStereo1.setLayout(new GridLayout(1, 0, 5, 5));
		addButton("Intercam correlations",panelStereo1);
		addButton("Test Line",panelStereo1);
		addButton("Linear Features",panelStereo1);
		addButton("Intercam rectangular", panelStereo1);
		addButton("Thershold Disparity", panelStereo1);
		addButton("Disparity Tiles0",panelStereo1);
		addButton("Disparity Tiles",panelStereo1);
		addButton("Disparity Section",panelStereo1);
		addButton("Disparity Map",panelStereo1);
		addButton("Create Ambiguity",panelStereo1);
		addButton("Initial Resolve",panelStereo1);
		addButton("Ambiguity Resolve",panelStereo1);
		
		add(panelStereo1);
		
		pack();
		GUI.center(this);
		setVisible(true);
		JP4_INSTANCE=       new JP46_Reader_camera(false);
		SDFA_INSTANCE=      new showDoubleFloatArrays();
// main loop
		while (true){
			synchronized (this.SYNC_COMMAND) {
				try {
					this.SYNC_COMMAND.wait();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				this.SYNC_COMMAND.isRunning=true;
			}
			if (this.SYNC_COMMAND.stopRequested.get()==0){
				try {
					runMenuCommand(this.SYNC_COMMAND.buttonLabel);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
//					System.out.println(stack2string(e));
					IJ.showMessage("Exception",stack2string(e));
				}
			}
			this.SYNC_COMMAND.isRunning=false;
			this.SYNC_COMMAND.stopRequested.set(0);
		}
		
	}
	public static String stack2string(Exception e) {
		try {
			StringWriter sw = new StringWriter();
			PrintWriter pw = new PrintWriter(sw);
			e.printStackTrace(pw);
			return sw.toString();
		}
		catch(Exception e2) {
			return "bad stack2string";
		}
	}
	void addButton(String label, Panel panel,Color color) {
		Button b = new Button(label);
		b.setBackground(color);
		b.addActionListener(this);
		b.addKeyListener(IJ.getInstance());
		panel.add(b);
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
		String label = e.getActionCommand();
		if        (label.equals("Abort")) {
			this.SYNC_COMMAND.stopRequested.set(1);
			return;
		} else if (label.equals("Stop")) {
			this.SYNC_COMMAND.stopRequested.set(2);
			return;
		}
//		System.out.println("actionPerformed: SwingUtilities.isEventDispatchThread()="+SwingUtilities.isEventDispatchThread());
		// executed in EDT, so no need to synchronize additionally?
		if (this.SYNC_COMMAND.isRunning){
			GenericDialog gd=new GenericDialog("Busy");
			gd.addMessage("Command \""+this.SYNC_COMMAND.buttonLabel+"\" is running. Ask it to stop "+((this.SYNC_COMMAND.stopRequested.get()!=0)?"(again!) ":"")+"(when possible)?");
			gd.enableYesNoCancel("ASAP", "When convenient");
			gd.showDialog();
			if (gd.wasCanceled()) return;
			this.SYNC_COMMAND.stopRequested.set(gd.wasOKed()?1:2);
			return;
		}
		synchronized (this.SYNC_COMMAND) {
			this.SYNC_COMMAND.buttonLabel=label;
			this.SYNC_COMMAND.notify();
		}
		//		matchSimulatedPattern.FFT_SIZE=FFT_SIZE;
	}
	public void runMenuCommand(String label){
		for (COMPONENTS.referenceComponent=5;(COMPONENTS.referenceComponent>=0) && (!COMPONENTS.colorsToCorrect[COMPONENTS.referenceComponent]); COMPONENTS.referenceComponent--);
		Runtime runtime = Runtime.getRuntime();
	    runtime.gc();
		if (DEBUG_LEVEL>0) System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");

		if (label==null) return;
/* ======================================================================== */
		if       (label.equals("Configure Globals")) {
			showConfigureGlobalsDialog();
			return;
/* ======================================================================== */
		} else if (label.equals("Conf. Components")) {
			showcolorComponentsDialog(COMPONENTS);
			return;
/* ======================================================================== */
		} else if (label.equals("Conf. Multifile")) {
			showMultiFilePSFDialog(MULTIFILE_PSF);
			return;
/* ======================================================================== */
		} else if (label.equals("Conf. Simulation")) {
			showSimulParametersDialog(SIMUL);
			return;
/* ======================================================================== */
		} else if (label.equals("Conf. Pattern Detection")) {
			showPatternDetectParametersDialog(PATTERN_DETECT);
			return;
/* ======================================================================== */
		} else if (label.equals("Conf. PSF")) {
			showPSFParametersDialog(PSF_PARS);
			return;
/* ======================================================================== */
		} else if (label.equals("Conf. OTF Filter")) {
			showOTFFilterParametersDialog(OTF_FILTER);
			return;
/* ======================================================================== */
		} else if (label.equals("Conf. Interpolation")) {
			showInterpolateParametersDialog(INTERPOLATE);
			return;
/* ======================================================================== */
		} else if (label.equals("Conf. Inversion")) {
			showInverseParametersDialog(INVERSE);
			return;
/* ======================================================================== */
		} else if (label.equals("Conf. Results")) {
			showShowResultsDialog(SHOW_RESULTS);
			return;
/* ======================================================================== */
		} else if (label.equals("Select source directory")) {
	    	String fileName= selectSourceDirectory(PROCESS_PARAMETERS.sourceSuperDirectory);
	        if (fileName!=null) PROCESS_PARAMETERS.sourceSuperDirectory=fileName;
			return;
/* ======================================================================== */
		} else if (label.equals("Select intermediate directory")) {
	    	String fileName= selectPartialKernelsDirectory(PROCESS_PARAMETERS.partialKernelsSuperDirectory);
	        if (fileName!=null) PROCESS_PARAMETERS.partialKernelsSuperDirectory=fileName;
			return;
/* ======================================================================== */
		} else if (label.equals("Select results directory")) {
	    	String fileName= selectKernelsDirectory(PROCESS_PARAMETERS.kernelsDirectory);
	        if (fileName!=null) PROCESS_PARAMETERS.kernelsDirectory=fileName;
			return;
/* ======================================================================== */
		     
	    } else if (label.equals("Save")) {
	    	saveProperties(null,PROCESS_PARAMETERS.kernelsDirectory,PROCESS_PARAMETERS.useXML, PROPERTIES);
	    	return;
/* ======================================================================== */
	    } else if (label.equals("Save Selected")) {
	    	Properties selectedProperties=new Properties();
	    	selectedProperties.setProperty("selected", "true");
	    	saveProperties(null,PROCESS_PARAMETERS.kernelsDirectory, PROCESS_PARAMETERS.useXML, selectedProperties);
	    	return;
	    	
/* ======================================================================== */
	        
	    } else if (label.equals("Restore") || label.equals("Restore no autoload")) {
	    	boolean noAuto=label.equals("Restore no autoload");
	    	ABERRATIONS_PARAMETERS.autoRestore=false;
	    	loadProperties(null,PROCESS_PARAMETERS.kernelsDirectory,PROCESS_PARAMETERS.useXML, PROPERTIES);
	    	if (ABERRATIONS_PARAMETERS.autoRestore && !noAuto){
	    		if (DEBUG_LEVEL>0)System.out.println("Auto-loading configuration files");
	    		if (LENS_DISTORTIONS==null) {
	    			LENS_DISTORTIONS=new Distortions(LENS_DISTORTION_PARAMETERS,PATTERN_PARAMETERS,REFINE_PARAMETERS,this.SYNC_COMMAND.stopRequested);
	    		}
	    		LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
	    		boolean dcdUpdated=autoLoadFiles(
	    				ABERRATIONS_PARAMETERS,
	    				LENS_DISTORTIONS, // should be initialized, after update DISTORTION_CALIBRATION_DATA from this
	    				PATTERN_PARAMETERS,
	    				EYESIS_CAMERA_PARAMETERS, //Distortions.EyesisCameraParameters eyesisCameraParameters,
	    				UPDATE_STATUS,
	    				DEBUG_LEVEL
	    		);
	    		if (dcdUpdated) DISTORTION_CALIBRATION_DATA=LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData;
		    	if (ABERRATIONS_PARAMETERS.autoReCalibrate){
					if (LENS_DISTORTIONS.fittingStrategy==null) {
						IJ.showMessage("Can not recalibrate grids - LENS_DISTORTION.fittingStrategy is not set");
						return;
					}
		    		if (DEBUG_LEVEL>0)System.out.println("=== Re-calibrating grids ===");
		    		int numMatched=LENS_DISTORTIONS.applyHintedGrids(
		    				LASER_POINTERS, // MatchSimulatedPattern.LaserPointer laserPointer, // LaserPointer object that specifies actual laser poiners on the target
		    				DISTORTION_PROCESS_CONFIGURATION.removeOutOfGridPointers, // boolean removeOutOfGridPointers,
		    				5.0,                   //double  hintGridTolerance, // alllowed mismatch (fraction of period) or 0 - orientation only
		    				true, //boolean processAll, // if true - process all images, false - only disabeld
		    				false, //?
		    				false, //processBlind,
		    				-1,    //       imageNumber,
		    				true, // useSetData - use imageSets data if available (false - use camera data)
		    				THREADS_MAX,                 //int threadsMax,
		    				UPDATE_STATUS,               // boolean updateStatus,
		    				DISTORTION.loop_debug_level, // int mspDebugLevel,
		    				MASTER_DEBUG_LEVEL,          //int global_debug_level, // DEBUG_LEVEL
		    				MASTER_DEBUG_LEVEL           //int debug_level // debug level used inside loops
		    		);
		    		System.out.println("Number of matched images: "+numMatched);
		    		
					LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
					LENS_DISTORTIONS.updateStatus=UPDATE_STATUS;
					LENS_DISTORTIONS.fittingStrategy.debugLevel=DEBUG_LEVEL;
					LENS_DISTORTIONS.markBadNodces(
							-1,// series - all images
							DEBUG_LEVEL);

		    		
		    	}	    		
	    	}
	    	return;
		
/* ======================================================================== */
		} else if (label.equals("Process Calibration Files")) {
			if (!showProcessCalibrationFilesDialog(PROCESS_PARAMETERS)) return;
			boolean noMessageBoxes=false; //TODO: move to configuration 
			long 	  startTime=System.nanoTime();
			long      tmpTime;
			String resultPath;
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			int loop_debug_level=1;
			String [][][] filePaths=null;
			int dirNum, fileNum;
			ImageStack stack;
			ImagePlus imp_psf;
			File file;
			String fileName;
			if (PROCESS_PARAMETERS.processSourceImages) {
				if ((PROCESS_PARAMETERS.sourceSuperDirectory==null) || (PROCESS_PARAMETERS.sourceSuperDirectory.length()==0)) {
			    	fileName= selectSourceDirectory(PROCESS_PARAMETERS.sourceSuperDirectory);
			        if (fileName!=null) PROCESS_PARAMETERS.sourceSuperDirectory=fileName;
				}
				if ((PROCESS_PARAMETERS.partialKernelsSuperDirectory==null) || (PROCESS_PARAMETERS.partialKernelsSuperDirectory.length()==0)) {
			    	fileName= selectPartialKernelsDirectory(PROCESS_PARAMETERS.partialKernelsSuperDirectory);
			        if (fileName!=null) PROCESS_PARAMETERS.partialKernelsSuperDirectory=fileName;
				}
				filePaths=prepareCalibrationFilesList(PROCESS_PARAMETERS);
				if (filePaths==null) {
					IJ.showMessage("Error","No files to process\nProcess canceled");
					return;
				}
				startTime=System.nanoTime(); // restart timer after possible interactive dialogs
				for (dirNum=0;dirNum<filePaths.length; dirNum++) {
					if (filePaths[dirNum]!=null){
						if (DEBUG_LEVEL>1) System.out.println("======= directory number="+dirNum+", number of files="+filePaths[dirNum].length+" ========");
						for (fileNum=0;fileNum<filePaths[dirNum].length; fileNum++) {
							if (DEBUG_LEVEL>1) System.out.println(filePaths[dirNum][fileNum][0] +" ==> "+filePaths[dirNum][fileNum][1]);
							imp_sel=JP4_INSTANCE.open(
									"", // path,
									filePaths[dirNum][fileNum][0],
									"",  //arg - not used in JP46 reader
									true, // un-apply camera color gains
									imp_sel); // reuse the same image window
// Remove for old method?
							matchSimulatedPattern= new MatchSimulatedPattern(DISTORTION.FFTSize);
//							   matchSimulatedPattern.invalidateFlatFieldForGrid(); //It is already reset, no need to do it again
//							   matchSimulatedPattern.invalidateFocusMask();
							
							matchSimulatedPattern.calculateDistortions(
									DISTORTION, //
									PATTERN_DETECT,
									SIMUL,
									COMPONENTS.equalizeGreens,
									imp_sel,
									null, // LaserPointer laserPointer, // LaserPointer object or null
									true, // don't care -removeOutOfGridPointers
									null, //   double [][][] hintGrid, // predicted grid array (or null)
									0,    //   double  hintGridTolerance, // alllowed mismatch (fraction of period) or 0 - orientation only

									THREADS_MAX,
									UPDATE_STATUS,
									DEBUG_LEVEL,
									DISTORTION.loop_debug_level, // debug level
									noMessageBoxes);
							
/*
  							SIM_ARRAY=	simulateGridAll (
 									imp_sel.getWidth(),
									imp_sel.getHeight(),
									matchSimulatedPattern.getDArray(),
									2, // gridFrac, // number of grid steps per pattern full period
									SIMUL,
									THREADS_MAX,
									UPDATE_STATUS,
									DISTORTION.loop_debug_level); // debug level
*/
							SIM_ARRAY=	(new SimulationPattern(SIMUL)).simulateGridAll (
									imp_sel.getWidth(),
									imp_sel.getHeight(),
									matchSimulatedPattern,
									2, // gridFrac, // number of grid steps per pattern full period
									SIMUL,
									THREADS_MAX,
									UPDATE_STATUS,
									DEBUG_LEVEL,
									DISTORTION.loop_debug_level); // debug level
							
							createPSFMap(
									matchSimulatedPattern,
									matchSimulatedPattern.applyFlatField (imp_sel), // if grid is flat-field calibrated, apply it
//									imp_sel, // linearized Bayer mosaic image form the camera, GR/BG
									null,     //  int [][][] sampleList, // optional (or null) 2-d array: list of coordinate pairs (2d - to match existent  PSF_KERNEL_MAP structure)  
									MULTIFILE_PSF.overexposedMaxFraction,
									SIMUL, //simulation parameters
									MAP_FFT_SIZE, // scanImageForPatterns:FFT size
									PATTERN_DETECT,
									FFT_OVERLAP, // scanImageForPatterns:high-pass gaussian filter sigma when correlating power spectrum
									FFT_SIZE, // maximal distance between maximum on spectrum and predicted maximum on autocorrelation of gamma(|spectrum|)
									COMPONENTS,
									PSF_SUBPIXEL, // maximal iterations when looking for local maximum
									OTF_FILTER,
									PSF_PARS, // step of the new map (should be multiple of map step)
									PSF_PARS.minDefinedArea,
									INVERSE.dSize, // size of square used in the new map (should be multiple of map step)
									THREADS_MAX,
									UPDATE_STATUS,
									loop_debug_level);// debug level used inside loops
							DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
							stack=mergeKernelsToStack(PSF_KERNEL_MAP);
							if (stack!=null) {
								imp_psf = new ImagePlus(filePaths[dirNum][fileNum][1], stack);
								//							if (DEBUG_LEVEL>1) imp_psf.show();
								if (DEBUG_LEVEL>1) System.out.println("Saving result to"+filePaths[dirNum][fileNum][1]+ " at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
								FileSaver fs=new FileSaver(imp_psf);
								fs.saveAsTiffStack(filePaths[dirNum][fileNum][1]);
							} else {
								System.out.println("File "+filePaths[dirNum][fileNum][1]+ " has no useful PSF kernels - at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));    			
							}
						}
					}
				}
			}
			
			if (PROCESS_PARAMETERS.combinePSFfiles) {
				tmpTime=System.nanoTime();
				if ((PROCESS_PARAMETERS.partialKernelsSuperDirectory==null) || (PROCESS_PARAMETERS.partialKernelsSuperDirectory.length()==0)) {
			    	fileName= selectPartialKernelsDirectory(PROCESS_PARAMETERS.partialKernelsSuperDirectory);
			        if (fileName!=null) PROCESS_PARAMETERS.partialKernelsSuperDirectory=fileName;
				}
				filePaths=preparePartialKernelsFilesList(PROCESS_PARAMETERS);
				startTime-=(System.nanoTime()-tmpTime); // do not count time used for selection of files
				if (filePaths==null) {
					IJ.showMessage("Error","No partila kernel files to process, finished in "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+ " seconds");
					return;
				}
				String [] filenames;
				for (dirNum=0;dirNum<filePaths.length; dirNum++) {
					if ((filePaths[dirNum]!=null) && (filePaths[dirNum].length!=0)){
						resultPath=filePaths[dirNum][0][1];
						filenames=new String[filePaths[dirNum].length];
						for (fileNum=0;fileNum<filenames.length;fileNum++) filenames[fileNum]=filePaths[dirNum][fileNum][0];
						if (!combinePSFKernels (
								INTERPOLATE,
								MULTIFILE_PSF,
								filenames,
								resultPath,
								SDFA_INSTANCE,
								imp_sel,         // re-use global
								true,            // saveResult,
								false,           // showResult,
								UPDATE_STATUS, 
								DEBUG_LEVEL)) continue; // return; // no overlap, bad result kernel
					}

				}
			}
			if (PROCESS_PARAMETERS.interpolatePSFkernel) {
				tmpTime=System.nanoTime();
				if ((PROCESS_PARAMETERS.partialKernelsSuperDirectory==null) || (PROCESS_PARAMETERS.partialKernelsSuperDirectory.length()==0)) {
			    	fileName= selectPartialKernelsDirectory(PROCESS_PARAMETERS.partialKernelsSuperDirectory);
			        if (fileName!=null) PROCESS_PARAMETERS.partialKernelsSuperDirectory=fileName;
				}
				IJ.showMessage("Notice","partialKernelsSuperDirectory="+PROCESS_PARAMETERS.partialKernelsSuperDirectory);
				startTime-=(System.nanoTime()-tmpTime); // do not count time used for selection of files
				filePaths=prepareInterpolateKernelsList(PROCESS_PARAMETERS); 
				for (fileNum=0;fileNum<filePaths.length;fileNum++) if ((filePaths[fileNum]!=null) && (filePaths[fileNum].length>0)) {
					file=new File(filePaths[fileNum][0][0]);
					if (!file.exists()) {
						if (DEBUG_LEVEL>1) System.out.println("Raw PSF kernel stack file "+filePaths[fileNum][0][0]+" does not exist");
						continue;
					}
					imp_psf=new ImagePlus(filePaths[fileNum][0][0]);
					if (imp_psf.getStackSize()<3) {
						System.out.println("Need a 3-layer stack with raw PSF kernels");
						continue;
					}
					stack= interpolateKernelStack(imp_psf.getStack(), // Image stack, each slice consists of square kernels of one channel
							INTERPOLATE,
							UPDATE_STATUS); // update status info

					imp_psf = new ImagePlus(filePaths[fileNum][0][1], stack);
					if (DEBUG_LEVEL>2) {
					  imp_psf.getProcessor().resetMinAndMax(); // imp_psf will be reused
					  imp_psf.show();
					}  
					if (DEBUG_LEVEL>1) System.out.println("Saving interpolation result to"+filePaths[fileNum][0][1]+ " at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
					FileSaver fs=new FileSaver(imp_psf);
					fs.saveAsTiffStack(filePaths[fileNum][0][1]);
				}
			}	
			if (PROCESS_PARAMETERS.invertKernels) {
				tmpTime=System.nanoTime();
				if ((PROCESS_PARAMETERS.partialKernelsSuperDirectory==null) || (PROCESS_PARAMETERS.partialKernelsSuperDirectory.length()==0)) {
			    	fileName= selectPartialKernelsDirectory(PROCESS_PARAMETERS.partialKernelsSuperDirectory);
			        if (fileName!=null) PROCESS_PARAMETERS.partialKernelsSuperDirectory=fileName;
				}
				if ((PROCESS_PARAMETERS.kernelsDirectory==null) || (PROCESS_PARAMETERS.kernelsDirectory.length()==0)) {
			    	fileName= selectKernelsDirectory(PROCESS_PARAMETERS.kernelsDirectory);
			        if (fileName!=null) PROCESS_PARAMETERS.kernelsDirectory=fileName;
				}
				startTime-=(System.nanoTime()-tmpTime); // do not count time used fro selection of files
				filePaths=prepareInvertGaussianKernelsList(PROCESS_PARAMETERS, PROCESS_PARAMETERS.rpsfPrefix ); 

				for (fileNum=0;fileNum<filePaths.length;fileNum++) if ((filePaths[fileNum]!=null) && (filePaths[fileNum].length>0)) {
					file=new File(filePaths[fileNum][0][0]);
					if (!file.exists()) {
						if (DEBUG_LEVEL>0) System.out.println("Interpolated PSF kernel stack file "+filePaths[fileNum][0][0]+" does not exist");
						continue;
					}
					imp_psf=new ImagePlus(filePaths[fileNum][0][0]);
					if (imp_psf.getStackSize()<3) {
						System.out.println("Need a 3-layer stack with interpolated PSF kernels");
						continue;
					}
					stack= reversePSFKernelStack(imp_psf.getStack(), //  stack of 3 32-bit (float) images, made of square kernels
							INVERSE,
							THREADS_MAX,   // size (side of square) of reverse PSF kernel
							UPDATE_STATUS);         // update status info
					imp_psf = new ImagePlus(filePaths[fileNum][0][1], stack);
					if (DEBUG_LEVEL>2) {
					  imp_psf.getProcessor().resetMinAndMax(); // imp_psf will be reused
					  imp_psf.show();
					}  
					if (DEBUG_LEVEL>0) System.out.println("Saving PSF inversion result to"+filePaths[fileNum][0][1]+ " at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
					FileSaver fs=new FileSaver(imp_psf);
					fs.saveAsTiffStack(filePaths[fileNum][0][1]);
				}
			}
			if (PROCESS_PARAMETERS.gaussianKernels) {
				tmpTime=System.nanoTime();
				if ((PROCESS_PARAMETERS.partialKernelsSuperDirectory==null) || (PROCESS_PARAMETERS.partialKernelsSuperDirectory.length()==0)) {
			    	fileName= selectPartialKernelsDirectory(PROCESS_PARAMETERS.partialKernelsSuperDirectory);
			        if (fileName!=null) PROCESS_PARAMETERS.partialKernelsSuperDirectory=fileName;
				}
				if ((PROCESS_PARAMETERS.kernelsDirectory==null) || (PROCESS_PARAMETERS.kernelsDirectory.length()==0)) {
			    	fileName= selectKernelsDirectory(PROCESS_PARAMETERS.kernelsDirectory);
			        if (fileName!=null) PROCESS_PARAMETERS.kernelsDirectory=fileName;
				}
				startTime-=(System.nanoTime()-tmpTime); // do not count time used fro selection of files
				filePaths=prepareInvertGaussianKernelsList(PROCESS_PARAMETERS, PROCESS_PARAMETERS.gaussianPrefix); 
				for (fileNum=0;fileNum<filePaths.length;fileNum++) if ((filePaths[fileNum]!=null) && (filePaths[fileNum].length>0)) {
					file=new File(filePaths[fileNum][0][0]);
					if (!file.exists()) {
						if (DEBUG_LEVEL>0) System.out.println("Interpolated PSF kernel stack file "+filePaths[fileNum][0][0]+" does not exist");
						continue;
					}
					imp_psf=new ImagePlus(filePaths[fileNum][0][0]);
					if (imp_psf.getStackSize()<3) {
						System.out.println("Need a 3-layer stack with interpolated PSF kernels");
						continue;
					}
					stack= generateGaussianStackFromDirect(imp_psf.getStack(), // stack of 3 32-bit (float) images, made of square kernels
							INVERSE,
							UPDATE_STATUS);  // update status info
					imp_psf = new ImagePlus(filePaths[fileNum][0][1], stack);
					if (DEBUG_LEVEL>2) {
					  imp_psf.getProcessor().resetMinAndMax(); // imp_psf will be reused
					  imp_psf.show();
					}  
					if (DEBUG_LEVEL>0) System.out.println("Saving Gaussian kernels result to"+filePaths[fileNum][0][1]+ " at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
					FileSaver fs=new FileSaver(imp_psf);
					fs.saveAsTiffStack(filePaths[fileNum][0][1]);
				}
			}
			System.out.println("Processing done in "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+ " seconds");
			if (PROCESS_PARAMETERS.saveSettings) saveProperties(PROCESS_PARAMETERS.kernelsDirectory+Prefs.getFileSeparator()+"calibration_settings",null,PROCESS_PARAMETERS.useXML, PROPERTIES);
			return;
/* ======================================================================== */

/* ======================================================================== */
		} else if (label.equals("Configure")) {
			if (!showConfigureDialog(OTF_FILTER, PATTERN_DETECT,COMPONENTS,INVERSE,MULTIFILE_PSF)) return;
			return;
/* ======================================================================== */
		} else if (label.equals("Configure Simul.")) {
			if (!showSimulDialog(SIMUL)) return;
			return;

/* ======================================================================== */

		} else if (label.equals("Map PSF")) {
			if (!showPSFDialog(PSF_PARS,INVERSE,OTF_FILTER,SIMUL)) return;
			int loop_debug_level=DISTORTION.loop_debug_level;
//			int loop_debug_level=1;
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			imp_sel = WindowManager.getCurrentImage();
			if (imp_sel==null){
				IJ.showMessage("Error","There are no images open\nProcess canceled");
				return;
			}
			// Updates global PSF_KERNEL_MAP
//			update selected roi to include at least one PSF kernel
			int psf_overlap_sensor=FFT_OVERLAP*PSF_SUBPIXEL/2;
			int psf_fft_sensor=    FFT_SIZE*   PSF_SUBPIXEL/2;
			Roi roi= imp_sel.getRoi();
			Rectangle selection;
			if (roi!=null){
				selection=roi.getBounds();
				int selectionCenterX=selection.x+selection.width/2;
				int selectionCenterY=selection.y+selection.height/2;
				// round selection width, height
				selection.width=psf_overlap_sensor*((selection.width+psf_overlap_sensor/2)/psf_overlap_sensor);
				selection.height=psf_overlap_sensor*((selection.height+psf_overlap_sensor/2)/psf_overlap_sensor);
				// enforce minimal size
				if (selection.width <psf_fft_sensor) selection.width=psf_fft_sensor;
				if (selection.height<psf_fft_sensor) selection.height=psf_fft_sensor;
				// recalculate center
				selection.x=selectionCenterX-selection.width/2;
				selection.y=selectionCenterY-selection.width/2;
				if (selection.x<0) selection.x=0;
				if (selection.y<0) selection.y=0;
				// add half step for later rounding by truncating;
				selection.x+=psf_overlap_sensor/2;
				selection.y+=psf_overlap_sensor/2;
				//limit by right and bottom edges 
				if ((selection.width+selection.x)>imp_sel.getWidth()) selection.x=imp_sel.getWidth()-selection.width;
				if ((selection.height+selection.y)>imp_sel.getHeight()) selection.y=imp_sel.getHeight()-selection.height;
				// round to the grid (already added psf_overlap_sensor/2, so just truncate)
				selection.x=psf_overlap_sensor*(selection.x/psf_overlap_sensor);
				selection.y=psf_overlap_sensor*(selection.y/psf_overlap_sensor);
				
				imp_sel.setRoi(selection);
			}
/*
  			SIM_ARRAY=	simulateGridAll (
 					imp_sel.getWidth(),
					imp_sel.getHeight(),
					matchSimulatedPattern.getDArray(),
					2, // gridFrac, // number of grid steps per pattern full period
					SIMUL,
					THREADS_MAX,
					UPDATE_STATUS,
					DISTORTION.loop_debug_level); // debug level
*/			
			SIM_ARRAY=	(new SimulationPattern(SIMUL)).simulateGridAll (
					imp_sel.getWidth(),
					imp_sel.getHeight(),
					matchSimulatedPattern,
					2, // gridFrac, // number of grid steps per pattern full period
					SIMUL,
					THREADS_MAX,
					UPDATE_STATUS,
					DEBUG_LEVEL,
					DISTORTION.loop_debug_level); // debug level
			
			//(new showDoubleFloatArrays()).showArrays(kernels, "***kernels-"+nTX+"-"+nTY);
			createPSFMap(
					matchSimulatedPattern,
					matchSimulatedPattern.applyFlatField (imp_sel), // if grid is flat-field calibrated, apply it
//					imp_sel, // linearized Bayer mosaic image form the camera, GR/BG
					
					null,     //  int [][][] sampleList, // optional (or null) 2-d array: list of coordinate pairs (2d - to match existent  PSF_KERNEL_MAP structure)  
					MULTIFILE_PSF.overexposedMaxFraction,
					SIMUL, // simulation parameters
					MAP_FFT_SIZE, // scanImageForPatterns:FFT size
					PATTERN_DETECT, // pattern detection parameters
					FFT_OVERLAP,
					FFT_SIZE,  
					COMPONENTS,
					PSF_SUBPIXEL, // maximal iterations when looking for local maximum
					OTF_FILTER,
					PSF_PARS, // step of the new map (should be multiple of map step)
					PSF_PARS.minDefinedArea,
					INVERSE.dSize, // size of square used in the new map (should be multiple of map step)
					THREADS_MAX,
					UPDATE_STATUS,
					loop_debug_level);// debug level used inside loops
//if (DEBUG_LEVEL>0) return;			
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			ImageStack mergedStack=mergeKernelsToStack(PSF_KERNEL_MAP);
			if (mergedStack==null) {
				IJ.showMessage("Error","No useful PSF kernels to show");
				return;
			}
			ImagePlus imp_psf=SDFA_INSTANCE.showImageStack(mergedStack, imp_sel.getTitle()+"-PSF_KERNEL");
			if (PSF_SAVE_FILE) {
				String path=imp_sel.getOriginalFileInfo().directory+"PSF-"+imp_sel.getTitle();
				if (DEBUG_LEVEL>1) {
					System.out.println("Saving result to "+path);
				}
				IJ.saveAs(imp_psf,"tif",path);
			}
			return;
/* ======================================================================== */
			// calculate a single PSF      
/* ======================================================================== */
		} else if (label.equals("Interpolate Stack")) {
			//      int loop_debug_level=1;
			if (!showInterpolateKernelsDialog(INTERPOLATE)) return;
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
					INTERPOLATE,
					UPDATE_STATUS); // update status info

			ImagePlus imp_interpolated_stack = new ImagePlus(imp_kernels.getTitle()+"-"+INTERPOLATE.step+ "X-interpolated", interpolatedStack);
			imp_interpolated_stack.getProcessor().resetMinAndMax();
			imp_interpolated_stack.show();

			String	result_path=imp_kernels.getOriginalFileInfo().directory+imp_interpolated_stack.getTitle();
			if (PSF_SAVE_FILE && (result_path!=null)) {
				if (DEBUG_LEVEL>1) {
					System.out.println("Saving result to "+result_path);
				}
				IJ.saveAs(imp_interpolated_stack,"tif",result_path);
			}
			return;

/* ======================================================================== */
		} else if (label.equals("Scan and Map")) {
			boolean noMessageBoxes=false;
			int loop_debug_level=1;
			if (!showPSFDialog(PSF_PARS,INVERSE,OTF_FILTER,SIMUL)) return;
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			JFileChooser fc= new JFileChooser();
			fc.setMultiSelectionEnabled(true);

			if (DIR==null) { // global
				String sdir = OpenDialog.getDefaultDirectory();
				if (sdir!=null)
					DIR = new File(sdir);
			}
			if (DIR!=null)
				fc.setCurrentDirectory(DIR);
			int returnVal = fc.showOpenDialog(IJ.getInstance());
			if (returnVal!=JFileChooser.APPROVE_OPTION)
				return;
			File[] files = fc.getSelectedFiles();
			if (files.length==0) { // getSelectedFiles does not work on some JVMs
				files = new File[1];
				files[0] = fc.getSelectedFile();
			}
			long startTime = System.nanoTime();

			String path = fc.getCurrentDirectory().getPath()+Prefs.getFileSeparator();
			DIR = fc.getCurrentDirectory();
			System.out.println("path= "+path+", files:");
			String [] filenames=new String[files.length]; // global
			int nFile;
			ImageStack mergedStack;
			ImagePlus imp_psf;
			String outPath;
			for (nFile=0;nFile<files.length;nFile++) {
				filenames[nFile]= files[nFile].getName();
				System.out.println(nFile+": "+filenames[nFile]);
				imp_sel=JP4_INSTANCE.open(
						"", // path,
						path+filenames[nFile],
						"",  //arg - not used in JP46 reader
						true, // un-apply camera color gains
						imp_sel); // reuse the same image window
				matchSimulatedPattern= new MatchSimulatedPattern(DISTORTION.FFTSize);
//				   matchSimulatedPattern.invalidateFlatFieldForGrid(); //It is already reset, no need to do it again
//				   matchSimulatedPattern.invalidateFocusMask();
				matchSimulatedPattern.calculateDistortions(
						DISTORTION, //
						PATTERN_DETECT,
						SIMUL,
						COMPONENTS.equalizeGreens,
						imp_sel,
						null, // LaserPointer laserPointer, // LaserPointer object or null
						true, // don't care -removeOutOfGridPointers
						null, //   double [][][] hintGrid, // predicted grid array (or null)
						0,    //   double  hintGridTolerance, // alllowed mismatch (fraction of period) or 0 - orientation only
						THREADS_MAX,
						UPDATE_STATUS,
						DEBUG_LEVEL,
						DISTORTION.loop_debug_level, // debug level
						noMessageBoxes);
				/*
				SIM_ARRAY=	simulateGridAll (
						imp_sel.getWidth(),
						imp_sel.getHeight(),
						matchSimulatedPattern.getDArray(),
						2, // gridFrac, // number of grid steps per pattern full period
						SIMUL,
						THREADS_MAX,
						UPDATE_STATUS,
						DISTORTION.loop_debug_level); // debug level
				*/
				SIM_ARRAY=	(new SimulationPattern(SIMUL)).simulateGridAll (
						imp_sel.getWidth(),
						imp_sel.getHeight(),
						matchSimulatedPattern,
						2, // gridFrac, // number of grid steps per pattern full period
						SIMUL,
						THREADS_MAX,
						UPDATE_STATUS,
						DEBUG_LEVEL,
						DISTORTION.loop_debug_level); // debug level
				
				createPSFMap(
						matchSimulatedPattern,
						matchSimulatedPattern.applyFlatField (imp_sel), // if grid is flat-field calibrated, apply it
//						imp_sel, // linearized Bayer mosaic image form the camera, GR/BG
						null,     //  int [][][] sampleList, // optional (or null) 2-d array: list of coordinate pairs (2d - to match existent  PSF_KERNEL_MAP structure)  
						MULTIFILE_PSF.overexposedMaxFraction,
						SIMUL, //simulation parameters
						MAP_FFT_SIZE, // scanImageForPatterns:FFT size
						PATTERN_DETECT,
						FFT_OVERLAP, // scanImageForPatterns:high-pass gaussian filter sigma when correlating power spectrum
						FFT_SIZE, // maximal distance between maximum on spectrum and predicted maximum on autocorrelation of gamma(|spectrum|)
						COMPONENTS,
						PSF_SUBPIXEL, // maximal iterations when looking for local maximum
						OTF_FILTER,
						PSF_PARS, // step of the new map (should be multiple of map step)
						PSF_PARS.minDefinedArea,
						INVERSE.dSize, // size of square used in the new map (should be multiple of map step)
						THREADS_MAX,
						UPDATE_STATUS,
						loop_debug_level);// debug level used inside loops
				DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
				//				outPath=path+"PSF-"+filenames[nFile];
				outPath=updateExtension(path+"PSF-"+filenames[nFile],".tif");
				mergedStack=mergeKernelsToStack(PSF_KERNEL_MAP);
				if (mergedStack!=null) {
					imp_psf = new ImagePlus(imp_sel.getTitle()+"-PSF_KERNEL", mergedStack);
					if (DEBUG_LEVEL>1) System.out.println("Saving result to"+outPath+ " at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
					IJ.saveAs(imp_psf,"tif",outPath);
				} else {
					System.out.println("File "+outPath+ " has no useful PSF kernels - at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));    			
				}
			}

			return;

/* ======================================================================== */
		} else if (label.equals("Combine PSF files")) {
			if (!showValidateKernelsDialog(INTERPOLATE, MULTIFILE_PSF)) return;
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			JFileChooser fc= new JFileChooser();
			fc.setMultiSelectionEnabled(true);

			if (DIR==null) { // global
				String sdir = OpenDialog.getDefaultDirectory();
				if (sdir!=null)
					DIR = new File(sdir);
			}
			if (DIR!=null)
				fc.setCurrentDirectory(DIR);
			int returnVal = fc.showOpenDialog(IJ.getInstance());
			if (returnVal!=JFileChooser.APPROVE_OPTION)
				return;
			File[] files = fc.getSelectedFiles();
			if (files.length==0) { // getSelectedFiles does not work on some JVMs
				files = new File[1];
				files[0] = fc.getSelectedFile();
			}
			// 	  long startTime = System.nanoTime();

			String path = fc.getCurrentDirectory().getPath()+Prefs.getFileSeparator();
			DIR = fc.getCurrentDirectory();
			if (DEBUG_LEVEL>1) System.out.println("path= "+path+", files:");
			String [] filenames=new String[files.length]; // global
			double [][][][][] kernelsElllipsePars = new double[files.length][][][][];
			int nFile;
			Opener opener=new Opener();
			for (nFile=0;nFile<files.length;nFile++) {
				filenames[nFile]= files[nFile].getName();
				if (UPDATE_STATUS) IJ.showStatus("Scanning file "+(nFile+1)+" (of "+(files.length)+"): "+filenames[nFile]);
				if (DEBUG_LEVEL>1) System.out.println((nFile+1)+": "+filenames[nFile]);
				imp_sel=opener.openImage(path, filenames[nFile]);  // or (path+filenames[nFile])
				kernelsElllipsePars[nFile]= kernelStackToEllipseCoefficients(
						imp_sel.getStack(), // Image stack, each slice consists of square kernels of one channel
						INTERPOLATE.size, // size of each kernel (should be square)
						MULTIFILE_PSF.validateThreshold);               //      threshold) // to find ellipse
			}

			// Visualize the array as stacks
			int nFiles=kernelsElllipsePars.length;
			int kHeight=kernelsElllipsePars[0].length;
			int kWidth=kernelsElllipsePars[0][0].length;
			int kLength=kHeight*kWidth;
			int nChn=imp_sel.getStack().getSize();
			int numResults=7;
			double [][][][] c= new double[numResults][nChn][nFiles+1][kLength];
			double [][][] numVals=new double[numResults][nChn][kLength];
			int chn, tileY,tileX;
			boolean [] channels=new boolean[nChn];
			double a;
			if (DEBUG_LEVEL>1) { 
				System.out.println("nFiles="+nFiles);
				System.out.println("kWidth="+kWidth);
				System.out.println("kHeight="+kHeight);
				System.out.println("nChn="+nChn);
			}
			Double D;
			int nOut;
			for (chn=0;chn<nChn;chn++) {
				channels[chn]=false;		
				for (nFile=0;nFile<nFiles;nFile++) for (tileY=0;tileY<kHeight;tileY++) for (tileX=0;tileX<kWidth;tileX++) {
					//   			  System.out.println("nChn="+nChn+" nFile="+nFile+" tileY="+tileY+" tileX="+tileX);
					if (kernelsElllipsePars[nFile][tileY][tileX][chn]!=null) {
						channels[chn]=true;
						c[0][chn][nFile+1][tileY*kWidth+tileX]=kernelsElllipsePars[nFile][tileY][tileX][chn][0];
						c[1][chn][nFile+1][tileY*kWidth+tileX]=kernelsElllipsePars[nFile][tileY][tileX][chn][1];
						c[2][chn][nFile+1][tileY*kWidth+tileX]=kernelsElllipsePars[nFile][tileY][tileX][chn][2];
						c[3][chn][nFile+1][tileY*kWidth+tileX]=kernelsElllipsePars[nFile][tileY][tileX][chn][3];
						c[4][chn][nFile+1][tileY*kWidth+tileX]=kernelsElllipsePars[nFile][tileY][tileX][chn][4];
						a=1/Math.sqrt(kernelsElllipsePars[nFile][tileY][tileX][chn][2]*kernelsElllipsePars[nFile][tileY][tileX][chn][3]-
								kernelsElllipsePars[nFile][tileY][tileX][chn][4]*kernelsElllipsePars[nFile][tileY][tileX][chn][4]/4);
						c[5][chn][nFile+1][tileY*kWidth+tileX]= Math.sqrt(a);
						c[6][chn][nFile+1][tileY*kWidth+tileX]=kernelsElllipsePars[nFile][tileY][tileX][chn][5];

					} else {
						c[0][chn][nFile+1][tileY*kWidth+tileX]=Double.NaN;
						c[1][chn][nFile+1][tileY*kWidth+tileX]=Double.NaN;
						c[2][chn][nFile+1][tileY*kWidth+tileX]=Double.NaN;
						c[3][chn][nFile+1][tileY*kWidth+tileX]=Double.NaN;
						c[4][chn][nFile+1][tileY*kWidth+tileX]=Double.NaN;
						c[5][chn][nFile+1][tileY*kWidth+tileX]=Double.NaN;
						c[6][chn][nFile+1][tileY*kWidth+tileX]=Double.NaN;
					}

				}
			}
			/* 
			 * Combine files - now just average all that are not NaN
			 */
			int [][] dirs={{-1,-1},{-1,0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1}};
			int yn,xn,index;
			// remove any tiles that are not OK in all channels
			double [][] weights=new double[nFiles+1][kLength];
			for (nFile=0;nFile<nFiles;nFile++) {
				for (int i=0;i<kLength;i++){
					weights[nFile+1][i]=1.0;
					for (chn=0;chn<nChn;chn++) {
						D=c[0][chn][nFile+1][i];
						if (D.isNaN()) weights[nFile+1][i]=0.0;
					}
				}
				// Set weight to 0.5 if it has zero cells around        	
				for (tileY=0;tileY<kHeight;tileY++) for (tileX=0;tileX<kWidth;tileX++) {
					index=tileY*kWidth+tileX;
					if ( weights[nFile+1][index]>0.0){
						for (int i=0;i<dirs.length;i++) {
							yn=tileY+dirs[i][1];
							xn=tileX+dirs[i][0];
							if ((yn>=0) && (yn<kHeight) && (xn>=0) && (xn<kWidth) && (weights[nFile+1][yn*kWidth+xn]==0.0)){
								weights[nFile+1][index]=0.5;
							}
						}
					}  
					weights[0][index]+=weights[nFile+1][index]; 
				}
			}
			//    	double [][][][] c= new double[numResults][nChn][nFiles+1][kLength];
			//     	double [][][] numVals=new double[numResults][nChn][kLength];
			if (MULTIFILE_PSF.validateShowEllipse) {
				for (chn=0;chn<nChn;chn++) {

					for (nOut=0;nOut<c.length;nOut++) {
						c[nOut][chn][0]=null;
						for (int i=0;i<kLength;i++) {
							numVals[nOut][chn][i]=0.0;
						}
					}
					if (channels[chn]) {
						for (nOut=0;nOut<c.length;nOut++) {
							c[nOut][chn][0]=new double [kLength];
							for (nFile=0;nFile<nFiles;nFile++) {
								for (int i=0;i<kLength;i++){
									D=c[nOut][chn][nFile+1][i];
									if (!D.isNaN()){
										numVals[nOut][chn][i]+=1.0;
										c[nOut][chn][0][i]+=D*weights[nFile+1][i]/weights[0][i];
									}
								}

							}
							for (int i=0;i<kLength;i++){
								if (numVals[nOut][chn][i]==0.0 )c[nOut][chn][0][i]=Double.NaN;
								//    			  else c[nOut][chn][0][i]/=numVals[nOut][chn][i];
							}       	    	
						}
						SDFA_INSTANCE.showArrays(c[0][chn], kWidth, kHeight,  true, "x-shift-"+chn);
						SDFA_INSTANCE.showArrays(c[1][chn], kWidth, kHeight,  true, "y-shift-"+chn);
						SDFA_INSTANCE.showArrays(c[5][chn],kWidth, kHeight,  true, "radius-"+chn);
						if (DEBUG_LEVEL>1) {
							SDFA_INSTANCE.showArrays(c[2][chn], kWidth, kHeight,  true, "x2-"+chn);
							SDFA_INSTANCE.showArrays(c[3][chn], kWidth, kHeight,  true, "y2-"+chn);
							SDFA_INSTANCE.showArrays(c[4][chn], kWidth, kHeight,  true, "xy-"+chn);
							SDFA_INSTANCE.showArrays(c[6][chn], kWidth, kHeight,  true, "area-"+chn);
						}  
					}

				}
				//    		SDFA_INSTANCE.showArrays(weights, kWidth, kHeight,  true, "weights");
			}
			SDFA_INSTANCE.showArrays(weights, kWidth, kHeight,  true, "weights");
			//    	double [][] weights=new double[nFiles+1][kLength];
			for (int i=0;i<kLength;i++) weights[0][i]=0.0;
			PSF_KERNEL_MAP=new double [kHeight][kWidth][nChn][];
			for (tileY=0;tileY<kHeight;tileY++) for (tileX=0;tileX<kWidth;tileX++) for (chn=0;chn<nChn;chn++){
				PSF_KERNEL_MAP[tileY][tileX][chn]=null;
			}
			String [] originalSliceLabels=null;
			String result_path=null;
			for (nFile=0;nFile<nFiles;nFile++) {
				if (UPDATE_STATUS) IJ.showStatus("Accumulating file "+(nFile+1)+" (of "+nFiles+"): "+filenames[nFile]);
				if (DEBUG_LEVEL>1) System.out.println("Accumulating file "+nFile+": "+filenames[nFile]);
				imp_sel=opener.openImage(path, filenames[nFile]);  // or (path+filenames[nFile])
				if (originalSliceLabels==null) {
					originalSliceLabels=imp_sel.getStack().getSliceLabels();
					result_path=imp_sel.getOriginalFileInfo().directory+"PSF-kernels";
				}
				accumulatePartialKernelStack(
						imp_sel.getStack(), // Image stack with partial array of kernels, each slice consists of square kernels of one channel
						INTERPOLATE.size, // size of each kernel (should be square)
						weights[nFile+1], // weights of the kernel tiles in the current stack
						weights[0]);// weights of the kernel tiles already accumulated (will be updated)

			}
			//Finalize accumulated kernels - transform them from frequency to space domain
			inverseTransformKernels();
			int numMissing=0;
			ImageStack mergedStack= mergeKernelsToStack(PSF_KERNEL_MAP,originalSliceLabels);
			System.out.println("mergedStack.getSize()= "+mergedStack.getSize());
			System.out.println("mergedStack.getWidth()= "+mergedStack.getWidth()  );
			System.out.println("mergedStack.getHeight()= "+mergedStack.getHeight()  );
			System.out.println("PSF_KERNEL_MAP.length= "+PSF_KERNEL_MAP.length  );
			System.out.println("PSF_KERNEL_MAP[0].length= "+PSF_KERNEL_MAP[0].length  );
			System.out.println("mergedStack= "+((mergedStack==null)?"null":"not null"));

			if (mergedStack.getSize()==0) {
				System.out.println("*** Error - result is empty");
				return;
			}

			for (tileY=0;tileY<kHeight;tileY++) for (tileX=0;tileX<kWidth;tileX++) if ((PSF_KERNEL_MAP[tileY][tileX]==null) || (PSF_KERNEL_MAP[tileY][tileX][0]==null)) numMissing++;
			ImagePlus imp_psf=SDFA_INSTANCE.showImageStack(mergedStack,imp_sel.getTitle()+"-PSF_KERNEL");
			if (PSF_SAVE_FILE && (result_path!=null)) {
				if (DEBUG_LEVEL>1) {
					System.out.println("Saving result to "+result_path);
				}
				IJ.saveAs(imp_psf,"tif",result_path);
			}
			if (numMissing>0) {
				System.out.println("*** Error "+numMissing+" kernel tiles are missing from the results (insufficient overlap)");
				IJ.showMessage("Error",numMissing+" kernel tiles are missing from the results\n (insufficient overlap)");

			}
			return;   
/* ======================================================================== */
		} else if (label.equals("Invert Stack")) {
			if (!showInverseStackDialog(INVERSE)) return;
			long startTime = System.nanoTime();
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
			ImageStack invertedStack=null;
			invertedStack= reversePSFKernelStack(imp_kernels.getStack(), //  stack of 3 32-bit (float) images, made of square kernels
					INVERSE,
					THREADS_MAX,   // size (side of square) of reverse PSF kernel
					UPDATE_STATUS);         // update status info
			ImagePlus imp_invertedStack = new ImagePlus(imp_kernels.getTitle()+"-rPSF", invertedStack);
			imp_invertedStack.getProcessor().resetMinAndMax();
			imp_invertedStack.show();
			String	result_path=imp_kernels.getOriginalFileInfo().directory+imp_invertedStack.getTitle();
			if (PSF_SAVE_FILE && (result_path!=null)) {
				if (DEBUG_LEVEL>1) {
					System.out.println("Saving inverted result to "+result_path);
				}
				IJ.saveAs(imp_invertedStack,"tif",result_path);
			}
			if (DEBUG_LEVEL>1)  System.out.println("Kernel inversion done in "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+" seconds");
			return;

/* ======================================================================== */
		} else if (label.equals("Gaussian Stack")) {
			if (!showGaussianStackDialog(INVERSE)) return;
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
			ImageStack gaussianStack=  generateGaussianStackFromDirect(imp_kernels.getStack(), // stack of 3 32-bit (float) images, made of square kernels
					INVERSE,
					UPDATE_STATUS);  // update status info
			ImagePlus imp_gaussianStack = new ImagePlus(imp_kernels.getTitle()+"-Gaussian", gaussianStack);
			imp_gaussianStack.getProcessor().resetMinAndMax();
			imp_gaussianStack.show();
			String	result_path=imp_kernels.getOriginalFileInfo().directory+imp_gaussianStack.getTitle();
			if (PSF_SAVE_FILE && (result_path!=null)) {
				if (DEBUG_LEVEL>1) {
					System.out.println("Saving gaussian stack to "+result_path);
				}
				IJ.saveAs(imp_gaussianStack,"tif",result_path);
			}
			return;

/* ======================================================================== */
		} else if (label.equals("Correlation tests")) {
			if (!showDistortionDialog(DISTORTION)) return;
			long 	  startTime=System.nanoTime();
			
//			if ((PATTERN_GRID==null) || (PATTERN_GRID.length!=distortionParameters.gridSize)) {
//				PATTERN_GRID=setPatternGridArray(distortionParameters.gridSize);
//			}
			
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			imp_sel = WindowManager.getCurrentImage();
			if (imp_sel==null){
				IJ.showMessage("Error","There are no images open\nProcess canceled");
				return;
			}
			matchSimulatedPattern= new MatchSimulatedPattern(DISTORTION.FFTSize);
			matchSimulatedPattern.debugLevel=DEBUG_LEVEL;
			matchSimulatedPattern.distortionsTest (
					DISTORTION, //
					PATTERN_DETECT,
					SIMUL,
					COMPONENTS.equalizeGreens,
					imp_sel,
					THREADS_MAX,
					UPDATE_STATUS,
					DISTORTION.loop_debug_level); // debug level
			if (DEBUG_LEVEL>0) System.out.println("Pattern correlation test done at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));

            return;
            /* ======================================================================== */
		} else if (label.equals("Configure Distortion")) {
			if (!showDistortionDialog(DISTORTION)) return;
            return;
            /* ======================================================================== */
		} else if (label.equals("Configure Pointers")) {
			if (!LASER_POINTERS.showDialog("Laser Pointers Setup ")) return;
            return;
/* ======================================================================== */
		} else if (label.equals("Distortion")) {
//			if (!showDistortionDialog(DISTORTION)) return;
			boolean noMessageBoxes=false;
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			imp_sel = WindowManager.getCurrentImage();
//			imp_distortions=imp_sel;
			matchSimulatedPattern= new MatchSimulatedPattern(DISTORTION.FFTSize);
//			   matchSimulatedPattern.invalidateFlatFieldForGrid(); //It is already reset, no need to do it again
//			   matchSimulatedPattern.invalidateFocusMask();
			int numAbsolutePoints=matchSimulatedPattern.calculateDistortions(
					DISTORTION, //
					PATTERN_DETECT,
					SIMUL,
					COMPONENTS.equalizeGreens,
					imp_sel,
					LASER_POINTERS, // LaserPointer laserPointer, // LaserPointer object or null
					true, // don't care -removeOutOfGridPointers
					null, //   double [][][] hintGrid, // predicted grid array (or null)
					0,    //   double  hintGridTolerance, // alllowed mismatch (fraction of period) or 0 - orientation only
					THREADS_MAX,
					UPDATE_STATUS,
					DEBUG_LEVEL,
					DISTORTION.loop_debug_level, // debug level
					noMessageBoxes);
			if (DEBUG_LEVEL>0) System.out.println("Mapped with "+numAbsolutePoints+" absolute points (laser pointers)");
			//			if (numAbsolutePoints>0) {//  no lasers but still crashed in getCalibratedPatternCurvatureAsImage
			ImagePlus imp_curv=matchSimulatedPattern.getCalibratedPatternCurvatureAsImage("curv-"+imp_sel.getTitle());
			imp_curv.show(); // null pointer
			//			if (numAbsolutePoints>0){
			//Calculate grid contrast and brightness for each color component 
			matchSimulatedPattern.calcGridIntensities (
					DISTORTION, //final DistortionParameters distortionParameters, //
					COMPONENTS.equalizeGreens,
					imp_sel, //distortions, // image to process
					THREADS_MAX);				
			ImagePlus imp_calibrated0=matchSimulatedPattern.getCalibratedPatternAsImage(imp_sel,numAbsolutePoints);
			imp_calibrated0.show();
			//			}
			//			} else {
			//				if (DEBUG_LEVEL>0) System.out.println("No absolute mapping info/pointers availabe for this image");
			//			}
			int distWidth=matchSimulatedPattern.getDArrayWidth();
			int distHeight=matchSimulatedPattern.getDArrayHeight();
			if ((distWidth==0) || (distHeight==0)) return;
			double [][] dist=new double [(DEBUG_LEVEL>1)?9:4][distWidth*distHeight];
			double [] distAverage=new double [dist.length];
			int j,k, index;
			for (k=0;k<distAverage.length;k++) distAverage[k]=0.0;
			int nDefined=0;
		    double [][] wv=new double [2][2];
		    
			for (int i=0;i<distHeight;i++) for (j=0;j<distWidth;j++){
				index=j+i*distWidth;
//				for (k=0;k<dist.length;k++) dist[k][index]=0.0;
				if ((matchSimulatedPattern.getDArray(i,j)!=null) && (matchSimulatedPattern.getDArray(i,j,0)!=null)
						&& (matchSimulatedPattern.getDArray(i,j,1)!=null)) {
					wv[0][0]=matchSimulatedPattern.getDArray(i,j,1,0);
					wv[0][1]=matchSimulatedPattern.getDArray(i,j,1,1);
					wv[1][0]=matchSimulatedPattern.getDArray(i,j,2,0);
					wv[1][1]=matchSimulatedPattern.getDArray(i,j,2,1);
					wv=matrix2x2_invert(wv);
					dist[0][index]=wv[0][0];
					dist[1][index]=wv[0][1];
					dist[2][index]=wv[1][0];
					dist[3][index]=wv[1][1];
					if (DEBUG_LEVEL>1) {
						dist[4][index]=Math.sqrt(dist[0][index]*dist[0][index]+dist[1][index]*dist[1][index]);
						dist[5][index]=Math.sqrt(dist[2][index]*dist[2][index]+dist[3][index]*dist[3][index]);
						dist[6][index]=Math.sqrt(dist[4][index]*dist[4][index]+dist[5][index]*dist[5][index]);
						dist[7][index]=matchSimulatedPattern.getDArray(i,j,0,0);
						dist[8][index]=matchSimulatedPattern.getDArray(i,j,0,1);
					}
					for (k=0;k<dist.length;k++) distAverage[k]+=dist[k][index];
					nDefined++;
				}
				
			}
			for (k=0;k<dist.length;k++) distAverage[k]/=nDefined;
	
			for (int i=0;i<distHeight;i++) for (j=0;j<distWidth;j++){
				index=j+i*distWidth;
				if ((matchSimulatedPattern.getDArray(i,j)==null) || (matchSimulatedPattern.getDArray(i,j,0)==null)
						|| (matchSimulatedPattern.getDArray(i,j,1)==null))
				for (k=0;k<dist.length;k++) dist[k][index]=distAverage[k];
			}			
			if (DEBUG_LEVEL>1) SDFA_INSTANCE.showArrays(dist, distWidth, distHeight, true, "distortions");
			if (DEBUG_LEVEL>4) SDFA_INSTANCE.showArrays(dist, distWidth, distHeight, false, "distortions");
			if (imp_sel.getStackSize()>1) {
// make calibration by laser pointers
				calcLaser(imp_sel,noMessageBoxes);
				if (DEBUG_LEVEL>1) {
					ImagePlus imp_calibrated=matchSimulatedPattern.getCalibratedPatternAsImage(imp_sel.getTitle()+"-grid",0);
					imp_calibrated.getProcessor().resetMinAndMax();
					imp_calibrated.show();
				}
			}
			
            return;
/* ======================================================================== */
		} else if (label.equals("Grid Brightness")) {
//			if (!showDistortionDialog(DISTORTION)) return;
//			boolean noMessageBoxes=false;
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (matchSimulatedPattern==null) {
				IJ.showMessage("matchSimulatedPattern is not defined. Run \"Distortion\"");
				return;
			}
//			imp_distortions;
			imp_sel = WindowManager.getCurrentImage();
			GenericDialog gd=new GenericDialog("Grid brightness calculation");
			double agp=matchSimulatedPattern.averageGridPeriod(matchSimulatedPattern.PATTERN_GRID);
			gd.addMessage("Average grid.period="+ agp);
			gd.addNumericField ("Scale average grid period when averaging brightness",DISTORTION.averagingAreaScale,2,4,"x");
    		gd.showDialog();
    		if (gd.wasCanceled()) return;
    		DISTORTION.averagingAreaScale= gd.getNextNumber();
    		String [] titles={"Contrast","Red","Green","Blue", "R/G", "B/G","B/Rr"};
    		double [][][] intensities =matchSimulatedPattern.calcGridIntensities (
 				   DISTORTION, //final DistortionParameters distortionParameters, //
 				   COMPONENTS.equalizeGreens,
 				   imp_sel, //imp_distortions, // image to process
 				   THREADS_MAX);
			double [][] showGridIintensity=new double [titles.length][intensities[0].length*intensities[0][0].length];
			for (int n=0;n<intensities.length;n++){
				int index=0;
				for (int i=0;i<intensities[n].length;i++) for (int j=0;j<intensities[n][i].length;j++) showGridIintensity[n][index++]=intensities[n][i][j];
				
			}
			showGridIintensity[4]=showGridIintensity[1].clone();
			showGridIintensity[5]=showGridIintensity[3].clone();
			showGridIintensity[6]=showGridIintensity[3].clone();
			
			for (int i=0;i<showGridIintensity[4].length;i++) {
				showGridIintensity[4][i]/=showGridIintensity[2][i]; // r/g
				showGridIintensity[5][i]/=showGridIintensity[2][i]; // b/g
				showGridIintensity[6][i]/=showGridIintensity[1][i]; // b/r
			}
			SDFA_INSTANCE.showArrays(showGridIintensity, intensities[0][0].length,intensities[0].length, true,"gridIntensities", titles);
			return;
/* ======================================================================== */
		} else if (label.equals("List Pointers")) {
			boolean noMessageBoxes=false;
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
//			if (!LASER_POINTERS.showDialog("Laser Pointers Setup ")) return;
			ImagePlus imp_pointers = WindowManager.getCurrentImage();
			if (imp_pointers==null){
				IJ.showMessage("Error","There is no stack with the laser pointer images to process");
				return;
			}
			if (imp_pointers.getStackSize()<2) {
				IJ.showMessage("Error","Need at least a 2-layer stack with the laser pointer images");
				return;
			}
// make calibration by laser pointers
			calcLaser(imp_pointers,noMessageBoxes);
			if (DEBUG_LEVEL>1) {
				ImagePlus imp_calibrated=matchSimulatedPattern.getCalibratedPatternAsImage(imp_sel.getTitle()+"-grid",0);
				imp_calibrated.getProcessor().resetMinAndMax();
				imp_calibrated.show();
			}

            return;
            
/* ======================================================================== */
		} else if (label.equals("Simulate Full")) {
			if ((matchSimulatedPattern==null) || (matchSimulatedPattern.patternOK())){
				IJ.showMessage("Error","Distortion data is not yet calculated\nProcess canceled");
				return;
			}
/*			
			SIM_ARRAY=	simulateGridAll (
					matchSimulatedPattern.getImageWidth(),
					matchSimulatedPattern.getImageHeight(),
					matchSimulatedPattern.getDArray(),
					2, // gridFrac, // number of grid steps per pattern full period
					SIMUL,
					THREADS_MAX,
					UPDATE_STATUS,
					DISTORTION.loop_debug_level); // debug level
*/			
			SIM_ARRAY=	(new SimulationPattern(SIMUL)).simulateGridAll (
					matchSimulatedPattern.getImageWidth(),
					matchSimulatedPattern.getImageHeight(),
					matchSimulatedPattern,
					2, // gridFrac, // number of grid steps per pattern full period
					SIMUL,
					THREADS_MAX,
					UPDATE_STATUS,
					DEBUG_LEVEL,
					DISTORTION.loop_debug_level); // debug level
			return;
		}
/* ======================================================================== */
		if       (label.equals("Configure Flat Field")) {
			if (showFlatFieldDialog(FLATFIELD_PARAMETERS,true,false)) showFlatFieldDialog(FLATFIELD_PARAMETERS,false,true);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Remove FF Source Dirs")) {
			selectDirectoriesToRemove(FLATFIELD_PARAMETERS);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Add FF Source Dirs")) {
			addFlatFieldSources(FLATFIELD_PARAMETERS);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Accumulate Flat Fields")) {
			if (!showFlatFieldDialog(FLATFIELD_PARAMETERS, true,false)) return;
	        DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			flatFieldAccumulate(FLATFIELD_PARAMETERS,PROCESS_PARAMETERS); // returns true if did accumulate
			if (FLATFIELD_PARAMETERS.saveSettings)
				saveProperties(FLATFIELD_PARAMETERS.flatFieldDirectory+Prefs.getFileSeparator()+"flatfield_settings",null,FLATFIELD_PARAMETERS.useXML, PROPERTIES);

			return;
		}
/* ======================================================================== */
		if       (label.equals("Process Flat Field")) {
			if (!showFlatFieldDialog(FLATFIELD_PARAMETERS,false)) return;
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			imp_sel = WindowManager.getCurrentImage();
			if (imp_sel==null){
				IJ.showMessage("Error","There are no images open\nProcess canceled");
				return;
			}
            ImageProcessor ip_ff=imp_sel.getProcessor();
            float [] pixels_ff= (float[]) ip_ff.getPixels();
            processFlatField (FLATFIELD_PARAMETERS,pixels_ff, imp_sel.getWidth(), imp_sel.getTitle());
			return;
		}
/* ======================================================================== */
		if       (label.equals("Configure Pattern Dimensions")) {
			PATTERN_PARAMETERS.showDialog();
			return;
		}
/* ======================================================================== */
		if       (label.equals("Configure Distortion/Location")) {
			LENS_DISTORTION_PARAMETERS.showDialog();
			return;
		}
/* ======================================================================== */
		if       (label.equals("Configure Eyesis4pi")) {
			EYESIS_CAMERA_PARAMETERS.showDialog("Configure Eyesis4pi camera"); // returns false
			if (LENS_DISTORTIONS!=null) LENS_DISTORTIONS.resetGridImageMasks(); // will force recalculation
			return;
		}
/* ======================================================================== */
		if       (label.equals("List Eyesis4pi")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			System.out.println("DISTORTION_CALIBRATION_DATA="+((DISTORTION_CALIBRATION_DATA==null)?"null":" NOT null"));
			System.out.println("+++++++++++ EYESIS_CAMERA_PARAMETERS.numStations="+EYESIS_CAMERA_PARAMETERS.numStations+
					" +EYESIS_CAMERA_PARAMETERS.goniometerHorizontal.length="+EYESIS_CAMERA_PARAMETERS.goniometerHorizontal.length);
			if (DISTORTION_CALIBRATION_DATA!=null) {
			if (DEBUG_LEVEL>1) System.out.println("+++++++++++ DISTORTION_CALIBRATION_DATA.eyesisCameraParameters.numStations="+DISTORTION_CALIBRATION_DATA.eyesisCameraParameters.numStations+
					" +DISTORTION_CALIBRATION_DATA.eyesisCameraParameters.goniometerHorizontal.length="+DISTORTION_CALIBRATION_DATA.eyesisCameraParameters.goniometerHorizontal.length);
			}

			Distortions.DistortionCalibrationData dcd=(DISTORTION_CALIBRATION_DATA!=null)?DISTORTION_CALIBRATION_DATA:
				new Distortions.DistortionCalibrationData(EYESIS_CAMERA_PARAMETERS);
			if (DEBUG_LEVEL>1) System.out.println("+++++++++++ dcd.eyesisCameraParameters.numStations="+dcd.eyesisCameraParameters.numStations+
					" +dcd.eyesisCameraParameters.goniometerHorizontal.length="+dcd.eyesisCameraParameters.goniometerHorizontal.length);
			dcd.listCameraParameters();
			return;
		}
		
/* ======================================================================== */
		if       (label.equals("Process Lens Distortion")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			LENS_DISTORTIONS=new Distortions(LENS_DISTORTION_PARAMETERS,PATTERN_PARAMETERS,REFINE_PARAMETERS,this.SYNC_COMMAND.stopRequested);
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			double[][] gridDisplay=LENS_DISTORTIONS.prepareDisplayGrid();
			SDFA_INSTANCE.showArrays(gridDisplay,LENS_DISTORTIONS.getGridWidth(),LENS_DISTORTIONS.getGridHeight(),  true, "Grid",
					LENS_DISTORTIONS.displayGridTitles());
			LENS_DISTORTIONS.calcGridOnSensor(0.0);
			double[][] gridDisplayOnSensor=LENS_DISTORTIONS.prepareDisplayGridOnSensor(false);
			SDFA_INSTANCE.showArrays(gridDisplayOnSensor,LENS_DISTORTIONS.getGridWidth(),LENS_DISTORTIONS.getGridHeight(),  true, "GridOnSensor",
					LENS_DISTORTIONS.displayGridOnSensorTitles());
			double[][] gridDisplayOnSensorAll=LENS_DISTORTIONS.prepareDisplayGridOnSensor(true);
			SDFA_INSTANCE.showArrays(gridDisplayOnSensorAll,LENS_DISTORTIONS.getGridWidth(),LENS_DISTORTIONS.getGridHeight(),  true, "GridOnSensor-all",
					LENS_DISTORTIONS.displayGridOnSensorTitles());
			LENS_DISTORTIONS.calcGridOnSensor(0.01);
			double[][] gridDisplayOnSensorAllDelta=LENS_DISTORTIONS.prepareDisplayGridOnSensor(true);
			SDFA_INSTANCE.showArrays(gridDisplayOnSensorAllDelta,LENS_DISTORTIONS.getGridWidth(),LENS_DISTORTIONS.getGridHeight(),  true, "GridOnSensor-all-0.01",
					LENS_DISTORTIONS.displayGridOnSensorTitles());
			return;
		}
/* ======================================================================== */
		if       (label.equals("Select Grid Files")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			LENS_DISTORTIONS=new Distortions(LENS_DISTORTION_PARAMETERS,PATTERN_PARAMETERS,REFINE_PARAMETERS,this.SYNC_COMMAND.stopRequested);
			String [] extensions={".tif",".tiff"};
			CalibrationFileManagement.MultipleExtensionsFileFilter gridFilter = 
				new CalibrationFileManagement.MultipleExtensionsFileFilter("grid",extensions,"Calibrated grid files");
		    GenericDialog gd = new GenericDialog("Setup Goniometer/Camera Stations");
		    gd.addMessage("Setting up calibration that includes multiple camera tripod or goniometer positions.");
		    gd.addMessage("File selection dialog will open for each station separateley.");
		    gd.addNumericField("Number of goniometer/camera stations", 1,0);
		    gd.showDialog();
		    if (gd.wasCanceled()) return;
		    int numStations= (int) gd.getNextNumber();
	    	String [][] gridFiles= new String [numStations][];
		    for (int numStation=0;numStation<numStations;numStation++){
		    	gridFiles[numStation]=CalibrationFileManagement.selectFiles(false,
						"Select Calibrated Grid Files for Station "+numStation+ "("+(numStation+1)+" of "+numStations+")",
						"Select",
						gridFilter,
						null); // String [] defaultPaths);
		       	if ((gridFiles[numStation]==null) || (gridFiles[numStation].length==0)) {
	        		if (!IJ.showMessageWithCancel("No files selected","Retry? (Cancel will abort the command)")) return;
	        		numStation--;
	        	}
		    }
/*			
			String [] gridFiles=CalibrationFileManagement.selectFiles(false,
					"Select Calibrated Grid Files",
					"Select",
					gridFilter,
					null); // String [] defaultPaths);
	       	if ((gridFiles==null) || (gridFiles.length==0)) {
        		IJ.showMessage("No files selected");
        		return;
        	}
*/        	
			PATTERN_PARAMETERS.debugLevel=MASTER_DEBUG_LEVEL;
			EYESIS_CAMERA_PARAMETERS.updateNumstations (numStations);
//if (MASTER_DEBUG_LEVEL==0) return; //TODO: Remove - just debugging
			DISTORTION_CALIBRATION_DATA=new Distortions.DistortionCalibrationData(
					gridFiles,
					PATTERN_PARAMETERS,
					EYESIS_CAMERA_PARAMETERS,
					MASTER_DEBUG_LEVEL);
			LENS_DISTORTIONS.initImageSet(
 					DISTORTION_CALIBRATION_DATA,
					EYESIS_CAMERA_PARAMETERS
					);
			// set initial orientation of the cameras from the sensors that see mpst of the matching pointers
			// just for the LMA to start
			DISTORTION_CALIBRATION_DATA.setInitialOrientation(true);
			return;
		}
		
/* ======================================================================== */
/*
  		if       (label.equals("Use Current Images")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			
			LENS_DISTORTIONS=new Distortions(LENS_DISTORTION_PARAMETERS,PATTERN_PARAMETERS,REFINE_PARAMETERS);
			DISTORTION_CALIBRATION_DATA=new Distortions.DistortionCalibrationData(gridFiles,PATTERN_PARAMETERS,EYESIS_CAMERA_PARAMETERS);
			LENS_DISTORTIONS.initImageSet(
 					DISTORTION_CALIBRATION_DATA,
					EYESIS_CAMERA_PARAMETERS
					);
			return;
		}
*/		
/* ======================================================================== */
		if       (label.equals("Edit Calibration")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (DISTORTION_CALIBRATION_DATA==null) {
				IJ.showMessage("Distortion calibration data required");
				return;
			}
            int numImage=0;
            while (numImage>=0) {
            	numImage=DISTORTION_CALIBRATION_DATA.editImageParameters(numImage);
            }
// numImage =-1 - "done", <-1 - "canceled"            
			return;
		}
/* ======================================================================== */
		if       (label.equals("Save Calibration")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (DISTORTION_CALIBRATION_DATA==null) return;
			DISTORTION_CALIBRATION_DATA.selectAndSaveToXML(false,"");
			return;
		}
/* ======================================================================== */
		if       (label.equals("Restore Calibration")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			String defaultPath= (DISTORTION_CALIBRATION_DATA!=null)?DISTORTION_CALIBRATION_DATA.pathName:"";
//    public static Distortions.DistortionCalibrationData DISTORTION_CALIBRATION_DATA=null; 
			Distortions.DistortionCalibrationData oldDISTORTION_CALIBRATION_DATA=DISTORTION_CALIBRATION_DATA;
			DISTORTION_CALIBRATION_DATA=new Distortions.DistortionCalibrationData(
					false,
					defaultPath,
					PATTERN_PARAMETERS,
					EYESIS_CAMERA_PARAMETERS,
					null); // gridImages null - use specified files
			if (DISTORTION_CALIBRATION_DATA.pathName== null){ // failed to select/open the file
				DISTORTION_CALIBRATION_DATA=oldDISTORTION_CALIBRATION_DATA;
				return;
			}
			// Update referenced data
			if ((LENS_DISTORTIONS!=null) && (LENS_DISTORTIONS.fittingStrategy!=null)) {
				LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData=DISTORTION_CALIBRATION_DATA;
			}
			return;
		}
/* ======================================================================== */
		if       (label.equals("New Strategy")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (DISTORTION_CALIBRATION_DATA==null) {
				IJ.showMessage("Distortion calibration data required");
				return;
			}
			if (LENS_DISTORTIONS==null) {
				if (LENS_DISTORTION_PARAMETERS==null){
					IJ.showMessage("LENS_DISTORTION_PARAMETERS is not set");
					return;
				}
				if (PATTERN_PARAMETERS==null){
					IJ.showMessage("PATTERN_PARAMETERS is not set");
					return;
				}
				LENS_DISTORTIONS=new Distortions(LENS_DISTORTION_PARAMETERS,PATTERN_PARAMETERS,REFINE_PARAMETERS,this.SYNC_COMMAND.stopRequested);
			}
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			if (DEBUG_LEVEL>2) System.out.println("New Strategy");
			LENS_DISTORTIONS.fittingStrategy=new Distortions.FittingStrategy(DISTORTION_CALIBRATION_DATA);
			LENS_DISTORTIONS.fittingStrategy.debugLevel=DEBUG_LEVEL;
			IJ.showMessage("Empty new strategy initialized");
			return;
		}
/* ======================================================================== */
		
		if       (label.equals("Edit Strategy")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				if (LENS_DISTORTION_PARAMETERS==null){
					IJ.showMessage("LENS_DISTORTION_PARAMETERS is not set");
					return;
				}
				if (PATTERN_PARAMETERS==null){
					IJ.showMessage("PATTERN_PARAMETERS is not set");
					return;
				}
				LENS_DISTORTIONS=new Distortions(LENS_DISTORTION_PARAMETERS,PATTERN_PARAMETERS,REFINE_PARAMETERS,this.SYNC_COMMAND.stopRequested);
			}
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			if (LENS_DISTORTIONS.fittingStrategy==null) {
				if (DISTORTION_CALIBRATION_DATA==null) return;
				LENS_DISTORTIONS.fittingStrategy=new Distortions.FittingStrategy(DISTORTION_CALIBRATION_DATA);
			}
			LENS_DISTORTIONS.fittingStrategy.debugLevel=DEBUG_LEVEL;
			LENS_DISTORTIONS.fittingStrategy.selectStrategy(LENS_DISTORTIONS.seriesNumber);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Save Strategy")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			if (LENS_DISTORTIONS.fittingStrategy==null) return;
			LENS_DISTORTIONS.fittingStrategy.debugLevel=DEBUG_LEVEL;
			LENS_DISTORTIONS.fittingStrategy.selectAndSaveToXML(false,null);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Restore Strategy")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (DISTORTION_CALIBRATION_DATA==null) {
				IJ.showMessage("Distortion calibration data required");
				return;
			}
			if (LENS_DISTORTIONS==null) {
				if (LENS_DISTORTION_PARAMETERS==null){
					IJ.showMessage("LENS_DISTORTION_PARAMETERS is not set");
					return;
				}
				if (PATTERN_PARAMETERS==null){
					IJ.showMessage("PATTERN_PARAMETERS is not set");
					return;
				}
				LENS_DISTORTIONS=new Distortions(LENS_DISTORTION_PARAMETERS,PATTERN_PARAMETERS,REFINE_PARAMETERS,this.SYNC_COMMAND.stopRequested);
			}
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			Distortions.FittingStrategy fs=LENS_DISTORTIONS.fittingStrategy; // save old value
			String defaultPath= ((fs!=null) && (fs.pathName != null) && (fs.pathName.length()>0)) ? fs.pathName : "";
			LENS_DISTORTIONS.fittingStrategy=new Distortions.FittingStrategy(
					false,
					defaultPath,
					DISTORTION_CALIBRATION_DATA);
			if (LENS_DISTORTIONS.fittingStrategy.pathName== null){ // failed to select/open the file
				LENS_DISTORTIONS.fittingStrategy=null;
				IJ.showMessage("Nothing selected");
				LENS_DISTORTIONS.fittingStrategy=fs; // restore old strategy
				return;
			}
			LENS_DISTORTIONS.fittingStrategy.debugLevel=DEBUG_LEVEL;
			LENS_DISTORTIONS.fittingStrategy.adjustNumberOfImages(DISTORTION_CALIBRATION_DATA.gIP.length);
			//		lensDistortions.fittingStrategy.adjustNumberOfImages(imp_calibrated.length);

			return;
		}
/* ======================================================================== */
		if       (label.equals("Run LMA")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			LENS_DISTORTIONS.updateStatus=UPDATE_STATUS;
			if (LENS_DISTORTIONS.fittingStrategy==null) return;
			LENS_DISTORTIONS.fittingStrategy.debugLevel=DEBUG_LEVEL;
			
			LENS_DISTORTIONS.LevenbergMarquardt(true); // open dialog
			return;
		}
/* ======================================================================== */
		if       (label.equals("Bad Nodes")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			LENS_DISTORTIONS.updateStatus=UPDATE_STATUS;
			if (LENS_DISTORTIONS.fittingStrategy==null) return;
			LENS_DISTORTIONS.fittingStrategy.debugLevel=DEBUG_LEVEL;
			LENS_DISTORTIONS.dialogMarkBadNodes(DEBUG_LEVEL);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Debug deriv")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			if (LENS_DISTORTIONS.fittingStrategy==null) return;
			LENS_DISTORTIONS.fittingStrategy.debugLevel=DEBUG_LEVEL;
			LENS_DISTORTIONS.compareDerivatives();
			return;
		}

/* ======================================================================== */
		if       (label.equals("Configure Lasers")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			LASERS.debugLevel=DEBUG_LEVEL;
			LASERS.showDialog("Configure laser pointers hardware", true);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Manual laser pointers")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			LASERS.debugLevel=DEBUG_LEVEL;
			for (;	LASERS.manualSetLasers(););
			return;
		}
/* ======================================================================== */
		if       (label.equals("Configure Cameras")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			CAMERAS.debugLevel=DEBUG_LEVEL;
			CAMERAS.showDialog("Configure cameras interface", 0, true);
			if (!FOCUS_MEASUREMENT_PARAMETERS.getLensSerial()) return;
			// reset histories
			MOTORS.clearPreFocus();
			MOTORS.clearHistory();
			return;
		}
/* ======================================================================== */
		if       (label.equals("Cameras Settings")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			CAMERAS.debugLevel=DEBUG_LEVEL;
			if (!CAMERAS.editCameraSettings("Camera acquisition settings")) return;
			
			CAMERAS.probeCameraState(); // testing detection
			CAMERAS.setupCameraAcquisition();
			return;
		}
/* ======================================================================== */
		if       (label.equals("Test Cameras")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			CAMERAS.debugLevel=DEBUG_LEVEL;
			CAMERAS.setNumberOfThreads(THREADS_MAX);
			CAMERAS.test1(true);
			
			return;
		}
/* ======================================================================== */
		if       (label.equals("Test No Lasers")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			CAMERAS.debugLevel=DEBUG_LEVEL;
			CAMERAS.setNumberOfThreads(THREADS_MAX);
			CAMERAS.test1(false);
			return;
		}
		
/* ======================================================================== */
		if       (label.equals("Configure Process Distortions")) {
			DISTORTION_PROCESS_CONFIGURATION.debugLevel=MASTER_DEBUG_LEVEL;
			if (DISTORTION_PROCESS_CONFIGURATION.showDialog("")){
				MASTER_DEBUG_LEVEL=DISTORTION_PROCESS_CONFIGURATION.debugLevel;
			}
			return;
		}
/* ======================================================================== */
		if       (label.equals("Quick get&show")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			CAMERAS.debugLevel=DEBUG_LEVEL;
			checkSerialAndRestore(); // returns true if did not change or was restored 
			long 	  startTime=System.nanoTime();
			if (!FOCUS_MEASUREMENT_PARAMETERS.cameraIsConfigured) {
				if (CAMERAS.showDialog("Configure cameras interface", 1, true)){
					FOCUS_MEASUREMENT_PARAMETERS.cameraIsConfigured=true;
					if (!FOCUS_MEASUREMENT_PARAMETERS.getLensSerial()) return;
					// reset histories
					MOTORS.clearPreFocus();
					MOTORS.clearHistory();
				} else {
					IJ.showMessage("Error","Camera is not configured\nProcess canceled");
					return;
				}
			}
			// acquire camera image here, no lasers.
			FOCUS_MEASUREMENT_PARAMETERS.sensorTemperature=CAMERAS.getSensorTemperature(0,FOCUS_MEASUREMENT_PARAMETERS.EEPROM_channel);
			imp_sel= CAMERAS.acquireSingleImage (
					false, //boolean useLasers,
					UPDATE_STATUS);
			if (imp_sel==null){
				IJ.showMessage("Error","Failed to get camera image\nProcess canceled");
				return;
			}
			imp_sel.show();
			imp_sel.updateAndDraw();
			if (DEBUG_LEVEL>0) System.out.println("Image acquisition (@"+FOCUS_MEASUREMENT_PARAMETERS.sensorTemperature+"C) done at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
			return;
		}
		
		
		
/* ======================================================================== */
		if       (label.equals("Acquire")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			CAMERAS.debugLevel=DEBUG_LEVEL;
			CAMERAS.reportTiming=DEBUG_LEVEL>0;
			String src_dir=DISTORTION_PROCESS_CONFIGURATION.selectSourceDirectory(true, DISTORTION_PROCESS_CONFIGURATION.sourceDirectory, true);
			if (src_dir==null) {
		    	IJ.showMessage("Error", "Nothing selected\nProcess canceled");
				return;
			}
			DISTORTION_PROCESS_CONFIGURATION.sourceDirectory=src_dir;
			CAMERAS.setNumberOfThreads(THREADS_MAX);
			long startTime=System.nanoTime();
			CAMERAS.acquire(DISTORTION_PROCESS_CONFIGURATION.sourceDirectory,true, UPDATE_STATUS); // true - use lasers
			System.out.println("\"Acquire\" command finished in ("+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+" sec");
			return;
		}
/* ======================================================================== */
		if       (label.equals("Grid candidate")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			imp_sel = WindowManager.getCurrentImage();
			if (imp_sel==null){
				IJ.showMessage("Error","There is no image selected");
				return;
			}
			LASER_POINTERS.showFilterDialog("Pattern filter parameters");
			
			float [] fpixels= (float[]) imp_sel.getProcessor().getPixels();
			int imp_sel_width=imp_sel.getWidth();
			double [] dpixels=new double [fpixels.length];
			for (int i=0;i<fpixels.length;i++) dpixels[i]=fpixels[i];
			boolean [] patternMask=LASER_POINTERS.getPatternMask(
					dpixels,
					imp_sel_width
					);
			if (patternMask!=null) {
				SDFA_INSTANCE.showArrays(patternMask, imp_sel_width, dpixels.length/imp_sel_width,  "pattern_mask");
			}
			return;
		}
//		
		
/* ======================================================================== */
		if       (label.equals("Calculate grids")) {
		    if ((LASER_POINTERS==null) || (LASER_POINTERS.laserUVMap.length==0)){
		    	IJ.showMessage("Laser pointer data needed for this function is not provided");
		    	return;
		    }
			long 	  startTime=System.nanoTime();
		    boolean noMessageBoxes=true;
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			DISTORTION_PROCESS_CONFIGURATION.debugLevel=MASTER_DEBUG_LEVEL;
            if (matchSimulatedPattern==null) matchSimulatedPattern= new MatchSimulatedPattern(DISTORTION.FFTSize);
            matchSimulatedPattern.debugLevel=MASTER_DEBUG_LEVEL;
            String [] sourceFilesList=DISTORTION_PROCESS_CONFIGURATION.selectSourceFiles(); // select files - with/without dialog
            boolean saveGrids=DISTORTION_PROCESS_CONFIGURATION.saveGridImages;
            if (sourceFilesList==null) return;
            showPatternMinMaxPeriodDialog(PATTERN_DETECT);
            for (int numFile=0;numFile<sourceFilesList.length;numFile++){
    			long 	  startFileTime=System.nanoTime();
            	if (DEBUG_LEVEL>0){
            		System.out.println(IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+"s: Processing file # "+(numFile+1)+ " (of "+ sourceFilesList.length+"): "+sourceFilesList[numFile]);
            	}
            	imp_sel=new ImagePlus(sourceFilesList[numFile]); // read source file
            	JP4_INSTANCE.decodeProperiesFromInfo(imp_sel);
            	matchSimulatedPattern= new MatchSimulatedPattern(DISTORTION.FFTSize); // TODO: is it needed each time?
            	if (!DISTORTION_PROCESS_CONFIGURATION.useNoPonters && (matchSimulatedPattern.getPointersXY(imp_sel,LASER_POINTERS.laserUVMap.length)==null)) {
    				if (this.SYNC_COMMAND.stopRequested.get()>0) {
    					System.out.println("User requested stop");
    					break;
    				}
            		continue; // no pointers in this image
            	}
            	// /getPointersXY(ImagePlus imp, int numPointers){               if
            	// calculate distortion grid for it

            	matchSimulatedPattern.invalidateFlatFieldForGrid(); //Reset Flat Field calibration - different image. 
            	matchSimulatedPattern.invalidateFocusMask();
            	int numAbsolutePoints=matchSimulatedPattern.calculateDistortions(
            			DISTORTION, //
            			PATTERN_DETECT,
            			SIMUL,
            			COMPONENTS.equalizeGreens,
            			imp_sel,
            			LASER_POINTERS, // LaserPointer laserPointer, // LaserPointer object or null
            			DISTORTION_PROCESS_CONFIGURATION.removeOutOfGridPointers, // 
            			null, //   double [][][] hintGrid, // predicted grid array (or null)
            			0,    //   double  hintGridTolerance, // allowed mismatch (fraction of period) or 0 - orientation only
            			THREADS_MAX,
            			UPDATE_STATUS,
            			DEBUG_LEVEL,
            			DISTORTION.loop_debug_level, // debug level
            			noMessageBoxes);
            	if (DEBUG_LEVEL>1) System.out.println("numAbsolutePoints="+numAbsolutePoints);
            	if ((numAbsolutePoints==DISTORTION.errPatternNotFound) || (numAbsolutePoints==DISTORTION.errTooFewCells)) {
    				if (DEBUG_LEVEL>0) System.out.println("Grid "+(numFile+1)+" not found or too small ("+numAbsolutePoints+"), wasted "+
    						IJ.d2s(0.000000001*(System.nanoTime()-startFileTime),3)+" seconds )\n");
    				if (this.SYNC_COMMAND.stopRequested.get()>0) {
    					System.out.println("User requested stop");
    					break;
    				}
            		continue; // too few cells detected
            	}
            	if (DISTORTION_PROCESS_CONFIGURATION.useNoPonters || (numAbsolutePoints>0)){
        			//Calculate grid contrast and brightness for each color component 
        			matchSimulatedPattern.calcGridIntensities (
        					DISTORTION, //final DistortionParameters distortionParameters, //
        					COMPONENTS.equalizeGreens,
        					imp_sel, // image to process
        					THREADS_MAX);				
            		ImagePlus imp_calibrated=matchSimulatedPattern.getCalibratedPatternAsImage(imp_sel,numAbsolutePoints);
            		if (DISTORTION_PROCESS_CONFIGURATION.showGridImages) imp_calibrated.show();
            		if (saveGrids){
            			FileSaver fs=new FileSaver(imp_calibrated);
            			String srcDir=DISTORTION_PROCESS_CONFIGURATION.selectGridFileDirectory(true,DISTORTION_PROCESS_CONFIGURATION.gridDirectory,true);
            			if (srcDir==null){
            				saveGrids=false; // do not ask about the next ones too
            			} else {
            				String path=DISTORTION_PROCESS_CONFIGURATION.gridDirectory+Prefs.getFileSeparator()+imp_calibrated.getTitle();
            				if (UPDATE_STATUS) IJ.showStatus("Saving "+path);
            				if (DEBUG_LEVEL>0) System.out.println("-->>> Saving "+path+" - using "+numAbsolutePoints+" laser pointer references");
            				fs.saveAsTiffStack(path);
            			}
            		}
            	}
				if (DEBUG_LEVEL>0) System.out.println("Grid "+(numFile+1)+" calculation done at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+" (in "+
						IJ.d2s(0.000000001*(System.nanoTime()-startFileTime),3)+"s )\n");
				
// 				
				if (this.SYNC_COMMAND.stopRequested.get()>0) {
					System.out.println("User requested stop");
					break;
				}

            }
			if (DEBUG_LEVEL>0) System.out.println(((this.SYNC_COMMAND.stopRequested.get()>0)?"Partial (interrupted by user) set of grids":"All")+ " grids calculation done at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
            return;
		}
/* ======================================================================== */
		if       (label.equals("Calculate Sensor Masks")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
		    if (DISTORTION_CALIBRATION_DATA==null){
		    	IJ.showMessage("Distortion Calibration data is not initialized - create it with"+
		    			"\"SelectGrid Files\" or \"Restore Calibration\"");
		    	return;
		    }
			DISTORTION_CALIBRATION_DATA.debugLevel=DEBUG_LEVEL;
		    DISTORTION_CALIBRATION_DATA.updateStatus=UPDATE_STATUS;
			DISTORTION_CALIBRATION_DATA.calculateSensorMasks();
			if (LENS_DISTORTIONS!=null) LENS_DISTORTIONS.updateSensorMasks();
			return;
		}
/* ======================================================================== */
		if       (label.equals("Save Sensor Masks")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
		    if (DISTORTION_CALIBRATION_DATA==null){
		    	IJ.showMessage("Distortion Calibration data is not initialized - create it with"+
		    			"\"SelectGrid Files\" or \"Restore Calibration\"");
		    	return;
		    }
			DISTORTION_CALIBRATION_DATA.debugLevel=DEBUG_LEVEL;
		    DISTORTION_CALIBRATION_DATA.updateStatus=UPDATE_STATUS;
		    if (DISTORTION_CALIBRATION_DATA.sensorMasks==null){
		    	IJ.showMessage("Sensor masks are not initialized - create it with"+
		    			"\"Calculate Sensor Masks\" or \"Restore Sensor Masks\"");
		    	return;
		    }
			String [] extensions={".mask-tiff","-masks.tiff"};
			CalibrationFileManagement.MultipleExtensionsFileFilter parFilter = new CalibrationFileManagement.MultipleExtensionsFileFilter("",extensions,"*.mask-tiff files");
			String pathname=CalibrationFileManagement.selectFile(true,
					"Save Sensor Masks file",
					"Save",
					parFilter,
					""); //String defaultPath
			DISTORTION_CALIBRATION_DATA.saveMaskAsImageStack("Sensor Masks", pathname);
			
			return;
		}
/* ======================================================================== */
		if       (label.equals("Restore Sensor Masks")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
		    if (DISTORTION_CALIBRATION_DATA==null){
		    	IJ.showMessage("Distortion Calibration data is not initialized - create it with"+
		    			"\"SelectGrid Files\" or \"Restore Calibration\"");
		    	return;
		    }
		    DISTORTION_CALIBRATION_DATA.debugLevel=DEBUG_LEVEL;
		    DISTORTION_CALIBRATION_DATA.updateStatus=UPDATE_STATUS;
			String [] extensions={".mask-tiff","-masks.tiff"};
			
			CalibrationFileManagement.MultipleExtensionsFileFilter parFilter = new CalibrationFileManagement.MultipleExtensionsFileFilter("",extensions,"Sensor masks *.mask-tiff files");
			String pathname=CalibrationFileManagement.selectFile(false,
					"Restore Sensor Masks",
					"Restore",
					parFilter,
					""); //String defaultPath
			if ((pathname==null) || (pathname=="")) return;
		    DISTORTION_CALIBRATION_DATA.setMaskFromImageStack(pathname);
			if (LENS_DISTORTIONS!=null) LENS_DISTORTIONS.updateSensorMasks();
			return;
		}
/* ======================================================================== */
		if       (label.equals("Reset Grid")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
		    if (PATTERN_PARAMETERS==null){
		    	IJ.showMessage("PATTERN_PARAMETERS IS NULL");
		    	return;
		    }
		    PATTERN_PARAMETERS.debugLevel=DEBUG_LEVEL;
		    PATTERN_PARAMETERS.updateStatus=UPDATE_STATUS;
		    
		    GenericDialog gd = new GenericDialog("Reset grid options");
		    gd.addCheckbox("Reset geometry", false);
		    gd.addNumericField("Update pattern views to number of sensors (26 - Eyesis4pi, 0 - do not change)", 0 ,0);
		    gd.addMessage ("OK will at least reset pattern margins/alpha");
		    gd.showDialog();
		    if (gd.wasCanceled()) return;
		    boolean resetGeometry=gd.getNextBoolean();
		    int numSensors= (int) gd.getNextNumber();
		    if (numSensors>0){
		    	PATTERN_PARAMETERS.initDefaultChannels(numSensors);
		    	PATTERN_PARAMETERS.setPhotometric();
		    }
			PATTERN_PARAMETERS.calculateGridGeometryAndPhotometric(resetGeometry);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Reset Margins")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
		    if (PATTERN_PARAMETERS==null){
		    	IJ.showMessage("PATTERN_PARAMETERS IS NULL");
		    	return;
		    }
		    PATTERN_PARAMETERS.debugLevel=DEBUG_LEVEL;
		    PATTERN_PARAMETERS.updateStatus=UPDATE_STATUS;
			PATTERN_PARAMETERS.calculateGridGeometryAndPhotometric(false);
			return;
		}
		
/* ======================================================================== */
		if       (label.equals("Restore Grid")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
		    if (PATTERN_PARAMETERS==null){
		    	IJ.showMessage("PATTERN_PARAMETERS IS NULL");
		    	return;
		    }
		    PATTERN_PARAMETERS.debugLevel=DEBUG_LEVEL;
		    PATTERN_PARAMETERS.updateStatus=UPDATE_STATUS;
			if (DISTORTION_CALIBRATION_DATA==null){
				IJ.showMessage("Distortion Calibration data is not initialized - create it with"+
				"\"SelectGrid Files\" or \"Restore Calibration\"");
				return;
			}
		    PATTERN_PARAMETERS.selectAndRestore(
		    		false,
		    		null,
					DISTORTION_CALIBRATION_DATA.eyesisCameraParameters.numStations);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Correct Grid")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			PATTERN_PARAMETERS.debugLevel=DEBUG_LEVEL;
			if (LENS_DISTORTIONS.fittingStrategy==null){
				IJ.showMessage("Distortion Fitting strategy is not initialized - create it with" +
				"\"New Strategy\", \"Edit Strategy\" or \"Restore Strategy\"");
				return;
			}
			if (DISTORTION_CALIBRATION_DATA==null){
				IJ.showMessage("Distortion Calibration data is not initialized - create it with"+
				"\"SelectGrid Files\" or \"Restore Calibration\"");
				return;
			}
			DISTORTION_CALIBRATION_DATA.debugLevel=DEBUG_LEVEL;
			LENS_DISTORTIONS.fittingStrategy.debugLevel=DEBUG_LEVEL;
			if (DISTORTION_CALIBRATION_DATA.sensorMasks==null){
				IJ.showMessage("Sensor mask(s) are not initialized - create them with"+
				"\"Calculate Sensor Masks\" or \"Restore Sensor Masks\"");
				return;
			}
			LENS_DISTORTIONS.modifyGrid(
					DISTORTION_CALIBRATION_DATA,
					THREADS_MAX,
					UPDATE_STATUS);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Reset Variations")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			PATTERN_PARAMETERS.debugLevel=DEBUG_LEVEL;
			LENS_DISTORTIONS.patternParameters.resetStationZCorr();
			return;
		}
		
/* ======================================================================== */
		if       (label.equals("Correct Grid0")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			PATTERN_PARAMETERS.debugLevel=DEBUG_LEVEL;
			if (LENS_DISTORTIONS.fittingStrategy==null){
				IJ.showMessage("Distortion Fitting strategy is not initialized - create it with" +
				"\"New Strategy\", \"Edit Strategy\" or \"Restore Strategy\"");
				return;
			}
			if (DISTORTION_CALIBRATION_DATA==null){
				IJ.showMessage("Distortion Calibration data is not initialized - create it with"+
				"\"SelectGrid Files\" or \"Restore Calibration\"");
				return;
			}
			DISTORTION_CALIBRATION_DATA.debugLevel=DEBUG_LEVEL;
			LENS_DISTORTIONS.fittingStrategy.debugLevel=DEBUG_LEVEL;
			if (DISTORTION_CALIBRATION_DATA.sensorMasks==null){
				IJ.showMessage("Sensor mask(s) are not initialized - create them with"+
				"\"Calculate Sensor Masks\" or \"Restore Sensor Masks\"");
				return;
			}
			LENS_DISTORTIONS.modifyGrid0(DISTORTION_CALIBRATION_DATA);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Grid Diffs")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			PATTERN_PARAMETERS.debugLevel=DEBUG_LEVEL;
			if (LENS_DISTORTIONS.fittingStrategy==null){
				IJ.showMessage("Distortion Fitting strategy is not initialized - create it with" +
				"\"New Strategy\", \"Edit Strategy\" or \"Restore Strategy\"");
				return;
			}
			if (DISTORTION_CALIBRATION_DATA==null){
				IJ.showMessage("Distortion Calibration data is not initialized - create it with"+
				"\"SelectGrid Files\" or \"Restore Calibration\"");
				return;
			}
			DISTORTION_CALIBRATION_DATA.debugLevel=DEBUG_LEVEL;
			LENS_DISTORTIONS.fittingStrategy.debugLevel=DEBUG_LEVEL;
			/*
			if (DISTORTION_CALIBRATION_DATA.sensorMasks==null){
				IJ.showMessage("Sensor mask(s) are not initialized - create them with"+
				"\"Calculate Sensor Masks\" or \"Restore Sensor Masks\"");
				return;
			}
			*/
			LENS_DISTORTIONS.patternErrors(
					THREADS_MAX,
					UPDATE_STATUS,
					DEBUG_LEVEL);
			return;
		}
		
		
/* ======================================================================== */
		if       (label.equals("Save Grid")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
		    if (PATTERN_PARAMETERS==null){
		    	IJ.showMessage("PATTERN_PARAMETERS IS NULL");
		    	return;
		    }
		    PATTERN_PARAMETERS.debugLevel=DEBUG_LEVEL;
		    PATTERN_PARAMETERS.updateStatus=UPDATE_STATUS;
		    PATTERN_PARAMETERS.selectAndSave(false, null);
			return;
		}
		
/* ======================================================================== */
		if       (label.equals("Reset Sensor")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			LENS_DISTORTIONS.resetSensorCorrection();
			LENS_DISTORTIONS.initSensorCorrection(); // set zero corerctions (to be able to save sesnor correction files)
			return;
		}
/* ======================================================================== */
		if       (label.equals("Restore Sensor")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			String [] extensions={".calib-tiff"};
			CalibrationFileManagement.MultipleExtensionsFileFilter parFilter = new CalibrationFileManagement.MultipleExtensionsFileFilter("",extensions,"distortion calibration .calib-tiff files");
			String pathname=CalibrationFileManagement.selectFile(false,
					"Restore distortion calibration for sensor",
					"Restore",
					parFilter,
					LENS_DISTORTIONS.getSensorPath(-1)); //String defaultPath
			if ((pathname==null) || (pathname=="")) return;
		    GenericDialog gd = new GenericDialog("Select sensor number");
		    gd.addCheckbox("Restore all sensors", true);
			gd.addNumericField("Number of sensor/channel to apply calibration (if \"all\" is not selected)", 0,0);
			gd.addCheckbox("Overwrite SFE position/orientation from the sensor calibration data", true);
		    gd.showDialog();
		    if (gd.wasCanceled()) return;
		    boolean allFiles=gd.getNextBoolean();
		    int numSensor= (int) gd.getNextNumber();
		    if (allFiles) numSensor=-1;
			boolean overwriteExtrinsic=gd.getNextBoolean();
		    if (numSensor<0) LENS_DISTORTIONS.setDistortionFromImageStack(pathname,overwriteExtrinsic); // requires fitting strategy to be set?
		    else LENS_DISTORTIONS.setDistortionFromImageStack(pathname, numSensor,true,overwriteExtrinsic); // report missing files
			return;
		}
/* ======================================================================== */
		if       (label.equals("Correct Sensor Old")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			PATTERN_PARAMETERS.debugLevel=DEBUG_LEVEL;
			if (LENS_DISTORTIONS.fittingStrategy==null){
				IJ.showMessage("Distortion Fitting strategy is not initialized - create it with" +
				"\"New Strategy\", \"Edit Strategy\" or \"Restore Strategy\"");
				return;
			}
			if (DISTORTION_CALIBRATION_DATA==null){
				IJ.showMessage("Distortion Calibration data is not initialized - create it with"+
				"\"SelectGrid Files\" or \"Restore Calibration\"");
				return;
			}
			DISTORTION_CALIBRATION_DATA.debugLevel=DEBUG_LEVEL;
			LENS_DISTORTIONS.fittingStrategy.debugLevel=DEBUG_LEVEL;
			/*
			if (DISTORTION_CALIBRATION_DATA.sensorMasks==null){
				IJ.showMessage("Sensor mask(s) are not initialized - create them with"+
				"\"Calculate Sensor Masks\" or \"Restore Sensor Masks\"");
				return;
			}
			*/
			LENS_DISTORTIONS.modifyPixelCorrection(DISTORTION_CALIBRATION_DATA);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Correct Sensor")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			PATTERN_PARAMETERS.debugLevel=DEBUG_LEVEL;
			if (LENS_DISTORTIONS.fittingStrategy==null){
				IJ.showMessage("Distortion Fitting strategy is not initialized - create it with" +
				"\"New Strategy\", \"Edit Strategy\" or \"Restore Strategy\"");
				return;
			}
			if (DISTORTION_CALIBRATION_DATA==null){
				IJ.showMessage("Distortion Calibration data is not initialized - create it with"+
				"\"SelectGrid Files\" or \"Restore Calibration\"");
				return;
			}
			DISTORTION_CALIBRATION_DATA.debugLevel=DEBUG_LEVEL;
			LENS_DISTORTIONS.fittingStrategy.debugLevel=DEBUG_LEVEL;
	    	if (LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData.eyesisCameraParameters==null){
	    		String msg="Eyesis camera parameters (and sensor dimensions) are not defined";
	    		IJ.showMessage("Error",msg);
	    		throw new IllegalArgumentException (msg);
	    	}
//	    	int series=refineParameters.showDialog("Select Lens Distortion Residual Compensation Parameters", 0x1efff, (this.seriesNumber>=0)?this.seriesNumber:0);
	    	int series=LENS_DISTORTIONS.refineParameters.showDialog(
	    			"Select Lens Distortion Residual Compensation Parameters",
//	    			0x846f1,
	    			0x94ff1,
	    			((LENS_DISTORTIONS.seriesNumber>=0)?LENS_DISTORTIONS.seriesNumber:0),
	    			null); // averageRGB - only for target flat-field correction
	    	if (series<0) return;
	    	LENS_DISTORTIONS.seriesNumber=series;
			LENS_DISTORTIONS.modifyPixelCorrection(
					true, // enable show
					THREADS_MAX,
					UPDATE_STATUS,
					DEBUG_LEVEL);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Save Sensor")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			String [] extensions={".calib-tiff"};
			CalibrationFileManagement.MultipleExtensionsFileFilter parFilter = new CalibrationFileManagement.MultipleExtensionsFileFilter("",extensions,"distortion calibration .calib-tiff files");
			String pathname=CalibrationFileManagement.selectFile(true,
					"Save distortion calibration for sensor (will add channel number when saving all)",
					"Save",
					parFilter,
					LENS_DISTORTIONS.getSensorPath(-1)); //String defaultPath
			if ((pathname==null) || (pathname=="")) return;
		    GenericDialog gd = new GenericDialog("Select sensor number");
		    gd.addCheckbox("Save all sensors", true);
			gd.addNumericField("Number of sensor/channel to save calibration (if \"all\" is not selected)", 0,0);
			gd.addCheckbox("Save non-calibared sensors", false);
		    gd.showDialog();
		    if (gd.wasCanceled()) return;
		    boolean allFiles=gd.getNextBoolean();
		    int numSensor= (int) gd.getNextNumber();
		    if (allFiles) numSensor=-1;
		    boolean saveNonCalibrated=gd.getNextBoolean();
		    if (numSensor<0){
		    	LENS_DISTORTIONS.saveDistortionAsImageStack(
		    			CAMERAS, // to save channel map
		    			"sensor_calibratrion" , //String title,
		    			pathname,
		    			saveNonCalibrated); // boolean emptyOK) if false will throw for non-calibrated sensors 
		    } else {
		    	LENS_DISTORTIONS.saveDistortionAsImageStack(
		    			CAMERAS, // to save channel map
		    			"sensor_calibratrion" , //String title,
		    			pathname,
		    			numSensor,
		    			saveNonCalibrated); // boolean emptyOK) if false will throw for non-calibrated sensors
		    }
			return;
		}
/* ======================================================================== */
/*		if       (label.equals("TestIpl")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			LENS_DISTORTIONS.testExtrapolateSensorCorrection();
			return;
		}
*/		
/* ======================================================================== */
		if       (label.equals("Convert X/Y slices to color")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			imp_sel = WindowManager.getCurrentImage();
			if (imp_sel==null){
				IJ.showMessage("Error","There is no image selected");
				return;
			}
			ImagePlus imp_flow=SDFA_INSTANCE.showFlowFromSlices(imp_sel);
			imp_flow.getProcessor().resetMinAndMax();
			imp_flow.show();

			return;
		}
/* ======================================================================== */
		if       (label.equals("List Calibration")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			if (LENS_DISTORTIONS.fittingStrategy==null) return;
			LENS_DISTORTIONS.fittingStrategy.debugLevel=DEBUG_LEVEL;
			LENS_DISTORTIONS.listImageParameters(false); // not silent
			return;
		}
/* ======================================================================== */
		if       (label.equals("Reproject")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			if (LENS_DISTORTIONS.fittingStrategy==null) return;
			LENS_DISTORTIONS.fittingStrategy.debugLevel=DEBUG_LEVEL;
			LENS_DISTORTIONS.showImageReprojectionErrorsDialog(DEBUG_LEVEL);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Configure Focusing")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (FOCUS_MEASUREMENT_PARAMETERS.showDialog("Focus Measurement Parameters")) return;
			MOTORS.setHysteresis(FOCUS_MEASUREMENT_PARAMETERS.motorHysteresis);
			MOTORS.setCalmMotors(FOCUS_MEASUREMENT_PARAMETERS.motorCalm);
			MOTORS.setDebug(FOCUS_MEASUREMENT_PARAMETERS.motorDebug);

			if (FOCUS_MEASUREMENT_PARAMETERS.configureCamera) {
				if (CAMERAS.showDialog("Configure cameras interface", 1, true)) FOCUS_MEASUREMENT_PARAMETERS.cameraIsConfigured=true;
			}
			return;
		}
/* ======================================================================== */
		if       (label.equals("Reset Focusing")) {
			matchSimulatedPattern=null;
			return;
		}
/* ======================================================================== */
///Get Focusing Grid		
		if       (label.equals("Get Focusing Grid")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (matchSimulatedPattern==null) return;
			matchSimulatedPattern.debugLevel=DEBUG_LEVEL;
			if (DEBUG_LEVEL>0) { // calulate/print number of defined nodes ina grid
				System.out.println("number of defined grid cells (before maskFocus()) = "+matchSimulatedPattern.numDefinedCells());
			}
			if (LENS_DISTORTION_PARAMETERS==null){
				IJ.showMessage("LENS_DISTORTION_PARAMETERS is not set");
				return;
			}
			matchSimulatedPattern.maskFocus(
			       	LENS_DISTORTION_PARAMETERS.px0, // pixel coordinate of the the optical center
			       	LENS_DISTORTION_PARAMETERS.py0, // pixel coordinate of the the optical center
					FOCUS_MEASUREMENT_PARAMETERS);
			if (DEBUG_LEVEL>0) { // calulate/print number of defined nodes ina grid
				System.out.println("number of defined grid cells (after maskFocus()) = "+matchSimulatedPattern.numDefinedCells());
			}
			
			if (DEBUG_LEVEL>0){
				SDFA_INSTANCE.showArrays(
						matchSimulatedPattern.flatFieldForGrid,
						matchSimulatedPattern.getImageWidth(),
						matchSimulatedPattern.getImageHeight(),
						"flatFieldForGrid");
			}
/*			SIM_ARRAY=	simulateGridAll (
					matchSimulatedPattern.getImageWidth(),
					matchSimulatedPattern.getImageHeight(),
					matchSimulatedPattern.getDArray(),
					2, // gridFrac, // number of grid steps per pattern full period
					SIMUL,
					THREADS_MAX,
					UPDATE_STATUS,
					DISTORTION.loop_debug_level); // debug level
*/			
			SIM_ARRAY=	(new SimulationPattern(SIMUL)).simulateGridAll (
					matchSimulatedPattern.getImageWidth(),
					matchSimulatedPattern.getImageHeight(),
					matchSimulatedPattern,
					2, // gridFrac, // number of grid steps per pattern full period
					SIMUL,
					THREADS_MAX,
					UPDATE_STATUS,
					DEBUG_LEVEL,
					DISTORTION.loop_debug_level); // debug level
			
			if (DEBUG_LEVEL>0){
				SDFA_INSTANCE.showArrays(
						SIM_ARRAY,
						matchSimulatedPattern.getImageWidth()*SIMUL.subdiv/2,
						matchSimulatedPattern.getImageHeight()*SIMUL.subdiv/2,
						"focus-simulation");
			}
			return;
		}
/* ======================================================================== */
//		addButton("Update Focusing Grid",panelFocusing);
		if  (label.equals("Update Focusing Grid")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (matchSimulatedPattern==null) {
				matchSimulatedPattern= new MatchSimulatedPattern(DISTORTION.FFTSize); // new instance, all reset
			}
			matchSimulatedPattern.debugLevel=DEBUG_LEVEL;
			imp_sel = WindowManager.getCurrentImage();
			if (imp_sel==null){
				IJ.showMessage("Error","There are no images open\nProcess canceled");
				return;
			}
			if (LENS_DISTORTION_PARAMETERS==null){
				IJ.showMessage("LENS_DISTORTION_PARAMETERS is not set");
				return;
			}
			LENS_ADJUSTMENT.updateFocusGrid(
					LENS_DISTORTION_PARAMETERS.px0, // pixel coordinate of the the optical center
					LENS_DISTORTION_PARAMETERS.py0, // pixel coordinate of the the optical center
					imp_sel,
					matchSimulatedPattern,
					DISTORTION,
					FOCUS_MEASUREMENT_PARAMETERS,
					PATTERN_DETECT,
					null, //LASER_POINTERS
					SIMUL,
					true, // remove non-PSF areas
					COMPONENTS.equalizeGreens,
					THREADS_MAX,
					UPDATE_STATUS,
//					DEBUG_LEVEL);
			DISTORTION.loop_debug_level);
// just to display
			if (DEBUG_LEVEL>1){
/*
				SIM_ARRAY=	simulateGridAll (
						matchSimulatedPattern.getImageWidth(),
						matchSimulatedPattern.getImageHeight(),
						matchSimulatedPattern.getDArray(),
						2, // gridFrac, // number of grid steps per pattern full period
						SIMUL,
						THREADS_MAX,
						UPDATE_STATUS,
						DISTORTION.loop_debug_level); // debug level
*/				
				SIM_ARRAY=	(new SimulationPattern(SIMUL)).simulateGridAll (
						matchSimulatedPattern.getImageWidth(),
						matchSimulatedPattern.getImageHeight(),
						matchSimulatedPattern,
						2, // gridFrac, // number of grid steps per pattern full period
						SIMUL,
						THREADS_MAX,
						UPDATE_STATUS,
						DEBUG_LEVEL,
						DISTORTION.loop_debug_level); // debug level

				
				if (DEBUG_LEVEL>0){
					SDFA_INSTANCE.showArrays(
							SIM_ARRAY,
							matchSimulatedPattern.getImageWidth()*SIMUL.subdiv/2,
							matchSimulatedPattern.getImageHeight()*SIMUL.subdiv/2,
					"focus-simulation");
				}
			}
			return;
		}
		
/* ======================================================================== */
		if       (label.equals("Focusing PSF")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (matchSimulatedPattern==null) return;
			matchSimulatedPattern.debugLevel=DEBUG_LEVEL;
			if (LENS_DISTORTION_PARAMETERS==null){
				IJ.showMessage("LENS_DISTORTION_PARAMETERS is not set");
				return;
			}
			matchSimulatedPattern.maskFocus(
					LENS_DISTORTION_PARAMETERS.px0, // pixel coordinate of the the optical center
					LENS_DISTORTION_PARAMETERS.py0, // pixel coordinate of the the optical center
					FOCUS_MEASUREMENT_PARAMETERS);
			imp_sel = WindowManager.getCurrentImage();
			if (imp_sel==null){
				IJ.showMessage("Error","There are no images open\nProcess canceled");
				return;
			}

			focusPSF(
					LENS_DISTORTION_PARAMETERS.px0, // pixel coordinate of the the optical center
					LENS_DISTORTION_PARAMETERS.py0, // pixel coordinate of the the optical center
					imp_sel,
					matchSimulatedPattern,
					FOCUS_MEASUREMENT_PARAMETERS,
					SIMUL,
					0.0, // double overexposedMaxFraction, ( MULTIFILE_PSF.overexposedMaxFraction, )
					COMPONENTS,
//					PSF_SUBPIXEL,
					OTF_FILTER,
					PSF_PARS, // step of the new map (should be multiple of map step)
					THREADS_MAX,
					UPDATE_STATUS,
					DISTORTION.loop_debug_level); // debug level
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			ImageStack mergedStack=mergeKernelsToStack(PSF_KERNEL_MAP);
			if (mergedStack==null) {
				IJ.showMessage("Error","No PSF kernels to show");
				return;
			}
			ImagePlus imp_psf=SDFA_INSTANCE.showImageStack(mergedStack, imp_sel.getTitle()+"-FOCUS-PSF_KERNEL");
			if (PSF_SAVE_FILE) {
				String dir="";
				if (imp_sel.getOriginalFileInfo()!=null) dir=imp_sel.getOriginalFileInfo().directory;
				String path=dir+"FOCUS-PSF-"+imp_sel.getTitle();
				if (DEBUG_LEVEL>1) {
					System.out.println("Saving result to "+path);
				}
				IJ.saveAs(imp_psf,"tif",path);
			}
			
			
			return;
		}
/* ======================================================================== */
		if      ((label.equals("Focusing New PSF")) || (label.equals("Focusing Acquire PSF"))) {
			long 	  startTime=System.nanoTime();
			if (label.equals("Focusing Acquire PSF")) {
//				FOCUS_MEASUREMENT_PARAMETERS.showDialog("Focus Measurement Parameters");
				if (!FOCUS_MEASUREMENT_PARAMETERS.cameraIsConfigured) {
					if (CAMERAS.showDialog("Configure cameras interface", 1, true)){
						FOCUS_MEASUREMENT_PARAMETERS.cameraIsConfigured=true;
						if (!FOCUS_MEASUREMENT_PARAMETERS.getLensSerial()) return;
//						IJ.showMessage("Notice","Make sure camera is in TRIG=4 mode, JP4, correct exposure/white balance...");
						// reset histories
						MOTORS.clearPreFocus();
						MOTORS.clearHistory();
					} else {
						IJ.showMessage("Error","Camera is not configured\nProcess canceled");
						return;
					}
				}
// acquire camera image here, no lasers.
				FOCUS_MEASUREMENT_PARAMETERS.sensorTemperature=CAMERAS.getSensorTemperature(0,FOCUS_MEASUREMENT_PARAMETERS.EEPROM_channel);
				imp_sel= CAMERAS.acquireSingleImage (
						false, //boolean useLasers,
						UPDATE_STATUS);
				if (imp_sel==null){
					IJ.showMessage("Error","Failed to get camera image\nProcess canceled");
					return;
				}
				if (FOCUS_MEASUREMENT_PARAMETERS.showAcquiredImages){
					imp_sel.show();
					imp_sel.updateAndDraw();
				}
				if (DEBUG_LEVEL>0) System.out.println("Image acquisition (@"+FOCUS_MEASUREMENT_PARAMETERS.sensorTemperature+"C) done at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
			} else {
				imp_sel = WindowManager.getCurrentImage();
				if (imp_sel==null){
					IJ.showMessage("Error","There are no images open\nProcess canceled");
					return;
				}
				
			}
	//"Focusing Acquire PSF"		
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (matchSimulatedPattern==null) {
				matchSimulatedPattern= new MatchSimulatedPattern(DISTORTION.FFTSize); // new instance, all reset
			}
			matchSimulatedPattern.debugLevel=DEBUG_LEVEL;
			if (LENS_DISTORTION_PARAMETERS==null){
				IJ.showMessage("LENS_DISTORTION_PARAMETERS is not set");
				return;
			}
			double [][][][][] rFullResults=new double [1][][][][];
			double [][] metrics=measurePSFMetrics(
					imp_sel,
					LENS_DISTORTION_PARAMETERS,
					matchSimulatedPattern,
					FOCUS_MEASUREMENT_PARAMETERS,
					PATTERN_DETECT,
					DISTORTION,
					SIMUL,
					COMPONENTS,
					OTF_FILTER,
					PSF_PARS,
					rFullResults,
					THREADS_MAX,
					UPDATE_STATUS,
					DEBUG_LEVEL,
					DISTORTION.loop_debug_level);
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL; // Needed!
			if (DEBUG_LEVEL>0){
            	String [] compColors={"green","red","blue","green","green","green","AVERAGE"};
            	for (int c=0;c<metrics.length;c++) if (metrics[c]!=null){
            		System.out.println(compColors[c]+": Far/Near="+IJ.d2s(metrics[c][0],3)+
            				"  Tilt X="+IJ.d2s(metrics[c][1],3)+
            				"  Tilt Y="+IJ.d2s(metrics[c][2],3)+
            				"  R50% average="+IJ.d2s(metrics[c][3],3)+" sensor pixels,"+
            				"  A50% average="+IJ.d2s(metrics[c][4],3)+" sensor pixels,"+
            				"  B50% average="+IJ.d2s(metrics[c][5],3)+" sensor pixels,"+
    					    "  R50%Center="+IJ.d2s(metrics[c][6],3)+" sensor pixels,"+
           				"  component weight="+IJ.d2s(100*metrics[c][7],1)+"%");
            	}
            }
			if (DEBUG_LEVEL>0) System.out.println("Calculating PSF parameters done at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
			return;
		}
/* ======================================================================== */
		if      (label.equals("Select WOI")) {
			imp_sel = WindowManager.getCurrentImage();
			if ((imp_sel==null) || (imp_sel.getRoi()==null)){
				IJ.showMessage("Error","Image with selection is required");
				return;
			}
			FOCUS_MEASUREMENT_PARAMETERS.margins=imp_sel.getRoi().getBounds();
			return;
		}
//		addButton("Head Orientation",panelFocusing);
/* ======================================================================== */
		if      (label.equals("Head Orientation")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			UV_LED_LASERS.debugLevel=DEBUG_LEVEL;
			UV_LED_LASERS.setParameters(FOCUS_MEASUREMENT_PARAMETERS);
			if (!FOCUS_MEASUREMENT_PARAMETERS.cameraIsConfigured) {
				if (CAMERAS.showDialog("Configure cameras interface", 1, true)){
					FOCUS_MEASUREMENT_PARAMETERS.cameraIsConfigured=true;
					if (!FOCUS_MEASUREMENT_PARAMETERS.getLensSerial()) return;
//					IJ.showMessage("Notice","Make sure camera is in TRIG=4 mode, JP4, correct exposure/white balance...");
					// reset histories
					MOTORS.clearPreFocus();
					MOTORS.clearHistory();
				} else {
					IJ.showMessage("Error","Camera is not configured\nProcess canceled");
					return;
				}
			}
//			long 	  startTime=System.nanoTime();
			FOCUS_MEASUREMENT_PARAMETERS.sensorTemperature=CAMERAS.getSensorTemperature(0,FOCUS_MEASUREMENT_PARAMETERS.EEPROM_channel);
			CAMERAS.debugLevel=DEBUG_LEVEL;
			imp_sel= CAMERAS.acquireSingleImage (
					UV_LED_LASERS,
					UPDATE_STATUS);
			if (imp_sel==null){
				IJ.showMessage("Error","Failed to get camera image\nProcess canceled");
				return;
			}
			if (FOCUS_MEASUREMENT_PARAMETERS.showAcquiredImages){
				imp_sel.show();
				imp_sel.updateAndDraw();
			}
			double [][] headPointers=CAMERAS.getHeadPointers(imp_sel);
			if (headPointers==null) {
				System.out.println("Failed to locate optical head laser pointers");
				return;
			}
			if (DEBUG_LEVEL>1){
				for (int n=0;n<headPointers.length;n++) if (headPointers[n]!=null){
					System.out.println("Head pointer "+n+": X="+IJ.d2s(headPointers[n][0],2)+", Y="+IJ.d2s(headPointers[n][1],2));
				}
			}
			double headPointersTilt=Double.NaN;
			if ((headPointers[0]!=null) && (headPointers[1]!=null)){
				headPointersTilt=180.0/Math.PI*Math.atan2(headPointers[1][1]-headPointers[0][1], headPointers[1][0]-headPointers[0][0])-LASER_POINTERS.headLasersTilt;
				if (DEBUG_LEVEL>0){
					System.out.println("SFE is rotated by "+IJ.d2s(headPointersTilt,3)+" degrees according to optical head laser pointers (clockwise positive)");
				}
			}
	   		FOCUS_MEASUREMENT_PARAMETERS.result_ROT=headPointersTilt;
			return;
		}
		
/* ======================================================================== */
//TODO: - reduce length of the laser sequence when focusing?		
		if      ((label.equals("Lens Center")) || (label.equals("Find Grid"))) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			checkSerialAndRestore(); // returns true if did not change or was restored 
			boolean findCenter=label.equals("Lens Center");
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (!FOCUS_MEASUREMENT_PARAMETERS.cameraIsConfigured) {
				if (CAMERAS.showDialog("Configure cameras interface", 1, true)){
					FOCUS_MEASUREMENT_PARAMETERS.cameraIsConfigured=true;
//					IJ.showMessage("Notice","Make sure camera is in TRIG=4 mode, JP4, correct exposure/white balance...");
					if (!FOCUS_MEASUREMENT_PARAMETERS.getLensSerial()) return;
					DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
// reset histories
					MOTORS.clearPreFocus();
					MOTORS.clearHistory();

				} else {
					IJ.showMessage("Error","Camera is not configured\nProcess canceled");
					return;
				}
			}
			long 	  startTime=System.nanoTime();
			FOCUS_MEASUREMENT_PARAMETERS.sensorTemperature=CAMERAS.getSensorTemperature(0,FOCUS_MEASUREMENT_PARAMETERS.EEPROM_channel);
			imp_sel= CAMERAS.acquireSingleImage (
					true, //boolean useLasers,
					UPDATE_STATUS);
			if (imp_sel==null){
				IJ.showMessage("Error","Failed to get camera image\nProcess canceled");
				return;
			}
			if (FOCUS_MEASUREMENT_PARAMETERS.showAcquiredImages){
				imp_sel.show();
				imp_sel.updateAndDraw();
			}
			if (LENS_DISTORTION_PARAMETERS==null){
				IJ.showMessage("LENS_DISTORTION_PARAMETERS is not set");
				return;
			}
			double headPointersTilt=Double.NaN;
			// measure rotation from the optical head lasers
			if (findCenter &&FOCUS_MEASUREMENT_PARAMETERS.useHeadLasers){
				ImagePlus imp_headLasers= CAMERAS.acquireSingleImage (
						UV_LED_LASERS,
						UPDATE_STATUS);
				if (imp_headLasers==null){
					IJ.showMessage("Error","Failed to get camera image\nProcess canceled");
					return;
				}
				if (FOCUS_MEASUREMENT_PARAMETERS.showAcquiredImages){
					imp_headLasers.show();
					imp_headLasers.updateAndDraw();
				}
				double [][] headPointers=CAMERAS.getHeadPointers(imp_headLasers);
				if (headPointers==null) {
					System.out.println("Failed to locate optical head laser pointers");
					return;
				}
				if (DEBUG_LEVEL>1){
					for (int n=0;n<headPointers.length;n++) if (headPointers[n]!=null){
						System.out.println("Head pointer "+n+": X="+IJ.d2s(headPointers[n][0],2)+", Y="+IJ.d2s(headPointers[n][1],2));
					}
				}
				if ((headPointers[0]!=null) && (headPointers[1]!=null)){
					headPointersTilt=180.0/Math.PI*Math.atan2(headPointers[1][1]-headPointers[0][1], headPointers[1][0]-headPointers[0][0])-LASER_POINTERS.headLasersTilt;
					if (DEBUG_LEVEL>0){
						System.out.println("SFE is rotated by "+IJ.d2s(headPointersTilt,3)+" degrees according to optical head laser pointers (clockwise positive)");
					}
				}
				FOCUS_MEASUREMENT_PARAMETERS.result_ROT=headPointersTilt;
			}

			if (DEBUG_LEVEL>0) System.out.println("Image acquisition (@"+FOCUS_MEASUREMENT_PARAMETERS.sensorTemperature+"C) done at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
			// reset matchSimulatedPattern, so it will start from scratch
			matchSimulatedPattern= new MatchSimulatedPattern(DISTORTION.FFTSize); // new instance, all reset
			// next 2 lines are not needed for the new instance, but can be used alternatively if keeipg it
			   matchSimulatedPattern.invalidateFlatFieldForGrid(); //Reset Flat Field calibration - different image. 
			   matchSimulatedPattern.invalidateFocusMask();

			
			if (matchSimulatedPattern.getPointersXY(imp_sel,LASER_POINTERS.laserUVMap.length)==null) { // no pointers in this image
				IJ.showMessage("Error","No laser pointers detected - they are needed for absolute grid positioning\nProcess canceled");
				return;
			}

			matchSimulatedPattern.debugLevel=DEBUG_LEVEL;
			int numAbsolutePoints=LENS_ADJUSTMENT.updateFocusGrid(
					LENS_DISTORTION_PARAMETERS.px0, // pixel coordinate of the the optical center
					LENS_DISTORTION_PARAMETERS.py0, // pixel coordinate of the the optical center
					imp_sel,
					matchSimulatedPattern,
					DISTORTION,
					FOCUS_MEASUREMENT_PARAMETERS,
					PATTERN_DETECT,
					LASER_POINTERS,
					SIMUL,
					false, // keep (not remove) non-PSF areas
					COMPONENTS.equalizeGreens,
					THREADS_MAX,
					UPDATE_STATUS,
//					DEBUG_LEVEL);
			DISTORTION.loop_debug_level);
			if (numAbsolutePoints<=0) { // no pointers in this image
				IJ.showMessage("Error","No laser pointers matched - they are needed for absolute grid positioning\nProcess canceled");
				return;
			}
			if (DEBUG_LEVEL>0) System.out.println("Matched "+numAbsolutePoints+" laser pointers, grid generated at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
			ImagePlus [] imp_calibrated={matchSimulatedPattern.getCalibratedPatternAsImage(imp_sel,numAbsolutePoints)};
			if (FOCUS_MEASUREMENT_PARAMETERS.showAcquiredImages) imp_calibrated[0].show(); // DISTORTION_PROCESS_CONFIGURATION.showGridImages
			if (findCenter){
				// Read required calibration files
				// initial calibration			
				DISTORTION_CALIBRATION_DATA=new Distortions.DistortionCalibrationData(
						true, // skip dialog if file exists
						FOCUS_MEASUREMENT_PARAMETERS.initialCalibrationFile,
						PATTERN_PARAMETERS,
						EYESIS_CAMERA_PARAMETERS, // is it null or 1?
						imp_calibrated); // gridImages null - use specified files - single image
				if (DISTORTION_CALIBRATION_DATA.pathName== null){ // failed to select/open the file
					DISTORTION_CALIBRATION_DATA=null;
					IJ.showMessage("Error","Failed to open initial calibration data file");
					return;
				}
				// fitting strategy			
				FOCUS_MEASUREMENT_PARAMETERS.initialCalibrationFile=DISTORTION_CALIBRATION_DATA.pathName;
				if (LENS_DISTORTIONS==null) {
					if (LENS_DISTORTION_PARAMETERS==null){
						IJ.showMessage("LENS_DISTORTION_PARAMETERS is not set");
						return;
					}
					if (PATTERN_PARAMETERS==null){
						IJ.showMessage("PATTERN_PARAMETERS is not set");
						return;
					}
					LENS_DISTORTIONS=new Distortions(LENS_DISTORTION_PARAMETERS,PATTERN_PARAMETERS,REFINE_PARAMETERS,this.SYNC_COMMAND.stopRequested);
//				} else {
//					LENS_DISTORTION_PARAMETERS=LENS_DISTORTIONS.lensDistortionParameters;  // after working on goniometer LENS_DISTORTION_PARAMETERS stopped updating?
				}
				LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;

				LENS_DISTORTIONS.fittingStrategy=new Distortions.FittingStrategy(
						true,
						FOCUS_MEASUREMENT_PARAMETERS.strategyFile,
						DISTORTION_CALIBRATION_DATA);
				if (LENS_DISTORTIONS.fittingStrategy.pathName== null){ // failed to select/open the file
					LENS_DISTORTIONS.fittingStrategy=null;
					IJ.showMessage("Error","Failed to open fitting strategy file");
					return;
				}
				LENS_DISTORTIONS.fittingStrategy.debugLevel=DEBUG_LEVEL;
				FOCUS_MEASUREMENT_PARAMETERS.strategyFile=LENS_DISTORTIONS.fittingStrategy.pathName;
				// Grid file
				PATTERN_PARAMETERS.debugLevel=DEBUG_LEVEL;
				PATTERN_PARAMETERS.updateStatus=UPDATE_STATUS;
				String gridPathname=PATTERN_PARAMETERS.selectAndRestore(
						true, // skip dialog if not needed
						FOCUS_MEASUREMENT_PARAMETERS.gridGeometryFile,
						DISTORTION_CALIBRATION_DATA.eyesisCameraParameters.numStations);
				if (gridPathname== null){ // failed to select/open the file
					IJ.showMessage("Error","Failed to open grid geometry file");
					return;
				}
				FOCUS_MEASUREMENT_PARAMETERS.gridGeometryFile=gridPathname;

				// Calculate Sensor Masks			
				DISTORTION_CALIBRATION_DATA.debugLevel=DEBUG_LEVEL;
				DISTORTION_CALIBRATION_DATA.updateStatus=UPDATE_STATUS;
				DISTORTION_CALIBRATION_DATA.calculateSensorMasks();

				if (DEBUG_LEVEL>0) System.out.println("Starting LMA at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));

				LENS_DISTORTIONS.seriesNumber=   0; // start from 0;
				LENS_DISTORTIONS.stopEachStep=   false;
				LENS_DISTORTIONS.stopEachSeries= false;

				// TODO: configure through FOCUS_MEASUREMENT_PARAMETERS 			
				LENS_DISTORTIONS.thresholdFinish=FOCUS_MEASUREMENT_PARAMETERS.thresholdFinish;
				LENS_DISTORTIONS.numIterations=FOCUS_MEASUREMENT_PARAMETERS.numIterations;

				LENS_DISTORTIONS.LevenbergMarquardt(false); //  skip dialog
				if (DEBUG_LEVEL>0) System.out.println("Finished LMA at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
				int stationNumber=0;
				// Read camera parameters
				Distortions.EyesisCameraParameters camPars=
					LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData.eyesisCameraParameters;
				double threadPitch=0.35; // M1.6
				double dPx0=camPars.eyesisSubCameras[stationNumber][0].px0-(camPars.sensorWidth/2)-FOCUS_MEASUREMENT_PARAMETERS.centerDeltaX;
				double dPy0=camPars.eyesisSubCameras[stationNumber][0].py0-(camPars.sensorHeight/2)-FOCUS_MEASUREMENT_PARAMETERS.centerDeltaY;
			   	double psi=camPars.eyesisSubCameras[stationNumber][0].psi;     // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
			   	
			   	FOCUS_MEASUREMENT_PARAMETERS.result_PX0=camPars.eyesisSubCameras[stationNumber][0].px0;
			   	FOCUS_MEASUREMENT_PARAMETERS.result_PY0=camPars.eyesisSubCameras[stationNumber][0].py0;
			   	FOCUS_MEASUREMENT_PARAMETERS.result_PSI=-camPars.eyesisSubCameras[stationNumber][0].psi; // "-" to show the current rotation (error), not where to rotate to correct
			   	FOCUS_MEASUREMENT_PARAMETERS.result_FocalLength=camPars.eyesisSubCameras[stationNumber][0].focalLength;
// Use rotation from head lasers
			   	if (FOCUS_MEASUREMENT_PARAMETERS.useHeadLasers){
			   		psi=-headPointersTilt; // "-" - correction, instead of the target tilt
			   	}
			   	
			   	double diffY=12*Math.sin(Math.PI/180.0*psi);
			   	double diffYTurns=diffY/threadPitch;
				double dPx0Mm=dPx0*camPars.eyesisSubCameras[stationNumber][0].pixelSize/1000;
				double dPy0Mm=dPy0*camPars.eyesisSubCameras[stationNumber][0].pixelSize/1000;
				double dPx0Turns=dPx0Mm/threadPitch;
				double dPy0Turns=dPy0Mm/threadPitch;

				double dPy0Mm1=dPy0Mm+diffY;
				double dPy0Mm2=dPy0Mm-diffY;
				double dPy0Turns1=dPy0Turns+diffYTurns;
				double dPy0Turns2=dPy0Turns-diffYTurns;
				String dPYDir1=(dPy0Mm1>0)?"AWAY":"TO";
				String dPYDir2=(dPy0Mm2>0)?"AWAY":"TO";
				String lensCenter=
					(((FOCUS_MEASUREMENT_PARAMETERS.centerDeltaX!=0.0) || (FOCUS_MEASUREMENT_PARAMETERS.centerDeltaX!=0.0))?
							("Required lens center shift "+FOCUS_MEASUREMENT_PARAMETERS.centerDeltaX+"/"+FOCUS_MEASUREMENT_PARAMETERS.centerDeltaY+" pix\n"):"")+
					"Lens center: "+IJ.d2s(camPars.eyesisSubCameras[stationNumber][0].px0,1)+":"+IJ.d2s(camPars.eyesisSubCameras[stationNumber][0].py0,1)+"\n"+ // yes, is updated
					//camPars.eyesisSubCameras[0]
					"dX="+IJ.d2s(dPx0,1)+"pix,  "+IJ.d2s(dPx0Mm,2)+" mm,  "+IJ.d2s(dPx0Turns,2)+" turns."+
					((dPx0>0)?"As it is >0, move sensor LEFT (looking to the target)":" As it is <0, move sensor RIGHT (looking to the target)")+"\n"+
					"dY="+IJ.d2s(dPy0,1)+"pix,  "+IJ.d2s(dPy0Mm,2)+" mm,  "+IJ.d2s(dPy0Turns,2)+" turns."+
					((dPy0>0)?" As it is >0, move sensor away from the target":" As it is <0, move sensor to the target")+"\n"+
					"Sensor rotation (determined from "+(FOCUS_MEASUREMENT_PARAMETERS.useHeadLasers?"optical head lasers":"target tilt")+
					"): "+IJ.d2s(psi,3)+" degrees, compensate by differential of two screws by "+IJ.d2s(diffY,3)+" mm ("+IJ.d2s(diffYTurns,3)+" turns\n"+
					"Combined movement of Y screws (right/left looking to the target): "+IJ.d2s(dPy0Mm1,3)+" "+dPYDir1+"/"+IJ.d2s(dPy0Mm2,3)+" "+dPYDir2+" mm, "+
					IJ.d2s(dPy0Turns1,3)+" "+dPYDir1+"/"+IJ.d2s(dPy0Turns2,3)+" "+dPYDir2+" turns";
			
				if (DEBUG_LEVEL>0) System.out.println(
						"\n\n=====================================================================\n"+
						lensCenter+"\n=====================================================================");
				IJ.showMessage("Instruction",lensCenter);

				if (DEBUG_LEVEL>0) System.out.println("Focusing Distortions finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
				if (FOCUS_MEASUREMENT_PARAMETERS.saveResults) {
					String dir=getResultsPath(FOCUS_MEASUREMENT_PARAMETERS);
					File dFile=new File(dir);
					if (!dFile.isDirectory() &&  !dFile.mkdirs()) {
						String msg="Failed to create directory "+dir;
						IJ.showMessage(msg);
						throw new IllegalArgumentException (msg);
					}
					String lensPrefix="";
					if (FOCUS_MEASUREMENT_PARAMETERS.includeLensSerial && (FOCUS_MEASUREMENT_PARAMETERS.lensSerial.length()>0)){
//						lensPrefix=String.format("LENS%S-",FOCUS_MEASUREMENT_PARAMETERS.lensSerial);
						lensPrefix=String.format("LENS%S-S%02d-",FOCUS_MEASUREMENT_PARAMETERS.lensSerial,FOCUS_MEASUREMENT_PARAMETERS.manufacturingState);
					}
					String path=dFile+Prefs.getFileSeparator()+lensPrefix+CAMERAS.getLastTimestampUnderscored()+".dcal-xml";
					if (MASTER_DEBUG_LEVEL>0) System.out.println ("Saving distortion parameters (including lens center) to "+path);
					DISTORTION_CALIBRATION_DATA.saveToXML(path, FOCUS_MEASUREMENT_PARAMETERS.comment); // just save
					saveCurrentConfig();
				}
				
			} else {
				System.out.println("TODO: add automatic WOI from grid, till then use \"Select WOI\" command");
			}
			return;
		}
/* ======================================================================== */
		if       (label.equals("Reset Histories")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			MOTORS.clearPreFocus();
			MOTORS.clearHistory();
			return;
		}
/* ======================================================================== */
		if       (label.equals("Motors Home")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			MOTORS.setDebug(FOCUS_MEASUREMENT_PARAMETERS.motorDebug);

			int [] positions={0,0,0};
			MOTORS.moveElphel10364Motors(true, positions, 0, true, "moving motors to home position",true);

			MOTORS.clearPreFocus();
			MOTORS.clearHistory();
			return;
		}
		//
/* ======================================================================== */
		if       (label.equals("List Pre-focus")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			MOTORS.listPreFocus(null,FOCUS_MEASUREMENT_PARAMETERS.comment);
			return;
		}
/* ======================================================================== */
		if       ((label.equals("Manual Pre-focus")) || (label.equals("Auto Pre-focus"))) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			checkSerialAndRestore(); // returns true if did not change or was restored 
			MOTORS.setDebug(FOCUS_MEASUREMENT_PARAMETERS.motorDebug);
			boolean autoMove=label.equals("Auto Pre-focus");
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (!FOCUS_MEASUREMENT_PARAMETERS.cameraIsConfigured) {
				if (CAMERAS.showDialog("Configure cameras interface", 1, true)){
					FOCUS_MEASUREMENT_PARAMETERS.cameraIsConfigured=true;
					if (!FOCUS_MEASUREMENT_PARAMETERS.getLensSerial()) return;
//					IJ.showMessage("Notice","Make sure camera is in TRIG=4 mode, JP4, correct exposure/white balance...");
					// reset histories
					MOTORS.clearPreFocus();
					MOTORS.clearHistory();
				} else {
					IJ.showMessage("Error","Camera is not configured\nProcess canceled");
					return;
				}
			}
			if (matchSimulatedPattern==null) {
				matchSimulatedPattern= new MatchSimulatedPattern(DISTORTION.FFTSize); // new instance, all reset
			}
			matchSimulatedPattern.debugLevel=DEBUG_LEVEL;
			if (LENS_DISTORTION_PARAMETERS==null){
				IJ.showMessage("LENS_DISTORTION_PARAMETERS is not set");
				return;
			}
			int mode=autoMove? (FOCUS_MEASUREMENT_PARAMETERS.confirmFirstAuto?2:3):1;
// not to forget to reset history after forcusing lens by thread 			
			if (autoMove) MOTORS.clearPreFocus();

			int [] newPos=null;
			boolean isAdjusted=false;
			while (true) {
				if ((mode==2)|| (mode==3)){
					newPos=preFocusingStepsAuto(
							(mode==2)?0:FOCUS_MEASUREMENT_PARAMETERS.maxAutoIterations, // maximal number of iterations (0 - suggest only, do not move). When calling from the button - first time single iteration, second time - as specified
									MOTORS,
									CAMERAS,
									LENS_DISTORTION_PARAMETERS,
									matchSimulatedPattern, // should not bee null
									FOCUS_MEASUREMENT_PARAMETERS,
									UPDATE_STATUS,
									DEBUG_LEVEL);
					isAdjusted=((mode==3) && (newPos!=null));
					
					if ((mode==2) && (newPos==null)){
						String message="Failed to suggest automatic adjustment. Exiting";
						System.out.println(message);
						IJ.showMessage(message);
						break;
					}
				} else {
					newPos=null; // will use current position - asigned below
					
				}
				if (newPos==null) newPos=MOTORS.readElphel10364Motors().clone();
				int result=MOTORS.dialogFocusingSharpness(
						isAdjusted,
						newPos, // null is not OK (as it may be be modified!)
		    			mode, // 1 - manual, 2 auto first, 3 - auto normal
		    			FOCUS_MEASUREMENT_PARAMETERS);
//			    System.out.println("*** mode="+mode+" result="+result);
				if (result<=0) break; // -1 - cancel, 0 - OK after adjustment - do nothing
				if ((result==1) || (result==3)) {
					moveAndMeasureSharpness(
							newPos, // move and measure
							MOTORS,
							CAMERAS,
							LENS_DISTORTION_PARAMETERS,
							matchSimulatedPattern, // should not bee null
							FOCUS_MEASUREMENT_PARAMETERS,
							UPDATE_STATUS,
							DEBUG_LEVEL);
				}
				if (result==1) {
					mode = 1; // manual
				} else if (result==2) {
					mode = 2; // auto with confirmation
				} else if (result==3) {
					mode = 3;
				} else if (result==4) {
					mode = 3;
				}
//			    System.out.println("*** new mode="+mode);

			} // while (true), with break;
			if (FOCUS_MEASUREMENT_PARAMETERS.saveResults) {
				String dir=getResultsPath(FOCUS_MEASUREMENT_PARAMETERS);
				File dFile=new File(dir);
				if (!dFile.isDirectory() &&  !dFile.mkdirs()) {
					String msg="Failed to create directory "+dir;
					IJ.showMessage(msg);
					throw new IllegalArgumentException (msg);
				}
				String lensPrefix="";
				if (FOCUS_MEASUREMENT_PARAMETERS.includeLensSerial && (FOCUS_MEASUREMENT_PARAMETERS.lensSerial.length()>0)){
//					lensPrefix=String.format("LENS%S-",FOCUS_MEASUREMENT_PARAMETERS.lensSerial);
					lensPrefix=String.format("LENS%S-S%02d-",FOCUS_MEASUREMENT_PARAMETERS.lensSerial,FOCUS_MEASUREMENT_PARAMETERS.manufacturingState);
				}
				String path=dFile+Prefs.getFileSeparator()+lensPrefix+CAMERAS.getLastTimestampUnderscored()+"-prefocus.csv";
				if (MASTER_DEBUG_LEVEL>0) System.out.println ("Saving pre-focusing log data to "+path);
				MOTORS.listPreFocus(path,FOCUS_MEASUREMENT_PARAMETERS.comment);
				saveCurrentConfig();
			}

			return;
		}
/* ======================================================================== */
		if       (label.equals("Scan Calib")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			checkSerialAndRestore(); // returns true if did not change or was restored
			if (FOCUS_MEASUREMENT_PARAMETERS.showScanningSetup("Setup scanning parameters")) return;
			MOTORS.setHysteresis(FOCUS_MEASUREMENT_PARAMETERS.motorHysteresis);
			MOTORS.setDebug(FOCUS_MEASUREMENT_PARAMETERS.motorDebug);
//			int historyFrom=MOTORS.historySize()-1; // before scanning
			int historyFrom=MOTORS.historySize();  // first during scanning (not yet exist)
			double[] range= ScanFocus(
					null, // center at current position
					MOTORS,
					CAMERAS,
					LENS_DISTORTION_PARAMETERS,
					matchSimulatedPattern, // should not bee null
					FOCUS_MEASUREMENT_PARAMETERS,
					PATTERN_DETECT,
					DISTORTION,
					SIMUL,
					COMPONENTS,
					OTF_FILTER,
					PSF_PARS,
					THREADS_MAX,
					UPDATE_STATUS,
					MASTER_DEBUG_LEVEL,
					DISTORTION.loop_debug_level);
			if (range==null ){
				String msg="Scanning failed";
				System.out.println(msg);
				IJ.showMessage(msg);
				return;
			}
			int historyTo=MOTORS.historySize()-2; // skip last movement 
			if (FOCUS_MEASUREMENT_PARAMETERS.scanHysteresis){
				historyTo-=FOCUS_MEASUREMENT_PARAMETERS.scanHysteresisNumber; // exclude hysteresis measurement
			}
			MOTORS.setLinearReductionRatio(FOCUS_MEASUREMENT_PARAMETERS.linearReductionRatio);
			MOTORS.focusingHistory.setCalibrationHistoryFromTo(historyFrom, historyTo);
			if (FOCUS_MEASUREMENT_PARAMETERS.scanHysteresis) {
				MOTORS.focusingHistory.setHysteresisHistoryFromTo(historyTo+1,MOTORS.historySize()-2);
			}
			String path=null;
			if (FOCUS_MEASUREMENT_PARAMETERS.saveResults) {
				String dir=getResultsPath(FOCUS_MEASUREMENT_PARAMETERS);
				File dFile=new File(dir);
				if (!dFile.isDirectory() &&  !dFile.mkdirs()) {
					String msg="Failed to create directory "+dir;
					IJ.showMessage(msg);
					throw new IllegalArgumentException (msg);
				}
				
				String lensPrefix="";
				if (FOCUS_MEASUREMENT_PARAMETERS.includeLensSerial && (FOCUS_MEASUREMENT_PARAMETERS.lensSerial.length()>0)){
//					lensPrefix=String.format("LENS%S-",FOCUS_MEASUREMENT_PARAMETERS.lensSerial);
					lensPrefix=String.format("LENS%S-S%02d-",FOCUS_MEASUREMENT_PARAMETERS.lensSerial,FOCUS_MEASUREMENT_PARAMETERS.manufacturingState);
				}
				path=dFile+Prefs.getFileSeparator()+lensPrefix+CAMERAS.getLastTimestampUnderscored()+"-scan-calib.csv";
			}			

			int [] newPos=MOTORS.focusingHistory.calibrateLensDistance(
					path,
					FOCUS_MEASUREMENT_PARAMETERS,
	    			true, //centerOnly,
	    			MOTORS.getMicronsPerStep(), //double micronsPerStep,
	    			DEBUG_LEVEL); //  int debugLevel
			if (FOCUS_MEASUREMENT_PARAMETERS.scanHysteresis) {
				FOCUS_MEASUREMENT_PARAMETERS.measuredHysteresis=MOTORS.focusingHistory.getMeasuredHysteresis();
				// report measured and configured hysteresis, suggest to modify and/or check hardware.
				if (MASTER_DEBUG_LEVEL>0) System.out.println("Measured hysteresis is "+
						IJ.d2s(FOCUS_MEASUREMENT_PARAMETERS.measuredHysteresis,0)+" steps, configured hysteresis is set to "+
						FOCUS_MEASUREMENT_PARAMETERS.motorHysteresis+" motor steps");
			}
			if (FOCUS_MEASUREMENT_PARAMETERS.lensDistanceMoveToGoal && (newPos!=null)){
				moveAndMaybeProbe(
						true, // just move, not probe
						newPos, // position maybe updated by the dialog
						MOTORS,
						CAMERAS,
						LENS_DISTORTION_PARAMETERS,
						matchSimulatedPattern, // should not bee null
						FOCUS_MEASUREMENT_PARAMETERS,
						PATTERN_DETECT,
						DISTORTION,
						SIMUL,
						COMPONENTS,
						OTF_FILTER,
						PSF_PARS,
						THREADS_MAX,
						UPDATE_STATUS,
						MASTER_DEBUG_LEVEL,
						DISTORTION.loop_debug_level);
			}
			saveCurrentConfig();
			return;
		}
/* ======================================================================== */
		if       (label.equals("Scan Calib LMA")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			checkSerialAndRestore(); // returns true if did not change or was restored
			if (!FOCUS_MEASUREMENT_PARAMETERS.showScanningSetup("Setup scanning parameters for LMA")) return;
			MOTORS.setHysteresis(FOCUS_MEASUREMENT_PARAMETERS.motorHysteresis);
			MOTORS.setDebug(FOCUS_MEASUREMENT_PARAMETERS.motorDebug);
			double[] range= ScanFocusTilt(
					null, // center at current position
					MOTORS,
					CAMERAS,
					LENS_DISTORTION_PARAMETERS,
					matchSimulatedPattern, // should not bee null
					FOCUS_MEASUREMENT_PARAMETERS,
					PATTERN_DETECT,
					DISTORTION,
					SIMUL,
					COMPONENTS,
					OTF_FILTER,
					PSF_PARS,
					THREADS_MAX,
					UPDATE_STATUS,
					MASTER_DEBUG_LEVEL,
					DISTORTION.loop_debug_level);
			if (range==null ){
				String msg="Scanning failed";
				System.out.println(msg);
				IJ.showMessage(msg);
				return;
			}

			double pX0=FOCUS_MEASUREMENT_PARAMETERS.result_PX0;
			double pY0=FOCUS_MEASUREMENT_PARAMETERS.result_PY0;
			double [][][] sampleCoord=FOCUS_MEASUREMENT_PARAMETERS.sampleCoordinates( //{x,y,r}
					pX0,   // lens center on the sensor
					pY0);
			// set file path
			String path=null;
			String dir=getResultsPath(FOCUS_MEASUREMENT_PARAMETERS);
			File dFile=new File(dir);
			if (!dFile.isDirectory() &&  !dFile.mkdirs()) {
				String msg="Failed to create directory "+dir;
				IJ.showMessage(msg);
				throw new IllegalArgumentException (msg);
			}

			String lensPrefix="";
			if (FOCUS_MEASUREMENT_PARAMETERS.includeLensSerial && (FOCUS_MEASUREMENT_PARAMETERS.lensSerial.length()>0)){
				lensPrefix=String.format("LENS%S-S%02d-",FOCUS_MEASUREMENT_PARAMETERS.lensSerial,FOCUS_MEASUREMENT_PARAMETERS.manufacturingState);
			}
			path=dFile+Prefs.getFileSeparator()+lensPrefix+CAMERAS.getLastTimestampUnderscored()+".history-xml";
			FOCUSING_FIELD= new FocusingField(
					EYESIS_CAMERA_PARAMETERS.getSensorWidth(),
					EYESIS_CAMERA_PARAMETERS.getSensorHeight(),
					0.001*EYESIS_CAMERA_PARAMETERS.getPixelSize(0), //subCamera_0.pixelSize,
					FOCUS_MEASUREMENT_PARAMETERS.serialNumber,
					FOCUS_MEASUREMENT_PARAMETERS.lensSerial, // String lensSerial, // if null - do not add average
					FOCUS_MEASUREMENT_PARAMETERS.comment, // String comment,
					pX0,
					pY0,
					sampleCoord,
					this.SYNC_COMMAND.stopRequested);
			FOCUSING_FIELD.setDebugLevel(DEBUG_LEVEL);
			FOCUSING_FIELD.setAdjustMode(false);
			if (PROPERTIES!=null) FOCUSING_FIELD.getProperties("FOCUSING_FIELD.", PROPERTIES, true);
			MOTORS.addCurrentHistoryToFocusingField(FOCUSING_FIELD);
			System.out.println("Saving measurement history to "+path);
			FOCUSING_FIELD.saveXML(path);
			saveCurrentConfig();
// for now just copying from "Restore History". TODO: Make both more automatic (move number of parameters outside?)			
			if (!FOCUSING_FIELD.configureDataVector(
					true, //boolean silent
					"Configure curvature - TODO: fix many settings restored from properties", // String title,
					true, //boolean forcenew
					true) // boolean enableReset
					) return;
			System.out.println("TODO: fix many settings restored from properties, overwriting entered fields. Currently run \"Modify LMA\" to re-enter values");
			System.out.println("TODO: Probably need to make a separate dialog that enters number of parameters.");
	    	double [] sv=          FOCUSING_FIELD.fieldFitting.createParameterVector(FOCUSING_FIELD.sagittalMaster);
			FOCUSING_FIELD.setDataVector(
					true, // calibrate mode
					FOCUSING_FIELD.createDataVector());
			double [] focusing_fx= FOCUSING_FIELD.createFXandJacobian(sv, false);
			double rms=            FOCUSING_FIELD.calcErrorDiffY(focusing_fx, false);
			double rms_pure=       FOCUSING_FIELD.calcErrorDiffY(focusing_fx, true);
			System.out.println("rms="+rms+", rms_pure="+rms_pure+" - with old parameters may be well off.");
			return;
		}
		
		
		
/* ======================================================================== */
		if       ((label.equals("Manual Focus/Tilt")) || (label.equals("Auto Focus/Tilt"))|| (label.equals("Fine Focus"))) {
			checkSerialAndRestore(); // returns true if did not change or was restored 
			MOTORS.setDebug(FOCUS_MEASUREMENT_PARAMETERS.motorDebug);
			// will try to apply absolute calibration to the center lens distance
//			int [] mmm=	 MOTORS.focusingHistory.getPosition();

//System.out.println("@@"+MOTORS.historySize()+": "+MOTORS.curpos[0]+", "+MOTORS.curpos[1]+", "+MOTORS.curpos[2]+" --- "+mmm[0]+", "+mmm[1]+", "+mmm[2]);

			MOTORS.focusingHistory.optimalMotorPosition( // recalculate calibration to estimate current distamnce from center PSF
					FOCUS_MEASUREMENT_PARAMETERS,
	    			MOTORS.getMicronsPerStep(), //double micronsPerStep,
	    			DEBUG_LEVEL);
//System.out.println("@@@"+MOTORS.historySize()+": "+MOTORS.curpos[0]+", "+MOTORS.curpos[1]+", "+MOTORS.curpos[2]+" --- "+mmm[0]+", "+mmm[1]+", "+mmm[2]);
			boolean autoMove= label.equals("Auto Focus/Tilt");
			boolean fineFocus=label.equals("Fine Focus");
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (!FOCUS_MEASUREMENT_PARAMETERS.cameraIsConfigured) {
				if (CAMERAS.showDialog("Configure cameras interface", 1, true)){
					FOCUS_MEASUREMENT_PARAMETERS.cameraIsConfigured=true;
					if (!FOCUS_MEASUREMENT_PARAMETERS.getLensSerial()) return;
//					IJ.showMessage("Notice","Make sure camera is in TRIG=4 mode, JP4, correct exposure/white balance...");
					// reset histories
					MOTORS.clearPreFocus();
					MOTORS.clearHistory();

				} else {
					IJ.showMessage("Error","Camera is not configured\nProcess canceled");
					return;
				}
			}
			if (matchSimulatedPattern==null) {
				matchSimulatedPattern= new MatchSimulatedPattern(DISTORTION.FFTSize); // new instance, all reset
			}
			matchSimulatedPattern.debugLevel=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTION_PARAMETERS==null){
				IJ.showMessage("LENS_DISTORTION_PARAMETERS is not set");
				return;
			}
			int mode=autoMove? (FOCUS_MEASUREMENT_PARAMETERS.confirmFirstAuto?2:3):1;
			if (fineFocus) mode=FOCUS_MEASUREMENT_PARAMETERS.confirmFirstAuto?4:5;
			int [] newPos=null;
			boolean isAdjusted=false;
			int result=-1;
			while (true) {
				if ((mode==2)|| (mode==3)){
					newPos=focusingStepsAuto(
							(mode==2)?0:FOCUS_MEASUREMENT_PARAMETERS.maxAutoIterations, // maximal number of iterations (0 - suggest only, do not move). When calling from the button - first time single iteration, second time - as specified
									MOTORS,
									CAMERAS,
									LENS_DISTORTION_PARAMETERS,
									matchSimulatedPattern, // should not bee null
									FOCUS_MEASUREMENT_PARAMETERS,
									PATTERN_DETECT,
									DISTORTION,
									SIMUL,
									COMPONENTS,
									OTF_FILTER,
									PSF_PARS,
									THREADS_MAX,
									UPDATE_STATUS,
									MASTER_DEBUG_LEVEL,
									DISTORTION.loop_debug_level);
					if ((mode==2) && (newPos==null)){
						String message="Failed to suggest automatic adjustment. Exiting";
						System.out.println(message);
						IJ.showMessage(message);
						result=-1;
						break;
					}
					isAdjusted=((mode==3) && (newPos!=null));
				} else	if ((mode==4)|| (mode==5)){
					newPos=fineFocusingStepsAuto(
							(mode==4)?0:FOCUS_MEASUREMENT_PARAMETERS.maxAutoIterations, // maximal number of iterations (0 - suggest only, do not move). When calling from the button - first time single iteration, second time - as specified
									MOTORS,
									CAMERAS,
									LENS_DISTORTION_PARAMETERS,
									matchSimulatedPattern, // should not bee null
									FOCUS_MEASUREMENT_PARAMETERS,
									PATTERN_DETECT,
									DISTORTION,
									SIMUL,
									COMPONENTS,
									OTF_FILTER,
									PSF_PARS,
									THREADS_MAX,
									UPDATE_STATUS,
									MASTER_DEBUG_LEVEL,
									DISTORTION.loop_debug_level);
					if ((mode==4) && (newPos==null)){
						String message="Failed to suggest automatic fine focus adjustment. Exiting";
						System.out.println(message);
						IJ.showMessage(message);
						result=-1;
						break;
					}
	
				} else {
					newPos=null; // will be changed to current position
				}
				if (newPos==null) newPos=MOTORS.readElphel10364Motors().clone();
				result=MOTORS.dialogFocusing(
						isAdjusted,
						newPos, // null not OK - not here, it can be modified for the next move
		    			mode, // 1 - manual, 2 auto first, 3 - auto normal, 4 fine focus (first), 5 fine focus (auto)
		    			FOCUS_MEASUREMENT_PARAMETERS);
//			    System.out.println("*** mode="+mode+" result="+result);
				if (result<=0) break; //-1 - canceled, 0 - finished
				if ((result==1) || (result==3)  || (result==5) || (result==7)) {
					moveAndMaybeProbe(
							true, // just move, not probe
							newPos, // position maybe updated by the dialog
							MOTORS,
							CAMERAS,
							LENS_DISTORTION_PARAMETERS,
							matchSimulatedPattern, // should not bee null
							FOCUS_MEASUREMENT_PARAMETERS,
							PATTERN_DETECT,
							DISTORTION,
							SIMUL,
							COMPONENTS,
							OTF_FILTER,
							PSF_PARS,
							THREADS_MAX,
							UPDATE_STATUS,
							MASTER_DEBUG_LEVEL,
							DISTORTION.loop_debug_level);
				}
				if (result==1) {
					mode = 1; // manual
				} else if (result==2) {
					mode = 2; // auto with confirmation
				} else if (result==3) {
					mode = 3;
				} else if (result==4) {
					mode = 3;
				} else if (result==6) {
					mode = 4;
				} else if (result==7) {
					mode = 5;
				} else if (result==5) {
					int [] scanPos=newPos.clone();
					boolean scanOK= (ScanFocus(
							scanPos, // position maybe updated by the dialog
							MOTORS,
							CAMERAS,
							LENS_DISTORTION_PARAMETERS,
							matchSimulatedPattern, // should not bee null
							FOCUS_MEASUREMENT_PARAMETERS,
							PATTERN_DETECT,
							DISTORTION,
							SIMUL,
							COMPONENTS,
							OTF_FILTER,
							PSF_PARS,
							THREADS_MAX,
							UPDATE_STATUS,
							MASTER_DEBUG_LEVEL,
							DISTORTION.loop_debug_level)!=null);
					mode = 	1;
					if (!scanOK ){
						String msg="Scanning failed";
						System.out.println(msg);
						IJ.showMessage(msg);
						result=-1;
						break;
					}
				}
			}
			if (FOCUS_MEASUREMENT_PARAMETERS.saveResults) {
				String dir=getResultsPath(FOCUS_MEASUREMENT_PARAMETERS);
				File dFile=new File(dir);
				if (!dFile.isDirectory() &&  !dFile.mkdirs()) {
					String msg="Failed to create directory "+dir;
					IJ.showMessage(msg);
					throw new IllegalArgumentException (msg);
				}
				String lensPrefix="";
				if (FOCUS_MEASUREMENT_PARAMETERS.includeLensSerial && (FOCUS_MEASUREMENT_PARAMETERS.lensSerial.length()>0)){
//					lensPrefix=String.format("LENS%S-",FOCUS_MEASUREMENT_PARAMETERS.lensSerial);
					lensPrefix=String.format("LENS%S-S%02d-",FOCUS_MEASUREMENT_PARAMETERS.lensSerial,FOCUS_MEASUREMENT_PARAMETERS.manufacturingState);
				}
				String path=dFile+Prefs.getFileSeparator()+lensPrefix+CAMERAS.getLastTimestampUnderscored()+"-focus.csv";
				if (MASTER_DEBUG_LEVEL>0) System.out.println ("Saving focusing log data to "+path);
				MOTORS.listHistory(
						path, // on screen, path - to csv
						FOCUS_MEASUREMENT_PARAMETERS.lensSerial,
						FOCUS_MEASUREMENT_PARAMETERS.comment,
						FOCUS_MEASUREMENT_PARAMETERS.showHistoryDetails,
						FOCUS_MEASUREMENT_PARAMETERS.showHistorySamples,
						FOCUS_MEASUREMENT_PARAMETERS.weightRatioRedToGreen,
						FOCUS_MEASUREMENT_PARAMETERS.weightRatioBlueToGreen,
						FOCUS_MEASUREMENT_PARAMETERS.lensDistanceWeightK,
						FOCUS_MEASUREMENT_PARAMETERS.lensDistanceWeightY
						);
			}
			saveCurrentConfig();
			return;
		}
/* ======================================================================== */
		if       (label.equals("Probe around")) {
			checkSerialAndRestore(); // returns true if did not change or was restored 
			MOTORS.setDebug(FOCUS_MEASUREMENT_PARAMETERS.motorDebug);
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (!FOCUS_MEASUREMENT_PARAMETERS.cameraIsConfigured) {
				if (CAMERAS.showDialog("Configure cameras interface", 1, true)){
					FOCUS_MEASUREMENT_PARAMETERS.cameraIsConfigured=true;
					if (!FOCUS_MEASUREMENT_PARAMETERS.getLensSerial()) return;
//					IJ.showMessage("Notice","Make sure camera is in TRIG=4 mode, JP4, correct exposure/white balance...");
					// reset histories
					MOTORS.clearPreFocus();
					MOTORS.clearHistory();

				} else {
					IJ.showMessage("Error","Camera is not configured\nProcess canceled");
					return;
				}
			}
			if (matchSimulatedPattern==null) {
				matchSimulatedPattern= new MatchSimulatedPattern(DISTORTION.FFTSize); // new instance, all reset
			}
			matchSimulatedPattern.debugLevel=DEBUG_LEVEL;
			if (LENS_DISTORTION_PARAMETERS==null){
				IJ.showMessage("LENS_DISTORTION_PARAMETERS is not set");
				return;
			}
			MOTORS.focusingHistory.optimalMotorPosition( // recalculate calibration to estimate current distamnce from center PSF
					FOCUS_MEASUREMENT_PARAMETERS,
	    			MOTORS.getMicronsPerStep(), //double micronsPerStep,
	    			DEBUG_LEVEL);

			moveAndMaybeProbe(
					false,
					MOTORS.readElphel10364Motors().clone(), // null OK
					MOTORS,
					CAMERAS,
					LENS_DISTORTION_PARAMETERS,
					matchSimulatedPattern, // should not be null
					FOCUS_MEASUREMENT_PARAMETERS,
					PATTERN_DETECT,
					DISTORTION,
					SIMUL,
					COMPONENTS,
					OTF_FILTER,
					PSF_PARS,
					THREADS_MAX,
					UPDATE_STATUS,
					MASTER_DEBUG_LEVEL+1,
					DISTORTION.loop_debug_level);
			return;
		}
/* ======================================================================== */
		if       (label.equals("List History")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			MOTORS.listHistory(
					null, // on screen, path - to csv
					FOCUS_MEASUREMENT_PARAMETERS.lensSerial,
					FOCUS_MEASUREMENT_PARAMETERS.comment,
					FOCUS_MEASUREMENT_PARAMETERS.showHistoryDetails,
					FOCUS_MEASUREMENT_PARAMETERS.showHistorySamples,
					FOCUS_MEASUREMENT_PARAMETERS.weightRatioRedToGreen,
					FOCUS_MEASUREMENT_PARAMETERS.weightRatioBlueToGreen,
					FOCUS_MEASUREMENT_PARAMETERS.lensDistanceWeightK,
					FOCUS_MEASUREMENT_PARAMETERS.lensDistanceWeightY
					);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Save History")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			double pX0=FOCUS_MEASUREMENT_PARAMETERS.result_PX0;
			double pY0=FOCUS_MEASUREMENT_PARAMETERS.result_PY0;
			double [][][] sampleCoord=FOCUS_MEASUREMENT_PARAMETERS.sampleCoordinates( //{x,y,r}
					pX0,   // lens center on the sensor
					pY0);
			// set file path
			String path=null;
			String dir=getResultsPath(FOCUS_MEASUREMENT_PARAMETERS);
			File dFile=new File(dir);
			if (!dFile.isDirectory() &&  !dFile.mkdirs()) {
				String msg="Failed to create directory "+dir;
				IJ.showMessage(msg);
				throw new IllegalArgumentException (msg);
			}

			String lensPrefix="";
			if (FOCUS_MEASUREMENT_PARAMETERS.includeLensSerial && (FOCUS_MEASUREMENT_PARAMETERS.lensSerial.length()>0)){
				lensPrefix=String.format("LENS%S-S%02d-",FOCUS_MEASUREMENT_PARAMETERS.lensSerial,FOCUS_MEASUREMENT_PARAMETERS.manufacturingState);
			}
			path=dFile+Prefs.getFileSeparator()+lensPrefix+CAMERAS.getLastTimestampUnderscored()+".history-xml";
			FOCUSING_FIELD= new FocusingField(
					EYESIS_CAMERA_PARAMETERS.getSensorWidth(),
					EYESIS_CAMERA_PARAMETERS.getSensorHeight(),
					0.001*EYESIS_CAMERA_PARAMETERS.getPixelSize(0), //subCamera_0.pixelSize,
					FOCUS_MEASUREMENT_PARAMETERS.serialNumber,
					FOCUS_MEASUREMENT_PARAMETERS.lensSerial, // String lensSerial, // if null - do not add average
					FOCUS_MEASUREMENT_PARAMETERS.comment, // String comment,
					pX0,
					pY0,
					sampleCoord,
					this.SYNC_COMMAND.stopRequested);

			System.out.println("Saving measurement history to "+path);
			MOTORS.addCurrentHistoryToFocusingField(FOCUSING_FIELD);
			FOCUSING_FIELD.saveXML(path);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Restore History")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			FOCUSING_FIELD=new FocusingField(
					true, // boolean smart,       // do not open dialog if default matches 
					"",//); //String defaultPath); //			AtomicInteger stopRequested
					this.SYNC_COMMAND.stopRequested);
			FOCUSING_FIELD.setDebugLevel(DEBUG_LEVEL);
			FOCUSING_FIELD.setAdjustMode(false);
			if (PROPERTIES!=null) FOCUSING_FIELD.getProperties("FOCUSING_FIELD.", PROPERTIES,true); // keep distortions center from history
			System.out.println("Loaded FocusingField");
			if (!FOCUSING_FIELD.configureDataVector(
					true, // boolean silent (maybe add option with false to change number of parameters?)
					"Configure curvature - TODO: fix many settings restored from properties", // String title
					true, // boolean forcenew,
					true) // boolean enableReset
					) return;
			System.out.println("TODO: fix many settings restored from properties, overwriting entered fields. Currently run \"Modify LMA\" to re-enter values");
			System.out.println("TODO: Probably need to make a separate dialog that enters number of parameters.");
	    	double [] sv=          FOCUSING_FIELD.fieldFitting.createParameterVector(FOCUSING_FIELD.sagittalMaster);
			FOCUSING_FIELD.setDataVector(
					true, // calibrate mode
					FOCUSING_FIELD.createDataVector());
			double [] focusing_fx= FOCUSING_FIELD.createFXandJacobian(sv, false);
			double rms=            FOCUSING_FIELD.calcErrorDiffY(focusing_fx, false);
			double rms_pure=       FOCUSING_FIELD.calcErrorDiffY(focusing_fx, true);
			System.out.println("rms="+rms+", rms_pure="+rms_pure);
			return;
		}
/* ======================================================================== */
		if       (label.equals("History RMS")) {
			if (FOCUSING_FIELD==null) return;
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			FOCUSING_FIELD.setDebugLevel(DEBUG_LEVEL);
//			FOCUSING_FIELD.setAdjustMode(false);
	    	double [] sv=          FOCUSING_FIELD.fieldFitting.createParameterVector(FOCUSING_FIELD.sagittalMaster);
			FOCUSING_FIELD.setDataVector(
					true, // calibrate mode
					FOCUSING_FIELD.createDataVector());
			double [] focusing_fx= FOCUSING_FIELD.createFXandJacobian(sv, false);
			double rms=            FOCUSING_FIELD.calcErrorDiffY(focusing_fx, false);
			double rms_pure=       FOCUSING_FIELD.calcErrorDiffY(focusing_fx, true);
			System.out.println("rms="+rms+", rms_pure="+rms_pure);
			FOCUSING_FIELD.printSetRMS(focusing_fx);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Modify LMA")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (FOCUSING_FIELD==null) return;
			FOCUSING_FIELD.setDebugLevel(DEBUG_LEVEL);
			FOCUSING_FIELD.setAdjustMode(false);
			if (!FOCUSING_FIELD.configureDataVector(
					false, // boolean silent,
					"Re-configure curvature parameters", // String title
					false, // boolean forcenew
					true)  // boolean enableReset
					) return;
			FOCUSING_FIELD.setDataVector(
					true, // calibrate mode
					FOCUSING_FIELD.createDataVector());
			return;
		}
/* ======================================================================== */
		if       (label.equals("Load strategies")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (FOCUSING_FIELD==null) return;
			FOCUSING_FIELD.setDebugLevel(DEBUG_LEVEL);
			FOCUSING_FIELD.fieldFitting.fieldStrategies.loadStrategies(null,PROCESS_PARAMETERS.kernelsDirectory);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Save strategies")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (FOCUSING_FIELD==null) return;
			FOCUSING_FIELD.setDebugLevel(DEBUG_LEVEL);
			FOCUSING_FIELD.fieldFitting.fieldStrategies.saveStrategies(null,PROCESS_PARAMETERS.kernelsDirectory);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Organize strategies")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (FOCUSING_FIELD==null) return;
			FOCUSING_FIELD.setDebugLevel(DEBUG_LEVEL);
			int resp=0;
			while (resp==0){
				resp=FOCUSING_FIELD.organizeStrategies("Organize LMA strategies");
			}
			return;
		}
///organizeStrategies(String title)		
/* ======================================================================== */
		if       (label.equals("LMA History")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (FOCUSING_FIELD==null) return;
			FOCUSING_FIELD.setDebugLevel(DEBUG_LEVEL);
			FOCUSING_FIELD.setAdjustMode(false);
			FOCUSING_FIELD.LevenbergMarquardt(
					null, // measurement
					true, // open dialog
//					false, // filterZ
					DEBUG_LEVEL); //boolean openDialog, int debugLevel){
			return;
		}
/* ======================================================================== */
		if       (label.equals("List curv pars")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (FOCUSING_FIELD==null) return;
			FOCUSING_FIELD.setDebugLevel(DEBUG_LEVEL);
			FOCUSING_FIELD.listParameters("Field curvature measurement parameters",null); // to screen
			return;
		}
/* ======================================================================== */
		if       (label.equals("List curv data")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (FOCUSING_FIELD==null) return;
			FOCUSING_FIELD.setDebugLevel(DEBUG_LEVEL);
			FOCUSING_FIELD.listData("Field curvature measurement data",null); // to screen
			return;
		}

/* ======================================================================== */
		if       (label.equals("List qualB")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (FOCUSING_FIELD==null) return;
			FOCUSING_FIELD.setDebugLevel(DEBUG_LEVEL);
			FOCUSING_FIELD.listScanQB(); // to screen
			return;
		}
/* ======================================================================== */
		if       (label.equals("List curv")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (FOCUSING_FIELD==null) return;
			FOCUSING_FIELD.setDebugLevel(DEBUG_LEVEL);
			FOCUSING_FIELD.listCombinedResults(); // to screen
			return;
		}
/* ======================================================================== */
		if       (label.equals("Show curv corr")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (FOCUSING_FIELD==null) return;
			FOCUSING_FIELD.setDebugLevel(DEBUG_LEVEL);
			FOCUSING_FIELD.showCurvCorr(); // to screen
			return;
		}
/* ======================================================================== */
		if       (label.equals("Test measurement")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (FOCUSING_FIELD==null) return;
			FOCUSING_FIELD.setDebugLevel(DEBUG_LEVEL);
			FOCUSING_FIELD.testMeasurement();
			return;
		}
/* ======================================================================== */
		if       (label.equals("Optimize qualB")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (FOCUSING_FIELD==null) return;
			FOCUSING_FIELD.setDebugLevel(DEBUG_LEVEL);
			FOCUSING_FIELD.testQualB(true); //     public double[] testQualB(boolean interactive)
			return;
		}
		
				
		
/* ======================================================================== */
		if       (label.equals("Focus/Tilt LMA")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (FOCUSING_FIELD==null) return;
			FOCUSING_FIELD.setDebugLevel(DEBUG_LEVEL);
			if (!FOCUS_MEASUREMENT_PARAMETERS.cameraIsConfigured) {
				if (CAMERAS.showDialog("Configure cameras interface", 1, true)){
					FOCUS_MEASUREMENT_PARAMETERS.cameraIsConfigured=true;
					if (!FOCUS_MEASUREMENT_PARAMETERS.getLensSerial()) return;
					// reset histories
					MOTORS.clearPreFocus();
					MOTORS.clearHistory();
				} else {
					IJ.showMessage("Error","Camera is not configured\nProcess canceled");
					return;
				}
			}
			// Just for old focal distance calculation
			MOTORS.focusingHistory.optimalMotorPosition( // recalculate calibration to estimate current distance from center PSF
					FOCUS_MEASUREMENT_PARAMETERS,
	    			MOTORS.getMicronsPerStep(), //double micronsPerStep,
	    			DEBUG_LEVEL);
			
			while (adjustFocusTiltLMA());
			return;
		}
//		
/* ======================================================================== */
		if       (label.equals("Show PSF")) {
			
//TODO: change name to have a history step there			
			if (PSF_KERNEL_MAP==null){
				IJ.showMessage("Warning","PSF_KERNEL_MAP is null, nothing to show" );
				return;
			}
			double [][][][] psfRGB=new double [PSF_KERNEL_MAP.length][PSF_KERNEL_MAP[0].length][][];
			int [] rgbChn={1,5,2};
			String [] rgbNames={"Red","Green","Blue"};
			for (int tileY=0;tileY< PSF_KERNEL_MAP.length;tileY++) for (int tileX=0;tileX< PSF_KERNEL_MAP[0].length;tileX++){
				if (PSF_KERNEL_MAP[tileY][tileX]!=null){
					psfRGB[tileY][tileX]=new double [3][];
					for (int rgbi=0;rgbi<3;rgbi++) psfRGB[tileY][tileX][rgbi]=PSF_KERNEL_MAP[tileY][tileX][rgbChn[rgbi]];
					
				} else psfRGB[tileY][tileX]=null;
			}
			ImageStack mergedStack=mergeKernelsToStack(psfRGB,rgbNames);
			if (mergedStack==null) {
				IJ.showMessage("Error","No PSF kernels to show");
				return;
			}
			ImagePlus imp_psf=SDFA_INSTANCE.showImageStack(mergedStack, imp_sel.getTitle()+"-FOCUS-PSF");
			if (PSF_SAVE_FILE) {
				String dir="";
				if (imp_sel.getOriginalFileInfo()!=null) dir=imp_sel.getOriginalFileInfo().directory;
				String path=dir+"FOCUS-PSF-"+imp_sel.getTitle();
				if (DEBUG_LEVEL>1) {
					System.out.println("Saving result to "+path);
				}
				IJ.saveAs(imp_psf,"tif",path);
			}
			return;
		}
/* ======================================================================== */
//	private  ImageStack mergeKernelsToStack(double [][][][] kernels,String [] names) { // use oldStack.getSliceLabels() to get names[]
//
		if       (label.equals("Show Grid")) {
			if (matchSimulatedPattern==null){
				IJ.showMessage("Warning","matchSimulatedPattern is null, nothing to show" );
				return;
			}
			matchSimulatedPattern.showFlatFieldForGrid();
			matchSimulatedPattern.showFFCorrectedGrid();
			matchSimulatedPattern.showFocusMask();
			matchSimulatedPattern.showUVIndex();
			return;

		}
		
/* ======================================================================== */
		if       (label.equals("Calibrate Distance")) {
			MOTORS.setDebug(FOCUS_MEASUREMENT_PARAMETERS.motorDebug);
			checkSerialAndRestore(); // returns true if did not change or was restored 
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
//			System.out.println("microns per step="+MOTORS.getMicronsPerStep());
			MOTORS.setLinearReductionRatio(FOCUS_MEASUREMENT_PARAMETERS.linearReductionRatio);
//TODO: - save/use hysteresis calibration data			
			int [] newPos=MOTORS.focusingHistory.calibrateLensDistance(
					null, // do not save
//					0,
//					MOTORS.historySize()-2, // all but last?
					FOCUS_MEASUREMENT_PARAMETERS,
//	    			true, //boolean sameTiltOnly,
	    			true, //centerOnly,
//	    			-5000, //double xMin,
//	    			200, //double xMax,
	    			MOTORS.getMicronsPerStep(), //double micronsPerStep,
	    			DEBUG_LEVEL); //  int debugLevel
			int  measuredHysteresis=0;
			if (FOCUS_MEASUREMENT_PARAMETERS.scanHysteresis) {
				measuredHysteresis=(int) MOTORS.focusingHistory.getMeasuredHysteresis();
				// report measured and configured hysteresis, suggest to modify and/or check hardware.
				if (MASTER_DEBUG_LEVEL>0) System.out.println("Measured hysteresis is "+measuredHysteresis+" steps, configured hysteresis is set to "+
						FOCUS_MEASUREMENT_PARAMETERS.motorHysteresis+" motor steps");
			}

			if (FOCUS_MEASUREMENT_PARAMETERS.lensDistanceMoveToGoal && (newPos!=null)){
				moveAndMaybeProbe(
						true, // just move, not probe
						newPos, // position maybe updated by the dialog
						MOTORS,
						CAMERAS,
						LENS_DISTORTION_PARAMETERS,
						matchSimulatedPattern, // should not bee null
						FOCUS_MEASUREMENT_PARAMETERS,
						PATTERN_DETECT,
						DISTORTION,
						SIMUL,
						COMPONENTS,
						OTF_FILTER,
						PSF_PARS,
						THREADS_MAX,
						UPDATE_STATUS,
						MASTER_DEBUG_LEVEL,
						DISTORTION.loop_debug_level);

			}
			saveCurrentConfig();
			return;
		}
/* ======================================================================== */
		if       (label.equals("Lasers Toggle")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			checkSerialAndRestore(); // returns true if did not change or was restored 
			UV_LED_LASERS.debugLevel=DEBUG_LEVEL;
			UV_LED_LASERS.lasersToggle(FOCUS_MEASUREMENT_PARAMETERS);
			return;
		}		
/* ======================================================================== */
		if       (label.equals("UV on")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			UV_LED_LASERS.debugLevel=DEBUG_LEVEL;
			if (UV_LED_LASERS.uvControl(FOCUS_MEASUREMENT_PARAMETERS)) {
				MOTORS.setEnable(false); // disable motors if LED is/was on
			}
			return;
		}		
/* ======================================================================== */
		if       (label.equals("UV off")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			checkSerialAndRestore(); // returns true if did not change or was restored 
			UV_LED_LASERS.debugLevel=DEBUG_LEVEL;
			UV_LED_LASERS.uvOff(FOCUS_MEASUREMENT_PARAMETERS);
			// need update amp-sec!!
			// get final image in any case, even if not asked for
			getAndSaveImage(
					false, // boolean alwaysShow, // true overwrites focusMeasurementParameters.showResults
					false, //boolean alwaysSave, // true overwrites focusMeasurementParameters.saveResults
					CAMERAS,
					LENS_DISTORTION_PARAMETERS,
					FOCUS_MEASUREMENT_PARAMETERS,
					UPDATE_STATUS,
					MASTER_DEBUG_LEVEL);
			saveCurrentConfig();
			return;
		}	
/* ======================================================================== */
		if       (label.equals("Camera Power Cycled")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			checkSerialAndRestore(); // Compare S/N before it was will be reset 
			MOTORS.setDebug(FOCUS_MEASUREMENT_PARAMETERS.motorDebug);
			MOTORS.resetInitialization();
			CAMERAS.resetInitialization();
			FOCUS_MEASUREMENT_PARAMETERS.cameraIsConfigured=false;
//			FOCUS_MEASUREMENT_PARAMETERS.serialNumber="";
//			updateSerial(FOCUS_MEASUREMENT_PARAMETERS); // is it still needed?
			return;
		}	
/* ======================================================================== */
		if       (label.equals("Restore SFE Latest")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			restoreSFELatest();
			return;
		}	
/* ======================================================================== */
		if       (label.equals("List SFE")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
//			MOTORS.setDebug(FOCUS_MEASUREMENT_PARAMETERS.motorDebug);
//			MOTORS.resetInitialization();
//			CAMERAS.resetInitialization();
//			FOCUS_MEASUREMENT_PARAMETERS.serialNumber="";
			String dir= CalibrationFileManagement.selectDirectory(
					false, //true, // smart,
					false, //true, // newAllowed, // save  
					"Superdirectory to read SFE focusing results (with SFE serial subdirs)", // title
					"Select results superdirectory (with SFE serial subdirs)", // button
					null, // filter
					FOCUS_MEASUREMENT_PARAMETERS.resultsSuperDirectory);
			if (dir!=null) FOCUS_MEASUREMENT_PARAMETERS.resultsSuperDirectory=dir;
			else {
				String msg="SFE directory is not selected";
				IJ.showMessage(msg);
				System.out.println("Error: "+msg);
				return;
			}
			File dFile=new File(dir);
			File[] sfeList=dFile.listFiles(); // all files
			if ((sfeList==null) || (sfeList.length==0)){
//				String msg="No SFE subdirectories in "+dir;
				String msg="Empty directory: "+dir;
				IJ.showMessage(msg);
				System.out.println("Error: "+msg);
				return;
			}
			boolean stageResults=false;
			boolean multiLensPerSFE=false;
			boolean blankRepeats=true;
			int minThermalState=50;
			boolean decodeManufacturingState=true;
			GenericDialog gd = new GenericDialog("Listing measured SFE");
			gd.addCheckbox("Report all lenses for the same DFE (false - only the latest)",multiLensPerSFE);
			gd.addCheckbox("Include each manufacturing state result (false - only final)",stageResults);
			gd.addCheckbox("Blank lens numbers if they do not change",blankRepeats);
			gd.addCheckbox("Decode manufacturing state",decodeManufacturingState);
			int [] manufacturingIndexMod=FOCUS_MEASUREMENT_PARAMETERS.getManufacturingIndexMod(minThermalState);
			gd. addChoice("Minimal manufacturing state to show thermal parameters",
					FOCUS_MEASUREMENT_PARAMETERS.manufacturingStateNames,
					FOCUS_MEASUREMENT_PARAMETERS.manufacturingStateNames[manufacturingIndexMod[0]]);
			int maxMod=9;
			if (manufacturingIndexMod[0]<(FOCUS_MEASUREMENT_PARAMETERS.manufacturingStateValues.length-1)){
				maxMod=FOCUS_MEASUREMENT_PARAMETERS.manufacturingStateValues[manufacturingIndexMod[0]+1]-
				FOCUS_MEASUREMENT_PARAMETERS.manufacturingStateValues[manufacturingIndexMod[0]]-1;
			}
			gd.addNumericField("Optional manufacturing state modifier for thermal parameters (0.."+maxMod+")",      manufacturingIndexMod[1], 0,1,"");
    		gd.showDialog();
    		if (gd.wasCanceled()) return;
    		multiLensPerSFE=gd.getNextBoolean();
			stageResults=gd.getNextBoolean();
			blankRepeats=gd.getNextBoolean();
			decodeManufacturingState=gd.getNextBoolean();
			int manIndex=                    gd.getNextChoiceIndex();
			int manMod=                (int) gd.getNextNumber();
			if (manMod<0)           manMod=0;
			else if (manMod>maxMod) manMod=maxMod;
			minThermalState=FOCUS_MEASUREMENT_PARAMETERS.manufacturingStateValues[manIndex]+manMod;
			
			Pattern upperHex6Format=Pattern.compile("[0-9A-F]{6}");
			Map<Integer,File> sortedSFE = new TreeMap<Integer,File>();
			for (int iSubDir=0;iSubDir<sfeList.length;iSubDir++) if (sfeList[iSubDir].isDirectory()){
				Matcher matcher=upperHex6Format.matcher(sfeList[iSubDir].getName());
				if (matcher.find()){ // only 6-character serial numbers accepted
					sortedSFE.put(new Integer(Integer.parseInt(matcher.group(0),16)), sfeList[iSubDir]);
				}
			}
			sfeList=new File [sortedSFE.size()];
	        int numSFEDir = 0;
	        for (Map.Entry<Integer, File> e : sortedSFE.entrySet()) {
	        	sfeList[numSFEDir++] = e.getValue();
	        }
			LensAdjustment.FocusMeasurementParameters[][][] sfeParameters=new LensAdjustment.FocusMeasurementParameters [sfeList.length][][];
			LensAdjustment.FocusMeasurementParameters currentSFEParameters=FOCUS_MEASUREMENT_PARAMETERS.clone(); // save current parameters to restore later
			int [][]sensorDimensions=new int [sfeList.length][2];
			Pattern fileNameFormat0=Pattern.compile("LENS(\\d+)");
			Pattern fileNameFormat= Pattern.compile("LENS(\\d+)-S(\\d+)");
			Pattern timestampFormat= Pattern.compile("-(\\d+)_(\\d+)");
			class timestampFile{
				double timestamp;
				File file;
				public timestampFile(double timestamp, File file){
					setFile(file);
					setTimestamp(timestamp);
				}
				File getFile () {return this.file;}
				double getTimestamp() {return this.timestamp;}
				void setFile(File file) {this.file=file;}
				void setTimestamp(double timestamp) {this.timestamp=timestamp;}
			}
			for (int numSFE=0;numSFE<sfeList.length;numSFE++){
				sfeParameters[numSFE]=null;
				String sfeDirPath=sfeList[numSFE].getAbsolutePath();
				File[] fileList=sfeList[numSFE].listFiles(new Filter(".conf-xml"));
				if ((fileList==null) || (fileList.length==0)) {
					String msg="No configuration filers in "+sfeDirPath;
					if (DEBUG_LEVEL>1) System.out.println(msg);
					continue;
				}
				// find all different manufacturing states
				Map<Integer,Map<Integer,timestampFile>> lensStateMap = new TreeMap<Integer,Map<Integer,timestampFile>>();
				for (int i=0;i<fileList.length;i++) {
					Integer lens=-1;
					Integer manState=-1;
					Matcher matcher=fileNameFormat.matcher(fileList[i].getName());
					if (matcher.find()){
						if ((matcher.group(1)!=null) && multiLensPerSFE) lens= Integer.parseInt(matcher.group(1));
						if ((matcher.group(2)!=null) && stageResults && multiLensPerSFE) manState=Integer.parseInt(matcher.group(2));
					} else {
						matcher=fileNameFormat0.matcher(fileList[i].getName());
						if (matcher.find()){
							if ((matcher.group(1)!=null) && multiLensPerSFE) lens=Integer.parseInt(matcher.group(1));
						}
					}
					Matcher tsMatcher=timestampFormat.matcher(fileList[i].getName());
					double ts=0.0;
					if (tsMatcher.find()){
						String sts="";
						if (tsMatcher.group(1)!=null) sts=tsMatcher.group(1);
						if (tsMatcher.group(2)!=null) sts+="."+tsMatcher.group(2);
						ts=Double.parseDouble(sts);
					} else {
						System.out.println("Failed to find timestamp in "+fileList[i].getName()+", lens="+lens+" manState="+manState);
					}
					if (DEBUG_LEVEL>1){
						System.out.println (i+": "+fileList[i].getName()+ " lens="+lens+" manState="+manState+" timestamp="+ts);
					}
					if (lensStateMap.get(lens)==null) {
						lensStateMap.put(lens,new TreeMap<Integer,timestampFile>());
					}
					if (lensStateMap.get(lens).get(manState)==null) {
						lensStateMap.get(lens).put(manState,new timestampFile(ts, fileList[i]));
					} else {
						if (ts>lensStateMap.get(lens).get(manState).getTimestamp()){
							lensStateMap.get(lens).get(manState).setTimestamp(ts);
							lensStateMap.get(lens).get(manState).setFile(fileList[i]);
						}
					}
				}
				sfeParameters[numSFE]=new LensAdjustment.FocusMeasurementParameters[lensStateMap.size()][];
				int iLens=0;
				for (Map.Entry<Integer,Map<Integer,timestampFile>> eLens : lensStateMap.entrySet()) {
					Map<Integer,timestampFile> stateMap= eLens.getValue();
					sfeParameters[numSFE][iLens]=new LensAdjustment.FocusMeasurementParameters[stateMap.size()];
					int iState=0;
					for (Map.Entry<Integer,timestampFile> eState : stateMap.entrySet()) {
						File file = eState.getValue().getFile();
						String path=file.getAbsolutePath();
						if (DEBUG_LEVEL>1) {
							String msg=sfeList[numSFE].getName()+": found latest configuration file "+path;
							System.out.println(msg);
						}
						FOCUS_MEASUREMENT_PARAMETERS.resetResults();
						PROPERTIES=new Properties(); // reset properties
						FOCUS_MEASUREMENT_PARAMETERS.manufacturingState=0; // new SFE - reset for old format
						loadProperties(path, null, true, PROPERTIES);
						if (DEBUG_LEVEL>2) System.out.println("numSFE="+numSFE+" iLens="+iLens+" iState="+iState+" sfeParameters.length="+sfeParameters.length);
						if (DEBUG_LEVEL>2) System.out.println("sfeParameters["+numSFE+"].length="+sfeParameters[numSFE].length);
						if (DEBUG_LEVEL>2) System.out.println("sfeParameters["+numSFE+"]["+iLens+"].length="+sfeParameters[numSFE][iLens].length);
						sfeParameters[numSFE][iLens][iState]=FOCUS_MEASUREMENT_PARAMETERS.clone();
						sensorDimensions[numSFE][0]=EYESIS_CAMERA_PARAMETERS.sensorWidth; // should be the same for all lenses/states
						sensorDimensions[numSFE][1]=EYESIS_CAMERA_PARAMETERS.sensorHeight;
						iState++;
					}			        	
					iLens++;
				}
			}
			if (DEBUG_LEVEL>0) System.out.println("Found "+sfeList.length+" SFE directories with configuration files");
    		String header="#\tS/N"+
    		"\tLens"+
    		"\tSHFT"+
    		"\tMS"+
	    	"\tlastKT"+ // focal distance temperature coefficient (um/C), measured from last run 
	    	"\tlastFD20"+ // focal distance for 20C, measured from last run 
	    	"\tallHistoryKT"+   // focal distance temperature coefficient (um/C), measured from all the measurement histgory 
	    	"\tallHistoryFD20"+ // focal distance for 20C, measured from  all the measurement histgory
	    	"\tfDistance"+ // last measured focal distance
	    	"\ttiltX"+ // last measured tilt X
	    	"\ttiltY"+ // last measured tilt Y
	    	"\tR50"+   // last measured R50 (average PSF 50% level radius, pixels - somewhat larged than actual because of measurement settings)
	    	"\tA50"+   // last measured A50 (simailar, but R^2 are averaged) 
	    	"\tB50"+   // last measured B50 (simailar, but R^4 are averaged)
	    	"\tRC50"+  // last measured RC50(R50 calculated only for the 2 center samples)
//	    	"\tPX0"+ // lens center shift, X
//	    	"\tPY0"+ // lens center shift, Y
	    	"\tdPX"+ // lens center shift, X
	    	"\tdPY"+ // lens center shift, Y
//	    	"\tPSI"+ // SFE rotation (from grid)
	    	"\tROT"+ // SFE rotation (from head lasers)
	    	"\tFocalLength"+ // lens focal length
	    	"\tComment"; // Comment
    		StringBuffer sb = new StringBuffer();
    		for ( int iSFE=0;iSFE<sfeParameters.length;iSFE++)
    			for (int iLens=0;iLens<sfeParameters[iSFE].length;iLens++)
    				for (int iState=0;iState<sfeParameters[iSFE][iLens].length;iState++){
    					if (!blankRepeats || (iState==0)) sb.append(iSFE); // skip numbers for the same lens
    					sb.append("\t"+sfeParameters[iSFE][iLens][iState].serialNumber);
    					if (!blankRepeats || (iState==0)) { // to visually separate different lenses
    						sb.append("\t"+sfeParameters[iSFE][iLens][iState].lensSerial);
    					} else {
    						sb.append("\t");
    					}
    					sb.append("\t"+sfeParameters[iSFE][iLens][iState].centerDeltaX);
    					sb.append("\t"+sfeParameters[iSFE][iLens][iState].manufacturingState);
    					if (decodeManufacturingState){
    						int [] manIndexMod=FOCUS_MEASUREMENT_PARAMETERS.getManufacturingIndexMod(sfeParameters[iSFE][iLens][iState].manufacturingState);
    						sb.append(" - "+ FOCUS_MEASUREMENT_PARAMETERS.manufacturingStateNames[manIndexMod[0]]);
    						if (manIndexMod[1]>0) sb.append(" +"+manIndexMod[1]); 
    					}
    					if (sfeParameters[iSFE][iLens][iState].manufacturingState >=minThermalState) {
    						sb.append("\t"+IJ.d2s(sfeParameters[iSFE][iLens][iState].result_lastKT,3));
    						sb.append("\t"+IJ.d2s(sfeParameters[iSFE][iLens][iState].result_lastFD20,3)); 
    						sb.append("\t"+IJ.d2s(sfeParameters[iSFE][iLens][iState].result_allHistoryKT,3)); 
    						sb.append("\t"+IJ.d2s(sfeParameters[iSFE][iLens][iState].result_allHistoryFD20,3));
    					} else {
    						sb.append("\t\t\t\t");
    					}
    					sb.append("\t"+IJ.d2s(sfeParameters[iSFE][iLens][iState].result_fDistance,3));
    					sb.append("\t"+IJ.d2s(sfeParameters[iSFE][iLens][iState].result_tiltX,3));
    					sb.append("\t"+IJ.d2s(sfeParameters[iSFE][iLens][iState].result_tiltY,3));
    					sb.append("\t"+IJ.d2s(sfeParameters[iSFE][iLens][iState].result_R50,3));
    					sb.append("\t"+IJ.d2s(sfeParameters[iSFE][iLens][iState].result_A50,3)); 
    					sb.append("\t"+IJ.d2s(sfeParameters[iSFE][iLens][iState].result_B50,3));
    					sb.append("\t"+IJ.d2s(sfeParameters[iSFE][iLens][iState].result_RC50,3));
    					//    	    	if (iState==0) { // to visually separate different lenses
    					sb.append("\t"+IJ.d2s(sfeParameters[iSFE][iLens][iState].result_PX0-0.5*sensorDimensions[iSFE][0],3));
    					sb.append("\t"+IJ.d2s(sfeParameters[iSFE][iLens][iState].result_PY0-0.5*sensorDimensions[iSFE][1],3));
    					//    	    	} else {
    					//    	    		sb.append("\t\t");
    					//    	    	}
//    					sb.append("\t"+IJ.d2s(sfeParameters[iSFE][iLens][iState].result_PSI,3));
    					sb.append("\t"+IJ.d2s(sfeParameters[iSFE][iLens][iState].result_ROT,3));
    					sb.append("\t"+IJ.d2s(sfeParameters[iSFE][iLens][iState].result_FocalLength,3));
    					sb.append("\t"+       sfeParameters[iSFE][iLens][iState].comment);
    					sb.append("\n");
    				}
    		new TextWindow("SFE_calibration_results", header, sb.toString(), 1600,1000);
    		FOCUS_MEASUREMENT_PARAMETERS=currentSFEParameters; // restore original parameters
    		return;
		}	
		
/* ======================================================================== */
		if       (label.equals("Acquire&Save")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			getAndSaveImage(
					false, // boolean alwaysShow, // true overwrites focusMeasurementParameters.showResults
					true, //boolean alwaysSave, // true overwrites focusMeasurementParameters.saveResults
					CAMERAS,
					LENS_DISTORTION_PARAMETERS,
					FOCUS_MEASUREMENT_PARAMETERS,
					UPDATE_STATUS,
					MASTER_DEBUG_LEVEL);
			return;
		}	
/* ======================================================================== */
		if       (label.equals("No-move measure")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (!FOCUS_MEASUREMENT_PARAMETERS.cameraIsConfigured) {
				if (CAMERAS.showDialog("Configure cameras interface", 1, true)){
					FOCUS_MEASUREMENT_PARAMETERS.cameraIsConfigured=true;
					if (!FOCUS_MEASUREMENT_PARAMETERS.getLensSerial()) return;
					// reset histories
					MOTORS.clearPreFocus();
					MOTORS.clearHistory();

//					IJ.showMessage("Notice","Make sure camera is in TRIG=4 mode, JP4, correct exposure/white balance...");
				} else {
					IJ.showMessage("Error","Camera is not configured\nProcess canceled");
					return;
				}
			}
			MOTORS.focusingHistory.optimalMotorPosition( // recalculate calibration to estimate current distance from center PSF
					FOCUS_MEASUREMENT_PARAMETERS,
	    			MOTORS.getMicronsPerStep(), //double micronsPerStep,
	    			DEBUG_LEVEL);
			moveAndMaybeProbe(
					true, // just move, not probe
					null, // no move, just measure
					MOTORS,
					CAMERAS,
					LENS_DISTORTION_PARAMETERS,
					matchSimulatedPattern, // should not bee null - is null after grid center!!!
					FOCUS_MEASUREMENT_PARAMETERS,
					PATTERN_DETECT,
					DISTORTION,
					SIMUL,
					COMPONENTS,
					OTF_FILTER,
					PSF_PARS,
					THREADS_MAX,
					UPDATE_STATUS,
					MASTER_DEBUG_LEVEL,
					DISTORTION.loop_debug_level);
			return;
		}	
/* ======================================================================== */
		if       (label.equals("Temp. Scan") || label.equals("Focus Average")) {
			checkSerialAndRestore(); // returns true if did not change or was restored 
			boolean modeAverage=label.equals("Focus Average");
			boolean noTiltScan=true;
			boolean useLMA=true;
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			MOTORS.focusingHistory.optimalMotorPosition( // recalculate calibration to estimate current distance from center PSF
					FOCUS_MEASUREMENT_PARAMETERS,
	    			MOTORS.getMicronsPerStep(), //double micronsPerStep,
	    			DEBUG_LEVEL);
			GenericDialog gd = new GenericDialog("Temperature Scan");
			if (modeAverage) {
				gd.addMessage("This program will repetitively measure focal distance for specified time and average (and record) results.");
			} else {
				gd.addMessage("This program will repetitively measure focal distance and temperature, recording the results.");
				gd.addMessage("Temperature has to be varied separately.");
				gd.addMessage("Use \"List History\" to see the results.");
				gd.addCheckbox("Turn lasers off to protect from overheating",true);
			}
			gd.addCheckbox("Erase previous measurement history",modeAverage);
			gd.addCheckbox("Allow tilt scan when looking for the best fit",!noTiltScan);
			gd.addCheckbox("Use LMA calculations for focus/tilt",useLMA);
			double scanMinutes=modeAverage?1.0:30.0;
    		gd.addNumericField("Measure for ",   scanMinutes , 1,5," minutes");
    		gd.showDialog();
    		if (gd.wasCanceled()) return;
    		if (!modeAverage && gd.getNextBoolean()){
    			UV_LED_LASERS.debugLevel=DEBUG_LEVEL;
    			UV_LED_LASERS.lasersOff(FOCUS_MEASUREMENT_PARAMETERS);
				if (MASTER_DEBUG_LEVEL>0) System.out.println ("Turned laser pointers off to protect from overheating");
    		}
			if (gd.getNextBoolean()) MOTORS.clearHistory();
			noTiltScan=!gd.getNextBoolean();
			useLMA=gd.getNextBoolean();
//			int startHistPos=MOTORS.historySize();
			scanMinutes=gd.getNextNumber();
			long startTime=System.nanoTime();
			long endTime=startTime+(long) (6E10*scanMinutes);
			if (MASTER_DEBUG_LEVEL>0) System.out.println(" startTime= "+startTime+", endTime="+endTime);
			int runs=0;
			while (System.nanoTime()<endTime){


				moveAndMaybeProbe(
						true, // just move, not probe
						null, // no move, just measure
						MOTORS,
						CAMERAS,
						LENS_DISTORTION_PARAMETERS,
						matchSimulatedPattern, // should not be null
						FOCUS_MEASUREMENT_PARAMETERS,
						PATTERN_DETECT,
						DISTORTION,
						SIMUL,
						COMPONENTS,
						OTF_FILTER,
						PSF_PARS,
						THREADS_MAX,
						UPDATE_STATUS,
						MASTER_DEBUG_LEVEL,
						DISTORTION.loop_debug_level);
				runs++;
				long secondsLeft=(long) (0.000000001*(endTime-System.nanoTime()));
				if (secondsLeft<0) secondsLeft=0;
				if (MASTER_DEBUG_LEVEL>0) System.out.println(" Measured "+runs+", "+secondsLeft+" seconds left");
				if (this.SYNC_COMMAND.stopRequested.get()>0) {
					System.out.println("User requested stop");
					break;
				}

			}
			// LMA version
			
			double pX0=FOCUS_MEASUREMENT_PARAMETERS.result_PX0;
			double pY0=FOCUS_MEASUREMENT_PARAMETERS.result_PY0;
			double [][][] sampleCoord=FOCUS_MEASUREMENT_PARAMETERS.sampleCoordinates( //{x,y,r}
					pX0,   // lens center on the sensor
					pY0);
			FocusingField ff= null;
			if (useLMA){
				ff= new FocusingField(
						EYESIS_CAMERA_PARAMETERS.getSensorWidth(),
						EYESIS_CAMERA_PARAMETERS.getSensorHeight(),
						0.001*EYESIS_CAMERA_PARAMETERS.getPixelSize(0), //subCamera_0.pixelSize,
						FOCUS_MEASUREMENT_PARAMETERS.serialNumber,
						FOCUS_MEASUREMENT_PARAMETERS.lensSerial, // String lensSerial, // if null - do not add average
						FOCUS_MEASUREMENT_PARAMETERS.comment, // String comment,
						pX0,
						pY0,
						sampleCoord,
						this.SYNC_COMMAND.stopRequested);
				ff.setDebugLevel(DEBUG_LEVEL);
				ff.setAdjustMode(false);
				if (PROPERTIES!=null) ff.getProperties("FOCUSING_FIELD.", PROPERTIES,true);
				MOTORS.addCurrentHistoryToFocusingField(
						ff,
						(runs==0)?0:(MOTORS.historySize()-runs),
								MOTORS.historySize()); // all newly acquired
				if (modeAverage && (FOCUSING_FIELD!=null)){ // calculate/show average over the last run - only in "average" mode
					double [] ZTM=FOCUSING_FIELD.averageZTM(noTiltScan,ff); // no tilt scan - faster
					if (MASTER_DEBUG_LEVEL>0) {
						String msg="Failed to calulate average focus/tilt";
						if (ZTM!=null) msg="Average:\n"+
								"Relative focal shift "+IJ.d2s(ZTM[0],3)+"um (absolute - "+IJ.d2s(ZTM[0]+FOCUSING_FIELD.qualBOptimizationResults[0],3)+"um)\n"+
								"Relative horizontal tilt "+IJ.d2s(ZTM[1],3)+"um/mm (absolute - "+IJ.d2s(ZTM[1]+FOCUSING_FIELD.qualBOptimizationResults[1],3)+"um.mm)\n"+
								"Relative vertical tilt "+IJ.d2s(ZTM[2],3)+"um/mm (absolute - "+IJ.d2s(ZTM[2]+FOCUSING_FIELD.qualBOptimizationResults[2],3)+"um.mm)\n"+
								"Suggested M1 "+IJ.d2s(ZTM[3],0)+"steps\n"+
								"Suggested M2 "+IJ.d2s(ZTM[4],0)+"steps\n"+
								"Suggested M3 "+IJ.d2s(ZTM[5],0)+"steps";
						System.out.println(msg);
						IJ.showMessage(msg);
					}
				}
			}
			if (FOCUS_MEASUREMENT_PARAMETERS.saveResults) {
				String dir=getResultsPath(FOCUS_MEASUREMENT_PARAMETERS);
				File dFile=new File(dir);
				if (!dFile.isDirectory() &&  !dFile.mkdirs()) {
					String msg="Failed to create directory "+dir;
					IJ.showMessage(msg);
					throw new IllegalArgumentException (msg);
				}
				String lensPrefix="";
				if (FOCUS_MEASUREMENT_PARAMETERS.includeLensSerial && (FOCUS_MEASUREMENT_PARAMETERS.lensSerial.length()>0)){
//					lensPrefix=String.format("LENS%S-",FOCUS_MEASUREMENT_PARAMETERS.lensSerial);
					lensPrefix=String.format("LENS%S-S%02d-",FOCUS_MEASUREMENT_PARAMETERS.lensSerial,FOCUS_MEASUREMENT_PARAMETERS.manufacturingState);
				}
				String path=dFile+Prefs.getFileSeparator()+lensPrefix+CAMERAS.getLastTimestampUnderscored()+
				(modeAverage?"-summary.csv":"-tempscan.csv");
				if (MASTER_DEBUG_LEVEL>0) System.out.println ((modeAverage?"Saving averaged measurements to ":"Saving temperature measurement log data to ")+path);
				int sensorWidth=2992,sensorHeight=1936;
				if ((LENS_DISTORTIONS!=null) && (LENS_DISTORTIONS.fittingStrategy!=null) && (LENS_DISTORTIONS.fittingStrategy!=null)&&
						(LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData!=null) &&
						(LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData.eyesisCameraParameters!=null)){
					sensorWidth=LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.sensorWidth;
					sensorHeight=LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.sensorHeight;
				}
				if (FOCUSING_FIELD!=null){
					sensorWidth=FOCUSING_FIELD.sensorWidth;
					sensorHeight=FOCUSING_FIELD.sensorHeight;
				}
				MOTORS.listHistory(
						path, // on screen, path - to csv
						FOCUS_MEASUREMENT_PARAMETERS.serialNumber,
						FOCUS_MEASUREMENT_PARAMETERS.lensSerial,
						FOCUS_MEASUREMENT_PARAMETERS.comment,
						FOCUS_MEASUREMENT_PARAMETERS.showHistoryDetails,
						FOCUS_MEASUREMENT_PARAMETERS.showHistorySamples, // separate them?
						FOCUS_MEASUREMENT_PARAMETERS.weightRatioRedToGreen,
						FOCUS_MEASUREMENT_PARAMETERS.weightRatioBlueToGreen,
						FOCUS_MEASUREMENT_PARAMETERS.lensDistanceWeightK,
						FOCUS_MEASUREMENT_PARAMETERS.lensDistanceWeightY,
						label.equals("Focus Average"),
						FOCUS_MEASUREMENT_PARAMETERS.result_PX0-sensorWidth/2,
						FOCUS_MEASUREMENT_PARAMETERS.result_PY0-sensorHeight/2,
						FOCUS_MEASUREMENT_PARAMETERS.result_lastKT,
						FOCUS_MEASUREMENT_PARAMETERS.result_allHistoryKT
						);
				if (useLMA && (ff!=null)){
					String focusingPath=dFile+Prefs.getFileSeparator()+lensPrefix+CAMERAS.getLastTimestampUnderscored()+".history-xml";
					System.out.println("Saving measurement history to "+focusingPath);
					ff.saveXML(focusingPath);
				}
			}
// Calculate and show average distances and tilts for measured history			
			
			if (!modeAverage) {
				
				double [] lastKT=MOTORS.focusingHistory.temperatureLinearApproximation(
						runs,             // number of last samples from history to use, 0 - use all
						FOCUS_MEASUREMENT_PARAMETERS.lensDistanceWeightK,
						FOCUS_MEASUREMENT_PARAMETERS.lensDistanceWeightY
				);
				double [] allKT=MOTORS.focusingHistory.temperatureLinearApproximation(
						0,             // number of last samples from history to use, 0 - use all
						FOCUS_MEASUREMENT_PARAMETERS.lensDistanceWeightK,
						FOCUS_MEASUREMENT_PARAMETERS.lensDistanceWeightY
				);
				FOCUS_MEASUREMENT_PARAMETERS.result_lastKT=lastKT[1];   // focal distance temperature coefficient (um/C), measured from last run 
				FOCUS_MEASUREMENT_PARAMETERS.result_lastFD20=lastKT[0]+20*lastKT[1]; // focal distance for 20C, measured from last run 
				FOCUS_MEASUREMENT_PARAMETERS.result_allHistoryKT=allKT[1];   // focal distance temperature coefficient (um/C), measured from all the measurement histgory 
				FOCUS_MEASUREMENT_PARAMETERS.result_allHistoryFD20=allKT[0]+20*allKT[1]; // focal distance for 20C, measured from  all the measurement histgory
				if (MASTER_DEBUG_LEVEL>0){
					String msg=
						"Focal distance temperature expansion is "+IJ.d2s(FOCUS_MEASUREMENT_PARAMETERS.result_lastKT,2)+
						" ("+IJ.d2s(FOCUS_MEASUREMENT_PARAMETERS.result_allHistoryKT,2)+")"+
						" microns per C, measured from this run (measured from all the history)\n"+
						"Focal distance at 20C is "+IJ.d2s(FOCUS_MEASUREMENT_PARAMETERS.result_lastFD20,2)+
						" ("+IJ.d2s(FOCUS_MEASUREMENT_PARAMETERS.result_allHistoryFD20,2)+")"+
						" microns, measured from this run (measured from all the history)";
					System.out.println(msg);
					IJ.showMessage("Info",msg);
				}
// Now in LMA mode - recalculate and overwrite
				if (useLMA){
					ff= new FocusingField(
							EYESIS_CAMERA_PARAMETERS.getSensorWidth(),
							EYESIS_CAMERA_PARAMETERS.getSensorHeight(),
							0.001*EYESIS_CAMERA_PARAMETERS.getPixelSize(0), //subCamera_0.pixelSize,
							FOCUS_MEASUREMENT_PARAMETERS.serialNumber,
							FOCUS_MEASUREMENT_PARAMETERS.lensSerial, // String lensSerial, // if null - do not add average
							FOCUS_MEASUREMENT_PARAMETERS.comment, // String comment,
							pX0,
							pY0,
							sampleCoord,
							this.SYNC_COMMAND.stopRequested);
					ff.setDebugLevel(DEBUG_LEVEL);
					ff.setAdjustMode(false);
					if (PROPERTIES!=null) ff.getProperties("FOCUSING_FIELD.", PROPERTIES,true);
					MOTORS.addCurrentHistoryToFocusingField(ff); // all, not just newly acquired
					if (MASTER_DEBUG_LEVEL>0){
						System.out.println ("*** Calculating focal distance shift for each of "+MOTORS.historySize()+" recorded measurments ***");
					}
					double [][] allZTM=FOCUSING_FIELD.getAllZTM(noTiltScan,ff); // no tilt scan (supposed to be adjusted
					lastKT=MOTORS.focusingHistory.temperatureLinearApproximation(
							allZTM,
							runs
							);
					allKT=MOTORS.focusingHistory.temperatureLinearApproximation(
							allZTM,
							0
							);
					FOCUS_MEASUREMENT_PARAMETERS.result_lastKT=lastKT[1];   // focal distance temperature coefficient (um/C), measured from last run 
					FOCUS_MEASUREMENT_PARAMETERS.result_lastFD20=lastKT[0]+20*lastKT[1]; // focal distance for 20C, measured from last run 
					FOCUS_MEASUREMENT_PARAMETERS.result_allHistoryKT=allKT[1];   // focal distance temperature coefficient (um/C), measured from all the measurement histgory 
					FOCUS_MEASUREMENT_PARAMETERS.result_allHistoryFD20=allKT[0]+20*allKT[1]; // focal distance for 20C, measured from  all the measurement histgory
					if (MASTER_DEBUG_LEVEL>0){
						String msg=
								"Determined by LMA: Focal distance temperature expansion is "+IJ.d2s(FOCUS_MEASUREMENT_PARAMETERS.result_lastKT,2)+
								" ("+IJ.d2s(FOCUS_MEASUREMENT_PARAMETERS.result_allHistoryKT,2)+")"+
								" microns per C, measured from this run (measured from all the history)\n"+
								"Focal distance at 20C is "+IJ.d2s(FOCUS_MEASUREMENT_PARAMETERS.result_lastFD20,2)+
								" ("+IJ.d2s(FOCUS_MEASUREMENT_PARAMETERS.result_allHistoryFD20,2)+")"+
								" microns, measured from this run (measured from all the history)";
						System.out.println(msg);
						IJ.showMessage("Info",msg);
					}
				}
			}
			saveCurrentConfig();
			return;
		}	
///	"No-move measure"	
		
/* ======================================================================== */
		if       (label.equals("Configure Goniometer")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (GONIOMETER_PARAMETERS.showDialog("Goniometer  Parameters")) return;
/*
			MOTORS.setHysteresis(FOCUS_MEASUREMENT_PARAMETERS.motorHysteresis);
			MOTORS.setCalmMotors(FOCUS_MEASUREMENT_PARAMETERS.motorCalm);
			MOTORS.setDebug(FOCUS_MEASUREMENT_PARAMETERS.motorDebug);

			if (FOCUS_MEASUREMENT_PARAMETERS.configureCamera) {
				if (CAMERAS.showDialog("Configure cameras interface", 1, true)) FOCUS_MEASUREMENT_PARAMETERS.cameraIsConfigured=true;
			}
*/			
			return;
		}
/* ======================================================================== */
		if       (label.equals("Goniometer Scan")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			
			CAMERAS.setNumberOfThreads(THREADS_MAX);
			CAMERAS.debugLevel=DEBUG_LEVEL;

			if (GONIOMETER==null) {
				GONIOMETER= new Goniometer(
						CAMERAS, // CalibrationHardwareInterface.CamerasInterface cameras,
						DISTORTION, //MatchSimulatedPattern.DistortionParameters distortion,
						PATTERN_DETECT, //MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
						EYESIS_CAMERA_PARAMETERS, //Distortions.EyesisCameraParameters eyesisCameraParameters,
						LASER_POINTERS, // MatchSimulatedPattern.LaserPointer laserPointers
						SIMUL,                       //SimulationPattern.SimulParameters  simulParametersDefault,
						GONIOMETER_PARAMETERS, //LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
						DISTORTION_PROCESS_CONFIGURATION
				);
				if (DEBUG_LEVEL>1){
					System.out.println("Initiaslizing Goniometer class");
				}
			} else if (DEBUG_LEVEL>1){
				System.out.println("GONIOMETER was initialized");
			}
			
			// calculate angular size of the target as visible from the camera
			/*
			Distortions.DistortionCalibrationData dcd=(DISTORTION_CALIBRATION_DATA!=null)?DISTORTION_CALIBRATION_DATA:
				new Distortions.DistortionCalibrationData(EYESIS_CAMERA_PARAMETERS);*/
//			double distanceToTarget=dcd.eyesisCameraParameters.GXYZ[2];
			double distanceToTarget=GONIOMETER_PARAMETERS.targetDistance;
			double patternWidth= PATTERN_PARAMETERS.patternWidth;
			double patternHeight=PATTERN_PARAMETERS.patternHeight;
			double targetAngleHorizontal=360*Math.atan(patternWidth/2/distanceToTarget)/Math.PI;
			double targetAngleVertical=  360*Math.atan(patternHeight/2/distanceToTarget)/Math.PI;
			if (DEBUG_LEVEL>0) System.out.println(
					"Using:\n"+
					"Distance from target:          "+IJ.d2s(distanceToTarget,1)+" mm\n"+
					"         Taget width:          "+IJ.d2s(patternWidth,1)+" mm\n"+
					"         Taget height:         "+IJ.d2s(patternHeight,1)+" mm\n"+
					"Taget angular size horizontal: "+IJ.d2s(targetAngleHorizontal,1)+" degrees\n"+
					"Taget angular size vertical:   "+IJ.d2s(targetAngleVertical,1)+" degrees\n"
			);
			GONIOMETER.debugLevel=DEBUG_LEVEL;
			boolean goniometerScanOK=GONIOMETER.scanAndAcquire(
					targetAngleHorizontal,
					targetAngleVertical,
					this.SYNC_COMMAND.stopRequested,
					UPDATE_STATUS);
			System.out.println ("GONIOMETER.scanAndAcquireI() "+(goniometerScanOK?"finished OK":"failed"));
			return;
		}		
		
/* ======================================================================== */
		if       (label.equals("Filter Grids")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (DISTORTION_CALIBRATION_DATA==null) {
				IJ.showMessage("Distortion calibration data required");
				return;
			}
			boolean overwriteAll=(DISTORTION_CALIBRATION_DATA.gIS==null);
			if ((DEBUG_LEVEL>0) && (DISTORTION_CALIBRATION_DATA.gIS!=null)){
				System.out.println("There are "+DISTORTION_CALIBRATION_DATA.getNumberOfEstimated(true)+ "("+DISTORTION_CALIBRATION_DATA.getNumberOfEstimated(false)+") images with estimated orientation");
			}
			GenericDialog gd=new GenericDialog ("Configure grid filter");
    		gd.addCheckbox    ("Reset calibration by predicted grids", false);
    		gd.addNumericField("Minimal number of laser pointers (w/o calibration from predicted)", 2,0);
 //   		gd.addCheckbox    ("Organize all images (false - only enabled)", false);
    		gd.addNumericField("Minimal registered grid period as a fraction of maximal (to filter reflections)", 0.4,2); //was 0.7
    		gd.addCheckbox    ("Reset orientation from the image with most pointers (false - only if it was not set yet)", overwriteAll);
    		gd.addCheckbox    ("Disable (old) images without vignetting information", true);
    		gd.addNumericField("Minimal number of grids in no-pointer images and estimated orientation", 1000,0);
    		gd.showDialog();
    		if (gd.wasCanceled()) return;
    		boolean resetHinted=       gd.getNextBoolean();
    		int minPointers=      (int)gd.getNextNumber();
//    		boolean organizeAllImages= gd.getNextBoolean();
    		double minGridPeriod=      gd.getNextNumber();
    		overwriteAll=      gd.getNextBoolean();
    		boolean disableNoVignetting=gd.getNextBoolean();
    		int minGridsNoPointer= (int) gd.getNextNumber();
    		int [] numImages=DISTORTION_CALIBRATION_DATA.filterImages(
    				resetHinted,
    				minPointers,
    				minGridPeriod,
    				disableNoVignetting,
    				minGridsNoPointer);
    		System.out.println("Number of enabled grid images: "+numImages[0]+
    				", of them new: "+numImages[1]+
    				", disabled without vignetting info: "+numImages[2]+
    				", disabled having less than "+minGridsNoPointer+" nodes and no matched pointers: "+numImages[3]);
    		if (DISTORTION_CALIBRATION_DATA.gIS==null) {
    			int numImageSets=DISTORTION_CALIBRATION_DATA.buildImageSets(false); // from scratch
    			if (DEBUG_LEVEL>0) System.out.println("Image set was empty, built a new one with "+numImageSets+" image sets (\"panoramas\"): ");
    			DISTORTION_CALIBRATION_DATA.updateSetOrientation(null); // restore orientation from (enabled) image files
    			if (DEBUG_LEVEL>0) System.out.println("Setting sets orientation from per-grid image data");
    		}
    		if (overwriteAll) DISTORTION_CALIBRATION_DATA.setInitialOrientation(overwriteAll);  // needed here? modify?
			if ((DEBUG_LEVEL>0) && (DISTORTION_CALIBRATION_DATA.gIS!=null)){
				System.out.println("There are now "+DISTORTION_CALIBRATION_DATA.getNumberOfEstimated(true)+ "("+DISTORTION_CALIBRATION_DATA.getNumberOfEstimated(false)+") images with estimated orientation");
			}
			//	minGridsNoPointer
			
			return; 
		}
/* ======================================================================== */
		if       (label.equals("Update Image Set")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (DISTORTION_CALIBRATION_DATA==null) {
				IJ.showMessage("Distortion calibration data required");
				return;
			}
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			LENS_DISTORTIONS.interactiveUpdateImageSet(
 					DISTORTION_CALIBRATION_DATA,
					EYESIS_CAMERA_PARAMETERS
					);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Remove Outlayers")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			if (LENS_DISTORTIONS.fittingStrategy==null) {
				IJ.showMessage("LENS_DISTORTION.fittingStrategy is not set");
				return;
			}
			LENS_DISTORTIONS.removeOutLayers(
					-1, //int series, (<0 - ask)
					-1  //int numOutLayers  (<0 - ask)
					);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Remove Sets")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			if (LENS_DISTORTIONS.fittingStrategy==null) {
				IJ.showMessage("LENS_DISTORTION.fittingStrategy is not set");
				return;
			}
			LENS_DISTORTIONS.removeOutLayerSets(
					-1  //int numOutLayers  (<0 - ask)
					);
			return;
		}
		//
/* ======================================================================== */
		if  (label.equals("Update Sets Orientation")) { // update "panoramas" orientation after LMA
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (DISTORTION_CALIBRATION_DATA==null) {
				IJ.showMessage("Distortion calibration data required");
				return;
			}
			DISTORTION_CALIBRATION_DATA.updateSetOrientation(null); // update from any of the enabled images - should not be needed, runs after each successful LMA
			return;
		}

/* ======================================================================== */
		if       (label.equals("List Image Sets")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			if (LENS_DISTORTIONS.fittingStrategy==null) {
				IJ.showMessage("LENS_DISTORTION.fittingStrategy is not set");
				return;
			}
			LENS_DISTORTIONS.listImageSets();
			return;
		}
/* ======================================================================== */
		if       (label.equals("Re-calibrate Grids")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			if (LENS_DISTORTIONS.fittingStrategy==null) {
				IJ.showMessage("LENS_DISTORTION.fittingStrategy is not set");
				return;
			}
			GenericDialog gd=new GenericDialog ("Parameters of the filter by predicted grid");
    		gd.addNumericField("Mismatch tolerance of match between the predicted and acquired grid", 5.0, 1,4,"pix");
    		gd.addCheckbox    ("Calibrate by translation (false - orientation only)", true);
    		gd.addCheckbox    ("Process all images (false - only not yet enabled)", false);
    		gd.addCheckbox    ("Ignore laser pointers", false);
    		gd.addCheckbox    ("Process images with estimated orientation and no matched laser pointers", false);
    		gd.addCheckbox    ("Use image sets data if available (false - use camera data)", true);
    		gd.addNumericField("Process only one specified image (<0 - all)", -1, 0);
 
    		gd.showDialog();
    		if (gd.wasCanceled()) return;
    		double hintGridTolerance=       gd.getNextNumber();
    		boolean useHintTolerance=       gd.getNextBoolean();
    		boolean processAll=             gd.getNextBoolean();
    		boolean ignoreLaserPointers=    gd.getNextBoolean();
    		boolean processBlind=           gd.getNextBoolean();
    		boolean useSetsData=            gd.getNextBoolean();
    		int imageNumber=          (int) gd.getNextNumber();
    		if (imageNumber>=0) processAll=true; // no additional filtering
    		
    		int numMatched=LENS_DISTORTIONS.applyHintedGrids(
    				LASER_POINTERS, // MatchSimulatedPattern.LaserPointer laserPointer, // LaserPointer object that specifies actual laser poiners on the target
    				DISTORTION_PROCESS_CONFIGURATION.removeOutOfGridPointers, // boolean removeOutOfGridPointers,
    				(useHintTolerance?hintGridTolerance:0.0),                   //double  hintGridTolerance, // alllowed mismatch (fraction of period) or 0 - orientation only
    				processAll, //boolean processAll, // if true - process all images, false - only disabeld
    				ignoreLaserPointers,
    				processBlind,
    				imageNumber,
    				useSetsData,
    				THREADS_MAX,                 //int threadsMax,
    				UPDATE_STATUS,               // boolean updateStatus,
    				DISTORTION.loop_debug_level, // int mspDebugLevel,
    				MASTER_DEBUG_LEVEL,          //int global_debug_level, // DEBUG_LEVEL
    				MASTER_DEBUG_LEVEL           //int debug_level // debug level used inside loops
    		);
    		System.out.println("Number of (new) matched images: "+numMatched);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Re-calibrate Set")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			if (LENS_DISTORTIONS.fittingStrategy==null) {
				IJ.showMessage("LENS_DISTORTION.fittingStrategy is not set");
				return;
			}
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			GenericDialog gd=new GenericDialog ("Select image set to re-calibrate");
    		gd.addNumericField("Image set number", 0, 0);
    		gd.addCheckbox    ("Set image set parameters from closest, (re-)estimate orientation", true);
    		gd.showDialog();
    		if (gd.wasCanceled()) return;
    		int imageSetNumber=(int) gd.getNextNumber();
    		boolean reEstimate=gd.getNextBoolean();
    		if (imageSetNumber<0) return;
    		if (reEstimate){
//    			boolean OK=    
    					LENS_DISTORTIONS.setSetFromClosestAndEstimateOrientation(
    					imageSetNumber,
    		    		null, //boolean [] selectedImages,
    		    		null, //boolean [] parameterMask,
    		    		LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData,
    		    		LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData.eyesisCameraParameters);
    		}
    		double tiltCenter=LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData.gIS[imageSetNumber].goniometerTilt;
    		double axialCenter=LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData.gIS[imageSetNumber].goniometerAxial;
    		int tiltNumSteps=1;
    		int axialNumSteps=21;
    		double tiltStep=0.25;
    		double axialStep=0.25;
    		boolean processAllImages=true;
    		boolean ignoreLaserPointers=false;
    		double hintGridTolerance=200.0;
    		boolean useSetsData=  true;
			gd=new GenericDialog ("Image set # "+imageSetNumber+" re-calibration without laser pointers");
			gd.addMessage("Strategy 0 should have all parameters but 2 goniometer axes disabled");
			gd.addMessage("Imgages belonging to the set will be selected, possible to check with \"Remove Outlayers\" for strategy 0");
    		gd.addNumericField("Mismatch tolerance of match between the predicted and acquired grid", hintGridTolerance, 1,4,"pix");
    		gd.addCheckbox("Ignore laser pointers", ignoreLaserPointers);
    		gd.addNumericField("Image set tilt", tiltCenter, 2,6,"degrees");
    		gd.addNumericField("Image set axial", axialCenter, 2,6,"degrees");
    		gd.addNumericField("Tilt number of steps", tiltNumSteps, 0,3,"");
    		gd.addNumericField("Axial number of steps", axialNumSteps, 0,3,"");
    		gd.addNumericField("Tilt scan step", tiltStep, 2,5,"degrees");
    		gd.addNumericField("Axial scan step", axialStep, 2,5,"degrees");
    		gd.addCheckbox    ("Use image sets data if available (false - use camera data)", useSetsData);
    		gd.addCheckbox("process all images (false - enabled only)", processAllImages);
    		gd.showDialog();
    		if (gd.wasCanceled()) return;
    		hintGridTolerance=        gd.getNextNumber();
    		ignoreLaserPointers=      gd.getNextBoolean();
    		tiltCenter=               gd.getNextNumber();
    		axialCenter=              gd.getNextNumber();
    		tiltNumSteps=       (int) gd.getNextNumber();
    		axialNumSteps=      (int) gd.getNextNumber();
    		tiltStep=                 gd.getNextNumber();
    		axialStep=                gd.getNextNumber();
    		useSetsData=              gd.getNextBoolean();
    		processAllImages=         gd.getNextBoolean();
    		double [] initialTilt =   new double [tiltNumSteps];
    		double [] initialAxial=   new double [axialNumSteps];
    		double [][] finalTilt =   new double [tiltNumSteps][axialNumSteps];
    		double [][] finalAxial=   new double [tiltNumSteps][axialNumSteps];
    		double [][] finalError=   new double [tiltNumSteps][axialNumSteps];
    		for (int tiltIndex=0; tiltIndex<tiltNumSteps;tiltIndex++) initialTilt[tiltIndex]=tiltCenter+tiltStep*(tiltIndex-0.5*(tiltNumSteps-1));
    		for (int axialIndex=0;axialIndex<axialNumSteps;axialIndex++) initialAxial[axialIndex]=axialCenter+axialStep*(axialIndex-0.5*(axialNumSteps-1));
    		for (int tiltIndex=0; tiltIndex<tiltNumSteps;tiltIndex++) for (int axialIndex=0;axialIndex<axialNumSteps;axialIndex++) { 
    			LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData.gIS[imageSetNumber].goniometerTilt=initialTilt[tiltIndex];
    			LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData.gIS[imageSetNumber].goniometerAxial=initialAxial[axialIndex];
    			if (DEBUG_LEVEL>0) System.out.println("Image Set #"+imageSetNumber+" Initial tilt="+initialTilt[tiltIndex]+" Initial axial="+initialAxial[axialIndex]);
    			int [][] imageSets=LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData.listImages(!processAllImages);
    			int [] imageSet=imageSets[imageSetNumber];
    			if (imageSet==null){
    				IJ.showMessage("Image set #"+imageSetNumber+" is empty");
    				return;
    			}
    			for (int i=0;i<imageSet.length;i++){
    				int imageNumber=imageSet[i];
    				LENS_DISTORTIONS.applyHintedGrids(
    						LASER_POINTERS, // MatchSimulatedPattern.LaserPointer laserPointer, // LaserPointer object that specifies actual laser poiners on the target
    						DISTORTION_PROCESS_CONFIGURATION.removeOutOfGridPointers, // boolean removeOutOfGridPointers,
    						hintGridTolerance,                   //double  hintGridTolerance, // alllowed mismatch (fraction of period) or 0 - orientation only
    						true, //processAll, //boolean processAll, // if true - process all images, false - only disabled
    						ignoreLaserPointers, //true, //ignoreLaserPointers,
    						true, //processBlind,
    						imageNumber,
    						useSetsData,
    						THREADS_MAX,                 //int threadsMax,
    						UPDATE_STATUS,               // boolean updateStatus,
    						DISTORTION.loop_debug_level, // int mspDebugLevel,
    						MASTER_DEBUG_LEVEL,          //int global_debug_level, // DEBUG_LEVEL
    						MASTER_DEBUG_LEVEL           //int debug_level // debug level used inside loops
    				);

    			}
    			// set series 0 to this set images
    			boolean [] selection =LENS_DISTORTIONS.fittingStrategy.selectAllImages(0); // enable all images in series 0
    			for (int i=0;i<selection.length;i++) selection[i]=false;
    			for (int i=0;i<imageSet.length;i++) selection[imageSet[i]]=true;
    			LENS_DISTORTIONS.fittingStrategy.setImageSelection(0,selection);
				LENS_DISTORTIONS.seriesNumber=   0; // start from 0;
				LENS_DISTORTIONS.initFittingSeries(false,LENS_DISTORTIONS.filterForAll,0); // will set this.currentVector
				//this.stopAfterThis[numSeries]
				
    			LENS_DISTORTIONS.fittingStrategy.stopAfterThis[0]=true;
				LENS_DISTORTIONS.stopEachStep=   false;
				LENS_DISTORTIONS.stopEachSeries= false; // will not ask for confirmation after done
				LENS_DISTORTIONS.stopOnFailure=false;
				LENS_DISTORTIONS.lambda=LENS_DISTORTIONS.fittingStrategy.lambdas[0]; // 0.001; // why it does not use fitting series lambda?
				boolean LMA_OK=LENS_DISTORTIONS.LevenbergMarquardt(false); //  skip dialog
				if (LMA_OK) {
					finalTilt[tiltIndex][axialIndex]=LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData.gIS[imageSetNumber].goniometerTilt;
					finalAxial[tiltIndex][axialIndex]=LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData.gIS[imageSetNumber].goniometerAxial;
					finalError[tiltIndex][axialIndex]=LENS_DISTORTIONS.currentRMS;
				} else {
					finalTilt[tiltIndex][axialIndex]=initialTilt[tiltIndex];
					finalAxial[tiltIndex][axialIndex]=initialAxial[axialIndex];
					finalError[tiltIndex][axialIndex]=Double.NaN;
					if (DEBUG_LEVEL>0) System.out.println("----------------- LMA FAILED -------------------------");
				}
    		}
    		int bestTiltIndex=0;
    		int bestAxialIndex=0;
    		double bestRMS=finalError[bestTiltIndex][bestAxialIndex];
    		for (int tiltIndex=0; tiltIndex<tiltNumSteps;tiltIndex++) for (int axialIndex=0;axialIndex<axialNumSteps;axialIndex++) {
    			if (Double.isNaN(bestRMS) || (bestRMS>finalError[tiltIndex][axialIndex])){
    				bestTiltIndex=tiltIndex;
    				bestAxialIndex=axialIndex;
    				bestRMS=finalError[bestTiltIndex][bestAxialIndex];
    			}
    		} 
			if (DEBUG_LEVEL>0) System.out.println("================ Image Set #"+imageSetNumber+" rms="+IJ.d2s(bestRMS, 6)+
					" final tilt="+finalTilt[bestTiltIndex][bestAxialIndex]+" ("+tiltCenter+": "+initialTilt[bestTiltIndex]+") " +
					" final axial="+finalAxial[bestTiltIndex][bestAxialIndex]+" ("+axialCenter+": "+initialAxial[bestAxialIndex]+")");

// repeat with the best indices			
			LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData.gIS[imageSetNumber].goniometerTilt=finalTilt[bestTiltIndex][bestAxialIndex]; //initialTilt[bestTiltIndex];
			LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData.gIS[imageSetNumber].goniometerAxial=finalAxial[bestTiltIndex][bestAxialIndex]; //initialAxial[bestAxialIndex];
			if (DEBUG_LEVEL>0) System.out.println("Image Set #"+imageSetNumber+" Initial tilt="+finalTilt[bestTiltIndex][bestAxialIndex]+" Initial axial="+finalAxial[bestTiltIndex][bestAxialIndex]);
			int [][] imageSets=LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData.listImages(!processAllImages);
			int [] imageSet=imageSets[imageSetNumber];
			if (imageSet==null){
				IJ.showMessage("Image set #"+imageSetNumber+" is empty");
				return;
			}
			for (int i=0;i<imageSet.length;i++){
				int imageNumber=imageSet[i];
				LENS_DISTORTIONS.applyHintedGrids(
						LASER_POINTERS, // MatchSimulatedPattern.LaserPointer laserPointer, // LaserPointer object that specifies actual laser poiners on the target
						DISTORTION_PROCESS_CONFIGURATION.removeOutOfGridPointers, // boolean removeOutOfGridPointers,
						hintGridTolerance,                   //double  hintGridTolerance, // alllowed mismatch (fraction of period) or 0 - orientation only
						true, //processAll, //boolean processAll, // if true - process all images, false - only disabeld
						ignoreLaserPointers, //ignoreLaserPointers,
						true, //processBlind,
						imageNumber,
						useSetsData,
						THREADS_MAX,                 //int threadsMax,
						UPDATE_STATUS,               // boolean updateStatus,
						DISTORTION.loop_debug_level, // int mspDebugLevel,
						MASTER_DEBUG_LEVEL,          //int global_debug_level, // DEBUG_LEVEL
						MASTER_DEBUG_LEVEL           //int debug_level // debug level used inside loops
				);

			}
			// set series 0 to this set images
			boolean [] selection =LENS_DISTORTIONS.fittingStrategy.selectAllImages(0); // enable all images in series 0
			for (int i=0;i<selection.length;i++) selection[i]=false;
			for (int i=0;i<imageSet.length;i++) selection[imageSet[i]]=true;
			LENS_DISTORTIONS.fittingStrategy.setImageSelection(0,selection);
			LENS_DISTORTIONS.seriesNumber=   0; // start from 0;
			LENS_DISTORTIONS.initFittingSeries(false,LENS_DISTORTIONS.filterForAll,0); // will set this.currentVector
			//this.stopAfterThis[numSeries]
			
			LENS_DISTORTIONS.fittingStrategy.stopAfterThis[0]=true;
			LENS_DISTORTIONS.stopEachStep=   false;
			LENS_DISTORTIONS.stopEachSeries= false; // will not ask for confirmation after done
			LENS_DISTORTIONS.lambda=LENS_DISTORTIONS.fittingStrategy.lambdas[0];
			LENS_DISTORTIONS.LevenbergMarquardt(false); //  skip dialog

// save safe settings to run LMA manually			
			LENS_DISTORTIONS.seriesNumber=   0; // start from 0;
			LENS_DISTORTIONS.initFittingSeries(false,LENS_DISTORTIONS.filterForAll,0); // will set this.currentVector
			LENS_DISTORTIONS.stopEachSeries= true; // will not ask for confirmation after done
			LENS_DISTORTIONS.stopOnFailure=true;
			LENS_DISTORTIONS.lambda=LENS_DISTORTIONS.fittingStrategy.lambdas[0];

			if (DEBUG_LEVEL>0) System.out.println("================ Image Set #"+imageSetNumber+" rms="+IJ.d2s(bestRMS, 6)+
					" final tilt="+finalTilt[bestTiltIndex][bestAxialIndex]+" ("+tiltCenter+": "+initialTilt[bestTiltIndex]+") " +
					" final axial="+finalAxial[bestTiltIndex][bestAxialIndex]+" ("+axialCenter+": "+initialAxial[bestAxialIndex]+")");

			return;
		}

		
/* ======================================================================== */
		if       (label.equals("Get Orientation")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (DEBUG_LEVEL>0){
				IJ.showMessage("disabled option");
				return;
			}
//			CAMERAS.setNumberOfThreads(THREADS_MAX);
			if (GONIOMETER==null) {
				GONIOMETER= new Goniometer(
						CAMERAS, // CalibrationHardwareInterface.CamerasInterface cameras,
						DISTORTION, //MatchSimulatedPattern.DistortionParameters distortion,
						PATTERN_DETECT, //MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
						EYESIS_CAMERA_PARAMETERS, //Distortions.EyesisCameraParameters eyesisCameraParameters,
						LASER_POINTERS, // MatchSimulatedPattern.LaserPointer laserPointers
						SIMUL,                       //SimulationPattern.SimulParameters  simulParametersDefault,
						GONIOMETER_PARAMETERS, //LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
						DISTORTION_PROCESS_CONFIGURATION
				);
				if (DEBUG_LEVEL>1){
					System.out.println("Initiaslizing Goniometer class");
				}
			} else if (DEBUG_LEVEL>1){
				System.out.println("GONIOMETER was initialized");
			}
			// initialize needed classes
			DISTORTION_CALIBRATION_DATA=new Distortions.DistortionCalibrationData( // images are not setup yet
	        		EYESIS_CAMERA_PARAMETERS); //EyesisCameraParameters eyesisCameraParameters
			if ((LENS_DISTORTIONS!=null) && (LENS_DISTORTIONS.fittingStrategy!=null)) {
				LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData=DISTORTION_CALIBRATION_DATA;
			}
			
			if (DEBUG_LEVEL>1){
				System.out.println("Initializing DistortionCalibrationData class");
			}
			if (LENS_DISTORTIONS==null) {
				if (LENS_DISTORTION_PARAMETERS==null){
					String msg="LENS_DISTORTION_PARAMETERS is not set";
					System.out.println("Error" + msg);
					IJ.showMessage("Error",msg);
					return;
				}
				if (PATTERN_PARAMETERS==null){
					String msg="PATTERN_PARAMETERS is not set";
					System.out.println("Error" + msg);
					IJ.showMessage("Error",msg);
					return;
				}
				LENS_DISTORTIONS=new Distortions(LENS_DISTORTION_PARAMETERS,PATTERN_PARAMETERS,REFINE_PARAMETERS,this.SYNC_COMMAND.stopRequested);
				if (DEBUG_LEVEL>1){
					System.out.println("Initiaslizing Distortions class");
				}
			} else if (DEBUG_LEVEL>1){
				System.out.println("LENS_DISTORTIONS was initialized");
			}
 
			//	Also - do it once
			// Grid file
			PATTERN_PARAMETERS.debugLevel=DEBUG_LEVEL;
			PATTERN_PARAMETERS.updateStatus=UPDATE_STATUS;
			String gridPathname=PATTERN_PARAMETERS.selectAndRestore(
					true, // skip dialog if not needed
					GONIOMETER_PARAMETERS.gridGeometryFile,
					1); // numStations
			if (gridPathname== null){ // failed to select/open the file
				String msg="Failed to open grid geometry file";
				System.out.println("Error" + msg);
				IJ.showMessage("Error",msg);
				return;
			}
			GONIOMETER_PARAMETERS.gridGeometryFile=gridPathname;

// Find curernt orientation
			double [] currentOrientation=GONIOMETER.estimateOrientation (
					CAMERAS.getImages(1), // last acquired images with number of pointers detected>0
//					DISTORTION, //MatchSimulatedPattern.DistortionParameters distortionParametersDefault,
//					GONIOMETER_PARAMETERS, //LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
//					PATTERN_DETECT, //MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
//					LASER_POINTERS, //MatchSimulatedPattern.LaserPointer laserPointer, // null OK
//					SIMUL,                       //SimulationPattern.SimulParameters  simulParametersDefault,
					DISTORTION_CALIBRATION_DATA, // Distortions.DistortionCalibrationData distortionCalibrationData,
					PATTERN_PARAMETERS,          //Distortions.PatternParameters patternParameters, // should not be null
					LENS_DISTORTIONS,            //Distortions lensDistortions, // should not be null
					COMPONENTS.equalizeGreens,   //boolean equalizeGreens,
					THREADS_MAX,                 // int       threadsMax,
					UPDATE_STATUS,               //boolean   updateStatus,
					DEBUG_LEVEL);                // debug level used inside loops

			if (MASTER_DEBUG_LEVEL>0){
				String msg="Goniometer horizontal="+currentOrientation[0]+"\nGoniometer axial="+currentOrientation[1];
				System.out.println(msg);
				IJ.showMessage("Info",msg);
			}
			return;
		}
/* ======================================================================== */
		if       (label.equals("Test Hinted Grid")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			ImagePlus [] images={WindowManager.getCurrentImage()};
			if (images[0]==null){
				IJ.showMessage("Error","There are no images open\nProcess canceled");
				return;
			}
			   if ((images[0].getProperty("timestamp")==null) || (((String) images[0].getProperty("timestamp")).length()==0)) {
				   (new JP46_Reader_camera(false)).decodeProperiesFromInfo(images[0]);
			   }

			if (GONIOMETER==null) {
				GONIOMETER= new Goniometer(
						CAMERAS, // CalibrationHardwareInterface.CamerasInterface cameras,
						DISTORTION, //MatchSimulatedPattern.DistortionParameters distortion,
						PATTERN_DETECT, //MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
						EYESIS_CAMERA_PARAMETERS, //Distortions.EyesisCameraParameters eyesisCameraParameters,
						LASER_POINTERS, // MatchSimulatedPattern.LaserPointer laserPointers
						SIMUL,                       //SimulationPattern.SimulParameters  simulParametersDefault,
						GONIOMETER_PARAMETERS, //LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
						DISTORTION_PROCESS_CONFIGURATION
				);
				if (DEBUG_LEVEL>1){
					System.out.println("Initiaslizing Goniometer class");
				}
			} else if (DEBUG_LEVEL>1){
				System.out.println("GONIOMETER was initialized");
			}
			if (LENS_DISTORTIONS==null) {
				String msg="LENS_DISTORTIONS is not set";
				System.out.println("Error" + msg);
				IJ.showMessage("Error",msg);
				return;
			}
			if (DISTORTION_CALIBRATION_DATA==null) {
				String msg="DISTORTION_CALIBRATION_DATA is not set";
				System.out.println("Error" + msg);
				IJ.showMessage("Error",msg);
				return;
			}

			GONIOMETER.testHintedTarget (
					images, //CAMERAS.getImages(0), // last acquired images with number of pointers detected>=0
					LENS_DISTORTIONS,            //Distortions lensDistortions, // should not be null
					DISTORTION_CALIBRATION_DATA, // Distortions.DistortionCalibrationData distortionCalibrationData,
					PATTERN_PARAMETERS,          //Distortions.PatternParameters patternParameters, // should not be null
					COMPONENTS.equalizeGreens,   //boolean equalizeGreens,
					THREADS_MAX,                 // int       threadsMax,
					UPDATE_STATUS,               //boolean   updateStatus,
					DEBUG_LEVEL);                // debug level used inside loops
					
					
			return;
		}
/* ======================================================================== */
		if       (label.equals("Test Hinted Grid Cameras")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (GONIOMETER==null) {
				GONIOMETER= new Goniometer(
						CAMERAS, // CalibrationHardwareInterface.CamerasInterface cameras,
						DISTORTION, //MatchSimulatedPattern.DistortionParameters distortion,
						PATTERN_DETECT, //MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
						EYESIS_CAMERA_PARAMETERS, //Distortions.EyesisCameraParameters eyesisCameraParameters,
						LASER_POINTERS, // MatchSimulatedPattern.LaserPointer laserPointers
						SIMUL,                       //SimulationPattern.SimulParameters  simulParametersDefault,
						GONIOMETER_PARAMETERS, //LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
						DISTORTION_PROCESS_CONFIGURATION
				);
				if (DEBUG_LEVEL>1){
					System.out.println("Initiaslizing Goniometer class");
				}
			} else if (DEBUG_LEVEL>1){
				System.out.println("GONIOMETER was initialized");
			}
			if (LENS_DISTORTIONS==null) {
				String msg="LENS_DISTORTIONS is not set";
				System.out.println("Error" + msg);
				IJ.showMessage("Error",msg);
				return;
			}
			if (DISTORTION_CALIBRATION_DATA==null) {
				String msg="DISTORTION_CALIBRATION_DATA is not set";
				System.out.println("Error" + msg);
				IJ.showMessage("Error",msg);
				return;
			}

			GONIOMETER.testHintedTarget (
					CAMERAS.getImages(0), // last acquired images with number of pointers detected>=0
					LENS_DISTORTIONS,            //Distortions lensDistortions, // should not be null
					DISTORTION_CALIBRATION_DATA, // Distortions.DistortionCalibrationData distortionCalibrationData,
					PATTERN_PARAMETERS,          //Distortions.PatternParameters patternParameters, // should not be null
					COMPONENTS.equalizeGreens,   //boolean equalizeGreens,
					THREADS_MAX,                 // int       threadsMax,
					UPDATE_STATUS,               //boolean   updateStatus,
					DEBUG_LEVEL);                // debug level used inside loops
					
					
			return;
		}
/* ======================================================================== */
		if       (label.equals("Simulate Grid View")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			if (LENS_DISTORTIONS.fittingStrategy==null) {
				IJ.showMessage("LENS_DISTORTION.fittingStrategy is not set");
				return;
			}
			int stationNumber=0;
			GenericDialog gd=new GenericDialog ("Parameters for simulating pattern grid on sensor");
			if (LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.numStations>1){
	    		gd.addNumericField("Station number (0.."+(LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.numStations-1), stationNumber, 0);
			}
    		gd.addNumericField("Number of sub-camera (starting with 0)", 0, 0);
    		gd.addNumericField("Camera tilt (0 - vertical, >0 looking above horizon on the target", 0.0, 1,6,"degrees");
    		gd.addNumericField("Camera axial (0 - subcamera 0 looking to the target, >0 - rotated clockwise", 0.0, 1,6,"degrees");
    		gd.showDialog();
    		if (gd.wasCanceled()) return;
			if (LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.numStations>1){
	    		stationNumber=        (int) gd.getNextNumber();
			}
    		int channelNumber=        (int) gd.getNextNumber();
    		double tilt=       gd.getNextNumber();
    		double axial=      gd.getNextNumber();
    		ImagePlus imp_simulatePatternOnSensor=LENS_DISTORTIONS.simulatePatternOnSensor(
    				stationNumber,
    				channelNumber,
    				tilt,
    				axial,
    				SIMUL,
    				THREADS_MAX,                 //int threadsMax,
    				UPDATE_STATUS,               // boolean updateStatus,
    				DISTORTION.loop_debug_level, // int mspDebugLevel,
    				MASTER_DEBUG_LEVEL,          //int global_debug_level, // DEBUG_LEVEL
    				MASTER_DEBUG_LEVEL           //int debug_level // debug level used inside loops
    		);
    		imp_simulatePatternOnSensor.show();	
			return;
		}
/* ======================================================================== */
		if       (label.equals("Show grid/hint")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			if (LENS_DISTORTIONS.fittingStrategy==null) {
				IJ.showMessage("LENS_DISTORTION.fittingStrategy is not set");
				return;
			}
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			LENS_DISTORTIONS.showGridAndHint();
			return;
		}
		
		//"Simulate Grid View"
//		/"Test Progress"
/* ======================================================================== */
		if       (label.equals("Test Progress")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			long 	  endTime=System.nanoTime();
			for (int ii=0;ii<=100; ii++){
				endTime+=100000000;
				while ((endTime-System.nanoTime())>0);
				IJ.showProgress(0.01*ii);
				IJ.showStatus(ii+"%");
			}
			return;
		}
/* ======================================================================== */
		if       (label.equals("Load Pixel Mapping")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			PIXEL_MAPPING=new PixelMapping((String)null,DEBUG_LEVEL);
			return;
		}
/* ======================================================================== */
		if       (label.equals("List Mapping Parameters")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) PIXEL_MAPPING=new PixelMapping((String)null,DEBUG_LEVEL);
			PIXEL_MAPPING.listParameters();
			return;
		}
/* ======================================================================== */
		if       (label.equals("Test Direct Mapping")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) PIXEL_MAPPING=new PixelMapping((String)null,DEBUG_LEVEL);
			GenericDialog gd=new GenericDialog("Select parameters for sensor->equirectangular pixel mapping");
			gd.addNumericField("Channel number (0..."+PIXEL_MAPPING.sensors.length,0,0);
			gd.addNumericField("Output frame width", 2592,0,4,"output pix");
			gd.addNumericField("Output frame height", 1936,0,4,"output pix");
			gd.addNumericField("Resolution scale (2.0 - twice resolution, 0.5 - half)", 1.0, 4,6,"x" );
			gd.addNumericField("Left margin", 1.0, 1,6,"sensor pix" );
			gd.addNumericField("Top margin", 1.0, 1,6,"sensor pix" );
			gd.addCheckbox("Convert to equirectangular (azimuth/elevation), false - to mm in focal plane",true);
			gd.showDialog();
			if (gd.wasCanceled()) return;
			int channelNumber= (int) gd.getNextNumber();
			int outputWidth=   (int) gd.getNextNumber();
			int outputHeight=  (int) gd.getNextNumber();
			double resolutionScale=  gd.getNextNumber();
			double x0=               gd.getNextNumber();
			double y0=               gd.getNextNumber();
			boolean toAngles=        gd.getNextBoolean();
			String [] flatTitles={"X-focal","Y-focal","mask"};
			String [] angleTitles={"Azimuth","Elevation","mask"};
			String [] titles=toAngles?angleTitles:flatTitles;
			PIXEL_MAPPING.sensors[channelNumber].combineDistortionsSensorToEquirectangular(
					outputWidth, //int width,
					outputHeight, //int height,
		    		x0, //double x0,
		    		y0, //double y0,
		    		1.0/resolutionScale, //double pixelStep,
		    		!toAngles); //boolean flat)
			this.SDFA_INSTANCE.showArrays(PIXEL_MAPPING.sensors[channelNumber].directMap.map, outputWidth, outputHeight,  true, "DPM-"+channelNumber, titles);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Test Equirectangular Mapping")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) PIXEL_MAPPING=new PixelMapping((String)null,DEBUG_LEVEL);
			GenericDialog gd=new GenericDialog("Select parameters for equirectangular->sensor pixel mapping");
			gd.addMessage("Equirectangular area");
			gd.addNumericField("Longitude left", -180.0, 1,6,"degrees" );
			gd.addNumericField("Longitude right", 180.0, 1,6,"degrees" );
			gd.addNumericField("Latitude top", 90.0, 1,6,"degrees" );
			gd.addNumericField("Latitude bottom", -90.0, 1,6,"degrees" );
			gd.addNumericField("Pixels horizontal ", 14000,0,5,"image pix");
			
			gd.addMessage("Source image parameters");
			gd.addNumericField("Channel number (0..."+PIXEL_MAPPING.sensors.length,0,0);
			gd.addNumericField("Input image width", 2592,0,4,"image pix");
			gd.addNumericField("Input image height", 1936,0,4,"image pix");
			gd.addNumericField("Input image resolution scale (2.0 - twice resolution, 0.5 - half)", 1.0, 4,6,"x" );
			gd.addNumericField("Input image left margin", 0.0, 1,6,"sensor pix" );
			gd.addNumericField("Input image top margin", 0.0, 1,6,"sensor pix" );

			gd.addCheckbox    ("Show result", false);

			gd.showDialog();
			if (gd.wasCanceled()) return;
			
			double longitudeLeft=       gd.getNextNumber();
			double longitudeRight=      gd.getNextNumber();
			double latitudeTop=         gd.getNextNumber();
			double latitudeBottom=      gd.getNextNumber();
			int pixelsHorizontal= (int) gd.getNextNumber();
			int channelNumber=    (int) gd.getNextNumber();
			int imageWidth=      (int) gd.getNextNumber();
			int imageHeight=     (int) gd.getNextNumber();
			double resolutionScale=     gd.getNextNumber();
			double x0=                  gd.getNextNumber();
			double y0=                  gd.getNextNumber();
			boolean showResult=         gd.getNextBoolean();
			String [] titles={"X-pixel","Y-pixel","mask"};
			PIXEL_MAPPING.sensors[channelNumber].combineDistortionsEquirectangularToSensor(
					channelNumber,
					longitudeLeft,
					longitudeRight,
					latitudeTop,
					latitudeBottom,
					pixelsHorizontal,
					imageWidth, //int width,
					imageHeight, //int height,
		    		x0, //double x0,
		    		y0, //double y0,
		    		1.0/resolutionScale, //double pixelStep,
		    		THREADS_MAX);
			if (showResult) this.SDFA_INSTANCE.showArraysSparse(
					PIXEL_MAPPING.sensors[channelNumber].equirectangularMap.map,
					PIXEL_MAPPING.sensors[channelNumber].equirectangularMap.pixelsHorizontal,
					PIXEL_MAPPING.sensors[channelNumber].equirectangularMap.pixelsVertical,
					true,
					"DPM-"+channelNumber, titles);
			  if (DEBUG_LEVEL>0) System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");

			return;
		}
/* ======================================================================== */
		if       (label.equals("Crop Equirectangular Mapping")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) PIXEL_MAPPING=new PixelMapping((String) null,DEBUG_LEVEL);
			GenericDialog gd=new GenericDialog("Select parameters cropping equirectangular->sensor pixel mapping");
			gd.addNumericField("Channel number (0..."+PIXEL_MAPPING.sensors.length,0,0);
			gd.addNumericField("Map width",    2400,0,4,"image pix");
			gd.addCheckbox    ("Clear full map",    false);
			gd.showDialog();
			if (gd.wasCanceled()) return;
			int channelNumber=    (int) gd.getNextNumber();
			int outputWidth=      (int) gd.getNextNumber();
			boolean clearFullMap=gd.getNextBoolean();
			if (PIXEL_MAPPING.sensors[channelNumber].equirectangularMap.map==null){
				IJ.showMessage("Equirectangular map for channel "+channelNumber+" does not exist.");
				return;
			}
			int numPixels=PIXEL_MAPPING.sensors[channelNumber].equirectangularMap.createPartialMap(outputWidth);
			String [] titles={"X-pixel","Y-pixel","mask"};
			this.SDFA_INSTANCE.showArrays(PIXEL_MAPPING.sensors[channelNumber].equirectangularMap.partialMap,
					PIXEL_MAPPING.sensors[channelNumber].equirectangularMap.mapWOI.width,
					PIXEL_MAPPING.sensors[channelNumber].equirectangularMap.mapWOI.height,
					true,
					"DPM-"+channelNumber, titles);
			  if (DEBUG_LEVEL>0) System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");
			  if (DEBUG_LEVEL>0) System.out.println ("Number of defined pixels ="+numPixels+" ("+
					  IJ.d2s(100.0*numPixels/PIXEL_MAPPING.sensors[channelNumber].equirectangularMap.partialMap[0].length)+"%)");
			  if (clearFullMap){
				  PIXEL_MAPPING.sensors[channelNumber].equirectangularMap.map=null;
				  if (DEBUG_LEVEL>0) System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");
			  }
			return;
		}
/* ======================================================================== */
		if       (label.equals("Generate & Save Equirectangular")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) PIXEL_MAPPING=new PixelMapping((String)null,DEBUG_LEVEL); // ask for and load sensor calibration files
			String [] extensions={".eqr-tiff",".eqrect-tiff"};
			CalibrationFileManagement.MultipleExtensionsFileFilter parFilter = new CalibrationFileManagement.MultipleExtensionsFileFilter("",extensions,"equirectangular map "+extensions[0]+" files");
			String pathname=CalibrationFileManagement.selectFile(true,
					"Save equirectangular calibration file (will add channel number suffix)",
					"Save",
					parFilter,
					""); //String defaultPath
			if ((pathname==null) || (pathname=="")) return;
			GenericDialog gd=new GenericDialog("Select parameters for equirectangular->sensor pixel mapping");
			gd.addMessage("Equirectangular area");
			gd.addNumericField("Longitude left", -180.0, 1,6,"degrees" );
			gd.addNumericField("Longitude right", 180.0, 1,6,"degrees" );
			gd.addNumericField("Latitude top", 90.0, 1,6,"degrees" );
			gd.addNumericField("Latitude bottom", -90.0, 1,6,"degrees" );
			gd.addNumericField("Pixels horizontal ", 14268,0,5,"image pix");
			
			gd.addMessage("Source image parameters");
			gd.addNumericField("Input image width", 2592,0,4,"image pix");
			gd.addNumericField("Input image height", 1936,0,4,"image pix");
			gd.addNumericField("Input image resolution scale (2.0 - twice resolution, 0.5 - half)", 1.0, 4,6,"x" );
			gd.addNumericField("Input image left margin", 0.0, 1,6,"sensor pix" );
			gd.addNumericField("Input image top margin", 0.0, 1,6,"sensor pix" );
			gd.addMessage("Reduction of memory usage");
			gd.addNumericField("Crop files horizontally to ", 3000,0,4,"longitude pix");
			gd.addCheckbox    ("Clear full map",    true);
			gd.addCheckbox    ("Clear all data",    true);
			gd.showDialog();
			if (gd.wasCanceled()) return;
			
			double longitudeLeft=       gd.getNextNumber();
			double longitudeRight=      gd.getNextNumber();
			double latitudeTop=         gd.getNextNumber();
			double latitudeBottom=      gd.getNextNumber();
			int pixelsHorizontal= (int) gd.getNextNumber();
			int imageWidth=      (int) gd.getNextNumber();
			int imageHeight=     (int) gd.getNextNumber();
			double resolutionScale=     gd.getNextNumber();
			double x0=                  gd.getNextNumber();
			double y0=                  gd.getNextNumber();
			int longitudeWidth=   (int) gd.getNextNumber();
			boolean clearFullMap=gd.getNextBoolean();
			boolean clearAllMaps=gd.getNextBoolean();

			PIXEL_MAPPING.generateAndSaveEquirectangularMaps(
					pathname,
					longitudeLeft,
					longitudeRight,
					latitudeTop,
					latitudeBottom,
					pixelsHorizontal,
					imageWidth, //int width,
					imageHeight, //int height,
		    		x0, //double x0,
		    		y0, //double y0,
		    		1.0/resolutionScale, //double pixelStep,
		    		longitudeWidth,
		    		clearFullMap,
		    		clearAllMaps,
		    		THREADS_MAX);
			  if (DEBUG_LEVEL>0) System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");

			return;
		}
/* ======================================================================== */
		if       (label.equals("Load Equirectangular Maps")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null){
				if (DEBUG_LEVEL>1) System.out.println("Creating new PixelMapping from equirectangular maps");
				PIXEL_MAPPING=new PixelMapping((String)null,DEBUG_LEVEL); // ask for and load sensor calibration files
//				PIXEL_MAPPING=new PixelMapping(null,true,DEBUG_LEVEL); // will not ask for sensor calibration files
				PIXEL_MAPPING.loadChannelMaps (null); // ask
			}
			else {
				if (DEBUG_LEVEL>1) System.out.println("Loadidng maps, preserving sensor data");
				PIXEL_MAPPING.loadChannelMaps (null); // ask
			}
			return;
		}
/* ======================================================================== */
		if       (label.equals("Show Maps Overlap")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) PIXEL_MAPPING=new PixelMapping(null,true,DEBUG_LEVEL);
			float [] pixels= PIXEL_MAPPING.generateOverlapMap();
			this.SDFA_INSTANCE.showArrays(pixels,
					PIXEL_MAPPING.getPanoWidth(),
					PIXEL_MAPPING.getPanoHeight(),
					"Overlap");
			GenericDialog gd=new GenericDialog("Overlap alpha threshold");
			gd.addNumericField("pixel alpha threshold ",50,1,5,"%");
			gd.showDialog();
			if (gd.wasCanceled()) return;
			double threshold=0.01*gd.getNextNumber();
			double [][] overlap=PIXEL_MAPPING.overlapPairsAreas(threshold);
			PIXEL_MAPPING.listOverlap(overlap);
			for (int i=0;i<(overlap.length-1);i++) for (int j=i+1;j<overlap.length;j++){
				overlap[i][j]/=Math.sqrt(overlap[i][i]*overlap[j][j])/100.0;
			}
			PIXEL_MAPPING.listOverlap(overlap);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Test Lanczos")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) PIXEL_MAPPING=new PixelMapping(null,true,DEBUG_LEVEL);
			GenericDialog gd=new GenericDialog("Test Lanczos Kernel stack");
			gd.addNumericField("Lanczos A",  3, 0 );
			gd.addNumericField("Oversampled image",  2, 0 );
			gd.addNumericField("Bins per half-pixel", 50, 0 );
			gd.showDialog();
			if (gd.wasCanceled()) return;
    		int lacroszosA       = (int)  gd.getNextNumber();
    		int oversampled      = (int)  gd.getNextNumber();
    		int binsPerHalfPixel = (int)  gd.getNextNumber();

			double [][][][] lanczos=PIXEL_MAPPING.generateLanczosStack(
					lacroszosA,
					oversampled,
					binsPerHalfPixel);
			//   public double [] testLanczosCenter(double [][][][] lanczos){
		    double [] lanczosTestCenter= PIXEL_MAPPING.testLanczosCenter(lanczos);
			int sizeCenter=(int) Math.sqrt(lanczosTestCenter.length);
			this.SDFA_INSTANCE.showArrays(
					lanczosTestCenter,
					sizeCenter,
					sizeCenter,
					"LanczosCenter_a"+lacroszosA+"_"+oversampled+"_"+binsPerHalfPixel);

		    double [] lanczosTest0= PIXEL_MAPPING.testLanczosStack(lanczos);
			int size=(int) Math.sqrt(lanczosTest0.length);
			this.SDFA_INSTANCE.showArrays(
					lanczosTest0,
					size,
					size,
					"Lanczos_a"+lacroszosA+"_"+oversampled+"_"+binsPerHalfPixel);
		    double [] lanczosTestNormalize0= PIXEL_MAPPING.testLanczosStackNormalization(lanczos);
			int sizeNorm=(int) Math.sqrt(lanczosTestNormalize0.length);
			this.SDFA_INSTANCE.showArrays(
					lanczosTestNormalize0,
					sizeNorm,
					sizeNorm,
					"Lanczos-norm_a"+lacroszosA+"_"+oversampled+"_"+binsPerHalfPixel);
			PIXEL_MAPPING.normalizeLanczosStack(lanczos);
			
		    double [] lanczosTest1= PIXEL_MAPPING.testLanczosStack(lanczos);
			this.SDFA_INSTANCE.showArrays(
					lanczosTest1,
					size,
					size,
					"Lanczos1_a"+lacroszosA+"_"+oversampled+"_"+binsPerHalfPixel);
		    double [] lanczosTestNormalize1= PIXEL_MAPPING.testLanczosStackNormalization(lanczos);
			this.SDFA_INSTANCE.showArrays(
					lanczosTestNormalize1,
					sizeNorm,
					sizeNorm,
					"Lanczos1-norm_a"+lacroszosA+"_"+oversampled+"_"+binsPerHalfPixel);
			return;
		}

/* ======================================================================== */
		if       (label.equals("Warp Image")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) PIXEL_MAPPING=new PixelMapping(null,true,DEBUG_LEVEL);
			imp_sel = WindowManager.getCurrentImage();
			if (imp_sel==null) {
				IJ.showMessage("No image selected");
				return;
			}
			GenericDialog gd=new GenericDialog("Test warping images");
			gd.addNumericField("Channel number (0..."+PIXEL_MAPPING.sensors.length,0,0);
			gd.addNumericField("Source image pixel scale",2,0,1,"x");
			gd.showDialog();
			if (gd.wasCanceled()) return;
    		int channel          = (int)  gd.getNextNumber();
    		int sourceImageScale = (int)  gd.getNextNumber();
			
			ImagePlus imp_warped= PIXEL_MAPPING.resampleToEquirectangular(
					imp_sel,
					channel,
					sourceImageScale,
					THREADS_MAX);
			imp_warped.show();
			return;
		}
 
/* ======================================================================== */
		if       (label.equals("Warp Files")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) PIXEL_MAPPING=new PixelMapping(null,true,DEBUG_LEVEL);
	    	String [] extensions={".tiff",".jpeg"};
	    	CalibrationFileManagement.MultipleExtensionsFileFilter parFilter =
	    		new CalibrationFileManagement.MultipleExtensionsFileFilter("",extensions,"channel images, may be scaled");
	    	String [] files=CalibrationFileManagement.selectFiles(false,
	    			"Select channel images to warp to equirectangular",
	    			"Select",
	    			parFilter,
	    			null); // String [] defaultPaths);
	    	if ((files==null) || (files.length==0)) {
	    		IJ.showMessage("No files selected");
	    		return;
	    	}
			GenericDialog gd=new GenericDialog("Test warping images");
			gd.addNumericField("Source image pixel scale",2,0,1,"x");
			gd.addCheckbox("Show warped images", false);
			gd.addCheckbox("Save warped images", true);
			gd.showDialog();
			if (gd.wasCanceled()) return;
    		int sourceImageScale = (int)  gd.getNextNumber();
    		boolean showResults =gd.getNextBoolean();
    		boolean saveResults =gd.getNextBoolean();
			PIXEL_MAPPING.resampleToEquirectangular(
					files,
					sourceImageScale,
		    		showResults,
		    		saveResults,
					THREADS_MAX);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Create Plane Map")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) PIXEL_MAPPING=new PixelMapping(null,true,DEBUG_LEVEL);
//    		if (useSelectedChannels && (ABERRATIONS_PARAMETERS!=null))selectedChannels=ABERRATIONS_PARAMETERS.selectedChannels;
			boolean [] selectedChannels=ABERRATIONS_PARAMETERS.selectedChannels;
			int numUsedChannels=0;
			int [] usedChannels;
			String [] extensions={".plane-tiff",".intermap.tiff"};
			String extList="";
			for (int i=0;i<extensions.length;i++){
				if (i>0) extList+=", ";
				extList+="*"+extensions[i];
			}
			if (selectedChannels!=null) for (int i=0;i<selectedChannels.length;i++) if (selectedChannels[i]) numUsedChannels++;
			if (numUsedChannels==0){
				IJ.showMessage("Error","No channels selected to create maps, using all loaded");
				numUsedChannels=PIXEL_MAPPING.sensors.length;
				usedChannels=new int [numUsedChannels];
				for (int i=0;i<numUsedChannels;i++)usedChannels[i]=i;
			} else {
				usedChannels=new int [numUsedChannels];
				int index=0;
				for (int i=0;i<selectedChannels.length;i++) if (selectedChannels[i]) usedChannels[index++]=i;
			}
			if (numUsedChannels==0){
				IJ.showMessage("Error","No channels to create maps");
				return;
			}
			boolean strictHorizontalDisparity=true;
			int channelRight = 1;
			int channelLeft  = 2;
			double nominalHorizontalDisparity=60.0; // mm
			GenericDialog gd;
			if (numUsedChannels>=3) {
				gd=new GenericDialog("Select generation of sensor pixel maps to the common plane mode");
				gd.addMessage     ("Select projection plane rotation - align to a pair of horizontallty spaced cameras or minimize roll");
				gd.addCheckbox    ("Strict horizontal disparity (unchecked - minimize view roll)", strictHorizontalDisparity);
				gd.addNumericField("Right camera channel number (0..2)",channelRight,0);
				gd.addNumericField("Left  camera channel number (0..2)",channelLeft, 0);
				gd.addNumericField("Nominal distance between left and right cameras",nominalHorizontalDisparity,2,6,"mm");
				gd.showDialog();
				if (gd.wasCanceled()) return;
				strictHorizontalDisparity=gd.getNextBoolean();
				channelRight=  (int) gd.getNextNumber();
				channelLeft=   (int) gd.getNextNumber();
				nominalHorizontalDisparity=     gd.getNextNumber();
			}
			gd=new GenericDialog("Generation of sensor pixel maps to the common plane");
			String sChannels="";
			for (int i=0;i<usedChannels.length;i++) sChannels+=" "+usedChannels[i];
			gd.addMessage("Using camera channels: "+sChannels);
			double [] suggestedYER={0.0,0.0,0.0};
			if (numUsedChannels>=3) suggestedYER=PIXEL_MAPPING.suggestProjectionPlaneYawElevationRoll(
					usedChannels,
					channelRight, //1, //int channel1,
					channelLeft,  //2, //int channel2,
					strictHorizontalDisparity,
		    		DEBUG_LEVEL
				);
			double projectionElevation=     suggestedYER[1];
			double projectionYaw=           suggestedYER[0];
			double projectionRoll=          suggestedYER[2];
			double projectionPixelSize=0.00044036902;
			int    projectionWidth=   2920;
			int     projectionHeight= 2220;
			double projectionCenterX=       0.5*projectionWidth;
			double projectionCenterY=       0.5*projectionHeight;
			
			gd.addStringField ("Result map filename ("+extList+")", "",100);
			gd.addCheckbox    ("Select result file", false);
			gd.addNumericField("View axis elevation (orthogonal to projection plane)",projectionElevation,5,9,"degrees");
			gd.addNumericField("View axis heading   (orthogonal to projection plane)",projectionYaw,5,9,"degrees");
			gd.addNumericField("View plane rotation (roll) around the view axis",     projectionRoll,5,9,"degrees");
			gd.addNumericField("Projection pixel size (relative)     ",1000*projectionPixelSize,5,9,"x1/1000");
			gd.addNumericField("Projection plane width", projectionWidth,0,5,"pix");
			gd.addNumericField("Projection plane height",projectionHeight,0,5,"pix");
			gd.addNumericField("Projection plane Center X (point orthogonal to the view axis), right",   projectionCenterX,3,9,"pix");
			gd.addNumericField("Projection plane Center Y (point orthogonal to the view axis), down",    projectionCenterY,3,9,"pix");
			gd.showDialog();
			if (gd.wasCanceled()) return;
			String resultPath=gd.getNextString();
			if (gd.getNextBoolean() || (resultPath==null) || (resultPath.length()==0)){
				CalibrationFileManagement.MultipleExtensionsFileFilter filter= new CalibrationFileManagement.MultipleExtensionsFileFilter("",extensions,"sensors map file ("+extList+")");
				resultPath=CalibrationFileManagement.selectFile(
						false, // smart
						true,  // save
						"Select sensors map file",
						"Save",
						filter,
						null); // String [] defaultPaths);
		       	if ((resultPath==null) || (resultPath.length()==0)) return;
			}
			projectionElevation=     gd.getNextNumber();
			projectionYaw=           gd.getNextNumber();
			projectionRoll=          gd.getNextNumber();
			projectionPixelSize=      0.001*gd.getNextNumber();
			projectionWidth=   (int) gd.getNextNumber();
			projectionHeight=  (int) gd.getNextNumber();
			projectionCenterX=       gd.getNextNumber();
			projectionCenterY=       gd.getNextNumber();
			String title="Projection_plane_map";

			ImagePlus imp_pixelmap= PIXEL_MAPPING.getPlaneToSensorsMap(
					usedChannels, // int [] channels,
		    		projectionElevation, // Latitude (in degrees) of the normal to the projection plane
		    		projectionYaw,       // Longitude (in degrees) of the normal to the projection plane
		    		projectionRoll,      // Rotation of the projection plane around the perpendicular from the lens centers
		    		projectionPixelSize, // Ratio of the plane pixel size to the distance from the lens center to the projection plane
		    		projectionWidth,     // Width of the projection rectangle
		    		projectionHeight,    // Height of the projection rectangle
		    		projectionCenterX,   // X-coordinate (along the projection plane X - right) of the intersection of the projection plane with the perpendicular from the lens center
		    		projectionCenterY,   // Y-coordinate (along the projection plane Y - down) of the intersection of the projection plane with the perpendicular from the lens center
		    		nominalHorizontalDisparity, // nominal distance between horizontal cameras
		    		title,
		    		DEBUG_LEVEL
		    		);
        	if (imp_pixelmap!=null) {
        		if (DEBUG_LEVEL>2) {
        			imp_pixelmap.getProcessor().resetMinAndMax(); // imp_psf will be reused
        			imp_pixelmap.show();
        		}
        		FileSaver fs=new FileSaver(imp_pixelmap);
        		String msg="Saving pixel map to a plane for sensors "+sChannels+": "+resultPath;
        		if (UPDATE_STATUS) IJ.showStatus(msg);
        		if (DEBUG_LEVEL>0) System.out.println(msg);
        		fs.saveAsTiffStack(resultPath);
        	} else {
           	 System.out.println("Failed to create pixel map for sensors "+sChannels);
        	}
			return;
		}
		
		
	    /**
	     * Create 3xNumber of sensors (channels) map from plane pixels to individual sensor images pixels (pixel-X, pixel-Y, alpha(mask) )
	     * Uses intermediate equirectangular map that should be defined
	     * @param channels - array of sensor numbers
	     * @param projectionElevation - Latitude (in degrees) of the normal to the projection plane
	     * @param projectionYaw - Longitude (in degrees) of the normal to the projection plane
	     * @param projectionRoll - Rotation of the projection plane around the perpendicular from the lens centers
	     * @param projectionPixelSize ratio of the plane pixel size to the distance from the lens center to the projection plane
	     * @param projectionWidth - width of the projection rectangle 
	     * @param projectionHeight - height of the projection rectangle 
	     * @param projectionCenterX - X-coordinate (along the projection plane X - right) of the intersection of the projection plane with the perpendicular from the lens center  
	     * @param projectionCenterY - Y-coordinate (along the projection plane Y - down) of the intersection of the projection plane with the perpendicular from the lens center  
	     * @param title - Image title
	     * @param debugLevel - debug level
	     * @return image with 3*N slices - {pixelX1, pixelY1, alpha1, pixelX2, pixelY2, alpha2, ...} and metadata needed to map images 
	     */
		
//		new CalibrationFileManagement.MultipleExtensionsFileFilter("grid",extensions,"Calibrated grid files");
/* ======================================================================== */
		if       (label.equals("Test Plane Map")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) PIXEL_MAPPING=new PixelMapping(null,true,DEBUG_LEVEL);
			GenericDialog gd=new GenericDialog("Testing of inter-sensor overlap maps");
			gd.addStringField ("Map file", "/home/andrey/goniometer/triclope01/calibration/plane_maps/plane04.plane-tiff",100);
//			gd.addStringField ("Result directory", "/home/andrey/goniometer/corrections01/equirectangular_maps/test1",80);
			gd.addStringField ("Result directory", "/home/andrey/goniometer/triclope01/correction01/results02/plane01",100);
//			gd.addStringField ("Source format", "/home/andrey/goniometer/corrections01/equirectangular_maps/1327916658_133555-%02d-DECONV-RGB24.tiff",80);
			gd.addStringField ("Source format", "/home/andrey/goniometer/triclope01/correction01/results02/1350200009_096256-%02d-DECONV-RGB24.tiff",100);
			gd.showDialog();
			if (gd.wasCanceled()) return;
			String mapFile=         gd.getNextString();
			String resultDirectory= gd.getNextString();
			String fileNameFormat=  gd.getNextString();
			PIXEL_MAPPING.debugLevel= DEBUG_LEVEL;
			PIXEL_MAPPING.applyPlaneMap(
		    		mapFile, //String path, //
		    		fileNameFormat, //String imgPathFormat,
		    		resultDirectory, //String resultDirectory,
		    		2, //int sourceImageScale, // 2.0
		    		true, // boolean saveTiff,
					THREADS_MAX,
					DEBUG_LEVEL);
			return;
		}
		
/* ======================================================================== */
		if       (label.equals("Create Intermaps")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) PIXEL_MAPPING=new PixelMapping(null,true,DEBUG_LEVEL);
			GenericDialog gd=new GenericDialog("Generation of inter-sensor overlap maps");
			gd.addStringField ("Result directory", "",60);
			gd.addStringField ("Result maps format",  "intersensor-%02d-%02d.intermap.tiff",40);
			gd.addNumericField("Pixel alpha threshold",50,1,5,"%");
			gd.addNumericField("Minimal overlap area", 2.0,1,5,"%");
			gd.addNumericField("Scan overlap angle, direction between sensors", 20.0,1,5,"degrees");
			gd.addNumericField("Scan overlap angle, direction perpendicular to the sensors line", 70.0,1,5,"degrees");
			gd.addNumericField("Overlap extra pixels in inter-sensor direction", 200,0,4,"pixels");
			gd.showDialog();
			if (gd.wasCanceled()) return;
			String directory=gd.getNextString();
			String fileNameFormat=gd.getNextString();
			double alphaThreshold=     0.01*gd.getNextNumber();
			double overlapThreshold=   0.01*gd.getNextNumber();
			double angleSizeX=              gd.getNextNumber();
			double angleSizeY=              gd.getNextNumber();
			int    overlapExtraPixels=(int) gd.getNextNumber();
			PIXEL_MAPPING.debugLevel= DEBUG_LEVEL;
			PIXEL_MAPPING.generateAndSaveInterSensorMaps(
		    		directory, 
		    		fileNameFormat, //prefix%02d-%02dsuffix
		    		alphaThreshold,
		    		overlapThreshold,
		    		angleSizeX,
		    		angleSizeY,
		    		overlapExtraPixels,
		    		UPDATE_STATUS,
		    		DEBUG_LEVEL
		    		);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Test Intermaps")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) PIXEL_MAPPING=new PixelMapping(null,true,DEBUG_LEVEL);
			GenericDialog gd=new GenericDialog("Testing of inter-sensor overlap maps");
			gd.addStringField ("Map file", "/home/andrey/goniometer/corrections01/equirectangular_maps/intersensor-designedlocations-01-09.intermap.tiff",80);
//			gd.addStringField ("Result directory", "/home/andrey/goniometer/corrections01/equirectangular_maps/test1",80);
			gd.addStringField ("Result directory", "",80);
			gd.addStringField ("Source format", "/home/andrey/goniometer/corrections01/equirectangular_maps/1327916658_133555-%02d-DECONV-RGB24.tiff",80);
			gd.showDialog();
			if (gd.wasCanceled()) return;
			String mapFile=         gd.getNextString();
			String resultDirectory= gd.getNextString();
			String fileNameFormat=  gd.getNextString();
			PIXEL_MAPPING.debugLevel= DEBUG_LEVEL;
			double [][] overlapStack=PIXEL_MAPPING.applyOverlapMap(
					mapFile,
					fileNameFormat,
					resultDirectory,
					2, // image scale
		    		true, //boolean saveTiff,
		    		true, //boolean convertToDouble,
					THREADS_MAX,
					DEBUG_LEVEL);
			String [] titles={"Alpha1","Y1","Cb1","Cr1","Alpha2","Y2","Cb2","Cr2"};
			this.SDFA_INSTANCE.showArrays(
					overlapStack,
					PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
					PIXEL_MAPPING.lastUsedInterSensor.mapHeight,
					true,
			"overlapAYCbCr",
			titles);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Show Sobel")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) {
				String msg="PIXEL_MAPPING is not initialized";
				System.out.println(msg);
				IJ.showMessage(msg);
				return;
			}
			if ((PIXEL_MAPPING.lastUsedInterSensor==null) || (PIXEL_MAPPING.lastUsedInterSensor.overlapImages==null)) {
				String msg="Inter-sensor data is not initialized (use \"Test Intermaps\"";
				System.out.println(msg);
				IJ.showMessage(msg);
				return;
			}
			GenericDialog gd=new GenericDialog("Edge detection");
			gd.addNumericField("Blur sigma",PIXEL_MAPPING.lastUsedInterSensor.preSobelSigma,3,5,"pix");
			gd.addNumericField("Spread portion of edge pixel in vertical/horizontal directions",PIXEL_MAPPING.lastUsedInterSensor.edgeSpreadOrtho,3,5,"pix");
			gd.addNumericField("Spread portion of edge pixel in diagonal directions",PIXEL_MAPPING.lastUsedInterSensor.edgeSpreadDiagonal,3,5,"pix");
			gd.addCheckbox("Threshold edges", true);
			gd.addNumericField("Threshold high",PIXEL_MAPPING.lastUsedInterSensor.edgeThresholdHigh,3,5,"pix");
			gd.addNumericField("Threshold low",PIXEL_MAPPING.lastUsedInterSensor.edgeThresholdLow,3,5,"pix");
			gd.addCheckbox("Apply alpha", true);
			gd.addCheckbox("Show progress", true);
			gd.showDialog();
			if (gd.wasCanceled()) return;
			PIXEL_MAPPING.lastUsedInterSensor.preSobelSigma=      gd.getNextNumber();
			PIXEL_MAPPING.lastUsedInterSensor.edgeSpreadOrtho=    gd.getNextNumber();
			PIXEL_MAPPING.lastUsedInterSensor.edgeSpreadDiagonal= gd.getNextNumber();
			boolean thresholdEdges=                               gd.getNextBoolean();
			PIXEL_MAPPING.lastUsedInterSensor.edgeThresholdHigh=gd.getNextNumber();
			PIXEL_MAPPING.lastUsedInterSensor.edgeThresholdLow=gd.getNextNumber();
			
			PIXEL_MAPPING.lastUsedInterSensor.initSobelY(true); // recalculate even if it existed
			
			boolean applyAlpha=gd.getNextBoolean();
			boolean showProgress=gd.getNextBoolean();
			String [] titles={"Left","Right"};
			this.SDFA_INSTANCE.showArrays(
					PIXEL_MAPPING.lastUsedInterSensor.sobelY,
					PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
					PIXEL_MAPPING.lastUsedInterSensor.mapHeight,
					true,
			"Sobel"+((""+PIXEL_MAPPING.lastUsedInterSensor.preSobelSigma).replace('.','_')),
			titles);
			if (!thresholdEdges) return;
			float [][]fEdges=new float [PIXEL_MAPPING.lastUsedInterSensor.sobelY.length][];
			PIXEL_MAPPING.lastUsedInterSensor.booleanEdges=new boolean[2][];
			for (int side=0;side<PIXEL_MAPPING.lastUsedInterSensor.sobelY.length;side++) {
				if (PIXEL_MAPPING.lastUsedInterSensor.sobelY[side]!=null){
					//    				preImage=this.overlapImages[1+this.overlapImages.length/2];

					PIXEL_MAPPING.lastUsedInterSensor.booleanEdges[side]=PIXEL_MAPPING.lastUsedInterSensor.thresholdEdges(
							PIXEL_MAPPING.lastUsedInterSensor.sobelY[side],
							applyAlpha?PIXEL_MAPPING.lastUsedInterSensor.overlapImages[side*(PIXEL_MAPPING.lastUsedInterSensor.overlapImages.length/2)]:null, //alpha
							PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
							PIXEL_MAPPING.lastUsedInterSensor.edgeThresholdHigh,
							PIXEL_MAPPING.lastUsedInterSensor.edgeThresholdLow,
							showProgress,
							DEBUG_LEVEL);
					fEdges[side]=new float[PIXEL_MAPPING.lastUsedInterSensor.booleanEdges[side].length];
					for (int i=0;i<PIXEL_MAPPING.lastUsedInterSensor.booleanEdges[side].length;i++)fEdges[side][i]=PIXEL_MAPPING.lastUsedInterSensor.booleanEdges[side][i]?1.0f:0.0f;
				} else fEdges[side]=null;
			}
			this.SDFA_INSTANCE.showArrays(
					fEdges,
					PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
					PIXEL_MAPPING.lastUsedInterSensor.mapHeight,
					true,
			"Edges-TH"+((""+PIXEL_MAPPING.lastUsedInterSensor.edgeThresholdHigh).replace('.','_'))+"-TL"+((""+PIXEL_MAPPING.lastUsedInterSensor.edgeThresholdLow).replace('.','_')),
			titles);
			return;
//		addButton("Show Sobel",panelStereo);
		}	
/* ======================================================================== */
		if       (label.equals("Edge Thinning")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) {
				String msg="PIXEL_MAPPING is not initialized";
				System.out.println(msg);
				IJ.showMessage(msg);
				return;
			}
			if ((PIXEL_MAPPING.lastUsedInterSensor==null) || (PIXEL_MAPPING.lastUsedInterSensor.booleanEdges==null)) {
				String msg="booleanEdges is not initialized (use \"Show Sobel\"";
				System.out.println(msg);
				IJ.showMessage(msg);
				return;
			}
			String [] choices={"left image","right image","both images"};
			GenericDialog gd=new GenericDialog("Edge Thinning Parameters");
			gd.addChoice("Thin edges for",choices, choices[2]);

			gd.addNumericField("Minimal cycle area",40,0,5,"pix");
			gd.addNumericField("Ends test radius",50,0,5,"pix");
			gd.addNumericField("Minimal alpha",     0.5,3,5,"");
			gd.showDialog();
			if (gd.wasCanceled()) return;
			int sides=               gd.getNextChoiceIndex()+1; // 1-2-3
			int minCycleArea= (int) gd.getNextNumber();
			int endRadius= (int) gd.getNextNumber();
			double minAlpha=gd.getNextNumber();
			PIXEL_MAPPING.lastUsedInterSensor.thinEdges=new boolean[2][];
			for (int side=0;side<2;side++) if ((sides & (1<<side))!=0){
				PIXEL_MAPPING.lastUsedInterSensor.thinEdges[side]=PIXEL_MAPPING.lastUsedInterSensor.EdgeThinning(
						PIXEL_MAPPING.lastUsedInterSensor.sobelY[side],
						PIXEL_MAPPING.lastUsedInterSensor.booleanEdges[side],
						PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
						minCycleArea, // minimal area of the closed loop to keep
						endRadius,
						UPDATE_STATUS,
						DEBUG_LEVEL
				);
				if (minAlpha>0){
					double [] alphaLayer=PIXEL_MAPPING.lastUsedInterSensor.overlapImages[side*(PIXEL_MAPPING.lastUsedInterSensor.overlapImages.length/2)];
					for (int i=0;i<PIXEL_MAPPING.lastUsedInterSensor.thinEdges[side].length;i++)
						if (alphaLayer[i]<minAlpha)PIXEL_MAPPING.lastUsedInterSensor.thinEdges[side][i]=false;
				}
			}
			float [][]fEdges=new float [PIXEL_MAPPING.lastUsedInterSensor.sobelY.length][];
			String [] titles={"Left","Right"};
			for (int side=0;side<2;side++) if ((sides & (1<<side))!=0){
				fEdges[side]=new float[PIXEL_MAPPING.lastUsedInterSensor.thinEdges[side].length];
				for (int i=0;i<PIXEL_MAPPING.lastUsedInterSensor.thinEdges[side].length;i++)
					fEdges[side][i]=PIXEL_MAPPING.lastUsedInterSensor.thinEdges[side][i]?1.0f:0.0f;
			} else fEdges[side]=null;
			this.SDFA_INSTANCE.showArrays(
					fEdges,
					PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
					PIXEL_MAPPING.lastUsedInterSensor.mapHeight,
					true,
					"ThinEdges-MC"+minCycleArea+"-ER"+endRadius,
					titles);
			return;
		}	
/* ======================================================================== */
		if       (label.equals("Distance from Edges")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) {
				String msg="PIXEL_MAPPING is not initialized";
				System.out.println(msg);
				IJ.showMessage(msg);
				return;
			}
			if ((PIXEL_MAPPING.lastUsedInterSensor==null) || (PIXEL_MAPPING.lastUsedInterSensor.thinEdges==null)) {
				String msg="thinEdges is not initialized (use \"Edge Thinning\"";
				System.out.println(msg);
				IJ.showMessage(msg);
				return;
			}
			String [] choices={"left image","right image","both images"};
			GenericDialog gd=new GenericDialog("Edge Thinning Parameters");
			gd.addChoice("Distance from Edges for",choices, choices[2]);
			gd.addCheckbox("Margins as edges", false);
			gd.addNumericField("Alpha threshold",50.0,1,5,"%");
			gd.addNumericField("Number of passes",2,0,5,"");
			gd.addNumericField("Directions mask",63,0,5,"pix");
			

			gd.showDialog();
			if (gd.wasCanceled()) return;
			int sides=               gd.getNextChoiceIndex()+1; // 1-2-3
			boolean marginsAsEdges=gd.getNextBoolean();
			double alphaThreshold= 0.01*gd.getNextNumber();
			int numPasses=  (int) gd.getNextNumber();
			int directions= (int) gd.getNextNumber();
			PIXEL_MAPPING.lastUsedInterSensor.distanceFromEdges =new float [2][];
			for (int side=0;side<2;side++) if (((sides & (1<<side))!=0) && (PIXEL_MAPPING.lastUsedInterSensor.thinEdges[side]!=null)){
				boolean [] alphaMask= null;
				if (marginsAsEdges){
					alphaMask=new boolean[PIXEL_MAPPING.lastUsedInterSensor.thinEdges[side].length];
					for (int i=0;i<alphaMask.length;i++)alphaMask[i]=
						PIXEL_MAPPING.lastUsedInterSensor.overlapImages[side*(PIXEL_MAPPING.lastUsedInterSensor.overlapImages.length/2)][i]>alphaThreshold;
				}
				PIXEL_MAPPING.lastUsedInterSensor.distanceFromEdges[side]=PIXEL_MAPPING.lastUsedInterSensor.edgesSegmeniting.distanceFromEdges(
						PIXEL_MAPPING.lastUsedInterSensor.thinEdges[side],
						alphaMask,
						PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
						numPasses,
						directions,
		    			DEBUG_LEVEL);
			}
			String [] titles={"Left","Right"};
			this.SDFA_INSTANCE.showArrays(
					PIXEL_MAPPING.lastUsedInterSensor.distanceFromEdges,
					PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
					PIXEL_MAPPING.lastUsedInterSensor.mapHeight,
					true,
					"DistanceFromEdges",
					titles);
			return;
		}	
		/* ======================================================================== */
		if       (label.equals("Edge Areas")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) {
				String msg="PIXEL_MAPPING is not initialized";
				System.out.println(msg);
				IJ.showMessage(msg);
				return;
			}
			if ((PIXEL_MAPPING.lastUsedInterSensor==null) || (PIXEL_MAPPING.lastUsedInterSensor.distanceFromEdges==null)) {
				String msg="distanceFromEdges is not initialized (use \"Distance from Edges\"";
				System.out.println(msg);
				IJ.showMessage(msg);
				return;
			}
			String [] choices={"left image","right image","both images"};
			GenericDialog gd=new GenericDialog("Edge Thinning Parameters");
			gd.addChoice      ("Edge Areas for",choices, choices[2]);
			gd.addCheckbox    ("Gaps not Edges", false);
			gd.addNumericField("Initial distance threshold (higher)",10.0,1,5,"");
			gd.addNumericField("Final   distance threshold (lower)", 0.0,1,5,"");
			gd.addNumericField("Distance threshold step between passes)", 1.0,1,5,"");
			gd.addCheckbox    ("Repeat each step", false);
			gd.addCheckbox    ("Repeat last step", false);
			
			gd.showDialog();
			if (gd.wasCanceled()) return;
			int sides=               gd.getNextChoiceIndex()+1; // 1-2-3
			boolean gapsNotEdges=         gd.getNextBoolean();
			double thresholdStart=        gd.getNextNumber();
			double thresholdFinal=        gd.getNextNumber();
			double thresholdStep=         gd.getNextNumber();
			boolean repeatEachStep=       gd.getNextBoolean();
			boolean repeatLastStep=       gd.getNextBoolean();
			PIXEL_MAPPING.lastUsedInterSensor.edgeAreas =new int [2][];
			Rectangle [][][] bounds=new Rectangle[2][][];
			for (int side=0;side<2;side++) if (((sides & (1<<side))!=0) && (PIXEL_MAPPING.lastUsedInterSensor.thinEdges[side]!=null)){
				bounds[side]=new Rectangle[1][0];
				PIXEL_MAPPING.lastUsedInterSensor.edgeAreas[side]=PIXEL_MAPPING.lastUsedInterSensor.edgesSegmeniting.extractEdgesAreas(
						gapsNotEdges,
						thresholdStart,
	    				PIXEL_MAPPING.lastUsedInterSensor.distanceFromEdges[side],
	    				PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
	    				bounds[side],  // null or one element Rectangle[] - will return bounding rectangle for each of the edge areas
	    				UPDATE_STATUS,
	    				DEBUG_LEVEL);
				
				
			}
			float [][]fEdgeAreas=new float [PIXEL_MAPPING.lastUsedInterSensor.edgeAreas.length][];
			for (int side=0;side<2;side++) if ((sides & (1<<side))!=0){
				fEdgeAreas[side]=new float[PIXEL_MAPPING.lastUsedInterSensor.edgeAreas[side].length];
				for (int i=0;i<PIXEL_MAPPING.lastUsedInterSensor.edgeAreas[side].length;i++)
					fEdgeAreas[side][i]=PIXEL_MAPPING.lastUsedInterSensor.edgeAreas[side][i];
			} else fEdgeAreas[side]=null;
			String [] titles={"Left","Right"};
			this.SDFA_INSTANCE.showArrays(
					fEdgeAreas,
					PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
					PIXEL_MAPPING.lastUsedInterSensor.mapHeight,
					true,
					"Edge_areas-TH"+((""+thresholdStart).replace('.','_')),
					titles);
//			PIXEL_MAPPING.lastUsedInterSensor.edgeAreas =new int [2][];
			if (gapsNotEdges) return;
			if (thresholdFinal>thresholdStart) return; // don't vacuum
			for (int side=0;side<2;side++) if ((sides & (1<<side))!=0){
				for (int edgeAreaNumber=0; edgeAreaNumber<bounds[side][0].length; edgeAreaNumber++){
					if (DEBUG_LEVEL>0) System.out.println("Processing edge area "+(edgeAreaNumber+1)+" of "+bounds[side][0].length);
					ArrayList<ArrayList<Integer>> boundariesList=PIXEL_MAPPING.lastUsedInterSensor.edgesSegmeniting.createInitialEdgeAreaBorder(
							PIXEL_MAPPING.lastUsedInterSensor.edgeAreas[side],
							edgeAreaNumber+1, // starting with 1
							//						int imageWidth,
							bounds[side][0][edgeAreaNumber], // bounding rectangle to speed up search. Should not have y, x<2 and > appropriate dimensions 
							UPDATE_STATUS,
							DEBUG_LEVEL);
					if (DEBUG_LEVEL>0){
						System.out.println("Created "+boundariesList.size()+ " boundaries for edge area #"+edgeAreaNumber);
						for (int i=0;i<boundariesList.size();i++)
							System.out.println("---Boundary segment "+i+": - length= "+boundariesList.get(i).size());
					}
					int numNewCells=PIXEL_MAPPING.lastUsedInterSensor.edgesSegmeniting.vacuumEdgeAreaBorder(
		    				thresholdStart, // the same as for initial segmenting
		    				thresholdFinal, // threshold (lower) for the shrank edge areas
		    				thresholdStep,  // decrease threshold between passes
		    				repeatEachStep,
		    				repeatLastStep,

		    				PIXEL_MAPPING.lastUsedInterSensor.distanceFromEdges[side], // array of distances from the nearest edge
		    				PIXEL_MAPPING.lastUsedInterSensor.edgeAreas[side],   // or null. If not null, will zero out removed areas
		    				edgeAreaNumber, // any if edgeAreas is null
		    				boundariesList, // will be modified
		    				PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
		    				bounds[side][0][edgeAreaNumber], // bounding rectangle to speed up search. Should not have y, x<2 and > appropriate dimensions 
		    				UPDATE_STATUS,
		    				DEBUG_LEVEL);
					if (DEBUG_LEVEL>0){
						System.out.println("After vacuuming "+boundariesList.size()+ " boundaries for edge area #"+edgeAreaNumber+" New cells: "+numNewCells);
						for (int i=0;i<boundariesList.size();i++)
							System.out.println("---Boundary segment "+i+": - length= "+boundariesList.get(i).size());
					}
					
				}
			}
			float [][]fEdgeAreasVacuumed=new float [PIXEL_MAPPING.lastUsedInterSensor.edgeAreas.length][];
			for (int side=0;side<2;side++) if ((sides & (1<<side))!=0){
				fEdgeAreasVacuumed[side]=new float[PIXEL_MAPPING.lastUsedInterSensor.edgeAreas[side].length];
				if (fEdgeAreasVacuumed[side]==null) System.out.println("fEdgeAreasVacuumed["+side+"]==null");
				if (PIXEL_MAPPING.lastUsedInterSensor.edgeAreas[side]==null) System.out.println("PIXEL_MAPPING.lastUsedInterSensor.edgeAreas["+side+"]==null");
				for (int i=0;i<PIXEL_MAPPING.lastUsedInterSensor.edgeAreas[side].length;i++)
					fEdgeAreasVacuumed[side][i]=PIXEL_MAPPING.lastUsedInterSensor.edgeAreas[side][i];
			} else fEdgeAreasVacuumed[side]=null;
			this.SDFA_INSTANCE.showArrays(
					fEdgeAreasVacuumed,
					PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
					PIXEL_MAPPING.lastUsedInterSensor.mapHeight,
					true,
					"Vacuumed-THS"+((""+thresholdStart).replace('.','_'))+
					"-THF"+((""+thresholdFinal).replace('.','_'))+
					"-STP"+((""+thresholdStep).replace('.','_')),
					titles);
			return;
		}	

		
/* ======================================================================== */
		if       (label.equals("Test Edges")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) {
				String msg="PIXEL_MAPPING is not initialized";
				System.out.println(msg);
				IJ.showMessage(msg);
				return;
			}
			if ((PIXEL_MAPPING.lastUsedInterSensor==null) || (PIXEL_MAPPING.lastUsedInterSensor.booleanEdges==null)) {
				String msg="booleanEdges is not initialized (use \"Show Sobel\"";
				System.out.println(msg);
				IJ.showMessage(msg);
				return;
			}
			GenericDialog gd=new GenericDialog("Edge detection");
//			gd.addNumericField("Start X",226,0,5,"pix");
//			gd.addNumericField("Start Y",436,0,5,"pix");
			gd.addNumericField("Start X",120,0,5,"pix");
			gd.addNumericField("Start Y",52,0,5,"pix");
			gd.addNumericField("Minimal cycle area",0,0,5,"pix");
			gd.addNumericField("Ends test radius",10,0,5,"pix");
			gd.addNumericField("debugXMin",PIXEL_MAPPING.lastUsedInterSensor.debugXMin,0,5,"pix");
			gd.addNumericField("debugXMax",PIXEL_MAPPING.lastUsedInterSensor.debugXMax,0,5,"pix");
			gd.addNumericField("debugYMin",PIXEL_MAPPING.lastUsedInterSensor.debugYMin,0,5,"pix");
			gd.addNumericField("debugYMax",PIXEL_MAPPING.lastUsedInterSensor.debugYMax,0,5,"pix");
			gd.showDialog();
			if (gd.wasCanceled()) return;
			int startX= (int) gd.getNextNumber();
			int startY= (int) gd.getNextNumber();
			int minCycleArea= (int) gd.getNextNumber();
			int endRadius= (int) gd.getNextNumber();
			PIXEL_MAPPING.lastUsedInterSensor.debugXMin=(int) gd.getNextNumber();
			PIXEL_MAPPING.lastUsedInterSensor.debugXMax=(int) gd.getNextNumber();
			PIXEL_MAPPING.lastUsedInterSensor.debugYMin=(int) gd.getNextNumber();
			PIXEL_MAPPING.lastUsedInterSensor.debugYMax=(int) gd.getNextNumber();

			String [] titles={"analog","time", "area","length","parents","trimmed","parents"};
			int side=0;
			float [][]fTest=PIXEL_MAPPING.lastUsedInterSensor.testEdgeThinning(
					PIXEL_MAPPING.lastUsedInterSensor.sobelY[side],
					PIXEL_MAPPING.lastUsedInterSensor.booleanEdges[side],
					PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
					startX, // start looking for the edge to try
					startY,
					minCycleArea, // minimal area of the closed loop to keep
					endRadius,
					//TODO: - other parameters to trim short branches and detect ends?
					true, //UPDATE_STATUS,
					DEBUG_LEVEL);
			this.SDFA_INSTANCE.showArrays(
					fTest,
					PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
					PIXEL_MAPPING.lastUsedInterSensor.mapHeight,
					true,
			"Test-TH"+((""+PIXEL_MAPPING.lastUsedInterSensor.edgeThresholdHigh).replace('.','_'))+"-TL"+((""+PIXEL_MAPPING.lastUsedInterSensor.edgeThresholdLow).replace('.','_')),
			titles);
			return;
//		addButton("Show Sobel",panelStereo);
		}	

/* ======================================================================== */
		
		if       (label.equals("Intercam correlations")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) {
				String msg="PIXEL_MAPPING is not initialized";
				System.out.println(msg);
				IJ.showMessage(msg);
			}
			if ((PIXEL_MAPPING.lastUsedInterSensor==null) || (PIXEL_MAPPING.lastUsedInterSensor.overlapImages==null)) {
				String msg="Inter-sensor data is not initialized (use \"Test Intermaps\"";
				System.out.println(msg);
				IJ.showMessage(msg);
			}
			
			GenericDialog gd=new GenericDialog("Testing of inter-sensor correlation");
			gd.addCheckbox    ("Autocorrelation",                   false);
			gd.addNumericField("Correlation area size (power of 2)",256,0);
			gd.addNumericField("Selection area center horizontal/right (XC)",351,0,4,"pix");
			gd.addNumericField("Selection area center vertical/down  (YC)",978,0,4,"pix");
			gd.addNumericField("Selection area shift between images in a pair",28,0,4,"pix");
			gd.addNumericField("Use phase correlation (0.0 - pure normal correlation, 1.0 - pure phase one)",0.5,3,5,"");// 0.25?
			gd.addNumericField("Weight of Cb (relative to Y) in correlation",0.5,3,5,"");
			gd.addNumericField("Weight of Cr (relative to Y) in correlation",0.5,3,5,"");
			gd.addNumericField("correlationHighPassSigma, - pixels in frequency domain",1.5,3,5,"");
			gd.addNumericField("correlationLowPassSigma, - fraction of the frequency range",0.6,3,5,"");
			gd.addNumericField("Noise normalization sigma for Y-component",3.0,2,5,"pix");
			gd.addNumericField("Noise normalization sigma for Cb, Cr components",5.0,2,5,"pix");
			gd.addNumericField("Contrast threshold",1.5,1,5,"");
			gd.addCheckbox    ("Use binary alpha)", true);
			gd.addCheckbox    ("Enable negative disparity)", false);
			gd.showDialog();
			if (gd.wasCanceled()) return;
			boolean autocorrelation=gd.getNextBoolean();
			int corrFFTSize= (int) gd.getNextNumber();
			int corrXC=      (int) gd.getNextNumber();
			int corrYC=      (int) gd.getNextNumber();
			int corrShift=   (int) gd.getNextNumber();
			double corrPhaseFraction=gd.getNextNumber();
			double corrCbWeight=gd.getNextNumber();
			double corrCrWeight=gd.getNextNumber();
			double correlationHighPassSigma=gd.getNextNumber();
			double correlationLowPassSigma=gd.getNextNumber();
			double noiseNormalizationSignaY=   gd.getNextNumber();
			double noiseNormalizationSignaCbCr=gd.getNextNumber();
			double contrastThreshold=gd.getNextNumber();
			boolean useBinaryAlpha=gd.getNextBoolean();
			boolean enableNegativeDisparity=gd.getNextBoolean(); 
			double [][] corr=PIXEL_MAPPING.lastUsedInterSensor.correlate(
					autocorrelation,
					corrFFTSize,
					corrXC,
					corrYC,
					corrShift,
					corrPhaseFraction,
					corrCbWeight,
					corrCrWeight,
					correlationHighPassSigma,
					correlationLowPassSigma,
					noiseNormalizationSignaY,
					noiseNormalizationSignaCbCr,
					contrastThreshold,
					useBinaryAlpha,
					enableNegativeDisparity,
					THREADS_MAX,
					DEBUG_LEVEL);
			String [] titles={"Combo","Y","Cb","Cr"};
			this.SDFA_INSTANCE.showArrays(
					corr,
					corrFFTSize,
					corrFFTSize,
					true,
			"corr-x"+corrXC+"-y"+corrYC+"-SHFT"+corrShift+"-PC"+((""+corrPhaseFraction).replace('.','_')),
			titles);
			return;
		}

/* ======================================================================== */
		if       (label.equals("Test Line")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) {
				String msg="PIXEL_MAPPING is not initialized";
				System.out.println(msg);
				IJ.showMessage(msg);
			}
			if ((PIXEL_MAPPING.lastUsedInterSensor==null) || (PIXEL_MAPPING.lastUsedInterSensor.overlapImages==null)) {
				String msg="Inter-sensor data is not initialized (use \"Test Intermaps\"";
				System.out.println(msg);
				IJ.showMessage(msg);
			}
			GenericDialog gd=new GenericDialog("Testing line detection");
			gd.addNumericField("Correlation area size (power of 2)",32,0);
			gd.addNumericField("Selection area center horizontal/right (XC)",376,0,4,"pix");
			gd.addNumericField("Selection area center vertical/down  (YC)",826,0,4,"pix");
			gd.addNumericField("Use phase correlation (0.0 - pure normal correlation, 1.0 - pure phase one)",0.5,3,5,"");// 0.25?
			gd.addNumericField("correlationHighPassSigma, - pixels in frequency domain",1.5,3,5,"");
			gd.addNumericField("correlationLowPassSigma, - fraction of the frequency range",0.6,3,5,"");
			gd.addNumericField("Threshold (relative to RMS) ",5.0,3,5,"");
			
			gd.addCheckbox    ("Remove non-connected \"islands\"",true);
			gd.addCheckbox    ("Adjust direction",false);
			//gd.addCheckbox    ("Add PI to phase",true);
			gd.addNumericField("Phase integration width ",3.5,3,5,"frequency samples");
			gd.addNumericField("Result high-pass filter",1.0,3,5,"frequency samples");
			gd.addNumericField("Phase dispersion cost when finding direction (0.0 - only amplitudes) ",0.5,3,5,"");
			gd.addNumericField("Zero area half size (for phase continuity) ",0.5,3,5,"pix");
			gd.addNumericField("Feature filter 0 - any, +1 white lines, +2 - black lines, + 4 black-to-white in the direction, 8 - white-to-black", 0,0);
			gd.addNumericField("Diff phase - RMS fraction ",2.0,3,5,"");
			gd.addNumericField("Diff phase - absolute threshold ",0.0,3,5,"");
			gd.addNumericField("Diff phase - maximal phase mismatch",45.0,3,5,"degrees");
			gd.addNumericField("Simulation phase mismatch tolerance",45.0,3,5,"degrees");

			gd.showDialog();
			if (gd.wasCanceled()) return;
			int corrFFTSize= (int) gd.getNextNumber();
			int corrXC=      (int) gd.getNextNumber();
			int corrYC=      (int) gd.getNextNumber();
			double corrPhaseFraction=gd.getNextNumber();
			double correlationHighPassSigma=gd.getNextNumber();
			double correlationLowPassSigma= gd.getNextNumber();
			double threshold=               gd.getNextNumber();
			boolean removeIslands=          gd.getNextBoolean();
			boolean adjustDirection=        gd.getNextBoolean();
			double phaseIntegrationWidth=   gd.getNextNumber();
			double resultHighPass=          gd.getNextNumber();
			double dispertionCost=          gd.getNextNumber();
//			boolean addPI=                  gd.getNextBoolean();
			double zeroBinHalfSize=         gd.getNextNumber();
			int featureFilter=        (int) gd.getNextNumber();
			double minRMSFrac=              gd.getNextNumber();
			double minAbs=                  gd.getNextNumber();
			double maxPhaseMismatch=        Math.PI/180.0*gd.getNextNumber();
			double phaseTolerance=          Math.PI/180.0*gd.getNextNumber();

			PIXEL_MAPPING.lastUsedInterSensor.testLineDetect(
					corrFFTSize,
					corrXC,
					corrYC,
					corrPhaseFraction,
					correlationHighPassSigma,
					correlationLowPassSigma,
					removeIslands,
					threshold,
					adjustDirection,
					//addPI,
					phaseIntegrationWidth,
					resultHighPass,
					dispertionCost,
					featureFilter,
					zeroBinHalfSize,
					minRMSFrac,
					minAbs,
					maxPhaseMismatch,
					phaseTolerance,
					DEBUG_LEVEL);
			return;
		}
/* ======================================================================== */
		if       (label.equals("Linear Features")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) {
				String msg="PIXEL_MAPPING is not initialized";
				System.out.println(msg);
				IJ.showMessage(msg);
			}
			if ((PIXEL_MAPPING.lastUsedInterSensor==null) || (PIXEL_MAPPING.lastUsedInterSensor.overlapImages==null)) {
				String msg="Inter-sensor data is not initialized (use \"Test Intermaps\"";
				System.out.println(msg);
				IJ.showMessage(msg);
			}
			GenericDialog gd=new GenericDialog("Testing line detection");
			gd.addCheckbox    ("(Re)-extract linear features", true);
			gd.addCheckbox    ("Ignore feature phase", false);
			gd.addCheckbox    ("Preserve feature DC", false);
			gd.addNumericField("Scale features mode: -2.0 - no scale, -1.0 - phase strength, 0.0 - absolute strength, 1.0 - relative strength",0.5,3,5,"");
			
			String [] choices={"left image","right image","both images"};
			gd.addChoice("Extract linear features for",choices, choices[2]);
			gd.addNumericField("Tile size (power of 2)",32,0);
			gd.addNumericField("Overlap as fraction of the tile width (FFT size)",8,0);
			gd.addNumericField("Minimal fractional tile to use",50,1,5,"%");
			gd.addCheckbox    ("Use binary alpha (all >0.0 treat as 1.0)", false);
			
			gd.addNumericField("correlationHighPassSigma, - pixels in frequency domain",1.5,3,5,"");
			gd.addNumericField("correlationLowPassSigma, - fraction of the frequency range",0.6,3,5,"");
			gd.addNumericField("Phase integration width ",3.5,3,5,"frequency samples");
			gd.addNumericField("Result high-pass filter",1.0,3,5,"frequency samples");
			
//			gd.addNumericField("Zero area half size (for phase continuity) ",0.5,3,5,"pix");
			gd.addNumericField("Feature filter 0 - any, +1 white lines, +2 - black lines, + 4 black-to-white in the direction, 8 - white-to-black", 0,0);
			gd.addNumericField("Minimal frequency sample value relative to RMS to be considered when looking for the linear phase feature (0.0 - skip test) ",2.0,3,5,"");
			gd.addNumericField("Minimal frequency sample absolute value to be considered when looking for the linear phase feature  (0.0 - skip test)",0.0,3,5,"");
			gd.addNumericField("Maximal phase mismatch (between diagonals in a 2x2 sqaure) to use frequaency sample",45.0,3,5,"degrees");
			gd.addNumericField("Phase dispersion cost when finding direction (0.0 - only amplitudes) ",0.5,3,5,"");
//			gd.addNumericField("Simulation phase mismatch tolerance",45.0,3,5,"degrees");
			gd.addCheckbox    ("Calculate features strengths", true);
			
			gd.addMessage     ("Merging features:");
			gd.addNumericField("Direction tolerance",20.0,3,5,"degrees");
			gd.addNumericField("Normal distance  tolerance",2.0,3,5,"pixels");
			gd.addNumericField("Tangential distance  tolerance",4.0,3,5,"pixels");
			gd.addNumericField("Hosts cells tolerance (if two neighbors claim same point as it's own",1.0,3,5,"pixels");
			gd.addNumericField("Direction sigma - fraction of direction tolerance",0.6,3,5,"x");
			gd.addNumericField("Normal distance  sigma - fraction of normal distance  tolerance",0.6,3,5,"x");
			gd.addNumericField("Tangential distance  sigma - fraction of normal distance  tolerance",0.6,3,5,"x");
			gd.addNumericField("Minimal number of merged (0 - including single with off-cell features, 1 - including same cell no-merge, >1 - merged)",2,0);
			gd.addNumericField("Scale distance (compensate for apparent decrease of measured distance to feature caused by windowing)",1.0,3,5,"x");
			
			gd.addNumericField("Maximal tangential shift while moving features to the closer cells",4.0,3,5,"x");
			gd.addNumericField("Number of cells around feature center to search for the new host cell",3,0,5," cells each direction");
			
			gd.addCheckbox    ("Multiply by cell usage", true);
			gd.addNumericField("Cell usage shift",1.0,3,5,"x");

			gd.addCheckbox    ("Show cell usage", false);
			gd.addNumericField("Cell usage display scale (if enabled))",0.1,3,5,"x");
			
			gd.addNumericField("Debug row",-10,0);
			gd.addNumericField("Debug column",-10,0);
			WindowTools.addScrollBars(gd);
			gd.showDialog();
			if (gd.wasCanceled()) return;
			boolean extractFeatures=        gd.getNextBoolean();
			boolean ignorePhase=            gd.getNextBoolean();
			boolean preserveDC=              gd.getNextBoolean();
			double strengthMode=             gd.getNextNumber();
			int sides=                      gd.getNextChoiceIndex()+1; // 1-2-3
			int corrFFTSize=          (int) gd.getNextNumber();
			int overlapFraction=      (int) gd.getNextNumber();
			double alphaThreshold=     0.01*gd.getNextNumber();
			boolean useBinaryAlpha=         gd.getNextBoolean();
			double correlationHighPassSigma=gd.getNextNumber();
			double correlationLowPassSigma= gd.getNextNumber();
			double phaseIntegrationWidth=   gd.getNextNumber();
			double resultHighPass=          gd.getNextNumber();

//			double zeroBinHalfSize=         gd.getNextNumber();
			int featureFilter=        (int) gd.getNextNumber();
			double minRMSFrac=              gd.getNextNumber();
			double minAbs=                  gd.getNextNumber();
			double maxPhaseMismatch=        Math.PI/180.0*gd.getNextNumber();
			double dispertionCost=          gd.getNextNumber();
			
			boolean calculateStrength=      gd.getNextBoolean();
			
			double directionTolerance=      Math.PI/180.0*gd.getNextNumber();
			double normalDistanceTolerance= gd.getNextNumber();
			double tangentialDistanceTolerance= gd.getNextNumber();
			double hostsTolerance=          gd.getNextNumber();
			double directionFracSigma=      gd.getNextNumber(); // make a fixed fraction of directionTolerance?
			double normalDistanceFracSigma= gd.getNextNumber(); // make a fixed fraction of distanceTolerance?
			double tangentialDistanceFracSigma= gd.getNextNumber(); // make a fixed fraction of distanceTolerance?
			int minMerged=            (int) gd.getNextNumber();
			double scaleDistances=          gd.getNextNumber();
			
			double swapTangentialTolerance= gd.getNextNumber();
			int    swapSearchRange=   (int) gd.getNextNumber();
			
			boolean multiplyByCellUsage=    gd.getNextBoolean();
			double cellUsageShift=          gd.getNextNumber();
			if (!multiplyByCellUsage) cellUsageShift=Double.NaN;

			boolean enableDisplayUsage=     gd.getNextBoolean();
			double displayUsageScale=       gd.getNextNumber()*(enableDisplayUsage?1.0:0.0);
    		int debugRow=             (int) gd.getNextNumber();
    		int debugColumn=          (int) gd.getNextNumber();


			double directionSigma=     directionTolerance*directionFracSigma;  // make a fixed fraction of directionTolerance?
			double normalDistanceSigma=      normalDistanceTolerance*normalDistanceFracSigma;// make a fixed fraction of distanceTolerance?
			double tangentialDistanceSigma=  tangentialDistanceTolerance*tangentialDistanceFracSigma;// make a fixed fraction of distanceTolerance?
			
			
//			double phaseTolerance=          Math.PI/180.0*gd.getNextNumber();
			
			if (PIXEL_MAPPING.lastUsedInterSensor.linearFeatures==null) {
				PIXEL_MAPPING.lastUsedInterSensor.linearFeatures=new double [2][][][];
				for (int i=0;i<PIXEL_MAPPING.lastUsedInterSensor.linearFeatures.length;i++) PIXEL_MAPPING.lastUsedInterSensor.linearFeatures[i]=null;
			}
			long startTime = System.nanoTime();

			for (int side=0;side<2;side++) if ((sides & (1<<side))!=0) {
				boolean reExtract=extractFeatures || (PIXEL_MAPPING.lastUsedInterSensor.linearFeatures[side]==null);
				if (reExtract) PIXEL_MAPPING.lastUsedInterSensor.linearFeatures[side]= PIXEL_MAPPING.lastUsedInterSensor.detectLinearFeatures(
		    			(side>0),
		    			corrFFTSize,
		    			overlapFraction, // default 8
		    			alphaThreshold, // minFracArea
		    			useBinaryAlpha,
						correlationHighPassSigma,
						correlationLowPassSigma,
						phaseIntegrationWidth,
						resultHighPass, //cutting "ribbon" near zero frequency - influences length of the detected line (used only for strength calculation)
						dispertionCost,
						featureFilter,
						minRMSFrac,
						minAbs,
						maxPhaseMismatch,
						calculateStrength,
		        		debugRow,
		        		debugColumn,
						THREADS_MAX,
						DEBUG_LEVEL);
				if (DEBUG_LEVEL>1) System.out.println("Linear Features for "+((side>0)?"right":"left")+" side finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));


				// startTime	
			}
			// startTime
			double [][] null2 = {null,null};
			double [][] null4 = {null,null,null,null};
			
			double [][] linearFeatures =(displayUsageScale>0)?null4:null2;
			for (int side=0;side<2;side++) if ((sides & (1<<side))!=0) {
				if (PIXEL_MAPPING.lastUsedInterSensor.linearFeatures[side]!=null){
					linearFeatures[side<<((displayUsageScale>0)?1:0)]=PIXEL_MAPPING.lastUsedInterSensor.reconstructImageFeature(
							corrFFTSize,
							overlapFraction, // default 8
							PIXEL_MAPPING.lastUsedInterSensor.linearFeatures[side],
							ignorePhase,
							preserveDC,
							strengthMode, // 0.0 - absolute, 1.0 - relative, 0.5 - "balanced", "-1" - use phaseStrength instead 
							phaseIntegrationWidth, // use the same as during extraction?
							resultHighPass, //cutting "ribbon" near zero frequency - influences length of the detected line
			        		debugRow,
			        		debugColumn,
							THREADS_MAX,
							DEBUG_LEVEL);
					if (displayUsageScale>0){
						linearFeatures[(side<<1)+1]=PIXEL_MAPPING.lastUsedInterSensor.displayUsage(
								PIXEL_MAPPING.lastUsedInterSensor.linearFeatures[side],
				    			corrFFTSize,
				    			overlapFraction, // default 8
				    			displayUsageScale);
					}
					if (DEBUG_LEVEL>1) System.out.println("Linear Features simulated for "+((side>0)?"right":"left")+" side finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
				}
			}			
//TODO: just show linearFeatures array and measure time			
			String [] titlesLR={"Left","Right"};
			String [] titlesLRU={"Left","Usage Left","Right","Usage Right"};
			String [] titles=(displayUsageScale>0)?titlesLRU:titlesLR;
			//			double displayUsageScale=       gd.getNextNumber();

			this.SDFA_INSTANCE.showArrays(
					linearFeatures,
					PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
					PIXEL_MAPPING.lastUsedInterSensor.mapHeight,
					true,
					"Linear features",
					titles);
			double [][][][] filteredFeatures={null,null};
			double [][] filteredFeaturesImage =(displayUsageScale>0)?null4:null2;;
			for (int side=0;side<2;side++) if ((sides & (1<<side))!=0) {
				
				filteredFeatures[side]=PIXEL_MAPPING.lastUsedInterSensor.mergeLinearFeatures(
		    			corrFFTSize,
		    			overlapFraction, // default 8
		    			PIXEL_MAPPING.lastUsedInterSensor.linearFeatures[side],
		    			strengthMode, // 0.0 - absolute, 1.0 - relative, 0.5 - "balanced", "-1" - use phaseStrength instead 
		    			directionTolerance,
		    			normalDistanceTolerance,
		    			tangentialDistanceTolerance,
		     			hostsTolerance,
		    			directionSigma, // make a fixed fraction of directionTolerance?
		    			normalDistanceSigma, // make a fixed fraction of distanceTolerance?
		    			tangentialDistanceSigma, // make a fixed fraction of distanceTolerance?
		    			minMerged,
		    			cellUsageShift,
		    			scaleDistances,
		    			swapTangentialTolerance,
		    			swapSearchRange,
		        		debugRow,
		        		debugColumn,
		    			DEBUG_LEVEL);
				if (PIXEL_MAPPING.lastUsedInterSensor.linearFeatures[side]!=null){
					filteredFeaturesImage[side<<((displayUsageScale>0)?1:0)]=PIXEL_MAPPING.lastUsedInterSensor.reconstructImageFeature(
							corrFFTSize,
							overlapFraction, // default 8
							filteredFeatures[side],
							ignorePhase,
							preserveDC,
							strengthMode, // 0.0 - absolute, 1.0 - relative, 0.5 - "balanced", "-1" - use phaseStrength instead 
							phaseIntegrationWidth, // use the same as during extraction?
							resultHighPass, //cutting "ribbon" near zero frequency - influences length of the detected line
			        		debugRow,
			        		debugColumn,
							THREADS_MAX,
							DEBUG_LEVEL);
					if (displayUsageScale>0){
						filteredFeaturesImage[(side<<1)+1]=PIXEL_MAPPING.lastUsedInterSensor.displayUsage(
								filteredFeatures[side],
				    			corrFFTSize,
				    			overlapFraction, // default 8
				    			displayUsageScale);
					}

					if (DEBUG_LEVEL>1) System.out.println("Filtered linear Features simulated for "+((side>0)?"right":"left")+" side finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
				}
			}			
//TODO: just show linearFeatures array and measure time			
			this.SDFA_INSTANCE.showArrays(
					filteredFeaturesImage,
					PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
					PIXEL_MAPPING.lastUsedInterSensor.mapHeight,
					true,
					"Filtered Linear features",
					titles);
			
			return;
		}
		
/* ======================================================================== */
		
		if       (label.equals("Intercam rectangular")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) {
				String msg="PIXEL_MAPPING is not initialized";
				System.out.println(msg);
				IJ.showMessage(msg);
			}
			if ((PIXEL_MAPPING.lastUsedInterSensor==null) || (PIXEL_MAPPING.lastUsedInterSensor.overlapImages==null)) {
				String msg="Inter-sensor data is not initialized (use \"Test Intermaps\"";
				System.out.println(msg);
				IJ.showMessage(msg);
			}
			
			GenericDialog gd=new GenericDialog("Testing of inter-sensor correlation");
			gd.addCheckbox    ("Autocorrelation",                   false);
			gd.addNumericField("Correlation area size (power of 2)",32,0);
			gd.addNumericField("Overlap as fraction of the tile width (FFT size)",8,0);
			gd.addNumericField("Selection area center horizontal/right (XC)",326,0,4,"pix");
			gd.addNumericField("Selection area center vertical/down  (YC)",488,0,4,"pix");
			gd.addNumericField("Maximal disparity",50,0,4,"pix");
			gd.addNumericField("Use phase correlation (0.0 - pure normal correlation, 1.0 - pure phase one)",0.5,3,5,"");// 0.25?
			gd.addNumericField("Weight of Cb (relative to Y) in correlation",0.5,3,5,"");
			gd.addNumericField("Weight of Cr (relative to Y) in correlation",0.5,3,5,"");
			gd.addNumericField("correlationHighPassSigma, - pixels in frequency domain",1.5,3,5,"");
			gd.addNumericField("correlationLowPassSigma, - fraction of the frequency range",0.6,3,5,"");
			gd.addNumericField("Noise normalization sigma for Y-component",3.0,2,5,"pix");
			gd.addNumericField("Noise normalization sigma for Cb, Cr components",5.0,2,5,"pix");
			gd.addNumericField("Contrast threshold",1.5,1,5,"");
			gd.addCheckbox    ("Use binary alpha", true);
			gd.addCheckbox    ("Enable negative disparity)", false);
			gd.addMessage     ("---- For testing correlateOneRow(): ---");
			gd.addCheckbox    ("Calculate for right (second) image", false);
			gd.addNumericField("Minimal area fraction",50,1,5,"%");
			gd.addMessage     ("---- For testing correlateImage(): ---");
			gd.addCheckbox    ("Correlate full image", true);
			//correlateImage
			
			gd.showDialog();
			if (gd.wasCanceled()) return;
			boolean autocorrelation=gd.getNextBoolean();
			int corrFFTSize=     (int) gd.getNextNumber();
			int overlapFraction= (int) gd.getNextNumber();
			int corrXC=          (int) gd.getNextNumber();
			int corrYC=          (int) gd.getNextNumber();
			int maxDisparity=    (int) gd.getNextNumber();
			double corrPhaseFraction=gd.getNextNumber();
			double corrCbWeight=gd.getNextNumber();
			double corrCrWeight=gd.getNextNumber();
			double correlationHighPassSigma=gd.getNextNumber();
			double correlationLowPassSigma=gd.getNextNumber();
			double noiseNormalizationSignaY=   gd.getNextNumber();
			double noiseNormalizationSignaCbCr=gd.getNextNumber();
			double contrastThreshold=gd.getNextNumber();
			boolean useBinaryAlpha=gd.getNextBoolean();
			boolean enableNegativeDisparity=gd.getNextBoolean(); 
			
			boolean secondSide=gd.getNextBoolean();
			double  minFracArea=0.01*gd.getNextNumber();
			
			boolean correlateFullImage=gd.getNextBoolean();
			double [] corr=PIXEL_MAPPING.lastUsedInterSensor.correlateRectangular(
					autocorrelation,
					corrFFTSize,
					overlapFraction,
					corrXC,
					corrYC,
					maxDisparity,
					corrPhaseFraction,
					corrCbWeight,
					corrCrWeight,
					correlationHighPassSigma,
					correlationLowPassSigma,
					noiseNormalizationSignaY,
					noiseNormalizationSignaCbCr,
					contrastThreshold,
					useBinaryAlpha,
					enableNegativeDisparity,
					THREADS_MAX,
					DEBUG_LEVEL);
//			String [] titles={"Combo","Y","Cb","Cr"};
			this.SDFA_INSTANCE.showArrays(
					corr,
					corr.length/corrFFTSize,
					corrFFTSize,
			"corr-x"+corrXC+"-y"+corrYC+"-MAXDSP"+maxDisparity+"-PC"+((""+corrPhaseFraction).replace('.','_')));
			
			double [][] rowCorr=PIXEL_MAPPING.lastUsedInterSensor.correlateOneRow(
					null,             // DoubleFHT doubleFHT
					secondSide,
					autocorrelation,
					corrFFTSize,
					overlapFraction,
					corrYC,
					maxDisparity,
					corrPhaseFraction,
					corrCbWeight,
					corrCrWeight,
					correlationHighPassSigma,
					correlationLowPassSigma,
					noiseNormalizationSignaY,
					noiseNormalizationSignaCbCr,
					minFracArea,
					useBinaryAlpha,
					DEBUG_LEVEL);
			int rowCorrLen=0;
			for (int i=0;i<rowCorr.length;i++) if (rowCorr[i]!=null) {rowCorrLen=rowCorr[i].length;break;}
			
			double [] rowCorrImage=new double [rowCorr.length*rowCorrLen];
			for (int i=0;i<rowCorrImage.length;i++) rowCorrImage[i]=Double.NaN; 
			for (int i=0;i<rowCorr.length;i++) if (rowCorr[i]!=null) for (int j=0;j<rowCorr[i].length;j++) rowCorrImage[j*rowCorr.length+i]=rowCorr[i][j];
			this.SDFA_INSTANCE.showArrays(
					rowCorrImage,
					rowCorr.length,
					rowCorrLen,
			"corrRow"+corrYC+(secondSide?"R":"L")+(autocorrelation?"A":"")+"-MAXDSP"+maxDisparity+"-PC"+((""+corrPhaseFraction).replace('.','_')));
			if (correlateFullImage) {
				if (PIXEL_MAPPING.lastUsedInterSensor.disparityMap==null) PIXEL_MAPPING.lastUsedInterSensor.disparityMap=new double [2][2][][][];
				PIXEL_MAPPING.lastUsedInterSensor.disparityMap[secondSide?1:0][autocorrelation?1:0]=PIXEL_MAPPING.lastUsedInterSensor.correlateImage(
						secondSide,
						autocorrelation,
						corrFFTSize,
						overlapFraction,
						maxDisparity,
						corrPhaseFraction,
						corrCbWeight,
						corrCrWeight,
						correlationHighPassSigma,
						correlationLowPassSigma,
						noiseNormalizationSignaY,
						noiseNormalizationSignaCbCr,
						minFracArea,
						useBinaryAlpha,
						THREADS_MAX,
						DEBUG_LEVEL);
				PIXEL_MAPPING.lastUsedInterSensor.disparityMapFFTSize=corrFFTSize;
				PIXEL_MAPPING.lastUsedInterSensor.disparityMapOverlapFraction=overlapFraction;
			}
//PIXEL_MAPPING.lastUsedInterSensor
			return;
		}
/* ======================================================================== */
		if       (label.equals("Thershold Disparity")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) {
				String msg="PIXEL_MAPPING is not initialized";
				System.out.println(msg);
				IJ.showMessage(msg);
			}
			if ((PIXEL_MAPPING.lastUsedInterSensor==null) || (PIXEL_MAPPING.lastUsedInterSensor.overlapImages==null)) {
				String msg="Inter-sensor data is not initialized (use \"Test Intermaps\"";
				System.out.println(msg);
				IJ.showMessage(msg);
			}
			if (PIXEL_MAPPING.lastUsedInterSensor.disparityMap==null) {
				String msg="disparityMap[][][] not initialized (use \"Intercam rectangular\"";
				System.out.println(msg);
				IJ.showMessage(msg);
			}
			
			GenericDialog gd=new GenericDialog("Testing of inter-sensor correlation");
			gd.addCheckbox    ("Calculate for right (second) image", false);
			gd.addCheckbox    ("Autocorrelation",                   false);
			gd.addNumericField("Contrast threshold",1.5,1,5,"");
			//correlateImage
			gd.showDialog();
			if (gd.wasCanceled()) return;
			boolean secondSide=gd.getNextBoolean();
			boolean autocorrelation=gd.getNextBoolean();
			double contrastThreshold=gd.getNextNumber();
			double [] thresholdedDisparity=PIXEL_MAPPING.lastUsedInterSensor.showThresholdedDisparityPixels(
					PIXEL_MAPPING.lastUsedInterSensor.disparityMap[secondSide?1:0][autocorrelation?1:0],
					PIXEL_MAPPING.lastUsedInterSensor.disparityMapFFTSize/2, // 1/2 corrFFTSize
					contrastThreshold);
			this.SDFA_INSTANCE.showArrays(
					thresholdedDisparity,
//					PIXEL_MAPPING.lastUsedInterSensor.disparityMap[secondSide?1:0][autocorrelation?1:0][0].length,
//					PIXEL_MAPPING.lastUsedInterSensor.disparityMap[secondSide?1:0][autocorrelation?1:0].length,
					PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
					PIXEL_MAPPING.lastUsedInterSensor.mapHeight,
			"D"+((""+contrastThreshold).replace('.','_'))+(secondSide?"R":"L")+(autocorrelation?"A":""));
			return;
		}
/* ======================================================================== */
		
		if       (label.equals("Disparity Tiles0")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) {
				String msg="PIXEL_MAPPING is not initialized";
				System.out.println(msg);
				IJ.showMessage(msg);
			}
			if ((PIXEL_MAPPING.lastUsedInterSensor==null) || (PIXEL_MAPPING.lastUsedInterSensor.overlapImages==null)) {
				String msg="Inter-sensor data is not initialized (use \"Test Intermaps\"";
				System.out.println(msg);
				IJ.showMessage(msg);
			}
			
			GenericDialog gd=new GenericDialog("Disparity tiles");
			gd.addNumericField("Top tile size (half correlation FFT size , power of 2)",128,0);
			gd.addNumericField("Number of levels (each has twice smaller tile size",4,0);
			gd.addNumericField("Use phase correlation (0.0 - pure normal correlation, 1.0 - pure phase one)",0.5,3,5,"");// 0.25?
			gd.addNumericField("Weight of Cb (relative to Y) in correlation",0.5,3,5,"");
			gd.addNumericField("Weight of Cr (relative to Y) in correlation",0.5,3,5,"");
			gd.addNumericField("correlationHighPassSigma, - pixels in frequency domain",1.5,3,5,"");
			gd.addNumericField("correlationLowPassSigma, - fraction of the frequency range",0.6,3,5,"");
			gd.addNumericField("Noise normalization sigma for Y-component",3.0,2,5,"pix");
			gd.addNumericField("Noise normalization sigma for Cb, Cr components",5.0,2,5,"pix");
			gd.addNumericField("Contrast threshold",1.5,1,5,"");
			gd.addNumericField("Decrease threshold for smaller tiles as 1/2^(k*level), k=",0.0,1,5,"");
			gd.addCheckbox    ("Use binary alpha)", true);
			gd.addNumericField("Minimal overlap area",30.0,1,6,"%");
			gd.addCheckbox    ("Enable negative disparity)", false);
			gd.showDialog();
			if (gd.wasCanceled()) return;
			int tileSize=   (int) gd.getNextNumber();
			int numLevels=  (int) gd.getNextNumber();
			double corrPhaseFraction=gd.getNextNumber();
			double corrCbWeight=gd.getNextNumber();
			double corrCrWeight=gd.getNextNumber();
			double correlationHighPassSigma=gd.getNextNumber();
			double correlationLowPassSigma=gd.getNextNumber();
			double noiseNormalizationSignaY=   gd.getNextNumber();
			double noiseNormalizationSignaCbCr=gd.getNextNumber();
			double contrastThreshold=gd.getNextNumber();
			double contrastThresholdDecrease=gd.getNextNumber();
			boolean useBinaryAlpha=gd.getNextBoolean();
			double minTileFraction=0.01*gd.getNextNumber();
			boolean enableNegativeDisparity=gd.getNextBoolean();
			// setup tiles arrays
			int [] tilesStats=PIXEL_MAPPING.lastUsedInterSensor.setupTiles0(
					tileSize,
					numLevels); // number of tile levels, each having half linear size (4..5)
			if (DEBUG_LEVEL>0) System.out.println("Using "+numLevels+" levels of tiles, top level is "+tilesStats[0]+"x"+tilesStats[1]+", largest is "+tileSize+
					" pixels square, smallest - "+(tileSize>>(numLevels-1))+" pixels");
			// generate disparity tiles data
			int [] tilesGenStats= PIXEL_MAPPING.lastUsedInterSensor.generateTileLevels0(
					corrPhaseFraction, // double phaseCorrelationFraction,
					corrCbWeight,
					corrCrWeight,
					correlationHighPassSigma, //double correlationHighPassSigma,
					correlationLowPassSigma, //double correlationLowPassSigma,
					noiseNormalizationSignaY, //double noiseNormalizationSignaY,
					noiseNormalizationSignaCbCr, //double noiseNormalizationSignaCbCr,
					contrastThreshold, //double contrastThreshold,
					contrastThresholdDecrease, //double contrastThresholdDecrease, //Decrease threshold for smaller tiles as 1/2^(k*level), k=
					useBinaryAlpha, //boolean binaryAlpha,
					enableNegativeDisparity, //boolean enableNegativeDisparity,
					minTileFraction, // double minTileFraction,
					null, //DoubleFHT doubleFHT,
					THREADS_MAX, //int threadsMax,
					DEBUG_LEVEL); //int debugLevel)			
			if (DEBUG_LEVEL>0) System.out.println("Generated "+tilesGenStats[0]+" non-empty tiles, "+tilesGenStats[1]+" of them have foreground");
			double [][] disparity=PIXEL_MAPPING.lastUsedInterSensor.collectDisparityFromTiles(DEBUG_LEVEL);
			String [] titles={"Best-shift","Best-contrast","BG-shift","BG-contrast","FG-shift","FG-contrast"};
			this.SDFA_INSTANCE.showArrays(
					disparity,
					PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
					PIXEL_MAPPING.lastUsedInterSensor.mapHeight,
					true,
			"disparity-PC"+((""+corrPhaseFraction).replace('.','_'))+"-THR"+((""+contrastThreshold).replace('.','_')),
			titles);
			return;
		}
/* ======================================================================== */
		
		if       (label.equals("Disparity Tiles")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) {
				String msg="PIXEL_MAPPING is not initialized";
				System.out.println(msg);
				IJ.showMessage(msg);
				return;
			}
			if ((PIXEL_MAPPING.lastUsedInterSensor==null) || (PIXEL_MAPPING.lastUsedInterSensor.overlapImages==null)) {
				String msg="Inter-sensor data is not initialized (use \"Test Intermaps\"";
				System.out.println(msg);
				IJ.showMessage(msg);
				return;
			}
			
			GenericDialog gd=new GenericDialog("Disparity tiles");
			//gd.addString("Generate disparity for ")
			String [] choices={"left image","right image","both images"};
			String [] choices1={"cross-correlate","auto-correlate","cross and auto"};
			gd.addChoice("Generate disparity for",choices, choices[2]);
			gd.addChoice("Cross/auto-correlate",  choices1,choices1[2]);
			gd.addNumericField("Top tile size (half correlation FFT size , power of 2)",128,0);
			gd.addNumericField("Number of levels (each has twice smaller tile size",4,0);
			gd.addNumericField("Use phase correlation (0.0 - pure normal correlation, 1.0 - pure phase one)",0.5,3,5,"");// 0.25?
			gd.addNumericField("Weight of Cb (relative to Y) in correlation",0.5,3,5,"");
			gd.addNumericField("Weight of Cr (relative to Y) in correlation",0.5,3,5,"");
			gd.addNumericField("correlationHighPassSigma, - pixels in frequency domain",1.5,3,5,"");
			gd.addNumericField("correlationLowPassSigma, - fraction of the frequency range",0.6,3,5,"");
			gd.addNumericField("Noise normalization sigma for Y-component",3.0,2,5,"pix");
			gd.addNumericField("Noise normalization sigma for Cb, Cr components",5.0,2,5,"pix");
			gd.addNumericField("Contrast threshold",0.6,1,5,"");
			gd.addNumericField("Decrease threshold for smaller tiles as 1/2^(k*level), k=",0.0,1,5,"");
			gd.addCheckbox    ("Use binary alpha)", true);

			gd.addNumericField("Minimal overlap area",30.0,1,6,"%");
			gd.addNumericField("Minimal disparity to accumulate (-1)", -1,0);
			gd.addNumericField("Maximal disparity to accumulate (tile size/2)", 64,0);
			gd.addNumericField("Relative (to tile size) correlation shift step", 0.05,3,4,"");
			gd.addCheckbox    ("Probe at limits (calculate correlation with shifts close to smallest/largest disparity above thershold)", true);
			gd.addMessage("Filtering parameters");
			gd.addNumericField("Fraction of tile parent to combine ",0.0,3,5,"");
			gd.addNumericField("Fraction of tile horizontal neighbors to combine ",0.0,3,5,"");
			gd.addNumericField("Neighbos power coefficient higher the more close to max(neighbors) ",4.0,3,5,"");
			gd.showDialog();
			if (gd.wasCanceled()) return;
			int sides=               gd.getNextChoiceIndex()+1; // 1-2-3
			int selves=              gd.getNextChoiceIndex()+1; // 1-2-3
			
			int tileSize=   (int) gd.getNextNumber();
			int numLevels=  (int) gd.getNextNumber();
			double corrPhaseFraction=gd.getNextNumber();
			double corrCbWeight=gd.getNextNumber();
			double corrCrWeight=gd.getNextNumber();
			double correlationHighPassSigma=gd.getNextNumber();
			double correlationLowPassSigma=gd.getNextNumber();
			double noiseNormalizationSignaY=   gd.getNextNumber();
			double noiseNormalizationSignaCbCr=gd.getNextNumber();
			double contrastThreshold=gd.getNextNumber();
			double contrastThresholdDecrease=gd.getNextNumber();
			boolean useBinaryAlpha=     gd.getNextBoolean();

			double minTileFraction=0.01*gd.getNextNumber();
			int minDisparity=      (int)gd.getNextNumber();
			int maxDisparity=      (int)gd.getNextNumber();
			double relativeStep=        gd.getNextNumber();
			boolean probeAtLimits=      gd.getNextBoolean();
			
			double parentFraction=      gd.getNextNumber();
			double neighborsFraction=   gd.getNextNumber();
			double neighborsPower=      gd.getNextNumber();
			if (DEBUG_LEVEL>1) System.out.println("sides="+sides+" selves="+selves);

			long startTime = System.nanoTime();
			// setup tiles arrays
			int [] tilesStats=PIXEL_MAPPING.lastUsedInterSensor.setupTiles(
					tileSize,
					numLevels, // number of tile levels, each having half linear size (4..5)
	    			minDisparity, // -1
	    			maxDisparity);  // =size/2
	 			if (DEBUG_LEVEL>0) System.out.println("Using "+numLevels+" levels of tiles, top level is "+tilesStats[0]+"x"+tilesStats[1]+", largest is "+tileSize+
					" pixels square, smallest - "+(tileSize>>numLevels)+" pixels");
			// generate disparity tiles data
			if (DEBUG_LEVEL>1) System.out.println("Tiles initialized at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));

			int  tilesGenStats= PIXEL_MAPPING.lastUsedInterSensor.generateTileLevels(
					sides,
					selves,
					corrPhaseFraction, // double phaseCorrelationFraction,
					corrCbWeight,
					corrCrWeight,
					correlationHighPassSigma, //double correlationHighPassSigma,
					correlationLowPassSigma, //double correlationLowPassSigma,
					noiseNormalizationSignaY, //double noiseNormalizationSignaY,
					noiseNormalizationSignaCbCr, //double noiseNormalizationSignaCbCr,
					contrastThreshold, //double contrastThreshold,
					contrastThresholdDecrease, //double contrastThresholdDecrease, //Decrease threshold for smaller tiles as 1/2^(k*level), k=
					relativeStep,
					probeAtLimits,
					useBinaryAlpha, //boolean binaryAlpha,
					minTileFraction, // double minTileFraction,
					null, //DoubleFHT doubleFHT,
					THREADS_MAX, //int threadsMax,
					DEBUG_LEVEL); //int debugLevel)		
			
			if (DEBUG_LEVEL>0) System.out.println("Generated "+tilesGenStats+" non-empty tiles at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
//TODO:  Add modification based on neighbors -amplify by the layer above and tiles around 			
//TODO:  Generate disparity map by filtering out disparity values that contradict individual pixels 
// Maybe - a two-pass? First generate a range?
// Or process a tile, accumulating disparity arrays for each pixel and using neighbors - same as for tiles?			
/*			
			double [][] disparity=PIXEL_MAPPING.lastUsedInterSensor.collectDisparityFromTiles(DEBUG_LEVEL);
			String [] titles={"Best-shift","Best-contrast","BG-shift","BG-contrast","FG-shift","FG-contrast"};
			this.SDFA_INSTANCE.showArrays(
					disparity,
					PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
					PIXEL_MAPPING.lastUsedInterSensor.mapHeight,
					true,
			"disparity-PC"+((""+corrPhaseFraction).replace('.','_'))+"-THR"+((""+contrastThreshold).replace('.','_')),
			titles);
*/			
			// just temporary to visualize tiles - actual disparity will be processed diffirently
			
	    	double [][] disparity=PIXEL_MAPPING.lastUsedInterSensor.showDisparityFromTiles(
	    			DEBUG_LEVEL);
			String [] titles={"left disparity","right disparity","left contrast","right contrast",
					"left auto","right auto","left auto contrast","right auto contrast"};
			this.SDFA_INSTANCE.showArrays(
					disparity,
					PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
					PIXEL_MAPPING.lastUsedInterSensor.mapHeight,
					true,
			"disparity-PC"+((""+corrPhaseFraction).replace('.','_'))+"-THR"+((""+contrastThreshold).replace('.','_')),
			titles);
			if (DEBUG_LEVEL>0) System.out.println("Unfiltered shown at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
			if ((parentFraction>0) || (neighborsFraction>0)){
				PIXEL_MAPPING.lastUsedInterSensor.applyParentAndNeighbors(
						parentFraction,
						neighborsFraction,
						neighborsPower,
						DEBUG_LEVEL);
				if (DEBUG_LEVEL>0) System.out.println("Applied parent and neighbors at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));

				disparity=PIXEL_MAPPING.lastUsedInterSensor.showDisparityFromTiles(
						DEBUG_LEVEL);
				this.SDFA_INSTANCE.showArrays(
						disparity,
						PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
						PIXEL_MAPPING.lastUsedInterSensor.mapHeight,
						true,
						"filtered-PF"+(""+parentFraction).replace('.','_')+
						"-NF"+(""+neighborsFraction).replace('.','_')+
						"-NP"+(""+neighborsFraction).replace('.','_')+
						"-PC"+((""+corrPhaseFraction).replace('.','_'))+"-THR"+((""+contrastThreshold).replace('.','_')),
						titles);
				if (DEBUG_LEVEL>0) System.out.println("Filtered shown at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
			}
			if (DEBUG_LEVEL>0) System.out.println("Finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
			return;
		}
/* ======================================================================== */
		if       (label.equals("Disparity Section")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) {
				String msg="PIXEL_MAPPING is not initialized";
				System.out.println(msg);
				IJ.showMessage(msg);
				return;
			}
			if ((PIXEL_MAPPING.lastUsedInterSensor==null) || (PIXEL_MAPPING.lastUsedInterSensor.overlapImages==null)) {
				String msg="Inter-sensor data is not initialized (use \"Test Intermaps\"";
				System.out.println(msg);
				IJ.showMessage(msg);
				return;
			}
			
			GenericDialog gd=new GenericDialog("Disparity section");
			//gd.addString("Generate disparity for ")
			String [] choices= {"left image","right image","both images"};
			String [] choices1={"cross-correlate","auto-correlate","cross and auto"};
			gd.addChoice("Generate disparity for",choices, choices[2]);
			gd.addChoice("Croos/auto-correlate",  choices1,choices1[2]);
			gd.addNumericField("Section Y",780,0);
			gd.addNumericField("Weight edge-detected Y (0.0 ... 1.0)",0.0,3,5,"");
			gd.addNumericField("Weight of Cb (relative to Y) in correlation",0.5,3,5,"");
			gd.addNumericField("Weight of Cr (relative to Y) in correlation",0.5,3,5,"");
			gd.addNumericField("Noise level",5,3,5,"%");
			gd.addNumericField("Mix with edge disparity",2.0,3,5,"");
			
			gd.addMessage("parameters to remove artifacts caused by periodic patterns");
			gd.addNumericField("Minimal period to consider",3.0,3,5,"pix");
			gd.addNumericField("Minimal autocorrelation contrast at zero",0.5,3,5,"x");
			gd.addNumericField("Minimal period contrast",0.15,3,5,"x");
			gd.addNumericField("Minimal period contrast (local max) relative to absolute maximum",0.5,3,5,"x");
			gd.addNumericField("Blur rejection mask sigma (or 0 to skip)",1.0,3,5,"pix");
			gd.addNumericField("Rejection mask scale",10.0,3,5,"x");
			gd.showDialog();
			if (gd.wasCanceled()) return;
			int sides=               gd.getNextChoiceIndex()+1; // 1-2-3
			int selves=              gd.getNextChoiceIndex()+1; // 1-2-3
			int sectionY=      (int) gd.getNextNumber();
			double corrSobelWeight=  gd.getNextNumber();
			double corrCbWeight=     gd.getNextNumber();
			double corrCrWeight=     gd.getNextNumber();
			double noiseLevel=  0.01*gd.getNextNumber();
			double edgeMix=          gd.getNextNumber();
			
			double minPeriod=        gd.getNextNumber();
			double minZeroAuto=      gd.getNextNumber();
			double minPeriodContrast=gd.getNextNumber();
			double minAbsoluteFraction=gd.getNextNumber();
			double blurMaskSigma=    gd.getNextNumber();
			double scale=            gd.getNextNumber();

			long startTime = System.nanoTime();
			double [][][][] sections = new double [2][2][5][];
			for (int side=0;side<2;side++) if ((sides & (1<<side))!=0) for (int self=0;self<2;self++) if ((selves & (1<<self))!=0){
				
				sections[side][self][1]= PIXEL_MAPPING.lastUsedInterSensor.tilesSection(
		    			sectionY,
		    			side,
		    			self,
		    			DEBUG_LEVEL);
			}
 			if (DEBUG_LEVEL>0) System.out.println("tilesSection() done at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));

			for (int side=0;side<2;side++) if ((sides & (1<<side))!=0) for (int self=0;self<2;self++) if ((selves & (1<<self))!=0){
				sections[side][self][2]= PIXEL_MAPPING.lastUsedInterSensor.sectionPixelDifference(
		    			sectionY,
		    			side,
		    			self,
		    			corrSobelWeight,
		    			corrCbWeight,
		    			corrCrWeight,
		    			noiseLevel,
		    			DEBUG_LEVEL);
			}
 			if (DEBUG_LEVEL>0) System.out.println("sectionPixelDifference() done at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
			for (int side=0;side<2;side++) if ((sides & (1<<side))!=0) for (int self=0;self<2;self++) if ((selves & (1<<self))!=0){
				sections[side][self][3]= PIXEL_MAPPING.lastUsedInterSensor.alphaSection(
		    			sectionY,
		    			side,
		    			self,
		    			DEBUG_LEVEL);
			}
 			if (DEBUG_LEVEL>0) System.out.println("alphaSection() done at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));

			for (int side=0;side<2;side++) if ((sides & (1<<side))!=0) for (int self=0;self<2;self++) if ((selves & (1<<self))!=0){
				sections[side][self][4]= PIXEL_MAPPING.lastUsedInterSensor.sectionPixelSobel(
		    			sectionY,
		    			side,
		    			self,
		    			DEBUG_LEVEL);
			}
 			if (DEBUG_LEVEL>0) System.out.println("sectionPixelSobel() done at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));

    		int numSamples= PIXEL_MAPPING.lastUsedInterSensor.maxDisparity-PIXEL_MAPPING.lastUsedInterSensor.minDisparity+1;

			for (int side=0;side<2;side++) if ((sides & (1<<side))!=0) for (int self=0;self<2;self++) if ((selves & (1<<self))!=0){
				sections[side][self][0]=new double[sections[side][self][1].length];
				if (DEBUG_LEVEL>1) System.out.println("Combining layers for side="+side+" self="+self);
				double em=(self==0)?edgeMix:0.0;
				for (int i=0;i<sections[side][self][0].length;i++){
//					sections[side][self][4][i]=(sections[side][self][4][i]>0.0)?Math.sqrt(sections[side][self][4][i]):0.0;
//					sections[side][0][i]=sections[side][1][i]*sections[side][2][i]*sections[side][3][i]*
//					(1+	edgeMix*sections[side][4][i]*sections[side][3][i]);
					int iX=i/numSamples;
					int index=sectionY*PIXEL_MAPPING.lastUsedInterSensor.mapWidth+iX;
					double thisEdge=PIXEL_MAPPING.lastUsedInterSensor.sobelY[side][index];
					sections[side][self][0][i]=sections[side][self][1][i]*sections[side][self][2][i]*sections[side][self][3][i]*
					(1+	em*thisEdge); // just increase weight of this image edges
				}
			}
			double [][] noPeriodic={null,null};
			if (selves==3) {
				for (int side=0;side<2;side++) if ((sides & (1<<side))!=0){

					if (DEBUG_LEVEL>1) System.out.println("Removing periodic ambiguity for image"+side);
					noPeriodic[side]=PIXEL_MAPPING.lastUsedInterSensor.removePeriodicSection(
							sections[side][0][0],
							sections[side][1][0],
							minPeriod,
							minZeroAuto,
							minPeriodContrast,
							minAbsoluteFraction,
							blurMaskSigma,
							scale,
							DEBUG_LEVEL	);
				}
				String [] titlesLR={"Left","Right"};
				this.SDFA_INSTANCE.showArrays(
						noPeriodic,
						PIXEL_MAPPING.lastUsedInterSensor.maxDisparity-PIXEL_MAPPING.lastUsedInterSensor.minDisparity+1,
						PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
						true,
						"Noperiodic-Y"+sectionY+"-BS"+((""+blurMaskSigma).replace('.','_'))+
						"-SCALE"+((""+scale).replace('.','_'))+"-MPC"+((""+minPeriodContrast).replace('.','_')),
						titlesLR);
			}
			
			
			
			String [] titles={"Product section","Tiles section","pixdiff section","alpha","Sobel"};
			for (int side=0;side<2;side++) if ((sides & (1<<side))!=0) for (int self=0;self<2;self++) if ((selves & (1<<self))!=0){
				if (DEBUG_LEVEL>1) System.out.println("side="+side+" self="+self);
				if (sections[side][self]==null){
					System.out.println("sections["+side+"]["+self+"]=null");
					
				} else {
					if (DEBUG_LEVEL>1) {
						System.out.println("sections["+side+"]["+self+"].length="+sections[side][self].length);
						for (int i=0;i<sections[side][self].length;i++){
							if (sections[side][self][i]==null){
								System.out.println("sections["+side+"]["+self+"]["+i+"]=null");
								
							} else {
								System.out.println("sections["+side+"]["+self+"]["+i+"].length="+sections[side][self][i].length);
							}							
						}
					}
				}
				this.SDFA_INSTANCE.showArrays(
						sections[side][self],
						PIXEL_MAPPING.lastUsedInterSensor.maxDisparity-PIXEL_MAPPING.lastUsedInterSensor.minDisparity+1,
						PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
						true,
						((side==0)?"Left":"Right")+((self==0)?"-cross":"-auto")+"-Y"+sectionY+"-N"+((""+noiseLevel).replace('.','_'))+
						"-Cb"+((""+corrCbWeight).replace('.','_'))+"-Cr"+((""+corrCrWeight).replace('.','_')),
						titles);
			}
			if (DEBUG_LEVEL>0) System.out.println("Finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
			return;
		}
/* ======================================================================== */
		if       (label.equals("Disparity Map")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) {
				String msg="PIXEL_MAPPING is not initialized";
				System.out.println(msg);
				IJ.showMessage(msg);
				return;
			}
			if ((PIXEL_MAPPING.lastUsedInterSensor==null) || (PIXEL_MAPPING.lastUsedInterSensor.overlapImages==null)) {
				String msg="Inter-sensor data is not initialized (use \"Test Intermaps\"";
				System.out.println(msg);
				IJ.showMessage(msg);
				return;
			}
			
			GenericDialog gd=new GenericDialog("Disparity Map");
			//gd.addString("Generate disparity for ")
			String [] choices={"left image","right image","both images"};
			String [] choices1={"cross-correlate","auto-correlate","cross and auto"};
			gd.addChoice("Generate disparity for",choices, choices[2]);
			gd.addChoice("Croos/auto-correlate",  choices1,choices1[2]);
			gd.addNumericField("Weight edge-detected Y (0.0 ... 1.0)",0.5,3,5,"");
			gd.addNumericField("Weight of Cb (relative to Y) in correlation",0.5,3,5,"");
			gd.addNumericField("Weight of Cr (relative to Y) in correlation",0.5,3,5,"");
			gd.addNumericField("Noise level",2,3,5,"%");
			gd.showDialog();
			if (gd.wasCanceled()) return;
			int sides=               gd.getNextChoiceIndex()+1; // 1-2-3
			int selves=              gd.getNextChoiceIndex()+1; // 1-2-3
			double corrSobelWeight=  gd.getNextNumber();
			double corrCbWeight=     gd.getNextNumber();
			double corrCrWeight=     gd.getNextNumber();
			double noiseLevel=  0.01*gd.getNextNumber();

			long startTime = System.nanoTime();
			double [][][] disparityMap = new double [2][2][];

			for (int side=0;side<2;side++) if ((sides & (1<<side))!=0) for (int self=0;self<2;self++) if ((selves & (1<<self))!=0){
				disparityMap[side][self]= PIXEL_MAPPING.lastUsedInterSensor.getFrameDisparity(
		    			side,
		    			self,
		    			corrSobelWeight,
		    			corrCbWeight,
		    			corrCrWeight,
		    			noiseLevel,
		    			DEBUG_LEVEL);
			}
 			if (DEBUG_LEVEL>0) System.out.println("getFrameDisparity() done at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
			for (int side=0;side<2;side++) if ((sides & (1<<side))!=0) for (int self=0;self<2;self++) if ((selves & (1<<self))!=0){
 				this.SDFA_INSTANCE.showArrays(
 						disparityMap[side][self],
 						PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
 						PIXEL_MAPPING.lastUsedInterSensor.mapHeight,
 						((side==0)?"Left":"Right")+((self==0)?"-cross":"-auto")+"-N"+((""+noiseLevel).replace('.','_'))+
 						"-Cb"+((""+corrCbWeight).replace('.','_'))+"-Cr"+((""+corrCrWeight).replace('.','_')));
 			}
 			if (DEBUG_LEVEL>0) System.out.println("Finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
 			return;
		}

/* ======================================================================== */
		if       (label.equals("Create Ambiguity")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) {
				String msg="PIXEL_MAPPING is not initialized";
				System.out.println(msg);
				IJ.showMessage(msg);
				return;
			}
			if ((PIXEL_MAPPING.lastUsedInterSensor==null) || (PIXEL_MAPPING.lastUsedInterSensor.overlapImages==null)) {
				String msg="Inter-sensor data is not initialized (use \"Test Intermaps\"";
				System.out.println(msg);
				IJ.showMessage(msg);
				return;
			}
			
			GenericDialog gd=new GenericDialog("Create ambiguity map");
			//gd.addString("Generate disparity for ")
			String [] choices={"left image","right image","both images"};
			gd.addChoice("Generate disparity for",choices,choices[2]);
			gd.addNumericField("Weight edge-detected Y (0.0 ... 1.0)",0.5,3,5,"");
			gd.addNumericField("Weight of Cb (relative to Y) in correlation",0.5,3,5,"");
			gd.addNumericField("Weight of Cr (relative to Y) in correlation",0.5,3,5,"");
			gd.addNumericField("Noise level",5,3,5,"%");
			gd.addMessage("Parameters to remove artifacts caused by periodic patterns");
			gd.addCheckbox("Filter periodic patterns artifacts",true);
			gd.addNumericField("Minimal period to consider",3.0,3,5,"pix");
			gd.addNumericField("Minimal autocorrelation contrast at zero",0.5,3,5,"x");
			gd.addNumericField("Minimal period contrast",0.15,3,5,"x");
			gd.addNumericField("Minimal period contrast (local max) relative to absolute maximum",0.5,3,5,"x");
			gd.addNumericField("Blur rejection mask sigma (or 0 to skip)",2.0,3,5,"pix");
			gd.addNumericField("Rejection mask scale",50.0,3,5,"x");

			gd.showDialog();
			if (gd.wasCanceled()) return;
			int sides=               gd.getNextChoiceIndex()+1; // 1-2-3
			double corrSobelWeight=  gd.getNextNumber();
			double corrCbWeight=     gd.getNextNumber();
			double corrCrWeight=     gd.getNextNumber();
			double noiseLevel=  0.01*gd.getNextNumber();

			boolean removePeriodic=  gd.getNextBoolean();
			double minPeriod=        gd.getNextNumber();
			double minZeroAuto=        gd.getNextNumber();
			double minPeriodContrast=gd.getNextNumber();
			double minAbsoluteFraction=gd.getNextNumber();
			double blurMaskSigma=    gd.getNextNumber();
			double scale=            gd.getNextNumber();


			long startTime = System.nanoTime();

			for (int side=0;side<2;side++) if ((sides & (1<<side))!=0){
				// TODO: - mayebe provide 3 arrays: auro (from fft only), cross (from fft only) - use them to create a mask, then apply mask to full cross
				PIXEL_MAPPING.lastUsedInterSensor.createAmbiguityMaps( // result will be stored in class instance
		    			side,
		    			corrSobelWeight,
		    			corrCbWeight,
		    			corrCrWeight,
		    			noiseLevel,
		    			removePeriodic,
						minPeriod,
						minZeroAuto,
						minPeriodContrast,
						minAbsoluteFraction,
						blurMaskSigma,
						scale,
		    			DEBUG_LEVEL);
			}
 			if (DEBUG_LEVEL>0) System.out.println("createAmbiguityMaps() done at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
 			String [] titles = {"Disparity","Strength","Second Disparity","Second Strength"};
			for (int side=0;side<2;side++) if ((sides & (1<<side))!=0){
				this.SDFA_INSTANCE.showArrays(
						PIXEL_MAPPING.lastUsedInterSensor.ambiguityMap[side],
						PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
						PIXEL_MAPPING.lastUsedInterSensor.mapHeight,
						true,
						"Amb"+((side==0)?"Left":"Right")+"-N"+((""+noiseLevel).replace('.','_'))+
						"-Cb"+((""+corrCbWeight).replace('.','_'))+"-Cr"+((""+corrCrWeight).replace('.','_')),
						titles);
			}
			if (DEBUG_LEVEL>0) System.out.println("Finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));

			return;
		}
/* ======================================================================== */
		if       (label.equals("Initial Resolve")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) {
				String msg="PIXEL_MAPPING is not initialized";
				System.out.println(msg);
				IJ.showMessage(msg);
				return;
			}
			if ((PIXEL_MAPPING.lastUsedInterSensor==null) || (PIXEL_MAPPING.lastUsedInterSensor.ambiguityMap==null)) {
				String msg="Ambiguitty maps do not exist (use \"Create Ambiguity\"";
				System.out.println(msg);
				IJ.showMessage(msg);
				return;
			}
			
			GenericDialog gd=new GenericDialog("Initial resolve ambiguous maps");
			//gd.addString("Generate disparity for ")
			String [] choices={"left image","right image","both images"};
			gd.addChoice("Generate disparity for",choices,choices[2]);
			gd.addNumericField("Contrast threshold",0.5,3,5,"");
			gd.addNumericField("Maximal allowed second maximum as a fraction of the absolute maximum",20.0,3,5,"%");
			gd.addNumericField("Edges bonus (multiply strength by 1+edgesBonus*edges",10.0,3,5,"");

			gd.showDialog();
			if (gd.wasCanceled()) return;
			int sides=               gd.getNextChoiceIndex()+1; // 1-2-3
			double contrastThreshold=     gd.getNextNumber();
			double minSecondFrac=     0.01*gd.getNextNumber();
			double edgesBonus=             gd.getNextNumber();

			long startTime = System.nanoTime();
			float [][][] mapAndContrast=new float [2][2][];
			for (int i=0;i<mapAndContrast.length;i++) for (int j=0;j<mapAndContrast[i].length;j++) mapAndContrast[i][j]=null;
			
			for (int side=0;side<2;side++) if ((sides & (1<<side))!=0){
				mapAndContrast[side][0]= PIXEL_MAPPING.lastUsedInterSensor.initialResolveMaps( // result will be stored in class instance
		    			side,
		    			contrastThreshold,
		    			minSecondFrac,
		    			edgesBonus,
		    			DEBUG_LEVEL);
				mapAndContrast[side][1]= PIXEL_MAPPING.lastUsedInterSensor.getResolvedState( // result will be stored in class instance
		    			side,
		    			contrastThreshold,
		    			minSecondFrac,
		    			edgesBonus,
		    			32.0, // so the resultys will be 0/32/64 matching disparity values 
		    			DEBUG_LEVEL);
				
			}
 			if (DEBUG_LEVEL>0) System.out.println("initialResolveMaps() and getResolvedState done at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
 			String [] titles = {"Disparity","Ambiguity"};
			for (int side=0;side<2;side++) if ((sides & (1<<side))!=0){
				this.SDFA_INSTANCE.showArrays(
						mapAndContrast[side],
						PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
						PIXEL_MAPPING.lastUsedInterSensor.mapHeight,
						true,
				((side==0)?"Left":"Right")+"-CT"+((""+contrastThreshold).replace('.','_'))+
				"-MSF"+((""+minSecondFrac).replace('.','_')),
				titles);
			}
			if (DEBUG_LEVEL>0) System.out.println("Finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
			return;
		}
/* ======================================================================== */
		if       (label.equals("Ambiguity Resolve")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (PIXEL_MAPPING==null) {
				String msg="PIXEL_MAPPING is not initialized";
				System.out.println(msg);
				IJ.showMessage(msg);
				return;
			}
			if ((PIXEL_MAPPING.lastUsedInterSensor==null) || (PIXEL_MAPPING.lastUsedInterSensor.resolvedMap==null)) {
				String msg="Ambiguitty maps do not exist (use \"Initial Resolve\"";
				System.out.println(msg);
				IJ.showMessage(msg);
				return;
			}

			double minExistentAlpha=0.5;
			double minExistentStrength=1.0;
			double weightCb=0.5; 
			double weightCr=0.5;
			double maxToneDiff=0.02;
			double maxDisparityDifference=0.5;
			double minNewStrength=1.0;
			boolean repeatWhilePossible=false;
			boolean showResult=true;
			int []  applicationStats={-1,-1};
			int sides=3;
			while (true) {
				GenericDialog gd=new GenericDialog("Initial resolve ambiguous maps");
				if (applicationStats[0]>=0){
					gd.addMessage("Resolved "+applicationStats[0]+" ambiguous pixels");
					gd.addMessage("Filtered out "+applicationStats[1]+" ambiguous pixels because of constraints");
				}
				//gd.addString("Generate disparity for ")
				String [] choices={"left image","right image","both images"};
				gd.addChoice("Generate disparity for",choices,choices[sides-1]);
				gd.addNumericField("Minimal combined alpha for parent (already resolved) pixels",minExistentAlpha,3,5,"");
				gd.addNumericField("Minimal strength for parent (already resolved) pixels",minExistentStrength,3,5,"");
				gd.addNumericField("Weight of Cb (relative to Y) in correlation",weightCb,3,5,"");
				gd.addNumericField("Weight of Cr (relative to Y) in correlation",weightCr,3,5,"");
				gd.addNumericField("Maximal weighted tone difference between the parent and new pixels",100.0*maxToneDiff,3,5,"%");
				gd.addNumericField("Maximal disparity difference between the parent and new pixels",maxDisparityDifference,3,5,"pix");
				gd.addNumericField("Minimal strength of the new pixels",minNewStrength,3,5,"");
				gd.addCheckbox    ("Repeat resolution while it produces new points", repeatWhilePossible);
				gd.addCheckbox    ("Show result image(s)", showResult);
				gd.showDialog();
				if (gd.wasCanceled()) return;
				sides=               gd.getNextChoiceIndex()+1; // 1-2-3
				minExistentAlpha=gd.getNextNumber();
				minExistentStrength=gd.getNextNumber();
				weightCb=gd.getNextNumber(); 
				weightCr=gd.getNextNumber();
				maxToneDiff=0.01*gd.getNextNumber();
				maxDisparityDifference=gd.getNextNumber();
				minNewStrength=gd.getNextNumber();
				repeatWhilePossible=gd.getNextBoolean();
				showResult=gd.getNextBoolean();
				long startTime = System.nanoTime();
				int totalApplied=0;
				while (true) {
					applicationStats[0]=0;
					applicationStats[1]=0;
					for (int side=0;side<2;side++) if ((sides & (1<<side))!=0){
						int [] applicationSubStats=PIXEL_MAPPING.lastUsedInterSensor.resolveAmbiguityStep(
								side,
								minExistentAlpha,
								minExistentStrength,
								weightCb,
								weightCr,
								maxToneDiff, // using weighted distance)
								maxDisparityDifference,
								minNewStrength,
								DEBUG_LEVEL);
						if (DEBUG_LEVEL>0) {
							System.out.println(((side==0)?"Left":"Right")+" image: " + applicationSubStats[0]+" pixels resolved, "+
									applicationSubStats[1]+" pixels filtered at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
						}
						applicationStats[0]+=applicationSubStats[0];
						applicationStats[1]+=applicationSubStats[1];
					}
					totalApplied+=applicationStats[0];
					if (!repeatWhilePossible || (applicationStats[0]==0)) break;
				}
				applicationStats[0]=totalApplied;
				if (showResult){
//					String [] titles = {"Disparity","Ambiguity"};
					for (int side=0;side<2;side++) if ((sides & (1<<side))!=0){
						this.SDFA_INSTANCE.showArrays(
								PIXEL_MAPPING.lastUsedInterSensor.resolvedMap[side],
								PIXEL_MAPPING.lastUsedInterSensor.mapWidth,
								PIXEL_MAPPING.lastUsedInterSensor.mapHeight,
								//							true,
								((side==0)?"Left":"Right")+"-resolved");
					}
					//				if (DEBUG_LEVEL>0) System.out.println("Finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
				}
			}

		}
/* ======================================================================== */
		if       (label.equals("Pattern Flat-Field")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			
			int series=LENS_DISTORTIONS.refineParameters.showDialog(
	    			"Select Lens Distortion Residual Compensation Parameters",
	    			0x100000,
	    			((LENS_DISTORTIONS.seriesNumber>=0)?LENS_DISTORTIONS.seriesNumber:0),
	    			LENS_DISTORTIONS.patternParameters.averageRGB);
			if (series<0) return;
			LENS_DISTORTIONS.correctPatternFlatField(true); // boolean enableShow
			return;
		}	
/* ======================================================================== */
		if       (label.equals("Remove Specular")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			int series=LENS_DISTORTIONS.refineParameters.showDialog(
					"Removing specular reflections from the target",
	    			0x400000,
	    			((LENS_DISTORTIONS.seriesNumber>=0)?LENS_DISTORTIONS.seriesNumber:0),
	    			LENS_DISTORTIONS.patternParameters.averageRGB);
			if (series<0) return;
/*
    		double highPassSigma=10.0;
			double diffFromAverageThreshold=0.1;
			int numIter=2;
			boolean applyNewWeights=true;
			boolean positiveDiffOnly=true;
			GenericDialog gd=new GenericDialog("Removing specular reflections from the target");
			gd.addCheckbox    ("Process only positive difference fr4om average",               positiveDiffOnly);
			gd.addNumericField("High-pass sigma for difference from average (to detect specular)", highPassSigma,2,6,"nodes");
			gd.addNumericField("Difference from average threshold", diffFromAverageThreshold,3,5,"");
			gd.addNumericField("Number of iterations for calculating average", numIter,0);
			gd.addCheckbox    ("Apply new weights",                            applyNewWeights);
			
			gd.showDialog();
			if (gd.wasCanceled()) return;
			positiveDiffOnly=         gd.getNextBoolean();
			highPassSigma=            gd.getNextNumber();
			diffFromAverageThreshold= gd.getNextNumber();
			numIter=            (int) gd.getNextNumber();
			applyNewWeights=          gd.getNextBoolean();
			*/

			LENS_DISTORTIONS.removeSpecular (
					LENS_DISTORTIONS.refineParameters.specularPositiveDiffOnly,
					LENS_DISTORTIONS.refineParameters.specularHighPassSigma,
					LENS_DISTORTIONS.refineParameters.specularLowPassSigma,
					LENS_DISTORTIONS.refineParameters.specularDiffFromAverageThreshold,
					LENS_DISTORTIONS.refineParameters.specularNumIter,
					LENS_DISTORTIONS.refineParameters.specularApplyNewWeights,
					LENS_DISTORTIONS.refineParameters.specularShowDebug);
			return;
		}	
/* ======================================================================== */
		if       (label.equals("Flat-Field")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			PATTERN_PARAMETERS.debugLevel=DEBUG_LEVEL;
			if (LENS_DISTORTIONS.fittingStrategy==null){
				IJ.showMessage("Distortion Fitting strategy is not initialized - create it with" +
				"\"New Strategy\", \"Edit Strategy\" or \"Restore Strategy\"");
				return;
			}
			if (DISTORTION_CALIBRATION_DATA==null){
				IJ.showMessage("Distortion Calibration data is not initialized - create it with"+
				"\"SelectGrid Files\" or \"Restore Calibration\"");
				return;
			}
			DISTORTION_CALIBRATION_DATA.debugLevel=DEBUG_LEVEL;
			LENS_DISTORTIONS.fittingStrategy.debugLevel=DEBUG_LEVEL;
	    	if (LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData.eyesisCameraParameters==null){
	    		String msg="Eyesis camera parameters (and sensor dimensions) are not defined";
	    		IJ.showMessage("Error",msg);
	    		throw new IllegalArgumentException (msg);
	    	}
//	    	int series=refineParameters.showDialog("Select Lens Distortion Residual Compensation Parameters", 0x1efff, (this.seriesNumber>=0)?this.seriesNumber:0);
	    	int series=LENS_DISTORTIONS.refineParameters.showDialog(
	    			"Select Sensor and Target Flat-Field Correction Parameters",
	    			0x794ff1,
//	    			/0x94ff1,
	    			((LENS_DISTORTIONS.seriesNumber>=0)?LENS_DISTORTIONS.seriesNumber:0),
	    			LENS_DISTORTIONS.patternParameters.averageRGB); // averageRGB - only for target flat-field correction
	    	if (series<0) return;
	    	LENS_DISTORTIONS.seriesNumber=series;
			long 	  startTime=System.nanoTime();    	
	    	for (int nIteration = 0; nIteration<LENS_DISTORTIONS.refineParameters.repeatFlatFieldSensor;nIteration++){
	    		if (this.SYNC_COMMAND.stopRequested.get()>0){
	    			System.out.println("Stop requested, command aborted, some changes may already be applied");
	    			break;
	    		}
	    		if(DEBUG_LEVEL>0) System.out.println("Calculating sensors flat-field correction, iteration #"+
	    				nIteration+" ("+(nIteration+1)+" of "+LENS_DISTORTIONS.refineParameters.repeatFlatFieldSensor+
	    				", started at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
				boolean OK=LENS_DISTORTIONS.modifyPixelCorrection(
						(nIteration>=(LENS_DISTORTIONS.refineParameters.repeatFlatFieldSensor-1)), // enable show - only on last iteration
						THREADS_MAX,
						UPDATE_STATUS,
						DEBUG_LEVEL);
				if (!OK){
	    			System.out.println("Stop requested in iternal loop, command aborted, some changes may already be applied");
	    			break;
					
				}
	    		if(DEBUG_LEVEL>0) System.out.println("Calculating target flat-field correction, iteration #"+
	    				nIteration+" ("+(nIteration+1)+" of "+LENS_DISTORTIONS.refineParameters.repeatFlatFieldSensor+
	    				", started at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
				LENS_DISTORTIONS.correctPatternFlatField(
						(nIteration>=(LENS_DISTORTIONS.refineParameters.repeatFlatFieldSensor-1))); // enable show - only on last iteration
				if (LENS_DISTORTIONS.refineParameters.specularApplyNewWeights){
					LENS_DISTORTIONS.removeSpecular (
							LENS_DISTORTIONS.refineParameters.specularPositiveDiffOnly,
							LENS_DISTORTIONS.refineParameters.specularHighPassSigma,
							LENS_DISTORTIONS.refineParameters.specularLowPassSigma,
							LENS_DISTORTIONS.refineParameters.specularDiffFromAverageThreshold,
							LENS_DISTORTIONS.refineParameters.specularNumIter,
							true, // LENS_DISTORTIONS.refineParameters.specularApplyNewWeights);
							LENS_DISTORTIONS.refineParameters.specularShowDebug);
				}
	    	}
    		if(DEBUG_LEVEL>0) System.out.println("Flat-field correction ("+LENS_DISTORTIONS.refineParameters.repeatFlatFieldSensor+
    				" iterations) finished in "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+" seconds");
			return;
		}		
		
		
/* ======================================================================== */
		if       (label.equals("Pattern Flat-Field0")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			LENS_DISTORTIONS.debugLevel=DEBUG_LEVEL;
			
			
			GenericDialog gd=new GenericDialog("Pattern Flat Field parameters");
			gd.addNumericField("Fitting series number (to select images), negative - use all enabled images", -1,0);
			gd.addNumericField("Reference station number (unity target brightness)", 0,0);
			gd.addNumericField("Shrink sensor mask",    100.0, 1,6,"sensor pix");
			gd.addNumericField("Non-vignetted radius", 1000.0, 1,6,"sensor pix");
			gd.addNumericField("Minimal alpha",           1.0, 3,7,"%");
			gd.addNumericField("Minimal contrast (occlusion detection)", 0.1, 3,7,"(0 .. ~0.8");
			gd.addNumericField("Minimal alpha for accumulation", 1.0, 3,7,"%");
			gd.addNumericField("Shrink pattern for matching", 2.0, 3,7,"grid nodes");
			gd.addNumericField("Maximal relative difference between nodes", 10.0, 3,7,"%");
			gd.addNumericField("Shrink pattern border", 2, 0,3,"grid nodes");
			gd.addNumericField("Fade pattern border", 2.0, 3,7,"grid nodes");
			gd.addMessage("Update pattern white balance (if the illumination is yellowish, increase red and green here)");
			
			gd.addNumericField("Average grid RED   (1.0 for white)",  LENS_DISTORTIONS.patternParameters.averageRGB[0], 3,5,"x"); //
			gd.addNumericField("Average grid GREEN (1.0 for white)",  LENS_DISTORTIONS.patternParameters.averageRGB[1], 3,5,"x"); //
			gd.addNumericField("Average grid BLUE  (1.0 for white)",  LENS_DISTORTIONS.patternParameters.averageRGB[2], 3,5,"x"); //
			gd.addCheckbox("Reset pattern mask",               true);
			gd.addCheckbox("Show non-vignetting sensor masks", false);
			gd.addCheckbox("Show per-sensor patterns",         false);
			gd.addCheckbox("Show result mask",                 true);
			gd.addCheckbox("Apply pattern flat field and mask",true);
			gd.addCheckbox("Use interpolation for sensor correction",true);
			gd.addNumericField("Suspect occlusion only if grid is missing in the area where sensor mask is above this threshold",15, 3,7,"%");
			gd.addNumericField("Expand suspected occlusion  area", 2, 0,3,"grid nodes");
			gd.addNumericField("Fade grid on image (occlusion handling)", 2.0, 3,7,"grid nodes");
			gd.addCheckbox("Ignore existent sensor flat-field calibration",false);
			gd.addCheckbox("Use only selected channels",false);
			
			gd.showDialog();
			if (gd.wasCanceled()) return;
			int serNumber=            (int) gd.getNextNumber();
			int referenceStation=     (int) gd.getNextNumber();
    		double shrink=                  gd.getNextNumber();
    		double nonVignettedRadius =     gd.getNextNumber();
    		double minimalAlpha =      0.01*gd.getNextNumber();
    		
    		double minimalContrast =        gd.getNextNumber();
    		double minimalAccumulate = 0.01*gd.getNextNumber();
    		double shrinkForMatching =      gd.getNextNumber();
    		double maxRelDiff =        0.01*gd.getNextNumber();
    		int shrinkMask=           (int) gd.getNextNumber();
    		double fadeBorder =             gd.getNextNumber();
    		LENS_DISTORTIONS.patternParameters.averageRGB[0]=gd.getNextNumber();
    		LENS_DISTORTIONS.patternParameters.averageRGB[1]=gd.getNextNumber();
    		LENS_DISTORTIONS.patternParameters.averageRGB[2]=gd.getNextNumber();
    		boolean resetMask=      gd.getNextBoolean();
    		boolean showSensorMasks=gd.getNextBoolean();
    		boolean showIndividual= gd.getNextBoolean();
    		boolean showResult=     gd.getNextBoolean();
    		boolean applyResult=    gd.getNextBoolean();
    		boolean useInterpolate= gd.getNextBoolean();
    		double maskThresholdOcclusion=0.01*gd.getNextNumber();
    		int shrinkOcclusion= (int)gd.getNextNumber();
    		double fadeOcclusion=   gd.getNextNumber();

    		boolean ignoreSensorFlatField= gd.getNextBoolean();
    		boolean useSelectedChannels= gd.getNextBoolean();
    		boolean [] selectedChannels=null;
    		
    		LENS_DISTORTIONS.patternParameters.updateNumStations(LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getNumStations());
    		if (useSelectedChannels && (ABERRATIONS_PARAMETERS!=null))selectedChannels=ABERRATIONS_PARAMETERS.selectedChannels;
    	    double [][] masks= LENS_DISTORTIONS.nonVignettedMasks(
    	    		shrink,
    	    		nonVignettedRadius,
    	    		minimalAlpha);
    	    if (selectedChannels!=null){
    	    	for (int nChn=0;nChn<masks.length;nChn++) if ((nChn<selectedChannels.length)&&!selectedChannels[nChn]) masks[nChn]=null;
    	    }
    	    
			if (showSensorMasks) this.SDFA_INSTANCE.showArrays( //java.lang.ArrayIndexOutOfBoundsException: 313632
					masks,
					LENS_DISTORTIONS.pixelCorrectionWidth/ LENS_DISTORTIONS.pixelCorrectionDecimation,
					LENS_DISTORTIONS.pixelCorrectionHeight/LENS_DISTORTIONS.pixelCorrectionDecimation,
					true,
					"nonVinetting masks");
			double [][][][] sensorGrids=LENS_DISTORTIONS.calculateGridFlatField(
					serNumber,
					masks,
					minimalContrast,
					minimalAccumulate,
					useInterpolate,
				    maskThresholdOcclusion, // suspect occlusion only if grid is missing in the area where sensor mask is above this threshold
		    		shrinkOcclusion,
		    		fadeOcclusion,
					ignoreSensorFlatField);
			double [][][] geometry= LENS_DISTORTIONS.patternParameters.getGeometry();
			if (showIndividual){
				for (int station=0;station<sensorGrids.length;station++) if (sensorGrids[station]!=null){
					for (int i=0;i<sensorGrids[station].length;i++) if (sensorGrids[station][i]!=null){
						this.SDFA_INSTANCE.showArrays(
								sensorGrids[station][i],
								geometry[0].length,
								geometry.length,
								true,
								"chn"+i+":"+station+"-pattern");
					}
				}
			}
			double [][][][] patternArray= LENS_DISTORTIONS.combineGridFlatField(
					referenceStation,
					sensorGrids,
					shrinkForMatching,
					resetMask,
					maxRelDiff,
				    shrinkMask,
					fadeBorder);
			if (showResult) {
				String [] titles={"Alpha","Red","Green","Blue","Number of images used"};
				for (int station=0;station<patternArray.length;station++) if (patternArray[station]!=null){
					for (int nView=0;nView<patternArray[station].length;nView++) if (patternArray[station][nView]!=null){
						this.SDFA_INSTANCE.showArrays(
								patternArray[station][nView],
								geometry[0].length,
								geometry.length,
								true,
								"St"+station+"_V"+nView+"_Pattern_Colors "+maxRelDiff,
								titles);
					}
				}
			}
			if (applyResult) LENS_DISTORTIONS.applyGridFlatField(patternArray); // {alpha, red,green,blue, number of images used}[pixel_index]
			return;
		}		
//"Pattern Flat-Field"
//"Generate & Save Equirectangular" 
/* ======================================================================== */
		if       (label.equals("Configure aberrations")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			ABERRATIONS_PARAMETERS.showDialog("Aberration processing parameters",LENS_DISTORTIONS); //LENS_DISTORTIONS==null OK
			return;
		}
/* ======================================================================== */
		if       (label.equals("Select Channels")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			ABERRATIONS_PARAMETERS.selectChannelsToProcess("Select channels to process",LENS_DISTORTIONS); //LENS_DISTORTIONS==null OK
			return;
		}
/* ======================================================================== */
		if       (label.equals("Partial Kernels")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			EYESIS_ABERRATIONS.setDistortions(LENS_DISTORTIONS);
			if (EYESIS_ABERRATIONS.aberrationParameters.partialKernelDirectory.length()>0){
				File dFile=new File(EYESIS_ABERRATIONS.aberrationParameters.partialKernelDirectory);
				if (!dFile.isDirectory() &&  !dFile.mkdirs()) {
					String msg="Failed to create directory "+EYESIS_ABERRATIONS.aberrationParameters.partialKernelDirectory;
					IJ.showMessage("Warning",msg);
		    		System.out.println("Warning: "+msg);
		    		EYESIS_ABERRATIONS.aberrationParameters.partialKernelDirectory=""; // start over with selecting directory
				}
			}
			String configPath=EYESIS_ABERRATIONS.aberrationParameters.selectPartialKernelDirectory(
					true,
					EYESIS_ABERRATIONS.aberrationParameters.partialKernelDirectory,
					true);
			if (configPath==null){
	    		String msg="No partial kernel directory selected, command aborted";
	    		System.out.println("Warning: "+msg);
	    		IJ.showMessage("Warning",msg);
	    		return;
			}
			configPath+=Prefs.getFileSeparator()+"config-partial";
			try {
				saveTimestampedProperties(
						configPath,      // full path or null
						null, // use as default directory if path==null 
						true,
						PROPERTIES);
				
			} catch (Exception e){
	    		String msg="Failed to save configuration to "+configPath+", command aborted";
	    		System.out.println("Error: "+msg);
	    		IJ.showMessage("Error",msg);
	    		return;
			}
			
			
			EYESIS_ABERRATIONS.createPartialKernels(
					this.SYNC_COMMAND.stopRequested,
					MAP_FFT_SIZE, // scanImageForPatterns:FFT size //int             mapFFTsize, // scanImageForPatterns:FFT size
					FFT_OVERLAP, ////int            fft_overlap,
					FFT_SIZE, //int               fft_size,
					PSF_SUBPIXEL, // //int           PSF_subpixel, 
					OTF_FILTER, // //OTFFilterParameters otfFilterParameters,
					PSF_PARS, //PSFParameters psfParameters,
					INVERSE.dSize, //int          PSFKernelSize, // size of square used in the new map (should be multiple of map step)

					GAUSS_WIDTH, // double       gaussWidth,  // ** NEW
					MULTIFILE_PSF, //multiFilePSF
					DISTORTION, //MatchSimulatedPattern.DistortionParameters distortionParameters, //
					PATTERN_DETECT, // MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
					SIMUL, //SimulationPattern.SimulParameters  simulParameters,
					COMPONENTS, //boolean equalizeGreens,
					THREADS_MAX, // int threadsMax,
					UPDATE_STATUS,
					DISTORTION.loop_debug_level, // int loopDebugLevel, // debug level used inside loops
					DEBUG_LEVEL //int debugLevel
					);
///getPath()
			saveTimestampedProperties( // save config again
					configPath,      // full path or null
					null, // use as default directory if path==null 
					true,
					PROPERTIES);
			LENS_DISTORTIONS.fittingStrategy.distortionCalibrationData.
			saveTimestampedToXML(EYESIS_ABERRATIONS.aberrationParameters.partialKernelDirectory+Prefs.getFileSeparator()+"calib_", null); 

			return;
		}
		
/* ======================================================================== */
		if       (label.equals("Combine Kernels")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			EYESIS_ABERRATIONS.setDistortions(LENS_DISTORTIONS);
			if (EYESIS_ABERRATIONS.aberrationParameters.psfKernelDirectory.length()>0){
				File dFile=new File(EYESIS_ABERRATIONS.aberrationParameters.psfKernelDirectory);
				if (!dFile.isDirectory() &&  !dFile.mkdirs()) {
					String msg="Failed to create directory "+EYESIS_ABERRATIONS.aberrationParameters.psfKernelDirectory;
					IJ.showMessage("Warning",msg);
					System.out.println("Warning: "+msg);
					EYESIS_ABERRATIONS.aberrationParameters.psfKernelDirectory=""; // start over with selecting directory
				}
			}
			String configPath=EYESIS_ABERRATIONS.aberrationParameters.selectPSFKernelDirectory(
					true,
					EYESIS_ABERRATIONS.aberrationParameters.psfKernelDirectory,
					true);
			if (configPath==null){
				String msg="No PSF kernel directory selected, command aborted";
				System.out.println("Warning: "+msg);
				IJ.showMessage("Warning",msg);
				return;
			}
			configPath+=Prefs.getFileSeparator()+"config-combined";
			try {
				saveTimestampedProperties(
						configPath,      // full path or null
						null, // use as default directory if path==null 
						true,
						PROPERTIES);

			} catch (Exception e){
				String msg="Failed to save configuration to "+configPath+", command aborted";
				System.out.println("Error: "+msg);
				IJ.showMessage("Error",msg);
				return;
			}
			boolean combinedOK=	EYESIS_ABERRATIONS.combinePSFKernels ( // save configuration to combined kernels directory before calling this methode
					this.SYNC_COMMAND.stopRequested,
					INTERPOLATE,          // INTERPOLATE
					MULTIFILE_PSF ,         // MULTIFILE_PSF = new EyesisAberrations.MultiFilePSF(
					true, //              saveResult,
					false, //              showResult,
					UPDATE_STATUS,          // UPDATE_STATUS
					DEBUG_LEVEL,
					DEBUG_LEVEL);

			if (combinedOK) {
				///getPath()
				saveTimestampedProperties( // save config again
						configPath,      // full path or null
						null, // use as default directory if path==null 
						true,
						PROPERTIES);
			}
			return;
		}
		//		addButton("Combine Kernels",panelAberrations);

/* ======================================================================== */
		if       (label.equals("Interpolate Kernels")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			EYESIS_ABERRATIONS.setDistortions(LENS_DISTORTIONS);
			String configPath=EYESIS_ABERRATIONS.aberrationParameters.selectPSFKernelDirectory(
					true,
					EYESIS_ABERRATIONS.aberrationParameters.psfKernelDirectory,
					true);
			if (configPath==null){
				String msg="No PSF kernel directory selected, command aborted";
				System.out.println("Warning: "+msg);
				IJ.showMessage("Warning",msg);
				return;
			}
			configPath+=Prefs.getFileSeparator()+"config-interpolated";
			try {
				saveTimestampedProperties(
						configPath,      // full path or null
						null, // use as default directory if path==null 
						true,
						PROPERTIES);

			} catch (Exception e){
				String msg="Failed to save configuration to "+configPath+", command aborted";
				System.out.println("Error: "+msg);
				IJ.showMessage("Error",msg);
				return;
			}
			boolean interpolatedOK=	EYESIS_ABERRATIONS.interpolateKernels ( // save configuration to combined kernels directory before calling this methode
					this.SYNC_COMMAND.stopRequested,
					INTERPOLATE,          // INTERPOLATE
					MULTIFILE_PSF ,         // MULTIFILE_PSF = new EyesisAberrations.MultiFilePSF(
					true, //              saveResult,
					false, //              showResult,
					UPDATE_STATUS,          // UPDATE_STATUS
					DEBUG_LEVEL);

			if (interpolatedOK) {
				///getPath()
				saveTimestampedProperties( // save config again
						configPath,      // full path or null
						null, // use as default directory if path==null 
						true,
						PROPERTIES);
			}
			return;
		}
		//		addButton("Combine Kernels",panelAberrations);
//		"Invert Kernels"
/* ======================================================================== */
		if       (label.equals("Invert Kernels")) {
			DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
			if (LENS_DISTORTIONS==null) {
				IJ.showMessage("LENS_DISTORTION is not set");
				return;
			}
			EYESIS_ABERRATIONS.setDistortions(LENS_DISTORTIONS);
			String configPath=EYESIS_ABERRATIONS.aberrationParameters.selectAberrationsKernelDirectory(
					true,
					EYESIS_ABERRATIONS.aberrationParameters.aberrationsKernelDirectory,
					true);
			if (configPath==null){
				String msg="No aberrations (inverted) kernel directory selected, command aborted";
				System.out.println("Warning: "+msg);
				IJ.showMessage("Warning",msg);
				return;
			}
			configPath+=Prefs.getFileSeparator()+"config-inverted";
			try {
				saveTimestampedProperties(
						configPath,      // full path or null
						null, // use as default directory if path==null 
						true,
						PROPERTIES);

			} catch (Exception e){
				String msg="Failed to save configuration to "+configPath+", command aborted";
				System.out.println("Error: "+msg);
				IJ.showMessage("Error",msg);
				return;
			}
			boolean invertedOK=	EYESIS_ABERRATIONS.reverseKernels(
					this.SYNC_COMMAND.stopRequested,
		  			INVERSE,      // inverseParameters, // size (side of square) of direct PSF kernel
		   			true,         // saveResult,
		   			false,        // showResult,
		   			UPDATE_STATUS,// updateStatus,          // UPDATE_STATUS
		   			THREADS_MAX,  // threadsMax,
		   			DEBUG_LEVEL); // globalDebugLevel

			if (invertedOK) {
				///getPath()
				saveTimestampedProperties( // save config again
						configPath,      // full path or null
						null, // use as default directory if path==null 
						true,
						PROPERTIES);
			}
			return;
		}

/* ======================================================================== */

		IJ.showMessage("Not yet implemented");
		DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
		return;
	}
	
/* ===== Other methods ==================================================== */
	public boolean adjustFocusTiltLMA(){
		// just for reporting distance old way
/*		
		MOTORS.focusingHistory.optimalMotorPosition( // recalculate calibration to estimate current distance from center PSF
				FOCUS_MEASUREMENT_PARAMETERS,
    			MOTORS.getMicronsPerStep(), //double micronsPerStep,
    			DEBUG_LEVEL);
*/    			
		// No-move measure, add to history
		moveAndMaybeProbe(
				true, // just move, not probe
				null, // no move, just measure
				MOTORS,
				CAMERAS,
				LENS_DISTORTION_PARAMETERS,
				matchSimulatedPattern, // should not bee null - is null after grid center!!!
				FOCUS_MEASUREMENT_PARAMETERS,
				PATTERN_DETECT,
				DISTORTION,
				SIMUL,
				COMPONENTS,
				OTF_FILTER,
				PSF_PARS,
				THREADS_MAX,
				UPDATE_STATUS,
				MASTER_DEBUG_LEVEL,
				DISTORTION.loop_debug_level);
		//get measurement
		FocusingField.FocusingFieldMeasurement fFMeasurement=MOTORS.getThisFFMeasurement(FOCUSING_FIELD);	
		// calculate z, tx, ty, m1,m2,m3
	    double [] zTxTyM1M2M3 = FOCUSING_FIELD.adjustLMA(false,fFMeasurement,false); // allow tilt scan
		// show dialog: Apply, re-calculate, exit
    	int [] currentMotors=fFMeasurement.motors;
    	int [] newMotors=currentMotors.clone();
    	double [] zTxTy={Double.NaN,Double.NaN,Double.NaN};
    	if (zTxTyM1M2M3!=null){
    		newMotors[0]=(int) Math.round(zTxTyM1M2M3[3]);
    		newMotors[1]=(int) Math.round(zTxTyM1M2M3[4]);
    		newMotors[2]=(int) Math.round(zTxTyM1M2M3[5]);
    		zTxTy[0]=zTxTyM1M2M3[0];
    		zTxTy[1]=zTxTyM1M2M3[1];
    		zTxTy[2]=zTxTyM1M2M3[2];
    	}
    	double [] targetTilts={0.0,0.0};
    	double [] manualScrewsCW=null;
    	if (zTxTyM1M2M3!=null){
    		manualScrewsCW=FOCUSING_FIELD.fieldFitting.mechanicalFocusingModel.getManualScrews(
    				zTxTy[0]-FOCUSING_FIELD.targetRelFocalShift, //double zErr, // positive - away from lens
    				zTxTy[1]-targetTilts[0],                     // double tXErr,// positive - 1,2 away from lens, 3 - to the lens
    				zTxTy[2]-targetTilts[1]);                    // double tYErr);
    	}
    	double scaleMovement=1.0; // calculate automatically - reduce when close
    	boolean parallelMove=false;
    	if (MASTER_DEBUG_LEVEL>0){
    		System.out.println("----- Optimal (for qualB) focus/tilt -----");
    		System.out.println("Optimal absolute Zc="+FOCUSING_FIELD.qualBOptimizationResults[0]);
    		System.out.println("Optimal Tx="+FOCUSING_FIELD.qualBOptimizationResults[1]);
    		System.out.println("Optimal Ty="+FOCUSING_FIELD.qualBOptimizationResults[2]);

    		System.out.println("----- Focus/tilt measurement results -----");
    		System.out.println("Relative to optimal focal shift "+IJ.d2s(zTxTy[0],3)+" um ("+IJ.d2s(FOCUSING_FIELD.targetRelFocalShift,3)+"um)");
    		System.out.println("Relative to optimal horizontal tilt "+IJ.d2s(zTxTy[1],3)+" um/mm ("+IJ.d2s(targetTilts[0],3)+"um/mm)");
    		System.out.println("Relative to optimal vertical tilt "+IJ.d2s(zTxTy[2],3)+" um/mm ("+IJ.d2s(targetTilts[1],3)+"um/mm)");
    		for (int i=0;i<newMotors.length;i++){
        		System.out.println("Suggested for motor "+(i+1)+" "+newMotors[i]+" ("+currentMotors[i]+")");
    		}
    		if (manualScrewsCW!=null) for (int i=0;i<manualScrewsCW.length;i++){
    			double deg=360*Math.abs(manualScrewsCW[i]);
    			if (manualScrewsCW[i]>=0) System.out.println("Suggested rotation for screw # "+(i+1)+
    					" "+IJ.d2s(manualScrewsCW[i],3)+" ("+IJ.d2s(deg,0)+"\u00b0 CW)");
    			else  System.out.println("Suggested rotation for screw # "+(i+1)+
    					" "+IJ.d2s(manualScrewsCW[i],3)+" ("+IJ.d2s(deg,0)+"\u00b0 CCW)");
    		}
    		System.out.println("----- end of Focus/tilt measurement results -----");
    		
    		if (MASTER_DEBUG_LEVEL>0) System.out.println(FOCUSING_FIELD.showSamples());
    	}
    	GenericDialog gd = new GenericDialog("Adjusting focus/tilt");
    	if (zTxTyM1M2M3==null){
    		gd.addMessage("**** Failed to determine focus/tilt, probably too far out of focus. ****");
    		gd.addMessage("**** You may cancel the command and try \"Auto pre-focus\" first. ****");
    	}
        gd.addNumericField("Target focus (relative to best composite)",FOCUSING_FIELD.targetRelFocalShift,2,5,"um ("+IJ.d2s(zTxTy[0],3)+")");
        gd.addNumericField("Target horizontal tilt relative to optimal (normally 0)",targetTilts[0],2,5,"um/mm ("+IJ.d2s(zTxTy[1],3)+")");
        gd.addNumericField("Target vertical tilt  relative to optimal (normally 0)",targetTilts[1],2,5,"um/mm ("+IJ.d2s(zTxTy[2],3)+")");

        gd.addMessage("Optimal absolute Zc="+FOCUSING_FIELD.qualBOptimizationResults[0]);
        gd.addMessage("Optimal Tx="+FOCUSING_FIELD.qualBOptimizationResults[1]);
        gd.addMessage("Optimal Ty="+FOCUSING_FIELD.qualBOptimizationResults[2]);
        
        gd.addCheckbox("Optimize focal distance",(FOCUSING_FIELD.qualBOptimizeMode & 1) != 0);
		gd.addCheckbox("Optimize tiltX",         (FOCUSING_FIELD.qualBOptimizeMode & 2) != 0);
		gd.addCheckbox("Optimize tiltY",         (FOCUSING_FIELD.qualBOptimizeMode & 4) != 0);
        
        
		gd.addNumericField("Motor 1",newMotors[0],0,5,"steps ("+currentMotors[0]+")");
		gd.addNumericField("Motor 2",newMotors[1],0,5,"steps ("+currentMotors[1]+")");
		gd.addNumericField("Motor 3",newMotors[2],0,5,"steps ("+currentMotors[2]+")");
		gd.addMessage("Suggested rotation of the top screws, use if motor positions are out of limits - outside of +/-25,000");
		if (manualScrewsCW!=null)  for (int i=0;i<manualScrewsCW.length;i++){
			double deg=360*Math.abs(manualScrewsCW[i]);
			if (manualScrewsCW[i]>=0) gd.addMessage("Screw # "+(i+1)+" "+IJ.d2s(manualScrewsCW[i],3)+" ("+IJ.d2s(deg,0)+"\u00b0 CW)");
			else                      gd.addMessage("Screw # "+(i+1)+" "+IJ.d2s(manualScrewsCW[i],3)+" ("+IJ.d2s(deg,0)+"\u00b0 CCW)");
		}
		gd.addNumericField("Scale movement",scaleMovement,3,5,"x");
        gd.addCheckbox("Recalculate and apply parallel move only",parallelMove); // should be false after manual movement
		
        gd.addCheckbox("Filter samples/channels by Z",FOCUSING_FIELD.filterZ); // should be false after manual movement
        gd.addCheckbox("Filter by value (leave lower than maximal fwhm used in focal scan mode)",FOCUSING_FIELD.filterByScanValue);
        gd.addNumericField("Filter by value (remove samples above scaled best FWHM for channel/location)",FOCUSING_FIELD.filterByValueScale,2,5,"x");

        gd.addNumericField("Z min",FOCUSING_FIELD.zMin,2,5,"um");
        gd.addNumericField("Z max",FOCUSING_FIELD.zMax,2,5,"um");
        gd.addNumericField("Z step",FOCUSING_FIELD.zStep,2,5,"um");

        gd.addNumericField("Tilt min",FOCUSING_FIELD.tMin,2,5,"um/mm");
        gd.addNumericField("Tilt max",FOCUSING_FIELD.tMax,2,5,"um/mm");
        gd.addNumericField("Tilt step",FOCUSING_FIELD.tStep,2,5,"um/mm");
		
		gd.addNumericField("Motor anti-hysteresis travel (last measured was "+IJ.d2s(FOCUS_MEASUREMENT_PARAMETERS.measuredHysteresis,0)+")", FOCUS_MEASUREMENT_PARAMETERS.motorHysteresis, 0,7,"motors steps");
		gd.addNumericField("Debug Level:",                                MASTER_DEBUG_LEVEL, 0);
		gd.enableYesNoCancel("Apply movement","Re-measure"); // default OK (on enter) - "Apply"
    	WindowTools.addScrollBars(gd);
    	gd.showDialog();
		if (gd.wasCanceled()) return false;
		FOCUSING_FIELD.targetRelFocalShift=gd.getNextNumber();
		targetTilts[0]=                    gd.getNextNumber();
		targetTilts[1]=                    gd.getNextNumber();
		
		FOCUSING_FIELD.qualBOptimizeMode=0;
		FOCUSING_FIELD.qualBOptimizeMode+= gd.getNextBoolean()?1:0;
		FOCUSING_FIELD.qualBOptimizeMode+= gd.getNextBoolean()?2:0;
		FOCUSING_FIELD.qualBOptimizeMode+= gd.getNextBoolean()?4:0;
		
		newMotors[0]=               (int)  gd.getNextNumber();
		newMotors[1]=               (int)  gd.getNextNumber();
		newMotors[2]=               (int)  gd.getNextNumber();
		scaleMovement=                     gd.getNextNumber();
		parallelMove=                      gd.getNextBoolean();
        FOCUSING_FIELD.filterZ=            gd.getNextBoolean();
        FOCUSING_FIELD.filterByScanValue=  gd.getNextBoolean();
        FOCUSING_FIELD.filterByValueScale= gd.getNextNumber();

        FOCUSING_FIELD.zMin=               gd.getNextNumber();
        FOCUSING_FIELD.zMax=               gd.getNextNumber();
        FOCUSING_FIELD.zStep=              gd.getNextNumber();

        FOCUSING_FIELD.tMin=               gd.getNextNumber();
        FOCUSING_FIELD.tMax=               gd.getNextNumber();
        FOCUSING_FIELD.tStep=              gd.getNextNumber();
		FOCUS_MEASUREMENT_PARAMETERS.motorHysteresis= (int) gd.getNextNumber();
		MASTER_DEBUG_LEVEL=(         int)  gd.getNextNumber();
		DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
		FOCUSING_FIELD.setDebugLevel(DEBUG_LEVEL);
		if (parallelMove){ // ignore/recalculate newMotors data 
		    zTxTyM1M2M3 = FOCUSING_FIELD.adjustLMA(false,fFMeasurement,true); // recalculate with parallel move only, allow tilt scan
	    	newMotors=currentMotors.clone();
	    	if (zTxTyM1M2M3!=null){
	    		newMotors[0]=(int) Math.round(zTxTyM1M2M3[3]);
	    		newMotors[1]=(int) Math.round(zTxTyM1M2M3[4]);
	    		newMotors[2]=(int) Math.round(zTxTyM1M2M3[5]);
	    	}			
    		System.out.println("Parallel move position for motor 1 "+newMotors[0]+" ("+currentMotors[0]+")");
    		System.out.println("Parallel move position for motor 2 "+newMotors[1]+" ("+currentMotors[1]+")");
    		System.out.println("Parallel move position for motor 3 "+newMotors[2]+" ("+currentMotors[2]+")");
		}
		
//	Scale motor movement	
		newMotors[0]=currentMotors[0]+((int) Math.round((newMotors[0]-currentMotors[0])*scaleMovement));
		newMotors[1]=currentMotors[1]+((int) Math.round((newMotors[1]-currentMotors[1])*scaleMovement));
		newMotors[2]=currentMotors[2]+((int) Math.round((newMotors[2]-currentMotors[2])*scaleMovement));
		
		if (gd.wasOKed()){
			// Move, no measure
			MOTORS.moveElphel10364Motors(
					true, //boolean wait,
					newMotors,
					0.0, //double sleep,
					true, //boolean showStatus,
					"",   //String message,
					false); //!noHysteresis);			
		}
		return true;
	}
	
	
	
	public boolean checkSerialAndRestore(){
		// wait for camera
		CAMERAS.debugLevel=DEBUG_LEVEL;
		int frameNumber=CAMERAS.getCurrentFrameNumberWithTimeout(0, true,this.SYNC_COMMAND.stopRequested); // will throw on timeout
		if ((frameNumber<LAST_FRAME_NUMBER) || (frameNumber<0)){
			if (DEBUG_LEVEL>0) System.out.println("Camera was power cycled/rebooted (frameNumber="+frameNumber+", LAST_FRAME_NUMBER="+LAST_FRAME_NUMBER+")");
			MOTORS.setDebug(FOCUS_MEASUREMENT_PARAMETERS.motorDebug);
			MOTORS.resetInitialization();
			CAMERAS.resetInitialization();
			FOCUS_MEASUREMENT_PARAMETERS.cameraIsConfigured=false;
		} else {
			if (DEBUG_LEVEL>0) System.out.println("frameNumber ("+frameNumber+")>=LAST_FRAME_NUMBER ("+LAST_FRAME_NUMBER+")");
		}
		LAST_FRAME_NUMBER=frameNumber;
		if (FOCUS_MEASUREMENT_PARAMETERS.serialNumber==null) FOCUS_MEASUREMENT_PARAMETERS.serialNumber="";
		String currentSerial=CAMERAS.getSerialNumber(0,FOCUS_MEASUREMENT_PARAMETERS.EEPROM_channel);
		if ((currentSerial==null) || currentSerial.equals("")) {
			String msg="Failed to get SFE serial number";
			IJ.showStatus(msg);
			System.out.println (msg);
			if (MASTER_DEBUG_LEVEL>0) System.out.println(msg);
			FOCUS_MEASUREMENT_PARAMETERS.serialNumber="UNDEFINED"; //
			return false;
		}
		if (FOCUS_MEASUREMENT_PARAMETERS.serialNumber.equals(currentSerial)) return true; // serial did not change, OK to proceed
		GenericDialog gd=new GenericDialog("New SFE Serial Detected");
		gd.addMessage("Detected camera SFE serial ("+currentSerial+") differs from the one used ("+FOCUS_MEASUREMENT_PARAMETERS.serialNumber+")");
		gd.addMessage("Do you want to try to read previous data for the current SFE?");
		gd.addMessage("'Cancel' will keep the current lens S/N ("+FOCUS_MEASUREMENT_PARAMETERS.lensSerial+ ") and comments: "+FOCUS_MEASUREMENT_PARAMETERS.comment);
		gd.enableYesNoCancel("Yes, try to restore", "No, it is a new SFE");
		gd.showDialog();
		FOCUS_MEASUREMENT_PARAMETERS.serialNumber=currentSerial; // do that in any case, even if canceled
		
		if (gd.wasCanceled()) return false; // no change to comments, lens serial
		boolean restored=gd.wasOKed();
		FOCUS_MEASUREMENT_PARAMETERS.manufacturingState=0; // new SFE - reset for old format
		if (restored) restored=restoreSFELatest(); // may fail and return false
		if (restored) return true; // OK, restored
		// erase comments
		MOTORS.setDebug(FOCUS_MEASUREMENT_PARAMETERS.motorDebug);
		MOTORS.resetInitialization();
		CAMERAS.resetInitialization();
		FOCUS_MEASUREMENT_PARAMETERS.cameraIsConfigured=false;
		FOCUS_MEASUREMENT_PARAMETERS.serialNumber="";
		FOCUS_MEASUREMENT_PARAMETERS.comment="no comments"; // Comment to add to the results
		FOCUS_MEASUREMENT_PARAMETERS.lensSerial="????"; // Lens serial number
		FOCUS_MEASUREMENT_PARAMETERS.manufacturingState=0; // new SFE
		updateSerial(FOCUS_MEASUREMENT_PARAMETERS);
		return false;
	}
	
	public boolean restoreSFELatest(){
		MOTORS.setDebug(FOCUS_MEASUREMENT_PARAMETERS.motorDebug);
		MOTORS.resetInitialization();
		CAMERAS.resetInitialization();
		FOCUS_MEASUREMENT_PARAMETERS.serialNumber="";
		updateSerial(FOCUS_MEASUREMENT_PARAMETERS);
//		String dir=getResultsPath(FOCUS_MEASUREMENT_PARAMETERS);
		String dir=getResultsPath(FOCUS_MEASUREMENT_PARAMETERS,false); // always ask for directory
		File dFile=new File(dir);
		if (!dFile.isDirectory()) {
			String msg="Could not find directory with saved SFE configuration: "+dir;
			IJ.showMessage(msg);
			System.out.println("Error: "+msg);
			return false;
		}
		File[] fileList=dFile.listFiles(new Filter(".conf-xml"));
		if ((fileList==null) || (fileList.length==0)){
			String msg="Could not find any configuration files (*.conf-xml) in the directory "+dir;
			IJ.showMessage(msg);
			System.out.println("Error: "+msg);
			return false;
		}
		long lastTime=0;
		int index=0;
		for (int i=0;i<fileList.length;i++) if (fileList[i].lastModified()>lastTime) {
			lastTime=fileList[i].lastModified();
			index=i;
		} 
		String configPath=fileList[index].getAbsolutePath();
    	loadProperties(configPath, null, true, PROPERTIES);
    	return true;
	}
	
	
	public boolean autoLoadFiles(
			EyesisAberrations.AberrationParameters aberrationParameters,
			Distortions distortions, // should be initialized, after update DISTORTION_CALIBRATION_DATA from this
			Distortions.PatternParameters patternParameters,
    		Distortions.EyesisCameraParameters eyesisCameraParameters,
    		boolean updateStstus,
    		int debugLevel
			){
		if (distortions==null){
			return false;
		}
		distortions.debugLevel=debugLevel;
		String [] configPaths=aberrationParameters.autoLoadPaths();
		if (configPaths[0]==null) return false;
		System.out.println("+++++++++++ autoLoadFiles() eyesisCameraParameters.numStations="+eyesisCameraParameters.numStations+
				" +eyesisCameraParameters.goniometerHorizontal.length="+eyesisCameraParameters.goniometerHorizontal.length);
		Distortions.DistortionCalibrationData dcd=new Distortions.DistortionCalibrationData(
				true,
				configPaths[0],
				patternParameters,
				eyesisCameraParameters,
				null); // gridImages null - use specified files
		if (distortions.fittingStrategy!=null) {
			distortions.fittingStrategy.distortionCalibrationData=dcd;
		} else if (configPaths[1]==null) return false; // fitting strategy was null ind is not specified
		
		System.out.println("+++++++++++ autoLoadFiles() dcd.eyesisCameraParameters.numStations="+dcd.eyesisCameraParameters.numStations+
				" +dcd.eyesisCameraParameters.goniometerHorizontal.length="+dcd.eyesisCameraParameters.goniometerHorizontal.length);

		
		if (configPaths[1]!=null) {
			Distortions.FittingStrategy fs=distortions.fittingStrategy; // save old value
			distortions.fittingStrategy=new Distortions.FittingStrategy(
					true, // do not ask if specified
					configPaths[1],
					dcd);
			if (distortions.fittingStrategy.pathName== null){ // failed to select/open the file
				IJ.showMessage("No strategy selected");
				distortions.fittingStrategy=fs; // restore old strategy
				return false;
			}
			distortions.fittingStrategy.debugLevel=debugLevel;
			distortions.fittingStrategy.adjustNumberOfImages(dcd.gIP.length);
		}
		if (configPaths[2] !=null){ // load grid file
			patternParameters.debugLevel=debugLevel;
			patternParameters.updateStatus=updateStstus;
			if (DEBUG_LEVEL>0) System.out.println("Autolading grid file "+configPaths[2]);
			patternParameters.selectAndRestore(true,configPaths[2],dcd.eyesisCameraParameters.numStations); // returns path or null on failure
		}
		if (configPaths[3] !=null){ // load sensor
			if (distortions.fittingStrategy==null) return false;
			if (DEBUG_LEVEL>0) System.out.println("Autoloading sensor calibration files "+configPaths[3]);
			distortions.setDistortionFromImageStack(configPaths[3],aberrationParameters.autoRestoreSensorOverwriteOrientation);
		}
		return true;
	}
	

	
	
	
	
//	UV_LED_LASERS
	

	public void getAndSaveImage(
			boolean alwaysShow, // true overwrites focusMeasurementParameters.showResults
			boolean alwaysSave, // true overwrites focusMeasurementParameters.saveResults
			CalibrationHardwareInterface.CamerasInterface camerasInterface,
			Distortions.LensDistortionParameters lensDistortionParameters,
			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
			boolean updateStatus,
			int debugLevel
	){
		long 	  startTime=System.nanoTime();
		if (!focusMeasurementParameters.cameraIsConfigured) {
			if (camerasInterface.showDialog("Configure cameras interface", 1, true)){
				focusMeasurementParameters.cameraIsConfigured=true;
//				IJ.showMessage("Notice","Make sure camera is in TRIG=4 mode, JP4, correct exposure/white balance...");
			} else {
				IJ.showMessage("Error","Camera is not configured\nProcess canceled");
				return;
			}
		}
		// acquire camera image here, no lasers.
		focusMeasurementParameters.sensorTemperature=camerasInterface.getSensorTemperature(0,FOCUS_MEASUREMENT_PARAMETERS.EEPROM_channel);
		ImagePlus imp= camerasInterface.acquireSingleImage (
				true, //boolean useLasers,
				updateStatus);
		if (imp==null){
			IJ.showMessage("Error","Failed to get camera image\nProcess canceled");
			return;
		}
		// set motors, timestamp, ...
		imp.show();
		imp.updateAndDraw();
		if (debugLevel>0) System.out.println("Image acquisition done at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
// all are read as "null"
		imp.setProperty("comment",focusMeasurementParameters.comment);
		if (!Double.isNaN(focusMeasurementParameters.sensorTemperature))
			imp.setProperty("sensorTemperature", ""+focusMeasurementParameters.sensorTemperature);
		imp.setProperty("px0", ""+lensDistortionParameters.px0);
		imp.setProperty("py0", ""+lensDistortionParameters.py0);
		imp.setProperty("motor1", ""+focusMeasurementParameters.motorPos[0]);
		imp.setProperty("motor2", ""+focusMeasurementParameters.motorPos[1]);
		imp.setProperty("motor3", ""+focusMeasurementParameters.motorPos[2]);
    	(new JP46_Reader_camera(false)).encodeProperiesToInfo(imp);
		if (alwaysShow || focusMeasurementParameters.showResults) imp.show();
		if (alwaysSave || focusMeasurementParameters.saveResults) {
			String dir=getResultsPath(focusMeasurementParameters);
			File dFile=new File(dir);
			if (!dFile.isDirectory() &&  !dFile.mkdirs()) {
				String msg="Failed to create directory "+dir;
				IJ.showMessage(msg);
				throw new IllegalArgumentException (msg);
			}
			String lensPrefix="";
			if (FOCUS_MEASUREMENT_PARAMETERS.includeLensSerial && (FOCUS_MEASUREMENT_PARAMETERS.lensSerial.length()>0)){
//				lensPrefix=String.format("LENS%S-",FOCUS_MEASUREMENT_PARAMETERS.lensSerial);
				lensPrefix=String.format("LENS%S-S%02d-",FOCUS_MEASUREMENT_PARAMETERS.lensSerial,FOCUS_MEASUREMENT_PARAMETERS.manufacturingState);
			}
			String path=dFile+Prefs.getFileSeparator()+lensPrefix+
			((String)imp.getProperty("timestamp")).replace('.','_')+".tiff";
			if (debugLevel>0) System.out.println ("Saving annotated registered image to "+path);
			if (updateStatus) IJ.showStatus("Saving annotated registered image to "+path);
			(new FileSaver(imp)).saveAsTiff(path);
		}
	}
	
	
	public boolean saveCurrentConfig(){
		if (!FOCUS_MEASUREMENT_PARAMETERS.saveResults) return false;
		String dir=getResultsPath(FOCUS_MEASUREMENT_PARAMETERS);
		File dFile=new File(dir);
		if (!dFile.isDirectory() &&  !dFile.mkdirs()) {
			String msg="Failed to create directory "+dir;
			IJ.showMessage(msg);
			throw new IllegalArgumentException (msg);
		}
		String lensPrefix="";
		if (FOCUS_MEASUREMENT_PARAMETERS.includeLensSerial && (FOCUS_MEASUREMENT_PARAMETERS.lensSerial.length()>0)){
//			lensPrefix=String.format("LENS%S-",FOCUS_MEASUREMENT_PARAMETERS.lensSerial);
			lensPrefix=String.format("LENS%S-S%02d-",FOCUS_MEASUREMENT_PARAMETERS.lensSerial,FOCUS_MEASUREMENT_PARAMETERS.manufacturingState);

		}
		String path=dFile+Prefs.getFileSeparator()+lensPrefix+CAMERAS.getLastTimestampUnderscored(); // +".conf-xml" will be added in saveProperties
		if (MASTER_DEBUG_LEVEL>0) System.out.println ("Saving current configuration parameters to "+path);
		saveProperties(
				path,
				null,
				true, // use xml
				PROPERTIES);
		return true;
	}
	public String getResultsPath(
			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters){
		return getResultsPath(focusMeasurementParameters, true);
	}
	
	public String getResultsPath(
			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
			boolean smart){
		updateSerial(focusMeasurementParameters);
		String dir= CalibrationFileManagement.selectDirectory(
				smart, //true, // smart,
				true, // newAllowed, // save  
				"Superdirectory to save/restore SFE focusing results", // title
				"Select results superdirectory (having SFE subdirs)", // button
				null, // filter
				focusMeasurementParameters.resultsSuperDirectory);
		if (dir!=null) focusMeasurementParameters.resultsSuperDirectory=dir;
		dir+=Prefs.getFileSeparator()+FOCUS_MEASUREMENT_PARAMETERS.serialNumber;
		return dir;
	}
	
	public void updateSerial(
			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters){
		if ((FOCUS_MEASUREMENT_PARAMETERS.serialNumber==null)|| (FOCUS_MEASUREMENT_PARAMETERS.serialNumber=="")) {
			FOCUS_MEASUREMENT_PARAMETERS.serialNumber=CAMERAS.getSerialNumber(0,FOCUS_MEASUREMENT_PARAMETERS.EEPROM_channel);
			if ((FOCUS_MEASUREMENT_PARAMETERS.serialNumber==null)|| (FOCUS_MEASUREMENT_PARAMETERS.serialNumber=="")) {
				String msg="Failed to get SFE serial number";
				IJ.showStatus(msg);
				if (MASTER_DEBUG_LEVEL>0) System.out.println(msg);
				FOCUS_MEASUREMENT_PARAMETERS.serialNumber="UNDEFINED"; // 
			} else {
				String msg="SFE serial number="+FOCUS_MEASUREMENT_PARAMETERS.serialNumber;
				IJ.showStatus(msg);
				if (MASTER_DEBUG_LEVEL>0) System.out.println(msg);
			}
		}

	}
	
	public int [] focusingCenterStepsAuto(
			int numIterations, // maximal number of iterations (0 - suggest only, do not move). When calling from the button - first time single iteration, second time - as specified
			double focusTolerance, // will exit after whatever comes first tolearance or number of iterations
			CalibrationHardwareInterface.FocusingMotors focusingMotors,
			CalibrationHardwareInterface.CamerasInterface camerasInterface,
			Distortions.LensDistortionParameters lensDistortionParameters,
			MatchSimulatedPattern matchSimulatedPattern, // should not bee null
			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			MatchSimulatedPattern.DistortionParameters distortionParameters,
			SimulationPattern.SimulParameters  simulParameters,
			EyesisAberrations.ColorComponents colorComponents,
			EyesisAberrations.OTFFilterParameters otfFilterParameters,
			EyesisAberrations.PSFParameters psfParameters,
			int       threadsMax,
			boolean   updateStatus,
			int       debugLevel,
			int       loopDebugLevel){
		int debugThreshold=0;
		if (focusingMotors.historySize()==0){
			if (debugLevel>debugThreshold) {
				System.out.println("focusingCenterStepsAuto(): history is empty, measuring the first sample");
			}
			moveMeasureAndSave( // measure "here" and put to history
					false,
					null, // first time - will repeat measurement at the current position
					focusingMotors,
					camerasInterface,
					lensDistortionParameters,
					matchSimulatedPattern, // should not bee null
					focusMeasurementParameters,
					patternDetectParameters,
					distortionParameters,
					simulParameters,
					colorComponents,
					otfFilterParameters,
					psfParameters,
					threadsMax,
					updateStatus,
					debugLevel,
					loopDebugLevel);
		
		}
//		for (int nIter=0;nIter<numIterations;nIter++){
		int nIter=0;
		int [] newPos=null;
		while (true) {
			newPos=focusingMotors.focusingHistory.focusReadjustStep(
					focusMeasurementParameters.targetMicrons, //double targetMicrons, // target focal distance
					focusMeasurementParameters.goodDistanceSigma, // double micronsFade, // reduce influence of far points
					focusMeasurementParameters.lensDistanceWeightK, // double lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
					focusMeasurementParameters.lensDistanceWeightY, // double lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
					debugLevel+1); //int debugLevel
			if (newPos==null) {
				System.out.println("focusingCenterStepsAuto() failed");
				return null;
			}
			if (debugLevel>debugThreshold) System.out.println("focusingCenterStepsAuto(): step #"+(nIter+1)+" - moving to ["+newPos[0]+","+newPos[1]+","+newPos[2]+"]");
			if (numIterations<=0) break;

			moveMeasureAndSave( // measure "here" and put to history
					false,
					newPos, // first time - will repeat measurement at the current position
					focusingMotors,
					camerasInterface,
					lensDistortionParameters,
					matchSimulatedPattern, // should not bee null
					focusMeasurementParameters,
					patternDetectParameters,
					distortionParameters,
					simulParameters,
					colorComponents,
					otfFilterParameters,
					psfParameters,
					threadsMax,
					updateStatus,
					debugLevel,
					loopDebugLevel);
			// measure focal distance and comparte to threshold
			double fDist=    	focusingMotors.focusingHistory.getLensDistance(
	    			true, // return absolutely calibrated data
	    			focusMeasurementParameters.lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
	    			focusMeasurementParameters.lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
	    			debugLevel
	    			);
			if (Math.abs(fDist-focusMeasurementParameters.targetMicrons)<focusTolerance){
				if (debugLevel>debugThreshold) System.out.println("focusingCenterStepsAuto(): matched tolerance, distance="+fDist +" um ( target= "+
						focusMeasurementParameters.targetMicrons+" um, tolerance= +/-"+focusTolerance+" um)");
				break;
			}
			if (++nIter>=numIterations) {
				if (debugLevel>debugThreshold) System.out.println("focusingCenterStepsAuto(): exhausted iterations, distance="+fDist +" um ( target= "+
						focusMeasurementParameters.targetMicrons+" um, tolerance= +/-"+focusTolerance+" um)");
				break;
			}
			
		}
		return newPos;
		
	}
	public int [] focusingStepsAuto(
			int numIterations, // maximal number of iterations (0 - suggest only, do not move). When calling from the button - first time single iteration, second time - as specified
			CalibrationHardwareInterface.FocusingMotors focusingMotors,
			CalibrationHardwareInterface.CamerasInterface camerasInterface,
			Distortions.LensDistortionParameters lensDistortionParameters,
			MatchSimulatedPattern matchSimulatedPattern, // should not bee null
			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			MatchSimulatedPattern.DistortionParameters distortionParameters,
			SimulationPattern.SimulParameters  simulParameters,
			EyesisAberrations.ColorComponents colorComponents,
			EyesisAberrations.OTFFilterParameters otfFilterParameters,
			EyesisAberrations.PSFParameters psfParameters,
			int       threadsMax,
			boolean   updateStatus,
			int       debugLevel,
			int       loopDebugLevel){
		long 	  startTime=System.nanoTime();
        int [] initialMotorPos=focusingMotors.readElphel10364Motors().clone();
		int [] position=null;
		// just to move before probing around, otherwise may be removed (will do it anyway if center is far from the goal)
		if ((focusingMotors.historySize()>0) &&focusMeasurementParameters.lensDistanceMoveToGoal && !focusMeasurementParameters.useRadialTangential) {
			// try to move center first (if calibrated)
			// skip if history is empty (will start with probing around)
			
			position=focusingMotors.focusingHistory.focusReadjustStep(
					focusMeasurementParameters.targetMicrons, //double targetMicrons, // target focal distance
					focusMeasurementParameters.goodDistanceSigma, // double micronsFade, // reduce influence of far points
					focusMeasurementParameters.lensDistanceWeightK, // double lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
					focusMeasurementParameters.lensDistanceWeightY, // double lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
					debugLevel//int debugLevel
	     			);			
			if (position==null) {
				String msg="Failure: Calibration does not exist or requested position is out of the calibrated range. Did you run \"Scan Calib\"?";
	        	if (debugLevel>0) System.out.println(msg);
	        	IJ.showMessage(msg);
	        	return null; // need calibration first
			}
		} else {
			position=initialMotorPos;
		}

		moveAndMaybeProbe(
				false,
				position, //null, // first time - will repeat measurement at the current position
				focusingMotors,
				camerasInterface,
				lensDistortionParameters,
				matchSimulatedPattern, // should not be null
				focusMeasurementParameters,
				patternDetectParameters,
				distortionParameters,
				simulParameters,
				colorComponents,
				otfFilterParameters,
				psfParameters,
				threadsMax,
				updateStatus,
				debugLevel,
				loopDebugLevel);
		int leftFinal=focusMeasurementParameters.numFinalCorr;
		int [] result=null;
		int iterations=0;
		double micronsFade=(focusMeasurementParameters.filterGoodDistance)?focusMeasurementParameters.goodDistanceSigma:0.0;
		double tiltFade=(focusMeasurementParameters.filterGoodDistance)?focusMeasurementParameters.goodTiltSigma:0.0;
			
        for (int iterNumber=0;(iterNumber==0) || (iterNumber<numIterations);iterNumber++) {// normally will break out earlier, it is just to prevent endless loops
        	if (debugLevel>0) System.out.println("focusingStepsAuto(): Tilt/Focus iteration="+(iterNumber+1)+" of maximum "+numIterations);
    		int [] lastMotorPos=focusingMotors.focusingHistory.getPosition();
    		// check focal distance error, if too big - readjust parallel
    		double lastDist=focusingMotors.focusingHistory.getLensDistance(
    				focusingMotors.focusingHistory.getCenterResolutions(), // {R-sharpness,G-sharpness,B-sharpness}
	    			true, // boolean absolute, // return absolutely calibrated data
	    			focusMeasurementParameters.lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
	    			focusMeasurementParameters.lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
	    			1 //debugLevel
	    			);
    		int [] newMotorPos=null;
    		// focal distance in the center too far, readjust before correcting tilts
    		if ((focusMeasurementParameters.parallelAdjustThreshold>0) && (Math.abs(lastDist-focusMeasurementParameters.targetMicrons)>focusMeasurementParameters.parallelAdjustThreshold)){
    			if (debugLevel>0) System.out.println("Focal distance for the center  ("+IJ.d2s(lastDist,2)+"um is too far from the requested "+
    					focusMeasurementParameters.targetMicrons+"um ( threshold="+micronsFade+"um), readjusting distance by moving 3 motors together instead of full adjustment");
    			
    			newMotorPos=focusingMotors.focusingHistory.focusReadjustStep(
    					focusMeasurementParameters.targetMicrons, //double targetMicrons, // target focal distance
    					focusMeasurementParameters.goodDistanceSigma, // double micronsFade, // reduce influence of far points
    					focusMeasurementParameters.lensDistanceWeightK, // double lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
    					focusMeasurementParameters.lensDistanceWeightY, // double lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
    					debugLevel//int debugLevel
    	     			);			
    			if (position==null) {
    				String msg="Failure: Calibration does not exist or requested position is out of the calibrated range. Did you run \"Scan Calib\"?";
    	        	if (debugLevel>0) System.out.println(msg);
    	        	IJ.showMessage(msg);
    	        	return null; // need calibration first
    			}
    			
    		} else {
    			newMotorPos=focusingMotors.focusingHistory.focusTiltStep( 
    					focusMeasurementParameters.useTheBest, // start from the best sample, (false - from the last)
    					focusMeasurementParameters.targetMicrons, //double targetMicrons, // target focal distance
    					micronsFade,
    					tiltFade,
    					focusMeasurementParameters.probe_M1M2M3*FOCUS_MEASUREMENT_PARAMETERS.sigmaToProbe, // sigma for decay for average of the 3 motors: (M1+M2+M3)/3
    					focusMeasurementParameters.probe_M3_M1M2*FOCUS_MEASUREMENT_PARAMETERS.sigmaToProbe, // sigma for decay for M3 opposite to M1 and M2: M3-(M1+M2)/2
    					focusMeasurementParameters.probe_M2_M1*FOCUS_MEASUREMENT_PARAMETERS.sigmaToProbe,// sigma for decay for  M2 opposite to M1:    M2-M1
    					focusMeasurementParameters.lensDistanceWeightK, // double lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
    					focusMeasurementParameters.lensDistanceWeightY, // double lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
    					focusMeasurementParameters.weightRatioRedToGreen,  // used for tilt averaging
    					focusMeasurementParameters.weightRatioBlueToGreen, // used for tilt averaging
    					focusMeasurementParameters.toleranceMicrons,
    					focusMeasurementParameters.toleranceTilt,
    					focusMeasurementParameters.toleranceThreshold,
    					focusMeasurementParameters.maxStep,        // Maximal motor move for each step
    					debugLevel); //+1); //int debugLevel
    		}
// if no move, that means that correction was under threshold and was zeroed-out - exit with success. That would happen after
// previous movement (i.e. readjustment)			
			if ((newMotorPos[0]==lastMotorPos[0]) && (newMotorPos[1]==lastMotorPos[1]) &&	(newMotorPos[2]==lastMotorPos[2])){
	        	if (debugLevel>0) {
	        		System.out.println("Suggested position did not change - probably residuals are under tolerances, m1="+lastMotorPos[0]+"m2="+lastMotorPos[1]+"m3="+lastMotorPos[2]);
	        	}
				result=newMotorPos; 
				break; 
			}
			if (numIterations==0){ // do not move, just suggest 
				result=newMotorPos;
				break; 
			}
			double overallMoveDistance=Math.sqrt(
					(newMotorPos[0]-initialMotorPos[0])*(newMotorPos[0]-initialMotorPos[0])+
					(newMotorPos[1]-initialMotorPos[1])*(newMotorPos[1]-initialMotorPos[1])+
					(newMotorPos[2]-initialMotorPos[2])*(newMotorPos[2]-initialMotorPos[2]));
			if (overallMoveDistance>focusMeasurementParameters.maxAutoDistance){
				String msg="Exceeded overall travel for adjustment, it is calculated "+IJ.d2s(overallMoveDistance,0)+" steps, limit is set to "+focusMeasurementParameters.maxAutoDistance+
				" steps.\nYou may verify the motors are safe to continue and just re-run adjustment from this point, or probably change some shims first.";
				System.out.println(msg);
				IJ.showMessage(msg);
				break; // failure - too far
			}
			
			boolean probeIsNeeded=!(focusingMotors.distFromProbed()<focusMeasurementParameters.reProbeDistance); // distFromProbed() may return NaN, and will - first time
				moveAndMaybeProbe(
						!probeIsNeeded,
						newMotorPos, // will first measure around destination, then go to the destination and measure there
						focusingMotors,
						camerasInterface,
						lensDistortionParameters,
						matchSimulatedPattern, // should not bee null
						focusMeasurementParameters,
						patternDetectParameters,
						distortionParameters,
						simulParameters,
						colorComponents,
						otfFilterParameters,
						psfParameters,
						threadsMax,
						updateStatus,
						debugLevel,
						loopDebugLevel);
			iterations++;
			// See if it is time to exit
			// is it a success?
			double lastMoveDistance=Math.sqrt(
					(newMotorPos[0]-lastMotorPos[0])*(newMotorPos[0]-lastMotorPos[0])+
					(newMotorPos[1]-lastMotorPos[1])*(newMotorPos[1]-lastMotorPos[1])+
					(newMotorPos[2]-lastMotorPos[2])*(newMotorPos[2]-lastMotorPos[2]));
			if 	(lastMoveDistance<focusMeasurementParameters.minCorr){
				if (--leftFinal<=0){
					result=newMotorPos;
					break;
				}
			} else leftFinal=focusMeasurementParameters.numFinalCorr;
			if (numIterations==1)result=newMotorPos; // if the requested was just a single iteration and no other failures - it is OK (returns true)
		}
        if (iterations==0){
        	if (debugLevel>0) System.out.println("Automatic focus adjustment "+((result!=null)?"succeeded":"failed")+
    				" to suggest new position at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
        	
        } else {
        	if (debugLevel>0) System.out.println("Automatic focus adjustment "+((result!=null)?"succeeded":"failed")+
    				" after "+iterations+" iterations at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
        }
        
		return result;
	}
	
	public int [] focusingStepsAutoOld(
			int numIterations, // maximal number of iterations (0 - suggest only, do not move). When calling from the button - first time single iteration, second time - as specified
			CalibrationHardwareInterface.FocusingMotors focusingMotors,
			CalibrationHardwareInterface.CamerasInterface camerasInterface,
			Distortions.LensDistortionParameters lensDistortionParameters,
			MatchSimulatedPattern matchSimulatedPattern, // should not bee null
			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			MatchSimulatedPattern.DistortionParameters distortionParameters,
			SimulationPattern.SimulParameters  simulParameters,
			EyesisAberrations.ColorComponents colorComponents,
			EyesisAberrations.OTFFilterParameters otfFilterParameters,
			EyesisAberrations.PSFParameters psfParameters,
			int       threadsMax,
			boolean   updateStatus,
			int       debugLevel,
			int       loopDebugLevel){
		long 	  startTime=System.nanoTime();
        int [] initialMotorPos=focusingMotors.readElphel10364Motors().clone();
		int [] position=null;
		if (focusMeasurementParameters.lensDistanceMoveToGoal && !focusMeasurementParameters.useRadialTangential) {
			// try to move center first (if calibrated)
			position=focusingMotors.focusingHistory.parallelPositionMove(focusMeasurementParameters.targetMicrons,
        			focusMeasurementParameters.lensDistanceWeightY, 
        			focusMeasurementParameters.lensDistanceWeightK,
        			debugLevel);
			if (position==null) {
				String msg="Failure: Calibration does not exist or requested position is out of the calibrated range. Did you run \"Scan Calib\"?";
	        	if (debugLevel>0) System.out.println(msg);
	        	IJ.showMessage(msg);
	        	return null; // need calibration first
			}
		} else {
			position=initialMotorPos;
		}
//		if (!(focusingMotors.distFromProbed(position)<focusMeasurementParameters.reProbeDistance)){ // distFromProbed() may return NaN, and will - first time
			moveAndMaybeProbe(
					false,
					position, //null, // first time - will repeat measurement at the current position
					focusingMotors,
					camerasInterface,
					lensDistortionParameters,
					matchSimulatedPattern, // should not bee null
					focusMeasurementParameters,
					patternDetectParameters,
					distortionParameters,
					simulParameters,
					colorComponents,
					otfFilterParameters,
					psfParameters,
					threadsMax,
					updateStatus,
					debugLevel,
					loopDebugLevel);
//		}
		int leftFinal=focusMeasurementParameters.numFinalCorr;
		int [] result=null;
		int iterations=0;

        for (int iterNumber=0;(iterNumber==0) || (iterNumber<numIterations);iterNumber++) {// normally will break out earlier, it is just to prevent endless loops
        	if (debugLevel>0) System.out.println("focusingStepsAuto(): Tilt/Focus iteration="+(iterNumber+1)+" of maximum "+numIterations);
			//		if ((focusingMotors.historySize()<4){
//            int [] lastMotorPos=focusingMotors.readElphel10364Motors().clone();
    		int [] lastMotorPos=focusingMotors.focusingHistory.getPosition();
			int [] newMotorPos=focusingMotors.focusingHistory.solveFocusing(
					focusMeasurementParameters.weightRatioRedToGreen,
   	    			focusMeasurementParameters.weightRatioBlueToGreen,
	    			focusMeasurementParameters.lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
	    			focusMeasurementParameters.lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
					focusMeasurementParameters.useRadialTangential,
					focusMeasurementParameters.targetFarNear,
					focusMeasurementParameters.targetMicrons,
					focusMeasurementParameters.toleranceMicrons,
					focusMeasurementParameters.toleranceTilt,
					focusMeasurementParameters.toleranceThreshold,
					focusMeasurementParameters.motorsSigma,
					focusMeasurementParameters.believeLast,   // 0 - use 'honest' best fit, 1.0 - make each plane go through the last sample 
					focusMeasurementParameters.maxStep,        // Maximal motor move for each step
					(focusMeasurementParameters.filterGoodDistance)?focusMeasurementParameters.goodDistanceSigma:0.0, // to use only samples with small distance errors
					debugLevel+1
			);
// if no move, that means that correction was under threshold and was zeroed-out - exit with success. That would happen after
// previous movent (i.e. readjustment)			
			if ((newMotorPos[0]==lastMotorPos[0]) && (newMotorPos[1]==lastMotorPos[1]) &&	(newMotorPos[2]==lastMotorPos[2])){
	        	if (debugLevel>0) {
	        		System.out.println("Suggested position did not change - probably residuals are under tolerances, m1="+lastMotorPos[0]+"m2="+lastMotorPos[1]+"m3="+lastMotorPos[2]);
	        	}
				result=newMotorPos; 
				break; 
			}
			if (numIterations==0){ // do not move, just suggest 
				result=newMotorPos;
				break; 
			}
			double overallMoveDistance=Math.sqrt(
					(newMotorPos[0]-initialMotorPos[0])*(newMotorPos[0]-initialMotorPos[0])+
					(newMotorPos[1]-initialMotorPos[1])*(newMotorPos[1]-initialMotorPos[1])+
					(newMotorPos[2]-initialMotorPos[2])*(newMotorPos[2]-initialMotorPos[2]));
			if (overallMoveDistance>focusMeasurementParameters.maxAutoDistance){
				break; // failure - too far
			}
			if (focusMeasurementParameters.parallelAdjustThreshold>0){ // 
				//Re-adjust focus by rotating 3 motors in parallel - focus is much more sensitive than tilt
	        	if (debugLevel>0) System.out.println("focusingStepsAuto() - before readjusting focus, move-and-probe");
				moveAndMaybeProbe(
						true, // just move
						newMotorPos, // will first measure around destination, then go to the destination and measure there
						focusingMotors,
						camerasInterface,
						lensDistortionParameters,
						matchSimulatedPattern, // should not bee null
						focusMeasurementParameters,
						patternDetectParameters,
						distortionParameters,
						simulParameters,
						colorComponents,
						otfFilterParameters,
						psfParameters,
						threadsMax,
						updateStatus,
						debugLevel,
						loopDebugLevel);
				double worstOver=focusingMotors.focusingHistory.worstOverTolerance (
						focusMeasurementParameters.weightRatioRedToGreen,
						focusMeasurementParameters.weightRatioBlueToGreen,
						focusMeasurementParameters.lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
						focusMeasurementParameters.lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
						focusMeasurementParameters.useRadialTangential,
						focusMeasurementParameters.targetFarNear,
						focusMeasurementParameters.targetMicrons,
						focusMeasurementParameters.toleranceMicrons,
						focusMeasurementParameters.toleranceTilt,
		    			debugLevel);
				if (worstOver<1.0){

					result=focusingMotors.focusingHistory.getPosition();
		        	if (debugLevel>0) {
		        		System.out.println("-- worstOver= "+IJ.d2s(worstOver,4)+", m1="+result[0]+"m2="+result[1]+"m3="+result[2]);
		        	}
					break;
				}
				int [] focusReAdjusted=focusingMotors.focusingHistory.parallelPositionMove(focusMeasurementParameters.targetMicrons,
	        			focusMeasurementParameters.lensDistanceWeightY, 
	        			focusMeasurementParameters.lensDistanceWeightK,
	        			debugLevel);
				if (focusReAdjusted!=null){
		        	if (debugLevel>0) System.out.println("focusingStepsAuto() - readjusting focus: Tilt/Focus iteration="+(iterNumber+1)+" of maximum "+numIterations);
					newMotorPos=focusReAdjusted;
					overallMoveDistance=Math.sqrt(
							(newMotorPos[0]-initialMotorPos[0])*(newMotorPos[0]-initialMotorPos[0])+
							(newMotorPos[1]-initialMotorPos[1])*(newMotorPos[1]-initialMotorPos[1])+
							(newMotorPos[2]-initialMotorPos[2])*(newMotorPos[2]-initialMotorPos[2]));
					if (overallMoveDistance>focusMeasurementParameters.maxAutoDistance){
						break; // failure - too far
					}
				}
			}
			boolean probeIsNeeded=!(focusingMotors.distFromProbed()<focusMeasurementParameters.reProbeDistance); // distFromProbed() may return NaN, and will - first time
				moveAndMaybeProbe(
						!probeIsNeeded,
						newMotorPos, // will first measure around destination, then go to the destination and measure there
						focusingMotors,
						camerasInterface,
						lensDistortionParameters,
						matchSimulatedPattern, // should not bee null
						focusMeasurementParameters,
						patternDetectParameters,
						distortionParameters,
						simulParameters,
						colorComponents,
						otfFilterParameters,
						psfParameters,
						threadsMax,
						updateStatus,
						debugLevel,
						loopDebugLevel);
			iterations++;
			// See if it is time to exit
			// is it a success?
			double lastMoveDistance=Math.sqrt(
					(newMotorPos[0]-lastMotorPos[0])*(newMotorPos[0]-lastMotorPos[0])+
					(newMotorPos[1]-lastMotorPos[1])*(newMotorPos[1]-lastMotorPos[1])+
					(newMotorPos[2]-lastMotorPos[2])*(newMotorPos[2]-lastMotorPos[2]));
			if 	(lastMoveDistance<focusMeasurementParameters.minCorr){
				if (--leftFinal<=0){
					result=newMotorPos;
					break;
				}
			} else leftFinal=focusMeasurementParameters.numFinalCorr;
			if (numIterations==1)result=newMotorPos; // if the requested was just a single iteration and no other failures - it is OK (returns true)
		}
        if (iterations==0){
        	if (debugLevel>0) System.out.println("Automatic focus adjustment "+((result!=null)?"succeeded":"failed")+
    				" to suggest new position at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
        	
        } else {
        	if (debugLevel>0) System.out.println("Automatic focus adjustment "+((result!=null)?"succeeded":"failed")+
    				" after "+iterations+" iterations at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
        }
        
		return result;
	}

	
	public int [] fineFocusingStepsAuto(
			int numIterations, // maximal number of iterations (0 - suggest only, do not move). When calling from the button - first time single iteration, second time - as specified
			CalibrationHardwareInterface.FocusingMotors focusingMotors,
			CalibrationHardwareInterface.CamerasInterface camerasInterface,
			Distortions.LensDistortionParameters lensDistortionParameters,
			MatchSimulatedPattern matchSimulatedPattern, // should not bee null
			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			MatchSimulatedPattern.DistortionParameters distortionParameters,
			SimulationPattern.SimulParameters  simulParameters,
			EyesisAberrations.ColorComponents colorComponents,
			EyesisAberrations.OTFFilterParameters otfFilterParameters,
			EyesisAberrations.PSFParameters psfParameters,
			int       threadsMax,
			boolean   updateStatus,
			int       debugLevel,
			int       loopDebugLevel){
		long 	  startTime=System.nanoTime();
		if (focusingMotors.sizeFocus()==0){
			moveAndMaybeProbe(
					true, // just move, no "probing" in this focusing mode
					null, // first time - will repeat measurement at the current position
					focusingMotors,
					camerasInterface,
					lensDistortionParameters,
					matchSimulatedPattern, // should not bee null
					focusMeasurementParameters,
					patternDetectParameters,
					distortionParameters,
					simulParameters,
					colorComponents,
					otfFilterParameters,
					psfParameters,
					threadsMax,
					updateStatus,
					debugLevel,
					loopDebugLevel);
		}

		int leftFinal=focusMeasurementParameters.numFinalCorr;
		int [] result=null;
		int iterations=0;
        for (int iterNumber=0;(iterNumber==0) || (iterNumber<numIterations);iterNumber++) {// normally will break out earlier, it is just to prevent endless loops
        	if (debugLevel>0) System.out.println("fineFocusingStepsAuto(): iteration="+(iterNumber+1)+" of maximum "+numIterations);
			//		if ((focusingMotors.historySize()<4){
            int [] lastMotorPos=focusingMotors.readElphel10364Motors().clone();
            int [] newMotorPos=focusingMotors.solveSinglePoly(
            		(iterNumber==0), // first call
            		focusMeasurementParameters.weightRatioRedToGreen,
            		1.0,
            		focusMeasurementParameters.weightRatioBlueToGreen,
            		true, //boolean sameTiltOnly, // us only history where the difference between motors was the same as in the last sample
            		true, //boolean centerOnly, // only use center samples
            		2, //int polyDegree,      // polynomial order (just 2 always)

            		focusMeasurementParameters.motorsSigma3,    // another sigma (when all 3 together, smaller than when tilt is involved!)
            		focusMeasurementParameters.maxAutoDistance, // Maximal motor move for each step
            		focusMeasurementParameters.maxLinearStep,
            		focusMeasurementParameters.motorsMinSigma, 
            		focusMeasurementParameters.motorsVarSigmaToTravel,
            		focusMeasurementParameters.motorsFadeSigma, 
            		focusMeasurementParameters.motorsOverShootToBalance,

            		debugLevel+1
            );
			if (numIterations==0){ // do not move, just suggest 
				result=newMotorPos;
				break; 
			}
			moveAndMaybeProbe(
					true, // just move, no "probing" in this focusing mode
					newMotorPos, // first time - will repeat measurement at the current position
					focusingMotors,
					camerasInterface,
					lensDistortionParameters,
					matchSimulatedPattern, // should not bee null
					focusMeasurementParameters,
					patternDetectParameters,
					distortionParameters,
					simulParameters,
					colorComponents,
					otfFilterParameters,
					psfParameters,
					threadsMax,
					updateStatus,
					debugLevel,
					loopDebugLevel);
			iterations++;
			// See if it is time to exit
			// is it a success?
			double lastMoveDistance=Math.abs( // here - center movement
					(newMotorPos[0]-lastMotorPos[0])+
					(newMotorPos[1]-lastMotorPos[1])+
					(newMotorPos[2]-lastMotorPos[2]))/3.0;
			if 	(lastMoveDistance<focusMeasurementParameters.minCorr){
				if (--leftFinal<=0){
					result=newMotorPos;
					break;
				}
			} else leftFinal=focusMeasurementParameters.numFinalCorr;
			if (numIterations==1) result=newMotorPos; // if the requested was just a single iteration and no other failures - it is OK (returns true)
		}
		if (DEBUG_LEVEL>0) System.out.println("Automatic fine-focus adjustment "+((result!=null)?"succeeded":"failed")+
				" after "+iterations+" iterations at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		return result;
	}

	public int [] preFocusingStepsAuto(
			int numIterations, // maximal number of iterations (0 - suggest only, do not move). When calling from the button - first time single iteration, second time - as specified
			CalibrationHardwareInterface.FocusingMotors focusingMotors,
			CalibrationHardwareInterface.CamerasInterface camerasInterface,
			Distortions.LensDistortionParameters lensDistortionParameters,
			MatchSimulatedPattern matchSimulatedPattern, // should not bee null
			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
			boolean   updateStatus,
			int       debugLevel){
		long 	  startTime=System.nanoTime();
		if (focusingMotors.sizePreFocus()==0){
			moveAndMeasureSharpness(
					null, // just add current to pre-focus history
					focusingMotors,
					camerasInterface,
					lensDistortionParameters,
					matchSimulatedPattern, // should not bee null
					focusMeasurementParameters,
					updateStatus,
					debugLevel);

		}
		int leftFinal=focusMeasurementParameters.numFinalCorrPre;
		int [] result=null;
		int iterations=0;
        for (int iterNumber=0;(iterNumber==0) || (iterNumber<numIterations);iterNumber++) {// normally will break out earlier, it is just to prevent endless loops
        	if (debugLevel>0) System.out.println("preFocusingStepsAuto(): iteration="+(iterNumber+1)+" of maximum "+numIterations);
            int [] lastMotorPos=focusingMotors.readElphel10364Motors().clone();
			int [] newMotorPos=focusingMotors.solvePreFocus(
					focusMeasurementParameters.motorsSigma,
					focusMeasurementParameters.maxAutoDistance, // Maximal motor move for each step
					focusMeasurementParameters.maxLinearStep,
            		focusMeasurementParameters.motorsOverShootToBalance,
					debugLevel
			);
			if (numIterations==0){ // do not move, just suggest 
				result=newMotorPos;
				break; 
			}
			moveAndMeasureSharpness(
					newMotorPos, // move, measure and save
					focusingMotors,
					camerasInterface,
					lensDistortionParameters,
					matchSimulatedPattern, // should not be null
					focusMeasurementParameters,
					updateStatus,
					debugLevel);
			iterations++;
			// See if it is time to exit
			// is it a success?
			double lastMoveDistance=Math.abs( // here - center movement
					(newMotorPos[0]-lastMotorPos[0])+
					(newMotorPos[1]-lastMotorPos[1])+
					(newMotorPos[2]-lastMotorPos[2]))/3.0;
			if 	(lastMoveDistance<focusMeasurementParameters.minCorrPre){
				if (DEBUG_LEVEL>1) System.out.println("lastMoveDistance="+lastMoveDistance+" (<"+focusMeasurementParameters.minCorrPre+"), leftFinal="+leftFinal);
				if (--leftFinal<=0){
					result=newMotorPos;
					break;
				}
			} else leftFinal=focusMeasurementParameters.numFinalCorrPre;
			if (numIterations==1) result=newMotorPos; // if the requested was just a single iteration and no other failures - it is OK (returns true)
		}
		if (DEBUG_LEVEL>0) System.out.println("Automatic pre-focus adjustment "+((result!=null)?"succeeded ":"failed ")+
				((iterations==0)?"to suggest move": "after "+iterations+" iterations")+" at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		return result;
	}
	
	// returns {Xmin-1,Xmax+1} - average motors

	
	public double[] ScanFocusTilt(
			//			boolean scanHysteresis, // after scanning forward, go in reverse (different number of steps to measure hysteresis
			int [] centerMotorPos, // null OK
			CalibrationHardwareInterface.FocusingMotors focusingMotors,
			CalibrationHardwareInterface.CamerasInterface camerasInterface,
			Distortions.LensDistortionParameters lensDistortionParameters,
			MatchSimulatedPattern matchSimulatedPattern, // should not bee null
			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			MatchSimulatedPattern.DistortionParameters distortionParameters,
			SimulationPattern.SimulParameters  simulParameters,
			EyesisAberrations.ColorComponents colorComponents,
			EyesisAberrations.OTFFilterParameters otfFilterParameters,
			EyesisAberrations.PSFParameters psfParameters,
			int       threadsMax,
			boolean   updateStatus,
			int       debugLevel,
			int       loopDebugLevel){
		if (centerMotorPos==null) centerMotorPos=focusingMotors.readElphel10364Motors().clone();

		boolean allOK=true;
		boolean aborted=false;
		long 	  startTime=System.nanoTime();
		if (debugLevel>0) System.out.println("Starting scanning focus in the center, number of steps="+ focusMeasurementParameters.scanNumber+
				", of them "+focusMeasurementParameters.scanNumberNegative+" in negative direction"+
				", step size="+focusMeasurementParameters.scanStep);
		int scanNegative=focusMeasurementParameters.scanNumberNegative;
		int scanPositive=focusMeasurementParameters.scanNumber-focusMeasurementParameters.scanNumberNegative-1; // not including zero
		int [] scanPos={
				centerMotorPos[0]-scanNegative * focusMeasurementParameters.scanStep,
				centerMotorPos[1]-scanNegative * focusMeasurementParameters.scanStep,
				centerMotorPos[2]-scanNegative * focusMeasurementParameters.scanStep
		};
		double centerAverage=(centerMotorPos[0]+centerMotorPos[1]+centerMotorPos[2])/3.0;
		double [] range={
				Math.floor(centerAverage- scanNegative * focusMeasurementParameters.scanStep)-1.0,
				Math.ceil(centerAverage+ scanPositive * focusMeasurementParameters.scanStep)+1.0
		};
		int [] scanPosLast=null;
		for (int numStep=0;numStep<focusMeasurementParameters.scanNumber;numStep++){
			if (debugLevel>0) System.out.println("Scanning focus in the center, step#"+(numStep+1)+" (of "+ focusMeasurementParameters.scanNumber+
					") at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
			allOK &=moveAndMaybeProbe(
					true,
					scanPos, // null OK
					focusingMotors,
					camerasInterface,
					lensDistortionParameters,
					matchSimulatedPattern, // should not bee null
					focusMeasurementParameters,
					patternDetectParameters,
					distortionParameters,
					simulParameters,
					colorComponents,
					otfFilterParameters,
					psfParameters,
					threadsMax,
					updateStatus,
					debugLevel,
					loopDebugLevel);
			// do not advance position after last measurement
			if (numStep<(focusMeasurementParameters.scanNumber-1)) for (int nm=0;nm<3;nm++) scanPos[nm]+=focusMeasurementParameters.scanStep;
			scanPosLast=scanPos.clone();
    		if (this.SYNC_COMMAND.stopRequested.get()>0){
    			aborted=true;
    			allOK=false;
    			System.out.println("Stop requested, command aborted, returning motors to initial position");
    			break;
    		}
			if (!allOK) break; // failed
		}
		if (allOK && focusMeasurementParameters.scanTiltEnable) {
			if (focusMeasurementParameters.scanTiltStepsX >1 ) { // 0 or 1 STOPS - do not scan
				double scanStepX=1.0*focusMeasurementParameters.scanTiltRangeX/(focusMeasurementParameters.scanTiltStepsX-1);
				if (debugLevel>0) System.out.println("Starting scanning tilt in X direction, number of stops="+ focusMeasurementParameters.scanTiltStepsX+
						", step size="+IJ.d2s(scanStepX,0));
				for (int numStep=0;numStep<focusMeasurementParameters.scanTiltStepsX;numStep++){
					int delta=(int) Math.round(focusMeasurementParameters.scanTiltRangeX*
							(1.0*numStep/(focusMeasurementParameters.scanTiltStepsX-1) -0.5));
					scanPos[0]=centerMotorPos[0]-delta;
					scanPos[1]=centerMotorPos[1]-delta;
					scanPos[2]=centerMotorPos[2]+delta;
					if (debugLevel>0) System.out.println("Scanning tilt in X direction, step#"+(numStep+1)+" (of "+
					focusMeasurementParameters.scanTiltStepsX+") at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
					allOK &=moveAndMaybeProbe(
							true,
							scanPos, // null OK
							focusingMotors,
							camerasInterface,
							lensDistortionParameters,
							matchSimulatedPattern, // should not bee null
							focusMeasurementParameters,
							patternDetectParameters,
							distortionParameters,
							simulParameters,
							colorComponents,
							otfFilterParameters,
							psfParameters,
							threadsMax,
							updateStatus,
							debugLevel,
							loopDebugLevel);
		    		if (this.SYNC_COMMAND.stopRequested.get()>0){
		    			aborted=true;
		    			allOK=false;
		    			System.out.println("Stop requested, command aborted, returning motors to initial position");
		    			break;
		    		}
					if (!allOK) break; // failed
				}
				if (focusMeasurementParameters.scanTiltReverse) {
					if (debugLevel>0) System.out.println("Starting reverse scanning tilt in X direction, number of stops="+ focusMeasurementParameters.scanTiltStepsX+
							", step size="+IJ.d2s(scanStepX,0));
					for (int numStep=0;numStep<focusMeasurementParameters.scanTiltStepsX;numStep++){
						int delta=(int) Math.round(focusMeasurementParameters.scanTiltRangeX*
								(1.0*numStep/(focusMeasurementParameters.scanTiltStepsX-1) -0.5));
						scanPos[0]=centerMotorPos[0]+delta;
						scanPos[1]=centerMotorPos[1]+delta;
						scanPos[2]=centerMotorPos[2]-delta;
						if (debugLevel>0) System.out.println("Reverse scanning tilt in X direction, step#"+(numStep+1)+" (of "+
						focusMeasurementParameters.scanTiltStepsX+") at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
						allOK &=moveAndMaybeProbe(
								true,
								scanPos, // null OK
								focusingMotors,
								camerasInterface,
								lensDistortionParameters,
								matchSimulatedPattern, // should not bee null
								focusMeasurementParameters,
								patternDetectParameters,
								distortionParameters,
								simulParameters,
								colorComponents,
								otfFilterParameters,
								psfParameters,
								threadsMax,
								updateStatus,
								debugLevel,
								loopDebugLevel);
			    		if (this.SYNC_COMMAND.stopRequested.get()>0){
			    			aborted=true;
			    			allOK=false;
			    			System.out.println("Stop requested, command aborted, returning motors to initial position");
			    			break;
			    		}
						if (!allOK) break; // failed
					}
				}
			}
			if (allOK && (focusMeasurementParameters.scanTiltStepsY >1 )) { // 0 or 1 STOPS - do not scan
				double scanStepY=1.0*focusMeasurementParameters.scanTiltRangeY/(focusMeasurementParameters.scanTiltStepsY-1);
				if (debugLevel>0) System.out.println("Starting scanning tilt in Y direction, number of stops="+ focusMeasurementParameters.scanTiltStepsY+
						", step size="+IJ.d2s(scanStepY,0));
				for (int numStep=0;numStep<focusMeasurementParameters.scanTiltStepsY;numStep++){
					int delta=(int) Math.round(focusMeasurementParameters.scanTiltRangeY*
							(1.0*numStep/(focusMeasurementParameters.scanTiltStepsY-1) -0.5));
					scanPos[0]=centerMotorPos[0]+delta;
					scanPos[1]=centerMotorPos[1]-delta;
					scanPos[2]=centerMotorPos[2]+0;
					if (debugLevel>0) System.out.println("Scanning tilt in Y direction, step#"+(numStep+1)+" (of "+
					focusMeasurementParameters.scanTiltStepsY+") at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
					allOK &=moveAndMaybeProbe(
							true,
							scanPos, // null OK
							focusingMotors,
							camerasInterface,
							lensDistortionParameters,
							matchSimulatedPattern, // should not bee null
							focusMeasurementParameters,
							patternDetectParameters,
							distortionParameters,
							simulParameters,
							colorComponents,
							otfFilterParameters,
							psfParameters,
							threadsMax,
							updateStatus,
							debugLevel,
							loopDebugLevel);
		    		if (this.SYNC_COMMAND.stopRequested.get()>0){
		    			aborted=true;
		    			allOK=false;
		    			System.out.println("Stop requested, command aborted, returning motors to initial position");
		    			break;
		    		}
					if (!allOK) break; // failed
				}
				if (focusMeasurementParameters.scanTiltReverse) {
					if (debugLevel>0) System.out.println("Starting reverse scanning tilt in Y direction, number of stops="+ focusMeasurementParameters.scanTiltStepsY+
							", step size="+IJ.d2s(scanStepY,0));
					for (int numStep=0;numStep<focusMeasurementParameters.scanTiltStepsY;numStep++){
						int delta=(int) Math.round(focusMeasurementParameters.scanTiltRangeY*
								(1.0*numStep/(focusMeasurementParameters.scanTiltStepsY-1) -0.5));
						scanPos[0]=centerMotorPos[0]-delta;
						scanPos[1]=centerMotorPos[1]+delta;
						scanPos[2]=centerMotorPos[2]+0;
						if (debugLevel>0) System.out.println("Reverse scanning tilt in Y direction, step#"+(numStep+1)+" (of "+
						focusMeasurementParameters.scanTiltStepsY+") at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
						allOK &=moveAndMaybeProbe(
								true,
								scanPos, // null OK
								focusingMotors,
								camerasInterface,
								lensDistortionParameters,
								matchSimulatedPattern, // should not bee null
								focusMeasurementParameters,
								patternDetectParameters,
								distortionParameters,
								simulParameters,
								colorComponents,
								otfFilterParameters,
								psfParameters,
								threadsMax,
								updateStatus,
								debugLevel,
								loopDebugLevel);
			    		if (this.SYNC_COMMAND.stopRequested.get()>0){
			    			aborted=true;
			    			allOK=false;
			    			System.out.println("Stop requested, command aborted, returning motors to initial position");
			    			break;
			    		}
						if (!allOK) break; // failed
					}
				}
			}
		}

		if (allOK && focusMeasurementParameters.scanHysteresis && (scanPosLast!=null)){
			focusingMotors.moveElphel10364Motors( // return to last direct scan position
					true, //boolean wait,
					scanPosLast,
					0.0, //double sleep,
					true, //boolean showStatus,
					"",   //String message,
					false); //focusMeasurementParameters.compensateHysteresis); //boolean hysteresis)

			double hystStep=((double) focusMeasurementParameters.scanStep*(focusMeasurementParameters.scanNumber-1))/focusMeasurementParameters.scanHysteresisNumber; // errors will accumulate, but that's OK
			for (int numStep=0;numStep<focusMeasurementParameters.scanHysteresisNumber;numStep++){
				if (debugLevel>0) System.out.println("Scanning focus in reverse direction for hysteresis estimation, step#"+(numStep+1)+" (of "+ focusMeasurementParameters.scanHysteresisNumber+
						") at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
				int [] revPosition=new int[scanPos.length];
				for (int nm=0;nm<3;nm++) revPosition[nm]=scanPosLast[nm]-(int) Math.round(hystStep*(numStep+1));
				allOK &=moveAndMaybeProbe(
						true, // boolean noHysteresis,
						true,
						revPosition, // null OK
						focusingMotors,
						camerasInterface,
						lensDistortionParameters,
						matchSimulatedPattern, // should not bee null
						focusMeasurementParameters,
						patternDetectParameters,
						distortionParameters,
						simulParameters,
						colorComponents,
						otfFilterParameters,
						psfParameters,
						threadsMax,
						updateStatus,
						debugLevel,
						loopDebugLevel);
	    		if (this.SYNC_COMMAND.stopRequested.get()>0){
	    			aborted=true;
	    			allOK=false;
	    			System.out.println("Stop requested, command aborted, returning motors to initial position");
	    			break;
	    		}
				if (!allOK) break; // failed
			}		
		}
		if (aborted) {
			System.out.println("Returning motors to initial position due to command aborted at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
			focusingMotors.moveElphel10364Motors( // return to last direct scan position
					true, //boolean wait,
					centerMotorPos,
					0.0, //double sleep,
					true, //boolean showStatus,
					"",   //String message,
					false); //focusMeasurementParameters.compensateHysteresis); //boolean hysteresis)
			return null;
		}
		if (allOK) {
			if (focusMeasurementParameters.scanMeasureLast) {
			allOK &= moveAndMaybeProbe(
					true,
					centerMotorPos, // null OK
					focusingMotors,
					camerasInterface,
					lensDistortionParameters,
					matchSimulatedPattern, // should not bee null
					focusMeasurementParameters,
					patternDetectParameters,
					distortionParameters,
					simulParameters,
					colorComponents,
					otfFilterParameters,
					psfParameters,
					threadsMax,
					updateStatus,
					debugLevel,
					loopDebugLevel);
			} else { // just move, no measuring
				System.out.println("Returning motors to initial position at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
				focusingMotors.moveElphel10364Motors( // return to last direct scan position
						true, //boolean wait,
						centerMotorPos,
						0.0, //double sleep,
						true, //boolean showStatus,
						"",   //String message,
						false); //focusMeasurementParameters.compensateHysteresis); //boolean hysteresis)
			}
		}

		if (debugLevel>0) System.out.println("Scanning focus in the center, number of steps="+ focusMeasurementParameters.scanNumber+
				", step size="+focusMeasurementParameters.scanStep+" finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+
				", allOK="+allOK);
		// focusMeasurementParameters.scanStep	
		// focusMeasurementParameters.scanNumber
		return allOK?range:null;
	}	

	
	
	public double[] ScanFocus(
//			boolean scanHysteresis, // after scanning forward, go in reverse (different number of steps to measure hysteresis
			int [] centerMotorPos, // null OK
			CalibrationHardwareInterface.FocusingMotors focusingMotors,
			CalibrationHardwareInterface.CamerasInterface camerasInterface,
		    Distortions.LensDistortionParameters lensDistortionParameters,
			MatchSimulatedPattern matchSimulatedPattern, // should not bee null
			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			MatchSimulatedPattern.DistortionParameters distortionParameters,
			SimulationPattern.SimulParameters  simulParameters,
			EyesisAberrations.ColorComponents colorComponents,
			EyesisAberrations.OTFFilterParameters otfFilterParameters,
			EyesisAberrations.PSFParameters psfParameters,
			int       threadsMax,
			boolean   updateStatus,
			int       debugLevel,
			int       loopDebugLevel){
		if (centerMotorPos==null) centerMotorPos=focusingMotors.readElphel10364Motors().clone();

		boolean allOK=true;
		long 	  startTime=System.nanoTime();
		if (debugLevel>0) System.out.println("Starting scanning focus in the center, number of steps="+ focusMeasurementParameters.scanNumber+
				", of them "+focusMeasurementParameters.scanNumberNegative+" in negative direction"+
				", step size="+focusMeasurementParameters.scanStep);
		int scanNegative=focusMeasurementParameters.scanNumberNegative;
		int scanPositive=focusMeasurementParameters.scanNumber-focusMeasurementParameters.scanNumberNegative-1; // not including zero

//		int range=(focusMeasurementParameters.scanNumber-1) * focusMeasurementParameters.scanStep;
		int [] scanPos={
//				centerMotorPos[0]-((focusMeasurementParameters.scanNumber-1) * focusMeasurementParameters.scanStep)/2,
//				centerMotorPos[1]-((focusMeasurementParameters.scanNumber-1) * focusMeasurementParameters.scanStep)/2,
//				centerMotorPos[2]-((focusMeasurementParameters.scanNumber-1) * focusMeasurementParameters.scanStep)/2

				centerMotorPos[0]-scanNegative * focusMeasurementParameters.scanStep,
				centerMotorPos[1]-scanNegative * focusMeasurementParameters.scanStep,
				centerMotorPos[2]-scanNegative * focusMeasurementParameters.scanStep
				};
		double centerAverage=(centerMotorPos[0]+centerMotorPos[1]+centerMotorPos[2])/3.0;
		double [] range={
//				Math.floor(centerAverage-((focusMeasurementParameters.scanNumber-1) * focusMeasurementParameters.scanStep)/2)-1.0,
//				Math.ceil(centerAverage+((focusMeasurementParameters.scanNumber-1) * focusMeasurementParameters.scanStep)/2)+1.0
				Math.floor(centerAverage- scanNegative * focusMeasurementParameters.scanStep)-1.0,
				Math.ceil(centerAverage+ scanPositive * focusMeasurementParameters.scanStep)+1.0
		};
	for (int numStep=0;numStep<focusMeasurementParameters.scanNumber;numStep++){
		if (debugLevel>0) System.out.println("Scanning focus in the center, step#"+(numStep+1)+" (of "+ focusMeasurementParameters.scanNumber+
				") at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		allOK &=moveAndMaybeProbe(
				true,
				scanPos, // null OK
				focusingMotors,
				camerasInterface,
			    lensDistortionParameters,
				matchSimulatedPattern, // should not bee null
				focusMeasurementParameters,
				patternDetectParameters,
				distortionParameters,
				simulParameters,
				colorComponents,
				otfFilterParameters,
				psfParameters,
				threadsMax,
				updateStatus,
				debugLevel,
				loopDebugLevel);
		// do not advance position after last measurement
		if (numStep<(focusMeasurementParameters.scanNumber-1)) for (int nm=0;nm<3;nm++) scanPos[nm]+=focusMeasurementParameters.scanStep;
		if (!allOK) break; // failed 
	}
	if (focusMeasurementParameters.scanHysteresis){
		double hystStep=((double) focusMeasurementParameters.scanStep*(focusMeasurementParameters.scanNumber-1))/focusMeasurementParameters.scanHysteresisNumber; // errors will accumulate, but that's OK
		for (int numStep=0;numStep<focusMeasurementParameters.scanHysteresisNumber;numStep++){
			if (debugLevel>0) System.out.println("Scanning focus in reverse direction for hysteresis estimation, step#"+(numStep+1)+" (of "+ focusMeasurementParameters.scanHysteresisNumber+
					") at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
			int [] revPosition=new int[scanPos.length];
			for (int nm=0;nm<3;nm++) revPosition[nm]=scanPos[nm]-(int) Math.round(hystStep*(numStep+1));
			allOK &=moveAndMaybeProbe(
					true, // boolean noHysteresis,
					true,
					revPosition, // null OK
					focusingMotors,
					camerasInterface,
				    lensDistortionParameters,
					matchSimulatedPattern, // should not bee null
					focusMeasurementParameters,
					patternDetectParameters,
					distortionParameters,
					simulParameters,
					colorComponents,
					otfFilterParameters,
					psfParameters,
					threadsMax,
					updateStatus,
					debugLevel,
					loopDebugLevel);
			if (!allOK) break; // failed 
		}		
	}
	allOK &= moveAndMaybeProbe(
			true,
			centerMotorPos, // null OK
			focusingMotors,
			camerasInterface,
		    lensDistortionParameters,
			matchSimulatedPattern, // should not bee null
			focusMeasurementParameters,
			patternDetectParameters,
			distortionParameters,
			simulParameters,
			colorComponents,
			otfFilterParameters,
			psfParameters,
			threadsMax,
			updateStatus,
			debugLevel,
			loopDebugLevel);

	
	if (debugLevel>0) System.out.println("Scanning focus in the center, number of steps="+ focusMeasurementParameters.scanNumber+
			", step size="+focusMeasurementParameters.scanStep+" finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
// focusMeasurementParameters.scanStep	
		// focusMeasurementParameters.scanNumber
	return allOK?range:null;
	}	
	
	public boolean moveAndMaybeProbe(
			boolean justMove, // just move, no probing around
			int [] newMotorPos, // null OK
			CalibrationHardwareInterface.FocusingMotors focusingMotors,
			CalibrationHardwareInterface.CamerasInterface camerasInterface,
		    Distortions.LensDistortionParameters lensDistortionParameters,
			MatchSimulatedPattern matchSimulatedPattern, // should not bee null
			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			MatchSimulatedPattern.DistortionParameters distortionParameters,
			SimulationPattern.SimulParameters  simulParameters,
			EyesisAberrations.ColorComponents colorComponents,
			EyesisAberrations.OTFFilterParameters otfFilterParameters,
			EyesisAberrations.PSFParameters psfParameters,
			int       threadsMax,
			boolean   updateStatus,
			int       debugLevel,
			int       loopDebugLevel){
		return  moveAndMaybeProbe(
				false, // apply hysteresis correction as specified in focusMeasurementParameters
				justMove, // just move, no probing around
				newMotorPos, // null OK
				focusingMotors,
				camerasInterface,
			    lensDistortionParameters,
				matchSimulatedPattern, // should not bee null
				focusMeasurementParameters,
				patternDetectParameters,
				distortionParameters,
				simulParameters,
				colorComponents,
				otfFilterParameters,
				psfParameters,
				threadsMax,
				updateStatus,
				debugLevel,
				loopDebugLevel);
	}
	public boolean moveAndMaybeProbeOld(
			boolean noHysteresis,
			boolean justMove, // just move, no probing around
			int [] newMotorPos, // null OK
			CalibrationHardwareInterface.FocusingMotors focusingMotors,
			CalibrationHardwareInterface.CamerasInterface camerasInterface,
		    Distortions.LensDistortionParameters lensDistortionParameters,
			MatchSimulatedPattern matchSimulatedPattern, // should not bee null
			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			MatchSimulatedPattern.DistortionParameters distortionParameters,
			SimulationPattern.SimulParameters  simulParameters,
			EyesisAberrations.ColorComponents colorComponents,
			EyesisAberrations.OTFFilterParameters otfFilterParameters,
			EyesisAberrations.PSFParameters psfParameters,
			int       threadsMax,
			boolean   updateStatus,
			int       debugLevel,
			int       loopDebugLevel){
		boolean noMove=false;
		if (newMotorPos==null) {
			newMotorPos=focusingMotors.readElphel10364Motors().clone();
			justMove=true;
			noMove=true;
		}
		
		if (focusingMotors.distFromProbed(newMotorPos)<focusMeasurementParameters.reProbeDistance){ // distFromProbed() may return NaN, and will - first time
			justMove=true; // no probing is needed
		}
		
		int [][] seqMove={{0,0,0}};
		int [][] seqShort={
				{(int) focusMeasurementParameters.probeStep,0,0},
				{0,(int) focusMeasurementParameters.probeStep,0},
				{0,0,(int) focusMeasurementParameters.probeStep},
				{-(int) (focusMeasurementParameters.probeStep/3),-(int) (focusMeasurementParameters.probeStep/3),-(int) (focusMeasurementParameters.probeStep/3)},
				{0,0,0}};
		int [][] seqSymm={
				{  (int) focusMeasurementParameters.probeStep, 0,0},
				{-((int) focusMeasurementParameters.probeStep),0,0},
				{0,  (int) focusMeasurementParameters.probeStep,0},
				{0,-((int) focusMeasurementParameters.probeStep),0},
				{0,0,  (int) focusMeasurementParameters.probeStep},
				{0,0,-((int) focusMeasurementParameters.probeStep)},
				{0,0,0}};
		int [][]seq= justMove?seqMove:(focusMeasurementParameters.probeSymmetrical?seqSymm:seqShort);
		ImagePlus imp;
//System.out.println("target position: m1="+newMotorPos[0]+" m2="+newMotorPos[1]+" m3="+newMotorPos[2]);		
		for (int seqNum=0;seqNum<seq.length;seqNum++){
        	if (!justMove && (debugLevel>0)) System.out.println("moveAndMaybeProbe() iteration="+(seqNum+1)+" of total "+seq.length);
			for (int i=0;i<seq[0].length;i++) seq[seqNum][i]+=newMotorPos[i];
//System.out.println("# "+seqNum+":  m1="+seq[seqNum][0]+" m2="+seq[seqNum][1]+" m3="+seq[seqNum][2]);
			if (!noMove)focusingMotors.moveElphel10364Motors(
				  true, //boolean wait,
				  seq[seqNum],
				  0.0, //double sleep,
				  true, //boolean showStatus,
				  "",   //String message,
				  !noHysteresis); //false); //focusMeasurementParameters.compensateHysteresis); //boolean hysteresis)
//			System.out.println("A#"+focusingMotors.historySize()+": "+seq[seqNum][0]+", "+seq[seqNum][1]+", "+seq[seqNum][2]);
			focusMeasurementParameters.sensorTemperature=camerasInterface.getSensorTemperature(0,focusMeasurementParameters.EEPROM_channel);
			imp= camerasInterface.acquireSingleImage (
					false, //boolean useLasers,
					updateStatus);
			String ts=(String) imp.getProperty("timestamp");
			if (debugLevel>0) System.out.println("Image timestamp="+((ts==null)?"null":ts));
		
/*
  			if (imp==null){
 				   String msg="Failed to get camera image\nProcess canceled";
				   IJ.showMessage("Error",msg);
				   throw new IllegalArgumentException (msg);
			}
*/			
			if (matchSimulatedPattern==null) {
				   String msg="matchSimulatedPattern is null - it should be initialized before calling this method";
				   IJ.showMessage("Error",msg);
				   throw new IllegalArgumentException (msg);
			}
			matchSimulatedPattern.debugLevel=debugLevel;
			if (debugLevel>1){
				System.out.println("Samples Map:\n"+
						focusMeasurementParameters.showSamplesMap(
								lensDistortionParameters.px0, // pixel coordinate of the the optical center
								lensDistortionParameters.py0, // pixel coordinate of the the optical center
								focusMeasurementParameters.numInCenter));
			}
			double [][][][][] rFullResults=new double [1][][][][];
			double [][] metrics=measurePSFMetrics(
					imp,
					lensDistortionParameters,
					matchSimulatedPattern,
					focusMeasurementParameters,
					patternDetectParameters,
					distortionParameters,
					simulParameters,
					colorComponents,
					otfFilterParameters,
					psfParameters,
					rFullResults,
					threadsMax,
					updateStatus,
					debugLevel,
					loopDebugLevel);
			focusingMotors.addToHistory(ts,focusMeasurementParameters.sensorTemperature,metrics,rFullResults[0]);
			
// end of measure and save to history			
			
//			if (Double.isNaN(metrics[6][3])){ //????
			if (Double.isNaN(metrics[6][6])){ // average colors, sharpness in the center (others might be undefined)
				int ca=6;
				String msg="#"+focusingMotors.historySize()+": "+seq[seqNum][0]+", "+seq[seqNum][1]+", "+seq[seqNum][2]+
						": Far/Near="+IJ.d2s(metrics[ca][0],3)+
						"  Tilt X="+IJ.d2s(metrics[ca][1],3)+
						"  Tilt Y="+IJ.d2s(metrics[ca][2],3)+
						"  R50% average="+IJ.d2s(metrics[ca][3],3)+" sensor pixels,"+
						"  A50% average="+IJ.d2s(metrics[ca][4],3)+" sensor pixels,"+
						"  B50% average="+IJ.d2s(metrics[ca][5],3)+" sensor pixels,"+
					    "  R50%Center="+IJ.d2s(metrics[ca][6],3)+" sensor pixels";
    			IJ.showMessage("Errror",msg);
    			throw new IllegalArgumentException (msg);

			}
//			System.out.println("B#"+focusingMotors.historySize()+": "+seq[seqNum][0]+", "+seq[seqNum][1]+", "+seq[seqNum][2] + ", debugLevel="+debugLevel);
			if (debugLevel>0){
// see if lens is calibrated
				double [] resolutions={1.0/metrics[1][6],1.0/metrics[5][6],1.0/metrics[2][6]}; // R,G,B
				double fDistance=focusingMotors.focusingHistory.getLensDistance(
						resolutions, // {R-sharpness,G-sharpness,B-sharpness}
		    			true, // boolean absolute, // return absolutely calibrated data
		    			focusMeasurementParameters.lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
		    			focusMeasurementParameters.lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
		    			debugLevel
		    			);
				
				int ca=6;
				int [] actualPosition=focusingMotors.focusingHistory.getPosition();
				System.out.println("#"+focusingMotors.historySize()+": "+actualPosition[0]+", "+actualPosition[1]+", "+actualPosition[2]+
						": fDistance="+IJ.d2s(fDistance,3)+"um"+
						"  Far/Near="+IJ.d2s(metrics[ca][0],3)+
						"  Tilt X="+IJ.d2s(metrics[ca][1],3)+
						"  Tilt Y="+IJ.d2s(metrics[ca][2],3)+
						"  R50% average="+IJ.d2s(metrics[ca][3],3)+" sensor pixels,"+
						"  A50% average="+IJ.d2s(metrics[ca][4],3)+" sensor pixels,"+
						"  B50% average="+IJ.d2s(metrics[ca][5],3)+" sensor pixels,"+
					    "  R50%Center="+IJ.d2s(metrics[ca][6],3)+" sensor pixels");
				if (debugLevel>1){
					String [] compColors={"green","red","blue","green","green","green","AVERAGE"};
					for (int c=0;c<metrics.length-1;c++) if (metrics[c]!=null){
						System.out.println(compColors[c]+": Far/Near="+IJ.d2s(metrics[c][0],3)+
								"  Tilt X="+IJ.d2s(metrics[c][1],3)+
								"  Tilt Y="+IJ.d2s(metrics[c][2],3)+
								"  R50% average="+IJ.d2s(metrics[c][3],3)+" sensor pixels,"+
								"  A50% average="+IJ.d2s(metrics[c][4],3)+" sensor pixels,"+
								"  B50% average="+IJ.d2s(metrics[c][5],3)+" sensor pixels,"+
							    "  R50%Center="+IJ.d2s(metrics[c][6],3)+" sensor pixels,"+
								"  component weight="+IJ.d2s(100*metrics[c][7],1)+"%");
					}
				}
			}
		}
		if (!justMove) focusingMotors.setLastProbed();
		return true;
	}
	public boolean moveAndMaybeProbe(
			boolean noHysteresis,
			boolean justMove, // just move, no probing around
			int [] newMotorPos, // null OK
			CalibrationHardwareInterface.FocusingMotors focusingMotors,
			CalibrationHardwareInterface.CamerasInterface camerasInterface,
		    Distortions.LensDistortionParameters lensDistortionParameters,
			MatchSimulatedPattern matchSimulatedPattern, // should not bee null
			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			MatchSimulatedPattern.DistortionParameters distortionParameters,
			SimulationPattern.SimulParameters  simulParameters,
			EyesisAberrations.ColorComponents colorComponents,
			EyesisAberrations.OTFFilterParameters otfFilterParameters,
			EyesisAberrations.PSFParameters psfParameters,
			int       threadsMax,
			boolean   updateStatus,
			int       debugLevel,
			int       loopDebugLevel){
		boolean noMove=false;
		if (newMotorPos==null) {
			newMotorPos=focusingMotors.readElphel10364Motors().clone();
			justMove=true;
			noMove=true;
		}
		
		if (focusingMotors.distFromProbed(newMotorPos)<focusMeasurementParameters.reProbeDistance){ // distFromProbed() may return NaN, and will - first time
			justMove=true; // no probing is needed
		}
		
		if (!justMove && focusMeasurementParameters.parallelBeforeProbing){
			moveMeasureAndSave(
					noHysteresis,
					newMotorPos,
					focusingMotors,
					camerasInterface,
					lensDistortionParameters,
					matchSimulatedPattern, // should not bee null
					focusMeasurementParameters,
					patternDetectParameters,
					distortionParameters,
					simulParameters,
					colorComponents,
					otfFilterParameters,
					psfParameters,
					threadsMax,
					updateStatus,
					1, //debugLevel,
					loopDebugLevel);
			if (debugLevel>0) System.out.println("Re-adjusting focus by parallel move before probing around");
			newMotorPos=focusingMotors.focusingHistory.focusReadjustStep(
					focusMeasurementParameters.targetMicrons, //double targetMicrons, // target focal distance
					focusMeasurementParameters.goodDistanceSigma, // double micronsFade, // reduce influence of far points
					focusMeasurementParameters.lensDistanceWeightK, // double lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
					focusMeasurementParameters.lensDistanceWeightY, // double lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
					debugLevel//int debugLevel
	     			);			

		}
		int [][] seq={newMotorPos};
		if (!justMove) seq=focusingMotors.focusingHistory.probeSequence(
				newMotorPos, // should not be null
				focusMeasurementParameters.probeSymmetrical,
				focusMeasurementParameters.probe_M1M2M3, // sigma for decay for average of the 3 motors: (M1+M2+M3)/3
				focusMeasurementParameters.probe_M3_M1M2, // sigma for decay for M3 opposite to M1 and M2: M3-(M1+M2)/2
				focusMeasurementParameters.probe_M2_M1);
//System.out.println("target position: m1="+newMotorPos[0]+" m2="+newMotorPos[1]+" m3="+newMotorPos[2]);		
		for (int seqNum=0;seqNum<seq.length;seqNum++){
			if (!justMove && (debugLevel>0)) System.out.println("moveAndMaybeProbe() iteration="+(seqNum+1)+" of total "+seq.length);
			moveMeasureAndSave(
					noHysteresis,
					(noMove?null:seq[seqNum]), // null OK
					focusingMotors,
					camerasInterface,
					lensDistortionParameters,
					matchSimulatedPattern, // should not bee null
					focusMeasurementParameters,
					patternDetectParameters,
					distortionParameters,
					simulParameters,
					colorComponents,
					otfFilterParameters,
					psfParameters,
					threadsMax,
					updateStatus,
					debugLevel, //debugLevel, // 1,
					loopDebugLevel);
		}
		if (!justMove) focusingMotors.setLastProbed();
		return true;
	}

	public void moveMeasureAndSave(
			boolean noHysteresis,
			int [] newMotorPos, // null OK
			CalibrationHardwareInterface.FocusingMotors focusingMotors,
			CalibrationHardwareInterface.CamerasInterface camerasInterface,
		    Distortions.LensDistortionParameters lensDistortionParameters,
			MatchSimulatedPattern matchSimulatedPattern, // should not bee null
			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			MatchSimulatedPattern.DistortionParameters distortionParameters,
			SimulationPattern.SimulParameters  simulParameters,
			EyesisAberrations.ColorComponents colorComponents,
			EyesisAberrations.OTFFilterParameters otfFilterParameters,
			EyesisAberrations.PSFParameters psfParameters,
			int       threadsMax,
			boolean   updateStatus,
			int       debugLevel,
			int       loopDebugLevel){
//		System.out.println(">"+focusingMotors.historySize()+": "+focusingMotors.curpos[0]+", "+focusingMotors.curpos[1]+", "+focusingMotors.curpos[2]);
		boolean noMove=false;
		if (newMotorPos==null) {
			newMotorPos=focusingMotors.readElphel10364Motors().clone();
			noMove=true;
		}
		if (!noMove)focusingMotors.moveElphel10364Motors(
				true, //boolean wait,
				newMotorPos,
				0.0, //double sleep,
				true, //boolean showStatus,
				"",   //String message,
				!noHysteresis); //false); //focusMeasurementParameters.compensateHysteresis); //boolean hysteresis)
		focusMeasurementParameters.sensorTemperature=camerasInterface.getSensorTemperature(0,focusMeasurementParameters.EEPROM_channel);
		ImagePlus imp= camerasInterface.acquireSingleImage (
				false, //boolean useLasers,
				updateStatus);
		String ts=(String) imp.getProperty("timestamp");
		if (debugLevel>0) System.out.println("Image timestamp="+((ts==null)?"null":ts));

/*
  		if (imp==null){
 			String msg="Failed to get camera image\nProcess canceled";
			IJ.showMessage("Error",msg);
			throw new IllegalArgumentException (msg);
		}
*/		
		if (matchSimulatedPattern==null) {
			String msg="matchSimulatedPattern is null - it should be initialized before calling this method";
			IJ.showMessage("Error",msg);
			throw new IllegalArgumentException (msg);
		}
		matchSimulatedPattern.debugLevel=debugLevel;
		if (debugLevel>2){
			System.out.println("Samples Map:\n"+
					focusMeasurementParameters.showSamplesMap(
							lensDistortionParameters.px0, // pixel coordinate of the the optical center
							lensDistortionParameters.py0, // pixel coordinate of the the optical center
							focusMeasurementParameters.numInCenter));
		}
		double [][][][][] rFullResults=new double [1][][][][];
		double [][] metrics=measurePSFMetrics(
				imp,
				lensDistortionParameters,
				matchSimulatedPattern,
				focusMeasurementParameters,
				patternDetectParameters,
				distortionParameters,
				simulParameters,
				colorComponents,
				otfFilterParameters,
				psfParameters,
				rFullResults,
				threadsMax,
				updateStatus,
				debugLevel,
				loopDebugLevel);
		//			System.out.println(">>"+focusingMotors.historySize()+": "+focusingMotors.curpos[0]+", "+focusingMotors.curpos[1]+", "+focusingMotors.curpos[2]);

		focusingMotors.addToHistory(ts,focusMeasurementParameters.sensorTemperature,metrics,rFullResults[0]);
//		System.out.println("focusMeasurementParameters.lensDistanceWeightK="+focusMeasurementParameters.lensDistanceWeightK);
//		System.out.println("focusMeasurementParameters.lensDistanceWeightY="+focusMeasurementParameters.lensDistanceWeightY);
		if ((debugLevel>0) && (metrics!=null)){
			// see if lens is calibrated
			double [] resolutions={1.0/metrics[1][6],1.0/metrics[5][6],1.0/metrics[2][6]}; // R,G,B
			double fDistance=focusingMotors.focusingHistory.getLensDistance(
					resolutions, // {R-sharpness,G-sharpness,B-sharpness}
					true, // boolean absolute, // return absolutely calibrated data
					focusMeasurementParameters.lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
					focusMeasurementParameters.lensDistanceWeightY, // R-frac, B-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
					debugLevel
			);

			int ca=6;
			int [] actualPosition=focusingMotors.focusingHistory.getPosition();
			System.out.println("##"+focusingMotors.historySize()+": "+actualPosition[0]+", "+actualPosition[1]+", "+actualPosition[2]+
					": fDistance="+IJ.d2s(fDistance,3)+"um"+
					"  Far/Near="+IJ.d2s(metrics[ca][0],3)+
					"  Tilt X="+IJ.d2s(metrics[ca][1],3)+
					"  Tilt Y="+IJ.d2s(metrics[ca][2],3)+
					"  R50% average="+IJ.d2s(metrics[ca][3],3)+" sensor pixels,"+
					"  A50% average="+IJ.d2s(metrics[ca][4],3)+" sensor pixels,"+
					"  B50% average="+IJ.d2s(metrics[ca][5],3)+" sensor pixels,"+
					"  R50%Center="+IJ.d2s(metrics[ca][6],3)+" sensor pixels"+
					" temp="+focusMeasurementParameters.sensorTemperature);
			if (debugLevel>1){
				String [] compColors={"green","red","blue","green","green","green","AVERAGE"};
				for (int c=0;c<metrics.length-1;c++) if (metrics[c]!=null){
					System.out.println(compColors[c]+": Far/Near="+IJ.d2s(metrics[c][0],3)+
							"  Tilt X="+IJ.d2s(metrics[c][1],3)+
							"  Tilt Y="+IJ.d2s(metrics[c][2],3)+
							"  R50% average="+IJ.d2s(metrics[c][3],3)+" sensor pixels,"+
							"  A50% average="+IJ.d2s(metrics[c][4],3)+" sensor pixels,"+
							"  B50% average="+IJ.d2s(metrics[c][5],3)+" sensor pixels,"+
							"  R50%Center="+IJ.d2s(metrics[c][6],3)+" sensor pixels,"+
							"  component weight="+IJ.d2s(100*metrics[c][7],1)+"%");
				}
			}
			focusMeasurementParameters.result_fDistance=fDistance; // last measured focal distance
			focusMeasurementParameters.result_tiltX=metrics[ca][1]; // last measured tilt X
			focusMeasurementParameters.result_tiltY=metrics[ca][2]; // last measured tilt Y
			focusMeasurementParameters.result_R50=metrics[ca][3];   // last measured R50 (average PSF 50% level radius, pixels - somewhat larged than actual because of measurement settings)
			focusMeasurementParameters.result_A50=metrics[ca][4];   // last measured A50 (simailar, but R^2 are averaged) 
			focusMeasurementParameters.result_B50=metrics[ca][5];   // last measured B50 (simailar, but R^4 are averaged)
			focusMeasurementParameters.result_RC50=metrics[ca][6];  // last measured RC50(R50 calculated only for the 2 center samples)
		}

	}
	
	
	
	public boolean moveAndMeasureSharpness(
			int [] newMotorPos, // null OK
			CalibrationHardwareInterface.FocusingMotors focusingMotors,
			CalibrationHardwareInterface.CamerasInterface camerasInterface,
			Distortions.LensDistortionParameters lensDistortionParameters,
			MatchSimulatedPattern matchSimulatedPattern, // should not bee null
			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
			boolean   updateStatus,
			int       debugLevel){
		if (newMotorPos==null){
			newMotorPos=focusingMotors.readElphel10364Motors().clone();
			if (debugLevel>0) System.out.println("moveAndMeasureSharpness(): Measuring at "+newMotorPos[0]+", "+newMotorPos[1]+", "+newMotorPos[2]);
		}else {
			if (debugLevel>0) System.out.println("moveAndMeasureSharpness(): Moving to "+newMotorPos[0]+", "+newMotorPos[1]+", "+newMotorPos[2]);
			focusingMotors.moveElphel10364Motors(
				true, //boolean wait,
				newMotorPos,
				0.0, //double sleep,
				true, //boolean showStatus,
				"",   //String message,
				focusMeasurementParameters.compensateHysteresis); //boolean hysteresis)
		}
		focusMeasurementParameters.sensorTemperature=camerasInterface.getSensorTemperature(0,FOCUS_MEASUREMENT_PARAMETERS.EEPROM_channel);
		ImagePlus imp= camerasInterface.acquireSingleImage (
				false, //boolean useLasers,
				updateStatus);
		if (imp==null){
			String msg="Failed to get camera image\nProcess canceled";
			IJ.showMessage("Error",msg);
			throw new IllegalArgumentException (msg);
		}
		if (matchSimulatedPattern==null) {
			String msg="matchSimulatedPattern is null - it should be initialized before calling this method";
			IJ.showMessage("Error",msg);
			throw new IllegalArgumentException (msg);
		}
		matchSimulatedPattern.debugLevel=debugLevel;

		double sharpnessOld=matchSimulatedPattern.focusQualityOld(
				imp,
				focusMeasurementParameters.sampleSize, // will be twice the regualr FFT size
				5,
				lensDistortionParameters.px0,
				lensDistortionParameters.py0,
				debugLevel);
		if (debugLevel>0) System.out.println("Focus qualityOld="+sharpnessOld);
		double sharpnessOld1=matchSimulatedPattern.focusQualityOld1(
				imp,
				focusMeasurementParameters.sampleSize, // will be twice the regualr FFT size
				5,
				lensDistortionParameters.px0,
				lensDistortionParameters.py0,
				debugLevel);
		if (debugLevel>0) System.out.println("Focus qualityOld1="+sharpnessOld1);
		double sharpness=matchSimulatedPattern.focusQuality(
				imp,
				focusMeasurementParameters.sampleSize, // will be twice the regualr FFT size
				5,
				lensDistortionParameters.px0,
				lensDistortionParameters.py0,
				debugLevel);
		if (debugLevel>0) System.out.println("Focus quality="+sharpness);
		focusingMotors.addPreFocus(sharpness);
		return true;
	}

	public double [][] measurePSFMetrics(
			ImagePlus imp_sel,
		    Distortions.LensDistortionParameters lensDistortionParameters,
			MatchSimulatedPattern matchSimulatedPattern,
			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			MatchSimulatedPattern.DistortionParameters distortionParameters,
			SimulationPattern.SimulParameters  simulParameters,
			EyesisAberrations.ColorComponents colorComponents,
			EyesisAberrations.OTFFilterParameters otfFilterParameters,
			EyesisAberrations.PSFParameters psfParameters,
			double [][][][][] returnFullResults, // null OK =- will return [0][i][j][0,1,2(r/g/b)][0,1]
			int       threadsMax,
			boolean   updateStatus,
			int       debugLevel,
			int       loopDebugLevel){
		//"Focusing Acquire PSF"		
		ImagePlus imp_src=imp_sel;
		LENS_ADJUSTMENT.updateFocusGrid(
				lensDistortionParameters.px0, // pixel coordinate of the the optical center
				lensDistortionParameters.py0, // pixel coordinate of the the optical center
				imp_src,
				matchSimulatedPattern,
				distortionParameters,
				focusMeasurementParameters,
				patternDetectParameters,
				null, //LASER_POINTERS
				simulParameters,
				true, // remove non-PSF areas
				colorComponents.equalizeGreens,
				threadsMax,
				updateStatus,
				loopDebugLevel);
		double [][][][] psf_kernels= focusPSF(
				lensDistortionParameters.px0, // pixel coordinate of the the optical center
				lensDistortionParameters.py0, // pixel coordinate of the the optical center
				imp_src,
				matchSimulatedPattern,
				focusMeasurementParameters,
				simulParameters,
				0.0, // double overexposedMaxFraction, ( MULTIFILE_PSF.overexposedMaxFraction, )
				colorComponents,
				otfFilterParameters,
				psfParameters, // step of the new map (should be multiple of map step)
				threadsMax,
				updateStatus,
				loopDebugLevel); // debug level
//		DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
    	double [][][] sampleCoord=focusMeasurementParameters.sampleCoordinates(
				lensDistortionParameters.px0, // pixel coordinate of the the optical center
				lensDistortionParameters.py0); // pixel coordinate of the the optical center
    	
		double [][][] anglePerPixel=pixToAngles(
				sampleCoord,
				lensDistortionParameters.px0, // pixel coordinate of the the optical center
				lensDistortionParameters.py0, // pixel coordinate of the the optical center
				lensDistortionParameters.pixelSize, // in microns
				lensDistortionParameters.focalLength, // in mm
				lensDistortionParameters.distortionRadius, //=  2.8512; // mm - half width of the sensor
				lensDistortionParameters.distortionA5, //r^5 (normalized to focal length or to sensor half width?)
				lensDistortionParameters.distortionA, // r^4 (normalized to focal length or to sensor half width?)
				lensDistortionParameters.distortionB, // r^3
				lensDistortionParameters.distortionC); // r^2
		double [] componentWeights={
				1.0,
				focusMeasurementParameters.weightRatioRedToGreen,
				focusMeasurementParameters.weightRatioBlueToGreen,
				1.0,
				1.0,
				1.0};
		double [][][][] fullResults=new double [sampleCoord.length][sampleCoord[0].length][][];
		double [][]metrics= extractPSFMetrics(
				sampleCoord,
				lensDistortionParameters.px0, // pixel coordinate of the the optical center
				lensDistortionParameters.py0, // pixel coordinate of the the optical center
				psf_kernels,
				focusMeasurementParameters,
				componentWeights,
				fullResults,// double [][][][] fullResults, // [y][x][color],should be null or match samples, last dimension {R50,  
				debugLevel // +2 // +0
				);
		int [] rgbChn={1,5,2};
		// modify fullResults to use just 3 R,G,B channels 
		if (returnFullResults!=null) {
			returnFullResults[0]=new double [fullResults.length][fullResults[0].length][][];
			for (int i=0;i<fullResults.length;i++) for (int j=0;j<fullResults[i].length;j++){
				if (fullResults[i][j]!=null){
					returnFullResults[0][i][j]=new double [rgbChn.length][];
					for (int c=0;c<rgbChn.length;c++){
						if (fullResults[i][j].length>rgbChn[c]){
							returnFullResults[0][i][j][c]=fullResults[i][j][rgbChn[c]];
						} else returnFullResults[0][i][j][c]=null;
					}
				} else returnFullResults[0][i][j]=null;
			}
		}
		boolean showPSF=debugLevel>1;
		
		if (metrics==null) return null;
		
		if (showPSF || (focusMeasurementParameters.saveResults)){
			double [][][][] psfRGB=new double [psf_kernels.length][psf_kernels[0].length][][];
			// reorder and assign names to kernel color channels
			String [] rgbNames={"Red","Green","Blue"};
			for (int tileY=0;tileY< psf_kernels.length;tileY++) for (int tileX=0;tileX< psf_kernels[0].length;tileX++){
				if (psf_kernels[tileY][tileX]!=null){
					psfRGB[tileY][tileX]=new double [3][];
					for (int rgbi=0;rgbi<3;rgbi++) psfRGB[tileY][tileX][rgbi]=psf_kernels[tileY][tileX][rgbChn[rgbi]];
				} else psfRGB[tileY][tileX]=null;
			}
			ImageStack mergedStack=mergeKernelsToStack(psfRGB,rgbNames);
//			ImageStack mergedStack=mergeKernelsToStack(psf_kernels);
			if (mergedStack==null) {
				IJ.showMessage("Error","No PSF kernels- it is a bug");
				return null;
			}
			String psfTitle=((String) imp_sel.getProperty("timestamp")).replace('.','_')+".psf-tiff";
			if (focusMeasurementParameters.includeLensSerial && (focusMeasurementParameters.lensSerial.length()>0)){
//				psfTitle=String.format("LENS%S-",focusMeasurementParameters.lensSerial)+psfTitle;
				psfTitle=String.format("LENS%S-S%02d-",FOCUS_MEASUREMENT_PARAMETERS.lensSerial,FOCUS_MEASUREMENT_PARAMETERS.manufacturingState)+psfTitle;
			}
			
			ImagePlus imp_psf = new ImagePlus(psfTitle, mergedStack);
			imp_psf.getProcessor().resetMinAndMax();
			imp_psf.setProperty("comment",focusMeasurementParameters.comment);
			imp_psf.setProperty("timestamp", (String) imp_sel.getProperty("timestamp"));
			if (!Double.isNaN(focusMeasurementParameters.sensorTemperature))
				imp_psf.setProperty("sensorTemperature", ""+focusMeasurementParameters.sensorTemperature);
			imp_psf.setProperty("px0", ""+lensDistortionParameters.px0);
			imp_psf.setProperty("py0", ""+lensDistortionParameters.py0);
			imp_psf.setProperty("motor1", ""+focusMeasurementParameters.motorPos[0]);
			imp_psf.setProperty("motor2", ""+focusMeasurementParameters.motorPos[1]);
			imp_psf.setProperty("motor3", ""+focusMeasurementParameters.motorPos[2]);
			imp_psf.setProperty("lensSerial", ""+focusMeasurementParameters.lensSerial);
			//x0,y0,
			for (int ii=0;ii<sampleCoord.length;ii++){
				for (int jj=0;jj<sampleCoord[0].length;jj++){
					int index=ii*sampleCoord[0].length+jj;
					imp_psf.setProperty("pX_"+index, ""+sampleCoord[ii][jj][0]);
					imp_psf.setProperty("pY_"+index, ""+sampleCoord[ii][jj][1]);
					if (fullResults[ii][jj]!=null){
						for (int cc=0;cc<fullResults[ii][jj].length;cc++) if (fullResults[ii][jj][cc]!=null){
//							imp_psf.setProperty("R50_"+cc+"_"+index, ""+fullResults[ii][jj][cc][0]);
//							imp_psf.setProperty("RTratio_"+cc+"_"+index, ""+fullResults[ii][jj][cc][1]);
							imp_psf.setProperty("R50_x2_"+cc+"_"+index, ""+fullResults[ii][jj][cc][0]);
							imp_psf.setProperty("R50_y2_"+cc+"_"+index, ""+fullResults[ii][jj][cc][1]);
							imp_psf.setProperty("R50_xy_"+cc+"_"+index, ""+fullResults[ii][jj][cc][2]);
						}
					}
				}
			}
	    	(new JP46_Reader_camera(false)).encodeProperiesToInfo(imp_psf);
			if (showPSF) imp_psf.show();
			if (focusMeasurementParameters.saveResults) {
				String dir=getResultsPath(focusMeasurementParameters);
				File dFile=new File(dir);
				if (!dFile.isDirectory() &&  !dFile.mkdirs()) {
					String msg="Failed to create directory "+dir;
					IJ.showMessage(msg);
					throw new IllegalArgumentException (msg);
				}
				String path=dFile+Prefs.getFileSeparator()+psfTitle;
				if (debugLevel>0) System.out.println ("Saving annotated PSF kernels to "+path);
				if (updateStatus) IJ.showStatus("Saving annotated PSF kernels to "+path);
				(new FileSaver(imp_psf)).saveAsTiffStack(path);
			}
		}
		
		if (debugLevel>1) {
			for (int ii=0;ii<sampleCoord.length;ii++){
				for (int jj=0;jj<sampleCoord[0].length;jj++){
					System.out.println((ii*sampleCoord[0].length+jj+1)+": pX="+IJ.d2s(sampleCoord[ii][jj][0],1)+
							" pY="+IJ.d2s(sampleCoord[ii][jj][1],1)+
							" radial "+IJ.d2s(anglePerPixel[ii][jj][0]*180*60/Math.PI,3)+"'/pix"+
							" tangential "+IJ.d2s(anglePerPixel[ii][jj][1]*180*60/Math.PI,3)+"'/pix");
				}

			}
		}
		return metrics;
	}
	
	
	
	
	
	public double [][] extractPSFMetrics(
			double [][][] sampleCoord,
			double x0,
			double y0,
			double [][][][] psf,
			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
			double [] componentWeights, // for averaging metrics
			double [][][][] fullResults, // [y][x][color],should be null or match samples, last dimension {R50, BA} 
			int       debugLevel
	){
		if (fullResults!=null){
			for (int i=0;i<fullResults.length;i++) for (int j=0;j<fullResults[0].length;j++) fullResults[i][j]=null;
		}
		boolean [] squared={false,false,false,false,true,true,false,false}; // which metrics elements are equavalent to squared pixels
		double tiltScale=1000; //assymmetry per 1000 pixels
		int numMetrics=7;
		boolean [][] centerMask=focusMeasurementParameters.getCenterMask(
				x0,   // lens center on the sensor
    			y0,  // lens center on the sensor
    			focusMeasurementParameters.numInCenter);
		int numColors=0;
		for (int i=0;i<psf.length;i++) for (int j=0;j<psf[0].length;j++)if ((psf[i][j]!=null) && (psf[i][j].length>numColors))numColors=psf[i][j].length;
		double [][] metrics= new double[numColors+1][]; // includes weighted average as last line
		for (int color=0;color<numColors;color++){
			double S0=  0.0;
			double SX=  0.0;
			double SX2= 0.0;
			double SY=  0.0;
			double SY2= 0.0;
			double SXY= 0.0;
			double SF=  0.0;
			double SFX= 0.0;
			double SFY= 0.0;
			double SR=  0.0; // sum of PSF radiuses (or use something else) 
			double SA=  0.0; // sum of PSF areas (no PI - just square radiuses) 
			double SB=  0.0; // sum of modified PSF areas, attemprt to estimate "quality" 

			double S0center=  0.0;
			double SFcenter=  0.0;
			
			for (int i=0;i<psf.length;i++) for (int j=0;j<psf[0].length;j++){
				if ((psf[i][j]!=null) && (psf[i][j][color]!=null)){
					if ((fullResults!=null) && (fullResults[i][j]==null)){ //first non-null color for this sample
						fullResults[i][j]=new double [numColors][];
						for (int cc=0;cc<numColors;cc++) fullResults[i][j][cc]=null;
					}
//					if (fullResults!=null) fullResults[i][j][color]=new double[2];
					double x=(sampleCoord[i][j][0]-x0);
					double y=(sampleCoord[i][j][1]-y0);
					double r=Math.sqrt(x*x+y*y);			
					double ca=(sampleCoord[i][j][0]-x0)/r;
					double sa=(sampleCoord[i][j][1]-y0)/r;
//					System.out.println("extractPSFMetrics.. color="+color+" i="+i+" j="+j+" cos="+ca+" sin="+sa+
//							" psf_cutoffEnergy="+focusMeasurementParameters.psf_cutoffEnergy+
//							" psf_cutoffLevel="+focusMeasurementParameters.psf_cutoffLevel+
//							" psf_minArea="+focusMeasurementParameters.psf_minArea+
//							" psf_blurSigma="+focusMeasurementParameters.psf_blurSigma);
/*					
					double [] tanRad=		   matchSimulatedPattern.tangetRadialSizes(
							   ca, // cosine of the center to sample vector
							   sa, // sine of the center to sample vector
								psf[i][j][color],     // PSF function, square array, nominally positive
								focusMeasurementParameters.psf_cutoffEnergy,     // fraction of energy in the pixels to be used
								focusMeasurementParameters.psf_cutoffLevel,      // minimal level as a fraction of maximal
								focusMeasurementParameters.psf_minArea,      // minimal selected area in pixels
								focusMeasurementParameters.psf_blurSigma,    // optionally blur the selection
								0.1, //maskCutOff,
								debugLevel-2, // debug level
								i+":"+j+"["+color+"]"); //	   String        title) {    // prefix used for debug images
					if (tanRad==null){
						tanRad=new double[2];
						tanRad[0]=Double.NaN;
						tanRad[1]=Double.NaN;
					} else {
					
// debug calculations
						double [] x2y2xy= matchSimulatedPattern.x2y2xySizes(
									psf[i][j][color],     // PSF function, square array, nominally positive
									focusMeasurementParameters.psf_cutoffEnergy,     // fraction of energy in the pixels to be used
									focusMeasurementParameters.psf_cutoffLevel,      // minimal level as a fraction of maximal
									focusMeasurementParameters.psf_minArea,      // minimal selected area in pixels
									focusMeasurementParameters.psf_blurSigma,    // optionally blur the selection
									0.1, //maskCutOff,
									debugLevel-2, // debug level
									i+":"+j+"["+color+"]"); //	   String        title) {    // prefix used for debug images
						double [] tanRad1={
								Math.sqrt(sa*sa*x2y2xy[0]+ca*ca*x2y2xy[1]-2*ca*sa*x2y2xy[2]), 
								Math.sqrt(ca*ca*x2y2xy[0]+sa*sa*x2y2xy[1]+2*ca*sa*x2y2xy[2])};
						
						System.out.println("i="+i+" j="+j+" tanRad[0]="+tanRad[0]+ "("+tanRad1[0]+")"+
						" tanRad[1]="+tanRad[1]+ "("+tanRad1[1]+") ### "+x2y2xy[0]+", "+x2y2xy[1]+", "+x2y2xy[2]);
						
					}
*/
					double [] x2y2xy=null;
					try{
					x2y2xy= matchSimulatedPattern.x2y2xySizes(
							psf[i][j][color],     // PSF function, square array, nominally positive
							focusMeasurementParameters.psf_cutoffEnergy,     // fraction of energy in the pixels to be used
							focusMeasurementParameters.psf_cutoffLevel,      // minimal level as a fraction of maximal
							focusMeasurementParameters.psf_minArea,      // minimal selected area in pixels
							focusMeasurementParameters.psf_blurSigma,    // optionally blur the selection
							0.1, //maskCutOff,
							debugLevel-2, // debug level
							i+":"+j+"["+color+"]"); //	   String        title) {    // prefix used for debug images
					} catch (Exception e){
						System.out.println("Failed to get PSF size for i="+i+", j="+j);
					}
					if (x2y2xy==null){
						x2y2xy=new double[3];
						for (int ii=0;ii<3;ii++) x2y2xy[ii]=Double.NaN;
					} else {
						for (int ii=0;ii<x2y2xy.length;ii++) x2y2xy[ii]/=focusMeasurementParameters.subdiv*focusMeasurementParameters.subdiv/4;
					}
/*					
					double [] tanRad={
							Math.sqrt(sa*sa*x2y2xy[0]+ca*ca*x2y2xy[1]-2*ca*sa*x2y2xy[2]), 
							Math.sqrt(ca*ca*x2y2xy[0]+sa*sa*x2y2xy[1]+2*ca*sa*x2y2xy[2])};
*/
					// tanRad[0]/=(focusMeasurementParameters.subdiv/2); // to sesnor pixels
					// tanRad[1]/=(focusMeasurementParameters.subdiv/2);
					double [][]quadCoeff=null;
					try {
						quadCoeff=matchSimulatedPattern.approximatePSFQuadratic(
								psf[i][j][color],     // PSF function, square array, nominally positive
								focusMeasurementParameters.psf_cutoffEnergy,     // fraction of energy in the pixels to be used
								focusMeasurementParameters.psf_cutoffLevel,      // minimal level as a fraction of maximal
								focusMeasurementParameters.psf_minArea,      // minimal selected area in pixels
								focusMeasurementParameters.psf_blurSigma,    // optionally blur the selection
								0.1, //maskCutOff,
								debugLevel-2, // debug level
								i+":"+j+"["+color+"]"); //	   String        title) {    // prefix used for debug images
					} catch (Exception e) {
						System.out.println("Failed to get approximatePSFQuadratic(...) for i="+i+", j="+j);
						continue;
					}
/*
 *  f(x,y)=A*x^2+B*y^2+C*x*y+D*x+E*y+F
 *  
 *  Xc=(2*B*D-C*E)/(C^2-4*A*B)
 *  Yc=(2*A*E-C*D)/(C^2-4*A*B)
 *  top=A*Xc^2+B*Yc^2+C*Xc*Yc+D*Xc+E*Yc+F;
 *  x1=x-Xc
 *  y1=y-Yc
 *  f(x,y)=A*x1^2+B*y1^2+C*x1*y1+top
 *  
 *  Ellipse at half height:
 *  f(x,y)=top/2,
 *  A*x1^2+B*y1^2+C*x1*y1+top/2=0
 *  
 *  x2=cos(alpha)*x1+sin(alpha)*y1
 *  y2=-sin(alpha)*x1+cos(alpha)*y1
 *  
 *  x1=cos(alpha)*x2-sin(alpha)*y2
 *  y1=sin(alpha)*x2+cos(alpha)*y2
 *  
 *  A*(ca*x2-sa*y2)^2+B*(sa*x2+ca*y2)^2+C*(ca*x2-sa*y2)*(sa*x2+ca*y2)+top/2=0; // 50% height
 *  
   A*ca^2*x^2+ A*sa^2*y^2-2*A*ca*sa*x*y+
   B*sa^2*x^2+ B*ca^2*y^2+2*B*ca*sa*x*y+
   C*ca*sa*x^2-C*sa*ca*y^2+C*ca*ca*x*y-C*sa*sa*x*y=
   
   x^2*(A*ca^2+B*sa^2+C*ca*sa)+
   y^2*(A*sa^2+B*ca^2-C*ca*sa)+
//   x*y*((A+B)*ca*sa+C*(ca^2-sa^2))
   x*y*(2*(B-A)*ca*sa+C*(ca^2-sa^2))
   
   Ar=-2*(A*ca^2+B*sa^2+C*ca*sa)/top
   Br=-2*(A*sa^2+B*ca^2-C*ca*sa)/top
//   Cr=-((A+B)*ca*sa+C*(ca^2-sa^2))/top
   Cr=-2*(2*(B-A)*ca*sa+C*(ca^2-sa^2))/top
   Ar*Xr^2+Br*Yr^2+Cr*Xr*Yr=1
   
 *  x2=cos(beta)*x3-sin(beta)*y3
 *  y2=sin(beta)*x3+cos(beta)*y3
   
Ar*x2^2+Br*y2^2+Cr*x2*y2=1

area=2*pi/sqrt(Ar*Br-(Cr/2)^2) 
Reff=sqrt(1/sqrt(Ar*Br-(Cr/2)^2))
   
   
   
 */
					if (quadCoeff[0].length<6){
						System.out.println("extractPSFMetrics() - trying to approximate bad psf["+i+"]["+j+"]["+color+"], probably far out of focus or too low otf deconvolution invert parameter - that should not happen...");
						continue;
					}
					
					double xc=(2*quadCoeff[0][1]*quadCoeff[0][3]-quadCoeff[0][2]*quadCoeff[0][4])/(quadCoeff[0][2]*quadCoeff[0][2]-4*quadCoeff[0][0]*quadCoeff[0][1]);
/* 3944 was above line					
					Using linear approximation, M.det()=5.3717030823464884E-8 normMatix(mAarrayQ)=7.7524357668672154E17
					Exception in thread "Run$_AWT-EventQueue-0" java.lang.ArrayIndexOutOfBoundsException: 3
					        at Aberration_Calibration.extractPSFMetrics(Aberration_Calibration.java:3944)
					        at Aberration_Calibration.measurePSFMetrics(Aberration_Calibration.java:3826)
					        at Aberration_Calibration.moveAndMaybeProbe(Aberration_Calibration.java:3616)
					        at Aberration_Calibration.fineFocusingStepsAuto(Aberration_Calibration.java:3392)
					        at Aberration_Calibration.actionPerformed(Aberration_Calibration.java:2894)
*/					        
					 					
					double yc=(2*quadCoeff[0][0]*quadCoeff[0][4]-quadCoeff[0][2]*quadCoeff[0][3])/(quadCoeff[0][2]*quadCoeff[0][2]-4*quadCoeff[0][0]*quadCoeff[0][1]);

					double top=quadCoeff[0][0]*xc*xc+quadCoeff[0][1]*yc*yc+quadCoeff[0][2]*xc*yc+quadCoeff[0][3]*xc+quadCoeff[0][4]*yc+quadCoeff[0][5];
					double Ar= -2*(quadCoeff[0][0]*ca*ca + quadCoeff[0][1]*sa*sa + quadCoeff[0][2]*ca*sa)/top;
					double Br= -2*(quadCoeff[0][0]*sa*sa + quadCoeff[0][1]*ca*ca - quadCoeff[0][2]*ca*sa)/top;
					//						double Cr=((quadCoeff[0][0]+quadCoeff[0][1])*ca*sa+ quadCoeff[0][2]*(ca*ca-sa*sa))/top;					
					double Cr=-2*(2*(quadCoeff[0][1]-quadCoeff[0][1])*ca*sa+ quadCoeff[0][2]*(ca*ca-sa*sa))/top;
					double R50=Math.sqrt(1.0/Math.sqrt(Ar*Br-Cr*Cr/4));
//					double A50=Math.PI*R50*R50;
					double A50=R50*R50; // radius squared (to apply sqrt in the end)
					double BA=Br/Ar;
					double B50=(A50/2)*(BA*BA+1.0/(BA*BA));// for circles will be the same as A50, will be higher for long of the same area
					double w=1.0; // change later?
					double logBA=Math.log(Br/Ar);
					S0+=w;
					SX+=w*x;
					SX2+=w*x*x;
					SY+=w*y;
					SY2+=w*y*y;
					SXY+=w*x*y;
					SF+=w*logBA;
					SFX+=w*logBA*x;
					SFY+=w*logBA*y;
					SR+=R50;
					SA+=A50;
					SB+=B50;
					// calculate this only from the center samples
					if (centerMask[i][j]){
						S0center+=w;
						SFcenter+=R50; // only center samples
					}
					if (fullResults!=null){
						
//						fullResults[i][j][color][0]=Math.sqrt((tanRad[0]*tanRad[0]+tanRad[1]*tanRad[1])/2); // sensor pixels
//						fullResults[i][j][color][1]=tanRad[1]/tanRad[0]; // radial-to-tangential ratio
						fullResults[i][j][color]=x2y2xy.clone();
						
					}

					if (debugLevel>1){
						String debugStr="extractPSFMetrics: "+i+":"+j+"["+color+"]";
						for (int k=0;k<quadCoeff[0].length;k++) debugStr+=" "+quadCoeff[0][k];
						if (debugLevel>3)System.out.println(debugStr);
						System.out.println(i+":"+j+"["+color+"] x="+IJ.d2s(x,1)+" y="+IJ.d2s(y,1)+" xc="+IJ.d2s(xc,2)+" yc="+IJ.d2s(yc,2)+
//								" tan="+IJ.d2s(tanRad[0],3)+" rad="+IJ.d2s(tanRad[1],3)+
								"x2="+IJ.d2s(x2y2xy[0],3)+"y2="+IJ.d2s(x2y2xy[1],3)+"xy="+IJ.d2s(x2y2xy[2],3)+
								" Ar="+IJ.d2s(Ar,3)+" Br="+IJ.d2s(Br,3)+" Cr="+IJ.d2s(Cr,3)+ " **** Br/Ar="+IJ.d2s(Br/Ar,3)+" R50="+IJ.d2s(R50,3)+" subpixels" +" A50="+IJ.d2s(Math.sqrt(A50),2)+" subpixels" +" B50=" + IJ.d2s(Math.sqrt(B50),3)+" subpixels");
					}
				}
			}
			if (S0>0){
				metrics[color]=new double [numMetrics+1]; // last - weights of components
				metrics[color][3]=SR/S0/(focusMeasurementParameters.subdiv/2);
				metrics[color][4]=Math.sqrt(SA/S0)/(focusMeasurementParameters.subdiv/2);
				metrics[color][5]=Math.sqrt(SB/S0)/(focusMeasurementParameters.subdiv/2);
				if (S0center>0.0){
					metrics[color][6]=SFcenter/S0center/(focusMeasurementParameters.subdiv/2);
				} else {
					metrics[color][6]=0.0;
				}
				double [][] aM={
						{SX2,SXY,SX},
						{SXY,SY2,SY},
						{SX, SY, S0}};
				double [][] aB={{SFX},{SFY},{SF}};
				Matrix M=new Matrix(aM);
				Matrix B=new Matrix(aB);
				Matrix V=M.solve(B);
				if (debugLevel>3) {
					M.print(9,5);
					B.print(9,5);
					V.print(9,5);
				}
				metrics[color][0]=V.get(2,0);
				metrics[color][1]=V.get(0,0)*tiltScale;
				metrics[color][2]=V.get(1,0)*tiltScale;
			} else 	metrics[color]=null;
		}
		//componentWeights
		metrics[numColors]=new double [numMetrics+1];
		for (int n=0;n<numMetrics;n++) metrics[numColors][n]=0.0;
		double metricsSumWeight=0.0;
    	for (int c=0;c<numColors;c++) if (metrics[c]!=null){
    		for (int n=0;n<numMetrics;n++) if (squared[n]){
    			metrics[numColors][n]+=metrics[c][n]*metrics[c][n]*componentWeights[c];
    		} else metrics[numColors][n]+=metrics[c][n]*componentWeights[c];
    		metricsSumWeight+=componentWeights[c];
    	}
    	if (metricsSumWeight==0) {
    		System.out.println("extractPSFMetrics(): not a single color available for this PSF!");
    		return null;  
    	}
       	for (int c=0;c<numColors;c++) if (metrics[c]!=null){
       		metrics[c][numMetrics]=componentWeights[c]/metricsSumWeight;
       	}
		for (int n=0;n<numMetrics;n++) {
			metrics[numColors][n]/=metricsSumWeight;
			if (squared[n]) metrics[numColors][n]=Math.sqrt(metrics[numColors][n]);
		}
		metrics[numColors][numMetrics]=1.0;
		return metrics;
	}
	
	
	
	//returns number of matched laser pointers
/*	
	public int updateFocusGrid(
			double x0,   // lens center on the sensor
			double y0,  // lens center on the sensor
			ImagePlus imp,
			MatchSimulatedPattern matchSimulatedPattern,
			MatchSimulatedPattern.DistortionParameters distortionParametersDefault,
			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			MatchSimulatedPattern.LaserPointer laserPointer, // null OK
			SimulationPattern.SimulParameters  simulParametersDefault,
			boolean maskNonPSF, // mask out areas not needed for focusing PSF measurements
			boolean equalizeGreens,
			int       threadsMax,
			boolean   updateStatus,
			int debug_level){// debug level used inside loops
		long 	  startTime=System.nanoTime();

		SimulationPattern.SimulParameters  simulParameters=simulParametersDefault.clone();
		simulParameters.smallestSubPix=             focusMeasurementParameters.smallestSubPix;
		simulParameters.bitmapNonuniforityThreshold=focusMeasurementParameters.bitmapNonuniforityThreshold;
		simulParameters.subdiv=                     focusMeasurementParameters.subdiv;

		MatchSimulatedPattern.DistortionParameters distortionParameters=  distortionParametersDefault.clone();
		distortionParameters.refineInPlace=false;

		distortionParameters.correlationMaxOffset=focusMeasurementParameters.maxCorr;

		distortionParameters.correlationSize=focusMeasurementParameters.correlationSize;
		distortionParameters.correlationGaussWidth=focusMeasurementParameters.correlationGaussWidth;
		distortionParameters.refineCorrelations=false;
		distortionParameters.fastCorrelationOnFirstPass=true;
		distortionParameters.fastCorrelationOnFinalPass=true;

		distortionParameters.correlationAverageOnRefine=false;
		distortionParameters.minUVSpan=focusMeasurementParameters.minUVSpan;
		distortionParameters.flatFieldCorrection=focusMeasurementParameters.flatFieldCorrection;
		distortionParameters.flatFieldExpand=focusMeasurementParameters.flatFieldExpand;
		
		if (maskNonPSF) {
			distortionParameters.numberExtrapolated=0; //1; //3; // measuring PSF - extrapolate
		} else {
			distortionParameters.numberExtrapolated=1; // measuring distortions - do not extrapolate
		}
		
		
//System.out.println("distortionParameters.correlationSize="+distortionParameters.correlationSize);
		// add more overwrites
		boolean updating=(matchSimulatedPattern.PATTERN_GRID!=null) &&
		(!distortionParameters.flatFieldCorrection || (matchSimulatedPattern.flatFieldForGrid!=null));

		ImagePlus imp_eq=matchSimulatedPattern.applyFlatField (imp); // will throw if image size mismatch
		if (updating) {
			double maxActualCorr= matchSimulatedPattern.refineDistortionCorrelation (
					distortionParameters, //
					patternDetectParameters,
					simulParameters,
					equalizeGreens,
					imp_eq, // if grid is flat-field calibrated, apply it
					focusMeasurementParameters.maxCorr, // maximal allowed correction, in pixels (0.0) - any
					threadsMax,
					updateStatus,
					debug_level);// debug level used inside loops
			matchSimulatedPattern.recalculateWaveVectors (
					updateStatus,
					debug_level);// debug level used inside loops

			if (debug_level>1) System.out.println("refineDistortionCorrelation() finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));

			//Bad:	NaN (no cells),	maxActualCorr<0 -number of new empty nodes, >focusMeasurementParameters.maxCorr - also bad (no correction)
			if (!((maxActualCorr>=0) && (maxActualCorr<=focusMeasurementParameters.maxCorr))) {
				if (debug_level>0) System.out.println("updateFocusGrid() failed, refineDistortionCorrelation() ->"+maxActualCorr+ " (maxCorr="+focusMeasurementParameters.maxCorr+")");
				// Do full 			


				updating=false;
			} else {
				if (debug_level>1) System.out.println("updateFocusGrid() ->"+maxActualCorr+ " (maxCorr="+focusMeasurementParameters.maxCorr+")");
			}
		}
		int numAbsolutePoints=0;
		if (updating) {		
			// add new nodes if the appeared after shift of the pattern
			if (debug_level>1) { // calulate/print number of defined nodes in the grid
				System.out.println("updateFocusGrid(), number of defined grid cells (before distortions()) = "+matchSimulatedPattern.numDefinedCells());
			}
			matchSimulatedPattern.distortions(
					distortionParameters, //
					patternDetectParameters,
					simulParameters,
					equalizeGreens,
					imp_eq, // image to process
					threadsMax,
					updateStatus,
					debug_level);// debug level used inside loops
			if (debug_level>1) System.out.println("distortions() finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));

			matchSimulatedPattern.recalculateWaveVectors (
					updateStatus,
					debug_level);// debug level used inside loops
			if (debug_level>1) System.out.println("recalculateWaveVectors() finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
			if (debug_level>1) { // calulate/print number of defined nodes ina grid
				System.out.println("updateFocusGrid(), number of defined grid cells (after distortions()) = "+matchSimulatedPattern.numDefinedCells());
			}
		} else {
//			   matchSimulatedPattern.invalidateFlatFieldForGrid(); //Keep these!
//			   matchSimulatedPattern.invalidateFocusMask();
			numAbsolutePoints=matchSimulatedPattern.calculateDistortions( // allow more of grid around pointers?
					distortionParameters, //
					patternDetectParameters,
					simulParameters,
					equalizeGreens,
					imp_eq,
					laserPointer, //null, //LASER_POINTERS, // LaserPointer laserPointer, // LaserPointer object or null
					true, // don't care -removeOutOfGridPointers
					threadsMax,
					updateStatus,
					DEBUG_LEVEL,
					distortionParameters.loop_debug_level); // debug level

		}
		if (maskNonPSF) {
			matchSimulatedPattern.maskFocus(
	    			x0,   // lens center on the sensor
	    			y0,  // lens center on the sensor
					focusMeasurementParameters);
			if (debug_level>1) System.out.println("maskFocus() finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
			if (debug_level>1) { // calulate/print number of defined nodes ina grid
				System.out.println("number of defined grid cells (after maskFocus()) = "+matchSimulatedPattern.numDefinedCells());
			}
			if(debug_level>2) {
			       double [] test_masked=new double [matchSimulatedPattern.focusMask.length];
			       float [] pixels_eq=(float []) imp_eq.getProcessor().getPixels();
				   for (int i=0;i<test_masked.length;i++) test_masked[i]=matchSimulatedPattern.focusMask[i]?pixels_eq[i]:0.0;
				   SDFA_INSTANCE.showArrays(test_masked,matchSimulatedPattern.getImageWidth(), matchSimulatedPattern.getImageHeight(), "MASKED");
			}

		}
		if (debug_level>1) System.out.println("updateFocusGrid() finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		if (debug_level>2) {
	       double [] test_uv=new double [matchSimulatedPattern.UV_INDEX.length];
		   for (int i=0;i<matchSimulatedPattern.UV_INDEX.length;i++) test_uv[i]=matchSimulatedPattern.UV_INDEX[i];
		   SDFA_INSTANCE.showArrays(test_uv,matchSimulatedPattern.getImageWidth(), matchSimulatedPattern.getImageHeight(), "UV_INDEX");
		}
		return numAbsolutePoints;

	}
*/
	//====================================================
	public double [][][] pixToAngles(
			double [][][] pXYArray,
			double pX0, // pixel coordinate of the the optical center
			double pY0, // pixel coordinate of the the optical center
			double pixelSize, // in microns
			double focalLength, // in mm
			double r0, // in mm - half sensor width, for wich distortion is calculated
			double a4, // polynomial coefficients i Rpix/Rph=a4*r^4+a3*r^3+a2*r^2+a1*r+ (1-a4-a3-a2-a1); r=Rph/r0
			double a3,
			double a2,
			double a1) {
		double [][][] p2A=new double[pXYArray.length][pXYArray[0].length][2];
		for (int i=0;i<p2A.length;i++) for (int j=0;j<p2A[0].length;j++){
			p2A[i][j]=pixToAngles(
					pXYArray[i][j][0], // pixel X for which the data is requested
					pXYArray[i][j][1], // pixel Y for which the data is requested
					pX0, // pixel coordinate of the the optical center
					pY0, // pixel coordinate of the the optical center
					pixelSize, // in microns
					focalLength, // in mm
					r0, // in mm - half sensor width, for wich distortion is calculated
					a4, // polynomial coefficients i Rpix/Rph=a4*r^4+a3*r^3+a2*r^2+a1*r+ (1-a4-a3-a2-a1); r=Rph/r0
					a3,
					a2,
					a1);

		}
		return p2A;
	}
	/**
	 *  Calculates angle (radians) per sensor pixel in radial and tangential directions
	 */
	public double [] pixToAngles(
			double pX, // pixel X for which the data is requested
			double pY, // pixel Y for which the data is requested
			double pX0, // pixel coordinate of the the optical center
			double pY0, // pixel coordinate of the the optical center
			double pixelSize, // in microns
			double focalLength, // in mm
			double r0, // in mm - half sensor width, for wich distortion is calculated
			double a4, // polynomial coefficients i Rpix/Rph=a4*r^4+a3*r^3+a2*r^2+a1*r+ (1-a4-a3-a2-a1); r=Rph/r0
			double a3,
			double a2,
			double a1) {
		double rPix=0.001*pixelSize*Math.sqrt((pX-pX0)*(pX-pX0)+(pY-pY0)*(pY-pY0)); // distance in mm
		double rPinHole=rPix; // distance for non-distorted (pin hole model), calculate in iterations
		double iterationThreshold=0.0001; // relative (to r0) difference between approximation and rPix to end iterations
		int maxNumIter=100;
		double rPixTorPinHole=1.0; // rPix/rPinHole ratio, should not be NaN when both are 0.0
		double drPixdrPh=1.0; // derivative drPix/drPinHole as a function of rPinHole
		for (int numIter=0;numIter<maxNumIter;numIter++){
			// calculate approximated rPix from rPinHole
			double r=rPinHole/r0;
			double r2=r*r;
			double r3=r2*r; 
			double r4=r3*r; 
			rPixTorPinHole=a4*r4+a3*r3+a2*r2+a1*r+(1.0-a4-a3-a2-a1); // rPix/rPinHole
			double rPixApprox=rPinHole*rPixTorPinHole;
			drPixdrPh=(1/r0)* (4*a4*r3+3*a3*r2+2*a2*+a1)* rPinHole  + rPixTorPinHole;
			if (DEBUG_LEVEL>1) System.out.println("--pX="+IJ.d2s(pX,1)+" pY="+IJ.d2s(pY,1)+" rPix="+IJ.d2s(rPix,4)+", rPinHole="+IJ.d2s(rPinHole,4));
			if (DEBUG_LEVEL>1) System.out.println("r0="+IJ.d2s(r0,3)+" a4="+IJ.d2s(a4,4)+" a34="+IJ.d2s(a3,4)+" a2="+IJ.d2s(a2,4)+" a1="+IJ.d2s(a1,4));
			if (Math.abs(rPixApprox-rPix)/r0<iterationThreshold) break;
			rPinHole+=(rPix-rPixApprox)/drPixdrPh;
		}
		
		double cosAlpha=focalLength/Math.sqrt(focalLength*focalLength+rPinHole*rPinHole);
		if (DEBUG_LEVEL>1) System.out.println("pX="+IJ.d2s(pX,1)+" pY="+IJ.d2s(pY,1)+" rPix="+IJ.d2s(rPix,4)+
				", rPinHole="+IJ.d2s(rPinHole,4)+", cos(alpha)="+IJ.d2s(cosAlpha,4)+", drPixdrPh="+IJ.d2s(drPixdrPh,4)+
				", rPixTorPinHole="+IJ.d2s(rPixTorPinHole,4));
        double [] result={
				0.001*pixelSize*cosAlpha*cosAlpha/focalLength/drPixdrPh,
				0.001*pixelSize*cosAlpha/focalLength/rPixTorPinHole};
		
		return result;
	}
	
	public double [][][][] focusPSF (
			double x0,   // lens center on the sensor
			double y0,   // lens center on the sensor
			ImagePlus imp,
			MatchSimulatedPattern matchSimulatedPattern,
			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
			SimulationPattern.SimulParameters  simulParametersDefault,
			double overexposedMaxFraction,
			EyesisAberrations.ColorComponents colorComponents,
//			int           PSF_subpixel, 
			EyesisAberrations.OTFFilterParameters otfFilterParametersDefault,
			EyesisAberrations.PSFParameters psfParametersDefault,
			int       threadsMax,
			boolean   updateStatus,
			int debug_level){// debug level used inside loops
		long 	  startTime=System.nanoTime();
		
// combine with 			matchSimulatedPattern.maskFocus(FOCUS_MEASUREMENT_PARAMETERS); ?
		SimulationPattern.SimulParameters  simulParameters=simulParametersDefault.clone();

		simulParameters.smallestSubPix=             focusMeasurementParameters.smallestSubPix;
		simulParameters.bitmapNonuniforityThreshold=focusMeasurementParameters.bitmapNonuniforityThreshold;
		simulParameters.subdiv=                     focusMeasurementParameters.subdiv;
		
		EyesisAberrations.PSFParameters psfParameters=  psfParametersDefault.clone();
		psfParameters.approximateGrid=focusMeasurementParameters.approximateGrid;
		psfParameters.centerPSF=focusMeasurementParameters.centerPSF;

		psfParameters.mask1_sigma=focusMeasurementParameters.mask1_sigma;
		psfParameters.mask1_threshold=focusMeasurementParameters.mask1_threshold;
		psfParameters.gaps_sigma=focusMeasurementParameters.gaps_sigma;
		psfParameters.mask_denoise=focusMeasurementParameters.mask_denoise;
		
		EyesisAberrations.OTFFilterParameters otfFilterParameters=otfFilterParametersDefault.clone();
		otfFilterParameters.deconvInvert=focusMeasurementParameters.deconvInvert;
		if (!psfParameters.approximateGrid) {
/*			SIM_ARRAY=	simulateGridAll (
				matchSimulatedPattern.getImageWidth(),
				matchSimulatedPattern.getImageHeight(),
				matchSimulatedPattern.getDArray(),
				2, // gridFrac, // number of grid steps per pattern full period
				simulParameters,  // TODO: copy some parameters from focusMeasurementParameters (faster, not so precise)??
				threadsMax,
				updateStatus,
				debug_level); // debug level*/
			SIM_ARRAY=	(new SimulationPattern(simulParameters)).simulateGridAll (
					matchSimulatedPattern.getImageWidth(),
					matchSimulatedPattern.getImageHeight(),
					matchSimulatedPattern,
					2, // gridFrac, // number of grid steps per pattern full period
					simulParameters,  // TODO: copy some parameters from focusMeasurementParameters (faster, not so precise)??
					threadsMax,
					updateStatus,
					DEBUG_LEVEL,
					debug_level); // debug level
		}
		int [][][] sampleMap=new int [focusMeasurementParameters.numSamples[1]][focusMeasurementParameters.numSamples[0]][2];
		int halfSize=focusMeasurementParameters.sampleSize/2;
    	double [][][] sampleCoord= focusMeasurementParameters.sampleCoordinates(
    			x0,   // lens center on the sensor
    			y0);  // lens center on the sensor

		for (int i=0;i<focusMeasurementParameters.numSamples[1];i++){
			for (int j=0;j<focusMeasurementParameters.numSamples[0];j++){
				sampleMap[i][j][0]=(int) sampleCoord[i][j][0]- halfSize;
				sampleMap[i][j][1]=(int) sampleCoord[i][j][1]- halfSize;
			}
		}
		double [][][][] psf_kernels=createPSFMap(
				matchSimulatedPattern,
				matchSimulatedPattern.applyFlatField (imp), // if grid is flat-field calibrated, apply it
				sampleMap, //final int [][][]        sampleList, // optional (or null) 2-d array: list of coordinate pairs (2d - to match existent  PSF_KERNEL_MAP structure)  
				overexposedMaxFraction,
				simulParameters, // simulation parameters
				0, // MAP_FFT_SIZE, // scanImageForPatterns:FFT size
				PATTERN_DETECT, // needed for testing! null, //PATTERN_DETECT, // pattern detection parameters
				0, //FFT_OVERLAP,
				focusMeasurementParameters.sampleSize/2, // FFT_SIZE,  
				colorComponents,
				focusMeasurementParameters.subdiv, //PSF_subpixel, // maximal iterations when looking for local maximum
				otfFilterParameters,
				psfParameters, // step of the new map (should be multiple of map step)
				focusMeasurementParameters.minDefinedArea, /******************************************* CHANGE */
				focusMeasurementParameters.PSFKernelSize, // size of square used in the new map (should be multiple of map step)
				threadsMax,
				updateStatus,
				debug_level);// debug level used inside loops
		if (debug_level>1) System.out.println("focusPSF() finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		return psf_kernels;
		
	}
	
	
	
	
	public void calcLaser(ImagePlus imp_pointers, boolean noMessageBoxes){
		double [][] pointersXY=new double [imp_pointers.getStackSize()-1][];
		int dbgLev=1;
		double[][] backgroundBayer =splitBayer (imp_pointers,  1, null, true); // full window
		double[][] pointedBayer;
		for (int pn=0;pn<pointersXY.length;pn++) {
			pointedBayer=splitBayer (imp_pointers,  pn+2, null, true);
			pointersXY[pn]=LASER_POINTERS.getPointerXY( // returns x,y pair or null if pointer not detected
					backgroundBayer,        // Bayer array of the background (lasers off) image 0,3 -G, 1-R,2-B
					pointedBayer,           // Bayer array of the background (laser on) image
					imp_pointers.getWidth(),// image width in pixels
	    			(pn==0),                //boolean modBackground,       // modify background array (on the first pass)
	    			"laser"+pn,             // String title,
					DEBUG_LEVEL             // debug level (normal == 1)
			);
		}
		if ((matchSimulatedPattern!=null) && (matchSimulatedPattern.patternOK())) {
			dbgLev=matchSimulatedPattern.debugLevel;
			matchSimulatedPattern.debugLevel=DEBUG_LEVEL;
//			pointersUV=matchSimulatedPattern.uvFromXY(pointersXY);
			int acalibrated=matchSimulatedPattern.calibrateGrid(LASER_POINTERS,
					pointersXY,
					true,
					-1, // hinted rotation undefined
					null, // hinted translation undefined
					noMessageBoxes,
					DEBUG_LEVEL); // remove out-of-grid pointers
			matchSimulatedPattern.debugLevel=dbgLev;
			if (DEBUG_LEVEL>0) {
				System.out.println("matchSimulatedPattern.laserCalibrateGrid() returned "+acalibrated+
				((acalibrated>0)?" laser points used.":" - error code"));
			}
		}
		
	}
	
/* ======================================================================== */
//Flat-field related
    public void processFlatField (FlatFieldParameters flatFieldParameters, float [] fpixels, int origWidth, String title){
        int numSect=5; // assuming 5 sections
    	int origHeight=fpixels.length/origWidth; 
    	int width= (origWidth-1)/flatFieldParameters.decimate+1;
    	int height= ((fpixels.length/origWidth)-1)/flatFieldParameters.decimate+1;
    	double[] dpixels= new double[width*height];
    	for (int i=0;i<dpixels.length;i++)dpixels[i]=0.0;
    	double scale=1.0/flatFieldParameters.decimate/flatFieldParameters.decimate;
    	for (int i=0;i<fpixels.length;i++) dpixels[((i/origWidth/flatFieldParameters.decimate)*width)+((i%origWidth)/flatFieldParameters.decimate)]+=scale*Math.log(fpixels[i]);
    	int [] dmargins=new int [flatFieldParameters.margins.length];
    	for (int i=0;i<dmargins.length;i++) dmargins[i]=flatFieldParameters.margins[i]/flatFieldParameters.decimate;
    	
    	if (DEBUG_LEVEL>1) SDFA_INSTANCE.showArrays(dpixels,width,height, title+"-downsampled");
    	double tilt=calculateFlatFieldTilt(flatFieldParameters, dpixels, width, title);
	    if (DEBUG_LEVEL>1) System.out.println("Optimal tilt="+tilt);
    	
    	double [] coeff=null;
//    	if      (flatFieldParameters.functionType==0) coeff=estimateVignettingPoly(flatFieldParameters,dpixels, width, tilt, title);
//    	else if (flatFieldParameters.functionType==1) coeff=estimateVignettingExp(flatFieldParameters,dpixels, width, tilt, title);
    	coeff=estimateVignetting(flatFieldParameters,dpixels, width, tilt, title);
    	if (coeff==null) return;
        if (DEBUG_LEVEL>1) {
        	double [] diffPixels=buildFFModel(flatFieldParameters,numSect,width,height,coeff);
        	for (int i=0;i<dpixels.length;i++) diffPixels[i]=dpixels[i]-diffPixels[i];
        	SDFA_INSTANCE.showArrays(diffPixels,width,height, title+"-diff"+coeff[3]+"_"+coeff[4]+"_"+coeff[5]+"_"+coeff[6]);
        }
        double [] coeff0=coeff.clone();
        for (int i=0;i<numSect;i++) coeff0[i]=0.0;
        if  (flatFieldParameters.functionType==0) {
//            coeff0[10]*=flatFieldParameters.decimate;
//            coeff0[11]*=flatFieldParameters.decimate;
            coeff0[ 9]*=flatFieldParameters.decimate;
            coeff0[10]*=flatFieldParameters.decimate;
        } else {
            coeff0[7]*=flatFieldParameters.decimate;
            coeff0[8]*=flatFieldParameters.decimate;
        }
        double [] dVignetting=buildFFModel(flatFieldParameters,numSect,origWidth,origHeight,coeff0);
        for (int i=0;i<dVignetting.length;i++) dVignetting[i]=Math.exp(dVignetting[i]);
    	SDFA_INSTANCE.showArrays(dVignetting ,origWidth, origHeight, title+"-vignetting");
   	
    }
// from real center (not from clear), 1.- <-> half full height
    private double [][] calcFF1DSections(FlatFieldParameters flatFieldParameters,
    		                             double [] dpixels,
    		                             int width,
    		                             double tilt,
    		                             String title,
    		                             double [] y01234) { // should be initialized to double[3], will be set here
    	int numSect=y01234.length;
    	int height= dpixels.length/width;
    	int [] dmargins=new int [flatFieldParameters.margins.length];
    	for (int i=0;i<dmargins.length;i++) dmargins[i]=flatFieldParameters.margins[i]/flatFieldParameters.decimate;
    	int dSampleWidth=flatFieldParameters.sampleWidth/flatFieldParameters.decimate;
    	if (dSampleWidth<1) dSampleWidth=1;
    	double max=0.0;
    	double min=0.0;
    	double c0=1.0/dSampleWidth;
    	int clearHeight=height-dmargins[2]-dmargins[3];
    	double [] tilts=new double[numSect];
    	for (int n=0;n<numSect;n++) tilts[n]=tilt;
    	if (flatFieldParameters.noTiltEdges) {
    		tilts[1]=0.0;
    		tilts[2]=0.0;
    	}
    	double [][] sections=new double [numSect][clearHeight];
    	int [] j0=new int [numSect];
    	int [] jt=new int [numSect];
    	jt[0]=width/2;
    	jt[1]=dmargins[0];
    	jt[2]=width-1-dmargins[1];
    	if (numSect>3){
    		if (flatFieldParameters.noTiltEdges) {
      		  jt[3]=(int) (jt[0]+flatFieldParameters.section34*(jt[1]-jt[0]) - (height/2-dmargins[2])*Math.abs(tilt));
    		  jt[4]=(int) (jt[0]+flatFieldParameters.section34*(jt[2]-jt[0]) + (height/2-dmargins[2])*Math.abs(tilt));
    		  if (jt[3]<(jt[1]+dSampleWidth)) jt[3]=jt[1]+dSampleWidth;
    		  if (jt[4]>(jt[2]-dSampleWidth)) jt[4]=jt[2]-dSampleWidth;
    		} else {
    		  jt[3]=(int) (jt[0]+flatFieldParameters.section34*(jt[1]-jt[0]));
    		  jt[4]=(int) (jt[0]+flatFieldParameters.section34*(jt[2]-jt[0]));
    		}  
    	}
    	
    	j0[0]=(int) ((width+(height-2*dmargins[2])*tilts[0])/2); // center band
    	if (tilt>0){
    		j0[1]= jt[1]+(int) (clearHeight*tilts[1]);
    		j0[2]= jt[2];
        	if (numSect>3){
        		j0[3]=jt[3]+(int) (clearHeight*tilts[3]);
        		j0[4]=jt[4];
        	}
    	
    	} else {
    		j0[1]= jt[1];
    		j0[2]= jt[2]+(int) (clearHeight*tilts[2]);
        	if (numSect>3){
        		j0[3]=jt[3];
        		j0[4]=jt[4]+(int) (clearHeight*tilts[4]);
        	}
    	}
    	for (int n=0;n<numSect;n++) {
    		y01234[n]=j0[n]+0.5*dSampleWidth-(height/2-dmargins[2])*tilts[n];
    	}
    	for (int i=0;i<clearHeight;i++) {
    		int ia=i+dmargins[2];
    		for (int n=0;n<numSect;n++) sections[n][i]=0.0;
    		for (int j=0;j<dSampleWidth;j++) {
        		for (int n=0;n<numSect;n++) {
        			int ja=j- (int) (i*tilts[n]); // actual j after tilt
        			sections[n][i]+=c0*dpixels[width*ia+ (j0[n]+ja)];
        		}
    		}
       		for (int n=0;n<3;n++) {
         		  if ((i==0) || (sections[n][i]<min)) min=sections[n][i];
         		  if ((i==0) || (sections[n][i]>max)) max=sections[n][i];
     		}
    	}
	    if ((DEBUG_LEVEL>1) && (title!=null)) {
	    	double [] dp=dpixels.clone();
	    	double bump=0.1*(max-min);
	    	for (int i=0;i<clearHeight;i++) {
	    		int ia=i+dmargins[2];
	    		for (int j=0;j<dSampleWidth;j++) {
	        		for (int n=0;n<numSect;n++) {
		    			int ja=j- (int) (i*tilts[n]); // actual j after tilt
	        			dp[width*ia+ (j0[n]+ja)]+=bump;
	        		}
	    		}
	    	}
        	SDFA_INSTANCE.showArrays(dp,width,height, title+"-selections");
	    	
	    }
        return sections;   	
    }
    
    private double [] calcFF1DWeights(FlatFieldParameters flatFieldParameters,
    		                          double [][] sections, // data rows (vertical, tilted, CW - positive) - just dimensions used
    		                          int margin,           // margin for the beginning of each data row
    		                          int height,           // image height (radius==1.0 is half height)
    		                          double [] y01234,     // horizontal positions of the center (half height!) of each row
    		                          double centerHor,     // hor. location of the center (for weight calculations) from top left corner of full image (no margins) 
    		                          double centerVert,    // same for vertical.
    		                          double tilt) {        // rows tilt (hor pixels per pixel in vertical direction, tilted CW - positive, CCW - negative
    	int numSect=y01234.length;
    	double [] tilts=new double[numSect];
    	for (int n=0;n<numSect;n++) tilts[n]=tilt;
    	if (flatFieldParameters.noTiltEdges) {
    		tilts[1]=0.0;
    		tilts[2]=0.0;
    	}
    	double [] weights=new double [sections.length*sections[0].length];
    	double x; //vertical, down
    	double y; //horizontal, right
    	double rr02=4.0/(height*height);
//    	double [] tilts={tilt,(flatFieldParameters.noTiltEdges)?0.0:tilt,(flatFieldParameters.noTiltEdges)?0.0:tilt};
    	
    	if (flatFieldParameters.centerWeight<0) {
    		for (int i=0; i<weights.length;i++) weights[i]=1.0;
    		return weights;
    	}
    	for (int n=0;n<sections.length;n++) for (int i=0;i<sections[n].length;i++) {
    		x=i+margin-centerVert;
    		y=y01234[n] - tilts[n]* (i+margin -0.5*height);
    		weights[n*sections[n].length+i]=flatFieldParameters.centerWeight+rr02*(x*x+y*y);
    	}
        return weights;   	
    }

    private double [] estimateVignetting(FlatFieldParameters flatFieldParameters, double [] dpixels, int width, double tilt, String title) {
    	int height= dpixels.length/width;
	    String [] coefficientNamesExp={"a00", // shift for the center tilted from vertical band
	    		                    "a01", // shift for the left tilted from vertical band
	    		                    "a02", // shift for the right tilted from vertical band
	    		                    "a03", // shift for the middle-left tilted from vertical band
	    		                    "a04", // shift for the middle-right tilted from vertical band
	    		                     "a1", // 1-a1*r2^a2
	    		                     "a2",
	    		                     "x0",
	    		                     "y0"};
	    String [] coefficientNamesPoly={
	    		"a00", // shift for the center tilted from vertical band
                "a01", // shift for the left tilted from vertical band
                "a02", // shift for the rightr tilted from vertical band
                "a03", // shift for the middle-left tilted from vertical band
                "a04", // shift for the middle-right tilted from vertical band
                 "a1",
                 "a2",
                 "a3",
                 "a4",
//                 "a5",
                 "x0",
                 "y0"};
// make them programmable strategies (list of numbers?)	    
    	boolean [][] masksExp=
    	       {{true,true,true,true,true,true,true,false,false},
    			{true,true,true,true,true,true,true,true,true}};
//    	boolean [] maskPoly={true,true,true,true,true,true,true,true,false,false,false,false};
    	boolean [][] masksPoly={
    			{true,true,true,true,true,false,false,false,true,true,true},
//    			{true,true,true,true,true,true, false,false,true,false,false}};
		{true,true,true,true,true,
    				(flatFieldParameters.functionModifier == 0), //false, // r^2
    				(flatFieldParameters.functionModifier == 1), //true, // r^4
    				(flatFieldParameters.functionModifier == 2), //false, // r^6
    				true, // r^8
    				false,
    				false}};
        if (flatFieldParameters.functionModifier == 1){
        	
        }
	    String [] coefficientNames=(flatFieldParameters.functionType==1)?coefficientNamesExp:coefficientNamesPoly;
	    boolean [][] bmasks=           (flatFieldParameters.functionType==1)?masksExp:masksPoly;

    	int [] dmargins=new int [flatFieldParameters.margins.length];
    	for (int i=0;i<dmargins.length;i++) dmargins[i]=flatFieldParameters.margins[i]/flatFieldParameters.decimate;
    	int dSampleWidth=flatFieldParameters.sampleWidth/flatFieldParameters.decimate;
    	if (dSampleWidth<1) dSampleWidth=1;
    	int clearHeight=height-dmargins[2]-dmargins[3];
//    	double [] y012=new double[3];
    	double [] y01234=new double[5];
    	double [][] sections=calcFF1DSections(flatFieldParameters, dpixels, width, tilt, title,y01234);
    	double max=sections[0][0];
    	for (int i=0;i<sections.length;i++) for (int j=0;j<sections[i].length;j++) if (max<sections[i][j]) max=sections[i][j]; 
    	
    	if (DEBUG_LEVEL>1) {
		  for (int l=0;l<y01234.length;l++) {
			  System.out.println(">>>>   y01234["+l+"]="+y01234[l]);
			
		  }
    	}
    	
    	double [] errorWeight= calcFF1DWeights(flatFieldParameters,
                sections,      // data rows (vertical, tilted, CW - positive) - just dimensions used
                dmargins[2],   // margin for the beginning of each data row
                height,        // image height (radius==1.0 is half height)
                y01234,        // horizontal positions of the center (half height!) of each row
                0.5*width,     // hor. location of the center (for weight calculations) from top left corner of full image (no margins) 
                0.5*height,    // same for vertical.
                tilt);         // rows tilt (hor pixels per pixel in vertical direction, tilted CW - positive, CCW - negative 
   	    int numSect=y01234.length;
    	double [] coeffExp=  {max,max,max,max,max,0.0,2.0,0.5*height,0.5*width};
    	double [] coeffPoly= {max,max,max,max,max,0.0,0.0,0.0,0.0,0.5*height,0.5*width};
    	double [] coefficients= (flatFieldParameters.functionType==1)?coeffExp:coeffPoly;

    	ApproximationState approximationState = new ApproximationState(coefficients,flatFieldParameters.LM_lambdaInitial,bmasks);
    	
    	if (!decideFF1DIteration(flatFieldParameters, approximationState, coefficientNames) ) return null;    	
    	
    	double [][] jacobian=null;
    	double [][] JtJ=null;
    	double [] rightSide=null;
    	double [] diffVector;
    	double [] deltas;
    	double r0=0.5*height; // height, not clearHeight here
//    	double [] modelPixels=null;
//     	double [][] modelSections=null;
    	while (!approximationState.finished) { // iterate
//calculate rms error
    		approximationState.rms=
    			vectorRMS(calcFF1DDiff( flatFieldParameters,sections, dmargins[2],r0, y01234, tilt, approximationState.coeff), errorWeight);
        	if (!decideFF1DIteration(flatFieldParameters,approximationState, coefficientNames) ) break; // end iterations    	        	
        	while (true) {
// calculate deltas for the coefficients        		
        		jacobian= calcFF1DJacobian( flatFieldParameters,clearHeight,dmargins[2],r0,y01234,tilt,approximationState.coeff);
        		JtJ=jacobianByJacobian(jacobian,errorWeight, approximationState.getMask(),approximationState.lambda);
        		diffVector=calcFF1DDiff( flatFieldParameters,sections, dmargins[2],r0, y01234, tilt, approximationState.coeff);
        		approximationState.rms= vectorRMS(diffVector, errorWeight);
        		rightSide=jacobianByVector(jacobian, diffVector, errorWeight, approximationState.getMask());
        		if (DEBUG_LEVEL>1) {
        			System.out.println("Old coefficients: ");
//coefficientNames
        			for (int ii=0;ii<approximationState.coeff.length;ii++)
        				System.out.println((approximationState.getMask()[ii]?"* ":"- ")+"    "+coefficientNames[ii]+"="+approximationState.coeff[ii]);
        			for (int ii=0;ii<rightSide.length;ii++) System.out.println("   rightSide["+ii+"]="+rightSide[ii]);
        			System.out.println("RMS="+approximationState.rms);
        		}

// solve JtJ*delta=      rightSide;
        		Matrix M=new Matrix(JtJ);
        		Matrix mRightSide=new Matrix(rightSide,rightSide.length); // 1 column
        		if (DEBUG_LEVEL>2) {
        			double [][] debug_M=M.getArray();
        			double [][] debug_b=mRightSide.getArray();
        			System.out.println("M["+debug_M.length+"]["+debug_M[0].length+"]");
        			System.out.println("mRightSide["+debug_b.length+"]["+debug_b[0].length+"]");
        			
        		}
        		Matrix mDeltas = M.solve(mRightSide);
        		
        		deltas=mDeltas.getColumnPackedCopy();
        		if (DEBUG_LEVEL>1) {
        			System.out.println("Deltas for lambda="+approximationState.lambda);
        			for (int i=0;i<deltas.length;i++)
        			  System.out.println("  delta["+i+"]="+deltas[i]);
        		}
// save original values
  // approximationState.save() // should be already saved?
// apply deltas 
        		int j=0;
        		for (int i=0;i<approximationState.coeff.length;i++) if (approximationState.getMask()[i]) approximationState.coeff[i]+=deltas[j++];
// calculate new rms;
        		approximationState.rms=
        			vectorRMS(calcFF1DDiff( flatFieldParameters,sections, dmargins[2],r0, y01234, tilt, approximationState.coeff), errorWeight);
        	    
        	    if (!decideFF1DIteration(flatFieldParameters,approximationState, coefficientNames) ) {  	    	
// roll back
        	    	approximationState.restore(); // already done inside?
            		break; // undo last step, continue outer loop
            	}
                if (approximationState.finished) break; // finished iterations
                if ((DEBUG_LEVEL>1) && (!flatFieldParameters.LM_auto)) {
                  SDFA_INSTANCE.showArrays(buildFFModel(flatFieldParameters,numSect,width,height,approximationState.coeff),
                		  width,height, title+"-rms-"+approximationState.rms);
                }
// just continue with improving results
                //buildFFModel(int width,int height,double [] coeff)
        	}
    	}
        if (DEBUG_LEVEL>1) {
        	SDFA_INSTANCE.showArrays(buildFFModel(flatFieldParameters,numSect,width,height,approximationState.coeff),
        			width,height, title+"-final-"+approximationState.coeff[3]+"_"+approximationState.coeff[4]+"_"+approximationState.coeff[5]+"_"+approximationState.coeff[6]);
        }
    	return approximationState.coeff;
    }
   
	private class ApproximationState {
		public double lambda;
		public int    step;
		public int    pass;
		public int [] masks;
		public double rms;
		public double old_rms;
		public int    maskLength;
		public double [] coeff;
		public double old_lambda;
		public double [] old_coeff=null;
		public double lastImprovement=0;
		public boolean finished=false;
		public ApproximationState(
				double [] coeff,
				double lambda,
				boolean [][] bmasks){
			this.coeff=coeff.clone();
			this.lambda=lambda;
			this.masks= new int[bmasks.length+1];
			this.maskLength=bmasks[0].length;
			for (int i=0;i<bmasks.length;i++){
				masks[i]=0;
				for (int j=0;j<bmasks[i].length;j++) if (bmasks[i][j]) masks[i] |= (1<<j);
			}
			masks[masks.length-1]=0;
		}
		public int getIMask(){return this.masks[pass];}
		public boolean[] getMask(){
			boolean [] mask = new boolean [this.maskLength];
			for (int i=0;i<this.maskLength;i++) mask[i]=((masks[pass] & (1<<i))!=0);
			return mask;
		}
//		public void setMask(int imask){this.masks[pass]=imask;}
		public void setMask(boolean[] mask){
			this.masks[pass]=0;
			for (int i=0;i<this.maskLength;i++) if (mask[i]) this.masks[pass] |= (1<<i);
		}
		public void save() {
			this.old_rms=this.rms;
			this.old_lambda=this.lambda;
			this.old_coeff=this.coeff.clone();
		}
		public boolean restore() {
			if (this.old_coeff==null) return false;
			this.lambda=this.old_lambda;
			this.coeff=this.old_coeff;
			this.rms=this.old_rms;
			this.old_coeff=null;
			return true;
		}
	}

	private double []  calcFF1DDiff(
			FlatFieldParameters flatFieldParameters,
			double [][] sections,
			int margin,
			double r0,
			double [] y01234,
			double tilt,
			double [] coeff) {
		return (flatFieldParameters.functionType==1)?
				calcFF1DDiffExp(flatFieldParameters, sections, margin, r0, y01234, tilt,coeff):
					calcFF1DDiffPoly(flatFieldParameters, sections, margin, r0, y01234, tilt,coeff);
	}
    private double []  calcFF1DDiffPoly(
    		FlatFieldParameters flatFieldParameters,
    		double [][] sections,
    		int margin,
    		double r0,
    		double [] y01234,
    		double tilt,
    		double [] coeff) {
    	int numSect=y01234.length;
    	int width=sections[0].length;
    	double [] result = new double [sections.length*width];
    	double [] r2=new double[sections.length];
    	double rr02=1.0/(r0*r0); // width is vertical,
    	double x2;
//    	double [] a0={coeff[0],coeff[1],coeff[2]};
    	double [] a0=new double[numSect];
    	for (int i=0;i<numSect;i++) a0[i]=coeff[i];
    	double a1= coeff[numSect+0];
    	double a2= coeff[numSect+1];
    	double a3= coeff[numSect+2];
    	double a4= coeff[numSect+3];
//    	double a5= coeff[numSect+4];
//    	double x0= coeff[numSect+5];
//    	double y0= coeff[numSect+6];
    	double x0= coeff[numSect+4];
    	double y0= coeff[numSect+5];
    	double y,y2;
    	double [] tilts=new double[numSect];
    	for (int n=0;n<numSect;n++) tilts[n]=tilt;
    	if (flatFieldParameters.noTiltEdges) {
    		tilts[1]=0.0;
    		tilts[2]=0.0;
    	}
    	for (int i=0;i<width;i++) {
    		int ia=i+margin;
    		x2=(ia-x0)*(ia-x0);
    		for (int n=0;n<sections.length;n++) {
//    			y=y012[n]-(i-r0)*tilt-y0;
    			y=y01234[n]-(ia-r0)*tilts[n]-y0;
    			y2=y*y;
        		r2[n]=rr02*(x2+y2);
//        		result[i+n*width]=   sections[n][i]- a0[n]+ r2[n]* (a1 + r2[n]*(a2 + r2[n]*(a3 + r2[n]*(a4 + r2[n]*a5))));
        		result[i+n*width]=   sections[n][i]- a0[n]+ r2[n]* (a1 + r2[n]*(a2 + r2[n]*(a3 + r2[n]*(a4           ))));
    		}
    	}
    	return result;
    }
  
    private double []  calcFF1DDiffExp(
    		FlatFieldParameters flatFieldParameters,
    		double [][] sections,
    		int margin,
    		double r0,
    		double [] y01234,
    		double tilt,
    		double [] coeff) {
    	int numSect=y01234.length;
    	int width=sections[0].length;
    	double [] result = new double [sections.length*width];
    	double [] r2=new double[sections.length];
    	double rr02=1.0/(r0*r0); // width is vertical,
    	double x2;
    	double [] a0=new double[numSect];
    	for (int i=0;i<numSect;i++) a0[i]=coeff[i];
//    	double [] a0={coeff[0],coeff[1],coeff[2]};
    	double a1= coeff[numSect+0];
    	double a2= coeff[numSect+1];
    	double x0= coeff[numSect+2];
    	double y0= coeff[numSect+3];
    	double y,y2;
    	double [] tilts=new double[numSect];
    	for (int n=0;n<numSect;n++) tilts[n]=tilt;
    	if (flatFieldParameters.noTiltEdges) {
    		tilts[1]=0.0;
    		tilts[2]=0.0;
    	}

    	for (int i=0;i<width;i++) {
    		int ia=i+margin;
    		x2=(ia-x0)*(ia-x0);
    		for (int n=0;n<sections.length;n++) {
//    			y=y012[n]-(i-r0)*tilt-y0;
    			y=y01234[n]-(ia-r0)*tilts[n]-y0;
    			y2=y*y;
        		r2[n]=rr02*(x2+y2);
        		result[i+n*width]=   sections[n][i]- a0[n]+ a1*Math.pow(r2[n],a2);
    		}
    	}
    	return result;
    }

    private double [][] calcFF1DJacobian(
    		FlatFieldParameters flatFieldParameters,
    		int width,
    		int margin,
    		double r0,
    		double [] y01234,
    		double tilt,
    		double [] coeff){
    	return calcFF1DJacobian(
    			flatFieldParameters,
        		0, //int img_width,
        		0, //int img_height,
        		width,
        		margin,
        		r0,
        		y01234,
        		tilt,
        		coeff);
    }
    private double [][] calcFF1DJacobian(FlatFieldParameters flatFieldParameters,
    		int img_width,
    		int img_height,
    		int width,
    		int margin,
    		double r0,
    		double [] y01234,
    		double tilt,
    		double [] coeff){
    	
		return (flatFieldParameters.functionType==1)?
       			calcFF1DJacobianExp( flatFieldParameters,img_width,img_height,width,margin,r0,y01234,tilt,coeff):
       			calcFF1DJacobianPoly(flatFieldParameters,img_width,img_height,width,margin,r0,y01234,tilt,coeff);
    	
    }
    
    
    private double [][] calcFF1DJacobianPoly(FlatFieldParameters flatFieldParameters,
    		int img_width,
    		int img_height,
    		int width,
    		int margin,
    		double r0,
    		double [] y01234,
    		double tilt,
    		double [] coeff){
    	int numSect=y01234.length;

//    	f1(x)=(a01- a1* ((x-x0)^2+(y1-y0)^2) - a2* ((x-x0)^2+(y1-y0)^2)^2 -  a3* ((x-x0)^2+(y1-y0)^2)^3) // first band at y1
//    	f2(x)=(a02- a1* ((x-x0)^2+(y2-y0)^2) - a2* ((x-x0)^2+(y2-y0)^2)^2 -  a3* ((x-x0)^2+(y2-y0)^2)^3) // second band at y2
    	double [][] jacobian=new double[numSect*width][coeff.length];
    	double [] r2=new double[numSect];
    	double rr02=1.0/(r0*r0);
    	double x2;
    	double a1= coeff[numSect+0];
    	double a2= coeff[numSect+1];
    	double a3= coeff[numSect+2];
    	double a4= coeff[numSect+3];
//    	double a5= coeff[numSect+4];
//    	double x0= coeff[numSect+5];
//    	double y0= coeff[numSect+6];
    	double x0= coeff[numSect+4];
    	double y0= coeff[numSect+5];

    	double y,y2;
    	double [] tilts=new double[numSect];
    	for (int n=0;n<numSect;n++) tilts[n]=tilt;
    	if (flatFieldParameters.noTiltEdges) {
    		tilts[1]=0.0;
    		tilts[2]=0.0;
    	}
    	double [] ddr2=new double[numSect]; // -df/dr2
    	for (int i=0; i<width;i++)	 {
    		int ia=i+margin;
    		x2=(ia-x0)*(ia-x0);
    		for (int n=0;n<numSect;n++) {
//    			y=y012[n]-(i-r0)*tilt-y0;
    			y=y01234[n]-(ia-r0)*tilts[n]-y0;
    			y2=y*y;
        		r2[n]=rr02*(x2+y2);
        		for (int l=0;l<numSect;l++)
        			jacobian[i+n*width][l]=(l==n)?1.0:0.0;
        		jacobian[i+n*width][numSect+0]=      -r2[n];
        		jacobian[i+n*width][numSect+1]=      -r2[n]*r2[n];
        		jacobian[i+n*width][numSect+2]=      jacobian[i+n*width][numSect+1]*r2[n];
        		jacobian[i+n*width][numSect+3]=      jacobian[i+n*width][numSect+2]*r2[n];
        		jacobian[i+n*width][numSect+4]=      jacobian[i+n*width][numSect+3]*r2[n];
//        		ddr2[n]=2*rr02*(a1 + r2[n]*(2*a2  + r2[n]*3*a3));
//        		ddr2[n]=2*rr02*(a1 + r2[n]*(2*a2  + r2[n]*(3*a3+r2[n]*(4*a4+r2[n]*a5))));
        		ddr2[n]=2*rr02*(a1 + r2[n]*(2*a2  + r2[n]*(3*a3+r2[n]*(4*a4         ))));
//        		jacobian[i+n*width][numSect+5]=ddr2[n] * (ia-x0);
//        		jacobian[i+n*width][numSect+6]=ddr2[n] * y;
        		jacobian[i+n*width][numSect+4]=ddr2[n] * (ia-x0);
        		jacobian[i+n*width][numSect+5]=ddr2[n] * y;
    		}   		
    	}
    	return jacobian;
    }
    
    private double [][] calcFF1DJacobianExp(FlatFieldParameters flatFieldParameters,
    		int img_width,
    		int img_height,
    		int width,
    		int margin,
    		double r0,
    		double [] y01234,
    		double tilt,
    		double [] coeff){
    	int numSect=y01234.length;
    	double [][] jacobian=new double[numSect*width][coeff.length];
    	double [] r2=new double[numSect];
    	double rr02=1.0/(r0*r0);
    	double x2;
//    	double [] a0={coeff[0],coeff[1],coeff[2]};
    	double a1= coeff[numSect+0];
    	double a2= coeff[numSect+1];
    	double x0= coeff[numSect+2];
    	double y0= coeff[numSect+3];
       	double y,y2;
    	double [] tilts=new double[numSect];
    	for (int n=0;n<numSect;n++) tilts[n]=tilt;
    	if (flatFieldParameters.noTiltEdges) {
    		tilts[1]=0.0;
    		tilts[2]=0.0;
    	}
    	double [] ddr2=new double[numSect]; // -df/dr2
    	double [] dpixels =null;
    	if ((img_width>0) && (img_height>0)) {
    		dpixels = new double [img_width*img_height];
    		for (int i=0;i<dpixels.length;i++) dpixels[i]=0.0;
    	}
    	for (int i=0; i<width;i++)	 {
    		int ia=i+margin;
    		x2=(ia-x0)*(ia-x0);
    		for (int n=0;n<numSect;n++) {
//    			y=y012[n]-(i-r0)*tilt-y0;
    			y=y01234[n]-(ia-r0)*tilts[n]-y0;
    			y2=y*y;
        		r2[n]=(x2+y2);
        		if (r2[n]<1.0) r2[n]=1.0;
        		r2[n]*=rr02;
        		for (int l=0;l<numSect;l++)
        			jacobian[i+n*width][l]=(l==n)?1.0:0.0;
        		jacobian[i+n*width][numSect+0]=      -Math.pow(r2[n],a2);
        		jacobian[i+n*width][numSect+1]=      -a1*Math.log(r2[n])*Math.pow(r2[n], a2);
        		ddr2[n]=2*rr02*(a1 * a2 * Math.pow(r2[n], a2-1));
        		jacobian[i+n*width][numSect+2]=ddr2[n] * (ia-x0);
        		jacobian[i+n*width][numSect+3]=ddr2[n] * y;
        		if (dpixels!=null) dpixels[ia*img_width+((int) (y+y0+0.5))]=r2[n];
    		}   		
    	}
    	if (dpixels!=null) SDFA_INSTANCE.showArrays(dpixels,img_width,img_height,"calcFF1DJacobianExp-r2");
		System.out.println("    x0= "+x0);
		System.out.println("    y0= "+y0);
		for (int l=0;l<numSect;l++) {
			System.out.println("    y01234["+l+"]="+y01234[l]);
			
		}

    	return jacobian;
    }
    private double [] buildFFModel(
    		FlatFieldParameters flatFieldParameters,
    		int numSect,
    		int width,
    		int height,
    		double [] coeff){ // width is horizontal, long 
    	if (flatFieldParameters.functionType==1) return buildFFModelExp(numSect,width,height,coeff);    	
    	else                                     return buildFFModelPoly(numSect,width,height,coeff);    	
    }	 
     private double [] buildFFModelPoly(int numSect, int width,int height,double [] coeff){ // width is horizontal, long 
    	double rr02=4.0/(height*height);
    	double [] result = new double [width*height];
    	double x2,r2;
//    	double a0=(coeff[0]+coeff[1]+coeff[2])/3.0;
     	double a0=0.0;
     	for (int i=0;i<numSect;i++) a0+=coeff[i];
     	a0/=numSect;
    	double a1= coeff[numSect+0];
    	double a2= coeff[numSect+1];
    	double a3= coeff[numSect+2];
    	double a4= coeff[numSect+3];
//    	double a5= coeff[numSect+4];
//    	double x0= coeff[numSect+5];
//    	double y0= coeff[numSect+6];
    	double x0= coeff[numSect+4];
    	double y0= coeff[numSect+5];
    	for (int x=0;x<height;x++) {
    		x2=(x-x0)*(x-x0);
        	for (int y=0;y<width;y++) {
        		r2=rr02*(x2+(y-y0)*(y-y0));
//        		result[x*width+y]=a0-r2*(a1+r2*(a2+r2*(a3+r2*(a4+r2*a5))));
        		result[x*width+y]=a0-r2*(a1+r2*(a2+r2*(a3+r2*(a4      ))));
        	}
    	}
    	return result;
    }

     private double [] buildFFModelExp (int numSect, int width,int height,double [] coeff){ // width is horizontal, long 
     	double rr02=4.0/(height*height);
     	double [] result = new double [width*height];
     	double x2,r2;
//     	double a0=(coeff[0]+coeff[1]+coeff[2])/3.0;
     	double a0=0.0;
     	for (int i=0;i<numSect;i++) a0+=coeff[i];
     	a0/=numSect;
     	double a1= coeff[numSect+0];
     	double a2= coeff[numSect+1];
     	double x0= coeff[numSect+2];
     	double y0= coeff[numSect+3];
     	for (int x=0;x<height;x++) {
     		x2=(x-x0)*(x-x0);
         	for (int y=0;y<width;y++) {
         		r2=rr02*(x2+(y-y0)*(y-y0));
         		result[x*width+y]=a0-a1*Math.pow(r2, a2);
         	}
     	}
     	return result;
     }
     
     //jacobian[0].length==mask.length
    //jacobian.length==vector.length
/*     
    private double [] jacobianByVector(double [][] jacobian, double[] vector, boolean [] mask) {
    	return jacobianByVector(jacobian, vector,null,  mask);
    }
*/    
    private double [] jacobianByVector(double [][] jacobian, double[] vector, double [] weights, boolean [] mask) {
    	int size=0;
    	int [] indices=     new int [mask.length];
    	for (int i=0;i<mask.length;i++) if (mask[i]) indices[size++]=i;
    	double [] result = new double [size];
        for (int i=0;i<size;i++){
            result[i]=0.0;
        	if (weights==null)  for (int j=0;j<jacobian.length;j++) result[i]+=jacobian[j][indices[i]]*vector[j];
        	else for (int j=0;j<jacobian.length;j++) result[i]+=jacobian[j][indices[i]]*vector[j]*weights[j];
        }
    	return result;
    }

    //jacobian[0].length==mask.length
   @SuppressWarnings("unused")
private double [][] jacobianByJacobian(double [][] jacobian, boolean [] mask) {
    	return jacobianByJacobian(jacobian, null,mask);
    }
    
    private double [][] jacobianByJacobian(double [][] jacobian, double [] weights, boolean [] mask) {
    	if (jacobian[0].length!=mask.length){
    		System.out.println("Error: jacobian[0].length="+jacobian[0].length+", mask.length="+mask.length);
		    IJ.showMessage("Error","jacobian[0].length="+jacobian[0].length+", mask.length="+mask.length);
    		return null;
    	}
    	int size=0;
    	int [] indices=     new int [mask.length];
    	for (int i=0;i<mask.length;i++) if (mask[i]) indices[size++]=i;
//   System.out.println("size="+size);
//   for (int ii=0;ii<mask.length;ii++) 	System.out.println("mask["+ii+"]="+mask[ii]);
//   for (int ii=0;ii<size;ii++) 	System.out.println("indices["+ii+"]="+indices[ii]);
    	double [][] result = new double [size][size];
    	if (weights==null) { 
    		for (int i=0;i<size;i++) for (int j=0; j<size; j++){
    			result[i][j]=0.0; // 8 - oob
    			for (int k=0;k<jacobian.length;k++) result[i][j]+=jacobian[k][indices[i]]*jacobian[k][indices[j]];
    		}
    	} else {
    		for (int i=0;i<size;i++) for (int j=0; j<size; j++){
    			result[i][j]=0.0; // 8 - oob
    			for (int k=0;k<jacobian.length;k++) result[i][j]+=jacobian[k][indices[i]]*jacobian[k][indices[j]]*weights[k];
    		}
    	}
    	return result;
    }
    
    @SuppressWarnings("unused")
	private double [][] jacobianByJacobian(double [][] jacobian, boolean [] mask, double lambda) {
        return jacobianByJacobian(jacobian, null, mask, lambda);
    }
    
    private double [][] jacobianByJacobian(double [][] jacobian, double [] weights, boolean [] mask, double lambda) {
    	double [][]result=jacobianByJacobian(jacobian, weights, mask);
    	for (int i=0;i<result.length;i++) result[i][i]+=lambda; // this is just Levenberg
    	for (int i=0;i<result.length;i++) result[i][i]*=(1.0+lambda); // this is LMA
    	return result;
    }
/*    
   private double vectorRMS(double [] vector) {
   	return vectorRMS(vector, null);
   }
*/   
    private double vectorRMS(double [] vector, double [] weight) {
    	double l2=0.0;
    	if (weight==null) for (int i=0;i<vector.length;i++) l2+= vector[i]* vector[i];
    	else              for (int i=0;i<vector.length;i++) l2+= vector[i]* vector[i]*weight[i];
    	return Math.sqrt(l2/vector.length);
    }
    /*
    f1(x)=a01- a1* ((x-x0)^2+(y1-y0)^2) - a2* ((x-x0)^2+(y1-y0)^2)^2 -  a3* ((x-x0)^2+(y1-y0)^2)^3
    r2=((x-x0)^2+(y1-y0)^2)
    f1(x)=a01- a1* r2 - a2* r2^2 -  a3* r2^3
    df1/dx0= df1/dr2 * dr2/dx0
    dr2/dx0=-2*(x-x0)
    df1/dr2= -(a1 +2*a2*r2 +3*a3*r2^2)
    df1/dx0= 2*(a1 +2*a2*r2 +3*a3*r2^2)*(x-x0)
     */

	private boolean decideFF1DIteration(FlatFieldParameters flatFieldParameters,
			                            ApproximationState approximationState,
			                            String [] coefficientNames) {
		if (!flatFieldParameters.LM_auto) return calcFF1DIterationDialog(
                approximationState,
				coefficientNames);
//		if (flatFieldParameters.functionType!=1) return false; // not automated (yet?)
		double improvement=-1.0;
		if ((approximationState.old_coeff!=null) && (approximationState.old_rms>0))
			improvement=(approximationState.old_rms-approximationState.rms)/approximationState.old_rms;
		boolean thresholdFired=(approximationState.lastImprovement>0.0) &&
		                       (approximationState.lastImprovement<flatFieldParameters.LM_thresholdFinish) &&
		                        (improvement>-0.001*flatFieldParameters.LM_thresholdFinish) &&
		                        (improvement<flatFieldParameters.LM_thresholdFinish);
		approximationState.lastImprovement=improvement;
		approximationState.step++;
		if (thresholdFired) { // this will happen only after improvement, so we do not need to modify lambda
			if (improvement<0.0) approximationState.restore();
			approximationState.save();
			approximationState.pass++;
			if (approximationState.getIMask()==0){
				if (DEBUG_LEVEL>1) System.out.println("Iteration #"+approximationState.step+" - finished with rms="+approximationState.rms);
				approximationState.finished=true;
				return true;
			} else {
				if (DEBUG_LEVEL>1) System.out.println("Iteration #"+approximationState.step+" - finished "+approximationState.pass+" pass with rms="+approximationState.rms);
				return true;
			}
		}
	    if(approximationState.step > flatFieldParameters.LM_numIterations) {
		  if (DEBUG_LEVEL>1) System.out.println("Iteration #"+approximationState.step+" - ** TOO MANY ITERATIONS **, finished with rms="+approximationState.rms);
		  approximationState.finished=true;
		  return true;
	    }

		
		if (approximationState.old_coeff==null) { // first of first after roll-back call
			if (approximationState.rms>0) approximationState.save();
			return true; // nothing yet to compare 
		}
		if (improvement>0) {
			approximationState.lambda*=flatFieldParameters.LM_lambdaStepDown; //was improvement, reduce lambda
    		if (DEBUG_LEVEL>1) System.out.println("Iteration #"+approximationState.step+", improvement "+(100*improvement)+
    				"%, reducing lambda to "+approximationState.lambda);
			approximationState.save();
			return true;
		} else {
			approximationState.restore();
			approximationState.lambda*=flatFieldParameters.LM_lambdaStepUp; //worsened,increase lambda
			if (approximationState.lambda>1000000.0) approximationState.lambda =1000000.0;
			if (DEBUG_LEVEL>1) System.out.println("Iteration #"+approximationState.step+", worsened by "+(-100*improvement)+
					"%, rolling back, increasing lambda to "+approximationState.lambda);
			return false; //rolled back
		}
	}
	
	private boolean calcFF1DIterationDialog(
            ApproximationState approximationState,
			String [] coefficientNames) {
		    int i;
            boolean [] mask=approximationState.getMask();
		    GenericDialog gd = new GenericDialog("Flat-field estimation iteration step");
		    
	    	gd.addMessage("Average error - "+approximationState.rms+
	    			((approximationState.old_coeff!=null)?( " (was "+approximationState.old_rms+") , delta="+
	    					(100.0*(approximationState.rms-approximationState.old_rms)/approximationState.old_rms)+"%"):""));
    		gd.addMessage("  New coefficients");
    		for (i=0;i<coefficientNames.length;i++) {
    			if (approximationState.old_coeff!=null)
    				gd.addNumericField(i+" ("+coefficientNames[i]+"): ",approximationState.coeff[i], 5,9,
    						(mask[i]?"*":" ")+"("+coefficientNames[i]+"): "+approximationState.old_coeff[i]);
    			else gd.addNumericField(i+" ("+coefficientNames[i]+"): ",     approximationState.coeff[i], 5);
    		}
    		gd.addMessage("  Will be adjusted");
    		for (i=0;i<coefficientNames.length;i++) {
    			gd.addCheckbox(i+" ("+coefficientNames[i]+"): ",     mask[i]);
    		}
			gd.addNumericField("Lambda ",     approximationState.lambda, 5);
			gd.addNumericField("Pass ",     approximationState.pass+1, 0);
			gd.addCheckbox("Exit iterations",     false);
		    gd.showDialog();
		    if (gd.wasCanceled()) {
		    	approximationState.restore();
		    	return false;
		    }
    		for (i=0;i<coefficientNames.length;i++) {
    			approximationState.coeff[i]= gd.getNextNumber();
    		}
    		for (i=0;i<coefficientNames.length;i++) {
    			mask[i]=  gd.getNextBoolean();
    		}
    		approximationState.lambda=         gd.getNextNumber();
    		approximationState.setMask(mask);
    		approximationState.save();
    		approximationState.pass=     (int) gd.getNextNumber()-1;
    		if (gd.getNextBoolean() || (approximationState.lambda<0)) approximationState.finished=true; // exit
		    return true;
	   }
    
    
    // positive - clockwise, negative - CCW
    private double calculateFlatFieldTilt(FlatFieldParameters flatFieldParameters, double [] dpixels, int width, String title) {
    	int height= dpixels.length/width;
    	int [] dmargins=new int [flatFieldParameters.margins.length];
    	for (int i=0;i<dmargins.length;i++) dmargins[i]=flatFieldParameters.margins[i]/flatFieldParameters.decimate;
    	int clearWidth= width- dmargins[0]-dmargins[1];
    	int clearheight=height-dmargins[2]-dmargins[3];
    	double highPassSigma=flatFieldParameters.highPassSigma/flatFieldParameters.decimate;
		double [] highPassPixels= dpixels.clone();
// 		fill right/left margins before using 1-d gaussian blur
    	for (int i=dmargins[2];i<height-dmargins[3];i++) {
    		int indx=i*width+dmargins[0];
    		for (int j=0;j<dmargins[0];j++) highPassPixels[i*width+j]=highPassPixels[indx]; 
    		indx=i*width+width-1-dmargins[1];
    		for (int j=width-1-dmargins[1];j<width;j++) highPassPixels[i*width+j]=highPassPixels[indx]; 
    	}
 //   	SDFA_INSTANCE.showArrays(dpixels,width,height, title+"-downsampled");
		DoubleGaussianBlur gb=new DoubleGaussianBlur();
//		gb.blurDouble(highPassPixels, width, height, highPassSigma, highPassSigma, 0.01);
		gb.blur1Direction(highPassPixels, width, height, highPassSigma, 0.01,true);
		
/*
	    public void blur1Direction(double [] pixels,
	    		                   int        width,
	    		                   int       height,
	    		                   double     sigma,
	    		                   double   accuracy,
	                               boolean xDirection
	                               ) {
		
 */
		
		for (int i=0;i<dpixels.length;i++) highPassPixels[i]=dpixels[i]-highPassPixels[i];
		if (DEBUG_LEVEL>2) SDFA_INSTANCE.showArrays(highPassPixels,width,height, title+"-highpass1d");
// calculate camera average tilt to horizon (in this mode horizon is vertical)
    	int tiltRange= (int) ((clearheight-1)*flatFieldParameters.maxTilt);
    	double [] tiltCorr=new double[2*tiltRange+1];
    	int lastRow=clearheight-1;
    	int lastPixel=tiltCorr.length-1;
    	double dMax=0.0;
    	int iMax=-1;
    	for (int i=0;i<tiltCorr.length;i++) {
    		tiltCorr[i]=0.0;
    		for (int j=0;j<clearWidth-tiltCorr.length;j++) {
        		double lineSum=0.0;
                for (int k=0;k<clearheight;k++)  {
                	int i1=((lastRow-k)*i + k*(lastPixel-i))/lastRow;
                	 lineSum+=highPassPixels[width*(dmargins[2]+k) + (i1+j+dmargins[0])];
                }
    			tiltCorr[i]+=lineSum*lineSum;
    		}
    		if (tiltCorr[i]>dMax) {
    			dMax=tiltCorr[i];
    			iMax=i;
    		}
//    	    if (DEBUG_LEVEL>1) System.out.println("tilt="+(i-tiltRange)+", corr="+tiltCorr[i]);
    	}
    	double tilt=(2.0*(iMax-tiltRange))/lastRow;
	    if (DEBUG_LEVEL>2) System.out.println("Optimal tilt="+(tilt)+" ("+(iMax-tiltRange)+")");
	    if (DEBUG_LEVEL>1) {
	    	int j=(clearWidth-tiltCorr.length)/2;
            for (int k=0;k<clearheight;k++)  {
            	int i1=((lastRow-k)*iMax + k*(lastPixel-iMax))/lastRow;
            	 highPassPixels[width*(dmargins[2]+k) + (i1+j+dmargins[0])]=0.0;
            }
           	SDFA_INSTANCE.showArrays(highPassPixels,width,height, title+"-tilt");
	    }
	    return tilt;
    }
    private void initFlatFieldArrays (int width, double [][] masks, int [] ranges, double [][] weights) {
    	for (int i=0;i<weights.length;i++) {
    		weights[i]=new double[width];
    		for (int j=0;j<width;j++) weights[i][j]=0.0; 
    	}
    	ranges[0]=0;
    	ranges[1]=width/4;
    	ranges[2]=3*width/4;
    	ranges[3]=width;
    	for (int i=0;i<8;i++) masks[i]=new double[width];
    	double k=2*Math.PI/width;
    	for (int j=0;j<width/2;j++) {
    		masks[1][j]=0.5*(1.0+Math.cos(j*k));
    		masks[2][j]=1.0-masks[1][j];
    		masks[4][j]=0.0;
    	}
    	for (int j=width/2;j<width;j++) {
    		masks[1][j]=0.0;
    		masks[4][j]=0.5*(1.0+Math.cos(j*k));
    		masks[2][j]=1.0-masks[4][j];
    	}
    	for (int j=0;j<width;j++) {
    		masks[0][j]=0.0;
    		masks[3][j]=masks[1][j]+masks[2][j];
    		masks[5][j]=masks[1][j]+masks[4][j];
    		masks[6][j]=masks[2][j]+masks[4][j];
    		masks[7][j]=1.0;
    	}
    }
// images cap[ture in portrait mode, opened in landscape mode (rotated CCW90 from original)
// returns valid images for sub-bands: top, middle and bottom:   mask bits 1- top (left) valid, 2 - middle valid, 4 - bottom (right) vlaid     
    private int calcValidFlatFieldMask(FlatFieldParameters flatFieldParameters, int [] ranges, ImagePlus imp){
    	int mask=0, maskO=0, maskU=0;
    	int [] mranges=ranges.clone();
    	mranges[0]=flatFieldParameters.margins[0];
    	mranges[3]-=flatFieldParameters.margins[1];
    	double [] oMap=JP4_INSTANCE.overexposedMap (imp, flatFieldParameters.overExpValue);
    	double [] uMap=JP4_INSTANCE.overexposedMap (imp, flatFieldParameters.underExpValue);
    	int width=imp.getWidth();
    	int height=imp.getHeight()-flatFieldParameters.margins[2]-flatFieldParameters.margins[3];
    	for (int i=0;i<3;i++) {
           if  (JP4_INSTANCE.fracOverExposed(oMap,   // map of overexposed pixels 0.0 - 0K, >0 (==1.0) - overexposed
        		   width,    // width of the map
        		   ranges[i],          // X of the top left corner of the selection
        		   flatFieldParameters.margins[2],          // Y of the top left corner of the selection
        		   ranges[i+1]-ranges[i],  // selection width
        		   height)<flatFieldParameters.overExpFrac) // selection height
        	   maskO |= (1<<i);    		
        	   
			if (JP4_INSTANCE.fracOverExposed(uMap,   // map of overexposed pixels 0.0 - 0K, >0 (==1.0) - overexposed
						width,    // width of the map
						ranges[i],          // X of the top left corner of the selection
						flatFieldParameters.margins[2],          // Y of the top left corner of the selection
						ranges[i+1]-ranges[i],  // selection width
						height) >(1.0-flatFieldParameters.underExpFrac)) // selection height
        	   maskU |= (1<<i);    		
    	}
    	mask = (maskU & maskO) | (maskO<<4) | (maskU<<8) ; // only 3 LSB will be finally used, other bits - ju
    	return mask;
    }
    private boolean flatFieldAccumulate(FlatFieldParameters flatFieldParameters,
			ProcessCalibrationFilesParameters processCalibrationFilesParameters){
		if ((flatFieldParameters.flatFieldDirectory==null) || (flatFieldParameters.flatFieldDirectory.length()==0))
			flatFieldParameters.flatFieldDirectory=selectFlatFieldDirectory(flatFieldParameters.flatFieldDirectory);
		if ((flatFieldParameters.flatFieldDirectory==null) || (flatFieldParameters.flatFieldDirectory.length()==0)) return false; // nothing to do
	    if (DEBUG_LEVEL>1) System.out.println("flatFieldAccumulate(): Do flat-field accumulation");
// initialize result arrays
	    int numSensors=flatFieldParameters.eyesisMode?(flatFieldParameters.numEyesisChannels*flatFieldParameters.numEyesisSubChannels):1;
	    float [][] flatFieldPixels= new float[numSensors][];
	    int   [][] numAveraged=       new int[numSensors][3];
	    int   [][] imageDimensions= new int[numSensors][2];
	    for (int i=0;i<numSensors;i++){
	    	flatFieldPixels[i]=null;
	    	numAveraged[i][0]=0;  // top (left)
	    	numAveraged[i][1]=0;  // middle
	    	numAveraged[i][2]=0;  // bottom (right)
	    	imageDimensions[i][0]=0;
	    	imageDimensions[i][1]=0;
	    }
	    double [][] weightMasks = new double [8][];
	    weightMasks[0]=null; 
	    double [][] weights = new double [numSensors][];
	    int [] ranges=new int [4];
	    File srcDir;
		File [] fileList;
		boolean [] usedChannels=null;
//		String [] paths=null;
		String [] names=null;
		int numSubChannels=1;
		boolean [][] subChannels=null;
        int [] fileChannels=null; // channel number of each file, <0 - not used
  	    ImagePlus imp_composite=null;
  	    ImagePlus imp_single=null;
  	    ImageProcessor ip;
  	    float [] pixels;
  	    Runtime runtime = Runtime.getRuntime();
		if (flatFieldParameters.eyesisMode){
			numSubChannels=flatFieldParameters.numEyesisSubChannels;
			usedChannels=new boolean[flatFieldParameters.numEyesisChannels];
			subChannels=new boolean[flatFieldParameters.numEyesisChannels][flatFieldParameters.numEyesisSubChannels];
			for (int i=0;i<flatFieldParameters.numEyesisChannels; i++) for (int j=0;j<flatFieldParameters.numEyesisSubChannels; j++) subChannels[i][j]=false; 
			for (int i=0;i<usedChannels.length; i++) usedChannels[i]=false; 
			for (int i=0;i<flatFieldParameters.processChannels.length; i++) if (flatFieldParameters.processChannels[i] || flatFieldParameters.processAllChannels){
				usedChannels[i/flatFieldParameters.numEyesisSubChannels] = true; 
				subChannels[i/flatFieldParameters.numEyesisSubChannels][i%flatFieldParameters.numEyesisSubChannels]=true;
			}
		} else {
			subChannels=new boolean[1][1];
			subChannels[0][0]=true;
		}
	    for (int numDir=0;numDir<flatFieldParameters.sourceDirPaths.length;numDir++) if (flatFieldParameters.sourceDirPathsEn[numDir]){
//            if (DEBUG_LEVEL>1) System.out.println("??Processing directory "+flatFieldParameters.sourceDirPaths[numDir]);
	    	srcDir= new File (flatFieldParameters.sourceDirPaths[numDir]);
			if (!srcDir.exists()) continue;
//            if (DEBUG_LEVEL>1) System.out.println("++Processing directory "+flatFieldParameters.sourceDirPaths[numDir]+" filter="+("."+flatFieldParameters.sourceFileExtension));
			fileList=srcDir.listFiles(new Filter("."+flatFieldParameters.sourceFileExtension));
//            if (DEBUG_LEVEL>1) System.out.println("--Processing directory "+flatFieldParameters.sourceDirPaths[numDir]+" ("+fileList.length+" files)");
			if ((fileList==null) || (fileList.length==0)) continue; // null pointer
            if (DEBUG_LEVEL>1) System.out.println("Processing directory "+flatFieldParameters.sourceDirPaths[numDir]+" ("+fileList.length+" files)");
//            paths=new String[fileList.length];
            names=new String[fileList.length];
            fileChannels= new int[fileList.length];
            for (int i=0;i<fileList.length;i++) {
//            	paths[i]=fileList[i].getPath();
            	names[i]=fileList[i].getName();
            	if (flatFieldParameters.eyesisMode) {
            	  fileChannels[i]=-1;
            	  for (int j=0;j<flatFieldParameters.numEyesisChannels;j++) if  (usedChannels[j] &&(names[i].indexOf("_"+(j+1)+"-")>=0)) {
            		  fileChannels[i]=j;
            		  break;
            	  }
            	} else fileChannels[i]=0; 
            }
            int numFile=0;
            for (int i=0;i<names.length;i++) if (fileChannels[i]>=0) { 
   			  int channel=fileChannels[i];
                if (DEBUG_LEVEL>1) System.out.println("  Processing file "+(numFile+1)+": "+names[i]);
  			  imp_composite=JP4_INSTANCE.open(
  					flatFieldParameters.sourceDirPaths[numDir]+Prefs.getFileSeparator(), // directory
					  names[i],
					  "",  //arg - not used in JP46 reader
					  true, // un-apply camera color gains
					  null, // new window
					  false); // do not show
			  if (imp_composite==null) continue;
              for (int subChannel=0;subChannel<numSubChannels;subChannel++) if (subChannels[channel][subChannel]){
                if (DEBUG_LEVEL>1) System.out.println("    Processing subchannel "+(subChannel+1));
            	int sensorNumber=channel*numSubChannels+subChannel;
          		if (flatFieldParameters.eyesisMode){
          			imp_single=JP4_INSTANCE.demuxImage(imp_composite, subChannel);
          		} else {
          			imp_single=imp_composite;
          		}
          		ip=imp_single.getChannelProcessor();
// currently all images in all channels should be the same width          		
        	    if (weightMasks[0]==null)initFlatFieldArrays (imp_single.getWidth(), weightMasks, ranges, weights);
        	    int diagMask= calcValidFlatFieldMask(flatFieldParameters, ranges,  imp_single);
        	    int validBandsMask = diagMask & 7;
        	    String imageExpStatus="";
        	    for (int n=0;n<3;n++) {
        	    	if ((diagMask &      (1     << n))!=0) imageExpStatus+=" *";
        	    	else if ((diagMask & (0x10  << n))!=0) imageExpStatus+=" -";
        	    	else if ((diagMask & (0x100 << n))!=0) imageExpStatus+=" +";
        	    	else                                   imageExpStatus+=" X";
        	    }
                if (DEBUG_LEVEL>1) System.out.println("    Sensor # "+sensorNumber+" - image status: "+imageExpStatus);
                if (validBandsMask==0) {
//                  if (DEBUG_LEVEL>1) System.out.println("    no usable areas in the image, skipping");
                	continue;
                }
          		pixels=(float[]) ip.getPixels();
          		if (flatFieldPixels[sensorNumber]==null) {
          			flatFieldPixels[sensorNumber]= new float [pixels.length];
              		for (int n=0;n<pixels.length;n++) flatFieldPixels[sensorNumber][n]=0.0f;
          			imageDimensions[sensorNumber][0]=ip.getWidth();
          			imageDimensions[sensorNumber][1]=ip.getHeight();
           if (DEBUG_LEVEL>1) System.out.println("    New image "+imageDimensions[sensorNumber][0]+"*"+imageDimensions[sensorNumber][1]);
          		}
      			if (flatFieldPixels[sensorNumber].length!=pixels.length) { // TODO: change to width,height?
      				System.out.println("Image size ("+pixels.length+") does not match that of the previous ("+flatFieldPixels[sensorNumber].length+") - skipping");
      				continue;
      			}
          		if (flatFieldParameters.normalize) { // actually - logarithm
//              		for (int n=0;n<pixels.length;n++) flatFieldPixels[sensorNumber][n]+=Math.log(pixels[n]+flatFieldParameters.fatZero);
              		for (int n=0;n<pixels.length;n++) flatFieldPixels[sensorNumber][n]+=Math.log(pixels[n]+flatFieldParameters.fatZero)*weightMasks[validBandsMask][n%imageDimensions[sensorNumber][0]];
          		} else {
//              		for (int n=0;n<pixels.length;n++) flatFieldPixels[sensorNumber][n]+=pixels[n];
          			for (int n=0;n<pixels.length;n++) flatFieldPixels[sensorNumber][n]+=pixels[n]*weightMasks[validBandsMask][n%imageDimensions[sensorNumber][0]];
          		}
          		for (int n=0;n<imageDimensions[sensorNumber][0];n++) weights[sensorNumber][n]+=weightMasks[validBandsMask][n];
// TODO: add normalization per-color? here
//      		numAveraged[sensorNumber]++;
          		for (int n=0;n<3;n++) if ((validBandsMask & (1<<n))!=0)  numAveraged[sensorNumber][n]++;
                if (DEBUG_LEVEL>1) System.out.println("    Number averaged for the sensor "+sensorNumber+": "+
                		numAveraged[sensorNumber][0]+" : "+
                		numAveraged[sensorNumber][1]+" : "+
                		numAveraged[sensorNumber][2]);
//	  ImagePlus [][][] result=null;

//subChannels[fileChannels[i]] - which subchannels to process;
//TODO: do actual processing
//		ask for FlatFieldDirectory if it is not defined
//		  private String selectFlatFieldDirectory(String defaultPath);
              }
              numFile++;
            }
		  runtime.gc();
		  if (DEBUG_LEVEL>1) System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");
            
	    }
	    for (int sensorNumber=0; sensorNumber<flatFieldPixels.length;sensorNumber++) if ((numAveraged[sensorNumber][0]>0) &&
	    		(numAveraged[sensorNumber][1]>0) &&	(numAveraged[sensorNumber][2]>0)) {
//	    	for (int n=0;n<flatFieldPixels[sensorNumber].length;n++) flatFieldPixels[sensorNumber][n]/=numAveraged[sensorNumber];
	    	for (int n=0;n<flatFieldPixels[sensorNumber].length;n++) flatFieldPixels[sensorNumber][n]/=weights[sensorNumber][n%imageDimensions[sensorNumber][0]];

      		if (flatFieldParameters.normalize) { // actually - logarithm
      			for (int n=0;n<flatFieldPixels[sensorNumber].length;n++) flatFieldPixels[sensorNumber][n]=(float) (Math.exp(flatFieldPixels[sensorNumber][n])-flatFieldParameters.fatZero);
      		}
	    }
// just display for now	    
	    for (int sensorNumber=0; sensorNumber<flatFieldPixels.length;sensorNumber++) {
	    	SDFA_INSTANCE.showArrays(flatFieldPixels[sensorNumber],
	    			                 imageDimensions[sensorNumber][0],
	    			                 imageDimensions[sensorNumber][1],
	    			                 "flat-field-"+sensorNumber+"-"+numAveraged[sensorNumber][0]+"-"+
	    		                		numAveraged[sensorNumber][1]+"-"+
	    		                		numAveraged[sensorNumber][2]); 
            if (DEBUG_LEVEL>1) System.out.println("    Number averaged for sensor "+sensorNumber+"(top,middle,bottom): "+
            		numAveraged[sensorNumber][0]+" : "+
            		numAveraged[sensorNumber][1]+" : "+
            		numAveraged[sensorNumber][2]);

	    
	    }
		return true;
	}
	
	private void addFlatFieldSources(FlatFieldParameters flatFieldParameters){
//	  String defaultDirectory=	((flatFieldParameters.sourceDirPaths==null) || (flatFieldParameters.sourceDirPaths.length==0))?"":flatFieldParameters.sourceDirPaths[flatFieldParameters.sourceDirPaths.length-1];
	// will start with already selected ones - is it OK?	
      String [] newDirectories=selectSourceDirectories(flatFieldParameters.sourceDirPaths);
      if (newDirectories==null) return;
      int numNew=0;
      int numOld=(flatFieldParameters.sourceDirPaths==null)?0:flatFieldParameters.sourceDirPaths.length;
      boolean [] newDirs=new boolean[newDirectories.length];
      for (int i=0;i<newDirectories.length;i++) {
    	  newDirs[i]=true;
    	  for (int j=0;j<numOld;j++) if (newDirectories[i].equals(flatFieldParameters.sourceDirPaths[j])) {
    		  newDirs[i]=false;
    		  break;
    	  }
    	  // duplicates in the new? probably impossible, just in case
    	  if (newDirs[i]) for (int j=0;j<i;j++) if (newDirs[j] && (newDirectories[i].equals(newDirectories[j]))) {
    		  newDirs[i]=false;
    		  break;
    	  }
    	  if (newDirs[i]) numNew++;
      }
      if (numNew==0) return; // nothing new
      String [] oldPaths=null;
      boolean [] oldPathsEn=null;
      if (numOld>0) {
        oldPaths=flatFieldParameters.sourceDirPaths.clone();
        oldPathsEn=flatFieldParameters.sourceDirPathsEn.clone();
      }  
      flatFieldParameters.sourceDirPaths=  new String [numOld+numNew];
      flatFieldParameters.sourceDirPathsEn=new boolean[numOld+numNew];
      for (int i=0;i<numOld;i++){
    	  flatFieldParameters.sourceDirPaths[i]=  oldPaths[i];
    	  flatFieldParameters.sourceDirPathsEn[i]=oldPathsEn[i];
      }
      int j=numOld;
      for (int i=0;i<newDirectories.length;i++) if(newDirs[i]){
    	  flatFieldParameters.sourceDirPaths[j]=  newDirectories[i];
    	  flatFieldParameters.sourceDirPathsEn[j++]=true;
      }
	}
	private boolean selectDirectoriesToRemove(FlatFieldParameters flatFieldParameters) {
		GenericDialog gd = new GenericDialog("Remove source directories for flat-field calcualtion");
		gd.addMessage("Select source directories to remove");
		int i;
		int numDirs=(flatFieldParameters.sourceDirPaths==null)?0:flatFieldParameters.sourceDirPaths.length;
		for (i=0;i<numDirs;i++) gd.addCheckbox((flatFieldParameters.sourceDirPathsEn[i]?"+":"-")+" "+flatFieldParameters.sourceDirPaths[i], false); 
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		int numDelete=0;
		boolean [] selToDel=new boolean[numDirs];
		for (i=0;i<numDirs;i++) {
			selToDel[i]=gd.getNextBoolean();
			if (selToDel[i]) numDelete++;
		}
		if (numDelete==0) return false; // nothing to delete, same as cancel
		if (numDelete==numDirs) {
			flatFieldParameters.sourceDirPaths=  null;
			flatFieldParameters.sourceDirPathsEn=null;
			return true; // nothing left
		}
		String [] oldPaths=flatFieldParameters.sourceDirPaths.clone();
		boolean [] oldPathsEn=flatFieldParameters.sourceDirPathsEn.clone();
		flatFieldParameters.sourceDirPaths=  new String [numDirs-numDelete];
		flatFieldParameters.sourceDirPathsEn=new boolean[numDirs-numDelete];
		int j=0;
		for (i=0;i<numDirs;i++) if (!selToDel[i]) {
			flatFieldParameters.sourceDirPaths[j]=oldPaths[i];
			flatFieldParameters.sourceDirPathsEn[j]=oldPathsEn[i];
			j++;
		}
		return true;
}
	private boolean showFlatFieldDialog(FlatFieldParameters flatFieldParameters, boolean dirSelect) {
		return 	showFlatFieldDialog(flatFieldParameters, dirSelect, true);

	}
	private boolean showFlatFieldDialog(FlatFieldParameters flatFieldParameters, boolean dirSelect, boolean showApproxPars) {
		   int i;
		    GenericDialog gd = new GenericDialog("Flat-field parameters");
		    gd.addCheckbox    ("Normalize before averaging",                                                flatFieldParameters.normalize);
		    
			gd.addNumericField("Over-exposure relative value (0..1)",                                       flatFieldParameters.overExpValue, 5);
			gd.addNumericField("Over-exposure allowed fraction (0..1)",                                     flatFieldParameters.overExpFrac, 5);
			gd.addNumericField("Under-exposure relative value (0..1)",                                      flatFieldParameters.underExpValue, 5);
			gd.addNumericField("Under-exposure allowed fraction (0..1)",                                    flatFieldParameters.underExpFrac, 5);
			gd.addNumericField("Fat zero (before calculating ln):",                                         flatFieldParameters.fatZero, 2);
		    
		    if (showApproxPars) {
		      gd.addCheckbox    ("Do not use tilt in edge sections (to get closer to the corners)",           flatFieldParameters.noTiltEdges);
			  gd.addNumericField("Approximation function type (0 - polynomial, 1 - power)",                   flatFieldParameters.functionType, 0);
			  gd.addNumericField("Approximation function modifier",                                           flatFieldParameters.functionModifier, 0);
			  
			  gd.addNumericField("relative distance from the center to the  middle sections 3,4 (1-3-0-4-2)", flatFieldParameters.section34, 3);
			  
			  gd.addNumericField("Weight for the error function in the center (to be added to r^2",           flatFieldParameters.centerWeight, 5);

		      gd.addCheckbox    ("Iterate automatically",                                                     flatFieldParameters.LM_auto);
			  gd.addNumericField("Levenberg–Marquardt initial lambda",                                        flatFieldParameters.LM_lambdaInitial, 5);
			  gd.addNumericField("Levenberg–Marquardt lambda step up (if worsened)",                          flatFieldParameters.LM_lambdaStepUp,  1);
			  gd.addNumericField("Levenberg–Marquardt lambda step down (if improved)",                        flatFieldParameters.LM_lambdaStepDown,2);
			  gd.addNumericField("Levenberg–Marquardt final relative improvement",                            flatFieldParameters.LM_thresholdFinish,6);
			  gd.addNumericField("Levenberg–Marquardt maximal number of iterations",                          flatFieldParameters.LM_numIterations,0);
			

			  gd.addNumericField("Margin left, pixels:",                                                      flatFieldParameters.margins[0], 0);
			  gd.addNumericField("Margin right, pixels:",                                                     flatFieldParameters.margins[1], 0);
			  gd.addNumericField("Margin top, pixels:",                                                       flatFieldParameters.margins[2], 0);
			  gd.addNumericField("Margin bottom, pixels:",                                                    flatFieldParameters.margins[3], 0);
			  gd.addNumericField("Downsacale flat-field image:",                                              flatFieldParameters.decimate, 0);
			  gd.addNumericField("Sample width for initial estimation:",                                      flatFieldParameters.sampleWidth, 0);
			  gd.addNumericField("High-pass filtering for tilt calculation",                                  flatFieldParameters.highPassSigma, 1);
			  gd.addNumericField("MaximalTilt",                                                               flatFieldParameters.maxTilt, 3);
		    }			
			
//			true,   // eyesisMode
		    gd.addCheckbox    ("Eyesis mode (multiple channels)",                                           flatFieldParameters.eyesisMode);
		    int numDirs=0;
		    if (dirSelect) {
		    	gd.addStringField ("Source files extension (no dot):    ",                                      flatFieldParameters.sourceFileExtension, 5);
		    	gd.addMessage("Select channels to process (applicable in eyesis mode only)");
		    	gd.addCheckbox    ("Process all channels (if false will use only individually enabled, below)", flatFieldParameters.processAllChannels);
		    	gd.addCheckbox    ("Process channel 1-1",                                                       flatFieldParameters.processChannels[0]);
		    	gd.addCheckbox    ("Process channel 1-2",                                                       flatFieldParameters.processChannels[1]);
		    	gd.addCheckbox    ("Process channel 1-3",                                                       flatFieldParameters.processChannels[2]);
		    	gd.addCheckbox    ("Process channel 2-1",                                                       flatFieldParameters.processChannels[3]);
		    	gd.addCheckbox    ("Process channel 2-2",                                                       flatFieldParameters.processChannels[4]);
		    	gd.addCheckbox    ("Process channel 2-3",                                                       flatFieldParameters.processChannels[5]);
		    	gd.addCheckbox    ("Process channel 3-1",                                                       flatFieldParameters.processChannels[6]);
		    	gd.addCheckbox    ("Process channel 3-2",                                                       flatFieldParameters.processChannels[7]);
		    	gd.addCheckbox    ("Process channel 3-3",                                                       flatFieldParameters.processChannels[8]);
		    	gd.addMessage("Select directories to process");
		    	numDirs=(flatFieldParameters.sourceDirPaths==null)?0:flatFieldParameters.sourceDirPaths.length;
		    	for (i=0;i<numDirs;i++) gd.addCheckbox(flatFieldParameters.sourceDirPaths[i], flatFieldParameters.sourceDirPathsEn[i]);
		    	gd.addStringField ("Result directory for flat-field calibration files:", flatFieldParameters.flatFieldDirectory, 50);
		    }
		    gd.addCheckbox    ("Save current settings with results",                                        flatFieldParameters.saveSettings);
	        gd.addCheckbox    ("Use XML format to save/restore settings",                                   flatFieldParameters.useXML);
			gd.addNumericField("Debug Level:",                                                              MASTER_DEBUG_LEVEL, 0);
		    gd.showDialog();
		    if (gd.wasCanceled()) return false;
		    flatFieldParameters.normalize=                                   gd.getNextBoolean();
			flatFieldParameters.overExpValue=                                gd.getNextNumber();
			flatFieldParameters.overExpFrac=                                 gd.getNextNumber();
			flatFieldParameters.underExpValue=                               gd.getNextNumber();
			flatFieldParameters.underExpFrac=                                gd.getNextNumber();
		    flatFieldParameters.fatZero=                                     gd.getNextNumber();
		    if (showApproxPars) {
		      flatFieldParameters.noTiltEdges=                                 gd.getNextBoolean();
		      flatFieldParameters.functionType=                          (int) gd.getNextNumber();
		      flatFieldParameters.functionModifier=                      (int) gd.getNextNumber();
		      
		      flatFieldParameters.section34=                                   gd.getNextNumber();
		      flatFieldParameters.centerWeight=                                gd.getNextNumber();
		      flatFieldParameters.LM_auto=                                     gd.getNextBoolean();
		      flatFieldParameters.LM_lambdaInitial=                            gd.getNextNumber();
		      flatFieldParameters.LM_lambdaStepUp=                             gd.getNextNumber();
		      flatFieldParameters.LM_lambdaStepDown=                           gd.getNextNumber();
		      flatFieldParameters.LM_thresholdFinish=                          gd.getNextNumber();
		      flatFieldParameters.LM_numIterations=                      (int) gd.getNextNumber();

		    
			  flatFieldParameters.margins[0]=                            (int) gd.getNextNumber();
			  flatFieldParameters.margins[1]=                            (int) gd.getNextNumber();
			  flatFieldParameters.margins[2]=                            (int) gd.getNextNumber();
			  flatFieldParameters.margins[3]=                            (int) gd.getNextNumber();
			  flatFieldParameters.decimate=                              (int) gd.getNextNumber();
			  flatFieldParameters.sampleWidth=                           (int) gd.getNextNumber();
			  flatFieldParameters.highPassSigma=                               gd.getNextNumber();
			  flatFieldParameters.maxTilt=                                     gd.getNextNumber();
		    }
		    flatFieldParameters.eyesisMode=                                  gd.getNextBoolean();
		    if (dirSelect) {
		    	flatFieldParameters.sourceFileExtension=                         gd.getNextString();
		    	flatFieldParameters.processAllChannels=                          gd.getNextBoolean();
		    	for (i=0; i<flatFieldParameters.processChannels.length;i++)
		    		flatFieldParameters.processChannels[i]=                      gd.getNextBoolean();
		    	for (i=0;i<numDirs;i++) flatFieldParameters.sourceDirPathsEn[i] =gd.getNextBoolean();
		    	flatFieldParameters.flatFieldDirectory=                          gd.getNextString();
		    }
		    flatFieldParameters.saveSettings=                                gd.getNextBoolean();
		    flatFieldParameters.useXML=		                                 gd.getNextBoolean();
			MASTER_DEBUG_LEVEL=                                        (int) gd.getNextNumber();
		    return true;
	   }
	   
	
// end of Flat-field related
	/* ======================================================================== */
/*
	public float[][] simulateGridAll (
			int width, // extend to full image, width, height - original (not scaled) image size
			int height, 
			double [][][][] patternGrid, // should be aligned to gridFrac 
			int gridFrac, // number of grid steps per pattern full period
			SimulationPattern.SimulParameters  simulParameters,
			int       threadsMax,
			boolean   updateStatus,
			int debug_level){// debug level used inside loops
		SimulationPattern simulationPattern=new SimulationPattern(simulParameters);
		float [][] simArray0=simulateGridAll (
				patternGrid, // should be aligned to gridFrac 
				gridFrac, // number of grid steps per pattern full period
				simulParameters,
				simulationPattern,
				threadsMax,
				updateStatus,
				debug_level);
		Rectangle woi=matchSimulatedPattern.getWOI();
		if ((woi.x==0) && (woi.y==0) && (woi.width==width) && (woi.height==height)) return simArray0;
		int k=simulParameters.subdiv/2;
		Rectangle scaledWoi=new Rectangle(k*woi.x, k*woi.y, k*woi.width, k*woi.height);
		float [][] simArray=new float [2][];
		simArray[0]=(new SimulationPattern(simulParameters)).combineWithCanvas(0.0,  k*width, k*height, scaledWoi,simArray0[0]);
		simArray[1]=(new SimulationPattern(simulParameters)).combineWithCanvas(0.0,  k*width, k*height, scaledWoi,simArray0[1]);
		if (DEBUG_LEVEL>1) SDFA_INSTANCE.showArrays(simArray,width*k,height*k,true, "full-simulation");
		return simArray;
	}
	
	public float[][] simulateGridAll (
			double [][][][] patternGrid, // should be aligned to gridFrac 
			int gridFrac, // number of grid steps per pattern full period
			SimulationPattern.SimulParameters  simulParameters,
			SimulationPattern simulationPattern, // or null
			int       threadsMax,
			boolean   updateStatus,
			int debug_level){// debug level used inside loops
		long 	  startTime=System.nanoTime();
		double [][] xy0={{simulParameters.offsetX,simulParameters.offsetY},{simulParameters.offsetX-0.5,simulParameters.offsetY-0.5}} ;
		if (simulationPattern==null) simulationPattern=new SimulationPattern(simulParameters);
		float[][] simArray=new float[2][];
		simArray[0]=  simulationPattern.simulateGrid (
				matchSimulatedPattern.getDArray(),
				2, // gridFrac, // number of grid steps per pattern full period
				simulParameters,
				matchSimulatedPattern.getWOI(),
				simulParameters.subdiv/2,
				xy0[0],    // add to patterGrid xy
				threadsMax,
				updateStatus,
				debug_level); // debug level
		simArray[1]=  simulationPattern.simulateGrid (
				matchSimulatedPattern.getDArray(),
				2, // gridFrac, // number of grid steps per pattern full period
				simulParameters,
				matchSimulatedPattern.getWOI(),
				simulParameters.subdiv/2,
				xy0[1],    // add to patterGrid xy
				threadsMax,
				updateStatus,
				debug_level); // debug level
		if (DEBUG_LEVEL>1) SDFA_INSTANCE.showArrays(simArray,matchSimulatedPattern.getWOI().width*simulParameters.subdiv/2,matchSimulatedPattern.getWOI().height*simulParameters.subdiv/2,true, "a-simulation");
		if (DEBUG_LEVEL>0) System.out.println("Finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		return simArray;
	}	

*/	
	/* ======================================================================== */
	public void saveTimestampedProperties(
			String path,      // full path or null
			String directory, // use as default directory if path==null 
			boolean useXML,
			Properties properties ){
		if (path==null){
   			String msg="path should not be null";
   			IJ.showMessage("Error",msg); 
			throw new IllegalArgumentException (msg);
		}
		saveProperties(
				path+"_"+IJ.d2s(0.000001*(System.nanoTime()/1000),6).replace('.', '_'),      // full path or null
				directory, // use as default directory if path==null 
				useXML,
				properties );
		
	}			

	public void saveProperties(
			String path,      // full path or null
			String directory, // use as default directory if path==null 
			boolean useXML,
			Properties properties ){
   	    String [] XMLPatterns= {".conf-xml",".xml"};
   	    String [] confPatterns={".conf"};
   	    String [] patterns=useXML?XMLPatterns:confPatterns;
     if (path==null) {
	    path= CalibrationFileManagement.selectFile(true, // save  
			  "Save configuration selection", // title
			  "Select configuration file", // button
			  new CalibrationFileManagement.MultipleExtensionsFileFilter(patterns, (useXML?"XML ":"")+"Configuration files ("+(useXML?"*.conf-xml":"*.conf")+")"), // filter
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
		setAllProperties(properties); // NEW. Setting properties from current parameters so missing (in the file) parameters will not cause an error
   	    String [] XMLPatterns= {".conf-xml",".xml"};
	    String [] confPatterns={".conf"};
	    String [] patterns=useXML?XMLPatterns:confPatterns;
	     if (path==null) {
	 	    path= CalibrationFileManagement.selectFile(false, // save  
	 			  "Configuration file selection", // title
	 			  "Read configuration file", // button
	 			  new CalibrationFileManagement.MultipleExtensionsFileFilter(patterns,(useXML?"XML ":"")+"Configuration files ("+(useXML?"*.conf-xml":"*.conf")+")"), // filter
	 			  directory); // may be ""
	      }  else {
	    	  // do not add extension if it already exists
	    	  if ((path.length()<patterns[0].length()) || (!path.substring(path.length()-patterns[0].length()).equals(patterns[0]))){
	    		  path+=patterns[0];
	    	  }
	      }
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
    	boolean select= (properties.getProperty("selected")!=null);
    	boolean select_MASTER_DEBUG_LEVEL=!select;
    	boolean select_SHOW_AS_STACKS=!select;
    	boolean select_FFT_SIZE=!select;
    	boolean select_MAP_FFT_SIZE=!select;
    	boolean select_GAUSS_WIDTH=!select;
    	boolean select_FFT_OVERLAP=!select;
    	boolean select_PSF_SUBPIXEL=!select;
    	boolean select_PSF_SAVE_FILE=!select;
    	boolean select_THREADS_MAX=!select;
    	boolean select_UPDATE_STATUS=!select;
    	boolean select_SIMUL=!select;
    	boolean select_INTERPOLATE=!select;
    	boolean select_INVERSE=!select;
    	boolean select_PSF_PARS=!select;
    	boolean select_OTF_FILTER=!select;
    	boolean select_PATTERN_DETECT=!select;
    	boolean select_COMPONENTS=!select;
    	boolean select_SHOW_RESULTS=!select;
    	boolean select_MULTIFILE_PSF=!select;
    	boolean select_PROCESS_PARAMETERS=!select;
    	boolean select_DISTORTION=!select;
    	boolean select_LASER_POINTERS=!select;
    	boolean select_FLATFIELD_PARAMETERS=!select;
    	boolean select_PATTERN_PARAMETERS=!select;
    	boolean select_LENS_DISTORTION_PARAMETERS=!select;
    	boolean select_EYESIS_CAMERA_PARAMETERS=!select;
    	boolean select_LASERS=!select;
    	boolean select_CAMERAS=!select;
    	boolean select_DISTORTION_PROCESS_CONFIGURATION=!select;
    	boolean select_REFINE_PARAMETERS=!select;
    	boolean select_FOCUS_MEASUREMENT_PARAMETERS=!select;
    	boolean select_MOTORS_focusingHistory=!select;
    	boolean select_GONIOMETER_PARAMETERS=!select;
    	boolean select_ABERRATIONS_PARAMETERS=!select;
    	boolean select_FOCUSING_FIELD=!select;
    	if (select) {
    		GenericDialog gd = new GenericDialog("Select parameters to save");
    		gd.addMessage("===== Individual parameters ======");
        	gd.addCheckbox("MASTER_DEBUG_LEVEL", select_MASTER_DEBUG_LEVEL);
        	gd.addCheckbox("SHOW_AS_STACKS", select_SHOW_AS_STACKS);
        	gd.addCheckbox("FFT_SIZE",select_FFT_SIZE);
        	gd.addCheckbox("MAP_FFT_SIZE",select_MAP_FFT_SIZE);
        	gd.addCheckbox("GAUSS_WIDTH",select_GAUSS_WIDTH);
        	gd.addCheckbox("FFT_OVERLAP",select_FFT_OVERLAP);
        	gd.addCheckbox("PSF_SUBPIXEL",select_PSF_SUBPIXEL);
        	gd.addCheckbox("PSF_SAVE_FILE",select_PSF_SAVE_FILE);
        	gd.addCheckbox("THREADS_MAX",select_THREADS_MAX);
        	gd.addCheckbox("UPDATE_STATUS",select_UPDATE_STATUS);
    		gd.addMessage("===== Parameter classes ======");
        	gd.addCheckbox("SIMUL",select_SIMUL);
        	gd.addCheckbox("INTERPOLATE",select_INTERPOLATE);
        	gd.addCheckbox("INVERSE",select_INVERSE);
        	gd.addCheckbox("PSF_PARS",select_PSF_PARS);
        	gd.addCheckbox("OTF_FILTER",select_OTF_FILTER);
        	gd.addCheckbox("PATTERN_DETECT",select_PATTERN_DETECT);
        	gd.addCheckbox("COMPONENTS",select_COMPONENTS);
        	gd.addCheckbox("SHOW_RESULTS",select_SHOW_RESULTS);
        	gd.addCheckbox("MULTIFILE_PSF",select_MULTIFILE_PSF);
        	gd.addCheckbox("PROCESS_PARAMETERS",select_PROCESS_PARAMETERS);
        	gd.addCheckbox("DISTORTION",select_DISTORTION);
        	gd.addCheckbox("LASER_POINTERS",select_LASER_POINTERS);
        	gd.addCheckbox("FLATFIELD_PARAMETERS",select_FLATFIELD_PARAMETERS);
        	gd.addCheckbox("PATTERN_PARAMETERS",select_PATTERN_PARAMETERS);
        	gd.addCheckbox("LENS_DISTORTION_PARAMETERS",select_LENS_DISTORTION_PARAMETERS);
        	gd.addCheckbox("EYESIS_CAMERA_PARAMETERS",select_EYESIS_CAMERA_PARAMETERS);
        	gd.addCheckbox("LASERS",select_LASERS);
        	gd.addCheckbox("CAMERAS",select_CAMERAS);
        	gd.addCheckbox("DISTORTION_PROCESS_CONFIGURATION",select_DISTORTION_PROCESS_CONFIGURATION);
        	gd.addCheckbox("REFINE_PARAMETERS",select_REFINE_PARAMETERS);
        	gd.addCheckbox("FOCUS_MEASUREMENT_PARAMETERS",select_FOCUS_MEASUREMENT_PARAMETERS);
        	gd.addCheckbox("MOTORS_focusingHistory",select_MOTORS_focusingHistory);
        	gd.addCheckbox("GONIOMETER_PARAMETERS",select_GONIOMETER_PARAMETERS);
        	gd.addCheckbox("ABERRATIONS_PARAMETERS",select_ABERRATIONS_PARAMETERS);
        	gd.addCheckbox("FOCUSING_FIELD",select_FOCUSING_FIELD);
            WindowTools.addScrollBars(gd);
            gd.showDialog();
            if (gd.wasCanceled()) return;
        	select_MASTER_DEBUG_LEVEL=gd.getNextBoolean();
        	select_SHOW_AS_STACKS=gd.getNextBoolean();
        	select_FFT_SIZE=gd.getNextBoolean();
        	select_MAP_FFT_SIZE=gd.getNextBoolean();
        	select_GAUSS_WIDTH=gd.getNextBoolean();
        	select_FFT_OVERLAP=gd.getNextBoolean();
        	select_PSF_SUBPIXEL=gd.getNextBoolean();
        	select_PSF_SAVE_FILE=gd.getNextBoolean();
        	select_THREADS_MAX=gd.getNextBoolean();
        	select_UPDATE_STATUS=gd.getNextBoolean();

        	select_SIMUL=gd.getNextBoolean();
        	select_INTERPOLATE=gd.getNextBoolean();
        	select_INVERSE=gd.getNextBoolean();
        	select_PSF_PARS=gd.getNextBoolean();
        	select_OTF_FILTER=gd.getNextBoolean();
        	select_PATTERN_DETECT=gd.getNextBoolean();
        	select_COMPONENTS=gd.getNextBoolean();
        	select_SHOW_RESULTS=gd.getNextBoolean();
        	select_MULTIFILE_PSF=gd.getNextBoolean();
        	select_PROCESS_PARAMETERS=gd.getNextBoolean();
        	select_DISTORTION=gd.getNextBoolean();
        	select_LASER_POINTERS=gd.getNextBoolean();
        	select_FLATFIELD_PARAMETERS=gd.getNextBoolean();
        	select_PATTERN_PARAMETERS=gd.getNextBoolean();
        	select_LENS_DISTORTION_PARAMETERS=gd.getNextBoolean();
        	select_EYESIS_CAMERA_PARAMETERS=gd.getNextBoolean();
        	select_LASERS=gd.getNextBoolean();
        	select_CAMERAS=gd.getNextBoolean();
        	select_DISTORTION_PROCESS_CONFIGURATION=gd.getNextBoolean();
        	select_REFINE_PARAMETERS=gd.getNextBoolean();
        	select_FOCUS_MEASUREMENT_PARAMETERS=gd.getNextBoolean();
        	select_MOTORS_focusingHistory=gd.getNextBoolean();
        	select_GONIOMETER_PARAMETERS=gd.getNextBoolean();
        	select_ABERRATIONS_PARAMETERS=gd.getNextBoolean();
        	select_FOCUSING_FIELD=gd.getNextBoolean();
    	}
    	
       	if (select_MASTER_DEBUG_LEVEL) properties.setProperty("MASTER_DEBUG_LEVEL", MASTER_DEBUG_LEVEL+"");
       	if (select_SHOW_AS_STACKS) properties.setProperty("SHOW_AS_STACKS",     SHOW_AS_STACKS+"");//  Show debug images as stacks (false - individual)
       	if (select_FFT_SIZE) properties.setProperty("FFT_SIZE",           FFT_SIZE+"");
       	if (select_MAP_FFT_SIZE) properties.setProperty("MAP_FFT_SIZE",       MAP_FFT_SIZE+"");// used to find where grid covers the image
       	if (select_GAUSS_WIDTH) properties.setProperty("GAUSS_WIDTH",        GAUSS_WIDTH+""); // 0.4 (0 - use Hamming window)
       	if (select_FFT_OVERLAP) properties.setProperty("FFT_OVERLAP",        FFT_OVERLAP+""); // createPSFMap()
       	if (select_PSF_SUBPIXEL) properties.setProperty("PSF_SUBPIXEL",       PSF_SUBPIXEL+""); // sub-pixel decimation
       	if (select_PSF_SAVE_FILE) properties.setProperty("PSF_SAVE_FILE",      PSF_SAVE_FILE+"");// save PSF array to a multi-slice TIFF file
       	if (select_THREADS_MAX) properties.setProperty("THREADS_MAX",        THREADS_MAX+"");  // 100, testing multi-threading, limit maximal number of threads
       	if (select_UPDATE_STATUS) properties.setProperty("UPDATE_STATUS",      UPDATE_STATUS+""); // update ImageJ status info

       	if (select_SIMUL) SIMUL.setProperties(             "SIMUL.", properties);
        if (select_INTERPOLATE) INTERPOLATE.setProperties(       "INTERPOLATE.", properties);
        if (select_INVERSE) INVERSE.setProperties(           "INVERSE.", properties);
        if (select_PSF_PARS) PSF_PARS.setProperties(          "PSF_PARS.", properties);
        if (select_OTF_FILTER) OTF_FILTER.setProperties(        "OTF_FILTER.", properties);
        if (select_PATTERN_DETECT) PATTERN_DETECT.setProperties(    "PATTERN_DETECT.", properties);
        if (select_COMPONENTS) COMPONENTS.setProperties(        "COMPONENTS.", properties);
        if (select_SHOW_RESULTS) SHOW_RESULTS.setProperties(      "SHOW_RESULTS.", properties);
        if (select_MULTIFILE_PSF) MULTIFILE_PSF.setProperties(     "MULTIFILE_PSF.", properties);
        if (select_PROCESS_PARAMETERS) PROCESS_PARAMETERS.setProperties("PROCESS_PARAMETERS.", properties);
        if (select_DISTORTION) DISTORTION.setProperties(        "DISTORTION.", properties);
        if (select_LASER_POINTERS) LASER_POINTERS.setProperties("LASER.", properties);
        if (select_FLATFIELD_PARAMETERS) FLATFIELD_PARAMETERS.setProperties("FLATFIELD_PARAMETERS.", properties);
        if (select_PATTERN_PARAMETERS) PATTERN_PARAMETERS.setProperties  ("PATTERN_PARAMETERS.", properties);
        if (select_LENS_DISTORTION_PARAMETERS) LENS_DISTORTION_PARAMETERS.setProperties("LENS_DISTORTION_PARAMETERS.", properties);
        if (select_EYESIS_CAMERA_PARAMETERS) EYESIS_CAMERA_PARAMETERS.setProperties("EYESIS_CAMERA_PARAMETERS.", properties);
        if (select_LASERS) LASERS.setProperties("LASERS.", properties);
        if (select_CAMERAS) CAMERAS.setProperties("CAMERAS.", properties);
        if (select_DISTORTION_PROCESS_CONFIGURATION) DISTORTION_PROCESS_CONFIGURATION.setProperties("DISTORTION_PROCESS_CONFIGURATION.", properties);
        if (select_REFINE_PARAMETERS) REFINE_PARAMETERS.setProperties("REFINE_PARAMETERS.", properties);
        if (select_FOCUS_MEASUREMENT_PARAMETERS) FOCUS_MEASUREMENT_PARAMETERS.setProperties("FOCUS_MEASUREMENT_PARAMETERS.", properties);
        if (select_MOTORS_focusingHistory) MOTORS.focusingHistory.setProperties("FOCUSING_HISTORY.", properties);
        if (select_GONIOMETER_PARAMETERS) GONIOMETER_PARAMETERS.setProperties("GONIOMETER_PARAMETERS.", properties);
        if (select_ABERRATIONS_PARAMETERS) ABERRATIONS_PARAMETERS.setProperties("ABERRATIONS_PARAMETERS.", properties);
        if ((select_FOCUSING_FIELD) && (FOCUSING_FIELD!=null)) FOCUSING_FIELD.setProperties("FOCUSING_FIELD.", properties);
    	if (select) properties.remove("selected");
    }
/* ======================================================================== */
    public void getAllProperties(Properties properties){
       MASTER_DEBUG_LEVEL = Integer.parseInt(properties.getProperty("MASTER_DEBUG_LEVEL"));
   	   SHOW_AS_STACKS =     Boolean.parseBoolean(properties.getProperty("SHOW_AS_STACKS")); //  Show debug images as stacks (false - individual)
	   FFT_SIZE =           Integer.parseInt(properties.getProperty("FFT_SIZE"));
	   MAP_FFT_SIZE =       Integer.parseInt(properties.getProperty("MAP_FFT_SIZE")); // used to find where grid covers the image
	   GAUSS_WIDTH =        Double.parseDouble(properties.getProperty("GAUSS_WIDTH"));
	   FFT_OVERLAP =        Integer.parseInt(properties.getProperty("FFT_OVERLAP")); // createPSFMap()
	   PSF_SUBPIXEL =       Integer.parseInt(properties.getProperty("PSF_SUBPIXEL")); // sub-pixel decimation 
	   PSF_SAVE_FILE =      Boolean.parseBoolean(properties.getProperty("PSF_SAVE_FILE")); // save PSF array to a multi-slice TIFF file
	   THREADS_MAX =        Integer.parseInt(properties.getProperty("THREADS_MAX"));
       UPDATE_STATUS=       Boolean.parseBoolean(properties.getProperty("UPDATE_STATUS"));
       
       SIMUL.getProperties("SIMUL.", properties);
       INTERPOLATE.getProperties("INTERPOLATE.", properties);
       INVERSE.getProperties("INVERSE.", properties);
       PSF_PARS.getProperties("PSF_PARS.", properties);
       OTF_FILTER.getProperties("OTF_FILTER.", properties);
       PATTERN_DETECT.getProperties("PATTERN_DETECT.", properties);
       COMPONENTS.getProperties("COMPONENTS.", properties);
       SHOW_RESULTS.getProperties("SHOW_RESULTS.", properties);
       MULTIFILE_PSF.getProperties("MULTIFILE_PSF.", properties);
       PROCESS_PARAMETERS.getProperties("PROCESS_PARAMETERS.", properties);
       DISTORTION.getProperties(        "DISTORTION.", properties);
       LASER_POINTERS.getProperties("LASER.", properties);
       FLATFIELD_PARAMETERS.getProperties("FLATFIELD_PARAMETERS.", properties);
       PATTERN_PARAMETERS.getProperties("PATTERN_PARAMETERS.", properties);
       LENS_DISTORTION_PARAMETERS.getProperties("LENS_DISTORTION_PARAMETERS.", properties);
       EYESIS_CAMERA_PARAMETERS.getProperties("EYESIS_CAMERA_PARAMETERS.", properties);
       LASERS.getProperties("LASERS.", properties);
       CAMERAS.getProperties("CAMERAS.", properties);
       DISTORTION_PROCESS_CONFIGURATION.getProperties("DISTORTION_PROCESS_CONFIGURATION.", properties);
       REFINE_PARAMETERS.getProperties("REFINE_PARAMETERS.", properties);
       FOCUS_MEASUREMENT_PARAMETERS.getProperties("FOCUS_MEASUREMENT_PARAMETERS.", properties);
       MOTORS.focusingHistory.getProperties("FOCUSING_HISTORY.", properties);
       GONIOMETER_PARAMETERS.getProperties("GONIOMETER_PARAMETERS.", properties);
       ABERRATIONS_PARAMETERS.getProperties("ABERRATIONS_PARAMETERS.", properties);
       if (FOCUSING_FIELD!=null) FOCUSING_FIELD.getProperties("FOCUSING_FIELD.", properties,false); // false -> overwrite distortions center
    }
	
	  private String selectSourceDirectory(String defaultPath) {
		  return CalibrationFileManagement.selectDirectory(false, // save  
				  "Source calibration superdirectory selection (should have 1-1...3-3 subdirectories)", // title
				  "Select source calibration super-directory", // button
				  null, // filter
				  defaultPath);
	  }
	  private String selectPartialKernelsDirectory(String defaultPath) {
		  return  CalibrationFileManagement.selectDirectory(
				  true, // save
				  "Partial kernels superdirectory selection (should have 1-1...3-3 subdirectories)", // title
				  "Select partial kernels super-directory", // button
				  null, // filter
				  defaultPath);
	  }
	  private String selectKernelsDirectory(String defaultPath) {
		  return  CalibrationFileManagement.selectDirectory(
				  true, // save
				  "Result kernels directory selection",
				  "Select result kernels",
				  null, // filter
				  defaultPath);
	  }
	  private String selectFlatFieldDirectory(String defaultPath) {
		  return  CalibrationFileManagement.selectDirectory(
				  true, // save
				  "Flat-field results directory selection",
				  "Select flat-field results directory",
				  null, // filter
				  defaultPath);
	  }
	  private String [] selectSourceDirectories(String [] defaultPaths) {
		  return CalibrationFileManagement.selectDirectories(false, // save  
				  "Source directories for flat-field averaging selection", // title
				  "Select source directory(ies) for flat-field calculation", // button
				  null, // filter
				  defaultPaths);
	  }
	  
/* ======================================================================== */
/*	  public String selectDirectory(boolean save, String title, String button, FileFilter filter, String defaultPath) {
		  return selectDirectoryOrFile(save,true, title, button, filter,defaultPath);
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
*/
/* ======================================================================== */
/*
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
//		  for (int nFile=0;nFile<files.length;nFile++) filenames[nFile]= files[nFile].getName();
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
*/
/* ======================================================================== */
/*
	  class MultipleExtensionsFileFilter extends FileFilter {
		  protected String [] patterns; // case insensitive
		  protected String    description="JP4 files";
		  protected String    prefix=""; // case sensitive
		  
		  public MultipleExtensionsFileFilter (String prefix, String [] patterns,String description) {
			  this.prefix=     prefix;
			  this.description=description;
			  this.patterns=   patterns.clone();
		  }
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
			  if (!name.startsWith(this.prefix)) return false; // empty prefix OK
			  for (i=0;i<patterns.length;i++) {
				  if (name.toLowerCase().endsWith(patterns[i].toLowerCase())) return true;
			  }
			  return false; 
		  }
		  public String getDescription() {
			  return description;
		  }
	  }
*/
/* ======================================================================== */
	public boolean combinePSFKernels (
			EyesisAberrations.InterpolateParameters  interpolateParameters, // INTERPOLATE
			EyesisAberrations.MultiFilePSF           multiFilePSF ,         // MULTIFILE_PSF = new EyesisAberrations.MultiFilePSF(
			String []              filenames,
			String                 resultPath,
			showDoubleFloatArrays  sdfa_instance,        // SDFA_INSTANCE
			ImagePlus              imp_sel,
			boolean                saveResult,
			boolean                showResult,
			boolean                updateStatus,          // UPDATE_STATUS
			int                    thisDebugLevel

	){	
//		double [][][][] psfKernelMap=null;
		double [][][][][] kernelsElllipsePars = new double[filenames.length][][][][];
		int i;
		int nFile;
		Opener opener=new Opener();
		for (nFile=0;nFile<filenames.length;nFile++) {
			if (updateStatus) IJ.showStatus("Scanning file "+(nFile+1)+" (of "+(filenames.length)+"): "+filenames[nFile]);
			if (thisDebugLevel>1) System.out.println((nFile+1)+": "+filenames[nFile]);
			imp_sel=opener.openImage("", filenames[nFile]);  // or (path+filenames[nFile])
			kernelsElllipsePars[nFile]= kernelStackToEllipseCoefficients(
					imp_sel.getStack(), // Image stack, each slice consists of square kernels of one channel
					interpolateParameters.size, // size of each kernel (should be square)
					multiFilePSF.validateThreshold);               //      threshold) // to find ellipse
		}

		// Visualize the array as stacks
		int nFiles=kernelsElllipsePars.length;
		int kHeight=kernelsElllipsePars[0].length;
		int kWidth=kernelsElllipsePars[0][0].length;
		int kLength=kHeight*kWidth;
		int nChn=imp_sel.getStack().getSize();
		int numResults=7;
		double [][][][] c= new double[numResults][nChn][nFiles+1][kLength];
		double [][][] numVals=new double[numResults][nChn][kLength];
		int chn, tileY,tileX;
		boolean [] channels=new boolean[nChn];
		double a;
		if (thisDebugLevel>1) { 
			System.out.println("nFiles="+nFiles);
			System.out.println("kWidth="+kWidth);
			System.out.println("kHeight="+kHeight);
			System.out.println("nChn="+nChn);
		}
		Double D;
		int nOut;
		for (chn=0;chn<nChn;chn++) {
			channels[chn]=false;		
			for (nFile=0;nFile<nFiles;nFile++) for (tileY=0;tileY<kHeight;tileY++) for (tileX=0;tileX<kWidth;tileX++) {
				//   			  System.out.println("nChn="+nChn+" nFile="+nFile+" tileY="+tileY+" tileX="+tileX);
				if (kernelsElllipsePars[nFile][tileY][tileX][chn]!=null) {
					channels[chn]=true;
					c[0][chn][nFile+1][tileY*kWidth+tileX]=kernelsElllipsePars[nFile][tileY][tileX][chn][0];
					c[1][chn][nFile+1][tileY*kWidth+tileX]=kernelsElllipsePars[nFile][tileY][tileX][chn][1];
					c[2][chn][nFile+1][tileY*kWidth+tileX]=kernelsElllipsePars[nFile][tileY][tileX][chn][2];
					c[3][chn][nFile+1][tileY*kWidth+tileX]=kernelsElllipsePars[nFile][tileY][tileX][chn][3];
					c[4][chn][nFile+1][tileY*kWidth+tileX]=kernelsElllipsePars[nFile][tileY][tileX][chn][4];
					a=1/Math.sqrt(kernelsElllipsePars[nFile][tileY][tileX][chn][2]*kernelsElllipsePars[nFile][tileY][tileX][chn][3]-
							kernelsElllipsePars[nFile][tileY][tileX][chn][4]*kernelsElllipsePars[nFile][tileY][tileX][chn][4]/4);
					c[5][chn][nFile+1][tileY*kWidth+tileX]= Math.sqrt(a);
					c[6][chn][nFile+1][tileY*kWidth+tileX]=kernelsElllipsePars[nFile][tileY][tileX][chn][5];

				} else {
					c[0][chn][nFile+1][tileY*kWidth+tileX]=Double.NaN;
					c[1][chn][nFile+1][tileY*kWidth+tileX]=Double.NaN;
					c[2][chn][nFile+1][tileY*kWidth+tileX]=Double.NaN;
					c[3][chn][nFile+1][tileY*kWidth+tileX]=Double.NaN;
					c[4][chn][nFile+1][tileY*kWidth+tileX]=Double.NaN;
					c[5][chn][nFile+1][tileY*kWidth+tileX]=Double.NaN;
					c[6][chn][nFile+1][tileY*kWidth+tileX]=Double.NaN;
				}

			}
		}
		/* 
		 * Combine files - now just average all that are not NaN
		 */
		int [][] dirs={{-1,-1},{-1,0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1}};
		int yn,xn,index;
		// remove any tiles that are not OK in all channels
		double [][] weights=new double[nFiles+1][kLength];
		for (nFile=0;nFile<nFiles;nFile++) {
			for (i=0;i<kLength;i++){
				weights[nFile+1][i]=1.0;
				for (chn=0;chn<nChn;chn++) {
					D=c[0][chn][nFile+1][i];
					if (D.isNaN()) weights[nFile+1][i]=0.0;
				}
			}
			// Set weight to 0.5 if it has zero cells around        	
			for (tileY=0;tileY<kHeight;tileY++) for (tileX=0;tileX<kWidth;tileX++) {
				index=tileY*kWidth+tileX;
				if ( weights[nFile+1][index]>0.0){
					for (i=0;i<dirs.length;i++) {
						yn=tileY+dirs[i][1];
						xn=tileX+dirs[i][0];
						if ((yn>=0) && (yn<kHeight) && (xn>=0) && (xn<kWidth) && (weights[nFile+1][yn*kWidth+xn]==0.0)){
							weights[nFile+1][index]=0.5;
						}
					}
				}  
				weights[0][index]+=weights[nFile+1][index]; 
			}
		}
		//    	double [][][][] c= new double[numResults][nChn][nFiles+1][kLength];
		//     	double [][][] numVals=new double[numResults][nChn][kLength];
		if (multiFilePSF.validateShowEllipse) {
			for (chn=0;chn<nChn;chn++) {

				for (nOut=0;nOut<c.length;nOut++) {
					c[nOut][chn][0]=null;
					for (i=0;i<kLength;i++) {
						numVals[nOut][chn][i]=0.0;
					}
				}
				if (channels[chn]) {
					for (nOut=0;nOut<c.length;nOut++) {
						c[nOut][chn][0]=new double [kLength];
						for (nFile=0;nFile<nFiles;nFile++) {
							for (i=0;i<kLength;i++){
								D=c[nOut][chn][nFile+1][i];
								if (!D.isNaN()){
									numVals[nOut][chn][i]+=1.0;
									c[nOut][chn][0][i]+=D*weights[nFile+1][i]/weights[0][i];
								}
							}

						}
						for (i=0;i<kLength;i++){
							if (numVals[nOut][chn][i]==0.0 )c[nOut][chn][0][i]=Double.NaN;
							//    			  else c[nOut][chn][0][i]/=numVals[nOut][chn][i];
						}       	    	
					}
					sdfa_instance.showArrays(c[0][chn], kWidth, kHeight,  true, "x-shift-"+chn);
					sdfa_instance.showArrays(c[1][chn], kWidth, kHeight,  true, "y-shift-"+chn);
					sdfa_instance.showArrays(c[5][chn],kWidth, kHeight,  true, "radius-"+chn);
					if (thisDebugLevel>1) {
						sdfa_instance.showArrays(c[2][chn], kWidth, kHeight,  true, "x2-"+chn);
						sdfa_instance.showArrays(c[3][chn], kWidth, kHeight,  true, "y2-"+chn);
						sdfa_instance.showArrays(c[4][chn], kWidth, kHeight,  true, "xy-"+chn);
						sdfa_instance.showArrays(c[6][chn], kWidth, kHeight,  true, "area-"+chn);
					}  
				}

			}
		}
		if (multiFilePSF.showWeights) sdfa_instance.showArrays(weights, kWidth, kHeight,  true, "weights");
		//    	double [][] weights=new double[nFiles+1][kLength];
		for (i=0;i<kLength;i++) weights[0][i]=0.0;
		PSF_KERNEL_MAP=new double [kHeight][kWidth][nChn][];
		for (tileY=0;tileY<kHeight;tileY++) for (tileX=0;tileX<kWidth;tileX++) for (chn=0;chn<nChn;chn++){
			PSF_KERNEL_MAP[tileY][tileX][chn]=null;
		}
		String [] originalSliceLabels=null;
		for (nFile=0;nFile<nFiles;nFile++) {
			if (updateStatus) IJ.showStatus("Accumulating file "+(nFile+1)+" (of "+nFiles+"): "+filenames[nFile]);
			if (thisDebugLevel>1) System.out.println("Accumulating file "+nFile+": "+filenames[nFile]);
			imp_sel=opener.openImage("", filenames[nFile]);  // or (path+filenames[nFile])
			if (originalSliceLabels==null) {
				originalSliceLabels=imp_sel.getStack().getSliceLabels();
			}
			accumulatePartialKernelStack(
					imp_sel.getStack(), // Image stack with partial array of kernels, each slice consists of square kernels of one channel
					interpolateParameters.size, // size of each kernel (should be square)
					weights[nFile+1], // weights of the kernel tiles in the current stack
					weights[0]);// weights of the kernel tiles already accumulated (will be updated)

		}
// optionally fill in blanks from nearest neighbors
		int filledMissing=0;
		//Finalize accumulated kernels - transform them from frequency to space domain
		inverseTransformKernels();
// should be done after inversion, because filled in kernels are just pointers to original ones		
		if (multiFilePSF.fillMissing) filledMissing=fillMissingKernels (PSF_KERNEL_MAP);
		int numMissing=0;
		ImageStack mergedStack= mergeKernelsToStack(PSF_KERNEL_MAP,originalSliceLabels);
		System.out.println("mergedStack.getSize()= "+mergedStack.getSize());
		System.out.println("mergedStack.getWidth()= "+mergedStack.getWidth()  );
		System.out.println("mergedStack.getHeight()= "+mergedStack.getHeight()  );
		System.out.println("PSF_KERNEL_MAP.length= "+PSF_KERNEL_MAP.length  );
		System.out.println("PSF_KERNEL_MAP[0].length= "+PSF_KERNEL_MAP[0].length  );
		System.out.println("mergedStack= "+((mergedStack==null)?"null":"not null"));

		if (mergedStack.getSize()==0) {
			System.out.println("*** Error - result is empty");
			return false;
		}

		for (tileY=0;tileY<kHeight;tileY++) for (tileX=0;tileX<kWidth;tileX++) if ((PSF_KERNEL_MAP[tileY][tileX]==null) || (PSF_KERNEL_MAP[tileY][tileX][0]==null)) numMissing++;
        ImagePlus imp_psf = new ImagePlus(resultPath, mergedStack);
        if (showResult) {
        	imp_psf.getProcessor().resetMinAndMax();
        	imp_psf.show();
        }
		if (saveResult) {
			if (numMissing==0) {
			  if (thisDebugLevel>1) System.out.println("Saving result to "+resultPath);
			  FileSaver fs=new FileSaver(imp_psf);
			  fs.saveAsTiffStack(resultPath);
			  if (multiFilePSF.fillMissing && (filledMissing>0)) {
					System.out.println("*** Warning "+filledMissing+" kernel tiles are missing from the results (insufficient overlap) \n"+
					"You may disable filling missing kernels from neighbors in Conf. Multifile");
					IJ.showMessage("Warning",filledMissing+" kernel tiles were missing from the results\n"+
							"(i.e.insufficient overlap) and filled from neighbors, it is OK only for the fisheye lens.\n"+
					        "You may disable filling missing kernels from neighbors in Conf. Multifile");
			  }
			  return true;
			} else {
				System.out.println("*** Error "+numMissing+" kernel tiles are missing from the results (insufficient overlap), result is not saved\n"+
				"You may enable filling missing kernels from neighbors if it is a fisheye lens (in Conf. Multifile)");

				IJ.showMessage("Error",numMissing+" kernel tiles are missing from the results\n (insufficient overlap), result file is not saved\n"+
				"You may enable filling missing kernels from neighbors if it is a fisheye lens (in Conf. Multifile)");

				if (!showResult) { // not yet shown
		        	imp_psf.getProcessor().resetMinAndMax();
		        	imp_psf.show();
				}
				return false;
			}
		}
		if (numMissing>0) {
			System.out.println("*** Error "+numMissing+" kernel tiles are missing from the results (insufficient overlap) \n"+
					"You may enable filling missing kernels from neighbors if it is a fisheye lens (in Conf. Multifile)");
			IJ.showMessage("Error",numMissing+" kernel tiles are missing from the results\n (insufficient overlap)\n"+
					"You may enable filling missing kernels from neighbors if it is a fisheye lens (in Conf. Multifile)");
			return false;

		} else if (multiFilePSF.fillMissing && (filledMissing>0)) {
			System.out.println("*** Warning "+filledMissing+" kernel tiles are missing from the results (insufficient overlap) \n"+
			"You may disable filling missing kernels from neighbors in Conf. Multifile");
			IJ.showMessage("Warning",filledMissing+" kernel tiles were missing from the results\n"+
					"(i.e.insufficient overlap) and filled from neighbors, it is OK only for the fisheye lens.\n"+
			        "You may disable filling missing kernels from neighbors in Conf. Multifile");
		}
		return true;
	}
	private int fillMissingKernels(double [][][][] kernels){
		int [][] dirs={{-1,0},{1,0},{0,-1},{0,1}};
		List <Integer> kernelList=new ArrayList<Integer>(100);
		Integer Index;
		kernelList.clear();
		int tileY,tileX,newTileY,newTileX,nDir,numMissing=0;
		int width= kernels[0].length;
		int height=kernels.length;
		for (tileY=0;tileY<height;tileY++) for (tileX=0;tileX<width;tileX++) {
			if ((kernels[tileY][tileX]==null) || (kernels[tileY][tileX][0]==null)) {
				Index=tileY*width+tileX;
				for (nDir=0;nDir<dirs.length;nDir++) {
					newTileX=tileX+dirs[nDir][0];
					newTileY=tileY+dirs[nDir][1];
					if ((newTileX>=0) && (newTileY>=0) && (newTileX<width) && (newTileY<height) &&
							(kernels[newTileY][newTileX]!=null) && (kernels[newTileY][newTileX][0]!=null) ) {
						kernelList.add(Index);
					}
				}				
				numMissing++;
			}
		}
		System.out.println("fillMissingKernels: numMissing="+numMissing);
		System.out.println("fillMissingKernels: kernelList.size()="+kernelList.size());

		while (kernelList.size()>0) {
			Index=kernelList.get(0);
			kernelList.remove(0);
			tileY=Index/width;
			tileX=Index%width;
			if ((kernels[tileY][tileX]==null) || (kernels[tileY][tileX][0]==null)) {// may be duplicates (added several times)
//TODO: - change order of directions?
				for (nDir=0;nDir<dirs.length;nDir++) {
					newTileX=tileX+dirs[nDir][0];
					newTileY=tileY+dirs[nDir][1];
					if ((newTileX>=0) && (newTileY>=0) && (newTileX<width) && (newTileY<height)) {
						if ((kernels[newTileY][newTileX]==null) || (kernels[newTileY][newTileX][0]==null)) {
							Index=newTileY*width+newTileX;
							kernelList.add(Index);
						} else if ((kernels[tileY][tileX]==null) || (kernels[tileY][tileX][0]==null)) { // may be already added
// need to copy - they will be subject to reverse fht	?						
							kernels[tileY][tileX]=kernels[newTileY][newTileX];
							System.out.println("fillMissingKernels: filled "+tileX+"/"+tileY);
						}
					}
				}
			}
		}
		return numMissing;
		
	}
/*
 * 				while (listIndex<pixelList.size() ) {
					Index=pixelList.get(listIndex++);
					pixelList.remove(i);

				pixelList.clear();
				pixelList.add (Index);
	PSF_KERNEL_MAP=new double [kHeight][kWidth][nChn][];
	for (tileY=0;tileY<kHeight;tileY++) for (tileX=0;tileX<kWidth;tileX++) for (chn=0;chn<nChn;chn++){
		PSF_KERNEL_MAP[tileY][tileX][chn]=null;
	}
*/	
	
/* ======================================================================== */
// returns [*][0][*] - only one pair src/dest per channel
/* ======================================================================== */
	// returns [*][0][*] - only one pair src/dest per channel
	public String [][][] prepareInvertGaussianKernelsList(ProcessCalibrationFilesParameters processCalibrationFilesParameters,
			String prefix) { // one of rpsfPrefix and gaussianPrefix
		String [][][] fileNames = new String [processCalibrationFilesParameters.subdirNames.length][][];
		int dirNum;

		File dir,destDir,file;
		String dirName,destDirName,fileName;
		dirName=processCalibrationFilesParameters.partialKernelsSuperDirectory+Prefs.getFileSeparator()+processCalibrationFilesParameters.combinedSubDirectory;
		dir= new File (dirName);
		if (!dir.exists()) {
			IJ.showMessage("Error","Combined PSF kernels directory "+dirName+" does not exist");
			return null;
		}
		destDirName=processCalibrationFilesParameters.kernelsDirectory;
		destDir= new File (destDirName);
		if (!destDir.exists()){
			if (!destDir.mkdirs()) {
				IJ.showMessage("Error","Failed to create result kernels directory "+destDirName);
				return null;
			}
		}
		for (dirNum=0;dirNum<fileNames.length;dirNum++) {
			fileNames[dirNum]=null;
			// process only selected channels (or all)			
			if (processCalibrationFilesParameters.processAllChannels || processCalibrationFilesParameters.processChannels[dirNum]) {
				fileName=dirName+Prefs.getFileSeparator()+
				processCalibrationFilesParameters.psfInterpoaltedPrefix+
				processCalibrationFilesParameters.suffixes[dirNum]+
				"."+processCalibrationFilesParameters.kernelFileExtension;
				file=new File(fileName);
				if (!file.exists()) {
					if (DEBUG_LEVEL>1) System.out.println("PSF kernel file "+fileName+" does not exist.");
					continue;
				}
				fileNames[dirNum]=new String[1][2];
				fileNames[dirNum][0][0]=fileName;
				fileNames[dirNum][0][1]=destDirName+Prefs.getFileSeparator()+
				prefix+
				processCalibrationFilesParameters.suffixes[dirNum]+
				"."+processCalibrationFilesParameters.kernelFileExtension;
			}
		}

		return fileNames;
	}
/* ======================================================================== */
	// returns [*][0][*] - only one pair src/dest per channel
	public String [][][] prepareInterpolateKernelsList(ProcessCalibrationFilesParameters processCalibrationFilesParameters) {
		String [][][] fileNames = new String [processCalibrationFilesParameters.subdirNames.length][][];
		int dirNum;
		File dir,file;
		String dirName,fileName;
		dirName=processCalibrationFilesParameters.partialKernelsSuperDirectory+Prefs.getFileSeparator()+processCalibrationFilesParameters.combinedSubDirectory;
		dir= new File (dirName);
		if (!dir.exists()) {
			IJ.showMessage("Error","Combined PSF kernels directory "+dirName+" does not exist");
            return null;
		}
		if (DEBUG_LEVEL>1) System.out.println("fileNames.length= "+fileNames.length);
		for (dirNum=0;dirNum<fileNames.length;dirNum++) {
			fileNames[dirNum]=null;
// process only selected channels (or all)			
			if (processCalibrationFilesParameters.processAllChannels || processCalibrationFilesParameters.processChannels[dirNum]) {
				fileName=dirName+Prefs.getFileSeparator()+
			       processCalibrationFilesParameters.psfRawPrefix+
			       processCalibrationFilesParameters.suffixes[dirNum]+
			       "."+processCalibrationFilesParameters.kernelFileExtension;
				file=new File(fileName);
				if (!file.exists()) {
					if (DEBUG_LEVEL>1) System.out.println("PSF kernel file "+fileName+" does not exist.");
		            continue;
				}
				fileNames[dirNum]=new String[1][2];
				fileNames[dirNum][0][0]=fileName;
				fileNames[dirNum][0][1]=dirName+Prefs.getFileSeparator()+
			       processCalibrationFilesParameters.psfInterpoaltedPrefix+
			       processCalibrationFilesParameters.suffixes[dirNum]+
			       "."+processCalibrationFilesParameters.kernelFileExtension;
				if (DEBUG_LEVEL>1) System.out.println("fileNames[dirNum][0][0]= "+fileNames[dirNum][0][0]);
				if (DEBUG_LEVEL>1) System.out.println("fileNames[dirNum][0][1]= "+fileNames[dirNum][0][1]);

			}
		}
		
		return fileNames;
	}
/* ======================================================================== */
	public String [][][] preparePartialKernelsFilesList(ProcessCalibrationFilesParameters processCalibrationFilesParameters) {
		String [][][] fileNames = new String [processCalibrationFilesParameters.subdirNames.length][][];
		int dirNum,fileNum,numFiles,i;
		File dir,destDir;
		File [] fileList;
		String dirName,destDirName;
		String [] names;
		dirName=processCalibrationFilesParameters.partialKernelsSuperDirectory;
		dir= new File (dirName);
		if (!dir.exists()) {
			IJ.showMessage("Error","Partial PSF kernels directory "+dirName+" does not exist");
            return null;
		}
		destDirName=dirName+Prefs.getFileSeparator()+processCalibrationFilesParameters.combinedSubDirectory;
		destDir= new File (destDirName);
		if (!destDir.exists()){
			if (!destDir.mkdirs()) {
				IJ.showMessage("Error","Failed to create combined PSF kernels directory "+destDirName);
	            return null;
			}
		}
		for (dirNum=0;dirNum<fileNames.length;dirNum++) {
			fileNames[dirNum]=null;
// process only selected channels (or all)			
			if (processCalibrationFilesParameters.processAllChannels || processCalibrationFilesParameters.processChannels[dirNum]) {
				dirName=processCalibrationFilesParameters.partialKernelsSuperDirectory+Prefs.getFileSeparator()+processCalibrationFilesParameters.subdirNames[dirNum];
				dir= new File (dirName);
				if (!dir.exists()) {
					if (DEBUG_LEVEL>1) System.out.println("Partial PSF kerneldirectory "+dirName+" does not exist.");
					continue;
				}
// get list of all tiff files starting with specified prefix 	
				fileList=dir.listFiles(new Filter(processCalibrationFilesParameters.kernelFilePrefix,"."+processCalibrationFilesParameters.kernelFileExtension));
				if (DEBUG_LEVEL>1) System.out.println("Calibration files directory has "+fileList.length+" files.");
				if (fileList.length==0) continue;

				names = new String[fileList.length];
				numFiles=names.length;
				for (fileNum=0;fileNum<names.length;fileNum++){
					names[fileNum]=fileList[fileNum].getName(); // here names include extensions, keepOld is not used
					if (DEBUG_LEVEL>1) System.out.println(fileNum+": "+names[fileNum]);

				}
				if (numFiles==0) {
					if (DEBUG_LEVEL>1) System.out.println("No files to be processed in "+dirName);
					continue;
				}
				// manually select files to process				
				if (processCalibrationFilesParameters.selectFiles) {
					if (!selectFilesToProcess (names,
							null, // extension is already included
							"Selection "+processCalibrationFilesParameters.subdirNames[dirNum],
							"Select partial PSF kernels to process in "+dirName)){
						if (DEBUG_LEVEL>1) System.out.println("Operation canceled");
						return null; // operation canceled
					}
				}
				numFiles=0;
				for (i=0;i<names.length;i++) if (names[i]!=null) numFiles++;
				if (numFiles==0) {
					if (DEBUG_LEVEL>1) System.out.println("No files to be processed in "+dirName);
					continue;
				}
				// create source path, remove null-ed files (destination path will be in the first "pair"				
				fileNames[dirNum]=new String [numFiles][2];
				fileNum=0;
				for (i=0;i<names.length;i++) if (names[i]!=null) {
					fileNames[dirNum][fileNum  ][0]=    dirName+Prefs.getFileSeparator()+names[i];
					fileNames[dirNum][fileNum++][1]=    null;
				}
				if (numFiles>0) {
					fileNames[dirNum][0][1]=destDirName+Prefs.getFileSeparator()+
					processCalibrationFilesParameters.psfRawPrefix+
					processCalibrationFilesParameters.suffixes[dirNum]+
					"."+processCalibrationFilesParameters.kernelFileExtension;
				}
				
			}			
		}		
		return fileNames;
	}
/* ======================================================================== */
	public String [][][] prepareCalibrationFilesList(ProcessCalibrationFilesParameters processCalibrationFilesParameters) {
		String [][][] fileNames = new String [processCalibrationFilesParameters.subdirNames.length][][];
		int dirNum,fileNum,numFiles,i;
		File dir,destDir,file;
		File [] fileList;
		String dirName,destDirName,destName;
		String [] names;
		dirName=processCalibrationFilesParameters.sourceSuperDirectory;
		dir= new File (dirName);
		if (!dir.exists()) {
			IJ.showMessage("Error","Source calibration files directory "+dirName+" does not exist");
            return null;
		}
		destDirName=processCalibrationFilesParameters.partialKernelsSuperDirectory;
		destDir= new File (destDirName);
		if (!destDir.exists()){
			if (!destDir.mkdirs()) {
				IJ.showMessage("Error","Failed to create partial PSF kernels directory "+destDirName);
	            return null;
			}
		}

// see if result super-directory exists, if not create it		
		for (dirNum=0;dirNum<fileNames.length;dirNum++) {
			fileNames[dirNum]=null;
			if (processCalibrationFilesParameters.processAllChannels || processCalibrationFilesParameters.processChannels[dirNum]) {
				dirName=processCalibrationFilesParameters.sourceSuperDirectory+Prefs.getFileSeparator()+processCalibrationFilesParameters.subdirNames[dirNum];
				dir= new File (dirName);
				if (!dir.exists()) {
					if (DEBUG_LEVEL>1) System.out.println("Calibration files directory "+dirName+" does not exist.");
					continue;
				}
// get list of all jp4 files	
				fileList=dir.listFiles(new Filter("."+processCalibrationFilesParameters.sourceFileExtension));
				if (DEBUG_LEVEL>1) System.out.println("Calibration files directory has "+fileList.length+" files.");
				if (fileList.length==0) continue;
//				see if result directory exists, if not - create it
				destDirName=processCalibrationFilesParameters.partialKernelsSuperDirectory+Prefs.getFileSeparator()+processCalibrationFilesParameters.subdirNames[dirNum];
				destDir= new File (destDirName);
				if (!destDir.exists()){
					if (!destDir.mkdirs()) {
						IJ.showMessage("Error","Failed to create partial PSF kernels directory "+destDirName);
			            return null;
					}
				}
// for each file in the list, see if the result file exists and keepOld is true.
				names = new String[fileList.length];
				numFiles=names.length;
				for (fileNum=0;fileNum<names.length;fileNum++){
					names[fileNum]=fileList[fileNum].getName();
// remove extension (directory is already removed by getName() itself
					names[fileNum]=names[fileNum].substring(0, names[fileNum].length()-processCalibrationFilesParameters.sourceFileExtension.length()-1);
					if (processCalibrationFilesParameters.keepOld) {
						destName=destDirName+Prefs.getFileSeparator()+processCalibrationFilesParameters.kernelFilePrefix+names[fileNum]+"."+processCalibrationFilesParameters.kernelFileExtension;
						file=new File(destName);
						if (file.exists()) {
							if (DEBUG_LEVEL>1) System.out.println("Skipping existent PSF file re-calculation: "+destName);
							names[fileNum]=null;
							numFiles--;
						}
					}
				}
				if (numFiles==0) {
					if (DEBUG_LEVEL>1) System.out.println("No files to be processed in "+dirName);
					continue;
				}
// manually select files to process				
				if (processCalibrationFilesParameters.selectFiles) {
					if (!selectFilesToProcess (names,
							processCalibrationFilesParameters.sourceFileExtension,
							"Selection "+processCalibrationFilesParameters.subdirNames[dirNum],
							"Select calibration files to process in "+dirName)){
					  if (DEBUG_LEVEL>1) System.out.println("Operation canceled");
					  return null; // operation canceled
					}
					numFiles=0;
					for (i=0;i<names.length;i++) if (names[i]!=null) numFiles++;
					if (numFiles==0) {
						if (DEBUG_LEVEL>1) System.out.println("No files to be processed in "+dirName);
						continue;
					}
				}
// create source/destination path pairs, remove null-ed files				
			    fileNames[dirNum]=new String [numFiles][2];
				fileNum=0;
				for (i=0;i<names.length;i++) if (names[i]!=null) {
					fileNames[dirNum][fileNum  ][0]=    dirName+Prefs.getFileSeparator()+names[i]+"."+processCalibrationFilesParameters.sourceFileExtension;
					fileNames[dirNum][fileNum++][1]=destDirName+Prefs.getFileSeparator()+processCalibrationFilesParameters.kernelFilePrefix+names[i]+"."+processCalibrationFilesParameters.kernelFileExtension;
				}
			
			}
		}
		return fileNames;
	}
	
	
	private boolean selectFilesToProcess (String []  names,
			                           String extension, // does not include "."
			                           String     title,
                                       String     message) {
		GenericDialog gd = new GenericDialog(title);
		gd.addMessage(message);
		int i;
		for (i=0;i<names.length;i++) if (names[i]!=null) gd.addCheckbox(names[i]+(((extension==null) || extension.equals(""))?"":("."+extension)), false); 
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		for (i=0;i<names.length;i++) if (names[i]!=null) if (!gd.getNextBoolean()) names[i]=null; 
		return true;
	}
	
	class Filter implements FilenameFilter {
		  protected String pattern;
		  protected String prefix=null;
		  public Filter (String str) {
		    pattern = str;
		  }
		  public Filter (String pref, String str) {
			    prefix=   pref;
			    pattern = str;
		  }
		  public boolean accept (File dir, String name) {
			if (prefix!=null) return (name.startsWith(prefix)) && (name.toLowerCase().endsWith(pattern.toLowerCase()));
		    return name.toLowerCase().endsWith(pattern.toLowerCase());
		  }
	}

	
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */


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



/* ======================================================================== */

	public ImageStack  generateGaussianStackFromDirect(ImageStack PSFStack, // stack of 3 32-bit (float) images, made of square kernels
			EyesisAberrations.InverseParameters inverseParameters, // size (side of square) of direct PSF kernel
			//		  int dSize, // size (side of square) of direct PSF kernel
			//		  int size, // size (side of square) of the Gaussian kernel
			//		  double psf_cutoff_energy,  // OTF_cutoff_energy
			//		  double [] sigmas, // array of sigmas in the center, matching stacks sequence. Null if no blurring is needed
			boolean updateStatus){  // update status info
		if (PSFStack==null) return null;
		double [] sigmas={inverseParameters.gaussianSigmaIndividual,inverseParameters.gaussianSigmaIndividual,inverseParameters.gaussianSigmaChecker};

		int inWidth=PSFStack.getWidth();
		int inHeight=PSFStack.getHeight();
		int tilesX=inWidth/inverseParameters.dSize;
		int tilesY=inHeight/inverseParameters.dSize;
		int nChn=PSFStack.getSize();
		float [][] outPixels=new float[nChn][tilesX*inverseParameters.rSize*tilesY*inverseParameters.rSize];
		float [] pixels;
		double [] kernel= new double[inverseParameters.dSize*inverseParameters.dSize];
		int  [][]selection;
		double [] ellipse_coeff;
		int chn,tileY,tileX;
		double nSigma2;
		double x,y,x2,y2;
		int       length=   inverseParameters.rSize*inverseParameters.rSize;
		double [] gaussian= new double[length];
		double [] gaussianX=new double[inverseParameters.rSize];
		double [] gaussianY=new double[inverseParameters.rSize];
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
					selection=    findClusterOnPSF(kernel,  inverseParameters.psfCutoffEnergy, "");
					ellipse_coeff=findEllipseOnPSF(kernel,  selection, ""); // coefficients for direct PSF, for rPSF [0] and [1] need to be opposite size
					nSigma2=(4*sigmas[chn])*(4*sigmas[chn]);
					k=(sigmas[chn]<minsigma)?(0.5/(minsigma*minsigma)):(0.5/(sigmas[chn]*sigmas[chn]));
					for ( i=0; i<inverseParameters.rSize;i++) {
						x=i-inverseParameters.rSize/2+ellipse_coeff[0];
						x2=x*x;
						if (x2>nSigma2) gaussianX[i]=0.0;
						else gaussianX[i]=Math.exp(-x2*k);
						y=i-inverseParameters.rSize/2+ellipse_coeff[1];
						y2=y*y;
						if (y2>nSigma2) gaussianY[i]=0.0;
						else gaussianY[i]=Math.exp(-y2*k);
					}
					sum=0.0;
					for ( i=0; i<inverseParameters.rSize;i++) for (j=0;j<inverseParameters.rSize;j++) {
						d=gaussianX[j]*gaussianY[i];
						sum+=d;
						gaussian[i*inverseParameters.rSize+j]=d;
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
		ImageStack outStack=new ImageStack(tilesX*inverseParameters.rSize,tilesY*inverseParameters.rSize);
		for (chn=0;chn<nChn;chn++) {
			outStack.addSlice(PSFStack.getSliceLabel(chn+1), outPixels[chn]);
		}
		return outStack;
	}



/* ======================================================================== */

	//getSliceLabels()
	// Will build global PSF_KERNEL_MAP (each [][][]element should be set to null?
	// kernels are supposed to be normalized?
	public void accumulatePartialKernelStack(
			ImageStack   kernelStack, // Image stack with partial array of kernels, each slice consists of square kernels of one channel
			int                 size, // size of each kernel (should be square)
			double []   theseWeights, // weights of the kernel tiles in the current stack
			double []   accumWeights){// weights of the kernel tiles already accumulated (will be updated)
		DoubleFHT fht_instance =new DoubleFHT(); // provide DoubleFHT instance to save on initializations (or null)
		int tilesX=kernelStack.getWidth()/size;
		int tilesY=kernelStack.getHeight()/size;
		int nChn=kernelStack.getSize();
		int tileY,tileX, chn; //,subTileY,subTileX;
		float [] pixels;
		int length=size*size;
		double [] kernel=new double[length];
		int index;
		boolean lastChn;
		for (chn=0;chn<nChn;chn++) {
			pixels=(float[]) kernelStack.getPixels(chn+1);
			lastChn= (chn==(nChn-1));
			for (tileY=0;tileY<tilesY;tileY++) for (tileX=0;tileX<tilesX;tileX++) {
				index=tileY*tilesX+tileX;
				if (theseWeights[index]>0.0){
					extractOneKernel(
							pixels, //  array of combined square kernels, each 
							kernel, // will be filled, should have correct size before call
							tilesX, // number of kernels in a row
							tileX, // horizontal number of kernel to extract
							tileY); // vertical number of kernel to extract
					// convert to frequency domain (interpolation is for FHT)
					fht_instance.swapQuadrants(kernel);
					fht_instance.transform(    kernel);
					if (!(accumWeights[index]>0.0)) { // nothing yet in this tile
						PSF_KERNEL_MAP[tileY][tileX][chn]=kernel.clone();
						if (lastChn) accumWeights[index]=theseWeights[index];
					} else { // "accumulate" - interpolate between existent and new kernel, using/updating weights
						if (DEBUG_LEVEL>5) {
							System.out.println("tileY="+tileY+" tileX="+tileX+" chn= "+chn);
							System.out.println("PSF_KERNEL_MAP[tileY][tileX][chn].length= "+PSF_KERNEL_MAP[tileY][tileX][chn].length);
							System.out.println("kernel.length= "+kernel.length);
						}

						kernel=fht_instance.interpolateFHT (
								PSF_KERNEL_MAP[tileY][tileX][chn],    // first FHT array
								kernel,    // second FHT array
								theseWeights[index]/accumWeights[index]);    //interpolation ratio - 0.0 - fht0, 1.0 - fht1
						if (lastChn) accumWeights[index]+=theseWeights[index];
					}
				}
			}
		}
	}

/* ======================================================================== */
	//Finalize accumulated kernels - transfrom them from frequency to space domain
	public void inverseTransformKernels(){
		DoubleFHT fht_instance =new DoubleFHT(); // provide DoubleFHT instance to save on initializations (or null)
		int tilesX=PSF_KERNEL_MAP[0].length;
		int tilesY=PSF_KERNEL_MAP.length;
		int tileY,tileX, chn; //,subTileY,subTileX;
		for (tileY=0;tileY<tilesY;tileY++) for (tileX=0;tileX<tilesX;tileX++) for (chn=0; chn<PSF_KERNEL_MAP[tileY][tileX].length; chn++){
			if (PSF_KERNEL_MAP[tileY][tileX][chn]!=null) {
				fht_instance.inverseTransform(PSF_KERNEL_MAP[tileY][tileX][chn]);
				fht_instance.swapQuadrants   (PSF_KERNEL_MAP[tileY][tileX][chn]);
			}
		}
	}



/* ======================================================================== */

	public double [][][][] kernelStackToEllipseCoefficients(
			ImageStack kernelStack, // Image stack, each slice consists of square kernels of one channel
			int               size, // size of each kernel (should be square)
			double       threshold) // to find ellipse
	// update status info
	{
		//	  DoubleFHT fht_instance =new DoubleFHT(); // provide DoubleFHT instance to save on initializations (or null)
		if (kernelStack==null) return null;
		int tilesX=kernelStack.getWidth()/size;
		int tilesY=kernelStack.getHeight()/size;
		int nChn=kernelStack.getSize();
		float [] pixels;
		int i,j;
		int tileY,tileX, chn; //,subTileY,subTileX;
		double [][][][] ellipseCoeffs=new double [tilesY][tilesX][nChn][];
		int length=size*size;
		double [] kernel=new double[length]; 
		double max;
		int  [][]selection;
		double [] ec;
		int l;
		for (chn=0;chn<nChn;chn++) {
			pixels=(float[]) kernelStack.getPixels(chn+1);
			for (tileY=0;tileY<tilesY;tileY++) for (tileX=0;tileX<tilesX;tileX++) {
				extractOneKernel(
						pixels, //  array of combined square kernels, each 
						kernel, // will be filled, should have correct size before call
						tilesX, // number of kernels in a row
						tileX, // horizontal number of kernel to extract
						tileY); // vertical number of kernel to extract
				max=0.0;
				for (i=0;i<length;i++) if (max<kernel[i]) max=kernel[i];
				if (max<=0.0) ellipseCoeffs[tileY][tileX][chn]=null;
				else {
					selection= findClusterOnPSF(
							kernel, // PSF function, square array
							threshold, // fraction of energy in the pixels to be used
					"");
					//				  ellipseCoeffs[tileY][tileX][chn]=findEllipseOnPSF(kernel,  selection,   "");
					ec=findEllipseOnPSF(kernel,  selection,   "");
					l=ec.length;
					ellipseCoeffs[tileY][tileX][chn]=new double[l+1];
					for (i=0;i<ec.length;i++) ellipseCoeffs[tileY][tileX][chn][i]=ec[i];
					ellipseCoeffs[tileY][tileX][chn][l]=0;
					for (i=0;i<selection.length;i++) for (j=0;j<selection[0].length;j++) ellipseCoeffs[tileY][tileX][chn][l]+=selection[i][j];
				}
			}
		}
		return ellipseCoeffs;

	} 
/* ======================================================================== */

/* interpolate kernels minimizing memory image - use directly the image stack (32-bit, float) with kernels.
  Add kernels around by either replication or extrapolation to compensate for "margins" in the original; kernels */
	public ImageStack interpolateKernelStack(
			ImageStack kernelStack, // Image stack, each slice consists of square kernels of one channel
			EyesisAberrations.InterpolateParameters interpolateParameters,
			//		  int               size, // size of each kernel (should be square)
			//		  int             subdiv, // number of subdivisions form input to output
			//		  int            addLeft, // add this number of kernel columns to the output on the left of the existent/interpolated
			//		  int             addTop, // add this number of kernel rows to the output above the existent/interpolated
			//		  int           addRight, // add this number of kernel columns to the output on the right of the existent/interpolated
			//		  int          addBottom, // add this number of kernel rows to the output below the existent/interpolated
			//		  double     extrapolate, // 0 - duplicate, 1.0 - extrapolate outside of the known kernels
			boolean   updateStatus) // update status info
	{
		DoubleFHT fht_instance =new DoubleFHT(); // provide DoubleFHT instance to save on initializations (or null)
		if (kernelStack==null) return null;
		int inTilesX=kernelStack.getWidth()/interpolateParameters.size;
		int inTilesY=kernelStack.getHeight()/interpolateParameters.size;
		int outTilesX= (inTilesX-1)*interpolateParameters.step +1 + interpolateParameters.add_left + interpolateParameters.add_right;
		int outTilesY= (inTilesY-1)*interpolateParameters.step +1 + interpolateParameters.add_top + interpolateParameters.add_bottom;
		int nChn=kernelStack.getSize();
		float [][] outPixels=new float[nChn][outTilesX*interpolateParameters.size*outTilesY*interpolateParameters.size];
		float [] pixels;
		int i,j,chn;
		int xTile0=(interpolateParameters.add_left>0)?-1:0;
		int xTile1=inTilesX+((interpolateParameters.add_right>0)?0:-1);
		int yTile0=(interpolateParameters.add_top>0)?-1:0;
		int yTile1=inTilesY+((interpolateParameters.add_bottom>0)?0:-1);
		int tileY,tileX; //,subTileY,subTileX;

		int tileWidth, tileHeight; // for inner cells (interpolateParameters.step+1)*(interpolateParameters.step+1), for outer includes exte row/column fro extrapolation
		//  int maxTileWidth= Math.max(interpolateParameters.step,1+Math.max(interpolateParameters.add_right,interpolateParameters.add_left));
		//  int maxTileHeight=Math.max(interpolateParameters.step,1+Math.max(interpolateParameters.add_bottom,interpolateParameters.add_top));
		boolean lastColumn=false;  //last column - inverse convert and copy the last column of rectangleFHT to the result array
		boolean lastRow=false;     //last row - interpolate, inverse convert and copy the last row of rectangleFHT to the result array

		double [] pointsVert;
		double [] pointsHor;
		double [][] fhtLine;
		double extraScale=interpolateParameters.extrapolate/interpolateParameters.step;
		int [] outTopLeft=new int [2]; // top left kernel in the output array
		int [] inTopLeft=new int [2]; // top left kernel in the input array
		double [][] firstFHTColumn=null;
		double [][] secondFHTColumn=null;
		double [][][] cornerFHT=new double[2][2][interpolateParameters.size*interpolateParameters.size]; //[y][x][pixel] 
		double [] swapArray=null;

		for (chn=0;chn<nChn;chn++) {
			pixels=(float[]) kernelStack.getPixels(chn+1);
			for (tileY=yTile0;tileY<yTile1;tileY++) {
				if (updateStatus) IJ.showStatus("Interpolating kernels, channel "+kernelStack.getSliceLabel(chn+1)+", row "+(tileY-yTile0+1)+" of "+(yTile1-yTile0));
				lastRow=(tileY==(yTile1-1));
				if (tileY<0) {
					inTopLeft[1]=0;
					tileHeight=interpolateParameters.add_top;
					outTopLeft[1]=0;
					pointsVert=new double[tileHeight];
					for (i=0;i<tileHeight;i++)  pointsVert[i]=(i-tileHeight)*extraScale; // negative values
				} else if (tileY>=(inTilesY-1)){
					inTopLeft[1]=tileY-1;
					tileHeight=interpolateParameters.add_bottom+1; // always last row, if got here at all (interpolateParameters.add_bottom>0)
					outTopLeft[1]=interpolateParameters.add_top+interpolateParameters.step*tileY;
					pointsVert=new double[tileHeight];
					for (i=0;i<tileHeight;i++)  pointsVert[i]=1.0+i*extraScale;
				} else {
					inTopLeft[1]=tileY;
					tileHeight=interpolateParameters.step+ (lastRow?1:0); // last tile row includes bottom outpout kernel row
					outTopLeft[1]=interpolateParameters.add_top+interpolateParameters.step*tileY;
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
						tileWidth=interpolateParameters.add_left;
						outTopLeft[0]=0;
						pointsHor=new double[tileWidth];
						for (i=0;i<tileWidth;i++)  pointsHor[i]=(i-tileWidth)*extraScale; // negative values
					} else if (tileX>=(inTilesX-1)){
						inTopLeft[0]=tileX-1;
						tileWidth=interpolateParameters.add_right+1; // always last columnw, if got here at all (interpolateParameters.add_right>0)
						outTopLeft[0]=interpolateParameters.add_left+interpolateParameters.step*tileX;
						pointsHor=new double[tileWidth];
						for (i=0;i<tileWidth;i++)  pointsHor[i]=1.0+ i*extraScale;
						// else keep both firstFHTColumn and secondFHTColumn
						if (DEBUG_LEVEL>2)  System.out.println("last column: tileX="+tileX);
					} else {
						inTopLeft[0]=tileX;
						tileWidth=interpolateParameters.step+ (lastColumn?1:0); // last tile column includes rightmost outpout kernel column
						outTopLeft[0]=interpolateParameters.add_left+interpolateParameters.step*tileX;
						pointsHor=new double[tileWidth];
						for (i=1;i<tileWidth;i++)  pointsHor[i]=(1.0*i)/tileWidth;
						//  if (DEBUG_LEVEL>2)  System.out.println("else: tileX="+tileX);
						if (tileX!=0) {
							firstFHTColumn=secondFHTColumn;
							secondFHTColumn=null; // invalidate
							//  if (DEBUG_LEVEL>2)  System.out.println(" secondFHTColumn==null");
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
						firstFHTColumn=fht_instance.interpolateFHT (cornerFHT[0][0],    // first FHT array
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
						secondFHTColumn=fht_instance.interpolateFHT (cornerFHT[0][1],    // first FHT array
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

						fhtLine=fht_instance.interpolateFHT ( firstFHTColumn[i],    // first FHT array
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
		ImageStack outStack=new ImageStack(outTilesX*interpolateParameters.size,outTilesY*interpolateParameters.size);
		for (chn=0;chn<nChn;chn++) {
			outStack.addSlice(kernelStack.getSliceLabel(chn+1), outPixels[chn]);
		}
		return outStack;
	}
/* ======================================================================== */
/* Used in interpolateKernelStack() */  
	private void storeOneKernel(
			float []  pixels, // float [] array of combined square kernels - will be filled
			double [] kernel, // square kernel to store
			int       numHor, // number of kernels in a row
			int        xTile, // horizontal number of kernel to store
			int        yTile) { // vertical number of kernel to store
		int length=kernel.length;
		int size=(int) Math.sqrt(length);
		int i,j;
		int pixelsWidth=numHor*size;
		int base=(yTile*pixelsWidth+xTile)*size;
		for (i=0;i<size;i++) for (j=0;j<size;j++) pixels[base+i*pixelsWidth+j]= (float) kernel[i*size+j];
	}

/* ======================================================================== */
	private void extractOneKernel(
			float []  pixels, //  array of combined square kernels, each 
			double [] kernel, // will be filled, should have correct size before call
			int       numHor, // number of kernels in a row
			int        xTile, // horizontal number of kernel to extract
			int        yTile) { // vertical number of kernel to extract
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




/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
	 * @param inverseParameters TODO*/

	public ImageStack  reversePSFKernelStack(
			final ImageStack            PSFStack, // stack of 3 32-bit (float) images, made of square kernels
			final EyesisAberrations.InverseParameters inverseParameters, // size (side of square) of direct PSF kernel
			final int                 threadsMax, // size (side of square) of reverse PSF kernel
			final boolean           updateStatus){  // update status info
		if (PSFStack==null) return null;
		final int tilesX=PSFStack.getWidth()/inverseParameters.dSize;
		final int tilesY=PSFStack.getHeight()/inverseParameters.dSize;
		final int nChn=PSFStack.getSize();
		final double [] sigmas={inverseParameters.blurIndividual,inverseParameters.blurIndividual,inverseParameters.blurChecker};
		final float [][] outPixels=new float[nChn][tilesX*inverseParameters.rSize*tilesY*inverseParameters.rSize];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final int numberOfKernels=     tilesY*tilesX*nChn;
		final int numberOfKernelsInChn=tilesY*tilesX;
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					float [] pixels=null;
					double [] kernel= new double[inverseParameters.dSize*inverseParameters.dSize];
					double [] rKernel=new double[inverseParameters.rSize*inverseParameters.rSize];
					int  [][]selection;
					double [] ellipse_coeff;
					double [] variableSigmas;
					int chn,tileY,tileX;
					int chn0=-1;
					DoubleFHT fht_instance =new DoubleFHT(); // provide DoubleFHT instance to save on initializations (or null)
					for (int nTile = ai.getAndIncrement(); nTile < numberOfKernels; nTile = ai.getAndIncrement()) {
						chn=nTile/numberOfKernelsInChn;
						tileY =(nTile % numberOfKernelsInChn)/tilesX;
						tileX = nTile % tilesX;
						if (updateStatus) IJ.showStatus("Invertinging PSF, channel "+(chn+1)+" of "+nChn+", row "+(tileY+1)+" of "+tilesY);
						if (chn!=chn0) {
							pixels=(float[]) PSFStack.getPixels(chn+1);
							chn0=chn;
						}
						extractOneKernel( pixels, //  array of combined square kernels, each 
								kernel, // will be filled, should have correct size before call
								tilesX, // number of kernels in a row
								tileX, // horizontal number of kernel to extract
								tileY); // vertical number of kernel to extract
						/* Find direct kernel approximation ellipse, increase it, mirror center around 0,0 and use it as a mask for the reversed kernel */
						selection=    findClusterOnPSF(kernel,  inverseParameters.psfCutoffEnergy, "");
						ellipse_coeff=findEllipseOnPSF(kernel,  selection, ""); // coefficients for direct PSF, for rPSF [0] and [1] need to be opposite size

						rKernel=resizeForFFT(kernel,inverseParameters.rSize);
/* Apply variable blur to direct kernel using it's center X,Y */
						if (inverseParameters.filterDirect) {
							variableSigmas= createSigmasFromCenter(inverseParameters.rSize, // side of square
									inverseParameters.sigmaToRadiusDirect, // variable blurring - sigma will be proportional distance from the center
									sigmas[chn]*inverseParameters.sigmaScaleDirect, //blurring in the center sigma(r)=sqrt((sigma_to_radius*r)^2)+center_sigma^2)
									ellipse_coeff[0], // coordinates of the center (0:0 - size/2: size/2)
									ellipse_coeff[1]);
							rKernel=variableGaussBlurr(          rKernel, // input square pixel array, preferably having many exact zeros (they will be skipped)
									variableSigmas, // array of sigmas to be used for each pixel, matches pixels[]
									3.5, // drop calculatin if farther then nSigma
									0, // int WOICenterX, // window of interest in pixels[] array - do not generate data outside it
									0, // int WOICenterY, // 
									inverseParameters.rSize, //int WOIWidth, reduce later
									inverseParameters.rSize); //int WOIHeight)
						}
						
/* reverse PSF kernel */
						rKernel= cleanupAndReversePSF (rKernel,  // input pixels
								inverseParameters,
								//    						  false,  // fold high frequency into low, when false - use Hamming to cut off high frequencies
								fht_instance,
						""); // just for the plot names
/*  mask  the reversed kernel */
						rKernel= maskReversePSFKernel(rKernel, // reversed psf, square array
								ellipse_coeff, // ellipse coefficients from _direct_ kernel
								inverseParameters.psfEllipseScale,
								inverseParameters.rpsfMinMaskThreshold); // zero output element if elliptical Gauss mask is below this threshold

						normalizeKernel(rKernel); // in-place
/* Apply variable blur to inversed kernel, using (and reversing sign) the center X,Y from the direct kernel */
						if (inverseParameters.filter) {
							variableSigmas= createSigmasFromCenter(inverseParameters.rSize, // side of square
									inverseParameters.sigmaToRadius, // variable blurring - sigma will be proportional distance from the center
									sigmas[chn]*inverseParameters.sigmaScale, //blurring in the center sigma(r)=sqrt((sigma_to_radius*r)^2)+center_sigma^2)
									-ellipse_coeff[0], // coordinates of the center (0:0 - size/2: size/2)
									-ellipse_coeff[1]);
							rKernel=variableGaussBlurr(          rKernel, // input square pixel array, preferrably having many exact zeros (they will be skipped)
									variableSigmas, // array of sigmas to be used for each pixel, matches pixels[]
									3.5, // drop calculation if farther then nSigma
									0, // int WOICenterX, // window of interest in pixels[] array - do not generate data outside it
									0, // int WOICenterY, // 
									inverseParameters.rSize, //int WOIWidth, reduce later
									inverseParameters.rSize); //int WOIHeight)

						}
						//TODO: verify if it is OK that sum changed (was 10.5) after variableGaussBlurr(), for now - just re-calibrate
						normalizeKernel(rKernel); // in-place

						storeOneKernel( outPixels[chn], // float [] array of combined square kernels - will be filled
								rKernel, // square kernel to store
								tilesX, // number of kernels in a row
								tileX, // horizontal number of kernel to store
								tileY); // vertical number of kernel to store

					}
				}
			};
		}
		startAndJoin(threads);
		//	  System.out.println("Threads done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
/* prepare result stack to return */
		final ImageStack outStack=new ImageStack(tilesX*inverseParameters.rSize,tilesY*inverseParameters.rSize);
		for (int chn=0;chn<nChn;chn++) {
			outStack.addSlice(PSFStack.getSliceLabel(chn+1), outPixels[chn]);
		}
		return outStack;
	}



/* ======================================================================== */


/* calculates and updates global array PSF_KERNEL_MAP */
	public double [][][][] createPSFMap(
			final MatchSimulatedPattern commonMatchSimulatedPattern, // to be cloned in therads, using common data
			final ImagePlus          imp_sel, // linearized Bayer mosaic image form the camera, GR/BG
			final int [][][]        sampleList, // optional (or null) 2-d array: list of coordinate pairs (2d - to match existent  PSF_KERNEL_MAP structure)  
			final double  overexposedAllowed, // fraction of pixels OK to be overexposed
			final SimulationPattern.SimulParameters simulParameters,
			final int             mapFFTsize, // scanImageForPatterns:FFT size
			final MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			final int            fft_overlap,
			final int               fft_size,
			final EyesisAberrations.ColorComponents colorComponents,
			final int           PSF_subpixel, 
			final EyesisAberrations.OTFFilterParameters otfFilterParameters,
			final EyesisAberrations.PSFParameters psfParameters,
			final double       minDefinedArea,
			final int          PSFKernelSize, // size of square used in the new map (should be multiple of map step)
			final int          threadsMax,
			final boolean      updateStatus,          // UPDATE_STATUS
			final int          debug_level){// debug level used inside loops
		final long startTime = System.nanoTime();
		  Runtime runtime = Runtime.getRuntime();
		  runtime.gc();
		  if (DEBUG_LEVEL>1) System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");

		// Generate hi-res pattern bitmap (one cell)	
		SimulationPattern simulationPattern= new SimulationPattern();
		simulationPattern.debugLevel=DEBUG_LEVEL;
		final double [] bitmaskPattern= simulationPattern.patternGenerator(simulParameters);
		int nTileX,nTileY;
		int numPatternCells=0;
/* Filter results based on correlation with the actual pattern */
		boolean [][]   PSFBooleanMap; // map of 2*fft_size x 2*fft_size squares with 2*fft_overlap step, true means that that tile can be used for PSF
		final boolean useOld=(matchSimulatedPattern==null) || (matchSimulatedPattern.getUVIndex()==null);
		if (sampleList==null){
			if (useOld){
				IJ.showMessage("Error","MASK_ARRAY does not exist, process distortions first\nUsing old method");
				double [][][][] patternMap=scanImageForPatterns(imp_sel,
						mapFFTsize,
						patternDetectParameters,
						bitmaskPattern,
						simulParameters,
						//			  simulFill,
						1); // debug level to be used while scanning cells

				for (nTileY=0;nTileY<patternMap.length;nTileY++) for (nTileX=0;nTileX<patternMap[0].length;nTileX++) if (patternMap[nTileY][nTileX]!=null) numPatternCells++;
				if (DEBUG_LEVEL>1) {
					System.out.println("Finished mapping, created array["+patternMap.length+"]["+patternMap[0].length+"][][], "+
							numPatternCells+" cells (of "+(patternMap.length*patternMap[0].length)+") with pattern detected");
				}
				if (numPatternCells==0) {
					IJ.showMessage("Error","There are no cells with patternn\nProcess canceled");
					return null;
				}
				PSFBooleanMap= remapSquares(patternMap, // [][]map of either null or 2 wave vectors
						mapFFTsize/2, // step of initial map
						mapFFTsize, // size of square used in scanning of initial map (should be multiple of map step)
						fft_overlap, // step of the new map (should be multiple of map step)
						fft_size); // size of square used in the new map (should be multiple of map step)


			} else {
				PSFBooleanMap= mapFromPatternMask ( // count number of defined cells
						matchSimulatedPattern.getUVIndex(), // int array, >=0 - uv exist, <0 - rmpty 
						//					matchSimulatedPattern.getWOI().width, // image (mask) width
						imp_sel.getWidth(), // image (mask) width
						fft_size*2,
						fft_overlap*2,
						fft_size,   // backward compatibility margin==tileSize/2
						GAUSS_WIDTH,
						//	psfParameters.minDefinedArea);
						minDefinedArea);
			}
		} else {
			PSFBooleanMap= new boolean[sampleList.length][sampleList[0].length];
			for (int i=0;i<sampleList.length;i++) for (int j=0;j<sampleList[0].length;j++) PSFBooleanMap[i][j]=(sampleList[i][j][0]>=0); // all with positive X
		}
		numPatternCells=0;
		for (nTileY=0;nTileY<PSFBooleanMap.length;nTileY++) for (nTileX=0;nTileX<PSFBooleanMap[0].length;nTileX++) if (PSFBooleanMap[nTileY][nTileX]) numPatternCells++;
		if (DEBUG_LEVEL>1) {
			System.out.println("Remapped for PSF measurment- converted to an array["+PSFBooleanMap.length+"]["+PSFBooleanMap[0].length+"], "+
					numPatternCells+" cells (of "+(PSFBooleanMap.length*PSFBooleanMap[0].length)+") with pattern detected");
		}
		PSF_KERNEL_MAP=new double[PSFBooleanMap.length][PSFBooleanMap[0].length][][]; //PSF_KERNEL_MAP - global (or final)
		int saved_DEBUG_LEVEL=DEBUG_LEVEL;
		DEBUG_LEVEL=debug_level;
		simulationPattern.debugLevel=DEBUG_LEVEL;
		int ncell=0;
/* Create array of coordinates of cells to process, fill result array with zeros (to be actually written by threads */       
		final int [][] tilesToProcessXY=new int [numPatternCells][4];

		for (nTileY=0;nTileY<PSFBooleanMap.length;nTileY++) for (nTileX=0;nTileX<PSFBooleanMap[0].length;nTileX++){
			if (PSFBooleanMap[nTileY][nTileX]) {
				tilesToProcessXY[ncell  ][0]=nTileX;
				tilesToProcessXY[ncell  ][1]=nTileY;
				tilesToProcessXY[ncell  ][2]=(sampleList==null)?(fft_overlap*2*nTileX):sampleList[nTileY][nTileX][0];
				tilesToProcessXY[ncell++][3]=(sampleList==null)?(fft_overlap*2*nTileY):sampleList[nTileY][nTileX][1];
				PSF_KERNEL_MAP[nTileY][nTileX]=new double[colorComponents.colorsToCorrect.length][];
			} else PSF_KERNEL_MAP[nTileY][nTileX]=null;
		}
		final Thread[] threads = newThreadArray(threadsMax);
		 if (DEBUG_LEVEL>1) System.out.println("Starting "+threads.length+" threads: "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		final AtomicInteger ai = new AtomicInteger(0);
		final int patternCells=numPatternCells;
		//	  final double []   overexposedMap, // map of overexposed pixels in the image (may be null)
		final double [] overexposed=(overexposedAllowed>0)?JP4_INSTANCE.overexposedMap (imp_sel):null;
		final int mapWidth=imp_sel.getWidth();
		for (int ithread = 0; ithread < threads.length; ithread++) {
			// Concurrently run in as many threads as CPUs
			threads[ithread] = new Thread() {
				public void run() {

					// Each thread processes a few items in the total list
					// Each loop iteration within the run method has a unique 'i' number to work with
					// and to use as index in the results array:
					//	double [] sum_kern_el=new double[6]; // just testing					
					int x0,y0,nTX,nTY,nChn;
					double [][] kernels;
					double [] windowFFTSize=    initWindowFunction(fft_size); //=initHamming( fft_size) calculate once
					double [] windowFullFFTSize=initWindowFunction(fft_size*PSF_subpixel); //=initHamming( fft_size*subpixel);
					DoubleFHT fht_instance =new DoubleFHT(); // provide DoubleFHT instance to save on initializations (or null)
					double over;
// individual per-thread - will be needed when converted to doubleFHT					
//				    MatchSimulatedPattern matchSimulatedPattern=new MatchSimulatedPattern(FFT_SIZE);
				    MatchSimulatedPattern matchSimulatedPattern=commonMatchSimulatedPattern.clone();
				    matchSimulatedPattern.debugLevel=DEBUG_LEVEL;
					SimulationPattern simulationPattern= new SimulationPattern(bitmaskPattern);
					simulationPattern.debugLevel=DEBUG_LEVEL;
					for (int nTile = ai.getAndIncrement(); nTile < patternCells; nTile = ai.getAndIncrement()) {
						nTX=tilesToProcessXY[nTile][0];
						nTY=tilesToProcessXY[nTile][1];
//						y0=fft_overlap*2*nTY;
//						x0=fft_overlap*2*nTX;
						y0=tilesToProcessXY[nTile][3];
						x0=tilesToProcessXY[nTile][2];
						if (updateStatus) IJ.showStatus("Processing tile["+nTY+"]["+nTX+"] ("+(nTile+1)+" of "+patternCells+")");
						if (MASTER_DEBUG_LEVEL>1) System.out.println("Processing tile["+nTY+"]["+nTX+"] ("+(nTile+1)+" of "+patternCells+") : "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
						if (overexposed!=null){
							over=JP4_INSTANCE.fracOverExposed(overexposed,   // map of overexposed pixels 0.0 - 0K, >0 (==1.0) - overexposed
									mapWidth,    // width of the map
									x0,          // X of the top left corner of the selection
									y0,          // Y of the top left corner of the selection
									2*fft_size,  // selection width
									2*fft_size); // selection height
							//						  if ((DEBUG_LEVEL>0) && (over>0.0)) System.out.println("Overexposed fraction of "+over+" at x0="+x0+" y0="+y0+" width"+(2*fft_size));
						} else over=-1.0;
						if ( over > overexposedAllowed) {
							PSF_KERNEL_MAP[nTY][nTX]=null;
							if (DEBUG_LEVEL>0) System.out.println("Overexposed fraction of "+over+" at x0="+x0+" y0="+y0+" width"+(2*fft_size));
						} else {
							kernels=getPSFKernels(imp_sel,
									useOld?null:SIM_ARRAY, //simulation image, scaled PSF_subpixel/2
									2*fft_size,       // size in pixels (twice fft_size)
									x0,               // top left corner X (pixels)
									y0,               // top left corner Y (pixels)
									simulationPattern,
									matchSimulatedPattern,
									patternDetectParameters,
									windowFFTSize,    //=initHamming( fft_size) calculate once
									windowFullFFTSize,//=initHamming( fft_size*subpixel);
									PSF_subpixel,     // use finer grid than actual pixels 
									simulParameters,
									colorComponents,  // color channels to process, equalizeGreens
									otfFilterParameters,
									5,                // int referenceComp, // number of color component to reference lateral chromatic aberration to (now 4 - checkerboard greens)
									psfParameters,
									fht_instance,      // provide DoubleFHT instance to save on initializations (or null)
									debug_level // ((x0<512)&& (y0<512))?3:debug_level DEBUG during "focusing"
							);
							if (kernels!=null) {
								if (kernelLength(kernels)>(PSFKernelSize*PSFKernelSize)) kernels=resizeForFFT(kernels,PSFKernelSize); // shrink before normalizing
								normalizeKernel(kernels); // in-place
								if (kernelLength(kernels)<(PSFKernelSize*PSFKernelSize)) kernels=resizeForFFT(kernels,PSFKernelSize); // expand after normalizing
								for (nChn=0;nChn<kernels.length;nChn++) if (kernels[nChn]!=null){
									PSF_KERNEL_MAP[nTY][nTX][nChn]=kernels[nChn]; // not .clone()?
								}
//(new showDoubleFloatArrays()).showArrays(kernels, "***kernels-"+nTX+"-"+nTY);
							} else {
								if (MASTER_DEBUG_LEVEL>1) System.out.println("Empty kernel for tile["+nTY+"]["+nTX+"]");
							}
							//save results into common array
							//PSF_KERNEL_MAP[nTY][nTX]
						}
					}
				}
			};
		}
		startAndJoin(threads);
		 if (DEBUG_LEVEL>1) System.out.println("Threads done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		DEBUG_LEVEL=saved_DEBUG_LEVEL;
		return PSF_KERNEL_MAP;
	}
	public void createPSFMapOld(
			final ImagePlus          imp_sel, // linearized Bayer mosaic image form the camera, GR/BG
			final double  overexposedAllowed, // fraction of pixels OK to be overexposed
			final SimulationPattern.SimulParameters simulParameters,
			final int             mapFFTsize, // scanImageForPatterns:FFT size
			final MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			final int            fft_overlap,
			final int               fft_size,
			final EyesisAberrations.ColorComponents colorComponents,
			final int           PSF_subpixel, 
			final EyesisAberrations.OTFFilterParameters otfFilterParameters,
			final EyesisAberrations.PSFParameters psfParameters,
			final int          PSFKernelSize, // size of square used in the new map (should be multiple of map step)
			final int          threadsMax,
			final int               debug_level){// debug level used inside loops
		final long startTime = System.nanoTime();

		// Generate hi-res pattern bitmap (one cell)	
		SimulationPattern simulationPattern= new SimulationPattern();
		simulationPattern.debugLevel=DEBUG_LEVEL;
		final double [] bitmaskPattern= simulationPattern.patternGenerator(simulParameters);
		double [][][][] patternMap=scanImageForPatterns(imp_sel,
				mapFFTsize,
				patternDetectParameters,
				bitmaskPattern,
				simulParameters,
				//			  simulFill,
				1); // debug level to be used while scanning cells

		int nTileX,nTileY;
		int numPatternCells=0;
		for (nTileY=0;nTileY<patternMap.length;nTileY++) for (nTileX=0;nTileX<patternMap[0].length;nTileX++) if (patternMap[nTileY][nTileX]!=null) numPatternCells++;
		if (DEBUG_LEVEL>1) {
			System.out.println("Finished mapping, created array["+patternMap.length+"]["+patternMap[0].length+"][][], "+
					numPatternCells+" cells (of "+(patternMap.length*patternMap[0].length)+") with pattern detected");
		}
		if (numPatternCells==0) {
			IJ.showMessage("Error","There are no cells with patternn\nProcess canceled");
			return;
		}
/* Filter results based on correlation with the actual pattern */
		boolean [][]   PSFBooleanMap; // map of 2*fft_size x 2*fft_size squares with 2*fft_overlap step, true means that that tile can be used for PSF
		PSFBooleanMap= remapSquares(patternMap, // [][]map of either null or 2 wave vectors
				mapFFTsize/2, // step of initial map
				mapFFTsize, // size of square used in scanning of initial map (should be multiple of map step)
				fft_overlap, // step of the new map (should be multiple of map step)
				fft_size); // size of square used in the new map (should be multiple of map step)
		numPatternCells=0;
		for (nTileY=0;nTileY<PSFBooleanMap.length;nTileY++) for (nTileX=0;nTileX<PSFBooleanMap[0].length;nTileX++) if (PSFBooleanMap[nTileY][nTileX]) numPatternCells++;
		if (DEBUG_LEVEL>1) {
			System.out.println("Remapped for PSF measurment- converted to an array["+PSFBooleanMap.length+"]["+PSFBooleanMap[0].length+"], "+
					numPatternCells+" cells (of "+(PSFBooleanMap.length*PSFBooleanMap[0].length)+") with pattern detected");
		}
		PSF_KERNEL_MAP=new double[PSFBooleanMap.length][PSFBooleanMap[0].length][][]; //PSF_KERNEL_MAP - global (or final)
		DEBUG_LEVEL=debug_level;
		simulationPattern.debugLevel=DEBUG_LEVEL;
		int ncell=0;
/* Create array of coordinates of cells to process, fill result array with zeros (to be actually written by threads */       
		final int [][] tilesToProcessXY=new int [numPatternCells][2];
		for (nTileY=0;nTileY<PSFBooleanMap.length;nTileY++) for (nTileX=0;nTileX<PSFBooleanMap[0].length;nTileX++){
			if (PSFBooleanMap[nTileY][nTileX]) {
				tilesToProcessXY[ncell][0]=nTileX;
				tilesToProcessXY[ncell++][1]=nTileY;
				PSF_KERNEL_MAP[nTileY][nTileX]=new double[colorComponents.colorsToCorrect.length][];
			} else PSF_KERNEL_MAP[nTileY][nTileX]=null;
		}
		final Thread[] threads = newThreadArray(threadsMax);
		System.out.println("Starting "+threads.length+" threads: "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		final AtomicInteger ai = new AtomicInteger(0);
		final int patternCells=numPatternCells;
		//	  final double []   overexposedMap, // map of overexposed pixels in the image (may be null)
		final double [] overexposed=JP4_INSTANCE.overexposedMap (imp_sel);
		final int mapWidth=imp_sel.getWidth();
		for (int ithread = 0; ithread < threads.length; ithread++) {
			// Concurrently run in as many threads as CPUs
			threads[ithread] = new Thread() {
				public void run() {

					// Each thread processes a few items in the total list
					// Each loop iteration within the run method has a unique 'i' number to work with
					// and to use as index in the results array:
					//	double [] sum_kern_el=new double[6]; // just testing					
					int x0,y0,nTX,nTY,nChn;
					double [][] kernels;
					double [] windowFFTSize=    initWindowFunction(fft_size); //=initHamming( fft_size) calculate once
					double [] windowFullFFTSize=initWindowFunction(fft_size*PSF_subpixel); //=initHamming( fft_size*subpixel);
					DoubleFHT fht_instance =new DoubleFHT(); // provide DoubleFHT instance to save on initializations (or null)
					double over;
// individual per-thread - will be needed when converted to doubleFHT					
				    MatchSimulatedPattern matchSimulatedPattern=new MatchSimulatedPattern(FFT_SIZE);
				    matchSimulatedPattern.debugLevel=DEBUG_LEVEL;
					SimulationPattern simulationPattern= new SimulationPattern(bitmaskPattern);
					simulationPattern.debugLevel=DEBUG_LEVEL;
					for (int nTile = ai.getAndIncrement(); nTile < patternCells; nTile = ai.getAndIncrement()) {
						nTX=tilesToProcessXY[nTile][0];
						nTY=tilesToProcessXY[nTile][1];
						y0=fft_overlap*2*nTY;
						x0=fft_overlap*2*nTX;
						if (MASTER_DEBUG_LEVEL>0) IJ.showStatus("Processing tile["+nTY+"]["+nTX+"] ("+(nTile+1)+" of "+patternCells+")");
						if (MASTER_DEBUG_LEVEL>0) System.out.println("Processing tile["+nTY+"]["+nTX+"] ("+(nTile+1)+" of "+patternCells+") : "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
						if (overexposed!=null){
							over=JP4_INSTANCE.fracOverExposed(overexposed,   // map of overexposed pixels 0.0 - 0K, >0 (==1.0) - overexposed
									mapWidth,    // width of the map
									x0,          // X of the top left corner of the selection
									y0,          // Y of the top left corner of the selection
									2*fft_size,  // selection width
									2*fft_size); // selection height
							//						  if ((DEBUG_LEVEL>0) && (over>0.0)) System.out.println("Overexposed fraction of "+over+" at x0="+x0+" y0="+y0+" width"+(2*fft_size));
						} else over=-1.0;
						if ( over > overexposedAllowed) {
							PSF_KERNEL_MAP[nTY][nTX]=null;
							if (DEBUG_LEVEL>0) System.out.println("Overexposed fraction of "+over+" at x0="+x0+" y0="+y0+" width"+(2*fft_size));
						} else {
							kernels=getPSFKernels (imp_sel,
									null,
									2*fft_size,       // size in pixels (twice fft_size)
									x0,               // top left corner X (pixels)
									y0,               // top left corner Y (pixels)
									simulationPattern,
									matchSimulatedPattern,
									patternDetectParameters,
									windowFFTSize,    //=initHamming( fft_size) calculate once
									windowFullFFTSize,//=initHamming( fft_size*subpixel);
									PSF_subpixel,     // use finer grid than actual pixels 
									simulParameters,
									colorComponents,  // color channels to process, equalizeGreens
									otfFilterParameters,
									5,                // int referenceComp, // number of color component to reference lateral chromatic aberration to (now 4 - checkerboard greens)
									psfParameters,
									fht_instance,      // provide DoubleFHT instance to save on initializations (or null)
									debug_level // ((x0<512)&& (y0<512))?3:debug_level
							);

							if (kernelLength(kernels)>(PSFKernelSize*PSFKernelSize)) kernels=resizeForFFT(kernels,PSFKernelSize); // shrink before normalizing
							normalizeKernel(kernels); // in-place
							if (kernelLength(kernels)<(PSFKernelSize*PSFKernelSize)) kernels=resizeForFFT(kernels,PSFKernelSize); // expand after normalizing
							for (nChn=0;nChn<kernels.length;nChn++) if (kernels[nChn]!=null){
								PSF_KERNEL_MAP[nTY][nTX][nChn]=kernels[nChn]; // not .clone()?
							}
							//save results into common array
							//PSF_KERNEL_MAP[nTY][nTX]
						}

					}
				}
			};
		}
		startAndJoin(threads);
		System.out.println("Threads done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		return;
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
	public double [][][][] scanImageForPatterns(ImagePlus imp,
			int                                        size, // FFT size
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			double [] bPattern,
			SimulationPattern.SimulParameters  simulParameters,
//			double simulFill,
			int debug_level) {  // debug level to use while iterating through steps
		if (imp==null){
			IJ.showMessage("Error","No image specified\nProcess canceled");
			return null;
		}
	    MatchSimulatedPattern matchSimulatedPattern=new MatchSimulatedPattern(size); // here size does not matter?
	    matchSimulatedPattern.debugLevel=DEBUG_LEVEL;
	    SimulationPattern     simulationPattern=    new SimulationPattern (bPattern); //  reuse bPattern
		simulationPattern.debugLevel=DEBUG_LEVEL;
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
/* title is the same - window, not file */
			System.out.println("Mapping image with squares of "+(size*2)+"x"+(size*2)+"pixels, with 50% overlap, covering total area of "+
					((mapWidth+1)*size)+"x"+((mapHeight+1)*size));
		}
		double [][][][] patternMap=new double [mapHeight][mapWidth][][];
		Rectangle mapCell=new Rectangle (0,0,size*2,size*2);
		int nTileX,nTileY,yc,xc;
		int wasDebug=DEBUG_LEVEL;
		double [][] pixels;
		double [] hamming=initWindowFunction(size);
/* Prepare for filtering - common part outside of the iteration */
		double [][]  sim_pix;
		SimulationPattern.SimulParameters  thisSimulParameters=simulParameters.clone();
//		int this_simul_subdiv=2;
		thisSimulParameters.subdiv=2;
		double [] model_corr;
		double [][] convMatrix= {{1.0,-1.0},{1.0,1.0}}; // from greens2 to pixel WV
		double [][] invConvMatrix= matrix2x2_scale(matrix2x2_invert(convMatrix),2.0);
		double [][] WVgreens;
		double contrast;
		DEBUG_LEVEL=debug_level; // modify DEBUG_LEVEL to mute it while scanning many cells
//		DoubleFHT fht_instance =new DoubleFHT(); // provide DoubleFHT instance to save on initializations (or null)
		matchSimulatedPattern.debugLevel=DEBUG_LEVEL;
		simulationPattern.debugLevel=DEBUG_LEVEL;

		for (nTileY=0;nTileY<mapHeight;nTileY++) {
			if (DEBUG_LEVEL>0) IJ.showStatus("Mapping row "+(nTileY+1)+" (of "+mapHeight+")");
			for (nTileX=0;nTileX<mapHeight;nTileX++) patternMap[nTileY][nTileX]=null;
			yc=(nTileY+1)*size;
//			if ((yc<selection.y) || (yc>=selection.y+selection.height)) continue;
			if ((yc<(selection.y+size)) || (yc>(selection.y+selection.height-size))) continue;
			
			mapCell.y=nTileY*size;
			for (nTileX=0;nTileX<mapWidth;nTileX++) {
				xc=(nTileX+1)*size;
//				if ((xc<selection.x) || (xc>=selection.x+selection.width)) continue;
				if ((xc<(selection.x+size)) || (xc>(selection.x+selection.width-size))) continue;
				mapCell.x=nTileX*size;
				pixels=splitBayer(imp, mapCell,COMPONENTS.equalizeGreens);
				pixels[4]= normalizeAndWindow (pixels[4], hamming);
				patternMap[nTileY][nTileX]=matchSimulatedPattern.findPattern(pixels[4],
						size,
						patternDetectParameters,
						true,
						title); // title - will not be used
/* Now verify by correlating with the actual pattern */
				if ((patternMap[nTileY][nTileX]!=null) && (patternDetectParameters.minCorrContrast>0.0)) {
					simulationPattern.simulatePatternFullPattern(
							patternMap[nTileY][nTileX][0][0],
							patternMap[nTileY][nTileX][0][1],
							patternMap[nTileY][nTileX][0][2],
							patternMap[nTileY][nTileX][1][0],
							patternMap[nTileY][nTileX][1][1],
							patternMap[nTileY][nTileX][1][2],
							null, // no mesh distortion here
							thisSimulParameters.subdiv,// SIMUL.subdiv, - do not need high quality here
							size,
							true); // center for greens
					sim_pix= simulationPattern.extractSimulPatterns (
							thisSimulParameters,
							1,
							size,    // number of Bayer cells in width of the square selection (half number of pixels)
							0.0,    // selection center, X (in pixels)
							0.0);   // selection center, y (in pixels)
					sim_pix[4]= normalizeAndWindow (sim_pix[4], hamming);
//TODO: test and replace 	matchSimulatedPattern.correlateWithModel with 	fht_instance.correlate !!!			
					model_corr=matchSimulatedPattern.correlateWithModel (pixels[4],  // measured pixel array
							sim_pix[4],  // simulated (model) pixel array)
							0.0,  // double sigma,   // Sigma for high pass filtering
							imp.getTitle());
//					model_corr=fht_instance.correlate(pixels[4],sim_pix[4],0); // destroys operands
					WVgreens=matrix2x2_mul(patternMap[nTileY][nTileX],invConvMatrix);
					contrast= matchSimulatedPattern.correlationContrast (model_corr,    // square pixel array
							WVgreens,    // wave vectors (same units as the pixels array)
							patternDetectParameters.corrRingWidth,   // ring (around r=0.5 dist to opposite corr) width
							0.0,    //  x0,              // center coordinates
							0.0,    //y0,
							title);   // title base for optional plots names
					//      System.out.println("Pattern correlation contrast= "+IJ.d2s(contrast,3)+ ", threshold is "+PATTERN_DETECT.minCorrContrast);
					if (!(contrast >= patternDetectParameters.minCorrContrast)) patternMap[nTileY][nTileX]=null; // still getting NaN sometimes
				}
			}
		}
		DEBUG_LEVEL=wasDebug; // restore original debug level
		//COMPONENT_COLOR_NAMES
		//  public double[][] splitBayer (ImagePlus imp, Rectangle r) {
		//imp_src.getWidth(), imp_src.getHeight(),imp_src.getTitle()
		return patternMap;
		//
	}

/* ======================================================================== */
/* build a map of available overlapping squares of different size than initially mapped. New size and step should be multiple of initial step */
	private boolean [][] remapSquares (double [][][][] map, // [][]map of either null or 2 wave vectors
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
			System.out.println("mapStep= "+mapStep+" mapSquare="+mapSquare+" newStep="+newStep+" newSquare="+newSquare);
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

	private boolean [][] mapFromPatternMask (
			int [] uvIndex, // grid detected >=0, no grid -1
			int width, // image (mask) width
			int tileSize,
			int tileStep,
			int margin,   // backward compatibility margin==tileSize/2
			double gaussWidth,
			double threshold){
		double[] windowFunction= initWindowFunction(tileSize, gaussWidth);
		int height =uvIndex.length/width;
		int tileHeight=(height-2*margin)/tileStep+1;
		int tileWidth= (width- 2*margin)/tileStep+1;
		boolean [][] result = new boolean [tileHeight][tileWidth];
		int index;
		int len=tileSize*tileSize;
		double absThresh=0.0, sum;
		for (index=0;index<len;index++) absThresh+=windowFunction[index];
		absThresh*=threshold;
		if (DEBUG_LEVEL>1) System.out.println(" threshold="+threshold+" absThresh="+absThresh);
		
		int y,x,y0,x0;
		for (int tileY=0;tileY<tileHeight;tileY++) for (int tileX=0;tileX<tileWidth;tileX++) {
			y0=-tileSize/2+margin+tileStep*tileY;
			x0=-tileSize/2+margin+tileStep*tileX;
			sum=0;
			for (index=0;index<len;index++) {
				y=index/tileSize+y0;
				x=index%tileSize+x0;
				if ((y>=0) && (x>=0) && (y<height) && (x<width) && (uvIndex[y*width+x]>=0)) sum+= windowFunction[index];
//				if ((DEBUG_LEVEL>0) && (tileY==22) && (tileX==32)) System.out.println(" x="+x+" y="+y);
//				if ((DEBUG_LEVEL>0) && (tileY==22) && (tileX==32) && (y>=0) && (x>=0) && (y<height) && (x<width))System.out.println(" uvIndex["+(y*width+x)+"]="+uvIndex[y*width+x]);
			}
			result[tileY][tileX]=(sum>absThresh);
			if (DEBUG_LEVEL>1) System.out.println("mapFromPatternMask(..,"+threshold +" tileY="+tileY+" tileX="+tileX+" x0="+x0+" y0="+y0+" sum="+sum);
			
		}
		return result;
	}
	
	
/* ======================================================================== */
	public double [][] matrix2x2_invert(double [][] m ){
		double det=m[0][0]*m[1][1]-m[0][1]*m[1][0];
		double [][] rslt= {{ m[1][1]/det,  -m[0][1]/det},
				{-m[1][0]/det,   m[0][0]/det}};
		return rslt;
	}
	public double [][] matrix2x2_mul(double [][] a, double [][] b ){
		double [][] rslt={{a[0][0]*b[0][0]+a[0][1]*b[1][0], a[0][0]*b[0][1]+a[0][1]*b[1][1]},
				{a[1][0]*b[0][0]+a[1][1]*b[1][0], a[1][0]*b[0][1]+a[1][1]*b[1][1]}};
		return rslt;
	}
	public double [] matrix2x2_mul(double [][] a, double [] b ){
		double [] rslt={a[0][0]*b[0]+a[0][1]*b[1],
				a[1][0]*b[0]+a[1][1]*b[1]};
		return rslt;
	}
	public double [][] matrix2x2_scale(double [][] a, double  b ){
		double [][] rslt={{a[0][0]*b, a[0][1]*b},
				{a[1][0]*b, a[1][1]*b}};
		return rslt;
	}

	public double [][] matrix2x2_add(double [][] a, double [][] b ){
		double [][] rslt={{a[0][0]+b[0][0], a[0][1]+b[0][1]},
		         		  {a[1][0]+b[1][0], a[1][1]+b[1][1]}};
		return rslt;
	}

	public double [] matrix2x2_add(double [] a, double [] b ){
		double [] rslt={a[0]+b[0], a[1]+b[1]};
		return rslt;
	}
	
	public double [][] matrix2x2_transp(double [][] m ){
		double [][] rslt= {{ m[0][0],  m[1][0]},
		            	   { m[0][1],  m[1][1]}};
		return rslt;
	}
	
/* ======================================================================== */
/* Use ROI */
/* Supply rectangle */
/* Accept slices, default to first slice */
/*	
	private double[][] splitBayer (ImagePlus imp, Rectangle r, boolean equalize_greens) {
		return splitBayer ((float[])imp.getProcessor().getPixels(), imp.getWidth(),  r,  equalize_greens);
	}
*/
	private double[][] splitBayer (ImagePlus imp, Rectangle r, boolean equalize_greens) {
		return splitBayer (imp, 1, r, equalize_greens);
	}
	private double[][] splitBayer (ImagePlus imp,  int sliceNumber, Rectangle r, boolean equalize_greens) {
		if (imp.getStackSize()>1){
		  return splitBayer ((float[])imp.getStack().getProcessor(sliceNumber).getPixels(), imp.getWidth(),  r,  equalize_greens);
		} else {
		  return splitBayer ((float[])imp.getProcessor().getPixels(), imp.getWidth(),  r,  equalize_greens);
		}  
    }    

	private double[][] splitBayer (float [] pixels, int full_width, Rectangle r, boolean equalize_greens) {
		int full_height=pixels.length/full_width; // full image height
		if (r==null) r=new Rectangle(0,0,full_width,full_height);
		if (DEBUG_LEVEL>10) IJ.showMessage("splitBayer","r.width="+r.width+
				"\nr.height="+r.height+
				"\nr.x="+r.x+
				"\nr.y="+r.y+
				"\nlength="+pixels.length);
		if ((DEBUG_LEVEL>2) && ((r.x<0) || (r.y<0) || ((r.x+r.width)>=full_width) || ((r.y+r.height)>=full_height))) System.out.println("r.width="+r.width+
				" r.height="+r.height+
				" r.x="+r.x+
				" r.y="+r.y);
		int x,y,base,base_b,bv,i,j;
		int half_height=r.height>>1;
		int half_width=r.width>>1;
        // make them all 0 if not a single pixel falls into the image        
        int numColors=(half_height==half_width)?5:4;
        int pixX,pixY;
        double [][] bayer_pixels=new double[numColors][half_height * half_width];
        if ((r.x>=full_width) || (r.y>=full_height) || ((r.x+r.width)<0)  || ((r.y+r.height)<0)) {
        	for (i=0;i<bayer_pixels.length;i++) for (j=0;j<bayer_pixels[i].length;j++) bayer_pixels[i][j]=0.0;
            return bayer_pixels;
        }
        //      base=r.width*((y<<1)+bv);
        for (y=0; y<half_height; y++) for (bv=0;bv<2;bv++){
        	pixY=(y*2)+bv+r.y;
        	base_b=half_width*y;
        	//					if ((pixY>=0)

        	if (pixY<0) {
        		pixY=bv;
        	} else if (pixY>=full_height){
        		pixY=full_height-2+bv;
        	}
        	base=full_width*pixY+((r.x>0)?r.x:0);
        	//						base=full_width*((y*2)+bv+r.y)+r.x;
        	pixX=r.x;
        	if (bv==0) for (x=0; x<half_width; x++) {
        		if ((pixX<0) || (pixX>=(full_width-2))) {
        			bayer_pixels[0][base_b]= pixels[base];
        			bayer_pixels[1][base_b]= pixels[base+1];
        		} else {
        			bayer_pixels[0][base_b]= pixels[base++];
        			bayer_pixels[1][base_b]= pixels[base++];
        		}
        		base_b++;
        		pixX+=2;
        	} else  for (x=0; x<half_width; x++) {
        		if ((pixX<0) || (pixX>=(full_width-2))) {
        			bayer_pixels[2][base_b]= pixels[base];
        			bayer_pixels[3][base_b]= pixels[base+1];
        		} else {
        			bayer_pixels[2][base_b]= pixels[base++];
        			bayer_pixels[3][base_b]= pixels[base++];
        		}
        		base_b++;
        		pixX+=2;
        	}
        }
        if (equalize_greens) {
        	double g0=0.0,g3=0.0,g02=0.0,g32=0.0,a0,a3,b0,b3;
        	int n=bayer_pixels[0].length;
        	for (i=0;i<bayer_pixels[0].length;i++) {
        		g0 +=bayer_pixels[0][i];
        		g02+=bayer_pixels[0][i]*bayer_pixels[0][i];
        		g3 +=bayer_pixels[3][i];
        		g32+=bayer_pixels[3][i]*bayer_pixels[3][i];
        	}
        	g0/=n; // mean value
        	g3/=n; // meran value
        	g02=g02/n-g0*g0;
        	g32=g32/n-g3*g3;
        	b0=Math.sqrt(Math.sqrt(g32/g02));
        	b3 = 1.0/b0;
        	a0= (g0+g3)/2 -b0*g0;
        	a3= (g0+g3)/2 -b3*g3;
        	if (DEBUG_LEVEL>2) {
        		System.out.println("g0= "+g0+ ", g3= "+g3);
        		System.out.println("g02="+g02+", g32="+g32);
        		System.out.println("a0="+a0+", b0="+b0);
        		System.out.println("a3="+a3+", b3="+b3);
        	}
        	for (i=0;i<bayer_pixels[0].length;i++) {
        		bayer_pixels[0][i]=a0+bayer_pixels[0][i]*b0;
        		bayer_pixels[3][i]=a3+bayer_pixels[3][i]*b3;
        	}

        }

        if (numColors>4) bayer_pixels[4]=combineDiagonalGreens (bayer_pixels[0], bayer_pixels[3],  half_width, half_height);
        return bayer_pixels;
	}

	public double[][] splitBayerZero (ImagePlus imp, Rectangle r, boolean equalize_greens) {
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
			if (equalize_greens) {
				double g0=0.0,g3=0.0,g02=0.0,g32=0.0,a0,a3,b0,b3;
				int n=bayer_pixels[0].length;
				for (i=0;i<bayer_pixels[0].length;i++) {
					g0 +=bayer_pixels[0][i];
					g02+=bayer_pixels[0][i]*bayer_pixels[0][i];
					g3 +=bayer_pixels[3][i];
					g32+=bayer_pixels[3][i]*bayer_pixels[3][i];
				}
				g0/=n; // mean value
				g3/=n; // meran value
				g02=g02/n-g0*g0;
				g32=g32/n-g3*g3;
				b0=Math.sqrt(Math.sqrt(g32/g02));
				b3 = 1.0/b0;
				a0= (g0+g3)/2 -b0*g0;
				a3= (g0+g3)/2 -b3*g3;
				if (DEBUG_LEVEL>2) {
					System.out.println("g0= "+g0+ ", g3= "+g3);
					System.out.println("g02="+g02+", g32="+g32);
					System.out.println("a0="+a0+", b0="+b0);
					System.out.println("a3="+a3+", b3="+b3);
				}
				for (i=0;i<bayer_pixels[0].length;i++) {
					bayer_pixels[0][i]=a0+bayer_pixels[0][i]*b0;
					bayer_pixels[3][i]=a3+bayer_pixels[3][i]*b3;
				}

			}

			if (numColors>4) bayer_pixels[4]=combineDiagonalGreens (bayer_pixels[0], bayer_pixels[3],  half_width, half_height);
			return bayer_pixels;
	}

	private  double [] combineDiagonalGreens (double [] green0, double []green3, int half_width, int half_height) {
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

/* ======================================================================== */
	private  void normalizeKernel(double [][] kernel) {
		int i;
		for (i=0;i<kernel.length;i++) if (kernel[i]!=null) normalizeKernel(kernel[i]);
	}

	private  void normalizeKernel(double [] kernel) {
		//	    if (kernel==null) return null;
		int i;
		double s=0;
		for (i=0;i<kernel.length;i++) s+= kernel[i];
		s=1.0/s;
		for (i=0;i<kernel.length;i++) kernel[i]*=s;
	}

/* ======================================================================== */
/* extends/shrinks image to make it square for FFT */
	private double[][] resizeForFFT (double[][]kernels, int size) {
		if (kernels==null) return null;
		double [][]result=new double [kernels.length][];
		for (int i=0;i<kernels.length;i++) {
			if (kernels[i]!=null) result[i]=resizeForFFT(kernels[i],size);
			else result[i]=null;
		}
		return result;
	}

	private double[] resizeForFFT (double[]kernel, int size) {
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
/* ======================================================================== */

	private  int kernelLength(double[][]kernels) {
		if (kernels==null) return 0;
		for (int i=0; i<kernels.length;i++) if (kernels[i]!=null) return kernels[i].length;
		return 0;
	}

/* ======================================================================== */
	private  ImageStack mergeKernelsToStack(double [][][][] kernels) {
		return mergeKernelsToStack(kernels,null);

	}
	private  ImageStack mergeKernelsToStack(double [][][][] kernels,String [] names) { // use oldStack.getSliceLabels() to get names[]
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
		int [] channelsMask = new int [kernel.length];
		for (i=0;i<kernel.length;i++) channelsMask[i]=0;
		for (i=0;i<tilesY ;i++)  for (j=0;j<tilesX;j++) if (kernels[i][j]!=null) {
			for (k=0;(k<kernel.length)&& (k<channelsMask.length);k++) if (kernels[i][j][k]!=null) {
				channelsMask[k]=1;
				if (kernels[i][j][k].length>length) length=kernels[i][j][k].length;
			}
		}

		nChn=0;
		for (i=0;i<channelsMask.length;i++) if (channelsMask[i]!=0) nChn++;
		int [] channels = new int [nChn];
		nChn=0;
		for (i=0;i<channelsMask.length;i++) if (channelsMask[i]!=0) channels[nChn++]=i;

		//	    for (i=0;i<kernel.length;i++) if (kernel[i]!=null)  if (nChn<channels.length) channels[nChn++]=i;
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
					if ((kernels[i][j]==null || (kernels[i][j][k]==null))) fpixels[index]=0.0f;
					else {
						fpixels[index]= (float) kernels[i][j][k][y*size+x];
					}
				}
			}
			if (names==null) stack.addSlice("channel"+k, fpixels);
			else             stack.addSlice(names[chn], fpixels);
		}
		return stack;
	}

/* ======================================================================== */
	public double [][] getPSFKernels ( ImagePlus imp,
			float [][] simArray, //simulation image, scaled PSF_subpixel/2 (or null), [0] - main pixels, [1] - shifted diagonally by 0.5 pixels (for checker greens)
			int size,   // size in pixels (twice FFT_SIZE)
			int x0,      // top left corner X (pixels)
			int y0,      // top left corner Y (pixels)
			SimulationPattern     simulationPattern,
		    MatchSimulatedPattern matchSimulatedPattern,
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			double [] Hamming, //=initHamming( fft_size) calculate once
			double [] fullHamming, //=initHamming( fft_size*subpixel);
			int subpixel, // use finer grid than actual pixels 
			SimulationPattern.SimulParameters  simulParameters,
			EyesisAberrations.ColorComponents colorComponents,
			EyesisAberrations.OTFFilterParameters otfFilterParameters,
			int referenceComp, // number of color component to reference lateral chromatic aberration to (now 4 - checkerboard greens)
			EyesisAberrations.PSFParameters psfParameters,
			DoubleFHT fht_instance, // provide DoubleFHT instance to save on initializations (or null)
			int debug
	){
		if (DEBUG_LEVEL>1){
			System.out.println("getPSFKernels(), simArray is "+((simArray==null)?"":"not ")+"null"); 
		}

		if (imp==null) return null; // Maybe convert to double pixel array once to make it faster?
		if (fht_instance==null) fht_instance=new DoubleFHT(); // move upstream to reduce number of initializations
		double [][] kernels=         new double[6][];  // was global
		String title=imp.getTitle()+"X"+x0+"Y"+y0;
		Rectangle PSFCell=new Rectangle (x0,y0,size,size);
		int fft_size=size/2;
		double [][] input_bayer=splitBayer (imp,PSFCell,colorComponents.equalizeGreens);
		//int greensToProcess=4;
		int i,j,l;
		double [][] simul_pixels;
		double [][]wVectors=new double[2][2];
		int imgWidth=imp.getWidth();
		
		double [][] dbgSimPix=null;
		
		if ((simArray==null) || (psfParameters.approximateGrid)){ // just for testing
			/* Calculate pattern parameters, including distortion */
			if (matchSimulatedPattern.PATTERN_GRID==null) {
				double[][] distortedPattern= matchSimulatedPattern.findPatternDistorted(input_bayer, // pixel array to process (no windowing!)
						patternDetectParameters,
						true, //(greensToProcess==4), // boolean greens, // this is a pattern for combined greens (diagonal), adjust results accordingly
						title); // title prefix to use for debug  images

				if (distortedPattern==null) return null;
				if (DEBUG_LEVEL>3){
					System.out.println(
							" W0x="+     IJ.d2s(distortedPattern[0][0],4)+
							" W0y="+     IJ.d2s(distortedPattern[0][1],4)+
							" W0_phase="+IJ.d2s(distortedPattern[0][2],2)+
							" W1x="+     IJ.d2s(distortedPattern[1][0],4)+
							" W1y="+     IJ.d2s(distortedPattern[1][1],4)+
							" W1_phase="+IJ.d2s(distortedPattern[1][2],2));

				}
				simulationPattern.simulatePatternFullPattern(
						distortedPattern[0][0],
						distortedPattern[0][1],
						distortedPattern[0][2],
						distortedPattern[1][0],
						distortedPattern[1][1],
						distortedPattern[1][2],
						distortedPattern[2], //
						simulParameters.subdiv,
						fft_size,
						simulParameters.center_for_g2);
				wVectors[0][0]=2.0*distortedPattern[0][0]/subpixel;
				wVectors[0][1]=2.0*distortedPattern[0][1]/subpixel;
				wVectors[1][0]=2.0*distortedPattern[1][0]/subpixel;
				wVectors[1][1]=2.0*distortedPattern[1][1]/subpixel;
			} else { // approximate pattern grid and simulate it
				double[][] distPatPars=null;
				try {
				distPatPars= matchSimulatedPattern.findPatternFromGrid(
						x0, // top-left pixel of the square WOI
						y0,
						size, // size of square (pixels)
						Hamming, // only half-window!
						false,  // use linear approximation (instead of quadratic)
						1.0E-10,  // thershold ratio of matrix determinant to norm for linear approximation (det too low - fail)
						1.0E-20);  // thershold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
				} catch (Exception e){
					System.out.println("Failed findPatternFromGrid("+x0+","+y0+",...)");
				}
				if (distPatPars==null) return null;
				if ((distPatPars[0].length<6) || (distPatPars[1].length<6)){
					System.out.println("Failure: findPatternFromGrid("+x0+","+y0+",...) returned only linear coefficients");
					return null;
				}
				int [] iUV={(int) Math.floor(distPatPars[0][2]),(int) Math.floor(distPatPars[1][2])};
				boolean negative=((iUV[0]^iUV[1])&1)!=0;
				double [] simCorr={
						distPatPars[0][3]/4,
						distPatPars[0][4]/4,
						distPatPars[0][5]/4,
						distPatPars[1][3]/4,
						distPatPars[1][4]/4,
						distPatPars[1][5]/4,
						0.0,0.0,0.0,0.0};
				double [] phases={
						1.0*Math.PI*(distPatPars[0][2]-iUV[0]+(negative?(-0.5):0.5)), // measured from the center of white
						1.0*Math.PI*(distPatPars[1][2]-iUV[1]+0.5)};
//				double [][]wVectors={{distPatPars[0][0],distPatPars[0][1]},{distPatPars[1][0],distPatPars[1][1]}};
				wVectors[0][0]=distPatPars[0][0];
				wVectors[0][1]=distPatPars[0][1];
				wVectors[1][0]=distPatPars[1][0];
				wVectors[1][1]=distPatPars[1][1];
				simulationPattern.simulatePatternFullPattern(
						wVectors[0][0],
						wVectors[0][1],
						phases[0],
						wVectors[1][0],
						wVectors[1][1],
						phases[1],
						simCorr, //
						simulParameters.subdiv,
						fft_size,
						simulParameters.center_for_g2);
			}
			simul_pixels= simulationPattern.extractSimulPatterns (
					simulParameters,
					subpixel, // subdivide pixels
					fft_size*subpixel, // number of Bayer cells in width of the square selection (half number of pixels)
					0.0,    // selection center, X (in pixels)
					0.0);   // selection center, y (in pixels)
			if (subpixel>1) {
				if (colorComponents.colorsToCorrect[5])  simul_pixels=combineCheckerGreens (simul_pixels,   // pixel arrays after oversampleFFTInput() or extractSimulPatterns())
						subpixel); // same as used in oversampleFFTInput() - oversampling ratio
			}
			for (i=0;i<simul_pixels.length; i++) {
				if (!colorComponents.colorsToCorrect[i]) simul_pixels[i]=null; // removed unused
			}
			simul_pixels= normalizeAndWindow (simul_pixels, fullHamming);
		} else {
			Rectangle PSFCellSim=new Rectangle (x0*subpixel/2,y0*subpixel/2,size*subpixel/2,size*subpixel/2);
			
			simul_pixels=new double[6][];
// simulationPattern.debugLevel=DEBUG_LEVEL;
			for (i=0;i<simul_pixels.length; i++) {
				if (colorComponents.colorsToCorrect[i]) simul_pixels[i]=simulationPattern.extractBayerSim (
						SIM_ARRAY, // [0] - regular pixels, [1] - shifted by 1/2 diagonally, for checker greens
						imgWidth*subpixel/2,
						PSFCellSim,
						subpixel, // 4
						i);
				else simul_pixels[i]=null;
			}
//System.out.println("PSFCell.y="+PSFCell.y+" PSFCell.height="+PSFCell.height+" imgWidth="+imgWidth+" PSFCell.x="+PSFCell.x+" PSFCell.width="+PSFCell.width+" matchSimulatedPattern.UV_INDEX.length="+matchSimulatedPattern.UV_INDEX.length);
			int index=matchSimulatedPattern.getUVIndex((PSFCell.y+PSFCell.height/2)*imgWidth+(PSFCell.x+PSFCell.width/2));
//			int index=matchSimulatedPattern.getUVIndex((PSFCell.y+PSFCell.height/2)*matchSimulatedPattern.getWOI().width+(PSFCell.x+PSFCell.width/2));
			
			if (index<0) {
				System.out.println ("Error, No UV pattern @ x="+(PSFCell.x+PSFCell.width/2)+", y="+(PSFCell.y+PSFCell.height/2));
				return null;
			}
			//			int [] iUV={index % matchSimulatedPattern.getDArrayHeight(), index / matchSimulatedPattern.getDArrayHeight()}; // TODO: make sure it is correct?
			int [] iUV={index % matchSimulatedPattern.getDArrayWidth(), index / matchSimulatedPattern.getDArrayWidth()}; // TODO: make sure it is correct?
			if (matchSimulatedPattern.getDArray(iUV[1],iUV[0])==null) {
				if (DEBUG_LEVEL>0){
					System.out.println ( "Tried to extract wave vectors from non-existent node "+iUV[0]+"/"+iUV[1]);
					System.out.println ( "index="+index+"  matchSimulatedPattern.getDArrayHeight()"+ matchSimulatedPattern.getDArrayHeight());
					System.out.println("PSFCell.y="+PSFCell.y+" PSFCell.height="+PSFCell.height+" imgWidth="+imgWidth+" PSFCell.x="+PSFCell.x+" PSFCell.width="+PSFCell.width+
							" matchSimulatedPattern.UV_INDEX.length="+matchSimulatedPattern.UV_INDEX.length);
				}
				return null;
			}
			if (matchSimulatedPattern.getDArray(iUV[1],iUV[0],1)==null) {
				if (DEBUG_LEVEL>0) System.out.println ( "Tried to extract non-existent wave vectors from "+iUV[0]+"/"+iUV[1]);
				return null;
			}
			//TODO:  Need to define wave vectors here - how?
			wVectors[0]=matchSimulatedPattern.getDArray(iUV[1],iUV[0],1); //null pointer
			wVectors[1]=matchSimulatedPattern.getDArray(iUV[1],iUV[0],2);
			// should it be averaged WV?			
			if (DEBUG_LEVEL>2) System.out.println ( " x0="+x0+" y0="+y0);
			if (DEBUG_LEVEL>2) SDFA_INSTANCE.showArrays(input_bayer, true, title+"-in");
			if (DEBUG_LEVEL>2) SDFA_INSTANCE.showArrays(simul_pixels, true, title+"-S");
			//if (DEBUG_LEVEL>2) System.out.println (SIM_ARRAY[0][-1]); // cause error
			if (MASTER_DEBUG_LEVEL>1){
				dbgSimPix=new double[simul_pixels.length][];
				for (int ii=0;ii<dbgSimPix.length;ii++)
					if (simul_pixels[ii]!=null) dbgSimPix[ii]=simul_pixels[ii].clone();
					else dbgSimPix[ii]=null;

			}
			simul_pixels= normalizeAndWindow (simul_pixels, fullHamming);
		}
		
		input_bayer= normalizeAndWindow (input_bayer, Hamming);
		if (subpixel>1) {
			input_bayer= oversampleFFTInput (input_bayer,subpixel);
			if (colorComponents.colorsToCorrect[5])  input_bayer=combineCheckerGreens (input_bayer,   // pixel arrays after oversampleFFTInput() or extractSimulPatterns())
					subpixel); // same as used in oversampleFFTInput() - oversampling ratio
		}
		for (i=0;i<4;i++) if (!colorComponents.colorsToCorrect[i]) input_bayer[i]=null; // leave composite greens even if disabled
		if (DEBUG_LEVEL>3) {
			SDFA_INSTANCE.showArrays(input_bayer, fft_size*subpixel, fft_size*subpixel, title);
		}
		if (DEBUG_LEVEL>2) System.out.println ( " input_bayer.length="+input_bayer.length+" simul_pixels.length="+simul_pixels.length+" fft_size*subpixel="+fft_size*subpixel);
		for (i=0;(i<input_bayer.length) && (i<simul_pixels.length);i++) if ((colorComponents.colorsToCorrect[i]) && (input_bayer[i]!=null)){
			if (DEBUG_LEVEL>2) System.out.println ( "input_bayer["+i+"].length="+input_bayer[i].length+" simul_pixels["+i+"].length="+simul_pixels[i].length);
		}
		if (DEBUG_LEVEL>2) SDFA_INSTANCE.showArrays(input_bayer, true, title+"-input");
		if (DEBUG_LEVEL>2) SDFA_INSTANCE.showArrays(simul_pixels, true, title+"-SIM");

if (DEBUG_LEVEL>2)DEBUG_LEVEL=0; //************************************************************
		double [][] inverted=new double[colorComponents.colorsToCorrect.length][];
		double wvAverage=Math.sqrt(0.5*(wVectors[0][0]*wVectors[0][0]+wVectors[0][1]*wVectors[0][1]+
				wVectors[1][0]*wVectors[1][0]+wVectors[1][1]*wVectors[1][1]));

		for (i=0;(i<input_bayer.length) && (i<simul_pixels.length);i++) if ((colorComponents.colorsToCorrect[i]) && (input_bayer[i]!=null)){
			if (DEBUG_LEVEL>2) System.out.println ( "Color "+COMPONENT_COLOR_NAMES[i]+" is re-calculated into bayer pixels ");
			if (DEBUG_LEVEL>2) System.out.println ( "input_bayer["+i+"].length="+input_bayer[i].length+" simul_pixels["+i+"].length="+simul_pixels[i].length);
			inverted[i]=limitedInverseOfFHT(input_bayer[i],
					simul_pixels[i],
					fft_size*subpixel,
					(i==5),     //    boolean checker // checkerboard pattern in the source file (use when filtering)
					true, //      forwardOTF,
					subpixel,
					otfFilterParameters,
					fht_instance,
					psfParameters.mask1_sigma*size*wvAverage,      // normalize to wave vectors!
					psfParameters.mask1_threshold,
					psfParameters.gaps_sigma*size*wvAverage,
					psfParameters.mask_denoise,
					debug,
					title+"-"+i);
		}
		if (DEBUG_LEVEL>2) SDFA_INSTANCE.showArrays(inverted, fft_size*subpixel, fft_size*subpixel, title+"_Combined-PSF");
/* correct composite greens */
/* Here we divide wave vectors by subpixel as the pixels are already added */
		double [][] wVrotMatrix= {{0.5,0.5},{-0.5,0.5}};
		double [][]wVectors4= new double [2][2];
		for (i=0;i<2;i++) for (j=0;j<2;j++) {
			wVectors4[i][j]=0.0;
			for (l=0;l<2;l++) wVectors4[i][j]+=wVectors[i][l]*wVrotMatrix[l][j];
		}
		double [][] PSF_shifts=         new double [input_bayer.length][];    // X/Y shift of the PSF array, in Bayer component pixel coordinates (same as PSF arrays)
		double [][] PSF_centroids=      new double [input_bayer.length][];    // X/Y coordinates of the centroids of PSF in Bayer component pioxel coordinates (same as PSF arrays) (after they were optionally shifted)
		double [][] lateralChromatic=   new double [input_bayer.length][]; // X/Y coordinates of the centroids of Bayer component PSF in sensor pixel coordinates
		double [][] kernelsForFFT=      new double [input_bayer.length][];
		double [][] psf_inverted=       new double [input_bayer.length][];
		double [][] psf_inverted_masked=new double [input_bayer.length][];
		double [] lateralChromaticAbs=new double [input_bayer.length];
		double [] zeroVector={0.0,0.0};
		for (i=input_bayer.length-1;i>=0;i--) {
			if (colorComponents.colorsToCorrect[i]) {
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
		//int [][]  clusterMask;
/* Start with referenceComp */
		i= referenceComp;
		if (DEBUG_LEVEL>2) {
			System.out.println("1-PSF_shifts.length= "+PSF_shifts.length+" i="+i+" input_bayer.length="+input_bayer.length);
			System.out.println("Before: color Component "+i+" PSF_shifts["+i+"][0]="+IJ.d2s(PSF_shifts[i][0],3)+
					" PSF_shifts["+i+"][1]="+IJ.d2s(PSF_shifts[i][1],3));
		}  


		kernels[i]=combinePSF (inverted[i], // Square array of pixels with multiple repeated PSF (alternating sign)
				true, // master, force ignoreChromatic
				PSF_shifts[i],  // centerXY[] - will be modified inside combinePSF() if PSF_PARS.ignoreChromatic is true
				PSF_centroids[i], // will return array of XY coordinates of the result centroid
				(i==4)?wVectors4:wVectors, // two wave vectors, lengths in cycles/pixel (pixels match pixel array)
						psfParameters,
						fht_instance,
						title+"_"+i,    // reduce the PSF cell size to this part of the area connecting first negative clones
						(DEBUG_LEVEL>4));
		if (DEBUG_LEVEL>2)     System.out.println("After-1: color Component "+i+"    PSF_shifts["+i+"][0]="+IJ.d2s(PSF_shifts   [i][0],3)+"    PSF_shifts["+i+"][1]="+IJ.d2s(   PSF_shifts[i][1],3));
		if (DEBUG_LEVEL>2)     System.out.println("After-1: color Component "+i+" PSF_centroids["+i+"][0]="+IJ.d2s(PSF_centroids[i][0],3)+" PSF_centroids["+i+"][1]="+IJ.d2s(PSF_centroids[i][1],3));

		if (!psfParameters.ignoreChromatic) { /* Recalculate center to pixels from greens (diagonal)) and supply it to other colors (lateral chromatic aberration correction) */
			for (j=0;j<input_bayer.length;j++) if ((colorComponents.colorsToCorrect[j]) && (j!=referenceComp)) {
				PSF_shifts[j]=shiftSensorToBayer (shiftBayerToSensor(PSF_shifts[referenceComp],referenceComp,subpixel),j,subpixel);
				if (DEBUG_LEVEL>2)       System.out.println("After-2 (recalc): color Component "+j+" PSF_shifts["+j+"][0]="+IJ.d2s(PSF_shifts[j][0],3)+" PSF_shifts["+j+"][1]="+IJ.d2s(PSF_shifts[j][1],3));
			}
		}

		lateralChromatic[i]=shiftBayerToSensor ( PSF_shifts[i][0]+PSF_centroids[i][0],
				PSF_shifts[i][1]+PSF_centroids[i][1],
				i,
				subpixel);
		lateralChromaticAbs[i]=Math.sqrt(lateralChromatic[i][0]*lateralChromatic[i][0]+lateralChromatic[i][1]*lateralChromatic[i][1]);
/* Now process all the other components */
		for (i=0; i<input_bayer.length;i++) if ((i!=referenceComp) && (colorComponents.colorsToCorrect[i])) {
			if (DEBUG_LEVEL>2) {
				System.out.println("2-PSF_shifts.length= "+PSF_shifts.length+" i="+i+" input_bayer.length="+input_bayer.length);

				System.out.println("Before: color Component "+i+" PSF_shifts["+i+"][0]="+IJ.d2s(PSF_shifts[i][0],3)+
						" PSF_shifts["+i+"][1]="+IJ.d2s(PSF_shifts[i][1],3));
			}  
			kernels[i]=combinePSF (inverted[i], // Square array of pixels with multiple repeated PSF (alternating sign)
					false, // !master, use ignoreChromatic
					PSF_shifts[i],  // centerXY[] - will be modified inside combinePSF() if psfParameters.ignoreChromatic is true
					PSF_centroids[i], // will return array of XY coordinates of the result centroid
					(i==4)?wVectors4:wVectors, // two wave vectors, lengths in cycles/pixel (pixels match pixel array)
							psfParameters,
							fht_instance,
							title+"_"+i,    // reduce the PSF cell size to this part of the area connecting first negative clones
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
		if (DEBUG_LEVEL>1) { //1
			for (i=0;i<PSF_shifts.length;i++) if (colorComponents.colorsToCorrect[i]){
				if (DEBUG_LEVEL>2) { //2
					System.out.println("Color Component "+i+" subpixel="+subpixel+
							" psfParameters.ignoreChromatic="+psfParameters.ignoreChromatic+
							" psfParameters.symm180="+psfParameters.symm180);
					System.out.println(                     " PSF_shifts["+i+"][0]="+IJ.d2s(PSF_shifts[i][0],3)+
							" PSF_shifts["+i+"][1]="+IJ.d2s(PSF_shifts[i][1],3)+
							" PSF_centroids["+i+"][0]="+IJ.d2s(PSF_centroids[i][0],3)+
							" PSF_centroids["+i+"][1]="+IJ.d2s(PSF_centroids[i][1],3));
					System.out.println("  lateralChromatic["+i+"][0]="+IJ.d2s(lateralChromatic[i][0],3)+
							"  lateralChromatic["+i+"][1]="+IJ.d2s(lateralChromatic[i][1],3));
				}
			}
			if (colorComponents.colorsToCorrect[referenceComp]) for (i=0;i<colorComponents.colorsToCorrect.length;i++) if ((colorComponents.colorsToCorrect[i])&& (i!=referenceComp)){
				System.out.println(COMPONENT_COLOR_NAMES[i]+" lateral chromatic (from green) "+IJ.d2s(lateralChromaticAbs[i],3)+"pix(sensor):  ["+i+"][0]="+IJ.d2s(lateralChromatic[i][0]-lateralChromatic[referenceComp][0],3)+
						"  ["+i+"][1]="+IJ.d2s(lateralChromatic[i][1]-lateralChromatic[referenceComp][1],3));
			}
			System.out.println("Lateral shift green from simulation "+IJ.d2s(lateralChromaticAbs[referenceComp],3)+"pix(sensor):  ["+referenceComp+"][0]="+IJ.d2s(lateralChromatic[referenceComp][0],3)+
					"  ["+referenceComp+"][1]="+IJ.d2s(lateralChromatic[referenceComp][1],3));
		}
		return kernels;
	}




/* ======================================================================== */
/* Combine both greens as a checkerboard pattern (after oversampleFFTInput()) */
	private  double [][] combineCheckerGreens (double[][] input_pixels,   // pixel arrays after oversampleFFTInput() or extractSimulPatterns())
			int ratio) { // same as used in oversampleFFTInput() - oversampling ratio
		int width=(int) Math.sqrt(input_pixels[0].length);
		return combineCheckerGreens (input_pixels,   // pixel arrays after oversampleFFTInput() or extractSimulPatterns())
				width,   // width of the image
				ratio);
	}

	private  double [][] combineCheckerGreens (double[][] input_pixels,   // pixel arrays after oversampleFFTInput() or extractSimulPatterns())
			int width,   // width of the image
			int ratio) { // same as used in oversampleFFTInput() - oversampling ratio
		if (DEBUG_LEVEL>5) System.out.println ("combineCheckerGreens(), ratio="+ratio+" input_pixels.length="+input_pixels.length);

		if ((ratio<2) ||
				(input_pixels==null) ||
				((input_pixels.length>5) && (input_pixels[5]!=null)) ||
				(input_pixels.length<4) ||
				(input_pixels[0]==null) ||
				(input_pixels[3]==null)) return input_pixels;
		int height=input_pixels[0].length/width;
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

/* ======================================================================== */
/* inserts zeros between pixels */ 
	private  double [][] oversampleFFTInput (double[][] input_pixels,
			int ratio) {
		double [][] pixels=new double[input_pixels.length][];
		int i;
		for (i=0;i<pixels.length;i++) pixels[i]= oversampleFFTInput (input_pixels[i], ratio);
		return pixels;
	}


	private  double [] oversampleFFTInput (double[] input_pixels, int ratio) {
		if (input_pixels==null) return null;
		int width=(int) Math.sqrt(input_pixels.length);
		return oversampleFFTInput (input_pixels,
				width,   // width of the image
				ratio);
	}

	private  double [] oversampleFFTInput (double[] input_pixels,
			int width,   // width of the image
			int ratio) {
		if (input_pixels==null) return null;
		if (DEBUG_LEVEL>2) System.out.println ("oversampleFFTInput(), width="+width+" ratio="+ratio+" input_pixels.length="+input_pixels.length);
		double [] pixels=new double[input_pixels.length*ratio*ratio];
		int i,j,x,y;
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

/* ======================================================================== */

	private double[] limitedInverseOfFHT(double [] measuredPixels,  // measured pixel array
			double [] modelPixels,  // simulated (model) pixel array)
			int size,  // FFT size
			boolean checker,  // checkerboard pattern in the source file (use when filtering)
			boolean forward_OTF,  // divide measured by simulated when true, simulated by measured - when false
			int oversample,  // measured array is sampled at 1/oversample frequency than model (will add more parameters later)
			EyesisAberrations.OTFFilterParameters filterOTFParameters,  //  fraction of the maximal value to be used to limit zeros
			DoubleFHT fht_instance,  // add rejection of zero frequency (~2-3pix)
			double mask1_sigma,
			double mask1_threshold,
			double gaps_sigma,
			double mask_denoise,
			int debug,
			String title){ // title base for optional plots names
		return limitedInverseOfFHT(measuredPixels,
				modelPixels,
				size,
				checker,
				forward_OTF,
				oversample,  // measured array is sampled at 1/oversample frequency than model (will add more parameters later)
				filterOTFParameters.deconvInvert,
				filterOTFParameters.zerofreqSize,  // add rejection of zero frequency (~2-3pix)
				filterOTFParameters.smoothPS,       // 0 - none, otherwise Gauss width
				filterOTFParameters.thresholdHigh,  // reject completely if energy is above this part of maximal
				filterOTFParameters.thresholdLow,  // leave intact if energy is below this part of maximal
				-1.0, // if 0 use normalize amplitude, if 0..1 - make binary: 1.0 if > threshold, 0.0 - otherwise -1 - disable mask
				0.0, // low-pass result with low pass filter (should be later defined automatically)
				fht_instance,
				mask1_sigma,
				mask1_threshold,
				gaps_sigma,
				mask_denoise,
				debug,
				title);
	}
// TODO: It now selects a single PSF, so combinePSF() and binPSF() can be simplified and eliminated
	private double[] limitedInverseOfFHT(double [] measuredPixels,  // measured pixel array
			double [] modelPixels,  // simulated (model) pixel array)
			int size,  // FFT size
			boolean checker,  // checkerboard pattern in the source file (use when filtering)
			boolean forward_OTF,  // divide measured by simulated when true, simulated by measured - when false
			int oversample,  // measured array is sampled at 1/oversample frequency than model (will add more parameters later)
			double deconvInvert,  //  fraction of the maximal value to be used to limit zeros
			double zerofreq_size,  // add rejection of zero frequency (~2-3pix)
			double smoothPS,       // 0 - none, otherwise Gauss width = FFT size/2/smoothPS
			double threshold_high,  // reject completely if energy is above this part of maximal
			double threshold_low,  // leave intact if energy is below this part of maximal
			double threshold, // if 0 use normalize amplitude, if 0..1 - make binary: 1.0 if > threshold, 0.0 - otherwise -1 - disable mask
			double radius, // low-pass result with low pass filter (should be later defined automatically)
			DoubleFHT fht_instance, // provide DoubleFHT instance to save on initializations (or null)
			double mask1_sigma,
			double mask1_threshold,
			double gaps_sigma,
			double mask_denoise,
			int debug,
			String title){
        
		double [] denominatorPixels= forward_OTF? modelPixels.clone():    measuredPixels.clone();
		double [] nominatorPixels=   forward_OTF? measuredPixels.clone(): modelPixels.clone();
		if (fht_instance==null) fht_instance=new DoubleFHT(); // move upstream to reduce number of initializations
		int i;
		fht_instance.swapQuadrants(denominatorPixels);
		fht_instance.transform(denominatorPixels);
		double [] mask= null;
		double [] mask1=null;
		DoubleGaussianBlur gb=new DoubleGaussianBlur();
		if ((oversample>1) && (threshold_low<1.0)) {
			double [] ps=fht_instance.calculateAmplitude2(denominatorPixels);
/* create mask */
			mask= maskAliases (denominatorPixels,   // complex spectrum, [size/2+1][size]
					checker, // checkerboard pattern in the source file (use when filtering)
					oversample,   // measured array is sampled at 1/oversample frequency than model (will add more parameters later)
					zerofreq_size,   // add rejection of zero frequency (~2-3pix)
					smoothPS,
					deconvInvert,
					threshold_high,   // reject completely if energy is above this part of maximal
					threshold_low,  // leave intact if energy is below this part of maximal
					fht_instance);
/* debug show the mask */
			if ((debug>1) ||((DEBUG_LEVEL>2) && (title!=""))) { /* Increase debug level later */ // was 3
				SDFA_INSTANCE.showArrays(mask, title+"-MASK");
			}
			for (int ii=0;ii<ps.length;ii++) ps[ii]=Math.log(ps[ii]); // can be twice faster
			if ((debug>1) ||((DEBUG_LEVEL>2) && (title!=""))) { /* Increase debug level later */ // was 3
				SDFA_INSTANCE.showArrays(ps, "LOG-"+title);
			}
			double [] ps_smooth=ps.clone();
			gb.blurDouble(ps_smooth, size, size, mask1_sigma, mask1_sigma, 0.01);
			if ((debug>1) ||((DEBUG_LEVEL>2) && (title!=""))) { /* Increase debug level later */ // was 3
				SDFA_INSTANCE.showArrays(ps_smooth, "SM-"+title);
			}
			double threshold1=Math.log(2.0*mask1_threshold);
			mask1=new double [ps.length];
			for (int ii=0;ii<ps.length;ii++) mask1[ii]= ps[ii]-ps_smooth[ii]-threshold1;
			if ((debug>1) ||((DEBUG_LEVEL>2) && (title!=""))) { /* Increase debug level later */ // was 3
				SDFA_INSTANCE.showArrays(mask1, "M1-"+title);
			}
			fht_instance.swapQuadrants(mask1); // zero in the corner
			for (int ii=0;ii<mask1.length;ii++){
				if (mask1[ii]<0) {
					//				mask[ii]=0.0;
					mask1[ii]=0.0;
				}
				mask1[ii]*=mask[ii];
			}
			if ((debug>1) ||((DEBUG_LEVEL>2) && (title!=""))) { /* Increase debug level later */ // was 3
				SDFA_INSTANCE.showArrays(mask1, "M1A-"+title);
			}
		}
/* Mask already includes zeros on ps, so we can just use divisions of FHT*/		  
		//Swapping quadrants of the nominator, so the center will be 0,0
		fht_instance.swapQuadrants(nominatorPixels);
		//get to frequency domain
		fht_instance.transform(nominatorPixels);
		if ((debug>1) ||((DEBUG_LEVEL>2) && (title!=""))) { /* Increase debug evel later */ // was 3
			SDFA_INSTANCE.showArrays(nominatorPixels, title+"-NOM-FHT");
			SDFA_INSTANCE.showArrays(denominatorPixels, title+"-DENOM-FHT");
		}			    
		double [] pixels=fht_instance.divide(nominatorPixels,denominatorPixels);
		if ((debug>1) ||((DEBUG_LEVEL>2) && (title!=""))) { /* Increase debug evel later */ // was 3
			SDFA_INSTANCE.showArrays(pixels, title+"-DECONV");
		}			    
		for (i=0;i<pixels.length;i++) {
			if (mask[i]==0.0) pixels[i]=0.0; // preventing NaN*0.0
			else pixels[i]*=mask[i];
		}
		if ((debug>1) ||((DEBUG_LEVEL>2) && (title!=""))) { /* Increase debug level later */ // was 3
			SDFA_INSTANCE.showArrays(pixels, title+"-MASKED");
			double [][] aphase=fht_instance.fht2AmpHase(pixels,true);
			SDFA_INSTANCE.showArrays(aphase, true,"AP="+title+"-MASKED");
			
		}
		if (gaps_sigma>0.0){
			double [][] fft_reIm_centered=fht_instance.fht2ReIm(pixels, true); //0 in the center, full square
			fht_instance.swapQuadrants(mask1); // zero in the center
			for (int ii=0;ii<2;ii++) for (int jj=0;jj<mask1.length;jj++) fft_reIm_centered[ii][jj]*=mask1[jj];
			gb.blurDouble(mask1, size, size, gaps_sigma, gaps_sigma, 0.01);
			gb.blurDouble(fft_reIm_centered[0], size, size, gaps_sigma, gaps_sigma, 0.01);
			gb.blurDouble(fft_reIm_centered[1], size, size, gaps_sigma, gaps_sigma, 0.01);
			for (int ii=0;ii<2;ii++) for (int jj=0;jj<mask1.length;jj++)
				if (mask1[jj]>mask_denoise) fft_reIm_centered[ii][jj]/=mask1[jj];
				else if (mask1[jj]>=0.0) fft_reIm_centered[ii][jj]/=mask_denoise;
				else  fft_reIm_centered[ii][jj]=0.0;
			if ((debug>1) ||((DEBUG_LEVEL>2) && (title!=""))) { /* Increase debug level later */ // was 3
				SDFA_INSTANCE.showArrays(fft_reIm_centered, true,"ReIm-"+title);
			}
			fht_instance.swapQuadrants(fft_reIm_centered[0]); // zero in the corner
			fht_instance.swapQuadrants(fft_reIm_centered[1]); // zero in the corner
			pixels=fht_instance.FFTHalf2FHT(fft_reIm_centered, size);
		//mask_denoise	
		}
		/// transform to space
		fht_instance.inverseTransform(pixels);
		fht_instance.swapQuadrants(pixels);
		if ((debug>1) ||((DEBUG_LEVEL>2) && (title!=""))) { /* Increase debug level later */ // was 3
			SDFA_INSTANCE.showArrays(pixels, "PSF-"+title);
		}
		return pixels;
	}


/* ======================================================================== */
/* shift (like lateral chromatic aberration) in Bayer component to sensor pixels */

	private  double [] shiftBayerToSensor ( double [] dxy,
			int color,
			int subPixel) {
		return shiftBayerToSensor (dxy[0], dxy[1], color, subPixel);
	}

	private  double [] shiftBayerToSensor ( double dx,
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

	private  double [] shiftSensorToBayer ( double [] dxy,
			int color,
			int subPixel) {
		return shiftSensorToBayer (dxy[0], dxy[1], color, subPixel);
	}
	private  double [] shiftSensorToBayer ( double dx,
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



/* ======================================================================== 
/**
 * Mostly done, need to move where szis\
 * TODO: currently the shift of the PSF during binning is done with the integer steps. If ignoreChromatic - to all colors
 * independently, if it is false - all components are moved in sync, but again - with integer steps. That causes
 * mis-match between the PSF calculated in nearly identical runs (i.e. use the data shifted by 2 pixels) caused by 1 pixel shift.
 * That can be improved if PSF are shifted smoothly (not so easy though). It is probably already handled when averaging PSF - 
 * amplitude and phase is handled separately so shift should be OK.
 * 	
 */
	
	
	double [] combinePSF (double []pixels,         // Square array of pixels with multiple repeated PSF (alternating sign)
			boolean   master,          // force ignoreChromatic
			double[] centerXY,         // coordinates (x,y) of the center point (will update if ignoreChromatic is true)
			double [] centroid_xy,    // RETURNS centroid of the result array (should be small) if ignoreChromatic is true
			double [][] wVectors,    // two wave vectors, lengths in cycles/pixel (pixels match pixel array)
			EyesisAberrations.PSFParameters psfParameters,    // minimal instance contrast to use in binning
			DoubleFHT fht_instance, // provide DoubleFHT instance to save on initializations (or null) // used for sub-pixel shift, null OK
			String title,     // reduce the PSF cell size to this part of the area connecting first negative clones
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
				psfParameters.minContrast,  // minimal contrast of PSF clones
				centerXY,  //  coordinates (x,y) of the center point
				null,  // coordinates of the center of symmetry - not applicable
				1, // pass 1
				title,
				debug);
		//                   true);
		
		if (!master && !psfParameters.ignoreChromatic && psfParameters.centerPSF && (centerXY!=null)){
//			System.out.println("1:pixelsPSF.length="+pixelsPSF.length+" outSize+"+outSize);

			// TODO: Shift +/- 0.5 Pix here {centerXY[0]-Math.round(centerXY[0]),centerXY[1]-Math.round(centerXY[1])}	
			if (fht_instance==null) fht_instance=new DoubleFHT();
//			fht_instance.debug=(centerXY[0]-Math.round(centerXY[0]))<-0.4; // just reducing number
//			double dx=centerXY[0]-Math.round(centerXY[0]);
//			double dy=centerXY[1]-Math.round(centerXY[1]);
//			if (dx<-0.4) SDFA_INSTANCE.showArrays(pixelsPSF.clone(), "before:"+dx+":"+dy);
        
			pixelsPSF=fht_instance.translateSubPixel (
					 pixelsPSF,
					 -(centerXY[0]-Math.round(centerXY[0])),
					 -(centerXY[1]-Math.round(centerXY[1])));
//			fht_instance.debug=false;
//			if (dx<-0.4) SDFA_INSTANCE.showArrays(pixelsPSF.clone(), "after:"+dx+":"+dy);

		}

		double distToNegativeClones=0.5*Math.sqrt(((g[0][0]+g[1][0])*(g[0][0]+g[1][0])+
				(g[0][1]+g[1][1])*(g[0][1]+g[1][1])+
				(g[0][0]-g[1][0])*(g[0][0]-g[1][0])+
				(g[0][1]-g[1][1])*(g[0][1]-g[1][1]))/2.0);
		if (DEBUG_LEVEL>2) {
			System.out.println("distToNegativeClones="+distToNegativeClones+ " gaussWidth="+ distToNegativeClones*psfParameters.smoothSeparate);
		}
		double smoothSigma=distToNegativeClones*psfParameters.smoothSeparate;

		//	double [] smoothPixelsPSF= lowPassGauss(pixelsPSF, smoothSigma, true);
		double [] smoothPixelsPSF= pixelsPSF.clone();
		DoubleGaussianBlur gb=new DoubleGaussianBlur();
		gb.blurDouble(smoothPixelsPSF, outSize, outSize, smoothSigma, smoothSigma, 0.01);

/* find amplitude of smoothed pixel array */
		double smoothMin=0.0;
		double smoothMax=0.0;
		for (i=0;i<smoothPixelsPSF.length;i++) {
			if      (smoothPixelsPSF[i] > smoothMax) smoothMax=smoothPixelsPSF[i];
			else if (smoothPixelsPSF[i] < smoothMin) smoothMin=smoothPixelsPSF[i];
		}
		int [][]  clusterMask = findClusterOnPSF(smoothPixelsPSF, // PSF function, square array (use smooth array)
				-psfParameters.topCenter, // fraction of energy in the pixels to be used (or minimal level if it is negative)
				outSize/2,  // location of a start point, x-coordinate
				outSize/2,  // location of a start point, y-coordinate
				title);
		double [] centroidXY=       calcCentroidFromCenter(pixelsPSF, // use original array (mask from the smoothed one)
				//--centroidXY is in function call arguments
				//centroidXY=            calcCentroidFromCenter(pixelsPSF, // use original array (mask from the smoothed one)
				clusterMask, // integer mask -0 - don't use this pixel, 1 - use it
				psfParameters.topCenter);// subtract level below topCenter*max
		double [] centroidXY_smooth=calcCentroidFromCenter(smoothPixelsPSF, // use smooth - not final, just for clones rejection
				clusterMask, // integer mask -0 - don't use this pixel, 1 - use it
				psfParameters.topCenter);// subtract level below topCenter*max

		if (DEBUG_LEVEL>2) System.out.println("Centroid after first binPSF: x="+IJ.d2s(centroidXY[0],3)+" y="+IJ.d2s(centroidXY[1],3)+" center was at x="+IJ.d2s(centerXY[0],3)+" y="+IJ.d2s(centerXY[1],3));

/* Re-bin results with the new center if ignoreChromatic is true, update centerXY[](shift of the result PSF array) and centroidXY[] (center of the optionally shifted PDF array) */
		if (master || psfParameters.ignoreChromatic) {
			if (centerXY!=null) {
				centerXY[0]+=centroidXY[0];
				centerXY[1]+=centroidXY[1];
			}
			pixelsPSF= binPSF(   pixels,
					g,
					outSize,
					psfParameters.minContrast,  // minimal contrast of PSF clones
					centerXY,  // now includes centroid from the pass 1
					psfParameters.symm180?centroidXY:null,
							2, // pass2
							title,
							debug);
			if (psfParameters.centerPSF && (centerXY!=null)){
//				System.out.println("2:pixelsPSF.length="+pixelsPSF.length+" outSize+"+outSize);
				// TODO: Shift +/- 0.5 Pix here {centerXY[0]-Math.round(centerXY[0]),centerXY[1]-Math.round(centerXY[1])}	
				if (fht_instance==null) fht_instance=new DoubleFHT();
//				fht_instance.debug=(centerXY[0]-Math.round(centerXY[0]))<-0.4; // just reducing number
				pixelsPSF=fht_instance.translateSubPixel (
						 pixelsPSF,
						 -(centerXY[0]-Math.round(centerXY[0])),
						 -(centerXY[1]-Math.round(centerXY[1])));
//				fht_instance.debug=false;
			}
/*  recalculate centroids  */
			smoothPixelsPSF= pixelsPSF.clone();
			gb.blurDouble(smoothPixelsPSF, outSize, outSize, smoothSigma, smoothSigma, 0.01);
			smoothMin=0.0;
			smoothMax=0.0;
			for (i=0;i<smoothPixelsPSF.length;i++) {
				if      (smoothPixelsPSF[i] > smoothMax) smoothMax=smoothPixelsPSF[i];
				else if (smoothPixelsPSF[i] < smoothMin) smoothMin=smoothPixelsPSF[i];
			}
			clusterMask = findClusterOnPSF(smoothPixelsPSF, // PSF function, square array (use smooth array)
					-psfParameters.topCenter, // fraction of energy in the pixels to be used (or minimal level if it is negative)
					outSize/2,  // location of a start point, x-coordinate
					outSize/2,  // location of a start point, y-coordinate
					title);
			centroidXY= calcCentroidFromCenter(pixelsPSF, // use original array (mask from the smoothed one)
					clusterMask, // integer mask -0 - don't use this pixel, 1 - use it
					psfParameters.topCenter);// subtract level below topCenter*max
			// seems it is not used anymore
			centroidXY_smooth=calcCentroidFromCenter(smoothPixelsPSF, // use smooth - not final, just for clones rejection
					clusterMask, // integer mask -0 - don't use this pixel, 1 - use it
					psfParameters.topCenter);// subtract level below topCenter*max
			if (DEBUG_LEVEL>2) System.out.println("Centroid after second binPSF: x="+IJ.d2s(centroidXY[0],3)+" y="+IJ.d2s(centroidXY[1],3)+" center was at x="+IJ.d2s(centerXY[0],3)+" y="+IJ.d2s(centerXY[1],3));

		}
		
		
		
/* compensate center point and/or add center-symmetrical points if enabled */
		double [] rejectedClonesPixels=null;
		double [][] modelPSFVectors={{0.5*(g[0][0]+g[1][0]),0.5*(g[0][1]+g[1][1])},
				{0.5*(g[0][0]-g[1][0]),0.5*(g[0][1]-g[1][1])}};
/********* removed subtraction of clones *****************************************************************/		
		rejectedClonesPixels=pixelsPSF; // Maybe fo the opposite?
		maskClonesPSF(rejectedClonesPixels, // square pixel array where the model PSF is added
				psfParameters.windowFrac, // multiply window by this value
				centroidXY[0], // Center of the remaining single PSF
				centroidXY[1], // same for Y
				modelPSFVectors, // vectors that connect center of PSF with two oppositre sign clones
				psfParameters.useWindow);  // use Hamming window, if false - just cut sharp

		if (psfParameters.wingsEnergy>0.0) {
			rejectedClonesPixels=cutPSFWings (rejectedClonesPixels, // direct PSF function, square array, may be proportionally larger than reversed
					psfParameters.wingsEnergy, // fraction of energy in the pixels to be used
					psfParameters.wingsEllipseScale,
					0.003, // wings_min_mask_threshold, // zero output element if elliptical Gauss mask is below this threshold
					title+"-w");
		}
		double [] sigmas=createSigmasRadius(rejectedClonesPixels, // input square pixel array, preferrably having many exact zeros (they will be skipped)
				psfParameters.sigmaToRadius, // sigma is proportional to the distance from the center
				centroidXY[0], // model PSF center X-coordinate (in pixels[] units, from the center of the array )
				centroidXY[1], // same for Y
				0, // int WOICenterX, // window of interest in pixels[] array - do not generate data outside it
				0, // int WOICenterY, // 
				outSize, //int WOIWidth, reduce later
				outSize); //int WOIHeight)

		double max1=0;
		for (i=0;i<smoothPixelsPSF.length;i++) if (smoothPixelsPSF[i]>max1) max1=smoothPixelsPSF[i];
		double minSigma=0.5;
		double varSigmaTop=1.0 ; //0.7;
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
			if (psfParameters.sigmaToRadius>0.0) {
				float [] floatPixelsSigmas=new float[sigmas.length];
				for (j=0;j<sigmas.length;j++) floatPixelsSigmas[j]=(float) sigmas[j];
				ImageProcessor ip_Sigmas=new FloatProcessor(outSize,outSize);
				ip_Sigmas.setPixels(floatPixelsSigmas);
				ip_Sigmas.resetMinAndMax();
				ImagePlus imp_Sigmas=  new ImagePlus(title+"_Sigmas", ip_Sigmas);
				imp_Sigmas.show();
			}

			System.out.println("title="+title+" center X(pix)="+centroidXY_smooth[0]+"(smooth) center Y(pix)="+centroidXY_smooth[1]+"(smooth)");
			System.out.println("title="+title+" center X(pix)="+centroidXY[0]+"          center Y(pix)="+centroidXY[1]);
		}
		centroid_xy[0]=centroidXY[0];
		centroid_xy[1]=centroidXY[1];
		return  varFilteredPSF;
	}

/* ======================================================================== */

/* ======================================================================== */
/* Trying to remove aliasing artifacts when the decimated (pixel resolution) image is deconvolved with full resolution (sub-pixel resolution)
model pattern. This effect is also easily visible if the decimated model is deconvolved with the same one art full resolution.
Solution is to clone the power spectrum of the full resolution model with the shifts to match oversampling (15 clones for the 4x oversampling),
And add them together (adding also zero frequerncy point - it might be absent on the model) but not include the original (true one) and
use the result to create a rejectiobn mask - if the energy was high, (multiplicative) mask should be zero at those points. */

	private double [] maskAliases (double [] fht, // complex spectrum, [size/2+1][size]
			boolean checker, // checkerboard pattern in the source file (use when filtering)
			int oversample,  // measured array is sampled at 1/oversample frequency than model (will add more parameters later)
			double zerofreq_size,  // add rejection of zero frequency (~2-3pix)
			double sigma,
			double deconvInvert,
			double threshold_high,  // reject completely if energy is above this part of maximal
			double threshold_low,  // leave intact if energy is below this part of maximal
			DoubleFHT fht_instance){ // provide DoubleFHT instance to save on initializations (or null)
		double th=threshold_high*threshold_high;
		double tl=threshold_low*threshold_low;

		int length=fht.length;
		int size=(int) Math.sqrt(fht.length);
		//	double [][] ps=new double [size/2+1][size];
		int i,ix,iy, cloneNx, cloneNy, cloneX, cloneY;
		int cloneStep=size/oversample;
/* generating power spectrum for the high-res complex spectrum, find maximum value and normalize */
		if (fht_instance==null) fht_instance=new DoubleFHT(); // move upstream to reduce number of initializations
		double [] ps=fht_instance.calculateAmplitude2(fht);
		double psMax=0.0;
		for (i=0;i<length; i++) if (psMax<ps[i]) psMax=ps[i];
		double k=1.0/psMax;
		for (i=0;i<length; i++) ps[i]*=k;
		if (DEBUG_LEVEL>2) SDFA_INSTANCE.showArrays(ps, "PS");
/* Add maximum at (0,0) */
		double [] psWithZero=ps;
		if (zerofreq_size>0.0) {
			psWithZero=ps.clone();
			int zs=(int) (4*zerofreq_size);
			int base=size*(size+1)/2;
			k=0.5/(zerofreq_size*zerofreq_size);
			if (zs>=size/2) zs =size/2;
			for (iy=-zs;iy<=zs;iy++) for (ix=-zs; ix <= zs; ix++) {
				psWithZero[base+iy*size+ix]+=Math.exp(-k*(iy*iy+ix*ix));
			}
		}
/* put zero in the center */
		double [] mask=new double [length];
		for (i=0;i<length; i++) mask[i]=0.0;
/* clone spectrums */
		for (iy=0;iy<size;iy++) for (ix=0;ix<size;ix++){
			for (cloneNy=0;cloneNy<oversample;cloneNy++) for (cloneNx=0;cloneNx<oversample;cloneNx++)
				if (((cloneNy!=0) || (cloneNx!=0)) && // not a zero point
						(!checker ||                      // use all if it is not a checkerboard pattren
								(((cloneNx ^ cloneNy) & 1)==0) )) { // remove clones in a checker pattern
					cloneY=(iy+cloneNy*cloneStep)%size;
					cloneX=(ix+cloneNx*cloneStep)%size;
					mask[cloneY*size+cloneX]+=psWithZero[iy*size+ix];
				}			
		}
/* debug show the mask */
		if (DEBUG_LEVEL>2) SDFA_INSTANCE.showArrays(mask, "PS-cloned");
		if (sigma>0) {
			DoubleGaussianBlur gb = new DoubleGaussianBlur();
			gb.blurDouble(mask,size,size,sigma,sigma, 0.01);
			if (DEBUG_LEVEL>2) SDFA_INSTANCE.showArrays(mask, "PS-smooth");
		}

/* make mask of cloned power spectrums */
		double a;
		double k2=deconvInvert*deconvInvert;
		double min=0.01*k2; // less than 1/10 of that value - mask=0.0
		if (DEBUG_LEVEL>2) System.out.println("maskAliases() threshold_high="+threshold_high+" threshold_low="+threshold_low+" th="+th+" tl="+tl+" k2="+k2+" min="+min);
		for (i=0;i<length;i++) {
			if      (mask[i]<tl)  mask[i]=1.0;
			else if (mask[i]>th) mask[i]=0.0;
			else { // make smooth transition
				a=(2.0 * mask[i] - th - tl)/(th - tl);
				mask[i]=0.5*(1.0-a*a*a);
			}
			// now mask out zeros on the ps		
			if (ps[i]<min) mask[i]=0.0;
			else {
				mask[i]*=ps[i]/(ps[i]+k2);
			}
		}
		if (DEBUG_LEVEL>2) SDFA_INSTANCE.showArrays(mask, "mask-all");
		/* zeros are now for FHT - in the top left corner */	
		fht_instance.swapQuadrants(mask);
		return mask;
	}

/* ======================================================================== */
	private double [] binPSF(double [] pixels,
			double [][] g,
			int outSize,
			//		int      decimate,     // sub-pixel decimation 
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
		double threshold=minContrast*contrastAtXY(1, pixels, pixelSize, 0.0, 0.0,  g, contrastCache);
		if (debug)  {
			System.out.println("binPSF title="+title+" g[0][0]="+IJ.d2s(g[0][0],4)+" g[0][1]="+IJ.d2s(g[0][1],4));
			System.out.println("binPSF title="+title+" g[1][0]="+IJ.d2s(g[1][0],4)+" g[1][1]="+IJ.d2s(g[1][1],4));
			System.out.println("  center[0]="+center[0]+"  center[1]="+center[1]);
			//		System.out.println("  decimate="+decimate+"  threshold="+threshold);
			System.out.println("  threshold="+threshold);
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
			//if (debug)  debugPixels[indx]=(int)Math.floor((1+q)/2);
/* see if the point is in the cell of positive or negative OTF instance */
			PSF_sign= (((np + nq) & 1)==0)?1:-1;
/* find x,y coordinates of the center of the cell */
			uc=0.5*(np+nq);
			vc=0.5*(np-nq);
			//xc=g[0][0]*uc + g[1][0]*vc;
			//yc=g[0][1]*uc + g[1][1]*vc;

			yc=-g[0][0]*uc - g[1][0]*vc;
			xc= g[0][1]*uc + g[1][1]*vc;


			//if (debug) debugPixels[indx]=p/2-Math.round(p/2);

/* See if this cell has enough contrast */
			overThreshold=contrastAtXY(PSF_sign,pixels, pixelSize, xc,yc,  g, contrastCache);
			//if (debug) debugPixels[indx]=overThreshold;
			if (overThreshold<threshold) {
				if (debug) debugPixels[indx]=0.0;
				//if (debug) debugPixels[indx]=yc;
			} else {
				//if (debug) debugPixels[indx]=yc;

/* Do binning itself here */
				d=PSF_sign*PSFAtXY(pixels, pixelSize, x,y);

/* map to the segment around 0,0 */        
				dp=p/2-Math.round(p/2);
				dq=q/2-Math.round(q/2);
/* dp, dq are between +/- 0.5 - use them for Hamming windowing -NOT HERE, moved later*/
				du=(dp+dq)/2;
				dv=(dp-dq)/2;

/* bin this point to the center and some (positive) duplicates if enabled */
				for (i=-(multiple/2); i<=(multiple/2); i++) for (j=-(multiple/2); j<=(multiple/2); j++) {
					xr= uv2xy[0][0]*(j+du) + uv2xy[0][1]*(i+dv);
					yr= uv2xy[1][0]*(j+du) + uv2xy[1][1]*(i+dv);
					xr= Math.round(xr-center[0]);
					yr= Math.round(yr-center[1]);
/* does it fit into output array ? */
					if ((yr>=-halfOutSize) && (yr<halfOutSize) && (xr>=-halfOutSize) && (xr<halfOutSize)) {
						outIndex=outSize*(outSize/2+ ((int) yr))+(outSize/2)+((int) xr);
						pixelsPSFCount[outIndex]++;
						pixelsPSF[outIndex]+=d*overThreshold;
						pixelsPSFWeight[outIndex]+=overThreshold;
					}
				}
/* bin this to center-symmetrical point if enabled */
				if (symmXY!=null) {
					for (i=-(multiple/2); i<=(multiple/2); i++) for (j=-(multiple/2); j<=(multiple/2); j++) {
						xr= uv2xy[0][0]*(j+du) + uv2xy[0][1]*(i+dv);
						yr= uv2xy[1][0]*(j+du) + uv2xy[1][1]*(i+dv);
						xr= Math.round(symmXY[0]*2.0-xr-center[0]);
						yr= Math.round(symmXY[1]*2.0-yr-center[1]);
						//does it fit into output array ?
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
					xr= Math.round(xr-center[0]);
					yr= Math.round(yr-center[1]);
					//does it fit into output array ?
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
						xr= Math.round(symmXY[0]*2.0-xr-center[0]);
						yr= Math.round(symmXY[1]*2.0-yr-center[1]);
						//does it fit into output array ?
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
		if ((DEBUG_LEVEL>2) && (pass==1))  SDFA_INSTANCE.showArrays(pixelsPSF, title+"_Used-PSF");
		if (debug) {
			SDFA_INSTANCE.showArrays(debugPixels, title+"_mask_PSF");
			double [] doublePixelsPSFCount=new double [pixelsPSF.length];
			for (j=0;j<doublePixelsPSFCount.length;j++) doublePixelsPSFCount[j]=(double)pixelsPSFCount[j];
			SDFA_INSTANCE.showArrays(doublePixelsPSFCount, title+"_PSF_bin_count");
			SDFA_INSTANCE.showArrays(pixelsPSFWeight,      title+"_PSF_bin_weight");
			double [] doubleContrastCache=new double [contrastCache.length];
			for (j=0;j<doubleContrastCache.length;j++) doubleContrastCache[j]=(double)((contrastCache[j]>=0.0)?contrastCache[j]:-0.00001);
			SDFA_INSTANCE.showArrays(doubleContrastCache,  title+"_ContrastCache");
		}
		return pixelsPSF;
	}


/* ======================================================================== */
/* pixels should be a square array, zero is in the center (/center+0.5 for even dimensions) */
//	private  double [] calcCentroidFromCenter(double [] pixels) {return calcCentroidFromCenter(pixels, (int[]) null, 0.0);}
	private  double [] calcCentroidFromCenter(double [] pixels, // square pixel array
			int[][] mask, // integer mask -0 - don't use this pixel, 1 - use it
			double refLevel) { // subtract this fraction of maximal level from all pixels
		return calcCentroidFromCenter(pixels, convert2d_1d(mask), refLevel);
	}
	private  double [] calcCentroidFromCenter(double [] pixels, // square pixel array
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
		return result;
	}


/* ======================================================================== */
/* zeroes out area outside of the area bound by 4 negative clones (or a fraction of it), either sharp or with Hamming */
	private double [] maskClonesPSF(double [] pixels, // square pixel array where the model PSF is added
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

/* ======================================================================== */
/* create aray (to be used with variableGaussBlurr() ) of per-pixel sigma values for gauss blur, proportional to distance from the specified center */
	private double [] createSigmasRadius (double []pixels, // input square pixel array, preferrably having many exact zeros (they will be skipped)
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
				//        sigma=r*sigmaToRadius;
				//        sigma=r*r/radiusSigma;
				//        sigmas[iy*size+ix]=(r*sigmaToRadius)+1;
				sigmas[iy*size+ix]=(r*sigmaToRadius);
			}
		}


		return sigmas;
	}

/* ======================================================================== */
	private double [] variableGaussBlurr (double []pixels, // input square pixel array, preferrably having many exact zeros (they will be skipped)
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

						for (iky=yk0;iky<=yk1;iky++) {
							gy=gauss[Math.abs(iky)]/2; // Extra /2 because we'll calculate the convolution twice from the [ix,iy] and to [ix,iy]
							for (ikx=xk0;ikx<=xk1;ikx++) {
								index=(iy+iky)*size+ix+ikx;
								g=gy*gauss[Math.abs(ikx)];
								result[index]+=d*g;
								result[iy*size+ix]+=pixels[index]*g;

							}
						}
					}
				}
			}
		}
		return result;
	}
/* ======================================================================== */
/* find ellipse approximating section of the PSF, scale ellipse and use it as a mask to remove PSF far wings */
	private double [] cutPSFWings (double [] psf_pixels, // direct PSF function, square array, may be proportionally larger than reversed
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


/* ======================================================================== */
	private double PSFAtXY(double [] pixels, int size, double x, double y) {
		int ix=(int) Math.round(x);
		int iy=(int) Math.round(y);
		if      (ix <  -size/2) ix=-size/2;
		else if (ix >=  size/2) ix= size/2-1;
		if      (iy <  -size/2) iy=-size/2;
		else if (iy >=  size/2) iy= size/2-1;
		int index=size* (size/2 + iy)+ size/2 + ix;
		if ((index<0) || (index > pixels.length)) {
			System.out.println("PSFAtXY error, x="+IJ.d2s(x,0)+" y="+IJ.d2s(y,0)+ " index="+(size*(size/2 + (int) Math.round(y))+ size/2 + (int) Math.round(x))+ " pixels.length="+pixels.length);
		}
		return pixels[index];
	}
/* ======================================================================== */

	private double contrastAtXY(int sign, double [] pixels, int size, double x, double y, double [][] g, double [] cache) {
		int ir= (int) Math.round(0.2*Math.min(Math.max(Math.abs(g[0][0]),Math.abs(g[1][0])),Math.max(Math.abs(g[0][1]),Math.abs(g[1][1])))); // sample at square 1 1/2x1/2 of the grid "square"

		int ix=(int) Math.round(x);
		int iy=(int) Math.round(y);
		if      (ix <  -size/2) ix=-size/2;
		else if (ix >=  size/2) ix= size/2-1;
		if      (iy <  -size/2) iy=-size/2;
		else if (iy >=  size/2) iy= size/2-1;
		int index= size* (size/2 + iy)+ size/2 + ix;
		//  if ((cache!=null) && (cache[index]>=0)) return sign*cache[index];
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
		return rslt/ir/ir;
	}

/* ======================================================================== */
/* calculates 2x2 matrix that converts two pairs of vectors: u2=M*u1, v2=M*v1*/

/* ======================================================================== */

/* ======================================================================== */

	private  int [] convert2d_1d(int [][] pixels){
		int i,j;
		int width=pixels[0].length;
		int [] rslt=new int[pixels.length*pixels[0].length];
		for (i=0;i<pixels.length;i++) for (j=0;j<width;j++) rslt[i*width+j]=pixels[i][j];
		return rslt;
	}


/* ======================================================================== */
	public double [] cleanupAndReversePSF (double []   psf_pixels,  // input pixels
			EyesisAberrations.InverseParameters inverseParameters, // size (side of square) of direct PSF kernel
			DoubleFHT fht_instance,  // provide DoubleFHT instance to save on initializations (or null)
			String           title   // just for the plot names
	) {
		int size=(int) Math.sqrt(psf_pixels.length);
		double[][][] fft_complex;
		int i,j,ix,iy;
		double a,k,r,r2,k2;

		double [] cpixels=psf_pixels.clone();
		if (fht_instance==null) fht_instance=new DoubleFHT(); // move upstream to reduce number of initializations
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
		int  [][] clusterPS = findClusterOnPS(fft_energy, inverseParameters.otfCutoffEnergy,title);
		double [] ellipse_coeff = findEllipseOnPS(fft_energy, clusterPS, title);
/* create ellipse window using Hamming */
/* TODO: scale radius */
		double [][] ellipseMask=new double [size/2+1][size];
		k2=1/inverseParameters.otfEllipseScale/inverseParameters.otfEllipseScale;
		for (i=0;i<(size/2+1);i++) for (j=0;j<size;j++) {
			iy=(i==size/2)?-i:i;
			ix=(j>=(size/2))?(j-size):j;
			if (iy<0) ix=-ix;
			r2=ellipse_coeff[0]*ix*ix+ellipse_coeff[1]*iy*iy+ellipse_coeff[2]*ix*iy;
			if (inverseParameters.otfEllipseGauss){
				ellipseMask[i][j]=Math.exp(-k2*r2);
			} else {
				r=Math.sqrt(r2)/inverseParameters.otfEllipseScale;
				ellipseMask[i][j]=(r>1.0)?0.0:(0.54+0.46*Math.cos(r*Math.PI));
			}
		}
/* optionally display selection */
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
			ImagePlus imp_ellipse= new ImagePlus(title+"_EL-MASK_"+ inverseParameters.otfCutoffEnergy+"-"+inverseParameters.otfEllipseScale, ip_ellipse);
			imp_ellipse.show();
		}

/* inverse fft_complex */
		if (inverseParameters.invertRange>0.0) {
			/// Invert Z for large values, but make them Z - for small ones. So it will be a mixture of correlation and deconvolution
			//here the targets are round, but what will th\be the colrrect way fo assymmetrical ones?
			/// First - find maximal value
			double fft_max=0;
			for (i=0;i<fft_complex.length; i++) for (j=0;j<fft_complex[0].length;j++) {
				r2=fft_complex[i][j][0]*fft_complex[i][j][0]+fft_complex[i][j][1]*fft_complex[i][j][1];
				if (r2>fft_max) fft_max=r2;
			}
			k=Math.sqrt(fft_max)*inverseParameters.invertRange;
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
/* Convert fft array back to fht array and 
    set fht pixels with new values */
	    pixels=FFTHalf2FHT (fft_complex,size);
/* optionally show the result FHT*/
/* transform to space */
		fht_instance.inverseTransform(pixels);
		fht_instance.swapQuadrants(pixels);
/*   return inverted psf pixels */
		return pixels;
	}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
	private double [] maskReversePSFKernel( double []rpsf_pixels, // reversed psf, square array
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

/* ======================================================================== */
	private  double [] createSigmasFromCenter(
			int               size, // side of square
			double sigma_to_radius, // variable blurring - sigma will be proportional distance from the center
			double    center_sigma, //blurring in the center sigma(r)=sqrt((sigma_to_radius*r)^2)+center_sigma^2)
			double         centerX, // coordinates of the center (0:0 - size/2: size/2)
			double         centerY) {
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

/* ======================================================================== */
/* finds cluster on the PSF (with the center at specidfied point)  by flooding from the specified center, so total energy is cutoff_energy fraction
returns integer array (same dimensions as input) with 1 - selected, 0 - not selected
cutoff_energy: if positive - specifies fraction of total energy, if negative -cutoff_energy is the minimal value of the pixel to be included 
UPDATE: follows gradient from the start point to a local maximum if "cutoff_energy" is negative" */
	private int [][] findClusterOnPSF(
			double []        psf, // PSF function, square array
			double cutoff_energy, // fraction of energy in the pixels to be used
			String         title) {
		int size=(int) Math.sqrt(psf.length);
		return findClusterOnPSF(psf,          // PSF function, square array
				cutoff_energy, // fraction of energy in the pixels to be used
				size/2,        // X0
				size/2,        // Y0
				title);
	}



	private int [][] findClusterOnPSF(
			double []        psf, // PSF function, square array
			double cutoff_energy, // fraction of energy in the pixels to be used (or minimal level if it is negative)
			int           startX,  // location of a start point, x-coordinate
			int           startY,  // location of a start point, y-coordinate
			String         title) {
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
		boolean noThreshold=(cutoff_energy<=0);
		double threshold=full_energy*((cutoff_energy>0)?cutoff_energy:1.0); // no limit for negative values of cutoff_energy
		double minValue=0.0; // no limit if total energy is controlled
		double cluster_energy=0.0;
		int clusterSize=0;
		boolean noNew=true;
		if (cutoff_energy<=0) { // find nearest local maximum following gradient
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
//
if (DEBUG_LEVEL>1)		System.out.println("findClusterOnPSF: full_energy="+full_energy+" minValue="+minValue+" maxValue="+maxValue);
if (DEBUG_LEVEL>1)		System.out.println("findClusterOnPSF: ix="+ix+" iy="+iy);
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
		while ((pixelList.size()>0) &&  (noThreshold || (cluster_energy<threshold) )) { // will break from the loop if  (psf[Index] <minValue)
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
				if (!noThreshold) System.out.println("findClusterOnPSF: - should not get here - no points around >0, and threshold is not reached yet.");
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

/* calculates ellipse (with the center at DC) that interpolates area of the points defined by flooding from the initial center,
so total energy is cutoff_energy fraction
returns {x0,y0,a,b,c} , where a*x^2+b*y^2 + c*x*y=r^2 , so r^2 can be used for a window that removes high far pixels
distribute the whol mass at the ends of short and long ellipse axis
u^2/Ru^2+V^2/Rv^2=1, u=cos(a)*x+sin(a)*y, v=-sin(a)*x+cos(a)*y
c=cos(a), s=sin(a), S0=sum(f(x,y), SX2=sum(f(x,y)*(x-x0)*(x-x0)),SY2=sum(f(x,y)*(y-y0)*(y-y0)), SXY=sum(f(x,y)*(x-x0)*(y-y0))
"effective" squared radius (to be used in Gaussian)
r2= u^2/Ru^2+V^2/Rv^2
r2= 1/Ru^2 * 1/Rv^2 * (x^2*(c^2*Rv^2+s^2*Ru^2)+y^2*(c^2*Ru^2+s^2*Rv^2)+2*x*y*c*s*(Rv^2-Ru^2)

SX2/S0=1/2* ((c*Ru)^2 + (s*Rv)^2)         =1/2*(c^2*Ru^2 + s^2*Rv^2)
SY2/S0=1/2* ((s*Ru)^2 + (c*Rv)^2)         =1/2*(c^2*Rv^2 + s^2*Ru^2)
SXY/S0=1/2* ((c*Ru)*(s*Ru)-(c*Rv)*(s*rv)) =1/2*(c*s*(Ru^2 -Rv^2))

r2= 1/Ru^2 * 1/Rv^2 * (x^2*(2*SY2/S0))+y^2*(2*SX2/S0)-2*2*x*y*(SXY/S0)

SX2/S0+SY2/S0= 1/2*(Ru^2 + Rv^2)
Ru^2+Rv^2= 2*(SX2+SY2)/S0
Ru^2-Rv^2= 2* SXY /S0

Ru^2=(SX2+SY2+SXY)/S0
Rv^2=(SX2+SY2-SXY)/S0

r2= a* x^2*+b*y^2+c*x*y
a=  1/Ru^2 * 1/Rv^2 * (2*SY2/S0)
b=  1/Ru^2 * 1/Rv^2 * (2*SX2/S0)
c= -1/Ru^2 * 1/Rv^2 * (4*SXY/S0)
	 */
	private double [] findEllipseOnPSF(
			double []         psf,   // Point Spread Function (may be off-center)
			int    [][] selection, // 0/1 - selected/not selected
			String          title) {
		int i,j;
		double x,y;
		int size=(int) Math.sqrt(psf.length);
		double SX=0.0;
		double SY=0.0;
		double SX2=0.0;
		double SY2=0.0;
		double SXY=0.0;
		double S0=0.0;
		double d; //,k;
		//	double area=0; // selection area
/* find centyer */

		for (i=0;i<size;i++) {
			y=i-size/2;
			for (j=0;j<size;j++) if (selection[i][j]>0){
				x=j-size/2;
				d=psf[i*size+j];
				S0+=d;
				SX+=x*d;
				SY+=y*d;
				//			area+=1.0;
			}
		}
		double centerX=SX/S0;
		double centerY=SY/S0;
		if (DEBUG_LEVEL>5) {
			//		System.out.println("findEllipseOnPSF: title="+title+" area="+area+" S0="+S0+" SX="+SX+" SY="+SY+" centerX="+centerX+" centerY="+centerY);
			System.out.println("findEllipseOnPSF: title="+title+" S0="+S0+" SX="+SX+" SY="+SY+" centerX="+centerX+" centerY="+centerY);
		}

/* second pass (could all be done in a single) */
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
		/*
Ru^2=(SX2+SY2+SXY)/S0
Rv^2=(SX2+SY2-SXY)/S0

r2= a* x^2*+b*y^2+c*x*y
a=  1/Ru^2 * 1/Rv^2 * (2*SY2/S0)
b=  1/Ru^2 * 1/Rv^2 * (2*SX2/S0)
c= -1/Ru^2 * 1/Rv^2 * (4*SXY/S0)
		 */
		double Ru2=(SX2+SY2+SXY)/S0;
		double Rv2=(SX2+SY2-SXY)/S0;
		double [] result = {centerX,
				centerY,
				1/Ru2 * 1/Rv2 * (2*SY2/S0),
				1/Ru2 * 1/Rv2 * (2*SX2/S0),
				-1/Ru2 * 1/Rv2 * (4*SXY/S0)};
		//	k=Math.PI*Math.PI/(2.0*S0*area*area);
		//	double [] result = {centerX,centerY,k*SY2,k*SX2,-2*k*SXY};
		if (DEBUG_LEVEL>3) {
			System.out.println("findEllipseOnPS: title="+title+" x0="+result[0]+" y0="+result[1]+" a="+result[2]+" b="+result[3]+" c="+result[4]);
		}
		return result;
	}

	// Old version, seem wrong


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* finds cluster (with the center at DC)  by flooding from DC, so total energy is cutoff_energy fraction
returns integer array (same dimensions as input) with 1 - selected, 0 - not selected */
	private int [][] findClusterOnPS(
			double [][]       ps, // half power spectrum, starting from 0.0 (DC)
			double cutoff_energy, // fraction of energy in the pixels to be used
			String         title) {
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

	private double [] findEllipseOnPS(
			double [][]        ps,   // half power spectrum, starting from 0.0 (DC)
			int    [][] selection, // 0/1 - selected/not selected
			String          title) {
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
		//k=Math.PI*Math.PI/(2.0*S0*S0*area*area);
		//double [] result = {k*SY2,k*SX2,2*k*SXY};
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


/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
	public boolean showSimulParametersDialog(SimulationPattern.SimulParameters simulParameters) {
		GenericDialog gd = new GenericDialog("Simulated pattern parameters");
		gd.addNumericField ("Pattern bitmap size",                 simulParameters.patternSize,0); // 512
		gd.addNumericField ("Blur pattern bitmap (sigma):       ", simulParameters.bPatternSigma, 3,5,"pattern cell");
		gd.addNumericField ("Pattern type (0-linear, 1 - curved, ...)", simulParameters.pattern_type, 0);
		gd.addNumericField ("Pattern modifier (i.e. scale)      ", simulParameters.pattern_modifier, 3); // 1.0
		gd.addNumericField ("Frequency 1, X component:          ", simulParameters.freq_x1, 4);
		gd.addNumericField ("Frequency 1, Y component:          ", simulParameters.freq_y1, 4);
		gd.addNumericField ("Phase 1, (radians):                ", simulParameters.phase1, 4);
		gd.addNumericField ("Frequency 2, X component:          ", simulParameters.freq_x2, 4);
		gd.addNumericField ("Frequency 2, Y component:          ", simulParameters.freq_y2, 4);
		gd.addNumericField ("Phase 2, (radians):                ", simulParameters.phase2, 4);
		gd.addNumericField ("Subdivide pixels by:               ", simulParameters.subdiv, 0);
		gd.addNumericField ("Blur pattern (sigma):              ", simulParameters.barraySigma, 3,5,"sensor pix");
		gd.addNumericField ("Photosensitive center part of each simulated pixel", simulParameters.fill, 4); //0.5;  part of the (center) pixel area being "phptosensitive"
		gd.addCheckbox     ("Center pattern for combined greens",  simulParameters.center_for_g2); // true;  // Align pattern to phases for the diagonal (both greens) sub-array
		
		gd.addNumericField ("Smallest fraction to subdivide pixels at simulation", simulParameters.smallestSubPix, 3,5,"sensor pix");
		gd.addNumericField ("Maximal difference of the pattern value in the corners that triggers subdivision", simulParameters.bitmapNonuniforityThreshold, 3);
		gd.addNumericField ("Pattern horizontal (right) offset (debug):", simulParameters.offsetX, 3,5,"sensor pix");
		gd.addNumericField ("Pattern vertical (down) offset (debug):", simulParameters.offsetY, 3,5,"sensor pix");
		
		
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		simulParameters.patternSize= (int) gd.getNextNumber();
		simulParameters.bPatternSigma=     gd.getNextNumber();
		simulParameters.pattern_type=(int) gd.getNextNumber();
		simulParameters.pattern_modifier=  gd.getNextNumber();
		simulParameters.freq_x1=           gd.getNextNumber();
		simulParameters.freq_y1=           gd.getNextNumber();
		simulParameters.phase1=            gd.getNextNumber();
		simulParameters.freq_x2=           gd.getNextNumber();
		simulParameters.freq_y2=           gd.getNextNumber();
		simulParameters.phase2=            gd.getNextNumber();
		simulParameters.subdiv=      (int) gd.getNextNumber();
		simulParameters.barraySigma=       gd.getNextNumber();
		simulParameters.fill=              gd.getNextNumber();
		simulParameters.center_for_g2=     gd.getNextBoolean();

		simulParameters.smallestSubPix=    gd.getNextNumber();
		simulParameters.bitmapNonuniforityThreshold=gd.getNextNumber();
		simulParameters.offsetX=           gd.getNextNumber();
		simulParameters.offsetY=           gd.getNextNumber();
		return true;
	}
/* ======================================================================== */
	public boolean showInverseParametersDialog(EyesisAberrations.InverseParameters inverseParameters) {
		int i;
		GenericDialog gd = new GenericDialog("PSF stack inversion parameters");
		gd.addNumericField("Direct kernel size:",                                             inverseParameters.dSize,                  0); // 32
		gd.addNumericField("Inverted kernel size:",                                           inverseParameters.rSize,                  0); //64
		gd.addNumericField("OTF deconvolution parameter ",                                    inverseParameters.invertRange,            3); //0.01; // when FFT component is less than this fraction of the maximal value, replace 1/z with Z
		gd.addNumericField("OTF cutoff energy (used to determine bounding ellipse)",          inverseParameters.otfCutoffEnergy,        3); //0.6; use frequency points that have inverseParameters.otfCutoffEnergy of the total to determine ellipse for limiting frequency responce
		gd.addNumericField("OTF size of elliptical window relative to cluster size",          inverseParameters.otfEllipseScale,        3); //1.5;  // size of elliptical window relative to the cluster defined by inverseParameters.otfCutoffEnergy
		gd.addCheckbox    ("Use Gauss window (instead of Hamming) for Ellipse mask",          inverseParameters.otfEllipseGauss          ); //   // size of elliptical window relative to the cluster defined by inverseParameters.otfCutoffEnergy
		gd.addNumericField("PSF cutoff energy (used to determine bounding ellipse)",          inverseParameters.psfCutoffEnergy,        3); //0.9;  // Limit result kernel to proportional of the PSF, calculate initial cluster shape by this cutoff energy
		gd.addNumericField("Reversed PSF elliptical mask relative to the direct PSF ellipse", inverseParameters.psfEllipseScale,        3); //1.5;  // size of elliptical window to limuit reverse PSF as proportional to direct one
		gd.addNumericField("Reversed PSF mask thershold (completely ignore far pixels)",      inverseParameters.rpsfMinMaskThreshold,   3); //0.01; // completely zero reversed kernel elements where elliptical mask is below this threshold
		gd.addNumericField("Gaussian blur for individual colors",                             inverseParameters.blurIndividual,         3); //1.8
		gd.addNumericField("Gaussian blur for diagonal greens",                               inverseParameters.blurDiagonal,           3); //1.4
		gd.addNumericField("Gaussian blur for checkerboard greens",                           inverseParameters.blurChecker,            3); //1.4
		gd.addCheckbox    ("Filter direct PSF (blur proportional to radius)",                 inverseParameters.filterDirect);
		gd.addNumericField("Direct PSF: scale variable sigma (in the center) from the uniform one",  inverseParameters.sigmaScaleDirect, 3); //=0.8; // reduce variable sigma in the center from uniuform one
		gd.addNumericField("Direct PSF: sigma to radius ratio (0 to disable)",                inverseParameters.sigmaToRadiusDirect, 3); //0.2- variable blurring - sigma will be proportional distance from the center
		gd.addCheckbox    ("Filter inverted PSF",                                             inverseParameters.filter);
		gd.addNumericField("Reversed PSF: scale variable sigma (in the center) from the uniform one",  inverseParameters.sigmaScale, 3); //=0.8; // reduce variable sigma in the center from uniuform one
		gd.addNumericField("Reversed PSF sigma to radius ratio (0 to disable)",               inverseParameters.sigmaToRadius, 3); //0.2- variable blurring - sigma will be proportional distance from the center
		gd.addMessage("Parameters for Gaussian kernels (used as low resolution/low noise for smooth areas, has the same lateral chromatic correction");
		gd.addNumericField("Gaussian kernels sigma for individual colors",                    inverseParameters.gaussianSigmaIndividual,3); //1.8
		gd.addNumericField("Gaussian kernels sigma for diagonal greens",                      inverseParameters.gaussianSigmaDiagonal,  3); //1.8
		gd.addNumericField("Gaussian kernels sigma for checkerboard greens",                  inverseParameters.gaussianSigmaChecker,   3); //1.4
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		inverseParameters.dSize=1;
		for (i=(int) gd.getNextNumber(); i >1; i>>=1) inverseParameters.dSize <<=1; /* make it to be power of 2 */
		inverseParameters.rSize=1;
		for (i=(int) gd.getNextNumber(); i >1; i>>=1) inverseParameters.rSize <<=1; /* make it to be power of 2 */
		inverseParameters.invertRange=                 gd.getNextNumber();
		inverseParameters.otfCutoffEnergy=             gd.getNextNumber();
		inverseParameters.otfEllipseScale=             gd.getNextNumber();
		inverseParameters.otfEllipseGauss=             gd.getNextBoolean();
		inverseParameters.psfCutoffEnergy=             gd.getNextNumber();
		inverseParameters.psfEllipseScale=             gd.getNextNumber();
		inverseParameters.rpsfMinMaskThreshold=        gd.getNextNumber();
		inverseParameters.blurIndividual=              gd.getNextNumber();
		inverseParameters.blurDiagonal=                gd.getNextNumber();
		inverseParameters.blurChecker=                 gd.getNextNumber();
		inverseParameters.filterDirect=                gd.getNextBoolean();
		inverseParameters.sigmaScaleDirect=            gd.getNextNumber();
		inverseParameters.sigmaToRadiusDirect=         gd.getNextNumber();
		inverseParameters.filter=                      gd.getNextBoolean();
		inverseParameters.sigmaScale=                  gd.getNextNumber();
		inverseParameters.sigmaToRadius=               gd.getNextNumber();
		inverseParameters.gaussianSigmaIndividual=     gd.getNextNumber();
		inverseParameters.gaussianSigmaDiagonal=       gd.getNextNumber();
		inverseParameters.gaussianSigmaChecker=        gd.getNextNumber();
		return true;
	}
/* ======================================================================== */
	public boolean showPSFParametersDialog(EyesisAberrations.PSFParameters psfParameters) {

		GenericDialog gd = new GenericDialog("PSF Parameters");
		gd.addNumericField("Minimal PSF contrast to use",                                        psfParameters.minContrast, 3); // 0.1 - minimal instance contrast to use in binning (compared to the one at [0,0]
		gd.addNumericField("PSF cell reduction from centers of negative clones",                 psfParameters.windowFrac, 3); // 0.9 reduce the PSF cell size to this part of the area connecting first negative clones
		gd.addCheckbox    ("Multiply PSF cell by Hamming window",                                psfParameters.useWindow); //  true;
		gd.addCheckbox    ("Force PSF center- symmetrical (around centroid)",                    psfParameters.symm180); //  true; // make OTF center-symmetrical (around centroid that is defined by lateral chromatic aberration)
		gd.addCheckbox     ("Ignore lateral chromatic aberrations, center PSF",                  psfParameters.ignoreChromatic); // true; // ignore lateral chromatic aberration (center OTF to 0,0)
		gd.addNumericField("PSF separation: low-pass filter width (to PSF half-period) ",        psfParameters.smoothSeparate,  3); // 0.125 low pass filter width (relative to PSF pitch) when separation individual PSF
		gd.addNumericField("PSF separation: threshold to find the PSF maximum",                  psfParameters.topCenter,  3); // 0.75 consider only points above this fraction of the peak to find the centroid
		gd.addNumericField("PSF variable Gauss blurring (farther from center, higher the sigma", psfParameters.sigmaToRadius,3); // 0.4 variable-sigma blurring to reduce high frequencies more for the pixels farther from the PSF center
		gd.addNumericField("PSF wings energy (searching for ellipse approximation)",             psfParameters.wingsEnergy, 3); //  0.8 fraction of energy in the pixels to be used
		gd.addNumericField("PSF wings ellipse scale (multiply PSF by elliptical gaussian)",      psfParameters.wingsEllipseScale, 3);// 2.0 increase wings cutoff ellipse by this from one defined by the  cutoff energy
		gd.addNumericField("Min fraction of the FFT square (weighted) to have defined pattern",  psfParameters.minDefinedArea, 3);
		gd.addCheckbox    ("Approximate pattern grid with a polynomial",                         psfParameters.approximateGrid); // true; // ignore lateral chromatic aberration (center OTF to 0,0)
		gd.addCheckbox    ("Center PSF by modifying phase",                                      psfParameters.centerPSF); // true; // ignore lateral chromatic aberration (center OTF to 0,0)
				
		gd.addNumericField("Bluring power spectrum to remove pattern grid (in pattern base freq)",  psfParameters.mask1_sigma, 3);
		gd.addNumericField("Threshold to supress spectral points not present in the pattern ",  psfParameters.mask1_threshold, 3);
		gd.addNumericField("Sigma for filling the OTF ",                                         psfParameters.gaps_sigma, 3);
		gd.addNumericField("Denoise mask ",                                                      psfParameters.mask_denoise, 3);
		
		
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		psfParameters.minContrast=               gd.getNextNumber();
		psfParameters.windowFrac=                gd.getNextNumber();
		psfParameters.useWindow=                 gd.getNextBoolean();
		psfParameters.symm180=                   gd.getNextBoolean();
		psfParameters.ignoreChromatic=           gd.getNextBoolean();
		psfParameters.smoothSeparate=            gd.getNextNumber();
		psfParameters.topCenter=                 gd.getNextNumber();
		psfParameters.sigmaToRadius=             gd.getNextNumber();
		psfParameters.wingsEnergy=               gd.getNextNumber();
		psfParameters.wingsEllipseScale=         gd.getNextNumber();
		psfParameters.minDefinedArea=            gd.getNextNumber();
		psfParameters.approximateGrid=           gd.getNextBoolean();
		psfParameters.centerPSF=                 gd.getNextBoolean();
		psfParameters.mask1_sigma=               gd.getNextNumber();
		psfParameters.mask1_threshold=           gd.getNextNumber();
		psfParameters.gaps_sigma=                gd.getNextNumber();
		psfParameters.mask_denoise=              gd.getNextNumber();
		return true;
	}
	
/* ======================================================================== */
/* ======================================================================== */
	public boolean showInterpolateParametersDialog(EyesisAberrations.InterpolateParameters interpolateParameters) {
		GenericDialog gd = new GenericDialog("Interpolate kernels parameters");
		gd.addNumericField("Input kernel size",                                                            interpolateParameters.size,     0); // 64
		gd.addNumericField("Interpolation step between original kernels",                                  interpolateParameters.step,       0); //4
		gd.addNumericField("Interpolation add on the top (in output, subdivided steps)",                   interpolateParameters.add_top,     0); //16
		gd.addNumericField("Interpolation add on the left (in output, subdivided steps)",                  interpolateParameters.add_left,    0); //16
		gd.addNumericField("Interpolation add on the right (in output, subdivided steps)",                 interpolateParameters.add_right,   0); //16
		gd.addNumericField("Interpolation add on the bottom (in output, subdivided steps)",                interpolateParameters.add_bottom,  0); //16
		gd.addNumericField("Interpolation: extrapolate margins - 0.0 - duplicate, 1.0 - full extrapolate", interpolateParameters.extrapolate,3); //1.0
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		interpolateParameters.size=       (int) gd.getNextNumber();
		interpolateParameters.step=       (int) gd.getNextNumber();
		interpolateParameters.add_top=    (int) gd.getNextNumber();
		interpolateParameters.add_left=   (int) gd.getNextNumber();
		interpolateParameters.add_right=  (int) gd.getNextNumber();
		interpolateParameters.add_bottom= (int) gd.getNextNumber();
		interpolateParameters.extrapolate=      gd.getNextNumber();
		return true;
	}  
/* ======================================================================== */
/* ======================================================================== */
/**TODO: add variable gaussian filter to direct psf */
	public boolean showOTFFilterParametersDialog(EyesisAberrations.OTFFilterParameters otfFilterParameters) {
		GenericDialog gd = new GenericDialog("OTF Filter parameters");
		gd.addNumericField("Invert deconvolution if less than",                           otfFilterParameters.deconvInvert, 3);
		gd.addNumericField("OTF zero frequency size on power spectrum ",                  otfFilterParameters.zerofreqSize, 3); //2.0;
		gd.addNumericField("OTF smouth PS to generate alias rejection mask (0 - none)",   otfFilterParameters.smoothPS,      3); //2.5 - smooth model PS for rejecting aliases (0 - no smouth, >0 additional Gauss )
		gd.addNumericField("OTF relative high value of PS for rejection mask ",           otfFilterParameters.thresholdHigh, 3); //0.1
		gd.addNumericField("OTF relative low  value of PS for rejection mask ",           otfFilterParameters.thresholdLow,  3); //0.01; // when FFT component is less than this fraction of the maximal value, replace 1/z with Z
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		otfFilterParameters.deconvInvert=        gd.getNextNumber();
		otfFilterParameters.zerofreqSize=        gd.getNextNumber();
		otfFilterParameters.smoothPS=            gd.getNextNumber();
		otfFilterParameters.thresholdHigh=       gd.getNextNumber();
		otfFilterParameters.thresholdLow=        gd.getNextNumber();
		return true;
	}
/* ======================================================================== */
	public boolean showPatternDetectParametersDialog(MatchSimulatedPattern.PatternDetectParameters patternDetectParameters) {
///gaussWidth
		GenericDialog gd = new GenericDialog("Parameters");
		gd.addNumericField("Gaussian width for the window function (<=0 - use Hamming):", patternDetectParameters.gaussWidth, 3);
		gd.addNumericField("Gamma value for pattern frequency measurement:",              patternDetectParameters.corrGamma, 3);
		gd.addNumericField("Sigma value for gauss high-pass pattern filtering:",          patternDetectParameters.corrSigma, 3);// 1.5; // pattern detection: high-pass filter (0.0 - none) gamma(PS)
		gd.addNumericField("Location diff between spectrum and correlation:",             patternDetectParameters.diffSpectrCorr, 0); //2
		gd.addNumericField("Shrink clusters after initial separation (0.0 for sqrt(5*s), negative - absolute size):",   patternDetectParameters.shrinkClusters, 3); // 0.0; //  Shrink clusters by this ratio (remove lowest) after initial separation
		gd.addNumericField("Number secondary maximums to try:",                           patternDetectParameters.multiplesToTry, 0); //4 -  try this number of maximums proportionally farther from 0,0 than the two closest (increase precision)
		gd.addNumericField("Deviation",                                                   patternDetectParameters.deviation, 3); // 1.0 - when looking for maximums - maximal distance from predicted from the lower order one
		gd.addNumericField("Deviation steps",                                             patternDetectParameters.deviationSteps, 0); // 6; maximal iterations when looking for local maximum
		gd.addNumericField("High-pass filter for correlation with the model",             patternDetectParameters.highpass, 3); //1.5 model correlation high-pass filter (relative to pattern fundamental frequency - average of 2)
		gd.addNumericField("Model correlation relative ring width",                       patternDetectParameters.corrRingWidth, 3);   //0.4, ring (around r=0.5 dist to opposite corr) width , center circle r=0.5*patternDetectParameters.corrRingWidth
		gd.addNumericField("Minimal pattern correlation contrast" ,                       patternDetectParameters.minCorrContrast, 3);   //5.0; // Discrimination threshold between good and bad pattern correleation
		gd.addNumericField("Minimal pattern grid period (<=0.0 - do not check)" ,         patternDetectParameters.minGridPeriod, 2,5,"pix");
		gd.addNumericField("Maximal pattern grid period (<=0.0 - do not check)" ,         patternDetectParameters.maxGridPeriod, 2,5,"pix");
		
		
		//			0.0, // minGridPeriod
//		0.0  // maxGridPeriod

		gd.showDialog();
		if (gd.wasCanceled()) return false;
		patternDetectParameters.gaussWidth=               gd.getNextNumber();  //0.4
		patternDetectParameters.corrGamma=               gd.getNextNumber();  //0.2; lower the value - higher harmonics will participate in pattern frequency measurements
		patternDetectParameters.corrSigma=               gd.getNextNumber();  // 1.5; // pattern detection: high-pass filter (0.0 - none) gamma(PS)
		patternDetectParameters.diffSpectrCorr=     (int) gd.getNextNumber();
		patternDetectParameters.shrinkClusters=           gd.getNextNumber(); // 0.5; //  Shrink clusters by this ratio (remove lowest) after initial separation
		patternDetectParameters.multiplesToTry=     (int) gd.getNextNumber();
		patternDetectParameters.deviation=                gd.getNextNumber();
		patternDetectParameters.deviationSteps=    (int) gd.getNextNumber();
		patternDetectParameters.highpass=           gd.getNextNumber();
		patternDetectParameters.corrRingWidth=            gd.getNextNumber();
		patternDetectParameters.minCorrContrast=          gd.getNextNumber();
		patternDetectParameters.minGridPeriod=            gd.getNextNumber();
		patternDetectParameters.maxGridPeriod=            gd.getNextNumber();
		return true;
	}
	public boolean showPatternMinMaxPeriodDialog(MatchSimulatedPattern.PatternDetectParameters patternDetectParameters) {
		///gaussWidth
				GenericDialog gd = new GenericDialog("Min/Max opattern grid period");
				gd.addNumericField("Minimal pattern grid period (<=0.0 - do not check)" ,         patternDetectParameters.minGridPeriod, 2,5,"pix");
				gd.addNumericField("Maximal pattern grid period (<=0.0 - do not check)" ,         patternDetectParameters.maxGridPeriod, 2,5,"pix");
				gd.addMessage ("TODO: Calculate min/max from pattern data and distance");
				gd.showDialog();
				if (gd.wasCanceled()) return false;
				patternDetectParameters.minGridPeriod=            gd.getNextNumber();
				patternDetectParameters.maxGridPeriod=            gd.getNextNumber();
				return true;
			}

	
/* ======================================================================== */
	public boolean showcolorComponentsDialog(EyesisAberrations.ColorComponents colorComponents) {
		GenericDialog gd = new GenericDialog("Color components parameters");
		gd.addCheckbox     ("Calculate correction for RED component",               colorComponents.colorsToCorrect[1]); //  true;
		gd.addCheckbox     ("Calculate correction for BLUE component",              colorComponents.colorsToCorrect[2]); //  true;
		gd.addCheckbox     ("Calculate correction for GREEN components (diagonal)", colorComponents.colorsToCorrect[4]); // false;
		gd.addCheckbox     ("Calculate correction for GREEN components (checker)",  colorComponents.colorsToCorrect[5]); //  true;
		gd.addCheckbox     ("Calculate correction for GREEN1 component (individual)",colorComponents.colorsToCorrect[0]); //  false;
		gd.addCheckbox     ("Calculate correction for GREEN2 component (individual)",colorComponents.colorsToCorrect[3]); //  false;
		gd.addNumericField ("Color component used as base for lateral chromatic (4,5)" ,colorComponents.referenceComponent, 0);   //4; // (change to 5 later) component to calculate lateral chromatic from (0 - G1, 1 - R, 2 - B, 3 - G2,4 - diagonal greens, 5 - checker greens)
		gd.addCheckbox     ("Equalize average values of the two greens in Bayer mosaic",colorComponents.equalizeGreens); //  true; Equalize average values of the two greens in Bayer mosaic
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		colorComponents.colorsToCorrect[1]=          gd.getNextBoolean();
		colorComponents.colorsToCorrect[2]=         gd.getNextBoolean();
		colorComponents.colorsToCorrect[4]=     gd.getNextBoolean();
		colorComponents.colorsToCorrect[5]=      gd.getNextBoolean();
		colorComponents.colorsToCorrect[0]=       gd.getNextBoolean();
		colorComponents.colorsToCorrect[3]=       gd.getNextBoolean();
		colorComponents.referenceComponent=(int)   gd.getNextNumber();
		colorComponents.equalizeGreens=           gd.getNextBoolean();
		if (!colorComponents.colorsToCorrect[colorComponents.referenceComponent]) {
		  for (colorComponents.referenceComponent=5;(colorComponents.referenceComponent>=0) && (!colorComponents.colorsToCorrect[colorComponents.referenceComponent]); colorComponents.referenceComponent--);
		}
		return true;
	}
	
/* ======================================================================== */
	public boolean showShowResultsDialog(ShowResults showResults) {
		GenericDialog gd = new GenericDialog("Show MTF and other results");
		gd.addCheckbox    ("Show PSF",                                                        showResults.showPSF); // show combined PSF kernels (per-color and/or composite, as defined in SHOW_INDIVIDUAL, SHOW_COMPOSITE)
		gd.addCheckbox    ("Show MTF",                                                        showResults.showMTF);  // calculate/show MTF (see notes to showResults.showPSF)
		gd.addCheckbox    ("Show inverted PSF kernels",                                       showResults.showInverted);        // show inverted kernels (unfiltered), same notes
		gd.addCheckbox    ("Show filtered inverted PSF kernels",                              showResults.showFiltered);        // filter and show inverted kernels
		gd.addCheckbox    ("Show gaussian kernels (same centers as in filtered)",             showResults.showGaussians);       // create gaussian kernels with the same centers as inverted ones (low noise, use in low details areas)
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		showResults.showPSF=                           gd.getNextBoolean();
		showResults.showMTF=                           gd.getNextBoolean();
		showResults.showInverted=                      gd.getNextBoolean();
		showResults.showFiltered=                      gd.getNextBoolean();
		showResults.showGaussians=                     gd.getNextBoolean();
		return true;
	}
/* ======================================================================== */
	public boolean showMultiFilePSFDialog(EyesisAberrations.MultiFilePSF multiFilePSF) {
		GenericDialog gd = new GenericDialog("Validate kernels parameters");
		gd.addNumericField("Allowed overexposed pixels (fraction of the area) ",multiFilePSF.overexposedMaxFraction,3); //  0.005; // allowed fraction of the overexposed pixels in the PSF kernel measurement area
		gd.addNumericField("Weight of the PSF from partial file if it has empty neighbors",multiFilePSF.weightOnBorder,3); //  0.5; weight of PSF instances on the border (neglected if there is overlap with not border one)

		gd.addNumericField("Remove partial kernel cell if radius differs from average by this fraction (if there are other files with this cell)",multiFilePSF.radiusDiffLow,3); //  0.01; weight of PSF instances on the border (neglected if there is overlap with not border one)
		gd.addNumericField("Remove partial kernel cell if radius differs from average by this fraction (even if there are no other files with this cell)",multiFilePSF.radiusDiffHigh,3); //  0.01; weight of PSF instances on the border (neglected if there is overlap with not border one)
		
		gd.addNumericField("Center shift (in pixels) addition to the difference relative to radius difference (in pixels)",multiFilePSF.shiftToRadiusContrib,3); //  1.0; // Center shift (in pixels) addition to the difference relative to radius difference (in pixels)
		gd.addNumericField("Increase weight of the \"sharp\" kernels by dividing weights by the radius to this power",multiFilePSF.sharpBonusPower,3); //  1.0; // Center shift (in pixels) addition to the difference relative to radius difference (in pixels)
		
		gd.addNumericField("Discard up to this fraction of samples that have larger radius (i.e. falling on the target seam that may only make PSF larger)",100.0*multiFilePSF.maxFracDiscardWorse,1,5,"%"); // 0.1, discard up to this fraction of samples that have larger radius (i.e. falling on the target seam that may only make PSF larger)
		gd.addNumericField("Continue removing outlayers (combined radius and shift), removing not more that this fraction (including maxFracDiscardWorse)",100.0*multiFilePSF.maxFracDiscardAll,1,5,"%"); // 0.5, continue removing outlayers (combined radius and shift), removing not more that this fraction (including maxFracDiscardWorse)
		
		gd.addNumericField("Cell having 8 around will have thresholds above twice higher than having none when this is set to 1.0, 0.0 - all same thresholds",multiFilePSF.internalBonus,3);
		
		// cell having 8 around will "seem" twice better than having none (radiusDiff* twice higher)
		gd.addCheckbox    ("Fill missing kernels from nearest defined ones",    multiFilePSF.fillMissing);
		gd.addCheckbox    ("Show ellipse parameters",                 multiFilePSF.validateShowEllipse); // 
		gd.addNumericField("Ellipse threshold (fraction of full energy)",       multiFilePSF.validateThreshold,3); //1.0
		gd.addCheckbox    ("Show images indicating frame coverage",             multiFilePSF.showWeights);
		
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		multiFilePSF.overexposedMaxFraction=  gd.getNextNumber();
		multiFilePSF.weightOnBorder=          gd.getNextNumber();
		multiFilePSF.radiusDiffLow=           gd.getNextNumber();
		multiFilePSF.radiusDiffHigh=          gd.getNextNumber();
		
		multiFilePSF.shiftToRadiusContrib=    gd.getNextNumber();
		multiFilePSF.sharpBonusPower=         gd.getNextNumber();
		
		multiFilePSF.maxFracDiscardWorse=0.01*gd.getNextNumber();
		multiFilePSF.maxFracDiscardAll=  0.01*gd.getNextNumber();

		
		
		multiFilePSF.internalBonus=           gd.getNextNumber();
		multiFilePSF.fillMissing=             gd.getNextBoolean();
		multiFilePSF.validateShowEllipse=     gd.getNextBoolean();
		multiFilePSF.validateThreshold=       gd.getNextNumber();
		multiFilePSF.showWeights=             gd.getNextBoolean();
		return true;
	}  

/* ======================================================================== */
	public boolean showConfigureGlobalsDialog() {
		int i;
		GenericDialog gd = new GenericDialog("Parameters");
		gd.addStringField ("Filename prefix:                   ",         JP4_INSTANCE.getTitle(), 20);
		gd.addStringField ("Camera address:                    ",         JP4_INSTANCE.getURL(), 20);
		gd.addNumericField("FFT_Size:",                                   FFT_SIZE, 0);
		gd.addNumericField("FFT_Size for mapping image (areas covered):", MAP_FFT_SIZE, 0); //64 used to find where grid covers the image
		gd.addNumericField("Gauss to window ratio (0 - use Hamming:",     GAUSS_WIDTH, 3);//0.4; //0 - use Hamming window
		gd.addNumericField("FFT Overlap (=FFT_Size recommended):",        FFT_OVERLAP, 0);
		gd.addNumericField("Decimate PSF before binning:",                PSF_SUBPIXEL,    0); // OTF sub-pixel decimation
		gd.addCheckbox    ("Save result to file system",                  PSF_SAVE_FILE); // true;
		gd.addCheckbox    ("Use XML to save/restore parameters",          PROCESS_PARAMETERS.useXML);
		
		gd.addNumericField("Maximal number of concurrent threads",        THREADS_MAX, 0); //   100
		gd.addNumericField("Debug Level:",                                MASTER_DEBUG_LEVEL, 0);
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		JP4_INSTANCE.setTitle(gd.getNextString());
		JP4_INSTANCE.setURL(gd.getNextString());
		FFT_SIZE=4;
		for (i=(int) gd.getNextNumber(); i >4; i>>=1) FFT_SIZE <<=1; /* make FFT_SIZE to be power of 2*/
		MAP_FFT_SIZE=4;
		for (i=(int) gd.getNextNumber(); i >4; i>>=1) MAP_FFT_SIZE <<=1; /* make FFT_SIZE to be power of 2*/
		GAUSS_WIDTH=               gd.getNextNumber();
		FFT_OVERLAP=         (int) gd.getNextNumber();
		PSF_SUBPIXEL=        (int) gd.getNextNumber();
		PSF_SAVE_FILE=             gd.getNextBoolean();
		PROCESS_PARAMETERS.useXML=gd.getNextBoolean();
		THREADS_MAX=         (int) gd.getNextNumber();
		MASTER_DEBUG_LEVEL=  (int) gd.getNextNumber();
		return true;
	}

/* ======================================================================== */
	   public boolean showProcessCalibrationFilesDialog(ProcessCalibrationFilesParameters processCalibrationFilesParameters) {
		    int i;
		    GenericDialog gd = new GenericDialog("Process parameters");
		    gd.addMessage("Paths:");
			gd.addStringField ("Source calibration files directory (with 1-1...3-3 subdirs):    ",          processCalibrationFilesParameters.sourceSuperDirectory, 50);
			gd.addStringField ("Partial PSF kernels directory (with 1-1...3-3 subdirs) - will be created:", processCalibrationFilesParameters.partialKernelsSuperDirectory, 50);
			gd.addStringField ("Result kernels directory (needed for image correction program) - will be created:", processCalibrationFilesParameters.kernelsDirectory, 50);
		    gd.addMessage("Parameters for extacting partial PSF arrays from the source images:");
		    gd.addCheckbox    ("Process source calibration images ",                                        processCalibrationFilesParameters.processSourceImages);
		    gd.addCheckbox    ("Process all channels (if false will use only individually enabled, below)", processCalibrationFilesParameters.processAllChannels);
		    gd.addCheckbox    ("Process channel 1-1",                                                       processCalibrationFilesParameters.processChannels[0]);
		    gd.addCheckbox    ("Process channel 1-2",                                                       processCalibrationFilesParameters.processChannels[1]);
		    gd.addCheckbox    ("Process channel 1-3",                                                       processCalibrationFilesParameters.processChannels[2]);
		    gd.addCheckbox    ("Process channel 2-1",                                                       processCalibrationFilesParameters.processChannels[3]);
		    gd.addCheckbox    ("Process channel 2-2",                                                       processCalibrationFilesParameters.processChannels[4]);
		    gd.addCheckbox    ("Process channel 2-3",                                                       processCalibrationFilesParameters.processChannels[5]);
		    gd.addCheckbox    ("Process channel 3-1",                                                       processCalibrationFilesParameters.processChannels[6]);
		    gd.addCheckbox    ("Process channel 3-2",                                                       processCalibrationFilesParameters.processChannels[7]);
		    gd.addCheckbox    ("Process channel 3-3",                                                       processCalibrationFilesParameters.processChannels[8]);
		    gd.addCheckbox    ("Keep PSF kernels if they exist (false - overwrite)",                        processCalibrationFilesParameters.keepOld);
		    gd.addCheckbox    ("Select individual files to process",                                        processCalibrationFilesParameters.selectFiles);
		    gd.addMessage("Parameters for combining partial PSF kernels and processing them:");
		    
		    gd.addCheckbox    ("Combine partilal PSF kernels",                                              processCalibrationFilesParameters.combinePSFfiles);
		    gd.addCheckbox    ("Interpolate combined PSF kernels",                                          processCalibrationFilesParameters.interpolatePSFkernel);
		    gd.addCheckbox    ("Invert interpolated kernels",                                               processCalibrationFilesParameters.invertKernels);
		    gd.addCheckbox    ("Create gaussian kernels (low-res, same lateral chromatic correction)",      processCalibrationFilesParameters.gaussianKernels);
		    
		    gd.addMessage("Filename extensions and prefixes:");
			gd.addStringField ("Source calibration files extension (no dot):    ",                          processCalibrationFilesParameters.sourceFileExtension, 5);
			gd.addStringField ("PSF kernel files extension (no dot):    ",                                  processCalibrationFilesParameters.kernelFileExtension, 5);
			gd.addStringField ("PSF kernel files prefix:    ",                                              processCalibrationFilesParameters.kernelFilePrefix, 15);
			gd.addStringField ("Combined raw PSF kernel files prefix:    ",                                 processCalibrationFilesParameters.psfRawPrefix, 15);
			gd.addStringField ("Interpolated PSF kernel files prefix:    ",                                 processCalibrationFilesParameters.psfInterpoaltedPrefix, 15);
			gd.addStringField ("Inversed PSF kernel files prefix:    ",                                     processCalibrationFilesParameters.rpsfPrefix, 15);
			gd.addStringField ("Gaussian kernel files prefix:    ",                                         processCalibrationFilesParameters.gaussianPrefix, 15);
		    gd.addMessage("");
		    gd.addCheckbox    ("Save current settings with results",                                        processCalibrationFilesParameters.saveSettings);
	        gd.addCheckbox    ("Use XML format to save/restore settings",                                   processCalibrationFilesParameters.useXML);

			gd.addNumericField("Debug Level:",                                                              MASTER_DEBUG_LEVEL, 0);
		    gd.showDialog();
		    if (gd.wasCanceled()) return false;
		    processCalibrationFilesParameters.sourceSuperDirectory=         gd.getNextString();
		    processCalibrationFilesParameters.partialKernelsSuperDirectory= gd.getNextString();
		    processCalibrationFilesParameters.kernelsDirectory=             gd.getNextString();
		    processCalibrationFilesParameters.processSourceImages=          gd.getNextBoolean();
		    processCalibrationFilesParameters.processAllChannels=           gd.getNextBoolean();
		    for (i=0; i<processCalibrationFilesParameters.processChannels.length;i++)
			     processCalibrationFilesParameters.processChannels[i]=      gd.getNextBoolean();
		    processCalibrationFilesParameters.keepOld=                      gd.getNextBoolean();
		    processCalibrationFilesParameters.selectFiles=                  gd.getNextBoolean();

		    processCalibrationFilesParameters.combinePSFfiles=              gd.getNextBoolean(); 
		    processCalibrationFilesParameters.interpolatePSFkernel=         gd.getNextBoolean();
		    processCalibrationFilesParameters.invertKernels=                gd.getNextBoolean();
		    processCalibrationFilesParameters.gaussianKernels=              gd.getNextBoolean();
		    
			processCalibrationFilesParameters.sourceFileExtension=          gd.getNextString();
			processCalibrationFilesParameters.kernelFileExtension=          gd.getNextString();
			processCalibrationFilesParameters.kernelFilePrefix=             gd.getNextString();
			processCalibrationFilesParameters.psfRawPrefix=                 gd.getNextString();
			processCalibrationFilesParameters.psfInterpoaltedPrefix=        gd.getNextString();
			processCalibrationFilesParameters.rpsfPrefix=                   gd.getNextString();
			processCalibrationFilesParameters.gaussianPrefix=               gd.getNextString();
			processCalibrationFilesParameters.saveSettings=                 gd.getNextBoolean();
			processCalibrationFilesParameters.useXML=		                gd.getNextBoolean();
			MASTER_DEBUG_LEVEL=                                       (int) gd.getNextNumber();
		  return true;
	  }
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
		public boolean showDistortionDialog(MatchSimulatedPattern.DistortionParameters distortionParameters) {
			int i;
			GenericDialog gd = new GenericDialog("Distrortion parameters");
			gd.addNumericField("FFTSize (Initial pattern detection only):",        distortionParameters.FFTSize, 0); // 128
			gd.addNumericField("FFT Gaussian width (relative):",                   distortionParameters.fftGaussWidth, 3);
			gd.addNumericField("Correlation size:",                                distortionParameters.correlationSize, 0); // 64
			gd.addNumericField("Correlation Gauss width (relative):",              distortionParameters.correlationGaussWidth, 3);
			gd.addCheckbox("Keep Gaussian width absolute when increasing FFT size",distortionParameters.absoluteCorrelationGaussWidth);
//			/phaseCorrelationFraction
			//// leave this number of zeros on teh margins of the window (toatal from both sides). If correlationGaussWidth>0 will 
	        // additionally multiply by Hamming
			gd.addNumericField("Leave zeros on the window margins (toatal numbedr from both sides)", distortionParameters.zeros, 0);
			
			gd.addNumericField("Phase correlation modifier (1.0 - phase corr., 0 - just corr.)", distortionParameters.phaseCorrelationFraction, 5);
			gd.addNumericField("Correlation high-pass sigma:",          distortionParameters.correlationHighPassSigma, 3);
			
			gd.addNumericField("Correlation low-pass sigma (fraction of sqrt(2)*Nyquist, lower - more filtering, 0 -none):",distortionParameters.correlationLowPassSigma, 3);
			gd.addNumericField("Correlation maximal offset from predicted:",distortionParameters.correlationMaxOffset, 3);
			gd.addNumericField("Detection ring width (fraction):",      distortionParameters.correlationRingWidth, 3);
			gd.addNumericField("Correlation minimal contrast:",         distortionParameters.correlationMinContrast, 3);
			gd.addNumericField("Correlation minimal contrast for initial search:", distortionParameters.correlationMinInitialContrast, 3);

			gd.addNumericField("Minimal initial pattern cluster size (0 - disable retries)", distortionParameters.minimalPatternCluster, 0);
			gd.addNumericField("Scale minimal contrast if the initial cluster is nonzero but smaller", distortionParameters.scaleMinimalInitialContrast, 3);
			gd.addNumericField("Overlap of FFT areas when searching for pattern", distortionParameters.searchOverlap, 3);
			
			gd.addNumericField("Pattern subdivision:",                  distortionParameters.patternSubdiv, 0);
			gd.addNumericField("Blur pattern bitmap (sigma):       ",   distortionParameters.bPatternSigma, 3,5,"pattern cell");
			gd.addNumericField("Blur pattern (sigma):                ", distortionParameters.barraySigma, 3,5,"sensor pix");
			gd.addNumericField("Correlation weights (around maximum):", distortionParameters.correlationWeightSigma, 3,5,"nodes");
			gd.addNumericField("Correlation radius scale (0 - sharp sigma)", distortionParameters.correlationRadiusScale, 1,3,"sigmas");
			
			gd.addNumericField("Correlation maximal radius to use",      distortionParameters.correlationRadius, 0,1,"pix");
			gd.addNumericField("Correlation maximum calculation threshold", distortionParameters.correlationThreshold*100, 2,5,"%");
			gd.addNumericField("Interpolate correlation (FFT*linear)",  distortionParameters.correlationSubdiv, 0,1,"x");
			gd.addNumericField("Interpolate correlation with FFT",      distortionParameters.correlationFFTSubdiv, 0,1,"x");
			
			gd.addNumericField("Correlation dx (debug)",                distortionParameters.correlationDx, 3);
			gd.addNumericField("Correlation dy (debug)",                distortionParameters.correlationDy, 3);
			gd.addNumericField("Maximal size of the pattern grid (square)", distortionParameters.gridSize, 0);
			gd.addCheckbox    ("Refine correlations",                   distortionParameters.refineCorrelations);
			gd.addCheckbox    ("Use fast correlation on first pass",    distortionParameters.fastCorrelationOnFirstPass);
			gd.addCheckbox    ("Use fast correlation on refine pass",   distortionParameters.fastCorrelationOnFinalPass);

			gd.addCheckbox    ("Average correlation measurements between neighbors (on refine)", distortionParameters.correlationAverageOnRefine);
			gd.addCheckbox    ("Update coordinates of the grid points as they are recalculated (false - then update all at once)", distortionParameters.refineInPlace);
			
			gd.addNumericField("Distance to ortho neighbors (for averaging)", distortionParameters.averageOrthoDist, 3,5,"sensor pix");
			gd.addNumericField("Combined weight of ortho neighbors (fraction of 1.0)", distortionParameters.averageOrthoWeight, 3);
			gd.addNumericField("Distance to diagonal neighbors (for averaging)", distortionParameters.averageDiagDist, 3,5,"sensor pix");
			gd.addNumericField("Combined weight of diagonal neighbors (fraction of 1.0)", distortionParameters.averageDiagWeight, 3);
			gd.addCheckbox    ("Use quadratic extrapolation (false - force linear)", distortionParameters.useQuadratic);
			
			gd.addCheckbox    ("Remove outer (unreliable) layer before extrapolation", distortionParameters.removeLast);
			gd.addNumericField("Number of extrapolated layers of nodes (final stage)", distortionParameters.numberExtrapolated, 0);
			gd.addNumericField("Sigma during final extrapolation stage", distortionParameters.extrapolationSigma, 3,5,"nodes");
			gd.addNumericField("Minimal UV span in correlation window to trigger FFT size increase", distortionParameters.minUVSpan, 3);
			
			gd.addCheckbox    ("Compensate uneven pattern intensity", distortionParameters.flatFieldCorrection);
			gd.addNumericField("Extrapolate pattern intensity map (relative to pattern period)", distortionParameters.flatFieldExtarpolate, 3);
			gd.addNumericField("Blur pattern intensity map (relative to pattern period)", distortionParameters.flatFieldBlur, 3);
			gd.addNumericField("Do not use areas where intensity map is below this part of maximal", distortionParameters.flatFieldMin, 3);
			
			gd.addNumericField("Shrink before extrapolating intensity map (relative to the average grid period)", distortionParameters.flatFieldShrink, 3);
			gd.addNumericField("Expand during extrapolation (relative to the average grid period)", distortionParameters.flatFieldExpand, 3);
			gd.addNumericField("Extrapolation weight effective radius (relative to the average grid period)", distortionParameters.flatFieldSigmaRadius, 3);
			gd.addNumericField("Consider pixels in a square with the side twice this (relative to flatFieldSigmaRadius)", distortionParameters.flatFieldExtraRadius, 3);
			gd.addNumericField("Multiply the average grid period to determine the area for averaging the grig brightness", distortionParameters.averagingAreaScale, 3);
		
			gd.addCheckbox    ("Legacy mode (deprecated)", distortionParameters.legacyMode);
			
			gd.addNumericField("Debug level inside the loop", distortionParameters.loop_debug_level, 0);
			
			
			gd.addNumericField("Debug Level:",                          MASTER_DEBUG_LEVEL, 0);
		    WindowTools.addScrollBars(gd);
			gd.showDialog();
			if (gd.wasCanceled()) return false;

			distortionParameters.FFTSize=1;
			for (i=(int) gd.getNextNumber(); i >1; i>>=1) distortionParameters.FFTSize <<=1; /* make it to be power of 2 */
			distortionParameters.fftGaussWidth=           gd.getNextNumber();
			distortionParameters.correlationSize=1;
			for (i=(int) gd.getNextNumber(); i >1; i>>=1) distortionParameters.correlationSize <<=1; /* make it to be power of 2 */
			distortionParameters.correlationGaussWidth=   gd.getNextNumber();
			distortionParameters.absoluteCorrelationGaussWidth=gd.getNextBoolean();
			distortionParameters.zeros=             (int) gd.getNextNumber();
			distortionParameters.phaseCorrelationFraction=gd.getNextNumber();
			distortionParameters.correlationHighPassSigma=gd.getNextNumber();
			distortionParameters.correlationLowPassSigma= gd.getNextNumber();
			distortionParameters.correlationMaxOffset=    gd.getNextNumber();
			distortionParameters.correlationRingWidth=    gd.getNextNumber();
			distortionParameters.correlationMinContrast=  gd.getNextNumber();
			distortionParameters.correlationMinInitialContrast=  gd.getNextNumber();
			
			distortionParameters.minimalPatternCluster=(int) gd.getNextNumber();
			distortionParameters.scaleMinimalInitialContrast=gd.getNextNumber();
			distortionParameters.searchOverlap=           gd.getNextNumber();
			
			distortionParameters.patternSubdiv=     (int) gd.getNextNumber();
			distortionParameters.bPatternSigma=           gd.getNextNumber();
			distortionParameters.barraySigma=             gd.getNextNumber();
			distortionParameters.correlationWeightSigma=  gd.getNextNumber();
			distortionParameters.correlationRadiusScale=  gd.getNextNumber();
			distortionParameters.correlationRadius= (int) gd.getNextNumber();
			distortionParameters.correlationThreshold= 0.01*gd.getNextNumber();
			distortionParameters.correlationSubdiv= (int) gd.getNextNumber();
			distortionParameters.correlationFFTSubdiv=1;
			for (i=(int) gd.getNextNumber(); i >1; i>>=1) distortionParameters.correlationFFTSubdiv <<=1; /* make it to be power of 2 */
			distortionParameters.correlationDx=           gd.getNextNumber();
			distortionParameters.correlationDy=           gd.getNextNumber();
			distortionParameters.gridSize=          (int) gd.getNextNumber();
			distortionParameters.refineCorrelations=      gd.getNextBoolean();
			distortionParameters.fastCorrelationOnFirstPass=gd.getNextBoolean();
			distortionParameters.fastCorrelationOnFinalPass=gd.getNextBoolean();

			distortionParameters.correlationAverageOnRefine=gd.getNextBoolean();
			distortionParameters.refineInPlace=           gd.getNextBoolean();
			distortionParameters.averageOrthoDist=        gd.getNextNumber();
			distortionParameters.averageOrthoWeight=      gd.getNextNumber();
			distortionParameters.averageDiagDist=         gd.getNextNumber();
			distortionParameters.averageDiagWeight=       gd.getNextNumber();
			distortionParameters.useQuadratic=            gd.getNextBoolean();

			distortionParameters.removeLast=              gd.getNextBoolean();
			distortionParameters.numberExtrapolated=(int) gd.getNextNumber();
			distortionParameters.extrapolationSigma=      gd.getNextNumber();
			distortionParameters.minUVSpan=               gd.getNextNumber();
			distortionParameters.flatFieldCorrection=     gd.getNextBoolean();
			distortionParameters.flatFieldExtarpolate=    gd.getNextNumber();
			distortionParameters.flatFieldBlur=           gd.getNextNumber();
			distortionParameters.flatFieldMin=            gd.getNextNumber();
			distortionParameters.flatFieldShrink=         gd.getNextNumber();
			distortionParameters.flatFieldExpand=         gd.getNextNumber();
			distortionParameters.flatFieldSigmaRadius=    gd.getNextNumber();
			distortionParameters.flatFieldExtraRadius=    gd.getNextNumber();
			distortionParameters.averagingAreaScale=      gd.getNextNumber();
			distortionParameters.legacyMode=              gd.getNextBoolean();
			distortionParameters.loop_debug_level=  (int) gd.getNextNumber();
			MASTER_DEBUG_LEVEL=                     (int) gd.getNextNumber();
			return true;
		}
/* ======================================================================== */
/* ======================================================================== */
	public boolean showGaussianStackDialog(EyesisAberrations.InverseParameters inverseParameters) {
		int i;
		GenericDialog gd = new GenericDialog("Gaussian stack generation parameters");
		gd.addNumericField("Direct kernel size:",                                             inverseParameters.dSize, 0); // 32
		gd.addNumericField("Inverted kernel size:",                                           inverseParameters.rSize, 0); //64

		gd.addNumericField("PSF cutoff energy (used to determine bounding ellipse)",          inverseParameters.psfCutoffEnergy, 3); //0.9;  // Limit result kernel to proportional of the PSF, calculate initial cluster shape by this cutoff energy

		gd.addNumericField("Gaussian blur for individual colors",                             inverseParameters.gaussianSigmaIndividual, 3); //1.8
		gd.addNumericField("Gaussian blur for checkerboard greens",                           inverseParameters.gaussianSigmaChecker,    3); //1.4

		gd.addCheckbox    ("Update ImageJ status",                                            UPDATE_STATUS);
		gd.addNumericField("Debug Level:",                                                    MASTER_DEBUG_LEVEL, 0);
		gd.showDialog();
		if (gd.wasCanceled()) return false;

		inverseParameters.dSize=1;
		for (i=(int) gd.getNextNumber(); i >1; i>>=1) inverseParameters.dSize <<=1; /* make it to be power of 2 */
		inverseParameters.rSize=1;
		for (i=(int) gd.getNextNumber(); i >1; i>>=1) inverseParameters.rSize <<=1; /* make it to be power of 2 */
		inverseParameters.psfCutoffEnergy=             gd.getNextNumber();
		inverseParameters.gaussianSigmaIndividual=       gd.getNextNumber();
		inverseParameters.gaussianSigmaChecker=          gd.getNextNumber();
		UPDATE_STATUS=                 gd.getNextBoolean();
		MASTER_DEBUG_LEVEL=      (int) gd.getNextNumber();
		return true;
	}
/* ======================================================================== */
	public boolean showInverseStackDialog(EyesisAberrations.InverseParameters inverseParameters) {
		int i;
		GenericDialog gd = new GenericDialog("PSF stack inversion parameters");
		gd.addNumericField("Direct kernel size:",                                             inverseParameters.dSize, 0); // 32
		gd.addNumericField("Inverted kernel size:",                                           inverseParameters.rSize, 0); //64

		gd.addNumericField("OTF cutoff energy (used to determine bounding ellipse)",          inverseParameters.otfCutoffEnergy, 3); //0.6; use frequency points that have inverseParameters.otfCutoffEnergy of the total to determine ellipse for limiting frequency responce
		gd.addNumericField("OTF size of elliptical window relative to cluster size",          inverseParameters.otfEllipseScale, 3); //1.5;  // size of elliptical window relative to the cluster defined by inverseParameters.otfCutoffEnergy
		gd.addCheckbox    ("Use Gauss window (instead of Hamming) for Ellipse mask",          inverseParameters.otfEllipseGauss); //   // size of elliptical window relative to the cluster defined by inverseParameters.otfCutoffEnergy
		gd.addNumericField("OTF deconvolution parameter ",                                    inverseParameters.invertRange, 3); //0.01; // when FFT component is less than this fraction of the maximal value, replace 1/z with Z
		gd.addNumericField("PSF cutoff energy (used to determine bounding ellipse)",          inverseParameters.psfCutoffEnergy, 3); //0.9;  // Limit result kernel to proportional of the PSF, calculate initial cluster shape by this cutoff energy
		gd.addNumericField("Reversed PSF elliptical mask relative to the direct PSF ellipse", inverseParameters.psfEllipseScale, 3); //1.5;  // size of elliptical window to limuit reverse PSF as proportional to direct one
		gd.addNumericField("Reversed PSF mask thershold (completely ignore far pixels)",      inverseParameters.rpsfMinMaskThreshold, 3); //0.01; // completely zero reversed kernel elements where elliptical mask is below this threshold

		gd.addCheckbox    ("Filter inverted PSF",                                             inverseParameters.filter);

		gd.addNumericField("Reversed PSF sigma to radius ratio (0 to disable)",               inverseParameters.sigmaToRadius, 3); //0.2- variable blurring - sigma will be proportional distance from the center
		gd.addNumericField("Reversed PSF: scale variable sigma (in the center) from the uniform one",  inverseParameters.sigmaScale, 3); //=0.8; // reduce variable sigma in the center from uniuform one

		gd.addNumericField("Gaussian blur for individual colors",                             inverseParameters.blurIndividual, 3); //1.8
		gd.addNumericField("Gaussian blur for checkerboard greens",                           inverseParameters.blurChecker,    3); //1.4


		gd.addCheckbox    ("Update ImageJ status",                                            UPDATE_STATUS);
		gd.addNumericField("Maximal number of concurrent threads",                            THREADS_MAX, 0); //   100
		gd.addNumericField("Debug Level:",                                                    MASTER_DEBUG_LEVEL, 0);
		gd.showDialog();
		if (gd.wasCanceled()) return false;

		inverseParameters.dSize=1;
		for (i=(int) gd.getNextNumber(); i >1; i>>=1) inverseParameters.dSize <<=1; /* make it to be power of 2 */
		inverseParameters.rSize=1;
		for (i=(int) gd.getNextNumber(); i >1; i>>=1) inverseParameters.rSize <<=1; /* make it to be power of 2 */
		inverseParameters.otfCutoffEnergy=             gd.getNextNumber();
		inverseParameters.otfEllipseScale=             gd.getNextNumber();
		inverseParameters.otfEllipseGauss=             gd.getNextBoolean();
		inverseParameters.invertRange=                 gd.getNextNumber();
		inverseParameters.psfCutoffEnergy=             gd.getNextNumber();
		inverseParameters.psfEllipseScale=             gd.getNextNumber();
		inverseParameters.rpsfMinMaskThreshold=        gd.getNextNumber();
		inverseParameters.filter=                      gd.getNextBoolean();
		inverseParameters.sigmaToRadius=               gd.getNextNumber();
		inverseParameters.sigmaScale=                  gd.getNextNumber();

		inverseParameters.blurIndividual=              gd.getNextNumber();
		inverseParameters.blurChecker=                 gd.getNextNumber();
		UPDATE_STATUS=                                 gd.getNextBoolean();
		THREADS_MAX=                             (int) gd.getNextNumber();
		MASTER_DEBUG_LEVEL=                      (int) gd.getNextNumber();
		return true;
	}

/* ======================================================================== */
	//public static EyesisAberrations.InterpolateParameters INTERPOLATE= new EyesisAberrations.InterpolateParameters (

	public boolean showValidateKernelsDialog(
			EyesisAberrations.InterpolateParameters interpolateParameters,
			EyesisAberrations.MultiFilePSF multiFilePSF
			) {
		GenericDialog gd = new GenericDialog("Validate kernels parameters");
		gd.addNumericField("Kernel size",                                 interpolateParameters.size,     0); // 64
		gd.addCheckbox    ("Calculate/show ellipse parameters",           multiFilePSF.validateShowEllipse); // 
		gd.addNumericField("Ellipse threshold (fraction of full energy)", multiFilePSF.validateThreshold,3); //1.0
		gd.addCheckbox    ("Save result to file system",                  PSF_SAVE_FILE); // true;
		gd.addCheckbox    ("Update ImageJ status",                        UPDATE_STATUS);
		gd.addNumericField("Debug Level:",                                MASTER_DEBUG_LEVEL,      0);
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		interpolateParameters.size= (int) gd.getNextNumber();
		multiFilePSF.validateShowEllipse= gd.getNextBoolean();
		multiFilePSF.validateThreshold=   gd.getNextNumber();
		PSF_SAVE_FILE=                    gd.getNextBoolean();
		UPDATE_STATUS=                    gd.getNextBoolean();
		MASTER_DEBUG_LEVEL=         (int) gd.getNextNumber();
		return true;
	}  
/* ======================================================================== */
	public boolean showInterpolateKernelsDialog(EyesisAberrations.InterpolateParameters interpolateParameters) {
		GenericDialog gd = new GenericDialog("Interpolate kernels parameters");
		gd.addNumericField("Input kernel size",                                                            interpolateParameters.size,     0); // 64
		gd.addNumericField("Interpolation step between original kernels",                                  interpolateParameters.step,       0); //4
		gd.addNumericField("Interpolation add on the top (in output, subdivided steps)",                   interpolateParameters.add_top,     0); //16
		gd.addNumericField("Interpolation add on the left (in output, subdivided steps)",                  interpolateParameters.add_left,    0); //16
		gd.addNumericField("Interpolation add on the right (in output, subdivided steps)",                 interpolateParameters.add_right,   0); //16
		gd.addNumericField("Interpolation add on the bottom (in output, subdivided steps)",                interpolateParameters.add_bottom,  0); //16
		gd.addNumericField("Interpolation: extrapolate margins - 0.0 - duplicate, 1.0 - full extrapolate", interpolateParameters.extrapolate,3); //1.0
		gd.addCheckbox    ("Update ImageJ status",                                                         UPDATE_STATUS);
		gd.addCheckbox    ("Save result to file system",                                                    PSF_SAVE_FILE); // true;

		gd.addNumericField("Debug Level:",                                                                 MASTER_DEBUG_LEVEL,      0);
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		interpolateParameters.size=       (int) gd.getNextNumber();
		interpolateParameters.step=       (int) gd.getNextNumber();
		interpolateParameters.add_top=    (int) gd.getNextNumber();
		interpolateParameters.add_left=   (int) gd.getNextNumber();
		interpolateParameters.add_right=  (int) gd.getNextNumber();
		interpolateParameters.add_bottom= (int) gd.getNextNumber();
		interpolateParameters.extrapolate=      gd.getNextNumber();
		PSF_SAVE_FILE=                          gd.getNextBoolean();
		UPDATE_STATUS=                          gd.getNextBoolean();
		MASTER_DEBUG_LEVEL=               (int) gd.getNextNumber();
		return true;
	}  
/* ======================================================================== */

	public boolean showSimulDialog(SimulationPattern.SimulParameters simulParameters) {
		GenericDialog gd = new GenericDialog("Simulated pattern parameters");
		gd.addNumericField ("Pattern type (0-linear, 1 - curved, ...)", simulParameters.pattern_type, 0);
		gd.addNumericField ("Pattern modifier (i.e. scale)      ", simulParameters.pattern_modifier, 3); // 1.0
		gd.addNumericField ("Frequency 1, X component:          ", simulParameters.freq_x1, 4);
		gd.addNumericField ("Frequency 1, Y component:          ", simulParameters.freq_y1, 4);
		gd.addNumericField ("Phase 1, (radians):                ", simulParameters.phase1, 4);
		gd.addNumericField ("Frequency 2, X component:          ", simulParameters.freq_x2, 4);
		gd.addNumericField ("Frequency 2, Y component:          ", simulParameters.freq_y2, 4);
		gd.addNumericField ("Phase 2, (radians):                ", simulParameters.phase2, 4);

		gd.addCheckbox     ("Center pattern for combined greens",  simulParameters.center_for_g2); // true;  // Align pattern to phases for the diagonal (both greens) sub-array
		gd.addNumericField ("Subdivide pixels by:            ",    simulParameters.subdiv, 0);
		gd.addNumericField ("Photosensitive center part of each simulated pixel", simulParameters.fill, 4); //0.5;  part of the (center) pixel area being "phptosensitive"
		gd.addNumericField("Debug Level:",      MASTER_DEBUG_LEVEL, 0);
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		simulParameters.pattern_type=(int) gd.getNextNumber();
		simulParameters.pattern_modifier=  gd.getNextNumber();
		simulParameters.freq_x1=           gd.getNextNumber();
		simulParameters.freq_y1=           gd.getNextNumber();
		simulParameters.phase1=            gd.getNextNumber();
		simulParameters.freq_x2=           gd.getNextNumber();
		simulParameters.freq_y2=           gd.getNextNumber();
		simulParameters.phase2=            gd.getNextNumber();
		simulParameters.center_for_g2=     gd.getNextBoolean();
		simulParameters.subdiv=      (int) gd.getNextNumber();
		simulParameters.fill=              gd.getNextNumber();
		MASTER_DEBUG_LEVEL=          (int) gd.getNextNumber();
		return true;
	}


/* ======================================================================== */
	public boolean showRPSFDialog(
			EyesisAberrations.InverseParameters inverseParameters,
			ShowResults showResults) {
		int i;
		GenericDialog gd = new GenericDialog("PSF reversal parameters");
		gd.addNumericField("OTF cutoff energy (used to determine bounding ellipse)",          inverseParameters.otfCutoffEnergy, 3); //0.6; use frequency points that have inverseParameters.otfCutoffEnergy of the total to determine ellipse for limiting frequency responce
		gd.addNumericField("OTF size of elliptical window relative to cluster size",          inverseParameters.otfEllipseScale, 3); //1.5;  // size of elliptical window relative to the cluster defined by inverseParameters.otfCutoffEnergy
		gd.addCheckbox    ("Use Gauss window (instead of Hamming) for Ellipse mask",          inverseParameters.otfEllipseGauss); //   // size of elliptical window relative to the cluster defined by inverseParameters.otfCutoffEnergy
		gd.addNumericField("OTF deconvolution parameter ",                                    inverseParameters.invertRange, 3); //0.01; // when FFT component is less than this fraction of the maximal value, replace 1/z with Z
		gd.addNumericField("Direct PSF kernel size (side of square array to be stored) ",     inverseParameters.dSize,  0); // 64 -kernel (to be stored) size (per color component)
		gd.addNumericField("PSF cutoff energy (used to determine bounding ellipse)",          inverseParameters.psfCutoffEnergy, 3); //0.9;  // Limit result kernel to proportional of the PSF, calculate initial cluster shape by this cutoff energy
		gd.addNumericField("Reversed PSF elliptical mask relative to the direct PSF ellipse", inverseParameters.psfEllipseScale, 3); //1.5;  // size of elliptical window to limuit reverse PSF as proportional to direct one
		gd.addNumericField("Reversed PSF mask thershold (completely ignore far pixels)",      inverseParameters.rpsfMinMaskThreshold, 3); //0.01; // completely zero reversed kernel elements where elliptical mask is below this threshold
		gd.addNumericField("Reversed PSF sigma to radius ratio (0 to disable)",               inverseParameters.sigmaToRadius, 3); //0.2- variable blurring - sigma will be proportional distance from the center
		gd.addNumericField("Reversed PSF: scale variable sigma (in the center) from the uniform one",  inverseParameters.sigmaScale, 3); //=0.8; // reduce variable sigma in the center from uniuform one
		gd.addNumericField("Gaussian blur for individual colors",                             inverseParameters.blurIndividual, 3); //1.8
		gd.addNumericField("Gaussian blur for diagonal greens",                               inverseParameters.blurDiagonal,   3); //2.0
		gd.addNumericField("Gaussian blur for checkerboard greens",                           inverseParameters.blurChecker,    3); //1.4
		gd.addNumericField("Reversed PSF  kernel size",                                       inverseParameters.rSize, 0); //   32 - size of deconvolution kernel
		gd.addCheckbox    ("Update ImageJ status",                                            UPDATE_STATUS);
		gd.addCheckbox    ("Show PSF",                                                        showResults.showPSF); // show combined PSF kernels (per-color and/or composite, as defined in SHOW_INDIVIDUAL, SHOW_COMPOSITE)
		gd.addCheckbox    ("Show MTF",                                                        showResults.showMTF);  // calculate/show MTF (see notes to showResults.showPSF)
		gd.addCheckbox    ("Show inverted PSF kernels",                                       showResults.showInverted);        // show inverted kernels (unfiltered), same notes
		gd.addCheckbox    ("Show filtered inverted PSF kernels",                              showResults.showFiltered);        // filter and show inverted kernels
		gd.addCheckbox    ("Show gaussian kernels (same centers as in filtered)",             showResults.showGaussians);       // create gaussian kernels with the same centers as inverted ones (low noise, use in low details areas)
		gd.addCheckbox    ("Show debug images as stacks",                                     SHOW_AS_STACKS);       // true
		gd.addNumericField("Debug Level:",                                                    MASTER_DEBUG_LEVEL, 0);
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		inverseParameters.otfCutoffEnergy=             gd.getNextNumber();
		inverseParameters.otfEllipseScale=             gd.getNextNumber();
		inverseParameters.otfEllipseGauss=             gd.getNextBoolean();
		inverseParameters.invertRange=                 gd.getNextNumber();
		inverseParameters.dSize=                 (int) gd.getNextNumber();
		inverseParameters.psfCutoffEnergy=             gd.getNextNumber();
		inverseParameters.psfEllipseScale=             gd.getNextNumber();
		inverseParameters.rpsfMinMaskThreshold=        gd.getNextNumber();
		inverseParameters.sigmaToRadius=               gd.getNextNumber();
		inverseParameters.sigmaScale=                  gd.getNextNumber();
		inverseParameters.blurIndividual=              gd.getNextNumber();
		inverseParameters.blurDiagonal=                gd.getNextNumber();
		inverseParameters.blurChecker=                 gd.getNextNumber();
		inverseParameters.rSize=1;
		for (i=(int) gd.getNextNumber(); i >1; i>>=1) inverseParameters.rSize <<=1; /* make it to be power of 2 */
		UPDATE_STATUS=                                 gd.getNextBoolean();
		showResults.showPSF=                           gd.getNextBoolean();
		showResults.showMTF=                           gd.getNextBoolean();
		showResults.showInverted=                      gd.getNextBoolean();
		showResults.showFiltered=                      gd.getNextBoolean();
		showResults.showGaussians=                     gd.getNextBoolean();
		SHOW_AS_STACKS=                                gd.getNextBoolean();
		MASTER_DEBUG_LEVEL=                      (int) gd.getNextNumber();
		return true;
	}


	public ImagePlus acquireAndDisplay(ImagePlus imp_src) {
		ImagePlus imp=JP4_INSTANCE.openURL(imp_src);
		//     if (imp!=null) measureAndDisplay(imp);
		return imp;
	}
/* ======================================================================== */
	public boolean showPSFDialog(EyesisAberrations.PSFParameters psfParameters,
			                     EyesisAberrations.InverseParameters inverseParameters,
			                     EyesisAberrations.OTFFilterParameters otfFilterParameters,
			                     SimulationPattern.SimulParameters  simulParameters
			                     ) {

		GenericDialog gd = new GenericDialog("PSF measurement parameters");
		gd.addNumericField("Subdivide simulation pixels during generation by:",           simulParameters.subdiv, 0);
		gd.addNumericField("Photosensitive center part of each simulated pixel",          simulParameters.fill, 3); //0.5;  part of the (center) pixel area being "phptosensitive"
		gd.addNumericField("Decimate PSF before binning:",                                PSF_SUBPIXEL,    0); // OTF sub-pixel decimation
		gd.addNumericField("Minimal PSF contrast to use",                                 psfParameters.minContrast, 3); // 0.1 - minimal instance contrast to use in binning (compared to the one at [0,0]
		gd.addNumericField("PSF cell reduction from centers of negative clones",          psfParameters.windowFrac, 3); // 0.9 reduce the PSF cell size to this part of the area connecting first negative clones
		gd.addCheckbox    ("Multiply PSF cell by Hamming window",                         psfParameters.useWindow); //  true;
		gd.addCheckbox    ("Force PSF center- symmetrical (around centroid)",             psfParameters.symm180); //  true; // make OTF center-symmetrical (around centroid that is defined by lateral chromatic aberration)
		gd.addNumericField("PSF separation: low-pass filter width (to PSF half-period) ", psfParameters.smoothSeparate,  3); // 0.125 low pass filter width (relative to PSF pitch) when separation individual PSF
		gd.addNumericField("PSF separation: threshold to find the PSF maximum",           psfParameters.topCenter,  3); // 0.75 consider only points above this fraction of the peak to find the centroid
		gd.addNumericField("PSF variable Gauss blurring (farther from center, higher the sigma",    psfParameters.sigmaToRadius,3); // 0.4 variable-sigma blurring to reduce high frequencies more for the pixels farther from the PSF center
		gd.addNumericField("PSF wings energy (searching for ellipse approximation)",      psfParameters.wingsEnergy, 3); //  0.8 fraction of energy in the pixels to be used
		gd.addNumericField("PSF wings ellipse scale (multiply PSF by elliptical gaussian)",psfParameters.wingsEllipseScale, 3);// 2.0 increase wings cutoff ellipse by this from one defined by the  cutoff energy
		gd.addCheckbox     ("Ignore lateral chromatic aberrations, center PSF",           psfParameters.ignoreChromatic); // true; // ignore lateral chromatic aberration (center OTF to 0,0)

		gd.addNumericField("OTF cutoff energy (used to determine bounding ellipse)",      inverseParameters.otfCutoffEnergy, 3); //0.6; use frequency points that have inverseParameters.otfCutoffEnergy of the total to determine ellipse for limiting frequency responce
		gd.addNumericField("OTF size of elliptical window relative to cluster size",      inverseParameters.otfEllipseScale, 3); //1.5;  // size of elliptical window relative to the cluster defined by inverseParameters.otfCutoffEnergy
		gd.addCheckbox    ("Use Gauss window (instead of Hamming) for Ellipse mask",      inverseParameters.otfEllipseGauss); //   // size of elliptical window relative to the cluster defined by inverseParameters.otfCutoffEnergy
		gd.addNumericField("OTF deconvolution parameter ",                                inverseParameters.invertRange, 3); //0.01; // when FFT component is less than this fraction of the maximal value, replace 1/z with Z
		gd.addNumericField("Invert deconvolution if less than",                           otfFilterParameters.deconvInvert, 3);
		gd.addNumericField("OTF zero frequency size on power spectrum ",                  otfFilterParameters.zerofreqSize, 3); //2.0;
		gd.addNumericField("OTF smouth PS to generate alias rejection mask (0 - none)",   otfFilterParameters.smoothPS,      3); //2.5 - smooth model PS for rejecting aliases (0 - no smouth, >0 additional Gauss )
		gd.addNumericField("OTF relative high value of PS for rejection mask ",           otfFilterParameters.thresholdHigh, 3); //0.1
		gd.addNumericField("OTF relative low  value of PS for rejection mask ",           otfFilterParameters.thresholdLow,  3); //0.01; // when FFT component is less than this fraction of the maximal value, replace 1/z with Z


		gd.addNumericField("PSF kernel size (side of square array to be stored) ",            inverseParameters.dSize,  0); // 64 -kernel (to be stored) size (per color component)
		gd.addNumericField("PSF cutoff energy (used to determine bounding ellipse)",          inverseParameters.psfCutoffEnergy, 3); //0.9;  // Limit result kernel to proportional of the PSF, calculate initial cluster shape by this cutoff energy
		gd.addNumericField("Reversed PSF elliptical mask relative to the direct PSF ellipse", inverseParameters.psfEllipseScale, 3); //1.5;  // size of elliptical window to limuit reverse PSF as proportional to direct one
		gd.addNumericField("Reversed PSF mask thershold (completely ignore far pixels)",      inverseParameters.rpsfMinMaskThreshold, 3); //0.01; // completely zero reversed kernel elements where elliptical mask is below this threshold
		gd.addCheckbox     ("Save result to file system",                                     PSF_SAVE_FILE); // true;
		gd.addNumericField("Debug Level:",                                                    MASTER_DEBUG_LEVEL, 0);
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		simulParameters.subdiv=            (int) gd.getNextNumber();
		simulParameters.fill=                    gd.getNextNumber();
		PSF_SUBPIXEL=                      (int) gd.getNextNumber();
		psfParameters.minContrast=               gd.getNextNumber();
		psfParameters.windowFrac=                gd.getNextNumber();
		psfParameters.useWindow=                 gd.getNextBoolean();
		psfParameters.symm180=                   gd.getNextBoolean();
		psfParameters.smoothSeparate=            gd.getNextNumber();
		psfParameters.topCenter=                 gd.getNextNumber();
		psfParameters.sigmaToRadius=             gd.getNextNumber();
		psfParameters.wingsEnergy=               gd.getNextNumber();
		psfParameters.wingsEllipseScale=         gd.getNextNumber();
		psfParameters.ignoreChromatic=           gd.getNextBoolean();
		inverseParameters.otfCutoffEnergy=       gd.getNextNumber();
		inverseParameters.otfEllipseScale=       gd.getNextNumber();
		inverseParameters.otfEllipseGauss=       gd.getNextBoolean();
		inverseParameters.invertRange=           gd.getNextNumber();
		otfFilterParameters.deconvInvert=        gd.getNextNumber();
		otfFilterParameters.zerofreqSize=        gd.getNextNumber();
		otfFilterParameters.smoothPS=            gd.getNextNumber();
		otfFilterParameters.thresholdHigh=       gd.getNextNumber();
		otfFilterParameters.thresholdLow=        gd.getNextNumber();
		inverseParameters.dSize=           (int) gd.getNextNumber();
		inverseParameters.psfCutoffEnergy=       gd.getNextNumber();
		inverseParameters.psfEllipseScale=       gd.getNextNumber();
		inverseParameters.rpsfMinMaskThreshold=  gd.getNextNumber();
		PSF_SAVE_FILE=                 gd.getNextBoolean();
		MASTER_DEBUG_LEVEL=      (int) gd.getNextNumber();
		return true;
	}

	public boolean showConfigureDialog(
			EyesisAberrations.OTFFilterParameters     otfFilterParameters,
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			EyesisAberrations.ColorComponents         colorComponents,
			EyesisAberrations.InverseParameters       inverseParameters,
			EyesisAberrations.MultiFilePSF            multiFilePSF) {
		int i;
		GenericDialog gd = new GenericDialog("Parameters");
		gd.addStringField ("Filename prefix:                   ", JP4_INSTANCE.getTitle(), 20);
		gd.addStringField ("Camera address:                    ", JP4_INSTANCE.getURL(), 20);
		gd.addNumericField("FFT_Size:",                           FFT_SIZE, 0);
		gd.addNumericField("FFT_Size for mapping image (areas covered):", MAP_FFT_SIZE, 0); //64 used to find where grid covers the image
		gd.addNumericField("Gauss to window ratio (0 - use Hamming:",GAUSS_WIDTH, 3);//0.4; //0 - use Hamming window
		gd.addNumericField("FFT Overlap (=FFT_Size recommended):", FFT_OVERLAP, 0);
		gd.addNumericField("Invert deconvolution if less than",   otfFilterParameters.deconvInvert, 3);
		gd.addNumericField("Gamma value for pattern frequency measurement:",   patternDetectParameters.corrGamma, 3);
		gd.addNumericField("Sigma value for gauss high-pass pattern filtering:",   patternDetectParameters.corrSigma, 3);// 1.5; // pattern detection: high-pass filter (0.0 - none) gamma(PS)
		gd.addNumericField("Location diff between spectrum and correlation:",  patternDetectParameters.diffSpectrCorr, 0); //2
		gd.addNumericField("Shrink clusters after initial separation (0.0 for sqrt(5*s), negative - absolute size):",   patternDetectParameters.shrinkClusters, 3); // 0.0; //  Shrink clusters by this ratio (remove lowest) after initial separation
		gd.addNumericField("Number secondary maximums to try:",   patternDetectParameters.multiplesToTry, 0); //4 -  try this number of maximums proportionally farther from 0,0 than the two closest (increase precision)
		gd.addNumericField("Deviation",                           patternDetectParameters.deviation, 3); // 1.0 - when looking for maximums - maximal distance from predicted from the lower order one
		gd.addNumericField("Deviation steps",                     patternDetectParameters.deviationSteps, 0); // 6; maximal iterations when looking for local maximum
		gd.addNumericField("High-pass filter for correlation with the model",patternDetectParameters.highpass, 3); //1.5 model correlation high-pass filter (relative to pattern fundamental frequency - average of 2)
		gd.addNumericField("Model correlation relative ring width",patternDetectParameters.corrRingWidth, 3);   //0.4, ring (around r=0.5 dist to opposite corr) width , center circle r=0.5*patternDetectParameters.corrRingWidth
		gd.addNumericField("Minimal pattern correlation contrast" ,patternDetectParameters.minCorrContrast, 3);   //5.0; // Discrimination threshold between good and bad pattern correleation
		gd.addCheckbox     ("Calculate correction for RED component",               colorComponents.colorsToCorrect[1]); //  true;
		gd.addCheckbox     ("Calculate correction for BLUE component",              colorComponents.colorsToCorrect[2]); //  true;
		gd.addCheckbox     ("Calculate correction for GREEN components (diagonal)", colorComponents.colorsToCorrect[4]); // false;
		gd.addCheckbox     ("Calculate correction for GREEN components (checker)",  colorComponents.colorsToCorrect[5]); //  true;
		gd.addCheckbox     ("Calculate correction for GREEN1 component (individual)",colorComponents.colorsToCorrect[0]); //  false;
		gd.addCheckbox     ("Calculate correction for GREEN2 component (individual)",colorComponents.colorsToCorrect[3]); //  false;
		gd.addNumericField ("Color component used as base for lateral chromatic (4,5)" ,colorComponents.referenceComponent, 0);   //4; // (change to 5 later) component to calculate lateral chromatic from (0 - G1, 1 - R, 2 - B, 3 - G2,4 - diagonal greens, 5 - checker greens)
		gd.addCheckbox     ("Equalize average values of the two greens in Bayer mosaic",colorComponents.equalizeGreens); //  true; Equalize average values of the two greens in Bayer mosaic
		gd.addNumericField("Decimate PSF before binning:",                           PSF_SUBPIXEL,    0); // OTF sub-pixel decimation
		gd.addNumericField("Reversed PSF  kernel size",                              inverseParameters.rSize, 0); //   32 - size of deconvolution kernel
		gd.addNumericField("Allowed overexposed pixels (fraction of the area) ",multiFilePSF.overexposedMaxFraction,3); //  0.005; // allowed fraction of the overexposed pixels in the PSF kernel measurement area 
		gd.addNumericField("Maximal number of concurrent threads",                   THREADS_MAX, 0); //   100
		gd.addNumericField("Debug Level:",      MASTER_DEBUG_LEVEL, 0);
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		JP4_INSTANCE.setTitle(gd.getNextString());
		JP4_INSTANCE.setURL(gd.getNextString());
		FFT_SIZE=4;
		for (i=(int) gd.getNextNumber(); i >4; i>>=1) FFT_SIZE <<=1; /* make FFT_SIZE to be power of 2*/
		MAP_FFT_SIZE=4;
		for (i=(int) gd.getNextNumber(); i >4; i>>=1) MAP_FFT_SIZE <<=1; /* make FFT_SIZE to be power of 2*/
		GAUSS_WIDTH=              gd.getNextNumber();
		FFT_OVERLAP=         (int) gd.getNextNumber();
		otfFilterParameters.deconvInvert=             gd.getNextNumber();  //0.05; // when FFT component is lass than this fraction of the maximal value, replace 1/z with Z
		patternDetectParameters.corrGamma=               gd.getNextNumber();  //0.2; lower the value - higher harmonics will participate in pattern frequency measurements
		patternDetectParameters.corrSigma=               gd.getNextNumber();  // 1.5; // pattern detection: high-pass filter (0.0 - none) gamma(PS)
		patternDetectParameters.diffSpectrCorr=     (int) gd.getNextNumber();
		patternDetectParameters.shrinkClusters=           gd.getNextNumber(); // 0.5; //  Shrink clusters by this ratio (remove lowest) after initial separation
		patternDetectParameters.multiplesToTry=     (int) gd.getNextNumber();
		patternDetectParameters.deviation=                gd.getNextNumber();
		patternDetectParameters.deviationSteps=    (int) gd.getNextNumber();
		patternDetectParameters.highpass=           gd.getNextNumber();
		patternDetectParameters.corrRingWidth=            gd.getNextNumber();
		patternDetectParameters.minCorrContrast=          gd.getNextNumber();
		colorComponents.colorsToCorrect[1]=          gd.getNextBoolean();
		colorComponents.colorsToCorrect[2]=         gd.getNextBoolean();
		colorComponents.colorsToCorrect[4]=     gd.getNextBoolean();
		colorComponents.colorsToCorrect[5]=      gd.getNextBoolean();
		colorComponents.colorsToCorrect[0]=       gd.getNextBoolean();
		colorComponents.colorsToCorrect[3]=       gd.getNextBoolean();
		colorComponents.referenceComponent=(int)   gd.getNextNumber();
		colorComponents.equalizeGreens=           gd.getNextBoolean();
		PSF_SUBPIXEL=            (int) gd.getNextNumber();
		inverseParameters.rSize=        (int) gd.getNextNumber();
		for (colorComponents.referenceComponent=5;(colorComponents.referenceComponent>=0) && (!colorComponents.colorsToCorrect[colorComponents.referenceComponent]); colorComponents.referenceComponent--);
		multiFilePSF.overexposedMaxFraction=  gd.getNextNumber();
		THREADS_MAX=             (int) gd.getNextNumber();
		MASTER_DEBUG_LEVEL=      (int) gd.getNextNumber();
		return true;
	}
/* ======================================================================== */
/* TODO: REPLACE doubleFHT  */
/* converts FHT results (frequency space) to complex numbers of [fftsize/2+1][fftsize] */

	private double[][][] FHT2FFTHalf (double [] fht_pixels, int fftsize) {
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


	private double[] FFTHalf2FHT (double [][][] fft, int fftsize) {
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



/* ======================================================================== */
/* ======================================================================== */
	private  double[] initWindowFunction(int size) {
		return initWindowFunction(size,GAUSS_WIDTH);
	}
	private  double[] initWindowFunction(int size, double gaussWidth) {
		double [] windowFunction =new double [size*size];
		double [] windowFunction_line=new double [size];
		double a,k;
		int i,j;
		if (gaussWidth<=0) {
			for (i=0; i<size; i++) windowFunction_line[i]= (0.54-0.46*Math.cos((i*2.0*Math.PI)/size));
		} else {
			k=2.0/(size*gaussWidth);
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
	private  double[][] normalizeAndWindow (double [][] pixels, double [] windowFunction) {
		return normalizeAndWindow (pixels, windowFunction, true);
	}
	private  double[] normalizeAndWindow (double [] pixels, double [] windowFunction) {
		return normalizeAndWindow (pixels, windowFunction, true);
	}
	private  double[][] normalizeAndWindow (double [][] pixels, double [] windowFunction, boolean removeDC) {
		int i;
		for (i=0;i<pixels.length;i++)  if (pixels[i]!=null) pixels[i]=normalizeAndWindow (pixels[i],  windowFunction, removeDC);
		return pixels;
	}
	private  double[] normalizeAndWindow (double [] pixels, double [] windowFunction, boolean removeDC) {
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
//	copied from IJ.java 
	private String updateExtension(String path, String extension) {
		if (path==null) return null;
		int dotIndex = path.lastIndexOf(".");
		int separatorIndex = path.lastIndexOf(File.separator);
		if (dotIndex>=0 && dotIndex>separatorIndex && (path.length()-dotIndex)<=5) {
			if (dotIndex+1<path.length() && Character.isDigit(path.charAt(dotIndex+1)))
				path += extension;
			else
				path = path.substring(0, dotIndex) + extension;
		} else
			path += extension;
		return path;
	}
	
/* ======================================================================== */
/* === Parameter classes === */
	/*
	public static class MultiFilePSF {
		public double  overexposedMaxFraction; // allowed fraction of the overexposed pixels in the PSF kernel measurement area 
		public double  validateThreshold;      // fraction of full PSF "energy"
		public boolean validateShowEllipse;    // show ellipse parameters of partial PSF arrays
		public boolean showWeights;            // show image indicating frame coverage
		public boolean fillMissing;            // replace missing kernels with neighbors
		public MultiFilePSF (
				double  overexposedMaxFraction,
				double  validateThreshold,
				boolean validateShowEllipse,
				boolean showWeights,
				boolean fillMissing
		) {
			this.overexposedMaxFraction=overexposedMaxFraction;
			this.validateThreshold=validateThreshold;
			this.validateShowEllipse=validateShowEllipse;
			this.showWeights=showWeights;
			this.fillMissing=fillMissing;
		}

		public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"overexposedMaxFraction",this.overexposedMaxFraction+"");
			properties.setProperty(prefix+"validateThreshold",this.validateThreshold+"");
			properties.setProperty(prefix+"validateShowEllipse",this.validateShowEllipse+"");
			properties.setProperty(prefix+"showWeights",this.showWeights+"");
			properties.setProperty(prefix+"fillMissing",this.fillMissing+"");
		}
		public void getProperties(String prefix,Properties properties){
			if (properties.getProperty(prefix+"overexposedMaxFraction")!=null) this.overexposedMaxFraction=Double.parseDouble(properties.getProperty(prefix+"overexposedMaxFraction"));
			if (properties.getProperty(prefix+"validateThreshold")!=null)this.validateThreshold=Double.parseDouble(properties.getProperty(prefix+"validateThreshold"));
			if (properties.getProperty(prefix+"validateShowEllipse")!=null)this.validateShowEllipse=Boolean.parseBoolean(properties.getProperty(prefix+"validateShowEllipse"));
			if (properties.getProperty(prefix+"showWeights")!=null)this.showWeights=Boolean.parseBoolean(properties.getProperty(prefix+"showWeights"));
			if (properties.getProperty(prefix+"fillMissing")!=null)this.fillMissing=Boolean.parseBoolean(properties.getProperty(prefix+"fillMissing"));
			
		}

	}

	public static class ColorComponents {
		public boolean [] colorsToCorrect=    new boolean[6];
		public int        referenceComponent; // component to calculate lateral chromatic from (0 - G1, 1 - R, 2 - B, 3 - G2,4 - diagonal greens, 5 - checker greens)
		public boolean    equalizeGreens;   // equalize 2 greens in Bayer mosaic
		public ColorComponents (
				boolean green1,
				boolean red,
				boolean blue,
				boolean green2,
				boolean diagonal, // both greens combined in a 45-degree rotated array 
				boolean checker,   // both greens combined in a checkerboard pattern
				int        referenceComponent,
				boolean    equalizeGreens
		) {
			this.colorsToCorrect[0]=green1;
			this.colorsToCorrect[1]=red;
			this.colorsToCorrect[2]=blue;
			this.colorsToCorrect[3]=green2;
			this.colorsToCorrect[4]=diagonal;
			this.colorsToCorrect[5]=checker;
			this.referenceComponent=referenceComponent;
			this.equalizeGreens=equalizeGreens;
		}
		public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"green1",this.colorsToCorrect[0]+"");
			properties.setProperty(prefix+"red",this.colorsToCorrect[1]+"");
			properties.setProperty(prefix+"blue",this.colorsToCorrect[2]+"");
			properties.setProperty(prefix+"green2",this.colorsToCorrect[3]+"");
			properties.setProperty(prefix+"diagonal",this.colorsToCorrect[4]+"");
			properties.setProperty(prefix+"checker",this.colorsToCorrect[5]+"");
			properties.setProperty(prefix+"referenceComponent",this.referenceComponent+"");
			properties.setProperty(prefix+"equalizeGreens",this.equalizeGreens+"");
		}
		public void getProperties(String prefix,Properties properties){
			this.colorsToCorrect[0]=Boolean.parseBoolean(properties.getProperty(prefix+"green1"));
			this.colorsToCorrect[1]=Boolean.parseBoolean(properties.getProperty(prefix+"red"));
			this.colorsToCorrect[2]=Boolean.parseBoolean(properties.getProperty(prefix+"blue"));
			this.colorsToCorrect[3]=Boolean.parseBoolean(properties.getProperty(prefix+"green2"));
			this.colorsToCorrect[4]=Boolean.parseBoolean(properties.getProperty(prefix+"diagonal"));
			this.colorsToCorrect[5]=Boolean.parseBoolean(properties.getProperty(prefix+"checker"));
			this.referenceComponent=Integer.parseInt(properties.getProperty(prefix+"referenceComponent"));
			this.equalizeGreens=Boolean.parseBoolean(properties.getProperty(prefix+"equalizeGreens"));
		}
	}
*/
	public static class ShowResults {
		public boolean       showPSF; // show combined PSF kernels (per-color and/or composite, as defined in SHOW_INDIVIDUAL, SHOW_COMPOSITE)
		public boolean       showMTF; // calculate/show MTF (see notes to SHOW_RESULTS.showPSF)
		public boolean  showInverted; // show inverted kernels (unfiltered), same notes
		public boolean  showFiltered; // filter and show inverted kernels
		public boolean showGaussians; // create Gaussian kernels with the same centers as inverted ones (low noise, use in low details areas)
		public ShowResults (
				boolean showPSF,
				boolean showMTF,
				boolean showInverted,
				boolean showFiltered,
				boolean showGaussians
		){
			this.showPSF=showPSF;
			this.showMTF=showMTF;
			this.showInverted=showInverted;
			this.showFiltered=showFiltered;
			this.showGaussians=showGaussians;
		}

		public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"showPSF",this.showPSF+"");
			properties.setProperty(prefix+"showMTF",this.showMTF+"");
			properties.setProperty(prefix+"showInverted",this.showInverted+"");
			properties.setProperty(prefix+"showFiltered",this.showFiltered+"");
			properties.setProperty(prefix+"showGaussians",this.showGaussians+"");
		}
		public void getProperties(String prefix,Properties properties){
			this.showPSF=Boolean.parseBoolean(properties.getProperty(prefix+"showPSF"));
			this.showMTF=Boolean.parseBoolean(properties.getProperty(prefix+"showMTF"));
			this.showInverted=Boolean.parseBoolean(properties.getProperty(prefix+"showInverted"));
			this.showFiltered=Boolean.parseBoolean(properties.getProperty(prefix+"showFiltered"));
			this.showGaussians=Boolean.parseBoolean(properties.getProperty(prefix+"showGaussians"));
		}
		
	}
/*
	public static class OTFFilterParameters {
		public double deconvInvert;
		public double zerofreqSize;
		public double smoothPS;
		public double thresholdHigh;
		public double thresholdLow;

		public OTFFilterParameters(
				double deconvInvert,
				double zerofreqSize,
				double smoothPS,
				double thresholdHigh,
				double thresholdLow) {
			this.deconvInvert = deconvInvert;
			this.zerofreqSize = zerofreqSize;
			this.smoothPS = smoothPS;
			this.thresholdHigh = thresholdHigh;
			this.thresholdLow = thresholdLow;
		}
        public OTFFilterParameters clone() {
        	return new OTFFilterParameters(
        			this.deconvInvert,
        			this.zerofreqSize,
        			this.smoothPS,
        			this.thresholdHigh,
        			this.thresholdLow);
        }
		public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"deconvInvert",this.deconvInvert+"");
			properties.setProperty(prefix+"zerofreqSize",this.zerofreqSize+"");
			properties.setProperty(prefix+"smoothPS",this.smoothPS+"");
			properties.setProperty(prefix+"thresholdHigh",this.thresholdHigh+"");
			properties.setProperty(prefix+"thresholdLow",this.thresholdLow+"");
		}
		public void getProperties(String prefix,Properties properties){
			this.deconvInvert=Double.parseDouble(properties.getProperty(prefix+"deconvInvert"));
			this.zerofreqSize=Double.parseDouble(properties.getProperty(prefix+"zerofreqSize"));
			this.smoothPS=Double.parseDouble(properties.getProperty(prefix+"smoothPS"));
			this.thresholdHigh=Double.parseDouble(properties.getProperty(prefix+"thresholdHigh"));
			this.thresholdLow=Double.parseDouble(properties.getProperty(prefix+"thresholdLow"));
		}

	}


	public static class PSFParameters {
		public double minContrast;
		public double windowFrac;
		public boolean useWindow;
		public boolean symm180;
		public boolean ignoreChromatic;
		public double smoothSeparate;
		public double topCenter;
		public double sigmaToRadius;
		public double wingsEnergy;
		public double wingsEllipseScale;
		public double minDefinedArea;   // minimal (weighted) fraction of the defined patter pixels in the FFT area
		public boolean approximateGrid; // approximate grid with polynomial 
		public boolean centerPSF;       // Center PSF by modifying phase 
		public double mask1_sigma;
		public double mask1_threshold;
		public double gaps_sigma;
		public double mask_denoise;
		

		public PSFParameters(double minContrast,
				double windowFrac,
				boolean useWindow,
				boolean symm180,
				boolean ignoreChromatic,
				double smoothSeparate,
				double topCenter,
				double sigmaToRadius,
				double wingsEnergy,
				double wingsEllipseScale,
				double minDefinedArea, // minimal (weighted) fraction of the defined patter pixels in the FFT area
				boolean approximateGrid, // approximate grid with polynomial
				boolean centerPSF,       // Center PSF by modifying phase
				double mask1_sigma,
				double mask1_threshold,
				double gaps_sigma,
				double mask_denoise

		) {
			this.minContrast = minContrast;
			this.windowFrac = windowFrac;
			this.useWindow = useWindow;
			this.symm180 = symm180;
			this.ignoreChromatic = ignoreChromatic;
			this.smoothSeparate = smoothSeparate;
			this.topCenter = topCenter;
			this.sigmaToRadius = sigmaToRadius;
			this.wingsEnergy = wingsEnergy;
			this.wingsEllipseScale = wingsEllipseScale;
			this.minDefinedArea = minDefinedArea; // minimal (weighted) fraction of the defined patter pixels in the FFT area
			this.approximateGrid = approximateGrid; // approximate grid with polynomial 
			this.centerPSF = centerPSF; // approximate grid with polynomial 
			this.mask1_sigma = mask1_sigma;
			this.mask1_threshold = mask1_threshold;
			this.gaps_sigma=gaps_sigma;
			this.mask_denoise=mask_denoise;


		}
        public PSFParameters clone(){
        	return new PSFParameters(
        			this.minContrast,
        			this.windowFrac,
        			this.useWindow,
        			this.symm180,
        			this.ignoreChromatic,
        			this.smoothSeparate,
        			this.topCenter,
        			this.sigmaToRadius,
        			this.wingsEnergy,
        			this.wingsEllipseScale,
        			this.minDefinedArea, // minimal (weighted) fraction of the defined patter pixels in the FFT area
        			this.approximateGrid, // approximate grid with polynomial
        			this.centerPSF, // approximate grid with polynomial
        			this.mask1_sigma,
        			this.mask1_threshold,
        			this.gaps_sigma,
        			this.mask_denoise
                      );
        }
		public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"minContrast",this.minContrast+"");
			properties.setProperty(prefix+"windowFrac",this.windowFrac+"");
			properties.setProperty(prefix+"useWindow",this.useWindow+"");
			properties.setProperty(prefix+"symm180",this.symm180+"");
			properties.setProperty(prefix+"ignoreChromatic",this.ignoreChromatic+"");
			properties.setProperty(prefix+"smoothSeparate",this.smoothSeparate+"");
			properties.setProperty(prefix+"topCenter",this.topCenter+"");
			properties.setProperty(prefix+"sigmaToRadius",this.sigmaToRadius+"");
			properties.setProperty(prefix+"wingsEnergy",this.wingsEnergy+"");
			properties.setProperty(prefix+"wingsEllipseScale",this.wingsEllipseScale+"");
			properties.setProperty(prefix+"minDefinedArea",this.minDefinedArea+"");
			properties.setProperty(prefix+"approximateGrid",this.approximateGrid+"");
			properties.setProperty(prefix+"centerPSF",this.centerPSF+"");
			properties.setProperty(prefix+"mask1_sigma",this.mask1_sigma+"");
			properties.setProperty(prefix+"mask1_threshold",this.mask1_threshold+"");
			properties.setProperty(prefix+"gaps_sigma",this.gaps_sigma+"");
			properties.setProperty(prefix+"mask_denoise",this.mask_denoise+"");
		}
		public void getProperties(String prefix,Properties properties){
			if (properties.getProperty(prefix+"minContrast")!=null)       this.minContrast=Double.parseDouble(properties.getProperty(prefix+"minContrast"));
			if (properties.getProperty(prefix+"windowFrac")!=null)        this.windowFrac=Double.parseDouble(properties.getProperty(prefix+"windowFrac"));
			if (properties.getProperty(prefix+"useWindow")!=null)         this.useWindow=Boolean.parseBoolean(properties.getProperty(prefix+"useWindow"));
			if (properties.getProperty(prefix+"symm180")!=null)           this.symm180=Boolean.parseBoolean(properties.getProperty(prefix+"symm180"));
			if (properties.getProperty(prefix+"ignoreChromatic")!=null)   this.ignoreChromatic=Boolean.parseBoolean(properties.getProperty(prefix+"ignoreChromatic"));
			if (properties.getProperty(prefix+"smoothSeparate")!=null)    this.smoothSeparate=Double.parseDouble(properties.getProperty(prefix+"smoothSeparate"));
			if (properties.getProperty(prefix+"topCenter")!=null)         this.topCenter=Double.parseDouble(properties.getProperty(prefix+"topCenter"));
			if (properties.getProperty(prefix+"sigmaToRadius")!=null)     this.sigmaToRadius=Double.parseDouble(properties.getProperty(prefix+"sigmaToRadius"));
			if (properties.getProperty(prefix+"wingsEnergy")!=null)       this.wingsEnergy=Double.parseDouble(properties.getProperty(prefix+"wingsEnergy"));
			if (properties.getProperty(prefix+"wingsEllipseScale")!=null) this.wingsEllipseScale=Double.parseDouble(properties.getProperty(prefix+"wingsEllipseScale"));
			if (properties.getProperty(prefix+"minDefinedArea")!=null)    this.minDefinedArea=Double.parseDouble(properties.getProperty(prefix+"minDefinedArea"));
			if (properties.getProperty(prefix+"approximateGrid")!=null)   this.approximateGrid=Boolean.parseBoolean(properties.getProperty(prefix+"approximateGrid"));
			if (properties.getProperty(prefix+"centerPSF")!=null)         this.centerPSF=Boolean.parseBoolean(properties.getProperty(prefix+"centerPSF"));
			if (properties.getProperty(prefix+"mask1_sigma")!=null)       this.mask1_sigma=Double.parseDouble(properties.getProperty(prefix+"mask1_sigma"));
			if (properties.getProperty(prefix+"mask1_threshold")!=null)   this.mask1_threshold=Double.parseDouble(properties.getProperty(prefix+"mask1_threshold"));
			if (properties.getProperty(prefix+"gaps_sigma")!=null)        this.mask1_threshold=Double.parseDouble(properties.getProperty(prefix+"gaps_sigma"));
			if (properties.getProperty(prefix+"mask_denoise")!=null)        this.mask_denoise=Double.parseDouble(properties.getProperty(prefix+"mask_denoise"));
		}
	}

	public static class InverseParameters {
		public int dSize;
		public int rSize;
		public double invertRange;
		public double otfCutoffEnergy;
		public double otfEllipseScale;
		public boolean otfEllipseGauss;
		public double psfCutoffEnergy;
		public double psfEllipseScale;
		public double rpsfMinMaskThreshold;
		public boolean filter;
		public double blurIndividual;
		public double blurDiagonal;
		public double blurChecker;
		public double gaussianSigmaIndividual;
		public double gaussianSigmaDiagonal;
		public double gaussianSigmaChecker;
		public double sigmaScale;
		public double sigmaToRadius;
		public boolean filterDirect;
		public double sigmaScaleDirect;
		public double sigmaToRadiusDirect;

		public InverseParameters(int dSize, int rSize, double invertRange,
				double otfCutoffEnergy, double otfEllipseScale,
				boolean otfEllipseGauss, double psfCutoffEnergy,
				double psfEllipseScale, double rpsfMinMaskThreshold,
				boolean filter, double blurIndividual, double blurDiagonal, double blurChecker,
				double gaussianSigmaIndividual, double gaussianSigmaDiagonal, double gaussianSigmaChecker,
				double sigmaScale, double sigmaToRadius,
				boolean filterDirect, double sigmaScaleDirect, double sigmaToRadiusDirect
				) {
			this.dSize = dSize;
			this.rSize = rSize;
			this.invertRange = invertRange;
			this.otfCutoffEnergy = otfCutoffEnergy;
			this.otfEllipseScale = otfEllipseScale;
			this.otfEllipseGauss = otfEllipseGauss;
			this.psfCutoffEnergy = psfCutoffEnergy;
			this.psfEllipseScale = psfEllipseScale;
			this.rpsfMinMaskThreshold = rpsfMinMaskThreshold;
			this.filter = filter;
			this.blurIndividual = blurIndividual;
			this.blurDiagonal = blurDiagonal;
			this.blurChecker = blurChecker;
			this.gaussianSigmaIndividual = gaussianSigmaIndividual;
			this.gaussianSigmaDiagonal = gaussianSigmaDiagonal;
			this.gaussianSigmaChecker = gaussianSigmaChecker;
			this.sigmaScale = sigmaScale;
			this.sigmaToRadius = sigmaToRadius;
			this.filterDirect=filterDirect;
			this.sigmaScaleDirect=sigmaScaleDirect;
			this.sigmaToRadiusDirect=sigmaToRadiusDirect;
		}

		public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"dSize",this.dSize+"");
			properties.setProperty(prefix+"rSize",this.rSize+"");
			properties.setProperty(prefix+"invertRange",this.invertRange+"");
			properties.setProperty(prefix+"otfCutoffEnergy",this.otfCutoffEnergy+"");
			properties.setProperty(prefix+"otfEllipseScale",this.otfEllipseScale+"");
			properties.setProperty(prefix+"otfEllipseGauss",this.otfEllipseGauss+"");
			properties.setProperty(prefix+"psfCutoffEnergy",this.psfCutoffEnergy+"");
			properties.setProperty(prefix+"psfEllipseScale",this.psfEllipseScale+"");
			properties.setProperty(prefix+"rpsfMinMaskThreshold",this.rpsfMinMaskThreshold+"");
			properties.setProperty(prefix+"filter",this.filter+"");
			properties.setProperty(prefix+"blurIndividual",this.blurIndividual+"");
			properties.setProperty(prefix+"blurDiagonal",this.blurDiagonal+"");
			properties.setProperty(prefix+"blurChecker",this.blurChecker+"");
			properties.setProperty(prefix+"gaussianSigmaIndividual",this.gaussianSigmaIndividual+"");
			properties.setProperty(prefix+"gaussianSigmaDiagonal",this.gaussianSigmaDiagonal+"");
			properties.setProperty(prefix+"gaussianSigmaChecker",this.gaussianSigmaChecker+"");
			properties.setProperty(prefix+"sigmaScale",this.sigmaScale+"");
			properties.setProperty(prefix+"sigmaToRadius",this.sigmaToRadius+"");
			properties.setProperty(prefix+"filterDirect",this.filterDirect+"");
			properties.setProperty(prefix+"sigmaScaleDirect",this.sigmaScaleDirect+"");
			properties.setProperty(prefix+"sigmaToRadiusDirect",this.sigmaToRadiusDirect+"");
		}
		public void getProperties(String prefix,Properties properties){
			this.dSize=Integer.parseInt(properties.getProperty(prefix+"dSize"));
			this.rSize=Integer.parseInt(properties.getProperty(prefix+"rSize"));
			this.invertRange=Double.parseDouble(properties.getProperty(prefix+"invertRange"));
			this.otfCutoffEnergy=Double.parseDouble(properties.getProperty(prefix+"otfCutoffEnergy"));
			this.otfEllipseScale=Double.parseDouble(properties.getProperty(prefix+"otfEllipseScale"));
			this.otfEllipseGauss=Boolean.parseBoolean(properties.getProperty(prefix+"otfEllipseGauss"));
			this.psfCutoffEnergy=Double.parseDouble(properties.getProperty(prefix+"psfCutoffEnergy"));
			this.psfEllipseScale=Double.parseDouble(properties.getProperty(prefix+"psfEllipseScale"));
			this.rpsfMinMaskThreshold=Double.parseDouble(properties.getProperty(prefix+"rpsfMinMaskThreshold"));
			this.filter=Boolean.parseBoolean(properties.getProperty(prefix+"filter"));
			this.blurIndividual=Double.parseDouble(properties.getProperty(prefix+"blurIndividual"));
			this.blurDiagonal=Double.parseDouble(properties.getProperty(prefix+"blurDiagonal"));
			this.blurChecker=Double.parseDouble(properties.getProperty(prefix+"blurChecker"));
			this.gaussianSigmaIndividual=Double.parseDouble(properties.getProperty(prefix+"gaussianSigmaIndividual"));
			this.gaussianSigmaDiagonal=Double.parseDouble(properties.getProperty(prefix+"gaussianSigmaDiagonal"));
			this.gaussianSigmaChecker=Double.parseDouble(properties.getProperty(prefix+"gaussianSigmaChecker"));
			this.sigmaScale=Double.parseDouble(properties.getProperty(prefix+"sigmaScale"));
			this.sigmaToRadius=Double.parseDouble(properties.getProperty(prefix+"sigmaToRadius"));
			this.filterDirect=Boolean.parseBoolean(properties.getProperty(prefix+"filterDirect"));
			this.sigmaScaleDirect=Double.parseDouble(properties.getProperty(prefix+"sigmaScaleDirect"));
			this.sigmaToRadiusDirect=Double.parseDouble(properties.getProperty(prefix+"sigmaToRadiusDirect"));
		}
	}
	
	
	public static class InterpolateParameters {
		public int    size;        // size of each kernel (should be square)
		public int    step;        // number of subdivisions from input to output
		public int    add_top;     // add this number of kernel rows to the output above the existent/interpolated
		public int    add_left;    // add this number of kernel columns to the output on the left of the existent/interpolated
		public int    add_right;   // add this number of kernel columns to the output on the right of the existent/interpolated
		public int    add_bottom;  // add this number of kernel rows to the output below the existent/interpolated
		public double extrapolate; // 0 - duplicate, 1.0 - extrapolate outside of the known kernels

		public InterpolateParameters(
				int    size,
				int    step,
				int    add_top,
				int    add_left,
				int    add_right,
				int    add_bottom,
				double extrapolate
		) {
			this.size=size;
			this.step=step;
			this.add_top=add_top;
			this.add_left=add_left;
			this.add_right=add_right;
			this.add_bottom=add_bottom;
			this.extrapolate=extrapolate;
		}
		public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"size",this.size+"");
			properties.setProperty(prefix+"step",this.step+"");
			properties.setProperty(prefix+"add_top",this.add_top+"");
			properties.setProperty(prefix+"add_left",this.add_left+"");
			properties.setProperty(prefix+"add_right",this.add_right+"");
			properties.setProperty(prefix+"add_bottom",this.add_bottom+"");
			properties.setProperty(prefix+"extrapolate",this.extrapolate+"");
		}
		public void getProperties(String prefix,Properties properties){
			this.size=Integer.parseInt(properties.getProperty(prefix+"size"));
			this.step=Integer.parseInt(properties.getProperty(prefix+"step"));
			this.add_top=Integer.parseInt(properties.getProperty(prefix+"add_top"));
			this.add_left=Integer.parseInt(properties.getProperty(prefix+"add_left"));
			this.add_right=Integer.parseInt(properties.getProperty(prefix+"add_right"));
			this.add_bottom=Integer.parseInt(properties.getProperty(prefix+"add_bottom"));
			this.extrapolate=Double.parseDouble(properties.getProperty(prefix+"extrapolate"));
		}
		
	}
*/	
    public static class FlatFieldParameters {
	    public int numEyesisChannels=3;
	    public int numEyesisSubChannels=3;

	    public String []  sourceDirPaths=null;
    	public boolean [] sourceDirPathsEn=null;
    	public boolean    normalize=true; // normalize each color component before averaging
///			0.025, // overexposedMaxFraction -  allowed fraction of the overexposed pixels in the PSF kernel measurement area 
// to filter good/bad images in each sub-band    	
    	public double overExpValue= 0.99;
    	public double overExpFrac=  0.025;
    	public double underExpValue=0.25;
    	public double underExpFrac= 0.5;
    	
    	public boolean    noTiltEdges=true; // do not apply tilt to corners (to get closer)
    	public int        functionType=0; //0 polynomial, 1 - power
    	public int        functionModifier=0; //additional modifier for the approxiamtion function
    	public double     section34=0.5; // location of 4-th and 5-th section (ratio from o to 1 and from 0 to 2 (1-3-0-4-2) 
///Levenberg–Marquardt parameters
    	public double     centerWeight=0.5; // weight for the error function will be proportional to r^2 (r - half smallest dimension) plus weight in the center
    	public boolean    LM_auto=true;   // automatically iterate (false open - dialogs)
    	public double     LM_lambdaInitial=0.001;
    	public double     LM_lambdaStepUp=   8.0; // multiply lambda by this if result is worse
    	public double     LM_lambdaStepDown= 0.5; // multiply lambda by this if result is better
    	public double     LM_thresholdFinish=0.0001; // stop iterations if 2 last steps had less improvement (but not worsening ) 
    	public int        LM_numIterations=  100; // maximal number of iterations 
    	
    	
    	public double     fatZero=3; // normalize each color component before averaging
    	
    	
    	public int []     margins= {100,100,100,100};
    	public int        decimate=2;
    	public int        sampleWidth=16;
    	public double     highPassSigma=512.0;
    	public double     maxTilt=0.25; // real life 0.16
		public String     flatFieldDirectory;            // results with flat field calibration files
		public boolean    eyesisMode;                   // multiple channels
		public String     sourceFileExtension="jp4";
		public boolean    processAllChannels;            // if true - process all channels, otherwise only enabled in processChannels[]       
		public boolean [] processChannels=new boolean [9];
		public boolean    useXML;                        // save/restore settings as xml file
		public boolean    saveSettings;                  // save current settings in results directory
	    public  FlatFieldParameters (
	    		    String  sourceFileExtension,
	    		    boolean normalize,
	    	    	double overExpValue,
	    	    	double overExpFrac,
	    	    	double underExpValue,
	    	    	double underExpFrac,
	    	    	boolean    noTiltEdges, // do not apply tilt to corners (to get closer)
	    		    int     functionType,
	    		    int        functionModifier, //additional modifier for the approxiamtion function
	    	    	double     section34, // location of 4-th and 5-th section (ratio from o to 1 and from 0 to 2 (1-3-0-4-2) 
	    	    	double     centerWeight, // weight for the error function will be proportional to r^2 (r - half smallest dimension) plus weight in the center
	    		  ///Levenberg–Marquardt parameters
	    	    	boolean    LM_auto,   // automatically iterate (false open - dialogs)
	    	    	double     LM_lambdaInitial,
	    	    	double     LM_lambdaStepUp,
	    	    	double     LM_lambdaStepDown,
	    	    	double     LM_thresholdFinish, 
	    	    	int        LM_numIterations, 
	    		    double  fatZero,
	    		    int     margin_left,
	    		    int     margin_right,
	    		    int     margin_top,
	    		    int     margin_bottom,
	    		    int     decimate,
	    		    int     sampleWidth,
	    		    double  highPassSigma,
	    	    	double  maxTilt,       // real life 0.16
					String  flatFieldDirectory,           // results with flat field calibration files
					boolean    eyesisMode,                   // multiple channels
					boolean processAllChannels,           // if true - process all channels, otherwise only enabled in processChannels[]       
					boolean processChannels11,
					boolean processChannels12,
					boolean processChannels13,
					boolean processChannels21,
					boolean processChannels22,
					boolean processChannels23,
					boolean processChannels31,
					boolean processChannels32,
					boolean processChannels33,
	  		        boolean useXML,
					boolean saveSettings
			) {
	    	this.sourceFileExtension=  sourceFileExtension;
			this.normalize=            normalize;
			this.overExpValue=         overExpValue;
			this.overExpFrac=          overExpFrac;
			this.underExpValue=        underExpValue;
			this.underExpFrac=         underExpFrac;
			this.noTiltEdges=          noTiltEdges;
			this.functionType=         functionType;
			this.functionModifier=     functionModifier;
			this.section34=            section34;
			this.centerWeight=         centerWeight;
			this.LM_auto=              LM_auto;
			this.LM_lambdaInitial=     LM_lambdaInitial;
			this.LM_lambdaStepUp=      LM_lambdaStepUp;
			this.LM_lambdaStepDown=    LM_lambdaStepDown;
			this.LM_thresholdFinish=   LM_thresholdFinish;
			this.LM_numIterations=     LM_numIterations;
			this.fatZero=              fatZero;
			this.margins[0]=           margin_left;
			this.margins[1]=           margin_right;
			this.margins[2]=           margin_top;
			this.margins[3]=           margin_bottom;
			this.decimate=             decimate;
			this.sampleWidth=          sampleWidth;
			this.highPassSigma=        highPassSigma;
			this.maxTilt=              maxTilt;
			this.flatFieldDirectory=   flatFieldDirectory;
			this.eyesisMode=           eyesisMode;
			this.processAllChannels=   processAllChannels;
			this.processChannels[0]=   processChannels11;
			this.processChannels[1]=   processChannels12;
			this.processChannels[2]=   processChannels13;
			this.processChannels[3]=   processChannels21;
			this.processChannels[4]=   processChannels22;
			this.processChannels[5]=   processChannels23;
			this.processChannels[6]=   processChannels31;
			this.processChannels[7]=   processChannels32;
			this.processChannels[8]=   processChannels33;
			this.useXML=useXML;
			this.saveSettings=saveSettings;
	    }

		public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"sourceFileExtension",this.sourceFileExtension);
			properties.setProperty(prefix+"numEyesisChannels",this.numEyesisChannels+"");
			properties.setProperty(prefix+"numEyesisSubChannels",this.numEyesisSubChannels+"");
			properties.setProperty(prefix+"normalize",this.normalize+"");
			properties.setProperty(prefix+"overExpValue",this.overExpValue+"");
			properties.setProperty(prefix+"overExpFrac",this.overExpFrac+"");
			properties.setProperty(prefix+"underExpValue",this.underExpValue+"");
			properties.setProperty(prefix+"underExpFrac",this.underExpFrac+"");
			properties.setProperty(prefix+"noTiltEdges",this.noTiltEdges+"");
			
			
			properties.setProperty(prefix+"functionType",this.functionType+"");
			properties.setProperty(prefix+"functionModifier",this.functionModifier+"");
			properties.setProperty(prefix+"section34",this.section34+"");
			
			properties.setProperty(prefix+"centerWeight",this.centerWeight+"");
			properties.setProperty(prefix+"LM_auto",this.LM_auto+"");
			properties.setProperty(prefix+"LM_lambdaInitial",this.LM_lambdaInitial+"");
			properties.setProperty(prefix+"LM_lambdaStepUp",this.LM_lambdaStepUp+"");
			properties.setProperty(prefix+"LM_lambdaStepDown",this.LM_lambdaStepDown+"");
			properties.setProperty(prefix+"LM_thresholdFinish",this.LM_thresholdFinish+"");
			properties.setProperty(prefix+"LM_numIterations",this.LM_numIterations+"");
			
			properties.setProperty(prefix+"fatZero",this.fatZero+"");
			properties.setProperty(prefix+"margin_left",this.margins[0]+"");
			properties.setProperty(prefix+"margin_right",this.margins[1]+"");
			properties.setProperty(prefix+"margin_top",this.margins[2]+"");
			properties.setProperty(prefix+"margin_bottom",this.margins[3]+"");
			properties.setProperty(prefix+"decimate",this.decimate+"");
			properties.setProperty(prefix+"sampleWidth",this.sampleWidth+"");
			properties.setProperty(prefix+"highPassSigma",this.highPassSigma+"");
			properties.setProperty(prefix+"maxTilt",this.maxTilt+"");
			
			properties.setProperty(prefix+"flatFieldDirectory",this.flatFieldDirectory);
			properties.setProperty(prefix+"eyesisMode",        this.eyesisMode+"");
			properties.setProperty(prefix+"processAllChannels",this.processAllChannels+"");
			properties.setProperty(prefix+"processChannels11",this.processChannels[0]+"");
			properties.setProperty(prefix+"processChannels12",this.processChannels[1]+"");
			properties.setProperty(prefix+"processChannels13",this.processChannels[2]+"");
			properties.setProperty(prefix+"processChannels21",this.processChannels[3]+"");
			properties.setProperty(prefix+"processChannels22",this.processChannels[4]+"");
			properties.setProperty(prefix+"processChannels23",this.processChannels[5]+"");
			properties.setProperty(prefix+"processChannels31",this.processChannels[6]+"");
			properties.setProperty(prefix+"processChannels32",this.processChannels[7]+"");
			properties.setProperty(prefix+"processChannels33",this.processChannels[8]+"");
			properties.setProperty(prefix+"useXML",this.useXML+"");
			properties.setProperty(prefix+"saveSettings",this.saveSettings+"");
			int j=(this.sourceDirPaths==null)?0:this.sourceDirPaths.length;
			properties.setProperty(prefix+"sourceDirPaths_length",j+"");
			for (int i=0;i<j;i++) {
				properties.setProperty(prefix+"sourceDirPaths_"+i,this.sourceDirPaths[i]);
				properties.setProperty(prefix+"sourceDirPathsEn_"+i,this.sourceDirPathsEn[i]+"");
			}	
		}
		
		public void getProperties(String prefix,Properties properties){
			if (properties.getProperty(prefix+"sourceFileExtension")!=null) this.sourceFileExtension=   properties.getProperty(prefix+"sourceFileExtension");
			if (properties.getProperty(prefix+"numEyesisChannels")!=null) this.numEyesisChannels=Integer.parseInt(properties.getProperty(prefix+"numEyesisChannels"));
			if (properties.getProperty(prefix+"numEyesisSubChannels")!=null) this.numEyesisSubChannels=Integer.parseInt(properties.getProperty(prefix+"numEyesisSubChannels"));
			if (properties.getProperty(prefix+"normalize")!=null)          this.normalize=          Boolean.parseBoolean(properties.getProperty(prefix+"normalize"));

			if (properties.getProperty(prefix+"overExpValue")!=  null)    this.overExpValue=        Double.parseDouble(properties.getProperty(prefix+"overExpValue"));
			if (properties.getProperty(prefix+"overExpFrac")!=  null)    this.overExpFrac=        Double.parseDouble(properties.getProperty(prefix+"overExpFrac"));
			if (properties.getProperty(prefix+"underExpValue")!=  null)    this.underExpValue=        Double.parseDouble(properties.getProperty(prefix+"underExpValue"));
			if (properties.getProperty(prefix+"underExpFrac")!=  null)    this.underExpFrac=        Double.parseDouble(properties.getProperty(prefix+"underExpFrac"));

			if (properties.getProperty(prefix+"functionType")!=  null)     this.functionType=       Integer.parseInt(properties.getProperty(prefix+"functionType"));
			if (properties.getProperty(prefix+"functionModifier")!=  null) this.functionModifier=   Integer.parseInt(properties.getProperty(prefix+"functionModifier"));
			
			if (properties.getProperty(prefix+"section34")!=null)          this.section34=          Double.parseDouble(properties.getProperty(prefix+"section34"));
			
			if (properties.getProperty(prefix+"noTiltEdges")!=null)          this.noTiltEdges=          Boolean.parseBoolean(properties.getProperty(prefix+"noTiltEdges"));
			if (properties.getProperty(prefix+"centerWeight")!=null)       this.centerWeight=       Double.parseDouble(properties.getProperty(prefix+"centerWeight"));
			if (properties.getProperty(prefix+"LM_auto")!=null)            this.LM_auto=            Boolean.parseBoolean(properties.getProperty(prefix+"LM_auto"));
			if (properties.getProperty(prefix+"LM_lambdaInitial")!=null)   this.LM_lambdaInitial=   Double.parseDouble(properties.getProperty(prefix+"LM_lambdaInitial"));
			if (properties.getProperty(prefix+"LM_lambdaStepUp")!=null)    this.LM_lambdaStepUp=    Double.parseDouble(properties.getProperty(prefix+"LM_lambdaStepUp"));
			if (properties.getProperty(prefix+"LM_lambdaStepDown")!=null)  this.LM_lambdaStepDown=  Double.parseDouble(properties.getProperty(prefix+"LM_lambdaStepDown"));
			if (properties.getProperty(prefix+"LM_thresholdFinish")!=null) this.LM_thresholdFinish= Double.parseDouble(properties.getProperty(prefix+"LM_thresholdFinish"));
			if (properties.getProperty(prefix+"LM_numIterations")!=  null) this.LM_numIterations=   Integer.parseInt(properties.getProperty(prefix+"LM_numIterations"));
			
			if (properties.getProperty(prefix+"fatZero")!=null)            this.fatZero=             Double.parseDouble(properties.getProperty(prefix+"fatZero"));

			if (properties.getProperty(prefix+"margin_left")!=  null)      this.margins[0]=           Integer.parseInt(properties.getProperty(prefix+"margin_left"));
			if (properties.getProperty(prefix+"margin_right")!= null)      this.margins[1]=           Integer.parseInt(properties.getProperty(prefix+"margin_right"));
			if (properties.getProperty(prefix+"margin_top")!=   null)      this.margins[2]=           Integer.parseInt(properties.getProperty(prefix+"margin_bottom"));
			if (properties.getProperty(prefix+"margin_bottom")!=null)      this.margins[3]=           Integer.parseInt(properties.getProperty(prefix+"margin_left"));
			if (properties.getProperty(prefix+"decimate")!=     null)      this.decimate=             Integer.parseInt(properties.getProperty(prefix+"decimate"));
			if (properties.getProperty(prefix+"sampleWidth")!=  null)      this.sampleWidth=          Integer.parseInt(properties.getProperty(prefix+"sampleWidth"));
			if (properties.getProperty(prefix+"highPassSigma")!=  null)    this.highPassSigma=        Double.parseDouble(properties.getProperty(prefix+"highPassSigma"));
			if (properties.getProperty(prefix+"maxTilt")!=  null)          this.maxTilt=              Double.parseDouble(properties.getProperty(prefix+"maxTilt"));
			
			if (properties.getProperty(prefix+"flatFieldDirectory")!=null) this.flatFieldDirectory=   properties.getProperty(prefix+"flatFieldDirectory");
			if (properties.getProperty(prefix+"eyesisMode")!=null)         this.eyesisMode=           Boolean.parseBoolean(properties.getProperty(prefix+"eyesisMode"));
			if (properties.getProperty(prefix+"processAllChannels")!=null) this.processAllChannels=   Boolean.parseBoolean(properties.getProperty(prefix+"processAllChannels"));
			if (properties.getProperty(prefix+"processChannels11")!=null)  this.processChannels[0]=   Boolean.parseBoolean(properties.getProperty(prefix+"processChannels11"));
			if (properties.getProperty(prefix+"processChannels12")!=null)  this.processChannels[1]=   Boolean.parseBoolean(properties.getProperty(prefix+"processChannels12"));
			if (properties.getProperty(prefix+"processChannels13")!=null)  this.processChannels[2]=   Boolean.parseBoolean(properties.getProperty(prefix+"processChannels13"));
			if (properties.getProperty(prefix+"processChannels21")!=null)  this.processChannels[3]=   Boolean.parseBoolean(properties.getProperty(prefix+"processChannels21"));
			if (properties.getProperty(prefix+"processChannels22")!=null)  this.processChannels[4]=   Boolean.parseBoolean(properties.getProperty(prefix+"processChannels22"));
			if (properties.getProperty(prefix+"processChannels23")!=null)  this.processChannels[5]=   Boolean.parseBoolean(properties.getProperty(prefix+"processChannels23"));
			if (properties.getProperty(prefix+"processChannels31")!=null)  this.processChannels[6]=   Boolean.parseBoolean(properties.getProperty(prefix+"processChannels31"));
			if (properties.getProperty(prefix+"processChannels32")!=null)  this.processChannels[7]=   Boolean.parseBoolean(properties.getProperty(prefix+"processChannels32"));
			if (properties.getProperty(prefix+"processChannels33")!=null)  this.processChannels[8]=   Boolean.parseBoolean(properties.getProperty(prefix+"processChannels33"));
			if (properties.getProperty(prefix+"useXML")!=null)             this.useXML=               Boolean.parseBoolean(properties.getProperty(prefix+"useXML"));
			if (properties.getProperty(prefix+"saveSettings")!=null)       this.saveSettings=         Boolean.parseBoolean(properties.getProperty(prefix+"saveSettings"));
			int j=0;
			if (properties.getProperty(prefix+"sourceDirPaths_length")!=null) j=Integer.parseInt(properties.getProperty(prefix+"sourceDirPaths_length"));
			if (j==0) this.sourceDirPaths=null;
			else {
			  this.sourceDirPaths=  new String[j];
			  this.sourceDirPathsEn=new boolean[j];
			  for (int i=0;i<j;i++) {
				this.sourceDirPaths[i]=properties.getProperty(prefix+"sourceDirPaths_"+i);
				if (properties.getProperty(prefix+"sourceDirPathsEn_"+i)!=null)
					this.sourceDirPathsEn[i]=Boolean.parseBoolean(properties.getProperty(prefix+"sourceDirPathsEn_"+i));
				else
					this.sourceDirPathsEn[i]=true;
			  }	
			}

		}

    }

    
    public static class ProcessCalibrationFilesParameters {
		public String [] subdirNames={"1-1","1-2","1-3","2-1","2-2","2-3","3-1","3-2","3-3"};
		public String [] suffixes=   {"11","12","13","21","22","23","31","32","33"};
		public String combinedSubDirectory="combined";
		public String sourceFileExtension="jp4"; 
		public String kernelFileExtension="tiff"; 
		public String kernelFilePrefix;
		public String psfRawPrefix;
		public String psfInterpoaltedPrefix;
		public String rpsfPrefix;
		public String gaussianPrefix;
		
		public String  sourceSuperDirectory; // having 1-1, 1-2,.. 3-3 subdirs
		public String  partialKernelsSuperDirectory;  // having 1-1, 1-2,.. 3-3 subdirs
		public String  kernelsDirectory;              // results with Gaussian and deconvolution kernels for all channels
		
		public boolean processAllChannels;            // if true - process all channels, otherwise only enabled in processChannels[]       
		public boolean [] processChannels=new boolean [9];
		public boolean keepOld;                       // do not re-calculate existent partial kernels, only the new ones 
		public boolean selectFiles;                   // select individual files to process
		
		public boolean processSourceImages;
		public boolean combinePSFfiles;               // combine partial PSF kernels 
		public boolean interpolatePSFkernel;          // interpolate PSF kernels (fail if missing??)
		public boolean invertKernels;                 // invert interpolated kernels
		public boolean gaussianKernels;               // create Gaussian kernels
		public boolean useXML;                        // save/restore settings as xml file
		public boolean saveSettings;                  // save current settings in results directory

		
		public ProcessCalibrationFilesParameters(
				String sourceFileExtension, 
				String kernelFileExtension,
				String kernelFilePrefix,
				
				String psfRawPrefix,
				String psfInterpoaltedPrefix,
				String rpsfPrefix,
				String gaussianPrefix,

				
				String  sourceSuperDirectory, // having 1-1, 1-2,.. 3-3 subdirs
				String  partialKernelsSuperDirectory, // having 1-1, 1-2,.. 3-3 subdirs
				String  kernelsDirectory,              // results with Gaussian and deconvolution kernels for all channels
				boolean processAllChannels,           // if true - process all channels, otherwise only enabled in processChatnnels[]       
				boolean processChannels11,
				boolean processChannels12,
				boolean processChannels13,
				boolean processChannels21,
				boolean processChannels22,
				boolean processChannels23,
				boolean processChannels31,
				boolean processChannels32,
				boolean processChannels33,
				boolean keepOld,                       // do not re-calculate existent partial kernels, only the new ones 
				boolean selectFiles,                   // select individual files to process
				
				boolean processSourceImages,           // process source calibration files		
				boolean combinePSFfiles,               // combine partial PSF kernels 
				boolean interpolatePSFkernel,          // interpolate PSF kernels (fail if missing??)
				boolean invertKernels,                 // invert interpolated kernels
				boolean gaussianKernels,                // create Gaussian kernels
  		        boolean useXML,
				boolean saveSettings
		) {
			this.sourceFileExtension=  sourceFileExtension;
			this.kernelFileExtension=  kernelFileExtension;
			this.kernelFilePrefix=     kernelFilePrefix;

			this.psfRawPrefix=         psfRawPrefix;
			this.psfInterpoaltedPrefix=psfInterpoaltedPrefix;
			this.rpsfPrefix=           rpsfPrefix;
			this.gaussianPrefix=       gaussianPrefix;
			
			this.sourceSuperDirectory= sourceSuperDirectory;
			this.partialKernelsSuperDirectory=partialKernelsSuperDirectory;
			this.kernelsDirectory=     kernelsDirectory;
			
			this.processAllChannels=   processAllChannels;
			this.processChannels[0]=   processChannels11;
			this.processChannels[1]=   processChannels12;
			this.processChannels[2]=   processChannels13;
			this.processChannels[3]=   processChannels21;
			this.processChannels[4]=   processChannels22;
			this.processChannels[5]=   processChannels23;
			this.processChannels[6]=   processChannels31;
			this.processChannels[7]=   processChannels32;
			this.processChannels[8]=   processChannels33;
			this.keepOld=keepOld;
			this.selectFiles=selectFiles;
			this.processSourceImages=processSourceImages;
			this.combinePSFfiles=combinePSFfiles; 
			this.interpolatePSFkernel=interpolatePSFkernel;
			this.invertKernels=invertKernels;
			this.gaussianKernels=gaussianKernels;
			this.useXML=useXML;
			this.saveSettings=saveSettings;
		}

		public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"sourceFileExtension",this.sourceFileExtension);
			properties.setProperty(prefix+"kernelFileExtension",this.kernelFileExtension);
			properties.setProperty(prefix+"kernelFilePrefix",this.kernelFilePrefix);
			properties.setProperty(prefix+"psfRawPrefix",this.psfRawPrefix);
			properties.setProperty(prefix+"psfInterpoaltedPrefix",this.psfInterpoaltedPrefix);
			properties.setProperty(prefix+"rpsfPrefix",this.rpsfPrefix);
			properties.setProperty(prefix+"gaussianPrefix",this.gaussianPrefix);
			properties.setProperty(prefix+"sourceSuperDirectory",this.sourceSuperDirectory);
			properties.setProperty(prefix+"partialKernelsSuperDirectory",this.partialKernelsSuperDirectory);
			properties.setProperty(prefix+"kernelsDirectory",this.kernelsDirectory);

			properties.setProperty(prefix+"processAllChannels",this.processAllChannels+"");
			properties.setProperty(prefix+"processChannels11",this.processChannels[0]+"");
			properties.setProperty(prefix+"processChannels12",this.processChannels[1]+"");
			properties.setProperty(prefix+"processChannels13",this.processChannels[2]+"");
			properties.setProperty(prefix+"processChannels21",this.processChannels[3]+"");
			properties.setProperty(prefix+"processChannels22",this.processChannels[4]+"");
			properties.setProperty(prefix+"processChannels23",this.processChannels[5]+"");
			properties.setProperty(prefix+"processChannels31",this.processChannels[6]+"");
			properties.setProperty(prefix+"processChannels32",this.processChannels[7]+"");
			properties.setProperty(prefix+"processChannels33",this.processChannels[8]+"");
			properties.setProperty(prefix+"keepOld",this.keepOld+"");
			properties.setProperty(prefix+"selectFiles",this.selectFiles+"");
			properties.setProperty(prefix+"processSourceImages",this.processSourceImages+"");
			properties.setProperty(prefix+"combinePSFfiles",this.combinePSFfiles+"");
			properties.setProperty(prefix+"interpolatePSFkernel",this.interpolatePSFkernel+"");
			properties.setProperty(prefix+"invertKernels",this.invertKernels+"");
			properties.setProperty(prefix+"gaussianKernels",this.gaussianKernels+"");
			properties.setProperty(prefix+"useXML",this.useXML+"");
			properties.setProperty(prefix+"saveSettings",this.saveSettings+"");
		}
		public void getProperties(String prefix,Properties properties){
			this.sourceFileExtension=         properties.getProperty(prefix+"sourceFileExtension");
			this.kernelFileExtension=         properties.getProperty(prefix+"kernelFileExtension");
			this.kernelFilePrefix=            properties.getProperty(prefix+"kernelFilePrefix");
			this.psfRawPrefix=                properties.getProperty(prefix+"psfRawPrefix");
			this.psfInterpoaltedPrefix=       properties.getProperty(prefix+"psfInterpoaltedPrefix");
			this.rpsfPrefix=                  properties.getProperty(prefix+"rpsfPrefix");
			this.gaussianPrefix=              properties.getProperty(prefix+"gaussianPrefix");
			this.sourceSuperDirectory=        properties.getProperty(prefix+"sourceSuperDirectory");
			this.partialKernelsSuperDirectory=properties.getProperty(prefix+"partialKernelsSuperDirectory");
			this.kernelsDirectory=            properties.getProperty(prefix+"kernelsDirectory");
			
			this.processAllChannels=   Boolean.parseBoolean(properties.getProperty(prefix+"processAllChannels"));
			this.processChannels[0]=   Boolean.parseBoolean(properties.getProperty(prefix+"processChannels11"));
			this.processChannels[1]=   Boolean.parseBoolean(properties.getProperty(prefix+"processChannels12"));
			this.processChannels[2]=   Boolean.parseBoolean(properties.getProperty(prefix+"processChannels13"));
			this.processChannels[3]=   Boolean.parseBoolean(properties.getProperty(prefix+"processChannels21"));
			this.processChannels[4]=   Boolean.parseBoolean(properties.getProperty(prefix+"processChannels22"));
			this.processChannels[5]=   Boolean.parseBoolean(properties.getProperty(prefix+"processChannels23"));
			this.processChannels[6]=   Boolean.parseBoolean(properties.getProperty(prefix+"processChannels31"));
			this.processChannels[7]=   Boolean.parseBoolean(properties.getProperty(prefix+"processChannels32"));
			this.processChannels[8]=   Boolean.parseBoolean(properties.getProperty(prefix+"processChannels33"));
			this.keepOld=              Boolean.parseBoolean(properties.getProperty(prefix+"keepOld"));
			this.selectFiles=          Boolean.parseBoolean(properties.getProperty(prefix+"selectFiles"));
			this.processSourceImages=  Boolean.parseBoolean(properties.getProperty(prefix+"processSourceImages"));
			this.combinePSFfiles=      Boolean.parseBoolean(properties.getProperty(prefix+"combinePSFfiles")); 
			this.interpolatePSFkernel= Boolean.parseBoolean(properties.getProperty(prefix+"interpolatePSFkernel"));
			this.invertKernels=        Boolean.parseBoolean(properties.getProperty(prefix+"invertKernels"));
			this.gaussianKernels=      Boolean.parseBoolean(properties.getProperty(prefix+"gaussianKernels"));
			this.useXML=               Boolean.parseBoolean(properties.getProperty(prefix+"useXML"));
			this.saveSettings=         Boolean.parseBoolean(properties.getProperty(prefix+"saveSettings"));
			
		}

	
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

}

