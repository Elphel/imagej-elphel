package com.elphel.imagej.calibration;
/**
 **
 ** LensAdjustment.jave - processing related to focus measurement/adjustment machine 
 **
 ** Copyright (C) 2011 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  LensAdjustment.java is free software: you can redistribute it and/or modify
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
import java.awt.Rectangle;
import java.util.Properties;

import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.dp.SimulationPattern;
import com.elphel.imagej.dp.WindowTools;
import com.elphel.imagej.dp.SimulationPattern.SimulParameters;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;

// Started to move methods from Aberration_Calibration
public class LensAdjustment {
	private ShowDoubleFloatArrays sdfaInstance =new ShowDoubleFloatArrays(); // just for debugging?
//	public int debugLevel=2;
	
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
        boolean noMessageBoxes=true;
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
		
		distortionParameters.correlationMinContrast=focusMeasurementParameters.correlationMinContrast;
		distortionParameters.correlationMinInitialContrast=focusMeasurementParameters.correlationMinInitialContrast;
		distortionParameters.correlationMinAbsoluteContrast=focusMeasurementParameters.correlationMinAbsoluteContrast;
		distortionParameters.correlationMinAbsoluteInitialContrast=focusMeasurementParameters.correlationMinAbsoluteInitialContrast;
		
		if (maskNonPSF) {
			distortionParameters.numberExtrapolated=0; //1; //3; // measuring PSF - extrapolate
		} else {
			distortionParameters.numberExtrapolated=1; // measuring distortions - do not extrapolate
		}
		
//System.out.println("distortionParameters.correlationSize="+distortionParameters.correlationSize);
		// add more overwrites
//		boolean updating=(matchSimulatedPattern.PATTERN_GRID!=null) &&
		if ((matchSimulatedPattern.PATTERN_GRID!=null) && (matchSimulatedPattern.PATTERN_GRID.length==0)){
			System.out.println("\n**** LensAdjusting.updateFocusGrid(): possible BUG - (matchSimulatedPattern.PATTERN_GRID!=null) && (matchSimulatedPattern.PATTERN_GRID.length==0)\n");
		}
		boolean updating=((matchSimulatedPattern.PATTERN_GRID!=null) && (matchSimulatedPattern.PATTERN_GRID.length>0)) &&
		(!distortionParameters.flatFieldCorrection || (matchSimulatedPattern.flatFieldForGrid!=null));
		if (debug_level>1) System.out.println("updateFocusGrid(): updating="+updating+" numberExtrapolated="+distortionParameters.numberExtrapolated);

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
			if (debug_level>0) { // calculate/print number of defined nodes in the grid
				System.out.println("updateFocusGrid(), number of already defined grid cells (before distortions()) = "+matchSimulatedPattern.numDefinedCells());
			}
			matchSimulatedPattern.distortions(
					null, // is not used in update mode
					distortionParameters, //
					patternDetectParameters,
					simulParameters,
					equalizeGreens,
					imp_eq, // image to process
					threadsMax,
					updateStatus,
					debug_level,// debug level used inside loops
					debug_level); // global debug level ?
			if (debug_level>1) System.out.println("distortions() finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));

			matchSimulatedPattern.recalculateWaveVectors (
					updateStatus,
					debug_level);// debug level used inside loops
			if (debug_level>1) System.out.println("recalculateWaveVectors() finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
			if (debug_level>1) { // calculate/print number of defined nodes ina grid
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
					null, // double [][][] hintGrid
					0,    // hintTolerance
					threadsMax,
					updateStatus,
					debug_level,
					distortionParameters.loop_debug_level, // debug level
					noMessageBoxes);
			if (numAbsolutePoints<0){
				System.out.println ("updateFocusGrid() numAbsolutePoints="+numAbsolutePoints);
			}

		}
		if (maskNonPSF) {
			matchSimulatedPattern.maskFocus(
	    			x0,   // lens center on the sensor
	    			y0,  // lens center on the sensor
					focusMeasurementParameters);
			if (debug_level>1) System.out.println("maskFocus() finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
			if (debug_level>1) { // calculate/print number of defined nodes ina grid
				System.out.println("number of defined grid cells (after maskFocus()) = "+matchSimulatedPattern.numDefinedCells());
			}
			if(debug_level>2) {
			       double [] test_masked=new double [matchSimulatedPattern.focusMask.length];
			       float [] pixels_eq=(float []) imp_eq.getProcessor().getPixels();
				   for (int i=0;i<test_masked.length;i++) test_masked[i]=matchSimulatedPattern.focusMask[i]?pixels_eq[i]:0.0;
				   sdfaInstance.showArrays(test_masked,matchSimulatedPattern.getImageWidth(), matchSimulatedPattern.getImageHeight(), "MASKED");
			}

		}
		if (debug_level>1) System.out.println("updateFocusGrid() finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		if (debug_level>2) {
	       double [] test_uv=new double [matchSimulatedPattern.UV_INDEX.length];
		   for (int i=0;i<matchSimulatedPattern.UV_INDEX.length;i++) test_uv[i]=matchSimulatedPattern.UV_INDEX[i];
		   sdfaInstance.showArrays(test_uv,matchSimulatedPattern.getImageWidth(), matchSimulatedPattern.getImageHeight(), "UV_INDEX");
		}
		return numAbsolutePoints;

	}
	
    public static class FocusMeasurementParameters {
    	public String gridGeometryFile="";
    	public String initialCalibrationFile="";
    	public String focusingHistoryFile="";
    	public boolean useLMAMetrics=true; // measure/report focal distance and tilts using lens model/LMA (when available)
    	public String strategyFile="";
    	public String resultsSuperDirectory=""; // directory with subdirectories named as serial numbers to stro results
    	public int EEPROM_channel=1; // EEPROM channel to read serial number from
    	public boolean saveResults=true; // save focusing results
    	public boolean showResults=true; // show focusing (including intermediate) results

    	public String serialNumber=""; // camera serial number string
    	public double sensorTemperature=Double.NaN; // last measured sensor temperature
//other summary results to be saved with parameters
    	public double result_lastKT=Double.NaN;   // focal distance temperature coefficient (um/C), measured from last run 
    	public double result_lastFD20=Double.NaN; // focal distance for 20C, measured from last run 
    	public double result_allHistoryKT=Double.NaN;   // focal distance temperature coefficient (um/C), measured from all the measurement histgory 
    	public double result_allHistoryFD20=Double.NaN; // focal distance for 20C, measured from  all the measurement histgory
    	public double result_fDistance=Double.NaN; // last measured focal distance
    	public double result_tiltX=Double.NaN; // last measured tilt X
    	public double result_tiltY=Double.NaN; // last measured tilt Y
    	public double result_R50=Double.NaN;   // last measured R50 (average PSF 50% level radius, pixels - somewhat larged than actual because of measurement settings)
    	public double result_A50=Double.NaN;   // last measured A50 (simailar, but R^2 are averaged) 
    	public double result_B50=Double.NaN;   // last measured B50 (simailar, but R^4 are averaged)
    	public double result_RC50=Double.NaN;  // last measured RC50(R50 calculated only for the 2 center samples)
    	public double result_PX0=Double.NaN; // lens center shift, X
    	public double result_PY0=Double.NaN; // lens center shift, Y
    	public double result_PSI=Double.NaN; // SFE rotation (from grid)
    	public double result_ROT=Double.NaN; // SFE rotation (from head lasers)
    	public double result_FocalLength=Double.NaN; // lens focal length

    	public String comment="no comments"; // Comment to add to the results
    	public int lensSerialLength=4;
    	public String lensSerial=""; // Lens serial number
    	public boolean askLensSerial=true; // Ask lens serial on camera power cycle
    	public double reportTemperature=50;      // temperature to report focal length  
    	public boolean showLegacy=false; // Show old focusing parameters (most are not supported anyway)
    	public boolean includeLensSerial=true; // add lens serial to config filename
    	public double centerDeltaX=0.0; // required X-difference between lens center and sensor center 
    	public double centerDeltaY=0.0; // required Y-difference between lens center and sensor center
    	// with the seam in the middle - make even # of samples horizontally
    	// to be made private or removed
//    	public Rectangle margins=new Rectangle (100,100,2392,1736) ; // maximal height 1816 (1936-120)
    	private Rectangle margins=new Rectangle (100,100,2392,1736) ; // maximal height 1816 (1936-120)
    	public int [] numSamples={8,5};       // number of samples in x and y directions
    	public int    sampleSize=256;// 512;       // size of square (2^n), in sensor pixels (twice FFT size)
    	public int    numInCenter=(numSamples[0]-2)*(numSamples[1]-2);// 2 - number of "center" samples
    	public boolean centerSamples=true; // Select samples in the WOI symmetrical around the lens center

    	public double maxCorr=5;     // maximal grid correction between images allowed (larger will trigger full grid rebuild)
    	public boolean showHistoryDetails=false;   // show color info
    	public boolean showHistorySamples=true;   // show individual samples
    	public boolean showHistorySingleLine=true; // all parameters in a single line (easier to copy to spreadsheet)
        public boolean showAcquiredImages=false;
        public boolean showROI=true;
        public boolean showSamples=true;
        
        public boolean showFittedParameters=true;
        public boolean useHeadLasers=true;
        
        // when approximating PSF with a second degree polynomial:
        public double psf_cutoffEnergy=0.98; //0.5; // disregard pixels outside of this fraction of the total energy
        public double psf_cutoffLevel= 0.2; // disregard pixels below this fraction of the maximal value
        public int    psf_minArea    = 10;  // continue increasing the selected area, even if beyond psf_cutoffEnergy and psf_cutoffLevel,
                                             // if the selected area is smaller than this (so approximation would work)
        public double psf_blurSigma  = 0.0; // optionally blur the calculated mask
        
        public double weightRatioRedToGreen=0.7;  // Use this data when combining defocusing data from different color PSF
        public double weightRatioBlueToGreen=0.2;
        public double targetFarNear=0.0;         // OBSOLETE target logarithm of average tangential-to-radial resolution
        public boolean useRadialTangential=false; // Use targetFarNear (radial/tangential resolution)  as a proxy for the distance  
        public double targetMicrons=0.0;          // target lens center distance (away from "best focus"
        public double toleranceMicrons=0.5; // microns
        public double toleranceTilt=0.02; // 
        public double toleranceThreshold=3.0; // When each error is under scaled threshold, reduce correction step twice
//        public boolean parallelAdjust=true;   // move 3 motors parallel after each 3-motor focus/tilt adjustment   
        public double  parallelAdjustThreshold=4.0;   // adjust 3 motors parallel if focal distance error in the center exceeds this    
        
        
        public double motorsSigma=896.0;   // when fitting planes for far/near, tiltX and tiltY the weights of the samples decay with this sigma 1/4 turn
        public double motorsSigma3=896.0;  // all 3 motors together (focusing center)
        public double motorsMinSigma=200.0;// sigma will not drop below this value when fitting walk is getting smaller 
        public double motorsVarSigmaToTravel=2.0;  // when walk is getting smaller, sigma will keep going down proportionally
        public double motorsFadeSigma=0.3;         // after each step new sigma will have this part of the calculated from the travel 
        public double motorsOverShootToBalance=0.9; //For quadratic maximum the correction will be increased by 1+motorsOverShootToBalance if there are less samples on teh other side
        public boolean filterGoodDistance=true; // when measuring tilt, use those with good center with higher weight
        public double goodDistanceSigma=5.0;    // sigma for the weight function of tilt measurements, depending on the center distance error
        public double goodTiltSigma=0.05;       // weight decay for heavily tilted samples
        public double maxStep=2000.0;       // maximal allowed single-step focusing adjustment
        public double probeStep=896.0;     // how far to go to probe around the current point to measure derivatives
        public double probe_M1M2M3= 448.0;    // how far to move average of the 3 motors: (M1+M2+M3)/3
        public double probe_M3_M1M2=3584.0;  // how far to move M3 opposite to M1 and M2: M3-(M1+M2)/2
        public double probe_M2_M1=  3584.0;  // how far to move M2 opposite to M1:    M2-M1
        public double sigmaToProbe=1.0;      // data from far samples decay proportionally to the probe distances 
        public boolean useTheBest=true;      // adjust from the best known position (false - from the last)  
        
        
        public boolean probeSymmetrical=false; // if true, probe 6 measurements), if false - only 4  (tetrahedron)  
        public boolean parallelBeforeProbing=true; // move 3 motors before probing around  
        public double reProbeDistance=4000.0; // re-run probing around in orthoganal directions, if the current position move farther from the last probing one
    	public double believeLast=0.0;       // coefficient 0.. 1.0. When each ot the 3 parameters is linearized, add shift, so 1.0 the planes will go throug the last sample
    	
    	public boolean compensateHysteresis=true; // move motors in the same direction to compensate for hysteresis
    	public double minCorr=100.0;              // minimal correction movement to initiate final numFinalCorr moves
    	public int    numFinalCorr=5;             // exit if this number of last corrections where below  minCorr
    	public double minCorrPre=200.0;           // minimal correction movement to initiate final numFinalCorr moves
    	public int    numFinalCorrPre=2;          // exit if this number of last corrections where below  minCorr
    	public int    maxAutoIterations=25;      // exit if history grows above this
    	public double maxAutoDistance=  10000;     // Maximal allowed automatic correction, motor steps
    	public boolean confirmFirstAuto=true;     // ask confirmation after first automatic adjustment step (before moving)
    	
        public double motorsPreSigma=3584.0; // when fitting parabola for focusing sharpness in the center, far measurements decay with this sigma
        public double maxLinearStep= 3584.0; // If there are insufficient measurements to fit parabola - make this step
        
        public int scanStep=320; // 200;             // motor step (all 3 motors) in scan focus mode (signed value)
        public int scanNumber=50;            // number of scanStep steps to run
        public int scanNumberNegative=20; // 15;    // number of scanStep steps negative from the start 
        public boolean scanHysteresis=false; // true;  // scan both ways
        public int scanHysteresisNumber=5;   // number of test points for the Hysteresis measurement
        
        public boolean scanTiltEnable=false; //true;  // enable scanning tilt
        public boolean scanTiltReverse=false;  // enable scanning tilt in both directions
        public boolean scanMeasureLast=false;  // Calculate PSF after last move (to original position)
        public boolean scanRunLMA=true;  // Calculate PSF after last move (to original position)
        public int scanTiltRangeX=14336;    // 4 periods
        public int scanTiltRangeY=14336;    // 4 periods
        public int scanTiltStepsX=24;
        public int scanTiltStepsY=24;
        
        
        public int motorHysteresis=300;
        public double measuredHysteresis=0; // actually measured (will be saved/restored)
        public double motorCalm=2; // wait (seconds) after motors reached final position (for the first time) before acquiring image   
        public double linearReductionRatio=4.0/38.0; // sensor travel to motors travel (all 3 together), By design it is 4/38~=0.105
        public int motorDebug=0; // 1 show  motor moves, 2 - show hysteresis back-ups too
        
        // parameters for appoximating sensor center position
        public int    lensDistanceNumPoints=1000; // number of points to tabulate center focus parameters vs. focal distance 
        public int    lensDistancePolynomialDegree=8; // polynomial degree to approximate center focus parameters vs. focal distance 
        public double lensDistanceWeightY=0.5;     // normalize overall sharpness (that depends on the lens quality and/or PSF parameters) to differential ones 
        public double lensDistanceWeightK=0.0;     // 0.0 use 3 distances (to frac-R, frac-B and Y) with teh same weight when fitting, 1.0 - proportional to squared derivatives
        public boolean lensDistanceInteractive=true; // Open dialog when calibrating focal distance
        public boolean lensDistanceShowResults=true; // show results window from foca
        public boolean lensDistanceMoveToGoal=true;  // Move to targetMicrons
        
        public boolean powerControlEnable=true;
        public double powerControlMaximalTemperature=60.0;
        public double powerControlHeaterOnMinutes=10.0;
        public double powerControlNeitherOnMinutes=5.0;
        public double powerControlFanOnMinutes=15.0;
        
        public String uvLasersIP="192.168.0.236"; // IP address of the camera with UV LEDs and aiming lasers are connected
        public int    uvLasersBus=0;              // 0 if 103641 board is connected to the sensor port (through 10-359), 1 - to 10369
        public double [] uvLasersCurrents={40.0,40.0,40.0,40.0}; // default LED on currents (mA)
        
    	// the following 3 overwrite SimulParameters members
    	public double smallestSubPix=0.3; // subdivide pixels down to that fraction when simulating
    	public double bitmapNonuniforityThreshold=0.1	; // subdivide pixels until difference between the corners is below this value
    	public int    subdiv=4; 
    	// overwrites  	public static class MultiFilePSF.overexposedMaxFraction
    	public double overexposedMaxFraction=0.1; // allowed fraction of the overexposed pixels in the PSF kernel measurement area 
    	// overwrites	public static class PSFParameters.minDefinedArea
    	public double minDefinedArea=0.75;    // minimal (weighted) fraction of the defined patter pixels in the FFT area
        public int PSFKernelSize=32;          // size of the detected PSF kernel
		public boolean approximateGrid=true; // approximate grid with polynomial 
		public boolean centerPSF=true;       // Center PSF by modifying phase
		
		public double mask1_sigma=    1.0;
		public double mask1_threshold=0.25;
		public double gaps_sigma=     1.0;
		public double mask_denoise=   0.25;
    //  OTFFilterParameters
		public double deconvInvert=0.03; // with good focus can go to 0.015 or smaller
    //  DistortionParameters
        public int correlationSize=32;
        public double correlationGaussWidth=0.75;
        public double minUVSpan;           // Minimal u/v span in correlation window that triggers increase of the correlation FFT size

		public double correlationMinContrast=1.2;   // minimal contrast for the pattern to pass
		public double correlationMinInitialContrast=1.5;   // minimal contrast for the pattern of the center (initial point)
		public double correlationMinAbsoluteContrast=0.5;   // minimal contrast for the pattern to pass, does not compensate for low ligt
		public double correlationMinAbsoluteInitialContrast=0.5;   // minimal contrast for the pattern of the center (initial point)
        
		public double [] postUVscrewSensitivity={-300,-300,-120}; //{-2.908571735,-3.8198374024,-2.4491867448};
        
        public boolean flatFieldCorrection=true;
        public double flatFieldExpand=4.0;
        
        public double thresholdFinish=0.001; // Stop iterations if 2 last steps had less improvement (but not worsening ) 
        public int    numIterations=  100; // maximal number of iterations
        
        public boolean cameraIsConfigured=false;
        public boolean useExtraSensor = true;
        public boolean enRoundOff = true; // create round mask to refine lens center position (second LMA run) 
        public boolean keepCircularMask=false; // keep sensor 0 round mask after running LMA, false - recalculate actual mask
        public boolean configureCamera=false; // only valid after dialog
        public int [] motorPos=null; // will point to
    	public double [] ampsSeconds={0.0,0.0,0.0,0.0}; // cumulative Amps*seconds (read only, but will be saved/restored)
    	public int manufacturingState=0;
    	
    	// temporary - will re-calculate according to required lens center
    	public Rectangle getMargins(){
    		int iXC = (int) Math.round(centerDeltaX);
    		int iYC = (int) Math.round(centerDeltaY);
    		Rectangle margins=new Rectangle(this.margins);
    		if ((iXC !=0) || (iYC!=0)){
    			int hw=this.margins.width/2;
    			int hh=this.margins.height/2;
    			int ix0=this.margins.x+hw;
    			int iy0=this.margins.y+hh;
    			hw += (iXC < 0)?iXC:(-iXC);
    			hh += (iYC < 0)?iYC:(-iYC);
    			margins=new Rectangle(ix0+iXC-hw,iy0+iYC-hh,2*hw,2*hh);
    		}
    		
    		return margins;
    	}
    	
    	public String [] manufacturingStateNames={
    			"New SFE",
    			"UV cured (not released)",
    			"UV cured (released)",
    			"UV cured, re-tested later",
    			"Epoxy, not cured",
    			"Epoxy, cured by thermo-cycling (not released)",
    			"Epoxy, cured by thermo-cycling (released)",
    			"Epoxy, cured at room temperature",
    			"Epoxy cured, re-testing later"
    			};
    			
    	public int [] manufacturingStateValues={
    			00, // "New SFE",
    			10, //"UV cured (not released)",
    			20, //"UV cured (released)",
    			30, //"UV cured, re-tested later",
    			40, //"Epoxy, not cured",
    			50, //"Epoxy, cured by thermo-cycling (not released)",
    			60, //"Epoxy, cured by thermo-cycling (released)",
    			70, //"Epoxy, cured at room temperature",
    			80  //"Epoxy cured, re-testing later"
    	};
    	
    	public int [] getManufacturingIndexMod(int value){
    		if (value<this.manufacturingStateValues[0]) value= this.manufacturingStateValues[0];
    		for (int i=this.manufacturingStateValues.length-1;i>=0;i--){
    			if (manufacturingStateValues[i]<=value) {
    				int [] result = {i,value-manufacturingStateValues[i]};
    				return result;
    			}
    		}
    		return null; // should not get here
    	}
    	
    	
    	public FocusMeasurementParameters(int []motorPos){
            this.motorPos=motorPos;
    	}
    	public void resetResults(){
	    	this.result_lastKT=Double.NaN;   // focal distance temperature coefficient (um/C), measured from last run 
	    	this.result_lastFD20=Double.NaN; // focal distance for 20C, measured from last run 
	    	this.result_allHistoryKT=Double.NaN;   // focal distance temperature coefficient (um/C), measured from all the measurement histgory 
	    	this.result_allHistoryFD20=Double.NaN; // focal distance for 20C, measured from  all the measurement histgory
	    	this.result_fDistance=Double.NaN; // last measured focal distance
	    	this.result_tiltX=Double.NaN; // last measured tilt X
	    	this.result_tiltY=Double.NaN; // last measured tilt Y
	    	this.result_R50=Double.NaN;   // last measured R50 (average PSF 50% level radius, pixels - somewhat larged than actual because of measurement settings)
	    	this.result_A50=Double.NaN;   // last measured A50 (simailar, but R^2 are averaged) 
	    	this.result_B50=Double.NaN;   // last measured B50 (similar, but R^4 are averaged)
	    	this.result_RC50=Double.NaN;  // last measured RC50(R50 calculated only for the 2 center samples)
	    	this.result_PX0=Double.NaN; // lens center shift, X
	    	this.result_PY0=Double.NaN; // lens center shift, Y
	    	this.result_PSI=Double.NaN; // SFE rotation (from grid)
	    	this.result_ROT=Double.NaN; // SFE rotation (from head lasers)
	    	this.result_FocalLength=Double.NaN; // lens focal length
    	}
    	public FocusMeasurementParameters(
    	    	String gridGeometryFile,
    	    	String initialCalibrationFile,
    	    	String strategyFile,
    	    	String resultsSuperDirectory, // directory with subdirectories named as serial numbers to stro results
    	    	String focusingHistoryFile,
    	    	boolean useLMAMetrics, // measure/report focal distance and tilts using lens model/LMA (when available)
    	    	int EEPROM_channel, // EEPROM channel to read serial number from
    	    	boolean saveResults, // save focusing results
    	    	boolean showResults, // show focusing (includingh intermediate) results
    	    	String serialNumber, // camera serial number string
    	    	double sensorTemperature, // last measured sensor temperature
    	    	// other summary results to be saved with parameters
    	    	double result_lastKT,   // focal distance temperature coefficient (um/C), measured from last run 
    	    	double result_lastFD20, // focal distance for 20C, measured from last run 
    	    	double result_allHistoryKT,   // focal distance temperature coefficient (um/C), measured from all the measurement histgory 
    	    	double result_allHistoryFD20, // focal distance for 20C, measured from  all the measurement histgory
    	    	double result_fDistance, // last measured focal distance
    	    	double result_tiltX, // last measured tilt X
    	    	double result_tiltY, // last measured tilt Y
    	    	double result_R50,   // last measured R50 (average PSF 50% level radius, pixels - somewhat larged than actual because of measurement settings)
    	    	double result_A50,   // last measured A50 (simailar, but R^2 are averaged) 
    	    	double result_B50,   // last measured B50 (simailar, but R^4 are averaged)
    	    	double result_RC50,  // last measured RC50(R50 calculated only for the 2 center samples)
    	    	double result_PX0, // lens center shift, X
    	    	double result_PY0, // lens center shift, Y
    	    	double result_PSI, // SFE rotation (from grid)
    	    	double result_ROT, // SFE rotation (from head lasers)
    	    	double result_FocalLength, // lens focal length (currently not calculated!)
    	    	String comment, // Comment to add to the results
    	    	String lensSerial, // Lens serial number
    	    	int manufacturingState, // SFE manufacturing state
    	    	boolean askLensSerial,
    	    	boolean includeLensSerial, // add lens serial to config filename
    	    	double centerDeltaX, // required X-difference between lens center and sensor center 
    	    	double centerDeltaY, // required Y-difference between lens center and sensor center
    	    	Rectangle margins,
    			int [] numSamples,
    			int sampleSize,
    	    	int numInCenter,// 2 - number of "center" samples
    	    	boolean centerSamples, // selsct samples in WOI symmetrical around the lens center
    			double maxCorr,     // maximal grid correction between images allowed (larger will trigger full grid rebuild)
    	    	boolean showHistoryDetails,
    	    	boolean showHistorySamples,
    	    	boolean showHistorySingleLine, // all parameters in a single line (easier to copy to spreadsheet)
    	    	boolean showAcquiredImages,
    	        boolean showROI,
    	        boolean showSamples,
    	    	
    	    	boolean showFittedParameters,
    	        boolean useHeadLasers,
                double psf_cutoffEnergy, // disregard pixels outside of this fraction of the total energy
                double psf_cutoffLevel,  // disregard pixels below this fraction of the maximal value
                int    psf_minArea,      // continue increasing the selected area, even if beyound psf_cutoffEnergy and psf_cutoffLevel,
                                         // if the selected area is smaller than this (so approximation wpuld work)
                double psf_blurSigma,    // optionally blur the calculated mask
                
                double weightRatioRedToGreen,  // Use this data when combining defocusing data from different color PSF
                double weightRatioBlueToGreen,
                double targetFarNear,         // target logariphm of average tangential-to-radial resolution
                boolean useRadialTangential, // Use targetFarNear (radial/tangential resolution)  as a proxy for the distance  
                double targetMicrons,          // target lens center distance (away from "best focus"
                double toleranceMicrons, // microns
                double toleranceTilt, // 
                double toleranceThreshold, // When each error is under swcaled thereshold, reduce correxction step twice
                double  parallelAdjustThreshold,   // adjust 3 motors parallel if focal distance error in the center exceeds this    
                double motorsSigma,   // when fitting planes for far/near, tiltX and tiltY the weights of the samples decay with this sigma
                double motorsSigma3,  // all 3 motors together (focusing center)
                double motorsMinSigma,// sigma will not drop below this value when fitting walk is getting smaller 
                double motorsVarSigmaToTravel,  // when walk is getting smaller, sigma will keep going down proportionally
                double motorsFadeSigma,       // after each step new sigma will have this part of the calculated from the travel 
                double motorsOverShootToBalance, //For quadratic maximum the correction will be increased by 1+motorsOverShootToBalance if there are less samples on teh other side
                boolean filterGoodDistance, // when measuring tilt, use those with good center with higher weight
                double goodDistanceSigma,    // sigma for the weight function of tilt measurements, depending on the center distance error
                double goodTiltSigma,       // weight decay for heavily tilted samples
                double maxStep,          // maximal allowed single-step focusing adjustment
                double probeStep,        // how far to go to probe around the current point to measure derivatives
                double probe_M1M2M3,     // how far to move average of the 3 motors: (M1+M2+M3)/3
                double probe_M3_M1M2,    // how far to move M3 opposite to M1 and M2: M3-(M1+M2)/2
                double probe_M2_M1,      // how far to move M2 opposite to M1:    M2-M1
                double sigmaToProbe,     // data from far samples decay proportionally to the probe distances 
                boolean useTheBest,      // adjust from the best known position (false - from the last)  	                
                boolean probeSymmetrical, // if true, probe 6 measurements), if false - only 4  (tetrahedron)
                boolean parallelBeforeProbing, // move 3 motors before probing around
                double reProbeDistance,   // re-run probing around in orthoganal directions, if the current position move farther from the last probing one
    	    	double believeLast,       // coefficient 0.. 1.0. When each ot the 3 parameters is linearized, add shift, so 1.0 the planes will go throug the last sample
    	    	boolean compensateHysteresis,// move motors in the same direction to compensate fro hysteresis
    	    	double minCorr,              // minimal correction movement to initiate final numFinalCorr moves
    	    	int    numFinalCorr,         // exit if this number of last corrections where below  minCorr
    	    	double minCorrPre,              // minimal correction movement to initiate final numFinalCorr moves
    	    	int    numFinalCorrPre,         // exit if this number of last corrections where below  minCorr
    	    	int    maxAutoIterations,   // exit if history grows above this
    	    	double maxAutoDistance,     // Maximal allowed automatic correction, motor steps
    	    	boolean confirmFirstAuto,     // ask confirmation after first automatic adjustment step (before moving)
                double motorsPreSigma, // when fitting parabola for focusing sharpness in the center, far measurements decay with this sigma
                double maxLinearStep,  // If there are insufficient measurements to fit parabola - make this step
                int scanStep,             // motor step (all 3 motors) in scan focus mode (signed value)
                int scanNumber,            // number of scanStep steps to run
                int scanNumberNegative,    // of them negative
                boolean scanHysteresis,  // scan both ways
                int scanHysteresisNumber,   // number of test points for the Hysteresis measurement
                
                boolean scanTiltEnable, //=true;  // enable scanning tilt
                boolean scanTiltReverse,
                boolean scanMeasureLast,
                boolean scanRunLMA,
                int scanTiltRangeX, //=14336;    // 4 periods
                int scanTiltRangeY, //=14336;    // 4 periods
                int scanTiltStepsX, //=24;
                int scanTiltStepsY, //=24;

                
                int motorHysteresis,
                double measuredHysteresis, // actually measured (will be saved/restored)
                double motorCalm, // wait (seconds) after motors reached final position (for the first time) before acquiring image   
                double linearReductionRatio, // sensor travel to motors travel (all 3 together), By design it is 4/38~=0.105
                int motorDebug,// 1 show  motor moves, 2 - show hysteresis back-ups too
                int    lensDistanceNumPoints, // number of points to tabulate center focus parameters vs. focal distance 
                int    lensDistancePolynomialDegree, // polynomial degree to approximate center focus parameters vs. focal distance 
                double lensDistanceWeightY,     // normalize overall sharpness (that depends on the lens quality and/or PSF parameters) to differential ones 
                double lensDistanceWeightK,     // 0.0 use 3 distances (to frac-R, frac-B and Y) with teh same weight when fitting, 1.0 - proportional to squared derivatives
                boolean lensDistanceInteractive, // Open dialog when calibrating focal distance
                boolean lensDistanceShowResults, // show results window from foca
                boolean lensDistanceMoveToGoal,  // Move to targetMicrons
                boolean powerControlEnable,
                double powerControlMaximalTemperature,
                double powerControlHeaterOnMinutes,
                double powerControlNeitherOnMinutes,
                double powerControlFanOnMinutes,
                String uvLasersIP, // IP address of the camera with UV LEDs and aiming lasers are connected
                int    uvLasersBus,             // 0 if 103641 board is connected to the sensor port (through 10-359), 1 - to 10369
                double [] uvLasersCurrents, // default LED on currents (mA)

    	    	double smallestSubPix, // subdivide pixels down to that fraction when simulating
    	    	double bitmapNonuniforityThreshold, // subdivide pixels until difference between the corners is below this value
    	    	int    subdiv, 
    	    	double overexposedMaxFraction, // allowed fraction of the overexposed pixels in the PSF kernel measurement area 
    	    	double minDefinedArea, // minimal (weighted) fraction of the defined patter pixels in the FFT area
    	    	int PSFKernelSize,
				boolean approximateGrid, // approximate grid with polynomial
				boolean centerPSF,       // Center PSF by modifying phase
				double mask1_sigma,
				double mask1_threshold,
				double gaps_sigma,
				double mask_denoise,
				double deconvInvert,
                int correlationSize,
                double correlationGaussWidth,
                double minUVSpan,           // Minimal u/v span in correlation window that triggers increase of the correlation FFT size
        		double correlationMinContrast,   // minimal contrast for the pattern to pass
        		double correlationMinInitialContrast,   // minimal contrast for the pattern of the center (initial point)
        		double correlationMinAbsoluteContrast,   // minimal contrast for the pattern to pass, does not compensate for low ligt
        		double correlationMinAbsoluteInitialContrast,   // minimal contrast for the pattern of the center (initial point)
        		double [] postUVscrewSensitivity, // microns/turn for 3 post-UV fixture adjustment screws
                boolean flatFieldCorrection,
                double flatFieldExpand,
                double thresholdFinish,// (copied from series) stop iterations if 2 last steps had less improvement (but not worsening ) 
                int    numIterations, // maximal number of iterations
				boolean cameraIsConfigured,
				boolean useExtraSensor,
				boolean enRoundOff,
				boolean keepCircularMask,
				int [] motorPos,
	        	double [] ampsSeconds, // cumulative Amps*seconds (read only, but will be saved/restored)
	        	double reportTemperature,      // temperature to report focal length
	        	boolean showLegacy
    			){
    		this.gridGeometryFile=gridGeometryFile;
    		this.initialCalibrationFile=initialCalibrationFile;
    		this.strategyFile=strategyFile;
    		this.resultsSuperDirectory=resultsSuperDirectory; // directory with subdirectories named as serial numbers to stro results
    		this.focusingHistoryFile=focusingHistoryFile;
    		this.useLMAMetrics=useLMAMetrics; // measure/report focal distance and tilts using lens model/LMA (when available)
    		this.EEPROM_channel=EEPROM_channel; // EEPROM channel to read serial number from
    		this.saveResults=saveResults; // save focusing results
    		this.showResults=showResults; // show focusing (includingh intermediate) results
    		this.serialNumber=serialNumber; // camera serial number string
    		this.sensorTemperature=sensorTemperature; // last measured sensor temperature
    		this.result_lastKT=result_lastKT;   // focal distance temperature coefficient (um/C), measured from last run 
    		this.result_lastFD20=result_lastFD20; // focal distance for 20C, measured from last run 
    		this.result_allHistoryKT=result_allHistoryKT;   // focal distance temperature coefficient (um/C), measured from all the measurement histgory 
    		this.result_allHistoryFD20=result_allHistoryFD20; // focal distance for 20C, measured from  all the measurement histgory
    		this.result_fDistance=result_fDistance; // last measured focal distance
    		this.result_tiltX=result_tiltX; // last measured tilt X
    		this.result_tiltY=result_tiltY; // last measured tilt Y
    		this.result_R50=result_R50;   // last measured R50 (average PSF 50% level radius, pixels - somewhat larged than actual because of measurement settings)
    		this.result_A50=result_A50;   // last measured A50 (simailar, but R^2 are averaged) 
    		this.result_B50=result_B50;   // last measured B50 (simailar, but R^4 are averaged)
    		this.result_RC50=result_RC50;  // last measured RC50(R50 calculated only for the 2 center samples)
    		this.result_PX0=result_PX0; // lens center shift, X
    		this.result_PY0=result_PY0; // lens center shift, Y
    		this.result_PSI=result_PSI; // SFE rotation (from grid)
    		this.result_ROT=result_ROT; // SFE rotation (from head lasers)
    		this.result_FocalLength=result_FocalLength; // lerns focal length
    		this.comment=comment; // Comment to add to the results
    		this.lensSerial=lensSerial; // Lens serial number
    		this.manufacturingState=manufacturingState;
    		this.askLensSerial=askLensSerial;
    		this.includeLensSerial=includeLensSerial; // add lens serial to config filename
        	this.centerDeltaX=centerDeltaX; // required X-difference between lens center and sensor center 
        	this.centerDeltaY=centerDeltaY; // required Y-difference between lens center and sensor center
    		this.margins=(Rectangle) margins.clone();
//			this.origin=origin.clone(); // top left corner
//			this.periods=periods.clone(); // period {x,y}of the sample points
			this.numSamples=numSamples.clone();
			this.sampleSize=sampleSize;
			this.numInCenter=numInCenter;// 2 - number of "center" samples
			this.centerSamples=centerSamples;
			this.maxCorr=maxCorr;
			this.showHistoryDetails=showHistoryDetails;
			this.showHistorySamples=showHistorySamples;
			this.showHistorySingleLine=showHistorySingleLine; // all parameters in a single line (easier to copy to spreadsheet)
			this.showAcquiredImages=showAcquiredImages;
			this.showROI=showROI;
			this.showSamples=showSamples;
			this.showFittedParameters=showFittedParameters;
			this.useHeadLasers=useHeadLasers;
			this.psf_cutoffEnergy=psf_cutoffEnergy;
			this.psf_cutoffLevel= psf_cutoffLevel;
			this.psf_minArea=     psf_minArea;
			this.psf_blurSigma=   psf_blurSigma;
			this.weightRatioRedToGreen=weightRatioRedToGreen;  // Use this data when combining defocusing data from different color PSF
			this.weightRatioBlueToGreen=weightRatioBlueToGreen;
			this.targetFarNear=targetFarNear;                  // target logariphm of average tangential-to-radial resolution
			this.useRadialTangential=useRadialTangential;      // Use targetFarNear (radial/tangential resolution)  as a proxy for the distance  
			this.targetMicrons=targetMicrons;                  // target lens center distance (away from "best focus"
			this.toleranceMicrons=toleranceMicrons; // microns
			this.toleranceTilt=toleranceTilt; // 
			this.toleranceThreshold=toleranceThreshold; // When each error is under swcaled thereshold, reduce correxction step twice
			this.parallelAdjustThreshold=parallelAdjustThreshold; // adjust 3 motors parallel if focal distance error in the center exceeds this    
			this.motorsSigma=motorsSigma;                      // when fitting planes for far/near, tiltX and tiltY the weights of the samples decay with this sigma
			this.motorsSigma3=motorsSigma3;                    // same when 3 motors move together
			this.motorsMinSigma=motorsMinSigma;                // sigma will not drop below this value when fitting walk is getting smaller 
			this.motorsVarSigmaToTravel=motorsVarSigmaToTravel;// when walk is getting smaller, sigma will keep going down proportionally
			this.motorsFadeSigma=motorsFadeSigma;              // after each step new sigma will have this part of the calculated from the travel 
			this.motorsOverShootToBalance=motorsOverShootToBalance;//For quadratic maximum the correction will be increased by 1+motorsOverShootToBalance if there are less samples on teh other side
			this.filterGoodDistance=filterGoodDistance; // when measuring tilt, use those with good center with higher weight
			this.goodDistanceSigma=goodDistanceSigma;    // sigma for the weight function of tilt measurements, depending on the center distance error
			this.goodTiltSigma=goodTiltSigma;            // weight decay for heavily tilted samples
			this.maxStep=maxStep;       // maximal allowed single-step focusing adjustment
			this.probeStep=probeStep;     // how far to go to probe around the current point to measure derivatives
			this.probe_M1M2M3=probe_M1M2M3;     // how far to move average of the 3 motors: (M1+M2+M3)/3
			this.probe_M3_M1M2=probe_M3_M1M2;    // how far to move M3 opposite to M1 and M2: M3-(M1+M2)/2
			this.probe_M2_M1=probe_M2_M1;      // how far to move M2 opposite to M1:    M2-M1
			this.sigmaToProbe=sigmaToProbe;     // data from far samples decay proportionally to the probe distances
			this.useTheBest=useTheBest;      // adjust from the best known position (false - from the last)  	                
			this.probeSymmetrical=probeSymmetrical; // if true, probe 6 measurements), if false - only 4  (tetrahedron)
			this.parallelBeforeProbing=parallelBeforeProbing; // move 3 motors before probing around
			this.reProbeDistance=reProbeDistance;   // re-run probing around in orthoganal directions, if the current position move farther from the last probing one
			this.believeLast=believeLast;       // coefficient 0.. 1.0. When each ot the 3 parameters is linearized, add shift, so 1.0 the planes will go throug the last sample
			this.compensateHysteresis=compensateHysteresis;// move motors in the same direction to compensate fro hysteresis
			this.minCorr=minCorr;              // minimal correction movement to initiate final numFinalCorr moves
			this.numFinalCorr=numFinalCorr;    // exit if this number of last corrections where below  minCorr
			this.minCorrPre=minCorrPre;              // minimal correction movement to initiate final numFinalCorr moves
			this.numFinalCorrPre=numFinalCorrPre;    // exit if this number of last corrections where below  minCorr
			this.maxAutoIterations=maxAutoIterations;           // exit if history grows above this
			this.maxAutoDistance=maxAutoDistance;     // Maximal allowed automatic correction, motor steps
			this.confirmFirstAuto=confirmFirstAuto;   // ask confirmation after first automatic adjustment step (before moving)
			this.motorsPreSigma=motorsPreSigma;       // when fitting parabola for focusing sharpness in the center, far measurements decay with this sigma
			this.maxLinearStep= maxLinearStep;        // If there are insufficient measurements to fit parabola - make this step
			this.scanStep=scanStep;                   // motor step (all 3 motors) in scan focus mode (signed value)
			this.scanNumber=scanNumber;               // number of scanStep steps to run
            this.scanNumberNegative=scanNumberNegative;    // of them negative
			this.scanHysteresis=scanHysteresis;       // scan both ways
			this.scanHysteresisNumber=scanHysteresisNumber; // number of test points for the Hysteresis measurement
			
			this.scanTiltEnable=scanTiltEnable; //=true;  // enable scanning tilt
			this.scanTiltReverse=scanTiltReverse;
			this.scanMeasureLast=scanMeasureLast;
			this.scanRunLMA=scanRunLMA;
			this.scanTiltRangeX=scanTiltRangeX; //, //=14336;    // 4 periods
			this.scanTiltRangeY=scanTiltRangeY; //, //=14336;    // 4 periods
			this.scanTiltStepsX=scanTiltStepsX; //=24;
			this.scanTiltStepsY=scanTiltStepsY; //=24;
			
			this.motorHysteresis=motorHysteresis;
			this.measuredHysteresis=measuredHysteresis; // actually measured (will be saved/restored)
			this.motorCalm=motorCalm; // wait (seconds) after motors reached final position (for the first time) before acquiring image   
			this.linearReductionRatio=linearReductionRatio; // sensor travel to motors travel (all 3 together), By design it is 4/38~=0.105
			this.motorDebug=motorDebug;// 1 show  motor moves, 2 - show hysteresis back-ups too
			this.lensDistanceNumPoints=lensDistanceNumPoints; // number of points to tabulate center focus parameters vs. focal distance 
			this.lensDistancePolynomialDegree=lensDistancePolynomialDegree; // polynomial degree to approximate center focus parameters vs. focal distance 
			this.lensDistanceWeightY=lensDistanceWeightY;     // normalize overall sharpness (that depends on the lens quality and/or PSF parameters) to differential ones 
			this.lensDistanceWeightK=lensDistanceWeightK;     // 0.0 use 3 distances (to frac-R, frac-B and Y) with teh same weight when fitting, 1.0 - proportional to squared derivatives
			this.lensDistanceInteractive=lensDistanceInteractive; // Open dialog when calibrating focal distance
			this.lensDistanceShowResults=lensDistanceShowResults; // show results window from foca
			this.lensDistanceMoveToGoal=lensDistanceMoveToGoal;  // Move to targetMicrons
			this.powerControlEnable=powerControlEnable;
			this.powerControlMaximalTemperature=powerControlMaximalTemperature;
			this.powerControlHeaterOnMinutes=powerControlHeaterOnMinutes;
			this.powerControlNeitherOnMinutes=powerControlNeitherOnMinutes;
			this.powerControlFanOnMinutes=powerControlFanOnMinutes;
            this.uvLasersIP=new String(uvLasersIP); // IP address of the camera with UV LEDs and aiming lasers are connected
            this.uvLasersBus=uvLasersBus; // 0 if 103641 board is connected to the sensor port (through 10-359), 1 - to 10369
            this.uvLasersCurrents=uvLasersCurrents.clone(); // default LED on currents (mA)
			this.smallestSubPix=smallestSubPix;
			this.bitmapNonuniforityThreshold=bitmapNonuniforityThreshold;
			this.subdiv=subdiv; 
			this.overexposedMaxFraction=overexposedMaxFraction; 
			this.minDefinedArea=minDefinedArea;
			this.PSFKernelSize=PSFKernelSize;
			this.approximateGrid = approximateGrid; // approximate grid with polynomial
			this.centerPSF = centerPSF; // approximate grid with polynomial 
			this.mask1_sigma = mask1_sigma;
			this.mask1_threshold = mask1_threshold;
			this.gaps_sigma=gaps_sigma;
			this.mask_denoise=mask_denoise;
			this.deconvInvert = deconvInvert;
			this.correlationSize=correlationSize;
			this.correlationGaussWidth=correlationGaussWidth;
			this.minUVSpan=minUVSpan;
			this.correlationMinContrast=correlationMinContrast;   // minimal contrast for the pattern to pass
			this.correlationMinInitialContrast=correlationMinInitialContrast;   // minimal contrast for the pattern of the center (initial point)
			this.correlationMinAbsoluteContrast=correlationMinAbsoluteContrast;   // minimal contrast for the pattern to pass, does not compensate for low ligt
			this.correlationMinAbsoluteInitialContrast=correlationMinAbsoluteInitialContrast;   // minimal contrast for the pattern of the center (initial point)
			this.postUVscrewSensitivity=postUVscrewSensitivity.clone(); // microns/turn for 3 post-UV fixture adjustment screws
			this.flatFieldCorrection=flatFieldCorrection;
			this.flatFieldExpand=flatFieldExpand;
			this.thresholdFinish=thresholdFinish;// (copied from series) stop iterations if 2 last steps had less improvement (but not worsening ) 
			this.numIterations=numIterations; // maximal number of iterations
            this.cameraIsConfigured=cameraIsConfigured;
            this.useExtraSensor=useExtraSensor;
            this.enRoundOff=enRoundOff;
            this.keepCircularMask=keepCircularMask;
            this.motorPos=motorPos;
        	this.ampsSeconds=ampsSeconds; // cumulative Amps*seconds (read only, but will be saved/restored)
        	this.reportTemperature=reportTemperature;
        	this.showLegacy=showLegacy;
    	}
    	public FocusMeasurementParameters clone(){
    		return new FocusMeasurementParameters(
    	    		this.gridGeometryFile,
    	    		this.initialCalibrationFile,
    	    		this.strategyFile,
    	    		this.resultsSuperDirectory, // directory with subdirectories named as serial numbers to stro results
    	    		this.focusingHistoryFile,
    	    		this.useLMAMetrics, // measure/report focal distance and tilts using lens model/LMA (when available)
    	    		this.EEPROM_channel,// EEPROM channel to read serial number from
    	    		this.saveResults, // save focusing results
    	    		this.showResults, // show focusing (includingh intermediate) results
    	    		this.serialNumber, // camera serial number string
    	    		this.sensorTemperature, // last measured sensor temperature
    	    		this.result_lastKT,   // focal distance temperature coefficient (um/C), measured from last run 
    	    		this.result_lastFD20, // focal distance for 20C, measured from last run 
    	    		this.result_allHistoryKT,   // focal distance temperature coefficient (um/C), measured from all the measurement histgory 
    	    		this.result_allHistoryFD20, // focal distance for 20C, measured from  all the measurement histgory
    	    		this.result_fDistance, // last measured focal distance
    	    		this.result_tiltX, // last measured tilt X
    	    		this.result_tiltY, // last measured tilt Y
    	    		this.result_R50,   // last measured R50 (average PSF 50% level radius, pixels - somewhat larged than actual because of measurement settings)
    	    		this.result_A50,   // last measured A50 (simailar, but R^2 are averaged) 
    	    		this.result_B50,   // last measured B50 (simailar, but R^4 are averaged)
    	    		this.result_RC50,  // last measured RC50(R50 calculated only for the 2 center samples)
    	    		this.result_PX0, // lens center shift, X
    	    		this.result_PY0, // lens center shift, Y
    	    		this.result_PSI, // SFE rotation (from grid)
    	    		this.result_ROT, // SFE rotation (from head lasers)
    	    		this.result_FocalLength,
    	    		this.comment, // Comment to add to the results
    	    		this.lensSerial, // Lens serial number
    	    		this.manufacturingState,
    	    		this.askLensSerial,
    	    		this.includeLensSerial, // add lens serial to config filename
    	        	this.centerDeltaX, // required X-difference between lens center and sensor center 
    	        	this.centerDeltaY, // required Y-difference between lens center and sensor center
    	    		this.margins,
        			this.numSamples,
        			this.sampleSize,
        			this.numInCenter,// 2 - number of "center" samples
        			this.centerSamples,
        			this.maxCorr,
    				this.showHistoryDetails,
    				this.showHistorySamples,
    				this.showHistorySingleLine, // all parameters in a single line (easier to copy to spreadsheet)
    				this.showAcquiredImages,
    				this.showROI,
    				this.showSamples,
    				this.showFittedParameters,
    				this.useHeadLasers,
    				this.psf_cutoffEnergy,
    				this.psf_cutoffLevel,
    				this.psf_minArea,
    				this.psf_blurSigma,
    				this.weightRatioRedToGreen,  // Use this data when combining defocusing data from different color PSF
    				this.weightRatioBlueToGreen,
    				this.targetFarNear,          // target logariphm of average tangential-to-radial resolution
    				this.useRadialTangential,      // Use targetFarNear (radial/tangential resolution)  as a proxy for the distance  
    				this.targetMicrons,          // target lens center distance (away from "best focus"
    				this.toleranceMicrons, // microns
    				this.toleranceTilt, // 
    				this.toleranceThreshold, // When each error is under swcaled thereshold, reduce correxction step twice
    				this.parallelAdjustThreshold, // adjust 3 motors parallel if focal distance error in the center exceeds this    
    				
    				this.motorsSigma,            // when fitting planes for far/near, tiltX and tiltY the weights of the samples decay with this sigma
    				this.motorsSigma3,           // same when 3 motors move together
    				this.motorsMinSigma,           // sigma will not drop below this value when fitting walk is getting smaller 
    				this.motorsVarSigmaToTravel,   // when walk is getting smaller, sigma will keep going down proportionally
    				this.motorsFadeSigma,          // after each step new sigma will have this part of the calculated from the travel 
    				this.motorsOverShootToBalance, //For quadratic maximum the correction will be increased by 1+motorsOverShootToBalance if there are less samples on teh other side
    				this.filterGoodDistance,     // when measuring tilt, use those with good center with higher weight
    				this.goodDistanceSigma,      // sigma for the weight function of tilt measurements, depending on the center distance error
    				this.goodTiltSigma,            // weight decay for heavily tilted samples
    				this.maxStep,                // maximal allowed single-step focusing adjustment
    				this.probeStep,              // how far to go to probe around the current point to measure derivatives
    				this.probe_M1M2M3,     // how far to move average of the 3 motors: (M1+M2+M3)/3
    				this.probe_M3_M1M2,    // how far to move M3 opposite to M1 and M2: M3-(M1+M2)/2
    				this.probe_M2_M1,      // how far to move M2 opposite to M1:    M2-M1
    				this.sigmaToProbe,     // data from far samples decay proportionally to the probe distances 
    				this.useTheBest,      // adjust from the best known position (false - from the last)  	                
    				this.probeSymmetrical,       // if true, probe 6 measurements), if false - only 4  (tetrahedron)
    				this.parallelBeforeProbing, // move 3 motors before probing around
    				this.reProbeDistance,        // re-run probing around in orthoganal directions, if the current position move farther from the last probing one
    				this.believeLast,            // coefficient 0.. 1.0. When each ot the 3 parameters is linearized, add shift, so 1.0 the planes will go throug the last sample
    				this.compensateHysteresis,   // move motors in the same direction to compensate fro hysteresis
    				this.minCorr,                // minimal correction movement to initiate final numFinalCorr moves
    				this.numFinalCorr,           // exit if this number of last corrections where below  minCorr
    				this.minCorrPre,                // minimal correction movement to initiate final numFinalCorr moves
    				this.numFinalCorrPre,           // exit if this number of last corrections where below  minCorr
    				this.maxAutoIterations,      // exit if history grows above this
    				this.maxAutoDistance,        // Maximal allowed automatic correction, motor steps
    				this.confirmFirstAuto,       // ask confirmation after first automatic adjustment step (before moving)
    				this.motorsPreSigma,         // when fitting parabola for focusing sharpness in the center, far measurements decay with this sigma
    				this.maxLinearStep,          // If there are insufficient measurements to fit parabola - make this step
    				this.scanStep,               // motor step (all 3 motors) in scan focus mode (signed value)
    				this.scanNumber,             // number of scanStep steps to run
    	            this.scanNumberNegative,    // of them negative
    				this.scanHysteresis,         // scan both ways
    				this.scanHysteresisNumber,   // number of test points for the Hysteresis measurement
    				
    				this.scanTiltEnable,  // enable scanning tilt
    				this.scanTiltReverse,
    				this.scanMeasureLast,
    				this.scanRunLMA,
    	    		this.scanTiltRangeX,    // 4 periods
    	    		this.scanTiltRangeY,    // 4 periods
    	    		this.scanTiltStepsX,
    	    		this.scanTiltStepsY,
    				this.motorHysteresis,
    				
    				this.measuredHysteresis,     // actually measured (will be saved/restored)
    				this.motorCalm,              // wait (seconds) after motors reached final position (for the first time) before acquiring image   
    				this.linearReductionRatio,   // sensor travel to motors travel (all 3 together), By design it is 4/38~=0.105
    				this.motorDebug,             // 1 show  motor moves, 2 - show hysteresis back-ups too
    				this.lensDistanceNumPoints,  // number of points to tabulate center focus parameters vs. focal distance 
    				this.lensDistancePolynomialDegree, // polynomial degree to approximate center focus parameters vs. focal distance 
    				this.lensDistanceWeightY,    // normalize overall sharpness (that depends on the lens quality and/or PSF parameters) to differential ones 
    				this.lensDistanceWeightK,    // 0.0 use 3 distances (to frac-R, frac-B and Y) with teh same weight when fitting, 1.0 - proportional to squared derivatives
    				this.lensDistanceInteractive, // Open dialog when calibrating focal distance
    				this.lensDistanceShowResults, // show results window from foca
    				this.lensDistanceMoveToGoal,  // Move to targetMicrons
    				this.powerControlEnable,
    				this.powerControlMaximalTemperature,
    				this.powerControlHeaterOnMinutes,
    				this.powerControlNeitherOnMinutes,
    				this.powerControlFanOnMinutes,
                    this.uvLasersIP,              // IP address of the camera with UV LEDs and aiming lasers are connected
                    this.uvLasersBus,             // 0 if 103641 board is connected to the sensor port (through 10-359), 1 - to 10369
                    this.uvLasersCurrents,        // default LED on currents (mA)

    				this.smallestSubPix,
    				this.bitmapNonuniforityThreshold,
    				this.subdiv, 
    				this.overexposedMaxFraction, 
    				this.minDefinedArea,
    				this.PSFKernelSize,
    				this.approximateGrid,
        			this.centerPSF,
        			this.mask1_sigma,
        			this.mask1_threshold,
        			this.gaps_sigma,
        			this.mask_denoise,
        			this.deconvInvert,
    				this.correlationSize,
    				this.correlationGaussWidth,
    				this.minUVSpan,
    				this.correlationMinContrast,
    				this.correlationMinInitialContrast,
    				this.correlationMinAbsoluteContrast,
    				this.correlationMinAbsoluteInitialContrast,
    				this.postUVscrewSensitivity,
    				this.flatFieldCorrection,
    				this.flatFieldExpand,
    				this.thresholdFinish, 
    				this.numIterations,
    				this.cameraIsConfigured,
    				this.useExtraSensor,
    				this.enRoundOff,
    				this.keepCircularMask,
    				this.motorPos,
    	        	this.ampsSeconds, // cumulative Amps*seconds (read only, but will be saved/restored)
    	        	this.reportTemperature,
    	        	this.showLegacy
        			);
    	}
		public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"gridGeometryFile",this.gridGeometryFile+"");
			properties.setProperty(prefix+"initialCalibrationFile",this.initialCalibrationFile+"");
			properties.setProperty(prefix+"strategyFile",this.strategyFile+"");
			properties.setProperty(prefix+"resultsSuperDirectory",this.resultsSuperDirectory+"");
			properties.setProperty(prefix+"focusingHistoryFile",this.focusingHistoryFile+"");
			properties.setProperty(prefix+"useLMAMetrics",this.useLMAMetrics+"");
			properties.setProperty(prefix+"serialNumber",this.serialNumber+"");
			if (!Double.isNaN(this.sensorTemperature))properties.setProperty(prefix+"sensorTemperature",this.sensorTemperature+"");
			if (!Double.isNaN(this.result_lastKT))properties.setProperty(prefix+"result_lastKT",this.result_lastKT+"");
			if (!Double.isNaN(this.result_lastFD20))properties.setProperty(prefix+"result_lastFD20",this.result_lastFD20+"");
			if (!Double.isNaN(this.result_allHistoryKT))properties.setProperty(prefix+"result_allHistoryKT",this.result_allHistoryKT+"");
			if (!Double.isNaN(this.result_allHistoryFD20))properties.setProperty(prefix+"result_allHistoryFD20",this.result_allHistoryFD20+"");
			if (!Double.isNaN(this.result_fDistance))properties.setProperty(prefix+"result_fDistance",this.result_fDistance+"");
			if (!Double.isNaN(this.result_tiltX))properties.setProperty(prefix+"result_tiltX",this.result_tiltX+"");
			if (!Double.isNaN(this.result_tiltY))properties.setProperty(prefix+"result_tiltY",this.result_tiltY+"");
			if (!Double.isNaN(this.result_R50))properties.setProperty(prefix+"result_R50",this.result_R50+"");
			if (!Double.isNaN(this.result_A50))properties.setProperty(prefix+"result_A50",this.result_A50+"");
			if (!Double.isNaN(this.result_B50))properties.setProperty(prefix+"result_B50",this.result_B50+"");
			if (!Double.isNaN(this.result_RC50))properties.setProperty(prefix+"result_RC50",this.result_RC50+"");
			
			if (!Double.isNaN(this.result_PX0))properties.setProperty(prefix+"result_PX0",this.result_PX0+"");
			if (!Double.isNaN(this.result_PY0))properties.setProperty(prefix+"result_PY0",this.result_PY0+"");
			if (!Double.isNaN(this.result_PSI))properties.setProperty(prefix+"result_PSI",this.result_PSI+"");
			if (!Double.isNaN(this.result_ROT))properties.setProperty(prefix+"result_ROT",this.result_ROT+"");
			if (!Double.isNaN(this.result_FocalLength))properties.setProperty(prefix+"result_FocalLength",this.result_FocalLength+"");
			properties.setProperty(prefix+"EEPROM_channel",this.EEPROM_channel+"");
			properties.setProperty(prefix+"saveResults",this.saveResults+"");
			properties.setProperty(prefix+"showResults",this.showResults+"");
			properties.setProperty(prefix+"comment","<![CDATA["+this.comment+ "]]>");
			properties.setProperty(prefix+"lensSerial",this.lensSerial);
			properties.setProperty(prefix+"manufacturingState",this.manufacturingState+"");
			properties.setProperty(prefix+"askLensSerial",this.askLensSerial+"");
			properties.setProperty(prefix+"includeLensSerial",this.includeLensSerial+"");
			properties.setProperty(prefix+"centerDeltaX",this.centerDeltaX+"");
			properties.setProperty(prefix+"centerDeltaY",this.centerDeltaY+"");
			properties.setProperty(prefix+"margins_x",this.margins.x+"");
			properties.setProperty(prefix+"margins_y",this.margins.y+"");
			properties.setProperty(prefix+"margins_width",this.margins.width+"");
			properties.setProperty(prefix+"margins_height",this.margins.height+"");
			properties.setProperty(prefix+"numSamples_0",this.numSamples[0]+"");
			properties.setProperty(prefix+"numSamples_1",this.numSamples[1]+"");
			properties.setProperty(prefix+"sampleSize",this.sampleSize+"");
			properties.setProperty(prefix+"numInCenter",this.numInCenter+"");
			properties.setProperty(prefix+"centerSamples",this.centerSamples+"");
			properties.setProperty(prefix+"maxCorr",this.maxCorr+"");
			properties.setProperty(prefix+"showHistoryDetails",this.showHistoryDetails+"");
			properties.setProperty(prefix+"showHistorySamples",this.showHistorySamples+"");
			
			properties.setProperty(prefix+"showHistorySingleLine",this.showHistorySingleLine+"");
			properties.setProperty(prefix+"showAcquiredImages",this.showAcquiredImages+"");
			properties.setProperty(prefix+"showROI",this.showROI+"");
			properties.setProperty(prefix+"showSamples",this.showSamples+"");
			properties.setProperty(prefix+"showFittedParameters",this.showFittedParameters+"");
			properties.setProperty(prefix+"useHeadLasers",this.useHeadLasers+"");
			properties.setProperty(prefix+"useExtraSensor",this.useExtraSensor+"");
			properties.setProperty(prefix+"enRoundOff",this.enRoundOff+"");
			properties.setProperty(prefix+"keepCircularMask",this.keepCircularMask+"");
			properties.setProperty(prefix+"psf_cutoffEnergy",this.psf_cutoffEnergy+"");
			properties.setProperty(prefix+"psf_cutoffLevel",this.psf_cutoffLevel+"");
			properties.setProperty(prefix+"psf_minArea",this.psf_minArea+"");
			properties.setProperty(prefix+"psf_blurSigma",this.psf_blurSigma+"");
			properties.setProperty(prefix+"weightRatioRedToGreen",this.weightRatioRedToGreen+"");
			properties.setProperty(prefix+"weightRatioBlueToGreen",this.weightRatioBlueToGreen+"");
			properties.setProperty(prefix+"targetFarNear",this.targetFarNear+"");
			properties.setProperty(prefix+"useRadialTangential",this.useRadialTangential+"");
			properties.setProperty(prefix+"targetMicrons",this.targetMicrons+"");
			properties.setProperty(prefix+"toleranceMicrons",this.toleranceMicrons+"");
			properties.setProperty(prefix+"toleranceTilt",this.toleranceTilt+"");
			properties.setProperty(prefix+"toleranceThreshold",this.toleranceThreshold+"");
			properties.setProperty(prefix+"parallelAdjustThreshold",this.parallelAdjustThreshold+"");
			properties.setProperty(prefix+"motorsSigma",this.motorsSigma+"");
			properties.setProperty(prefix+"motorsSigma3",this.motorsSigma3+"");
			properties.setProperty(prefix+"motorsMinSigma",this.motorsMinSigma+"");
			properties.setProperty(prefix+"motorsVarSigmaToTravel",this.motorsVarSigmaToTravel+"");
			properties.setProperty(prefix+"motorsFadeSigma",this.motorsFadeSigma+"");
			properties.setProperty(prefix+"motorsOverShootToBalance",this.motorsOverShootToBalance+"");
			
			properties.setProperty(prefix+"filterGoodDistance",this.filterGoodDistance+"");
			properties.setProperty(prefix+"goodDistanceSigma",this.goodDistanceSigma+"");
			properties.setProperty(prefix+"goodTiltSigma",this.goodTiltSigma+"");
			properties.setProperty(prefix+"maxStep",this.maxStep+"");
			properties.setProperty(prefix+"probeStep",this.probeStep+"");
			
			properties.setProperty(prefix+"probe_M1M2M3",this.probe_M1M2M3+"");
			properties.setProperty(prefix+"probe_M3_M1M2",this.probe_M3_M1M2+"");
			properties.setProperty(prefix+"probe_M2_M1",this.probe_M2_M1+"");
			properties.setProperty(prefix+"sigmaToProbe",this.sigmaToProbe+"");
			properties.setProperty(prefix+"useTheBest",this.useTheBest+"");
			properties.setProperty(prefix+"probeSymmetrical",this.probeSymmetrical+"");
			properties.setProperty(prefix+"parallelBeforeProbing",this.parallelBeforeProbing+"");
			properties.setProperty(prefix+"reProbeDistance",this.reProbeDistance+"");
			properties.setProperty(prefix+"believeLast",this.believeLast+"");
			properties.setProperty(prefix+"compensateHysteresis",this.compensateHysteresis+"");
			properties.setProperty(prefix+"minCorr",this.minCorr+"");
			properties.setProperty(prefix+"numFinalCorr",this.numFinalCorr+"");
			properties.setProperty(prefix+"minCorrPre",this.minCorrPre+"");
			properties.setProperty(prefix+"numFinalCorrPre",this.numFinalCorrPre+"");
			properties.setProperty(prefix+"maxAutoIterations",this.maxAutoIterations+"");
			properties.setProperty(prefix+"maxAutoDistance",this.maxAutoDistance+"");
			properties.setProperty(prefix+"confirmFirstAuto",this.confirmFirstAuto+"");
			properties.setProperty(prefix+"motorsPreSigma",this.motorsPreSigma+"");
			properties.setProperty(prefix+"maxLinearStep",this.maxLinearStep+"");
			properties.setProperty(prefix+"scanStep",this.scanStep+"");
			properties.setProperty(prefix+"scanNumber",this.scanNumber+"");
			properties.setProperty(prefix+"scanNumberNegative",this.scanNumberNegative+"");
			properties.setProperty(prefix+"scanHysteresis",this.scanHysteresis+"");
			properties.setProperty(prefix+"scanHysteresisNumber",this.scanHysteresisNumber+"");
			properties.setProperty(prefix+"scanTiltEnable",this.scanTiltEnable+"");  // enable scanning tilt
			properties.setProperty(prefix+"scanTiltReverse",this.scanTiltReverse+"");
			properties.setProperty(prefix+"scanMeasureLast",this.scanMeasureLast+"");
			properties.setProperty(prefix+"scanRunLMA",this.scanRunLMA+"");
			properties.setProperty(prefix+"scanTiltRangeX",this.scanTiltRangeX+"");    // 4 periods
			properties.setProperty(prefix+"scanTiltRangeY",this.scanTiltRangeY+"");    // 4 periods
			properties.setProperty(prefix+"scanTiltStepsX",this.scanTiltStepsX+"");
			properties.setProperty(prefix+"scanTiltStepsY",this.scanTiltStepsY+"");
			properties.setProperty(prefix+"motorHysteresis",this.motorHysteresis+"");
			properties.setProperty(prefix+"measuredHysteresis",this.measuredHysteresis+"");
			properties.setProperty(prefix+"motorCalm",this.motorCalm+"");
			properties.setProperty(prefix+"linearReductionRatio",this.linearReductionRatio+"");
			properties.setProperty(prefix+"motorDebug",this.motorDebug+"");
			properties.setProperty(prefix+"lensDistanceNumPoints",this.lensDistanceNumPoints+"");
			properties.setProperty(prefix+"lensDistancePolynomialDegree",this.lensDistancePolynomialDegree+"");
			properties.setProperty(prefix+"lensDistanceWeightY",this.lensDistanceWeightY+"");
			properties.setProperty(prefix+"lensDistanceWeightK",this.lensDistanceWeightK+"");
			properties.setProperty(prefix+"lensDistanceInteractive",this.lensDistanceInteractive+"");
			properties.setProperty(prefix+"lensDistanceShowResults",this.lensDistanceShowResults+"");
			properties.setProperty(prefix+"lensDistanceMoveToGoal",this.lensDistanceMoveToGoal+"");
			
			properties.setProperty(prefix+"powerControlEnable",this.powerControlEnable+"");
			properties.setProperty(prefix+"powerControlMaximalTemperature",this.powerControlMaximalTemperature+"");
			properties.setProperty(prefix+"powerControlHeaterOnMinutes",this.powerControlHeaterOnMinutes+"");
			properties.setProperty(prefix+"powerControlNeitherOnMinutes",this.powerControlNeitherOnMinutes+"");
			properties.setProperty(prefix+"powerControlFanOnMinutes",this.powerControlFanOnMinutes+"");
			
			properties.setProperty(prefix+"uvLasersIP",this.uvLasersIP);
			properties.setProperty(prefix+"uvLasersBus",this.uvLasersBus+"");
			properties.setProperty(prefix+"uvLasersCurrents_0",this.uvLasersCurrents[0]+"");
			properties.setProperty(prefix+"uvLasersCurrents_1",this.uvLasersCurrents[1]+"");
			properties.setProperty(prefix+"uvLasersCurrents_2",this.uvLasersCurrents[2]+"");
			properties.setProperty(prefix+"uvLasersCurrents_3",this.uvLasersCurrents[3]+"");
			
			properties.setProperty(prefix+"smallestSubPix",this.smallestSubPix+"");
			properties.setProperty(prefix+"bitmapNonuniforityThreshold",this.bitmapNonuniforityThreshold+"");
			properties.setProperty(prefix+"subdiv",this.subdiv+"");
			properties.setProperty(prefix+"overexposedMaxFraction",this.overexposedMaxFraction+"");
			properties.setProperty(prefix+"minDefinedArea",this.minDefinedArea+"");
			properties.setProperty(prefix+"PSFKernelSize",this.PSFKernelSize+"");
			properties.setProperty(prefix+"approximateGrid",this.approximateGrid+"");
			properties.setProperty(prefix+"centerPSF",this.centerPSF+"");
			properties.setProperty(prefix+"mask1_sigma",this.mask1_sigma+"");
			properties.setProperty(prefix+"mask1_threshold",this.mask1_threshold+"");
			properties.setProperty(prefix+"gaps_sigma",this.gaps_sigma+"");
			properties.setProperty(prefix+"mask_denoise",this.mask_denoise+"");
			properties.setProperty(prefix+"deconvInvert",this.deconvInvert+"");
			properties.setProperty(prefix+"correlationSize",this.correlationSize+"");
			properties.setProperty(prefix+"correlationGaussWidth",this.correlationGaussWidth+"");
			properties.setProperty(prefix+"minUVSpan",this.minUVSpan+"");
			
			properties.setProperty(prefix+"correlationMinContrast",this.correlationMinContrast+"");
			properties.setProperty(prefix+"correlationMinInitialContrast",this.correlationMinInitialContrast+"");
			properties.setProperty(prefix+"correlationMinAbsoluteContrast",this.correlationMinAbsoluteContrast+"");
			properties.setProperty(prefix+"correlationMinAbsoluteInitialContrast",this.correlationMinAbsoluteInitialContrast+"");
			for (int i=0;i<this.postUVscrewSensitivity.length;i++){
				properties.setProperty(prefix+"postUVscrewSensitivity_"+i,this.postUVscrewSensitivity[i]+"");
			}
			properties.setProperty(prefix+"flatFieldCorrection",this.flatFieldCorrection+"");
			properties.setProperty(prefix+"flatFieldExpand",this.flatFieldExpand+"");
			properties.setProperty(prefix+"thresholdFinish",this.thresholdFinish+"");
			properties.setProperty(prefix+"numIterations",this.numIterations+"");
			for (int i=0;i<this.ampsSeconds.length;i++) 
				properties.setProperty(prefix+"ampsSeconds_"+i,this.ampsSeconds[i]+"");
			properties.setProperty(prefix+"reportTemperature",this.reportTemperature+"");
			properties.setProperty(prefix+"showLegacy",this.showLegacy+"");
		}    	
		public void getProperties(String prefix,Properties properties){
			if (properties.getProperty(prefix+"gridGeometryFile")!=null)
				this.gridGeometryFile=properties.getProperty(prefix+"gridGeometryFile");
			if (properties.getProperty(prefix+"initialCalibrationFile")!=null)
				this.initialCalibrationFile=properties.getProperty(prefix+"initialCalibrationFile");
			if (properties.getProperty(prefix+"strategyFile")!=null)
				this.strategyFile=properties.getProperty(prefix+"strategyFile");
			if (properties.getProperty(prefix+"resultsSuperDirectory")!=null)
				this.resultsSuperDirectory=properties.getProperty(prefix+"resultsSuperDirectory");
			if (properties.getProperty(prefix+"focusingHistoryFile")!=null)
				this.focusingHistoryFile=properties.getProperty(prefix+"focusingHistoryFile");

			if (properties.getProperty(prefix+"useLMAMetrics")!=null)
				this.useLMAMetrics=Boolean.parseBoolean(properties.getProperty(prefix+"useLMAMetrics"));
			
			if (properties.getProperty(prefix+"serialNumber")!=null)
				this.serialNumber=properties.getProperty(prefix+"serialNumber");
			//	this.serialNumber is only written, but never read from the configuration file (only from devivce)
			
			if (properties.getProperty(prefix+"sensorTemperature")!=null) this.sensorTemperature=Double.parseDouble(properties.getProperty(prefix+"sensorTemperature"));
			else this.sensorTemperature=Double.NaN;
			if (properties.getProperty(prefix+"result_lastKT")!=null) this.result_lastKT=Double.parseDouble(properties.getProperty(prefix+"result_lastKT"));
			else this.result_lastKT=Double.NaN;
			if (properties.getProperty(prefix+"result_lastFD20")!=null) this.result_lastFD20=Double.parseDouble(properties.getProperty(prefix+"result_lastFD20"));
			else this.result_lastFD20=Double.NaN;
			if (properties.getProperty(prefix+"result_allHistoryKT")!=null) this.result_allHistoryKT=Double.parseDouble(properties.getProperty(prefix+"result_allHistoryKT"));
			else this.result_allHistoryKT=Double.NaN;
			if (properties.getProperty(prefix+"result_allHistoryFD20")!=null) this.result_allHistoryFD20=Double.parseDouble(properties.getProperty(prefix+"result_allHistoryFD20"));
			else this.result_allHistoryFD20=Double.NaN;
			if (properties.getProperty(prefix+"result_fDistance")!=null) this.result_fDistance=Double.parseDouble(properties.getProperty(prefix+"result_fDistance"));
			else this.result_fDistance=Double.NaN;
			if (properties.getProperty(prefix+"result_tiltX")!=null) this.result_tiltX=Double.parseDouble(properties.getProperty(prefix+"result_tiltX"));
			else this.result_tiltX=Double.NaN;
			if (properties.getProperty(prefix+"result_tiltY")!=null) this.result_tiltY=Double.parseDouble(properties.getProperty(prefix+"result_tiltY"));
			else this.result_tiltY=Double.NaN;
			if (properties.getProperty(prefix+"result_R50")!=null) this.result_R50=Double.parseDouble(properties.getProperty(prefix+"result_R50"));
			else this.result_R50=Double.NaN;
			if (properties.getProperty(prefix+"result_A50")!=null) this.result_A50=Double.parseDouble(properties.getProperty(prefix+"result_A50"));
			else this.result_A50=Double.NaN;
			if (properties.getProperty(prefix+"result_B50")!=null) this.result_B50=Double.parseDouble(properties.getProperty(prefix+"result_B50"));
			else this.result_B50=Double.NaN;
			if (properties.getProperty(prefix+"result_RC50")!=null) this.result_RC50=Double.parseDouble(properties.getProperty(prefix+"result_RC50"));
			else this.result_RC50=Double.NaN;
			if (properties.getProperty(prefix+"result_PX0")!=null) this.result_PX0=Double.parseDouble(properties.getProperty(prefix+"result_PX0"));
			else this.result_PX0=Double.NaN;
			if (properties.getProperty(prefix+"result_PY0")!=null) this.result_PY0=Double.parseDouble(properties.getProperty(prefix+"result_PY0"));
			else this.result_PY0=Double.NaN;
			if (properties.getProperty(prefix+"result_PSI")!=null) this.result_PSI=Double.parseDouble(properties.getProperty(prefix+"result_PSI"));
			else this.result_PSI=Double.NaN;
			if (properties.getProperty(prefix+"result_ROT")!=null) this.result_ROT=Double.parseDouble(properties.getProperty(prefix+"result_ROT"));
			else this.result_ROT=Double.NaN;
			if (properties.getProperty(prefix+"result_FocalLength")!=null) this.result_FocalLength=Double.parseDouble(properties.getProperty(prefix+"result_FocalLength"));
			else this.result_FocalLength=Double.NaN;
			if (properties.getProperty(prefix+"EEPROM_channel")!=null)
				this.EEPROM_channel=Integer.parseInt(properties.getProperty(prefix+"EEPROM_channel"));
			if (properties.getProperty(prefix+"saveResults")!=null)
				this.saveResults=Boolean.parseBoolean(properties.getProperty(prefix+"saveResults"));
			if (properties.getProperty(prefix+"showResults")!=null)
				this.showResults=Boolean.parseBoolean(properties.getProperty(prefix+"showResults"));
			if (properties.getProperty(prefix+"comment")!=null)
				this.comment=properties.getProperty(prefix+"comment");
			if ((this.comment.length()>10) && this.comment.substring(0,9).equals("<![CDATA[")) this.comment=this.comment.substring(9,this.comment.length()-3);
			if (properties.getProperty(prefix+"lensSerial")!=null)
				this.lensSerial= properties.getProperty(prefix+"lensSerial");
			if (properties.getProperty(prefix+"manufacturingState")!=null)
				this.manufacturingState=Integer.parseInt(properties.getProperty(prefix+"manufacturingState"));
			if (properties.getProperty(prefix+"askLensSerial")!=null)
				this.askLensSerial=Boolean.parseBoolean(properties.getProperty(prefix+"askLensSerial"));
			if (properties.getProperty(prefix+"includeLensSerial")!=null)
				this.includeLensSerial=Boolean.parseBoolean(properties.getProperty(prefix+"includeLensSerial"));
			if (properties.getProperty(prefix+"centerDeltaX")!=null)
				this.centerDeltaX=Double.parseDouble(properties.getProperty(prefix+"centerDeltaX"));
			if (properties.getProperty(prefix+"centerDeltaY")!=null)
				this.centerDeltaY=Double.parseDouble(properties.getProperty(prefix+"centerDeltaY"));
			if (properties.getProperty(prefix+"margins_x")!=null)
				this.margins.x=Integer.parseInt(properties.getProperty(prefix+"margins_x"));
			if (properties.getProperty(prefix+"margins_y")!=null)
				this.margins.y=Integer.parseInt(properties.getProperty(prefix+"margins_y"));
			if (properties.getProperty(prefix+"margins_width")!=null)
				this.margins.width=Integer.parseInt(properties.getProperty(prefix+"margins_width"));
			if (properties.getProperty(prefix+"margins_height")!=null)
				this.margins.height=Integer.parseInt(properties.getProperty(prefix+"margins_height"));
			if (properties.getProperty(prefix+"numSamples_0")!=null)
				this.numSamples[0]=Integer.parseInt(properties.getProperty(prefix+"numSamples_0"));
			if (properties.getProperty(prefix+"numSamples_1")!=null)
				this.numSamples[1]=Integer.parseInt(properties.getProperty(prefix+"numSamples_1"));
			if (properties.getProperty(prefix+"sampleSize")!=null)
				this.sampleSize=Integer.parseInt(properties.getProperty(prefix+"sampleSize"));
			if (properties.getProperty(prefix+"numInCenter")!=null)
				this.numInCenter=Integer.parseInt(properties.getProperty(prefix+"numInCenter"));
			if (properties.getProperty(prefix+"centerSamples")!=null)
				this.centerSamples=Boolean.parseBoolean(properties.getProperty(prefix+"centerSamples"));
			if (properties.getProperty(prefix+"maxCorr")!=null)
				this.maxCorr=Double.parseDouble(properties.getProperty(prefix+"maxCorr"));
			if (properties.getProperty(prefix+"showHistoryDetails")!=null)
				this.showHistoryDetails=Boolean.parseBoolean(properties.getProperty(prefix+"showHistoryDetails"));
			if (properties.getProperty(prefix+"showHistorySamples")!=null)
				this.showHistorySamples=Boolean.parseBoolean(properties.getProperty(prefix+"showHistorySamples"));
			if (properties.getProperty(prefix+"showHistorySingleLine")!=null)
				this.showHistorySingleLine=Boolean.parseBoolean(properties.getProperty(prefix+"showHistorySingleLine"));
			if (properties.getProperty(prefix+"showAcquiredImages")!=null)
				this.showAcquiredImages=Boolean.parseBoolean(properties.getProperty(prefix+"showAcquiredImages"));
			if (properties.getProperty(prefix+"showROI")!=null)
				this.showROI=Boolean.parseBoolean(properties.getProperty(prefix+"showROI"));
			if (properties.getProperty(prefix+"showSamples")!=null)
				this.showSamples=Boolean.parseBoolean(properties.getProperty(prefix+"showSamples"));
			if (properties.getProperty(prefix+"showFittedParameters")!=null)
				this.showFittedParameters=Boolean.parseBoolean(properties.getProperty(prefix+"showFittedParameters"));
			if (properties.getProperty(prefix+"useHeadLasers")!=null)
				this.useHeadLasers=Boolean.parseBoolean(properties.getProperty(prefix+"useHeadLasers"));
			if (properties.getProperty(prefix+"useExtraSensor")!=null)
				this.useExtraSensor=Boolean.parseBoolean(properties.getProperty(prefix+"useExtraSensor"));
			if (properties.getProperty(prefix+"enRoundOff")!=null)
				this.enRoundOff = Boolean.parseBoolean(properties.getProperty(prefix+"enRoundOff"));
			
			if (properties.getProperty(prefix+"keepCircularMask")!=null)
				this.keepCircularMask=Boolean.parseBoolean(properties.getProperty(prefix+"keepCircularMask"));
			if (properties.getProperty(prefix+"psf_cutoffEnergy")!=null)
				this.psf_cutoffEnergy=Double.parseDouble(properties.getProperty(prefix+"psf_cutoffEnergy"));
			if (properties.getProperty(prefix+"psf_cutoffLevel")!=null)
				this.psf_cutoffLevel=Double.parseDouble(properties.getProperty(prefix+"psf_cutoffLevel"));
			if (properties.getProperty(prefix+"psf_minArea")!=null)
				this.psf_minArea=Integer.parseInt(properties.getProperty(prefix+"psf_minArea"));
			if (properties.getProperty(prefix+"psf_blurSigma")!=null)
				this.psf_blurSigma=Double.parseDouble(properties.getProperty(prefix+"psf_blurSigma"));
			if (properties.getProperty(prefix+"weightRatioRedToGreen")!=null)
				this.weightRatioRedToGreen=Double.parseDouble(properties.getProperty(prefix+"weightRatioRedToGreen"));
			if (properties.getProperty(prefix+"weightRatioBlueToGreen")!=null)
				this.weightRatioBlueToGreen=Double.parseDouble(properties.getProperty(prefix+"weightRatioBlueToGreen"));
			if (properties.getProperty(prefix+"targetFarNear")!=null)
				this.targetFarNear=Double.parseDouble(properties.getProperty(prefix+"targetFarNear"));
			if (properties.getProperty(prefix+"useRadialTangential")!=null)
				this.useRadialTangential=Boolean.parseBoolean(properties.getProperty(prefix+"useRadialTangential"));
			if (properties.getProperty(prefix+"targetMicrons")!=null)
				this.targetMicrons=Double.parseDouble(properties.getProperty(prefix+"targetMicrons"));
			if (properties.getProperty(prefix+"toleranceMicrons")!=null)
				this.toleranceMicrons=Double.parseDouble(properties.getProperty(prefix+"toleranceMicrons"));
			if (properties.getProperty(prefix+"toleranceTilt")!=null)
				this.toleranceTilt=Double.parseDouble(properties.getProperty(prefix+"toleranceTilt"));
			if (properties.getProperty(prefix+"toleranceThreshold")!=null)
				this.toleranceThreshold=Double.parseDouble(properties.getProperty(prefix+"toleranceThreshold"));
			if (properties.getProperty(prefix+"parallelAdjustThreshold")!=null)
				this.parallelAdjustThreshold=Double.parseDouble(properties.getProperty(prefix+"parallelAdjustThreshold"));
			if (properties.getProperty(prefix+"motorsSigma")!=null)
				this.motorsSigma=Double.parseDouble(properties.getProperty(prefix+"motorsSigma"));
			if (properties.getProperty(prefix+"motorsSigma3")!=null)
				this.motorsSigma3=Double.parseDouble(properties.getProperty(prefix+"motorsSigma3"));
			if (properties.getProperty(prefix+"motorsMinSigma")!=null)
				this.motorsMinSigma=Double.parseDouble(properties.getProperty(prefix+"motorsMinSigma"));
			if (properties.getProperty(prefix+"motorsVarSigmaToTravel")!=null)
				this.motorsVarSigmaToTravel=Double.parseDouble(properties.getProperty(prefix+"motorsVarSigmaToTravel"));
			if (properties.getProperty(prefix+"motorsFadeSigma")!=null)
				this.motorsFadeSigma=Double.parseDouble(properties.getProperty(prefix+"motorsFadeSigma"));
			if (properties.getProperty(prefix+"motorsOverShootToBalance")!=null)
				this.motorsOverShootToBalance=Double.parseDouble(properties.getProperty(prefix+"motorsOverShootToBalance"));
			
			if (properties.getProperty(prefix+"filterGoodDistance")!=null)
				this.filterGoodDistance=Boolean.parseBoolean(properties.getProperty(prefix+"filterGoodDistance"));
			if (properties.getProperty(prefix+"goodDistanceSigma")!=null)
				this.goodDistanceSigma=Double.parseDouble(properties.getProperty(prefix+"goodDistanceSigma"));
			
			if (properties.getProperty(prefix+"goodTiltSigma")!=null)
				this.goodTiltSigma=Double.parseDouble(properties.getProperty(prefix+"goodTiltSigma"));
			
			if (properties.getProperty(prefix+"maxStep")!=null)
				this.maxStep=Double.parseDouble(properties.getProperty(prefix+"maxStep"));
			if (properties.getProperty(prefix+"probeStep")!=null)
				this.probeStep=Double.parseDouble(properties.getProperty(prefix+"probeStep"));
			
			if (properties.getProperty(prefix+"probe_M1M2M3")!=null)
				this.probe_M1M2M3=Double.parseDouble(properties.getProperty(prefix+"probe_M1M2M3"));
			if (properties.getProperty(prefix+"probe_M3_M1M2")!=null)
				this.probe_M3_M1M2=Double.parseDouble(properties.getProperty(prefix+"probe_M3_M1M2"));
			if (properties.getProperty(prefix+"probe_M2_M1")!=null)
				this.probe_M2_M1=Double.parseDouble(properties.getProperty(prefix+"probe_M2_M1"));
			if (properties.getProperty(prefix+"sigmaToProbe")!=null)
				this.sigmaToProbe=Double.parseDouble(properties.getProperty(prefix+"sigmaToProbe"));
			if (properties.getProperty(prefix+"useTheBest")!=null)
				this.useTheBest=Boolean.parseBoolean(properties.getProperty(prefix+"useTheBest"));
			if (properties.getProperty(prefix+"probeSymmetrical")!=null)
				this.probeSymmetrical=Boolean.parseBoolean(properties.getProperty(prefix+"probeSymmetrical"));
			if (properties.getProperty(prefix+"parallelBeforeProbing")!=null)
				this.parallelBeforeProbing=Boolean.parseBoolean(properties.getProperty(prefix+"parallelBeforeProbing"));
			if (properties.getProperty(prefix+"reProbeDistance")!=null)
				this.reProbeDistance=Double.parseDouble(properties.getProperty(prefix+"reProbeDistance"));
			if (properties.getProperty(prefix+"believeLast")!=null)
				this.believeLast=Double.parseDouble(properties.getProperty(prefix+"believeLast"));
			if (properties.getProperty(prefix+"compensateHysteresis")!=null)
				this.compensateHysteresis=Boolean.parseBoolean(properties.getProperty(prefix+"compensateHysteresis"));
			if (properties.getProperty(prefix+"minCorr")!=null)
				this.minCorr=Double.parseDouble(properties.getProperty(prefix+"minCorr"));
			if (properties.getProperty(prefix+"numFinalCorr")!=null)
				this.numFinalCorr=Integer.parseInt(properties.getProperty(prefix+"numFinalCorr"));
			if (properties.getProperty(prefix+"minCorrPre")!=null)
				this.minCorrPre=Double.parseDouble(properties.getProperty(prefix+"minCorrPre"));
			if (properties.getProperty(prefix+"numFinalCorrPre")!=null)
				this.numFinalCorrPre=Integer.parseInt(properties.getProperty(prefix+"numFinalCorrPre"));
			if (properties.getProperty(prefix+"maxAutoIterations")!=null)
				this.maxAutoIterations=Integer.parseInt(properties.getProperty(prefix+"maxAutoIterations"));
			if (properties.getProperty(prefix+"maxAutoDistance")!=null)
				this.maxAutoDistance=Double.parseDouble(properties.getProperty(prefix+"maxAutoDistance"));
			if (properties.getProperty(prefix+"confirmFirstAuto")!=null)
				this.confirmFirstAuto=Boolean.parseBoolean(properties.getProperty(prefix+"confirmFirstAuto"));
			if (properties.getProperty(prefix+"motorsPreSigma")!=null)
				this.motorsPreSigma=Double.parseDouble(properties.getProperty(prefix+"motorsPreSigma"));
			if (properties.getProperty(prefix+"maxLinearStep")!=null)
				this.maxLinearStep=Double.parseDouble(properties.getProperty(prefix+"maxLinearStep"));
			if (properties.getProperty(prefix+"scanStep")!=null)
				this.scanStep=Integer.parseInt(properties.getProperty(prefix+"scanStep"));
			if (properties.getProperty(prefix+"scanNumber")!=null)
				this.scanNumber=Integer.parseInt(properties.getProperty(prefix+"scanNumber"));
			if (properties.getProperty(prefix+"scanNumberNegative")!=null)
				this.scanNumberNegative=Integer.parseInt(properties.getProperty(prefix+"scanNumberNegative"));
			if (properties.getProperty(prefix+"scanHysteresis")!=null)
				this.scanHysteresis=Boolean.parseBoolean(properties.getProperty(prefix+"scanHysteresis"));
			if (properties.getProperty(prefix+"scanHysteresisNumber")!=null)
				this.scanHysteresisNumber=Integer.parseInt(properties.getProperty(prefix+"scanHysteresisNumber"));

			if (properties.getProperty(prefix+"scanTiltEnable")!=null)
				this.scanTiltEnable=Boolean.parseBoolean(properties.getProperty(prefix+"scanTiltEnable"));
			if (properties.getProperty(prefix+"scanTiltReverse")!=null)
				this.scanTiltReverse=Boolean.parseBoolean(properties.getProperty(prefix+"scanTiltReverse"));
			
			
			if (properties.getProperty(prefix+"scanMeasureLast")!=null)
				this.scanMeasureLast=Boolean.parseBoolean(properties.getProperty(prefix+"scanMeasureLast"));

			if (properties.getProperty(prefix+"scanRunLMA")!=null)
				this.scanRunLMA=Boolean.parseBoolean(properties.getProperty(prefix+"scanRunLMA"));
			
			if (properties.getProperty(prefix+"scanTiltRangeX")!=null)
				this.scanTiltRangeX=Integer.parseInt(properties.getProperty(prefix+"scanTiltRangeX"));
			if (properties.getProperty(prefix+"scanTiltRangeY")!=null)
				this.scanTiltRangeY=Integer.parseInt(properties.getProperty(prefix+"scanTiltRangeY"));
			if (properties.getProperty(prefix+"scanTiltStepsX")!=null)
				this.scanTiltStepsX=Integer.parseInt(properties.getProperty(prefix+"scanTiltStepsX"));
			if (properties.getProperty(prefix+"scanTiltStepsY")!=null)
				this.scanTiltStepsY=Integer.parseInt(properties.getProperty(prefix+"scanTiltStepsY"));
			
			if (properties.getProperty(prefix+"motorHysteresis")!=null)
				this.motorHysteresis=Integer.parseInt(properties.getProperty(prefix+"motorHysteresis"));
			if (properties.getProperty(prefix+"measuredHysteresis")!=null)
				this.measuredHysteresis=Double.parseDouble(properties.getProperty(prefix+"measuredHysteresis"));
			if (properties.getProperty(prefix+"motorCalm")!=null)
				this.motorCalm=Double.parseDouble(properties.getProperty(prefix+"motorCalm"));
			if (properties.getProperty(prefix+"linearReductionRatio")!=null)
				this.linearReductionRatio=Double.parseDouble(properties.getProperty(prefix+"linearReductionRatio"));
			if (properties.getProperty(prefix+"motorDebug")!=null)
				this.motorDebug=Integer.parseInt(properties.getProperty(prefix+"motorDebug"));
			if (properties.getProperty(prefix+"lensDistanceNumPoints")!=null)
				this.lensDistanceNumPoints=Integer.parseInt(properties.getProperty(prefix+"lensDistanceNumPoints"));
			if (properties.getProperty(prefix+"lensDistancePolynomialDegree")!=null)
				this.lensDistancePolynomialDegree=Integer.parseInt(properties.getProperty(prefix+"lensDistancePolynomialDegree"));
			if (properties.getProperty(prefix+"lensDistanceWeightY")!=null)
				this.lensDistanceWeightY=Double.parseDouble(properties.getProperty(prefix+"lensDistanceWeightY"));
			if (properties.getProperty(prefix+"lensDistanceWeightK")!=null)
				this.lensDistanceWeightK=Double.parseDouble(properties.getProperty(prefix+"lensDistanceWeightK"));
			if (properties.getProperty(prefix+"lensDistanceInteractive")!=null)
				this.lensDistanceInteractive=Boolean.parseBoolean(properties.getProperty(prefix+"lensDistanceInteractive"));
			if (properties.getProperty(prefix+"lensDistanceShowResults")!=null)
				this.lensDistanceShowResults=Boolean.parseBoolean(properties.getProperty(prefix+"lensDistanceShowResults"));
			if (properties.getProperty(prefix+"lensDistanceMoveToGoal")!=null)
				this.lensDistanceMoveToGoal=Boolean.parseBoolean(properties.getProperty(prefix+"lensDistanceMoveToGoal"));

			
			if (properties.getProperty(prefix+"powerControlEnable")!=null)
				this.powerControlEnable=Boolean.parseBoolean(properties.getProperty(prefix+"powerControlEnable"));
			if (properties.getProperty(prefix+"powerControlMaximalTemperature")!=null)
				this.powerControlMaximalTemperature=Double.parseDouble(properties.getProperty(prefix+"powerControlMaximalTemperature"));
			if (properties.getProperty(prefix+"powerControlHeaterOnMinutes")!=null)
				this.powerControlHeaterOnMinutes=Double.parseDouble(properties.getProperty(prefix+"powerControlHeaterOnMinutes"));
			if (properties.getProperty(prefix+"powerControlNeitherOnMinutes")!=null)
				this.powerControlNeitherOnMinutes=Double.parseDouble(properties.getProperty(prefix+"powerControlNeitherOnMinutes"));
			if (properties.getProperty(prefix+"powerControlFanOnMinutes")!=null)
				this.powerControlFanOnMinutes=Double.parseDouble(properties.getProperty(prefix+"powerControlFanOnMinutes"));
			
			if (properties.getProperty(prefix+"uvLasersIP")!=null)
				this.uvLasersIP=properties.getProperty(prefix+"uvLasersIP");
			if (properties.getProperty(prefix+"uvLasersBus")!=null)
				this.uvLasersBus=Integer.parseInt(properties.getProperty(prefix+"uvLasersBus"));
			if (properties.getProperty(prefix+"uvLasersCurrents_0")!=null)
				this.uvLasersCurrents[0]=Double.parseDouble(properties.getProperty(prefix+"uvLasersCurrents_0"));
			if (properties.getProperty(prefix+"uvLasersCurrents_1")!=null)
				this.uvLasersCurrents[1]=Double.parseDouble(properties.getProperty(prefix+"uvLasersCurrents_1"));
			if (properties.getProperty(prefix+"uvLasersCurrents_2")!=null)
				this.uvLasersCurrents[2]=Double.parseDouble(properties.getProperty(prefix+"uvLasersCurrents_2"));
			if (properties.getProperty(prefix+"uvLasersCurrents_3")!=null)
				this.uvLasersCurrents[3]=Double.parseDouble(properties.getProperty(prefix+"uvLasersCurrents_3"));
			
			if (properties.getProperty(prefix+"smallestSubPix")!=null)
				this.smallestSubPix=Double.parseDouble(properties.getProperty(prefix+"smallestSubPix"));
			if (properties.getProperty(prefix+"bitmapNonuniforityThreshold")!=null)
				this.bitmapNonuniforityThreshold=Double.parseDouble(properties.getProperty(prefix+"bitmapNonuniforityThreshold"));
			if (properties.getProperty(prefix+"subdiv")!=null)
				this.subdiv=Integer.parseInt(properties.getProperty(prefix+"subdiv"));
			if (properties.getProperty(prefix+"overexposedMaxFraction")!=null)
				this.overexposedMaxFraction=Double.parseDouble(properties.getProperty(prefix+"overexposedMaxFraction"));
			if (properties.getProperty(prefix+"minDefinedArea")!=null)
				this.minDefinedArea=Double.parseDouble(properties.getProperty(prefix+"minDefinedArea"));
			if (properties.getProperty(prefix+"PSFKernelSize")!=null)
				this.PSFKernelSize=Integer.parseInt(properties.getProperty(prefix+"PSFKernelSize"));
			if (properties.getProperty(prefix+"approximateGrid")!=null)
				this.approximateGrid=Boolean.parseBoolean(properties.getProperty(prefix+"approximateGrid"));
			if (properties.getProperty(prefix+"centerPSF")!=null)
				this.centerPSF=Boolean.parseBoolean(properties.getProperty(prefix+"centerPSF"));
			if (properties.getProperty(prefix+"mask1_sigma")!=null)
				this.mask1_sigma=Double.parseDouble(properties.getProperty(prefix+"mask1_sigma"));
			if (properties.getProperty(prefix+"mask1_threshold")!=null)
				this.mask1_threshold=Double.parseDouble(properties.getProperty(prefix+"mask1_threshold"));
			if (properties.getProperty(prefix+"gaps_sigma")!=null)
				this.mask1_threshold=Double.parseDouble(properties.getProperty(prefix+"gaps_sigma"));
			if (properties.getProperty(prefix+"mask_denoise")!=null)
				this.mask_denoise=Double.parseDouble(properties.getProperty(prefix+"mask_denoise"));
			if (properties.getProperty(prefix+"deconvInvert")!=null)
				this.deconvInvert=Double.parseDouble(properties.getProperty(prefix+"deconvInvert"));
			if (properties.getProperty(prefix+"correlationSize")!=null)
				this.correlationSize=Integer.parseInt(properties.getProperty(prefix+"correlationSize"));
			if (properties.getProperty(prefix+"correlationGaussWidth")!=null)
				this.correlationGaussWidth=Double.parseDouble(properties.getProperty(prefix+"correlationGaussWidth"));
			if (properties.getProperty(prefix+"minUVSpan")!=null)
				this.minUVSpan=Double.parseDouble(properties.getProperty(prefix+"minUVSpan"));
			if (properties.getProperty(prefix+"correlationMinContrast")!=null)
				this.correlationMinContrast=Double.parseDouble(properties.getProperty(prefix+"correlationMinContrast"));
			if (properties.getProperty(prefix+"correlationMinInitialContrast")!=null)
				this.correlationMinInitialContrast=Double.parseDouble(properties.getProperty(prefix+"correlationMinInitialContrast"));
			if (properties.getProperty(prefix+"correlationMinAbsoluteContrast")!=null)
				this.correlationMinAbsoluteContrast=Double.parseDouble(properties.getProperty(prefix+"correlationMinAbsoluteContrast"));
			if (properties.getProperty(prefix+"correlationMinAbsoluteInitialContrast")!=null)
				this.correlationMinAbsoluteInitialContrast=Double.parseDouble(properties.getProperty(prefix+"correlationMinAbsoluteInitialContrast"));
			for (int i=0;i<this.postUVscrewSensitivity.length;i++){
				if (properties.getProperty(prefix+"postUVscrewSensitivity_"+i)!=null)
					this.postUVscrewSensitivity[i]=Double.parseDouble(properties.getProperty(prefix+"postUVscrewSensitivity_"+i));
			}
			if (properties.getProperty(prefix+"flatFieldCorrection")!=null)
				this.flatFieldCorrection=Boolean.parseBoolean(properties.getProperty(prefix+"flatFieldCorrection"));
			if (properties.getProperty(prefix+"flatFieldExpand")!=null)
				this.flatFieldExpand=Double.parseDouble(properties.getProperty(prefix+"flatFieldExpand"));
			if (properties.getProperty(prefix+"thresholdFinish")!=null)
				this.thresholdFinish=Double.parseDouble(properties.getProperty(prefix+"thresholdFinish"));
			if (properties.getProperty(prefix+"numIterations")!=null)
				this.numIterations=Integer.parseInt(properties.getProperty(prefix+"numIterations"));
			for (int i=0;i<this.ampsSeconds.length;i++) if (properties.getProperty(prefix+"ampsSeconds_"+i)!=null)
				this.ampsSeconds[i]=Double.parseDouble(properties.getProperty(prefix+"ampsSeconds_"+i));
			if (properties.getProperty(prefix+"reportTemperature")!=null)
				this.reportTemperature=Double.parseDouble(properties.getProperty(prefix+"reportTemperature"));
			if (properties.getProperty(prefix+"showLegacy")!=null)
				this.showLegacy=Boolean.parseBoolean(properties.getProperty(prefix+"showLegacy"));
		}
		public boolean getLensSerial(){
			while (true) { // loop until OK-ed
				GenericDialog gd = new GenericDialog("Enter lens serial number");
				gd.addMessage("Sensor board serial number is "+(((this.serialNumber==null)||(this.serialNumber==""))?"not specified":this.serialNumber));
				gd.addStringField  ("Comment to add to the result files",                 this.comment,80);
				gd.addStringField  ("Lens serial number",this.lensSerial);

				int [] manufacturingIndexMod=getManufacturingIndexMod(this.manufacturingState);
				gd. addChoice("Manufacturing state", this.manufacturingStateNames, this.manufacturingStateNames[manufacturingIndexMod[0]]);
				int maxMod=9;
				if (manufacturingIndexMod[0]<(this.manufacturingStateValues.length-1)){
					maxMod=this.manufacturingStateValues[manufacturingIndexMod[0]+1]-this.manufacturingStateValues[manufacturingIndexMod[0]]-1;
				}
				gd.addNumericField("Optional manufacturing state modifier (0.."+maxMod+")",      manufacturingIndexMod[1], 0,1,"");

				gd.addCheckbox     ("Ask lens serial number on each camera power cycle",this.askLensSerial);
				gd.addNumericField("Required X-shift between the lens axis and the sensor center",      this.centerDeltaX, 0,4,"pix (180 for tilted)");
				gd.addNumericField("Required Y-shift between the lens axis and the sensor center",      this.centerDeltaY, 0,4,"pix");
				gd.showDialog();
				if (gd.wasCanceled()) return false;
				this.comment=                    gd.getNextString();
				this.lensSerial=                 gd.getNextString();
				if (this.lensSerial.length()>0){
					while (this.lensSerial.length()<lensSerialLength) this.lensSerial="0"+this.lensSerial;
				}
				int manIndex=                    gd.getNextChoiceIndex();
				int manMod=                (int) gd.getNextNumber();
				this.askLensSerial=              gd.getNextBoolean();
	    		this.centerDeltaX=               gd.getNextNumber(); 
		    	this.centerDeltaY=               gd.getNextNumber();
				if (manMod<0)           manMod=0;
				else if (manMod>maxMod) manMod=maxMod;
				if (manIndex<manufacturingIndexMod[0]){
					gd = new GenericDialog("Confirm lower manufacturing state");
					gd.addMessage("Selected manufacturing state is \""+ this.manufacturingStateNames[manIndex]+"\" - it is lower than ");
					gd.addMessage("current manufacturing state (\""+ this.manufacturingStateNames[manufacturingIndexMod[0]]+"\".");
					gd.enableYesNoCancel("Yes, update state to lower", "No, re-enter state");
					if (gd.wasCanceled()) return false;
					if (gd.wasOKed()){
						this.manufacturingState=this.manufacturingStateValues[manIndex]+manMod;
						return true;
					}
				} else {
					this.manufacturingState=this.manufacturingStateValues[manIndex]+manMod;
					return true;
				}
			}
		}
// subset of showDialog() - only set parameters realated to scanning		
	   	public boolean showScanningSetup(String title) {
    		GenericDialog gd = new GenericDialog(title);
    		gd.addNumericField("Motor single movement (all 3 motors) in scan focus mode (signed value)",         this.scanStep, 0,7,"motors steps");
    		gd.addNumericField("Number of scan steps during (center) focus scanning",                            this.scanNumber,        0);
    		gd.addNumericField("... of them - in the negative direction (closer lens to sensor)",                this.scanNumberNegative,        0);

    		gd.addCheckbox    ("Scan focus in 2 directions, after the calibration estimate hysteresis (play)",   this.scanHysteresis);
    		gd.addNumericField("Number of scan steps during hysteresis (play) measurement",                      this.scanHysteresisNumber, 0);

    		gd.addCheckbox    ("Scan for tilt measurement (approximately preserving center)",                    this.scanTiltEnable);
    		gd.addCheckbox    ("Scan for tilt measurement in both directions",                                   this.scanTiltReverse);
    		gd.addCheckbox    ("Calculate PSF after returning to the initial position",                          this.scanMeasureLast);
    		gd.addCheckbox    ("Calculate model parameters after scanning",                                      this.scanRunLMA);
    		
    		
    		gd.addNumericField("Full range of scanning motors tilting in X-direction",                           this.scanTiltRangeX, 0,7,"motors steps");
    		gd.addNumericField("Full range of scanning motors tilting in Y-direction",                           this.scanTiltRangeY, 0,7,"motors steps");
    		gd.addNumericField("Number of stops measurements when tilting in X-deirection",                      this.scanTiltStepsX, 0);
    		gd.addNumericField("Number of stops measurements when tilting in Y-deirection",                      this.scanTiltStepsY, 0);
    		gd.addMessage("");
    		gd.addNumericField("Motor anti-hysteresis travel (last measured was "+IJ.d2s(this.measuredHysteresis,0)+")", this.motorHysteresis, 0,7,"motors steps");


    		WindowTools.addScrollBars(gd);
    		gd.showDialog();
    		if (gd.wasCanceled()) return false;

			this.scanStep=             (int) gd.getNextNumber();
			this.scanNumber=           (int) gd.getNextNumber();
            this.scanNumberNegative=   (int) gd.getNextNumber();
			this.scanHysteresis=             gd.getNextBoolean();
			this.scanHysteresisNumber= (int) gd.getNextNumber();

    		this.scanTiltEnable=             gd.getNextBoolean();
    		this.scanTiltReverse=            gd.getNextBoolean();
            this.scanMeasureLast=            gd.getNextBoolean();
            this.scanRunLMA=                 gd.getNextBoolean();
    		
    		this.scanTiltRangeX=       (int) gd.getNextNumber();
    		this.scanTiltRangeY=       (int) gd.getNextNumber();
    		this.scanTiltStepsX=       (int) gd.getNextNumber();
    		this.scanTiltStepsY=       (int) gd.getNextNumber();
    		this.motorHysteresis=      (int) gd.getNextNumber();
    		return true;  
    	}

		
    	public boolean showDialog(String title) { 
    		GenericDialog gd = new GenericDialog(title);
    		//	    		this.serialNumber, // camera serial number string
    		gd.addMessage("Sensor board serial number is "+(((this.serialNumber==null)||(this.serialNumber==""))?"not specified":this.serialNumber));
			gd.addStringField  ("Comment to add to the result files",                 this.comment,80);
    		gd.addStringField  ("Lens serial number",this.lensSerial);
    		
			int [] manufacturingIndexMod=getManufacturingIndexMod(this.manufacturingState);
			gd. addChoice("Manufacturing state", this.manufacturingStateNames, this.manufacturingStateNames[manufacturingIndexMod[0]]);
			int maxMod=9;
			if (manufacturingIndexMod[0]<(this.manufacturingStateValues.length-1)){
				maxMod=this.manufacturingStateValues[manufacturingIndexMod[0]+1]-this.manufacturingStateValues[manufacturingIndexMod[0]]-1;
			}
			gd.addNumericField("Optional manufacturing state modifier (0.."+maxMod+")",      manufacturingIndexMod[1], 0,1,"");

    		
    		gd.addCheckbox     ("Ask lens serial number on each camera power cycle",this.askLensSerial);
    		gd.addCheckbox     ("Add lens serial number to filenames",this.includeLensSerial);
			gd.addNumericField("Required X-shift between the lens axis and the sensor center",      this.centerDeltaX, 0,4,"pix (180 for tilted)");
			gd.addNumericField("Required Y-shift between the lens axis and the sensor center",      this.centerDeltaY, 0,4,"pix");
    		gd.addStringField  ("Grid geometry file",                                 this.gridGeometryFile,40);
			gd.addStringField  ("Initial camera intrinsic/extrinsic parametres file", this.initialCalibrationFile,40);
			gd.addStringField  ("Levenberg-Marquardt algorithm strategy file",        this.strategyFile,40);
			gd.addStringField  ("Focusing results superdirectory (individual will be named by serial numbers)", this.resultsSuperDirectory,40);
			gd.addStringField  ("Measurement history (acquired during \"Scan Calib LMA\") file", this.focusingHistoryFile,80);
			gd.addCheckbox     ("Use lens aberration model (if available) for focal distance and tilts", this.useLMAMetrics);
			gd.addNumericField("EEPROM channel to read sensor serial number from",    this.EEPROM_channel, 0,4,"");
			gd.addCheckbox    ("Save SFE focusing results (including intermediate) ", this.saveResults);
			gd.addCheckbox    ("Show SFE focusing results (including intermediate) ", this.showResults);
			gd.addCheckbox    ("Configure camera",                       !this.cameraIsConfigured);
			gd.addCheckbox    ("Use extra sensor for 3d-location of the optical head", this.useExtraSensor);
			gd.addCheckbox    ("Use circular sensor mask to refine lens center (second LMA run)", this.enRoundOff);
			gd.addCheckbox    ("Keep circular sensor mask after finding lens center",  this.keepCircularMask);
			gd.addNumericField("Camera FOV left margin (clear WOI)",      this.margins.x, 0,4,"pix");
			gd.addNumericField("Camera FOV top margin (clear WOI)",       this.margins.y, 0,4,"pix");
			gd.addNumericField("Camera FOV width (clear WOI)",            this.margins.width, 0,4,"pix");
			gd.addNumericField("Camera FOV height (clear WOI)",           this.margins.height, 0,4,"pix");
    		gd.addNumericField("Number of samples in X direction",        this.numSamples[0], 0,2,"");
    		gd.addNumericField("Number of samples in Y direction",        this.numSamples[1], 0,2,"");
    		gd.addNumericField("Size of sample square (power of 2)",      this.sampleSize, 0,4,"pix");
    		gd.addNumericField("Number of the samples closest to the optical center for \"center focus\"", this.numInCenter, 0);
    		gd.addCheckbox    ("Select samples in the WOI symmetrical around the lens center",this.centerSamples);
    		gd.addNumericField("Maximal grid correction between images",  this.maxCorr, 3,5,"pix");
    		gd.addCheckbox    ("Show history details (per color info)",   this.showHistoryDetails);
    		gd.addCheckbox    ("Show history details for each FOV sample",   this.showHistorySamples);
    		gd.addCheckbox    ("Show history details in a single line (for spreadheets)",  this.showHistorySingleLine);
    		gd.addCheckbox    ("Show acquired images",                    this.showAcquiredImages); // false;
    		gd.addCheckbox    ("Mark ROI on \"Quick get&show\"",          this.showROI); // false;
    		gd.addCheckbox    ("Mark individual samples on \"Quick get&show\"",this.showSamples); // false;
    		gd.addCheckbox    ("Show LMA fitted parameters",              this.showFittedParameters); // true;
    		gd.addCheckbox    ("Use optical head lasers to determine SFE rotation",this.useHeadLasers); // true;
//    		gd.addCheckbox    ("Measure SFE rotation with optical head lasers",  this.useHeadLasers); // true;
    		gd.addMessage("When approximating measured PSF for different areas/colors:");
    		gd.addNumericField("Disregard pixels outside of this fraction of the total energy",           100*this.psf_cutoffEnergy, 2,6,"%");
    		gd.addNumericField("Disregard pixels below this fraction of the maximal value",               100*this.psf_cutoffLevel,  2,6,"%");
    		gd.addNumericField("Minimal selection size (will continue even if previous conditions matched)",  this.psf_minArea,      0,3,"sub-pix");
    		gd.addNumericField("Optionally blur the calculated selection mask",                               this.psf_blurSigma,    2,6,"sub-pix");
    		gd.addNumericField("Weight ratio  red/green (use when combining defocusing data)",                    this.weightRatioRedToGreen,    2,6,"");
    		gd.addNumericField("Weight ratio blue/green (use when combining defocusing data)",                    this.weightRatioBlueToGreen,    2,6,"");
    		gd.addNumericField("Far/Near Focusing target (logariphm of average tangential-to-radial resolution)", this.targetFarNear,    2,6,"");
    		gd.addCheckbox    ("Use targetFarNear (radial/tangential resolution)  as a proxy for the distance",  this.useRadialTangential); // true; // ignore lateral chromatic aberration (center OTF to 0,0)
    		gd.addNumericField("Needed lens center distance (away from the \"best focus\"",                   this.targetMicrons,    2,6,"um");
    		gd.addNumericField("Tolerance of the center focal distance",                                      this.toleranceMicrons,    2,6,"um");
    		gd.addNumericField("Tolerance of the tilt ",                                                      this.toleranceTilt,    2,6,"");
    		gd.addNumericField("When under this scaled tolerance, correction step is scaled down twice",      this.toleranceThreshold,    2,6,"");
    		gd.addNumericField("Adjust 3 motors by parallel movement if focal distance error in the center exceeds this value", this.parallelAdjustThreshold,    2,6,"um");
    		gd.addNumericField("When fitting far/near and tilts, samples 'decay' with this sigma",                this.motorsSigma,      1,7,"motors steps");
    		gd.addNumericField("Same for the 3 motors moving together (center focusing)",                         this.motorsSigma3,      1,7,"motors steps");
    		gd.addNumericField("Low limit for motor sigma when fitting is in the final stage",                this.motorsMinSigma,      1,7,"motors steps");
    		gd.addNumericField("When walk is getting smaller, sigma be reduced proportionally",                this.motorsVarSigmaToTravel,      1,7,"motors steps");
    		gd.addNumericField("After each step new sigma will have this part of the calculated from the travel",                this.motorsFadeSigma,      1,7,"x");
    		gd.addNumericField("For quadratic maximum the correction will be increased by (1+motorsOverShootToBalance) if there are less samples on the other side",                this.motorsOverShootToBalance,      1,7,"motors steps");
    		gd.addCheckbox    ("When measuring tilt, use those with good center with higher weight",       this.filterGoodDistance);
    		gd.addNumericField("Sigma for the weight function of tilt measurements, depending on the center distance error", this.goodDistanceSigma,      1,7,"microns");

    		gd.addNumericField("Weight decay for heavily tilted samples ", this.goodTiltSigma,      3,7,"");

    		gd.addNumericField("Maximal allowed single-step focusing adjustment",                                   this.maxStep,      1,7,"motors steps");
    		gd.addNumericField("How far to go to probe around the current point to measure derivatives",            this.probeStep,      1,7,"motors steps");
    		
			
    		gd.addNumericField("How far too move average of the 3 motors: (M1+M2+M3)/3 during probing",             this.probe_M1M2M3,      1,7,"motors steps");
    		gd.addNumericField("How far to move M3 opposite to M1 and M2: M3-(M1+M2)/2 during probing",             this.probe_M3_M1M2,      1,7,"motors steps");
    		gd.addNumericField("How far to move M2 opposite to M1: M2-M1 during probing",                           this.probe_M2_M1,      1,7,"motors steps");
    		gd.addNumericField("Data from far samples decays proportionally to the probe distances with this scale",this.sigmaToProbe,      1,7,"x");
    		gd.addCheckbox    ("Adjust from the best known position (false - from the last)",                       this.useTheBest);
    		gd.addCheckbox    ("Probe 6 measurements for each motor, if unchecked - only 4 (tetrahedron)",       this.probeSymmetrical);
    		gd.addCheckbox    ("Readjust focus before probing around",       this.parallelBeforeProbing);
    		gd.addNumericField("Re-run probing around in orthogonal directions, if the current position moved farther from the last probing one", this.reProbeDistance, 1,7,"motors steps");
    		gd.addNumericField("Bias towards the last measurement: 0% use best fit for all, 100% make planes throgh the last point",100*this.believeLast,      1,5,"%");
    		gd.addCheckbox    ("Move motors in the same direction to compensate for the hysteresis",             this.compensateHysteresis);
    		gd.addNumericField("Motor anti-hysteresis travel (last measured was "+IJ.d2s(this.measuredHysteresis,0)+")", this.motorHysteresis, 0,7,"motors steps");
    		gd.addNumericField("wait after motors reached final position before acquiring image ",               this.motorCalm, 0,7,"seconds");
    		gd.addNumericField("Sensor travel to motors travel (all 3 together), by design it is 4/38~=0.105",   this.linearReductionRatio, 5,7,"");
    		gd.addNumericField("Motors debug (1 - show moves, 2 show moves+hysteresis mioves)",                  this.motorDebug,        0);
    		gd.addMessage("Parameters to calculate lens distance measurements from 3-color PSF measurements");
    		gd.addNumericField("Number of points to tabulate center focus parameters vs. focal distance",        this.lensDistanceNumPoints,        0);
    		gd.addNumericField("polynomial degree to approximate center focus parameters vs. focal distance ",   this.lensDistancePolynomialDegree,        0);
    		gd.addNumericField("Normalize overall sharpness (that depends on the lens quality and/or PSF parameters) to differential ones",            this.lensDistanceWeightY,        3,5,"");
    		gd.addNumericField("Use derivartives-dependent weights for components (1.0), equal weight - 0.0",    this.lensDistanceWeightK,        3,5,"");
    		gd.addCheckbox    ("Open dialog when calibrating focal distance",                                    this.lensDistanceInteractive);
    		gd.addCheckbox    ("Show results window from focal distance calibration",                            this.lensDistanceShowResults);
    		gd.addCheckbox    ("Move motors together to the requested microns from the \"best focus\"",          this.lensDistanceMoveToGoal);
    		
    		gd.addCheckbox    ("Enable power control for heater and fan",                                        this.powerControlEnable);
    		gd.addNumericField("Maximal allowed temperature",                                                    this.powerControlMaximalTemperature,  3,5,"C");
    		gd.addNumericField("Heater ON time",                                                                 this.powerControlHeaterOnMinutes,  1,5,"min");
    		gd.addNumericField("Both heater and fan OFF time",                                                   this.powerControlNeitherOnMinutes,  1,5,"min");
    		gd.addNumericField("Fan ON time",                                                                    this.powerControlFanOnMinutes,  1,5,"min");
    		
			gd.addStringField  ("IP address of the camera with 103641 board (UV LEDs and lasers) are attached",  this.uvLasersIP,40);
    		gd.addNumericField("I2C bus where LED/laser board is attached (0 - through 10359, 1 - through 10369)",this.uvLasersBus,        0);
    		gd.addNumericField("UV LED1 \"on\" current (left/near  when looking from the target)",               this.uvLasersCurrents[0],  3,5,"mA");
    		gd.addNumericField("UV LED2 \"on\" current (right/near when looking from the target)",               this.uvLasersCurrents[0],  3,5,"mA");
    		gd.addNumericField("UV LED3 \"on\" current (right/far  when looking from the target)",               this.uvLasersCurrents[0],  3,5,"mA");
    		gd.addNumericField("UV LED4 \"on\" current (left/far   when looking from the target)",               this.uvLasersCurrents[0],  3,5,"mA");
			
    		gd.addMessage("");
    		gd.addNumericField("Minimal correction movement to initiate final series of corrections - focus/tilt mode",            this.minCorr,        1,5,"motors steps");
    		gd.addNumericField("Finish if this number of last corrections where below  minimum (previous input) - focus/tilt mode",this.numFinalCorr,     0);
    		gd.addNumericField("Minimal correction movement to initiate final series of corrections - pre-focus mode",            this.minCorrPre,        1,5,"motors steps");
    		gd.addNumericField("Finish if this number of last corrections where below  minimum (previous input) - pre-focus mode",this.numFinalCorrPre,     0);
    		gd.addNumericField("Unconditionally exit focusing adjustment after these number of itereations",     this.maxAutoIterations,        0);
    		gd.addNumericField("Maximal allowed total motors travel for automatic correction",                   this.maxAutoDistance, 1,7,"motors steps");
    		gd.addCheckbox    ("Ask confirmation after first automatic adjustment step (before moving)",         this.confirmFirstAuto);
    		gd.addNumericField("When fitting parabola for focusing sharpness in the center, far measurements decay with this sigma", this.motorsPreSigma, 1,7,"motors steps");
    		gd.addNumericField("If there are insufficient measurements to fit parabola - make this step",        this.maxLinearStep, 1,7,"motors steps");
    		gd.addNumericField("Motor single movement (all 3 motors) in scan focus mode (signed value)",         this.scanStep, 0,7,"motors steps");
    		gd.addNumericField("Number of scan steps during (center) focus scanning",                            this.scanNumber,        0);
    		gd.addNumericField("... of them - in the negative direction (closer lens to sensor)",                this.scanNumberNegative,        0);
   		
    		gd.addCheckbox    ("Scan focus in 2 directions, after the calibration estimate hysteresis (play)",   this.scanHysteresis);
    		gd.addNumericField("Number of scan steps during hysteresis (play) measurement",                      this.scanHysteresisNumber, 0);

    		gd.addCheckbox    ("Scan for tilt measurement (approximately preserving center)",                    this.scanTiltEnable);
    		gd.addCheckbox    ("Scan for tilt measurement in both directions",                                   this.scanTiltReverse);
    		gd.addCheckbox    ("Calculate PSF after returning to the initial position",                          this.scanMeasureLast);
    		
    		
    		gd.addNumericField("Full range of scanning motors tilting in X-direction",                           this.scanTiltRangeX, 0,7,"motors steps");
    		gd.addNumericField("Full range of scanning motors tilting in Y-direction",                           this.scanTiltRangeY, 0,7,"motors steps");
    		gd.addNumericField("Number of stops measurements when tilting in X-deirection",                      this.scanTiltStepsX, 0);
    		gd.addNumericField("Number of stops measurements when tilting in Y-deirection",                      this.scanTiltStepsY, 0);
    		
    		gd.addMessage     ("The following parameters overwrite some defined for aberration measurements in other dialogs");
    		gd.addNumericField("Smallest fraction to subdivide pixels at simulation", this.smallestSubPix, 3,5,"sensor pix");
    		gd.addNumericField("Maximal difference of the pattern value in the corners that triggers subdivision", this.bitmapNonuniforityThreshold, 3);
    		gd.addNumericField("Subdivide simulated pattern by:",         this.subdiv, 0);
    		gd.addNumericField("Allowed overexposed pixels (fraction of the area) ",this.overexposedMaxFraction,3); //  0.005; // allowed fraction of the overexposed pixels in the PSF kernel measurement area 
    		gd.addNumericField("Min fraction of the FFT square (weighted) to have defined pattern",  this.minDefinedArea, 3);
    		gd.addNumericField ("PSF kernel size",                        this.PSFKernelSize, 0);
    		gd.addCheckbox    ("Approximate pattern grid with a polynomial",this.approximateGrid); // true; // ignore lateral chromatic aberration (center OTF to 0,0)
    		gd.addCheckbox    ("Center PSF by modifying phase",           this.centerPSF); // true; // ignore lateral chromatic aberration (center OTF to 0,0)
    		gd.addNumericField("Bluring power spectrum to remove pattern grid (in pattern base freq)",  this.mask1_sigma, 3);
    		gd.addNumericField("Threshold to supress spectral points not present in the pattern ",  this.mask1_threshold, 3);
    		gd.addNumericField("Sigma for filling the OTF ",               this.gaps_sigma, 3);
    		gd.addNumericField("Denoise mask ",                            this.mask_denoise, 3);
    		gd.addNumericField("Invert deconvolution if less than",        this.deconvInvert, 3);
    		gd.addMessage(     "--- Pattern correlation parameters ---");
			gd.addNumericField("Correlation size:",                        this.correlationSize, 0); // 64
			gd.addNumericField("Correlation Gauss width (relative):",      this.correlationGaussWidth, 3);
			gd.addNumericField("Minimal UV span in correlation window to trigger FFT size increase",this.minUVSpan, 3);
			gd.addNumericField("Correlation minimal contrast (normalized)",                    this.correlationMinContrast, 3);
			gd.addNumericField("Correlation minimal contrast for initial search (normalized)", this.correlationMinInitialContrast, 3);
			gd.addNumericField("Correlation minimal contrast (absolute)",                      this.correlationMinAbsoluteContrast, 3);
			gd.addNumericField("Correlation minimal contrast for initial search (absolute)",   this.correlationMinAbsoluteInitialContrast, 3);
			gd.addCheckbox    ("Compensate uneven pattern intensity",      this.flatFieldCorrection);
			gd.addNumericField("Expand during extrapolation (relative to the average grid period)", this.flatFieldExpand, 3);
			gd.addNumericField("Threshold RMS to exit LMA",                this.thresholdFinish, 7,9,"pix");
			gd.addNumericField("Maximal number of LMA iterations per series",this.numIterations, 0);
    		gd.addMessage(     "--- Post-UV fixture screws sensitivity ---");
			for (int i=0;i<this.postUVscrewSensitivity.length;i++){
				gd.addNumericField("Screw "+i+" sensitivity",              this.postUVscrewSensitivity[i], 4,6,"um/turn CW");
			}
    		
    		gd.addMessage("-----");
    		gd.addNumericField("Report focal length at this temperature", this.reportTemperature, 1,5,"C");
			gd.addCheckbox    ("Show legacy focusing parameters (most are already not supported anyway)",      this.showLegacy);
			
    		if (!Double.isNaN(this.sensorTemperature)) gd.addMessage("Last measured sensor temperature is "+this.sensorTemperature+" C");

    		if (!Double.isNaN(this.result_lastKT)) gd.addMessage("Temperature focal distance coefficient measured in last run is "+this.result_lastKT+"microns/C");
    		if (!Double.isNaN(this.result_lastFD20)) gd.addMessage("Focal distance @20C measured at last run is "+this.result_lastFD20+" microns");
    		if (!Double.isNaN(this.result_lastKT) && !Double.isNaN(this.result_lastFD20)){
    			gd.addMessage("Focal distance @"+this.reportTemperature+"C measured at last run is "+
    					(this.result_lastFD20+(this.reportTemperature-20.0)*this.result_lastKT)+" microns");	
    		}
    		if (!Double.isNaN(this.result_allHistoryKT)) gd.addMessage("Temperature focal distance coefficient calculated from all measurements is "+this.result_allHistoryKT+" microns");
    		if (!Double.isNaN(this.result_allHistoryFD20)) gd.addMessage("Focal distance @20C calculated from all measurements is "+this.result_allHistoryFD20+" microns");
    		if (!Double.isNaN(this.result_allHistoryKT) && !Double.isNaN(this.result_allHistoryFD20)){
    			gd.addMessage("Focal distance @"+this.reportTemperature+"C calculated from all measurements is "+
    					(this.result_allHistoryFD20+(this.reportTemperature-20.0)*this.result_allHistoryKT)+" microns");	
    		}
    		if (!Double.isNaN(this.result_fDistance)) gd.addMessage("Focal distance is "+this.result_fDistance+" microns");
    		if (!Double.isNaN(this.result_tiltX)) gd.addMessage("Horizontal angular/tangential asymmetry "+this.result_tiltX);
    		if (!Double.isNaN(this.result_tiltY)) gd.addMessage("Vertical angular/tangential asymmetry "+this.result_tiltY);
    		if (!Double.isNaN(this.result_R50)) gd.addMessage("Average PSF radius at 50% level is (higher than actual) "+this.result_R50+" pixels");
    		if (!Double.isNaN(this.result_A50)) gd.addMessage("Same, but with averaged r^2: "+this.result_A50+" C");
    		if (!Double.isNaN(this.result_B50)) gd.addMessage("Same, but with averaged r^4: "+this.result_B50+" C");
    		if (!Double.isNaN(this.result_RC50)) gd.addMessage("Average PSF 50% radius for center samples "+this.result_RC50+" pixels");
    		
    		if (!Double.isNaN(this.result_PX0)) gd.addMessage("Lens center X-coordinate on the sensor "+this.result_PX0+" pixels");
    		if (!Double.isNaN(this.result_PY0)) gd.addMessage("Lens center Y-coordinate on the sensor  "+this.result_PY0+" pixels");
    		if (!Double.isNaN(this.result_PSI)) gd.addMessage("SFE rotation relative to target (clockwise - positive) "+this.result_PSI+" degrees");
    		if (!Double.isNaN(this.result_ROT)) gd.addMessage("SFE rotation relative to optical head lasers (clockwise - positive) "+this.result_ROT+" degrees");
    		if (!Double.isNaN(this.result_FocalLength)) gd.addMessage("Lens focal length "+this.result_FocalLength+" mm");
			gd.addMessage("Cumulative currents that ran through UV LEDs:");
			for (int i=0;i<this.ampsSeconds.length;i++) gd.addMessage("UV LED "+(i+1)+":"+IJ.d2s(this.ampsSeconds[i],3)+" coulombs  (amp-seconds)");
			
    		WindowTools.addScrollBars(gd);
    		gd.showDialog();
    		if (gd.wasCanceled()) return false;
    		this.comment=                    gd.getNextString();
    		this.lensSerial=                 gd.getNextString();
    		if (this.lensSerial.length()>0){
    			while (this.lensSerial.length()<lensSerialLength) this.lensSerial="0"+this.lensSerial;
    		}
    		
			int manIndex=                    gd.getNextChoiceIndex();
			int manMod=                (int) gd.getNextNumber();
			if (manMod<0)           manMod=0;
			else if (manMod>maxMod) manMod=maxMod;
			this.manufacturingState=this.manufacturingStateValues[manIndex]+manMod; // here no restriction on the order - can go backwards
			
    		this.askLensSerial=              gd.getNextBoolean();
    		this.includeLensSerial=          gd.getNextBoolean();

    		this.centerDeltaX=               gd.getNextNumber(); 
	    	this.centerDeltaY=               gd.getNextNumber();

			this.gridGeometryFile=           gd.getNextString();
			this.initialCalibrationFile=     gd.getNextString();
			this.strategyFile=               gd.getNextString();
    		this.resultsSuperDirectory=      gd.getNextString();
    		this.focusingHistoryFile=        gd.getNextString();
			this.useLMAMetrics =             gd.getNextBoolean();

    		this.EEPROM_channel=       (int) gd.getNextNumber();
			this.saveResults=                gd.getNextBoolean();
			this.showResults=                gd.getNextBoolean();
//    		this.comment=                    gd.getNextString().replace(' ','_'); //TODO: - add escape
    		this.configureCamera=            gd.getNextBoolean();
    		this.useExtraSensor=             gd.getNextBoolean();
    		this.enRoundOff=                 gd.getNextBoolean();
    		this.keepCircularMask=           gd.getNextBoolean();
    		this.margins.x=            (int) gd.getNextNumber();
    		this.margins.y=            (int) gd.getNextNumber();
    		this.margins.width=        (int) gd.getNextNumber();
    		this.margins.height=       (int) gd.getNextNumber();
    		this.numSamples[0]=        (int) gd.getNextNumber();
    		this.numSamples[1]=        (int) gd.getNextNumber();
    		int inputNumber=           (int) gd.getNextNumber();
    		for (this.sampleSize=16; this.sampleSize<inputNumber; this.sampleSize<<=1);
    		this.numInCenter=           (int) gd.getNextNumber();
    		this.centerSamples=              gd.getNextBoolean();
    		this.maxCorr=                    gd.getNextNumber();
    		this.showHistoryDetails=         gd.getNextBoolean();
    		this.showHistorySamples=         gd.getNextBoolean();
    		this.showHistorySingleLine=      gd.getNextBoolean();
    		this.showAcquiredImages=         gd.getNextBoolean();
			this.showROI=                    gd.getNextBoolean();
			this.showSamples=                gd.getNextBoolean();
    		this.showFittedParameters=       gd.getNextBoolean();
			this.useHeadLasers=              gd.getNextBoolean();
    		this.psf_cutoffEnergy=      0.01*gd.getNextNumber();
    		this.psf_cutoffLevel=       0.01*gd.getNextNumber();
    		this.psf_minArea=          (int) gd.getNextNumber();
    		this.psf_blurSigma=              gd.getNextNumber();
			this.weightRatioRedToGreen=      gd.getNextNumber();
			this.weightRatioBlueToGreen=     gd.getNextNumber();
			this.targetFarNear=              gd.getNextNumber();
			this.useRadialTangential=        gd.getNextBoolean();  
    		this.targetMicrons=              gd.getNextNumber();
    		this.toleranceMicrons=           gd.getNextNumber();
    		this.toleranceTilt=              gd.getNextNumber();
    		this.toleranceThreshold=         gd.getNextNumber();
    		this.parallelAdjustThreshold=    gd.getNextNumber();
			this.motorsSigma=                gd.getNextNumber();
			this.motorsSigma3=               gd.getNextNumber();
			this.motorsMinSigma=             gd.getNextNumber(); 
			this.motorsVarSigmaToTravel=     gd.getNextNumber();
			this.motorsFadeSigma=            gd.getNextNumber(); 
			this.motorsOverShootToBalance=   gd.getNextNumber();
			this.filterGoodDistance=         gd.getNextBoolean();
			this.goodDistanceSigma=          gd.getNextNumber();
			this.goodTiltSigma=              gd.getNextNumber();
			this.maxStep=                    gd.getNextNumber();
			this.probeStep=                  gd.getNextNumber();
			this.probe_M1M2M3=               gd.getNextNumber();
			this.probe_M3_M1M2=              gd.getNextNumber();
			this.probe_M2_M1=                gd.getNextNumber();
			this.sigmaToProbe=               gd.getNextNumber(); 
			this.useTheBest=                 gd.getNextBoolean();  	                
			this.probeSymmetrical=           gd.getNextBoolean();
			this.parallelBeforeProbing=      gd.getNextBoolean();
			this.reProbeDistance=            gd.getNextNumber();
			this.believeLast=           0.01*gd.getNextNumber(); if (this.believeLast>1.0) this.believeLast=1.0; else if (this.believeLast<0) this.believeLast=0.0;
    		this.compensateHysteresis=       gd.getNextBoolean();
    		this.motorHysteresis=      (int) gd.getNextNumber();
    		this.motorCalm=                   gd.getNextNumber();
			this.linearReductionRatio=       gd.getNextNumber();
    		this.motorDebug=           (int) gd.getNextNumber();
    		this.lensDistanceNumPoints=(int) gd.getNextNumber();
    		this.lensDistancePolynomialDegree=(int) gd.getNextNumber();
    		this.lensDistanceWeightY=        gd.getNextNumber();
    		this.lensDistanceWeightK=        gd.getNextNumber();
			this.lensDistanceInteractive=    gd.getNextBoolean();
			this.lensDistanceShowResults=    gd.getNextBoolean();
			this.lensDistanceMoveToGoal=     gd.getNextBoolean();

			this.powerControlEnable=            gd.getNextBoolean();
			this.powerControlMaximalTemperature=gd.getNextNumber();
			this.powerControlHeaterOnMinutes=   gd.getNextNumber();
			this.powerControlNeitherOnMinutes=  gd.getNextNumber();
			this.powerControlFanOnMinutes=      gd.getNextNumber();
			
			this.uvLasersIP=                 gd.getNextString();
			this.uvLasersBus=          (int) gd.getNextNumber();
			this.uvLasersCurrents[0]=        gd.getNextNumber();
			this.uvLasersCurrents[1]=        gd.getNextNumber();
			this.uvLasersCurrents[2]=        gd.getNextNumber();
			this.uvLasersCurrents[3]=        gd.getNextNumber();
			
			this.minCorr=                    gd.getNextNumber();
			this.numFinalCorr=         (int) gd.getNextNumber();
			this.minCorrPre=                 gd.getNextNumber();
			this.numFinalCorrPre=      (int) gd.getNextNumber();
			this.maxAutoIterations=    (int) gd.getNextNumber();
			this.maxAutoDistance=            gd.getNextNumber();
			this.confirmFirstAuto=           gd.getNextBoolean();
			this.motorsPreSigma=             gd.getNextNumber();
			this.maxLinearStep=              gd.getNextNumber();
			this.scanStep=             (int) gd.getNextNumber();
			this.scanNumber=           (int) gd.getNextNumber();
            this.scanNumberNegative=   (int) gd.getNextNumber();
			this.scanHysteresis=             gd.getNextBoolean();
			this.scanHysteresisNumber= (int) gd.getNextNumber();
			
    		this.scanTiltEnable=             gd.getNextBoolean();
    		this.scanTiltReverse=            gd.getNextBoolean();
    		this.scanMeasureLast=            gd.getNextBoolean();
    		this.scanTiltRangeX=       (int) gd.getNextNumber();
    		this.scanTiltRangeY=       (int) gd.getNextNumber();
    		this.scanTiltStepsX=       (int) gd.getNextNumber();
    		this.scanTiltStepsY=       (int) gd.getNextNumber();
			
			
    		this.smallestSubPix=             gd.getNextNumber();
    		this.bitmapNonuniforityThreshold=gd.getNextNumber();
    		this.subdiv=               (int) gd.getNextNumber();
    		this.overexposedMaxFraction=     gd.getNextNumber(); 
    		this.minDefinedArea=             gd.getNextNumber();
    		this.PSFKernelSize=        (int) gd.getNextNumber();
    		this.approximateGrid=            gd.getNextBoolean();
    		this.centerPSF=                  gd.getNextBoolean();
    		this.mask1_sigma=                gd.getNextNumber();
    		this.mask1_threshold=            gd.getNextNumber();
    		this.gaps_sigma=                 gd.getNextNumber();
    		this.mask_denoise=               gd.getNextNumber();
    		this.deconvInvert=               gd.getNextNumber();
			this.correlationSize=      (int) gd.getNextNumber();
			this.correlationGaussWidth=      gd.getNextNumber();
			this.minUVSpan=                  gd.getNextNumber();
			this.correlationMinContrast=  gd.getNextNumber();
			this.correlationMinInitialContrast=  gd.getNextNumber();
			this.correlationMinAbsoluteContrast=  gd.getNextNumber();
			this.correlationMinAbsoluteInitialContrast=  gd.getNextNumber();
			this.flatFieldCorrection=        gd.getNextBoolean();
			this.flatFieldExpand=            gd.getNextNumber();
			this.thresholdFinish=            gd.getNextNumber();
			this.numIterations=        (int) gd.getNextNumber();
			for (int i=0;i<this.postUVscrewSensitivity.length;i++){
				this.postUVscrewSensitivity[i]=            gd.getNextNumber();
			}
    		this.reportTemperature=          gd.getNextNumber();
			this.showLegacy=                 gd.getNextBoolean();
    		return true;
    	}
/* ======================================================================== */
    	//returns triads - x,y,distance from the lens center 
    	public double [][][] sampleCoordinates(
    			double x0,   // lens center on the sensor
    			double y0){  // lens center on the sensor
    		int ix0=(int) Math.round(x0);
    		int iy0=(int) Math.round(y0);
    		
//    		Rectangle woi=new Rectangle(this.margins);
    		Rectangle woi=this.getMargins();
//    		System.out.println("Selection Rectangle("+woi.x+", "+woi.y+", "+woi.width+", "+woi.height+ ")");
    		int extra=this.sampleSize+2*(this.numSamples[1]+this.numSamples[0]);
    		if ((this.centerSamples) &&
    				(woi.x<(x0-extra)) &&
    				(woi.y<(y0-extra)) &&
    				((woi.x+woi.width)>(x0+extra)) &&
    				((woi.y+woi.height)>(y0+extra))){
    			int halfWidth= ((woi.x+woi.width/2)>ix0)? (ix0-woi.x) : (woi.x+woi.width-ix0);
    			int halfHeight=((woi.y+woi.height/2)>iy0)?(iy0-woi.y) : (woi.y+woi.height-iy0);
    			woi.x=(ix0-halfWidth);
    			woi.y=(iy0-halfHeight);
    			woi.width=2*halfWidth;
    			woi.height=2*halfHeight;
    //			System.out.println("ix0="+ix0+" iy0="+iy0+", halfWidth="+halfWidth+" halfHeight="+halfHeight+"  Modified selection= Rectangle("+woi.x+", "+woi.y+", "+woi.width+", "+woi.height+")");
    //			System.out.println("Original selection was= Rectangle("+this.margins.x+", "+this.margins.y+", "+this.margins.width+", "+this.margins.height+")");
    		} else {
//    			System.out.println("ix0="+ix0+" iy0="+iy0+", Modified selection= Rectangle("+woi.x+", "+woi.y+", "+woi.width+", "+woi.height+")");
//    			System.out.println("sampleCoordinates() - NO CENTERING!, extra="+extra);
    		}
    		if ((woi.width<extra) || (woi.height<extra)){
    			String msg="Selection Rectangle("+woi.x+", "+woi.y+", "+woi.width+", "+woi.height+ ")+is too small for the meashurements";
    			IJ.showMessage("Error",msg);
    			throw new IllegalArgumentException (msg);
    		}
    		woi.x+=this.sampleSize/2;
    		woi.y+=this.sampleSize/2;
    		woi.width-=this.sampleSize;
    		woi.height-=this.sampleSize;
//    		System.out.println("Selection - sample centers: Rectangle("+woi.x+", "+woi.y+", "+woi.width+", "+woi.height+ ")");
    		
    		double [][][] sampleCoord=new double[this.numSamples[1]][this.numSamples[0]][3];
    		int y,x;
    		for (int i=0;i<this.numSamples[1];i++){
    			if (numSamples[1]<=1) {
    			  if ((iy0>=woi.y) && (iy0<(woi.y+woi.height))) y=iy0;
    			  else y=woi.y+woi.height/2;
    			} else y= (woi.y+(woi.height*i)/(numSamples[1]-1));
    			y&=~1; // make it even
    			for (int j=0;j<this.numSamples[0];j++){
    				if (numSamples[0]<=1) {
    					if ((ix0>=woi.x) && (ix0<(woi.x+woi.width))) x=ix0;
    					else x=woi.x+woi.width/2;
    				} else x= (woi.x+(woi.width*j)/(numSamples[0]-1));
    				x&=~1; // make it even
    				sampleCoord[i][j][0]=x;
    				sampleCoord[i][j][1]=y;
    				sampleCoord[i][j][2]=Math.sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)); // distance from tyhe center
    			}
    		}
/*	    		
    		for (int i=0;i<this.numSamples[1];i++){
    			for (int j=0;j<this.numSamples[0];j++){
    				System.out.println ("i="+i+" j="+j+" sampleCoord[i][j][0]="+sampleCoord[i][j][0]+" sampleCoord[i][j][1]="+sampleCoord[i][j][1]);
    			}
    		}
*/
    		return sampleCoord;
    	}
		public boolean [][] getCenterMask(
				double x0,   // lens center on the sensor
    			double y0,  // lens center on the sensor
                int size){
			double [][][] sampleCoord=sampleCoordinates(x0,y0);
			boolean [][] mask =new boolean[sampleCoord.length][sampleCoord[0].length];
			for (int i=0;i<mask.length;i++) for (int j=0;j<mask.length;j++) mask[i][j]=false;
			if (size>(mask.length*mask[0].length)) size=(mask.length*mask[0].length);
			for (int n=0;n<size;n++){
				double minDist=Double.NaN;
				int i0=0,j0=0;
    			for (int i=0;i<mask.length;i++) for (int j=0;j<mask.length;j++) if (!mask[i][j]) {
    				if (!(minDist<sampleCoord[i][j][2])) {
    					minDist=sampleCoord[i][j][2];
    					i0=i;
    					j0=j;
    				}
    			}
				mask[i0][j0]=true;
			}
			return mask;
		}
		public String showSamplesMap(
				double x0,   // lens center on the sensor
    			double y0,  // lens center on the sensor
                int size){
			double [][][] sampleCoord=sampleCoordinates(x0,y0);
			boolean [][] mask=getCenterMask(x0,y0,size);
			StringBuffer sb=new StringBuffer();
//   			String s="";
			for (int i=0;i<mask.length;i++) {
				for (int j=0;j<mask[0].length;j++) {
					if (mask[i][j])	sb.append("["+IJ.d2s(sampleCoord[i][j][0],1)+":"+IJ.d2s(sampleCoord[i][j][1],1)+" ("+IJ.d2s(sampleCoord[i][j][2],1)+")] ");
					else         	sb.append(" "+IJ.d2s(sampleCoord[i][j][0],1)+":"+IJ.d2s(sampleCoord[i][j][1],1)+" ("+IJ.d2s(sampleCoord[i][j][2],1)+")  ");
//   					if (mask[i][j])	s+=("["+IJ.d2s(sampleCoord[i][j][2],1)+"] ");
//  					else         	s+=(" "+IJ.d2s(sampleCoord[i][j][2],1)+"  ");
				}
//   				s+="\n";
				sb.append("\n");
			}
//			System.out.println(s);
//   			System.out.println(sb.toString());
			return sb.toString();
			
		}

    }

	
	
	
}
