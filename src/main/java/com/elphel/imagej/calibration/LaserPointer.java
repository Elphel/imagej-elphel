package com.elphel.imagej.calibration;

import java.util.Properties;

import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.PolynomialApproximation;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.common.WindowTools;

import ij.IJ;
import ij.gui.GenericDialog;

public class LaserPointer{
	public double headLasersTilt=  1.06; // degrees, right laser lower than left laser
	public double minimalIntensity=0.05; // of scaled saturation when laser is on
	public double maximalIntensity=1.5;  // of scaled saturation when laser is off
	public int    overexposedRadius = 30;   // no pointers closer than this to overexposed areas
	public double lowpassSigma=1.0; // 0.8;    // low pass sigma, in pixels
	public double highpassSigma=20;    // high pass sigma, in pixels
	public double headLowpassSigma=0.8;    // low pass sigma, in pixels for optical head lasers
	public double quadraticScaleSigma=1.0; // find local maximum by quadratic intrepolating pixels around maximal value (relative to lkow pass sigma)
	public int    algorithmNumber=4;
	public int    closestOffender=3;
	public int    fartherstOffender=200;
	public double fatZero =0.05;
	public double greenFloor=0.6;      // when dividing by green, add this fraction of maximal value (decrease green accordingly)
	public boolean useOther=true; // when true - use red and other color, when false - only red
	public boolean otherGreen=true; // other color is green (false - blue)
	public double threshold=0.1;
	// default grid orientation, used if not enough pointers visible (modified when more visible)
	public boolean swapUV=false; // first
	public boolean flipU=false;
	public boolean flipV=false;
	public boolean whiteOnly=true; // verify laser is on the white pattern cell
	public double  maxOffsetFromCenter=0.6; // maximal offset of the laser spot from the center, relative to cell radius
	public double [][] laserUVMap; // first index - number of pointer points
	// new variables TODO: add handling (all linear dimensions in sensor pixels)
	public double laserSignalToNoise=1.5; // Minimal signal-to-noise ratio for laser pointers
	public double localMaxRadius=10; // sensor pix. currently uses just square (2*localMaxRadius+1)**2
	public boolean usePatternFilter=true; // Filter laser positions by likely pattern white cells
	public int decimatePatternFilter=2; // reduce resolution for pattern filter
	public double localContrastSigma=40; // use to calculate local level and contrast
	public double localToGlobalContrast=0.8; // 0 - same contrast normalization for the whole image, 1.0 - pure local
	public double patternLowPassSigma=4.0; // filter normalized patetrn before thresholding
	public double patternThreshold=0.2; // fraction of dispersion (same positive for white cells, negative for black ones)
	public double maximalCellSize=30.0; // White cells should have black pixels in all 4 quadrants not farther than this
	public int numPasses=3; // number of black/white alternations of the surrounding cells to use in quadrant filtering
	public boolean bordersOK=false; // frame border as good cell for quadrant filter
	public double blurredMaskThreshold=0.1; // select only areas with multiple pattern white cells
	public double maskGrow=2.0; // grow final mask (pixels)
	public int    debugLevel=1;
	private long startTime=System.nanoTime();
	private long lastTime=startTime;
	private long thisTime=startTime;
	private void printTiming(String title){
		this.thisTime=System.nanoTime();
		System.out.println(title+ " calculated at "+IJ.d2s(0.000000001*(this.thisTime-this.startTime),3)+
				" (+"+IJ.d2s(0.000000001*(this.thisTime-this.lastTime),3)+") sec");
		this.lastTime=this.thisTime;
	}
	private void printTimingInit(){
		this.startTime=System.nanoTime();
		this.lastTime=this.startTime;
		this.thisTime=this.startTime;
	}
	public LaserPointer(
			double headLasersTilt, // degrees, right laser lower than left laser
			double minimalIntensity,
			double maximalIntensity,
			int    overexposedRadius,
			double lowpassSigma,    // low pass sigma, in pixels
			double highpassSigma,
			double headLowpassSigma,    // low pass sigma, in pixels for optical head lasers
			double quadraticScaleSigma, // find local maximum by quadratic intrepolating pixels around maximal value (relative to lkow pass sigma)
			int    algorithmNumber,
			int    closestOffender,
			int    fartherstOffender,
			double fatZero,
			double greenFloor,      // when dividing by green, add this fraction of maximal value (decrease green accordingly)
			boolean useOther, // when true - use red and other color, when false - only red
			boolean otherGreen, // other color is green (false - blue)
			double threshold,
			boolean swapUV, // first
			boolean flipU,
			boolean flipV,
			boolean whiteOnly, // verify laser is on the white pattern cell
			double  maxOffsetFromCenter, // maximal offset of the laser spot from the center (<0.5)
			double [][] laserUVMap, // first index - number of pointer points
			double laserSignalToNoise, // Minimal signal-to-noise ratio for laser pointers
			double localMaxRadius, // sensor pix. currently uses just square (2*localMaxRadius+1)**2
			boolean usePatternFilter, // Filter laser positions by likely pattern white cells
			int decimatePatternFilter,
			double localContrastSigma,
			double localToGlobalContrast,
			double patternLowPassSigma,
			double patternThreshold,
			double maximalCellSize,
			int numPasses, // number of black/white alternations of the surrounding cells to use in quadrant filtering
			boolean bordersOK, // frame boreder as good cell for quadrant filter
			double blurredMaskThreshold, // select only areas with multiple pattern white cells
			double maskGrow, // grow final mask (pixels)
			int    debugLevel
			) {
		this.headLasersTilt=headLasersTilt; // degrees, right laser lower than left laser
		this.minimalIntensity=minimalIntensity;
		this.maximalIntensity=maximalIntensity;
		this.overexposedRadius=overexposedRadius;
		this.lowpassSigma=lowpassSigma;
		this.highpassSigma=highpassSigma;
		this.headLowpassSigma=headLowpassSigma;    // low pass sigma, in pixels for optical head lasers
		this.quadraticScaleSigma= quadraticScaleSigma; // find local maximum by quadratic intrepolating pixels around maximal value (relative to lkow pass sigma)

		this.algorithmNumber=algorithmNumber;
		this.closestOffender=closestOffender;
		this.fartherstOffender=fartherstOffender;
		this.fatZero=fatZero;

		this.greenFloor=  greenFloor;
		this.useOther=useOther;
		this.otherGreen=otherGreen;
		this.threshold=   threshold;
		this.swapUV=swapUV; // first
		this.flipU=flipU;
		this.flipV=flipV;
		this.whiteOnly=whiteOnly;
		this.maxOffsetFromCenter=maxOffsetFromCenter;
		this.laserUVMap=new double[laserUVMap.length][2];
		for (int i=0;i<laserUVMap.length;i++) {
			this.laserUVMap[i][0]=laserUVMap[i][0];
			this.laserUVMap[i][1]=laserUVMap[i][1];
		}
		this.laserSignalToNoise=laserSignalToNoise; // Minimal signal-to-noise ratio for laser pointers
		this.localMaxRadius=localMaxRadius; // sensor pix. currently uses just square (2*localMaxRadius+1)**2
		this.usePatternFilter=usePatternFilter; // Filter laser positions by likely pattern white cells
		this.decimatePatternFilter=decimatePatternFilter;
		this.localContrastSigma=localContrastSigma;
		this.localToGlobalContrast=localToGlobalContrast;
		this.patternLowPassSigma=patternLowPassSigma;
		this.patternThreshold=patternThreshold;
		this.maximalCellSize=maximalCellSize;
		this.numPasses=numPasses; // number of black/white alternations of the surrounding cells to use in quadrant filtering
		this.bordersOK=bordersOK; // frame boreder as good cell for quadrant filter
		this.blurredMaskThreshold=blurredMaskThreshold; // select only areas with multiple pattern white cells
		this.maskGrow=maskGrow; // grow final mask (pixels)
		this.debugLevel=debugLevel;
	}
	public int getNumberOfPointers(){
		return this.laserUVMap.length;
	}
	@Override
	public LaserPointer clone(){
		return new LaserPointer(
				this.headLasersTilt,// degrees, right laser lower than left laser
				this.minimalIntensity,
				this.maximalIntensity,
				this.overexposedRadius,
				this.lowpassSigma,    // low pass sigma, in pixels
				this.highpassSigma,
				this.quadraticScaleSigma, // find local maximum by quadratic intrepolating pixels around maximal value (relative to lkow pass sigma)
				this.headLowpassSigma,
				this.algorithmNumber,
				this.closestOffender,
				this.fartherstOffender,
				this.fatZero,
				this.greenFloor,      // when dividing by green, add this fraction of maximal value (decrease green accordingly)
				this.useOther,
				this.otherGreen,
				this.threshold,
				this.swapUV,
				this.flipU,
				this.flipV,
				this.whiteOnly,
				this.maxOffsetFromCenter,
				this.laserUVMap, // first index - number of pointer points
				this.laserSignalToNoise, // Minimal signal-to-noise ratio for laser pointers
				this.localMaxRadius, // sensor pix. currently uses just square (2*localMaxRadius+1)**2
				this.usePatternFilter, // Filter laser positions by likely pattern white cells
				this.decimatePatternFilter,
				this.localContrastSigma,
				this.localToGlobalContrast,
				this.patternLowPassSigma,
				this.patternThreshold,
				this.maximalCellSize,
				this.numPasses, // number of black/white alternations of the surrounding cells to use in quadrant filtering
				this.bordersOK,
				this.blurredMaskThreshold,// select only areas with multiple pattern white cells
				this.maskGrow,  // grow final mask (pixels)
				this.debugLevel
				);
	}
	public void setProperties(String prefix,Properties properties){
		properties.setProperty(prefix+"headLasersTilt",this.headLasersTilt+"");
		properties.setProperty(prefix+"minimalIntensity",this.minimalIntensity+"");
		properties.setProperty(prefix+"maximalIntensity",this.maximalIntensity+"");
		properties.setProperty(prefix+"overexposedRadius",this.overexposedRadius+"");
		properties.setProperty(prefix+"lowpassSigma",this.lowpassSigma+"");
		properties.setProperty(prefix+"highpassSigma",this.highpassSigma+"");
		properties.setProperty(prefix+"headLowpassSigma",this.headLowpassSigma+"");
		properties.setProperty(prefix+"quadraticScaleSigma",this.quadraticScaleSigma+"");
		properties.setProperty(prefix+"algorithmNumber",this.algorithmNumber+"");
		properties.setProperty(prefix+"closestOffender",this.closestOffender+"");
		properties.setProperty(prefix+"fartherstOffender",this.fartherstOffender+"");
		properties.setProperty(prefix+"fatZero",this.fatZero+"");
		properties.setProperty(prefix+"greenFloor",this.greenFloor+"");
		properties.setProperty(prefix+"useOther",this.useOther+"");
		properties.setProperty(prefix+"otherGreen",this.otherGreen+"");
		properties.setProperty(prefix+"threshold",this.threshold+"");
		properties.setProperty(prefix+"swapUV",this.swapUV+"");
		properties.setProperty(prefix+"flipU",this.flipU+"");
		properties.setProperty(prefix+"flipV",this.flipV+"");
		properties.setProperty(prefix+"whiteOnly",this.whiteOnly+"");
		properties.setProperty(prefix+"maxOffsetFromCenter",this.maxOffsetFromCenter+"");
		properties.setProperty(prefix+"numberOfLaserPoints",this.laserUVMap.length+"");
		for (int i=0;i<this.laserUVMap.length;i++) {
			properties.setProperty(prefix+"laserUVMap_"+i+"u",this.laserUVMap[i][0]+"");
			properties.setProperty(prefix+"laserUVMap_"+i+"v",this.laserUVMap[i][1]+"");
		}
		properties.setProperty(prefix+"laserSignalToNoise",this.laserSignalToNoise+"");
		properties.setProperty(prefix+"localMaxRadius",this.localMaxRadius+"");
		properties.setProperty(prefix+"usePatternFilter",this.usePatternFilter+"");
		properties.setProperty(prefix+"decimatePatternFilter",this.decimatePatternFilter+"");
		properties.setProperty(prefix+"localContrastSigma",this.localContrastSigma+"");
		properties.setProperty(prefix+"localToGlobalContrast",this.localToGlobalContrast+"");
		properties.setProperty(prefix+"patternLowPassSigma",this.patternLowPassSigma+"");
		properties.setProperty(prefix+"patternThreshold",this.patternThreshold+"");
		properties.setProperty(prefix+"maximalCellSize",this.maximalCellSize+"");
		properties.setProperty(prefix+"numPasses",this.numPasses+"");
		properties.setProperty(prefix+"bordersOK",this.bordersOK+"");
		properties.setProperty(prefix+"blurredMaskThreshold",this.blurredMaskThreshold+"");
		properties.setProperty(prefix+"maskGrow",this.maskGrow+"");
		properties.setProperty(prefix+"debugLevel",this.debugLevel+"");
	}
	public void getProperties(String prefix,Properties properties){
		int numberOfLaserPoints=0;
		if (properties.getProperty(prefix+"headLasersTilt")!=null)
			this.headLasersTilt=Double.parseDouble(properties.getProperty(prefix+"headLasersTilt"));
		if (properties.getProperty(prefix+"minimalIntensity")!=null)
			this.minimalIntensity=Double.parseDouble(properties.getProperty(prefix+"minimalIntensity"));
		if (properties.getProperty(prefix+"maximalIntensity")!=null)
			this.maximalIntensity=Double.parseDouble(properties.getProperty(prefix+"maximalIntensity"));
		if (properties.getProperty(prefix+"overexposedRadius")!=null)
			this.overexposedRadius=Integer.parseInt(properties.getProperty(prefix+"overexposedRadius"));
		if (properties.getProperty(prefix+"lowpassSigma")!=null)
			this.lowpassSigma=Double.parseDouble(properties.getProperty(prefix+"lowpassSigma"));
		if (properties.getProperty(prefix+"highpassSigma")!=null)
			this.highpassSigma=Double.parseDouble(properties.getProperty(prefix+"highpassSigma"));
		if (properties.getProperty(prefix+"headLowpassSigma")!=null)
			this.headLowpassSigma=Double.parseDouble(properties.getProperty(prefix+"headLowpassSigma"));
		if (properties.getProperty(prefix+"quadraticScaleSigma")!=null)
			this.quadraticScaleSigma=Double.parseDouble(properties.getProperty(prefix+"quadraticScaleSigma"));
		if (properties.getProperty(prefix+"algorithmNumber")!=null)
			this.algorithmNumber=Integer.parseInt(properties.getProperty(prefix+"algorithmNumber"));
		if (properties.getProperty(prefix+"closestOffender")!=null)
			this.closestOffender=Integer.parseInt(properties.getProperty(prefix+"closestOffender"));
		if (properties.getProperty(prefix+"fartherstOffender")!=null)
			this.fartherstOffender=Integer.parseInt(properties.getProperty(prefix+"fartherstOffender"));
		if (properties.getProperty(prefix+"fatZero")!=null)
			this.fatZero=Double.parseDouble(properties.getProperty(prefix+"fatZero"));
		if (properties.getProperty(prefix+"greenFloor")!=null)
			this.greenFloor=Double.parseDouble(properties.getProperty(prefix+"greenFloor"));
		if (properties.getProperty(prefix+"useOther")!=null)
			this.useOther=Boolean.parseBoolean(properties.getProperty(prefix+"useOther"));
		if (properties.getProperty(prefix+"otherGreen")!=null)
			this.otherGreen=Boolean.parseBoolean(properties.getProperty(prefix+"otherGreen"));
		if (properties.getProperty(prefix+"threshold")!=null)
			this.threshold=Double.parseDouble(properties.getProperty(prefix+"threshold"));
		if (properties.getProperty(prefix+"swapUV")!=null)
			this.swapUV=Boolean.parseBoolean(properties.getProperty(prefix+"swapUV"));
		if (properties.getProperty(prefix+"flipU")!=null)
			this.flipU=Boolean.parseBoolean(properties.getProperty(prefix+"flipU"));
		if (properties.getProperty(prefix+"flipV")!=null)
			this.flipV=Boolean.parseBoolean(properties.getProperty(prefix+"flipV"));
		if (properties.getProperty(prefix+"whiteOnly")!=null)
			this.whiteOnly=Boolean.parseBoolean(properties.getProperty(prefix+"whiteOnly"));
		if (properties.getProperty(prefix+"maxOffsetFromCenter")!=null)
			this.maxOffsetFromCenter=Double.parseDouble(properties.getProperty(prefix+"maxOffsetFromCenter"));
		if (properties.getProperty(prefix+"numberOfLaserPoints")!=null) {
			numberOfLaserPoints=Integer.parseInt(properties.getProperty(prefix+"numberOfLaserPoints"));
			this.laserUVMap=new double[numberOfLaserPoints][2];
			for (int i=0;i<this.laserUVMap.length;i++) {
				this.laserUVMap[i][0]=Double.parseDouble(properties.getProperty(prefix+"laserUVMap_"+i+"u"));
				this.laserUVMap[i][1]=Double.parseDouble(properties.getProperty(prefix+"laserUVMap_"+i+"v"));
			}
		}
		if (properties.getProperty(prefix+"laserSignalToNoise")!=null)
			this.laserSignalToNoise=Double.parseDouble(properties.getProperty(prefix+"laserSignalToNoise"));
		if (properties.getProperty(prefix+"localMaxRadius")!=null)
			this.localMaxRadius=Double.parseDouble(properties.getProperty(prefix+"localMaxRadius"));
		if (properties.getProperty(prefix+"usePatternFilter")!=null)
			this.usePatternFilter=Boolean.parseBoolean(properties.getProperty(prefix+"usePatternFilter"));
		if (properties.getProperty(prefix+"decimatePatternFilter")!=null)
			this.decimatePatternFilter=Integer.parseInt(properties.getProperty(prefix+"decimatePatternFilter"));
		if (properties.getProperty(prefix+"localContrastSigma")!=null)
			this.localContrastSigma=Double.parseDouble(properties.getProperty(prefix+"localContrastSigma"));
		if (properties.getProperty(prefix+"localToGlobalContrast")!=null)
			this.localToGlobalContrast=Double.parseDouble(properties.getProperty(prefix+"localToGlobalContrast"));
		if (properties.getProperty(prefix+"patternLowPassSigma")!=null)
			this.patternLowPassSigma=Double.parseDouble(properties.getProperty(prefix+"patternLowPassSigma"));
		if (properties.getProperty(prefix+"patternThreshold")!=null)
			this.patternThreshold=Double.parseDouble(properties.getProperty(prefix+"patternThreshold"));
		if (properties.getProperty(prefix+"maximalCellSize")!=null)
			this.maximalCellSize=Double.parseDouble(properties.getProperty(prefix+"maximalCellSize"));
		if (properties.getProperty(prefix+"numPasses")!=null)
			this.numPasses=Integer.parseInt(properties.getProperty(prefix+"numPasses"));
		if (properties.getProperty(prefix+"bordersOK")!=null)
			this.bordersOK=Boolean.parseBoolean(properties.getProperty(prefix+"bordersOK"));
		if (properties.getProperty(prefix+"blurredMaskThreshold")!=null)
			this.blurredMaskThreshold=Double.parseDouble(properties.getProperty(prefix+"blurredMaskThreshold"));
		if (properties.getProperty(prefix+"maskGrow")!=null)
			this.maskGrow=Double.parseDouble(properties.getProperty(prefix+"maskGrow"));
		if (properties.getProperty(prefix+"debugLevel")!=null)
			this.debugLevel=Integer.parseInt(properties.getProperty(prefix+"debugLevel"));


	}
	public boolean showFilterDialog(String title){
		GenericDialog gd = new GenericDialog(title);
		gd.addNumericField("Decimate image for Pattern filter", this.decimatePatternFilter, 0,1,"x");
		gd.addNumericField("Sigma to calculate local level and contrast", this.localContrastSigma, 2,5,"sensor pixels");
		gd.addNumericField("Local/global contrast normalization", 100*this.localToGlobalContrast, 1,5,"%");
		gd.addNumericField("Filter sigma to apply to the pattern before thresholding",       this.patternLowPassSigma, 1,5,"sensor pix");
		gd.addNumericField("Pattern cell threshold as a fraction of dispersion",             100*this.patternThreshold, 1,5,"%");
		gd.addNumericField("Maximal pattern cell size (for discrimination)",  this.maximalCellSize, 1,5,"pix");

		gd.addNumericField("Number of black/white alternations of the surrounding cells to use in quadrant filtering", this.numPasses, 0);
		gd.addCheckbox   ("Borders as good cells for quadrant filter of possible pattern",  this.bordersOK);
		gd.addNumericField("Blurred white pattern cells mask (to remove separate white cells)", 100*this.blurredMaskThreshold, 1,5,"%");
		gd.addNumericField("Grow final detected pattern white cells mask", this.maskGrow, 2,5,"sensor pixels");
		gd.addNumericField("Debug level",   this.debugLevel, 0);
		gd.showDialog();
		if (gd.wasCanceled()) return false;

		this.decimatePatternFilter= (int) gd.getNextNumber();
		this.localContrastSigma=          gd.getNextNumber();
		this.localToGlobalContrast=  0.01*gd.getNextNumber();
		this.patternLowPassSigma=         gd.getNextNumber();
		this.patternThreshold=       0.01*gd.getNextNumber();
		this.maximalCellSize=             gd.getNextNumber();
		this.numPasses=             (int) gd.getNextNumber();
		this.bordersOK=                   gd.getNextBoolean();
		this.blurredMaskThreshold=   0.01*gd.getNextNumber();
		this.maskGrow=                    gd.getNextNumber();
		this.debugLevel=            (int) gd.getNextNumber();
		return true;
	}



	public boolean showDialog(String title){
		return showDialog(title, -1);
	}
	public boolean showDialog(String title,
			int numberOfPoints) { // >0 - ask only UV, <=0 all, use same number of points
		GenericDialog gd = new GenericDialog(title);
		if (numberOfPoints<=0) {

			gd.addNumericField("Angle between 2 laser spots and horizontal (right lower than left - positive)", this.headLasersTilt, 3,7,"degrees");
			gd.addNumericField("Minimal intensity at the expected pointer as a fraction of saturation", 100*this.minimalIntensity, 1,5,"%");
			gd.addNumericField("Maximal expected intensity at pointer (when it is off as a fraction of scaled saturation)", 100*this.maximalIntensity, 1,5,"%");
			gd.addNumericField("Do not look for the pointer closer than this distance from overexposed areas", this.overexposedRadius, 0,5,"pix");

			gd.addNumericField("Target laser spot detection low pass filter (4 spots)",         this.lowpassSigma, 1,5,"pix");
			gd.addNumericField("Target spot detection high pass filter  (4 spots)",             this.highpassSigma, 1,5,"pix");

			gd.addNumericField("Optical head laser spot detection low pass filter  (2 spots)",  this.headLowpassSigma, 1,5,"pix");
			gd.addNumericField("Scale low pass sigma when quadratic interpolate for maximum (0 - no interpolation)",   this.quadraticScaleSigma, 1,5,"x");

			gd.addNumericField("Algorithm number to detect pointers ",                          this.algorithmNumber, 0,1,"");
			gd.addNumericField("Do not check for above-threshold closer to the current point",  this.closestOffender, 0,5,"pix");
			gd.addNumericField("Prohibit above-threshold points closer to each other than",     this.fartherstOffender, 0,5,"pix");
			gd.addNumericField("Fat zero for combining differences",                            this.fatZero, 3,5,"");

			gd.addNumericField("Normalization to green, floor (100% - no normalization)", 100.0*this.greenFloor,  1,5,"%");

			gd.addCheckbox("Compare red to other color (green or blue)", this.useOther);
			gd.addCheckbox("If compare, compare to green (unchecked - blue)", this.otherGreen);

			gd.addNumericField("Red/Green difference to R/G average to be a laser spot",  100.0*this.threshold,  1,7,"%");
			gd.addMessage("Default grid orientation, used if not enough pointers are visible (auto-modified when more appear)");
			gd.addCheckbox("Swap U avd V",this.swapUV); // first
			gd.addCheckbox("Flip U direction",this.flipU);
			gd.addCheckbox("Flip V direction",this.flipV);

			gd.addCheckbox("Allow laser pointers on the white cells only",this.whiteOnly);
			gd.addNumericField("Maximal relative distance of the laser spot from the pattern cell center ",  100.0*this.maxOffsetFromCenter,  1,7,"%");
			gd.addNumericField("Minimal pointer S/N ratio", this.laserSignalToNoise, 2,5,"x");
			gd.addNumericField("Radius for finding local maximums", this.localMaxRadius, 2,5,"sensor pixels");
			gd.addMessage("==== Pattern Filter Parameters ====");
			gd.addCheckbox    ("Use pattern filter for laser pointer detection (and apply next settings)",  this.usePatternFilter);
			gd.addNumericField("Decimate image for Pattern filter", this.decimatePatternFilter, 0,1,"x");
			gd.addNumericField("Sigma to calculate local level and contrast", this.localContrastSigma, 2,5,"sensor pixels");
			gd.addNumericField("Local/global contrast normalization", 100*this.localToGlobalContrast, 1,5,"%");
			gd.addNumericField("Filter sigma to apply to the pattern before thresholding",       this.patternLowPassSigma, 1,5,"sensor pix");
			gd.addNumericField("Pattern cell threshold as a fraction of dispersion",             100*this.patternThreshold, 1,5,"%");
			gd.addNumericField("Maximal pattern cell size (for discrimination)",  this.maximalCellSize, 1,5,"pix");
			gd.addNumericField("Number of black/white alternations of the surrounding cells to use in quadrant filtering", this.numPasses, 0);
			gd.addCheckbox    ("Borders as good cells for quadrant filter of possible pattern",  this.bordersOK);
			gd.addNumericField("Blurred white pattern cells mask (to remove separate white cells)", 100*this.blurredMaskThreshold, 1,5,"%");
			gd.addNumericField("Grow final detected pattern white cells mask", this.maskGrow, 2,5,"sensor pixels");
			gd.addNumericField("Pattern filter debug level",   this.debugLevel, 0);

		} else {
			double [][] newLaserUVMap = new double [numberOfPoints][2];
			for (int i=0; i<numberOfPoints;i++) {
				newLaserUVMap[i][0]=(i>(laserUVMap.length-1))?0.0:this.laserUVMap[i][0];
				newLaserUVMap[i][1]=(i>(laserUVMap.length-1))?0.0:this.laserUVMap[i][1];
			}
			this.laserUVMap=newLaserUVMap;
		}
		for (int i=0;i<this.laserUVMap.length;i++) {
			gd.addMessage("Laser point "+i+" grid coordinates:");
			gd.addNumericField("Grid U "+i,                   this.laserUVMap[i][0], 2,7,"grid periods");
			gd.addNumericField("Grid V "+i,                   this.laserUVMap[i][1], 2,7,"grid periods");
		}
		if (numberOfPoints<=0) {
			gd.addNumericField("Number of laser points (will ask  for U V if modified)",  this.laserUVMap.length,  0,2,"");
		}

		WindowTools.addScrollBars(gd);
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		if (numberOfPoints<=0) {
			this.headLasersTilt=gd.getNextNumber();
			this.minimalIntensity=   0.01*gd.getNextNumber();
			this.maximalIntensity=   0.01*gd.getNextNumber();
			this.overexposedRadius= (int) gd.getNextNumber();
			this.lowpassSigma=            gd.getNextNumber();
			this.highpassSigma=           gd.getNextNumber();
			this.headLowpassSigma=        gd.getNextNumber();
			this.quadraticScaleSigma=     gd.getNextNumber();
			this.algorithmNumber=   (int) gd.getNextNumber();
			this.closestOffender=   (int) gd.getNextNumber();
			this.fartherstOffender=  (int) gd.getNextNumber();
			this.fatZero=                 gd.getNextNumber();
			this.greenFloor=         0.01*gd.getNextNumber();
			this.useOther=                gd.getNextBoolean();
			this.otherGreen=              gd.getNextBoolean();
			this.threshold=          0.01*gd.getNextNumber();
			this.swapUV=                  gd.getNextBoolean();
			this.flipU=                   gd.getNextBoolean();
			this.flipV=                   gd.getNextBoolean();
			this.whiteOnly=               gd.getNextBoolean();
			this.maxOffsetFromCenter=0.01*gd.getNextNumber();
			this.laserSignalToNoise=      gd.getNextNumber();
			this.localMaxRadius=          gd.getNextNumber();
			this.usePatternFilter=            gd.getNextBoolean();
			this.decimatePatternFilter= (int) gd.getNextNumber();
			this.localContrastSigma=          gd.getNextNumber();
			this.localToGlobalContrast=  0.01*gd.getNextNumber();
			this.patternLowPassSigma=         gd.getNextNumber();
			this.patternThreshold=       0.01*gd.getNextNumber();
			this.maximalCellSize=             gd.getNextNumber();
			this.numPasses=             (int) gd.getNextNumber();
			this.bordersOK=                   gd.getNextBoolean();
			this.blurredMaskThreshold=   0.01*gd.getNextNumber();
			this.maskGrow=                    gd.getNextNumber();
			this.debugLevel=            (int) gd.getNextNumber();
		}
		for (int i=0;i<this.laserUVMap.length;i++) {
			this.laserUVMap[i][0]=     gd.getNextNumber();
			this.laserUVMap[i][1]=     gd.getNextNumber();
		}
		if (numberOfPoints<=0) {
			numberOfPoints= (int) gd.getNextNumber();
			if ((numberOfPoints > 0) && (numberOfPoints!=this.laserUVMap.length)) showDialog(title,numberOfPoints);
		}
		return true;
	}
	/*
				ponterXY=laserPointers.laserPointer.getPointerXY( // returns x,y pair or null if pointer not detected
						greens,        // combined Bayer greens for each image, starting with no-laser
						this.laserPointers.laserWasOn(nPointer), // array specifying which image should have pointer on
						imp_pointed.getWidth(),// image width in pixels
		    			imp_pointed.getTitle()+"-"+nPointer,             // String title,
						this.debugLevel             // debug level (normal == 1)
				);

	 */
	public boolean [] localMaximum(
			double [] pixels,
			int width,
			int radius,
			int debugLevel){
		int height=pixels.length/width;
		ShowDoubleFloatArrays sdfra_instance= null;
		if (debugLevel>1) {
			sdfra_instance= new ShowDoubleFloatArrays(); // just for debugging?
		}

		boolean [] bmax=new boolean[pixels.length];
		// horizontal pass
		double []hmax = new double [pixels.length];
		for (int i=0;i<height;i++){
			int j0,j1;
			double max=0.0;
			for (int j=0;j<width;j++){
				j0=j-radius;
				if (j0<0) j0=0;
				j1=j+radius;
				if (j1>=width) j1=width-1;
				if (j>0) {
					max=Math.max(max,pixels[width*i+j1]);
					if ((j0>0) && (pixels[width*i+j0-1]<max)) {
						hmax[i*width+j]=max;
						continue; // first (to be removed) pixel was not max
					}
				}
				max=pixels[width*i+j0];
				for (int k=width*i+j0+1;k<=width*i+j1;k++) if (pixels[k]>max) max= pixels[k];
				hmax[i*width+j]=max;
			}
		}
		if (debugLevel>2) sdfra_instance.showArrays(hmax,   width, height, "hmax-"+radius);

		//vertical pass
		for (int j=0;j<width;j++){
			int i0,i1;
			double max=0.0;
			for (int i=0;i<height;i++){
				int index=i*width+j;
				i0=i-radius;
				if (i0<0) i0=0;
				i1=i+radius;
				if (i1>=height) i1=height-1;
				if (i>0) {
					max=Math.max(max,hmax[width*i1+j]);
					if ((i0>0) && (hmax[width*(i0-1)+j]<max)) {
						bmax[index]= (pixels[index]==max);
						continue; // first (to be removed) pixel was not max
					}
				}
				max=pixels[width*i0+j];
				for (int k=width*i0+j+width;k<=width*i1+j;k+=width) if (hmax[k]>max) max= hmax[k];
				bmax[index]= (pixels[index]==max);
			}
		}
		return bmax;
	}

	public boolean [] getPatternMask(
			double [] pixels,
			int width
			){
		ShowDoubleFloatArrays sdfra_instance= null;
		if (this.debugLevel>1) sdfra_instance= new ShowDoubleFloatArrays(); // just for debugging?
		DoubleGaussianBlur gb=new DoubleGaussianBlur();
		double initialScale=0.5;
		int height=pixels.length/width;
		double scale=initialScale/this.decimatePatternFilter;
		double [] dpixels;
		int dheight,dwidth;
		int extraWidth=0,extraHeight=0;
		double k=1.0/(this.decimatePatternFilter*this.decimatePatternFilter);
		double k1=0.0,k2=0.0,k3=0.0;
		dwidth=width/this.decimatePatternFilter;
		if (dwidth*this.decimatePatternFilter<width){
			extraWidth=width-(dwidth*this.decimatePatternFilter);
			dwidth++;
			k1=1.0/(this.decimatePatternFilter*extraWidth);
		}
		dheight=height/this.decimatePatternFilter;
		if (dheight*this.decimatePatternFilter<height) {
			extraHeight=height- (dheight*this.decimatePatternFilter);
			dheight++;
			k2=1.0/(this.decimatePatternFilter*extraHeight);
		}
		if ((dheight>0) && (dwidth>0)) k3=1.0/(extraWidth*extraHeight);
		double d;
		// decimate original image (to speedup calculations)
		if (this.decimatePatternFilter>1){
			int i,j;
			dpixels=new double [dwidth*dheight];
			for (i=0;i<height/this.decimatePatternFilter;i++){
				for (j=0;j<width/this.decimatePatternFilter;j++){
					d=0;
					int index=this.decimatePatternFilter*(i*width+j);
					for (int m=0;m<this.decimatePatternFilter;m++){
						for (int n=0;n<this.decimatePatternFilter;n++){
							d+=pixels[index+n];
						}
						index+=width;
					}
					/*
						if (i*dwidth+j>=dpixels.length){
							System.out.println(
									"dpixels.length="+dpixels.length+"\n"+
									"width="+width+"\n"+
									"height="+height+"\n"+
									"dwidth="+dwidth+"\n"+
									"dheight="+dheight+"\n"+
									"i="+i+"\n"+
									"j="+j+"\n");
						}
					 */
					dpixels[i*dwidth+j]=d*k;
				}
				if (extraWidth>0){ // j==dWidth-1
					d=0;
					int index=this.decimatePatternFilter*(i*width+j);
					for (int m=0;m<extraWidth;m++){
						for (int n=0;n<this.decimatePatternFilter;n++){
							d+=pixels[index+n];
						}
						index+=width;
					}
					dpixels[i*dwidth+j]=d*k1;
				}
			}
			if (extraHeight>0){ // i==dHeight-1
				for (j=0;j<width/this.decimatePatternFilter;j++){
					d=0;
					int index=this.decimatePatternFilter*(i*width+j);
					for (int m=0;m<extraHeight;m++){
						for (int n=0;n<this.decimatePatternFilter;n++){
							d+=pixels[index+n];
						}
						index+=width;
					}
					dpixels[i*dwidth+j]=d*k2;
				}
				if (extraWidth>0){ // j==dWidth-1
					d=0;
					int index=this.decimatePatternFilter*(i*width+j);
					for (int m=0;m<extraWidth;m++){
						for (int n=0;n<this.decimatePatternFilter;n++){
							d+=pixels[index+n];
						}
						index+=width;
					}
					dpixels[i*dwidth+j]=d*k3;
				}
			}
		} else dpixels=pixels.clone();
		if (this.debugLevel>2) sdfra_instance.showArrays(dpixels, dwidth, dheight,  "decimated");
		double [] dpixels_lp=dpixels.clone();
		gb.blurDouble(dpixels_lp, dwidth, dheight, scale*this.localContrastSigma, scale*this.localContrastSigma, 0.01);
		if (this.debugLevel>2) sdfra_instance.showArrays(dpixels_lp, dwidth, dheight,  "dpixels_lp");
		double sum=0.0;
		for (int i=0;i<dpixels.length;i++){
			dpixels[i]-=dpixels_lp[i];
			dpixels_lp[i]=dpixels[i]*dpixels[i];
			sum+=dpixels_lp[i];
		}
		double corr=(1.0-this.localToGlobalContrast)*Math.sqrt(sum/(dwidth*dheight));
		if (this.debugLevel>2) sdfra_instance.showArrays(dpixels_lp, dwidth, dheight,  "dpixels_var");
		gb.blurDouble(dpixels_lp, dwidth, dheight, scale*this.localContrastSigma, scale*this.localContrastSigma, 0.01);
		if (this.debugLevel>2) sdfra_instance.showArrays(dpixels_lp, dwidth, dheight,  "dpixels_var_blur");
		for (int i=0;i<dpixels.length;i++){
			dpixels_lp[i]=this.localToGlobalContrast*Math.sqrt(dpixels_lp[i])+corr;
			dpixels[i]/=dpixels_lp[i];
		}
		if (this.debugLevel>2) sdfra_instance.showArrays(dpixels_lp, dwidth, dheight,  "dpixels_denom");
		if (this.patternLowPassSigma>0) {
			if (this.debugLevel>2) sdfra_instance.showArrays(dpixels, dwidth, dheight,    "dpixels_normalized");
			gb.blurDouble(dpixels, dwidth, dheight, scale*this.patternLowPassSigma, scale*this.patternLowPassSigma, 0.01);
		}
		if (this.debugLevel>1) sdfra_instance.showArrays(dpixels, dwidth, dheight,    "dpixels_norm_blur");
		int [] ipixels=new int [dpixels.length];
		for (int i=0;i<dpixels.length;i++){
			if (dpixels[i]>this.patternThreshold) ipixels[i]=1;
			else if (dpixels[i]<-this.patternThreshold) ipixels[i]=-1;
			else ipixels[i]=0;
		}
		if (this.debugLevel>1)	sdfra_instance.showArrays(ipixels, dwidth, dheight,    "ipixels_threshold");

		int maxDist = (int) Math.round(scale*this.maximalCellSize);
		if (maxDist<1) maxDist=1;
		// TODO: save intermediate ipixels and use it with blurred version from more passes?
		for (int pass=this.numPasses-1;pass>=0;pass--) {
			filterQuadrant(
					ipixels,
					dwidth,
					(pass & 1)==0,  // false, // each black cell should have white in all 4 quadrants
					maxDist,
					this.bordersOK
					);
			if (this.debugLevel>2)	sdfra_instance.showArrays(ipixels, dwidth, dheight,    "pass-"+pass);
		}
		//blur-threshold-grow
		// convert to Double (white cells - 1.0, black and none - 0.0)
		for (int i=0;i<dpixels.length;i++) dpixels[i]= (ipixels[i]>0)?1.0:0.0;
		// blur result with sigma > cell period to find continuous pattern cell areas
		gb.blurDouble(dpixels, dwidth, dheight, scale*this.localContrastSigma, scale*this.localContrastSigma, 0.01);
		if (this.debugLevel>2)	sdfra_instance.showArrays(dpixels, dwidth, dheight,    "white_blurred");
		// Threshold to boolean mask
		boolean [] dmask=new boolean[dpixels.length];
		for (int i=0;i<dpixels.length;i++) dmask[i]= (dpixels[i]>=this.blurredMaskThreshold);
		if (this.debugLevel>1)	sdfra_instance.showArrays(dmask, dwidth, dheight,    "white_threshold");
		int iBlurMaskGrow= (int) Math.round(scale*this.localContrastSigma); // does in need separate coefficient?
		growMask( dmask, //boolean [] pixels,
				dwidth, // int width,
				iBlurMaskGrow); //int grow);
		if (this.debugLevel>1)	sdfra_instance.showArrays(dmask, dwidth, dheight,    "white_threshold_grown"+iBlurMaskGrow);
		// Mask out white cells outside of the compact areas just found
		for (int i=0;i<dpixels.length;i++) if (!dmask[i])  ipixels[i]=0;
		if (this.debugLevel>1)	sdfra_instance.showArrays(ipixels, dwidth, dheight,    "ipixels_masked");
		boolean [] mask=new boolean[pixels.length];
		if (this.decimatePatternFilter>1){
			for (int i=0;i<pixels.length;i++){
				int y= (i/width)/this.decimatePatternFilter;
				int x= (i%width)/this.decimatePatternFilter;
				mask[i]= (ipixels[x+dwidth*y]==1); // "good" white cells
			}

		} else {
			for (int i=0;i<pixels.length;i++) mask[i]= (ipixels[i]==1); // "good" white cells
		}
		int finalMaskGrow=(int) Math.round(this.maskGrow*initialScale);
		if (this.maskGrow>0.0) {
			growMask(mask, //boolean [] pixels,
					width, // int width,
					finalMaskGrow);// //int grow);
		}
		return mask;
	}
	void growMask(
			boolean [] pixels,
			int width,
			int grow){
		int height=pixels.length/width;
		int [] distLeft= new int [pixels.length];
		int [] distRight=new int [pixels.length]; // also used for down
		int [] distUp=  new int [pixels.length];
		int v1,v2;
		int index1,index2;
		int initialValue=width+height;
		for (int i=0;i<height;i++){
			v1=initialValue;
			v2=initialValue;
			index1=i*width;
			index2=index1+width-1;
			for (int j=0;j<width;j++){
				if (pixels[index1]) v1=0; else v1++;
				if (pixels[index2]) v2=0; else v2++;
				distLeft[index1++]=v1;
				distRight[index2--]=v2;
			}
		}
		// combine two (min distance)
		for (int i=0;i<pixels.length;i++) if (distRight[i]<distLeft[i])distLeft[i]=distRight[i];
		// Down
		for (int j=0;j<width;j++)distRight[j]= distLeft[j]; //  very top line (now used for down)
		index1=0;
		for (int i=width;i<pixels.length;i++){
			v1=distRight[index1++]+1;
			distRight[i]= (v1>distLeft[i])?distLeft[i]:v1;
		}
		// Up
		for (int j=pixels.length-width;j<pixels.length;j++) distUp[j]= distLeft[j]; // very bottom line
		index1=pixels.length-1;
		for (int i=pixels.length-1-width;i>=0;i--){
			v1=distUp[index1--]+1;
			distUp[i]= (v1>distLeft[i])?distLeft[i]:v1;
		}
		// combine two (min distance)
		for (int i=0;i<pixels.length;i++) if (distUp[i]<distRight[i])distRight[i]=distUp[i];
		// Now distRight contains the shortest distance from the nearest enabled pixels
		// update the original pixels to include tyhe new ones
		for (int i=0;i<pixels.length;i++) pixels[i] |= (distRight[i]<=grow);
	}

	void filterQuadrant(
			int [] pixels,
			int width,
			boolean fromBlack,
			int maxDist,
			boolean bordersOK
			){
		int height=pixels.length/width;
		ShowDoubleFloatArrays sdfra_instance= null;
		if (this.debugLevel>1) sdfra_instance= new ShowDoubleFloatArrays(); // just for debugging?

		int [] distLeft= new int [pixels.length];
		int [] distRight=new int [pixels.length]; // also used for down
		int [] distUp=  new int [pixels.length];
		int sign=fromBlack?-1:1;
		int initialValue=bordersOK?0:(width+height);
		int v1,v2;
		int index1,index2;
		for (int i=0;i<height;i++){
			v1=initialValue;
			v2=initialValue;
			index1=i*width;
			index2=index1+width-1;
			for (int j=0;j<width;j++){
				if (pixels[index1]==sign) v1=0; else v1++;
				if (pixels[index2]==sign) v2=0; else v2++;
				distLeft[index1++]=v1;
				distRight[index2--]=v2;
			}
		}
		if (this.debugLevel>3){
			sdfra_instance.showArrays(distLeft, width, height,    "distLeft-"+fromBlack);
			sdfra_instance.showArrays(distRight, width, height,   "distRight-"+fromBlack);

		}
		// combine two (max distance)
		for (int i=0;i<pixels.length;i++) if (distRight[i]>distLeft[i])distLeft[i]=distRight[i];
		if (this.debugLevel>3) sdfra_instance.showArrays(distLeft, width, height,    "comboLeftRight-"+fromBlack);


		// Down
		for (int j=0;j<width;j++)distRight[j]= bordersOK?0:distLeft[j]; // now used for down
		index1=0;
		for (int i=width;i<pixels.length;i++){
			//					v1=distLeft[index1++]+1;
			v1=distRight[index1++]+1;
			distRight[i]= (v1>distLeft[i])?distLeft[i]:v1;
		}
		if (this.debugLevel>3) sdfra_instance.showArrays(distRight, width, height,    "distDown-"+fromBlack);
		// Up
		for (int j=pixels.length-width;j<pixels.length;j++) distUp[j]= bordersOK?0:distLeft[j];
		index1=pixels.length-1;
		for (int i=pixels.length-1-width;i>=0;i--){
			//					v1=distLeft[index1--]+1;
			v1=distUp[index1--]+1;
			distUp[i]= (v1>distLeft[i])?distLeft[i]:v1;
		}
		if (this.debugLevel>3) sdfra_instance.showArrays(distUp, width, height,    "distUp-"+fromBlack);
		// combine two  (max distance)
		for (int i=0;i<pixels.length;i++) if (distUp[i]>distRight[i])distRight[i]=distUp[i];
		if (this.debugLevel>3) sdfra_instance.showArrays(distRight, width, height,    "combo-"+fromBlack);
		// Now distRight contains the longest of 4 quadrants distance from "good" pixels - remove bad opposite side ones
		for (int i=0;i<pixels.length;i++) if ((distRight[i]>maxDist) && (pixels[i]!=sign)) pixels[i]=0; // if opposite signe (or 0) make 0
	}


	public double [] getPointerXY( // returns x,y pair or null if pointer not detected
			double [][] backgroundBayer, // Bayer array of the background (lasers off) image 0,3 -G, 1-R,2-B
			double [][] pointedBayer,    // Bayer array of the (laser on) image
			int width,                   // image width in pixels
			boolean modBackground,       // modify background array (on the first pass)
			String title,
			int debugLevel   // debug level (normal == 1)
			){
		ShowDoubleFloatArrays sdfra_instance= null;
		if (debugLevel>1) sdfra_instance= new ShowDoubleFloatArrays(); // just for debugging?
		// As high precision is not needed we can map Bayer pixels to the same grid of half resolution of the image
		// 0,3 - green 1 - red (laser)
		int bayerG1=0;
		int bayerG2=3;
		int bayerR=1;
		double avrgGreenB=0.0;
		//	    		double avrgGreenP=0.0;
		int len=backgroundBayer[bayerG1].length;
		int halfWidth=width/2;
		int halfHeight=len/halfWidth;
		if (debugLevel>2){
			String subtitles[] ={"green1","red","blue","green2","combo","5"};
			sdfra_instance.showArrays(pointedBayer.clone(), halfWidth, halfHeight, true, title+"-bayer", subtitles);
		}

		for (int i=0; i<len;i++){
			if (modBackground) backgroundBayer[bayerG1][i]=0.5*(backgroundBayer[bayerG1][i]+backgroundBayer[bayerG2][i]);
			avrgGreenB+=backgroundBayer[bayerG1][i];
			pointedBayer[bayerG1][i]=0.5*(pointedBayer[bayerG1][i]+pointedBayer[bayerG2][i]);
			//	    			avrgGreenP+=pointedBayer[bayerG1][i];
		}
		avrgGreenB/=len;
		//	    		avrgGreenP/=len;
		for (int i=0; i<len;i++){
			if (modBackground) {
				backgroundBayer[bayerR][i]/=(backgroundBayer[bayerG1][i]*(1.0-this.greenFloor)+avrgGreenB*this.greenFloor);
			}
			pointedBayer[bayerR][i] = pointedBayer[bayerR][i]/
					(backgroundBayer[bayerG1][i]*(1.0-this.greenFloor)+avrgGreenB*this.greenFloor)-backgroundBayer[bayerR][i];

		}
		if (debugLevel>3) sdfra_instance.showArrays(pointedBayer[bayerR].clone(), halfWidth, halfHeight, title);
		// low pass filter, 2-d
		DoubleGaussianBlur gb=new DoubleGaussianBlur();
		gb.blurDouble(pointedBayer[bayerR], halfWidth, halfHeight, this.lowpassSigma, this.lowpassSigma, 0.01);
		if (debugLevel>2) sdfra_instance.showArrays(pointedBayer[bayerR].clone(), halfWidth, halfHeight, title+"-smooth");
		// Finding just maximum, to centroid here
		int indx=0;
		double max=pointedBayer[bayerR][indx];
		for (int i=0; i<len;i++){
			if (pointedBayer[bayerR][i]>max) {
				max=pointedBayer[bayerR][i];
				indx=i;
			}
		}
		double [] result={2*(indx%halfWidth)+1.0, 2*(indx/halfWidth)-1.0};
		if (debugLevel>1) System.out.println("Max="+max+"(>"+this.threshold+"?), x="+result[0]+", y="+result[1]);
		if (max>=this.threshold) return result;
		return null;
	}

	public double [][] getPointerXY( // returns x,y pair or null if pointer not detected
			boolean headLaserMode,
			double saturationRed, // maximal red intensity scaled to reduced exposure
			double scaleExposureForLasers,
			//	    			boolean skipFirst, // do not use no-pointer image (may have different exposure time)
			double [][] pre_greens, // combined Bayer greens for each image, starting with no-laser (maybe Blue or none)
			double [][] reds,   // red Bayer component for each image, starting with no-laser
			boolean [][] whichOn, // array specifying which image should have pointer on
			int width,                   // image width in pixels
			String title,
			int debugLevel   // debug level (normal == 1)
			){
		int debugTiming=1;
		if (debugLevel>debugTiming) printTimingInit();
		boolean skipFirst=scaleExposureForLasers>0.0; // do not use no-pointer image (may have different exposure time)
		///								if (scaleExposureForLasers>0) saturationRed*=scaleExposureForLasers; // scaled to reduced exposure time
		double noLaserMaxRed=(scaleExposureForLasers>0)?(saturationRed/scaleExposureForLasers):saturationRed;
		if (this.maximalIntensity<1.0) noLaserMaxRed*=this.maximalIntensity;


		double [][] greens=headLaserMode?null:pre_greens;
		int startIndex=skipFirst?1:0;
		ShowDoubleFloatArrays sdfra_instance= null;
		String [] subtitles_all= null;
		String [] subtitles= null;
		if (debugLevel>1) {
			sdfra_instance= new ShowDoubleFloatArrays(); // just for debugging?
			subtitles_all= new String [reds.length];
			for (int i=0;i<reds.length;i++){
				subtitles_all[i]="";
				for (int j=0;j<whichOn.length;j++) subtitles_all[i]+=whichOn[j][i]?"+":"-";
			}
			subtitles= new String [whichOn.length];
			for (int i=0;i<whichOn.length;i++){
				subtitles[i]="p-"+i;
			}

		}
		// As high precision is not needed we can map Bayer pixels to the same grid of half resolution of the image
		// 0,3 - green 1 - red (laser)
		double avrgGreen;
		double avrgRed;
		int len=reds[0].length;
		int halfWidth=width/2;
		int halfHeight=len/halfWidth;

		double minimalIntensity=saturationRed*this.minimalIntensity;
		double maximalIntensity=saturationRed*this.maximalIntensity;
		int [] thresholds = new int [reds[0].length];
		boolean [] overexposed = new boolean [len]; /* somewhat duplicates thersholds */
		for (int i=0;i<len;i++) overexposed[i]=reds[0][i] > noLaserMaxRed;
		for (int i=0;i<thresholds.length;i++) thresholds[i]=0;

		if (debugLevel>2){
			if (greens!=null) sdfra_instance.showArrays(greens, halfWidth, halfHeight, true, "green-"+title, subtitles_all);
			sdfra_instance.showArrays(reds,   halfWidth, halfHeight, true, "red-"+title, subtitles_all);
		}
		if (debugLevel>debugTiming) printTiming("red-"+title);
		boolean[]  patternMask=null;
		if (this.usePatternFilter){
			double [] ppixels=(greens!=null)?greens[0]:reds[0];
			patternMask=getPatternMask(	ppixels, halfWidth );
			if (debugLevel>2)	sdfra_instance.showArrays(patternMask, halfWidth, halfHeight,  "patternMask-"+title);
			if (debugLevel>debugTiming) printTiming("patternMask-"+title);
		}
		for (int i=startIndex;i<reds.length;i++) {
			avrgGreen=0.0;
			avrgRed=0.0;
			for (int j=0; j<len;j++) avrgRed+=  reds[i][j];
			if (greens!=null) for (int j=0; j<len;j++) avrgGreen+=greens[i][j];
			avrgGreen/=len;
			avrgRed/=len;
			if (greens!=null){
				for (int j=0; j<len;j++){
					if ((reds[i][j]<minimalIntensity)) thresholds[j] |= 1; // under
					if ((reds[i][j]>maximalIntensity)) thresholds[j] |= 2; // over
					reds[i][j]=reds[i][j]/avrgRed- greens[i][j]/avrgGreen*(1-this.greenFloor); // always
				}
			} else {
				for (int j=0; j<len;j++){
					if ((reds[i][j]<minimalIntensity)) thresholds[j] |= 1; // under
					if ((reds[i][j]>maximalIntensity)) thresholds[j] |= 2; // over
					reds[i][j]/=avrgRed;
				}
			}
		}
		DoubleGaussianBlur gb=new DoubleGaussianBlur();
		greens=null; // don't need it anymore
		if (debugLevel>2){
			System.out.println("saturationRed="+saturationRed+" minimalIntensity="+minimalIntensity+" maximalIntensity="+maximalIntensity);
			sdfra_instance.showArrays(reds,   halfWidth, halfHeight, true, "normalized-"+title, subtitles_all);
		}
		if (debugLevel>debugTiming) printTiming("normalized-"+title);
		double [] dThresholds=new double [len];
		for (int i=0;i<len;i++){
			dThresholds[i]=((thresholds[i] & 1)!=0)?(-1): ( ((thresholds[i] & 2)!=0)?1.0:0.0 );
		}
		if (debugLevel>2)	sdfra_instance.showArrays(dThresholds,   halfWidth, halfHeight, "intensity_limits-"+title);
		if (debugLevel>debugTiming) printTiming("intensity_limits-"+title);
		// high-pass images
		if (!headLaserMode && (this.highpassSigma>0)){
			for (int numImg=startIndex;numImg<reds.length;numImg++){
				double [] redBlured=reds[numImg].clone();
				gb.blurDouble(redBlured, halfWidth, halfHeight, this.highpassSigma, this.highpassSigma, 0.01);
				for (int i=0;i<redBlured.length;i++)reds[numImg][i]-=redBlured[i];
			}
		}
		if (debugLevel>2)sdfra_instance.showArrays(reds,   halfWidth, halfHeight, true, "highpass-"+title, subtitles_all);
		if (debugLevel>debugTiming) printTiming("highpass-"+title);
		// low pass filter, 2-d
		double threshold=headLaserMode?this.threshold:this.threshold; // so far - the same?
		double lpSigma=headLaserMode?this.headLowpassSigma:this.lowpassSigma;

		double [][] diffOnOff=new double[whichOn.length][len];
		double [][] noise=null; // used in algorithm 4
		//	    		int debugX=952,debugY=665,debugAround=3,debugPointer=1;

		//	    		if (headLaserMode || (this.algorithmNumber==0) || (this.algorithmNumber==2)){
		if (headLaserMode || (this.algorithmNumber==0)){
			for (int numPointer=0; numPointer <whichOn.length;numPointer++) for (int i=0;i<len;i++){
				double on=-1.0,off=-1.0;
				for (int nImg=startIndex;nImg<whichOn[numPointer].length; nImg++){
					if (whichOn[numPointer][nImg]){
						if ((on<0)  || (on>reds[nImg][i])) on=reds[nImg][i]; // smallest among "on"
					} else {
						if ((off<0) || (off<reds[nImg][i])) off=reds[nImg][i]; // largest among "off"
					}
				}
				diffOnOff[numPointer][i]=on-off; // difference between lowest "on" and largest "off"
			}
		} else if (this.algorithmNumber==1){
			for (int numPointer=0; numPointer <whichOn.length;numPointer++) {
				//	    				int numSamples=0;
				int numPositiveSamples=0;
				for (int nImg=startIndex;nImg<whichOn[numPointer].length; nImg++){
					//	    					numSamples++;
					if (whichOn[numPointer][nImg]) numPositiveSamples++;
				}
				for (int i=0;i<len;i++) {
					diffOnOff[numPointer][i]=1.0;
					for (int nImg=startIndex;nImg<whichOn[numPointer].length; nImg++) {
						if (whichOn[numPointer][nImg]) diffOnOff[numPointer][i] +=reds[nImg][i];
						else  diffOnOff[numPointer][i] -=reds[nImg][i];
					}
					diffOnOff[numPointer][i]/=numPositiveSamples;
				}
			}

		} else 	if (this.algorithmNumber==2){
			for (int numPointer=0; numPointer <whichOn.length;numPointer++) {
				int numSamples=0;
				//	    				int numPositiveSamples=0;
				for (int nImg=startIndex;nImg<whichOn[numPointer].length; nImg++){
					numSamples++;
					//	    					if (whichOn[numPointer][nImg]) numPositiveSamples++;
				}
				double pwr=1.0/numSamples;
				for (int i=0;i<len;i++) {
					double on=-1.0,off=-1.0;
					for (int nImg=startIndex;nImg<whichOn[numPointer].length; nImg++){
						if (whichOn[numPointer][nImg]){
							if ((on<0)  || (on>reds[nImg][i])) on=reds[nImg][i]; // smallest among "on"
						} else {
							if ((off<0) || (off<reds[nImg][i])) off=reds[nImg][i]; // largest among "off"
						}
					}
					if (on<off){
						diffOnOff[numPointer][i]=0.0;
						continue;
					}
					double average=0.5*(on+off); // middle between highest "off" and lowest "on"
					diffOnOff[numPointer][i]=1.0;
					for (int nImg=startIndex;nImg<whichOn[numPointer].length; nImg++) {
						double diff=reds[nImg][i]-average+this.fatZero;
						if (!whichOn[numPointer][nImg]) diff=-diff;
						if (diff>0.0) diffOnOff[numPointer][i]*=diff;
						else {
							diffOnOff[numPointer][i]=0;
							break;
						}
					}
					if (diffOnOff[numPointer][i]>0.0) {
						diffOnOff[numPointer][i]=Math.pow(diffOnOff[numPointer][i],pwr);
					}
					diffOnOff[numPointer][i]-=this.fatZero;
				}
			}
		} else 	if (this.algorithmNumber==3){
			for (int numPointer=0; numPointer <whichOn.length;numPointer++){
				for (int i=0;i<len;i++){
					double d=0.0;
					for (int nImg=1;nImg<whichOn[numPointer].length; nImg++){ // skip first image with all 0ff
						if (whichOn[numPointer][nImg]) d+=reds[nImg][i];
						else                          d-=reds[nImg][i];
					}
					diffOnOff[numPointer][i]=d/(whichOn[numPointer].length-1);
				}
			}
		} else 	if (this.algorithmNumber==4){ // sometimes long - 11.391 sec, other - 0.558 s
			double avOn,avOff;
			noise=new double [whichOn.length][];
			if (debugLevel>debugTiming) printTiming("alg_"+this.algorithmNumber+"-start-"+title);
			for (int numPointer=0; numPointer <whichOn.length;numPointer++){
				if (debugLevel>debugTiming) printTiming("numPointer-"+numPointer+"-"+whichOn[numPointer].length+"-"+len+"-"+title);
				noise[numPointer]=new double[len];
				for (int i=0;i<len;i++){
					//	    					double d=0.0;
					avOn= 0.0;
					avOff=0.0;
					for (int nImg=1;nImg<whichOn[numPointer].length; nImg++){
						if (whichOn[numPointer][nImg]) {
							//	    							d+=reds[nImg][i];
							avOn+=reds[nImg][i];
						} else {
							//	    							d-=reds[nImg][i];
							avOff+=reds[nImg][i];
						}
					}
					double d=avOn-avOff;
					diffOnOff[numPointer][i]=d/(whichOn[numPointer].length-1);
					avOn/=((whichOn[numPointer].length-1)/2);
					avOff/=((whichOn[numPointer].length-1)/2);
					double s2=0.0;
					for (int nImg=1;nImg<whichOn[numPointer].length; nImg++){
						d=reds[nImg][i]-((whichOn[numPointer][nImg])?avOn:avOff);
						s2+=d*d;
					}
					//	    					s2n[numPointer][i]=diffOnOff[numPointer][i]/(Math.sqrt(s2/(whichOn[numPointer].length-1)));
					noise[numPointer][i]=(Math.sqrt(s2/(whichOn[numPointer].length-1)));
				}
			}
		}

		// mask out too dim/too bright - TODO: still need to grow overexposed mask
		/* Mask out pixels closer than overexposedRadius from overesposed areas*/
		//overexposedRadius
		if (debugLevel>debugTiming) printTiming("alg_"+this.algorithmNumber+"-done-"+title);
		if (!headLaserMode) {
			int discardBorder=5;
			if (this.overexposedRadius>0){
				boolean [] overexposedTmp;
				for (int n=0;n<this.overexposedRadius;n++){
					overexposedTmp=overexposed;
					overexposed=new boolean[len];
					int i=0;
					//	    					System.out.println("halfWidth="+halfWidth+"  halfHeight="+halfHeight+"  len="+len);
					for (int iy=0;iy<halfHeight;iy++) for (int ix=0;ix<halfWidth;ix++){
						/*
    						if (i>=len){
    	    					System.out.println("  iy="+iy+"  ix="+ix+"  i="+i);
    							break;
    						}
    						overexposed[iy*halfWidth+ix]  = ((iy>0) && overexposedTmp[i-halfWidth]);
    						overexposed[iy*halfWidth+ix] |= ((ix>0) && overexposedTmp[i-1]);
    						overexposed[iy*halfWidth+ix] |= ((iy<(halfHeight-1)) && overexposedTmp[i+halfWidth]);
    						overexposed[iy*halfWidth+ix] |= ((ix<(halfWidth-1)) && overexposedTmp[i+1]);


    						overexposed[i++]=
    								((iy>0) && overexposedTmp[i-halfWidth]) ||
    								((ix>0) && overexposedTmp[i-1]) ||
    								((iy<(halfHeight-1)) && overexposedTmp[i+halfWidth]) ||
    								((ix<(halfWidth-1)) && overexposedTmp[i+1]);
						 */
						overexposed[i]  = overexposedTmp[i];
						overexposed[i] |= ((iy>0) && overexposedTmp[i-halfWidth]);
						overexposed[i] |= ((ix>0) && overexposedTmp[i-1]);
						overexposed[i] |= ((iy<(halfHeight-1)) && overexposedTmp[i+halfWidth]);
						overexposed[i] |= ((ix<(halfWidth-1)) && overexposedTmp[i+1]);
						i++;
					}
				}
			}
			if (debugLevel>debugTiming) printTiming("grown_by_"+this.overexposedRadius+"-"+title);
			int debugNumBorder=0;
			for (int iy=0;iy<halfHeight;iy++) for (int ix=0;ix<halfWidth;ix++){
				if ((iy<discardBorder) || (ix<discardBorder) || (ix>=(halfWidth-discardBorder)) || (iy>=(halfHeight-discardBorder))){
					overexposed[iy*halfWidth+ix]=true;
					debugNumBorder++;
				}
			}
			int debugNumOver=0;
			for (int i=0;i<len;i++){
				if (overexposed[i]) debugNumOver++;
			}
			if (debugLevel>debugTiming) printTiming("Number of overexposed/border pixels="+debugNumOver+" border pixels="+debugNumBorder+" for "+title);

			for (int numPointer=0; numPointer <whichOn.length;numPointer++) {
				for (int i=0;i<len;i++) if (overexposed[i]) diffOnOff[numPointer][i]=0.0; // too close to overexposed in no-lasers image
			}
			for (int numPointer=0; numPointer <whichOn.length;numPointer++) for (int nImg=startIndex;nImg<whichOn[numPointer].length; nImg++){
				if (whichOn[numPointer][nImg]){
					//	    					for (int i=0;i<len;i++) if ((thresholds[i]& (1<<nImg))!=0) diffOnOff[numPointer][i]=0.0; // too dim for on
					for (int i=0;i<len;i++) if ((thresholds[i]& 1)!=0) diffOnOff[numPointer][i]=0.0; // too dim for on
				}
				//	    				 else {
				//	    					for (int i=0;i<len;i++) if ((thresholds[i]& (1<<(nImg+numImages)))!=0) diffOnOff[numPointer][i]=0.0; // too bright for off
				//	    					for (int i=0;i<len;i++) if ((thresholds[i]& 2)!=0) diffOnOff[numPointer][i]=0.0; // too bright for off
				//	    				}
			}
		}
		if (debugLevel>2){
			sdfra_instance.showArrays(diffOnOff,   halfWidth, halfHeight, true, "diffOnOff-"+title, subtitles);
			if (noise!=null) sdfra_instance.showArrays(noise,   halfWidth, halfHeight, true, "noise-"+title, subtitles);
		}
		if (debugLevel>debugTiming) printTiming("diffOnOff-"+title);
		if (lpSigma>0.0) {
			for (int numPointer=0; numPointer <whichOn.length;numPointer++){
				gb.blurDouble(diffOnOff[numPointer], halfWidth, halfHeight, lpSigma, lpSigma, 0.01);
				// TODO: use different (2x?) sigma for noise? Otherwise
				if (noise!=null) gb.blurDouble(noise[numPointer], halfWidth, halfHeight, lpSigma, lpSigma, 0.01);
			}
		}
		if (debugLevel>2){
			sdfra_instance.showArrays(diffOnOff,   halfWidth, halfHeight, true, "diffSmooth-"+(threshold)+"-"+title, subtitles);
			if (noise!=null) sdfra_instance.showArrays(noise,   halfWidth, halfHeight, true, "noise_Smooth-"+(threshold)+"-"+title, subtitles);
		}
		if (debugLevel>debugTiming) printTiming("diffSmooth-"+title);
		if (noise!=null){
			for (int numPointer=0; numPointer <whichOn.length;numPointer++) for (int i=0;i<len;i++){
				noise[numPointer][i]=diffOnOff[numPointer][i]/noise[numPointer][i]; // now - s/n
			}
			if (debugLevel>2)	sdfra_instance.showArrays(noise,   halfWidth, halfHeight, true, "s2n_Smooth-"+(threshold)+"-"+title, subtitles);
			if (debugLevel>debugTiming) printTiming("s2n_Smooth-"+title);
		}


		// zero out non-local max pixels
		double [][] localMaxPixels=new double [diffOnOff.length][diffOnOff[0].length];
		for (int numPointer=0; numPointer <whichOn.length;numPointer++){
			localMaxPixels[numPointer]=diffOnOff[numPointer].clone();
			boolean [] isLocalMax=localMaximum(
					localMaxPixels[numPointer], //double [] pixels,
					halfWidth, // int width,
					(int) this.localMaxRadius/2, //); //int radius)
					debugLevel-1);
			if (debugLevel>3) sdfra_instance.showArrays(isLocalMax,   halfWidth, halfHeight, "isLocalMax-"+(numPointer)+"-"+title);

			for (int i=0;i<localMaxPixels[numPointer].length;i++) if (!isLocalMax[i]) localMaxPixels[numPointer][i]=0.0;
		}

		// remove close to overexposed and border pixels
		if (!headLaserMode) {
			for (int numPointer=0; numPointer <whichOn.length;numPointer++) {
				for (int i=0;i<len;i++) if (overexposed[i]) localMaxPixels[numPointer][i]=0.0; // too close to overexposed in no-lasers image
			}
		}

		if (debugLevel>2) sdfra_instance.showArrays(localMaxPixels,   halfWidth, halfHeight, true, "localmax-"+(threshold)+"-"+title, subtitles);
		if (debugLevel>debugTiming) printTiming("localmax-"+title);
		// filter by S/N ratio
		if (noise!=null){
			for (int numPointer=0; numPointer <whichOn.length;numPointer++){
				for (int i=0;i<localMaxPixels[numPointer].length;i++) if (noise[numPointer][i]<this.laserSignalToNoise) localMaxPixels[numPointer][i]=0.0;
			}
			if (debugLevel>2)sdfra_instance.showArrays(localMaxPixels,   halfWidth, halfHeight, true, "localmax_s2n-"+(this.laserSignalToNoise)+"-"+title, subtitles);
			if (debugLevel>debugTiming) printTiming("localmax_s2n-"+title);
		}


		double [] firstMax=  new double [whichOn.length];
		int [] firstMaxIndex=new int [whichOn.length];
		int [] numberOfMax=new int [whichOn.length];
		// final filter by a global threshold
		for (int numPointer=0; numPointer <whichOn.length;numPointer++){
			firstMax[numPointer]=0.0;
			numberOfMax[numPointer]=0;
			firstMaxIndex[numPointer]=0;
			for (int i=0;i<localMaxPixels[numPointer].length;i++) {
				if ((localMaxPixels[numPointer][i]> 0.0) && (localMaxPixels[numPointer][i]>=threshold)){ // is threshold>0.0 always
					if (localMaxPixels[numPointer][i]>firstMax[numPointer]) {
						firstMax[numPointer]=localMaxPixels[numPointer][i];
						firstMaxIndex[numPointer]=i;
					}
					numberOfMax[numPointer]++;
				} else{
					localMaxPixels[numPointer][i]=0.0;
				}
			}
		}
		// zero out all out of pattern maximums
		if (patternMask!=null){
			for (int numPointer=0; numPointer <whichOn.length;numPointer++){
				for (int i=0;i<localMaxPixels[numPointer].length;i++) if (!patternMask[i]) localMaxPixels[numPointer][i]=0.0;
			}
			if (debugLevel>2)	sdfra_instance.showArrays(localMaxPixels,   halfWidth, halfHeight, true, "in-patt-max-"+(threshold)+"-"+title, subtitles);
			if (debugLevel>debugTiming) printTiming("in-patt-max-"+title);
		}
		if (debugLevel>2)sdfra_instance.showArrays(localMaxPixels,   halfWidth, halfHeight, true, "localmax_threshold-"+(this.laserSignalToNoise)+"-"+title, subtitles);
		if (debugLevel>debugTiming) printTiming("localmax_threshold-"+(this.laserSignalToNoise)+"-"+title);
		if (debugLevel>1) {
			String dbgStr="";
			for (int numPointer=0; numPointer <numberOfMax.length;numPointer++) if (numberOfMax[numPointer]>0){
				dbgStr+=" "+numPointer+":"+IJ.d2s(firstMax[numPointer],4)+" ["+(firstMaxIndex[numPointer]%halfWidth)+":"+
						(firstMaxIndex[numPointer]/halfWidth)+"]("+numberOfMax[numPointer]+")";
			}
			if (dbgStr.length()>0) System.out.println("Maximums found: "+dbgStr);
		}
		// remove smaller maximum if it is too close to larger one, and if there are some left - find the next replacement

		while (this.fartherstOffender>0){
			int n1=-1,n2=-1;
			boolean tooClose=false;
			double co2=0.25*this.fartherstOffender*this.fartherstOffender;

			for (n1=0;n1<numberOfMax.length;n1++) if (numberOfMax[n1]>0){
				for (n2=n1+1;n2<numberOfMax.length;n2++) if (numberOfMax[n2]>0){
					double dx=(firstMaxIndex[n1]%halfWidth)-(firstMaxIndex[n2]%halfWidth);
					double dy=(firstMaxIndex[n1]/halfWidth)-(firstMaxIndex[n2]/halfWidth);
					if (dx*dx+dy*dy < co2) {
						tooClose=true;
						if (debugLevel>0) { // rare events, output details to verify program
							System.out.println ("Two detected pointers in image "+title+": "+n1+"("+numberOfMax[n1]+") and "+n2+"("+numberOfMax[n2]+") are too close,");
							System.out.println ("Distance="+IJ.d2s(2*Math.sqrt(dx*dx+dy*dy),3)+" sensor pixels is smaller than configured limit of "+this.fartherstOffender);
						}
						break;
					}

				}
				if (tooClose) break;
			}
			if (!tooClose) break;
			int numPointer= (firstMax[n1]>firstMax[n2])?n2:n1;
			if (debugLevel>0) { // rare events, output details to verify program
				int numFirstPointer=(firstMax[n1]>firstMax[n2])?n1:n2;
				System.out.println ("More intense pointer is #"+numFirstPointer+" ("+firstMax[numFirstPointer]+"), x="
						+(2*(firstMaxIndex[numFirstPointer]%halfWidth))+" y="+(2*(firstMaxIndex[numFirstPointer]/halfWidth)));
				System.out.println ("Will remove other pointer #"+numPointer+" ("+firstMax[numPointer]+"), x="
						+(2*(firstMaxIndex[numPointer]%halfWidth))+" y="+(2*(firstMaxIndex[numPointer]/halfWidth)));
			}

			localMaxPixels[numPointer][firstMaxIndex[numPointer]]=0.0;
			numberOfMax[numPointer]--;


			if (numberOfMax[numPointer]>0){ // look for the next maximum
				firstMax[numPointer]=localMaxPixels[numPointer][0];
				firstMaxIndex[numPointer]=0;
				for (int i=0;i<localMaxPixels[numPointer].length;i++) if ((localMaxPixels[numPointer][i]>0) && (localMaxPixels[numPointer][i]>firstMax[numPointer])){
					firstMax[numPointer]=localMaxPixels[numPointer][i];
					firstMaxIndex[numPointer]=i;
				}
				if (debugLevel>0) { // rare events, output details to verify program
					System.out.println ("Found replacement pointer for #"+numPointer+" ("+firstMax[numPointer]+"), x??="
							+(2*(firstMaxIndex[numPointer]%halfWidth))+" y??="+(2*(firstMaxIndex[numPointer]/halfWidth)));
					System.out.println (numberOfMax[numPointer]+ " pointer candidate"+((numberOfMax[numPointer]>1)?"s":"")+" remain"+((numberOfMax[numPointer]>1)?"s":"")+".");
				}
			}
		}

		/*
    		// remove points close (but not too close) to bright/overexposed OBSOLETE?
    		if (!headLaserMode && (this.fartherstOffender>0)) {
    			for (int numPointer=0; numPointer <whichOn.length;numPointer++){ //
//    				for (int i=0;i<diffOnOff[numPointer].length;i++) if (diffOnOff[numPointer][i]>=threshold){
    				for (int i=0;i<diffOnOff[numPointer].length;i++) if (localMaxPixels[numPointer][i]>=threshold){
    					int xc=i%halfWidth;
    					int yc=i/halfWidth;
    					int xMin=xc-this.fartherstOffender;
    					int xMax=xc+this.fartherstOffender;
    					int yMin=yc-this.fartherstOffender;
    					int yMax=yc+this.fartherstOffender;
    					int xMinC=xc-this.closestOffender;
    					int xMaxC=xc+this.closestOffender;
    					int yMinC=yc-this.closestOffender;
    					int yMaxC=yc+this.closestOffender;
    					if (xMin<0) xMin=0;
    					if (xMax>=halfWidth) xMax=halfWidth-1;
    					if (yMin<0) yMin=0;
    					if (yMax>=halfHeight) yMax=halfHeight-1;
    					for (int y=yMin;y<=yMax;y++)
    						for (int x=xMin;x<=xMax;x++)
    							if (!((x>xMinC) && (x<xMaxC) && (y>yMinC) && (y<yMaxC)) && ((thresholds[y*halfWidth+x]& 2)!=0)) {
    								diffOnOff[numPointer][i]=0.0;
    								break;
    							}
    				}
    			}
	    		if (debugLevel>2)	sdfra_instance.showArrays(diffOnOff,   halfWidth, halfHeight, true, "filteredSmooth-"+title, subtitles);
    		}
		 */
		// Finding just maximums, no centroid here
		double [][] pointers = new double[whichOn.length][];
		for (int numPointer=0; numPointer <whichOn.length;numPointer++){
			if (numberOfMax[numPointer]>0) {

				pointers[numPointer]=new double[2];
				int [] localMaxXY={firstMaxIndex[numPointer]%halfWidth,firstMaxIndex[numPointer]/halfWidth};
				double [] dMaxXY= {localMaxXY[0], localMaxXY[1]};
				if (debugLevel>2) System.out.println("Pointer "+numPointer+" max X="+dMaxXY[0]+" max Y="+dMaxXY[1]);
				if (this.quadraticScaleSigma>0.0){
					double quadSigma=this.quadraticScaleSigma*((lpSigma>0.0)?lpSigma:1.0);
					double k=0.5/(quadSigma*quadSigma);
					int range=(int) Math.round(quadSigma*3.0);
					if (range<1) range=1;
					int minX=localMaxXY[0]-range;
					int maxX=localMaxXY[0]+range;
					int minY=localMaxXY[1]-range;
					int maxY=localMaxXY[1]+range;
					if (minX<0)minX=0;
					if (minY<0)minY=0;
					if (maxX>= halfWidth) maxX= halfWidth-1;
					if (maxY>=halfHeight) maxY=halfHeight-1;
					double [][][] data =new double [(maxX-minX+1)*(maxY-minY+1)][3][];
					int index=0;
					for (int y=minY;y<=maxY;y++) for (int x=minX;x<=maxX;x++){
						data[index][0] = new double [2];
						data[index][0][0] = x-localMaxXY[0];
						data[index][0][1] = y-localMaxXY[1];
						data[index][1]=     new double[1];
						data[index][1][0]=  diffOnOff[numPointer][x+halfWidth*y];
						data[index][2]=     new double[1];
						data[index][2][0]=  Math.exp(-k*(data[index][0][0]*data[index][0][0] + data[index][0][1]*data[index][0][1]));
						index++;
					}
					double [] corrXY=(new PolynomialApproximation()).quadraticMax2d (data);
					if (corrXY!=null) {
						dMaxXY[0]+=corrXY[0];
						dMaxXY[1]+=corrXY[1];
					}
					if (debugLevel>2) { // something wrong, (E-17)
						if (corrXY!=null) System.out.println("Pointer "+numPointer+" corr X="+corrXY[0]+" corr Y="+corrXY[1]);
						else System.out.println("Pointer "+numPointer+": failed to find correction by quadratic approximation of the vicinity");
					}
				}
				// bayer shift
				pointers[numPointer][0]=2*dMaxXY[0]+1.0;
				pointers[numPointer][1]=2*dMaxXY[1]+0.0; //-1.0;
			} else pointers[numPointer]=null;
			if (debugLevel>1) {
				if (pointers[numPointer]!=null) 	System.out.println("Pointer "+numPointer+": Max="+firstMax[numPointer]+"(>"+this.threshold+"?), x="+
						pointers[numPointer][0]+", y="+pointers[numPointer][1]);
				else if (debugLevel>2) System.out.println("Pointer "+numPointer+" is not detected");
			}
		}
		if (debugLevel>debugTiming) printTiming("getPointerXY()-done for "+title);
		return pointers;
	}
}
