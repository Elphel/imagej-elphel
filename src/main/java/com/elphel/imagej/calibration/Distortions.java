package com.elphel.imagej.calibration;
/*
 **
 ** Distortions.java - Calculate lens distortion parameters from the pattern image
 **
 ** Copyright (C) 2011-2020 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  Distortions.java is free software: you can redistribute it and/or modify
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

import java.awt.Point;
import java.awt.Rectangle;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Arrays;
//import java.util.Arrays;
//import java.io.StringWriter;
import java.util.List;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import javax.swing.SwingUtilities;

import com.elphel.imagej.calibration.DistortionCalibrationData.GridImageParameters;
import com.elphel.imagej.calibration.hardware.CamerasInterface;
import com.elphel.imagej.cameras.EyesisCameraParameters;
import com.elphel.imagej.cameras.EyesisSubCameraParameters;
import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.PolynomialApproximation;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.common.WindowTools;
import com.elphel.imagej.jp4.JP46_Reader_camera;

import Jama.LUDecomposition;
import Jama.Matrix;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
//import ij.process.*;
import ij.gui.GenericDialog;
import ij.gui.PointRoi;
import ij.io.FileSaver;
import ij.io.Opener;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.text.TextWindow;
//import src.java.org.apache.commons.configuration.*;
// to work both in Eclipse and ImageJ:
// 1 - put commons-configuration-1.7.jar under ImageJ plugins directory (I used ImageJ-Elphel)
// 2 - in Eclipse project properties -> Build Path -> Libraries -> Add External jar
public class Distortions {
//	final public double hintedMaxRelativeRadius=1.2; // make adjustable?
	final public double hintedMaxRelativeRadiusToDiagonal= 1.3; // 0.96; // make adjustable?
//    int numInputs=27; // with A8...// 24;   // parameters in subcamera+...
//    int numOutputs=16; // with A8...//13;  // parameters in a single camera
	public PatternParameters patternParameters;
	public LensDistortionParameters lensDistortionParameters;
	public RefineParameters refineParameters= new RefineParameters(); //create with default values
	public FittingStrategy fittingStrategy=null;
    public double [][][][] gridOnSensor =null; // [v][u][px,py][0-value, 1..14 - derivative]
    public double [][] interParameterDerivatives=null; //new double[this.numInputs][]; //partial derivative matrix from subcamera-camera-goniometer to single camera (12x21)
    public double []   currentVector; // current variable parameter vector
    public double []   Y=null; // array of "y" - for each grid image, each defined grid node - 2 elements
    public int    []   imageStartIndex=null; // elements containing index of the start point of the selected image, first element 0, last - total number of points.
    public double []   weightFunction=null; //  array of weights for pixels (to fade values near borders), corresponding to Y array
	public double [][] dTA_dUV = null; // null or double [2][2] to return averaged {{dU/dAz,dU/dTl}{dV/dAz,dV/dTl}}.inverse

    public double      sumWeights;
    public double [][] targetXYZ=null; // array of target {x,y,z} matching each image each grid point
    public double [][] jacobian=null; // partial derivatives of fX (above) by parameters to be adjusted (rows)
    public double []   nextVector; // next variable parameter vector
    public double []   currentfX=null; // array of "f(x)" - simulated data for all images, combining pixel-X and pixel-Y (odd/even)
    public double []   nextfX=null; // array of "f(x)" - simulated data for all images, combining pixel-X and pixel-Y (odd/even)

    public double      currentRMS=-1.0; // calculated RMS for the currentVector->currentfX
    public double      nextRMS=-1.0; // calculated RMS for the nextVector->nextfX
    public double      firstRMS=-1.0; // RMS before current series of LMA started

    public double      currentRMSPure=-1.0; // calculated RMS for the currentVector->currentfX
    public double      nextRMSPure=-1.0; // calculated RMS for the nextVector->nextfX
    public double      firstRMSPure=-1.0; // RMS before current series of LMA started

    public double lambdaStepUp=   8.0; // multiply lambda by this if result is worse
    public double lambdaStepDown= 0.5; // multiply lambda by this if result is better
    public double thresholdFinish=0.001; // (copied from series) stop iterations if 2 last steps had less improvement (but not worsening )
    public int    numIterations=  100; // maximal number of iterations
    public double maxLambda=      100.0;  // max lambda to fail

    public double lambda=0.001;        // copied from series
    public double [] lastImprovements= {-1.0,-1.0}; // {last improvement, previous improvement}. If both >0 and < thresholdFinish - done
    public int    iterationStepNumber=0;
    public boolean stopEachStep=  true;  // open dialog after each fitting step
    public boolean stopEachSeries=true;  // open dialog when each fitting series finished
    public boolean stopOnFailure= true;  // open dialog when fitting series failed
    public boolean showParams=   false;   // show modified parameters
    public boolean showThisImages=false; // show debug images for the current ("this" state,before correction) state of parameters
    public boolean showNextImages=false; // show debug images for the current (after correction) state of parameters
    public boolean askFilter=     false; // show debug images for the current (after correction) state of parameters

 //   public boolean showGridCorr=  true;  // show grid correction
 //   public boolean showIndividual=true;  // show individual image residuals
 //   public double  corrScale=     1.0;   // scale grid correction before applying

    public int     seriesNumber=0; // just for the dialog
    public boolean saveSeries=false;   // just for the dialog
    public double [][][] pixelCorrection=null; // for each sensor: corr-X, corr-Y, mask, flat-field-Red, flat-field-Green, flat-field-Blue
    public String []  pathNames=null;

    // Will have to chage for different resolution
//    public int [][]   pixelCorrectionWHD= null; // For each sensor -width, height, decimation
//    public int        defaultPixelCorrectionDecimation=   1;
//   public int        defaultPixelCorrectionWidth=     2592;
//    public int        defaultPixelCorrectionHeight=    1936;

//    @Deprecated
//    public int        pixelCorrectionDecimation=   1;
//    @Deprecated
//    public int        pixelCorrectionWidth=     2592;
//    @Deprecated
//    public int        pixelCorrectionHeight=    1936;




    public double     RMSscale=Math.sqrt(2.0); // errors for x and y are calculated separately, so actual error is larger

    public boolean  showIndex=true;
    public boolean  showRMS=true;
    public boolean  showPoints=true;
    public boolean  showLensLocation=true;
    public boolean  showEyesisParameters=true;
    public boolean  showIntrinsicParameters=true;
    public boolean  showExtrinsicParameters=true;
    public int      extraDecimals=0;

    public boolean   threadedLMA=true; // use threaded/partial method to solve LMA
    public LMAArrays lMAArrays=null;
    public LMAArrays  savedLMAArrays=null;
    public long startTime=0;
    public int debugLevel=2;
    public boolean updateStatus=true;
    public int threadsMax=100;
    public AtomicInteger stopRequested=null; // 1 - stop now, 2 - when convenient

    public String [] status ={"",""};
    public double [] lastUsedManualGridHint_UV = {0.5, 0.5};

    public int getSensorWidth(int subCam) { return fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getSensorWidth(subCam);} // for the future? different sensors
    public int getSensorHeight(int subCam) { return fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getSensorHeight(subCam);}// for the future? different sensors
    public int getDecimateMasks(int subCam) { return fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getDecimateMasks(subCam);}// for the future? different sensors

    public int getSensorWidth() { return fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getSensorWidth();} // for the future? different sensors
    public int getSensorHeight() { return fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getSensorHeight();}// for the future? different sensors
    public int getDecimateMasks() { return fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getDecimateMasks();}// for the future? different sensors

    public int getSensorCorrWidth(int subCam) { return(getSensorWidth(subCam)-1)/getDecimateMasks(subCam)+1;}
    public int getSensorCorrWidth() { return(getSensorWidth()-1)/getDecimateMasks()+1;}

    public void setSensorWidth(int subCam, int v)  {
    	fittingStrategy.distortionCalibrationData.eyesisCameraParameters.setSensorWidth(subCam, v);
    }
    public void setSensorHeight(int subCam, int v) {fittingStrategy.distortionCalibrationData.eyesisCameraParameters.setSensorHeight(subCam, v);}
    public void setDecimateMasks(int subCam, int v){fittingStrategy.distortionCalibrationData.eyesisCameraParameters.setDecimateMasks(subCam, v);}
    
    public double [][] getDtaDuv(){
    	return dTA_dUV;
    }


    public class LMAArrays {
        public double [][] jTByJ=  null; // jacobian multiplied by Jacobian transposed
        public double []   jTByDiff=null; // jacobian multiplied difference vector
        @Override
		public LMAArrays clone() {
        	LMAArrays lma=new LMAArrays();
        	lma.jTByJ = this.jTByJ.clone();
        	for (int i=0;i<this.jTByJ.length;i++) lma.jTByJ[i]=this.jTByJ[i].clone();
        	lma.jTByDiff=this.jTByDiff.clone();
        	return lma;
        }
    }
    public Distortions (){}
	public Distortions (
			LensDistortionParameters lensDistortionParameters,
			PatternParameters patternParameters,
			RefineParameters refineParameters,
			AtomicInteger stopRequested
	){
//		this.patternParameters=patternParameters.clone();  // why clone here?
//		this.lensDistortionParameters=lensDistortionParameters.clone();
		this.patternParameters=patternParameters;  // why clone here?
		this.lensDistortionParameters=lensDistortionParameters;
		this.refineParameters=refineParameters;
		this.stopRequested=stopRequested;
		if (this.lensDistortionParameters!=null) {
			interParameterDerivatives=new double[this.lensDistortionParameters.getNumInputs()][];
		}

	}
//	public int getNumInputs(){return numInputs;}
//	public int getNumOutputs(){return numOutputs;}
/**
 * Prerequisites:
 * this.patternParameters, this.fittingStrategy are already initialized
 *
 */
	/*
	private void initImageSetAndGrids(){  // never used??
// Calculate patter x,y,z==0 and alpha (1.0 - inside, 0.0 - outside) for the grid
// TODO: and save/restore to file to account for non-perfect grid
		patternParameters.calculateGridGeometry();
//  Read all grid data files (4-slice TIFF images) and create  pixelsXY and  pixelsUV arrays
		fittingStrategy.distortionCalibrationData.readAllGrids(patternParameters);
		if (this.debugLevel>3) {
			for (int n=0;n<fittingStrategy.distortionCalibrationData.pixelsXY.length;n++) {
				for (int i=0;i<fittingStrategy.distortionCalibrationData.pixelsXY[n].length;i++){
					System.out.println(n+":"+i+"  "+
							fittingStrategy.distortionCalibrationData.pixelsUV[n][i][0]+"/"+
							fittingStrategy.distortionCalibrationData.pixelsUV[n][i][1]+"  "+
							IJ.d2s(fittingStrategy.distortionCalibrationData.pixelsXY[n][i][0], 2)+"/"+
							IJ.d2s(fittingStrategy.distortionCalibrationData.pixelsXY[n][i][1], 2)
					);
				}
			}
		}
	}
	*/

	public DistortionCalibrationData getDistortionCalibrationData() {
		return (fittingStrategy == null)?null:fittingStrategy.distortionCalibrationData;
	}

	public void resetGridImageMasks(){
		int numImg=fittingStrategy.distortionCalibrationData.getNumImages();
		System.out.println("resetGridImageMasks()");
		for (int imgNum=0;imgNum<numImg;imgNum++){
			fittingStrategy.distortionCalibrationData.gIP[imgNum].resetMask();
		}
	}
	// TODO - make station-dependent? Pass sensor mask and combine it?
	public void calculateGridImageMasks(
			final boolean proportional,
			final double gridMarginScale,
			final double minContrast,
			final double shrinkBlurSigma,
			final double shrinkBlurLevel,
			final int threadsMax,
			final boolean updateStatus
			){
		final int numImg=fittingStrategy.distortionCalibrationData.getNumImages();
		final  DistortionCalibrationData.GridImageParameters [] distortionCalibrationData=this.fittingStrategy.distortionCalibrationData.gIP;
		if (updateStatus) IJ.showStatus("Calculating grid image masks...");
		System.out.print("Calculating grid image masks...");
		System.out.print(" minContrast="+minContrast+" shrinkBlurSigma="+shrinkBlurSigma+" shrinkBlurLevel="+shrinkBlurLevel);

   		final AtomicInteger imageNumberAtomic = new AtomicInteger(0);
   		final AtomicInteger imageFinishedAtomic = new AtomicInteger(0);
   		final Thread[] threads = newThreadArray(threadsMax);
   		for (int ithread = 0; ithread < threads.length; ithread++) {
   			threads[ithread] = new Thread() {
   				@Override
				public void run() {
   					for (int imgNum=imageNumberAtomic.getAndIncrement(); imgNum<numImg;imgNum=imageNumberAtomic.getAndIncrement()){
//   						if (imgNum == 488) {
//   							System.out.println("calculateGridImageMasks(), imgNum="+imgNum);
//   						}
   						distortionCalibrationData[imgNum].calculateMask(
   								proportional,
   								gridMarginScale,
   			        			minContrast,
   			        			shrinkBlurSigma,
   			        			shrinkBlurLevel);
							final int numFinished=imageFinishedAtomic.getAndIncrement();
   							SwingUtilities.invokeLater(new Runnable() {
   								@Override
								public void run() {
   									if (updateStatus) IJ.showProgress(numFinished,numImg);
   								}
   							});

   					} // for (int numImage=imageNumberAtomic.getAndIncrement(); ...
   				} // public void run() {
   			};
   		}
   		startAndJoin(threads);
		if (updateStatus) IJ.showProgress(0);
		if (updateStatus) IJ.showStatus("Calculating grid image masks... DONE");
		System.out.println("  Done");

	}


/**
 * once per fitting strategy series:
 *   1) repeat for each image/point patternParameters.getXYZM(int u, int v) and create
 *      this.targetXYZ;
 *
 *   2)fittingStrategy.buildParameterMap (int numSeries)
 *   Creates map from the parameter vector index to the {grid image number, parameter number}
 *   When the parameter is shared by several images, the map points to the one which value will be used
 *   (they might be different). Timestamp of the masterImages[] is used to determine which image to use.
 *   Simultaneously creates this.reverseParameterMap that maps each of the image/parameter to the parameter vector
 *   Needs to be run for each new strategy series
 *
 * 	 3)this.currentVector=fittingStrategy.getSeriesVector(); // and save it in the class instance
 *   Calculate vector of the parameters used in LMA algorithm, extracted from the
 *   individual data, using parameter map (calculated once after changing series)
 *
 *    public double []   currentVector; // current variable parameter vector
 *
 */
	final public int filterMulti=            1;
	final public int filterContrast=         2;
	final public int filterSensor=           4;
	final public int filterTargetMask=       8;
	final public int filterTargetAlpha=     16;
	final public int filterTargetErrors=    32;
	final public int filterMaskBadNodes=    64;
	final public int filterDiameter=       128; // use measured grid "diameter" to change image weight
	final public int filterChannelWeights= 256; // different weights for channels (higher weight for bottom sensors)
	final public int filterYtoX=           512; // different weights for channels (higher weight for bottom sensors)

	final public int filterForAll=             filterMulti+filterContrast+filterSensor+filterTargetMask+filterTargetAlpha+filterTargetErrors+filterMaskBadNodes+
	filterDiameter+filterChannelWeights+filterYtoX;
	final public int filterForSensor=          filterMulti+filterContrast             +filterTargetMask+filterTargetAlpha+filterTargetErrors+filterMaskBadNodes+
	filterDiameter+filterChannelWeights+filterYtoX;
	final public int filterForTargetGeometry=  filterMulti+filterContrast+filterSensor+filterMaskBadNodes+filterDiameter+filterChannelWeights+filterYtoX;
	final public int filterForTargetFlatField= filterMulti+filterContrast+filterSensor+filterMaskBadNodes+filterDiameter+filterChannelWeights+filterYtoX;

	public int selectFilter(int dfltFilter){
		GenericDialog gd = new GenericDialog("Select series to process");
		int filter=    dfltFilter;
		gd.addCheckbox("filterMulti",         (filterForAll & filterMulti)!=0);
		gd.addCheckbox("filterContrast",      (filterForAll & filterContrast)!=0);
		gd.addCheckbox("filterSensor",        (filterForAll & filterSensor)!=0);
		gd.addCheckbox("filterTargetMask",    (filterForAll & filterTargetMask)!=0);
		gd.addCheckbox("filterTargetAlpha",   (filterForAll & filterTargetAlpha)!=0);
		gd.addCheckbox("filterTargetErrors",  (filterForAll & filterTargetErrors)!=0);
		gd.addCheckbox("filterMaskBadNodes",  (filterForAll & filterMaskBadNodes)!=0);
		gd.addCheckbox("filterDiameter",      (filterForAll & filterDiameter)!=0);
		gd.addCheckbox("filterChannelWeights",(filterForAll & filterChannelWeights)!=0);
		gd.addCheckbox("filterYtoX",          (filterForAll & filterYtoX)!=0);



		gd.showDialog();
		if (gd.wasCanceled()) return filter;
		filter=0;
		if (gd.getNextBoolean()) filter |= filterMulti;
		if (gd.getNextBoolean()) filter |= filterContrast;
		if (gd.getNextBoolean()) filter |= filterSensor;
		if (gd.getNextBoolean()) filter |= filterTargetMask;
		if (gd.getNextBoolean()) filter |= filterTargetAlpha;
		if (gd.getNextBoolean()) filter |= filterTargetErrors;
		if (gd.getNextBoolean()) filter |= filterMaskBadNodes;
		if (gd.getNextBoolean()) filter |= filterDiameter;
		if (gd.getNextBoolean()) filter |= filterChannelWeights;
		if (gd.getNextBoolean()) filter |= filterYtoX;
		if (this.debugLevel>1) System.out.println("Using filter bitmap: "+filter);
		return filter;
    }

	public void initFittingSeries(
			boolean justSelection, // use series to get selection only
			int filter,
			int numSeries) {
		if (initFittingSeries(
			justSelection, // use series to get selection only
			filter,
			numSeries,
			1)){
			initFittingSeries(
					justSelection, // use series to get selection only
					filter,
					numSeries,
					2);
		}

	}
	//returns true if some images were disabled and re-calculation is needed
	public boolean initFittingSeries(
			boolean justSelection, // use series to get selection only
			int filter,
			int numSeries,
			int pass) {
		boolean skipMinVal=this.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.minimalValidNodes<0;
		if ((pass>1) && skipMinVal){ System.out.println("initFittingSeries("+justSelection+","+filter+","+numSeries+"), skipMinVal="+skipMinVal); return false;} // debug - skipping new functionality
		System.out.println("initFittingSeries("+justSelection+","+filter+","+numSeries+"), pass="+pass);
		//TODO: ********* Implement comments above ************
		  // calculate total number of x/y pairs in the selected images
		if (numSeries <0 ) justSelection=true;
		if ((pass==1) && (numSeries>=0)) fittingStrategy.invalidateSelectedImages(numSeries); // next selectedImages() will select all, including empty
		if (!justSelection) {
			fittingStrategy.buildParameterMap (numSeries); // also sets currentSeriesNumber
		} else{
			fittingStrategy.currentSeriesNumber=numSeries;
		}
		int numXYPairs=0;
		int numImg=fittingStrategy.distortionCalibrationData.getNumImages();
		if (this.debugLevel>3)	System.out.println("initFittingSeries("+numSeries+"), numImg="+numImg);
		if ((pass==1) && (numSeries>=0) && !skipMinVal) fittingStrategy.initSelectedValidImages(numSeries); // copy from selected images

		boolean [] selectedImages=fittingStrategy.selectedImages(numSeries); // -1 OK, will select all
		if (this.debugLevel>3)	System.out.println("initFittingSeries("+numSeries+"), selectedImages.length="+selectedImages.length);
		int [] dbg_indices = new int[numImg];
		for (int imgNum=0;imgNum<numImg;imgNum++) if (selectedImages[imgNum]) {
			dbg_indices[imgNum] = numXYPairs;
			numXYPairs+=fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV.length;
		}
		this.targetXYZ=new double[numXYPairs][3];
		this.Y= new double[numXYPairs*2];
		this.weightFunction=new double[numXYPairs*2];
		this.sumWeights=0.0;
		this.imageStartIndex=new int [numImg+1];
		// added here, was using pixelCorrectionDecimation==1
///		this.pixelCorrectionDecimation=fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getDecimateMasks();
///		this.pixelCorrectionWidth=   fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getSensorWidth();
///		this.pixelCorrectionHeight=  fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getSensorHeight();
///		int sensorCorrWidth= (this.pixelCorrectionWidth-1)/this.pixelCorrectionDecimation+1;

		double [] multiWeight=new double [numImg];
		for (int imgNum=0;imgNum<numImg;imgNum++) multiWeight[imgNum]=0.0;
        double minimalGridContrast=fittingStrategy.distortionCalibrationData.eyesisCameraParameters.minimalGridContrast;
        double shrinkBlurSigma=fittingStrategy.distortionCalibrationData.eyesisCameraParameters.shrinkBlurSigma;
        double shrinkBlurLevel=fittingStrategy.distortionCalibrationData.eyesisCameraParameters.shrinkBlurLevel;
        double gridMarginScale = fittingStrategy.distortionCalibrationData.eyesisCameraParameters.gridMarginScale; // apply -scaled maximal to grid margins (_extra) for masks
        boolean proportional = false;
        calculateGridImageMasks(
        		proportional,
        		gridMarginScale,
        		minimalGridContrast, // final double minContrast,
        		shrinkBlurSigma, //final double shrinkBlurSigma,
        		shrinkBlurLevel, //final double shrinkBlurLevel,
    			100, //final int threadsMax,
    			true //final boolean updateStatus
    			);
//        this.imageSetWeight=new double[this.fittingStrategy.distortionCalibrationData.gIS.length];

        if ((filter & this.filterChannelWeights)!=0) calculateChannelsWeights(
        		this.fittingStrategy.currentSeriesNumber,
        		fittingStrategy.distortionCalibrationData.eyesisCameraParameters.balanceChannelWeightsMode);

        for (int imgSet=0;imgSet<this.fittingStrategy.distortionCalibrationData.gIS.length;imgSet++){
        	this.fittingStrategy.distortionCalibrationData.gIS[imgSet].setWeight=0.0;
        	int numUsed=0;
        	int stationNumber=0;
        	int numInSet=((this.fittingStrategy.distortionCalibrationData.gIS[imgSet].imageSet!=null)?
        			this.fittingStrategy.distortionCalibrationData.gIS[imgSet].imageSet.length:0);
        	for (int i=0;i<numInSet;i++){
        		if (this.fittingStrategy.distortionCalibrationData.gIS[imgSet].imageSet[i]!=null) {
        			stationNumber=this.fittingStrategy.distortionCalibrationData.gIS[imgSet].imageSet[i].getStationNumber(); // should be the same for all images
        			int imgNum=this.fittingStrategy.distortionCalibrationData.gIS[imgSet].imageSet[i].imgNumber;
        			if ((imgNum>=0) && selectedImages[imgNum]) numUsed++; // counting only selected in this fitting series, not all enabled !
        		}
        	}
        	if (numUsed>0) {
        		double d;
        		switch (fittingStrategy.distortionCalibrationData.eyesisCameraParameters.weightMultiImageMode){
        		case 0: d=1.0; break;
        		case 1: d=Math.pow(numUsed,fittingStrategy.distortionCalibrationData.eyesisCameraParameters.weightMultiExponent);
        		break;
        		case 2: d=(numUsed>1)?(Math.pow(numUsed,fittingStrategy.distortionCalibrationData.eyesisCameraParameters.weightMultiExponent)):0.001; break; // virtually eliminate single-image sets, but prevent errors
        		case 3: d=numUsed*numUsed; break;
        		default: d=1.0;
        		}
        		d*=fittingStrategy.distortionCalibrationData.eyesisCameraParameters.stationWeight[stationNumber];
//        		set weight will be calculated as sum of all points weights
//        		this.fittingStrategy.distortionCalibrationData.gIS[imgSet].setWeight=d;
        		for (int i=0;i<numInSet;i++){
        			if (this.fittingStrategy.distortionCalibrationData.gIS[imgSet].imageSet[i]!=null) {
        				int imgNum=this.fittingStrategy.distortionCalibrationData.gIS[imgSet].imageSet[i].imgNumber;
        				if ((imgNum>=0) && selectedImages[imgNum]) multiWeight[imgNum]= d;
        			}
        		}
        	}
        }
        int patternMaskIndex=3;
        int patternAlphaIndex=7;
        int patternErrorMaskIndex=8;
		int index=0;
		double weightScaleX=1.0,weightScaleY=1.0;
		if (((filter & this.filterYtoX)!=0) && (this.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.weightYtoX!=1.0)) {
			weightScaleX/=Math.sqrt(this.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.weightYtoX);
			weightScaleY*=Math.sqrt(this.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.weightYtoX);
		}
		double weightSumXY=weightScaleX+weightScaleY;
		for (int imgNum=0;imgNum<numImg;imgNum++){
			this.imageStartIndex[imgNum]=index;
			if (selectedImages[imgNum]) {
				int chnNum=this.fittingStrategy.distortionCalibrationData.gIP[imgNum].channel; // number of sub-camera
				int station=this.fittingStrategy.distortionCalibrationData.gIP[imgNum].getStationNumber(); // number of sub-camera
				int setNumber=this.fittingStrategy.distortionCalibrationData.gIP[imgNum].getSetNumber();
				double [] gridWeight=fittingStrategy.distortionCalibrationData.gIP[imgNum].getGridWeight();
				double gridImageWeight=1.0;
				if (((filter & this.filterDiameter)!=0) && (fittingStrategy.distortionCalibrationData.eyesisCameraParameters.weightDiameterExponent>0.0)) {
					gridImageWeight*=Math.pow(setImageDiameter(imgNum),fittingStrategy.distortionCalibrationData.eyesisCameraParameters.weightDiameterExponent);
				}
				if ((filter & this.filterChannelWeights)!=0) {
					gridImageWeight*=this.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[station][chnNum].getChannelWeightCurrent();
				}
				for (int pointNumber=0;pointNumber<fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV.length;pointNumber++){

					double [] XYZMP=patternParameters.getXYZMPE(
							fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV[pointNumber][0],
							fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV[pointNumber][1],
							station,
							chnNum,
							false);
//		 * @return null if out of grid, otherwise X,Y,Z,mask (binary),R (~0.5..1.2),G,B,alpha (0.0..1.0)
/*					double [] XYZM=patternParameters.getXYZM( // will throw if outside or masked out
							fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV[pointNumber][0],
							fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV[pointNumber][1]);*/
					if ((targetXYZ[index]==null) || (XYZMP==null)) {
						System.out.println("Null problem in imgNum="+imgNum+", point "+pointNumber);
						continue;
					}
					this.targetXYZ[index][0]=XYZMP[0];
					this.targetXYZ[index][1]=XYZMP[1];
					this.targetXYZ[index][2]=XYZMP[2];
					double weight=1.0;
					if ((filter & this.filterSensor)!=0) {
						weight*=fittingStrategy.distortionCalibrationData.getMask( // returns 1.0 if sensor mask is not available
								chnNum,
								fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsXY[pointNumber][0],
								fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsXY[pointNumber][1]);
					}
//					Individual image mask is needed as some parts can be obscured by moving parts - not present on  all images.
//					grid "contrast" may be far from 1.0 but probably should work OK
///					double gridContrast= fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsXY[pointNumber][2]-minimalGridContrast;//minimalGridContrast\
					double dbg;
//					if (weight > 0) {
					if ((weight > 0) && (imgNum == 244)) {
						dbg = weight;
					}
					if ((filter & this.filterContrast)!=0) {
						double gridContrast= gridWeight[pointNumber];
						weight*=gridContrast;
						if (Double.isNaN(gridContrast) && (this.debugLevel>1)) System.out.println("gridContrast=NaN, imgNum="+imgNum);
					}
//					if (weight > 0) {
					if ((weight > 0) && (imgNum == 244)) {
						dbg = weight;
					}

					if ((filter & this.filterTargetMask)!=0) {
						weight*=XYZMP[patternMaskIndex];//DONE: Use grid mask also (fade out outer grid nodes?)
						if (Double.isNaN(XYZMP[patternMaskIndex]) && (this.debugLevel>1)) System.out.println("XYZMP["+patternMaskIndex+"]=NaN, imgNum="+imgNum);
					}
//					if (weight > 0) {
					if ((weight > 0) && (imgNum == 244)) {
						dbg = weight;
					}
					if ((filter & this.filterTargetAlpha)!=0) {
						weight*=XYZMP[patternAlphaIndex];//DONE: Use grid mask also (fade out outer grid nodes?)
						if (Double.isNaN(XYZMP[patternAlphaIndex]) && (this.debugLevel>1)) System.out.println("XYZMP["+patternAlphaIndex+"]=NaN, imgNum="+imgNum);
					}
//					if (weight > 0) {
					if ((weight > 0) && (imgNum == 244)) {
						dbg = weight;
					}
					if ((filter & this.filterTargetErrors)!=0) {
						weight*=XYZMP[patternErrorMaskIndex];//DONE: Use grid mask also (fade out outer grid nodes?)
						if (Double.isNaN(XYZMP[patternErrorMaskIndex]) && (this.debugLevel>1)) System.out.println("XYZMP["+patternErrorMaskIndex+"]=NaN, imgNum="+imgNum);
					}
//					if (weight > 0) {
					if ((weight > 0) && (imgNum == 244)) {
						dbg = weight;
					}
					if ((filter & this.filterMulti)!=0) {
						weight*=multiWeight[imgNum];
						if (Double.isNaN(multiWeight[imgNum]) && (this.debugLevel>1)) System.out.println("multiWeight["+imgNum+"]=NaN, imgNum="+imgNum);
					}
//					if (weight > 0) {
					if ((weight > 0) && (imgNum == 244)) {
						dbg = weight;
					}
					if ((filter & this.filterMaskBadNodes)!=0) {
						if (fittingStrategy.distortionCalibrationData.gIP[imgNum].isNodeBad(pointNumber)) weight=0.0;
					}
//					if (weight > 0) {
					if ((weight > 0) && (imgNum == 244)) {
						dbg = weight; // got here
					}

					//fittingStrategy.distortionCalibrationData.eyesisCameraParameters.weightMultiExponent)
					if (((filter & this.filterDiameter)!=0) && (fittingStrategy.distortionCalibrationData.eyesisCameraParameters.weightDiameterExponent>0.0)) {
						weight*=gridImageWeight;
						if (Double.isNaN(gridImageWeight) && (this.debugLevel>1)) System.out.println("gridImageWeight=NaN, imgNum="+imgNum);
					}
//					if (weight > 0) {
					if ((weight > 0) && (imgNum == 244)) {
						dbg = weight;
					}
					if (Double.isNaN(weight)) {
						weight=0.0; // find who makes it NaN
						if (Double.isNaN(multiWeight[imgNum])) System.out.println("weight is null, imgNum="+imgNum);
					}


					this.weightFunction[2*index]=  weight*weightScaleX;
					this.weightFunction[2*index+1]=weight*weightScaleY;
					this.sumWeights+=              weight*weightSumXY;
	        		this.fittingStrategy.distortionCalibrationData.gIS[setNumber].setWeight+=2.0*weight;  // used for variances - proportional to the set weight
					if ((this.pixelCorrection==null) || (this.pixelCorrection[chnNum] == null)){
						this.Y[2*index]=  fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsXY[pointNumber][0];
						this.Y[2*index+1]=fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsXY[pointNumber][1];
						if (Double.isNaN(this.Y[2*index]) || Double.isNaN(this.Y[2*index+1])) {
							System.out.println("Bug 1 in initFittingSeries(): NaN! distortionCalibrationData.gIP["+imgNum+"].pixelsXY["+pointNumber+"][0]="+
						          fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsXY[pointNumber][0]+
						          ", distortionCalibrationData.gIP["+imgNum+"].pixelsXY["+pointNumber+"][1]="+
						          fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsXY[pointNumber][1]);
							this.Y[2*index]=      0.0;
							this.Y[2*index + 1]=  0.0;
						}
					} else {
// TODO: remove and use new code (if tested OK)
						double [] pXY={
								fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsXY[pointNumber][0],
								fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsXY[pointNumber][1]
						};
// TODO: Should it be interpolated? Correction is normally small/smooth, so it may be not important
						int indexXY=((int) Math.floor(pXY[0]/getDecimateMasks(chnNum))) +
						((int) Math.floor(pXY[1]/getDecimateMasks(chnNum)))*getSensorCorrWidth(chnNum);
						if (this.pixelCorrection[chnNum][0].length<=indexXY){
							System.out.println("initFittingSeries("+numSeries+") bug:");
							System.out.println("this.pixelCorrection["+chnNum+"][0].length="+this.pixelCorrection[chnNum][0].length);
							System.out.println("indexXY="+indexXY+" pXY[0]="+pXY[0]+", pXY[1]="+pXY[1]+" sensorCorrWidth="+getSensorCorrWidth(chnNum));

						} else {
							this.Y[2*index]=  pXY[0]-this.pixelCorrection[chnNum][0][indexXY]; //java.lang.ArrayIndexOutOfBoundsException: 3204663
							this.Y[2*index+1]=pXY[1]-this.pixelCorrection[chnNum][1][indexXY];

							if (Double.isNaN(this.Y[2*index]) || Double.isNaN(this.Y[2*index+1])) {
								System.out.println("Bug 2 in initFittingSeries(): NaN! this.pixelCorrection["+chnNum+"][0]["+indexXY+"]="+
										this.pixelCorrection[chnNum][0][indexXY]+
							          ", this.pixelCorrection["+chnNum+"][1]["+indexXY+"]="+
							          this.pixelCorrection[chnNum][1][indexXY]+
							          ", pXY[0]="+pXY[0]+", pXY[1]="+pXY[1]+
							          ", imgNum="+imgNum+", pointNumber="+pointNumber);
								this.Y[2*index]=      0.0;
								this.Y[2*index + 1]=  0.0;
							}
						}
// TODO: remove above and un-comment below	(after testing)
/*
						double [] vector=interpolateCorrectionVector(
								chnNum,
								fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsXY[pointNumber][0],
								fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsXY[pointNumber][1]);
						this.Y[2*index]=  pXY[0]-vector[0];
						this.Y[2*index+1]=pXY[1]-vector[1];
*/
					}
					index++;
				}
//				numXYPairs+=fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV.length; ??
			}
		}
		this.imageStartIndex[numImg]=index; // one after last
		if ((pass==1) && (numSeries>=0) && !skipMinVal){
    		// count non-zero weight nodes for each image, disable image if this number is less than
    		int needReCalc=0;
    		for (int imgNum=0;imgNum<numImg;imgNum++) if (selectedImages[imgNum]) {
    			index=this.imageStartIndex[imgNum];
    			int numValidNodes=0;
    			for (int pointNumber=0;pointNumber<fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV.length;pointNumber++){
    				if (2*(index+pointNumber)>=this.weightFunction.length){
    					System.out.println("BUG@535: this.weightFunction.length="+this.weightFunction.length+" index="+index+
    							" pointNumber="+pointNumber+" imgNum="+imgNum+" pixelsUV.length="+
    							fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV.length+
    							" numXYPairs="+numXYPairs);

    					continue;
    				}
    				if (this.weightFunction[2*(index+pointNumber)]>0.0) {
    					numValidNodes++; //OOB 5064
    				}
    			}
    			if (numValidNodes<this.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.minimalValidNodes){
    				this.fittingStrategy.invalidateSelectedImage(numSeries,imgNum);
    				needReCalc++;
    				if (this.debugLevel>1){
    					System.out.println("Number of valid nodes in image #"+imgNum+" is "+numValidNodes+" < "+
    							this.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.minimalValidNodes+
    							", this image will be temporarily disabled");
    				}
    			}
    		}
    		if (needReCalc>0) {
    			if (this.debugLevel>1) System.out.println("Number of temporarily disabled images="+needReCalc );
    			return true; // will need a second pass
    		} else {
    			if (this.debugLevel>1) System.out.println("No images disabled, no need for pass #2");
    		}
    	}
		// Normalize set weights
		int numSetsUsed=0+0;
		double totalSetWeight=0.0;
        for (int imgSet=0;imgSet<this.fittingStrategy.distortionCalibrationData.gIS.length;imgSet++){
        	if (this.fittingStrategy.distortionCalibrationData.gIS[imgSet].setWeight>0){
        		numSetsUsed++;
        		totalSetWeight+=this.fittingStrategy.distortionCalibrationData.gIS[imgSet].setWeight;
        	}
        }
        double setWeightScale=numSetsUsed/totalSetWeight;
        if (numSetsUsed>0){
            for (int imgSet=0;imgSet<this.fittingStrategy.distortionCalibrationData.gIS.length;imgSet++){
            	if (this.fittingStrategy.distortionCalibrationData.gIS[imgSet].setWeight>0){
            		numSetsUsed++;
            		this.fittingStrategy.distortionCalibrationData.gIS[imgSet].setWeight*=setWeightScale;
            	}
            }
        }
// last? not here! numSetsUsed counted twice (should be = 1, is 2)
//		this.imageStartIndex[numImg]=index;
		if (justSelection) {
			this.currentVector = null;
			this.lambda=0.0;

		} else {
			this.currentVector =fittingStrategy.getSeriesVector(); // here?
			// for now - use common parameters, later maybe restore /add individual
			//    	this.lambda=fittingStrategy.getLambda();
			//    	was commented out???

			this.lambda=fittingStrategy.getLambda();
		   	if ((this.fittingStrategy.varianceModes!=null)
		   			&& (this.fittingStrategy.varianceModes[numSeries]!=this.fittingStrategy.varianceModeDisabled)) fittingStrategy.buildVariancesMaps (numSeries); // return value lost
		}
//    	this.thresholdFinish=fittingStrategy.getStepDone();
    	this.iterationStepNumber=0;
    	// should be calculated after series weights are set
//    public int    []   imageStartIndex=null; // elements containing index of the start point of the selected image, first element 0, last - total number of points.
// TODO: add copying  lambdaStepUp,lambdaStepDown?
    	return false;
	}


	public void calculateChannelsWeights(
			int numSeries,
			double balanceChannelWeightsMode){
		if (balanceChannelWeightsMode==0) return; // keep current weights
		int numImg=fittingStrategy.distortionCalibrationData.getNumImages();
		int numStations=fittingStrategy.distortionCalibrationData.getNumStations();
		int numChannels=fittingStrategy.distortionCalibrationData.getNumChannels();

		if (balanceChannelWeightsMode<0) { //copy specified defaults to current values
			for (int station=0;station<numStations;station++){
				for (int chn=0;chn<numChannels;chn++){
					this.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[station][chn].setChannelWeightCurrent(
							this.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[0][chn].getChannelWeightDefault()); // from station 0
				}
			}
		} else {
			double exp=balanceChannelWeightsMode;
			double [][] sumChnWeights=new double [numStations][numChannels];
			double [] avgWeights=new double [numStations];
			int [] numNonzeroChannels=new int [numStations];
			for (int station=0;station<numStations;station++){
				avgWeights[station] =0.0;
				numNonzeroChannels[station] =0;
				for (int chn=0;chn<numChannels;chn++) sumChnWeights[station][chn] =0.0;
			}
			boolean [] selectedImages=fittingStrategy.selectedImages(numSeries); // -1 OK, will select all
			for (int imgNum=0;imgNum<numImg;imgNum++)if (selectedImages[imgNum]) {
					int chn=this.fittingStrategy.distortionCalibrationData.gIP[imgNum].channel; // number of sub-camera
					int station=this.fittingStrategy.distortionCalibrationData.gIP[imgNum].getStationNumber(); // number of sub-camera
					sumChnWeights[station][chn]+=this.fittingStrategy.distortionCalibrationData.gIP[imgNum].getNumContrastNodes(
							this.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.minimalGridContrast);
			}
			for (int station=0;station<numStations;station++){
				for (int chn=0;chn<numChannels;chn++) if (sumChnWeights[station][chn]>0) {
					avgWeights[station]+=sumChnWeights[station][chn];
					numNonzeroChannels[station]++;
				}
				if (numNonzeroChannels[station]>0) avgWeights[station]/=numNonzeroChannels[station];
			}
			for (int station=0;station<numStations;station++){
				for (int chn=0;chn<numChannels;chn++) if (sumChnWeights[station][chn]>0) {
					double weight=(sumChnWeights[station][chn]>0.0)?Math.pow(avgWeights[station]/sumChnWeights[station][chn],exp):0.0;
					this.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[station][chn].setChannelWeightCurrent(
							weight);
				}
			}

		}

	}

	public double setImageDiameter(int imgNum){
		int debugThreshold=2;
		int chnNum=this.fittingStrategy.distortionCalibrationData.gIP[imgNum].channel; // number of sub-camera
        double minimalGridContrast=this.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.minimalGridContrast;
		int station=this.fittingStrategy.distortionCalibrationData.gIP[imgNum].getStationNumber(); // number of sub-camera
		EyesisSubCameraParameters esp=this.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[station][chnNum];
        double r0pix=1000.0*esp.distortionRadius/esp.pixelSize;
        this.fittingStrategy.distortionCalibrationData.gIP[imgNum].setImageDiameter( // need to get image center px,py. Maybe r0 - use to normalize result diameter
    			esp.px0, // double xc,
    			esp.py0, // double yc,
    			r0pix,   // double r0,
    			minimalGridContrast,//  double minContrast
    			(this.debugLevel>debugThreshold)?imgNum:-1);
        return this.fittingStrategy.distortionCalibrationData.gIP[imgNum].getGridDiameter();
	}

	public void listImageSets(int mode){ // TODO: use series -1 - should work now
//		boolean [] oldSelection=this.fittingStrategy.selectAllImages(0); // enable all images in series 0
		if (this.fittingStrategy.distortionCalibrationData.gIS!=null){
			if (this.debugLevel>2){
				System.out.println("listImageSets() 1: ");
				for (int is=0;is<this.fittingStrategy.distortionCalibrationData.gIS.length;is++){
					System.out.println("listImageSets() 1: "+is+
							": tilt="+    this.fittingStrategy.distortionCalibrationData.gIS[is].goniometerTilt+
							" axial="+    this.fittingStrategy.distortionCalibrationData.gIS[is].goniometerAxial+
							" interAxis="+this.fittingStrategy.distortionCalibrationData.gIS[is].interAxisAngle+
							" estimated="+this.fittingStrategy.distortionCalibrationData.gIS[is].orientationEstimated);
				}
			}
		}

		int filter=this.filterForAll;
		if (this.askFilter) filter=selectFilter(filter);
		initFittingSeries(false,filter,-1); // first step in series
		if (this.fittingStrategy.distortionCalibrationData.gIS!=null){
			if (this.debugLevel>2){
				System.out.println("listImageSets() 2: ");
				for (int is=0;is<this.fittingStrategy.distortionCalibrationData.gIS.length;is++){
					System.out.println("listImageSets() 2: "+is+
							": tilt="+    this.fittingStrategy.distortionCalibrationData.gIS[is].goniometerTilt+
							" axial="+    this.fittingStrategy.distortionCalibrationData.gIS[is].goniometerAxial+
							" interAxis="+this.fittingStrategy.distortionCalibrationData.gIS[is].interAxisAngle+
							" estimated="+this.fittingStrategy.distortionCalibrationData.gIS[is].orientationEstimated);
				}
			}
		}
//	    initFittingSeries(true,this.filterForAll,0); // will set this.currentVector
		this.currentfX=calculateFxAndJacobian(this.currentVector, false); // is it always true here (this.jacobian==null)
		double [] errors=calcErrors(calcYminusFx(this.currentfX));
		int [] numPairs=calcNumPairs();

	    int [][] imageSets=this.fittingStrategy.distortionCalibrationData.listImages(
	    		false, // true - only enabled images
	    		null);    // do not filter eo, lwir
	    boolean hasLWIR = this.fittingStrategy.distortionCalibrationData.hasSmallSensors();

	    int [] numSetPoints=new int [imageSets.length*(hasLWIR?2:1)];
	    double [] rmsPerSet=new double[imageSets.length*(hasLWIR?2:1)];
	    int [][] numImgPoints=new int [imageSets.length][this.fittingStrategy.distortionCalibrationData.getNumSubCameras()];
	    double [][] rmsPerImg=new double[imageSets.length][this.fittingStrategy.distortionCalibrationData.getNumSubCameras()];


	    boolean [] hasNaNInSet=new boolean[imageSets.length*(hasLWIR?2:1)];
	    if (hasLWIR) {
	    	for (int setNum=0;setNum<imageSets.length;setNum++){
	    		double [] error2= {0.0,0.0};
	    		int [] numInSet= {0,0};
	    		hasNaNInSet[2*setNum]=false;
	    		hasNaNInSet[2*setNum+1]=false;
	    		for (int imgInSet=0;imgInSet<imageSets[setNum].length;imgInSet++) { // upper limit depends (now 4/20)
	    			int imgNum=imageSets[setNum][imgInSet]; // image number
	    			int chn =    fittingStrategy.distortionCalibrationData.gIP[imgNum].getChannel();
	    			int isLwir = fittingStrategy.distortionCalibrationData.isSmallSensor(imgNum)?1:0;
	    			int num=numPairs[imgNum];
	    			rmsPerImg[setNum][chn] = errors[imgNum];
	    			numImgPoints[setNum][chn] = num;
	    			if (Double.isNaN(errors[imgNum])){
	    				hasNaNInSet[2 * setNum + isLwir]=true;
	    			} else {
	    				error2[isLwir]+=errors[imgNum]*errors[imgNum]*num;
	    				numInSet[isLwir]+=num;
	    			}
	    		}
	    		numSetPoints[2 * setNum + 0]= numInSet[0];
	    		rmsPerSet   [2 * setNum + 0]= Math.sqrt(error2[0]/numInSet[0]);
	    		numSetPoints[2 * setNum + 1]= numInSet[1];
	    		rmsPerSet   [2 * setNum + 1]= Math.sqrt(error2[1]/numInSet[1]);
	    	}

	    } else {
	    	for (int setNum=0;setNum<imageSets.length;setNum++){
	    		double error2=0.0;
	    		int numInSet=0;
	    		hasNaNInSet[setNum]=false;
	    		for (int imgInSet=0;imgInSet<imageSets[setNum].length;imgInSet++) {
	    			int imgNum=imageSets[setNum][imgInSet];
	    			int chn =    fittingStrategy.distortionCalibrationData.gIP[imgNum].getChannel();
	    			int num=numPairs[imgNum];
	    			rmsPerImg[setNum][chn] = errors[imgNum];
	    			numImgPoints[setNum][chn] = num;
	    			if (Double.isNaN(errors[imgNum])){
	    				hasNaNInSet[setNum]=true;
	    			} else {
	    				error2+=errors[imgNum]*errors[imgNum]*num;
	    				numInSet+=num;
	    			}
	    		}
	    		numSetPoints[setNum]=numInSet;
	    		rmsPerSet[setNum]=Math.sqrt(error2/numInSet);
	    	}
	    }
	    this.fittingStrategy.distortionCalibrationData.listImageSet(
	    		mode,
	    		numSetPoints,
	    		rmsPerSet,
	    		hasNaNInSet,
    			numImgPoints,
    			rmsPerImg
	    		);
//		this.fittingStrategy.setImageSelection(0, oldSelection); // restore original selection in series 0
	}


	public void updateSensorMasks(){
		int alphaIndex=2;
		if (this.pixelCorrection==null){
			System.out.println("Sensor data is null, can not update sensor masks");
			return;
		}
		if (this.debugLevel>0) System.out.println("Updating sensor masks in sensor data");
		for (int i=0;(i<this.fittingStrategy.distortionCalibrationData.sensorMasks.length) && (i<this.pixelCorrection.length);i++){
			this.pixelCorrection[i][alphaIndex]=this.fittingStrategy.distortionCalibrationData.sensorMasks[i].clone();
		}

	}

	public boolean correctPatternFlatField(boolean enableShow){
		if (this.debugLevel>0) System.out.println("=== Performing pattern flat field correction");
		this.patternParameters.updateNumStations(this.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getNumStations());
		double [][] masks= nonVignettedMasks(
				this.refineParameters);

		boolean same_size = true;
		for (int nChn=1; nChn < masks.length; nChn++) same_size &= (masks[nChn].length == masks[0].length);


		if (enableShow && this.refineParameters.flatFieldShowSensorMasks) {
			if (same_size) {
			ShowDoubleFloatArrays.showArrays( //java.lang.ArrayIndexOutOfBoundsException: 313632
				masks,
				getSensorWidth(0)/ getDecimateMasks(0),
				getSensorHeight(0)/getDecimateMasks(0),
				true,
		"nonVinetting masks");
			} else {
				System.out.println ("Can not display different size masks in a stack");
			}
		}


		double [][][][] sensorGrids=calculateGridFlatField(
				this.refineParameters.flatFieldSerNumber,
				masks,
				this.refineParameters);
//				this.refineParameters.flatFieldMinimalContrast,
//				this.refineParameters.flatFieldMinimalAccumulate,
//				this.refineParameters.flatFieldUseInterpolate,
//				this.refineParameters.flatFieldMaskThresholdOcclusion, // suspect occlusion only if grid is missing in the area where sensor mask is above this threshold
//				this.refineParameters.flatFieldShrinkOcclusion,
//				this.refineParameters.flatFieldFadeOcclusion,
//				this.refineParameters.flatFieldIgnoreSensorFlatField);
		double [][][] geometry= patternParameters.getGeometry();
		if (enableShow && this.refineParameters.flatFieldShowIndividual){
			for (int station=0;station<sensorGrids.length;station++) if (sensorGrids[station]!=null){
				for (int i=0;i<sensorGrids[station].length;i++) if (sensorGrids[station][i]!=null){
					ShowDoubleFloatArrays.showArrays(
							sensorGrids[station][i],
							geometry[0].length,
							geometry.length,
							true,
							"chn"+i+":"+station+"-pattern");
				}
			}
		}


		double [][][][] patternArray= combineGridFlatField(
//				this.refineParameters,
				this.refineParameters.flatFieldReferenceStation,
				sensorGrids,
				this.refineParameters.flatFieldShrinkForMatching,
				this.refineParameters.flatFieldResetMask,
				this.refineParameters.flatFieldMaxRelDiff,
				this.refineParameters.flatFieldShrinkMask,
				this.refineParameters.flatFieldFadeBorder);
		if (enableShow && this.refineParameters.flatFieldShowResult) {
			String [] titles={"Alpha","Red","Green","Blue","Number of images used"};
			for (int station=0;station<patternArray.length;station++) if (patternArray[station]!=null){
				for (int nView=0;nView<patternArray[station].length;nView++) if (patternArray[station][nView]!=null){
					ShowDoubleFloatArrays.showArrays(
							patternArray[station][nView],
							geometry[0].length,
							geometry.length,
							true,
							"St"+station+"_V"+nView+"_Pattern_Colors "+this.refineParameters.flatFieldMaxRelDiff,
							titles);
				}
			}
		}
		if (this.refineParameters.flatFieldApplyResult) applyGridFlatField(patternArray); // {alpha, red,green,blue, number of images used}[pixel_index]
		return true;
	}

	public boolean modifyPixelCorrection(
			boolean   enableShow,
			int       threadsMax,
			boolean   updateStatus,
			int debugLevel
	){
		int filter=this.filterForSensor;
		if (this.askFilter) filter=selectFilter(filter);
    	initFittingSeries(true,filter,this.seriesNumber); // first step in series now uses pattern alpha
//    	initFittingSeries(true,this.filterForSensor,this.seriesNumber); // first step in series now uses pattern alpha
    	this.currentfX=calculateFxAndJacobian(this.currentVector, false);
    	//        	this.currentRMS= calcError(calcYminusFx(this.currentfX));
    	if (this.debugLevel>2) {
    		System.out.println("this.currentVector");
    		for (int i=0;i<this.currentVector.length;i++){
    			System.out.println(i+": "+ this.currentVector[i]);
    		}
    	}
		boolean [] selectedImages=fittingStrategy.selectedImages();
		double [][][] sensorXYRGBCorr=  allImagesCorrectionMapped(
				selectedImages,
				enableShow, //  && this.refineParameters.showPerImage,
				this.refineParameters, // .showIndividualNumber,
				threadsMax,
				updateStatus,
				debugLevel);
		String [] titles={"X-corr(pix)","Y-corr(pix)","weight","Red","Green","Blue"};
		for (int numChn=0;numChn<sensorXYRGBCorr.length;numChn++) if (sensorXYRGBCorr[numChn]!=null){
//			boolean small_sensor = fittingStrategy.distortionCalibrationData.isSmallSensor(numChn);
			boolean small_sensor = fittingStrategy.distortionCalibrationData.getSmallSensors()[numChn];
			RefineParameters rp = small_sensor ? refineParameters.refineParametersSmall : refineParameters;
			int decimate=getDecimateMasks(numChn); // Reduce for LWIR? Make form sensor width?
			int sWidth= (getSensorWidth(numChn)-1)/decimate+1;
			if (rp.showUnfilteredCorrection &&  enableShow) { //  && this.refineParameters.showUnfilteredCorrection) {
				showWithRadialTangential(
						titles,
						"chn_"+numChn+"_visible_correction",
						sensorXYRGBCorr[numChn], // [0] - dx, [1] - dy
						sWidth,
						decimate,
						fittingStrategy.distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[0][numChn].px0, // using station 0
						fittingStrategy.distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[0][numChn].py0);
			}
			smoothSensorPolar(
					titles,
					"chn_"+numChn+"_polar_correction",
					sensorXYRGBCorr[numChn], // [0] - dx, [1] - dy
					sWidth,
					decimate,
					fittingStrategy.distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[0][numChn].px0, // using station 0
					fittingStrategy.distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[0][numChn].py0,
					rp.center_fract,       // 0.5, // double center_fract, // 0.5 of half-height
					rp.transit_fract,      // 0.2, // double transit_fract, // 0.2 of half-height - transition from center ortho to outer polar
					rp.gaus_ang,           // 0.2, // 1, // 0.2, // double gaus_ang,  // in radians
					rp.gaus_rad,           // 0.05, // 0.1, // double gaus_rad   // in fractions of the full radius
					rp.max_diff_err_geom,  // 0.25, // double max_diff_err_geom,   // before second pass linearly fade R/T and RGB where high-frequency error nears thios value
					rp.max_diff_err_photo, // 0.25, // double max_diff_err_photo
					rp.showExtrapolationCorrection &&  enableShow); // boolean showDebugImages);

			if (rp.showThisCorrection &&  enableShow) { //  && this.refineParameters.showUnfilteredCorrection) {
				showWithRadialTangential(
						titles,
						"after-chn_"+numChn+"_after_polar",
						sensorXYRGBCorr[numChn], // [0] - dx, [1] - dy
						sWidth,
						decimate,
						fittingStrategy.distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[0][numChn].px0, // using station 0
						fittingStrategy.distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[0][numChn].py0);
			}
		}
		//eyesisSubCameras
		//extrapolate
		// TODO: different extrapolation for FF - not circular where shades are in effect (top/bottom)
		//    	if (!this.refineParameters.sensorExtrapolateDiff) { // add current correction BEFORE extrapolating/blurring
		addOldXYCorrectionToCurrent(
				//    				true, // if (!this.refineParameters.sensorExtrapolateDiff)
				this.refineParameters, // .correctionScale,
				sensorXYRGBCorr
				);
		/*
		allSensorsExtrapolationMapped(
				0, //final int stationNumber, // has to be selected
				sensorXYRGBCorr, //final double [][][] gridPCorr,
				this.refineParameters,
				threadsMax,
				updateStatus,
				enableShow && this.refineParameters.showExtrapolationCorrection, //final boolean showDebugImages,
				debugLevel
				);
		 */
 //   	if (this.refineParameters.smoothCorrection) {
		/*
    	boolean [] whichBlur={true,true,false,true,true,true}; // all but weight
    	IJ.showStatus("Bluring sensor corrections...");
    	for (int numChn=0;numChn<sensorXYRGBCorr.length;numChn++) if (sensorXYRGBCorr[numChn]!=null){
    		boolean small_sensor = fittingStrategy.distortionCalibrationData.isSmallSensor(numChn);
    		RefineParameters rp = small_sensor ? refineParameters.refineParametersSmall : refineParameters;
    		if (rp.smoothCorrection) {
    			int decimate=getDecimateMasks(numChn);
    			int sWidth= (getSensorWidth(numChn)-1)/decimate+1;
    			int sHeight=(getSensorHeight(numChn)-1)/decimate+1;
    			DoubleGaussianBlur gb=new DoubleGaussianBlur();
    			for (int m=0;m<whichBlur.length;m++) if (whichBlur[m]){
    				gb.blurDouble(
    						sensorXYRGBCorr[numChn][m],
    						sWidth,
    						sHeight,
    						rp.smoothSigma/decimate,
    						rp.smoothSigma/decimate,
    						0.01);
    			}
    			IJ.showProgress(numChn+1, sensorXYRGBCorr.length);
    		}
    	}
    	IJ.showProgress(1.0);
 //   	}
		 */
//    	if (enableShow && this.refineParameters.showThisCorrection ) {
		/*
    	if (enableShow) {
    		for (int numChn=0;numChn<sensorXYRGBCorr.length;numChn++) if (sensorXYRGBCorr[numChn]!=null){
    			boolean small_sensor = fittingStrategy.distortionCalibrationData.isSmallSensor(numChn);
    			RefineParameters rp = small_sensor ? refineParameters.refineParametersSmall : refineParameters;
    			if (rp.showThisCorrection) {
    				int decimate=getDecimateMasks(numChn);
    				int sWidth= (getSensorWidth(numChn)-1)/decimate+1;
    				//    		   ShowDoubleFloatArrays.showArrays(sensorXYRGBCorr[numChn], sWidth, sHeight,  true, "chn_"+numChn+"_filtered", titles);
    				showWithRadialTangential(
    						titles,
    						"chn_"+numChn+"_filtered",
    						sensorXYRGBCorr[numChn], // [0] - dx, [1] - dy
    						sWidth,
    						decimate,
    						fittingStrategy.distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[0][numChn].px0, // using station 0
    						fittingStrategy.distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[0][numChn].py0);
    			}
    		}
    	}
    	*/
		/*
//    	if (this.refineParameters.sensorExtrapolateDiff) { // add current correction AFTER extrapolationg/bluring
    		addOldXYCorrectionToCurrent(
    				false, // if (this.refineParameters.sensorExtrapolateDiff) { // add current correction AFTER extrapolationg/bluring
    				this.refineParameters, // .correctionScale,
    				sensorXYRGBCorr
    		);
//    	}
	    */
//   	if (!selectCorrectionScale()) return false;

		IJ.showStatus("Applying corrections:"+((!this.refineParameters.applyCorrection && !this.refineParameters.applyFlatField)?
				"none ":((this.refineParameters.applyCorrection?"geometry ":"")+(this.refineParameters.applyFlatField?"flat field":""))));
    	boolean result=applySensorCorrection(
    			this.refineParameters,
    			sensorXYRGBCorr, //sensorXYCorr, // modified to accept both 7(old) and 6(new) entries
    			fittingStrategy.distortionCalibrationData);
    	if (enableShow && this.refineParameters.showCumulativeCorrection) {
    		for (int numChn=0;numChn<sensorXYRGBCorr.length;numChn++) if (sensorXYRGBCorr[numChn]!=null){
				int decimate=getDecimateMasks(numChn);
				int sWidth= (getSensorWidth(numChn)-1)/decimate+1;
//    		   ShowDoubleFloatArrays.showArrays(sensorXYRGBCorr[numChn], sWidth, sHeight,  true, "Cumulative_chn_"+numChn+"_corrections", titles);
				showWithRadialTangential(
						titles,
						"Cumulative_chn_"+numChn+"_corrections",
						sensorXYRGBCorr[numChn], // [0] - dx, [1] - dy
						sWidth,
						decimate,
						fittingStrategy.distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[0][numChn].px0, // using station 0
						fittingStrategy.distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[0][numChn].py0);
    		}
    	}
    	if (result) {
			updateCameraParametersFromCalculated(false); // update camera parameters from enabled only images (may overwrite some of the above)

    	}
		IJ.showStatus("");
    	return result;
    }

	public void showWithRadialTangential(
			String [] preTitles,
			String title,
			double [][] preData, // [0] - dx, [1] - dy
			int width,
			int decimate,
			double x0,
			double y0){
		int indexDx=0;
		int indexDy=1;
		int indexDr=0;
		int indexDt=1;
		int indexDa=2;
		String [] extraTitles={"R-corr(pix)","T-corr{pix)","A-corr(pix)"};
		int newImages=extraTitles.length;
		int length=preData[0].length;
		int height=length/width;
		double [][] data= new double [preData.length+newImages] [length];
		String [] titles= new String [preTitles.length+newImages];
		for (int i=0;i<preData.length;i++){
			data[i+newImages]=preData[i];
			titles[i+newImages]=preTitles[i];
		}
		for (int i=0;i<newImages;i++){
			titles[i]=extraTitles[i];
			data[i]=new double[length];
		}
		Point2D Z=new Point2D.Double(0.0,0.0);
		for (int i=0;i<length;i++){
			Point2D R=new Point2D.Double((decimate*(i%width))-x0,(decimate*(i/width))-y0);
			double r=R.distance(Z);
			Point2D uR=new Point2D.Double(1.0,0.0);
			if (r>0) uR.setLocation(R.getX()/r,R.getY()/r);
			Point2D dXY=new Point2D.Double(preData[indexDx][i],preData[indexDy][i]);
			data[indexDr][i]= dXY.getX()*uR.getX()+dXY.getY()*uR.getY();
			data[indexDt][i]=-dXY.getX()*uR.getY()+dXY.getY()*uR.getX();
			data[indexDa][i]=dXY.distance(Z);
		}
	   ShowDoubleFloatArrays.showArrays(data, width, height,  true, title, titles);
	}





	public static double [][] DXYtoDRT(
			double [][] dataXY, // [0] - dx, [1] - dy
			int width,
			int decimate,
			double x0,
			double y0){
		int length=dataXY[0].length;
//		int height=length/width;
		double [][] dataRT = new double[2][length];
		Point2D Z=new Point2D.Double(0.0,0.0);
		for (int i=0;i<length;i++){
			Point2D R=new Point2D.Double((decimate*(i%width))-x0,(decimate*(i/width))-y0);
			double r=R.distance(Z);
			Point2D uR=new Point2D.Double(1.0,0.0);
			if (r>0) uR.setLocation(R.getX()/r,R.getY()/r);
			Point2D dXY=new Point2D.Double(dataXY[0][i],dataXY[1][i]);
			dataRT[0][i]= dXY.getX()*uR.getX()+dXY.getY()*uR.getY();
			dataRT[1][i]=-dXY.getX()*uR.getY()+dXY.getY()*uR.getX();
		}
		return dataRT;
	}

	public double [][] DRTtoDXY(
			double [][] dataRT, // [0] - dx, [1] - dy
			int width,
			int decimate,
			double x0,
			double y0){

		int length=dataRT[0].length;
//		int height=length/width;
		double [][] dataXY = new double[2][length];
		Point2D Z=new Point2D.Double(0.0,0.0);
		for (int i=0;i<length;i++){
			Point2D R=new Point2D.Double((decimate*(i%width))-x0,(decimate*(i/width))-y0);
			double r=R.distance(Z);
			Point2D uR=new Point2D.Double(1.0,0.0);
			if (r>0) uR.setLocation(R.getX()/r,R.getY()/r); // cos, sin
			Point2D dRT=new Point2D.Double(dataRT[0][i],dataRT[1][i]);
			dataXY[0][i]= dRT.getX()*uR.getX()-dRT.getY()*uR.getY();
			dataXY[1][i]= dRT.getX()*uR.getY()+dRT.getY()*uR.getX();
		}
		return dataXY;
	}




	public double [][] correctionToRadial(
			double [][] preData, // variable number of layers, [0] - weight dRT or RGB
			int width,
			int decimate,
			double x0,
			double y0,
			int inner_rad_d, // margins (minimal from all sides)
			int lines_overlap, // add linesto the bottom (corresponding to the angle)
			int indexWeight) {
		int length=preData[0].length;
		int height=length/width;
		int [] ints = radialInts(
				length,
				width,
				decimate,
				x0,
				y0,
				inner_rad_d); // margins (minimal from all sides)
//		int inner_rad_d =    ints[0];
		int polar_height0 = ints[0];
		int polar_width =   ints[1];

		int polar_height =  polar_height0 + lines_overlap;


		double [][] polar = new double[preData.length][polar_width * polar_height];
		for (int i = 0; i < polar.length; i++) if (i != indexWeight){
 			for (int j = 0; j < polar[i].length; j++) {
 				polar[i][j] = Double.NaN;
 			}
 		}
 		double k = 2* Math.PI / polar_height0;
 		for (int iangle = 0; iangle < polar_height; iangle++) {
 			double angle = k * iangle;
 			double s = Math.sin(angle);
 			double c = Math.cos(angle);
 			for (int ir = 0; ir < polar_width; ir++) {
 				double r = decimate*(ir + inner_rad_d);
 				// calculate x,y - no interpolation here
 				int ix = (int) Math.round((x0 + r * c)/decimate);
 				int iy = (int) Math.round((y0 + r * s)/decimate);
 				if ((ix >= 0) && (iy >= 0) && (ix < width) && (iy < height)) {
 					int indx = ix+ iy* width;
 					int indx_polar = ir + iangle* polar_width;
 					for (int l = 0; l < polar.length; l++) {
 						polar[l][indx_polar] = preData[l][indx];
 					}
 				}
 			}
 		}
 		return polar;
	}

	public void correctionToOrtho(
			double [][] data_ortho, // variable number of layers, will be modifier from the overlapping radial)
			double [][] data_polar, // same layers as ortho, will be used as a source
			int width,
			int decimate,
			double x0,
			double y0,
			int inner_rad_d, // margins (minimal from all sides)
			int lines_overlap, // add lines to the bottom (corresponding to the angle)
			int transit_rad) {
		int length=data_ortho[0].length;
		int height=length/width;

		int [] ints = radialInts(
				length,
				width,
				decimate,
				x0,
				y0,
				inner_rad_d); // margins (minimal from all sides)

		double [] transition = new double[transit_rad];
		for (int i = 0; i < transit_rad; i++) {
			transition[i] = 0.5*(1.0-Math.cos(Math.PI * (i +0.5)/transit_rad));
		}
		int polar_height0 = ints[0];
		int polar_width =   ints[1];
// 		double k = 2 * Math.PI / polar_height0;
 		int iangle0 = lines_overlap / 2; // polar was created with overlap (pver2PI) for angular blurring, use middle portion

 		for (int iy = 0; iy < height; iy++) {
 			double rs = iy * decimate - y0;
 			for (int ix = 0; ix < width; ix++) {
 				int indx = ix+ iy* width;
 	 			double rc = ix * decimate - x0;
 	 			double r = Math.sqrt(rc*rc + rs*rs);
 	 			double angle = Math.atan2(rs, rc);
 				if (angle < 0.0) angle += 2* Math.PI;
 				int iangle = (int) Math.round(polar_height0 * angle/ (2 * Math.PI));
 				if (iangle < iangle0) {
 					iangle += polar_height0;
 				}
 				int ir = (int) Math.round(r/decimate - inner_rad_d);
 				if ((ir >=0) && (ir < polar_width)) {
 					int indx_polar = ir + iangle* polar_width;
 					if (ir >= transit_rad){
 						for (int l = 0; l < data_ortho.length; l++) {
 							data_ortho[l][indx] = data_polar[l][indx_polar];
 						}
 					} else {
 						for (int l = 0; l < data_ortho.length; l++) {
 							data_ortho[l][indx] = transition[ir]*data_polar[l][indx_polar] + (1.0-transition[ir])*data_ortho[l][indx];
 						}
 					}
 				}
 			}
 		}
	}


	public void correctionToOrtho0XXX(
			double [][] data_ortho, // variable number of layers, will be modifier from the overlapping radial)
			double [][] data_polar, // same layers as ortho, will be used as a source
			int width,
			int decimate,
			double x0,
			double y0,
			int inner_rad_d, // margins (minimal from all sides)
			int lines_overlap){ // add lines to the bottom (corresponding to the angle)
		int length=data_ortho[0].length;
		int height=length/width;

		int [] ints = radialInts(
				length,
				width,
				decimate,
				x0,
				y0,
				inner_rad_d); // margins (minimal from all sides)

		int polar_height0 = ints[0];
		int polar_width =   ints[1];
 		double k = 2 * Math.PI / polar_height0;
 		int iangle0 = lines_overlap / 2; // polar was created with overlap (pver2PI) for angular blurring, use middle portion

 		for (int iangle = iangle0; iangle < (iangle0 + polar_height0); iangle++) {
 			double angle = k * iangle;
 			double s = Math.sin(angle);
 			double c = Math.cos(angle);
 			for (int ir = 0; ir < polar_width; ir++) {
 				double r = decimate*(ir + inner_rad_d);
 				// calculate x,y - no interpolation here
 				int ix = (int) Math.round((x0 + r * c)/decimate);
 				int iy = (int) Math.round((y0 + r * s)/decimate);
 				if ((ix >= 0) && (iy >= 0) && (ix < width) && (iy < height)) {
 					int indx = ix+ iy* width;
 					int indx_polar = ir + iangle* polar_width;
 					for (int l = 0; l < data_ortho.length; l++) {
 						data_ortho[l][indx] = data_polar[l][indx_polar];
 					}
 				}
 			}
 		}
	}

	public int [] radialInts( // helper to correctionToRadial() and correctionToOrtho()
			int length,
			int width,
			int decimate,
			double x0,
			double y0,
			int inner_rad_d) { // margins (minimal from all sides)
		int height=length/width;

		//find longest distance from the center to the corners
		Point2D C=new Point2D.Double(x0, y0);
		Point2D [] corners = {
				new Point2D.Double(-decimate,     -decimate),
				new Point2D.Double(decimate*(width+1),-decimate),
				new Point2D.Double(decimate*(width+1), decimate*(height+1)),
				new Point2D.Double(-decimate,      decimate*(height+1)),
		};
		double radius = 0.0;
		for (Point2D p :corners) {
			radius = Math.max(radius, C.distance(p));
		}
//		double x0d = x0/decimate;
//		double y0d = y0/decimate;
//		double min_rad_d_h = Math.min(x0d, width- 1-x0d);
//		double min_rad_d_v = Math.min(y0d, height-1-y0d);

		int imax_rad_d = (int) Math.ceil(radius/decimate);

//		int inner_rad_d = (int) Math.floor(Math.min(min_rad_d_h,min_rad_d_v)) -1 - margins;
		int polar_height0 = (int) Math.ceil(2*Math.PI*radius/decimate);
		int polar_width =  imax_rad_d - inner_rad_d + 1;
		int [] rslt = {polar_height0, polar_width};
		return rslt;
	}

	public void smoothSensorPolar(
			String [] preTitles,
			String title,
			double [][] preData, // [0] - dx, [1] - dy
			int width,
			int decimate,
			double x0,
			double y0,
			double center_fract, // 0.5 of half-height
			double transit_fract, // 0.2 of half-height - transition from center ortho to outer polar
			double gaus_ang_rel,  // in radians 0.2
			double gaus_rad_rel,  // in fractions of the full radius 0.1
			double max_diff_err_geom,   // before second pass linearly fade R/T and RGB where high-frequency error nears thios value
			double max_diff_err_photo,
			boolean showDebugImages
			) {
		int indexDx=    0;
		int indexDy=    1;
		int indexWeight=2;
		int indexR=     3;
		int indexG=     4;
		int indexB=     5;
		int length=preData[0].length;
		int height=length/width;
		double sensor_radius_d = Math.sqrt(width*width+height*height)/2;
		double gaus_ang =   gaus_ang_rel * sensor_radius_d;
		double gaus_rad =   gaus_rad_rel *  sensor_radius_d;
		int lines_overlap = (int) Math.round(4 * gaus_ang);
		int inner_rad_d = (int) Math.round(height/2 * center_fract);
		int transit_rad = (int) Math.round(height/2 * transit_fract);

		int [] ints = radialInts(
				length,
				width,
				decimate,
				x0,
				y0,
				inner_rad_d); // margins (minimal from all sides)
		int polar_height0 = ints[0];
		int polar_width =   ints[1];
		int polar_height =  polar_height0 + lines_overlap;
		int polar_length = polar_width * polar_height;

// convert dXY to dRT
		double [][] data_xy = {preData[indexDx].clone(), preData[indexDy].clone()}; // clone to later blur (actually src is not needed)
 		double [][] data_rt = DXYtoDRT(
				data_xy, // [0] - dx, [1] - dy
				width,
				decimate,
				x0,
				y0);
		double [][] data_ortho= {
				preData[indexWeight],
				data_rt[0],
				data_rt[1],
				preData[indexR],
				preData[indexG],
				preData[indexB]};

		String [] dbg_titles= {"Weight","dR", "dT","Red", "Green", "Blue"};
		if (showDebugImages) {
			ShowDoubleFloatArrays.showArrays(data_ortho, width, height,  true, "ortho-"+title, dbg_titles);
		}
		double [][] data_polar = correctionToRadial(
				data_ortho, // variable number of layers, [0] - weight dRT or RGB
				width,
				decimate,
				x0,
				y0,
				inner_rad_d,    // distorted inner radis
				lines_overlap, // add lines to the bottom (corresponding to the angle)
				0); // indexWeight);

		if (showDebugImages) {
			ShowDoubleFloatArrays.showArrays(data_polar, polar_width, polar_height,  true, "polar-"+title, dbg_titles);
		}
		double [][] data_polar_blur = new double[data_polar.length][];
		data_polar_blur[0] = data_polar[0].clone();
		for (int i = 1; i < data_polar_blur.length; i++) {
			data_polar_blur[i] = data_polar[i].clone();
			for (int j = 0; j < polar_length; j++ ) {
				data_polar_blur[i][j]*=data_polar[0][j]; // weight (was *= data_polar_blur[0][j] restore?
				if (Double.isNaN(data_polar_blur[i][j])) data_polar_blur[i][j]= 0.0;
			}
			(new DoubleGaussianBlur() ).blurDouble(data_polar_blur[i], polar_width, polar_height, gaus_rad, gaus_ang, 0.01);
		}
		(new DoubleGaussianBlur() ).blurDouble(data_polar_blur[0], polar_width, polar_height, gaus_rad, gaus_ang, 0.01);
		for (int i = 1; i < data_polar_blur.length; i++) {
			for (int j = 0; j < polar_length; j++ ) {
				if (data_polar_blur[0][j] == 0.0) {
					data_polar_blur[i][j] = Double.NaN;
				} else {
					data_polar_blur[i][j] /= data_polar_blur[0][j]; // weight
				}
			}

		}
		if (showDebugImages) {
			ShowDoubleFloatArrays.showArrays(data_polar_blur, polar_width, polar_height,  true, "polar-bl-"+title, dbg_titles);
		}
		double [][] data_polar_diff = new double[data_polar.length][data_polar[0].length];
		for (int i = 0; i < data_polar_diff.length; i++) {
			for (int j = 0; j < polar_length; j++ ) {
				data_polar_diff[i][j] = data_polar[i][j] - data_polar_blur[i][j];
			}
		}
		if (showDebugImages) {
			ShowDoubleFloatArrays.showArrays(data_polar_diff, polar_width, polar_height,  true, "polar-diff-"+title, dbg_titles);
		}

		double [][] errors_gr = new double[5][data_polar[0].length];
		for (int j = 0; j < polar_length; j++) {
			for (int i = 1; i < 3; i++) {
				errors_gr[0][j] += data_polar_diff[i][j]*data_polar_diff[i][j];
			}
			errors_gr[0][j] = Math.sqrt(errors_gr[0][j]);

			for (int i = 3; i < 6; i++) {
				errors_gr[1][j] += data_polar_diff[i][j]*data_polar_diff[i][j];
			}
			errors_gr[1][j] = Math.sqrt(errors_gr[1][j]);
		}
		errors_gr[2] = data_polar[0].clone(); // weight

		// reduce weight where high frequency difference was high, re-run blur
		double [] weight_geom =  data_polar[0].clone();
		double [] weight_photo = data_polar[0].clone();
		for (int i = 0; i <polar_length; i++) {
			if (Double.isNaN(errors_gr[0][i]) || (errors_gr[0][i] > max_diff_err_geom)) {
				weight_geom[i] = 0.0;
			} else {
				weight_geom[i] *= (max_diff_err_geom - errors_gr[0][i]) / max_diff_err_geom;
			}
			if (Double.isNaN(errors_gr[1][i]) || (errors_gr[1][i] > max_diff_err_photo)) {
				weight_photo[i] = 0.0;
			} else {
				weight_photo[i] *= (max_diff_err_photo - errors_gr[1][i])/max_diff_err_photo;
			}
		}
		errors_gr[3] = weight_geom.clone(); // weight
		errors_gr[4] = weight_photo.clone(); // weight
		if (showDebugImages) {
			String [] titles_err = {"geom","photo","weight", "weight_geom","weight_photo"};
			ShowDoubleFloatArrays.showArrays(errors_gr, polar_width, polar_height,  true, "polar-err-"+title, titles_err);
		}

		// geometry
		for (int i = 1; i < 3; i++) { // geom
			data_polar_blur[i] = data_polar[i].clone();
			for (int j = 0; j < polar_length; j++ ) {
				data_polar_blur[i][j] *= weight_geom[j]; // weight
				if (Double.isNaN(data_polar_blur[i][j])) data_polar_blur[i][j]= 0.0;
			}
			(new DoubleGaussianBlur() ).blurDouble(data_polar_blur[i], polar_width, polar_height, gaus_rad, gaus_ang, 0.01);
		}
		(new DoubleGaussianBlur() ).blurDouble(weight_geom, polar_width, polar_height, gaus_rad, gaus_ang, 0.01);
		for (int i = 1; i < 3; i++) { // geom
			for (int j = 0; j < polar_length; j++ ) {
				if (weight_geom[j] == 0.0) {
					data_polar_blur[i][j] = Double.NaN;
				} else {
					data_polar_blur[i][j] /= weight_geom[j]; // weight
				}
			}
		}

		// photo
		for (int i = 3; i < 6; i++) { // photo
			data_polar_blur[i] = data_polar[i].clone();
			for (int j = 0; j < polar_length; j++ ) {
				data_polar_blur[i][j] *= weight_photo[j]; // weight
				if (Double.isNaN(data_polar_blur[i][j])) data_polar_blur[i][j]= 0.0;
			}
			(new DoubleGaussianBlur() ).blurDouble(data_polar_blur[i], polar_width, polar_height, gaus_rad, gaus_ang, 0.01);
		}
		(new DoubleGaussianBlur() ).blurDouble(weight_photo, polar_width, polar_height, gaus_rad, gaus_ang, 0.01);
		for (int i = 3; i < 6; i++) { // photo
			for (int j = 0; j < polar_length; j++ ) {
				if (weight_photo[j] == 0.0) {
					data_polar_blur[i][j] = Double.NaN;
				} else {
					data_polar_blur[i][j] /= weight_photo[j]; // weight
				}
			}
		}
		data_polar_blur[0] = weight_geom.clone();
		if (showDebugImages) {
			ShowDoubleFloatArrays.showArrays(data_polar_blur, polar_width, polar_height,  true, "polar-bl2-"+title, dbg_titles);
		}
		data_polar_diff = new double[data_polar.length][polar_length];
		for (int i = 0; i < data_polar_diff.length; i++) {
			for (int j = 0; j < polar_length; j++ ) {
				data_polar_diff[i][j] = data_polar[i][j] - data_polar_blur[i][j];
			}
		}
		if (showDebugImages) {
			ShowDoubleFloatArrays.showArrays(data_polar_diff, polar_width, polar_height,  true, "polar-diff2-"+title, dbg_titles);
			ShowDoubleFloatArrays.showArrays(data_ortho, width, height,  true, "ortho-back-pre-"+title, dbg_titles);
		}

		// put back to ortho (so far no transition, no ortho blur)
		// Blur dXY, not dRT!

		for (int l = 0; l < 2; l++) {
			(new DoubleGaussianBlur() ).blurDouble(data_xy[l], width, height, gaus_rad, gaus_rad, 0.01);
		}
		data_rt = DXYtoDRT(
				data_xy, // [0] - dx, [1] - dy
				width,
				decimate,
				x0,
				y0);
		for (int l = 0; l < 2; l++) {
			data_ortho[l+1] = data_rt[l];
		}


		double [][] polar_back = {
				data_polar[0], // weight
				data_polar_blur[1],
				data_polar_blur[2],
				data_polar_blur[3],
				data_polar_blur[4],
				data_polar_blur[5]};


		(new DoubleGaussianBlur() ).blurDouble(data_ortho[0], width, height, gaus_rad, gaus_rad, 0.01); // is it needed - blur alpha?

		for (int l = 3; l < data_ortho.length; l++) { // only R,G,B
			(new DoubleGaussianBlur() ).blurDouble(data_ortho[l], width, height, gaus_rad, gaus_rad, 0.01);
		}

//		if (showDebugImages) {
//			ShowDoubleFloatArrays.showArrays(polar_back, polar_width, polar_height,  true, "polar-back-"+title, dbg_titles);
//		}

		correctionToOrtho(
				data_ortho, // variable number of layers, will be modifier from the overlapping radial)
				polar_back, // same layers as ortho, will be used as a source
				width,
				decimate,
				x0,
				y0,
				inner_rad_d, // margins (minimal from all sides)
				lines_overlap, // add lines to the bottom (corresponding to the angle)
				transit_rad); // smooth transition between ortho and polar width
		if (showDebugImages) {
			ShowDoubleFloatArrays.showArrays(data_ortho, width, height,  true, "ortho-back"+title, dbg_titles);
		}

		double [][] data_rt_back = {data_ortho[1], data_ortho[2]};
		double [][] data_xy_back = DRTtoDXY(
				data_rt_back, // [0] - dx, [1] - dy
				width,
				decimate,
				x0,
				y0);

		preData[indexWeight] = data_ortho[0];
		preData[indexDx] =     data_xy_back[0]; // data_ortho[1];
		preData[indexDy] =     data_xy_back[1]; // data_ortho[2];
		preData[indexR] =      data_ortho[3];
		preData[indexG] =      data_ortho[4];
		preData[indexB] =      data_ortho[5];

	}

	public void addOldXYCorrectionToCurrent(
//			boolean invert,
    		RefineParameters refineParameters,
    		double [][][] sensorXYCorr
			){
		if (this.pixelCorrection==null) return; // no modifications are needed
		for (int i=0;i<sensorXYCorr.length;i++) if ((sensorXYCorr[i]!=null) && (this.pixelCorrection[i]!=null)) {
//			boolean small_sensor = fittingStrategy.distortionCalibrationData.isSmallSensor(i);
			boolean small_sensor = fittingStrategy.distortionCalibrationData.getSmallSensors()[i];
			
			RefineParameters rp = small_sensor ? refineParameters.refineParametersSmall : refineParameters;
//			if (rp.sensorExtrapolateDiff ^ invert) { // add current correction AFTER extrapolationg/bluring
				double scale = rp.correctionScale;
				for (int j=0;j<sensorXYCorr[i][0].length;j++){
					sensorXYCorr[i][0][j]=this.pixelCorrection[i][0][j]+scale*sensorXYCorr[i][0][j];
					sensorXYCorr[i][1][j]=this.pixelCorrection[i][1][j]+scale*sensorXYCorr[i][1][j];
				}
			}
//		}
	}



	public void patternErrors(
			final int       threadsMax,
			final boolean   updateStatus,
			final int debugLevel
			){
		GenericDialog gd=new GenericDialog("Setup pattern errors map");
		gd.addNumericField("Series number", this.seriesNumber, 0,2,"");
		gd.addCheckbox    ("Show map", true);

		gd.addNumericField("Minimal RMS", .07, 3,6,"pix");
		gd.addNumericField("Maximal RMS", 0.12, 3,6,"pix");
		gd.addNumericField("Expand EMS mask", 1, 0,2,"nodes");
		gd.addCheckbox    ("Update pattern weights", false);
		gd.addCheckbox    ("Reset error-based target map", false);
		gd.showDialog();
		if (gd.wasCanceled()) return;
		this.seriesNumber =      (int) gd.getNextNumber();
		boolean showMap=               gd.getNextBoolean();
 		double minRMS =                gd.getNextNumber();
 		double maxRMS =                gd.getNextNumber();
		int expandMask =         (int) gd.getNextNumber();

		boolean updateMap=              gd.getNextBoolean();
		boolean resetMap=              gd.getNextBoolean();

		if (resetMap){
			this.patternParameters.resetPatternErrorMask();
			return;
		} else {
			double [] worstImageNumber=calculatePatterErrorRMS(
					this.seriesNumber,
					threadsMax,
					updateStatus,
					debugLevel);
			this.patternParameters.savePatternErrorMask();
			double [] savedMask=this.patternParameters.getSavedPatternErrorMask();
			this.patternParameters.calculatePatternErrorMask(maxRMS,minRMS);
			for (int i=0;i<expandMask;i++)this.patternParameters.expandPatternErrorMask();
			if (showMap){
				String [] titles={"mask","rms","worst image number","savedMask"};
				double [][] debugData={
						this.patternParameters.getPatternErrorMask(),
						this.patternParameters.getPatternErrors(),
						worstImageNumber,
						savedMask};
				 Rectangle gridDimensions=patternParameters.getUVDimensions();
				ShowDoubleFloatArrays.showArrays(
						debugData,
						gridDimensions.width,
						gridDimensions.height,
						true,
						"TM_"+maxRMS+":"+minRMS,
						titles);
			}
			if (!updateMap) {
				System.out.println("Restoring mask to the previous state");
				this.patternParameters.restorePatternErrorMask();
			}
		}
	}



	public double []  calculatePatterErrorRMS( // returns worst image number array
			final int       series,
			final int       threadsMax,
			final boolean   updateStatus,
			final int debugLevel

	){
    	if (fittingStrategy==null) {
    		String msg="Fitting strategy does not exist, exiting";
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
    	}
    	if (fittingStrategy.distortionCalibrationData.eyesisCameraParameters==null){
    		String msg="Eyesis camera parameters (and sensor dimensions) are not defined";
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
    	}
    	//	fittingStrategy.distortionCalibrationData.readAllGrids();
//    	if (! selectGridEnhanceParameters()) return false;
//    	if (series<0) return null; // false; // make "all " later?
    	this.seriesNumber=series;

    	initFittingSeries(true,this.filterForTargetGeometry,this.seriesNumber); // first step in series now uses pattern alpha
    	this.currentfX=calculateFxAndJacobian(this.currentVector, false);
    	//        	this.currentRMS= calcError(calcYminusFx(this.currentfX));
    	if (this.debugLevel>2) {
    		System.out.println("this.currentVector");
    		for (int i=0;i<this.currentVector.length;i++){
    			System.out.println(i+": "+ this.currentVector[i]);
    		}
    	}
		final boolean [] selectedImages=fittingStrategy.selectedImages();
		final Rectangle gridDimensions=patternParameters.getUVDimensions();
		final int width=  gridDimensions.width;
		final int height= gridDimensions.height;
//		final int U0=     gridDimensions.x;
//		final int V0=     gridDimensions.y;
		final double [][] gridErrors=new double [4][width*height]; // added debug features - worst image number
		for (int n=0;n<gridErrors.length;n++) for (int i=0;i<gridErrors[n].length;i++) gridErrors[n][i]=0.0;
		int numSelected=0;
		for (int imgNum=0;imgNum<selectedImages.length;imgNum++) if (selectedImages[imgNum]) numSelected++;
		final int finalSelected=numSelected;
		if (updateStatus) IJ.showStatus("Calculating pattern grid errors...");
   		final AtomicInteger imageNumberAtomic = new AtomicInteger(0);
   		final AtomicInteger imageFinishedAtomic = new AtomicInteger(0);
   		final Thread[] threads = newThreadArray(threadsMax);
   		for (int ithread = 0; ithread < threads.length; ithread++) {
   			threads[ithread] = new Thread() {
   				@Override
				public void run() {
   					double [][] partialGridErrors=new double [4][width*height];
   					for (int n=0;n<partialGridErrors.length;n++) for (int i=0;i<partialGridErrors[n].length;i++) partialGridErrors[n][i]=0.0;
   					for (int imgNum=imageNumberAtomic.getAndIncrement(); imgNum<selectedImages.length;imgNum=imageNumberAtomic.getAndIncrement()){
   						if (selectedImages[imgNum]){
   							accumulatePatternErrors(
   									partialGridErrors,
   									imgNum,
   									gridDimensions);
   							final int numFinished=imageFinishedAtomic.getAndIncrement();
   							SwingUtilities.invokeLater(new Runnable() {
   								@Override
								public void run() {
   									if (updateStatus) IJ.showProgress(numFinished,finalSelected);
   								}
   							});
   						} //if (selectedImages[numImage]){
   					} // for (int numImage=imageNumberAtomic.getAndIncrement(); ...
   					combinePatternErrors(partialGridErrors,gridErrors);
   				} // public void run() {
   			};
   		}
   		startAndJoin(threads);
   		for (int i=0;i<gridErrors[0].length;i++){
   			gridErrors[0][i]=(gridErrors[0][i]>0.0)?Math.sqrt(gridErrors[0][i]/gridErrors[1][i]):Double.NaN;

   		}
   		patternParameters.setPatternErrors(gridErrors[0]);
   		return gridErrors[2]; // worst image number for target grid nodes
	}

	public void accumulatePatternErrors(
			double [][] errorMap,
			int imgNum,
			Rectangle gridDimensions){
		int width=  gridDimensions.width;
//		int height= gridDimensions.height;
		int U0=     gridDimensions.x; // location of the grid center (U==0,V==0)
		int V0=     gridDimensions.y;
		double [] diff=calcYminusFx(this.currentfX, 2*this.imageStartIndex[imgNum],2*this.imageStartIndex[imgNum+1]);
		int [][] imgUV=	  this.fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV;
		for (int i=0;i<imgUV.length;i++){
			int index=width*(imgUV[i][1]+V0) + (imgUV[i][0]+U0);
			double w=this.weightFunction[2*(this.imageStartIndex[imgNum]+i)];
			double dX=diff[2*i];
			double dY=diff[2*i+1];
			double e2w=w*(dX*dX+dY*dY);
			errorMap[0][index]+=e2w;
			errorMap[1][index]+=w;
			if (e2w>errorMap[3][index]){
				errorMap[3][index]=e2w;    // worst error for this node
				errorMap[2][index]=imgNum; // worst (for that particular grig node) image number
			}
		}
	}

	public synchronized void combinePatternErrors(
			double [][] partialErrorMap,
			double [][] fullErrorMap ){
//		for (int n=0;n<fullErrorMap.length;n++) for (int i=0;i<fullErrorMap[n].length;i++) fullErrorMap[n][i]+=partialErrorMap[n][i];
		for (int i=0;i<fullErrorMap[0].length;i++){
			fullErrorMap[0][i]+=partialErrorMap[0][i];
			fullErrorMap[1][i]+=partialErrorMap[1][i];
			if (fullErrorMap[3][i]<partialErrorMap[3][i]){
				fullErrorMap[2][i]=partialErrorMap[2][i];
				fullErrorMap[3][i]=partialErrorMap[3][i];
			}

		}

	}




	/**
	 * Calculate each sensor correction increment for geometry and photometry contributed by all images selected in a series
	 * @param selectedImages process only selected images
	 * @param showIndividual show per-image intermediate results
	 * @param threadsMax maximal number of concurrent threads
	 * @param updateStatus update IJ status/progress
	 * @param debugLevel debug level
	 * @return [sensor]{dpX,dpY,alpha,R,G,B}[pixelIndex] . dpX, dpY - correction to previous, RGB - total FF, not increment!
	 */

	public double [][][]  allImagesCorrectionMapped(
			final boolean []       selectedImages,
			final boolean          si,               // showIndividual,
			final RefineParameters refineParameters, // final int showIndividualNumber,
			final int              threadsMax,
			final boolean          updateStatus,
			final int debugLevel
			){
		int numChannels=  fittingStrategy.distortionCalibrationData.getNumChannels(); // number of used channels
		final double [][][] gridPCorr=new double [numChannels][][];
		for (int chnNum=0;chnNum<gridPCorr.length;chnNum++) gridPCorr[chnNum]=null;
		int numSelected=0;
		for (int imgNum=0;imgNum<selectedImages.length;imgNum++) if (selectedImages[imgNum]) numSelected++;
		final int finalSelected=numSelected;
		if (updateStatus) IJ.showStatus("Calculating sensor corrections...");
   		final AtomicInteger imageNumberAtomic = new AtomicInteger(0);
   		final AtomicInteger imageFinishedAtomic = new AtomicInteger(0);
   		final Thread[] threads = newThreadArray(threadsMax);
   		final AtomicInteger stopRequested=this.stopRequested;
		final AtomicBoolean interruptedAtomic=new AtomicBoolean();
		final int alphaIndex=2;

   		for (int ithread = 0; ithread < threads.length; ithread++) {
   			threads[ithread] = new Thread() {
   				@Override
				public void run() {
   					for (int imgNum=imageNumberAtomic.getAndIncrement(); (imgNum<selectedImages.length) && !interruptedAtomic.get();imgNum=imageNumberAtomic.getAndIncrement()){
   						if (selectedImages[imgNum]){
   							int chnNum=fittingStrategy.distortionCalibrationData.gIP[imgNum].channel; // number of sub-camera
   	   						//boolean small_sensor = fittingStrategy.distortionCalibrationData.isSmallSensor(chnNum);
   	   						boolean small_sensor = fittingStrategy.distortionCalibrationData.getSmallSensors()[chnNum];
   	   						RefineParameters rp = small_sensor ? refineParameters.refineParametersSmall : refineParameters;
   	   						boolean   showIndividual = si && rp.showPerImage;
   	   						int showIndividualNumber = rp.showIndividualNumber;
   							double [][] singleCorr=
   								singleImageCorrectionMapped(
   									imgNum, // image number
   									showIndividual && ((showIndividualNumber<0) || (showIndividualNumber==chnNum)),
   									debugLevel);
   							combineImageCorrection(
   									chnNum,
   									gridPCorr,
   									singleCorr
   							);
   							final int numFinished=imageFinishedAtomic.getAndIncrement();
   							SwingUtilities.invokeLater(new Runnable() {
   								@Override
								public void run() {
   									if (updateStatus) IJ.showProgress(numFinished,finalSelected);
   								}
   							});
   	   						if (stopRequested.get()==1){ // ASAP
   	   							interruptedAtomic.set(true);
   	   						}
   						} //if (selectedImages[numImage]){
   					} // for (int numImage=imageNumberAtomic.getAndIncrement(); ...
   				} // public void run() {
   			};
   		}
   		startAndJoin(threads);
   		// divide by weight;
   		for (int nChn=0;nChn<gridPCorr.length;nChn++) if (gridPCorr[nChn]!=null){
   			for (int i=0;i<gridPCorr[nChn].length;i++) {
   				if (i!=alphaIndex){
   					for (int j=0; j<gridPCorr[nChn][i].length;j++){
   						if (gridPCorr[nChn][alphaIndex][j]>0) gridPCorr[nChn][i][j]/=gridPCorr[nChn][alphaIndex][j];
   					}
   				}
   			}
   		}

		if (updateStatus) IJ.showProgress(0);

   		if (interruptedAtomic.get()) {
   			System.out.println("allImagesCorrection() aborted by user request");
   			return null;
   		}
   		return gridPCorr;
	}
	@Deprecated
	public void allSensorsExtrapolationMapped(
			final int stationNumber, // has to be selected
			final double [][][] gridPCorr,
			final RefineParameters refineParameters, //
			final int       threadsMax,
			final boolean   updateStatus,
			final boolean showDebugImages,
			final int debugLevel
			){
		if (updateStatus) IJ.showStatus("Extrapolating sensor corrections...");
   		final AtomicInteger sensorNumberAtomic = new AtomicInteger(0);
   		final AtomicInteger sensorFinishedAtomic = new AtomicInteger(0);
   		final Thread[] threads = newThreadArray(threadsMax);
   		final AtomicInteger stopRequested=this.stopRequested;
		final AtomicBoolean interruptedAtomic=new AtomicBoolean();
		final EyesisSubCameraParameters [] eyesisSubCameras = this.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[stationNumber];
		final double [][] sensorMasks=this.fittingStrategy.distortionCalibrationData.sensorMasks;

		final int alphaIndex=2;

		final boolean extraShowDebug=showDebugImages&& (debugLevel>2);

   		for (int ithread = 0; ithread < threads.length; ithread++) {
   			threads[ithread] = new Thread() {
   				@Override
				public void run() {
   					DoubleGaussianBlur gb=null;
   					double [][] debugMasks1=null;
   					double [][] debugMasks2=null;
   					String [] debugMaskTitles={"original","blured"};
   					if (extraShowDebug){
   						debugMasks1=new double[2][];
   						debugMasks2=new double[2][];
   					}
//   					if (shrinkBlurComboSigma>0.0) gb=new DoubleGaussianBlur();
   					gb=new DoubleGaussianBlur();
   					for (int sensorNum=sensorNumberAtomic.getAndIncrement(); (sensorNum<gridPCorr.length) && !interruptedAtomic.get();sensorNum=sensorNumberAtomic.getAndIncrement()){
//   						boolean small_sensor = fittingStrategy.distortionCalibrationData.isSmallSensor(sensorNum);
   						boolean small_sensor = fittingStrategy.distortionCalibrationData.getSmallSensors()[sensorNum];
   						RefineParameters rp = small_sensor ? refineParameters.refineParametersSmall : refineParameters;
   						if (rp.extrapolate) {
   							int sensorWidth=   fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getSensorWidth(sensorNum);
   							int sensorHeight=  fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getSensorHeight(sensorNum);
   							int decimation=fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getDecimateMasks(sensorNum);
   							int width= (sensorWidth-1)/decimation+1; // decimated width (648)
   							int height= (sensorHeight-1)/decimation+1; // decimated width (648)

   							if (gridPCorr[sensorNum]!=null){
   								final double [] centerPXY={
   										eyesisSubCameras[sensorNum].px0,
   										eyesisSubCameras[sensorNum].py0
   								};
   								if (rp.sensorShrinkBlurComboSigma>0.0){
   									double sigma=rp.sensorShrinkBlurComboSigma/decimation;
   									int margin=(int) (2*sigma);
   									int width1=width+2*margin;
   									int height1=height+2*margin;
   									if (extraShowDebug) debugMasks2[0]=gridPCorr[sensorNum][alphaIndex].clone();
   									double [] mask= addMarginsThreshold(
   											gridPCorr[sensorNum][alphaIndex], // double [] data,
   											0.0, // double threshold,
   											width,
   											height,
   											margin);
   									if (extraShowDebug) debugMasks1[0]=mask.clone();
   									gb.blurDouble(
   											mask,
   											width1,
   											height1,
   											sigma,
   											sigma,
   											0.01);

   									double k=1.0/(1.0-rp.sensorShrinkBlurComboLevel);
   									for (int i=0;i<mask.length;i++) {
   										mask[i]=k*(mask[i]-rp.sensorShrinkBlurComboLevel);
   										mask[i]=(mask[i]>0.0)?(mask[i]*mask[i]):0.0;
   									}
   									if (extraShowDebug) debugMasks1[1]=mask.clone();
   									gridPCorr[sensorNum][alphaIndex]=removeMargins(
   											mask, //double [] data,
   											width, // w/o margins
   											height,
   											margin); //mask; // replace with 0.0 .. 1.0 mask
   									if (extraShowDebug) debugMasks2[1]=gridPCorr[sensorNum][alphaIndex].clone();
   									if (extraShowDebug) {
   										ShowDoubleFloatArrays.showArrays(
   												debugMasks1,
   												width1,
   												height1,
   												true,
   												"M1-"+sensorNum,
   												debugMaskTitles);
   										ShowDoubleFloatArrays.showArrays(
   												debugMasks2,
   												width,
   												height,
   												true,
   												"M2-"+sensorNum,
   												debugMaskTitles);
   									}

   								}
   								singleSensorExtrapolationMapped(
   										sensorNum,
   										gridPCorr[sensorNum],
   										sensorMasks[sensorNum],
   										width,
   										decimation,
   										rp.sensorAlphaThreshold,
   										rp.sensorStep,
   										centerPXY,
   										rp.sensorInterpolationSigma,
   										rp.sensorTangentialRadius,
   										rp.sensorScanDistance,
   										rp.sensorResultDistance,
   										rp.sensorInterpolationDegree,
   										(rp.sensorShrinkBlurComboSigma > 0.0),
   										showDebugImages,
   										debugLevel);
   								final int numFinished=sensorFinishedAtomic.getAndIncrement();
   								SwingUtilities.invokeLater(new Runnable() {
   									@Override
   									public void run() {
   										if (updateStatus) IJ.showProgress(numFinished,gridPCorr.length);
   									}
   								});
   								if (stopRequested.get()==1){ // ASAP
   									interruptedAtomic.set(true);
   								}
   							}
   						} //if (rp.extrapolate); // if (selectedImages[numImage]){
   					} // for (int numImage=imageNumberAtomic.getAndIncrement(); ...
   				} // public void run() {
   			};
   		}
   		startAndJoin(threads);
		if (updateStatus) IJ.showProgress(0);
   		return;
	}
 @Deprecated
 	public double [] addMargins(
			double [] data,
			double marginData,
			int width,
			int height,
			int margin){
		int width1= width+ 2*margin;
		int height1=height+2*margin;
		int length1=width1*height1;
		double [] result = new double [length1];
		for (int i=0;i<length1;i++) result[i] = marginData;
		int indexDest=margin*(width1+1);
		int indexSrc=0;
		for (int y=0;y<height;y++){
			for (int x=0;x<width;x++){
				result[indexDest++]=data[indexSrc++];
			}
			indexDest+=2*margin;
		}
		return result;
	}
    @Deprecated
	public double [] addMarginsThreshold(
			double [] data,
			double threshold,
			int width,
			int height,
			int margin){
		int width1= width+ 2*margin;
		int height1=height+2*margin;
		int length1=width1*height1;
		double [] result = new double [length1];
		for (int i=0;i<length1;i++) result[i] = -1.0;
		int indexDest=margin*(width1+1);
		int indexSrc=0;
		for (int y=0;y<height;y++){
			for (int x=0;x<width;x++){
				result[indexDest++]=(data[indexSrc++]>threshold)?1.0:-1.0;
			}
			indexDest+=2*margin;
		}
		return result;
	}
    @Deprecated
	public double [] removeMargins(
			double [] data,
			int width, // w/o margins
			int height,
			int margin){
		int width1= width+ 2*margin;
//		int height1=height+2*margin;
		int length=width*height;
		double [] result = new double [length];
		int indexSrc=margin*(width1+1);
		int indexDest=0;
		for (int y=0;y<height;y++){
			for (int x=0;x<width;x++){
				result[indexDest++]=data[indexSrc++];
			}
			indexSrc+=2*margin;
		}
		return result;
	}
    @Deprecated
	public void singleSensorExtrapolationMapped(
			int sensoNum,
			double [][] gridPCorr,
			double [] sensorMask,
			int width,
			int decimation,
			double alphaThreshold,
			double step,
			double [] centerPXY,
			double interpolationSigma, // sensor pixels
			double tangentialRadius,
			int    scanDistance,       // sensor pixels
			int resultDistance,
			int interpolationDegree,
			boolean useAlpha, // false - sensor mask
			boolean showDebugImages,
			int debugLevel
			){
		int dxIndex=0;
		int alphaIndex=2;
		int rIndex=3;
		int height=gridPCorr[0].length/width;
		double gaussianK=-0.5/(interpolationSigma*interpolationSigma);
		double tangR0=tangentialRadius*Math.sqrt(width*height)*decimation/2; // sigma in tangential direction is interpolationSigma*(1+r/tangR0), in radial - interpolationSigma
		PolynomialApproximation polynomialApproximation =new PolynomialApproximation(0);// no debug
		int length=gridPCorr[0].length;
		DirInc dirInc= new DirInc(width,height);
		int [] iMap = new int[length];
		for (int i=0;i<length;i++) iMap[i]= (gridPCorr[alphaIndex][i]>=alphaThreshold)?1:0;
		List <Integer>waveList=new ArrayList<Integer>(1000);
		for (int index0=0;index0<length;index0++) if (iMap[index0]==0){
			for (int iDir=0;iDir<8;iDir+=2){
				int index=dirInc.newIndex(index0,iDir);
				if ((index>=0) && (iMap[index]==1)){
					iMap[index0]=2;
					waveList.add(new Integer(index0));
					break;
				}
			}
		}
// decimate the wave list
		List <Integer> seedList=new ArrayList<Integer>(1000);
		int oldIndex=0; // find better start?
		int s2= (int) Math.floor(step*step);
		while (waveList.size()>0){
			int oldX=oldIndex%width;
			int oldY=oldIndex/width;
			int bestD2=height*height+width*width;
			int nBest=-1;
			for (int n=0;n<waveList.size();n++){
				int index=waveList.get(n);
				int dx=index%width-oldX;
				int dy=index/width-oldY;
				int d2=dx*dx+dy*dy;
				if (d2<bestD2){
					bestD2=d2;
					nBest=n;
				}
			}
			oldIndex=waveList.remove(nBest);
			seedList.add(new Integer(oldIndex));
			oldX=oldIndex%width;
			oldY=oldIndex/width;
			// remove all closer than step
			for (int n=0;n<waveList.size();n++){ // size will change
				int index=waveList.get(n);
				int dx=index%width-oldX;
				int dy=index/width-oldY;
				int d2=dx*dx+dy*dy;
				if (d2<s2){
					waveList.remove(n);
				}
			}

		} //while (waveList.size()>0)
		// debug show waves?
		Rectangle full=new Rectangle (0,0,width,height);
		double [][] extrapolated=new double [gridPCorr.length][length];
		for (int n=0;n<extrapolated.length;n++) for (int i=0;i<extrapolated[n].length;i++) extrapolated[n][i]=0.0;
		int halfScanSize=scanDistance/decimation+1;
		int halfInterpolteSize=resultDistance/decimation+1;
		for (int n=0; n<seedList.size();n++) {
			int index0=seedList.get(n);
			int x0=index0%width;
			int y0=index0/width;
			double [] dCxy0={
					x0*decimation-centerPXY[0],
					y0*decimation-centerPXY[1]
			};
			double r0=Math.sqrt(dCxy0[0]*dCxy0[0]+dCxy0[1]*dCxy0[1]);
			final Rectangle scan =full.intersection(new Rectangle (x0-halfScanSize,y0-halfScanSize,2*halfScanSize+1,2*halfScanSize+1));
			waveList.clear();
			for (int y=scan.y;y<(scan.y+scan.height);y++) for (int x=scan.x;x<(scan.x+scan.width);x++) {
				int index=y*width+x;
				if (iMap[index]==1)	waveList.add(new Integer(index));
			}
			double [][][] data = new double [5][waveList.size()][3]; // x,y,w
			double sumWeights=0.0;
			double rScaleTangSigma=1.0/(1.0+r0/tangR0); //
			for (int i=0;i<data[0].length;i++){
				int index=waveList.get(i);
				int x=index%width;
				int y=index/width;
				double [] dCxy={
						x*decimation-centerPXY[0],
						y*decimation-centerPXY[1]
				};
				double [] ddCxy={
						dCxy[0]-dCxy0[0],
						dCxy[1]-dCxy0[1]
				};
				double rc=Math.sqrt(dCxy[0]*dCxy[0]+dCxy[1]*dCxy[1]); // distance from lens center (in sensor pixels)
				double rDiff=rc-r0;
				double [] uRadVect={(rc>0.0)?(dCxy[0]/rc):0.0, (rc>0.0)?(dCxy[1]/rc):0.0};

				double distRad= ddCxy[0]*uRadVect[0]+ddCxy[1]*uRadVect[1]; // radial distance form the center (seed point)
				double distTan=-ddCxy[0]*uRadVect[1]+ddCxy[1]*uRadVect[0]; // tangential distance form the center (seed point)
				distTan*=rScaleTangSigma; // // for the center (seed point). was  distTan/=(1.0+rc/tangR0);
				double w=Math.exp(gaussianK*(distRad*distRad+distTan*distTan))*gridPCorr[alphaIndex][index];
				sumWeights+=w;

				double dRad= gridPCorr[dxIndex+0][index]*uRadVect[0]+gridPCorr[dxIndex+1][index]*uRadVect[1]; // radial component
				double dTan=-gridPCorr[dxIndex+0][index]*uRadVect[1]+gridPCorr[dxIndex+1][index]*uRadVect[0]; // tangential component
				data[0][i][1]=dRad;
				data[1][i][1]=dTan;

				data[2][i][1]=gridPCorr[rIndex+0][index]; // R
				data[3][i][1]=gridPCorr[rIndex+1][index]; // G
				data[4][i][1]=gridPCorr[rIndex+2][index]; // B
				for (int j=0;j<data.length;j++){
					data[j][i][0]=rDiff;
					data[j][i][2]=w;
				}
			}
			sumWeights*=rScaleTangSigma; // normalize for expanded in one dimension gaussian
			double [][] poly=new double [data.length][];
			for (int j=0;j<poly.length;j++) {
				poly[j]=polynomialApproximation.polynomialApproximation1d(data[j],interpolationDegree);
			}
			if (poly[0]==null) { // all will be either null, or not - [0] testing is enough
				System.out.println("singleSensorExtrapolationMapped() BUG - poly[0]==null");
//				stageReprojPXY[index0]=null;
				continue;
			}
			final Rectangle rInterpolate =full.intersection(new Rectangle (x0-halfInterpolteSize,y0-halfInterpolteSize,2*halfInterpolteSize+1,2*halfInterpolteSize+1));
			for (int y=rInterpolate.y;y<(rInterpolate.y+rInterpolate.height);y++) for (int x=rInterpolate.x;x<(rInterpolate.x+rInterpolate.width);x++) {
				int index=y*width+x;
				double [] dCxy={
						x*decimation-centerPXY[0],
						y*decimation-centerPXY[1]
				};
				double [] ddCxy={
						dCxy[0]-dCxy0[0],
						dCxy[1]-dCxy0[1]
				};
				double rc=Math.sqrt(dCxy[0]*dCxy[0]+dCxy[1]*dCxy[1]); // distance from lens center (in sensor pixels)
				double rDiff=rc-r0;
				double [] uRadVect={(rc>0.0)?(dCxy[0]/rc):0.0, (rc>0.0)?(dCxy[1]/rc):0.0};

				double distRad= ddCxy[0]*uRadVect[0]+ddCxy[1]*uRadVect[1]; // radial distance form the center (seed point)
				double distTan=-ddCxy[0]*uRadVect[1]+ddCxy[1]*uRadVect[0]; // tangential distance form the center (seed point)
				distTan*=rScaleTangSigma;
				double w=Math.exp(gaussianK*(distRad*distRad+distTan*distTan)); //*gridPCorr[alphaIndex][index];
				w*=sumWeights; // more points were used in coefficients calculation, more trust to that extrapolation
				// extrapolate each value using polynomial coefficients
				double [] results= new double [poly.length];
				for (int nPar=0;nPar<results.length;nPar++){
					double rN=1.0;
					results[nPar]=0.0;
					for (int dgr=0;dgr<poly[nPar].length;dgr++){
						results[nPar]+=poly[nPar][dgr]*rN;
						rN*=rDiff;
					}
				}
				// restore dX, dY from radial/tangential
				double [] diffPXY={
						results[0]*uRadVect[0]-results[1]*uRadVect[1],
						results[0]*uRadVect[1]+results[1]*uRadVect[0]};
                //accumulate
				extrapolated[dxIndex+0][index]+=diffPXY[0]*w;
				extrapolated[dxIndex+1][index]+=diffPXY[1]*w;
				extrapolated[rIndex+0][index]+=results[2]*w;
				extrapolated[rIndex+1][index]+=results[3]*w;
				extrapolated[rIndex+2][index]+=results[4]*w;
				extrapolated[alphaIndex][index]+=w;
			}
		} // for (int n=0; n<seedList.size();n++) {
		// divide by weight
		for (int index=0;index<length;index++) if (extrapolated[alphaIndex][index]>0.0){
			for (int i=0;i<extrapolated.length;i++) if (i!=alphaIndex){
				extrapolated[i][index]/=extrapolated[alphaIndex][index];
			}
		}
		// debug show here extrapolated
		if (showDebugImages){
			String [] debugTiles={"dX","dY","alpha","R","G","B","mask"};
			double [] debugMask=new double[length];
			for (int i=0;i<length;i++) debugMask[i]=iMap[i];
			for (int n=0; n<seedList.size();n++) {
				int index=seedList.get(n);
				debugMask[index]+=3.0;
			}
			//iMap[index0]
			double [][] debugData={
					extrapolated[0],
					extrapolated[1],
					extrapolated[2],
					extrapolated[3],
					extrapolated[4],
					extrapolated[5],
					debugMask};
			ShowDoubleFloatArrays.showArrays(
					debugData,
					width,
					height,
					true,
					"EX-"+sensoNum,
					debugTiles);

		}
		// mix interpolated with original data
// double [] sensorMask,
//gridPCorr
		for (int index=0;index<length;index++) if (extrapolated[alphaIndex][index]>0.0){
			for (int i=0;i<extrapolated.length;i++) if (i!=alphaIndex){
				double w=useAlpha?(gridPCorr[alphaIndex][index]):((gridPCorr[alphaIndex][index]>0.0)?sensorMask[index]:0.0);
				gridPCorr[i][index]=gridPCorr[i][index]*w+extrapolated[i][index]*(1.0-w);
			}
		}
	}
    // used 2020
	public synchronized void combineImageCorrection(
			int chnNum,
			double [][][] gridPCorr,
			double [][] singleCorr
	){
		int alphaIndex=2;

		if (gridPCorr[chnNum]==null){
			gridPCorr[chnNum]=new double [singleCorr.length][singleCorr[0].length];
			for (int i=0;i<singleCorr.length;i++) for (int j=0; j<singleCorr[i].length;j++){
				gridPCorr[chnNum][i][j]=0.0;
			}
		}
		for (int i=0;i<singleCorr.length;i++) {
			if (i==alphaIndex){
				for (int j=0; j<singleCorr[i].length;j++) gridPCorr[chnNum][i][j]+=singleCorr[i][j];
			} else {
				for (int j=0; j<singleCorr[i].length;j++) gridPCorr[chnNum][i][j]+=singleCorr[i][j]*singleCorr[alphaIndex][j];
			}
		}
	}

	/**
	 * Calculate sensor correction increment for geometry and photometry contributed by a single image
	 * @param imgNum  number of image
	 * @param maxSensorMask maximal value of the sensor mask for this sensor to start extrapolating
	 * @param minContrast minimal measured grid contrast to seed extrapolating  - to prevent expansion in the areas where this particular sensor has bad data
	 * @param minTargetAlpha - minimal alpha of the target node
	 * @param useTargetAlpha   false - only use contrast of the detected grid, true - multiply contrast by grid alpha
	 * @param showIntermediate - show intermediate data as images
	 * @param debugLevel debug level
	 * @return scan-line pixels additional correction arrays {dpX,dpY,alpha,R,G,B}[pixelIndex]
	 */
	public double [][]  singleImageCorrectionMapped(
			int imgNum, // image number
			boolean showIntermediate,
			int debugLevel
			){
		CorrectionInNodes correctionInNodes=extractNodeCorrections(
				imgNum, // image number
				showIntermediate,
				debugLevel);
		if (showIntermediate) correctionInNodes.show("finNode-"+imgNum);
		int chnNum=fittingStrategy.distortionCalibrationData.gIP[imgNum].channel;
		int sensorWidth=   fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getSensorWidth(chnNum);
		int sensorHeight=  fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getSensorHeight(chnNum);
		int decimation=fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getDecimateMasks(chnNum);
		double [][] additionalCorrection=correctionInNodes.mapToPixels(
				decimation,
				sensorWidth,
				sensorHeight,
				debugLevel);
		if (showIntermediate){
			String [] dbgTitles={"dPX","dPY","alpha","R","G","B"};
			ShowDoubleFloatArrays.showArrays(
					additionalCorrection,
					sensorWidth/decimation,
					sensorHeight/decimation,
					true,
					"AC-"+imgNum,
					dbgTitles);
		}

		return additionalCorrection;
	}


	/**
	 * @param imgNum  number of image
	 * @param showIntermediate - show intermediate images
	 * @param debugLevel debug level
	 * @return CorrectionInNodes data correction, image and grid data for some target grid nodes
	 */

	public CorrectionInNodes extractNodeCorrections(
			int imgNum, // image number
			boolean showIntermediate,
			int debugLevel
			){
//		int debugThreshold=2;
    	int imgRGBIndex=   3;
		int chnNum=fittingStrategy.distortionCalibrationData.gIP[imgNum].channel; // number of sub-camera
		int station=fittingStrategy.distortionCalibrationData.gIP[imgNum].getStationNumber(); // number of sub-camera
		LensDistortionParameters lensDistortionParameters= setupLensDistortionParameters(
				imgNum,
				debugLevel);     // Axial - may be Double.NaN

		int [][] imgUV=	  fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV;
		double [][] imgXY=fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsXY; // for each image, each grid node - a set of of {px,py,contrast,vignR,vignG,vignB} vign* is in the 0..1.0 range
		if ((imgUV==null) || (imgUV.length==0)) {
			System.out.println("expandMeasuredGrid("+imgNum+",..) empty image");
			return null;
		}
		int minU=imgUV[0][0];
		int minV=imgUV[0][1];
		int maxU=minU;
		int maxV=minV;
		for (int i=1;i<imgUV.length;i++){
			if (minU>imgUV[i][0]) minU=imgUV[i][0];
			if (minV>imgUV[i][1]) minV=imgUV[i][1];
			if (maxU<imgUV[i][0]) maxU=imgUV[i][0];
			if (maxV<imgUV[i][1]) maxV=imgUV[i][1];
		}
		int extraMargins=1;
		int [] uv0= {minU-extraMargins,minV-extraMargins}; // target U,V at the stageXYA[0]
		int width= maxU-minU+1+2*extraMargins;
		int height=maxV-minV+1+2*extraMargins;
		double [][] stagePXY= new double [width*height][]; //reprojected {px,py}
		double [][] stageDiffPXY=   new double [width*height][]; // difference between corrected measured and reprojected (to add to correction)
		double [][] stageDiffRGB=   new double [width*height][]; // difference (measured RGB)/(grid RGB) and current correction RGB (pixel sensitivity RGB)
		double [] stageMask=        new double [width*height];   // weight
		for (int i=0;i<stagePXY.length;i++) {
			stagePXY[i]=null;
			stageDiffPXY[i]  =null;
			stageDiffRGB[i] = null;
		}
//		int vignRIndex=3; //  in measured data
		int corrRIndex=3; // in correction vector
//		int gridRIndex=3; // in reprojected vector
		double [] diff=calcYminusFx(this.currentfX, 2*this.imageStartIndex[imgNum],2*this.imageStartIndex[imgNum+1]);
		double [][] photometrics=patternParameters.getPhotometricBySensor(station,chnNum); // head/bottom grid intensity/alpha
		int targetGridWidth=getGridWidth();
		double [][] debugRGB=null;
		if (showIntermediate){
			debugRGB = new double [12][width*height];
			for (int n=0;n<debugRGB.length;n++) for (int i=0;i<debugRGB[n].length;i++) debugRGB[n][i]=0.0;
		}
		for (int i=0;i<imgUV.length;i++){
			int index=width*(imgUV[i][1]-uv0[1]) + (imgUV[i][0]-uv0[0]);
			int targetGridIndex=targetGridWidth*(imgUV[i][1]+patternParameters.V0) +(imgUV[i][0]+patternParameters.U0); // index in photometrics[][]
			int doublePairIndex=2*(this.imageStartIndex[imgNum]+i); // number of a pair in a full vector
			stageMask[index]=this.weightFunction[doublePairIndex];
			stagePXY[index]=null;
			double [] debugCorrVector=null;
			if (showIntermediate) {
				debugCorrVector=interpolateCorrectionVector ( //  vector of {corrX, corrY, alpha, flatfield_red, flatfield_green, flatfield_blue}
						chnNum,
						imgXY[i][0], //double px, measured
						imgXY[i][1]); //double py, measured);
			}
			double [] reprojectedNode= reprojectGridNode( //{pX,pY,grid mask (binary), grid R, grid G, grid B, alpha}
					lensDistortionParameters,
					imgNum,
					imgUV[i][0], //int u, // grid signed u,v
					imgUV[i][1]); //int v
			if (reprojectedNode==null) {
				continue; // out of grid - should not happen here (now - also: target point behind the camera sensor)?
			}
//			double [] reprojPXY={reprojectedNode[0],reprojectedNode[1]};
			double [] nodePXY={this.Y[doublePairIndex],this.Y[doublePairIndex+1]};
			stagePXY[index]=nodePXY;// measured pixels Px,Py with correction applied  // reprojPXY;
//			double [] diffPXY= {imgXY[i][0]-corrVector[0]-reprojectedNode[0],imgXY[i][1]-debugCorrVector[1]-reprojectedNode[1]};
			double [] diffPXY= {diff[2*i],diff[2*i+1]};
			stageDiffPXY[index]=diffPXY;
			//{px,py,contrast,vignR,vignG,vignB}
			double [] diffRGB={0.0,0.0,0.0};
			for (int c=0;c<diffRGB.length;c++){
				double gridPhotometrics=photometrics[c][targetGridIndex];
//				if (gridPhotometrics>0.0) diffRGB[c]=imgXY[i][imgRGBIndex+c]/gridPhotometrics-debugCorrVector[corrRIndex+c];
				if (gridPhotometrics>0.0) diffRGB[c]=imgXY[i][imgRGBIndex+c]/gridPhotometrics; // don't use old correction at all!
			}
			stageDiffRGB[index]=diffRGB;
			stageMask[index]=this.weightFunction[2*(this.imageStartIndex[imgNum]+i)];

			if (showIntermediate) for (int c=0;c<3;c++){
				debugRGB[4*c+0][index]=photometrics[c][targetGridIndex];
				debugRGB[4*c+1][index]=imgXY[i][imgRGBIndex+c];
				debugRGB[4*c+2][index]=debugCorrVector[corrRIndex+c];
				debugRGB[4*c+3][index]=imgXY[i][imgRGBIndex+c]/photometrics[c][targetGridIndex];
			}

		}
		if (showIntermediate){
			double [][] debugData = new double [8][width*height];
			String [] dbgTitles={"rep-X","rep-Y","dX","dY","R","G","B","Weight"};//
			for (int i=0;i<debugData[0].length;i++){
				if (stagePXY[i]==null){
					for (int j=0;j<debugData.length;j++) {
						debugData[j][i]=Double.NaN; // 0.0?
					}
				} else {
					debugData[0][i]=stagePXY[i][0];
					debugData[1][i]=stagePXY[i][1];
					debugData[2][i]=  stageDiffPXY[i][0];
					debugData[3][i]=  stageDiffPXY[i][1];
					debugData[4][i]=  stageDiffRGB[i][0];
					debugData[5][i]=  stageDiffRGB[i][1];
					debugData[6][i]=  stageDiffRGB[i][2];
					debugData[7][i]=     stageMask[i];
				}

			}
			ShowDoubleFloatArrays.showArrays(
					debugData,
					width,
					height,
					true,
					"PRE_EXP-"+imgNum+"-"+chnNum,
					dbgTitles);
			String [] dbgTitles1={"R-tar","R-grid","R-corr","R-FF","G-tar","G-grid","G-corr","G-FF","B-tar","B-grid","B-corr","B-FF",};//
			ShowDoubleFloatArrays.showArrays(
					debugRGB,
					width,
					height,
					true,
					"CORR-RGB-"+imgNum+"-"+chnNum,
					dbgTitles1);

		}
		CorrectionInNodes correctionInNodes=new CorrectionInNodes(
				imgNum,
				uv0[0],
				uv0[1],
				width,
				height,
				stagePXY,
				stageDiffPXY,
				stageDiffRGB,
				stageMask
				);
		return correctionInNodes;
	}

	class CorrectionInNodes{
		public int numImg;
		public Rectangle uv0;
		public double [][] reprojPXY; //= new double [width*height][]; //reprojected {px,py}
		public double [][] diffPXY; //=   new double [width*height][]; // difference between corrected measured and reprojected (to add to correction)
		public double [][] diffRGB; //=   new double [width*height][]; // difference (measured RGB)/(grid RGB) and current correction RGB (pixel sensitivity RGB)
		public double []   mask;
//		public int stageMasksSensor=0, stageMasksTarget=1, stageMasksContrast=2;
		public CorrectionInNodes (
				int numImg,
				int u0,
				int v0,
				int width,
				int height,
				double [][] reprojPXY, //= new double [width*height][]; //reprojected {px,py}
				double [][] diffPXY, //=   new double [width*height][]; // difference between corrected measured and reprojected (to add to correction)
				double [][] diffRGB, //=   new double [width*height][]; // difference (measured RGB)/(grid RGB) and current correction RGB (pixel sensitivity RGB)
				double [] mask //=     new double [width*height][]; // {sensor mask, target mask, measured contrast}
		){
			this.numImg=numImg;
			this.uv0=new Rectangle(u0,v0,width,height);
			this.reprojPXY=reprojPXY; //= new double [width*height][]; //reprojected {px,py}
			this.diffPXY=diffPXY; //=   new double [width*height][]; // difference between corrected measured and reprojected (to add to correction)
			this.diffRGB=diffRGB; //=   new double [width*height][]; // difference (measured RGB)/(grid RGB) and current correction RGB (pixel sensitivity RGB)
			this.mask=mask; //=     new double [width*height][]; // {sensor mask, target mask, measured contrast}
		}


		public void show(
				String title
				){
				double [][] debugData = new double [8][this.uv0.width*this.uv0.height];
				String [] dbgTitles={"rep-X","rep-Y","dX","dY","R","G","B","Weight"};
				for (int i=0;i<debugData[0].length;i++){
					if (this.reprojPXY[i]==null){
						for (int j=0;j<debugData.length;j++) {
							debugData[j][i]=Double.NaN; // 0.0?
						}
					} else {
						debugData[0][i]=this.reprojPXY[i][0];
						debugData[1][i]=this.reprojPXY[i][1];
						debugData[2][i]=  this.diffPXY[i][0];
						debugData[3][i]=  this.diffPXY[i][1];
						debugData[4][i]=  this.diffRGB[i][0];
						debugData[5][i]=  this.diffRGB[i][1];
						debugData[6][i]=  this.diffRGB[i][2];
						debugData[7][i]=     this.mask[i];
					}
				}
				ShowDoubleFloatArrays.showArrays(
						debugData,
						this.uv0.width,
						this.uv0.height,
						true,
						title,
						dbgTitles);
		}
		/**
		 * Convert correction for grid nodes (detected and extrapolated) into decimated pixel array
		 * result should be added to the current (prior) correction. Use alpha as weight when accumulating for multiple images
		 * @param decimation decimate correction pixels from sensor pixels
		 * @param sensorWidth sensor width in pixels (2592)
		 * @param sensorHeight sensor height in pixels (1936)
		 * @param debugLevel debug level (verbose if >3)
		 * @return scan-line pixels correction arrays {dpX,dpY,alpha,R,G,B}[pixelIndex]
		 */
		public double [][] mapToPixels(
				int decimation,
				int sensorWidth,
				int sensorHeight,
				int debugLevel){
			int debugThreshold=2;
			int sWidth= (sensorWidth-1)/decimation+1; // decimated width (648)
			int sHeight=(sensorHeight-1)/decimation+1; // decimated height (484)

			int [] uvInc={0,1,this.uv0.width,this.uv0.width+1}; // four corners as vu index
			int [][] cycles={ // counter-clockwise corners bounding the area  (only orthogonal sides?)
					{1,0,2},
					{2,3,1},
					{0,2,3},
					{3,1,0}};

			double [][] thisPCorr=  new double [6][sWidth*sHeight]; // calculate for a single (this) image, accumulate in the end
			int    []   thisCounted=new    int    [sWidth*sHeight]; // some pixels accumulated twice - divide in the end
			for (int n=0;n<thisPCorr.length;n++) for (int i=0;i<thisPCorr[0].length;i++) thisPCorr[n][i]=0.0;
			for (int i=0;i<thisCounted.length;i++) thisCounted[i]=0;

			// now use imgData array to fill thisPCorr by linear interpolation
			for (int v=0;v<(this.uv0.height-1); v++) for (int u=0; u<(this.uv0.width-1);u++){
				int vu=u+this.uv0.width*v;
                double [][] cornerXY =new double[4][];
                for (int i=0;i<uvInc.length;i++){
                	int vu1=vu+uvInc[i];
                	cornerXY[i]=null;
                	if (this.reprojPXY[vu1]!=null){
                		double w=this.mask[vu1];
                		if (w>0.0) {
                			cornerXY[i]=new double[2];
                			cornerXY[i][0]=this.reprojPXY[vu1][0];
                			cornerXY[i][1]=this.reprojPXY[vu1][1];
                		}
                	}
                }
                boolean [] cycleFits=new boolean[cycles.length];
                boolean anyFits=false;
                for (int i=0;i<cycles.length;i++){
                	cycleFits[i]=true;
                	for (int j=0;j<cycles[i].length;j++) if (cornerXY[cycles[i][j]]==null) {
                		cycleFits[i]=false;
                		break;
                	}
                	anyFits |=cycleFits[i];
                }
                if (!anyFits) continue; // not a single cycle
				if (debugLevel>debugThreshold) {
					String debugString="cycleFits ";
					for (int i =0;i<cycleFits.length; i++) debugString+=" "+cycleFits[i];
					System.out.println(debugString);
				}
                if (cycleFits[0]&&cycleFits[1]){ // remove overlaps
                	cycleFits[2]=false;
                	cycleFits[3]=false;
                }
                boolean minMaxUndefined=true;
				double minX=0,maxX=0,minY=0,maxY=0;
				// find bounding rectangle;
				for (int nCycle=0;nCycle<cycles.length;nCycle++) if (cycleFits[nCycle]){
					int [] cycle=cycles[nCycle];
					for (int corner=0; corner<cycle.length;corner++){
						if (minMaxUndefined || (minX>cornerXY[cycle[corner]][0])) minX=cornerXY[cycle[corner]][0];
						if (minMaxUndefined || (maxX<cornerXY[cycle[corner]][0])) maxX=cornerXY[cycle[corner]][0];
						if (minMaxUndefined || (minY>cornerXY[cycle[corner]][1])) minY=cornerXY[cycle[corner]][1];
						if (minMaxUndefined || (maxY<cornerXY[cycle[corner]][1])) maxY=cornerXY[cycle[corner]][1];
						minMaxUndefined=false;
					}
				}
				int iMinX=(int) Math.floor(minX/decimation);
				int iMinY=(int) Math.floor(minY/decimation);
				int iMaxX=(int) Math.ceil(maxX/decimation);
				int iMaxY=(int) Math.ceil(maxY/decimation);
				// not sure if these checks are needed, got out of bounds wheriDy was =484=sHeight
				if (iMinX<0) iMinX=0;
				if (iMaxX>=sWidth) iMaxX=sWidth-1;
				if (iMinY<0) iMinY=0;
				if (iMaxY>=sHeight) iMaxY=sHeight-1;
				double [] originXY=new double [2];
				double [] endXY=new double [2];
				boolean debugHadPixels=false;
//TODO: scan X,Y in this rectangle, for points in defined squares/triangles find if the point is inside (accurate not to loose any).
				for (int idY=iMinY; idY<=iMaxY;idY++){

					double pY=idY*decimation; // in sensor pixels
					for (int idX=iMinX; idX<=iMaxX;idX++){
						double pX=idX*decimation; // in sensor pixels
						// scan allowed triangles, usually 2
						for (int nCycle=0;nCycle<cycles.length;nCycle++) if (cycleFits[nCycle]){
							int [] cycle=cycles[nCycle];
							// is this point inside?
							boolean inside=true;
							for (int nEdge=0;nEdge<cycle.length;nEdge++){
								int nextNEdge=(nEdge==(cycle.length-1))?0:(nEdge+1);

								originXY[0]=this.reprojPXY[vu+uvInc[cycle[nEdge]]][0];     // imgData[2][vu+uvInc[cycle[nEdge]]];
								originXY[1]=this.reprojPXY[vu+uvInc[cycle[nEdge]]][1];     // imgData[3][vu+uvInc[cycle[nEdge]]];
								endXY[0]=   this.reprojPXY[vu+uvInc[cycle[nextNEdge]]][0]; // imgData[2][vu+uvInc[cycle[nextNEdge]]];
								endXY[1]=   this.reprojPXY[vu+uvInc[cycle[nextNEdge]]][1]; // imgData[3][vu+uvInc[cycle[nextNEdge]]];
								if (((pX-originXY[0])*(endXY[1]-originXY[1]) - (pY-originXY[1])*(endXY[0]-originXY[0]))<0.0){
									inside=false;
									break;
								}
							}
							if (!inside) continue; // point is outside of the interpolation area, try next triangle (if any)
							if (debugLevel>debugThreshold) {
								System.out.println("idX="+idX+" idY="+idY+" nCycle="+nCycle);
								String debugString1="cycle:";
								for (int i =0;i<cycle.length; i++) debugString1+=" "+cycle[i];
								System.out.println(debugString1);
							}

							/* interpolate:
							1. taking cycles[0] as origin and two (non co-linear) edge vectors - V1:from 0 to 1 and V2 from 1 to 2
							    find a1 and a2  so that vector V  (from 0  to pXY) = a1*V1+ a2*V2
							2. if F0 is the value of the interpolated function at cycles[0], F1 and F2 - at cycles[1] and cycles2
							   then F=F0+(F1-F0)*a1 +(F2-F1)*a2
							 */
							double [] XY0={this.reprojPXY[vu+uvInc[cycle[0]]][0],this.reprojPXY[vu+uvInc[cycle[0]]][1]};
							double [] XY1={this.reprojPXY[vu+uvInc[cycle[1]]][0],this.reprojPXY[vu+uvInc[cycle[1]]][1]};
							double [] XY2={this.reprojPXY[vu+uvInc[cycle[2]]][0],this.reprojPXY[vu+uvInc[cycle[2]]][1]};
							double [] V= {pX-XY0[0],pY-XY0[1]};
							double [][] M={
									{XY1[0]-XY0[0],XY2[0]-XY1[0]},
									{XY1[1]-XY0[1],XY2[1]-XY1[1]}};
							double det=M[0][0]*M[1][1]-M[1][0]*M[0][1];
							double [][] MInverse={
									{ M[1][1]/det,-M[0][1]/det},
									{-M[1][0]/det, M[0][0]/det}};
							double [] a12={
									MInverse[0][0]*V[0]+MInverse[0][1]*V[1],
									MInverse[1][0]*V[0]+MInverse[1][1]*V[1]};
							int pCorrIndex=idY*sWidth+idX;
// some points may be accumulated multiple times - thisPCorr[3] will take care of this
							if (debugLevel>debugThreshold) {
								System.out.println("XY0="+IJ.d2s(XY0[0],3)+":"+IJ.d2s(XY0[1],3));
								System.out.println("XY1="+IJ.d2s(XY1[0],3)+":"+IJ.d2s(XY1[1],3));
								System.out.println("XY2="+IJ.d2s(XY2[0],3)+":"+IJ.d2s(XY2[1],3));
								System.out.println("M00="+IJ.d2s(M[0][0],3)+" M01="+IJ.d2s(M[0][1],3));
								System.out.println("M10="+IJ.d2s(M[1][0],3)+" M11="+IJ.d2s(M[1][1],3));
								System.out.println("MInverse00="+IJ.d2s(MInverse[0][0],5)+" MInverse01="+IJ.d2s(MInverse[0][1],5));
								System.out.println("MInverse10="+IJ.d2s(MInverse[1][0],5)+" MInverse11="+IJ.d2s(MInverse[1][1],5));
								System.out.println("a12="+IJ.d2s(a12[0],3)+":"+IJ.d2s(a12[1],3));
								System.out.println("this.diffPXY[vu+uvInc[cycle[0]]][0]="+IJ.d2s(this.diffPXY[vu+uvInc[cycle[0]]][0],3)+
										"this.diffPXY[vu+uvInc[cycle[0]]][1]="+IJ.d2s(this.diffPXY[vu+uvInc[cycle[0]]][1],3));
								System.out.println("this.diffPXY[vu+uvInc[cycle[1]]][0]="+IJ.d2s(this.diffPXY[vu+uvInc[cycle[1]]][0],3)+
										"this.diffPXY[vu+uvInc[cycle[1]]][1]="+IJ.d2s(this.diffPXY[vu+uvInc[cycle[1]]][1],3));
								System.out.println("this.diffPXY[vu+uvInc[cycle[2]]][0]="+IJ.d2s(this.diffPXY[vu+uvInc[cycle[2]]][0],3)+
										"this.diffPXY[vu+uvInc[cycle[2]]][1]="+IJ.d2s(this.diffPXY[vu+uvInc[cycle[2]]][1],3));
							}

							double [] corr={
									 this.diffPXY[vu+uvInc[cycle[0]]][0]+ // dPx
									(this.diffPXY[vu+uvInc[cycle[1]]][0]-this.diffPXY[vu+uvInc[cycle[0]]][0])*a12[0]+
									(this.diffPXY[vu+uvInc[cycle[2]]][0]-this.diffPXY[vu+uvInc[cycle[1]]][0])*a12[1],

									 this.diffPXY[vu+uvInc[cycle[0]]][1]+ // dPy
									(this.diffPXY[vu+uvInc[cycle[1]]][1]-this.diffPXY[vu+uvInc[cycle[0]]][1])*a12[0]+
									(this.diffPXY[vu+uvInc[cycle[2]]][1]-this.diffPXY[vu+uvInc[cycle[1]]][1])*a12[1],

									 this.mask[vu+uvInc[cycle[0]]]+ // alpha
									(this.mask[vu+uvInc[cycle[1]]]-this.mask[vu+uvInc[cycle[0]]])*a12[0]+
									(this.mask[vu+uvInc[cycle[2]]]-this.mask[vu+uvInc[cycle[1]]])*a12[1],

									 this.diffRGB[vu+uvInc[cycle[0]]][0]+ // Red measured/pattern
									(this.diffRGB[vu+uvInc[cycle[1]]][0]-this.diffRGB[vu+uvInc[cycle[0]]][0])*a12[0]+
									(this.diffRGB[vu+uvInc[cycle[2]]][0]-this.diffRGB[vu+uvInc[cycle[1]]][0])*a12[1],

									 this.diffRGB[vu+uvInc[cycle[0]]][1]+ // Red measured/pattern
									(this.diffRGB[vu+uvInc[cycle[1]]][1]-this.diffRGB[vu+uvInc[cycle[0]]][1])*a12[0]+
									(this.diffRGB[vu+uvInc[cycle[2]]][1]-this.diffRGB[vu+uvInc[cycle[1]]][1])*a12[1],

									 this.diffRGB[vu+uvInc[cycle[0]]][2]+ // Red measured/pattern
									(this.diffRGB[vu+uvInc[cycle[1]]][2]-this.diffRGB[vu+uvInc[cycle[0]]][2])*a12[0]+
									(this.diffRGB[vu+uvInc[cycle[2]]][2]-this.diffRGB[vu+uvInc[cycle[1]]][2])*a12[1]};
							if (debugLevel>debugThreshold) {
								System.out.println("corr="+IJ.d2s(corr[0],3)+" "+IJ.d2s(corr[1],3)+" "+IJ.d2s(corr[2],3));
							}
 if (pCorrIndex>thisPCorr[0].length) {
//	 System.out.println("imgNum=" + imgNum+": "+	fittingStrategy.distortionCalibrationData.gIP[imgNum].path);
	 System.out.println("thisPCorr[0].length="+thisPCorr[0].length+" pCorrIndex="+pCorrIndex+" sWidth="+sWidth+" idY="+idY+" idX="+idX);
 }
                            for (int i=0;i<corr.length;i++) {
                            	thisPCorr[i][pCorrIndex]+= corr[i]; // OOB: -8, -1433
                            }
							thisCounted[pCorrIndex]++;

							if (debugLevel>debugThreshold) {
								debugHadPixels=true;
							}
						}
					} // idX
					// use same order in calculations, make sure no gaps
				} // idY
				if ((debugLevel>debugThreshold) && (debugHadPixels)){
//					if (!debugExit) {
						System.out.println(
								" minX="+IJ.d2s(minX,1)+
								" maxX="+IJ.d2s(maxX,1));
						System.out.println(
								" minY="+IJ.d2s(minY,1)+
								" maxY="+IJ.d2s(maxY,1));
						System.out.println(
								" iMinX="+iMinX+
								" iMaxX="+iMaxX);
						System.out.println(
								" iMinY="+iMinY+
								" iMaxY="+iMaxY);
//					}
//					if (!debugExit) debugCntr--;
//					if (debugCntr==0) debugExit=true; // exit after first non-empty tile

				}
			} //for (int v=0;v<(this.uv0.height-1); v++) for (int u=0; u<(this.uv0.width-1);u++){
            for (int i=0;i<thisCounted.length;i++) if (thisCounted[i]>1) {
            	for (int j=0;j<thisPCorr[i].length;j++)	thisPCorr[j][i]/= thisCounted[i];
            }
            return thisPCorr;
		}
	}

	class DirInc{
		private int top=   1 | 2 | 4 | 8 | 16;
		private int bottom=1 |             16 | 32 | 64 | 128;
		private int left=  1 | 2 |                   64 | 128;
		private int right=         4 | 8 | 16 | 32 | 64;
		private int [] inc=null;
		private int [] validDirs=null;
		private double [][] unityVector=null;
		public int dirs=8;
		public DirInc(int width, int height){
//			int [] dirs8={1,1+width,width,-1+width,-1,-1-width,-width,1-width};
			int [][] incXY8={{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}};
			this.inc=new int [incXY8.length];
			this.unityVector=new double [incXY8.length][2];
			for (int i=0;i<incXY8.length;i++){
				this.inc[i]=incXY8[i][0]+width*incXY8[i][1];
				double len=Math.sqrt(incXY8[i][0]*incXY8[i][0]+incXY8[i][1]*incXY8[i][1]);
				this.unityVector[i][0]=incXY8[i][0]/len;
				this.unityVector[i][1]=incXY8[i][1]/len;
			}
//			this.inc=dirs8;
			this.validDirs=new int [width*height];
			for (int i=0;i<this.validDirs.length;i++) this.validDirs[i]=0xff;
			for (int i=0;i<width;i++){
				this.validDirs[                 i]&=top;
				this.validDirs[(height-1)*width+i]&=bottom;
			}
			for (int i=0;i<height;i++){
				this.validDirs[i*width]&=left;
				this.validDirs[i*width + width-1]&=right;
			}

		}
		public int newIndex(int oldIndex, int dir){
			if ((validDirs[oldIndex] & (1<<dir))==0) return -1; // invalid dir for this location (border)
			return oldIndex+this.inc[dir];
		}
		public double [] unity(int dir) {
			return this.unityVector[(dir+this.unityVector.length)%this.unityVector.length];
		}
	}

	public class PixXYUV{
		double [][]xy=null;
		int [][]uv=null;
		double [] alpha=null;
		double [][]dxy=null;
		public PixXYUV(){}
		public PixXYUV(int len){
			this.uv=new int [len][2];
			this.xy=new double [len][2];
			this.alpha=new double [len];
			this.dxy=new double [len][2];
		}
	}




	/**
	 * Interpolate (bi-linear) X/Y corrections and flat-field data for the sensor
	 * @param chnNum - sensor (channel) number
	 * @param px     - pixel X coordinate (non-decimated)
	 * @param py     - pixel Y coordinate (non-decimated)
	 * @return       - vector of {corrX, corrY, alpha, flatfield_red, flatfield_green, flatfield_blue}
	 */
	public double [] interpolateCorrectionVector (
			int chnNum,
			double px,
			double py){
		if (this.pixelCorrection==null){
			double [] vector={0.0,0.0,1.0,1.0,1.0,1.0};
			return vector;
		}
//		this.pixelCorrectionDecimation=fittingStrategy.distortionCalibrationData.eyesisCameraParameters.decimateMasks;
//		this.pixelCorrectionWidth=   fittingStrategy.distortionCalibrationData.eyesisCameraParameters.sensorWidth;
//		this.pixelCorrectionHeight=  fittingStrategy.distortionCalibrationData.eyesisCameraParameters.sensorHeight;

		int sensorCorrWidth= getSensorCorrWidth(chnNum);
		int sensorCorrHeight=this.pixelCorrection[chnNum][0].length/sensorCorrWidth;

		int [] ix={(int) Math.floor(px/getDecimateMasks(chnNum)), (int) Math.floor(px/getDecimateMasks(chnNum))+1};
		int [] iy={(int) Math.floor(py/getDecimateMasks(chnNum)),(int) Math.floor(py/getDecimateMasks(chnNum))+1};
		for (int i=0;i<2;i++){
			if (ix[i]<0) ix[i]=0;
			else if (ix[i]>=sensorCorrWidth) ix[i]=sensorCorrWidth-1;
			if (iy[i]<0) iy[i]=0;
			else if (iy[i]>=sensorCorrHeight) iy[i]=sensorCorrHeight-1;
		}
		int index00=ix[0] + iy[0]*sensorCorrWidth;
		int indexX0=ix[1] + iy[0]*sensorCorrWidth;
		int index0Y=ix[0] + iy[1]*sensorCorrWidth;
		int indexXY=ix[1] + iy[1]*sensorCorrWidth;

		double corrDX=0,corrDY=0;
		if ((px>ix[0])&& (px<ix[1])) corrDX=px-ix[0];
		if ((py>iy[0])&& (py<iy[1])) corrDY=py-iy[0];
		double [] vector=new double [this.pixelCorrection[chnNum].length];
		for (int n=0;n<vector.length;n++){
			// bilinear interpolation
			vector[n]=
				(1-corrDX)* (1-corrDY)* this.pixelCorrection[chnNum][n][index00]+
				corrDX * (1-corrDY)* this.pixelCorrection[chnNum][n][indexX0]+
				(1-corrDX)*    corrDY * this.pixelCorrection[chnNum][n][index0Y]+
				corrDX *    corrDY * this.pixelCorrection[chnNum][n][indexXY];
		}
		return vector;
	}
	/**
	 * Bilinear interpolate sensor mask array
	 * @param mask decimated mask data
	 * @param px     - pixel X coordinate (non-decimated)
	 * @param py     - pixel Y coordinate (non-decimated)
	 * @return interpolated mask data at specified fractional pixel
	 */
	public double interpolateMask (
			int       chnNum,
			double [] mask,
			double px,
			double py){
///		this.pixelCorrectionDecimation=fittingStrategy.distortionCalibrationData.eyesisCameraParameters.decimateMasks;
///		this.pixelCorrectionWidth=   fittingStrategy.distortionCalibrationData.eyesisCameraParameters.sensorWidth;
///		this.pixelCorrectionHeight=  fittingStrategy.distortionCalibrationData.eyesisCameraParameters.sensorHeight;

		int sensorCorrWidth= getSensorCorrWidth(chnNum); // (this.pixelCorrectionWidth-1)/this.pixelCorrectionDecimation+1;
		int sensorCorrHeight=mask.length/sensorCorrWidth;
		int [] ix={(int) Math.floor(px/getDecimateMasks(chnNum)), (int) Math.floor(px/getDecimateMasks(chnNum))+1};
		int [] iy={(int) Math.floor(py/getDecimateMasks(chnNum)), (int) Math.floor(py/getDecimateMasks(chnNum))+1};
		for (int i=0;i<2;i++){
			if (ix[i]<0) ix[i]=0;
			else if (ix[i]>=sensorCorrWidth) ix[i]=sensorCorrWidth-1;
			if (iy[i]<0) iy[i]=0;
			else if (iy[i]>=sensorCorrHeight) iy[i]=sensorCorrHeight-1;
		}
		int index00=ix[0] + iy[0]*sensorCorrWidth;
		int indexX0=ix[1] + iy[0]*sensorCorrWidth;
		int index0Y=ix[0] + iy[1]*sensorCorrWidth;
		int indexXY=ix[1] + iy[1]*sensorCorrWidth;

		double corrDX=0,corrDY=0;
		if ((px>ix[0])&& (px<ix[1])) corrDX=px-ix[0];
		if ((py>iy[0])&& (py<iy[1])) corrDY=py-iy[0];
		double result=
				(1-corrDX)* (1-corrDY)* mask[index00]+
				corrDX * (1-corrDY)* mask[indexX0]+
				(1-corrDX)*    corrDY * mask[index0Y]+
				corrDX *    corrDY * mask[indexXY];
		return result;
	}


/**
 *   after fitting finished and accepted - 	fittingStrategy.saveSeriesVector(double [] vector)
 */
	public void saveFittingSeries() {
		fittingStrategy.saveSeriesVector(this.currentVector);
	}
	/*
	 * For each image in the series:

    	public double [] fittingStrategy.getImageParametersVector(int numImg, double [] parameterVector);
    	 * Calculates current values of all parameters for the particular sensor - some ("fixed")
    	 * are taken from the data stored for this individual image, others - from the parameter
    	 * vector (used in fitting)
    	 * @param numImg number of image
    	 * @param vector parameters vector
    	 * @return vector used for the current image (parameters influencing the acquired grid
    	 * on the sensor (common parameters and those of the sensor's subchannel)

   public void calcInterParamers(
    		double [] parVect,
    		boolean [] mask, // calculate only selected derivatives (all parVect values are still
    		boolean calculateDerivatives // calculate this.interParameterDerivatives -derivatives array (false - just this.values)
    		){
     * Calculate/set  this.lensDistortionParameters and this.interParameterDerivatives
     * @param parVect 21-element vector for eyesis sub-camera, including common and individual parameters
     * @param mask -mask - which partial derivatives are needed to be calculated (others will be null)
     * @param calculateDerivatives calculate array of partial derivatives, if false - just the values


For each point in the image
      public double [][] lensDistortionParameters.reorderPartialDerivatives (double [][] srcDerivatives){
      double [][] lensDistortionParameters.calcPartialDerivatives(
        		double xp, // target point horizontal, positive - right,  mm
        		double yp, // target point vertical,   positive - down,  mm
        		double zp, // target point horizontal, positive - away from camera,  mm
        		boolean calculateAll){ // calculate derivatives, false - values only

    public double [][] interParameterDerivatives=null; //partial derivative matrix from subcamera-camera-goniometer to single camera (12x21)
    public double []   currentVector; // current variable parameter vector
    public double []   Y=null; // array of "y" - for each grid image, each defined grid node - 2 elements
    public double [][] targetXYZ=null; // array of target {x,y,z} matching each image each grid point
    public double []   fX=null; // array of "f(x)" - simulated data for all images, combining pixel-X and pixel-Y (odd/even)
    public double [][] jacobian=null; // partial derivatives of fX (above) by parameters to be adjusted (rows)
	 */

	public ImagePlus simulatePatternOnSensor(
			int stationNumber,
			int subCam,
			double goniometerTilt,
			double goniometerAxial,
			double goniometerInterAxis,
			SimulationPattern.SimulParameters simulParametersDefault,
			int threadsMax,
			boolean updateStatus,
			int mspDebugLevel,
			int global_debug_level, // DEBUG_LEVEL
			int debug_level // debug level used inside loops
	){
		MatchSimulatedPattern matchSimulatedPattern = new MatchSimulatedPattern(64); // new instance, all reset, FFTSize=64 will not be used
		matchSimulatedPattern.debugLevel = mspDebugLevel;
		//		MatchSimulatedPattern.DistortionParameters distortionParameters = modifyDistortionParameters();
		//		SimulationPattern.SimulParameters simulParameters = modifySimulParameters();
		int sensorWidth=  getSensorWidth(subCam);
		int sensorHeight= getSensorHeight(subCam);

		double [][][] hintGrid=estimateGridOnSensor(
				stationNumber,
				subCam,
				goniometerTilt, // Tilt, goniometerHorizontal
				goniometerAxial,  // Axial,goniometerAxial
				goniometerInterAxis, // inter-axis angle
				-1, // use camera parameters, not imageSet
				true // filter border
		);
		if (hintGrid==null){
			String msg="Grid is not visible for subcamera="+subCam+",  tilt="+goniometerTilt+", axial="+goniometerAxial;
			IJ.showMessage("Error",msg);
			System.out.println("Error: "+msg);
			return null;
		}
		if (global_debug_level>1){
			double [][] pixels=new double[4][hintGrid.length*hintGrid[0].length];
			int index=0;
			String [] titles={"pixel-X","pixel-Y","grid-U","grid-V"};
			for (int v=0; v<hintGrid.length;v++) for (int u=0;u<hintGrid[v].length;u++){
				if (hintGrid[v][u]!=null){
					for (int i=0; i<4;i++)	pixels[i][index]=hintGrid[v][u][i];
				} else {
					for (int i=0; i<4;i++)	pixels[i][index]=-1;
				}
				index++;
			}
			ShowDoubleFloatArrays.showArrays(pixels, hintGrid[0].length, hintGrid.length,  true, "hintGrid", titles);

		}

		if (global_debug_level>0){
			System.out.println("simulatePatternOnSensor(): subcamera="+subCam+",  tilt="+goniometerTilt+", axial="+goniometerAxial);
		}
		int numCells=matchSimulatedPattern.restoreSimulatedPatternGridFromHint(hintGrid, sensorWidth, sensorHeight);
		matchSimulatedPattern.recalculateWaveVectors (
				   updateStatus,
				   debug_level);// debug level used inside loops

		if (global_debug_level>0){
			System.out.println("simulatePatternOnSensor(): "+numCells+" grid cells");
		}
		SimulationPattern.SimulParameters simulParameters = simulParametersDefault.clone();
		SimulationPattern simulationPattern=new SimulationPattern(simulParameters);
		double [][] xy0={{simulParameters.offsetX,simulParameters.offsetY},{simulParameters.offsetX-0.5,simulParameters.offsetY-0.5}} ;
// TODO: add marks for the laser pointers when visible?
		float[] simPixels=simulationPattern.simulateGrid (
				matchSimulatedPattern.getDArray(),
				2, // gridFrac, // number of grid steps per pattern full period
				simulParameters,
				matchSimulatedPattern.getWOI(),
				1, // simulParameters.subdiv/2,
				xy0[0],    // add to patternGrid xy
				threadsMax,
				updateStatus,
				(debug_level>1)?1:0); //debug_level); // debug level
		if (global_debug_level>0){
			System.out.println("simulatePatternOnSensor(): simPixels.length="+simPixels.length+" sensorWidth="+sensorWidth+" sensorHeight="+sensorHeight);
		}
		for (int i=0;i<simPixels.length;i++) simPixels[i]*=255.0;
		ImageProcessor ip_simGrid = new FloatProcessor(sensorWidth, sensorHeight);
		ip_simGrid.setPixels(simPixels);
		ip_simGrid.resetMinAndMax();
		ImagePlus imp_simGrid= new ImagePlus("Simulated_Grid_CHN"+subCam+"_TILT"+goniometerTilt+"_AXIAL"+goniometerAxial, ip_simGrid);
		return imp_simGrid;
	}

//TODO: add additional parameter - process all, but with matched pointers less than 2
	public int applyHintedGrids(
			LaserPointer laserPointer, // LaserPointer object that specifies actual laser pointers on the target
			boolean removeOutOfGridPointers,
			double  hintGridTolerance, // allowed mismatch (fraction of period) or 0 - orientation only
			boolean processAll, // if true - process all images, false - only disabled
			boolean ignoreLaserPointers, // ignore laser pointers, rely on hints only
			boolean processBlind, // try to match without known orientation and no laser pointers
			int     imageNumber, // <0 - all, >=0 only this image
			boolean useSetData,
			int threadsMax,
			boolean updateStatus,
			int mspDebugLevel,
			int global_debug_level, // DEBUG_LEVEL
			int debug_level // debug level used inside loops
	){
		return applyHintedGrids(
				laserPointer, // LaserPointer object that specifies actual laser pointers on the target
				removeOutOfGridPointers,
				hintGridTolerance, // allowed mismatch (fraction of period) or 0 - orientation only
				processAll, // if true - process all images, false - only disabled
				ignoreLaserPointers, // ignore laser pointers, rely on hints only
				processBlind, // try to match without known orientation and no laser pointers
				imageNumber, // <0 - all, >=0 only this image
				0, // int     start_set,
				this.fittingStrategy.distortionCalibrationData.getNumSets()-1, // int     end_set,
				useSetData,
				threadsMax,
				updateStatus,
				mspDebugLevel,
				global_debug_level, // DEBUG_LEVEL
				debug_level // debug level used inside loops
		);
	}


	public int applyHintedGrids(
			LaserPointer laserPointer, // LaserPointer object that specifies actual laser pointers on the target
			boolean removeOutOfGridPointers,
			double  hintGridTolerance, // allowed mismatch (fraction of period) or 0 - orientation only
			boolean processAll, // if true - process all images, false - only disabled
			boolean ignoreLaserPointers, // ignore laser pointers, rely on hints only
			boolean processBlind, // try to match without known orientation and no laser pointers
			int     imageNumber, // <0 - all, >=0 only this image
			int     start_set,
			int     end_set,
			boolean useSetData,
			int threadsMax,
			boolean updateStatus,
			int mspDebugLevel,
			int global_debug_level, // DEBUG_LEVEL
			int debug_level // debug level used inside loops
	){
		boolean invert = false;
		int debugThreshold0=0;
		int debugThreshold=2;
		MatchSimulatedPattern matchSimulatedPattern = new MatchSimulatedPattern(64); // new instance, all reset, FFTSize=64 will not be used
		// next 2 lines are not needed for the new instance, but can be
		// used alternatively if keeping it
		//		matchSimulatedPattern.invalidateFlatFieldForGrid(); // Reset Flat Filed calibration - different image.
		//		matchSimulatedPattern.invalidateFocusMask();
		matchSimulatedPattern.debugLevel = mspDebugLevel;
		//		ImagePlus imp_eq = matchSimulatedPattern.applyFlatField(images[nImg]); // current image with grid flat-field  correction

		//		if (debug_level > 0){
		//			System.out.println("\n   ======= Looking for grid, matching pointers in image " +images[nImg].getTitle()+
		//					", initial number of pointers was "+numPointers);
		//		}
		//matchSimulatedPatterns[numSensor].getChannel(images[numSensor])+" ");
		//		MatchSimulatedPattern.DistortionParameters distortionParameters = modifyDistortionParameters();
		//		SimulationPattern.SimulParameters simulParameters = modifySimulParameters();

		boolean noMessageBoxes=true;
		double [] xy0={0.0,0.0} ; //(old) debug only
		int numSuccess=0;
		DistortionCalibrationData dcd=fittingStrategy.distortionCalibrationData;
		for (int numGridImage=0;numGridImage<dcd.gIP.length;numGridImage++) {
			/*
			if (numGridImage >= 1680)	{
				System.out.println("Processing debug image "+numGridImage);
				System.out.println("Processing debug image "+numGridImage);
			}
			if ((set_number >= start_set) &&
					(set_number <= end_set) &&
					(((imageNumber<0) || ((imageNumber==numGridImage)) && (processAll) ||
					(!dcd.gIP[numGridImage].enabled &&
							((hintGridTolerance>0.0) ||
							 ((dcd.gIP[numGridImage].matchedPointers>0)) && !ignoreLaserPointers))))){ // skip no-pointers if only orientation is hinted
			*/
			int set_number = dcd.gIP[numGridImage].getSetNumber();
			if ((set_number >= start_set) && (set_number <= end_set) && // correct set range
					((imageNumber < 0) || (imageNumber==numGridImage)) && // either all images or selected image
					(processAll || !dcd.gIP[numGridImage].enabled) && // "process all" (including disabled) or this is disabled
					((hintGridTolerance > 0.0) || ((dcd.gIP[numGridImage].matchedPointers>0) && !ignoreLaserPointers)) // hint tolerance is provided, or there are lasers not disabled
					){ // skip no-pointers if only orientation is hinted
				
				if (((dcd.gIP[numGridImage].matchedPointers==0) || ignoreLaserPointers)&&
						(dcd.gIS[dcd.get_gIS_index(numGridImage)].orientationEstimated)) {
					if ( !processBlind) {
						if (this.debugLevel>0) {
							System.out.println("\n**** Orientation is not known exactly for image # "+numGridImage+" - "+dcd.gIP[numGridImage].path+
									", and there are no laser pointer references (processBlind==false) - skipping");
						}
						continue;
					} else {
						if (this.debugLevel>0) {
							System.out.println("\n**** Orientation is not known exactly for image # "+numGridImage+" - "+dcd.gIP[numGridImage].path+
									", and there are no laser pointer references, but processBlind is enabled, proceeding");
						}
					}
				}
				if ((dcd.gIP[numGridImage].matchedPointers > 0) && !ignoreLaserPointers) { // just re-enable with the same shifts (will fail if pointers were just added, but it failed anyway) 
					if (!dcd.gIP[numGridImage].enabled) {
						if (this.debugLevel>0) {
							System.out.println("Re-enabling grid #"+numGridImage+" that has pointer(s) with the previously set UVShiftRot =={0,0,0}");
						}
						dcd.gIP[numGridImage].enabled = true;
						dcd.gIP[numGridImage].newEnabled = true;
					}
					continue;
				}
				
				if (this.debugLevel>debugThreshold0) {
					System.out.println("\n---- applyHintedGrids() image #"+numGridImage+" (imageNumber="+imageNumber+") "+
							" dcd.gIP["+numGridImage+"].pixelsXY.length="+dcd.gIP[numGridImage].pixelsXY.length+
							" dcd.gIP["+numGridImage+"].pixelsXY_extra.length="+dcd.gIP[numGridImage].pixelsXY_extra.length+
							" grid period="+dcd.gIP[numGridImage].getGridPeriod()+
							" enabled="+dcd.gIP[numGridImage].enabled+
							" hintedMatch="+dcd.gIP[numGridImage].hintedMatch
							);
					if (this.debugLevel>(debugThreshold)){
						for (int i=0;i<dcd.gIP[numGridImage].pixelsXY.length;i++){
							System.out.println(i+": dcd.gIP["+numGridImage+"].pixelsXY={"+dcd.gIP[numGridImage].pixelsXY[i][0]+
									","+dcd.gIP[numGridImage].pixelsXY[i][1]+"}"+
									" uv={"+dcd.gIP[numGridImage].pixelsUV[i][0]+
									","+dcd.gIP[numGridImage].pixelsUV[i][1]+"}");
						}
						for (int i=0;i<dcd.gIP[numGridImage].pixelsXY_extra.length;i++){
							System.out.println(i+": dcd.gIP["+numGridImage+"].pixelsXY_extra={"+dcd.gIP[numGridImage].pixelsXY_extra[i][0]+
									","+dcd.gIP[numGridImage].pixelsXY_extra[i][1]+"}"+
									" uv={"+dcd.gIP[numGridImage].pixelsUV_extra[i][0]+
									","+dcd.gIP[numGridImage].pixelsUV_extra[i][1]+"}");
						}
					}
				}

				double [][][] pixelsXYSet={
						dcd.gIP[numGridImage].pixelsXY,
						dcd.gIP[numGridImage].pixelsXY_extra};
				int   [][][] pixelsUVSet={
						dcd.gIP[numGridImage].pixelsUV,
						dcd.gIP[numGridImage].pixelsUV_extra};
				// shifts pixelsUV to have minimal u,v of 0 (stores shift in this.minUV), sets PATTERN_GRID
				matchSimulatedPattern.restorePatternGridFromGridList(
						pixelsXYSet, //double [][] pixelsXY,
						pixelsUVSet, // int [][] pixelsUV,
						dcd.gIP[numGridImage].intensityRange
						); // width and height will be calculated from maximal of pixelsXY
				boolean OK=matchSimulatedPattern.createUV_INDEX( /// **** fails here
						null, //imp, // or null - just to determine WOI (when getWOI matches image size)
						xy0, // add to patterGrid xy, null OK
						threadsMax,
						updateStatus,
						global_debug_level, // DEBUG_LEVEL
						debug_level); // debug level used inside loops
				if (!OK) {
					System.out.println("++++++ BUG: in applyHintedGrids() failed in createUV_INDEX()");
					continue;
				}
				double [] goniometerTiltAxial=dcd.getImagesetTiltAxial(numGridImage);
				if ((goniometerTiltAxial==null) || Double.isNaN(goniometerTiltAxial[0])  || Double.isNaN(goniometerTiltAxial[1])){
					if (this.debugLevel>0) {
						System.out.println("No goniometer orientation is available for image # "+numGridImage+" - "+dcd.gIP[numGridImage].path);
					}
				} else {
					if ((numGridImage >= 234) && (numGridImage< 245)) {
						System.out.println("debug numGridImage="+numGridImage);
						System.out.println();
					}
					int station=dcd.getImageStation(numGridImage);
					int setNumber=dcd.gIP[numGridImage].getSetNumber();
					double [][][] hintGrid=estimateGridOnSensor(
							station, // station number
							dcd.gIP[numGridImage].channel,
							goniometerTiltAxial[0], // Tilt, goniometerHorizontal
							goniometerTiltAxial[1],  // Axial,goniometerAxial
							goniometerTiltAxial[2],  // inter-axis angle
							setNumber, // -1 or specific image set
							true // filter border
							);
					if (global_debug_level>0){
						System.out.println("\n**** applyHintedGrids(): processing grid image # "+numGridImage+", path="+dcd.gIP[numGridImage].path);
					}
					if (hintGrid==null){
						if (global_debug_level>0){
							System.out.println("estimateGridOnSensor() failed - skipping");
						}
						dcd.gIP[numGridImage].hintedMatch =0;
						continue;
					}
					int rslt= matchSimulatedPattern.combineGridCalibration(
							laserPointer, // LaserPointer object or null
							ignoreLaserPointers?null:dcd.gIP[numGridImage].laserPixelCoordinates, //pointersXY,
							removeOutOfGridPointers, //
							hintGrid, // predicted grid array (or null)
							hintGridTolerance, // allowed mismatch (fraction of period) or 0 - orientation only
							invert,
							global_debug_level, // DEBUG_LEVEL
							noMessageBoxes );
					if (global_debug_level>0){
						System.out.println("applyHintedGrids(): rslt="+rslt);
					}
					if (rslt<0) { // failed hinting
						dcd.gIP[numGridImage].hintedMatch =0;
					} else {
						// re-create pixelsXY, pixelsXY_extra, pixelsUV, pixelsUV_extra
						int size=0;
						int size_extra=0;
						/*	            		System.out.println("numGridImage="+numGridImage+" matchSimulatedPattern.getHeight()="+matchSimulatedPattern.getHeight()+
	            				" matchSimulatedPattern.getWidth()="+matchSimulatedPattern.getWidth()+
	            				" matchSimulatedPattern.targetUV is "+((matchSimulatedPattern.targetUV==null)?"null":"not null")+
	            				" matchSimulatedPattern.pixelsUV is "+((matchSimulatedPattern.pixelsUV==null)?"null":"not null")
	            				);
	            		System.out.println(
	            				" matchSimulatedPattern.targetUV[0] is "+((matchSimulatedPattern.targetUV[0]==null)?"null":"not null")+
	            				" matchSimulatedPattern.pixelsUV[0] is "+((matchSimulatedPattern.pixelsUV[0]==null)?"null":"not null")
	            				);*/
						for (int v=0;v<matchSimulatedPattern.getHeight();v++) for (int u=0;u<matchSimulatedPattern.getWidth();u++) {
							/*		            		System.out.println("v="+v+", u="+u);
		            		System.out.println(" matchSimulatedPattern.targetUV[v][u] is "+((matchSimulatedPattern.targetUV[v][u]==null)?"null":"not null"));
		            		System.out.println(" matchSimulatedPattern.pixelsUV[v][u] is "+((matchSimulatedPattern.pixelsUV[v][u]==null)?"null":"not null"));*/
							if ((matchSimulatedPattern.targetUV[v][u]!=null) && (matchSimulatedPattern.pXYUV [v][u]!=null)){

								if ((matchSimulatedPattern.targetUV[v][u]!=null) && (matchSimulatedPattern.pXYUV [v][u]!=null) &&
										(matchSimulatedPattern.pXYUV[v][u][0]>=0.0) || (matchSimulatedPattern.pXYUV[v][u][1]>=0.0)) { // disregard negative sensor pixels
									//				            		System.out.println(" matchSimulatedPattern.targetUV[v][u] is "+((matchSimulatedPattern.targetUV[v][u]==null)?"null":"not null"));
									//				            		System.out.println(" matchSimulatedPattern.targetUV[v][u][0]= "+matchSimulatedPattern.targetUV[v][u][0]);
									//				            		System.out.println(" matchSimulatedPattern.targetUV[v][u][1]= "+matchSimulatedPattern.targetUV[v][u][1]); //********
									//				            		System.out.println(" patternParameters is "+((patternParameters==null)?"null":"not null"));
									//				            		int tu=matchSimulatedPattern.targetUV[v][u][0];
									//				            		int tv=matchSimulatedPattern.targetUV[v][u][1];
									//
									if (patternParameters.getXYZM(matchSimulatedPattern.targetUV[v][u][0],matchSimulatedPattern.targetUV[v][u][1],false,station)!=null) {
										size++;
									} else {
										size_extra++;
									}
								}
							}
						}
						// Move to DCD?
						dcd.gIP[numGridImage].resetMask();
						dcd.gIP[numGridImage].pixelsXY=new double [size][6];
						dcd.gIP[numGridImage].pixelsUV=new int    [size][2];
						dcd.gIP[numGridImage].pixelsXY_extra=new double [size_extra][6];
						dcd.gIP[numGridImage].pixelsUV_extra=new int    [size_extra][2];
						int index=0;
						int index_extra=0;
						for (int v=0;v<matchSimulatedPattern.getHeight();v++) for (int u=0;u<matchSimulatedPattern.getWidth();u++) {
							/*		            		System.out.println("+ v="+v+", u="+u);
		            		System.out.println(" + matchSimulatedPattern.targetUV[v][u] is "+((matchSimulatedPattern.targetUV[v][u]==null)?"null":"not null"));
		            		System.out.println(" + matchSimulatedPattern.pixelsUV[v][u] is "+((matchSimulatedPattern.pixelsUV[v][u]==null)?"null":"not null"));*/
							if ((matchSimulatedPattern.targetUV[v][u]!=null) &&(matchSimulatedPattern.pXYUV[v][u]!=null) ) {
								//			            		System.out.println("++ v="+v+", u="+u+" index="+index+" ("+size+"), index_extra="+index_extra+" ("+size_extra+")");

								if ((matchSimulatedPattern.targetUV[v][u]!=null) &&(matchSimulatedPattern.pXYUV[v][u]!=null) &&
										(matchSimulatedPattern.pXYUV[v][u][0]>=0.0) || (matchSimulatedPattern.pXYUV[v][u][1]>=0.0)) { // disregard negative sensor pixels
									if (
											(v>=matchSimulatedPattern.gridContrastBrightness[0].length) ||
											(u>=matchSimulatedPattern.gridContrastBrightness[0][0].length)){
										System.out.println(
												" matchSimulatedPattern.gridContrastBrightness[0].length="+matchSimulatedPattern.gridContrastBrightness[0].length+
												" matchSimulatedPattern.gridContrastBrightness[0][0].length="+matchSimulatedPattern.gridContrastBrightness[0][0].length+
												" v="+v+" u="+u);
									}
								}
								// setting dcd.gIP[numGridImage].pixelsUV[index] with rotated/shifted
								if (patternParameters.getXYZM(matchSimulatedPattern.targetUV[v][u][0],matchSimulatedPattern.targetUV[v][u][1],false,station)!=null) {
									dcd.gIP[numGridImage].pixelsXY[index][0]=matchSimulatedPattern.pXYUV[v][u][0];
									dcd.gIP[numGridImage].pixelsXY[index][1]=matchSimulatedPattern.pXYUV[v][u][1];
									dcd.gIP[numGridImage].pixelsUV[index][0]=matchSimulatedPattern.targetUV[v][u][0];
									dcd.gIP[numGridImage].pixelsUV[index][1]=matchSimulatedPattern.targetUV[v][u][1];
									dcd.gIP[numGridImage].pixelsXY[index][2]=matchSimulatedPattern.gridContrastBrightness[0][v][u]; // grid contrast
									dcd.gIP[numGridImage].pixelsXY[index][3]=matchSimulatedPattern.gridContrastBrightness[1][v][u]/dcd.gIP[numGridImage].intensityRange[0]; // red
									dcd.gIP[numGridImage].pixelsXY[index][4]=matchSimulatedPattern.gridContrastBrightness[2][v][u]/dcd.gIP[numGridImage].intensityRange[1]; // green
									dcd.gIP[numGridImage].pixelsXY[index][5]=matchSimulatedPattern.gridContrastBrightness[3][v][u]/dcd.gIP[numGridImage].intensityRange[2]; // blue
									index++;
								} else {
									dcd.gIP[numGridImage].pixelsXY_extra[index_extra][0]=matchSimulatedPattern.pXYUV[v][u][0];
									dcd.gIP[numGridImage].pixelsXY_extra[index_extra][1]=matchSimulatedPattern.pXYUV[v][u][1];
									dcd.gIP[numGridImage].pixelsUV_extra[index_extra][0]=matchSimulatedPattern.targetUV[v][u][0];
									dcd.gIP[numGridImage].pixelsUV_extra[index_extra][1]=matchSimulatedPattern.targetUV[v][u][1];
									dcd.gIP[numGridImage].pixelsXY_extra[index_extra][2]=matchSimulatedPattern.gridContrastBrightness[0][v][u]; // grid contrast
									dcd.gIP[numGridImage].pixelsXY_extra[index_extra][3]=matchSimulatedPattern.gridContrastBrightness[1][v][u]/dcd.gIP[numGridImage].intensityRange[0]; // red
									dcd.gIP[numGridImage].pixelsXY_extra[index_extra][4]=matchSimulatedPattern.gridContrastBrightness[2][v][u]/dcd.gIP[numGridImage].intensityRange[1]; // green
									dcd.gIP[numGridImage].pixelsXY_extra[index_extra][5]=matchSimulatedPattern.gridContrastBrightness[3][v][u]/dcd.gIP[numGridImage].intensityRange[2]; // blue
									index_extra++;
								}
							}
						}
						dcd.gIP[numGridImage].hintedMatch =(hintGridTolerance>0.0)?2:1; // orientation or both orientation and translation
						dcd.gIP[numGridImage].matchedPointers=rslt; // update number of matched pointers
						if ((dcd.gIP[numGridImage].hintedMatch>1) || (dcd.gIP[numGridImage].matchedPointers>0)) numSuccess++;
						// Update rotation/shift
						//matchSimulatedPattern
						int [] fileUVShiftRot=dcd.gIP[numGridImage].getUVShiftRot();
						int [] extraUVShiftRot=matchSimulatedPattern.getUVShiftRot(true); // last shift/rotation during matching pattern, correct for zero shift
						//						int [] extraDbg=matchSimulatedPattern.getUVShiftRot(false);
						int [] combinedUVShiftRot=MatchSimulatedPattern.combineUVShiftRot(fileUVShiftRot,extraUVShiftRot);
						dcd.gIP[numGridImage].setUVShiftRot(combinedUVShiftRot);
						System.out.println("applyHintedGrids(): dcd.gIP["+numGridImage+"].hintedMatch="+dcd.gIP[numGridImage].hintedMatch+
								" dcd.gIP["+numGridImage+"].matchedPointers="+dcd.gIP[numGridImage].matchedPointers+ " points:"+index+" extra points:"+index_extra);
						// testing rot/shift:
						String nonzero=((extraUVShiftRot[0]==0)&&(extraUVShiftRot[1]==0)&&(extraUVShiftRot[2]==0))?" ":"*";
						System.out.println("applyHintedGrids(): fileUVShiftRot=    "+fileUVShiftRot[0]+"/"+fileUVShiftRot[1]+":"+fileUVShiftRot[2]);
						System.out.println("                   "+nonzero+"extraUVShiftRot=   "+extraUVShiftRot[0]+"/"+extraUVShiftRot[1]+":"+extraUVShiftRot[2]);
						System.out.println("                    combinedUVShiftRot="+combinedUVShiftRot[0]+"/"+combinedUVShiftRot[1]+":"+combinedUVShiftRot[2]);
						//						System.out.println("                    extraDbg="+extraDbg[0]+"/"+extraDbg[1]+":"+extraDbg[2]);
					}
				}
			}
		}
		return numSuccess;
	}
	public void showSourceImage(int numGridImage){
		String source_path=fittingStrategy.distortionCalibrationData.gIP[numGridImage].source_path;
		if (source_path != null) {
			ImagePlus imp = new ImagePlus(source_path);
			imp.show();
		}
	}

	public int [][] getImageMarkers(int numGridImage){
		String source_path=fittingStrategy.distortionCalibrationData.gIP[numGridImage].source_path;
		if (source_path != null) {
			ImagePlus imp = new ImagePlus(source_path);
			imp.show();
			/*
			Thread msg_box_thread  = new Thread() {
   				@Override
				public void run() {
   					IJ.showMessage("Please place point markers on the "+imp.getTitle());
   				}
   			};
   			msg_box_thread.setPriority(Thread.MIN_PRIORITY);
   			msg_box_thread.start();
   			try {
   				msg_box_thread.join();
   			} catch (InterruptedException ie) {
   				throw new RuntimeException(ie);
   			}
	*/


//			IJ.showMessage("Please place point markers on the "+imp.getTitle());
			System.out.println("got it");
			PointRoi pointRoi = null;

			if (imp.getRoi() instanceof PointRoi) {
				pointRoi =  (PointRoi) imp.getRoi();
			} else {
				System.out.println("This image does not have point marks - please mark it in "+source_path);
				IJ.showMessage("This image does not have point marks - please mark it in "+source_path);
				return null;
				/*
				
				boolean mark_and_continue = IJ.showMessageWithCancel("Mark and Continue", "This image does not have point marks - please mark it in "+source_path);
				
				if (mark_and_continue) {
					imp = new ImagePlus(source_path);
					System.out.println("got it again!");
					pointRoi = null;
					if (imp.getRoi() instanceof PointRoi) {
						pointRoi =  (PointRoi) imp.getRoi();
					} else {
						System.out.println("This image does not have point marks - please mark it in "+source_path);
						IJ.showMessage("This image does not have point marks - please mark it in "+source_path);
						return null;
					}
				} else {
					return null;
				}
				*/
			}
			Point [] points = pointRoi.getContainedPoints();
			int [][] ipoints = new int [points.length][2];
			for (int n = 0; n < ipoints.length; n++) {
				ipoints[n][0] = points[n].x;
				ipoints[n][1] = points[n].y;
			}
			return ipoints;
		}
		return null;
	}


	public void showGridImage(int numGridImage){
		DistortionCalibrationData.GridImageParameters grid=fittingStrategy.distortionCalibrationData.gIP[numGridImage];
		boolean valid=false;
		int minU=0,maxU=0,minV=0,maxV=0;
		for (int i=0;i<grid.pixelsUV.length;i++){
			if (!valid){
				minU=grid.pixelsUV[i][0];
				minV=grid.pixelsUV[i][1];
				maxU=minU;
				maxV=minV;
				valid=true;
			} else {
				if (minU>grid.pixelsUV[i][0]) minU=grid.pixelsUV[i][0];
				if (minV>grid.pixelsUV[i][1]) minV=grid.pixelsUV[i][1];
				if (maxU<grid.pixelsUV[i][0]) maxU=grid.pixelsUV[i][0];
				if (maxV<grid.pixelsUV[i][1]) maxV=grid.pixelsUV[i][1];
			}
		}
		for (int i=0;i<grid.pixelsUV_extra.length;i++){
			if (!valid){
				minU=grid.pixelsUV_extra[i][0];
				minV=grid.pixelsUV_extra[i][1];
				maxU=minU;
				maxV=minV;
				valid=true;
			} else {
				if (minU>grid.pixelsUV_extra[i][0]) minU=grid.pixelsUV_extra[i][0];
				if (minV>grid.pixelsUV_extra[i][1]) minV=grid.pixelsUV_extra[i][1];
				if (maxU<grid.pixelsUV_extra[i][0]) maxU=grid.pixelsUV_extra[i][0];
				if (maxV<grid.pixelsUV_extra[i][1]) maxV=grid.pixelsUV_extra[i][1];
			}
		}
		String [] titles={"X","Y","U","V","valid","extra"};
		int height=maxV-minV+1;
		int width= maxU-minU+1;
//		System.out.println("showGridImage(): minU="+minU+" maxU="+maxU+" minV="+minV+" maxV="+maxV+" width="+width+" height="+height);
//		System.out.println("showGridImage(): grid.pixelsXY.length="+grid.pixelsXY.length+" grid.pixelsXY.length="+grid.pixelsXY.length);
		double [][] pixels=new double [titles.length][width*height];
		for (int i=0;i<pixels[0].length;i++) {
			pixels[0][i]=-1.0; // x
			pixels[1][i]=-1.0; // y
			pixels[2][i]= 0.0; // u
			pixels[3][i]= 0.0; // v
			pixels[4][i]=-1000.0; // valid
			pixels[5][i]=-1000.0; // extra
		}
		for (int i=0;i<grid.pixelsUV.length;i++){
			int u=grid.pixelsUV[i][0]-minU;
			int v=grid.pixelsUV[i][1]-minV;
			int index=u+width*v;
			pixels[0][index]=grid.pixelsXY[i][0];
			pixels[1][index]=grid.pixelsXY[i][1];
			pixels[2][index]=grid.pixelsUV[i][0];
			pixels[3][index]=grid.pixelsUV[i][1];
			pixels[4][index]=1000.0;
		}
		for (int i=0;i<grid.pixelsUV_extra.length;i++){
			int u=grid.pixelsUV_extra[i][0]-minU;
			int v=grid.pixelsUV_extra[i][1]-minV;
			int index=u+width*v;
			pixels[0][index]=grid.pixelsXY_extra[i][0];
			pixels[1][index]=grid.pixelsXY_extra[i][1];
			pixels[2][index]=grid.pixelsUV_extra[i][0];
			pixels[3][index]=grid.pixelsUV_extra[i][1];
			pixels[4][index]=1000.0;
		}
		ShowDoubleFloatArrays.showArrays(pixels, width, height,  true, "grid-"+numGridImage, titles);
	}



	public void manualGridHint(int imgNumber) {
		int [][] markers = getImageMarkers(imgNumber);
		if ((markers != null) && (markers.length > 0)) {
			if (markers.length > 1) {
				System.out.println("This image has multiple point marks - please remove extra");
				IJ.showMessage("This image has multiple point marks - please remove extra");
				return;
			}
			double [][] xyuv = new double [markers.length][4];
			for (int i =0; i < markers.length; i++) {
				xyuv[i][0] = markers[i][0];
				xyuv[i][1] = markers[i][1];
				xyuv[i][2] = lastUsedManualGridHint_UV[0]; // 16.5; // 15.5; // 0.5;
				xyuv[i][3] = lastUsedManualGridHint_UV[1]; //  0.5; // -8.5;//0.5;
			}
			GenericDialog gd=new GenericDialog("Specify U,V coordinates of the marker(s)");
			gd.addMessage("Center white (LWIR black) cell U=0.5, V=0.5");
			for (int n = 0; n < markers.length; n++) {
				String label = "Marker "+(n+1)+" (x="+markers[n][0]+", y="+markers[n][1];
				gd.addNumericField(label+" U", xyuv[n][2], 1, 5, "");
				gd.addNumericField(label+" V", xyuv[n][3], 1, 5, "");
			}
			gd.showDialog();
			if (gd.wasCanceled()) return;
			for (int i =0; i < markers.length; i++) {
				xyuv[i][2] = gd.getNextNumber();
				xyuv[i][3] = gd.getNextNumber();
			}

		    lastUsedManualGridHint_UV[0] = xyuv[0][2];
		    lastUsedManualGridHint_UV[1] = xyuv[0][3];
			
			
			// read grid image

			String grid_path=fittingStrategy.distortionCalibrationData.gIP[imgNumber].path;
			if (grid_path != null) {
				ImagePlus imp = new ImagePlus(grid_path);
				JP46_Reader_camera jp4_instance= new JP46_Reader_camera(false);
				jp4_instance.decodeProperiesFromInfo(imp);
				MatchSimulatedPattern.setPointersXYUV(imp, xyuv);
				updateGridToPointer(imp, xyuv);
				jp4_instance.encodeProperiesToInfo(imp);
				System.out.println("Updated "+grid_path);
				(new FileSaver(imp)).saveAsTiff(grid_path);
//				imp.show();
			}
			return;
		}
	}
	public void updateGridToPointer(ImagePlus imp_grid, double[][] xyuv) {
		ImageStack stack=imp_grid.getStack();
		if ((stack==null) || (stack.getSize()<4)) {
			String msg="Expected a 8-slice stack in "+imp_grid.getTitle();
			IJ.showMessage("Error",msg);
			throw new IllegalArgumentException (msg);
		}
		float [][] pixels=new float[stack.getSize()][]; // now - 8 (x,y,u,v,contrast, vignR,vignG,vignB
		for (int i=0;i<pixels.length;i++) pixels[i]= (float[]) stack.getPixels(i+1); // pixel X : negative - no grid here
		int width = imp_grid.getWidth();
		int height = imp_grid.getHeight();
		// start with translation only using xyuv[0][], may use full matching - same as laser pointers later
		int    indx_best = -1;
		double d2_best = Double.NaN;
		for (int indx = 0; indx < pixels[0].length; indx++) {
			double dx = pixels[0][indx] - xyuv[0][0];
			double dy = pixels[1][indx] - xyuv[0][1] ;
			double d2 = dx*dx + dy*dy;
			if (Double.isNaN(d2_best) || (d2 < d2_best)) {
				indx_best = indx;
				d2_best = d2;
			}
		}
		int ix0 = indx_best % width;
		int iy0 = indx_best / width;
		int half_range = 2; // was 1;
		PolynomialApproximation polynomialApproximation =new PolynomialApproximation(0);// no debug
//		double [][][] data = new double[9][3][];
		double [][][] data = new double[(2*half_range+1)*(2*half_range+1)][3][];
		int indx = 0;
		for (int idy = -half_range; idy <=half_range; idy++) {
			int iy = iy0+idy;
			for (int idx = -half_range; idx <= half_range; idx++) {
				int ix = ix0 + idx;
				data[indx][0] = new double[2];
				data[indx][1] = new double[2];
				data[indx][2] = new double[1];
				data[indx][0][0] = idx;
				data[indx][0][1] = idy;
				data[indx][2][0] = 0.0;
				if ((ix >= 0) && (ix < width) && (iy >= 0) && (iy < height)) {
					int offs = iy * width + ix;
					if ((pixels[0][offs] >= 0) && (pixels[1][offs] >= 0)) {
						data[indx][1][0] = pixels[0][offs] - xyuv[0][0];
						data[indx][1][1] = pixels[1][offs] - xyuv[0][1];
						data[indx][2][0] = 1.0;
					}
				}
				indx++;
			}
		}
		double [][] coeff = polynomialApproximation.quadraticApproximation(
				data,
				true); // force linear
		double [][] aA = {{coeff[0][0],coeff[0][1]},{coeff[1][0],coeff[1][1]}};
		double [][] aB = {{-coeff[0][2]},{-coeff[1][2]}};
		Matrix A = new Matrix(aA);
		Matrix B = new Matrix(aB);
		Matrix V = A.solve(B);
		double [] av = V.getColumnPackedCopy();
		double u, v; //  = xyuv[0][2]-()
		if (av[0] < 0) {
			av[0] += 1.0;
			ix0 -= 1;
		}
		if (av[1] < 0) {
			av[1] += 1.0;
			iy0 -= 1;
		}
		u = xyuv[0][2] - (
				(1-av[0])*(1-av[1]) * pixels[2][(iy0 + 0) * width + ix0 + 0]+
				(  av[0])*(1-av[1]) * pixels[2][(iy0 + 0) * width + ix0 + 1]+
				(1-av[0])*(  av[1]) * pixels[2][(iy0 + 1) * width + ix0 + 0]+
				(  av[0])*(  av[1]) * pixels[2][(iy0 + 1) * width + ix0 + 1]);
		v = xyuv[0][3] - (
				(1-av[0])*(1-av[1]) * pixels[3][(iy0 + 0) * width + ix0 + 0]+
				(  av[0])*(1-av[1]) * pixels[3][(iy0 + 0) * width + ix0 + 1]+
				(1-av[0])*(  av[1]) * pixels[3][(iy0 + 1) * width + ix0 + 0]+
				(  av[0])*(  av[1]) * pixels[3][(iy0 + 1) * width + ix0 + 1]);
		int idu = (int)Math.round(u);
		int idv = (int)Math.round(v);
		// Verify that idy+idv - even number
		if (((idu + idv) & 1) != 0) {
			String msg = "Incorrect shift - u="+u+", v="+v+", idu="+idu+", idv="+idv+", idu+idv="+(idu+idv)+" SHOULD BE EVEN!";
			System.out.println(msg);
			IJ.showMessage(msg);
		}

		for (int i = 0; i < pixels[2].length; i++) {
			pixels[2][i] += idu;
			pixels[3][i] += idv;
		}
	}

	public void showGridAndHint(){
		GenericDialog gd=new GenericDialog("Show selected grid and/or hint grid");
		gd.addNumericField("Grid Image index", 0,0);
		gd.addCheckbox("Show source image (if available)", true);
		gd.addCheckbox("Show grid image", true);
		gd.addCheckbox("Show hint grid", true);
		gd.addCheckbox("Use imageSet data if available (unchecked - camera data)", true);

		gd.showDialog();
		if (gd.wasCanceled()) return;
		int numGridImage= (int) gd.getNextNumber();
		boolean showGrid=gd.getNextBoolean();
		boolean showSource=gd.getNextBoolean();
		boolean showHint=gd.getNextBoolean();
		boolean useSetData=gd.getNextBoolean();
		IJ.showStatus("grid: "+((fittingStrategy.distortionCalibrationData.gIP[numGridImage].path==null)?"":fittingStrategy.distortionCalibrationData.gIP[numGridImage].path));
//		showStatus("grid: "+((fittingStrategy.distortionCalibrationData.gIP[numGridImage].path==null)?"":fittingStrategy.distortionCalibrationData.gIP[numGridImage].path),0);

        if (showGrid)	showGridImage(numGridImage);
        if (showSource)	showSourceImage(numGridImage);
//        if (showSource)	getImageMarkers(numGridImage);
        if (showHint)	calcAndShowHintGrid(numGridImage,useSetData);
	}



	public void calcAndShowHintGrid(int numGridImage, boolean useSetData){
		double [] goniometerTiltAxial=fittingStrategy.distortionCalibrationData.getImagesetTiltAxial(numGridImage);
		if ((goniometerTiltAxial==null) || Double.isNaN(goniometerTiltAxial[0])  || Double.isNaN(goniometerTiltAxial[1])){
			if (this.debugLevel>0)System.out.println("No goniometer orientation is available for image # "+numGridImage+" - "+fittingStrategy.distortionCalibrationData.gIP[numGridImage].path);
			GenericDialog gd=new GenericDialog("Specify camera orientation (channel"+fittingStrategy.distortionCalibrationData.gIP[numGridImage].channel+")");
			gd.addMessage("No goniometer orientation is available for image # "+numGridImage+" - "+fittingStrategy.distortionCalibrationData.gIP[numGridImage].path+
			", please specify orientation manually");
			gd.addNumericField("Camera tilt (0 - vertical, >0 looking above horizon on the target)", 0.0, 1,6,"degrees");
			gd.addNumericField("Camera axial (0 - subcamera 0 looking to the target, >0 - rotated clockwise)", 0.0, 1,6,"degrees");
			gd.addNumericField("Camera inter-axis angle (from 90) ", 0.0, 1,6,"degrees");
			gd.showDialog();
			if (gd.wasCanceled()) return;
			goniometerTiltAxial=new double[3];
			goniometerTiltAxial[0]=      gd.getNextNumber();
			goniometerTiltAxial[1]=      gd.getNextNumber();
			goniometerTiltAxial[2]=      gd.getNextNumber();
		}
		double [][][] hintGrid=estimateGridOnSensor(
				fittingStrategy.distortionCalibrationData.getImageStation(numGridImage), // station number
				fittingStrategy.distortionCalibrationData.gIP[numGridImage].channel,
				goniometerTiltAxial[0], // Tilt, goniometerHorizontal
				goniometerTiltAxial[1],  // Axial,goniometerAxial
				goniometerTiltAxial[2],  // inter-axis angle
				(useSetData?fittingStrategy.distortionCalibrationData.gIP[numGridImage].getSetNumber():-1),
				true // filter border
				);
		if (hintGrid == null) {
			String msg = "hintGrid is null";
			IJ.showMessage("Error",msg);
			System.out.println(msg);
			return;
		}
		showHintGrid(hintGrid,"hint-"+numGridImage);

	}
	public void showHintGrid(double [][][] hintGrid){
		showHintGrid(hintGrid,"hintGrid");
	}

	public void showHintGrid(double [][][] hintGrid, String title){
		double [][] pixels=new double[4][hintGrid.length*hintGrid[0].length];
		int index=0;
		String [] titles={"pixel-X","pixel-Y","grid-U","grid-V"};
		for (int v=0; v<hintGrid.length;v++) for (int u=0;u<hintGrid[v].length;u++){
			if (hintGrid[v][u]!=null){
				for (int i=0; i<4;i++)	pixels[i][index]=hintGrid[v][u][i];
			} else {
				for (int i=0; i<4;i++)	pixels[i][index]=0;
			}
			index++;
		}
		ShowDoubleFloatArrays.showArrays(pixels, hintGrid[0].length, hintGrid.length,  true, title, titles);
	}

	/**
	 * Calculate grid on sensor using current camera parameters (including goniometer angles), sub-camera number
	 * @param subCamera
	 * @return grid array [v][u][0- x,  1 - y, 2 - u, 3 - v]
	 */
	/*
	 // wrong, orientation depends on timestamp
	public double [][][] estimateGridOnSensor(
			int subCamera){
		double [] parVector=fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getParametersVector(subCamera);
		return estimateGridOnSensor(
				subCamera,
				parVector[fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getGoniometerHorizontalIndex()],
				parVector[fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getGoniometerAxialIndex()]);
	}
	*/


	public LensDistortionParameters setupLensDistortionParameters(
			int numImg,
			int debugLevel){     // Axial - may be Double.NaN

		LensDistortionParameters lensDistortionParameters = new LensDistortionParameters (
				this.fittingStrategy.distortionCalibrationData.isTripod(),
				this.fittingStrategy.distortionCalibrationData.isCartesian(),
	    		this.fittingStrategy.distortionCalibrationData.getPixelSize(numImg),
	    		this.fittingStrategy.distortionCalibrationData.getDistortionRadius(numImg),
	            null, //double [][] interParameterDerivatives, //partial derivative matrix from subcamera-camera-goniometer to single camera (12x21) if null - just values, no derivatives
	            this.fittingStrategy.distortionCalibrationData.getParameters(numImg), //parVector,
	    		null, //boolean [] mask, // calculate only selected derivatives (all parVect values are still
	    		debugLevel
				);
		return lensDistortionParameters;
	}

	public LensDistortionParameters setupLensDistortionParameters(
			int stationNumber,
			int subCamera,
			double goniometerHorizontal, // Tilt - may be Double.NaN
			double goniometerAxial,
			int debugLevel){     // Axial - may be Double.NaN
		double [] parVector=fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getParametersVector(stationNumber,subCamera);
		int goniometerHorizontalIndex=fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getGoniometerHorizontalIndex();
		int goniometerAxialIndex=fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getGoniometerAxialIndex();
		if (!Double.isNaN(goniometerHorizontal))parVector[goniometerHorizontalIndex]=goniometerHorizontal;
		if (!Double.isNaN(goniometerAxial))parVector[goniometerAxialIndex]=goniometerAxial;
		LensDistortionParameters lensDistortionParameters = new LensDistortionParameters (
				this.fittingStrategy.distortionCalibrationData.isTripod(),
				this.fittingStrategy.distortionCalibrationData.isCartesian(),
	    		this.fittingStrategy.distortionCalibrationData.getPixelSize(stationNumber, subCamera),
	    		this.fittingStrategy.distortionCalibrationData.getDistortionRadius(stationNumber, subCamera),
	            null, //double [][] interParameterDerivatives, //partial derivative matrix from subcamera-camera-goniometer to single camera (12x21) if null - just values, no derivatives
	    		parVector,
	    		null, //boolean [] mask, // calculate only selected derivatives (all parVect values are still
	    		debugLevel
				);
		return lensDistortionParameters;
	}

	/**
	 * Calculate grid projection to pixel X, Y (not counting sensor correction (add?) and grid photometrics
	 * @param lensDistortionParameters LensDistortionParameters instance created for particular image with setupLensDistortionParameters()
	 * @param numImg image number
	 * @param u grid U (signed, 0 in the center)
	 * @param v grid V (signed, 0 in the center)
	 * @return [7] {pX,pY,grid mask (binary), grid R, grid G, grid B, alpha}
	 */
	public double [] reprojectGridNode(
			LensDistortionParameters lensDistortionParameters,
			int numImg,
			int u, // grid signed u,v
			int v){

		
		int subCamera=   this.fittingStrategy.distortionCalibrationData.gIP[numImg].channel;
		int sensorWidth=fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getSensorWidth(subCamera);
		int sensorHeight=fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getSensorHeight(subCamera);
//		double maxRelativeRadius=this.hintedMaxRelativeRadius; // make adjustable
		double maxRelativeRadius=hintedMaxRelativeRadiusToDiagonal * Math.sqrt(sensorWidth * sensorWidth + sensorHeight*sensorHeight)/ sensorWidth;
		
		return  reprojectGridNode(
				lensDistortionParameters,
				numImg,
				u, // grid signed u,v
				v,
		       	maxRelativeRadius);
	}
	public double [] reprojectGridNode(
			LensDistortionParameters lensDistortionParameters,
			int numImg,
			int u, // grid signed u,v
			int v,
	       	double maxRelativeRadius //=2.0;
	){
		int debugThreshold=1;
		int nChn=   this.fittingStrategy.distortionCalibrationData.gIP[numImg].channel;
		int station=this.fittingStrategy.distortionCalibrationData.gIP[numImg].getStationNumber();
//		if (!lensDistortionParameters.isTargetVisible(false)) return null; // camera is looking away from the target (does not mean target is in FOV)
//		double [][][] patternGeometry=this.patternParameters.getGeometry(); // [v][u]{x,y,z,alpha} - no photometric
		double [] result= new double[7];
			double [] XYZMP=this.patternParameters.getXYZMP( // null pointer
					u,
					v,
					station,
					nChn,
					false);
			if (XYZMP==null) return null;
			// project the target point to this sensor
			double [][]pXY=  lensDistortionParameters.calcPartialDerivatives(
					XYZMP[0], // target point horizontal, positive - right,  mm
					XYZMP[1], // target point vertical,   positive - down,  mm
					XYZMP[2], // target point horizontal, positive - away from camera,  mm
					maxRelativeRadius, //
					false); // calculate derivatives, false - values only (NaN for behind points - only when false here)
			if (Double.isNaN(pXY[0][0])) {
				if (this.debugLevel>debugThreshold){
					System.out.println("reprojectGridNode(...,"+numImg+","+u+","+"v"+") - point behind the sensor");
				}
				return null; // point behind camera
			}
			result[0]=pXY[0][0];
			result[1]=pXY[0][1];
			result[2]=XYZMP[3]; // binary mask
			result[3]=XYZMP[4]; // R
			result[4]=XYZMP[5]; // G
			result[5]=XYZMP[6]; // B
			result[6]=XYZMP[7]; // alpha
// get photometrics here

		return result;
	}

	/**
	 * Apply sensor correction to the projected grid (generated by estimateGridOnSensor())
	 * @param gridOnSensor array [v][u][0- x,  1 - y, 2 - targetAbsolute-u, 3 - targetAbsolute-v]
	 * @param subCamera channel number
	 * @return true if the correction was applied (in-place) false if no correction is available
	 */
	public boolean correctGridOnSensor(
			double [][][] gridOnSensor,
			int subCamera){
		if (this.pixelCorrection==null) return false;
		for (double [][] row:gridOnSensor) for (double [] cell:row) if ((cell!=null) && (cell.length>1)){
			double [] corrXYARGB=interpolateCorrectionVector ( // vector of {corrX, corrY, alpha, flatfield_red, flatfield_green, flatfield_blue}
					subCamera, //int chnNum,
					cell[0], //double px,
					cell[1]); //double py)
			cell[0]+=corrXYARGB[0]; // measured-> corrected : subtract, projected->simulated:add;
			cell[1]+=corrXYARGB[1]+0.0; // Debugging by adding +1.0!!
		}
//		System.out.println("================== Added +0.0 to pixel y for debugging purposes! =====================");
		return true;
	}

	/**
	 * Calculate grid on sensor using current Camera parameters, sub-camera number and the two goniometer angles
	 * @param stationNumber
	 * @param subCamera
	 * @param goniometerHorizontal
	 * @param goniometerAxial
	 * @param imageSet - if >=0 - use this set number data  instead of the camera data)
	 * @return grid array [v][u][0- x,  1 - y, 2 - u, 3 - v]
	 */

	// TODO:calcInterParamers() -> lensDistortionParameters.lensCalcInterParamers
	public double [][][] estimateGridOnSensor( // not yet thread safe
			int stationNumber,
			int subCamera,
			double goniometerHorizontal, // Tilt
			double goniometerAxial,     // Axial
			double goniometerInterAxis,     // interAxisAngle
			int  imageSet,
			boolean filterBorder){
		int sensorWidth=fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getSensorWidth(subCamera);
		int sensorHeight=fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getSensorHeight(subCamera);
//		double maxRelativeRadius=this.hintedMaxRelativeRadius; // make adjustable
		double maxRelativeRadius=hintedMaxRelativeRadiusToDiagonal * Math.sqrt(sensorWidth * sensorWidth + sensorHeight*sensorHeight)/ sensorWidth;
		// 1.1 is sufficient
//		double maxRelativeRadius= 2.0*Math.sqrt(sensorWidth * sensorWidth + sensorHeight*sensorHeight)/ sensorWidth;
		int debugThreshold=2;
		// Get parameter vector (22) for the selected sensor, current Eyesisparameters and specified orientation angles
		double [] parVector=fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getParametersVector(stationNumber,subCamera);
		if ((imageSet>=0) &&
				(this.fittingStrategy.distortionCalibrationData.gIS!=null) &&
				(this.fittingStrategy.distortionCalibrationData.gIS[imageSet]!=null)){
			this.fittingStrategy.distortionCalibrationData.gIS[imageSet].updateParameterVectorFromSet(parVector);
		}
		if (!Double.isNaN(goniometerHorizontal)) {
			int goniometerHorizontalIndex=fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getGoniometerHorizontalIndex();
			parVector[goniometerHorizontalIndex]=goniometerHorizontal;
		}
		if (!Double.isNaN(goniometerAxial)) {
			int goniometerAxialIndex=fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getGoniometerAxialIndex();
			parVector[goniometerAxialIndex]=     goniometerAxial;
		}
		if (!Double.isNaN(goniometerInterAxis)) {
			int goniometerInterAxisAngleIndex=fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getInterAxisAngleIndex();
			parVector[goniometerInterAxisAngleIndex]=  goniometerInterAxis;
		}
//		/interAxis
		System.out.println("estimateGridOnSensor(): subCamera="+subCamera+", goniometerHorizontal="+goniometerHorizontal+", goniometerAxial="+goniometerAxial);
		this.lensDistortionParameters.lensCalcInterParamers(
				this.lensDistortionParameters, // 22-long parameter vector for the image
				this.fittingStrategy.distortionCalibrationData.isTripod(),
				this.fittingStrategy.distortionCalibrationData.isCartesian(),
	    		this.fittingStrategy.distortionCalibrationData.getPixelSize(stationNumber, subCamera),
	    		this.fittingStrategy.distortionCalibrationData.getDistortionRadius(stationNumber, subCamera),
				null, // this.interParameterDerivatives, // [22][]
				parVector,
				null); // if no derivatives, null is OK


		if (!lensDistortionParameters.isTargetVisible(this.debugLevel>0)) {
			if (this.debugLevel>debugThreshold) System.out.println("Camera is looking away from the target");
//			return null; // camera is looking away from the target (does not mean target is in FOV)
		}
		double [][][] patternGeometry=this.patternParameters.getGeometry(); // [v][u]{x,y,z,alpha} - no photometric
		double [][][] result= new double[patternGeometry.length][patternGeometry[0].length][4];
		int visibleCells=0;
		double [][] debugPixels=null;
		String [] debugTitles={"pX","pY","X","Y","Z","mask"};
		if (this.debugLevel>debugThreshold){
			debugPixels=new double [6][patternGeometry.length*patternGeometry[0].length];
			for (int c=0;c<debugPixels.length;c++) for (int i=0;i<debugPixels[c].length;i++) debugPixels[c][i]=Double.NaN;
		}
		// was bug cased by +/- infinity (and sometimes numbers falling into the sensor range) when the image plane intersected target
		// simple fix - remove pixels with too few neighbors (maybe just all border pixels?


		for (int v=0;v<result.length;v++) {
			for (int u=0;u<result[v].length;u++){
				int [] iUV=this.patternParameters.uvIndicesToUV (u, v);
				if (iUV==null) {
					result[v][u]=null;
				} else {
					double [] XYZM=this.patternParameters.getXYZM(iUV[0],iUV[1],stationNumber);
					// project the target point to this sensor
					double [][]pXY=  this.lensDistortionParameters.calcPartialDerivatives(
							XYZM[0], // target point horizontal, positive - right,  mm
							XYZM[1], // target point vertical,   positive - down,  mm
							XYZM[2], // target point horizontal, positive - away from camera,  mm
							maxRelativeRadius,
							false); // calculate derivatives, false - values only (NaN for behind points - only when false here)
					// verify the grid is inside the sensor area (may use sensor mask later too? probably not needed)
					// Now NaN if point is behind the sensor
					if (Double.isNaN(pXY[0][0]) || (pXY[0][0]<0) || (pXY[0][0]>=sensorWidth) || (pXY[0][1]<0) || (pXY[0][1]>=sensorHeight)){
						if (this.debugLevel>debugThreshold){
							System.out.println("--- estimateGridOnSensor():v="+v+" u="+u+" X="+XYZM[0]+" Y="+XYZM[1]+" Z="+XYZM[2]+" M="+XYZM[3]+
									" pXY[0][0]="+pXY[0][0]+", pXY[0][1]="+pXY[0][1]+", iUV[0]="+iUV[0]+", iUV[1]="+iUV[1]);
						}
						result[v][u]=null;
					} else {
						double [] resultCell={pXY[0][0],pXY[0][1],iUV[0],iUV[1]};
						result[v][u]=resultCell;
						if (this.debugLevel>debugThreshold){
							System.out.println("+++ estimateGridOnSensor():v="+v+" u="+u+" X="+XYZM[0]+" Y="+XYZM[1]+" Z="+XYZM[2]+" M="+XYZM[3]+
									" pXY[0][0]="+pXY[0][0]+", pXY[0][1]="+pXY[0][1]+", iUV[0]="+iUV[0]+", iUV[1]="+iUV[1]);
						}
						visibleCells++;
					}
					if (this.debugLevel>debugThreshold){
						int uv=u+v*result[v].length;
						debugPixels[0][uv]=pXY[0][0];
						debugPixels[1][uv]=pXY[0][1];
						debugPixels[2][uv]=XYZM[0];
						debugPixels[3][uv]=XYZM[1];
						debugPixels[4][uv]=XYZM[2];
					}
				}
			}
		}
		if (filterBorder){
			// now filter border nodes
			boolean [] mask= new boolean [patternGeometry.length*patternGeometry[0].length];
			int index=0;
			for (int v=0;v<result.length;v++) for (int u=0;u<result[v].length;u++){
				mask [index++]=(result[v][u]!=null) &&
						((v==0) || (result[v-1][u]!=null)) &&
						((v==(result.length-1)) || (result[v+1][u]!=null)) &&
						((u==0) || (result[v][u-1]!=null))&&
						((u==(result[v].length-1)) || (result[v][u+1]!=null));
			}
			index=0;
			for (int v=0;v<result.length;v++) for (int u=0;u<result[v].length;u++){
				if (!mask[index++]) result[v][u]=null;
			}
		}
		if (this.debugLevel>debugThreshold){
			for (int v=0;v<result.length;v++) for (int u=0;u<result[v].length;u++){
				int uv=u+v*result[v].length;
				debugPixels[5][uv]=(result[v][u]!=null)?3000:-3000; // masked
			}
			ShowDoubleFloatArrays.showArrays(
					debugPixels,
					result[0].length,
					result.length,
					true,
					"Hinted-All",
					debugTitles);
		}
		if (this.debugLevel>0) {
			System.out.println("Grid in the FOV of the subcamera "+subCamera+
					" tilt="+goniometerHorizontal+" axial="+goniometerAxial+" has "+visibleCells+" cells");
		}
		if (visibleCells==0) return null; // no grid cells in FOV
		return result;
	}

    public void debugCompareInterparameterDerivatives(
    		double [] vector,
    		int imgNum,
    		double delta){
		if (this.debugLevel>1) {
			System.out.println("debugCompareInterparameterDerivatives(vector, imgNum="+imgNum+", delta="+delta+")");
			for (int ii=0;ii<vector.length;ii++) System.out.println(ii+": "+vector[ii]);
		}
		int numImg=fittingStrategy.distortionCalibrationData.getNumImages();
		if (imgNum<0){ // find first selected image
			boolean [] selectedImages=fittingStrategy.selectedImages();
			imgNum=0;
			while ((imgNum<numImg) && (!selectedImages[imgNum])) imgNum++;
		}
		if (imgNum>=numImg){
			IJ.showMessage("No images found for this fitting strategy");
			return; // no images found
		}
		double [] imgVector=fittingStrategy.getImageParametersVector(imgNum, vector); //this.currentVector);
		boolean [] imgMask= new boolean[imgVector.length];
		for (int i=0;i<imgMask.length;i++) imgMask[i]=true;
		this.lensDistortionParameters.lensCalcInterParamers(
				this.lensDistortionParameters,
				this.fittingStrategy.distortionCalibrationData.isTripod(),
				this.fittingStrategy.distortionCalibrationData.isCartesian(),
	    		this.fittingStrategy.distortionCalibrationData.getPixelSize(imgNum),
	    		this.fittingStrategy.distortionCalibrationData.getDistortionRadius(imgNum),
				this.interParameterDerivatives, // [22][]
				imgVector,
				imgMask); // calculate only selected derivatives (all parVect values are still
//				true); // calculate this.interParameterDerivatives -derivatives array (false - just this.values)
// reorder derivatives to match lensDistortionParameters.getExtrinsicVector(); (dist,x0,y0,yaw,pitch,roll)
//					double [] parameterVector0=lensDistortionParameters.getAllVector();
		double [] values=lensDistortionParameters.getExtrinsicVector();
		double [][] derivatives_true = new double [this.lensDistortionParameters.getNumInputs()][6];
		for (int i=0;i<this.lensDistortionParameters.getNumInputs();i++){
			derivatives_true[i][0]=this.interParameterDerivatives[i][2]; // d distance /d vector[i]
			derivatives_true[i][1]=this.interParameterDerivatives[i][0]; // d x0 /d vector[i]
			derivatives_true[i][2]=this.interParameterDerivatives[i][1]; // d y0 /d vector[i]
			derivatives_true[i][3]=this.interParameterDerivatives[i][3]; // d jaw /d vector[i]
			derivatives_true[i][4]=this.interParameterDerivatives[i][4]; // d pitch /d vector[i]
			derivatives_true[i][5]=this.interParameterDerivatives[i][5]; // d roll /d vector[i]
		}
		double [][] derivatives_delta = new double [this.lensDistortionParameters.getNumInputs()][values.length];
		for (int i=0;i<this.lensDistortionParameters.getNumInputs();i++){
			double [] vector_delta=imgVector.clone();
			vector_delta[i]+=delta;
			this.lensDistortionParameters.lensCalcInterParamers(
					this.lensDistortionParameters,
					this.fittingStrategy.distortionCalibrationData.isTripod(),
					this.fittingStrategy.distortionCalibrationData.isCartesian(),
		    		this.fittingStrategy.distortionCalibrationData.getPixelSize(imgNum),
		    		this.fittingStrategy.distortionCalibrationData.getDistortionRadius(imgNum),
					null, // this.interParameterDerivatives, // just values, no derivatives
					vector_delta,
					imgMask);
//					false); // just values, no derivatives
			double [] values_delta=lensDistortionParameters.getExtrinsicVector();
			for (int j=0;j<derivatives_delta[i].length;j++) derivatives_delta[i][j]=(values_delta[j]-values[j])/delta;
		}
		String [] lensParNames = lensDistortionParameters.getExtrinsicNames();
	    String header="#\tphysical/lens\t ";
	    for (int i=0;i<lensParNames.length;i++)header+="\t"+lensParNames[i];
	    StringBuffer sb = new StringBuffer();
	    for (int parNum=0;parNum<this.lensDistortionParameters.getNumInputs();parNum++){
	    	sb.append(parNum+"\t"+fittingStrategy.distortionCalibrationData.descrField(parNum,0)+"\tderivative");
	    	for (int i=0;i<lensParNames.length;i++) sb.append("\t"+derivatives_true[parNum][i]);
	    	sb.append("\n");
	    	sb.append("\t \tdelta");
	    	for (int i=0;i<lensParNames.length;i++) sb.append("\t"+derivatives_delta[parNum][i]);
	    	sb.append("\n");
	    	sb.append("\t \tdifference");
	    	for (int i=0;i<lensParNames.length;i++) sb.append("\t"+(derivatives_true[parNum][i]-derivatives_delta[parNum][i]));
	    	sb.append("\n");
	    	sb.append("---\t---\t---");
	    	for (int i=0;i<lensParNames.length;i++) sb.append("\t---");
	    	sb.append("\n");
	    }
	    new TextWindow("Comparisison of the interparameter dcerivatives (true and compared as deltas)", header, sb.toString(), 500,900);
    }
// after stepping back - no need to rerun calculateFxAndJacobian(false), just keep
	/**
	 *  Calculates f(X) and optionally Jacobian for the current parameters
	 *  @parameter vector - parameter vector to be used
	 *  @parameter calcJacobian  if true, calculates Jacobian as this.jacobian
	 *  @return  f(X) - pixel coordinates for each (visible) grid pattern node for current parameters this.currentVector
	 *   as a 1-d array that alternates pixel-X and pixel-Y for all images
	 *   NOTE: this one is not thread safe
	 */
	public double [] calculateFxAndJacobian(
			double [] vector,
			boolean calcJacobian){ // when false, modifies only this.lensDistortionParameters.*
		if (vector==null) {
			calcJacobian=false;
//			vector = new double[0];
		}
		// TODO: verify classes/arrays exist?
        int doubleNumAllPoints=this.Y.length; // all points in all images multiplied by 2 (x and y error are separate)
		int fittedParNumber=(vector==null)?0:vector.length; //this.currentVector.length;
		double [] vectorFX=new double[doubleNumAllPoints];
//		this.fX=new double[doubleNumAllPoints];
		if (this.debugLevel>2) {
			System.out.println("calculateFxAndJacobian(), calcJacobian="+calcJacobian+" D3304 + this.debugLevel="+this.debugLevel);
			if (vector!=null) {
			  for (int ii=0;ii<vector.length;ii++) System.out.println(ii+": "+vector[ii]);
			} else {
				System.out.println("calculateFxAndJacobian() : vector==null");
			}
		}
		if (calcJacobian) {
			this.jacobian=new double[fittedParNumber][doubleNumAllPoints];
			for (int i=0;i<fittedParNumber;i++) for (int j=0;j<doubleNumAllPoints;j++) this.jacobian[i][j]=0.0;
		}
		int numImg=fittingStrategy.distortionCalibrationData.getNumImages();
		boolean [] selectedImages=fittingStrategy.selectedImages();
		int index=0;
		IJ.showProgress(0);
		for (int imgNum=0;imgNum<numImg;imgNum++) if (selectedImages[imgNum]) {
//			initialize arrays for parameters and derivatives conversion
			double [] imgVector=fittingStrategy.getImageParametersVector(imgNum, vector); // null is OK now
			boolean [] imgMask=null;
			int []     imgMap= null;
			if (calcJacobian) {
				imgMask= fittingStrategy.getImageParametersVectorMask(imgNum);
				int []     imgRMap=  fittingStrategy.getImageParametersVectorReverseMap(imgNum);
				imgMap=new int[vector.length];
				for (int i=0;i<imgMap.length;i++) imgMap[i]=-1;
				for (int i=0;i<imgRMap.length;i++) if (imgRMap[i]>=0)imgMap[imgRMap[i]]=i;
			}

// Calculate/set  this.lensDistortionParameters class, so it will calculate values/derivatives correctly)
// and this.interParameterDerivatives
//			if (this.debugLevel>1) {
			if (this.debugLevel>2) {
				System.out.println("calculateFxAndJacobian(), imgNum="+imgNum+" calcInterParamers(): (D3336)");
			}
			this.lensDistortionParameters.debugLevel=this.debugLevel;
			this.lensDistortionParameters.lensCalcInterParamers(
					this.lensDistortionParameters,
					this.fittingStrategy.distortionCalibrationData.isTripod(),
					this.fittingStrategy.distortionCalibrationData.isCartesian(),
		    		this.fittingStrategy.distortionCalibrationData.getPixelSize(imgNum),
		    		this.fittingStrategy.distortionCalibrationData.getDistortionRadius(imgNum),
					calcJacobian?this.interParameterDerivatives:null, // [22][]
					imgVector,
					imgMask); // imgMask may be null if no derivativescalculate only selected derivatives (all parVect values are still
			int numPoints=fittingStrategy.distortionCalibrationData.getImageNumPoints(imgNum);
			if (this.debugLevel>2) {
				System.out.println("calculateFxAndJacobian(), numPoints="+numPoints+" (imgNum="+imgNum+")");
			}
// iterate through points, for each calculate pixelx, pixely and derivatives
			for (int pointNum=0;pointNum<numPoints;pointNum++){
				int fullIndex=index+pointNum;
				if (fullIndex>=this.targetXYZ.length){
					System.out.println("BUG: calculateFxAndJacobian() imgNum="+imgNum+" pointNum="+pointNum+" fullIndex="+fullIndex+" this.targetXYZ.length="+this.targetXYZ.length);
				}
				double [][]derivatives15=  lensDistortionParameters.calcPartialDerivatives( // [NaN, NaN]
						this.targetXYZ[fullIndex][0], // target point horizontal, positive - right,  mm
						this.targetXYZ[fullIndex][1], // target point vertical,   positive - down,  mm
						this.targetXYZ[fullIndex][2], // target point horizontal, positive - away from camera,  mm
						calcJacobian); // calculate derivatives, false - values only
	       		if (this.debugLevel>3) {
	    			System.out.println(fullIndex+": calculateFxAndJacobian->calcPartialDerivatives("+IJ.d2s(targetXYZ[fullIndex][0],2)+","+
	    					IJ.d2s(targetXYZ[fullIndex][1],2)+","+
	    					IJ.d2s(targetXYZ[fullIndex][2],2)+" ("+calcJacobian+") -> "+
	    					IJ.d2s(derivatives15[0][0],2)+"/"+IJ.d2s(derivatives15[0][1],2));
	    			String all="derivatives15: D3365";
	    			for (int ii=0;ii<derivatives15.length;ii++) all+=" "+ii+":"+IJ.d2s(derivatives15[ii][0],3)+"/"+IJ.d2s(derivatives15[ii][1],3);
	    			System.out.println(all);
	    		}
				vectorFX[2*fullIndex]=  derivatives15[0][0];
				vectorFX[2*fullIndex+1]=derivatives15[0][1];
				if (calcJacobian) {
					double [][]derivatives = lensDistortionParameters.reorderPartialDerivatives(derivatives15);
		       		if (this.debugLevel>3) {
		    			String all="derivatives:";
		    			for (int ii=0;ii<derivatives.length;ii++) all+=" "+ii+":"+IJ.d2s(derivatives[ii][0],3)+"/"+IJ.d2s(derivatives[ii][1],3);
		    			System.out.println(all);
		    		}
					for (int i=0;i<this.jacobian.length;i++) if (imgMap[i]>=0){
						double sX=0,sY=0;
						for (int k=0;k<derivatives.length;k++){
							sX+=this.interParameterDerivatives[imgMap[i]][k]*derivatives[k][0];
							sY+=this.interParameterDerivatives[imgMap[i]][k]*derivatives[k][1];
						}
						this.jacobian[i][2*fullIndex]=  sX;
						this.jacobian[i][2*fullIndex+1]=sY;
					}
				}
			}
			index+=numPoints;
			IJ.showProgress(imgNum, numImg-1);
		}
//		IJ.showProgress(0); not needed, will turn off automatically

		return vectorFX;
	}

	/**
	 * Calculate FX and (optionally) Jacobian for one image. FX is a single vector for all images, jacobian - only for one (to save on memory usage)
	 * @param numImage    number of image being processed
	 * @param vector      parameters vector
	 * @param patternXYZ  X,Y,Z of the physical target for each node of each image (TODO: memory may be reduced)
	 * @param vectorFX    Vector to be filled here , twice length as patternXYZ (x and y alternating)
	 * @param imageStartIndex  start index in patternXYZ array (length - difference to the next, includes extra last element)
	 * @param lensDistortionParameters LensDistortionParameters class instance (may be reused between calls)
	 * @param calcJacobian calculate Jacobian matrix (if false - only FX)
	 * @return partial Jacobian matrix, number of rows= vector.length, number of columns - 2*indexCount
	 *   NOTE: this one is thread safe
	 */

	public double [][] calculatePartialFxAndJacobian(
			final int numImage,      // number of grid image
			final double [] vector,  // parameters vector
			final double [][] patternXYZ, // this.targetXYZ
			final double [] vectorFX,     // non-overlapping segments will be filled
			final int []  imageStartIndex, // start index in patternXYZ array (length - difference to the next, includes extra last element)
			final LensDistortionParameters lensDistortionParameters, // initialize one per each thread? Or for each call?
			boolean calcJacobian){ // when false, modifies only this.lensDistortionParameters.*
		final int    indexStart=imageStartIndex[numImage];      // start index in patternXYZ array
		final int    indexCount=imageStartIndex[numImage+1]-imageStartIndex[numImage]; // number of nodes in the current grid image
		int fittedParNumber=vector.length; //this.currentVector.length;
		if (this.debugLevel>3) {
			System.out.println("calculatePartialFxAndJacobian(), calcJacobian="+calcJacobian+" indexStart="+indexStart+" indexCount="+indexCount);
			for (int ii=0;ii<vector.length;ii++) System.out.println("vector["+ii+"]: "+vector[ii]);
		}
		boolean [] imgMask= fittingStrategy.getImageParametersVectorMask(numImage);        // thread safe
		int []     imgRMap=  fittingStrategy.getImageParametersVectorReverseMap(numImage); // thread safe
		int []     imgMap=new int[vector.length];
		for (int i=0;i<imgMap.length;i++) imgMap[i]=-1;
		for (int i=0;i<imgRMap.length;i++) if (imgRMap[i]>=0)imgMap[imgRMap[i]]=i;
		double [][] jacobian=null;
		if (calcJacobian) {
//			jacobian=new double[fittedParNumber][indexCount*2];
//			for (int i=0;i<fittedParNumber;i++) for (int j=0;j<jacobian[0].length;j++) jacobian[i][j]=0.0;
			jacobian=new double[fittedParNumber][];
			// TODO: verify that only small number of rows is calculated
			for (int i=0;i<fittedParNumber;i++) {
				if (imgMap[i]>=0) {
					jacobian[i]=new double [indexCount*2];
					for (int j=0;j<jacobian[i].length;j++) jacobian[i][j]=0.0;
				} else {
					jacobian[i]=null;
				}
			}
		}
		double [][] interParameterDerivatives=new double [this.lensDistortionParameters.getNumInputs()][];
		//			initialize arrays for parameters and derivatives conversion
		double [] imgVector=fittingStrategy.getImageParametersVector(numImage, vector);     // thread safe
		if (this.debugLevel>3) {
			String all="imgVector: ";
			for (int jj=0;jj<imgVector.length;jj++) all+=" "+imgVector[jj];
			System.out.println(all);
		}
		// Calculate/set  this.lensDistortionParameters class, so it will calculate values/derivatives correctly)
		// and this.interParameterDerivatives
		if (this.debugLevel>3) System.out.println("calculatePartialFxAndJacobian(), numImage="+numImage+" calcInterParamers():");
		lensDistortionParameters.lensCalcInterParamers(
				lensDistortionParameters,
				this.fittingStrategy.distortionCalibrationData.isTripod(),
				this.fittingStrategy.distortionCalibrationData.isCartesian(),
	    		this.fittingStrategy.distortionCalibrationData.getPixelSize(numImage),
	    		this.fittingStrategy.distortionCalibrationData.getDistortionRadius(numImage),
				calcJacobian?interParameterDerivatives:null, // [22][]
						imgVector,
						imgMask); // calculate only selected derivatives (all parVect values are still

		// iterate through points, for each calculate pixelx, pixely and derivatives
		for (int pointNum=0;pointNum<indexCount;pointNum++){
			int fullIndex=indexStart+pointNum;
			double [][]derivatives15=  lensDistortionParameters.calcPartialDerivatives(
					patternXYZ[fullIndex][0], // target point horizontal, positive - right,  mm
					patternXYZ[fullIndex][1], // target point vertical,   positive - down,  mm
					patternXYZ[fullIndex][2], // target point horizontal, positive - away from camera,  mm
					calcJacobian); // calculate derivatives, false - values only
			if (this.debugLevel>3) {
				System.out.println(fullIndex+": calculateFxAndJacobian->calcPartialDerivatives("+IJ.d2s(patternXYZ[fullIndex][0],2)+","+
						IJ.d2s(patternXYZ[fullIndex][1],2)+","+
						IJ.d2s(patternXYZ[fullIndex][2],2)+" ("+calcJacobian+") -> "+
						IJ.d2s(derivatives15[0][0],2)+"/"+IJ.d2s(derivatives15[0][1],2));
				String all="derivatives15: D3476";
				for (int ii=0;ii<derivatives15.length;ii++) all+=" "+ii+":"+IJ.d2s(derivatives15[ii][0],3)+"/"+IJ.d2s(derivatives15[ii][1],3);
				System.out.println(all);
			}
			vectorFX[2*fullIndex]=  derivatives15[0][0];
			vectorFX[2*fullIndex+1]=derivatives15[0][1];
			if (calcJacobian) {
				double [][]derivatives = lensDistortionParameters.reorderPartialDerivatives(derivatives15);
				if (this.debugLevel>3) {
					String all="derivatives:";
					for (int ii=0;ii<derivatives.length;ii++) all+=" "+ii+":"+IJ.d2s(derivatives[ii][0],3)+"/"+IJ.d2s(derivatives[ii][1],3);
					System.out.println(all);
				}
				for (int i=0;i<jacobian.length;i++) if (imgMap[i]>=0){
					double sX=0,sY=0;
					for (int k=0;k<derivatives.length;k++){
						sX+=interParameterDerivatives[imgMap[i]][k]*derivatives[k][0];
						sY+=interParameterDerivatives[imgMap[i]][k]*derivatives[k][1];
					}
					jacobian[i][2*pointNum]=  sX;
					jacobian[i][2*pointNum+1]=sY;
				}
			}
		}
		return jacobian;
	}


	/**
	 *
	 * @param vector - parameter vector to be used
	 * @param imgNumber - number of image to process or -1 - use the first of selected in this strategy
	 * @return return Jacobian matrix for the selected image and individual parameters
	 *   NOTE: this one is not thread safe (used this.lensDistortionParameters)
	 */
	// used only to debug derivatives (delta==0 - real derivatives, delta>0 - difference)
	public double [][] calculateJacobian16(
			double [] vector,
			int imgNumber,
			double delta){ // these parameters can work for one image only
        int doubleNumAllPoints=this.Y.length; // all points in all images multiplied by 2 (x and y error are separate)
		double [][] jacobian16=new double[lensDistortionParameters.getNumOutputs()][doubleNumAllPoints];
		double []   values=    new double[doubleNumAllPoints];

		for (int i=0;i<jacobian16.length;i++) for (int j=0;j<doubleNumAllPoints;j++) jacobian16[i][j]=0.0;
		int numImg=fittingStrategy.distortionCalibrationData.getNumImages();
		boolean [] selectedImages=fittingStrategy.selectedImages();
		int index=0;
		for (int imgNum=0;imgNum<numImg;imgNum++) if (selectedImages[imgNum]) {
			int numPoints=fittingStrategy.distortionCalibrationData.getImageNumPoints(imgNum);
			if (imgNumber<0) imgNumber=imgNum; // -1 - use the first image in the list
			if (imgNum==imgNumber) {
				double [] imgVector=fittingStrategy.getImageParametersVector(imgNum, vector); //this.currentVector);
				boolean [] imgMask= fittingStrategy.getImageParametersVectorMask(imgNum);
				int []     imgRMap=  fittingStrategy.getImageParametersVectorReverseMap(imgNum);
				int []     imgMap=new int[vector.length];
				for (int i=0;i<imgMap.length;i++) imgMap[i]=-1;
				for (int i=0;i<imgRMap.length;i++) if (imgRMap[i]>=0)imgMap[imgRMap[i]]=i;
				// Calculate/set  this.lensDistortionParameters class, so it will calculate values/derivatives correctly)
				// and this.interParameterDerivatives
				if (this.debugLevel>2) {
					System.out.println("calculateJacobian15(), imgNum="+imgNum+" calcInterParamers():");
				}
				this.lensDistortionParameters.lensCalcInterParamers(
						this.lensDistortionParameters,
						this.fittingStrategy.distortionCalibrationData.isTripod(),
						this.fittingStrategy.distortionCalibrationData.isCartesian(),
			    		this.fittingStrategy.distortionCalibrationData.getPixelSize(imgNum),
			    		this.fittingStrategy.distortionCalibrationData.getDistortionRadius(imgNum),
						null, //this.interParameterDerivatives, // [22][]
						imgVector,
						imgMask); // calculate only selected derivatives (all parVect values are still
//						false); // probably can use false
				if (this.debugLevel>2) {
					System.out.println("calculateJacobian16(), numPoints="+numPoints+" (imgNum="+imgNum+")");
				}
				if (delta<=0) {
					// iterate through points, for each calculate pixelx, pixely and derivatives
					for (int pointNum=0;pointNum<numPoints;pointNum++){
						int fullIndex=index+pointNum;
						double [][]derivatives15=  lensDistortionParameters.calcPartialDerivatives(
								targetXYZ[fullIndex][0], // target point horizontal, positive - right,  mm
								targetXYZ[fullIndex][1], // target point vertical,   positive - down,  mm
								targetXYZ[fullIndex][2], // target point horizontal, positive - away from camera,  mm
								true); // calculate derivatives, false - values only
						if (this.debugLevel>3) {
							System.out.println(fullIndex+": calculateFxAndJacobian->calcPartialDerivatives("+IJ.d2s(targetXYZ[fullIndex][0],2)+","+
									IJ.d2s(targetXYZ[fullIndex][1],2)+","+
									IJ.d2s(targetXYZ[fullIndex][2],2)+" -> "+
									IJ.d2s(derivatives15[0][0],2)+"/"+IJ.d2s(derivatives15[0][1],2));
							String all="derivatives15: D3563";
							for (int ii=0;ii<derivatives15.length;ii++) all+=" "+ii+":"+IJ.d2s(derivatives15[ii][0],3)+"/"+IJ.d2s(derivatives15[ii][1],3);
							System.out.println(all);
						}
						double [][]derivatives = lensDistortionParameters.reorderPartialDerivativesAsNames(derivatives15);
						if (this.debugLevel>3) {
							String all="derivatives:";
							for (int ii=0;ii<derivatives.length;ii++) all+=" "+ii+":"+IJ.d2s(derivatives[ii][0],3)+"/"+IJ.d2s(derivatives[ii][1],3);
							System.out.println(all);
						}
						for (int i=0;i<derivatives.length;i++){
							jacobian16[i][2*fullIndex]= derivatives[i][0]; // oob 16
							jacobian16[i][2*fullIndex+1]=derivatives[i][1];
						}

					}
				} else {
					double [] parameterVector0=lensDistortionParameters.getAllVector();
					for (int pointNum=0;pointNum<numPoints;pointNum++){
						int fullIndex=index+pointNum;
						double [][]values2=  lensDistortionParameters.calcPartialDerivatives(
								targetXYZ[fullIndex][0], // target point horizontal, positive - right,  mm
								targetXYZ[fullIndex][1], // target point vertical,   positive - down,  mm
								targetXYZ[fullIndex][2], // target point horizontal, positive - away from camera,  mm
								false); // calculate derivatives, false - values only
						values[2*fullIndex]= values2[0][0];
						values[2*fullIndex+1]=values2[0][1];
					}
					for (int nPar=0;nPar<jacobian16.length;nPar++){
						double [] parameterVector=parameterVector0.clone();
						parameterVector[nPar]+=delta;
						lensDistortionParameters.setAllVector(parameterVector);
						for (int pointNum=0;pointNum<numPoints;pointNum++){
							int fullIndex=index+pointNum;
							double [][] values2=lensDistortionParameters.calcPartialDerivatives(
									targetXYZ[fullIndex][0], // target point horizontal, positive - right,  mm
									targetXYZ[fullIndex][1], // target point vertical,   positive - down,  mm
									targetXYZ[fullIndex][2], // target point horizontal, positive - away from camera,  mm
									false); // calculate derivatives, false - values only
							jacobian16[nPar][2*fullIndex]=  (values2[0][0]- values[2*fullIndex])/delta;
							jacobian16[nPar][2*fullIndex+1]=(values2[0][1]- values[2*fullIndex+1])/delta;
						}
					}
				}
				return jacobian16;
			}
			index+=numPoints;
		}
		return null; // should normally return from inside the for loop
	}


/*
List calibration
 */
    public boolean listImageParameters(boolean silent){
    	if (this.fittingStrategy==null) {
        		String msg="Fitting strategy does not exist, exiting";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
    	}
    	int numSeries=fittingStrategy.getNumSeries();
    	if (silent) {
    		this.seriesNumber=0;
    	} else {
    		GenericDialog gd = new GenericDialog("Settings for the parameter list");
    		gd.addNumericField("Iteration number to start (0.."+(numSeries-1)+")", this.seriesNumber, 0);
    		gd.addCheckbox("Show image number (from 0)",                           this.showIndex);
    		gd.addCheckbox("Show per-image RMS",                                   this.showRMS);
    		gd.addCheckbox("Show number of grid points",                           this.showPoints);
    		gd.addCheckbox("Show lens coordinates (relative to target)",           this.showLensLocation);

    		gd.addCheckbox("Show physical camera parameters",                      this.showEyesisParameters);
    		gd.addCheckbox("Show intrinsic lens/sensor parameters ",               this.showIntrinsicParameters);
    		gd.addCheckbox("Show extrinsic lens/sensor parameters",                this.showExtrinsicParameters);
    		gd.addNumericField("Extra decimal places (precision) in the list",     this.extraDecimals, 0);
    		gd.showDialog();
    		if (gd.wasCanceled()) return false;
    		this.seriesNumber=          (int) gd.getNextNumber();
    		this.showIndex=                   gd.getNextBoolean();
    		this.showRMS=                     gd.getNextBoolean();
    		this.showPoints=                  gd.getNextBoolean();
    		this.showLensLocation=            gd.getNextBoolean();
    		this.showEyesisParameters=        gd.getNextBoolean();
    		this.showIntrinsicParameters=     gd.getNextBoolean();
    		this.showExtrinsicParameters=     gd.getNextBoolean();
    		this.extraDecimals=         (int) gd.getNextNumber();
    	}
// need to select strategy
	    initFittingSeries(true,this.filterForAll,this.seriesNumber); // will set this.currentVector
		this.currentfX=calculateFxAndJacobian(this.currentVector, false); // is it always true here (this.jacobian==null)
		double [] errors=calcErrors(calcYminusFx(this.currentfX));
		double    rms=   calcError (calcYminusFx(this.currentfX));
		int [] numPairs=calcNumPairs();


		boolean [] selectedImages=fittingStrategy.selectedImages();
//TODO: add display of per-image RMS
	    listImageParameters (
	    		selectedImages,
	    		rms,
	    		errors,
	    		numPairs,
	    		this.showIndex,
	    		true, // grid match
	    		this.showRMS,
	    		this.showPoints,
	    		this.showLensLocation,
	    		this.showEyesisParameters,
	    		this.showIntrinsicParameters,
	    		this.showExtrinsicParameters,
	    		this.extraDecimals);
    	return true;
    }

    public void markBadNodces(int seriesNumber,
    		int debugLevel){
    	int oldSeries=this.seriesNumber;
    	this.seriesNumber=seriesNumber;
    	int totalBadNodes=markBadNodes(
    			fittingStrategy.distortionCalibrationData.eyesisCameraParameters.removeOverRMS,
    			fittingStrategy.distortionCalibrationData.eyesisCameraParameters.removeOverRMSNonweighted,
    			false,
    			debugLevel
    	);
    	if (debugLevel>0) {
    		System.out.println("Marked "+totalBadNodes+" nodes as bad (excessive errors, used fitting series #"+this.seriesNumber+")");
    	}
    	this.seriesNumber=oldSeries;
    }

    public boolean dialogMarkBadNodes(int debugLevel){
    	int numSeries=fittingStrategy.getNumSeries();
    	boolean verbose=false;
    	GenericDialog gd = new GenericDialog("Select parameters for marking bad nodes");
    	gd.addNumericField("Series number to use for selection (0.."+(numSeries-1)+")", this.seriesNumber, 0);
    	gd.addNumericField("Remove nodes with error greater than scaled RMS in that image, weighted",  fittingStrategy.distortionCalibrationData.eyesisCameraParameters.removeOverRMS, 2,6,"xRMS");
    	gd.addNumericField("Same, not weghted (not more permissive near the borders with low weight)", fittingStrategy.distortionCalibrationData.eyesisCameraParameters.removeOverRMSNonweighted, 2,6,"xRMS");
    	gd.addCheckbox    ("Verbose (report number of bad nodes per image)",verbose);
    	gd.showDialog();
    	if (gd.wasCanceled()) return false;
    	this.seriesNumber=          (int) gd.getNextNumber();
    	fittingStrategy.distortionCalibrationData.eyesisCameraParameters.removeOverRMS=            gd.getNextNumber();
    	fittingStrategy.distortionCalibrationData.eyesisCameraParameters.removeOverRMSNonweighted= gd.getNextNumber();
    	verbose=                                                                                   gd.getNextBoolean();
    	int totalBadNodes=markBadNodes(
    			fittingStrategy.distortionCalibrationData.eyesisCameraParameters.removeOverRMS,
    			fittingStrategy.distortionCalibrationData.eyesisCameraParameters.removeOverRMSNonweighted,
    			verbose,
    			debugLevel
    	);
    	if (debugLevel>0) {
    		System.out.println("Marked "+totalBadNodes+" nodes as bad (excessive errors");
    	}
    	return true;
    }

    public boolean removeOutLiers(
    		int series,
    		int numOutLiers,
    		boolean [] selectedChannels){
    	int numSeries=fittingStrategy.getNumSeries();
    	boolean removeEmpty=false;
    	boolean recalculate=false;
    	boolean applyChannelFilter=false;
		int filter=filterForAll;
    	if ((series<0) || (numOutLiers<0)) {
    		GenericDialog gd = new GenericDialog("Select series to process");
    		gd.addNumericField("Iteration number to start (0.."+(numSeries-1)+")", this.seriesNumber, 0);
    		if (selectedChannels != null) {
    			String s="";
    			for (boolean b:selectedChannels)s+=b?"+":"-";
    			gd.addCheckbox("Filter by channel selection ("+s+")", applyChannelFilter);
    		}
    		gd.addCheckbox("Recalculate parameters vector from selected strategy",recalculate);
    		gd.addNumericField("Number of outliers to show", 10, 0);
    		gd.addCheckbox("Remove empty (rms==NaN) images", removeEmpty);
    		gd.addCheckbox("Ask filter (current filter="+filter+")",    this.askFilter);
    		gd.showDialog();
    		if (gd.wasCanceled()) return false;
    		this.seriesNumber=                           (int) gd.getNextNumber();
    		if (selectedChannels != null) applyChannelFilter=  gd.getNextBoolean();
    		recalculate=                                       gd.getNextBoolean();
    		numOutLiers=                                (int) gd.getNextNumber();
    		removeEmpty=                                       gd.getNextBoolean();
    		this.askFilter=                                    gd.getNextBoolean();
    		if (this.askFilter) filter=  selectFilter(filter);
    		filter=0;
    	} else {
    		this.seriesNumber=series;
    	}
    	if (!applyChannelFilter) selectedChannels=null;
    	if (recalculate) {
    		resetGridImageMasks(); // FIXME: move elsewhere?
    	}
//	    initFittingSeries(!recalculate,this.filterForAll,this.seriesNumber); // will set this.currentVector
	    initFittingSeries(!recalculate,this.filterForAll,this.seriesNumber); // will set this.currentVector
		this.currentfX=calculateFxAndJacobian(this.currentVector, false); // is it always true here (this.jacobian==null)
		double [] errors=calcErrors(calcYminusFx(this.currentfX)); // seem to have errors? - now may return NaN!
		double    rms=   calcError (calcYminusFx(this.currentfX));
		boolean [] selectedImages=fittingStrategy.selectedImages();
		if (selectedChannels!=null){
			selectedImages=selectedImages.clone(); // disconnect from original for modification
			for (int i=0;i<selectedImages.length;i++) if (selectedImages[i]){
				int chn=this.fittingStrategy.distortionCalibrationData.gIP[i].channel;
				if ((chn<0) || (chn>=selectedChannels.length) || !selectedChannels[chn]){
					selectedImages[i]=false;
				}
			}
		}
		int numSelectedNotNaNImages=0;
		int numNaN=0;
		for (int i=0;i<selectedImages.length;i++) if (selectedImages[i]) {
			if (!Double.isNaN(errors[i])) numSelectedNotNaNImages++;
			else numNaN++;
		}
		int [] imgIndices=new int[numSelectedNotNaNImages];
		int index=0;
		for (int i=0;i<selectedImages.length;i++) if ( selectedImages[i] && !Double.isNaN(errors[i])) imgIndices[index++]=i; // OOB 2389

		if (numOutLiers>numSelectedNotNaNImages) numOutLiers=numSelectedNotNaNImages;
		int [] indices=new int [numOutLiers];
		boolean [] availableImages=selectedImages.clone();
		for (int i=0;i<selectedImages.length;i++) if (selectedImages[i] && Double.isNaN(errors[i])) availableImages[i]=false;


		if ((this.debugLevel>0) && (numNaN>0)){
			System.out.println("removeOutLiers(): Number of empty (rms=NaN) images="+numNaN+":");
			int n=0;
			for (int i=0;i<selectedImages.length;i++) if (selectedImages[i] && Double.isNaN(errors[i])){
				n++;
				System.out.println(n+": "+i+": "+this.fittingStrategy.distortionCalibrationData.gIP[i].path);
			}
		}
		if (removeEmpty){
			int n=0;
			for (int i=0;i<selectedImages.length;i++) if (selectedImages[i] && Double.isNaN(errors[i])){
				n++;
				if (this.debugLevel>0) System.out.println(n+"removing empty image #"+i+": "+this.fittingStrategy.distortionCalibrationData.gIP[i].path);
				this.fittingStrategy.distortionCalibrationData.gIP[i].enabled=false;
				this.fittingStrategy.distortionCalibrationData.gIP[i].hintedMatch=-1; // so can be re-calibrated again w/o others
			}

		}

		System.out.println("removeOutLiers(): availableImages.length="+availableImages.length+" numSelectedNotNaNImages="+numSelectedNotNaNImages);
		for (int n=0;n<numOutLiers;n++){
			double maxRMS=-1.0;
			indices[n]=-1;
			for (int i=0;i<availableImages.length;i++)if (availableImages[i] && (Double.isNaN(errors[i]) || (errors[i]>maxRMS))){ // Double.NaN will be greater
					maxRMS=errors[i];
					indices[n]=i;
			}
			if (indices[n]<0){
				System.out.println("removeOutLiers(): indices["+n+"]="+indices[n]);
				continue;
			}
			availableImages[indices[n]]=false; // java.lang.ArrayIndexOutOfBoundsException: -1
		}

		GenericDialog gd = new GenericDialog("Select images to remove (RMS="+IJ.d2s(rms,3)+")");
		if (this.debugLevel>0) System.out.println("Listing "+numOutLiers+" worst images:");
		for (int n=0;n<indices.length;n++){
			String msg=n+" ("+indices[n]+" / "+ this.fittingStrategy.distortionCalibrationData.gIP[indices[n]].getSetNumber()+"): "+
			IJ.d2s(errors[indices[n]],3)+" "+
			this.fittingStrategy.distortionCalibrationData.gIP[indices[n]].path+
			" ("+this.fittingStrategy.distortionCalibrationData.gIP[indices[n]].pixelsXY.length+
			" points) "+selectedImages[indices[n]];

			if (this.debugLevel>0) System.out.println(
					msg);
			gd.addCheckbox(msg, true);
		}
		WindowTools.addScrollBars(gd);
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		if (this.debugLevel>0) System.out.println("Removing outliers:");
		for (int n=0;n<indices.length;n++){
			if (gd.getNextBoolean()) {
				if (this.debugLevel>0) System.out.println(n+" :"+IJ.d2s(errors[indices[n]],3)+" "+this.fittingStrategy.distortionCalibrationData.gIP[indices[n]].path);
				this.fittingStrategy.distortionCalibrationData.gIP[indices[n]].enabled=false;
				this.fittingStrategy.distortionCalibrationData.gIP[indices[n]].hintedMatch=-1; // so can be re-calibrated again w/o others
			}
		}
		fittingStrategy.distortionCalibrationData.updateSetOrientation(null); // remove orientation information from the image set if none is enabled
		return true;
    }

    public boolean removeOutLierSets(int numOutLiers){
    	boolean removeEmptySets=true; // false;
    	if (numOutLiers<0) {
    		GenericDialog gd = new GenericDialog("Select sets to process");
    		gd.addNumericField("Series number (<0 - all images)", -1, 0);
    		gd.addNumericField("Number of outliers to show", 5, 0);
    		gd.addCheckbox("Remove empty sets", removeEmptySets);
    		gd.addCheckbox("Ask for weight function filter",     this.askFilter);
    		gd.showDialog();
    		if (gd.wasCanceled()) return false;
    		this.seriesNumber= (int) gd.getNextNumber();
    		numOutLiers=      (int) gd.getNextNumber();
    		removeEmptySets=         gd.getNextBoolean();
    		this.askFilter=         gd.getNextBoolean();
    	}
//		boolean [] oldSelection=this.fittingStrategy.selectAllImages(0); // enable all images in series 0
		int filter=this.filterForAll;
		if (this.askFilter) filter=selectFilter(filter);

    	initFittingSeries(true,filter, -1); // this.seriesNumber); // will set this.currentVector
		this.currentfX=calculateFxAndJacobian(this.currentVector, false); // is it always true here (this.jacobian==null)
		double [] errors=calcErrors(calcYminusFx(this.currentfX)); // for all images
		double    rms=   calcError (calcYminusFx(this.currentfX));
		int []    numPairs=calcNumPairs(); // for all images, not only selected

// re-init for the selected series
    	initFittingSeries(true,filter,this.seriesNumber); // will set this.currentVector
//    	initFittingSeries(true,this.filterForAll,this.seriesNumber); // will set this.currentVector
		this.currentfX=calculateFxAndJacobian(this.currentVector, false); // is it always true here (this.jacobian==null)
//		double [] errors=calcErrors(calcYminusFx(this.currentfX)); // error - always for -1?
		rms=   calcError (calcYminusFx(this.currentfX)); // for selected series



		boolean [] selectedImages=fittingStrategy.selectedImages(this.seriesNumber);
//		int [] numPairs=calcNumPairs();
		int [][] imageSets=this.fittingStrategy.distortionCalibrationData.listImages(
				false, // true - only enabled images
				null);    // do not filter eo, lwir
    
	    int [] numSetPoints=new int [imageSets.length];
	    double [] rmsPerSet=new double[imageSets.length];
	    boolean [] hasNaNInSet=new boolean[imageSets.length];
	    boolean [] allNaNInSet=new boolean[imageSets.length];
	    for (int setNum=0;setNum<imageSets.length;setNum++){
	    	double error2=0.0;
	    	int numInSet=0;
	    	int numInSetOther=0; // not selected in this
    		hasNaNInSet[setNum]=false;
    		boolean has_selected = false; // remove selected with all nan in selected, and unselected with all nans
    		for (int imgInSet=0;imgInSet<imageSets[setNum].length;imgInSet++){
    			int imgNum=imageSets[setNum][imgInSet];
				int num=numPairs[imgNum];
    			if (selectedImages[imgNum]) {
    				has_selected = true;
    				if (Double.isNaN(errors[imgNum])){
    					hasNaNInSet[setNum]=true;
    				} else {
    					error2+=errors[imgNum]*errors[imgNum]*num;
    					numInSet+=num;

    				}
    			} else {
    				if (!Double.isNaN(errors[imgNum])){
    					numInSetOther += num;
    				}
    			}
	    	}
//    		allNaNInSet[setNum]= hasNaNInSet[setNum] && (numInSet==0);
    		allNaNInSet[setNum]= has_selected? (numInSet == 0) : ((numInSet + numInSetOther) == 0);
	    	numSetPoints[setNum]=numInSet;
	    	rmsPerSet[setNum]=(numInSet>0)?Math.sqrt(error2/numInSet) : Double.NaN; // only count selected images (i.e. only eo, ignore lwir)
	    }
//		int numSelectedNotNaNSets=0;
		int numSelectedSets=0;
		int numNaN=0;
		for (int i=0;i<imageSets.length;i++)  {
//			if (!Double.isNaN(rmsPerSet[i])) numSelectedSets++;
			if (!allNaNInSet[i]) numSelectedSets++;
			else numNaN++;
		}
//		int [] imgIndices=new int[numSelectedNotNaNSets];
//		int index=0;
//		for (int i=0;i<imageSets.length;i++) if ( selectedImages[i]) imgIndices[index++]=i;

		if (numOutLiers>numSelectedSets) numOutLiers=numSelectedSets;
		int [] indices=new int [numOutLiers];
		boolean [] availableSets= new boolean  [imageSets.length];
		for (int i=0;i<imageSets.length;i++) availableSets[i]= !allNaNInSet[i]; //!Double.isNaN(rmsPerSet[i]);
/*
		// Remove all empty, not just selected by strategy. Now errors are calculated for all images, not juet selected
		if (removeEmptySets  && (numNaN>0)){ //(this.debugLevel>0)
			if (this.debugLevel>-1) System.out.println("removeOutLierSets(): Number of empty (rms=NaN) sets="+numNaN+":");
//			int n=0;
			for (int setNum=0;setNum<imageSets.length;setNum++) if (!availableSets[setNum]){
//				n++;
				if (this.debugLevel>0) System.out.println("Set "+setNum);
	    		for (int imgInSet=0;imgInSet<imageSets[setNum].length;imgInSet++){
					int numImg=imageSets[setNum][imgInSet];
					if (this.debugLevel>0) System.out.println(setNum+":"+imgInSet+" #"+ numImg+" "+IJ.d2s(errors[numImg],3)+" "+
							this.fittingStrategy.distortionCalibrationData.gIP[numImg].path);
					this.fittingStrategy.distortionCalibrationData.gIP[numImg].enabled=false;
					this.fittingStrategy.distortionCalibrationData.gIP[numImg].hintedMatch=-1; // so can be re-calibrated again w/o others
	    		}
			}
		}
*/
		System.out.println("removeOutLierSets(): availableSets.length="+availableSets.length+" numSelectedSets="+numSelectedSets);
		for (int n=0;n<numOutLiers;n++){
			double maxRMS=-1.0;
			indices[n]=-1;
			for (int i=0;i<availableSets.length;i++)if (availableSets[i] && (rmsPerSet[i]>maxRMS)){ // NaN are already skipped
					maxRMS=rmsPerSet[i];
					indices[n]=i;
			}
			if (indices[n]<0){
				System.out.println("removeOutLierSets(): indices["+n+"]="+indices[n]);
				continue;
			}
			availableSets[indices[n]]=false; // java.lang.ArrayIndexOutOfBoundsException: -1
		}

		GenericDialog gd = new GenericDialog("Select image Sets to remove (RMS="+IJ.d2s(rms,3)+")");
		if (this.debugLevel>0) System.out.println("Listing "+numOutLiers+" worst image sets");
		for (int n=0;n<indices.length;n++){
			int numSet=indices[n];
			if (numSet >= 0) {
				double setWeight=this.fittingStrategy.distortionCalibrationData.gIS[numSet].setWeight; //-1
				if (this.debugLevel>0) System.out.println(n+" ("+numSet+"): "+(hasNaNInSet[numSet]?"* ":"")+IJ.d2s(rmsPerSet[numSet],3)+
						" points: "+numSetPoints[numSet]+" weight:"+setWeight);
				gd.addCheckbox(n+": "+numSet+": "+(hasNaNInSet[numSet]?"* ":"")+IJ.d2s(rmsPerSet[numSet],3)+" weight:"+setWeight, true);
				for (int i=0;i<imageSets[numSet].length;i++){
					int numImg=imageSets[numSet][i];
					double diameter=this.fittingStrategy.distortionCalibrationData.gIP[numImg].getGridDiameter();
					gd.addMessage(i+":"+numImg+": "+IJ.d2s(errors[numImg],3)+" "+
							" ("+this.fittingStrategy.distortionCalibrationData.gIP[numImg].pixelsXY.length+" points, diameter="+diameter+") "+
							this.fittingStrategy.distortionCalibrationData.gIP[numImg].path);
					if (this.debugLevel>0) System.out.println("  --- "+numImg+": "+IJ.d2s(errors[numImg],3)+" "+
							" ("+this.fittingStrategy.distortionCalibrationData.gIP[numImg].pixelsXY.length+" points, diameter="+diameter+") "+
							this.fittingStrategy.distortionCalibrationData.gIP[numImg].path);
				}
			}
		}
		WindowTools.addScrollBars(gd);
		gd.showDialog();
		if (gd.wasCanceled()){
//			this.fittingStrategy.setImageSelection(0, oldSelection); // restore original selection in series 0
			return false;
		}
		if (this.debugLevel>0) System.out.println("Removing outliers:");
		for (int n=0;n<indices.length;n++){
			int numSet=indices[n];
			if (numSet >= 0) {
				if (gd.getNextBoolean()) {
					if (this.debugLevel>0) System.out.println(" Removing imgages in image set "+numSet);
					for (int i=0;i<imageSets[numSet].length;i++){
						int numImg=imageSets[numSet][i];
						if (this.debugLevel>0) System.out.println(n+":"+i+"("+numImg+")"+IJ.d2s(errors[numImg],3)+" "+
								this.fittingStrategy.distortionCalibrationData.gIP[numImg].path);
						this.fittingStrategy.distortionCalibrationData.gIP[numImg].enabled=false;
						this.fittingStrategy.distortionCalibrationData.gIP[numImg].hintedMatch=-1; // so can be re-calibrated again w/o others
					}
				}
			}
		}

		// Remove all empty, not just selected by strategy. Now errors are calculated for all images, not juet selected
		if (removeEmptySets  && (numNaN>0)){ //(this.debugLevel>0)
			if (this.debugLevel>-1) System.out.println("removeOutLierSets(): Number of empty (rms=NaN) sets="+numNaN+":");
//			int n=0;
			for (int setNum=0;setNum<imageSets.length;setNum++) if (!availableSets[setNum]){
//				n++;
				if (this.debugLevel>0) System.out.println("Set "+setNum);
	    		for (int imgInSet=0;imgInSet<imageSets[setNum].length;imgInSet++){
					int numImg=imageSets[setNum][imgInSet];
					if (this.debugLevel>0) System.out.println(setNum+":"+imgInSet+" #"+ numImg+" "+IJ.d2s(errors[numImg],3)+" "+
							this.fittingStrategy.distortionCalibrationData.gIP[numImg].path);
					this.fittingStrategy.distortionCalibrationData.gIP[numImg].enabled=false;
					this.fittingStrategy.distortionCalibrationData.gIP[numImg].hintedMatch=-1; // so can be re-calibrated again w/o others
	    		}
			}
		}


		// next is not needed
		fittingStrategy.distortionCalibrationData.updateSetOrientation(null); //selectedImages); // null); // remove orientation information from the image set if none is enabled
//		this.fittingStrategy.setImageSelection(0, oldSelection); // restore original selection in series 0
		return true;
    }

    public boolean removeOutLiersJunk(int series, int numOutLiers){
    	int numSeries=fittingStrategy.getNumSeries();
    	if ((series<0) || (numOutLiers<0)) {
    		GenericDialog gd = new GenericDialog("Select series to process");
    		gd.addNumericField("Iteration number to start (0.."+(numSeries-1)+")", this.seriesNumber, 0);
    		gd.addNumericField("Number of outliers to show", 10, 0);
    		gd.showDialog();
    		if (gd.wasCanceled()) return false;
    		this.seriesNumber=          (int) gd.getNextNumber();
    		numOutLiers=               (int) gd.getNextNumber();
    	} else {
    		this.seriesNumber=series;
    	}
	    initFittingSeries(true,this.filterForAll,this.seriesNumber); // will set this.currentVector
		this.currentfX=calculateFxAndJacobian(this.currentVector, false); // is it always true here (this.jacobian==null)
		double [] errors=calcErrors(calcYminusFx(this.currentfX));
		double    rms=   calcError (calcYminusFx(this.currentfX));
		boolean [] selectedImages=fittingStrategy.selectedImages();
		int numSelectedImages=0;
		for (int i=0;i<selectedImages.length;i++) if (selectedImages[i]) numSelectedImages++;
		int [] imgIndices=new int[numSelectedImages];
		int index=0;
		for (int i=0;i<selectedImages.length;i++) if ( selectedImages[i]) imgIndices[index++]=i;

		if (numOutLiers>numSelectedImages) numOutLiers=numSelectedImages;
		int [] indices=new int [numOutLiers];
		int [] indicesSelected=new int [numOutLiers];
		boolean [] availableImages=new boolean[numSelectedImages];
		for (int i=0;i<availableImages.length;i++)availableImages[i]=true;
		for (int n=0;n<numOutLiers;n++){
			double maxRMS=0;
			indices[n]=-1;
			indicesSelected[n]=-1;
			int imgIndex=0;
			for (int i=0;i<selectedImages.length;i++)if (selectedImages[i]){
				if (availableImages[imgIndex] && (errors[imgIndex]>maxRMS)){
					maxRMS=errors[imgIndex];
					indicesSelected[n]=imgIndex;
					indices[n]=i;
				}
				imgIndex++;
			}
			availableImages[indicesSelected[n]]=false;
		}
		GenericDialog gd = new GenericDialog("Select images to remove (RMS="+IJ.d2s(rms,3)+")");
		for (int n=0;n<indices.length;n++){
			gd.addCheckbox(indices[n]+": "+IJ.d2s(errors[indicesSelected[n]],3)+" "+this.fittingStrategy.distortionCalibrationData.gIP[indices[n]].path, true);
		}
		WindowTools.addScrollBars(gd);
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		for (int n=0;n<indices.length;n++){
			if (gd.getNextBoolean()) this.fittingStrategy.distortionCalibrationData.gIP[indices[n]].enabled=false;
		}
		return true;
    }

	/**
	 * Opens a text window with the parameter table
	 * @param imageSelection which images include in the output
	 * @param showEyesisParameters show physical location/attitude based on Eyesis
	 * @param showIntrinsicParameters show lens distortion/alignment parameters)
	 * @param showExtrinsicParameters show position/attitude of the individual cameras
	 * @param extraDecimals add this many decimals to data
	 */
	public void listImageParameters (boolean [] imageSelection,
    		double rms,
    		double [] errors,
    		int    [] numPairs,
    		boolean showIndex,
    		boolean showGridMatch,
    		boolean showErrors,
    	    boolean showPoints,
    		boolean showLensCoordinates,
			boolean showEyesisParameters,
			boolean showIntrinsicParameters,
			boolean showExtrinsicParameters,
			int extraDecimals){
		int numImages=0;

		for (int i=0;i<fittingStrategy.distortionCalibrationData.getNumImages();i++) {
			if((imageSelection==null) || ((i<imageSelection.length) && imageSelection[i])) numImages++;
		}
		double  [][] intrinsic=new double [numImages][];
		double  [][] extrinsic=new double [numImages][];
		double  [][] lensCoordinates=new double [numImages][];
		int [] imgIndices=new int[numImages];
		int index=0;
		for (int i=0;i<fittingStrategy.distortionCalibrationData.getNumImages();i++) {
			if((imageSelection==null) || ((i<imageSelection.length) && imageSelection[i])) imgIndices[index++]=i;
		}
		for (int imgIndex=0;imgIndex<numImages;imgIndex++) {
			int imgNum=imgIndices[imgIndex]; // image number
			if (this.debugLevel>2) {
				System.out.println("listImageParameters(), imgNum="+imgNum+" calcInterParamers():");
			}
			this.lensDistortionParameters.lensCalcInterParamers(
					this.lensDistortionParameters,
					this.fittingStrategy.distortionCalibrationData.isTripod(),
					this.fittingStrategy.distortionCalibrationData.isCartesian(),
		    		this.fittingStrategy.distortionCalibrationData.getPixelSize(imgNum),
		    		this.fittingStrategy.distortionCalibrationData.getDistortionRadius(imgNum),
					null, //this.interParameterDerivatives, // [22][]
//					fittingStrategy.distortionCalibrationData.pars[imgNum], // 22-long parameter vector for the image
					fittingStrategy.distortionCalibrationData.getParameters(imgNum), // 22-long parameter vector for the image
					null); // if no derivatives, null is OK
//					false); // calculate this.interParameterDerivatives -derivatives array (false - just this.values)
			intrinsic[imgIndex]=      lensDistortionParameters.getIntrinsicVector().clone();
			extrinsic[imgIndex]=      lensDistortionParameters.getExtrinsicVector().clone();
			lensCoordinates[imgIndex]=lensDistortionParameters.getLensCenterCoordinates();
		}

	    String header="Name\tUnits";
		for (int imgIndex=0;imgIndex<numImages;imgIndex++)
			header+="\t"+IJ.d2s(fittingStrategy.distortionCalibrationData.getImageTimestamp(imgIndices[imgIndex]),6);
	    StringBuffer sb = new StringBuffer();
	    if (showIndex) {
			sb.append("Station \t");
			for (int imgIndex=0;imgIndex<numImages;imgIndex++){
				sb.append("\t"+fittingStrategy.distortionCalibrationData.gIP[imgIndices[imgIndex]].getStationNumber());
			}
			sb.append("\n");
			sb.append("Set \t");
			for (int imgIndex=0;imgIndex<numImages;imgIndex++){
				sb.append("\t"+fittingStrategy.distortionCalibrationData.gIP[imgIndices[imgIndex]].getSetNumber());
			}
			sb.append("\n");
			sb.append("Index \t");
			for (int imgIndex=0;imgIndex<numImages;imgIndex++){
				sb.append("\t"+imgIndices[imgIndex]);
			}

	    }
	    if (showGridMatch){
			sb.append("Grid Match"+"\tX/Y:ROT");
			for (int imgIndex=0;imgIndex<numImages;imgIndex++){
				int imgNum=imgIndices[imgIndex]; // image number
				int [] shiftRot=fittingStrategy.distortionCalibrationData.getUVShiftRot(imgNum);
				sb.append("\t"+shiftRot[0]+"/"+shiftRot[1]+":"+shiftRot[2]);
			}
			sb.append("\n");

			sb.append("Lasers(matched)"+"\t");
			for (int imgIndex=0;imgIndex<numImages;imgIndex++){
				int imgNum=imgIndices[imgIndex]; // image number

    			int numPointers=0; // count number of laser pointers
    	        DistortionCalibrationData.GridImageParameters gip=fittingStrategy.distortionCalibrationData.getGridImageParameters(imgNum);
    			if (gip.laserPixelCoordinates!=null){
    				for (int j=0;j<gip.laserPixelCoordinates.length;j++) if (gip.laserPixelCoordinates[j]!=null) numPointers++;
    			}
    			sb.append("\t");
    			if (!gip.enabled) sb.append("(");
    			sb.append(numPointers+"("+gip.matchedPointers+"):"+gip.hintedMatch +
    					" "+IJ.d2s(gip.getGridPeriod(),1));
    			if (!gip.enabled) sb.append(")");
			}
			sb.append("\n");
	    }
		if (showErrors) {
				sb.append("--- RMS "+IJ.d2s(rms,3+extraDecimals)+"\tpix");
				for (int imgIndex=0;imgIndex<numImages;imgIndex++){
					int imgNum=imgIndices[imgIndex]; // image number
					sb.append("\t"+IJ.d2s(errors[imgNum],3+extraDecimals));
				}
				sb.append("\n");
		}
		if (showPoints) {
			int totalPoints=0;
			for (int imgIndex=0;imgIndex<numImages;imgIndex++){
				totalPoints+=numPairs[imgIndices[imgIndex]];
			}
			sb.append(" points "+totalPoints+"\t");
			for (int imgIndex=0;imgIndex<numImages;imgIndex++){
				sb.append("\t"+numPairs[imgIndices[imgIndex]]);
			}
			sb.append("\n");
			sb.append(" Diameter\trel");
			for (int imgIndex=0;imgIndex<numImages;imgIndex++){
				int imgNum=imgIndices[imgIndex]; // image number
				sb.append("\t"+IJ.d2s(this.fittingStrategy.distortionCalibrationData.gIP[imgNum].getGridDiameter(),2));
			}
			sb.append("\n");
		}
		if (showEyesisParameters) {
//			getImageSubcamera
			sb.append("Sub-camera\t");
			for (int imgIndex=0;imgIndex<numImages;imgIndex++){
				int imgNum=imgIndices[imgIndex]; // image number
				sb.append("\t"+fittingStrategy.distortionCalibrationData.getImageSubcamera(imgNum));
			}
			sb.append("\n");

			for (int parNumber=0;parNumber<fittingStrategy.distortionCalibrationData.getNumDescriptions();parNumber++){
				sb.append(
						fittingStrategy.distortionCalibrationData.descrField(parNumber,0)+"\t"+
						fittingStrategy.distortionCalibrationData.descrField(parNumber,2));
				for (int imgIndex=0;imgIndex<numImages;imgIndex++){
					int imgNum=imgIndices[imgIndex]; // image number
//					sb.append("\t"+IJ.d2s(fittingStrategy.distortionCalibrationData.pars[imgNum][parNumber],3+extraDecimals)); // TODO: make an array of decimals per parameter
					sb.append("\t"+IJ.d2s(fittingStrategy.distortionCalibrationData.getParameterValue(imgNum,parNumber),3+extraDecimals)); // TODO: make an array of decimals per parameter
				}
				sb.append("\n");
			}
		}
		if (showIntrinsicParameters) {
			sb.append("--- Intrinsic\t");
			for (int imgIndex=0;imgIndex<numImages;imgIndex++) sb.append("\t---");
			sb.append("\n");
			for (int parNumber=0;parNumber<lensDistortionParameters.getIntrinsicNames().length;parNumber++){
				sb.append(
						lensDistortionParameters.getIntrinsicNames()[parNumber]+"\t"+
						lensDistortionParameters.getIntrinsicUnits()[parNumber]);
				for (int imgIndex=0;imgIndex<numImages;imgIndex++){
					sb.append("\t"+IJ.d2s(intrinsic[imgIndex][parNumber],3+extraDecimals)); // TODO: make an array of decimals per parameter
				}
				sb.append("\n");
			}
		}
		if (showExtrinsicParameters ||  showLensCoordinates) {
			sb.append("--- Extrinsic\t");
			for (int imgIndex=0;imgIndex<numImages;imgIndex++) sb.append("\t---");
			sb.append("\n");
	        if (showLensCoordinates){
				sb.append("Lens X(right)\tmm");
				for (int imgIndex=0;imgIndex<numImages;imgIndex++){
//					int imgNum=imgIndices[imgIndex]; // image number
					sb.append("\t"+IJ.d2s(lensCoordinates[imgIndex][0],3+extraDecimals));
				}
				sb.append("\n");
				sb.append("Lens Y(down)\tmm");
				for (int imgIndex=0;imgIndex<numImages;imgIndex++){
//					int imgNum=imgIndices[imgIndex]; // image number
					sb.append("\t"+IJ.d2s(lensCoordinates[imgIndex][1],3+extraDecimals));
				}
				sb.append("\n");
				sb.append("Lens Z(into)\tmm");
				for (int imgIndex=0;imgIndex<numImages;imgIndex++){
//					int imgNum=imgIndices[imgIndex]; // image number
					sb.append("\t"+IJ.d2s(lensCoordinates[imgIndex][2],3+extraDecimals));
				}
				sb.append("\n");
	        }
			if (showExtrinsicParameters){
				for (int parNumber=0;parNumber<lensDistortionParameters.getExtrinsicNames().length;parNumber++){
					sb.append(
							lensDistortionParameters.getExtrinsicNames()[parNumber]+"\t"+
							lensDistortionParameters.getExtrinsicUnits()[parNumber]);
					for (int imgIndex=0;imgIndex<numImages;imgIndex++){
						sb.append("\t"+IJ.d2s(extrinsic[imgIndex][parNumber],3+extraDecimals)); // TODO: make an array of decimals per parameter
					}
					sb.append("\n");
				}
			}
		}
	    new TextWindow("Camera/lens parameters", header, sb.toString(), 500,900);
	}

	/**
	 * Calculate differences vector
	 * @param fX vector of calculated pixelX,pixelY on the sensors
	 * @return same dimension vector of differences from this.Y (measured grid pixelxX, pixelY)
	 */
	public double [] calcYminusFx(double [] fX){
		double [] result=this.Y.clone();
		for (int i=0;i<result.length;i++) result[i]-=fX[i];
	    return result;
	}
	/**
	 * Calcualte partial differences vector
	 * @param fX vector of reprojected pixelX,pixelY on the sensors (number of elements - double number of points
	 * @param startIndex start index to extract (even number, twice point index)
	 * @param endIndex end index (1 greater than the last to extract)
	 * @return partial differences (measured/corrected -reprojected), twice number of points long
	 */
	public double [] calcYminusFx(double [] fX, int startIndex, int endIndex){
		double [] result=new double [endIndex-startIndex];
		for (int i=0;i<result.length;i++) {
			int index=startIndex+i;
			result[i]=this.Y[index]-fX[index];
		}
		return result;
	}


	/**
	 * Calculate the RMS from the differences vector
	 * @param diff - differences vector
	 * @return RMS for the mean error (in sensor pixels)
	 */
	public double calcError(double [] diff){
		double result=0;
		double sumw = 0;
		if (this.weightFunction!=null) {
			for (int i=0;i<diff.length;i++) if (!Double.isNaN(diff[i])){
				result+=diff[i]*diff[i]*this.weightFunction[i];
				sumw += this.weightFunction[i];
			}
			result/=sumw; // this.sumWeights;
		} else {
			for (int i=0;i<diff.length;i++)  if (!Double.isNaN(diff[i])){
				result+=diff[i]*diff[i];
				sumw += 1.0;
			}
			result/=sumw; // diff.length;
		}
		return Math.sqrt(result)*this.RMSscale;
	}


	public double calcErrorDiffY(double [] fX){
		double result=0;
		double dbg_maxdiff = 0.0+0.0;
		if (this.weightFunction!=null) {
			for (int i=0;i<fX.length;i++)if (this.weightFunction[i] != 0.0){
				double diff=this.Y[i]-fX[i];
				if ((this.debugLevel>1) &&(Math.abs(diff) > dbg_maxdiff)) {
					System.out.println("i="+i+", diff="+diff+" y="+this.Y[i]+" fX="+fX[i]+ " w="+this.weightFunction[i]);
					dbg_maxdiff = Math.abs(diff);
				}
				result+=diff*diff*this.weightFunction[i];
			}
			result/=this.sumWeights;
		} else {
			for (int i=0;i<fX.length;i++){
				double diff=this.Y[i]-fX[i];
//				if (Math.abs(diff) > 1.0) {
//					System.out.print("");
//				}
				result+=diff*diff;
			}
			result/=fX.length;
		}
		return Math.sqrt(result)*this.RMSscale;
	}
	public double calcErrorDiffY(
			double [] fX,
			double [] extraWeightedErrors,
			double [] extraWeights){
		double result=0;
		double effectiveWeight;
		if (this.weightFunction!=null) {
			effectiveWeight=this.sumWeights;
			for (int i=0;i<fX.length;i++){
				double diff=this.Y[i]-fX[i];
//				if (Math.abs(diff) > 1.0) {
//					System.out.print("");
//				}
				result+=diff*diff*this.weightFunction[i];
			}
		} else {
			effectiveWeight=fX.length;
			for (int i=0;i<fX.length;i++){
				double diff=this.Y[i]-fX[i];
//				if (Math.abs(diff) > 1.0) {
//					System.out.print("");
//				}
				result+=diff*diff;
			}
		}
		if ((extraWeightedErrors!=null) && (extraWeights!=null)) {
			for (int i=0;i<extraWeightedErrors.length;i++){
				result+=extraWeightedErrors[i];
				effectiveWeight+=extraWeights[i];
			}
		}
		result/=effectiveWeight;
		return Math.sqrt(result)*this.RMSscale;
	}

	public void resetBadNodes(){
		for (int imgNum=0;imgNum<fittingStrategy.distortionCalibrationData.gIP.length;imgNum++) if (fittingStrategy.distortionCalibrationData.gIP[imgNum]!=null){
			fittingStrategy.distortionCalibrationData.gIP[imgNum].resetBadNodes();
		}
	}

    public int markBadNodes(
//    		int numSeries,
    		double removeOverRMS,
    		double removeOverRMSNonweighted,
    		boolean verbose,
    		int debugLevel){
    	int debugThreshold=2;
//		this.seriesNumber=series;
    	resetBadNodes(); // before calculating weight function
	    initFittingSeries(false,this.filterForAll,this.seriesNumber); // recalculate
		this.currentfX=calculateFxAndJacobian(this.currentVector, false); // is it always true here (this.jacobian==null)
    	int totalBadNodes=0;
		boolean [] selectedImages=fittingStrategy.selectedImages();
		double [] diff=calcYminusFx(this.currentfX);
		int index=0;
		for (int imgNum=0;imgNum<selectedImages.length;imgNum++) {
			double e2 =0.0;
//			errors[imgNum]=0.0;
			if (selectedImages[imgNum]) {
				int numThisRemoved=0;
				int len=fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV.length*2;
				double w=0;
				int nw=0;
				if (this.weightFunction!=null) {
					for (int i=index;i<index+len;i++) {
						e2+=diff[i]*diff[i]*this.weightFunction[i];
						w+=this.weightFunction[i];
						if (this.weightFunction[i]>0.0) nw++;
					}
				} else {
					for (int i=index;i<index+len;i++) {
						e2+=diff[i]*diff[i];
						w+=1.0;
						nw++;
					}
				}
				if (w>0.0) {
//					e2/=w;
					double threshold2Weighted=   2.0*removeOverRMS*removeOverRMS*e2/nw; // 2.0 because x^2+y^2
					double threshold2NonWeighted=2.0*removeOverRMSNonweighted*removeOverRMSNonweighted*e2/w; // 2.0 because x^2+y^2
					if (debugLevel>debugThreshold){
						boolean someRemoved=false;
						for (int i=index;i<index+len;i+=2) {
							e2=(diff[i]*diff[i]+ diff[i+1]*diff[i+1]);
							if ((e2>threshold2NonWeighted) || ((this.weightFunction!=null) && ((e2*this.weightFunction[i]) > threshold2Weighted )) ) {
								double ww=(this.weightFunction==null)?1.0:(this.weightFunction[i]);
								if (ww>0.0) someRemoved=true;

							}
						}
						if (someRemoved || (debugLevel>2)) System.out.println("imgNum="+imgNum+" len="+len+" e2/w="+(e2/w)+" w="+w+" e2/nw="+(e2/nw)+
								" threshold2Weighted="+threshold2Weighted+" threshold2NonWeighted="+threshold2NonWeighted);
					}

					for (int i=index;i<index+len;i+=2) {
						e2=(diff[i]*diff[i]+ diff[i+1]*diff[i+1]);
						if ((e2>threshold2NonWeighted) || ((this.weightFunction!=null) && ((e2*this.weightFunction[i]) > threshold2Weighted )) ) {
							double ww=(this.weightFunction==null)?1.0:(this.weightFunction[i]);
							int pointIndex=(i-index)/2;
							if (ww>0.0) {
								fittingStrategy.distortionCalibrationData.gIP[imgNum].setBadNode(pointIndex);
								numThisRemoved++;
								if (debugLevel>debugThreshold){
									int iu=fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV[pointIndex][0];
									int iv=fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV[pointIndex][1];
									System.out.println(numThisRemoved+": "+pointIndex +
											" uv="+iu+":"+iv+
											" e2="+e2+" ww="+ww+" e2w="+
											(e2*ww)+" ["+diff[i]+":"+diff[i+1]+"]");
								}
							}
						}
					}
					if (verbose && (numThisRemoved>0)) {
						System.out.println("Image "+imgNum+": removed "+numThisRemoved+" nodes over threshold");
					}
					totalBadNodes+=numThisRemoved;
				}
				index+=len;
			}
		}
    	return totalBadNodes;
    }
	public boolean showImageReprojectionErrorsDialog( int debugLevel){
		boolean eachImageInSet=false;
		boolean showGrids = false;
	    GenericDialog gd = new GenericDialog("Show Reprojection errors for image/image set/image selection");
		gd.addNumericField("Series number for image selection (-1 - all enabled images)", -1, 0);
		gd.addNumericField("Single image number to show (<0 - do not select)", -1,0);
		gd.addNumericField("Image set number to show (<0 - do not select)", -1,0);
		gd.addCheckbox("Open each image in the set",     eachImageInSet);
		gd.addCheckbox("Ask for weight function filter",     this.askFilter);
		gd.addCheckbox("Show grid images",                   showGrids);
//		gd.addNumericField("Weight function filter (-1 - use default for all )",-1,0);
	    gd.showDialog();
	    if (gd.wasCanceled()) return false;
	    this.seriesNumber=        (int) gd.getNextNumber();
	    int singleImageNumber=    (int) gd.getNextNumber();
	    int imageSetNumber=       (int) gd.getNextNumber();
	    eachImageInSet=                 gd.getNextBoolean();
	    this.askFilter=                 gd.getNextBoolean();
//	    int weightFunctionFilter= (int) gd.getNextNumber();
	    showGrids =                  gd.getNextBoolean();
		int filter=this.filterForAll;
		if (this.askFilter) filter=selectFilter(filter);
	    int [] imageNumbers = null;
	    if (singleImageNumber>=0){
	    	imageNumbers=new int [1];
	    	imageNumbers[0]=singleImageNumber;
	    	if (showGrids) {
	    		this.fittingStrategy.distortionCalibrationData.gIP[singleImageNumber].showGridImage();
	    	}
	    } else if (imageSetNumber>=0){
	    	int numInSet=0;
	    	for (int nChn=0;nChn<this.fittingStrategy.distortionCalibrationData.gIS[imageSetNumber].imageSet.length;nChn++){
	    		if (this.fittingStrategy.distortionCalibrationData.gIS[imageSetNumber].imageSet[nChn]!=null) numInSet++;
	    	}
	    	imageNumbers=new int [numInSet];
	    	numInSet=0;
	    	for (int nChn=0;nChn<this.fittingStrategy.distortionCalibrationData.gIS[imageSetNumber].imageSet.length;nChn++){
	    		if (this.fittingStrategy.distortionCalibrationData.gIS[imageSetNumber].imageSet[nChn]!=null) {
	    			imageNumbers[numInSet++]=this.fittingStrategy.distortionCalibrationData.gIS[imageSetNumber].imageSet[nChn].imgNumber;
	    		}
	    	}

	    	if (showGrids) {
	    		for (int nChn=0;nChn<imageNumbers.length;nChn++){
	    			int img_num= imageNumbers[nChn];
	    			this.fittingStrategy.distortionCalibrationData.gIP[img_num].showGridImage();
	    		}
	    	}
	    	if (eachImageInSet){
	    		for (int nChn=0;nChn<imageNumbers.length;nChn++){
	    			int [] imageNumber={imageNumbers[nChn]};
	    			showImageReprojectionErrors(
	    		    		imageNumber, // if null - use all images in a series
	    		    		filter, //weightFunctionFilter,
	    		    		debugLevel);
	    		}
	    		// Do not exit, continue and show combine reprojection errors for all set
	    	}
	    }
	    showImageReprojectionErrors(
	    		imageNumbers, // if null - use all images in a series
	    		filter, //weightFunctionFilter,
	    		debugLevel);
	    return true;

	}

    public void showImageReprojectionErrors(
    		int [] imageNumbers, // if null - use all images in a series
    		int filter,
    		int debugLevel){
    	if (filter<0) filter=this.filterForAll;
    	if (debugLevel>1) {
    		System.out.print("showImageReprojectionErrors: ");
    		if (imageNumbers!=null){
    			for (int i=0;i<imageNumbers.length;i++) System.out.print(" "+imageNumbers[i]);
    		} else {
        		System.out.println("imageNumbers is NULL");
    		}
    	}
	    initFittingSeries(false, filter,this.seriesNumber); // recalculate
		this.currentfX=calculateFxAndJacobian(this.currentVector, false); // is it always true here (this.jacobian==null)
		boolean [] tmpSelectedImages=fittingStrategy.selectedImages();
		boolean [] selectedImages;
		double [] diff=calcYminusFx(this.currentfX);
		if ((imageNumbers!=null) && (imageNumbers.length>0)){
			selectedImages=new boolean[tmpSelectedImages.length];
			for (int i=0;i<selectedImages.length;i++) selectedImages[i]=false;
			for (int i=0;i<imageNumbers.length;i++) if ((imageNumbers[i]>=0) && (imageNumbers[i]<=selectedImages.length)){
				selectedImages[imageNumbers[i]]=tmpSelectedImages[imageNumbers[i]];
			}
		} else {
			selectedImages=tmpSelectedImages;
			int numImg = 0;
			for (int i=0;i<selectedImages.length;i++) if (selectedImages[i]) numImg++;
			imageNumbers = new int [numImg];
			numImg = 0;
			for (int i=0;i<selectedImages.length;i++)  if (selectedImages[i]) {
				imageNumbers[numImg++] = i;
			}
		}
		int width= getGridWidth();
		int height=getGridHeight();
		double [][] imgData=new double[5][height * width]; // dPX, dPY, err
		String [] titles={"dX","dY", "Err","W_Err","Weight"};
		for (int i=0;i<(width*height);i++){
			for (int c=0;c<imgData.length;c++) imgData[c][i]=Double.NaN;
		}
		for (int imgNum=0;imgNum<selectedImages.length;imgNum++) if (selectedImages[imgNum]){
			int len=fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV.length;
			int index=this.imageStartIndex[imgNum]; // pair index
			for (int i=0;i<len;i++){
				int u=fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV[i][0]+this.patternParameters.U0;
				int v=fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV[i][1]+this.patternParameters.V0;
				int vu=u+width*v;
				double w=this.weightFunction[2*(index+i)];
				double dx=diff[2*(index+i)];
				double dy=diff[2*(index+i)+1];
				double e2=dx*dx+dy*dy;
				if (w>0.0){
					if (Double.isNaN(imgData[0][vu])) for (int c=0;c<imgData.length;c++) imgData[c][vu]=0.0;
					imgData[0][vu]+=dx*w;
					imgData[1][vu]+=dy*w;
					imgData[2][vu]+=e2*w;
					imgData[4][vu]+=w;
				}
			}
		}
		int nonEmpty=0;
	    double sumWeights=0.0;
	    for (int vu=0;vu<(width*height);vu++) if (!Double.isNaN(imgData[0][vu])){
	    	nonEmpty++;
	    	sumWeights+=imgData[4][vu];
	    }
		if ((nonEmpty==0) || (sumWeights==0.0)){
			System.out.println("showImageReprojectionErrors():  No non-empty points");
			return;
		}
		double averageWeight=sumWeights/nonEmpty;
	    for (int vu=0;vu<(width*height);vu++) if (!Double.isNaN(imgData[0][vu])){
	    	imgData[0][vu]/=imgData[4][vu];
	    	imgData[1][vu]/=imgData[4][vu];
	    	imgData[2][vu] =Math.sqrt(imgData[2][vu]/imgData[4][vu]);
	    	imgData[3][vu] =Math.sqrt(imgData[2][vu]/averageWeight);
	    }
	    String title="RPRJ";
	    int maxNumInTitle=10;
	    for (int i=0;(i<imageNumbers.length) && (i<maxNumInTitle); i++) title+="-"+imageNumbers[i];
		ShowDoubleFloatArrays.showArrays(
				imgData,
				width,
				height,
				true,
				title,
				titles);
    }




	public double [] calcErrors(double [] diff){
		boolean [] selectedImages=fittingStrategy.selectedImages();
		double [] errors=new double [selectedImages.length];
		int index=0;
		for (int imgNum=0;imgNum<selectedImages.length;imgNum++) {
			errors[imgNum]=Double.NaN; //0.0;
			if (selectedImages[imgNum]) {
				errors[imgNum]=0.0;
				int len=fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV.length*2;
				double w=0;
				if (this.weightFunction!=null) {
					for (int i=index;i<index+len;i++) {
						errors[imgNum]+=diff[i]*diff[i]*this.weightFunction[i];
						w+=this.weightFunction[i];
					}
				} else {
					for (int i=index;i<index+len;i++) {
						errors[imgNum]+=diff[i]*diff[i];
						w+=1.0;
					}
				}
				if (w>0.0) {
					errors[imgNum]/=w;
					errors[imgNum]=Math.sqrt(errors[imgNum])*this.RMSscale;
				} else {
					errors[imgNum]=Double.NaN;
				}
				index+=len;
			}
		}
		return errors;
	}
	/**
	 * Calculate number of used grid points (x/y pairs) for each image in the current fitting series
	 * @return
	 */
	public int [] calcNumPairs(){
		boolean [] selectedImages=fittingStrategy.selectedImages();
		int [] numPairs=new int [selectedImages.length];
		for (int imgNum=0;imgNum<selectedImages.length;imgNum++) {
			numPairs[imgNum]=0;
			if (selectedImages[imgNum]) {
				numPairs[imgNum]=fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV.length;
			}
		}
		return numPairs;
	}
	/**
	 * Calculate corrections to the current parameter values
	 * @param fX calculated grid pixelX, PixelY for current parameter values
	 * @param lambda damping parameter
	 * @return array of deltas to be applied to the coefficients
	 */
	public double [] solveLevenbergMarquardtOldNotUsed(double [] fX, double lambda){
		// calculate JtJ
		double [] diff=calcYminusFx(fX);
		int numPars=this.jacobian.length; // number of parameters to be adjusted
		int length=diff.length; // should be the same as this.jacobian[0].length
	    double [][] JtByJmod=new double [numPars][numPars]; //Transposed Jacobian multiplied by Jacobian
	    double [] JtByDiff=new double [numPars];
	    for (int i=0;i<numPars;i++) for (int j=i;j<numPars;j++){
	    	JtByJmod[i][j]=0.0;
	    	if (this.weightFunction!=null)
	    		for (int k=0;k<length;k++) JtByJmod[i][j]+=this.jacobian[i][k]*this.jacobian[j][k]*this.weightFunction[k];
	    	else
	    		for (int k=0;k<length;k++) JtByJmod[i][j]+=this.jacobian[i][k]*this.jacobian[j][k];
	    }
	    for (int i=0;i<numPars;i++) { // subtract lambda*diagonal , fill the symmetrical half below the diagonal
	    	JtByJmod[i][i]+=lambda*JtByJmod[i][i]; //Marquardt mod
	    	for (int j=0;j<i;j++) JtByJmod[i][j]= JtByJmod[j][i]; // it is symmetrical matrix, just copy
	    }
	    for (int i=0;i<numPars;i++) {
	    	JtByDiff[i]=0.0;
	    	if (this.weightFunction!=null)
		    	for (int k=0;k<length;k++) JtByDiff[i]+=this.jacobian[i][k]*diff[k]*this.weightFunction[k];
	    	else
		    	for (int k=0;k<length;k++) JtByDiff[i]+=this.jacobian[i][k]*diff[k];

	    }
//	    M*Ma=Mb
	    Matrix M=new Matrix(JtByJmod);
//  public Matrix (double vals[], int m) {
/*
		if (this.debugLevel>2) {
			for (int n=0;n<fittingStrategy.distortionCalibrationData.pixelsXY.length;n++) {
				for (int i=0;i<fittingStrategy.distortionCalibrationData.pixelsXY[n].length;i++){
					System.out.println(n+":"+i+"  "+
							fittingStrategy.distortionCalibrationData.pixelsUV[n][i][0]+"/"+
							fittingStrategy.distortionCalibrationData.pixelsUV[n][i][1]+"  "+
							IJ.d2s(fittingStrategy.distortionCalibrationData.pixelsXY[n][i][0], 2)+"/"+
							IJ.d2s(fittingStrategy.distortionCalibrationData.pixelsXY[n][i][1], 2)
					);
				}
			}
		}

 */

		if (this.debugLevel>3) {
//		if (this.debugLevel>1) {
			System.out.println("Jt*J -lambda* diag(Jt*J), lambda="+lambda+":");
			M.print(10, 5);
		}

	    Matrix Mb=new Matrix(JtByDiff,numPars); // single column
	    if (!(new LUDecomposition(M)).isNonsingular()){
	    	double [][] arr=M.getArray();
			System.out.println("Singular Matrix "+arr.length+"x"+arr[0].length);
			// any rowsx off all 0.0?
			for (int n=0;n<arr.length;n++){
				boolean zeroRow=true;
				for (int i=0;i<arr[n].length;i++) if (arr[n][i]!=0.0){
					zeroRow=false;
					break;
				}
				if (zeroRow){
					System.out.println("Row of all zeros: "+n);
				}
			}
//			M.print(10, 5);
	    	return null;
	    }

//	    Matrix Ma=M.solve(Mb); // singular
	    if (this.debugLevel>0) System.out.print("Running Cholesky decomposition...");
	    long decompositionTime=System.nanoTime();
	    Matrix Ma=M.chol().solve(Mb); // singular
	    decompositionTime=System.nanoTime()-decompositionTime;
	    if (this.debugLevel>0) System.out.println("done in "+(decompositionTime/1E9)+" sec");

	    return Ma.getColumnPackedCopy();
	}



	public LMAArrays calculateJacobianArrays(double [] fX){
		// calculate JtJ
		double [] diff=calcYminusFx(fX);
		int numPars=this.jacobian.length; // number of parameters to be adjusted
		int length=diff.length; // should be the same as this.jacobian[0].length
	    double [][] JtByJmod=new double [numPars][numPars]; //Transposed Jacobian multiplied by Jacobian
	    double [] JtByDiff=new double [numPars];
	    for (int i=0;i<numPars;i++) for (int j=i;j<numPars;j++){
	    	JtByJmod[i][j]=0.0;
	    	if (this.weightFunction!=null)
	    		for (int k=0;k<length;k++) JtByJmod[i][j]+=this.jacobian[i][k]*this.jacobian[j][k]*this.weightFunction[k];
	    	else
	    		for (int k=0;k<length;k++) JtByJmod[i][j]+=this.jacobian[i][k]*this.jacobian[j][k];
	    }
	    for (int i=0;i<numPars;i++) { // subtract lambda*diagonal , fill the symmetrical half below the diagonal
//	    	JtByJmod[i][i]+=lambda*JtByJmod[i][i]; //Marquardt mod
	    	for (int j=0;j<i;j++) JtByJmod[i][j]= JtByJmod[j][i]; // it is symmetrical matrix, just copy
	    }
	    for (int i=0;i<numPars;i++) {
	    	JtByDiff[i]=0.0;
	    	if (this.weightFunction!=null)
		    	for (int k=0;k<length;k++) JtByDiff[i]+=this.jacobian[i][k]*diff[k]*this.weightFunction[k];
	    	else
		    	for (int k=0;k<length;k++) JtByDiff[i]+=this.jacobian[i][k]*diff[k];

	    }

   		LMAArrays lMAArrays = new LMAArrays();
   		lMAArrays.jTByJ=JtByJmod;
   		lMAArrays.jTByDiff=JtByDiff;
   		return lMAArrays;
/*

//	    M*Ma=Mb
	    Matrix M=new Matrix(JtByJmod);
//  public Matrix (double vals[], int m) {

		if (this.debugLevel>2) {
//		if (this.debugLevel>1) {
			System.out.println("Jt*J -lambda* diag(Jt*J), lambda="+lambda+":");
			M.print(10, 5);
		}

	    Matrix Mb=new Matrix(JtByDiff,numPars); // single column
	    if (!(new LUDecomposition(M)).isNonsingular()){
			System.out.println("Singular Matrix");
			M.print(10, 5);
	    	return null;
	    }
	    Matrix Ma=M.solve(Mb); // singular
	    return Ma.getColumnPackedCopy();
*/
	}



	public LMAArrays calculateJacobianArrays (
			final boolean [] selectedImages, // selected images to process
			final double [] Y,  // should be initialized
			final double [] fX, // should be initialized to correct length, data is not needed
			final double [] vector,  // parameters vector
			final int    [] imageStartIndex, // index of the first point of each image (including extra element in the end so n+1 is always valid)
			final double [][] patternXYZ, // this.targetXYZ
			final double [] weightFunction, // may be null - make it twice smaller? - same for X and Y?
			final LensDistortionParameters lensDistortionParametersProto,
			final double [][] dTA_dUV, // null or double [][] to return averaged per-image {{dU/dAz,dU/dTl}{dV/dAz,dV/dTl}}
			int threadsMax,
			boolean updateStatus){

		// calculate JtJ
//		double [] diff=calcYminusFx(fX);
//		int numPars=this.jacobian.length; // number of parameters to be adjusted
		final int numPars=vector.length; // number of parameters to be adjusted
//		int length=diff.length; // should be the same as this.jacobian[0].length
//		final int length=fX.length; // should be the same as this.jacobian[0].length
		final double [][] JtByJmod=new double [numPars][numPars]; //Transposed Jacobian multiplied by Jacobian
		final double [] JtByDiff=new double [numPars];
	    for (int i=0;i<numPars;i++){
	    	JtByDiff[i]=0.0;
	    	for (int j=0;j<numPars;j++) JtByJmod[i][j]=0.0;
	    }
	    final int debugLevel=this.debugLevel;
   		final Thread[] threads = newThreadArray(threadsMax);
   		final AtomicInteger imageNumberAtomic = new AtomicInteger(0);
   		final AtomicInteger imageFinishedAtomic = new AtomicInteger(0);
   		final double [] progressValues=new double [selectedImages.length];
   		int numSelectedImages=0;
   		for (int i=0;i<selectedImages.length;i++) if (selectedImages[i]) numSelectedImages++;
   		int selectedIndex=0;
   		for (int i=0;i<selectedImages.length;i++) {
   			progressValues[i]=(selectedIndex+1.0)/numSelectedImages;
   			if (selectedImages[i]) selectedIndex++;
   			if (selectedIndex>=numSelectedImages) selectedIndex--;
   		}
   		final double [][][] dUV_image =   (dTA_dUV != null) ? new double [selectedImages.length][][]: null;
   		final double []     dUV_weights = (dTA_dUV != null) ? new double [selectedImages.length]: null;
   		final AtomicInteger stopRequested=this.stopRequested;
		final AtomicBoolean interruptedAtomic=new AtomicBoolean();
   		for (int ithread = 0; ithread < threads.length; ithread++) {
   			threads[ithread] = new Thread() {
   				@Override
				public void run() {
   					LensDistortionParameters lensDistortionParameters=lensDistortionParametersProto.clone(); // see - if that is needed - maybe new is OK
   					//   					LensDistortionParameters lensDistortionParameters= new LensDistortionParameters();
   					for (int numImage=imageNumberAtomic.getAndIncrement(); (numImage<selectedImages.length) && !interruptedAtomic.get();numImage=imageNumberAtomic.getAndIncrement()){
   						int length=2*(imageStartIndex[numImage+1]-imageStartIndex[numImage]);
   						if (length == 0) {
   							continue;
   						}
   						int start= 2*imageStartIndex[numImage];
   						double [][] partialJacobian= calculatePartialFxAndJacobian(
   								numImage,      // number of grid image
   								vector,  // parameters vector
   								patternXYZ, // this.targetXYZ
   								fX,     // non-overlapping segments will be filled
   								imageStartIndex, // start index in patternXYZ array (length - difference to the next, includes extra last element)
   								lensDistortionParameters, // initialize one per each tread? Or for each call?
   								true); // when false, modifies only this.lensDistortionParameters.*

//   						int length=2*(imageStartIndex[numImage+1]-imageStartIndex[numImage]);
//   						int start= 2*imageStartIndex[numImage];

   						double [][] partialJtByJmod=new double [numPars][numPars]; // out of heap space
   						double []   partialJtByDiff=new double [numPars];

   						for (int i=0;i<numPars;i++) if (partialJacobian[i]!=null) {
   							for (int j=i;j<numPars;j++) if (partialJacobian[j]!=null) {
   								partialJtByJmod[i][j]=0.0;
   								if (weightFunction!=null) {
   									for (int k=0;k<length;k++) partialJtByJmod[i][j]+=partialJacobian[i][k]*partialJacobian[j][k]*weightFunction[start+k];
   								} else {
   									for (int k=0;k<length;k++) partialJtByJmod[i][j]+=partialJacobian[i][k]*partialJacobian[j][k];
   								}
   							}
   						}

   						double [] partialDiff=new double[length];
   						for (int k=0;k<length;k++) 	partialDiff[k]=Y[start+k]-fX[start+k];

   						for (int i=0;i<numPars;i++) if (partialJacobian[i]!=null) {
   							partialJtByDiff[i]=0.0;
   							if (weightFunction!=null)
   								for (int k=0;k<length;k++) {
   									if (Double.isNaN(partialDiff[k])) {
   										System.out.println("calculateJacobianArrays() BUG1:partialDiff["+k+"]=NaN, i="+i);
   									} else {
   										partialJtByDiff[i]+=partialJacobian[i][k]*partialDiff[k]*weightFunction[start+k];
   									}
   								}
   							else
   								for (int k=0;k<length;k++) {
   									if (Double.isNaN(partialDiff[k])) {
   										System.out.println("calculateJacobianArrays() BUG2:partialDiff["+k+"]=NaN, i="+i);
   									} else {
   										partialJtByDiff[i]+=partialJacobian[i][k]*partialDiff[k];
   									}
   								}

   						}
   						int par_gh = -1;
   						int par_ga = -1;
   						if (dUV_image != null) {
   							dUV_image[numImage] = null; //new double [2][2];
   							for (int i = 0; i < numPars; i++) {
   								int parNum=fittingStrategy.parameterMap[i][1];
   								if (parNum == fittingStrategy.distortionCalibrationData.index_gh) {
   									par_gh = i;
   								} else if (parNum == fittingStrategy.distortionCalibrationData.index_ga) {
   									par_ga = i;
   								}
   							}
   							if ((par_gh >= 0) && (par_ga >= 0)) { // both defined
   								// partialJacobian[par_gh][2*k  ] dPx/dGh
   								// partialJacobian[par_gh][2*k+1] dPy/dGh
   								// partialJacobian[par_gh][2*k  ] dPx/dGa
   								// partialJacobian[par_gh][2*k+1] dPy/dGa
   								// d3780:	public void updateGridToPointer(ImagePlus imp_grid, double[][] xyuv) {
   								// find average dX/dU, dY/dU, dX/dV, dY/dV 
   								
   								
   								PolynomialApproximation polynomialApproximation =new PolynomialApproximation(0);// no debug
   								GridImageParameters gip = fittingStrategy.distortionCalibrationData.gIP[numImage];
  								int np = gip.pixelsXY.length;
  								double wsx=0.0, wsy=0.0, ws_dpx_dgh= 0.0, ws_dpy_dgh= 0.0, ws_dpx_dga= 0.0, ws_dpy_dga= 0.0;
   								for (int k = 0; k < length; k+=1) {
   									double w = (weightFunction!=null)? weightFunction[start+k] : 1.0;
   									wsx +=        w;
   									ws_dpx_dgh += w * partialJacobian[par_gh][k];
   									ws_dpx_dga += w * partialJacobian[par_ga][k];
   									k++;
   									w = (weightFunction!=null)? weightFunction[start+k] : 1.0;
   									wsy +=        w;
   									ws_dpy_dgh += w * partialJacobian[par_gh][k];
   									ws_dpy_dga += w * partialJacobian[par_ga][k];
   								}
   								if ((wsx == 0.0) || (wsy == 0.0)) {
   									if (debugLevel>2) {
   										System.out.println("Not enough data for dpx/dgh, dpy/dgh, dpx/dga, dpy/dga for image #"+numImage);
   									}
   									continue;
   								}
   								ws_dpx_dgh /= wsx;
   								ws_dpx_dga /= wsx;
   								ws_dpy_dgh /= wsy;
   								ws_dpy_dga /= wsy;
   								
   								Matrix mXY_HA = new Matrix(new double[][] {{ws_dpx_dgh, ws_dpx_dga},{ws_dpy_dgh, ws_dpy_dga}});
  								
   								double [][][] data = new double[np][3][];
   								for (int indx = 0; indx < np; indx++) {
   									data[indx][0] = new double[2];
   									data[indx][1] = new double[2];
   									data[indx][2] = new double[1];
   									data[indx][0][0] = gip.pixelsUV[indx][0]; // U
   									data[indx][0][1] = gip.pixelsUV[indx][1]; // V
   									data[indx][1][0] = gip.pixelsXY[indx][0]; // pX
   									data[indx][1][1] = gip.pixelsXY[indx][1]; // pY
   									data[indx][2][0] = gip.pixelsMask[indx];  // weighth
   								}
   								double [][] coeff = polynomialApproximation.quadraticApproximation(
   										data,
   										true); // force linear
   								Matrix mXY_UV=new Matrix(new double[][] {{coeff[0][0],coeff[0][1]},{coeff[1][0],coeff[1][1]}});
   								if (!(new LUDecomposition(mXY_UV)).isNonsingular()){
   									if (debugLevel>2) {
   										System.out.println("Skipping singular matrix for image #"+numImage);
   									}
   									continue;
   								}
   								Matrix mUV_XY= mXY_UV.inverse();
   								Matrix mUV_HA = mUV_XY.times(mXY_HA);
   								dUV_image  [numImage] = mUV_HA.getArray();
   								dUV_weights[numImage] = wsx + wsy;
   							} else {
   								System.out.println ("Bug: ga/gh are not among the parameters, par_gh="+par_gh+", par_ga="+par_ga);
   							}
   						}
   						
   						synchronizedCombinePartialJacobians(
   								JtByJmod, //Transposed Jacobian multiplied by Jacobian
   								JtByDiff,
   								partialJacobian,
   								partialJtByDiff,
   								partialJtByJmod,
   								numPars	);

   						final int numFinished=imageFinishedAtomic.getAndIncrement();
   						SwingUtilities.invokeLater(new Runnable() {
   							@Override
							public void run() {
   								// Here, we can safely update the GUI
   								// because we'll be called from the
   								// event dispatch thread
   								IJ.showProgress(progressValues[numFinished]);
   							}
   						});
   						if (stopRequested.get()==1){ // ASAP
   							interruptedAtomic.set(true);
   						}
   					}
   				}
   			};
   		}
   		startAndJoin(threads);
   		if (interruptedAtomic.get()) {
   			System.out.println("calculateJacobianArrays() aborted by user request");
   			return null;
   		}
   		if (debugLevel>3){
   			String msg="calculateJacobianArrays() ALL_trace=";
   			for (int ii=0;ii<numPars;ii++) msg+=IJ.d2s(JtByJmod[ii][ii],5);
   			System.out.println(msg);

   		}

   		for (int i=0;i<numPars;i++) { // subtract lambda*diagonal , fill the symmetrical half below the diagonal
   			for (int j=0;j<i;j++) JtByJmod[i][j]= JtByJmod[j][i]; // it is symmetrical matrix, just copy
   		}
   		LMAArrays lMAArrays = new LMAArrays();
   		lMAArrays.jTByJ=JtByJmod;
   		lMAArrays.jTByDiff=JtByDiff;
   		if (debugLevel>3){
   			String msg="calculateJacobianArrays() lMAArrays.jTByJ trace=";
   			for (int ii=0;ii<numPars;ii++) msg+=IJ.d2s(lMAArrays.jTByJ[ii][ii],5)+" ";
   			System.out.println(msg);

   		}
   		if (dTA_dUV != null) {
   			double [][] dUV_average = new double [2][2];
   			for (int i = 0; i < 2; i++) for (int j=0; j<2; j++) dUV_average[i][j] = 0.0;
   			double sw = 0.0;
   			for (int ni = 0; ni < dUV_weights.length; ni++) {
   				double w = dUV_weights[ni]; 
   				if (w > 0.0) {
   					for (int i = 0; i < 2; i++) for (int j=0; j<2; j++) {
   						dUV_average[i][j] += w *dUV_image[ni][i][j] ;
   					}
   					sw += w;
   				}
   			}
   			if (sw > 0) {
   				for (int i = 0; i < 2; i++) for (int j=0; j<2; j++) {
   					dUV_average[i][j] /= sw;
   				}
   				Matrix mdTA_dUV = (new Matrix(dUV_average)).inverse();
   				this.dTA_dUV = mdTA_dUV.getArray();
   			}
   		}
   		return lMAArrays;
	}

	public synchronized void synchronizedCombinePartialJacobians(
			double [][] JtByJmod, //Transposed Jacobian multiplied by Jacobian
			double []   JtByDiff,
			double [][] partialJacobian,
			double []   partialJtByDiff,
			double [][] partialJtByJmod,
			int numPars
	){
		for (int i=0;i<numPars;i++) if (partialJacobian[i]!=null){
			JtByDiff[i]+=partialJtByDiff[i];
			for (int j=i;j<numPars;j++) JtByJmod[i][j]+=partialJtByJmod[i][j];
		}
	}




	public double [] solveLMA(
			LMAArrays lMAArrays,
			double lambda){
		double [][] JtByJmod= lMAArrays.jTByJ.clone();
		int numPars=JtByJmod.length;
		for (int i=0;i<numPars;i++){
			JtByJmod[i]=lMAArrays.jTByJ[i].clone();
   			JtByJmod[i][i]+=lambda*JtByJmod[i][i]; //Marquardt mod
		}
//	    M*Ma=Mb
		System.out.println("JtByJmod.length="+JtByJmod.length+" numPars="+numPars);
		if (numPars==0) {
			return null;
		}
	    Matrix M=new Matrix(JtByJmod);
//  public Matrix (double vals[], int m) {

		if (this.debugLevel>2) {
//		if (this.debugLevel>1) {
			System.out.println("Jt*J -lambda* diag(Jt*J), lambda="+lambda+":");
			M.print(10, 5);
		}

	    Matrix Mb=new Matrix(lMAArrays.jTByDiff,numPars); // single column {NaN,NaN}
	    if (!(new LUDecomposition(M)).isNonsingular()){
	    	double [][] arr=M.getArray();
			System.out.println("Singular Matrix "+arr.length+"x"+arr[0].length);
			// any rowsx off all 0.0?
			for (int n=0;n<arr.length;n++){
				boolean zeroRow=true;
				for (int i=0;i<arr[n].length;i++) if (arr[n][i]!=0.0){
					zeroRow=false;
					break;
				}
				if (zeroRow){
					System.out.println("Row of all zeros: "+n);
				}
			}
//			M.print(10, 5);
	    	return null;
	    }
	    Matrix Ma=M.solve(Mb); // singular
	    return Ma.getColumnPackedCopy();

	}


	/**
	 * Calculates  next parameters vector, holds some arrays
	 * @param numSeries
	 * @return array of two booleans: { improved, finished}
	 */
	public boolean [] stepLevenbergMarquardtFirst(
			int     numSeries,
			boolean calc_dUV){
		double [] deltas=null;
		if ((this.currentVector==null) || (this.currentVector.length==0)) { // length==0 was debugging
			int filter=this.filterForAll;
			if (this.askFilter) filter=selectFilter(filter);
			initFittingSeries(false,filter,numSeries); // first step in series
			this.currentRMS=-1;
			this.currentRMSPure=-1;
			this.currentfX=null; // invalidate
			this.jacobian=null;  // invalidate
			this.lMAArrays=null;
			lastImprovements[0]=-1.0;
			lastImprovements[1]=-1.0;
		}
		// calculate  this.currentfX, this.jacobian if needed
		if (this.debugLevel>2) {
			System.out.println("this.currentVector");
			for (int i=0;i<this.currentVector.length;i++){
				System.out.println(i+": "+ this.currentVector[i]);
			}
		}
		//    	if ((this.currentfX==null)|| ((this.jacobian==null) && !this.threadedLMA )) {
		if ((this.currentfX==null)|| (this.lMAArrays==null)) {
			if (this.updateStatus){
//				IJ.showStatus(this.seriesNumber+": "+"Step #"+this.iterationStepNumber+" RMS="+IJ.d2s(this.currentRMS,8)+ " ("+IJ.d2s(this.firstRMS,8)+")");
				IJ.showStatus(this.seriesNumber+": initial Jacobian matrix calculation. Points:"+this.Y.length+" Parameters:"+this.currentVector.length);
			}
			if (this.debugLevel >1) {
				System.out.println(this.seriesNumber+": initial Jacobian matrix calculation. Points:"+this.Y.length+" Parameters:"+this.currentVector.length);
			}
    		if (this.threadedLMA) {
    			this.currentfX=new double[this.Y.length];
   				this.dTA_dUV = 	calc_dUV ? (new double [2][2]): null;
    			//    			deltas=solveLevenbergMarquardtThreaded(
    			this.lMAArrays=calculateJacobianArrays(
    					this.fittingStrategy.selectedImages(), // selected images to process
    					this.Y,  // should be initialized
    					this.currentfX, // should be initialized to correct length, data is not needed
    					this.currentVector,  // parameters vector
    					this.imageStartIndex, // index of the first point of each image (including extra element in the end so n+1 is always valid)
    					this.targetXYZ, // this.targetXYZ
    					this.weightFunction, // may be null - make it twice smaller? - same for X and Y?
    					this.lensDistortionParameters,
    					this.dTA_dUV,// final double [][][] dUV_average, // null or double [selectedImages.length][][] to return per-image {{dU/dAz,dU/dTl}{dV/dAz,dV/dTl}}
    					//    					this.lambda,
    					this.threadsMax,
    					this.updateStatus);
    			if (this.lMAArrays == null) {
    				return null ; // aborted
    			}
    		} else {
    			this.currentfX=calculateFxAndJacobian(this.currentVector, true); // is it always true here (this.jacobian==null)
    			this.lMAArrays=calculateJacobianArrays(this.currentfX);
//    			deltas=solveLevenbergMarquardt(this.currentfX,this.lambda);
    		}
    		// add termes that push selected extrinsic parameters towards average (global, per station, per tilt-station)
    		this.currentRMSPure= calcErrorDiffY(this.currentfX) +0.0;
    	   	if ((this.fittingStrategy.varianceModes!=null) && (this.fittingStrategy.varianceModes[numSeries]!=this.fittingStrategy.varianceModeDisabled)) {
    	   		this.fittingStrategy.addVarianceToLMA(
    	    			numSeries,
    	    			this.currentVector,
    	    			this.lMAArrays.jTByJ, // jacobian multiplied by Jacobian transposed (or null)
    	    			this.lMAArrays.jTByDiff);
    	   		this.currentRMS= calcErrorDiffY(
    	   				this.currentfX,
    	   				this.fittingStrategy.getVarianceError2(), //double [] extraWeightedErrors,
    	   				this.fittingStrategy.getWeights()); //double [] extraWeights);
    			if (this.updateStatus){
    				IJ.showStatus(this.seriesNumber+": initial RMS="+IJ.d2s(this.currentRMS,8)+
    						" ("+IJ.d2s(this.currentRMSPure,8)+")"+
    						". Calculating next Jacobian. Points:"+this.Y.length+" Parameters:"+this.currentVector.length);
    			}
    			if ((this.debugLevel>0) && ((this.debugLevel>1) || ((System.nanoTime()-this.startTime)>10000000000.0))) {
    				System.out.println(this.seriesNumber+": initial RMS="+IJ.d2s(this.currentRMS,8)+
    						" ("+IJ.d2s(this.currentRMSPure,8)+")"+
    						". Calculating next Jacobian. Points:"+this.Y.length+" Parameters:"+this.currentVector.length);
    			}

    	   	} else {
        		this.currentRMS= this.currentRMSPure;
    			if (this.updateStatus){
    				IJ.showStatus(this.seriesNumber+": initial RMS="+IJ.d2s(this.currentRMS,8)+
    						". Calculating next Jacobian. Points:"+this.Y.length+" Parameters:"+this.currentVector.length);
    			}
    			if (this.debugLevel>1) {
    				System.out.println(this.seriesNumber+": initial RMS="+IJ.d2s(this.currentRMS,8)+
    						". Calculating next Jacobian. Points:"+this.Y.length+" Parameters:"+this.currentVector.length);
    			}
    	   	}
    	} else {
    		this.currentRMSPure= calcErrorDiffY(this.currentfX);
    	   	if ((this.fittingStrategy.varianceModes!=null) && (this.fittingStrategy.varianceModes[numSeries]!=this.fittingStrategy.varianceModeDisabled)) {
    	   		this.fittingStrategy.addVarianceToLMA(// recalculating as this may keep from nextVector (or just being restored)
    	    			numSeries,
    	    			this.currentVector,
    	    			null, //this.lMAArrays.jTByJ, // jacobian multiplied by Jacobian transposed (or null)
    	    			null); //this.lMAArrays.jTByDiff);

    	   		this.currentRMS= calcErrorDiffY(
    	   				this.currentfX,
    	   				this.fittingStrategy.getVarianceError2(), //double [] extraWeightedErrors,
    	   				this.fittingStrategy.getWeights()); //double [] extraWeights);
    	   	} else {
        		this.currentRMS= this.currentRMSPure;
    	   	}

    	}
//		this.currentRMS= calcError(calcYminusFx(this.currentfX));
    	if (this.firstRMS<0) {
    		this.firstRMS=this.currentRMS;
    		this.firstRMSPure=this.currentRMSPure;
    	}
// calculate deltas
//    	double [] deltas=solveLevenbergMarquardt(this.currentfX,fittingStrategy.getLambda());

		deltas=solveLMA(this.lMAArrays,	this.lambda	);

    	boolean matrixNonSingular=true;
    	if (deltas==null) {
    		deltas=new double[this.currentVector.length];
    		for (int i=0;i<deltas.length;i++) deltas[i]=0.0;
    		matrixNonSingular=false;
    	}
		if (this.debugLevel>2) {
			System.out.println("deltas");
			for (int i=0;i<deltas.length;i++){
				System.out.println(i+": "+ deltas[i]);
			}
		}
// apply deltas
    	this.nextVector=this.currentVector.clone();
    	for (int i=0;i<this.nextVector.length;i++) this.nextVector[i]+=deltas[i];
// another option - do not calculate J now, just fX. and late - calculate both if it was improvement
//    	save current Jacobian
		if (this.debugLevel>2) {
			System.out.println("this.nextVector");
			for (int i=0;i<this.nextVector.length;i++){
				System.out.println(i+": "+ this.nextVector[i]);
			}
		}

//        this.savedJacobian=this.jacobian;
        this.savedLMAArrays=lMAArrays.clone();
        this.jacobian=null; // not needed, just to catch bugs
// calculate next vector and Jacobian  (this.jacobian)
//    	this.nextfX=calculateFxAndJacobian(this.nextVector,true); //=========== OLD

		if (this.threadedLMA) {
			this.nextfX=new double[this.Y.length];
			//    			deltas=solveLevenbergMarquardtThreaded(
			this.lMAArrays=calculateJacobianArrays(
					this.fittingStrategy.selectedImages(), // selected images to process
					this.Y,  // should be initialized
					this.nextfX, // should be initialized to correct length, data is not needed
					this.nextVector,  // parameters vector
					this.imageStartIndex, // index of the first point of each image (including extra element in the end so n+1 is always valid)
					this.targetXYZ, // this.targetXYZ
					this.weightFunction, // may be null - make it twice smaller? - same for X and Y?
					this.lensDistortionParameters,
					this.dTA_dUV,// final double [][][] dUV_average, // null or double [selectedImages.length][][] to return per-image {{dU/dAz,dU/dTl}{dV/dAz,dV/dTl}}
					//    					this.lambda,
					this.threadsMax,
					this.updateStatus);
			if (this.lMAArrays == null) {
				return null ; // aborted
			}
		} else {
	    	this.nextfX=calculateFxAndJacobian(this.nextVector,true);
			this.lMAArrays=calculateJacobianArrays(this.nextfX);
		}
//		this.nextRMS=calcErrorDiffY(this.nextfX);

		this.nextRMSPure= calcErrorDiffY(this.nextfX);
	   	if ((this.fittingStrategy.varianceModes!=null) && (this.fittingStrategy.varianceModes[numSeries]!=this.fittingStrategy.varianceModeDisabled)) {
	   		this.fittingStrategy.addVarianceToLMA(
	    			numSeries,
	    			this.nextVector,
	    			this.lMAArrays.jTByJ, // jacobian multiplied by Jacobian transposed (or null)
	    			this.lMAArrays.jTByDiff);
	   		this.nextRMS= calcErrorDiffY(
	   				this.nextfX,
	   				this.fittingStrategy.getVarianceError2(), //double [] extraWeightedErrors,
	   				this.fittingStrategy.getWeights()); //double [] extraWeights);
	   	} else {
	   		this.nextRMS= this.nextRMSPure;
	   	}

		this.lastImprovements[1]=this.lastImprovements[0];
		this.lastImprovements[0]=this.currentRMS-this.nextRMS;
		if (this.debugLevel>2) {
			System.out.println("stepLMA this.currentRMS="+this.currentRMS+
					", this.currentRMSPure="+this.currentRMSPure+
					", this.nextRMS="+this.nextRMS+
					", this.nextRMSPure="+this.nextRMSPure+
					", delta="+(this.currentRMS-this.nextRMS)+
					", deltaPure="+(this.currentRMSPure-this.nextRMSPure));
		}
		boolean [] status={matrixNonSingular && (this.nextRMS<=this.currentRMS),!matrixNonSingular};
		// additional test if "worse" but the difference is too small, it was be caused by computation error, like here:
		//stepLevenbergMarquardtAction() step=27, this.currentRMS=0.17068403807026408,   this.nextRMS=0.1706840380702647

		if (!status[0] && matrixNonSingular) {
			if (this.nextRMS<(this.currentRMS+this.currentRMS*this.thresholdFinish*0.01)) {
				this.nextRMS=this.currentRMS;
				this.nextRMSPure=this.currentRMSPure;
				status[0]=true;
				status[1]=true;
				this.lastImprovements[0]=0.0;
				if (this.debugLevel>1) {
					System.out.println("New RMS error is larger than the old one, but the difference is too small to be trusted ");
					System.out.println(
							"stepLMA this.currentRMS="+this.currentRMS+
							", this.currentRMSPure="+this.currentRMSPure+
							", this.nextRMS="+this.nextRMS+
							", this.nextRMSPure="+this.nextRMSPure+
							", delta="+(this.currentRMS-this.nextRMS)+
							", deltaPure="+(this.currentRMSPure-this.nextRMSPure));
				}

			}
		}
    	if (status[0] && matrixNonSingular) { //improved
    		status[1]=(this.iterationStepNumber>this.numIterations) || ( // done
    				(this.lastImprovements[0]>=0.0) &&
    				(this.lastImprovements[0]<this.thresholdFinish*this.currentRMS) &&
    				(this.lastImprovements[1]>=0.0) &&
    				(this.lastImprovements[1]<this.thresholdFinish*this.currentRMS));
    	} else if (matrixNonSingular){
//    		this.jacobian=this.savedJacobian;// restore saved Jacobian
    		this.lMAArrays=this.savedLMAArrays; // restore Jt*J and Jt*diff

    		status[1]=(this.iterationStepNumber>this.numIterations) || // failed
    		((this.lambda*this.lambdaStepUp)>this.maxLambda);
    	}
///this.currentRMS
//TODO: add other failures leading to result failure?
		if (this.debugLevel>2) {
			System.out.println("stepLevenbergMarquardtFirst("+numSeries+")=>"+status[0]+","+status[1]);
		}
		return status;
    }
    /**
     * Apply fitting step
     */
    public void stepLevenbergMarquardtAction(){//
    	this.iterationStepNumber++;
// apply/revert,modify lambda
		if (this.debugLevel>1) {
			System.out.println(
					"stepLevenbergMarquardtAction() step="+this.iterationStepNumber+
					", this.currentRMS="+this.currentRMS+
					", this.currentRMSPure="+this.currentRMSPure+
					", this.nextRMS="+this.nextRMS+
					", this.nextRMSPure="+this.nextRMSPure+
					" lambda="+this.lambda+" at "+IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec");
		}
    	if (this.nextRMS<this.currentRMS) { //improved
    		this.lambda*=this.lambdaStepDown;
    		this.currentRMS=this.nextRMS;
    		this.currentRMSPure=this.nextRMSPure;
    		this.currentfX=this.nextfX;
    		this.currentVector=this.nextVector;
    	} else {
    		this.lambda*=this.lambdaStepUp;
//    		this.jacobian=this.savedJacobian;// restore saved Jacobian
    		this.lMAArrays=this.savedLMAArrays; // restore Jt*J and Jt*diff
    	}
    }

    /**
     * Dialog to select Levenberg-Marquardt algorithm and related parameters
     * @return true if OK, false if canceled
     */
    public boolean selectLMAParameters(){
    	int numSeries=fittingStrategy.getNumSeries();
	    GenericDialog gd = new GenericDialog("Levenberg-Marquardt algorithm parameters for cameras distortions/locations");
		gd.addNumericField("Iteration number to start (0.."+(numSeries-1)+")", this.seriesNumber, 0);
		gd.addNumericField("Initial LMA Lambda ",            this.lambda, 5);
		gd.addNumericField("Multiply lambda on success",     this.lambdaStepDown, 5);
		gd.addNumericField("Threshold RMS to exit LMA",      this.thresholdFinish, 7,9,"pix");
		gd.addNumericField("Multiply lambda on failure",     this.lambdaStepUp, 5);
		gd.addNumericField("Threshold lambda to fail",       this.maxLambda, 5);
		gd.addNumericField("Maximal number of iterations",   this.numIterations, 0);

		gd.addCheckbox("Dialog after each iteration step",   this.stopEachStep);
		gd.addCheckbox("Dialog after each iteration series", this.stopEachSeries);
		gd.addCheckbox("Dialog after each failure",          this.stopOnFailure);
		gd.addCheckbox("Ask for weight function filter",     this.askFilter);

		gd.addCheckbox("Show modified parameters",           this.showParams);
		gd.addCheckbox("Show debug images before correction",this.showThisImages);
		gd.addCheckbox("Show debug images after correction", this.showNextImages);
		gd.addNumericField("Maximal number of threads",   this.threadsMax, 0);
		gd.addCheckbox("Use memory-saving/multithreaded version", this.threadedLMA);
	    gd.showDialog();
	    if (gd.wasCanceled()) return false;
	    this.seriesNumber=     (int) gd.getNextNumber();
		this.lambda=                 gd.getNextNumber();
		this.lambdaStepDown=         gd.getNextNumber();
		this.thresholdFinish=        gd.getNextNumber();
		this.lambdaStepUp=           gd.getNextNumber();
		this.maxLambda=              gd.getNextNumber();
		this.numIterations=    (int) gd.getNextNumber();
		this.stopEachStep=           gd.getNextBoolean();
		this.stopEachSeries=         gd.getNextBoolean();
		this.stopOnFailure=          gd.getNextBoolean();
		this.askFilter=              gd.getNextBoolean();
		this.showParams=             gd.getNextBoolean();
		this.showThisImages=         gd.getNextBoolean();
		this.showNextImages=         gd.getNextBoolean();
		this.threadsMax=       (int) gd.getNextNumber();
		this.threadedLMA=            gd.getNextBoolean();
	    return true;
    }


    public boolean dialogLMAStep(boolean [] state){
    	String [] states={
    			"Worse, increase lambda",
    			"Better, decrease lambda",
    			"Failed to fit",
    			"Fitting Successful"};
    	int iState=(state[0]?1:0)+(state[1]?2:0);

	    GenericDialog gd = new GenericDialog("Levenberg-Marquardt algorithm step");
//    	String [][] parameterDescriptions=fittingStrategy.distortionCalibrationData.parameterDescriptions;
    	gd.addMessage("Current state="+states[iState]);
    	gd.addMessage("Iteration step="+this.iterationStepNumber);

    	gd.addMessage("Initial RMS="+IJ.d2s(this.firstRMS,6)+", Current RMS="+IJ.d2s(this.currentRMS,6)+", new RMS="+IJ.d2s(this.nextRMS,6));
    	gd.addMessage("Pure initial RMS="+IJ.d2s(this.firstRMSPure,6)+", Current RMS="+IJ.d2s(this.currentRMSPure,6)+", new RMS="+IJ.d2s(this.nextRMSPure,6));
    	if (this.showParams) {
    		for (int i=0;i<this.currentVector.length;i++){
    			int parNum=fittingStrategy.parameterMap[i][1];
    			int imgNum=fittingStrategy.parameterMap[i][0];
    			double delta= this.nextVector[i] - this.currentVector[i];
//    			gd.addMessage(i+": "+parameterDescriptions[parNum][0]+
//    					"["+imgNum+"]("+parameterDescriptions[parNum][2]+") "+IJ.d2s(this.currentVector[i],3)+
//    					" + "+IJ.d2s(delta,3)+" = "+IJ.d2s(this.nextVector[i],3));
    			gd.addMessage(i+": "+fittingStrategy.distortionCalibrationData.descrField(parNum,0)+
    					"["+imgNum+"]("+fittingStrategy.distortionCalibrationData.descrField(parNum,2)+") "+IJ.d2s(this.currentVector[i],3)+
    					" + "+IJ.d2s(delta,3)+" = "+IJ.d2s(this.nextVector[i],3));
    		}
    	}
		gd.addNumericField("Lambda ",                        this.lambda, 5);
		gd.addNumericField("Multiply lambda on success",     this.lambdaStepDown, 5);
		gd.addNumericField("Threshold RMS to exit LMA",      this.thresholdFinish, 7,9,"pix");
		gd.addNumericField("Multiply lambda on failure",     this.lambdaStepUp, 5);
		gd.addNumericField("Threshold lambda to fail",       this.maxLambda, 5);
		gd.addNumericField("Maximal number of iterations",   this.numIterations, 0);

		gd.addCheckbox("Dialog after each iteration step",   this.stopEachStep);
		gd.addCheckbox("Dialog after each iteration series", this.stopEachSeries);
		gd.addCheckbox("Dialog after each failure",          this.stopOnFailure);
		gd.addCheckbox("Show modified parameters",           this.showParams);
		gd.addCheckbox("Show debug images before correction",this.showThisImages);
		gd.addCheckbox("Show debug images after correction", this.showNextImages);
		gd.addMessage("Done will save the current (not new!) state and exit, Continue will proceed according to LMA");
		gd.enableYesNoCancel("Continue", "Done");
		WindowTools.addScrollBars(gd);

	    gd.showDialog();
	    if (gd.wasCanceled()) {
	    	this.saveSeries=false;
	    	return false;
	    }
		this.lambda=                 gd.getNextNumber();
		this.lambdaStepDown=         gd.getNextNumber();
		this.thresholdFinish=        gd.getNextNumber();
		this.lambdaStepUp=           gd.getNextNumber();
		this.maxLambda=              gd.getNextNumber();
		this.numIterations=    (int) gd.getNextNumber();
		this.stopEachStep=           gd.getNextBoolean();
		this.stopEachSeries=         gd.getNextBoolean();
		this.stopOnFailure=          gd.getNextBoolean();
		this.showParams=             gd.getNextBoolean();
		this.showThisImages=         gd.getNextBoolean();
		this.showNextImages=         gd.getNextBoolean();
	    this.saveSeries=true;
	    return gd.wasOKed();
    }

    public boolean modifyGrid(
    		DistortionCalibrationData distortionCalibrationData,
			int threadsMax,
			boolean updateStatus){
    	if (fittingStrategy==null) {
    		String msg="Fitting strategy does not exist, exiting";
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
    	}
    	if (distortionCalibrationData.sensorMasks==null){
    		String msg="Sensor mask(s) are not defined";
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
    	}
    	if (distortionCalibrationData.eyesisCameraParameters==null){
    		String msg="Eyesis camera parameters (and sensor dimensions) are not defined";
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
    	}
//    	if (! selectGridEnhanceParameters()) return false;
//    	int series=refineParameters.showDialog("Select Grid Tuning Parameters", 0xdc0, fittingStrategy.getNumSeries());
//    	int series=refineParameters.showDialog("Select Grid Tuning Parameters", 0x31cc0, (this.seriesNumber>=0)?this.seriesNumber:0); // 0x1dco with show result, but we can not show it easily
// todo - add and implement 0x10000 to show just one individual image
//    	int series=refineParameters.showDialog("Select Grid Tuning Parameters", 0x21cc0, (this.seriesNumber>=0)?this.seriesNumber:0); // 0x1dco with show result, but we can not show it easily
    	int series=refineParameters.showDialog(
    			"Select Grid Tuning Parameters",
    			0x61000,
    			((this.seriesNumber>=0)?this.seriesNumber:0),
    			null, // averageRGB - only for target flat-field correction
    			false); // no difference for LWIR sensors?


    	if (series<0) return false;
    	this.seriesNumber=series;

		int filter=this.filterForTargetGeometry;
		if (this.askFilter) filter=selectFilter(filter);
    	initFittingSeries(true,filter,this.seriesNumber); // first step in series
//    	initFittingSeries(true,this.filterForTargetGeometry, this.seriesNumber); // first step in series
    	this.currentfX=calculateFxAndJacobian(this.currentVector, false);
    	//        	this.currentRMS= calcError(calcYminusFx(this.currentfX));
    	if (this.debugLevel>2) {
    		System.out.println("this.currentVector");
    		for (int i=0;i<this.currentVector.length;i++){
    			System.out.println(i+": "+ this.currentVector[i]);
    		}
    	}
    	if (this.showThisImages) showDiff (this.currentfX, "residual-series-"+this.seriesNumber);
    	if (this.refineParameters.resetVariations) {
    		this.patternParameters.resetStationZCorr();
    	}
    	double [][][] correctionCombo= calculateGridXYZCorr3D(
    			this.refineParameters.variationPenalty,
    			this.refineParameters.fixXY,
                this.refineParameters.useVariations?(this.fittingStrategy.zGroups[this.seriesNumber]):null, //stationGroups,
				this.refineParameters.grid3DCorrection,
				this.refineParameters.rotateCorrection,
				this.refineParameters.grid3DMaximalZCorr, //20.0,
				this.refineParameters.noFallBack,
				this.refineParameters.targetShowPerImage,
				threadsMax,
				updateStatus);
    	double [][] gridXYZCorr=correctionCombo[0];
		double [][] gridZCorr3d =correctionCombo[1];
		double [][] gridZCorr3dWeight =correctionCombo[2];
		String [] titles={"X-correction(mm)","Y-correction(mm)","Z-correction","Weight"};
    	String [] titlesStations=new String [2*gridZCorr3d.length];
    	for (int i=0;i<gridZCorr3d.length;i++){
    		titlesStations[i]="Z_"+i;
    		titlesStations[i+gridZCorr3d.length]="W_"+i;
    	}
    	if (this.refineParameters.targetShowThisCorrection) {
    		if (this.debugLevel>1){
    			double [][] debugData=new double [2*gridZCorr3d.length][];
    	    	for (int i=0;i<gridZCorr3d.length;i++){
    	    		debugData[i]=gridZCorr3d[i];
    	    		debugData[i+gridZCorr3d.length]=gridZCorr3dWeight[i];
    	    	}
        		ShowDoubleFloatArrays.showArrays(debugData, getGridWidth(), getGridHeight(),  true, "Z corrections", titlesStations);
    		}
    	}


// TODO: make configurable and optional
		shrinkExtrapolateGridCorrection(
				gridXYZCorr, // dx,dy,dz, mask >0
				gridZCorr3d,
				getGridWidth(),
				1, //preShrink,
				5, // expand,
				3.0, //  sigma,
				2.0); //double ksigma

    	if (this.refineParameters.targetShowThisCorrection) {
    		ShowDoubleFloatArrays.showArrays(gridXYZCorr, getGridWidth(), getGridHeight(),  true, "Grid corrections", titles);
    		if (this.debugLevel>1){

    		}
    	}
    	if (!this.refineParameters.targetApplyCorrection) return false;
    	patternParameters.applyGridCorrection(gridXYZCorr, this.refineParameters.targetCorrectionScale);
    	patternParameters.applyZGridCorrection(gridZCorr3d, this.refineParameters.targetCorrectionScale);
    	return true;
    }

    public boolean modifyGrid0(DistortionCalibrationData distortionCalibrationData){
    	if (fittingStrategy==null) {
    		String msg="Fitting strategy does not exist, exiting";
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
    	}
    	if (distortionCalibrationData.sensorMasks==null){
    		String msg="Sensor mask(s) are not defined";
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
    	}
    	if (distortionCalibrationData.eyesisCameraParameters==null){
    		String msg="Eyesis camera parameters (and sensor dimensions) are not defined";
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
    	}
//    	if (! selectGridEnhanceParameters()) return false;
//    	int series=refineParameters.showDialog("Select Grid Tuning Parameters", 0xdc0, fittingStrategy.getNumSeries());
//    	int series=refineParameters.showDialog("Select Grid Tuning Parameters", 0x31cc0, (this.seriesNumber>=0)?this.seriesNumber:0); // 0x1dco with show result, but we can not show it easily
// todo - add and implement 0x10000 to show just one individual image
//    	int series=refineParameters.showDialog("Select Grid Tuning Parameters", 0x21cc0, (this.seriesNumber>=0)?this.seriesNumber:0); // 0x1dco with show result, but we can not show it easily
    	int series=refineParameters.showDialog(
    			"Select Grid Tuning Parameters",
    			0x61000,
    			((this.seriesNumber>=0)?this.seriesNumber:0),
    			null, // averageRGB - only for target flat-field correction
    			false); // no difference for LWIR?


    	if (series<0) return false;
    	this.seriesNumber=series;



		int filter=this.filterForTargetGeometry;
		if (this.askFilter) filter=selectFilter(filter);
    	initFittingSeries(true,filter,this.seriesNumber); // first step in series
//    	initFittingSeries(true,this.filterForTargetGeometry, this.seriesNumber); // first step in series
    	this.currentfX=calculateFxAndJacobian(this.currentVector, false);
    	//        	this.currentRMS= calcError(calcYminusFx(this.currentfX));
    	if (this.debugLevel>2) {
    		System.out.println("this.currentVector");
    		for (int i=0;i<this.currentVector.length;i++){
    			System.out.println(i+": "+ this.currentVector[i]);
    		}
    	}
    	if (this.showThisImages) showDiff (this.currentfX, "residual-series-"+this.seriesNumber);

    	double [][] gridXYZCorr=null;
		gridXYZCorr=	calculateGridXYZCorr3D(
//				distortionCalibrationData,
				this.refineParameters.grid3DCorrection,
				this.refineParameters.rotateCorrection,
				this.refineParameters.grid3DMaximalZCorr, //20.0,
				this.refineParameters.targetShowPerImage);

// TODO: make configurable and optional
		shrinkExtrapolateGridCorrection(
				gridXYZCorr, // dx,dy,dz, mask >0
				null,
				getGridWidth(),
				1, //preShrink,
				5, // expand,
				3.0, //  sigma,
				2.0); //double ksigma

    	String [] titles={"X-correction(mm)","Y-correction(mm)","Z-correction","Weight"};
    	if (this.refineParameters.targetShowThisCorrection) {
    		ShowDoubleFloatArrays.showArrays(gridXYZCorr, getGridWidth(), getGridHeight(),  true, "Grid corrections", titles);
    	}
    	if (!this.refineParameters.targetApplyCorrection) return false;
    	patternParameters.applyGridCorrection(gridXYZCorr, this.refineParameters.targetCorrectionScale);
    	return true;
    }

//    public boolean modifyPixelCorrectionOld(DistortionCalibrationData distortionCalibrationData){ // old removed

    public void resetSensorCorrection(){
    	this.pixelCorrection=null;
    	this.pathNames=null;
    }
    public void resetSensorCorrection(int sensorNum){
    	if ((this.pixelCorrection!=null) && (sensorNum<this.pixelCorrection.length) && (sensorNum>=0)) {
    		this.pixelCorrection[sensorNum]=null;
    		this.pathNames[sensorNum]=null;
    	}
    }
    public void initSensorCorrection(){
    	int numLayers=7;
    	int numChannels=this.fittingStrategy.distortionCalibrationData.getNumChannels(); // number of used channels
    	this.pixelCorrection=new double [numChannels][][];
    	this.pathNames=new String[numChannels];
    	double [][] masks=this.fittingStrategy.distortionCalibrationData.calculateSensorMasks();
    	for (int i=0;i<this.pixelCorrection.length;i++){
    		this.pixelCorrection[i]=new double [numLayers][];
    		this.pathNames[i]=null;
    		for (int n=0;n<numLayers;n++) this.pixelCorrection[i][n]=new double [masks[i].length];
    		for (int j=0;j<masks[i].length;j++) {
        		this.pixelCorrection[i][0][j]=0.0;
        		this.pixelCorrection[i][1][j]=0.0;
        		this.pixelCorrection[i][2][j]=masks[i][j];
    			this.pixelCorrection[i][3][j]=1.0;
    			this.pixelCorrection[i][4][j]=1.0;
    			this.pixelCorrection[i][5][j]=1.0;
    		}
    	}
    }

    public void initSensorCorrection(int sensorNum){
    	int numLayers=7;
    	if ((this.pixelCorrection!=null) && (sensorNum<this.pixelCorrection.length) && (sensorNum>=0)) {
    		this.pixelCorrection[sensorNum]=null;
    		this.pathNames[sensorNum]=null;
    	}
    	double [] mask = this.fittingStrategy.distortionCalibrationData.calculateSensorMasks(sensorNum);
    	this.pixelCorrection[sensorNum]=new double [numLayers][];
    	this.pathNames[sensorNum]=null;
    	for (int n=0;n<numLayers;n++) this.pixelCorrection[sensorNum][n]=new double [mask.length];
    	for (int j=0;j<mask.length;j++) {
    		this.pixelCorrection[sensorNum][0][j]=0.0;
    		this.pixelCorrection[sensorNum][1][j]=0.0;
    		this.pixelCorrection[sensorNum][2][j]=mask[j];
    		this.pixelCorrection[sensorNum][3][j]=1.0;
    		this.pixelCorrection[sensorNum][4][j]=1.0;
    		this.pixelCorrection[sensorNum][5][j]=1.0;
    	}
    }


    /*
     * Adds new correction to the current one with the result to the new one. If update, the old arrays are also modified/created
     */
    public boolean applySensorCorrection(
    		RefineParameters refineParameters,
//    		boolean update,
//    		boolean updateFlatField,
//    		double scale,
    		double [][][] sensorXYCorr,
    		DistortionCalibrationData distortionCalibrationData){
		int numLayers=6;
///		int decimate=distortionCalibrationData.eyesisCameraParameters.decimateMasks;
///		int width= distortionCalibrationData.eyesisCameraParameters.sensorWidth;
///		int height=distortionCalibrationData.eyesisCameraParameters.sensorHeight;
///    	if ((this.pixelCorrection!=null) && (this.pixelCorrectionDecimation!=decimate)){
///    		IJ.showMessage("Error","Can not apply correction as the current correction and the new one have different decimations");
///    		return false;
///    	}
//    	if ((this.pixelCorrection==null) && !update && !updateFlatField) return true;
//    	if (update){
///    		this.pixelCorrectionDecimation=decimate;
///    		this.pixelCorrectionWidth=width;
///    		this.pixelCorrectionHeight=height;
//    	}
        if (this.pixelCorrection==null) {
        	if (this.debugLevel>1) System.out.println("Initializing pixelCorrection array...");
        	this.pixelCorrection=new double [sensorXYCorr.length][][];
        	this.pathNames=new String[sensorXYCorr.length];
        	for (int i=0;i<this.pixelCorrection.length;i++){
        		this.pixelCorrection[i]=null;
        		this.pathNames[i]=null;
        	}
        }

        if (this.pixelCorrection.length<sensorXYCorr.length){ // OK to update even if !update
        	if (this.debugLevel>1) System.out.println("Increasing number of sensors in pixelCorrection array");
        	double [][][] tmp=new double[sensorXYCorr.length][][];
        	String [] tmpPaths=new String[sensorXYCorr.length];
        	for (int i=0;i<tmp.length;i++){
        		if (i<this.pixelCorrection.length){
        			tmp[i]=this.pixelCorrection[i];
        			tmpPaths[i]=this.pathNames[i];
        		}
        		else {
        			tmp[i]=null;
        			tmpPaths[i]=null;
        		}
        	}
        	this.pixelCorrection=tmp;
        	this.pathNames=tmpPaths;
        }
        for (int i=0;i<sensorXYCorr.length;i++) if (sensorXYCorr[i]!=null){
//			boolean small_sensor = fittingStrategy.distortionCalibrationData.isSmallSensor(i);
			boolean small_sensor = fittingStrategy.distortionCalibrationData.getSmallSensors()[i];
			RefineParameters rp = small_sensor ? refineParameters.refineParametersSmall : refineParameters;
    		boolean update =          rp.applyCorrection;
    		boolean updateFlatField = rp.applyFlatField;
    		double scale =            rp.correctionScale;


        	boolean in6=sensorXYCorr[i].length==6; // was - 7
        	int indxR=in6?3:4;
        	int indxG=in6?4:5;
        	int indxB=in6?5:6;
        	double [] sensorMask=in6?((fittingStrategy.distortionCalibrationData.sensorMasks==null)?null:fittingStrategy.distortionCalibrationData.sensorMasks[i]):sensorXYCorr[i][2];
        	if (this.pixelCorrection[i]==null) {
        		if (update || updateFlatField) {
        			this.pixelCorrection[i]=new double [numLayers][];
        			this.pixelCorrection[i][0]=sensorXYCorr[i][0];
        			this.pixelCorrection[i][1]=sensorXYCorr[i][1];
        			if (sensorMask!=null){
        				this.pixelCorrection[i][2]=sensorMask;
        			} else {
        				this.pixelCorrection[i][2]= new double[this.pixelCorrection[i][0].length];
    					for (int j=0;j<this.pixelCorrection[i][2].length;j++) this.pixelCorrection[i][2][j]=1.0;
        			}
        			if (sensorXYCorr[i].length>=7){
        				this.pixelCorrection[i][3]=sensorXYCorr[i][indxR];
        				this.pixelCorrection[i][4]=sensorXYCorr[i][indxG];
        				this.pixelCorrection[i][5]=sensorXYCorr[i][indxB];
        			} else {
        				for (int n=3;n<numLayers;n++){
//        					this.pixelCorrection[i][n]=new double[this.pixelCorrection[0].length];
        					this.pixelCorrection[i][n]=new double[this.pixelCorrection[i][0].length]; // number of pixels
        					for (int j=0;j<this.pixelCorrection[i][0].length;j++) {
        						if ((i >= pixelCorrection.length) || (n >= pixelCorrection[i].length) || (j >= pixelCorrection[i][n].length)){
        							System.out.println("i="+i+", n="+n+", j="+j);
        							continue;
        						}
        						this.pixelCorrection[i][n][j]=1.0; // java.lang.ArrayIndexOutOfBoundsException: Index 6 out of bounds for length 6
        					}
        				}

        			}
        		}
        	} else  {
        		for (int j=0;j<sensorXYCorr[i][0].length;j++){
        			// removed - now it is already done
///        			sensorXYCorr[i][0][j]=this.pixelCorrection[i][0][j]+scale*sensorXYCorr[i][0][j];
///        			sensorXYCorr[i][1][j]=this.pixelCorrection[i][1][j]+scale*sensorXYCorr[i][1][j];
        			if (scale==1.0) { // recovering from Double.NaN in old values - still do not know where it came from in the first place
        			} else {
        				if (!in6){
        					sensorXYCorr[i][2][j]=this.pixelCorrection[i][2][j]+scale*(sensorXYCorr[i][2][j]-this.pixelCorrection[i][2][j]);
        				}
            			sensorXYCorr[i][indxR][j]=this.pixelCorrection[i][3][j]+scale*(sensorXYCorr[i][indxR][j]-this.pixelCorrection[i][3][j]);
            			sensorXYCorr[i][indxG][j]=this.pixelCorrection[i][4][j]+scale*(sensorXYCorr[i][indxG][j]-this.pixelCorrection[i][4][j]);
            			sensorXYCorr[i][indxB][j]=this.pixelCorrection[i][5][j]+scale*(sensorXYCorr[i][indxB][j]-this.pixelCorrection[i][5][j]);

        			}
        		}
        		if (update){
        			this.pixelCorrection[i][0]=sensorXYCorr[i][0];
        			this.pixelCorrection[i][1]=sensorXYCorr[i][1];
        		}
        		if (updateFlatField){
        			if (!in6){
        				this.pixelCorrection[i][2]=sensorXYCorr[i][2];
        			}
        			this.pixelCorrection[i][3]=sensorXYCorr[i][indxR];
        			this.pixelCorrection[i][4]=sensorXYCorr[i][indxG];
        			this.pixelCorrection[i][5]=sensorXYCorr[i][indxB];
        		}
        	}
        }
        return true;

    }
    public String getSensorPath(int numSensor){ //<0 - first available;
    	if ((this.pathNames == null) || (numSensor>=this.pathNames.length)) return null;
    	if (numSensor>=0) return this.pathNames[numSensor];
    	for (int i=0;i<this.pathNames.length;i++) if ((this.pathNames[i]!=null) && (this.pathNames[i].length()>0)) return this.pathNames[i];
    	return null;
    }
    public void saveDistortionAsImageStack(
    		DistortionCalibrationData distortionCalibrationData, // null OK
    		CamerasInterface camerasInterface, // to save channel map
    		String title,
    		String path,
    		boolean emptyOK){
    	int indexPeriod=path.indexOf('.',path.lastIndexOf(Prefs.getFileSeparator()));
    	int indexSuffix=indexPeriod;
    	String digits="0123456789";
    	for (int i=1;i<=2;i++) if (digits.indexOf(path.charAt(indexSuffix-1))>=0) indexSuffix--; // remove 1 or 2 digits before period
    	boolean hadSuffix= (path.charAt(indexSuffix-1)=='-');
    	int numSubCameras=fittingStrategy.distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[0].length;
    	for (int chNum=0;chNum<numSubCameras;chNum++) if (emptyOK  || ((this.pixelCorrection!=null) && (chNum<this.pixelCorrection.length) && (this.pixelCorrection[chNum]!=null)))  {
    		String channelPath= (hadSuffix?path.substring(0,indexSuffix):(path.substring(0,indexPeriod)+"-"))+
    		String.format("%02d",chNum)+path.substring(indexPeriod);
    		saveDistortionAsImageStack(
    				distortionCalibrationData,
    				camerasInterface, // to save channel map
    				title,
    				channelPath,
    				chNum,
    				emptyOK);
    	}
    }

    public ImagePlus saveDistortionAsImageStack(
    		DistortionCalibrationData distortionCalibrationData, // null OK
    		CamerasInterface camerasInterface, // to save channel map
    		String title,
    		String path,
    		int numSensor,
    		boolean emptyOK){
    	ImagePlus imp=getDistortionAsImageStack(
    			distortionCalibrationData,
    			camerasInterface, // to save channel map
    			title,
    			numSensor,
    			emptyOK);
    	if (imp==null) return null;
    	boolean realData= (this.pixelCorrection!=null) && (this.pixelCorrection[numSensor]!=null);
    	FileSaver fs=new FileSaver(imp);
    	String msg="Saving "+(realData?"":"EMPTY")+" sensor distortions to "+path;
    	if (updateStatus) IJ.showStatus(msg);
    	if (this.debugLevel>0) System.out.println(msg);
    	fs.saveAsTiffStack(path);
    	if (this.pathNames==null){
    		this.pathNames=new String[this.fittingStrategy.distortionCalibrationData.getNumChannels()];
    		for (int i=0;i<this.pathNames.length;i++) this.pathNames[i]=null;
    	}
    	this.pathNames[numSensor]=path;
    	return imp;
    }

//  /    	int numChannels=this.fittingStrategy.distortionCalibrationData.getNumChannels(); // number of used channels
// TODO: Currently saves data from Station 0
    public ImagePlus getDistortionAsImageStack(
    		DistortionCalibrationData distortionCalibrationData, // null OK - will use old way from fittingStrategy
    		CamerasInterface camerasInterface, // to save channel map
    		String title,
    		int numSensor,
    		boolean emptyOK){
    	if (distortionCalibrationData == null) {
    		distortionCalibrationData = this.fittingStrategy.distortionCalibrationData;
    	}
    	int stationNumber=0;
    	String [] titles={"X-corr","Y-corr","mask","R-vign","G-vign","B-vign"};
    	double [][] pixelCorr=null;
    	if (!emptyOK &&((this.pixelCorrection==null) ||
    			(numSensor<0) ||
    		(numSensor>=this.pixelCorrection.length) ||
    		(this.pixelCorrection[numSensor]==null)))
    			{
    		String msg="Sensor correction data is not available";
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
    	}
    	if ((this.pixelCorrection!=null) && (numSensor>=0) && (numSensor<this.pixelCorrection.length))
    		pixelCorr=this.pixelCorrection[numSensor];

    	int width =  distortionCalibrationData.eyesisCameraParameters.getSensorWidth(numSensor) /
    			distortionCalibrationData.eyesisCameraParameters.getDecimateMasks(numSensor);
    	int height = distortionCalibrationData.eyesisCameraParameters.getSensorHeight(numSensor) /
    			distortionCalibrationData.eyesisCameraParameters.getDecimateMasks(numSensor);

//    	int length=this.pixelCorrection[numSensor][0].length; // should be == width*height
    	int length=width*height;

    	float [][]pixels=new float [titles.length][length]; // dx, dy, sensor mask,v-r,v-g,v-b
    	// assuming all sensors have the same dimension
    	double [] mask=null;
    	if (distortionCalibrationData.sensorMasks == null) {
        	if (this.debugLevel>0) {
        		System.out.println("sensorMasks are null, calculating");
        	}

    		distortionCalibrationData.calculateSensorMasks();
    	}
    	
    	if (distortionCalibrationData.sensorMasks.length<=numSensor) return null; // no data (make if distortionCalibrationData.sensorMasks==null
    	if ((distortionCalibrationData.sensorMasks!=null) &&
    			(distortionCalibrationData.sensorMasks[numSensor]!=null)){
    		mask=distortionCalibrationData.sensorMasks[numSensor];
    	}

    	for (int index=0;index<length;index++){
    		if (pixelCorr==null){
        		pixels[0][index]=  0.0f;
        		pixels[1][index]=  0.0f;
        		for (int n=3;n<pixels.length;n++) pixels[n][index]= 1.0f; // normalize?
    		} else {
    			pixels[0][index]=  (float) pixelCorr[0][index];
    			pixels[1][index]=  (float) pixelCorr[1][index];
        		for (int n=3;n<pixels.length;n++) {
        			if (pixelCorr[n] != null) {
        				if ((n>=pixels.length) || (index >= pixels[n].length) || (n>=pixelCorr.length) || (index >= pixelCorr[n].length)) {
        					System.out.println (" Bug: n="+n+", index="+index);
        					continue;
        				}
        				pixels[n][index]= (float) pixelCorr[n][index]; // java.lang.NullPointerException
        			}
        		}
    		}
    		// get sensor mask here
    		pixels[2][index]=  (mask==null)? 1.0f:((float) mask[index]);

    	}
    	ImagePlus imp=null;
  		ImageStack stack=new ImageStack(width,height);
   		for (int n=0;n<pixels.length;n++)  stack.addSlice(titles[n],    pixels[n]);
   		imp = new ImagePlus(title, stack);
        // set properties sufficient to un-apply distortions to the image
   		// First - corrections
    	EyesisSubCameraParameters subCam=distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[stationNumber][numSensor];
    	subCam.updateCartesian(); // recalculate other parameters
    	double entrancePupilForward=distortionCalibrationData.eyesisCameraParameters.entrancePupilForward[stationNumber];
    	imp.setProperty("VERSION",  "1.0");
    	imp.setProperty("comment_arrays",  "Array corrections from acquired image to radially distorted, in pixels");
    	imp.setProperty("arraysSet",  ""+(pixelCorr!=null)); // per-pixel arrays are not set, using 0.0
    	imp.setProperty("maskSet",     ""+(mask!=null)); // per-pixel masks is not set, using 1.0
    	imp.setProperty("pixelCorrectionWidth",  ""+distortionCalibrationData.eyesisCameraParameters.getSensorWidth(numSensor)); // this.pixelCorrectionWidth);
    	imp.setProperty("pixelCorrectionHeight", ""+distortionCalibrationData.eyesisCameraParameters.getSensorHeight(numSensor));
    	imp.setProperty("pixelCorrectionDecimation", ""+distortionCalibrationData.eyesisCameraParameters.getDecimateMasks(numSensor));
    	imp.setProperty("comment_decimation", "when decimation use integer divide to find the index, corection values are in non-decimated pixels");
    	imp.setProperty("distortion_formula",  "(normalized by distortionRadius in mm) Rdist/R=A8*R^7+A7*R^6+A6*R^5+A5*R^4+A*R^3+B*R^2+C*R+(1-A6-A7-A6-A5-A-B-C)");
    	imp.setProperty("distortionRadius", ""+subCam.distortionRadius);
    	imp.setProperty("distortionRadius_unuts", "mm");
    	imp.setProperty("focalLength", ""+subCam.focalLength);
    	imp.setProperty("focalLength_units", "mm");
    	imp.setProperty("pixelSize", ""+subCam.pixelSize);
    	imp.setProperty("lineTime",  ""+subCam.lineTime);
    	imp.setProperty("pixelSize_units", "um");
    	imp.setProperty("distortionA8", ""+subCam.distortionA8);
    	imp.setProperty("distortionA7", ""+subCam.distortionA7);
    	imp.setProperty("distortionA6", ""+subCam.distortionA6);
    	imp.setProperty("distortionA5", ""+subCam.distortionA5);
    	imp.setProperty("distortionA", ""+subCam.distortionA);
    	imp.setProperty("distortionB", ""+subCam.distortionB);
    	imp.setProperty("distortionC", ""+subCam.distortionC);
    	imp.setProperty("comment_px0_py0", "lens center on the sensor, in pixels");
    	imp.setProperty("px0", ""+subCam.px0);
    	imp.setProperty("py0", ""+subCam.py0);
    	imp.setProperty("comment_azimuth", "lens center azimuth, CW from top, degrees");
    	imp.setProperty("height",  ""+subCam.height);
    	imp.setProperty("comment_elevation", "lens elevation from horizontal, positive - above horizon, degrees");
    	imp.setProperty("elevation",  ""+subCam.theta);
    	imp.setProperty("comment_roll", "lens rotation around the lens axis. Positive - CW looking to the target, degrees");
    	imp.setProperty("roll",  ""+subCam.psi);

    	imp.setProperty("comment_cartesian", "Use cartesian coordinates for the sensor in the camera CS (forward, right,aheading), instead of (radius, azimuth, heading)");
    	imp.setProperty("cartesian",  ""+subCam.cartesian);
// cartesian parameters
    	imp.setProperty("comment_forward", "lens forward (towards target) displacement in the camera CS");
    	imp.setProperty("forward",  ""+subCam.forward);
    	imp.setProperty("comment_right", "lens right (looking towards target) displacement in the camera CS");
    	imp.setProperty("right",  ""+subCam.right);
    	imp.setProperty("comment_aheading", "lens axis horizontal direction, degrees. Positive - CW from the target (looking from top)");
    	imp.setProperty("aheading",  ""+subCam.heading);
// cylindrical parameters
    	imp.setProperty("azimuth", ""+subCam.azimuth);
    	imp.setProperty("comment_radius", "lens center distance from the camera vertical axis, mm");
    	imp.setProperty("radius",  ""+subCam.radius);
    	imp.setProperty("comment_height", "lens center vertical position from the head center, mm");
    	imp.setProperty("comment_heading", "lens heading - added to azimuth");
    	imp.setProperty("heading",  ""+subCam.phi);

    	imp.setProperty("comment_channel", "number of the sensor (channel) in the camera");
    	imp.setProperty("channel",  ""+numSensor);
    	imp.setProperty("comment_subcamera", "number of the subcamera with individual IP, starting with 0");
    	if (camerasInterface != null) {
    		subCam.subcamera =   camerasInterface.getSubCamera(numSensor);
    		subCam.sensor_port = camerasInterface.getSensorPort(numSensor);
    		subCam.subchannel =  camerasInterface.getSubChannel(numSensor);
    	}

    	imp.setProperty("subcamera",  ""+subCam.subcamera);
		imp.setProperty("sensor_port",""+subCam.sensor_port);

    	imp.setProperty("comment_subchannel", "number of the sensor port on a subcamera (0..2)");
    	imp.setProperty("subchannel",  ""+subCam.subchannel);
    	imp.setProperty("comment_entrancePupilForward",  "entrance pupil distance from the azimuth/radius/height, outwards in mm");
    	imp.setProperty("entrancePupilForward",  ""+entrancePupilForward); // currently global, decoders will use per-sensor
       	imp.setProperty("comment_defects", "Sensor hot/cold pixels list as x:y:difference");

       	imp.setProperty("comment_lensDistortionModel", "Integer specifying lens distrotion model (0 - radial)");
       	imp.setProperty("lensDistortionModel", ""+subCam.lensDistortionModel);

		for (int i=0;i<subCam.r_xy.length;i++){
			imp.setProperty("r_xy_"+i+"_x",subCam.r_xy[i][0]+"");
			imp.setProperty("r_xy_"+i+"_y",subCam.r_xy[i][1]+"");
		}
		for (int i=0;i<subCam.r_od.length;i++){
			imp.setProperty("r_od_"+i+"_o",subCam.r_od[i][0]+"");
			imp.setProperty("r_od_"+i+"_d",subCam.r_od[i][1]+"");
		}
       	if (subCam.defectsXY!=null){
    		StringBuffer sb = new StringBuffer();
    		for (int i=0;i<subCam.defectsXY.length;i++){
    			if (sb.length()>0) sb.append(" ");
    			sb.append(subCam.defectsXY[i][0]+":"+subCam.defectsXY[i][1]+":"+subCam.defectsDiff[i]);
    		}
    		imp.setProperty("defects", sb.toString());
//       	} else {
//    		imp.setProperty("defects", null);
       	}

    	//camerasInterface, numSensor
    	(new JP46_Reader_camera(false)).encodeProperiesToInfo(imp);
    	imp.getProcessor().resetMinAndMax();
    	return imp;
    }


    public void setDistortionFromImageStack(
    		DistortionCalibrationData distortionCalibrationData, // null OK
    		EyesisCameraParameters eyesisCameraParameters, // null OK
    		String path,
    		boolean overwriteExtrinsic,
    		boolean overwriteDistortion){
    	int indexPeriod=path.indexOf('.',path.lastIndexOf(Prefs.getFileSeparator()));
    	if (eyesisCameraParameters == null) {
    		eyesisCameraParameters = fittingStrategy.distortionCalibrationData.eyesisCameraParameters;
    	}
    	if (distortionCalibrationData == null) {
    		distortionCalibrationData = fittingStrategy.distortionCalibrationData;
    	}
    	int numSubCameras=distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[0].length;
    	for (int chNum=0;chNum<numSubCameras;chNum++){
    		String channelPath=path.substring(0,indexPeriod-2)+String.format("%02d",chNum)+path.substring(indexPeriod);
    		try { // disable here for now
    			setDistortionFromImageStack(
    					distortionCalibrationData,
    					channelPath,
    					chNum,
    					false,
    					overwriteExtrinsic,
    					overwriteDistortion);
    		} catch (Exception e) {
    			System.out.println("setDistortionFromImageStack(): " + e.toString());
    			e.printStackTrace();
    		}
    	}
    }

    public void setDistortionFromImageStack(
    		DistortionCalibrationData distortionCalibrationData,
    		String path,
    		int numSensor,
    		boolean reportProblems,
    		boolean overwriteExtrinsic,
    		boolean overwriteDistortion){
    	Opener opener=new Opener();
    	ImagePlus imp=opener.openImage("", path);
    	if (imp==null) {
    		if (!reportProblems) return;
    		String msg="Failed to read sensor calibration data file "+path;
    		IJ.showMessage("Error",msg);
    		System.out.println(msg);
    		throw new IllegalArgumentException (msg);
    	}
    	if (this.debugLevel>0) System.out.println("Read "+path+" as a sensor calibration data");
    	(new JP46_Reader_camera(false)).decodeProperiesFromInfo(imp);
    	setDistortionFromImageStack(distortionCalibrationData,
    			imp,
    			numSensor,
    			overwriteExtrinsic,
    			overwriteDistortion);
    	this.pathNames[numSensor]=path;
    }

    //TODO: look more after testing. Currently all station parameters are set from the sensor images, may be minor differences
    public void setDistortionFromImageStack(
    		DistortionCalibrationData distortionCalibrationData,
    		ImagePlus imp,
    		int numSensor,
    		boolean overwriteExtrinsic,
    		boolean overwriteDistortion){
    	if (distortionCalibrationData == null) {
    		distortionCalibrationData = this.fittingStrategy.distortionCalibrationData;
    	}
//    	int corrX=0,corrY=1,
    	int corrMask=2;
    	if (numSensor<0) {
    		System.out.println("setDistortionFromImageStack(): Tried to read negative channel");
    		return;
    	}
//		System.out.println("setDistortionFromImageStack(): processing channel channel "+numSensor);
    	if (imp == null){
    		String msg="Distortions image is null";
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
    	}
        String [] requiredProperties={
        		"pixelCorrectionWidth",
        		"pixelCorrectionHeight",
        		"pixelCorrectionDecimation",
        		"distortionRadius",
        		"focalLength",
        		"pixelSize",
//        		"distortionA8",
//        		"distortionA7",
//        		"distortionA6",
        		"distortionA5",
        		"distortionA",
        		"distortionB",
        		"distortionC",
        		"px0",
        		"py0"};
        for (int i=0; i<requiredProperties.length;i++) if (imp.getProperty(requiredProperties[i])==null){
    		String msg="Required property "+requiredProperties[i]+" is not defined in "+imp.getTitle();
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
        }

    	if (imp.getStackSize()<3){
    		String msg="Expecting >=3 slices, got "+imp.getStackSize();
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
    	}

		ImageStack stack = imp.getStack();
		float [][] pixels =new float[stack.getSize()][];
    	for (int i=0;i<pixels.length;i++) pixels[i]= (float[]) stack.getPixels(i+1);


        int numSubCameras=distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[0].length;
        if (numSensor>=numSubCameras){
    		String msg="Loaded calibration channel number "+numSensor+"is higher than maximal in the system "+(numSubCameras-1);
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
        }
//		System.out.println("setDistortionFromImageStack(): processing channel channel "+numSensor);

    	EyesisSubCameraParameters subCam;
    	EyesisCameraParameters cam=     distortionCalibrationData.eyesisCameraParameters;
        if ((distortionCalibrationData!=null) && (distortionCalibrationData.eyesisCameraParameters!=null)){
        	// Now it is the same
///    		distortionCalibrationData.eyesisCameraParameters.decimateMasks=this.pixelCorrectionDecimation;
///    		distortionCalibrationData.eyesisCameraParameters.sensorWidth=  this.pixelCorrectionWidth;
///    		distortionCalibrationData.eyesisCameraParameters.sensorHeight=this.pixelCorrectionHeight;
        }
        for (int stationNumber=0; stationNumber < distortionCalibrationData.eyesisCameraParameters.numStations; stationNumber++){
        	subCam=distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[stationNumber][numSensor];
        	distortionCalibrationData.eyesisCameraParameters.setSensorWidth(  numSensor, Integer.parseInt ((String) imp.getProperty("pixelCorrectionWidth")));
        	distortionCalibrationData.eyesisCameraParameters.setSensorHeight( numSensor, Integer.parseInt ((String) imp.getProperty("pixelCorrectionHeight")));
        	distortionCalibrationData.eyesisCameraParameters.setDecimateMasks(numSensor, Integer.parseInt ((String) imp.getProperty("pixelCorrectionDecimation")));

        	subCam.distortionRadius=        Double.parseDouble((String) imp.getProperty("distortionRadius"));
        	subCam.focalLength=             Double.parseDouble((String) imp.getProperty("focalLength"));
        	subCam.pixelSize=               Double.parseDouble((String) imp.getProperty("pixelSize"));
        	if (imp.getProperty("lineTime") != null) {
            	subCam.pixelSize=               Double.parseDouble((String) imp.getProperty("lineTime"));
        	} else { // fix older saved files
        		if (subCam.pixelSize < 5.0) {
        			subCam.pixelSize=3.638E-5;
	        	} else if (distortionCalibrationData.eyesisCameraParameters.getSensorWidth(numSensor) == 640){ // Boson
	        		subCam.lineTime = 2.7778e-05; // 12um pixel, Boson
	        	} else {
	        		subCam.lineTime = 7.8e-05; // 12um pixel, Lepton (may7 be wrong)
	        	}
        	}
        	if (imp.getProperty("distortionA8")!=null) {
        		subCam.distortionA8=            Double.parseDouble((String) imp.getProperty("distortionA8"));
        	} else subCam.distortionA8=0.0;
        	if (imp.getProperty("distortionA7")!=null) {
        		subCam.distortionA7=            Double.parseDouble((String) imp.getProperty("distortionA7"));
        	} else subCam.distortionA7=0.0;
        	if (imp.getProperty("distortionA6")!=null) {
        		subCam.distortionA6=            Double.parseDouble((String) imp.getProperty("distortionA6"));
        	} else subCam.distortionA6=0.0;
        	subCam.distortionA5=            Double.parseDouble((String) imp.getProperty("distortionA5"));
        	subCam.distortionA=             Double.parseDouble((String) imp.getProperty("distortionA"));
        	subCam.distortionB=             Double.parseDouble((String) imp.getProperty("distortionB"));
        	subCam.distortionC=             Double.parseDouble((String) imp.getProperty("distortionC"));
        	subCam.px0=                     Double.parseDouble((String) imp.getProperty("px0"));
        	subCam.py0=                     Double.parseDouble((String) imp.getProperty("py0"));
        	if (imp.getProperty("azimuth")  !=null) subCam.azimuth= Double.parseDouble((String) imp.getProperty("azimuth"));
        	if (imp.getProperty("radius")   !=null) subCam.radius=  Double.parseDouble((String) imp.getProperty("radius"));
        	if (imp.getProperty("height")   !=null) subCam.height=  Double.parseDouble((String) imp.getProperty("height"));
        	if (imp.getProperty("entrancePupilForward")!=null) cam.entrancePupilForward[stationNumber]= Double.parseDouble((String) imp.getProperty("entrancePupilForward"));
        	if (imp.getProperty("heading")  !=null) subCam.phi=     Double.parseDouble((String) imp.getProperty("heading"));
        	if (imp.getProperty("elevation")!=null) subCam.theta=   Double.parseDouble((String) imp.getProperty("elevation"));
        	if (imp.getProperty("roll")!=null) subCam.psi=          Double.parseDouble((String) imp.getProperty("roll"));

        	if (imp.getProperty("forward")  !=null) subCam.forward=  Double.parseDouble((String) imp.getProperty("forward"));
        	if (imp.getProperty("right")    !=null) subCam.right=    Double.parseDouble((String) imp.getProperty("right"));
        	if (imp.getProperty("aheading") !=null) subCam.heading= Double.parseDouble((String) imp.getProperty("aheading"));

        	if (imp.getProperty("cartesian") !=null) {
        		subCam.cartesian= Boolean.parseBoolean((String) imp.getProperty("cartesian"));
        		subCam.updateCartesian(); // recalculate other parameters (they may or may nort be provided
        	} else {
        		subCam.cartesian = false;
        	}

        	// Update intrinsic image parameters
        	this.lensDistortionParameters.pixelSize=subCam.pixelSize;
        	this.lensDistortionParameters.distortionRadius=subCam.distortionRadius;
        	if (imp.getProperty("defects")!=null) {
        		String sDefects=(String) imp.getProperty("defects");
        		String [] asDefects=sDefects.trim().split(" ");
        		subCam.defectsXY=new int [asDefects.length][2];
        		subCam.defectsDiff=new double [asDefects.length];
        		for (int i=0;i<asDefects.length;i++) {
        			String [] stDefect=asDefects[i].split(":");
        			subCam.defectsXY[i][0]=Integer.parseInt(stDefect[0]);
        			subCam.defectsXY[i][1]=Integer.parseInt(stDefect[1]);
        			subCam.defectsDiff[i]=Double.parseDouble(stDefect[2]);
        		}
        	} else {
        		subCam.defectsXY=null;
        		subCam.defectsDiff=null;
        	}
 // non-radial
        	if (imp.getProperty("lensDistortionModel")  !=null) subCam.lensDistortionModel= Integer.parseInt((String) imp.getProperty("lensDistortionModel"));
        	subCam.setDefaultNonRadial();
			for (int i=0;i<subCam.r_xy.length;i++) {
				if (imp.getProperty("r_xy_"+i+"_x")  !=null) subCam.r_xy[i][0]= Double.parseDouble((String) imp.getProperty("r_xy_"+i+"_x"));
				if (imp.getProperty("r_xy_"+i+"_y")  !=null) subCam.r_xy[i][1]= Double.parseDouble((String) imp.getProperty("r_xy_"+i+"_y"));
			}
			for (int i=0;i<subCam.r_od.length;i++) {
				if (imp.getProperty("r_od_"+i+"_o")  !=null) subCam.r_od[i][0]= Double.parseDouble((String) imp.getProperty("r_od_"+i+"_o"));
				if (imp.getProperty("r_od_"+i+"_d")  !=null) subCam.r_od[i][1]= Double.parseDouble((String) imp.getProperty("r_od_"+i+"_d"));
			}
			if (imp.getProperty("subcamera")   !=null) subCam.subcamera=   Integer.parseInt((String) imp.getProperty("subcamera"));
			if (imp.getProperty("sensor_port") !=null) subCam.sensor_port= Integer.parseInt((String) imp.getProperty("sensor_port"));
			if (imp.getProperty("subchannel")  !=null) subCam.subchannel=   Integer.parseInt((String) imp.getProperty("subchannel"));
        }


        for (int imgNum=0;imgNum < distortionCalibrationData.getNumImages();imgNum++){
        	int imageSubCam=  distortionCalibrationData.getImageSubcamera(imgNum);
        	int stationNumber=distortionCalibrationData.getImageStation(imgNum);
        	if (imageSubCam==numSensor){
        		// vector from the data we just set
        		double [] parVector=           distortionCalibrationData.eyesisCameraParameters.getParametersVector(stationNumber,imageSubCam);
        		if       (overwriteExtrinsic)  distortionCalibrationData.setSubcameraParameters(parVector,imgNum);
        		else  if (overwriteDistortion) distortionCalibrationData.setIntrinsicParameters(parVector,imgNum);
        	}
        }

        // now read the calibration data and mask
    	if (this.pixelCorrection==null) {
    		this.pixelCorrection=new double [numSubCameras][][];
    		this.pathNames=new String [numSubCameras];
    		for (int i=0;i<this.pixelCorrection.length;i++){
    			this.pixelCorrection[i]=null;
    			this.pathNames[i]=null;
    		}
    	}
        if (numSensor>=this.pixelCorrection.length){ // increase number of elements
        	double [][][] tmp=this.pixelCorrection.clone();
        	String [] tmpPaths=this.pathNames.clone();
        	this.pixelCorrection=new double[numSensor+1][][];
        	this.pathNames=new String[numSensor+1];
        	for (int i=0;i<this.pixelCorrection.length;i++)
        		if (i<tmp.length){
        			this.pixelCorrection[i]=tmp[i];
        			this.pathNames[i]=tmpPaths[i];
        		}else {
        			this.pixelCorrection[i]=null;
        			this.pathNames[i]=null;
        		}
        }
        int numLayers=6; //corr-x, corr-y,mask, ff-R, ff-G, ff-b
        if (numLayers<pixels.length) numLayers=pixels.length; // for the future?
//        this.pixelCorrection[numSensor]=new double [pixels.length] [pixels[0].length];
        this.pixelCorrection[numSensor]=new double [numLayers][pixels[0].length];
        for (int i= 0;i<this.pixelCorrection[numSensor][0].length;i++){
        	for (int n=0;n<pixels.length;n++)	this.pixelCorrection[numSensor][n][i]=pixels[n][i]; // mask will go to two places
        }
        if (pixels.length<numLayers){
            for (int i= 0;i<this.pixelCorrection[numSensor][0].length;i++){
            	for (int n=pixels.length;n<numLayers;n++)	this.pixelCorrection[numSensor][n][i]=1.0; // default ff if no data is available
            }
        }
        // now mask
        boolean defined=false;
		for (int i=0;i<pixels[2].length;i++) if ((pixels[2][i]!=0.0) && (pixels[2][i]!=1.0)){
			defined=true;
			break;
		}
//    	System.out.println("setDistortionFromImageStack(): defined="+defined );
		if (defined) {
	        if (distortionCalibrationData.sensorMasks==null) {
	        	distortionCalibrationData.sensorMasks=new double [numSubCameras][];
	        	for (int i=0;i<distortionCalibrationData.sensorMasks.length;i++)
	        		distortionCalibrationData.sensorMasks[i]=null;
//	        	System.out.println("setDistortionFromImageStack(): created distortionCalibrationData.sensorMasks["+numSubCameras+"] of null-s" );
	        }
	        if (numSensor>=distortionCalibrationData.sensorMasks.length){ // increase number of elements
	        	double [][] tmp=distortionCalibrationData.sensorMasks;
	        	distortionCalibrationData.sensorMasks=new double[numSensor+1][];
	        	for (int i=0;i<distortionCalibrationData.sensorMasks.length;i++)
	        		if (i<tmp.length)distortionCalibrationData.sensorMasks[i]=tmp[i];
	        		else distortionCalibrationData.sensorMasks[i]=null;
	        }
	        if (distortionCalibrationData.sensorMasks[numSensor]==null){
	        	distortionCalibrationData.sensorMasks[numSensor]=new double[pixels[corrMask].length];
//	        	System.out.println("setDistortionFromImageStack(): created distortionCalibrationData.sensorMasks["+numSensor+"] of ["+pixels[corrMask].length+"]" );
	        }
	        for (int i= 0;i<distortionCalibrationData.sensorMasks[numSensor].length;i++) // null pointer
	        	distortionCalibrationData.sensorMasks[numSensor][i]=pixels[corrMask][i];
		}
    }
    /**
     * Accumulate per-sensor grid R,G,B intensities using current sensor flat-field values
     * @param serNumber - fitting series number to select images (-1 - all enabled)
     * @param sensorMasks "pessimistic" masks to use only center (low-vignetting) part of each sensor (at least on the first runs?)
     * @param minContrast - minimal contrast to consider a node
     * @param threshold - not yet used - disregard grid nodes with low data - in the end
     * @param interpolate - interpolate sensor data
     * @param maskThresholdOcclusion suspect occlusion only if grid is missing in the area where sensor mask is above this threshold
     * @param expandOcclusion - shrink defined grid on image by this steps - to handle occlusion by rollers
     * @param fadeOcclusion - fade shrank occlusion border
     * @param ignoreSensorFlatField - ignorfe previously calculated sensors flat-field calibration
     * @return
     */

    public double [][][][] calculateGridFlatField(
    		int serNumber,
    		double [][] sensorMasks,
    		RefineParameters refineParameters){
   // TODO: add standard weight function used elsewhere.
    	int indexContrast=2;
    	boolean [] selectedImages=fittingStrategy.selectedImages(serNumber); // negative series number OK - will select all enabled
    	int gridHeight=this.patternParameters.gridGeometry.length;
    	int gridWidth=this.patternParameters.gridGeometry[0].length;

    	int maxChannel=0;
    	int numStations=this.patternParameters.getNumStations();
       	for (int numImg=0;numImg<fittingStrategy.distortionCalibrationData.gIP.length;numImg++) if (selectedImages[numImg]){
    		if (fittingStrategy.distortionCalibrationData.gIP[numImg].channel>maxChannel) maxChannel=fittingStrategy.distortionCalibrationData.gIP[numImg].channel;
    	}

    	double [][][][] sensorGrids=new double [numStations][maxChannel+1][][]; //{alpha, red,green, blue}
    	for (int ns=0;ns<sensorGrids.length;ns++) for (int n=0;n<sensorGrids[ns].length;n++) sensorGrids[ns][n]=null;
    	// For each sensor separately accumulate grid intensity using current sensor flat field calibration
    	for (int numImg=0;numImg<fittingStrategy.distortionCalibrationData.gIP.length;numImg++) if (selectedImages[numImg]) {
    		int channel=fittingStrategy.distortionCalibrationData.gIP[numImg].channel;
//    		boolean small_sensor = fittingStrategy.distortionCalibrationData.isSmallSensor(channel);
			boolean small_sensor = fittingStrategy.distortionCalibrationData.getSmallSensors()[channel];
    		RefineParameters rp = small_sensor ? refineParameters.refineParametersSmall : refineParameters;
    		double minContrast = rp.flatFieldMinimalContrast;
//    		double threshold = rp.flatFieldMinimalContrast;
    		boolean interpolate = rp.flatFieldUseInterpolate;
    		double maskThresholdOcclusion = rp.flatFieldMaskThresholdOcclusion;
    		int expandOcclusion = rp.flatFieldShrinkOcclusion;
    		double fadeOcclusion = rp.flatFieldFadeOcclusion;
    		boolean ignoreSensorFlatField = rp.flatFieldIgnoreSensorFlatField;

    		int station=fittingStrategy.distortionCalibrationData.gIP[numImg].getStationNumber();
    		if (sensorMasks[channel]==null) continue;
    		if (sensorGrids[station][channel]==null){ // null pointer
    			sensorGrids[station][channel]=new double [4][gridHeight*gridWidth]; //{alpha, red,green, blue}
    			for (int c=0;c<sensorGrids[station][channel].length;c++){
    				for (int i=0;i<sensorGrids[station][channel][0].length;i++) sensorGrids[station][channel][c][i]=0.0;
    			}
    		}
    		double [][] pixelsXY=fittingStrategy.distortionCalibrationData.gIP[numImg].pixelsXY;
    		if ((pixelsXY.length<1) || (pixelsXY[0].length<6)){
    			if (this.debugLevel>0) System.out.println("No flat-field data in image #"+numImg+
    					" - "+fittingStrategy.distortionCalibrationData.gIP[numImg].path+
    					" pixelsXY.length="+pixelsXY.length+
    					" pixelsXY[0].length="+((pixelsXY.length==0)?"nan": pixelsXY[0].length));
    			continue;
    		}
    		int [][]    pixelsUV=fittingStrategy.distortionCalibrationData.gIP[numImg].pixelsUV;
    		double [] defaultVector={0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    		// detect if there is any occlusion (i.e. by goniometer rollers)
    		boolean [] bMask=new boolean [gridHeight*gridWidth];
    		double [] mask=new double[bMask.length];
    		for (int i=0;i<bMask.length;i++){
    			bMask[i]=false;
    			mask[i]=0.0;
    		}
    		for (int i=0;i<pixelsXY.length;i++){
    			double [] xyzmrgb=patternParameters.getXYZM(
    					pixelsUV[i][0],
    					pixelsUV[i][1],
    					false,
    					station);
    			if (xyzmrgb!=null){
    	   			int index=patternParameters.getGridIndex(pixelsUV[i][0], pixelsUV[i][1]);
    	   			bMask[index]=(pixelsXY[i][indexContrast]>=minContrast);
    	   			mask[index]=interpolateMask (
    	   					channel,
        					sensorMasks[channel],
        					pixelsXY[i][0],
        					pixelsXY[i][1]);
    			}
    		}
    		boolean [] occlusionMask=new boolean[bMask.length];
    		for (int i=0;i<occlusionMask.length;i++){
    			occlusionMask[i]=false;
    		}
    		boolean occlusion=false;
    		for (int i=1;i<(gridHeight-1);i++){
    			for (int j=1;j<(gridWidth-1);j++){
    				int index=i*gridWidth+j;
    				if (bMask[index]){
    					if ((   !bMask[(i-1)*gridWidth+j] ||
    							!bMask[(i+1)*gridWidth+j] ||
    							!bMask[i*    gridWidth+j-1] ||
    							!bMask[i*    gridWidth+j+1]) &&
    							(mask[index]>=maskThresholdOcclusion)
    					){
    						occlusionMask[index]=true;
    						occlusion=true;
    					}
    				}
    			}
    		}
    		if (occlusion){
    			for (int n=0;n<expandOcclusion;n++){ // expand
    				boolean [] bMaskPrevious=occlusionMask.clone();
    				for (int i=1;i<(gridHeight-1);i++){
    					for (int j=1;j<(gridWidth-1);j++){
    						if (!occlusionMask[i*gridHeight+j]){
    							if (
    									bMaskPrevious[(i-1)*gridWidth+j] ||
    									bMaskPrevious[(i+1)*gridWidth+j] ||
    									bMaskPrevious[i*    gridWidth+j-1] ||
    									bMaskPrevious[i*    gridWidth+j+1]){
    								occlusionMask[i*gridWidth+j]=true;
    							}
    						}
    					}
    				}
    			}
        		double [] maskNonOccluded=new double [occlusionMask.length];
        		for (int i=0;i<maskNonOccluded.length;i++) maskNonOccluded[i]=occlusionMask[i]?0.0:1.0;
        		if (fadeOcclusion>0.0){
        			(new DoubleGaussianBlur() ).blurDouble(
        					maskNonOccluded,
        					gridWidth,
        					gridHeight,
        					fadeOcclusion,
        					fadeOcclusion,
        					0.01);
        		}
        		if (fadeOcclusion>=0.0 )for (int i=0;i<mask.length;i++){
    				double d=2.0*(maskNonOccluded[i]-0.5);
    				mask[i]*=(!occlusionMask[i] && (d>0))?(d*d):0.0;
    			}
    		}
    		for (int i=0;i<pixelsXY.length;i++){
    			double [] xyzmrgb=patternParameters.getXYZM(
    					pixelsUV[i][0],
    					pixelsUV[i][1],
    					false,
    					station);
    			if (xyzmrgb==null) continue; // out of grid
    			double [] vector=ignoreSensorFlatField?defaultVector:
    				((interpolate)?
    					interpolateCorrectionVector (
    							channel,
    							pixelsXY[i][0],
    							pixelsXY[i][1]):
 						getCorrectionVector (
 								channel,
 								pixelsXY[i][0],
 								pixelsXY[i][1]))	;
    			int index=patternParameters.getGridIndex(pixelsUV[i][0], pixelsUV[i][1]);
    			double weight=mask[index];

    			sensorGrids[station][channel][0][index]+=weight;
    			for (int c=0;c<3;c++)if (vector[c+3]>0.0){
    				sensorGrids[station][channel][c+1][index]+=weight*pixelsXY[i][c+3]/vector[c+3];
    			}
    		}
    	}
    	for (int station=0;station<sensorGrids.length;station++){
    		for (int channel=0;channel<sensorGrids[station].length; channel++) if (sensorGrids[station][channel]!=null){
//        		boolean small_sensor = fittingStrategy.distortionCalibrationData.isSmallSensor(channel);
				boolean small_sensor = fittingStrategy.distortionCalibrationData.getSmallSensors()[channel];
        		RefineParameters rp = small_sensor ? refineParameters.refineParametersSmall : refineParameters;
        		double threshold = rp.flatFieldMinimalContrast;
    			if (this.pixelCorrection[channel]==null) {
    				sensorGrids[station][channel]=null;
    			} else {
    				for (int i=0;i<sensorGrids[station][channel][0].length;i++){
    					if (sensorGrids[station][channel][0][i]<threshold) {
    						for (int j=0;j<sensorGrids[station][channel].length;j++) sensorGrids[station][channel][j][i]=0.0;
    					} else {
    						for (int j=1;j<sensorGrids[station][channel].length;j++) sensorGrids[station][channel][j][i]/=sensorGrids[station][channel][0][i];
    					}
    				}
    			}
    		}
    	}
    	return sensorGrids;
    }

    public double [] getCorrectionVector(
			int chnNum,
			double px,
			double py){
		int sensorCorrWidth= getSensorCorrWidth(chnNum);
		int indexXY=((int) Math.floor(px/getDecimateMasks(chnNum))) +
		((int) Math.floor(py/getDecimateMasks(chnNum)))*sensorCorrWidth;
		double []vector=new double[this.pixelCorrection[chnNum].length];
		for (int i=0;i<vector.length;i++) vector[i]=this.pixelCorrection[chnNum][i][indexXY];
        return vector;
    }

    /**
     * Calculate color flat-field data for the pattern grid, calculate pattern grid mask (alpha)
     * @param referenceStation - station number for unity target brightness
     * @param flatFields partial, per-station, per-sensor pattern flat-field data
     * @param shrinkForMatching shrink pattern mask for calculating pattern average (removing unreliable borders)
     * @param resetMask reset pattern mask to default before (re)-calculating mask
     * @param maxDiffNeighb maximal relative difference between neghbor nodes (ignoring off-grid)
     * @param shrinkMask shrink result mask
     * @param fadeMask smooth fade the alpha on the pattern edge, keep zeros zeros
     * @return {alpha, r,g,b,number of images used} for each view group separately
     */

    public double [][][][] combineGridFlatField(
//    		RefineParameters refineParameters,
    		int referenceStation,
    		double [][][][] flatFields,
    		double shrinkForMatching,
    		boolean resetMask,
    		double maxDiffNeighb,  // maximal allowed relative difference between neighbor nodes (relative), 0 - do not  filter any
    		int shrinkMask, // shrink result mask
    		double fadeMask
    		){
    	// So far all parameters common for large/small sensors!
//    	boolean resetMask = refineParameters.flatFieldResetMask;                // common for all sensors!
//    	int referenceStation = refineParameters.flatFieldReferenceStation;      // common for all sensors!
//    	double shrinkForMatching = refineParameters.flatFieldShrinkForMatching; // common for all sensors!
//    	double maxDiffNeighb  = refineParameters.flatFieldMaxRelDiff;           // common for all sensors!
//		int shrinkMask  = refineParameters.flatFieldShrinkMask;                 // common for all sensors!
//		double fadeMask  = refineParameters.flatFieldFadeBorder;                // common for all sensors!

    	int maskIndex=3;
    	if (resetMask) patternParameters.calculateGridGeometryAndPhotometric(false);
    	double [][][] gridGeometry= patternParameters.getGeometry();
    	int [] viewMap=	patternParameters.getViewMap();

    	int gridHeight=gridGeometry.length;
    	int gridWidth=gridGeometry[0].length;
    	int numStations=patternParameters.getNumStations();
    	int numViews=patternParameters.getNumViews();
    	double [][][][] viewPatterns=new double [numStations][numViews][][];
		double [][][] gridMask= new double[numStations][numViews][gridWidth*gridHeight];
		double [][] scaleIndividual=new double[flatFields[referenceStation].length][3]; // scale individual sensor patters before averaging
    	for (int station=0;station<numStations;station++){
    		for (int numView=0;numView<numViews;numView++){
    			viewPatterns[station][numView]=null;
//    			double [] gridMask= new double[gridWidth*gridHeight];
    			for (int v=0;v<gridHeight;v++) for (int u=0;u<gridWidth;u++) gridMask[station][numView][u+v*gridWidth]=(gridGeometry[v][u]!=null)?gridGeometry[v][u][maskIndex]:0.0;
    			if (shrinkForMatching>0){
    				(new DoubleGaussianBlur() ).blurDouble(gridMask[station][numView], gridWidth, gridHeight, shrinkForMatching, shrinkForMatching, 0.01);
    				for (int i=0;i<gridMask[station][numView].length;i++){
    					double d=2.0*(gridMask[station][numView][i]-0.5);
    					gridMask[station][numView][i]=(d>0)?(d*d):(0.0);
    				}
    			}
    			for (int v=0;v<gridHeight;v++) for (int u=0;u<gridWidth;u++){
    				if ((gridGeometry[v][u]==null) || (gridGeometry[v][u][maskIndex]<=0.0)) gridMask[station][numView][u+v*gridWidth]=0.0;
    			}
    			if (this.debugLevel>2){
    				ShowDoubleFloatArrays.showArrays(gridMask[station][numView], gridWidth, gridHeight,   "MATCH_MASK"+numView);
    			}
//    			double [][] scaleIndividual=new double[flatFields[station].length][3]; // scale individual sensor patters before averaging
    			//    		for (int numSensor=0;numSensor<flatFields.length; numSensor++ ) if (flatFields[numSensor]!=null){
    			// process only sensors from the same view of the target (i.e. 0 - eyesis head, 1 - eyesis bottom)
    			if (station==referenceStation) {
    				int numUsedSensors=0;
    				for (int numSensor=0;numSensor<flatFields[station].length; numSensor++ ) if ((flatFields[station][numSensor]!=null) && (viewMap[numSensor]==numView)){
    					numUsedSensors++;
    					double [] weightedSums={0.0,0.0,0.0};
    					double sumWeights=0;
    					for (int i=0;i<flatFields[station][numSensor][0].length;i++){
    						if ((gridMask[station][numView][i]>0.0) && (flatFields[station][numSensor][0][i]>1.0)){ // more than one overlapping image
    							double weight=flatFields[station][numSensor][0][i]*gridMask[station][numView][i];
    							sumWeights+=weight;
    							for (int c=0;c<weightedSums.length;c++) weightedSums[c]+=weight*flatFields[station][numSensor][c+1][i];

    						}
    					}
    					for (int c=0;c<weightedSums.length;c++){
    						scaleIndividual[numSensor][c]=patternParameters.averageRGB[c]*sumWeights/weightedSums[c];
    						if (this.debugLevel>2){
    							System.out.println("combineGridFlatField(): scaleIndividual["+numSensor+"]["+c+"]="+scaleIndividual[numSensor][c]);
    						}
    					}

    				}
    				if (numUsedSensors==0){
    					System.out.println("No data for target view #"+numView+" reference station ="+referenceStation);
    					continue;
    				}
    			}

    		}
    	}
    	for (int station=0;station<numStations;station++){
    		for (int numView=0;numView<numViews;numView++){
    			//    		double [][] combinedPattern=new double [5][gridWidth*gridHeight];
    			viewPatterns[station][numView]=new double [5][gridWidth*gridHeight];
    			double [][] combinedPattern=viewPatterns[station][numView];
    			for (int i=0;i<combinedPattern[0].length;i++){
    				double sumWeights=0;
    				double [] weightedSums={0.0,0.0,0.0};
    				for (int numSensor=0;numSensor<flatFields[station].length; numSensor++ ) if ((flatFields[station][numSensor]!=null) && (viewMap[numSensor]==numView)){
    					double weight=flatFields[station][numSensor][0][i];
    					sumWeights+=weight;
    					for (int c=0;c<weightedSums.length;c++) weightedSums[c]+=weight*flatFields[station][numSensor][c+1][i]*scaleIndividual[numSensor][c];
    				}
    				combinedPattern[4][i]=sumWeights; // just for debugging - no, actually used? - number of images used for this grid node
    				for (int c=0;c<weightedSums.length;c++){
    					combinedPattern[c+1][i]=(sumWeights>0.0)?(weightedSums[c]/sumWeights):0.0;
    				}

    			}
    			/*
    		}
    	}
    	for (int station=0;station<viewPatterns.length;station++){
    		for (int numView=0;numView<viewPatterns[station].length;numView++){

    			double [][] combinedPattern=viewPatterns[station][numView];
    	*/
    			//    	double [] gridMask[station][numView]= new double[gridWidth*gridHeight];
    			// calculate final mask
    			for (int v=0;v<gridHeight;v++) for (int u=0;u<gridWidth;u++) gridMask[station][numView][u+v*gridWidth]=(gridGeometry[v][u]!=null)?gridGeometry[v][u][maskIndex]:0.0;
    			if (maxDiffNeighb>0.0) { // throw away bad (having sharp gradients) nodes
    				int expWidth=gridWidth+2;
    				int expHeight=gridHeight+2;
    				double [] expandedGrid=new double [expWidth*expHeight];
    				boolean [] enabled=new boolean[expandedGrid.length];
    				for (int v=0;v<expHeight;v++) for (int u=0;u<expWidth;u++){
    					int index=u+expWidth*v;
    					if ((u==0) || (v==0) || (u==(expWidth-1)) || (v==(expHeight-1))){
    						expandedGrid[index]=0.0;
    						enabled[index]=false;
    					} else {
    						int indexSrc=(u-1)+gridWidth*(v-1);
    						expandedGrid[index]=(combinedPattern[1][indexSrc]+combinedPattern[2][indexSrc]+combinedPattern[3][indexSrc])/3.0; // average value;
    						enabled[index]=gridMask[station][numView][indexSrc]>0.0;
    					}
    				}
    				boolean [] badNodes=enabled.clone();
    				int [] dirs={
    						-expWidth-1,-expWidth,-expWidth+1,  1,
    						expWidth+1, expWidth, expWidth-1, -1};
    				int numBadOnTheBorder=1; // just to make while(true) happy
    				int minNeighb=3; // remove nodes with less than 3 neighbors
    				while (numBadOnTheBorder>0){
    					// build/update badNodes array
    					numBadOnTheBorder=0;
    					int numBad=0;
    					double [] diffs = new double [8];
    					int [] indices=new int [8];
    					for (int i=0;i<8;i++) {
    						diffs[i]=  -1.0; // diff==0 on isolated pair?
    						indices[i]=-1;
    					}
    					for (int index=0;index<badNodes.length;index++) if (badNodes[index]){
    						int numNeighb=0;
    						double maxDiff=0.0;
    						for (int dir=0;dir<dirs.length;dir++) {
    							int index1=index+dirs[dir];
    							if (enabled[index1]) {
    								numNeighb++;
    								double d=2.0*Math.abs((expandedGrid[index1]-expandedGrid[index])/(expandedGrid[index]+expandedGrid[index1]));
    								if (maxDiff<d) maxDiff=d;
    							}
    						}
    						if ((maxDiff<((maxDiffNeighb*numNeighb)/minNeighb)) && (numNeighb>=minNeighb)){ //more neighbors - more likely to keep
    							badNodes[index]=false; // rehabilitate node
    						} else {
    							numBad++;
    							if (numNeighb<8) { // do nothing if bad node is inside - it may be removed in the next passes
    								numBadOnTheBorder++;
    								if (maxDiff>diffs[numNeighb]){
    									diffs[numNeighb]=maxDiff;
    									indices[numNeighb]=index;
    								}
    							}
    						}
    					}
    					if (this.debugLevel>1) System.out.println("combineGridFlatField(): "+numBad+" bad nodes, "+numBadOnTheBorder+" of them on the border");
    					if (numBadOnTheBorder==0) break; // nothing to remove - break from the while(true) loop
    					// find bad node with least enabled neighbors - there will be one at least
    					for (int n=0;n<8;n++){
    						if (indices[n]>=0){
    							enabled[indices[n]]=false;   // disable this node
    							badNodes[indices[n]]=false;  // and remove from bad nodes (it is dead now)
    							// Any orphans around (were not bad, but now have to few neighbors)
    							for (int dir=0;dir<dirs.length;dir++) {
    								int index1=indices[n]+dirs[dir];
    								if (enabled[index1]) {
    									badNodes[index1]=true;
    								}
    							}
    							break;
    						}
    					}
    				}
    				// shrink enabled cells by shrinkMask
    				for (int n=0;n<shrinkMask;n++){
    					for (int i=0;i<badNodes.length;i++) badNodes[i]=false;
    					for (int v=1;v<(expHeight-1);v++) for (int u=1;u<(expWidth-1);u++){
    						int index=u+expWidth*v;
    						badNodes[index]=!enabled[index+1] || !enabled[index-1] || !enabled[index+expWidth] || !enabled[index-expWidth];
    					}
    					for (int i=0;i<badNodes.length;i++) if (badNodes[i]) enabled[i]=false;
    				}
    				// copy back to the gridMask[station][numView]
    				for (int v=1;v<(expHeight-1);v++) for (int u=1;u<(expWidth-1);u++){
    					int index=u+expWidth*v;
    					int indexSrc=(u-1)+gridWidth*(v-1);
    					if (!enabled[index]) gridMask[station][numView][indexSrc]=0.0;
    				}
    				for (int i=0;i<gridMask[station][numView].length;i++) if (gridMask[station][numView][i]==0.0){
    					combinedPattern[1][i]=0.0;
    					combinedPattern[2][i]=0.0;
    					combinedPattern[3][i]=0.0;
    				}

    			}
    			// fade mask on the borders, keep zeros - zeros
    			if (fadeMask>0.0){
    				double [] gridMask1=gridMask[station][numView].clone();
    				(new DoubleGaussianBlur() ).blurDouble(gridMask[station][numView], gridWidth, gridHeight, fadeMask, fadeMask, 0.01);
    				for (int i=0;i<gridMask[station][numView].length;i++){
    					double d=2.0*(gridMask[station][numView][i]-0.5);
    					gridMask[station][numView][i]=(gridMask1[i]>0)?((d>0)?(d*d):(0.0)):0.0;
    					if (combinedPattern[4][i]==0.0) gridMask[station][numView][i]=0.0; // how can it be zero combinedPattern[4][i] with non-zero gridMask[i]?
    				}
    			}
        		combinedPattern[0]=gridMask[station][numView];
    		}
    	}
    	return viewPatterns;
    }
    /**
     * Applies calculated [][] pattern alpha, r,g,b to the current grid geometry
     * @param patternFlatField
     */
    public void applyGridFlatField(
    		double [][][][] patternFlatField // {alpha, red,green,blue, number of images used}[pixel_index] for each view of the pattern
    ){
    	for (int station=0;station<patternParameters.getNumStations();station++){
    		for (int nView=0;nView<patternParameters.getNumViews();nView++)
    			if ((patternFlatField[station]!=null) && (patternFlatField[station][nView]!=null)) {
    				double [][] photometrics=patternParameters.getPhotometricByView(station,nView);
    				photometrics[0]=patternFlatField[station][nView][1].clone(); // red
    				photometrics[1]=patternFlatField[station][nView][2].clone(); // green
    				photometrics[2]=patternFlatField[station][nView][3].clone(); // blue
    				photometrics[3]=patternFlatField[station][nView][0].clone(); // alpha
    			}
    	}
    	/*
    	double [][][] gridGeometry= patternParameters.getGeometry();
    	int gridHeight=gridGeometry.length;
    	int gridWidth=gridGeometry[0].length;
    	for (int v=0;v<gridHeight;v++) for (int u=0;u<gridWidth;u++) {
    		int index=u+v*gridWidth;
    		gridGeometry[v][u][3]=patternFlatField[0][index];
    		gridGeometry[v][u][4]=patternFlatField[1][index];
    		gridGeometry[v][u][5]=patternFlatField[2][index];
    		gridGeometry[v][u][6]=patternFlatField[3][index];
    	}
    	 */
    }
    /**
     * Remove areas on the target flat-field data with specular reflections of the lamps by matching different views
     * @param highPassSigma - subtract weighted average of the difference with this
     * @param thershold mismatch causing 50% drop of the weight function
     * @param numIterations number of iterations of comparing to the weighted/masked average
     * @param apply apply changes
     * @param debugShowMode 0 - do not show debug images, 1 show only during last iteration, 2 - show always
     */
    public void removeSpecular (
    		boolean positiveDiffOnly,
    		double highPassSigma,
    		double lowPassSigma,
    		double thershold,
    		int numIterations,
    		boolean apply,
    		int debugShowMode){ // 0 - none, 1 - last iteration, 2 - all iterations
    	int debugThreshold=1;
    	int length=0;

    	double [][][] weights=new double [patternParameters.getNumStations()][patternParameters.getNumViews()][];
		double [][][][] photometrics=new double [patternParameters.getNumStations()][patternParameters.getNumViews()][][];
		double [][][][] highPassDiff=new double [patternParameters.getNumStations()][patternParameters.getNumViews()][][];
		double [][][][] lowPassDiff=new double [patternParameters.getNumStations()][patternParameters.getNumViews()][][];

		int width=  patternParameters.gridGeometry[0].length;
		int height= patternParameters.gridGeometry.length;
    	for (int station=0;station<patternParameters.getNumStations();station++){
    		for (int nView=0;nView<patternParameters.getNumViews();nView++) {
    			photometrics[station][nView]=patternParameters.getPhotometricByView(station,nView);
    			if (photometrics[station][nView]!=null){
    				length=photometrics[0][0][3].length; // should all be the same length (or null)
    				weights[station][nView]=new double [length];
    				for (int nPix=0;nPix<length;nPix++) weights[station][nView][nPix]=(photometrics[station][nView][3][nPix]>0.0)?1.0:0.0;
    			} else {
    				weights[station][nView]=null;
    			}
    		}
    	}
    	double threshold23=9.0*thershold*thershold;
    	for (int nIter=0;nIter<numIterations;nIter++){
    		boolean showDebug=(debugShowMode>1) || ((debugShowMode>0) && (nIter== (numIterations-1)));
        	// Calculate weighted average among different stations/views.
    		double [][] average=new double [4][length];
    		for (int nPix=0;nPix<length;nPix++){
    			double w0=0.0;
    	    	for (int station=0;station<patternParameters.getNumStations();station++){
    	    		for (int nView=0;nView<patternParameters.getNumViews();nView++) {
    	    			if (photometrics[station][nView]!=null){
    	    				double w=weights[station][nView][nPix]*photometrics[station][nView][3][nPix];
    	    				average[0][nPix]+=w*photometrics[station][nView][0][nPix];
    	    				average[1][nPix]+=w*photometrics[station][nView][1][nPix];
    	    				average[2][nPix]+=w*photometrics[station][nView][2][nPix];
    	    				w0+=w;
    	    			}
    	    		}
    	    	}
    	    	double k= (w0>0.0)?(1.0/w0):0.0;
    	    	average[0][nPix]*=k;
    	    	average[1][nPix]*=k;
    	    	average[2][nPix]*=k;
    	    	average[3][nPix]=w0/(patternParameters.getNumStations()*patternParameters.getNumViews());
    		}
    		double [][][][] diffFromAverage=new double [photometrics.length][photometrics[0].length][4][length];
        	// Scale each station/view for best fit
	    	for (int station=0;station<patternParameters.getNumStations();station++){
	    		for (int nView=0;nView<patternParameters.getNumViews();nView++) {
	    			double scale=0.0;
	    			if (photometrics[station][nView]!=null){
	    				double [] weightsHighLowPass=new double[length];
	    				double sf2=0.0,sfg=0.0;
	    				for (int nPix=0;nPix<length;nPix++){
    	    				double w=weights[station][nView][nPix]*photometrics[station][nView][3][nPix];
    	    				weightsHighLowPass[nPix]=w;
    	    				sf2+=w*(photometrics[station][nView][0][nPix]*photometrics[station][nView][0][nPix]+
    	    						photometrics[station][nView][1][nPix]*photometrics[station][nView][1][nPix]+
    	    						photometrics[station][nView][2][nPix]*photometrics[station][nView][2][nPix]);
    	    				sfg+=w*(photometrics[station][nView][0][nPix]*average[0][nPix]+
    	    						photometrics[station][nView][1][nPix]*average[1][nPix]+
    	    						photometrics[station][nView][2][nPix]*average[2][nPix]);

	    				}
	    				scale=sfg/sf2;
	    				if ((this.debugLevel>=debugThreshold) && showDebug){
	    					System.out.println("removeSpecular(), pass"+nIter+" scale["+station+"]["+nView+"]="+scale);
	    				}
//	    				Calculate difference from average
	    				for (int nPix=0;nPix<length;nPix++){
	    					if (photometrics[station][nView][3][nPix]>0.0){
	    						for (int c=0;c<3;c++){
	    							double d=scale*photometrics[station][nView][c][nPix]-average[c][nPix];
	    							diffFromAverage[station][nView][c][nPix]=d;
	    						}
	    					}
	    				}
	    				if (highPassSigma>0.0){
	    					double [] weightsHighPass=weightsHighLowPass.clone();
	    					(new DoubleGaussianBlur()).blurDouble(
	    							weightsHighPass,
	    							width,
	    							height,
	    							highPassSigma,
	    							highPassSigma,
	    							0.01);
	    					highPassDiff[station][nView]=new double [3][];
	    					for (int c=0;c<3;c++){
	    	    				highPassDiff[station][nView][c]=diffFromAverage[station][nView][c].clone(); // deep

	    						for (int nPix=0;nPix<length;nPix++){
	    							highPassDiff[station][nView][c][nPix]*=weightsHighLowPass[nPix];
	    						}
		    					(new DoubleGaussianBlur()).blurDouble(
		    							highPassDiff[station][nView][c],
		    							width,
		    							height,
		    							highPassSigma,
		    							highPassSigma,
		    							0.01);
	    						for (int nPix=0;nPix<length;nPix++){
	    							highPassDiff[station][nView][c][nPix]=
	    								diffFromAverage[station][nView][c][nPix]-
	    								((weightsHighPass[nPix]>0)?(highPassDiff[station][nView][c][nPix]/weightsHighPass[nPix]):0.0);
	    						}
	    					}
	    				} else {
    	    				highPassDiff[station][nView]=diffFromAverage[station][nView].clone(); // shallow
	    				}

	    				if (lowPassSigma>0.0){
	    					double [] weightsLowPass=weightsHighLowPass.clone();
	    					(new DoubleGaussianBlur()).blurDouble(
	    							weightsLowPass,
	    							width,
	    							height,
	    							lowPassSigma,
	    							lowPassSigma,
	    							0.01);
	    					lowPassDiff[station][nView]=new double [3][];
	    					for (int c=0;c<3;c++){
	    						lowPassDiff[station][nView][c]=highPassDiff[station][nView][c].clone();
	    						for (int nPix=0;nPix<length;nPix++){
	    							lowPassDiff[station][nView][c][nPix]*=weightsHighLowPass[nPix];
	    						}
		    					(new DoubleGaussianBlur()).blurDouble(
		    							lowPassDiff[station][nView][c],
		    							width,
		    							height,
		    							lowPassSigma,
		    							lowPassSigma,
		    							0.01);
	    						for (int nPix=0;nPix<length;nPix++){
	    							lowPassDiff[station][nView][c][nPix]=
	    								(weightsLowPass[nPix]>0)?(lowPassDiff[station][nView][c][nPix]/weightsLowPass[nPix]):0.0;
	    						}
	    					}
	    				} else {
    	    				lowPassDiff[station][nView]=highPassDiff[station][nView].clone(); // shallow
    					}



// TODO: display, calculate new weight from filtered difference.
// Calculate new weight
	    				for (int nPix=0;nPix<length;nPix++){
	    					if (photometrics[station][nView][3][nPix]>0.0){
	    						double e2=0.0;
	    						for (int c=0;c<3;c++){
	    							double d=lowPassDiff[station][nView][c][nPix];
//	    							double d=scale*photometrics[station][nView][c][nPix]-average[c][nPix];
//	    							diffFromAverage[station][nView][c][nPix]=d;
	    							if (!positiveDiffOnly || (d>0)) e2+=d*d;
	    						}
	    						weights[station][nView][nPix]=1.0/(e2/threshold23+1.0);
	    					} else {
	    						weights[station][nView][nPix]=0.0;
	    					}
	    				}

	    			}
	    		}
	    	}
	    	if ((this.debugLevel>=debugThreshold) && showDebug) {
	    		String [] titles = new String [weights.length*weights[0].length];
	    		double [][] debugData= new double [weights.length*weights[0].length][];
	        	for (int station=0;station<patternParameters.getNumStations();station++){
	        		for (int nView=0;nView<patternParameters.getNumViews();nView++) {
	        			int n=station*weights[0].length+nView;
	        			titles[n]="S"+station+" V"+nView;
	        			if (photometrics[station][nView]!=null){
	        				debugData[n]=weights[station][nView];
	        			}
	        		}
	        	}
	        	ShowDoubleFloatArrays.showArrays(
	        			debugData,
	        			width,
	        			height,
	        			true,
	        			"GridWeights"+nIter,
	        			titles);
	        	double [][] debugDiffGreen= new double [weights.length*weights[0].length][];
	        	double [][] debugHighpassDiffGreen= new double [weights.length*weights[0].length][];
	        	double [][] debugLowpassDiffGreen= new double [weights.length*weights[0].length][];
	        	for (int station=0;station<patternParameters.getNumStations();station++){
	        		for (int nView=0;nView<patternParameters.getNumViews();nView++) {
	        			int n=station*weights[0].length+nView;
	        			debugDiffGreen[n]=diffFromAverage[station][nView][1];
	        			debugHighpassDiffGreen[n]=highPassDiff[station][nView][1];
	        			debugLowpassDiffGreen[n]=lowPassDiff[station][nView][1];
	        		}
	        	}
	        	if (this.debugLevel>=(debugThreshold+1)) ShowDoubleFloatArrays.showArrays(
	        			debugDiffGreen,
	        			width,
	        			height,
	        			true,
	        			"DiffGreen"+nIter,
	        			titles);
	        	if (this.debugLevel>=(debugThreshold+1)) ShowDoubleFloatArrays.showArrays(
	        			debugHighpassDiffGreen,
	        			width,
	        			height,
	        			true,
	        			"HighpassGreen"+nIter,
	        			titles);

	        	ShowDoubleFloatArrays.showArrays(
	        			debugLowpassDiffGreen,
	        			width,
	        			height,
	        			true,
	        			"LowpassGreen"+nIter,
	        			titles);

	        	String [] averageTitles={"Red","Green","Blue","Weight"};
	        	ShowDoubleFloatArrays.showArrays(
	        			average,
	        			width,
	        			height,
	        			true,
	        			"Average-"+nIter,
	        			averageTitles);

	    	}
    	} // for (int nIter=0;nIter<numIterations;nIter++){
    	// Apply new weights
    	if (apply) {
    		for (int station=0;station<patternParameters.getNumStations();station++){
    			for (int nView=0;nView<patternParameters.getNumViews();nView++) {
    				if (photometrics[station][nView]!=null){
    					for (int nPix=0;nPix<length;nPix++) photometrics[station][nView][3][nPix]*=weights[station][nView][nPix];
    				}
    			}
    		}
    	}

    }

    /**
     *
     * @param shrink sensor mask by this amount (sensor, non-decimated pixels)
     * @param radius radial mask - zero if farther than radius, 0.5*(cos(pi*r/radius)+1.0) if less
     * @param minimalAlpha - zero mask below this threshold
     * @return returns arrray with the same size as sensorMask that corresponds to low-vignetting areas of the sensor/lens
     */
    // Using station 0 - should be not much difference
    // older version with one type of sensors
    public double [][] nonVignettedMasks(
    		double shrink,
    		double radius,
    		double minimalAlpha){
    	if (this.pixelCorrection==null){
    		initSensorCorrection();
    	}
    	double [][]masks=new double [this.pixelCorrection.length][];
    	int maskIndex=2;
    	for (int numSensor=0;numSensor<masks.length;numSensor++){
//    		boolean small_sensor = fittingStrategy.distortionCalibrationData.isSmallSensor(numSensor);
			boolean small_sensor = fittingStrategy.distortionCalibrationData.getSmallSensors()[numSensor];
    		RefineParameters rp = small_sensor ? refineParameters.refineParametersSmall : refineParameters;

    		if (this.pixelCorrection[numSensor]==null) masks[numSensor] = null;
    		else {
    			masks[numSensor] = fittingStrategy.distortionCalibrationData.nonVignettedMask(
    					this.pixelCorrection[numSensor][maskIndex],
    					getSensorWidth(numSensor), // this.pixelCorrectionWidth,
    					getSensorHeight(numSensor), // this.pixelCorrectionHeight,
    					fittingStrategy.distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[0][numSensor].px0,     // lens center X (sensor, non-decimated pix)
    					fittingStrategy.distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[0][numSensor].py0,     // lens center Y (sensor, non-decimated pix)
    					shrink,
    					radius,
    					minimalAlpha);
    			//    			System.out.println("nonVignettedMasks(), masks["+numSensor+"].length="+masks[numSensor].length);
    		}
    	}
    	return masks;
    }

// modified for multiple types of sensors:
    public double [][] nonVignettedMasks(
    		RefineParameters refineParameters){
    	if (this.pixelCorrection==null){
    		initSensorCorrection();
    	}
    	double [][]masks=new double [this.pixelCorrection.length][];
    	int maskIndex=2;
    	for (int numSensor=0;numSensor<masks.length;numSensor++){
//				boolean small_sensor = fittingStrategy.distortionCalibrationData.isSmallSensor(numSensor);
				boolean small_sensor = fittingStrategy.distortionCalibrationData.getSmallSensors()[numSensor];
				RefineParameters rp = small_sensor ? refineParameters.refineParametersSmall : refineParameters;

    		double shrink = rp.flatFieldShrink;
    		double radius = rp.flatFieldNonVignettedRadius;
    		double minimalAlpha=rp.flatFieldMinimalAlpha;

    		if (this.pixelCorrection[numSensor]==null) masks[numSensor] = null;
    		else {
    			masks[numSensor] = fittingStrategy.distortionCalibrationData.nonVignettedMask(
    					this.pixelCorrection[numSensor][maskIndex],
    					getSensorWidth(numSensor), // this.pixelCorrectionWidth,
    					getSensorHeight(numSensor), // this.pixelCorrectionHeight,
    	        		fittingStrategy.distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[0][numSensor].px0,     // lens center X (sensor, non-decimated pix)
    	        		fittingStrategy.distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[0][numSensor].py0,     // lens center Y (sensor, non-decimated pix)
    	        		shrink,
    	        		radius,
    	        		minimalAlpha);
//    			System.out.println("nonVignettedMasks(), masks["+numSensor+"].length="+masks[numSensor].length);
    		}
    	}
    	return masks;
    }




    // Use series0 to find grid mismatch (and set it correctly). Uses that pattern in the world coordinate system is
    // approximately in XY plane, so by freezing all other parameters but GXY0 and GXY1 it is possible to find
    // the pattern grid match.
    public int [] findImageGridOffset(
    		int     num_img,
    		int     ser_num, // number of series to reprogram (< 0 - from the last one)
    		boolean adjust_attitude, // true for eo, false for lwir (uses exact attitude from eo)
    		boolean even,
    		PatternParameters patternParameters,
    		double [] stats) { // result quality, null OK
    	if (ser_num < 0) {
    		ser_num = fittingStrategy.parameterMode.length + ser_num; // use from the last one
    	}
		if (this.debugLevel > 0) {
			System.out.println("Will use/modify fitting series "+ser_num+" to adjust image "+num_img);
			if (adjust_attitude) {
				System.out.println("Will adjust GXYZ0, GXYZ1, goniometerHorizontal, and goniometerAxial");
			} else {
				System.out.println("Will adjust only GXYZ0 and GXYZ1.");
			}
		}

    	int was_seriesNumber = seriesNumber;
		// set series ser_num (was 0) to this set images
		DistortionCalibrationData dcd = fittingStrategy.distortionCalibrationData;
		int par_index_goniometerHorizontal = dcd.getParameterIndexByName("goniometerHorizontal");
		int par_index_goniometerAxial =      dcd.getParameterIndexByName("goniometerAxial");
		int par_index_GXYZ0 =                dcd.getParameterIndexByName("GXYZ0");
		int par_index_GXYZ1 =                dcd.getParameterIndexByName("GXYZ1");
		for (int i = 0; i < fittingStrategy.parameterMode[ser_num].length; i++) {
			 fittingStrategy.parameterMode[ser_num][i] = FittingStrategy.modeFixed; // 0
		}
		fittingStrategy.parameterMode[ser_num][par_index_GXYZ0] =                    FittingStrategy.modeIndividual; // 3
		fittingStrategy.parameterMode[ser_num][par_index_GXYZ1] =                    FittingStrategy.modeIndividual; // 3
		if (adjust_attitude) {
			fittingStrategy.parameterMode[ser_num][par_index_goniometerHorizontal] = FittingStrategy.modeIndividual; // 3
			fittingStrategy.parameterMode[ser_num][par_index_goniometerAxial] =      FittingStrategy.modeIndividual; // 3
		}

		boolean [] selection = fittingStrategy.selectAllImages(ser_num); // enable all images in series 0
		for (int i=0;i<selection.length;i++) selection[i]=false;
		selection[num_img]=true;
		fittingStrategy.setImageSelection(ser_num,selection);
		seriesNumber=   ser_num; // start from 0;
		initFittingSeries(false, filterForAll,ser_num); // will set this.currentVector, will build parameter map too
		//this.stopAfterThis[numSeries]
		fittingStrategy.stopAfterThis[ser_num]=true;
		stopEachStep=      false;
		stopEachSeries=    false; // will not ask for confirmation after done
		lambda =           fittingStrategy.lambdas[ser_num];
		boolean   LMA_OK = false;
		try {
			LMA_OK = 	LevenbergMarquardt (false, true); //  skip dialog, dry run (no updates)
		} catch (Exception e) {
			// LMA failed - e.g. not enough points (Singular Matrix)
		}
		if (!LMA_OK) {
			return null; // LMA did not converge
		}

//		DistortionCalibrationData dcd = this.fittingStrategy.distortionCalibrationData;
// find indices of GXYZ0 and GXYZ1 in the vector
		int index_GXYZ0 = this.fittingStrategy.reverseParameterMap[num_img][par_index_GXYZ0];
		int index_GXYZ1 = this.fittingStrategy.reverseParameterMap[num_img][par_index_GXYZ1];

		double [] new_XY = {this.currentVector[index_GXYZ0],this.currentVector[index_GXYZ1]};
//		DistortionCalibrationData dcd = this.fittingStrategy.distortionCalibrationData;
//		int num_set = dcd.gIP[num_img].getSetNumber();
		double [] 	ref_XYZ = dcd.getXYZ(num_img);
		double []   diff_XY = {
				new_XY[0] -ref_XYZ[0],
				new_XY[1] -ref_XYZ[1]};

//save safe settings to run LMA manually (or save what was set)
		seriesNumber=      was_seriesNumber; // start from 0;
		initFittingSeries  (false,filterForAll,seriesNumber); // will set this.currentVector
		stopEachSeries=    true; // will not ask for confirmation after done
		stopOnFailure=     true;
		lambda=            fittingStrategy.lambdas[seriesNumber];
		double [] errs = new double[2];

		int [] grid_offset = dcd.suggestOffset (
        		num_img,
        		diff_XY, // This XYZ minus reference XYZ  z is not used, may be just[2]
        		even,
        		patternParameters,
        		errs);
		if (stats != null) {
			stats[0] = this.currentRMS;
			stats[1] = errs[0]; // dU
			stats[2] = errs[1]; // dV
		}
		return grid_offset;
    }

    public boolean adjustAttitudeAfterOffset(
    		int     num_img,
    		int     ser_num, // number of series to reprogram
    		PatternParameters patternParameters) {
    	int [] num_imgs = {num_img};
    	return  adjustAttitudeAfterOffset(
        		num_imgs,
        		ser_num, // number of series to reprogram
        		patternParameters);
    }

    public boolean adjustAttitudeAfterOffset(
    		int  []   num_imgs,
    		int     ser_num, // number of series to reprogram
    		PatternParameters patternParameters) {
    	if (ser_num < 0) {
    		ser_num = fittingStrategy.parameterMode.length + ser_num; // use from the last one
    	}
		if (this.debugLevel > 0) {
			System.out.print("Will use/modify fitting series "+ser_num+" for to adjust  az, tilt of images: ");
			for (int num_img: num_imgs) {
				System.out.print(" "+num_img);
			}
			System.out.println();
			System.out.println("Will adjust goniometerHorizontal, and goniometerAxial");
		}
    	int was_seriesNumber = seriesNumber;
		// set series ser_num (was 0) to this set images
		DistortionCalibrationData dcd = fittingStrategy.distortionCalibrationData;
		int par_index_goniometerHorizontal = dcd.getParameterIndexByName("goniometerHorizontal");
		int par_index_goniometerAxial =      dcd.getParameterIndexByName("goniometerAxial");
		for (int i = 0; i < fittingStrategy.parameterMode[ser_num].length; i++) {
			 fittingStrategy.parameterMode[ser_num][i] = FittingStrategy.modeFixed; // 0
		}
		fittingStrategy.parameterMode[ser_num][par_index_goniometerHorizontal] = FittingStrategy.modeIndividual; // 3
		fittingStrategy.parameterMode[ser_num][par_index_goniometerAxial] =      FittingStrategy.modeIndividual; // 3
		boolean [] selection = fittingStrategy.selectAllImages(ser_num); // enable all images in series 0
		for (int i=0;i<selection.length;i++) selection[i]=false;
		for (int num_img:num_imgs) {
			selection[num_img]=true;
		}
		fittingStrategy.setImageSelection(ser_num,selection);
		seriesNumber=   ser_num; // start from 0;
		initFittingSeries(false, filterForAll,ser_num); // will set this.currentVector, will build parameter map too
		fittingStrategy.stopAfterThis[ser_num]=true;
		stopEachStep=      false;
		stopEachSeries=    false; // will not ask for confirmation after done
		lambda =           fittingStrategy.lambdas[ser_num];
		boolean   LMA_OK = false;
		try {
			LMA_OK = 	LevenbergMarquardt (false, false); // true); //  skip dialog, dry run (no updates)
		} catch (Exception e) {
			// LMA failed - e.g. not enough points (Singular Matrix)
		}
		if (!LMA_OK) {
			return false; // LMA did not converge
		}
//save safe settings to run LMA manually (or save what was set)
		seriesNumber=      was_seriesNumber; // start from 0;
		initFittingSeries  (false,filterForAll,seriesNumber); // will set this.currentVector
		stopEachSeries=    true; // will not ask for confirmation after done
		stopOnFailure=     true;
		lambda=            fittingStrategy.lambdas[seriesNumber];
		return true;
    }

    
    public boolean LevenbergMarquardt(
    		boolean openDialog,
    		boolean dry_run){ // do not save results
    	return LevenbergMarquardt(
        		openDialog,
        		dry_run,
        		false);
    }

    public boolean LevenbergMarquardt(
    		boolean openDialog,
    		boolean dry_run,  // do not save results
    		boolean calc_dUV){
    	if (this.fittingStrategy==null) {
        		String msg="Fitting strategy does not exist, exiting";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
    	}
//    	fittingStrategy.distortionCalibrationData.readAllGrids();
    	if (openDialog && !selectLMAParameters()) return false;
    	this.startTime=System.nanoTime();
//    	while (this.seriesNumber<fittingStrategy.getNumSeries()){ // TODO: Add "stop" tag to series
    	this.firstRMS=-1; //undefined
    	this.fittingStrategy.invalidateSelectedImages(this.seriesNumber); // undo any filters, only user selection of the  images will be in effect
    	while (this.fittingStrategy.isSeriesValid(this.seriesNumber)){ // TODO: Add "stop" tag to series
    		this.currentVector=null; // invalidate for the new series
    		boolean wasLastSeries=false;
    		while (true) { // loop for the same series

    			boolean [] state=stepLevenbergMarquardtFirst(this.seriesNumber, calc_dUV);
    			if (!this.fittingStrategy.isSeriesValid(this.seriesNumber)){
    				System.out.println("Series "+this.seriesNumber+" is invalid when weight function filters are applied (probably removed some images)");
    				return false;
    			}
    			if (state==null) {
    				String msg="Calculation aborted by user request";
    				IJ.showMessage(msg);
    				System.out.println(msg);
    				return false;
    			}

    			if (this.debugLevel>1) System.out.println(this.seriesNumber+":"+this.iterationStepNumber+": stepLevenbergMarquardtFirst("+this.seriesNumber+")==>"+state[1]+":"+state[0]);
				boolean cont=true;
				// Make it success if this.currentRMS<this.firstRMS even if LMA failed to converge
				if (state[1] && !state[0] && (this.firstRMS>this.currentRMS)){
					if (this.debugLevel>1) System.out.println("LMA failed to converge, but RMS improved from the initial value ("+this.currentRMS+" < "+this.firstRMS+")");
					state[0]=true;
				}
    			if (
    					(this.stopRequested.get()>0) || // graceful stop requested
    					(this.stopEachStep) ||
    					(this.stopEachSeries && state[1]) ||
    					(this.stopOnFailure && state[1] && !state[0])){

    				if (this.debugLevel>0){
    					if (this.stopRequested.get()>0) System.out.println("User requested stop");
    					System.out.println("LevenbergMarquardt(): series:step ="+this.seriesNumber+":"+this.iterationStepNumber+
    						", RMS="+IJ.d2s(this.currentRMS,8)+
    						" ("+IJ.d2s(this.firstRMS,8)+") "+
    						", RMSPure="+IJ.d2s(this.currentRMSPure,8)+
    						" ("+IJ.d2s(this.firstRMSPure,8)+
    						") at "+ IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3));
    				}
    				long startDialogTime=System.nanoTime();
    				cont=dialogLMAStep(state);
    				this.stopRequested.set(0); // Will not stop each run
    				this.startTime+=(System.nanoTime()-startDialogTime); // do not count time used by the User.
    				if (this.showThisImages) showDiff (this.currentfX, "fit-"+this.iterationStepNumber);
    				if (this.showNextImages) showDiff (this.nextfX,    "fit-"+(this.iterationStepNumber+1));
    			} else if (this.debugLevel>1){
					System.out.println("==> LevenbergMarquardt(): before action series:step ="+this.seriesNumber+":"+this.iterationStepNumber+
    						", RMS="+IJ.d2s(this.currentRMS,8)+
    						" ("+IJ.d2s(this.firstRMS,8)+") "+
    						", RMSPure="+IJ.d2s(this.currentRMSPure,8)+
    						" ("+IJ.d2s(this.firstRMSPure,8)+
    						") at "+ IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3));
    			}
				stepLevenbergMarquardtAction(); // apply step - in any case?
				if (this.updateStatus){
   	   				IJ.showStatus(this.seriesNumber+": "+"Step #"+this.iterationStepNumber+
   	   						" RMS="+IJ.d2s(this.currentRMS,8)+
   	   						" ("+IJ.d2s(this.firstRMS,8)+")"+
   	   						" RMSPure="+IJ.d2s(this.currentRMSPure,8)+
   	   						" ("+IJ.d2s(this.firstRMSPure,8)+")"+
   	   						" ");
//   	   				showStatus(this.seriesNumber+": "+"Step #"+this.iterationStepNumber+" RMS="+IJ.d2s(this.currentRMS,8)+ " ("+IJ.d2s(this.firstRMS,8)+")",0);
				}
				if ((this.debugLevel>0) && ((this.debugLevel>1) || ((System.nanoTime()-this.startTime)>10000000000.0))){ // > 10 sec
					System.out.println("--> LevenbergMarquardt(): series:step ="+this.seriesNumber+":"+this.iterationStepNumber+
							", RMS="+IJ.d2s(this.currentRMS,8)+
							" ("+IJ.d2s(this.firstRMS,8)+") "+
							", RMSPure="+IJ.d2s(this.currentRMSPure,8)+
							" ("+IJ.d2s(this.firstRMSPure,8)+
							") at "+ IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3));
				}
				if (!cont){
					if (this.saveSeries && !dry_run) {
						saveFittingSeries(); // will save series even if it ended in failure, vector will be only updated
						updateCameraParametersFromCalculated(true); // update camera parameters from all (even disabled) images
						updateCameraParametersFromCalculated(false); // update camera parameters from enabled only images (may overwrite some of the above)
					}
					// if RMS was decreased. this.saveSeries==false after dialogLMAStep(state) only if "cancel" was pressed
					return this.saveSeries; // TODO: Maybe change result?
				}
//stepLevenbergMarquardtAction();
    			if (state[1]) {
    				if (!state[0]) return false; // sequence failed
    				if (dry_run) {
    					wasLastSeries= true; // always just one series
    				} else {
    					saveFittingSeries();
    					updateCameraParametersFromCalculated(true); // update camera parameters from all (even disabled) images
    					updateCameraParametersFromCalculated(false); // update camera parameters from enabled only images (may overwrite some of the above)
    					wasLastSeries=this.fittingStrategy.isLastSeries(this.seriesNumber);
    					this.seriesNumber++;
    				}
    				break; // while (true), proceed to the next series
    			}
    		}
//    		if (this.fittingStrategy.isLastSeries(this.seriesNumber)) break;
    		if (wasLastSeries) break;
//    		this.seriesNumber++;
    	} // while (this.fittingStrategy.isSeriesValid(this.seriesNumber)){ // TODO: Add "stop" tag to series
		if (this.debugLevel>0) System.out.println("LevenbergMarquardt(): series="+this.seriesNumber+
				", RMS="+this.currentRMS+
				" ("+this.firstRMS+") "+
				", RMSPure="+this.currentRMSPure+
				" ("+this.firstRMSPure+
				") at "+ IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3));
    	return true; // all series done

    }


    /**
     * Show debug image (see showDiff (int imgNumber, double [] fX, String title ) above)
     * for each image used in the current fitting series
     * @param fX - calculated data for all images (use with this.Y)
     * @param title - Image title
     */
    	public void showDiff (double [] fX, String title ){
    		boolean [] selectedImages=fittingStrategy.selectedImages();
    		for (int imgNum=0;imgNum<selectedImages.length;imgNum++) if (selectedImages[imgNum]) {
    			showDiff (imgNum, fX, title+"-"+imgNum);
    		}
    	}

	/**
	 * Shows a 7-slice image for provided f(X) array (this.Y is also used):
	 * 1 - distance   - sqrt (dx^2+dy^2)
	 * 2 - difference for pixel-X
	 * 3 - difference for pixel-Y
	 * 4 - calculated pixel-X
	 * 5 - calculated pixel-Y
	 * 6 - measured   pixel-X
	 * 7 - measured   pixel-Y
	 * @param imgNumber - number of image
	 * @param fX - calculated data for all images (use with this.Y)
	 * @param title - Image title
	 */

	public void showDiff (int imgNumber, double [] fX, String title ){
		String [] titles={"distance","diff-X","diff-Y","f(x)-X","f(x)-Y","y-X","y-Y"};
		double [] diff=calcYminusFx(fX);
// find data range for the selected image
		int index=0;
		int numImg=fittingStrategy.distortionCalibrationData.getNumImages();
		boolean [] selectedImages=fittingStrategy.selectedImages();
		for (int imgNum=0;(imgNum<imgNumber) && (imgNum<numImg) ;imgNum++) if (selectedImages[imgNum])
			index+=fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV.length;
		int width=getGridWidth();
		if (this.debugLevel>1) {
			System.out.println("showDiff("+imgNumber+",...): fX.length="+fX.length+" this image index="+index);
		}
		double [][] imgData=new double[7][getGridHeight() * width];
		for (int i=0;i<imgData.length;i++) for (int j=0;j<imgData[i].length;j++)imgData[i][j]=0.0;

		for (int i=0;i<fittingStrategy.distortionCalibrationData.gIP[imgNumber].pixelsUV.length;i++){
			int u=fittingStrategy.distortionCalibrationData.gIP[imgNumber].pixelsUV[i][0]+patternParameters.U0;
			int v=fittingStrategy.distortionCalibrationData.gIP[imgNumber].pixelsUV[i][1]+patternParameters.V0;
			int vu=u+width*v;
			imgData[0][vu]=   Math.sqrt(diff[2*(index+i)]*diff[2*(index+i)] + diff[2*(index+i)+1]*diff[2*(index+i)+1]);
			imgData[1][vu]=   diff[2*(index+i)]; // out of bound 1410
			imgData[2][vu]=   diff[2*(index+i)+1];
			imgData[3][vu]=     fX[2*(index+i)];
			imgData[4][vu]=     fX[2*(index+i)+1];
			imgData[5][vu]= this.Y[2*(index+i)];
			imgData[6][vu]= this.Y[2*(index+i)+1];
		}
		ShowDoubleFloatArrays.showArrays(imgData, width, getGridHeight(),  true, title, titles);
	}
	/**
	 * Calculates corrections to X and Y coordinates of the grid nodes
	 * @param variationPenalty - cost of different Z for different stations
	 * @param fixXY - if true, do not adjust X,Y - only Z
	 * @param stationGroupsIn - consider some stations have the same pattern - assign them the same number. Negative - do not process the station
	 * @param grid3DCorrection - if true - calculate 3d correction, false - slow 3d (2d perpendicular to view)
	 * @param maxZCorr - maximal allowed correction in Z-direction (if wants more, will fall back to 2-d correction (perpendicular to the view)
	 * @param showIndividual - show individual images
	 * @return  combination of 3 arrays: 1 (original) - first index - 0 - correction x (mm), 1 - correction y(mm), 2 - correction z(mm)  3 - weight (number of images used)
	 * 2 - gridZCorr3d - per station differential Z correction
	 * 3 - gridZCorr3dWeight - per station weight of Z-corrections
	 */

	public double [][][] calculateGridXYZCorr3D(
			double variationPenalty,
            boolean fixXY,
            int [] stationGroupsIn,
			boolean grid3DCorrection,
			boolean rotateCorrection,
			double maxZCorr,
			boolean noFallBack,
			boolean showIndividual,
			int threadsMax,
			boolean updateStatus
			){
		int debugThreshold=2;
		// Normalize stationGroups
		int numStations=fittingStrategy.distortionCalibrationData.getNumStations();
		int [] stationGroups=new int [numStations];
		int [] stationGroupsTmp=(stationGroupsIn==null)?(new int [0]):stationGroupsIn.clone();
		for (int i=0;i<numStations;i++) stationGroups[i]=-1;
		int numberOfZGroups=0;
		for (int i=0;i<stationGroupsTmp.length;i++) if (stationGroupsTmp[i]>=0){
			for (int j=i;j<stationGroupsTmp.length;j++) if (stationGroupsTmp[j]==stationGroupsTmp[i]){
				stationGroups[j]=numberOfZGroups;
				if (j>i) stationGroupsTmp[j]=-1;
			}
			numberOfZGroups++;
		}
		if (numberOfZGroups==0) {
			System.out.println ("calculateGridXYZCorr3D(), no groups defined - using a single group for all stations");
			numberOfZGroups=1;
			for (int i=0;i<numStations;i++) stationGroups[i]=0;
		}
		if (this.debugLevel>1) {
			System.out.println ("calculateGridXYZCorr3D(), groups: "+numberOfZGroups);
			for (int i=0;i<stationGroups.length;i++) if (stationGroups[i]>=0){
				System.out.println ("  station "+i+": group "+stationGroups[i]);
			}
		}

		int width=getGridWidth();
		int height=getGridHeight();
		boolean [] selectedImages=fittingStrategy.selectedImages();
		double [][] cameraXYZ=new double [selectedImages.length][];

		double [][][] gridCorr2d=calculateGridXYZCorr2D(
				width,
				height,
				stationGroups,
				selectedImages,
				cameraXYZ,
				this.lensDistortionParameters,
				showIndividual,
				threadsMax,
				updateStatus);


		IJ.showStatus("Calculating pattern 3d correction...");
// now using gridCorr2d[imgNum], cameraXYZ[imgNum] and patternParameters.gridGeometry[v][u] find the 3d correction     public double [][][] gridGeometry=null; // [v][u]{x,y,z,"alpha"} alpha=0 - no ghrid, 1 - grid
		double [][] gridCorr3d=  new double [4][width*height];
		double [][] gridZCorr3d =new double [numStations][width*height];
		double [][] gridZCorr3dWeight =new double [numStations][width*height];
		for (int n=0;n<gridCorr3d.length;n++) for (int i=0;i<gridCorr3d[0].length;i++) gridCorr3d[n][i]=0.0;
		for (int n=0;n<gridZCorr3d.length;n++) for (int i=0;i<gridZCorr3d[0].length;i++){
			gridZCorr3d[n][i]=0.0;
			gridZCorr3dWeight[n][i]=0.0;
		}
		double Cx,Cy,Cz,Cxy,Cxz,Cyz;
		double [] V= new double[3];
		double [] V2= new double[3];
		int debugIndex=(height/2)*width+ (width/2);
		int debugIndex1=(height/2)*width+ (width/4);
		double [] alphaStation=new double [numStations];
		int zIndex=fixXY?0:2;
		int numVariables=numberOfZGroups+zIndex;
		double [][] aM=new double [numVariables][numVariables];
		double [][] aB=new double [numVariables][1];
		double []   zPerStation= new double [numStations];

		for (int v=0;v<height;v++) for (int u=0;u<width; u++){
			int index=u+v*width;
			boolean thisDebug=(this.debugLevel>debugThreshold) && ((index==debugIndex) || (index==debugIndex1));
			if (thisDebug) System.out.println("calculateGridXYZCorr3D() debug("+this.debugLevel+"): index="+index+" v="+v+" u="+u);
			for (int i=0;i<numVariables;i++){
				aB[i][0]=0.0;
				for (int j=0;j<numVariables;j++) aM[i][j]=0.0;
			}
			for (int i=0;i<numStations;i++) alphaStation[i]=0.0;
			double alpha=0.0;
			boolean fallBack2D=true;
			if (grid3DCorrection) {
				for (int imgNum=0;imgNum<selectedImages.length;imgNum++) {
					int station=fittingStrategy.distortionCalibrationData.gIP[imgNum].getStationNumber(); // number of sub-camera
					int zGroup=stationGroups[station];
					if ((gridCorr2d[imgNum]!=null)  && (gridCorr2d[imgNum][3][index]>0.0) && (zGroup>=0)) {
						zPerStation[station]=gridCorr2d[imgNum][2][index]; // should all be the same for the same station
						// calculate unity vector from the camera lens to the grid point
						double absV=0.0;
						for (int i=0;i<V.length;i++){
							V[i]=patternParameters.gridGeometry[v][u][i]+gridCorr2d[imgNum][i][index]-cameraXYZ[imgNum][i]; // corrected value, including zCorr
							absV+=V[i]*V[i];
						}
						absV=Math.sqrt(absV);
						if (absV>0) for (int i=0;i<V.length;i++) V[i]/=absV;
						for (int i=0;i<V.length;i++) V2[i]=V[i]*V[i];
						if (thisDebug) System.out.println(" imgNum="+imgNum+" V[0]="+IJ.d2s(V[0],4)+" V[1]="+IJ.d2s(V[1],4)+" V[2]="+IJ.d2s(V[2],4)+
								" V2[0]="+IJ.d2s(V2[0],4)+" V2[1]="+IJ.d2s(V2[1],4)+" V2[2]="+IJ.d2s(V2[2],4));
						// When performin 3-d correction (X,Y,Z) the result point has to have minimal weighted sum of squared distances to all rays
// when summing for different stations, multiply W by sign(image belongs to station)
/*
Px, Py - calculated correction for individual image
V={Vx,Vy,Vz} unity vector from the camera lens center to the {Px,Py,0}
A - vector from the {Px,Py,0} to {X,Y,Z} = {X-Px,Y-Py,Z}
Projection of A on V will have length of A(.)V, Vector B=V*(A(.)V)
Vector D=A-B = A - V*(A(.)V)
D2=D(.)D= A(.)A - 2* (A(.)V ) * (A(.)V ) + (A(.)V ) * (A(.)V ) = A(.)A -  (A(.)V ) * (A(.)V )
D2=A(.)A -  (A(.)V )^2

A(.)A=(X-Px)^2 + (Y-Py)^2 + Z^2 =X^2 -2*X*Px +Px^2 +Y^2 -2*Y*Py +Py^2 +Z^2
A(.)A=X^2 -2*X*Px +Px^2 +Y^2 -2*Y*Py +Py^2 +Z^2
A(.)V=      (X-Px)*Vx + (Y-Py)*Vy + Z*Vz
(A(.)V)^2= ((X-Px)*Vx + (Y-Py)*Vy + Z*Vz)^2 = ((X-Px)*Vx)^2 + ((Y-Py)*Vy)^2 + (Z*Vz)^2 + 2*((X-Px)*Vx)*((Y-Py)*Vy)+ 2*((X-Px)*Vx)*(Z*Vz)+2*((Y-Py)*Vy)*(Z*Vz)
(A(.)V)^2= X^2*Vx^2 +Px^2*Vx^2 - 2*X*Px*Vx^2 +Y^2*Vy^2+Py^2*Vy^2-2*Y*Py*Vy^2 +Z^2*Vz^2 +2*X*Y*Vx*Vy +2*Px*Py*Vx*Vy - 2*X*Py*Vx*Vy - 2*Y*Px*Vx*Vy +2*X*Z*Vx*Vz - 2*Z*Px*Vx*Vz +2*Y*Z*Vy*Vz -2*z*Py*Vy*Vz

D2=
  +X^2 - X^2*Vx^2
  +Y^2 - Y^2*Vy^2
  +Z^2 - Z^2*Vz^2
-2*X*Y* Vx*Vy
-2*X*Z* Vx*Vz
-2*Y*Z* Vy*Vz
-2*X*Px +2*X*Px*Vx^2+ 2*X*Py*Vx*Vy
-2*Y*Py +2*Y*Py*Vy^2+ 2*Y*Px*Vx*Vy
+2*Z*Px*Vx*Vz   +2*Z*Py*Vy*Vz
+Px^2  +Py^2 -Px^2*Vx^2   -Py^2*Vy^2    -2*Px*Py*Vx*Vy

0= dD2/dX/2= X*(1-Vx^2) - Y* Vx*Vy - Z* Vx*Vz -Px + Px*Vx^2  + Py*Vx*Vy
0= dD2/dY/2= Y*(1-Vy^2) - X* Vx*Vy - Z* Vy*Vz -Py + Py*Vy^2  + Px*Vx*Vy
0= dD2/dZ/2= Z*(1-Vz^2) - X* Vx*Vz - Y* Vy*Vz     + Px*Vx*Vz + Py*Vy*Vz


 X*(Vx^2-1) + Y* (Vx*Vy)  + Z* (Vx*Vz)   =  Px * (Vx^2-1)  + Py* (Vx*Vy)
 X*(Vx*Vy)  + Y* (Vy^2-1) + Z* (Vy*Vz)   =  Px * (Vx*Vy)   + Py * (Vy^2-1)
 X*(Vx*Vz)  + Y* (Vy*Vz)  + Z* (Vz^2-1)  =  Px * (Vx*Vz)   + Py* (Vy*Vz)


 */

//   | sum(Wi*Cxi),  sum(Wi*Cxyi), sum(Wi*Cxzi) |
//M= | sum(Wi*Cxyi), sum(Wi*Cyi ), sum(Wi*Cyzi) |
//   | sum(Wi*Cxzi), sum(Wi*Cyzi), sum(Wi*Czi ) |

//   | sum(Wi*(P0xi*Cxi + P0yi*Cxyi + P0zi*Cxzi)) |
//B= | sum(Wi*(P0yi*Cyi + P0xi*Cyxi + P0zi*Cyzi)) |
//   | sum(Wi*(P0zi*Czi + P0yi*Czyi + P0xi*Czxi)) |
/*
	X*(Vxi^2-1) + Y*(Vxi*Vyi) + Z*(Vxi*Vzi) = P0xi*(Vxi^2-1) +P0yi*(Vxi*Vyi) + P0zi*(Vxi*Vzi)
	X*(Vxi*Vyi) + Y*(Vyi^2-1) + Z*(Vyi*Vzi) = P0xi*(Vxi*Vyi) +P0yi*(Vyi^2-1) + P0zi*(Vyi*Vzi)
	X*(Vxi*Vzi) + Y*(Vxi*Vyi) + Z*(Vzi^2-1) = P0xi*(Vxi*Vzi) +P0yi*(Vxi*Vyi) + P0zi*(Vzi^2-1)

	X*Cx  + Y*Cxy + Z*Cxz = P0xi*Cx  +P0yi*Cxy + P0zi*Cxz
	X*Cxy + Y*Cy  + Z*Cyz = P0xi*Cxy +P0yi*Cy  + P0zi*Cyz
	X*Cxz + Y*Cyz + Z*Cz  = P0xi*Cxz +P0yi*Cyz + P0zi*Cz
	P0zi==0.0, so - now we'll use P0zi - difference from this station to average

	X*Cx  + Y*Cxy + Z*Cxz = P0xi*Cx  +P0yi*Cxy
	X*Cxy + Y*Cy  + Z*Cyz = P0xi*Cxy +P0yi*Cy
	X*Cxz + Y*Cyz + Z*Cz  = P0xi*Cxz +P0yi*Cyz

*/
						Cx=V2[0]-1.0;
						Cy=V2[1]-1.0;
						Cz=V2[2]-1.0;
						Cxy= V[0]*V[1];
						Cxz= V[0]*V[2];
						Cyz= V[1]*V[2];
						if (thisDebug) System.out.println(" Cx="+IJ.d2s(Cx,6)+" Cy="+IJ.d2s(Cy,6)+" Cz="+IJ.d2s(Cz,6)+
								" Cxy="+IJ.d2s(Cxy,6)+" Cxz="+IJ.d2s(Cxz,6)+" Cyz="+IJ.d2s(Cyz,6));


						double W=gridCorr2d[imgNum][3][index];
						double Px=gridCorr2d[imgNum][0][index];
						double Py=gridCorr2d[imgNum][1][index];
						double Pz=gridCorr2d[imgNum][2][index];
						alpha+=W;
						alphaStation[station]+=W;
						if (thisDebug) System.out.println(imgNum+": Px="+IJ.d2s(Px,6)+" Py="+IJ.d2s(Py,6)+" W="+IJ.d2s(W,6));
						if (zIndex>0){ // X,Y correction is enabled, not only Z
							aM[0][0]+=W*Cx;
							aM[0][1]+=W*Cxy;
							aM[1][1]+=W*Cy;
							aM[0][2+zGroup]+=W*Cxz;
							aM[1][2+zGroup]+=W*Cyz;
							aB[0][0]+=W*(Px*Cx  + Py*Cxy + Pz*Cxz);
							aB[1][0]+=W*(Px*Cxy + Py*Cy  + Pz*Cyz);
						}
						aM[zIndex+zGroup][zIndex+zGroup]+=W*(Cz-variationPenalty); // -1>>Cz<0
						aB[zIndex+zGroup][0]+=W*(Px*Cxz + Py*Cyz + Pz*Cz);
					}
				}
				if (zIndex>0){// X,Y correction is enabled, not only Z
					aM[1][0]+=aM[0][1]; // why "+=" - just "="
					for (int zGroup=0;zGroup<numberOfZGroups;zGroup++){
						aM[zIndex+zGroup][0]+=aM[0][zIndex+zGroup];
						aM[zIndex+zGroup][1]+=aM[1][zIndex+zGroup];
					}
				}
				Matrix M=new Matrix(aM);
				Matrix B=new Matrix(aB);
				if (thisDebug) {
					System.out.println(" M:");
					M.print(8, 6);
					System.out.println(" B:");
					B.print(8, 6);
				}

				//			boolean fallBack2D=true;
				if ((new LUDecomposition(M)).isNonsingular()){
					double [] dXYZ=M.solve(B).getRowPackedCopy();
//// Now save per station group (with weights)
					if (zIndex>0){// X,Y correction is enabled, not only Z
						for (int i=0;i<2;i++) gridCorr3d[i][index]=dXYZ[i];
					}
					double zAverage=0.0;
					double sumW=0;
					for (int station=0;station<numStations;station++){
						double w=alphaStation[station];
						sumW+=w;
						gridZCorr3dWeight[station][index]=w;
						int zGroup=stationGroups[station];
						zAverage+=w*dXYZ[zIndex+zGroup];
					}
					if (sumW>0.0) {
						zAverage/=sumW;
						gridCorr3d[2][index]=zAverage; // weighted average of grid Z correction (from current pattern Z)
						gridCorr3d[3][index]=alpha; // same as sumW?
// second pass - calculate per-station Z corrections - referenced to existent current values
//zPerStation[station]
						for (int station=0;station<numStations;station++){
							int zGroup=stationGroups[station];
//							gridZCorr3d[station][index]=dXYZ[zIndex+zGroup]-zPerStation[station]; // differential from the current pattern geometry
							gridZCorr3d[station][index]=dXYZ[zIndex+zGroup]-zAverage; // differential from the current pattern geometry
						}
					}
					fallBack2D=false; //TODO:  make sure delta Z (Math.abs(gridCorr3d[2][index])) is not too big!!

					if (Math.abs(gridCorr3d[2][index])>maxZCorr) {
						fallBack2D=true; // temporary limit
					}
					if (thisDebug) System.out.println(" dX="+IJ.d2s(gridCorr3d[0][index],6)+" dY="+IJ.d2s(gridCorr3d[1][index],6)+" dZ="+IJ.d2s(gridCorr3d[2][index],6));
				}
			}
			if(fallBack2D && !(grid3DCorrection && noFallBack)) { // make a 2d averaging of weighted dx, dy correction - separately for each station group
				double [] gridZcorrPerGroup=      new double [numberOfZGroups];
				double [] gridZcorrAddPerGroup=      new double [numberOfZGroups];
				double [] gridZcorrWeightPerGroup=new double [numberOfZGroups];
				for (int i=0;i<numberOfZGroups;i++){
					gridZcorrPerGroup[i]=0.0;
					gridZcorrWeightPerGroup[i]=0.0;
					gridZcorrAddPerGroup[i]=0.0;
				}
				for (int imgNum=0;imgNum<selectedImages.length;imgNum++) {
					int station=fittingStrategy.distortionCalibrationData.gIP[imgNum].getStationNumber(); // number of sub-camera
					int zGroup=stationGroups[station];
					if ((gridCorr2d[imgNum]!=null)  && (gridCorr2d[imgNum][3][index]>0.0) && (zGroup>=0)) {
						double w=gridCorr2d[imgNum][3][index];
						double z=gridCorr2d[imgNum][2][index]; // difference from average Z
						gridZcorrPerGroup[zGroup]+=w*z;
						gridZcorrWeightPerGroup[zGroup]+=w;
					}
				}
				for (int i=0;i<numberOfZGroups;i++) if (gridZcorrWeightPerGroup[i]>0.0) gridZcorrPerGroup[i]/=gridZcorrWeightPerGroup[i];
				for (int i=0;i<gridCorr3d.length;i++) gridCorr3d[i][index]=0.0;

				double s=0;
				for (int imgNum=0;imgNum<selectedImages.length;imgNum++) {
					int station=fittingStrategy.distortionCalibrationData.gIP[imgNum].getStationNumber(); // number of sub-camera
					int zGroup=stationGroups[station];
					if ((gridCorr2d[imgNum]!=null)  && (gridCorr2d[imgNum][3][index]>0.0) && (zGroup>=0)) {
						double z=patternParameters.gridGeometry[v][u][2]+gridZcorrPerGroup[zGroup];
						double [] cv={
								patternParameters.gridGeometry[v][u][0]-cameraXYZ[imgNum][0],
								patternParameters.gridGeometry[v][u][1]-cameraXYZ[imgNum][1],
								z-cameraXYZ[imgNum][2]};
						double cv2=cv[0]*cv[0]+cv[1]*cv[1]+cv[2]*cv[2];
						double acv=Math.sqrt(cv2);
						for (int i=0;i<3;i++)cv[i]/=acv; // make unity vector;
						// intersection of the corrected view ray with the average taget plane
						double [] dXYplane0={-gridZcorrPerGroup[zGroup]/cv[2]*cv[0],-gridZcorrPerGroup[zGroup]/cv[2]*cv[1]};
						double [] modCorrXY={gridCorr2d[imgNum][0][index]+dXYplane0[0], gridCorr2d[imgNum][1][index]+dXYplane0[1]};
						double kv=(modCorrXY[0]*cv[0]+modCorrXY[1]*cv[1])/cv2;
						double w=gridCorr2d[imgNum][3][index];
						gridCorr3d[0][index]+=w*(gridCorr2d[imgNum][0][index]-cv[0]*kv);
						gridCorr3d[1][index]+=w*(gridCorr2d[imgNum][1][index]-cv[1]*kv);
						gridZcorrAddPerGroup[zGroup]+=w*(                            -cv[2]*kv);
// not finished per station/per group 2d correction, will just use corerction average
						gridCorr3d[2][index]+=w*(                            -cv[2]*kv);
						s+=w;
					}
				}
				for (int i=0;i<numberOfZGroups;i++) if (gridZcorrWeightPerGroup[i]>0.0) gridZcorrAddPerGroup[i]/=gridZcorrWeightPerGroup[i];
				for (int station=0;station<numStations;station++){
					int zGroup=stationGroups[station];
					gridZCorr3d[station][index]=gridZcorrAddPerGroup[zGroup]; // differential from the current pattern geometry
				}
				if (s>0){
					gridCorr3d[0][index]/=s;
					gridCorr3d[1][index]/=s;
					gridCorr3d[2][index]/=s;
				} else {
					gridCorr3d[0][index]=0.0;
					gridCorr3d[1][index]=0.0;
					gridCorr3d[2][index]=0.0;
				}
				gridCorr3d[3][index]=s;
				if (thisDebug) System.out.println(" Using 2d averaging: dX="+IJ.d2s(gridCorr3d[0][index],6)+
						" dY="+IJ.d2s(gridCorr3d[1][index],6)+" dZ="+IJ.d2s(gridCorr3d[2][index],6));
			}
		}
		// Make average correction zero is it needed?
		// create "reliable" mask for averaging/tilting - disregard the outmost grid pixels
		boolean [] reliable=new boolean [width*height];
		double wThreshold=0.0;
		for (int v=0;v<height;v++) for (int u=0;u<width;u++){
			int index=u+v*width;
			reliable[index]=false;
			if ((v>0) && (u>0) && (v<(height-1)) && (u<(width-1)) &&
					(gridCorr3d[3][index]>wThreshold) &&
					(gridCorr3d[3][index-1]>wThreshold) &&
					(gridCorr3d[3][index+1]>wThreshold) &&
					(gridCorr3d[3][index-width]>wThreshold) &&
					(gridCorr3d[3][index+width]>wThreshold) ){
				reliable[index]=true;
			}

		}
		double corrAverage;
		for (int c=0;c<3;c++){
			corrAverage=0.0;
			double s=0.0;
			for (int i=0;i<gridCorr3d[0].length;i++) if (reliable[i]) {
				corrAverage+=gridCorr3d[c][i]*gridCorr3d[3][i];
				s+=gridCorr3d[3][i];
			}
			corrAverage/=s;
			//			System.out.println("zCorrAverage["+c+"="+corrAverage);
			for (int i=0;i<gridCorr3d[c].length;i++) gridCorr3d[c][i]-=corrAverage;
		}
		// for Z correction compensate for x/y tilts
		String [] titles={"X-correction(mm)","Y-correction(mm)","Z-correction","Weight"};
		if (rotateCorrection) {
			double SX=0.0,SX2=0.0,SZ=0.0,SXY=0.0,SXZ=0.0,S0=0.0,SY=0.0,SY2=0.0,SYZ=0.0;
			double [][] gridGeom=new double [3][gridCorr3d[0].length];
			for (int c=0;c<gridGeom.length;c++) for (int i=0;i<gridGeom[c].length;i++)gridGeom[c][i]=0.0;

			for (int v=0;v<height;v++) for (int u=0;u<width; u++){
				int index=u+v*width;
				double W=gridCorr3d[3][index];
				gridGeom[0][index]=patternParameters.gridGeometry[v][u][0];
				gridGeom[1][index]=patternParameters.gridGeometry[v][u][1];
				gridGeom[2][index]=W;
				if ((reliable[index]) && (W>0.0)){
					S0+=W;
					SX+=  W*patternParameters.gridGeometry[v][u][0];
					SX2+= W*patternParameters.gridGeometry[v][u][0]*patternParameters.gridGeometry[v][u][0];
					SY+=  W*patternParameters.gridGeometry[v][u][1];
					SY2+= W*patternParameters.gridGeometry[v][u][1]*patternParameters.gridGeometry[v][u][1];
					SXY+= W*patternParameters.gridGeometry[v][u][0]*patternParameters.gridGeometry[v][u][1];
					SZ+=  W*gridCorr3d[2][index];
					SXZ+= W*gridCorr3d[2][index]*patternParameters.gridGeometry[v][u][0];
					SYZ+= W*gridCorr3d[2][index]*patternParameters.gridGeometry[v][u][1];
				}
			}
			double [][] aM1= {
					{SX2, SXY, SX},
					{SXY, SY2, SY},
					{SX,  SY,  S0}};
			double [][] aB1= {{SXZ},{SYZ},{SZ}};
			Matrix M=new Matrix(aM1);
			Matrix B=new Matrix(aB1);
			if (this.debugLevel>2) {
				System.out.println(" M:");
				M.print(8, 6);
				System.out.println(" B:");
				B.print(8, 6);
				System.out.println(" Ax,Ay,B:");
				M.solve(B).print(8, 6);
			}
			double [] tilts=M.solve(B).getRowPackedCopy(); // singular ???
			if (this.debugLevel>2) {
				if (this.refineParameters.showThisCorrection) {
					ShowDoubleFloatArrays.showArrays(gridCorr3d, getGridWidth(), getGridHeight(),  true, "before tilt:", titles);
				}
			}
			for (int v=0;v<height;v++) for (int u=0;u<width; u++){
				int index=u+v*width;
				gridCorr3d[2][index]-=tilts[0]*patternParameters.gridGeometry[v][u][0]+tilts[1]*patternParameters.gridGeometry[v][u][1]+tilts[2];
			}
		}
    	if (this.debugLevel>2) {
    		if (this.refineParameters.showThisCorrection) {
    			double [][] gridCorr3dClone=new double [4][width*height];
    			for (int c=0;c<gridCorr3dClone.length;c++) for (int i=0;i<gridCorr3dClone[c].length;i++)
    				gridCorr3dClone[c][i]=reliable[i]? gridCorr3d[c][i]:0.0;
    			ShowDoubleFloatArrays.showArrays(gridCorr3dClone, getGridWidth(), getGridHeight(),  true, "after tilt:", titles);
    		}
    	}
    	IJ.showStatus("");

// combine in a single array?

    	double [][][] result={gridCorr3d,gridZCorr3d,gridZCorr3dWeight};
		return  result;
	}
	public double [][][] calculateGridXYZCorr2D(
			final int width,
			final int height,
			final int [] stationGroups,
			final boolean [] selectedImages,
			final double [][] cameraXYZ,
			final LensDistortionParameters lensDistortionParametersProto,
			final boolean showIndividual,
			final int threadsMax,
			final boolean updateStatus
			){
//		final boolean isTripod=this.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.is_tripod;
//		final boolean cartesian=this.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.cartesian;
		final int [][] dirs=            {{0,0},{-1,0},{1,0},{0,-1},{0,1}}; // possible to make 8 directions
		final double [][][] derivatives={ // for of /du, /dv 3 variants, depending on which neighbors are available
				{
					{ 0.0,-0.5, 0.5, 0.0, 0.0},
					{ 1.0,-1.0, 0.0, 0.0, 0.0},
					{-1.0, 0.0, 1.0, 0.0, 0.0}
				},
				{
					{ 0.0, 0.0, 0.0,-0.5, 0.5},
					{ 1.0, 0.0, 0.0,-1.0, 0.0},
					{-1.0, 0.0, 0.0, 0.0, 1.0}}};
		final double [][][] gridCorr2d=new double [selectedImages.length][][]; // per-image grid {dx,dy,weight} corrections
		for (int i=0;i<gridCorr2d.length;i++) {
			gridCorr2d[i]=null;
			cameraXYZ[i]=null;
		}
		// Should it be just once - common for all images? (removed from the "for" loop)
		final double [] diff=calcYminusFx(this.currentfX);
		final int debugLevel=this.debugLevel;
		final int [] imageStartIndex=this.imageStartIndex;
		final double [] Y=this.Y;
		final double [] weightFunction= this.weightFunction;
   		final Thread[] threads = newThreadArray(threadsMax);
   		final AtomicInteger imageNumberAtomic = new AtomicInteger(0);
   		final AtomicInteger imageFinishedAtomic = new AtomicInteger(0);
   		final double [] progressValues=new double [selectedImages.length];
   		int numSelectedImages=0;
   		for (int i=0;i<selectedImages.length;i++) if (selectedImages[i]) numSelectedImages++;
   		int selectedIndex=0;
   		for (int i=0;i<selectedImages.length;i++) {
   			progressValues[i]=(selectedIndex+1.0)/numSelectedImages;
   			if (selectedImages[i]) selectedIndex++;
   			if (selectedIndex>=numSelectedImages) selectedIndex--;
   		}
		IJ.showStatus("Calculating pattern geometry correction...");
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					LensDistortionParameters lensDistortionParameters=lensDistortionParametersProto.clone(); // see - if that is needed - maybe new is OK
					//   					LensDistortionParameters lensDistortionParameters= new LensDistortionParameters();
					for (int imgNum=imageNumberAtomic.getAndIncrement(); imgNum<selectedImages.length; imgNum=imageNumberAtomic.getAndIncrement()) if (selectedImages[imgNum]){
						//		for (int imgNum=0;imgNum<selectedImages.length;imgNum++) if (selectedImages[imgNum]) {
						int station=fittingStrategy.distortionCalibrationData.gIP[imgNum].getStationNumber(); // number of sub-camera
						if (stationGroups[station]<0) continue; // do not process images that do not belong to selected stations
						gridCorr2d[imgNum]=new double [4][width*height]; // dx, dy only - added zCorr per station
						for (int n=0;n<gridCorr2d[imgNum].length;n++) for (int i=0;i<gridCorr2d[imgNum][0].length;i++) gridCorr2d[imgNum][n][i]=0.0;
						//		int chnNum=fittingStrategy.distortionCalibrationData.gIP[imgNum].channel; // number of sub-camera
						cameraXYZ[imgNum]=new double[3];
						// The following method sets this.lensDistortionParameters and invokes this.lensDistortionParameters.recalcCommons();
						lensDistortionParameters.lensCalcInterParamers(
								lensDistortionParameters,
								fittingStrategy.distortionCalibrationData.isTripod(),
								fittingStrategy.distortionCalibrationData.isCartesian(),
					    		fittingStrategy.distortionCalibrationData.getPixelSize(imgNum),
					    		fittingStrategy.distortionCalibrationData.getDistortionRadius(imgNum),
								null, //this.interParameterDerivatives, // [22][]
								fittingStrategy.distortionCalibrationData.getParameters(imgNum), // 22-long parameter vector for the image
								null); // if no derivatives, null is OK
						//				false); // calculate this.interParameterDerivatives -derivatives array (false - just this.values)
						cameraXYZ[imgNum]=lensDistortionParameters.getLensCenterCoordinates();
						if (debugLevel>2) {
							System.out.println("calculateGridXYZCorr(): imgNum="+imgNum+" lens coordinates (mm)={"+
									IJ.d2s(cameraXYZ[imgNum][0],3)+", "+IJ.d2s(cameraXYZ[imgNum][1],3)+", "+IJ.d2s(cameraXYZ[imgNum][2],3)+"}");
						}
						//		double [] diff=calcYminusFx(this.currentfX); // removed from the loop
						// find data range for the selected image
						int index=imageStartIndex[imgNum]; // set when fitting series is init
						double [][] imgData=new double[showIndividual?7:5][getGridHeight() * width]; // dPX, dPY, Px, Py, alpha
						for (int i=0;i<imgData.length;i++) for (int j=0;j<imgData[i].length;j++)imgData[i][j]=0.0;
						// first pass - prepare [v][u]arrays
						for (int i=0;i<fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV.length;i++){
							int u=fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV[i][0]+patternParameters.U0;
							int v=fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV[i][1]+patternParameters.V0;
							int vu=u+width*v;
							imgData[0][vu]= diff[2*(index+i)]; // out of bound 1410
							imgData[1][vu]= diff[2*(index+i)+1];
							imgData[2][vu]= Y[2*(index+i)];  // measured pixel x
							imgData[3][vu]= Y[2*(index+i)+1];// measured pixel y

							//				imgData[4][vu]= fittingStrategy.distortionCalibrationData.getMask(chnNum, imgData[2][vu], imgData[3][vu]);

							if (weightFunction!=null) {
								imgData[4][vu]= weightFunction[2*(index+i)];
							} else {
								imgData[4][vu]= 1.0;
							}
						}
						// second pass - calculate derivatives, and residuals in mm
						for (int i=0;i<fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV.length;i++){
							int u=fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV[i][0]+patternParameters.U0;
							int v=fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV[i][1]+patternParameters.V0;
							int vu=u+width*v;
							gridCorr2d[imgNum][0][vu]=0.0;
							gridCorr2d[imgNum][1][vu]=0.0;
							gridCorr2d[imgNum][2][vu]=patternParameters.getZCorr(vu,station); // per-station Z correction from average
							gridCorr2d[imgNum][3][vu]=0.0; // weight
							double [][] gXY =new double[dirs.length][3];
							//			double [][] gpXY=new double[dirs.length][2];
							double [][] gpXY=new double[dirs.length][3];
							boolean [] dirMask=new boolean [dirs.length];
							for (int dir=0;dir<dirs.length;dir++){
								int u1=u+dirs[dir][0];
								int v1=v+dirs[dir][1];
								int vu1=u1+width*v1;
								dirMask[dir] = (u1>=0) && (v1>=0) && (u1<width) && (v1<height) && (imgData[4][vu1]>0);
								if (dirMask[dir]){
									gXY[dir][0]= patternParameters.gridGeometry[v1][u1][0];
									gXY[dir][1]= patternParameters.gridGeometry[v1][u1][1];
									gXY[dir][2]= patternParameters.gridGeometry[v1][u1][2]; // Here - average Z
									gpXY[dir][0]=imgData[2][vu1];
									gpXY[dir][1]=imgData[3][vu1];
								} else {
									gXY[dir][0]= 0.0;
									gXY[dir][1]= 0.0;
									gXY[dir][2]= 0.0;
									gpXY[dir][0]=0.0;
									gpXY[dir][1]=0.0;

								}
							}
							int [] variants={-1,-1}; // {horizontal, vertical}
							boolean variantsExist=true;
							for (int duv=0;duv<2;duv++){ // 0 - horizontal, 1 - vertical
								for (int variant=0;variant<derivatives[duv].length;variant++) { // variants: 0 half of right/left, 1 left deriv, 2 - right deriv
									boolean fit=true;
									for (int dir=0;dir<dirs.length;dir++) if ((derivatives[duv][variant][dir]!=0) && !dirMask[dir]){
										fit=false;
										break;
									}
									if (fit) {
										variants[duv]=variant;
										break;
									}
								}
								if (variants[duv]<0) { // could not find any variant to calculate derivatives for this direction
									variantsExist=false;
									break;
								}
							}
							if (!variantsExist){
								imgData[4][vu]=0.0;
								continue;
							}
							double [][] dXY_dUV= new double [2][2];
							double [][] dpXY_dUV=new double [2][2];
							for (int nom=0;nom<2;nom++) { // 0-x, 1 - y
								for (int denom=0;denom<2;denom++) { //0 - du, 1 - dv
									dXY_dUV [nom][denom]=0.0;
									dpXY_dUV[nom][denom]=0.0;
									for (int dir=0;dir<dirs.length;dir++){
										dXY_dUV [nom][denom]+=gXY [dir][nom]*derivatives[denom][variants[denom]][dir];
										dpXY_dUV[nom][denom]+=gpXY[dir][nom]*derivatives[denom][variants[denom]][dir];
									}
								}
							}
							double [] dpXY={imgData[0][vu],imgData[1][vu]};
							Matrix MdpXY=    new Matrix(dpXY,2); // 2 rows
							Matrix MdXY_dUV= new Matrix(dXY_dUV);
							Matrix MdpXY_dUV=new Matrix(dpXY_dUV);
							if ((new LUDecomposition(MdpXY_dUV)).isNonsingular()){
								/*
								 * MdpXY= MdpXY_dUV* MdUV
								 * MdXY=  MdXY_dUV * MdUV
								 * MdUV=  MdpXY_dUV.solve(MdpXY);
								 * MdXY=  MdXY_dUV * MdpXY_dUV.solve(MdpXY);
								 */
								Matrix MdXY=MdXY_dUV.times(MdpXY_dUV.solve(MdpXY));
								double [] dXY=MdXY.getRowPackedCopy();
								gridCorr2d[imgNum][0][vu]=dXY[0];
								gridCorr2d[imgNum][1][vu]=dXY[1];
								gridCorr2d[imgNum][3][vu]=imgData[4][vu]; // weight
							}
						} // end scanning pixels
						if (showIndividual) {
							String [] titles={"diff-X","diff-Y","pX","pY","alpha","X-correction(mm)","Y-correction(mm)","Z-correction(mm)"};
							ShowDoubleFloatArrays.showArrays(imgData, width, height,  true, "Grid"+imgNum, titles);
						}
   						final int numFinished=imageFinishedAtomic.getAndIncrement();
//						IJ.showProgress(progressValues[numFinished]);
						SwingUtilities.invokeLater(new Runnable() {
							@Override
							public void run() {
								// Here, we can safely update the GUI
								// because we'll be called from the
								// event dispatch thread
								IJ.showProgress(progressValues[numFinished]);
							}
						});
					}
				}
			};
		}
		startAndJoin(threads);

		IJ.showProgress(1.0);
		return gridCorr2d;
	}





	/**
	 * Calculates corrections to X and Y coordinates of the grid nodes
	 * //@param distortionCalibrationData - used to receive sensor mask(s)
	 * @param grid3DCorrection - if true - calculate 3d correction, false - slow 3d (2d perpendicular to view)
	 * @param maxZCorr - maximal allowed correction in Z-direction (if wants more, will fall back to 2-d correction (perpendicular to the view)
	 * @param showIndividual - show individual images
	 * @return first index - 0 - correction x (mm), 1 - correction y(mm), 2 - correction z(mm)  3 - weight (number of images used)
	 */

	public double [][] calculateGridXYZCorr3D( // old version
//			DistortionCalibrationData distortionCalibrationData,
			boolean grid3DCorrection,
			boolean rotateCorrection,
			double maxZCorr,
			boolean showIndividual){
		int width=getGridWidth();
		int height=getGridHeight();
		int [][] dirs=            {{0,0},{-1,0},{1,0},{0,-1},{0,1}}; // possible to make 8 directions
		double [][][] derivatives={ // for of /du, /dv 3 variants, depending on which neighbors are available
				{
					{ 0.0,-0.5, 0.5, 0.0, 0.0},
					{ 1.0,-1.0, 0.0, 0.0, 0.0},
					{-1.0, 0.0, 1.0, 0.0, 0.0}
				},
				{
					{ 0.0, 0.0, 0.0,-0.5, 0.5},
					{ 1.0, 0.0, 0.0,-1.0, 0.0},
					{-1.0, 0.0, 0.0, 0.0, 1.0}}};
		boolean [] selectedImages=fittingStrategy.selectedImages();
		double [][][] gridCorr2d=new double [selectedImages.length][][]; // per-image grid {dx,dy,weight} corrections

		double [][] cameraXYZ=new double [selectedImages.length][];
		for (int i=0;i<gridCorr2d.length;i++) {
			gridCorr2d[i]=null;
			cameraXYZ[i]=null;
		}
		int numSelected=0;
		for (int imgNum=0;imgNum<selectedImages.length;imgNum++) if (selectedImages[imgNum]) numSelected++;
		int numProcessed=0;
		IJ.showStatus("Calculating pattern geometry correction...");
		for (int imgNum=0;imgNum<selectedImages.length;imgNum++) if (selectedImages[imgNum]) {
			gridCorr2d[imgNum]=new double [3][width*height]; // dx, dy only
			for (int n=0;n<gridCorr2d[imgNum].length;n++) for (int i=0;i<gridCorr2d[imgNum][0].length;i++) gridCorr2d[imgNum][n][i]=0.0;
//			int chnNum=fittingStrategy.distortionCalibrationData.gIP[imgNum].channel; // number of sub-camera
			cameraXYZ[imgNum]=new double[3];
			// The following method sets this.lensDistortionParameters and invokes this.lensDistortionParameters.recalcCommons();
			this.lensDistortionParameters.lensCalcInterParamers(
					this.lensDistortionParameters,
					fittingStrategy.distortionCalibrationData.isTripod(),
					fittingStrategy.distortionCalibrationData.isCartesian(),
		    		fittingStrategy.distortionCalibrationData.getPixelSize(imgNum),
		    		fittingStrategy.distortionCalibrationData.getDistortionRadius(imgNum),
					null, //this.interParameterDerivatives, // [22][]
//					fittingStrategy.distortionCalibrationData.pars[imgNum], // 22-long parameter vector for the image
					fittingStrategy.distortionCalibrationData.getParameters(imgNum), // 22-long parameter vector for the image
					null); // if no derivatives, null is OK
//					false); // calculate this.interParameterDerivatives -derivatives array (false - just this.values)
			cameraXYZ[imgNum]=lensDistortionParameters.getLensCenterCoordinates();
			if (this.debugLevel>2) {
				System.out.println("calculateGridXYZCorr(): imgNum="+imgNum+" lens coordinates (mm)={"+
						IJ.d2s(cameraXYZ[imgNum][0],3)+", "+IJ.d2s(cameraXYZ[imgNum][1],3)+", "+IJ.d2s(cameraXYZ[imgNum][2],3)+"}");
			}
			double [] diff=calcYminusFx(this.currentfX);
			// find data range for the selected image
			int index=this.imageStartIndex[imgNum]; // set when fitting series is init
/*
			int index=0;
			int numImg=fittingStrategy.distortionCalibrationData.getNumImages();

			for (int iNum=0;(iNum<imgNum) && (iNum<numImg) ;iNum++) if (selectedImages[iNum]) //
				index+=fittingStrategy.distortionCalibrationData.gIP[iNum].pixelsUV.length;
			//System.out.println ("+++++++++++++imgNum="+imgNum+" index="+index);
*/
			if (this.debugLevel>2) {
				System.out.println("calculateGridXYZCorr(): fX.length="+this.currentfX.length+" this image index="+index);
			}
			double [][] imgData=new double[showIndividual?7:5][getGridHeight() * width]; // dPX, dPY, Px, Py, alpha
			for (int i=0;i<imgData.length;i++) for (int j=0;j<imgData[i].length;j++)imgData[i][j]=0.0;
			// first pass - prepare [v][u]arrays
			for (int i=0;i<fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV.length;i++){
				int u=fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV[i][0]+patternParameters.U0;
				int v=fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV[i][1]+patternParameters.V0;
				int vu=u+width*v;
				imgData[0][vu]=   diff[2*(index+i)]; // out of bound 1410
				imgData[1][vu]=   diff[2*(index+i)+1];
				imgData[2][vu]= this.Y[2*(index+i)];  // measured pixel x
				imgData[3][vu]= this.Y[2*(index+i)+1];// measured pixel y

				//				imgData[4][vu]= fittingStrategy.distortionCalibrationData.getMask(chnNum, imgData[2][vu], imgData[3][vu]);

				if (this.weightFunction!=null) {
					imgData[4][vu]= this.weightFunction[2*(index+i)];
				} else {
					imgData[4][vu]= 1.0;
				}
				//				if (imgNum==1) System.out.println ("---index="+index+" i="+i+" vu="+vu+ " v="+v+" u="+u+" x="+IJ.d2s(this.Y[2*(index+i)],1)+" y="+IJ.d2s(this.Y[2*(index+i)+1],1));
			}
			// second pass - calculate derivatives, and residuals in mm
			for (int i=0;i<fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV.length;i++){
				int u=fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV[i][0]+patternParameters.U0;
				int v=fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV[i][1]+patternParameters.V0;
				int vu=u+width*v;
				gridCorr2d[imgNum][0][vu]=0.0;
				gridCorr2d[imgNum][1][vu]=0.0;
				gridCorr2d[imgNum][2][vu]=0.0; // weight

				double [][] gXY =new double[dirs.length][3];
				double [][] gpXY=new double[dirs.length][2];
				boolean [] dirMask=new boolean [dirs.length];
				for (int dir=0;dir<dirs.length;dir++){
					int u1=u+dirs[dir][0];
					int v1=v+dirs[dir][1];
					int vu1=u1+width*v1;
					dirMask[dir] = (u1>=0) && (v1>=0) && (u1<width) && (v1<height) && (imgData[4][vu1]>0);
					gXY[dir][0]= dirMask[dir]?patternParameters.gridGeometry[v1][u1][0]:0.0;
					gXY[dir][1]= dirMask[dir]?patternParameters.gridGeometry[v1][u1][1]:0.0;
					gXY[dir][2]= dirMask[dir]?patternParameters.gridGeometry[v1][u1][2]:0.0; // Add per-station (optionally)
					gpXY[dir][0]=dirMask[dir]?imgData[2][vu1]:0.0;
					gpXY[dir][1]=dirMask[dir]?imgData[3][vu1]:0.0;
				}
				int [] variants={-1,-1}; // {horizontal, vertical}
				boolean variantsExist=true;
				for (int duv=0;duv<2;duv++){ // 0 - horizontal, 1 - vertical
					for (int variant=0;variant<derivatives[duv].length;variant++) { // variants: 0 half of right/left, 1 left deriv, 2 - right deriv
						boolean fit=true;
						for (int dir=0;dir<dirs.length;dir++) if ((derivatives[duv][variant][dir]!=0) && !dirMask[dir]){
							fit=false;
							break;
						}
						if (fit) {
							variants[duv]=variant;
							break;
						}
					}
					if (variants[duv]<0) { // could not find any variant to calculate derivatives for this direction
						variantsExist=false;
						break;
					}
				}
				if (!variantsExist){
					imgData[4][vu]=0.0;
					continue;
				}
				double [][] dXY_dUV= new double [2][2];
				double [][] dpXY_dUV=new double [2][2];
				for (int nom=0;nom<2;nom++) { // 0-x, 1 - y
					for (int denom=0;denom<2;denom++) { //0 - du, 1 - dv
						dXY_dUV [nom][denom]=0.0;
						dpXY_dUV[nom][denom]=0.0;
						for (int dir=0;dir<dirs.length;dir++){
							dXY_dUV [nom][denom]+=gXY [dir][nom]*derivatives[denom][variants[denom]][dir];
							dpXY_dUV[nom][denom]+=gpXY[dir][nom]*derivatives[denom][variants[denom]][dir];
						}
					}

				}
				double [] dpXY={imgData[0][vu],imgData[1][vu]};
				Matrix MdpXY=    new Matrix(dpXY,2); // 2 rows
				Matrix MdXY_dUV= new Matrix(dXY_dUV);
				Matrix MdpXY_dUV=new Matrix(dpXY_dUV);
				if ((new LUDecomposition(MdpXY_dUV)).isNonsingular()){
					/*
					 * MdpXY= MdpXY_dUV* MdUV
					 * MdXY=  MdXY_dUV * MdUV
					 * MdUV=  MdpXY_dUV.solve(MdpXY);
					 * MdXY=  MdXY_dUV * MdpXY_dUV.solve(MdpXY);
					 */
					Matrix MdXY=MdXY_dUV.times(MdpXY_dUV.solve(MdpXY));
					double [] dXY=MdXY.getRowPackedCopy();
					gridCorr2d[imgNum][0][vu]=dXY[0];
					gridCorr2d[imgNum][1][vu]=dXY[1];
					gridCorr2d[imgNum][2][vu]=imgData[4][vu]; // weight
				}
			} // end scanning pixels
			if (showIndividual) {
		        String [] titles={"diff-X","diff-Y","pX","pY","alpha","X-correction(mm)","Y-correction(mm)","Z-correction(mm)"};
				ShowDoubleFloatArrays.showArrays(imgData, width, height,  true, "Grid"+imgNum, titles);
			}
			IJ.showProgress(++numProcessed,numSelected);
		}
		IJ.showProgress(1.0);
		IJ.showStatus("Calculating pattern 3d correction...");
// now using gridCorr2d[imgNum], cameraXYZ[imgNum] and patternParameters.gridGeometry[v][u] find the 3d correction     public double [][][] gridGeometry=null; // [v][u]{x,y,z,"alpha"} alpha=0 - no ghrid, 1 - grid
		double [][] gridCorr3d=new double [4][width*height];
		for (int n=0;n<gridCorr3d.length;n++) for (int i=0;i<gridCorr3d[0].length;i++) gridCorr3d[n][i]=0.0;
		double Cx,Cy,Cz,Cxy,Cxz,Cyz;
		double [] V= new double[3];
		double [] V2= new double[3];
		int debugIndex=(height/2)*width+ (width/2);
		int debugIndex1=(height/2)*width+ (width/4);
		for (int v=0;v<height;v++) for (int u=0;u<width; u++){
			int index=u+v*width;
			boolean thisDebug=(this.debugLevel>1) && ((index==debugIndex) || (index==debugIndex1));
			if (thisDebug) System.out.println("calculateGridXYZCorr3D() debug("+this.debugLevel+"): index="+index+" v="+v+" u="+u);
			double [][] aM={{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
			double [][] aB ={{0.0},{0.0},{0.0}};
			double alpha=0.0;
			boolean fallBack2D=true;
			if (grid3DCorrection) {
				for (int imgNum=0;imgNum<selectedImages.length;imgNum++)
					if ((gridCorr2d[imgNum]!=null)  && (gridCorr2d[imgNum][2][index]>0.0)) {
						// calculate unity vector from the camera lens to the grid point
						double absV=0.0;
						for (int i=0;i<V.length;i++){
							V[i]=patternParameters.gridGeometry[v][u][i]-cameraXYZ[imgNum][i];
							absV+=V[i]*V[i];
						}
						absV=Math.sqrt(absV);
						if (absV>0) for (int i=0;i<V.length;i++) V[i]/=absV;
						for (int i=0;i<V.length;i++) V2[i]=V[i]*V[i];
						if (thisDebug) System.out.println(" imgNum="+imgNum+" V[0]="+IJ.d2s(V[0],4)+" V[1]="+IJ.d2s(V[1],4)+" V[2]="+IJ.d2s(V[2],4)+
								" V2[0]="+IJ.d2s(V2[0],4)+" V2[1]="+IJ.d2s(V2[1],4)+" V2[2]="+IJ.d2s(V2[2],4));
						// When performin 3-d correction (X,Y,Z) the result point has to have minimal weighted sum of squared distances to all rays
// when summing for different stations, multiply W by sign(image belongs to station)
/*
Px, Py - calculated correction for individual image
V={Vx,Vy,Vz} unity vector from the camera lens center to the {Px,Py,0}
A - vector from the {Px,Py,0} to {X,Y,Z} = {X-Px,Y-Py,Z}
Projection of A on V will have length of A(.)V, Vector B=V*(A(.)V)
Vector D=A-B = A - V*(A(.)V)
D2=D(.)D= A(.)A - 2* (A(.)V ) * (A(.)V ) + (A(.)V ) * (A(.)V ) = A(.)A -  (A(.)V ) * (A(.)V )
D2=A(.)A -  (A(.)V )^2

A(.)A=(X-Px)^2 + (Y-Py)^2 + Z^2 =X^2 -2*X*Px +Px^2 +Y^2 -2*Y*Py +Py^2 +Z^2
A(.)A=X^2 -2*X*Px +Px^2 +Y^2 -2*Y*Py +Py^2 +Z^2
A(.)V=      (X-Px)*Vx + (Y-Py)*Vy + Z*Vz
(A(.)V)^2= ((X-Px)*Vx + (Y-Py)*Vy + Z*Vz)^2 = ((X-Px)*Vx)^2 + ((Y-Py)*Vy)^2 + (Z*Vz)^2 + 2*((X-Px)*Vx)*((Y-Py)*Vy)+ 2*((X-Px)*Vx)*(Z*Vz)+2*((Y-Py)*Vy)*(Z*Vz)
(A(.)V)^2= X^2*Vx^2 +Px^2*Vx^2 - 2*X*Px*Vx^2 +Y^2*Vy^2+Py^2*Vy^2-2*Y*Py*Vy^2 +Z^2*Vz^2 +2*X*Y*Vx*Vy +2*Px*Py*Vx*Vy - 2*X*Py*Vx*Vy - 2*Y*Px*Vx*Vy +2*X*Z*Vx*Vz - 2*Z*Px*Vx*Vz +2*Y*Z*Vy*Vz -2*z*Py*Vy*Vz

D2=
  +X^2 - X^2*Vx^2
  +Y^2 - Y^2*Vy^2
  +Z^2 - Z^2*Vz^2
-2*X*Y* Vx*Vy
-2*X*Z* Vx*Vz
-2*Y*Z* Vy*Vz
-2*X*Px +2*X*Px*Vx^2+ 2*X*Py*Vx*Vy
-2*Y*Py +2*Y*Py*Vy^2+ 2*Y*Px*Vx*Vy
+2*Z*Px*Vx*Vz   +2*Z*Py*Vy*Vz
+Px^2  +Py^2 -Px^2*Vx^2   -Py^2*Vy^2    -2*Px*Py*Vx*Vy

0= dD2/dX/2= X*(1-Vx^2) - Y* Vx*Vy - Z* Vx*Vz -Px + Px*Vx^2  + Py*Vx*Vy
0= dD2/dY/2= Y*(1-Vy^2) - X* Vx*Vy - Z* Vy*Vz -Py + Py*Vy^2  + Px*Vx*Vy
0= dD2/dZ/2= Z*(1-Vz^2) - X* Vx*Vz - Y* Vy*Vz     + Px*Vx*Vz + Py*Vy*Vz


 X*(Vx^2-1) + Y* (Vx*Vy)  + Z* (Vx*Vz)   =  Px * (Vx^2-1)  + Py* (Vx*Vy)
 X*(Vx*Vy)  + Y* (Vy^2-1) + Z* (Vy*Vz)   =  Px * (Vx*Vy)   + Py * (Vy^2-1)
 X*(Vx*Vz)  + Y* (Vy*Vz)  + Z* (Vz^2-1)  =  Px * (Vx*Vz)   + Py* (Vy*Vz)


 */

//   | sum(Wi*Cxi),  sum(Wi*Cxyi), sum(Wi*Cxzi) |
//M= | sum(Wi*Cxyi), sum(Wi*Cyi ), sum(Wi*Cyzi) |
//   | sum(Wi*Cxzi), sum(Wi*Cyzi), sum(Wi*Czi ) |

//   | sum(Wi*(P0xi*Cxi + P0yi*Cxyi + P0zi*Cxzi)) |
//B= | sum(Wi*(P0yi*Cyi + P0xi*Cyxi + P0zi*Cyzi)) |
//   | sum(Wi*(P0zi*Czi + P0yi*Czyi + P0xi*Czxi)) |
/*
	X*(Vxi^2-1) + Y*(Vxi*Vyi) + Z*(Vxi*Vzi) = P0xi*(Vxi^2-1) +P0yi*(Vxi*Vyi) + P0zi*(Vxi*Vzi)
	X*(Vxi*Vyi) + Y*(Vyi^2-1) + Z*(Vyi*Vzi) = P0xi*(Vxi*Vyi) +P0yi*(Vyi^2-1) + P0zi*(Vyi*Vzi)
	X*(Vxi*Vzi) + Y*(Vxi*Vyi) + Z*(Vzi^2-1) = P0xi*(Vxi*Vzi) +P0yi*(Vxi*Vyi) + P0zi*(Vzi^2-1)

	X*Cx  + Y*Cxy + Z*Cxz = P0xi*Cx  +P0yi*Cxy + P0zi*Cxz
	X*Cxy + Y*Cy  + Z*Cyz = P0xi*Cxy +P0yi*Cy  + P0zi*Cyz
	X*Cxz + Y*Cyz + Z*Cz  = P0xi*Cxz +P0yi*Cyz + P0zi*Cz
	P0zi==0.0, so - now we'll use P0zi - difference from this station to average

	X*Cx  + Y*Cxy + Z*Cxz = P0xi*Cx  +P0yi*Cxy
	X*Cxy + Y*Cy  + Z*Cyz = P0xi*Cxy +P0yi*Cy
	X*Cxz + Y*Cyz + Z*Cz  = P0xi*Cxz +P0yi*Cyz

*/
						Cx=V2[0]-1.0;
						Cy=V2[1]-1.0;
						Cz=V2[2]-1.0;
						Cxy= V[0]*V[1];
						Cxz= V[0]*V[2];
						Cyz= V[1]*V[2];
						if (thisDebug) System.out.println(" Cx="+IJ.d2s(Cx,6)+" Cy="+IJ.d2s(Cy,6)+" Cz="+IJ.d2s(Cz,6)+
								" Cxy="+IJ.d2s(Cxy,6)+" Cxz="+IJ.d2s(Cxz,6)+" Cyz="+IJ.d2s(Cyz,6));


						double W=gridCorr2d[imgNum][2][index];
						double Px=gridCorr2d[imgNum][0][index];
						double Py=gridCorr2d[imgNum][1][index];
						alpha+=W;
						if (thisDebug) System.out.println(imgNum+": Px="+IJ.d2s(Px,6)+" Py="+IJ.d2s(Py,6)+" W="+IJ.d2s(W,6));
						aM[0][0]+=W*Cx;
						aM[0][1]+=W*Cxy;
						aM[0][2]+=W*Cxz;
						aM[1][1]+=W*Cy;
						aM[1][2]+=W*Cyz;
						aM[2][2]+=W*Cz;
						aB[0][0]+=W*(Px*Cx  + Py*Cxy);// Pz==0.0
						aB[1][0]+=W*(Px*Cxy + Py*Cy);// Pz==0.0
						aB[2][0]+=W*(Px*Cxz + Py*Cyz);// Pz==0.0
					}
				aM[1][0]+=aM[0][1];
				aM[2][0]+=aM[0][2];
				aM[2][1]+=aM[1][2];
				Matrix M=new Matrix(aM);
				Matrix B=new Matrix(aB);
				if (thisDebug) {
					System.out.println(" M:");
					M.print(8, 6);
					System.out.println(" B:");
					B.print(8, 6);
				}

				//			boolean fallBack2D=true;
				if ((new LUDecomposition(M)).isNonsingular()){
					double [] dXYZ=M.solve(B).getRowPackedCopy();
					for (int i=0;i<3;i++) gridCorr3d[i][index]=dXYZ[i];
					gridCorr3d[3][index]=alpha;
					fallBack2D=false; //TODO:  make sure delta Z (Math.abs(gridCorr3d[2][index])) is not too big!!
					if (Math.abs(gridCorr3d[2][index])>maxZCorr) {

						fallBack2D=true; // temporary limit
					}
					if (thisDebug) System.out.println(" dX="+IJ.d2s(gridCorr3d[0][index],6)+" dY="+IJ.d2s(gridCorr3d[1][index],6)+" dZ="+IJ.d2s(gridCorr3d[2][index],6));
				}
			}
			if(fallBack2D) { // make a 2d averaging of weighted dx, dy correction
				for (int i=0;i<gridCorr3d.length;i++) gridCorr3d[i][index]=0.0;
				double s=0;
				for (int imgNum=0;imgNum<selectedImages.length;imgNum++) if ((gridCorr2d[imgNum]!=null) &&(gridCorr2d[imgNum][2][index]>0.0)) {
					double W=gridCorr2d[imgNum][2][index];
					//			V[i]=patternParameters.gridGeometry[v][u][i]-cameraXYZ[imgNum][i];

					double [] cv={
							patternParameters.gridGeometry[v][u][0]-cameraXYZ[imgNum][0],
							patternParameters.gridGeometry[v][u][1]-cameraXYZ[imgNum][1],
							patternParameters.gridGeometry[v][u][2]-cameraXYZ[imgNum][2]};
					double cv2=cv[0]*cv[0]+cv[1]*cv[1]+cv[2]*cv[2];
					double kv=(gridCorr2d[imgNum][0][index]*cv[0]+gridCorr2d[imgNum][1][index]*cv[1])/cv2;
					gridCorr3d[0][index]+=W*(gridCorr2d[imgNum][0][index]-cv[0]*kv);
					gridCorr3d[1][index]+=W*(gridCorr2d[imgNum][1][index]-cv[1]*kv);
					gridCorr3d[2][index]+=W*(                            -cv[2]*kv);
					s+=W;
				}
				if (s>0){
					gridCorr3d[0][index]/=s;
					gridCorr3d[1][index]/=s;
					gridCorr3d[2][index]/=s;
				} else {
					gridCorr3d[0][index]=0.0;
					gridCorr3d[1][index]=0.0;
					gridCorr3d[2][index]=0.0;
				}
				gridCorr3d[3][index]=s;
				if (thisDebug) System.out.println(" Using 2d averaging: dX="+IJ.d2s(gridCorr3d[0][index],6)+
						" dY="+IJ.d2s(gridCorr3d[1][index],6)+" dZ="+IJ.d2s(gridCorr3d[2][index],6));
			}
		}
		// Make average correction zero is it needed?
		// create "reliable" mask for averaging/tilting - disregard the outmost grid pixels
		boolean [] reliable=new boolean [width*height];
		double wThreshold=0.0;
		for (int v=0;v<height;v++) for (int u=0;u<width;u++){
			int index=u+v*width;
			reliable[index]=false;
			if ((v>0) && (u>0) && (v<(height-1)) && (u<(width-1)) &&
					(gridCorr3d[3][index]>wThreshold) &&
					(gridCorr3d[3][index-1]>wThreshold) &&
					(gridCorr3d[3][index+1]>wThreshold) &&
					(gridCorr3d[3][index-width]>wThreshold) &&
					(gridCorr3d[3][index+width]>wThreshold) ){
				reliable[index]=true;
			}

		}
		double corrAverage;
		for (int c=0;c<3;c++){
			corrAverage=0.0;
			double s=0.0;
			for (int i=0;i<gridCorr3d[0].length;i++) if (reliable[i]) {
				corrAverage+=gridCorr3d[c][i]*gridCorr3d[3][i];
				s+=gridCorr3d[3][i];
			}
			corrAverage/=s;
//			System.out.println("zCorrAverage["+c+"="+corrAverage);
			for (int i=0;i<gridCorr3d[c].length;i++) gridCorr3d[c][i]-=corrAverage;
		}
// for Z correction compensate for x/y tilts
		String [] titles={"X-correction(mm)","Y-correction(mm)","Z-correction","Weight"};
		if (rotateCorrection) {
			double SX=0.0,SX2=0.0,SZ=0.0,SXY=0.0,SXZ=0.0,S0=0.0,SY=0.0,SY2=0.0,SYZ=0.0;
			double [][] gridGeom=new double [3][gridCorr3d[0].length];
			for (int c=0;c<gridGeom.length;c++) for (int i=0;i<gridGeom[c].length;i++)gridGeom[c][i]=0.0;

			for (int v=0;v<height;v++) for (int u=0;u<width; u++){
				int index=u+v*width;
				double W=gridCorr3d[3][index];
				gridGeom[0][index]=patternParameters.gridGeometry[v][u][0];
				gridGeom[1][index]=patternParameters.gridGeometry[v][u][1];
				gridGeom[2][index]=W;
				if ((reliable[index]) && (W>0.0)){
					S0+=W;
					SX+=  W*patternParameters.gridGeometry[v][u][0];
					SX2+= W*patternParameters.gridGeometry[v][u][0]*patternParameters.gridGeometry[v][u][0];
					SY+=  W*patternParameters.gridGeometry[v][u][1];
					SY2+= W*patternParameters.gridGeometry[v][u][1]*patternParameters.gridGeometry[v][u][1];
					SXY+= W*patternParameters.gridGeometry[v][u][0]*patternParameters.gridGeometry[v][u][1];
					SZ+=  W*gridCorr3d[2][index];
					SXZ+= W*gridCorr3d[2][index]*patternParameters.gridGeometry[v][u][0];
					SYZ+= W*gridCorr3d[2][index]*patternParameters.gridGeometry[v][u][1];
				}
			}
			double [][] aM= {
					{SX2, SXY, SX},
					{SXY, SY2, SY},
					{SX,  SY,  S0}};
			double [][] aB= {{SXZ},{SYZ},{SZ}};
			Matrix M=new Matrix(aM);
			Matrix B=new Matrix(aB);
			if (this.debugLevel>2) {
				System.out.println(" M:");
				M.print(8, 6);
				System.out.println(" B:");
				B.print(8, 6);
				System.out.println(" Ax,Ay,B:");
				M.solve(B).print(8, 6);
			}
			double [] tilts=M.solve(B).getRowPackedCopy();
			if (this.debugLevel>2) {
				if (this.refineParameters.showThisCorrection) {
					ShowDoubleFloatArrays.showArrays(gridCorr3d, getGridWidth(), getGridHeight(),  true, "before tilt:", titles);
				}
			}
			for (int v=0;v<height;v++) for (int u=0;u<width; u++){
				int index=u+v*width;
				gridCorr3d[2][index]-=tilts[0]*patternParameters.gridGeometry[v][u][0]+tilts[1]*patternParameters.gridGeometry[v][u][1]+tilts[2];
			}
		}
    	if (this.debugLevel>2) {
    		if (this.refineParameters.showThisCorrection) {
    			double [][] gridCorr3dClone=new double [4][width*height];
    			for (int c=0;c<gridCorr3dClone.length;c++) for (int i=0;i<gridCorr3dClone[c].length;i++)
    				gridCorr3dClone[c][i]=reliable[i]? gridCorr3d[c][i]:0.0;
    			ShowDoubleFloatArrays.showArrays(gridCorr3dClone, getGridWidth(), getGridHeight(),  true, "after tilt:", titles);
    		}
    	}
    	IJ.showStatus("");
		return  gridCorr3d;
	}
/**
 *
 * @param gridCorr3D Array of grid corrections (1-st index: dx, dy, dz, mask (>0 - valid point)
 * @param gridZCorr Optional per-station z-correction (or null)
 * @param width // width of the grid array
 * @param preShrink // shrink input array by this number of pixels (hor/vert) befere extrapolating (remove bad border nodes)
 * @param expand    // expand/extrapolate this number of steps after shrinking (or until no pixels left
 * @param sigma     // effective radius for fitting the extrapolation plane, in nodes
 * @param ksigma    // size if square to consider (measured in ksigma-s). 2.0 means square is 4*sigma by 4*sigma
 * @return true if OK, false if error
 */
	public boolean shrinkExtrapolateGridCorrection(
			double [][] gridCorr3D, // dx,dy,dz, mask >0
			double [][] gridZCorr, // per-station additional Z-correction (or null)
			int width,
			int preShrink,
			int expand,
			double sigma,
			double ksigma){
		int length=gridCorr3D[0].length;
        int height=	length/width;
//		int decimate=fittingStrategy.distortionCalibrationData.eyesisCameraParameters.decimateMasks;
//		int sWidth= (fittingStrategy.distortionCalibrationData.eyesisCameraParameters.sensorWidth-1)/decimate+1;
//		int sHeight=(fittingStrategy.distortionCalibrationData.eyesisCameraParameters.sensorHeight-1)/decimate+1;
//		double sigma=nsigma/decimate;
		boolean [] fMask=new boolean[length];
		for (int i=0;i<fMask.length;i++)
			fMask[i]=gridCorr3D[3][i]>0;
		int len= (int) Math.ceil(sigma*ksigma);
		double [] gaussian=new double[len+1];
		double k=0.5/sigma/sigma;
		for (int i=0;i<=len;i++) gaussian[i]=Math.exp(-i*i*k);
		int [][] dirs={{-1,0},{1,0},{0,-1},{0,1}};
		List <Integer> extList=new ArrayList<Integer>(1000);
		Integer Index,Index2;
		extList.clear();
		// create initial wave
		int debugThreshold=2;
		if (this.debugLevel>debugThreshold) System.out.println("shrinkExtrapolateGridCorrection width="+width+" height="+height);

		for (int iy=0;iy<height;iy++) for (int ix=0;ix<width;ix++) {
			Index=iy*width+ix;
			if (fMask[Index]) {
				int numNew=0;
				for (int dir=0;dir<dirs.length;dir++){
					int ix1=ix+dirs[dir][0];
					int iy1=iy+dirs[dir][1];
					// Will not shrink from the array border!
					if ((ix1>=0) && (iy1>=0) && (ix1<width) && (iy1<height)) {
						if (!fMask[iy1*width+ix1]) numNew++;
					}
					if (numNew>0) extList.add(Index); // neighbor will have non-singular matrix
				}
			}
		}
		if (this.debugLevel>debugThreshold) System.out.println("Initial wave length="+extList.size());
		// now shrink
		// unmask current wave
		for (int i=extList.size()-1; i>=0;i--) fMask[extList.get(i)]=false;
		if (extList.size()==0) return false; // no points
		for (int nShrink=0;nShrink<preShrink;nShrink++){
			int size=extList.size();
			if (this.debugLevel>debugThreshold) System.out.println("shrinking, size="+size);
			if (size==0) return false; // no points
			// wave step, unmasking
			for (int i=0; i<size;i++) {
				Index=extList.get(0);
				extList.remove(0);
				int iy=Index/width;
				int ix=Index%width;
				for (int dir=0;dir<dirs.length;dir++){
					int ix1=ix+dirs[dir][0];
					int iy1=iy+dirs[dir][1];
					if ((ix1>=0) && (iy1>=0) && (ix1<width) && (iy1<height)){
						Index=iy1*width+ix1;
						if (fMask[Index]){
							extList.add(Index);
							fMask[Index]=false; // restore later?
						}
					}
				}
			}
		}
		// restore mask on the front
		for (int i=extList.size()-1; i>=0;i--) fMask[extList.get(i)]=true;

		   // repeat with the wave until there is place to move, but not more than "expand" steps
		   int [] dirs2=new int [2];
		   for (int n=0; (n<expand) && (extList.size()>0); n++ ){
			   if (this.updateStatus) IJ.showStatus("Expanding, step="+(n+1)+" (of "+expand+"), extList.size()="+extList.size());
//			   if (this.updateStatus) showStatus("Expanding, step="+(n+1)+" (of "+expand+"), extList.size()="+extList.size(),0);
			   if (this.debugLevel>debugThreshold) System.out.println("Expanding, step="+n+", extList.size()="+extList.size());
			   // move wave front 1 pixel hor/vert
			   for (int i=extList.size();i>0;i--){ // repeat current size times
				   Index=extList.get(0);
				   extList.remove(0);
				   int iy=Index/width;
				   int ix=Index%width;
				   for (int dir=0;dir<dirs.length;dir++){
					   int ix1=ix+dirs[dir][0];
					   int iy1=iy+dirs[dir][1];
					   if ((ix1>=0) && (iy1>=0) && (ix1<width) && (iy1<height)){
						   Index=iy1*width+ix1;
						   if (!fMask[Index]){
							   // verify it has neighbors in the perpendicular direction to dir
							   dirs2[0]=(dir+2) & 3;
							   dirs2[1]=dirs2[0] ^ 1;
							   for (int dir2=0;dir2<dirs2.length;dir2++){
								   int ix2=ix+dirs[dirs2[dir2]][0]; // from the old, not the new point!
								   int iy2=iy+dirs[dirs2[dir2]][1];
								   if ((ix2>=0) && (iy2>=0) && (ix2<width) && (iy2<height)){
									   Index2=iy2*width+ix2;
									   if (fMask[Index2]){ // has orthogonal neighbor, OK to add
										   extList.add(Index);
										   fMask[Index]=true; // remove later
										   break;
									   }
								   }
							   }
						   }
					   }
				   }
			   }
			   // now un-mask the pixels in new list new
			   for (int i =0;i<extList.size();i++){
				   Index=extList.get(i);
				   fMask[Index]=false; // now mask is only set for known pixels
			   }
	// Calculate values (extrapolate) for the pixels in the list
			/*
Err = sum (W(x,y)*(f(x,y)-F0-Ax*(x-X0)-Ay*(y-Y0))^2)=
sum (Wxy*(Fxy^2+F0^2+Ax^2*(x-X0)^2+Ay^2*(y-Y0)^2
-2*Fxy*F0 -2*Fxy*Ax*(x-X0) - 2*Fxy*Ay*(y-Y0)
+2*F0*Ax*(x-X0) + 2*F0*Ay*(y-Y0)
+2*Ax*(x-X0)*Ay*(y-Y0))
(1)0=dErr/dF0= 2*sum (Wxy*(F0-Fxy+Ax*(x-X0)+Ay(y-Y0)))
(2)0=dErr/dAx= 2*sum (Wxy*(Ax*(x-X0)^2-Fxy*(x-X0) +F0*(x-X0)+Ay*(x-x0)*(y-Y0)))
(3)0=dErr/dAy= 2*sum (Wxy*(Ay*(y-y0)^2-Fxy*(y-Y0) +F0*(y-Y0)+Ax*(x-x0)*(y-Y0)))

S0 = sum(Wxy)
SF=  sum(Wxy*Fxy)
SX=  sum(Wxy*(x-X0)
SY=  sum(Wxy*(y-Y0)
SFX= sum(Wxy*Fxy*(x-X0)
SFY= sum(Wxy*Fxy*(y-Y0)
SX2= sum(Wxy*(x-X0)^2
SY2= sum(Wxy*(y-Y0)^2
SXY= sum(Wxy*(x-X0)*(y-Y0)

(1) F0*S0 - SF + Ax*SX +Ay*Sy = 0
(2) Ax*SX2-SFX+F0*SX+Ay*SXY = 0
(3) Ay*Sy2 -SFY + F0*SY +Ax*SXY = 0

(1) F0*S0  + Ax*SX +Ay*SY = SF
(2) Ax*SX2+F0*SX+Ay*SXY = SFX
(3) Ay*Sy2  + F0*SY +Ax*SXY = SFY


   | F0 |
V= | Ax |
   | Ay |

     | SF  |
B =  | SFX |
     | SFY |

     | S0  SX   SY  |
M =  | SX  SX2  SXY |
     | SY  SXY  SY2 |

M * V = B
			 */
			   int numOriginalComponents=3;
			   boolean useExtra=gridZCorr!=null;
			for (int i =0;i<extList.size();i++){
        		Index=extList.get(i);
        		int iy=Index/width;
        		int ix=Index%width;
        		double [] S0=new double [3+(useExtra?gridZCorr.length:0)];
        		for (int ii=0;ii<S0.length;ii++) S0[ii]=0.0;
//        		double [] S0= {0.0,0.0,0.0};
        		double [] SF= S0.clone();
        		double [] SX= S0.clone();
        		double [] SY= S0.clone();
        		double [] SFX=S0.clone();
        		double [] SFY=S0.clone();
        		double [] SX2=S0.clone();
        		double [] SY2=S0.clone();
        		double [] SXY=S0.clone();
        		int iYmin=iy-len; if (iYmin<0) iYmin=0;
        		int iYmax=iy+len; if (iYmax>=height) iYmax=height-1;
        		int iXmin=ix-len; if (iXmin<0) iXmin=0;
        		int iXmax=ix+len; if (iXmax>=width) iXmax=width-1;
        		for (int iy1=iYmin;iy1<=iYmax;iy1++) for (int ix1=iXmin;ix1<=iXmax;ix1++) {
        			int ind=ix1+iy1*width;
        			if (fMask[ind]){
        				double w=gaussian[(iy1>=iy)?(iy1-iy):(iy-iy1)]*gaussian[(ix1>=ix)?(ix1-ix):(ix-ix1)];
        				for (int m=0;m<S0.length;m++){
        					double d=(m<numOriginalComponents)?gridCorr3D[m][ind]:gridZCorr[m-numOriginalComponents][ind];
        					S0[m]+= w;
        					SF[m]+= w*d;
        					SX[m]+= w*(ix1-ix);
        					SY[m]+= w*(iy1-iy);
        					SFX[m]+=w*d*(ix1-ix);
        					SFY[m]+=w*d*(iy1-iy);
        					SX2[m]+=w*(ix1-ix)*(ix1-ix);
        					SY2[m]+=w*(iy1-iy)*(iy1-iy);
        					SXY[m]+=w*(ix1-ix)*(iy1-iy);
        				}
        			}

        		}
        		for (int m=0;m<S0.length;m++){
        			double [][] aB={{SF[m]},{SFX[m]},{SFY[m]}};
        			double [][] aM={
        					{S0[m],SX[m], SY[m]},
        					{SX[m],SX2[m],SXY[m]},
        					{SY[m],SXY[m],SY2[m]}
        			};
        			Matrix B=new Matrix(aB);
        			Matrix M=new Matrix(aM);
        			Matrix V=M.solve(B);
        			if (m<numOriginalComponents) gridCorr3D[m][Index]=V.get(0,0);
        			else gridZCorr[m-numOriginalComponents][Index]=V.get(0,0);
        		}
    			if (this.debugLevel>debugThreshold) System.out.println("updated v="+(Index/width)+" u="+(Index%width)+" {"+
    					IJ.d2s(gridCorr3D[0][Index],2)+","+IJ.d2s(gridCorr3D[1][Index],2)+","+IJ.d2s(gridCorr3D[2][Index],2)+"}");
			}

// set mask again for the new calculated layer of pixels
			for (int i =0;i<extList.size();i++){
        		Index=extList.get(i);
				fMask[Index]=true;
			}
        }
	   return true;
	}

	public void logScale(
			double [] data,
			double fatZero){
		for (int i=0;i<data.length;i++){
			double d=((data[i]>=0)?data[i]:0.0);
			data[i]=(fatZero>0)?(Math.log(fatZero+d)):d;
		}
	}
	public void unLogScale(
			double [] data,
			double fatZero){
		for (int i=0;i<data.length;i++){
			if (fatZero>0.0) data[i]=Math.exp(data[i])-fatZero;
			if (data[i]<0.0) data[i]=0.0;
		}
	}

	/**
	 * Extrapolates sensor correction beyond known data (in-place)
	 * @param fieldXY [2][nPixels] vector field to extrapolate
	 * @param sMask [nPixels] alpha (0.0 .. 1.0) "reliability" mask to apply to vector field
	 * @param alphaThreshold start with pixels with alpha above this value (disregard border unreliable pixels)
	 * @param nsigma when fitting plane through new point use Gaussian weight function for the neighbors
	 *  (normalized to non-decimated points)
	 * @param ksigma Process pixels in a square with the side 2*sigma*ksigma
	 * @return false if nothing to extrapolate (too small mask)?
	 */
	public boolean extrapolateSensorCorrection(
			int    numChn,
			boolean [] whichExtrapolate,
			double [][] fieldXY,
			double []sMask,
			double alphaThreshold,
			double nsigma,
			double ksigma){

		int decimate = getDecimateMasks(numChn);
		int sWidth =   (getSensorWidth(numChn)-1)/decimate+1;
		int sHeight =  (getSensorHeight(numChn)-1)/decimate+1;
		double sigma=nsigma/decimate;
		boolean [] fMask=new boolean[fieldXY[0].length];
		for (int i=0;i<fMask.length;i++)
			fMask[i]=sMask[i]>=alphaThreshold;
		int len= (int) Math.ceil(sigma*ksigma);
		double [] gaussian=new double[len+1];
		double k=0.5/sigma/sigma;
		for (int i=0;i<=len;i++) gaussian[i]=Math.exp(-i*i*k);
		int [][] dirs={{-1,0},{1,0},{0,-1},{0,1}};
		List <Integer> extList=new ArrayList<Integer>(1000);
		Integer Index;
		extList.clear();
		// create initial wave
		if (this.debugLevel>2) System.out.println("extrapolateSensorCorrection() decimate="+decimate+", sWidth="+sWidth+" sHeight="+sHeight);
		for (int iy=0;iy<sHeight;iy++) for (int ix=0;ix<sWidth;ix++) {
			Index=iy*sWidth+ix;
			if (fMask[Index]) {
				int numOld=0;
				int numNew=0;
				for (int dir=0;dir<dirs.length;dir++){
					int ix1=ix+dirs[dir][0];
					int iy1=iy+dirs[dir][1];
					if ((ix1>=0) && (iy1>=0) && (ix1<sWidth) && (iy1<sHeight)) {
						if (fMask[iy1*sWidth+ix1]) numOld++;
						else numNew++;
					}
					if ((numNew>0) && (numOld>1)) extList.add(Index); // neighbor will have non-singular matrix
				}
			}
		}
		if (extList.size()==0) return false;
        while (extList.size()>0){
    		if (this.debugLevel>2) System.out.println("extList.size()="+extList.size());

        	// move wave front 1 pixel hor/vert
        	for (int i=extList.size();i>0;i--){ // repeat current size times
        		Index=extList.get(0);
        		extList.remove(0);
        		int iy=Index/sWidth;
        		int ix=Index%sWidth;
				for (int dir=0;dir<dirs.length;dir++){
					int ix1=ix+dirs[dir][0];
					int iy1=iy+dirs[dir][1];
					if ((ix1>=0) && (iy1>=0) && (ix1<sWidth) && (iy1<sHeight)){
						Index=iy1*sWidth+ix1;
						if (!fMask[Index]){
							extList.add(Index);
							fMask[Index]=true; // remove later
						}
					}
				}
        	}
			// now un-mask the pixels in new list new
			for (int i =0;i<extList.size();i++){
        		Index=extList.get(i);
				fMask[Index]=false; // now mask is only set for known pixels
			}
// Calculate values (extrapolate) for the pixels in the list
			/*
Err = sum (W(x,y)*(f(x,y)-F0-Ax*(x-X0)-Ay*(y-Y0))^2)=
sum (Wxy*(Fxy^2+F0^2+Ax^2*(x-X0)^2+Ay^2*(y-Y0)^2
-2*Fxy*F0 -2*Fxy*Ax*(x-X0) - 2*Fxy*Ay*(y-Y0)
+2*F0*Ax*(x-X0) + 2*F0*Ay*(y-Y0)
+2*Ax*(x-X0)*Ay*(y-Y0))
(1)0=dErr/dF0= 2*sum (Wxy*(F0-Fxy+Ax*(x-X0)+Ay(y-Y0)))
(2)0=dErr/dAx= 2*sum (Wxy*(Ax*(x-X0)^2-Fxy*(x-X0) +F0*(x-X0)+Ay*(x-x0)*(y-Y0)))
(3)0=dErr/dAy= 2*sum (Wxy*(Ay*(y-y0)^2-Fxy*(y-Y0) +F0*(y-Y0)+Ax*(x-x0)*(y-Y0)))

S0 = sum(Wxy)
SF=  sum(Wxy*Fxy)
SX=  sum(Wxy*(x-X0)
SY=  sum(Wxy*(y-Y0)
SFX= sum(Wxy*Fxy*(x-X0)
SFY= sum(Wxy*Fxy*(y-Y0)
SX2= sum(Wxy*(x-X0)^2
SY2= sum(Wxy*(y-Y0)^2
SXY= sum(Wxy*(x-X0)*(y-Y0)

(1) F0*S0 - SF + Ax*SX +Ay*Sy = 0
(2) Ax*SX2-SFX+F0*SX+Ay*SXY = 0
(3) Ay*Sy2 -SFY + F0*SY +Ax*SXY = 0

(1) F0*S0  + Ax*SX +Ay*SY = SF
(2) Ax*SX2+F0*SX+Ay*SXY = SFX
(3) Ay*Sy2  + F0*SY +Ax*SXY = SFY


   | F0 |
V= | Ax |
   | Ay |

     | SF  |
B =  | SFX |
     | SFY |

     | S0  SX   SY  |
M =  | SX  SX2  SXY |
     | SY  SXY  SY2 |

M * V = B
			 */
			double [] zeros= new double[whichExtrapolate.length];
			for (int i=0;i<zeros.length;i++)zeros[i]=0.0;
			for (int i =0;i<extList.size();i++){
        		Index=extList.get(i);
        		int iy=Index/sWidth;
        		int ix=Index%sWidth;
        		double [] S0= zeros.clone();
        		double [] SF= zeros.clone();
        		double [] SX= zeros.clone();
        		double [] SY= zeros.clone();
        		double [] SFX=zeros.clone();
        		double [] SFY=zeros.clone();
        		double [] SX2=zeros.clone();
        		double [] SY2=zeros.clone();
        		double [] SXY=zeros.clone();
        		int iYmin=iy-len; if (iYmin<0) iYmin=0;
        		int iYmax=iy+len; if (iYmax>=sHeight) iYmax=sHeight-1;
        		int iXmin=ix-len; if (iXmin<0) iXmin=0;
        		int iXmax=ix+len; if (iXmax>=sWidth) iXmax=sWidth-1;
        		for (int iy1=iYmin;iy1<=iYmax;iy1++) for (int ix1=iXmin;ix1<=iXmax;ix1++) {
        			int ind=ix1+iy1*sWidth;
        			if (fMask[ind]){
        				double w=gaussian[(iy1>=iy)?(iy1-iy):(iy-iy1)]*gaussian[(ix1>=ix)?(ix1-ix):(ix-ix1)];
        				for (int m=0;m<whichExtrapolate.length;m++) if(whichExtrapolate[m]){
        					S0[m]+= w;
        					SF[m]+= w*fieldXY[m][ind];
        					SX[m]+= w*(ix1-ix);
        					SY[m]+= w*(iy1-iy);
        					SFX[m]+=w*fieldXY[m][ind]*(ix1-ix);
        					SFY[m]+=w*fieldXY[m][ind]*(iy1-iy);
        					SX2[m]+=w*(ix1-ix)*(ix1-ix);
        					SY2[m]+=w*(iy1-iy)*(iy1-iy);
        					SXY[m]+=w*(ix1-ix)*(iy1-iy);
        				}
        			}

        		}
        		for (int m=0;m<whichExtrapolate.length;m++) if(whichExtrapolate[m]){
        			double [][] aB={{SF[m]},{SFX[m]},{SFY[m]}};
        			double [][] aM={
        					{S0[m],SX[m], SY[m]},
        					{SX[m],SX2[m],SXY[m]},
        					{SY[m],SXY[m],SY2[m]}
        			};
        			Matrix B=new Matrix(aB);
        			Matrix M=new Matrix(aM);
        			Matrix V=M.solve(B);
        			fieldXY[m][Index]=V.get(0,0);
        		}

			}

// set mask again for the new calculated layer of pixels
			for (int i =0;i<extList.size();i++){
        		Index=extList.get(i);
				fMask[Index]=true;
			}
        }
		return true;

	}


	/**
	 * Calculates residual correction from  the measured sensor pX, pY to the calculated {pixel X, pixel Y}
	 * @param distortionCalibrationData
	 * @param showIndividual - show individual images
	 * @param showIndividualNumber - which image to show (-1 - all)
	 * @return
	 */
	public double [][][] calculateSensorXYCorr(
			DistortionCalibrationData distortionCalibrationData,
			boolean showIndividual,
			int showIndividualNumber, // which image to show (-1 - all)
			boolean useGridAlpha // use grid alpha, false - use old calculations
			){
		int numChannels=distortionCalibrationData.getNumChannels(); // number of used channels
		int width=getGridWidth();
		int height=getGridHeight();
    	int imgRGBIndex=   3;

		int [] uvInc={0,1,width,width+1}; // four corners as vu index
		int [][] cycles={ // counter-clockwise corners bounding the area  (only orthogonal sides?)
				{1,0,2},
				{2,3,1},
				{0,2,3},
				{3,1,0}};
		double [][][] gridPCorr=new double [numChannels][][];
		for (int chnNum=0;chnNum<gridPCorr.length;chnNum++) gridPCorr[chnNum]=null;
		boolean [] selectedImages=fittingStrategy.selectedImages();
		boolean debugExit=false;
		int debugCntr=2;
		int numSelected=0;
		for (int imgNum=0;imgNum<selectedImages.length;imgNum++) if (selectedImages[imgNum]) numSelected++;
		int numProcessed=0;
		IJ.showStatus("Calculating sensor corrections...");
		for (int imgNum=0;imgNum<selectedImages.length;imgNum++) if (selectedImages[imgNum]) {
			if (debugExit) break;
			int chnNum=fittingStrategy.distortionCalibrationData.gIP[imgNum].channel; // number of sub-camera

			int decimate=getDecimateMasks(chnNum);
			int sWidth= (getSensorWidth(chnNum)-1)/decimate+1;
			int sHeight=(getSensorHeight(chnNum)-1)/decimate+1;

			int station=fittingStrategy.distortionCalibrationData.gIP[imgNum].getStationNumber(); // number of sub-camera
			double [][] photometrics=patternParameters.getPhotometricBySensor(station,chnNum); // head/bottom grid intensity/alpha

			if (showIndividual && ((showIndividualNumber<0) || (showIndividualNumber==chnNum))) {
				String [] titles={"R","G","B","A"};
				ShowDoubleFloatArrays.showArrays(photometrics, width, height,  true, "Photometrics"+chnNum+"-"+imgNum, titles);
			}


			// initialize this array if it is needed, leave unused null
			if (gridPCorr[chnNum]==null){
				 gridPCorr[chnNum]=new double [7][sWidth*sHeight];
				for (int n=0;n<gridPCorr[chnNum].length;n++) for (int i=0;i<gridPCorr[chnNum][0].length;i++) gridPCorr[chnNum][n][i]=0.0;
			}
			double [][] thisPCorr=null;

			thisPCorr=new double [7][sWidth*sHeight]; // calculate for a single (this) image, accumulate in the end
			for (int n=0;n<thisPCorr.length;n++) for (int i=0;i<thisPCorr[0].length;i++) thisPCorr[n][i]=0.0;
			double [] diff=calcYminusFx(this.currentfX);
			// find data range for the selected image
			int index=0;
			int numImg=fittingStrategy.distortionCalibrationData.getNumImages();
			for (int iNum=0;(iNum<imgNum) && (iNum<numImg) ;iNum++) if (selectedImages[iNum])
				index+=fittingStrategy.distortionCalibrationData.gIP[iNum].pixelsUV.length;
			if (this.debugLevel>2) {
				System.out.println("calculateGridXYCorr(): fX.length="+this.currentfX.length+" this image index="+index);
			}
			double [][] imgData=new double[8][height * width]; // dPX, dPY, Px, Py, alpha,R,G,B
			for (int i=0;i<imgData.length;i++) for (int j=0;j<imgData[i].length;j++)imgData[i][j]=0.0;

			for (int i=0;i<fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV.length;i++){
				int u=fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV[i][0]+patternParameters.U0; // starting from 0
				int v=fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV[i][1]+patternParameters.V0; // starting from 0
				int vu=u+width*v;
				imgData[0][vu]=   diff[2*(index+i)];
				imgData[1][vu]=   diff[2*(index+i)+1];
				imgData[2][vu]= this.Y[2*(index+i)];  // measured pixel x
				imgData[3][vu]= this.Y[2*(index+i)+1];// measured pixel y
				imgData[4][vu]= this.weightFunction[2*(index+i)];

				for (int c=0;c<3;c++){
//					double g=gridGeometry[v][u][gridRGBIndex+c];
					double g=photometrics[c][vu];
					imgData[5+c][vu]=(g>0)?(fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsXY[i][imgRGBIndex+c]/g):g;
				}
			}
			if (showIndividual && ((showIndividualNumber<0) || (showIndividualNumber==chnNum))) {
				String [] titles={"dPx","dPy","Px","Py","A","R","G","B"};// dPX, dPY, Px, Py, alpha,R,G,B - rgb - full, not incremental
				ShowDoubleFloatArrays.showArrays(imgData, width, height,  true, "imgData"+imgNum, titles);

			}

			// now use imgData array to fill thisPCorr by linear interpolation
			for (int v=0;v<(height-1); v++) for (int u=0; u<(width-1);u++){
				if (debugExit) break;
				int vu=u+width*v;
                double [][] cornerXY =new double[4][];
                for (int i=0;i<uvInc.length;i++){
                	int vu1=vu+uvInc[i];
                	if (imgData[4][vu1]>0.0){
                		cornerXY[i]=new double[2];
                		cornerXY[i][0]=imgData[2][vu1];
                		cornerXY[i][1]=imgData[3][vu1];
                	} else cornerXY[i]=null;
                }
                boolean [] cycleFits=new boolean[cycles.length];
                boolean anyFits=false;
                for (int i=0;i<cycles.length;i++){
                	cycleFits[i]=true;
                	for (int j=0;j<cycles[i].length;j++) if (cornerXY[cycles[i][j]]==null) {
                		cycleFits[i]=false;
                		break;
                	}
                	anyFits |=cycleFits[i];
                }
                if (!anyFits) continue; // not a single cycle
				if ((this.debugLevel>3) && !debugExit) {
					String debugString="cycleFits ";
					for (int i =0;i<cycleFits.length; i++) debugString+=" "+cycleFits[i];
					System.out.println(debugString);
				}
                if (cycleFits[0]&&cycleFits[1]){ // remove overlaps
                	cycleFits[2]=false;
                	cycleFits[3]=false;
                }
                boolean minMaxUndefined=true;
				double minX=0,maxX=0,minY=0,maxY=0;
				// find bounding rectangle;
				for (int nCycle=0;nCycle<cycles.length;nCycle++) if (cycleFits[nCycle]){
					int [] cycle=cycles[nCycle];
					for (int corner=0; corner<cycle.length;corner++){
						if (minMaxUndefined || (minX>cornerXY[cycle[corner]][0])) minX=cornerXY[cycle[corner]][0];
						if (minMaxUndefined || (maxX<cornerXY[cycle[corner]][0])) maxX=cornerXY[cycle[corner]][0];
						if (minMaxUndefined || (minY>cornerXY[cycle[corner]][1])) minY=cornerXY[cycle[corner]][1];
						if (minMaxUndefined || (maxY<cornerXY[cycle[corner]][1])) maxY=cornerXY[cycle[corner]][1];
						minMaxUndefined=false;
					}
				}
				int iMinX=(int) Math.floor(minX/decimate);
				int iMinY=(int) Math.floor(minY/decimate);
				int iMaxX=(int) Math.ceil(maxX/decimate);
				int iMaxY=(int) Math.ceil(maxY/decimate);
				// not sure if these checks are needed, got out of bounds wheriDy was =484=sHeight
				if (iMinX<0) iMinX=0;
				if (iMaxX>=sWidth) iMaxX=sWidth-1;
				if (iMinY<0) iMinY=0;
				if (iMaxY>=sHeight) iMaxY=sHeight-1;
				double [] originXY=new double [2];
				double [] endXY=new double [2];
				boolean debugHadPixels=false;
//TODO: scan X,Y in this rectangle, for points in defined squares/triangles find if the point is inside (accurate not to loose any).
				for (int idY=iMinY; idY<=iMaxY;idY++){

					double pY=idY*decimate; // in sensor pixels
					for (int idX=iMinX; idX<=iMaxX;idX++){
						double pX=idX*decimate; // in sensor pixels
						// scan allowed triangles, usually 2
						for (int nCycle=0;nCycle<cycles.length;nCycle++) if (cycleFits[nCycle]){
							int [] cycle=cycles[nCycle];
							// is this point inside?
							if (debugExit) {
								for (int nEdge=0;nEdge<cycle.length;nEdge++){
									int nextNEdge=(nEdge==(cycle.length-1))?0:(nEdge+1);
									System.out.println("nEdge="+nEdge+" nextNEdge"+nextNEdge);

									originXY[0]=imgData[2][vu+uvInc[cycle[nEdge]]];
									originXY[1]=imgData[3][vu+uvInc[cycle[nEdge]]];
									endXY[0]=   imgData[2][vu+uvInc[cycle[nextNEdge]]];
									endXY[1]=   imgData[3][vu+uvInc[cycle[nextNEdge]]];
									System.out.println("--- pX="+IJ.d2s(pX,1)+" originXY[0]="+IJ.d2s(originXY[0],1)+
											" endXY[1]="+IJ.d2s(endXY[1],1)+" originXY[1]="+IJ.d2s(originXY[1],1));
									System.out.println("--- pY="+IJ.d2s(pY,1)+" originXY[1]="+IJ.d2s(originXY[1],1)+
											" endXY[0]="+IJ.d2s(endXY[0],1)+" originXY[0]="+IJ.d2s(originXY[0],1));
									System.out.println("Cross-product="+IJ.d2s(((pX-originXY[0])*(endXY[1]-originXY[1]) - (pY-originXY[1])*(endXY[0]-originXY[0])),1));

								}
							}

							boolean inside=true;
							for (int nEdge=0;nEdge<cycle.length;nEdge++){
								int nextNEdge=(nEdge==(cycle.length-1))?0:(nEdge+1);

								originXY[0]=imgData[2][vu+uvInc[cycle[nEdge]]];
								originXY[1]=imgData[3][vu+uvInc[cycle[nEdge]]];
								endXY[0]=   imgData[2][vu+uvInc[cycle[nextNEdge]]];
								endXY[1]=   imgData[3][vu+uvInc[cycle[nextNEdge]]];
								if (((pX-originXY[0])*(endXY[1]-originXY[1]) - (pY-originXY[1])*(endXY[0]-originXY[0]))<0.0){
									inside=false;
									break;
								}
							}
							if (!inside) continue; // point is outside of the interpolation area, try next triangle (if any)
//							if ((this.debugLevel>3) && !debugExit) {
							if (this.debugLevel>3) {
								System.out.println("idX="+idX+" idY="+idY+" nCycle="+nCycle);
								String debugString1="cycle:";
								for (int i =0;i<cycle.length; i++) debugString1+=" "+cycle[i];
								System.out.println(debugString1);
							}

							/* interpolate:
							1. taking cycles[0] as origin and two (non co-linear) edge vectors - V1:from 0 to 1 and V2 from 1 to 2
							    find a1 and a2  so that vector V  (from 0  to pXY) = a1*V1+ a2*V2
							2. if F0 is the value of the interpolated function at cycles[0], F1 and F2 - at cycles[1] and cycles2
							   then F=F0+(F1-F0)*a1 +(F2-F1)*a2
							 */
							double [] XY0={imgData[2][vu+uvInc[cycle[0]]],imgData[3][vu+uvInc[cycle[0]]]};
							double [] XY1={imgData[2][vu+uvInc[cycle[1]]],imgData[3][vu+uvInc[cycle[1]]]};
							double [] XY2={imgData[2][vu+uvInc[cycle[2]]],imgData[3][vu+uvInc[cycle[2]]]};
							double [] V= {pX-XY0[0],pY-XY0[1]};
							double [][] M={
									{XY1[0]-XY0[0],XY2[0]-XY1[0]},
									{XY1[1]-XY0[1],XY2[1]-XY1[1]}};
							double det=M[0][0]*M[1][1]-M[1][0]*M[0][1];
							double [][] MInverse={
									{ M[1][1]/det,-M[0][1]/det},
									{-M[1][0]/det, M[0][0]/det}};
							double [] a12={
									MInverse[0][0]*V[0]+MInverse[0][1]*V[1],
									MInverse[1][0]*V[0]+MInverse[1][1]*V[1]};
							int pCorrIndex=idY*sWidth+idX;
// some points may be accumulated multiple times - thisPCorr[3] will take care of this
							if (this.debugLevel>3) {
								System.out.println("XY0="+IJ.d2s(XY0[0],3)+":"+IJ.d2s(XY0[1],3));
								System.out.println("XY1="+IJ.d2s(XY1[0],3)+":"+IJ.d2s(XY1[1],3));
								System.out.println("XY2="+IJ.d2s(XY2[0],3)+":"+IJ.d2s(XY2[1],3));
								System.out.println("M00="+IJ.d2s(M[0][0],3)+" M01="+IJ.d2s(M[0][1],3));
								System.out.println("M10="+IJ.d2s(M[1][0],3)+" M11="+IJ.d2s(M[1][1],3));
								System.out.println("MInverse00="+IJ.d2s(MInverse[0][0],5)+" MInverse01="+IJ.d2s(MInverse[0][1],5));
								System.out.println("MInverse10="+IJ.d2s(MInverse[1][0],5)+" MInverse11="+IJ.d2s(MInverse[1][1],5));
								System.out.println("a12="+IJ.d2s(a12[0],3)+":"+IJ.d2s(a12[1],3));
								System.out.println("imgData[0][vu+uvInc[cycle[0]]]="+IJ.d2s(imgData[0][vu+uvInc[cycle[0]]],3)+
										"imgData[1][vu+uvInc[cycle[0]]]="+IJ.d2s(imgData[1][vu+uvInc[cycle[0]]],3));
								System.out.println("imgData[0][vu+uvInc[cycle[1]]]="+IJ.d2s(imgData[0][vu+uvInc[cycle[1]]],3)+
										"imgData[1][vu+uvInc[cycle[1]]]="+IJ.d2s(imgData[1][vu+uvInc[cycle[1]]],3));
								System.out.println("imgData[0][vu+uvInc[cycle[2]]]="+IJ.d2s(imgData[0][vu+uvInc[cycle[2]]],3)+
										"imgData[1][vu+uvInc[cycle[2]]]="+IJ.d2s(imgData[1][vu+uvInc[cycle[2]]],3));
							}

							double [] corr={
									 imgData[0][vu+uvInc[cycle[0]]]+ // dPx
									(imgData[0][vu+uvInc[cycle[1]]]-imgData[0][vu+uvInc[cycle[0]]])*a12[0]+
									(imgData[0][vu+uvInc[cycle[2]]]-imgData[0][vu+uvInc[cycle[1]]])*a12[1],

									 imgData[1][vu+uvInc[cycle[0]]]+ // dPy
									(imgData[1][vu+uvInc[cycle[1]]]-imgData[1][vu+uvInc[cycle[0]]])*a12[0]+
									(imgData[1][vu+uvInc[cycle[2]]]-imgData[1][vu+uvInc[cycle[1]]])*a12[1],

									 imgData[4][vu+uvInc[cycle[0]]]+ // alpha
									(imgData[4][vu+uvInc[cycle[1]]]-imgData[4][vu+uvInc[cycle[0]]])*a12[0]+
									(imgData[4][vu+uvInc[cycle[2]]]-imgData[4][vu+uvInc[cycle[1]]])*a12[1],
									 imgData[5][vu+uvInc[cycle[0]]]+ // Red measured/pattern
									(imgData[5][vu+uvInc[cycle[1]]]-imgData[5][vu+uvInc[cycle[0]]])*a12[0]+
									(imgData[5][vu+uvInc[cycle[2]]]-imgData[5][vu+uvInc[cycle[1]]])*a12[1],
									 imgData[6][vu+uvInc[cycle[0]]]+ // Green measured/pattern
									(imgData[6][vu+uvInc[cycle[1]]]-imgData[6][vu+uvInc[cycle[0]]])*a12[0]+
									(imgData[6][vu+uvInc[cycle[2]]]-imgData[6][vu+uvInc[cycle[1]]])*a12[1],
									 imgData[7][vu+uvInc[cycle[0]]]+ // Blue  measured/pattern
									(imgData[7][vu+uvInc[cycle[1]]]-imgData[7][vu+uvInc[cycle[0]]])*a12[0]+
									(imgData[7][vu+uvInc[cycle[2]]]-imgData[7][vu+uvInc[cycle[1]]])*a12[1]
									};
							if (this.debugLevel>3) {
								System.out.println("corr="+IJ.d2s(corr[0],3)+" "+IJ.d2s(corr[1],3)+" "+IJ.d2s(corr[2],3));
							}
 if (pCorrIndex>thisPCorr[0].length) {
	 System.out.println("imgNum=" + imgNum+": "+	fittingStrategy.distortionCalibrationData.gIP[imgNum].path);
	 System.out.println("thisPCorr[0].length="+thisPCorr[0].length+" pCorrIndex="+pCorrIndex+" sWidth="+sWidth+" idY="+idY+" idX="+idX);
 }
							thisPCorr[0][pCorrIndex]+= corr[0];// dPx
							thisPCorr[1][pCorrIndex]+= corr[1];// dPy
							thisPCorr[2][pCorrIndex]+= corr[2];// alpha
							thisPCorr[3][pCorrIndex]+= 1.0;    // number of times accumulated
							thisPCorr[4][pCorrIndex]+= corr[3];// red
							thisPCorr[5][pCorrIndex]+= corr[4];// green
							thisPCorr[6][pCorrIndex]+= corr[5];// blue

							if (this.debugLevel>3) {
								debugHadPixels=true;
//								if (!debugExit) debugCntr--;
//								if (debugCntr==0) debugExit=true; // exit after first non-empty tile
							}

//gridPCorr[chnNum]
						}
					} // idX
					// use same order in calculations, make sure no gaps
				} // idY
				if ((this.debugLevel>3) && (debugHadPixels)){
					if (!debugExit) {
						System.out.println(
								" minX="+IJ.d2s(minX,1)+
								" maxX="+IJ.d2s(maxX,1));
						System.out.println(
								" minY="+IJ.d2s(minY,1)+
								" maxY="+IJ.d2s(maxY,1));
						System.out.println(
								" iMinX="+iMinX+
								" iMaxX="+iMaxX);
						System.out.println(
								" iMinY="+iMinY+
								" iMaxY="+iMaxY);
					}
					if (!debugExit) debugCntr--;
					if (debugCntr==0) debugExit=true; // exit after first non-empty tile

				}
			} // finished image


/*			if (showIndividual) {
				String [] titles={"dPx","dPy","alpha","Multiple","Red","Green","Blue"};
				ShowDoubleFloatArrays.showArrays(thisPCorr, sWidth, sHeight,  true, "thisPCorr_pre"+imgNum, titles);
			}
*/
			// some points may be calculated multiple times
			for (int i=0;i<gridPCorr[chnNum][0].length;i++) if (thisPCorr[3][i]>=1.0){
				thisPCorr[0][i]/=thisPCorr[3][i]; // dPx
				thisPCorr[1][i]/=thisPCorr[3][i]; // dPy
				thisPCorr[2][i]/=thisPCorr[3][i]; // alpha
				thisPCorr[4][i]/=thisPCorr[3][i]; // r
				thisPCorr[5][i]/=thisPCorr[3][i]; // g
				thisPCorr[6][i]/=thisPCorr[3][i]; // b
			}

			if (showIndividual && ((showIndividualNumber<0) || (showIndividualNumber==chnNum))) {
				String [] titles={"dPx","dPy","alpha","Multiple","Red","Green","Blue"};
				ShowDoubleFloatArrays.showArrays(thisPCorr, sWidth, sHeight,  true, "thisPCorr"+imgNum, titles);
			}
			for (int i=0;i<gridPCorr[chnNum][0].length;i++) if (thisPCorr[2][i]>0){
				gridPCorr[chnNum][0][i]+=thisPCorr[0][i]*thisPCorr[2][i];
				gridPCorr[chnNum][1][i]+=thisPCorr[1][i]*thisPCorr[2][i];
				/**TODO: not used anyway - just for debugging? see if just the sensor mask should go here? Or when saving?*/
				if (gridPCorr[chnNum][2][i]<thisPCorr[2][i]) gridPCorr[chnNum][2][i]=thisPCorr[2][i]; // best alpha
				gridPCorr[chnNum][3][i]+=                thisPCorr[2][i]; // sum of weights from all images
				gridPCorr[chnNum][4][i]+=thisPCorr[4][i]*thisPCorr[2][i];
				gridPCorr[chnNum][5][i]+=thisPCorr[5][i]*thisPCorr[2][i];
				gridPCorr[chnNum][6][i]+=thisPCorr[6][i]*thisPCorr[2][i];
			}
			IJ.showProgress(++numProcessed, numSelected);
		}
/*
		if (showIndividual) {
			String [] titles={"dPx","dPy","alpha","Multiple","Red","Green","Blue"};
			for (int chnNum=0;chnNum<gridPCorr.length;chnNum++) if (gridPCorr[chnNum]!=null) ShowDoubleFloatArrays.showArrays(gridPCorr[chnNum], sWidth, sHeight,  true, "gridPCorr1"+chnNum, titles);
		}
*/
		for (int chnNum=0;chnNum<gridPCorr.length;chnNum++) if (gridPCorr[chnNum]!=null){
			for (int i=0;i<gridPCorr[chnNum][0].length;i++) if (gridPCorr[chnNum][2][i]>0){ //null pointer
				gridPCorr[chnNum][0][i]/=gridPCorr[chnNum][3][i];
				gridPCorr[chnNum][1][i]/=gridPCorr[chnNum][3][i];
				gridPCorr[chnNum][4][i]/=gridPCorr[chnNum][3][i];
				gridPCorr[chnNum][5][i]/=gridPCorr[chnNum][3][i];
				gridPCorr[chnNum][6][i]/=gridPCorr[chnNum][3][i];
			}
		}
/*
		if (showIndividual) {
			String [] titles={"dPx","dPy","alpha","Multiple","Red","Green","Blue"};
			for (int chnNum=0;chnNum<gridPCorr.length;chnNum++) if (gridPCorr[chnNum]!=null) ShowDoubleFloatArrays.showArrays(gridPCorr[chnNum], sWidth, sHeight,  true, "gridPCorr2"+chnNum, titles);
		}
*/
		return gridPCorr;
	}

/**
 * Calculate partial derivative analytically (as the Jacobian calculation) and by difference divided by delta and compare
 * Done to debug derivatives calculation
 */
	public void compareDerivatives(){
		if (fittingStrategy==null) {
			String msg="Fitting strategy does not exist, exiting";
			IJ.showMessage("Error",msg);
			throw new IllegalArgumentException (msg);
		}
    	int numSeries=fittingStrategy.getNumSeries();
	    GenericDialog gd = new GenericDialog("Debug: verifying partial derivatives calculation, select series number");
		gd.addNumericField("Series number to show (0.."+(numSeries-1), this.seriesNumber, 0);
		gd.addCheckbox("Show actual parameters (false: X0,Y0,distance, angles)", true);
		gd.addCheckbox("Apply sensor mask (fade near edges)", true);
		gd.addCheckbox("Debug derivatives (show analytic/difference match)",true);
	    gd.showDialog();
	    if (gd.wasCanceled()) return;
	    this.seriesNumber=     (int) gd.getNextNumber();
	    boolean useActualParameters=gd.getNextBoolean();
	    boolean applySensorMask=gd.getNextBoolean();
	    boolean debugDerivatives=gd.getNextBoolean();
	    // currently not possible to debug "internal" parameters, so
//	    debugDerivatives&=useActualParameters; //*******************
		initFittingSeries(false,filterForAll,this.seriesNumber);
		int numPars=this.currentVector.length;
    	String [] parameterNames;
    	String [] parameterUnits;

    	if (useActualParameters) {
    		parameterNames=new String[fittingStrategy.distortionCalibrationData.getNumDescriptions()];
    		parameterUnits=new String[fittingStrategy.distortionCalibrationData.getNumDescriptions()];
    		for (int i=0;i<parameterNames.length;i++){
    			// TODO: move to DdistortionCalibrationData methods()
    			parameterNames[i]=fittingStrategy.distortionCalibrationData.descrField(i,0);
    			parameterUnits[i]=fittingStrategy.distortionCalibrationData.descrField(i,2);
    		}
    	} else {
    		parameterNames=lensDistortionParameters.getAllNames();
    		parameterUnits=lensDistortionParameters.getAllUnits();
    	}
	    gd = new GenericDialog((debugDerivatives?"Debug: verifying partial derivatives calculation,":"Showing partial derivatives,") +" select parameter number");
	    if (useActualParameters) {
	    	for (int i=0;i<this.currentVector.length;i++){
	    		int parNum=fittingStrategy.parameterMap[i][1];
	    		int imgNum=fittingStrategy.parameterMap[i][0];
	    		gd.addMessage(i+": "+parameterNames[parNum]+
	    				"["+imgNum+"]("+parameterUnits[parNum]+") "+IJ.d2s(this.currentVector[i],3));
	    	}
			gd.addNumericField("Select parameter number (0.."+(numPars-1)+") from above", 0, 0);
	    } else {
	    	for (int i=0;i<parameterNames.length;i++){
	    		gd.addMessage(i+": "+parameterNames[i]+"("+parameterUnits[i]+") ");
	    	}
			gd.addNumericField("Select parameter number (0.."+(parameterNames.length-1)+") from above", 0, 0);

	    }
		if (debugDerivatives) gd.addNumericField("Select delta to increment selected parameter", .001, 5);
		if (debugDerivatives) gd.addCheckbox("Show inter-parameter derivatives matrix", true);
		WindowTools.addScrollBars(gd);
	    gd.showDialog();
	    if (gd.wasCanceled()) return;
	    int selectedParameter=     (int) gd.getNextNumber();
	    double delta=0;
	    if (debugDerivatives) delta=     gd.getNextNumber();
	    boolean showInterparameterDerivatives=false;
	    if (debugDerivatives) showInterparameterDerivatives=gd.getNextBoolean();
		double [] this_currentfX=null;
	    double [] d_derivative;
	    double [] d_delta=null;
	    String title;
	    if (useActualParameters) {
	    	this_currentfX=calculateFxAndJacobian(this.currentVector, true); // is it always true here (this.jacobian==null)
	    	d_derivative=this.jacobian[selectedParameter].clone(); //  wrong?
	    	if (debugDerivatives) {
	    		double[] modVector=this.currentVector.clone();
	    		modVector[selectedParameter]+=delta;
	    		d_delta=calculateFxAndJacobian(modVector, true);
	    		if (this.debugLevel>3) {
	    			for (int i=0;i<d_delta.length;i++) {
	    				System.out.println(i+": "+IJ.d2s(d_delta[i],3)+" - "+IJ.d2s(this_currentfX[i],3)+" = "
	    						+ IJ.d2s(d_delta[i]-this_currentfX[i],3));
	    			}

	    		}
	    		for (int i=0;i<d_delta.length;i++) d_delta[i]= (d_delta[i]-this_currentfX[i])/delta;
	    	}
		    int parNum=fittingStrategy.parameterMap[selectedParameter][1];
			int imgNum=fittingStrategy.parameterMap[selectedParameter][0];
			title=parameterNames[parNum]+"_derivatives:"+imgNum;
	    } else {
	    	d_derivative=calculateJacobian16(this.currentVector, -1,0.0)[selectedParameter].clone();
	    	if (debugDerivatives) d_delta=     calculateJacobian16(this.currentVector, -1,delta)[selectedParameter].clone();
			title=parameterNames[selectedParameter]+"_derivatives";
	    }
	    if (this.debugLevel>3) {
		    for (int i=0;i<d_delta.length;i++) {
		    	System.out.println(i+":: "+IJ.d2s(d_delta[i],3)+" - "+IJ.d2s(d_derivative[i],3));
		    }
	    }
	    double [] sumWeight=showCompareDerivatives (d_derivative, d_delta, applySensorMask, !useActualParameters,  title ); // d_delta==null - no debug
	    if (showInterparameterDerivatives && (delta>0)) {
	    debugCompareInterparameterDerivatives(
	    		this.currentVector.clone(),
	    		-1, //int imgNum,
	    		delta);
	    }
	    for (int i=0;i<sumWeight.length; i++) if (sumWeight[i]>0.0){
	    	System.out.println("Image "+i+", "+title+"derivative RMS="+sumWeight[i]);
	    }
	}


	/**
	 * Show comparison of the calculated partial derivatives in Jacobian and approximated by difference
	 * for incremented parameters
	 * @param imgNumber - number of image in series to show
	 * @param d_derivative vector array of "true" derivatives (from Jacobian)
	 * @param d_delta approximated derivatives from varying parameter
	 * @param title image title
	 * @return rms
	 */
	public double showCompareDerivatives(int imgNumber, double [] d_derivative, double [] d_delta, boolean applySensorMask, String title ){
		String [] titlesDebug={"dX-derivative","dY-derivative","abs-derivative","diff-X (should be 0)","diff-Y (should be 0)","dX-delta/delta","dY-delta/delta","dX-delta","dY-delta"};
		String [] titlesNoDebug={"dX-derivative","dY-derivative","abs-derivative"};
		String [] titles= (d_delta==null)? titlesNoDebug:titlesDebug;
		double [] d_diff=new double [d_derivative.length];
		double [] r_diff=new double [d_derivative.length];
		double [] aDeriv=new double [d_derivative.length/2];

		if (d_delta!=null) for (int i=0;i<d_diff.length;i++){
			d_diff[i]=d_derivative[i]-d_delta[i];
			r_diff[i]=d_diff[i]/d_delta[i];
		}
// find data range for the selected image
		int index=0;
		int numImg=fittingStrategy.distortionCalibrationData.getNumImages();
		boolean [] selectedImages=fittingStrategy.selectedImages();
		for (int imgNum=0;(imgNum<imgNumber) && (imgNum<numImg) ;imgNum++) if (selectedImages[imgNum])
			index+=fittingStrategy.distortionCalibrationData.gIP[imgNum].pixelsUV.length;
		double sumWeights=0.0;
		double sumDerivatives2=0.0;
		double w,sqrtW;
		for (int i=2*index;i<2*(2*index+fittingStrategy.distortionCalibrationData.gIP[imgNumber].pixelsUV.length);i++){
			w=applySensorMask?this.weightFunction[i]:1.0;
			if (w<0.0) w=0.0;
			sumWeights+=w;
			sumDerivatives2+=d_derivative[i]*d_derivative[i]*w;
			sqrtW=Math.sqrt(w);
			d_derivative[i]*=sqrtW; // for display
			if (d_delta!=null) d_delta[i]*=sqrtW;
			if ((i&1)==0) aDeriv[i>>1]=Math.sqrt(d_derivative[i]*d_derivative[i]+d_derivative[i+1]*d_derivative[i+1]);
		}
		sumDerivatives2=Math.sqrt(sumDerivatives2/sumWeights*2.0); // 2.0 because x,y pair should not be averaged, just added
		titles[2]+=":rms="+sumDerivatives2;
		int width=getGridWidth();
		double [][] imgData=new double[titles.length][getGridHeight() * width];
		for (int i=0;i<imgData.length;i++) for (int j=0;j<imgData[i].length;j++)imgData[i][j]=0.0;

		for (int i=0;i<fittingStrategy.distortionCalibrationData.gIP[imgNumber].pixelsUV.length;i++){
			int u=fittingStrategy.distortionCalibrationData.gIP[imgNumber].pixelsUV[index+i][0]+patternParameters.U0;
			int v=fittingStrategy.distortionCalibrationData.gIP[imgNumber].pixelsUV[index+i][1]+patternParameters.V0;
			int vu=u+width*v;
			imgData[0][vu]=   d_derivative[2*(index+i)];
			imgData[1][vu]=   d_derivative[2*(index+i)+1];
			imgData[2][vu]=   aDeriv[index+i];
			if (d_delta!=null) {
				imgData[3][vu]=   d_diff[2*(index+i)];
				imgData[4][vu]=   d_diff[2*(index+i)+1];
				imgData[5][vu]=   r_diff[2*(index+i)];
				imgData[6][vu]=   r_diff[2*(index+i)+1];
				imgData[7][vu]=   d_delta[2*(index+i)];
				imgData[8][vu]=   d_delta[2*(index+i)+1];
			}
		}
		ShowDoubleFloatArrays.showArrays(imgData, width, getGridHeight(),  true, title, titles);
		return sumDerivatives2;
	}
	/**
	 * Show comparison of the calculated partial derivatives in Jacobian and approximated by difference
	 * for incremented parameters (for all selected images in the series)
	 * @param d_derivative vector array of "true" derivatives (from Jacobian)
	 * @param d_delta approximated derivatives form varying parameter
	 * @param applySensorMask Multiply by sensor mask (fade near edges)
	 * @param single calculate just the first selected image
	 * @param title image title
	 * @return array of rms
	 */

	public double[] showCompareDerivatives (double [] d_derivative, double [] d_delta, boolean applySensorMask, boolean single, String title ){
		boolean [] selectedImages=fittingStrategy.selectedImages();
		double [] diffs= new double [selectedImages.length];
		for (int imgNum=0;imgNum<diffs.length;imgNum++) diffs[imgNum]=0.0;
		for (int imgNum=0;imgNum<selectedImages.length;imgNum++) if (selectedImages[imgNum]) {
			diffs[imgNum] =showCompareDerivatives(imgNum, d_derivative, d_delta, applySensorMask, title+"-"+imgNum);
			if (single) break;
		}
		return diffs;
	}

	/**
	 *
	 * @param delta if 0 - actual derivatives, >0 - approximate derivatives by deltas
	 * @return for each v,u - values and derivatives
	 */
	public double [][][][] calcGridOnSensor( double delta) {
    	int gridHeight=patternParameters.gridGeometry.length;
    	int gridWidth=patternParameters.gridGeometry[0].length;
    	this.gridOnSensor=new double[gridHeight][gridWidth][2][15];
    	double [][] node;
 //   	double [][][] nodes=new double [15][][];
    	boolean dMode=delta>0;
        if (this.debugLevel>2){
        	System.out.println("calcGridOnSensor()");
        	System.out.println("this.lensDistortionParameters.distance="+IJ.d2s(this.lensDistortionParameters.distance, 3));
        	System.out.println("this.lensDistortionParameters.x0="+      IJ.d2s(this.lensDistortionParameters.x0, 3));
        	System.out.println("this.lensDistortionParameters.y0="+      IJ.d2s(this.lensDistortionParameters.y0, 3));
        	System.out.println("this.lensDistortionParameters.z0="+      IJ.d2s(this.lensDistortionParameters.z0, 3));
        	System.out.println("this.lensDistortionParameters.pitch="+   IJ.d2s(this.lensDistortionParameters.pitch, 3));
        	System.out.println("this.lensDistortionParameters.yaw="+IJ.d2s(this.lensDistortionParameters.yaw, 3));
        	System.out.println("this.lensDistortionParameters.roll="+IJ.d2s(this.lensDistortionParameters.roll, 3));
        	System.out.println("this.lensDistortionParameters.focalLength="+IJ.d2s(this.lensDistortionParameters.focalLength, 3));
        	System.out.println("this.lensDistortionParameters.px0="+IJ.d2s(this.lensDistortionParameters.px0, 3));
        	System.out.println("this.lensDistortionParameters.py0="+IJ.d2s(this.lensDistortionParameters.py0, 3));
        	System.out.println("this.lensDistortionParameters.distortionA8="+IJ.d2s(this.lensDistortionParameters.distortionA8, 5));
        	System.out.println("this.lensDistortionParameters.distortionA7="+IJ.d2s(this.lensDistortionParameters.distortionA7, 5));
        	System.out.println("this.lensDistortionParameters.distortionA6="+IJ.d2s(this.lensDistortionParameters.distortionA6, 5));
        	System.out.println("this.lensDistortionParameters.distortionA5="+IJ.d2s(this.lensDistortionParameters.distortionA5, 5));
        	System.out.println("this.lensDistortionParameters.distortionA="+IJ.d2s(this.lensDistortionParameters.distortionA, 5));
        	System.out.println("this.lensDistortionParameters.distortionB="+IJ.d2s(this.lensDistortionParameters.distortionB, 5));
        	System.out.println("this.lensDistortionParameters.distortionC="+IJ.d2s(this.lensDistortionParameters.distortionC, 5));
        	System.out.println("this.lensDistortionParameters.lensDistortionModel="+this.lensDistortionParameters.lensDistortionModel);
        	for (int i=0;i<this.lensDistortionParameters.r_xy.length;i++){
            	System.out.println("this.lensDistortionParameters.r_xy["+i+"][0]="+IJ.d2s(this.lensDistortionParameters.r_xy[i][0], 5));
            	System.out.println("this.lensDistortionParameters.r_xy["+i+"][1]="+IJ.d2s(this.lensDistortionParameters.r_xy[i][1], 5));
        	}
        	for (int i=0;i<this.lensDistortionParameters.r_od.length;i++){
            	System.out.println("this.lensDistortionParameters.r_od["+i+"][0]="+IJ.d2s(this.lensDistortionParameters.r_od[i][0], 5));
            	System.out.println("this.lensDistortionParameters.r_od["+i+"][1]="+IJ.d2s(this.lensDistortionParameters.r_od[i][1], 5));
        	}
        }
        LensDistortionParameters ldp=this.lensDistortionParameters.clone();
        // 06/2019 - need to update distortionRadius, pixelSize)

//		public void setLensDistortionParameters(LensDistortionParameters ldp

        for (int v=0; v<gridHeight; v++) for (int u=0; u<gridWidth; u++) if (patternParameters.gridGeometry[v][u][3]>0) {
        	this.lensDistortionParameters.setLensDistortionParameters(ldp); // restore
        	node=this.lensDistortionParameters.calcPartialDerivatives(
        			patternParameters.gridGeometry[v][u][0],//double xp, // target point horizontal, positive - right,  mm
        			patternParameters.gridGeometry[v][u][1],//double yp, // target point vertical,   positive - down,  mm
        			patternParameters.gridGeometry[v][u][2],//double zp, // target point horizontal, positive - away from camera,  mm
        			!dMode);//boolean calculateAll){ // calculate derivatives, false - values only
        	if (this.debugLevel>3) {
        		System.out.println("calcPartialDerivatives("+
        				IJ.d2s(patternParameters.gridGeometry[v][u][0],2)+","+
        				IJ.d2s(patternParameters.gridGeometry[v][u][1],2)+","+
        				IJ.d2s(patternParameters.gridGeometry[v][u][2],2)+" ("+true+") -> "+
        				IJ.d2s(node[0][0],2)+"/"+IJ.d2s(node[0][1],2));
        	}
        	if (dMode) {
//        		double []pXY=node[0]; // px,py values
        		this.gridOnSensor[v][u][0][0]=node[0][0];
        		this.gridOnSensor[v][u][1][0]=node[0][1];
        		for (int j=1;j<15;j++) {  // was 14
        			this.lensDistortionParameters.setLensDistortionParameters(ldp, j, delta); // set one of the parameters (j) with added delta to ldp
                	node=this.lensDistortionParameters.calcPartialDerivatives(
                			patternParameters.gridGeometry[v][u][0],//double xp, // target point horizontal, positive - right,  mm
                			patternParameters.gridGeometry[v][u][1],//double yp, // target point vertical,   positive - down,  mm
                			patternParameters.gridGeometry[v][u][2],//double zp, // target point horizontal, positive - away from camera,  mm
                			false);
            		this.gridOnSensor[v][u][0][j]=(node[0][0]-this.gridOnSensor[v][u][0][0])/delta;
            		this.gridOnSensor[v][u][1][j]=(node[0][1]-this.gridOnSensor[v][u][1][0])/delta;
        		}

        	} else for (int i=0;i<2;i++) for (int j=0;j<15;j++){ // was 14
        		this.gridOnSensor[v][u][i][j]=node[j][i];
        	}

        } else {
        	this.gridOnSensor[v][u]=null;
        }
        return this.gridOnSensor;
    }
    public int getGridWidth() {
    	return patternParameters.gridGeometry[0].length;
    }
    public int getGridHeight() {
    	return patternParameters.gridGeometry.length;
    }

    public double [][] prepareDisplayGrid(){
    	int gridHeight=this.patternParameters.gridGeometry.length;
    	int gridWidth=this.patternParameters.gridGeometry[0].length;
    	double [][] dgrid=new double[3][gridHeight*gridWidth];
    	double average;
    	int num,index;
    	for (int i=0;i<dgrid.length;i++){
    		average=0.0;
    		num=0;
    		for (int v=0; v<gridHeight; v++) for (int u=0; u<gridWidth; u++) if (this.patternParameters.gridGeometry[v][u][3]>0) {
    			average+=this.patternParameters.gridGeometry[v][u][i];
    			num++;
    		}
    		average/=num;
    		index=0;
    		for (int v=0; v<gridHeight; v++) for (int u=0; u<gridWidth; u++) if (this.patternParameters.gridGeometry[v][u][3]>0) {
    			dgrid[i][index++]=this.patternParameters.gridGeometry[v][u][i];
    		} else {
    			dgrid[i][index++]=average;
    		}
    	}
    	return dgrid;
    }
    public String [] displayGridTitles() {
    	String [] titles={"Grid-X","Grid-Y","Grid-Z"};
    	return titles;
    }
    public String [] displayGridOnSensorTitles() {
    	String [] titles={
    			"PX","PY",
    			"dPX/dphi","dPY/dphi",
    			"dPX/dtheta","dPY/dtheta",
    			"dPX/dpsi","dPY/dpsi",
    			"dPX/dX0","dPY/dX0",
    			"dPX/dY0","dPY/dY0",
    			"dPX/dZ0","dPY/dZ0",
    			"dPX/df","dPY/df",
    			"dPX/ddist","dPY/dist",
    			"dPX/dDa","dPY/dDa",
    			"dPX/dDb","dPY/dDb",
    			"dPX/dDc","dPY/dDc",
    			"dPX/dPX0","dPY/dPX0",
    			"dPX/dPY0","dPY/dPY0"
    	};
    	return titles;
    }
    public double [][] prepareDisplayGridOnSensor(boolean showAll){
    	int gridHeight=this.patternParameters.gridGeometry.length;
    	int gridWidth=this.patternParameters.gridGeometry[0].length;
//    	double [][] dgrid=new double[showAll?28:2][gridHeight*gridWidth];
    	double [][] dgrid=new double[showAll?(2*15):2][gridHeight*gridWidth];
    	double average;
    	int num,index;
    	for (int i=0;i<dgrid.length/2;i++) for (int j=0;j<2;j++){
    		int ii=i*2+j;
    		average=0.0;
    		num=0;
    		for (int v=0; v<gridHeight; v++) for (int u=0; u<gridWidth; u++) if (this.patternParameters.gridGeometry[v][u][3]>0) {
    			average+=this.gridOnSensor[v][u][j][i];
    			num++;
    		}
    		average/=num;
    		index=0;
    		for (int v=0; v<gridHeight; v++) for (int u=0; u<gridWidth; u++) if (this.patternParameters.gridGeometry[v][u][3]>0) {
    			dgrid[ii][index++]=this.gridOnSensor[v][u][j][i];
    		} else {
    			dgrid[ii][index++]=average;
    		}
    	}
    	return dgrid;
    }
    /**
     * initialize image data with camera defaults
     * @param distortionCalibrationData grid distortionCalibrationData
     * @param eyesisCameraParameters deafault camera parameters
     */

    // Used in Aberration_Calibration
    public void initImageSet(
    		DistortionCalibrationData distortionCalibrationData,
    		EyesisCameraParameters eyesisCameraParameters) {
//    	DistortionCalibrationData distortionCalibrationData= new DistortionCalibrationData(filenames);
    	for (int i=0;i<distortionCalibrationData.getNumImages();i++){
    		int stationNumber=distortionCalibrationData.getImageStation(i);
    		int subCam=distortionCalibrationData.getImageSubcamera(i);
    		distortionCalibrationData.setParameters(eyesisCameraParameters.getParametersVector(stationNumber,subCam), i);
    		this.lensDistortionParameters.pixelSize=eyesisCameraParameters.getPixelSize(subCam);
    		this.lensDistortionParameters.distortionRadius=eyesisCameraParameters.getDistortionRadius(subCam);
    	}
    }
    public void copySensorConstants(EyesisCameraParameters eyesisCameraParameters) { // copy from the first channel
    		this.lensDistortionParameters.pixelSize=eyesisCameraParameters.getPixelSize(0);
    		this.lensDistortionParameters.distortionRadius=eyesisCameraParameters.getDistortionRadius(0);
    }

    /**
     * Update per-image parameters from those of the camera and those that have the same timestamp. Usually needed after adding or
     * enabling new images.
     * @param distortionCalibrationData grid distortionCalibrationData
     * @param eyesisCameraParameters - camera parameters (common and per sub-camera)
     * @return true if dialog was not canceled and programs ran
     */

    public boolean interactiveUpdateImageSet(
    		DistortionCalibrationData distortionCalibrationData,
    		EyesisCameraParameters eyesisCameraParameters
    ){
    	boolean resetParametersToZero=false;
    	boolean updateAllSubcameras = false;
    	boolean [] parameterMask= new boolean[distortionCalibrationData.getNumParameters()];
    	boolean [] channelMask=   new boolean[distortionCalibrationData.getNumSubCameras()];
    	boolean [] stationMask=   new boolean[distortionCalibrationData.getNumStations()];
    	String [] source_stations = new String[stationMask.length+1];
    	source_stations[0] = "---";
    	for (int i = 0; i < stationMask.length; i++) {
    		source_stations[i+1]=""+i;
    	}
    	for (int i=0;i<parameterMask.length;i++) parameterMask[i]=false;
    	for (int i=0;i<channelMask.length;i++)   channelMask[i]=  true;
    	for (int i=0;i<stationMask.length;i++)   stationMask[i]=  true;
    	GenericDialog gd=new GenericDialog("Update (new) image settings from known data");
    	//
    	gd.addCheckbox("Reset selected parameters to zero (false - update from camera parameters)", resetParametersToZero);
    	if (stationMask.length > 1) {
    		gd.addChoice("Copy selected parameters from tis station to all other stations" , source_stations, source_stations[0]);
    	}
    	gd.addMessage("Select which individual image parameters to be updated from the camera parameters (or reset to 0)");
    	
    	gd.addCheckbox("Update all subcamera parameters", updateAllSubcameras);
    	
    	for (int i=0;i<parameterMask.length;i++) {
    		gd.addCheckbox(i+": "+distortionCalibrationData.getParameterName(i), parameterMask[i]);
    	}
    	gd.addMessage("----------");
    	gd.addMessage("Select which channels (sub-cameras) to update");
    	for (int i=0;i<channelMask.length;i++) gd.addCheckbox("Subcamera "+i, channelMask[i]);
    	if (stationMask.length>1) {
        	gd.addMessage("----------");
        	gd.addMessage("Select which stations (camera/goniometer locations) to update");
        	for (int i=0;i<stationMask.length;i++) gd.addCheckbox("Station "+i, stationMask[i]);
    	}
    	gd.addMessage("----------");
    	gd.addCheckbox("Applying known extrinsic parameters to the same timestamp images", true);
    	gd.addCheckbox("Use closest (by motor steps) image if none for the same timestamp is enabled", true);
    	gd.addMessage("==== Note: The following correction will be applied to all subcameras, use selection above to specify which heights should be averaged" );
    	gd.addCheckbox("Vertically center the camera head by calculateing center above horizontal", false);
//    	gd.addCheckbox("Update currently disabled images", true);
	    WindowTools.addScrollBars(gd);
    	gd.showDialog();
    	if (gd.wasCanceled()) return false;
    	resetParametersToZero=gd.getNextBoolean();
    	int source_station = -1;
    	if (stationMask.length > 1) {
    		source_station = gd.getNextChoiceIndex() - 1;
    	}

    	updateAllSubcameras = gd.getNextBoolean();
    	for (int i=0;i<parameterMask.length;i++) parameterMask[i]= gd.getNextBoolean();
    	for (int i=0;i<channelMask.length;i++)   channelMask[i]=   gd.getNextBoolean();
    	if (stationMask.length>1) {
    		for (int i=0;i<stationMask.length;i++) stationMask[i]= gd.getNextBoolean();
    	}
    	boolean updateFromTimestamps= gd.getNextBoolean();
    	boolean allowClosest=         gd.getNextBoolean();
    	boolean reCenterVertically=   gd.getNextBoolean();

    	if (updateAllSubcameras) {
    		resetParametersToZero = false; // just for safety
        	for (int i=0;i<parameterMask.length;i++) {
        		parameterMask[i] |= distortionCalibrationData.isSubcameraParameter(i);
        	}
    	}
    	
    	if (source_station >= 0) {
    		updateOtherStations(
    	    		eyesisCameraParameters,
    	    		source_station,
    	    		parameterMask,
    	    		channelMask,
    	    		stationMask);    		
    	}
    	
    	
    	
    	if (reCenterVertically){
    		eyesisCameraParameters.recenterVertically(channelMask, stationMask);
    		for (int i=0;i<channelMask.length;i++) channelMask[i]= true;
    		parameterMask[distortionCalibrationData.getParameterIndexByName("subcamHeight")] = true;
    	}


//		boolean updateDisabled=       gd.getNextBoolean();
    	updateImageSetFromCamera(
    			resetParametersToZero,
    			distortionCalibrationData,
    			eyesisCameraParameters,
    			parameterMask, //boolean [] parameterMask,
    			channelMask, // copy X,Y,Z (usually true)
    			stationMask // copy 2 goniometer angles (usually false)
    	);
    	if (updateFromTimestamps) {
    		updateImageSetFromSameTimestamps(
    				distortionCalibrationData,
    				eyesisCameraParameters,
    				null, // boolean [] selectedImages,
    				null, //boolean [] parameterMask,
    				allowClosest
    				//,updateDisabled
    		);
    		distortionCalibrationData.updateSetOrientation(null); // update orientation of image sets
    	}
    	return true;
    }

    public boolean setSetFromClosestAndEstimateOrientation(
    		int numSet,
    		boolean [] selectedImages,
    		boolean [] parameterMask,
    		DistortionCalibrationData distortionCalibrationData,
    		EyesisCameraParameters eyesisCameraParameters){
    	if (selectedImages==null) {
    		selectedImages= new boolean[distortionCalibrationData.getNumImages()];
    		for (int i=0;i<selectedImages.length;i++) selectedImages[i]=distortionCalibrationData.gIP[i].enabled;
    	}
    	if (parameterMask==null) {
    		parameterMask= new boolean[distortionCalibrationData.getNumParameters()];
    		for (int i=0;i<parameterMask.length;i++) parameterMask[i]=true;
    	}
    	for (int i=0;i<parameterMask.length;i++) {
    		if (distortionCalibrationData.isSubcameraParameter(i))    	parameterMask[i]=false;
    	}

    	int enabledImage=getClosestImage( // {numEnabledSet,enabledChannel,enabledImage};
	    		distortionCalibrationData,
	    		selectedImages,
	    		numSet);
    	if (enabledImage<0) return false; // failed to find closest
		updateSetFromClosest(
				numSet,
				enabledImage,
				parameterMask,
				distortionCalibrationData);
		// invalidate current angles
		distortionCalibrationData.gIS[numSet].goniometerAxial=Double.NaN;
		distortionCalibrationData.gIS[numSet].goniometerTilt= Double.NaN;
		// re-estimate orientation
		double [] ta=distortionCalibrationData.getImagesetTiltAxial(distortionCalibrationData.gIS[numSet].timeStamp); // updates tilt/axial (now interAxis too!)
	    if ((ta==null) || Double.isNaN(ta[0]) || Double.isNaN(ta[1])) return false;
	    return true;
    }


    public boolean interactiveUpdateImageSetOld(
    		DistortionCalibrationData distortionCalibrationData,
    		EyesisCameraParameters eyesisCameraParameters
    ){
    	GenericDialog gd=new GenericDialog("Update (new) image settings from known data");
    	gd.addCheckbox("Update per-image parameters from those of the camera", true);
    	gd.addCheckbox("Copy location of the camera (X,Y,Z)", true);
    	gd.addCheckbox("Copy orientation of the camera (tilt and axial)", false);
    	gd.addMessage("");
    	gd.addCheckbox("Update per-image parameters from those with the same timestamp", true);
    	gd.addCheckbox("Use closest (by motor steps) image if none for the same timestamp is enabled", true);
//    	gd.addCheckbox("Update currently disabled images", true);

    	gd.showDialog();
    	if (gd.wasCanceled()) return false;
    	boolean updateFromCamera=     gd.getNextBoolean();
    	boolean copyLocation=         gd.getNextBoolean();
    	boolean copyOrientation=      gd.getNextBoolean();
    	boolean updateFromTimestamps= gd.getNextBoolean();
    	boolean allowClosest=         gd.getNextBoolean();
//		boolean updateDisabled=       gd.getNextBoolean();

    	boolean [] parameterMask= new boolean[distortionCalibrationData.getNumParameters()];
    	for (int i=0;i<parameterMask.length;i++) {
    		parameterMask[i]=true;
    		if (distortionCalibrationData.isLocationParameter(i)    && !copyLocation)    	parameterMask[i]=false;
    		if (distortionCalibrationData.isOrientationParameter(i) && !copyOrientation)	parameterMask[i]=false;
    	}

    	if (updateFromCamera) updateImageSetFromCamera(
    			false, //resetParametersToZero
    			distortionCalibrationData,
    			eyesisCameraParameters,
    			parameterMask, //boolean [] parameterMask,
    			null,
    			null
    	);
    	if (updateFromTimestamps) {
    		updateImageSetFromSameTimestamps(
    				distortionCalibrationData,
    				eyesisCameraParameters,
    				null, // boolean [] selectedImages,
    				null, //boolean [] parameterMask,
    				allowClosest
//    				,updateDisabled
    		);
    		distortionCalibrationData.updateSetOrientation(null); // update orientation of image sets
    	}
    	return true;
    }

    
    /**
     * Update selected parameters from sourceStation to selected (stationMask) stations, filtered by channelMask
     * @param eyesisCameraParameters
     * @param sourceStation
     * @param parameterMask
     * @param channelMask
     * @param stationMask
     */
    public void updateOtherStations(
    		EyesisCameraParameters eyesisCameraParameters,
    		int                    sourceStation,
    		boolean []             parameterMask,
    		boolean []             channelMask,
    		boolean []             stationMask
    		) {
    	for (int stationNumber = 0; stationNumber <  stationMask.length; stationNumber++) if (stationMask[stationNumber]) {
    		for (int subCam=0; subCam < channelMask.length; subCam++) if (channelMask[subCam]) {
        		double [] oldVector=eyesisCameraParameters.getParametersVector(stationNumber,subCam);
        		double [] newVector=eyesisCameraParameters.getParametersVector(sourceStation,subCam);
        		for (int j=0;j<oldVector.length;j++) if (parameterMask[j]){
        			oldVector[j]=newVector[j];
        		}
    			eyesisCameraParameters.setParametersVector(
    					newVector,
    					parameterMask,
    					stationNumber,
    					subCam);
    		}
    	}
    }
    

    /**
     * Copies selected parameters from the camera parameters to per-image parameters (i.e. for new/previously disabled images)
     * @param distortionCalibrationData grid distortionCalibrationData
     * @param eyesisCameraParameters - camera parameters (common and per sub-camera)
     * @param parameterMask when element is true - copy parameters, false - keep current value. Null - selects all (filtered by the next parameters)
     * @param copyLocation copy location (x,Y,Z) of the camera , normally should be true
     * @param copyOrientation copy 2 goniometer angles, normally should be false
     */
    public void updateImageSetFromCamera(
    		boolean resetParametersToZero, // reset to 0 instead of camera parameters
    		DistortionCalibrationData distortionCalibrationData,
    		EyesisCameraParameters eyesisCameraParameters,
    		boolean [] parameterMask,
    		boolean [] channelMask,
    		boolean [] stationMask
    		) {
//    	DistortionCalibrationData distortionCalibrationData= new DistortionCalibrationData(filenames);
    	for (int i=0;i<distortionCalibrationData.getNumImages();i++){
    		int stationNumber=distortionCalibrationData.getImageStation(i);
    		int subCam=distortionCalibrationData.getImageSubcamera(i);
    		if ((channelMask!=null) && !channelMask[subCam])        continue;
    		if ((stationMask!=null) && !stationMask[stationNumber]) continue;
    		double [] oldVector=distortionCalibrationData.getParameters(i);
    		double [] newVector=eyesisCameraParameters.getParametersVector(stationNumber,subCam);
    		for (int j=0;j<oldVector.length;j++) if (parameterMask[j]){
    			if (resetParametersToZero) newVector[j]=0.0;
    			oldVector[j]=newVector[j];
    		}
    		if (resetParametersToZero){
    			eyesisCameraParameters.setParametersVector(
    					newVector,
    					parameterMask,
    					stationNumber,
    					subCam);
    		}
    		distortionCalibrationData.setParameters(oldVector, i);
    		this.lensDistortionParameters.pixelSize=eyesisCameraParameters.getPixelSize(subCam);
    		this.lensDistortionParameters.distortionRadius=eyesisCameraParameters.getDistortionRadius(subCam);
    	}
    }
    /**
     * Copies selected (normally all) parameters from the selected images with the same timestamp (i.e. for new/previously disabled images)
     * @param distortionCalibrationData grid distortionCalibrationData
     * @param eyesisCameraParameters - camera parameters (common and per sub-camera)
     * @param selectedImages Use only selected images (null - all enabled)
     * @param parameterMask when element is true - copy parameters, false - keep current value. Null - selects all (and should be normally null)
     * @param allowClosest If there is no enabled image for the current timestamp, find the closest selected using motor coordinates
     * @param updateDisabled update disable images also
     */
    public void updateImageSetFromSameTimestamps(
    		DistortionCalibrationData distortionCalibrationData,
    		EyesisCameraParameters eyesisCameraParameters,
    		boolean [] selectedImages,
    		boolean [] parameterMask,
    		boolean allowClosest
    		){
		System.out.println("updateImageSetFromSameTimestamps(), allowClosest="+allowClosest); //+" updateDisabled="+updateDisabled);
    	if (selectedImages==null) {
    		selectedImages= new boolean[distortionCalibrationData.getNumImages()];
    		for (int i=0;i<selectedImages.length;i++) selectedImages[i]=distortionCalibrationData.gIP[i].enabled;
//    		for (int i=0;i<selectedImages.length;i++) selectedImages[i]=distortionCalibrationData.gIP[i].enabled || updateDisabled;
    	}
    	if (parameterMask==null) {
    		parameterMask= new boolean[distortionCalibrationData.getNumParameters()];
    		for (int i=0;i<parameterMask.length;i++) parameterMask[i]=true;
    	}
    	for (int i=0;i<parameterMask.length;i++) {
    		if (distortionCalibrationData.isSubcameraParameter(i))    	parameterMask[i]=false;
    	}
    	for (int numSet=0; numSet<distortionCalibrationData.gIS.length;numSet++){
// find enabled image for this set
    		int enabledImage=-1;
    		// look for enabled image in the same imageSet

    		for (int nChn=0;nChn<distortionCalibrationData.gIS[numSet].imageSet.length;nChn++) if (distortionCalibrationData.gIS[numSet].imageSet[nChn]!=null){
    			int img=distortionCalibrationData.gIS[numSet].imageSet[nChn].imgNumber;
    			if (selectedImages[img]){
        			enabledImage=img;
    				break;
    			}
    		}
    		// look for closest in the other imageSet
    		if ((enabledImage<0) && (allowClosest)){
    			enabledImage=getClosestImage( // {numEnabledSet,enabledChannel,enabledImage};
    			    		distortionCalibrationData,
    			    		selectedImages,
    			    		numSet);
    		}
    		if (enabledImage>=0){
    			updateSetFromClosest(
    					numSet,
    					enabledImage,
    					parameterMask,
    					distortionCalibrationData);
    		}
    	}

    }

    public int getClosestImage(
    		DistortionCalibrationData distortionCalibrationData,
    		boolean [] selectedImages,
    		int numSet
    ){
    	int enabledChannel=-1;
    	int enabledImage=-1;
    	if (distortionCalibrationData.gIS[numSet].motors==null ){
    		if (this.debugLevel>0) System.out.println("getClosestSetChannelImage(): No motor data for timestamp "+distortionCalibrationData.gIS[numSet].timeStamp);
    		return -1;
    	}
    	double d2Min=-1;
    	for (int numOtherSet=0;numOtherSet<distortionCalibrationData.gIS.length;numOtherSet++)
    		if ((numOtherSet!=numSet) &&
    				(distortionCalibrationData.gIS[numOtherSet].stationNumber==distortionCalibrationData.gIS[numSet].stationNumber) &&
    				(distortionCalibrationData.gIS[numOtherSet].motors!=null) &&
    				(distortionCalibrationData.gIS[numOtherSet].imageSet!=null)
    		) {
    			enabledChannel=-1;
    			int otherImage=-1;
    			for (int nChn=0;nChn<distortionCalibrationData.gIS[numOtherSet].imageSet.length;nChn++)
    				if (distortionCalibrationData.gIS[numOtherSet].imageSet[nChn]!=null){
    					otherImage=distortionCalibrationData.gIS[numOtherSet].imageSet[nChn].imgNumber;
    					if (selectedImages[otherImage]){
    						enabledChannel=nChn;
    						break;
    					}
    				}
    			if (enabledChannel>=0){
    				double d2=0;
    				for (int k=0;k<distortionCalibrationData.gIS[numOtherSet].motors.length;k++){
    					d2+=1.0*(distortionCalibrationData.gIS[numOtherSet].motors[k]-distortionCalibrationData.gIS[numSet].motors[k])*
    					(distortionCalibrationData.gIS[numOtherSet].motors[k]-distortionCalibrationData.gIS[numSet].motors[k]);
    				}
    				if ((d2Min<0) || (d2Min>d2)) {
    					d2Min=d2;
    					enabledImage=otherImage;
    				}
    			}
    		}
    	return enabledImage;
    }

    public void updateSetFromClosest(
    		int numSet,
    		int enabledImage,
    		boolean [] parameterMask,
    		DistortionCalibrationData distortionCalibrationData
    		){
		int numEnabledSet=distortionCalibrationData.gIP[enabledImage].getSetNumber();
		distortionCalibrationData.gIS[numSet].setSetVector(distortionCalibrationData.gIS[numEnabledSet].getSetVector());
		System.out.println("getClosestSetChannelImage(): imageSet "+numSet+" set orientationEstimated=true, updated from imageSet "+numEnabledSet);
		distortionCalibrationData.gIS[numSet].orientationEstimated=(numSet!=numEnabledSet);
		double [] newVector=distortionCalibrationData.getParameters(enabledImage);
		for (int nChn=0;nChn<distortionCalibrationData.gIS[numSet].imageSet.length;nChn++)
			if (distortionCalibrationData.gIS[numSet].imageSet[nChn]!=null){ // will copy back to itself, OK
				int targetImage=distortionCalibrationData.gIS[numSet].imageSet[nChn].imgNumber;
				double [] oldVector=distortionCalibrationData.getParameters(targetImage);
				for (int j=0;j<oldVector.length;j++) if (parameterMask[j]) oldVector[j]=newVector[j];
				distortionCalibrationData.setParameters(oldVector, targetImage);
			}
    }

    // TODO: Add updating to all Stations depending on type of adjustment. Initially only teh same station as image will be updated
    // Not needed - for "super" unselected images are also updated
    /**
     * Update camera/subcamera parameters from the currently selected set of images
     * several images may have different values for the same parameter, in that case
     * these parameters will have the value of the last image
     */
    public void updateCameraParametersFromCalculated(
    		boolean allImages ){
    	int numSeries=allImages?(-1):this.fittingStrategy.currentSeriesNumber;
		boolean [] selectedImages=fittingStrategy.selectedImages(numSeries); // all enabled
		boolean [] selectedImagesDebug=null;
		boolean include_disabled = allImages;
		if (include_disabled) {
			Arrays.fill(selectedImages, true);
		}
		boolean debugThis=false;
		int maxDebugImages=10;
		if (this.debugLevel>0) System.out.println("updateCameraParametersFromCalculated("+allImages+")");
		if (this.debugLevel>2){
			int numSel=0;
			for (int i=0;i<selectedImages.length;i++) if (selectedImages[i]) numSel++;
			if (numSel<=maxDebugImages) debugThis=true;
			else {
				System.out.println ("Too many images ("+numSel+">"+ maxDebugImages +") to debug, skipping console println.");
				selectedImagesDebug=fittingStrategy.selectedImages(this.fittingStrategy.currentSeriesNumber); // all enabled

			}
		}
		for (int numImg=0;numImg<selectedImages.length; numImg++) if (selectedImages[numImg]){ // here only adjusted images should participate
			int subCam=fittingStrategy.distortionCalibrationData.getImageSubcamera(numImg);
//			double [] par=fittingStrategy.distortionCalibrationData.pars[numImg];
			double [] par=fittingStrategy.distortionCalibrationData.getParameters(numImg);
    		boolean [] update=new boolean[par.length];
    		for (int i=0;i<update.length;i++) update[i]=true;
    		int stationNumber=fittingStrategy.distortionCalibrationData.getImageStation(numImg);
    		// TODO: maybe determine - which parameters to be updated, not all - i.e. "super-common", or having the same value, etc.
    		// but all those intrinsic are required to match calibration files saved
			fittingStrategy.distortionCalibrationData.eyesisCameraParameters.setParametersVector(par, update, stationNumber, subCam);
			if (debugThis || ((selectedImagesDebug!=null) && selectedImagesDebug[numImg])){
				System.out.println ("Updating from image #"+numImg+" (subCam="+subCam+" stationNumber="+stationNumber+"):");
//getParameterName
				for (int i=0;i<par.length;i++){
					System.out.println(i+": "+fittingStrategy.distortionCalibrationData.getParameterName(i)+" = "+par[i]);
				}
			}
//			System.out.println(numImg+"[21]: "+fittingStrategy.distortionCalibrationData.getParameterName(21)+" = "+par[21]);
		}
		if (this.debugLevel>1) System.out.println("updateCameraParametersFromCalculated("+allImages+") for series="+numSeries);
		// Next line is not needed anymore (will harm as will set orientationEstimated for all unselected sets)
//		if (!allImages) fittingStrategy.distortionCalibrationData.updateSetOrientation(selectedImages); // only for selected images (not all enabled), OK
    }
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
	private static void startAndJoin(Thread[] threads)
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
}

