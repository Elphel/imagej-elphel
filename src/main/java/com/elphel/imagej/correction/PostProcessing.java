package com.elphel.imagej.correction;
/**
** -----------------------------------------------------------------------------**
** PostProcessing.java
**
** Experimenting with different postprocessing methods
**
**
** Copyright (C) 2012 Elphel, Inc.
**
** -----------------------------------------------------------------------------**
**
**  PostProcessing.java is free software: you can redistribute it and/or modify
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

//import java.awt.Rectangle;
import java.awt.Rectangle;
import java.io.File;
//import java.io.IOException;
//import java.util.ArrayList;
//import java.util.List;
import java.util.Properties;
//import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.calibration.CalibrationFileManagement;
import com.elphel.imagej.calibration.MultipleExtensionsFileFilter;
import com.elphel.imagej.calibration.PixelMapping;
import com.elphel.imagej.common.DoubleFHT;
import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.common.WindowTools;
import com.elphel.imagej.jp4.JP46_Reader_camera;

//import javax.swing.SwingUtilities;

//import loci.common.services.DependencyException;
//import loci.common.services.ServiceException;
//import loci.formats.FormatException;

//import Jama.Matrix;
//import PixelMapping.InterSensor;
import ij.IJ;
import ij.ImagePlus;
//import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
//import ij.io.FileSaver;
import ij.io.Opener;
//import ij.process.ColorProcessor;
//import ij.text.TextWindow;



public class PostProcessing {
	public int debugLevel=1;
	public PostProcessingParameters       postProcessingParameters= new PostProcessingParameters();
	public LinearFeaturesParameters       linearFeaturesParameters=new LinearFeaturesParameters();
	public DisparityCorrelationParameters disparityCorrelationParameters=new DisparityCorrelationParameters();
//	public PixelMapping pixelMapping=new PixelMapping(null,true,debugLevel);
	public PixelMapping                   pixelMapping=new PixelMapping(debugLevel); // just to access PixelMapping.InterSensor
//	public double [][][] aYCbCr=null;
	public int [] channels;
	public PixelMapping.InterSensor       interSensor=null;
	public void setProperties(String prefix,Properties properties){ // currently prefix ==""
    	this.postProcessingParameters.setProperties(prefix+"POST_PROCESSING_PARAMETERS.", properties);
    	this.linearFeaturesParameters.setProperties(prefix+"LINEAR_FEATURES_PARAMETERS.", properties);
    	this.disparityCorrelationParameters.setProperties(prefix+"DISPARITY_CORRELATION_PARAMETERS.", properties);
	}
	public void getProperties(String prefix,Properties properties){ // currently prefix ==""
    	this.postProcessingParameters.getProperties(prefix+"POST_PROCESSING_PARAMETERS.", properties);
    	this.linearFeaturesParameters.getProperties(prefix+"LINEAR_FEATURES_PARAMETERS.", properties);
    	this.disparityCorrelationParameters.getProperties(prefix+"DISPARITY_CORRELATION_PARAMETERS.", properties);
	}
	public void convertSourceFiles(int debugLevel){
   		String mapsPath=postProcessingParameters.equirectangularDirectory+
		Prefs.getFileSeparator()+postProcessingParameters.planeMapPrefix+postProcessingParameters.planeMapSuffix;

		Opener opener=new Opener();
		ImagePlus impMap=opener.openImage("", mapsPath);
    	if (impMap==null) {
    		String msg="Failed to read inter-sensor overlap map file "+mapsPath;
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
    	}
		if (debugLevel>0) System.out.println("Read "+mapsPath+" as an inter-sensor overlap map");
    	(new JP46_Reader_camera(false)).decodeProperiesFromInfo(impMap);
    	this.interSensor=pixelMapping.new InterSensor(impMap,debugLevel);
		System.out.println("Converting source files");
		this.debugLevel= debugLevel;
		this.pixelMapping.debugLevel= debugLevel;
		if (postProcessingParameters.sourcePaths==null) {
			System.out.println("No files selected");
			return;
		}
		this.channels=new int [postProcessingParameters.sourcePaths.length];
		interSensor.overlapImages=new double [4*this.interSensor.numOverlapChannels][];
		for (int i=0;i<interSensor.overlapImages.length;i++) interSensor.overlapImages[i]=null;
		double [][] aYCbCr=null;
		int width=0,height=0;
		for (int iFile=0;iFile<postProcessingParameters.sourcePaths.length;iFile++){
			String path=this.postProcessingParameters.sourcePaths[iFile];
			this.channels[iFile]=postProcessingParameters.getChannelFromSourceTiff(path);
			if (debugLevel>0) System.out.println(iFile+": processing "+path+", channel="+this.channels[iFile]);
			  ImagePlus imp_src=new ImagePlus(path);
			  aYCbCr=argb24ToAYCbCr(imp_src);
			  for (int i=0;i<aYCbCr.length;i++){
				  interSensor.overlapImages[4*this.channels[iFile]+i]=aYCbCr[i];
			  }
			  width=imp_src.getWidth();
			  height=imp_src.getHeight();

		}
		if (debugLevel>1) {
			String [] titles0={"Alpha","Y","Cb","Cr"};
			String [] titles=new String[interSensor.overlapImages.length];
			for (int i=0;i<titles.length;i++){
				titles[i]=titles0[i%titles0.length]+(i/titles0.length);
			}
			ShowDoubleFloatArrays.showArrays(
					interSensor.overlapImages,
					width,
					height,
					true,
					"aYCbCr",
					titles);
		}

	}
	//public double [][][][][]
	public PixelMapping.InterSensor.DisparityTiles getDisparityTiles(
			String externalTitle,
			double [][] externalData, // if not null - process instead of channel data
			int channelMask,
			int disparitySteps,
			Rectangle resultWindow,
			int threadsMax,
			boolean updateStatus,
			int debugLevel){
		  Runtime runtime = Runtime.getRuntime();
		  long 	  startTime=System.nanoTime();

//		this.disparityCorrelationParameters.xx
		int size=this.disparityCorrelationParameters.corrFFTSize;
		int numSensors=this.interSensor.channel.length;
		// cover the whole image with tiles, result window will just select active ones from all
		// this.interSensor.mapWidth
		int width= this.interSensor.mapWidth;
		int height=this.interSensor.mapHeight;
		int overlapStep=size/this.disparityCorrelationParameters.tileOverlapFraction;
		int numTilesX=width/overlapStep;  if (width>numTilesX*overlapStep) numTilesX++;
		int numTilesY=height/overlapStep; if (height>numTilesY*overlapStep) numTilesY++;
		int tileXMin=resultWindow.x/overlapStep;                        // tile including leftmost pixel
		int tileXMax=(resultWindow.x+resultWindow.width-1)/overlapStep; // tile including rightmost pixel
		int tileYMin=resultWindow.y/overlapStep;                        // tile including top pixel
		int tileYMax=(resultWindow.y+resultWindow.height-1)/overlapStep; // tile including bottom pixel
		// create list of image pairs
		int [][] otherImage = new int [numSensors][numSensors-1]; // here we need both {1,2} and {2,1}
		for (int i=0;i< numSensors;i++)  for (int j=0;j<(numSensors-1);j++) otherImage[i][j]=(j<i)?j:(j+1);
		// for each image pair create a list of displacements (in discrete tiles) needed for calculation of disparity sections over the specified range.
		// This list is needed to maintain a cache of direct FHT tiles deleting those that will not be needed anymore to reduce needed memory
		int [][][][] displacementTilesList = new int [numSensors][numSensors-1][][];
		// these show how far FFT tiles are needed . Or should we save more memory and group per each target image?
		int [] minDispalcementTileY=new int [numSensors];
		int [] maxDispalcementTileY=new int [numSensors];
		int [] minDispalcementTileX=new int [numSensors];
		int [] maxDispalcementTileX=new int [numSensors];
		for (int i=0;i<numSensors;i++){
			minDispalcementTileY[i]=0;
			maxDispalcementTileY[i]=0;
			minDispalcementTileX[i]=0;
			maxDispalcementTileX[i]=0;

		}
		int extraPixels=this.disparityCorrelationParameters.interpolationSize/2;
		int [][][][][][] correlationMap=new int [numSensors][][][][][]; // [image number][other image index][pixelY][pixelX]{{tiledX,tilledY,pixelIndex}, ...}
//		int [][][][] bilinearPixelIndex=new int [numSensors][][][];     // [image number][other image index]{{pixelX,pixelY},...}

		for (int i=0;i< numSensors;i++)  {
			correlationMap[i]=    new int [otherImage[i].length][][][][];
//			bilinearPixelIndex[i]=new int [otherImage[i].length][][];
//			int [][] tileDirs8={{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}};
//			int [][] tileDirs9={{0,0},{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}};
			int [][] tileDirs={{0,0},{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}}; // includes center
			for (int j=0;j<otherImage[i].length;j++) {
				int other=otherImage[i][j]; // number of image in a pair
				double maxDX=this.disparityCorrelationParameters.disparityRange*
				(this.disparityCorrelationParameters.disparityScales[other][0]-this.disparityCorrelationParameters.disparityScales[i][0]);
				double maxDY=this.disparityCorrelationParameters.disparityRange*
				(this.disparityCorrelationParameters.disparityScales[other][1]-this.disparityCorrelationParameters.disparityScales[i][1]);
				int maxDist = (int) Math.ceil(Math.sqrt(maxDX*maxDX+maxDY*maxDY)/overlapStep);
				maxDist++; // needed more?
				if (debugLevel>1){
					System.out.println("this="+i+" other="+other+" maxDX="+maxDX+" maxDY="+maxDY+" maxDist="+maxDist);
				}

				boolean [][] neededTilesMap=new boolean[2*maxDist+1][2*maxDist+1]; // zero displacement in the center, [maxDist][maxDist]
				for (int idy=0;idy<neededTilesMap.length;idy++) for (int idx=0;idx<neededTilesMap[0].length;idx++) neededTilesMap[idy][idx]=false;
				for (int di=0;di<=this.disparityCorrelationParameters.disparitySteps;di++){
					double dy=(maxDY*di)/this.disparityCorrelationParameters.disparitySteps;
					double dx=(maxDX*di)/this.disparityCorrelationParameters.disparitySteps;
					int dxtMin=(int) Math.floor((dx-extraPixels)/overlapStep);
					int dxtMax=(int) Math.ceil ((dx+extraPixels)/overlapStep);
					int dytMin=(int) Math.floor((dy-extraPixels)/overlapStep);
					int dytMax=(int) Math.ceil ((dy+extraPixels)/overlapStep);
					if (dxtMin==dxtMax){
						dxtMin--;
						dxtMax++;
					}
					if (dytMin==dytMax){
						dytMin--;
						dytMax++;
					}
					for (int dyt=dytMin;dyt<=dytMax;dyt++) for (int dxt=dxtMin;dxt<=dxtMax;dxt++) neededTilesMap[maxDist+dyt][maxDist+dxt]=true;
				}
				int numNeeded=0;
				for (int idy=0;idy<neededTilesMap.length;idy++) for (int idx=0;idx<neededTilesMap[0].length;idx++) if (neededTilesMap[idy][idx]) numNeeded++;
				int indexNeeded=0;
				displacementTilesList[i][j]=new int [numNeeded][2];
				for (int idy=0;idy<neededTilesMap.length;idy++) for (int idx=0;idx<neededTilesMap[0].length;idx++) if (neededTilesMap[idy][idx]){
					int dx=idx-maxDist;
					int dy=idy-maxDist;
					displacementTilesList[i][j][indexNeeded][0]=dx;
					displacementTilesList[i][j][indexNeeded][1]=dy;
					if (dy<minDispalcementTileY[other]) minDispalcementTileY[other]=dy;
					if (dy>maxDispalcementTileY[other]) maxDispalcementTileY[other]=dy;
					if (dx<minDispalcementTileX[other]) minDispalcementTileX[other]=dx;
					if (dx>maxDispalcementTileX[other]) maxDispalcementTileX[other]=dx;
					indexNeeded++;
				}
				 if (debugLevel>2){
				for (int idy=0;idy<neededTilesMap.length;idy++) {
					for (int idx=0;idx<neededTilesMap[0].length;idx++){
						if ((idy==maxDist) && (idx==maxDist)) System.out.print("* ");
						else System.out.print(neededTilesMap[idy][idx]?"+ ":"  ");
					}
					System.out.println();
				}
				 }
				correlationMap[i][j]=new int [neededTilesMap.length*overlapStep][neededTilesMap[0].length*overlapStep][][];
				for (int py=0;py<correlationMap[i][j].length;py++) for (int px=0;px<correlationMap[i][j][py].length;px++) correlationMap[i][j][py][px]=null;
				for (int tileY=0;tileY<neededTilesMap.length;tileY++){
					for (int tileX=0;tileX<neededTilesMap[tileY].length;tileX++) if (neededTilesMap[tileY][tileX]){
						int numNeib=0;
						for (int iDir=0;iDir<tileDirs.length;iDir++){
							int tileXn=tileX+tileDirs[iDir][0];
							int tileYn=tileY+tileDirs[iDir][1];
							if ((tileXn>=0) && (tileXn<neededTilesMap[0].length) && (tileYn>=0) && (tileYn<neededTilesMap.length) && neededTilesMap[tileYn][tileXn]) numNeib++;
						}
						for (int pdY=0;pdY<overlapStep;pdY++) {
							int pY=tileY*overlapStep+pdY;
							for (int pdX=0;pdX<overlapStep;pdX++) {
								int pX=tileX*overlapStep+pdX;
								correlationMap[i][j][pY][pX]=new int [numNeib][3];
								int indexNeib=0;
								int marginCenter=(size-overlapStep)/2;
								for (int iDir=0;iDir<tileDirs.length;iDir++){
									int tileXn=tileX+tileDirs[iDir][0];
									int tileYn=tileY+tileDirs[iDir][1];
									if ((tileXn>=0) && (tileXn<neededTilesMap[0].length) && (tileYn>=0) && (tileYn<neededTilesMap.length) && neededTilesMap[tileYn][tileXn]) {
										correlationMap[i][j][pY][pX][indexNeib][0]=tileDirs[iDir][0]; // tileDX
										correlationMap[i][j][pY][pX][indexNeib][0]=tileDirs[iDir][1]; // tileDY
										correlationMap[i][j][pY][pX][indexNeib][0]=marginCenter*(size+1)-overlapStep*(tileDirs[iDir][0]+tileDirs[iDir][1]*size); // tilePixelIndex
										indexNeib++;
									}
								}
							}
						}
					}
				}

			}
		}
		// debug show correlationMap, 	bilinearPixelIndex
		if (debugLevel>1){
			System.out.println("getDisparityTiles(): needed FHT tiles from each image:");
			for (int i=0;i<numSensors;i++){
				System.out.println("Image "+i+
						" minDispalcementTileY="+minDispalcementTileY[i]+
						" maxDispalcementTileY="+maxDispalcementTileY[i]+
						" minDispalcementTileX="+minDispalcementTileX[i]+
						" maxDispalcementTileX="+maxDispalcementTileX[i] );
			}
		}
       	double [][][][][] disparityTilesDouble= interSensor.calculateDisparityTiles(
    			externalTitle,
    			externalData, // if not null - process instead of channel data
    			channelMask,
    			otherImage,
    			displacementTilesList,
    			disparitySteps,
    			this.disparityCorrelationParameters.disparityRange,
    			this.disparityCorrelationParameters.disparityScales,
    			this.disparityCorrelationParameters.interpolationSize,
    			this.disparityCorrelationParameters.interpolationUpSample,
    			correlationMap, // =new int [numSensors][][][][][]; // [image number][other image index][pixelY][pixelX]{{tiledX,tilledY,pixelIndex}, ...}
 //   			bilinearPixelIndex, //=new int [numSensors][][][];     // [image number][other image index]{{pixelX,pixelY},...}
    			minDispalcementTileX,
    			maxDispalcementTileX,
    			minDispalcementTileY,
    			maxDispalcementTileY,
    			this.disparityCorrelationParameters.corrPhaseFraction,
    			this.disparityCorrelationParameters.correlationHighPassSigma,
    			this.disparityCorrelationParameters.correlationLowPassSigma,
    			this.disparityCorrelationParameters.corrCbWeight,
    			this.disparityCorrelationParameters.corrCrWeight,
    			this.disparityCorrelationParameters.subpixAMax,
    			this.disparityCorrelationParameters.subpixRMax,
    			this.disparityCorrelationParameters.subpixNMax,
    			this.disparityCorrelationParameters.corrFFTSize,
				overlapStep,
				numTilesX,
				numTilesY,
				tileXMin,
				tileXMax,
				tileYMin,
				tileYMax,
    			threadsMax,
    			updateStatus,
				debugLevel);

		  runtime.gc();
		  if (debugLevel>0) System.out.println("getDisparityTiles() done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+
				  " seconds  Free memory="+IJ.d2s(runtime.freeMemory()/(1024.0*1024.0*1024.0),3)+" GB (of "+
				  IJ.d2s(runtime.totalMemory()/(1024.0*1024.0*1024.0),3)+" GB), used "+
				  IJ.d2s((runtime.totalMemory()-runtime.freeMemory())/(1024.0*1024.0*1024.0),3)+" GB");

		  // 			if (this.disparityCorrelationParameters.showCorrelationToImage){
		  String title;
		  title=(externalData==null)?"Encoded_disparity":((externalTitle==null)?"external":externalTitle);
		  title+="_"+this.disparityCorrelationParameters.corrFFTSize+"-"+overlapStep+".corr-tiff";

		  PixelMapping.InterSensor.DisparityTiles dispatityTiles=interSensor.new DisparityTiles (
				  this.disparityCorrelationParameters.disparityScales,
				  this.disparityCorrelationParameters.corrFFTSize,
				  overlapStep, // to properties
				  this.disparityCorrelationParameters.disparityRange/disparitySteps, //disparity increment (in pix) per array element
				  title,
				  disparityTilesDouble,
				  this.disparityCorrelationParameters.showCorrelationToImage);
		  return dispatityTiles; //Double; // just temporarily
	}
	/*
	public void  setupZMap(
			Rectangle zMapWOI,
			double zMapMinFirst,
			double zMapMinAbsolute,
			double zMapMinRelative,
			double zMapMergeMax,
			int zMapMaxNumber,
			int zMapOverlap,
    		int threadsMax, //final int          threadsMax,  // maximal number of threads to launch
    		boolean updateStatus,
    		int debugLevel){
		 //int debugLevel
	}
*/

	public void matchTest(
			double blurVarianceSigma,
			int refineTilePeriod,
			double refinePhaseCoeff,
			double refineHighPassSigma,
			double refineLowPassSigma,
			double refineCorrMaxDistance,
			double refineCorrThreshold,
			int refineSubPixel,
			double [][] externalData,
			double disparity, //POST_PROCESSING.disparityCorrelationParameters.corrShift,
			double xc, //POST_PROCESSING.disparityCorrelationParameters.corrXC,
			double yc, //POST_PROCESSING.disparityCorrelationParameters.corrYC,
			int channelMask,
			int firstImageNumber,
			int secondImageNumber,
			boolean doubleSizeOutput,
			int size0,
			int debugLevel){
		DoubleGaussianBlur gb=new DoubleGaussianBlur();
		double dX=disparity*(this.disparityCorrelationParameters.disparityScales[secondImageNumber][0]-
				this.disparityCorrelationParameters.disparityScales[firstImageNumber][0]);
		double dY=disparity*(this.disparityCorrelationParameters.disparityScales[secondImageNumber][1]-
				this.disparityCorrelationParameters.disparityScales[firstImageNumber][1]);
		double dXYlen=Math.sqrt(dX*dX+dY*dY);
		double [] udXY={dX/dXYlen,dY/dXYlen};
		double [][] slices=interSensor.getShiftedSlices(
				(int) Math.round(xc),
				(int) Math.round(yc),
				size0,
				firstImageNumber, //int first,
				secondImageNumber, //int second,
				dX, //double dxA,
				dY, //double dyA,
				channelMask,
				externalData,
				doubleSizeOutput,
				null, // double [] window, // should be 4*size*size long;
				null, //DoubleFHT doubleFHT, // to reuse tables
				debugLevel);
//
		int [] img= {firstImageNumber,secondImageNumber};
		int numImg=img.length;
		String [] channelNames={"Y","Cb","Cr","Ext"};
	    String [] resultNames={"corr","bdiff","diff","inv_diff2","var","bdiff/var","diff-dX","diff-dY"};
		int numLayers=0;
		for (int chm=channelMask; chm!=0; chm>>=1) if ((chm & 1)!=0) numLayers++;
		if (externalData!=null) numLayers++;
//		int [] img= {first,second};
		int [] channels=new int [numLayers];
		int c=0;
		int l=0;
		for (int chm=channelMask; chm!=0; chm>>=1) {
			if ((chm & 1)!=0) channels[l++]=c;
			c++;
		}
		if (externalData!=null) channels[l]=3;
		int length=slices[0].length;
		int size=(int) Math.sqrt(length);
		String [] titlesSlices=new String[2*numLayers];
		String titleSlices="SH";
		for (int i=0;i<numImg;i++) titleSlices+="-"+img[i];
		double dd=Math.sqrt(dX*dX+dY*dY);
		titleSlices+="_D"+IJ.d2s(dd,2)+"_X"+IJ.d2s(dX,2)+"_Y"+IJ.d2s(dY,2);
		for (int i=0;i<titlesSlices.length;i++){
			titlesSlices[i]=img[i%img.length]+channelNames[channels[i/img.length]];
		}
		if (debugLevel>1){
			ShowDoubleFloatArrays.showArrays(
					slices,
					size,
					size,
					true,
					titleSlices,
					titlesSlices);
		}
		double [][] results=new double[(slices.length/numImg)*resultNames.length][slices[0].length];
		String [] titles= new String[results.length];
		int [] dirs1={1,size+1,size,size-1,-1,-size-1,-size,-size+1,0};
		int [][] dirsDxy={
				{1,   -1,   size+1, size-1,-size+1,-size-1},
				{size,-size,size+1,-size+1, size-1,-size-1}};
		double [] weightsDxy={0.25,-0.25,0.125,-0.125,0.125,-0.125};
		int refineTilePerRow=size/2/refineTilePeriod;
		double [] refineWindow=new double [4*refineTilePeriod*refineTilePeriod];
		refineWindow[0]=Double.NaN; // mark that it needs to be initialized
//		double [] udXY={dX/dXYlen,dY/dXYlen};
		double [][] disparityArrays=new double[refineTilePerRow*refineTilePerRow][];
		DoubleFHT refineFHT= new DoubleFHT();

		for (int chn=0;chn<slices.length/numImg;chn++){
			double [][] pair ={slices[numImg*chn],slices[numImg*chn+1]};
			int indexCorr=    chn*resultNames.length+0;
			int indexBDiff=   chn*resultNames.length+1;
			int indexDiff=    chn*resultNames.length+2;
			int indexInvDiff= chn*resultNames.length+3;
			int indexVar=     chn*resultNames.length+4;
			int indexDiffVar= chn*resultNames.length+5;
			int indexDiffDx=  chn*resultNames.length+6;
			int indexDiffDy=  chn*resultNames.length+7;

			double [] refineSmooth=interSensor.refineCorrelation(
					refineFHT, // DoubleFHT doubleFHT, // will generate if null, can be used to reuse initialization
					pair, //double [][] data,   // difines size
					refineWindow, //double [] window, // defines tile size
					udXY, // double [] dXY,    // unity vector defines disparity direction
					refinePhaseCoeff, //phaseCoeff, // phase correlation fraction (0 - normal, 1.0 - pure phase)
					refineHighPassSigma,
					refineLowPassSigma,
					disparityArrays, // double [][] disparityArrays, // top level should be initialized
					refineCorrMaxDistance, //double corrMaxDist, // distance from rhe center to look for the correlation maximum for disparity correction
					refineCorrThreshold, //double corrThreshold, // relative correlation contrast to correct disparity
					refineSubPixel, //int subPixel,       // subdivide pixels when updating disparity
	    			true, //boolean generateSmooth,
	    			debugLevel
	    			);
			for (int nTyp=0;nTyp<resultNames.length;nTyp++){
				int index=chn*resultNames.length+nTyp;
				titles[index]=resultNames[nTyp]+"-"+channelNames[channels[chn]];
				switch (nTyp){
				case 0:	// corr
					for (int i=0;i<length;i++) results[indexCorr][i]=refineSmooth[i];
					break;

				case 1:	// diff
				//	for (int i=0;i<length;i++) results[indexDiff][i]=pair[1][i]-pair[0][i];
				// bestDiff
					for (int i=0;i<length;i++) {
						boolean less=false;
						boolean more=false;
						double bdiff=0;
						for (int d=0;d<dirs1.length;d++){
							int i1=(i+dirs1[d]+length)%length;
							double df=pair[1][i1]-pair[0][i1];
							less |= df<=0;
							more |= df>=0;
							df=Math.abs(df);
							if ((d==0) || (df<d)) bdiff=df;
							if (more && less){
								bdiff=0;
								break;
							}
						}
						results[indexBDiff][i]=bdiff;
					}
					break;
				case 2:	// diff
						for (int i=0;i<length;i++) results[indexDiff][i]=pair[1][i]-pair[0][i];
						break;
				case 3:	// 1/diff2
//					for (int i=0;i<length;i++) results[index][i]=1.0/(results[indexDiff][i]*results[indexDiff][i]);
					for (int i=0;i<length;i++) results[indexInvDiff][i]=Math.abs(results[indexDiff][i]);
					break;
				case 4:	// variance
					for (int i=0;i<length;i++) {
						double [] S1={0.0,0.0};
						double [] S2=S1.clone();
						double  v2=0;
						for (int n=0;n<2;n++) {
							for (int d=0;d<dirs1.length;d++){
								int i1=(i+dirs1[d]+length)%length;
								S1[n]+=pair[n][i1];
								S2[n]+=(pair[n][i1]*pair[n][i1]);
							}
							v2+=(S2[n]-(S1[n]*S1[n])/dirs1.length)/dirs1.length;
						}
						results[indexVar][i]=Math.sqrt(v2);
					}
					if (blurVarianceSigma>0.0){
			        	   gb.blurDouble(
			        			   results[indexVar],
			        			   size,
			        			   size,
			        			   blurVarianceSigma,
			        			   blurVarianceSigma,
			        			   0.01);
					}
					break;
				case 5:	// variance/|diff| (squared?)
//					for (int i=0;i<length;i++) results[indexDiffVar][i]=results[indexVar][i]/Math.abs(results[indexDiff][i]); // use square?
//					for (int i=0;i<length;i++) results[indexDiffVar][i]=Math.abs(results[indexDiff][i])/results[indexVar][i]; // use square?
					for (int i=0;i<length;i++) results[indexDiffVar][i]=results[indexBDiff][i]/results[indexVar][i]; //
					break;
				case 6:	// diff/(d/dx)
					for (int i=0;i<length;i++) {
						double ddx=0;
						for (int n=0;n<2;n++) {
							for (int d=0;d<dirsDxy[0].length;d++){
								int i1=(i+dirsDxy[0][d]+length)%length;
								ddx+=weightsDxy[d]*pair[n][i1];
							}
						}
						ddx*=0.5; // was a sum for 2 images, need average
						results[indexDiffDx][i]=results[indexDiff][i]/ddx;
					}
					break;
				case 7:	// diff/(d/dy)
					for (int i=0;i<length;i++) {
						double ddy=0;
						for (int n=0;n<2;n++) {
							for (int d=0;d<dirsDxy[1].length;d++){
								int i1=(i+dirsDxy[1][d]+length)%length;
								ddy+=weightsDxy[d]*pair[n][i1];
							}
						}
						ddy*=0.5; // was a sum for 2 images, need average
						results[indexDiffDy][i]=results[indexDiff][i]/ddy;
					}
					break;
				}

			}

		}
		if (debugLevel>0){
			String title="DSH";
			for (int i=0;i<numImg;i++) titleSlices+="-"+img[i];
			titleSlices+="_D"+IJ.d2s(dd,2)+"_X"+IJ.d2s(dX,2)+"_Y"+IJ.d2s(dY,2);
			ShowDoubleFloatArrays.showArrays(
					results,
					size,
					size,
					true,
					title,
					titles);
		}
	}

	public double [][] correlationTest(
			double disparity,
			double xc,
			double yc,
			int threadsMax,
			int debugLevel){
		if (this.disparityCorrelationParameters.disparityScales==null){
			System.out.println ("disparityCorrelationParameters.disparityScales==null");
			return null;
		}
		int [][] iCenterXY=new int [this.disparityCorrelationParameters.disparityScales.length][];

//    		    int [] zeroShift={0,0};
//	    int [][] iShifts=new int [iShifts0.length][];
//	    for (int i=0;i<iShifts.length;i++) iShifts[i]= (autoCorrelation || (iShifts0[i]==null))?zeroShift.clone():iShifts0[i].clone();


		// Calculate rounded integer shifts per channel


		for (int i=0;i<iCenterXY.length;i++){
			if (this.disparityCorrelationParameters.disparityScales[i]!=null){
				iCenterXY[i]=new int[2];
				iCenterXY[i][0]= (int) Math.round(xc+disparity*this.disparityCorrelationParameters.disparityScales[i][0]);
				iCenterXY[i][1]= (int) Math.round(yc+disparity*this.disparityCorrelationParameters.disparityScales[i][1]);
			} else iCenterXY[i]=null;
		}

//		int numSensors=this.interSensor.channel.length;
//		int numLayers=this.interSensor.overlapImages.length/numSensors;
		double [][] corr=interSensor.correlate(
    			iCenterXY, // for each image - centerX, centerY
    			this.disparityCorrelationParameters.autocorrelation,
    			this.disparityCorrelationParameters.corrFFTSize,
    			this.disparityCorrelationParameters.corrPhaseFraction,
    			this.disparityCorrelationParameters.corrCbWeight,
    			this.disparityCorrelationParameters.corrCrWeight,
    			this.disparityCorrelationParameters.correlationHighPassSigma,
    			this.disparityCorrelationParameters.correlationLowPassSigma,
    			this.disparityCorrelationParameters.noiseNormalizationSignaY,
    			this.disparityCorrelationParameters.noiseNormalizationSignaCbCr,
    			this.disparityCorrelationParameters.contrastThreshold,
    			this.disparityCorrelationParameters.useBinaryAlpha,
				threadsMax,
				debugLevel);

//		int numPairs=corr.length/numLayers;
	  return corr;



	}

	public void linearFeatures(
			int threadsMax,
			int debugLevel){
		if (pixelMapping==null) {
			String msg="pixelMapping is not initialized";
			System.out.println(msg);
			IJ.showMessage(msg);
		}
		if ((interSensor==null) || (interSensor.overlapImages==null)) {
			String msg="Inter-sensor data is not initialized (use \"Source PP\" command";
			System.out.println(msg);
			IJ.showMessage(msg);
		}
//		boolean multiplyByCellUsage=    gd.getNextBoolean();
//		double cellUsageShift=          gd.getNextNumber();
//		if (!multiplyByCellUsage) cellUsageShift=Double.NaN;

		if (interSensor.linearFeatures==null) {
			interSensor.linearFeatures=new double [this.interSensor.numOverlapChannels][][][];
			for (int i=0;i<interSensor.linearFeatures.length;i++) interSensor.linearFeatures[i]=null;
		}
		long startTime = System.nanoTime();
		//this.interSensor.numOverlapChannels
		for (int numImage=0;numImage<this.interSensor.numOverlapChannels;numImage++) if (interSensor.overlapImages[4*numImage]!=null) {
			boolean reExtract=this.linearFeaturesParameters.extractFeatures || (interSensor.linearFeatures[numImage]==null);
//			/this.overlapImages[i]
			if (debugLevel>1){
				for (int i=0;i<interSensor.overlapImages.length;i++){
					System.out.println ("overlapImages["+i+"] is "+((interSensor.overlapImages[i]==null)?"null":"not null"));

				}
			}
			if (reExtract) interSensor.linearFeatures[numImage]= interSensor.detectLinearFeatures(
					numImage, //(side>0),
					this.linearFeaturesParameters.corrFFTSize,
					this.linearFeaturesParameters.overlapFraction, // default 8
					this.linearFeaturesParameters.alphaThreshold, // minFracArea
					this.linearFeaturesParameters.useBinaryAlpha,
					this.linearFeaturesParameters.correlationHighPassSigma,
					this.linearFeaturesParameters.correlationLowPassSigma,
					this.linearFeaturesParameters.phaseIntegrationWidth,
					this.linearFeaturesParameters.resultHighPass, //cutting "ribbon" near zero frequency - influences length of the detected line (used only for strength calculation)
					this.linearFeaturesParameters.dispertionCost,
					this.linearFeaturesParameters.featureFilter,
					this.linearFeaturesParameters.minRMSFrac,
					this.linearFeaturesParameters.minAbs,
					this.linearFeaturesParameters.maxPhaseMismatch,
					this.linearFeaturesParameters.calculateStrength,
					this.linearFeaturesParameters.debugRow,
					this.linearFeaturesParameters.debugColumn,
	        		threadsMax,
					debugLevel);
			if (debugLevel>1) System.out.println("Linear Features for "+numImage+" finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
			// startTime
		}
		// startTime
//		double [][] null2 = {null,null};
//		double [][] null4 = {null,null,null,null};

		double [][] linearFeatures = new double [((this.linearFeaturesParameters.displayUsageScale>0)?2:1)*this.interSensor.numOverlapChannels][];

		for (int numImage=0;numImage<this.interSensor.numOverlapChannels;numImage++) if (interSensor.overlapImages[4*numImage]!=null){
			linearFeatures[numImage<<((this.linearFeaturesParameters.displayUsageScale>0)?1:0)]=interSensor.reconstructImageFeature(
					this.linearFeaturesParameters.corrFFTSize,
					this.linearFeaturesParameters.overlapFraction, // default 8
					this.interSensor.linearFeatures[numImage],
					this.linearFeaturesParameters.ignorePhase,
					this.linearFeaturesParameters.preserveDC,
					this.linearFeaturesParameters.strengthMode, // 0.0 - absolute, 1.0 - relative, 0.5 - "balanced", "-1" - use phaseStrength instead
					this.linearFeaturesParameters.phaseIntegrationWidth, // use the same as during extraction?
					this.linearFeaturesParameters.resultHighPass, //cutting "ribbon" near zero frequency - influences length of the detected line
					this.linearFeaturesParameters.debugRow,
					this.linearFeaturesParameters.debugColumn,
					threadsMax,
					debugLevel);
				if (this.linearFeaturesParameters.displayUsageScale>0){
					linearFeatures[(numImage<<1)+1]=interSensor.displayUsage(
							interSensor.linearFeatures[numImage],
							this.linearFeaturesParameters.corrFFTSize,
							this.linearFeaturesParameters.overlapFraction, // default 8
							this.linearFeaturesParameters.displayUsageScale);
				}
				if (debugLevel>1) System.out.println("Linear Features simulated for "+(numImage)+" finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));

		}
//TODO: just show linearFeatures array and measure time
//		String [] titlesLR={"Left","Right"};
//		String [] titlesLRU={"Left","Usage Left","Right","Usage Right"};
//		String [] titles=(this.linearFeaturesParameters.displayUsageScale>0)?titlesLRU:titlesLR;
		//			double displayUsageScale=       gd.getNextNumber();

		String [] titles = new String[((this.linearFeaturesParameters.displayUsageScale>0)?2:1)*this.interSensor.numOverlapChannels];
		if (this.linearFeaturesParameters.displayUsageScale>0){
			for( int i=0;i<this.interSensor.numOverlapChannels;i++){
				titles[2*i]="Chn "+i;
				titles[2*i+1]="Usage "+i;
			}
		} else {
			for( int i=0;i<this.interSensor.numOverlapChannels;i++){
				titles[i]="Chn "+i;
			}
		}

		ShowDoubleFloatArrays.showArrays(
				linearFeatures,
				interSensor.mapWidth,
				interSensor.mapHeight,
				true,
				"Linear features",
				titles);
//		double [][] linearFeatures = new double [((displayUsageScale>0)?2:1)*this.interSensor.numOverlapChannels][];
//		double [][][][] filteredFeatures={null,null};
//		double [][] filteredFeaturesImage =(displayUsageScale>0)?null4:null2;;
		double [][][][] filteredFeatures=new double [this.interSensor.numOverlapChannels][][][];
		double [][] filteredFeaturesImage =  new double [((this.linearFeaturesParameters.displayUsageScale>0)?2:1)*this.interSensor.numOverlapChannels][];
		for (int i=0;i<filteredFeatures.length;i++)filteredFeatures[i]=null;
		for (int i=0;i<filteredFeaturesImage.length;i++)filteredFeaturesImage[i]=null;
//		for (int side=0;side<2;side++) if ((sides & (1<<side))!=0) {
			for (int numImage=0;numImage<this.interSensor.numOverlapChannels;numImage++) if (interSensor.overlapImages[4*numImage]!=null){
			filteredFeatures[numImage]=interSensor.mergeLinearFeatures(
					this.linearFeaturesParameters.corrFFTSize,
					this.linearFeaturesParameters.overlapFraction, // default 8
	    			this.interSensor.linearFeatures[numImage],
	    			this.linearFeaturesParameters.strengthMode, // 0.0 - absolute, 1.0 - relative, 0.5 - "balanced", "-1" - use phaseStrength instead
	    			this.linearFeaturesParameters.directionTolerance,
	    			this.linearFeaturesParameters.normalDistanceTolerance,
	    			this.linearFeaturesParameters.tangentialDistanceTolerance,
	    			this.linearFeaturesParameters.hostsTolerance,
	    			this.linearFeaturesParameters.directionFracSigma, // make a fixed fraction of directionTolerance?
	    			this.linearFeaturesParameters.normalDistanceFracSigma, // make a fixed fraction of distanceTolerance?
	    			this.linearFeaturesParameters.tangentialDistanceFracSigma, // make a fixed fraction of distanceTolerance?
	    			this.linearFeaturesParameters.minMerged,
	    			this.linearFeaturesParameters.multiplyByCellUsage?this.linearFeaturesParameters.cellUsageShift:Double.NaN,
	    			this.linearFeaturesParameters.scaleDistances,
	    			this.linearFeaturesParameters.swapTangentialTolerance,
	    			this.linearFeaturesParameters.swapSearchRange,
	    			this.linearFeaturesParameters.debugRow,
	    			this.linearFeaturesParameters.debugColumn,
	    			debugLevel);
			if (interSensor.linearFeatures[numImage]!=null){
				filteredFeaturesImage[numImage<<((this.linearFeaturesParameters.displayUsageScale>0)?1:0)]=interSensor.reconstructImageFeature(
						this.linearFeaturesParameters.corrFFTSize,
						this.linearFeaturesParameters.overlapFraction, // default 8
						filteredFeatures[numImage],
						this.linearFeaturesParameters.ignorePhase,
						this.linearFeaturesParameters.preserveDC,
						this.linearFeaturesParameters.strengthMode, // 0.0 - absolute, 1.0 - relative, 0.5 - "balanced", "-1" - use phaseStrength instead
						this.linearFeaturesParameters.phaseIntegrationWidth, // use the same as during extraction?
						this.linearFeaturesParameters.resultHighPass, //cutting "ribbon" near zero frequency - influences length of the detected line
						this.linearFeaturesParameters.debugRow,
						this.linearFeaturesParameters.debugColumn,
						threadsMax,
						debugLevel);
				if (this.linearFeaturesParameters.displayUsageScale>0){
					filteredFeaturesImage[(numImage<<1)+1]=interSensor.displayUsage(
							filteredFeatures[numImage],
							this.linearFeaturesParameters.corrFFTSize,
							this.linearFeaturesParameters.overlapFraction, // default 8
							this.linearFeaturesParameters.displayUsageScale);
				}

				if (debugLevel>1) System.out.println("Filtered linear Features simulated for "+(numImage)+" finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
			}
		}
//TODO: just show linearFeatures array and measure time
			ShowDoubleFloatArrays.showArrays(
				filteredFeaturesImage,
				interSensor.mapWidth,
				interSensor.mapHeight,
				true,
				"Filtered Linear features",
				titles);
		return;


	}

	public double [][] argb24ToAYCbCr(
			ImagePlus imp){
		int [] iPixels = (int []) imp.getProcessor().getPixels();
		double KR=0.299, KB=0.114;
		double KG=1.0-KR-KB;
		double K255=1.0/255;
		double KPB=0.5/(1-KB);
		double KPR=0.5/(1-KR);
		double [][] aYCbCr=new double [4][iPixels.length];
		for (int i=0;i<iPixels.length;i++){
			aYCbCr[0][i]=K255*((iPixels[i]>>24) & 0xff); // alpha channel
			double r=K255*((iPixels[i]>>16) & 0xff);
			double g=K255*((iPixels[i]>> 8) & 0xff);
			double b=K255*( iPixels[i]      & 0xff);
			double y=KR*r+KG*g+KB*b;
			aYCbCr[1][i]=y;
			aYCbCr[2][i]=KPB*(b-y);
			aYCbCr[3][i]=KPR*(r-y);
		}
		return aYCbCr;
	}



	public class PostProcessingParameters{
  		public boolean saveSettings =          true;
		public String [] sourcePaths={};
		public String sourceDirectory="";
		public String sourcePrefix="";
		public String sourceSuffix=".tiff"; //".jp4"
    	public String equirectangularDirectory="";
    	public String planeMapPrefix="";
    	public String planeMapSuffix=".plane-proj-tiff";
    	public String resultsDirectory="";

    	public void setProperties(String prefix,Properties properties){
  			properties.setProperty(prefix+"saveSettings",this.saveSettings+"");
    		properties.setProperty(prefix+"sourceDirectory",this.sourceDirectory);
    		properties.setProperty(prefix+"sourcePrefix",this.sourcePrefix);
    		properties.setProperty(prefix+"sourceSuffix",this.sourceSuffix);
    		properties.setProperty(prefix+"equirectangularDirectory",this.equirectangularDirectory);
    		properties.setProperty(prefix+"planeMapPrefix",this.planeMapPrefix+"");
    		properties.setProperty(prefix+"planeMapSuffix",this.planeMapSuffix+"");

    		properties.setProperty(prefix+"resultsDirectory",this.resultsDirectory);
    		if (this.sourcePaths!=null) {
        		properties.setProperty(prefix+"sourcePaths",this.sourcePaths.length+"");
        		for (int i=0;i<this.sourcePaths.length;i++){
        			properties.setProperty(prefix+"sourcePath"+i,this.sourcePaths[i]);
        		}
    		}

    	}

    	public void getProperties(String prefix,Properties properties){
  		    if (properties.getProperty(prefix+"saveSettings")!=null) this.saveSettings=Boolean.parseBoolean(properties.getProperty(prefix+"saveSettings"));
			if (properties.getProperty(prefix+"sourceDirectory")!=      null) this.sourceDirectory=properties.getProperty(prefix+"sourceDirectory");
			if (properties.getProperty(prefix+"sourcePrefix")!=         null) this.sourcePrefix=properties.getProperty(prefix+"sourcePrefix");
			if (properties.getProperty(prefix+"sourceSuffix")!=         null) this.sourceSuffix=properties.getProperty(prefix+"sourceSuffix");
			if (properties.getProperty(prefix+"equirectangularDirectory")!=null) this.equirectangularDirectory=properties.getProperty(prefix+"equirectangularDirectory");
			if (properties.getProperty(prefix+"planeMapPrefix")!=null) this.planeMapPrefix=properties.getProperty(prefix+"planeMapPrefix");
			if (properties.getProperty(prefix+"planeMapSuffix")!=null) this.planeMapSuffix=properties.getProperty(prefix+"planeMapSuffix");

			if (properties.getProperty(prefix+"resultsDirectory")!=     null) this.resultsDirectory=properties.getProperty(prefix+"resultsDirectory");
			if (properties.getProperty(prefix+"sourcePaths")!=   null){
				int numFiles=Integer.parseInt(properties.getProperty(prefix+"sourcePaths"));
				this.sourcePaths=new String[numFiles];
				for (int i=0;i<this.sourcePaths.length;i++){
					this.sourcePaths[i]=properties.getProperty(prefix+"sourcePath"+i);
        		}
			}
    	}

    	public boolean showDialog(String title) {
    		GenericDialog gd = new GenericDialog(title);
    		gd.addCheckbox ("Save current settings with results",  this.saveSettings);
    		gd.addStringField("Source files directory",            this.sourceDirectory, 60);
    		gd.addCheckbox("Select source directory",              false);
    		gd.addStringField("Pixel maps directory",              this.equirectangularDirectory, 60);
    		gd.addCheckbox("Select pixel maps directory",          false);
    		gd.addStringField("Results directory",                 this.resultsDirectory, 40);
    		gd.addCheckbox("Select results directory",             false);
    		gd.addStringField("Source files prefix",               this.sourcePrefix, 40);
    		gd.addStringField("Source files suffix",               this.sourceSuffix, 40);
    		gd.addStringField("Plane maps prefix",                 this.planeMapPrefix, 40);
    		gd.addStringField("Plane maps suffix",                 this.planeMapSuffix, 40);



    		WindowTools.addScrollBars(gd);
    		gd.showDialog();
    		if (gd.wasCanceled()) return false;
    		this.saveSettings=              gd.getNextBoolean();
    		this.sourceDirectory=           gd.getNextString(); if (gd.getNextBoolean()) selectSourceDirectory(false, false);
    		this.equirectangularDirectory=  gd.getNextString(); if (gd.getNextBoolean()) selectMapsDirectory(false, false);
    		this.resultsDirectory=          gd.getNextString(); if (gd.getNextBoolean()) selectResultsDirectory(false, true);
    		this.sourcePrefix=              gd.getNextString();
    		this.sourceSuffix=              gd.getNextString();
    		this.planeMapPrefix=  gd.getNextString();
    		this.planeMapSuffix=  gd.getNextString();

    		return true;
    	}


// TODO: extract timestamnp from JP4 or, at least combine movie timestamp+frame into a single filename string
    	public String [] getSourcePaths(){
    		String [] empty={};
    		return (this.sourcePaths!=null)?this.sourcePaths:empty;
    	}
    	public String selectSourceDirectory(boolean smart, boolean newAllowed) { // normally newAllowed=false
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save
    				"Source (\"de-warped\") image directory", // title
    				"Select source directory", // button
    				null, // filter
    				this.sourceDirectory); // this.sourceDirectory);
    		if (dir!=null) this.sourceDirectory=dir;
    		return dir;
    	}
    	public String selectMapsDirectory(boolean smart, boolean newAllowed) {
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save
    				"Maps directory", // title
    				"Select maps directory", // button
    				null, // filter
    				this.equirectangularDirectory);
    		if (dir!=null) this.equirectangularDirectory=dir;
    		return dir;
    	}

    	public String selectResultsDirectory(boolean smart, boolean newAllowed) {
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save
    				"Results directory", // title
    				"Select results directory", // button
    				null, // filter
    				this.resultsDirectory); //this.sourceDirectory);
    		if (dir!=null) this.resultsDirectory=dir;
    		return dir;
    	}
    	public boolean selectSourceFiles(boolean allFiles) {
    		return selectSourceFiles(allFiles, 1); // debug level 1 - modify here
    	}
    	public boolean selectSourceFiles(boolean allFiles, int debugLevel) {
    		String [] defaultPaths=this.sourcePaths;
    		if ((defaultPaths==null) || (defaultPaths.length==0)){
    			defaultPaths = new String[1];
    			if ((this.sourceDirectory==null) || (this.sourceDirectory.length()==0)){
    				defaultPaths[0]="";
    			} else {
    				defaultPaths[0]=this.sourceDirectory+Prefs.getFileSeparator();
    			}
    		}
    		String [] extensions={this.sourceSuffix};
			MultipleExtensionsFileFilter sourceFilter =
				new MultipleExtensionsFileFilter(this.sourcePrefix,extensions,"Source files");
			String [] sourceFiles=null;
    		if (allFiles){
				File dir= new File (this.sourceDirectory);
				if (debugLevel>1) System.out.println("selectSourceFiles, dir="+this.sourceDirectory);
				if (!dir.exists()) {
					String error="Source directory "+this.sourceDirectory+" does not exist.";
					if (debugLevel>1) System.out.println("selectSourceFiles() ERROR:"+error);
					if (debugLevel>1) IJ.showMessage("No files selected");
					return false;
				}
				File [] fileList=dir.listFiles(sourceFilter);
				if (debugLevel>1) System.out.println("Source directory "+this.sourceDirectory+" has "+fileList.length+" files.");
				sourceFiles = new String[fileList.length];
				for (int i=0;i<sourceFiles.length;i++) sourceFiles[i]=fileList[i].getPath();
    		} else {
    			sourceFiles=CalibrationFileManagement.selectFiles(false,
    					"Select Source files, saved as TIFF",
    					"Select",
    					sourceFilter,
    					defaultPaths); // String [] defaultPaths); //this.sourceDirectory // null
    	       	if ((sourceFiles==null) || (sourceFiles.length==0)) {
					if (debugLevel>1) System.out.println("selectSourceFiles() ERROR: No files selected");
					if (debugLevel>1) IJ.showMessage("No files selected");
            		return false;
            	}
    		}
    		this.sourcePaths=sourceFiles;
    		if ((this.sourcePaths!=null) && (this.sourcePaths.length>0)){
    			this.sourceDirectory=this.sourcePaths[0].substring(0, this.sourcePaths[0].lastIndexOf(Prefs.getFileSeparator()));
    			//sourceNames
    		}
			return true;
    	}
    	public int getChannelFromSourceTiff(String path){ return getChannelFromTiff(path, this.sourceSuffix);	}
    	public int getChannelFromTiff(String path, String suffix){
    		int indexSuffix=path.length()-suffix.length();
    		int indexLastDash=indexSuffix-1; // in jp4 it will be underscore, not dash? Or should we change that?
    		while ((indexLastDash>0) &&
    				(indexLastDash>(indexSuffix-3)) &&
    				(path.charAt(indexLastDash)!='_') &&
    				(path.charAt(indexLastDash)!='-')) indexLastDash--;
    		return Integer.parseInt(path.substring(indexLastDash+1,indexSuffix));

    	}


	}
	public class DisparityCorrelationParameters{
    	public double [][] disparityScales=       null; // for each channel - a pair of {scaleX, scaleY} or null if undefined (interSensor has the same)
    	public boolean autocorrelation=           false;
    	public int corrFFTSize=                 128;
    	public double corrXC=                  1280.0;
    	public double corrYC=                   960.0;
    	public double corrShift=                 30.0; // was int corrShift->disparity
    	public int    firstImage=                 0;
    	public int    secondImage=                1;

    	public double corrPhaseFraction=          0.5;
    	public double corrCbWeight=               0.5;
    	public double corrCrWeight=               0.5;
    	public double correlationHighPassSigma=   1.5;
    	public double correlationLowPassSigma=    0.6;
    	public int    normalizeCorrelation=         1; // 0 - none, 1 - "old", 2 - new
    	public double noiseNormalizationSignaY=   3.0;
    	public double noiseNormalizationSignaCbCr=5.0;
    	public double contrastThreshold=          1.5;
    	public boolean useBinaryAlpha=            true;

    	public int    tileOverlapFraction=        4; // calculate FHT tiles with corrFFTSize/tileOverlapFraction overlap
    	public double disparityRange=           100.0;  // maximal allowed disparity in pixels
    	public int    disparitySteps=           100;    // calculate Z with disparitySteps [disparitySteps] over disparityRange
    	public boolean useFullWindow=             true;
    	public Rectangle resultWindow= new Rectangle(0,0,2592,1936); // only process for the result pixels in the specified rectangle
    	public int    channelMask=                1; // 1 - process Y, +2 process Cb, +4 - process Cr
    	public boolean useFileData=             false; // use multi-slice TIFF file (1 layer per image) instead of channels (alfa will still be from channels)
    	public String  dataPathName=              "";  // external data image path
    	public String  corrPathName=              "";  // external data image path

    	public int    interpolationSize=          8;     // use interpolationSizexinterpolationSize square to up-sample center 2x2 pixel area
    	public int    interpolationUpSample=     16;    // subdivide pixels near local maximums on correlation vs. disparity
    	public double[]  subpixAMax=null;
    	public double subpixRMax=                 0.05;
    	public int    subpixNMax=                 4;
    	public boolean showCorrelationToImage=  true;
    	public boolean saveCorrelationToImage=  true;
    	public DisparityCorrelationParameters(){
    		this.subpixAMax=new double[4];
    		this.subpixAMax[0]=0.0;
    		this.subpixAMax[1]=0.0;
    		this.subpixAMax[2]=0.0;
    		this.subpixAMax[3]=0.0;
    	}
    	public double bgFraction=              0.5; // try the background plane if it is not less than this fraction of maximal correlation plane
    	public double blurVarianceSigma=       2.0;
    	public int    refineTilePeriod=        4;
    	public double refinePhaseCoeff=        0.5;
    	public double refineHighPassSigma=     1.5;
    	public double refineLowPassSigma=      0.6;
    	public double refineCorrMaxDistance=   1.5;
    	public double refineCorrThreshold=     1.0; //??
    	public int    refineSubPixel=          8;
 // during setup
    	public Rectangle  zMapWOI=             null;
    	public int    zMapMaxNumber=           3;
    	public double zMapMinFirst  =          0.015;
    	public double zMapMinAbsolute  =       0.002;
    	public double zMapMinRelative  =       0.15;
    	public double zMapMergeMax  =          1.5;
    	public int    zMapOverlap  =           2;
    	public double zMapFatZero=             0.05; // actually - for combining to center
// during re-run/filter
    	public double zMapMinForeground =      0.001;
    	public int    zMapVarMask=             7; // +1 - Y, +2 - Cb, +4 - Cr, +8 - Ext
    	public double [] zMapVarThresholds=    {2.5,2.5,2.5,2.5};
    	public double [] zMapVarWeights=       {1.0,0.5,0.5,0.5};//
    	public int    auxVarMode=              0; // 0 multiply all, add pairs
    	public int    zMapCorrMask=            7; // +1 - Y, +2 - Cb, +4 - Cr, +8 - Ext
    	public double [] zMapCorrThresholds=   {0.0,0.0,0.0,0.0}; // 0 - positive correlation, >+1.0 - strong, negative - loosen
    	public double [] zMapCorrThresholdsRel={0.5,0.5,0.5,0.5}; // cross correlation to sqrt(auto1*auto2) ratio
    	public double [] zMapCorrWeights=      {1.0,0.5,0.5,0.5};//
    	public double fillFgGapMin=            0.001;
    	public double fillFgGapDiff=           1.5;
    	public int    fillFgGapNeib=           3;

    	public double    filter2DisparityMax =100.0;
    	public double    filter2disparityMin = 50.0;
    	public boolean   filter2UpdateMax=     true;
    	public double    filter2MinAbsolute=    0.01;
    	public double    filter2MinRelative=    0.1;
    	public boolean   filter2ByForeground= true; // apply known certain masks
    	public double    filter2ByForegroundMargin=2.0;
    	public boolean   filter2ByDisabled=   true;
    	public double    filter2DisparityTolearnce=1.0;
    	public double    filter2MaskBlurSigma=   2.0;
    	public double    filter2corrHighPassSigma=2.0; // subtract blurred version to minimize correlation caused by masks

    	//photometric
    	public double    photometricIgnoreFraction=        0.001;
    	public int       photometricSubdivAverage=       256;
    	public int       photometricSubdivHalfDifference=128;
    	public double    photometricSmoothVarianceSigma=  10.0;
    	public double    photometricScaleVariance=         3.0;

    	public double    matchStatVarianceBlurScale=1.0;
    	public double    matchStatKLocal=          0.5; // 0 use global variance, 1.0 - local
    	public int       matchStatMode=              0; // 0 - multiply likelyhood for all pairs, 1 - add, 2 - min, 3 - max


		public int debugRow=                   -10;   // "Debug row"
		public int debugColumn=                -10;   // "Debug column"

 //   	public boolean enableNegativeDisparity=   false;
    	public void setProperties(String prefix,Properties properties){
  			properties.setProperty(prefix+"autocorrelation",this.autocorrelation+"");
  			properties.setProperty(prefix+"corrFFTSize",this.corrFFTSize+"");
  			properties.setProperty(prefix+"corrXC",this.corrXC+"");
  			properties.setProperty(prefix+"corrYC",this.corrYC+"");
  			properties.setProperty(prefix+"corrShift",this.corrShift+"");
  			properties.setProperty(prefix+"firstImage",this.firstImage+"");
  			properties.setProperty(prefix+"secondImage",this.secondImage+"");
  			properties.setProperty(prefix+"corrPhaseFraction",this.corrPhaseFraction+"");
  			properties.setProperty(prefix+"corrCbWeight",this.corrCbWeight+"");
  			properties.setProperty(prefix+"corrCrWeight",this.corrCrWeight+"");
  			properties.setProperty(prefix+"correlationHighPassSigma",this.correlationHighPassSigma+"");
  			properties.setProperty(prefix+"correlationLowPassSigma",this.correlationLowPassSigma+"");
  			properties.setProperty(prefix+"normalizeCorrelation",this.normalizeCorrelation+"");
  			properties.setProperty(prefix+"noiseNormalizationSignaY",this.noiseNormalizationSignaY+"");
  			properties.setProperty(prefix+"noiseNormalizationSignaCbCr",this.noiseNormalizationSignaCbCr+"");
  			properties.setProperty(prefix+"contrastThreshold",this.contrastThreshold+"");
  			properties.setProperty(prefix+"useBinaryAlpha",this.useBinaryAlpha+"");
  			if (disparityScales!=null) {
  	  			properties.setProperty(prefix+"disparityScales",this.disparityScales.length+"");
  				for (int i=0;i<this.disparityScales.length;i++) if (this.disparityScales[i]!=null){
  		  			properties.setProperty(prefix+"disparityScales_"+i+"_X",this.disparityScales[i][0]+"");
  		  			properties.setProperty(prefix+"disparityScales_"+i+"_Y",this.disparityScales[i][1]+"");
  				}
  			}
  			properties.setProperty(prefix+"tileOverlapFraction",this.tileOverlapFraction+"");
  			properties.setProperty(prefix+"disparityRange",this.disparityRange+"");
  			properties.setProperty(prefix+"disparitySteps",this.disparitySteps+"");
  			properties.setProperty(prefix+"useFullWindow",this.useFullWindow+"");
  			properties.setProperty(prefix+"resultWindow_x",this.resultWindow.x+"");
  			properties.setProperty(prefix+"resultWindow_y",this.resultWindow.y+"");
  			properties.setProperty(prefix+"resultWindow_width",this.resultWindow.width+"");
  			properties.setProperty(prefix+"resultWindow_height",this.resultWindow.height+"");
  			properties.setProperty(prefix+"channelMask",this.channelMask+"");
  			properties.setProperty(prefix+"useFileData",this.useFileData+"");
  			properties.setProperty(prefix+"dataPathName",this.dataPathName);
  			properties.setProperty(prefix+"corrPathName",this.corrPathName);
  			properties.setProperty(prefix+"interpolationSize",this.interpolationSize+"");
  			properties.setProperty(prefix+"interpolationUpSample",this.interpolationUpSample+"");
  			properties.setProperty(prefix+"subpixAMax_Y",this.subpixAMax[0]+"");
  			properties.setProperty(prefix+"subpixAMax_Cb",this.subpixAMax[1]+"");
  			properties.setProperty(prefix+"subpixAMax_Cr",this.subpixAMax[2]+"");
  			properties.setProperty(prefix+"subpixAMax_Ext",this.subpixAMax[3]+"");
  			properties.setProperty(prefix+"subpixRMax",this.subpixRMax+"");
  			properties.setProperty(prefix+"subpixNMax",this.subpixNMax+"");
  			properties.setProperty(prefix+"showCorrelationToImage",this.showCorrelationToImage+"");
  			properties.setProperty(prefix+"saveCorrelationToImage",this.saveCorrelationToImage+"");

  			properties.setProperty(prefix+"bgFraction",this.bgFraction+"");
  			properties.setProperty(prefix+"blurVarianceSigma",this.blurVarianceSigma+"");
  			properties.setProperty(prefix+"refineTilePeriod",this.refineTilePeriod+"");
  			properties.setProperty(prefix+"refinePhaseCoeff",this.refinePhaseCoeff+"");
  			properties.setProperty(prefix+"refineHighPassSigma",this.refineHighPassSigma+"");
  			properties.setProperty(prefix+"refineLowPassSigma",this.refineLowPassSigma+"");
  			properties.setProperty(prefix+"refineCorrMaxDistance",this.refineCorrMaxDistance+"");
  			properties.setProperty(prefix+"refineCorrThreshold",this.refineCorrThreshold+"");
  			properties.setProperty(prefix+"refineSubPixel",this.refineSubPixel+"");

  			if (this.zMapWOI!=null) {
  				properties.setProperty(prefix+"zMapWOI_x",this.zMapWOI.x+"");
  				properties.setProperty(prefix+"zMapWOI_y",this.zMapWOI.y+"");
  				properties.setProperty(prefix+"zMapWOI_width",this.zMapWOI.width+"");
  				properties.setProperty(prefix+"zMapWOI_height",this.zMapWOI.height+"");
  			}
  			properties.setProperty(prefix+"zMapMaxNumber",this.zMapMaxNumber+"");
  			properties.setProperty(prefix+"zMapMinFirst",this.zMapMinFirst+"");
  			properties.setProperty(prefix+"zMapMinAbsolute",this.zMapMinAbsolute+"");
  			properties.setProperty(prefix+"zMapMinRelative",this.zMapMinRelative+"");
  			properties.setProperty(prefix+"zMapMergeMax",this.zMapMergeMax+"");
  			properties.setProperty(prefix+"zMapOverlap",this.zMapOverlap+"");
  			properties.setProperty(prefix+"zMapFatZero",this.zMapFatZero+"");

  			properties.setProperty(prefix+"zMapMinForeground",this.zMapMinForeground+"");
  			properties.setProperty(prefix+"zMapVarMask",this.zMapVarMask+"");
  			properties.setProperty(prefix+"zMapVarThresholds_length",this.zMapVarThresholds.length+"");
  			for (int i=0;i<zMapVarThresholds.length;i++){
  				properties.setProperty(prefix+"zMapVarThresholds_"+i,this.zMapVarThresholds[i]+"");
  			}
  			properties.setProperty(prefix+"auxVarMode",this.auxVarMode+"");

  			properties.setProperty(prefix+"zMapVarWeights_length",this.zMapVarWeights.length+"");
  			for (int i=0;i<zMapVarWeights.length;i++){
  				properties.setProperty(prefix+"zMapVarWeights_"+i,this.zMapVarWeights[i]+"");
  			}
  			properties.setProperty(prefix+"zMapCorrMask",this.zMapCorrMask+"");
  			properties.setProperty(prefix+"zMapCorrThresholds_length",this.zMapCorrThresholds.length+"");
  			for (int i=0;i<zMapCorrThresholds.length;i++){
  				properties.setProperty(prefix+"zMapCorrThresholds_"+i,this.zMapCorrThresholds[i]+"");
  			}
  			properties.setProperty(prefix+"zMapCorrThresholdsRel_length",this.zMapCorrThresholdsRel.length+"");
  			for (int i=0;i<zMapCorrThresholdsRel.length;i++){
  				properties.setProperty(prefix+"zMapCorrThresholdsRel_"+i,this.zMapCorrThresholdsRel[i]+"");
  			}
  			properties.setProperty(prefix+"zMapCorrWeights_length",this.zMapCorrWeights.length+"");
  			for (int i=0;i<zMapCorrWeights.length;i++){
  				properties.setProperty(prefix+"zMapCorrWeights_"+i,this.zMapCorrWeights[i]+"");
  			}
  			properties.setProperty(prefix+"fillFgGapMin",this.fillFgGapMin+"");
  			properties.setProperty(prefix+"fillFgGapDiff",this.fillFgGapDiff+"");
  			properties.setProperty(prefix+"fillFgGapNeib",this.fillFgGapNeib+"");

  			properties.setProperty(prefix+"filter2DisparityMax",this.filter2DisparityMax+"");
  			properties.setProperty(prefix+"filter2disparityMin",this.filter2disparityMin+"");
  			properties.setProperty(prefix+"filter2UpdateMax",this.filter2UpdateMax+"");
  			properties.setProperty(prefix+"filter2MinAbsolute",this.filter2MinAbsolute+"");
  			properties.setProperty(prefix+"filter2MinRelative",this.filter2MinRelative+"");
  			properties.setProperty(prefix+"filter2ByForeground",this.filter2ByForeground+"");
  			properties.setProperty(prefix+"filter2ByForegroundMargin",this.filter2ByForegroundMargin+"");
  			properties.setProperty(prefix+"filter2ByDisabled",this.filter2ByDisabled+"");
  			properties.setProperty(prefix+"filter2DisparityTolearnce",this.filter2DisparityTolearnce+"");
  			properties.setProperty(prefix+"filter2MaskBlurSigma",this.filter2MaskBlurSigma+"");
  			properties.setProperty(prefix+"filter2corrHighPassSigma",this.filter2corrHighPassSigma+"");

  			properties.setProperty(prefix+"photometricIgnoreFraction",this.photometricIgnoreFraction+"");
  			properties.setProperty(prefix+"photometricSubdivAverage",this.photometricSubdivAverage+"");
  			properties.setProperty(prefix+"photometricSubdivHalfDifference",this.photometricSubdivHalfDifference+"");
  			properties.setProperty(prefix+"photometricSmoothVarianceSigma",this.photometricSmoothVarianceSigma+"");
  			properties.setProperty(prefix+"photometricScaleVariance",this.photometricScaleVariance+"");

  			properties.setProperty(prefix+"matchStatVarianceBlurScale",this.matchStatVarianceBlurScale+"");
  			properties.setProperty(prefix+"matchStatKLocal",this.matchStatKLocal+"");
  			properties.setProperty(prefix+"matchStatMode",this.matchStatMode+"");

    		properties.setProperty(prefix+"debugRow",this.debugRow+"");
    		properties.setProperty(prefix+"debugColumn",this.debugColumn+"");

    	}
    	public void getProperties(String prefix,Properties properties){
    		if (this.resultWindow==null) this.resultWindow=new Rectangle();

    		if (properties.getProperty(prefix+"autocorrelation")!=null)             this.autocorrelation=Boolean.parseBoolean(properties.getProperty(prefix+"autocorrelation"));
  		    if (properties.getProperty(prefix+"corrFFTSize")!=null)                 this.corrFFTSize=Integer.parseInt(properties.getProperty(prefix+"corrFFTSize"));
  		    if (properties.getProperty(prefix+"corrXC")!=null)                      this.corrXC=Double.parseDouble(properties.getProperty(prefix+"corrXC"));
  		    if (properties.getProperty(prefix+"corrYC")!=null)                      this.corrYC=Double.parseDouble(properties.getProperty(prefix+"corrYC"));
  		    if (properties.getProperty(prefix+"corrShift")!=null)                   this.corrShift=Double.parseDouble(properties.getProperty(prefix+"corrShift"));
  		    if (properties.getProperty(prefix+"firstImage")!=null)                  this.firstImage=Integer.parseInt(properties.getProperty(prefix+"firstImage"));
  		    if (properties.getProperty(prefix+"secondImage")!=null)                 this.secondImage=Integer.parseInt(properties.getProperty(prefix+"secondImage"));
  		    if (properties.getProperty(prefix+"corrPhaseFraction")!=null)           this.corrPhaseFraction=Double.parseDouble(properties.getProperty(prefix+"corrPhaseFraction"));
  		    if (properties.getProperty(prefix+"corrCbWeight")!=null)                this.corrCbWeight=Double.parseDouble(properties.getProperty(prefix+"corrCbWeight"));
  		    if (properties.getProperty(prefix+"corrCrWeight")!=null)                this.corrCrWeight=Double.parseDouble(properties.getProperty(prefix+"corrCrWeight"));
  		    if (properties.getProperty(prefix+"correlationHighPassSigma")!=null)    this.correlationHighPassSigma=Double.parseDouble(properties.getProperty(prefix+"correlationHighPassSigma"));
  		    if (properties.getProperty(prefix+"correlationLowPassSigma")!=null)     this.correlationLowPassSigma=Double.parseDouble(properties.getProperty(prefix+"correlationLowPassSigma"));
    		if (properties.getProperty(prefix+"normalizeCorrelation")!=null)        this.normalizeCorrelation=Integer.parseInt(properties.getProperty(prefix+"normalizeCorrelation"));
  		    if (properties.getProperty(prefix+"noiseNormalizationSignaY")!=null)    this.noiseNormalizationSignaY=Double.parseDouble(properties.getProperty(prefix+"noiseNormalizationSignaY"));
  		    if (properties.getProperty(prefix+"noiseNormalizationSignaCbCr")!=null) this.noiseNormalizationSignaCbCr=Double.parseDouble(properties.getProperty(prefix+"noiseNormalizationSignaCbCr"));
  		    if (properties.getProperty(prefix+"contrastThreshold")!=null)           this.contrastThreshold=Double.parseDouble(properties.getProperty(prefix+"contrastThreshold"));
    		if (properties.getProperty(prefix+"useBinaryAlpha")!=null)              this.useBinaryAlpha=Boolean.parseBoolean(properties.getProperty(prefix+"useBinaryAlpha"));
    		if (properties.getProperty(prefix+"disparityScales")!=null){
    			this.disparityScales=new double [Integer.parseInt(properties.getProperty(prefix+"disparityScales"))][];
    			for (int i=0;i<this.disparityScales.length;i++) {
    				if ((properties.getProperty(prefix+"disparityScales_"+i+"_X")!=null) && (properties.getProperty(prefix+"disparityScales_"+i+"_Y")!=null)){
    					this.disparityScales[i]=new double[2];
    					this.disparityScales[i][0]=Double.parseDouble(properties.getProperty(prefix+"disparityScales_"+i+"_X"));
    					this.disparityScales[i][1]=Double.parseDouble(properties.getProperty(prefix+"disparityScales_"+i+"_Y"));
    				} else this.disparityScales[i]=null;
    			}
    		}
  		    if (properties.getProperty(prefix+"tileOverlapFraction")!=null)  this.tileOverlapFraction=Integer.parseInt(properties.getProperty(prefix+"tileOverlapFraction"));
  		    if (properties.getProperty(prefix+"disparityRange")!=null)       this.disparityRange=Double.parseDouble(properties.getProperty(prefix+"disparityRange"));
  		    if (properties.getProperty(prefix+"disparitySteps")!=null)       this.disparitySteps=Integer.parseInt(properties.getProperty(prefix+"disparitySteps"));
    		if (properties.getProperty(prefix+"useFullWindow")!=null)        this.useFullWindow=Boolean.parseBoolean(properties.getProperty(prefix+"useFullWindow"));
  		    if (properties.getProperty(prefix+"resultWindow_x")!=null)       this.resultWindow.x=Integer.parseInt(properties.getProperty(prefix+"resultWindow_x"));
  		    if (properties.getProperty(prefix+"resultWindow_y")!=null)       this.resultWindow.y=Integer.parseInt(properties.getProperty(prefix+"resultWindow_y"));
  		    if (properties.getProperty(prefix+"resultWindow_width")!=null)   this.resultWindow.width= Integer.parseInt(properties.getProperty(prefix+"resultWindow_width"));
  		    if (properties.getProperty(prefix+"resultWindow_height")!=null)  this.resultWindow.height=Integer.parseInt(properties.getProperty(prefix+"resultWindow_height"));
		    if (properties.getProperty(prefix+"channelMask")!=null) this.channelMask=Integer.parseInt(properties.getProperty(prefix+"channelMask"));
    		if (properties.getProperty(prefix+"useFileData")!=null)          this.useFileData=Boolean.parseBoolean(properties.getProperty(prefix+"useFileData"));
    		if (properties.getProperty(prefix+"dataPathName")!=null)         this.dataPathName=(properties.getProperty(prefix+"dataPathName"));
    		if (properties.getProperty(prefix+"corrPathName")!=null)         this.corrPathName=(properties.getProperty(prefix+"corrPathName"));
		    if (properties.getProperty(prefix+"interpolationSize")!=null)    this.interpolationSize=Integer.parseInt(properties.getProperty(prefix+"interpolationSize"));
		    if (properties.getProperty(prefix+"interpolationUpSample")!=null)this.interpolationUpSample=Integer.parseInt(properties.getProperty(prefix+"interpolationUpSample"));
  		    if (properties.getProperty(prefix+"subpixAMax_Y")!=null)         this.subpixAMax[0]=Double.parseDouble(properties.getProperty(prefix+"subpixAMax_Y"));
  		    if (properties.getProperty(prefix+"subpixAMax_Cb")!=null)        this.subpixAMax[1]=Double.parseDouble(properties.getProperty(prefix+"subpixAMax_Cb"));
  		    if (properties.getProperty(prefix+"subpixAMax_Cr")!=null)        this.subpixAMax[2]=Double.parseDouble(properties.getProperty(prefix+"subpixAMax_Cr"));
  		    if (properties.getProperty(prefix+"subpixAMax_Ext")!=null)       this.subpixAMax[3]=Double.parseDouble(properties.getProperty(prefix+"subpixAMax_Ext"));
  		    if (properties.getProperty(prefix+"subpixRMax")!=null)           this.subpixRMax=Double.parseDouble(properties.getProperty(prefix+"subpixRMax"));
		    if (properties.getProperty(prefix+"subpixNMax")!=null)           this.subpixNMax=Integer.parseInt(properties.getProperty(prefix+"subpixNMax"));
    		if (properties.getProperty(prefix+"showCorrelationToImage")!=null)this.showCorrelationToImage=Boolean.parseBoolean(properties.getProperty(prefix+"showCorrelationToImage"));
    		if (properties.getProperty(prefix+"saveCorrelationToImage")!=null)this.saveCorrelationToImage=Boolean.parseBoolean(properties.getProperty(prefix+"saveCorrelationToImage"));

    		if (properties.getProperty(prefix+"bgFraction")!=null)           this.bgFraction=Double.parseDouble(properties.getProperty(prefix+"bgFraction"));
    		if (properties.getProperty(prefix+"blurVarianceSigma")!=null)    this.blurVarianceSigma=Double.parseDouble(properties.getProperty(prefix+"blurVarianceSigma"));
		    if (properties.getProperty(prefix+"refineTilePeriod")!=null)     this.refineTilePeriod=Integer.parseInt(properties.getProperty(prefix+"refineTilePeriod"));
  		    if (properties.getProperty(prefix+"refinePhaseCoeff")!=null)     this.refinePhaseCoeff=Double.parseDouble(properties.getProperty(prefix+"refinePhaseCoeff"));
  		    if (properties.getProperty(prefix+"refineHighPassSigma")!=null)  this.refineHighPassSigma=Double.parseDouble(properties.getProperty(prefix+"refineHighPassSigma"));
  		    if (properties.getProperty(prefix+"refineLowPassSigma")!=null)   this.refineLowPassSigma=Double.parseDouble(properties.getProperty(prefix+"refineLowPassSigma"));
  		    if (properties.getProperty(prefix+"refineCorrMaxDistance")!=null)this.refineCorrMaxDistance=Double.parseDouble(properties.getProperty(prefix+"refineCorrMaxDistance"));
  		    if (properties.getProperty(prefix+"refineCorrThreshold")!=null)  this.refineCorrThreshold=Double.parseDouble(properties.getProperty(prefix+"refineCorrThreshold"));
		    if (properties.getProperty(prefix+"refineSubPixel")!=null)       this.refineSubPixel=Integer.parseInt(properties.getProperty(prefix+"refineSubPixel"));
		    if (properties.getProperty(prefix+"zMapWOI_x")!=null) { // should be defined all or none
		    	this.zMapWOI=new Rectangle();
		    	this.zMapWOI.x=Integer.parseInt(properties.getProperty(prefix+"zMapWOI_x"));
		    	this.zMapWOI.y=Integer.parseInt(properties.getProperty(prefix+"zMapWOI_y"));
		    	this.zMapWOI.width= Integer.parseInt(properties.getProperty(prefix+"zMapWOI_width"));
		    	this.zMapWOI.height=Integer.parseInt(properties.getProperty(prefix+"zMapWOI_height"));
		    }
		    if (properties.getProperty(prefix+"zMapMaxNumber")!=null) this.zMapMaxNumber=Integer.parseInt(properties.getProperty(prefix+"zMapMaxNumber"));
  		    if (properties.getProperty(prefix+"zMapMinFirst")!=null) this.zMapMinFirst=Double.parseDouble(properties.getProperty(prefix+"zMapMinFirst"));
  		    if (properties.getProperty(prefix+"zMapMinAbsolute")!=null) this.zMapMinAbsolute=Double.parseDouble(properties.getProperty(prefix+"zMapMinAbsolute"));
  		    if (properties.getProperty(prefix+"zMapMinRelative")!=null) this.zMapMinRelative=Double.parseDouble(properties.getProperty(prefix+"zMapMinRelative"));
  		    if (properties.getProperty(prefix+"zMapMergeMax")!=null) this.zMapMergeMax=Double.parseDouble(properties.getProperty(prefix+"zMapMergeMax"));
		    if (properties.getProperty(prefix+"zMapOverlap")!=null) this.zMapOverlap=Integer.parseInt(properties.getProperty(prefix+"zMapOverlap"));
  		    if (properties.getProperty(prefix+"zMapFatZero")!=null) this.zMapFatZero=Double.parseDouble(properties.getProperty(prefix+"zMapFatZero"));


  		    if (properties.getProperty(prefix+"zMapMinForeground")!=null) this.zMapMinForeground=Double.parseDouble(properties.getProperty(prefix+"zMapMinForeground"));
		    if (properties.getProperty(prefix+"zMapVarMask")!=null) this.zMapVarMask=Integer.parseInt(properties.getProperty(prefix+"zMapVarMask"));
		    if (properties.getProperty(prefix+"zMapVarThresholds_length")!=null){
		    	this.zMapVarThresholds= new double[Integer.parseInt(properties.getProperty(prefix+"zMapVarThresholds_length"))];
	  			for (int i=0;i<zMapVarThresholds.length;i++){
	  	  		    if (properties.getProperty(prefix+"zMapVarThresholds_"+i)!=null) this.zMapVarThresholds[i]=Double.parseDouble(properties.getProperty(prefix+"zMapVarThresholds_"+i));
	  			}
		    }
		    if (properties.getProperty(prefix+"auxVarMode")!=null) this.auxVarMode=Integer.parseInt(properties.getProperty(prefix+"auxVarMode"));

		    if (properties.getProperty(prefix+"zMapVarWeights_length")!=null){
		    	this.zMapVarWeights= new double[Integer.parseInt(properties.getProperty(prefix+"zMapVarWeights_length"))];
	  			for (int i=0;i<zMapVarWeights.length;i++){
	  	  		    if (properties.getProperty(prefix+"zMapVarWeights_"+i)!=null) this.zMapVarWeights[i]=Double.parseDouble(properties.getProperty(prefix+"zMapVarWeights_"+i));
	  			}
		    }

		    if (properties.getProperty(prefix+"zMapCorrMask")!=null) this.zMapCorrMask=Integer.parseInt(properties.getProperty(prefix+"zMapCorrMask"));
		    if (properties.getProperty(prefix+"zMapCorrThresholds_length")!=null){
		    	this.zMapCorrThresholds= new double[Integer.parseInt(properties.getProperty(prefix+"zMapCorrThresholds_length"))];
	  			for (int i=0;i<zMapCorrThresholds.length;i++){
	  	  		    if (properties.getProperty(prefix+"zMapCorrThresholds_"+i)!=null) this.zMapCorrThresholds[i]=Double.parseDouble(properties.getProperty(prefix+"zMapCorrThresholds_"+i));
	  			}
		    }
		    if (properties.getProperty(prefix+"zMapCorrThresholdsRel_length")!=null){
		    	this.zMapCorrThresholdsRel= new double[Integer.parseInt(properties.getProperty(prefix+"zMapCorrThresholdsRel_length"))];
	  			for (int i=0;i<zMapCorrThresholdsRel.length;i++){
	  	  		    if (properties.getProperty(prefix+"zMapCorrThresholdsRel_"+i)!=null) this.zMapCorrThresholdsRel[i]=Double.parseDouble(properties.getProperty(prefix+"zMapCorrThresholdsRel_"+i));
	  			}
		    }
		    if (properties.getProperty(prefix+"zMapCorrWeights_length")!=null){
		    	this.zMapCorrWeights= new double[Integer.parseInt(properties.getProperty(prefix+"zMapCorrWeights_length"))];
	  			for (int i=0;i<zMapCorrWeights.length;i++){
	  	  		    if (properties.getProperty(prefix+"zMapCorrWeights_"+i)!=null) this.zMapCorrWeights[i]=Double.parseDouble(properties.getProperty(prefix+"zMapCorrWeights_"+i));
	  			}
		    }
  		    if (properties.getProperty(prefix+"fillFgGapMin")!=null)  this.fillFgGapMin=Double.parseDouble(properties.getProperty(prefix+"fillFgGapMin"));
  		    if (properties.getProperty(prefix+"fillFgGapDiff")!=null) this.fillFgGapDiff=Double.parseDouble(properties.getProperty(prefix+"fillFgGapDiff"));
		    if (properties.getProperty(prefix+"fillFgGapNeib")!=null) this.fillFgGapNeib=Integer.parseInt(properties.getProperty(prefix+"fillFgGapNeib"));

		    if (properties.getProperty(prefix+"filter2DisparityMax")!=null) this.filter2DisparityMax=Double.parseDouble(properties.getProperty(prefix+"filter2DisparityMax"));
  		    if (properties.getProperty(prefix+"filter2disparityMin")!=null) this.filter2disparityMin=Double.parseDouble(properties.getProperty(prefix+"filter2disparityMin"));
    		if (properties.getProperty(prefix+"filter2UpdateMax")!=null)this.filter2UpdateMax=Boolean.parseBoolean(properties.getProperty(prefix+"filter2UpdateMax"));
  		    if (properties.getProperty(prefix+"filter2MinAbsolute")!=null) this.filter2MinAbsolute=Double.parseDouble(properties.getProperty(prefix+"filter2MinAbsolute"));
  		    if (properties.getProperty(prefix+"filter2MinRelative")!=null) this.filter2MinRelative=Double.parseDouble(properties.getProperty(prefix+"filter2MinRelative"));
    		if (properties.getProperty(prefix+"filter2ByForeground")!=null)this.filter2ByForeground=Boolean.parseBoolean(properties.getProperty(prefix+"filter2ByForeground"));
  		    if (properties.getProperty(prefix+"filter2ByForegroundMargin")!=null) this.filter2ByForegroundMargin=Double.parseDouble(properties.getProperty(prefix+"filter2ByForegroundMargin"));
    		if (properties.getProperty(prefix+"filter2ByDisabled")!=null)this.filter2ByDisabled=Boolean.parseBoolean(properties.getProperty(prefix+"filter2ByDisabled"));
  		    if (properties.getProperty(prefix+"filter2DisparityTolearnce")!=null) this.filter2DisparityTolearnce=Double.parseDouble(properties.getProperty(prefix+"filter2DisparityTolearnce"));
  		    if (properties.getProperty(prefix+"filter2MaskBlurSigma")!=null) this.filter2MaskBlurSigma=Double.parseDouble(properties.getProperty(prefix+"filter2MaskBlurSigma"));
  		    if (properties.getProperty(prefix+"filter2corrHighPassSigma")!=null) this.filter2corrHighPassSigma=Double.parseDouble(properties.getProperty(prefix+"filter2corrHighPassSigma"));

  		    if (properties.getProperty(prefix+"photometricIgnoreFraction")!=null) this.photometricIgnoreFraction=Double.parseDouble(properties.getProperty(prefix+"photometricIgnoreFraction"));
		    if (properties.getProperty(prefix+"photometricSubdivAverage")!=null) this.photometricSubdivAverage=Integer.parseInt(properties.getProperty(prefix+"photometricSubdivAverage"));
		    if (properties.getProperty(prefix+"photometricSubdivHalfDifference")!=null) this.photometricSubdivHalfDifference=Integer.parseInt(properties.getProperty(prefix+"photometricSubdivHalfDifference"));
  		    if (properties.getProperty(prefix+"photometricSmoothVarianceSigma")!=null) this.photometricSmoothVarianceSigma=Double.parseDouble(properties.getProperty(prefix+"photometricSmoothVarianceSigma"));
  		    if (properties.getProperty(prefix+"photometricScaleVariance")!=null) this.photometricScaleVariance=Double.parseDouble(properties.getProperty(prefix+"photometricScaleVariance"));

  		    if (properties.getProperty(prefix+"matchStatVarianceBlurScale")!=null) this.matchStatVarianceBlurScale=Double.parseDouble(properties.getProperty(prefix+"matchStatVarianceBlurScale"));
  		    if (properties.getProperty(prefix+"matchStatKLocal")!=null) this.matchStatKLocal=Double.parseDouble(properties.getProperty(prefix+"matchStatKLocal"));
		    if (properties.getProperty(prefix+"matchStatMode")!=null) this.matchStatMode=Integer.parseInt(properties.getProperty(prefix+"matchStatMode"));

  		    if (properties.getProperty(prefix+"debugRow")!=null) this.debugRow=Integer.parseInt(properties.getProperty(prefix+"debugRow"));
  		    if (properties.getProperty(prefix+"debugColumn")!=null) this.debugColumn=Integer.parseInt(properties.getProperty(prefix+"debugColumn"));

    	}
    	public void combineDisparityScales(double [][] actualScales, boolean newPriority){
    		if (actualScales==null) return; // do nothing
    		if ((this.disparityScales==null) || (this.disparityScales.length!=actualScales.length)){
    			this.disparityScales=actualScales.clone(); // only top level
    			return;
    		}
    		for (int i=0;i<this.disparityScales.length;i++){
    			if (((this.disparityScales[i]==null) || newPriority) && (actualScales[i]!=null)){
    				this.disparityScales[i]=actualScales[i];
    			} else if (((actualScales[i]==null) || !newPriority) && (this.disparityScales[i]!=null)){
    				actualScales[i]=this.disparityScales[i];
    			}
    		}
    	}

    	public boolean showDialog(String title, double [][] actualScales, int mode, int fullWidth, int fullHeight){
    		if (actualScales!=null) combineDisparityScales(actualScales, true);
			GenericDialog gd=new GenericDialog("title");
			if ((mode & 1)!=0) gd.addCheckbox    ("Autocorrelation",                   this.autocorrelation); // false
			gd.addNumericField("Correlation area size (power of 2)",this.corrFFTSize,0); // 256
			if ((mode & 16)!=0) gd.addNumericField("Tile overlap fraction",this.tileOverlapFraction,0,2,""); //4
			if ((mode & 1)!=0) gd.addNumericField("Selection area center horizontal/right (XC)",this.corrXC,2,6,"pix"); //1280
			if ((mode & 1)!=0) gd.addNumericField("Selection area center vertical/down  (YC)",this.corrYC,2,6,"pix"); // 960
			if ((mode & 1)!=0) gd.addNumericField("Shift from the center",this.corrShift,2,6,"pix"); // 30.0 - was int 28

			if ((mode & 8)!=0)gd.addNumericField("First image in a pair (0-2)",this.firstImage,0); // 0
			if ((mode & 8)!=0)gd.addNumericField("Second image in a pair (0-2)",this.secondImage,0); // 1
			String [] componentNames={"Y","Cb","Cr","Aux."};
			gd.addNumericField("Use phase correlation (0.0 - pure normal correlation, 1.0 - pure phase one)",this.corrPhaseFraction,3,5,"");// 0.5?
			gd.addNumericField("Weight of Cb (relative to Y) in correlation",this.corrCbWeight,3,5,""); //0.5
			gd.addNumericField("Weight of Cr (relative to Y) in correlation",this.corrCrWeight,3,5,""); //0.5
			gd.addNumericField("correlationHighPassSigma, - pixels in frequency domain",this.correlationHighPassSigma,3,5,""); // 1.5
			gd.addNumericField("correlationLowPassSigma, - fraction of the frequency range",this.correlationLowPassSigma,3,5,""); // 0.6
			gd.addNumericField("Noise normalization sigma for Y-component",this.noiseNormalizationSignaY,2,5,"pix"); // 3.0
			gd.addNumericField("Noise normalization sigma for Cb, Cr components",this.noiseNormalizationSignaCbCr,2,5,"pix"); // 5.0
			gd.addNumericField("Contrast threshold",this.contrastThreshold,1,5,""); //1.5
			gd.addCheckbox    ("Use binary alpha)", this.useBinaryAlpha); // true;
//   	public double [][] disparityScales=       null; // for each channel - a pair of {scaleX, scaleY} or null if undefined (interSensor has the same)
			if (this.disparityScales!=null) {
				double [] disparityZero={0.0,0.0};
				double [][] disparityPair={{-0.5,0.0},{0.5,0.0}};
				double [][] disparityTriplet={{0.0,Math.sqrt(3)/3}, {-0.5,-Math.sqrt(3)/6},{0.5,-Math.sqrt(3)/6}};
				gd.addMessage("Channel shift X,Y per unit of disparity (distance between images in a pair)");
				for (int nChn=0;nChn<this.disparityScales.length;nChn++){
					if (this.disparityScales[nChn]==null){
						if      (this.disparityScales.length==2) this.disparityScales[nChn]=disparityPair[nChn].clone();
						else if (this.disparityScales.length==3) this.disparityScales[nChn]=disparityTriplet[nChn].clone();
						else this.disparityScales[nChn]=disparityZero.clone();
					}
					gd.addMessage("Channel "+nChn+":");
					gd.addNumericField("   Channel "+nChn+" shift X (right)",this.disparityScales[nChn][0],6,8,"x");
					gd.addNumericField("   Channel "+nChn+" shift Y (down)",this.disparityScales[nChn][1],6,8,"x");
				}
			}
			if ((mode & 2)!=0) {
				gd.addNumericField("Disparity range",this.disparityRange,2,6,"pix"); //100.0
				gd.addNumericField("Disparity steps (over the range)",this.disparitySteps,0,3,""); //100
				gd.addNumericField("Size of square used for sub-pixel disparity interpolation near local maximums",this.interpolationSize,0,8,"pix"); //8
				gd.addNumericField("Subdivide pixels fraction near correlation local maximums",this.interpolationUpSample,0,3,""); //16
				gd.addMessage("Full image width: "+fullWidth+" pix");
				gd.addMessage("Full image height: "+fullHeight+" pix");
				gd.addCheckbox    ("Use full image window (overwrite next fields)", this.useFullWindow); // true;
				gd.addNumericField("Result window left margin",this.resultWindow.x,0,4,""); //0
				gd.addNumericField("Result window top  margin",this.resultWindow.y,0,4,""); //0
				gd.addNumericField("Result window width", this.resultWindow.width,0,4,""); //2592
				gd.addNumericField("Result window height",this.resultWindow.height,0,4,""); //1936
			}
			if ((mode & 4)!=0) {
				gd.addCheckbox    ("Process Y  channel", ((this.channelMask & 1)!=0)); // true;
				gd.addCheckbox    ("Process Cb channel", ((this.channelMask & 2)!=0)); // true;
				gd.addCheckbox    ("Process Cr channel", ((this.channelMask & 4)!=0)); // true;
				gd.addCheckbox    ("Use external data file (plus/instead of YCbCr)", this.useFileData); // false;
				gd.addStringField ("External data file path (if selected plus/instead of YCbCr)",this.dataPathName,100);
				gd.addStringField ("Calculated correlation vs. disparity data file",this.corrPathName,100);
			}
			if ((mode & 2)!=0) {
				gd.addNumericField("Absolute correlation maximum threshold for subpixel resolution, Y",this.subpixAMax[0],3,7,"");
				gd.addNumericField("Absolute correlation maximum threshold for subpixel resolution, Cb",this.subpixAMax[1],3,7,"");
				gd.addNumericField("Absolute correlation maximum threshold for subpixel resolution, Cr",this.subpixAMax[2],3,7,"");
				gd.addNumericField("Absolute correlation maximum threshold for subpixel resolution, External file",this.subpixAMax[3],3,7,"");
				gd.addNumericField("Relative correlation maximum threshold for subpixel resolution, all channels",100*this.subpixRMax,2,6,"%");
				gd.addNumericField("Number of correlation maximums to enable subpixel resolution, all channels",this.subpixNMax,0,2,"");
				gd.addCheckbox    ("Show correlation results as a multi-slice image", this.showCorrelationToImage); // true;
				gd.addCheckbox    ("Save correlation results as a multi-slice image", this.saveCorrelationToImage); // true;

			}
			if ((mode & 32)!=0){
				gd.addNumericField("Try background plane not less that this fraction of maximal correlation plane",100*this.bgFraction,1,5,"%");
				gd.addNumericField("Blur variance sigma",this.blurVarianceSigma,2,5,"pix");
				gd.addNumericField("Tile period for correlation refinement (power of 2)",this.refineTilePeriod,0,2,"pix"); // 4
				gd.addNumericField("Refine correlation phase correlation coefficient (0.0 - correlation, 1.0 - phase correlation)",this.refinePhaseCoeff,2,4,"");
				gd.addNumericField("refineHighPassSigma, - pixels in frequency domain",this.refineHighPassSigma,3,5,"pix"); // 1.5
				gd.addNumericField("refineLowPassSigma, - fraction of the frequency range",this.refineLowPassSigma,3,5,""); // 0.6
				gd.addNumericField("Maximal distance from the center to correlation maximum to refine disparity",this.refineCorrMaxDistance,2,5,"pix");
				gd.addNumericField("Relative correlation strength to enable disparity refinement",this.refineCorrThreshold,2,5,"");
				gd.addNumericField("Subdivide pixels to calculate disparity correction (power of 2)",this.refineSubPixel,0,2,"");
			}
			if ((mode & 64)!=0){
				if (this.zMapWOI==null) this.zMapWOI=new Rectangle(0,0,fullWidth,fullHeight);
				gd.addNumericField("Minimal absolute correlation strength for strongest maximum",this.zMapMinFirst,4,6,"");
				gd.addNumericField("Minimal absolute correlation strength for any maximum",this.zMapMinAbsolute,4,6,"");
				gd.addNumericField("Minimal correlation strength as a fraction of the strongest",100*this.zMapMinRelative,2,6,"%");
				gd.addNumericField("Merge maximums with disparity difference less than",this.zMapMergeMax,2,5,"pix");
				gd.addNumericField("Maximal number of z-planes to consider",this.zMapMaxNumber,0,2,"");
				gd.addNumericField("Store overlap padding to z-tiles, on each side",this.zMapOverlap,0,3,"pix");
				gd.addNumericField("Render window left margin",this.zMapWOI.x,0,4,"");
				gd.addNumericField("Result window top  margin",this.zMapWOI.y,0,4,"");
				gd.addNumericField("Result window width", this.zMapWOI.width,0,4,"");
				gd.addNumericField("Result window height",this.zMapWOI.height,0,4,"");
				gd.addNumericField("Fat zero for combining to center",this.zMapFatZero,4,6,"");
			}
			if ((mode & 128)!=0){
				gd.addNumericField("Minimal absolute correlation strength for the foreground",this.zMapMinForeground,4,6,"");
				gd.addMessage(" -------------- filter foreground by tone matching --------------");
				for (int i=0;i<componentNames.length;i++){
					gd.addCheckbox    ("Filter by "+componentNames[i]+" component matching", (this.zMapVarMask & (1<<i))!=0); // true;
					gd.addNumericField("Threshold for "+componentNames[i]+" component relative difference",this.zMapVarThresholds[i],4,6,"");
					gd.addNumericField("Weight of  "+componentNames[i]+" component",this.zMapVarWeights[i],4,6,"");
				}
				gd.addNumericField("Aux channel variance mode (0 - all multiply, 1 - add pairs)",this.auxVarMode,0,1,"");
				gd.addMessage(" -------------- filter foreground by local correlation --------------");
				gd.addNumericField    ("Normalize correlation mode ", this.normalizeCorrelation,0,1,""); // true;
				for (int i=0;i<componentNames.length;i++){
					gd.addCheckbox    ("Filter by "+componentNames[i]+" component correlation", (this.zMapCorrMask & (1<<i))!=0); // true;
					gd.addNumericField("Threshold for "+componentNames[i]+" component correlation (0.0 - just positive, >1.0 - strong, <0 - loosen",this.zMapCorrThresholds[i],4,6,"");
					gd.addNumericField("Relative (to autocorrelation) threshold for "+componentNames[i]+" component correlation (0.0 - just positive, >1.0 - strong, <0 - loosen",100*this.zMapCorrThresholdsRel[i],4,6,"%");
					gd.addNumericField("Weight of  "+componentNames[i]+" component",this.zMapCorrWeights[i],4,6,"");
				}
			}
			if ((mode & 256)!=0){
				gd.addNumericField("Minimal absolute correlation for the new foreground",this.fillFgGapMin,4,6,"");
				gd.addNumericField("Maximal difference between the new foreground disparity and that of the surrounding tiles",this.fillFgGapDiff,2,5,"pix");
				gd.addNumericField("Minimal number of neighboring tiles (of 8) having the similar disparity as the new foreground",this.fillFgGapNeib,0,1,"");
			}
			if ((mode & 512)!=0){
				gd.addMessage(" -------------- filter2 parameters --------------");
				gd.addNumericField("Maximal disparity to process",this.filter2DisparityMax,2,6,"pix");
				gd.addNumericField("Minimal disparity to process",this.filter2disparityMin,2,6,"pix");
				gd.addCheckbox    ("Update abolute/relative maximums", this.filter2UpdateMax); // true;
				gd.addNumericField("Absolute threshold for correlation vs. disparity maximums",this.filter2MinAbsolute,4,5,"");
				gd.addNumericField("Relative threshold for correlation vs. disparity maximums",100*this.filter2MinRelative,2,6,"%");
				gd.addCheckbox    ("Filter by occluding foreground (closer than this) planes", this.filter2ByForeground); // true;
				gd.addNumericField("Minimal disparity difference to closer occluding objects",this.filter2ByForegroundMargin,2,5,"");
				gd.addCheckbox    ("Filter by disabled pixels at the same distance", this.filter2ByDisabled); // true;
				gd.addNumericField("Disparity tolerance to consider objects to be at the same distance",this.filter2DisparityTolearnce,2,5,"");
				gd.addNumericField("Blur enabling masks with this Gaussian sigma ",this.filter2MaskBlurSigma,2,5,"pix");
				gd.addNumericField("HIgh-pass image data before multiplying by the enabling masks",this.filter2corrHighPassSigma,2,5,"pix");
			}

			if ((mode & 1024)!=0){
				gd.addMessage(" -------------- photometric parameters --------------");
				gd.addNumericField("Ignore fraction of high/low pixels",100*this.photometricIgnoreFraction,4,8,"%");
				gd.addNumericField("Number of intensity value columns in the table",this.photometricSubdivAverage,0,4,"samples");
				gd.addNumericField("Half number of intensity difference rows in the table",this.photometricSubdivHalfDifference,0,4,"samples");
				gd.addNumericField("Blur variance between different value samples",this.photometricSmoothVarianceSigma,2,6,"samples");
				gd.addNumericField("Scale measured variance to create result probability data",this.photometricScaleVariance,2,6,"x");
			}
			if ((mode & 2048)!=0){
				gd.addMessage(" -------------- tone-matching statistical parameters --------------");
				gd.addNumericField("Use scaled varinace as a gaussian sigma for a 2-d histogram",this.matchStatVarianceBlurScale,1,4,"");
				gd.addNumericField("Normalize differences by global variance (0.0) or local (1.0)",100.0*this.matchStatKLocal,1,5,"%");
				gd.addNumericField("Likelyhood mode (0 - mul, 1 - add, 2 - min, 3 - max)",this.matchStatMode,0,1,"");

			}
//			/matchStatVarianceBlurScale
			gd.addNumericField("Debug row",this.debugRow,0);
			gd.addNumericField("Debug column",this.debugColumn,0);

			WindowTools.addScrollBars(gd);
			gd.showDialog();
			if (gd.wasCanceled()) return false;
			if ((mode & 1)!=0) this.autocorrelation=gd.getNextBoolean();
			this.corrFFTSize=       (int) gd.getNextNumber();
			if ((mode & 16)!=0) this.tileOverlapFraction=  (int) gd.getNextNumber(); //4
			if ((mode & 1)!=0) this.corrXC=                 gd.getNextNumber();
			if ((mode & 1)!=0) this.corrYC=                 gd.getNextNumber();
			if ((mode & 1)!=0) this.corrShift=              gd.getNextNumber();
			if ((mode & 8)!=0) this.firstImage=       (int) gd.getNextNumber();
			if ((mode & 8)!=0)this.secondImage=       (int) gd.getNextNumber();
			this.corrPhaseFraction=          gd.getNextNumber();
			this.corrCbWeight=               gd.getNextNumber();
			this.corrCrWeight=               gd.getNextNumber();
			this.correlationHighPassSigma=   gd.getNextNumber();
			this.correlationLowPassSigma=    gd.getNextNumber();
			this.noiseNormalizationSignaY=   gd.getNextNumber();
			this.noiseNormalizationSignaCbCr=gd.getNextNumber();
			this.contrastThreshold=gd.getNextNumber();
			this.useBinaryAlpha=gd.getNextBoolean();
    		if (actualScales!=null) combineDisparityScales(actualScales, false);
			if (this.disparityScales!=null) {
				for (int nChn=0;nChn<this.disparityScales.length;nChn++){
					this.disparityScales[nChn][0]=gd.getNextNumber();
					this.disparityScales[nChn][1]=gd.getNextNumber();
				}
			}
			if ((mode & 2)!=0) {
				this.disparityRange=             gd.getNextNumber();
				this.disparitySteps=       (int) gd.getNextNumber();
				this.interpolationSize=    (int) gd.getNextNumber();
				this.interpolationUpSample=(int) gd.getNextNumber();
				this.resultWindow=new Rectangle();
				this.useFullWindow=              gd.getNextBoolean();
				this.resultWindow.x=       (int) gd.getNextNumber();
				this.resultWindow.y=       (int) gd.getNextNumber();
				this.resultWindow.width=   (int) gd.getNextNumber();
				this.resultWindow.height=  (int) gd.getNextNumber();
				if (this.useFullWindow){
					this.resultWindow.x=       0;
					this.resultWindow.y=       0;
					this.resultWindow.width=   fullWidth;
					this.resultWindow.height=  fullHeight;
				}
			}
			if ((mode & 4)!=0) {
				this.channelMask=0;
				if (gd.getNextBoolean()) this.channelMask+=1;
				if (gd.getNextBoolean()) this.channelMask+=2;
				if (gd.getNextBoolean()) this.channelMask+=4;
				this.useFileData=               gd.getNextBoolean();
				this.dataPathName=              gd.getNextString();
				this.corrPathName=              gd.getNextString();
			}
			if ((mode & 2)!=0) {
				this.subpixAMax[0]=             gd.getNextNumber();
				this.subpixAMax[1]=             gd.getNextNumber();
				this.subpixAMax[2]=             gd.getNextNumber();
				this.subpixAMax[3]=             gd.getNextNumber();
				this.subpixRMax=           0.01*gd.getNextNumber();
				this.subpixNMax=          (int) gd.getNextNumber();
				this.showCorrelationToImage=    gd.getNextBoolean();
				this.saveCorrelationToImage=    gd.getNextBoolean();
			}
			if ((mode & 32)!=0){
				this.bgFraction=           0.01*gd.getNextNumber();
				this.blurVarianceSigma=         gd.getNextNumber();
				this.refineTilePeriod=    (int) gd.getNextNumber();
				this.refinePhaseCoeff=          gd.getNextNumber();
				this.refineHighPassSigma=       gd.getNextNumber();
				this.refineLowPassSigma=        gd.getNextNumber();
				this.refineCorrMaxDistance=     gd.getNextNumber();
				this.refineCorrThreshold=       gd.getNextNumber();
				this.refineSubPixel=      (int) gd.getNextNumber();
			}
			if ((mode & 64)!=0){
				this.zMapMinFirst=              gd.getNextNumber();
				this.zMapMinAbsolute=           gd.getNextNumber();
				this.zMapMinRelative=      0.01*gd.getNextNumber();
				this.zMapMergeMax=              gd.getNextNumber();
				this.zMapMaxNumber=       (int) gd.getNextNumber();
				this.zMapOverlap=         (int) gd.getNextNumber();
				this.zMapWOI.x=           (int) gd.getNextNumber();
				this.zMapWOI.y=           (int) gd.getNextNumber();
				this.zMapWOI.width=       (int) gd.getNextNumber();
				this.zMapWOI.height=      (int) gd.getNextNumber();
				this.zMapFatZero=               gd.getNextNumber();
			}
			if ((mode & 128)!=0){
				this.zMapMinForeground=               gd.getNextNumber();
				this.zMapVarMask=0;
				for (int i=0;i<componentNames.length;i++){
					if (gd.getNextBoolean()) this.zMapVarMask |= (1<<i);
					this.zMapVarThresholds[i]=               gd.getNextNumber();
					this.zMapVarWeights[i]=                  gd.getNextNumber();
				}
				this.auxVarMode=          (int) gd.getNextNumber();
				this.normalizeCorrelation=(int) gd.getNextNumber();
				this.zMapCorrMask=0;
				for (int i=0;i<componentNames.length;i++){
					if (gd.getNextBoolean()) this.zMapCorrMask |= (1<<i);
					this.zMapCorrThresholds[i]=              gd.getNextNumber();
					this.zMapCorrThresholdsRel[i]=      0.01*gd.getNextNumber();
					this.zMapCorrWeights[i]=                 gd.getNextNumber();
				}
			}
			if ((mode & 256)!=0){
				this.fillFgGapMin=              gd.getNextNumber();
				this.fillFgGapDiff=             gd.getNextNumber();
				this.fillFgGapNeib=       (int) gd.getNextNumber();
			}
			if ((mode & 512)!=0){
				this.filter2DisparityMax=       gd.getNextNumber();
				this.filter2disparityMin=       gd.getNextNumber();
				this.filter2UpdateMax=          gd.getNextBoolean();
				this.filter2MinAbsolute=        gd.getNextNumber();
				this.filter2MinRelative=   0.01*gd.getNextNumber();
				this.filter2ByForeground=       gd.getNextBoolean();
				this.filter2ByForegroundMargin= gd.getNextNumber();
				this.filter2ByDisabled=         gd.getNextBoolean();
				this.filter2DisparityTolearnce= gd.getNextNumber();
				this.filter2MaskBlurSigma=      gd.getNextNumber();
				this.filter2corrHighPassSigma=  gd.getNextNumber();
			}
			if ((mode & 1024)!=0){
				this.photometricIgnoreFraction=       0.01*gd.getNextNumber();
				this.photometricSubdivAverage=       (int) gd.getNextNumber();
				this.photometricSubdivHalfDifference=(int) gd.getNextNumber();
				this.photometricSmoothVarianceSigma=       gd.getNextNumber();
				this.photometricScaleVariance=             gd.getNextNumber();
			}
			if ((mode & 2048)!=0){
				this.matchStatVarianceBlurScale=           gd.getNextNumber();
				this.matchStatKLocal=                 0.01*gd.getNextNumber();
				this.matchStatMode=                  (int) gd.getNextNumber();
			}

			this.debugRow=                (int) gd.getNextNumber();
			this.debugColumn=             (int) gd.getNextNumber();
			return true;
    	}
	}

	public class LinearFeaturesParameters{
		public boolean extractFeatures=        true;  // "(Re)-extract linear features"
		public boolean ignorePhase=            false; // "Ignore feature phase"
		public boolean preserveDC=             false; // "Preserve DC when reconstructing features"
		public double strengthMode=            0.5;   // "Scale features mode: -2.0 - no scale, -1.0 - phase strength, 0.0 - absolute strength, 1.0 - relative strength"
		public int corrFFTSize=                32;    // "Tile size (power of 2)"
		public int overlapFraction=            8;     // "Overlap as fraction of the tile width (FFT size)"
		public double alphaThreshold=          0.5;   // "Minimal fractional tile to use"
		public boolean useBinaryAlpha=         false; // "Use binary alpha (all >0.0 treat as 1.0)"
		public double correlationHighPassSigma=1.5;   // "correlationHighPassSigma, - pixels in frequency domain"
		public double correlationLowPassSigma= 0.6;   // "correlationLowPassSigma, - fraction of the frequency range"
		public double phaseIntegrationWidth=   3.5;   // "Phase integration width in frequency samples"
		public double resultHighPass=          1.0;   // "Result high-pass filter, frequency samples"
		public int featureFilter=              0;     // "Feature filter 0 - any, +1 white lines, +2 - black lines, + 4 black-to-white in the direction, 8 - white-to-black
		// seems 1 - white, 4 black, 2,8 - step
		public double minRMSFrac=              2.0;   // "Minimal frequency sample value relative to RMS to be considered when looking for the linear phase feature (0.0 - skip test) "
		public double minAbs=                  0.0;   // "Minimal frequency sample absolute value to be considered when looking for the linear phase feature  (0.0 - skip test)"
		public double maxPhaseMismatch=        45.0*Math.PI/180.0;  // "Maximal phase mismatch (between diagonals in a 2x2 sqaure) to use frequaency sample",45.0,3,5,"degrees"
		public double dispertionCost=          0.5;   // "Phase dispersion cost when finding direction (0.0 - only amplitudes) ",0.5
		public boolean calculateStrength=      true;  // "Calculate features strengths"
		public double directionTolerance=      20.0*Math.PI/180.0;  // "Direction tolerance, degrees"
		public double normalDistanceTolerance= 2.0;   // "Normal distance  tolerance, pixels"
		public double tangentialDistanceTolerance= 4.0; // "Tangential distance  tolerance, pixels"
		public double hostsTolerance=          1.0;   //"Hosts cells tolerance (if two neighbors claim same point as it's own, pixels"
		public double directionFracSigma=      0.6;   // "Direction sigma - fraction of direction tolerance"
		public double normalDistanceFracSigma= 0.6;   // "Normal distance  sigma - fraction of normal distance  tolerance"
		public double tangentialDistanceFracSigma= 0.6; // "Tangential distance  sigma - fraction of normal distance  tolerance"
		public int    minMerged=               2;     // "Minimal number of merged (0 - including single with off-cell features, 1 - including same cell no-merge, >1 - merged)"
		public double scaleDistances=          1.0;   //
		public double swapTangentialTolerance= 4.0;   // "Maximal tangential shift while moving features to the closer cells"
		public int    swapSearchRange=         3;     // "Number of cells around feature center to search for the new host cell",3,0,5," cells each direction"
		public boolean multiplyByCellUsage=    true;  // "Multiply by cell usage"
		public double cellUsageShift=          1.0;   // "Cell usage shift"
		public boolean enableDisplayUsage=     false; // "Show cell usage"
		public double displayUsageScale=       0.1;   // "Cell usage display scale (if enabled))"
		public int debugRow=                   -10;   // "Debug row"
		public int debugColumn=                -10;   // "Debug column"

    	public void setProperties(String prefix,Properties properties){
  			properties.setProperty(prefix+"extractFeatures",this.extractFeatures+"");
    		properties.setProperty(prefix+"ignorePhase",this.ignorePhase+"");
    		properties.setProperty(prefix+"preserveDC", this.preserveDC+"");
    		properties.setProperty(prefix+"strengthMode",this.strengthMode+"");
    		properties.setProperty(prefix+"corrFFTSize",this.corrFFTSize+"");
    		properties.setProperty(prefix+"overlapFraction",this.overlapFraction+"");
    		properties.setProperty(prefix+"alphaThreshold",this.alphaThreshold+"");
    		properties.setProperty(prefix+"useBinaryAlpha",this.useBinaryAlpha+"");
    		properties.setProperty(prefix+"correlationHighPassSigma",this.correlationHighPassSigma+"");
    		properties.setProperty(prefix+"correlationLowPassSigma",this.correlationLowPassSigma+"");
    		properties.setProperty(prefix+"phaseIntegrationWidth",this.phaseIntegrationWidth+"");
    		properties.setProperty(prefix+"resultHighPass",this.resultHighPass+"");
    		properties.setProperty(prefix+"featureFilter",this.featureFilter+"");
    		properties.setProperty(prefix+"minRMSFrac",this.minRMSFrac+"");
    		properties.setProperty(prefix+"minAbs",this.minAbs+"");
    		properties.setProperty(prefix+"maxPhaseMismatch",this.maxPhaseMismatch+"");
    		properties.setProperty(prefix+"dispertionCost",this.dispertionCost+"");
    		properties.setProperty(prefix+"calculateStrength",this.calculateStrength+"");
    		properties.setProperty(prefix+"directionTolerance",this.directionTolerance+"");
    		properties.setProperty(prefix+"normalDistanceTolerance",this.normalDistanceTolerance+"");
    		properties.setProperty(prefix+"tangentialDistanceTolerance",this.tangentialDistanceTolerance+"");
    		properties.setProperty(prefix+"hostsTolerance",this.hostsTolerance+"");
    		properties.setProperty(prefix+"directionFracSigma",this.directionFracSigma+"");
    		properties.setProperty(prefix+"normalDistanceFracSigma",this.normalDistanceFracSigma+"");
    		properties.setProperty(prefix+"tangentialDistanceFracSigma",this.tangentialDistanceFracSigma+"");
    		properties.setProperty(prefix+"minMerged",this.minMerged+"");
    		properties.setProperty(prefix+"scaleDistances",this.scaleDistances+"");
    		properties.setProperty(prefix+"swapTangentialTolerance",this.swapTangentialTolerance+"");
    		properties.setProperty(prefix+"swapSearchRange",this.swapSearchRange+"");
    		properties.setProperty(prefix+"multiplyByCellUsage",this.multiplyByCellUsage+"");
    		properties.setProperty(prefix+"cellUsageShift",this.cellUsageShift+"");
    		properties.setProperty(prefix+"enableDisplayUsage",this.enableDisplayUsage+"");
    		properties.setProperty(prefix+"displayUsageScale",this.displayUsageScale+"");
    		properties.setProperty(prefix+"debugRow",this.debugRow+"");
    		properties.setProperty(prefix+"debugColumn",this.debugColumn+"");
    	}

    	public void getProperties(String prefix,Properties properties){
  		    if (properties.getProperty(prefix+"extractFeatures")!=null) this.extractFeatures=Boolean.parseBoolean(properties.getProperty(prefix+"extractFeatures"));
  		    if (properties.getProperty(prefix+"ignorePhase")!=null) this.ignorePhase=Boolean.parseBoolean(properties.getProperty(prefix+"ignorePhase"));
  		    if (properties.getProperty(prefix+"preserveDC")!=null)  this.preserveDC=Boolean.parseBoolean(properties.getProperty(prefix+"preserveDC"));
  		    if (properties.getProperty(prefix+"strengthMode")!=null) this.strengthMode=Double.parseDouble(properties.getProperty(prefix+"strengthMode"));
  		    if (properties.getProperty(prefix+"corrFFTSize")!=null) this.corrFFTSize=Integer.parseInt(properties.getProperty(prefix+"corrFFTSize"));
  		    if (properties.getProperty(prefix+"overlapFraction")!=null) this.overlapFraction=Integer.parseInt(properties.getProperty(prefix+"overlapFraction"));
  		    if (properties.getProperty(prefix+"alphaThreshold")!=null) this.alphaThreshold=Double.parseDouble(properties.getProperty(prefix+"alphaThreshold"));
  		    if (properties.getProperty(prefix+"useBinaryAlpha")!=null) this.useBinaryAlpha=Boolean.parseBoolean(properties.getProperty(prefix+"useBinaryAlpha"));
  		    if (properties.getProperty(prefix+"correlationHighPassSigma")!=null) this.correlationHighPassSigma=Double.parseDouble(properties.getProperty(prefix+"correlationHighPassSigma"));
  		    if (properties.getProperty(prefix+"correlationLowPassSigma")!=null) this.correlationLowPassSigma=Double.parseDouble(properties.getProperty(prefix+"correlationLowPassSigma"));
  		    if (properties.getProperty(prefix+"phaseIntegrationWidth")!=null) this.phaseIntegrationWidth=Double.parseDouble(properties.getProperty(prefix+"phaseIntegrationWidth"));
  		    if (properties.getProperty(prefix+"resultHighPass")!=null) this.resultHighPass=Double.parseDouble(properties.getProperty(prefix+"resultHighPass"));
  		    if (properties.getProperty(prefix+"featureFilter")!=null) this.featureFilter=Integer.parseInt(properties.getProperty(prefix+"featureFilter"));
  		    if (properties.getProperty(prefix+"minRMSFrac")!=null) this.minRMSFrac=Double.parseDouble(properties.getProperty(prefix+"minRMSFrac"));
  		    if (properties.getProperty(prefix+"minAbs")!=null) this.minAbs=Double.parseDouble(properties.getProperty(prefix+"minAbs"));
  		    if (properties.getProperty(prefix+"maxPhaseMismatch")!=null) this.maxPhaseMismatch=Double.parseDouble(properties.getProperty(prefix+"maxPhaseMismatch"));
  		    if (properties.getProperty(prefix+"dispertionCost")!=null) this.dispertionCost=Double.parseDouble(properties.getProperty(prefix+"dispertionCost"));
  		    if (properties.getProperty(prefix+"calculateStrength")!=null) this.calculateStrength=Boolean.parseBoolean(properties.getProperty(prefix+"calculateStrength"));
  		    if (properties.getProperty(prefix+"directionTolerance")!=null) this.directionTolerance=Double.parseDouble(properties.getProperty(prefix+"directionTolerance"));
  		    if (properties.getProperty(prefix+"normalDistanceTolerance")!=null) this.normalDistanceTolerance=Double.parseDouble(properties.getProperty(prefix+"normalDistanceTolerance"));
  		    if (properties.getProperty(prefix+"tangentialDistanceTolerance")!=null) this.tangentialDistanceTolerance=Double.parseDouble(properties.getProperty(prefix+"tangentialDistanceTolerance"));
  		    if (properties.getProperty(prefix+"hostsTolerance")!=null) this.hostsTolerance=Double.parseDouble(properties.getProperty(prefix+"hostsTolerance"));
  		    if (properties.getProperty(prefix+"directionFracSigma")!=null) this.directionFracSigma=Double.parseDouble(properties.getProperty(prefix+"directionFracSigma"));
  		    if (properties.getProperty(prefix+"normalDistanceFracSigma")!=null) this.normalDistanceFracSigma=Double.parseDouble(properties.getProperty(prefix+"normalDistanceFracSigma"));
  		    if (properties.getProperty(prefix+"tangentialDistanceFracSigma")!=null) this.tangentialDistanceFracSigma=Double.parseDouble(properties.getProperty(prefix+"tangentialDistanceFracSigma"));
  		    if (properties.getProperty(prefix+"minMerged")!=null) this.minMerged=Integer.parseInt(properties.getProperty(prefix+"minMerged"));
  		    if (properties.getProperty(prefix+"scaleDistances")!=null) this.scaleDistances=Double.parseDouble(properties.getProperty(prefix+"scaleDistances"));
  		    if (properties.getProperty(prefix+"swapTangentialTolerance")!=null) this.swapTangentialTolerance=Double.parseDouble(properties.getProperty(prefix+"swapTangentialTolerance"));
  		    if (properties.getProperty(prefix+"swapSearchRange")!=null) this.swapSearchRange=Integer.parseInt(properties.getProperty(prefix+"swapSearchRange"));
  		    if (properties.getProperty(prefix+"multiplyByCellUsage")!=null) this.multiplyByCellUsage=Boolean.parseBoolean(properties.getProperty(prefix+"multiplyByCellUsage"));
  		    if (properties.getProperty(prefix+"cellUsageShift")!=null) this.cellUsageShift=Double.parseDouble(properties.getProperty(prefix+"cellUsageShift"));
  		    if (properties.getProperty(prefix+"enableDisplayUsage")!=null) this.enableDisplayUsage=Boolean.parseBoolean(properties.getProperty(prefix+"enableDisplayUsage"));
  		    if (properties.getProperty(prefix+"displayUsageScale")!=null) this.displayUsageScale=Double.parseDouble(properties.getProperty(prefix+"displayUsageScale"));
  		    if (properties.getProperty(prefix+"debugRow")!=null) this.debugRow=Integer.parseInt(properties.getProperty(prefix+"debugRow"));
  		    if (properties.getProperty(prefix+"debugColumn")!=null) this.debugColumn=Integer.parseInt(properties.getProperty(prefix+"debugColumn"));
    	}

		public boolean showDialog(String title) {
			GenericDialog gd = new GenericDialog(title);
			gd.addCheckbox    ("(Re)-extract linear features", this.extractFeatures);
			gd.addCheckbox    ("Ignore feature phase", this.ignorePhase);
			gd.addCheckbox    ("Preserve DC when reconstructing features", this.preserveDC);
			gd.addNumericField("Scale features mode: -2.0 - no scale, -1.0 - phase strength, 0.0 - absolute strength, 1.0 - relative strength",this.strengthMode,3,5,"");
			gd.addNumericField("Tile size (power of 2)",this.corrFFTSize,0);
			gd.addNumericField("Overlap as fraction of the tile width (FFT size)",this.overlapFraction,0);
			gd.addNumericField("Minimal fractional tile to use",100*this.alphaThreshold,1,5,"%");
			gd.addCheckbox    ("Use binary alpha (all >0.0 treat as 1.0)", this.useBinaryAlpha);
			gd.addNumericField("correlationHighPassSigma, - pixels in frequency domain",this.correlationHighPassSigma,3,5,"");
			gd.addNumericField("correlationLowPassSigma, - fraction of the frequency range",this.correlationLowPassSigma,3,5,"");
			gd.addNumericField("Phase integration width ",this.phaseIntegrationWidth,3,5,"frequency samples");
			gd.addNumericField("Result high-pass filter",this.resultHighPass,3,5,"frequency samples");
			gd.addNumericField("Feature filter 0 - any, +1 white lines, +4 - black lines, + 2 black-to-white in the direction, 8 - white-to-black",this.featureFilter,0);
			gd.addNumericField("Minimal frequency sample value relative to RMS to be considered when looking for the linear phase feature (0.0 - skip test) ",this.minRMSFrac,3,5,"");
			gd.addNumericField("Minimal frequency sample absolute value to be considered when looking for the linear phase feature  (0.0 - skip test)",this.minAbs,3,5,"");
			gd.addNumericField("Maximal phase mismatch (between diagonals in a 2x2 sqaure) to use frequaency sample",180.0/Math.PI*this.maxPhaseMismatch,3,5,"degrees");
			gd.addNumericField("Phase dispersion cost when finding direction (0.0 - only amplitudes) ",this.dispertionCost,3,5,"");
			gd.addCheckbox    ("Calculate features strengths", this.calculateStrength);
			gd.addMessage     ("Merging features:");
			gd.addNumericField("Direction tolerance",180.0/Math.PI*this.directionTolerance,3,5,"degrees");
			gd.addNumericField("Normal distance  tolerance",this.normalDistanceTolerance,3,5,"pixels");
			gd.addNumericField("Tangential distance  tolerance",this.tangentialDistanceTolerance,3,5,"pixels");
			gd.addNumericField("Hosts cells tolerance (if two neighbors claim same point as it's own",this.hostsTolerance,3,5,"pixels");
			gd.addNumericField("Direction sigma - fraction of direction tolerance",this.directionFracSigma,3,5,"x");
			gd.addNumericField("Normal distance  sigma - fraction of normal distance  tolerance",this.normalDistanceFracSigma,3,5,"x");
			gd.addNumericField("Tangential distance  sigma - fraction of normal distance  tolerance",this.tangentialDistanceFracSigma,3,5,"x");
			gd.addNumericField("Minimal number of merged (0 - including single with off-cell features, 1 - including same cell no-merge, >1 - merged)",this.minMerged,0);
			gd.addNumericField("Scale distance (compensate for apparent decrease of measured distance to feature caused by windowing)",this.scaleDistances,3,5,"x");
			gd.addNumericField("Maximal tangential shift while moving features to the closer cells",this.swapTangentialTolerance,3,5,"x");
			gd.addNumericField("Number of cells around feature center to search for the new host cell",this.swapSearchRange,0,5," cells each direction");
			gd.addCheckbox    ("Multiply by cell usage", this.multiplyByCellUsage);
			gd.addNumericField("Cell usage shift",this.cellUsageShift,3,5,"x");
			gd.addCheckbox    ("Show cell usage", this.enableDisplayUsage);
			gd.addNumericField("Cell usage display scale (if enabled))",this.displayUsageScale,3,5,"x");
			gd.addNumericField("Debug row",this.debugRow,0);
			gd.addNumericField("Debug column",this.debugColumn,0);
			WindowTools.addScrollBars(gd);
			gd.showDialog();
			if (gd.wasCanceled()) return false;
			this.extractFeatures=        gd.getNextBoolean();
			this.ignorePhase=            gd.getNextBoolean();
			this.preserveDC=             gd.getNextBoolean();

			this.strengthMode=             gd.getNextNumber();
			this.corrFFTSize=          (int) gd.getNextNumber();
			this.overlapFraction=      (int) gd.getNextNumber();
			this.alphaThreshold=     0.01*gd.getNextNumber();
			this.useBinaryAlpha=         gd.getNextBoolean();
			this.correlationHighPassSigma=gd.getNextNumber();
			this.correlationLowPassSigma= gd.getNextNumber();
			this.phaseIntegrationWidth=   gd.getNextNumber();
			this.resultHighPass=          gd.getNextNumber();
			this.featureFilter=        (int) gd.getNextNumber();
			this.minRMSFrac=              gd.getNextNumber();
			this.minAbs=                  gd.getNextNumber();
			this.maxPhaseMismatch=        Math.PI/180.0*gd.getNextNumber();
			this.dispertionCost=          gd.getNextNumber();
			this.calculateStrength=      gd.getNextBoolean();
			this.directionTolerance=      Math.PI/180.0*gd.getNextNumber();
			this.normalDistanceTolerance= gd.getNextNumber();
			this.tangentialDistanceTolerance= gd.getNextNumber();
			this.hostsTolerance=          gd.getNextNumber();
			this.directionFracSigma=      gd.getNextNumber(); // make a fixed fraction of directionTolerance?
			this.normalDistanceFracSigma= gd.getNextNumber(); // make a fixed fraction of distanceTolerance?
			this.tangentialDistanceFracSigma= gd.getNextNumber(); // make a fixed fraction of distanceTolerance?
			this.minMerged=            (int) gd.getNextNumber();
			this.scaleDistances=          gd.getNextNumber();
			this.swapTangentialTolerance= gd.getNextNumber();
			this.swapSearchRange=   (int) gd.getNextNumber();
			this.multiplyByCellUsage=    gd.getNextBoolean();
			this.cellUsageShift=          gd.getNextNumber();
// !!!!!!			if (!multiplyByCellUsage) cellUsageShift=Double.NaN;
			this.enableDisplayUsage=     gd.getNextBoolean();
			this.displayUsageScale=       gd.getNextNumber()*(enableDisplayUsage?1.0:0.0);
			this.debugRow=             (int) gd.getNextNumber();
			this.debugColumn=          (int) gd.getNextNumber();
			return true;
		}
	}
}
