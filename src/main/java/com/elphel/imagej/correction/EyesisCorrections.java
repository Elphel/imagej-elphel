package com.elphel.imagej.correction;
/**
** -----------------------------------------------------------------------------**
** EyesisCorrections.java
**
** Aberration correction for Eyesis4pi
**
**
** Copyright (C) 2012 Elphel, Inc.
**
** -----------------------------------------------------------------------------**
**
**  EyesisCorrections.java is free software: you can redistribute it and/or modify
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

import java.io.IOException;
import java.util.concurrent.atomic.AtomicInteger;

import javax.swing.SwingUtilities;

import com.elphel.imagej.calibration.DebayerScissorsClass;
import com.elphel.imagej.calibration.PixelMapping;
import com.elphel.imagej.cameras.EyesisCorrectionParameters;
import com.elphel.imagej.common.DoubleFHT;
import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.jp4.JP46_Reader_camera;
import com.elphel.imagej.readers.EyesisTiff;

import ij.CompositeImage;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.io.FileInfo;
import ij.io.FileSaver;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;


public class EyesisCorrections {
	public JP46_Reader_camera JP4_INSTANCE=       new JP46_Reader_camera(false);
	ShowDoubleFloatArrays SDFA_INSTANCE=   new ShowDoubleFloatArrays();
	DebayerScissorsClass debayerScissors=null;
    public AtomicInteger stopRequested=null; // 1 - stop now, 2 - when convenient
	public PixelMapping pixelMapping=null;
	public EyesisCorrectionParameters.CorrectionParameters correctionsParameters=null;
	public boolean [] usedChannels;
	public float [][] channelVignettingCorrection=null;
	public int [][][] defectsXY=null; // per each channel: pixel defects coordinates list (starting with worst)
	public double [][] defectsDiff=null; // per each channel: pixel defects value (diff from average of neighbors), matching defectsXY

	public int   [][] channelWidthHeight=null;
	public ImagePlus [] imageNoiseGains=null;
	public String [] sharpKernelPaths=null;
	public String [] smoothKernelPaths=null;
	public int debugLevel;
	public String [] stackColorNames= {"Red","Green","Blue"};
	public int psfSubpixelShouldBe4=4;         // sub-pixel decimation
	public long   startTime=0;

//	public boolean BUG_subchannel=true; // top channel - 1, middle - 0, bottom - 2 (should be -0-1-2)
//	public boolean BUG_subchannel=false; // top channel - 1, middle - 0, bottom - 2 (should be -0-1-2)



	public EyesisCorrections (
			AtomicInteger stopRequested,
			EyesisCorrectionParameters.CorrectionParameters correctionsParameters
			){
		this.correctionsParameters=correctionsParameters;
		this.stopRequested=stopRequested;
	}
	public void setDebug(int debugLevel){
		this.debugLevel=debugLevel;
	}

	public int getNumChannels(){return (this.usedChannels!=null)?this.usedChannels.length:0;}


	// TODO: preserve some data when re-running with new source files
	// FIXME: Make forgiving alien files
	public void initSensorFiles(int debugLevel){
		initSensorFiles(debugLevel, false);
	}
	public void initSensorFiles(int debugLevel, boolean missing_ok){
		this.sharpKernelPaths=null;
		this.smoothKernelPaths=null;
		String [] sensorPaths=correctionsParameters.selectSensorFiles(this.debugLevel);
		this.pixelMapping=new PixelMapping(sensorPaths,debugLevel);
		this.usedChannels= usedChannels(correctionsParameters.getSourcePaths(),missing_ok);
		// TODO: Combine with additional channel map to be able to select single image (of all 3)
		if (correctionsParameters.removeUnusedSensorData){
			for (int nChn=0;nChn< this.usedChannels.length; nChn++) if (!this.usedChannels[nChn]) this.pixelMapping.removeChannel(nChn);
		}
		int numUsedChannels=0;
		for (int nChn=0;nChn< this.usedChannels.length; nChn++) if (this.usedChannels[nChn]) numUsedChannels++;
		if (this.debugLevel>0) {
			String sChannels="";
			for (int nChn=0;nChn< this.usedChannels.length; nChn++) if (this.usedChannels[nChn]) sChannels+=" "+nChn;
			System.out.println ("Number of used channels: "+numUsedChannels+" ("+sChannels+" )");
		}
		createChannelVignetting();
		if ((this.debugLevel>101) && (correctionsParameters.sourcePaths!=null) && (correctionsParameters.sourcePaths.length>0)) {
			testFF(correctionsParameters.sourcePaths[0]);
		}
		if (debugLevel > -2) {
			System.out.println();
		}
	}

	public void initSensorFilesAux(int debugLevel){
//		this.sharpKernelPaths=null;
//		this.smoothKernelPaths=null;
		String [] sensorPaths=correctionsParameters.selectSensorFiles(this.debugLevel);
		this.pixelMapping=new PixelMapping(sensorPaths,debugLevel);
		this.usedChannels= usedChannels(correctionsParameters.getSourcePaths());
		// TODO: Combine with additional channel map to be able to select single image (of all 3)
		if (correctionsParameters.removeUnusedSensorData){
			for (int nChn=0;nChn< this.usedChannels.length; nChn++) if (!this.usedChannels[nChn]) this.pixelMapping.removeChannel(nChn);
		}
		int numUsedChannels=0;
		for (int nChn=0;nChn< this.usedChannels.length; nChn++) if (this.usedChannels[nChn]) numUsedChannels++;
		if (this.debugLevel>0) {
			String sChannels="";
			for (int nChn=0;nChn< this.usedChannels.length; nChn++) if (this.usedChannels[nChn]) sChannels+=" "+nChn;
			System.out.println ("Number of used channels: "+numUsedChannels+" ("+sChannels+" )");
		}
		createChannelVignetting();
		//java.lang.ClassCastException: [B cannot be cast to [F
		//if ((this.debugLevel>1) && (correctionsParameters.sourcePaths!=null) && (correctionsParameters.sourcePaths.length>0)) {
		if ((this.debugLevel>101) && (correctionsParameters.sourcePaths!=null) && (correctionsParameters.sourcePaths.length>0)) {
			testFF(correctionsParameters.sourcePaths[0]);
		}
	}

	public double [] calcReferenceExposures(int debugLevel){
		String [] paths=this.correctionsParameters.getSourcePaths();
		double [] exposures=new double [paths.length];
		if (this.correctionsParameters.exposureCorrectionMode<2){
			for (int nFile=0;nFile<paths.length;nFile++) {
				exposures[nFile]=(this.correctionsParameters.exposureCorrectionMode>0)?this.correctionsParameters.referenceExposure:Double.NaN;
			}
		} else {
			ImagePlus imp; // using that composite image has same exposure
			for (int nFile=0;nFile<paths.length;nFile++){
				if (this.correctionsParameters.isJP4()){
					imp=JP4_INSTANCE.open(
							"", // path,
							paths[nFile],
							"",  //arg - not used in JP46 reader
							true, // un-apply camera color gains
							null, // new window
							false); // do not show
				} else {
					imp=new ImagePlus(paths[nFile]);
					//				  (new JP46_Reader_camera(false)).decodeProperiesFromInfo(imp_src); // decode existent properties from info
					JP4_INSTANCE.decodeProperiesFromInfo(imp); // decode existent properties from info
				}
				if (imp.getProperty("EXPOSURE")!=null){
					exposures[nFile]=Double.parseDouble((String)imp.getProperty("EXPOSURE"));
				} else {
					exposures[nFile]=Double.NaN;
				}
				if (debugLevel>1){
					System.out.println((nFile+1)+": "+paths[nFile]+" - exposure="+exposures[nFile]);
				}
			}
			String [] names=new String [paths.length];
			for (int nFile=0;nFile<paths.length;nFile++) names[nFile]=this.correctionsParameters.getNameFromSourceTiff(paths[nFile]);
			int [] firstImageIndex=new int [paths.length];
			for (int nFile=0;nFile<paths.length;nFile++) {
				firstImageIndex[nFile]=nFile;
				for (int j=0;j<nFile;j++) if (names[j].equals(names[nFile])){
					firstImageIndex[nFile]=j;
					break;
				}
			}
			double [][] minMaxExposure=new double[paths.length][2];
			for (int nFile=0;nFile<paths.length;nFile++) {
				minMaxExposure[nFile][0]=Double.NaN;
				minMaxExposure[nFile][1]=Double.NaN;
			}
			for (int nFile=0;nFile<paths.length;nFile++) if (!Double.isNaN(exposures[nFile])){
				int j=firstImageIndex[nFile];
				if (Double.isNaN(minMaxExposure[j][0]) || (minMaxExposure[j][0]>exposures[nFile])) minMaxExposure[j][0]=exposures[nFile];
				if (Double.isNaN(minMaxExposure[j][1]) || (minMaxExposure[j][1]<exposures[nFile])) minMaxExposure[j][1]=exposures[nFile];
			}
			for (int nFile=0;nFile<paths.length;nFile++) if (!Double.isNaN(exposures[nFile])){
				int j=firstImageIndex[nFile];
				exposures[nFile]=(1.0-this.correctionsParameters.relativeExposure)*minMaxExposure[j][0]+
				this.correctionsParameters.relativeExposure*minMaxExposure[j][1];
			}
		}
		// apply modes
		return exposures;
	}

	public void rebuildEquirectangularMaps(
			EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
			int          threadsMax,  // maximal number of threads to launch
			boolean    updateStatus,
			int        debugLevel){
		this.sharpKernelPaths=null;
		this.smoothKernelPaths=null;
		String [] sensorPaths=correctionsParameters.selectSensorFiles(this.debugLevel);
		String directory= correctionsParameters.selectEquirectangularDirectory(true,true);
		if (directory==null) {
			System.out.println ("No directory selected for equirectangular maps to save");
			return;
		}
//		boolean processPlaneProjection= equirectangularParameters.generateCommonPlane &&
//			equirectangularParameters.selectChannelsToProcess("Select channels for plane projection", this.pixelMapping.sensors.length);

		this.pixelMapping=new PixelMapping(sensorPaths,debugLevel);
		pixelMapping.generateAndSaveEquirectangularMaps(
				correctionsParameters.equirectangularDirectory+
				Prefs.getFileSeparator()+
				correctionsParameters.equirectangularPrefix+
				correctionsParameters.equirectangularSuffix,
				equirectangularParameters.longitudeLeft,
				equirectangularParameters.longitudeRight,
				equirectangularParameters.latitudeTop,
				equirectangularParameters.latitudeBottom,
				equirectangularParameters.pixelsHorizontal,
				equirectangularParameters.imageWidth, //int width,
				equirectangularParameters.imageHeight, //int height,
				equirectangularParameters.x0, //double x0,
				equirectangularParameters.y0, //double y0,
	    		1.0/equirectangularParameters.resolutionScale, //double pixelStep,
	    		equirectangularParameters.longitudeWidth,
	    		equirectangularParameters.clearFullMap,
	    		equirectangularParameters.clearAllMaps,
	    		threadsMax);
		boolean processPlaneProjection= equirectangularParameters.generateCommonPlane &&
		equirectangularParameters.selectChannelsToProcess("Select channels for plane projection", this.pixelMapping.sensors.length);

		if (processPlaneProjection){
//			equirectangularParameters.projectionPixelSize=pixelMapping.degreesPerPixel*Math.PI/180.0;
			boolean [] channelMask= equirectangularParameters.getChannelsToProcess();
			int numUsedChannels=0;
			for (int i=0;i<channelMask.length;i++) if (channelMask[i]) numUsedChannels++;
			int [] channelList=new int [numUsedChannels];
			int iChannel=0;
			for (int i=0;i<channelMask.length;i++) if (channelMask[i]) channelList[iChannel++]=i;
			String sChannels="";
			for (int i=0;i<channelList.length;i++) sChannels+=" "+channelList[i];

			for (int i=0;i<channelList.length;i++) {
				int channel=channelList[i];
				if (!pixelMapping.isChannelAvailable(channel)){
					String msg="No sensor data for channel "+channel;
					System.out.println("Error "+msg);
					IJ.showMessage("Error",msg);
					return;
				}
				if (!pixelMapping.isEquirectangularMapAvailable(channel)){
					String path=correctionsParameters.selectEquirectangularMapFile(
							channel,
							debugLevel);

					if (path==null) {
						String msg="No equirectangular map found for channel "+channel;
						System.out.println("Error "+msg);
						IJ.showMessage("Error",msg);
						return;
					}
					if (debugLevel>1) System.out.println("rebuildEquirectangularMaps(): channel="+channel+" path="+path);

					pixelMapping.loadChannelEquirectangularMap(
							channel,
							path);

					if (!this.pixelMapping.isEquirectangularMapAvailable(channel)){
						String msg="Failed to load equirectangular map for channel "+channel;
						System.out.println("Error "+msg);
						IJ.showMessage("Error",msg);
						return;
					}
				}
			}



			String title="Projection_plane_map";
			ImagePlus imp_pixelmap= pixelMapping.getPlaneToSensorsMap( // need to re-load equirectangular maps?
					channelList, // int [] channels,
					equirectangularParameters.projectionElevation, // Latitude (in degrees) of the normal to the projection plane
					equirectangularParameters.projectionYaw,       // Longitude (in degrees) of the normal to the projection plane
					equirectangularParameters.projectionRoll,      // Rotation of the projection plane around the perpendicular from the lens centers
					equirectangularParameters.matchPixelSize?0.0:	equirectangularParameters.projectionPixelSize, // Ratio of the plane pixel size to the distance from the lens center to the projection plane
					equirectangularParameters.projectionWidth,     // Width of the projection rectangle
					equirectangularParameters.projectionHeight,    // Height of the projection rectangle
					equirectangularParameters.projectionCenterX,   // X-coordinate (along the projection plane X - right) of the intersection of the projection plane with the perpendicular from the lens center
					equirectangularParameters.projectionCenterY,   // Y-coordinate (along the projection plane Y - down) of the intersection of the projection plane with the perpendicular from the lens center
					equirectangularParameters.nominalHorizontalDisparity, // 60.0 - nominal distance between horizontal cameras, mm
		    		title,
		    		debugLevel
		    		);
			if (equirectangularParameters.matchPixelSize) {
				equirectangularParameters.projectionPixelSize=Math.PI/180.0*pixelMapping.panoDegreesPerPixel;
	    		if (debugLevel>0) System.out.println("rebuildEquirectangularMaps(): Setting  equirectangularParameters.projectionPixelSize="+equirectangularParameters.projectionPixelSize);
			}
        	if (imp_pixelmap!=null) {
        		if (debugLevel>2) {
        			imp_pixelmap.getProcessor().resetMinAndMax(); // imp_psf will be reused
        			imp_pixelmap.show();
        		}
        		FileSaver fs=new FileSaver(imp_pixelmap);
        		String resultPath=correctionsParameters.equirectangularDirectory+
				Prefs.getFileSeparator()+correctionsParameters.planeMapPrefix+correctionsParameters.planeMapSuffix;
        		String msg="Saving pixel map to a common plane for sensors "+sChannels+": "+resultPath;
        		if (updateStatus) IJ.showStatus(msg);
        		if (debugLevel>0) System.out.println(msg);
        		fs.saveAsTiffStack(resultPath);
        	} else {
           	 System.out.println("Failed to create pixel map for sensors "+sChannels);
        	}

		}
	}


	public boolean updateImageNoiseGains(
			EyesisCorrectionParameters.NonlinParameters nonlinParameters,
			int          fftSize, // 128 - fft size, kernel size should be size/2
			int          threadsMax,  // maximal number of threads to launch
			boolean    updateStatus,
			int        globalDebugLevel){
		boolean removeUnused=this.correctionsParameters.removeUnusedSensorData;
		int numChannels=this.usedChannels.length;
		if (this.imageNoiseGains==null){
			this.imageNoiseGains= new ImagePlus[0];
		}
		if (this.imageNoiseGains.length!=numChannels){
			ImagePlus [] tmp=this.imageNoiseGains.clone();
			this.imageNoiseGains=new ImagePlus[numChannels];
			for (int chn=0;chn<numChannels;chn++) this.imageNoiseGains[chn]=(chn<tmp.length)?tmp[chn]:null;
		}
		this.sharpKernelPaths=correctionsParameters.selectKernelChannelFiles(
				0,  // 0 - sharp, 1 - smooth
				numChannels, // number of channels
				this.debugLevel);
		if (this.sharpKernelPaths==null) return false;
		if (nonlinParameters.useDiffNoiseGains) {
			this.smoothKernelPaths=correctionsParameters.selectKernelChannelFiles(
					1,  // 0 - sharp, 1 - smooth
					numChannels, // number of channels
					this.debugLevel);
			if (this.smoothKernelPaths==null) return false;
		}
		for (int chn=0;chn<this.usedChannels.length;chn++){
			if (this.usedChannels[chn] && (this.sharpKernelPaths[chn]!=null) && (!nonlinParameters.useDiffNoiseGains ||(this.smoothKernelPaths[chn]!=null))){
				if (
						(this.imageNoiseGains[chn]==null) ||
						(!this.sharpKernelPaths[chn].equals(this.imageNoiseGains[chn].getProperty("sharpKernelPath"))) ||
						(!this.smoothKernelPaths[chn].equals(this.imageNoiseGains[chn].getProperty("smoothKernelPath")))){

					ImagePlus imp_kernel_sharp=new ImagePlus(this.sharpKernelPaths[chn]);
					  if (imp_kernel_sharp.getStackSize()<3) {
						  System.out.println("Need a 3-layer stack with kernels");
						  this.sharpKernelPaths[chn]=null;
						  continue;
					  }
					  ImageStack kernel_sharp_stack= imp_kernel_sharp.getStack();
					  ImageStack kernel_smooth_stack=null;
					  if (nonlinParameters.useDiffNoiseGains) {
							ImagePlus imp_kernel_smooth=new ImagePlus(this.smoothKernelPaths[chn]);
							  if (imp_kernel_smooth.getStackSize()<3) {
								  System.out.println("Need a 3-layer stack with kernels");
								  this.smoothKernelPaths[chn]=null;
								  continue;
							  }
							  kernel_smooth_stack= imp_kernel_smooth.getStack();
					  }
					ImageStack kernelsNoise=
						calculateKernelsNoiseGains (
								kernel_sharp_stack, //final ImageStack kernelStack1, // first stack with 3 colors/slices convolution kernels
								  kernel_smooth_stack, //final ImageStack kernelStack2, // second stack with 3 colors/slices convolution kernels (or null)
								  fftSize, //size, // 128 - fft size, kernel size should be size/2
								  nonlinParameters.blurSigma,
								  threadsMax,  // maximal number of threads to launch
								  updateStatus,
								  globalDebugLevel);
					  kernel_sharp_stack= null; // TODO: - maybe keep one set to speed-up single-channel processing?
					  kernel_smooth_stack=null;
					  Runtime.getRuntime().gc();
			     	  String title="noiseGains_"+(nonlinParameters.useDiffNoiseGains?"diff_":"")+String.format("%02d",chn);
					  imageNoiseGains[chn]= new ImagePlus(title, kernelsNoise);
					  imageNoiseGains[chn].setProperty("sharpKernelPath", this.sharpKernelPaths[chn]);
					  imageNoiseGains[chn].setProperty("smoothKernelPath", nonlinParameters.useDiffNoiseGains?this.smoothKernelPaths[chn]:"");
					  if (this.correctionsParameters.saveNoiseGains || this.correctionsParameters.showNoiseGains) {
						  saveAndShow(this.imageNoiseGains[chn],
								  this.correctionsParameters,
								  this.correctionsParameters.saveNoiseGains,
								  this.correctionsParameters.showNoiseGains
						  );
					  }
				}
			} else {
				if (removeUnused) this.imageNoiseGains[chn]=null;
			}
			if (this.stopRequested.get()>0) {
				System.out.println("User requested stop");
				return false;
			}

		}
		return true;
	}

	public void createChannelVignetting(){
		this.channelWidthHeight=new int [this.usedChannels.length][];
		this.channelVignettingCorrection=new float [this.usedChannels.length][];
		this.defectsXY=new int [this.usedChannels.length][][];
		this.defectsDiff=new double [this.usedChannels.length][];

		for (int nChn=0;nChn< this.usedChannels.length; nChn++){
			this.channelWidthHeight[nChn]=null;
			this.channelVignettingCorrection[nChn]=null;
			this.defectsXY[nChn]=null;
			this.defectsDiff[nChn]=null;
		}
		int [][] bayer={{1,0},{2,1}}; // GR/BG
		ImagePlus imp=null,imp_composite=null;
		for (int nFile=0;nFile<correctionsParameters.getSourcePaths().length;nFile++){
			int [] channels={correctionsParameters.getChannelFromSourceTiff(correctionsParameters.getSourcePaths()[nFile])};
			if (correctionsParameters.isJP4()){
				int subCamera= channels[0]- correctionsParameters.firstSubCamera; // to match those in the sensor files
				channels=this.pixelMapping.channelsForSubCamera(subCamera);
			}

			if (channels!=null) {
				imp=null;
				imp_composite=null;
				if (correctionsParameters.isJP4()){
					  imp_composite=JP4_INSTANCE.open(
							  "", // path,
							  correctionsParameters.getSourcePaths()[nFile],
							  "",  //arg - not used in JP46 reader
							  true, // un-apply camera color gains
							  null, // new window
							  false); // do not show
				} else imp=new ImagePlus(correctionsParameters.getSourcePaths()[nFile]);
				if ((imp==null) && (imp_composite==null)) {
					if (this.debugLevel>0) System.out.println("createChannelVignetting(): can not open "+correctionsParameters.getSourcePaths()[nFile]+
							" as "+(correctionsParameters.isJP4()?"JP4":"TIFF")+" file");
					continue;
				}
				for (int chn=0;chn<channels.length;chn++) {
					int srcChannel=channels[chn];
					if ((this.channelWidthHeight[srcChannel]==null) && this.pixelMapping.isChannelAvailable(srcChannel)){
						int subChannel=this.pixelMapping.getSubChannel(srcChannel);
						  if (this.correctionsParameters.swapSubchannels01) {
							  switch (subChannel){
							  case 0: subChannel=1; break;
							  case 1: subChannel=0; break;
							  }
						  }

						if (subChannel<0){
							System.out.println("BUG in createChannelVignetting(): chn="+chn+
									" subChannel="+subChannel+
									" nFile="+nFile+" : "+correctionsParameters.getSourcePaths()[nFile]+
									" channels.length="+channels.length);
							for (int i=0;i<channels.length;i++) System.out.print(" "+channels[i]);
							System.out.println();
							for (int i=0;i<this.usedChannels.length;i++) if (this.usedChannels[i]) {
								System.out.println(i+": subCamera="+this.pixelMapping.sensors[i].subcamera);
							}

						}
						if (correctionsParameters.isJP4()) imp=JP4_INSTANCE.demuxImage(imp_composite, subChannel);
						if (imp==null) imp=imp_composite; // not a composite image

//						int [] widthHeight={imp.getWidth(),imp.getHeight()};

//						this.channelWidthHeight[srcChannel]=widthHeight;
						this.channelVignettingCorrection[srcChannel]=this.pixelMapping.getBayerFlatFieldFloat(
								srcChannel,
//								this.channelWidthHeight[srcChannel][0],
//								this.channelWidthHeight[srcChannel][1],
								bayer);
						if (this.debugLevel>0){
							System.out.println("Created vignetting info for channel "+srcChannel+
									" subchannel="+subChannel+" ("+
									correctionsParameters.getSourcePaths()[nFile]+")");
							System.out.println("imageWidth= "+this.channelWidthHeight[srcChannel][0]+" imageHeight="+this.channelWidthHeight[srcChannel][1]);
						}
						this.defectsXY[srcChannel]=this.pixelMapping.getDefectsXY(srcChannel);
						this.defectsDiff[srcChannel]=this.pixelMapping.getDefectsDiff(srcChannel);
						if (this.debugLevel>0){
							if (this.defectsXY[srcChannel]==null){
								System.out.println("No pixel defects info is available for channel "+srcChannel);
							} else {
								System.out.println("Extracted "+this.defectsXY[srcChannel].length+" pixel outlayers for channel "+srcChannel+
										" (x:y:difference");
								int numInLine=8;
								for (int i=0;i<this.defectsXY[srcChannel].length;i++){
									System.out.print(this.defectsXY[srcChannel][0]+":"+this.defectsXY[srcChannel][1]);
									if ((this.defectsDiff[srcChannel]!=null) && (this.defectsDiff[srcChannel].length>i)){
										System.out.print(":"+IJ.d2s(this.defectsDiff[srcChannel][i],3)+" ");
									}
									if (((i%numInLine)==(numInLine-1)) || (i == (this.defectsXY[srcChannel].length-1))) System.out.println();
								}
							}
						}
					}
				}
			}
		}
	}

	boolean [] usedChannels(String [] paths){
		return usedChannels(paths, false);
	}

	boolean [] usedChannels(String [] paths, boolean missing_ok){
		if (paths==null) paths=new String[0];
		int numChannels=this.pixelMapping.getNumChannels();
		boolean [] usedChannels=new boolean[numChannels];
		for (int i=0;i<numChannels;i++) usedChannels[i]= false; // this.pixelMapping.isChannelAvailable(i);
		for (int i=0;i<paths.length;i++){
			int srcChannel=correctionsParameters.getChannelFromSourceTiff(paths[i]); // different for JP4
			if (correctionsParameters.isJP4()){
				int subCamera= srcChannel- correctionsParameters.firstSubCamera; // to match those in the sensor files
				int [] channels=this.pixelMapping.channelsForSubCamera(subCamera);
				if (channels!=null) for (int j=0;j<channels.length;j++) usedChannels[channels[j]]=true;
			} else {
				if (!this.pixelMapping.isChannelAvailable(srcChannel)){
					if ((debugLevel>0) && !missing_ok) {
						System.out.println("No sensor data for channel "+srcChannel+", needed for source file "+paths[i]);
					}
				} else usedChannels[srcChannel] = true;
			}
		}
		return usedChannels;
	}

	public void testFF(String path){
		ImagePlus imp=new ImagePlus(path);
		imp.getProcessor().resetMinAndMax(); // imp_psf will be reused
		imp.show();
		int srcChannel=correctionsParameters.getChannelFromSourceTiff(path);
		int [] channels={srcChannel};
		if (correctionsParameters.isJP4()){
			int subCamera= srcChannel- correctionsParameters.firstSubCamera; // to match those in the sensor files
			channels=this.pixelMapping.channelsForSubCamera(subCamera);
		}
		/*		int [][] bayer={{1,0},{2,1}}; // GR/BG
//		double [] corrFF=this.pixelMapping.getBayerFlatField(
		float [] corrFF=this.pixelMapping.getBayerFlatFieldFloat(
				srcChannel,
				imp.getWidth(),
				imp.getHeight(),
				bayer);
				*/
		for (int chn=0;chn<channels.length;chn++){
			srcChannel=channels[chn];
			float [] corrFF=this.channelVignettingCorrection[srcChannel];
			if (corrFF==null) return;
			float [] pixels=(float[]) imp.getProcessor().getPixels();
			double [] pixelsFlat=new double [corrFF.length];
			for (int i=0;i<corrFF.length;i++) pixelsFlat[i]=pixels[i]*corrFF[i];
			SDFA_INSTANCE.showArrays(corrFF, imp.getWidth(), imp.getHeight(), srcChannel+"-FF-correction");
			SDFA_INSTANCE.showArrays(pixelsFlat, imp.getWidth(), imp.getHeight(), srcChannel+"-flat-"+imp.getTitle());
		}
	}

	public boolean isChannelEnabled(int channel){
		return ((channel>=0) && (channel<this.usedChannels.length) && this.usedChannels[channel]);
	}
	public void processChannelImages(
			EyesisCorrectionParameters.SplitParameters         splitParameters,
			EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			EyesisCorrectionParameters.RGBParameters             rgbParameters,
			EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
			int          convolveFFTSize, // 128 - fft size, kernel size should be size/2
			final int          threadsMax,  // maximal number of threads to launch
			final boolean    updateStatus,
			final int        debugLevel){
		  this.startTime=System.nanoTime();
		  String [] sourceFiles=correctionsParameters.getSourcePaths();
		  boolean [] enabledFiles=new boolean[sourceFiles.length];
		  for (int i=0;i<enabledFiles.length;i++) enabledFiles[i]=false;
		  int numFilesToProcess=0;
		  int numImagesToProcess=0;
		  for (int nFile=0;nFile<enabledFiles.length;nFile++){
			  if ((sourceFiles[nFile]!=null) && (sourceFiles[nFile].length()>1)) {
				  int [] channels={correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile])};
				  if (correctionsParameters.isJP4()){
						int subCamera= channels[0]- correctionsParameters.firstSubCamera; // to match those in the sensor files
						channels=this.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (isChannelEnabled(channels[i])){
						  if (!enabledFiles[nFile]) numFilesToProcess++;
						  enabledFiles[nFile]=true;
						  numImagesToProcess++;
					  }
				  }
			  }
		  }
		  if (numFilesToProcess==0){
			  System.out.println("No files to process (of "+sourceFiles.length+")");
			  return;
		  } else {
			  if (debugLevel>0) System.out.println(numFilesToProcess+ " files to process (of "+sourceFiles.length+"), "+numImagesToProcess+" images to process");
		  }
		  double [] referenceExposures=calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		  int [][] fileIndices=new int [numImagesToProcess][2]; // file index, channel number
		  int index=0;
		  for (int nFile=0;nFile<enabledFiles.length;nFile++){
			  if ((sourceFiles[nFile]!=null) && (sourceFiles[nFile].length()>1)) {
				  int [] channels={correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile])};
				  if (correctionsParameters.isJP4()){
						int subCamera= channels[0]- correctionsParameters.firstSubCamera; // to match those in the sensor files
						channels=this.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (isChannelEnabled(channels[i])){
						  fileIndices[index  ][0]=nFile;
						  fileIndices[index++][1]=channels[i];
					  }
				  }
			  }
		  }
		  for (int iImage=0;iImage<fileIndices.length;iImage++){
			  int nFile=fileIndices[iImage][0];
			  ImagePlus imp_src=null;
//			  int srcChannel=correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile]);
			  int srcChannel=fileIndices[iImage][1];
			  if (correctionsParameters.isJP4()){
				  int subchannel=this.pixelMapping.getSubChannel(srcChannel);
				  if (this.correctionsParameters.swapSubchannels01) {
					  switch (subchannel){
					  case 0: subchannel=1; break;
					  case 1: subchannel=0; break;
					  }
				  }
				  if (debugLevel>0) System.out.println("Processing channel "+fileIndices[iImage][1]+" - subchannel "+subchannel+" of "+sourceFiles[nFile]);
				  ImagePlus imp_composite=JP4_INSTANCE.open(
								  "", // path,
								  sourceFiles[nFile],
								  "",  //arg - not used in JP46 reader
								  true, // un-apply camera color gains
								  null, // new window
								  false); // do not show
				  imp_src=JP4_INSTANCE.demuxImage(imp_composite, subchannel);
				  if (imp_src==null) imp_src=imp_composite; // not a composite image

// do we need to add any properties?
			  } else {
				  imp_src=new ImagePlus(sourceFiles[nFile]);
//				  (new JP46_Reader_camera(false)).decodeProperiesFromInfo(imp_src); // decode existent properties from info
				  JP4_INSTANCE.decodeProperiesFromInfo(imp_src); // decode existent properties from info
				  if (debugLevel>0) System.out.println("Processing "+sourceFiles[nFile]);
			  }
			  double scaleExposure=1.0;
			  if (!Double.isNaN(referenceExposures[nFile]) && (imp_src.getProperty("EXPOSURE")!=null)){
				  scaleExposure=referenceExposures[nFile]/Double.parseDouble((String) imp_src.getProperty("EXPOSURE"));
//				  imp_src.setProperty("scaleExposure", scaleExposure); // it may already have channel
				  if (debugLevel>0) System.out.println("Will scale intensity (to compensate for exposure) by  "+scaleExposure);
			  }
			  imp_src.setProperty("name",    correctionsParameters.getNameFromSourceTiff(sourceFiles[nFile]));
			  imp_src.setProperty("channel", srcChannel); // it may already have channel
			  imp_src.setProperty("path",    sourceFiles[nFile]); // it may already have channel
//			  ImagePlus result=processChannelImage( // returns ImagePlus, but it already should be saved/shown
			  processChannelImage( // returns ImagePlus, but it already should be saved/shown
					  imp_src, // should have properties "name"(base for saving results), "channel","path"
					  splitParameters,
					  debayerParameters,
					  nonlinParameters,
					  colorProcParameters,
					  channelGainParameters,
					  rgbParameters,
					  convolveFFTSize, // 128 - fft size, kernel size should be size/2
					  scaleExposure,
					  threadsMax,  // maximal number of threads to launch
					  updateStatus,
					  debugLevel);
// warp result (add support for different color modes)
			  if (this.correctionsParameters.equirectangular){
				   if (equirectangularParameters.clearFullMap) pixelMapping.deleteEquirectangularMapFull(srcChannel); // save memory? //removeUnusedSensorData - no, use equirectangular specific settings
				   if (equirectangularParameters.clearAllMaps) pixelMapping.deleteEquirectangularMapAll(srcChannel); // save memory? //removeUnusedSensorData - no, use equirectangular specific settings
			  }
			  //pixelMapping
			  Runtime.getRuntime().gc();
			  if (debugLevel>0) System.out.println("Processing image "+(iImage+1)+" (of "+fileIndices.length+") finished at "+
					  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
			  if (this.stopRequested.get()>0) {
				  System.out.println("User requested stop");
				  return;
			  }
		  }
	}

	public void saveTiffWithAlpha(
			ImagePlus imp,
			EyesisCorrectionParameters.CorrectionParameters correctionsParameters)
	throws IOException, FormatException, ServiceException, DependencyException{
		int fullWidth,x0;
		boolean cutTile=!this.correctionsParameters.usePlaneProjection && correctionsParameters.equirectangularCut;
		if (cutTile){
			fullWidth=Integer.parseInt((String) imp.getProperty("ImageFullWidth"));
			x0=       Integer.parseInt((String) imp.getProperty("XPosition"));
			cutTile=(x0+imp.getWidth()>fullWidth) ;
		}
		if (cutTile  ) {
			if (this.debugLevel>0) System.out.println("Cutting result image in two parts to prevent roll-over");
			saveTiffWithAlpha(
					cropEquirectangular(imp, false),
					correctionsParameters);
			saveTiffWithAlpha(
					cropEquirectangular(imp, true),
					correctionsParameters);
			return;
		} else {
			String path= correctionsParameters.selectResultsDirectory(
					true,  // smart,
					true);  //newAllowed, // save
			if (path!=null){
				path+=Prefs.getFileSeparator()+imp.getTitle()+".tiff";
	 			 if (this.debugLevel>0) System.out.println("Saving equirectangular result to "+path);
				(new EyesisTiff(correctionsParameters.tiffCompression)).saveTiff(
						imp,
						path,
						correctionsParameters.equirectangularFormat,
						((correctionsParameters.equirectangularFormat==3)?correctionsParameters.outputRangeFP:correctionsParameters.outputRangeInt),
						correctionsParameters.imageJTags,
						debugLevel);
			}
		}
	}

	public ImagePlus cropEquirectangular(
			ImagePlus imp,
			boolean right){ // right - rolled over to the left of the result, adds "-LEFT", "-RIGHT" to the title
		  int fullWidth=Integer.parseInt((String) imp.getProperty("ImageFullWidth"));
		  int x0=       Integer.parseInt((String) imp.getProperty("XPosition"));
		  int height=   imp.getHeight();
		  int width=    imp.getWidth();
		  int dx=(fullWidth-x0);
		  int cropedWidth=right?(imp.getWidth()+-dx):dx;
		  ImagePlus imp_croped;
		  int addX=right?dx:0;
		  if (imp.getType()==ImagePlus.COLOR_RGB) {
			  int [] allPixels= (int []) imp.getProcessor().getPixels();
			  int [] cropedPixels=new int[cropedWidth*height];
			  for (int y=0;y<height;y++) for (int x=0;x< cropedWidth; x++){
				  cropedPixels[y*cropedWidth+x] = allPixels[y*width+(x+addX)];
			  }
			  ColorProcessor colorProcessor=new ColorProcessor(cropedWidth,height);
			  colorProcessor.setPixels(cropedPixels);
			  imp_croped= new ImagePlus(imp.getTitle()+(right?"-RIGHT":"-LEFT"),colorProcessor);
		  } else if (imp.getStackSize()==4) {
			  float [][] allPixels= new float [4][];
			  for (int c=0;c<allPixels.length;c++) allPixels[c]= (float []) imp.getStack().getPixels(c+1);
			  float [][] cropedPixels=new float [4][];
			  for (int c=0;c<cropedPixels.length;c++) cropedPixels[c]= new float [cropedWidth*height];
			  for (int c=0;c<cropedPixels.length;c++) for (int y=0;y<height;y++) for (int x=0;x< cropedWidth; x++){
				  cropedPixels[c][y*cropedWidth+x] = allPixels[c][y*width+(x+addX)];
			  }
			  ImageStack croppedStack=new ImageStack(cropedWidth,height);
			  for (int c=0;c<allPixels.length;c++){
				  croppedStack.addSlice(imp.getStack().getSliceLabel(c+1),  cropedPixels[c]);
			  }
			  imp_croped= new ImagePlus(imp.getTitle()+(right?"-RIGHT":"-LEFT"),croppedStack);
		  } else{
				String msg="cropEquirectangular(): Unsupported image format";
				System.out.println("Error "+msg);
				IJ.showMessage("Error",msg);
				return null;
		  }
		  (new JP46_Reader_camera(false)).copyProperties (imp, imp_croped);
		  if (right) imp_croped.setProperty("XPosition", "0");
		  (new JP46_Reader_camera(false)).encodeProperiesToInfo(imp_croped);
		  return imp_croped;
	}


	public ImagePlus applyEquirectangular(
			int channel,
			ImagePlus imp,
			int threadsMax,
			int debugLevel){
		if (!pixelMapping.isChannelAvailable(channel)){
			String msg="No sensor data for channel "+channel;
			System.out.println("Error "+msg);
			IJ.showMessage("Error",msg);
			return null;
		}
		if (!pixelMapping.isEquirectangularMapAvailable(channel)){
			String path=correctionsParameters.selectEquirectangularMapFile(
					channel,
					debugLevel);

			if (path==null) {
				String msg="No equirectangular map found for channel "+channel;
				System.out.println("Error "+msg);
				IJ.showMessage("Error",msg);
				return null;
			}
			if (debugLevel>1) System.out.println("applyEquirectangular(): channel="+channel+" path="+path);

			pixelMapping.loadChannelEquirectangularMap(
					channel,
					path);

			if (!this.pixelMapping.isEquirectangularMapAvailable(channel)){
				String msg="Failed to load equirectangular map for channel "+channel;
				System.out.println("Error "+msg);
				IJ.showMessage("Error",msg);
				return null;
			}
		}
		// apply warping here
		//		  double sourceImageScale=2.0*this.correctionsParameters.JPEG_scale;
		int sourceImageScale=2; // *this.correctionsParameters.JPEG_scale;
		ImagePlus imp_warped= pixelMapping.resampleToEquirectangular( // will Add "_EQR"
				imp,
				channel,
				sourceImageScale,
				threadsMax);
 		return 	imp_warped;
	}

	public ImagePlus applyCommonPlane(
			int channel,
			ImagePlus imp,
			int threadsMax,
			int debugLevel){
		if (!pixelMapping.isChannelAvailable(channel)){
			String msg="No sensor data for channel "+channel;
			System.out.println("Error "+msg);
			IJ.showMessage("Error",msg);
			return null;
		}
		if (!pixelMapping.isPlaneMapMapAvailable(channel)){
			String path=correctionsParameters.selectPlaneMapFile(
//					channel,
					debugLevel);

			if (path==null) {
				String msg="No common plane projection map found";
				System.out.println("Error "+msg);
				IJ.showMessage("Error",msg);
				return null;
			}
			if (debugLevel>1) System.out.println("applyCommonPlane(): channel="+channel+" path="+path);
			pixelMapping.loadPlaneMap(
//					channel,
					path,
					debugLevel);

			if (!this.pixelMapping.isPlaneMapMapAvailable(channel)){
				String msg="Failed to load a common plane projection map for channel "+channel+", or that file does not have this sensor data";
				System.out.println("Error "+msg);
				IJ.showMessage("Error",msg);
				return null;
			}
		}
		// apply warping here
		//		  double sourceImageScale=2.0*this.correctionsParameters.JPEG_scale;
		int sourceImageScale=2; // *this.correctionsParameters.JPEG_scale;
		ImagePlus imp_warped= pixelMapping.applyPlaneMap(
	    		imp, //ImagePlus impSrc,
	    		channel,
	    		sourceImageScale, // 2.0
	    		threadsMax,
	    		debugLevel
	    		);
 		return 	imp_warped;
	}

	public int correctDefects(
			ImagePlus imp,
			int channel,
			int debugLevel){
		int numApplied=0;
		if (this.correctionsParameters.pixelDefects && (this.defectsXY!=null)&& (this.defectsXY[channel]!=null)){
			// apply pixel correction
			float [] pixels=(float []) imp.getProcessor().getPixels();
			int width=imp.getWidth();
			int height=pixels.length/width;
			int [] dirsRB={2,2*width,-2,-2*width};
			int [] dirsG={width+1,width-1,-width-1,-width+1};
			for (int i=0;i<this.defectsXY[channel].length;i++){
				if ( // difference provided and is smaller than threshold
						(this.defectsDiff != null) &&
						(this.defectsDiff[channel]!=null) &&
						(this.defectsDiff[channel].length>i) &&
						(Math.abs(this.defectsDiff[channel][i])<this.correctionsParameters.pixelDefectsThreshold)) break;
				int x=this.defectsXY[channel][i][0];
				int y=this.defectsXY[channel][i][1];
				int index=x+y*width;
				int [] dirs=(((x^y)&1)==0)?dirsG:dirsRB;
				// do not bother to correct border pixels
				if ((x<2) || (y<2) || (x>(width-3)) || (y>height-3)) continue;
				double s=0.0;
				for (int dir=0;dir<dirs.length;dir++) s+=pixels[index+dirs[dir]];
				pixels[index]=(float) (s/dirs.length);
				numApplied++;
			}
		}
		return numApplied;
	}


	public ImagePlus processChannelImage(
			ImagePlus imp_src, // should have properties "name"(base for saving results), "channel","path"
			EyesisCorrectionParameters.SplitParameters         splitParameters,
			EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			EyesisCorrectionParameters.RGBParameters             rgbParameters,
			int          convolveFFTSize, // 128 - fft size, kernel size should be size/2
			double 		     scaleExposure,
			final int        threadsMax,  // maximal number of threads to launch
			final boolean    updateStatus,
			final int        debugLevel){
		boolean advanced=this.correctionsParameters.zcorrect || this.correctionsParameters.equirectangular;
		boolean crop=      advanced? true: this.correctionsParameters.crop;
		boolean rotate=    advanced? false: this.correctionsParameters.rotate;
		double JPEG_scale= advanced? 1.0: this.correctionsParameters.JPEG_scale;
		boolean toRGB=     advanced? true: this.correctionsParameters.toRGB;

		// may use this.StartTime to report intermediate steps execution times
		String name=(String) imp_src.getProperty("name");
		//		int channel= Integer.parseInt((String) imp_src.getProperty("channel"));
		int channel= (Integer) imp_src.getProperty("channel");
		String path= (String) imp_src.getProperty("path");
		if (this.correctionsParameters.pixelDefects && (this.defectsXY!=null)&& (this.defectsXY[channel]!=null)){
			// apply pixel correction
			int numApplied=	correctDefects(
					imp_src,
					channel,
					debugLevel);
			if ((debugLevel>0) && (numApplied>0)) { // reduce verbosity after verified defect correction works
				System.out.println("Corrected "+numApplied+" pixels in "+path);
			}
		}
		if (this.correctionsParameters.vignetting){
			if ((this.channelVignettingCorrection==null) || (channel<0) || (channel>=this.channelVignettingCorrection.length) || (this.channelVignettingCorrection[channel]==null)){
				System.out.println("No vignetting data for channel "+channel);
				return null;
			}
			float [] pixels=(float []) imp_src.getProcessor().getPixels();
			if (pixels.length!=this.channelVignettingCorrection[channel].length){
				System.out.println("Vignetting data for channel "+channel+" has "+this.channelVignettingCorrection[channel].length+" pixels, image "+path+" has "+pixels.length);
				return null;
			}
			for (int i=0;i<pixels.length;i++){
				pixels[i]*=this.channelVignettingCorrection[channel][i];
			}
		}
		String title=name+"-"+String.format("%02d", channel);
		ImagePlus result=imp_src;
		if (debugLevel>1) System.out.println("processing: "+path);
		result.setTitle(title+"RAW");
		if (!this.correctionsParameters.split){
			saveAndShow(result, this.correctionsParameters);
			return result;
		}
		// Split into Bayer components, oversample, increase canvas
		ImageStack stack= bayerToStack(
				result, // source Bayer image, linearized, 32-bit (float))
				splitParameters);
		String titleFull=title+"-SPLIT";
		if (!this.correctionsParameters.debayer) {
			result= new ImagePlus(titleFull, stack);
			saveAndShow(result, this.correctionsParameters);
			return result;
		}
		// Demosaic image
		if (debayerScissors==null) debayerScissors=new DebayerScissorsClass(this.stopRequested);
		debayerScissors.setDebug(debugLevel);
		stack= debayerScissors.aliasScissorsStack(stack,  // stack with 3 colors/slices with the image
				debayerParameters,
				(this.correctionsParameters.saveDebayerEnergy || this.correctionsParameters.showDebayerEnergy),
				threadsMax, // number of image pixels/ sensor pixels (each direction) == 2
				updateStatus,// update status info
				debugLevel);
		if (this.correctionsParameters.saveDebayerEnergy || this.correctionsParameters.showDebayerEnergy) {
			if (debayerScissors.getDebayerEnergy()!=null) {
				ImagePlus debayerMask=SDFA_INSTANCE.makeArrays (debayerScissors.getDebayerEnergy(),
						debayerScissors.getDebayerEnergyWidth(),
						debayerScissors.getDebayerEnergy().length/debayerScissors.getDebayerEnergyWidth(),
						title+"-DEBAYER-ENERGY");
				saveAndShow(debayerMask,
						this.correctionsParameters,
						this.correctionsParameters.saveDebayerEnergy,
						this.correctionsParameters.showDebayerEnergy
				);
			}
		}
		titleFull=title+"-DEMOSAIC";
		CorrectionDenoise correctionDenoise=new CorrectionDenoise(stopRequested);
		CorrectionColorProc correctionColorProc=new CorrectionColorProc(stopRequested);
		result= new ImagePlus(titleFull, stack);
		if (this.correctionsParameters.deconvolve) {
			//Ask for the kernel directory if it is undefined
			if (this.sharpKernelPaths==null){ // make sure the paths list is reset after changing parameters
				this.sharpKernelPaths=correctionsParameters.selectKernelChannelFiles(
						0,  // 0 - sharp, 1 - smooth
						this.usedChannels.length, // number of channels
						this.debugLevel);
			}
			if ((this.sharpKernelPaths==null) || (this.sharpKernelPaths[channel]==null)){
				System.out.println("Sharp kernel path does not exist");
				return null;
			}
			// Read deconvolution kernels
			ImagePlus imp_sharp_kernels=new ImagePlus(this.sharpKernelPaths[channel]);
			if (imp_sharp_kernels.getStackSize()<3) {
				System.out.println("Need a 3-layer stack with kernels - file "+this.sharpKernelPaths[channel]);
				return null;
			}
			ImageStack convolutionSharpKernelStack=imp_sharp_kernels.getStack();
			if (debugLevel>1) System.out.println("Using kernel stack "+this.sharpKernelPaths[channel]+" for convolution with "+result.getTitle());
			ImageStack stackDeconvolvedSharp= convolveStackWithKernelStack( //  stack_d
					stack,  // stack with 3 colors/slices with the image
					convolutionSharpKernelStack, // stack with 3 colors/slices convolution kernels
					convolveFFTSize, // 128 - fft size, kernel size should be size/2
					threadsMax,
					updateStatus, // update status info
					debugLevel);
			imp_sharp_kernels=null; // free memory
			convolutionSharpKernelStack=null;
			Runtime.getRuntime().gc();
			titleFull=title+"-DECONV";
			if (this.correctionsParameters.combine) {
				// Read "smooth" kernels
				if (this.smoothKernelPaths==null){ // make sure the paths list is reset after changing parameters
					this.smoothKernelPaths=correctionsParameters.selectKernelChannelFiles(
							1,  // 0 - sharp, 1 - smooth
							this.usedChannels.length, // number of channels
							this.debugLevel);
				}
				if ((this.smoothKernelPaths==null) || (this.smoothKernelPaths[channel]==null)){
					System.out.println("Smooth kernel path does not exist");
					return null;
				}
				ImagePlus imp_smooth_kernels=new ImagePlus(this.smoothKernelPaths[channel]);
				if (imp_smooth_kernels.getStackSize()<3) {
					System.out.println("Need a 3-layer stack with kernels - file "+this.smoothKernelPaths[channel]);
					return null;
				}
				ImageStack convolutionSmoothKernelStack=imp_smooth_kernels.getStack();
				if (debugLevel>1) System.out.println("Using smooth kernel stack "+this.smoothKernelPaths[channel]+" for convolution with "+result.getTitle());
				ImageStack stackDeconvolvedSmooth = convolveStackWithKernelStack( //stack_g
						stack,  // stack with 3 colors/slices with the image
						convolutionSmoothKernelStack, // stack with 3 colors/slices convolution kernels
						convolveFFTSize, // 128 - fft size, kernel size should be size/2
						threadsMax,
						updateStatus, // update status info
						debugLevel);
				imp_smooth_kernels=null; // free memory
				convolutionSmoothKernelStack=null;
				Runtime.getRuntime().gc();
				// Combine Smooth and Sharp images
				double [][] noiseMask= extractNoiseMask(
						this.imageNoiseGains[channel],// contains 3-slice stack (r,b,g)
						nonlinParameters.noiseGainWeights[0], // coefficient for slice 0 (r)
						nonlinParameters.noiseGainWeights[1], // coefficient for slice 1 (b)
						nonlinParameters.noiseGainWeights[2], // coefficient for slice 2 (g)
						1,     // decimate result (not yet supported)
						nonlinParameters.noiseGainPower
				);

				// show noise mask here?
				nonlinParameters.showMask=this.correctionsParameters.showDenoiseMask;
				//		          if (DEBUG_LEVEL>1) System.out.println ( " noiseMask.length="+((noiseMask==null)?"null":(noiseMask.length+" noiseMask[0].length="+noiseMask[0].length)));
//				CorrectionDenoise correctionDenoise=new CorrectionDenoise(stopRequested);
				correctionDenoise.setDebug(debugLevel); // not yet used
				stack=  correctionDenoise.combineLoHiStacks(
						stackDeconvolvedSharp, // ImageStack with the image, convolved with the reversed PSF (sharp but with high noise)
						stackDeconvolvedSmooth,  // ImageStack with the image, convolved with the Gaussian (just lateral compensated)  (blurred, but low noise)
						channel,
						nonlinParameters, // show mask generated and used
						noiseMask, // 2-d array of kernelsNoiseGain (divide mask by it)
						32,        // linear pixels per noiseMask pixels (32)
						threadsMax,
						updateStatus, // update status info
						debugLevel);
				if (this.correctionsParameters.saveDenoiseMask || this.correctionsParameters.showDenoiseMask) {
					ImagePlus denoiseMask=SDFA_INSTANCE.makeArrays (
							correctionDenoise.getDenoiseMask(),
							correctionDenoise.getDenoiseMaskWidth(),
							correctionDenoise.getDenoiseMask().length/correctionDenoise.getDenoiseMaskWidth(),
							title+"-MASK");
					if (this.correctionsParameters.jpeg) {
						//crop Mask to original image size
						if (this.correctionsParameters.crop){
							denoiseMask=cropImage32(denoiseMask,splitParameters);
						}
						//rotate the result
						if (this.correctionsParameters.rotate){
							denoiseMask=rotateImage32CW(denoiseMask);
						}
						//scale the result
						if (this.correctionsParameters.JPEG_scale!=1.0){
							ImageProcessor ip=denoiseMask.getProcessor();
							ip.setInterpolationMethod(ImageProcessor.BICUBIC);
							ip=ip.resize((int)(ip.getWidth()*this.correctionsParameters.JPEG_scale),(int) (ip.getHeight()*this.correctionsParameters.JPEG_scale));
							denoiseMask= new ImagePlus(denoiseMask.getTitle(),ip);
							denoiseMask.updateAndDraw();
						}
						if (this.correctionsParameters.showDenoiseMask) denoiseMask.show();
						//public ImagePlus Image32toGreyRGB24(ImagePlus imp);
						if (this.correctionsParameters.saveDenoiseMask) {
							ImagePlus denoiseMaskRGB24=Image32toGreyRGB24(denoiseMask);
							saveAndShow(denoiseMaskRGB24,
									this.correctionsParameters,
									this.correctionsParameters.saveDenoiseMask,
									false, //processParameters.showDenoiseMask,
									this.correctionsParameters.JPEG_quality);
							denoiseMaskRGB24=null;
						}
					} else {
						saveAndShow(denoiseMask,
								this.correctionsParameters,
								this.correctionsParameters.saveDenoiseMask,
								this.correctionsParameters.showDenoiseMask
						);
					}
				}

			} else { // end of if (this.correctionsParameters.combine)
				stack=stackDeconvolvedSharp;
			} // end of else if (this.correctionsParameters.combine)
		}  else if (this.correctionsParameters.combine) { // "combine" w/o "deconvolve" - just use convolution with smooth kernels
			// Read smooth kernels
			// Read "smooth" kernels
			if (this.smoothKernelPaths==null){ // make sure the paths list is reset after changing parameters
				this.smoothKernelPaths=correctionsParameters.selectKernelChannelFiles(
						1,  // 0 - sharp, 1 - smooth
						this.usedChannels.length, // number of channels
						this.debugLevel);
			}
			if ((this.smoothKernelPaths==null) || (this.smoothKernelPaths[channel]==null)){
				System.out.println("Smooth kernel path does not exist");
				return null;
			}
			ImagePlus imp_smooth_kernels=new ImagePlus(this.smoothKernelPaths[channel]);
			if (imp_smooth_kernels.getStackSize()<3) {
				System.out.println("Need a 3-layer stack with kernels - file "+this.smoothKernelPaths[channel]);
				return null;
			}
			ImageStack convolutionSmoothKernelStack=imp_smooth_kernels.getStack();
			if (debugLevel>1) System.out.println("Using smooth kernel stack "+this.smoothKernelPaths[channel]+" for convolution with "+result.getTitle());
			ImageStack stackDeconvolvedSmooth = convolveStackWithKernelStack( // stack_g
					stack,  // stack with 3 colors/slices with the image
					convolutionSmoothKernelStack, // stack with 3 colors/slices convolution kernels
					convolveFFTSize, // 128 - fft size, kernel size should be size/2
					threadsMax,
					updateStatus, // update status info
					debugLevel);
			imp_smooth_kernels=null; // free memory
			stack=stackDeconvolvedSmooth;
			convolutionSmoothKernelStack=null;
			Runtime.getRuntime().gc();
			titleFull=title+"-LOWRES";
		}// end of if (this.correctionsParameters.deconvolve)
		//stack now has the result, titleFull - correct title for the image
		  if (!this.correctionsParameters.colorProc){
			  result= new ImagePlus(titleFull, stack);
			  saveAndShow(
					  result,
					  this.correctionsParameters);
			  return result;
		  }
		  //Processing colors - changing stack sequence to r-g-b (was r-b-g)
		  if (!fixSliceSequence(
				  stack,
				  debugLevel)){
			  if (debugLevel>0) System.out.println("fixSliceSequence() returned false");
			  return null;
		  }
//		  if (debugLevel>2){
		  if (debugLevel>1){
			  ImagePlus imp_dbg=new ImagePlus(imp_src.getTitle()+"-"+channel+"-preColors",stack);
			  saveAndShow(
					  imp_dbg,
					  this.correctionsParameters);
		  }

		  correctionColorProc.processColorsWeights(stack,
//				  255.0/this.psfSubpixelShouldBe4/this.psfSubpixelShouldBe4, //  double scale,     // initial maximal pixel value (16))
				  255.0/this.psfSubpixelShouldBe4/this.psfSubpixelShouldBe4/scaleExposure, //  double scale,     // initial maximal pixel value (16))
				  colorProcParameters,
				  channelGainParameters,
				  channel,
				  correctionDenoise.getDenoiseMask(),
				  this.correctionsParameters.blueProc,
				  debugLevel);
		  if (debugLevel>1) System.out.println("Processed colors to YPbPr, total number of slices="+stack.getSize());
		  if (debugLevel>2){
			  ImagePlus imp_dbg=new ImagePlus("procColors",stack);
			  saveAndShow(
					  imp_dbg,
					  this.correctionsParameters);
		  }

		// Show/save color denoise mask
		  if ((this.correctionsParameters.saveChromaDenoiseMask || this.correctionsParameters.showChromaDenoiseMask) && (correctionColorProc.getDenoiseMaskChroma()!=null)) {
			  ImagePlus chromaDenoiseMask=SDFA_INSTANCE.makeArrays (correctionColorProc.getDenoiseMaskChroma(),
					  correctionColorProc.getDenoiseMaskChromaWidth(),
					  correctionColorProc.getDenoiseMaskChroma().length/correctionColorProc.getDenoiseMaskChromaWidth(),
					  title+"-MASK_CHROMA");
			  if (this.correctionsParameters.jpeg) {
//crop Mask to original image size
				  if (this.correctionsParameters.crop){
					  chromaDenoiseMask=cropImage32(chromaDenoiseMask,splitParameters);
				  }
//rotate the result
				  if (this.correctionsParameters.rotate){
					  chromaDenoiseMask=rotateImage32CW(chromaDenoiseMask);
				  }
//scale the result
				  if (this.correctionsParameters.JPEG_scale!=1.0){
					  ImageProcessor ip=chromaDenoiseMask.getProcessor();
					  ip.setInterpolationMethod(ImageProcessor.BICUBIC);
					  ip=ip.resize((int)(ip.getWidth()*this.correctionsParameters.JPEG_scale),(int) (ip.getHeight()*this.correctionsParameters.JPEG_scale));
					  chromaDenoiseMask= new ImagePlus(chromaDenoiseMask.getTitle(),ip);
					  chromaDenoiseMask.updateAndDraw();
				  }
				  if (this.correctionsParameters.showChromaDenoiseMask) chromaDenoiseMask.show();
//public ImagePlus Image32toGreyRGB24(ImagePlus imp);
				  if (this.correctionsParameters.saveChromaDenoiseMask) {
					  ImagePlus chromaDenoiseMaskRGB24=Image32toGreyRGB24(chromaDenoiseMask);
					  saveAndShow(chromaDenoiseMaskRGB24,
							  this.correctionsParameters,
							  this.correctionsParameters.saveChromaDenoiseMask,
							  false, //processParameters.showChromaDenoiseMask,
							  this.correctionsParameters.JPEG_quality);
					  chromaDenoiseMaskRGB24=null;
				  }
			  } else {
				  saveAndShow(chromaDenoiseMask,
						  this.correctionsParameters,
						  this.correctionsParameters.saveChromaDenoiseMask,
						  this.correctionsParameters.showChromaDenoiseMask
				  );
			  }
		  }
		  if (toRGB) {
			  correctionColorProc.YPrPbToRGB(stack,
					  colorProcParameters.kr,        // 0.299;
					  colorProcParameters.kb,        // 0.114;
					  colorProcParameters.useFirstY?9:8,        //  int sliceY,
							  6, // int slicePr,
							  7// int slicePb
			  );
			  title=titleFull; // including "-DECONV" or "-COMBO"
			  titleFull=title+"-RGB-float";
			  //Trim stack to just first 3 slices
			  if (debugLevel>2){
				  ImagePlus imp_dbg=new ImagePlus("YPrPbToRGB",stack);
				  saveAndShow(
						  imp_dbg,
						  this.correctionsParameters);
			  }

			  while (stack.getSize()>3) stack.deleteLastSlice();
			  if (debugLevel>1) System.out.println("Trimming color stack");
		  } else {
			  title=titleFull; // including "-DECONV" or "-COMBO"
			  titleFull=title+"-YPrPb"; // including "-DECONV" or "-COMBO"
			  if (debugLevel>1) System.out.println("Using full stack, including YPbPr");
		  }
		  result= new ImagePlus(titleFull, stack);
		  // Crop image to match original one (scaled to oversampling)
		  if (crop){ // always crop if equirectangular
			  stack=cropStack32(stack,splitParameters);
			  if (debugLevel>2){
				  ImagePlus imp_dbg=new ImagePlus("cropped",stack);
				  saveAndShow(
						  imp_dbg,
						  this.correctionsParameters);
			  }

		  }
		  // rotate the result
		  if (rotate){ // never rotate for equirectangular
			  stack=rotateStack32CW(stack);
		  }
		  if (!toRGB && !this.correctionsParameters.jpeg){ // toRGB set for equirectangular
			  saveAndShow(result, this.correctionsParameters);
			  return result;
		  } else { // that's not the end result, save if required
			  saveAndShow(result, this.correctionsParameters, this.correctionsParameters.save32, false,this.correctionsParameters.JPEG_quality); // save, no show
		  }
		  // convert to RGB48 (16 bits per color component)
		  ImagePlus imp_RGB;
		  if (this.correctionsParameters.equirectangularFormat==0){
			  stack=convertRGB32toRGB16Stack(
					  stack,
					  rgbParameters);

			  titleFull=title+"-RGB48";
			  result= new ImagePlus(titleFull, stack);
//			  ImagePlus imp_RGB24;
			  CompositeImage compositeImage=convertToComposite(result);
			  if (!this.correctionsParameters.jpeg && !advanced){ // RGB48 was the end result
				  saveAndShow(compositeImage, this.correctionsParameters);
				  return result;
			  } else { // that's not the end result, save if required
				  saveAndShow(compositeImage, this.correctionsParameters, this.correctionsParameters.save16, false); // save, no show
			  }
			  imp_RGB=convertRGB48toRGB24(
					  stack,
					  title+"-RGB24",
					  0, 65536, // r range 0->0, 65536->256
					  0, 65536, // g range
					  0, 65536,// b range
					  0, 65536);// alpha range
			  if (JPEG_scale!=1.0){
				  ImageProcessor ip=imp_RGB.getProcessor();
				  ip.setInterpolationMethod(ImageProcessor.BICUBIC);
				  ip=ip.resize((int)(ip.getWidth()*JPEG_scale),(int) (ip.getHeight()*JPEG_scale));
				  imp_RGB= new ImagePlus(imp_RGB.getTitle(),ip);
				  imp_RGB.updateAndDraw();
			  }
		  } else {
			  switch (this.correctionsParameters.equirectangularFormat){
			   case 0:
				   titleFull=title+"-INT8";
				   break;
			   case 1:
				   titleFull=title+"-INT16";
				   break;
			   case 2:
				   titleFull=title+"-INT32";
				   break;
			   case 3:
				   titleFull=title+"-FLOAT32";
				   break;
			   case 4:
				   titleFull=title+"-IJSTACK";
				   break;
			  }
			  result= new ImagePlus(titleFull, stack);
			  imp_RGB=result;
		  }
		  if (advanced){
			  if (this.correctionsParameters.zcorrect){ // save double resolution image for distance measurements
				  saveAndShow(
						  imp_RGB,
						  this.correctionsParameters,
						  true,  // save
						  false, // show
						  0);    // force Tiff
			  }
			  if (this.correctionsParameters.equirectangular) {
				  ImagePlus impWarped=null;
				  if (this.correctionsParameters.usePlaneProjection){
					  impWarped=	applyCommonPlane(
							  channel,
							  imp_RGB,
							  threadsMax,
							  debugLevel);

				  } else {
					  impWarped=	applyEquirectangular(
							  channel,
							  imp_RGB,
							  threadsMax,
							  debugLevel);
				  }
				  if ((impWarped.getFileInfo().fileType== FileInfo.RGB) && // only for 8-bit RGB (this.correctionsParameters.equirectangularFormat==0)
						  this.correctionsParameters.planeAsJPEG &&
						  this.correctionsParameters.usePlaneProjection){ // equirectangular - always TIFF
					  saveAndShow(
							  impWarped,
							  this.correctionsParameters,
							  this.correctionsParameters.save,
							  this.correctionsParameters.show,
							  this.correctionsParameters.JPEG_quality);

				  } else {
					  if (this.correctionsParameters.equirectangularFormat<4){
						  try {
							  saveTiffWithAlpha(impWarped,this.correctionsParameters);
						  } catch (IOException e) {
							  // TODO Auto-generated catch block
							  e.printStackTrace();
						  } catch (FormatException e) {
							  // TODO Auto-generated catch block
							  e.printStackTrace();
						  } catch (ServiceException e) {
							  // TODO Auto-generated catch block
							  e.printStackTrace();
						  } catch (DependencyException e) {
							  // TODO Auto-generated catch block
							  e.printStackTrace();
						  }
					  } else { // temporarily
						  saveAndShow(impWarped, this.correctionsParameters);
					  }
				  }
			  }
		  } else {
			  saveAndShow(imp_RGB, this.correctionsParameters);
		  }
		  return result;
	}

	 /* ======================================================================== */


//	  private boolean fixSliceSequence (
	public boolean fixSliceSequence ( // for EyesisDCT
			  ImageStack stack,
			  int debugLevel){
		  int i,j;
		  int [] rgbNumbers= {0,0,0};
		  for (j=0;j<3;j++) {
			  for (i=1;i<=3;i++) if (stack.getSliceLabel(i).toLowerCase().equals(this.stackColorNames[j].toLowerCase())){
	// fix case (capitalized)
//				  System.out.println ( "stack.getSliceLabel("+i+")="+stack.getSliceLabel(i));
//				  System.out.println ( "stackColorNames["+j+"]="+stackColorNames[j]);
				  stack.setSliceLabel(this.stackColorNames[j],i);
				  rgbNumbers[j]=i;
//				  System.out.println ( "rgbNumbers["+j+"]="+rgbNumbers[j]);
				  break;
			  }
		  }
		  if (debugLevel > 2) {
			  System.out.println ( "Input file color slice numbers:");
			  System.out.println ( "  Red -   slice "+((rgbNumbers[0]>0)?rgbNumbers[0]:"missing"));
			  System.out.println ( "  Green - slice "+((rgbNumbers[1]>0)?rgbNumbers[1]:"missing"));
			  System.out.println ( "  Blue -  slice "+((rgbNumbers[2]>0)?rgbNumbers[2]:"missing"));
		  }

		  for (i=0;i<3;i++) if (rgbNumbers[i]<=0) {
			  System.out.println ( this.stackColorNames[i]+ "  slice is missing in the input file. Please check slice names");
			  return false;
		  }
		  while ((rgbNumbers[0]!=1) || (rgbNumbers[1]!=2) ||(rgbNumbers[2]!=3)) {
			  if      (rgbNumbers[0]==1) swapStackSlices(stack,2,3);
			  else if (rgbNumbers[2]==3) swapStackSlices(stack,1,2);
			  else                       swapStackSlices(stack,1,3);
			  for (j=0;j<3;j++) {
				  for (i=1;i<=3;i++) if (stack.getSliceLabel(i).equals(this.stackColorNames[j])){
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
			  EyesisCorrectionParameters.SplitParameters splitParameters) {
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
			  EyesisCorrectionParameters.SplitParameters splitParameters) {
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
		  if (imp.isComposite()) return null;
		  if (imp.getNChannels()>1) {
			  return null; // number of channels should be just 1
		  }
		  int c = imp.getStackSize();
		  imp.setDimensions(c, 1, 1);
		  CompositeImage ci = new CompositeImage(imp, CompositeImage.COMPOSITE);
		  return ci;
	  }

	/* ======================================================================== */
	  public ImageStack convertRGB32toRGB16Stack(
			  ImageStack stack32,
			  EyesisCorrectionParameters.RGBParameters rgbParameters) {
		  ImageStack stack16 = new ImageStack(stack32.getWidth(), stack32.getHeight());
		  int length=stack32.getWidth()*stack32.getHeight();
		  int i,j;
		  float [] fpixels;
		  short [] spixels;
		  int numSlices  = stack32.getSize();
		  if (numSlices > 4) numSlices = 4;
		  if (numSlices < 3) return null;
		  double [] mins= {rgbParameters.r_min, rgbParameters.g_min, rgbParameters.b_min, rgbParameters.alpha_min};
		  double [] maxs= {rgbParameters.r_max, rgbParameters.g_max, rgbParameters.b_max, rgbParameters.alpha_max};
		  double value;
		  double scale;
		  for (i = 0; i < numSlices; i++) {
			 fpixels= (float[])stack32.getPixels(i + 1);
			 scale=65535.0 / (maxs[i] - mins[i]);
			 spixels=new short [length];
			 for (j = 0; j < length; j++) {
				 value=(fpixels[j]-mins[i])*scale;
				 if      (value<0.0)     value = 0.0;
				 else if (value>65535.0) value = 65535.0;
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
		  int [] sliceSeq = new int [numSlices];
		  for (int j = 0; j < numSlices; j++) sliceSeq[j] = (j + ((numSlices > 3)? 3:0)) % 4;

		  int [] pixels=new int[length];
		  int c,d;

		  double [] scale=new double[numSlices];
		  for (c = 0; c < numSlices; c++) {
			  scale[c]=256.0/(maxs[c]-mins[c]);
		  }
		  for (i = 0; i < numSlices; i++) spixels[i]= (short[])stack16.getPixels(i+1);
		  for (i = 0; i < length; i++) {
			  pixels[i]=0;
			  for (int j=0; j < numSlices; j++) {
				  c = sliceSeq[j];
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

	  public ImageStack convertRGB48toRGBA24Stack(
			  ImageStack stack16,
			  double [] dalpha, // alpha pixel array 0..1.0 or null
//			  String title,
			  int r_min,
			  int r_max,
			  int g_min,
			  int g_max,
			  int b_min,
			  int b_max){
		  ImageStack stack8 = new ImageStack(stack16.getWidth(), stack16.getHeight());
		  int [] mins= {r_min,g_min,b_min};
		  int [] maxs= {r_max,g_max,b_max};
	      int i;
		  int length=stack16.getWidth()*stack16.getHeight();
		  short [][] spixels=new short[3][];
		  byte [][] bpixels=new byte [(dalpha != null)? 4: 3][length];
		  int c,d;
		  double [] scale=new double[3];
		  for (c=0;c<3;c++) {
			  scale[c]=256.0/(maxs[c]-mins[c]);
		  }
		  for (i = 0; i <3; i++) spixels[i]= (short[])stack16.getPixels(i+1);
		  for (c = 0; c <3; c++) {
			  for (i = 0; i < length; i++){
				  d=(int)(((spixels[c][i]& 0xffff)-mins[c])*scale[c]);
				  if (d>255) d=255;
				  else if (d<0) d=0;
				  bpixels[c][i]= (byte) d;
			  }
			  stack8.addSlice(stack16.getSliceLabel(c+1), bpixels[c]);
		  }
		  if (dalpha != null) {
			  for (i = 0; i < length; i++){
				  bpixels[3][i] = (byte) (255*dalpha[i]);
			  }
			  stack8.addSlice("alpha", bpixels[3]);
		  }
		  return stack8;
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
    public double [] getSlidingMask(int size) { // duplicate with DebayerScissors
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
	  /* convolve image stack with the kernel stack using FHT. kernels should be (size/2)*(size/2) - currently 64x64, then image will be split into same
	      (size/2)*(size/2) overlapping by step=size/4 segments. Both are zero-padded to size x size, so after convolution the result will not roll over, and
	      processed 128x128 result arrays are accumulated in the output stack.
	      The input image should be properly extended by size/4 in each direction (and so the kernel arrays should match it) - that would minimize border effects.*/

	  /* ======================================================================== */
	  public ImageStack convolveStackWithKernelStack (
			  final ImageStack  imageStack,  // stack with 3 colors/slices with the image
			  final ImageStack kernelStack, // stack with 3 colors/slices convolution kernels
			  final int               size, // 128 - fft size, kernel size should be size/2
			  final int          threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus, // update status info
			  final int globalDebugLevel)
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

		  int ichn,indx,dx,dy,tx,ty,li;
		  final int [] nonOverlapSeq = new int[numberOfKernels];
		  int [] nextFirstFindex=new int[16*nChn];
		  indx = 0;
		  li=0;
		  for (ichn=0;ichn<nChn;ichn++) for (dy=0;dy<4;dy++) for (dx=0;dx<4;dx++) {
			  for (ty=dy; ty < tilesY; ty+=4) for (tx=dx; tx < tilesX; tx+=4){
				  nonOverlapSeq[indx++] = ichn*numberOfKernelsInChn+ ty*tilesX + tx;
			  }
			  nextFirstFindex[li++] = indx;
		  }
		  final AtomicInteger aStopIndex = new AtomicInteger(0);
		  final AtomicInteger tilesFinishedAtomic = new AtomicInteger(1); // first finished will be 1

		  if (globalDebugLevel>1)
			  System.out.println("Eyesis_Corrections:convolveStackWithKernelStack :\n"+
				  "globalDebugLevel="+globalDebugLevel+"\n"+
				  "imgWidth="+imgWidth+"\n"+
				  "imgHeight="+imgHeight+"\n"+
				  "tilesX="+tilesX+"\n"+
				  "tilesY="+tilesY+"\n"+
				  "nChn="+nChn+"\n"+
				  "step="+step+"\n"+
				  "size="+size+"\n"+
				  "kernelSize="+kernelSize+"\n"+
				  "kernelWidth="+kernelWidth+"\n"+
				  "kernelNumHor="+kernelNumHor+"\n"+
				  "numberOfKernelsInChn="+numberOfKernelsInChn+"\n");

		  if (updateStatus) IJ.showStatus("Convolving image with kernels, "+nChn+" channels, "+tilesY+" rows");
		  final long startTime = System.nanoTime();
		  for (li = 0; li < nextFirstFindex.length; li++){
			  aStopIndex.set(nextFirstFindex[li]);
			  if (li>0){
				  ai.set(nextFirstFindex[li-1]);
			  }
//			  System.out.println("\n=== nextFirstFindex["+li+"] =" + nextFirstFindex[li]+" === ");
			  for (int ithread = 0; ithread < threads.length; ithread++) {
				  threads[ithread] = new Thread() {
					  @Override
					public void run() {
						  float [] pixels=null;       // will be initialized at first use
						  float [] kernelPixels=null; // will be initialized at first use
						  double [] kernel=       new double[kernelSize*kernelSize];
						  double [] inTile=       new double[kernelSize*kernelSize];
						  double [] outTile=      new double[size * size];
						  double [] doubleKernel= new double[size * size];
						  int chn,tileY,tileX;
						  int chn0=-1;
						  //					  double debug_sum;
						  //					  int i;
						  DoubleFHT fht_instance =new DoubleFHT(); // provide DoubleFHT instance to save on initializations (or null)
						  //					  for (int nTile0 = ai.getAndIncrement(); nTile0 < numberOfKernels; nTile0 = ai.getAndIncrement()) {
						  for (int nTile0 = ai.getAndIncrement(); nTile0 < aStopIndex.get(); nTile0 = ai.getAndIncrement()) {
							  //aStopIndex
							  int nTile = nonOverlapSeq[nTile0];
							  chn=nTile/numberOfKernelsInChn;
							  tileY =(nTile % numberOfKernelsInChn)/tilesX;
							  tileX = nTile % tilesX;
							  if (tileX < 4) {
	  							  int trow=(tileY+ ((tileY & 3) * tilesY))/4;
								  if (updateStatus) IJ.showStatus("Convolving image with kernels, channel "+(chn+1)+" of "+nChn+", row "+(trow+1)+" of "+tilesY);
								  if (globalDebugLevel>2) System.out.println("Processing kernels, channel "+(chn+1)+" of "+nChn+", row "+(trow+1)+" of "+tilesY+" : "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
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
							  //						  debug_sum=0;
							  //						  for (i=0;i<doubleKernel.length;i++) debug_sum+=doubleKernel[i];
							  //						  if (globalDebugLevel>1) System.out.println("kernel sum="+debug_sum);

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
							   * takes just small fraction of several FHTs, it should be OK - reasonable number of threads will spread and not "stay in line"
							   */

							  // Add smart synchronization - wait only if is too far ahead. First test - no synchronization at all

							  //accumulateSquareTile(
//							  System.out.print(tileY+":"+tileX+"/"+chn+"("+nTile0+"/"+nTile+") ");
//							  if (tileX < 4)System.out.println();
							  nonSyncAccumulateSquareTile(
									  outPixels[chn], //  float pixels array to accumulate tile
									  outTile, // data to accumulate to the pixels array
									  imgWidth, // width of pixels array
									  (tileX-1)*step, // left corner X
									  (tileY-1)*step); // top corner Y
							  final int numFinished=tilesFinishedAtomic.getAndIncrement();
							  if (numFinished % (numberOfKernels/100+1) == 0) {
								  SwingUtilities.invokeLater(new Runnable() {
									  @Override
									public void run() {
										  IJ.showProgress(numFinished,numberOfKernels);
									  }
								  });
							  }

							  //numberOfKernels
						  }
					  }
				  };
			  }
			  startAndJoin(threads);
		  }
		  if (updateStatus) IJ.showStatus("Convolution DONE");
		  if (globalDebugLevel > 1) System.out.println("Threads done in "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		  IJ.showProgress(1.0);
		  /* prepare result stack to return */
		  ImageStack outStack=new ImageStack(imgWidth,imgHeight);
		  for (i=0;i<nChn;i++) {
			  outStack.addSlice(imageStack.getSliceLabel(i+1), outPixels[i]);
		  }
		  if (globalDebugLevel > 0) System.out.println("Convolution done in "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		  return outStack;
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



// duplicates with DebayerScissors

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

	     void  nonSyncAccumulateSquareTile(
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

// end of duplicates with DebayerScissors
/* Convert source Bayer pattern (GR/BG) image to higher resolution, add margins by duplicating pattern around */
	  public ImageStack  bayerToStack(ImagePlus imp, // source bayer image, linearized, 32-bit (float))
			  EyesisCorrectionParameters.SplitParameters splitParameters){

	    if (imp==null) return null;
//	    String [] chnNames={"red","blue","green"};
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

	  // return 2-d double array instead of the stack
	  public double [][]bayerToDoubleStack(ImagePlus imp, // source bayer image, linearized, 32-bit (float))
			  EyesisCorrectionParameters.SplitParameters splitParameters){ // if null - no margins, no oversample
		  if (imp==null) return null;
		  boolean adv = splitParameters != null;
		  int oversample = adv? splitParameters.oversample : 1;
		  int addTop=      adv?splitParameters.addTop:       0;
		  int addLeft=     adv?splitParameters.addLeft:      0;
		  int addBottom=   adv?splitParameters.addBottom:    0;
		  int addRight=    adv?splitParameters.addRight:     0;
		  String [] chnNames={"Red","Blue","Green"}; //Different sequence than RGB!!
		  int nChn=chnNames.length;
		  ImageProcessor ip=imp.getProcessor();
		  int inWidth=imp.getWidth();
		  int inHeight=imp.getHeight();
		  int outHeight=inHeight*oversample + addTop +  addBottom;
		  int outWidth= inWidth* oversample + addLeft + addRight;
		  int outLength=outWidth*outHeight;

		  double [][] outPixels=new double[nChn][outLength];
		  float [] pixels = (float[]) ip.getPixels();
		  int chn,y,x,i,index;
		  int bayerPeriod=2*oversample;
		  int ovrWidth=   inWidth * oversample;
		  int ovrHeight= inHeight *  oversample;
		  for (chn=0;chn<nChn;chn++) for (i=0;i<outPixels[chn].length;i++) outPixels[chn][i]=0.0f;
		  /* Can be optimized - now it calculate input address for all those 0-es */
		  for (index=0; index<outLength; index++) {
			  y=(index / outWidth) - addTop;
			  x=(index % outWidth) - addLeft;
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
		  return outPixels;
	  }



//double []  DENOISE_MASK=null;


// TODO: do similar for JP4, using "subcamera" to "use" all channels for it
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
			  final boolean    updateStatus,
			  final int        globalDebugLevel) // update status info
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
				  @Override
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
							  if (globalDebugLevel>2) System.out.println("Processing kernels, channel "+(chn+1)+" of "+nChn+", row "+(tileY+1)+" of "+kernelNumVert+" : "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
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
//						  System.out.println("Processing kernels, channel "+(chn+1)+" of "+nChn+", row "+(tileY+1)+" of "+kernelNumVert+" sum="+sum);
					  }
				  }
			  };
		  }
		  startAndJoin(threads);
		  if (globalDebugLevel > 1) System.out.println("Threads done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		  /* prepare result stack to return */
		  ImageStack outStack=new ImageStack(kernelNumHor,kernelNumVert);
		  for (i=0;i<nChn;i++) {
			  outStack.addSlice(kernelStack1.getSliceLabel(i+1), outPixles[i]);
		  }
		  return outStack;
	  }

	  public void extractOneKernel(float [] pixels, //  array of combined square kernels, each
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
//	   private void saveAndShow(
	   public void saveAndShow(
			   ImagePlus             imp,
			   EyesisCorrectionParameters.CorrectionParameters  correctionsParameters){
		   saveAndShowEnable( imp,    correctionsParameters , true, true);
	   }

	   public void saveAndShowEnable(
			   ImagePlus             imp,
			   EyesisCorrectionParameters.CorrectionParameters  correctionsParameters,
			   boolean               enableSave,
			   boolean               enableShow){
		   saveAndShow(
				   imp,
				   correctionsParameters,
				   correctionsParameters.save && enableSave,
				   correctionsParameters.show && enableShow,
				   correctionsParameters.JPEG_quality);
	   }

	   public void saveAndShow(
//	   public void saveAndShow(
			   ImagePlus             imp,
			   EyesisCorrectionParameters.CorrectionParameters  correctionsParameters,
			   boolean               save,
			   boolean               show){
		   saveAndShow(imp, correctionsParameters,  save,  show, -1);
	   }

	   public void saveAndShow(
			   ImagePlus             imp,
			   EyesisCorrectionParameters.CorrectionParameters  correctionsParameters,
			   boolean               save,
			   boolean               show,
			   int                   jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
		   String path=null;
		   if (save) path= correctionsParameters.selectResultsDirectory(
				   true,  // smart,
				   true);  //newAllowed, // save

		   saveAndShow(
				   imp,
				   path,
				   correctionsParameters.png,
				   show,
				   jpegQuality);
	   }

	   public void saveAndShow(
			   ImagePlus             imp,
			   String                path,
			   boolean               png,
			   boolean               show,
			   int                   jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
		   saveAndShow(
				   imp,
				   path,
				   png,
				   show,
				   jpegQuality,//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
				   this.debugLevel);
	   }


	   public void saveAndShow(
			   ImagePlus             imp,
			   String                path,
			   boolean               png,
			   boolean               show,
			   int                   jpegQuality,//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
			   int                   debugLevel){
		   if (path!=null) {
			   path+=Prefs.getFileSeparator()+imp.getTitle();
			   boolean hasAlphaHighByte = false;
			   if ((imp.getStackSize()==1) && (imp.getFileInfo().fileType== FileInfo.RGB)) {
				   int [] pixels = (int []) imp.getProcessor().getPixels();
				   for (int i = 0; i < pixels.length; i++){
					   if ((pixels[i] & 0xff000000) != 0){
						   hasAlphaHighByte = true;
						   break;
					   }
				   }
			   }

			   if (hasAlphaHighByte && (jpegQuality <= 0)){
				   if (png){
					   if (debugLevel > 0) System.out.println("Saving RGBA result to "+path+".png");
					   (new EyesisTiff()).savePNG_ARGB32(
							   imp,
							   path+".png"
							   );

				   } else {
					   if (debugLevel > 0) System.out.println("Saving RGBA result to "+path+".tiff");
					   try {
						   (new EyesisTiff()).saveTiffARGB32(
								   imp,
								   path+".tiff",
								   false, // correctionsParameters.imageJTags,
								   debugLevel);
					   } catch (IOException e) {
						   e.printStackTrace();
					   } catch (FormatException e) {
						   e.printStackTrace();
					   } catch (ServiceException e) {
						   e.printStackTrace();
					   } catch (DependencyException e) {
						   e.printStackTrace();
					   }
				   }

			   } else if (((imp.getStackSize()==1)) && (jpegQuality!=0) && ((imp.getFileInfo().fileType== FileInfo.RGB) || (jpegQuality>0))) {
				   if (debugLevel>0) System.out.println("Saving result to "+path+".jpeg");
				   FileSaver fs=new FileSaver(imp);
				   if (jpegQuality>0) FileSaver.setJpegQuality(jpegQuality);
				   fs.saveAsJpeg(path+".jpeg");
			   } else {
				   if (debugLevel>0) System.out.println("Saving result to "+path+".tiff");
				   FileSaver fs=new FileSaver(imp);
				   // Save as RGBA if it is RGBA
				   int bytesPerPixel = imp.getBytesPerPixel();
				   int mode = 0;
				   if (bytesPerPixel > 1) mode = 1;
				   if (bytesPerPixel > 3) mode = 3;// float, mode 2 not supported

				   if ((imp.getStackSize() == 4) && imp.getStack().getSliceLabel(4).equals("alpha")){
					   try {
						   (new EyesisTiff()).saveTiff(
								   imp,
								   path+".tiff",
								   mode, //
								   1.0, // full scale, absolute
								   false, // correctionsParameters.imageJTags,
								   debugLevel);
					   } catch (IOException e) {
						   e.printStackTrace();
					   } catch (FormatException e) {
						   e.printStackTrace();
					   } catch (ServiceException e) {
						   e.printStackTrace();
					   } catch (DependencyException e) {
						   e.printStackTrace();
					   }
				   } else if (imp.getStackSize()>1)  fs.saveAsTiffStack(path+".tiff");
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
	 		  EyesisCorrectionParameters.CorrectionParameters correctionsParameters,
	 		  boolean               save,
	 		  boolean               show){
	 	  String path=null;
	 	  if (save)  path= correctionsParameters.selectResultsDirectory(
	 		    				true,  // smart,
	 		    				true);  //newAllowed, // save

	 	  if (path!=null) {
	 		  path+=Prefs.getFileSeparator()+compositeImage.getTitle();
	 		  if (debugLevel>0) System.out.println("Saving result to "+path+".tiff");
	 		  FileSaver fs=new FileSaver(compositeImage);
	 		  if (compositeImage.getStackSize()>1)  fs.saveAsTiffStack(path+".tiff");
	 		  else fs.saveAsTiff(path+".tiff");
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


}
