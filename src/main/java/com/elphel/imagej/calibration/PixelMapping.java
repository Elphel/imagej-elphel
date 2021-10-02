package com.elphel.imagej.calibration;
/**
** -----------------------------------------------------------------------------**
** PixelMapping.java
**
** Using camera calibration files to resample images for equirectangular projection
**
**
** Copyright (C) 2012 Elphel, Inc.
**
** -----------------------------------------------------------------------------**
**
**  PixelMapping.java is free software: you can redistribute it and/or modify
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
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import javax.swing.SwingUtilities;

import com.elphel.imagej.cameras.EyesisCorrectionParameters;
import com.elphel.imagej.common.DoubleFHT;
import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.jp4.JP46_Reader_camera;
import com.elphel.imagej.readers.EyesisTiff;

import Jama.Matrix;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
//import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.io.Opener;
import ij.process.ColorProcessor;
import ij.text.TextWindow;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;


public class PixelMapping {
	public SensorData [] sensors;
	public int panoWidth=0;
	public int panoHeight=0;
	public double panoDegreesPerPixel=Double.NaN; // needs to be initialized
	public double panoLongitudeLeft=Double.NaN;
	public double panoLongitudeRight=Double.NaN;
	public double panoLatitudeTop=Double.NaN;
	public double panoLatitudeBottom=Double.NaN;

	public int debugLevel=1;
	public int lanczosA=3;
	public int oversampled=2;
	public int binsPerHalfPixel=50;
	public double [][][][] lanczos=null;
	public int maxSensors=100;
	InterSensor lastUsedInterSensor=null;
//	public enum cellTypes {EMPTY_CELL, BAD_CELL, OLD_CELL,NEW_CELL}


    public PixelMapping (String defaultPath, int debugLevel){
    	this.debugLevel=debugLevel;
		String [] extensions={".calib-tiff"};
		String [] defaultPaths=((defaultPath==null) || defaultPath.equals(""))?null:(new String[1]);
		if (defaultPaths!=null) defaultPaths[0]=defaultPath;
		MultipleExtensionsFileFilter parFilter = new MultipleExtensionsFileFilter("",extensions,"distortion calibration .calib-tiff files");
		String [] calibFiles=CalibrationFileManagement.selectFiles(false,
				"Select camera sub-modules calibration files",
				"Select",
				parFilter,
				defaultPaths); // String [] defaultPaths);
       	if ((calibFiles==null) || (calibFiles.length==0)) {
    		IJ.showMessage("No files selected");
       	}
       	int maxChannel=0;
       	for (int i=0;i<calibFiles.length;i++){
			int indexPeriod=calibFiles[i].indexOf('.',calibFiles[i].lastIndexOf(Prefs.getFileSeparator()));
	    	int channel=Integer.parseInt(calibFiles[i].substring(indexPeriod-2,indexPeriod));
      		if (channel>maxChannel) maxChannel=channel;
       	}
       	this.sensors=new SensorData[maxChannel+1];
       	for (int i=0;i<calibFiles.length;i++){
			int indexPeriod=calibFiles[i].indexOf('.',calibFiles[i].lastIndexOf(Prefs.getFileSeparator()));
	    	int channel=Integer.parseInt(calibFiles[i].substring(indexPeriod-2,indexPeriod));
	    	String channelPath=calibFiles[i].substring(0,indexPeriod-2)+String.format("%02d",channel)+calibFiles[i].substring(indexPeriod);
	    	this.sensors[channel]=new SensorData (channelPath, this.debugLevel);
       	}
    }

    public PixelMapping (
    		String [] calibFiles,
    		int first_channel, // 0 - old way
    		int num_channels,  // 0 - any
    		boolean update_channel, // false (replace file channel with effective channel (subtract first_channel)
    		int debugLevel){
    	this.debugLevel=debugLevel;
    	if (calibFiles==null) calibFiles=new String[0];
       	this.sensors=new SensorData[(num_channels > 0)? num_channels: this.maxSensors];
       	for (int i=0;i<this.sensors.length;i++) this.sensors[i]=null;
       	int maxChannel=0;
       	for (int i=0;i<calibFiles.length;i++){
       		SensorData sensorData=new SensorData (calibFiles[i],this.debugLevel);
       		int channel=sensorData.getChannel(); // from Properties
       		// Temporary to use old sensor files
       		//    	public static int getChannelFromTiff(String path, String [] suffixes){

       		int channel_file = EyesisCorrectionParameters.CorrectionParameters.getChannelFromTiff(calibFiles[i], new String []{".calib-tiff"}); 
//       		if ((i>=16) && (channel < 16)) {
       		if (channel_file != channel) {
       			System.out.println ("****Temporary fix to use old sensor files, COMMENT OUT when done debugging****");
       			System.out.println ("i="+i+", channel_file="+channel_file+", channel ="+channel+" -> "+channel_file);
       			channel = channel_file; // i;
       		}
       		
       		
       		
       		if ((channel < first_channel) || ((num_channels > 0) && (channel >= (first_channel + num_channels)))) {
       			continue; // wrong channels
       		}
   			channel -= first_channel;
       		if (update_channel) {
       			sensorData.setChannel(channel);
       		}
       		this.sensors[channel]=sensorData;
      		if (channel>maxChannel) maxChannel=channel;
       	}
       	SensorData [] tmpSensors=this.sensors.clone();
       	this.sensors=new SensorData[maxChannel+1];
       	for (int i=0;i<this.sensors.length;i++) this.sensors[i]= tmpSensors[i];
    }
    public int getNumChannels(){
    	return (this.sensors==null)?0:this.sensors.length;
    }

    public String getPath(int channel ){
    	if ((channel<0) || (channel>=this.sensors.length)) return null;
    	return this.sensors[channel].getPath();
    }

    public int getSubChannel(int channel ){
    	if ((channel<0) || (channel>=this.sensors.length)) {
    		System.out.println("ERROR: getSubChannel("+channel+"): channel out of range");
    		return -1;
    	}
//		System.out.println("getSubChannel("+channel+") => "+this.sensors[channel].getSubChannel());
    	return this.sensors[channel].getSubChannel();
    }
    public int getSubChannelSilent(int channel ){
        	if ((channel<0) || (channel>=this.sensors.length)) {
        		return -1;
        	}
//    		System.out.println("getSubChannel("+channel+") => "+this.sensors[channel].getSubChannel());
        	return this.sensors[channel].getSubChannel();
    }
    public int getSubCamera(int channel ){
    	if ((channel<0) || (channel>=this.sensors.length)) return -1;
    	return this.sensors[channel].getSubCamera();
    }

// Enumerated calibration file

    public boolean isChannelAvailable(int channel){
    	return (this.sensors != null) && (channel>=0)  && (channel<this.sensors.length) && (this.sensors[channel]!=null);
    }

    public boolean isChannelAvailable(int [] channels){
    	if (channels == null) return false;
    	for (int c:channels) {
    		if (isChannelAvailable(c)) return true;
    	}
    	return false;
    }

    // Updating for nc393. subCamera here is 0..9 for Eyesis4pi393 - 0-based index of the file, so it combines physical camera (separate IP)
    // as stored in "subcamera" field of the calibration file and "sensor_port". sensor_port may start from non-0, so we need to count all combinations
    //removeUnusedSensorData should be off!


    public boolean subcamerasUsed() {
    	for (int i=0;i<this.sensors.length;i++)  if ((this.sensors[i]!=null) && (this.sensors[i].subcamera >= 0)) {
    		return  true;
    	}
    	return false;
    }

    // if subcameras are not used (-1), as in Quads, return {channel index} or null if alien file

	public int [] getSensorWH(int chn) {
		if ((sensors == null) || (chn < 0) || (chn >= sensors.length)) {
			return null;
		}
		return sensors[chn].getSensorWH();
	}



    public int [] channelsForSubCamera(int subCamera){
    	if (subCamera < 0) {
    		return null;
    	}
    	if (this.sensors == null) return null;

    	if (this.debugLevel>1) {
    		System.out.print("channelsForSubCamera("+subCamera+"),this.sensors.length="+this.sensors.length+": ");
    	}
    	if (!subcamerasUsed()) {
        	for (int i=0;i<this.sensors.length;i++)  if ((this.sensors[i]!=null) && (this.sensors[i].channel == subCamera)) {
        		int [] rslt = {i};
        		return rslt;
        	}
        	return null;
    	}

    	// Eyesis version
    	ArrayList<Point> cam_port = new ArrayList<Point>();
    	for (int i=0;i<this.sensors.length;i++)  if (this.sensors[i]!=null) {
    		Point cp = new Point(this.sensors[i].subcamera, this.sensors[i].sensor_port);
    		if (!cam_port.contains(cp)){
    			cam_port.add(cp);
    		}
    	}
    	Point [] cam_port_arr = cam_port.toArray(new Point[0]);
		Arrays.sort(cam_port_arr, new Comparator<Point>() {
			@Override
			public int compare(Point o1, Point o2) {
				return (o1.x>o2.x)? 1:((o1.x < o2.x)?-1:(o1.y > o2.y)? 1:((o1.y < o2.y)?-1:0));
			}
		});
		// debugging:
		if (subCamera >= cam_port_arr.length) {
			if (this.debugLevel>1) {
				System.out.println("Error: Subcamera "+subCamera+" > that total number of sensor ports in the system = "+cam_port_arr.length);
			}
			return null;
		}
		if (this.debugLevel>1) {
			System.out.println("----- This filename subcamera "+subCamera+": physical camera "+cam_port_arr[subCamera].x+", sensor_port "+cam_port_arr[subCamera].y);
		} else if (this.debugLevel > -2) {
    		System.out.print(".");
		}
    	int numChannels=0;
    	for (int i=0;i<this.sensors.length;i++) if (this.sensors[i]!=null) {
    		if (new Point(this.sensors[i].subcamera, this.sensors[i].sensor_port).equals(cam_port_arr[subCamera])) numChannels++;
    	}
    	int [] result=new int [numChannels];
    	numChannels=0;
    	for (int i=0;i<this.sensors.length;i++) if (this.sensors[i]!=null){
    		if (new Point(this.sensors[i].subcamera, this.sensors[i].sensor_port).equals(cam_port_arr[subCamera])) result[numChannels++]=i;
    	}
    	return result;
    }



    public void removeChannel(int channel){
    	if ((this.sensors != null) && (channel>=0)  && (channel<this.sensors.length)) this.sensors[channel]=null;
    }

	public float [] getBayerFlatFieldFloat(
			int channel,
			int [][] bayer, //{{1,0},{2,1}} GR/BG
			double range) { // max/min
		if ((this.sensors == null) || (channel<0)  && (channel>=this.sensors.length))return null;
		return this.sensors[channel].getBayerFlatFieldFloat(
				bayer,
				range);
	}

	public double [] getBayerFlatField(
			int channel,
			int [][] bayer){ //{{1,0},{2,1}} GR/BG
		if ((this.sensors == null) || (channel<0)  && (channel>=this.sensors.length))return null;
		return this.sensors[channel].getBayerFlatField(
//				width,
//				height,
				bayer);
	}

	public int [][] getDefectsXY(
			int channel){
		if ((this.sensors == null) || (channel<0)  && (channel>=this.sensors.length))return null;
		return this.sensors[channel].getDefectsXY();
	}

	public double[] getDefectsDiff(
			int channel){
		if ((this.sensors == null) || (channel<0)  && (channel>=this.sensors.length))return null;
		return this.sensors[channel].getDefectsDiff();
	}

///SensorData
       	/*
		public float [] getBayerFlatFieldFloat(
				int width,
				int height,
				int [][] bayer){ //{{1,0},{2,1}} GR/BG
       	 *
       	int maxChannel=0;
       	for (int i=0;i<calibFiles.length;i++){
			int indexPeriod=calibFiles[i].indexOf('.',calibFiles[i].lastIndexOf(Prefs.getFileSeparator()));
	    	int channel=Integer.parseInt(calibFiles[i].substring(indexPeriod-2,indexPeriod));
      		if (channel>maxChannel) maxChannel=channel;
       	}
       	this.sensors=new SensorData[maxChannel+1];
       	for (int i=0;i<calibFiles.length;i++){
			int indexPeriod=calibFiles[i].indexOf('.',calibFiles[i].lastIndexOf(Prefs.getFileSeparator()));
	    	int channel=Integer.parseInt(calibFiles[i].substring(indexPeriod-2,indexPeriod));
	    	String channelPath=calibFiles[i].substring(0,indexPeriod-2)+String.format("%02d",channel)+calibFiles[i].substring(indexPeriod);
	    	this.sensors[channel]=new SensorData (channelPath);
       	}
    }
       	*/


    public PixelMapping (int debugLevel){ // boolean just to make a different constructor
    	this.debugLevel=debugLevel;
    }


    public PixelMapping (String defaultPath,boolean ok,int debugLevel){ // boolean just to make a different constructor
    	this.debugLevel=debugLevel;
    	loadChannelMaps (defaultPath);
    }
    public void loadChannelMaps (String defaultPath){
    	String [] extensions={".eqr-tiff", ".eqrect-tiff"};
    	String [] defaultPaths=((defaultPath==null) || defaultPath.equals(""))?null:(new String[1]);
    	if (defaultPaths!=null) defaultPaths[0]=defaultPath;
    	MultipleExtensionsFileFilter parFilter = new MultipleExtensionsFileFilter("",extensions,"equirectangular map "+extensions[0]+" files");
    	String [] eqrectFiles=CalibrationFileManagement.selectFiles(false,
    			"Select equirectangular calibration files",
    			"Select",
    			parFilter,
    			defaultPaths); // String [] defaultPaths);
    	if ((eqrectFiles==null) || (eqrectFiles.length==0)) {
    		IJ.showMessage("No files selected");
    		return;
    	}
    	int maxChannel=0;
    	for (int i=0;i<eqrectFiles.length;i++){
    		int indexPeriod=eqrectFiles[i].indexOf('.',eqrectFiles[i].lastIndexOf(Prefs.getFileSeparator()));
    		int channel=Integer.parseInt(eqrectFiles[i].substring(indexPeriod-2,indexPeriod));
    		if (channel>maxChannel) maxChannel=channel;
    	}
    	if (this.sensors==null)	{
    		this.sensors=new SensorData[maxChannel+1];
    		for (int i=0;i<this.sensors.length;i++) this.sensors[i]=null;
    	}
    	if (this.sensors.length<(maxChannel+1)){
    		SensorData[] tmp=this.sensors.clone();
    		this.sensors=new SensorData[maxChannel+1];
    		for (int i=0;i<this.sensors.length;i++) this.sensors[i]=(i<tmp.length)?tmp[i]:null;
    	}

    	for (int i=0;i<eqrectFiles.length;i++){
    		int indexPeriod=eqrectFiles[i].indexOf('.',eqrectFiles[i].lastIndexOf(Prefs.getFileSeparator()));
    		int channel=Integer.parseInt(eqrectFiles[i].substring(indexPeriod-2,indexPeriod));
    		String channelPath=eqrectFiles[i].substring(0,indexPeriod-2)+String.format("%02d",channel)+eqrectFiles[i].substring(indexPeriod);

    		String msg= "Loading "+channelPath;
    		IJ.showStatus(msg);
    		if (this.debugLevel>0)	System.out.println(msg);
    		if (this.sensors[channel]==null) this.sensors[channel] =new SensorData(channelPath,true);
    		else this.sensors[channel].createEquirectangularMap(channelPath);
    	}
    }

    public void loadChannelEquirectangularMap(
    		int channel,
    		String path){
    	if ((this.sensors==null) || (this.sensors[channel]==null)){
    		String msg="Sensor "+channel+" data is not initialized";
   			IJ.showMessage("Error",msg);
			throw new IllegalArgumentException (msg);
    	}
//    	this.sensors[channel] =new SensorData(path,true);
    	this.sensors[channel].createEquirectangularMap(path);
    }

    public void resampleToEquirectangular(
    		String [] paths,
    		int sourceImageScale,
    		boolean showResults,
    		boolean saveResults,
    		int maxThreads){
    	for (int i=0;i<paths.length;i++){
    		String path=paths[i];
        	int indexPeriod=path.indexOf('.',path.lastIndexOf(Prefs.getFileSeparator()));
    		String outPath=path.substring(0,indexPeriod)+".eqr-tiff";
    		ImagePlus imp_eqr=resampleToEquirectangular(
    	    		path,
    	    		sourceImageScale,
    	    		maxThreads);
    		if (showResults) imp_eqr.show();
    		if (saveResults) {
   				FileSaver fs=new FileSaver(imp_eqr);
   				fs.saveAsTiff(outPath);

    		}
    	}
    }

    public ImagePlus resampleToEquirectangular(
    		String path,
    		int sourceImageScale,
    		int maxThreads){
    	String regexChannel=".*-\\d\\d-.*";
    	if (!path.matches(regexChannel)){
    		String msg="Can not determine channel number - path "+path+" does not contain \"-NN-\"";
    		IJ.showMessage (msg);
    		return null;
    	}
    	int indexPeriod=path.indexOf('.',path.lastIndexOf(Prefs.getFileSeparator()));

		int channel=-1;
    	for (int i=path.lastIndexOf(Prefs.getFileSeparator());i<(indexPeriod-4);i++){
    		if (path.substring(i,i+4).matches(regexChannel)){
    			channel=Integer.parseInt(path.substring(i+1,i+3));
    			break;
    		}
    	}
    	if (channel<0){
    		String msg="Can not determine channel number - path "+path+" does not contain \"-NN-\"";
    		IJ.showMessage (msg);
    		return null;
    	}
    	return  resampleToEquirectangular(
        		path,
        		channel,
        		sourceImageScale,
        		maxThreads);
    }
    public ImagePlus resampleToEquirectangular(
    		String path,
    		int channel,
    		int sourceImageScale,
    		int maxThreads){
    	Opener opener=new Opener();
    	ImagePlus imp=opener.openImage("", path);
    	if (imp==null) {
    		String msg="Failed to  "+path;
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
    	}
    	(new JP46_Reader_camera(false)).decodeProperiesFromInfo(imp);
    	return resampleToEquirectangular(
        		imp,
        		channel,
        		sourceImageScale,
        		maxThreads);
    }
    public void  deleteEquirectangularMapAll(){
    	if (this.sensors!=null){
    	for (int channel=0;channel<this.sensors.length;channel++) deleteEquirectangularMapAll(channel);
    	}
    }
    public void  deleteEquirectangularMapFull(){
    	if (this.sensors!=null){
    	for (int channel=0;channel<this.sensors.length;channel++) deleteEquirectangularMapFull(channel);
    	}
    }

    public void  deleteEquirectangularMapAll(int channel){
    	if ((this.sensors!=null) && (this.sensors.length>channel) && (this.sensors[channel]!=null)){
    		this.sensors[channel].equirectangularMap=null;
    	}
    }
    public void  deleteEquirectangularMapFull(int channel){
    	if ((this.sensors!=null) && (this.sensors.length>channel) && (this.sensors[channel]!=null) && (this.sensors[channel].equirectangularMap!=null)){
    		this.sensors[channel].equirectangularMap.map=null;
    	}
    }

    public boolean  isEquirectangularMapAvailable(int channel){
    	if ((this.sensors!=null) && (this.sensors.length>channel) && (this.sensors[channel]!=null)){
    		return((this.sensors[channel].equirectangularMap!=null) && (this.sensors[channel].equirectangularMap.partialMap!=null));
    	}
    	return false;
    }
    public boolean  isPlaneMapMapAvailable(int channel){
    	if ((this.sensors!=null) && (this.sensors.length>channel) && (this.sensors[channel]!=null)){
    		return((this.sensors[channel].interSensor!=null));
    	}
    	return false;
    }

    public ImagePlus resampleToEquirectangular(
    		ImagePlus imp,
    		int channel,
    		int sourceImageScale,
    		int maxThreads){
    	if (imp==null) return null;
		if ((this.sensors==null) || (this.sensors.length<(channel+1) || (this.sensors[channel]==null))){
			IJ.showMessage("No sensor data for channel "+channel);
			return null;
		}
		if (( this.sensors[channel].equirectangularMap==null) || ( this.sensors[channel].equirectangularMap.partialMap==null)){
			IJ.showMessage("No equirectangular map for channel "+channel);
			return null;
		}
		if (this.lanczos==null)   generateLanczosStack();

		if (imp.getType()==ImagePlus.COLOR_RGB) {
			return resampleToEquirectangularRGB24(imp, channel, sourceImageScale, maxThreads);
		} else if (imp.getStackSize()>=3) {
			return resampleToEquirectangularRGBFP32(imp, channel, sourceImageScale, maxThreads);
		}
		IJ.showMessage("Not yet implemented for this image type");
    	return null;
    }
/**
 * Warp 3-slice 32-bit FP stack (RGB) to 4-slice stack RGBA
 * @param imp 3-slice image stack with straight image
 * @param channel sensor channel number
 * @param scale source image scale (normally 2x)
 * @param maxThreads
 * @return 4-slice RGBA image, floating point 32 bits per sample per color channel
 */

    public ImagePlus resampleToEquirectangularRGBFP32(
    		ImagePlus imp,
    		int channel,
    		final int scale,
    		int maxThreads){
    	final SensorData.EquirectangularMap erm=this.sensors[channel].equirectangularMap;
    	final int width=imp.getWidth();
    	final int height=imp.getHeight();
		if (scale==0) {
			IJ.showMessage("resampleToEquirectangularRGBFP32(): Image has less resolution than the map, can not proceed");
			return null;
		}
		if (imp.getStackSize()<3) {
			IJ.showMessage("resampleToEquirectangularRGBFP32(): Expecting 3-layer stack(R,G,B)");
			return null;
		}
		final float [][] imagePixels= new float [3][];
		for (int c=0;c<imagePixels.length;c++) imagePixels[c]=	(float []) imp.getStack().getPixels(c+1);
		final float [][] outPixels=new float [4][erm.mapWOI.width*erm.mapWOI.height];
		System.out.println("resampleToEquirectangularRGBFP32(): imagePixels[0].length.length="+imagePixels[0].length);
		if (this.oversampled!=scale) generateLanczosStack(
        		this.lanczosA,
        		scale,
        		this.binsPerHalfPixel);
		final double [][][][] lanczos=this.lanczos;
		final int center=this.lanczos[0][0][0].length/2;
   		final Thread[] threads = newThreadArray(maxThreads);
   		final AtomicInteger opyAtomic     = new AtomicInteger(0);
   		final AtomicInteger opyDoneAtomic = new AtomicInteger(1);
   		final int binsPerHalfPixel=this.binsPerHalfPixel;
//   		for (int opy=0;opy<erm.mapWOI.height;opy++){
   		IJ.showStatus("Warping image channel "+channel);
   		for (int ithread = 0; ithread < threads.length; ithread++) {
   			threads[ithread] = new Thread() {
   				@Override
				public void run() {
   					for (int opy=opyAtomic.getAndIncrement(); opy<erm.mapWOI.height;opy=opyAtomic.getAndIncrement()){

   						for (int opx=0;opx<erm.mapWOI.width;opx++){

   							int oIndex=opy*erm.mapWOI.width+opx;
   							double [] RGBA={erm.partialMap[2][oIndex],0.0,0.0,0.0};
   							outPixels[3][oIndex]=(float) RGBA[0];
   							if (outPixels[3][oIndex]>0){ // do not convolve pixels with alpha=0;
   								double x=scale*erm.partialMap[0][oIndex];
   								double y=scale*erm.partialMap[1][oIndex];
   								int ix= (int) Math.round(x);
   								int iy= (int) Math.round(y);
   								double dx=x-ix;
   								double dy=y-iy;
   								int indxX= (int) Math.round((2*dx+1)*binsPerHalfPixel);
   								int indxY= (int) Math.round((2*dy+1)*binsPerHalfPixel);
   								double [][] lk=lanczos[indxY][indxX];
   								for (int i=0;i<lk.length;i++) {
   									int ipy=iy+i-center;
   									if (ipy<0) ipy=0;
   									else if (ipy>=height) ipy=height-1;
   									int baseI=ipy*width;
   									for (int j=0;j<lk[0].length;j++){
   										int ipx=ix+j-center;
   										if (ipx<0) ipx=0;
   										else if (ipx>=width) ipx=width-1;
   										/*
   										int pix=imagePixels[baseI+ipx];
   										RGBA[1]+=((pix>>16) & 0xff)*lk[i][j]; // R
   										RGBA[2]+=((pix>>8)  & 0xff)*lk[i][j]; // G
   										RGBA[3]+=( pix      & 0xff)*lk[i][j]; // B
   										*/
   										if ((baseI+ipx)>imagePixels[0].length){
   											System.out.println("baseI+ipx="+(baseI+ipx)+" baseI="+baseI+" ipx="+ipx+" i="+i+" j="+j+
   													" erm.mapWOI.width="+erm.mapWOI.width+" erm.mapWOI.height="+erm.mapWOI.height);
   										}
   										RGBA[1]+=imagePixels[0][baseI+ipx]*lk[i][j]; // R
   										RGBA[2]+=imagePixels[1][baseI+ipx]*lk[i][j]; // G
   										RGBA[3]+=imagePixels[2][baseI+ipx]*lk[i][j]; // B

   									}

   								}
   								outPixels[0][oIndex]=(float) RGBA[1];
   								outPixels[1][oIndex]=(float) RGBA[2];
   								outPixels[2][oIndex]=(float) RGBA[3];
   								/*
   								int c=0;
   								for (int i=0;i<4;i++){
   									if      (RGBA[i]< 0.0)    RGBA[i]=0.0;
   									else if (RGBA[i] > 255.0) RGBA[i]=255.0;
   									c=(c<<8) | ((int) RGBA[i]);
   									outPixels [oIndex]=c;
   								}
   								*/
   							} else {
   								outPixels[0][oIndex]=0.0F;
   								outPixels[1][oIndex]=0.0F;
   								outPixels[2][oIndex]=0.0F;
//   								outPixels[oIndex]=0;
   							}

   						}
    					final int numFinished=opyDoneAtomic.getAndIncrement();
    					//					IJ.showProgress(progressValues[numFinished]);
    					SwingUtilities.invokeLater(new Runnable() {
    						@Override
							public void run() {
    							// Here, we can safely update the GUI
    							// because we'll be called from the
    							// event dispatch thread
    							IJ.showProgress(numFinished,erm.mapWOI.height);
    						}
    					});
   					}
   				}
   			};
   		}
   		startAndJoin(threads);
   		ImageStack outStack=new ImageStack(erm.mapWOI.width,erm.mapWOI.height);
   		outStack.addSlice("Red",   outPixels[0]);
   		outStack.addSlice("Green", outPixels[1]);
   		outStack.addSlice("Blue",  outPixels[2]);
   		outStack.addSlice("Alpha", outPixels[3]);
//		ColorProcessor cp=new ColorProcessor(erm.mapWOI.width,erm.mapWOI.height);
//		cp.setPixels(outPixels);
		ImagePlus impOut=new ImagePlus(imp.getTitle()+"_EQR",outStack);
		impOut.setProperty("channel",   ""+erm.channel);  // image channel
		impOut.setProperty("XPosition", ""+erm.mapWOI.x); // real XPosition as tiff tag 0x11e is rational (5)
		impOut.setProperty("YPosition", ""+erm.mapWOI.y); // real XPosition as tiff tag 0x11f is rational (5)
		impOut.setProperty("ImageFullWidth", "" +erm.pixelsHorizontal);// real ImageFullWidth as tiff tag 0x8214 is rational (5)
		impOut.setProperty("ImageFullLength", ""+erm.pixelsVertical);   // real ImageFullLength as tiff tag 0x8215 is rational (5)
		(new JP46_Reader_camera(false)).encodeProperiesToInfo(impOut);
		impOut.getProcessor().resetMinAndMax();
		return impOut;
    }
    public ImagePlus resampleToEquirectangularRGB24(
    		ImagePlus imp,
    		int channel,
    		final int scale,
    		int maxThreads){
    	final SensorData.EquirectangularMap erm=this.sensors[channel].equirectangularMap;
    	final int width=imp.getWidth();
    	final int height=imp.getHeight();
		if (scale==0) {
			IJ.showMessage("Image has less resolution than the map, can not proceed");
			return null;
		}
		final int [] imagePixels=(int []) imp.getProcessor().getPixels();
		final int [] outPixels=new int [erm.mapWOI.width*erm.mapWOI.height];
		if (this.oversampled!=scale) generateLanczosStack(
        		this.lanczosA,
        		scale,
        		this.binsPerHalfPixel);
		final double [][][][] lanczos=this.lanczos;
		final int center=this.lanczos[0][0][0].length/2;
   		final Thread[] threads = newThreadArray(maxThreads);
   		final AtomicInteger opyAtomic     = new AtomicInteger(0);
   		final AtomicInteger opyDoneAtomic = new AtomicInteger(1);
   		final int binsPerHalfPixel=this.binsPerHalfPixel;
//   		for (int opy=0;opy<erm.mapWOI.height;opy++){
   		IJ.showStatus("Warping image channel "+channel);
   		for (int ithread = 0; ithread < threads.length; ithread++) {
   			threads[ithread] = new Thread() {
   				@Override
				public void run() {
   					for (int opy=opyAtomic.getAndIncrement(); opy<erm.mapWOI.height;opy=opyAtomic.getAndIncrement()){

   						for (int opx=0;opx<erm.mapWOI.width;opx++){

   							int oIndex=opy*erm.mapWOI.width+opx;
   							double [] RGBA={255*erm.partialMap[2][oIndex],0.0,0.0,0.0};
   							if (RGBA[0]>0){ // do not convolve pixels with alpha=0;
   								double x=scale*erm.partialMap[0][oIndex];
   								double y=scale*erm.partialMap[1][oIndex];
   								int ix= (int) Math.round(x);
   								int iy= (int) Math.round(y);
   								double dx=x-ix;
   								double dy=y-iy;
   								int indxX= (int) Math.round((2*dx+1)*binsPerHalfPixel);
   								int indxY= (int) Math.round((2*dy+1)*binsPerHalfPixel);
   								double [][] lk=lanczos[indxY][indxX];
   								for (int i=0;i<lk.length;i++) {
   									int ipy=iy+i-center;
   									if (ipy<0) ipy=0;
   									else if (ipy>=height) ipy=height-1;
   									int baseI=ipy*width;
   									for (int j=0;j<lk[0].length;j++){
   										int ipx=ix+j-center;
   										if (ipx<0) ipx=0;
   										else if (ipx>=width) ipx=width-1;
   										int pix=imagePixels[baseI+ipx];
   										RGBA[1]+=((pix>>16) & 0xff)*lk[i][j]; // R
   										RGBA[2]+=((pix>>8)  & 0xff)*lk[i][j]; // G
   										RGBA[3]+=( pix      & 0xff)*lk[i][j]; // B
   									}

   								}
   								int c=0;
   								for (int i=0;i<4;i++){
   									if      (RGBA[i]< 0.0)    RGBA[i]=0.0;
   									else if (RGBA[i] > 255.0) RGBA[i]=255.0;
   									c=(c<<8) | ((int) RGBA[i]);
   									outPixels [oIndex]=c;
   								}
   							} else {
   								outPixels[oIndex]=0;
   							}

   						}
    					final int numFinished=opyDoneAtomic.getAndIncrement();
    					//					IJ.showProgress(progressValues[numFinished]);
    					SwingUtilities.invokeLater(new Runnable() {
    						@Override
							public void run() {
    							// Here, we can safely update the GUI
    							// because we'll be called from the
    							// event dispatch thread
    							IJ.showProgress(numFinished,erm.mapWOI.height);
    						}
    					});
   					}
   				}
   			};
   		}
   		startAndJoin(threads);
		ColorProcessor cp=new ColorProcessor(erm.mapWOI.width,erm.mapWOI.height);
		cp.setPixels(outPixels);
		ImagePlus impOut=new ImagePlus(imp.getTitle()+"_EQR",cp);
		impOut.setProperty("channel",   ""+erm.channel);  // image channel
		impOut.setProperty("XPosition", ""+erm.mapWOI.x); // real XPosition as tiff tag 0x11e is rational (5)
		impOut.setProperty("YPosition", ""+erm.mapWOI.y); // real XPosition as tiff tag 0x11f is rational (5)
		impOut.setProperty("ImageFullWidth", "" +erm.pixelsHorizontal);// real ImageFullWidth as tiff tag 0x8214 is rational (5)
		impOut.setProperty("ImageFullLength", ""+erm.pixelsVertical);   // real ImageFullLength as tiff tag 0x8215 is rational (5)
		(new JP46_Reader_camera(false)).encodeProperiesToInfo(impOut);
		impOut.getProcessor().resetMinAndMax();
		return impOut;
    }




    /**
     * Generate a stack of Lanczos kernels to for re-sampling
     * @param lacrososA - value of "a" in sinc(x)*sinc(x/a)
     * @param oversampled reduce resolution of the image if it was over-sampled (2 for current Eeyesis software)
     * @param binsPerHalfPixel calculate for this number of fractional pixels (oversampled pixels, not sensor ones)
     * @return [fractional pixels vertical] [fractional pixels horizontal][delta pixel horizontal][ delta pixel vertical]
     */
    public double  [][][][] generateLanczosStack(){
    	return generateLanczosStack(
        		this.lanczosA,
        		this.oversampled,
        		this.binsPerHalfPixel);
    }


    public double  [][][][] generateLanczosStack(
    		int lanczosA,
    		int oversampled,
    		int binsPerHalfPixel){
    	this.lanczosA=lanczosA;
    	this.oversampled=oversampled;
    	this.binsPerHalfPixel=binsPerHalfPixel;
    	if (this.debugLevel>0) System.out.println("Generating Lanczos kernel stack with A="+ this.lanczosA+", oversampled="+oversampled+", binsPerHalfPixel="+binsPerHalfPixel);
    	double [][][][] lanczos= new double [2*binsPerHalfPixel+1][2*binsPerHalfPixel+1][2*lanczosA*oversampled+1][2*lanczosA*oversampled+1];
    	double step=0.5/binsPerHalfPixel;
    	double centerPix=lanczosA*oversampled;
    	double pi2=Math.PI*Math.PI;
    	double s=lanczosA/pi2/oversampled;
    	double sync0=1.0/oversampled;
    	for (int iy=0;iy<lanczos.length;iy++){
    		double yc=-0.5+iy*step;
        	for (int ix=0;ix<lanczos[0].length;ix++){
        		double xc=-0.5+ix*step;
        		for (int ipy=0;ipy<lanczos[0][0].length;ipy++){
        			double py=(ipy-centerPix-yc)/oversampled;
        			double lanczosY=(Math.abs(py)>=lanczosA)?0.0:((py==0.0)?sync0:(s*Math.sin(py*Math.PI)*Math.sin(py*Math.PI/lanczosA)/(py*py)));
            		for (int ipx=0;ipx<lanczos[0][0][0].length;ipx++){
            			if (lanczosY==0.0) lanczos[iy][ix][ipy][ipx]=0.0;
            			else {
            				double px=(ipx-centerPix-xc)/oversampled;
            				lanczos[iy][ix][ipy][ipx]=lanczosY*
            				((Math.abs(px)>=lanczosA)?0.0:((px==0.0)?sync0:(s*Math.sin(px*Math.PI)*Math.sin(px*Math.PI/lanczosA)/(px*px))));
            			}
            		}
        		}
        	}
    	}
    	this.lanczos=lanczos;
    	if (this.debugLevel>0) System.out.println("Generated Lanczos kernel stack with A="+ this.lanczosA+", oversampled="+oversampled+", binsPerHalfPixel="+binsPerHalfPixel);
    	return lanczos;
    }


    public double [] testLanczosCenter(double [][][][] lanczos){
    	int size=    lanczos.length;
    	int center=	lanczos[0][0].length/2;
    	double [] testArr=new double [size*size];
    	for (int i=0;i<testArr.length;i++) testArr[i]=0.0;
    	for (int iy=0;iy<lanczos.length;iy++){
        	for (int ix=0;ix<lanczos[0].length;ix++){
            			testArr[iy*size+ix]=lanczos[iy][ix][center][center];
        	}
    	}
    	return testArr;
    }

    public double [] testLanczosStack(double [][][][] lanczos){
    	int halfDirSamples= lanczos.length/2;
    	int halfKernelSize= lanczos[0][0].length/2;
    	int halfSize=halfDirSamples + halfKernelSize;
//    	int size=    2*halfSize+1;
    	int size=    lanczos.length*lanczos[0][0].length+1;
    	System.out.println(
    			"testLanczosStack() lanczos.size="+lanczos.length+"\n"+
    			"  lanczos[0].length="+lanczos[0].length+"\n"+
    			"  lanczos[0][0].length="+lanczos[0][0].length+"\n"+
    			"  lanczos[0][0][0].length="+lanczos[0][0][0].length+"\n"+
    			"  halfDirSamples="+halfDirSamples+"\n"+
    			"  halfKernelSize="+halfKernelSize+"\n"+
    			"  halfSize="+halfSize+"\n"+
    			"  size="+size+"\n");

    	double [] testArr=new double [size*size];
    	for (int i=0;i<testArr.length;i++) testArr[i]=0.0;
    	for (int iy=0;iy<(lanczos.length-1);iy++){
        	for (int ix=0;ix<(lanczos[0].length-1);ix++){
        		for (int ipy=0;ipy<lanczos[0][0].length;ipy++){
        			int fullIY=iy+(2*halfDirSamples)*(lanczos[0][0].length-1-ipy);
            		for (int ipx=0;ipx<lanczos[0][0][0].length;ipx++){
            			int fullIX=ix+(2*halfDirSamples)*(lanczos[0][0][0].length-1-ipx);
            			testArr[fullIY*size+fullIX]=lanczos[iy][ix][ipy][ipx];
            		}
        		}
        	}
    	}
    	return testArr;
    }


    public void normalizeLanczosStack(double [][][][] lanczos){
    	for (int iy=0;iy< lanczos.length;iy++){
        	for (int ix=0;ix<lanczos[0].length;ix++){
        		double sum=0;
        		for (int ipy=0;ipy<lanczos[0][0].length;ipy++){
            		for (int ipx=0;ipx<lanczos[0][0][0].length;ipx++){
            			sum+=lanczos[iy][ix][ipy][ipx];
            		}
        		}
        		double k=1.0/sum;
        		for (int ipy=0;ipy<lanczos[0][0].length;ipy++){
            		for (int ipx=0;ipx<lanczos[0][0][0].length;ipx++){
            			lanczos[iy][ix][ipy][ipx]*=k;
            		}
        		}
        	}
    	}
    }





    public double [] testLanczosStackNormalization(double [][][][] lanczos){
    	int size=    lanczos.length;
    	double [] testArr=new double [size*size];
    	for (int i=0;i<testArr.length;i++) testArr[i]=0.0;
    	for (int iy=0;iy< lanczos.length;iy++){
        	for (int ix=0;ix<lanczos[0].length;ix++){
        		for (int ipy=0;ipy<lanczos[0][0].length;ipy++){
            		for (int ipx=0;ipx<lanczos[0][0][0].length;ipx++){
            			testArr[iy*size+ix]+=lanczos[iy][ix][ipy][ipx];
            		}
        		}
        	}
    	}
    	return testArr;
    }


    public float[] generateOverlapMap(){
    	float [] pixels=null;
    	for (int channelNumber=0;channelNumber<this.sensors.length;channelNumber++) if ((this.sensors[channelNumber]!=null) && (this.sensors[channelNumber].equirectangularMap!=null)){
    		if (this.debugLevel>0) System.out.println("Accumulating channel "+channelNumber);
    		SensorData.EquirectangularMap erm=this.sensors[channelNumber].equirectangularMap;
    		if (pixels==null){
    			this.panoWidth=erm.pixelsHorizontal;
    			this.panoHeight=erm.pixelsVertical;
    			pixels= new float [this.panoWidth*this.panoHeight];
    			for (int i=0;i<pixels.length;i++) pixels[i]=0.0f;
    		}
    		int tileIndex=0;
    		Rectangle woi=erm.mapWOI;
    		for (int iTileLat=0;iTileLat<erm.mapWOI.height;iTileLat++){
    			for (int iTileLong=0;iTileLong<erm.mapWOI.width;iTileLong++){
    				if (erm.partialMap[2][tileIndex++]>0){ // alpha
    					pixels[(iTileLat+woi.y)*this.panoWidth+((iTileLong+woi.x) % this.panoWidth)]+=erm.partialMap[2][tileIndex-1];
    				}
    			}
    		}
    	}
    	return pixels;
    }
    /**
     *  For up to 32 sensors, calculates in which of them current lat/long is available
     * @param threshold - count only pixels with alpha > this threshold (0.0 - all non-zero)
     * @return bit masks in scan-line (long, then lat) order
     * sets this.panoDegreesPerPixel
     */

    public int [] generateOverlapBits( double threshold, int [] sensorList ){
    	if (sensorList==null){ // all sensors
    		sensorList=new int [this.sensors.length];
    		for (int i=0;i<sensorList.length;i++) sensorList[i]=i;
    	}
    	int [] usedSensors=null;
//    	for (int channelNumber=0;channelNumber<this.sensors.length;channelNumber++) if ((this.sensors[channelNumber]!=null) && (this.sensors[channelNumber].equirectangularMap!=null)){
    	for (int channelIndex=0;channelIndex<sensorList.length;channelIndex++) {
    		int channelNumber=sensorList[channelIndex];
    		if ((this.sensors[channelNumber]!=null) && (this.sensors[channelNumber].equirectangularMap!=null)){


    			if (this.debugLevel>0) System.out.println("Accumulating channel "+channelNumber);
    			SensorData.EquirectangularMap erm=this.sensors[channelNumber].equirectangularMap;
    			if (usedSensors==null){
    				this.panoWidth=erm.pixelsHorizontal;
    				this.panoHeight=erm.pixelsVertical;
    				usedSensors= new int [this.panoWidth*this.panoHeight];
    				for (int i=0;i<usedSensors.length;i++) usedSensors[i]=0;
    				this.panoLongitudeLeft=erm.longitudeLeft;
    				this.panoLongitudeRight=erm.longitudeRight;
    				this.panoLatitudeTop=erm.latitudeTop;
    				this.panoLatitudeBottom=erm.latitudeBottom;
    				this.panoDegreesPerPixel=(this.panoLongitudeRight-this.panoLongitudeLeft)/(this.panoWidth-0);
    			}
    			int tileIndex=0;
    			Rectangle woi=erm.mapWOI;
    			for (int iTileLat=0;iTileLat<erm.mapWOI.height;iTileLat++){
    				for (int iTileLong=0;iTileLong<erm.mapWOI.width;iTileLong++){
    					if (erm.partialMap[2][tileIndex++]>threshold){ // alpha
    						usedSensors[(iTileLat+woi.y)*this.panoWidth+((iTileLong+woi.x) % this.panoWidth)] |= (1<<channelNumber);
    					}
    				}
    			}
    		}
    	}
    	return usedSensors;
    }
    /**
     * Calculate overlap areas for pairs of sensors, measured in angular pixels (equirectangular at equator)
     * @param threshold
     * @return square symmetrical 2-d array [numSesnsors][numSesnsors]
     */
    public double [][] overlapPairsAreas(double threshold){
    	int [] usedSensors=generateOverlapBits(threshold,null); // also sets this.panoWidth, this.panoHeight
    	int numChannels=this.sensors.length;
    	double [][] overlap=new double [numChannels][numChannels];
    	for (int i=0;i<numChannels;i++) for (int j=0;j<numChannels;j++) overlap[i][j]=0.0;
    	for (int iLat=0; iLat<this.panoHeight;iLat++){
    		double lat=this.panoLatitudeTop-this.panoDegreesPerPixel*iLat; // top: +90, bottom: -90
    		double pixelWeight=Math.cos(Math.PI*lat/180.0);
    		for (int iLong=0; iLong<this.panoWidth;iLong++){
        		int us=usedSensors[iLat*this.panoWidth+iLong];
    			for (int i=0;i<numChannels;i++) if ((us & (1<<i))!=0){
        			for (int j=i;j<numChannels;j++) if ((us & (1<<j))!=0){
        				overlap[i][j]+=pixelWeight; // will include [i][i]
        			}
    			}
    		}
    	}
    	// make it symmetrical
    	for (int i=0;i<(numChannels-1);i++) for (int j=i+1;j<numChannels;j++) overlap[j][i]=overlap[i][j];
    	return overlap;
    }

    /**
     * Generate a list of sensor pairs that have sufficient overlap areas
     * @param alphaThreshold minimal alpha value to consider (i.e. 0.0 - all with non-zero alpha)
     * @param overlapThreshold minimal fraction of the overlap area to square root of the product of the sensor areas
     * @return list of sensor number pairs
     */

    public int [][] findSensorPairs(
    		double alphaThreshold,
    		double overlapThreshold){
    	int numPairs=0;
    	double [][] overlap=overlapPairsAreas(alphaThreshold);
		for (int i=0;i<(overlap.length-1);i++) for (int j=i+1;j<overlap.length;j++){
			overlap[i][j]/=Math.sqrt(overlap[i][i]*overlap[j][j])/100.0;
			if (overlap[i][j]>=overlapThreshold) numPairs++;
		}
		int [][] pairs = new int [numPairs][2];
		numPairs=0;
		for (int i=0;i<(overlap.length-1);i++) for (int j=i+1;j<overlap.length;j++){
			if (overlap[i][j]>=overlapThreshold) {
				pairs[numPairs  ][0]=i;
				pairs[numPairs++][1]=j;
			}
		}
		return pairs;
    }
  /*
			double threshold=0.01*gd.getNextNumber();
			double [][] overlap=PIXEL_MAPPING.overlapPairsAreas(threshold);
			PIXEL_MAPPING.listOverlap(overlap);
			for (int i=0;i<(overlap.length-1);i++) for (int j=i+1;j<overlap.length;j++){
				overlap[i][j]/=Math.sqrt(overlap[i][i]*overlap[j][j])/100.0;
			}
			PIXEL_MAPPING.listOverlap(overlap);

   */

    public void listOverlap(double [][] overlap){
    	StringBuffer sb=new StringBuffer();
    	String header="Sensor ";
    	int numChannels=this.sensors.length;
    	for (int i=0;i<numChannels;i++) header +="\t"+i;
    	for (int i=0;i<numChannels;i++){
    		sb.append (i+"");
    		for (int j=0;j<numChannels;j++){
    			sb.append ("\t"+IJ.d2s(overlap[i][j],1));
    		}
    		sb.append ("\n");
    	}
	    new TextWindow("sensor overlaps", header, sb.toString(), 900,900);
    }

    public int getPanoWidth() {return this.panoWidth;}
    public int getPanoHeight() {return this.panoHeight;}

    public void listParameters(){
        int numSubCameras=getNumSubCameras();
        if (numSubCameras==0) return;

        double [][] cameraPars=new double [numSubCameras][];
        for (int i=0;i<numSubCameras;i++) cameraPars[i]=getParametersVector(i);
        // parameters same order as in this
 	    String header="Name\tUnits";
 	    for (int i=0;i<numSubCameras;i++) header+="\t"+i;
 	    StringBuffer sb = new StringBuffer();
 	    for (int n=0;n<cameraPars[0].length;n++) {
 	    	sb.append(getParameterName(n)+"\t"+getParameterUnits(n));
 		    for (int i=0;i<numSubCameras;i++) sb.append("\t"+IJ.d2s(cameraPars[i][n],3));
 		    sb.append("\n");
 	    }
 	    new TextWindow("Camera parameters", header, sb.toString(), 85*(numSubCameras+2),350);
    }

    public void generateAndSaveEquirectangularMaps(
    		String path,
    		double longitudeLeft,
    		double longitudeRight,
    		double latitudeTop,
    		double latitudeBottom,
    		int    pixelsHorizontal,
    		int    imageWidth,
    		int    imageHeight,
    		double x0,
    		double y0,
    		double pixelStep,
    		int    longitudeCrop,
    		boolean deleteFull,
    		boolean deleteAll,
    		int maxThreads){

    	if (this.sensors==null) {
    		return;
    	}
		long 	  startTime=System.nanoTime();
    	int indexPeriod=path.indexOf('.',path.lastIndexOf(Prefs.getFileSeparator()));
    	int indexSuffix=indexPeriod;
    	String digits="0123456789";
    	for (int i=0;i<2;i++) if (digits.indexOf(path.charAt(indexSuffix-1))>=0) indexSuffix--; // remove 1 or 2 digits before period

    	for (int channelNumber=0;channelNumber<this.sensors.length;channelNumber++) if (this.sensors[channelNumber]!=null){
    		String msg="Creating equirectangular map, channel "+(channelNumber+1)+" of "+this.sensors.length;
    		IJ.showStatus(msg);
    		if (this.debugLevel>0) System.out.println(msg);
    		sensors[channelNumber].combineDistortionsEquirectangularToSensor(
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
    				pixelStep, //double pixelStep,
    				maxThreads);
    		sensors[channelNumber].equirectangularMap.createPartialMap(longitudeCrop);
    		if (deleteFull) sensors[channelNumber].equirectangularMap.map=null;
    		String channelPath=path.substring(0,indexSuffix)+String.format("%02d",channelNumber)+path.substring(indexPeriod);
    		msg="Saving equirectangular map to "+channelPath;
    		IJ.showStatus(msg);
    		if (this.debugLevel>0) System.out.println(msg);
    		sensors[channelNumber].equirectangularMap.saveMapAsImageStack("EQRCT_MAP"+channelNumber, channelPath);
    		if (deleteAll) sensors[channelNumber].equirectangularMap=null;
			if (this.debugLevel >1) System.out.println("Saved equirectangular map "+channelNumber+" at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));

    	}
		if (this.debugLevel >0) System.out.println("Finished equirectangular maps at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		IJ.showStatus("Finished equirectangular maps");
    }

    public int getNumSubCameras(){
    	return (this.sensors==null)?0:this.sensors.length;
    }
    public double [] getParametersVector(int chn){
    	return this.sensors[chn].getParametersVector();
    }

    public String getParameterName(int index){
    	return (new SensorData()).getParameterName(index);
    }
    public String getParameterDescription(int index){
    	return (new SensorData()).getParameterDescription(index);
    }
    public String getParameterUnits(int index){
    	return (new SensorData()).getParameterUnits(index);
    }

    public double[][] applyOverlapMap(
    		String path, //
    		String imgPathFormat,
    		String resultDirectory,
    		int sourceImageScale, // 2.0
    		boolean saveTiff,
    		boolean convertToDouble,
    		int maxThreads,
    		int debugLevel
    		){
		Opener opener=new Opener();
		ImagePlus impMap=opener.openImage("", path);
    	if (impMap==null) {
    		String msg="Failed to read inter-sensor overlap map file "+path;
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
    	}
		if (debugLevel>0) System.out.println("Read "+path+" as an inter-sensor overlap map");
    	(new JP46_Reader_camera(false)).decodeProperiesFromInfo(impMap);
    	InterSensor interSensor=new InterSensor(impMap,debugLevel);
    	this.lastUsedInterSensor=interSensor;
    	double [][] aYCbCr=convertToDouble?(new double [8][]):null;
    	for (int iChn=0;iChn<2;iChn++){
    		int channel=interSensor.channel[iChn];
    		String sourcePath=String.format(imgPathFormat,channel);
    		ImagePlus impSrc=opener.openImage("", sourcePath);
        	if (impSrc==null) {
        		String msg="Failed to read sensor image file "+sourcePath;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
    		if (debugLevel>0) System.out.println("Read "+sourcePath+" as a sensor image");
    		if (this.lanczos==null)   generateLanczosStack();
    		if (this.oversampled!=sourceImageScale) generateLanczosStack(
            		this.lanczosA,
            		sourceImageScale,
            		this.binsPerHalfPixel);
    		if (impSrc.getType()!=ImagePlus.COLOR_RGB) {
        		String msg="Not yet implemented for this image type, only RGB24 is supported";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
    		}
    		ImagePlus impOverlap=interSensor.resampleToOverlapRGB24(
    				true, // it is overlap of a pair of images
    				impSrc,
    				channel,
    				sourceImageScale,
        			this.lanczos,
        			this.binsPerHalfPixel,
        			maxThreads);
    		if (impOverlap==null){
        		String msg="Failed to apply overlap map to "+sourcePath;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
    		}
    		if (debugLevel>1){
    			impOverlap.getProcessor().resetMinAndMax(); // imp_psf will be reused
    			impOverlap.show();
    		}
    		//resultDirectory
    		if (saveTiff){
    			String outPath=((resultDirectory.length()<1)?"":(resultDirectory+Prefs.getFileSeparator()))+impOverlap.getTitle()+".tiff";
    			try {
    				(new EyesisTiff()).saveTiff(
    						impOverlap,
    						outPath,
    						3, // float32
    						255.0,
    						true,
    						debugLevel);
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
    		}
    		if (convertToDouble) interSensor.argb24ToAYCbCr(
    				impOverlap,
    				iChn, // 0 - first image, 1 - second image in the pair
        			aYCbCr);
    	}
    	if (convertToDouble) interSensor.overlapImages=aYCbCr;
    	return aYCbCr;
    }

    public InterSensor loadPlaneMap(
    		String path, //
    		int debugLevel
    		){
		Opener opener=new Opener();
		ImagePlus impMap=opener.openImage("", path);
    	if (impMap==null) {
    		String msg="Failed to read plane-to-sensor map file "+path;
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
    	}
		if (debugLevel>0) System.out.println("Read "+path+" as a plane-to-sensor map");
    	(new JP46_Reader_camera(false)).decodeProperiesFromInfo(impMap);
    	InterSensor interSensor=new InterSensor(impMap,debugLevel);
    	this.lastUsedInterSensor=interSensor;
    	for (int iChn=0;iChn<interSensor.channel.length;iChn++){
    		int channel=interSensor.channel[iChn];
    		if (this.sensors[channel]==null){
    			if (debugLevel>1) System.out.println("loadPlaneMap():this.sensors["+channel+"]==null");
    		} else {
    		  this.sensors[channel].interSensor=interSensor; // OK if some sensors are not used for the selected source files
    		}
    	}
    	return interSensor;
    }


    public void applyPlaneMap(
    		String path, //
    		String imgPathFormat,
    		String resultDirectory,
    		int sourceImageScale, // 2.0
    		boolean saveTiff,
    		int maxThreads,
    		int debugLevel
    		){
    	InterSensor interSensor=loadPlaneMap(
        		path, //
        		debugLevel
        		);
		Opener opener=new Opener();
    	for (int iChn=0;iChn<interSensor.channel.length;iChn++){
    		int channel=interSensor.channel[iChn];
    		String sourcePath=String.format(imgPathFormat,channel);
    		ImagePlus impSrc=opener.openImage("", sourcePath);
        	if (impSrc==null) {
        		if (debugLevel>0) System.out.println("Failed to read sensor image file "+sourcePath);
        		continue;
        	}
    		if (debugLevel>0) System.out.println("Read "+sourcePath+" as a sensor image");
    		if (this.lanczos==null)   generateLanczosStack();
    		if (this.oversampled!=sourceImageScale) generateLanczosStack(
            		this.lanczosA,
            		sourceImageScale,
            		this.binsPerHalfPixel);
    		if (impSrc.getType()!=ImagePlus.COLOR_RGB) {
        		String msg="Not yet implemented for this image type, only RGB24 is supported";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
    		}
    		ImagePlus impOverlap=interSensor.resampleToOverlapRGB24(
    				false, // it is plane, not overlap of a pair of images
    				impSrc,
    				channel,
    				sourceImageScale,
        			this.lanczos,
        			this.binsPerHalfPixel,
        			maxThreads);
    		if (impOverlap==null){
        		String msg="Failed to apply plane map to "+sourcePath;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
    		}
    		if (debugLevel>1){
    			impOverlap.getProcessor().resetMinAndMax(); // imp_psf will be reused
    			impOverlap.show();
    		}
    		//resultDirectory
    		if (saveTiff){
    			String outPath=((resultDirectory.length()<1)?"":(resultDirectory+Prefs.getFileSeparator()))+impOverlap.getTitle()+".tiff";
    			try {
    				(new EyesisTiff()).saveTiff(
    						impOverlap,
    						outPath,
    						3, // float32
    						255.0,
    						true,
    						debugLevel);
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
    		}
    	}
    }

    public ImagePlus applyPlaneMap(
    		ImagePlus impSrc,
    		int channel,
    		int sourceImageScale, // 2.0
    		int maxThreads,
    		int debugLevel
    		){
    	InterSensor interSensor=this.sensors[channel].interSensor;
    		if (this.lanczos==null)   generateLanczosStack();
    		if (this.oversampled!=sourceImageScale) generateLanczosStack(
            		this.lanczosA,
            		sourceImageScale,
            		this.binsPerHalfPixel);
    		if (impSrc.getType()!=ImagePlus.COLOR_RGB) {
        		String msg="Not yet implemented for this image type, only RGB24 is supported";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
    		}
    		ImagePlus impPlane=interSensor.resampleToOverlapRGB24(
    				false, // it is plane, not overlap of a pair of images
    				impSrc,
    				channel,
    				sourceImageScale,
        			this.lanczos,
        			this.binsPerHalfPixel,
        			maxThreads);
    		if (impPlane==null){
        		String msg="Failed to apply plane map to channel "+channel;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
    		}
    		if (debugLevel>2){
    			impPlane.getProcessor().resetMinAndMax(); // imp_psf will be reused
    			impPlane.show();
    		}
    		return impPlane;
    }



    public void generateAndSaveInterSensorMaps(
    		String directory,
    		String fileNameFormat, //prefix%02d-%02dsuffix
    		double alphaThreshold,
    		double overlapThreshold,
    		double angleSizeX,
    		double angleSizeY,
    		int overlapExtraPixels,
    		boolean updateStatus,
    		int debugLevel
    		){
        int [][] pairs =findSensorPairs(
        		alphaThreshold,
        		overlapThreshold);
        if (debugLevel>0) {
        	System.out.println("Detected "+pairs.length+" overlapping sensor pairs");
        	if (debugLevel>1){
                for (int numPair=0;numPair<pairs.length;numPair++){
                	 System.out.println("Sensor pair "+(numPair+1)+ " ( of "+pairs.length+") - sensor "+
                     		pairs[numPair][0]+" and sensor "+pairs[numPair][1]);
                }
        	}
        }
        for (int numPair=0;numPair<pairs.length;numPair++){
            if (debugLevel>1) System.out.println("Processing sensor pair "+(numPair+1)+ " ( of "+pairs.length+") - sensor "+
            		pairs[numPair][0]+" and sensor "+pairs[numPair][1]);
            String title=String.format(fileNameFormat,pairs[numPair][0],pairs[numPair][1]);
        	ImagePlus imp=interSensorOverlapMap(
        			pairs[numPair][0],
        			pairs[numPair][1],
            		angleSizeX,
            		angleSizeY,
            		overlapExtraPixels,
            		title,
            		debugLevel);
        	if (imp!=null) {
        		if (debugLevel>2) {
        			imp.getProcessor().resetMinAndMax(); // imp_psf will be reused
        			imp.show();
        		}

        		String path=((directory!=null) && (directory.length()>0)?(directory+Prefs.getFileSeparator()):"")+title;
        		FileSaver fs=new FileSaver(imp);
        		String msg="Saving ovelap map for sensors "+pairs[numPair][0]+" and "+pairs[numPair][1]+": "+path;
        		if (updateStatus) IJ.showStatus(msg);
        		if (debugLevel>0) System.out.println(msg);
        		fs.saveAsTiffStack(path);
        	} else { // imp==null
            	if (debugLevel>0){
                    	 System.out.println("No overlap detected for sensors "+pairs[numPair][0]+" and "+pairs[numPair][1]);
            	}
        	}
        }
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

    // TODO: Make NaN for projectionPixelSize,projectionCenterX,projectionCenterY and 0 for projectionWidth - then find automatically
    // also save as parameters
    public ImagePlus getPlaneToSensorsMap(
    		int [] channels,
    		double projectionElevation,
    		double projectionYaw,
    		double projectionRoll,
    		double projectionPixelSize,
    		int projectionWidth,
    		int projectionHeight,
    		double projectionCenterX,
    		double projectionCenterY,
    		double nominalHorizontalDisparity, // nominal distance between horizontal cameras
    		String title,
    		int debugLevel
    		){

        generateOverlapBits(0.1, channels ); // sets  this.panoDegreesPerPixel among other things
    	// determine projectionPixelSize (should be the same for all sensors
    	if (Double.isNaN(projectionPixelSize) || (projectionPixelSize<=0.0)) {
    		projectionPixelSize=Math.PI/180.0*this.panoDegreesPerPixel;
    		if (debugLevel>0) System.out.println("getPlaneToSensorsMap(): Using projectionPixelSize="+projectionPixelSize);

//        	SensorData.EquirectangularMap erm=this.sensors[channels[0]].equirectangularMap;
//        	projectionPixelSize=erm.degreesPerPixel*Math.PI/180;
    	}
    	// create lat/long pair for each pixel of the projection plane

        double projPsi=Math.PI/180*projectionRoll;
    	double projTheta=Math.PI/180*projectionElevation;
    	double projPhi=Math.PI/180*projectionYaw;
       	// rotate by - projPsi around projection plane z (normal)
    	double [][] aR10={
    			{ Math.cos(projPsi), Math.sin(projPsi), 0.0},
    			{-Math.sin(projPsi), Math.cos(projPsi), 0.0},
    			{      0.0,                 0.0,        1.0}};
    	Matrix R10=new Matrix(aR10,3,3);

       	// rotate by - projTheta around projection plane X
    	double [][] aR21={
    			{1.0,    0.0,               0.0        },
    			{0.0, Math.cos(projTheta), Math.sin(projTheta)},
    			{0.0,-Math.sin(projTheta), Math.cos(projTheta)}};
    	Matrix R21=new Matrix(aR21,3,3);
    	// rotate by -projPhi around projection plane Y
    	double [][] aR32={
    			{ Math.cos(projPhi), 0.0, Math.sin(projPhi)},
    			{       0.0,         1.0,      0.0        },
    			{-Math.sin(projPhi), 0.0, Math.cos(projPhi)}};
    	Matrix MA=((new Matrix(aR32,3,3)).times(R21)).times(R10);
    	float [][][] projLatLong = new float [projectionHeight][projectionWidth][2];

    	double [][] V={{0.0},{0.0},{1.0}};
    	Matrix MV=new Matrix(V,3,1);
    	for (int py=0;py<projectionHeight;py++) {
    		V[1][0]=-projectionPixelSize*(py-projectionCenterY);
    		for (int px=0;px<projectionWidth;px++){
        		V[0][0]=projectionPixelSize*(px-projectionCenterX);
        		double [] projPointUnity=unityVector((MA.times(MV)).getRowPackedCopy());
        		projLatLong[py][px][0]=(float) ((180.0/Math.PI)*Math.asin(projPointUnity[1])); // latitude v[1] (y) - up -90..+90
        		projLatLong[py][px][1]=(float) ((180.0/Math.PI)*Math.atan2(projPointUnity[0],projPointUnity[2])); // longitude v[0](x) right, V[2](z) to the target +/-180
    		}
    	}
    	if (debugLevel>2) {
    		int i0=0;
    		double [][] dbgVect=new double [3][projectionHeight*projectionWidth];
        	for (int py=0;py<projectionHeight;py++) {
        		V[1][0]=-projectionPixelSize*(py-projectionCenterY);
        		for (int px=0;px<projectionWidth;px++){
            		V[0][0]=projectionPixelSize*(px-projectionCenterX);
            		double [] projPointUnity=unityVector((MA.times(MV)).getRowPackedCopy());
            		dbgVect[0][i0]=projPointUnity[0];
            		dbgVect[1][i0]=projPointUnity[1];
            		dbgVect[2][i0]=projPointUnity[2];
            		i0++;
        		}
        	}
   			String [] titles0={"Vx","Vy","Vz"};
			(new ShowDoubleFloatArrays()).showArrays(
					dbgVect,
					projectionWidth,
					projectionHeight,
					true,
					"vectors",
					titles0);
    		double [][] dbgLatLong=new double [2][projectionHeight*projectionWidth];
    		int i=0;
        	for (int py=0;py<projectionHeight;py++) {
        		for (int px=0;px<projectionWidth;px++){
        			dbgLatLong[0][i]=projLatLong[py][px][0];
        			dbgLatLong[1][i]=projLatLong[py][px][1];
        			i++;
        		}
        	}
   			String [] titles={"lat","long"};
			(new ShowDoubleFloatArrays()).showArrays(
					dbgLatLong,
					projectionWidth,
					projectionHeight,
					true,
					"LatLong_map",
					titles);
    	}

    	float [][] pixels=new float [3*channels.length][];
    	String [] titles=new String [3*channels.length];
    	for (int iSens=0;iSens<channels.length;iSens++){
    		int chn=channels[iSens];
    		float [][] channelPixels=mapRectangleToSensor(
            		chn,
            		projLatLong, //(=/-90, +/- 180)
            		debugLevel);
    		for (int i=0;i<3;i++)pixels[3*iSens+i]=channelPixels[i];
    		titles[3*iSens+0]="Channel "+chn+" pX";
    		titles[3*iSens+1]="Channel "+chn+" pY";
    		titles[3*iSens+2]="Channel "+chn+" alpha";
    	}
        ImageStack stack=new ImageStack(projectionWidth,projectionHeight);
        for (int i=0;i<pixels.length;i++) stack.addSlice(titles[i], pixels[i]);
        ImagePlus imp = new ImagePlus(title, stack);
        imp.setProperty("version", "1.1");
    	imp.setProperty("comment_projectionYaw", "Projection plane normal yaw in degrees, positive - CW looking from top");
    	imp.setProperty("comment_projectionElevation", "Projection plane normal elevation in degrees (up - positive)");
    	imp.setProperty("comment_projectionRoll", "Projection plane rotation around the normal, positive CW looking from the camera center");

    	imp.setProperty("projectionYaw", ""+projectionYaw);
    	imp.setProperty("projectionElevation", ""+projectionElevation);
    	imp.setProperty("projectionRoll", ""+projectionRoll);

    	double [][] projectedXYZ=new double [channels.length][];
        double [] projectedXYZCenter={0.0,0.0,0.0};
    	for (int iSens=0;iSens<channels.length;iSens++){
    		projectedXYZ[iSens]=projectPupil(
    	    		iSens,
    	    		projectionYaw,
    	    		projectionElevation,
    	    		projectionRoll);
    		for (int i=0;i<projectedXYZCenter.length;i++) projectedXYZCenter[i]+=projectedXYZ[iSens][i];
    	}
		for (int i=0;i<projectedXYZCenter.length;i++) projectedXYZCenter[i]/=channels.length;
        imp.setProperty("comment_channel","sensor number");
    	for (int iSens=0;iSens<channels.length;iSens++){
    		int chn=channels[iSens];
            imp.setProperty("channel"+(iSens+1),chn+"");
            imp.setProperty("path"+iSens,this.sensors[chn].path);
    	}
        imp.setProperty("comment_azimuth", "lens center azimuth, CW from top, degrees");
    	imp.setProperty("comment_radius", "lens center distance from the camera vertical axis, mm");
    	imp.setProperty("comment_height", "lens center vertical position from the head center, mm");
    	imp.setProperty("comment_heading", "lens heading - added to azimuth");
    	imp.setProperty("comment_elevation", "lens elevation from horizontal, positive - above horizon, degrees");
    	imp.setProperty("comment_roll",    "lens roll, positive - CW looking away from the camera");

    	imp.setProperty("comment_ppaXYZ", "Projection coordinates of each lens center on the common plane passing through the camera center");
    	imp.setProperty("comment_ppXYZ", "Projection coordinates of each lens center on the common plane passing through the sensors center of gravity");


    	for (int iSens=0;iSens<channels.length;iSens++){
    		int chn=channels[iSens];
            imp.setProperty("azimuth"+(iSens+1),  ""+this.sensors[chn].azimuth);
            imp.setProperty("radius"+(iSens+1),   ""+this.sensors[chn].radius);
            imp.setProperty("height"+(iSens+1),   ""+this.sensors[chn].height);
            imp.setProperty("heading"+(iSens+1),  ""+this.sensors[chn].phi);
            imp.setProperty("elevation"+(iSens+1),""+this.sensors[chn].theta);
            imp.setProperty("roll"+   (iSens+1),  ""+this.sensors[chn].psi);

            imp.setProperty("ppaX"+(iSens+1),""+projectedXYZ[chn][0]);
            imp.setProperty("ppaY"+(iSens+1),""+projectedXYZ[chn][1]);
            imp.setProperty("ppaZ"+(iSens+1),""+projectedXYZ[chn][2]);

            imp.setProperty("ppX"+(iSens+1),""+(projectedXYZ[chn][0]-projectedXYZCenter[0]));
            imp.setProperty("ppY"+(iSens+1),""+(projectedXYZ[chn][1]-projectedXYZCenter[1]));
            imp.setProperty("ppZ"+(iSens+1),""+(projectedXYZ[chn][2]-projectedXYZCenter[2]));
    	}
// save image plane vectors
    	imp.setProperty("comment_projectionCenter", "Pixel X and Y corresponding to the center of the image plane (closest to the camera center)");
    	imp.setProperty("projectionCenterX",  ""+projectionCenterX); // currently integer
    	imp.setProperty("projectionCenterY",  ""+projectionCenterY); // currently integer

    	imp.setProperty("comment_projectionPixelSize", "Projection pixel size relative to the distance from the camera center to the image plane");
    	imp.setProperty("projectionPixelSize",  ""+projectionPixelSize);

    	double [][] VZ={{0.0},{0.0},{1.0}};
    	double [][] VX={{1.0},{0.0},{0.0}};

		double [] imagePlane= unityVector((MA.times(new Matrix(VZ))).getRowPackedCopy());
		double [] imagePlaneX=unityVector((MA.times(new Matrix(VX))).getRowPackedCopy());

    	imp.setProperty("comment_imagePlane", "X,Y,Z components of the unity vector normal to the image plane. Plane X axis is defined by provection of the interVector");
    	imp.setProperty("imagePlaneX",  ""+imagePlane[0]);
    	imp.setProperty("imagePlaneY",  ""+imagePlane[1]);
    	imp.setProperty("imagePlaneZ",  ""+imagePlane[2]);

    	imp.setProperty("comment_interVector", "X,Y,Z components of the vector from camera1 to camera2 centers, in camera coord system, mm");
    	imp.setProperty("interVectorX",  ""+imagePlaneX[0]);
    	imp.setProperty("interVectorY",  ""+imagePlaneX[1]);
    	imp.setProperty("interVectorZ",  ""+imagePlaneX[2]);
    	// data for triclops viewer
    	imp.setProperty("comment_dispScales", "Disparity data for triclops viewer");
    	for (int iSens=0;iSens<channels.length;iSens++){
    		int chn=channels[iSens];
            imp.setProperty("dispScales_"+iSens+"_x",""+(-1.0*(projectedXYZ[chn][0]-projectedXYZCenter[0])/nominalHorizontalDisparity));
            imp.setProperty("dispScales_"+iSens+"_y",""+(     (projectedXYZ[chn][1]-projectedXYZCenter[1])/nominalHorizontalDisparity));
    	}

/* from older calibration, manual:
    	var dispScales_0_x = 0.0103208521;
    	var dispScales_0_y = 0.5763058147;
    	var dispScales_1_x = -0.5059333518;
    	var dispScales_1_y = -0.2864608745;
    	var dispScales_2_x = 0.4956124998;
    	var dispScales_2_y = -0.2898449402;
*/

    	(new JP46_Reader_camera(false)).encodeProperiesToInfo(imp);
    	return imp;
    }

    /**
     * Suggest common projection plane in camera coordinates perpendicular to average view vector,
     *  X-axis along projection of channlel1-> channel2 vector
     * @param channels array of channels (sensors) to use
     * @param channel1 right sub-camera
     * @param channel2 left sub-camera
     * @param debugLevel debug level
     * @param strictHorizontalDisparity - set roll to make two cameras have pure horizontal disparity (false - minimize overall roll)
     * @return (Yaw, elevation, roll) of the projection plane in degrees
     */

    public double [] suggestProjectionPlaneYawElevationRoll(
    		int [] channels,
    		int channel1,
    		int channel2,
    		boolean strictHorizontalDisparity,
    		int debugLevel
		){
    	int debugThreshold=1;
    	double [][] VZ={{0.0},{0.0},{1.0}};
    	double [][] VX={{1.0},{0.0},{0.0}};
    	double [][] VY={{0.0},{1.0},{0.0}}; // up?
    	double [] averageViewVector={0.0,0.0,0.0};
    	double [] averageViewXVector={0.0,0.0,0.0};
    	for (int iSens=0;iSens<channels.length;iSens++){
    		int chn=channels[iSens];
            double psi=Math.PI/180*this.sensors[chn].psi;
        	double theta=Math.PI/180*this.sensors[chn].theta; // radians
        	double phi=Math.PI/180*(this.sensors[chn].azimuth+this.sensors[chn].phi); // radians

        	if (debugLevel>=debugThreshold) {
        		System.out.println("iSens="+iSens+", chn="+chn);
        		System.out.println("this.sensors["+chn+"].theta"+this.sensors[chn].theta+", theta="+theta);
        		System.out.println("this.sensors["+chn+"].azimuth"+this.sensors[chn].azimuth+
        				" this.sensors["+chn+"].phi"+this.sensors[chn].phi+" (sum="+
        				(this.sensors[chn].azimuth+this.sensors[chn].phi)+
        				", phi="+phi);
        	}

           	// rotate by - projPsi around projection plane z (normal)
        	double [][] aR10={
        			{ Math.cos(psi), Math.sin(psi), 0.0},
        			{-Math.sin(psi), Math.cos(psi), 0.0},
        			{      0.0,           0.0,      1.0}};
        	Matrix R10=new Matrix(aR10,3,3);
           	// rotate by - projTheta around projection plane X
        	double [][] aR21={
        			{1.0,    0.0,               0.0        },
        			{0.0, Math.cos(theta), Math.sin(theta)},
        			{0.0,-Math.sin(theta), Math.cos(theta)}};
        	Matrix R21=new Matrix(aR21,3,3);
        	// rotate by -projPhi around projection plane Y
        	double [][] aR32={
        			{ Math.cos(phi), 0.0, Math.sin(phi)},
        			{       0.0,         1.0,      0.0        },
        			{-Math.sin(phi), 0.0, Math.cos(phi)}};
        	Matrix MA=(new Matrix(aR32,3,3)).times(R21.times(R10));
        	if (debugLevel>=debugThreshold) {
        		System.out.println("R21=");
        		R21.print(10,5);
        		System.out.println("R32=");
        		(new Matrix(aR32,3,3)).print(10,5);
        		System.out.println("MA=");
        		MA.print(10,5);
        	}
        	double [] sensorViwewVector= unityVector((MA.times(new Matrix(VZ))).getRowPackedCopy());
        	double [] sensorViwewXVector=unityVector((MA.times(new Matrix(VX))).getRowPackedCopy());
        	for (int i=0;i<averageViewVector.length;i++){
            	if (debugLevel>=debugThreshold) {
            		System.out.println("Sensor view vector=("+sensorViwewVector[0]+","+sensorViwewVector[1]+","+sensorViwewVector[2]+")");
            		System.out.println("Sensor view X vector=("+sensorViwewXVector[0]+","+sensorViwewXVector[1]+","+sensorViwewXVector[2]+")");
            	}
        		averageViewVector[i]+= sensorViwewVector[i];
        		averageViewXVector[i]+=sensorViwewXVector[i];
        	}
    	}
    	double [] planeZ=unityVector(averageViewVector);
    	double [] planeX=unityVector(crossProduct3(planeZ,crossProduct3(averageViewXVector,planeZ)));
    	double [] planeXM={-planeX[0],-planeX[1],-planeX[2]};
    	if (debugLevel>=debugThreshold) {
    		System.out.println("suggestProjectionPlaneYawElevationRoll(channels[],"+channel1+","+channel2+")");
    		System.out.println("planeZ=("+planeZ[0]+","+planeZ[1]+","+planeZ[2]+")");
    		System.out.println("planeX=("+planeX[0]+","+planeX[1]+","+planeX[2]+")");
    		System.out.println("planeXM=("+planeXM[0]+","+planeXM[1]+","+planeXM[2]+")");
    	}
    	double [] v1=lensCenterVector(channel1);
    	double [] v2=lensCenterVector(channel2);
    	double [] interVector={0.0,0.0,0.0};
    	for (int i=0;i<interVector.length;i++){
    		interVector[i]=v2[i]-v1[i];
    	}
    	double [] interProjectionUnity=unityVector(crossProduct3(planeZ,crossProduct3(interVector,planeZ)));
    	double [] horProjectionUnity=unityVector(crossProduct3(planeZ,new Matrix(VY).getRowPackedCopy())); //horizontal in projection plane
    	double [] vertProjectionUnity=unityVector(crossProduct3(horProjectionUnity, planeZ)); //"vertical" (closest to up) in projection plane
//    	double cosPPRoll=dotProduct(strictHorizontalDisparity?interProjectionUnity:planeXM,horProjectionUnity);
//    	double sinPPRoll=dotProduct(strictHorizontalDisparity?interProjectionUnity:planeXM,vertProjectionUnity);
    	double cosPPRoll=dotProduct(strictHorizontalDisparity?interProjectionUnity:planeX,horProjectionUnity);
    	double sinPPRoll=dotProduct(strictHorizontalDisparity?interProjectionUnity:planeX,vertProjectionUnity);
    	double ppRoll=180/Math.PI*Math.atan2(sinPPRoll, cosPPRoll);
    	double ppYaw= 180/Math.PI*Math.atan2(planeZ[0], planeZ[2]); // x/z
    	double ppElevation= 180/Math.PI*Math.atan2(planeZ[1], Math.sqrt(planeZ[0]*planeZ[0]+planeZ[2]*planeZ[2])); // x/z
    	double [] result={ppYaw,ppElevation,ppRoll};
    	if (debugLevel>=debugThreshold) {
    		System.out.println("interProjectionUnity=("+interProjectionUnity[0]+","+interProjectionUnity[1]+","+interProjectionUnity[2]+")");
    		System.out.println("horProjectionUnity=("+horProjectionUnity[0]+","+horProjectionUnity[1]+","+horProjectionUnity[2]+")");
    		System.out.println("vertProjectionUnity=("+vertProjectionUnity[0]+","+vertProjectionUnity[1]+","+vertProjectionUnity[2]+")");

    		System.out.println("cosPPRoll="+cosPPRoll+" sinPPRoll="+sinPPRoll);
    		System.out.println("ppYaw="+ppYaw+" ppElevation="+ppElevation+" ppRoll="+ppRoll);
    	}
    	return result;
    }
    /**
     * Project lens center (entrance pupil) on the projection plane (going through the camera center)
     * @param channel subcamera number
     * @param ppYaw   projection plane normal yaw in degrees, positive - CW looking from top
     * @param ppElevation projection plane normal elevation in degrees (up - positive)
     * @param ppRoll projection plane rotation around the normal, positive CW looking from the camera center
     * @param debugLevel debug level
     * @return triplet of coordinates aligned to the projection plane (z - away from the camera center)
     */
    double [] projectPupil(
    		int channel,
    		double ppYaw,
    		double ppElevation,
    		double ppRoll){
    	double psi=Math.PI/180*ppRoll;
    	double theta=Math.PI/180*ppElevation; // radians
    	double phi=Math.PI/180*ppYaw; // radians
    	double [] v=lensCenterVector(channel);
    	double [][] V={{v[0]},{v[1]},{v[2]}};

     	/*1) rotate by -phi) around
     	 | cos(phi)    0  -sin(phi)   |
     	 |     0       1         0    |
     	 |  sin(phi)   0   cos(phi)   |
     	 */
     	double [][] aR0={
     			{Math.cos(phi),0.0,-Math.sin(phi)},
     			{0.0,          1.0, 0.0},
     			{Math.sin(phi),0.0, Math.cos(phi)}};
     	Matrix R0=new Matrix(aR0);

     	/*
     	2) rotate by theta around
     	|    1         0         0        |
     	|    0    cos(theta)  -sin(theta) |
     	|    0    sin(theta)   cos(theta) |
     	*/
     	double [][] aR1={
     			{1.0,  0.0,             0.0},
     			{0.0,  Math.cos(theta),-Math.sin(theta)},
     			{0.0,  Math.sin(theta), Math.cos(theta)}};
     	Matrix R1=new Matrix(aR1);
       	// 3) rotate by  psi around projection plane z (normal)
    	double [][] aR2={
    			{ Math.cos(psi), -Math.sin(psi), 0.0},
    			{ Math.sin(psi),  Math.cos(psi), 0.0},
    			{ 0.0,            0.0,           1.0}};
    	Matrix R2=new Matrix(aR2);
    	double [] xyz= R2.times(R1.times(R0.times(new Matrix(V)))).getRowPackedCopy();
    	return xyz;
    }

//    double [] crossProduct3(double [] a,double [] b){
    /**
     * Calculate sensor entrance pupil coordinates in camera coordinate system
     * @param chn sensor number
     * @return (x,y,z)
     */
    public double [] lensCenterVector(int chn){
    	double theta=Math.PI/180*this.sensors[chn].theta; // radians
    	double phi=Math.PI/180*this.sensors[chn].phi; // radians
    	double azimuth=Math.PI/180*this.sensors[chn].azimuth; // radians
    	/* 0) Translate by distance to entrance pupil (lens center)
     	| Xc0 |   | 0                     |   |Xc|
     	| Yc0 | = | 0                     | + |Yc|
     	| Zc0 |   | entrancePupilForward  |   |Zc|
     	*/
     	double [][] aT0={{0.0},{0.0},{this.sensors[chn].entrancePupilForward}};
     	Matrix T0=new Matrix(aT0);
     	/*
     	2) rotate by - theta around C1X:Vc2= R2*Vc1
     	| Xc2 |   |    1         0         0        |   |Xc1|
     	| Yc2 | = |    0    cos(theta)   sin(theta) | * |Yc1|
     	| Zc2 |   |    0   -sin(theta)   cos(theta) |   |Zc1|
     	*/
     	double [][] aR2={
     			{1.0,0.0,0.0},
     			{0.0,Math.cos(theta),Math.sin(theta)},
     			{0.0,-Math.sin(theta),Math.cos(theta)}};
     	Matrix R2=new Matrix(aR2);
     	/*
     	3) rotate by -(azimuth+phi) around C2Y:Vc3= R3*Vc2
     	| Xc3 |   | cos(azimuth+phi)    0   sin(azimuth+phi)   |   |Xc2|
     	| Yc3 | = |     0               1         0            | * |Yc2|
     	| Zc3 |   | -sin(azimuth+phi)   0   cos(azimuth+phi)   |   |Zc2|
     	 */
     	double [][] aR3={
     			{Math.cos(phi+azimuth),0.0,Math.sin(phi+azimuth)},
     			{0.0,1.0,0.0},
     			{-Math.sin(phi+azimuth),0.0,Math.cos(phi+azimuth)}};
     	Matrix R3=new Matrix(aR3);
     	/*
     	4) Now axes are aligned, just translate to get to eyesis coordinates: Vey= T1+Vc3
     	| Xey |   |      r * sin (azimuth)       |   |Xc3|
     	| Yey | = | height+centerAboveHorizontal | + |Yc3|
     	| Zey |   |      r * cos (azimuth)       |   |Zc3|
     	 */
     	double [][] aT1={
     			{this.sensors[chn].radius*Math.sin(azimuth)},
     			{this.sensors[chn].height},
     			{this.sensors[chn].radius*Math.cos(azimuth)}};
     	Matrix T1=new Matrix(aT1);
     // MA=R3*R2;
     // MB=T1+R3*R2*T0;
//     	 Matrix MA=R3.times(R2);
//     	 Matrix MB=T1.plus(R3.times(R2.times(T0)));
     	 double [] aB=T1.plus(R3.times(R2.times(T0))).getRowPackedCopy();
     	 return aB;
    }


    /**
     * Create a 6-slice map of the flat rectangular intersection area of two sensors
     * to the sensor pixel coordinates. The map is centered around the mid-point of the two
     * sensors on the sphere, horizontal dimension (x) from sensor1 to sensor 2.
     * 	this.panoDegreesPerPixel, this.panoLongitudeLeft, this.panoLongitudeRight,
     *  this.panoLatitudeTop and this.panoLatitudeBottom should be already set (i.e. by
     *  overlapPairsAreas() ).
     * @param channel1 first sensor channel number
     * @param channel2 second sensor channel number
     * @param angleSizeX overlap dimension in degrees (usually narrow)
     * @param angleSizeY overlap dimension in degrees (usually high)
     * @param extraPixels trim result to have only this number of non-overlapping pixels along X-axis
     * @param debugLevel debug level
     * @return image with 6 slices - {pixelX1, pixelY1, alpha1, pixelX2, pixelY2, alpha2}
     */
    public ImagePlus interSensorOverlapMap(
    		int channel1,
    		int channel2,
    		double angleSizeX,
    		double angleSizeY,
    		int extraPixels,
    		String title,
    		int debugLevel
    		){
//    	double [] interDistance={0.0};
    	double [][] vectors = new double [2][];
        float [][][] latLong= mapInterSensorToLatLong(
        		channel1,
        		channel2,
        		angleSizeX,
        		angleSizeY,
        		vectors,
        		debugLevel);
        int height=latLong.length;
        int width=latLong[0].length;
        float [][] pixels1=mapRectangleToSensor(
        		channel1,
        		latLong, //(=/-90, +/- 180)
        		debugLevel);
        float [][] pixels2=mapRectangleToSensor(
        		channel2,
        		latLong, //(=/-90, +/- 180)
        		debugLevel);
        // find first and last non-empty lines
        // Trim  result arrays
    	int minX=width,maxX=-1,minY=height,maxY=-1;
    	for (int iy=0;iy<height;iy++) for (int ix=0;ix<width;ix++){
    		int index=ix+width*iy;
    		if ((pixels1[2][index]>0.0) && (pixels2[2][index]>0.0)){ // overlap
    			if (ix>maxX) maxX=ix;
    			if (ix<minX) minX=ix;
    			if (iy>maxY) maxY=iy;
    			if (iy<minY) minY=iy;
    		}
    	}
    	if (minX>maxX){
        	if (debugLevel>0){
        		System.out.println ("interSensorOverlapMap(): No overlap width="+width+" height="+height+
        				" minX="+minX+" maxX="+maxX+" minY="+minY+" maxY="+maxY);
        	}
        	return null;
    	}
    	if (debugLevel>1){
    		System.out.println ("interSensorOverlapMap(): width="+width+" height="+height+
    				" minX="+minX+" maxX="+maxX+" minY="+minY+" maxY="+maxY);
    	}
    	minX-=extraPixels;
    	if (minX<0) minX=0;
    	maxX+=extraPixels;
    	if (maxX>=width) maxX=width-1;
    	int projectionCenterX=(width-1)/2-minX; // width is odd, symmetrical around the center
    	int projectionCenterY=(height-1)/2-minY; // height is odd, symmetrical around the center
    	int oWidth=maxX-minX+1;
    	int oHeight=maxY-minY+1;
    	float [][] pixels=new float [6][oWidth*oHeight];
    	for (int iy=0;iy<oHeight;iy++) for (int ix=0;ix<oWidth;ix++){
    		int iIndex= (iy+minY)*width+(ix+minX);
    		int oIndex= iy*oWidth+ ix;
            for (int n=0;n<3;n++){
            	pixels[n  ][oIndex]=  pixels1[n][iIndex];
            	pixels[n+3][oIndex]=  pixels2[n][iIndex];
            	pixels[n+3]=pixels2[n]; //????
            }
    	}

        String [] titles={
        		"Channel "+channel1+" pX",
        		"Channel "+channel1+" pY",
        		"Channel "+channel1+" alpha",
        		"Channel "+channel2+" pX",
        		"Channel "+channel2+" pY",
        		"Channel "+channel2+" alpha"};
        ImageStack stack=new ImageStack(oWidth,oHeight);
        for (int i=0;i<pixels.length;i++) stack.addSlice(titles[i], pixels[i]);
        ImagePlus imp = new ImagePlus(title, stack);
        imp.setProperty("comment_channel1","first sensor number");
        imp.setProperty("channel1",channel1+"");
        imp.setProperty("comment_channel2","second sensor number");
        imp.setProperty("channel2",channel2+"");
        imp.setProperty("path1",this.sensors[channel1].path);
        imp.setProperty("path2",this.sensors[channel2].path);

    	imp.setProperty("comment_azimuth", "lens center azimuth, CW from top, degrees");
    	imp.setProperty("azimuth1", ""+this.sensors[channel1].azimuth);
    	imp.setProperty("azimuth2", ""+this.sensors[channel2].azimuth);
    	imp.setProperty("comment_radius", "lens center distance from the camera vertical axis, mm");
    	imp.setProperty("radius1",  ""+this.sensors[channel1].radius);
    	imp.setProperty("radius2",  ""+this.sensors[channel2].radius);
    	imp.setProperty("comment_height", "lens center vertical position from the head center, mm");
    	imp.setProperty("height1",  ""+this.sensors[channel1].height);
    	imp.setProperty("height2",  ""+this.sensors[channel2].height);
    	imp.setProperty("comment_heading", "lens heading - added to azimuth");
    	imp.setProperty("heading1",  ""+this.sensors[channel1].phi);
    	imp.setProperty("heading2",  ""+this.sensors[channel2].phi);
    	imp.setProperty("comment_elevation", "lens elevation from horizontal, positive - above horizon, degrees");
    	imp.setProperty("elevation1",  ""+this.sensors[channel1].theta);
    	imp.setProperty("elevation2",  ""+this.sensors[channel2].theta);
// save image plane vectors
    	imp.setProperty("comment_projectionCenter", "Pixel X and Y corresponding to the center of the image plane (closest to the camera center)");
    	imp.setProperty("projectionCenterX",  ""+projectionCenterX); // currently integer
    	imp.setProperty("projectionCenterY",  ""+projectionCenterY); // currently integer

    	imp.setProperty("comment_projectionPixelSize", "Projection pixel size relative to the distance from the camera center to the image plane");
    	imp.setProperty("projectionPixelSize",  ""+Math.PI/180.0*this.panoDegreesPerPixel);

    	imp.setProperty("comment_imagePlane", "X,Y,Z components of the unity vector normal to the image plane. Plane X axis is defined by projection of the interVector");
    	imp.setProperty("imagePlaneX",  ""+vectors[0][0]);
    	imp.setProperty("imagePlaneY",  ""+vectors[0][1]);
    	imp.setProperty("imagePlaneZ",  ""+vectors[0][2]);

    	imp.setProperty("comment_interVector", "X,Y,Z components of the vector from camera1 to camera2 centers, in camera coord system, mm");
    	imp.setProperty("interVectorX",  ""+vectors[1][0]);
    	imp.setProperty("interVectorY",  ""+vectors[1][1]);
    	imp.setProperty("interVectorZ",  ""+vectors[1][2]);

    	(new JP46_Reader_camera(false)).encodeProperiesToInfo(imp);
    	return imp;
    }

    /**
     * Map sensor overlap flat rectangular area to coordinates on the sphere. Adding multiple of 360 to parameters is OK here,
     * result is -90<=..<=+90 for latitude and -180>=..>+180 for longitude
     * @param channel1 - number of the first camera
     * @param channel1 - number of the second camera
     * @param angleSizeX angular overlap size along the line connecting two cameras, in degrees
     * @param angleSizeY angular overlap size perpendicular to the line connecting two cameras, in degrees
     * @param vectors  should be initialized as double [2][]. Returns 2 vectors:
     *  first - unity vector, normal to the image plane, second - vector from camera1 to camera 2, mm.
     *  Image plane axis X is a projection of this vector on the image plane
     * @param debugLevel debug level
     * @return [pixels perpendicular to line between cameras][pixels in the direction from camera 1 to camera 2]{lat, long}
     */

    public float [][][] mapInterSensorToLatLong(
    		int channel1,
    		int channel2,
    		double angleSizeX,
    		double angleSizeY,
    		double [][] vectors, // should be initialized as double[2][]
    		int debugLevel
    ){


		double lat1= this.sensors[channel1].theta; //latitude (elevation) of the first camera optical axis, degrees
		double long1 = this.sensors[channel1].azimuth+this.sensors[channel1].phi; //longitude (heading+azimuth) of the first camera (CW from top) adding multiple 360 is OK here
		double lat2 = this.sensors[channel2].theta; //latitude of the second camera
		double long2= this.sensors[channel2].azimuth+this.sensors[channel2].phi; // longitude of the second camera, dding multiple 360 is OK here
		double degreesPerPixel= this.panoDegreesPerPixel; //angular size of a pixel on equator, in degrees


    	// camera coordinates y - up, to the target (from the camera at zero angle/zero heading), x - right (looking to the target)
    	// converting sensor vector {0,0,1} to the camera coordinates. skipping roll

    	double theta1= lat1* Math.PI/180;
    	double phi1=   long1*Math.PI/180;
    	double theta2= lat2* Math.PI/180;
    	double phi2=   long2*Math.PI/180;
    	double [][] aSensorAxis={{0.0},{0.0},{1.0}};
    	Matrix mSensorAxis=new Matrix(aSensorAxis);

    	// rotate by - theta1 around sensor1 x:Vc2= R2*Vc1
    	double [][] aR21={
    			{1.0,    0.0,               0.0        },
    			{0.0, Math.cos(theta1), Math.sin(theta1)},
    			{0.0,-Math.sin(theta1), Math.cos(theta1)}};
    	Matrix R21=new Matrix(aR21);
    	// rotate by -(phi1) around C2Y:Vc3= R3*Vc2
    	double [][] aR31={
    			{ Math.cos(phi1), 0.0, Math.sin(phi1)},
    			{    0.0,         1.0,    0.0        },
    			{-Math.sin(phi1), 0.0, Math.cos(phi1)}};
    	Matrix MA1=(new Matrix(aR31)).times(R21);
    	double [] u1=(MA1.times(mSensorAxis)).getColumnPackedCopy(); // unity vector in the direction of the first camera axis

    	double [][] aR22={
    			{1.0,    0.0,               0.0        },
    			{0.0, Math.cos(theta2), Math.sin(theta2)},
    			{0.0,-Math.sin(theta2), Math.cos(theta2)}};
    	Matrix R22=new Matrix(aR22);
    	// rotate by -(phi1) around C2Y:Vc3= R3*Vc2
    	double [][] aR32={
    			{ Math.cos(phi2), 0.0, Math.sin(phi2)},
    			{    0.0,         1.0,     0.0       },
    			{-Math.sin(phi2), 0.0, Math.cos(phi2)}};
    	Matrix MA2=(new Matrix(aR32) ).times(R22);
    	double [] u2=(MA2.times(mSensorAxis)).getColumnPackedCopy(); // unity vector in the direction of the first camera axis

    	double [] rsltUZ={
    			u2[0]+u1[0],
    			u2[1]+u1[1],
    			u2[2]+u1[2]};
    	rsltUZ=unityVector(rsltUZ); // unity vector between the views of the two cameras (maybe - find center of the overlap?
    	if (debugLevel>1) System.out.println("mapInterSensorToLatLong(): rsltUZ={"+rsltUZ[0]+","+rsltUZ[1]+","+rsltUZ[2]+"}");
    	if (debugLevel>1) System.out.println("mapInterSensorToLatLong(): radius1="+this.sensors[channel1].radius+
    			", azimuth1="+this.sensors[channel1].azimuth+
    			", height1="+this.sensors[channel1].height);
    	if (debugLevel>1) System.out.println("mapInterSensorToLatLong(): radius2="+this.sensors[channel2].radius+
    			", azimuth2="+this.sensors[channel2].azimuth+
    			", height2="+this.sensors[channel2].height);
    	// lens centers locations:
    	double [] camera1XYZ={
    			this.sensors[channel1].radius*Math.sin(Math.PI/180.0* this.sensors[channel1].azimuth),
    			this.sensors[channel1].height,
    			this.sensors[channel1].radius*Math.cos(Math.PI/180.0* this.sensors[channel1].azimuth)};
    	double [] camera2XYZ={
    			this.sensors[channel2].radius*Math.sin(Math.PI/180.0* this.sensors[channel2].azimuth),
    			this.sensors[channel2].height,
    			this.sensors[channel2].radius*Math.cos(Math.PI/180.0* this.sensors[channel2].azimuth)};
    	if (debugLevel>1) System.out.println("mapInterSensorToLatLong(): camera1XYZ={"+camera1XYZ[0]+","+camera1XYZ[1]+","+camera1XYZ[2]+"}");
    	if (debugLevel>1) System.out.println("mapInterSensorToLatLong(): camera2XYZ={"+camera2XYZ[0]+","+camera2XYZ[1]+","+camera2XYZ[2]+"}");
    	double [] interCamerVector={
    			camera2XYZ[0]-camera1XYZ[0],
    			camera2XYZ[1]-camera1XYZ[1],
    			camera2XYZ[2]-camera1XYZ[2]};
    	if (debugLevel>1) System.out.println("mapInterSensorToLatLong(): interCamerVector={"+interCamerVector[0]+","+interCamerVector[1]+","+interCamerVector[2]+"}");
    	// here x,y,z - left!
    	double [] rsltUY=unityVector(crossProduct3(rsltUZ,interCamerVector));
//    	double [] rsltUX= crossProduct3(interCamerVector,rsltUY);
    	double [] rsltUX= crossProduct3(rsltUY,rsltUZ);
    	if (debugLevel>1) System.out.println("mapInterSensorToLatLong(): rsltUY={"+rsltUY[0]+","+rsltUY[1]+","+rsltUY[2]+"}");
    	if (debugLevel>1) System.out.println("mapInterSensorToLatLong(): rsltUX={"+rsltUX[0]+","+rsltUX[1]+","+rsltUX[2]+"}");
    	// calculate output dimensions
    	double outPixelSize=Math.PI/180.0*degreesPerPixel; // radius =1.0
    	double linearXHalfRange=Math.tan(angleSizeX*Math.PI/360.0)/outPixelSize; // using half-angle
    	double linearYHalfRange=Math.tan(angleSizeY*Math.PI/360.0)/outPixelSize; // using half-angle
    	int xHalfSize=(int) Math.round(linearXHalfRange);
    	int yHalfSize=(int) Math.round(linearYHalfRange);
    	float [][][] result = new float [2*yHalfSize+1][2*xHalfSize+1][2];
    	for (int iy=0;iy<=2*yHalfSize;iy++) for (int ix=0;ix<=2*xHalfSize;ix++){
//    		double y=(iy-yHalfSize)*outPixelSize;
    		double y=(yHalfSize-iy)*outPixelSize; // in pixel order - y-down
    		double x=(ix-xHalfSize)*outPixelSize;
    		double [] v={ // point on the result plane in camera coordinates
    				x*rsltUX[0]+y*rsltUY[0]+rsltUZ[0],
    				x*rsltUX[1]+y*rsltUY[1]+rsltUZ[1],
    				x*rsltUX[2]+y*rsltUY[2]+rsltUZ[2],
    		};
    		// convert to latitude/longitude
    		v=unityVector(v); // re-normalize to increase precision?
    		result[iy][ix][0]=(float) ((180.0/Math.PI)*Math.asin(v[1])); // latitude v[1] (y) - up -90..+90
    		result[iy][ix][1]=(float) ((180.0/Math.PI)*Math.atan2(v[0],v[2])); // longitude v[0](x) right, V[2](z) to the target +/-180
    	}

    	vectors[0]=rsltUZ;
    	vectors[1]=interCamerVector;

    	return result;
    }
    /**
     * Calculate cross-product of 2 vectors in 3d
     * @param a 3-element vector a
     * @param b 3-element vector b
     * @return 3-element cross product of vectors a and b
     */
    double [] crossProduct3(double [] a, double [] b){
    	double []cp3={
    			a[1]*b[2]-a[2]*b[1],
    			a[2]*b[0]-a[0]*b[2],
    			a[0]*b[1]-a[1]*b[0]};
    	return cp3;
    }

    /**
     * Calculate dot-product of 2 vectors
     * @param a vector a
     * @param b vector b
     * @return dot product of two vectors
     */
    double dotProduct(double [] a, double [] b){
    	double dp=0;
    	for (int i=0;i<a.length;i++) dp+=a[i]*b[i];
    	return dp;
    }

    /**
     * Make a unity vector by dividing each component by the vector length
     * @param a vector
     * @return unity vector in the direction of a
     */
    double [] unityVector(double [] a) {
    	double s=0.0;
    	for (int i=0;i<a.length;i++)s+=a[i]*a[i];
    	s=1.0/Math.sqrt(s);
    	double [] uvect=new double [a.length];
    	for (int i=0;i<a.length;i++) uvect[i]=s*a[i];
    	return uvect;
    }



    /**
     * Map rectangle of {lat,long} pairs to sensor pixels;
     * @param channel number of sensor (channel)
     * @param latlong 2-d array of lat, long pairs
     * @param debugLevel debug level
     * @return array of pixelx, pixely and alpha slices, indexed in scanline order
     */
    public float [][] mapRectangleToSensor(
    		int channel,
    		float [][][] latLong, //(=/-90, +/- 180)
    		int debugLevel){
    	SensorData.EquirectangularMap erm=this.sensors[channel].equirectangularMap;
    	Rectangle woi=erm.mapWOI;
    	if (debugLevel>1){
    		System.out.println("this.sensors["+channel+"].channel].equirectangularMap.maoWOI= {x="+
    				woi.x+", y="+woi.y+", width="+woi.width+", height="+woi.height+"}");
    		System.out.println("this.panoDegreesPerPixel="+this.panoDegreesPerPixel);
    		System.out.println("this.panoLatitudeTop="+this.panoLatitudeTop);
    		System.out.println("this.panoLongitudeLeft="+this.panoLongitudeLeft);
    	}
    	double minLat=this.panoLatitudeTop-(woi.y+woi.height)*this.panoDegreesPerPixel;
    	double maxLat=this.panoLatitudeTop-woi.y*this.panoDegreesPerPixel;
    	double minLong=this.panoLongitudeLeft+woi.x*this.panoDegreesPerPixel;
    	double maxLong=this.panoLongitudeLeft+(woi.x+woi.width)*this.panoDegreesPerPixel;
    	while (minLong<-180.0) minLong+=360.0;
    	while (minLong>=180.0) minLong-=360.0; // -180<=long<180
    	while (maxLong<-180.0) maxLong+=360.0;
    	while (maxLong>=180.0) maxLong-=360.0; // -180<=long<180
    	boolean rollover= (maxLong<=minLong);
    	if (debugLevel>1){
    		System.out.println("mapRectangleToSensor("+channel+",...): minLong="+minLong+", maxLong="+maxLong+
    				", minLat="+minLat+", maxLat="+maxLat+", rollover="+rollover);

    	}
    	float [][] pixelCoords =new float[3][latLong.length*latLong[0].length];
    	for (int n=0;n<pixelCoords.length;n++) for (int i=0;i<pixelCoords[0].length;i++) pixelCoords[n][i]=0.0f;
    	for (int iy=0;iy<latLong.length;iy++) for (int ix=0;ix<latLong[0].length;ix++){
    		double lat= latLong[iy][ix][0];
    		double lng= latLong[iy][ix][1];
        	if ((iy==(latLong.length/2)) && (ix==(latLong[0].length/2)) && (debugLevel>1)){
        		System.out.println("mapRectangleToSensor("+channel+",...): lat["+iy+"]["+ix+"]="+lat+", long["+iy+"]["+ix+"]="+lng);

        	}
    		if ((lat<minLat) || (lat>maxLat)) continue;
    		if (!rollover &&((lng<minLong) || (lng>maxLong))) continue;
    		if ( rollover &&((lng<minLong) && (lng>maxLong))) continue;
    		if (lng<minLong)lng+=360.0;
    		double apy=(maxLat - lat)/this.panoDegreesPerPixel;
    		double apx=(lng-minLong)/this.panoDegreesPerPixel;
        	if ((iy==(latLong.length/2)) && (ix==(latLong[0].length/2)) && (debugLevel>1)){
        		System.out.println("mapRectangleToSensor("+channel+",...): apy["+iy+"]["+ix+"]="+apy+", apx["+iy+"]["+ix+"]="+apx);
        	}
    		int iapy= (int) Math.floor(apy);
    		int iapx= (int) Math.floor(apx);
    		if ((iapy>= (woi.height-1)) || (iapy<0)) continue;
    		if ((iapx>= (woi.width -1)) || (iapx<0)) continue;
    		int index00=iapx+iapy*woi.width;
    		if (erm.partialMap[2][index00]<=0.0)  continue; // alpha is 0 - skip
			int indexX0=index00+1;
			int index0Y=index00+woi.width;
			int indexXY=index0Y+1;
    		double dapy=apy-iapy;
    		double dapx=apx-iapx;
        	if ((iy==(latLong.length/2)) && (ix==(latLong[0].length/2)) && (debugLevel>1)){
        		System.out.println("mapRectangleToSensor("+channel+",...): dapy["+iy+"]["+ix+"]="+dapy+", dapx["+iy+"]["+ix+"]="+dapx);
        	}

    		if ((dapx>0.0)                 && (erm.partialMap[2][indexX0]<=0.0)) continue;
    		if ((dapy>0.0)                 && (erm.partialMap[2][index0Y]<=0.0)) continue;
    		if (((dapy>0.0) || (dapx>0.0)) && (erm.partialMap[2][indexXY]<=0.0)) continue;

    		// bi-linear interpolate
    		int index=iy*latLong[0].length+ix;
    		for (int n=0;n<pixelCoords.length;n++) {
    			pixelCoords[n][index]= (float)(
    				(1-dapx)*(1-dapy)* erm.partialMap[n][index00]+
    				dapx *   (1-dapy)* erm.partialMap[n][indexX0]+
					(1-dapx)*   dapy * erm.partialMap[n][index0Y]+
					dapx   *    dapy * erm.partialMap[n][indexXY]);
    		}
    	}
    	return pixelCoords;
    }

    public class InterSensor{
    	public int numOverlapChannels=2;
    	public int debugLevel=1;
    	public int [] channel=new int[2];
    	public double [][] disparityScales; // for each channel - a pair of {scaleX, scaleY} or null if undefined;
    	public String path=null;
//    	String [] paths={null,null};
    	public double [] imagePlane=new double[3];  // X,Y,Z components of the unity vector normal to the image plane. Plane X axis is defined by projection of the interVector
//    	double [] interVector=null; //"X,Y,Z components of the vector from camera1 to camera2 centers, in camera coord system, mm"
    	public float [][] map; //{pixelx1,pixely1,alpha1,pixelX2,pixelY2,alpha2}[]  // may have now 3 (more) components
    	public double [][] overlapImages=null; // {alpha1, Y1, Cb1, Cr1, alpha2, Y2, Cb2, Cr2}[pixels]
    	public double [][] sobelY=null; // edge-detection ran on {Y1,Y2}, same format;
    	public boolean [][] booleanEdges=null;
    	public boolean [][] thinEdges=null;
    	public double preSobelSigma=1.4;
    	public double edgeSpreadOrtho=0.5;
    	public double edgeSpreadDiagonal=0.4;

    	public double edgeThresholdHigh=0.2;
    	public double edgeThresholdLow=0.08; //4;

    	public int debugXMin=144;
    	public int debugXMax=151;
    	public int debugYMin=65;
    	public int debugYMax=72;

    	public int mapWidth;
    	public int mapHeight;
    	public double [] projectionCenterXY=new double[2]; // pixel coordinates that correspond to the image plane normal vector (so far integer)
    	public double projectionPixelSize;
    	public int tilesSize; // size of largest (top) tiles, power of 2.  (FFT Size will be twice that
    	public int tilesX0;   // pixel X position of the top left corner of the top left tile
    	public int tilesY0;   // pixel Y position of the top left corner of the top left tile
    	public int tilesNHor; // number of tiles horizontally
    	public int tilesNVert;// number of tiles vertically
    	public double [][][][][] tiles0; // [layer][vert][hor][background,foreground]{position,strength,rmsY,rmsCb,rmsCr} // old, to be deleted
    	public double [][][][][][] tiles; // [side][self][layer][vert][hor][disparity]
    	public int minDisparity;
    	public int maxDisparity;

    	public double [][][][][] disparityMap=null; // [side][self][vert][hor][disparity]
    	public double [][][][]   linearFeatures=null; // [side][vert][hor][feature_parameter]
    	public int   disparityMapFFTSize;
    	public int   disparityMapOverlapFraction;

    	public float [][][] ambiguityMap=null; // [side][first pos, first contrast, second pos, second contrast][pixel
    	public float [][]   resolvedMap=null; //  [side][pixel]

    	//duplicates of the sensor data, needed later
    	public double [] azimuth=new double[2];
    	public double [] radius=new double[2];
    	public double [] height= new double[2];

    	public EdgesSegmeniting edgesSegmeniting=new EdgesSegmeniting();
    	public float [][] distanceFromEdges=null;
    	public int [][] edgeAreas=null;
 /*
    	public int [] channel=new int[2];
    	public double [][] disparityScales; // for each channel - a pair of {scaleX, scaleY} or null if undefined;

  */
    	public double [][] getDisparityScales(){
    		if (this.channel==null) return null; // even channels are not initialized
    		if ((this.disparityScales==null) || (this.disparityScales.length!=this.channel.length)){
    			this.disparityScales=new double [this.channel.length][];
    			for (int i=0;i<this.disparityScales.length;i++) this.disparityScales[i]=null;
    		}
    		return this.disparityScales;
    	}
    	private String[] generateRequiredProperties(int numChannels){
        	String [] requiredPropertiesCommon={
        			"imagePlaneX",
        			"imagePlaneY",
        			"imagePlaneZ",
        			"projectionCenterX",
        			"projectionCenterY",
        			"projectionPixelSize"
        	};
        	String [] requiredPropertiesSensor={
        			"channel",
        			"azimuth",
        			"radius",
        			"height"
        	};
        	String [] requiredProperties = new String [requiredPropertiesCommon.length+numChannels*requiredPropertiesSensor.length];
        	int index=0;
        	for (int i=0;i<requiredPropertiesCommon.length;i++)     requiredProperties[index++]=requiredPropertiesCommon[i];
        	for (int chn=0;chn<numChannels;chn++){
        		for (int i=0;i<requiredPropertiesSensor.length;i++) requiredProperties[index++]=requiredPropertiesSensor[i]+(chn+1);
        	}
    		return requiredProperties;
    	}
    	public InterSensor(ImagePlus imp, int debugLevel){
    		this.debugLevel=debugLevel;
    		int minNumPlanes=3; // was 6, trying for a single cameras
    		if (imp == null){
    			String msg="Sensor averlap area map image is null";
    			IJ.showMessage("Error",msg);
    			throw new IllegalArgumentException (msg);
    		}
    		if (imp.getStackSize()<minNumPlanes){
    			String msg="Expecting >="+minNumPlanes+" slices, got "+imp.getStackSize();
    			IJ.showMessage("Error",msg);
    			throw new IllegalArgumentException (msg);
    		}
    		if ((imp.getStackSize()%3)!=0){
    			String msg="Expecting multiple of 3 slices, got "+imp.getStackSize();
    			IJ.showMessage("Error",msg);
    			throw new IllegalArgumentException (msg);
    		}
    		this.numOverlapChannels=imp.getStackSize()/3;
    		this.channel=new int[this.numOverlapChannels];
    		this.projectionCenterXY=new double[2]; // pixel coordinates that correspond to the image plane normal vector (so far integer)
    		this.azimuth=new double[this.numOverlapChannels];
    		this.radius=new double[this.numOverlapChannels];
    		this.height= new double[this.numOverlapChannels];
    		String [] requiredProperties=generateRequiredProperties(this.numOverlapChannels);

    		// try to decode info if required properties are not set
    		if (imp.getProperty(requiredProperties[0])==null) (new JP46_Reader_camera(false)).decodeProperiesFromInfo(imp);

    		for (int i=0; i<requiredProperties.length;i++) if (imp.getProperty(requiredProperties[i])==null){
    			String msg="Required property "+requiredProperties[i]+" is not defined in "+imp.getTitle();
    			IJ.showMessage("Error",msg);
    			throw new IllegalArgumentException (msg);
    		}
    		ImageStack stack = imp.getStack();
    		this.map =new float[stack.getSize()][];
    		for (int i=0;i<this.map.length;i++) this.map[i]= (float[]) stack.getPixels(i+1);
    		this.mapWidth=imp.getWidth();
    		this.mapHeight=imp.getHeight();
    		if (imp.getProperty("path")!=null) this.path=(String) imp.getProperty("path");
    		this.imagePlane[0]= Double.parseDouble((String) imp.getProperty("imagePlaneX"));
    		this.imagePlane[1]= Double.parseDouble((String) imp.getProperty("imagePlaneY"));
    		this.imagePlane[2]= Double.parseDouble((String) imp.getProperty("imagePlaneZ"));
    		this.projectionCenterXY[0]=  Double.parseDouble((String) imp.getProperty("projectionCenterX"));
    		this.projectionCenterXY[1]=  Double.parseDouble((String) imp.getProperty("projectionCenterY"));
    		this.projectionPixelSize=    Double.parseDouble((String) imp.getProperty("projectionPixelSize"));

    		for (int chn=0;chn<this.numOverlapChannels;chn++){
        		this.channel[chn]=Integer.parseInt  ((String) imp.getProperty("channel"+(chn+1)));
        		this.azimuth[chn]= Double.parseDouble((String) imp.getProperty("azimuth"+(chn+1)));
        		this.radius[chn]=  Double.parseDouble((String) imp.getProperty("radius"+(chn+1)));
        		this.height[chn]=  Double.parseDouble((String) imp.getProperty("height"+(chn+1)));
    		}
    	}

    	public void argb24ToAYCbCr(
    			ImagePlus imp,
    			int imgNum, // 0 - first image, 1 - second image in teh pair
    			double [][] aYCbCr){ // should be initialized as double [8][]
    		double [][] aYCbCr1=argb24ToAYCbCr(imp);
    		for (int i=0;i<aYCbCr1.length;i++){
    			aYCbCr[i+imgNum*aYCbCr1.length]= aYCbCr1[i];
    		}
    	}

    	public double [][] argb24ToAYCbCr(
    			ImagePlus imp1,
    			ImagePlus imp2
    			){
    		double [][] aYCbCr1=argb24ToAYCbCr(imp1);
    		double [][] aYCbCr2=argb24ToAYCbCr(imp2);
    		double [][] aYCbCr=new double [2*aYCbCr1.length][];
    		for (int i=0;i<aYCbCr1.length;i++){
    			aYCbCr[i]=               aYCbCr1[i];
    			aYCbCr[i+aYCbCr1.length]=aYCbCr2[i];
    		}
    		return aYCbCr;
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


    	public double [] showThresholdedDisparity(
    			double [][][] disparityMap,
    			int zeroShift, // 1/2 corrFFTSize
    			double threshold
    			){
    		double [] result =new double [disparityMap.length*disparityMap[0].length];
    		for (int index=0;index<result.length;index++) {
    			result[index]=Double.NaN;
    			int y=index/disparityMap[0].length;
    			int x=index%disparityMap[0].length;
    			if (disparityMap[y][x]!=null){
    				int iMax=-1;
    				double max=0.0;
    				for (int i=zeroShift;i<disparityMap[y][x].length;i++) if (disparityMap[y][x][i]>max){
    					iMax=i;
    					max=disparityMap[y][x][i];
    				}
    				if (max>=threshold){
    					if ((iMax==zeroShift) || (iMax==(disparityMap[y][x].length-1))){
    						result[index]=iMax-zeroShift;
    					} else {
    						// y=ax^2+bx+c, x==0 for the middle point, c=y[0] already defined
    						double b=0.5*(disparityMap[y][x][iMax+1]-disparityMap[y][x][iMax-1]);
    						double a=0.5*(disparityMap[y][x][iMax+1]+disparityMap[y][x][iMax-1])-disparityMap[y][x][iMax]; // negative
//    						double x=-0.5*b/a;
    						result[index]=iMax-zeroShift+(-0.5*b/a);
    					}
    				}
    			}
    		}
    		return result;
    	}
    	public double [] showThresholdedDisparityPixels(
    			double [][][] disparityMap,
    			int zeroShift, // 1/2 corrFFTSize
    			double threshold
    			){
    		double [] map=showThresholdedDisparity(
        			disparityMap,
        			zeroShift, // 1/2 corrFFTSize
        			threshold);

    		double [] result =new double [this.mapWidth*this.mapHeight];
    		for (int index=0;index<result.length;index++) {
    			int y=index/this.mapWidth;
    			int x=index%this.mapWidth;
    			int ny=y*this.disparityMapOverlapFraction/this.disparityMapFFTSize;
    			int nx=x*this.disparityMapOverlapFraction/this.disparityMapFFTSize;
    			if ((ny<disparityMap.length)&&(nx<disparityMap[0].length)){
    				result[index]=map[ny*disparityMap[0].length+nx];
    			} else result[index]=Double.NaN;
    		}
    		return result;
    	}


    	public double [][][] correlateImage(
    			final boolean side,
    			final boolean autoCorrelation,
    			final int corrFFTSize,
    			final int overlapFraction,
//				int corrYC,
    			final int maxDisparity,
    			final double phaseCorrelationFraction,
    			final double corrCbWeight,
    			final double corrCrWeight,
    			final double correlationHighPassSigma,
    			final double correlationLowPassSigma,
    			final double noiseNormalizationSignaY,
    			final double noiseNormalizationSignaCbCr,
    			final double minFracArea,
    			final boolean useBinaryAlpha,
				int threadsMax,
				final int debugLevel){
    		final int numTileRows= this.mapHeight*overlapFraction/corrFFTSize+
    		(((this.mapHeight*overlapFraction/corrFFTSize)*corrFFTSize/overlapFraction<this.mapHeight)?1:0);
 		   final Thread[] threads = newThreadArray(threadsMax);
		   final AtomicInteger rowNum = new AtomicInteger(0);
		   final AtomicInteger rowsFinished = new AtomicInteger(0);
		   final double [][][] result=new double[numTileRows][][];

	   		for (int ithread = 0; ithread < threads.length; ithread++) {
	   			threads[ithread] = new Thread() {
	   				@Override
					public void run() {
	   					DoubleFHT doubleFHT=new DoubleFHT();
	   					for (int row=rowNum.getAndIncrement(); row<numTileRows;row=rowNum.getAndIncrement()) {
	   						int corrYC=row*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
	   						result[row]=correlateOneRow(
	   								doubleFHT,
	   				    			side,
	   				    			autoCorrelation,
	   								corrFFTSize,
	   								overlapFraction,
	   								corrYC,
	   								maxDisparity,
	   								phaseCorrelationFraction,
	   								corrCbWeight,
	   								corrCrWeight,
	   								correlationHighPassSigma,
	   								correlationLowPassSigma,
	   								noiseNormalizationSignaY,
	   								noiseNormalizationSignaCbCr,
	   								minFracArea,
	   								useBinaryAlpha,
	   								debugLevel);
    						final int numFinished=rowsFinished.getAndIncrement();
    						SwingUtilities.invokeLater(new Runnable() {
    							@Override
								public void run() {
    								IJ.showProgress(numFinished,numTileRows);
    							}
    						});

	   					}
	   				}

	   			};
	   		}
	   		startAndJoin(threads);
	   		IJ.showProgress(1.0);
	   		return result;
    	}


    	/**
    	 * Calculate correlation between the two images assuming strict horizontal disparity
    	 * @param doubleFHT instance of DoubleFHT to re-use, or null
    	 * @param side If true - calculate for right (second) image
    	 * @param autoCorrelation If true, calculate autocorrelation (to detect periodic patterns that lead to ambiguity
    	 * @param corrFFTSize FFT size (power of 2)
    	 * @param overlapFraction Shift by this fraction of the corrFFTSize between tiles
    	 * @param corrYC Center Y for the correlation
    	 * @param maxDisparity  Maximal disparity value (in pixels)
    	 * @param phaseCorrelationFraction 0.0 - normal correlation, 1.0 - phase correlation
    	 * @param corrCbWeight  Relative weight of the Cb component (<=0 - skip)
    	 * @param corrCrWeight  Relative weight of the Cr component (<=0 - skip)
    	 * @param correlationHighPassSigma High pass filter for correlation
    	 * @param correlationLowPassSigma Low-pass filter for correlation
    	 * @param noiseNormalizationSigmaY  Normalize Y component to the RMS value (do nothing if <=0.0)
    	 * @param noiseNormalizationSigmaCbCr Normalize color components to the RMS value (do nothing if <=0.0)
    	 * @param minFracArea use only tiles with the defined area (fraction of the full corrFFTSize*corrFFTSize) exceeds this value
    	 * @param useBinaryAlpha Accept any alpha>0.0 when calculating areas
    	 * @param debugLevel
    	 * @return For each tile (with the step corrFFTSize/overlapFraction) where defined area exceeds minimal limits return horizontal profile. For undefined tiles retuern null
    	 */


    	public double [][] correlateOneRow(
    			DoubleFHT doubleFHT,
    			boolean side,
    			boolean autoCorrelation,
				int corrFFTSize,
				int overlapFraction,
				int corrYC,
				int maxDisparity,
				double phaseCorrelationFraction,
				double corrCbWeight,
				double corrCrWeight,
				double correlationHighPassSigma,
				double correlationLowPassSigma,
				double noiseNormalizationSignaY,
				double noiseNormalizationSignaCbCr,
				double minFracArea,
				boolean useBinaryAlpha,
				int debugLevel){
    		int numLayers=this.overlapImages.length/2;
//    		int indexAlpha=0;
//    		int indexY= 1;
    		int indexCb=2;
    		int indexCr=3;

    		int numTiles=this.mapWidth*overlapFraction/corrFFTSize;
    		if (doubleFHT==null) doubleFHT=new DoubleFHT();
    		double[] hamming1d=doubleFHT.getHamming1d(corrFFTSize,0.0); // pure shifted cosine, not real Hamming
			double [] componentWeights={
					1.0/(1.0+corrCbWeight+corrCrWeight),
					corrCbWeight/(1.0+corrCbWeight+corrCrWeight),
					corrCrWeight/(1.0+corrCbWeight+corrCrWeight)};
    		double [][][] selection = new double[numTiles] [][];
    		double [] corr=new double[corrFFTSize*corrFFTSize];
    		double [] frequencyFilter=null;
    		for (int nTile=0;nTile<numTiles;nTile++){
    			int corrXC=corrFFTSize/2/overlapFraction+ corrFFTSize*nTile/overlapFraction;
    			selection[nTile] = getSelection(
        				autoCorrelation?(side?2:1):0,
        						corrFFTSize,
        						corrFFTSize,
        						corrXC,
        						corrYC,
        						0);// positive - shift second image, negative - shift first image
    			double [] stats= selectionStats(selection[nTile],useBinaryAlpha);
    			if (stats[0]<minFracArea) selection[nTile][0]=null;
    			if (stats[1]<minFracArea) selection[nTile][numLayers]=null;
    			if ((selection[nTile][0]!=null) || (selection[nTile][numLayers]!=null)) {
    				if (corrCbWeight<=0) {
    					selection[nTile][indexCb]=null;
    					selection[nTile][indexCb+numLayers]=null;
    				}
    				if (corrCrWeight<=0) {
    					selection[nTile][indexCr]=null;
    					selection[nTile][indexCr+numLayers]=null;
    				}
    				double [][] window={
    						(selection[nTile][0]!=null)?selection[nTile][0].clone():null,
    								(selection[nTile][numLayers]!=null)?selection[nTile][numLayers].clone():null};
    				for (int iy=0;iy<corrFFTSize;iy++) for (int ix=0;ix<corrFFTSize;ix++){
    					int index=iy*corrFFTSize+ix;
    					double h=hamming1d[iy]*hamming1d[ix];
    					if (window[0]!=null) window[0][index]*=h;
    					if (window[1]!=null) window[1][index]*=h;
    				}
    				// window and forward FHT
    				for (int nImage=0;nImage<2;nImage++) if (selection[nTile][nImage*numLayers]!=null) {
    					for (int nComp=1;nComp<=3;nComp++) if (selection[nTile][nComp+nImage*numLayers]!=null){
    		    			normalizeAndWindow(selection[nTile][nComp+nImage*numLayers], window[nImage],true);
    		    			doubleFHT.swapQuadrants(selection[nTile][nComp+nImage*numLayers]);
    		    			doubleFHT.transform(selection[nTile][nComp+nImage*numLayers],false);
    		    			if ((frequencyFilter==null) && ((correlationHighPassSigma>0.0) || (correlationLowPassSigma>0.0))){
    		    				frequencyFilter=doubleFHT.createFrequencyFilter(
    		    	    				selection[nTile][nComp+nImage*numLayers], // just for length
    		    	    				correlationHighPassSigma,
    		    	    				correlationLowPassSigma);
    		    			}
    					}
    				}
    			}
    		}
    		double [][] result = new double [numTiles][];
    		for (int nTile=0;nTile<numTiles;nTile++) result[nTile]=null;
    		double [] first;
    		double [] second;
    		int span= 1+maxDisparity*overlapFraction/corrFFTSize;
    		int outLength=corrFFTSize+(span-1)*corrFFTSize/overlapFraction;
    		double [][]debugCorr=null;
			if (debugLevel>2) debugCorr=new double[span][];
			for (int nTile=0;nTile<numTiles;nTile++) {
				if (debugLevel>2) for (int i=0;i<span;i++){
					debugCorr[i]=null;
				}
				for (int nDisp=0;nDisp<span;nDisp++){
					int numFirst= nTile+(side?nDisp:0);
					int numSecond=nTile-(side?0:nDisp);
					if (debugLevel>2) for (int i=0;i<corr.length;i++) corr[i]=0; // will combine color components
					if ((numFirst<numTiles) && (numSecond>=0) && (selection[numFirst][0]!=null) && (selection[numFirst][numLayers]!=null)){
						for (int i=0;i<corr.length;i++) corr[i]=0; // will combine color components
						if (result[nTile]==null){
							result[nTile]=new double[outLength];
							for (int i=0;i<outLength;i++) result[nTile][i]=0.0;
						}
						// Calculate and combine Y,Cb,Cr components
						for (int nComp=1;nComp<numLayers;nComp++){
							if ((selection[numFirst][nComp]!=null) && (selection[numFirst][nComp+numLayers]!=null)){
								first= selection[numFirst][nComp].clone(); // frequency domain
								second=selection[numSecond][nComp+numLayers]; // no need to clone here
								first= doubleFHT.phaseMultiply(first, second, phaseCorrelationFraction);
								if (frequencyFilter!=null)  doubleFHT.multiplyByReal(first, frequencyFilter);
								doubleFHT.transform(first,true) ; // inverse transform
								doubleFHT.swapQuadrants(first);
								//        					corr[nComp]=first;

								double sigma=(nComp==1)?noiseNormalizationSignaY:noiseNormalizationSignaCbCr;
								if (sigma>0){
									double rms=hFNoiseNormalize(
											first, // double [] data, will be modified
											sigma, // double filterSigma,
											true);  // boolean centerOnly);
									if (debugLevel>3){
										System.out.println("Correlation RMS for component "+((nComp==1)?"Y":((nComp==2)?"Cb":"Cr"))+ " was "+rms);
									}
								}
								for (int i=0;i<corr.length;i++) corr[i]+=componentWeights[nComp-1]*first[i]; // not all pixels are needed, just the centerline
							}
						}
						// debug-show partial correlation?
						// Accumulate result
						double scale=2.0/overlapFraction;
						int dstIndex=nDisp*corrFFTSize/overlapFraction;
						for (int i=0;i<corrFFTSize;i++){
							if ((dstIndex>=0) && (dstIndex<outLength)) result[nTile][dstIndex]+=scale*corr[corrFFTSize*corrFFTSize/2+i];
							dstIndex++;
						}
					}
					if (debugLevel>2) {
						debugCorr[nDisp]=corr.clone();
						System.out.println("correlateOneRow(): nTile="+nTile+" nDisp="+nDisp+" numFirst="+numFirst+" numSecond="+numSecond);
					}
				}
				if (debugLevel>2) {
        			(new ShowDoubleFloatArrays()).showArrays(
        					debugCorr,
        					corrFFTSize,
        					corrFFTSize,
        					true,
        			"corr"+nTile+"-y"+corrYC+"-MAXDSP"+maxDisparity+"-PC"+phaseCorrelationFraction);

				}
			}
    		return result;
    	}

    	/**
    	 * Merge several representations of the linear features (pointed from different cells)
    	 * @param corrFFTSize size of the square tile side
    	 * @param overlapFraction Shift by this fraction of the corrFFTSize between tiles
    	 * @param features array of [row][column][parameters] describing the linear feature (see tileLineDetect() method)
    	 * @param strengthMode scale normalized linear feature to (0.0): absolute strength, 1.0 - relative strength (0.0<strengthMode<1.0): intermediate, (-1): phaseStrength, -2 - no scaling
    	 * @param directionTolerance Maximal angle difference between merging features (radians)
    	 * @param normalDistanceTolerance Maximal distance between merging features (pixels) - orthogonal to the linear feature
    	 * @param tangentialDistanceTolerance Maximal distance between merging features (pixels) - parallel to the linear feature
    	 * @param hostsTolerance Merge neighbors claiming the same (close) point as belonging to each
    	 * @param directionSigma Influence of the merging features decreases as the Gaussian of angle with this sigma
    	 * @param normalDistanceSigma Influence of the merging features decreases as the Gaussian of distance between them with this sigma (orthogonal to the linear feature)
    	 * @param tangentialDistanceSigma Same in the direction of the linear feature
    	 * @param minMerged do not output features that have less number of merged cells (0 - non-merged feature outside of current cell, 1 nothing merged, feature inside cell, >1 - several merged)
    	 * @param scaleDistances // ~1.05? compensate for decrease in measured distances to the features (partially caused by windowing)
    	 * @param debugLevel Debug level
    	 * @return array of the same format as input "featured" with some features merged and some - removed
    	 */
     	public double [][][] mergeLinearFeatures(
    			int corrFFTSize,
    			int overlapFraction, // default 8
    			double [][][] features,
    			double strengthMode, // 0.0 - absolute, 1.0 - relative, 0.5 - "balanced", "-1" - use phaseStrength instead
    			double directionTolerance,
    			double normalDistanceTolerance,
    			double tangentialDistanceTolerance,
     			double hostsTolerance,
    			double directionSigma, // make a fixed fraction of directionTolerance?
    			double normalDistanceSigma, // make a fixed fraction of distanceTolerance?
    			double tangentialDistanceSigma, // make a fixed fraction of distanceTolerance?
    			int minMerged,
    			double cellUsageShift,
    			double scaleDistances,
    			double swapTangentialTolerance,
    			int    swapSearchRange,
        		int debugRow,
        		int debugColumn,
				int debugLevel){
//    		int debugRow=144; //106; //153; //148;
//    		int debugColumn=43; //60;// 49; //53;
//    		int modDebugLevel=debugLevel;

    		int height=features.length;
    		int width=features[0].length;

    		// now include not just same cell, but all other hosts
    		int [][] cellUsage=   getCellUsage(
        			corrFFTSize,
        			overlapFraction, // default 8
        			features,
        			directionTolerance,
        			normalDistanceTolerance,
        			tangentialDistanceTolerance,
         			hostsTolerance,
        			scaleDistances,
            		debugRow,
            		debugColumn,
    				debugLevel);
    		// find all cells where feature is in the same cell
    		if (debugLevel>1){
    			System.out.println("mergeLinearFeatures()1A:"+
    					" directionTolerance="+          directionTolerance+
    					" directionSigma="+              directionSigma+
  //  					" kAngle="+            kAngle+
    					" normalDistanceTolerance="+     normalDistanceTolerance+
    					" tangentialDistanceTolerance="+ tangentialDistanceTolerance+
    					" normalDistanceSigma="+         normalDistanceSigma+
    					" tangentialDistanceSigma="+     tangentialDistanceSigma
//    					+" kDist="+             kDist
    					);
    		}
    		double [][][] shares=new double [height][width][];
    		int  [][][][] hosts= new int    [height][width][][];
    		for (int row=0;row<height;row++) for (int col=0;col<width;col++){
    			shares[row][col]=null;
    			hosts [row][col]=null;
    		}
    		getFeatureShares(
         			hosts, // will be updated
         			shares,// will be updated
         			cellUsage,
        			corrFFTSize,
        			overlapFraction, // default 8
        			features,
        			0.0, //strengthMode, // 0.0 - absolute, 1.0 - relative, 0.5 - "balanced", "-1" - use phaseStrength instead
        			directionTolerance,
        			normalDistanceTolerance,
        			tangentialDistanceTolerance,
         			hostsTolerance,
        			directionSigma, // make a fixed fraction of directionTolerance?
        			normalDistanceSigma, // make a fixed fraction of distanceTolerance?
        			tangentialDistanceSigma, // make a fixed fraction of distanceTolerance?
        			scaleDistances,
            		debugRow,
            		debugColumn,
    				debugLevel);
    		// combine features
    		double [][][] filteredFeatures=combineFeatures(
         			hosts,
         			shares,
         			cellUsage, // will be updated to number of cells combined
        			corrFFTSize,
        			overlapFraction, // default 8
        			features,
        			0.0, // strengthMode, // 0.0 - absolute, 1.0 - relative, 0.5 - "balanced", "-1" - use phaseStrength instead
        			directionTolerance,
        			normalDistanceTolerance,
        			tangentialDistanceTolerance,
         			hostsTolerance,
        			directionSigma, // make a fixed fraction of directionTolerance?
        			normalDistanceSigma, // make a fixed fraction of distanceTolerance?
        			tangentialDistanceSigma, // make a fixed fraction of distanceTolerance?
        			scaleDistances,
            		debugRow,
            		debugColumn,
    				debugLevel);

    		for (int row=0;row<height;row++){
    			for (int col=0;col<width;col++){
    				if (cellUsage[row][col]<minMerged) filteredFeatures[row][col]=null;
    			}
    		}
    		if ((swapSearchRange>=0) && (swapTangentialTolerance>=0)){
    			for (int n=1;n<corrFFTSize;n++){
    			   int [] swapResults=swapHosts( // return num. swapped, max distance?
    		     			cellUsage,
    		    			corrFFTSize,
    		    			overlapFraction, // default 8
    		    			filteredFeatures,
    		    			0.0, //strengthMode, // 0.0 - absolute, 1.0 - relative, 0.5 - "balanced", "-1" - use phaseStrength instead
    		    			swapTangentialTolerance,
    		    			swapSearchRange,
    		        		debugRow,
    		        		debugColumn,
    						debugLevel);
    			   if (swapResults[0]==0) break;
    			}

    		}

    		if (!Double.isNaN(cellUsageShift)){
        		for (int row=0;row<height;row++){
        			for (int col=0;col<width;col++) if (filteredFeatures[row][col]!=null){
        				double k=(cellUsage[row][col]-cellUsageShift);
        				if (k<0) k=0;
        				filteredFeatures[row][col][0]*=k;
        			}
        		}
    		}
    		return filteredFeatures;
    	}

    	private int [] swapHosts( // return num. swapped, max distance?
     			int [][] cellUsage,
    			int corrFFTSize,
    			int overlapFraction, // default 8
    			double [][][] features,
    			double strengthMode, // 0.0 - absolute, 1.0 - relative, 0.5 - "balanced", "-1" - use phaseStrength instead
    			double tangentialDistanceTolerance,
    			int scanRange,
//    			double scaleDistances, // should be already applied!
        		int debugRow,
        		int debugColumn,
				int debugLevel){
    		int modDebugLevel=debugLevel;
    		int height=features.length;
    		int width=features[0].length;
    		int gridStep=corrFFTSize/overlapFraction;
    		int debugThreshold=2;
    		// find all cells where feature is in closer to the cell center, than matching feature - to the center of any other cell
    		int [][][] hosts=new int [height][width][];
    		double [][] hostStrengths=new double[height][width];
    		for (int row=0;row<height;row++){
        		for (int col=0;col<width;col++){
        			hosts[row][col]=null;
        			hostStrengths[row][col]=0.0;
        		}
    		}
    		int numOutsideCell=0;
    		int numWantMove=0;
    		double maximalDistance=0.0;
    		// find which cells want to move
    		for (int row=0;row<height;row++){
        		for (int col=0;col<width;col++) if (features[row][col]!=null){
//        			modDebugLevel=debugLevel+(((row==debugRow) && (col==debugColumn))?2:0);
        			modDebugLevel=debugLevel+(((Math.abs(row-debugRow)<=2) && (Math.abs(col-debugColumn)<=2))?2:0);
        			double distance=features[row][col][4];
        			double angle=   features[row][col][3];
        			double dX=distance*Math.cos(angle);
        			double dY=distance*Math.sin(angle);
        			int row0=row+ ((int) Math.round(dY/gridStep));
        			int col0=col+ ((int) Math.round(dX/gridStep));
        			double bestNewDist= distance;
        			int [] bestColRow={-1,-1};
        			if ((row0!=row) || (col0!=col)) {
    					if (modDebugLevel>debugThreshold){
							System.out.println("swapHosts()1: row="+row+" col="+col+
									" row0="+row0+" col0="+col0+
									" dX="+dX+" dY="+dY+
									" angle="+IJ.d2s(180/Math.PI*angle,2)+
									" distance="+distance);
    					}
        				numOutsideCell++;
        				if (distance<0.0){ // maybe it is not needed here at all?
        					distance=-distance;
        					angle+=Math.PI;
        					angle-=2*Math.PI*Math.round(angle/(2*Math.PI));
        				}
        				if ((distance>maximalDistance) && // just for debug/statistics, remove border elements
        						(row> overlapFraction/2) &&
        						(row< (height-overlapFraction/2)) &&
        						(col> overlapFraction/2) &&
        						(col< (width-overlapFraction/2))) {
        					maximalDistance=distance;
        					if (modDebugLevel>debugThreshold){
        						System.out.println("swapHosts()2: row="+row+" col="+col+	" distance="+distance+" maximalDistance"+maximalDistance);
        					}
        				}
        				double sinA=Math.sin(angle);
        				double cosA=Math.cos(angle);
        				for (int dRow=-scanRange;dRow<=scanRange;dRow++) {
    						int row1=row0+dRow;
    						if ((row1>=0) && (row1<height)) {
    							for (int dCol=-scanRange;dCol<=scanRange;dCol++){
    								int col1=col0+dCol;
    								if ((col1>=0) && (col1<width)) {
    									double xc=(col1-col)*gridStep;
    									double yc=(row1-row)*gridStep;
    									double newDist=xc*cosA+yc*sinA - distance;
    									if (features[row1][col1]== null){ // only considering empty cells
    										if ((Math.abs(newDist)<bestNewDist) && (Math.abs(yc*cosA-xc*sinA)<=tangentialDistanceTolerance)){ // improvement
    											bestNewDist=Math.abs(newDist);
    											bestColRow[0]=col1;
    											bestColRow[1]=row1;
    					    					numWantMove++;
        										if (modDebugLevel>debugThreshold){
        											System.out.println("swapHosts()3: row="+row+" col="+col+
        													" row0="+row0+" col0="+col0+
        													" row1="+row1+" col1="+col1+
        													" cell is "+((features[row1][col1]== null)?"FREE":"USED")+
        													" xc="+xc+" yc="+yc+
        													" xc*cosA+yc*sinA="+(xc*cosA+yc*sinA)+
        													" row1="+row1+" col1="+col1+
        													" angle="+IJ.d2s(180/Math.PI*angle,2)+
        													" distance="+distance+
        													" newDist="+newDist);
        										}
    										}
    									}
    								}
    							}
    						}
    					}
    					if (bestNewDist<distance) {
    						if (hosts[bestColRow[1]][bestColRow[0]]==null){ //-1
    							hosts[bestColRow[1]][bestColRow[0]]=new int[2];
    						}
    						double thisStrength=combinedStrength(strengthMode,features[row][col]);

    						if (thisStrength>hostStrengths[bestColRow[1]][bestColRow[0]]){
    							hostStrengths[bestColRow[1]][bestColRow[0]]=thisStrength;
    							hosts[bestColRow[1]][bestColRow[0]][0]=col;
    							hosts[bestColRow[1]][bestColRow[0]][1]=row;
    						}
				    		if (modDebugLevel>debugThreshold){
				    			System.out.print  ("swapHosts()4: bestColRow[1]="+bestColRow[1]+" bestColRow[0]="+bestColRow[0]);
				    			System.out.print  (
				    					" hosts["+bestColRow[1]+"]["+bestColRow[0]+"][0]="+hosts[bestColRow[1]][bestColRow[0]][0]+
				    					" hosts["+bestColRow[1]+"]["+bestColRow[0]+"][1]="+hosts[bestColRow[1]][bestColRow[0]][1]);
						    			System.out.println(" thisStrength="+thisStrength+" hostStrengths[bestColRow[1]][bestColRow[0]]="+hostStrengths[bestColRow[1]][bestColRow[0]]);
				    		}
    					}
        			}
        		}
    		}
    		if (debugLevel>2) System.out.println(">>>> swapHosts(): numWantMove="+numWantMove+" numOutsideCell="+numOutsideCell+" maximaDistance="+maximalDistance);
    		int [] result={numWantMove,numWantMove,numOutsideCell};
    		if (numWantMove==0) return result;
    		// now move linear features to closer to them cells
    		int numMoved=0;
    		for (int row=0;row<height;row++){
        		for (int col=0;col<width;col++) if (hosts[row][col]!=null){
        			modDebugLevel=debugLevel+(((Math.abs(row-debugRow)<=2) && (Math.abs(col-debugColumn)<=2))?2:0);
        			int row1=hosts[row][col][1];
        			int col1=hosts[row][col][0];
        			double [] feature=features[row1][col1]; // no need to clone - original will be null

        			features[row1][col1]=null;
        			double distance=feature[4];
        			double angle=   feature[3];
					if (distance<0.0){ // maybe it is not needed here at all?
						distance=-distance;
						angle+=Math.PI;
						angle-=2*Math.PI*Math.round(angle/(2*Math.PI));
						feature[5]=-feature[5];
					}
					double sinA=Math.sin(angle);
					double cosA=Math.cos(angle);
					double xc=(col-col1)*gridStep;
					double yc=(row-row1)*gridStep;
		    		if (modDebugLevel>debugThreshold){
						int corrYC= row*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
						int corrXC= col*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
						int corrYC1=row1*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
						int corrXC1=col1*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
		    			System.out.print("swapHosts()5: "+" YC="+corrYC+" XC="+corrXC+
		    					" row="+row+" col="+col+
							" row1="+row1+" col1="+col1+" YC1="+corrYC1+" XC1="+corrXC1+
							" angle="+IJ.d2s(180/Math.PI*angle,2)+
							" distance="+distance);
		    		}
//					distance=xc*cosA+yc*sinA - distance;
					distance-=xc*cosA+yc*sinA;
					if (distance<0.0){ // maybe it is not needed here at all?
						distance=-distance;
						angle+=Math.PI;
						angle-=2*Math.PI*Math.round(angle/(2*Math.PI));
						feature[5]=-feature[5];
					}
        			feature[3]=angle;
        			feature[4]=distance;
		    		if (modDebugLevel>debugThreshold){
		    			System.out.println(	" new angle="+IJ.d2s(180/Math.PI*angle,2)+
							" new distance="+distance);
		    		}
        			features[row][col]=feature;
        			cellUsage[row][col]=cellUsage[row1][col1];
        			cellUsage[row1][col1]=-1;
        			numMoved++;
        		}
    		}
    		result[0]=numMoved;
    		if (debugLevel>1) System.out.println(">>>> swapHosts(): numMoved="+numMoved+" numWantMove="+numWantMove+" numOutsideCell="+numOutsideCell+" maximaDistance="+maximalDistance);
    		return result;
     	}





     	private double getFeatureMatch(
     			int row10,
     			int col10,
     			int row20,
     			int col20,
     			int gridStep,
    			double [][][] features,
    			double directionTolerance,
    			double normalDistanceTolerance,
    			double tangentialDistanceTolerance,
    			double directionSigma, // make a fixed fraction of directionTolerance?
    			double normalDistanceSigma, // make a fixed fraction of distanceTolerance?
    			double tangentialDistanceSigma, // make a fixed fraction of distanceTolerance?
    			double scaleDistances,
				int debugLevel){
     		// swap ends to make sure the result absolutely does not depend on the order in a pair
     		int row1,row2,col1,col2;
     		if ((row20>row10) || ((row20==row10) && (col20> col10))){
     			row1=row10;
     			col1=col10;
     			row2=row20;
     			col2=col20;
     		} else {
     			row1=row20;
     			col1=col20;
     			row2=row10;
     			col2=col10;
     		}
     		int height=features.length;
     		int width= features[0].length;
     		if (	(row1<0) || (col1<0) || (row1>=height) || (col1>=width) ||
     				(row2<0) || (col2<0) || (row2>=height) || (col2>=width) ||
     				(features[row1][col1]==null) || (features[row2][col2]==null)) return 0;
     		double angle1= features[row1][col1][3];
     		angle1-=(2*Math.PI)*Math.round(angle1/(2*Math.PI));
     		double angle2= features[row2][col2][3];
     		angle2-=(2*Math.PI)*Math.round(angle2/(2*Math.PI));
			double dAngle=angle2-angle1;
			dAngle-=Math.PI*Math.round(dAngle/Math.PI);
			if (Math.abs(dAngle)>directionTolerance) return 0;
			if (Math.abs(angle1-angle1)>(Math.PI/2)){
				angle2+=Math.PI;
				angle2-=(2*Math.PI)*Math.round(angle2/(2*Math.PI));
			}
			double angle=(angle1+angle2)/2; // here ignoring weights
			double sinAngle=Math.sin(angle);
			double cosAngle=Math.cos(angle);


     		double distance1= scaleDistances*features[row1][col1][4];
     		double distance2= scaleDistances*features[row2][col2][4];
     		double x1=distance1*Math.cos(angle1);
     		double y1=distance1*Math.sin(angle1);
     		double x2=distance2*Math.cos(angle2)+gridStep*(col2-col1);
     		double y2=distance2*Math.sin(angle2)+gridStep*(row2-row1);

     		double dX=x2-x1;
     		double dY=y2-y1;
     		double dNormal=-sinAngle*dX+cosAngle*dY;
     		if (Math.abs(dNormal)>normalDistanceTolerance) return 0;
     		double dTangential=cosAngle*dX+sinAngle*dY;
     		if (Math.abs(dTangential)>tangentialDistanceTolerance) return 0;
     		if ((directionSigma<=0) || (normalDistanceSigma<=0) || (tangentialDistanceSigma<=0)) return 1.0;
     		return Math.exp(-(
     				dAngle*dAngle/(directionSigma*directionSigma)+
     				dNormal*dNormal/(normalDistanceSigma*normalDistanceSigma)+
     				dTangential*dTangential/(tangentialDistanceSigma*tangentialDistanceSigma)
     				)/2.0);
     	}

     	private double [][][] combineFeatures(
     			int [][][][] hosts,
     			double [][][] shares,
     			int [][] cellUsage,
    			int corrFFTSize,
    			int overlapFraction, // default 8
    			double [][][] features,
    			double strengthMode, // 0.0 - absolute, 1.0 - relative, 0.5 - "balanced", "-1" - use phaseStrength instead
    			double directionTolerance,
    			double normalDistanceTolerance,
    			double tangentialDistanceTolerance,
     			double hostsTolerance,
    			double directionSigma, // make a fixed fraction of directionTolerance?
    			double normalDistanceSigma, // make a fixed fraction of distanceTolerance?
    			double tangentialDistanceSigma, // make a fixed fraction of distanceTolerance?
    			double scaleDistances,
        		int debugRow,
        		int debugColumn,
				int debugLevel){
    		int modDebugLevel=debugLevel;
    		int height=features.length;
    		int width=features[0].length;
    		int scanAroundCells=overlapFraction/2;
    		int gridStep=corrFFTSize/overlapFraction;
       		double [][][] filteredFeatures = new double [height][width][];
    		double interCenter=corrFFTSize/overlapFraction;
    		// find all cells where feature is in closer to the cell center, than matching feature - to the center of any other cell
//    		int [][] theseHosts=  new int [(2*scanAroundCells+1)*(2*scanAroundCells+1)][2];
//    		double[] theseShares=new double [(2*scanAroundCells+1)*(2*scanAroundCells+1)];
    		for (int row=0;row<height;row++){
        		for (int col=0;col<width;col++) if (cellUsage[row][col]>0){ // "hosts" only
//        			modDebugLevel=debugLevel+(((row==debugRow) && (col==debugColumn))?2:0);
        			modDebugLevel=debugLevel+(((Math.abs(row-debugRow)<=1) && (Math.abs(col-debugColumn)<=1))?2:0);

        			double distance=scaleDistances*features[row][col][4];
        			double angle=   features[row][col][3];
					if (distance<0.0){ // maybe it is not needed here at all?
						distance=-distance;
						angle+=Math.PI;
						angle-=2*Math.PI*Math.round(angle/(2*Math.PI));
//						phaseAtZero=-phaseAtZero; // not used here
					}



        			double dX=distance*Math.cos(angle);
        			double dY=distance*Math.sin(angle);
        			int row0=row+ ((int) Math.round(dY/gridStep));
        			int col0=col+ ((int) Math.round(dX/gridStep));

					if (modDebugLevel>2){
						int corrYC=row*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
						int corrXC=col*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;

						System.out.println("combineFeatures()0:, row="+row+" col="+col+" YC="+corrYC+" XC="+corrXC+
								" row0="+row0+" col0="+col0+
								" angle="+IJ.d2s(180/Math.PI*angle,2)+
								" distance="+IJ.d2s(distance,2)+
								" phase zero="+IJ.d2s(180/Math.PI*features[row][col][5],2)
						);
					}

					double sumWeights=0.0, sumWDist=0.0, sumX=0.0, sumY=0.0, sumSinAngle=0.0, sumCosAngle=0.0, sumSinPhase=0.0, sumCosPhase=0.0,
				       sumAbsStrength=0.0, sumReferenceLevel=0.0, sumPhaseStrength=0.0;
					for (int dRow=-scanAroundCells;dRow<=scanAroundCells;dRow++) {
						int row1=row0+dRow;
						if ((row1>=0) && (row1<height)) {
							for (int dCol=-scanAroundCells;dCol<=scanAroundCells;dCol++){
								int col1=col0+dCol;
								if ((col1>=0) &&(col1<width) && (features[row1][col1]!=null) ){ //
    								// find share (if exists)
									double share=0.0;
									if (hosts[row1][col1]!=null) {
										int neibIndex=-1;
										for (int i=0;i<hosts[row1][col1].length;i++) if ((hosts[row1][col1][i][0]==col) && (hosts[row1][col1][i][1]==row)){
											neibIndex=i;
											break;
										}
										if (neibIndex>=0) share=shares[row1][col1][neibIndex]; // so the "energy" will be split between the hosts, preserving totals
									} else if ((row==row1) && (col==col1)){
										share=1.0; // self
									}
									if (share>0.0){
										double wDistance=getFeatureMatch(
												row,
												col,
												row1,
												col1,
												gridStep,
												features,
												directionTolerance,
												normalDistanceTolerance,
												tangentialDistanceTolerance,
												directionSigma, // make a fixed fraction of directionTolerance?
												normalDistanceSigma, // make a fixed fraction of distanceTolerance?
												tangentialDistanceSigma, // make a fixed fraction of distanceTolerance?
												scaleDistances,
												debugLevel);
										if (wDistance>0){
											double [] feature=features[row1][col1].clone();
											// compensate for negative distance
											if (feature[4]<0.0){
												feature[4]=-feature[4];
												feature[3]+=Math.PI;
												feature[3]-=2*Math.PI*Math.round(feature[3]/(2*Math.PI));
												feature[5]=-feature[5]; // not used here
											}

											double dX1=scaleDistances*feature[4]*Math.cos(feature[3])+(interCenter* (col1-col));//  dCol);
											double dY1=scaleDistances*feature[4]*Math.sin(feature[3])+(interCenter* (row1-row));// dRow);

											double dAngle=feature[3]-angle;
											dAngle-=2*Math.PI*Math.round(dAngle/(2*Math.PI));
											if (Math.abs(dAngle)>(Math.PI/2)){
												feature[3]+=Math.PI;
												feature[3]-=2*Math.PI*Math.round(feature[3]/(2*Math.PI));
												feature[5]=-feature[5]; //phaseAtZero
											}
											double wSharedDistance=wDistance * share;
											double w=combinedStrength(strengthMode,features[row][col])*wSharedDistance;
    	    								if (modDebugLevel>2){
    	    									int corrYC=row1*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
    	    									int corrXC=col1*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;

    	    									System.out.println("combineFeatures()1:, row1="+row1+" col1="+col1+" YC="+corrYC+" XC="+corrXC+
    	    											" strength="+features[row1][col1][0]+
    	    											" refLevel="+features[row1][col1][1]+
    	    											" angle="+IJ.d2s(180/Math.PI*features[row1][col1][3],2)+
    	    											" distance="+IJ.d2s(scaleDistances*features[row1][col1][4],2)+
    	    											" phase zero="+IJ.d2s(180/Math.PI*features[row1][col1][5],2)
    	    									);
    	    									System.out.println("combineFeatures()2:, row1="+row1+" col1="+col1+" YC="+corrYC+" XC="+corrXC+
    	    											" strength="+feature[0]+
    	    											" angle="+IJ.d2s(180/Math.PI*feature[3],2)+
    	    											" distance="+IJ.d2s(scaleDistances*feature[4],2)+
    	    											" phase zero="+IJ.d2s(180/Math.PI*feature[5],2)
    	    									);
    	    									System.out.println("combineFeatures()3:, w="+w+
    	    											" wDistance="+wDistance+
    	    											" wSharedDistance="+wSharedDistance+
    	    											" dX1="+dX1+
    	    											" dY1="+dY1+
    	    											" dX="+ dX+
    	    											" dY="+ dY+
    	    											" Math.sin(feature[3])="+Math.sin(feature[3])+
    	    											" Math.cos(feature[3])="+Math.cos(feature[3])+
    	    											" Math.sin(feature[5])="+Math.sin(feature[5])+
    	    											" Math.cos(feature[5])="+Math.cos(feature[5])
    	    									);
    	    								}

											sumWeights+=w;
											sumX+=w*dX1;
											sumY+=w*dY1;
											sumSinAngle+=w*Math.sin(feature[3]); // direction
											sumCosAngle+=w*Math.cos(feature[3]); // direction
											sumSinPhase+=w*Math.sin(feature[5]); // phaseAtZero
											sumCosPhase+=w*Math.cos(feature[5]); // phaseAtZero
											sumWDist+=          wSharedDistance;
											sumAbsStrength+=    wSharedDistance*feature[0];
											sumReferenceLevel+= wSharedDistance*feature[1];
											sumPhaseStrength+=  wSharedDistance*feature[2];
											cellUsage[row][col]++;
    	    								if (modDebugLevel>2){
    	    									System.out.println("combineFeatures()4:"+
    	    											" sumSinAngle="+sumSinAngle+
    	    											" sumCosAngle="+sumCosAngle+
    	    											" sumSinPhase="+sumSinPhase+
    	    											" sumCosPhase="+sumCosPhase+
    	    											" cellUsage["+row+"]["+col+"]="+cellUsage[row][col]);
    	    								}

										} else {
	    									int corrYC=row1*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
	    									int corrXC=col1*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
											System.out.println("*** BUG: row1="+row1+" col1="+col1+" YC="+corrYC+" XC="+corrXC);
											 //bug
										}
									}

								}
							}
						}
					}
					cellUsage[row][col]--; // host was counted twice
					angle=      Math.atan2(sumSinAngle,sumCosAngle);
					double phaseAtZero=Math.atan2(sumSinPhase,sumCosPhase);
					double x=sumX/sumWeights;
					double y=sumY/sumWeights;
//					double distance=x*Math.sin(angle)-y*Math.cos(angle);
					distance=x*Math.cos(angle)+y*Math.sin(angle);

					if (modDebugLevel>2){
						System.out.println("combineFeatures()5: "+
								" angle="+IJ.d2s(180/Math.PI*angle,2)+
								" distance="+distance+
								" phaseAtZero="+phaseAtZero+
								" x="+x+
								" y="+y);

					}

					if (distance<0.0){
						distance=-distance;
						angle+=Math.PI;
						angle-=2*Math.PI*Math.round(angle/(2*Math.PI));
						phaseAtZero=-phaseAtZero;
					}
					double [] feature={
							sumAbsStrength/sumWDist,
							sumReferenceLevel/sumWDist,
							sumPhaseStrength/sumWDist,
							angle,
							distance, // here - already scaled
							phaseAtZero,
							cellUsage[row][col]
							};
					if (modDebugLevel>2){
						System.out.println("combineFeatures()6: "+
								" angle="+IJ.d2s(180/Math.PI*angle,2)+
								" distance="+distance+
								" phaseAtZero="+phaseAtZero+
								" x="+x+
								" y="+y);
						System.out.println("combineFeatures()7:"+
								" strength="+feature[0]+
								" refLevel="+feature[1]+
								" angle="+IJ.d2s(180/Math.PI*feature[3],2)+
								" distance="+IJ.d2s(feature[4],2)+
								" phase zero="+IJ.d2s(180/Math.PI*feature[5],2)
						);
					}

					filteredFeatures[row][col]=feature;
        		}
    		}
    		return filteredFeatures;
     	}
     	public double [] displayUsage(
     			double [][][] features,
    			int corrFFTSize,
    			int overlapFraction, // default 8
    			double scale
     			){
			double [] usage=new double[this.mapWidth*this.mapHeight];
			int k=corrFFTSize/overlapFraction;
			for (int index=0;index<usage.length;index++){
				double d=-1.0;
				int x=(index%this.mapWidth)/k;
				int y=(index/this.mapWidth)/k;
				if ((x<features[0].length) && (y<features.length) && (features[y][x]!=null) ){
					if (features[y][x].length>6) d=features[y][x][6];
					else d=1.0;

				}
				usage[index]=d*scale;
			}
			return usage;

     	}


     	private void getFeatureShares(
     			int [][][][] hosts,
     			double [][][] shares,
     			int [][] cellUsage,
    			int corrFFTSize,
    			int overlapFraction, // default 8
    			double [][][] features,
    			double strengthMode, // 0.0 - absolute, 1.0 - relative, 0.5 - "balanced", "-1" - use phaseStrength instead
    			double directionTolerance,
    			double normalDistanceTolerance,
    			double tangentialDistanceTolerance,
     			double hostsTolerance,
    			double directionSigma, // make a fixed fraction of directionTolerance?
    			double normalDistanceSigma, // make a fixed fraction of distanceTolerance?
    			double tangentialDistanceSigma, // make a fixed fraction of distanceTolerance?
    			double scaleDistances,
        		int debugRow,
        		int debugColumn,
				int debugLevel){
    		int modDebugLevel=debugLevel;
    		int height=features.length;
    		int width=features[0].length;
    		int scanAroundCells=overlapFraction/2;
    		int gridStep=corrFFTSize/overlapFraction;
    		// find all cells where feature is in closer to the cell center, than matching feature - to the center of any other cell
    		int [][] theseHosts=  new int [(2*scanAroundCells+1)*(2*scanAroundCells+1)][2];
    		double[] theseShares=new double [(2*scanAroundCells+1)*(2*scanAroundCells+1)];
    		for (int row=0;row<height;row++){
        		for (int col=0;col<width;col++) if (features[row][col]!=null){
//        			modDebugLevel=debugLevel+(((row==debugRow) && (col==debugColumn))?2:0);
        			modDebugLevel=debugLevel+(((Math.abs(row-debugRow)<=2) && (Math.abs(col-debugColumn)<=2))?2:0);
        			double distance=scaleDistances*features[row][col][4];
        			double angle=   features[row][col][3];
        			double dX=distance*Math.cos(angle);
        			double dY=distance*Math.sin(angle);
        			int row0=row+ ((int) Math.round(dY/gridStep));
        			int col0=col+ ((int) Math.round(dX/gridStep));
        			int numHost=0;
        			double sumShares=0.0;
					if (modDebugLevel>2){
						int corrYC=row*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
						int corrXC=col*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
						System.out.println(">>>>++++ getFeatureShares()1: row="+row+" col="+col+" YC="+corrYC+" XC="+corrXC+
								" row0="+row0+" col0="+col0+
								" distance="+distance+
								" angle="+IJ.d2s(180/Math.PI*angle,2)+
								" dX="+dX+
								" dY="+dY);
					}
					for (int dRow=-scanAroundCells;dRow<=scanAroundCells;dRow++) {
						int row1=row0+dRow;
						if ((row1>=0) && (row1<height)) {
							for (int dCol=-scanAroundCells;dCol<=scanAroundCells;dCol++){
								int col1=col0+dCol;
								if (((row1!=row) || (col1!=col)) &&(col1>=0) && (col1<width)){
									if (modDebugLevel>3){
										int corrYC=row1*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
										int corrXC=col1*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
										System.out.println("getFeatureShares()1a: row1="+row1+" col1="+col1+" YC="+corrYC+" XC="+corrXC+
												" cellUsage[row1][col1]="+cellUsage[row1][col1]);
									}
									if (cellUsage[row1][col1]>0){ //cellUsage[row][col] - just to be sure, not needed
										double wDist=getFeatureMatch(
												row,
												col,
												row1,
												col1,
												gridStep,
												features,
												directionTolerance,
												normalDistanceTolerance,
												tangentialDistanceTolerance,
												directionSigma, // make a fixed fraction of directionTolerance?
												normalDistanceSigma, // make a fixed fraction of distanceTolerance?
												tangentialDistanceSigma, // make a fixed fraction of distanceTolerance?
												scaleDistances,
												debugLevel);
										if (wDist>0){
											theseHosts[numHost][0]=col1;
											theseHosts[numHost][1]=row1;
											double w=combinedStrength(strengthMode,features[row1][col1])*wDist;
											theseShares[numHost]=w;
											sumShares+=w;
											numHost++;
											if (modDebugLevel>2){
												int corrYC=row1*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
												int corrXC=col1*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
												System.out.println("getFeatureShares()2: row1="+row1+" col1="+col1+" YC="+corrYC+" XC="+corrXC+
														" wDist="+wDist+
														" w="+w+
														" sumShares="+sumShares+
														" numHost="+numHost);
											}
										} else {
											if (modDebugLevel>2){
												int corrYC=row1*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
												int corrXC=col1*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
												System.out.println("getFeatureShares()2: row1="+row1+" col1="+col1+" YC="+corrYC+" XC="+corrXC+
														" wDist="+wDist+
														" sumShares="+sumShares+
														" numHost="+numHost);
											}
										}

									}
								}
							}
						}
					}
					if (numHost>0) {
						hosts[row][col]=new int [numHost][2];
						shares[row][col]=new double [numHost];
						for (int i=0;i<numHost;i++){
							hosts[row][col][i]=theseHosts[i].clone();
							shares[row][col][i]=theseShares[i]/sumShares;
						}
					}
        		}
    		}
     	}

     	private int [][] getCellUsage(
    			int corrFFTSize,
    			int overlapFraction, // default 8
    			double [][][] features,
    			double directionTolerance,
    			double normalDistanceTolerance,
    			double tangentialDistanceTolerance,
     			double hostsTolerance,
    			double scaleDistances,
        		int debugRow,
        		int debugColumn,
				int debugLevel){
    		int modDebugLevel=debugLevel;
    		int height=features.length;
    		int width=features[0].length;
    		int [][] cellUsage=new int [height][width]; // 0 - unused, >0 - number of cells combined, <0: 1+y*width+x of the "host"
    		for (int row=0;row<height;row++){
        		for (int col=0;col<width;col++){
        			cellUsage[row][col]=0;
        		}
    		}
    		int gridStep=corrFFTSize/overlapFraction;
    		// find all cells where feature is in closer to the cell center, than matching feature - to the center of any other cell
    		int numHostCells=0;
    		for (int row=0;row<height;row++){
        		for (int col=0;col<width;col++) if (features[row][col]!=null){
        			modDebugLevel=debugLevel+(((Math.abs(row-debugRow)<=2) && (Math.abs(col-debugColumn)<=2))?2:0);
        			double distance=scaleDistances*features[row][col][4];
        			double angle=   features[row][col][3];
        			double dX=distance*Math.cos(angle);
        			double dY=distance*Math.sin(angle);
        			int row0=row+ ((int) Math.round(dY/gridStep));
        			int col0=col+ ((int) Math.round(dX/gridStep));
        			int range= (int) Math.ceil(distance/gridStep);
        			boolean best=true;
					if (modDebugLevel>2){
						int corrYC=row*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
						int corrXC=col*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
						System.out.println(">>>> getCellUsage(): row="+row+" col="+col+" YC="+corrYC+" XC="+corrXC+
								" row0="+row0+" col0="+col0+
								" dX="+IJ.d2s(dX,3)+" dY="+IJ.d2s(dY,3)+
								" range="+range+
								" angle="+IJ.d2s(180/Math.PI*angle,2)+
								" distance="+IJ.d2s(distance,3)+
								" strength="+IJ.d2s(features[row][col][0],4)+
								" refLev="+IJ.d2s(features[row][col][1],4));
					}
					for (int dRow=-range;(dRow<=range) && best;dRow++) {
						int row1=row0+dRow;
						if ((row1>=0) && (row1<height)) {
							for (int dCol=-range;dCol<=range;dCol++){
								int col1=col0+dCol;
								if (((row1!=row) || (col1!=col)) &&(col1>=0) && (col1<width) && (features[row1][col1]!=null)){ //cellUsage[row][col] - just to be sure, not needed
									double w=getFeatureMatch(
											row,
											col,
											row1,
											col1,
											gridStep,
											features,
											directionTolerance,
											normalDistanceTolerance,
											tangentialDistanceTolerance,
											-1, // just check tolerance
											-1, // just check tolerance
											-1, // just check tolerance
											scaleDistances,
											debugLevel);
									best&=(w<=0) || (Math.abs(distance)<=Math.abs(scaleDistances*features[row1][col1][4]));
    								if ((w>0) &&(modDebugLevel>2)){
    									int corrYC=row1*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
    									int corrXC=col1*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
    									System.out.println("getCellUsage(): row1="+row1+" col1="+col1+" YC="+corrYC+" XC="+corrXC+
    											" w="+w+
    											" best="+best+
    											" distance="+IJ.d2s(distance,3)+
    											" scaleDistances*features["+row1+"]["+col1+"]["+4+"]="+scaleDistances*features[row1][col1][4]);
    								}
									if (!best) break;
								}

							}
						}
					}
					if (best) {
        				cellUsage[row][col]=1;
        				numHostCells++;
					}
        		}
    		}
    		if (debugLevel>1){
    			System.out.println("mergeLinearFeatures(): Detected "+numHostCells+" host cells (common feature is closer to its center than to any other cell.");
    		}
    		return cellUsage;
     	}


    	/**
    	 *
    	 * @param corrFFTSize
    	 * @param overlapFraction
    	 * @param features
    	 * @param ignorePhase
    	 * @param strengthMode
    	 * @param phaseIntegrationWidth
    	 * @param resultHighPass
    	 * @param threadsMax
    	 * @param debugLevel
    	 * @return
    	 */
    	public double [] reconstructImageFeature(
    			final int corrFFTSize,
    			final int overlapFraction, // default 8
    			final double [][][] features,
    			final boolean ignorePhase,
    			final boolean preserveDC,
    			final double strengthMode, // 0.0 - absolute, 1.0 - relative, 0.5 - "balanced", "-1" - use phaseStrength instead
    			final double phaseIntegrationWidth, // use the same as during extraction?
    			final double resultHighPass, //cutting "ribbon" near zero frequency - influences length of the detected line
    			final int debugRow,
    			final int debugColumn,
				final int threadsMax,
				final int debugLevel){


//    		final int numTileRows=    this.mapHeight*overlapFraction/corrFFTSize+(((this.mapHeight*overlapFraction/corrFFTSize)*corrFFTSize/overlapFraction<this.mapHeight)?1:0);
//    		final int numTileColumns= this.mapWidth *overlapFraction/corrFFTSize+(((this.mapWidth *overlapFraction/corrFFTSize)*corrFFTSize/overlapFraction<this.mapWidth )?1:0);
    		final int mapHeight=this.mapHeight;
    		final int mapWidth= this.mapWidth;
    		final int numTileRows=    features.length;
    		final int numTileColumns= features[0].length;
    		final int numTiles= numTileRows*numTileColumns;
    		final Thread[] threads = newThreadArray(threadsMax);
    		final AtomicInteger tileNum =  new AtomicInteger(0);
    		final AtomicInteger tilesFinished = new AtomicInteger(0);
    		final double [] featuresPlot=new double [this.mapHeight*this.mapWidth];
    		for (int i=0;i<featuresPlot.length;i++) featuresPlot[i]=0.0;
    		final int hSize=corrFFTSize/2;
    		final AtomicInteger numNonEmpty =  new AtomicInteger(0);
    		for (int ithread = 0; ithread < threads.length; ithread++) {
    			threads[ithread] = new Thread() {
    				@Override
					public void run() {
    					DoubleFHT doubleFHT=new DoubleFHT();
    					double [] tilePlot;
    					for (int tile=tileNum.getAndIncrement(); tile<numTiles;tile=tileNum.getAndIncrement()) {
    						int row= tile/numTileColumns;
    						int col= tile%numTileColumns;
    						double [] feature=features[row][col];
    						if (feature!=null) {
    							numNonEmpty.getAndIncrement();
    							int corrYC=row*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
    							int corrXC=col*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
    							doubleFHT.debug=((row==debugRow) && (col==debugColumn));
    							tilePlot=reconstructTileFeature(
    									doubleFHT,
    									corrFFTSize,
    									feature,
    									ignorePhase,
    									preserveDC,
    									strengthMode, // 0.0 - absolute, 1.0 - relative, 0.5 - "balanced", "-1" - use phaseStrength instead
    									phaseIntegrationWidth, // use the same as during extraction?
    									resultHighPass, //cutting "ribbon" near zero frequency - influences length of the detected line
    									debugLevel + (doubleFHT.debug?1:0)
    							);
//								if ((row== 105) && (debugLevel>1)) {
//									System.out.println("reconstructImageFeature():  row="+row+" col="+col+" corrYC="+corrYC+" corrXC="+corrXC);
//								}
/*
    	   						if ((Math.abs(row-debugRow)<6) && (Math.abs(col-debugColumn)<6) && (feature!=null) && (debugLevel>1)){
    	   							System.out.println ("reconstructImageFeature(): tile="+tile+" ( of "+numTiles+"), row="+row+" col="+col+" corrYC="+corrYC+" corrXC="+corrXC+
    	   								" strength="+feature[0]+
    	   								" angle="+IJ.d2s(180/Math.PI*feature[3],2)+
    	   								" distance="+IJ.d2s(feature[4],2)+
    	   								" phase zero="+IJ.d2s(180/Math.PI*feature[5],2)
    	   								);
    	   						}
*/
    							synchronized(this){
    								int indexSrc=0;
    								for (int y=corrYC-hSize;y<(corrYC+hSize);y++){
    									if ((y>=0) && (y<mapHeight)) {
    										int indexDst=y*mapWidth+corrXC-hSize;
    										for (int x=corrXC-hSize;x<(corrXC+hSize);x++){
    											if ((x>=0) && (x<mapWidth)){
//    												if ((indexDst>featuresPlot.length) || (indexSrc>tilePlot.length)) {
//    													System.out.println("reconstructImageFeature(): x="+x+" y="+y+" row="+row+" col="+col+" corrYC="+corrYC+" corrXC="+corrXC+
//   															" indexDst="+indexDst+" featuresPlot.length="+featuresPlot.length+" indexSrc="+indexSrc+" tilePlot.length="+tilePlot.length);
//   												}
    												featuresPlot[indexDst]+=tilePlot[indexSrc];
    											}
    											indexSrc++;
    											indexDst++;
    										}
    									} else indexSrc+=corrFFTSize;
    								}
    							}
    							final int numFinished=tilesFinished.getAndIncrement();
    							SwingUtilities.invokeLater(new Runnable() {
    								@Override
									public void run() {
    									IJ.showProgress(numFinished,numTiles);
    								}
    							});
    						}
    					}
    				}
    			};
    		}
    		startAndJoin(threads);
    		IJ.showProgress(1.0);
    		if (debugLevel>1){
    			System.out.println("=== reconstructImageFeature(): number of non-empty cells="+numNonEmpty.get()+
    					" (of "+(numTileRows*numTileColumns)+"), "+100.0*numNonEmpty.get()/((numTileRows*numTileColumns))+"%");
    		}
    		return featuresPlot;
    	}


    	/**
    	 * Plot the extracted linear feature to a square tile
    	 * @param doubleFHT - instance of doubleFHT. if not null, will be reused and the cached data will not be re-calculated
    	 * @param corrFFTSize size of the square tile side
    	 * @param feature array of parameters describing the linear feature (see tileLineDetect() method)
    	 * @param ignorePhase  always plot feature with zero phase (white on black)
    	 * @param preserveDC  do not remove DC component from the plotted feature
    	 * @param strengthMode scale normalized linear feature to (0.0): absolute strength, 1.0 - relative strength (0.0<strengthMode<1.0): intermediate, (-1): phaseStrength, -2 - no scaling
    	 * @param phaseIntegrationWidth  use the same as during extraction?
    	 * @param resultHighPass (used only for strength calculation) - cut frequencies near zero for the linear phase band (shortens the result line);
    	 * @param debugLevel Debug level, turning on more text/graphic output when higher
    	 * @return feature plotted to a corrFFTSizexcorrFFTSize tile, linescan order
    	 */
    	public double [] reconstructTileFeature(
			DoubleFHT doubleFHT,
			int corrFFTSize,
			double [] feature,
			boolean ignorePhase,
			boolean preserveDC,
			double strengthMode, // 0.0 - absolute, 1.0 - relative, 0.5 - "balanced", "-1" - use phaseStrength instead
			double phaseIntegrationWidth, // use the same as during extraction?
			double resultHighPass, //cutting "ribbon" near zero frequency - influences length of the detected line
			int debugLevel
			){
    		double absoluteStrength=feature[0];
    		double referenceLevel=  feature[1];
    		double phaseStrength=   feature[2];
    		double angle=           feature[3];
    		double distance=        feature[4];
    		double phaseAtZero=     feature[5];
    		double scale=combinedStrength(strengthMode,feature);
//			(strengthMode<-1.0)?1.0:((strengthMode>0.0)?( absoluteStrength/Math.pow(referenceLevel,strengthMode)):phaseStrength);
    		if (ignorePhase)         phaseAtZero=0.0;
    		if (doubleFHT==null )doubleFHT= new DoubleFHT();
    		int length=corrFFTSize*corrFFTSize;
			int    [] bandIndices= doubleFHT.getBandIndices(
					preserveDC,
					phaseIntegrationWidth,
					angle);
			double [] bandMask= doubleFHT.getRibbonMask(
					preserveDC,
					bandIndices,
					phaseIntegrationWidth,
					resultHighPass,
					angle);

			double [] phaseStepCorr=doubleFHT.compensatePhaseStep(
					angle,
					phaseAtZero,
					distance,
					0, //phaseTolerance, // if >0 will zero amplitude if the phase is too off - - only used with non-zero modifiedAmpPhase
					null, // modifiedAmpPhase, // modifies phase and/or amplitude (if  phaseTolerance>0),
					bandIndices,
					0); // zeroBinHalfSize); - only used with non-zero modifiedAmpPhase
    		double[] hamming1d=doubleFHT.getHamming1d(corrFFTSize,0.0); // complete zero out on the border

    		double [] uniform=new double[length];
    		for (int iy=0;iy<corrFFTSize;iy++) for (int ix=0;ix<corrFFTSize;ix++){
    			int index=iy*corrFFTSize+ix;
    			uniform[index]=hamming1d[iy]*hamming1d[ix]*scale;
    		}
			doubleFHT.swapQuadrants(uniform);
			if (debugLevel>2){
				double [] debugBandMask=new double [corrFFTSize*corrFFTSize];
				for (int i=0;i<debugBandMask.length;i++) debugBandMask[i]=0;

				for (int i=0;i<bandIndices.length;i++){
					debugBandMask[bandIndices[i]]=bandMask[i];
				}
				(new ShowDoubleFloatArrays()).showArrays(
						debugBandMask,
						corrFFTSize,
						corrFFTSize,
						"bandMask");


				System.out.println("+++reconstructTileFeature(): "+
						" strengthMode="+strengthMode+
						" scale="+IJ.d2s(scale,3)+
						" absoluteStrength="+absoluteStrength+
						" referenceLevel="+referenceLevel+
						" phaseStrength="+phaseStrength+
						" angle="+IJ.d2s(180/Math.PI*angle,2)+
						" distance="+IJ.d2s(distance,2)+
						" phase zero="+IJ.d2s(180/Math.PI*phaseAtZero,2)

						);
				(new ShowDoubleFloatArrays()).showArrays(
						uniform,
						corrFFTSize,
						corrFFTSize,
						"uniform");

			}

			double [] reconstructedFeature=doubleFHT.reconstructLinearFeature(
					uniform,
					angle,
					distance,
					phaseIntegrationWidth,
					bandIndices,
					null, //modPhase,
					phaseStepCorr,
					bandMask);
			if (debugLevel>2){
				(new ShowDoubleFloatArrays()).showArrays(
						reconstructedFeature,
						corrFFTSize,
						corrFFTSize,
						"reconstructedFeature");

			}
			return reconstructedFeature;
    	}

    	private double combinedStrength(
    			double strengthMode,
    			double [] feature){
    		return (strengthMode<-1.0)?1.0:((strengthMode>0.0)?( feature[0]/Math.pow(feature[1],strengthMode)):feature[2]);
    	}

    	/**
    	 * Calculate parameters of the line features in the current image
    	 * @param side  0- left, 1 - right image
    	 * @param corrFFTSize size of the square tile side
    	 * @param overlapFraction overlap tiles by this fraction of the tile size (corrFFTSize)
    	 * @param alphaThreshold Minimal required overlap of the tile with image alpha (0 - do not calculate)
    	 * @param useBinaryAlpha consider all alpha>0.0 as 1.0
    	 * @param correlationHighPassSigma High pass filter before feature extraction (measured in frequency samples)
    	 * @param correlationLowPassSigma Low-pass filter for correlation (as fraction of the frequency range)
    	 * @param phaseIntegrationWidth - Width of the "band" when looking for the linear phase (actual will be smaller because of Hammimng window)
    	 * @param resultHighPass (used only for strength calculation) - cut frequencies near zero for the linear phase band (shortens the result line);
    	 * @param dispertionCost - >0 - use phase dispersion when calculating quality of phase match, 0 - only weights.
    	 * @param featureFilter Bitmask enabling different phase half-steps (+1 - 0, +2 - pi/2, +4 - pi, +8 - 3pi/2. Value zero allows arbitrary step
    	 *        step 0 corresponds to thin white line, pi - thin black line, +pi/2 - edge black-to-white in the direction of directionAngle,
    	 *        3pi/2 - white-to-black in the direction of directionAngle
    	 * @param minRMSFrac Minimal frequency sample value relative to RMS to be considered when looking for the linear phase feature (0.0 - skip test)
    	 * @param minAbs  Minimal frequency sample absolute value to be considered when looking for the linear phase feature  (0.0 - skip test)
    	 * @param maxPhaseMismatch calculate phase differences for each full (above threshold) set of 4 samples (in a 2x2 square) and measure the distance
    	 *        between two centers of diagonals, discard set if this distance is above threshold (0.0 - skip test)
    	 * @param calculateStrength If true, calculate feature strength by correlating the source image with the simulated linear phase band, if false (faster processing)
    	 *        these two fields will be returned equal 0.0
    	 * @param threadsMax maximal number of parallel threads to use
    	 * @param debugLevel Debug level, turning on more text/graphic output when higher
    	 * @return array [tileY][tileX] {feature parameters (or null if feature failed to be detected)}:
    	 * 		   0 - absolute strength of the feature (proportional to the image contrast)
    	 *		   1 - Strength reference level (RMS of the frequency samples amplitudes (in the selected band in the frequency domain)
    	 *		   2 - phaseStrength - composite value calculated when selecting the feature direction - this value is calculated even if calculateStrength is false
    	 *		   3 - angle from the tile center to the detected feature (perpendicular to the feature direction itself), zero - left, pi/2 - down (pixel coordiantes)
    	 *		   4 - distance in pixels from the tile center to the feature
    	 *		   5 - phaseAtZero Detected feature type (may be constrained by featureFilter - see that parameter for meaning of the different phase values
    	 */
    	public double [][][] detectLinearFeatures(
    			final boolean side,
    			final int corrFFTSize,
    			final int overlapFraction, // default 8
    			final double alphaThreshold, // minFracArea
    			final boolean useBinaryAlpha,
				final double correlationHighPassSigma,
				final double correlationLowPassSigma,
				final double phaseIntegrationWidth,
				final double resultHighPass, //cutting "ribbon" near zero frequency - influences length of the detected line (used only for strength calculation)
				final double dispertionCost,
				final int featureFilter,
				final double minRMSFrac,
				final double minAbs,
				final double maxPhaseMismatch,
				final boolean calculateStrength,
				final int debugRow,
				final int debugColumn,
				final int threadsMax,
				final int debugLevel) {
    		return detectLinearFeatures(
        			side?1:0,
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
    				threadsMax,
    				debugLevel);
    	}


    	public double [][][] detectLinearFeatures(
    			final int imageNumber, //boolean side,
    			final int corrFFTSize,
    			final int overlapFraction, // default 8
    			final double alphaThreshold, // minFracArea
    			final boolean useBinaryAlpha,
				final double correlationHighPassSigma,
				final double correlationLowPassSigma,
				final double phaseIntegrationWidth,
				final double resultHighPass, //cutting "ribbon" near zero frequency - influences length of the detected line (used only for strength calculation)
				final double dispertionCost,
				final int featureFilter,
				final double minRMSFrac,
				final double minAbs,
				final double maxPhaseMismatch,
				final boolean calculateStrength,
				final int debugRow,
				final int debugColumn,
				final int threadsMax,
				final int debugLevel){

    		if (debugLevel>1) System.out.println ("detectLinearFeatures("+imageNumber+",...)");
    		final int numTileRows=    this.mapHeight*overlapFraction/corrFFTSize+(((this.mapHeight*overlapFraction/corrFFTSize)*corrFFTSize/overlapFraction<this.mapHeight)?1:0);
    		final int numTileColumns= this.mapWidth *overlapFraction/corrFFTSize+(((this.mapWidth *overlapFraction/corrFFTSize)*corrFFTSize/overlapFraction<this.mapWidth )?1:0);
    		final int numTiles= numTileRows*numTileColumns;
    		final Thread[] threads = newThreadArray(threadsMax);
    		final AtomicInteger tileNum =  new AtomicInteger(0);
    		final AtomicInteger tilesFinished = new AtomicInteger(0);
    		final double [][][] result=new double[numTileRows][numTileColumns][];
	   		for (int ithread = 0; ithread < threads.length; ithread++) {
	   			threads[ithread] = new Thread() {
	   				@Override
					public void run() {
	   					DoubleFHT doubleFHT=new DoubleFHT();
	   					for (int tile=tileNum.getAndIncrement(); tile<numTiles;tile=tileNum.getAndIncrement()) {
	   						int row= tile/numTileColumns;
	   						int col= tile%numTileColumns;
	   						int corrYC=row*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
	   						int corrXC=col*corrFFTSize/overlapFraction+corrFFTSize/overlapFraction/2;
	   						if (debugLevel>2) System.out.println ("detectLinearFeatures(): tile="+tile+" ( of "+numTiles+"), row="+row+" col="+col+" corrYC="+corrYC+" corrXC="+corrXC);
	   						result[row][col]= tileLineDetect(
	   				    			doubleFHT,
	   				    			imageNumber, //side,
	   				    			alphaThreshold,
	   				    			useBinaryAlpha,
	   								corrFFTSize,
	   								corrXC,
	   								corrYC,
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
	   								debugLevel);
    						final int numFinished=tilesFinished.getAndIncrement();
    						SwingUtilities.invokeLater(new Runnable() {
    							@Override
								public void run() {
    								IJ.showProgress(numFinished,numTiles);
    							}
    						});
//	   						if ((result[row][col]!=null) && (debugLevel>1)) System.out.println ("detectLinearFeatures(): tile="+tile+" ( of "+numTiles+"), row="+row+" col="+col+" corrYC="+corrYC+" corrXC="+corrXC+
//	   								" strength="+result[row][col][0]);
	   						if ((row==debugRow) && ((col==debugColumn)) && (result[row][col]!=null) && (debugLevel>1)){
	   							System.out.println ("detectLinearFeatures(): tile="+tile+" ( of "+numTiles+"), row="+row+" col="+col+" corrYC="+corrYC+" corrXC="+corrXC+
	   								" strength="+result[row][col][0]+
	   								" angle="+IJ.d2s(180/Math.PI*result[row][col][3],2)+
	   								" distance="+IJ.d2s(result[row][col][4],2)+
	   								" phase zero="+IJ.d2s(180/Math.PI*result[row][col][5],2)+
	   								" featureFilter="+featureFilter
	   								);
	   						}

	   					}
	   				}

	   			};
	   		}
	   		startAndJoin(threads);
	   		IJ.showProgress(1.0);
	   		return result;
    	}


    	/**
    	 * Calculate parameters of the line feature in the current image tile
    	 * @param doubleFHT - instance of doubleFHT. if not null, will be reused and the cached data will not be re-calculated
    	 * @param side  0- left, 1 - right image
    	 * @param alphaThreshold Minimal required overlap of the tile with image alpha (0 - do not calculate)
    	 * @param useBinaryAlpha consider all alpha>0.0 as 1.0
    	 * @param corrFFTSize size of the square tile side
    	 * @param corrXC Tile center X
    	 * @param corrYC Tile center Y
    	 * @param correlationHighPassSigma High pass filter before feature extraction (measured in frequency samples)
    	 * @param correlationLowPassSigma Low-pass filter for correlation (as fraction of the frequency range)
    	 * @param phaseIntegrationWidth - Width of the "band" when looking for the linear phase (actual will be smaller because of Hammimng window)
    	 * @param resultHighPass (used only for strength calculation) - cut frequencies near zero for the linear phase band (shortens the result line);
    	 * @param dispertionCost - >0 - use phase dispersion when calculating quality of phase match, 0 - only weights.
    	 * @param featureFilter Bitmask enabling different phase half-steps (+1 - 0, +2 - pi/2, +4 - pi, +8 - 3pi/2. Value zero allows arbitrary step
    	 *        step 0 corresponds to thin white line, pi - thin black line, +pi/2 - edge black-to-white in the direction of directionAngle,
    	 *        3pi/2 - white-to-black in the direction of directionAngle
    	 * @param minRMSFrac Minimal frequency sample value relative to RMS to be considered when looking for the linear phase feature (0.0 - skip test)
    	 * @param minAbs  Minimal frequency sample absolute value to be considered when looking for the linear phase feature  (0.0 - skip test)
    	 * @param maxPhaseMismatch calculate phase differences for each full (above threshold) set of 4 samples (in a 2x2 square) and measure the distance
    	 *        between two centers of diagonals, discard set if this distance is above threshold (0.0 - skip test)
    	 * @param calculateStrength If true, calculate feature strength by correlating the source image with the simulated linear phase band, if false (faster processing)
    	 *        these two fields will be returned equal 0.0
    	 * @param debugLevel Debug level, turning on more text/graphic output when higher
    	 * @return array with the following result fields (or null if feature failed to be detected):
    	 * 		   0 - absolute strength of the feature (proportional to the image contrast)
    	 *		   1 - Strength reference level (RMS of the frequency samples amplitudes (in the selected band in the frequency domain)
    	 *		   2 - phaseStrength - composite value calculated when selecting the feature direction - this value is calculated even if calculateStrength is false
    	 *		   3 - angle from the tile center to the detected feature (perpendicular to the feature direction itself), zero - left, pi/2 - down (pixel coordiantes)
    	 *		   4 - distance in pixels from the tile center to the feature
    	 *		   5 - phaseAtZero Detected feature type (may be constrained by featureFilter - see that parameter for meaning of the different phase values
    	 */
    	public double [] tileLineDetect(
    			DoubleFHT doubleFHT,
    			int imageNumber, //boolean side,
    			double alphaThreshold,
    			final boolean useBinaryAlpha,
				int corrFFTSize,
				int corrXC,
				int corrYC,
				double correlationHighPassSigma,
				double correlationLowPassSigma,
				double phaseIntegrationWidth,
				double resultHighPass, //cutting "ribbon" near zero frequency - influences length of the detected line (used only for strength calculation)
				double dispertionCost,
				int featureFilter,
				double minRMSFrac,
				double minAbs,
				double maxPhaseMismatch,
				boolean calculateStrength,
				int debugLevel){
    		if (doubleFHT==null )doubleFHT= new DoubleFHT();
    		int length=corrFFTSize*corrFFTSize;
    		int numSensors=this.channel.length;
    		int numLayers=this.overlapImages.length/numSensors;
    		double [][] selection = getSelection(
    				0, // cross-correlation always
    				corrFFTSize,
    				corrXC,
    				corrYC,
    				0); // positive - shift second image, negative - shift first image
    		double[] hamming1d=doubleFHT.getHamming1d(corrFFTSize,0.0); // complete zero out on the border
    		double [] window=new double[length]; // TODO:cache window in doubleFHT instance?
    		for (int iy=0;iy<corrFFTSize;iy++) for (int ix=0;ix<corrFFTSize;ix++){
    			int index=iy*corrFFTSize+ix;
    			window[index]=hamming1d[iy]*hamming1d[ix];
    		}
//    		double [][] corr=new double[numLayers][length];
//    		for (int i=0;i<length;i++) corr[0][i]=0;
//    		int alphaIndex=side?numLayers:0;
    		int alphaIndex=imageNumber*numLayers;

    		double [] tileData=  selection[alphaIndex+1].clone();
    		double [] tileAlpha= selection[alphaIndex+0];
    		if (alphaThreshold>0.0){
    			double sumAlpha=0.0;
    			double sumWindow=0.0;
    			for (int index=0;index<length;index++){
    				double alpha=useBinaryAlpha?1.0:tileAlpha[index];
    				sumWindow+=window[index];
    				sumAlpha+=window[index]*alpha;
    			}
    			if ((sumAlpha/sumWindow) < alphaThreshold) return null;
    		}
//    		System.out.print(corrXC+"/"+corrYC+"    ");
    		normalizeAndWindow(tileData, window,true);
    		double [] frequencyFilter=null;
    		if ((correlationHighPassSigma>0.0) || (correlationLowPassSigma>0.0)){
				frequencyFilter=doubleFHT.createFrequencyFilter(
	    				tileData, // just for length
	    				correlationHighPassSigma,
	    				correlationLowPassSigma);
    		}
    		//calculateAmplitude
			doubleFHT.swapQuadrants(tileData);
			doubleFHT.transform(tileData,false);
			if (frequencyFilter!=null){
    			doubleFHT.multiplyByReal(tileData, frequencyFilter);
    		}
    		doubleFHT.debug=(debugLevel>3);
			double [] dirDistStrength= doubleFHT.calcPhaseApproximation(
					phaseIntegrationWidth,
					tileData,
					minRMSFrac,
					minAbs,
					maxPhaseMismatch,
					dispertionCost,
					featureFilter);
			if (dirDistStrength==null) return null;

			double angle=        dirDistStrength[0];
			double distance=     dirDistStrength[1];
			double phaseAtZero=  dirDistStrength[2];
			double phaseStrength=dirDistStrength[3];
			if (distance<0){
				angle+=Math.PI;
				angle+=(2*Math.PI)*Math.round(angle/(2*Math.PI));
				phaseAtZero=-phaseAtZero;
				distance=-distance;
			}

			double [] result= {
					Double.NaN, //strengths[0], // absolute strength
					Double.NaN, //strengths[1], // reference level
					phaseStrength, // phase strength (less computations than the first two, known anyway)
					angle,
					distance,
					phaseAtZero
			};
			if (!calculateStrength) return result;

			// needed data is already received, next just for debugging. Not only - also to measure the strength
			int    [] bandIndices= doubleFHT.getBandIndices(
					false,
					phaseIntegrationWidth,
					angle);
			double [] bandMask= doubleFHT.getRibbonMask(
					false,
					bandIndices,
					phaseIntegrationWidth,
					resultHighPass,
					angle);

			double [] phaseStepCorr=doubleFHT.compensatePhaseStep(
					angle,
					phaseAtZero,
					distance,
					0, //phaseTolerance, // if >0 will zero amplitude if the phase is too off - - only used with non-zero modifiedAmpPhase
					null, // modifiedAmpPhase, // modifies phase and/or amplitude (if  phaseTolerance>0),
					bandIndices,
					0); // zeroBinHalfSize); - only used with non-zero modifiedAmpPhase
			double [] strengths=doubleFHT.linearFeatureStrength(
					tileData,
					angle,
					distance,
					phaseIntegrationWidth,
					bandIndices,
					null, // modifiedAmpPhase,
					phaseStepCorr,
					bandMask);
			result[0]=strengths[0]; // absolute strength
			result[1]=strengths[1]; // reference level

			if (debugLevel>2){
				System.out.println("balanced strength="+IJ.d2s(Math.sqrt(strengths[0]*strengths[0]/strengths[1]),3));
				System.out.println("absolute strength="+IJ.d2s(strengths[0],3));
				System.out.println("relative strength="+IJ.d2s((strengths[0]/strengths[1]),3));
				System.out.println("   phase strength="+IJ.d2s(phaseStrength,3));
				System.out.println("            Angle="+IJ.d2s(180*angle/Math.PI,2)+" degrees");
				System.out.println("         distance="+IJ.d2s(distance,2)+" pixels");
				System.out.println("       zero phase="+IJ.d2s(180*phaseAtZero/Math.PI,2)+" degrees");
			}
			return result;
    	}


    	public void testLineDetect(
				int corrFFTSize,
				int corrXC,
				int corrYC,
				double phaseCorrelationFraction,
				double correlationHighPassSigma,
				double correlationLowPassSigma,
				boolean removeIslands,
				double threshold, // relative to RMS value
				boolean adjustDirection,
				double phaseIntegrationWidth,
				double resultHighPass,
				double dispertionCost,
				int featureFilter,
				double zeroBinHalfSize,
				double minRMSFrac,
				double minAbs,
				double maxPhaseMismatch,
				double phaseTolerance, //simulation phase tolerance
				int debugLevel){
    		int length=corrFFTSize*corrFFTSize;
    		DoubleFHT doubleFHT= new DoubleFHT();
    		int numLayers=this.overlapImages.length/2;
    		double [][] selection = getSelection(
    				1, // autocorrelation always
    				corrFFTSize,
    				corrXC,
    				corrYC,
    				0); // positive - shift second image, negative - shift first image
    		double[] hamming1d=(new DoubleFHT()).getHamming1d(corrFFTSize,0.0); // complete zero out on the border
    		double [] window=new double[length];
    		for (int iy=0;iy<corrFFTSize;iy++) for (int ix=0;ix<corrFFTSize;ix++){
    			int index=iy*corrFFTSize+ix;
    			window[index]=hamming1d[iy]*hamming1d[ix];
    		}
			double [] ones=window.clone(); // needed to visialize only
			doubleFHT.swapQuadrants(ones); // needed to visialize only

    		double [][] corr=new double[numLayers][length];
    		for (int i=0;i<length;i++) corr[0][i]=0;
    		double [] first= selection[1].clone();
    		normalizeAndWindow(first, window,true);
    		double [] debugWindowed=first.clone();
    		if (debugLevel>2){
    			(new ShowDoubleFloatArrays()).showArrays(
    					first,
    					corrFFTSize,
    					corrFFTSize,
    					"windowed-Y"+corrXC+"-y"+corrYC+"-PC"+phaseCorrelationFraction);
    		}
    		double [] frequencyFilter=null;
    		if ((correlationHighPassSigma>0.0) || (correlationLowPassSigma>0.0)){
				frequencyFilter=doubleFHT.createFrequencyFilter(
	    				first, // just for length
	    				correlationHighPassSigma,
	    				correlationLowPassSigma);
    		}
    		//calculateAmplitude
			doubleFHT.swapQuadrants(first);
			doubleFHT.transform(first,false);

			if (frequencyFilter!=null){
    			doubleFHT.multiplyByReal(first, frequencyFilter);
    			if (debugLevel>2){
    				double [] filteredAmplitude=doubleFHT.calculateAmplitude(first);
    				(new ShowDoubleFloatArrays()).showArrays(
    						filteredAmplitude,
    						corrFFTSize,
    						corrFFTSize,
    						"filteredAmp-Y"+corrXC+"-y"+corrYC+"-PC"+phaseCorrelationFraction);
    				for (int i=0;i<filteredAmplitude.length;i++) filteredAmplitude[i]=Math.log(filteredAmplitude[i]); // NaN OK
    				(new ShowDoubleFloatArrays()).showArrays(
    						filteredAmplitude,
    						corrFFTSize,
    						corrFFTSize,
    						"filteredLog-Y"+corrXC+"-y"+corrYC+"-PC"+phaseCorrelationFraction);
    			}
    		}
    		doubleFHT.debug=(debugLevel>1);
			double [] dirDistStrength= doubleFHT.calcPhaseApproximation(
					phaseIntegrationWidth,
					first,
					minRMSFrac,
					minAbs,
					maxPhaseMismatch,
					dispertionCost,
					featureFilter);
			if (dirDistStrength==null) return;
			double angle=        dirDistStrength[0];
			double distance=     dirDistStrength[1];
			double phaseAtZero=  dirDistStrength[2];
			double phaseStrength=dirDistStrength[3];

			// needed data is already received, next just for debugging
			int    [] bandIndices= doubleFHT.getBandIndices(
					false,
					phaseIntegrationWidth,
					angle);
			double [] bandMask= doubleFHT.getRibbonMask(
					false,
					bandIndices,
					phaseIntegrationWidth,
					resultHighPass,
					angle);

//			double [][] modifiedAmpPhase=doubleFHT.calcIndAmHPhase (first, bandIndices);
			double [] phaseStepCorr=doubleFHT.compensatePhaseStep(
					angle,
					phaseAtZero,
					distance,
					phaseTolerance, // if >0 will zero amplitude if the phase is too off
					null, // modifiedAmpPhase, // modifies phase and/or amplitude (if  phaseTolerance>0),
					bandIndices,
					zeroBinHalfSize);
			double [] strengths=doubleFHT.linearFeatureStrength(
					first,
					angle,
					distance,
					phaseIntegrationWidth,
					bandIndices,
					null, // modifiedAmpPhase,
					phaseStepCorr,
					bandMask);
			if (debugLevel>0){
				System.out.println("balanced strength="+IJ.d2s(Math.sqrt(strengths[0]*strengths[0]/strengths[1]),3));
				System.out.println("absolute strength="+IJ.d2s(strengths[0],3));
				System.out.println("relative strength="+IJ.d2s((strengths[0]/strengths[1]),3));
				System.out.println("   phase strength="+IJ.d2s(phaseStrength,3));
				System.out.println("            Angle="+IJ.d2s(180*angle/Math.PI,2)+" degrees");
				System.out.println("         distance="+IJ.d2s(distance,2)+" pixels");
				System.out.println("       zero phase="+IJ.d2s(180*phaseAtZero/Math.PI,2)+" degrees");
			}
			if (debugLevel>0){

				double [] reconstructed=doubleFHT.reconstructLinearFeature(
						ones,
						angle,
						distance,
						phaseIntegrationWidth,
						bandIndices,
						null, //modPhase,
						phaseStepCorr,
						bandMask);

				//debugWindowed
				double [][] reconstructionComparison={reconstructed, debugWindowed};
				String [] titles={"reconstructed","windowed original"};
    			(new ShowDoubleFloatArrays()).showArrays(
    					reconstructionComparison,
    					corrFFTSize,
    					corrFFTSize,
    					true,
    					"reconstructed-x"+corrXC+"-y"+corrYC,
    					titles);
			}
    	}




    	public void testLineDetectDebug(
				int corrFFTSize,
				int corrXC,
				int corrYC,
				double phaseCorrelationFraction,
				double correlationHighPassSigma,
				double correlationLowPassSigma,
				boolean removeIslands,
				double threshold, // relative to RMS value
				boolean adjustDirection,
				double phaseIntegrationWidth,
				double resultHighPass,
				double dispertionCost,
				int featureFilter,
				double zeroBinHalfSize,
				double minRMSFrac,
				double minAbs,
				double maxPhaseMismatch,
				double phaseTolerance, //simulation phase tolerance
				int debugLevel){
    		int length=corrFFTSize*corrFFTSize;
    		DoubleFHT doubleFHT= new DoubleFHT();

//    		double [][] selection = new double[this.overlapImages.length][length];
    		int numLayers=this.overlapImages.length/2;
    		double [][] selection = getSelection(
    				1, // autocorrelation always
    				corrFFTSize,
    				corrXC,
    				corrYC,
    				0); // positive - shift second image, negative - shift first image
// correlation stuff
//    		double[] hamming1d=(new DoubleFHT()).getHamming1d(corrFFTSize);
    		double[] hamming1d=(new DoubleFHT()).getHamming1d(corrFFTSize,0.0); // complete zero out on the border
    		double [] window=new double[length];
    		for (int iy=0;iy<corrFFTSize;iy++) for (int ix=0;ix<corrFFTSize;ix++){
    			int index=iy*corrFFTSize+ix;
    			window[index]=hamming1d[iy]*hamming1d[ix];
    		}
    		double [][] corr=new double[numLayers][length];
    		for (int i=0;i<length;i++) corr[0][i]=0;
    		double [] first= selection[1].clone();
    		normalizeAndWindow(first, window,true);
    		double [] debugWindowed=first.clone();
    		if (debugLevel>2){
    			(new ShowDoubleFloatArrays()).showArrays(
    					first,
    					corrFFTSize,
    					corrFFTSize,
    					"windowed-Y"+corrXC+"-y"+corrYC+"-PC"+phaseCorrelationFraction);
    		}
    		double [] frequencyFilter=null;
    		if ((correlationHighPassSigma>0.0) || (correlationLowPassSigma>0.0)){
				frequencyFilter=doubleFHT.createFrequencyFilter(
	    				first, // just for length
	    				correlationHighPassSigma,
	    				correlationLowPassSigma);
    		}
    		//calculateAmplitude
			doubleFHT.swapQuadrants(first);
			doubleFHT.transform(first,false);



    		if (debugLevel>2){
    			double [] unfilteredAmplitude=doubleFHT.calculateAmplitude(first);
    			for (int i=0;i<unfilteredAmplitude.length;i++) unfilteredAmplitude[i]=Math.log(unfilteredAmplitude[i]); // NaN OK
    			(new ShowDoubleFloatArrays()).showArrays(
    					unfilteredAmplitude,
    					corrFFTSize,
    					corrFFTSize,
    					"unfilteredLog-Y"+corrXC+"-y"+corrYC+"-PC"+phaseCorrelationFraction);
    		}

    		if (frequencyFilter!=null){
    			doubleFHT.multiplyByReal(first, frequencyFilter);
    			if (debugLevel>2){
    				double [] filteredAmplitude=doubleFHT.calculateAmplitude(first);
    				(new ShowDoubleFloatArrays()).showArrays(
    						filteredAmplitude,
    						corrFFTSize,
    						corrFFTSize,
    						"filteredAmp-Y"+corrXC+"-y"+corrYC+"-PC"+phaseCorrelationFraction);
    				for (int i=0;i<filteredAmplitude.length;i++) filteredAmplitude[i]=Math.log(filteredAmplitude[i]); // NaN OK
    				(new ShowDoubleFloatArrays()).showArrays(
    						filteredAmplitude,
    						corrFFTSize,
    						corrFFTSize,
    						"filteredLog-Y"+corrXC+"-y"+corrYC+"-PC"+phaseCorrelationFraction);
    			}
    		}
    		doubleFHT.debug=(debugLevel>1);
			double [] dirDistStrength= doubleFHT.calcPhaseApproximation(
					phaseIntegrationWidth,
					first,
					minRMSFrac,
					minAbs,
					maxPhaseMismatch,
					dispertionCost,
					featureFilter);
			if (dirDistStrength==null) return;
			double angle=        dirDistStrength[0];
			double distance=     dirDistStrength[1];
			double phaseAtZero=  dirDistStrength[2];
//			double phaseStrength=dirDistStrength[3];

			// needed data is already received, next just for debugging
			int    [] bandIndices= doubleFHT.getBandIndices(
					false,
					phaseIntegrationWidth,
					angle);
			double [] bandMask= doubleFHT.getRibbonMask(
					false,
					bandIndices,
					phaseIntegrationWidth,
					resultHighPass,
					angle);


			double [][] modifiedAmpPhase=doubleFHT.calcIndAmHPhase (first, bandIndices);


			double [] phaseStepCorr=doubleFHT.compensatePhaseStep(
					angle,
					phaseAtZero,
					distance,
					phaseTolerance, // if >0 will zero amplitude if the phase is too off
					modifiedAmpPhase, // modifies phase and/or amplitude (if  phaseTolerance>0),
					bandIndices,
					zeroBinHalfSize);


			if (debugLevel>0){

				double [][] dbgMask=new double [2][corrFFTSize*corrFFTSize/2];
				for (int i=0;i<dbgMask[0].length;i++){
					dbgMask[0][i]=0.0;
					dbgMask[1][i]=0.0;
				}
				for (int n=0;n<bandIndices.length;n++){
					int index=bandIndices[n];
					dbgMask[0][index]=bandMask[n];
					dbgMask[1][index]=1.0;
				}

    			(new ShowDoubleFloatArrays()).showArrays(
    					dbgMask,
    					corrFFTSize,
    					corrFFTSize/2,
    					true,
    					"bandMask");
//    			doubleFHT.debug=false;
				double [] reconstructed=doubleFHT.reconstructLinearFeature(
						first,
						angle,
						distance,
						phaseIntegrationWidth,
						bandIndices,
						modifiedAmpPhase,
						phaseStepCorr,
						bandMask);
				double [] reconstructedNoBand=doubleFHT.reconstructLinearFeature(
						first,
						angle,
						distance,
						phaseIntegrationWidth,
						bandIndices,
						modifiedAmpPhase,
						phaseStepCorr,
						null);
//				double [][] modPhase=new double[modifiedAmpPhase[0].length][2];
//				for (int i=0;i<modPhase[0].length;i++){
//					modPhase[i][0]=1.0;
//					modPhase[i][1]=modifiedAmpPhase[i][1];
//				}
//				double [] ones = new double[first.length];
//				for (int i=0;i<ones.length;i++)ones[i]=1.0;
//				System.out.println("ones.length="+ones.length);
//				System.out.println("first.length="+first.length);
 //   			doubleFHT.debug=true;
				double [] ones=window.clone();
				doubleFHT.swapQuadrants(ones);

				double [] reconstructed1=doubleFHT.reconstructLinearFeature(
						ones,
						angle,
						distance,
						phaseIntegrationWidth,
						bandIndices,
						null, //modPhase,
						phaseStepCorr,
						bandMask);
				double [] reconstructedNoBand1=doubleFHT.reconstructLinearFeature(
						ones,
						angle,
						distance,
						phaseIntegrationWidth,
						bandIndices,
						null, //modPhase,
						phaseStepCorr,
						null);


				//debugWindowed
				double [][] reconstructionComparison={reconstructed,reconstructedNoBand,reconstructed1,reconstructedNoBand1,debugWindowed};
				System.out.println("reconstructionComparison[0].length="+reconstructionComparison[0].length);
				System.out.println("reconstructionComparison[1].length="+reconstructionComparison[1].length);
				System.out.println("reconstructionComparison[2].length="+reconstructionComparison[2].length);
				System.out.println("reconstructionComparison[3].length="+reconstructionComparison[3].length);
				System.out.println("reconstructionComparison[4].length="+reconstructionComparison[4].length);
				String [] titles={"reconstructed","sharp-reconstructed","simulated","simulated-reconstructed","windowed original"};
    			(new ShowDoubleFloatArrays()).showArrays(
    					reconstructionComparison,
    					corrFFTSize,
    					corrFFTSize,
    					true,
    					"reconstructed-x"+corrXC+"-y"+corrYC,
    					titles);


			}
			double [] strengths=doubleFHT.linearFeatureStrength(
					first,
					angle,
					distance,
					phaseIntegrationWidth,
					bandIndices,
					modifiedAmpPhase,
					phaseStepCorr,
					bandMask);
			if (debugLevel>0){
				System.out.println(" ====== Linear feature strength: absolute="+strengths[0]+" relative="+(strengths[0]/strengths[1]));
			}
    	}

    	public double [] correlateRectangular(
    			boolean autoCorrelation,
				int corrFFTSize,
				int overlapFraction,
				int corrXC,
				int corrYC,
				int maxDisparity,
				double phaseCorrelationFraction,
				double corrCbWeight,
				double corrCrWeight,
				double correlationHighPassSigma,
				double correlationLowPassSigma,
				double noiseNormalizationSignaY,
				double noiseNormalizationSignaCbCr,
				double contrastThreshold,
				boolean useBinaryAlpha,
				boolean enableNegativeDisparity,
				int threadsMax,
				int debugLevel){
    		int numTiles=1+(maxDisparity/(corrFFTSize/overlapFraction));
    		int length=corrFFTSize*corrFFTSize;
    		int corrWidth=corrFFTSize + corrFFTSize*(numTiles-1)/overlapFraction;
    		int outLength=corrFFTSize*corrWidth;
    		double [] fullCorr=new double [outLength];
    		for (int i=0;i<outLength;i++) fullCorr[i]=0.0;
    		int numLayers=this.overlapImages.length/2;
    		double[] hamming1d=(new DoubleFHT()).getHamming1d(corrFFTSize,0.0); // pure shifted cosine, not real Hamming
    		double [][] corr=new double[numLayers][length];
			double [][] debugCorr=new double[3][];
			double [] componentWeights={
					1.0/(1.0+corrCbWeight+corrCrWeight),
					corrCbWeight/(1.0+corrCbWeight+corrCrWeight),
					corrCrWeight/(1.0+corrCbWeight+corrCrWeight)};
    		for (int nTile=0;nTile<numTiles;nTile++){
        		double [][] selection = getSelection(
        				autoCorrelation?1:0,
        						corrFFTSize,
        						corrFFTSize,
        						corrXC, // -corrFFTSize/2*(numTiles-1),
        						corrYC,
        						corrFFTSize*nTile/overlapFraction);// positive - shift second image, negative - shift first image
        		double [] stats= selectionStats(selection,useBinaryAlpha);
        		if (debugLevel>0){
        			System.out.println("Selection ("+(useBinaryAlpha?"binary":"analog")+" alpha ) xc="+corrXC+" yc="+corrYC+" size="+corrFFTSize+" maxDisparity="+maxDisparity+
        					" : first "+IJ.d2s(100*stats[0],3)+"%, second "+IJ.d2s(100*stats[1],3)+"%, overlap "+IJ.d2s(100*stats[2],3)+"%");
        		}
        		double [][] window={selection[0].clone(),selection[numLayers].clone()};
        		for (int iy=0;iy<corrFFTSize;iy++) for (int ix=0;ix<corrFFTSize;ix++){
        			int index=iy*corrFFTSize+ix;
        			double h=hamming1d[iy]*hamming1d[ix];
        			window[0][index]*=h;
        			window[1][index]*=h;
        		}
    			for (int i=0;i<length;i++) corr[0][i]=0;
    			for (int i=0;i<numLayers-1;i++){
    				double [] first= selection[i+1].clone();
    				double [] second=selection[i+1+numLayers].clone();
    				normalizeAndWindow(first, window[0],true);
    				normalizeAndWindow(second,window[1],true);
            		if ((i==0) &&(debugLevel>1)){
            			debugCorr[1]=first.clone();
            			debugCorr[2]=second.clone();
            		}
    				// TODO: use common (per-thread) DoubleFHT instance to use caching of sin,cos, filter tables
    				corr[i+1]=(new DoubleFHT()).correlate (
    						first,
    						second,
    						correlationHighPassSigma,
    						correlationLowPassSigma,
    						phaseCorrelationFraction);
    				double sigma=(i==0)?noiseNormalizationSignaY:noiseNormalizationSignaCbCr;
    				if (sigma>0){
    					double rms=hFNoiseNormalize(
    							corr[i+1], // double [] data,
    							sigma, // double filterSigma,
    							true);  // boolean centerOnly);
    					if (this.debugLevel>1){
    						System.out.println("Correlation RMS for component "+((i==0)?"Y":((i==1)?"Cb":"Cr"))+ " was "+rms);
    					}
    				}
    				//        		for (int j=0;j<length;j++) corr[0][j]+=corr[i][j];
    				for (int j=0;j<length;j++) corr[0][j]+=componentWeights[i]*corr[i+1][j];
    			}
        		if (debugLevel>2){
        			String [] titles={"corr","first","second"};
        			debugCorr[0]=corr[0];

        			(new ShowDoubleFloatArrays()).showArrays(
        					debugCorr,
        					corrFFTSize,
        					corrFFTSize,
        					true,
        			"corr"+nTile+"_selection-x"+corrXC+"-y"+corrYC+"-MAXDSP"+maxDisparity+"-PC"+phaseCorrelationFraction,
        			titles);
        		}

				int srcIndex=0;
				for (int y=0;y<corrFFTSize;y++){
					int dstIndex=corrWidth*y+ (nTile)*corrFFTSize/overlapFraction;
					for (int x=0;x<corrFFTSize;x++){
						fullCorr[dstIndex++]+=corr[0][srcIndex++];
					}
				}

    		}
    		return fullCorr;
    	}

    	/**
    	 * Create a stack of channel images (YCbCr and optionally with external data
    	 * @param xc first image center X
    	 * @param yc first image center Y
    	 * @param size square side, power of 2 (will use twice larger internally)
    	 * @param first image number of the first in pair
    	 * @param second image number of the second in pair
    	 * @param dx - how much to shift the second image from the first, horizontally, right
    	 * @param dy - how much to shift the second image from the first, vertically,down
    	 * @param channelMask - channles to use ! - Y, 2 - Cb, 4-Cr
    	 * @param externalData array of external image data, first index image number, secon - pixel number
    	 * @param doubleSizeOutput - output twice the required size tile (to be able to probe around needed pixels
    	 * @param window - 1d in linescan order, (2*size)*(2*size). Probably - flat==1.0 in the inner size*size area. May be null - will be generated internally
    	 * @param doubleFHT DoubleFHT or null (to save resources)
    	 * @param debugLevel debug level
    	 * @return 2d array of the shifted slices from 2 images (only enabled in the channelMask and externalData first image-first channel, second image fiurst channel, ...
    	 * first tile (and the result will be shifted half fractional distance to the sescond, second - all intyeger part and half fractional towards the first
    	 */
    	public double [][] getShiftedSlices(
    			int xc1,
    			int yc1,
    			int size,
    			int first,
    			int second,
    			double dxA,
    			double dyA,
    			int channelMask,
    			double [][] externalData,
    			boolean doubleSizeOutput,
    			double [] window, // should be 4*size*size long;
    			DoubleFHT doubleFHT, // to reuse tables
    			int debugLevel
    			){
    		int numImgChannels=4; // alpha-Y-Cb-Cr
    		if (doubleFHT==null) doubleFHT=new DoubleFHT();
    		int iDx= (int) Math.round(dxA); //
    		int iDy= (int) Math.round(dyA);
    		int [] xc={xc1,xc1+ iDx};
    		int [] yc={yc1,yc1+ iDy};
    		double [] dx={-(dxA-iDx)/2,(dxA-iDx)/2}; // sub-pixel shift of the first image, second image
    		double [] dy={-(dyA-iDy)/2,(dyA-iDy)/2};

    		boolean [] shiftImage= {(dx[0]!=0.0) || (dy[0]!=0.0), (dx[1]!=0.0) || (dy[1]!=0.0)};
    		int numLayers=0;
    		for (int chm=channelMask; chm!=0; chm>>=1) if ((chm & 1)!=0) numLayers++;
    		if (externalData!=null) numLayers++;
    		int [] img= {first,second};
    		int [] channels=new int [numLayers];
    		int c=0;
    		int l=0;
    		for (int chm=channelMask; chm!=0; chm>>=1) {
    			if ((chm & 1)!=0) channels[l++]=c;
    			c++;
    		}
    		if (externalData!=null) channels[l]=3;
    		if (debugLevel>1){
    			double d=Math.sqrt(dxA*dxA+dyA*dyA);
    			System.out.println(" getShiftedSlices(): xc1="+xc1+" yc1="+yc1+" size="+size+" first="+first+" second="+second);
    			System.out.println(" getShiftedSlices(): d="+IJ.d2s(d,2)+" dX="+IJ.d2s(dxA,2)+" dY="+IJ.d2s(dyA,2));
    			System.out.println(" getShiftedSlices(): shiftImage[0]="+shiftImage[0]+" shiftImage[1]="+shiftImage[1]);
    			System.out.println(" getShiftedSlices() first image:   dX="+IJ.d2s(dx[0],2)+" dY="+IJ.d2s(dy[0],2));
    			System.out.println(" getShiftedSlices() second image:  dX="+IJ.d2s(dx[1],2)+" dY="+IJ.d2s(dy[1],2));
    			System.out.println(" getShiftedSlices() second image: iDx="+iDx+" iDyY="+iDy);
    			System.out.println(" getShiftedSlices() second image:  xc="+(xc1+iDx)+" yc="+(yc1+iDy));
    		}

    		int doubleSize=2*size;
    		int doubleLength=doubleSize*doubleSize;
    		double [] slice= new double[doubleLength];
    		double [][] result=new double [2*numLayers][];
    		double [] imgSlice;
    		if (window==null) {
    			window = new double [doubleLength];
    			double[] window1d=doubleFHT.getHamming1d(size); // Hamming
    			int index=0;
        		int halfSize=size/2;
        		int size32=3*size/2;
    			for (int iy=0;iy<doubleSize;iy++) {
    				double wy=(iy<halfSize)?window1d[iy]:((iy>size32)?window1d[iy-size]:1.0);
    				for (int ix=0;ix<doubleSize;ix++) {
        				double wx=(ix<halfSize)?window1d[ix]:((ix>size32)?window1d[ix-size]:1.0);
    					window[index++]=wx*wy;
    				}
    			}
    		}
    		if (debugLevel>1){
 //   			double d=Math.sqrt(dxA*dxA+dyA*dyA);
    			System.out.print("\ngetShiftedSlices(): channelMask="+channelMask+" channels=");
    			for (int i=0;i<channels.length;i++) System.out.print(" "+channels[i]);
    			System.out.println();

    	   		if (debugLevel>2){
    					(new ShowDoubleFloatArrays()).showArrays(
    								window,
    								doubleSize,
    								doubleSize,
    								"window");
        		}

    		}
    		for (int iImg=0;iImg<img.length;iImg++){
    			int nImg=img[iImg];
    			for (int cN=0; cN<numLayers; cN++) {
    				int chn=channels[cN];
    				if ((chn>=0) && (chn<3)){ // starts from 0 (Y)
    					imgSlice=this.overlapImages[numImgChannels*nImg+chn+1]; // skip alpha
    					if (debugLevel>2) System.out.println("Using slice "+(numImgChannels*nImg+chn+1));
    				} else {
    					imgSlice=externalData[nImg];
    					if (debugLevel>2) System.out.println("Using externalData["+nImg+"]");
    				}
    				if (debugLevel>2) System.out.println("xc["+iImg+"]="+xc[iImg]+" yc["+iImg+"]="+yc[iImg]+" shiftImage["+iImg+"]="+shiftImage[iImg]);

    				slice= getSelection(
    						imgSlice, // source image/channel slice
    						slice,
    						doubleSize, //int width,
    						doubleSize, //int height,
    						xc[iImg] , //xc,
    						yc[iImg]); //yc);
    				if (debugLevel>2) {
    					(new ShowDoubleFloatArrays()).showArrays(
								slice,
								doubleSize,
								doubleSize,
								"slice");
    				}

    				// normalize and save DC
    				double dc=0.0;
    				if (shiftImage[iImg]){
    					dc=normalizeAndWindowGetDC (slice, window); //windowInterpolation
//    					doubleFHT.shift(slice, dx[iImg], dy[iImg]);
    					doubleFHT.shift(slice, -dx[iImg], -dy[iImg]);
    				}
					int oSliceNumber=img.length*cN+iImg;
    				if (doubleSizeOutput) {
    					if (debugLevel>2) System.out.println("doubleLength="+doubleLength+"doubleSize="+doubleSize);
    					result[oSliceNumber]=new double [doubleLength];
    					int oIndex=0;
    					for (int iY=0;iY<doubleSize;iY++) {
    						for (int iX=0;iX<doubleSize;iX++){
    							if ((oIndex>=result[oSliceNumber].length) || (oIndex>=slice.length) || (oIndex>=window.length)){
    								System.out.println("iX="+iX+" iY="+iY+" oIndex="+oIndex);

    							}
//    							result[oSliceNumber][oIndex++]=imgSlice[oIndex]/window[oIndex]+dc;// no NaN with Hamming
    							result[oSliceNumber][oIndex]=slice[oIndex]/window[oIndex]+dc;// no NaN with Hamming
    							oIndex++;
    						}
    					}
    				}else{
    					result[oSliceNumber]=new double [size*size];

    					int iIndex=size*size+size/2; // top left index of the inner square size*size
    					int oIndex=0;
    					for (int iY=0;iY<size;iY++) {
    						for (int iX=0;iX<size;iX++){
    							result[oSliceNumber][oIndex++]=slice[iIndex]/window[iIndex]+dc;
    							iIndex++;
    						}
    						iIndex+=size;
    					}
    				}
    			}
    		}
    		if (debugLevel>2){
    			String [] channelNames={"Y","Cb","Cr","Lin"};
    			String [] titles=new String[2*numLayers];
    			String title="SH";
    			for (int i=0;i<img.length;i++) title+="-"+img[i];
    			double ddA=Math.sqrt(dxA*dxA+dyA*dyA);
    			title+="_D"+IJ.d2s(ddA,2)+"_X"+IJ.d2s(dxA,2)+"_Y"+IJ.d2s(dyA,2);
    			for (int i=0;i<titles.length;i++){
    				titles[i]=img[i%img.length]+channelNames[channels[i/img.length]];
    			}
    			int resultSize=doubleSizeOutput?doubleSize:size;
					(new ShowDoubleFloatArrays()).showArrays(
								result,
								resultSize,
								resultSize,
								true,
								title,
								titles);
    		}
    		return result;
    	}


    	public class DisparityTiles{
        	public double [][] disparityScales=       null; // for each channel - a pair of {scaleX, scaleY} or null if undefined (interSensor has the same)
    		public ImagePlus impDisparity=null;
    		public int corrFFTSize; // to properties
    		public int overlapStep; // to properties
    		public int paddedSize;  // same as in zTile
    		private double disparityPerEntry; //disparity increment (in pix) per array element
    		private int tilesX;
    		private int tilesY;
    		public String title;
    		private float [][] pixels=null;
    		private float [] centerPixels=null; // disparity arrays combined for the center virtual image
    		private float [][] syntheticPixels=null; // disparity arrays for individual images restored from the centerPixels
    		private BitSet innerMask=null; // will be provided to zMap instances to quickly find that there are no inner (sans padding) pixels
    		private int  []   borderMask=null; // will be provided to zMap +1:top+2:bottom+8:left+16:right (to prevent roll over when iterating around
    		private double centerPixelsFatZero=0.0;
//    		private int [][] imagePairs=null;
    		private int [][] imagePairIndices=null;
    		private double [] doubleTileWindow=null;
    		private ZTile [][][] zMap=null;
    		private Rectangle zMapWOI=null; // full: 0,0,tilesX,tilesY
    		public  Photometric photometric=null;
    		public  void setPhotometric(
    				double [][][] images,
    				int imageWidth,
    				int margins,
    				double ignoreFraction,
    				int subdivAverage,
    				int subdivHalfDifference,
    				double smoothVarianceSigma,
    				double scaleVariance,
    				int debugLevel
    		){
    			this.photometric=new Photometric(
        				images,
        				imageWidth,
        				margins,
        				ignoreFraction,
        				subdivAverage,
        				subdivHalfDifference,
        				smoothVarianceSigma,
        				scaleVariance,
        				debugLevel
        		);
    		}
    		public void showPhotometric(){
    			this.photometric.showmatchingQuality();
    		}


    		public int getPadding() {return (this.paddedSize-this.overlapStep)/2;}
    		public DisparityTiles(
    				double [][] disparityScales,
        			int corrFFTSize, // to properties
        			int overlapStep, // to properties
        			double disparityPerEntry, //disparity increment (in pix) per array element
        			String title,
        			double [][][][][] disparityTiles,
        			boolean show){
    			this.corrFFTSize=corrFFTSize;
    			this.overlapStep=overlapStep;
    			this.disparityPerEntry=disparityPerEntry;
    			this.title=title;
    			this.disparityScales=disparityScales;
    			this.impDisparity=disparityTilesToImagePlus(
            			corrFFTSize, // to properties
            			overlapStep, // to properties
            			disparityPerEntry, //disparity increment (in pix) per array element
            			title,
            			disparityTiles,
            			show);
//    			this.pixels=(float [][]) impDisparity.getStack().getImageArray();
    			this.pixels=new float [this.impDisparity.getStack().getSize()][];
    			for (int i=0;i<this.pixels.length;i++)this.pixels[i]= (float []) this.impDisparity.getStack().getPixels(i+1);
    			initDoubleTileWindow();
    		}
    		public DisparityTiles(ImagePlus imp){
    			setTilesFromImagePlus(imp);
    			this.impDisparity=imp;
//    			this.pixels=(float [][]) impDisparity.getStack().getImageArray();
    			this.pixels=new float [this.impDisparity.getStack().getSize()][];
    			for (int i=0;i<this.pixels.length;i++)this.pixels[i]= (float []) this.impDisparity.getStack().getPixels(i+1);
    			initDoubleTileWindow();
    		}

    		public DisparityTiles(String path){
    			Opener opener=new Opener();
    			ImagePlus imp=opener.openImage("", path);
    	    	if (imp==null) {
    	    		String msg="Failed to read correlation disparity file "+path;
    	    		IJ.showMessage("Error",msg);
    	    		System.out.println(msg);
    	    		throw new IllegalArgumentException (msg);
    	    	}
    			setTilesFromImagePlus(imp);
    			this.impDisparity=imp;
    			System.out.println("stack size="+this.impDisparity.getStack().getSize());
//    			this.pixels=(float [][]) impDisparity.getStack().getImageArray();
    			this.pixels=new float [this.impDisparity.getStack().getSize()][];
    			for (int i=0;i<this.pixels.length;i++)this.pixels[i]= (float []) this.impDisparity.getStack().getPixels(i+1);
    			initDoubleTileWindow();
    		}
    		public boolean isCenterPixelsDefined(){return (this.centerPixels!=null);}
    		public boolean isZMapDefined(){return (this.zMap!=null);}
        	public void deletePairsCorrelation(){this.impDisparity=null; this.pixels=null;	Runtime.getRuntime().gc(); }
        	public void deleteCenterCorrelation(){this.centerPixels=null; Runtime.getRuntime().gc(); }
        	public void deleteSyntheticCorrelation(){this.syntheticPixels=null; Runtime.getRuntime().gc(); }
        	public double getDisparityPerEntry() {return disparityPerEntry;}
        	public int    getDisparityPoints() {return this.impDisparity.getWidth()/this.tilesX;}
        	public int [][] imagePairIndices(){return this.imagePairIndices;}
        	public int [][] imagePairNumbers(){
        		int [][] pairs=new int [this.imagePairIndices.length][2];
        		for (int i=0;i<pairs.length;i++){
        			pairs[i]=this.imagePairIndices[i].clone();
        			if (pairs[i][1]>=pairs[i][0]) pairs[i][1]++;
        		}
        		return pairs;
        	}

    		//       	public void setTilesFromImagePlus(ImagePlus imp){


        	/**
        	 * encode overlapping disparity tiles to multi-slice image, save metadata
        	 * @param corrFFTSize correlation FFT size
        	 * @param overlapStep tile period
        	 * @param disparityPerEntry disparity increment per one entry in the tile array
        	 * @param title image title
        	 * @param disparityTiles data array [tileY][tileX][image number][other image index][disparity index]
        	 * @param show show image
        	 * @return image with encoded data and metadata to be saved as TIFF file
        	 */
        	public ImagePlus disparityTilesToImagePlus(
        			int corrFFTSize, // to properties
        			int overlapStep, // to properties
        			double disparityPerEntry, //disparity increment (in pix) per array element
        			String title,
        			double [][][][][] disparityTiles,
        			boolean show){ //   double [numTilesY][numTilesX][][][]; // first [tileY][tileX][image number][other image index][disparity index]
        		this.tilesY=disparityTiles.length;
        		this.tilesX=-1;
        		int numImages=0;
        		int numSecondImages=0;
        		int disparityPoints=0;
        		for (int tY=0;tY<this.tilesY;tY++) if (disparityTiles[tY]!=null){
        			if (this.tilesX<0)tilesX=disparityTiles[tY].length;
        			for (int tX=0;tX<this.tilesX;tX++) if (disparityTiles[tY][tX]!=null){
        				if (disparityTiles[tY][tX].length>numImages) numImages=disparityTiles[tY][tX].length;
        				for (int nImg=0;nImg<numImages;nImg++){
        					if ((disparityTiles[tY][tX][nImg]!=null) && (disparityTiles[tY][tX][nImg].length>numSecondImages)){
        						numSecondImages=disparityTiles[tY][tX][nImg].length;
        					}
        				}
        			}
        		}
        		boolean [][] pairsMap=new boolean[numImages][numSecondImages];
        		for (int i=0;i<numImages;i++) for (int j=0;j<numSecondImages;j++) pairsMap[i][j]=false;
        		for (int tY=0;tY<this.tilesY;tY++) if (disparityTiles[tY]!=null){
        			for (int tX=0;tX<this.tilesX;tX++) if (disparityTiles[tY][tX]!=null){
        				for (int nImg=0;nImg<disparityTiles[tY][tX].length;nImg++) if (disparityTiles[tY][tX][nImg]!=null){
        					for (int nSImg=0;nSImg<disparityTiles[tY][tX][nImg].length;nSImg++) if (disparityTiles[tY][tX][nImg][nSImg]!=null){
        						pairsMap[nImg][nSImg]=true;
        						if (disparityTiles[tY][tX][nImg][nSImg].length>disparityPoints) disparityPoints=disparityTiles[tY][tX][nImg][nSImg].length;
        					}
        				}
        			}
        		}
        		int numPairs=0;
        		for (int i=0;i<numImages;i++) for (int j=0;j<numSecondImages;j++) if (pairsMap[i][j]) numPairs++;
        		this.imagePairIndices =new int [numPairs][2]; // {firstImage, second_image_INDEX}
//        		this.imagePairs =new int [numPairs][2]; // {firstImage, second_image_INDEX}
        		int index=0;
        		for (int i=0;i<numImages;i++) for (int j=0;j<numSecondImages;j++) if (pairsMap[i][j]) {
        			this.imagePairIndices[index  ][0]=i;
        			this.imagePairIndices[index++][1]=j;
//        			this.imagePairs[index  ][0]=i;
//        			this.imagePairs[index++][1]=(j>=i)?(j+1):j;
        		}

        		int width=this.tilesX*disparityPoints;
        		int height=this.tilesY;
        		int length=width*height;
        		float [][] pixels = new float [numPairs][length];
        		for (int nPair=0;nPair<numPairs;nPair++){
        			int nImg=this.imagePairIndices[nPair][0];
        			int nSImg=this.imagePairIndices[nPair][1];
        			for (int tY=0;tY<this.tilesY;tY++){
        				if (disparityTiles[tY]!=null) {
        					for (int tX=0;tX<this.tilesX;tX++){
        						if ((disparityTiles[tY][tX]!=null) &&
        								(disparityTiles[tY][tX].length>nImg) &&
        								(disparityTiles[tY][tX][nImg]!=null)&&
        								(disparityTiles[tY][tX][nImg].length>nSImg) &&
        								(disparityTiles[tY][tX][nImg][nSImg]!=null) ) {
        							// assuming the same number of disparity points for all image pairs
        							for (int i=0;i<disparityPoints;i++) pixels[nPair][width*tY+disparityPoints*tX+i]=(float) disparityTiles[tY][tX][nImg][nSImg][i];
        						} else for (int i=0;i<disparityPoints;i++) pixels[nPair][width*tY+disparityPoints*tX+i]=Float.NaN;
        					}
        				} else for (int i=0;i<width;i++) pixels[nPair][width*tY+i]=Float.NaN;
        			}
        		}
        		ImageStack stack=new ImageStack(width,height);
        		for (int nPair=0;nPair<numPairs;nPair++){
        			int nImg=this.imagePairIndices[nPair][0];
        			int nSImg=this.imagePairIndices[nPair][1];
        			if (nSImg>=nImg) nSImg++; // here - second image number, not index!
        			stack.addSlice(String.format("%03d-%03d",nImg,nSImg), pixels[nPair]);
        		}
        		ImagePlus imp_stack = new ImagePlus(title, stack);
        		// add and encode properties
                imp_stack.setProperty("corrFFTSize",  ""+corrFFTSize);
                imp_stack.setProperty("overlapStep",  ""+overlapStep);
                imp_stack.setProperty("disparityPerEntry",  ""+disparityPerEntry);

                imp_stack.setProperty("tilesY",  ""+this.tilesY);
                imp_stack.setProperty("tilesX",  ""+this.tilesX);

                imp_stack.setProperty("disparityPoints",  ""+disparityPoints);
                imp_stack.setProperty("numPairs",  ""+numPairs);
                for (int nPair=0;nPair<numPairs;nPair++){
        			int nImg=this.imagePairIndices[nPair][0];
        			int nSImg=this.imagePairIndices[nPair][1];
        			if (nSImg>=nImg) nSImg++; // here - second image number, not index!
                    imp_stack.setProperty("firstImage_"+nPair,  ""+nImg);
                    imp_stack.setProperty("secondImage_"+nPair,  ""+nSImg);
                }
      			if (this.disparityScales!=null) {
                    imp_stack.setProperty("numberOfImages", this.disparityScales.length+"");
                    for (int i=0;i<this.disparityScales.length;i++){
                    	imp_stack.setProperty("disparityScalesX_"+i, this.disparityScales[i][0]+"");
                    	imp_stack.setProperty("disparityScalesY_"+i, this.disparityScales[i][1]+"");
                    }
      			}
               	(new JP46_Reader_camera(false)).encodeProperiesToInfo(imp_stack);
        		imp_stack.getProcessor().resetMinAndMax();
        		if (show) imp_stack.show();
        		return imp_stack;
        	}

        	public void setTilesFromImagePlus(
        			ImagePlus imp){
                String [] requiredProperties={
                		"corrFFTSize",
                		"overlapStep",
                		"disparityPerEntry",
                		"tilesY",
                		"tilesX",
                		"numPairs"};
            	(new JP46_Reader_camera(false)).decodeProperiesFromInfo(imp);
                for (int i=0; i<requiredProperties.length;i++) if (imp.getProperty(requiredProperties[i])==null){
            		String msg="Required property "+requiredProperties[i]+" is not defined in "+imp.getTitle();
            		IJ.showMessage("Error",msg);
            		throw new IllegalArgumentException (msg);
                }
                int numPairs=      Integer.parseInt  ((String) imp.getProperty("numPairs"));
                int [][] pairs=new int [numPairs][2];
                String s,msg;
                for (int nPair=0;nPair<numPairs;nPair++){
                	s= (String) imp.getProperty("firstImage_"+nPair);
                	if (s== null) {
                		msg="Property "+"firstImage_"+nPair+" is not defined in "+imp.getTitle();
                		IJ.showMessage("Error",msg);
                		throw new IllegalArgumentException (msg);
                	}
                	pairs[nPair][0]=Integer.parseInt(s);
                	s= (String) imp.getProperty("secondImage_"+nPair);
                	if (s== null) {
                		msg="Property "+"secondImage_"+nPair+" is not defined in "+imp.getTitle();
                		IJ.showMessage("Error",msg);
                		throw new IllegalArgumentException (msg);
                	}
                	pairs[nPair][1]=Integer.parseInt(s);
                }
        		this.imagePairIndices =new int [numPairs][2]; // {firstImage, second_image_INDEX}
                for (int nPair=0;nPair<numPairs;nPair++){
                	this.imagePairIndices[nPair][0]=pairs[nPair][0];
                	this.imagePairIndices[nPair][1]=(pairs[nPair][1]>pairs[nPair][0])?(pairs[nPair][1]-1):pairs[nPair][1];
                }
                this.corrFFTSize=  Integer.parseInt  ((String) imp.getProperty("corrFFTSize"));
                this.overlapStep=  Integer.parseInt  ((String) imp.getProperty("overlapStep"));
                this.disparityPerEntry=  Double.parseDouble ((String) imp.getProperty("disparityPerEntry"));
                this.tilesY=  Integer.parseInt  ((String) imp.getProperty("tilesY"));
                this.tilesX=  Integer.parseInt  ((String) imp.getProperty("tilesX"));
                if (imp.getProperty("numberOfImages")==null){ // later - add to required properties, now just to use older files
                	int numImages=-1;
                	for (int i=0;i<pairs.length;i++){
                		if (pairs[i][0]>numImages) numImages=pairs[i][0];
                		if (pairs[i][1]>numImages) numImages=pairs[i][1];
                	}
                	numImages++;
    				double [][] disparityPair={{-0.5,0.0},{0.5,0.0}};
    				double [][] disparityTriplet={{0.0,Math.sqrt(3)/3}, {-0.5,-Math.sqrt(3)/6},{0.5,-Math.sqrt(3)/6}};
    				if ((numImages==2)|| (numImages==3)){
    					this.disparityScales=new double [numImages][2];
    					for (int i=0;i<numImages;i++) this.disparityScales[i]=((numImages==3)?disparityTriplet[i]:disparityPair[i]).clone();
    				} else {
                		msg="Only doublets and triplets are supported with default disparityScales - you seem to have "+numImages+
                		" images. Please generate files with disparityScales data included.";
                		IJ.showMessage("Error",msg);
                		throw new IllegalArgumentException (msg);
    				}
                } else {
                	this.disparityScales=new double [Integer.parseInt  ((String) imp.getProperty("numberOfImages"))][2];
                	for (int i=0;i<this.disparityScales.length;i++){
                    	this.disparityScales[i][0]=Double.parseDouble ((String) imp.getProperty("disparityScalesX_"+i));
                    	this.disparityScales[i][1]=Double.parseDouble ((String) imp.getProperty("disparityScalesY_"+i));
                	}
                }
        	}


           	public void combineDisparityToCenter(
        			final double fatZero,
        			final int threadsMax,
        			final boolean showProgress,
        			final int debugLevel){
        		if ((this.centerPixels!=null) && (fatZero==this.centerPixelsFatZero)) return; // already done
        		this.centerPixels=new float [this.pixels[0].length]; // all should be the same
        		// TODO - convert to multithread?
        		final int tiles=this.tilesY*this.tilesX;
        		final int tilesX=this.tilesX;
        		final double [][] disparityScales=this.disparityScales;
        		final int [][] imagePairIndices=this.imagePairIndices;
        		final Thread[] threads = newThreadArray(threadsMax);
        		final AtomicInteger tileIndexAtomic = new AtomicInteger(0);
        		final float [] centerPixels=this.centerPixels;
        		final float [][] pixels=this.pixels;
		   		if (showProgress) IJ.showProgress(0.0);
		   		if (showProgress) IJ.showStatus("Centering correlation data...");
		   		for (int ithread = 0; ithread < threads.length; ithread++) {
		   			threads[ithread] = new Thread() {
		   				@Override
						public void run() {
		   					// common for the tread
		   					double [][] partial=new double [imagePairIndices.length][];
           					for (int tile=tileIndexAtomic.getAndIncrement(); tile<tiles;tile=tileIndexAtomic.getAndIncrement()){
           						boolean debugThis= (debugLevel>2 ) && (tile== ( (tiles/tilesX/2)*tilesX +(tilesX/2)));
           						if (debugThis){
           							System.out.println("tile="+tile+" tiles="+tiles+" tilesX="+tilesX);
           						}

           						for (int nPair=0;nPair<partial.length;nPair++){
           							int nImg= imagePairIndices[nPair][0];
           							double [] dXY=disparityScales[nImg];
               						if (debugThis){
               							System.out.println("\n\nnImg="+nImg+" nPair="+nPair+" dXY[0]="+dXY[0]+" dXY[1]="+dXY[1]);
               						}
           							partial [nPair]=getRecenteredDisparity(
           									tile%tilesX, //int tileX,
           									tile/tilesX, //int tileY,
           									pixels[nPair],
           				        			dXY,
           				        			debugThis?(debugLevel+2):debugLevel);
           							if (partial [nPair]==null) System.out.println("combineDisparityToCenter(): partial["+nPair+"]==null, tile="+tile);
           						}
           						int numNotNull=0;
           						for (int i=0;i<partial.length;i++) if (partial[i]!=null) numNotNull++;
           						double pwr=1.0/numNotNull;
           						int disparityPoints=getDisparityPoints();
           						for (int i=0;i<disparityPoints;i++) {
           							double d=1.0;
           							for (int n=0;n<partial.length;n++) if (partial[n]!=null) d*=(partial[n][i]+fatZero);
           							partial[0][i]=((d>=0)?Math.pow(d,pwr):0.0)-fatZero;
           						}
           						if (debugThis) {
           							for (int i=0;i<partial[0].length;i++){
           								System.out.println(i+" "+partial[0][i]); // center tile

           							}
           						}

           				       	setTileDisparity(
           			        			partial[0],
           			        			tile%tilesX, //int itx,
           			        			tile/tilesX,
           			        			centerPixels); //int ity

           		   	   			if (showProgress){
           		   	   				final int fTile=tile;
           	    					SwingUtilities.invokeLater(new Runnable() {
           	    						@Override
										public void run() {
           	    							// Here, we can safely update the GUI
           	    							// because we'll be called from the
           	    							// event dispatch thread
           	    							IJ.showProgress(fTile,tiles);
           	    						}
           	    					});
           		   	   			}
           					}
		   				}
		   			};
		   		}
		   		startAndJoin(threads);
		   		if (showProgress) IJ.showStatus("");
		   		if (showProgress) IJ.showProgress(1.0);
        		this.centerPixelsFatZero=fatZero;
        		this.syntheticPixels=null; // invalidate
        		if (debugLevel>2){
        			(new ShowDoubleFloatArrays()).showArrays(
        					this.centerPixels,
        					this.impDisparity.getWidth(),
        					this.centerPixels.length/this.impDisparity.getWidth(),
        					"center-corr");
        		}
        	}


        	/**
        	 * Calculating combined correlation arrays for the actual cameras from the  center camera
        	 * Will do nothing if the data is already calcualted
        	 * @param threadsMax maximal number of therads to use
        	 * @param showProgress show ImageJ progress/status
        	 * @param debugLevel debug level
        	 */
        	public void moveDisparityFromCenter(
        			final int threadsMax,
        			final boolean showProgress,
        			final int debugLevel){
        		if (this.syntheticPixels!=null) return; // already done
        		if (this.centerPixels==null) {
    				String msg="Centered disparity data is not defined";
    				IJ.showMessage("Error",msg);
    				throw new IllegalArgumentException (msg);
        		}
        		final int numImages= this.disparityScales.length;
        		this.syntheticPixels=new float [numImages] [this.centerPixels.length]; // all should be the same
        		// TODO - convert to multithread?
        		final int tilesPerImage=this.tilesY*this.tilesX;
        		final int tiles=tilesPerImage*numImages;
        		final int tilesX=this.tilesX;
        		final double [][] disparityScales=this.disparityScales;
        		final Thread[] threads = newThreadArray(threadsMax);
        		final AtomicInteger tileIndexAtomic = new AtomicInteger(0);
        		final float [] centerPixels=this.centerPixels;
        		final float [][] syntheticPixels=this.syntheticPixels;
		   		if (showProgress) IJ.showProgress(0.0);
		   		if (showProgress) IJ.showStatus("Re-centering correlation data...");
		   		for (int ithread = 0; ithread < threads.length; ithread++) {
		   			threads[ithread] = new Thread() {
		   				@Override
						public void run() {
           					for (int tile=tileIndexAtomic.getAndIncrement(); tile<tiles;tile=tileIndexAtomic.getAndIncrement()){
           						int nImage=tile/tilesPerImage;
           						int tileX=tile%tilesX;
           						int tileY=(tile%tilesPerImage)/tilesX;

           						boolean debugThis= (debugLevel>2 ) && (nImage==0)&& (tile== ( (tiles/tilesX/2)*tilesX +(tilesX/2)));
           						if (debugThis){
           							System.out.println("nImage="+nImage+"tile="+tile+" tiles="+tiles+" tilesX="+tilesX);
           						}
           						double [] dXY={-disparityScales[nImage][0],-disparityScales[nImage][1]};
           						if (debugThis){
           							System.out.println("\n\nnImg="+nImage+" dXY[0]="+dXY[0]+" dXY[1]="+dXY[1]);
           						}
           						double [] tileData=getRecenteredDisparity(
           								tileX, //int tileX,
           								tileY, //int tileY,
           								centerPixels,
           								dXY,
           								debugThis?(debugLevel+2):debugLevel);
           						if (tileData==null) System.out.println("combineDisparityToCenter(): tileData==null, tile="+tile);
           						setTileDisparity(
           								tileData,
           								tileX,
           								tileY,
           								syntheticPixels[nImage]);
           		   	   			if (showProgress){
           		   	   				final int fTile=tile;
           	    					SwingUtilities.invokeLater(new Runnable() {
           	    						@Override
										public void run() {
           	    							// Here, we can safely update the GUI
           	    							// because we'll be called from the
           	    							// event dispatch thread
           	    							IJ.showProgress(fTile,tiles);
           	    						}
           	    					});
           		   	   			}
           					}
		   				}
		   			};
		   		}
		   		startAndJoin(threads);
		   		if (showProgress) IJ.showStatus("");
		   		if (showProgress) IJ.showProgress(1.0);
        	}

        	private void  initDoubleTileWindow(){
        		int size=this.overlapStep*2;
        		this.doubleTileWindow = new double [size*size];
    			double [] window1d=(new DoubleFHT()).getHamming1d(size,0); // pure cosine
    			for (int i=0;i<this.doubleTileWindow.length;i++) this.doubleTileWindow[i]=window1d[i/size]*window1d[i%size];

        	}
/*
        	public double [] getRecenteredDisparity(
        			int itx,
        			int ity,
        			double [] relDxy, // disparityScales new - disparityScales old
        			int debugLevel
        			){
        		return getRecenteredDisparity(itx,ity,this.centerPixels,relDxy,debugLevel);}
*/
        	/**
        	 * Calculate correlation tile for different center
        	 * @param tileX horizontal tile number
        	 * @param tileY vertical tile number
        	 * @param tiledCorrelation encoded tiled correlation ararys
        	 * @param relDxy differenctial disparity vector - x and y shift of the pair first image relative to the new center
        	 * @param debugLevel debug level
        	 * @return array of correlatin vs disparity for a virtual camera view
        	 */
        	public double [] getRecenteredDisparity(
        			int tileX,
        			int tileY,
        			float [] tiledCorrelation,
        			double [] relDxy,
        			int debugLevel
        			){
        		if (debugLevel>3) System.out.println("getRecenteredDisparity("+tileX+","+tileY+",pixels,...)");
        		int disparityPoints=getDisparityPoints();
        		double disparityPerEntry=getDisparityPerEntry();
        		double lDisp=Math.sqrt(relDxy[0]*relDxy[0]+relDxy[1]*relDxy[1])*(disparityPoints-1)*disparityPerEntry; // in Pixels
        		int maxTileDis=(int) Math.floor(lDisp/this.overlapStep)+2; // 1 is enough, maybe some rounding?
        		double [][][] tileCache=new double [2*maxTileDis+1][2*maxTileDis+1][];
        		for (int i=0;i<tileCache.length;i++) for (int j=0;j<tileCache[0].length;j++) tileCache[i][j]=null;
        		double [] result=new double [disparityPoints];
        		double [] kTXY={relDxy[0]*disparityPerEntry/this.overlapStep,relDxy[1]*disparityPerEntry/this.overlapStep};
        		int windowSize=2*this.overlapStep;
        		if (debugLevel>3) System.out.print(" maxTileDis="+maxTileDis);

        		for (int i=0;i<result.length;i++){
            		if (debugLevel>3) System.out.print("getRecenteredDisparity() i="+i);

        			double dTX=kTXY[0]*i; // tile period is 0.5 in this units, all tiles within -1.0...+1.0 may contribute - normally 4
        			double dTY=kTXY[1]*i;
        			double d=0.0;
        			double weight=0.0;
        			int idTX0=(int) Math.floor(dTX);
        			int idTY0=(int) Math.floor(dTY);
        			for (int iTY=0;iTY<2;iTY++) {
        				int iTYCache=idTY0+iTY+maxTileDis;
        				int pY= (int) Math.round ((dTY-(idTY0+iTY)+1.0)*this.overlapStep); // 0... 2*this.overlapStep
						int srcTY=tileY+idTY0+iTY;
        				if ((srcTY>=0) && (srcTY<this.tilesY) && (pY>=0) && (pY<windowSize)) for (int iTX=0;iTX<2;iTX++){
        					int iTXCache=idTX0+iTX+maxTileDis;
        					int pX=(int) Math.round ((dTX-(idTX0+iTX)+1.0)*this.overlapStep);  // 0... 2*this.overlapStep
							int srcTX=tileX+idTX0+iTX;
        					if ((srcTX>=0) && (srcTX<this.tilesX) && (pX>=0) && (pX<windowSize)){
        						if (tileCache[iTYCache][iTXCache]==null){
        							if ((debugLevel>3) || ((tileX+idTX0+iTX)<0) || ((tileY+idTY0+iTY)<0)){
        								System.out.println ("\n>>>getting tile to cache:  iTXCache="+iTXCache+" iTYCache="+iTYCache+" tileX="+tileX+" tileY="+tileY+" idTX0="+idTX0+" idTY0="+idTY0+ " iTX="+iTX+" iTY="+iTY);
        							}
        							tileCache[iTYCache][iTXCache]=getTileDisparity( // negative?
        									srcTX,
        									srcTY,
        									tiledCorrelation);
        						}
        						double w=this.doubleTileWindow[pY*windowSize+pX];
        						weight+=w;
        						d+=w*tileCache[iTYCache][iTXCache][i];
        		        		if (debugLevel>3) {
        		        			System.out.print(" pX="+pX+" pY="+pY+" srcTX="+srcTX+" srcTY="+srcTY+" d="+tileCache[iTYCache][iTXCache][i]+" w="+w);
        		        		}

        					}
        				}
        			}
        			result[i]=(weight>0)?(d/weight):0.0;
            		if (debugLevel>3) System.out.println(" weight="+weight+" result["+i+"]="+result[i]);

        		}
        		if (debugLevel>3) {
        			for (int i=0;i<result.length;i++){
        				System.out.println(i+" "+result[i]+" "+tileCache[maxTileDis][maxTileDis][i]); // center tile

        			}
        		}
         		return result;
        	}

        	public double [] getTileDisparity(
        			int itx,
        			int ity,
        			float [] data){ //npair <0 - use center data
        		int width=this.impDisparity.getWidth();
        		int disparityPoints=getDisparityPoints();
        		double [] result=new double [disparityPoints];
        		int index=ity*width+itx*disparityPoints;
        		for (int i=0;i<disparityPoints;i++) result[i]=data[index++]; //oob -146583
        		return result;
        	}

        	public void setTileDisparity(
        			double []disparityArray,
        			int itx,
        			int ity,
        			float [] data){ //npair <0 - use center data
        		int width=this.impDisparity.getWidth();
        		int disparityPoints=getDisparityPoints();
        		int index=ity*width+itx*disparityPoints;
        		for (int i=0;i<disparityPoints;i++) {
        			if ((index>=data.length) || (i>=disparityArray.length)){
        				System.out.println("setTileDisparity(disparityArray,"+itx+","+ity+", data): index="+index+
        						" data.length="+data.length+" i="+i+" disparityArray.length="+disparityArray.length+" disparityPoints="+disparityPoints+" width="+width);
        			}
        			data[index++]=(float) disparityArray[i]; //oob 20375037
        		}
        	}

        	public double [] getCorrelationArray(
        			double x,
        			double y,
        			double fatZero,
        			int threadsMax,
        			boolean showProgress,
        			int debugLevel ){
    			combineDisparityToCenter(
            			fatZero,
            			threadsMax,
            			showProgress,
            			debugLevel);
        		return getCorrelationArray(
        				this.centerPixels,
            			x,
            			y,
            			threadsMax,
            			showProgress,
            			debugLevel );
        	}

        	public double [] getCorrelationArray(
        			int imageNumber,
        			double x,
        			double y,
        			int threadsMax,
        			boolean showProgress,
        			int debugLevel ){
        		moveDisparityFromCenter(
            			threadsMax,
            			showProgress,
            			debugLevel);
        		float [] data=this.syntheticPixels[imageNumber];
        		return  getCorrelationArray(
            			data,
            			x,
            			y,
            			threadsMax,
            			showProgress,
            			debugLevel );
        	}


        	public double [] getCorrelationArray(
        			int firstImage,
        			int secondImageIndex,
//        			double fatZero,
        			double x,
        			double y,
        			int threadsMax,
        			boolean showProgress,
        			int debugLevel ){
        		int nPair;
        			nPair= imagesToNPair(firstImage, secondImageIndex);
        			if (nPair<0) {
        				System.out.println("getCorrelationArray(): No data for pair of images "+firstImage+" and second image index "+secondImageIndex+
        						" (second image = "+((secondImageIndex>=firstImage)?(secondImageIndex+1):(secondImageIndex))+")");
        				return null;
        			}
        		float [] data=(nPair>=0)?this.pixels[nPair]:this.centerPixels;

        		return  getCorrelationArray(
            			data,
            			x,
            			y,
            			threadsMax,
            			showProgress,
            			debugLevel );
        	}

        	private double [] getCorrelationArray(
        			float [] data,
        			double x,
        			double y,
        			int threadsMax,
        			boolean showProgress,
        			int debugLevel ){
        		double tx=x/this.overlapStep;
        		double ty=y/this.overlapStep;
        		int itx0=(int) Math.floor(tx);
        		int ity0=(int) Math.floor(ty);
        		tx-=itx0;
        		ty-=ity0;
        		int itx1=itx0+1;
        		int ity1=ity0+1;
        		if ((itx1<0) || (itx0>=this.tilesX) || (ity1<0) || (ity0>=this.tilesY)) {
        			System.out.println("getCorrelationArray(): this.overlapStep="+this.overlapStep+" this.tilesX="+this.tilesX+"this.tilesY="+this.tilesY);
        			System.out.println("getCorrelationArray(): itx0="+itx0+" ity0="+ity0+" out of bounds. x="+x+" y="+y);
        			return null;
        		}
        		if (itx0<0){itx0=itx1;tx=0;	}
        		if (itx1>=this.tilesX){itx1=itx0;tx=0;}
        		if (ity0<0){ity0=ity1;ty=0;	}
        		if (ity1>=this.tilesY){ity1=ity0;ty=0;}
        		int width=this.impDisparity.getWidth();
        		int disparityPoints=width/this.tilesX;
        		int index00=ity0*width+itx0*disparityPoints;
        		int index10=ity0*width+itx1*disparityPoints;
        		int index01=ity1*width+itx0*disparityPoints;
        		int index11=ity1*width+itx1*disparityPoints;
        		double [] result=new double [disparityPoints];
        		for (int i=0;i<disparityPoints;i++){
        			if ((i>result.length) ||
        					((index00+i)>=data.length) || ((index10+i)>=data.length) || ((index01+i)>=data.length) || ((index11+i)>=data.length)){
        				System.out.println("width="+width+" index00="+index00+" index10="+index10+" index01="+index01+" index11="+index11);
        				System.out.println("itx0="+itx0+" itx1="+itx1+" ity0="+ity0+" ity1="+ity1+" this.tilesX="+this.tilesX+" this.tilesY="+this.tilesY);

        			}
        			result[i]=  // bi-linear interpolation //TODO:java.lang.ArrayIndexOutOfBoundsException: 24760512at PixelMapping$InterSensor$DisparityTiles.getCorrelationArray(PixelMapping.java:4619)

        				(data[index00+i]*(1.0-tx)+data[index10+i]*tx)*(1.0-ty)+
        				(data[index01+i]*(1.0-tx)+data[index11+i]*tx)*ty;
        		}
        		return result;
        	}


        	public double [] getCorrelationArray(
        			int imageNumber,
        			double fatZero,
        			double x,
        			double y,
        			int threadsMax,
        			boolean showProgress,
        			int debugLevel ){
        		int numSecondImages=0;
        		for (int i=0;i<this.imagePairIndices.length;i++) if (this.imagePairIndices[i][0]==imageNumber) numSecondImages++;
        		if (numSecondImages==0){
        			System.out.println("getCorrelationArray(): No pairs for image "+imageNumber);
        			return null;
        		}
        		int [] secondImages=new int[numSecondImages];
        		int index=0;
        		for (int i=0;i<this.imagePairIndices.length;i++) if (this.imagePairIndices[i][0]==imageNumber) secondImages[index++]=this.imagePairIndices[i][1];
        		double [][] partial=new double[secondImages.length][];
        		for (int i=0;i<partial.length;i++) partial[i]=getCorrelationArray(
            			imageNumber,
            			secondImages[i],
            			x,
            			y,
            			threadsMax,
            			showProgress,
            			debugLevel );

        		double [] result=new double [partial[0].length];
        		double pwr=1.0/partial.length;  // sqrt for 2 other images
        		for (int nP=0;nP<result.length;nP++){
    				double d=1.0;
        			for (int nI=0;nI<partial.length;nI++){
        				if (partial[nI][nP]+fatZero>0) d*=(partial[nI][nP]+fatZero);
        				else d=0;
        			}
    				result[nP]=Math.pow(d, pwr)-fatZero;
        		}
        		return result;
        	}
        	private int imagesToNPair(int firstImage, int secondImageIndex){
        		for (int i=0;i<this.imagePairIndices.length;i++) if ((this.imagePairIndices[i][0]==firstImage)&& (this.imagePairIndices[i][1]==secondImageIndex)) return i;
        		return -1;
        	}
        	public void initZMap(){
        		for (int i=0;i<this.disparityScales.length;i++) initZMap(i);
        	}
        	public void initZMap(int nImg){
        		if (this.zMap==null) {
        			this.zMap=new ZTile[this.disparityScales.length][][];
        			for (int i=0;i<this.zMap.length;i++) this.zMap[i]=null;
        		}
        		if (this.zMap[nImg]==null){
        			this.zMap[nImg]=new ZTile[this.tilesY][this.tilesX];
        			for (int tileY=0;tileY<this.tilesY;tileY++) for (int tileX=0;tileX<this.tilesX;tileX++){
        				this.zMap[nImg][tileY][tileX]=null;
        			}
        		}
        	}
        	public Rectangle pixelsToTilesWOI(
        			Rectangle pixelWOI){
        		Rectangle woi=new Rectangle (
        				pixelWOI.x/this.overlapStep,
        				pixelWOI.y/this.overlapStep,
        				(pixelWOI.x+pixelWOI.width -1)/this.overlapStep - pixelWOI.x/this.overlapStep +1,
        				(pixelWOI.y+pixelWOI.height-1)/this.overlapStep - pixelWOI.y/this.overlapStep +1);
        		return woi;
        	}

        	public void setupZMap(
        			int maxNumber,
        			double minFirst,
        			double minAbsolute,
        			double minRelative,
        			double mergeMax,
        			int overlap,
        			double zMapMinForeground,
        			int threadsMax,
        			boolean showProgress,
        			int debugLevel){
        		for (int i=0;i<this.disparityScales.length;i++){
        			setupZMap(null, i, maxNumber, minFirst, minAbsolute, minRelative, mergeMax,overlap,zMapMinForeground, threadsMax,showProgress,debugLevel);
        		}
        	}

        	public void setupZMap(
        			int nImg,
        			int maxNumber,
        			double minFirst,
        			double minAbsolute,
        			double minRelative,
        			double mergeMax,
        			int overlap,
        			double zMapMinForeground,
        			int threadsMax,
        			boolean showProgress,
        			int debugLevel){
        		setupZMap(null, nImg, maxNumber, minFirst, minAbsolute, minRelative,mergeMax,overlap,zMapMinForeground,threadsMax,showProgress,debugLevel);
        	}


        	public void setupZMap(
        			Rectangle woi, // in tiles - may be
        			int maxNumber,
        			double minFirst,
        			double minAbsolute,
        			double minRelative,
        			double mergeMax,
        			final int overlap,
        			final double zMapMinForeground,
        			int threadsMax,
        			boolean showProgress,
        			int debugLevel){
        		for (int i=0;i<this.disparityScales.length;i++){
        			setupZMap(woi, i, maxNumber, minFirst, minAbsolute, minRelative,mergeMax,overlap,zMapMinForeground,threadsMax,showProgress,debugLevel);
        		}
        	}
           	public int fillForegroundGapsStep(
        			Rectangle woi, // in tiles - may be
        			int minNeib,
        			double maxDifference,
        			double foregroundThreshold,
        			int threadsMax,
        			boolean showProgress,
        			int debugLevel){
           		int totalFilled=0;
        		for (int i=0;i<this.disparityScales.length;i++){
        			totalFilled+=fillForegroundGapsStep(woi, i, minNeib, maxDifference, foregroundThreshold, threadsMax, showProgress,debugLevel);
        		}
        		return totalFilled;
           	}


           	public int fillForegroundGapsStep(
        			Rectangle woi, // in tiles - may be
        			final int nImg,
        			final int minNeib,
        			final double maxDifference,
        			final double foregroundThreshold,
        			final int threadsMax,
        			final boolean showProgress,
        			final int debugLevel){
        		if (woi==null) woi=new Rectangle (0, 0, this.tilesX,this.tilesY);
        		this.zMapWOI=new Rectangle();
        		this.zMapWOI.x=(woi.x>=0)?woi.x:0;
        		this.zMapWOI.y=(woi.y>=0)?woi.y:0;
        		this.zMapWOI.width= (((woi.x+woi.width) <=this.tilesX)?(woi.x+woi.width): this.tilesX)-this.zMapWOI.x;
        		this.zMapWOI.height=(((woi.y+woi.height)<=this.tilesY)?(woi.y+woi.height):this.tilesY)-this.zMapWOI.y;
        		final Rectangle zMapWOI=this.zMapWOI;
        		if ((this.zMap==null) || (this.zMap[nImg]==null)) initZMap(nImg);
        		final ZTile [][] thisZMap=this.zMap[nImg];
        		final int tilesX=this.zMapWOI.width;
        		final int tiles=this.zMapWOI.width*this.zMapWOI.height;
        		if (debugLevel>2) System.out.println("fillForegroundGaps() woi.x="+woi.x+" woi.y="+woi.y+
        				" woi.width="+woi.width+" woi.height="+woi.height);
        		if (debugLevel>2) System.out.println("fillForegroundGaps() zMapWOI.x="+zMapWOI.x+" zMapWOI.y="+zMapWOI.y+
        				" zMapWOI.width="+zMapWOI.width+" zMapWOI.height="+zMapWOI.height);

        		if (debugLevel>1) System.out.println("fillForegroundGapsStep() foregroundThreshold="+foregroundThreshold+" maxDifference="+maxDifference+" tiles="+tiles);

        		final Thread[] threads = newThreadArray(threadsMax);
		   		final AtomicInteger tileIndexAtomic     = new AtomicInteger(0);
		   		if (showProgress) IJ.showProgress(0.0);
		   		if (showProgress) IJ.showStatus("Setting up zMap for image "+nImg+" ("+(nImg+1)+" of "+this.disparityScales.length+") ...");
		   		final AtomicInteger numberUpdated = new AtomicInteger(0);
		   		for (int ithread = 0; ithread < threads.length; ithread++) {
		   			threads[ithread] = new Thread() {
		   				@Override
						public void run() {
		   					for (int tile=tileIndexAtomic.getAndIncrement(); tile<tiles;tile=tileIndexAtomic.getAndIncrement()){
		   						if (fillForegroundTile(
		   								zMapWOI.x+tile%tilesX, //int tileX,
		   								zMapWOI.y+tile/tilesX, //int tileY,
		   								thisZMap,
		   			        			minNeib,
		   			        			maxDifference,
		   			        			foregroundThreshold,
		   			        			debugLevel)) numberUpdated.getAndIncrement();
		   						if (showProgress){
		   							final int finalTile=tile;
		   							SwingUtilities.invokeLater(new Runnable() {
		   								@Override
										public void run() {
		   									IJ.showProgress(finalTile,tiles);
		   								}
		   							});
		   						}
		   					}
		   				}
		   			};
		   		}
		   		startAndJoin(threads);
		   		IJ.showProgress(1.0);
		   		return numberUpdated.get();
        	}
        	private boolean fillForegroundTile(
        			int tileX,
        			int tileY,
        			ZTile [][] zMap,
        			int minNeib,
        			double maxDifference,
        			double foregroundThreshold,
		        	int 	debugLevel){
        		ZTile zTile=zMap[tileY][tileX];
        		if (zTile==null) return false;
        		if (debugLevel>4) System.out.println("0:fillForegroundTile("+tileX+","+tileY+",zMap,"+minNeib+","+maxDifference+","+foregroundThreshold+","+debugLevel+"), zTile.foregroundIndex="+zTile.foregroundIndex);
        		if (zTile.foregroundIndex==0) return false; // already front
        		if (debugLevel>3) System.out.println("fillForegroundTile("+tileX+","+tileY+",zMap,"+minNeib+","+maxDifference+","+foregroundThreshold+","+debugLevel+"), zTile.foregroundIndex="+zTile.foregroundIndex);
        		int [][] dirs8={{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}};
        		for (int plane=0; plane <zTile.foregroundIndex;plane++) if (zTile.maximums[plane][1]>=foregroundThreshold){
        			// See if there are enough neighbors where current foreground is close enough to this one
        			int numMatched=0;
        			for (int iDir=0;iDir<dirs8.length;iDir++){
        				int tileX1=tileX+dirs8[iDir][0];
        				int tileY1=tileY+dirs8[iDir][1];
//        				if (tileY>=zMap.length){
//        					System.out.print("fillForegroundTile("+tileX+","+tileX+",...) zMap.length="+zMap.length+" zMap[0].length="+zMap[0].length);
//       				}
        				if ((tileY1>=0) && (tileY1<zMap.length) && (tileX1>=0) && (tileX1<zMap[0].length) && // other tile within limits
        						(zMap[tileY1][tileX1]!=null) &&                                              // and exists
 //       						(zMap[tileY1][tileX1].foregroundIndex<zMap[tileY][tileX].maximums.length) &&
        						(zMap[tileY1][tileX1].foregroundIndex<zMap[tileY1][tileX1].maximums.length) && // not at infinity
        						(Math.abs(zMap[tileY1][tileX1].maximums[zMap[tileY1][tileX1].foregroundIndex][0]-zTile.maximums[plane][0])<=maxDifference)){
        					numMatched++;
        				}
        			}
        			if (numMatched>=minNeib){
        				zTile.setForegroundPlane(plane);
        				return true;
        			}
        		}
        		return false;

        	}



    		public void foregroundByOcclusion(
    				double [][][] imageData, // [img]{Alpha, Y,Cb,Cr, Ext}. Alpha may be null - will handle later?
    				int imageFullWidth,
    				double bgFraction, // process maximal correlation or the first farther plane with correlation intensity >= this part of maximal
    				double blurVarianceSigma,
    				int refineTilePeriod,
    				double refinePhaseCoeff,
    				double refineHighPassSigma,
    				double refineLowPassSigma,
    				double refineCorrMaxDistance,
    				double refineCorrThreshold,
    				int refineSubPixel,
    				double zMapMinForeground,
    				int zMapVarMask,
    				double [] zMapVarThresholds,
    				double [] zMapVarWeights,
    				int     normalizeCorrelation,
    				int zMapCorrMask,
    				double [] zMapCorrThresholds,
    				double [] zMapCorrWeights,
    				int debugRow,
    				int debugColumn,
    				int       threadsMax,
    				boolean   showProgress,
    				int       debugLevel){
    			// process all image pairs
    			for (int nImg=0;nImg<this.disparityScales.length;nImg++){
    				int [] sImg=new int [this.disparityScales.length-1];
    				for (int j=0;j<sImg.length;j++)sImg[j]=(j<nImg)?j:(j+1);
    				foregroundByOcclusion(
    	    				nImg,
    	    				sImg, // list of second images
    	    				imageData, // [img]{Alpha, Y,Cb,Cr, Ext}. Alpha may be null - will handle later?
    	    				imageFullWidth,
    	    				bgFraction, // process maximal correlation or the first farther plane with correlation intensity >= this part of maximal
    	    				blurVarianceSigma,
    	    				refineTilePeriod,
    	    				refinePhaseCoeff,
    	    				refineHighPassSigma,
    	    				refineLowPassSigma,
    	    				refineCorrMaxDistance,
    	    				refineCorrThreshold,
    	    				refineSubPixel,
    	    				zMapMinForeground,
    	    				zMapVarMask,
    	    				zMapVarThresholds,
    	    				zMapVarWeights,
    	    				normalizeCorrelation,
    	    				zMapCorrMask,
    	    				zMapCorrThresholds,
    	    				zMapCorrWeights,
    	    				debugRow,
    	    				debugColumn,
    	    				threadsMax,
    	    				showProgress,
    	    				debugLevel);
    			}
    		}

        	public void foregroundByOcclusion( // remove unneeded parameters later
					final int nImg,
					final int[] sImg, // list of second images
    				final double [][][] imageData, // [img]{Alpha, Y,Cb,Cr, Ext}. Alpha may be null - will handle later?
    				final int imageFullWidth,
    				final double bgFraction, // process maximal correlation or the first farther plane with correlation intensity >= this part of maximal
    				final double blurVarianceSigma,
    				final int refineTilePeriod,
    				final double refinePhaseCoeff,
    				final double refineHighPassSigma,
    				final double refineLowPassSigma,
    				final double refineCorrMaxDistance,
    				final double refineCorrThreshold,
    				final int refineSubPixel,
    				final double zMapMinForeground,
    				final int zMapVarMask,
    				final double [] zMapVarThresholds,
    				final double [] zMapVarWeights,
    				final int     normalizeCorrelation,
    				final int zMapCorrMask,
    				final double [] zMapCorrThresholds,
    				final double [] zMapCorrWeights,
    				final int debugRow,
    				final int debugColumn,
    				final int threadsMax,
    				final boolean showProgress,
    				final int debugLevel){
        		final int debugThreshold=2;
        		final Rectangle zMapWOI=this.zMapWOI;
        		final ZTile [][] thisZMap=this.zMap[nImg];
        		final int tilesX=this.zMapWOI.width;
        		final int tiles=this.zMapWOI.width*this.zMapWOI.height;
        		if (debugLevel>2) System.out.println("setupZMap() zMapWOI.x="+zMapWOI.x+" zMapWOI.y="+zMapWOI.y+
        				" zMapWOI.width="+zMapWOI.width+" zMapWOI.height="+zMapWOI.height);
        		final Thread[] threads = newThreadArray(threadsMax);
		   		final AtomicInteger tileIndexAtomic     = new AtomicInteger(0);
		   		if (showProgress) IJ.showProgress(0.0);
		   		if (showProgress) IJ.showStatus("Filtering zMap foreground for image "+nImg+" ("+(nImg+1)+" of "+this.disparityScales.length+") ...");
				final int debugTile=debugRow*tilesX+debugColumn;
				//TODO: define here
				final double [] window=new double [4*this.overlapStep*this.overlapStep];
    			window[0]=Double.NaN;
				final double [] refineWindow=new double [refineTilePeriod*refineTilePeriod*4];
				refineWindow[0]=Double.NaN;
		   		for (int ithread = 0; ithread < threads.length; ithread++) {
		   			threads[ithread] = new Thread() {
		   				@Override
						public void run() {
		    				DoubleFHT doubleFHT = new DoubleFHT();
		    				DoubleFHT refineFHT= new DoubleFHT();
		   					for (int tile=tileIndexAtomic.getAndIncrement(); tile<tiles;tile=tileIndexAtomic.getAndIncrement()){
		   						boolean debugThis= (debugLevel>=debugThreshold) && (tile==debugTile);
		   						int thisDebugLevel=debugLevel+(debugThis?2:0);
//		   						while (!
		   								foregroundByOcclusionTile (
		   								zMapWOI.x+tile%tilesX, //int tileX,
		   								zMapWOI.y+tile/tilesX, //int tileY,
		   								thisZMap,
		   			    				nImg,
		   			    				sImg, // list of second images
		   			    				bgFraction, // process maximal correlation or the first farther plane with correlation intensity >= this part of maximal
		   			    				imageData, // [img]{Alpha, Y,Cb,Cr, Ext}. Alpha may be null - will handle later?
		   			    				imageFullWidth,
		   			    				blurVarianceSigma,
		   			    				refineTilePeriod,
		   			    				refinePhaseCoeff,
		   			    				refineHighPassSigma,
		   			    				refineLowPassSigma,
		   			    				refineCorrMaxDistance,
		   			    				refineCorrThreshold,
		   			    				refineSubPixel,
		   			    				zMapMinForeground,
		   			    				zMapVarMask,
		   			    				zMapVarThresholds,
		   			    				zMapVarWeights,
		   			    				normalizeCorrelation,
		   			    				zMapCorrMask,
		   			    				zMapCorrThresholds,
		   			    				zMapCorrWeights,
		   			    				window,
		   			    				doubleFHT,
		   			    				refineWindow,
		   			    				refineFHT,
		   			    				threadsMax,
		   			    				showProgress,
		   			    				thisDebugLevel); //);
		   						if (showProgress){
		   							final int finalTile=tile;
		   							SwingUtilities.invokeLater(new Runnable() {
		   								@Override
										public void run() {
		   									IJ.showProgress(finalTile,tiles);
		   								}
		   							});
		   						}
		   					}
		   				}
		   			};
		   		}
		   		startAndJoin(threads);
		   		IJ.showProgress(1.0);
        	}

    		private void foregroundByOcclusionTile ( // false if nothing left in the current foreground, may repeat
    				int tileX,
    				int tileY,
    				ZTile [][] thisZMap,
    				int nImg,
    				int [] sImgSet, // list of second images
    				double bgFraction, // process maximal correlation or the first farther plane with correlation intensity >= this part of maximal
    				double [][][] imageData, // [img]{Alpha, Y,Cb,Cr, Ext}. Alpha may be null - will handle later?
    				int imageFullWidth,
    				double blurVarianceSigma,
    				int refineTilePeriod,
    				double refinePhaseCoeff,
    				double refineHighPassSigma,
    				double refineLowPassSigma,
    				double refineCorrMaxDistance,
    				double refineCorrThreshold,
    				int refineSubPixel,
    				double zMapMinForeground,
    				int zMapVarMask,
    				double [] zMapVarThresholds,
    				double [] zMapVarWeights,
    				int     normalizeCorrelation,
    				int zMapCorrMask,
    				double [] zMapCorrThresholds,
    				double [] zMapCorrWeights,
    				double [] window,
    				DoubleFHT doubleFHT,
    				double [] refineWindow,
    				DoubleFHT refineFHT,
    				int threadsMax,
    				boolean showProgress,
    				int debugLevel){
    			double [] normVarWeights= normalizeWeights(zMapVarMask,zMapVarWeights);
    			double [] normCorrWeights=normalizeWeights(zMapCorrMask,zMapCorrWeights);
///    			final double [][] disparityScales=this.disparityScales;
// this filter will iterate through second images and AND filter for foreground pixels only
    			ZTile zTile=thisZMap[tileY][tileX];
    			if (zTile==null) return; // true; // tile is not defined, nothing left to do
    			// find the background plane to process
    			int bgPlane=0;
    			if (zTile.maximums.length>0){
    				for (int i=0;i<zTile.maximums.length;i++) if (zTile.maximums[i][1]>zTile.maximums[bgPlane][1]) bgPlane=i;
    				if (bgFraction>0){
    					int bgPlane1=bgPlane+1;
    					for (int i=bgPlane+1;i<zTile.maximums.length;i++) if (zTile.maximums[i][1]>zTile.maximums[bgPlane1][1]) bgPlane1=i;
    					if ((bgPlane1<zTile.maximums.length) && (zTile.maximums[bgPlane1][1]>=bgFraction*zTile.maximums[bgPlane][1])) bgPlane=bgPlane1;
    				} else if (bgFraction<0){ // compare with foreground
    					boolean hasFG=false;
    					for (int i=0;i<bgPlane;i++) if (zTile.maximums[i][1]>=(-bgFraction)*zTile.maximums[bgPlane][1]){
    						hasFG=true;
    						break;
    					}
    					if (!hasFG){
    						bgPlane++;
    						if (bgPlane<zTile.maximums.length){
    							for (int i=bgPlane+1;i<zTile.maximums.length;i++) if (zTile.maximums[i][1]>zTile.maximums[bgPlane][1]) bgPlane=i;
    						}
    					}
    				} else { // use gap-fileld FG
    					//zTile.foregroundIndex
    					if (zTile.foregroundIndex>=bgPlane){
    						bgPlane++;
    						if (bgPlane<zTile.maximums.length){
    							for (int i=bgPlane+1;i<zTile.maximums.length;i++) if (zTile.maximums[i][1]>zTile.maximums[bgPlane][1]) bgPlane=i;
    						}

    					}

    				}
    			}
    			// now bgPlane - background plane (if >=maximums.length - use infinity)
//				double bgFraction, // process maximal correlation or the first farther plane with correlation intensity >= this part of maximal

    			// is there any forground?
/*
    			if (!zTile.advanceForeground()) return true; // got to the end
    			if (zTile.enabledPixels[zTile.foregroundIndex]==null) {
    				zTile.initEnabledForegroundPixelsPlane();
    			}
*/
    			int tileOverlap=zTile.getOverlap();
    			int paddedSize=zTile.getPaddedSize();
    			int paddedLength=zTile.getPaddedLength();
//    			double disparity=zTile.maximums[zTile.foregroundIndex][0];
    			double disparity=(bgPlane>=zTile.maximums.length)?0.0:zTile.maximums[bgPlane][0];
    			int size=this.overlapStep*2;
    			int margin=size/4-tileOverlap;
    			int [] dirs1={1,size+1,size,size-1,-1,-size-1,-size,-size+1,0};
    			double [][][][] slices=new double [sImgSet.length][][][];
    			double [][][][] variance=new double [sImgSet.length][][][]; // [sIndex][pair 0/1][chn][pixel
    			int length=0; // initialize it to a correct value right here?
 //   			int numOther=sImgSet.length;
    			int [][] otherPairs=new int [(sImgSet.length*(sImgSet.length-1))/2][2];
    			// with 3 images there is only one other pair, but try to maki it work later with more images
    			{
    				int index=0;

    				for (int i=0;i<sImgSet.length-1;i++) for (int j=i+1;j<sImgSet.length;j++){
    					otherPairs[index][0]=i;
    					otherPairs[index++][1]=j;
    				}
    			}
    			int [] borderMask=zTile.getBorderMask();
    			if (debugLevel>4){
				(new ShowDoubleFloatArrays()).showArrays(
						borderMask,
						size,
						size,
						"borderMask");
    			}
//    			zTile.initAux(2);
//    			zTile.initAux(7);
    			zTile.initAux(11);
				if (debugLevel>3){ // +2 for selected tile
					System.out.println("otherPairs.length="+otherPairs.length+ "");
					for (int i=0;i<otherPairs.length;i++){
						System.out.println("--- nImg="+nImg+" otherPairs["+i+"][0]="+otherPairs[i][0]+" otherPairs["+i+"][1]="+otherPairs[i][1]);
					}
//					System.out.println("sImgSet.length="+sImgSet.length+ " numOther="+numOther);
					for (int i=0;i<sImgSet.length;i++){
						System.out.println("  nImg="+nImg+" sImgSet["+i+"]="+sImgSet[i]);
					}
				}
				float [] aux2=new float[paddedLength];
				for (int i=0;i<paddedLength;i++){ //tileOverlap
//					int iPix= (margin+(i/paddedSize))*size +(margin+ (i%paddedSize));
					aux2[i]=(float) disparity;
				}
				// debug
				zTile.setAux(2,aux2);

				for (int sIndex=0;sIndex<sImgSet.length;sIndex++){
					int sImg=sImgSet[sIndex];
					double dX=disparity*(this.disparityScales[sImg][0]-this.disparityScales[nImg][0]);
					double dY=disparity*(this.disparityScales[sImg][1]-this.disparityScales[nImg][1]);
//					double dXYlen=Math.sqrt(dX*dX+dY*dY);
//					double [] udXY={dX/dXYlen,dY/dXYlen};
					if (debugLevel>3){ // +2 for selected tile
						System.out.println("Debugging nImg="+nImg+" sImg="+sImg+" sIndex="+sIndex+" tileX="+tileX+", tileY="+tileY+" disparity="+disparity+" dX="+ dX+" dY="+dY);
					}
					slices[sIndex]=getShiftedSlices(
							tileX*this.overlapStep+this.overlapStep/2,
							tileY*this.overlapStep+this.overlapStep/2,
							this.overlapStep,
							nImg, //int first,
							sImg, //int second,
							dX, //double dxA,
							dY, //double dyA,
							zMapVarMask | zMapCorrMask, // used at least used in one filter,
							imageData,
							imageFullWidth,
							true, //doubleSizeOutput,
							window, // double [] window, // should be 4*size*size long;
							doubleFHT, //DoubleFHT doubleFHT, // to reuse tables
							debugLevel);

					if (debugLevel>3){ // +2 for selected tile
						System.out.println("Debugging tileX="+tileX+", tileY="+tileY+" disparity="+disparity+" dX="+ dX+" dY="+dY+
								" sIndex="+sIndex+" nImg="+nImg+" sImg="+sImg);
						int numActive=0;
						for (int i=0;i<slices[sIndex].length;i++) if (slices[sIndex][i]!=null) numActive++;
						String [] channelNames={"alpha","Y","Cb","Cr","Aux"};
						double [][] debugData=new double [2*numActive][];
						String [] debugTitles=new String [2*numActive];
						int index=0;
						for (int i=0;i<slices[sIndex].length;i++) if (slices[sIndex][i]!=null) {
							debugTitles[index]=channelNames[i/2+1]+"-"+nImg;
							debugData[index++]=slices[sIndex][i][0];
							debugTitles[index]=channelNames[i/2+1]+"-"+sImg;
							debugData[index++]=slices[sIndex][i][1];
						}
						(new ShowDoubleFloatArrays()).showArrays(
								debugData,
								size,
								size,
								true,
								"P"+nImg+"-"+sImg+"_sIndex-"+sIndex+"_tile"+tileX+"-"+tileY,
								debugTitles);
					}
				}
    			if (zMapVarMask!=0){
    				for (int sIndex=0;sIndex<sImgSet.length;sIndex++){
    					int sImg=sImgSet[sIndex];
    					variance [sIndex]= new double[2][slices[sIndex].length-1][];
    					for (int chn=0;chn<(slices[sIndex].length-1);chn++) if ((zMapVarMask& (1<<chn))!=0){
    						double [][] pair=slices[sIndex][chn+1]; // skip alpha
    						//        					int length=pair[0].length;
    						length=pair[0].length;
    						for (int n=0;n<2;n++) {
    							variance [sIndex][n][chn]= new double[length];
    							for (int i=0;i<length;i++) {
    								double  S1=0.0;
    								double  S2=0.0;
    								double  v2=0;
    								int numPix=0;
    								for (int d=0;d<dirs1.length;d++){
    									int i1=(i+dirs1[d]+length)%length;
    									if ((borderMask[i] & (borderMask[i1]<<1))==0) {
    										S1+=pair[n][i1];
    										S2+=(pair[n][i1]*pair[n][i1]);
    										numPix++;
    									}
    								}
    								v2=(S2-(S1*S1)/numPix)/numPix; //dirs1.length;
    								if ((debugLevel>4)&& (n==0) && (chn==0) && (numPix!=9) ){
    									System.out.println (" x="+(i%size)+" y="+(i/size)+" i="+i+" numPix="+numPix+" ="+borderMask[i]+" v2="+v2);
    								}
//    								v2+=(S2-(S1*S1)/numPix)/numPix; //dirs1.length;
    								variance[sIndex][n][chn][i]=Math.sqrt(v2); //java.lang.NullPointerException
    							}
    							if (blurVarianceSigma>0.0){
    								if (debugLevel>3) // +2 for selected tile
    		    						System.out.println("Bluring varaince["+sIndex+"]["+n+"]["+chn+"] with sigma="+blurVarianceSigma);
    								(new DoubleGaussianBlur()).blurDouble(
    										variance[sIndex][n][chn],
    										size,
    										size,
    										blurVarianceSigma,
    										blurVarianceSigma,
    										0.01);
    							}
    						}
    					}
    					if (debugLevel>3){ // +2 for selected tile
    						int numActive=0;
    						for (int i=0;i<slices[sIndex].length;i++) if (slices[sIndex][i]!=null) numActive++;
    						String [] channelNames={"alpha","Y","Cb","Cr","Aux"};
    						double [][] debugData=new double [2*numActive][];
    						String [] debugTitles=new String [2*numActive];
    						int index=0;
    						for (int i=0;i<slices[sIndex].length;i++) if (slices[sIndex][i]!=null) {
    							debugTitles[2*index]=channelNames[i/2+1]+"-"+nImg;
    							debugData[2*index]=variance[sIndex][0][index];
    							debugTitles[2*index+1]=channelNames[i/2+1]+"-"+sImg;
    							debugData[2*index+1]=variance[sIndex][1][index];
//    							if (variance[sIndex][0][index]==null) System.out.println("variance["+sIndex+"][0]["+index+"]==null");
//    							if (variance[sIndex][1][index]==null) System.out.println("variance["+sIndex+"][1]["+index+"]==null");
    							index++;
    						}
    						(new ShowDoubleFloatArrays()).showArrays(
    								debugData,
    								size,
    								size,
    								true,
    								"V"+nImg+"-"+sImg+"_sIndex-"+sIndex+"_tile"+tileX+"-"+tileY,
    								debugTitles);
    					}

    				}
    				double k=1.0;
    				//    					for (int iPix=0;iPix<length;iPix++){ // i - pixel index
    				float [] aux0=new float[paddedLength];
// debug layers
//    				float [] aux2=new float[paddedLength];
    				float [] aux3=new float[paddedLength];
    				float [] aux4=new float[paddedLength];
    				float [] aux5=new float[paddedLength];
    				float [] aux6=new float[paddedLength];
    				for (int i=0;i<paddedLength;i++){ //tileOverlap
    					int iPix= (margin+(i/paddedSize))*size +(margin+ (i%paddedSize));
    					// find pair with the smallest diff/diff_to_this ratio\
//    					int bestPair=0;
    					double occlVarVal=0.0;

						double pairVarDbgBest=0.0;
						double pairDiffDbgBest=0.0;
						double thisDiffDbgBest=0.0;
						double thisDiffVarBest=0.0;

    					for (int nPair=0;nPair<otherPairs.length;nPair++){
    						int sIndex0=otherPairs[nPair][0];
    						int sIndex1=otherPairs[nPair][1];
    						double occlVal=0.0;
// debug
    						double pairVarDbg=0.0;
    						double pairDiffDbg=0.0;
    						double thisDiffDbg=0.0;
    						double thisDiffVar=0.0;
    						for (int chn=0;chn<(slices[sIndex0].length-1);chn++) if ((zMapVarMask& (1<<chn))!=0){
 /*
  *    							if ((variance==null) ||
    									(variance[sIndex0]==null) ||
    									(variance[sIndex0][0]==null) ||
    									(variance[sIndex0][0][chn]==null) ||
    									(variance[sIndex0][1]==null) ||
    									(variance[sIndex0][1][chn]==null) ||
    									(variance[sIndex1]==null) ||
    									(variance[sIndex1][0]==null) ||
    									(variance[sIndex1][0][chn]==null) ||
    									(variance[sIndex1][1]==null) ||
    									(variance[sIndex1][1][chn]==null)){
    								System.out.println("foregroundByOcclusionTile("+tileX+","+tileY+",...): sIndex0="+sIndex0+"sIndex1="+sIndex1+" chn="+chn);
    								if (variance[sIndex0]==null) System.out.println("variance["+sIndex0+"]==null");
    								else {
    									if (variance[sIndex0][0]==null) System.out.println("variance["+sIndex0+"][0]==null");
    									else if (variance[sIndex0][0][chn]==null) System.out.println("variance["+sIndex0+"][0]["+chn+"]==null");
    									if (variance[sIndex0][1]==null) System.out.println("variance["+sIndex0+"][1]==null");
    									else if (variance[sIndex0][1][chn]==null) System.out.println("variance["+sIndex0+"][1]["+chn+"]==null");
    								}
    								if (variance[sIndex1]==null) System.out.println("variance["+sIndex1+"]==null");
    								else {
    									if (variance[sIndex1][0]==null) System.out.println("variance["+sIndex1+"][0]==null");
    									else if (variance[sIndex1][0][chn]==null) System.out.println("variance["+sIndex1+"][0]["+chn+"]==null");
    									if (variance[sIndex1][1]==null) System.out.println("variance["+sIndex1+"][1]==null");
    									else if (variance[sIndex1][1][chn]==null) System.out.println("variance["+sIndex1+"][1]["+chn+"]==null");
    								}
    							}
*/
    							double pairVar=Math.sqrt(0.25*(
    									variance[sIndex0][0][chn][iPix]*variance[sIndex0][0][chn][iPix]+
    									variance[sIndex1][0][chn][iPix]*variance[sIndex1][0][chn][iPix]+
    									variance[sIndex0][1][chn][iPix]*variance[sIndex0][1][chn][iPix]+
    									variance[sIndex1][1][chn][iPix]*variance[sIndex1][1][chn][iPix]));
    							if (debugLevel>4) {
    								System.out.print("chn="+chn+" nPair="+nPair+" i="+i+" iPix="+iPix+ " ix="+(i%paddedSize)+" iy="+(i/paddedSize)+
    										" ipx="+(iPix%size)+" ipy="+(iPix/size)+
    										" v["+sIndex0+"][0]["+chn+"]["+iPix+"]="+IJ.d2s(variance[sIndex0][0][chn][iPix],5)+
    										" v["+sIndex1+"][0]["+chn+"]["+iPix+"]="+IJ.d2s(variance[sIndex1][0][chn][iPix],5)+
    										" v["+sIndex0+"][1]["+chn+"]["+iPix+"]="+IJ.d2s(variance[sIndex0][1][chn][iPix],5)+
    										" v["+sIndex1+"][1]["+chn+"]["+iPix+"]="+IJ.d2s(variance[sIndex1][1][chn][iPix],5)+
    										" pairVar="+IJ.d2s(pairVar,5));
    							}
    							double pairDiff=Math.abs(slices[sIndex0][chn+1][1][iPix]-slices[sIndex1][chn+1][1][iPix])/pairVar/zMapVarThresholds[chn];
    							double thisDiff=Math.abs(0.5*(
    									slices[sIndex0][chn+1][0][iPix]+
    									slices[sIndex1][chn+1][0][iPix]-
    									slices[sIndex0][chn+1][1][iPix]-
    									slices[sIndex1][chn+1][1][iPix]));
    							occlVal+=normVarWeights[chn]*(thisDiff/pairVar/zMapVarThresholds[chn])/(pairDiff+k); // consider both - with /pairVar and without

        						pairVarDbg+=normVarWeights[chn]*pairVar;
        						pairDiffDbg+=normVarWeights[chn]*pairDiff;
        						thisDiffDbg+=normVarWeights[chn]*thisDiff;
        						double dd=slices[sIndex0][chn+1][1][iPix]-slices[sIndex1][chn+1][1][iPix];
        						thisDiffVar+=normVarWeights[chn]*Math.sqrt(0.5*(dd*dd+pairVar*pairVar));

    							// other variant - if pairDiff>1.0 - =>    occlVal=0.0;

    						}
    						if ((nPair==0) || (occlVarVal<occlVal)){
    							occlVarVal=occlVal;
//    							bestPair=nPair;

    							//debug
    							pairVarDbgBest=pairVarDbg;
    							pairDiffDbgBest=pairDiffDbg;
    							thisDiffDbgBest=thisDiffDbg;
    							thisDiffVarBest=thisDiffVar;

    						}
							if (debugLevel>4) {
								System.out.println(	" pairVarDbgBest="+IJ.d2s(pairVarDbgBest,5));
							}

    					}
    					aux0[i]=(float) occlVarVal;
// debug
//    					aux2[i]=(float) disparity;
    					aux3[i]=(float) pairVarDbgBest;
    					aux4[i]=(float) pairDiffDbgBest;
    					aux5[i]=(float) thisDiffDbgBest;
    					aux6[i]=(float) thisDiffVarBest;


    				}
    				zTile.setAux(0,aux0);
// debug
//    				zTile.setAux(2,aux2);
    				zTile.setAux(3,aux3);
    				zTile.setAux(4,aux4);
    				zTile.setAux(5,aux5);
    				zTile.setAux(6,aux6);

    				if (debugLevel>3){
						float [][] debugData={aux2,aux3,aux4,aux5,aux6,aux0};
						String [] debugTitles={"var","pair","this","occl","dvar","disp"};
						(new ShowDoubleFloatArrays()).showArrays(
								debugData,
								paddedSize,
								paddedSize,
								true,
								"R"+nImg+"_tile"+tileX+"-"+tileY,
								debugTitles);
    				}
    			}
    			if (zMapCorrMask!=0) { // 8x8 correlation to find occlusion
    				double [][][]refineSmoothThis= new double [sImgSet.length][][];
    				double [][][]refineSmoothOther=new double [otherPairs.length][][];
    				for (int sIndex=0;sIndex<sImgSet.length;sIndex++){
    					int sImg=sImgSet[sIndex];
    					double dX=disparity*(this.disparityScales[sImg][0]-this.disparityScales[nImg][0]);
    					double dY=disparity*(this.disparityScales[sImg][1]-this.disparityScales[nImg][1]);
    					double dXYlen=Math.sqrt(dX*dX+dY*dY);
    					double [] udXY={dX/dXYlen,dY/dXYlen};
    					refineSmoothThis[sIndex]=new double [slices[sIndex].length-1][];
    					for (int chn=0;chn<(slices.length-1);chn++) if ((zMapCorrMask& (1<<chn))!=0){
    						double [][] pair=slices[sIndex][chn+1];
    						refineSmoothThis[sIndex][chn]=localCorrelation(
    								refineFHT, // DoubleFHT doubleFHT, // will generate if null, can be used to reuse initialization
    								pair, //double [][] data,   // difines size
    								refineWindow, //double [] window, // defines tile size
    								udXY, // double [] dXY,    // unity vector defines disparity direction
    								refinePhaseCoeff, //phaseCoeff, // phase correlation fraction (0 - normal, 1.0 - pure phase)
    								refineHighPassSigma,
    								refineLowPassSigma,
    			    				normalizeCorrelation,
    								debugLevel,
    								(debugLevel>3)?("RST"+nImg+"-"+sImg+"_"+chn):null);
    					}
    				}
    				for (int nPair=0;nPair<otherPairs.length;nPair++){
    					int sIndex0=otherPairs[nPair][0];
    					int sIndex1=otherPairs[nPair][1];
    					int sImg0=sImgSet[sIndex0];
    					int sImg1=sImgSet[sIndex1];
    					double dX=disparity*(this.disparityScales[sImg1][0]-this.disparityScales[sImg0][0]);
    					double dY=disparity*(this.disparityScales[sImg1][1]-this.disparityScales[sImg0][1]);
    					double dXYlen=Math.sqrt(dX*dX+dY*dY);
    					double [] udXY={dX/dXYlen,dY/dXYlen};
    					refineSmoothOther[nPair]=new double [slices[sIndex0].length-1][];
    					for (int chn=0;chn<(slices.length-1);chn++) if ((zMapCorrMask& (1<<chn))!=0){
    						double [][] pair={slices[sIndex0][1][chn+1],slices[sIndex1][1][chn+1]};
    						refineSmoothOther[nPair][chn]=localCorrelation(
    								refineFHT, // DoubleFHT doubleFHT, // will generate if null, can be used to reuse initialization
    								pair, //double [][] data,   // difines size
    								refineWindow, //double [] window, // defines tile size
    								udXY, // double [] dXY,    // unity vector defines disparity direction
    								refinePhaseCoeff, //phaseCoeff, // phase correlation fraction (0 - normal, 1.0 - pure phase)
    								refineHighPassSigma,
    								refineLowPassSigma,
    			    				normalizeCorrelation,
    								debugLevel,
    								(debugLevel>3)?("RS-"+nImg+"-"+sImg0+"-"+sImg1+"_"+chn):null);
    					}
    				}
    				float [] aux1=new float[paddedLength];
//
    				float [] aux7=new float[paddedLength];
    				float [] aux8=new float[paddedLength];
    				float [] aux9=new float[paddedLength];
    				float [] aux10=new float[paddedLength];

    				for (int i=0;i<paddedLength;i++){ //tileOverlap
    					int iPix= (margin+(i/paddedSize))*size +(margin+ (i%paddedSize));
    					// zMapCorrThresholds - here weight? Maybe - add separate weights and use thersholds for thersholds?
    					double occlCorrVal=0.0;
    					double occlCorrOther=0.0;
    					double occlCorrFirst=0.0;
    					double occlCorrSecond=0.0;
    					for (int chn=0;chn<(slices.length-1);chn++) if ((zMapCorrMask& (1<<chn))!=0){
    						//								int bestPair=0;
    						double bestOcclCorrVal=0.0;
        					double bestOcclCorrOther=0.0;
        					double bestOcclCorrFirst=0.0;
        					double bestOcclCorrSecond=0.0;
    						for (int nPair=0;nPair<otherPairs.length;nPair++){
    							int sIndex0=otherPairs[nPair][0];
    							int sIndex1=otherPairs[nPair][1];
    							double thisOcclCorrVal=2*refineSmoothOther[nPair][chn][iPix]-refineSmoothThis[sIndex0][chn][iPix]-refineSmoothThis[sIndex1][chn][iPix];
    							if ((nPair==0) || (thisOcclCorrVal>bestOcclCorrVal)) {
    								bestOcclCorrVal=thisOcclCorrVal;
    								bestOcclCorrOther=refineSmoothOther[nPair][chn][iPix];
    	        					bestOcclCorrFirst=refineSmoothThis[sIndex0][chn][iPix];
    	        					bestOcclCorrSecond=refineSmoothThis[sIndex1][chn][iPix];
    								//										bestPair=nPair;
    							}
    						}
    						occlCorrVal+=normCorrWeights[chn]*bestOcclCorrVal;
        					occlCorrOther+=normCorrWeights[chn]*bestOcclCorrOther;
        					occlCorrFirst+=normCorrWeights[chn]*bestOcclCorrFirst;
        					occlCorrSecond+=normCorrWeights[chn]*bestOcclCorrSecond;
    					}
    					// Save result to aux[1]
    					aux1[i]=(float) occlCorrVal;
    					aux7[i]=(float) occlCorrOther;
    					aux8[i]=(float) (occlCorrFirst+occlCorrSecond)/2;
    					aux9[i]=(float) occlCorrFirst;
    					aux10[i]=(float) occlCorrSecond;
    				}
    				zTile.setAux(1,aux1);
    				zTile.setAux(7,aux7);
    				zTile.setAux(8,aux8);
    				zTile.setAux(9,aux9);
    				zTile.setAux(10,aux10);
    				if (debugLevel>3){
						float [][] debugData={aux1,aux7,aux8,aux9,aux10};
						String [] debugTitles={"rslt","other","average","first","second"};
						(new ShowDoubleFloatArrays()).showArrays(
								debugData,
								paddedSize,
								paddedSize,
								true,
								"C"+nImg+"_tile"+tileX+"-"+tileY,
								debugTitles);
    				}

    			}
    			return; // zTile.isForegroundValid(); // this current plane is not empty
    		}
    		public double [] normalizeWeights(int mask,double [] weights){
    			double [] maskedWeights=weights.clone();
    			for (int i=0;i<maskedWeights.length;i++) if ((mask & (1<<i))==0) maskedWeights[i]=0.0;
    			return normalizeWeights(maskedWeights);
    		}
    		public double [] normalizeWeights(double [] weights){
    			double [] normalizedWeights=new double [weights.length];
    			double sum=0.0;
    			for (int i=0;i<weights.length;i++) sum+=weights[i];
    			for (int i=0;i<weights.length;i++) normalizedWeights[i]=(sum==0.0)?0.0:(weights[i]/sum);
    			return normalizedWeights;
    		}


    		public void filterForegroundZMap(
//    				final int nImg,
//    				final int [] sImg, // list of second images
    				double [][][] imageData, // [img]{Alpha, Y,Cb,Cr, Ext}. Alpha may be null - will handle later?
    				int       imageFullWidth,
    				double    blurVarianceSigma,
    				int       refineTilePeriod,
    				double    refinePhaseCoeff,
    				double    refineHighPassSigma,
    				double    refineLowPassSigma,
    				double    refineCorrMaxDistance,
    				double    refineCorrThreshold,
    				int       refineSubPixel,
    				double    zMapMinForeground,
    				int       zMapVarMask,
    				double [] zMapVarThresholds,
    				double [] zMapVarWeights,
    				int       auxVarMode,
    				int       normalizeCorrelation,
    				int       zMapCorrMask,
    				double [] zMapCorrThresholds,
    				double [] zMapCorrWeights,
    				int       debugRow,
    				int       debugColumn,
    				int       threadsMax,
    				boolean   showProgress,
    				int       debugLevel){
    			// process all image pairs
    			for (int nImg=0;nImg<this.disparityScales.length;nImg++){
    				int [] sImg=new int [this.disparityScales.length-1];
    				for (int j=0;j<sImg.length;j++)sImg[j]=(j<nImg)?j:(j+1);
    				filterForegroundZMap(
    	    				nImg,
    	    				sImg, // list of second images
    	    				imageData, // [img]{Alpha, Y,Cb,Cr, Ext}. Alpha may be null - will handle later?
    	    				imageFullWidth,
    	    				blurVarianceSigma,
    	    				refineTilePeriod,
    	    				refinePhaseCoeff,
    	    				refineHighPassSigma,
    	    				refineLowPassSigma,
    	    				refineCorrMaxDistance,
    	    				refineCorrThreshold,
    	    				refineSubPixel,
    	    				zMapMinForeground,
    	    				zMapVarMask,
    	    				zMapVarThresholds,
    	    				zMapVarWeights,
    	    				auxVarMode,
    	    				normalizeCorrelation,
    	    				zMapCorrMask,
    	    				zMapCorrThresholds,
    	    				zMapCorrWeights,
    	    				debugRow,
    	    				debugColumn,
    	    				threadsMax,
    	    				showProgress,
    	    				debugLevel);
    			}
    		}

        	public void filterForegroundZMap(
    				final int nImg,
    				final int [] sImg, // list of second images
    				final double [][][] imageData, // [img]{Alpha, Y,Cb,Cr, Ext}. Alpha may be null - will handle later?
    				final int imageFullWidth,
    				final double blurVarianceSigma,
    				final int refineTilePeriod,
    				final double refinePhaseCoeff,
    				final double refineHighPassSigma,
    				final double refineLowPassSigma,
    				final double refineCorrMaxDistance,
    				final double refineCorrThreshold,
    				final int refineSubPixel,
    				final double zMapMinForeground,
    				final int zMapVarMask,
    				final double [] zMapVarThresholds,
    				final double [] zMapVarWeights,
    				final int       auxVarMode,
        			final int     normalizeCorrelation,
    				final int zMapCorrMask,
    				final double [] zMapCorrThresholds,
    				final double [] zMapCorrWeights,
    				final int debugRow,
    				final int debugColumn,
    				final int threadsMax,
    				final boolean showProgress,
    				final int debugLevel){
        		final int debugThreshold=2;
        		final Rectangle zMapWOI=this.zMapWOI;
        		final ZTile [][] thisZMap=this.zMap[nImg];
        		final int tilesX=this.zMapWOI.width;
        		final int tiles=this.zMapWOI.width*this.zMapWOI.height;
        		if (debugLevel>2) System.out.println("setupZMap() zMapWOI.x="+zMapWOI.x+" zMapWOI.y="+zMapWOI.y+
        				" zMapWOI.width="+zMapWOI.width+" zMapWOI.height="+zMapWOI.height);
        		final Thread[] threads = newThreadArray(threadsMax);
		   		final AtomicInteger tileIndexAtomic     = new AtomicInteger(0);
		   		if (showProgress) IJ.showProgress(0.0);
		   		if (showProgress) IJ.showStatus("Filtering zMap foreground for image "+nImg+" ("+(nImg+1)+" of "+this.disparityScales.length+") ...");
				final int debugTile=debugRow*tilesX+debugColumn;
				//TODO: define here
				final double [] window=new double [4*this.overlapStep*this.overlapStep];
    			window[0]=Double.NaN;
				final double [] refineWindow=new double [refineTilePeriod*refineTilePeriod*4];
				refineWindow[0]=Double.NaN;
		   		for (int ithread = 0; ithread < threads.length; ithread++) {
		   			threads[ithread] = new Thread() {
		   				@Override
						public void run() {
		    				DoubleFHT doubleFHT = new DoubleFHT();
		    				DoubleFHT refineFHT= new DoubleFHT();
		   					for (int tile=tileIndexAtomic.getAndIncrement(); tile<tiles;tile=tileIndexAtomic.getAndIncrement()){
		   						boolean debugThis= (debugLevel>=debugThreshold) && (tile==debugTile);
		   						int thisDebugLevel=debugLevel+(debugThis?2:0);
		   						while (!filterForegroundZMapTile (
		   								zMapWOI.x+tile%tilesX, //int tileX,
		   								zMapWOI.y+tile/tilesX, //int tileY,
		   								thisZMap,
		   			    				nImg,
		   			    				sImg, // list of second images
		   			    				imageData, // [img]{Alpha, Y,Cb,Cr, Ext}. Alpha may be null - will handle later?
		   			    				imageFullWidth,
		   			    				blurVarianceSigma,
		   			    				refineTilePeriod,
		   			    				refinePhaseCoeff,
		   			    				refineHighPassSigma,
		   			    				refineLowPassSigma,
		   			    				refineCorrMaxDistance,
		   			    				refineCorrThreshold,
		   			    				refineSubPixel,
		   			    				zMapMinForeground,
		   			    				zMapVarMask,
		   			    				zMapVarThresholds,
		   			    				zMapVarWeights,
		   			    				auxVarMode,
		   			    				normalizeCorrelation,
		   			    				zMapCorrMask,
		   			    				zMapCorrThresholds,
		   			    				zMapCorrWeights,
		   			    				window,
		   			    				doubleFHT,
		   			    				refineWindow,
		   			    				refineFHT,
		   			    				threadsMax,
		   			    				showProgress,
		   			    				thisDebugLevel));
		   						if (showProgress){
		   							final int finalTile=tile;
		   							SwingUtilities.invokeLater(new Runnable() {
		   								@Override
										public void run() {
		   									IJ.showProgress(finalTile,tiles);
		   								}
		   							});
		   						}
		   					}
		   				}
		   			};
		   		}
		   		startAndJoin(threads);
		   		IJ.showProgress(1.0);
        	}

    		private boolean filterForegroundZMapTile ( // false if nothing left in the current foreground, may repeat
    				int tileX,
    				int tileY,
    				ZTile [][] thisZMap,
    				int nImg,
    				int [] sImgSet, // list of second images
    				double [][][] imageData, // [img]{Alpha, Y,Cb,Cr, Ext}. Alpha may be null - will handle later?
    				int imageFullWidth,
    				double blurVarianceSigma,
    				int refineTilePeriod,
    				double refinePhaseCoeff,
    				double refineHighPassSigma,
    				double refineLowPassSigma,
    				double refineCorrMaxDistance,
    				double refineCorrThreshold,
    				int refineSubPixel,
    				double zMapMinForeground,
    				int zMapVarMask,
    				double [] zMapVarThresholds,
    				double [] zMapVarWeights,
    				int       auxVarMode,
    				int     normalize,
    				int zMapCorrMask,
    				double [] zMapCorrThresholds,
    				double [] zMapCorrWeights,
    				double [] window,
    				DoubleFHT doubleFHT,
    				double [] refineWindow,
    				DoubleFHT refineFHT,
    				int threadsMax,
    				boolean showProgress,
    				int debugLevel){
//    			double [] normVarWeights= normalizeWeights(zMapVarMask,zMapVarWeights);
//    			double [] normCorrWeights=normalizeWeights(zMapCorrMask,zMapCorrWeights);
                int auxChannelNumber=3;
///    			final double [][] disparityScales=this.disparityScales;
// this filter will iterate through second images and AND filter for foreground pixels only
    			ZTile zTile=thisZMap[tileY][tileX];
    			if (zTile==null) return true; // tile is not defined, nothing left to do
    			// is there any forground?
    			if (!zTile.advanceForeground()) return true; // got to the end
    			if (zTile.enabledPixels[zTile.foregroundIndex]==null) {
    				zTile.initEnabledForegroundPixelsPlane();
    			}
    			int tileOverlap=zTile.getOverlap();
    			int paddedSize=zTile.getPaddedSize();
    			int paddedLength=zTile.getPaddedLength();
    			double disparity=zTile.maximums[zTile.foregroundIndex][0];
    			int size=this.overlapStep*2;
    			int margin=size/4-tileOverlap;
    			int [] dirs1={1,size+1,size,size-1,-1,-size-1,-size,-size+1,0};
    			for (int sIndex=0;sIndex<sImgSet.length;sIndex++){
    				int sImg=sImgSet[sIndex];
    				double dX=disparity*(this.disparityScales[sImg][0]-this.disparityScales[nImg][0]);
    				double dY=disparity*(this.disparityScales[sImg][1]-this.disparityScales[nImg][1]);
    				double dXYlen=Math.sqrt(dX*dX+dY*dY);
    				double [] udXY={dX/dXYlen,dY/dXYlen};
    				double [][][] slices=getShiftedSlices(
    						tileX*this.overlapStep+this.overlapStep/2,
    						tileY*this.overlapStep+this.overlapStep/2,
    						this.overlapStep,
    						nImg, //int first,
    						sImg, //int second,
    						dX, //double dxA,
    						dY, //double dyA,
    						zMapVarMask | zMapCorrMask, // used at least used in one filter,
    						imageData,
    						imageFullWidth,
    						true, //doubleSizeOutput,
    						window, // double [] window, // should be 4*size*size long;
    						doubleFHT, //DoubleFHT doubleFHT, // to reuse tables
    						debugLevel);

    				if (debugLevel>3){ // +2 for selected tile
    					System.out.println("Debugging tileX="+tileX+", tileY="+tileY+" disparity="+disparity+" dX="+ dX+" dY="+dY);
    					int numActive=0;
    					for (int i=0;i<slices.length;i++) if (slices[i]!=null) numActive++;
    					String [] channelNames={"alpha","Y","Cb","Cr","Aux"};
    					double [][] debugData=new double [2*numActive][];
    					String [] debugTitles=new String [2*numActive];
    					int index=0;
    					for (int i=0;i<slices.length;i++) if (slices[i]!=null) {
    						debugTitles[index]=channelNames[i/2+1]+"-"+nImg;
    						debugData[index++]=slices[i][0];
    						debugTitles[index]=channelNames[i/2+1]+"-"+sImg;
    						debugData[index++]=slices[i][1];
    					}
    					String title="slice_"+tileX+"-"+tileY+"_"+nImg+"-"+sImg;
    					(new ShowDoubleFloatArrays()).showArrays(
    							debugData,
    							size,
    							size,
    							true,
    							title,
    							debugTitles);
    					System.out.println("\n\nProcessing slice "+title+": ");
    					System.out.print("Index\tX\tY");
    					for (int i=0;i<debugTitles.length;i++)System.out.print("\t"+debugTitles[i]);
    					System.out.println();
    					for (int i=0;i<(size*size);i++){
    						System.out.print(i+"\t"+(i%size)+"\t"+(i/size));
    						for (int j=0;j<debugData.length;j++){
    							System.out.print("\t"+debugData[j][i]);
    						}
    						System.out.println();
    					}
    				}
    				for (int chn=0;chn<(slices.length-1);chn++) if ((zMapVarMask& (1<<chn))!=0){
						double [][] pair=slices[chn+1]; // skip alpha
						int length=pair[0].length;
						double [] bDiff=new double[length];
    					if (chn!=auxChannelNumber){
    						// will be overkill here, later can be optimized
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
    							bDiff[i]=bdiff;
    						}
    						double [] variance=new double[length];
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
    							variance[i]=Math.sqrt(v2);
    						}
    						if (blurVarianceSigma>0.0){
    							(new DoubleGaussianBlur()) .blurDouble(
    									variance,
    									size,
    									size,
    									blurVarianceSigma,
    									blurVarianceSigma,
    									0.01);
    						}




    						for (int i=0;i<paddedLength;i++){ //tileOverlap
    							int i1= (margin+(i/paddedSize))*size +(margin+ (i%paddedSize));
    							if ((i1>=bDiff.length) || (i1>=variance.length)){
    								System.out.println("\nfilterForegroundZMapTile() i="+i+" margin="+margin+" paddedSize="+paddedSize+" bDiff.length="+bDiff.length+
    										" i1="+i1+" size="+size+" variance.length="+variance.length);
    							}
    							if ((bDiff[i1]/variance[i1])>zMapVarThresholds[chn]) { // 1024
    								zTile.disableForegroundPixel(i);
    							}
    						}
    						if (debugLevel>3){ // +2 for selected tile
    							String [] channelNames={"alpha","Y","Cb","Cr","Aux"};
    							String [] debugTitles={"bDiff","var","bdif/var","good","cumul"};
    							double [][] debugData=new double[5][length];
    							debugData[0]=bDiff;
    							debugData[1]=variance;
    							for (int i=0;i<length;i++){
    								debugData[2][i]=bDiff[i]/variance[i];
    								debugData[3][i]=0.0;
    								debugData[4][i]=0.0;
    							}

    							for (int i=0;i<paddedLength;i++){ //tileOverlap
    								int i1= (margin+(i/paddedSize))*size +(margin+ (i%paddedSize));
    								debugData[3][i]=((bDiff[i1]/variance[i1])>zMapVarThresholds[chn])?-1.0:1.0;
    								debugData[4][i]=zTile.isEnabledForegroundPixel(i)?-1.0:1.0;
    							}
    							(new ShowDoubleFloatArrays()).showArrays(
    									debugData,
    									size,
    									size,
    									true,
    									"diff_"+channelNames[chn+1]+"_"+tileX+"-"+tileY+"_"+nImg+"-"+sImg,
    									debugTitles);
    						}



    					} else { // aux channel
    						for (int i=0;i<length;i++) 	bDiff[i]=pair[0][i]*pair[1][i];
    						if (debugLevel>3){ // +2 for selected tile
    							(new ShowDoubleFloatArrays()).showArrays(
    									bDiff,
    									size,
    									size,
    									"AUX_"+nImg+"-"+sImg+"__"+tileX+"-"+tileY);

    						}
    						for (int i=0;i<paddedLength;i++){ //tileOverlap
    							int i1= (margin+(i/paddedSize))*size +(margin+ (i%paddedSize));
    							if (bDiff[i1]<zMapVarThresholds[chn]) { // 1024
    								zTile.disableForegroundPixel(i);
    							}
    						}
    					}
    				}
    				for (int chn=0;chn<(slices.length-1);chn++) if ((zMapCorrMask& (1<<chn))!=0){
    					double [][] pair=slices[chn+1];
//    	    			int length=pair[0].length;
    					double [] refineSmooth=localCorrelation(
    							refineFHT, // DoubleFHT doubleFHT, // will generate if null, can be used to reuse initialization
    							pair, //double [][] data,   // difines size
    							refineWindow, //double [] window, // defines tile size
    							udXY, // double [] dXY,    // unity vector defines disparity direction
    							refinePhaseCoeff, //phaseCoeff, // phase correlation fraction (0 - normal, 1.0 - pure phase)
    							refineHighPassSigma,
    							refineLowPassSigma,
    		        			normalize,
    			    			debugLevel,
    			    			(debugLevel>3)?("CORR"+nImg+"-"+sImg+"_"+chn):null);
    					for (int i=0;i<paddedLength;i++){ //tileOverlap
    						int i1= (margin+(i/paddedSize))*size +(margin+ (i%paddedSize));
    						if (refineSmooth[i1]<zMapCorrThresholds[chn]) {
    							zTile.disableForegroundPixel(i);
    						}
    					}
    				}

    			}
    			return zTile.isForegroundValid(); // this current plane is not empty
    		}

    		public void refinePlaneDisparity(
//    				final int nImg,
//    				final int [] sImg, // list of second images
    				double [][][] imageData, // [img]{Alpha, Y,Cb,Cr, Ext}. Alpha may be null - will handle later?
    				int       imageFullWidth,
//    				double    blurVarianceSigma,
//    				int       refineTilePeriod,
    				double    refinePhaseCoeff,
    				double    refineHighPassSigma,
    				double    refineLowPassSigma,
    				double    refineCorrMaxDistance,
    				double    refineCorrThreshold,
    				int       refineSubPixel,
//    				double    zMapMinForeground,
//    				int       zMapVarMask,
//    				double [] zMapVarThresholds,
//    				double [] zMapVarWeights,
//    				int       auxVarMode,
//    				int       normalizeCorrelation,
    				int       zMapCorrMask,
    				double [] zMapCorrThresholds,
    				double [] zMapCorrWeights,
    				int       debugRow,
    				int       debugColumn,

    				double    disparityMax,
    				double    disparityMin,
    				double    minAbsolute, // or NaN - will use enabled/disabled state of the tile
    				double    minRelative,
    				boolean   filterByForeground, // apply known certain masks
    				double    filterByForegroundMargin,
    				boolean   filterByDisabled,
    				double    disparityTolearnce,
    				double    maskBlurSigma,
    				double    corrHighPassSigma, // subtract blurred version to minimize correlation caused by masks

    				int       threadsMax,
    				boolean   showProgress,
    				int       debugLevel){
    			// process all image pairs
    			for (int nImg=0;nImg<this.disparityScales.length;nImg++){
    				int [] sImg=new int [this.disparityScales.length-1];
    				for (int j=0;j<sImg.length;j++)sImg[j]=(j<nImg)?j:(j+1);
    				refinePlaneDisparity(
    	    				nImg,
    	    				sImg, // list of second images
    	    				imageData, // [img]{Alpha, Y,Cb,Cr, Ext}. Alpha may be null - will handle later?
    	    				imageFullWidth,
//    	    				blurVarianceSigma,
//    	    				refineTilePeriod,
    	    				refinePhaseCoeff,
    	    				refineHighPassSigma,
    	    				refineLowPassSigma,
    	    				refineCorrMaxDistance,
    	    				refineCorrThreshold,
    	    				refineSubPixel,
//    	    				zMapMinForeground,
//    	    				zMapVarMask,
//    	    				zMapVarThresholds,
//    	    				zMapVarWeights,
//    	    				auxVarMode,
//    	    				normalizeCorrelation,
    	    				zMapCorrMask,
    	    				zMapCorrThresholds,
    	    				zMapCorrWeights,
    	    				debugRow,
    	    				debugColumn,
    	    				disparityMax,
    	    				disparityMin,
    	    				minAbsolute, // or NaN - will use enabled/disabled state of the tile
    	    				minRelative,
    	    				filterByForeground, // apply known certain masks
    	    				filterByForegroundMargin,
    	    				filterByDisabled,
    	    				disparityTolearnce,
    	    				maskBlurSigma,
    	    				corrHighPassSigma, // subtract blurred version to minimize correlation caused by masks

    	    				threadsMax,
    	    				showProgress,
    	    				debugLevel);
    			}
    		}

        	public void refinePlaneDisparity(
    				final int nImg,
    				final int [] sImg, // list of second images
    				final double [][][] imageData, // [img]{Alpha, Y,Cb,Cr, Ext}. Alpha may be null - will handle later?
    				final int imageFullWidth,
//    				final double blurVarianceSigma,
//    				final int refineTilePeriod,
    				final double refinePhaseCoeff,
    				final double refineHighPassSigma,
    				final double refineLowPassSigma,
    				final double refineCorrMaxDistance,
    				final double refineCorrThreshold,
    				final int refineSubPixel,
//    				final double zMapMinForeground,
//    				final int zMapVarMask,
//    				final double [] zMapVarThresholds,
//    				final double [] zMapVarWeights,
//    				final int       auxVarMode,
//        			final int     normalizeCorrelation,
    				final int zMapCorrMask,
    				final double [] zMapCorrThresholds,
    				final double [] zMapCorrWeights,
    				final int debugRow,
    				final int debugColumn,

    				final double disparityMax,
    				final double disparityMin,
    				final double minAbsolute, // or NaN - will use enabled/disabled state of the tile
    				final double minRelative,
    				final boolean filterByForeground, // apply known certain masks
    				final double filterByForegroundMargin,
    				final boolean filterByDisabled,
    				final double disparityTolearnce,
    				final double maskBlurSigma,
    				final double corrHighPassSigma, // subtract blurred version to minimize correlation caused by masks


    				final int threadsMax,
    				final boolean showProgress,
    				final int debugLevel){
        		final int debugThreshold=2;
        		final Rectangle zMapWOI=this.zMapWOI;
        		final ZTile [][][] zMapFinal=this.zMap;
        		final int tilesX=this.zMapWOI.width;
        		final int tiles=this.zMapWOI.width*this.zMapWOI.height;
        		if (debugLevel>2) System.out.println("setupZMap() zMapWOI.x="+zMapWOI.x+" zMapWOI.y="+zMapWOI.y+
        				" zMapWOI.width="+zMapWOI.width+" zMapWOI.height="+zMapWOI.height);
        		final Thread[] threads = newThreadArray(threadsMax);
		   		final AtomicInteger tileIndexAtomic     = new AtomicInteger(0);
		   		if (showProgress) IJ.showProgress(0.0);
		   		if (showProgress) IJ.showStatus("refinePlaneDisparity for image "+nImg+" ("+(nImg+1)+" of "+this.disparityScales.length+") ...");
				final int debugTile=debugRow*tilesX+debugColumn;
				//TODO: define here
//				final double [] window=new double [4*this.overlapStep*this.overlapStep];
//    			window[0]=Double.NaN;
//				final double [] refineWindow=new double [refineTilePeriod*refineTilePeriod*4];
//				refineWindow[0]=Double.NaN;
				final int size=2*this.overlapStep;
		   		for (int ithread = 0; ithread < threads.length; ithread++) {
		   			threads[ithread] = new Thread() {
		   				@Override
						public void run() {
		    				DoubleFHT doubleFHT = new DoubleFHT();
		    				final double [] window=new double [size*size];
		        			window[0]=Double.NaN;
//		    				DoubleFHT refineFHT= new DoubleFHT();
		   					for (int tile=tileIndexAtomic.getAndIncrement(); tile<tiles;tile=tileIndexAtomic.getAndIncrement()){
		   						boolean debugThis= (debugLevel>=debugThreshold) && (tile==debugTile);
		   						int thisDebugLevel=debugLevel+(debugThis?2:0);
		   						refinePlaneDisparityTile (
		   								zMapWOI.x+tile%tilesX, //int tileX,
		   								zMapWOI.y+tile/tilesX, //int tileY,
		   								zMapFinal,
		   			    				nImg,
		   			    				sImg, // list of second images
		   			    				imageData, // [img]{Alpha, Y,Cb,Cr, Ext}. Alpha may be null - will handle later?
		   			    				imageFullWidth,
//		   			    				blurVarianceSigma,
//		   			    				refineTilePeriod,
		   			    				refinePhaseCoeff,
		   			    				refineHighPassSigma,
		   			    				refineLowPassSigma,
		   			    				refineCorrMaxDistance,
		   			    				refineCorrThreshold,
		   			    				refineSubPixel,
//		   			    				zMapMinForeground,
//		   			    				zMapVarMask,
//		   			    				zMapVarThresholds,
//		   			    				zMapVarWeights,
//		   			    				auxVarMode,
//		   			    				normalizeCorrelation,
		   			    				zMapCorrMask,
		   			    				zMapCorrThresholds,
		   			    				zMapCorrWeights,
		   			    				window,
		   			    				doubleFHT,
//		   			    				refineWindow,
//		   			    				refineFHT,
		   			    				disparityMax,
		   			    				disparityMin,
		   			    				minAbsolute, // or NaN - will use enabled/disabled state of the tile
		   			    				minRelative,
		   			    				filterByForeground, // apply known certain masks
		   			    				filterByForegroundMargin,
		   			    				filterByDisabled,
		   			    				disparityTolearnce,
		   			    				maskBlurSigma,
		   			    				corrHighPassSigma, // subtract blurred version to minimize correlation caused by masks

		   			    				threadsMax,
		   			    				showProgress,
		   			    				thisDebugLevel);
		   						if (showProgress){
		   							final int finalTile=tile;
		   							SwingUtilities.invokeLater(new Runnable() {
		   								@Override
										public void run() {
		   									IJ.showProgress(finalTile,tiles);
		   								}
		   							});
		   						}
		   					}
		   				}
		   			};
		   		}
		   		startAndJoin(threads);
		   		IJ.showProgress(1.0);
        	}





       		private void refinePlaneDisparityTile ( // false if nothing left in the current foreground, may repeat
    				int tileX,
    				int tileY,
    				ZTile [][][] allZMap,
    				int nImg,
    				int [] sImgSet, // list of second images
    				double [][][] imageData, // [img]{Alpha, Y,Cb,Cr, Ext}. Alpha may be null - will handle later?
    				int imageFullWidth,
//    				double blurVarianceSigma,
//    				int refineTilePeriod,
    				double refinePhaseCoeff,
    				double refineHighPassSigma,
    				double refineLowPassSigma,
    				double refineCorrMaxDistance,
    				double refineCorrThreshold,
    				int refineSubPixel,
//    				double zMapMinForeground,
//    				int zMapVarMask,
//    				double [] zMapVarThresholds,
//    				double [] zMapVarWeights,
//    				int       auxVarMode,
//    				int     normalize,
    				int zMapCorrMask,
    				double [] zMapCorrThresholds,
    				double [] zMapCorrWeights,
    				double [] window,
    				DoubleFHT doubleFHT,
//    				double [] refineWindow,
//    				DoubleFHT refineFHT,
    				// new arguments
//    				int combineMode, // different image pairs - 0
    				double disparityMax,
    				double disparityMin,
    				double minAbsolute, // or NaN - will use enabled/disabled state of the tile
    				double minRelative,
    				boolean filterByForeground, // apply known certain masks
    				double filterByForegroundMargin,
    				boolean filterByDisabled,
    				double disparityTolearnce,
    				double maskBlurSigma,
    				double corrHighPassSigma, // subtract blurred version to minimize correlation caused by masks
    				int threadsMax,
    				boolean showProgress,
    				int debugLevel){
       			// TODO: also calculate "unlikely" - high autocorrelation, not occluded, low inter-correlation

    			double [] normCorrWeights=normalizeWeights(zMapCorrMask,zMapCorrWeights);
    			ZTile zTile=allZMap[nImg][tileY][tileX];
    			int tileOverlap=zTile.getOverlap();
    			int paddedSize=zTile.getPaddedSize();
    			int paddedLength=zTile.getPaddedLength();
    			int size=this.overlapStep*2;
    			int length=size*size;
				double[] zeros=new double [length];
				for (int i=0;i<length;i++) zeros[i]=0.0;

    			int margin=size/4-tileOverlap;
    			double [][][][] slices=new double [sImgSet.length][][][];
    			// Iterate through pairsd to this image - oter pairs will be processed separately and the result (likelyhood of belonguing to
    			// the particular plane can be evaluated  (they use the same "radar" data
    			zTile.setMinCorrelations(minAbsolute, minRelative, true); // Set here or before? include FG
    			if (debugLevel>3) {
    				System.out.println ("zTile.setMinCorrelations("+minAbsolute+","+ minRelative+")");
    				boolean [] enabled=zTile.enabledPlane;
    				for (int plane=0;plane<zTile.getNumberOfPlanes();plane++) {
    					System.out.println(plane+": disparity="+zTile.getPlaneDisparity(plane)+" strength="+zTile.getPlaneStrength(plane)+" enabled="+
    							enabled[plane]);

    				}
					System.out.println("zTile.minAbsolute="+zTile.minAbsolute+" zTile.minRelative="+zTile.minRelative);
    			}
        		double corrMaxDist2=refineCorrMaxDistance*refineCorrMaxDistance;
        		int halfDispRange=(int) Math.ceil(2*refineCorrMaxDistance*refineSubPixel);

        		if (Double.isNaN(window[0])) {
        			int index=0;
            		int quarterSize=size/4; //8
            		int halfSize=size/2; // 16
            		int size34=3*size/4; // 24
        			double[] window1d=doubleFHT.getHamming1d(halfSize); // Hamming
        			for (int iy=0;iy<size;iy++) {
        				double wy=(iy<quarterSize)?window1d[iy]:((iy>size34)?window1d[iy-halfSize]:1.0);
        				for (int ix=0;ix<size;ix++) {
            				double wx=(ix<quarterSize)?window1d[ix]:((ix>size34)?window1d[ix-halfSize]:1.0);
        					window[index++]=wx*wy;
        				}
        			}
					if (debugLevel>2){ // one per thread
						(new ShowDoubleFloatArrays()).showArrays(
								window,
								size,
								size,
								"Window_N"+nImg+"-X"+tileX+"-Y"+tileY
						);
					}

        		}
        		for (int plane=0;plane<zTile.getNumberOfPlanes();plane++) if (zTile.getPlaneEnabled(plane)) {
        			if (debugLevel>3) {
        				System.out.println ("refinePlaneDisparityTile("+tileX+","+tileY+", ...) plane="+plane+
        						" disparity="+zTile.getPlaneDisparity(plane)+" strength="+zTile.getPlaneStrength(plane));


        			}
        			if (!(zTile.getPlaneDisparity(plane)<disparityMin) && !(zTile.getPlaneDisparity(plane)>disparityMax)){ // NaN is OK
        				double disparity= zTile.getPlaneDisparity(plane);
        				if (debugLevel>3) {
//            				double testDisparityError=0.5;
        					System.out.println (" processing plane "+plane+" nImg="+nImg+
        							" disparity="+disparity+
        							" filterByForegroundMargin="+filterByForegroundMargin);
//        					disparity+=testDisparityError;
//        					System.out.println ("**** modified disparity for testing, new disparity="+disparity);
        				}
        				double [] thisEnabledMask;
        				if (filterByForeground){
        					boolean [] bMask=zTile.getEnabledNonOccluded(
        							disparity+filterByForegroundMargin,
        							disparity,
        							filterByDisabled?disparityTolearnce:Double.NaN,
        							debugLevel);
        					if (debugLevel>5){
        						(new ShowDoubleFloatArrays()).showArrays(
        								bMask,
        								paddedSize,
        								paddedSize,
        								"bMask_N"+nImg+"-X"+tileX+"-Y"+tileY
        						);
        					}

        					thisEnabledMask=zeros.clone();
        					for (int i=0;i<paddedLength;i++){ //tileOverlap
        						int i1= (margin+(i/paddedSize))*size +(margin+ (i%paddedSize));
        						thisEnabledMask[i1]=bMask[i]?1.0:0.0;
        					}
        					if (maskBlurSigma>0.0){
        						(new DoubleGaussianBlur()).blurDouble(
        								thisEnabledMask,
        								size,
        								size,
        								maskBlurSigma,
        								maskBlurSigma,
        								0.01);
        					}
        					if (debugLevel>5){
        						(new ShowDoubleFloatArrays()).showArrays(
        								thisEnabledMask,
        								size,
        								size,
        								"thisEnabledMask_N"+nImg+"-X"+tileX+"-Y"+tileY
        						);
        					}
        				} else {
        					thisEnabledMask=null;
        				}
        				double [][] pairEnabledMask=new double [sImgSet.length][];
        				double [] disparityArray=null;
        				//    				int disparityRange=(int) Math.ceil(refineCorrMaxDistance*refineSubPixel*2);
        				for (int sIndex=0;sIndex<sImgSet.length;sIndex++){
        					int sImg=sImgSet[sIndex];
        					double dX=disparity*(this.disparityScales[sImg][0]-this.disparityScales[nImg][0]);
        					double dY=disparity*(this.disparityScales[sImg][1]-this.disparityScales[nImg][1]);
        					double dXYlen=Math.sqrt(dX*dX+dY*dY);
        					double [] udXY={dX/dXYlen,dY/dXYlen};
        					if (debugLevel>3){ // +2 for selected tile
        						System.out.println("Iterating tileX="+tileX+", tileY="+tileY+" disparity="+disparity+" dX="+ dX+" dY="+dY+
        								" sIndex="+sIndex+" nImg="+nImg+" sImg="+sImg);
        					}
        					if (filterByForeground){
        						double [] otherNonOccluded=	getNonOccluded(
        								allZMap[sImg],//ZTile [][] thisZMap,
        								tileX*this.overlapStep+this.overlapStep/2+dX, //double xc,
        								tileY*this.overlapStep+this.overlapStep/2+dY, //double yc,
        								disparity, //double disparity,
        								filterByForegroundMargin,
        								filterByDisabled?disparityTolearnce:Double.NaN,
        								debugLevel);
            					if (debugLevel>5){
            						(new ShowDoubleFloatArrays()).showArrays(
            								otherNonOccluded,
            								paddedSize,
            								paddedSize,
            								"PEM_N"+nImg+"-S"+sImg+"-X"+tileX+"-Y"+tileY
            						);
            					}

        						pairEnabledMask[sIndex]=zeros.clone();
        						for (int i=0;i<paddedLength;i++){
        							int i1= (margin+(i/paddedSize))*size +(margin+ (i%paddedSize));
        							pairEnabledMask[sIndex][i1]=otherNonOccluded[i];
        						}

        						if (maskBlurSigma>0.0){
        							(new DoubleGaussianBlur()).blurDouble(
        									pairEnabledMask[sIndex],
        									size,
        									size,
        									maskBlurSigma,
        									maskBlurSigma,
        									0.01);
        						}
            					if (debugLevel>5){
            						(new ShowDoubleFloatArrays()).showArrays(
            								pairEnabledMask[sIndex],
            								size,
            								size,
            								"pairEnabledMask_N"+nImg+"-S"+sImg+"-X"+tileX+"-Y"+tileY
            						);
            					}
        						//multiply masks from both images
        						for (int i=0;i<pairEnabledMask[sIndex].length;i++) pairEnabledMask[sIndex][i]*=thisEnabledMask[i];
        					} else {
        						pairEnabledMask[sIndex]=window.clone(); // should be initialized to "sharp"?
        					}
        					double maskTotal=0.0;
        					for (int i=-0;i<pairEnabledMask[sIndex].length;i++) maskTotal+=pairEnabledMask[sIndex][i];
        					slices[sIndex]=getShiftedSlices(
        							tileX*this.overlapStep+this.overlapStep/2,
        							tileY*this.overlapStep+this.overlapStep/2,
        							this.overlapStep,
        							nImg, //int first,
        							sImg, //int second,
        							dX, //double dxA,
        							dY, //double dyA,
        							zMapCorrMask, // used at least used in one filter,
        							imageData,
        							imageFullWidth,
        							true, //doubleSizeOutput,
        							window, // double [] window, // should be 4*size*size long;
        							doubleFHT, //DoubleFHT doubleFHT, // to reuse tables
        							debugLevel);
        					if (debugLevel>3){ // +2 for selected tile
        						System.out.println("Debugging tileX="+tileX+", tileY="+tileY+" plane="+plane+" disparity="+disparity+" dX="+ dX+" dY="+dY+
        								" sIndex="+sIndex+" nImg="+nImg+" sImg="+sImg);
        						int numActive=0;
        						for (int i=0;i<slices[sIndex].length;i++) if (slices[sIndex][i]!=null) numActive++;
        						String [] channelNames={"alpha","Y","Cb","Cr","Aux"};
        						double [][] debugData=new double [2*numActive][];
        						String [] debugTitles=new String [2*numActive];
        						int index=0;
        						for (int i=0;i<slices[sIndex].length;i++) if (slices[sIndex][i]!=null) {
        							debugTitles[index]=channelNames[i]+"-"+nImg;
        							debugData[index++]=slices[sIndex][i][0];
        							debugTitles[index]=channelNames[i]+"-"+sImg;
        							debugData[index++]=slices[sIndex][i][1];
        						}
        						(new ShowDoubleFloatArrays()).showArrays(
        								debugData,
        								size,
        								size,
        								true,
        								"S"+plane+"-"+IJ.d2s(disparity,1)+"_N"+nImg+"-"+sImg+"_sIndex-"+sIndex+"_tile"+tileX+"-"+tileY,
        								debugTitles);
        					}

        					// process correlation, calculate auto-first, auto-second and cross for each pair,

        					// add them with weights and calculate disparity correction - at this stage - only for the whole tile
        					double [] centerCorr={0.0,0.0,0.0}; //00, 11, 01

        					double [][][] pairs=new double [slices[sIndex].length-1][][];
        					for (int chn=0;chn<(slices[sIndex].length-1);chn++) if ((zMapCorrMask& (1<<chn))!=0){
        						//    						int length=slices[sIndex][chn+1][0].length;
        						//    						size=(int) Math.sqrt(length)
        						pairs[chn]=new double[2][length]; //slices[sIndex][chn+1]; // skip alpha
        						for (int n=0;n<2;n++) {
        							pairs[chn][n]=slices[sIndex][chn+1][n].clone();
        							//high-pass each image in the pair (to reduce mask influence
        							if (corrHighPassSigma>0.0){
        								double [] loPass=pairs[chn][n].clone();
        								(new DoubleGaussianBlur()).blurDouble(
        										loPass,
        										size,
        										size,
        										corrHighPassSigma,
        										corrHighPassSigma,
        										0.01);
        								for (int i=0;i<length;i++) pairs[chn][n][i]-=loPass[i];
        							} else {
        								normalizeAndWindow (pairs[chn][n], null, true); // only remove DC
        							}
        							for (int i=0;i<length;i++) pairs[chn][n][i]*=pairEnabledMask[sIndex][i]; //
        						}
        						for (int i=0;i<length;i++){
        							centerCorr[0]+=normCorrWeights[chn]*pairs[chn][0][i]*pairs[chn][0][i];
        							centerCorr[1]+=normCorrWeights[chn]*pairs[chn][1][i]*pairs[chn][1][i];
        							centerCorr[2]+=normCorrWeights[chn]*pairs[chn][0][i]*pairs[chn][1][i];
        						}
        					}
        					if (debugLevel>3){ // +2 for selected tile
        						System.out.println("centerCorr[0]="+centerCorr[0]+" centerCorr[1]="+centerCorr[1]+" centerCorr[2]="+centerCorr[2]);
// show high-pass/masked pairs
        						int numActive=0;
        						for (int i=0;i<slices[sIndex].length;i++) if (slices[sIndex][i]!=null) numActive++;
        						String [] channelNames={"alpha","Y","Cb","Cr","Aux"};
        						double [][] debugData=new double [2*numActive+1][];
        						String [] debugTitles=new String [2*numActive+1];
        						int index=0;
        						debugTitles[index]="mask";
    							debugData[index++]=thisEnabledMask;
        						for (int i=0;i<slices[sIndex].length;i++) if (slices[sIndex][i]!=null) {
        							debugTitles[index]=channelNames[i]+"-"+nImg;
        							if ((index>=debugData.length) || ((i-1)>=pairs.length)){
        								System.out.println(
        										" index="+index+
        										" debugData.length="+debugData.length+
        										" i="+i+
        										" pairs.length="+pairs.length+
        										" slices.length="+slices.length+
        										" slices[sIndex].length="+slices[sIndex].length
        										);
        							}
        							debugData[index++]=pairs[i-1][0]; //ava.lang.ArrayIndexOutOfBoundsException: 3
        							debugTitles[index]=channelNames[i]+"-"+sImg;
        							debugData[index++]=pairs[i-1][1];
        						}
        						(new ShowDoubleFloatArrays()).showArrays(
        								debugData,
        								size,
        								size,
        								true,
        								"P"+plane+"_N"+nImg+"-"+sImg+"_sIndex-"+sIndex+"_tile"+tileX+"-"+tileY,
        								debugTitles);
        					}
        					// TODO: now look if the correlation is strong enough to perform correction (minimal of 3?)
        					double minCorr=1.0;
        					for (int i=0;i<centerCorr.length;i++) {
        						centerCorr[i]=Math.sqrt(centerCorr[i]/maskTotal);
        						minCorr*=centerCorr[i];
        					}
        					minCorr=Math.pow(minCorr,1.0/centerCorr.length); // Use other metrics?
        					//        				double refineCorrMaxDistance,
        					//        				double refineCorrThreshold,
        					if (debugLevel>3){ // +2 for selected tile
        						System.out.println("centerCorr[0]="+centerCorr[0]+" centerCorr[1]="+centerCorr[1]+" centerCorr[2]="+centerCorr[2]+
        								 " minCorr="+minCorr+" refineCorrThreshold="+refineCorrThreshold);
        					}

        					if (minCorr>=refineCorrThreshold){
        						double [] combinedChnCorr=new double [length];
        						for (int i=0;i<combinedChnCorr.length;i++) combinedChnCorr[i]=0.0;
        						for (int chn=0;chn<(slices[sIndex].length-1);chn++) if ((zMapCorrMask& (1<<chn))!=0){
        							double[]  corr=doubleFHT.correlate (
        									pairs[chn][1],
        									pairs[chn][0],
        									refineHighPassSigma,
        									refineLowPassSigma,
        									refinePhaseCoeff); // second will be modified
        							for (int i=0;i<combinedChnCorr.length;i++) combinedChnCorr[i]+=normCorrWeights[chn]*corr[i];
        						}
        						// See if maximum is close enough, then subpixel and calculate partial disparityArray - to be combined into a single one later
        						double max=0;
        						int iMax=0;
        						for (int i=0;i<length;i++) if (combinedChnCorr[i]>max){
        							max=combinedChnCorr[i];
        							iMax=i;
        						}
        						int ixc=iMax%size-size/2;
        						int iyc=iMax/size-size/2;
        						if (debugLevel>3) System.out.println("refineCorrelation(): max="+max+" iMax="+iMax+" refineCorrMaxDistance="+refineCorrMaxDistance+" corrMaxDist2="+corrMaxDist2+
        								" r2="+(ixc*ixc+iyc*iyc)+" r="+Math.sqrt(ixc*ixc+iyc*iyc));
        						if ((ixc*ixc+iyc*iyc)<=corrMaxDist2){ // maximum close enough
        							double [] upsampled=doubleFHT.upsample(combinedChnCorr,refineSubPixel);
        							int interpolatedSize=size*refineSubPixel;
        							int interpolatedCenter=(interpolatedSize+1)*interpolatedSize/2;
        							if (debugLevel>3) {
                						(new ShowDoubleFloatArrays()).showArrays(
                								upsampled,
                								interpolatedSize,
                								interpolatedSize,
                								"U"+plane+"_N"+nImg+"-"+sImg+"_sIndex-"+sIndex+"_tile"+tileX+"-"+tileY);

        							}
        							if (disparityArray==null) {
        								disparityArray=new double [2*halfDispRange+1];
        								for (int i=0;i<disparityArray.length;i++) disparityArray[i]=0.0;
        							}
        							for (int i=0;i<2*halfDispRange+1;i++){
        								double deltaX=udXY[0]*(i-halfDispRange);
        								double deltaY=udXY[1]*(i-halfDispRange);
        								int iX=(int) Math.floor(deltaX); // zero in the center
        								int iY=(int) Math.floor(deltaY);
        								deltaX-=iX;
        								deltaY-=iY;
        								int index00= interpolatedCenter+iY*interpolatedSize+iX;
        								int index01=index00+interpolatedSize;
        								int index10=index00+1;
        								int index11=index01+1;
        								disparityArray[i]+= // ACCUMULATE bi-linear interpolated data
        									(upsampled[index00]*(1.0-deltaX)+upsampled[index10]*deltaX)*(1.0-deltaY)+
        									(upsampled[index01]*(1.0-deltaX)+upsampled[index11]*deltaX)*     deltaY;
        								if (debugLevel>5){
        									System.out.println("disparityArray["+i+"]="+disparityArray[i]+
        											" deltaX="+deltaX+" deltaY="+deltaY+" iX="+iX+" iY="+iY+" index00="+index00+" index01="+index01+
        											" index10="+index10+" index11="+index11);
        								}
        							}
        						}

        					} // if (minCorr>=refineCorrThreshold)
        				} //for (int sIndex=0;sIndex<sImgSet.length;sIndex++){
    					if (disparityArray!=null){
    						int iMax=0;
    						for (int i=1;i<disparityArray.length;i++){
    							if (disparityArray[i]>disparityArray[iMax]) iMax=i;
    						}
    						double disparCorr=((double) (iMax-halfDispRange))/refineSubPixel;
    						if (debugLevel>5){
    							for (int i=0;i<disparityArray.length;i++){
    								System.out.println("combined_disparityArray["+i+"]="+disparityArray[i]);
    							}
    							System.out.println("\nDisparity correction="+disparCorr+" (iMax="+iMax+"), limit ="+refineCorrMaxDistance +" - VERIFY SIGN IS CORRECT!");
    						}
    						if (disparCorr<=refineCorrMaxDistance){
    							if (debugLevel>3){
    								System.out.println("Old disparity["+plane+"]="+disparity+" new disparity="+(disparity+disparCorr));
    							}
    							zTile.setPlaneDisparity(disparity+disparCorr,plane);
    						}
    					}

        			}
    			} // end for (int plane...)
       		}



    		public void planeLikely(
//    				int nImg,
//    				int [] sImg, // list of second images
    				double [][][] imageData, // [img]{Alpha, Y,Cb,Cr, Ext}. Alpha may be null - will handle later?
    				int imageFullWidth,
    				double blurVarianceSigma,
    				int refineTilePeriod,
    				double subTilePhaseCoeff,
    				double subTileHighPassSigma,
    				double subTileLowPassSigma,
//    				double refineCorrMaxDistance,
//    				double refineCorrThreshold,
    				int refineSubPixel,
    				double zMapMinForeground,
    				int zMapVarMask,
    				double [] zMapVarThresholds,
    				double [] zMapVarWeights,
    				int       auxVarMode,
        			int     normalizeCorrelation,
    				int zMapCorrMask,
    				double [] zMapCorrThresholds,
    				double [] zMapCorrThresholdsRel,
    				double [] zMapCorrWeights,
    				int debugRow,
    				int debugColumn,
    				int combineMode, // different image pairs - 0

    				double disparityMax,
    				double disparityMin,
    				double minAbsolute, // or NaN - will use enabled/disabled state of the tile
    				double minRelative,
    				boolean filterByForeground, // apply known certain masks
    				double filterByForegroundMargin,
    				boolean filterByDisabled,
    				double disparityTolearnce,
    				double maskBlurSigma,
    				double corrHighPassSigma, // subtract blurred version to minimize correlation caused by masks
    				double varianceBlurScale,
    				double kLocal,
    				int matchStatMode,
    				int threadsMax,
    				boolean showProgress,
    				int debugLevel){
    			// process all image pairs
    			for (int nImg=0;nImg<this.disparityScales.length;nImg++){
    				int [] sImg=new int [this.disparityScales.length-1];
    				for (int j=0;j<sImg.length;j++)sImg[j]=(j<nImg)?j:(j+1);
    				planeLikely(
    	    				nImg,
    	    				sImg, // list of second images
    	    				imageData, // [img]{Alpha, Y,Cb,Cr, Ext}. Alpha may be null - will handle later?
    	    				imageFullWidth,
    	    				blurVarianceSigma,
    	    				refineTilePeriod,
    	    				subTilePhaseCoeff,
    	    				subTileHighPassSigma,
    	    				subTileLowPassSigma,
//    	    				refineCorrMaxDistance,
//    	    				refineCorrThreshold,
    	    				refineSubPixel,
    	    				zMapMinForeground,
    	    				zMapVarMask,
    	    				zMapVarThresholds,
    	    				zMapVarWeights,
    	    				auxVarMode,
    	        			normalizeCorrelation,
    	    				zMapCorrMask,
    	    				zMapCorrThresholds,
    	    				zMapCorrThresholdsRel,
    	    				zMapCorrWeights,
    	    				debugRow,
    	    				debugColumn,
    	    				combineMode, // different image pairs - 0

    	    				disparityMax,
    	    				disparityMin,
    	    				minAbsolute, // or NaN - will use enabled/disabled state of the tile
    	    				minRelative,
    	    				filterByForeground, // apply known certain masks
    	    				filterByForegroundMargin,
    	    				filterByDisabled,
    	    				disparityTolearnce,
    	    				maskBlurSigma,
    	    				corrHighPassSigma, // subtract blurred version to minimize correlation caused by masks
    	    				varianceBlurScale,
    	    				kLocal,
    	    				matchStatMode,
    	    				threadsMax,
    	    				showProgress,
    	    				debugLevel);
    			}
    		}

        	public void planeLikely(
    				final int nImg,
    				final int [] sImg, // list of second images
    				final double [][][] imageData, // [img]{Alpha, Y,Cb,Cr, Ext}. Alpha may be null - will handle later?
    				final int imageFullWidth,
    				final double blurVarianceSigma,
    				final int refineTilePeriod,
    				final double subTilePhaseCoeff,
    				final double subTileHighPassSigma,
    				final double subTileLowPassSigma,
//    				final double refineCorrMaxDistance,
//    				final double refineCorrThreshold,
    				final int refineSubPixel,
    				final double zMapMinForeground,
    				final int zMapVarMask,
    				final double [] zMapVarThresholds,
    				final double [] zMapVarWeights,
    				final int       auxVarMode,
        			final int     normalizeCorrelation,
    				final int zMapCorrMask,
    				final double [] zMapCorrThresholds,
    				final double [] zMapCorrThresholdsRel,
    				final double [] zMapCorrWeights,
    				final int debugRow,
    				final int debugColumn,
    				final int combineMode, // different image pairs - 0

    				final double disparityMax,
    				final double disparityMin,
    				final double minAbsolute, // or NaN - will use enabled/disabled state of the tile
    				final double minRelative,
    				final boolean filterByForeground, // apply known certain masks
    				final double filterByForegroundMargin,
    				final boolean filterByDisabled,
    				final double disparityTolearnce,
    				final double maskBlurSigma,
    				final double corrHighPassSigma, // subtract blurred version to minimize correlation caused by masks
    				final double varianceBlurScale,
    				final double kLocal,
    				final int matchStatMode,
    				final int threadsMax,
    				final boolean showProgress,
    				final int debugLevel){
        		final int debugThreshold=2;
        		final Rectangle zMapWOI=this.zMapWOI;
        		final ZTile [][][] zMapFinal=this.zMap;
        		final int tilesX=this.zMapWOI.width;
        		final int tiles=this.zMapWOI.width*this.zMapWOI.height;
        		if (debugLevel>2) System.out.println("setupZMap() zMapWOI.x="+zMapWOI.x+" zMapWOI.y="+zMapWOI.y+
        				" zMapWOI.width="+zMapWOI.width+" zMapWOI.height="+zMapWOI.height);
        		final Thread[] threads = newThreadArray(threadsMax);
		   		final AtomicInteger tileIndexAtomic     = new AtomicInteger(0);
		   		if (showProgress) IJ.showProgress(0.0);
		   		if (showProgress) IJ.showStatus("planeLikely for image "+nImg+" ("+(nImg+1)+" of "+this.disparityScales.length+") ...");
				final int debugTile=debugRow*tilesX+debugColumn;
				//TODO: define here
//				final double [] window=new double [4*this.overlapStep*this.overlapStep];
//    			window[0]=Double.NaN;
//				final double [] refineWindow=new double [refineTilePeriod*refineTilePeriod*4];
//				refineWindow[0]=Double.NaN;
				final int size=2*this.overlapStep;
		   		for (int ithread = 0; ithread < threads.length; ithread++) {
		   			threads[ithread] = new Thread() {
		   				@Override
						public void run() {
		    				DoubleFHT doubleFHT = new DoubleFHT();
		    				final double [] window=new double [size*size];
		        			window[0]=Double.NaN;
		        			final double [] subTileWindow=new double [refineTilePeriod*refineTilePeriod*4];
		        			subTileWindow[0]=Double.NaN;
		    				DoubleFHT subTileFHT= new DoubleFHT();
		   					for (int tile=tileIndexAtomic.getAndIncrement(); tile<tiles;tile=tileIndexAtomic.getAndIncrement()){
		   						boolean debugThis= (debugLevel>=debugThreshold) && (tile==debugTile);
		   						int thisDebugLevel=debugLevel+(debugThis?2:0);
		   						planeLikelyTile (
		   								zMapWOI.x+tile%tilesX, //int tileX,
		   								zMapWOI.y+tile/tilesX, //int tileY,
		   								zMapFinal,
		   			    				nImg,
		   			    				sImg, // list of second images
		   			    				imageData, // [img]{Alpha, Y,Cb,Cr, Ext}. Alpha may be null - will handle later?
		   			    				imageFullWidth,
		   			    				blurVarianceSigma,
//		   			    				refineTilePeriod,
		   			    				subTilePhaseCoeff,
		   			    				subTileHighPassSigma,
		   			    				subTileLowPassSigma,
//		   			    				refineCorrMaxDistance,
//		   			    				refineCorrThreshold,
		   			    				refineSubPixel,
		   			    				zMapMinForeground,
		   			    				zMapVarMask,
		   			    				zMapVarThresholds,
		   			    				zMapVarWeights,
		   			    				auxVarMode,
		   			    				normalizeCorrelation,
		   			    				zMapCorrMask,
		   			    				zMapCorrThresholds,
		   			    				zMapCorrThresholdsRel,
		   			    				zMapCorrWeights,
		   			    				window,
		   			    				doubleFHT,
		   			    				subTileWindow,
		   			    				subTileFHT,
		   			    				combineMode, // different image pairs - 0
		   			    				disparityMax,
		   			    				disparityMin,
		   			    				minAbsolute, // or NaN - will use enabled/disabled state of the tile
		   			    				minRelative,
		   			    				filterByForeground, // apply known certain masks
		   			    				filterByForegroundMargin,
		   			    				filterByDisabled,
		   			    				disparityTolearnce,
		   			    				maskBlurSigma,
		   			    				corrHighPassSigma, // subtract blurred version to minimize correlation caused by masks
		   			    				varianceBlurScale,
		   			    				kLocal,
		   			    				matchStatMode,
		   			    				threadsMax,
		   			    				showProgress,
		   			    				thisDebugLevel);
		   						if (showProgress){
		   							final int finalTile=tile;
		   							SwingUtilities.invokeLater(new Runnable() {
		   								@Override
										public void run() {
		   									IJ.showProgress(finalTile,tiles);
		   								}
		   							});
		   						}
		   					}
		   				}
		   			};
		   		}
		   		startAndJoin(threads);
		   		IJ.showProgress(1.0);
        	}




       		private void planeLikelyTile ( // false if nothing left in the current foreground, may repeat
    				int tileX,
    				int tileY,
    				ZTile [][][] allZMap,
    				int nImg,
    				int [] sImgSet, // list of second images
    				double [][][] imageData, // [img]{Alpha, Y,Cb,Cr, Ext}. Alpha may be null - will handle later?
    				int imageFullWidth,
    				double blurVarianceSigma,
//    				int refineTilePeriod,
    				double subTilePhaseCoeff,
    				double subTileHighPassSigma,
    				double subTileLowPassSigma,
//    				double refineCorrMaxDistance,
//    				double refineCorrThreshold,
    				int refineSubPixel,
    				double zMapMinForeground,
    				int zMapVarMask,
    				double [] zMapVarThresholds,
    				double [] zMapVarWeights,
    				int       auxVarMode,
    				int     normalize,
    				int zMapCorrMask,
    				double [] zMapCorrThresholds,
    				double [] zMapCorrThresholdsRel,
    				double [] zMapCorrWeights,
    				double [] window,
    				DoubleFHT doubleFHT,
    				double [] subTileWindow,
    				DoubleFHT subTileFHT,
    				// new arguments
    				int combineMode, // different image pairs - 0
    				double disparityMax,
    				double disparityMin,
    				double minAbsolute, // or NaN - will use enabled/disabled state of the tile
    				double minRelative,
    				boolean filterByForeground, // apply known certain masks
    				double filterByForegroundMargin,
    				boolean filterByDisabled,
    				double disparityTolearnce,
    				double maskBlurSigma,
    				double corrHighPassSigma, // subtract blurred version to minimize correlation caused by masks
    				//tone-matching statistical parameters
    				double varianceBlurScale,
    				double kLocal,
    				int matchStatMode,
    				int threadsMax,
    				boolean showProgress,
    				int debugLevel){
       			// TODO: also calculate "unlikely" - high autocorrelation, not occluded, low inter-correlation

    			double [] normVarWeights= normalizeWeights(zMapVarMask,zMapVarWeights);
    			double [] normCorrWeights=normalizeWeights(zMapCorrMask,zMapCorrWeights);
//                int auxChannelNumber=3;
    			ZTile zTile=allZMap[nImg][tileY][tileX];
    			int tileOverlap=zTile.getOverlap();
    			int paddedSize=zTile.getPaddedSize();
    			int paddedLength=zTile.getPaddedLength();
    			int size=this.overlapStep*2;
    			int length=size*size;
				double[] zeros=new double [length];
				for (int i=0;i<length;i++) zeros[i]=0.0;
    			int margin=size/4-tileOverlap;
    			double [][][][] slices=new double [sImgSet.length][][][];
    			int [] dirs1={1,size+1,size,size-1,-1,-size-1,-size,-size+1,0};
//    			double [][][][] variance=new double [sImgSet.length][][][]; // [sIndex][pair 0/1][chn][pixel
//    			int length=4*this.overlapStep*this.overlapStep; // initialize it to a correct value right here?
    			int [] borderMask=zTile.getBorderMask();
    			// Iterate through pairsd to this image - oter pairs will be processed separately and the result (likelyhood of belonguing to
    			// the particular plane can be evaluated  (they use the same "radar" data
    			zTile.setMinCorrelations(minAbsolute, minRelative, true); // Set here or before? include FG
    			if (debugLevel>3) {
    				System.out.println ("zTile.setMinCorrelations("+minAbsolute+","+ minRelative+")");
    				boolean [] enabled=zTile.enabledPlane;
    				for (int plane=0;plane<zTile.getNumberOfPlanes();plane++) {
    					System.out.println(plane+": disparity="+zTile.getPlaneDisparity(plane)+" strength="+zTile.getPlaneStrength(plane)+" enabled="+
    							enabled[plane]);

    				}
					System.out.println("zTile.minAbsolute="+zTile.minAbsolute+" zTile.minRelative="+zTile.minRelative);
    			}
        		if (Double.isNaN(window[0])) {
        			int index=0;
            		int quarterSize=size/4; //8
            		int halfSize=size/2; // 16
            		int size34=3*size/4; // 24
        			double[] window1d=doubleFHT.getHamming1d(halfSize); // Hamming
        			for (int iy=0;iy<size;iy++) {
        				double wy=(iy<quarterSize)?window1d[iy]:((iy>size34)?window1d[iy-halfSize]:1.0);
        				for (int ix=0;ix<size;ix++) {
            				double wx=(ix<quarterSize)?window1d[ix]:((ix>size34)?window1d[ix-halfSize]:1.0);
        					window[index++]=wx*wy;
        				}
        			}
					if (debugLevel>2){ // one per thread
						(new ShowDoubleFloatArrays()).showArrays(
								window,
								size,
								size,
								"Window_N"+nImg+"-X"+tileX+"-Y"+tileY
						);
					}
        		}

				double [][][] planeStrength=new double [zTile.getNumberOfPlanes()][sImgSet.length][];
				for (int plane=0;plane<zTile.getNumberOfPlanes();plane++)  planeStrength[plane]=null;

        		for (int plane=0;plane<zTile.getNumberOfPlanes();plane++) if (zTile.getPlaneEnabled(plane)) {
        			if (debugLevel>3) {
        				System.out.println ("planeLikelyTile("+tileX+","+tileY+", ...) plane="+plane+
        						" disparity="+zTile.getPlaneDisparity(plane)+" strength="+zTile.getPlaneStrength(plane));
        			}

        			if (!(zTile.getPlaneDisparity(plane)<disparityMin) && !(zTile.getPlaneDisparity(plane)>disparityMax)){ // NaN is OK
        				planeStrength[plane]=new double [sImgSet.length][];
        				double disparity= zTile.getPlaneDisparity(plane);
        				if (debugLevel>3) {
        					System.out.println (" processing plane "+plane+" nImg="+nImg+
        							" disparity="+disparity+
        							" filterByForegroundMargin="+filterByForegroundMargin);
        				}
        				double [] thisEnabledMask;
        				if (filterByForeground){
        					boolean [] bMask=zTile.getEnabledNonOccluded(
        							disparity+filterByForegroundMargin,
        							disparity,
        							filterByDisabled?disparityTolearnce:Double.NaN,
        							debugLevel);
        					if (debugLevel>5){
        						(new ShowDoubleFloatArrays()).showArrays(
        								bMask,
        								paddedSize,
        								paddedSize,
        								"bMask_N"+nImg+"-X"+tileX+"-Y"+tileY
        						);
        					}

        					thisEnabledMask=zeros.clone();
        					for (int i=0;i<paddedLength;i++){ //tileOverlap
        						int i1= (margin+(i/paddedSize))*size +(margin+ (i%paddedSize));
        						thisEnabledMask[i1]=bMask[i]?1.0:0.0;
        					}
        					if (maskBlurSigma>0.0){
        						(new DoubleGaussianBlur()).blurDouble(
        								thisEnabledMask,
        								size,
        								size,
        								maskBlurSigma,
        								maskBlurSigma,
        								0.01);
        					}
        					if (debugLevel>5){
        						(new ShowDoubleFloatArrays()).showArrays(
        								thisEnabledMask,
        								size,
        								size,
        								"thisEnabledMask_N"+nImg+"-X"+tileX+"-Y"+tileY
        						);
        					}
        				} else {
        					thisEnabledMask=null;
        				}
        				double [][] pairEnabledMask=new double [sImgSet.length][];
//        				double [] disparityArray=null; // needed??
        				for (int sIndex=0;sIndex<sImgSet.length;sIndex++){
        					int sImg=sImgSet[sIndex];
        					double dX=disparity*(this.disparityScales[sImg][0]-this.disparityScales[nImg][0]);
        					double dY=disparity*(this.disparityScales[sImg][1]-this.disparityScales[nImg][1]);
//        					double dXYlen=Math.sqrt(dX*dX+dY*dY);
//        					double [] udXY={dX/dXYlen,dY/dXYlen};
        					if (debugLevel>3){ // +2 for selected tile
        						System.out.println("Iterating tileX="+tileX+", tileY="+tileY+" disparity="+disparity+" dX="+ dX+" dY="+dY+
        								" sIndex="+sIndex+" nImg="+nImg+" sImg="+sImg);
        					}
        					if (filterByForeground){
        						double [] otherNonOccluded=	getNonOccluded(
        								allZMap[sImg],//ZTile [][] thisZMap,
        								tileX*this.overlapStep+this.overlapStep/2+dX, //double xc,
        								tileY*this.overlapStep+this.overlapStep/2+dY, //double yc,
        								disparity, //double disparity,
        								filterByForegroundMargin,
        								filterByDisabled?disparityTolearnce:Double.NaN,
        								debugLevel);
            					if (debugLevel>5){
            						(new ShowDoubleFloatArrays()).showArrays(
            								otherNonOccluded,
            								paddedSize,
            								paddedSize,
            								"PEM_N"+nImg+"-S"+sImg+"-X"+tileX+"-Y"+tileY
            						);
            					}
        						pairEnabledMask[sIndex]=zeros.clone();
        						for (int i=0;i<paddedLength;i++){
        							int i1= (margin+(i/paddedSize))*size +(margin+ (i%paddedSize));
        							pairEnabledMask[sIndex][i1]=otherNonOccluded[i];
        						}

        						if (maskBlurSigma>0.0){
        							(new DoubleGaussianBlur()).blurDouble(
        									pairEnabledMask[sIndex],
        									size,
        									size,
        									maskBlurSigma,
        									maskBlurSigma,
        									0.01);
        						}
            					if (debugLevel>5){
            						(new ShowDoubleFloatArrays()).showArrays(
            								pairEnabledMask[sIndex],
            								size,
            								size,
            								"pairEnabledMask_N"+nImg+"-S"+sImg+"-X"+tileX+"-Y"+tileY
            						);
            					}
        						//multiply masks from both images
        						for (int i=0;i<pairEnabledMask[sIndex].length;i++) pairEnabledMask[sIndex][i]*=thisEnabledMask[i];
        					} else {
        						pairEnabledMask[sIndex]=window.clone(); // should be initialized to "sharp"?
        					}
        					slices[sIndex]=getShiftedSlices(
        							tileX*this.overlapStep+this.overlapStep/2,
        							tileY*this.overlapStep+this.overlapStep/2,
        							this.overlapStep,
        							nImg, //int first,
        							sImg, //int second,
        							dX, //double dxA,
        							dY, //double dyA,
        							zMapCorrMask | zMapVarMask, // used at least used in one filter,
        							imageData,
        							imageFullWidth,
        							true, //doubleSizeOutput,
        							window, // double [] window, // should be 4*size*size long;
        							doubleFHT, //DoubleFHT doubleFHT, // to reuse tables
        							debugLevel);
        					if (debugLevel>3){ // +2 for selected tile
        						System.out.println("Debugging tileX="+tileX+", tileY="+tileY+" plane="+plane+" disparity="+disparity+" dX="+ dX+" dY="+dY+
        								" sIndex="+sIndex+" nImg="+nImg+" sImg="+sImg);
        						int numActive=0;
        						for (int i=0;i<slices[sIndex].length;i++) if (slices[sIndex][i]!=null) numActive++;
        						String [] channelNames={"alpha","Y","Cb","Cr","Aux"};
        						double [][] debugData=new double [2*numActive][];
        						String [] debugTitles=new String [2*numActive];
        						int index=0;
        						for (int i=0;i<slices[sIndex].length;i++) if (slices[sIndex][i]!=null) {
        							debugTitles[index]=channelNames[i]+"-"+nImg;
        							debugData[index++]=slices[sIndex][i][0];
        							debugTitles[index]=channelNames[i]+"-"+sImg;
        							debugData[index++]=slices[sIndex][i][1];
        						}
        						(new ShowDoubleFloatArrays()).showArrays(
        								debugData,
        								size,
        								size,
        								true,
        								"S"+plane+"-"+IJ.d2s(disparity,1)+"_N"+nImg+"-"+sImg+"_sIndex-"+sIndex+"_tile"+tileX+"-"+tileY,
        								debugTitles);
        					}
        					// process correlation, calculate auto-first, auto-second and cross for each pair,
        					// add them with weights and calculate disparity correction - at this stage - only for the whole tile
        					double [][][] pairs=new double [slices[sIndex].length-1][][];
        					double [][][] corrs=null;
        					for (int chn=0;chn<(slices[sIndex].length-1);chn++) if ((zMapCorrMask& (1<<chn))!=0){ // should work with all channels disabled
        						//    						int length=slices[sIndex][chn+1][0].length;
        						//    						size=(int) Math.sqrt(length)
        						pairs[chn]=new double[2][length]; //slices[sIndex][chn+1]; // skip alpha
        						for (int n=0;n<2;n++) {
        							pairs[chn][n]=slices[sIndex][chn+1][n].clone();
        							//high-pass each image in the pair (to reduce mask influence
        							if (corrHighPassSigma>0.0){
        								double [] loPass=pairs[chn][n].clone();
        								(new DoubleGaussianBlur()).blurDouble(
        										loPass,
        										size,
        										size,
        										corrHighPassSigma,
        										corrHighPassSigma,
        										0.01);
        								for (int i=0;i<length;i++) pairs[chn][n][i]-=loPass[i];
        							} else {
        								normalizeAndWindow (pairs[chn][n], null, true); // only remove DC
        							}
        							for (int i=0;i<length;i++) pairs[chn][n][i]*=pairEnabledMask[sIndex][i]; //
        						}
        						double [][][] thisChnCorrs=subTileCorrelation(
        								subTileFHT, //DoubleFHT doubleFHT, // will generate if null, can be used to reuse initialization
        								pairs[chn], //double [][] data,   // difines size
        								subTileWindow, //double [] window, // defines tile size
        								subTilePhaseCoeff, //double phaseCoeff, // phase correlation fraction (0 - normal, 1.0 - pure phase)
        								subTileHighPassSigma, //double highPassSigma,
        								subTileLowPassSigma, //double lowPassSigma,
        			        			debugLevel);
        						if (corrs==null){
        							corrs=new double [thisChnCorrs.length][thisChnCorrs[0].length][3];
        							for (int tY=0;tY<corrs.length;tY++) for (int tX=0;tX<corrs.length;tX++) for(int n=0;n<3;n++){
        								corrs[tY][tX][n]=normCorrWeights[chn]*thisChnCorrs[tY][tX][n];
        							}
        						} else {
        							for (int tY=0;tY<corrs.length;tY++) for (int tX=0;tX<corrs.length;tX++) for(int n=0;n<3;n++){
        								corrs[tY][tX][n]+=normCorrWeights[chn]*thisChnCorrs[tY][tX][n];
        							}
        						}
        					}
        					if (debugLevel>3){ // +2 for selected tile
//        						System.out.println("centerCorr[0]="+centerCorr[0]+" centerCorr[1]="+centerCorr[1]+" centerCorr[2]="+centerCorr[2]);
// show high-pass/masked pairs
        						int numActive=0;
        						for (int i=0;i<slices[sIndex].length;i++) if (slices[sIndex][i]!=null) numActive++;
        						String [] channelNames={"alpha","Y","Cb","Cr","Aux"};
        						double [][] debugData=new double [2*numActive+1][];
        						String [] debugTitles=new String [2*numActive+1];
        						int index=0;
        						debugTitles[index]="mask";
    							debugData[index++]=thisEnabledMask;
        						for (int i=0;i<slices[sIndex].length;i++) if (slices[sIndex][i]!=null) {
        							debugTitles[index]=channelNames[i]+"-"+nImg;
        							if ((index>=debugData.length) || ((i-1)>=pairs.length)){
        								System.out.println(
        										" index="+index+
        										" debugData.length="+debugData.length+
        										" i="+i+
        										" pairs.length="+pairs.length+
        										" slices.length="+slices.length+
        										" slices[sIndex].length="+slices[sIndex].length
        										);
        							}
        							debugData[index++]=pairs[i-1][0]; //ava.lang.ArrayIndexOutOfBoundsException: 3 //Exception in thread "Thread-1033" java.lang.NullPointerException
        							debugTitles[index]=channelNames[i]+"-"+sImg;
        							debugData[index++]=pairs[i-1][1];
        						}
        						(new ShowDoubleFloatArrays()).showArrays(
        								debugData,
        								size,
        								size,
        								true,
        								"P"+plane+"_N"+nImg+"-"+sImg+"_sIndex-"+sIndex+"_tile"+tileX+"-"+tileY,
        								debugTitles);
        					}

//=================== several alternative modes. For likelyhood - find linear relation between both images using weights - proportional to cross-correlation
        					double [] linearMatchWeight=pairEnabledMask[sIndex].clone();
        					if (corrs!=null){ // and some mode?
        		        		int tileSize=(int) Math.sqrt(subTileWindow.length); // 8
        		        		int tileStep=tileSize/2; //4
        		        		int tileExtra=size-tileSize; // 24??
//        		        		int tileLength=tileSize*tileSize; // 64
//        		        		int tileCenter=(tileSize+1)*tileSize/2; // 36
        		        		int topLeft=(size+1)*(size/4-tileSize/4);
        		        		int numTilesRow=size/(2*tileStep); // 4
                				// accumulate result
        		        		double [] corrTile=new double [length];
        		        		for (int i=0;i<corrTile.length;i++) corrTile[i]=0.0;
        		        		for (int subTileY=0;subTileY<numTilesRow;subTileY++) for (int subTileX=0;subTileX<numTilesRow;subTileX++) {
        		        			int index=topLeft+ tileStep*(subTileY*size+ subTileX);
        		        			double corr=corrs[subTileY][subTileX][2];
        		        			double aCorr=Math.sqrt(corrs[subTileY][subTileX][0]*corrs[subTileY][subTileX][1]);
        		        			if ((aCorr<zMapCorrThresholds[0]) || (corr<aCorr*zMapCorrThresholdsRel[0])) corr=0.0; // use channel0 thresholds !!!
        		        			if (debugLevel>3){
        		        				System.out.println ("corrs["+subTileY+"]["+subTileX+"]={"+
        		        						corrs[subTileY][subTileX][0]+","+
        		        						corrs[subTileY][subTileX][1]+","+
        		        						corrs[subTileY][subTileX][2]+"} aCorr="+aCorr+
        		        						" corr="+corr);
        		        			}
        		        			for (int iy=0;iy<tileSize;iy++){
        		        				for (int ix=0;ix<tileSize;ix++) corrTile [index++]+=subTileWindow[iy*tileSize+ix]*corr; //java.lang.ArrayIndexOutOfBoundsException: 64
        		        				index+=tileExtra;
        		        			}
        		        			// corrtile fades at the tile edges similar to window[]
        		        		}
        		        		// TODO: corrTile may all be zero - handle this?

//        		        		double corrThreshold=0.0;
        		        		for (int i=0;i<length;i++) linearMatchWeight[i]*= ((corrTile[i]>0.0)?corrTile[i]:0.0);
        		        		double sumWindow=0.0;
        		        		double sumLMW=0.0;
        		        		for (int i=0;i<length;i++) {
        		        			sumLMW+=linearMatchWeight[i];
        		        			sumWindow+=window[i];
        		        		}
        		        		// Normalize for Bayesian
        		        		double [] linearMatchWeightNorm=new double[length];
        		        		double scaleLMW=sumWindow/sumLMW;
        		        		for (int i=0;i<length;i++) linearMatchWeightNorm[i]=scaleLMW*linearMatchWeight[i]-window[i];

//        		        		int numNonZero=0;
        		        		// calculate pair[1][i] ~= a*pair[0][i] + b using weight linearMatchWeight[i] - for each color channel !
// debug show linearMatchWeight
            					if (debugLevel>3){ // +2 for selected tile
            						double [][] debugData={linearMatchWeight,linearMatchWeightNorm};
            						String [] debugTitles={"LMW","Norm"};
            						(new ShowDoubleFloatArrays()).showArrays(
            								debugData,
            								size,
            								size,
            								true,
            								"LMW"+plane+"_N"+nImg+"-"+sImg+"_sIndex-"+sIndex+"_tile"+tileX+"-"+tileY,
            								debugTitles);
            					}
            					double [][]staging=new double [slices[sIndex].length-1][];
            					double [][][]variance=new double[slices[sIndex].length-1][][];
    							double [] avgVariance=new double[slices[sIndex].length-1];
    							double [] avgRange=   new double[slices[sIndex].length-1];

            					for (int chn=0;chn<slices[sIndex].length-1;chn++){
            						if ((zMapVarMask& (1<<chn))!=0) {
                						variance[chn]=new double[2][length];
            							staging[chn]=this.photometric.initStaging(); // combined for all channels
            							double sumWV2=0.0;
            							double sumW=0.0;
            							double [] S0tot={0.0,0.0};
            							double [] S1tot={0.0,0.0};
            							double [] S2tot={0.0,0.0};

            							for (int i=0;i<length;i++) if (linearMatchWeightNorm[i]!=0.0) {
            								for (int n=0;n<2;n++) {

                								double  S1=0.0;
                								double  S2=0.0;
                								double  v2=0;
                								int numPix=0;
                								for (int d=0;d<dirs1.length;d++){
                									int i1=(i+dirs1[d]+length)%length;
                									if ((borderMask[i] & (borderMask[i1]<<1))==0) {
                										S1+=slices[sIndex][chn+1][n][i1];
                										S2+=(slices[sIndex][chn+1][n][i1]*slices[sIndex][chn+1][n][i1]);
                										numPix++;
                									}
                								}

                								v2=(S2-(S1*S1)/numPix)/numPix;
                								variance[chn][n][i]=Math.sqrt(v2);
//                								sumWV2+=linearMatchWeight[i]*v2;
//                								sumW+=  linearMatchWeight[i];
                								sumWV2+=window[i]*v2;
                								sumW+=  window[i];
                								S0tot[n]+=window[i];
                								S1tot[n]+=window[i]*slices[sIndex][chn+1][n][i];
                								S2tot[n]+=window[i]*slices[sIndex][chn+1][n][i]*slices[sIndex][chn+1][n][i];

            								}

            								this.photometric.addStagingSample(
            										staging[chn], // double []staging,
            										chn, //int chn,
            										linearMatchWeightNorm[i], //double weight, average==0
            										slices[sIndex][chn+1][0][i], //double v1,
            										slices[sIndex][chn+1][1][i], //double v2
            				        				// reduce weight depending on difference, scale to measured variance
            										varianceBlurScale, //double scaleVariance, (if 0 - do not use variance at all)
            				        				kLocal, // double kLocal, // 0 - use global varaiance, 1.0 - use local
            				        				nImg, // int nImg1, // first image number
            				        				sImg, //int nImg2, // second image number
            				        				variance[chn][0][i], //double var1, // first image variance
            				        				variance[chn][1][i],
            				        				debugLevel>5); //double var2){ // second image variance
            							}
            							double [] fVar=new double[2];
            							for (int n=0;n<2;n++){
            								fVar[n]=S2tot[n]/S0tot[n]-S1tot[n]*S1tot[n]/S0tot[n]/S0tot[n];
            							}

            							avgVariance[chn]=Math.sqrt(sumWV2/sumW);
            							avgRange[chn]=   Math.sqrt(0.5*(fVar[0]+fVar[1]));
            							if (debugLevel>3){
            								System.out.println("Average variance for plane "+plane+" nImg="+nImg+" sImg="+sImg+
            										" chn="+chn+" is "+avgVariance[chn]+" total tile variance="+avgRange[chn]);
            							}
            							// blur stage with scaled variance
            							this.photometric.blurStaging(
            									staging[chn],
            									chn,
            									avgVariance[chn]*varianceBlurScale);
            							//varianceBlurScale
            						}else{
            							staging[chn]=null;
            							variance[chn]=null;
            							avgVariance[chn]=0.0;
            							avgRange[chn]=0.0;
            						}
            					}
            					double sumRelVar=0.0;
            					for (int chn=0;chn<slices[sIndex].length-1;chn++) if (avgVariance[chn]>0.0){
                					if (debugLevel>3){
                						System.out.println("avgVariance["+chn+"]="+avgVariance[chn]+
                								" avgRange["+chn+"]="+avgRange[chn]+
                								" getAverageVariance("+nImg+","+chn+")="+this.photometric.getAverageVariance(nImg,chn));
                					}
//            						avgVariance[chn]/=this.photometric.getAverageVariance(nImg,chn);
//            						avgVariance[chn]*=avgVariance[chn];
                					avgRange[chn]/=this.photometric.getAverageVariance(nImg,chn);
                					avgRange[chn]*=avgRange[chn];
            						sumRelVar+=avgVariance[chn];
            					}
            					for (int chn=0;chn<slices[sIndex].length-1;chn++) if (avgVariance[chn]>0.0){
            						avgVariance[chn]/=sumRelVar;
                					if (debugLevel>3){
                						System.out.println(" normalized avgVariance["+chn+"]="+avgVariance[chn]);
                					}
            					}
            					if (debugLevel>3){
            						int numActive=0;
            						for (int i=0;i<staging.length;i++) if (staging[i]!=null) numActive++;
            						String [] channelNames={"Y","Cb","Cr","Aux"};
            						double [][] debugData=new double [numActive][];
            						String [] debugTitles=new String [numActive];
            						int index=0;
            						for (int i=0;i<staging.length;i++) if (staging[i]!=null) {
            							debugTitles[index]=channelNames[i];
            							debugData[index++]=staging[i];
            						}
            						(new ShowDoubleFloatArrays()).showArrays( //java.lang.ArrayIndexOutOfBoundsException: 0
            								debugData,
            								this.photometric.subdivAverage,
            								debugData[0].length/this.photometric.subdivAverage,
            								true,
            								"STAG"+plane+"_N"+nImg+"-"+sImg+"_sIndex-"+sIndex+"_tile"+tileX+"-"+tileY,
            								debugTitles);
            					}
            					planeStrength[plane][sIndex]=new double [length];
            					for (int i=0;i<length;i++) planeStrength[plane][sIndex][i]=0.0;
            					double [][] planeStrengthChn=null;
            					if (debugLevel>3) planeStrengthChn=new double [staging.length][length];
            					for (int chn=0;chn<slices[sIndex].length-1;chn++){
            						if (staging[chn]!=null) {
            							if (debugLevel>3) {
            								planeStrengthChn[chn]=new double[length];
                        					for (int i=0;i<length;i++) planeStrengthChn[chn][i]=0.0;
            							}
            							for (int i=0;i<length;i++) if (linearMatchWeight[i]>0){
            								double strength=this.photometric.getStrength(
            				        				staging[chn],
            				        				chn,
            				        				slices[sIndex][chn+1][0][i], //double v1,
            				        				slices[sIndex][chn+1][1][i]); //double v2
            								strength*=linearMatchWeight[i]; // remove??
            								if (window[i]>0.0) strength/=window[i];
            								planeStrength[plane][sIndex][i]+=avgVariance[chn]*normVarWeights[chn]*strength;
                							if (debugLevel>3) planeStrengthChn[chn][i]=strength;
            							}
            						} else {
            							if (debugLevel>3) planeStrengthChn[chn]=null;
            						}
            					}

            					if (debugLevel>3){
            						int numActive=0;
            						for (int i=0;i<staging.length;i++) if (staging[i]!=null) numActive++;
            						String [] channelNames={"Y","Cb","Cr","Aux"};
            						double [][] debugData=new double [numActive+1][];
            						String [] debugTitles=new String [numActive+1];
            						int index=0;
        							debugTitles[index]="combined";
        							debugData[index++]=planeStrength[plane][sIndex];
            						for (int i=0;i<channelNames.length;i++) if (planeStrengthChn[i]!=null) {
            							debugTitles[index]=channelNames[i];
            							debugData[index++]=planeStrengthChn[i];
            						}
            						(new ShowDoubleFloatArrays()).showArrays(
            								debugData,
            								size,
            								size,
            								true,
            								"L"+plane+"_N"+nImg+"-"+sImg+"_sIndex-"+sIndex+"_tile"+tileX+"-"+tileY,
            								debugTitles);
            					}
            					// now - need to save "likely" back to the tiles or compare right now?
            					// "unlikely - just correlation"
            					//

        					}

        				} // for (int sIndex=0;sIndex<sImgSet.length;sIndex++){


        			}

        		} // end for (int plane...)
// quick test - use the "most likely " plane
				double [][] planeStrengthCombo=new double [zTile.getNumberOfPlanes()][];
        		for (int plane=0;plane<zTile.getNumberOfPlanes();plane++){
        			if (zTile.getPlaneEnabled(plane) && (planeStrength[plane]!=null) && (planeStrength[plane][0]!=null)) {
        				planeStrengthCombo[plane]=planeStrength[plane][0].clone(); // later only center?
        				for (int sIndex=1;sIndex<sImgSet.length;sIndex++){
        					switch (matchStatMode) {
        					case 0: // multiply
        						for (int i=0;i<planeStrengthCombo[plane].length;i++){
        							if ((planeStrengthCombo[plane][i]<0.0) || (planeStrength[plane][sIndex][i]<0.0)) planeStrengthCombo[plane][i]=0.0;
        							else planeStrengthCombo[plane][i]*=planeStrength[plane][sIndex][i];
        						}
        					break;
        					case 1: // add
        						for (int i=0;i<planeStrengthCombo[plane].length;i++){ // do not add negative values
        							if (planeStrength[plane][sIndex][i]>0.0) planeStrengthCombo[plane][i]+=planeStrength[plane][sIndex][i];
        						}
            					break;
        					case 2: // min
        						for (int i=0;i<planeStrengthCombo[plane].length;i++){
        							if (planeStrength[plane][sIndex][i]< planeStrengthCombo[plane][i]) planeStrengthCombo[plane][i]=planeStrength[plane][sIndex][i];
        						}
            					break;
        					case 3: // max
        						for (int i=0;i<planeStrengthCombo[plane].length;i++){
        							if (planeStrength[plane][sIndex][i]> planeStrengthCombo[plane][i]) planeStrengthCombo[plane][i]=planeStrength[plane][sIndex][i];
        						}
            					break;
        					}
        				}
						double a=1.0/planeStrength[plane].length;
    					switch (matchStatMode) {
    					case 0: // multiply
    						for (int i=0;i<planeStrengthCombo[plane].length;i++){
    							if (planeStrengthCombo[plane][i]>0.0) planeStrengthCombo[plane][i]=Math.pow(planeStrengthCombo[plane][i],a);
    						}
    					break;
    					case 1: // add
    						for (int i=0;i<planeStrengthCombo[plane].length;i++){ // do not add negative values
    							planeStrengthCombo[plane][i]*=a;
    						}
        					break;
    					case 2: // min
    					case 3: // max
        					break;
    					}
        			} else {
        				planeStrengthCombo[plane]=null;
        			}
        		}

        		zTile.initAux(1);
				float [] aux0=new float[paddedLength];
//				float [] aux1=new float[paddedLength];

        		int [] planeIndex=new int[paddedLength];
        		for (int i=0;i<paddedLength;i++)planeIndex[i]=-1;
        		for (int plane=0;plane<zTile.getNumberOfPlanes();plane++) if (planeStrengthCombo[plane]!=null) {
        			int numUsed_dbg=0;
        			for (int i=0;i<paddedLength;i++){
        				int iPix= (margin+(i/paddedSize))*size +(margin+ (i%paddedSize));
        				if (	(planeStrengthCombo[plane][iPix]>0) && // should be positive
        						((planeIndex[i]<0) ||((planeStrengthCombo[planeIndex[i]]!=null) &&
        						(planeStrengthCombo[plane]!=null) &&

        						(planeStrengthCombo[planeIndex[i]][iPix] < planeStrengthCombo[plane][iPix])))) { //java.lang.NullPointerException
        					planeIndex[i]=plane;
        					numUsed_dbg++;
        				}
        			}
        			if (debugLevel>3) {
        				System.out.println ("planeLikelyTile("+tileX+","+tileY+", ...) plane="+plane+
        						" disparity="+zTile.getPlaneDisparity(plane)+" strength="+zTile.getPlaneStrength(plane)+
        						" used "+numUsed_dbg+"pixels"); // only scan inner pixels?
        			}
        		}
    			if (debugLevel>3) {
    				int numActivePlanes=0;
            		for (int plane=0;plane<zTile.getNumberOfPlanes();plane++) if (planeStrengthCombo[plane]!=null) numActivePlanes++;
            		if (numActivePlanes>0){
            			int index=0;
            			int [] planeNumbers=new int [numActivePlanes];
            			for (int plane=0;plane<zTile.getNumberOfPlanes();plane++) if (planeStrengthCombo[plane]!=null) planeNumbers[index++]=plane;
            			double [][] debugData = new double[(sImgSet.length+1)*numActivePlanes][];
            			String [] debugTitles=new String [(sImgSet.length+1)*numActivePlanes];
            			for (int i=0;i<numActivePlanes;i++){
            				int plane=planeNumbers[i];
            				String sDisparity=IJ.d2s(zTile.getPlaneDisparity(plane),2);
            				debugTitles[i]=nImg+":"+sDisparity; // oob 7
            				debugData[i]=planeStrengthCombo[plane];
            				for (int sIndex=0;sIndex<sImgSet.length;sIndex++){
                				debugTitles[i+(sIndex+1)*numActivePlanes]=nImg+"-"+sImgSet[sIndex]+":"+sDisparity; // oob 10
                				debugData[i+(sIndex+1)*numActivePlanes]=planeStrength[plane][sIndex];
            				}
            			}
						(new ShowDoubleFloatArrays()).showArrays(
								debugData,
								size,
								size,
								true,
								"RSLT_"+nImg+"_tile"+tileX+"-"+tileY,
								debugTitles);
            		}
    			}
        		for (int i=0;i<paddedLength;i++)aux0[i]=(float) ((planeIndex[i]>=0)?zTile.getPlaneDisparity( planeIndex[i]):0.0); // make it NaN?
        		zTile.setAux(0,aux0);
       		}

       		/**
       		 * calculate auto/inter correaltion of the smaller subtiles
       		 * @param doubleFHT - DHT instance reused by the thread (should be initialized)
       		 * @param data first/second array - now [2][32*32]
       		 * @param window - subtile window (defines subtile size) - now 8x8=64 - has to be initialized to an array with [0] equal to NaN
       		 * @param phaseCoeff phase correlation fraction (0 - normal, 1.0 - pure phase)
       		 * @param highPassSigma high-pass sigma - fraction of the full range
       		 * @param lowPassSigma low-pass sigma, in pixels
       		 * @param debugLevel
       		 * @return for each row and column - value of the correlation in the center (auto 00, auto 11, inter 01): [vert][hor]{corr00, corr11, corr01}
       		 */
           	public double [][][] subTileCorrelation(
        			DoubleFHT doubleFHT, // will generate if null, can be used to reuse initialization
        			double [][] data,   // difines size
        			double [] window, // defines tile size
        			double phaseCoeff, // phase correlation fraction (0 - normal, 1.0 - pure phase)
        			double highPassSigma,
        			double lowPassSigma,
        			int debugLevel){
        		int size=(int) Math.sqrt(data[0].length); // 32
        		int tileSize=(int) Math.sqrt(window.length); // 8
        		int tileStep=tileSize/2; //4
        		int tileExtra=size-tileSize; // 24??
        		int tileLength=tileSize*tileSize; // 64
 //       		int tileCenter=(tileSize+1)*tileSize/2; // 36
        		int numTilesRow=size/(2*tileStep); // 4
        		double [][][]  corrs=new double[numTilesRow][numTilesRow][3];
        		if (Double.isNaN(window[0])){
        			double [] window1d=doubleFHT.getHamming1d(tileSize); // Hamming
        			for (int i=0;i<window.length;i++) window[i]=window1d[i/tileSize]*window1d[i%tileSize];
        		}
        		int topLeft=(size+1)*(size/4-tileSize/4);
        		double []first= new double[tileLength];
        		double []second=new double[tileLength];
        		for (int tileY=0;tileY<numTilesRow;tileY++){
        			for (int tileX=0;tileX<numTilesRow;tileX++){
        				int indexIn=topLeft+ tileStep*(tileY*size+ tileX);
        				int indexOut=0;
        				for (int iy=0;iy<tileSize;iy++){
        					for (int ix=0;ix<tileSize;ix++){
        						first [indexOut]=data[0][indexIn];
        						second[indexOut++]=data[1][indexIn++];
        					}
        					indexIn+=tileExtra;
        				}
        				normalizeAndWindowGetDC (first,  window); // DC not needed here
        				normalizeAndWindowGetDC (second, window); // DC not needed here
        					double a2=0.0,b2=0.0,ab=0.0;
        					for (int i=0;i<first.length;i++){
        						a2+=first[i]* first[i];
        						b2+=second[i]*second[i];
        						ab+=first[i]* second[i];
        					}
        					corrs[tileY][tileX][0]=a2;
        					corrs[tileY][tileX][1]=b2;
        					corrs[tileY][tileX][2]=ab;
        			}
        		}
        		return corrs;
        	}





       		/**
       		 * Calculate occlusion by closer planes to a pixel resolution (use fractions for partially occluded pixels)
       		 * @param thisZMap [tileY]tileX]
       		 * @param xc center of the WOI, X
       		 * @param yc center of the WOI, Y
       		 * @param disparity disparity only objects with higher disparity may block (Add margin to prevent blocking by the same plane?)
       		 * @param foregroundDisparityMargin only consider disparity+foregroundDisparityMargin as occluding
       		 * @param disparityTolerance - consider disparity+/-disparityTolerance as "this plane" to remove disabled pixels (in NaN - ignore)
       		 * @param debugLevel debug level
       		 * @return paddedLength (default 20x20=400) array whity 0.0 - completely blocked, 1.0 - completely not blocked
       		 */
       		public double [] getNonOccluded(
       				ZTile [][] thisZMap,
       				double xc,
       				double yc,
       				double disparity,
       				double foregroundDisparityMargin,
       				double disparityTolerance,
       				int debugLevel
       		){

       			int [] iXYtl={(int) Math.floor(xc)-this.overlapStep/2,(int) Math.floor(yc)-this.overlapStep/2};
       			int [] iXYbr={(int) Math.ceil(xc)+this.overlapStep/2-1,(int) Math.ceil(yc)+this.overlapStep/2-1};
       			double []dXY={xc-Math.floor(xc),yc-Math.floor(yc)};
       			boolean [][][] enMasks ={{null,null},{null,null}}; // cache for tile masks
/*       			boolean [] bothXY={
       					iXYbr[0]/this.overlapStep>iXYtl[0]/this.overlapStep,
       					iXYbr[1]/this.overlapStep>iXYtl[1]/this.overlapStep};
       			int tileX=iXYtl[0]/this.overlapStep;
       			int tileY=iXYtl[1]/this.overlapStep;
       			boolean [] bothXY={
       					((int) Math.floor(((double) iXYbr[0])/this.overlapStep))>((int) Math.floor(((double)iXYtl[0])/this.overlapStep)),
       					((int) Math.floor(((double) iXYbr[1])/this.overlapStep))>((int) Math.floor(((double)iXYtl[1])/this.overlapStep))};
       					*/
       			boolean [] bothXY={
       					Math.floor(((double) iXYbr[0])/this.overlapStep)>Math.floor(((double)iXYtl[0])/this.overlapStep),
       					Math.floor(((double) iXYbr[1])/this.overlapStep)>Math.floor(((double)iXYtl[1])/this.overlapStep)};
       			int tileX=(int)Math.floor(((double) iXYtl[0])/this.overlapStep);
       			int tileY=(int)Math.floor(((double) iXYtl[1])/this.overlapStep);
       			if (debugLevel>3) { //(tileX<0) || (tileY<0)){
       				System.out.println("getNonOccluded(thisZMap,"+xc+","+yc+","+disparity+
       						","+foregroundDisparityMargin+","+disparityTolerance+","+debugLevel+
       						") :iXYtl={"+iXYtl[0]+","+iXYtl[1]+"} :iXYbr={"+iXYbr[0]+","+iXYbr[1]+"} :dXY={"+dXY[0]+","+dXY[1]+"} tileX="+
       						tileX+" tileY="+tileY);
       			}
       			double disparityThresholdFg=disparity+foregroundDisparityMargin;
//       			enMasks[0][0]=thisZMap[tileY][tileX].getNonOccluded(disparityThresholdFg);
//       			enMasks[0][1]=(bothXY[0])?thisZMap[tileY][tileX+1].getNonOccluded(disparityThresholdFg):null;
//       			enMasks[1][0]=(bothXY[1])?thisZMap[tileY+1][tileX].getNonOccluded(disparityThresholdFg):null;
//       			enMasks[1][0]=(bothXY[0]&&bothXY[1])?thisZMap[tileY+1][tileX+1].getNonOccluded(disparityThresholdFg):null;
       			enMasks[0][0]=getEnabledNonOccluded(
       					thisZMap,
       					tileX,
       					tileY,
       					disparityThresholdFg,
       					disparity,
           				disparityTolerance,
           				debugLevel);
       			enMasks[0][1]=(bothXY[0])?
       					getEnabledNonOccluded(
       	       					thisZMap,
       	       					tileX+1,
       	       					tileY,
       	       					disparityThresholdFg,
       	       					disparity,
       	           				disparityTolerance,
       	           				debugLevel):null;
       			enMasks[1][0]=(bothXY[1])?
       					getEnabledNonOccluded(
       					thisZMap,
       					tileX,
       					tileY+1,
       					disparityThresholdFg,
       					disparity,
           				disparityTolerance,
           				debugLevel):null;
       			enMasks[1][1]=(bothXY[0]&&bothXY[1])?getEnabledNonOccluded(
       					thisZMap,
       					tileX+1,
       					tileY+1,
       					disparityThresholdFg,
       					disparity,
           				disparityTolerance,
           				debugLevel):null;
//    			int paddedSize=thisZMap[tileY][tileX].getPaddedSize();
//    			int margin=this.overlapStep/2-thisZMap[tileY][tileX].getOverlap();
//    			int padding=thisZMap[tileY][tileX].getOverlap();
    			int padding=getPadding();
    			int [][] whereX=new int [this.paddedSize+(bothXY[0]?1:0)][2];
    			int [][] whereY=new int [this.paddedSize+(bothXY[1]?1:0)][2];
    			int [] tile00tl= {tileX*this.overlapStep-padding,tileY*this.overlapStep-padding};
    			for (int i=0;i<whereX.length;i++){
					int px=iXYtl[0]-2+i; // absolute pixel X
					whereX[i][0]=px-tile00tl[0];
					whereX[i][1]=0;
					if (bothXY[0] && (px>=(tile00tl[0]+this.overlapStep))){
						whereX[i][0]-=this.overlapStep;
						whereX[i][1]=1;
					}
    			}
    			for (int i=0;i<whereY.length;i++){
					int py=iXYtl[1]-2+i; // absolute pixel X
					whereY[i][0]=py-tile00tl[1];
					whereY[i][1]=0;
					if (bothXY[1] && (py>=(tile00tl[1]+this.overlapStep))){
						whereY[i][0]-=this.overlapStep;
						whereY[i][1]=1;
					}
    			}
    			int longPaddedSize=this.paddedSize+1;
    			boolean [] extraBits=new boolean[longPaddedSize*longPaddedSize];
    			for (int y=0;y<whereY.length;y++) {
    				int index=y*longPaddedSize;
    				for (int x=0;x<whereX.length;x++){
    					if ((whereY[y][0]<0) || (whereX[x][0]<0) || (whereY[y][1]<0) || (whereX[x][1]<0) ||
    							(enMasks[whereY[y][1]][whereX[x][1]]==null)||
    							((whereY[y][0]*this.paddedSize + whereX[x][0]))>=enMasks[whereY[y][1]][whereX[x][1]].length){
    						System.out.println(
    								" xc="+xc+
    								" yc="+yc+
    								" tileX="+tileX+
    								" tileY="+tileY+
    								" iXYtl={"+iXYtl[0]+","+iXYtl[1]+"}"+
    								" iXYbr={"+iXYbr[0]+","+iXYbr[1]+"}"+
    								" bothXY={"+bothXY[0]+","+bothXY[1]+"}"+
    								" whereX.length="+whereX.length+
    								" whereY.length="+whereY.length+
    								" whereY["+y+"][0]="+whereY[y][0]+
    								" whereY["+y+"][1]="+whereY[y][1]+
    								" whereX["+x+"][0]="+whereX[x][0]+
    								" whereX["+x+"][1]="+whereX[x][1]+
    								" enMasks["+whereY[y][1]+"]["+whereX[x][1]+"].length="+enMasks[whereY[y][1]][whereX[x][1]].length+
    								" whereY[y][0]*this.paddedSize + whereX[x][0]="+(whereY[y][0]*this.paddedSize + whereX[x][0])+
    								" Math.floor(((double)iXYtl[0])/this.overlapStep)="+Math.floor(((double)iXYtl[0])/this.overlapStep)+
    								" Math.floor(((double)iXYbr[0])/this.overlapStep)="+Math.floor(((double)iXYbr[0])/this.overlapStep)+
    								" Math.floor(((double)iXYtl[1])/this.overlapStep)="+Math.floor(((double)iXYtl[1])/this.overlapStep)+
    								" Math.floor(((double)iXYbr[1])/this.overlapStep)="+Math.floor(((double)iXYbr[1])/this.overlapStep)+
    								" (int) Math.floor(((double)iXYtl[0])/this.overlapStep)="+((int) Math.floor(((double)iXYtl[0])/this.overlapStep))+
    								" (int) Math.floor(((double)iXYbr[0])/this.overlapStep)="+((int) Math.floor(((double)iXYbr[0])/this.overlapStep))+
    								" (int) Math.floor(((double)iXYtl[1])/this.overlapStep)="+((int) Math.floor(((double)iXYtl[1])/this.overlapStep))+
    								" (int) Math.floor(((double)iXYbr[1])/this.overlapStep)="+((int) Math.floor(((double)iXYbr[1])/this.overlapStep))/*+
    								" ((int) (Math.floor((double)iXYbr[1])/this.overlapStep))>((int) (Math.floor((double)iXYtl[1])/this.overlapStep))="+
    								(((int) (Math.floor((double)iXYbr[1])/this.overlapStep))>((int) (Math.floor((double)iXYtl[1])/this.overlapStep)))+
    								" ((Math.floor((double)iXYbr[1])/this.overlapStep))>((Math.floor((double)iXYtl[1])/this.overlapStep))="+
    								(((Math.floor((double)iXYbr[1])/this.overlapStep))>((Math.floor((double)iXYtl[1])/this.overlapStep)))+
    								" "+((int) (Math.floor((double)iXYbr[1])/this.overlapStep))+">"+((int) (Math.floor((double)iXYtl[1])/this.overlapStep))+"="+
    								(((int) (Math.floor((double)iXYbr[1])/this.overlapStep))>((int) (Math.floor((double)iXYtl[1])/this.overlapStep)))*/
    								);
    					}
/*
        					((int) (Math.floor((double)iXYbr[0])/this.overlapStep))>((int) (Math.floor((double)iXYtl[0])/this.overlapStep)),
       					((int) (Math.floor((double)iXYbr[1])/this.overlapStep))>((int) (Math.floor((double)iXYtl[1])/this.overlapStep))};

 */
    					extraBits[index++]=enMasks[whereY[y][1]][whereX[x][1]][whereY[y][0]*this.paddedSize + whereX[x][0]];
    				}
    			}
    			if (debugLevel>5){
					(new ShowDoubleFloatArrays()).showArrays(
							extraBits,
							longPaddedSize,
							longPaddedSize,
							"EB_-X"+tileX+"-Y"+tileY
					);


    			}
    			int index=0;
    			double [] result=new double[this.paddedSize*this.paddedSize];
    			double [][]k={
    					{(1-dXY[0])*(1-dXY[1]), (  dXY[0])*(1-dXY[1])},
    					{(1-dXY[0])*(  dXY[1]), (  dXY[0])*(  dXY[1])}};
    			for (int y=0;y<this.paddedSize;y++) {
    				int indexS=y*longPaddedSize;
    				for (int x=0;x<this.paddedSize;x++){
    					result[index]= extraBits[indexS]?k[0][0]:0.0;
    					if (bothXY[0] && extraBits[indexS+1]) result[index]+= k[0][1];
    					if (bothXY[1]){
    						if (extraBits[indexS+this.paddedSize]) result[index]+= k[1][0];
    						if (bothXY[0] && extraBits[indexS+this.paddedSize+1]) result[index]+= k[1][1];
    					}
    					index++;
    					indexS++;
    				}
    			}
    			return result;
       		}
    		boolean [] getEnabledNonOccluded(
    				ZTile [][] thisZMap,
    				int tileX,
    				int tileY,
    				double disparityThresholdFg,
       				double disparity,
       				double disparityTolerance,
       				int debugLevel){
    			if ((tileX<0) || (tileY<0) || (tileX>=this.tilesX) || (tileY>=this.tilesY) || (thisZMap[tileY][tileX]==null)){
    				boolean [] empty=new boolean[this.paddedSize*this.paddedSize];
    				for (int i=0;i<empty.length;i++) empty[i]=false;
    				return empty;
    			}
    			if (debugLevel>5){
    				boolean [] geno=thisZMap[tileY][tileX].getEnabledNonOccluded(
           					disparityThresholdFg,
           					disparity,
               				disparityTolerance,
               				debugLevel);
    				int dbgSize=(int)Math.sqrt(geno.length);
					(new ShowDoubleFloatArrays()).showArrays(
							geno,
							dbgSize,
							dbgSize,
							"GENO_-X"+tileX+"-Y"+tileY
					);


    			}
				return thisZMap[tileY][tileX].getEnabledNonOccluded(
       					disparityThresholdFg,
       					disparity,
           				disparityTolerance,
           				debugLevel);
    		}
        	/**
        	 * Calculate fine-step (4 pix?) correlation between the pair of square images (32x32) that are already shifted to compensate for known disparity
        	 * only center area (with 1/4 margins) will be calculated, if the correlation is strong enough, the updated disparity arrays will be calculated fro these tiles
        	 *
        	 * @param data a pair of square arrays to be correlated (32x32 pixels)
        	 * @param window cosine mask (8x8=64 pixels long). If the center element is 0.0 - will generate array
        	 * @param dXY unity vector in the direction of the disparity between the two images
        	 * @param phaseCoeff 0.0 - normal correlation, 1.0 - pure phase correlation
        	 * @param debugLevel debug level
        	 * @return smooth array, same dimension as data[0], with only center area trusted, proportional to correlation at zero
        	 */
        	public double [] localCorrelation(
        			DoubleFHT doubleFHT, // will generate if null, can be used to reuse initialization
        			double [][] data,   // difines size
        			double [] window, // defines tile size
        			double [] dXY,    // unity vector defines disparity direction
        			double phaseCoeff, // phase correlation fraction (0 - normal, 1.0 - pure phase)
        			double highPassSigma,
        			double lowPassSigma,
        			int     normalize,
        			int debugLevel
        			){
        		return localCorrelation(
        				doubleFHT, // will generate if null, can be used to reuse initialization
    			data,   // difines size
    			window, // defines tile size
    			dXY,    // unity vector defines disparity direction
    			phaseCoeff, // phase correlation fraction (0 - normal, 1.0 - pure phase)
    			highPassSigma,
    			lowPassSigma,
    			normalize,
    			debugLevel,
    			null);
        	}

        	public double [] localCorrelation(
        			DoubleFHT doubleFHT, // will generate if null, can be used to reuse initialization
        			double [][] data,   // difines size
        			double [] window, // defines tile size
        			double [] dXY,    // unity vector defines disparity direction
        			double phaseCoeff, // phase correlation fraction (0 - normal, 1.0 - pure phase)
        			double highPassSigma,
        			double lowPassSigma,
        			int     normalize,
        			int debugLevel,
        			String title
        			){
        		int size=(int) Math.sqrt(data[0].length); // 32
        		int tileSize=(int) Math.sqrt(window.length); // 8
        		int tileStep=tileSize/2; //4
        		int tileExtra=size-tileSize; // 24??
        		int tileLength=tileSize*tileSize; // 64
        		int tileCenter=(tileSize+1)*tileSize/2; // 36

        		int numTilesRow=size/(2*tileStep); // 4
        		if (doubleFHT==null) doubleFHT= new DoubleFHT();
        		if (Double.isNaN(window[0])){
        			double [] window1d=doubleFHT.getHamming1d(tileSize); // Hamming
        			for (int i=0;i<window.length;i++) window[i]=window1d[i/tileSize]*window1d[i%tileSize];
        		}
        		int topLeft=(size+1)*(size/4-tileSize/4);
        		double []first= new double[tileLength];
        		double []second=new double[tileLength];
        		double [] result=null;
        			result=new double [size*size];
        			for (int i=0;i<result.length;i++) result[i]=0.0;
        		double []debugArray=null;
        		if (title!=null) {
        			debugArray = new double [numTilesRow*tileSize*numTilesRow*tileSize];
        		}

        		for (int tileY=0;tileY<numTilesRow;tileY++){
        			for (int tileX=0;tileX<numTilesRow;tileX++){
        				int indexIn=topLeft+ tileStep*(tileY*size+ tileX);
        				int indexOut=0;
        				for (int iy=0;iy<tileSize;iy++){
        					for (int ix=0;ix<tileSize;ix++){
        						first [indexOut]=data[0][indexIn];
        						second[indexOut++]=data[1][indexIn++];
        					}
        					indexIn+=tileExtra;
        				}
        				normalizeAndWindowGetDC (first,  window); // DC not needed here
        				normalizeAndWindowGetDC (second, window); // DC not needed here
        				double sc=0.0;
        				if (normalize==2) {
        					double a2=0.0,b2=0.0,ab=0.0;
        					for (int i=0;i<first.length;i++){
        						a2+=first[i]* first[i];
        						b2+=second[i]*second[i];
        						ab+=first[i]* second[i];
        					}
        					sc=ab/Math.sqrt(a2*b2);

        				} else {
        					second=doubleFHT.correlate (second,first,highPassSigma,lowPassSigma,phaseCoeff); // second will be modified
        					if (normalize==1) {
        						for (int i=0;i<tileLength;i++)sc+=second[i]*second[i];
        						sc=Math.sqrt(sc/tileLength);
        						if ((debugLevel>3) && (tileY==0) && (tileX==0))	for (int i=0;i<second.length;i++) System.out.println("second["+i+"]="+second[i]);
        						if (debugLevel>2) System.out.println("localCorrelation(): phaseCoeff="+phaseCoeff+" sc="+sc+" second["+tileCenter+"]="+second[tileCenter]);

        						sc=second[tileCenter]/sc; // relative correlation strength at zero.
        					} else {
        						sc=second[tileCenter];
        					}
        				}
        				if (debugLevel>2) System.out.println("localCorrelation(): sc="+sc);

        				// accumulate result
        	    			int index=topLeft+ tileStep*(tileY*size+ tileX);
            				for (int iy=0;iy<tileSize;iy++){
            					for (int ix=0;ix<tileSize;ix++) result [index++]+=sc*window[iy*tileSize+ix];
            					index+=tileExtra;
            				}
            			if (debugArray!=null) {
            				int dbgIndex=0;
            				for (int i=0;i<tileSize;i++){
            					for (int j=0;j<tileSize;j++){
            						debugArray[(tileY*tileSize+i)* (numTilesRow*tileSize)+ (tileX*tileSize+j)]=(normalize==2)?sc:second[dbgIndex++];
            					}
            				}
            			}
        			}
        		}
    			if (debugArray!=null) {
    				(new ShowDoubleFloatArrays()).showArrays(
        					debugArray,
        					numTilesRow*tileSize,
        					numTilesRow*tileSize,
        					title);
    			}

/*
        		if (debugArray!=null){
        			String [] dbgTitles={"corr","first","second"};
        			(new showDoubleFloatArrays()).showArrays(
        					debugArray,
        					size*subPixel,
        					size*subPixel,
        					true,
        					"correlation-refine",
        					dbgTitles);
        		}
*/
        		return result;
        	}

        	/**
        	 * Create a stack of channel images (YCbCr and optionally with external data
        	 * @param xc1 first image center X
        	 * @param yc1 first image center Y
        	 * @param size square side, power of 2 (will use twice larger internally)
        	 * @param first image number of the first in pair
        	 * @param second image number of the second in pair
        	 * @param dxA - how much to shift the second image from the first, horizontally, right
        	 * @param dyA - how much to shift the second image from the first, vertically,down
        	 * @param channelMask - channels to use +1 - Y, +2 - Cb, +4-Cr +8 - Aux.
        	 * @param imageData [image][alpha,Y,Cb,Cr,ext][] array of source image data
        	 * @param doubleSizeOutput - output twice the required size tile (to be able to probe around needed pixels
        	 * @param window - 1d in linescan order, (2*size)*(2*size). Probably - flat==1.0 in the inner size*size area. May be null - will be generated internally
        	 * @param doubleFHT DoubleFHT or null (to save resources)
        	 * @param debugLevel debug level
        	 * @return 3d array of the shifted slices from 2 images [channel][image number 0,1][pixels]
        	 * first tile (and the result will be shifted half fractional distance to the sescond, second - all intyeger part and half fractional towards the first
        	 */
        	public double [][][] getShiftedSlices(
        			int xc1,
        			int yc1,
        			int size,
        			int first,
        			int second,
        			double dxA,
        			double dyA,
        			int channelMask,
        			double [][][] imageData,
        			int imageWidth,
        			boolean doubleSizeOutput, // make it always true?
        			double [] window, // should be 4*size*size long;
        			DoubleFHT doubleFHT, // to reuse tables
        			int debugLevel
        			){
        		int numImgChannels=5; // alpha-Y-Cb-Cr-Ext
        		if (doubleFHT==null) doubleFHT=new DoubleFHT();
        		int iDx= (int) Math.round(dxA); //
        		int iDy= (int) Math.round(dyA);
        		int [] xc={xc1,xc1+ iDx};
        		int [] yc={yc1,yc1+ iDy};
        		double [] dx={-(dxA-iDx)/2,(dxA-iDx)/2}; // sub-pixel shift of the first image, second image
        		double [] dy={-(dyA-iDy)/2,(dyA-iDy)/2};

        		boolean [] shiftImage= {(dx[0]!=0.0) || (dy[0]!=0.0), (dx[1]!=0.0) || (dy[1]!=0.0)};

        		double [][][] result=new double[numImgChannels][][];
        		result[0]=null; // alpha
        		for (int i=0;i<numImgChannels-1;i++) {
        			if ((channelMask & (1<<i))!=0) result[i+1]=new double [2][];
        			else result[i+1]=null;
        		}
        		if (debugLevel>4){
        			double d=Math.sqrt(dxA*dxA+dyA*dyA);
        			System.out.println("getShiftedSlices(): xc1="+xc1+" yc1="+yc1+" size="+size+" first="+first+" second="+second);
        			System.out.println("getShiftedSlices(): d="+IJ.d2s(d,2)+" dX="+IJ.d2s(dxA,2)+" dY="+IJ.d2s(dyA,2));
        			System.out.println("getShiftedSlices(): shiftImage[0]="+shiftImage[0]+" shiftImage[1]="+shiftImage[1]);
        			System.out.println("getShiftedSlices() first image:   dX="+IJ.d2s(dx[0],2)+" dY="+IJ.d2s(dy[0],2));
        			System.out.println("getShiftedSlices() second image:  dX="+IJ.d2s(dx[1],2)+" dY="+IJ.d2s(dy[1],2));
        			System.out.println("getShiftedSlices() second image: iDx="+iDx+" iDyY="+iDy);
        			System.out.println("getShiftedSlices() second image:  xc="+(xc1+iDx)+" yc="+(yc1+iDy));
        		}
        		int doubleSize=2*size;
        		int doubleLength=doubleSize*doubleSize;
        		double [] slice= new double[doubleLength];
        		double [] imgSlice;
        		int [] img= {first,second};
        		/*
        		if (window==null) {
        			window = new double [doubleLength];
        			window[0]=Double.NaN;
        		}
        		*/
        		if (Double.isNaN(window[0])) {
        			double[] window1d=doubleFHT.getHamming1d(size); // Hamming
        			int index=0;
            		int halfSize=size/2;
            		int size32=3*size/2;
        			for (int iy=0;iy<doubleSize;iy++) {
        				double wy=(iy<halfSize)?window1d[iy]:((iy>size32)?window1d[iy-size]:1.0);
        				for (int ix=0;ix<doubleSize;ix++) {
            				double wx=(ix<halfSize)?window1d[ix]:((ix>size32)?window1d[ix-size]:1.0);
        					window[index++]=wx*wy;
        				}
        			}
        		}
        		if (debugLevel>4){
        			System.out.print("getShiftedSlices(): channelMask="+channelMask);

        	   		if (debugLevel>3){
        					(new ShowDoubleFloatArrays()).showArrays(
        								window,
        								doubleSize,
        								doubleSize,
        								"window");
            		}
        		}
        		for (int iImg=0;iImg<img.length;iImg++){
        			int nImg=img[iImg];
        			for (int channel=0; channel<result.length; channel++) if (result[channel]!=null){
        				imgSlice=imageData[nImg][channel];  //java.lang.ArrayIndexOutOfBoundsException: 3
        				if (debugLevel>2) System.out.println("xc["+iImg+"]="+xc[iImg]+" yc["+iImg+"]="+yc[iImg]+" shiftImage["+iImg+"]="+shiftImage[iImg]);

        				slice= getSelection(
        						imgSlice, // source image/channel slice
        						slice,
        						doubleSize, //int width,
        						doubleSize, //int height,
        						imageWidth,
        						xc[iImg] , //xc,
        						yc[iImg]); //yc);
        				if (debugLevel>4) {
        					(new ShowDoubleFloatArrays()).showArrays(
    								slice,
    								doubleSize,
    								doubleSize,
    								"slice-getShiftedSlices");
        				}

        				// normalize and save DC
        				double dc=0.0;
        				if (shiftImage[iImg]){
        					dc=normalizeAndWindowGetDC (slice, window); //windowInterpolation
            				if (debugLevel>4) {
            					(new ShowDoubleFloatArrays()).showArrays(
        								slice,
        								doubleSize,
        								doubleSize,
        								"slice-normalized-dc"+dc);
            				}

//        					doubleFHT.shift(slice, dx[iImg], dy[iImg]);
        					doubleFHT.shift(slice, -dx[iImg], -dy[iImg]);
            				if (debugLevel>4) {
            					(new ShowDoubleFloatArrays()).showArrays(
        								slice,
        								doubleSize,
        								doubleSize,
        								"slice-shifted");
            				}
        				} else {
        					if (debugLevel>2) System.out.println("No shift is needed");
        				}
//    					int oSliceNumber=img.length*cN+iImg;
    					if (debugLevel>2) System.out.println("getShiftedSlices() dc="+dc);
        				if (doubleSizeOutput) {
        					if (debugLevel>2) System.out.println("doubleLength="+doubleLength+"doubleSize="+doubleSize);
        					result[channel][iImg]=new double [doubleLength];
        					int oIndex=0;
        					for (int iY=0;iY<doubleSize;iY++) {
        						for (int iX=0;iX<doubleSize;iX++){
        							if ((oIndex>=result[channel][iImg].length) || (oIndex>=slice.length) || (oIndex>=window.length)){
        								System.out.println("iX="+iX+" iY="+iY+" oIndex="+oIndex);

        							}
        							result[channel][iImg][oIndex]=shiftImage[iImg]?(slice[oIndex]/window[oIndex]+dc):slice[oIndex];// no NaN with Hamming
        							oIndex++;
        						}
        					}
        				}else{
        					result[channel][iImg]=new double [size*size];
        					int iIndex=size*size+size/2; // top left index of the inner square size*size
        					int oIndex=0;
        					for (int iY=0;iY<size;iY++) {
        						for (int iX=0;iX<size;iX++){
        							result[channel][iImg][oIndex++]=shiftImage[iImg]?(slice[iIndex]/window[iIndex]+dc):slice[iIndex];
        							iIndex++;
        						}
        						iIndex+=size;
        					}
        				}
        			}
        		}
        		return result;
        	}

           	public double [] getSelection(
        			double [] imageSlice, // one image/channel slice
        			double [] selection, // null or array to reuse
    				int width,
    				int height,
        			int fullWidth,
    				int xc,
    				int yc){
        		int length=width*height;
        		int fullHeight=imageSlice.length/fullWidth;
        		if (selection ==null) selection = new double[length];
    				int y0=yc-height/2;
    				int x0=xc-width/2;
    	    		for (int iy=0;iy<height;iy++) {
    	    			int srcY=iy+y0;
    	    			boolean oob=(srcY<0) || (srcY>=fullHeight);
    	    			for (int ix=0;ix<width;ix++){
            				int oIndex=iy*width+ix;
        					int srcX=x0+ix;
        					if (oob ||(srcX<0) || (srcX>=fullWidth)) {
        						if (oIndex>=selection.length) System.out.println("\ngetSelection(imageSlice["+imageSlice.length+"],selection["+selection.length+"]"+
        								","+fullWidth+","+width+","+height+","+xc+","+yc+") ix="+ix+" iy="+iy+" oIndex="+oIndex);
        						selection[oIndex]=0.0;
        					} else {
            						selection[oIndex]=imageSlice[srcY*fullWidth+srcX];
        					}
    	    			}
    	    		}
        		return selection;
        	}




        	public void setupZMap(
        			Rectangle woi, // in tiles - may be
        			final int nImg,
        			final int maxNumber,
        			final double minFirst,
        			final double minAbsolute,
        			final double minRelative,
        			final double mergeMax,
        			final int overlap,
        			final double zMapMinForeground,
        			final int threadsMax,
        			final boolean showProgress,
        			final int debugLevel){
        		if (this.centerPixels==null) {
    				String msg="Centered disparity data is not defined";
    				IJ.showMessage("Error",msg);
    				throw new IllegalArgumentException (msg);
        		}
        		if (this.syntheticPixels==null) moveDisparityFromCenter(threadsMax,showProgress,debugLevel);
        		final float [] thisSyntheticPixels= this.syntheticPixels[nImg];
        		if (woi==null) woi=new Rectangle (0, 0, this.tilesX,this.tilesY);
        		this.zMapWOI=new Rectangle();
        		this.zMapWOI.x=(woi.x>=0)?woi.x:0;
        		this.zMapWOI.y=(woi.y>=0)?woi.y:0;
        		this.zMapWOI.width= (((woi.x+woi.width) <=this.tilesX)?(woi.x+woi.width): this.tilesX)-this.zMapWOI.x;
        		this.zMapWOI.height=(((woi.y+woi.height)<=this.tilesY)?(woi.y+woi.height):this.tilesY)-this.zMapWOI.y;
        		final Rectangle zMapWOI=this.zMapWOI;
        		if ((this.zMap==null) || (this.zMap[nImg]==null)) initZMap(nImg);
        		final ZTile [][] thisZMap=this.zMap[nImg];
        		final int tilesX=this.zMapWOI.width;
        		final int tiles=this.zMapWOI.width*this.zMapWOI.height;
        		if (debugLevel>2) System.out.println("setupZMap() woi.x="+woi.x+" woi.y="+woi.y+
        				" woi.width="+woi.width+" woi.height="+woi.height);
        		if (debugLevel>2) System.out.println("setupZMap() zMapWOI.x="+zMapWOI.x+" zMapWOI.y="+zMapWOI.y+
        				" zMapWOI.width="+zMapWOI.width+" zMapWOI.height="+zMapWOI.height);
        		this.paddedSize=this.overlapStep+2*overlap;
        		this.innerMask=new BitSet(this.paddedSize*this.paddedSize);
        		for (int i=0;i<this.overlapStep;i++) for (int j=0;j<this.overlapStep;j++){
        			this.innerMask.set((overlap+i)*this.paddedSize+overlap+j);
        		}
//   		private int  []   borderMask=null; // will be provided to zMap +1:top+2:bottom+8:left+16:right (to prevent roll over when iterating around
        		int tileSize=2*this.overlapStep;
        		this.borderMask=new int [tileSize*tileSize];
        		for (int i=0;i<tileSize;i++) for (int j=0;j<tileSize;j++){
        			this.borderMask[i*tileSize+j]=((i==0)?5:0)+((i==(tileSize-1))?2:0)+((j==0)?0x50:0)+((j==(tileSize-1))?0x20:0);
        		}


        		final Thread[] threads = newThreadArray(threadsMax);
		   		final AtomicInteger tileIndexAtomic     = new AtomicInteger(0);
		   		if (showProgress) IJ.showProgress(0.0);
		   		if (showProgress) IJ.showStatus("Setting up zMap for image "+nImg+" ("+(nImg+1)+" of "+this.disparityScales.length+") ...");
		   		for (int ithread = 0; ithread < threads.length; ithread++) {
		   			threads[ithread] = new Thread() {
		   				@Override
						public void run() {
		   					for (int tile=tileIndexAtomic.getAndIncrement(); tile<tiles;tile=tileIndexAtomic.getAndIncrement()){
		   						setupZMapTile(
		   								zMapWOI.x+tile%tilesX, //int tileX,
		   								zMapWOI.y+tile/tilesX, //int tileY,
		   								thisSyntheticPixels,
		   								thisZMap,
		   			        			maxNumber,
		   			        			minFirst,
		   			        			minAbsolute,
		   			        			minRelative,
		   			        			mergeMax,
		   			        			overlap,
		   			        			zMapMinForeground,
		   			        			debugLevel);
		   						if (showProgress){
		   							final int finalTile=tile;
		   							SwingUtilities.invokeLater(new Runnable() {
		   								@Override
										public void run() {
		   									IJ.showProgress(finalTile,tiles);
		   								}
		   							});
		   						}
		   					}
		   				}
		   			};
		   		}
		   		startAndJoin(threads);
		   		IJ.showProgress(1.0);
        	}

        	public float [][] renderZMap(
        			Rectangle woi,
        			boolean overlap){
        		return renderZMap(woi,-1, overlap);
        	}
        	public float[] renderZMap(
        			int nImg,
        			Rectangle woi,
        			boolean overlap){
        		return overlap?renderZMapOverlap(nImg,woi,-1):renderZMap(nImg,woi,-1);
        	}
        	public float [][] renderZMap(
        			Rectangle woi,
        			int auxLayer,
        			boolean overlap
        			){
        		float [][] result = new float [this.zMap.length][];
        		for (int nImg=0;nImg<result.length;nImg++) result[nImg]=overlap?renderZMapOverlap(nImg,woi,auxLayer):renderZMap(nImg,woi,auxLayer);
        		return result;
        	}
        	public float [] renderZMap(
        			int nImg,
        			Rectangle woi,
        			int auxLayer){
        		float [] result = new float [woi.width*woi.height];
        		int yBottom=woi.y+woi.height;
        		int xRight= woi.x+woi.width;
        		int tileYmin=woi.y/this.overlapStep;
        		int tileYmax=(yBottom-1)/this.overlapStep+1; // +1 to the last needed
        		if (tileYmax>this.tilesY) tileYmax=this.tilesY;
        		int tileXmin=woi.x/this.overlapStep;
        		int tileXmax=(xRight-1)/this.overlapStep+1; // +1 to the last needed
        		if (tileXmax>this.tilesX) tileXmax=this.tilesX;

        		for (int i=0;i<result.length;i++) result[i]=Float.NaN;
        		for (int tileY=tileYmin;tileY<tileYmax;tileY++) for (int tileX=tileXmin;tileX<tileXmax; tileX++){
        			ZTile zTile=this.zMap[nImg][tileY][tileX];
        			if (zTile!=null){
        				float [] tile=(auxLayer<0)?zTile.getZmap():zTile.getAux(auxLayer);
        				if (tile!=null) {
        					int tileSize=this.overlapStep+2*zTile.overlap;
        					for (int i=0;i<this.overlapStep;i++) {
        						int oY=(tileY*this.overlapStep-woi.y)+i;
        						if ((oY>=0) && (oY<woi.height)){
        							int iIndex=tileSize*(zTile.overlap+i)+zTile.overlap;
        							int oIndex=woi.width*oY;
        							for (int j=0;j<this.overlapStep;j++){
        								int oX=(tileX*this.overlapStep-woi.x)+j;
        								if ((oX>=0) && (oX<woi.width)) result[oIndex+oX]=tile[iIndex];
        								iIndex++;
        							}
        						}
        					}
        				}
        			}
        		}
        		return result;
        	}
        	public float [] renderZMapOverlap (
        			int nImg,
        			Rectangle woi,
        			int auxLayer){
        		float [] result = new float [woi.width*woi.height];
        		int yBottom=woi.y+woi.height;
        		int xRight= woi.x+woi.width;
        		int tileYmin=woi.y/this.overlapStep;
        		int tileYmax=(yBottom-1)/this.overlapStep+1; // +1 to the last needed
        		if (tileYmax>this.tilesY) tileYmax=this.tilesY;
        		int tileXmin=woi.x/this.overlapStep;
        		int tileXmax=(xRight-1)/this.overlapStep+1; // +1 to the last needed
        		if (tileXmax>this.tilesX) tileXmax=this.tilesX;

        		for (int i=0;i<result.length;i++) result[i]=0.0f;
        		for (int tileY=tileYmin;tileY<tileYmax;tileY++) for (int tileX=tileXmin;tileX<tileXmax; tileX++){
        			ZTile zTile=this.zMap[nImg][tileY][tileX];
        			if (zTile!=null){
        				float [] tile=(auxLayer<0)?zTile.getZmap():zTile.getAux(auxLayer);
        				if (tile!=null) {
        					int tileSize=zTile.getPaddedSize(); //this.overlapStep+2*zTile.overlap;
            				int padding=(tileSize-zTile.getSize())/2;
        					for (int i=0;i<tileSize;i++) {
        						int oY=(tileY*this.overlapStep-woi.y)+i-padding;
        						if ((oY>=0) && (oY<woi.height)){
        							int iIndex=tileSize*i;
        							int oIndex=woi.width*oY;
        							for (int j=0;j<tileSize;j++){
        								int oX=(tileX*this.overlapStep-woi.x)+j-padding;
        								if ((oX>=0) && (oX<woi.width)) result[oIndex+oX]+=tile[iIndex];
        								iIndex++;
        							}
        						}
        					}
        				}
        			}
        		}
        		return result;
        	}


        	private void setupZMapTile(
        			int tileX,
        			int tileY,
        			float [] thisDisparityTiles,
        			ZTile [][] thisZMap,
        			//        			int nImg,
        			int maxNumber,
        			double minFirst,
        			double minAbsolute,
        			double minRelative,
        			double mergeMax,
        			int overlap,
        			double zMapMinForeground,
        			int debugLevel
        	){
        		if (debugLevel>2) System.out.println("setupZMapTile("+tileX+","+tileY+",...)");
        		boolean debugThisTile=debugLevel>3;
        		thisZMap[tileY][tileX]=new ZTile();
        		ZTile zTile=thisZMap[tileY][tileX];
        		double [] correlationSection = getTileDisparity(tileX, tileY, thisDisparityTiles);
        		//					double threshold=minAbsolute;
        		boolean [] isMax =new boolean[correlationSection.length];
        		int iMax=-1;
        		int numMax=0;
        		for (int i=0;i<correlationSection.length;i++){
        			if ((correlationSection[i]>=minAbsolute) &&
        					((i==0)||(correlationSection[i]>=correlationSection[i-1])) &&
        					((i==(correlationSection.length-1)) || (correlationSection[i]>=correlationSection[i+1]))){
        				isMax[i]=true; // local max above absolute threshold
        				if ((numMax==0) || (correlationSection[i]>correlationSection[iMax])) iMax=i;
        				numMax=1;
        			} else {
        				isMax[i]=false;
        			}
        		}
        		//						double mergeMax=1.5; //merge maximums for interpolation if closer than mergeMax;
        		int iMergeMax=(int) Math.round(mergeMax/this.disparityPerEntry);
        		if (debugThisTile){
        			System.out.println("mergeMax="+mergeMax+" iMergeMax="+iMergeMax);
        		}
        		if ((numMax>0) && (correlationSection[iMax]<minFirst)) numMax=0; // no maximums if the largest is below thershold
        		// for now - now maximums above thereshold - use infinity (disparity 0)
        		if (numMax==0) {
        			isMax[0]=true; // infinity
        			numMax=1;
        			iMax=0;
        		} else {
        			numMax=0;
//        			double threshold=correlationSection[iMax]*minRelative; // no relative filtering here - will be done later
        			for (int i=0;i<correlationSection.length;i++) if (isMax[i]){
//        				if (correlationSection[i]>=threshold) numMax++;
        				numMax++;
        				//else isMax[i]=false;
        			}
        		}
        		if (numMax>maxNumber) numMax=maxNumber; // limit to specified number of correlation maximums to subpixel
        		int [] maxIndices=new int [numMax];
        		maxIndices[0]=iMax;
        		isMax[iMax]=false;
        		int maxNum;
        		for (maxNum=1;maxNum<maxIndices.length;maxNum++){
        			// merge previous one if possible
        			int nearDown=-1, nearUp=-1;
        			for (int i=0;i<=iMergeMax;i++){
        				if ((maxIndices[maxNum-1]-i)<0) break;
        				if (isMax[maxIndices[maxNum-1]-i]){
        					nearDown = maxIndices[maxNum-1]-i;
        					break;
        				}
        			}
        			for (int i=0;i<=iMergeMax;i++){
        				if ((maxIndices[maxNum-1]+i)>=isMax.length) break;
        				if (isMax[maxIndices[maxNum-1]+i]){
        					nearUp = maxIndices[maxNum-1]+i;
        					break;
        				}
        			}
        			int n=1+((nearDown>=0)?1:0)+((nearUp>=0)?1:0);
        			if (n>1){
        				double s=maxIndices[maxNum-1]+((nearDown>=0)?nearDown:0)+((nearUp>=0)?nearUp:0);
        				int iMerged=(int) Math.round(s/n);

        				if (debugThisTile){
        					System.out.println("Merging close maximums: "+maxIndices[maxNum-1]+" with "+
        							((nearDown>=0)?nearDown:"")+" "+((nearUp>=0)?nearUp:"")+" to "+iMerged);
        				}
        				maxIndices[maxNum-1]=iMerged;
        				isMax[iMerged]=false;
        				if (nearDown>=0) isMax[nearDown]=false;
        				if (nearUp>=0)   isMax[nearUp]=false;
        			}
        			iMax=-1;
        			for (int i=0;i<correlationSection.length;i++) {
        				if (isMax[i]&&((iMax<0) || (correlationSection[i]>=correlationSection[iMax]))) iMax=i;
        			}
        			if (iMax<0) break; // no more maximums
        			isMax[iMax]=false;
        			maxIndices[maxNum]=iMax;
        		}
        		if (debugThisTile){
        			System.out.println("List maximums on correlation section");
        			for (int n=0;n<maxNum;n++) {
        				System.out.println(n+" "+maxIndices[n]+"("+(this.disparityPerEntry*maxIndices[n])+" pix) "+correlationSection[maxIndices[n]]);
        			}
        		}
        		//TODO: Always add infinity?
        		//TODO: use quadratic interpolation?
        		// reorder maximums by decreasing disparity (they are now ordered by strength
        		for (boolean ordered=false;!ordered;){
        			ordered=true;
        			for (int i=0;i<(maxIndices.length-1);i++) if (maxIndices[i]<maxIndices[i+1]){
        				ordered=false;
        				int tmp=maxIndices[i];
        				maxIndices[i]=maxIndices[i+1];
        				maxIndices[i+1]=tmp;
        			}
        		}
        		zTile.numPlanes=maxIndices.length;
        		zTile.maximums=new double [zTile.numPlanes][2];
        		for (int i=0;i<zTile.numPlanes;i++){
        			zTile.maximums[i][0]= (float) (maxIndices[i]*this.disparityPerEntry); // disparity in pixels
        			zTile.maximums[i][1]= (float) correlationSection[maxIndices[i]];      // strength
        		}
        		//overlap
        		zTile.overlap=overlap;
        		zTile.size=this.overlapStep;
        		zTile.reset(zMapMinForeground, minAbsolute, minRelative);
        		zTile.setInnerMask(this.innerMask);
        		zTile.setBorderMask(this.borderMask);
        	}

        	private class ZTile{
        		public int numPlanes; //=0;
        		public int size;
        		public double foregroundThreshold; // current foreground threshold (may reconsider if occlusion happens)
    	    	public double minAbsolute=0;
    	    	public double minRelative=0;

        		public int foregroundIndex; // current foreground index;
        		public double [][] maximums; //={}; // {position, strength} in the ordcer of decreasing disparity
        		public int overlap; //=0;
        		public float [] zmap; //=null; // square, centered in the center of the tile - may have margins - final z-value for each pixel
        		public boolean [] enabledPlane;
        		public BitSet [] enabledPixels; //=null;
        		public BitSet [] certainPixels; //=null;
        		public float [][] likely; //=null;
        		public float [][] unlikely; //=null;
        		public float [][] auxData=null;
        		private BitSet innerMask=null;
        		private int []borderMask=null;
        		public int getNumberOfPlanes(){return this.maximums.length;}
        		@SuppressWarnings("unused")
				public int getForegroundPlane(){return this.foregroundIndex;}
        		public double getPlaneDisparity (int plane){return this.maximums[plane][0];}
        		public void setPlaneDisparity (double disparity, int plane){this.maximums[plane][0]=disparity;}
        		public double getPlaneStrength (int plane){return this.maximums[plane][1];}
        		public boolean getPlaneEnabled (int plane){return this.enabledPlane[plane];}
        		public int getSize(){return this.size;}
        		public int getPaddedSize(){return this.size+2*this.overlap;}
        		public int getOverlap(){return this.overlap;}
        		public int getPaddedLength(){return (this.size+2*this.overlap)*(this.size+2*this.overlap);}
        		public void setInnerMask(BitSet mask){
        			this.innerMask=mask;
        		}
        		public void setBorderMask(int [] mask){
        			this.borderMask=mask;
        		}
        		public int [] getBorderMask(){return this.borderMask;}
        		public void setForegroundPlane(int plane){
    				this.foregroundIndex=plane;
    				if (this.maximums[plane][1]<this.foregroundThreshold) this.foregroundThreshold=this.maximums[plane][1];
    				this.enabledPlane[plane]=true;
        		}
        		public void disableForegroundPixel(int i){
        			if (this.enabledPixels[this.foregroundIndex]== null) {
            			int paddedSize=this.size+2*this.overlap;
            			int paddedLength=paddedSize*paddedSize;
        				this.enabledPixels[this.foregroundIndex]=new BitSet(paddedLength);
        				this.enabledPixels[this.foregroundIndex].set(0,paddedLength);
        			}
        			this.enabledPixels[this.foregroundIndex].clear(i);
        		}
        		public boolean isEnabledForegroundPixel(int i){
        			if (this.enabledPixels[this.foregroundIndex]== null) return true;
        			return this.enabledPixels[this.foregroundIndex].get(i);
        		}
        		public void initAux(int n){
        			this.auxData=new float[n][];
        			for (int i=0;i<n;i++) this.auxData[i]=null;
        		}
        		@SuppressWarnings("unused")
				public void resetAux(){
        			this.auxData=null;
        		}
        		public void setAux(int n, float [] data){
        			if (this.auxData==null) initAux(n+1);
        			this.auxData[n]=data;
        		}
        		@SuppressWarnings("unused")
				public float [][] getAux (){
        			return this.auxData;
        		}
        		public float [] getAux (int n){
        			if ((this.auxData==null) || (this.auxData.length<=n)) return null;
        			return this.auxData[n];
        		}

        		public void initLikely(){
        			this.likely=new float[this.numPlanes][];
        			for (int i=0;i<this.numPlanes;i++) this.likely[i]=null;
        		}
        		@SuppressWarnings("unused")
				public void resetLikely(){
        			this.likely=null;
        		}
        		@SuppressWarnings("unused")
				public void setLikely(int n, float [] data){
        			if (this.likely==null) initLikely();
        			this.likely[n]=data;
        		}
        		@SuppressWarnings("unused")
				public float [][] getLikely (){
        			return this.likely;
        		}
        		@SuppressWarnings("unused")
				public float [] getLikely (int n){
        			if ((this.likely==null) || (this.likely.length<=n)) return null;
        			return this.likely[n];
        		}

        		public void initUnlikely(){
        			this.unlikely=new float[this.numPlanes][];
        			for (int i=0;i<this.numPlanes;i++) this.unlikely[i]=null;
        		}
        		@SuppressWarnings("unused")
				public void resetUnlikely(){
        			this.unlikely=null;
        		}
        		@SuppressWarnings("unused")
				public void setUnlikely(int n, float [] data){
        			if (this.unlikely==null) initUnlikely();
        			this.unlikely[n]=data;
        		}
        		@SuppressWarnings("unused")
				public float [][] getUnlikely(){
        			return this.unlikely;
        		}
        		@SuppressWarnings("unused")
				public float [] getUnlikely(int n){
        			if ((this.unlikely==null) || (this.unlikely.length<=n)) return null;
        			return this.unlikely[n];
        		}





        		public void initEnabledForegroundPixelsPlane(){
        			if (this.foregroundThreshold>=this.enabledPixels.length) return; // or make it an error?
        			int paddedSize= this.size+this.overlap*2;
        			this.enabledPixels[this.foregroundIndex]=new BitSet(paddedSize*paddedSize);
        			this.enabledPixels[this.foregroundIndex].set(0,paddedSize*paddedSize);
        			this.enabledPlane[this.foregroundIndex]=true;
        		}
        		public boolean isForegroundValid(){
        			if (this.foregroundIndex>=this.maximums.length) return false; // no foreground
        			if (!this.enabledPlane[this.foregroundIndex]) return false; // disabled plane
        			if (this.enabledPixels[this.foregroundIndex]== null) return true; // all pixels enabled
        			BitSet bb=(BitSet) this.enabledPixels[this.foregroundIndex].clone();
        			bb.and(this.innerMask);
        			if (!bb.isEmpty()) return true;
        			return false;
        		}
        		public boolean advanceForeground(){
        			while (!isForegroundValid()){
        				if (this.foregroundIndex>=this.maximums.length) return false; // no foreground left
        				this.foregroundIndex++;
        			}
        			return true;
        		}
        		public float [] getZmap(){
        			if (this.zmap!=null) return this.zmap;
        			int paddedSize=this.size+2*this.overlap;
        			int paddedLength=paddedSize*paddedSize;
        			float [] map= new float [paddedLength];
        			BitSet needed=null;
        			BitSet newBits=null;
        			for (int plane=this.foregroundIndex; plane<this.maximums.length;plane++) if (this.enabledPlane[plane]){
        				if ((this.enabledPixels[plane]==null) && (needed==null)){
        					for (int i=0;i<paddedLength;i++) map[i]=(float) this.maximums[plane][0];
        					return map;
        				}
        				if (needed==null){
        					needed=new BitSet(paddedLength);
        					needed.set(0,paddedLength);
        				}
        				if (this.enabledPixels[plane]==null) { // null - all enabled
        					newBits= new BitSet(paddedLength);
        					newBits.set(0,paddedLength);
        				} else newBits=(BitSet) this.enabledPixels[plane].clone();
        				newBits.and(needed);
        				for (int i=newBits.nextSetBit(0); i>=0;i=newBits.nextSetBit(i+1)){
        					map[i]=(float) this.maximums[plane][0];
        				}
        				needed.andNot(newBits);
        				if (this.enabledPixels[plane]==null) return map;
        				if (needed.isEmpty()) return map;
        			}
        			if (needed==null){ // no maximums at all, probably
        				for (int i=0;i<paddedLength;i++) map[i]=0.0f;
        				return map;
        			}
        			if (!needed.isEmpty()) for (int i=needed.nextSetBit(0); i>=0;i=needed.nextSetBit(i+1))map[i]=0.0F;
        			return map;
        		}
        		/*
        		public boolean [] getNonOccluded(double disparity){
        			int paddedSize=this.size+2*this.overlap;
        			int paddedLength=paddedSize*paddedSize;
        			BitSet nonOccludedBits=new BitSet(paddedLength);
        			nonOccludedBits.set(0,paddedLength);
        			for (int plane=0;(plane<getNumberOfPlanes()) && (getPlaneDisparity(plane)>disparity);plane++) if (this.enabledPlane[plane]) {
        				if (this.certainPixels[plane]!=null) nonOccludedBits.andNot(this.certainPixels[plane]);
        			}
        			boolean [] nonOccluded=new boolean[paddedLength];
        			for (int i=0;i<paddedLength;i++) nonOccluded[i]=nonOccludedBits.get(i);
        			return nonOccluded;
        		}
        		*/
        		/**
        		 * Filter pixels by occlusion of the foreground and disabled in the current plane
        		 * @param disparityFG - threshold to consider pixel being in front
        		 * @param disparity - target disparity
        		 * @param disparityTolerance - consider "this" as being withing disparityTolerance of disparity. NaN - do not filter by this
        		 * @return
        		 */
        		@SuppressWarnings("unused")
				public boolean [] getEnabledNonOccluded(
        				double disparityFG,
        				double disparity,
        				double disparityTolerance){
            		return getEnabledNonOccluded(
            				disparityFG,
            				disparity,
            				disparityTolerance,
            				0);
        		}
        		public boolean [] getEnabledNonOccluded(
        				double disparityFG,
        				double disparity,
        				double disparityTolerance,
        				int debugLevel){
        			if (debugLevel>3) System.out.println(" ---- getEnabledNonOccluded("+disparityFG+","+
        				disparity+","+
        				disparityTolerance+","+
        				debugLevel+")");
        			int paddedSize=this.size+2*this.overlap;
        			int paddedLength=paddedSize*paddedSize;
        			BitSet nonOccludedBits=new BitSet(paddedLength);
        			nonOccludedBits.set(0,paddedLength);
        			for (int plane=0;(plane<getNumberOfPlanes()) && (getPlaneDisparity(plane)>disparityFG);plane++) if (this.enabledPlane[plane]) {
        				if (this.certainPixels[plane]!=null) nonOccludedBits.andNot(this.certainPixels[plane]);
            			if (debugLevel>3) System.out.println(" ------ plane="+plane+" cumul. nonOccludedBits.cardinality()="+nonOccludedBits.cardinality());
        			}
        			if (!Double.isNaN(disparityTolerance)){
        				BitSet enabledBits=new BitSet(paddedLength);
        				for (int plane=0;plane<getNumberOfPlanes();plane++)
        					if (this.enabledPlane[plane] && (Math.abs(getPlaneDisparity(plane)-disparity)<=disparityTolerance)) {
        						if (this.enabledPixels[plane]!=null) enabledBits.or(this.enabledPixels[plane]);
        						else enabledBits.set(0,paddedLength); // all enabled
                    			if (debugLevel>3) System.out.println(" ------ plane="+plane+" cumul. enabledBits.cardinality()="+enabledBits.cardinality());
        					}
        				nonOccludedBits.and(enabledBits);
            			if (debugLevel>3) System.out.println(" ------ result nonOccludedBitscardinality()="+nonOccludedBits.cardinality());
        			}
        			boolean [] nonOccluded=new boolean[paddedLength];
        			for (int i=0;i<paddedLength;i++) nonOccluded[i]=nonOccludedBits.get(i);
        			return nonOccluded;
        		}
            	public void reset(
            			double minForeground,
            	    	double minAbsolute,
            	    	double minRelative){
            		this.numPlanes=this.maximums.length;
            		this.enabledPixels=new BitSet[this.numPlanes];
            		this.certainPixels=new BitSet[this.numPlanes];
            		this.enabledPlane=new boolean[this.numPlanes];
            		for (int i=0;i<this.numPlanes;i++){
            			this.enabledPixels[i]= null;
            			this.certainPixels[i]=null;
            			this.enabledPlane[i]=true;
            		}
            		this.zmap=null;
            		setMinCorrelations(minAbsolute, minRelative,false);
            		setForeGround(minForeground);
            		this.likely=new float[this.numPlanes][]; // maybe use later
            		for (int i=0;i<this.likely.length;i++)this.likely[i]=null;
            		this.unlikely=new float[this.numPlanes][]; // maybe use later
            		for (int i=0;i<this.unlikely.length;i++)this.unlikely[i]=null;
            	}
            	public void setMinCorrelations(
            	    	double minAbsolute,
            	    	double minRelative,
            	    	boolean keepFG){
            		if (!Double.isNaN(minAbsolute))   this.minAbsolute=minAbsolute;
            		if (!Double.isNaN(minRelative))   this.minRelative=minRelative;
            		double aMax=0;
            		for (int i=0;i<this.numPlanes;i++) if (this.maximums[i][1]>aMax) aMax=this.maximums[i][1];
            		aMax=Math.max(aMax*this.minRelative,this.minAbsolute);
            		for (int i=0;i<this.numPlanes;i++)	this.enabledPlane[i]=this.maximums[i][1]>=aMax;
            		if (keepFG && (this.foregroundIndex<this.enabledPlane.length) && (this.foregroundIndex>=0)) {
            			this.enabledPlane[this.foregroundIndex]=true;
            		}
            	}
            	public void setForeGround(
            			double minForeground){
            		if (!Double.isNaN(minForeground)) this.foregroundThreshold=minForeground;
            		for (this.foregroundIndex=0;
          		  (this.foregroundIndex<this.maximums.length) && (!this.enabledPlane[this.foregroundIndex] || (this.maximums[this.foregroundIndex][1]<this.foregroundThreshold))
          		  ;this.foregroundIndex++);
            	}

        	} // end of private class ZTile

        	public class Photometric{
        		public String [] channelNames={"Y","Cb","Cr","Aux"};
        		public double  [][] valueLimits=null; // lo/high limits for each channel
        		public int subdivAverage=256;
        		public int subdivHalfDifference=128; // 257
//        		public int subdivVariance = 256;
        		public double smoothVarianceSigma=10.0;
        		public double scaleVariance=3.0; // sigma along difference (vertical axis) will be scaleVariance*variance for this average value
        		public int numImages=0;
        		public double [][][] standardVariance=null;// for each image each channel - average variance (among 9 neighbors) for different values
        		public double [][] averageVariance=null;
        		public double [][][][] matchingQuality=null; // for each image, each second image index, each channel [subdivAverage*subdivDifference] <=1.0 values
        		                                             // horizontal average, vertical - difference "2-1" (for "1-2" flip vertical)
        		public int histogramSize=1000;
        		public double getAverageVariance(int nImg,int chn){
        			return this.averageVariance[nImg][chn];
        		}
        		public double [] initStaging(){
        			double [] staging=new double[this.subdivAverage*(2*this.subdivHalfDifference+1)];
        			for (int i=0;i<staging.length;i++) staging[i]=0.0;
        			return staging;
        		}
        		public void addStagingSample(
        				double []staging,
        				int chn,
        				double weight,
        				double v1,
        				double v2){
        			double min=this.valueLimits[chn][0];
        			double step=(this.valueLimits[chn][1]-this.valueLimits[chn][0])/(this.subdivAverage-1);
        			double av=0.5*(v1+v2);
					int countAvg=(int) Math.round((av-min)/step);
					if (countAvg<0) countAvg=0;
					if (countAvg>=this.subdivAverage) countAvg=this.subdivAverage-1;
					int countDiff= (int)Math.round(0.5*(v2-v1))+subdivHalfDifference;
					if (countDiff<0) countDiff=0;
					if (countDiff>2*this.subdivHalfDifference) countDiff=2*this.subdivHalfDifference;
					staging[countDiff*this.subdivAverage+countAvg]+=weight;
        		}
        		public void addStagingSample(
        				double []staging,
        				int chn,
        				double weight,
        				double v1,
        				double v2,
        				// reduce weight depending on difference, scale to measured variance
        				double scaleVariance,
        				double kLocal, // 0 - use global varaiance, 1.0 - use local
        				int nImg1, // first image number
        				int nImg2, // second image number
        				double var1, // first image variance
        				double var2, // second image variance
        				boolean debug){
        			double min=this.valueLimits[chn][0];
        			double step=(this.valueLimits[chn][1]-this.valueLimits[chn][0])/(this.subdivAverage-1);
        			double av=0.5*(v1+v2);
					int countAvg=(int) Math.round((av-min)/step);
					if (countAvg<0) countAvg=0;
					if (countAvg>=this.subdivAverage) countAvg=this.subdivAverage-1;
					double diff=0.5*(v2-v1);
					int countDiff= (int)Math.round(diff/step)+subdivHalfDifference;
					if (countDiff<0) countDiff=0;
					if (countDiff>2*this.subdivHalfDifference) countDiff=2*this.subdivHalfDifference;
					if (debug){
						System.out.println("addStagingSample(...,"+chn+","+weight+","+v1+","+v2+","+scaleVariance+","+
								kLocal+","+nImg1+","+nImg2+","+var1+","+var2);
						System.out.println(" +++ min="+min+" step="+step+" av="+av+" diff="+diff+" countAvg="+countAvg+" countDiff="+countDiff);
					}
					if (scaleVariance>0.0){
						if (kLocal<1.0){ // mix local variance with average over all image
							int count1=(int) Math.round((v1-min)/step);
							if (count1<0) count1=0;
							if (count1>=this.subdivAverage) count1=this.subdivAverage-1;
							int count2=(int) Math.round((v2-min)/step);
							if (count2<0) count2=0;
							if (count2>=this.subdivAverage) count2=this.subdivAverage-1;
							if (kLocal<0) kLocal=0;
							var1=kLocal*var1+(1.0-kLocal)*this.standardVariance[nImg1][chn][count1];
							var2=kLocal*var2+(1.0-kLocal)*this.standardVariance[nImg2][chn][count2];
						}
						double sigma=scaleVariance*Math.sqrt(var1*var1);
						weight*=Math.exp(-diff*diff/(2.0*sigma*sigma))/sigma/Math.sqrt(2*Math.PI);
						if (debug)	System.out.println(" +++ sigma="+sigma+" weight="+weight);
					}
					staging[countDiff*this.subdivAverage+countAvg]+=weight;
					if (debug)	System.out.println("     staging["+(countDiff*this.subdivAverage+countAvg)+"]="+staging[countDiff*this.subdivAverage+countAvg]);

        		}
        		public double getStrength(
        				double []staging,
        				int chn,
        				double v1,
        				double v2){
        			double min=this.valueLimits[chn][0];
        			double step=(this.valueLimits[chn][1]-this.valueLimits[chn][0])/(this.subdivAverage-1);
        			double av=0.5*(v1+v2);
					int countAvg=(int) Math.floor((av-min)/step);
					int countAvg1=(int) Math.ceil((av-min)/step);
					double dx=(av-min)-step*countAvg;
					if (countAvg<0) {
						countAvg=0;
						countAvg1=0;
						dx=0.0;
					}
					if (countAvg1>=this.subdivAverage){
						countAvg=this.subdivAverage-1;
						countAvg1=this.subdivAverage-1;
						dx=0.0;
					}
					double diff=0.5*(v2-v1);
					int countDiff= (int)Math.floor(diff/step)+this.subdivHalfDifference;
					int countDiff1= (int)Math.ceil(diff/step)+this.subdivHalfDifference;
					double dy=diff-step*(countDiff-this.subdivHalfDifference);
					if (countDiff<0) {
						countDiff=0;
						countDiff1=0;
						dy=0.0;
					}
					if (countDiff1>2*this.subdivHalfDifference){
						countDiff=2*this.subdivHalfDifference;
						countDiff1=2*this.subdivHalfDifference;
						dy=0.0;
					}
					if ((countDiff1*this.subdivAverage+countAvg1)>=staging.length){
						System.out.println("BUG: getStrength() countAvg="+countAvg+" countAvg1="+countAvg1+
								" countDiff="+countDiff+" countDiff1="+countDiff1+" staging.length="+staging.length);
					}
					return
					staging[countDiff *this.subdivAverage+countAvg ]*(1-dy)*(1-dx)+
					staging[countDiff *this.subdivAverage+countAvg1]*(1-dy)*(  dx)+
					staging[countDiff1*this.subdivAverage+countAvg ]*(  dy)*(1-dx)+
					staging[countDiff1*this.subdivAverage+countAvg1]*(  dy)*(  dx);
        		}

        		public void blurStaging(
        				double []staging,
        				int chn,
        				double sigma){
        			if (sigma<=0) return;

//        			double min=this.valueLimits[chn][0];
        			double step=(this.valueLimits[chn][1]-this.valueLimits[chn][0])/(this.subdivAverage-1);

       				(new DoubleGaussianBlur()).blurDouble(
       						staging,
       						this.subdivAverage,
       						2*this.subdivHalfDifference+1,
    						sigma/step,
    						sigma/step,
    						0.01);
        		}
        		public  Photometric(
        				double [][][] images,
        				int imageWidth,
        				int margins,
        				double ignoreFraction,
                		int subdivAverage,
                		int subdivHalfDifference,
//                		int subdivVariance,
                		double smoothVarianceSigma,
                		double scaleVariance,
        				int debugLevel
        				){
        			if (margins<1) margins=1;
            		this.subdivAverage=subdivAverage;
            		this.subdivHalfDifference=subdivHalfDifference;
//            		this.subdivVariance = subdivVariance;
            		this.smoothVarianceSigma=smoothVarianceSigma;
            		this.scaleVariance=scaleVariance; // sigma along difference (vertical axis) will be scaleVariance*variance for this average value
            		this.numImages=images.length;
            		this.valueLimits=new double[channelNames.length][];
        			int [] dirs1={1,imageWidth+1,imageWidth,imageWidth-1,-1,-imageWidth-1,-imageWidth,-imageWidth+1,0};
        			this.standardVariance=new double[this.numImages][this.channelNames.length][];
        			this.averageVariance=new double[this.numImages][this.channelNames.length];
        			for (int i=0;i<this.standardVariance.length;i++) for (int j=0;j<this.standardVariance[0].length;j++){
        				this.standardVariance[i][j]=null;
        				this.averageVariance[i][j]=0.0;
        			}

            		for (int chn=0;chn<this.valueLimits.length;chn++){
            			this.valueLimits[chn]=null;
            			double min=Double.NaN, max=Double.NaN;
            			for (int nImg=0;nImg<this.numImages;nImg++) if ((images[nImg]!=null) && (images[nImg][chn]!=null)){
            				double [] image=images[nImg][chn];
            				int imageHeight=image.length/imageWidth;
            				if (Double.isNaN(min)){
            					min=image[margins*imageWidth+margins];
            					max=min;
            				}
                			if (debugLevel>2){
                				System.out.println("nImg="+nImg+" chn="+chn+
                						" min0="+IJ.d2s(min,3)+" max0="+IJ.d2s(max,3));
                			}
            				for (int y=margins;y<(imageHeight-margins);y++) for (int x=margins;x<(imageWidth-margins);x++){
            					if (image[y*imageWidth+x]<min) {
            						min=image[y*imageWidth+x];
                        			if (debugLevel>2)	System.out.println("y="+y+" x="+x+" min="+IJ.d2s(min,3));

            					}
            					if (image[y*imageWidth+x]>max) {
            						max=image[y*imageWidth+x];
                        			if (debugLevel>2)	System.out.println("y="+y+" x="+x+" max="+IJ.d2s(max,3));
            					}
            				}
                			if (debugLevel>2){
                				System.out.println("nImg="+nImg+" chn="+chn+
                						" min="+IJ.d2s(min,3)+" max="+IJ.d2s(max,3));
                			}

            			}
            			int [] histogram=new int [this.histogramSize];
            			for (int i=0;i<histogram.length;i++) histogram[i]=0;
            			double step=(max-min)/(this.histogramSize-0.0001);
            			if (debugLevel>2){
            				System.out.println(	" min="+IJ.d2s(min,3)+" max="+IJ.d2s(max,3)+" step="+IJ.d2s(step,6));
            			}

            			for (int nImg=0;nImg<this.numImages;nImg++) if ((images[nImg]!=null) && (images[nImg][chn]!=null)){
            				double [] image=images[nImg][chn];
            				int imageHeight=image.length/imageWidth;
            				for (int y=margins;y<(imageHeight-margins);y++) for (int x=margins;x<(imageWidth-margins);x++){
            					int index=(int) Math.floor((image[y*imageWidth+x]-min)/step);
            					if (index<0) index=0;
            					if (index>=this.histogramSize) index=this.histogramSize-1;
            					histogram[index]++; //java.lang.ArrayIndexOutOfBoundsException: 1005
            				}
            			}
            			int totalNum=0;
            			for (int i=0;i<histogram.length;i++)totalNum+= histogram[i];
            			int ignoreNumPix= (int)(Math.floor(totalNum*ignoreFraction));
            			if (debugLevel>2){
            				System.out.println(	" totalNum="+totalNum+" ignoreNumPix="+ignoreNumPix);
            			}

            			this.valueLimits[chn]=new double[2];
            			int num=0;
            			for (int i=0;i<histogram.length;i++){
            				num+=histogram[i];
                			if (debugLevel>2)	System.out.println("---  "+num+" min+"+i+"*step="+IJ.d2s(min+i*step,3));
            			    if (num>=ignoreNumPix){
            			    	this.valueLimits[chn][0]=min+i*step;
                    			if (debugLevel>2){
                    				System.out.println("i="+i+" min+i*step="+IJ.d2s(min+i*step,3));
                    			}
            			    	break;
            			    }
            			}
            			num=0;
            			for (int i=histogram.length-1;i>=0;i--){
            				num+=histogram[i];
                			if (debugLevel>2)	System.out.println("---  "+num+" min+"+i+"*step="+IJ.d2s(min+i*step,3));
            			    if (num>=ignoreNumPix){
            			    	this.valueLimits[chn][1]=min+i*step;
                    			if (debugLevel>2){
                    				System.out.println("i="+i+" min+i*step="+IJ.d2s(min+i*step,3));
                    			}
            			    	break;
            			    }
            			}
            			if (debugLevel>1){
            				System.out.println("Channel '"+this.channelNames[chn]+
            						"' min="+IJ.d2s(this.valueLimits[chn][0],3)+" max="+IJ.d2s(this.valueLimits[chn][1],3));
            			}
            			// calculatye variance for each emage/channel
            			for (int nImg=0;nImg<this.numImages;nImg++) if ((images[nImg]!=null) && (images[nImg][chn]!=null)){
            				double [] image=images[nImg][chn];
            				int imageHeight=image.length/imageWidth;
            				this.standardVariance[nImg][chn]=new double[this.subdivAverage];
            				int [] samples=new int [this.subdivAverage];
            				int totalSamples=0;
            				for (int i=0;i<this.subdivAverage;i++){
            					this.standardVariance[nImg][chn][i]=0.0;
            					samples[i]=0;
            				}
            				double stepVar=(this.valueLimits[chn][1]-this.valueLimits[chn][0])/(this.subdivAverage-1);
            				double minVar=this.valueLimits[chn][0];
            				for (int y=margins;y<(imageHeight-margins);y++) for (int x=margins;x<(imageWidth-margins);x++){
            					int index=y*imageWidth+x;
								double  S1=0.0;
								double  S2=0.0;
								for (int d=0;d<dirs1.length;d++){
									int indexOther=index+dirs1[d];
									S1+=image[indexOther];
									S2+=image[indexOther]*image[indexOther];
								}
								double v2=(S2-(S1*S1)/dirs1.length)/dirs1.length; //dirs1.length;
								double avg=S1/dirs1.length;
								int count=(int) Math.round((avg-minVar)/stepVar);
								if (count<0) count=0;
								if (count>=this.subdivAverage) count=this.subdivAverage-1;
								this.standardVariance[nImg][chn][count]+=v2; // squared
								samples[count]++;
								this.averageVariance[nImg][chn]+=v2;
								totalSamples++;
            				}
            				for (int i=0;i<this.subdivAverage;i++){
            					if (samples[i]>0) this.standardVariance[nImg][chn][i]=Math.sqrt(this.standardVariance[nImg][chn][i]/samples[i]);
            					else this.standardVariance[nImg][chn][i]=0.0;
//                    			if ((debugLevel>1) && (nImg==0) && (chn==0)){
                        			if ((debugLevel>2) && (nImg==0)){
                    				System.out.println("samples["+i+"]="+samples[i]+
                    						" this.standardVariance["+nImg+"]["+chn+"]["+i+"]="+IJ.d2s(this.standardVariance[nImg][chn][i],3));
                    			}

            				}
            				this.averageVariance[nImg][chn]=Math.sqrt(this.averageVariance[nImg][chn]/totalSamples);
            				for (int i=0;i<this.subdivAverage;i++) if (samples[i]==0){
            					int j=i+1;
            					for (;(j<this.subdivAverage) && (samples[j]==0);j++);
            					if (j==this.subdivAverage) j--;
            					int i0=i-1;
            					if (i0<0) i0=0;
            					if (samples[i0]==0) this.standardVariance[nImg][chn][i0]=this.standardVariance[nImg][chn][j];
            					if (samples[j]==0)  this.standardVariance[nImg][chn][j]=this.standardVariance[nImg][chn][i0];
            					if (i0<j){
            						double a=(this.standardVariance[nImg][chn][j]-this.standardVariance[nImg][chn][i0])/(j-i0);
            						for (int k=0;k<(j-i0);k++){
            							this.standardVariance[nImg][chn][i0+k]=this.standardVariance[nImg][chn][i0]+a*k;
            						}
            						i=j+1;
            					}
            				}
            				// fill gaps (if any)
            				if (this.smoothVarianceSigma>0.0){
            					(new DoubleGaussianBlur()).blur1Direction(
            							this.standardVariance[nImg][chn], //double [] pixels,
            							this.subdivAverage, //int        width,
            							1, //int       height,
            							this.smoothVarianceSigma, //double     sigma,
            							0.01, //double   accuracy,
            							true //boolean xDirection: true - horizontal
            					);
            				}
            			}
            		}
        			// now cteate matchingQuality arrays for each image pair/channel
//            		this.matchingQuality=new double [this.numImages][this.numImages-1][this.valueLimits.length][this.subdivAverage*this.subdivDifference];
            		this.matchingQuality=new double [this.numImages][this.numImages-1][][]; //[this.subdivAverage*this.subdivDifference];
            		int subdivDifference=2*this.subdivHalfDifference+1;
            		int matchingQualityLength=this.subdivAverage*subdivDifference;
            		for (int nImg=0;nImg<this.numImages;nImg++) for (int sIndex=0;sIndex<(this.numImages-1);sIndex++){
            			int sImg=(sIndex>=nImg)?(sIndex+1):sIndex;
            			if ((images[nImg]!=null) && (images[sImg]!=null)){
            				this.matchingQuality[nImg][sIndex]=new double [this.valueLimits.length][];
            				for (int chn=0;chn<this.valueLimits.length;chn++){
            					if ((images[nImg][chn]!=null) && (images[sImg][chn]!=null)){
            						this.matchingQuality[nImg][sIndex][chn]=new double [matchingQualityLength];
            						double diffStep=(this.valueLimits[chn][1]-this.valueLimits[chn][0])/this.subdivHalfDifference;
            						double k=0.5*diffStep*diffStep/(this.scaleVariance*this.scaleVariance);
            						for (int averageIndex=0;averageIndex<this.subdivAverage;averageIndex++){
            							double variance2=this.standardVariance[nImg][chn][averageIndex]*this.standardVariance[sImg][chn][averageIndex];
//        		public double scaleVariance=3.0; // sigma along difference (vertical axis) will be scaleVariance*variance for this average value
            							double a=k/variance2;
            							for (int i=0;i<=this.subdivHalfDifference;i++){
            								double d=Math.exp(-a*i*i);
            							    this.matchingQuality[nImg][sIndex][chn][(this.subdivHalfDifference-i)*this.subdivAverage+averageIndex]=d;
            							    this.matchingQuality[nImg][sIndex][chn][(this.subdivHalfDifference+i)*this.subdivAverage+averageIndex]=d;
            							}
            						}
            					} else {
            						this.matchingQuality[nImg][sIndex][chn]=null;
            					}
            				}
            			} else {
            				this.matchingQuality[nImg][sIndex]=null;
            			}
            		}
        		}
        		public void showmatchingQuality(){
            		int subdivDifference=2*this.subdivHalfDifference+1;
            		int matchingQualityLength=this.subdivAverage*subdivDifference;
            		double [] zero = new double [matchingQualityLength];
            		int numPairs=this.numImages*(this.numImages-1);
            		double [][] debugData=new double [numPairs*this.channelNames.length][];
            		String [] titles=new String [numPairs*this.channelNames.length];
            		int index=0;
            		for (int nImg=0;nImg<this.numImages;nImg++) for (int sIndex=0;sIndex<(this.numImages-1);sIndex++){
            			int sImg=(sIndex>=nImg)?(sIndex+1):sIndex;
            			for (int chn=0;chn<this.channelNames.length;chn++){
            				titles[index]=this.channelNames[chn]+"-"+nImg+"-"+sImg;
            				debugData[index++]=(this.matchingQuality[nImg][sIndex][chn]!=null)?this.matchingQuality[nImg][sIndex][chn]:zero;
            			}

            		}
        			(new ShowDoubleFloatArrays()).showArrays(
        					debugData,
        					this.subdivAverage,
        					subdivDifference,
        					true,
        					"MQ-"+IJ.d2s(this.scaleVariance,2)+"_"+IJ.d2s(this.smoothVarianceSigma,2),
        					titles);
        		}
        	} // end of private class Photometric
    	} // end of public class DisparityTiles

    	/**
    	 * Calculate fine-step (4 pix?) correlation between the pair of square images (32x32) that are already shifterd to compensate for known disparity
    	 * only center area (with 1/4 margins) will be calculated, if the correlation is strong enough, the updated disparity arrays will be calculated for these tiles
    	 *
    	 * @param data a pair of square arrays to be correlated (32x32 pixels)
    	 * @param window cosine mask (8x8=64 pixels long). If the center element is 0.0 - will generate array
    	 * @param dXY unity vector in the direction of the disparity between the two images
    	 * @param phaseCoeff 0.0 - normal correlation, 1.0 - pure phase correlation
    	 * @param disparityArrays tiles of disparity array, top level should be initialized to double [16] [] (for 32x32=1024 input, 8x8 window)
    	 * @param corrMaxDist maximal distance to the correlation maximum
    	 * @param corrThreshold minimal ratio of correlation(0) to avarage correlation to correct disparity
    	 * @param subPixel subdivide grid and result disparityArrays by this number (power of 2)
    	 * @param generateSmooth generate smooth output (otherwise return null)
    	 * @param debugLevel debug level
    	 * @return smooth array, same dimension as data[0], with only center area trusted, proportional to correlation at zero
    	 */
    	public double [] refineCorrelation(
    			DoubleFHT doubleFHT, // will generate if null, can be used to reuse initialization
    			double [][] data,   // difines size
    			double [] window, // defines tile size
    			double [] dXY,    // unity vector defines disparity direction
    			double phaseCoeff, // phase correlation fraction (0 - normal, 1.0 - pure phase)
    			double highPassSigma,
    			double lowPassSigma,
    			double [][] disparityArrays, // top level should be initialized
    			double corrMaxDist, // distance from rhe center to look for the correlation maximum for disparity correction
    			double corrThreshold, // relative correlation contrast to correct disparity
    			int subPixel,       // subdivide pixels when updating disparity
    			boolean generateSmooth,
    			int debugLevel
    			){
    		int size=(int) Math.sqrt(data[0].length);
    		int tileSize=(int) Math.sqrt(window.length);
    		int tileStep=tileSize/2;
    		int tileExtra=size-tileSize;
    		int tileLength=tileSize*tileSize;
    		int tileCenter=(tileSize+1)*tileSize/2;

    		int numTilesRow=size/(2*tileStep);
    		if (doubleFHT==null) doubleFHT= new DoubleFHT();
    		if (Double.isNaN(window[0])){
//    			double [] window1d=doubleFHT.getHamming1d(tileSize,0.0); // pure cosine
    			double [] window1d=doubleFHT.getHamming1d(tileSize); // Hamming
    			for (int i=0;i<window.length;i++) window[i]=window1d[i/tileSize]*window1d[i%tileSize];
    		}
    		int topLeft=(size+1)*(size/4-tileSize/4);
    		double []first= new double[tileLength];
    		double []second=new double[tileLength];
    		double [] corrZero=new double [numTilesRow*numTilesRow];
    		double corrMaxDist2=corrMaxDist*corrMaxDist;
    		int halfDispRange=(int) Math.ceil(corrMaxDist*subPixel);
    		double [] result=null;
    		if (generateSmooth) {
    			result=new double [size*size];
    			for (int i=0;i<result.length;i++) result[i]=0.0;
    		}
    		double [][] debugArray=null;
    		if (debugLevel>2) debugArray=new double[3][data[0].length*subPixel*subPixel];

    		for (int tileY=0;tileY<numTilesRow;tileY++){
    			for (int tileX=0;tileX<numTilesRow;tileX++){
    				int indexIn=topLeft+ tileStep*(tileY*size+ tileX);
    				int indexOut=0;
    				for (int iy=0;iy<tileSize;iy++){
    					for (int ix=0;ix<tileSize;ix++){
    						first [indexOut]=data[0][indexIn];
    						second[indexOut++]=data[1][indexIn++];
    					}
    					indexIn+=tileExtra;
    				}
    				normalizeAndWindowGetDC (first,  window); // DC not needed here
    				normalizeAndWindowGetDC (second, window); // DC not needed here
    				if (debugArray!=null){
    					int debugSize=size*subPixel;
    					int debugTileSize=tileSize*subPixel;
    					int debugTopLeft=debugTileSize*(tileY*debugSize+tileX);
    					for (int debug_i=0;debug_i<debugTileSize;debug_i++) {
    						for (int debug_j=0;debug_j<debugTileSize;debug_j++) {
    							debugArray[1][debugTopLeft+debug_i*debugSize+debug_j]=first[(debug_i/subPixel)*tileSize+(debug_j/subPixel)];
    							debugArray[2][debugTopLeft+debug_i*debugSize+debug_j]=second[(debug_i/subPixel)*tileSize+(debug_j/subPixel)];
    						}

    					}
    				}
    				second=doubleFHT.correlate (second,first,highPassSigma,lowPassSigma,phaseCoeff); // second will be modified
    				double sc=0.0;
    				for (int i=0;i<tileLength;i++)sc+=second[i]*second[i];
    				sc=Math.sqrt(sc/tileLength);
    				if ((debugLevel>3) && (tileY==0) && (tileX==0))	for (int i=0;i<second.length;i++) System.out.println("second["+i+"]="+second[i]);
    				if (debugLevel>2) System.out.println("refineCorrelation(): phaseCoeff="+phaseCoeff+" sc="+sc+" second["+tileCenter+"]="+second[tileCenter]);

    				sc=second[tileCenter]/sc; // realtive correlation strength at zero.
    				int tileIndex=tileY*numTilesRow+tileX;
    				corrZero[tileIndex]=sc; // realtive correlation strength at zero.
    				disparityArrays[tileIndex]=null;
    				if (debugArray!=null){
    					int debugSize=size*subPixel;
    					int debugTileSize=tileSize*subPixel;
    					int debugTopLeft=debugTileSize*(tileY*debugSize+tileX);
    					for (int debug_i=0;debug_i<debugTileSize;debug_i++) {
    						for (int debug_j=0;debug_j<debugTileSize;debug_j++) {
    							debugArray[0][debugTopLeft+debug_i*debugSize+debug_j]=second[(debug_i/subPixel)*tileSize+(debug_j/subPixel)];
    						}

    					}
    				}
    				if (debugLevel>2) System.out.println("refineCorrelation(): sc="+sc+" threshold="+corrThreshold);

    				if (sc>=corrThreshold) {
    					// find absolute maximum on correlation
    					double max=0;
    					int iMax=0;
    					for (int i=0;i<tileLength;i++) if (second[i]>max){
    						max=second[i];
    						iMax=i;
    					}
        				int ixc=iMax%tileSize-tileSize/2;
        				int iyc=iMax/tileSize-tileSize/2;
        				if (debugLevel>2) System.out.println("refineCorrelation(): max="+max+" iMax="+iMax+" corrMaxDist="+corrMaxDist+" corrMaxDist2="+corrMaxDist2+
        						" r2="+(ixc*ixc+iyc*iyc)+" r="+Math.sqrt(ixc*ixc+iyc*iyc));
        				if ((ixc*ixc+iyc*iyc)<=corrMaxDist2){ // maximum close enough
        					// upsample correlation
        					double [] upsampled=doubleFHT.upsample(second,subPixel);
        					int interpolatedSize=tileSize*subPixel;
 //       					int interpolatedLength=interpolatedSize*interpolatedSize;
        		    		int interpolatedCenter=(interpolatedSize+1)*interpolatedSize/2;

        					disparityArrays[tileIndex]=new double[2*halfDispRange+1];
        					for (int i=0;i<2*halfDispRange+1;i++){
        						double dX=dXY[0]*(i-halfDispRange);
        						double dY=dXY[1]*(i-halfDispRange);
        						int iX=(int) Math.floor(dX); // zero in the center
        						int iY=(int) Math.floor(dY);
        						dX-=iX;
        						dY-=iY;
        						int index00= interpolatedCenter+iY*interpolatedSize+iX;
        						int index01=index00+interpolatedSize;
        						int index10=index00+1;
        						int index11=index01+1;
        						disparityArrays[tileIndex][i]= // bi-linear interpolated data
        							(upsampled[index00]*(1.0-dX)+upsampled[index10]*dX)*(1.0-dY)+
        							(upsampled[index01]*(1.0-dX)+upsampled[index11]*dX)*     dY;
        						if (debugLevel>2){
        							System.out.println("disparityArrays["+tileIndex+"]["+i+"]="+disparityArrays[tileIndex][i]+
        									" dX="+dX+" dY="+dY+" iX="+iX+" iY="+iY+" index00="+index00+" index01="+index01+
        									" index10="+index10+" index11="+index11);
        						}
        					}


            				if (debugArray!=null){
            					int debugSize=size*subPixel;
            					int debugTileSize=tileSize*subPixel;
            					int debugTopLeft=debugTileSize*(tileY*debugSize+tileX);
            					for (int debug_i=0;debug_i<debugTileSize;debug_i++) {
            						for (int debug_j=0;debug_j<debugTileSize;debug_j++) {
            							debugArray[0][debugTopLeft+debug_i*debugSize+debug_j]=upsampled[debug_i*interpolatedSize+debug_j];
            						}
            					}
            				}
        					// TODO: update maximum at zero with maximum on the interpolation curve?

        				}
    				}
    				// accumulate result
    	    		if (generateSmooth) {
    	    			int index=topLeft+ tileStep*(tileY*size+ tileX);
        				for (int iy=0;iy<tileSize;iy++){
        					for (int ix=0;ix<tileSize;ix++) result [index++]+=sc*window[iy*tileSize+ix];
        					index+=tileExtra;
        				}
    	    		}
    			}
    		}

    		if (debugArray!=null){
    			String [] dbgTitles={"corr","first","second"};
    			(new ShowDoubleFloatArrays()).showArrays(
    					debugArray,
    					size*subPixel,
    					size*subPixel,
    					true,
    					"correlation-refine",
    					dbgTitles);
    		}
    		return result;
    	}




       	public double [][][][][] calculateDisparityTiles(
    			final String externalTitle,
    			final double [][] externalData, // if not null - process instead of channel data
    			final int channelMask,
    			final int [][] otherImage,
    			final int [][][][] displacementTilesList,
    			final int disparitySteps,
				final double disparityRange,
				final double [][] disparityScales,
				final int interpolationSize,
				final int interpolationUpSample,
    			final int [][][][][][] correlationMap, // [image number][other image index][pixelY][pixelX]{{tiledX,tilledY,pixelIndex}, ...}
    			final int [] minDispalcementTileX,
    			final int [] maxDispalcementTileX,
    			final int [] minDispalcementTileY,
    			final int [] maxDispalcementTileY,
    			final double corrPhaseFraction,
    			final double correlationHighPassSigma,
    			final double correlationLowPassSigma,
    			final double corrCbWeight,
    			final double corrCrWeight,
    			final double[] subpixAMax,
    			final double subpixRMax,
    			final int    subpixNMax,
				final int corrFFTSize,
				final int overlapStep,
				final int numTilesX,
				final int numTilesY,
				final int tileXMin,
				final int tileXMax,
				final int tileYMin,
				final int tileYMax,
    			final int threadsMax,
    			final boolean updateStatus,
				final int debugLevel){
				if (debugLevel>2){
					System.out.println("calculateDisparityTiles(), tileXMin="+tileXMin+", tileXMax="+tileXMax+", tileYMin="+tileYMin+", tileYMax="+tileYMax);
				}
    		final double [][][][][] disparityTiles= new double [numTilesY][numTilesX][][][]; // first [tileY][tileX][image number][other image index][disparity index]
    		final double [][][][] fhtCache=new double [numTilesY][][][];
    		for (int tileY=tileYMin;tileY<=tileYMax;tileY++){
    			for (int tileX=tileXMin;tileX<=tileXMax;tileX++){
    				disparityTiles[tileY][tileX]=null;
    			}
    			fhtCache[tileY]=null;
    		}
    		final int numSensors=this.channel.length;
			int minDispalcementTileXAll=minDispalcementTileX[0];
			int maxDispalcementTileXAll=maxDispalcementTileX[0];
			int minDispalcementTileYAll=minDispalcementTileY[0];
			int maxDispalcementTileYAll=maxDispalcementTileY[0];
			for (int i=1;i<numSensors;i++){
				if (minDispalcementTileXAll>minDispalcementTileX[i]) minDispalcementTileXAll=minDispalcementTileX[i];
				if (maxDispalcementTileXAll<maxDispalcementTileX[i]) maxDispalcementTileXAll=maxDispalcementTileX[i];
				if (minDispalcementTileYAll>minDispalcementTileY[i]) minDispalcementTileYAll=minDispalcementTileY[i];
				if (maxDispalcementTileYAll<maxDispalcementTileY[i]) maxDispalcementTileYAll=maxDispalcementTileY[i];
			}
			if (debugLevel>2){
				System.out.println("calculateDisparityTiles(), minDispalcementTileXAll="+minDispalcementTileXAll+", maxDispalcementTileXAll="+maxDispalcementTileXAll);
				System.out.println("calculateDisparityTiles(), minDispalcementTileXAll="+minDispalcementTileXAll+", maxDispalcementTileXAll="+maxDispalcementTileXAll);
				System.out.println("calculateDisparityTiles(), minDispalcementTileYAll="+minDispalcementTileYAll+", maxDispalcementTileYAll="+maxDispalcementTileYAll);
				System.out.println("calculateDisparityTiles(), minDispalcementTileYAll="+minDispalcementTileYAll+", maxDispalcementTileYAll="+maxDispalcementTileYAll);
				for (int i=0;i<numSensors;i++){
					System.out.println("calculateDisparityTiles(), minDispalcementTileX["+i+"]="+minDispalcementTileX[i]+", maxDispalcementTileX["+i+"]="+maxDispalcementTileX[i]);
					System.out.println("calculateDisparityTiles(), minDispalcementTileX["+i+"]="+minDispalcementTileX[i]+", maxDispalcementTileX["+i+"]="+maxDispalcementTileX[i]);
					System.out.println("calculateDisparityTiles(), minDispalcementTileY["+i+"]="+minDispalcementTileY[i]+", maxDispalcementTileY["+i+"]="+maxDispalcementTileY[i]);
					System.out.println("calculateDisparityTiles(), minDispalcementTileY["+i+"]="+minDispalcementTileY[i]+", maxDispalcementTileY["+i+"]="+maxDispalcementTileY[i]);
				}
			}

			int [] nextCachedY=new int [numSensors];
			for (int nImg=0;nImg<numSensors;nImg++) nextCachedY[nImg]=0;
			// Window geneartion  for correlation tiles
			double[] window1d=(new DoubleFHT()).getHamming1d(corrFFTSize,0.0); // pure elevated cosine
			final double [] window=new double [corrFFTSize*corrFFTSize];
    		for (int iy=0;iy<corrFFTSize;iy++) for (int ix=0;ix<corrFFTSize;ix++) window[iy*corrFFTSize+ix]=window1d[iy]*window1d[ix];
//  		Window to interpoalate subpixel resolution near correlation local maximums (original pixels)
			double[] windowInterpolation1d=(new DoubleFHT()).getHamming1d(interpolationSize); // Hamming
			final double [] windowInterpolation=new double [interpolationSize*interpolationSize];
    		for (int iy=0;iy<interpolationSize;iy++) for (int ix=0;ix<interpolationSize;ix++) windowInterpolation[iy*interpolationSize+ix]=windowInterpolation1d[iy]*windowInterpolation1d[ix];
//  		Window to interpoalate subpixel resolution near correlation local maximums (subdivided pixels) - divide by this window to cancel the previous one
    		int upSampleSize=interpolationSize*interpolationUpSample;
			double[] windowUpSample1d=(new DoubleFHT()).getHamming1d(upSampleSize); //Hamming
			final double [] iWindowUpSample=new double [upSampleSize*upSampleSize];  // inverse window
    		for (int iy=0;iy<upSampleSize;iy++) for (int ix=0;ix<upSampleSize;ix++) iWindowUpSample[iy*upSampleSize+ix]=1.0/(windowUpSample1d[iy]*windowUpSample1d[ix]);

    		//Create a flat-top window to switch before pixel and subpixel resolution correlation sections
    		double[] windowSample1d0=(new DoubleFHT()).getHamming1d(upSampleSize/2,0.0); //Pure cosine
    		double[] windowSample1d=new double [upSampleSize];
    		for (int i=0;i<upSampleSize;i++){
    			if      (i<upSampleSize/4)   windowSample1d[i]=windowSample1d0[i];
    			else if (i>3*upSampleSize/4) windowSample1d[i]=windowSample1d0[i-upSampleSize/2];
    			else windowSample1d[i]=1.0;
    		}
			final double [] windowSample=new double [upSampleSize*upSampleSize];  // inverse window
    		for (int iy=0;iy<upSampleSize;iy++) for (int ix=0;ix<upSampleSize;ix++) windowSample[iy*upSampleSize+ix]=windowSample1d[iy]*windowSample1d[ix];


    		final String [] channelName={"Y","Cb","Cr","Ext"};
    		int [] channelList=null;
    		if (externalData==null){
    			int numActiveChannels=0;
    			for (int i=0;i<3;i++)if ((channelMask& (1<<i))!=0) numActiveChannels++;
    			channelList=new int [numActiveChannels];
    			numActiveChannels=0;
    			for (int i=0;i<3;i++)if ((channelMask& (1<<i))!=0) channelList[numActiveChannels++]=i;
//    			String externalTitle,
//    			double [][] externalData, // if not null - process instead of channel data
    		} else {
    			channelList=new int [1];
    			channelList[0]=3;
    		}
    		int channelMaskMod=(externalData==null)?channelMask:8;
    		double [] maskedWeights={
    				((channelMaskMod & 1)!=0)?1.0:0.0,
    				((channelMaskMod & 2)!=0)?corrCbWeight:0.0,
    				((channelMaskMod & 4)!=0)?corrCrWeight:0.0,
    				((channelMaskMod & 8)!=0)?1.0:0.0};
    		final double [] componentWeights={
    				maskedWeights[0]/(maskedWeights[0]+maskedWeights[1]+maskedWeights[2]+maskedWeights[3]),
    				maskedWeights[1]/(maskedWeights[0]+maskedWeights[1]+maskedWeights[2]+maskedWeights[3]),
    				maskedWeights[2]/(maskedWeights[0]+maskedWeights[1]+maskedWeights[2]+maskedWeights[3]),
    				maskedWeights[3]/(maskedWeights[0]+maskedWeights[1]+maskedWeights[2]+maskedWeights[3])};

    		final Thread[] threads = newThreadArray(threadsMax);
    		int numPairs=0;
    		for (int i=0;i<otherImage.length;i++) numPairs+=otherImage[i].length;
    		final int [][] firstSecondIndex=new int [numPairs][2]; // first image numer, index  of the second image
    		int index=0;
    		for (int i=0;i<otherImage.length;i++) for (int j=0;j<otherImage[i].length;j++){
    			firstSecondIndex[index  ][0]=i;
    			firstSecondIndex[index++][1]=j;
    		}
			// main loop
//    		for (int chnYCbCr=0;chnYCbCr<3;chnYCbCr++) if ((channelMask& (1<<chnYCbCr))!=0){
        	for (int iChnYCbCr=0;iChnYCbCr<channelList.length;iChnYCbCr++){
    			for (int nImg=0;nImg<numSensors;nImg++) nextCachedY[nImg]=0; // reset before next channel run?
        		for (int tileY=tileYMin;tileY<=tileYMax;tileY++) fhtCache[tileY]=null;



    			final int finalChnYCbCr=channelList[iChnYCbCr];
    			if (updateStatus) IJ.showStatus("Correlating channel "+channelName[finalChnYCbCr]);
    			for (int iTileY=tileYMin;iTileY<=tileYMax;iTileY++){
    				final int tileY=iTileY;
    				if (debugLevel>2){
    					System.out.println("==== calculateDisparityTiles(), correlating channel "+channelName[finalChnYCbCr]+" tileY="+tileY+" ====");
    				}
    				// free memory - delete the full row of unneeded FHT tiles (all images)
    				if ((tileY+minDispalcementTileYAll)>0) {
        				if (debugLevel>2){
        					System.out.println("calculateDisparityTiles(), freeing tile row "+(tileY+minDispalcementTileYAll-1)+" for all images");
        				}
    					fhtCache[tileY+minDispalcementTileYAll-1]=null;
    					for (int nImg=0;nImg<numSensors;nImg++){
    						if (nextCachedY[nImg]<tileY+minDispalcementTileYAll) nextCachedY[nImg]= tileY+minDispalcementTileYAll; // for filling up cache initially only
    					}
    				}
    				// now free more - tiles for particular images that are "obsolete" earlier than others
    				for (int nImg=0;nImg<numSensors;nImg++){
    					for (int tileYDel=tileY+minDispalcementTileYAll;tileYDel<(tileY+minDispalcementTileY[nImg]);tileYDel++) if (tileYDel>=0){
    	    				if (debugLevel>2){
    	    					System.out.println("calculateDisparityTiles(), freeing tile row "+tileYDel+" for image "+nImg);
    	    				}
    						if (fhtCache[tileYDel]!=null) fhtCache[tileYDel][nImg]=null;
    						if (nextCachedY[nImg]<=tileYDel) nextCachedY[nImg]= tileYDel+1; // for filling up cache initially only
    					}
    				}
    				// calculate rows of FHT tiles
    				for (int nImg=0;nImg<numSensors;nImg++){
        				int correctedTileXMin=tileXMin+minDispalcementTileX[nImg];
        				int correctedTileXMax=tileXMax+maxDispalcementTileX[nImg];
        				if (correctedTileXMin<0) correctedTileXMin=0;
        				if (correctedTileXMax>=numTilesX)correctedTileXMax=numTilesX-1;

    					for (;nextCachedY[nImg]<=tileY+maxDispalcementTileY[nImg];nextCachedY[nImg]++) if (nextCachedY[nImg]<numTilesY){
    						getFHTRow(
    								(externalData==null)?null:externalData[nImg], // if not null - use instead of the AYCbCr stack
    								nextCachedY[nImg], // final int tileY, // tile to add
    								finalChnYCbCr, // channel to get
    								nImg,   // number of image to use
    								window, // 2-d window to use 0.5*(cos+1)
    								fhtCache,
    								corrFFTSize,
    								overlapStep,
    								numTilesX,
    								correctedTileXMin,
    								correctedTileXMax,
    								threadsMax,
    								debugLevel);

    					}
    				}
    				// multithreaded calculation of disparity array using fhtCache
    		   		final AtomicInteger rsltTileAtomic     = new AtomicInteger(0);
    		   		for (int ithread = 0; ithread < threads.length; ithread++) {
    		   			threads[ithread] = new Thread() {
    		   				@Override
							public void run() {
    		   					DoubleFHT doubleFHT=new DoubleFHT();
    		   					DoubleFHT interpolationFHT=null; // will be initialized if needed
    		   					boolean useFilter=(correlationHighPassSigma>=0.0) || (correlationLowPassSigma>=0.0);
    		   					if (useFilter) doubleFHT.createFrequencyFilter(
    		   							window, // just to specify size
    		   							correlationHighPassSigma,
    		   							correlationLowPassSigma);
    		   					double [][][] stage=null; // [pY][pX]{data, weight}
    		   					for (int rsltTile=rsltTileAtomic.getAndIncrement(); rsltTile<firstSecondIndex.length*(tileXMax-tileXMin+1);rsltTile=rsltTileAtomic.getAndIncrement()){
//    		   						boolean debugThisTile=(debugLevel>2) && (finalChnYCbCr==0) && (tileY==tileYMin) &&(rsltTile==0); // make
    		   						boolean debugThisTile=(debugLevel>2) && (finalChnYCbCr==0) && (tileY==tileYMin) &&(rsltTile==1); // make

    		   						int nImg=firstSecondIndex[rsltTile%firstSecondIndex.length][0];
    		   						int secondIndex=firstSecondIndex[rsltTile%firstSecondIndex.length][1];
    		   						int secondImage=otherImage[nImg][secondIndex];
    		   						int [][] thisDisplacementTilesList=displacementTilesList[nImg][secondIndex];
    		   						if ((stage==null) ||
    		   								(stage.length!=   correlationMap[nImg][secondIndex].length) ||
    		   								(stage[0].length!=correlationMap[nImg][secondIndex][0].length) ){
    		   							stage=new double [correlationMap[nImg][secondIndex].length][correlationMap[nImg][secondIndex][0].length][2];
    		   						}
    		   						double [] zero2={0.0,0.0};
    		   						// TODO - zero only needed pixels?
    		   						for (int pY=0;pY<stage.length;pY++) for (int pX=0;pX<stage[0].length;pX++){
    		   							stage[pY][pX]=zero2.clone();
    		   						}
    		   						int tileX=rsltTile/firstSecondIndex.length+tileXMin;
    		   						if (debugThisTile){
    		   							System.out.println("\n\n====Processing tileY="+tileY+" tileX="+tileX+" nImg="+nImg+" seconIndex="+secondIndex+" secondImage="+secondImage+"====");
    		   						}
    		   						if (disparityTiles[tileY][tileX]==null){
    		   							disparityTiles[tileY][tileX]=new double [otherImage.length][][];
    		   							for (int i=0;i<otherImage.length;i++)disparityTiles[tileY][tileX][i]=null;
    		   						}
    		   						if (disparityTiles[tileY][tileX][nImg]==null){
    		   							disparityTiles[tileY][tileX][nImg]=new double [otherImage[nImg].length][];
    		   							for (int i=0;i<otherImage[nImg].length;i++){
    		   								// where can this bug be?
    		   								if ((disparityTiles[tileY][tileX]==null)) {
    		   									System.out.println("disparityTiles["+tileY+"]["+tileX+"]==null");
    		   								} else if ((disparityTiles[tileY][tileX][nImg]==null)) {
    		   									System.out.println("disparityTiles["+tileY+"]["+tileX+"]["+nImg+"]==null");
    		   								}
    		   								disparityTiles[tileY][tileX][nImg][i]=null;      // TODO:java.lang.NullPointerException
    		   							}
    		   						}
    		   						if (disparityTiles[tileY][tileX][nImg][secondIndex]==null){
    		   							disparityTiles[tileY][tileX][nImg][secondIndex]=new double [disparitySteps+1];
    		   							for (int i=0;i<=disparitySteps;i++) disparityTiles[tileY][tileX][nImg][secondIndex][i]=0.0; // TODO:java.lang.NullPointerException- again
    		   						}
    		   						// ran with size =64, overlap=1/4 Exception in thread "Thread-8049" java.lang.NullPointerException   at PixelMapping$InterSensor$4.run(PixelMapping.java:4557)
    		   						// now calculate and accumulate correlation to stage
    		   						// sanity check
    		   						if ((tileY>=0) && (tileY<numTilesY) && (tileX>=0) && (tileX<numTilesX) &&
    		   								(fhtCache[tileY]!=null) && (fhtCache[tileY]!=null) && (fhtCache[tileY][secondImage]!=null) && (fhtCache[tileY][secondImage][tileX]!=null)){
    		   						} else {
    		   							System.out.println("reference correlation tile does not exist, tileX="+tileX+", tileY="+tileY+" channlel "+channelName[finalChnYCbCr]+" image "+nImg);
    		   							continue;
    		   						}
    		   						for (int tileNum=0;tileNum<thisDisplacementTilesList.length;tileNum++){
    		   							int tX=tileX+thisDisplacementTilesList[tileNum][0];
    		   							int tY=tileY+thisDisplacementTilesList[tileNum][1];
    		   							if ((tY>=0) && (tY<numTilesY) && (tX>=0) && (tX<numTilesX) &&
    		   									(fhtCache[tY]!=null) && (fhtCache[tY]!=null) && (fhtCache[tY][secondImage]!=null) && (fhtCache[tY][secondImage][tX]!=null)){
    		   								if (debugThisTile) System.out.println("Correlating with tile tY="+tY+" tX="+tX+" defined, proceeding...");
    		   								double [] second=fhtCache[tY][secondImage][tX].clone();
    		   								double [] first= fhtCache[tileY][nImg][tileX].clone(); // does it need to be cloned
    		   								double [] result=second; //alias
    		   								if (corrPhaseFraction<=0) doubleFHT.multiply(result,first,  true); // correlation, not convolution. Reverse order!
    		   								else result=doubleFHT.phaseMultiply(result,first,  corrPhaseFraction); // correlation, not convolution. Reverse order!
    		   								if (useFilter)  result=doubleFHT.applyFreqFilter(result,null); // use earlier calculated filter
    		   								doubleFHT.transform(result,true) ; // inverse transform
    		   								doubleFHT.swapQuadrants(result);
    		   								// accumulate result to stage area (window - also)
    		   								int topLeftY=(stage.length-corrFFTSize)/2 +overlapStep*thisDisplacementTilesList[tileNum][1];   // top left corner of the reference tile
    		   								int topLeftX=(stage[0].length-corrFFTSize)/2 +overlapStep*thisDisplacementTilesList[tileNum][0];
    		   								if (debugThisTile){
    		   									System.out.println("thisDisplacementTilesList["+tileNum+"][0]="+thisDisplacementTilesList[tileNum][0]+
    		   											" thisDisplacementTilesList["+tileNum+"][1]="+thisDisplacementTilesList[tileNum][1]);
    		   									System.out.println("stage.length="+stage.length+" stage[0].length="+stage[0].length);
    		   									System.out.println("corrFFTSize="+corrFFTSize+" overlapStep="+overlapStep);
    		   									System.out.println("topLeftY="+topLeftY+" topLeftX="+topLeftX);
    		   								}
    		   								if (debugThisTile && (debugLevel>3)) {
    		   									double [] second_dbg=fhtCache[tY][secondImage][tX].clone();
    		   									doubleFHT.transform(second_dbg,true) ; // inverse transform
    		   									doubleFHT.swapQuadrants(second_dbg);

    		   									double [] first_dbg= fhtCache[tileY][nImg][tileX].clone(); // does it need to be cloned
    		   									doubleFHT.transform(first_dbg,true) ; // inverse transform
    		   									doubleFHT.swapQuadrants(first_dbg);
    		   									double [][] result_dbg={first_dbg,second_dbg,result.clone()};
    		   									String [] titles_dbg={"first","second","corr"};

    		   									(new ShowDoubleFloatArrays()).showArrays(
    		   											result_dbg,
    		   											corrFFTSize,
    		   											corrFFTSize,
    		   											true,
    		   											"R"+tY+"-"+tX+"_"+tileY+"-"+tileX,
    		   											titles_dbg);
    		   								}
    		   								for (int index=0;index<result.length;index++){
    		   									int stageY=index/corrFFTSize+topLeftY;
    		   									int stageX=index%corrFFTSize+topLeftX;
    		   									if ((stageY>=0) && (stageY<stage.length) && (stageX>=0) && (stageX<stage[0].length)){
    		   										stage[stageY][stageX][0]+=result[index];
    		   										stage[stageY][stageX][1]+=window[index];
    		   										if (debugThisTile && (debugLevel>4)) {
    		   											System.out.println(tX+" "+tY+" "+index+" "+stageY+" "+stageX+" "+result[index]+" "+window[index]);
    		   										}
    		   									}
    		   								}
    		   							} else {
    		   								if (debugThisTile) System.out.println("!!!!! No FHT data for tile tY="+tY+" tX="+tX);
    		   								// possible debug here - matching tile FHT does not exist
    		   							}
    		   						}
    		   						if (debugThisTile && (debugLevel>4)){
    		   							for (int dbg_i=0;dbg_i<stage.length;dbg_i++){
    		   								System.out.print("\n==== "+dbg_i);

    		   								for (int dbg_j=0;dbg_j<stage[0].length;dbg_j++){
    		   									System.out.print(" "+dbg_j+":"+stage[dbg_i][dbg_j][0]+"/"+stage[dbg_i][dbg_j][1]);
    		   								}
    		   							}
    		   							System.out.println();
    		   						}
    		   						// debug center tile only
    		   						if (debugThisTile){
    		   							double [][] debugStage=new double [3][stage.length*stage[0].length];
    		   							for (int i=0;i<debugStage[0].length;i++){
    		   								int dbgY=i/stage[0].length;
    		   								int dbgX=i%stage[0].length;
    		   								debugStage[0][i]=(stage[dbgY][dbgX][1]>0)?(stage[dbgY][dbgX][0]/stage[dbgY][dbgX][1]):0.0;
    		   								debugStage[1][i]=stage[dbgY][dbgX][0];
    		   								debugStage[2][i]=stage[dbgY][dbgX][1];

    		   							}
    		   							String [] titles={"result","data","weight"};
    		   							(new ShowDoubleFloatArrays()).showArrays(
    		   									debugStage,
    		   									stage[0].length,
    		   									stage.length,
    		   									true,
    		   									"C"+nImg+"-"+secondImage+"_tY"+tileY+"_tX"+tileX+"-"+channelName[finalChnYCbCr],
    		   									titles);
    		   						}

    		   						// result accumulated data will need to be divided by weights. Maybe later do only on some?
    		   						double []correlationSection=new double [disparitySteps+1];
    		   						for (int pY=0;pY<stage.length;pY++) for (int pX=0;pX<stage[0].length;pX++){
    		   							if (stage[pY][pX][1]!=0) stage[pY][pX][0]/=stage[pY][pX][1];
    		   						}

    		   						// traverse disparity and linear interpolate (first pass)
    		   						double maxDX=disparityRange*(disparityScales[secondImage][0]-disparityScales[nImg][0]);
    		   						double maxDY=disparityRange*(disparityScales[secondImage][1]-disparityScales[nImg][1]);
    		   						double maxD=Math.sqrt(maxDX*maxDX+maxDY*maxDY);
    		   						int centerX=stage[0].length/2;
    		   						int centerY=stage.length/2;
    		   						for (int di=0;di<=disparitySteps;di++){
    		   							double dX=(maxDX*di)/disparitySteps+centerX;
    		   							double dY=(maxDY*di)/disparitySteps+centerY;
    		   							int iX=(int) Math.floor(dX);
    		   							int iY=(int) Math.floor(dY);
    		   							dX-=iX;
    		   							dY-=iY;
    		   							correlationSection[di]= // bi-linear interpolation
    		   								(stage[iY  ][iX][0]*(1.0-dX) + stage[iY  ][iX+1][0]*dX)*(1.0-dY)+
    		   								(stage[iY+1][iX][0]*(1.0-dX) + stage[iY+1][iX+1][0]*dX)*     dY;
    		   						}
    		   						if (debugThisTile){
    		   							System.out.println("\nPixel resolution bilinear disparity array for nImg="+nImg+" secondImage="+secondImage+" channel:"+channelName[finalChnYCbCr]);
    		   							for (int di=0;di<=disparitySteps;di++){
    		   								System.out.println(di+" "+nImg+" "+secondImage+" "+correlationSection[di]);
    		   							}
    		   						}
    		   						// find local max
    		   						double threshold=subpixAMax[finalChnYCbCr];
    		   						boolean [] isMax =new boolean[correlationSection.length];
    		   						int iMax=-1;
    		   						int numMax=0;
    		   						for (int i=0;i<correlationSection.length;i++){
    		   							if ((correlationSection[i]>=threshold) &&
    		   									((i==0)||(correlationSection[i]>=correlationSection[i-1])) &&
    		   									((i==(correlationSection.length-1)) || (correlationSection[i]>=correlationSection[i+1]))){
    		   								isMax[i]=true; // local max above absolute threshold
    		   								if ((numMax==0) || (correlationSection[i]>correlationSection[iMax])) iMax=i;
    		   								numMax=1;
    		   							} else {
    		   								isMax[i]=false;
    		   							}
    		   						}
    		   						double mergeMax=1.5; //merge maximums for interpolation if closer than mergeMax;
		   							int iMergeMax=(int) Math.round(mergeMax*disparitySteps/maxD);
    		   						if (debugThisTile){
    		   							System.out.println("mergeMax="+mergeMax+" iMergeMax="+iMergeMax);
    		   						}
    		   						if (numMax>0){
    		   							numMax=0;
    		   							threshold=correlationSection[iMax]*subpixRMax;
    		   							for (int i=0;i<correlationSection.length;i++) if (isMax[i]){
    		   								if (correlationSection[i]>=threshold) numMax++;
    		   								else isMax[i]=false;
    		   							}

    		   							if (numMax>subpixNMax) numMax=subpixNMax; // limit to specified number of correlation maximums to subpixel
    		   							int [] maxIndices=new int [numMax];
    		   							maxIndices[0]=iMax;
    		   							isMax[iMax]=false;
    		   							int maxNum;
    		   							for (maxNum=1;maxNum<maxIndices.length;maxNum++){
    		   								// merge previous one if possible
    		   								int nearDown=-1, nearUp=-1;
    		   								for (int i=0;i<=iMergeMax;i++){
    		   									if ((maxIndices[maxNum-1]-i)<0) break;
    		   									if (isMax[maxIndices[maxNum-1]-i]){
    		   										nearDown = maxIndices[maxNum-1]-i;
    		   										break;
    		   									}
    		   								}
    		   								for (int i=0;i<=iMergeMax;i++){
    		   									if ((maxIndices[maxNum-1]+i)>disparitySteps) break;
    		   									if (isMax[maxIndices[maxNum-1]+i]){
    		   										nearUp = maxIndices[maxNum-1]+i;
    		   										break;
    		   									}
    		   								}
    		   								int n=1+((nearDown>=0)?1:0)+((nearUp>=0)?1:0);
    		   								if (n>1){
    		   									double s=maxIndices[maxNum-1]+((nearDown>=0)?nearDown:0)+((nearUp>=0)?nearUp:0);
    		   									int iMerged=(int) Math.round(s/n);

    	        		   						if (debugThisTile){
    	        		   							System.out.println("Merging close maximums: "+maxIndices[maxNum-1]+" with "+
    	        		   									((nearDown>=0)?nearDown:"")+" "+((nearUp>=0)?nearUp:"")+" to "+iMerged);
    	        		   						}
    	        		   						maxIndices[maxNum-1]=iMerged;
    	        		   						isMax[iMerged]=false;
    	        		   						if (nearDown>=0) isMax[nearDown]=false;
    	        		   						if (nearUp>=0)   isMax[nearUp]=false;
    		   								}
    		   								iMax=-1;
    		   								for (int i=0;i<correlationSection.length;i++) {
    		   									if (isMax[i]&&((iMax<0) || (correlationSection[i]>=correlationSection[iMax]))) iMax=i;
    		   								}
    		   								if (iMax<0) break; // no more maximums
    		   								isMax[iMax]=false;
    		   								maxIndices[maxNum]=iMax;
    		   							}
    		   							//TODO: disable max closer than 1.5pix - they will be covered by previous interpolation anyway
        		   						if (debugThisTile){
        		   							System.out.println("List maximums on correlation section");
        		   							double d=Math.sqrt(maxDX*maxDX+maxDY*maxDY)/disparitySteps;
        		   							for (int n=0;n<maxNum;n++) {
        		   								System.out.println(n+" "+maxIndices[n]+"("+(d*maxIndices[n])+" pix) "+correlationSection[maxIndices[n]]);
        		   							}
        		   						}

    		   							// iterate through all maximums, subdivide pixels around them and replace data in correlationSection
    		   							// alternative way to get subpixel - correlate with proper centered tile (not on a grid), upsample result
    		   							for (int n=0;n<maxNum;n++){
    		   								double dXM=(maxDX*maxIndices[n])/disparitySteps; //+centerX;
    		   								double dYM=(maxDY*maxIndices[n])/disparitySteps; //+centerY;
    		   								int iXC=(int) Math.round(dXM);
    		   								int iYC=(int) Math.round(dYM);
    		   								int iXTL=centerX+iXC-interpolationSize/2; // interpolation top lext corner relative to TL of the stage array
    		   								int iYTL=centerY+iYC-interpolationSize/2;
    		   								double [] interpolationInput=new double [interpolationSize*interpolationSize];

            		   						if (debugThisTile){
            		   							System.out.println("iXC="+iXC+" iYC="+iYC+" iXTL="+iXTL+" iYTL="+iYTL);
            		   						}


    		   								for (int index=0;index<interpolationInput.length;index++){
    		   									int iiY=index/interpolationSize+iYTL;
    		   									int iiX=index%interpolationSize+iXTL; // no check for OOB - should be guaranteed
                		   						if (debugThisTile && (debugLevel>3)){
                		   							System.out.print("iiY="+iiY+" iiX="+iiX);
                		   							System.out.println(" stage["+iiY+"]["+iiX+"][0]="+stage[iiY][iiX][0]);
                		   						}
    		   									interpolationInput[index]=stage[iiY][iiX][0]; // data, already divided by weight TODO: OOB 272,293,288,325, ...
    		   								}
            		   						if (debugThisTile && (debugLevel>3)){
            		   							(new ShowDoubleFloatArrays()).showArrays(
            		   									interpolationInput,
            		   									interpolationSize,
            		   									interpolationSize,
            		   									"II"+n+"-"+iXTL+"-"+iYTL);

            		   						}
    		   								// normalize and save DC
    		   								double dc=normalizeAndWindowGetDC (interpolationInput, windowInterpolation);
            		   						if (debugThisTile && (debugLevel>3)){
            		   							(new ShowDoubleFloatArrays()).showArrays(
            		   									interpolationInput,
            		   									interpolationSize,
            		   									interpolationSize,
            		   									"IW"+n+"-"+iXTL+"-"+iYTL);

            		   						}

    		   								if (interpolationFHT==null) interpolationFHT=new DoubleFHT();
    		   								double [] upsampled=interpolationFHT.upsample(interpolationInput,interpolationUpSample); // add scalind to upsample()
            		   						if (debugThisTile && (debugLevel>3)){
            		   							(new ShowDoubleFloatArrays()).showArrays(
            		   									upsampled,
            		   									interpolationSize*interpolationUpSample,
            		   									interpolationSize*interpolationUpSample,
            		   									"RO"+n+"-"+iXTL+"-"+iYTL);

            		   						}

    		   								// divide by scaled window, add removed earlier DC level
    		   								for (int i=0;i<upsampled.length;i++){
    		   									upsampled[i]=upsampled[i]*iWindowUpSample[i]+dc;
    		   								}
            		   						if (debugThisTile){
            		   							(new ShowDoubleFloatArrays()).showArrays(
            		   									upsampled,
            		   									interpolationSize*interpolationUpSample,
            		   									interpolationSize*interpolationUpSample,
            		   									"IO"+n+"-"+iXTL+"-"+iYTL);

            		   						}
    		   								//windowSample
    		   								iXTL-=centerX; // relative to reference center
    		   								iYTL-=centerY; // relative to reference center
    		   								int iXBR=iXTL+interpolationSize; // interpolation bottom right corner relative to TL of the stage array
    		   								int iYBR=iYTL+interpolationSize;
    		   								int interpolatedSize=interpolationSize*interpolationUpSample;
            		   						if (debugThisTile){
            		   							System.out.println(" iXTL="+iXTL+" iYTL="+iYTL+" iXBR="+iXBR+" iYBR="+iYBR);
            		   						}

    		   								for (int di=0;di<=disparitySteps;di++){
    		   									double dX=(maxDX*di)/disparitySteps; //+centerX;
    		   									double dY=(maxDY*di)/disparitySteps; //+centerY;
    		   									int iX=(int) Math.floor(dX);
    		   									int iY=(int) Math.floor(dY);
    		   									if ((iY>=iYTL) && (iX>=iXTL) && (iY<iYBR) && (iX<iXBR)){
    		   										//dX-=iX-iXTL;
    		   										//dY-=iY-iYTL;
    		   										dX-=iXTL; // now relative to interpolation square
    		   										dY-=iYTL;
    		   										dX*=interpolationUpSample;
    		   										dY*=interpolationUpSample;
    		   										int isX=(int) Math.floor(dX);
    		   										int isY=(int) Math.floor(dY);
    		   										if ((isX>=(interpolatedSize-1)) || (isY>=(interpolatedSize-1))) continue; // no room for bilinear
    		   										dX-=isX;
    		   										dY-=isY;
    		   										//	        		   								bilinear on subpixels
    		   										int index00=isY*interpolatedSize+isX;
    		   										int index01=index00+interpolatedSize;
    		   										int index10=index00+1;
    		   										int index11=index01+1;
    	            		   						if (debugThisTile){
    	            		   							System.out.print(" Updating point="+di+" isX=" +isX+" isY=" +isY+
    	            		   									" old="+IJ.d2s(correlationSection[di],3));
    	            		   						}

    		   										double d= // bi-linear interpolated data
    		   											(upsampled   [index00]*(1.0-dX)+upsampled   [index10]*dX)*(1.0-dY)+
    		   											(upsampled   [index01]*(1.0-dX)+upsampled   [index11]*dX)*     dY;
    		   										double w= // bi-linear interpolated window
    		   											(windowSample[index00]*(1.0-dX)+windowSample[index10]*dX)*(1.0-dY)+
    		   											(windowSample[index01]*(1.0-dX)+windowSample[index11]*dX)*     dY;
    		   										correlationSection[di]=(1.0-w)*correlationSection[di]+w*d;
    	            		   						if (debugThisTile){
    	            		   							System.out.println(" new="+IJ.d2s(d,3)+" weight="+IJ.d2s(w,2)+" updated="+IJ.d2s(correlationSection[di],3));
    	            		   						}

    		   									}
    		   								}
    	    		   						if (debugThisTile){
    	    		   							System.out.println("\nUpdated subpixel resolution bilinear disparity array for nImg="+nImg+" secondImage="+secondImage+" channel:"+channelName[finalChnYCbCr]);
    	    		   							for (int di=0;di<=disparitySteps;di++){
    	    		   								System.out.println(di+" "+nImg+" "+secondImage+" "+correlationSection[di]);
    	    		   							}
    	    		   						}

    		   							} // for (int n=1;n<maxNum;n++){
    		   						} // if (numMax>0)
    		   						// now we have correlationSection[], need to combine multiple colors
    		   						for (int di=0;di<correlationSection.length;di++){
    		   							disparityTiles[tileY][tileX][nImg][secondIndex][di]+=componentWeights[finalChnYCbCr]*correlationSection[di];
    		   						}
    		   					}
    		   				}
    		   			};
    		   		}
    		   		startAndJoin(threads);
    		   		IJ.showProgress(iTileY-tileYMin+1,tileYMax-tileYMin);
    			}
    			IJ.showProgress(1.0);
    		}
    		return disparityTiles;
       	}


    	/**
    	 * Calculate one row of FHT tiles for 1 image/one YCbCr channel, save it in the cache for correlation calculation
    	 * @param altImage     alternative image array - use instead of the YCbCr slice
    	 * @param tileY        number of tile row to calculate
    	 * @param chnYCbCr     channel number 0 - Y, 1 - Cb, 2 - Cr
    	 * @param numImg       image number to use
    	 * @param window       elevated 2-d cosine (0.5*cos+1) window
    	 * @param fhtCache     [tilesY][numImg][tilesx][] array that will be filled with the direct [corrFFTSize*corrFFTSize] FHT data
    	 * @param corrFFTSize  size of the tile side, pixels
    	 * @param overlapStep  tile period, pixels (fraction of corrFFTSize)
    	 * @param numTilesX    number of hoprizontal tiles in a row
    	 * @param tileXMin     number of the first (lowest number) of the tile in a row to process
    	 * @param tileXMax     number of the last (highest number) of the tile in a row to process
    	 * @param threadsMax   numers of therads to use
    	 * @param debugLevel   debug level
    	 */

     	private void getFHTRow(
     			final double [] altImage, // alternative image array
				final int tileY, // tile to add
				final int chnYCbCr, // channel to get
				final int numImg,   // number of image to use
				final double [] window, // 2-d window to use 0.5*(cos+1)
				final double [][][][] fhtCache,
				final int corrFFTSize,
				final int overlapStep,
				final int numTilesX,
				final int tileXMin,
				final int tileXMax,
				final int threadsMax,
				final int debugLevel){
     		if (debugLevel>2){
     			System.out.println("getFHTRow(), tileY="+tileY);
     		}
    		    int numSensors=this.channel.length;
//    			if (fhtCache[tileY]==null) fhtCache[tileY]=new double [numSensors][numTilesX][];
    			if (fhtCache[tileY]==null) {
    				fhtCache[tileY]=new double [numSensors][][];
    				for (int i=0;i<fhtCache[tileY].length;i++) fhtCache[tileY][i]=null;
    			}
    			if (fhtCache[tileY][numImg]!=null){
    				System.out.println("getFHTRow() BUG: row "+tileY+" image "+numImg+" is already not null!!!!!");
    			}
    			fhtCache[tileY][numImg]=new double [numTilesX][];
    			for (int tileX=0;tileX<numTilesX;tileX++) fhtCache[tileY][numImg][tileX]=null; // null pointer
        		int numLayers=this.overlapImages.length/numSensors;
        		final double [] imageSliceAlpha=this.overlapImages[numImg*numLayers];
        		final double [] imageSlice=(altImage==null)?this.overlapImages[chnYCbCr+numImg*numLayers+1]:altImage;
        		final int length=corrFFTSize*corrFFTSize;
           		final Thread[] threads = newThreadArray(threadsMax);
           		final AtomicInteger tileXAtomic     = new AtomicInteger(tileXMin);
           		for (int ithread = 0; ithread < threads.length; ithread++) {
           			threads[ithread] = new Thread() {
           				@Override
						public void run() {
           	        		double [] alphaSlice=new double[length];
           	        		double [] channelSlice=new double[length];
           	        		DoubleFHT doubleFHT=new DoubleFHT();
           					for (int tileX=tileXAtomic.getAndIncrement(); tileX<=tileXMax;tileX=tileXAtomic.getAndIncrement()){
           			     		if (debugLevel>2){
           			     			System.out.println("getFHTRow(), selecting image "+numImg+" tileX="+tileX+" xc="+(overlapStep*tileX)+" yc="+(overlapStep*tileY));
           			     		}
           	        			getSelection(
           	        					imageSliceAlpha, // one image/channel slice
           	        					alphaSlice,
           	        					corrFFTSize, //int width,
           	        					corrFFTSize, //int height,
           	        					overlapStep*tileX , //xc,
           	        					overlapStep*tileY); //yc);
           	        			for (int i=0;i<length;i++) 	alphaSlice[i]*=window[i];
           	        			getSelection(
           	        					imageSlice, // one image/channel slice
           	        					channelSlice,
           	        					corrFFTSize, //int width,
           	        					corrFFTSize, //int height,
           	        					overlapStep*tileX , //xc,
           	        					overlapStep*tileY); //yc);
           	    				normalizeAndWindow(channelSlice, alphaSlice,true);
           	    				doubleFHT.swapQuadrants(channelSlice);
           	    				doubleFHT.transform(channelSlice,false); // direct FHT transform TODO: use frequency filter
           	    				fhtCache[tileY][numImg][tileX]=channelSlice.clone();
           					}
           				}
           			};
           		}
           		startAndJoin(threads);
    	}




    	public double [][] correlate(
    			int [][] iCenterXY, // for each image - centerX, centerY
    			boolean autoCorrelation,
				int corrFFTSize,
				double phaseCorrelationFraction,
				double corrCbWeight,
				double corrCrWeight,
				double correlationHighPassSigma,
				double correlationLowPassSigma,
				double noiseNormalizationSignaY,
				double noiseNormalizationSignaCbCr,
				double contrastThreshold,
				boolean useBinaryAlpha,
				int threadsMax,
				int debugLevel){

    		int length=corrFFTSize*corrFFTSize;
    		int numSensors=this.channel.length;
    		int numLayers=this.overlapImages.length/numSensors;
    		double [][] selection = getSelection(
    				corrFFTSize,
    				corrFFTSize,
    				iCenterXY);
    		/*
    		double [] stats= selectionStats(selection,useBinaryAlpha);
    		if (debugLevel>0){
    			System.out.println("Selection ("+(useBinaryAlpha?"binary":"analog")+" alpha ) "+
    					" iCenterXY[0][0]="+iCenterXY[0][0]+" iCenterXY[0][1]="+iCenterXY[0][1]+" size="+corrFFTSize+
    					" : first "+IJ.d2s(100*stats[0],3)+"%, second "+IJ.d2s(100*stats[1],3)+"%, overlap "+IJ.d2s(100*stats[2],3)+"%");
    		}
    		*/

    		if (debugLevel>1){
    				String [] titles0={"Alpha","Y","Cb","Cr"};
    				String [] titles=new String[this.overlapImages.length];
    				for (int i=0;i<titles.length;i++){
    					titles[i]=titles0[i%titles0.length]+(i/titles0.length);
    				}
    			(new ShowDoubleFloatArrays()).showArrays( // skip if any are zeros?
    					selection,
    					corrFFTSize,
    					corrFFTSize,
    					true,
//    			"selection-x"+corrXC+"-y"+corrYC+"-SHFT"+corrShift+"-PC"+phaseCorrelationFraction,
    			"selection-x"+iCenterXY[0][0]+"-y"+iCenterXY[0][1]+"-PC"+phaseCorrelationFraction,
    			titles);
    		}
// correlation stuff
    		double[] hamming1d=(new DoubleFHT()).getHamming1d(corrFFTSize);
    		double [][] window=new double [numSensors][];
    		for (int i=0;i<numSensors;i++) window[i]=selection[numLayers*i].clone();
    		for (int iy=0;iy<corrFFTSize;iy++) for (int ix=0;ix<corrFFTSize;ix++){
    			int index=iy*corrFFTSize+ix;
    			double h=hamming1d[iy]*hamming1d[ix];
    			for (int i=0;i<numSensors;i++) window[i][index]*=h;
    		}
    		double [] componentWeights={
    				1.0/(1.0+corrCbWeight+corrCrWeight),
    				corrCbWeight/(1.0+corrCbWeight+corrCrWeight),
    				corrCrWeight/(1.0+corrCbWeight+corrCrWeight)};
// create and iterate through channel pairs
    		int numValidImages=0;
    		for (int i=0;i<numSensors;i++) if (this.overlapImages[i*numLayers]!=null) numValidImages++;
    		int [][] pairs;
    		if (autoCorrelation) {
    			pairs=new int [numValidImages][2];
    			int index=0;
    			for (int i=0;i<numSensors;i++) if (this.overlapImages[i*numLayers]!=null){
    				pairs[index][0]	=i;
    				pairs[index][1]	=i;
    				index++;
    			}
    		} else {
    			pairs=new int [(numValidImages*(numValidImages-1))/2][2];
    			int index=0;
    			for (int i=0;i<numSensors;i++) if (this.overlapImages[i*numLayers]!=null){
    				for (int j=i+1;j<numSensors;j++) if (this.overlapImages[i*numLayers]!=null){
    					pairs[index][0]	=i;
    					pairs[index][1]	=j;
    					index++;
    				}
    			}
    		}
    		if (debugLevel>0){
    			for (int nPair=0;nPair<pairs.length;nPair++) if (pairs[nPair]!=null) {
    				System.out.println("pair["+nPair+"]={"+pairs[nPair][0]+","+pairs[nPair][1]+"}");
    			}
    		}
    		double [][] corr=new double[numLayers*pairs.length][length];
    		for (int nPair=0;nPair<pairs.length;nPair++) {
    			int iFirst=pairs[nPair][0];
    			int iSecond=pairs[nPair][1];
        		for (int i=0;i<length;i++) corr[nPair*numLayers][i]=0;
    			for (int l=1;l<numLayers;l++){
    				int layerFirst= numLayers*iFirst+l;
    				int layerSecond=numLayers*iSecond+l;
    				int layerResult0=numLayers*nPair;
    				int layerResult=layerResult0+l;
    				double [] first= selection[layerFirst].clone();
    				double [] second=selection[layerSecond].clone();
    				normalizeAndWindow(first, window[iFirst],true);
    				normalizeAndWindow(second,window[iSecond],true);
    				if ((debugLevel>2) && (l==1)){
    					double [][] firstSecond={first,second};
    					String [] firstSecondTitles={"first","second"};
    					(new ShowDoubleFloatArrays()).showArrays(
    							firstSecond,
    							corrFFTSize,
    							corrFFTSize,
    							true,
//    							"windowed-Y"+iCenterXY[0][0]+"-y"+iCenterXY[0][1]+"-PC"+phaseCorrelationFraction);
    							"windowed-Y"+iFirst+"-:"+iSecond,
    							firstSecondTitles);
    				}
    				// TODO: use common (per-thread) DoubleFHT instance to use caching of sin,cos, filter tables
    				corr[layerResult]=(new DoubleFHT()).correlate (
    						first,
    						second,
    						correlationHighPassSigma,
    						correlationLowPassSigma,
    						phaseCorrelationFraction);
    				double sigma=(l==1)?noiseNormalizationSignaY:noiseNormalizationSignaCbCr;
    				if (sigma>0){
    					double rms=hFNoiseNormalize(
    							corr[nPair*numLayers+l], // double [] data,
    							sigma, // double filterSigma,
    							true);  // boolean centerOnly);
    					if (this.debugLevel>1){
    						System.out.println("Correlation RMS for component "+((l==1)?"Y":((l==2)?"Cb":"Cr"))+ " was "+rms);
    					}
    				}
    				for (int j=0;j<length;j++) corr[layerResult0][j]+=componentWeights[l-1]*corr[layerResult][j];

    			}
    		}

// specific for disparity along horizontal line only
/*
    		double [] corrRes=correlationResults(
        			corr[0],
        			contrastThreshold,
        			corrShift,
        			enableNegativeDisparity,
        			true, // boolean orderMaximums,
        			debugLevel);
    		if (debugLevel>0){
    			if (corrRes==null){
    				System.out.println("Correlation results are below threshold");
    			} else {
    				System.out.println("Contrast threshold="+contrastThreshold);

    				System.out.println("Correlation lower disparity="+corrRes[0]+", higher disparity="+corrRes[1]+
    						" Best disparity="+IJ.d2s(corrRes[2],2)+" (contrast="+IJ.d2s(corrRes[3],3)+")");
    				if (corrRes.length>4){
    					for (int i=2;i<(corrRes.length/2-1);i++){
    						System.out.println("Correlation disparity max # "+i+": "+IJ.d2s(corrRes[2*i],2)+
    								" (contrast="+IJ.d2s(corrRes[2*i+1],3)+")");
    					}
    				}
    			}

    		}
*/
    		return corr;
    	}


    	public double [][] correlate( // old for image pairs
    			boolean autoCorrelation,
				int corrFFTSize,
				int corrXC,
				int corrYC,
				int corrShift,
				double phaseCorrelationFraction,
				double corrCbWeight,
				double corrCrWeight,
				double correlationHighPassSigma,
				double correlationLowPassSigma,
				double noiseNormalizationSignaY,
				double noiseNormalizationSignaCbCr,
				double contrastThreshold,
				boolean useBinaryAlpha,
				boolean enableNegativeDisparity,
				int threadsMax,
				int debugLevel){
    		int length=corrFFTSize*corrFFTSize;
//    		double [][] selection = new double[this.overlapImages.length][length];
    		int numLayers=this.overlapImages.length/2;
    		double [][] selection = getSelection(
    				autoCorrelation?1:0,
    				corrFFTSize,
    				corrXC,
    				corrYC,
    				corrShift); // positive - shift second image, negative - shift first image
    		double [] stats= selectionStats(selection,useBinaryAlpha);
    		if (debugLevel>0){
    			System.out.println("Selection ("+(useBinaryAlpha?"binary":"analog")+" alpha ) xc="+corrXC+" yc="+corrYC+" size="+corrFFTSize+" shift="+corrShift+
    					" : first "+IJ.d2s(100*stats[0],3)+"%, second "+IJ.d2s(100*stats[1],3)+"%, overlap "+IJ.d2s(100*stats[2],3)+"%");
    		}

    		if (debugLevel>1){
    			String [] titles={"Alpha1","Y1","Cb1","Cr1","Alpha2","Y2","Cb2","Cr2"};
    			(new ShowDoubleFloatArrays()).showArrays(
    					selection,
    					corrFFTSize,
    					corrFFTSize,
    					true,
    			"selection-x"+corrXC+"-y"+corrYC+"-SHFT"+corrShift+"-PC"+phaseCorrelationFraction,
    			titles);
    		}
// correlation stuff
    		double[] hamming1d=(new DoubleFHT()).getHamming1d(corrFFTSize);
    		double [][] window={selection[0].clone(),selection[numLayers].clone()};
    		for (int iy=0;iy<corrFFTSize;iy++) for (int ix=0;ix<corrFFTSize;ix++){
    			int index=iy*corrFFTSize+ix;
    			double h=hamming1d[iy]*hamming1d[ix];
    			window[0][index]*=h;
    			window[1][index]*=h;
    		}
    		double [][] corr=new double[numLayers][length];
    		for (int i=0;i<length;i++) corr[0][i]=0;
    		double [] componentWeights={
    				1.0/(1.0+corrCbWeight+corrCrWeight),
    				corrCbWeight/(1.0+corrCbWeight+corrCrWeight),
    				corrCrWeight/(1.0+corrCbWeight+corrCrWeight)};
    		for (int i=0;i<numLayers-1;i++){
    			double [] first= selection[i+1].clone();
    			double [] second=selection[i+1+numLayers].clone();
    			normalizeAndWindow(first, window[0],true);
        		if ((debugLevel>2) && (i==0)){
        			(new ShowDoubleFloatArrays()).showArrays(
        					first,
        					corrFFTSize,
        					corrFFTSize,
        			"windowed-Y"+corrXC+"-y"+corrYC+"-SHFT"+corrShift+"-PC"+phaseCorrelationFraction);
        		}


    			normalizeAndWindow(second,window[1],true);
// TODO: use common (per-thread) DoubleFHT instance to use caching of sin,cos, filter tables
    			corr[i+1]=(new DoubleFHT()).correlate (
    					first,
    					second,
						 correlationHighPassSigma,
						 correlationLowPassSigma,
						 phaseCorrelationFraction);
    			double sigma=(i==0)?noiseNormalizationSignaY:noiseNormalizationSignaCbCr;
    			if (sigma>0){
    				double rms=hFNoiseNormalize(
    						corr[i+1], // double [] data,
    		    			sigma, // double filterSigma,
    		    			true);  // boolean centerOnly);
    				if (this.debugLevel>1){
    					System.out.println("Correlation RMS for component "+((i==0)?"Y":((i==1)?"Cb":"Cr"))+ " was "+rms);
    				}
    			}
//        		for (int j=0;j<length;j++) corr[0][j]+=corr[i][j];
        		for (int j=0;j<length;j++) corr[0][j]+=componentWeights[i]*corr[i+1][j];

    		}

// specific for disparity along horizontal line only

    		double [] corrRes=correlationResults(
        			corr[0],
        			contrastThreshold,
        			corrShift,
        			enableNegativeDisparity,
        			true, // boolean orderMaximums,
        			debugLevel);
    		if (debugLevel>0){
    			if (corrRes==null){
    				System.out.println("Correlation results are below threshold");
    			} else {
    				System.out.println("Contrast threshold="+contrastThreshold);

    				System.out.println("Correlation lower disparity="+corrRes[0]+", higher disparity="+corrRes[1]+
    						" Best disparity="+IJ.d2s(corrRes[2],2)+" (contrast="+IJ.d2s(corrRes[3],3)+")");
    				if (corrRes.length>4){
    					for (int i=2;i<(corrRes.length/2-1);i++){
    						System.out.println("Correlation disparity max # "+i+": "+IJ.d2s(corrRes[2*i],2)+
    								" (contrast="+IJ.d2s(corrRes[2*i+1],3)+")");
    					}
    				}
    			}

    		}

    		return corr;
    	}

    	/**
    	 * Find min/max positions of above threshold values, positions and values of largest maximums above threshold.
    	 * All positions are relative to the center.
    	 * @param corr square correlation array (only horizontal middle line is used)
    	 * @param threshold minimal correlation value to consider
    	 * @param shift add this shift to positions
    	 * @param orderMaximums if true - order maximums in descending order, if false - leave them in increasing position (from back to foreground)
    	 * @param debugLevel debug level
    	 * @return variable length array or null, only defined values are output:null,
    	 *  {minPosition,maxPosition,firstMaxPosition,firstMaxValue},
    	 *  {minPosition,maxPosition,firstMaxPosition,firstMaxValue,secondMaxPosition,secondMaxValue}, ...
    	 */
    	double [] correlationResults(
    			double [] corr,
    			double threshold,
    			double shift,
    			boolean enableNegativeDisparity,
    			boolean orderMaximums,
    			int debugLevel
    			){
    		int size = (int) Math.round(Math.sqrt(corr.length));
    		int index=size*size/2; // start of the middle row;
    		int min=size, max=0,numMax=0;
    		boolean [] localMax=new boolean [size];
    		int scanMin=enableNegativeDisparity?0:(size/2-((int) shift));
    		if (scanMin<0) scanMin=0;
    		for (int i=0;i<size;i++) localMax[i]=false;
    		for (int i=scanMin;i<size;i++ ){
    			double d=corr[index+i];
    			if (d>threshold){

    				if (i>max) max=i;
    				if (i<min) min=i;
    				if ((i>0) && ((d>=corr[index+i-1]) || ((i==scanMin) && !enableNegativeDisparity)) && (i<(size-1)) && (d>corr[index+i+1]))  {
    					localMax[i]=true;
    					numMax++;
    				}
    				if (debugLevel>2) System.out.println("correlationResults(): d="+d+" threshold="+threshold+" i="+i+" (i+shift-size/2="+(i+shift-size/2)+") "+
    						" min="+min+" max="+max+ " numMax="+numMax);
    			}
    		}
    		if (numMax<1) return null;
    		double [] results= new double [2*(numMax+1)];
    		results[0]=min-(size/2)+shift;
    		results[1]=max-(size/2)+shift;
// first find all maximums, then reorder (relative maximum values may change
    		int iMax=0;
    		for (int i=scanMin;i<size;i++ ){
    			double c=corr[index+i];
    			if (c>threshold){
    				if (i>max) max=i;
    				if (i<min) min=i;

    				if (localMax[i]) {
    					if ((i==scanMin) && (corr[index+i-1]>corr[index+i])) {
    						if (debugLevel>2) System.out.println("correlationResults(): limiting negative disparity to 0");
    						results[2*(iMax+1)]=  i-size/2+shift;
    						results[2*(iMax+1)+1]=corr[index+i];
    					} else {
    						// y=ax^2+bx+c, x==0 for the middle point, c=y[0] already defined
    						double b=0.5*(corr[index+i+1]-corr[index+i-1]);
    						double a=0.5*(corr[index+i+1]+corr[index+i-1])-c; // negative
    						double x=-0.5*b/a;
    						results[2*(iMax+1)]=  i+x-size/2+shift;
    						results[2*(iMax+1)+1]=a*x*x+b*x+c;
    					}
    					iMax++;
    				}
    			}
    		}
    		// now sort the maximums
    		if (orderMaximums) {
    			if (debugLevel>2) System.out.println("correlationResults(): reordering maximums");
    			boolean ordered =false;
    			while (!ordered){
    				ordered=true;
    				for (int i=0; i<(numMax-1);i++) {
    	    			if (debugLevel>2) System.out.println("correlationResults(): results["+(2*i+3)+"]="+results[2*i+3]+
    	    					", results["+(2*i+5)+"]="+results[2*i+5]);

    					if (results[2*i+3]<results[2*i+5]){
    						ordered=false;
    						// swap maximums;
    						if (debugLevel>2) System.out.println("correlationResults(): swapping max "+i+" and "+(i+1));
    						double tmp=results[2*i+3];
    						results[2*i+3] = results[2*i+5];
    						results[2*i+5]=tmp;
    						tmp=results[2*i+2];
    						results[2*i+2] = results[2*i+4];
    						results[2*i+4]=tmp;
    					}

    				}
    			}
    		}
    		return results;
    	}

    	public double [][] collectDisparityFromTiles(
    			int debugLevel){
    		double [][] disparity=new double [6][this.mapWidth*this.mapHeight];
    		for (int n=0;n<disparity.length;n++) for (int i=0;i<disparity[0].length;i++)disparity[n][i]=Double.NaN;
    		int yMinLim=this.tilesY0;
    		int xMinLim=this.tilesX0;
    		int yMaxLim=this.tilesY0+this.tilesSize*this.tilesNVert;
    		int xMaxLim=this.tilesX0+this.tilesSize*this.tilesNHor;
    		if (yMinLim<0) yMinLim=0;
    		if (xMinLim<0) xMinLim=0;
    		if (yMaxLim>this.mapHeight) yMaxLim=this.mapHeight;
    		if (xMaxLim>this.mapWidth)  xMaxLim=this.mapWidth;
    		if (debugLevel>1) System.out.println("collectDisparityFromTiles(): this.tilesSize="+this.tilesSize+",yMinLim="+yMinLim+
    				", yMaxLim="+yMaxLim+", xMinLim="+xMinLim+", xMaxLim="+xMaxLim);
    		for (int iY=yMinLim;iY<yMaxLim;iY++){
    			int tY=iY-this.tilesY0;
    			for (int iX=xMinLim;iX<xMaxLim;iX++){
    				int tX=iX-this.tilesX0;
    				for (int level=this.tiles0.length-1;level>=0;level--){

    					int size = (this.tilesSize>>level);
						double [][] tile=this.tiles0[level][tY/size][tX/size];
    					if (tile!=null){
    						int index=iY*this.mapWidth+iX;
    						disparity[2][index]=tile[0][0];
    						disparity[3][index]=tile[0][1]; // contrast background
    						disparity[4][index]=tile[tile.length-1][0]; // last - FG if available, if not - BG
    						disparity[5][index]=tile[tile.length-1][1]; // contrast *-ground
    						disparity[0][index]=(disparity[3][index]>disparity[5][index])?disparity[2][index]:disparity[4][index];
    						disparity[1][index]=(disparity[3][index]>disparity[5][index])?disparity[3][index]:disparity[5][index];
    						break;
    					}
    				}
    			}
    		}
    		return disparity;
    	}

// just temporary to visualize tiles - actual disparity will be processed diffirently
    	public double [][] showDisparityFromTiles(
    			int debugLevel){
    		double [][] disparity=new double [8][this.mapWidth*this.mapHeight];
    		for (int n=0;n<disparity.length;n++) for (int i=0;i<disparity[0].length;i++)disparity[n][i]=Double.NaN;
    		int yMinLim=this.tilesY0;
    		int xMinLim=this.tilesX0;
    		int yMaxLim=this.tilesY0+this.tilesSize*this.tilesNVert;
    		int xMaxLim=this.tilesX0+this.tilesSize*this.tilesNHor;
    		if (yMinLim<0) yMinLim=0;
    		if (xMinLim<0) xMinLim=0;
    		if (yMaxLim>this.mapHeight) yMaxLim=this.mapHeight;
    		if (xMaxLim>this.mapWidth)  xMaxLim=this.mapWidth;
    		if (debugLevel>1) System.out.println("collectDisparityFromTiles(): this.tilesSize="+this.tilesSize+",yMinLim="+yMinLim+
    				", yMaxLim="+yMaxLim+", xMinLim="+xMinLim+", xMaxLim="+xMaxLim);
    		for (int side=0;side<this.tiles.length;side++) for (int self=0;self<this.tiles[side].length;self++) if (this.tiles[side][self]!=null){
    			for (int iY=yMinLim;iY<yMaxLim;iY++){
    				int tY=iY-this.tilesY0;
    				for (int iX=xMinLim;iX<xMaxLim;iX++){
    					int tX=iX-this.tilesX0;
    					for (int level=this.tiles[side][self].length-1;level>=0;level--){

    						int size = (this.tilesSize>>level);
    						double [] tile=this.tiles[side][self][level][tY/size][tX/size];
    						if (tile!=null){
// find single maximal value
    							double maxV=0;
    							int maxI=-1;
    							for (int i=0;i<tile.length;i++) if (tile[i]>maxV){
    								maxI=i;
    								maxV=tile[i];
    							}
    							if (maxI>=0){
    								double d=maxI;
    								double c=tile[maxI];
    								if ((maxI>0) && (maxI<(tile.length-1) && (tile[maxI-1]>0.0) && (tile[maxI+1]>0.0))){// do not interpolate if 0.0 (masked out)
    		    						// y=ax^2+bx+c, x==0 for the middle point, c=y[0] already defined
    		    						double b=0.5*(tile[maxI+1]-tile[maxI-1]);
    		    						double a=0.5*(tile[maxI+1]+tile[maxI-1])-c; // negative
    		    						double x=-0.5*b/a;
    		    						d+=x; // modify position
    		    						c+=a*x*x+b*x; // value at max
    								}
        							int index=iY*this.mapWidth+iX;
    								disparity[4*self+side  ][index]=d+this.minDisparity; // position
    								disparity[4*self+side+2][index]=c; // value (contrast)
    							}
    							break;
    						}
    					}
    				}
    			}
    		}
    		return disparity;
    	}

    	public float [][] createAmbiguityMaps(
    			int side,
    			double weightSobel,
    			double weightCb,
    			double weightCr,
    			double noiseLev,
    			boolean removePeriodic,
    			double minPeriod,
    			double minZeroAuto,
    			double minPeriodContrast,
    			double minAbsoluteFraction,
    			double blurMaskSigma,
    			double scale,

    			int debugLevel){

    		//TODO: *********** here combine with cross-correlation with auto correlation,
//			int self=0;


    		int numSamples= this.maxDisparity-this.minDisparity+1;
    		if (this.ambiguityMap==null){
    			this.ambiguityMap=new float[2][][];
    			for (int i=0;i<this.ambiguityMap.length;i++) this.ambiguityMap[i]=null;
    		}
//    		int width=this.mapWidth;
    		this.ambiguityMap[side]=new float [4][this.mapWidth*this.mapHeight];
    		for (int i=0;i<this.ambiguityMap[side].length;i++) for (int j=0;j<this.ambiguityMap[side][i].length;j++){
    			this.ambiguityMap[side][i][j]=Float.NaN;
    		}
    		for (int sectionY=1;sectionY<(this.mapHeight-1);sectionY++){
    			double [] fromTiles=tilesSection(
    					sectionY,
    					side,
    					0, //self,
    					debugLevel);
    			double [] fromPixelDiff=sectionPixelDifference(
    					sectionY,
    					side,
    					0, //self,
    	    			weightSobel,
    					weightCb,
    					weightCr,
    					noiseLev,
    					debugLevel);
    			double [] fromAlpha=tilesSection(
    					sectionY,
    					side,
    					0, //self,
    					debugLevel);
    			double [] fromTilesAuto=null;
    			double [] fromPixelDiffAuto=null;
    			double [] fromAlphaAuto=null;

    			if (removePeriodic) {
    				fromTilesAuto=tilesSection(
    						sectionY,
    						side,
    						1,
    						debugLevel);
    				fromPixelDiffAuto=sectionPixelDifference(
    						sectionY,
    						side,
    						1,
    						weightSobel,
    						weightCb,
    						weightCr,
    						noiseLev,
    						debugLevel);
    				fromAlphaAuto=tilesSection(
    						sectionY,
    						side,
    						1,
    						debugLevel);
    			}

    			for (int iX=1;iX<this.mapWidth;iX++){
    				double [] line =     new double[numSamples];
    				double [] autoLine = new double[numSamples];
    				double [] filteredLine=null;
    				boolean nonEmpty=false;
    				for (int iDisp=0;iDisp<numSamples;iDisp++){
    					int index=iX*numSamples+iDisp;
    					line[iDisp]=fromTiles[index]*fromPixelDiff[index]*fromAlpha[index]; // remove fromAlpha?
    					nonEmpty |= (line[iDisp]>0);
    					if (removePeriodic)autoLine[iDisp]=fromTilesAuto[index]*fromPixelDiffAuto[index]*fromAlphaAuto[index]; // remove fromAlpha? = new double[numSamples];
    				}
    				if (nonEmpty){
    					if (removePeriodic) {
    						filteredLine=removePeriodic(
    								line, // double [] imageCross,
    								autoLine, //double [] imageAuto,
    								minPeriod,
    								minZeroAuto,
    								minPeriodContrast,
    								minAbsoluteFraction,
    								blurMaskSigma,
    								scale,
    								debugLevel);
    						if (filteredLine!=null) line=filteredLine;
    					}
    					int [] iMax={-1,-1};
    					double [] max={0.0,0.0};
						double [] maxP={0.0,0.0};
						int index=sectionY*this.mapWidth+iX;
    					for (int nMax=0;nMax<iMax.length;nMax++){
    						if (nMax==0){ // find absolute maximum
    							for (int iDisp=0;iDisp<numSamples;iDisp++) if (line[iDisp]>max[nMax]){
    								max[nMax]=line[iDisp];
    								iMax[nMax]=iDisp;
    							}
    						} else { // find local maximum not equal to absolute one
    							for (int iDisp=0;iDisp<numSamples;iDisp++) if ((line[iDisp]>max[nMax]) && (iDisp!=iMax[0])) {
    								if (((iDisp==0) || (line[iDisp]>=line[iDisp-1])) && ((iDisp==(numSamples-1)) || (line[iDisp]>=line[iDisp+1]))){
    									max[nMax]=line[iDisp];
    									iMax[nMax]=iDisp;
    								}
    							}
    						}
    						if (max[nMax]<=0.0) break;
    							maxP[nMax]=iMax[nMax];
    							if ((iMax[nMax]>0) && (iMax[nMax]<(numSamples-1))){
    								// y=ax^2+bx+c, x==0 for the middle point, c=y[0] already defined
    								double b=0.5*(line[iMax[nMax]+1]-line[iMax[nMax]-1]);
    								double a=0.5*(line[iMax[nMax]+1]+line[iMax[nMax]-1])-max[nMax]; // negative
    								double x=-0.5*b/a;
    								maxP[nMax]+=x;
    								max[nMax]+= a*x*x+b*x; //+c)
    							}
    							maxP[nMax]+=this.minDisparity;
    							if (maxP[nMax]<0) maxP[nMax]=0.0; // limit disparity to 0.0
    							this.ambiguityMap[side][2*nMax][index]=  (float) maxP[nMax];
    							this.ambiguityMap[side][2*nMax+1][index]=(float) max[nMax];
    					}
    				}
    			}
    			IJ.showProgress(sectionY,this.mapHeight-2);
    		}
    		IJ.showProgress(1.0);
    		return this.ambiguityMap[side];
    	}

    	public float [] initialResolveMaps(
    			int side,
    			double threshold,
    			double minSecondFrac,
    			double edgesBonus,
    			int debugLevel){
    		if (this.resolvedMap==null){
    			this.resolvedMap=new float[2][];
    			for (int i=0;i<this.resolvedMap.length;i++) this.resolvedMap[i]=null;
    		}
    		initSobelY(false);
    		if (debugLevel>0) System.out.println("initialResolveMaps("+side+", "+threshold+", "+minSecondFrac+", "+edgesBonus+")");
    		this.resolvedMap[side]=new float[this.mapWidth*this.mapHeight];
    		int numLayers=this.overlapImages.length/2;
    		int iAlphaThis= (side==0)?0:numLayers;
    		double maxBonus=0;
    		for (int i=0;i<this.resolvedMap[side].length;i++) {
    			this.resolvedMap[side][i]=Float.NaN;
    			double bonus=(1.0+edgesBonus*this.sobelY[side][i])*this.overlapImages[iAlphaThis][i];
    			if (debugLevel>1){
    				if (maxBonus<bonus){
    					maxBonus=bonus;
    					System.out.println("initialResolveMaps(): x="+(i%this.mapWidth)+" y="+(i/this.mapWidth)+" sobel="+this.sobelY[side][i]+" bonus="+bonus);
    				}
    			}
    			if (!Float.isNaN(this.ambiguityMap[side][0][i]) && ((this.ambiguityMap[side][1][i]*bonus)>=threshold)) {
    				if ((Float.isNaN(this.ambiguityMap[side][2][i]) ||
    						(
//    								((this.ambiguityMap[side][3][i]*bonus)<threshold) &&
    								(this.ambiguityMap[side][3][i]<(minSecondFrac*this.ambiguityMap[side][1][i]))
    						))){
    					this.resolvedMap[side][i]=this.ambiguityMap[side][0][i];

    				} else {
    					this.resolvedMap[side][i]=-2; //ambiguous disparity
    				}
    			}
    		}
    		return this.resolvedMap[side];// (Float.NaN - undefined, >=0 - disparity, <-1 - ambiguity)
    	}

    	public float [] getResolvedState(
    			int side,
    			double threshold,
    			double minSecondFrac,
    			double edgesBonus,
    			double scale, // to match "brightness" to disparity values
    			int debugLevel){
    		if (debugLevel>0) System.out.println("getResolvedState("+side+", "+threshold+", "+minSecondFrac+", "+edgesBonus+")");
    		float [] resolvedState=new float[this.mapWidth*this.mapHeight];
    		int numLayers=this.overlapImages.length/2;
    		int iAlphaThis= (side==0)?0:numLayers;
    		double maxBonus=0;
    		for (int i=0;i<this.resolvedMap[side].length;i++) {
    			resolvedState[i]=0.0f;
    			double bonus=(1.0+edgesBonus*this.sobelY[side][i])*this.overlapImages[iAlphaThis][i];
    			if (debugLevel>1){
    				if (maxBonus<bonus){
    					maxBonus=bonus;
    					System.out.println("getResolvedState(): x="+(i%this.mapWidth)+" y="+(i/this.mapWidth)+" sobel="+this.sobelY[side][i]+" bonus="+bonus);
    				}
    			}

    			if (Float.isNaN(this.resolvedMap[side][i])){
        			if (
        					!Float.isNaN(this.ambiguityMap[side][2][i]) &&
        					((this.ambiguityMap[side][1][i]*bonus)>=threshold) &&
        					(
//        							(this.ambiguityMap[side][3][i]*bonus>=threshold) ||
        									(this.ambiguityMap[side][3][i]>=(minSecondFrac*this.ambiguityMap[side][1][i])))){
        				resolvedState[i]=(float) (2.0*scale); // ambiguity
        			}

    			} else {
    				resolvedState[i]=(float) (1.0*scale); // already resolved
    			}
    		}
    		return resolvedState;
    	}

    	public int [] fillVoidsStep(
    			int    side,
    			int    minNumNeib,
    			double minExistentAlpha,
    			double minExistentStrength,
    			double maxExistentEdge, //?? (do not cross edges
    			double weightCb,
    			double weightCr,
    			double maxToneDiff, // using weighted distance) if this disparity will be applied
    			double maxDisparityDifference, // among the neighbors
    			double minNewStrength,
    			int debugLevel){
 //   		int numLayers=this.overlapImages.length/2;
//    		int iAlphaThis= (side==0)?0:numLayers;
//    		int iAlphaOther=(side!=0)?0:numLayers;
			  int debugLevel0=debugLevel;
			  int debugX=280;
			  int debugY=970;
			  int debugDelta=10;
	    		int [] dirs={0,1,1+this.mapWidth,this.mapWidth,-1+this.mapWidth, -1, -1 -this.mapWidth, -this.mapWidth, 1-this.mapWidth};
	    		double [] weights={1.0/(1.0+weightCb+weightCr),weightCb/(1.0+weightCb+weightCr),weightCr/(1.0+weightCb+weightCr)};
				  if (debugLevel>1){
					  System.out.println("resolveAmbiguityStep():minExistentAlpha="+minExistentAlpha+" minExistentStrength="+minExistentStrength+
							  " weightCb="+weightCb+" weightCr="+weightCr+
							  " maxToneDiff="+maxDisparityDifference+" minNewStrength="+minNewStrength);
					  System.out.println("weights={"+weights[0]+", "+weights[1]+", "+weights[1]+"}");
				  }
	    		float [] newVals=new float [this.resolvedMap[side].length];
	    		for (int i=0;i<newVals.length;i++) newVals[i]=Float.NaN;
	    		for (int iY=1;iY<(this.mapHeight-1);iY++) { // to reduce number of checks for oob
	    			for (int iX=1;iX<=(this.mapWidth-1);iX++) { // to reduce number of checks for oob
	    				int index=iY*this.mapWidth+iX;
	    				debugLevel=debugLevel0+(((iY>=(debugY-debugDelta))&&(iY<=(debugY+debugDelta))&&
	    						(iX>=(debugX-debugDelta))&&(iX<=(debugX+debugDelta)))?1:0);
	    				if (Float.isNaN(this.resolvedMap[side][index])){ // looking for the undefined pixels
//	    					int numDefinedNeighbors=0;
//	    					double [] neighbors=new double [dirs.length];
	    					for (int dir=0;dir<dirs.length;dir++){

	    					}

	    				}
	    			}
	    		}
    	return null;
    	}




    	public int [] resolveAmbiguityStep(
    			int side,
    			double minExistentAlpha,
    			double minExistentStrength,
    			double weightCb,
    			double weightCr,
    			double maxToneDiff, // using weighted distance)
    			double maxDisparityDifference,
    			double minNewStrength,
    			int debugLevel){
    		int numLayers=this.overlapImages.length/2;
    		int iAlphaThis= (side==0)?0:numLayers;
    		int iAlphaOther=(side!=0)?0:numLayers;
			  int debugLevel0=debugLevel;
			  int debugX=280;
			  int debugY=970;
			  int debugDelta=10;


    		int [] dirs={0,1,1+this.mapWidth,this.mapWidth,-1+this.mapWidth, -1, -1 -this.mapWidth, -this.mapWidth, 1-this.mapWidth};
    		double [] weights={1.0/(1.0+weightCb+weightCr),weightCb/(1.0+weightCb+weightCr),weightCr/(1.0+weightCb+weightCr)};
			  if (debugLevel>1){
				  System.out.println("resolveAmbiguityStep():minExistentAlpha="+minExistentAlpha+" minExistentStrength="+minExistentStrength+
						  " weightCb="+weightCb+" weightCr="+weightCr+
						  " maxToneDiff="+maxDisparityDifference+" minNewStrength="+minNewStrength);
				  System.out.println("weights={"+weights[0]+", "+weights[1]+", "+weights[1]+"}");
			  }
    		byte [] mods=new byte [this.resolvedMap[side].length];
    		for (int i=0;i<mods.length;i++) mods[i]=0; // 0 - no change, 1 - use first value, 2 - use second value, 3 bad (does not fit)
    		for (int iY=2;iY<(this.mapHeight-2);iY++) { // to reduce number of checks for oob
    			for (int iX=2;iX<=(this.mapWidth-2);iX++) { // to reduce number of checks for oob
    				int index=iY*this.mapWidth+iX;
    				debugLevel=debugLevel0+(((iY>=(debugY-debugDelta))&&(iY<=(debugY+debugDelta))&&
    						(iX>=(debugX-debugDelta))&&(iX<=(debugX+debugDelta)))?1:0);

    				if (this.resolvedMap[side][index]>=0.0){
    					for (int iDir0=0;iDir0<dirs.length;iDir0++) {
    						int index1=index+dirs[iDir0];
    						if ((this.resolvedMap[side][index1]<-1.5) && (mods[index1]==0)){ // found ambiguous pixel, see if it matches and find the best parent
    							double [] existentAlphas=     new double [dirs.length];
    							double [] existentDisparities=new double [dirs.length];
    							int [] optionUsed=new int [dirs.length];
    							double [] existentStrength=   new double [dirs.length];
    							double [] weightredToneDiff=  new double [dirs.length];
    							double [] disparityDiff=      new double [dirs.length];
    							double [] newStrength=        new double [dirs.length];
    							int [] optionNew=new int [dirs.length]; // which of the two disparities is closer to the prototype
    							boolean [] newOptions={false,false};
    							for (int iDir=0;iDir<dirs.length;iDir++) {
    								int index2=index1+dirs[iDir];
    								if (this.resolvedMap[side][index2]>=0.0){
    									if (debugLevel>2){
    										System.out.println("resolveAmbiguityStep(): iY="+iY+" iX="+iX+" iDir0="+iDir0+
    												" index="+index+" index1="+index1+" index2="+index2);
    									}
    									existentDisparities[iDir]=this.resolvedMap[side][index2];
    									int dispar=((int) Math.round(existentDisparities[iDir]));
    									if (((index2 % this.mapWidth)+dispar)<this.mapWidth) {
    										existentAlphas[iDir]=
    											this.overlapImages[iAlphaThis][index2]*
    											this.overlapImages[iAlphaOther][index2+((side==0)?dispar:-dispar)];
    									} else {
    										continue; // 0 alpha
    									}
    									if (existentAlphas[iDir]<minExistentAlpha){
        									if (debugLevel>2){
        										System.out.println("resolveAmbiguityStep(): existentAlphas[iDir]="+existentAlphas[iDir]+" < "+minExistentAlpha);
        									}
    										continue; // too small: alpha
    									}
    									if (Float.isNaN(this.ambiguityMap[side][2][index2])) optionUsed[iDir]=0;
    									else {
    										optionUsed[iDir]=
    											(Math.abs(this.resolvedMap[side][index2]-this.ambiguityMap[side][0][index2])<
    													Math.abs(this.resolvedMap[side][index2]-this.ambiguityMap[side][2][index2]))?0:1;
    									}
    									existentStrength[iDir]=this.ambiguityMap[side][2*optionUsed[iDir]+1][index2];
    									if (existentStrength[iDir]<minExistentStrength){
        									if (debugLevel>2){
        										System.out.println("resolveAmbiguityStep(): existentStrength[iDir]="+existentStrength[iDir]+" < "+minExistentStrength);
        									}
    										continue; // too small: existentStrength
    									}
    									weightredToneDiff[iDir]=0.0;
    									for (int compNum=0;compNum<3;compNum++){
    										double d=this.overlapImages[iAlphaThis+compNum+1][index2]- // parent point
    										this.overlapImages[iAlphaThis+compNum+1][index1] ;// this, new point
    										weightredToneDiff[iDir]+=weights[compNum]*d*d;
    									}
    									weightredToneDiff[iDir]=Math.sqrt(weightredToneDiff[iDir]);
    									if (weightredToneDiff[iDir]>maxToneDiff){
        									if (debugLevel>2){
        										System.out.println("resolveAmbiguityStep(): weightredToneDiff[iDir]="+weightredToneDiff[iDir]+" > "+maxToneDiff);
        									}
    										continue; // too large: weightredToneDiff
    									}
    									// both disparities should be present for the new point
    									optionNew[iDir]=
    										(Math.abs(this.resolvedMap[side][index2]-this.ambiguityMap[side][0][index1])<
    												Math.abs(this.resolvedMap[side][index2]-this.ambiguityMap[side][2][index1]))?0:1;
    									disparityDiff[iDir]=Math.abs(this.resolvedMap[side][index2]-
    											this.ambiguityMap[side][2*optionUsed[iDir]][index1]);
    									if (disparityDiff[iDir]>maxDisparityDifference){
        									if (debugLevel>2){
        										System.out.println("resolveAmbiguityStep(): disparityDiff[iDir]="+disparityDiff[iDir]+" > "+maxDisparityDifference);
        									}
    										continue; // too large: disparityDiff
    									}
    									newStrength[iDir]= this.ambiguityMap[side][2*optionUsed[iDir]+1][index1];
    									if (newStrength[iDir]<minNewStrength){
        									if (debugLevel>2){
        										System.out.println("resolveAmbiguityStep(): newStrength[iDir]="+newStrength[iDir]+" > "+minNewStrength);
        									}
    										continue; // too small: newStrength
    									}
    									newOptions[optionNew[iDir]]=true;
    								}
    							} // for (int iDir=0;iDir<dirs.length;iDir++)
    							if (newOptions[0]&&newOptions[1]){

    								if (debugLevel>1){
    									System.out.println("resolveAmbiguityStep():iY="+iY+" iX="+iX+" index1="+index1+" iDir0="+iDir0+
    									" different parents suggested different disparities options, both match");
    								}
    								mods[index1]=3; // disable
    							} else if (newOptions[0]) mods[index1]=1; // use first disparity
    							else if (newOptions[1]) mods[index1]=2; // use second disparity
    							else mods[index1]=3; // disable second disparity
								if (debugLevel>2){
									System.out.println("resolveAmbiguityStep(): newOptions[0]="+newOptions[0]+
											" newOptions[1]="+newOptions[1]+
											" setting mod["+index1+"] (X="+(index1%this.mapWidth)+" Y="+(index1 / this.mapWidth)+
											" to "+mods[index1]);
								}
    						}
    					}
    				}
    			}

    		}
    		// Apply all changes;
    		int numFailed=0;
    		int numResol=0;
    		for (int i=0;i<mods.length;i++) {
    			if ((mods[i]>0) && (mods[i]<3)){
    		    			this.resolvedMap[side][i]=(mods[i]==1)?this.ambiguityMap[side][0][i]:this.ambiguityMap[side][2][i];
    		    			numResol++;
    			} else if (mods[i]==3){
    				numFailed++;
    			}
    		}
    		int [] result={numResol,numFailed};
    		return result;
    	}


    	/**
    	 * Combine pixel matching and correlation into disparity map
    	 * @param side 0 - left image, 1 - right
    	 * @param self 0 - with other image, 1 - with self (to detect periodic patterns)
    	 * @param weightCb weight of the image Cb component as fraction of Y
    	 * @param weightCr weight of the image Cr component as fraction of Y
    	 * @param noiseLev pixel noise as a fraction of the full scale (1.0)
    	 * @param debugLevel  Debug level
    	 * @return array of the same format as original image
    	 */
    	public double [] getFrameDisparity(
    			int side,
    			int self,
    			double weightSobel,
    			double weightCb,
    			double weightCr,
    			double noiseLev,
    			int debugLevel){
    		int numSamples= this.maxDisparity-this.minDisparity+1;
    		double [] map=new double [this.mapWidth*this.mapHeight];
    		for (int i=0;i<map.length;i++) map[i]=Double.NaN;
    		for (int sectionY=1;sectionY<(this.mapHeight-1);sectionY++){
    			double [] fromTiles=tilesSection(
    					sectionY,
    					side,
    					self,
    					debugLevel);
    			double [] fromPixelDiff=sectionPixelDifference(
    					sectionY,
    					side,
    					self,
    	    			weightSobel,
    					weightCb,
    					weightCr,
    					noiseLev,
    					debugLevel);
    			double [] fromAlpha=tilesSection(
    					sectionY,
    					side,
    					self,
    					debugLevel);
    			for (int iX=1;iX<this.mapWidth;iX++){
    				double [] line = new double[numSamples];
    				boolean nonEmpty=false;
    				for (int iDisp=0;iDisp<numSamples;iDisp++){
    					int index=iX*numSamples+iDisp;
    					line[iDisp]=fromTiles[index]*fromPixelDiff[index]*fromAlpha[index];
    					nonEmpty |= (line[iDisp]>0);
    				}
    				if (nonEmpty){
    					int iMax=-1;
    					double max=0.0;
    					for (int iDisp=0;iDisp<numSamples;iDisp++) if (line[iDisp]>max){
    						max=line[iDisp];
    						iMax=iDisp;
    					}
    					if (max>0.0){
    						double maxP=iMax;
    						if ((iMax>0) && (iMax<(numSamples-1))){
        						// y=ax^2+bx+c, x==0 for the middle point, c=y[0] already defined
        						double b=0.5*(line[iMax+1]-line[iMax-1]);
        						double a=0.5*(line[iMax+1]+line[iMax-1])-max; // negative
        						double x=-0.5*b/a;
        						maxP+=x;
    						}
    						maxP+=this.minDisparity;
    						if (maxP<0) maxP=0.0; // limit disparity to 0.0
    						map[sectionY*this.mapWidth+iX]=maxP;
    					}
    				}
    			}
    			IJ.showProgress(sectionY,this.mapHeight-2);
    		}
    		IJ.showProgress(1.0);
    		return map;
    	}



    	/**
    	 * Generate "cross-section" of the tiles data
    	 * @param sectionY Position of the section
    	 * @param side 0 - left image, 1 - right
    	 * @param debugLevel Debug level
    	 * @return 2d packed array - contrast of correlation for each disparity first,  each pixel - second
    	 */

    	public double[]  tilesSection(
    			int sectionY,
    			int side,
    			int self, // 1 - autocorrelation
    			int debugLevel){

    		int sectionLength=this.mapWidth;
    		int numSamples= this.maxDisparity-this.minDisparity+1;
    		double [] result = new double[sectionLength*numSamples];
    		for (int i=0;i<result.length;i++) result[i]=Double.NaN;
    		int xMinLim=this.tilesX0;
    		int xMaxLim=this.tilesX0+this.tilesSize*this.tilesNHor;
    		if (xMinLim<0) xMinLim=0;
    		if (xMaxLim>this.mapWidth)  xMaxLim=this.mapWidth;
    		int tY=sectionY-this.tilesY0;
    		for (int iX=xMinLim;iX<xMaxLim;iX++){
    			int tX=iX-this.tilesX0;
    			double [] tile=null;
    			int size;
    			for (int level=this.tiles[side][self].length-1;level>=0;level--){
    				size = (this.tilesSize>>level);
    				tile=this.tiles[side][self][level][tY/size][tX/size];
    				if (tile!=null) break;
    			}
    			if (tile!=null){
    				int index=numSamples*iX;
    				for (int i=0;i<numSamples;i++) result[index++]=tile[i];
    			}
    		}
    		return result;
    	}

    	/**
    	 * Calculate multiplied alpha of the two images along the section for the same array as in tilesSection
    	 * @param sectionY Position of the section
    	 * @param side 0 - left image, 1 - right
    	 * @param debugLevel Debug level
    	 * @return 2d packed array - contrast of correlation for each disparity first,  each pixel - second
    	 */
    	public double[]  alphaSection(
    			int sectionY,
    			int side,
    			int self,
    			int debugLevel){
    		int sectionLength=this.mapWidth;
    		int numSamples= this.maxDisparity-this.minDisparity+1;
    		double [] result = new double[sectionLength*numSamples];
    		for (int i=0;i<result.length;i++) result[i]=0.0;
    		int numLayers=this.overlapImages.length/2;
    		int iAlphaThis= (side==0)?0:numLayers;
    		int iAlphaOther=((side!=0)^(self!=0))?0:numLayers;
    		for (int iX=0;iX<sectionLength;iX++)for (int iDispar=0;iDispar<numSamples;iDispar++){
    			int iXOther=(side==0)?(iX-iDispar):(iX+iDispar);
    			if ((iXOther>=0) && (iXOther<sectionLength)){
    				result[iX*numSamples+iDispar]=
    					this.overlapImages[iAlphaThis][sectionY*sectionLength+iX]*
    					this.overlapImages[iAlphaOther][sectionY*sectionLength+iXOther];
    			}
    		}
    		return result;
    	}

    	public void initSobelY(
    			boolean uncoditionally){
    		if (uncoditionally || (this.sobelY==null)){
    			this.sobelY=new double[2][];
    			this.sobelY[0]=null;
    			this.sobelY[1]=null;
    		}
    		double [] preImage=null;
    		if (this.overlapImages[1]!=null){
    			if (this.preSobelSigma==0){
    				preImage=this.overlapImages[1];
    			}else {
    				preImage=this.overlapImages[1].clone();
    				(new DoubleGaussianBlur()).blurDouble(
    						preImage,
    						this.mapWidth,
    						this.mapHeight,
    						this.preSobelSigma,
    						this.preSobelSigma,
    						0.01);
    			}
    			this.sobelY[0]=sobelFilter(
    					preImage,
    					this.mapWidth,
    					this.edgeSpreadOrtho,
    					this.edgeSpreadDiagonal
    					);
    		}
    		if (this.overlapImages[1+this.overlapImages.length/2]!=null){
    			if (this.preSobelSigma==0){
    				preImage=this.overlapImages[1+this.overlapImages.length/2];
    			}else {
    				preImage=this.overlapImages[1+this.overlapImages.length/2].clone();
    				(new DoubleGaussianBlur()).blurDouble(
    						preImage,
    						this.mapWidth,
    						this.mapHeight,
    						this.preSobelSigma,
    						this.preSobelSigma,
    						0.01);
    			}

    			this.sobelY[1]=sobelFilter(
    					preImage,
    					this.mapWidth,
    					this.edgeSpreadOrtho,
    					this.edgeSpreadDiagonal
    					);
    		}
    	}
    	public boolean [] thresholdEdges(
    			double [] imageSrc,
    			double [] alpha, // or null
    			int imageWidth,
    			double thresholdHigh,
    			double thresholdLow,
    			boolean showProgress,
    			int debugLevel){
    		int updateCounts=1000;
    		int imageHeight=imageSrc.length/imageWidth;
   			List <Integer> edgeList=new ArrayList<Integer>(10000);
   			Integer Index;
//     		int [] dirsI={0,1,1+imageWidth,imageWidth,-1+imageWidth, -1, -1 -imageWidth, -imageWidth, 1-imageWidth};
			int [][]dirs={{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}}; // x,y
			double [] image=imageSrc;
			long startTime = System.nanoTime();
			if (alpha!=null){
				image=imageSrc.clone();
				for (int i=0;i<image.length;i++) image[i]*=alpha[i];
			}
			boolean [] edges=new boolean [image.length];
			for (int i=0;i<edges.length;i++) edges[i]=false;
			// create initial list
			for (Index=0;Index<image.length;Index++){
				if (image[Index]>=thresholdHigh){
					edges[Index]=true;
					edgeList.add(Index);
				}
			}
			if (debugLevel>0){
				System.out.println("Created initial list, length="+edgeList.size()+"  at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));

			}
			int maxLength=0;
			int updateCount=updateCounts;
   			while (edgeList.size()>0) {
   				Index=edgeList.remove(0); // faster will be to use get()?
   				int x=Index%imageWidth;
   				int y=Index/imageWidth;
   				for (int dir=0;dir<dirs.length;dir++){
   					int x1=x+dirs[dir][0];
   					int y1=y+dirs[dir][1];
   					if ((x1>=0) && (x1<imageWidth) && (y1>=0) && (y1<imageHeight)){
   						Index=x1+y1*imageWidth;
   						if (!edges[Index] && (image[Index]>=thresholdLow)){
   							edges[Index]=true;
   							edgeList.add(Index);
   						}
   					}
   				}
   	   			if (showProgress){
   	   				if (maxLength<edgeList.size()) maxLength=edgeList.size();
   	   				if (--updateCount<=0){
   	   				updateCount=updateCounts;
   	   				   IJ.showProgress(maxLength-edgeList.size(),maxLength+1);
   	   				   IJ.showStatus("length:"+edgeList.size()+" max length:"+maxLength);
   	   				}
   	   			}
   			}
   			if (showProgress)IJ.showProgress(1.0);
			if (debugLevel>0){
				System.out.println("Thresholded image at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));

			}

   			return edges;
    	}

//TODO: make an approximate version, that after making first N steps is looking when the parallel branches will come
    	// as close as 1 orthogonal (not diagonal) pixel apart. THen calculate lengths and areas (adding missing connection)

    	private int commonAncestor(
    			int index1,
    			int index2,
    			int [][] topology,
    			double [] timeToReach,
    			byte [] parents,
    			double [][][] node,// to find the best parent if there are many of them
    			int debugLevel){
    		double ttr1,ttr2;
    		int pDir;
    		if (debugLevel>2) System.out.println("\n---commonAncestor("+index1+", "+index2+", ...):");
    		while ((index1!=index2) && (index1>=0) && (index2>=0)) {
    			 ttr1=timeToReach[index1];
    			 ttr2=timeToReach[index2];
    			if (ttr1>ttr2){
    				if (node[index1]!=null) pDir=(int) node[index1][8][0];
    				else for (pDir=0;pDir<8;pDir++) if ((parents[index1]& (1<<pDir))!=0) break;
    	    		if (debugLevel>2) System.out.print (" parents["+index1+"]="+String.format("0x%02x",parents[index1])+" pDir="+pDir);
    				index1=topology[index1][pDir];
    	    		if (debugLevel>2) System.out.println (" index1="+index1);
    			} else {
    				if (node[index2]!=null) pDir=(int) node[index2][8][0];
    				else for (pDir=0;pDir<8;pDir++) if ((parents[index2]& (1<<pDir))!=0) break;
    	    		if (debugLevel>2) System.out.print (" parents["+index2+"]="+String.format("0x%02x",parents[index2])+" pDir="+pDir);
    				index2=topology[index2][pDir];
    	    		if (debugLevel>2) System.out.println (" index2="+index2);
    			}
    		}
    		if (debugLevel>2) System.out.println( ((index1==index2)?index1:-1)+ " => (index1="+index1+", index2="+index2);
    		return (index1==index2)?index1:-1;
    	}


    	//TODO: enhance nodes by "re-routing"  to the common point (first run wave from the approximate center for some distance, then go back and find the point
    	// with the smallest sum of times from each of the ends
    	public void thinEdge(
    			ArrayList<Integer> edgeList,
    			double [] image,
    			boolean [] bitmapAccumulator, // array were the edges segment will be drawn (over existent ones)
    			int imageWidth,
    			int minCycleArea, // minimal area of the closed loop to keep
    			int endRadius,
    			//TODO: - other parameters to trim short branches and detect ends?
    			boolean showProgress,
    			float [][] dbgImage, // should have length of 7 (or null)
    			int debugLevel0
    			){
//     		int updateCounts=10;

//     		int debugXMin=140;
//     		int debugXMax=155;
//     		int debugYMin=60;
//     		int debugYMax=75;
     		int debugLevel=debugLevel0;
 //    		int maxLength=0;
 //    		int updateCount=updateCounts;
    		int imageHeight=image.length/imageWidth;
    		int topologySize=edgeList.size();
    		int [][] topology = new int [topologySize][8];
			int [][]dirs={{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}}; // x,y
			double [] distance={1.0,Math.sqrt(2.0),1.0,Math.sqrt(2.0),1.0,Math.sqrt(2.0),1.0,Math.sqrt(2.0)};
			// Next code is extremely slow
			if (debugLevel>1) System.out.println("thinEdge():0");
			/*
			for (int i=0;i<topology.length;i++){
				Integer Index=edgeList.get(i);
   				int x=Index%imageWidth;
   				int y=Index/imageWidth;
   				for (int dir=0;dir<dirs.length;dir++){
   					int x1=x+dirs[dir][0];
   					int y1=y+dirs[dir][1];
   					if ((x1>=0) && (x1<imageWidth) && (y1>=0) && (y1<imageHeight)){
   						Index=x1+y1*imageWidth;
   	   					topology[i][dir]=edgeList.indexOf(Index);
   					} else topology[i][dir]=-1;
   				}
			}
			*/
			int [] topologyIndex=new int[image.length];
			for (int i=0;i<topologyIndex.length;i++) topologyIndex[i]=-1;
			for (int i=0;i<topology.length;i++) {
				Integer Index=edgeList.get(i);
				topologyIndex[Index]=i;
			}
			for (int i=0;i<topology.length;i++){
				Integer Index=edgeList.get(i);
   				int x=Index%imageWidth;
   				int y=Index/imageWidth;
   				for (int dir=0;dir<dirs.length;dir++){
   					int x1=x+dirs[dir][0];
   					int y1=y+dirs[dir][1];
   					if ((x1>=0) && (x1<imageWidth) && (y1>=0) && (y1<imageHeight)){
   	   					topology[i][dir]=topologyIndex[x1+y1*imageWidth];
   					} else topology[i][dir]=-1;
   				}
			}
			topologyIndex=null;



			double [] speed      = new double[topologySize]; // edge analog value, speed of wave
			double [] timeToReach= new double[topologySize]; // value on the wave front (when got from a single parent)
			int    [] pathArea   = new int   [topologySize]; // signed area of the path (YdX-XdY), when it comes from a single parent. Clockwise cycles positive, twice number of pixels
			byte   [] parents    = new byte  [topologySize]; // bit mask of parents (regular cells have just one)
			short  [] pixelX=      new short [topologySize]; // pixel X coordinate (just to calculate cycle area)
			short  [] pixelY=      new short [topologySize]; // pixel Y coordinate (just to calculate cycle area)
			short  [] pathLength=  new short [topologySize]; // length of the path, counting diagonal and orthogonal as 1.
			// for paths starting and anding in the same points,  pathLength1+pathlength2-2=abs(pathArea1-pathArea2) when there are no cycles
			if (debugLevel>1) System.out.println("thinEdge()");
			double [][][] node  = new double [topologySize][][]; // all cells but nodes have null here, nodes have {value, areas, lengths} triplets, the last row [8] - best parent
			for (int i=0;i<topologySize;i++){
				Integer Index=edgeList.get(i);
				pixelX[i]=(short) (Index%imageWidth);
				pixelY[i]=(short) (Index/imageWidth);

				speed[i]=image[Index];
				timeToReach[i]=0;
				pathArea[i]=0;
				pathLength[i]=0;
				parents[i]=0;
				node[i]=null;
			}
			if (debugLevel>1) System.out.println("thinEdge():1");
			if (debugLevel>2) {
				System.out.println("Topology:");
				for (int i=0;i<1000;i++) {
					System.out.print(String.format("%03d:" , i));

					for (int j=0;j<8;j++)   System.out.print(String.format(" %05d" , topology[i][j]));
					System.out.println(" X="+pixelX[i]+" Y="+pixelY[i]);
				}
			}
			//TODO: find a better start (i.e. run a wave from 0, but that still can result in a point on a circle)
			Integer start=0;
			timeToReach[0]=1.0; //
			ArrayList<Integer> waveList=new ArrayList<Integer>(topologySize);
			waveList.add(start);
			int processedCells=0;
   			while (waveList.size()>0) {
   				double min=timeToReach[waveList.get(0)];
   				int earliest=0;
   				for (int i=1;i<waveList.size();i++){
   					int wi=waveList.get(i);
   					if (min>timeToReach[wi]){
   						earliest=i;
   						min=timeToReach[wi];
   					}
   				}
   				Integer waveIndex=waveList.remove(earliest);
   				debugLevel=debugLevel0+(((pixelX[waveIndex]>=this.debugXMin) && (pixelX[waveIndex]<=this.debugXMax) &&
   						(pixelY[waveIndex]>=this.debugYMin) && (pixelY[waveIndex]<=this.debugYMax))?1:0);
   				if (debugLevel>2) System.out.println("\n>>>> pixelX["+waveIndex+"]="+pixelX[waveIndex]+" pixelY["+waveIndex+"]="+pixelY[waveIndex]);
   				for (int dir=0;dir<8;dir++){
   					Integer newIndex=topology[waveIndex][dir];
   					if (newIndex>=0){
   						double newValue=2.0*distance[dir]/(speed[waveIndex]+speed[newIndex])+timeToReach[waveIndex];
   						int newArea=pathArea[waveIndex]+
   						pixelY[waveIndex]*(pixelX[newIndex]-pixelX[waveIndex])-
   						pixelX[waveIndex]*(pixelY[newIndex]-pixelY[waveIndex]);
   						short newLength=(short) (pathLength[waveIndex]+1);
   						int dirToParent=((dir+4)%8);
   						byte parentMask=(byte)(1<<dirToParent); // direction opposite to dir
   		   				if (debugLevel>2) System.out.print("dir="+dir+" dirToParent="+dirToParent+String.format(" parentMask=0x%x",parentMask)+
   		   					 " newIndex: pixelX["+newIndex+"]="+pixelX[newIndex]+" pixelY["+newIndex+"]="+pixelY[newIndex]);
   						if (timeToReach[newIndex]==0.0){ // new cell - fill it
   							timeToReach[newIndex]=newValue;
   							parents[newIndex]=parentMask;
   							pathArea[newIndex]=newArea;
   							pathLength[newIndex]=newLength;
   							waveList.add(newIndex);
   			   	   			processedCells++;
   			    			if (showProgress)IJ.showProgress(processedCells,topologySize);
   			    			if (debugLevel>2) System.out.println( " --NEW CELL--");
   						} else {
   							boolean improve=timeToReach[newIndex]>newValue;
   							// as we take the smallest value from the list, all cells that can be improved do not have their children yet. But may have multiple parents
   							if (improve){
   								int cycleWithBestArea=newArea-pathArea[newIndex];
   								if (cycleWithBestArea<0) cycleWithBestArea=-cycleWithBestArea;
   								int commonIndex=commonAncestor(
   										waveIndex,
   						    			newIndex,
   						    			topology,
   						    			timeToReach,
   						    			parents,
   						    			node,// to find the best parent if there are many of them
   						    			debugLevel);

   								cycleWithBestArea-=(newLength+pathLength[newIndex]-2*(1+pathLength[commonIndex]));// should be 0 if the paths have no gaps - verify?
   	   			    			if (debugLevel>2) System.out.print( " -- OLD CELL, IMPROVE, cycleWithBestArea="+cycleWithBestArea+" newArea="+newArea+" pathArea["+newIndex+"]="+pathArea[newIndex]+
   	   			    				" newLength="+newLength+" pathLength["+newIndex+"]="+pathLength[newIndex]+" common: pathLength["+commonIndex+"]="+pathLength[commonIndex]);
   								if (cycleWithBestArea<=minCycleArea){ // replace the only or one of the older parents
   									if (node[newIndex]==null){ // just replace, no need to add as it is (still) in the list
   										parents[newIndex]=parentMask;
   									} else { // already has multiple parents - remove those that do not qualify cycles criteria
   										parents[newIndex]|=parentMask; // add new parent
   										node[newIndex][8][0]=dirToParent;
   									}
   								} else { // if (cycleWithBestArea<=minCycleArea) - add cycle
   	   	   			    			if (debugLevel>2) System.out.print( " -- OLD CELL, WORSE--");
   									if (node[newIndex]==null){
   		   	   			    			if (debugLevel>2) System.out.print( " -- not a node, setting a new node ");
   										node[newIndex]=new double [9][3];
   										int pDir=0;
   										for (int i=0;i<8;i++) if ((parents[newIndex] &(1<<i))!=0) { // should be done before parents[newIndex]|=parentMask;
   											pDir=i;
   											break;
   										}
   		   	   			    			if (debugLevel>2) System.out.print( " -- setting parent of a node, pDir= "+pDir+
   		   	   			    					", saving old values: ttr="+IJ.d2s(timeToReach[newIndex],3)+" area="+IJ.d2s(pathArea[newIndex],3)+" length="+IJ.d2s(pathLength[newIndex],3)  );
   										// Save the only data set to the parents data set
   										node[newIndex][pDir][0]=timeToReach[newIndex];
   										node[newIndex][pDir][1]=pathArea[newIndex];
   										node[newIndex][pDir][2]=pathLength[newIndex];
   										//  									} else { // if (node[newIndex]==null)

   									} // if (node[newIndex]==null)
   									parents[newIndex]|=parentMask;
   									node[newIndex][8][0]=dirToParent;
		   	   			    			if (debugLevel>2) System.out.println(String.format(" Set node[%d][8][0]=%d", newIndex,dirToParent));
   								}
   								// parents[newIndex] already has new (this) parent
   								if (debugLevel>2) System.out.println();
   								timeToReach[newIndex]=newValue;
   								pathArea[newIndex]=newArea;
   								pathLength[newIndex]=newLength;
   								if (node[newIndex]!=null){ // save the new data to node array, check/remove old parent(s
   									if (debugLevel>2) System.out.print( " -- node is not null, saving new values to node["+newIndex+"]["+dirToParent+"] ");
   									node[newIndex][dirToParent][0]=newValue;
   									node[newIndex][dirToParent][1]=newArea;
   									node[newIndex][dirToParent][2]=newLength;
   									for (int pDir=0;pDir<8;pDir++) if ((pDir!=dirToParent) && ((parents[newIndex] &(1<<pDir))<=0)){
   										//int parentIndex=topology[newIndex][pDir];
   										int cycleWithParentArea=newArea-(int) node[newIndex][pDir][1];
   										if (cycleWithParentArea<0) cycleWithParentArea=-cycleWithParentArea;
   										cycleWithParentArea-=(newLength+ (int) node[newIndex][pDir][2]-2);// should be 0 if the paths have no gaps - verify?

   		   	   			    			if (debugLevel>2) System.out.println( " Comparing paths in the node, cycleWithBestArea="+cycleWithBestArea+" newArea="+newArea+
   		   	   			    					" node["+newIndex+"]["+pDir+"][1]="+node[newIndex][pDir][1]+
   		    	   			    				" newLength="+newLength+" node["+newIndex+"]["+pDir+"][2]="+node[newIndex][pDir][2]);

   		   	   			    			if (cycleWithParentArea<minCycleArea){ // remove that parent
   		   	   			    		if (debugLevel>2) System.out.print( " Removed parent "+pDir+", was "+String.format("0x%02x",parents[newIndex]));
   											parents[newIndex] &= ~(1<<pDir);
   		   		   	   			    		if (debugLevel>2) System.out.print( " now: "+String.format("0x%02x",parents[newIndex]));
   										}
   									}
   									// That probably should not happen remove  - number of parents could get down to 1:
   									int numParents=0;
   									for (int pDir=0;pDir<8;pDir++) if ((parents[newIndex] &(1<<pDir))!=0) numParents++; // including this new one
   		   	   			    		if (debugLevel>2) System.out.println( " Now there are "+numParents+" parents fro the "+newIndex);
   									if (numParents<2) {
   										if (debugLevel>0){
   											System.out.println("Removing parents array from former node at X="+pixelX[newIndex]+", Y="+pixelY[newIndex]);
   										}
   										node[newIndex]=null;
   									}
   								}
   							}  else {// if (improve). Maybe we've got a valid cycle, even if not came here first?
		   	   			    		if (debugLevel>2) System.out.print( " No improvement for "+newIndex+" --> ");
   								if (node[newIndex]==null){
		   	   			    		if (debugLevel>2) System.out.print( " Node did not exist");
   	   								int cycleWithBestArea=newArea-pathArea[newIndex];
   	   								if (cycleWithBestArea<0) cycleWithBestArea=-cycleWithBestArea;
   	   								int commonIndex=commonAncestor(
   	   										waveIndex,
   	   						    			newIndex,
   	   						    			topology,
   	   						    			timeToReach,
   	   						    			parents,
   	   						    			node,// to find the best parent if there are many of them
   	   						    			debugLevel);
   	   								cycleWithBestArea-=(newLength+pathLength[newIndex]-2*(1+pathLength[commonIndex]));// should be 0 if the paths have no gaps - verify?
   	   	   			    			if (debugLevel>2) System.out.println( " cycleWithBestArea="+cycleWithBestArea+" newArea="+newArea+" pathArea["+newIndex+"]="+pathArea[newIndex]+
   	    	   			    				" newLength="+newLength+" pathLength["+newIndex+"]="+pathLength[newIndex]+" common: pathLength["+commonIndex+"]="+pathLength[commonIndex]);

   	   								if (cycleWithBestArea>minCycleArea){ // for most (no-node) pixels will fail
   	   	   	   			    			if (debugLevel>2) System.out.print( " cycleWithBestArea> "+minCycleArea+", old parents="+String.format("0x%02x ", parents[newIndex]));
   										int pDir=0;
   										for (int i=0;i<8;i++) if ((parents[newIndex] &(1<<i))!=0) {
   											pDir=i;
   											break;
   										}
   	   									parents[newIndex]|=parentMask; // add new parent - afre pDir is found
   	   	   	   			    			if (debugLevel>2) System.out.println( String.format(" new mask =0x%02x new parents[%d]=0x%02x ", parentMask, newIndex, parents[newIndex]));
   										node[newIndex]=new double [9][3];
   	   	   	   			    			if (debugLevel>2) System.out.println( String.format("Saved node[%d][%d] and node[%d][8][0]=%d", newIndex, pDir,newIndex, pDir));
   										// Save the only data set to the parents data set
   										node[newIndex][pDir][0]=timeToReach[newIndex];
   										node[newIndex][pDir][1]=pathArea[newIndex];
   										node[newIndex][pDir][2]=pathLength[newIndex];
   										node[newIndex][8][0]=pDir; // direction to the best (old) parent

   										if (debugLevel>2) System.out.println( String.format("Saved node[%d][%d]", newIndex, dirToParent));

   										// Save new parent data
   	   									node[newIndex][dirToParent][0]=newValue;
   	   									node[newIndex][dirToParent][1]=newArea;
   	   									node[newIndex][dirToParent][2]=newLength;
   	   								} // if (cycleWithBestArea>minCycleArea) - do not do anything

   								} else { // if (node[newIndex]==null){ - adding more to the same node - is that likely?
		   	   			    		if (debugLevel>2) System.out.print( " Node existed, adding data to direction "+dirToParent);

   									parents[newIndex]|=parentMask; // add new parent (we may remove it later
   									node[newIndex][dirToParent][0]=newValue;
   									node[newIndex][dirToParent][1]=newArea;
  									node[newIndex][dirToParent][2]=newLength;
									if (debugLevel>2) System.out.println( String.format("Saved node[%d][%d] , now parents[%d]=0x%02x", newIndex, dirToParent, newIndex,parents[newIndex]));

   									for (int pDir=0;pDir<8;pDir++) if ((pDir!=dirToParent) && ((parents[newIndex] & (1<<pDir))!=0)){
   										//int parentIndex=topology[newIndex][pDir];
   										int cycleWithParentArea=newArea-(int) node[newIndex][pDir][1];
   										if (cycleWithParentArea<0) cycleWithParentArea=-cycleWithParentArea;
   		   								int commonIndex=commonAncestor(
   		   										waveIndex,
   		   						    			newIndex,
   		   						    			topology,
   		   						    			timeToReach,
   		   						    			parents,
   		   						    			node,// to find the best parent if there are many of them
   		   						    			debugLevel);

   										cycleWithParentArea-=(newLength+ (int) node[newIndex][pDir][2]-2*(1+pathLength[commonIndex]));// should be 0 if the paths have no gaps - verify?
   		   	   			    			if (debugLevel>2) System.out.println( " Comparing paths in the node with all parents, cycleWithBestArea="+cycleWithParentArea+" newArea="+newArea+
   		   	   			    					" node["+newIndex+"]["+pDir+"][1]="+node[newIndex][pDir][1]+
   		    	   			    				" newLength="+newLength+" node["+newIndex+"]["+pDir+"][2]="+node[newIndex][pDir][2]+" common: pathLength["+commonIndex+"]="+pathLength[commonIndex]);


   										if (cycleWithParentArea<minCycleArea){ // remove that parent
   											if (node[newIndex][pDir][0]<newValue) { // remove new one
   	   											parents[newIndex] &=~parentMask;
   	   											break;
   											} else { // remove old one (that was worse).
   												parents[newIndex] &= ~(1<<pDir); // do not break - may happen that multiple will be removed ??
   											}
   										}
   									}
   								}
   							}
   						}

   					} //if (newIndex>=0)
   				}
   			}
   			if (showProgress)IJ.showProgress(1.0);
// May be a good idea to display wavefront image here for debugging?
//  			now we can remove extra pixels
   			if (dbgImage!=null){ //   			float [][]dbgImage=new float [7][image.length];
   				for (int n=0;n< dbgImage.length;n++) for (int i=0;i<image.length;i++){
   					dbgImage[n][i]=Float.NaN;
   				}
   				for (int i=0;i<topologySize;i++){
   					int index=pixelX[i]+pixelY[i]*imageWidth;
   					dbgImage[0][index]=(float) speed[i]; // same as source image, just masked
   					dbgImage[1][index]=(float) timeToReach[i]; //
   					dbgImage[2][index]=pathArea[i];    //
   					dbgImage[3][index]=pathLength[i];  //
   					if (node[i]!=null){
   						int numParents=0;
   						for (int dir=0;dir<8;dir++) if ((parents[i]& (1<<dir))!=0) numParents++;
   						dbgImage[4][index]=numParents;
   					}
   					dbgImage[6][index]=parents[i];  //


   				}

   				if (debugLevel>1) {
   					System.out.println("Topology:");
   					for (int i=0;i<1000;i++) {
   						System.out.print(String.format("%03d:" , i));

   						for (int j=0;j<8;j++)   System.out.print(String.format(" %05d" , topology[i][j]));
   						System.out.println(String.format(" X= %04d Y=%04d timeToReach=%f parents=0x%02x pathArea=%d pathLength=%d",
   								pixelX[i], pixelY[i],timeToReach[i],parents[i],pathArea[i],pathLength[i]));
   					}
   				}

   			}
// now trim - first pass
   			waveList=new ArrayList<Integer>(10000);
   			ArrayList<Integer>endsList=new ArrayList<Integer>(10000);
   			byte [] stageBytes=new byte[topologySize];
   			for (int i=0;i<topologySize;i++) stageBytes[i]=0;
   			for (int index=0;index<topologySize;index++){
   				if (parents[index]==0) continue; // already removed (not here)
   				int numParents=0;
   				for (int dir=0;dir<8;dir++) if ((parents[index]& (1<<dir))!=0) numParents++;
   				if ((debugLevel>1) &&(numParents>1)) {
   					System.out.println(String.format("Considering cell %05d with %d parents, parents[%05d]=0x%02x", index, numParents, index,parents[index]));
//   					continue; // multiple parents - keep it
   				}
   				boolean hasGreater=false;
   				boolean hasChildren=false;
   				for (int dir=0;dir<8;dir++){
   					int childIndex=topology[index][dir];
   					if (childIndex>=0){
   						int dirToParent=(dir+4)%8;
   						if ((parents[childIndex]& (1<<dirToParent))!=0) hasChildren=true;
   						if (timeToReach[childIndex]>timeToReach[index]) hasGreater=true;
   						if (hasGreater && hasChildren) break;
   					}
   				}
   				if (!hasGreater && !hasChildren){ // try farther
   					endsList.clear();
   					endsList.add(new Integer(index));
   					int listPos=0;
   					for (int step=0; !hasGreater && (step<endRadius);step++){
   						int listEnd=endsList.size();
   						for (int i=listPos;!hasGreater && (i<listEnd);i++){
   							Integer Index=endsList.get(i);
   							for (int dir=0;!hasGreater && (dir<8);dir++){
   								Integer childIndex=topology[Index][dir];
   								if ((childIndex>=0) && (stageBytes[childIndex]==0)){
   									if (timeToReach[childIndex]>timeToReach[index]){
   										hasGreater=true;
   										break;
   									} else {
   										endsList.add(childIndex);
   										stageBytes[childIndex]=1;
   									}

   								}
   							}

   						}

   					}
   					while (endsList.size()>0){
   						Integer Index=endsList.remove(0);
   						stageBytes[Index]=0;
   					}
   				}
   				if (hasGreater && !hasChildren){ // trim
   					waveList.add(new Integer(index));
   					if (debugLevel>1) System.out.println("Adding "+index+" for removal on the initial pass X="+pixelX[index]+", Y="+pixelY[index]+" numParents="+numParents);

   					if (parents[index]==0){
   						System.out.println("BUG at X="+pixelX[index]+", Y="+pixelY[index]+" - parents["+index+"]=0");
   					}
   				}
   			}
   			stageBytes=null;
   			while (waveList.size()>0){
   				int index=waveList.remove(0);
// multiple parents is not an excuse for survival
   				int parentIndex=-1;
   				if (parents[index]==0){
   					System.out.println("2: BUG at X="+pixelX[index]+", Y="+pixelY[index]+" - parents["+index+"]=0");
   				}
   				int hadParents=parents[index];
   				parents[index]=0;
   				for (int pDir=0;pDir<8;pDir++) if ((hadParents & (1<<pDir))!=0){ // should be just one - not true anymore, need to process all parents
   					parentIndex=topology[index][pDir];
   					//  				}
   					if (parentIndex<0){
   						System.out.println("BUG at X="+pixelX[index]+", Y="+pixelY[index]+" - parents["+index+"]="+hadParents+" parentIndex="+parentIndex);
   					}
   					boolean hasGreater=false;
   					boolean hasChildren=false;
   					for (int dir=0;dir<8;dir++){
   						int childIndex=topology[parentIndex][dir]; //-1
   						if (childIndex>=0) {
   							int dirToParent=(dir+4)%8;
   							if ((parents[childIndex]& (1<<dirToParent))!=0) hasChildren=true;
   							if (timeToReach[childIndex]>timeToReach[parentIndex]) hasGreater=true;
   							if (hasGreater && hasChildren) break;
   						}
   					}
   					if (hasGreater && !hasChildren){ // trim
   						waveList.add(new Integer(parentIndex));
   						if (debugLevel>1) System.out.println(String.format("ADDING for removal %05d X=%05d, Y=%05d child-X=%05d, child-Y=%05d parents[%05d]=0x%02x",
   								parentIndex, pixelX[parentIndex],pixelY[parentIndex],pixelX[index],pixelY[index],index, hadParents));
   					} else {
   						if (debugLevel>1) System.out.println(String.format("NOT adding for removal %05d X=%05d, Y=%05d child-X=%05d, child-Y=%05d parents[%05d]=0x%02x",
   								parentIndex, pixelX[parentIndex],pixelY[parentIndex],pixelX[index],pixelY[index],index, hadParents)+
   								" hasGreater="+hasGreater+" hasChildren="+hasChildren);
   					}
   				}

   			}
   			if (dbgImage!=null){ //   			float [][]dbgImage=new float [7][image.length];
   				for (int i=0;i<topologySize;i++){
   					int index=pixelX[i]+pixelY[i]*imageWidth;
   					dbgImage[5][index]=(parents[i]!=0)?1.0f:0.0f;
   				}
   			}
   			if (bitmapAccumulator!=null){
   				for (int i=0;i<topologySize;i++){
   					if (parents[i]!=0) bitmapAccumulator[pixelX[i]+pixelY[i]*imageWidth]=true;;
   				}
   			}
//   			return dbgImage;
    	}

    	public boolean [] EdgeThinning(
    			double [] image,
    			boolean [] edges,
    			int imageWidth,
    			int minCycleArea, // minimal area of the closed loop to keep
    			int endRadius,
    			boolean showProgress,
    			int debugLevel
    	){
    		int [][]dirs={{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}}; // x,y
    		int imageHeight=image.length/imageWidth;
    		boolean [] edgesStage=edges.clone();
    		boolean []bitmapAccumulator=new boolean [image.length];
    		for (int i=0;i<image.length;i++) bitmapAccumulator[i]=false;
			ArrayList<Integer> edgeList=new ArrayList<Integer>(10000);
			int numEdge=0;
    		for (int startIndex=0;startIndex<image.length;) {
    			if (showProgress) IJ.showStatus("Thinning edge "+(++numEdge)+"...");
    			while ((startIndex<image.length) && !edgesStage[startIndex]) startIndex++;
    			Integer Index=startIndex;
    			if (startIndex<image.length) {
    				edgeList.clear();
    				edgeList.add(Index);
    				edgesStage[Index]=false;
    				int position=0;
    	   			while (position<edgeList.size()) {
    	   				Index=edgeList.get(position++);
    	   				int x=Index%imageWidth;
    	   				int y=Index/imageWidth;
    	   				for (int dir=0;dir<dirs.length;dir++){
    	   					int x1=x+dirs[dir][0];
    	   					int y1=y+dirs[dir][1];
    	   					if ((x1>=0) && (x1<imageWidth) && (y1>=0) && (y1<imageHeight)){
    	   						Index=x1+y1*imageWidth;
    	   						if (edgesStage[Index]){
    	   							edgesStage[Index]=false;
    	   							edgeList.add(Index);
    	   						}
    	   					}
    	   				}
    	   			}
    	   			thinEdge(
    	   	    			edgeList,
    	   	    			image,
    	   	    			bitmapAccumulator,
    	   	    			imageWidth,
    	   	    			minCycleArea, // minimal area of the closed loop to keep
    	   	    			endRadius,
    	   	    			//TODO: - other parameters to trim short branches and detect ends?
    	   	    			showProgress,
    	   	    			null, //dbgImage,
    	   	    			debugLevel
    	   	    			);
    			} // if (startIndex<image.length)
    		}
			if (showProgress) IJ.showStatus("");
    		return bitmapAccumulator;
    	}



    	public float [][] testEdgeThinning(
    			double [] image,
    			boolean [] edges,
    			int imageWidth,
    			int startX, // start looking for the edge to try
    			int startY,
    			int minCycleArea, // minimal area of the closed loop to keep
    			int endRadius,
    			boolean showProgress,
    			//TODO: - other parameters to trim short branches and detect ends?
    			int debugLevel
    	){
    		int [][]dirs={{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}}; // x,y
    		int imageHeight=image.length/imageWidth;
    		boolean [] edgesStage=edges.clone();
    		Integer Index=-1;
    		for (int i=startX+imageWidth*startY;i<edgesStage.length;i++) if (edgesStage[i]){
    			Index=i;
    			break;
    		}
    		if (Index<0){
    			return null;
    		}
			ArrayList<Integer> edgeList=new ArrayList<Integer>(10000);
			edgeList.add(Index);
			edgesStage[Index]=false;
			int position=0;
   			while (position<edgeList.size()) {
   				Index=edgeList.get(position++);
   				int x=Index%imageWidth;
   				int y=Index/imageWidth;
   				for (int dir=0;dir<dirs.length;dir++){
   					int x1=x+dirs[dir][0];
   					int y1=y+dirs[dir][1];
   					if ((x1>=0) && (x1<imageWidth) && (y1>=0) && (y1<imageHeight)){
   						Index=x1+y1*imageWidth;
   						if (edgesStage[Index]){
   							edgesStage[Index]=false;
   							edgeList.add(Index);
   						}
   					}
   				}
   			}
   			if (debugLevel>0){
   				System.out.println("testEdgeThinning(): list size="+edgeList.size());
   			}
   			float [][]dbgImage=new float [7][image.length];
   			thinEdge(
   	    			edgeList,
   	    			image,
   	    			null,
   	    			imageWidth,
   	    			minCycleArea, // minimal area of the closed loop to keep
   	    			endRadius,
   	    			//TODO: - other parameters to trim short branches and detect ends?
   	    			showProgress,
   	    			dbgImage,
   	    			debugLevel
   	    			);
   			return dbgImage;
    	}
    	public class EdgesSegmeniting{
    		public byte emptyCell=0;
    		public byte badCell=1;
    		public byte oldCell=2;
    		public byte newCell=3;
    		public int imageWidth;
    		public int imageHeight;
    		public int imageLength;


    		public int [] iDirs4=null; // orthogonal only
    		public int [] iDirs8=null; // diagonoa
    		public int [][]dirs4={{1,0},      {0,1},       {-1,0},        {0,-1}}; // x,y
    		public int [][]dirs8={{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}}; // x,y

    		public EdgesSegmeniting(){}
    		public EdgesSegmeniting(int width,int height){setImageDimensions(width,height);}
    		public void setImageDimensions(
    				int width,
    				int height){
    			this.imageWidth=width;
    			this.imageHeight=height;
    			this.imageLength=width*height;
    			int [] iDirs4={1,                    this.imageWidth,                    -1,                     -this.imageWidth                    }; // orthogonal only
    			int [] iDirs8={1, this.imageWidth+1, this.imageWidth, this.imageWidth-1, -1, -this.imageWidth-1, -this.imageWidth, -this.imageWidth+1}; // diagonal
    			this.iDirs4=iDirs4;
    			this.iDirs8=iDirs8;
    		}

    		public float [] distanceFromEdges(
    				boolean [] edges,
    				boolean [] alphaMask,// or null
    				int imageWidth,
    				int numPasses,
    				int debugDirections,
    				int debugLevel
    		){
    			setImageDimensions(imageWidth,edges.length/imageWidth);
    			return distanceFromEdges(
    					edges,
    					alphaMask,
    					numPasses,
    					debugDirections,
    					debugLevel
    			);
    		}

    		public float [] distanceFromEdges(
    				boolean [] pureEdges,
    				boolean [] alphaMask,// or null
    				int numPasses,
    				int debugDirections,
    				int debugLevel
    		){
    			boolean [] edges=pureEdges;
    			if (alphaMask!=null){
    				edges=pureEdges.clone();
    				for (int i=0;i<edges.length;i++) edges[i]|= !alphaMask[i];
    			}
    			if (debugDirections==0) debugDirections=0xfffff;
    			float [] stage = new float [edges.length];
    			float maxDistance=(float) Math.sqrt(this.imageWidth*this.imageWidth+this.imageHeight*this.imageHeight);
    			for (int i=0;i<stage.length;i++) stage[i]=maxDistance;
    			for (int pass=0;pass<numPasses;pass++){
    				for (int y=0;y<this.imageHeight;y++) {
    					if ((debugDirections & 1)!=0) distanceFromEdgesLine(
    							edges,
    							stage, // should be initialized with large enough number
    							0, //int x0,
    							y, //int y0,
    							1, //int dx, // 0/=1
    							0, //int dy, // -1, 0, +1
    							debugLevel);
    					if ((debugDirections & 2)!=0) distanceFromEdgesLine(
    							edges,
    							stage, // should be initialized with large enough number
    							0, //int x0,
    							y, //int y0,
    							1, //int dx, // 0/=1
    							1, //int dy, // -1, 0, +1
    							debugLevel);

    					if ((debugDirections & 4)!=0) distanceFromEdgesLine(
    							edges,
    							stage, // should be initialized with large enough number
    							0, //int x0,
    							y, //int y0,
    							1, //int dx, // 0/=1
    							-1, //int dy, // -1, 0, +1
    							debugLevel);
    				}
    				for (int x=0;x<this.imageWidth;x++) {
    					if ((debugDirections & 8)!=0) distanceFromEdgesLine(
    							edges,
    							stage, // should be initialized with large enough number
    							x, //int x0,
    							0, //int y0,
    							0, //int dx, // 0/=1
    							1, //int dy, // -1, 0, +1
    							debugLevel
    					);
    					if (x>0){
    						if ((debugDirections & 16)!=0) distanceFromEdgesLine(
    								edges,
    								stage, // should be initialized with large enough number
    								x, //int x0,
    								0, //int y0,
    								1, //int dx, // 0/=1
    								1, //int dy, // -1, 0, +1
    								debugLevel);
    						if ((debugDirections & 32)!=0) distanceFromEdgesLine(
    								edges,
    								stage, // should be initialized with large enough number
    								x, //int x0,
    								this.imageHeight-1, //int y0,
    								1, //int dx, // 0/=1
    								-1, //int dy, // -1, 0, +1
    								debugLevel);
    					}
    				}
    			}
    			return stage;
    		}

    		private void distanceFromEdgesLine(
    				boolean [] edges,
    				float [] stage, // should be initialized with large enough number
    				int x0,
    				int y0,
    				int dx, // 0/=1
    				int dy, // -1, 0, +1
    				int debugLevel
    		){
    			if (((dx==0) && (dy==0)) || (dx<-1)  || (dy<-1)  || (dx>1)  || (dy>1)) {
    				String msg="Illegal values for dx="+dx+" and/or dy "+dy;
    				IJ.showMessage("Error",msg);
    				throw new IllegalArgumentException (msg);
    			}
    			if ((x0<0) || (y0<0) || (x0>=this.imageWidth) || (y0>=this.imageHeight)){
    				String msg=String.format("Illegal values for x0=%d (should be 0 <= x0 < %d) and/or y0=%d (should be 0 <= y0 <%d ",x0, this.imageWidth, y0, this.imageHeight);
    				IJ.showMessage("Error",msg);
    				throw new IllegalArgumentException (msg);
    			}
    			double stepDistance=Math.sqrt(dx*dx+dy*dy);
    			int length;
    			if (dx!=0){
    				if (dx>0) length= this.imageWidth-x0;
    				else length=x0;
    			} else length=this.imageHeight;
    			if (dy!=0){
    				if (dy>0) {
    					if ((this.imageHeight-y0)<length) length=(this.imageHeight-y0);
    				} else {
    					if ((y0)<length) length=y0;
    				}
    			}
    			if (length<1) return;
    			// direct pass;
    			int xStep=dx;
    			int yStep=this.imageWidth*dy;
    			int indexStep=xStep+yStep;
    			boolean isDiagonal=(dx!=0) && (dy!=0);
    			int index=x0+this.imageWidth*y0;
    			double edgeDistance=stage[index];
    			for (int i=0;i<length;i++){
    				if (edges[index]) edgeDistance=0.0;
    				else if ((i>0) && isDiagonal) {
    					if (edges[index-xStep] || edges[index-yStep])edgeDistance=0.5*stepDistance;
    				}
    				if (stage[index]>edgeDistance) stage[index]=(float) edgeDistance;
    				else edgeDistance=stage[index];
    				edgeDistance+=stepDistance;
    				index+=indexStep;
    			}
    			// reverse pass;
    			indexStep=-indexStep;
    			xStep=-xStep;
    			yStep=-yStep;
    			index+=indexStep; // just one step back
    			edgeDistance=stage[index];
    			for (int i=0;i<length;i++){
    				if (edges[index]) edgeDistance=0.0;
    				else if ((i>0) && isDiagonal) {
    					if (edges[index-xStep] || edges[index-yStep])edgeDistance=0.5*stepDistance;
    				}
    				if (stage[index]>edgeDistance) stage[index]=(float) edgeDistance;
    				else edgeDistance=stage[index];
    				edgeDistance+=stepDistance;
    				index+=indexStep;
    			}
    		}

    		public int [] extractEdgesAreas(
    				boolean gapsNotEdges,
    				double threshold,
    				float [] edgeDistance,
    				int imageWidth,
    				Rectangle [][] bounds0,  // null or one element Rectangle[] - will return bounding rectangle for each of the edge areas
    				boolean showProgress,
    				int debugLevel
    		) {
    			setImageDimensions(imageWidth,edgeDistance.length/imageWidth);
    			return extractEdgesAreas(
    					gapsNotEdges,
    					threshold,
    					edgeDistance,
    					bounds0,  // null or one element Rectangle[] - will return bounding rectangle for each of the edge areas
    					showProgress,
    					debugLevel);
    		}
    		public int [] extractEdgesAreas(
    				boolean gapsNotEdges,
    				double threshold,
    				float [] edgeDistance,
    				Rectangle [][] bounds0,  // null or one element Rectangle[] - will return bounding rectangle for each of the edge areas
    				boolean showProgress,
    				int debugLevel
    		){
    			int [] edgesAreas=new int [this.imageLength];
    			for (int i=0;i<this.imageLength;i++) edgesAreas[i]=0;
    			ArrayList<Integer> edgeList=new ArrayList<Integer>(10000);
    			int numEdge=0;
    			int xMin=2,xMax=this.imageWidth-3,yMin=2,yMax=this.imageHeight-3; // to avoid checking when following paths around (during edges vacuum)
    			for (int startIndex=0;startIndex<this.imageLength;) {
    				numEdge++;
    				if (showProgress) IJ.showStatus("Isolating edge group "+(numEdge)+"...");
    				if (debugLevel>1)   System.out.println("Isolating edge group "+(numEdge)+"...");
    				while ((startIndex<this.imageLength) &&
    						((edgesAreas[startIndex]!=0) ||
    								(gapsNotEdges?(edgeDistance[startIndex]<=threshold):(edgeDistance[startIndex]>threshold)) ||
    								((startIndex%this.imageWidth)<xMin) ||
    								((startIndex%this.imageWidth)>xMax) ||
    								((startIndex/this.imageWidth)<yMin) ||
    								((startIndex/this.imageWidth)>yMax) )) startIndex++;
    				Integer Index=startIndex;
    				if (debugLevel>1)System.out.println("Start point: x="+(startIndex%this.imageWidth)+", y="+(startIndex/this.imageWidth)+
    						" (xMin="+xMin+" xMax="+xMax+" (yMin="+yMin+" yMax="+yMax+")");
    				if (startIndex<this.imageLength) {
    					edgeList.clear();
    					edgeList.add(Index);
    					edgesAreas[Index]=numEdge;
//    					int position=0;
    					int dbgLength=0;
    					while (edgeList.size()>0) {
    						Index=edgeList.remove(0);
    						int x=Index%this.imageWidth;
    						int y=Index/this.imageWidth;
    						for (int dir=0;dir<this.dirs8.length;dir++){
    							int x1=x+this.dirs8[dir][0];
    							int y1=y+this.dirs8[dir][1];
    							if ((x1>=xMin) && (x1<=xMax) && (y1>=yMin) && (y1<=yMax)){
    								Index=x1+y1*this.imageWidth;
    								if ( (edgesAreas[Index]==0) && (gapsNotEdges?(edgeDistance[Index]>threshold):(edgeDistance[Index]<=threshold))){
    									edgesAreas[Index]=numEdge;
    									edgeList.add(Index);
    									dbgLength++;
    								}
    							}
    						}
    					}
        				if (debugLevel>1)System.out.println("Area size="+dbgLength);
    				}
    			}
    			if (bounds0!=null){
    				Rectangle [] bounds=new Rectangle[numEdge-1];
    				for (int i=0;i<bounds.length;i++) bounds[i]=new Rectangle(this.imageWidth, this.imageHeight,0,0);
    				for (int i=0;i<this.imageLength;i++)if (edgesAreas[i]>0){
    					int d= edgesAreas[i]-1;
    					int x=i%this.imageWidth;
    					int y=i/this.imageWidth;
    					if (x< bounds[d].x) bounds[d].x=x;
    					if (y< bounds[d].y) bounds[d].y=y;
    					if (x> bounds[d].width) bounds[d].width=x;
    					if (y> bounds[d].height) bounds[d].height=y;
    				}
    				for (int i=0;i<bounds.length;i++) {
    					bounds[i].width-= bounds[i].x-1;
    					bounds[i].height-=bounds[i].y-1;
    					if (debugLevel>1)System.out.println("bounds["+i+"]: x="+bounds[i].x+", y="+bounds[i].y+", width="+bounds[i].width+", height="+bounds[i].height);
    				}
    				bounds0[0]=bounds;
    			}
    			if (showProgress) IJ.showStatus("");
    			return edgesAreas;
    		}

    		public ArrayList<ArrayList<Integer>> createInitialEdgeAreaBorder(
    				int [] edgeAreas,
    				int edgeAreaNumber,
    				int imageWidth,
    				Rectangle bounds, // bounding rectangle to speed up search. Should not have y, x<2 and > appropriate dimensions
    				boolean showProgress,
    				int debugLevel
    		){
    			setImageDimensions(imageWidth,edgeAreas.length/imageWidth);
    			return createInitialEdgeAreaBorder(
    					edgeAreas,
    					edgeAreaNumber,
    					bounds, // bounding rectangle to speed up search. Should not have y, x<2 and > appropriate dimensions
    					showProgress,
    					debugLevel);
    		}

    		public ArrayList<ArrayList<Integer>> createInitialEdgeAreaBorder( // never returns
    				int [] edgeAreas,
    				int edgeAreaNumber,
    				Rectangle bounds, // bounding rectangle to speed up search. Should not have y, x<2 and > appropriate dimensions
    				boolean showProgress,
    				int debugLevel
    		){
    			int minDebugLevel=3;
    			int length=edgeAreas.length;
    			boolean [] stage=new boolean[length];
    			int margins=2;

    			int yMin=bounds.y-margins, yMax=bounds.y+bounds.height+margins,xMin=bounds.x-margins, xMax=bounds.x+bounds.width+margins;

    			for (int y=yMin;y<yMax;y++)	for (int x=xMin;x<xMax;x++){
    				//    			int index=x+y*imageWidth;
    				//    			stage[x+y*imageWidth]=edgeAreas[x+y*imageWidth]==edgeAreaNumber;
    				stage[x+y*this.imageWidth]=false; // just clear
    			}

    			ArrayList<ArrayList<Integer>> boundariesList=new ArrayList<ArrayList<Integer>>(100);
    			for (int startIndex=0;startIndex<(length-1);) {
    				//    			if (showProgress) IJ.showStatus("Isolating edge group "+(numEdge)+"...");
    				while ((startIndex<(length-1)) &&
    						(stage[startIndex] || // boundary already marked
    								(edgeAreas[startIndex+1]!=edgeAreaNumber) || // next does not fit
    								(edgeAreas[startIndex]==edgeAreaNumber))) startIndex++; // this fits
    				if (startIndex>=(length-1)) break;
    				// Now this current does not belong to edgeArea, and the east neighbor - does
    				int dir=3; // up
    				ArrayList<Integer> boundarysList=new ArrayList<Integer>(10000);
    				Integer Index=startIndex;
    				int startDir=dir;

					dir=(startDir+1)%4;
					if ((edgeAreas[Index+this.iDirs4[dir]]!=edgeAreaNumber)) startDir=dir;
					else if ((edgeAreas[Index+this.iDirs4[startDir]]==edgeAreaNumber)){
						dir=(startDir+3)%4;
						if ((edgeAreas[Index+this.iDirs4[dir]]!=edgeAreaNumber)) startDir=dir;
						else startDir=(startDir+2)%4; // turn back
					} // nothing matched - go on forward
    				dir=startDir;
    				if (debugLevel>=minDebugLevel) System.out.println("createInitialEdgeAreaBorder(,"+edgeAreaNumber+",..): startIndex="+startIndex+" startDir="+startDir);
    				int debugNum=0;
    				do {
    					boundarysList.add(Index);
    					stage[Index]=true;
    					Index+=this.iDirs4[dir];
    					if (edgeAreas[Index]==edgeAreaNumber) break; // one-point border
    					int dir1=(dir+1)%4;
    					if ((edgeAreas[Index+this.iDirs4[dir1]]!=edgeAreaNumber)) dir=dir1;
    					else if ((edgeAreas[Index+this.iDirs4[dir]]==edgeAreaNumber)){
    						dir1=(dir+3)%4;
    						if ((edgeAreas[Index+this.iDirs4[dir1]]!=edgeAreaNumber)) dir=dir1;
    						else dir=(dir+2)%4; // turn back
    					} // nothing matched - go on forward
        				if (debugLevel>=minDebugLevel){
        					System.out.println("createInitialEdgeAreaBorder(): Index="+Index+", dir="+dir+" x="+(Index%this.imageWidth)+" y=" +(Index/this.imageWidth)+
        							" (startIndex="+startIndex+" startDir="+startDir+" strat-x="+(startIndex%this.imageWidth)+" start-y=" +(startIndex/this.imageWidth)+")");
        				}
        				if (Index==startIndex) debugNum++;
        				if (debugNum>10){
        	        		String msg="Going circles 10 times";
        	        		IJ.showMessage("Error",msg);
        	        		throw new IllegalArgumentException (msg);
        				}
    				} while ((Index!=startIndex) || (dir!=startDir));

    				// follow border, mark boundary on stage
    				boundariesList.add(boundarysList);
    			}
    			return boundariesList;
    		}

    		public int vacuumEdgeAreaBorder(
    				double thresholdStart, // the same as for initial segmenting
    				double thresholdFinal, // threshold (lower) for the shrank edge areas
    				double thresholdStep,  // decrease threshold between passes
    				boolean repeatEachStep, // repeat each step until no new points
    				boolean repeatLastStep, // repeat last step until no new points

    				float [] edgeDistance, // array of distances from the nearest edge
    				int [] edgeAreas,   // or null. If not null, will zero out removed areas
    				int edgeAreaNumber, // any if edgeAreas is null
    				ArrayList<ArrayList<Integer>> borders, // will be modified
    				int imageWidth,
    				Rectangle bounds, // bounding rectangle to speed up search. Should not have y, x<2 and > appropriate dimensions
    				boolean showProgress,
    				int debugLevel
    		){
    			//    		byte emptyCell=0;
    			//    		byte badCell=1;
    			//    		byte oldCell=2;
    			//    		byte newCell=3;
    			///		private enum cellTypes {EMPTY_CELL, BAD_CELL, OLD_CELL,NEW_CELL}
    			int minDebugLevel=3;
    			int totalNewCells=0;
				if (debugLevel>=minDebugLevel){
					System.out.println("vacuumEdgeAreaBorder("+thresholdStart+","+thresholdFinal+","+thresholdStep+", [][] ,"+edgeAreaNumber+",<>,"+imageWidth+",,,"+debugLevel+")");
				}


    			int length=edgeAreas.length;
    			//    		int imageHeight=length/imageWidth;
    			byte [] stage=new byte[length];
    			int margins=2;

    			int yMin=bounds.y-margins, yMax=bounds.y+bounds.height+margins,xMin=bounds.x-margins, xMax=bounds.x+bounds.width+margins;
    			for (int y=yMin;y<yMax;y++)	for (int x=xMin;x<xMax;x++){
    				stage[x+y*imageWidth]=this.emptyCell; // just clear
    			}
    			// Mark initial border(s)
    			for (int numBorder=0;numBorder<borders.size();numBorder++) {
    				ArrayList<Integer> border= borders.get(numBorder);
    				for (int i=0;i<border.size();i++) stage[border.get(i)]=this.oldCell;
    			}

				if (debugLevel>=minDebugLevel){
					float [] fEdge=new float [length];
					for (int i=0;i<length;i++) fEdge[i]=stage[i];

					(new ShowDoubleFloatArrays()).showArrays(
							fEdge,
							imageWidth,
							length/imageWidth,
							"Edge="+edgeAreaNumber);
				}



    			if (thresholdFinal>thresholdStart) {
    	    		return -1; // error
    			}
    			thresholdStep=Math.abs(thresholdStep);
    			if (thresholdStep==0){
    				thresholdFinal=thresholdStart;
    				thresholdStep=1.0;
    			}
				if (debugLevel>=minDebugLevel){
					System.out.println("vacuumEdgeAreaBorder() thresholdStep="+thresholdStep+" thresholdStart="+thresholdStart+" thresholdFinal="+thresholdFinal);
				}

    			ArrayList<Integer> badCells= new ArrayList<Integer>(10000);
// TODO - make last pass until no new cells? Adjust new threshold to the maximal value of the new added cell
    			if ((thresholdStart-thresholdStep)>=thresholdFinal) thresholdStart-=thresholdStep; // thresholdStart was used fro initial border
    			for (double threshold=thresholdStart; threshold>=thresholdFinal;threshold-=thresholdStep){ // threshold is exclusive - allowedDistance should be > threshold
    				int thisNewCells=0;
    				for (int numBorder=0;numBorder<borders.size();numBorder++) {
    					ArrayList<Integer> border= borders.get(numBorder);
    					if (border.size()<2){
    						System.out.println("vacuumEdgeAreaBorder(): border.size()="+border.size()+"<2 - BUG???");
    						continue;
    					}
    					int index=-1;
    					int index1=-1,index2=border.get(0);
    					int startDir=-1;
    					for (int i=1;i<border.size();i++){
    						index1=index2;
    						index2=border.get(i);
    						startDir=(index2>index1)?((index2>(index1+1)?1:0)):((index2<(index1-1)?3:2));
    						index=index1+this.iDirs4[(startDir+1)%4];
    						if (stage[index]!=this.oldCell) break; // no newCells - only empty or bad
    					}
    					if (index<0){
    						String msg="all right-hand neighbors are OLD";
    						IJ.showMessage("Error",msg);
    						throw new IllegalArgumentException (msg);
    					}

    					// this.badCell - not currently used
    					// find the start point for the new border segment
    					if (stage[index]==this.badCell) {
    						index=index1; // start with old cell - will have to come here again
    					} else { // can only be a new cell here
    						if (edgeDistance[index]<=threshold) {
    							index=index1; // start with old cell - will have to come here again
    						} else {
    							int numNeibGroups=0;
    							for (int dir8=0;dir8<8;dir8++) if ((stage[index+this.iDirs8[dir8]]<this.oldCell) && (stage[index+this.iDirs8[(dir8+1)%8]]>=this.oldCell)) {
    								numNeibGroups++;
    								if (numNeibGroups>1) break; // no need to look more
    							}
    							if (numNeibGroups>1){
    								index=index1; // start with old cell - will have to come here again
    								// no, can not mark it "bad"
    							}
    						}
    					}
    					// find start direction for the new border
    					int dir=startDir;
    					int blockedCell=-1; // none
    					float blockedCellValue=Float.NaN;
    					//starting from a new cell (even as the first move will likely be to the old cell) ensures the border will shrink here.
    					// temoporarily blocking the blockedCell prevents "short-circuit" near the start cell.
    					// using edgeDistance[] - no new comparisons will be needed
    					if (index!=index1) {
    						for (dir=startDir;dir!=((startDir+3)%4); dir=(dir+1)%4) if (stage[index+this.iDirs4[dir]]!=this.oldCell) break;
    						if (dir==((startDir+3)%4)){
    							String msg="all old cells around index="+index+ " (x="+(index%this.imageWidth)+", y="+index/this.imageWidth+")";
    							IJ.showMessage("Error",msg);
    							throw new IllegalArgumentException (msg);
    						}
    						startDir=(dir+3)%4; // last old cell clockwise from the old start
    						blockedCell=index+this.iDirs4[dir];
    						blockedCellValue=edgeDistance[blockedCell];
    						edgeDistance[blockedCell]=(float)(threshold-1.0);
    						// temporarily block the first non-old cell by lovering threshold
    					}


    					// now we do not need old border segment list, can start again:
    					int startIndex=index;
    					border.clear();
    					badCells.clear();
    					Integer Index=startIndex;
    					dir=startDir;

    					if (debugLevel>=minDebugLevel){
    						System.out.println(" **** vacuumEdgeAreaBorder() startDir="+startDir+" startIndex="+startIndex+
    								" start-x="+(startIndex%imageWidth)+" start-y="+(startIndex/imageWidth)+
    								" blockedCell="+blockedCell+ " threshold="+threshold);
    					}

						if (stage[Index]!=this.oldCell) stage[Index]=this.newCell; //*** java.lang.ArrayIndexOutOfBoundsException: -612
    					do {
        					if (debugLevel>=minDebugLevel){
        						System.out.print(" --- dir="+dir+" Index="+Index+
        								" x="+(Index%imageWidth)+" y="+(Index/imageWidth)+
        								" stage[Index]="+stage[Index]+
        								" startDir="+startDir+" startIndex="+startIndex+
        								" start-x="+(startIndex%imageWidth)+" start-y="+(startIndex/imageWidth)+
        								" blockedCell="+blockedCell);
        					}
    						border.add(Index);
//    						if (stage[Index]!=this.oldCell) stage[Index]=this.newCell; //*** java.lang.ArrayIndexOutOfBoundsException: -612
    						Index+=this.iDirs4[dir];//
    						if (stage[Index]!=this.oldCell) stage[Index]=this.newCell; //*** java.lang.ArrayIndexOutOfBoundsException: -612

        					if (debugLevel>=minDebugLevel){
        						System.out.println(" === NEW  Index="+Index+
        								" x="+(Index%imageWidth)+" y="+(Index/imageWidth)+
        								" stage[Index]="+stage[Index]);
        					}

    						int dir1=(dir+1)%4;
    						for (dir1=(dir+1)%4;dir1!=(dir+2)%4; dir1=(dir1+3)%4){
    							// see if step in dir1 is permitted
    							index2=Index+this.iDirs4[dir1];
    							if (stage[index2]>=this.oldCell) break; // old /new cell - OK, use it
    							if ((stage[index2]==this.badCell) || // not yet used
    									(edgeDistance[index2]<=threshold)) continue; // does not fit threshold criteria - try next
    							// index2: still empty cell, fits criteria - make sure it doeas not break the link

    							int numNeibGroups=0;
    							for (int dir8=0;dir8<8;dir8++) if ((stage[index2+this.iDirs8[dir8]]<this.oldCell) && (stage[index2+this.iDirs8[(dir8+1)%8]]>=this.oldCell)) {
    								if (debugLevel>=minDebugLevel) System.out.print(" "+dir8+":"+((dir8+1)%8)+" ");
    								numNeibGroups++;
    								if (numNeibGroups>1) break; // no need to look more
    							}
            					if (debugLevel>=minDebugLevel){
            						System.out.println("\nnumNeibGroups="+numNeibGroups+" index2="+index2+
        								" x="+(index2%imageWidth)+" y="+(index2/imageWidth)+" stage[index2]="+stage[index2]+" edgeDistance[index2]="+edgeDistance[index2]);
            					}


//    							if (numNeibGroups>1) continue; // using this cell will break link between edges, do not use it.
    							if (numNeibGroups<2) break; // using this cell will not break link between edges, OK to use it.
    							badCells.add(new Integer(index2));
    							stage[index2]=this.badCell; // do niot try it again in this pass
    						}
    						// now index2 - where to go, dir1-direction to go
    						dir=dir1;
    					} while ((Index!=startIndex) || (dir!=startDir));
    					// restore original value of the blocked cell
    					if (blockedCell>=0) edgeDistance[blockedCell]=blockedCellValue;
    					// mark new cells as old
    					// restore bad cells
    					int newCells=0;
    					for (int i=0;i<badCells.size();i++) {
    						stage[badCells.get(i)]=this.emptyCell;
    					}
    					// count new cells, make them old
    					for (int i=0;i<border.size();i++) {
    						Index=border.get(i);
    						if (edgeAreas!=null) edgeAreas[Index]=0;
    						if (stage[Index]==this.newCell) {
    							newCells++;
        						stage[Index]=this.oldCell;
    						}
    					}

    					totalNewCells+=newCells;
    					thisNewCells+=newCells;
    					if (debugLevel>=minDebugLevel){
    						System.out.println("++++++++++++ New cells:"+newCells); // repeat untill no new cells?
    					}


    				}
    				if (thisNewCells>0){
    					if (repeatEachStep || (repeatLastStep && ((threshold-thresholdStep) < thresholdFinal))) threshold+=thresholdStep; // will repeat the pass
    				}
    			}
    			if (edgeAreas!=null){
    				for (int y=bounds.y; y<(bounds.y+bounds.height);y++) for (int x=bounds.x; x<(bounds.x+bounds.width);x++){
    					int index=y*this.imageWidth+x;
    					if ((edgeAreas[index]==edgeAreaNumber) && (stage[index]==this.oldCell)) edgeAreas[index]=0; // erase it here
    				}
    			}


        		return totalNewCells;
    		}
    	}



    	public double [] sobelFilter(
    			double [] srcImage,
    			int width,
    			double edgeSpreadOrtho, // not used
    			double edgeSpreadDiagonal // not used
    	){
    		int height=srcImage.length/width;
    		double p1,p2,p3,p4,p5,p6,p7,p8,p9,s1,s2;
    		double [] sobelResult=new double [srcImage.length];
    		int upRight,downRight;
    		int index;
    		for (int iY=0;iY<height;iY++){
    			index=iY*width;
    			p5=srcImage[index];p6=srcImage[index+1];
    			if (iY==0){
    				p2=p5;p3=p6;
            		upRight=  1;
    			} else {
    				p2=srcImage[index-width];p3=srcImage[index-width+1];
            		upRight=  -width+1;
    			}
    			if (iY==(height-1)){
    				p8=p5;p9=p6;
            		downRight= 1;
    			} else {
    				p8=srcImage[index+width];p9=srcImage[index+width+1];
            		downRight= width+1;
    			}
    			p1=p2;p4=p5;p7=p8;

    			for (int iX=0;iX<width;iX++){
    				s1=p1+2*p2+p3-p7-2*p8-p9;
    				s2=p1+2*p4+p7-p3-2*p6-p9;
    				sobelResult[index++]=Math.sqrt(s1*s1+s2*s2);
    				p1=p2;p2=p3;
    				p4=p5;p5=p6;
    				p7=p8;p8=p9;
    				if (iX<(width-2)){ // index already incremented, but iX - not yet
    					p3=srcImage[index+upRight];
    					p6=srcImage[index+1];
    					p9=srcImage[index+downRight];
    				}
    			}
    		}
    		return sobelResult;

    	}

    	public double [] sobelFilter1(
    			double [] srcImage,
    			int width,
    			double edgeSpreadOrtho,
    			double edgeSpreadDiagonal
    	){
    		int height=srcImage.length/width;
    		double p1,p2,p3,p4,p5,p6,p7,p8,p9,s1,s2;
    		double [] sobelResult=new double [srcImage.length];
    		int upRight,downRight;
    		int index;
    		edgeSpreadOrtho=Math.abs(edgeSpreadOrtho);
    		boolean spread=(edgeSpreadOrtho>0.0) || (edgeSpreadDiagonal>0.0);
//    		if (spread)
    			for (int i=0;i<sobelResult.length;i++) sobelResult[i]=0.0;
 //   		System.out.println ("kOrto0="+kOrto0+" kOrto1="+kOrto1+" kDiagonal0="+kDiagonal0+" kDiagonal="+kDiagonal1+" swap="+swap+", spread="+spread);
    		int upLeft,downLeft,down,up;
    		for (int iY=0;iY<height;iY++){
    			index=iY*width;
    			p5=srcImage[index];p6=srcImage[index+1];
    			if (iY==0){
    				p2=p5;p3=p6;
            		upRight=  1;
            		upLeft=   -1;
            		up=       0;
    			} else {
    				p2=srcImage[index-width];p3=srcImage[index-width+1];
            		upRight=  -width+1;
            		upLeft=   -width-1;
            		up=       -width;
    			}
    			if (iY==(height-1)){
    				p8=p5;p9=p6;
            		downRight= 1;
            		downLeft= -1;
            		down=      0;
    			} else {
    				p8=srcImage[index+width];p9=srcImage[index+width+1];
            		downRight= width+1;
            		downLeft=  width-1;
            		down=      width;
    			}
    			p1=p2;p4=p5;p7=p8;

    			for (int iX=0;iX<width;iX++){
    				s1=p1+2*p2+p3-p7-2*p8-p9;
    				s2=p1+2*p4+p7-p3-2*p6-p9;
    				if (!spread || (iX==0) || (iX==(width-1))) sobelResult[index++]+=Math.sqrt(s1*s1+s2*s2);
    				else {
    					double d=Math.sqrt(s1*s1+s2*s2);
    					double s12=Math.abs(s1)*edgeSpreadOrtho;
    					double s22=Math.abs(s2)*edgeSpreadOrtho;
    					double sp2=Math.abs(s1+s2)*edgeSpreadDiagonal;
    					double sm2=Math.abs(s1-s2)*edgeSpreadDiagonal;

						sobelResult[index+1]+=s12;
						sobelResult[index-1]+=s12;
						sobelResult[index+up]+=  s22;
						sobelResult[index+down]+=s22;
						sobelResult[index+upLeft]+=   sm2;
						sobelResult[index+downRight]+=sm2;
						sobelResult[index+upRight]+=  sp2;
						sobelResult[index+downLeft]+= sp2;
						sobelResult[index++]+=d;
    				}
    				p1=p2;p2=p3;
    				p4=p5;p5=p6;
    				p7=p8;p8=p9;
    				if (iX<(width-2)){ // index already incremented, but iX - not yet
    					p3=srcImage[index+upRight];
    					p6=srcImage[index+1];
    					p9=srcImage[index+downRight];
    				}
    			}
    		}
    		return sobelResult;

    	}
    	public double [] sobelFilter0(
    			double [] srcImage,
    			int width,
    			double edgeSpreadOrtho,
    			double edgeSpreadDiagonal
    	){
    		int height=srcImage.length/width;
    		double p1,p2,p3,p4,p5,p6,p7,p8,p9,s1,s2;
    		double [] sobelResult=new double [srcImage.length];
    		int upRight,downRight;
    		int index;
    		boolean swap=(edgeSpreadOrtho<0);
    		edgeSpreadOrtho=Math.abs(edgeSpreadOrtho);
    		boolean spread=(edgeSpreadOrtho>0.0) || (edgeSpreadDiagonal>0.0);
//    		if (spread)
    			for (int i=0;i<sobelResult.length;i++) sobelResult[i]=0.0;
    		double ka=4.0/Math.PI;
    		double kOrto0=1.0-edgeSpreadOrtho;
    		double kOrto1=0.5*edgeSpreadOrtho;
    		double kDiagonal0=1.0-edgeSpreadDiagonal;
    		double kDiagonal1=0.5*edgeSpreadDiagonal;
    		System.out.println ("kOrto0="+kOrto0+" kOrto1="+kOrto1+" kDiagonal0="+kDiagonal0+" kDiagonal="+kDiagonal1+" swap="+swap+", spread="+spread);
    		int upLeft,downLeft,down,up;
    		for (int iY=0;iY<height;iY++){
    			index=iY*width;
    			p5=srcImage[index];p6=srcImage[index+1];
    			if (iY==0){
    				p2=p5;p3=p6;
            		upRight=  1;
            		upLeft=   -1;
            		up=       0;
    			} else {
    				p2=srcImage[index-width];p3=srcImage[index-width+1];
            		upRight=  -width+1;
            		upLeft=   -width-1;
            		up=       -width;
    			}
    			if (iY==(height-1)){
    				p8=p5;p9=p6;
            		downRight= 1;
            		downLeft= -1;
            		down=      0;
    			} else {
    				p8=srcImage[index+width];p9=srcImage[index+width+1];
            		downRight= width+1;
            		downLeft=  width-1;
            		down=      width;
    			}
    			p1=p2;p4=p5;p7=p8;

    			for (int iX=0;iX<width;iX++){
    				s1=p1+2*p2+p3-p7-2*p8-p9;
    				s2=p1+2*p4+p7-p3-2*p6-p9;
//    				if (!spread || (iX==0) || (iX==(width-1))) sobelResult[index++]+=Math.sqrt(s1*s1+s2*s2);
//    				if (!spread || (iX==0) || (iX==(width-1))) sobelResult[index++]+=ka*Math.atan2(s1,s2);
    				if (!spread || (iX==0) || (iX==(width-1))) sobelResult[index++]+=(((int) Math.round(ka*Math.atan2(s1,s2)))+(swap?8:6))%4;
    				else {
    					double d=Math.sqrt(s1*s1+s2*s2);
    					int dir=(((int) Math.round(ka*Math.atan2(s1,s2)))+(swap?8:6))%4;
    					switch (dir){
    					case 0: // vertical gradient
    						sobelResult[index+1]+=d*kOrto1;
    						sobelResult[index-1]+=d*kOrto1;
    						sobelResult[index++]+=d*kOrto0;
    						break;
    					case 1: //up-right or down-left gradient
    						sobelResult[index+upLeft]+=   d*kDiagonal1;
    						sobelResult[index+downRight]+=d*kDiagonal1;
    						sobelResult[index++]+=        d*kDiagonal0;
    						break;
    					case 2: // horizontal gradient
    						sobelResult[index+up]+=  d*kOrto1;
    						sobelResult[index+down]+=d*kOrto1;
    						sobelResult[index++]+=   d*kOrto0;
    						break;
    					case 3:// down-right or up-left gradient
    						sobelResult[index+upRight]+=  d*kDiagonal1;
    						sobelResult[index+downLeft]+= d*kDiagonal1;
    						sobelResult[index++]+=        d*kDiagonal0;
    						break;
    						default:
    							System.out.print(" dir="+dir);
    					}
    				}
    				p1=p2;p2=p3;
    				p4=p5;p5=p6;
    				p7=p8;p8=p9;
    				if (iX<(width-2)){ // index already incremented, but iX - not yet
    					p3=srcImage[index+upRight];
    					p6=srcImage[index+1];
    					p9=srcImage[index+downRight];
    				}
    			}
    		}
    		return sobelResult;

    	}

    	/**
    	 * Calculate difference between pixels on the two images along a section line as (weighted_RMS+noiseLev)/(weighted distance+noiseLev)
    	 * @param sectionY Position of the section line
    	 * @param side 0 - left image, 1 - right
    	 * @param side 0 - correlate with other image, 1 - auto-correlate (to detect periodic)
    	 * @param weightSobel (0.0 ... 1.0) - weight of the image Y component processed by a Sobel edge-detection
    	 *  (Y will have 1-weightSobel) weight
    	 * @param weightCb weight of the image Cb component as fraction of Y
    	 * @param weightCr weight of the image Cr component as fraction of Y
    	 * @param noiseLev pixel noise as a fraction of the full scale (1.0)
    	 * @return 2d packed array - measure of similarity between matching pixels on the two images, same format as in tilesSection() and alphaSection()
    	 */

    	public double [] sectionPixelDifference(
    			int sectionY,
    			int side,
    			int self,
    			double weightSobel,
    			double weightCb,
    			double weightCr,
    			double noiseLev, // negative - do not normalize to the RMS
    			int debugLevel){
    		if (sectionY<1) sectionY=1; // to skip tests for the valid neighbors
    		if (sectionY> (this.mapHeight-2)) sectionY=this.mapHeight-2;
    		int sectionLength=this.mapWidth;
    		int numSamples= this.maxDisparity-this.minDisparity+1;
    		double [] result = new double[sectionLength*numSamples];
    		for (int i=0;i<result.length;i++) result[i]=0.0;
    		int numLayers=this.overlapImages.length/2;
    		int iYThis= 1+((side==0)?0:numLayers);
    		int iYOther=1+(((side!=0)^(self!=0))?0:numLayers);
    		int sideOther=(self==0)?(1-side):side;
    		boolean useSobel= (weightSobel>0);
    		if (useSobel) initSobelY(false); // only if it does not exist
    		double [][][] images={{
    				this.overlapImages[iYThis],   //Y
    				(useSobel)?this.sobelY[side]:this.overlapImages[iYThis], // Y-edges
    				this.overlapImages[iYThis+1], // Cb
    				this.overlapImages[iYThis+2]}, // Cr
    				{this.overlapImages[iYOther],   //Y
    				(useSobel)?this.sobelY[sideOther]:this.overlapImages[iYOther], // Y-edges
    				this.overlapImages[iYOther+1], // Cb
    				this.overlapImages[iYOther+2]}}; // Cr
    		int [] dirs={0,1,1+this.mapWidth,this.mapWidth,-1+this.mapWidth, -1, -1 -this.mapWidth, -this.mapWidth, 1-this.mapWidth};
    		double [] weights={
    				(1.0-weightSobel)/(1.0+weightSobel+weightCb+weightCr),
    				weightSobel/(1.0+weightSobel+weightCb+weightCr),
    				weightCb/(1.0+weightSobel+weightCb+weightCr),
    				weightCr/(1.0+weightSobel+weightCb+weightCr)};
    		for (int iX=1;iX<(sectionLength-1);iX++)for (int iDispar=0;iDispar<numSamples;iDispar++){
    			int iXOther=(side==0)?(iX-iDispar):(iX+iDispar); // same for self - right image in opposite direction
    			if ((iXOther>1) && (iXOther<(this.mapWidth-1))){
    				int indexThis= sectionY*this.mapWidth+iX;
    				int indexOther=sectionY*this.mapWidth+iXOther;
    				double weightedRMS=0.0;
    				double weightedDiff=0.0;
    				double [][] SX= new double[2][weights.length];
    				double [][] SX2=new double[2][weights.length];
    				for (int i=0;i<2;i++) for (int j=0;j<weights.length;j++){
    					SX[i][j]=0.0;
    					SX2[i][j]=0.0;
    				}
    				// Don't be too smart and use simultaneous calculation of the SX and SX2 - got square root of -2.220446049250313E-16,
    				// so first calculate and subtract average
    				for (int nComp=0;nComp<weights.length;nComp++){
                        double d= images[0][nComp][indexThis]-images[1][nComp][indexOther];
                        weightedDiff+=	0.5*weights[nComp]*d*d;
                        d= images[0][nComp][indexThis-1]-images[1][nComp][indexOther-1];
                        weightedDiff+=	0.25*weights[nComp]*d*d;
                        d= images[0][nComp][indexThis+1]-images[1][nComp][indexOther+1];
                        weightedDiff+=	0.25*weights[nComp]*d*d;

                        for (int iDir=0;iDir<dirs.length;iDir++){ // 9 points
    					  SX[0][nComp]+=images[0][nComp][indexThis+dirs[iDir]];
    					  SX[1][nComp]+=images[1][nComp][indexOther+dirs[iDir]];
    					}
    					SX[0][nComp]/=dirs.length; // make it average
    					SX[1][nComp]/=dirs.length; // make it average

    				}

    				for (int nComp=0;nComp<weights.length;nComp++){
    					double d;
    					for (int iDir=0;iDir<dirs.length;iDir++){ // 9 points
    					  d=images[0][nComp][indexThis+dirs[iDir]]-SX[0][nComp];
    					  SX2[0][nComp]+=d*d;
    					  d=images[1][nComp][indexOther+dirs[iDir]]-SX[1][nComp];
    					  SX2[1][nComp]+=d*d;
    					}
    					weightedRMS+=weights[nComp]*(SX2[0][nComp]+SX2[1][nComp])/(2*dirs.length);
/*
    					if (debugLevel>1) {
    						System.out.print("sectionPixelDifference(): iX="+iX+" iDispar="+iDispar+" iXOther="+iXOther+
    								" nComp="+nComp+ " SX[0][nComp]="+SX[0][nComp]+" SX2[0]="+SX2[0][nComp]+" SX[1]="+SX[1][nComp]+" SX2[1]="+SX2[1][nComp]);
    					}
    					if (debugLevel>1) {
    						System.out.println(" ##   weightedRMS="+weightedRMS+ " weightedDiff="+weightedDiff);
    					}
*/
    				}
    				weightedDiff=Math.sqrt(weightedDiff);
    				weightedRMS= Math.sqrt(weightedRMS);
    				if (noiseLev>=0) result[iX*numSamples+iDispar]=(weightedRMS+noiseLev)/(weightedDiff+noiseLev);
    				else  result[iX*numSamples+iDispar]=(-noiseLev)/(weightedDiff-noiseLev);
/*
					if (debugLevel>1) {
						System.out.println(">>> weightedRMS="+weightedRMS+ " weightedDiff="+weightedDiff+ " result[iX*numSamples+iDispar]="+result[iX*numSamples+iDispar]+
					" noiseLev="+noiseLev);
					}
*/
    			}
    		}
    		return result;
    	}

    	public double [] sectionPixelSobel(
    			int sectionY,
    			int side,
    			int self,
    			int debugLevel){
    		if (sectionY<1) sectionY=1; // to skip tests for the valid neighbors
    		if (sectionY> (this.mapHeight-2)) sectionY=this.mapHeight-2;
    		int sectionLength=this.mapWidth;
    		int numSamples= this.maxDisparity-this.minDisparity+1;
    		double [] result = new double[sectionLength*numSamples];
    		for (int i=0;i<result.length;i++) result[i]=0.0;
    		initSobelY(false); // only if it does not exist
    		int sideOther=(self==0)?(1-side):side;
    		for (int iX=1;iX<(sectionLength-1);iX++)for (int iDispar=0;iDispar<numSamples;iDispar++){
    			int iXOther=(side==0)?(iX-iDispar):(iX+iDispar);
    			if ((iXOther>1) && (iXOther<(this.mapWidth-1))){
    				int indexThis= sectionY*this.mapWidth+iX;
    				int indexOther=sectionY*this.mapWidth+iXOther;
    				result[iX*numSamples+iDispar]=this.sobelY[side][indexThis]*this.sobelY[sideOther][indexOther];
    			}
    		}
    		return result;
    	}






    	public double[]  removePeriodicSection(
    			double [] imageCross,
    			double [] imageAuto,
    			double minPeriod,
    			double minZeroAuto,
    			double minPeriodContrast,
    			double minAbsoluteFraction,
    			double blurMaskSigma,
    			double scale,
    			int debugLevel){
    		int dbgMin=360,dbgMax=380;

    		int sectionLength=this.mapWidth;
    		int numSamples= this.maxDisparity-this.minDisparity+1;
    		double [] result = new double[sectionLength*numSamples];
    		for (int i=0;i<result.length;i++) result[i]=0.0;
    		for (int iX=0;iX<sectionLength;iX++){
    			double [] cross=new double [numSamples];
    			double [] auto= new double [numSamples];
				int index=iX*numSamples;
    			for (int i=0;i<numSamples;i++){
    				cross[i]=imageCross[index];
    				auto[i]= imageAuto [index++];
    			}
    			int dbg=debugLevel;
    			if ((iX>=dbgMin) && (iX<=dbgMax)) dbg+=1;
    			if (dbg>1) System.out.print("iX="+iX+" >> ");
    			double [] combined=removePeriodic(
    	    			cross,
    	    			auto,
    	    			minPeriod,
    	    			minZeroAuto,
    	    			minPeriodContrast,
    	    			minAbsoluteFraction,
    	    			blurMaskSigma,
    	    			scale,
    	    			dbg); //debugLevel);
    			if (combined!=null) cross=combined;
				index=iX*numSamples;
    			for (int i=0;i<numSamples;i++){
    				result[index++]=cross[i];
    			}
    		}
    		return result;
    	}




    	/**
    	 * Combine correlation data with auto-correlation to remove ambiguity caused by periodic patterns
    	 * @param cross Array of correlation results (starting with this.minDisparity, usually -1)
    	 * @param auto Array of auto-correlation (same dimension as cross
    	 * @param minPeriod Minimal period to consider
    	 * @param minZeroAuto Minimal autocorrelation value at zero to process
    	 * @param minPeriodContrast Minimal ratio of the first auto-correlation maximum to autocorrelation at zero
    	 * @param minAbsoluteFraction Minimal ratio of the first auto-correlation maximum to the absolute autocorrelation maximum,
    	 * @param blurMaskSigma blur periodic rejection mask with this sigma
    	 * @param scale Scale rejection mask (scale==1 will reduce the false disparity contrast twice (if the auto had the same first amplitude as zero)
    	 * @param debugLevel
    	 * @return null if did not match criteria, otherwise the array to be used instead of cross.
    	 */

    	// TODO: What if we found not the fundamental period but a multiple?
    	// also - does not work (yet) if period is less than disparity

    	public double [] removePeriodic(
    			double [] cross,
    			double [] auto,
    			double minPeriod,
    			double minZeroAuto,
    			double minPeriodContrast,
    			double minAbsoluteFraction,
    			double blurMaskSigma,
    			double scale,
    			int debugLevel){
    		boolean dbg=(minPeriodContrast<0);
    		minPeriodContrast=Math.abs(minPeriodContrast);
    		int numSamples= this.maxDisparity-this.minDisparity+1;
    		if (debugLevel>2){
    			for (int i=0;i<numSamples;i++) System.out.print(" "+IJ.d2s(auto[i],4));
        		System.out.println("");
    		}
    		if (auto[-this.minDisparity]<=0.0){
    			System.out.println("removePeriodic(): Seems to be a BUG - autocorrelation at zero ("+auto[-this.minDisparity]+") <=0.0");
    			return null;
    		}
    		if (auto[-this.minDisparity]<=minZeroAuto){
    			if (debugLevel>2) {
    			   System.out.println("removePeriodic(): autocorrelation at zero ("+auto[-this.minDisparity]+") <= "+minZeroAuto);
    			}
    			return null;
    		}

    		int minIndex=(int) minPeriod-this.minDisparity;
    		int halfMinPeriod=(int) Math.ceil(minPeriod/2);
    		double max=0.0;
    		for (int i=minIndex; i<(numSamples-halfMinPeriod);i++) if (auto[i]>max) max=auto[i];
			double threshold=minPeriodContrast*auto[-this.minDisparity]; // minimal threshold
			if (threshold<(max*minAbsoluteFraction)) threshold=max*minAbsoluteFraction;
    		int iMax=-1;
    		if (minIndex<1) minIndex=1;
    		// looking for the first local max that fits
			boolean localMax;
    		for (int i=minIndex;i<(numSamples-halfMinPeriod);i++) if (auto[i]>threshold){
    			localMax=true;
    			for (int j=-halfMinPeriod; j<halfMinPeriod;j++) if (auto[i+j]>auto[i]){
    				localMax=false;
    				break;
    			}
    			if (localMax) { // also - above threshold
    				iMax=i;
    				break;
    			}
    		}

    		if (iMax<0) {
    			if (debugLevel>1) System.out.println("removePeriodic(): No qualifying period above "+minPeriod+
    					" detected, minIndex="+minIndex+" auto["+(-this.minDisparity)+"]="+auto[-this.minDisparity]);
    			return null; // no maximum above minPeriod
    		}

    		max=auto[iMax];
    		/*
// Make sure there are no local maximums between zero and this ?
    		if (!(max>(minPeriodContrast*auto[-this.minDisparity]))){
    			if (debugLevel>1) System.out.println("removePeriodic(): contrast too low ("+max+"<="+minPeriodContrast*auto[-this.minDisparity]);
    			return null; // contrast too low
    		}
    		*/
    		// find minimal value between zero max and this new found one
    		int iMin=0;
    		double min= auto[-this.minDisparity];
    		for (int i=-this.minDisparity+1;i<iMax;i++) if (auto[i]<min){
    			iMin=i;
    			min=auto[i];
    		}
    		int period=iMax-this.minDisparity;

			if (debugLevel>1) System.out.println("removePeriodic(): period ="+period+" auto["+(-this.minDisparity)+"]="+auto[-this.minDisparity]+
					" max="+max+" min="+min+" iMin="+iMin);
			if (min<0) min=0.0; // change floor

    		//double [] mask=cross.clone();
			double [] mask=new double [cross.length];
    		double k=scale/auto[-this.minDisparity];
    		for (int i=0;i<numSamples;i++){
    			int iSrc=i-period;
    			mask[i]=((iSrc>=0) && (cross[iSrc]>min))? (k*(cross[iSrc]-min)):0.0;
    		}
    		if (blurMaskSigma>0.0) {
    			(new DoubleGaussianBlur()).blur1Direction(
    					mask,
    					numSamples,
    					1,
    					blurMaskSigma,
    					0.01,
    					true);
    		}
    		double [] result=cross.clone();
    		if (dbg) for (int i=0;i<numSamples;i++) result[i]=1.0;
    		for (int i=period;i<numSamples;i++) result[i]/=1+mask[i];
    		return result;
    		/*
    		 double [] mask=auto.clone();
    		double k=scale/auto[-this.minDisparity];
    		for (int i=0;i<numSamples;i++) mask[i]=((i>iMin) && (auto[i]>min))? (k*(auto[i]-min)):0.0;
    		if (blurMaskSigma>0.0) {
    			(new DoubleGaussianBlur()).blur1Direction(
    					mask,
    					numSamples,
    					1,
    					blurMaskSigma,
    					0.01,
    					true);
    		}
    		double [] result=cross.clone();
    		for (int i=period;i<numSamples;i++) result[i]/=1+mask[i];
    		return dbg?mask:result;

    		 */
    	}


    	/**
    	 * Setup coverage of the overlap area with tiles of the specified size
    	 * @param size Tile size (half of the windowed FFT size)
    	 * @param levels Number of levels of progressively smaller tiles
    	 * @param minDisparity Minimal disparity value for accumulation (-1)
    	 * @param maxDisparity Maximal disparity value for accumulation (size/2)
    	 * @return number of tiles horizontally and vertically
    	 */
    	public int [] setupTiles(
    			int size,
    			int levels, // number of tile levels, each having half linear size (4..5)
    			int minDisparity, // -1
    			int maxDisparity  // =size/2
    			){
    		this.minDisparity=minDisparity;
    		this.maxDisparity=maxDisparity;

    		int numLayers=this.overlapImages.length/2;

    		int minY=this.mapHeight,minX=this.mapWidth,maxY=0,maxX=0;

    		for (int iy=0;iy<this.mapHeight;iy++) for (int ix=0;ix<this.mapWidth;ix++){
    			int index=iy*this.mapWidth+ix;
    			if ((this.overlapImages[0][index]>0) && (this.overlapImages[numLayers][index]>0)){
    				if (iy>maxY) maxY=iy;
    				if (iy<minY) minY=iy;
    				if (ix>maxX) maxX=ix;
    				if (ix<minX) minX=ix;
    			}
    		}
    		if (maxY<minY) {
    			if (this.debugLevel>0) System.out.println("No overlap");
    			return null;
    		}
    		int width= (maxX-minX+1);
    		int height=(maxY-minY+1);
    		this.tilesNHor=width/size;
    		if ((this.tilesNHor*size)<width) this.tilesNHor++;
    		this.tilesNVert=height/size;
    		if ((this.tilesNVert*size)<height) this.tilesNVert++;
    		// center tiles over overlap area
    		this.tilesX0=(maxX+minX+1)/2-(this.tilesNHor*size)/2;
    		this.tilesY0=(maxY+minY+1)/2-(this.tilesNVert*size)/2;
    		int [] nTiles={this.tilesNHor,this.tilesNVert};
    		this.tilesSize=size;
    		this.tiles = new double[2][2][levels][][][];
    		for (int side=0;side<2;side++)for (int self=0;self<2;self++){
    			for (int nLev=0;nLev<levels;nLev++){
    				this.tiles[side][self][nLev]=new double [this.tilesNVert<<nLev][this.tilesNHor<<nLev][];
    				for (int tileY=0;tileY<this.tiles[side][self][nLev].length;tileY++)
    					for (int tileX=0;tileX<this.tiles[side][self][nLev][0].length;tileX++) this.tiles[side][self][nLev][tileY][tileX]=null;
    			}
    		}
    		return nTiles;
    	}
    	/**
    	 * Combined tile disparity (from correlation) array with those of parent and neighbors
    	 * @param parentFraction Multiply by (1+max(parentFraction*parent(x),0))
    	 * @param neighborsFraction Multiply by (1+max(Multiply by (1+max(neighborsFraction*neighbors(x),0))[x),0))
    	 * @param neighborsPower amplify if any of the neighbors share the same disparity, so value is calculated as RMS bnut with variable power instead of 2
    	 * @param debugLevel Debug level
    	 */
		public void applyParentAndNeighbors(
				double parentFraction,
				double neighborsFraction,
				double neighborsPower,
				int debugLevel){
			int [][]dirs={{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1},{-1,0},{-1,1}};
			for (int side =0;side<this.tiles.length;side++){
				for (int self=0;self<2;self++){
					for (int level=0;level<this.tiles[side][self].length;level++){
						double [][][] tmpTiles= new double [this.tiles[side][self][level].length][][];
						for (int tileY=0;tileY<this.tiles[side][self][level].length; tileY++){
							tmpTiles[tileY]=new double[this.tiles[side][self][level][tileY].length][];
							for (int tileX=0;tileX<this.tiles[side][self][level][tileY].length; tileX++){
								tmpTiles[tileY][tileX]=null;
								if (this.tiles[side][self][level][tileY][tileX]!=null) {
									tmpTiles[tileY][tileX]=this.tiles[side][self][level][tileY][tileX].clone();
									double [] tile=tmpTiles[tileY][tileX];
									if ((parentFraction>0.0) && (level>0)){
										double [] parentTile=this.tiles[side][self][level-1][tileY/2][tileX/2]; // should never be null
										for (int i=0;i<tile.length;i++) if (parentTile[i]>0.0){
											tile[i]*=(1.0+parentFraction*parentTile[i]);
										}
									}
									if (neighborsFraction>0.0){
										double [] neighbors=new double[tile.length];
										int [] numUsed=  new int[tile.length];
										for (int i=0;i<neighbors.length;i++){
											neighbors[i]=0.0;
											numUsed[i]=  0;
										}

										for (int iDir=0;iDir<dirs.length;iDir++){
											int nTileX=tileX+dirs[iDir][0];
											int nTileY=tileY+dirs[iDir][1];
											if (debugLevel>10){
												System.out.println(
														" tmpTiles.length="+tmpTiles.length+
														" tmpTiles[0.length="+tmpTiles[0].length+
														" nTileY="+nTileY+
														" nTileX="+nTileX+
														" side="+side+
														" self="+self+
														" level="+level);
												if ((nTileX>=0) && (nTileX<tmpTiles[tileY].length)&&
														(nTileY>=0) && (nTileY<tmpTiles.length)){
													System.out.print("this.tiles.length=");System.out.print(this.tiles.length);
													System.out.print("this.tiles["+side+"].length=");System.out.print(this.tiles[side].length);
													System.out.print("this.tiles["+side+"]["+self+"].length=");System.out.print(this.tiles[side][self].length);
													System.out.print("this.tiles["+side+"]["+self+"]["+level+"].length=");System.out.print(this.tiles[side][self][level].length);
													System.out.print("this.tiles["+side+"]["+self+"]["+level+"]["+nTileY+"].length=");System.out.print(this.tiles[side][self][level][nTileY].length);
													System.out.print("this.tiles["+side+"]["+self+"]["+level+"]["+nTileY+"]["+nTileX+"].length=");System.out.print(this.tiles[side][self][level][nTileY][nTileX].length);
												}
											}
											if ( // oob 3
													(nTileX>=0) && (nTileX<tmpTiles[tileY].length)&&
													(nTileY>=0) && (nTileY<tmpTiles.length)&&
													(this.tiles[side][self][level][nTileY][nTileX]!=null)){
												double [] nTile=this.tiles[side][self][level][nTileY][nTileX];
												for (int i=0;i<nTile.length;i++) if (nTile[i]>0){
													neighbors[i]+=Math.pow(nTile[i],neighborsPower);
													numUsed[i]++;
												}
											}
										}
										double rPow=1.0/neighborsPower;
										for (int i=0;i<tile.length;i++) if (numUsed[i]>0){
											tile[i]*=(1.0+neighborsFraction* Math.pow(neighbors[i]/numUsed[i],rPow));

										}

									}
								}
							}

						}
						// copy layer back to orinal array
						this.tiles[side][self][level]=tmpTiles;
					}
				}
			}

		}



    	/**
    	 * Generate tile structure of the disparity values based on correlation between two images. Each (non-null) tile stores an
    	 * array of correlation contrasts for each integer disparity value from this.disparityMin to this.disparityMax (inclusive)
    	 * @param sides  1 - left image, 2 - right image - 3 both
    	 * @param selves 1 - only correlation with other, 2 - only autocorrelation,  3 - both
    	 * @param phaseCorrelationFraction 1.0 - use pure phase correlation, 0.0 - just correlation
    	 * @param corrCbWeight relative weight of the Cb component in correlation (Y has weight of 1.0)
    	 * @param corrCrWeight relative weight of the Cr component in correlation (Y has weight of 1.0)
    	 * @param correlationHighPassSigma Correlation high-pass filter sigma,  pixels in frequency domain
    	 * @param correlationLowPassSigma Correlation low-pass filter sigma,  fraction of the frequency range
    	 * @param noiseNormalizationSignaY Noise normalization sigma for Y-component
    	 * @param noiseNormalizationSignaCbCr Noise normalization sigma for Cb, Cr components
    	 * @param contrastThreshold Use correlation data with normalized values above this threshold only
    	 * @param contrastThresholdDecrease Decrease threshold for smaller tiles as 1/2^(contrastThresholdDecrease*level)
    	 * @param relativeStep Maximal distance between correlation shifts relative to size (if distance from the new
    	 *  shift and already calculated is less than that, new correlation will not be performed
    	 * @param probeAtLimits Add shift near the minimal and maximal shifts over threshold even if they are not local maximums
    	 * @param binaryAlpha consider all points with alpha>0 to have same alpha==1.0
    	 * @param minTileFraction minimal tile overlap fraction to process tile (will be squared in the corners of the tile map)
    	 * @param doubleFHT instance of DoubleFHT class (or null). Used to optimize by eliminating repetitive sin/cos/window tables generation
    	 * @param threadsMax Maximal number of concurrent threads to run
    	 * @param debugLevel Debug level
    	 * @return Number of non-null tiles
    	 */


    	public int generateTileLevels(
    			int sides,
    			int selves,
    			double phaseCorrelationFraction,
				double corrCbWeight,
				double corrCrWeight,
    			double correlationHighPassSigma,
    			double correlationLowPassSigma,
    			double noiseNormalizationSignaY,
    			double noiseNormalizationSignaCbCr,
    			double contrastThreshold,
    			double contrastThresholdDecrease, //Decrease threshold for smaller tiles as 1/2^(k*level), k=
    			double relativeStep,
    			boolean probeAtLimits,
    			boolean binaryAlpha,
    			double minTileFraction,
    			DoubleFHT doubleFHT,
    			int threadsMax,
    			int debugLevel){
    		if (doubleFHT==null) doubleFHT=new DoubleFHT();
    		int numNonEmptyTiles=0;
			if (debugLevel>1) System.out.println("generateTileLevels("+sides+", "+selves+"...), this.tiles.length="+this.tiles.length);
    		for (int side=0; side<this.tiles.length;side++) if ((sides & (1<<side))!=0) {
				if (debugLevel>1) System.out.println("side="+side+", this.tiles["+side+"].length="+this.tiles.length);
    			for (int self=0; self<this.tiles[side].length;self++) if ((selves & (1<<self))!=0) {
    				if (debugLevel>1) System.out.println("side="+side+" self="+self+", this.tiles["+side+"]["+self+"].length="+this.tiles.length);
    				for (int level=0;level<this.tiles[side][self].length;level++){
    					double contastThresholdMod=contrastThreshold*Math.pow(0.5,contrastThresholdDecrease*level);
    					int tileResult=generateTileLevel(
    							side,
    							self,
    							level,
    							phaseCorrelationFraction,
    							corrCbWeight,
    							corrCrWeight,
    							correlationHighPassSigma,
    							correlationLowPassSigma,
    							noiseNormalizationSignaY,
    							noiseNormalizationSignaCbCr,
    							contastThresholdMod,
    							relativeStep,
    							probeAtLimits,
    							binaryAlpha,
    							minTileFraction,
    							doubleFHT,
    							threadsMax,
    							debugLevel);
    					if (debugLevel>1){
    						System.out.println("On level "+level+": non-empty tiles="+tileResult);
    					}
    					numNonEmptyTiles+= tileResult;
    				}
    			}
    		}
    		return numNonEmptyTiles;
    	}

// Convert to multithread
    	public int  generateTileLevel(
    			int side,
    			int self,
    			int level,
				double phaseCorrelationFraction,
				double corrCbWeight,
				double corrCrWeight,
				double correlationHighPassSigma,
				double correlationLowPassSigma,
				double noiseNormalizationSignaY,
				double noiseNormalizationSignaCbCr,
				double contrastThreshold,
    			double relativeStep,
    			boolean probeAtLimits,
				boolean binaryAlpha,
				double minTileFraction,
				DoubleFHT doubleFHT,
				int threadsMax,
				int debugLevel){
    		int numNonEmptyTiles=0;
			if (debugLevel>1) System.out.println("generateTileLevel("+side+", "+self+"...)");

    		for (int tileY=0;tileY<this.tiles[side][self][level].length;tileY++) for (int tileX=0;tileX<this.tiles[side][self][level][0].length;tileX++){
    			double [] tileResult=generateTile(
    					side,
    					self,
    	    			level,
    	    			tileX,
    	    			tileY,
    					phaseCorrelationFraction,
    					corrCbWeight,
    					corrCrWeight,
    					correlationHighPassSigma,
    					correlationLowPassSigma,
    					noiseNormalizationSignaY,
    					noiseNormalizationSignaCbCr,
    					contrastThreshold,
    	    			relativeStep,
    	    			probeAtLimits,
    					binaryAlpha,
    					minTileFraction,
    					doubleFHT,
    					threadsMax,
    					debugLevel);
    			if (tileResult!=null) {
    				numNonEmptyTiles++;
    			}
    		}
    		return numNonEmptyTiles;
    	}

    	/**
    	 * Generate disparity tile by correlating the images, starting with info from the higher level tile (if available)
    	 * @param side 0 - left image, 1 - right image
    	 * @param self 0 - correlation with other, 1 - autocorrelation
    	 * @param level Level of the tile to generate
    	 * @param tileX Horizontal tile number
    	 * @param tileY Vertical tile number
    	 * @param phaseCorrelationFraction 1.0 - use pure phase correlation, 0.0 - just correlation
    	 * @param corrCbWeight relative weight of the Cb component in correlation (Y has weight of 1.0)
    	 * @param corrCrWeight relative weight of the Cr component in correlation (Y has weight of 1.0)
    	 * @param correlationHighPassSigma Correlation high-pass filter sigma,  pixels in frequency domain
    	 * @param correlationLowPassSigma Correlation low-pass filter sigma,  fraction of the frequency range
    	 * @param noiseNormalizationSignaY Noise normalization sigma for Y-component
    	 * @param noiseNormalizationSignaCbCr Noise normalization sigma for Cb, Cr components
    	 * @param contrastThreshold Use correlation data with normalized values above this threshold only
    	 * @param relativeStep Maximal distance between correlation shifts relative to size (if distance from the new
    	 *  shift and already calculated is less than that, new correlation will not be performed
    	 * @param probeAtLimits Add shift near the minimal and maximal shifts over threshold even if they are not local maximums
    	 * @param binaryAlpha consider all points with alpha>0 to have same alpha==1.0
    	 * @param minTileFraction minimal tile overlap fraction to process tile (will be squared in the corners of the tile map)
    	 * @param threadsMax Maximal number of concurrent threads to run
    	 * @param debugLevel Debug level
    	 * @return Generated tile (or null if no tile was generated)
    	 */
    	public double [] generateTile(
    			int side,
    			int self,
    			int level,
    			int tileX,
    			int tileY,
				double phaseCorrelationFraction,
				double corrCbWeight,
				double corrCrWeight,
				double correlationHighPassSigma,
				double correlationLowPassSigma,
				double noiseNormalizationSignaY,
				double noiseNormalizationSignaCbCr,
				double contrastThreshold,
    			double relativeStep,
    			boolean probeAtLimits,
				boolean binaryAlpha,
				double minTileFraction,
				DoubleFHT doubleFHT,
				int threadsMax,
				int debugLevel){
    		int size=(this.tilesSize>>level);
    		boolean isCorner=(
    				((tileY==0)||(tileY==(this.tiles[side][self][level].length))) &&
    				((tileX==0)||(tileX==(this.tiles[side][self][level][0].length))));
    		if (isCorner) minTileFraction*=minTileFraction;
    		int corrXC=this.tilesX0+tileX*size+size/2;
    		int corrYC=this.tilesY0+tileY*size+size/2;

    		double [] parentTile;
    		if (level==0)	parentTile=null;
    		else {
    			parentTile= this.tiles[side][self][level-1][tileY/2][tileX/2];
    			if (parentTile==null) return null; // no parent to this tile
    		}
    		if (debugLevel>2){
    			System.out.println("======> generateTile("+side+", "+self+", "+level+", "+tileX+", "+tileY+",...)");
    		}
    		// first series of correlations using shift from the higher level
    		int [] shifts =  generateShiftsList(
    				contrastThreshold, // double threshold,
        			size, // int size,
        			relativeStep, //double relativeStep,
        			probeAtLimits, //boolean probeAtLimits,
        			null, // int [] currentShifts,
        			parentTile, //double [] prevCorrelation,
        			debugLevel //int debugLevel
        			);
    		if (shifts==null){

    			return null;
    		}

    		if (debugLevel>2){
    			System.out.print("Accumulating - first pass: ");
    		}
    		double [][] accumulatedCorrs=correlateMultipleShifts(
    				side,
    				self,
        			size,
            		corrXC,
            		corrYC,
        			shifts,    // >=0;
    				phaseCorrelationFraction,
    				corrCbWeight,
    				corrCrWeight,
    				correlationHighPassSigma,
    				correlationLowPassSigma,
    				noiseNormalizationSignaY,
    				noiseNormalizationSignaCbCr,
    				binaryAlpha,
    				minTileFraction, // decrease for corners
    				doubleFHT,
    				threadsMax,
    				debugLevel);
    		if (accumulatedCorrs==null){
        		if (debugLevel>2){
        			System.out.println("First pass got null");
        		}
    			return null;
    		}
    		// additional correlations shifted according to the new results
    		int [] newShifts=generateShiftsList(
    				contrastThreshold, // double threshold,
        			size, // int size,
        			relativeStep, //double relativeStep,
        			probeAtLimits, //boolean probeAtLimits,
        			shifts, // int [] currentShifts,
        			accumulatedCorrs[0], //parentTile, //double [] prevCorrelation,
        			debugLevel //int debugLevel
        			);
    		if (newShifts==null){
        		if (debugLevel>2){
        			System.out.println("First pass - all below threshold");
        		}
    			return null;
    		}
    		if (newShifts.length>0){ // length==0 -> no new shifts to try
        		if (debugLevel>2){
        			System.out.print("Accumulating - second pass: ");
        		}
    			double [][] additionalCorrs=correlateMultipleShifts(
    					side,
    					self,
            			size,
                		corrXC,
                		corrYC,
            			newShifts,    // >=0;
        				phaseCorrelationFraction,
        				corrCbWeight,
        				corrCrWeight,
        				correlationHighPassSigma,
        				correlationLowPassSigma,
        				noiseNormalizationSignaY,
        				noiseNormalizationSignaCbCr,
        				binaryAlpha,
        				minTileFraction, // decrease for corners
        				doubleFHT,
        				threadsMax,
        				debugLevel);
    			if (additionalCorrs!=null){ // combine accumulatedCorrs with additionalCorrs
    				for (int i=0;i<accumulatedCorrs[0].length;i++){
    					double w=accumulatedCorrs[1][i]+additionalCorrs[1][i];
    					if (w>0.0) accumulatedCorrs[0][i]=(accumulatedCorrs[0][i]*accumulatedCorrs[1][i]+additionalCorrs[0][i]*additionalCorrs[1][i])/w;
    				}
    			}
    		}
    		this.tiles[side][self][level][tileY][tileX]=accumulatedCorrs[0];
    		return this.tiles[side][self][level][tileY][tileX];
    	}



    	/**
    	 * Create new or additional list of the correlation shifts to try
    	 * @param threshold Minimal correlation contrast to grant
    	 * @param size Tile size (half of the correlation FFT size
    	 * @param relativeStep Maximal distance between correlation shifts relative to size (if distance from the new
    	 *  shift and already calculated is less than that, new correlation will not be performed
    	 * @param probeAtLimits Add shift near the minimal and maximal shifts over threshold even if they are not local maximums
    	 * @param currentShifts Array of correlation shifts already calculated (or null)
    	 * @param prevCorrelation Array of previous correlation results (or null )
    	 * @param debugLevel Debug level
    	 * @return array of non-negative correlation shifts to try (may be zero length if none are needed), null - all below threshold
    	 */

    	int[] generateShiftsList(
    			double threshold,
    			int size,
    			double relativeStep,
    			boolean probeAtLimits,
    			int [] currentShifts,
    			double [] prevCorrelation,
    			int debugLevel
    			){
    		int minStep=1; // make it 2?
    		int step=(int) Math.round(relativeStep*size);
    		if (step<minStep) step=minStep;
    		if (prevCorrelation==null){
    			int [] initialShifts={0,size/2}; // try two shifts? or just one {0} ?
        		if (debugLevel>2){
        			System.out.println("generateShiftsList(): prevCorrelation was null, returning initial list");
        		}
    			return initialShifts;
    		}
    		int minShift=prevCorrelation.length,maxShift=0;
    		boolean [] localMax=new boolean [prevCorrelation.length];
    		double maxV=0.0;
    		for (int i=0;i< localMax.length;i++) {
    			if (prevCorrelation[i]>=maxV) maxV=prevCorrelation[i];
    			if (prevCorrelation[i]>=threshold){
    				if (i>maxShift) maxShift=i;
    				if (i<minShift) minShift=i;
    				localMax[i]= ((i>=-this.minDisparity) && (i>0) && (prevCorrelation[i]>=prevCorrelation[i-1]) &&
    						(i<(prevCorrelation.length-1)) && (prevCorrelation[i]>prevCorrelation[i+1]));
        			if ((i==-this.minDisparity) && (i<(prevCorrelation.length-1)) && (prevCorrelation[i]>prevCorrelation[i+1])) localMax[i]=true; // only if above limit
    			}
    		}
    		if (minShift>maxShift) {
        		if (debugLevel>2){
        			System.out.println("generateShiftsList(): All below threshold: maxV="+maxV+" threshold="+threshold+" (length="+localMax.length+")");
        			for (int i=0;i<prevCorrelation.length;i++){
        				System.out.print(" "+IJ.d2s(prevCorrelation[i],3));
        			}
        			System.out.println("");
        		}
    			return null; // all below threshold
    		}
    		if (debugLevel>2){
    			System.out.print("generateShiftsList(): before probeAtLimits: ");
    			for (int i=0;i<localMax.length;i++) System.out.print(localMax[i]?"+":".");
    			System.out.println("");
    		}
    		if (probeAtLimits) {
    			int lshift=minShift+step/2;
    			if (lshift>maxShift) lshift=(maxShift+minShift)/2;
    			if ((this.minDisparity+lshift)<0) lshift=-this.minDisparity;
    			if (lshift<localMax.length) localMax[lshift]=true;
    			lshift=maxShift-step/2;
    			if (lshift>=minShift) localMax[lshift]=true;
        		if (debugLevel>2){
        			System.out.print("generateShiftsList(): after  probeAtLimits: ");
        			for (int i=0;i<localMax.length;i++) System.out.print(localMax[i]?"+":".");
        			System.out.println("");
        		}

    		}

    		int numMax=0;
    		for (int i=0;i<localMax.length;i++) if (localMax[i]) numMax++;
    		int [] iShifts=new int [numMax];
    		int index=0;
    		for (int i=0;i<localMax.length;i++) if (localMax[i]) iShifts[index++]=i;
    		for (boolean ordered=false;!ordered; ) {
				ordered=true;
				for (int i=0; i<(numMax-1);i++) {
					if (prevCorrelation[iShifts[i]]<prevCorrelation[iShifts[i+1]]){
						ordered=false;
						int tmp=iShifts[i];
						iShifts[i]=iShifts[i+1];
						iShifts[i+1]=tmp;
					}
				}
			}
    		// mark around old shifts
    		boolean [] alreadyDone=new boolean [prevCorrelation.length];
    		for (int i=0;i<prevCorrelation.length;i++) alreadyDone[i]=false;
    		if (currentShifts!=null) {
    			for (int i=0;i<currentShifts.length;i++) {
    				for (int j=currentShifts[i]-this.minDisparity-step+1;j<(currentShifts[i]-this.minDisparity+step);j++){
    					if ((j>=0) && (j<prevCorrelation.length)) alreadyDone[j]=true;
    				}
    			}
    		}
    		if (debugLevel>2){
    			System.out.print("generateShiftsList(): alreadyDone from old=");
    			for (int i=0;i<alreadyDone.length;i++) System.out.print(alreadyDone[i]?"+":".");
    			System.out.println("");
    		}
    		int numShifts=0;
    		for (int i=0;i<iShifts.length;i++){
    			if (alreadyDone[iShifts[i]]) {
    				iShifts[i]=-1;
    			} else {
    				for (int j=iShifts[i]-step+1;j<(iShifts[i]+step);j++) if ((j>=0) && (j<prevCorrelation.length)) alreadyDone[j]=true;
    				numShifts++;
    			}
    		}
    		int [] shifts = new int [numShifts];
    		index=0;
    		for (int i=0;i<iShifts.length;i++) if (iShifts[i]>=0){
    			shifts[index++]=iShifts[i]+this.minDisparity;
    		}
    		return shifts; // may be zero length
    	}

    	/**
    	 * Correlate two images with different shifts (needed when shifts are large compared to FFT size)
    	 * and accumulate results into an array, starting with specified  this.minDisparity
    	 * @param side - 0 - left image, 1 - right one
    	 * @param self - 0 - cross-correlate, 1 - auto-correlate
    	 * @param size - Tile size (half of the FFT size
    	 * @param corrXC Center of the correlation area on the first image, X
    	 * @param corrYC Center of the correlation area on the first image, Y
    	 * @param shifts Array of shifts to use for correlation (>=0)
    	 * @param phaseCorrelationFraction 1.0 - use pure phase correlation, 0.0 - just correlation
    	 * @param corrCbWeight relative weight of the Cb component in correlation (Y has weight of 1.0)
    	 * @param corrCrWeight relative weight of the Cr component in correlation (Y has weight of 1.0)
    	 * @param correlationHighPassSigma Correlation high-pass filter sigma,  pixels in frequency domain
    	 * @param correlationLowPassSigma Correlation low-pass filter sigma,  fraction of the frequency range
    	 * @param noiseNormalizationSignaY Noise normalization sigma for Y-component
    	 * @param noiseNormalizationSignaCbCr Noise normalization sigma for Cb, Cr components
    	 * @param doubleFHT - DoubleFHT instance or null (reusing the same instance will use cache for sin/cos/window arrays)
    	 * @param binaryAlpha Consider all points with alpha>0 to have same alpha==1.0
    	 * @param minTileFraction Minimal tile overlap fraction to process tile (will be squared in the corners of the tile map)
    	 * @param threadsMax Maximal number of concurrent threads to run
    	 * @param debugLevel Debug level
    	 * @return accumulated Correlation results for the specified range of disparities as [0][], [1][] = number of samples averaged
    	 */

    	public double [][] correlateMultipleShifts(
    			int side,
    			int self,
    			int size,
        		int corrXC,
        		int corrYC,
    			int [] shifts,    // >=0;
				double phaseCorrelationFraction,
				double corrCbWeight,
				double corrCrWeight,
				double correlationHighPassSigma,
				double correlationLowPassSigma,
				double noiseNormalizationSignaY,
				double noiseNormalizationSignaCbCr,
				boolean binaryAlpha,
				double minTileFraction, // decrease for corners
				DoubleFHT doubleFHT,
				int threadsMax,
				int debugLevel){
    		if (debugLevel>2){
    			String dbgStr="correlateMultipleShifts("+side+","+size+","+ corrXC+","+corrYC+",{";
    			for (int nShft=0; nShft<shifts.length;nShft++){
    				dbgStr+=((nShft>0)?", ":"")+shifts[nShft];
    			}
    			dbgStr+="},...)";
    			System.out.println(dbgStr);
    		}
    		int corrFFTSize=2*size;
    		int numLayers=this.overlapImages.length/2;
    		int length=corrFFTSize*corrFFTSize;
    		if (doubleFHT==null) doubleFHT = new DoubleFHT();
    		double[] hamming1d=doubleFHT.getHamming1d(corrFFTSize);
    		double [] componentWeights={
    				1.0/(1.0+corrCbWeight+corrCrWeight),
    				corrCbWeight/(1.0+corrCbWeight+corrCrWeight),
    				corrCrWeight/(1.0+corrCbWeight+corrCrWeight)};
    		double [][] accumulatedCorrelation=new double[2][this.maxDisparity-this.minDisparity+1];
    		for (int n=0;n<accumulatedCorrelation.length;n++) for (int i=0;i<accumulatedCorrelation[0].length;i++) {
    			accumulatedCorrelation[n][i]=0.0;
    		}
    		int sign=(side>0)?-1:1;
    		int selectionType=(self!=0)?(side+1):0; // 0 - normal correlation, 1 - left auto, 2 - right auto
    		for (int nShift=0;nShift<shifts.length;nShift++){
        		double [][] selection = getSelection(
        				selectionType,
        				corrFFTSize,// twice larger than tile
        				corrXC,
        				corrYC,
        				sign*shifts[nShift]); // positive - shift second image, negative - shift first image
        		double [] stats= selectionStats(
        				selection,
        				binaryAlpha);
        		if (debugLevel>2){
        			System.out.println("correlateMultipleShifts(): Selection ("+(binaryAlpha?"binary":"original")+" alpha) xc="+corrXC+" yc="+corrYC+" corrFFTSize="+corrFFTSize+" shift="+shifts[nShift]+
        					" : first "+IJ.d2s(100*stats[0],3)+"%, second "+IJ.d2s(100*stats[1],3)+"%, overlap "+IJ.d2s(100*stats[2],3)+"%");
        		}
        		if (stats[2]<minTileFraction){
            		if (debugLevel>2){
            			System.out.println("Too little overlap - "+IJ.d2s(100*stats[2],3)+"%, abandoning tile XC="+corrXC+" YC="+corrXC);
            		}
        			continue;
        		}
        		double [][] window={selection[0].clone(),selection[numLayers].clone()};
        		for (int iy=0;iy<corrFFTSize;iy++) for (int ix=0;ix<corrFFTSize;ix++){
        			int index=iy*corrFFTSize+ix;
        			double h=hamming1d[iy]*hamming1d[ix];
        			window[0][index]*=h;
        			window[1][index]*=h;
        		}
        		int map0=shifts[nShift]-size; // actual disparity of the leftmost point in the correlation result (size=1/2 FFT square)
        		int index0=size*size*2; // index of the leftmost point, will be modified to be the first correlation point to accumulate
        		int resultIndex=map0-this.minDisparity;
        		int resultNumber=2*size; // number of result pixels to calculate;
        		if  (resultIndex<0){
        			index0+=-resultIndex; // increase
        			resultIndex=0;
        			resultNumber+=resultIndex; // decrease
        		}
        		if ((resultIndex+resultNumber)>accumulatedCorrelation[0].length) resultNumber=accumulatedCorrelation[0].length-resultIndex;
        		if (resultNumber<=0){
        			if (debugLevel>1){
        				System.out.println("Correlation results do not intersect with the accumulation array");
        			}
        			continue;
        		}
				if (debugLevel>2){
					System.out.println("resultNumber="+resultNumber+" resultIndex="+resultIndex+" index0="+index0);
				}

        		double [][] corr=new double[numLayers-1][length];
        		for (int i=0;i<numLayers-1;i++){
        			double [] first= selection[i+1].clone();
        			double [] second=selection[i+1+numLayers].clone();
        			normalizeAndWindow(first, window[0],true);
        			normalizeAndWindow(second,window[1],true);
        			corr[i]= doubleFHT.correlate (
        					first,
        					second,
    						 correlationHighPassSigma,
    						 correlationLowPassSigma,
    						 phaseCorrelationFraction);
        			double sigma=(i==0)?noiseNormalizationSignaY:noiseNormalizationSignaCbCr;
        			if (sigma>0){
        				double rms=hFNoiseNormalize(
        						corr[i], // double [] data,
        		    			sigma, // double filterSigma,
        		    			true);  // boolean centerOnly);
        				if (this.debugLevel>2){
        					System.out.println("Correlation RMS for component "+((i==0)?"Y":((i==1)?"Cb":"Cr"))+ " was "+rms);
        				}
        			}
        			// may only calculate for the middle line
        			for (int j=0;j<resultNumber;j++) {
        				accumulatedCorrelation[0][resultIndex+j]+=componentWeights[i]*corr[i][index0+j];
        				accumulatedCorrelation[1][resultIndex+j]+=1.0;
        			}
        		}

    		}
    		int numPoints=0;
    		for (int i=0; i<accumulatedCorrelation[0].length;i++){
    			if (accumulatedCorrelation[1][i]>0){
    			  accumulatedCorrelation[0][i]/=accumulatedCorrelation[1][i];
    			  numPoints++;
    			}
    		}

    		if (numPoints==0){
        		if (debugLevel>1){
        			System.out.println("Nothing accumulated");
        		}
    			return null;
    		}
    		return accumulatedCorrelation;
    	}

// old version (*0), delete
    	public double [][] collectDisparityFromTiles0(
    			int debugLevel){
    		double [][] disparity=new double [6][this.mapWidth*this.mapHeight];
    		for (int n=0;n<disparity.length;n++) for (int i=0;i<disparity[0].length;i++)disparity[n][i]=Double.NaN;
    		int yMinLim=this.tilesY0;
    		int xMinLim=this.tilesX0;
    		int yMaxLim=this.tilesY0+this.tilesSize*this.tilesNVert;
    		int xMaxLim=this.tilesX0+this.tilesSize*this.tilesNHor;
    		if (yMinLim<0) yMinLim=0;
    		if (xMinLim<0) xMinLim=0;
    		if (yMaxLim>this.mapHeight) yMaxLim=this.mapHeight;
    		if (xMaxLim>this.mapWidth)  xMaxLim=this.mapWidth;
    		if (debugLevel>1) System.out.println("collectDisparityFromTiles(): this.tilesSize="+this.tilesSize+",yMinLim="+yMinLim+
    				", yMaxLim="+yMaxLim+", xMinLim="+xMinLim+", xMaxLim="+xMaxLim);
    		for (int iY=yMinLim;iY<yMaxLim;iY++){
    			int tY=iY-this.tilesY0;
    			for (int iX=xMinLim;iX<xMaxLim;iX++){
    				int tX=iX-this.tilesX0;
    				for (int level=this.tiles0.length-1;level>=0;level--){

    					int size = (this.tilesSize>>level);
						double [][] tile=this.tiles0[level][tY/size][tX/size];
    					if (tile!=null){
    						int index=iY*this.mapWidth+iX;
    						disparity[2][index]=tile[0][0];
    						disparity[3][index]=tile[0][1]; // contrast background
    						disparity[4][index]=tile[tile.length-1][0]; // last - FG if available, if not - BG
    						disparity[5][index]=tile[tile.length-1][1]; // contrast *-ground
    						disparity[0][index]=(disparity[3][index]>disparity[5][index])?disparity[2][index]:disparity[4][index];
    						disparity[1][index]=(disparity[3][index]>disparity[5][index])?disparity[3][index]:disparity[5][index];
    						break;
    					}
    				}
    			}
    		}
    		return disparity;
    	}

    	public int [] generateTileLevels0(
    			double phaseCorrelationFraction,
				double corrCbWeight,
				double corrCrWeight,
    			double correlationHighPassSigma,
    			double correlationLowPassSigma,
    			double noiseNormalizationSignaY,
    			double noiseNormalizationSignaCbCr,
    			double contrastThreshold,
    			double contrastThresholdDecrease, //Decrease threshold for smaller tiles as 1/2^(k*level), k=
    			boolean binaryAlpha,
				boolean enableNegativeDisparity,
    			double minTileFraction,
    			DoubleFHT doubleFHT,
    			int threadsMax,
    			int debugLevel){
    		if (doubleFHT==null) doubleFHT=new DoubleFHT();
    		int numNonEmptyTiles=0;
    		int numWithForeground=0;
    		for (int level=0;level<this.tiles0.length;level++){
    			double contastThresholdMod=contrastThreshold*Math.pow(0.5,contrastThresholdDecrease*level);
    			int [] tileResult=generateTileLevel0(
    					level,
    					phaseCorrelationFraction,
    					corrCbWeight,
    					corrCrWeight,
    					correlationHighPassSigma,
    					correlationLowPassSigma,
    					noiseNormalizationSignaY,
    					noiseNormalizationSignaCbCr,
    					contastThresholdMod,
    					binaryAlpha,
    					enableNegativeDisparity,
    					minTileFraction,
    					doubleFHT,
    					threadsMax,
    					debugLevel);
    			if (debugLevel>1){
    				System.out.println("On level "+level+": non-empty tiles="+tileResult[0]+" tiles with foreground="+tileResult[1]);
    			}
    			numNonEmptyTiles+= tileResult[0];
    			numWithForeground+=tileResult[1];
    		}
    		int [] result={numNonEmptyTiles,numWithForeground};
    		return result;
    	}

// Convert to multithread
    	public int [] generateTileLevel0(
    			int level,
				double phaseCorrelationFraction,
				double corrCbWeight,
				double corrCrWeight,
				double correlationHighPassSigma,
				double correlationLowPassSigma,
				double noiseNormalizationSignaY,
				double noiseNormalizationSignaCbCr,
				double contrastThreshold,
				boolean binaryAlpha,
				boolean enableNegativeDisparity,
				double minTileFraction,
				DoubleFHT doubleFHT,
				int threadsMax,
				int debugLevel){
    		int numNonEmptyTiles=0;
    		int numWithForeground=0;
    		for (int tileY=0;tileY<this.tiles0[level].length;tileY++) for (int tileX=0;tileX<this.tiles0[level][0].length;tileX++){
    			double [][] tileResult=generateTile0(
    	    			level,
    	    			tileX,
    	    			tileY,
    					phaseCorrelationFraction,
    					corrCbWeight,
    					corrCrWeight,
    					correlationHighPassSigma,
    					correlationLowPassSigma,
    					noiseNormalizationSignaY,
    					noiseNormalizationSignaCbCr,
    					contrastThreshold,
    					binaryAlpha,
    					enableNegativeDisparity,
    					minTileFraction,
    					doubleFHT,
    					threadsMax,
    					debugLevel);
    			if (tileResult!=null) {
    				numNonEmptyTiles++;
    				if (tileResult.length>1){
    					numWithForeground++;
    				}
    			}
    		}
    		int [] result={numNonEmptyTiles,numWithForeground};
    		return result;
    	}


    	/**
    	 * Generate disparity tile by correlating the images, starting with info from the higher level tile (if available)
    	 * @param level level of the tile to generate
    	 * @param tileX Horizontal tile number
    	 * @param tileY Vertical tile number
    	 * @param phaseCorrelationFraction 1.0 - use pure phase correlation, 0.0 - just correlation
    	 * @param correlationHighPassSigma Correlation high-pass filter sigma,  pixels in frequency domain
    	 * @param correlationLowPassSigma Correlation low-pass filter sigma,  fraction of the frequency range
    	 * @param noiseNormalizationSignaY Noise normalization sigma for Y-component
    	 * @param noiseNormalizationSignaCbCr Noise normalization sigma for Cb, Cr components
    	 * @param contrastThreshold Use correlation data with normalized values above this threshold only
    	 * @param binaryAlpha consider all points with alpha>0 to have same alpha==1.0
    	 * @param minTileFraction minimal tile overlap fraction to process tile (will be squared in the corners of the tile map)
    	 * @param threadsMax Maximal number of concurrent threads to run
    	 * @param debugLevel Debug level
    	 * @return Generated tile (or null if no tile was generated)
    	 */
    	public double [][] generateTile0(
    			int level,
    			int tileX,
    			int tileY,
				double phaseCorrelationFraction,
				double corrCbWeight,
				double corrCrWeight,
				double correlationHighPassSigma,
				double correlationLowPassSigma,
				double noiseNormalizationSignaY,
				double noiseNormalizationSignaCbCr,
				double contrastThreshold,
				boolean binaryAlpha,
				boolean enableNegativeDisparity,
				double minTileFraction,
				DoubleFHT doubleFHT,
				int threadsMax,
				int debugLevel){
    		int tileDataLength=5; // only 2 used?
    		int size=(this.tilesSize>>level);
    		int corrFFTSize=2*size;
    		int numLayers=this.overlapImages.length/2;
    		boolean isCorner=(
    				((tileY==0)||(tileY==(this.tiles0[level].length))) &&
    				((tileX==0)||(tileX==(this.tiles0[level][0].length))));
    		if (isCorner) minTileFraction*=minTileFraction;
    		int corrXC=this.tilesX0+tileX*size+size/2;
    		int corrYC=this.tilesY0+tileY*size+size/2;
    		double [][] parentTile;
    		if ( level==0){
    			parentTile=new double[1][];
    			// {position,strength,rmsY,rmsCb,rmsCr}
    			parentTile[0]=new double [tileDataLength];
    			parentTile[0][0]=0; // center position, all the rest are not needed
    		} else {
    			parentTile= this.tiles0[level-1][tileY/2][tileX/2];
    		}
    		if (parentTile==null) return null; // no parent to this tile
    		double [][] corrResults ={null,null};
    		int [] shifts={-1,-1}; // undefined
    		double[] hamming1d=doubleFHT.getHamming1d(corrFFTSize);
    		double [] componentWeights={
    				1.0/(1.0+corrCbWeight+corrCrWeight),
    				corrCbWeight/(1.0+corrCbWeight+corrCrWeight),
    				corrCrWeight/(1.0+corrCbWeight+corrCrWeight)};


    		int length=corrFFTSize*corrFFTSize;
    		for (int numBgFg=0;numBgFg<parentTile.length;numBgFg++) if (parentTile[numBgFg]!=null){ // try for both for- and backgrounds
    			shifts[numBgFg]=(int) Math.round(parentTile[numBgFg][0]);
        		double [][] selection = getSelection(
        				corrFFTSize,// twice larger than tile
        				corrXC,
        				corrYC,
        				shifts[numBgFg]); // positive - shift second image, negative - shift first image
        		double [] stats= selectionStats(
        				selection,
        				binaryAlpha);
        		if (debugLevel>1){
        			System.out.println("Selection ("+(binaryAlpha?"binary":"original")+" alpha) xc="+corrXC+" yc="+corrYC+" corrFFTSize="+corrFFTSize+" shift="+shifts[numBgFg]+
        					" : first "+IJ.d2s(100*stats[0],3)+"%, second "+IJ.d2s(100*stats[1],3)+"%, overlap "+IJ.d2s(100*stats[2],3)+"%");
        		}
        		if (stats[2]<minTileFraction){
            		if (debugLevel>1){
            			System.out.println("Too little overlap - "+IJ.d2s(100*stats[2],3)+"%, abandoning tile XC="+corrXC+" YC="+corrXC);
            		}
        			continue;
        		}
        		if (doubleFHT==null) doubleFHT = new DoubleFHT();
        		double [][] window={selection[0].clone(),selection[numLayers].clone()};
        		for (int iy=0;iy<corrFFTSize;iy++) for (int ix=0;ix<corrFFTSize;ix++){
        			int index=iy*corrFFTSize+ix;
        			double h=hamming1d[iy]*hamming1d[ix];
        			window[0][index]*=h;
        			window[1][index]*=h;
        		}
        		double [][] corr=new double[numLayers][length];
        		for (int i=0;i<length;i++) corr[0][i]=0;
        		for (int i=0;i<numLayers-1;i++){
        			double [] first= selection[i+1].clone();
        			double [] second=selection[i+1+numLayers].clone();
        			normalizeAndWindow(first, window[0],true);
        			normalizeAndWindow(second,window[1],true);
        			corr[i+1]= doubleFHT.correlate (
        					first,
        					second,
    						 correlationHighPassSigma,
    						 correlationLowPassSigma,
    						 phaseCorrelationFraction);
        			double sigma=(i==0)?noiseNormalizationSignaY:noiseNormalizationSignaCbCr;
        			if (sigma>0){
        				double rms=hFNoiseNormalize(
        						corr[i+1], // double [] data,
        		    			sigma, // double filterSigma,
        		    			true);  // boolean centerOnly);
        				if (this.debugLevel>2){
        					System.out.println("Correlation RMS for component "+((i==0)?"Y":((i==1)?"Cb":"Cr"))+ " was "+rms);
        				}
        			}
            		for (int j=0;j<length;j++) corr[0][j]+=componentWeights[i]*corr[i+1][j];

        		}
        		corrResults[numBgFg]=correlationResults(
            			corr[0], // composite - Y, Cb, Cr
            			contrastThreshold,
            			shifts[numBgFg],
            			enableNegativeDisparity,
            			true, // boolean orderMaximums,
            			debugLevel);
				if (this.debugLevel>1){
					String dbgStr="corrResults["+numBgFg+"]";
					if (corrResults[numBgFg]==null) dbgStr+="=null";
					else {
						dbgStr+=": min="+IJ.d2s(corrResults[numBgFg][0],3)+
						", max="+IJ.d2s(corrResults[numBgFg][1],3)+
						", MAX1="+IJ.d2s(corrResults[numBgFg][2],3)+" ("+IJ.d2s(corrResults[numBgFg][3],3)+")";
						for (int j=2;j<(corrResults[numBgFg].length/2);j++){
							dbgStr+=", MAX"+j+"="+IJ.d2s(corrResults[numBgFg][2*j],3)+" ("+IJ.d2s(corrResults[numBgFg][2*j+1],3)+")";
						}
					}
					System.out.println(dbgStr);
				}

    		}
    		if ((corrResults[0]==null) && (corrResults[1]==null)) return null; // both are nulls
    		// combine results from both FG and BG
    		int minShift=this.tilesSize,maxShift=0;
    		for (int numBgFg=0;numBgFg<2;numBgFg++) if (corrResults[numBgFg]!=null){
    			for (int nMax=0;(nMax<2) && (nMax<(corrResults[numBgFg].length/2-1));nMax++){
    				int shift=(int) Math.round(corrResults[numBgFg][nMax*2+2]);
    				if (shift<minShift) minShift=shift;
    				if (shift>maxShift) maxShift=shift;
    			}
    		}
    		int [] oldShifts=shifts.clone(); // element = -1 - never tried
    		shifts=new int [(maxShift>minShift)?2:1];
    		shifts[0]=minShift;
    		shifts[shifts.length-1]=maxShift; // OK if minShift==maxShift
// use corrResults if the shifts are the same
    		double [][] oldResults=corrResults.clone();
    		corrResults=new double [shifts.length][];
    		boolean [] newShift={false,false};
    		for (int i=0;i<corrResults.length;i++){
    			if      (shifts[i]==oldShifts[0]) corrResults[i]=oldResults[0];
    			else if (shifts[i]==oldShifts[1]) corrResults[i]=oldResults[1];
    			else {
    				corrResults[i]=null;
    				newShift[i]=true;
    			}
    		}

    		// re-run correlation with selections centered at target shifts[i]
    		for (int numShift=0;numShift<shifts.length;numShift++) if (newShift[numShift]){
        		double [][] selection = getSelection(
        				corrFFTSize,// twice larger than tile
        				corrXC,
        				corrYC,
        				shifts[numShift]); // positive - shift second image, negative - shift first image
        		double [] stats= selectionStats(
        				selection,
        				binaryAlpha);
        		if (debugLevel>1){
        			System.out.println("Selection-second pass ("+(binaryAlpha?"binary":"original")+" alpha) xc="+corrXC+" yc="+corrYC+" corrFFTSize="+corrFFTSize+" shift="+shifts[numShift]+
        					" : first "+IJ.d2s(100*stats[0],3)+"%, second "+IJ.d2s(100*stats[1],3)+"%, overlap "+IJ.d2s(100*stats[2],3)+"%");
        		}
        		if (stats[2]<minTileFraction){
            		if (debugLevel>1){
            			System.out.println("Too little overlap - "+IJ.d2s(100*stats[2],3)+"%, abandoning tile XC="+corrXC+" YC="+corrXC);
            		}
        			continue;
        		}
        		double [][] window={selection[0].clone(),selection[numLayers].clone()};
        		for (int iy=0;iy<corrFFTSize;iy++) for (int ix=0;ix<corrFFTSize;ix++){
        			int index=iy*corrFFTSize+ix;
        			double h=hamming1d[iy]*hamming1d[ix];
        			window[0][index]*=h;
        			window[1][index]*=h;
        		}
        		double [][] corr=new double[numLayers][length];
        		for (int i=0;i<length;i++) corr[0][i]=0;
        		for (int i=0;i<numLayers-1;i++){
        			double [] first= selection[i+1].clone();
        			double [] second=selection[i+1+numLayers].clone();
        			normalizeAndWindow(first, window[0],true);
        			normalizeAndWindow(second,window[1],true);
        			corr[i+1]= doubleFHT.correlate (
        					first,
        					second,
    						 correlationHighPassSigma,
    						 correlationLowPassSigma,
    						 phaseCorrelationFraction);
        			double sigma=(i==0)?noiseNormalizationSignaY:noiseNormalizationSignaCbCr;
        			if (sigma>0){
        				double rms=hFNoiseNormalize(
        						corr[i+1], // double [] data,
        		    			sigma, // double filterSigma,
        		    			true);  // boolean centerOnly);
        				if (this.debugLevel>2){
        					System.out.println("--Correlation RMS for component "+((i==0)?"Y":((i==1)?"Cb":"Cr"))+ " was "+rms);
        				}
        			}
//            		for (int j=0;j<length;j++) corr[0][j]+=corr[i][j];
            		for (int j=0;j<length;j++) corr[0][j]+=componentWeights[i]*corr[i+1][j];

        		}
        		corrResults[numShift]=correlationResults(
            			corr[0], // composite - Y, Cb, Cr
            			contrastThreshold,
            			shifts[numShift],
            			enableNegativeDisparity,
            			true, // boolean orderMaximums,
            			debugLevel);
				if (this.debugLevel>1){
					String dbgStr="corrResults["+numShift+"]";
					if (corrResults[numShift]==null) dbgStr+="=null";
					else {
						dbgStr+=": min="+IJ.d2s(corrResults[numShift][0],3)+
						", max="+IJ.d2s(corrResults[numShift][1],3)+
						", MAX1="+IJ.d2s(corrResults[numShift][2],3)+" ("+IJ.d2s(corrResults[numShift][3],3)+")";
						for (int j=2;j<(corrResults[numShift].length/2);j++){
							dbgStr+=", MAX"+j+"="+IJ.d2s(corrResults[numShift][2*j],3)+" ("+IJ.d2s(corrResults[numShift][2*j+1],3)+")";
						}
					}
					System.out.println(dbgStr);
				}



    		}
  // now one or 2 correlations were run with shifts[] , each of the results is expected to be close to those (that)
    		double [][] rslt = {null,null};

    		if ((corrResults[0]==null) && ((corrResults.length==1) || (corrResults[1]==null))) return null;

    		if ((corrResults.length>1) && (corrResults[0]!=null) && (corrResults[1]!=null)){ // two results
    			for (int i=0;i<rslt.length;i++) {
    				rslt[i]=new double [2];
					boolean useFirst=(corrResults[i].length<=4) ||
					(Math.abs(corrResults[i][2]-shifts[i])<Math.abs(corrResults[i][4]-shifts[i])); // first result is closer to selection shift
    					rslt[i][0]=corrResults[i][useFirst?2:4];
    					rslt[i][1]=corrResults[i][useFirst?3:5];
    			}
    		} else {
    			int i=(corrResults[0]==null)?1:0;
    			rslt[0]=new double[2];
    			rslt[0][0]=corrResults[i][2];
    			rslt[0][1]=corrResults[i][3];
    			if (corrResults[i].length>4){
        			rslt[1]=new double[2];
        			rslt[1][0]=corrResults[i][4];
        			rslt[1][1]=corrResults[i][5];
    			}
    		}
    		// reorder to have background first
    		if ((rslt[1]!=null) && (rslt[1][0]<rslt[0][0])){
    			double [] tmp=rslt[1];
    			rslt[1]=rslt[0];
    			rslt[0]=tmp;
    		}
    		// combine 2 close ones into one
    		if ((rslt[1]!=null) && (Math.abs(rslt[1][0]-rslt[0][0])<1.0)){
    			rslt[0][0]=0.5*(rslt[0][0]+rslt[1][0]);
    			rslt[0][1]=0.5*(rslt[0][1]+rslt[1][1]);
    			rslt[1]=null;
    		}
    		// save result tile
    		this.tiles0[level][tileY][tileX]=new double[(rslt[1]!=null)?2:1][];
    		for (int i=0;i<this.tiles0[level][tileY][tileX].length;i++){
    			this.tiles0[level][tileY][tileX][i]=new double[tileDataLength];
    			this.tiles0[level][tileY][tileX][i][0]=rslt[i][0];
    			this.tiles0[level][tileY][tileX][i][1]=rslt[i][1];
    		}
    		return this.tiles0[level][tileY][tileX];
    	}




    	/**
    	 * Setup coverage of the overlap area with tiles of the specified size
    	 * @param size Tile size (half of the windowed FFT size)
    	 * @levels number of levels of progressively smaller tiles
    	 * @return number of tiles horizontally and vertically
    	 */
    	public int [] setupTiles0(
    			int size,
    			int levels){ // number of tile levels, each having half linear size (4..5)
    		int numLayers=this.overlapImages.length/2;

    		int minY=this.mapHeight,minX=this.mapWidth,maxY=0,maxX=0;

    		for (int iy=0;iy<this.mapHeight;iy++) for (int ix=0;ix<this.mapWidth;ix++){
    			int index=iy*this.mapWidth+ix;
    			if ((this.overlapImages[0][index]>0) && (this.overlapImages[numLayers][index]>0)){
    				if (iy>maxY) maxY=iy;
    				if (iy<minY) minY=iy;
    				if (ix>maxX) maxX=ix;
    				if (ix<minX) minX=ix;
    			}
    		}
    		if (maxY<minY) {
    			if (this.debugLevel>0) System.out.println("No overlap");
    			return null;
    		}
    		int width= (maxX-minX+1);
    		int height=(maxY-minY+1);
    		this.tilesNHor=width/size;
    		if ((this.tilesNHor*size)<width) this.tilesNHor++;
    		this.tilesNVert=height/size;
    		if ((this.tilesNVert*size)<height) this.tilesNVert++;
    		// center tiles over overlap area
    		this.tilesX0=(maxX+minX+1)/2-(this.tilesNHor*size)/2;
    		this.tilesY0=(maxY+minY+1)/2-(this.tilesNVert*size)/2;
    		int [] nTiles={this.tilesNHor,this.tilesNVert};
    		this.tilesSize=size;
    		this.tiles0 = new double[levels][][][][];
    		for (int nLev=0;nLev<levels;nLev++){
    			this.tiles0[nLev]=new double [this.tilesNVert<<nLev][this.tilesNHor<<nLev][][];
    			for (int tileY=0;tileY<this.tiles0[nLev].length;tileY++)
    				for (int tileX=0;tileX<this.tiles0[nLev][0].length;tileX++) this.tiles0[nLev][tileY][tileX]=null;
    		}
    		return nTiles;
    	}


    	/**
    	 * Get a square selection (for FFT-based correlation) from a pair of A-Y-Cb-Cr images. OK if selection extends beyond source image
    	 * @param size size of square side
    	 * @param xc selection center X
    	 * @param yc selection center Y (down)
    	 * @param shiftX if positive - shift second image right (selection - left), if negative - shift the first on left (selection - right)
    	 * @return array with half layers of the source (one set of A-Y-Cb-Cr), packed in line-scan order [size*size]
    	 */
    	public double [][] getSelection(
    			int size,
    			int xc,
    			int yc,
    			int shiftX // positive - shift second image, negative - shift first image
    			){
    		return getSelection(0,size,xc,yc,shiftX);
    	}

    	public double [][] getSelection(
    			int pairType, // 0 - norm, 1 - use first image for both, 2 - use second image for both
    			int size,
    			int xc,
    			int yc,
    			int shiftX // positive - shift second image, negative - shift first image
    			){
    		return  getSelection(
        			pairType, // 0 - norm, 1 - use first image for both, 2 - use second image for both
        			size,
        			size,
        			xc,
        			yc,
        			shiftX); // positive - shift second image, negative - shift first image
    	}

    	public double [][] getSelection(
    			int pairType, // 0 - norm, 1 - use first image for both, 2 - use second image for both // New - only applicable to 2 images
    			int width,
    			int height,
    			int xc,
    			int yc,
    			int shiftX // positive - shift second image, negative - shift first image
    			){
    		int [] x0={
    				xc-width/2-((shiftX<0)?shiftX:0),
    				xc-width/2-((shiftX>0)?shiftX:0)};
    		int y0=yc-height/2;
    		int length=width*height;
    		int numSensors=this.channel.length;
//    		double [][] selection = new double[this.overlapImages.length][length];
    		double [][] selection = new double[this.overlapImages.length][];
    		for (int i=0;i<selection.length;i++) {
    			if (this.overlapImages[i]!=null)	selection[i]=new double[length];
    			else selection[i]=null;
    		}
    		int numLayers=this.overlapImages.length/numSensors;
			x0=new int [numSensors];
			for (int i=0;i<numSensors;i++){
				switch (i) {
				case 0: x0[i]=xc-width/2-((shiftX<0)?shiftX:0); break;
				case 1: x0[i]=xc-width/2-((shiftX>0)?shiftX:0); break;
				default:x0[i]=xc-width/2;
				}
			}

    		if ((numSensors!=2) || (this.overlapImages[0]==null)  || (this.overlapImages[0]==null)) {
    			pairType=0;
    		}
    		for (int iy=0;iy<height;iy++) {
    			int srcY=iy+y0;
    			boolean oob=(srcY<0) || (srcY>=this.mapHeight);
    			for (int ix=0;ix<width;ix++){
    				for (int numImg=0;numImg<numSensors;numImg++) if (this.overlapImages[numImg*numLayers]!=null){
//    					System.out.print(iy+":"+ix);
    					int modNumImg=numImg;
    					if (pairType>0){
    						modNumImg=(pairType==1)?0:1;
    					}
        				int oIndex=iy*width+ix;
    					int srcX=x0[numImg]+ix;
//    					int numNonZero=0;
    					if (oob ||(srcX<0) || (srcX>=this.mapWidth)) {
    						for (int i=0;i<numLayers;i++) if (selection[i+numImg*numLayers]!=null)  selection[i+numImg*numLayers][oIndex]=0.0;
    					} else {
    						int iIndex=srcY*this.mapWidth+srcX;
        					for (int i=0;i<numLayers;i++) if (selection[i+numImg*numLayers]!=null) {
        						selection[i+numImg*numLayers][oIndex]=this.overlapImages[i+modNumImg*numLayers][iIndex];
//        						if ((selection[i+numImg*numLayers][oIndex]!=0.0) && (i>0)){
//            						System.out.print(iy+":"+ix+" i="+i+" i+modNumImg*numLayers="+(i+modNumImg*numLayers)+" ");
//        							numNonZero++;
//        						}
        					}
    					}
//    					if (numNonZero>1) {
//    						System.out.print(iy+":"+ix);
//    						System.out.print("="+numNonZero+"   ");
//    					}
    				}
    			}
    		}
    		return selection;
    	}
/*
    	public double [] getSelection(
				int width,
				int height,
				int numImg,
				int chn,
				int xc,
				int yc){
    		int length=width*height;
    		int numSensors=this.channel.length;
    		double [] selection = new double[length];
    		int numLayers=this.overlapImages.length/numSensors;
				int y0=yc-height/2;
				int x0=xc-width/2;
	    		for (int iy=0;iy<height;iy++) {
	    			int srcY=iy+y0;
	    			boolean oob=(srcY<0) || (srcY>=this.mapHeight);
	    			for (int ix=0;ix<width;ix++){
        				int oIndex=iy*width+ix;
    					int srcX=x0+ix;
    					if (oob ||(srcX<0) || (srcX>=this.mapWidth)) {
    						selection[oIndex]=0.0;
    					} else {
        						selection[oIndex]=this.overlapImages[chn+numImg*numLayers][srcY*this.mapWidth+srcX];
    					}
	    			}
	    		}
    		return selection;
    	}
  */
    	public double [] getSelection(
    			double [] imageSlice, // one image/channel slice
    			double [] selection, // null or array to reuse
				int width,
				int height,
				int xc,
				int yc){
    		int length=width*height;
    		if (selection ==null) selection = new double[length];
				int y0=yc-height/2;
				int x0=xc-width/2;
	    		for (int iy=0;iy<height;iy++) {
	    			int srcY=iy+y0;
	    			boolean oob=(srcY<0) || (srcY>=this.mapHeight);
	    			for (int ix=0;ix<width;ix++){
        				int oIndex=iy*width+ix;
    					int srcX=x0+ix;
    					if (oob ||(srcX<0) || (srcX>=this.mapWidth)) {
    						selection[oIndex]=0.0;
    					} else {
        						selection[oIndex]=imageSlice[srcY*this.mapWidth+srcX];
    					}
	    			}
	    		}
    		return selection;
    	}

    	public double [][] getSelection(
				int width,
				int height,
				int [][] iCenterXY){
    		int length=width*height;
    		int numSensors=this.channel.length;
    		double [][] selection = new double[this.overlapImages.length][];
    		for (int i=0;i<selection.length;i++) {
    			if (this.overlapImages[i]!=null)	selection[i]=new double[length];
    			else selection[i]=null;
    		}
    		int numLayers=this.overlapImages.length/numSensors;
			for (int numImg=0;numImg<numSensors;numImg++) if (this.overlapImages[numImg*numLayers]!=null){
				int y0=iCenterXY[numImg][1]-height/2;
				int x0=iCenterXY[numImg][0]-width/2;
	    		for (int iy=0;iy<height;iy++) {
	    			int srcY=iy+y0;
	    			boolean oob=(srcY<0) || (srcY>=this.mapHeight);
	    			for (int ix=0;ix<width;ix++){
        				int oIndex=iy*width+ix;
    					int srcX=x0+ix;
    					if (oob ||(srcX<0) || (srcX>=this.mapWidth)) {
    						for (int i=0;i<numLayers;i++) selection[i+numImg*numLayers][oIndex]=0.0;
    					} else {
    						int iIndex=srcY*this.mapWidth+srcX;
        					for (int i=0;i<numLayers;i++) if (selection[i+numImg*numLayers]!=null) {
        						int slice=i+numImg*numLayers;
        						selection[slice][oIndex]=this.overlapImages[slice][iIndex];
        					}
    					}
	    			}
	    		}
			}
    		return selection;
    	}



    	public double [] selectionStats(
    			double [][]selection,
    			boolean binaryAlpha){
    		int numLayers=selection.length/2;
    		double [] stats={0.0, 0.0, 0.0};
    		for (int i=0;i<selection[0].length;i++){
    			double a0=selection[0][i];
    			double a1=selection[numLayers][i];
    			if (binaryAlpha){
    				a0=(a0>0.0)?1.0:0.0;
    				a1=(a1>0.0)?1.0:0.0;
    			}
    			stats[0]+=a0;
    			stats[1]+=a1;
    			stats[2]+=a0*a1;
    		}
    		for (int i=0;i<stats.length;i++) stats[i]/=selection[0].length;
    		return stats;
    	}

    	/**
    	 * Normalize square data by dividing each element by the RMS of the high-frequency part of it
    	 * @param data square array (packed in 1d in line-scan order), will be modified
    	 * @param filterSigma filter sigma (in pixels), HF filter is subtraction of Gaussian- blurred copy
    	 * @return RMS of the original data
    	 */
    	public double hFNoiseNormalize(
    			double [] data,
    			double filterSigma,
    			boolean centerOnly){
    		int size = (int) Math.round(Math.sqrt(data.length));
    		double [] bluredData=data.clone();
    		(new DoubleGaussianBlur()).blurDouble(
    				bluredData,
					size,
					size,
					filterSigma,
					filterSigma,
					0.01);
    		int scanMin=centerOnly?(size/4):0;
    		int scanMax=centerOnly?(3*size/4):size;
    		int num=0;
    		double rms=0.0;
    		for (int iy=scanMin;iy<scanMax;iy++) for (int ix=scanMin;ix<scanMax;ix++){
    			int index=iy*size+ix;
    			double d=data[index]-bluredData[index];
    			rms+=d*d;
    			num++;
    		}
    		rms=Math.sqrt(rms/num);
    		if (rms>0){ // Prevent NaN if there is no data (i.e. oversaturated area)
    			double rrms=1.0/rms;
    			for (int i=0;i<data.length;i++){
    				data[i]*=rrms;
    			}
    		}
    		return rms;
    	}

    	public double[] normalizeAndWindow (double [] pixels, double [] windowFunction, boolean removeDC) {
    		int j;
    		if (pixels==null) return null;
    		double s=0.0,s0=0.0;
    		if (removeDC) {
    			for (j=0;j<pixels.length;j++){
    				s+=pixels[j]*windowFunction[j];
    				s0+=windowFunction[j];

    			}
    			s/=s0;
    		}
    		if (windowFunction!=null) for (j=0;j<pixels.length;j++) pixels[j]=(pixels[j]-s)*windowFunction[j];
    		else for (j=0;j<pixels.length;j++) pixels[j]=(pixels[j]-s);
    		return pixels;
    	}

    	public double normalizeAndWindowGetDC (double [] pixels, double [] windowFunction) {
    		int j;
    		if (pixels==null) return 0;
    		double s=0.0,s0=0.0;
    		for (j=0;j<pixels.length;j++){
    			s+=pixels[j]*windowFunction[j];
    			s0+=windowFunction[j];
    		}
    		s/=s0;
    		for (j=0;j<pixels.length;j++) pixels[j]=(pixels[j]-s)*windowFunction[j];
    		return s;
    	}


    	public ImagePlus resampleToOverlapRGB24(
    			boolean isOverlap, // false - justsome projection  plane image
    			ImagePlus imp,
    			int channel,
    			final int scale,
    			final double [][][][] lanczos,
    			final int binsPerHalfPixel,
    			int maxThreads){
//    		final SensorData.EquirectangularMap erm=this.sensors[channel].equirectangularMap;
    		final int width=imp.getWidth();
    		final int height=imp.getHeight();
    		if (scale==0) {
    			String msg="Image has less resolution than the map, can not proceed";
    			System.out.println(msg);
    			IJ.showMessage(msg);
    			return null;
    		}
//TODO - make for any number of channels available in the map
    		int index;
    		for (index=0;index<this.channel.length;index++) if (channel==this.channel[index]) break;

    		if (index>=this.channel.length){
    			String msg="This map is for channel"+((this.channel.length>1)?"s ":" ");
    			for (int i=0;i<this.channel.length;i++){
    				if (i>0) msg +=",";
    				msg+=" "+this.channel[i];
    			}
    			msg+=", and this image is for channel "+channel;
    			System.out.println(msg);
    			IJ.showMessage(msg);
    			return null;
    		}
//    		boolean second=(channel ==this.channel[1]);
    		final float [][] map={this.map[index*3],this.map[index*3+1],this.map[index*3+2]};
    		final int [] imagePixels=(int []) imp.getProcessor().getPixels();
    		final int mapWidth= this.mapWidth;
    		final int mapHeight=this.mapHeight;
    		final int [] outPixels=new int [mapWidth*mapHeight];
    		final int center=lanczos[0][0][0].length/2;
    		final Thread[] threads = newThreadArray(maxThreads);
    		final AtomicInteger opyAtomic     = new AtomicInteger(0);
    		final AtomicInteger opyDoneAtomic = new AtomicInteger(1);
    		//    	   		for (int opy=0;opy<erm.mapWOI.height;opy++){
    		IJ.showStatus("Warping image channel "+channel);
    		for (int ithread = 0; ithread < threads.length; ithread++) {
    			threads[ithread] = new Thread() {
    				@Override
					public void run() {
    					for (int opy=opyAtomic.getAndIncrement(); opy<mapHeight;opy=opyAtomic.getAndIncrement()){

    						for (int opx=0;opx<mapWidth;opx++){

    							int oIndex=opy*mapWidth+opx;
    							double [] RGBA={255*map[2][oIndex],0.0,0.0,0.0};
    							if (RGBA[0]>0){ // do not convolve pixels with alpha=0;
    								double x=scale*map[0][oIndex];
    								double y=scale*map[1][oIndex];
    								int ix= (int) Math.round(x);
    								int iy= (int) Math.round(y);
    								double dx=x-ix;
    								double dy=y-iy;
    								int indxX= (int) Math.round((2*dx+1)*binsPerHalfPixel);
    								int indxY= (int) Math.round((2*dy+1)*binsPerHalfPixel);
    								double [][] lk=lanczos[indxY][indxX];
    								for (int i=0;i<lk.length;i++) {
    									int ipy=iy+i-center;
    									if (ipy<0) ipy=0;
    									else if (ipy>=height) ipy=height-1;
    									int baseI=ipy*width;
    									for (int j=0;j<lk[0].length;j++){
    										int ipx=ix+j-center;
    										if (ipx<0) ipx=0;
    										else if (ipx>=width) ipx=width-1;
    										int pix=imagePixels[baseI+ipx];
    										RGBA[1]+=((pix>>16) & 0xff)*lk[i][j]; // R
    										RGBA[2]+=((pix>>8)  & 0xff)*lk[i][j]; // G
    										RGBA[3]+=( pix      & 0xff)*lk[i][j]; // B
    									}

    								}
    								int c=0;
    								for (int i=0;i<4;i++){
    									if      (RGBA[i]< 0.0)    RGBA[i]=0.0;
    									else if (RGBA[i] > 255.0) RGBA[i]=255.0;
    									c=(c<<8) | ((int) RGBA[i]);
    									outPixels [oIndex]=c;
    								}
    							} else {
    								outPixels[oIndex]=0;
    							}

    						}
    						final int numFinished=opyDoneAtomic.getAndIncrement();
    						//					IJ.showProgress(progressValues[numFinished]);
    						SwingUtilities.invokeLater(new Runnable() {
    							@Override
								public void run() {
    								// Here, we can safely update the GUI
    								// because we'll be called from the
    								// event dispatch thread
    								IJ.showProgress(numFinished,mapHeight);
    							}
    						});
    					}
    				}
    			};
    		}
    		startAndJoin(threads);
    		ColorProcessor cp=new ColorProcessor(mapWidth,mapHeight);
    		cp.setPixels(outPixels);
    		String imageTitle;
    		if (isOverlap && (this.channel.length>=2)) imageTitle=imp.getTitle()+String.format("OVLP_%02d-%02d",this.channel[0],this.channel[1]);
    		else  imageTitle=imp.getTitle()+"_PLNPRJ";
    			//isOverlap
    		ImagePlus impOut=new ImagePlus(imageTitle,cp);
    		impOut.setProperty("channel",   ""+channel);  // image channel
    		for (int i=0;i<this.channel.length;i++) impOut.setProperty("channel"+(i+1),   ""+this.channel[i]);  // image channel
    		(new JP46_Reader_camera(false)).encodeProperiesToInfo(impOut);
    		impOut.getProcessor().resetMinAndMax();
    		return impOut;
    	}
    }





	public class SensorData{
		public String path=null;
		public int    channel=   -1;
	    public int    subcamera= -1;
	    public int    subchannel=-1;
	    public int    sensor_port=-1;

// Following 4 parameters are not included in the vector, names, descriptions
	    public boolean cartesian = false; // cartesian coordinates mode (false - cylindrical)
    	public double right;   // distance to the right (radius*sin(azimuth))
    	public double forward; // distance forward (radius*cos(azimuth))
    	public double heading;   // absolute heading in degrees (used in cartesian mode)

    	public double azimuth; // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    	public double radius;  // mm, distance from the rotation axis
    	public double height;       // mm, up - from the origin point
    	public double phi;     // degrees, optical axis from azimuth/r vector, clockwise heading
    	public double theta;   // degrees, optical axis from the eyesis horizon, positive - up elevation
    	public double psi;     // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target roll
		public double focalLength=4.5;
		public double pixelSize=  2.2; //um
		public double lineTime =  3.638E-5; // 2.7778e-05 for Boson
		
		public boolean monochrome = false;
		public boolean lwir =       false;
		
		public double distortionRadius=  2.8512; // mm - half width of the sensor
		public double distortionA8=0.0; //r^8 (normalized to focal length or to sensor half width?)
		public double distortionA7=0.0; //r^7 (normalized to focal length or to sensor half width?)
		public double distortionA6=0.0; //r^6 (normalized to focal length or to sensor half width?)
		public double distortionA5=0.0; //r^5 (normalized to focal length or to sensor half width?)
		public double distortionA=0.0; // r^4 (normalized to focal length or to sensor half width?)
		public double distortionB=0.0; // r^3
		public double distortionC=0.0; // r^2
		public double px0=1296.0;          // center of the lens on the sensor, pixels
		public double py0=968.0;           // center of the lens on the sensor, pixels

	    public double [][] pixelCorrection=      null; // x,y, alpha, add flat, color, etc.
//	    public double []   sensorMask=           null;
	    public int    pixelCorrectionDecimation= 1;
	    public int    pixelCorrectionWidth=      -1; // 2592;
	    public int    pixelCorrectionHeight=     -1; // 1936;
	    public double entrancePupilForward=      0.0;


	    private double [] rByRDist=null;
	    private double    stepR=0.001;
	    private double    maxR=2.0; // calculate up to this*distortionRadius
	    public int    debugLevel=1;
	    public DirectMap directMap=null;
	    public EquirectangularMap equirectangularMap=null;

	    public int [][] defectsXY=null; // pixel defects coordinates list (starting with worst)
		public double [] defectsDiff=null; // pixel defects value (diff from average of neighbors), matching defectsXY

	    // TODO - add option to generate individual flat projections
		public InterSensor interSensor=null; // multiple sensors may have the same instance of the interSensor

		// new non-radial parameters
		final private boolean cummulativeCorrection=true; // r_xy, r_od for higher terms are relative to lower ones
		final private double [][] r_xy_dflt={{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0}}; // only 6, as for the first term delta x, delta y ==0
		final private double [][] r_od_dflt=   {{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0}}; // ortho

		public double [][] r_xy=null; // only 6, as for the first term delta x, delta y ==0
		public double [][] r_od=null; // ortho

		public double [][] r_xyod=null; //{x0,y0,ortho, diagonal}

		public int [] getSensorWH() {
			int [] wh = {this.pixelCorrectionWidth, this.pixelCorrectionHeight};
			return wh;
		}


	    public SensorData (String channelPath , boolean ok ){
	    	createEquirectangularMap(channelPath);
	    }

	    public void createEquirectangularMap(String channelPath){
	    	this.equirectangularMap= new EquirectangularMap(channelPath);

	    }

	    public void setDefaultNonRadial(){
			r_od=new double [r_od_dflt.length][2];
			for (int i=0;i<r_od.length;i++) r_od[i]=r_od_dflt[i].clone();
			r_xy=new double [r_xy_dflt.length][2];
			for (int i=0;i<r_xy.length;i++) r_xy[i]=r_xy_dflt[i].clone();
    	}


	    public double [][] getNonRadialCumulative(){
    		this.r_xyod=new double [this.r_od.length][4]; //{x0,y0,ortho, diagonal}
    		this.r_xyod[0][0]=0.0; // this.px0; //
    		this.r_xyod[0][1]=0.0; // this.py0;
    		this.r_xyod[0][2]=this.r_od[0][0];
    		this.r_xyod[0][3]=this.r_od[0][1];
    		if (cummulativeCorrection){
    			for (int i=1;i<this.r_xyod.length;i++){
    	    		this.r_xyod[i][0]=this.r_xyod[i-1][0]+this.r_xy[i-1][0];
    	    		this.r_xyod[i][1]=this.r_xyod[i-1][1]+this.r_xy[i-1][1];
    	    		this.r_xyod[i][2]=this.r_xyod[i-1][2]+this.r_od[i][0];
    	    		this.r_xyod[i][3]=this.r_xyod[i-1][3]+this.r_od[i][1];
    			}
    		} else {
    			for (int i=1;i<this.r_xyod.length;i++){
    	    		this.r_xyod[i][0]=this.r_xy[i-1][0]; // this.px0+this.r_xy[i-1][0];
    	    		this.r_xyod[i][1]=this.r_xy[i-1][1]; // this.py0+this.r_xy[i-1][1];
    	    		this.r_xyod[i][2]=this.r_od[i][0];
    	    		this.r_xyod[i][3]=this.r_od[i][1];
    			}
    		}
    		return this.r_xyod;
	    }


        public String [][] parameterDescriptions ={ // was 17, now 43
        		{"Azimuth",       "Subcamera azimuth, clockwise looking from top","degrees"},                                    //  0
        		{"Distance",      "Subcamera distance from the axis","mm"},                                                      //  1
        		{"Height",        "Subcamera height from the 'equator'","mm"},                                                   //  2
        		{"Heading",       "Optical axis heading (relative to azimuth)","degrees"},                                       //  3
        		{"Elevation",     "Optical axis elevation (up from equator)","degrees"},                                         //  4
        		{"Roll",          "Subcamera roll, positive CW looking to the target","degrees"},                                //  5
    			{"FocalLength",   "Lens focal length","mm","S","I"},                                                             //  6 15
    			{"PX0",           "Lens axis on the sensor (horizontal, from left edge)","pixels"},                              //  7 16
    			{"PY0",           "Lens axis on the sensor (vertical, from top edge)","pixels"},                                 //  8 17
    			{"DistortionA8",  "Distortion A8(r^8)","relative"},                                                              //  9 18
    			{"DistortionA7",  "Distortion A7(r^7)","relative"},                                                              // 10 19 18
    			{"DistortionA6",  "Distortion A6(r^6)","relative"},                                                              // 11 20 18
    			{"DistortionA5",  "Distortion A5(r^5)","relative"},                                                              // 12 21 18
    			{"DistortionA",   "Distortion A (r^4)","relative"},                                                              // 13 22 19
    			{"DistortionB",   "Distortion B (r^3)","relative"},                                                              // 14 23 20
    			{"DistortionC",   "Distortion C (r^2)","relative"},                                                              // 15 24 21
    			{"entrancePupilForward",   "Distance to actual entrance pupil from the lens center (positive - away from sensor)","mm"}, // 16 25

        		{"elong_C_o",   "Orthogonal elongation for r^2", "relative"},     //17
        		{"elong_C_d",   "Diagonal   elongation for r^2", "relative"},     //18

        		{"eccen_B_x",   "Distortion center shift X for r^3", "relative"}, //19
        		{"eccen_B_y",   "Distortion center shift Y for r^3", "relative"}, //20
        		{"elong_B_o",   "Orthogonal elongation for r^3", "relative"},     //21
        		{"elong_B_d",   "Diagonal   elongation for r^3", "relative"},     //22

        		{"eccen_A_x",   "Distortion center shift X for r^4", "relative"}, //23
        		{"eccen_A_y",   "Distortion center shift Y for r^4", "relative"}, //24
        		{"elong_A_o",   "Orthogonal elongation for r^4", "relative"},     //25
        		{"elong_A_d",   "Diagonal   elongation for r^4", "relative"},     //26

        		{"eccen_A5_x",  "Distortion center shift X for r^5", "relative"}, //27
        		{"eccen_A5_y",  "Distortion center shift Y for r^5", "relative"}, //28
        		{"elong_A5_o",  "Orthogonal elongation for r^5", "relative"},     //29
        		{"elong_A5_d",  "Diagonal   elongation for r^5", "relative"},     //30

        		{"eccen_A6_x",  "Distortion center shift X for r^6", "relative"}, //31
        		{"eccen_A6_y",  "Distortion center shift Y for r^6", "relative"}, //32
        		{"elong_A6_o",  "Orthogonal elongation for r^6", "relative"},     //33
        		{"elong_A6_d",  "Diagonal   elongation for r^6", "relative"},     //34

        		{"eccen_A7_x",  "Distortion center shift X for r^7", "relative"}, //35
        		{"eccen_A7_y",  "Distortion center shift Y for r^7", "relative"}, //36
        		{"elong_A7_o",  "Orthogonal elongation for r^7", "relative"},     //37
        		{"elong_A7_d",  "Diagonal   elongation for r^7", "relative"},     //38

        		{"eccen_A8_x",  "Distortion center shift X for r^8", "relative"}, //39
        		{"eccen_A8_y",  "Distortion center shift Y for r^8", "relative"}, //40
        		{"elong_A8_o",  "Orthogonal elongation for r^8", "relative"},     //41
        		{"elong_A8_d",  "Diagonal   elongation for r^8", "relative"}      //42
        };
	    public double [] getParametersVector(){ // was 17, now 43
	    	double [] vector={
	    			this.azimuth,      // 00
	    			this.radius,       // 01 mm, distance from the rotation axis
	    			this.height,       // 02 mm, up - from the origin point
	    			this.phi,          // 03 degrees, optical axis from azimuth/r vector, clockwise heading
	    			this.theta,        // 04 degrees, optical axis from the eyesis horizon, positive - up elevation
	    			this.psi,          // 05 degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target roll
	    			this.focalLength,  // 06
	    			this.px0,          // 07 center of the lens on the sensor, pixels
	    			this.py0,          // 08 center of the lens on the sensor, pixels
	    			this.distortionA8, // 09 r^8 (normalized to focal length or to sensor half width?)
	    			this.distortionA7, // 10 r^7 (normalized to focal length or to sensor half width?)
	    			this.distortionA6, // 11 r^6 (normalized to focal length or to sensor half width?)
	    			this.distortionA5, // 12 r^5 (normalized to focal length or to sensor half width?)
	    			this.distortionA,  // 13  r^4 (normalized to focal length or to sensor half width?)
	    			this.distortionB,  // 14  r^3
	    			this.distortionC,  // 15  r^2
	    			this.entrancePupilForward, // 16 Distance to actual entrance pupil from the lens center (positive - away from sensor)
	    			this.r_od[0][0],
	    			this.r_od[0][1],

	    			this.r_xy[0][0],
	    			this.r_xy[0][1],
	    			this.r_od[1][0],
	    			this.r_od[1][1],

	    			this.r_xy[1][0],
	    			this.r_xy[1][1],
	    			this.r_od[2][0],
	    			this.r_od[2][1],

	    			this.r_xy[2][0],
	    			this.r_xy[2][1],
	    			this.r_od[3][0],
	    			this.r_od[3][1],

	    			this.r_xy[3][0],
	    			this.r_xy[3][1],
	    			this.r_od[4][0],
	    			this.r_od[4][1],

	    			this.r_xy[4][0],
	    			this.r_xy[4][1],
	    			this.r_od[5][0],
	    			this.r_od[5][1],

	    			this.r_xy[5][0],
	    			this.r_xy[5][1],
	    			this.r_od[6][0],
	    			this.r_od[6][1]
	    	};
	    	return vector;
	    }
	    public double getParameter(int index){
	    	return getParametersVector()[index];
	    }
	    public String getParameterName(int index){
	    	return parameterDescriptions[index][0];
	    }
	    public String getParameterDescription(int index){
	    	return parameterDescriptions[index][1];
	    }
	    public String getParameterUnits(int index){
	    	return parameterDescriptions[index][2];
	    }
	    public String getPath(){return this.path;}

		/**
		 * Interpolate (bi-linear) X/Y corrections and flat-field data for the sensor
		 * @param rgbOnly - only flat field, no geometry
		 * @param px     - pixel X coordinate (non-decimated)
		 * @param py     - pixel Y coordinate (non-decimated)
		 * @return       - vector of {corrX, corrY, alpha, flatfield_red, flatfield_green, flatfield_blue}
		 */
		public double [] interpolateCorrectionVector (
				boolean rgbOnly,
				double px,
				double py){
			int start=rgbOnly?3:0;
			int sensorCorrWidth= (this.pixelCorrectionWidth-1)/this.pixelCorrectionDecimation+1;
			int sensorCorrHeight=this.pixelCorrection[0].length/sensorCorrWidth;
			int [] ix={(int) Math.floor(px/this.pixelCorrectionDecimation), (int) Math.floor(px/this.pixelCorrectionDecimation)+1};
			int [] iy={(int) Math.floor(py/this.pixelCorrectionDecimation),(int) Math.floor(py/this.pixelCorrectionDecimation)+1};
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
			double [] vector=new double [this.pixelCorrection.length-start];
			for (int n=0;n<vector.length;n++){
				// bilinear interpolation
				int n1=n+start;
				vector[n]=
					(1-corrDX)* (1-corrDY)* this.pixelCorrection[n1][index00]+
					corrDX * (1-corrDY)* this.pixelCorrection[n1][indexX0]+
					(1-corrDX)*    corrDY * this.pixelCorrection[n1][index0Y]+
					corrDX *    corrDY * this.pixelCorrection[n1][indexXY];
			}
			return vector;
		}

		public double [] getBayerFlatField(
//				int width,
//				int height,
				int [][] bayer){ //{{1,0},{2,1}} GR/BG
			int width = this.pixelCorrectionWidth; //  width,
			int height= this.pixelCorrectionHeight; // int height,

			double [] corrScale = new double [width*height];
			for (int y=0;y<height;y++) for (int x=0;x<width;x++){
				int iBayer=bayer[y&1][x&1];
				double [] v=interpolateCorrectionVector(
						true, // boolean rgbOnly,
						x,    //double px,
						y //double py)
						);

/*				double d=interpolateCorrectionVector(
						true, // boolean rgbOnly,
						x,    //double px,
						y //double py)
						)[iBayer];*/
				if (iBayer>=v.length){
					System.out.println("getBayerFlatField(): v.length="+v.length+" x="+x+" y="+y+" iBayer="+iBayer);
				}
				double d=v[iBayer];
				if (d!=0.0) d=1.0/d;
				corrScale[y*width+x]=d;
			}
			return corrScale;
		}
		public float [] getBayerFlatFieldFloat(
				int [][] bayer, //{{1,0},{2,1}} GR/BG
				double range) {
			return _getBayerFlatFieldFloat(
					this.pixelCorrectionWidth, //  width,
					this.pixelCorrectionHeight, // int height,
					bayer,
					range);
		}
		private float [] _getBayerFlatFieldFloat(
				int width,
				int height,
				int [][] bayer, //{{1,0},{2,1}} GR/BG
				double range) {

			float [] corrScale = new float [width*height];
			for (int y=0;y<height;y++) for (int x=0;x<width;x++){
				int iBayer=bayer[y&1][x&1];
				double [] v=interpolateCorrectionVector(
						true, // boolean rgbOnly,
						x,    //double px,
						y //double py)
						);
				if (iBayer>=v.length){
					System.out.println("getBayerFlatField(): v.length="+v.length+" x="+x+" y="+y+" iBayer="+iBayer);
				}
				double d=v[iBayer];
				if (d!=0.0) d=1.0/d;
				corrScale[y*width+x]= (float) d;
			}
			if (range > 0.0) {
				double [] mins = {Double.NaN,Double.NaN,Double.NaN,Double.NaN};
				int indx = 0;
				for (int y=0;y<height;y++) for (int x=0;x<width;x++){
					int iBayer= 2*(y&1) +(x&1);
					double d = corrScale[indx++];
					if (! ( d >=  mins[iBayer])) mins[iBayer] = d;
				}
				for (int i = 0; i < mins.length; i++) {
					mins[i] *= range;
				}
				indx = 0;
				for (int y=0;y<height;y++) for (int x=0;x<width;x++){
					int iBayer= 2*(y&1) +(x&1);
					if ( corrScale[indx] >=  mins[iBayer]) corrScale[indx] = (float) mins[iBayer];
					indx++;
				}
			}
			return corrScale;
		}
		public int [][] getDefectsXY(){
			return this.defectsXY;
		}
		public double [] getDefectsDiff(){
			return this.defectsDiff;
		}

	    public class  DirectMap{
	    	public int width;
	    	public int height;
	    	public double x0;
	    	public double y0;
	    	public double pixelStep;
	    	public double [][] map;
		    public DirectMap(
		    		int width,
		    		int height,
		    		double x0,
		    		double y0,
		    		double pixelStep,
		    		boolean flat,
		    		int numLayers){
		    	this.width=width;
		    	this.height=height;
		    	this.x0=x0;
		    	this.y0=y0;
		    	this.pixelStep=pixelStep;
		    	this.map=new double[numLayers][width*height];
		    }
	    }
	    public class  EquirectangularMap{
	    	public int    channel=          -1;
	    	public double longitudeLeft=    -180.0;
	    	public double longitudeRight=    180.0;
	    	public double latitudeTop=        90.0;
	    	public double latitudeBottom=    -90.0;
	    	public int    pixelsHorizontal=14000;
	    	public int    pixelsVertical;
	    	public double degreesPerPixel;

	    	public int width;
	    	public int height;
	    	public double x0;
	    	public double y0;
	    	public double pixelStep;
	    	public double minAlpha=0.1; // ignore pixels with mask below this value
	    	public int    numLayers;
	    	public float [][] map; // saving memory : first index - pixel #, empty -  null, otherwise {x,y,alpha, ...)
	    	public float [][] partialMap; //
	    	public Rectangle mapWOI=null;  // selection of the partialMap in the overall map. May

	    	public EquirectangularMap(String path){
		        setMapFromImageStack(path);
		    	this.degreesPerPixel=(this.longitudeRight-this.longitudeLeft)/(this.pixelsHorizontal-0); // with 360 - last pixel (this.longitudeRight) IS NOT included !

	    	}
//	    	public boolean flat=true;
		    public EquirectangularMap(
		    		int channel,
			    	double longitudeLeft,
			    	double longitudeRight,
			    	double latitudeTop,
			    	double latitudeBottom,
			    	int    pixelsHorizontal,
		    		int width,
		    		int height,
		    		double x0,
		    		double y0,
		    		double pixelStep,
		    		int numLayers){
		    	this.channel=         channel;
		    	this.longitudeLeft=   longitudeLeft;
		    	this.longitudeRight=  longitudeRight;
		    	this.latitudeTop=     latitudeTop;
		    	this.latitudeBottom=  latitudeBottom;
		    	this.pixelsHorizontal=pixelsHorizontal;
		    	while (this.longitudeRight<this.longitudeLeft)this.longitudeRight+=360.0;
		    	if (this.latitudeBottom>this.latitudeTop){
		    		double tmp=this.latitudeBottom;
		    		this.latitudeBottom=this.latitudeTop;
		    		this.latitudeTop=tmp;
		    	}
		    	this.degreesPerPixel=(this.longitudeRight-this.longitudeLeft)/(this.pixelsHorizontal-0); // with 360 - last pixel (this.longitudeRight) IS NOT included !
		    	this.pixelsVertical=(int) Math.round((this.latitudeTop-this.latitudeBottom)/this.degreesPerPixel) +1; // last pixel (this.latitudeBottom) IS included
		    	this.width=width;
		    	this.height=height;
		    	this.x0=x0;
		    	this.y0=y0;
		    	this.pixelStep=pixelStep;
		    	this.numLayers=numLayers;
		    	this.map=new float[this.pixelsVertical*this.pixelsHorizontal][];
		    }


	        public ImagePlus saveMapAsImageStack(String title, String path){
	        	ImagePlus imp=getMapAsImageStack(title);
	        	if (imp==null) return null;
		   				FileSaver fs=new FileSaver(imp);
		   				if (imp.getStackSize()>1)
		   					fs.saveAsTiffStack(path);
		   				else
		   					fs.saveAsTiff(path);
	        	return imp;
	        }

	        public ImagePlus getMapAsImageStack(String title){
	        	if (this.partialMap==null){
	        		String msg="Partial equirectangular map does not exist";
	        		IJ.showMessage("Error",msg);
	        		throw new IllegalArgumentException (msg);
	        	}
        		ImageStack stack=new ImageStack(this.mapWOI.width,this.mapWOI.height);
        		for (int n=0;n<this.partialMap.length;n++) {
        			if      (n==0) stack.addSlice("pixel-X", this.partialMap[n]);
        			else if (n==1) stack.addSlice("pixel-Y", this.partialMap[n]);
        			else if (n==2) stack.addSlice("mask",    this.partialMap[n]);
        			else            stack.addSlice("layer"+n,this.partialMap[n]);
        		}
        		ImagePlus imp = new ImagePlus(title, stack);
	// TODO: add more properties here (MAC+channel)? preserve other properties?
	        	imp.setProperty("channel",   ""+this.channel);  // image channel
	        	imp.setProperty("XPosition", ""+this.mapWOI.x); // real XPosition as tiff tag 0x11e is rational (5)
	        	imp.setProperty("YPosition", ""+this.mapWOI.y); // real XPosition as tiff tag 0x11f is rational (5)
	        	imp.setProperty("ImageFullWidth", "" +this.pixelsHorizontal);// real ImageFullWidth as tiff tag 0x8214 is rational (5)
	        	imp.setProperty("ImageFullLength", ""+this.pixelsVertical);   // real ImageFullLength as tiff tag 0x8215 is rational (5)

	        	imp.setProperty("longitudeLeft",  ""+this.longitudeLeft);
	        	imp.setProperty("longitudeRight", ""+this.longitudeRight);
	        	imp.setProperty("latitudeTop",    ""+this.latitudeTop);
	        	imp.setProperty("latitudeBottom", ""+this.latitudeBottom);

	        	(new JP46_Reader_camera(false)).encodeProperiesToInfo(imp);
	        	imp.getProcessor().resetMinAndMax();
	        	return imp;
	        }

	        public void setMapFromImageStack(String path){
	        	Opener opener=new Opener();
	        	ImagePlus imp=opener.openImage("", path);
	        	if (imp==null) {
	        		String msg="Failed to read map file "+path;
	        		IJ.showMessage("Error",msg);
	        		throw new IllegalArgumentException (msg);
	        	}
	        	(new JP46_Reader_camera(false)).decodeProperiesFromInfo(imp);
	        	this.mapWOI=new Rectangle();
	        	this.mapWOI.width=imp.getWidth();
	        	this.mapWOI.height=imp.getHeight();
	        	if (imp.getProperty("channel")!=null)
	        		this.channel=Integer.parseInt ((String) imp.getProperty("channel"));
	        	if (imp.getProperty("XPosition")!=null)
	        		this.mapWOI.x=Integer.parseInt ((String) imp.getProperty("XPosition"));
	        	if (imp.getProperty("YPosition")!=null)
	        		this.mapWOI.y=Integer.parseInt ((String) imp.getProperty("YPosition"));
	        	if (imp.getProperty("ImageFullWidth")!=null)
	        		this.pixelsHorizontal=Integer.parseInt ((String) imp.getProperty("ImageFullWidth"));
	        	if (imp.getProperty("ImageFullLength")!=null)
	        		this.pixelsVertical=Integer.parseInt ((String) imp.getProperty("ImageFullLength"));

	        	if (imp.getProperty("longitudeLeft")!=null)
	        		this.longitudeLeft=Double.parseDouble ((String) imp.getProperty("longitudeLeft"));
	        	if (imp.getProperty("longitudeRight")!=null)
	        		this.longitudeRight=Double.parseDouble ((String) imp.getProperty("longitudeRight"));
	        	if (imp.getProperty("latitudeTop")!=null)
	        		this.latitudeTop=Double.parseDouble ((String) imp.getProperty("latitudeTop"));
	        	if (imp.getProperty("latitudeBottom")!=null)
	        		this.latitudeBottom=Double.parseDouble ((String) imp.getProperty("latitudeBottom"));

        		ImageStack stack = imp.getStack();
            	if (stack==null) {
            		String msg="Expected a image stack with masks";
            		IJ.showMessage("Error",msg);
            		throw new IllegalArgumentException (msg);
            	}
	    		int numChannels=imp.getStackSize();
	    		this.partialMap =new float[numChannels][];
            	for (int i=0;i<numChannels;i++) this.partialMap[i]= (float[]) stack.getPixels(i+1);
//            	System.out.println("setMapFromImageStack(): loaded "+path);
	        }



		    public int createPartialMap(int width){ // does not delete the original map
		    	int [] histogram=new int [this.pixelsHorizontal];
		    	for (int i=0;i<this.pixelsHorizontal;i++) histogram[i]=0;
		    	int minILat=this.pixelsVertical;
		    	int maxILat=0;
		    	int index=0;
		    	int numLayers=-1;
		    	for (int iLat=0;iLat<pixelsVertical;iLat++) for (int iLong=0;iLong<pixelsHorizontal;iLong++) if (this.map[index++]!=null){
		    		if (numLayers<0) {
		    			numLayers=this.map[index-1].length;
		    		}
		    		histogram[iLong]++;
		    		if (minILat>iLat) minILat=iLat;
		    		if (maxILat<iLat) maxILat=iLat;
		    	}
		    	int maxHist=0;
		    	int minILong=-1;
		    	for (int iLong=0;iLong<this.pixelsHorizontal;iLong++) if (maxHist<histogram[iLong]) {
		    		maxHist=histogram[iLong];
		    		minILong=iLong;
		    	}
// grow horizontally (rollover allowed) to the specified width
//		    	for (int maxILong=minILong; maxILong<(minILong+width-1);){
		    	for (int i=1;i<=width;i++) {
		    		int maxILong=(minILong+i)% this.pixelsHorizontal;
		    		int minILong1=(minILong-1+this.pixelsHorizontal)%this.pixelsHorizontal;
		    		if ((histogram[maxILong]==0) && (histogram[minILong1]==0)){
		    			width=i;
		    			break;
		    		}
		    		if (histogram[maxILong]<histogram[minILong1]) minILong=minILong1; //>=0
		    	}
		    	this.mapWOI=new Rectangle(minILong,minILat,width,maxILat-minILat+1);
		    	this.partialMap=new float [numLayers][this.mapWOI.width*this.mapWOI.height];
		    	int indexDst=0;
		    	int numPix=0;
		    	for (int iLat=0;iLat<this.mapWOI.height;iLat++) {
		    		int iLatSrc=iLat+this.mapWOI.y;
		    		for (int iLong=0;iLong<this.mapWOI.width;iLong++){
		    			int iLongSrc=(iLong+this.mapWOI.x)%this.pixelsHorizontal;
		    			int indexSrc=iLatSrc*this.pixelsHorizontal+iLongSrc;
		    			if (this.map[indexSrc]!=null){
		    				for (int n=0;n<numLayers;n++) this.partialMap[n][indexDst]=this.map[indexSrc][n];
		    				numPix++;
		    			}	else                          for (int n=0;n<numLayers;n++) this.partialMap[n][indexDst]=0;
		    			indexDst++;
		    		}
		    	}
		    	return numPix;

		    }

	    }
	    public SensorData (String path, int debugLevel){
//			int indexPeriod=path.indexOf('.',path.lastIndexOf(Prefs.getFileSeparator())); // channel will be overwritten if available in properties
//	    	this.channel=Integer.parseInt(path.substring(indexPeriod-2,indexPeriod));
	    	this.debugLevel=debugLevel;
			Opener opener=new Opener();
			ImagePlus imp=opener.openImage("", path);
	    	if (imp==null) {
	    		String msg="Failed to read sensor calibration data file "+path;
	    		if (this.debugLevel>2) IJ.showMessage("Error",msg);
	    		if (this.debugLevel>0)System.out.println(msg);
	    		return;
	    	}
			if (this.debugLevel>1) System.out.println("Read "+path+" as a sensor calibration data");
	    	(new JP46_Reader_camera(false)).decodeProperiesFromInfo(imp);
	    	setSensorDataFromImageStack(imp);
	    }
	    public SensorData (){} // just to get parameter names

	    public void setChannel(int chn){this.channel = chn;}
	    public int getChannel(){return this.channel;}
	    public int getSubChannel(){return this.subchannel;}
	    public int getSubCamera(){return this.subcamera;}
	    public int getSensorPort(){return this.sensor_port;}


	    public void setSensorDataFromImageStack(ImagePlus imp){
//	    	int corrX=0,corrY=1,corrMask=2;
	    	boolean need_lineTime_fix =   false; // will use sensor default if not provided, 3.638E-5 for RGB and 2.7778e-05 for Boson
	    	boolean need_monochrome_fix = false; // will use sensor default if not provided, 3.638E-5 for RGB and 2.7778e-05 for Boson
	    	boolean need_lwir_fix =       false; // will use sensor default if not provided, 3.638E-5 for RGB and 2.7778e-05 for Boson
	    	if (imp == null){
	    		String msg="Sensor Calibration image is null";
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
	        		"distortionA5",
	        		"distortionA",
	        		"distortionB",
	        		"distortionC",
	        		"px0",
	        		"py0",
	    	        "azimuth",
	    	        "radius",
	    	        "height",
	    	        "heading",
	    	        "elevation",
	    	        "roll",
	    	        "channel",
	        		"lineTime",   // optional
	        		"monochrome", // optional
	        		"lwir"        // optional
	        		};
	        for (int i=0; i<requiredProperties.length;i++) if (imp.getProperty(requiredProperties[i])==null){
	    		String msg="Required property "+requiredProperties[i]+" is not defined in "+imp.getTitle();
	    		// is it 
	    		if (requiredProperties[i].equals("lineTime")) {
	    			System.out.println(requiredProperties[i]+" is not provided, will use default");
	    			need_lineTime_fix = true;
	    			continue;
	    		}
	    		if (requiredProperties[i].equals("monochrome")) {
	    			System.out.println(requiredProperties[i]+" is not provided, will use default");
	    			need_monochrome_fix = true;
	    			continue;
	    		}
	    		if (requiredProperties[i].equals("lwir")) {
	    			System.out.println(requiredProperties[i]+" is not provided, will use default");
	    			need_lwir_fix = true;
	    			continue;
	    		}
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



	        this.pixelCorrectionWidth=      Integer.parseInt  ((String) imp.getProperty("pixelCorrectionWidth"));
	        this.pixelCorrectionHeight=     Integer.parseInt  ((String) imp.getProperty("pixelCorrectionHeight"));
	        this.pixelCorrectionDecimation= Integer.parseInt  ((String) imp.getProperty("pixelCorrectionDecimation"));

	        this.distortionRadius= Double.parseDouble((String) imp.getProperty("distortionRadius"));
	        this.focalLength=      Double.parseDouble((String) imp.getProperty("focalLength"));
	        this.pixelSize=        Double.parseDouble((String) imp.getProperty("pixelSize"));
	        if (need_lineTime_fix) {
	        	if (this.pixelSize < 5.0) {
	        		this.lineTime = 3.638E-5;
	        	} else if (this.pixelCorrectionWidth == 640){ // Boson
	        		this.lineTime = 2.7778e-05; // 12um pixel, Boson
	        	} else {
	        		this.lineTime = 7.8e-05; // 12um pixel, Lepton (may7 be wrong)
	        	}
	        } else {
		        this.lineTime=        Double.parseDouble((String) imp.getProperty("lineTime"));
	        }

	        if (need_monochrome_fix) {
	        	if (this.pixelSize < 5.0) {
	        		this.monochrome = false;
	        	} else if (this.pixelCorrectionWidth == 640){ // Boson
	        		this.monochrome = true; // 12um pixel, Boson
	        	} else {
	        		this.monochrome = true; // 12um pixel, Lepton (may7 be wrong)
	        	}
	        } else {
		        this.monochrome=        Boolean.parseBoolean((String) imp.getProperty("monochrome"));
	        }
	        
	        
	        if (need_lwir_fix) {
	        	if (this.pixelSize < 5.0) {
	        		this.lwir = false;
	        	} else if (this.pixelCorrectionWidth == 640){ // Boson
	        		this.lwir = true; // 12um pixel, Boson
	        	} else {
	        		this.lwir = true; // 12um pixel, Lepton (may7 be wrong)
	        	}
	        } else {
		        this.lwir=        Boolean.parseBoolean((String) imp.getProperty("lwir"));
	        }
	        
	        if (imp.getProperty("distortionA8")!=null)	this.distortionA8=     Double.parseDouble((String) imp.getProperty("distortionA8"));
	        else this.distortionA8= 0.0;
	        if (imp.getProperty("distortionA7")!=null)	this.distortionA7=     Double.parseDouble((String) imp.getProperty("distortionA7"));
	        else this.distortionA7= 0.0;
	        if (imp.getProperty("distortionA6")!=null)	this.distortionA6=     Double.parseDouble((String) imp.getProperty("distortionA6"));
	        else this.distortionA6= 0.0;
	        this.distortionA5=     Double.parseDouble((String) imp.getProperty("distortionA5"));
	        this.distortionA=      Double.parseDouble((String) imp.getProperty("distortionA"));
	        this.distortionB=      Double.parseDouble((String) imp.getProperty("distortionB"));
	        this.distortionC=      Double.parseDouble((String) imp.getProperty("distortionC"));
	        this.px0=              Double.parseDouble((String) imp.getProperty("px0"));
	        this.py0=              Double.parseDouble((String) imp.getProperty("py0"));
	        this.azimuth=          Double.parseDouble((String) imp.getProperty("azimuth"));
	        this.radius=           Double.parseDouble((String) imp.getProperty("radius"));
	        this.height=           Double.parseDouble((String) imp.getProperty("height"));
	        if (imp.getProperty("entrancePupilForward")!=null) this.entrancePupilForward= Double.parseDouble((String) imp.getProperty("entrancePupilForward"));
	        this.phi=              Double.parseDouble((String) imp.getProperty("heading"));
	        this.theta=            Double.parseDouble((String) imp.getProperty("elevation"));
	        this.psi=              Double.parseDouble((String) imp.getProperty("roll"));
	        this.channel=          Integer.parseInt  ((String) imp.getProperty("channel"));
// older files do not have these properties
	        if (imp.getProperty("subcamera")!=null)  this.subcamera= Integer.parseInt((String) imp.getProperty("subcamera"));
	        if (imp.getProperty("subchannel")!=null) this.subchannel=Integer.parseInt((String) imp.getProperty("subchannel"));
	        if (imp.getProperty("sensor_port")!=null) this.sensor_port=Integer.parseInt((String) imp.getProperty("sensor_port"));



	        // now read the calibration data and mask
	        	this.pixelCorrection=null;
	        this.pixelCorrection=new double [pixels.length] [pixels[0].length];
	        for (int n=0;n<pixels.length;n++) for (int i= 0;i<this.pixelCorrection[0].length;i++){
	        	this.pixelCorrection[n][i]=pixels[n][i];
	        }

	        if (imp.getProperty("defects")!=null) {
        		String sDefects=(String) imp.getProperty("defects");
        		String [] asDefects=sDefects.trim().split(" ");
        		this.defectsXY=new int [asDefects.length][2];
        		this.defectsDiff=new double [asDefects.length];
        		for (int i=0;i<asDefects.length;i++) {
        			String [] stDefect=asDefects[i].split(":");
        			this.defectsXY[i][0]=Integer.parseInt(stDefect[0]);
        			this.defectsXY[i][1]=Integer.parseInt(stDefect[1]);
        			this.defectsDiff[i]=Double.parseDouble(stDefect[2]);
        		}
        	} else {
        		this.defectsXY=null;
        		this.defectsDiff=null;
        	}

        	setDefaultNonRadial();
			for (int i=0;i<this.r_xy.length;i++) {
				if (imp.getProperty("r_xy_"+i+"_x")  !=null) this.r_xy[i][0]= Double.parseDouble((String) imp.getProperty("r_xy_"+i+"_x"));
				if (imp.getProperty("r_xy_"+i+"_y")  !=null) this.r_xy[i][1]= Double.parseDouble((String) imp.getProperty("r_xy_"+i+"_y"));
			}
			for (int i=0;i<this.r_od.length;i++) {
				if (imp.getProperty("r_od_"+i+"_o")  !=null) this.r_od[i][0]= Double.parseDouble((String) imp.getProperty("r_od_"+i+"_o"));
				if (imp.getProperty("r_od_"+i+"_d")  !=null) this.r_od[i][1]= Double.parseDouble((String) imp.getProperty("r_od_"+i+"_d"));
			}

        	if (imp.getProperty("forward")  !=null) this.forward=  Double.parseDouble((String) imp.getProperty("forward"));
        	if (imp.getProperty("right")    !=null) this.right=    Double.parseDouble((String) imp.getProperty("right"));
        	if (imp.getProperty("aheading") !=null) this.heading=  Double.parseDouble((String) imp.getProperty("aheading"));

        	if (imp.getProperty("cartesian") !=null) {
        		this.cartesian= Boolean.parseBoolean((String) imp.getProperty("cartesian"));
        		updateCartesian(); // recalculate other parameters (they may or may nort be provided
        	} else {
        		this.cartesian = false;
        	}
	    }

    	public void updateCartesian(){ // set alternative parameters
    		if (cartesian) {
    			this.azimuth = Math.atan2(this.right, this.forward)*180.0/Math.PI;
    			this.radius = Math.sqrt(this.forward*this.forward + this.right*this.right);
    			this.phi = this.heading - this.azimuth;
    		} else {
    			this.forward = this.radius * Math.cos(Math.PI*this.azimuth/180.0);
    			this.right =   this.radius * Math.sin(Math.PI*this.azimuth/180.0);
    			this.heading = this.phi +  this.azimuth;
    		}
    	}



	    /**
	     * Calculate rotation matrix that converts sensor coordinates to world coordinates
	     * @return 3x3 rotation matrix
	     */

	    public Matrix rotateSensorCoordToPanoCoord(){
        	double sAZP=Math.sin((this.azimuth+this.phi)*Math.PI/180);
        	double cAZP=Math.cos((this.azimuth+this.phi)*Math.PI/180);
        	double sTH=Math.sin(this.theta*Math.PI/180);
        	double cTH=Math.cos(this.theta*Math.PI/180);
        	double sPS=Math.sin(this.psi*Math.PI/180);
        	double cPS=Math.cos(this.psi*Math.PI/180);
        	/*
        	Converting from the sub-camera coordinates to the target coordinates
        	1) rotate by -psi around CZ: Vc1= R1*Vc
        	| Xc1 |   | cos(psi)  sin(psi)    0         |   |Xc|
        	| Yc1 | = |-sin(psi)  cos(psi)    0         | * |Yc|
        	| Zc1 |   |    0         0        1         |   |Zc|
        	*/
        	    	double [][] aR1={{cPS,sPS,0.0},{-sPS,cPS,0.0},{0.0,0.0,1.0}};
        	    	Matrix R1=new Matrix(aR1);
        	/*
        	2) rotate by - theta around C1X:Vc2= R2*Vc1
        	| Xc2 |   |    1         0         0        |   |Xc1|
        	| Yc2 | = |    0    cos(theta)   sin(theta) | * |Yc1|
        	| Zc2 |   |    0   -sin(theta)   cos(theta) |   |Zc1|
        	*/
        	    	double [][] aR2={{1.0,0.0,0.0},{0.0,cTH,sTH},{0.0,-sTH,cTH}};
        	    	Matrix R2=new Matrix(aR2);
        	/*
        	3) rotate by -(azimuth+phi) around C2Y:Vc3= R3*Vc2
        	| Xc3 |   | cos(azimuth+phi)    0   sin(azimuth+phi)   |   |Xc2|
        	| Yc3 | = |     0               1         0            | * |Yc2|
        	| Zc3 |   | -sin(azimuth+phi)   0   cos(azimuth+phi)   |   |Zc2|
        	*/
        	    	double [][] aR3={{cAZP,0.0,sAZP},{0.0,1.0,0.0},{-sAZP,0.0,cAZP}};
        	    	Matrix R3=new Matrix(aR3);
        	    	// MA=R3*R2*R1;
        	    	// MB=T3+(R8*R7*R6*R5*(T2+R4*T1));
        	    	Matrix MA=R3.times(R2.times(R1));
        	    	return MA;

	    }
	    /**
	     * Build a pixel map from the (possibly decimated and shifted) pixel array to X,Y (in mm) from the lens axis. Will process alpha and flat field layers, if present
	     * @param width  width of the result map
	     * @param height height of the result map
	     * @param x0 top left corner of the map X. With different resolutions the common 0.0 is half pixels size (up-left) from the center of the pixel [0,0] in current resolution
	     * No, 0,0 is the center of the pixel 0,0
	     *  corrections use floor
	     * @param y0 top left corner of the map Y
	     * @param pixelStep step in sensor pixels (i.e step=0.5 will have twice number of pixels
	     */
	    public void combineDistortionsSensorToEquirectangular(
	    		int width,
	    		int height,
	    		double x0,
	    		double y0,
	    		double pixelStep,
	    		boolean flat){
	    	System.out.println("******* BUG: combineDistortionsSensorToEquirectangular() - not yet modified for non-radial !!!! ********");
	    	if (this.pixelCorrection==null){
	    		String msg="No pixel correction data";
	    		IJ.showMessage("Error",msg);
	    		throw new IllegalArgumentException (msg);
	    	}
	    	this.directMap=new DirectMap(
	    			width,
	    			height,
	    			x0,
	    			y0,
	    			pixelStep,
	    			flat,
	    			this.pixelCorrection.length);
	    	int sensorCorrWidth=  (this.pixelCorrectionWidth-1)/this.pixelCorrectionDecimation+1;
	    	int sensorCorrHeight= (this.pixelCorrectionHeight-1)/this.pixelCorrectionDecimation+1;
	    	Matrix MA=rotateSensorCoordToPanoCoord();
	    	for (int iy=0;iy<this.directMap.height;iy++){
	    		for (int ix=0;ix<this.directMap.width;ix++) {

	    			boolean debugThis= false; //(iy>=10) && (iy<=20) &&  (ix>=10) && (ix<=20);
	    			int outIndex=ix+this.directMap.width*iy;
	    			double x=x0+pixelStep*ix;
	    			double y=y0+pixelStep*iy;
	    			int corrX=((int) Math.floor(x/this.pixelCorrectionDecimation));
	    			int corrY=((int) Math.floor(y/this.pixelCorrectionDecimation));
	    			double corrDX=x/this.pixelCorrectionDecimation-corrX;
	    			double corrDY=y/this.pixelCorrectionDecimation-corrY;
	    			if (corrX<0)                      {corrX=0;corrDX=0;}
	    			else if (corrX>=sensorCorrWidth)  {corrX=sensorCorrWidth-1;corrDX=0;}
	    			if (corrY<0)                      {corrY=0;corrDY=0;}
	    			else if (corrY>=sensorCorrHeight) {corrY=sensorCorrHeight-1;corrDY=0;}
	    			int index00=corrX + corrY*sensorCorrWidth;
	    			int indexX0=index00+((corrX==(sensorCorrWidth-1))?0:1);
	    			int index0Y=index00+((corrY==(sensorCorrHeight-1))?0:sensorCorrWidth);
	    			int indexXY=index0Y+((corrX==(sensorCorrWidth-1))?0:1);
	    			if (debugThis){
	    				System.out.println("combineDistortionsSensorToFlat(): x="+x+" y="+y+" corrX="+corrX+ " corrY="+corrY+" corrDX="+corrDX+ " corrDY="+corrDY);
	    				System.out.println("combineDistortionsSensorToFlat(): index00="+index00+" indexX0="+indexX0+" index0Y="+index0Y+" indexXY="+indexXY);
	    			}
	    			// Linear interpolate all this.pixelCorrection[i][indexXY]
	    			double [] corr=new double [this.pixelCorrection.length];
	    			for (int n=0;n<this.pixelCorrection.length;n++){
	    				corr[n]=
	    					(1-corrDX)* (1-corrDY)* this.pixelCorrection[n][index00]+
	    					corrDX * (1-corrDY)* this.pixelCorrection[n][indexX0]+
	    					(1-corrDX)*    corrDY * this.pixelCorrection[n][index0Y]+
	    					corrDX *    corrDY * this.pixelCorrection[n][indexXY];
	    				if (debugThis){
	    					System.out.println("combineDistortionsSensorToFlat(): this.pixelCorrection["+n+"]["+index00+"]="+this.pixelCorrection[n][index00]);
	    					System.out.println("combineDistortionsSensorToFlat(): this.pixelCorrection["+n+"]["+indexX0+"]="+this.pixelCorrection[n][indexX0]);
	    					System.out.println("combineDistortionsSensorToFlat(): this.pixelCorrection["+n+"]["+index0Y+"]="+this.pixelCorrection[n][index0Y]);
	    					System.out.println("combineDistortionsSensorToFlat(): this.pixelCorrection["+n+"]["+indexXY+"]="+this.pixelCorrection[n][indexXY]);
	    				}
	    			}
	    			x-=corr[0]; // same sign as in Distortions.initFittingSeries()
	    			y-=corr[1];
	    			x-=this.px0; //
	    			y-=this.py0; // from the lens center
	    			x*=0.001*this.pixelSize;
	    			y*=0.001*this.pixelSize;
	    			double rD=Math.sqrt(x*x+y*y);
	    			if (debugThis){
	    				System.out.println("combineDistortionsSensorToFlat() 1: x="+x+" y="+y+" rD="+rD);
	    			}

	    			if (this.rByRDist==null){
	    				calcReverseDistortionTable();
	    			}
	    			double rND2R=getRByRDist(rD/this.distortionRadius,debugThis);
	    			x*= rND2R; // positive - right
	    			y*=-rND2R; // positive - up
	    			if (debugThis){
	    				System.out.println("combineDistortionsSensorToFlat() 2: x="+x+" y="+y+" rND2R="+rND2R);
	    			}
	    			this.directMap.map[0][outIndex]=x;
	    			this.directMap.map[1][outIndex]=y;
	    			for (int n=2; n<this.pixelCorrection.length;n++){
	    				this.directMap.map[n][outIndex]=corr[n];
	    			}
	    			if (!flat){ // replace x/y with azimuth/elevation
	    				double [][] aB={{x},{y},{this.focalLength}};
	    				Matrix MB=new Matrix(aB);
	    				double [] xyz=MA.times(MB).getRowPackedCopy();
	    				double dist=Math.sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
	    				for (int n=0; n<3;n++) xyz[n]/=dist;
	    				// now we have unity-length vector
	    				double r=Math.sqrt(xyz[0]*xyz[0]+xyz[2]*xyz[2]); // horizontal
	    				if (r>0){
	    					this.directMap.map[0][outIndex]=Math.atan2(xyz[0],xyz[2])*180.0/Math.PI;
	    					this.directMap.map[1][outIndex]=Math.atan2(xyz[1],r)*180.0/Math.PI;
	    				} else {
	    					this.directMap.map[0][outIndex]=0;
	    					this.directMap.map[1][outIndex]=(xyz[1]>0)?90.0:-90.0;
	    				}
	    			}
	    		}
	    		IJ.showProgress(iy, this.directMap.height-1);
	    	}
	    }

	    public double ipow(double a, int b){
	    	switch (b) {
	    	case 1: return a;
	    	case 2: return a*a;
	    	case 3: return a*a*a;
	    	case 0: return 1.0;
	    	default:
	    		if (b<0) {
	    			return 1.0/ipow(a,-b);
	    		} else {
	    			int b1=b>>1;
	    		return ipow(a,b1)*ipow(a,b-b1);
	    		}
	    	}
	    }


	    public void combineDistortionsEquirectangularToSensor(
	    		int    channel,
		    	double longitudeLeft0,
		    	double longitudeRight0,
		    	double latitudeTop0,
		    	double latitudeBottom0,
		    	int    pixelsHorizontal0,
	    		int width0,
	    		int height0,
	    		double x00,
	    		double y00,
	    		double pixelStep0,
	    		int maxThreads){
	    	final double maxKR=2.0; // maximal ratio to distortion radius to consider
	    	if (this.pixelCorrection==null){
	    		String msg="No pixel correction data";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
	    	}
	    	System.out.println("******* INFO: combineDistortionsEquirectangularToSensor() - already modified for non-radial !!!! ********");
	    	this.equirectangularMap=new EquirectangularMap(
	    			channel,
			    	longitudeLeft0,
			    	longitudeRight0,
			    	latitudeTop0,
			    	latitudeBottom0,
			    	pixelsHorizontal0,
	    			width0,
		    		height0,
		    		x00,
		    		y00,
		    		pixelStep0,
		    		this.pixelCorrection.length);
	    	final boolean debugThis=false;
			final int sensorCorrWidth=  (this.pixelCorrectionWidth-1)/this.pixelCorrectionDecimation+1;
			final int sensorCorrHeight= (this.pixelCorrectionHeight-1)/this.pixelCorrectionDecimation+1;
	    	final Matrix MA=rotateSensorCoordToPanoCoord().transpose(); // panorama coordinates to sensor coordinates (reversed rotation)
	   		final Thread[] threads = newThreadArray(maxThreads);
	   		final AtomicInteger ipLatAtomic     = new AtomicInteger(0);
	   		final AtomicInteger ipLatDoneAtomic = new AtomicInteger(1);

	    	final double [][] pixelCorrection=this.pixelCorrection;
	    	final double pixelCorrectionDecimation=this.pixelCorrectionDecimation;
	    	final double focalLength=     this.focalLength;
	    	final double px0=             this.px0;
	    	final double py0=             this.py0;
	    	final double distortionRadius=this.distortionRadius;
	    	final double pixelSize=       this.pixelSize;
	    	final double distortionA8=    this.distortionA8;
	    	final double distortionA7=    this.distortionA7;
	    	final double distortionA6=    this.distortionA6;
	    	final double distortionA5=    this.distortionA5;
	    	final double distortionA =    this.distortionA;
	    	final double distortionB =    this.distortionB;
	    	final double distortionC =    this.distortionC;
	    	final double [][] r_xyod=getNonRadialCumulative();

//	    	final double d=1.0-this.distortionA8-this.distortionA7-this.distortionA6-this.distortionA5-this.distortionA-this.distortionB-this.distortionC;
//	    	final boolean use8=(this.distortionA8!=0.0) || (this.distortionA7!=0.0) || (this.distortionA6!=0.0);
	    	final EquirectangularMap equirectangularMap=this.equirectangularMap;

//	    	for (int iLat=0;iLat<equirectangularMap.pixelsVertical;iLat++){
	    	for (int ithread = 0; ithread < threads.length; ithread++) {
	    		threads[ithread] = new Thread() {
	    			@Override
					public void run() {
	    				for (int iLat=ipLatAtomic.getAndIncrement(); iLat<equirectangularMap.pixelsVertical;iLat=ipLatAtomic.getAndIncrement()){
	    		        	double [] a={distortionC,distortionB,distortionA,distortionA5,distortionA6,distortionA7,distortionA8};
	    					double latitude=equirectangularMap.latitudeTop-equirectangularMap.degreesPerPixel*iLat;
	    					for (int iLong=0;iLong<equirectangularMap.pixelsHorizontal;iLong++){
	    						int outIndex=iLat*equirectangularMap.pixelsHorizontal+iLong;
	    						double [] pXY=null;
	    						double longitude=equirectangularMap.longitudeLeft+equirectangularMap.degreesPerPixel*iLong;
	    						double [][] xyz={ // unity vector in the direction
	    								{Math.cos(latitude*Math.PI/180.0)*Math.sin(longitude*Math.PI/180.0)},
	    								{Math.sin(latitude*Math.PI/180.0)},
	    								{Math.cos(latitude*Math.PI/180.0)*Math.cos(longitude*Math.PI/180.0)}};
	    						// get sensor coordinates
	    						double [] v=MA.times(new Matrix (xyz)).getRowPackedCopy();
	    						if (v[2]>0){ // potentially visible
	    							v[0]*=focalLength/v[2];
	    							v[1]*=focalLength/v[2];
	    							v[2]= focalLength;
	    							double rND=Math.sqrt(v[0]*v[0]+v[1]*v[1]);
	    							double r=rND/distortionRadius;
	    							if (r<maxKR){ // skip if they are too far from the center
	    					        	double deltaX=0.0,deltaY=0.0; // difference between distorted and non-distorted as a fraction of this.distortionRadius
	    					        	double xmmc=v[0]/distortionRadius; // non-distorted x relative to distortion radius
	    					        	double ymmc=v[1]/distortionRadius; // non-distorted y relative to distortion radius (math direction, not image)
	    					        	for (int i=0;i<r_xyod.length;i++){
	    					        		double x=xmmc-r_xyod[i][0]; // relative X-shift of this term center
	    					        		double y=ymmc-r_xyod[i][1]; // relative Y-shift of this term center
	    					        		double x2=x*x;   // relative squared X-shift from this term center
	    					        		double y2=y*y;   // relative squared Y-shift from this term center
	    					        		double r2=x2+y2; // relative squared distance from this term center
	    					        		// effective distance from this term center corrected for elongation
	    					        		double rr=Math.sqrt(r2+ r_xyod[i][2]*(y2-x2)+ 2.0*r_xyod[i][3]*x*y );
	    					                if (rr<0.00000001*distortionRadius) rr=0.00000001*distortionRadius; // avoid 1/0.0 -// Added 3 zeros (old comment)
	    					        		double rr_pow_i=ipow(rr,i);
	    					        		double ki=a[i]*(rr_pow_i*rr-1.0); // scaling of distorted distance from this term's center (does not include pinhole center shift itself)
	    					        		deltaX+=x*ki; // relative distorted distance from the center
	    					        		deltaY+=y*ki;
//	    					        		if ((debugLevel>2) && ((r_xyod[i][0]!=0.0) || (r_xyod[i][0]!=0.0))){
//	    					        			System.out.println ("i="+i+" r_xyod[i][0]="+r_xyod[i][0]+" r_xyod[i][1]="+r_xyod[i][1]+" a[0]="+a[0]+" a[1]="+a[1]+" a[2]="+a[2]+" a[3]="+a[3]);
//	    					        		}
	    					        	}
//	    					        	double [] xyDist={v[0]+distortionRadius*deltaX, v[1]+distortionRadius*deltaY};
	    								pXY=new double[2];
	    					        	// convert to sensor pixels coordinates
	    								pXY[0]= 1000.0/pixelSize*(v[0]+distortionRadius*deltaX) + px0;
	    								pXY[1]=-1000.0/pixelSize*(v[1]+distortionRadius*deltaY) + py0;
	    								/*
	    								double k;
	    								if (use8) k=((((((distortionA8*r+distortionA7)*r+distortionA6)*r+distortionA5)*r + distortionA)*r+distortionB)*r+distortionC)*r+d;
	    								else	  k=(((distortionA5*r + distortionA)*r+distortionB)*r+distortionC)*r+d;
	    								// calculate in sensor pixel
	    								pXY[0]= v[0]*k*1000.0/pixelSize+px0; // in pixels, right
	    								pXY[1]=-v[1]*k*1000.0/pixelSize+py0; // in pixels, down
	    								*/
	    								// un-apply residual correction, for now - consider it small/smooth to skip reverse mapping, make it a two step
	    								int corrX=((int) Math.floor(pXY[0]/pixelCorrectionDecimation));
	    								int corrY=((int) Math.floor(pXY[1]/pixelCorrectionDecimation));
	    								if (corrX<0)                      corrX=0;
	    								else if (corrX>=sensorCorrWidth)  corrX=sensorCorrWidth-1;
	    								if (corrY<0)                      corrY=0;
	    								else if (corrY>=sensorCorrHeight) corrY=sensorCorrHeight-1;
	    								// x1, y1 - one stage correction
	    								double x1=pXY[0]+pixelCorrection[0][corrX + corrY*sensorCorrWidth]; // positive correction
	    								double y1=pXY[1]+pixelCorrection[1][corrX + corrY*sensorCorrWidth];
	    								// second step correction
	    								corrX=((int) Math.floor(x1/pixelCorrectionDecimation));
	    								corrY=((int) Math.floor(y1/pixelCorrectionDecimation));
	    								double corrDX=x1/pixelCorrectionDecimation-corrX;
	    								double corrDY=y1/pixelCorrectionDecimation-corrY;
	    								if (corrX<0)                      {corrX=0;corrDX=0;}
	    								else if (corrX>=sensorCorrWidth)  {corrX=sensorCorrWidth-1;corrDX=0;}
	    								if (corrY<0)                      {corrY=0;corrDY=0;}
	    								else if (corrY>=sensorCorrHeight) {corrY=sensorCorrHeight-1;corrDY=0;}
	    								int index00=corrX + corrY*sensorCorrWidth;
	    								int indexX0=index00+((corrX==(sensorCorrWidth-1))?0:1);
	    								int index0Y=index00+((corrY==(sensorCorrHeight-1))?0:sensorCorrWidth);
	    								int indexXY=index0Y+((corrX==(sensorCorrWidth-1))?0:1);
	    								// Linear interpolate all this.pixelCorrection[i][indexXY]
	    								double [] corr=new double [pixelCorrection.length];
	    								for (int n=0;n<pixelCorrection.length;n++){
	    									corr[n]=
	    										(1-corrDX)* (1-corrDY)* pixelCorrection[n][index00]+
	    										corrDX * (1-corrDY)* pixelCorrection[n][indexX0]+
	    										(1-corrDX)*    corrDY * pixelCorrection[n][index0Y]+
	    										corrDX *    corrDY * pixelCorrection[n][indexXY];
	    									if (debugThis){
	    										System.out.println("combineDistortionsSensorToFlat(): this.pixelCorrection["+n+"]["+index00+"]="+pixelCorrection[n][index00]);
	    										System.out.println("combineDistortionsSensorToFlat(): this.pixelCorrection["+n+"]["+indexX0+"]="+pixelCorrection[n][indexX0]);
	    										System.out.println("combineDistortionsSensorToFlat(): this.pixelCorrection["+n+"]["+index0Y+"]="+pixelCorrection[n][index0Y]);
	    										System.out.println("combineDistortionsSensorToFlat(): this.pixelCorrection["+n+"]["+indexXY+"]="+pixelCorrection[n][indexXY]);
	    									}
	    								}
	    								pXY[0]+=corr[0]; // opposite sign as in Distortions.initFittingSeries()
	    								pXY[1]+=corr[1];
	    								// map to image pixels
	    								pXY[0]=(pXY[0]-equirectangularMap.x0)/equirectangularMap.pixelStep;
	    								pXY[1]=(pXY[1]-equirectangularMap.y0)/equirectangularMap.pixelStep;
	    								if (
	    										(pXY[0]<0) ||
	    										(pXY[0]>=equirectangularMap.width) ||
	    										(pXY[1]<0) ||
	    										(pXY[1]>=equirectangularMap.height) ||
	    										((corr.length>2) && (corr[2]<equirectangularMap.minAlpha))) {
	    									pXY=null;
	    								} else {
	    									equirectangularMap.map[outIndex]=new float [corr.length];
	    									for (int n=0;n<2;n++)           equirectangularMap.map[outIndex][n]=(float) pXY[n];
	    									for (int n=2;n<corr.length;n++) equirectangularMap.map[outIndex][n]=(float) corr[n];
	    								}
	    							}
	    						}
	    						if (pXY==null) equirectangularMap.map[outIndex]=null;
	    					}
	    					//	    		IJ.showProgress(iLat, equirectangularMap.pixelsVertical-1);
	    					final int numFinished=ipLatDoneAtomic.getAndIncrement();
	    					//					IJ.showProgress(progressValues[numFinished]);
	    					SwingUtilities.invokeLater(new Runnable() {
	    						@Override
								public void run() {
	    							// Here, we can safely update the GUI
	    							// because we'll be called from the
	    							// event dispatch thread
	    							IJ.showProgress(numFinished,equirectangularMap.pixelsVertical);
	    						}
	    					});
	    				}
	    			}
	    		};
	    	}
	    	startAndJoin(threads);
	    }


	    public void combineDistortionsEquirectangularToSensorRadial( // old version without non-radial terms
	    		int    channel,
		    	double longitudeLeft0,
		    	double longitudeRight0,
		    	double latitudeTop0,
		    	double latitudeBottom0,
		    	int    pixelsHorizontal0,
	    		int width0,
	    		int height0,
	    		double x00,
	    		double y00,
	    		double pixelStep0,
	    		int maxThreads){
	    	final double maxKR=2.0; // maximal ratio to distortion radius to consider
	    	if (this.pixelCorrection==null){
	    		String msg="No pixel correction data";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
	    	}
	    	this.equirectangularMap=new EquirectangularMap(
	    			channel,
			    	longitudeLeft0,
			    	longitudeRight0,
			    	latitudeTop0,
			    	latitudeBottom0,
			    	pixelsHorizontal0,
	    			width0,
		    		height0,
		    		x00,
		    		y00,
		    		pixelStep0,
		    		this.pixelCorrection.length);
	    	final boolean debugThis=false;
			final int sensorCorrWidth=  (this.pixelCorrectionWidth-1)/this.pixelCorrectionDecimation+1;
			final int sensorCorrHeight= (this.pixelCorrectionHeight-1)/this.pixelCorrectionDecimation+1;
	    	final Matrix MA=rotateSensorCoordToPanoCoord().transpose(); // panorama coordinates to sensor coordinates (reversed rotation)
	   		final Thread[] threads = newThreadArray(maxThreads);
	   		final AtomicInteger ipLatAtomic     = new AtomicInteger(0);
	   		final AtomicInteger ipLatDoneAtomic = new AtomicInteger(1);

	    	final double [][] pixelCorrection=this.pixelCorrection;
	    	final double pixelCorrectionDecimation=this.pixelCorrectionDecimation;
	    	final double focalLength=     this.focalLength;
	    	final double px0=             this.px0;
	    	final double py0=             this.py0;
	    	final double distortionRadius=this.distortionRadius;
	    	final double pixelSize=       this.pixelSize;
	    	final double distortionA8=    this.distortionA8;
	    	final double distortionA7=    this.distortionA7;
	    	final double distortionA6=    this.distortionA6;
	    	final double distortionA5=    this.distortionA5;
	    	final double distortionA =    this.distortionA;
	    	final double distortionB =    this.distortionB;
	    	final double distortionC =    this.distortionC;

	    	final double d=1.0-this.distortionA8-this.distortionA7-this.distortionA6-this.distortionA5-this.distortionA-this.distortionB-this.distortionC;
	    	final boolean use8=(this.distortionA8!=0.0) || (this.distortionA7!=0.0) || (this.distortionA6!=0.0);
	    	final EquirectangularMap equirectangularMap=this.equirectangularMap;

//	    	for (int iLat=0;iLat<equirectangularMap.pixelsVertical;iLat++){
	    	for (int ithread = 0; ithread < threads.length; ithread++) {
	    		threads[ithread] = new Thread() {
	    			@Override
					public void run() {
	    				for (int iLat=ipLatAtomic.getAndIncrement(); iLat<equirectangularMap.pixelsVertical;iLat=ipLatAtomic.getAndIncrement()){

	    					double latitude=equirectangularMap.latitudeTop-equirectangularMap.degreesPerPixel*iLat;
	    					for (int iLong=0;iLong<equirectangularMap.pixelsHorizontal;iLong++){
	    						int outIndex=iLat*equirectangularMap.pixelsHorizontal+iLong;
	    						double [] pXY=null;
	    						double longitude=equirectangularMap.longitudeLeft+equirectangularMap.degreesPerPixel*iLong;
	    						double [][] xyz={ // unity vector in the direction
	    								{Math.cos(latitude*Math.PI/180.0)*Math.sin(longitude*Math.PI/180.0)},
	    								{Math.sin(latitude*Math.PI/180.0)},
	    								{Math.cos(latitude*Math.PI/180.0)*Math.cos(longitude*Math.PI/180.0)}};
	    						// get sensor coordinates
	    						double [] v=MA.times(new Matrix (xyz)).getRowPackedCopy();
	    						if (v[2]>0){ // potentially visible
	    							v[0]*=focalLength/v[2];
	    							v[1]*=focalLength/v[2];
	    							v[2]= focalLength;
	    							double rND=Math.sqrt(v[0]*v[0]+v[1]*v[1]);
	    							double r=rND/distortionRadius;
	    							if (r<maxKR){ // skip if they are too far from the center
	    								double k;
	    								if (use8) k=((((((distortionA8*r+distortionA7)*r+distortionA6)*r+distortionA5)*r + distortionA)*r+distortionB)*r+distortionC)*r+d;
	    								else	  k=(((distortionA5*r + distortionA)*r+distortionB)*r+distortionC)*r+d;
	    								// calculate in sensor pixel
	    								pXY=new double[2];
	    								pXY[0]= v[0]*k*1000.0/pixelSize+px0; // in pixels, right
	    								pXY[1]=-v[1]*k*1000.0/pixelSize+py0; // in pixels, down
	    								// un-apply residual correction, for now - consider it small/smooth to skip reverse mapping, make it a two step
	    								int corrX=((int) Math.floor(pXY[0]/pixelCorrectionDecimation));
	    								int corrY=((int) Math.floor(pXY[1]/pixelCorrectionDecimation));
	    								if (corrX<0)                      corrX=0;
	    								else if (corrX>=sensorCorrWidth)  corrX=sensorCorrWidth-1;
	    								if (corrY<0)                      corrY=0;
	    								else if (corrY>=sensorCorrHeight) corrY=sensorCorrHeight-1;
	    								// x1, y1 - one stage correction
	    								double x1=pXY[0]+pixelCorrection[0][corrX + corrY*sensorCorrWidth]; // positive correction
	    								double y1=pXY[1]+pixelCorrection[1][corrX + corrY*sensorCorrWidth];
	    								// second step correction
	    								corrX=((int) Math.floor(x1/pixelCorrectionDecimation));
	    								corrY=((int) Math.floor(y1/pixelCorrectionDecimation));
	    								double corrDX=x1/pixelCorrectionDecimation-corrX;
	    								double corrDY=y1/pixelCorrectionDecimation-corrY;
	    								if (corrX<0)                      {corrX=0;corrDX=0;}
	    								else if (corrX>=sensorCorrWidth)  {corrX=sensorCorrWidth-1;corrDX=0;}
	    								if (corrY<0)                      {corrY=0;corrDY=0;}
	    								else if (corrY>=sensorCorrHeight) {corrY=sensorCorrHeight-1;corrDY=0;}
	    								int index00=corrX + corrY*sensorCorrWidth;
	    								int indexX0=index00+((corrX==(sensorCorrWidth-1))?0:1);
	    								int index0Y=index00+((corrY==(sensorCorrHeight-1))?0:sensorCorrWidth);
	    								int indexXY=index0Y+((corrX==(sensorCorrWidth-1))?0:1);
	    								// Linear interpolate all this.pixelCorrection[i][indexXY]
	    								double [] corr=new double [pixelCorrection.length];
	    								for (int n=0;n<pixelCorrection.length;n++){
	    									corr[n]=
	    										(1-corrDX)* (1-corrDY)* pixelCorrection[n][index00]+
	    										corrDX * (1-corrDY)* pixelCorrection[n][indexX0]+
	    										(1-corrDX)*    corrDY * pixelCorrection[n][index0Y]+
	    										corrDX *    corrDY * pixelCorrection[n][indexXY];
	    									if (debugThis){
	    										System.out.println("combineDistortionsSensorToFlat(): this.pixelCorrection["+n+"]["+index00+"]="+pixelCorrection[n][index00]);
	    										System.out.println("combineDistortionsSensorToFlat(): this.pixelCorrection["+n+"]["+indexX0+"]="+pixelCorrection[n][indexX0]);
	    										System.out.println("combineDistortionsSensorToFlat(): this.pixelCorrection["+n+"]["+index0Y+"]="+pixelCorrection[n][index0Y]);
	    										System.out.println("combineDistortionsSensorToFlat(): this.pixelCorrection["+n+"]["+indexXY+"]="+pixelCorrection[n][indexXY]);
	    									}
	    								}
	    								pXY[0]+=corr[0]; // opposite sign as in Distortions.initFittingSeries()
	    								pXY[1]+=corr[1];
	    								// map to image pixels
	    								pXY[0]=(pXY[0]-equirectangularMap.x0)/equirectangularMap.pixelStep;
	    								pXY[1]=(pXY[1]-equirectangularMap.y0)/equirectangularMap.pixelStep;
	    								if (
	    										(pXY[0]<0) ||
	    										(pXY[0]>=equirectangularMap.width) ||
	    										(pXY[1]<0) ||
	    										(pXY[1]>=equirectangularMap.height) ||
	    										((corr.length>2) && (corr[2]<equirectangularMap.minAlpha))) {
	    									pXY=null;
	    								} else {
	    									equirectangularMap.map[outIndex]=new float [corr.length];
	    									for (int n=0;n<2;n++)           equirectangularMap.map[outIndex][n]=(float) pXY[n];
	    									for (int n=2;n<corr.length;n++) equirectangularMap.map[outIndex][n]=(float) corr[n];
	    								}
	    							}
	    						}
	    						if (pXY==null) equirectangularMap.map[outIndex]=null;

	    					}
	    					//	    		IJ.showProgress(iLat, equirectangularMap.pixelsVertical-1);
	    					final int numFinished=ipLatDoneAtomic.getAndIncrement();
	    					//					IJ.showProgress(progressValues[numFinished]);
	    					SwingUtilities.invokeLater(new Runnable() {
	    						@Override
								public void run() {
	    							// Here, we can safely update the GUI
	    							// because we'll be called from the
	    							// event dispatch thread
	    							IJ.showProgress(numFinished,equirectangularMap.pixelsVertical);
	    						}
	    					});

	    					//ipLatDoneAtomic
	    				}
	    			}
	    		};
	    	}
	    	startAndJoin(threads);

	    }

	    /**
	     * Calculate reverse distortion table - from pixel radius to non-distorted radius
	     * Rdist/R=A5*R^4+A*R^3+B*R^2+C*R+(1-A5-A-B-C)
	     * @return false if distortion is too high
	     */
	    public boolean calcReverseDistortionTable(){
	    	boolean debugThis=false; //true;
	    	double delta=1E-8;
	    	double minDerivative=0.1;
	    	int numIterations=1000;
	    	double drDistDr=1.0;
//			public double distortionA5=0.0; //r^5 (normalized to focal length or to sensor half width?)
//			public double distortionA=0.0; // r^4 (normalized to focal length or to sensor half width?)
//			public double distortionB=0.0; // r^3
//			public double distortionC=0.0; // r^2
	    	boolean use8=(this.distortionA8!=0.0) || (this.distortionA7!=0.0) || (this.distortionA6!=0.0);
	    	double d=1.0-this.distortionA8-this.distortionA7-this.distortionA6-this.distortionA5-this.distortionA-this.distortionB-this.distortionC;
	    	double rPrev=0.0;
	    	this.rByRDist=new double [(int) Math.ceil(this.maxR/this.stepR)+1];
			for (int j=1;j<this.rByRDist.length;j++) this.rByRDist[j]=Double.NaN;
			this.rByRDist[0]=1.0/d;
	    	boolean bailOut=false;
			if (debugThis)	System.out.println("calcReverseDistortionTable()");

	    	for (int i=1;i<this.rByRDist.length;i++) {
	    		double rDist=this.stepR*i;
	    		double r=rPrev+this.stepR/drDistDr;
//				if (debugThis)	System.out.println("calcReverseDistortionTable() i="+i+" rDist="+rDist+" r="+r+" rPrev="+rPrev);

	    		for (int iteration=0;iteration<numIterations;iteration++){
	    			double k;
	    			if (use8){
	    				k=(((((((this.distortionA8)*r+this.distortionA7)*r+this.distortionA6)*r+this.distortionA5)*r + this.distortionA)*r+this.distortionB)*r+this.distortionC)*r+d;
		    			drDistDr=(((((((8*this.distortionA8)*r + 7*this.distortionA7)*r + 6*this.distortionA6)*r + 5*this.distortionA5)*r + 4*this.distortionA)*r+3*this.distortionB)*r+2*this.distortionC)*r+d;
	    			} else {
	    				k=(((this.distortionA5*r + this.distortionA)*r+this.distortionB)*r+this.distortionC)*r+d;
		    			drDistDr=(((5*this.distortionA5*r + 4*this.distortionA)*r+3*this.distortionB)*r+2*this.distortionC)*r+d;
	    			}
	    			double rD=r*k;
	    			if (drDistDr<minDerivative) {
	    				bailOut=true;
	    				break; // too high distortion
	    			}
	    			if (Math.abs(rD-rDist)<delta) break; // success
	    			r+=(rDist-rD)/drDistDr;
	    		}
	    		if (bailOut) {
					if (debugThis)	System.out.println("calcReverseDistortionTable() i="+i+" Bailing out, drDistDr="+drDistDr);
	    			return false;
	    		}
	    		rPrev=r;
	    		this.rByRDist[i]=r/rDist;
				if (debugThis)	System.out.println("calcReverseDistortionTable() i="+i+" rDist="+rDist+" r="+r+" rPrev="+rPrev+" this.rByRDist[i]="+this.rByRDist[i]);
	    	}
	    	return true;
	    }

	    public double getRByRDist(double rDist, boolean debug){
	    	// add exceptions;
	    	if (this.rByRDist==null) {
	    		if (debug)System.out.println("getRByRDist("+IJ.d2s(rDist,3)+"): this.rByRDist==null");
	    		return Double.NaN;
	    	}
	    	if (rDist<0) {
	    		if (debug)System.out.println("getRByRDist("+IJ.d2s(rDist,3)+"): rDist<0");
	    		return Double.NaN;
	    	}
	    	int index=(int) Math.floor(rDist/this.stepR);
	    	if (index>=(this.rByRDist.length-1)) {
	    		if (debug) System.out.println("getRByRDist("+IJ.d2s(rDist,3)+"): index="+index+">="+(this.rByRDist.length-1));
	    		return Double.NaN;
	    	}
	    	double result=this.rByRDist[index]+(this.rByRDist[index+1]-this.rByRDist[index])*(rDist/this.stepR-index);
	    	if (Double.isNaN(result)){
	    		if (debug) System.out.println("this.rByRDist["+index+"]="+this.rByRDist[index]);
	    		if (debug) System.out.println("this.rByRDist["+(index+1)+"]="+this.rByRDist[index+1]);
	    		if (debug) System.out.println("rDist="+rDist);
	    		if (debug) System.out.println("(rDist/this.stepR="+(rDist/this.stepR));

	    	}
	    	return result;

	    }

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
