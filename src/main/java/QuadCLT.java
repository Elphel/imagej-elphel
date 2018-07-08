/**
 **
 ** QuadCLT - Process images with CLT-based methods (code specific to ImageJ plugin)
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  QuadCLT.java is free software: you can redistribute it and/or modify
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

//import java.awt.Polygon;
import java.awt.Rectangle;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.attribute.PosixFilePermission;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;

import ij.CompositeImage;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
//import ij.gui.Overlay;
import ij.io.FileSaver;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;


public class QuadCLT {
	static String []                                       fine_corr_coeff_names = {"A","B","C","D","E","F"};
	static String []                                       fine_corr_dir_names = {"X","Y"};
	public static String                                   PREFIX =     "EYESIS_DCT.";    // change later (first on save)
	public static String                                   PREFIX_AUX = "EYESIS_DCT_AUX."; // change later (first on save)
	static int                                             QUAD =  4; // number of cameras
	public Properties                                      properties = null;
//	public String                                          properties_prefix = "EYESIS_DCT.";
	public EyesisCorrections                               eyesisCorrections = null;
	public EyesisCorrectionParameters.CorrectionParameters correctionsParameters=null;
	double [][][][][][]                                    clt_kernels = null;
	GeometryCorrection                                     geometryCorrection = null;
	double []                                              extrinsic_corr = new double [GeometryCorrection.CORR_NAMES.length]; // extrinsic corrections (needed from properties, before geometryCorrection
	public int                                             extra_items = 8; // number of extra items saved with kernels (center offset (partial, full, derivatives)
	public ImagePlus                                       eyesisKernelImage = null;
	public long                                            startTime;     // start of batch processing
	public long                                            startSetTime;  // start of set processing
	public long                                            startStepTime; // start of step processing

	public double [][][]                                   fine_corr  = new double [4][2][6]; // per port, per x/y, set of 6 coefficient for fine geometric corrections

	TileProcessor                                          tp = null;

	String                                                 image_name = null;
	double []                                              gps_lla =    null;
	double [][][]                                          image_data = null;
    boolean [][]                                           saturation_imp = null; // (near) saturated pixels or null
    boolean                                                is_aux = false;


// magic scale should be set before using  TileProcessor (calculated disparities depend on it)
    public void resetGroundTruthByRig() {
    	tp.rig_disparity_strength = null;
    }
    public double [][] getGroundTruthByRig(){
    	if (tp == null) return null;
    	return tp.rig_disparity_strength;
    }
	public void setTiles (ImagePlus imp, // set tp.tilesX, tp.tilesY
			EyesisCorrectionParameters.CLTParameters    clt_parameters,
			int threadsMax
			){
		setTiles(clt_parameters,
				imp.getWidth()/clt_parameters.transform_size,
				imp.getHeight()/clt_parameters.transform_size,
				threadsMax);
		/*
		if (tp == null){
			tp = new TileProcessor(imp.getWidth()/clt_parameters.transform_size,
					imp.getHeight()/clt_parameters.transform_size,
					clt_parameters.transform_size,
					clt_parameters.stSize,
					clt_parameters.corr_magic_scale,
					clt_parameters.grow_disp_trust,
					clt_parameters.max_overexposure, // double maxOverexposure,
					threadsMax);
		}
		*/
	}

	public void setTiles (
			EyesisCorrectionParameters.CLTParameters    clt_parameters,
			int tilesX,
			int tilesY,
			int threadsMax
			){
		if (tp == null){
			tp = new TileProcessor(
					tilesX,
					tilesY,
					clt_parameters.transform_size,
					clt_parameters.stSize,
					clt_parameters.corr_magic_scale,
					clt_parameters.grow_disp_trust,
					clt_parameters.max_overexposure, // double maxOverexposure,
					threadsMax);
		}
	}

// used for aux camera
	public boolean setupImageData(
			String image_name,
			String [] sourceFiles,
			EyesisCorrectionParameters.CLTParameters       clt_parameters,
			int threadsMax,
			int debugLevel) {
		  QuadCLT.SetChannels [] set_channels_aux =  setChannels(image_name, debugLevel);
		  if ((set_channels_aux == null) || (set_channels_aux.length==0)) {
			  System.out.println("No files for the auxiliary camera match series "+image_name);
			  return false;
		  }
		  double [] referenceExposures_aux =  eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		  int [] channelFiles_aux =  set_channels_aux[0].fileNumber();
		  // make single
		  boolean [][] saturation_aux =  (clt_parameters.sat_level > 0.0)? new boolean[channelFiles_aux.length][] : null;
		  double [] scaleExposures_aux =  new double[channelFiles_aux.length];
		  ImagePlus [] imp_srcs_aux = conditionImageSet(
				  clt_parameters,               // EyesisCorrectionParameters.CLTParameters  clt_parameters,
				  sourceFiles,                  // String []                                 sourceFiles,
				  image_name,      // set_channels_aux[0].name(), // String                                    set_name,
				  referenceExposures_aux,       // double []                                 referenceExposures,
				  channelFiles_aux,             // int []                                    channelFiles,
				  scaleExposures_aux,           //output  // double [] scaleExposures
				  saturation_aux,               //output  // boolean [][]                              saturation_imp,
				  debugLevel); // int                                       debugLevel);

		  double [][][] double_stacks_aux = new double [imp_srcs_aux.length][][];
		  for (int i = 0; i < double_stacks_aux.length; i++){
			  double_stacks_aux[i] = eyesisCorrections.bayerToDoubleStack(
					  imp_srcs_aux[i], // source Bayer image, linearized, 32-bit (float))
					  null); // no margins, no oversample
		  }

		  for (int i = 0; i < double_stacks_aux.length; i++){
			  for (int j =0 ; j < double_stacks_aux[i][0].length; j++){
				  double_stacks_aux[i][2][j]*=0.5; // Scale green 0.5 to compensate more pixels than R,B
			  }
		  }
		  setTiles (imp_srcs_aux[0], // set global tp.tilesX, tp.tilesY
				  clt_parameters,
				  threadsMax);
		  this.image_name =      image_name;
		  image_data =      double_stacks_aux;
		  saturation_imp =  saturation_aux;
		  tp.setTrustedCorrelation(clt_parameters.grow_disp_trust);
		  tp.resetCLTPasses();
		  return true;

	}


	public QuadCLT(
			String                                          prefix,
			Properties                                      properties,
			EyesisCorrections                               eyesisCorrections,
			EyesisCorrectionParameters.CorrectionParameters correctionsParameters
			){
		this.eyesisCorrections=      eyesisCorrections;
		this.correctionsParameters = correctionsParameters;
		this.properties =            properties;
		is_aux =                     prefix.equals(PREFIX_AUX);
//		this.properties_prefix =     prefix;
//		System.out.println("new QuadCLT(), prefix = "+prefix);
		getProperties(prefix);
	}

	// TODO:Add saving just calibration

//	public void setProperties(){
//		setProperties(this.properties_prefix);
//	}
	public void setProperties(String prefix, Properties properties){ // save
		if (properties == null) {
			properties = this.properties;
		}
//		System.out.println("setProperties("+prefix+")");
		for (int n = 0; n < fine_corr.length; n++){
			for (int d = 0; d < fine_corr[n].length; d++){
				for (int i = 0; i < fine_corr[n][d].length; i++){
					String name = prefix+"fine_corr_"+n+fine_corr_dir_names[d]+fine_corr_coeff_names[i];
					properties.setProperty(name,   this.fine_corr[n][d][i]+"");
				}
			}
		}
		GeometryCorrection gc = geometryCorrection;
		if (gc == null) { // if it was not yet created
			gc = new GeometryCorrection(this.extrinsic_corr);
		}
		for (int i = 0; i < GeometryCorrection.CORR_NAMES.length; i++){
			String name = prefix+"extrinsic_corr_"+GeometryCorrection.CORR_NAMES[i];
			properties.setProperty(name,  gc.getCorrVector().toArray()[i]+"");
//			System.out.println("setProperties():"+i+": setProperty("+name+","+gc.getCorrVector().toArray()[i]+"");
		}
		if (is_aux && (gc.rigOffset != null)) {
			gc.rigOffset.setProperties(prefix,properties);
		}
	}


	public void copyPropertiesFrom(Properties other_properties, String other_prefix, String this_prefix){ // save
//		System.out.println("copyPropertiesFrom(other_properties, "+other_prefix+", this_prefix"+")");
		for (int n = 0; n < fine_corr.length; n++){
			for (int d = 0; d < fine_corr[n].length; d++){
				for (int i = 0; i < fine_corr[n][d].length; i++){
					String other_name = other_prefix+"fine_corr_"+n+fine_corr_dir_names[d]+fine_corr_coeff_names[i];
					String this_name = this_prefix+"fine_corr_"+n+fine_corr_dir_names[d]+fine_corr_coeff_names[i];
		  			if (other_properties.getProperty(other_name)!=null) {
		  				this.fine_corr[n][d][i]=Double.parseDouble(other_properties.getProperty(other_name));
						properties.setProperty(this_name,   this.fine_corr[n][d][i]+"");
		  			}
				}
			}
		}
		GeometryCorrection gc = geometryCorrection;
		if (gc == null) { // if it was not yet created
			gc = new GeometryCorrection(this.extrinsic_corr);
		}
		for (int i = 0; i < GeometryCorrection.CORR_NAMES.length; i++){
			String other_name = other_prefix+"extrinsic_corr_"+GeometryCorrection.CORR_NAMES[i];
  			if (other_properties.getProperty(other_name)!=null) {
  				this.extrinsic_corr[i] = Double.parseDouble(other_properties.getProperty(other_name));
  				if (geometryCorrection != null){
  					geometryCorrection.getCorrVector().toArray()[i] = this.extrinsic_corr[i];
  				}
  			}
			String this_name =  this_prefix+"extrinsic_corr_"+GeometryCorrection.CORR_NAMES[i];
			properties.setProperty(this_name,  gc.getCorrVector().toArray()[i]+"");
//			System.out.println("copyPropertiesFrom():"+i+": setProperty("+this_name+","+gc.getCorrVector().toArray()[i]+"");
		}
//		System.out.println("Done copyPropertiesFrom");
	}

	public GeometryCorrection  getGeometryCorrection() {
		return geometryCorrection;
	}
	public double [][][][][][] getCLTKernels(){
		return clt_kernels;
	}
	public void listGeometryCorrection(boolean full){
		GeometryCorrection gc = geometryCorrection;
		if (gc == null) { // if it was not yet created
			gc = new GeometryCorrection(this.extrinsic_corr);
		}
		gc.listGeometryCorrection(full);
	}

	public void getProperties(String prefix){ // restore
//		System.out.println("getProperties("+prefix+")");
		for (int n = 0; n < fine_corr.length; n++){
			for (int d = 0; d < fine_corr[n].length; d++){
				for (int i = 0; i < fine_corr[n][d].length; i++){
					String name = prefix+"fine_corr_"+n+fine_corr_dir_names[d]+fine_corr_coeff_names[i];
		  			if (properties.getProperty(name)!=null) this.fine_corr[n][d][i]=Double.parseDouble(properties.getProperty(name));
				}
			}
		}
		// always set extrinsic_corr
		for (int i = 0; i < GeometryCorrection.CORR_NAMES.length; i++){
			String name = prefix+"extrinsic_corr_"+GeometryCorrection.CORR_NAMES[i];
  			if (properties.getProperty(name)!=null) {
  				this.extrinsic_corr[i] = Double.parseDouble(properties.getProperty(name));
//  				System.out.println("getProperties():"+i+": getProperty("+name+") -> "+properties.getProperty(name)+"");

  				if (geometryCorrection != null){
  					geometryCorrection.getCorrVector().toArray()[i] = this.extrinsic_corr[i];
  				}
  			}
		}
//		if (is_aux && (geometryCorrection != null)) {
//			geometryCorrection.setRigOffsetFromProperies(prefix, properties);
//		}
		if (geometryCorrection == null) {
			geometryCorrection = new GeometryCorrection(this.extrinsic_corr);
		}

		if (is_aux) {
			geometryCorrection.setRigOffsetFromProperies(prefix, properties);
		}
	}

	public void setKernelImageFile(ImagePlus img_kernels){
		eyesisKernelImage = img_kernels;
	}

	public boolean kernelImageSet(){
		return eyesisKernelImage != null;
	}

	public boolean CLTKernelsAvailable(){
		return clt_kernels != null;
	}
	public boolean geometryCorrectionAvailable(){
		return (geometryCorrection != null) && geometryCorrection.isInitialized();
	}
	public boolean initGeometryCorrection(int debugLevel){
		// keep rig offsets if edited
		if (geometryCorrection == null) {
			geometryCorrection = new GeometryCorrection(extrinsic_corr);
		}
		PixelMapping.SensorData [] sensors =  eyesisCorrections.pixelMapping.sensors;
		// verify that all sensors have the same distortion parameters
		int numSensors = sensors.length;
		for (int i = 1; i < numSensors; i++){
			if (	(sensors[0].focalLength !=           sensors[i].focalLength) ||
					(sensors[0].distortionC !=           sensors[i].distortionC) ||
					(sensors[0].distortionB !=           sensors[i].distortionB) ||
					(sensors[0].distortionA !=           sensors[i].distortionA) ||
					(sensors[0].distortionA5 !=          sensors[i].distortionA5) ||
					(sensors[0].distortionA6 !=          sensors[i].distortionA6) ||
					(sensors[0].distortionA7 !=          sensors[i].distortionA7) ||
					(sensors[0].distortionA8 !=          sensors[i].distortionA8) ||
					(sensors[0].distortionRadius !=      sensors[i].distortionRadius) ||
					(sensors[0].pixelCorrectionWidth !=  sensors[i].pixelCorrectionWidth) ||
					(sensors[0].pixelCorrectionHeight != sensors[i].pixelCorrectionHeight) ||
					(sensors[0].pixelSize !=             sensors[i].pixelSize)){
				System.out.println("initGeometryCorrection(): All sensors have to have the same distortion model, but channels 0 and "+i+" mismatch");
				return false;
			}
		}
		// set common distportion parameters
		geometryCorrection.setDistortion(
				sensors[0].focalLength,
				sensors[0].distortionC,
				sensors[0].distortionB,
				sensors[0].distortionA,
				sensors[0].distortionA5,
				sensors[0].distortionA6,
				sensors[0].distortionA7,
				sensors[0].distortionA8,
				sensors[0].distortionRadius,
				sensors[0].pixelCorrectionWidth,   // virtual camera center is at (pixelCorrectionWidth/2, pixelCorrectionHeight/2)
				sensors[0].pixelCorrectionHeight,
				sensors[0].pixelSize);
		// set other/individual sensor parameters
		for (int i = 1; i < numSensors; i++){
			if (	(sensors[0].theta !=                 sensors[i].theta) || // elevation
					(sensors[0].heading !=               sensors[i].heading)){
				System.out.println("initGeometryCorrection(): All sensors have to have the same elevation and heading, but channels 0 and "+i+" mismatch");
				return false;
			}
		}
		double []   forward = new double[numSensors];
		double []   right =   new double[numSensors];
		double []   height =  new double[numSensors];
		double []   roll =    new double[numSensors];
		double [][] pXY0 =    new double[numSensors][2];
		for (int i = 0; i < numSensors; i++){
			forward[i] = sensors[i].forward;
			right[i] =   sensors[i].right;
			height[i] =  sensors[i].height;
			roll[i] =    sensors[i].psi;
			pXY0[i][0] = sensors[i].px0;
			pXY0[i][1] = sensors[i].py0;
		}
		geometryCorrection.setSensors(
				numSensors,
				sensors[0].theta,
				sensors[0].heading,
				forward,
				right,
				height,
				roll,
				pXY0);
		geometryCorrection.planeProjectLenses(); // project all lenses to the common plane

		// calcualte reverse distortion as a table to be linear intr4epolated
		geometryCorrection.calcReverseDistortionTable();

		if (numSensors == 4){
			geometryCorrection.adustSquare();
			System.out.println("Adjusted camera to orient X Y along the sides of a square");
		} else {
			System.out.println("============= Cannot adustSquare() as it requires exactly 4 sensors, "+numSensors+" provided ==========");
			return false;
		}
		// Print parameters
		if (debugLevel > 0){
			geometryCorrection.listGeometryCorrection(debugLevel > 1);
		}

//listGeometryCorrection
		return true;
	}

//GeometryCorrection

	  public double [][][][][] calculateCLTKernel ( // per color/per tileY/ per tileX/per quadrant (plus offset as 5-th)/per pixel
			  final PixelMapping.SensorData sensor, // to calculate extra shift
			  final ImageStack kernelStack,  // first stack with 3 colors/slices convolution kernels
			  final int          kernelSize, // 64
			  final EyesisCorrectionParameters.CLTParameters clt_parameters,

			  final int          threadsMax,  // maximal number of threads to launch
			  final boolean      updateStatus,
			  final int          globalDebugLevel) // update status info
	  {
		  if (kernelStack==null) return null;
		  final int kernelWidth=kernelStack.getWidth();
		  final int kernelNumHor=kernelWidth/kernelSize;
		  final int kernelNumVert=kernelStack.getHeight()/kernelSize;
		  final int nChn=kernelStack.getSize();
		  final int dtt_size =      clt_parameters.transform_size;
		  final int dtt_len = dtt_size* dtt_size;
		  final double [][][][][] clt_kernels = new double [nChn][kernelNumVert][kernelNumHor][5][];
		  for (int chn = 0; chn < nChn; chn++){
			  for (int tileY = 0; tileY < kernelNumVert ; tileY++){
				  for (int tileX = 0; tileX < kernelNumHor ; tileX++){
					  for (int n = 0; n<4; n++){
						  clt_kernels[chn][tileY][tileX][n] = new double [dtt_len];
					  }
					  clt_kernels[chn][tileY][tileX][4] = new double [extra_items];
				  }
			  }
		  }
		  // currently each 64x64 kernel corresponds to 16x16 original pixels tile, 2 tiles margin each side
		  final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		  final AtomicInteger ai = new AtomicInteger(0);
		  final int numberOfKernels=     kernelNumHor*kernelNumVert*nChn;
		  final int numberOfKernelsInChn=kernelNumHor*kernelNumVert;
		  final double [] norm_sym_weights = clt_parameters.norm_kern ? new double [dtt_size*dtt_size]:null;
		  if (norm_sym_weights != null) {
			  for (int i = 0; i < dtt_size; i++){
				  for (int j = 0; j < dtt_size; j++){
					  norm_sym_weights[i*dtt_size+j] = Math.cos(Math.PI*i/(2*dtt_size))*Math.cos(Math.PI*j/(2*dtt_size));
				  }
			  }
		  }

		  final long startTime = System.nanoTime();
		  System.out.println("calculateCLTKernel():numberOfKernels="+numberOfKernels);
		  for (int ithread = 0; ithread < threads.length; ithread++) {
			  threads[ithread] = new Thread() {
				  @Override
				public void run() {
					  float [] kernelPixels= null; // will be initialized at first use NOT yet?
					  double [] kernel=      new double[kernelSize*kernelSize];
					  int centered_len = (2*dtt_size-1) * (2*dtt_size-1);
					  double [] kernel_centered = new double [centered_len + extra_items];
					  ImageDtt image_dtt = new ImageDtt();
					  int chn,tileY,tileX;
					  DttRad2 dtt = new DttRad2(dtt_size);
					  showDoubleFloatArrays sdfa_instance = null;
					  if (globalDebugLevel > -1) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?

					  for (int nTile = ai.getAndIncrement(); nTile < numberOfKernels; nTile = ai.getAndIncrement()) {
						  chn=nTile/numberOfKernelsInChn;
						  tileY =(nTile % numberOfKernelsInChn)/kernelNumHor;
						  tileX = nTile % kernelNumHor;
						  if (tileX==0) {
							  if (updateStatus) IJ.showStatus("Processing kernels, channel "+(chn+1)+" of "+nChn+", row "+(tileY+1)+" of "+kernelNumVert);
							  if (globalDebugLevel>2) System.out.println("Processing kernels, channel "+(chn+1)+" of "+nChn+", row "+(tileY+1)+" of "+kernelNumVert+" : "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
						  }
						  kernelPixels=(float[]) kernelStack.getPixels(chn+1);

						  /* read convolution kernel */
						  extractOneKernel(
								  kernelPixels, //  array of combined square kernels, each
								  kernel, // will be filled, should have correct size before call
								  kernelNumHor, // number of kernels in a row
								  tileX, // horizontal number of kernel to extract
								  tileY); // vertical number of kernel to extract
						  if ((globalDebugLevel > 0) && (tileY == clt_parameters.tileY/2)  && (tileX == clt_parameters.tileX/2)) {
							  int length=kernel.length;
							  int size=(int) Math.sqrt(length);
							  double s =0.0;
							  for (int i=0;i<kernel.length;i++) s+=kernel[i];
							  System.out.println("calculateCLTKernel(): sum(kernel_raw)="+s);
							  if (globalDebugLevel > 1) sdfa_instance.showArrays(
									  kernel,
									  size,
									  size,
									  "raw_kernel-"+chn+"-X"+(clt_parameters.tileX/2)+"-Y"+(clt_parameters.tileY/2));

						  }

						  // now has 64x64
						  image_dtt.clt_convert_double_kernel( // converts double resolution kernel
								  kernel,          // double []   src_kernel, //
								  kernel_centered, // double []   dst_kernel, // should be (2*dtt_size-1) * (2*dtt_size-1) + extra_items size - kernel and dx, dy to the nearest 1/2 pixels
								                   // also actual full center shifts in sensor pixels
								  kernelSize,      // int src_size, // 64
								  dtt_size);       // 8
						  if ((globalDebugLevel > 0) && (tileY == clt_parameters.tileY/2)  && (tileX == clt_parameters.tileX/2)) {
							  int length=kernel_centered.length;
							  int size=(int) Math.sqrt(length);
							  double s =0.0;
							  int klen = (2*dtt_size-1) * (2*dtt_size-1);
							  for (int i = 0; i < klen; i++) s += kernel_centered[i];
							  System.out.println("calculateCLTKernel(): sum(kernel_centered)="+s);
							  if (globalDebugLevel > 1) sdfa_instance.showArrays(
									  kernel_centered,
									  size,
									  size,
									  "kernel_centered-"+chn+"-X"+(clt_parameters.tileX/2)+"-Y"+(clt_parameters.tileY/2));
						  }

						  if (norm_sym_weights != null) {
							  image_dtt.clt_normalize_kernel( //
									  kernel_centered, // double []   kernel, // should be (2*dtt_size-1) * (2*dtt_size-1) + extra_items size (last (2*dtt_size-1) are not modified)
									  norm_sym_weights, // double []   window, // normalizes result kernel * window to have sum of elements == 1.0
									  dtt_size,
									  (globalDebugLevel > 0) && (tileY == clt_parameters.tileY/2)  && (tileX == clt_parameters.tileX/2)); // 8
							  if ((globalDebugLevel > 0) && (tileY == clt_parameters.tileY/2)  && (tileX == clt_parameters.tileX/2)) {
								  int length=kernel_centered.length;
								  int size=(int) Math.sqrt(length);
								  double s =0.0;
								  int klen = (2*dtt_size-1) * (2*dtt_size-1);
								  for (int i = 0; i < klen; i++) s += kernel_centered[i];
								  System.out.println("calculateCLTKernel(): sum(kernel_normalized)="+s);
								  if (globalDebugLevel > 1) sdfa_instance.showArrays(
										  kernel_centered,
										  size,
										  size,
										  "kernel_normalized-"+chn+"-X"+(clt_parameters.tileX/2)+"-Y"+(clt_parameters.tileY/2));

							  }
						  }
						  image_dtt.clt_symmetrize_kernel( //
								  kernel_centered, // double []     kernel,      // should be (2*dtt_size-1) * (2*dtt_size-1) +4 size (last 4 are not modified)
								  clt_kernels[chn][tileY][tileX], // 	double [][]   sym_kernels, // set of 4 SS, AS, SA, AA kdernels, each dtt_size * dtt_size (may have 5-th with center shift
								  dtt_size); // 8
						  for (int i = 0; i < extra_items; i++){
							  clt_kernels[chn][tileY][tileX][4][i] = kernel_centered [centered_len + i];
						  }
						  if ((globalDebugLevel > 0) && (tileY == clt_parameters.tileY/2)  && (tileX == clt_parameters.tileX/2)) {
							  double [][] dbg_clt = {
									  clt_kernels[chn][tileY][tileX][0],
									  clt_kernels[chn][tileY][tileX][1],
									  clt_kernels[chn][tileY][tileX][2],
									  clt_kernels[chn][tileY][tileX][3]};
							  String [] titles = {"CC", "SC", "CS", "SS"};
							  int length=dbg_clt[0].length;
							  int size=(int) Math.sqrt(length);
							  if (globalDebugLevel > 1) sdfa_instance.showArrays(
									  dbg_clt,
									  size,
									  size,
									  true,
									  "pre_clt_kernels-"+chn,
									  titles);
						  }
						  image_dtt.clt_dtt3_kernel( //
								  clt_kernels[chn][tileY][tileX], // double [][]   kernels, // set of 4 SS, AS, SA, AA kdernels, each dtt_size * dtt_size (may have 5-th with center shift
								  dtt_size, // 8
								  dtt);

						  if ((globalDebugLevel > 0) && (tileY == clt_parameters.tileY/2)  && (tileX == clt_parameters.tileX/2)) {
							  System.out.println("calculateCLTKernel() - before corr: chn="+chn+" "+
									  "tileX = "+clt_parameters.tileX+" ("+(clt_parameters.tileX/2)+") "+
									  "tileY = "+clt_parameters.tileY+" ("+(clt_parameters.tileY/2)+") "+
									  "center_x = "+clt_kernels[chn][tileY][tileX][4][0]+", "+
									  "center_y = "+clt_kernels[chn][tileY][tileX][4][1]+", "+
									  "full_dx =  "+clt_kernels[chn][tileY][tileX][4][2]+", "+
									  "full_dy =  "+clt_kernels[chn][tileY][tileX][4][3]);
						  }
						  // Add sensor geometry correction (optional?)
						  // Kernel center in pixels
						  double kpx0 = (tileX -1 +0.5) *  clt_parameters.kernel_step;
						  double kpy0 = (tileY -1 +0.5) *  clt_parameters.kernel_step;
						  double [] corrPxPy = sensor.interpolateCorrectionVector(false, kpx0, kpy0);
						  image_dtt.offsetKernelSensor(
								  clt_kernels[chn][tileY][tileX], // double [][] clt_tile, // clt tile, including [4] - metadata
								  //								  (corrPxPy[0] - kpx0),  //	double dx,
								  //								  (corrPxPy[1] - kpy0)); // double dy)
								  -corrPxPy[0],  //	double dx,
								  -corrPxPy[1]); // double dy)

						  if ((globalDebugLevel > 0) && (tileY == clt_parameters.tileY/2)  && (tileX == clt_parameters.tileX/2)) {
							  System.out.println("calculateCLTKernel() - after  corr: chn="+chn+" "+
									  "tileX = "+clt_parameters.tileX+" ("+(clt_parameters.tileX/2)+") "+
									  "tileY = "+clt_parameters.tileY+" ("+(clt_parameters.tileY/2)+") "+
									  "center_x = "+clt_kernels[chn][tileY][tileX][4][0]+", "+
									  "center_y = "+clt_kernels[chn][tileY][tileX][4][1]+", "+
									  "full_dx =  "+clt_kernels[chn][tileY][tileX][4][2]+", "+
									  "full_dy =  "+clt_kernels[chn][tileY][tileX][4][3]);
							  System.out.println("calculateCLTKernel() - after  corr: chn="+chn+" "+
									  "kpx0 = "+kpx0+
									  "kpy0 = "+kpy0+
									  "corrPxPy[0] = "+corrPxPy[0]+
									  "corrPxPy[1] = "+corrPxPy[1]);
						  }

						  if ((globalDebugLevel > 0) && (tileY == clt_parameters.tileY/2)  && (tileX == clt_parameters.tileX/2)) {
							  double [][] dbg_clt = {
									  clt_kernels[chn][tileY][tileX][0],
									  clt_kernels[chn][tileY][tileX][1],
									  clt_kernels[chn][tileY][tileX][2],
									  clt_kernels[chn][tileY][tileX][3]};
							  String [] titles = {"CC", "SC", "CS", "SS"};
							  int length=dbg_clt[0].length;
							  int size=(int) Math.sqrt(length);
							  if (globalDebugLevel > 1) sdfa_instance.showArrays(
									  dbg_clt,
									  size,
									  size,
									  true,
									  "dbg_clt_kernels-"+chn,
									  titles);
							  System.out.println("calculateCLTKernel() chn="+chn+" "+
									  "tileX = "+clt_parameters.tileX+" ("+(clt_parameters.tileX/2)+") "+
									  "tileY = "+clt_parameters.tileY+" ("+(clt_parameters.tileY/2)+") "+
									  "center_x = "+clt_kernels[chn][tileY][tileX][4][0]+", "+
									  "center_y = "+clt_kernels[chn][tileY][tileX][4][1]+", "+
									  "full_dx =  "+clt_kernels[chn][tileY][tileX][4][2]+", "+
									  "full_dy =  "+clt_kernels[chn][tileY][tileX][4][3]);
						  }
					  }
				  }
			  };
		  }
		  ImageDtt.startAndJoin(threads);
		  if (globalDebugLevel > 1) System.out.println("Threads done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		  System.out.println("1.Threads done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		  // Calculate differential offsets to interpolate for tiles between kernel centers
		  ImageDtt image_dtt = new ImageDtt();
		  image_dtt.clt_fill_coord_corr(
				  clt_parameters.kernel_step,  //  final int             kern_step, // distance between kernel centers, in pixels.
				  clt_kernels,                 // final double [][][][] clt_data,
				  threadsMax,                  // maximal number of threads to launch
				  globalDebugLevel);
		  return clt_kernels;
	  }

	  public double [][] flattenCLTKernels (      // per color, save 4 kernelas and displacement as (2*dtt_size+1)*(2*dtt_size) tiles in an image (last row - 4 values shift x,y)
			  final double [][][][][] kernels,    // per color/per tileY/ per tileX/per quadrant (plus offset as 5-th)/per pixel
			  final int          threadsMax,      // maximal number of threads to launch
			  final boolean      updateStatus,
			  final int          globalDebugLevel) // update status info
	  {
		  if (kernels==null) return null;
		  final int nChn = kernels.length;
		  final int kernelNumVert=kernels[0].length;
		  final int kernelNumHor=kernels[0][0].length;
		  final int dtt_len = kernels[0][0][0][0].length;
		  final int dtt_size =      (int) Math.sqrt(dtt_len);
		  final int tileWidth =  2 * dtt_size;
		  final int tileHeight = 2 * dtt_size + 1; // last row - shift with 0.5 pix steps
		  final int width =  tileWidth *  kernelNumHor;
		  final int height = tileHeight * kernelNumVert;
		  final double [][] clt_flat = new double [nChn][width * height];
		  // currently each 64x64 kernel corresponds to 16x16 original pixels tile, 2 tiles margin each side
		  final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		  final AtomicInteger ai = new AtomicInteger(0);
		  final int numberOfKernels=     kernelNumHor*kernelNumVert*nChn;
		  final int numberOfKernelsInChn=kernelNumHor*kernelNumVert;

		  final long startTime = System.nanoTime();
		  System.out.println("flattenCLTKernels():numberOfKernels="+numberOfKernels);
		  for (int ithread = 0; ithread < threads.length; ithread++) {
			  threads[ithread] = new Thread() {
				  @Override
				public void run() {
					  int chn,tileY,tileX;
					  for (int nTile = ai.getAndIncrement(); nTile < numberOfKernels; nTile = ai.getAndIncrement()) {
						  chn=nTile/numberOfKernelsInChn;
						  tileY =(nTile % numberOfKernelsInChn)/kernelNumHor;
						  tileX = nTile % kernelNumHor;
						  for (int i = 0; i < dtt_size; i++){
							System.arraycopy(
									kernels[chn][tileY][tileX][0],
									i * dtt_size,
									clt_flat[chn],
									(tileY*tileHeight + i)            * width + (tileX * tileWidth),
									dtt_size);
							System.arraycopy(
									kernels[chn][tileY][tileX][1],
									i * dtt_size,
									clt_flat[chn],
									(tileY*tileHeight + i)            * width + (tileX * tileWidth) + dtt_size,
									dtt_size);

							System.arraycopy(
									kernels[chn][tileY][tileX][2],
									i * dtt_size,
									clt_flat[chn],
									(tileY*tileHeight + i + dtt_size) * width + (tileX * tileWidth),
									dtt_size);
							System.arraycopy(
									kernels[chn][tileY][tileX][3],
									i * dtt_size,
									clt_flat[chn],
									(tileY*tileHeight + i + dtt_size) * width + (tileX * tileWidth) + 1 * dtt_size,
									dtt_size);
						  }
						  System.arraycopy(
								  kernels[chn][tileY][tileX][4], // just 2 values
								  0,
								  clt_flat[chn],
								  (tileY*tileHeight + 2 * dtt_size) * width + (tileX * tileWidth),
								  extra_items);
					  }
				  }
			  };
		  }
		  ImageDtt.startAndJoin(threads);
		  if (globalDebugLevel > 1) System.out.println("Threads done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		  System.out.println("1.Threads done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		  /* prepare result stack to return */
		  return clt_flat;
	  }

	  public void showCLTKernels(
			  final int          threadsMax,      // maximal number of threads to launch
			  final boolean      updateStatus,
			  final int          globalDebugLevel) // update status info
	  {
		  for (int chn=0;chn < clt_kernels.length; chn++){
			  if (clt_kernels[chn]!=null){
				  //					System.out.println("showKernels("+chn+")");
				  showCLTKernels(
						  chn,
						  threadsMax,
						  updateStatus,
						  globalDebugLevel);
			  }
		  }
	  }

	  public void showCLTKernels(
			  int chn,
			  final int          threadsMax,      // maximal number of threads to launch
			  final boolean      updateStatus,
			  final int          globalDebugLevel) // update status info
	  {
		  double [][] flat_kernels = flattenCLTKernels (      // per color, save 4 kernelas and displacement as (2*dtt_size+1)*(2*dtt_size) tiles in an image (last row - shift x,y)
				  clt_kernels[chn],     // per color/per tileY/ per tileX/per quadrant (plus offset as 5-th)/per pixel
				  threadsMax,      // maximal number of threads to launch
				  updateStatus,
				  globalDebugLevel); // update status info
		  int dtt_len = clt_kernels[chn][0][0][0][0].length;
		  int dtt_size= (int)Math.sqrt(dtt_len);
		  String [] titles = {"red", "blue", "green"};
		  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
		  sdfa_instance.showArrays(
				  flat_kernels,
				  clt_kernels[chn][0][0].length*(2*dtt_size),
				  clt_kernels[chn][0].length*(2*dtt_size+1),
				  true,
				  "clt_kernels-"+chn,
				  titles);
	  }


	  public double [][][][][] extractCLTKernels (      // per color, save 4 kernelas and displacement as (2*dtt_size+1)*(2*dtt_size) tiles in an image (last row - shift x,y)
			  final float [][]   flat_kernels,    // per color/per tileY/ per tileX/per quadrant (plus offset as 5-th)/per pixel
			  final int          width,
			  final int          dtt_size,
			  final int          threadsMax,      // maximal number of threads to launch
			  final boolean      updateStatus,
			  final int          globalDebugLevel) // update status info
	  {
		  if (flat_kernels==null) return null;
		  final int nChn =       flat_kernels.length;
		  final int height =     flat_kernels[0].length/width;
		  final int tileWidth =  2 * dtt_size;
		  final int tileHeight = 2 * dtt_size + 1; // last row - shift with 0.5 pix steps
		  final int kernelNumHor =  width / tileWidth;
		  final int kernelNumVert = height / tileHeight;
		  final int dtt_len = dtt_size*dtt_size;
		  final double [][][][][] clt_kernels = new double [nChn][kernelNumVert][kernelNumHor][5][];
		  for (int chn = 0; chn < nChn; chn++){
			  for (int tileY = 0; tileY < kernelNumVert ; tileY++){
				  for (int tileX = 0; tileX < kernelNumHor ; tileX++){
					  for (int n = 0; n<4; n++){
						  clt_kernels[chn][tileY][tileX][n] = new double [dtt_len];
					  }
					  clt_kernels[chn][tileY][tileX][4] = new double [extra_items];
				  }
			  }
		  }

		  // currently each 64x64 kernel corresponds to 16x16 original pixels tile, 2 tiles margin each side
		  final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		  final AtomicInteger ai = new AtomicInteger(0);
		  final int numberOfKernels=     kernelNumHor*kernelNumVert*nChn;
		  final int numberOfKernelsInChn=kernelNumHor*kernelNumVert;

		  final long startTime = System.nanoTime();
		  System.out.println("flattenCLTKernels():numberOfKernels="+numberOfKernels);
		  for (int ithread = 0; ithread < threads.length; ithread++) {
			  threads[ithread] = new Thread() {
				  @Override
				public void run() {
					  int chn,tileY,tileX;
					  for (int nTile = ai.getAndIncrement(); nTile < numberOfKernels; nTile = ai.getAndIncrement()) {
						  chn=nTile/numberOfKernelsInChn;
						  tileY =(nTile % numberOfKernelsInChn)/kernelNumHor;
						  tileX = nTile % kernelNumHor;
						  for (int i = 0; i < dtt_size; i++){
							  for (int j = 0; j<dtt_size; j++){
								  int indx = i*dtt_size+j;
								  int baddr = (tileY*tileHeight + i) * width + (tileX * tileWidth) + j;
								  clt_kernels[chn][tileY][tileX][0][indx] = flat_kernels[chn][baddr];
								  clt_kernels[chn][tileY][tileX][1][indx] = flat_kernels[chn][baddr + dtt_size];
								  clt_kernels[chn][tileY][tileX][2][indx] = flat_kernels[chn][baddr + dtt_size * width];
								  clt_kernels[chn][tileY][tileX][3][indx] = flat_kernels[chn][baddr + dtt_size * width + dtt_size];
							  }
						  }
						  for (int i = 0; i < extra_items; i ++) {
							  clt_kernels[chn][tileY][tileX][4][i] = flat_kernels[chn][(tileY*tileHeight + 2 * dtt_size) * width + (tileX * tileWidth) + i];
						  }
					  }
				  }
			  };
		  }
		  ImageDtt.startAndJoin(threads);
		  if (globalDebugLevel > 1) System.out.println("Threads done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		  System.out.println("1.Threads done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		  /* prepare result stack to return */
		  return clt_kernels;
	  }


	  public boolean createCLTKernels(
			  EyesisCorrectionParameters.CLTParameters clt_parameters,
			  int          srcKernelSize,
			  int          threadsMax,  // maximal number of threads to launch
			  boolean      updateStatus,
			  int          debugLevel
			  ){
		  // get sensor geometry correction to apply to kernels as extra shifts
		  PixelMapping.SensorData [] sensors =  eyesisCorrections.pixelMapping.sensors;

		  String [] sharpKernelPaths= correctionsParameters.selectKernelChannelFiles(
				  0,  // 0 - sharp, 1 - smooth
				  eyesisCorrections.usedChannels.length, // numChannels, // number of channels
				  eyesisCorrections.debugLevel);
		  if (sharpKernelPaths==null) return false;
		  for (int i=0;i<sharpKernelPaths.length;i++){
			  System.out.println(i+":"+sharpKernelPaths[i]);
		  }
			if (clt_kernels == null){
				clt_kernels = new double[eyesisCorrections.usedChannels.length][][][][][];
				for (int chn=0;chn<eyesisCorrections.usedChannels.length;chn++){
					clt_kernels[chn] = null;
				}
			}
		  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();

		  for (int chn=0;chn<eyesisCorrections.usedChannels.length;chn++){
			  if (eyesisCorrections.usedChannels[chn] && (sharpKernelPaths[chn]!=null) && (clt_kernels[chn]==null)){
				  ImagePlus imp_kernel_sharp=new ImagePlus(sharpKernelPaths[chn]);
				  if (imp_kernel_sharp.getStackSize()<3) {
					  System.out.println("Need a 3-layer stack with kernels");
					  sharpKernelPaths[chn]=null;
					  continue;
				  }
				  ImageStack kernel_sharp_stack= imp_kernel_sharp.getStack();
				  System.out.println("debugLevel = "+debugLevel+" kernel_sharp_stack.getWidth() = "+kernel_sharp_stack.getWidth()+
						  " kernel_sharp_stack.getHeight() = "+kernel_sharp_stack.getHeight());

				  double [][][][][] kernels = calculateCLTKernel ( // per color/per tileY/ per tileX/per quadrant (plus offset as 5-th)/per pixel
						  sensors[chn],                            // to calculate extra shift (kernels are centered around green)
						  kernel_sharp_stack,                      // final ImageStack kernelStack,  // first stack with 3 colors/slices convolution kernels
						  srcKernelSize,                           // final int          kernelSize, // 64
						  clt_parameters,                          // final EyesisCorrectionParameters.CLTParameters clt_parameters,
						  threadsMax,  // maximal number of threads to launch
						  updateStatus,
						  debugLevel); // update status info

				  double [][] flat_kernels = flattenCLTKernels (      // per color, save 4 kernels and displacement as (2*dtt_size+1)*(2*dtt_size) tiles in an image (last row - shift x,y)
						  kernels,    // per color/per tileY/ per tileX/per quadrant (plus offset as 5-th)/per pixel
						  threadsMax,  // maximal number of threads to launch
						  updateStatus,
						  debugLevel); // update status info
				  int kernelNumHor=kernels[0][0].length;
				  int dtt_len = kernels[0][0][0][0].length;
				  int dtt_size =      (int) Math.sqrt(dtt_len);
				  int tileWidth =  2 * dtt_size;
				  int width =  tileWidth *  kernelNumHor;
				  int height = flat_kernels[0].length/width;

				  String [] layerNames = {"red_clt_kernels","blue_clt_kernels","green_clt_kernels"};
				  ImageStack cltStack = sdfa_instance.makeStack(
						  flat_kernels,
						  width,
						  height,
						  layerNames);

				  String cltPath=correctionsParameters.cltKernelDirectory+
						  Prefs.getFileSeparator()+
						  correctionsParameters.cltKernelPrefix+
						  String.format("%02d",chn)+
						  correctionsParameters.cltSuffix;
				  String msg="Saving CLT convolution kernels to "+cltPath;
				  IJ.showStatus(msg);
				  if (debugLevel>0) System.out.println(msg);
				  ImagePlus imp_clt=new ImagePlus(imp_kernel_sharp.getTitle()+"-clt",cltStack);
				  if (debugLevel > 0) {
					  imp_clt.getProcessor().resetMinAndMax();
					  imp_clt.show();
				  }
				  FileSaver fs=new FileSaver(imp_clt);
				  fs.saveAsTiffStack(cltPath); // directory does not exist
			  }
		  }
		  return true;
	  }



	  public boolean readCLTKernels(
			  EyesisCorrectionParameters.CLTParameters clt_parameters,
			  int          threadsMax,  // maximal number of threads to launch
			  boolean      updateStatus,
			  int          debugLevel
			  ){
		  int dtt_size = clt_parameters.transform_size;
		  String [] cltKernelPaths = correctionsParameters.selectCLTChannelFiles(
				  //					0,  // 0 - sharp, 1 - smooth
				  eyesisCorrections.usedChannels.length, // numChannels, // number of channels
				  eyesisCorrections.debugLevel);
		  if (cltKernelPaths==null) return false;

		  for (int i=0;i<cltKernelPaths.length;i++){
			  System.out.println(i+":"+cltKernelPaths[i]); // some may be null!
		  }

		  if (clt_kernels == null){
			  clt_kernels = new double[eyesisCorrections.usedChannels.length][][][][][];
			  for (int chn=0;chn<eyesisCorrections.usedChannels.length;chn++){
				  clt_kernels[chn] = null;
			  }
		  }
		  showDoubleFloatArrays sdfa_instance = null;
		  if (debugLevel>0){
			  sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		  }

		  for (int chn=0;chn<eyesisCorrections.usedChannels.length;chn++){
			  if (eyesisCorrections.usedChannels[chn] && (cltKernelPaths[chn]!=null)){
				  ImagePlus imp_kernel_clt=new ImagePlus(cltKernelPaths[chn]);
				  if (imp_kernel_clt.getStackSize()<3) {
					  System.out.println("Need a 3-layer stack with symmetrical DCT kernels");
					  cltKernelPaths[chn]=null;
					  continue;
				  }

				  ImageStack kernel_clt_stack=  imp_kernel_clt.getStack();
				  if (debugLevel>0){
					  System.out.println(" kernel_clt_stack.getWidth() = "+kernel_clt_stack.getWidth()+
							  " kernel_clt_stack.getHeight() = "+kernel_clt_stack.getHeight());
				  }
				  int nColors = kernel_clt_stack.getSize();
				  float [][] flat_kernels = new float [nColors][];
				  for (int nc = 0; nc < nColors; nc++){
					  flat_kernels[nc]= (float[]) kernel_clt_stack.getPixels(nc + 1);
				  }
				  clt_kernels[chn] =  extractCLTKernels ( // per color, save 4 kernelas and displacement as (2*dtt_size+1)*(2*dtt_size) tiles in an image (last row - shift x,y)
						  flat_kernels,                   // per color/per tileY/ per tileX/per quadrant (plus offset as 5-th)/per pixel
						  kernel_clt_stack.getWidth(),    // final int          width,
						  dtt_size,
						  threadsMax,                     // maximal number of threads to launch
						  updateStatus,
						  debugLevel);                    // update status info


				  if (sdfa_instance != null){
					  for (int nc = 0; nc < clt_kernels[chn].length; nc++){
					  double [][] dbg_clt = {
							  clt_kernels[chn][nc][clt_parameters.tileY/2][clt_parameters.tileX/2][0],
							  clt_kernels[chn][nc][clt_parameters.tileY/2][clt_parameters.tileX/2][1],
							  clt_kernels[chn][nc][clt_parameters.tileY/2][clt_parameters.tileX/2][2],
							  clt_kernels[chn][nc][clt_parameters.tileY/2][clt_parameters.tileX/2][3]};
					  String [] titles = {"CC", "SC", "CS", "SS"};
					  int length=dbg_clt[0].length;
					  int size=(int) Math.sqrt(length);
					  sdfa_instance.showArrays(
							  dbg_clt,
							  size,
							  size,
							  true,
							  "dbg_clt-"+nc,
							  titles);
					  System.out.println("readCLTKernels() chn="+chn+", color="+nc+" "+
							  "tileX = "+clt_parameters.tileX+" ("+(clt_parameters.tileX/2)+") "+
							  "tileY = "+clt_parameters.tileY+" ("+(clt_parameters.tileY/2)+") "+
							  "center_x = "+clt_kernels[chn][nc][clt_parameters.tileY/2][clt_parameters.tileX/2][4][0]+", "+
							  "center_y = "+clt_kernels[chn][nc][clt_parameters.tileY/2][clt_parameters.tileX/2][4][1]);
					  }
				  }
			  }
		  }
		  return true;
	  }


// mostly for testing
//eyesisKernelImage
	  public double [] extractOneKernelFromStack(
			  final int          kernelSize, // 64
			  final int chn,
			  final int xTile, // horizontal number of kernel to extract
			  final int yTile)  // vertical number of kernel to extract
	  {
		  if (eyesisKernelImage == null) return null;
		  final ImageStack kernelStack = eyesisKernelImage.getStack();
		  return extractOneKernelFromStack(
				  kernelStack,  // first stack with 3 colors/slices convolution kernels
				  kernelSize, // 64
				  chn,
				  xTile, // horizontal number of kernel to extract
				  yTile);  // vertical number of kernel to extract
	  }

	  public double [] extractOneKernelFromStack(
			  final ImageStack kernelStack,  // first stack with 3 colors/slices convolution kernels
			  final int          kernelSize, // 64
			  final int chn,
			  final int xTile, // horizontal number of kernel to extract
			  final int yTile)  // vertical number of kernel to extract
	  {

		  final int kernelWidth=kernelStack.getWidth();
		  final int kernelNumHor=kernelWidth/kernelSize;
		  double [] kernel = new double [kernelSize*kernelSize];
		  extractOneKernel((float[]) kernelStack.getPixels(chn+1), //  array of combined square kernels, each
				  kernel,        // will be filled, should have correct size before call
				  kernelNumHor,  // number of kernels in a row
				  xTile,         // horizontal number of kernel to extract
				  yTile);
		  return kernel;
	  }
	  // to be used in threaded method
	  private void extractOneKernel(float [] pixels, //  array of combined square kernels, each
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

	  public double [] reformatKernel(
			  double [] src_kernel,// will be blured in-place
			  int       src_size,  // typical 64
			  int       dst_size,  // typical 15 // destination size
			  int       decimation,// typical 2
			  double    sigma)
	  {
		  double [] dst_kernel = new double [dst_size*dst_size];
		  DoubleGaussianBlur gb = null;
		  if (sigma > 0) gb = new DoubleGaussianBlur();
		  reformatKernel(
				  src_kernel,
				  dst_kernel,
				  src_size,
				  dst_size,
				  decimation,
				  sigma,
				  gb);
		  return dst_kernel;

	  }
	  // to be used in threaded method
	  private void reformatKernel(
			  double [] src_kernel, // will be blured in-place
			  double [] dst_kernel,
			  int src_size,
			  int dst_size,
			  int decimation,
			  double sigma,
			  DoubleGaussianBlur gb)
	  {
		  if (gb != null) gb.blurDouble(src_kernel, src_size, src_size, sigma, sigma, 0.01);
		  int src_center = src_size / 2; // 32
		  int dst_center = dst_size / 2; // 7
		  for (int i = 0; i< dst_size; i++){
			  int src_i = (i - dst_center)*decimation + src_center;
			  if ((src_i >= 0) && (src_i < src_size)) {
				  for (int j = 0; j< dst_size; j++) {
					  int src_j = (j - dst_center)*decimation + src_center;
					  if ((src_j >= 0) && (src_j < src_size)) {
						  dst_kernel[i*dst_size + j] = src_kernel[src_i*src_size + src_j];
					  } else {
						  dst_kernel[i*dst_size + j] = 0;
					  }
				  }
			  } else {
				  for (int j = 0; j< dst_size; j++) dst_kernel[i*dst_size + j] = 0;
			  }
		  }
	  }
	  public double []reformatKernel2( // averages by exactly 2 (decimate==2)
			  double [] src_kernel, //
			  int src_size,
			  int dst_size){
		  double [] dst_kernel = new double [dst_size*dst_size];
		  reformatKernel2(
				  src_kernel,
				  dst_kernel,
				  src_size,
				  dst_size);
		  return dst_kernel;
	  }

	  private void reformatKernel2( // averages by exactly 2 (decimate==2)
			  double [] src_kernel, //
			  double [] dst_kernel,
			  int src_size,
			  int dst_size)
	  {
		  int decimation = 2;
		  int [] indices = {0,-src_size,-1,1,src_size,-src_size-1,-src_size+1,src_size-1,src_size+1};
		  double [] weights = {0.25,0.125,0.125,0.125,0.125,0.0625,0.0625,0.0625,0.0625};
		  int src_center = src_size / 2; // 32
		  int dst_center = dst_size / 2; // 7
		  int src_len = src_size*src_size;
		  for (int i = 0; i< dst_size; i++){
			  int src_i = (i - dst_center)*decimation + src_center;
			  if ((src_i >= 0) && (src_i < src_size)) {
				  for (int j = 0; j< dst_size; j++) {
					  int src_j = (j - dst_center)*decimation + src_center;
					  int dst_index = i*dst_size + j;
					  dst_kernel[dst_index] = 0.0;
					  if ((src_j >= 0) && (src_j < src_size)) {
						  int src_index = src_i*src_size + src_j;
						  for (int k = 0; k < indices.length; k++){
							  int indx = src_index + indices[k]; // normally source kernel should be larger, these lines just to save from "out of bounds"
							  if      (indx < 0)        indx += src_len;
							  else if (indx >= src_len) indx -= src_len;
							  dst_kernel[dst_index] += weights[k]*src_kernel[indx];
						  }
					  }
				  }
			  } else {
				  for (int j = 0; j< dst_size; j++) dst_kernel[i*dst_size + j] = 0;
			  }
		  }
	  }

	  public void resetCLTKernels() // and geometry corection too
	  {
		  clt_kernels = null;
		  geometryCorrection=null;

	  }

	  public ImageStack  YPrPbToRGB(double [][] yPrPb,
			  double Kr,        // 0.299;
			  double Kb,        // 0.114;
			  int width
			  ) {
		  int length = yPrPb[0].length;
		  int height = length/width;

		  float [] fpixels_r= new float [length];
		  float [] fpixels_g= new float [length];
		  float [] fpixels_b= new float [length];
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
			  Pb=yPrPb[2][i];
			  Pr=yPrPb[1][i];
			  Y =yPrPb[0][i];
			  fpixels_r[i]=(float) (Y+ Pr*KPrR);
			  fpixels_b[i]=(float) (Y+ Pb*KPbB);
			  fpixels_g[i]=(float) (Y+ Pr*KPrG + Pb*KPbG);
		  }
		  ImageStack stack=new ImageStack(width,height);
		  stack.addSlice("red",   fpixels_r);
		  stack.addSlice("green", fpixels_g);
		  stack.addSlice("blue",  fpixels_b);
		  return stack;

	  }

	  public double [][]  YPrPbToRBG(double [][] yPrPb,
			  double Kr,        // 0.299;
			  double Kb,        // 0.114;
			  int width
			  ) {
		  int length = yPrPb[0].length;
		  //			int height = length/width;

		  double [][]rbg = new double[3][length];
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
			  Pb=yPrPb[2][i];
			  Pr=yPrPb[1][i];
			  Y =yPrPb[0][i];
			  rbg[0][i]=(float) (Y+ Pr*KPrR);
			  rbg[1][i]=(float) (Y+ Pb*KPbB);
			  rbg[2][i]=(float) (Y+ Pr*KPrG + Pb*KPbG);
		  }
		  return rbg;
	  }

	  public void debayer_rbg(
			  ImageStack stack_rbg){
		  debayer_rbg(stack_rbg, 1.0);
	  }

	  // Simple in-place debayer by (bi) linear approximation, assumes [0R/00], [00/B0], [G0/0G] slices
	  public void debayer_rbg(
			  ImageStack stack_rbg,
			  double scale)
	  {
		  int width = stack_rbg.getWidth();
		  int height = stack_rbg.getHeight();
		  float [] fpixels_r = (float[]) stack_rbg.getPixels(1);
		  float [] fpixels_b = (float[]) stack_rbg.getPixels(2);
		  float [] fpixels_g = (float[]) stack_rbg.getPixels(3);
		  int [][][] av_row = {
				  {{1,1},{-1,1},{-1,-1}},
				  {{1,1},{-1,1},{-1,-1}},
				  {{1,1},{-1,1},{-1,-1}}};
		  int [][][] av_col =  {
				  {{ width, width},{ width, width},{ width, width}},
				  {{-width, width},{-width, width},{-width, width}},
				  {{-width,-width},{-width,-width},{-width,-width}}};

		  int [][][] av_xcross =  {
				  {{ width+1, width+1, width+1, width+1}, { width-1, width+1, width-1, width+1}, { width-1, width-1, width-1, width-1}},
				  {{-width+1, width+1,-width+1, width+1}, {-width-1,-width+1, width-1, width+1}, {-width-1,-width-1, width-1, width-1}},
				  {{-width+1,-width+1,-width+1,-width+1}, {-width-1,-width+1,-width-1,-width+1}, {-width-1,-width-1,-width-1,-width-1}}};
		  int [][][] av_plus =  {
				  {{ width,         1,       1, width},   { width,        -1,       1, width  }, { width,        -1,      -1, width  }},
				  {{-width,         1,       1, width},   {-width,        -1,       1, width  }, {-width,        -1,      -1, width  }},
				  {{-width,         1,       1,-width},   {-width,        -1,       1,-width  }, {-width,        -1,      -1,-width  }}};
		  for (int y = 0; y < height; y++){
			  boolean odd_row = (y & 1) != 0;
			  int row_type = (y==0)? 0: ((y==(height-1))?2:1);
			  for (int x = 0; x < width; x++){
				  int indx = y*width+x;
				  boolean odd_col =   (x & 1) != 0;
				  int col_type = (x==0)? 0: ((x==(width-1))?2:1);
				  if (odd_row){
					  if (odd_col){ // GB site
						  fpixels_r[indx] = 0.5f*(
								  fpixels_r[indx+av_col[row_type][col_type][0]]+
								  fpixels_r[indx+av_col[row_type][col_type][1]]);
						  fpixels_b[indx] = 0.5f*(
								  fpixels_b[indx+av_row[row_type][col_type][0]]+
								  fpixels_b[indx+av_row[row_type][col_type][1]]);
					  } else { // !odd col  // B site
						  fpixels_r[indx] = 0.25f*(
								  fpixels_r[indx+av_xcross[row_type][col_type][0]]+
								  fpixels_r[indx+av_xcross[row_type][col_type][1]]+
								  fpixels_r[indx+av_xcross[row_type][col_type][2]]+
								  fpixels_r[indx+av_xcross[row_type][col_type][3]]);
						  fpixels_g[indx] = 0.25f*(
								  fpixels_g[indx+av_plus[row_type][col_type][0]]+
								  fpixels_g[indx+av_plus[row_type][col_type][1]]+
								  fpixels_g[indx+av_plus[row_type][col_type][2]]+
								  fpixels_g[indx+av_plus[row_type][col_type][3]]);
					  }

				  } else { // !odd_row
					  if (odd_col){  // R site
						  fpixels_b[indx] = 0.25f*(
								  fpixels_b[indx+av_xcross[row_type][col_type][0]]+
								  fpixels_b[indx+av_xcross[row_type][col_type][1]]+
								  fpixels_b[indx+av_xcross[row_type][col_type][2]]+
								  fpixels_b[indx+av_xcross[row_type][col_type][3]]);
						  fpixels_g[indx] = 0.25f*(
								  fpixels_g[indx+av_plus[row_type][col_type][0]]+
								  fpixels_g[indx+av_plus[row_type][col_type][1]]+
								  fpixels_g[indx+av_plus[row_type][col_type][2]]+
								  fpixels_g[indx+av_plus[row_type][col_type][3]]);
					  } else { // !odd col  // G site
						  fpixels_r[indx] = 0.5f*(
								  fpixels_r[indx+av_row[row_type][col_type][0]]+
								  fpixels_r[indx+av_row[row_type][col_type][1]]);
						  fpixels_b[indx] = 0.5f*(
								  fpixels_b[indx+av_col[row_type][col_type][0]]+
								  fpixels_b[indx+av_col[row_type][col_type][1]]);
					  }
				  }
			  }
		  }
		  if (scale !=1.0){
			  for (int i = 0; i< fpixels_r.length; i++){
				  fpixels_r[i] *= scale;
				  fpixels_b[i] *= scale;
				  fpixels_g[i] *= scale;
			  }
		  }
	  }

	  public void processCLTChannelImages(
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
//			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
//			  int          convolveFFTSize, // 128 - fft size, kernel size should be size/2
			  final int          threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel)
	  {
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
					  // removeUnusedSensorData should be off!?
					  channels=this.eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (eyesisCorrections.isChannelEnabled(channels[i])){
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
		  double [] referenceExposures=eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		  int [][] fileIndices=new int [numImagesToProcess][2]; // file index, channel number
		  int index=0;
		  for (int nFile=0;nFile<enabledFiles.length;nFile++){
			  if ((sourceFiles[nFile]!=null) && (sourceFiles[nFile].length()>1)) {
				  int [] channels={correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile])};
				  if (correctionsParameters.isJP4()){
					  int subCamera= channels[0]- correctionsParameters.firstSubCamera; // to match those in the sensor files
					  channels=eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (eyesisCorrections.isChannelEnabled(channels[i])){
						  fileIndices[index  ][0]=nFile;
						  fileIndices[index++][1]=channels[i];
					  }
				  }
			  }
		  }
		  for (int iImage=0;iImage<fileIndices.length;iImage++){
			  int nFile=fileIndices[iImage][0];
			  ImagePlus imp_src=null;
			  //				  int srcChannel=correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile]);
			  int srcChannel=fileIndices[iImage][1];
			  if (correctionsParameters.isJP4()){
				  int subchannel=eyesisCorrections.pixelMapping.getSubChannel(srcChannel);
				  if (this.correctionsParameters.swapSubchannels01) {
					  switch (subchannel){
					  case 0: subchannel=1; break;
					  case 1: subchannel=0; break;
					  }
				  }
				  if (debugLevel>0) System.out.println("Processing channel "+fileIndices[iImage][1]+" - subchannel "+subchannel+" of "+sourceFiles[nFile]);
				  ImagePlus imp_composite=eyesisCorrections.JP4_INSTANCE.open(
						  "", // path,
						  sourceFiles[nFile],
						  "",  //arg - not used in JP46 reader
						  true, // un-apply camera color gains
						  null, // new window
						  false); // do not show
				  imp_src=eyesisCorrections.JP4_INSTANCE.demuxImage(imp_composite, subchannel);
				  if (imp_src==null) imp_src=imp_composite; // not a composite image

				  // do we need to add any properties?
			  } else {
				  imp_src=new ImagePlus(sourceFiles[nFile]);
				  //					  (new JP46_Reader_camera(false)).decodeProperiesFromInfo(imp_src); // decode existent properties from info
				  eyesisCorrections.JP4_INSTANCE.decodeProperiesFromInfo(imp_src); // decode existent properties from info
				  if (debugLevel>0) System.out.println("Processing "+sourceFiles[nFile]);
			  }
			  double scaleExposure=1.0;
			  if (!Double.isNaN(referenceExposures[nFile]) && (imp_src.getProperty("EXPOSURE")!=null)){
				  scaleExposure=referenceExposures[nFile]/Double.parseDouble((String) imp_src.getProperty("EXPOSURE"));
				  //					  imp_src.setProperty("scaleExposure", scaleExposure); // it may already have channel
				  if (debugLevel>0) System.out.println("Will scale intensity (to compensate for exposure) by  "+scaleExposure);
			  }
			  imp_src.setProperty("name",    correctionsParameters.getNameFromSourceTiff(sourceFiles[nFile]));
			  imp_src.setProperty("channel", srcChannel); // it may already have channel
			  imp_src.setProperty("path",    sourceFiles[nFile]); // it may already have channel
			  //				  ImagePlus result=processChannelImage( // returns ImagePlus, but it already should be saved/shown
			  processCLTChannelImage( // returns ImagePlus, but it already should be saved/shown
					  imp_src, // should have properties "name"(base for saving results), "channel","path"
					  clt_parameters,
					  debayerParameters,
//					  nonlinParameters,
					  colorProcParameters,
					  channelGainParameters,
					  rgbParameters,
//					  convolveFFTSize, // 128 - fft size, kernel size should be size/2
					  scaleExposure,
					  threadsMax,  // maximal number of threads to launch
					  updateStatus,
					  debugLevel);
			  // warp result (add support for different color modes)
			  if (this.correctionsParameters.equirectangular){
				  if (equirectangularParameters.clearFullMap) eyesisCorrections.pixelMapping.deleteEquirectangularMapFull(srcChannel); // save memory? //removeUnusedSensorData - no, use equirectangular specific settings
				  if (equirectangularParameters.clearAllMaps) eyesisCorrections.pixelMapping.deleteEquirectangularMapAll(srcChannel); // save memory? //removeUnusedSensorData - no, use equirectangular specific settings
			  }
			  //pixelMapping
			  Runtime.getRuntime().gc();
			  if (debugLevel >-1) System.out.println("Processing image "+(iImage+1)+" (of "+fileIndices.length+") finished at "+
					  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
			  if (eyesisCorrections.stopRequested.get()>0) {
				  System.out.println("User requested stop");
				  return;
			  }
		  }
		  System.out.println("processCLTChannelImages(): Processing "+fileIndices.length+" files finished at "+
				  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
	  }

	  public ImagePlus processCLTChannelImage(
			  ImagePlus imp_src, // should have properties "name"(base for saving results), "channel","path"
//			  EyesisCorrectionParameters.DCTParameters           dct_parameters,
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
//			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
//			  int          convolveFFTSize, // 128 - fft size, kernel size should be size/2
			  double 		     scaleExposure,
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel){
		  boolean advanced=this.correctionsParameters.zcorrect || this.correctionsParameters.equirectangular;
		  boolean crop=      advanced? true: this.correctionsParameters.crop;
		  boolean rotate=    advanced? false: this.correctionsParameters.rotate;
		  double JPEG_scale= advanced? 1.0: this.correctionsParameters.JPEG_scale;
		  boolean toRGB=     advanced? true: this.correctionsParameters.toRGB;
		  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?

		  // may use this.StartTime to report intermediate steps execution times
		  String name=(String) imp_src.getProperty("name");
		  //		int channel= Integer.parseInt((String) imp_src.getProperty("channel"));
		  int channel= (Integer) imp_src.getProperty("channel");
		  String path= (String) imp_src.getProperty("path");
		  if (this.correctionsParameters.pixelDefects && (eyesisCorrections.defectsXY!=null)&& (eyesisCorrections.defectsXY[channel]!=null)){
			  // apply pixel correction
			  int numApplied=	eyesisCorrections.correctDefects(
					  imp_src,
					  channel,
					  debugLevel);
			  if ((debugLevel>0) && (numApplied>0)) { // reduce verbosity after verified defect correction works
				  System.out.println("Corrected "+numApplied+" pixels in "+path);
			  }
		  }

		  if (this.correctionsParameters.vignetting){
			  if ((eyesisCorrections.channelVignettingCorrection==null) || (channel<0) || (channel>=eyesisCorrections.channelVignettingCorrection.length) || (eyesisCorrections.channelVignettingCorrection[channel]==null)){
				  System.out.println("No vignetting data for channel "+channel);
				  return null;
			  }
			  float [] pixels=(float []) imp_src.getProcessor().getPixels();
			  if (pixels.length!=eyesisCorrections.channelVignettingCorrection[channel].length){
				  System.out.println("Vignetting data for channel "+channel+" has "+eyesisCorrections.channelVignettingCorrection[channel].length+" pixels, image "+path+" has "+pixels.length);
				  return null;
			  }
			  // TODO: Move to do it once:
			  double min_non_zero = 0.0;
			  for (int i=0;i<pixels.length;i++){
				  double d = eyesisCorrections.channelVignettingCorrection[channel][i];
				  if ((d > 0.0) && ((min_non_zero == 0) || (min_non_zero > d))){
					  min_non_zero = d;
				  }
			  }
			  double max_vign_corr = clt_parameters.vignetting_range*min_non_zero;

			  System.out.println("Vignetting data: channel="+channel+", min = "+min_non_zero);
			  for (int i=0;i<pixels.length;i++){
				  double d = eyesisCorrections.channelVignettingCorrection[channel][i];
				  if (d > max_vign_corr) d = max_vign_corr;
				  pixels[i]*=d;
			  }
			  // Scale here, combine with vignetting later?
			  int width =  imp_src.getWidth();
			  int height = imp_src.getHeight();
			  for (int y = 0; y < height-1; y+=2){
				  for (int x = 0; x < width-1; x+=2){
					  pixels[y*width+x        ] *= clt_parameters.scale_g;
					  pixels[y*width+x+width+1] *= clt_parameters.scale_g;
					  pixels[y*width+x      +1] *= clt_parameters.scale_r;
					  pixels[y*width+x+width  ] *= clt_parameters.scale_b;
				  }
			  }

		  } else { // assuming GR/BG pattern
			  System.out.println("Applying fixed color gain correction parameters: Gr="+
					  clt_parameters.novignetting_r+", Gg="+clt_parameters.novignetting_g+", Gb="+clt_parameters.novignetting_b);
			  float [] pixels=(float []) imp_src.getProcessor().getPixels();
			  int width =  imp_src.getWidth();
			  int height = imp_src.getHeight();
			  double kr = clt_parameters.scale_r/clt_parameters.novignetting_r;
			  double kg = clt_parameters.scale_g/clt_parameters.novignetting_g;
			  double kb = clt_parameters.scale_b/clt_parameters.novignetting_b;
			  for (int y = 0; y < height-1; y+=2){
				  for (int x = 0; x < width-1; x+=2){
					  pixels[y*width+x        ] *= kg;
					  pixels[y*width+x+width+1] *= kg;
					  pixels[y*width+x      +1] *= kr;
					  pixels[y*width+x+width  ] *= kb;
				  }
			  }
		  }
		  if (clt_parameters.gain_equalize){

		  }

		  String title=name+"-"+String.format("%02d", channel);
		  ImagePlus result=imp_src;
		  if (debugLevel>1) System.out.println("processing: "+path);
		  result.setTitle(title+"RAW");
		  if (!this.correctionsParameters.split){
			  eyesisCorrections.saveAndShow(result, this.correctionsParameters);
			  return result;
		  }
		  // Generate split parameters for DCT processing mode
		  EyesisCorrectionParameters.SplitParameters splitParameters = new EyesisCorrectionParameters.SplitParameters(
				                 1,  // oversample; // currently source kernels are oversampled
						  clt_parameters.transform_size/2, // addLeft
						  clt_parameters.transform_size/2, // addTop
						  clt_parameters.transform_size/2, // addRight
						  clt_parameters.transform_size/2  // addBottom
				  );

		  // Split into Bayer components, oversample, increase canvas

		  double [][] double_stack = eyesisCorrections.bayerToDoubleStack(
				  result, // source Bayer image, linearized, 32-bit (float))
				  null); // no margins, no oversample

//		  ImageStack stack= eyesisCorrections.bayerToStack(
//				  result, // source Bayer image, linearized, 32-bit (float))
//				  splitParameters);
		  String titleFull=title+"-SPLIT";
		  if (debugLevel > -1){
			  double [] chn_avg = {0.0,0.0,0.0};
			  int width =  imp_src.getWidth();
			  int height = imp_src.getHeight();
			  for (int c = 0; c < 3; c++){
				  for (int i = 0; i<double_stack[c].length; i++){
					  chn_avg[c] += double_stack[c][i];
				  }
			  }
			  chn_avg[0] /= width*height/4;
			  chn_avg[1] /= width*height/4;
			  chn_avg[2] /= width*height/2;
			  System.out.println("Split channels averages: R="+chn_avg[0]+", G="+chn_avg[2]+", B="+chn_avg[1]);
		  }
		  String [] rbg_titles = {"Red", "Blue", "Green"};
		  ImageStack stack;
		  if (!this.correctionsParameters.debayer) {
//			  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
//			  ImageStack
			  stack = sdfa_instance.makeStack(double_stack, imp_src.getWidth(), imp_src.getHeight(), rbg_titles);
			  result= new ImagePlus(titleFull, stack);
			  eyesisCorrections.saveAndShow(result, this.correctionsParameters);
			  return result;
		  }
		  // =================
		  if (debugLevel > 0) {
			  System.out.println("Showing image BEFORE_CLT_PROC");
			  sdfa_instance.showArrays(double_stack,  imp_src.getWidth(), imp_src.getHeight(), true, "BEFORE_CLT_PROC", rbg_titles);
		  }
		  if (this.correctionsParameters.deconvolve) { // process with DCT, otherwise use simple debayer
			  ImageDtt image_dtt = new ImageDtt();
/*
			  double [][][][][] clt_data = image_dtt.cltStack(
					  stack,
					  channel,
					  clt_parameters,
					  clt_parameters.ishift_x, //final int shiftX, // shift image horizontally (positive - right)
					  clt_parameters.ishift_y, //final int shiftY, // shift image vertically (positive - down)
					  threadsMax,
					  debugLevel,
					  updateStatus);
*/
			  for (int i =0 ; i < double_stack[0].length; i++){
//				  double_stack[0][i]*=2.0; // Scale red twice to compensate less pixels than green
//				  double_stack[1][i]*=2.0; // Scale blue twice to compensate less pixels than green
				  double_stack[2][i]*=0.5; // Scale blue twice to compensate less pixels than green
			  }
			  double [][][][][] clt_data = image_dtt.clt_aberrations(
					  double_stack,                 // final double [][]       imade_data,
					  imp_src.getWidth(),           //	final int               width,
					  clt_kernels[channel],         // final double [][][][][] clt_kernels, // [color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
					  clt_parameters.kernel_step,
					  clt_parameters.transform_size,
					  clt_parameters.clt_window,
					  clt_parameters.shift_x,       // final int               shiftX, // shift image horizontally (positive - right) - just for testing
					  clt_parameters.shift_y,       // final int               shiftY, // shift image vertically (positive - down)
					  clt_parameters.tileX,         // final int               debug_tileX,
					  clt_parameters.tileY,         // final int               debug_tileY,
					  (clt_parameters.dbg_mode & 64) != 0, // no fract shift
					  (clt_parameters.dbg_mode & 128) != 0, // no convolve
					  (clt_parameters.dbg_mode & 256) != 0, // transpose convolve
					  threadsMax,
					  debugLevel);
//					  updateStatus);


			  System.out.println("clt_data.length="+clt_data.length+" clt_data[0].length="+clt_data[0].length
					  +" clt_data[0][0].length="+clt_data[0][0].length+" clt_data[0][0][0].length="+
					  clt_data[0][0][0].length+" clt_data[0][0][0][0].length="+clt_data[0][0][0][0].length);
			  /*
			  if (dct_parameters.color_DCT){ // convert RBG -> YPrPb
				  dct_data = image_dtt.dct_color_convert(
						  dct_data,
						  colorProcParameters.kr,
						  colorProcParameters.kb,
						  dct_parameters.sigma_rb,        // blur of channels 0,1 (r,b) in addition to 2 (g)
						  dct_parameters.sigma_y,         // blur of Y from G
						  dct_parameters.sigma_color,     // blur of Pr, Pb in addition to Y
						  threadsMax,
						  debugLevel);
			  } else { // just LPF RGB
			  */
				  if (clt_parameters.corr_sigma > 0){ // no filter at all
					  for (int chn = 0; chn < clt_data.length; chn++) {
						  image_dtt.clt_lpf(
								  clt_parameters.corr_sigma,
								  clt_data[chn],
								  clt_parameters.transform_size,
								  threadsMax,
								  debugLevel);
					  }
				  }
/*
			  }
*/
			  int tilesY = imp_src.getHeight()/clt_parameters.transform_size;
			  int tilesX = imp_src.getWidth()/clt_parameters.transform_size;
			  if (debugLevel > 0){
				  System.out.println("--tp.tilesX="+tilesX);
				  System.out.println("--tp.tilesY="+tilesY);
			  }
			  if (debugLevel > 1){
			        double [][] clt = new double [clt_data.length*4][];
			        for (int chn = 0; chn < clt_data.length; chn++) {
			        	double [][] clt_set = image_dtt.clt_dbg(
			        			clt_data [chn],
								  threadsMax,
								  debugLevel);
			        	for (int ii = 0; ii < clt_set.length; ii++) clt[chn*4+ii] = clt_set[ii];
			        }

			        if (debugLevel > 0){
			        	sdfa_instance.showArrays(clt,
			        			tilesX*clt_parameters.transform_size,
			        			tilesY*clt_parameters.transform_size,
			        			true,
			        			result.getTitle()+"-CLT");
			        }
			  }
			  double [][] iclt_data = new double [clt_data.length][];
			  for (int chn=0; chn<clt_data.length;chn++){
				  iclt_data[chn] = image_dtt.iclt_2d(
						  clt_data[chn],                  // scanline representation of dcd data, organized as dct_size x dct_size tiles
						  clt_parameters.transform_size,  // final int
						  clt_parameters.clt_window,      // window_type
						  15,                             // clt_parameters.iclt_mask,       //which of 4 to transform back
						  0,                              // clt_parameters.dbg_mode,        //which of 4 to transform back
						  threadsMax,
						  debugLevel);

			  }
					  if (debugLevel > -1) sdfa_instance.showArrays(
							  iclt_data,
							  (tilesX + 1) * clt_parameters.transform_size,
							  (tilesY + 1) * clt_parameters.transform_size,
							  true,
							  result.getTitle()+"-rbg_sigma");
				/*
				  }
			  }
			 */
			  if (debugLevel > 0) sdfa_instance.showArrays(iclt_data,
					  (tilesX + 1) * clt_parameters.transform_size,
					  (tilesY + 1) * clt_parameters.transform_size,
					  true,
					  result.getTitle()+"-ICLT-RGB");

			  // convert to ImageStack of 3 slices
			  String [] sliceNames = {"red", "blue", "green"};
			  stack = sdfa_instance.makeStack(
					  iclt_data,
					  (tilesX + 1) * clt_parameters.transform_size,
					  (tilesY + 1) * clt_parameters.transform_size,
					  sliceNames); // or use null to get chn-nn slice names


		  } else { // if (this.correctionsParameters.deconvolve) - here use a simple debayer
			  System.out.println("Bypassing CLT-based aberration correction");
			  stack = sdfa_instance.makeStack(double_stack, imp_src.getWidth(), imp_src.getHeight(), rbg_titles);
			  debayer_rbg(stack, 0.25); // simple standard 3x3 kernel debayer
		  }
		  if (debugLevel > -1){
			  double [] chn_avg = {0.0,0.0,0.0};
			  float [] pixels;
			  int width =  stack.getWidth();
			  int height = stack.getHeight();

			  for (int c = 0; c <3; c++){
				  pixels = (float[]) stack.getPixels(c+1);
				  for (int i = 0; i<pixels.length; i++){
					  chn_avg[c] += pixels[i];
				  }
			  }
			  chn_avg[0] /= width*height;
			  chn_avg[1] /= width*height;
			  chn_avg[2] /= width*height;
			  System.out.println("Processed channels averages: R="+chn_avg[0]+", G="+chn_avg[2]+", B="+chn_avg[1]);
		  }

		  if (!this.correctionsParameters.colorProc){
			  result= new ImagePlus(titleFull, stack);
			  eyesisCorrections.saveAndShow(
					  result,
					  this.correctionsParameters);
			  return result;
		  }
		  if (debugLevel > 1) System.out.println("before colors.1");
		  //Processing colors - changing stack sequence to r-g-b (was r-b-g)
		  if (!eyesisCorrections.fixSliceSequence(
				  stack,
				  debugLevel)){
			  if (debugLevel > -1) System.out.println("fixSliceSequence() returned false");
			  return null;
		  }
		  if (debugLevel > 1) System.out.println("before colors.2");
		  if (debugLevel > 1){
			  ImagePlus imp_dbg=new ImagePlus(imp_src.getTitle()+"-"+channel+"-preColors",stack);
			  eyesisCorrections.saveAndShow(
					  imp_dbg,
					  this.correctionsParameters);
		  }
		  if (debugLevel > 1) System.out.println("before colors.3, scaleExposure="+scaleExposure+" scale = "+(255.0/eyesisCorrections.psfSubpixelShouldBe4/eyesisCorrections.psfSubpixelShouldBe4/scaleExposure));
		  CorrectionColorProc correctionColorProc=new CorrectionColorProc(eyesisCorrections.stopRequested);
		  double [][] yPrPb=new double [3][];
		  //			if (dct_parameters.color_DCT){
		  // need to get YPbPr - not RGB here
		  //			} else {
		  correctionColorProc.processColorsWeights(stack, // just gamma convert? TODO: Cleanup? Convert directly form the linear YPrPb
				  //					  255.0/this.psfSubpixelShouldBe4/this.psfSubpixelShouldBe4, //  double scale,     // initial maximal pixel value (16))
				  //					  255.0/eyesisCorrections.psfSubpixelShouldBe4/eyesisCorrections.psfSubpixelShouldBe4/scaleExposure, //  double scale,     // initial maximal pixel value (16))
				  //					  255.0/2/2/scaleExposure, //  double scale,     // initial maximal pixel value (16))
				  255.0/scaleExposure, //  double scale,     // initial maximal pixel value (16))
				  colorProcParameters,
				  channelGainParameters,
				  channel,
				  null, //correctionDenoise.getDenoiseMask(),
				  this.correctionsParameters.blueProc,
				  debugLevel);
		  if (debugLevel > 1) System.out.println("Processed colors to YPbPr, total number of slices="+stack.getSize());
		  if (debugLevel > 1) {
			  ImagePlus imp_dbg=new ImagePlus("procColors",stack);
			  eyesisCorrections.saveAndShow(
					  imp_dbg,
					  this.correctionsParameters);
		  }
		  float [] fpixels;
		  int [] slices_YPrPb = {8,6,7};
		  yPrPb=new double [3][];
		  for (int n = 0; n < slices_YPrPb.length; n++){
			  fpixels = (float[]) stack.getPixels(slices_YPrPb[n]);
			  yPrPb[n] = new double [fpixels.length];
			  for (int i = 0; i < fpixels.length; i++) yPrPb[n][i] = fpixels[i];
		  }

		  if (toRGB) {
			  if (debugLevel > 0){
				  System.out.println("correctionColorProc.YPrPbToRGB");
			  }
			  stack =  YPrPbToRGB(yPrPb,
					  colorProcParameters.kr,        // 0.299;
					  colorProcParameters.kb,        // 0.114;
					  stack.getWidth());

			  title=titleFull; // including "-DECONV" or "-COMBO"
			  titleFull=title+"-RGB-float";
			  //Trim stack to just first 3 slices
			  if (debugLevel > 1){ // 2){
				  ImagePlus imp_dbg=new ImagePlus("YPrPbToRGB",stack);
				  eyesisCorrections.saveAndShow(
						  imp_dbg,
						  this.correctionsParameters);
			  }
			  while (stack.getSize() > 3) stack.deleteLastSlice();
			  if (debugLevel > 1) System.out.println("Trimming color stack");
		  } else {
			  title=titleFull; // including "-DECONV" or "-COMBO"
			  titleFull=title+"-YPrPb"; // including "-DECONV" or "-COMBO"
			  if (debugLevel > 1) System.out.println("Using full stack, including YPbPr");
		  }

		  result= new ImagePlus(titleFull, stack);
		  // Crop image to match original one (scaled to oversampling)
		  if (crop){ // always crop if equirectangular
			  if (debugLevel > 1) System.out.println("cropping");
			  stack = eyesisCorrections.cropStack32(stack,splitParameters);
			  if (debugLevel > 2) { // 2){
				  ImagePlus imp_dbg=new ImagePlus("cropped",stack);
				  eyesisCorrections.saveAndShow(
						  imp_dbg,
						  this.correctionsParameters);
			  }
		  }
		  // rotate the result
		  if (rotate){ // never rotate for equirectangular
			  stack=eyesisCorrections.rotateStack32CW(stack);
		  }
		  if (!toRGB && !this.correctionsParameters.jpeg){ // toRGB set for equirectangular
			  if (debugLevel > 1) System.out.println("!toRGB && !this.correctionsParameters.jpeg");
			  eyesisCorrections.saveAndShow(result, this.correctionsParameters);
			  return result;
		  } else { // that's not the end result, save if required
			  if (debugLevel > 1) System.out.println("!toRGB && !this.correctionsParameters.jpeg - else");
			  eyesisCorrections.saveAndShow(result,
					  eyesisCorrections.correctionsParameters,
					  eyesisCorrections.correctionsParameters.save32,
					  false,
					  eyesisCorrections.correctionsParameters.JPEG_quality); // save, no show
		  }
		  // convert to RGB48 (16 bits per color component)
		  ImagePlus imp_RGB;
		  stack=eyesisCorrections.convertRGB32toRGB16Stack(
				  stack,
				  rgbParameters);

		  titleFull=title+"-RGB48";
		  result= new ImagePlus(titleFull, stack);
		  //			  ImagePlus imp_RGB24;
		  result.updateAndDraw();
		  if (debugLevel > 1) System.out.println("result.updateAndDraw(), "+titleFull+"-RGB48");

		  CompositeImage compositeImage=eyesisCorrections.convertToComposite(result);

		  if (!this.correctionsParameters.jpeg && !advanced){ // RGB48 was the end result
			  if (debugLevel > 1) System.out.println("if (!this.correctionsParameters.jpeg && !advanced)");
			  eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters);
			  return result;
		  } else { // that's not the end result, save if required
			  if (debugLevel > 1) System.out.println("if (!this.correctionsParameters.jpeg && !advanced) - else");
			  eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters, this.correctionsParameters.save16, false); // save, no show
			  //				  eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters, this.correctionsParameters.save16, true); // save, no show
		  }

		  imp_RGB=eyesisCorrections.convertRGB48toRGB24(
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
		  eyesisCorrections.saveAndShow(imp_RGB, this.correctionsParameters);

		  return result;
	  }

// Processing sets of 4 images together
	  public void processCLTSets(
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
//			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
//			  int          convolveFFTSize, // 128 - fft size, kernel size should be size/2
			  final int          threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel)
	  {
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
					  // removeUnusedSensorData should be off!?
					  channels=this.eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (eyesisCorrections.isChannelEnabled(channels[i])){
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
		  double [] referenceExposures=eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		  int [][] fileIndices=new int [numImagesToProcess][2]; // file index, channel number
		  int index=0;
		  for (int nFile=0;nFile<enabledFiles.length;nFile++){
			  if ((sourceFiles[nFile]!=null) && (sourceFiles[nFile].length()>1)) {
				  int [] channels={correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile])};
				  if (correctionsParameters.isJP4()){
					  int subCamera= channels[0]- correctionsParameters.firstSubCamera; // to match those in the sensor files
					  channels=eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (eyesisCorrections.isChannelEnabled(channels[i])){
						  fileIndices[index  ][0]=nFile;
						  fileIndices[index++][1]=channels[i];
					  }
				  }
			  }
		  }
		  ArrayList<String> setNames = new ArrayList<String>();
		  ArrayList<ArrayList<Integer>> setFiles = new ArrayList<ArrayList<Integer>>();

		  for (int iImage=0;iImage<fileIndices.length;iImage++){
			  int nFile=fileIndices[iImage][0];
			  String setName = correctionsParameters.getNameFromSourceTiff(sourceFiles[nFile]);
			  if (!setNames.contains(setName)) {
				  setNames.add(setName);
				  setFiles.add(new ArrayList<Integer>());
			  }
			  setFiles.get(setNames.indexOf(setName)).add(new Integer(nFile));
		  }
		  int iImage = 0;
		  for (int nSet = 0; nSet < setNames.size(); nSet++){
			  int maxChn = 0;
			  for (int i = 0; i < setFiles.get(nSet).size(); i++){
				  int chn = fileIndices[setFiles.get(nSet).get(i)][1];
				  if (chn > maxChn) maxChn = chn;
			  }
			  int [] channelFiles = new int[maxChn+1];
			  for (int i =0; i < channelFiles.length; i++) channelFiles[i] = -1;
			  for (int i = 0; i < setFiles.get(nSet).size(); i++){
				  channelFiles[fileIndices[setFiles.get(nSet).get(i)][1]] = setFiles.get(nSet).get(i);
			  }

			  ImagePlus [] imp_srcs = new ImagePlus[channelFiles.length];
			  double [] scaleExposure = new double[channelFiles.length];
			  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
				  int nFile=channelFiles[srcChannel];
				  imp_srcs[srcChannel]=null;
				  if (nFile >=0){
					  if (correctionsParameters.isJP4()){
						  int subchannel=eyesisCorrections.pixelMapping.getSubChannel(srcChannel);
						  if (this.correctionsParameters.swapSubchannels01) {
							  switch (subchannel){
							  case 0: subchannel=1; break;
							  case 1: subchannel=0; break;
							  }
						  }
						  if (debugLevel>0) System.out.println("Processing set " + setNames.get(nSet)+" channel "+srcChannel+" - subchannel "+subchannel+" of "+sourceFiles[nFile]);
						  ImagePlus imp_composite=eyesisCorrections.JP4_INSTANCE.open(
								  "", // path,
								  sourceFiles[nFile],
								  "",  //arg - not used in JP46 reader
								  true, // un-apply camera color gains
								  null, // new window
								  false); // do not show
						  imp_srcs[srcChannel]=eyesisCorrections.JP4_INSTANCE.demuxImage(imp_composite, subchannel);
						  if (imp_srcs[srcChannel] == null) imp_srcs[srcChannel] = imp_composite; // not a composite image
						  // do we need to add any properties?
					  } else {
						  imp_srcs[srcChannel]=new ImagePlus(sourceFiles[nFile]);
						  //					  (new JP46_Reader_camera(false)).decodeProperiesFromInfo(imp_src); // decode existent properties from info
						  eyesisCorrections.JP4_INSTANCE.decodeProperiesFromInfo(imp_srcs[srcChannel]); // decode existent properties from info
						  if (debugLevel>0) System.out.println("Processing "+sourceFiles[nFile]);
					  }
					  scaleExposure[srcChannel] = 1.0;
					  if (!Double.isNaN(referenceExposures[nFile]) && (imp_srcs[srcChannel].getProperty("EXPOSURE")!=null)){
						  scaleExposure[srcChannel] = referenceExposures[nFile]/Double.parseDouble((String) imp_srcs[srcChannel].getProperty("EXPOSURE"));
						  if (debugLevel>0) System.out.println("Will scale intensity (to compensate for exposure) by  "+scaleExposure);
					  }
					  imp_srcs[srcChannel].setProperty("name",    correctionsParameters.getNameFromSourceTiff(sourceFiles[nFile]));
					  imp_srcs[srcChannel].setProperty("channel", srcChannel); // it may already have channel
					  imp_srcs[srcChannel].setProperty("path",    sourceFiles[nFile]); // it may already have channel

					  if (this.correctionsParameters.pixelDefects && (eyesisCorrections.defectsXY!=null)&& (eyesisCorrections.defectsXY[srcChannel]!=null)){
						  // apply pixel correction
						  int numApplied=	eyesisCorrections.correctDefects(
								  imp_srcs[srcChannel],
								  srcChannel,
								  debugLevel);
						  if ((debugLevel>0) && (numApplied>0)) { // reduce verbosity after verified defect correction works
							  System.out.println("Corrected "+numApplied+" pixels in "+sourceFiles[nFile]);
						  }
					  }

					  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
					  int width =  imp_srcs[srcChannel].getWidth();
					  int height = imp_srcs[srcChannel].getHeight();

					  if (clt_parameters.sat_level > 0.0){
						  double [] saturations = {
								  Double.parseDouble((String) imp_srcs[srcChannel].getProperty("saturation_1")),
								  Double.parseDouble((String) imp_srcs[srcChannel].getProperty("saturation_0")),
								  Double.parseDouble((String) imp_srcs[srcChannel].getProperty("saturation_3")),
								  Double.parseDouble((String) imp_srcs[srcChannel].getProperty("saturation_2"))};
						  saturation_imp[srcChannel] = new boolean[width*height];
						  System.out.println(String.format("channel %d saturations = %6.2f %6.2f %6.2f %6.2f", srcChannel,
								  saturations[0],saturations[1],saturations[2],saturations[3]));
						  double [] scaled_saturations = new double [saturations.length];
						  for (int i = 0; i < scaled_saturations.length; i++){
							  scaled_saturations[i] = saturations[i] * clt_parameters.sat_level;
						  }
						  for (int y = 0; y < height-1; y+=2){
							  for (int x = 0; x < width-1; x+=2){
								  if (pixels[y*width+x        ] > scaled_saturations[0])  saturation_imp[srcChannel][y*width+x        ] = true;
								  if (pixels[y*width+x+      1] > scaled_saturations[1])  saturation_imp[srcChannel][y*width+x      +1] = true;
								  if (pixels[y*width+x+width  ] > scaled_saturations[2])  saturation_imp[srcChannel][y*width+x+width  ] = true;
								  if (pixels[y*width+x+width+1] > scaled_saturations[3])  saturation_imp[srcChannel][y*width+x+width+1] = true;
							  }
						  }
					  }



					  if (this.correctionsParameters.vignetting){
						  if ((eyesisCorrections.channelVignettingCorrection==null) || (srcChannel<0) || (srcChannel>=eyesisCorrections.channelVignettingCorrection.length) || (eyesisCorrections.channelVignettingCorrection[srcChannel]==null)){
							  System.out.println("No vignetting data for channel "+srcChannel);
							  return;
						  }
///						  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
						  if (pixels.length!=eyesisCorrections.channelVignettingCorrection[srcChannel].length){
							  System.out.println("Vignetting data for channel "+srcChannel+" has "+eyesisCorrections.channelVignettingCorrection[srcChannel].length+" pixels, image "+sourceFiles[nFile]+" has "+pixels.length);
							  return;
						  }
						  // TODO: Move to do it once:
						  double min_non_zero = 0.0;
						  for (int i=0;i<pixels.length;i++){
							  double d = eyesisCorrections.channelVignettingCorrection[srcChannel][i];
							  if ((d > 0.0) && ((min_non_zero == 0) || (min_non_zero > d))){
								  min_non_zero = d;
							  }
						  }
						  double max_vign_corr = clt_parameters.vignetting_range*min_non_zero;

						  System.out.println("Vignetting data: channel="+srcChannel+", min = "+min_non_zero);
						  for (int i=0;i<pixels.length;i++){
							  double d = eyesisCorrections.channelVignettingCorrection[srcChannel][i];
							  if (d > max_vign_corr) d = max_vign_corr;
							  pixels[i]*=d;
						  }
						  // Scale here, combine with vignetting later?
///						  int width =  imp_srcs[srcChannel].getWidth();
///						  int height = imp_srcs[srcChannel].getHeight();
						  for (int y = 0; y < height-1; y+=2){
							  for (int x = 0; x < width-1; x+=2){
								  pixels[y*width+x        ] *= clt_parameters.scale_g;
								  pixels[y*width+x+width+1] *= clt_parameters.scale_g;
								  pixels[y*width+x      +1] *= clt_parameters.scale_r;
								  pixels[y*width+x+width  ] *= clt_parameters.scale_b;
							  }
						  }

					  } else { // assuming GR/BG pattern
						  System.out.println("Applying fixed color gain correction parameters: Gr="+
								  clt_parameters.novignetting_r+", Gg="+clt_parameters.novignetting_g+", Gb="+clt_parameters.novignetting_b);
///						  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
///						  int width =  imp_srcs[srcChannel].getWidth();
///						  int height = imp_srcs[srcChannel].getHeight();
						  double kr = clt_parameters.scale_r/clt_parameters.novignetting_r;
						  double kg = clt_parameters.scale_g/clt_parameters.novignetting_g;
						  double kb = clt_parameters.scale_b/clt_parameters.novignetting_b;
						  for (int y = 0; y < height-1; y+=2){
							  for (int x = 0; x < width-1; x+=2){
								  pixels[y*width+x        ] *= kg;
								  pixels[y*width+x+width+1] *= kg;
								  pixels[y*width+x      +1] *= kr;
								  pixels[y*width+x+width  ] *= kb;
							  }
						  }
					  }
				  }
			  }
			  // may need to equalize gains between channels
			  // may need to equalize gains between channels
			  if (clt_parameters.gain_equalize || clt_parameters.colors_equalize){
				  channelGainsEqualize(
						  clt_parameters.gain_equalize,
						  clt_parameters.colors_equalize,
						  clt_parameters.nosat_equalize, // boolean nosat_equalize,
						  channelFiles,
						  imp_srcs,
						  saturation_imp, // boolean[][] saturated,
						  setNames.get(nSet), // just for debug messeges == setNames.get(nSet)
						  debugLevel);
			  }
			  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
				  int nFile=channelFiles[srcChannel];
				  if (nFile >=0){
					  processCLTSetImage( // returns ImagePlus, but it already should be saved/shown
							  imp_srcs[srcChannel], // should have properties "name"(base for saving results), "channel","path"
							  clt_parameters,
							  debayerParameters,
//							  nonlinParameters,
							  colorProcParameters,
							  channelGainParameters,
							  rgbParameters,
//							  convolveFFTSize, // 128 - fft size, kernel size should be size/2
							  scaleExposure[srcChannel],
							  threadsMax,  // maximal number of threads to launch
							  updateStatus,
							  debugLevel);
					  // warp result (add support for different color modes)
					  if (this.correctionsParameters.equirectangular){
						  if (equirectangularParameters.clearFullMap) eyesisCorrections.pixelMapping.deleteEquirectangularMapFull(srcChannel); // save memory? //removeUnusedSensorData - no, use equirectangular specific settings
						  if (equirectangularParameters.clearAllMaps) eyesisCorrections.pixelMapping.deleteEquirectangularMapAll(srcChannel); // save memory? //removeUnusedSensorData - no, use equirectangular specific settings
					  }
					  //pixelMapping
					  Runtime.getRuntime().gc();
					  if (debugLevel >-1) System.out.println("Processing image "+(iImage+1)+" (of "+fileIndices.length+") finished at "+
							  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
					  if (eyesisCorrections.stopRequested.get()>0) {
						  System.out.println("User requested stop");
						  return;
					  }
					  iImage++;
				  }
			  }
		  }
		  System.out.println("processCLTSets(): processing "+fileIndices.length+" files finished at "+
				  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
	  }

	  public ImagePlus processCLTSetImage(
			  ImagePlus imp_src, // should have properties "name"(base for saving results), "channel","path"
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
//			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
//			  int          convolveFFTSize, // 128 - fft size, kernel size should be size/2
			  double 		     scaleExposure,
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel){
		  boolean advanced=this.correctionsParameters.zcorrect || this.correctionsParameters.equirectangular;
//		  boolean crop=      advanced? true: this.correctionsParameters.crop;
		  boolean rotate=    advanced? false: this.correctionsParameters.rotate;
		  double JPEG_scale= advanced? 1.0: this.correctionsParameters.JPEG_scale;
		  boolean toRGB=     advanced? true: this.correctionsParameters.toRGB;
		  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?

		  // may use this.StartTime to report intermediate steps execution times
		  String name=(String) imp_src.getProperty("name");
		  //		int channel= Integer.parseInt((String) imp_src.getProperty("channel"));
		  int channel= (Integer) imp_src.getProperty("channel");
		  String path= (String) imp_src.getProperty("path");

		  String title=name+"-"+String.format("%02d", channel);
		  ImagePlus result=imp_src;
		  if (debugLevel>1) System.out.println("processing: "+path);
		  result.setTitle(title+"RAW");
		  if (!this.correctionsParameters.split){
			  eyesisCorrections.saveAndShow(result, this.correctionsParameters);
			  return result;
		  }
		  // Generate split parameters for DCT processing mode
		  // Split into Bayer components, oversample, increase canvas

		  double [][] double_stack = eyesisCorrections.bayerToDoubleStack(
				  result, // source Bayer image, linearized, 32-bit (float))
				  null); // no margins, no oversample

		  String titleFull=title+"-SPLIT";
		  if (debugLevel > -1){
			  double [] chn_avg = {0.0,0.0,0.0};
			  int width =  imp_src.getWidth();
			  int height = imp_src.getHeight();
			  for (int c = 0; c < 3; c++){
				  for (int i = 0; i<double_stack[c].length; i++){
					  chn_avg[c] += double_stack[c][i];
				  }
			  }
			  chn_avg[0] /= width*height/4;
			  chn_avg[1] /= width*height/4;
			  chn_avg[2] /= width*height/2;
			  System.out.println("Split channels averages: R="+chn_avg[0]+", G="+chn_avg[2]+", B="+chn_avg[1]);
		  }
		  String [] rbg_titles = {"Red", "Blue", "Green"};
		  ImageStack stack;
		  if (!this.correctionsParameters.debayer) {
//			  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
//			  ImageStack
			  stack = sdfa_instance.makeStack(double_stack, imp_src.getWidth(), imp_src.getHeight(), rbg_titles);
			  result= new ImagePlus(titleFull, stack);
			  eyesisCorrections.saveAndShow(result, this.correctionsParameters);
			  return result;
		  }
		  // =================
		  if (debugLevel > 0) {
			  System.out.println("Showing image BEFORE_CLT_PROC");
			  sdfa_instance.showArrays(double_stack,  imp_src.getWidth(), imp_src.getHeight(), true, "BEFORE_CLT_PROC", rbg_titles);
		  }
		  if (this.correctionsParameters.deconvolve) { // process with DCT, otherwise use simple debayer
			  ImageDtt image_dtt = new ImageDtt();
/*
			  double [][][][][] clt_data = image_dtt.cltStack(
					  stack,
					  channel,
					  clt_parameters,
					  clt_parameters.ishift_x, //final int shiftX, // shift image horizontally (positive - right)
					  clt_parameters.ishift_y, //final int shiftY, // shift image vertically (positive - down)
					  threadsMax,
					  debugLevel,
					  updateStatus);
*/
			  for (int i =0 ; i < double_stack[0].length; i++){
//				  double_stack[0][i]*=2.0; // Scale red twice to compensate less pixels than green
//				  double_stack[1][i]*=2.0; // Scale blue twice to compensate less pixels than green
				  double_stack[2][i]*=0.5; // Scale blue twice to compensate less pixels than green
			  }
			  double [][][][][] clt_data = image_dtt.clt_aberrations(
					  double_stack,                 // final double [][]       imade_data,
					  imp_src.getWidth(),           //	final int               width,
					  clt_kernels[channel],         // final double [][][][][] clt_kernels, // [color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
					  clt_parameters.kernel_step,
					  clt_parameters.transform_size,
					  clt_parameters.clt_window,
					  clt_parameters.shift_x,       // final int               shiftX, // shift image horizontally (positive - right) - just for testing
					  clt_parameters.shift_y,       // final int               shiftY, // shift image vertically (positive - down)
					  clt_parameters.tileX,         // final int               debug_tileX,
					  clt_parameters.tileY,         // final int               debug_tileY,
					  (clt_parameters.dbg_mode & 64) != 0, // no fract shift
					  (clt_parameters.dbg_mode & 128) != 0, // no convolve
					  (clt_parameters.dbg_mode & 256) != 0, // transpose convolve
					  threadsMax,
					  debugLevel);
//					  updateStatus);


			  System.out.println("clt_data.length="+clt_data.length+" clt_data[0].length="+clt_data[0].length
					  +" clt_data[0][0].length="+clt_data[0][0].length+" clt_data[0][0][0].length="+
					  clt_data[0][0][0].length+" clt_data[0][0][0][0].length="+clt_data[0][0][0][0].length);
			  /*
			  if (dct_parameters.color_DCT){ // convert RBG -> YPrPb
				  dct_data = image_dtt.dct_color_convert(
						  dct_data,
						  colorProcParameters.kr,
						  colorProcParameters.kb,
						  dct_parameters.sigma_rb,        // blur of channels 0,1 (r,b) in addition to 2 (g)
						  dct_parameters.sigma_y,         // blur of Y from G
						  dct_parameters.sigma_color,     // blur of Pr, Pb in addition to Y
						  threadsMax,
						  debugLevel);
			  } else { // just LPF RGB
			  */
				  if (clt_parameters.corr_sigma > 0){ // no filter at all
					  for (int chn = 0; chn < clt_data.length; chn++) {
						  image_dtt.clt_lpf(
								  clt_parameters.corr_sigma,
								  clt_data[chn],
								  clt_parameters.transform_size,
								  threadsMax,
								  debugLevel);
					  }
				  }
/*
			  }
*/
			  int tilesY = imp_src.getHeight()/clt_parameters.transform_size;
			  int tilesX = imp_src.getWidth()/clt_parameters.transform_size;
			  if (debugLevel > 0){
				  System.out.println("--tilesX="+tilesX);
				  System.out.println("--tilesY="+tilesY);
			  }
			  if (debugLevel > 1){
			        double [][] clt = new double [clt_data.length*4][];
			        for (int chn = 0; chn < clt_data.length; chn++) {
			        	double [][] clt_set = image_dtt.clt_dbg(
			        			clt_data [chn],
								  threadsMax,
								  debugLevel);
			        	for (int ii = 0; ii < clt_set.length; ii++) clt[chn*4+ii] = clt_set[ii];
			        }

			        if (debugLevel > 0){
			        	sdfa_instance.showArrays(clt,
			        			tilesX*clt_parameters.transform_size,
			        			tilesY*clt_parameters.transform_size,
			        			true,
			        			result.getTitle()+"-CLT");
			        }
			  }
			  double [][] iclt_data = new double [clt_data.length][];
			  for (int chn=0; chn<clt_data.length;chn++){
				  iclt_data[chn] = image_dtt.iclt_2d(
						  clt_data[chn],                  // scanline representation of dcd data, organized as dct_size x dct_size tiles
						  clt_parameters.transform_size,  // final int
						  clt_parameters.clt_window,      // window_type
						  15,                             // clt_parameters.iclt_mask,       //which of 4 to transform back
						  0,                              // clt_parameters.dbg_mode,        //which of 4 to transform back
						  threadsMax,
						  debugLevel);

			  }

//					  if (debugLevel > -1) System.out.println("Applyed LPF, sigma = "+dct_parameters.dbg_sigma);
					  if (debugLevel > 0) sdfa_instance.showArrays(
							  iclt_data,
							  (tilesX + 1) * clt_parameters.transform_size,
							  (tilesY + 1) * clt_parameters.transform_size,
							  true,
							  result.getTitle()+"-rbg_sigma");
				/*
				  }
			  }
			 */
			  if (debugLevel > 0) sdfa_instance.showArrays(iclt_data,
					  (tilesX + 0) * clt_parameters.transform_size,
					  (tilesY + 0) * clt_parameters.transform_size,
					  true,
					  result.getTitle()+"-ICLT-RGB");

			  // convert to ImageStack of 3 slices
			  String [] sliceNames = {"red", "blue", "green"};
			  stack = sdfa_instance.makeStack(
					  iclt_data,
					  (tilesX + 0) * clt_parameters.transform_size,
					  (tilesY + 0) * clt_parameters.transform_size,
					  sliceNames); // or use null to get chn-nn slice names


		  } else { // if (this.correctionsParameters.deconvolve) - here use a simple debayer
			  System.out.println("Bypassing CLT-based aberration correction");
			  stack = sdfa_instance.makeStack(double_stack, imp_src.getWidth(), imp_src.getHeight(), rbg_titles);
			  debayer_rbg(stack, 0.25); // simple standard 3x3 kernel debayer
		  }
		  if (debugLevel > -1){
			  double [] chn_avg = {0.0,0.0,0.0};
			  float [] pixels;
			  int width =  stack.getWidth();
			  int height = stack.getHeight();

			  for (int c = 0; c <3; c++){
				  pixels = (float[]) stack.getPixels(c+1);
				  for (int i = 0; i<pixels.length; i++){
					  chn_avg[c] += pixels[i];
				  }
			  }
			  chn_avg[0] /= width*height;
			  chn_avg[1] /= width*height;
			  chn_avg[2] /= width*height;
			  System.out.println("Processed channels averages: R="+chn_avg[0]+", G="+chn_avg[2]+", B="+chn_avg[1]);
		  }

		  if (!this.correctionsParameters.colorProc){
			  result= new ImagePlus(titleFull, stack);
			  eyesisCorrections.saveAndShow(
					  result,
					  this.correctionsParameters);
			  return result;
		  }
		  if (debugLevel > 1) System.out.println("before colors.1");
		  //Processing colors - changing stack sequence to r-g-b (was r-b-g)
		  if (!eyesisCorrections.fixSliceSequence(
				  stack,
				  debugLevel)){
			  if (debugLevel > -1) System.out.println("fixSliceSequence() returned false");
			  return null;
		  }
		  if (debugLevel > 1) System.out.println("before colors.2");
		  if (debugLevel > 1){
			  ImagePlus imp_dbg=new ImagePlus(imp_src.getTitle()+"-"+channel+"-preColors",stack);
			  eyesisCorrections.saveAndShow(
					  imp_dbg,
					  this.correctionsParameters);
		  }
		  if (debugLevel > 1) System.out.println("before colors.3, scaleExposure="+scaleExposure+" scale = "+(255.0/eyesisCorrections.psfSubpixelShouldBe4/eyesisCorrections.psfSubpixelShouldBe4/scaleExposure));
		  CorrectionColorProc correctionColorProc=new CorrectionColorProc(eyesisCorrections.stopRequested);
		  double [][] yPrPb=new double [3][];
		  //			if (dct_parameters.color_DCT){
		  // need to get YPbPr - not RGB here
		  //			} else {
		  correctionColorProc.processColorsWeights(stack, // just gamma convert? TODO: Cleanup? Convert directly form the linear YPrPb
				  //					  255.0/this.psfSubpixelShouldBe4/this.psfSubpixelShouldBe4, //  double scale,     // initial maximal pixel value (16))
				  //					  255.0/eyesisCorrections.psfSubpixelShouldBe4/eyesisCorrections.psfSubpixelShouldBe4/scaleExposure, //  double scale,     // initial maximal pixel value (16))
				  //					  255.0/2/2/scaleExposure, //  double scale,     // initial maximal pixel value (16))
				  255.0/scaleExposure, //  double scale,     // initial maximal pixel value (16))
				  colorProcParameters,
				  channelGainParameters,
				  channel,
				  null, //correctionDenoise.getDenoiseMask(),
				  this.correctionsParameters.blueProc,
				  debugLevel);
		  if (debugLevel > 1) System.out.println("Processed colors to YPbPr, total number of slices="+stack.getSize());
		  if (debugLevel > 1) {
			  ImagePlus imp_dbg=new ImagePlus("procColors",stack);
			  eyesisCorrections.saveAndShow(
					  imp_dbg,
					  this.correctionsParameters);
		  }
		  float [] fpixels;
		  int [] slices_YPrPb = {8,6,7};
		  yPrPb=new double [3][];
		  for (int n = 0; n < slices_YPrPb.length; n++){
			  fpixels = (float[]) stack.getPixels(slices_YPrPb[n]);
			  yPrPb[n] = new double [fpixels.length];
			  for (int i = 0; i < fpixels.length; i++) yPrPb[n][i] = fpixels[i];
		  }

		  if (toRGB) {
			  if (debugLevel > 0){
				  System.out.println("correctionColorProc.YPrPbToRGB");
			  }
			  stack =  YPrPbToRGB(yPrPb,
					  colorProcParameters.kr,        // 0.299;
					  colorProcParameters.kb,        // 0.114;
					  stack.getWidth());

			  title=titleFull; // including "-DECONV" or "-COMBO"
			  titleFull=title+"-RGB-float";
			  //Trim stack to just first 3 slices
			  if (debugLevel > 1){ // 2){
				  ImagePlus imp_dbg=new ImagePlus("YPrPbToRGB",stack);
				  eyesisCorrections.saveAndShow(
						  imp_dbg,
						  this.correctionsParameters);
			  }
			  while (stack.getSize() > 3) stack.deleteLastSlice();
			  if (debugLevel > 1) System.out.println("Trimming color stack");
		  } else {
			  title=titleFull; // including "-DECONV" or "-COMBO"
			  titleFull=title+"-YPrPb"; // including "-DECONV" or "-COMBO"
			  if (debugLevel > 1) System.out.println("Using full stack, including YPbPr");
		  }

		  result= new ImagePlus(titleFull, stack);
		  // Crop image to match original one (scaled to oversampling)
/*
		  if (crop){ // always crop if equirectangular
			  if (debugLevel > 1) System.out.println("cropping");
			  stack = eyesisCorrections.cropStack32(stack,splitParameters);
			  if (debugLevel > 2) { // 2){
				  ImagePlus imp_dbg=new ImagePlus("cropped",stack);
				  eyesisCorrections.saveAndShow(
						  imp_dbg,
						  this.correctionsParameters);
			  }
		  }
*/
		  // rotate the result
		  if (rotate){ // never rotate for equirectangular
			  stack=eyesisCorrections.rotateStack32CW(stack);
		  }
		  if (!toRGB && !this.correctionsParameters.jpeg){ // toRGB set for equirectangular
			  if (debugLevel > 1) System.out.println("!toRGB && !this.correctionsParameters.jpeg");
			  eyesisCorrections.saveAndShow(result, this.correctionsParameters);
			  return result;
		  } else { // that's not the end result, save if required
			  if (debugLevel > 1) System.out.println("!toRGB && !this.correctionsParameters.jpeg - else");
			  eyesisCorrections.saveAndShow(result,
					  eyesisCorrections.correctionsParameters,
					  eyesisCorrections.correctionsParameters.save32,
					  false,
					  eyesisCorrections.correctionsParameters.JPEG_quality); // save, no show
		  }
		  // convert to RGB48 (16 bits per color component)
		  ImagePlus imp_RGB;
		  stack=eyesisCorrections.convertRGB32toRGB16Stack(
				  stack,
				  rgbParameters);

		  titleFull=title+"-RGB48";
		  result= new ImagePlus(titleFull, stack);
		  //			  ImagePlus imp_RGB24;
		  result.updateAndDraw();
		  if (debugLevel > 1) System.out.println("result.updateAndDraw(), "+titleFull+"-RGB48");

		  CompositeImage compositeImage=eyesisCorrections.convertToComposite(result);

		  if (!this.correctionsParameters.jpeg && !advanced){ // RGB48 was the end result
			  if (debugLevel > 1) System.out.println("if (!this.correctionsParameters.jpeg && !advanced)");
			  eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters);
			  return result;
		  } else { // that's not the end result, save if required
			  if (debugLevel > 1) System.out.println("if (!this.correctionsParameters.jpeg && !advanced) - else");
			  eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters, this.correctionsParameters.save16, false); // save, no show
			  //				  eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters, this.correctionsParameters.save16, true); // save, no show
		  }

		  imp_RGB=eyesisCorrections.convertRGB48toRGB24(
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
		  eyesisCorrections.saveAndShow(imp_RGB, this.correctionsParameters);

		  return result;
	  }

	  public void processCLTQuads(
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
//			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
//			  int          convolveFFTSize, // 128 - fft size, kernel size should be size/2
			  final int          threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel)
	  {
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
					  // removeUnusedSensorData should be off!?
					  channels=this.eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (eyesisCorrections.isChannelEnabled(channels[i])){
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
		  double [] referenceExposures=eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		  int [][] fileIndices=new int [numImagesToProcess][2]; // file index, channel number
		  int index=0;
		  for (int nFile=0;nFile<enabledFiles.length;nFile++){
			  if ((sourceFiles[nFile]!=null) && (sourceFiles[nFile].length()>1)) {
				  int [] channels={correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile])};
				  if (correctionsParameters.isJP4()){
					  int subCamera= channels[0]- correctionsParameters.firstSubCamera; // to match those in the sensor files
					  channels=eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (eyesisCorrections.isChannelEnabled(channels[i])){
						  fileIndices[index  ][0]=nFile;
						  fileIndices[index++][1]=channels[i];
					  }
				  }
			  }
		  }
		  ArrayList<String> setNames = new ArrayList<String>();
		  ArrayList<ArrayList<Integer>> setFiles = new ArrayList<ArrayList<Integer>>();

		  for (int iImage=0;iImage<fileIndices.length;iImage++){
			  int nFile=fileIndices[iImage][0];
			  String setName = correctionsParameters.getNameFromSourceTiff(sourceFiles[nFile]);
			  if (!setNames.contains(setName)) {
				  setNames.add(setName);
				  setFiles.add(new ArrayList<Integer>());
			  }
			  setFiles.get(setNames.indexOf(setName)).add(new Integer(nFile));
		  }


		  for (int nSet = 0; nSet < setNames.size(); nSet++){
			  int maxChn = 0;
			  for (int i = 0; i < setFiles.get(nSet).size(); i++){
				  int chn = fileIndices[setFiles.get(nSet).get(i)][1];
				  if (chn > maxChn) maxChn = chn;
			  }
			  int [] channelFiles = new int[maxChn+1];
			  for (int i =0; i < channelFiles.length; i++) channelFiles[i] = -1;
			  for (int i = 0; i < setFiles.get(nSet).size(); i++){
				  channelFiles[fileIndices[setFiles.get(nSet).get(i)][1]] = setFiles.get(nSet).get(i);
			  }

			  ImagePlus [] imp_srcs = new ImagePlus[channelFiles.length];
			  double [] scaleExposures = new double[channelFiles.length];
			  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
				  int nFile=channelFiles[srcChannel];
				  imp_srcs[srcChannel]=null;
				  if (nFile >=0){
					  if (correctionsParameters.isJP4()){
						  int subchannel=eyesisCorrections.pixelMapping.getSubChannel(srcChannel);
						  if (this.correctionsParameters.swapSubchannels01) {
							  switch (subchannel){
							  case 0: subchannel=1; break;
							  case 1: subchannel=0; break;
							  }
						  }
						  if (debugLevel>0) System.out.println("Processing set " + setNames.get(nSet)+" channel "+srcChannel+" - subchannel "+subchannel+" of "+sourceFiles[nFile]);
						  ImagePlus imp_composite=eyesisCorrections.JP4_INSTANCE.open(
								  "", // path,
								  sourceFiles[nFile],
								  "",  //arg - not used in JP46 reader
								  true, // un-apply camera color gains
								  null, // new window
								  false); // do not show
						  imp_srcs[srcChannel]=eyesisCorrections.JP4_INSTANCE.demuxImage(imp_composite, subchannel);
						  if (imp_srcs[srcChannel] == null) imp_srcs[srcChannel] = imp_composite; // not a composite image
						  // do we need to add any properties?
					  } else {
						  imp_srcs[srcChannel]=new ImagePlus(sourceFiles[nFile]);
						  //					  (new JP46_Reader_camera(false)).decodeProperiesFromInfo(imp_src); // decode existent properties from info
						  eyesisCorrections.JP4_INSTANCE.decodeProperiesFromInfo(imp_srcs[srcChannel]); // decode existent properties from info
						  if (debugLevel>0) System.out.println("Processing "+sourceFiles[nFile]);
					  }
					  imp_srcs[srcChannel] =  padBayerToFullSize(
							  imp_srcs[srcChannel], // ImagePlus imp_src,
							  eyesisCorrections.pixelMapping.sensors[srcChannel].getSensorWH(),
							  true); // boolean replicate);
					  scaleExposures[srcChannel] = 1.0;
					  if (!Double.isNaN(referenceExposures[nFile]) && (imp_srcs[srcChannel].getProperty("EXPOSURE")!=null)){
						  scaleExposures[srcChannel] = referenceExposures[nFile]/Double.parseDouble((String) imp_srcs[srcChannel].getProperty("EXPOSURE"));
						  if (debugLevel>0) System.out.println("Will scale intensity (to compensate for exposure) by  "+scaleExposures[srcChannel]);
					  }
					  imp_srcs[srcChannel].setProperty("name",    correctionsParameters.getNameFromSourceTiff(sourceFiles[nFile]));
					  imp_srcs[srcChannel].setProperty("channel", srcChannel); // it may already have channel
					  imp_srcs[srcChannel].setProperty("path",    sourceFiles[nFile]); // it may already have channel

					  if (this.correctionsParameters.pixelDefects && (eyesisCorrections.defectsXY!=null)&& (eyesisCorrections.defectsXY[srcChannel]!=null)){
						  // apply pixel correction
						  int numApplied=	eyesisCorrections.correctDefects(
								  imp_srcs[srcChannel],
								  srcChannel,
								  debugLevel);
						  if ((debugLevel>0) && (numApplied>0)) { // reduce verbosity after verified defect correction works
							  System.out.println("Corrected "+numApplied+" pixels in "+sourceFiles[nFile]);
						  }
					  }

					  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
					  int width =  imp_srcs[srcChannel].getWidth();
					  int height = imp_srcs[srcChannel].getHeight();

					  if (clt_parameters.sat_level > 0.0){
						  double [] saturations = {
								  Double.parseDouble((String) imp_srcs[srcChannel].getProperty("saturation_1")),
								  Double.parseDouble((String) imp_srcs[srcChannel].getProperty("saturation_0")),
								  Double.parseDouble((String) imp_srcs[srcChannel].getProperty("saturation_3")),
								  Double.parseDouble((String) imp_srcs[srcChannel].getProperty("saturation_2"))};
						  saturation_imp[srcChannel] = new boolean[width*height];
						  System.out.println(String.format("channel %d saturations = %6.2f %6.2f %6.2f %6.2f", srcChannel,
								  saturations[0],saturations[1],saturations[2],saturations[3]));
						  double [] scaled_saturations = new double [saturations.length];
						  for (int i = 0; i < scaled_saturations.length; i++){
							  scaled_saturations[i] = saturations[i] * clt_parameters.sat_level;
						  }
						  for (int y = 0; y < height-1; y+=2){
							  for (int x = 0; x < width-1; x+=2){
								  if (pixels[y*width+x        ] > scaled_saturations[0])  saturation_imp[srcChannel][y*width+x        ] = true;
								  if (pixels[y*width+x+      1] > scaled_saturations[1])  saturation_imp[srcChannel][y*width+x      +1] = true;
								  if (pixels[y*width+x+width  ] > scaled_saturations[2])  saturation_imp[srcChannel][y*width+x+width  ] = true;
								  if (pixels[y*width+x+width+1] > scaled_saturations[3])  saturation_imp[srcChannel][y*width+x+width+1] = true;
							  }
						  }
					  }


					  if (this.correctionsParameters.vignetting){
						  if ((eyesisCorrections.channelVignettingCorrection==null) || (srcChannel<0) || (srcChannel>=eyesisCorrections.channelVignettingCorrection.length) || (eyesisCorrections.channelVignettingCorrection[srcChannel]==null)){
							  System.out.println("No vignetting data for channel "+srcChannel);
							  return;
						  }
///						  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
						  if (pixels.length!=eyesisCorrections.channelVignettingCorrection[srcChannel].length){
							  System.out.println("Vignetting data for channel "+srcChannel+" has "+eyesisCorrections.channelVignettingCorrection[srcChannel].length+" pixels, image "+sourceFiles[nFile]+" has "+pixels.length);
							  return;
						  }
						  // TODO: Move to do it once:
						  double min_non_zero = 0.0;
						  for (int i=0;i<pixels.length;i++){
							  double d = eyesisCorrections.channelVignettingCorrection[srcChannel][i];
							  if ((d > 0.0) && ((min_non_zero == 0) || (min_non_zero > d))){
								  min_non_zero = d;
							  }
						  }
						  double max_vign_corr = clt_parameters.vignetting_range*min_non_zero;

						  System.out.println("Vignetting data: channel="+srcChannel+", min = "+min_non_zero);
						  for (int i=0;i<pixels.length;i++){
							  double d = eyesisCorrections.channelVignettingCorrection[srcChannel][i];
							  if (d > max_vign_corr) d = max_vign_corr;
							  pixels[i]*=d;
						  }
						  // Scale here, combine with vignetting later?
///						  int width =  imp_srcs[srcChannel].getWidth();
///						  int height = imp_srcs[srcChannel].getHeight();
						  for (int y = 0; y < height-1; y+=2){
							  for (int x = 0; x < width-1; x+=2){
								  pixels[y*width+x        ] *= clt_parameters.scale_g;
								  pixels[y*width+x+width+1] *= clt_parameters.scale_g;
								  pixels[y*width+x      +1] *= clt_parameters.scale_r;
								  pixels[y*width+x+width  ] *= clt_parameters.scale_b;
							  }
						  }

					  } else { // assuming GR/BG pattern
						  System.out.println("Applying fixed color gain correction parameters: Gr="+
								  clt_parameters.novignetting_r+", Gg="+clt_parameters.novignetting_g+", Gb="+clt_parameters.novignetting_b);
///						  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
///						  int width =  imp_srcs[srcChannel].getWidth();
///						  int height = imp_srcs[srcChannel].getHeight();
						  double kr = clt_parameters.scale_r/clt_parameters.novignetting_r;
						  double kg = clt_parameters.scale_g/clt_parameters.novignetting_g;
						  double kb = clt_parameters.scale_b/clt_parameters.novignetting_b;
						  for (int y = 0; y < height-1; y+=2){
							  for (int x = 0; x < width-1; x+=2){
								  pixels[y*width+x        ] *= kg;
								  pixels[y*width+x+width+1] *= kg;
								  pixels[y*width+x      +1] *= kr;
								  pixels[y*width+x+width  ] *= kb;
							  }
						  }
					  }
				  }
			  }
			  // once per quad here
			  // may need to equalize gains between channels
			  // may need to equalize gains between channels
			  if (clt_parameters.gain_equalize || clt_parameters.colors_equalize){
				  channelGainsEqualize(
						  clt_parameters.gain_equalize,
						  clt_parameters.colors_equalize,
						  clt_parameters.nosat_equalize, // boolean nosat_equalize,
						  channelFiles,
						  imp_srcs,
						  saturation_imp, // boolean[][] saturated,
						  setNames.get(nSet), // just for debug messeges == setNames.get(nSet)
						  debugLevel);
			  }
			  // once per quad here
			  processCLTQuad( // returns ImagePlus, but it already should be saved/shown
					  imp_srcs, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
					  clt_parameters,
					  debayerParameters,
//					  nonlinParameters,
					  colorProcParameters,
					  channelGainParameters,
					  rgbParameters,
//					  convolveFFTSize, // 128 - fft size, kernel size should be size/2
					  scaleExposures,
					  threadsMax,  // maximal number of threads to launch
					  updateStatus,
					  debugLevel);
			  Runtime.getRuntime().gc();
			  if (debugLevel >-1) System.out.println("Processing set "+(nSet+1)+" (of "+setNames.size()+") finished at "+
					  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
			  if (eyesisCorrections.stopRequested.get()>0) {
				  System.out.println("User requested stop");
				  return;
			  }
		  }
		  System.out.println("processCLTQuads(): processing "+fileIndices.length+" files finished at "+
				  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
	  }

	  public ImagePlus [] processCLTQuad(
			  ImagePlus [] imp_quad, // should have properties "name"(base for saving results), "channel","path"
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
//			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
//			  int          convolveFFTSize, // 128 - fft size, kernel size should be size/2
			  double []	       scaleExposures, // probably not needed here
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel){
		  boolean advanced=this.correctionsParameters.zcorrect || this.correctionsParameters.equirectangular;
//		  boolean crop=      advanced? true: this.correctionsParameters.crop;
		  boolean rotate=    advanced? false: this.correctionsParameters.rotate;
		  double JPEG_scale= advanced? 1.0: this.correctionsParameters.JPEG_scale;
		  boolean toRGB=     advanced? true: this.correctionsParameters.toRGB;
		  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?

		  // may use this.StartTime to report intermediate steps execution times
		  String name=(String) imp_quad[0].getProperty("name");
		  //		int channel= Integer.parseInt((String) imp_src.getProperty("channel"));
		  int channel= (Integer) imp_quad[0].getProperty("channel");
		  String path= (String) imp_quad[0].getProperty("path");

		  ImagePlus [] results = new ImagePlus[imp_quad.length];
		  for (int i = 0; i < results.length; i++) {
			  results[i] = imp_quad[i];
			  results[i].setTitle(results[i].getTitle()+"RAW");
		  }
		  if (debugLevel>1) System.out.println("processing: "+path);
		  double [][][] double_stacks = new double [imp_quad.length][][];
		  for (int i = 0; i < double_stacks.length; i++){
			  double_stacks[i] = eyesisCorrections.bayerToDoubleStack(
					  imp_quad[i], // source Bayer image, linearized, 32-bit (float))
					  null); // no margins, no oversample
		  }

//		  String [] rbg_titles = {"Red", "Blue", "Green"};
		  ImageStack stack;
		  // =================
		  ImageDtt image_dtt = new ImageDtt();
		  for (int i = 0; i < double_stacks.length; i++){
			  for (int j =0 ; j < double_stacks[i][0].length; j++){
				  double_stacks[i][2][j]*=0.5; // Scale green 0.5 to compensate more pixels than R,B
			  }
		  }
		  double [][][][][][] clt_data = image_dtt.clt_aberrations_quad(
				  clt_parameters.disparity,     // final double            disparity,
				  double_stacks,                // final double [][][]      imade_data, // first index - number of image in a quad
				  imp_quad[0].getWidth(),       //	final int               width,
				  geometryCorrection,           // final GeometryCorrection  geometryCorrection,
				  clt_kernels,                  // final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				  clt_parameters.kernel_step,
				  clt_parameters.transform_size,
				  clt_parameters.clt_window,
				  clt_parameters.shift_x,       // final int               shiftX, // shift image horizontally (positive - right) - just for testing
				  clt_parameters.shift_y,       // final int               shiftY, // shift image vertically (positive - down)
				  clt_parameters.tileX,         // final int               debug_tileX,
				  clt_parameters.tileY,         // final int               debug_tileY,
				  (clt_parameters.dbg_mode & 64) != 0, // no fract shift
				  (clt_parameters.dbg_mode & 128) != 0, // no convolve
				  (clt_parameters.dbg_mode & 256) != 0, // transpose convolve
				  threadsMax,
				  debugLevel);

		  System.out.println("clt_data.length="+clt_data.length+" clt_data[0].length="+clt_data[0].length
				  +" clt_data[0][0].length="+clt_data[0][0].length+" clt_data[0][0][0].length="+
				  clt_data[0][0][0].length+" clt_data[0][0][0][0].length="+clt_data[0][0][0][0].length+
				  " clt_data[0][0][0][0][0].length="+clt_data[0][0][0][0][0].length);

		  for (int iQuad = 0; iQuad <clt_data.length; iQuad++){

			  String title=name+"-"+String.format("%02d", iQuad);
			  String titleFull=title+"-SPLIT";

			  if (clt_parameters.corr_sigma > 0){ // no filter at all
				  for (int chn = 0; chn < clt_data[iQuad].length; chn++) {
					  image_dtt.clt_lpf(
							  clt_parameters.corr_sigma,
							  clt_data[iQuad][chn],
							  clt_parameters.transform_size,
							  threadsMax,
							  debugLevel);
				  }
			  }

			  int tilesY = imp_quad[iQuad].getHeight()/clt_parameters.transform_size;
			  int tilesX = imp_quad[iQuad].getWidth()/clt_parameters.transform_size;
			  if (debugLevel > 0){
				  System.out.println("--tp.tilesX="+tilesX);
				  System.out.println("--tp.tilesY="+tilesY);
			  }
			  if (debugLevel > 1){
				  double [][] clt = new double [clt_data[iQuad].length*4][];
				  for (int chn = 0; chn < clt_data[iQuad].length; chn++) {
					  double [][] clt_set = image_dtt.clt_dbg(
							  clt_data [iQuad][chn],
							  threadsMax,
							  debugLevel);
					  for (int ii = 0; ii < clt_set.length; ii++) clt[chn*4+ii] = clt_set[ii];
				  }

				  if (debugLevel > 0){
					  sdfa_instance.showArrays(clt,
							  tilesX*clt_parameters.transform_size,
							  tilesY*clt_parameters.transform_size,
							  true,
							  results[iQuad].getTitle()+"-CLT");
				  }
			  }
			  double [][] iclt_data = new double [clt_data[iQuad].length][];
			  for (int chn=0; chn<iclt_data.length;chn++){
				  iclt_data[chn] = image_dtt.iclt_2d(
						  clt_data[iQuad][chn],           // scanline representation of dcd data, organized as dct_size x dct_size tiles
						  clt_parameters.transform_size,  // final int
						  clt_parameters.clt_window,      // window_type
						  15,                             // clt_parameters.iclt_mask,       //which of 4 to transform back
						  0,                              // clt_parameters.dbg_mode,        //which of 4 to transform back
						  threadsMax,
						  debugLevel);

			  }
			  if (debugLevel > 0) sdfa_instance.showArrays(
					  iclt_data,
					  (tilesX + 0) * clt_parameters.transform_size,
					  (tilesY + 0) * clt_parameters.transform_size,
					  true,
					  results[iQuad].getTitle()+"-rbg_sigma");
			  if (debugLevel > 0) sdfa_instance.showArrays(iclt_data,
					  (tilesX + 0) * clt_parameters.transform_size,
					  (tilesY + 0) * clt_parameters.transform_size,
					  true,
					  results[iQuad].getTitle()+"-ICLT-RGB");

			  // convert to ImageStack of 3 slices
			  String [] sliceNames = {"red", "blue", "green"};
			  stack = sdfa_instance.makeStack(
					  iclt_data,
					  (tilesX + 0) * clt_parameters.transform_size,
					  (tilesY + 0) * clt_parameters.transform_size,
					  sliceNames); // or use null to get chn-nn slice names

			  if (debugLevel > -1){
				  double [] chn_avg = {0.0,0.0,0.0};
				  float [] pixels;
				  int width =  stack.getWidth();
				  int height = stack.getHeight();

				  for (int c = 0; c <3; c++){
					  pixels = (float[]) stack.getPixels(c+1);
					  for (int i = 0; i<pixels.length; i++){
						  chn_avg[c] += pixels[i];
					  }
				  }
				  chn_avg[0] /= width*height;
				  chn_avg[1] /= width*height;
				  chn_avg[2] /= width*height;
				  System.out.println("Processed channels averages: R="+chn_avg[0]+", G="+chn_avg[2]+", B="+chn_avg[1]);
			  }

			  if (!this.correctionsParameters.colorProc){
				  results[iQuad]= new ImagePlus(titleFull, stack);
				  eyesisCorrections.saveAndShow(
						  results[iQuad],
						  this.correctionsParameters);
				  continue; // return results;
			  }
			  if (debugLevel > 1) System.out.println("before colors.1");
			  //Processing colors - changing stack sequence to r-g-b (was r-b-g)
			  if (!eyesisCorrections.fixSliceSequence(
					  stack,
					  debugLevel)){
				  if (debugLevel > -1) System.out.println("fixSliceSequence() returned false");
				  return null;
			  }
			  if (debugLevel > 1) System.out.println("before colors.2");
			  if (debugLevel > 1){
				  ImagePlus imp_dbg=new ImagePlus(imp_quad[iQuad].getTitle()+"-"+channel+"-preColors",stack);
				  eyesisCorrections.saveAndShow(
						  imp_dbg,
						  this.correctionsParameters);
			  }
			  if (debugLevel > 1) System.out.println("before colors.3, scaleExposure="+scaleExposures[iQuad]+" scale = "+(255.0/eyesisCorrections.psfSubpixelShouldBe4/eyesisCorrections.psfSubpixelShouldBe4/scaleExposures[iQuad]));
			  CorrectionColorProc correctionColorProc=new CorrectionColorProc(eyesisCorrections.stopRequested);
			  double [][] yPrPb=new double [3][];
			  //			if (dct_parameters.color_DCT){
			  // need to get YPbPr - not RGB here
			  //			} else {
			  correctionColorProc.processColorsWeights(stack, // just gamma convert? TODO: Cleanup? Convert directly form the linear YPrPb
					  //					  255.0/this.psfSubpixelShouldBe4/this.psfSubpixelShouldBe4, //  double scale,     // initial maximal pixel value (16))
					  //					  255.0/eyesisCorrections.psfSubpixelShouldBe4/eyesisCorrections.psfSubpixelShouldBe4/scaleExposure, //  double scale,     // initial maximal pixel value (16))
					  //					  255.0/2/2/scaleExposure, //  double scale,     // initial maximal pixel value (16))
					  255.0/scaleExposures[iQuad], //  double scale,     // initial maximal pixel value (16))
					  colorProcParameters,
					  channelGainParameters,
					  channel,
					  null, //correctionDenoise.getDenoiseMask(),
					  this.correctionsParameters.blueProc,
					  debugLevel);
			  if (debugLevel > 1) System.out.println("Processed colors to YPbPr, total number of slices="+stack.getSize());
			  if (debugLevel > 1) {
				  ImagePlus imp_dbg=new ImagePlus("procColors",stack);
				  eyesisCorrections.saveAndShow(
						  imp_dbg,
						  this.correctionsParameters);
			  }
			  float [] fpixels;
			  int [] slices_YPrPb = {8,6,7};
			  yPrPb=new double [3][];
			  for (int n = 0; n < slices_YPrPb.length; n++){
				  fpixels = (float[]) stack.getPixels(slices_YPrPb[n]);
				  yPrPb[n] = new double [fpixels.length];
				  for (int i = 0; i < fpixels.length; i++) yPrPb[n][i] = fpixels[i];
			  }

			  if (toRGB) {
				  if (debugLevel > 0){
					  System.out.println("correctionColorProc.YPrPbToRGB");
				  }
				  stack =  YPrPbToRGB(yPrPb,
						  colorProcParameters.kr,        // 0.299;
						  colorProcParameters.kb,        // 0.114;
						  stack.getWidth());

				  title=titleFull; // including "-DECONV" or "-COMBO"
				  titleFull=title+"-RGB-float";
				  //Trim stack to just first 3 slices
				  if (debugLevel > 1){ // 2){
					  ImagePlus imp_dbg=new ImagePlus("YPrPbToRGB",stack);
					  eyesisCorrections.saveAndShow(
							  imp_dbg,
							  this.correctionsParameters);
				  }
				  while (stack.getSize() > 3) stack.deleteLastSlice();
				  if (debugLevel > 1) System.out.println("Trimming color stack");
			  } else {
				  title=titleFull; // including "-DECONV" or "-COMBO"
				  titleFull=title+"-YPrPb"; // including "-DECONV" or "-COMBO"
				  if (debugLevel > 1) System.out.println("Using full stack, including YPbPr");
			  }

			  results[iQuad]= new ImagePlus(titleFull, stack);
			  // rotate the result
			  if (rotate){ // never rotate for equirectangular
				  stack=eyesisCorrections.rotateStack32CW(stack);
			  }
			  if (!toRGB && !this.correctionsParameters.jpeg){ // toRGB set for equirectangular
				  if (debugLevel > 1) System.out.println("!toRGB && !this.correctionsParameters.jpeg");
				  eyesisCorrections.saveAndShow(results[iQuad], this.correctionsParameters);
				  continue; // return result;
			  } else { // that's not the end result, save if required
				  if (debugLevel > 1) System.out.println("!toRGB && !this.correctionsParameters.jpeg - else");
				  eyesisCorrections.saveAndShow(results[iQuad],
						  eyesisCorrections.correctionsParameters,
						  eyesisCorrections.correctionsParameters.save32,
						  false,
						  eyesisCorrections.correctionsParameters.JPEG_quality); // save, no show
			  }
			  // convert to RGB48 (16 bits per color component)
			  ImagePlus imp_RGB;
			  stack=eyesisCorrections.convertRGB32toRGB16Stack(
					  stack,
					  rgbParameters);

			  titleFull=title+"-RGB48";
			  results[iQuad]= new ImagePlus(titleFull, stack);

			  //			  ImagePlus imp_RGB24;
			  results[iQuad].updateAndDraw();
			  if (debugLevel > 1) System.out.println("result.updateAndDraw(), "+titleFull+"-RGB48");

			  CompositeImage compositeImage=eyesisCorrections.convertToComposite(results[iQuad]);

			  if (!this.correctionsParameters.jpeg && !advanced){ // RGB48 was the end result
				  if (debugLevel > 1) System.out.println("if (!this.correctionsParameters.jpeg && !advanced)");
				  eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters);
				  continue; // return result;
			  } else { // that's not the end result, save if required
				  if (debugLevel > 1) System.out.println("if (!this.correctionsParameters.jpeg && !advanced) - else");
				  eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters, this.correctionsParameters.save16, false); // save, no show
				  //				  eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters, this.correctionsParameters.save16, true); // save, no show
			  }

			  imp_RGB=eyesisCorrections.convertRGB48toRGB24(
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
			  eyesisCorrections.saveAndShow(imp_RGB, this.correctionsParameters);
		  }
		  return results;
	  }

	  class SetChannels{
		  String set_name;    // set name (timestamp)
		  int [] file_number; // array of file numbers for channels
		  public SetChannels(String name, int[] fn){
			  set_name = name;
			  file_number = fn;
		  }
		  public String name() {
			  return set_name;
		  }
		  public int [] fileNumber() {
			  return file_number;
		  }
		  public int fileNumber(int i) {
			  return file_number[i];
		  }
	  }

	  SetChannels [] setChannels(
			  int debugLevel) {
		  return setChannels(null, debugLevel);
	  }



	  SetChannels [] setChannels(
			  String single_set_name, // process only files that contain specified series (timestamp) in the name
			  int debugLevel) {
		  String [] sourceFiles=correctionsParameters.getSourcePaths();
		  boolean [] enabledFiles=new boolean[sourceFiles.length];
		  for (int i=0;i<enabledFiles.length;i++) enabledFiles[i]=false;
		  int numFilesToProcess=0;
		  int numImagesToProcess=0;
		  for (int nFile=0;nFile<enabledFiles.length;nFile++){
			  if ((sourceFiles[nFile]!=null) && (sourceFiles[nFile].length()>1) && ((single_set_name == null) || (sourceFiles[nFile].contains(single_set_name)))) {
				  int [] channels={correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile])};
				  if (correctionsParameters.isJP4()){
					  int subCamera= channels[0]- correctionsParameters.firstSubCamera; // to match those in the sensor files
					  // removeUnusedSensorData should be off!?
					  channels=this.eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera); // limit here or disable Error
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (eyesisCorrections.isChannelEnabled(channels[i])){
						  if (!enabledFiles[nFile]) numFilesToProcess++;
						  enabledFiles[nFile]=true;
						  numImagesToProcess++;
					  }
				  }
			  }
		  }
		  if (numFilesToProcess==0){
			  System.out.println("No files to process (of "+sourceFiles.length+")");
			  return null;
		  } else {
			  if (debugLevel>0) System.out.println(numFilesToProcess+ " files to process (of "+sourceFiles.length+"), "+numImagesToProcess+" images to process");
		  }
//		  double [] referenceExposures=eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		  int [][] fileIndices=new int [numImagesToProcess][2]; // file index, channel number
		  int index=0;
		  for (int nFile=0;nFile<enabledFiles.length;nFile++){ // enabledFiles not used anymore?
			  if ((sourceFiles[nFile]!=null) && (sourceFiles[nFile].length()>1) && ((single_set_name == null) || (sourceFiles[nFile].contains(single_set_name)))) {
				  int [] channels={correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile])};
				  if (correctionsParameters.isJP4()){
					  int subCamera= channels[0]- correctionsParameters.firstSubCamera; // to match those in the sensor files
					  channels=eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (eyesisCorrections.isChannelEnabled(channels[i])){
						  fileIndices[index  ][0]=nFile;
						  fileIndices[index++][1]=channels[i];
					  }
				  }
			  }
		  }
		  ArrayList<String> setNames = new ArrayList<String>();
		  ArrayList<ArrayList<Integer>> setFiles = new ArrayList<ArrayList<Integer>>();

		  for (int iImage=0;iImage<fileIndices.length;iImage++){
			  int nFile=fileIndices[iImage][0];
			  String setName = correctionsParameters.getNameFromSourceTiff(sourceFiles[nFile]);
			  if (!setNames.contains(setName)) {
				  setNames.add(setName);
				  setFiles.add(new ArrayList<Integer>());
			  }
			  setFiles.get(setNames.indexOf(setName)).add(new Integer(iImage));
		  }
		  SetChannels [] sc = new SetChannels[setNames.size()];
		  for (int nSet = 0; nSet < setNames.size(); nSet++){
			  int maxChn = 0;
			  for (int i = 0; i < setFiles.get(nSet).size(); i++){
				  int chn = fileIndices[setFiles.get(nSet).get(i)][1];
				  if (chn > maxChn) maxChn = chn;
			  }
			  int [] channelFiles = new int[maxChn+1];
			  for (int i =0; i < channelFiles.length; i++) channelFiles[i] = -1;
			  for (int i = 0; i < setFiles.get(nSet).size(); i++){
				  channelFiles[fileIndices[setFiles.get(nSet).get(i)][1]] = fileIndices[setFiles.get(nSet).get(i)][0];
			  }
			  sc[nSet] = new SetChannels(setNames.get(nSet), channelFiles);
		  }
		  return sc;
	  }


	  int getTotalFiles(SetChannels [] sc) {
		  int nf = 0;
		  for (int i = 0; i < sc.length; i++) nf+=sc[i].fileNumber().length;
		  return nf;
	  }

	  /**
	   * Pad acquired Bayer image to the full sensor width/height. Used when optical center pixel coordinates do not match for channels
	   * and WOI is adjusted during image capture to avoid ERS mismatch between horizontal pairs
	   * @param imp_src source image with WOI specified as properties (sizes and offsets should be even)
	   * @param wh {sesnor_width, sensor_height} in pixels
	   * @param replicate fill gaps by replicating existing pixels
	   * @return full size image
	   */

	  ImagePlus padBayerToFullSize(
			  ImagePlus imp_src,
			  int [] wh,
			  boolean replicate) {
		  int woi_top =    Integer.parseInt((String) imp_src.getProperty("WOI_TOP")); // enforce even
		  int woi_left =   Integer.parseInt((String) imp_src.getProperty("WOI_LEFT"));
		  int woi_width =  imp_src.getWidth(); // Integer.parseInt((String) imp_src.getProperty("WOI_WIDTH"));
		  int woi_height = imp_src.getHeight(); // Integer.parseInt((String) imp_src.getProperty("WOI_HEIGHT"));
		  Properties properties = imp_src.getProperties();


		  if ((woi_top == 0) && (woi_left == 0) && (woi_width == wh[0])  && (woi_height == wh[1])){
			  return imp_src; // good as is
		  }
		  float [] full_pixels = new float [wh[0]*wh[1]];
		  float [] pixels=(float []) imp_src.getProcessor().getPixels();
		  int dst_col = woi_left;
		  int copy_width = woi_width;
		  if ((dst_col + copy_width) > wh[0]) {
			  copy_width = wh[0] - dst_col;
		  }
		  for (int src_row = 0; src_row < woi_height; src_row++) {
			  int dst_row = src_row + woi_top;
			  if (dst_row < wh[1]) {
				  System.arraycopy( pixels,   src_row * woi_width,  full_pixels, dst_row * wh[0] + dst_col,  copy_width);
			  }
		  }
		  if (replicate) {
			  // replicate top
			  for (int dst_row = 0; dst_row < woi_top; dst_row++) {
				  int src_row = woi_top + (dst_row & 1);
				  System.arraycopy( full_pixels,   src_row * wh[0] + dst_col,  full_pixels, dst_row * wh[0] + dst_col,  copy_width);
			  }
			  // replicate bottom
			  for (int dst_row = woi_top + woi_height; dst_row < wh[1]; dst_row++) {
				  int src_row = woi_top + woi_height - 2 + (dst_row & 1);
				  System.arraycopy( full_pixels,   src_row * wh[0] + dst_col,  full_pixels, dst_row * wh[0] + dst_col,  copy_width);
			  }
			  // right and left are not likely, as there is no need to use them - horizontal mismatch does not influence ERS
			  for (int col = 0; col < woi_left; col++) {
				  for (int row = 0; row < wh[1]; row++) {
					  full_pixels[row*wh[0] + col] = full_pixels[row*wh[0] + woi_left + (col & 1)];
				  }
			  }

			  for (int col = woi_left + woi_width; col < wh[0]; col++) {
				  for (int row = 0; row < wh[1]; row++) {
					  full_pixels[row*wh[0] + col] = full_pixels[row*wh[0] + woi_left + woi_width - 2 +(col & 1)];
				  }
			  }
		  }
		  ImageProcessor ip = new FloatProcessor(wh[0],wh[1]);
		  ip.setPixels(full_pixels);
		  ip.resetMinAndMax(); // is it needed here?
		  ImagePlus imp = new ImagePlus(imp_src.getTitle(),ip); // OK to have the same name?
		  for (Map.Entry<?, ?> entry: properties.entrySet()) {
			  String key = (String) entry.getKey();
			  String value = (String) entry.getValue();
			  imp.setProperty(key, value);
		  }
		  imp.setProperty("WOI_WIDTH", wh[0]+"");
		  imp.setProperty("WOI_HEIGHTH", wh[1]+"");
		  imp.setProperty("WOI_TOP", "0");
		  imp.setProperty("WOI_LEFT", "0");
		  return imp;
	  }

	  /**
	   * Conditions images for a single image set
	   * @param clt_parameters      various parameters
	   * @param sourceFiles         array of source file paths matching indices in channelFiles
	   * @param set_name            name of the current image set (normally timestamp with "_" for decimal point
	   * @param referenceExposures  array of per-channel reference exposures, data will be scaled using Exif exposure of each file
	   * @param channelFiles        array of file indices (in sourceFiles array) for the camera channels ([0] - index of the first channel file)
	   * @param scaleExposures      array of per-channel "brightening" of images (reference exposure/ actual exposure
	   * @param saturation_imp      per-channel bitmask of the saturated pixels or null. Should be initialized by the caller, will be filled here
	   * @param debugLevel          debug (verbosity) level
	   * @return array of per-channel ImagePlus objects to process (with saturation_imp)
	   */
	  public ImagePlus[] conditionImageSet(
			  EyesisCorrectionParameters.CLTParameters  clt_parameters,
			  String []                                 sourceFiles,
			  String                                    set_name,
			  double []                                 referenceExposures,
			  int []                                    channelFiles,
			  double []                                 scaleExposures,
			  boolean [][]                              saturation_imp,
			  int                                       debugLevel)
	  {
		  ImagePlus [] imp_srcs = new ImagePlus[channelFiles.length];
//		  double [] scaleExposures = new double[channelFiles.length]; //
		  double [][] dbg_dpixels = new double [channelFiles.length][];
//		  int [] fullWindowWH = geometryCorrection.getSensorWH();

		  for (int srcChannel=0; srcChannel < channelFiles.length; srcChannel++){
			  int nFile=channelFiles[srcChannel]; // channelFiles[srcChannel];

			  imp_srcs[srcChannel]=null;
			  if (nFile >=0){
				  if (correctionsParameters.isJP4()){
					  int subchannel=eyesisCorrections.pixelMapping.getSubChannel(srcChannel);
					  if (this.correctionsParameters.swapSubchannels01) {
						  switch (subchannel){
						  case 0: subchannel=1; break;
						  case 1: subchannel=0; break;
						  }
					  }
					  if (debugLevel>0) System.out.println("Processing set " + set_name+" channel "+srcChannel+" - subchannel "+subchannel+" of "+sourceFiles[nFile]);
					  ImagePlus imp_composite=eyesisCorrections.JP4_INSTANCE.open(
							  "", // path,
							  sourceFiles[nFile],
							  "",  //arg - not used in JP46 reader
							  true, // un-apply camera color gains
							  null, // new window
							  false); // do not show
					  imp_srcs[srcChannel]=eyesisCorrections.JP4_INSTANCE.demuxImage(imp_composite, subchannel);
					  if (imp_srcs[srcChannel] == null) imp_srcs[srcChannel] = imp_composite; // not a composite image
					  // do we need to add any properties?
				  } else {
					  imp_srcs[srcChannel]=new ImagePlus(sourceFiles[nFile]);
					  //					  (new JP46_Reader_camera(false)).decodeProperiesFromInfo(imp_src); // decode existent properties from info
					  eyesisCorrections.JP4_INSTANCE.decodeProperiesFromInfo(imp_srcs[srcChannel]); // decode existent properties from info
					  if (debugLevel>0) System.out.println("Processing "+sourceFiles[nFile]);
				  }
// imp_srcs[srcChannel].show(); // REMOVE ME!

				  imp_srcs[srcChannel] =  padBayerToFullSize(
						  imp_srcs[srcChannel], // ImagePlus imp_src,
						  eyesisCorrections.pixelMapping.sensors[srcChannel].getSensorWH(),
						  true); // boolean replicate);
				  scaleExposures[srcChannel] = 1.0;
				  if (!Double.isNaN(referenceExposures[nFile]) && (imp_srcs[srcChannel].getProperty("EXPOSURE")!=null)){
					  scaleExposures[srcChannel] = referenceExposures[nFile]/Double.parseDouble((String) imp_srcs[srcChannel].getProperty("EXPOSURE"));
					  if (debugLevel > -1) {
						  System.out.println("Will scale intensity (to compensate for exposure) by  "+scaleExposures[srcChannel]+
								  ", EXPOSURE = "+imp_srcs[srcChannel].getProperty("EXPOSURE"));
					  }
				  }
				  imp_srcs[srcChannel].setProperty("name",    correctionsParameters.getNameFromSourceTiff(sourceFiles[nFile]));
				  imp_srcs[srcChannel].setProperty("channel", srcChannel); // it may already have channel
				  imp_srcs[srcChannel].setProperty("path",    sourceFiles[nFile]); // it may already have channel

				  if (this.correctionsParameters.pixelDefects && (eyesisCorrections.defectsXY!=null)&& (eyesisCorrections.defectsXY[srcChannel]!=null)){
					  // apply pixel correction
					  int numApplied=	eyesisCorrections.correctDefects(
							  imp_srcs[srcChannel],
							  srcChannel,
							  debugLevel);
					  if ((debugLevel>0) && (numApplied>0)) { // reduce verbosity after verified defect correction works
						  System.out.println("Corrected "+numApplied+" pixels in "+sourceFiles[nFile]);
					  }
				  }
				  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
				  int width =  imp_srcs[srcChannel].getWidth();
				  int height = imp_srcs[srcChannel].getHeight();
				  if (debugLevel > -1) {
					  double [] max_pix= {0.0, 0.0, 0.0, 0.0};
//					  for (int y = 0; y < height-1; y+=2){
					  for (int y = 0; y < 499; y+=2){
//						  for (int x = 0; x < width-1; x+=2){
						  for (int x = width/2; x < width-1; x+=2){
							  if (pixels[y*width+x        ] > max_pix[0])  max_pix[0] = pixels[y*width+x        ];
							  if (pixels[y*width+x+      1] > max_pix[1])  max_pix[1] = pixels[y*width+x+      1];
							  if (pixels[y*width+x+width  ] > max_pix[2])  max_pix[2] = pixels[y*width+x+width  ];
							  if (pixels[y*width+x+width+1] > max_pix[3])  max_pix[3] = pixels[y*width+x+width+1];
						  }
					  }
					  System.out.println(String.format("channel %d max_pix[] = %6.2f %6.2f %6.2f %6.2f", srcChannel, max_pix[0], max_pix[1], max_pix[2], max_pix[3]));
					  dbg_dpixels[srcChannel] = new double [pixels.length];
					  for (int i = 0; i < pixels.length; i++) dbg_dpixels[srcChannel][i] = pixels[i];
					  //						  imp_srcs[srcChannel].show();
				  }
				  if (clt_parameters.sat_level > 0.0){
					  double [] saturations = {
							  Double.parseDouble((String) imp_srcs[srcChannel].getProperty("saturation_1")),
							  Double.parseDouble((String) imp_srcs[srcChannel].getProperty("saturation_0")),
							  Double.parseDouble((String) imp_srcs[srcChannel].getProperty("saturation_3")),
							  Double.parseDouble((String) imp_srcs[srcChannel].getProperty("saturation_2"))};
					  saturation_imp[srcChannel] = new boolean[width*height];
					  System.out.println(String.format("channel %d saturations = %6.2f %6.2f %6.2f %6.2f", srcChannel,
							  saturations[0],saturations[1],saturations[2],saturations[3]));
					  double [] scaled_saturations = new double [saturations.length];
					  for (int i = 0; i < scaled_saturations.length; i++){
						  scaled_saturations[i] = saturations[i] * clt_parameters.sat_level;
					  }
					  for (int y = 0; y < height-1; y+=2){
						  for (int x = 0; x < width-1; x+=2){
							  if (pixels[y*width+x        ] > scaled_saturations[0])  saturation_imp[srcChannel][y*width+x        ] = true;
							  if (pixels[y*width+x+      1] > scaled_saturations[1])  saturation_imp[srcChannel][y*width+x      +1] = true;
							  if (pixels[y*width+x+width  ] > scaled_saturations[2])  saturation_imp[srcChannel][y*width+x+width  ] = true;
							  if (pixels[y*width+x+width+1] > scaled_saturations[3])  saturation_imp[srcChannel][y*width+x+width+1] = true;
						  }
					  }
				  }


				  if (this.correctionsParameters.vignetting){
					  if ((eyesisCorrections.channelVignettingCorrection==null) || (srcChannel<0) || (srcChannel>=eyesisCorrections.channelVignettingCorrection.length) || (eyesisCorrections.channelVignettingCorrection[srcChannel]==null)){
						  System.out.println("No vignetting data for channel "+srcChannel);
						  return null;
					  }
					  ///						  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();


					  if (pixels.length!=eyesisCorrections.channelVignettingCorrection[srcChannel].length){
						  System.out.println("Vignetting data for channel "+srcChannel+" has "+eyesisCorrections.channelVignettingCorrection[srcChannel].length+" pixels, image "+sourceFiles[nFile]+" has "+pixels.length);
						  return null;
					  }
					  // TODO: Move to do it once:
					  double min_non_zero = 0.0;
					  for (int i=0;i<pixels.length;i++){
						  double d = eyesisCorrections.channelVignettingCorrection[srcChannel][i];
						  if ((d > 0.0) && ((min_non_zero == 0) || (min_non_zero > d))){
							  min_non_zero = d;
						  }
					  }
					  double max_vign_corr = clt_parameters.vignetting_range*min_non_zero;

					  System.out.println("Vignetting data: channel="+srcChannel+", min = "+min_non_zero);
					  for (int i=0;i<pixels.length;i++){
						  double d = eyesisCorrections.channelVignettingCorrection[srcChannel][i];
						  if (d > max_vign_corr) d = max_vign_corr;
						  pixels[i]*=d;
					  }
					  // Scale here, combine with vignetting later?
					  ///						  int width =  imp_srcs[srcChannel].getWidth();
					  ///						  int height = imp_srcs[srcChannel].getHeight();
					  for (int y = 0; y < height-1; y+=2){
						  for (int x = 0; x < width-1; x+=2){
							  pixels[y*width+x        ] *= clt_parameters.scale_g;
							  pixels[y*width+x+width+1] *= clt_parameters.scale_g;
							  pixels[y*width+x      +1] *= clt_parameters.scale_r;
							  pixels[y*width+x+width  ] *= clt_parameters.scale_b;
						  }
					  }

				  } else { // assuming GR/BG pattern
					  System.out.println("Applying fixed color gain correction parameters: Gr="+
							  clt_parameters.novignetting_r+", Gg="+clt_parameters.novignetting_g+", Gb="+clt_parameters.novignetting_b);
					  ///						  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
					  ///						  int width =  imp_srcs[srcChannel].getWidth();
					  ///						  int height = imp_srcs[srcChannel].getHeight();
					  double kr = clt_parameters.scale_r/clt_parameters.novignetting_r;
					  double kg = clt_parameters.scale_g/clt_parameters.novignetting_g;
					  double kb = clt_parameters.scale_b/clt_parameters.novignetting_b;
					  for (int y = 0; y < height-1; y+=2){
						  for (int x = 0; x < width-1; x+=2){
							  pixels[y*width+x        ] *= kg;
							  pixels[y*width+x+width+1] *= kg;
							  pixels[y*width+x      +1] *= kr;
							  pixels[y*width+x+width  ] *= kb;
						  }
					  }
				  }
			  }
		  }
		  // temporary applying scaleExposures[srcChannel] here, setting it to all 1.0
		  System.out.println("Temporarily applying scaleExposures[] here" );
		  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
			  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
			  for (int i = 0; i < pixels.length; i++){
				  pixels[i] *= scaleExposures[srcChannel];
			  }
			  scaleExposures[srcChannel] = 1.0;
		  }

		  if ((debugLevel > -1) && (saturation_imp != null)){
			  String [] titles = {"chn0","chn1","chn2","chn3"};
			  double [][] dbg_satur = new double [saturation_imp.length] [saturation_imp[0].length];
			  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
				  for (int i = 0; i < saturation_imp[srcChannel].length; i++){
					  dbg_satur[srcChannel][i] = saturation_imp[srcChannel][i]? 1.0 : 0.0;
				  }
			  }
			  int width =  imp_srcs[0].getWidth();
			  int height = imp_srcs[0].getHeight();
			  (new showDoubleFloatArrays()).showArrays(dbg_satur, width, height, true, "Saturated" , titles);

			  if (debugLevel > -1) { // 0){
				  double [][] dbg_dpixels_norm = new double [channelFiles.length][];
				  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
					  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
					  dbg_dpixels_norm[srcChannel] = new double[pixels.length];
					  for (int i = 0; i < pixels.length; i++){
						  dbg_dpixels_norm[srcChannel][i] = pixels[i];
					  }
				  }
				  (new showDoubleFloatArrays()).showArrays(dbg_dpixels, width, height, true, "dpixels" , titles);
				  (new showDoubleFloatArrays()).showArrays(dbg_dpixels_norm, width, height, true, "dpixels_norm" , titles);
				  double [][] dbg_dpixels_split = new double [4 * dbg_dpixels.length][dbg_dpixels[0].length / 4];
				  String [] dbg_titles = {"g1_0","r_0","b_0","g2_0","g1_2","r_1","b_1","g2_1","g1_2","r_2","b_2","g2_2","g1_3","r_3","b_3","g2_3"};
				  for (int srcChn = 0; srcChn < 4; srcChn++) {
					  for (int y = 0; y < height-1; y+=2){
						  for (int x = 0; x < width-1; x+=2){
							  dbg_dpixels_split[ 0 + 4 * srcChn][ y*width/4 +x/2 ] = dbg_dpixels_norm[srcChn][y * width + x             ];
							  dbg_dpixels_split[ 3 + 4 * srcChn][ y*width/4 +x/2 ] = dbg_dpixels_norm[srcChn][y * width + x  + width + 1];
							  dbg_dpixels_split[ 1 + 4 * srcChn][ y*width/4 +x/2 ] = dbg_dpixels_norm[srcChn][y * width + x  +         1];
							  dbg_dpixels_split[ 2 + 4 * srcChn][ y*width/4 +x/2 ] = dbg_dpixels_norm[srcChn][y * width + x  + width    ];
						  }
					  }
				  }
				  (new showDoubleFloatArrays()).showArrays(dbg_dpixels_split, width/2, height/2, true, "dpixels_split" , dbg_titles);
			  }

		  }
		  //			  Overlay ovl = imp_srcs[0].getOverlay();
		  // once per quad here
		  // may need to equalize gains between channels
		  if (clt_parameters.gain_equalize || clt_parameters.colors_equalize){
			  channelGainsEqualize(
					  clt_parameters.gain_equalize,
					  clt_parameters.colors_equalize,
					  clt_parameters.nosat_equalize, // boolean nosat_equalize,
					  channelFiles,
					  imp_srcs,
					  saturation_imp, // boolean[][] saturated,
					  set_name,       // setNames.get(nSet), // just for debug messages == setNames.get(nSet)
					  debugLevel);
		  }
		  return imp_srcs;
	  }


	  public void processCLTQuadCorrs(
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
//			  int          convolveFFTSize, // 128 - fft size, kernel size should be size/2
			  final boolean    apply_corr, // calculate and apply additional fine geometry correction
			  final boolean    infinity_corr, // calculate and apply geometry correction at infinity
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel)
	  {
		  if (infinity_corr && (clt_parameters.z_correction != 0.0)){
			  System.out.println(
					  "****************************************\n"+
					  "* Resetting manual infinity correction *\n"+
					  "****************************************\n");
			  clt_parameters.z_correction = 0.0;
		  }

		  this.startTime=System.nanoTime();
		  String [] sourceFiles=correctionsParameters.getSourcePaths();
		  SetChannels [] set_channels=setChannels(debugLevel);
		  if ((set_channels == null) || (set_channels.length==0)) {
			  System.out.println("No files to process (of "+sourceFiles.length+")");
			  return;
		  }
		  double [] referenceExposures=eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		  for (int nSet = 0; nSet < set_channels.length; nSet++){
			  int [] channelFiles = set_channels[nSet].fileNumber();
			  boolean [][] saturation_imp = (clt_parameters.sat_level > 0.0)? new boolean[channelFiles.length][] : null;
			  double [] scaleExposures = new double[channelFiles.length];

			  ImagePlus [] imp_srcs = conditionImageSet(
					  clt_parameters,             // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					  sourceFiles,                // String []                                 sourceFiles,
					  set_channels[nSet].name(),  // String                                    set_name,
					  referenceExposures,         // double []                                 referenceExposures,
					  channelFiles,               // int []                                    channelFiles,
					  scaleExposures,   //output  // double [] scaleExposures
					  saturation_imp,   //output  // boolean [][]                              saturation_imp,
					  debugLevel); // int                                       debugLevel);


			  // once per quad here
			  processCLTQuadCorr( // returns ImagePlus, but it already should be saved/shown
					  imp_srcs, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
					  saturation_imp, // boolean [][] saturation_imp, // (near) saturated pixels or null
					  clt_parameters,
					  debayerParameters,
					  colorProcParameters,
					  channelGainParameters,
					  rgbParameters,
//					  convolveFFTSize, // 128 - fft size, kernel size should be size/2
					  scaleExposures,
					  apply_corr, // calculate and apply additional fine geometry correction
					  infinity_corr, // calculate and apply geometry correction at infinity
					  threadsMax,  // maximal number of threads to launch
					  updateStatus,
					  debugLevel);
			  Runtime.getRuntime().gc();
			  if (debugLevel >-1) System.out.println("Processing set "+(nSet+1)+" (of "+set_channels.length+") finished at "+
					  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
			  if (eyesisCorrections.stopRequested.get()>0) {
				  System.out.println("User requested stop");
				  System.out.println("Processing "+(nSet + 1)+" file sets (of "+set_channels.length+") finished at "+
						  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
				  return;
			  }
		  }
		  System.out.println("processCLTQuadCorrs(): processing "+getTotalFiles(set_channels)+" files ("+set_channels.length+" file sets) finished at "+
				  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
	  }


	  public void channelGainsEqualize(
			  boolean gain_equalize,
			  boolean colors_equalize,
			  boolean nosat_equalize,
			  int [] channelFiles,
			  ImagePlus [] imp_srcs,
			  boolean[][] saturated,
			  String setName, // just for debug messages == setNames.get(nSet)
			  int debugLevel){
		  boolean use_new = true; // false;
		  if (use_new){
			  channelGainsEqualize_new(
					  gain_equalize,
					  colors_equalize,
					  nosat_equalize,
					  channelFiles,
					  imp_srcs,
					  saturated,
					  setName, // just for debug messages == setNames.get(nSet)
					  debugLevel);
		  } else {
			  channelGainsEqualize_old(
					  gain_equalize,
					  colors_equalize,
					  nosat_equalize,
					  channelFiles,
					  imp_srcs,
					  saturated,
					  setName, // just for debug messages == setNames.get(nSet)
					  debugLevel);
		  }
	  }
	  public void channelGainsEqualize_old(
			  boolean gain_equalize,
			  boolean colors_equalize,
			  boolean nosat_equalize,
			  int [] channelFiles,
			  ImagePlus [] imp_srcs,
			  boolean[][] saturated,
			  String setName, // just for debug messages == setNames.get(nSet)
			  int debugLevel){
		  double [][] avr_pix = new double [channelFiles.length][3];
		  double [] avr_RGB = {0.0,0.0,0.0};
		  int numChn = 0;
		  for (int srcChannel=0; srcChannel < channelFiles.length; srcChannel++){
			  int nFile=channelFiles[srcChannel];
			  if (nFile >=0){
				  for (int i = 0; i < avr_pix[srcChannel].length; i++) avr_pix[srcChannel][i] = 0;
//				  int [] num_nonsat = {0,0,0};
				  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
				  int width =  imp_srcs[srcChannel].getWidth();
				  int height = imp_srcs[srcChannel].getHeight();
				  for (int y = 0; y < height-1; y+=2){
					  for (int x = 0; x < width-1; x+=2){
						  avr_pix[srcChannel][0] += pixels[y*width+x      +1];
						  avr_pix[srcChannel][2] += pixels[y*width+x+width  ];
						  avr_pix[srcChannel][1] += pixels[y*width+x        ];
						  avr_pix[srcChannel][1] += pixels[y*width+x+width+1];
					  }
				  }
				  avr_pix[srcChannel][0] /= 0.25*width*height;
				  avr_pix[srcChannel][1] /= 0.5 *width*height;
				  avr_pix[srcChannel][2] /= 0.25*width*height;
				  for (int j=0; j < avr_RGB.length; j++) avr_RGB[j] += avr_pix[srcChannel][j];
				  numChn++;
				  if (debugLevel > -2) {
					  System.out.println("processCLTSets(): set "+ setName + " channel "+srcChannel+
							  " R"+avr_pix[srcChannel][0]+" G"+avr_pix[srcChannel][1]+" B"+avr_pix[srcChannel][2]);
				  }

			  }
		  }
		  for (int j=0; j < avr_RGB.length; j++) avr_RGB[j] /= numChn;
		  if (debugLevel > -2) {
			  System.out.println("processCLTSets(): set "+ setName + "average color values: "+
					  " R="+avr_RGB[0]+" G=" + avr_RGB[1]+" B=" + avr_RGB[2]);
		  }
		  for (int srcChannel=0; srcChannel < channelFiles.length; srcChannel++){
			  int nFile=channelFiles[srcChannel];
			  if (nFile >=0){
				  double [] scales = new double [avr_RGB.length];
				  for (int j=0;j < scales.length; j++){
					  scales[j] = 1.0;
					  if (gain_equalize){
						  scales[j] *=  avr_RGB[1]/avr_pix[srcChannel][1]; // 1 - index of green color
					  }
					  if (colors_equalize){
						  scales[j] *=  avr_RGB[j]/avr_pix[srcChannel][j] / (avr_RGB[1]/avr_pix[srcChannel][1]);
					  }
				  }
				  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
				  System.out.println("Channel "+srcChannel+ ":  scales[] = "+scales[0]+", "+scales[1]+", "+scales[2]);
				  int width =  imp_srcs[srcChannel].getWidth();
				  int height = imp_srcs[srcChannel].getHeight();
				  for (int y = 0; y < height-1; y+=2){
					  for (int x = 0; x < width-1; x+=2){
						  pixels[y*width+x        ] *= scales[1];
						  pixels[y*width+x+width+1] *= scales[1];
						  pixels[y*width+x      +1] *= scales[0];
						  pixels[y*width+x+width  ] *= scales[2];
					  }
				  }
			  }
		  }
	  }

	  public void channelGainsEqualize_new(
			  boolean gain_equalize,
			  boolean colors_equalize,
			  boolean nosat_equalize,
			  int [] channelFiles,
			  ImagePlus [] imp_srcs,
			  boolean[][] saturated,
			  String setName, // just for debug messages == setNames.get(nSet)
			  int debugLevel){
		  double [][] avr_pix = new double [channelFiles.length][4];
		  double [] avr_RGGB = {0.0,0.0,0.0,0.0};
		  int numChn = 0;
		  for (int srcChannel=0; srcChannel < channelFiles.length; srcChannel++){
			  int nFile=channelFiles[srcChannel];
			  if (nFile >=0){
				  for (int i = 0; i < avr_pix[srcChannel].length; i++) avr_pix[srcChannel][i] = 0;
				  int [] num_nonsat = {0,0,0,0};
				  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
				  int width =  imp_srcs[srcChannel].getWidth();
				  int height = imp_srcs[srcChannel].getHeight();
				  int indx;
				  for (int y = 0; y < height-1; y+=2){
					  for (int x = 0; x < width-1; x+=2){
						  indx = y * width + (x +1);  // R
						  if (!nosat_equalize || (saturated == null) || !saturated[srcChannel][indx]) {
							  avr_pix[srcChannel][0] += pixels[indx];
							  num_nonsat[0]++;
						  }

						  indx = (y + 1) * width + x;  // B
						  if (!nosat_equalize || (saturated == null) || !saturated[srcChannel][indx]) {
							  avr_pix[srcChannel][3] += pixels[indx];
							  num_nonsat[3]++;
						  }
						  indx = y * width + x;  // G1
						  if (!nosat_equalize || (saturated == null) || !saturated[srcChannel][indx]) {
							  avr_pix[srcChannel][1] += pixels[indx];
							  num_nonsat[1]++;
						  }
						  indx = (y + 1) * width + (x + 1);  // G2
						  if (!nosat_equalize || (saturated == null) || !saturated[srcChannel][indx]) {
							  avr_pix[srcChannel][2] += pixels[indx];
							  num_nonsat[2]++;
						  }

//						  avr_pix[srcChannel][0] += pixels[y*width+x      +1];
//						  avr_pix[srcChannel][2] += pixels[y*width+x+width  ];
//						  avr_pix[srcChannel][1] += pixels[y*width+x        ];
//						  avr_pix[srcChannel][1] += pixels[y*width+x+width+1];
					  }
				  }
				  for (int i = 0; i < num_nonsat.length; i++){
					  avr_pix[srcChannel][i] /= num_nonsat[i];
				  }
				  for (int j=0; j < avr_RGGB.length; j++) avr_RGGB[j] += avr_pix[srcChannel][j];
				  numChn++;
				  if (debugLevel > -2) {
					  System.out.println("processCLTSets(): set "+ setName + " channel "+srcChannel+
							  " R="+avr_pix[srcChannel][0]+" G1="+avr_pix[srcChannel][1]+" G2="+avr_pix[srcChannel][2]+" B="+avr_pix[srcChannel][3]);
				  }

			  }
		  }
		  for (int j=0; j < avr_RGGB.length; j++) avr_RGGB[j] /= numChn;
		  double avr_G = 0.5 * (avr_RGGB[1] + avr_RGGB[2]);
		  if (debugLevel > -2) {
			  System.out.println("processCLTSets(): set "+ setName + " average color values: "+
					  " R="+avr_RGGB[0]+" G=" + avr_G+" ("+avr_RGGB[1]+", "+avr_RGGB[2]+") B=" + avr_RGGB[3]);
		  }
		  for (int srcChannel=0; srcChannel < channelFiles.length; srcChannel++){
			  int nFile=channelFiles[srcChannel];
			  if (nFile >=0){
				  double [] scales = new double [avr_RGGB.length];
				  double avr_g = 0.5 * (avr_pix[srcChannel][1] + avr_pix[srcChannel][2]);
				  for (int j=0;j < scales.length; j++){
					  scales[j] = 1.0;
					  if (gain_equalize){
						  scales[j] *=  avr_G/avr_g;
					  }
					  if (colors_equalize){
						  switch (j) {
						  case 1: // G1
						  case 2: // G2
//							  scales[j] *=  avr_G/avr_pix[srcChannel][j] / (avr_G/avr_g);
							  scales[j] *=  avr_g/avr_pix[srcChannel][j];
							  break;
						  default: // R, B
							  scales[j] *=  avr_RGGB[j]/avr_pix[srcChannel][j] / (avr_G/avr_g);
						  }
//						  scales[j] *=  avr_RGGB[j]/avr_pix[srcChannel][j] / (avr_G/avr_g);
					  }
				  }
				  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
				  int width =  imp_srcs[srcChannel].getWidth();
				  int height = imp_srcs[srcChannel].getHeight();
				  System.out.println("Channel "+srcChannel+ ":  scales[] = "+scales[0]+", "+scales[1]+", "+scales[2]+", "+scales[3]);

				  for (int y = 0; y < height-1; y+=2){
					  for (int x = 0; x < width-1; x+=2){
						  pixels[y*width+x        ] *= scales[1];
						  pixels[y*width+x+width+1] *= scales[2];
						  pixels[y*width+x      +1] *= scales[0];
						  pixels[y*width+x+width  ] *= scales[3];
					  }
				  }
			  }
		  }
	  }



	  public ImagePlus [] processCLTQuadCorr(
			  ImagePlus [] imp_quad, // should have properties "name"(base for saving results), "channel","path"
			  boolean [][] saturation_imp, // (near) saturated pixels or null
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
//			  int              convolveFFTSize, // 128 - fft size, kernel size should be size/2
			  double []	       scaleExposures, // probably not needed here
			  final boolean    apply_corr, // calculate and apply additional fine geometry correction
			  final boolean    infinity_corr, // calculate and apply geometry correction at infinity
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel){
		  final boolean      batch_mode = clt_parameters.batch_run; //disable any debug images
		  boolean advanced=  this.correctionsParameters.zcorrect || this.correctionsParameters.equirectangular;
		  boolean toRGB=     advanced? true: this.correctionsParameters.toRGB;
		  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?

		  // may use this.StartTime to report intermediate steps execution times
		  String name=this.correctionsParameters.getModelName((String) imp_quad[0].getProperty("name"));
		  //		int channel= Integer.parseInt((String) imp_src.getProperty("channel"));
		  String path= (String) imp_quad[0].getProperty("path");

		  ImagePlus [] results = new ImagePlus[imp_quad.length];
		  for (int i = 0; i < results.length; i++) {
			  results[i] = imp_quad[i];
			  results[i].setTitle(results[i].getTitle()+"RAW");
		  }
		  if (debugLevel>1) System.out.println("processing: "+path);
		  double [][][] double_stacks = new double [imp_quad.length][][];
		  for (int i = 0; i < double_stacks.length; i++){
			  double_stacks[i] = eyesisCorrections.bayerToDoubleStack(
					  imp_quad[i], // source Bayer image, linearized, 32-bit (float))
					  null); // no margins, no oversample
		  }

		  ImageDtt image_dtt = new ImageDtt();
		  for (int i = 0; i < double_stacks.length; i++){
			  for (int j =0 ; j < double_stacks[i][0].length; j++){
				  double_stacks[i][2][j]*=0.5; // Scale green 0.5 to compensate more pixels than R,B
			  }
		  }

		  setTiles (imp_quad[0], // set global tp.tilesX, tp.tilesY
				  clt_parameters,
				  threadsMax);


		  // temporary setting up tile task file (one integer per tile, bitmask
		  // for testing defined for a window, later the tiles to process will be calculated based on previous passes results

		  int [][]    tile_op = tp.setSameTileOp(clt_parameters,  clt_parameters.tile_task_op, debugLevel);
		  double [][] disparity_array = tp.setSameDisparity(clt_parameters.disparity); // 0.0); // [tp.tilesY][tp.tilesX] - individual per-tile expected disparity

		  //TODO: Add array of default disparity - use for combining images in force disparity mode (no correlation), when disparity is predicted from other tiles

		  double [][][][]     clt_corr_combo =   null;
		  double [][][][][]   clt_corr_partial = null; // [tp.tilesY][tp.tilesX][pair][color][(2*transform_size-1)*(2*transform_size-1)]
		  double [][]         clt_mismatch =     null; // [3*4][tp.tilesY * tp.tilesX] // transpose unapplied
		  double [][][][]     texture_tiles =    null; // [tp.tilesY][tp.tilesX]["RGBA".length()][]; // tiles will be 16x16, 2 visualization mode full 16 or overlapped
		  // undecided, so 2 modes of combining alpha - same as rgb, or use center tile only
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();

		  if (clt_parameters.correlate){
			  //			  clt_corr_combo =    new double [2][tp.tilesY][tp.tilesX][];
			  clt_corr_combo =    new double [ImageDtt.TCORR_TITLES.length][tilesY][tilesX][];
			  texture_tiles =     new double [tilesY][tilesX][][]; // ["RGBA".length()][];
			  for (int i = 0; i < tilesY; i++){
				  for (int j = 0; j < tilesX; j++){
					  for (int k = 0; k<clt_corr_combo.length; k++){
						  clt_corr_combo[k][i][j] = null;
					  }
					  //					  clt_corr_combo[1][i][j] = null;
					  texture_tiles[i][j] = null;
				  }
			  }
			  if (!infinity_corr && clt_parameters.corr_keep){
				  clt_corr_partial = new double [tilesY][tilesX][][][];
				  for (int i = 0; i < tilesY; i++){
					  for (int j = 0; j < tilesX; j++){
						  clt_corr_partial[i][j] = null;
					  }
				  }
			  }
			  if (clt_parameters.corr_mismatch || apply_corr || infinity_corr){ // added infinity_corr
				  clt_mismatch = new double [12][];
			  }
		  }
		  double [][] disparity_map = new double [ImageDtt.DISPARITY_TITLES.length][]; //[0] -residual disparity, [1] - orthogonal (just for debugging) last 4 - max pixel differences

//		  double min_corr_selected = clt_parameters.corr_normalize? clt_parameters.min_corr_normalized: clt_parameters.min_corr;
		  double min_corr_selected = clt_parameters.min_corr;
		  double [][] shiftXY = new double [4][2];
		  if (!clt_parameters.fine_corr_ignore) {
			  double [][] shiftXY0 = {
					  {clt_parameters.fine_corr_x_0,clt_parameters.fine_corr_y_0},
					  {clt_parameters.fine_corr_x_1,clt_parameters.fine_corr_y_1},
					  {clt_parameters.fine_corr_x_2,clt_parameters.fine_corr_y_2},
					  {clt_parameters.fine_corr_x_3,clt_parameters.fine_corr_y_3}};
			  shiftXY = shiftXY0;
		  }
//		  final double disparity_corr = (clt_parameters.z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/clt_parameters.z_correction);
		  double z_correction =  clt_parameters.z_correction;
		  if (clt_parameters.z_corr_map.containsKey(name)){
			  z_correction +=clt_parameters.z_corr_map.get(name);
		  }
		  final double disparity_corr = (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);
		  double [][][][][][] clt_data = image_dtt.clt_aberrations_quad_corr(
				  clt_parameters.img_dtt,       // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
				  1,                            // final int  macro_scale, // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
				  tile_op,                      // per-tile operation bit codes
				  disparity_array,              // final double            disparity,
				  double_stacks,                // final double [][][]      imade_data, // first index - number of image in a quad
				  saturation_imp,               // boolean [][] saturation_imp, // (near) saturated pixels or null
				  // correlation results - final and partial
				  clt_corr_combo,               // [tp.tilesY][tp.tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				  clt_corr_partial,             // [tp.tilesY][tp.tilesX][pair][color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				  clt_mismatch,                 // [12][tp.tilesY * tp.tilesX] // transpose unapplied. null - do not calculate
				  disparity_map,                // [2][tp.tilesY * tp.tilesX]
				  texture_tiles,                // [tp.tilesY][tp.tilesX]["RGBA".length()][];
				  imp_quad[0].getWidth(),       // final int width,
				  clt_parameters.fat_zero,      // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
				  clt_parameters.corr_sym,
				  clt_parameters.corr_offset,
				  clt_parameters.corr_red,
				  clt_parameters.corr_blue,
				  clt_parameters.corr_sigma,
//				  clt_parameters.corr_mask,
				  clt_parameters.corr_normalize, // normalize correlation results by rms
				  min_corr_selected, // 0.0001; // minimal correlation value to consider valid
				  clt_parameters.max_corr_sigma,// 1.5;  // weights of points around global max to find fractional
				  clt_parameters.max_corr_radius,
//				  clt_parameters.enhortho_width,  // 2;    // reduce weight of center correlation pixels from center (0 - none, 1 - center, 2 +/-1 from center)
//				  clt_parameters.enhortho_scale,  // 0.2;  // multiply center correlation pixels (inside enhortho_width)
				  clt_parameters.max_corr_double, // Double pass when masking center of mass to reduce preference for integer values
				  clt_parameters.corr_mode,     // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
				  clt_parameters.min_shot,       // 10.0;  // Do not adjust for shot noise if lower than
				  clt_parameters.scale_shot,     // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
				  clt_parameters.diff_sigma,     // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_threshold, // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_gauss,     // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
				  clt_parameters.min_agree,      // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				  clt_parameters.dust_remove,    // Do not reduce average weight when only one image differes much from the average
				  clt_parameters.keep_weights,   // Add port weights to RGBA stack (debug feature)
				  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
				  null,                          // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
				  clt_kernels,                   // final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				  clt_parameters.kernel_step,
				  clt_parameters.transform_size,
				  clt_parameters.clt_window,
				  shiftXY, //
				  disparity_corr, // final double              disparity_corr, // disparity at infinity

				  (clt_parameters.fcorr_ignore? null: this.fine_corr),
				  clt_parameters.corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85

				  clt_parameters.shift_x,       // final int               shiftX, // shift image horizontally (positive - right) - just for testing
				  clt_parameters.shift_y,       // final int               shiftY, // shift image vertically (positive - down)
				  clt_parameters.tileX, // -1234, // clt_parameters.tileX,         // final int               debug_tileX,
				  clt_parameters.tileY,         // final int               debug_tileY, -1234 will cause port coordinates debug images
				  (clt_parameters.dbg_mode & 64) != 0, // no fract shift
				  (clt_parameters.dbg_mode & 128) != 0, // no convolve
				  //				  (clt_parameters.dbg_mode & 256) != 0, // transpose convolve
				  threadsMax,
				  debugLevel);
		  if (debugLevel > -1){
			  System.out.println("clt_data.length="+clt_data.length+" clt_data[0].length="+clt_data[0].length
					  +" clt_data[0][0].length="+clt_data[0][0].length+" clt_data[0][0][0].length="+
					  clt_data[0][0][0].length);
		  }
//		  +" clt_data[0][0][0][0].length="+clt_data[0][0][0][0].length+
//				  " clt_data[0][0][0][0][0].length="+clt_data[0][0][0][0][0].length);
		  // visualize texture tiles as RGBA slices
		  double [][] texture_nonoverlap = null;
		  double [][] texture_overlap = null;
		  String [] rgba_titles = {"red","blue","green","alpha"};
		  String [] rgba_weights_titles = {"red","blue","green","alpha","port0","port1","port2","port3","r-rms","b-rms","g-rms","w-rms"};
		  if (texture_tiles != null){
			  if (clt_parameters.show_nonoverlap){
				  texture_nonoverlap = image_dtt.combineRGBATiles(
						  texture_tiles,                 // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}
						  clt_parameters.transform_size,
						  false,                         // when false - output each tile as 16x16, true - overlap to make 8x8
						  clt_parameters.sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
						  threadsMax,                    // maximal number of threads to launch
						  debugLevel);
				  sdfa_instance.showArrays(
						  texture_nonoverlap,
						  tilesX * (2 * clt_parameters.transform_size),
						  tilesY * (2 * clt_parameters.transform_size),
						  true,
						  name + "-TXTNOL-D"+clt_parameters.disparity,
						  (clt_parameters.keep_weights?rgba_weights_titles:rgba_titles));

			  }
			  if (!infinity_corr && (clt_parameters.show_overlap || clt_parameters.show_rgba_color)){
				  int alpha_index = 3;
				  texture_overlap = image_dtt.combineRGBATiles(
						  texture_tiles,                 // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}
						  clt_parameters.transform_size,
						  true,                         // when false - output each tile as 16x16, true - overlap to make 8x8
						  clt_parameters.sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
						  threadsMax,                    // maximal number of threads to launch
						  debugLevel);
				  if (clt_parameters.alpha1 > 0){ // negative or 0 - keep alpha as it was
					  double scale = (clt_parameters.alpha1 > clt_parameters.alpha0) ? (1.0/(clt_parameters.alpha1 - clt_parameters.alpha0)) : 0.0;
					  for (int i = 0; i < texture_overlap[alpha_index].length; i++){
						  double d = texture_overlap[alpha_index][i];
						  if      (d >=clt_parameters.alpha1) d = 1.0;
						  else if (d <=clt_parameters.alpha0) d = 0.0;
						  else d = scale * (d- clt_parameters.alpha0);
						  texture_overlap[alpha_index][i] = d;
					  }
				  }

				  if (!batch_mode && clt_parameters.show_overlap) {
					  sdfa_instance.showArrays(
							  texture_overlap,
							  tilesX * clt_parameters.transform_size,
							  tilesY * clt_parameters.transform_size,
							  true,
							  name + "-TXTOL-D"+clt_parameters.disparity,
							  (clt_parameters.keep_weights?rgba_weights_titles:rgba_titles));
				  }
				  if (!batch_mode && clt_parameters.show_rgba_color) {
					  // for now - use just RGB. Later add option for RGBA
					  double [][] texture_rgb = {texture_overlap[0],texture_overlap[1],texture_overlap[2]};
					  double [][] texture_rgba = {texture_overlap[0],texture_overlap[1],texture_overlap[2],texture_overlap[3]};
//					  ImagePlus img_texture =
					  linearStackToColor(
							  clt_parameters,
							  colorProcParameters,
							  rgbParameters,
							  name+"-texture", // String name,
							  "-D"+clt_parameters.disparity, //String suffix, // such as disparity=...
							  toRGB,
							  !this.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
							  true, // boolean saveShowIntermediate, // save/show if set globally
							  true, // boolean saveShowFinal,        // save/show result (color image?)
							  ((clt_parameters.alpha1 > 0)? texture_rgba: texture_rgb),
							  tilesX *  clt_parameters.transform_size,
							  tilesY *  clt_parameters.transform_size,
							  1.0,         // double scaleExposure, // is it needed?
							  debugLevel );
				  }
			  }
		  }
		  // visualize correlation results
		  if (clt_corr_combo!=null){
			  if (disparity_map != null){
				  if (!batch_mode && clt_parameters.show_map &&  (debugLevel > -1)){
					  sdfa_instance.showArrays(
							  disparity_map,
							  tilesX,
							  tilesY,
							  true,
							  name+"-DISP_MAP-D"+clt_parameters.disparity,
							  ImageDtt.DISPARITY_TITLES);
				  }
			  }

			  if (infinity_corr && (disparity_map != null)){
				  System.out.println("=== applying geometry correction coefficients to correct disparity at infinity ===");
				  System.out.println("=== Set inf_repeat =0 to disable automatic correction and just generate a data image to correct in offline mode ===");
				  double [] mismatch_strength = disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX].clone();
				  int [] strength_indices = {2,5,8,11};
				  for (int nt = 0; nt < mismatch_strength.length; nt++) {
					  if (Double.isNaN(mismatch_strength[nt])) {
						  mismatch_strength[nt] = 0.0;
					  } else {
						  if (mismatch_strength[nt] > 0) {
							  for (int i = 0; i < strength_indices.length; i++) {
								  if (!(mismatch_strength[nt] < clt_mismatch[strength_indices[i]][nt])) {
									  mismatch_strength[nt] = clt_mismatch[strength_indices[i]][nt];
								  }
								  //								  if (nt == 37005) {
								  //									  System.out.println("breakpoint in processCLTQuadCorr()");
								  //								  }
							  }
						  }
					  }
				  }
				  double [][] inf_ds = { // when using with disparity, programmed disparity should be 0
						  disparity_map[ImageDtt.DISPARITY_INDEX_CM],
						  mismatch_strength}; // disparity_map[ ImageDtt.DISPARITY_STRENGTH_INDEX]};
				  String [] titles = {"disp_cm", "strength"};
				  if (!batch_mode && (clt_mismatch != null)){
					  double [][] inf_ds1 = {
							  disparity_map[ImageDtt.DISPARITY_INDEX_CM],
							  mismatch_strength, //disparity_map[ ImageDtt.DISPARITY_STRENGTH_INDEX],
							  clt_mismatch[0].clone(), // dx0
							  clt_mismatch[1].clone(), // dy0 +
							  clt_mismatch[3].clone(), // dx1
							  clt_mismatch[4].clone(), // dy1 +
							  clt_mismatch[6].clone(), // dx2 +
							  clt_mismatch[7].clone(), // dy2
							  clt_mismatch[9].clone(), // dx3 +
							  clt_mismatch[10].clone()};// dy3
					  // restore full d{xy}[i] with subtracted disparity - debugging (otherwise clone() above is not neded)
						// Add disparity to dx0, dx1, dy2, dy3 pairs
						if (clt_parameters.inf_restore_disp) {
							if (debugLevel > -2) {
								System.out.println("---- Adding disparity to  d{xy}[i] ---");
							}
							for (int nTile = 0; nTile < inf_ds1[0].length; nTile++) if (inf_ds1[1][nTile] > 0){ // strength
								for (int i = 0; i < AlignmentCorrection.INDICES_10_DISP.length; i++) {
									inf_ds1[AlignmentCorrection.INDICES_10_DISP[i]][nTile] += inf_ds1[AlignmentCorrection.INDEX_10_DISPARITY][nTile];
								}
							}
						} else {
							if (debugLevel > -2) {
								System.out.println("---- d{xy}[i] have disparity canceled, xy_mismatch will only reflect residualvalues---");
							}
						}

					  String [] titles1 = {"disp_cm", "strength", "dx0", "dy0", "dx1", "dy1", "dx2", "dy2", "dx3", "dy3"};
					  inf_ds = inf_ds1;
					  titles = titles1;
				  }
				  if (clt_parameters.inf_repeat < 1) {
					  System.out.println("=== Generating image to be saved and then used for correction with 'CLT ext infinity corr'===");
					  // This image can be saved and re-read with "CLT ext infinity corr" command
					  if (!batch_mode && (sdfa_instance != null)){
						  sdfa_instance.showArrays(
								  inf_ds,
								  tilesX,
								  tilesY,
								  true,
								  name + "-inf_corr",
								  titles );
					  }
				  } else { // calculate/apply coefficients
					  if (debugLevel + (clt_parameters.fine_dbg ? 1:0) > 0){
						  // still show image, even as it is not needed
						  if (!batch_mode && (sdfa_instance != null)) {
							  sdfa_instance.showArrays(
									  inf_ds,
									  tilesX,
									  tilesY,
									  true,
									  name + "-inf_corr",
									  titles );
						  }
					  }

					  AlignmentCorrection ac = new AlignmentCorrection(this);
					  // includes both infinity correction and mismatch correction for the same infinity tiles
					  double [][][] new_corr = ac.infinityCorrection(
							  clt_parameters.ly_poly,        // final boolean use_poly,
							  clt_parameters.fcorr_inf_strength, //  final double min_strenth,
							  clt_parameters.fcorr_inf_diff,     // final double max_diff,
							  clt_parameters.inf_iters,          // 20, // 0, // final int max_iterations,
							  clt_parameters.inf_final_diff,     // 0.0001, // final double max_coeff_diff,
							  clt_parameters.inf_far_pull,       // 0.0, // 0.25, //   final double far_pull, //  = 0.2; // 1; //  0.5;
							  clt_parameters.inf_str_pow,        // 1.0, //   final double     strength_pow,
							  clt_parameters.inf_smpl_side,      // 3, //   final int        smplSide, //        = 2;      // Sample size (side of a square)
							  clt_parameters.inf_smpl_num,       // 5, //   final int        smplNum,  //         = 3;      // Number after removing worst (should be >1)
							  clt_parameters.inf_smpl_rms,       // 0.1, // 0.05, //  final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
							  // histogram parameters
							  clt_parameters.ih_smpl_step,       // 8,    // final int        hist_smpl_side, // 8 x8 masked, 16x16 sampled
							  clt_parameters.ih_disp_min,        // -1.0, // final double     hist_disp_min,
							  clt_parameters.ih_disp_step,       // 0.05, // final double     hist_disp_step,
							  clt_parameters.ih_num_bins,        // 40,   // final int        hist_num_bins,
							  clt_parameters.ih_sigma,           // 0.1,  // final double     hist_sigma,
							  clt_parameters.ih_max_diff,        // 0.1,  // final double     hist_max_diff,
							  clt_parameters.ih_min_samples,     // 10,   // final int        hist_min_samples,
							  clt_parameters.ih_norm_center,     // true, // final boolean    hist_norm_center, // if there are more tiles that fit than min_samples, replace with
							  clt_parameters,                    // EyesisCorrectionParameters.CLTParameters           clt_parameters,
							  inf_ds,                            // double [][] disp_strength,
							  tilesX,                            // int         tilesX,
							  clt_parameters.corr_magic_scale,   // double      magic_coeff, // still not understood coefficent that reduces reported disparity value.  Seems to be around 8.5
							  debugLevel + (clt_parameters.fine_dbg ? 1:0));                   // int debugLevel)

					  if ((new_corr != null) && (debugLevel > -1)){
						  System.out.println("process_infinity_corr(): ready to apply infinity correction");
						  show_fine_corr(
								  new_corr, // double [][][] corr,
								  "");// String prefix)

					  }

					  if (clt_parameters.inf_disp_apply){
						  apply_fine_corr(
								  new_corr,
								  debugLevel + 2);
					  }
				  }
			  }


			  if (!batch_mode && !infinity_corr && clt_parameters.corr_show && (debugLevel > -1)){
				  double [][] corr_rslt = new double [clt_corr_combo.length][];
				  String [] titles = new String[clt_corr_combo.length]; // {"combo","sum"};
				  for (int i = 0; i< titles.length; i++) titles[i] = ImageDtt.TCORR_TITLES[i];
				  for (int i = 0; i<corr_rslt.length; i++) {
					  corr_rslt[i] = image_dtt.corr_dbg(
							  clt_corr_combo[i],
							  2*clt_parameters.transform_size - 1,
							  clt_parameters.corr_border_contrast,
							  threadsMax,
							  debugLevel);
				  }

				  sdfa_instance.showArrays(
						  corr_rslt,
						  tilesX*(2*clt_parameters.transform_size),
						  tilesY*(2*clt_parameters.transform_size),
						  true,
						  name + "-CORR-D"+clt_parameters.disparity,
						  titles );
			  }

			  if (!batch_mode && !infinity_corr && (clt_corr_partial!=null)){
				  if (debugLevel > -1){ // -1
					  String [] allColorNames = {"red","blue","green","combo"};
					  String [] titles = new String[clt_corr_partial.length];
					  for (int i = 0; i < titles.length; i++){
						  titles[i]=allColorNames[i % allColorNames.length]+"_"+(i / allColorNames.length);
					  }
					  double [][] corr_rslt_partial = image_dtt.corr_partial_dbg(
							  clt_corr_partial,
							  2*clt_parameters.transform_size - 1,	//final int corr_size,
							  4,	// final int pairs,
							  4,    // final int colors,
							  clt_parameters.corr_border_contrast,
							  threadsMax,
							  debugLevel);
					  sdfa_instance.showArrays(
							  corr_rslt_partial,
							  tilesX*(2*clt_parameters.transform_size),
							  tilesY*(2*clt_parameters.transform_size),
							  true,
							  name+"-PART_CORR-D"+clt_parameters.disparity,
							  titles);
				  }
			  }
		  }

		  if (!infinity_corr && (clt_parameters.gen_chn_img || clt_parameters.gen_4_img || clt_parameters.gen_chn_stacks)) {
			  ImagePlus [] imps_RGB = new ImagePlus[clt_data.length];
			  for (int iQuad = 0; iQuad < clt_data.length; iQuad++){

				  String title=name+"-"+String.format("%02d", iQuad);
				  //				  String titleFull=title+"-SPLIT-D"+clt_parameters.disparity;

				  if (clt_parameters.corr_sigma > 0){ // no filter at all
					  for (int chn = 0; chn < clt_data[iQuad].length; chn++) {
						  image_dtt.clt_lpf(
								  clt_parameters.corr_sigma,
								  clt_data[iQuad][chn],
								  clt_parameters.transform_size,
								  threadsMax,
								  debugLevel);
					  }
				  }

				  if (debugLevel > 0){
					  System.out.println("--tilesX="+tilesX);
					  System.out.println("--tilesY="+tilesY);
				  }
				  if (!batch_mode && (debugLevel > 0)){
					  double [][] clt = new double [clt_data[iQuad].length*4][];
					  for (int chn = 0; chn < clt_data[iQuad].length; chn++) {
						  double [][] clt_set = image_dtt.clt_dbg(
								  clt_data [iQuad][chn],
								  threadsMax,
								  debugLevel);
						  for (int ii = 0; ii < clt_set.length; ii++) clt[chn*4+ii] = clt_set[ii];
					  }

					  if (debugLevel > 0){
						  sdfa_instance.showArrays(clt,
								  tilesX*clt_parameters.transform_size,
								  tilesY*clt_parameters.transform_size,
								  true,
								  results[iQuad].getTitle()+"-CLT-D"+clt_parameters.disparity);
					  }
				  }
				  double [][] iclt_data = new double [clt_data[iQuad].length][];
				  for (int chn=0; chn<iclt_data.length;chn++){
					  iclt_data[chn] = image_dtt.iclt_2d(
							  clt_data[iQuad][chn],           // scanline representation of dcd data, organized as dct_size x dct_size tiles
							  clt_parameters.transform_size,  // final int
							  clt_parameters.clt_window,      // window_type
							  15,                             // clt_parameters.iclt_mask,       //which of 4 to transform back
							  0,                              // clt_parameters.dbg_mode,        //which of 4 to transform back
							  threadsMax,
							  debugLevel);

				  }

				  if (clt_parameters.gen_chn_stacks) sdfa_instance.showArrays(iclt_data,
						  (tilesX + 0) * clt_parameters.transform_size,
						  (tilesY + 0) * clt_parameters.transform_size,
						  true,
						  results[iQuad].getTitle()+"-ICLT-RGB-D"+clt_parameters.disparity);
				  if (!clt_parameters.gen_chn_img) continue;

				  imps_RGB[iQuad] = linearStackToColor(
						  clt_parameters,
						  colorProcParameters,
						  rgbParameters,
						  title, // String name,
						  "-D"+clt_parameters.disparity, //String suffix, // such as disparity=...
						  toRGB,
						  !this.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
						  !batch_mode, // true, // boolean saveShowIntermediate, // save/show if set globally
						  false, // boolean saveShowFinal,        // save/show result (color image?)
						  iclt_data,
						  tilesX *  clt_parameters.transform_size,
						  tilesY *  clt_parameters.transform_size,
						  scaleExposures[iQuad], // double scaleExposure, // is it needed?
						  debugLevel );

			  } // end of generating shifted channel images
			  if (clt_parameters.gen_chn_img) {
				  // combine to a sliced color image
				  int [] slice_seq = {0,1,3,2}; //clockwise
				  int width = imps_RGB[0].getWidth();
				  int height = imps_RGB[0].getHeight();
				  ImageStack array_stack=new ImageStack(width,height);
				  for (int i = 0; i<slice_seq.length; i++){
					  if (imps_RGB[slice_seq[i]] != null) {
						  array_stack.addSlice("port_"+slice_seq[i], imps_RGB[slice_seq[i]].getProcessor().getPixels());
					  } else {
						  array_stack.addSlice("port_"+slice_seq[i], results[slice_seq[i]].getProcessor().getPixels());
					  }
				  }
				  ImagePlus imp_stack = new ImagePlus(name+"-SHIFTED-D"+clt_parameters.disparity, array_stack);
				  imp_stack.getProcessor().resetMinAndMax();
				  if (!batch_mode) {
					  imp_stack.updateAndDraw();
				  }
				  //imp_stack.getProcessor().resetMinAndMax();
				  //imp_stack.show();
				  //				  eyesisCorrections.saveAndShow(imp_stack, this.correctionsParameters);
				  eyesisCorrections.saveAndShowEnable(
						  imp_stack,  // ImagePlus             imp,
						  this.correctionsParameters, // EyesisCorrectionParameters.CorrectionParameters  correctionsParameters,
						  true, // boolean               enableSave,
						  !batch_mode) ;// boolean               enableShow);
			  }
			  if (clt_parameters.gen_4_img) {
				  // Save as individual JPEG images in the model directory
				  String x3d_path= correctionsParameters.selectX3dDirectory(
						  name, // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
						  correctionsParameters.x3dModelVersion,
						  true,  // smart,
						  true);  //newAllowed, // save
				  for (int sub_img = 0; sub_img < 4; sub_img++){
					  eyesisCorrections.saveAndShow(
							  imps_RGB[sub_img],
							  x3d_path,
							  correctionsParameters.png && !clt_parameters.black_back,
							  !batch_mode && clt_parameters.show_textures,
							  correctionsParameters.JPEG_quality, // jpegQuality); // jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
							  (debugLevel > 0) ? debugLevel : 1); // int debugLevel (print what it saves)
				  }
			  }
		  }
		  return results;
	  }

	  double [][] resizeGridTexture(
			  double [][] imgData,
			  int tileSize,
			  int tilesX,
			  int tilesY,
			  Rectangle bounds){
		  int width =   tileSize*bounds.width;
		  int height =  tileSize*bounds.height;
		  // Adding row/column of all 0, assuming java zeroes arrays
		  int numSlices = imgData.length;
		  double [][] rslt = new double [numSlices][width*height];
		  int offset = (tileSize*bounds.y)*tileSize*tilesX + (tileSize*bounds.x);
		  for (int y = 0; y < height; y ++){
			  for (int x = 0; x < width; x ++){
				  int indx = width * y + x;
				  int indx_in = indx + offset;
				  for (int i = 0; i < numSlices; i++) {
					  rslt[i][indx] = Double.isNaN(imgData[i][indx_in])? 0.0:imgData[i][indx_in];
				  }
			  }
			  offset += tilesX * tileSize - width;

		  }
		  return rslt;
	  }

	  public ImagePlus linearStackToColor(
			  EyesisCorrectionParameters.CLTParameters         clt_parameters,
			  EyesisCorrectionParameters.ColorProcParameters   colorProcParameters,
			  EyesisCorrectionParameters.RGBParameters         rgbParameters,
			  String name,
			  String suffix, // such as disparity=...
			  boolean toRGB,
			  boolean bpp16, // 16-bit per channel color mode for result
			  boolean saveShowIntermediate, // save/show if set globally
			  boolean saveShowFinal,        // save/show result (color image?)
			  double [][] iclt_data,
			  int width, // int tilesX,
			  int height, // int tilesY,
			  double scaleExposure,
			  int debugLevel
			  )
	  {
		  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		  // convert to ImageStack of 3 slices
		  String [] sliceNames = {"red", "blue", "green"};
		  double []   alpha = null; // (0..1.0)
		  double [][] rgb_in = {iclt_data[0],iclt_data[1],iclt_data[2]};
		  if (iclt_data.length > 3) alpha = iclt_data[3];
		  ImageStack stack = sdfa_instance.makeStack(
				  rgb_in, // iclt_data,
				  width,  // (tilesX + 0) * clt_parameters.transform_size,
				  height, // (tilesY + 0) * clt_parameters.transform_size,
				  sliceNames,  // or use null to get chn-nn slice names
				  true); // replace NaN with 0.0
		  if (debugLevel > -1) { // 0){
			  double [] chn_avg = {0.0,0.0,0.0};
			  float [] pixels;
//			  int width =  stack.getWidth();
//			  int height = stack.getHeight();

			  for (int c = 0; c <3; c++){
				  pixels = (float[]) stack.getPixels(c+1);
				  for (int i = 0; i<pixels.length; i++){
					  chn_avg[c] += pixels[i];
				  }
			  }
			  chn_avg[0] /= width*height;
			  chn_avg[1] /= width*height;
			  chn_avg[2] /= width*height;
			  System.out.println("Processed channels averages: R="+chn_avg[0]+", G="+chn_avg[2]+", B="+chn_avg[1]);
		  }

		  if (debugLevel > 1) System.out.println("before colors.1");
//		  if (debugLevel > -1) System.out.println("before colors.1");
		  //Processing colors - changing stack sequence to r-g-b (was r-b-g)
		  if (!eyesisCorrections.fixSliceSequence(
				  stack,
				  debugLevel)){
			  if (debugLevel > -1) System.out.println("fixSliceSequence() returned false");
			  return null;
		  }
		  if (debugLevel > 1) System.out.println("before colors.2");
//		  if (debugLevel > -1) System.out.println("before colors.2");
		  if (saveShowIntermediate && (debugLevel > 1)){
//		  if (saveShowIntermediate && (debugLevel > -1)){
			  ImagePlus imp_dbg=new ImagePlus(name+"-preColors",stack);
			  eyesisCorrections.saveAndShow(
					  imp_dbg,
					  this.correctionsParameters);
		  }
//		  if (debugLevel > 1) System.out.println("before colors.3, scaleExposure="+scaleExposure+" scale = "+(255.0/eyesisCorrections.psfSubpixelShouldBe4/eyesisCorrections.psfSubpixelShouldBe4/scaleExposure));
		  if (debugLevel > -1) System.out.println("before colors.3, scaleExposure="+scaleExposure+" scale = "+(255.0/eyesisCorrections.psfSubpixelShouldBe4/eyesisCorrections.psfSubpixelShouldBe4/scaleExposure));
		  CorrectionColorProc correctionColorProc=new CorrectionColorProc(eyesisCorrections.stopRequested);
		  double [][] yPrPb=new double [3][];
		  //			if (dct_parameters.color_DCT){
		  // need to get YPbPr - not RGB here
		  //			} else {
		  correctionColorProc.processColorsWeights(stack, // just gamma convert? TODO: Cleanup? Convert directly form the linear YPrPb
				  //					  255.0/this.psfSubpixelShouldBe4/this.psfSubpixelShouldBe4, //  double scale,     // initial maximal pixel value (16))
				  //					  255.0/eyesisCorrections.psfSubpixelShouldBe4/eyesisCorrections.psfSubpixelShouldBe4/scaleExposure, //  double scale,     // initial maximal pixel value (16))
				  //					  255.0/2/2/scaleExposure, //  double scale,     // initial maximal pixel value (16))
				  255.0/scaleExposure, //  double scale,     // initial maximal pixel value (16))
				  colorProcParameters,
				  null, // channelGainParameters,
				  -1, // channel,
				  null, //correctionDenoise.getDenoiseMask(),
				  this.correctionsParameters.blueProc,
				  debugLevel);
		  if (debugLevel > 1) System.out.println("Processed colors to YPbPr, total number of slices="+stack.getSize());
		  if (saveShowIntermediate && (debugLevel > 1)) {
			  ImagePlus imp_dbg=new ImagePlus("procColors",stack);
			  eyesisCorrections.saveAndShow(
					  imp_dbg,
					  this.correctionsParameters);
		  }
		  float [] fpixels;
		  int [] slices_YPrPb = {8,6,7};
		  yPrPb=new double [3][];
		  for (int n = 0; n < slices_YPrPb.length; n++){
			  fpixels = (float[]) stack.getPixels(slices_YPrPb[n]);
			  yPrPb[n] = new double [fpixels.length];
			  for (int i = 0; i < fpixels.length; i++) yPrPb[n][i] = fpixels[i];
		  }
		  String titleFull = "";
		  if (toRGB) {
			  if (debugLevel > 0){
				  System.out.println("correctionColorProc.YPrPbToRGB");
			  }
			  stack =  YPrPbToRGB(yPrPb,
					  colorProcParameters.kr,        // 0.299;
					  colorProcParameters.kb,        // 0.114;
					  stack.getWidth());
			  titleFull=name+"-RGB-float"+suffix;
			  //Trim stack to just first 3 slices
			  if (saveShowIntermediate && (debugLevel > 1)){ // 2){
				  ImagePlus imp_dbg=new ImagePlus("YPrPbToRGB",stack);
				  eyesisCorrections.saveAndShow(
						  imp_dbg,
						  this.correctionsParameters);
			  }
			  while (stack.getSize() > 3) stack.deleteLastSlice();
			  if (debugLevel > 1) System.out.println("Trimming color stack");
		  } else {
			  titleFull=name+"-YPrPb"+suffix;
			  if (debugLevel > 1) System.out.println("Using full stack, including YPbPr");
		  }
		  if (alpha != null){
			  float [] alpha_pixels = new float [alpha.length];
			  for (int i = 0; i <alpha.length; i++){
				  alpha_pixels[i] = (float) alpha[i];
			  }
			  stack.addSlice("alpha",alpha_pixels);
		  }

		  ImagePlus result= new ImagePlus(titleFull, stack);
		  if (debugLevel> 1){
			  result.show();
		  }

		  if (!toRGB && !this.correctionsParameters.jpeg){ // toRGB set for equirectangular
			  if (debugLevel > 1) System.out.println("!toRGB && !this.correctionsParameters.jpeg");
			  if (saveShowIntermediate) eyesisCorrections.saveAndShow(result, this.correctionsParameters);
			  return result;
		  } else { // that's not the end result, save if required
			  if (debugLevel > 1) System.out.println("!toRGB && !this.correctionsParameters.jpeg - else");


			  if (saveShowIntermediate) eyesisCorrections.saveAndShow(result, // saves OK with alpha - 32-bit float
					  eyesisCorrections.correctionsParameters,
					  eyesisCorrections.correctionsParameters.save32,
					  false, // true, // false,
					  eyesisCorrections.correctionsParameters.JPEG_quality); // save, no show
		  }
		  stack=eyesisCorrections.convertRGB32toRGB16Stack(
				  stack,
				  rgbParameters);

		  titleFull=name+"-RGB48"+suffix;
		  result= new ImagePlus(titleFull, stack);
		  result.updateAndDraw();
		  if (debugLevel > 1) {
//		  if (debugLevel > -1) {
			  System.out.println("result.updateAndDraw(), "+titleFull+"-RGB48");
			  result.show();
		  }

		  CompositeImage compositeImage=eyesisCorrections.convertToComposite(result);

		  if (!this.correctionsParameters.jpeg && bpp16){ // RGB48 was the end result
			  if (debugLevel > 1) System.out.println("if (!this.correctionsParameters.jpeg && !advanced)");
			  if (saveShowIntermediate) eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters);
			  return compositeImage; // return result;
		  } else { // that's not the end result, save if required
			  if (debugLevel > 1) System.out.println("if (!this.correctionsParameters.jpeg && !advanced) - else");
			  if (saveShowIntermediate) eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters, this.correctionsParameters.save16, false); // save, no show
		  }
		  result = eyesisCorrections.convertRGB48toRGB24(
				  stack,
//				  name+"-RGB24"+suffix,
				  name+suffix,
				  0, 65536, // r range 0->0, 65536->256
				  0, 65536, // g range
				  0, 65536,// b range
				  0, 65536);// alpha range
		  // next will save either JPEG (if no alpha) or RGBA tiff (if alpha is present). ImageJ shows just RGB (no alpha)
		  if (saveShowFinal) eyesisCorrections.saveAndShow(result, this.correctionsParameters);
		  return result;

	  }

	  public void apply_fine_corr(
			  double [][][] corr,
			  int debugLevel)
	  {
		  if (debugLevel > 1){
			  if (debugLevel > 1){
				  show_fine_corr( this.fine_corr, "  was");
			  }
		  }
		  if (corr==null) {
			  System.out.println("New correction is null (only non-null for poly, not for infinity");
			  return;
		  } else {
			  if (debugLevel > 1){
				  show_fine_corr(corr, "added");
			  }
			  for (int n = 0; n< corr.length; n++){
				  for (int i = 0; i< corr[n].length; i++){
					  for (int j = 0; j< corr[n][i].length; j++){
						  this.fine_corr[n][i][j]+=corr[n][i][j];
					  }
				  }
			  }
			  if (debugLevel > 0){
				  show_fine_corr( this.fine_corr, "");
			  }
		  }
	  }
	  public void show_fine_corr()
	  {
		  show_fine_corr("");
	  }
	  public void show_fine_corr(String prefix)
	  {
		  show_fine_corr( this.fine_corr, prefix);
	  }

	  public void show_fine_corr(
			  double [][][] corr,
			  String prefix)
	  {
		  String sadd = (prefix.length() > 0)?(prefix+" "):"";

		  for (int n = 0; n< corr.length; n++){
			  for (int i = 0; i< corr[n].length; i++){
				  System.out.print(sadd+"port"+n+": "+fine_corr_dir_names[i]+": ");
				  for (int j = 0; j< corr[n][i].length; j++){
					  System.out.print(fine_corr_coeff_names[j]+"="+corr[n][i][j]+" ");
				  }
				  System.out.println();
			  }
		  }
		  System.out.println();
		  String name = (sadd.length() == 0)?"":("("+sadd+")");
		  showExtrinsicCorr(name);
	  }


	  public void reset_fine_corr()
	  {
		  this.fine_corr = new double [4][2][6]; // reset all coefficients to 0
	  }

	  public void showExtrinsicCorr(String name)
	  {
		  System.out.println("Extrinsic corrections "+name);
		  if (geometryCorrection == null){
			  System.out.println("are not set, will be:");
			  System.out.println(new GeometryCorrection(this.extrinsic_corr).getCorrVector().toString());
		  } else {
			  System.out.println(geometryCorrection.getCorrVector().toString());
		  }
	  }

	  public boolean editRig()
	  {
		  if (!is_aux) {
			  System.out.println("Rig offsets can only be edited for the auxiliary camera, not for the main one");
			  return false;
		  }
//		  GeometryCorrection gc = this.geometryCorrection;
		  if (this.geometryCorrection == null){
			  System.out.println("geometryCorrection is not set, creating one");
			  this.geometryCorrection = new GeometryCorrection(this.extrinsic_corr);
		  }
		  boolean edited = this.geometryCorrection.editRig();
//		  if (edited) {
//			  gc.rigOffset.setProperties(prefix,properties);
//		  }
		  return edited;
	  }
/*
		if (is_aux && (gc.rigOffset != null)) {
			gc.rigOffset.setProperties(prefix,properties);
		}

 */


	  public void resetExtrinsicCorr(
			  EyesisCorrectionParameters.CLTParameters           clt_parameters)
	  {
		  this.extrinsic_corr = new double [GeometryCorrection.CORR_NAMES.length];
		  if (geometryCorrection != null){
			  geometryCorrection.setCorrVector(null);
		  }
		  if (clt_parameters.fine_corr_apply){
			  clt_parameters.fine_corr_ignore = false;
		  }
	  }

	  public void cltDisparityScans(
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
//			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
//			  int          convolveFFTSize, // 128 - fft size, kernel size should be size/2
			  final int          threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel)
	  {
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
					  // removeUnusedSensorData should be off!?
					  channels=this.eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (eyesisCorrections.isChannelEnabled(channels[i])){
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
		  double [] referenceExposures=eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		  int [][] fileIndices=new int [numImagesToProcess][2]; // file index, channel number
		  int index=0;
		  for (int nFile=0;nFile<enabledFiles.length;nFile++){
			  if ((sourceFiles[nFile]!=null) && (sourceFiles[nFile].length()>1)) {
				  int [] channels={correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile])};
				  if (correctionsParameters.isJP4()){
					  int subCamera= channels[0]- correctionsParameters.firstSubCamera; // to match those in the sensor files
					  channels=eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (eyesisCorrections.isChannelEnabled(channels[i])){
						  fileIndices[index  ][0]=nFile;
						  fileIndices[index++][1]=channels[i];
					  }
				  }
			  }
		  }
		  ArrayList<String> setNames = new ArrayList<String>();
		  ArrayList<ArrayList<Integer>> setFiles = new ArrayList<ArrayList<Integer>>();

		  for (int iImage=0;iImage<fileIndices.length;iImage++){
			  int nFile=fileIndices[iImage][0];
			  String setName = correctionsParameters.getNameFromSourceTiff(sourceFiles[nFile]);
			  if (!setNames.contains(setName)) {
				  setNames.add(setName);
				  setFiles.add(new ArrayList<Integer>());
			  }
			  setFiles.get(setNames.indexOf(setName)).add(new Integer(nFile));
		  }
		  for (int nSet = 0; nSet < setNames.size(); nSet++){
			  int maxChn = 0;
			  for (int i = 0; i < setFiles.get(nSet).size(); i++){
				  int chn = fileIndices[setFiles.get(nSet).get(i)][1];
				  if (chn > maxChn) maxChn = chn;
			  }
			  int [] channelFiles = new int[maxChn+1];
			  for (int i =0; i < channelFiles.length; i++) channelFiles[i] = -1;
			  for (int i = 0; i < setFiles.get(nSet).size(); i++){
				  channelFiles[fileIndices[setFiles.get(nSet).get(i)][1]] = setFiles.get(nSet).get(i);
			  }

			  ImagePlus [] imp_srcs = new ImagePlus[channelFiles.length];
			  boolean [][] saturation_imp = (clt_parameters.sat_level > 0.0)? new boolean[channelFiles.length][] : null;
			  double [] scaleExposures = new double[channelFiles.length];
			  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
				  int nFile=channelFiles[srcChannel];
				  imp_srcs[srcChannel]=null;
				  if (nFile >=0){
					  if (correctionsParameters.isJP4()){
						  int subchannel=eyesisCorrections.pixelMapping.getSubChannel(srcChannel);
						  if (this.correctionsParameters.swapSubchannels01) {
							  switch (subchannel){
							  case 0: subchannel=1; break;
							  case 1: subchannel=0; break;
							  }
						  }
						  if (debugLevel>0) System.out.println("Processing set " + setNames.get(nSet)+" channel "+srcChannel+" - subchannel "+subchannel+" of "+sourceFiles[nFile]);
						  ImagePlus imp_composite=eyesisCorrections.JP4_INSTANCE.open(
								  "", // path,
								  sourceFiles[nFile],
								  "",  //arg - not used in JP46 reader
								  true, // un-apply camera color gains
								  null, // new window
								  false); // do not show
						  imp_srcs[srcChannel]=eyesisCorrections.JP4_INSTANCE.demuxImage(imp_composite, subchannel);
						  if (imp_srcs[srcChannel] == null) imp_srcs[srcChannel] = imp_composite; // not a composite image
						  // do we need to add any properties?
					  } else {
						  imp_srcs[srcChannel]=new ImagePlus(sourceFiles[nFile]);
						  //					  (new JP46_Reader_camera(false)).decodeProperiesFromInfo(imp_src); // decode existent properties from info
						  eyesisCorrections.JP4_INSTANCE.decodeProperiesFromInfo(imp_srcs[srcChannel]); // decode existent properties from info
						  if (debugLevel>0) System.out.println("Processing "+sourceFiles[nFile]);
					  }
					  imp_srcs[srcChannel] =  padBayerToFullSize(
							  imp_srcs[srcChannel], // ImagePlus imp_src,
							  eyesisCorrections.pixelMapping.sensors[srcChannel].getSensorWH(),
							  true); // boolean replicate);
					  scaleExposures[srcChannel] = 1.0;
					  if (!Double.isNaN(referenceExposures[nFile]) && (imp_srcs[srcChannel].getProperty("EXPOSURE")!=null)){
						  scaleExposures[srcChannel] = referenceExposures[nFile]/Double.parseDouble((String) imp_srcs[srcChannel].getProperty("EXPOSURE"));
						  if (debugLevel>0) System.out.println("Will scale intensity (to compensate for exposure) by  "+scaleExposures[srcChannel]);
					  }
					  imp_srcs[srcChannel].setProperty("name",    correctionsParameters.getNameFromSourceTiff(sourceFiles[nFile]));
					  imp_srcs[srcChannel].setProperty("channel", srcChannel); // it may already have channel
					  imp_srcs[srcChannel].setProperty("path",    sourceFiles[nFile]); // it may already have channel

					  if (this.correctionsParameters.pixelDefects && (eyesisCorrections.defectsXY!=null)&& (eyesisCorrections.defectsXY[srcChannel]!=null)){
						  // apply pixel correction
						  int numApplied=	eyesisCorrections.correctDefects(
								  imp_srcs[srcChannel],
								  srcChannel,
								  debugLevel);
						  if ((debugLevel>0) && (numApplied>0)) { // reduce verbosity after verified defect correction works
							  System.out.println("Corrected "+numApplied+" pixels in "+sourceFiles[nFile]);
						  }
					  }
					  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
					  int width =  imp_srcs[srcChannel].getWidth();
					  int height = imp_srcs[srcChannel].getHeight();

					  if (clt_parameters.sat_level > 0.0){
						  double [] saturations = {
								  Double.parseDouble((String) imp_srcs[srcChannel].getProperty("saturation_1")),
								  Double.parseDouble((String) imp_srcs[srcChannel].getProperty("saturation_0")),
								  Double.parseDouble((String) imp_srcs[srcChannel].getProperty("saturation_3")),
								  Double.parseDouble((String) imp_srcs[srcChannel].getProperty("saturation_2"))};
						  saturation_imp[srcChannel] = new boolean[width*height];
						  System.out.println(String.format("channel %d saturations = %6.2f %6.2f %6.2f %6.2f", srcChannel,
								  saturations[0],saturations[1],saturations[2],saturations[3]));
						  double [] scaled_saturations = new double [saturations.length];
						  for (int i = 0; i < scaled_saturations.length; i++){
							  scaled_saturations[i] = saturations[i] * clt_parameters.sat_level;
						  }
						  for (int y = 0; y < height-1; y+=2){
							  for (int x = 0; x < width-1; x+=2){
								  if (pixels[y*width+x        ] > scaled_saturations[0])  saturation_imp[srcChannel][y*width+x        ] = true;
								  if (pixels[y*width+x+      1] > scaled_saturations[1])  saturation_imp[srcChannel][y*width+x      +1] = true;
								  if (pixels[y*width+x+width  ] > scaled_saturations[2])  saturation_imp[srcChannel][y*width+x+width  ] = true;
								  if (pixels[y*width+x+width+1] > scaled_saturations[3])  saturation_imp[srcChannel][y*width+x+width+1] = true;
							  }
						  }
					  }



					  if (this.correctionsParameters.vignetting){
						  if ((eyesisCorrections.channelVignettingCorrection==null) || (srcChannel<0) || (srcChannel>=eyesisCorrections.channelVignettingCorrection.length) || (eyesisCorrections.channelVignettingCorrection[srcChannel]==null)){
							  System.out.println("No vignetting data for channel "+srcChannel);
							  return;
						  }
///						  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
						  if (pixels.length!=eyesisCorrections.channelVignettingCorrection[srcChannel].length){
							  System.out.println("Vignetting data for channel "+srcChannel+" has "+eyesisCorrections.channelVignettingCorrection[srcChannel].length+" pixels, image "+sourceFiles[nFile]+" has "+pixels.length);
							  return;
						  }
						  // TODO: Move to do it once:
						  double min_non_zero = 0.0;
						  for (int i=0;i<pixels.length;i++){
							  double d = eyesisCorrections.channelVignettingCorrection[srcChannel][i];
							  if ((d > 0.0) && ((min_non_zero == 0) || (min_non_zero > d))){
								  min_non_zero = d;
							  }
						  }
						  double max_vign_corr = clt_parameters.vignetting_range*min_non_zero;

						  System.out.println("Vignetting data: channel="+srcChannel+", min = "+min_non_zero);
						  for (int i=0;i<pixels.length;i++){
							  double d = eyesisCorrections.channelVignettingCorrection[srcChannel][i];
							  if (d > max_vign_corr) d = max_vign_corr;
							  pixels[i]*=d;
						  }
						  // Scale here, combine with vignetting later?
///						  int width =  imp_srcs[srcChannel].getWidth();
///						  int height = imp_srcs[srcChannel].getHeight();
						  for (int y = 0; y < height-1; y+=2){
							  for (int x = 0; x < width-1; x+=2){
								  pixels[y*width+x        ] *= clt_parameters.scale_g;
								  pixels[y*width+x+width+1] *= clt_parameters.scale_g;
								  pixels[y*width+x      +1] *= clt_parameters.scale_r;
								  pixels[y*width+x+width  ] *= clt_parameters.scale_b;
							  }
						  }

					  } else { // assuming GR/BG pattern
						  System.out.println("Applying fixed color gain correction parameters: Gr="+
								  clt_parameters.novignetting_r+", Gg="+clt_parameters.novignetting_g+", Gb="+clt_parameters.novignetting_b);
///						  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
///						  int width =  imp_srcs[srcChannel].getWidth();
///						  int height = imp_srcs[srcChannel].getHeight();
						  double kr = clt_parameters.scale_r/clt_parameters.novignetting_r;
						  double kg = clt_parameters.scale_g/clt_parameters.novignetting_g;
						  double kb = clt_parameters.scale_b/clt_parameters.novignetting_b;
						  for (int y = 0; y < height-1; y+=2){
							  for (int x = 0; x < width-1; x+=2){
								  pixels[y*width+x        ] *= kg;
								  pixels[y*width+x+width+1] *= kg;
								  pixels[y*width+x      +1] *= kr;
								  pixels[y*width+x+width  ] *= kb;
							  }
						  }
					  }
				  }
			  }

			  System.out.println("Temporarily applying scaleExposures[] here" );
			  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
				  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
				  for (int i = 0; i < pixels.length; i++){
					  pixels[i] *= scaleExposures[srcChannel];
				  }
				  scaleExposures[srcChannel] = 1.0;
			  }

			  // once per quad here
			  // may need to equalize gains between channels
			  if (clt_parameters.gain_equalize || clt_parameters.colors_equalize){
				  channelGainsEqualize(
						  clt_parameters.gain_equalize,
						  clt_parameters.colors_equalize,
						  clt_parameters.nosat_equalize, // boolean nosat_equalize,
						  channelFiles,
						  imp_srcs,
						  saturation_imp, // boolean[][] saturated,
						  setNames.get(nSet), // just for debug messeges == setNames.get(nSet)
						  debugLevel);
			  }
			  // once per quad here
			  cltDisparityScan( // returns ImagePlus, but it already should be saved/shown
					  imp_srcs,       // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
					  saturation_imp, // boolean [][] saturation_imp, // (near) saturated pixels or null
					  clt_parameters,
					  debayerParameters,
//					  nonlinParameters,
					  colorProcParameters,
					  channelGainParameters,
					  rgbParameters,
//					  convolveFFTSize, // 128 - fft size, kernel size should be size/2
					  scaleExposures,
					  threadsMax,  // maximal number of threads to launch
					  updateStatus,
					  debugLevel);
			  Runtime.getRuntime().gc();
			  if (debugLevel >-1) System.out.println("Processing set "+(nSet+1)+" (of "+fileIndices.length+") finished at "+
					  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
			  if (eyesisCorrections.stopRequested.get()>0) {
				  System.out.println("User requested stop");
				  return;
			  }
		  }
		  System.out.println("cltDisparityScans(): processing "+fileIndices.length+" files finished at "+
				  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
	  }


	  public ImagePlus [] cltDisparityScan(
			  ImagePlus [] imp_quad, // should have properties "name"(base for saving results), "channel","path"
			  boolean [][] saturation_imp, // (near) saturated pixels or null
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
//			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
//			  int              convolveFFTSize, // 128 - fft size, kernel size should be size/2
			  double []	       scaleExposures, // probably not needed here
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel){
		  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?

		  // may use this.StartTime to report intermediate steps execution times
		  String name=(String) imp_quad[0].getProperty("name");
		  //		int channel= Integer.parseInt((String) imp_src.getProperty("channel"));
		  String path= (String) imp_quad[0].getProperty("path");

		  ImagePlus [] results = new ImagePlus[imp_quad.length];
		  for (int i = 0; i < results.length; i++) {
			  results[i] = imp_quad[i];
			  results[i].setTitle(results[i].getTitle()+"RAW");
		  }
		  if (debugLevel>1) System.out.println("processing: "+path);
		  double [][][] double_stacks = new double [imp_quad.length][][];
		  for (int i = 0; i < double_stacks.length; i++){
			  double_stacks[i] = eyesisCorrections.bayerToDoubleStack(
					  imp_quad[i], // source Bayer image, linearized, 32-bit (float))
					  null); // no margins, no oversample
		  }
		  // =================
		  ImageDtt image_dtt = new ImageDtt();
		  for (int i = 0; i < double_stacks.length; i++){
			  for (int j =0 ; j < double_stacks[i][0].length; j++){
				  double_stacks[i][2][j]*=0.5; // Scale green 0.5 to compensate more pixels than R,B
			  }
		  }

		  setTiles (imp_quad[0], // set global tp.tilesX, tp.tilesY
				  clt_parameters,
				  threadsMax);
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();

		  // temporary setting up tile task file (one integer per tile, bitmask
		  // for testing defined for a window, later the tiles to process will be calculated based on previous passes results
		  int [][]    tile_op = tp.setSameTileOp(clt_parameters,  clt_parameters.tile_task_op, debugLevel);

		  //TODO: Add array of default disparity - use for combining images in force disparity mode (no correlation), when disparity is predicted from other tiles

		  // undecided, so 2 modes of combining alpha - same as rgb, or use center tile only
		  double [][][][]     clt_corr_combo =    new double [ImageDtt.TCORR_TITLES.length][tilesY][tilesX][]; // will only be used inside?
		  double min_corr_selected = clt_parameters.min_corr;

		  double [][][] disparity_maps = new double [clt_parameters.disp_scan_count][ImageDtt.DISPARITY_TITLES.length][]; //[0] -residual disparity, [1] - orthogonal (just for debugging)
		  double [][][] clt_mismatches = new double [clt_parameters.disp_scan_count][12][];
		  for (int scan_step = 0; scan_step < clt_parameters.disp_scan_count; scan_step++) {
			  double disparity = clt_parameters.disp_scan_start + scan_step * clt_parameters.disp_scan_step;
			  double [][] disparity_array = tp.setSameDisparity(disparity); // [tp.tilesY][tp.tilesX] - individual per-tile expected disparity

			  double [][] shiftXY = new double [4][2];
			  if (!clt_parameters.fine_corr_ignore) {
				  double [][] shiftXY0 = {
						  {clt_parameters.fine_corr_x_0,clt_parameters.fine_corr_y_0},
						  {clt_parameters.fine_corr_x_1,clt_parameters.fine_corr_y_1},
						  {clt_parameters.fine_corr_x_2,clt_parameters.fine_corr_y_2},
						  {clt_parameters.fine_corr_x_3,clt_parameters.fine_corr_y_3}};
				  shiftXY = shiftXY0;
			  }

//			  final double disparity_corr = (clt_parameters.z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/clt_parameters.z_correction);
			  double z_correction =  clt_parameters.z_correction;
			  if (clt_parameters.z_corr_map.containsKey(image_name)){
				  z_correction +=clt_parameters.z_corr_map.get(image_name);
			  }
			  final double disparity_corr = (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);

			  image_dtt.clt_aberrations_quad_corr(
					  clt_parameters.img_dtt,       // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					  1,                            // final int  macro_scale, // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
					  tile_op,                      // per-tile operation bit codes
					  disparity_array,              // clt_parameters.disparity,     // final double            disparity,
					  double_stacks,                // final double [][][]      imade_data, // first index - number of image in a quad
					  saturation_imp,               // boolean [][] saturation_imp, // (near) saturated pixels or null
					  // correlation results - final and partial
					  clt_corr_combo,               // [tp.tilesY][tp.tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					  null, // clt_corr_partial,    // [tp.tilesY][tp.tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					  clt_mismatches[scan_step], // null, [12][tilesY * tilesX] // transpose unapplied. null - do not calculate

//	Use it with disparity_maps[scan_step]?		  clt_mismatch,    // [tp.tilesY][tp.tilesX][pair]{dx,dy,weight}[(2*transform_size-1)*(2*transform_size-1)] // transpose unapplied. null - do not calculate
					  disparity_maps[scan_step],    // [2][tp.tilesY * tp.tilesX]
					  null, //texture_tiles,        // [tp.tilesY][tp.tilesX]["RGBA".length()][];
					  imp_quad[0].getWidth(),       // final int width,
					  clt_parameters.fat_zero,      // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
					  clt_parameters.corr_sym,
					  clt_parameters.corr_offset,
					  clt_parameters.corr_red,
					  clt_parameters.corr_blue,
					  clt_parameters.corr_sigma,
					  clt_parameters.corr_normalize, // normalize correlation results by rms
					  min_corr_selected, // 0.0001; // minimal correlation value to consider valid
					  clt_parameters.max_corr_sigma,// 1.5;  // weights of points around global max to find fractional
					  clt_parameters.max_corr_radius,
//					  clt_parameters.enhortho_width,  // 2;    // reduce weight of center correlation pixels from center (0 - none, 1 - center, 2 +/-1 from center)
//					  clt_parameters.enhortho_scale,  // 0.2;  // multiply center correlation pixels (inside enhortho_width)
					  clt_parameters.max_corr_double, // Double pass when masking center of mass to reduce preference for integer values
					  clt_parameters.corr_mode,     // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
					  clt_parameters.min_shot,       // 10.0;  // Do not adjust for shot noise if lower than
					  clt_parameters.scale_shot,     // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
					  clt_parameters.diff_sigma,     // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
					  clt_parameters.diff_threshold, // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
					  clt_parameters.diff_gauss,     // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
					  clt_parameters.min_agree,      // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
					  clt_parameters.dust_remove,    // Do not reduce average weight when only one image differes much from the average
					  clt_parameters.keep_weights,   // Add port weights to RGBA stack (debug feature)
					  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
					  null,                          // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
					  clt_kernels,                   // final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
					  clt_parameters.kernel_step,
					  clt_parameters.transform_size,
					  clt_parameters.clt_window,
					  shiftXY, //
					  disparity_corr, // final double              disparity_corr, // disparity at infinity
					  (clt_parameters.fcorr_ignore? null: this.fine_corr),
					  clt_parameters.corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
					  clt_parameters.shift_x,       // final int               shiftX, // shift image horizontally (positive - right) - just for testing
					  clt_parameters.shift_y,       // final int               shiftY, // shift image vertically (positive - down)
					  clt_parameters.tileX,         // final int               debug_tileX,
					  clt_parameters.tileY,         // final int               debug_tileY,
					  (clt_parameters.dbg_mode & 64) != 0, // no fract shift
					  (clt_parameters.dbg_mode & 128) != 0, // no convolve
					  //				  (clt_parameters.dbg_mode & 256) != 0, // transpose convolve
					  threadsMax,
					  debugLevel);
		  }

		  double [][] clt_mismatch = new double [12][tilesX*tilesY];
		  for (int pair = 0; pair < 4; pair++){
			  for (int scan_step = 0; scan_step < clt_parameters.disp_scan_count; scan_step++){
				  for (int i = 0; i < tilesX * tilesY; i++){
					  double w = clt_mismatches[scan_step][3 * pair + 2][i];
					  clt_mismatch[3 * pair + 0][i] += clt_mismatches[scan_step][3 * pair + 0][i] * w;
					  clt_mismatch[3 * pair + 1][i] += clt_mismatches[scan_step][3 * pair + 1][i] * w;
					  clt_mismatch[3 * pair + 2][i] += w;
				  }
			  }
		  }
		  for (int pair = 0; pair < 4; pair++){
			  for (int i = 0; i < tilesX * tilesY; i++){
				  double w = clt_mismatch[3 * pair + 2][i];
				  if (w != 0.0){
					  clt_mismatch[3 * pair + 0][i] /= w;
					  clt_mismatch[3 * pair + 1][i] /= w;
				  }
			  }
		  }

		  if (clt_mismatches != null) { // now always
			    AlignmentCorrection ac = new AlignmentCorrection(this);
			  double [][] scans = ac.combineCltMismatches(
					  clt_parameters, // EyesisCorrectionParameters.CLTParameters clt_parameters,
					  clt_mismatches, // double [][][]                            clt_mismatches,
					  disparity_maps, // double [][][]                            disparity_maps,
					  ImageDtt.DISPARITY_INDEX_CM, // int                                      disparity_index,
					  ImageDtt.DISPARITY_STRENGTH_INDEX); // int                                      strength_index)

			  if (clt_parameters.fine_dbg) {
				  ac.showCltMismatches(
						  "clt_mismatches", // String                                   title,
						  clt_parameters,   // EyesisCorrectionParameters.CLTParameters clt_parameters,
						  scans, // double [][]                              combo_data,
						  tp.getTilesX(),
						  tp.getTilesY());
				  if (debugLevel > 1) {
					  ac.showCltMismatch(
							  "clt_mismatch", // String                                   title,
							  clt_parameters,   // EyesisCorrectionParameters.CLTParameters clt_parameters,
							  clt_mismatch, // double [][][]                            clt_mismatch)
							  tp.getTilesX(),
							  tp.getTilesY());
				  }
			  }
			  // TODO: Add automatic run of the lazy eye here
			  if (true) {
				  System.out.println("=== Calculating lazy eye correction ===");
//				  final int nTiles =tilesX * tilesY;
				  double [][] disp_strength = ac.getFineCorrFromDoubleArray(
						  scans,  // double [][] data,
						  tilesX, // int         tilesX,
						  debugLevel); // int debugLevel)

				  double [][][] new_corr = ac.lazyEyeCorrection(
						  clt_parameters.ly_poly,        // final boolean use_poly,
						  true, // 	final boolean    restore_disp_inf, // Restore subtracted disparity for scan #0 (infinity)
						  clt_parameters.fcorr_radius,       // 	final double fcorr_radius,
						  clt_parameters.fcorr_inf_strength, //  final double min_strenth,
						  clt_parameters.fcorr_inf_diff,     // final double max_diff,
						  //								1.3, // final double comp_strength_var,
						  clt_parameters.inf_iters,          // 20, // 0, // final int max_iterations,
						  clt_parameters.inf_final_diff,     // 0.0001, // final double max_coeff_diff,
						  clt_parameters.inf_far_pull,       // 0.0, // 0.25, //   final double far_pull, //  = 0.2; // 1; //  0.5;
						  clt_parameters.inf_str_pow,        // 1.0, //   final double     strength_pow,
						  0.8*clt_parameters.disp_scan_step, // 1.5, // final double     lazyEyeCompDiff, // clt_parameters.fcorr_disp_diff
						  clt_parameters.ly_smpl_side,       // 3,   // final int        lazyEyeSmplSide, //        = 2;      // Sample size (side of a square)
						  clt_parameters.ly_smpl_num,        // 5,   // final int        lazyEyeSmplNum, //         = 3;      // Number after removing worst (should be >1)
						  clt_parameters.ly_smpl_rms,        // 0.1, // final double     lazyEyeSmplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
						  clt_parameters.ly_disp_var,        // 0.2, // final double     lazyEyeDispVariation, // 0.2, maximal full disparity difference between tgh tile and 8 neighborxs
						  clt_parameters.ly_disp_rvar,       // 0.2, // final double     lazyEyeDispRelVariation, // 0.02 Maximal relative full disparity difference to 8 neighbors
						  clt_parameters.ly_norm_disp,       //	final double     ly_norm_disp, //  =    5.0;     // Reduce weight of higher disparity tiles
						  clt_parameters.inf_smpl_side,      // 3, //   final int        smplSide, //        = 2;      // Sample size (side of a square)
						  clt_parameters.inf_smpl_num,       // 5, //   final int        smplNum,  //         = 3;      // Number after removing worst (should be >1)
						  clt_parameters.inf_smpl_rms,       // 0.1, // 0.05, //  final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
						  // histogram parameters
						  clt_parameters.ih_smpl_step,       // 8,    // final int        hist_smpl_side, // 8 x8 masked, 16x16 sampled
						  clt_parameters.ih_disp_min,        // -1.0, // final double     hist_disp_min,
						  clt_parameters.ih_disp_step,       // 0.05, // final double     hist_disp_step,
						  clt_parameters.ih_num_bins,        // 40,   // final int        hist_num_bins,
						  clt_parameters.ih_sigma,           // 0.1,  // final double     hist_sigma,
						  clt_parameters.ih_max_diff,        // 0.1,  // final double     hist_max_diff,
						  clt_parameters.ih_min_samples,     // 10,   // final int        hist_min_samples,
						  clt_parameters.ih_norm_center,     // true, // final boolean    hist_norm_center, // if there are more tiles that fit than min_samples, replace with
						  clt_parameters.ly_inf_frac,        // 0.5, // final double     inf_fraction,    // fraction of the weight for the infinity tiles
						  clt_parameters,  // EyesisCorrectionParameters.CLTParameters           clt_parameters,
						  disp_strength, // scans,   // double [][] disp_strength,
						  null,          // double [][]      target_disparity, // null or programmed disparity (1 per each 14 entries of scans_14)
						  tilesX, // int         tilesX,
						  clt_parameters.corr_magic_scale, // double      magic_coeff, // still not understood coefficent that reduces reported disparity value.  Seems to be around 8.5
						  debugLevel + (clt_parameters.fine_dbg ? 1:0)); // int debugLevel)
				  if (clt_parameters.ly_on_scan && clt_parameters.ly_poly && (new_corr != null)) {
					  if (debugLevel > -1){
						  System.out.println("=== Applying lazy eye correction === ");
					  }

					  apply_fine_corr(
							  new_corr,
							  debugLevel + 2);
				  } else {
					  if (debugLevel > -1){
						  System.out.println("=== SKIPPING application of the lazy eye correction because \"ly_on_scan\" is not set=== ");
					  }
				  }
			  }
		  }
		  int [] disp_indices = {
				  ImageDtt.DISPARITY_INDEX_INT,
				  ImageDtt.DISPARITY_INDEX_CM,
				  ImageDtt.DISPARITY_INDEX_POLY,
				  ImageDtt.DISPARITY_STRENGTH_INDEX,
				  ImageDtt.DISPARITY_VARIATIONS_INDEX};
		  String [] disparity_titles = new String [disp_indices.length];
		  for (int i = 0; i < disparity_titles.length; i++ ) disparity_titles[i] = ImageDtt.DISPARITY_TITLES[i];

		  //				  2,4,6,7};
		  String [] disparities_titles = new String [disparity_titles.length * clt_parameters.disp_scan_count];
		  double [][] disparities_maps = new double [disparity_titles.length * clt_parameters.disp_scan_count][];
		  int indx = 0;

		  for (int scan_step = 0; scan_step < clt_parameters.disp_scan_count; scan_step++) {
			  double disparity = clt_parameters.disp_scan_start + scan_step * clt_parameters.disp_scan_step;
			  for (int i = 0; i < disparity_titles.length; i++){
				  disparities_titles[indx] = disparity_titles[i]+"_"+disparity;
				  disparities_maps[indx++] = disparity_maps[scan_step][disp_indices[i]];
			  }
		  }

		  if (!clt_parameters.ly_on_scan) { // do not show if scan ran for the lazy eye correction
			  ImageStack array_stack = sdfa_instance.makeStack(
					  disparities_maps,
					  tilesX,
					  tilesY,
					  disparities_titles);

			  ImagePlus imp_stack = new ImagePlus( name+"-DISP_MAPS", array_stack);
			  imp_stack.getProcessor().resetMinAndMax();
			  imp_stack.updateAndDraw();
			  //imp_stack.getProcessor().resetMinAndMax();
			  //imp_stack.show();
			  eyesisCorrections.saveAndShow(imp_stack, this.correctionsParameters);
			  // process scan results
			  double [][] scan_trends = process_disparity_scan(
					  disparities_maps,
					  clt_parameters.disp_scan_step,
					  clt_parameters.disp_scan_start,
					  0.0); // min_corr_selected); // all what is remaining
			  String [] trend_titles={"b_cm","b_poly","a_cm","a_poly","rms_cm","rms_poly","strength","samples"};

			  ImageStack trends_stack = sdfa_instance.makeStack(
					  scan_trends,
					  tilesX,
					  tilesY,
					  trend_titles);
			  ImagePlus imp_stack_trends = new ImagePlus( name+"-DISP_TRENDS", trends_stack);
			  imp_stack_trends.getProcessor().resetMinAndMax();
			  imp_stack_trends.updateAndDraw();
			  eyesisCorrections.saveAndShow(imp_stack_trends, this.correctionsParameters);
		  }
		  return results;
	  }

	  public void process_infinity_corr( //from existing image
			  EyesisCorrectionParameters.CLTParameters clt_parameters,
			  int debugLevel
			  ) {
	        ImagePlus imp_src = WindowManager.getCurrentImage();
	        if (imp_src==null){
	            IJ.showMessage("Error","14*n-layer file with disparities/strengthspairs measured at infinity is required");
	            return;
	        }
	        ImageStack disp_strength_stack= imp_src.getStack();
		    final int tilesX = disp_strength_stack.getWidth(); // tp.getTilesX();
		    final int tilesY = disp_strength_stack.getHeight(); // tp.getTilesY();
		    final int nTiles =tilesX * tilesY;
		    final int num_scans =  disp_strength_stack.getSize()/AlignmentCorrection.NUM_ALL_SLICES;
		    final double [][] inf_disp_strength = new double [AlignmentCorrection.NUM_SLICES * num_scans][nTiles];
		    for (int n = 0; n <  disp_strength_stack.getSize(); n++){
	    		float [] fpixels = (float[]) disp_strength_stack.getPixels(n +1);
	    		for (int i = 0; i < nTiles; i++){
	    			inf_disp_strength[n][i] = fpixels[i];
	    		}
		    }
		    if (debugLevel > -1){
		    	System.out.println("process_infinity_corr(): proocessing "+num_scans+" disparity/strength pairs");
		    }
		    AlignmentCorrection ac = new AlignmentCorrection(this);
		    // includes both infinity correction and mismatch correction for the same infinity tiles

		    //FIXME: Here disparity now should be restored in dxy...


		    double [][][] new_corr = ac.infinityCorrection(
		    		clt_parameters.ly_poly,        // final boolean use_poly,
		    		clt_parameters.fcorr_inf_strength, //  final double min_strenth,
		    		clt_parameters.fcorr_inf_diff,     // final double max_diff,
		    		clt_parameters.inf_iters,          // 20, // 0, // final int max_iterations,
		    		clt_parameters.inf_final_diff,     // 0.0001, // final double max_coeff_diff,
		    		clt_parameters.inf_far_pull,       // 0.0, // 0.25, //   final double far_pull, //  = 0.2; // 1; //  0.5;
		    		clt_parameters.inf_str_pow,        // 1.0, //   final double     strength_pow,
		    		clt_parameters.inf_smpl_side,      // 3, //   final int        smplSide, //        = 2;      // Sample size (side of a square)
		    		clt_parameters.inf_smpl_num,       // 5, //   final int        smplNum,  //         = 3;      // Number after removing worst (should be >1)
		    		clt_parameters.inf_smpl_rms,       // 0.1, // 0.05, //  final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
					// histogram parameters
		    		clt_parameters.ih_smpl_step,       // 8,    // final int        hist_smpl_side, // 8 x8 masked, 16x16 sampled
		    		clt_parameters.ih_disp_min,        // -1.0, // final double     hist_disp_min,
		    		clt_parameters.ih_disp_step,       // 0.05, // final double     hist_disp_step,
		    		clt_parameters.ih_num_bins,        // 40,   // final int        hist_num_bins,
		    		clt_parameters.ih_sigma,           // 0.1,  // final double     hist_sigma,
		    		clt_parameters.ih_max_diff,        // 0.1,  // final double     hist_max_diff,
		    		clt_parameters.ih_min_samples,     // 10,   // final int        hist_min_samples,
		    		clt_parameters.ih_norm_center,     // true, // final boolean    hist_norm_center, // if there are more tiles that fit than min_samples, replace with
					clt_parameters,  // EyesisCorrectionParameters.CLTParameters           clt_parameters,
		    		inf_disp_strength,   // double [][] disp_strength,
		    		tilesX, // int         tilesX,
		    		clt_parameters.corr_magic_scale, // double      magic_coeff, // still not understood coefficent that reduces reported disparity value.  Seems to be around 8.5
		    		debugLevel + 1); // int debugLevel)




		    if (debugLevel > -1){
		    	System.out.println("process_infinity_corr(): ready to apply infinity correction");
		    	show_fine_corr(
		    			new_corr, // double [][][] corr,
		    			"");// String prefix)

		    }

		    if (debugLevel > -100){
		    	apply_fine_corr(
		    			new_corr,
		    			debugLevel + 2);
		    }
	  }


	  public void processLazyEye(
			  boolean dry_run,
			  EyesisCorrectionParameters.CLTParameters clt_parameters,
			  int debugLevel
			  ) {
	        ImagePlus imp_src = WindowManager.getCurrentImage();
	        if (imp_src==null){
	            IJ.showMessage("Error","2*n-layer file with disparities/strengthspairs measured at infinity is required");
	            return;
	        }
	        ImageStack disp_strength_stack= imp_src.getStack();
		    final int tilesX = disp_strength_stack.getWidth(); // tp.getTilesX();

		    AlignmentCorrection ac = new AlignmentCorrection(this);



			double [][] scans = ac.getDoubleFromImage(
					imp_src,
					debugLevel);
			double [][] disp_strength = ac.getFineCorrFromDoubleArray(
					scans,  // double [][] data,
					tilesX, // int         tilesX,
					debugLevel); // int debugLevel)


		    double [][][] new_corr = ac.lazyEyeCorrection(
		    		clt_parameters.ly_poly,        // final boolean use_poly,
					true, // final boolean    restore_disp_inf, // Restore subtracted disparity for scan #0 (infinity)
				    clt_parameters.fcorr_radius,       // 	final double fcorr_radius,
		    		clt_parameters.fcorr_inf_strength, //  final double min_strenth,
		    		clt_parameters.fcorr_inf_diff,     // final double max_diff,
		    		clt_parameters.inf_iters,          // 20, // 0, // final int max_iterations,
		    		clt_parameters.inf_final_diff,     // 0.0001, // final double max_coeff_diff,
		    		clt_parameters.inf_far_pull,       // 0.0, // 0.25, //   final double far_pull, //  = 0.2; // 1; //  0.5;
		    		clt_parameters.inf_str_pow,        // 1.0, //   final double     strength_pow,
		    		0.8*clt_parameters.disp_scan_step, // clt_parameters.ly_meas_disp,       // 1.5, // final double     lazyEyeCompDiff, // clt_parameters.fcorr_disp_diff
					clt_parameters.ly_smpl_side,       // 3,   // final int        lazyEyeSmplSide, //        = 2;      // Sample size (side of a square)
					clt_parameters.ly_smpl_num,        // 5,   // final int        lazyEyeSmplNum, //         = 3;      // Number after removing worst (should be >1)
					clt_parameters.ly_smpl_rms,        // 0.1, // final double     lazyEyeSmplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
					clt_parameters.ly_disp_var,        // 0.2, // final double     lazyEyeDispVariation, // 0.2, maximal full disparity difference between tgh tile and 8 neighbors
					clt_parameters.ly_disp_rvar,       // 0.2, // final double     lazyEyeDispRelVariation, // 0.02 Maximal relative full disparity difference to 8 neighbors
					clt_parameters.ly_norm_disp,       //	final double     ly_norm_disp, //  =    5.0;     // Reduce weight of higher disparity tiles
		    		clt_parameters.inf_smpl_side,      // 3, //   final int        smplSide, //        = 2;      // Sample size (side of a square)
		    		clt_parameters.inf_smpl_num,       // 5, //   final int        smplNum,  //         = 3;      // Number after removing worst (should be >1)
		    		clt_parameters.inf_smpl_rms,       // 0.1, // 0.05, //  final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
					// histogram parameters
		    		clt_parameters.ih_smpl_step,       // 8,    // final int        hist_smpl_side, // 8 x8 masked, 16x16 sampled
		    		clt_parameters.ih_disp_min,        // -1.0, // final double     hist_disp_min,
		    		clt_parameters.ih_disp_step,       // 0.05, // final double     hist_disp_step,
		    		clt_parameters.ih_num_bins,        // 40,   // final int        hist_num_bins,
		    		clt_parameters.ih_sigma,           // 0.1,  // final double     hist_sigma,
		    		clt_parameters.ih_max_diff,        // 0.1,  // final double     hist_max_diff,
		    		clt_parameters.ih_min_samples,     // 10,   // final int        hist_min_samples,
		    		clt_parameters.ih_norm_center,     // true, // final boolean    hist_norm_center, // if there are more tiles that fit than min_samples, replace with
		    		clt_parameters.ly_inf_frac,        // 0.5, // final double     inf_fraction,    // fraction of the weight for the infinity tiles
					clt_parameters,  // EyesisCorrectionParameters.CLTParameters           clt_parameters,
					disp_strength, // scans,   // double [][] disp_strength,
					null,          // double [][]      target_disparity, // null or programmed disparity (1 per each 14 entries of scans_14)
		    		tilesX, // int         tilesX,
		    		clt_parameters.corr_magic_scale, // double      magic_coeff, // still not understood coefficent that reduces reported disparity value.  Seems to be around 8.5
		    		debugLevel + (clt_parameters.fine_dbg ? 1:0)); // int debugLevel)

		    if (!dry_run  && clt_parameters.ly_poly && (new_corr != null)){
				  apply_fine_corr(
						  new_corr,
						  debugLevel + 2);

		    }
	  }


	  public double [][] process_disparity_scan(
			  double [][] disparities_maps,
			  double disp_step,
			  double disp_start,
			  double min_strength)
	  {
		  final int num_items = 5;
		  final int index_strength = 3;
		  final int index_cm = 1;
		  final int index_poly = 1;

		  final int ind_b_cm =   0;
		  final int ind_b_poly = 1;
		  final int ind_a_cm =   2;
		  final int ind_a_poly = 3;
		  final int ind_rms_cm =   4;
		  final int ind_rms_poly = 5;
		  final int ind_strength =   6;
		  final int ind_samples = 7;

		  final int num_steps = disparities_maps.length / num_items; // int, cm, poly, strength, variety
		  final int disp_len =  disparities_maps[0].length;
		  double [][] rslt = new double [8][disp_len];
		  for (int i = 0; i < disp_len; i++){
			  double s0 = 0.0, sx = 0.0, sx2 = 0.0, sy_cm = 0.0, sxy_cm = 0.0, sy_poly = 0.0, sxy_poly = 0.0;
			  int samples = 0;
			  for (int step = 0; step < num_steps; step++ ){
				  double wi = disparities_maps[num_items*step + index_strength][i];
				  if (wi > min_strength) {
					  double xi =      disp_start + step * disp_step;
					  double yi_cm =   disparities_maps[num_items*step + index_cm][i];
					  double yi_poly = disparities_maps[num_items*step + index_poly][i];
					  if (Double.isNaN(yi_cm) || Double.isNaN(yi_poly)) continue;
					  s0 +=     wi;
					  sx +=     wi*xi;
					  sx2 +=    wi*xi*xi;
					  sy_cm +=  wi*yi_cm;
					  sxy_cm += wi*xi*yi_cm;
					  sy_poly +=  wi*yi_poly;
					  sxy_poly += wi*xi*yi_poly;
					  samples++;
				  }
			  }
			  double denom = (s0*sx2 - sx*sx);
			  rslt[ind_strength][i] = s0;
			  rslt[ind_samples][i] =  samples;
			  if (denom != 0.0) {
				  rslt[ind_a_cm][i] =     (s0*sxy_cm - sx*sy_cm) /  denom;
				  rslt[ind_b_cm][i] =     (sy_cm*sx2 - sx*sxy_cm) / denom;
				  rslt[ind_a_poly][i] =   (s0*sxy_poly - sx*sy_poly) /  denom;
				  rslt[ind_b_poly][i] =   (sy_poly*sx2 - sx*sxy_poly) / denom;
				  rslt[ind_rms_cm][i] =   0.0;
				  rslt[ind_rms_poly][i] = 0.0;
				  for (int step = 0; step < num_steps; step++ ){
					  double wi = disparities_maps[num_items*step + index_strength][i];
					  if (wi > min_strength) {
						  double xi =      disp_start + step * disp_step;
						  double yi_cm =   disparities_maps[num_items*step + index_cm][i];
						  double yi_poly = disparities_maps[num_items*step + index_poly][i];
						  if (Double.isNaN(yi_cm) || Double.isNaN(yi_poly)) continue;
						  double d_cm =   yi_cm -   (rslt[ind_a_cm][i]*xi +rslt[ind_b_cm][i]);
						  double d_poly = yi_poly - (rslt[ind_a_poly][i]*xi +rslt[ind_b_poly][i]);
						  rslt[ind_rms_cm][i] +=   wi*d_cm*d_cm;
						  rslt[ind_rms_poly][i] += wi*d_poly*d_poly;
					  }
				  }
				  rslt[ind_rms_cm][i] =   Math.sqrt(rslt[ind_rms_cm][i]/s0);
				  rslt[ind_rms_poly][i] = Math.sqrt(rslt[ind_rms_cm][i]/s0);
			  } else {
				  rslt[ind_a_cm][i] =     Double.NaN;
				  rslt[ind_b_cm][i] =     Double.NaN;
				  rslt[ind_a_poly][i] =   Double.NaN;
				  rslt[ind_b_poly][i] =   Double.NaN;
				  rslt[ind_rms_cm][i] =   Double.NaN;
				  rslt[ind_rms_poly][i] = Double.NaN;
			  }

		  }
		  return rslt;
	  }

	  public void showCLTPlanes(
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  final int          threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel)
	  {
		  this.startStepTime=System.nanoTime();

		  if (tp == null){
			  System.out.println("showCLTPlanes(): tp is null");
			  return;
		  }
		  if (tp.clt_3d_passes == null){
			  System.out.println("showCLTPlanes(): tp.clt_3d_passes is null");
			  return;
		  }
		  tp.showPlanes(
				  clt_parameters,
				  geometryCorrection,
				  threadsMax,
				  updateStatus,
//				  false, // batch_mode
				  debugLevel);
		  Runtime.getRuntime().gc();
	      System.out.println("showCLTPlanes(): processing  finished at "+
			  IJ.d2s(0.000000001*(System.nanoTime()-this.startStepTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
	  }

	  public double [][]  assignCLTPlanes(
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  final int          threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel)
	  {
		  	if (tp == null){
		  		System.out.println("showCLTPlanes(): tp is null");
		  		return null;
		  	}
		  	if (tp.clt_3d_passes == null){
		  		System.out.println("showCLTPlanes(): tp.clt_3d_passes is null");
		  		return null;
		  	}
		  	this.startStepTime=System.nanoTime();
		  	double [][] assign_dbg = tp.assignTilesToSurfaces(
		  			clt_parameters,
		  			geometryCorrection,
		  			threadsMax,
		  			updateStatus,
//		  			false, //  boolean batch_mode,
		  			debugLevel);
		  	Runtime.getRuntime().gc();
		  	System.out.println("assignCLTPlanes(): processing  finished at "+
		  			IJ.d2s(0.000000001*(System.nanoTime()-this.startStepTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
		  	return assign_dbg;

	  }


	  public void out3d_old(
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  final int          threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel)
	  {
		  	if (tp == null){
		  		System.out.println("showCLTPlanes(): tp is null");
		  		return;
		  	}
		  	if (tp.clt_3d_passes == null){
		  		System.out.println("showCLTPlanes(): tp.clt_3d_passes is null");
		  		return;
		  	}
		  	tp.showPlanes(
		  			clt_parameters,
		  			geometryCorrection,
		  			threadsMax,
		  			updateStatus,
//		  			false, // batch_mode
		  			debugLevel);
//		  	CLTPass3d last_scan = tp.clt_3d_passes.get(tp.clt_3d_passes.size() -1); // get last one

	  }



	  public void processCLTQuads3d(
			  boolean adjust_extrinsics,
			  boolean adjust_poly,
			  TwoQuadCLT       twoQuadCLT, //maybe null in no-rig mode, otherwise may contain rig measurements to be used as infinity ground truth
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
			  final int          threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel)
	  {
		  this.startTime=System.nanoTime();
		  String [] sourceFiles=correctionsParameters.getSourcePaths();
		  SetChannels [] set_channels=setChannels(debugLevel);
		  if ((set_channels == null) || (set_channels.length==0)) {
			  System.out.println("No files to process (of "+sourceFiles.length+")");
			  return;
		  }
		  double [] referenceExposures=eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		  for (int nSet = 0; nSet < set_channels.length; nSet++){
			  int [] channelFiles = set_channels[nSet].fileNumber();
			  boolean [][] saturation_imp = (clt_parameters.sat_level > 0.0)? new boolean[channelFiles.length][] : null;
			  double [] scaleExposures = new double[channelFiles.length];

			  ImagePlus [] imp_srcs = conditionImageSet(
					  clt_parameters,             // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					  sourceFiles,                // String []                                 sourceFiles,
					  set_channels[nSet].name(),  // String                                    set_name,
					  referenceExposures,         // double []                                 referenceExposures,
					  channelFiles,               // int []                                    channelFiles,
					  scaleExposures,   //output  // double [] scaleExposures
					  saturation_imp,   //output  // boolean [][]                              saturation_imp,
					  debugLevel); // int                                       debugLevel);
//			  if (adjust_extrinsics && (debugLevel >-2)) {
//				  boolean tmp_exit = (debugLevel > -10); // == true;
//				  System.out.println("processCLTQuads3d(): adjust_extrinsics="+adjust_extrinsics);
//				  if (tmp_exit) {
//					  System.out.println("will now exit. To continue - change variable tmp_exit in debugger" );
//					  if (tmp_exit) {
//						  return;
//					  }
//				  }
//			  }
			  boolean use_rig = (twoQuadCLT != null) && (twoQuadCLT.getBiScan(0) != null);
			  if (!adjust_extrinsics || !use_rig) {
				  // Difficult to fix: adjust extrinsics for aux - when it is adjusted alone, it will not match tiles to those of a rig!
				  // can use only far tiles with small gradients?

				  // once per quad here
				  preExpandCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
						  imp_srcs, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
						  saturation_imp, // boolean [][] saturation_imp, // (near) saturated pixels or null
						  clt_parameters,
						  debayerParameters,
						  colorProcParameters,
						  rgbParameters,
						  threadsMax,  // maximal number of threads to launch
						  updateStatus,
						  debugLevel);

				  // adjust extrinsics here
			  }
			  if (adjust_extrinsics) {
				  if (use_rig) {
					  System.out.println("Adjust extrinsics using rig data here");
					  extrinsicsCLTfromGT(
							  twoQuadCLT,   // TwoQuadCLT       twoQuadCLT, //maybe null in no-rig mode, otherwise may contain rig measurements to be used as infinity ground truth
							  clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
							  adjust_poly,
							  threadsMax,  //final int        threadsMax,  // maximal number of threads to launch
							  updateStatus,// final boolean    updateStatus,
							  debugLevel + 2); // final int        debugLevel)

				  } else {
					  System.out.println("Adjust extrinsics here");
					  extrinsicsCLT(
							  //						  twoQuadCLT,   // TwoQuadCLT       twoQuadCLT, //maybe null in no-rig mode, otherwise may contain rig measurements to be used as infinity ground truth
							  clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
							  adjust_poly,
							  threadsMax,  //final int        threadsMax,  // maximal number of threads to launch
							  updateStatus,// final boolean    updateStatus,
							  debugLevel); // final int        debugLevel)
				  }

			  } else {
				  expandCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
						  clt_parameters,
						  debayerParameters,
						  colorProcParameters,
						  channelGainParameters,
						  rgbParameters,
						  threadsMax,  // maximal number of threads to launch
						  updateStatus,
						  debugLevel);
			  }

			  Runtime.getRuntime().gc();

//			  if (debugLevel >-1) System.out.println("Processing set "+(nSet+1)+" (of "+setNames.size()+") finished at "+
//					  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
			  if (debugLevel >-1) System.out.println("Processing set "+(nSet+1)+" (of "+set_channels.length+") finished at "+
					  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
			  if (eyesisCorrections.stopRequested.get()>0) {
				  System.out.println("User requested stop");
				  return;
			  }
		  }
//		  System.out.println("Processing "+fileIndices.length+" files finished at "+
//				  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
		  System.out.println("Processing "+getTotalFiles(set_channels)+" files finished at "+
				  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
	  }

	  public boolean preExpandCLTQuad3d(
			  ImagePlus []                                     imp_quad, // should have properties "name"(base for saving results), "channel","path"
			  boolean [][]                                     saturation_imp,   // (near) saturated pixels or null
			  EyesisCorrectionParameters.CLTParameters         clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  EyesisCorrectionParameters.ColorProcParameters   colorProcParameters,
			  EyesisCorrectionParameters.RGBParameters         rgbParameters,
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel)
	  {
		  this.startStepTime=System.nanoTime();
		  final boolean    show_init_refine = clt_parameters.show_init_refine;

		  //max_expand
		  String name = (String) imp_quad[0].getProperty("name");
		  double [][][] image_data = new double [imp_quad.length][][];
		  for (int i = 0; i < image_data.length; i++){
			  image_data[i] = eyesisCorrections.bayerToDoubleStack(
					  imp_quad[i], // source Bayer image, linearized, 32-bit (float))
					  null); // no margins, no oversample
		  }
		  for (int i = 0; i < image_data.length; i++){
			  for (int j =0 ; j < image_data[i][0].length; j++){
				  image_data[i][2][j]*=0.5; // Scale green 0.5 to compensate more pixels than R,B
			  }
		  }
		  setTiles (imp_quad[0], // set global tp.tilesX, tp.tilesY
				  clt_parameters,
				  threadsMax);

		  this.image_name =     name;
		  this.image_data =     image_data;
		  this.saturation_imp = saturation_imp;


		  tp.resetCLTPasses();
		  tp.setTrustedCorrelation(clt_parameters.grow_disp_trust);

		  CLTPass3d bgnd_data = CLTBackgroundMeas( // measure background
				  image_data, //
				  saturation_imp,    // boolean [][] saturation_imp, // (near) saturated pixels or null
				  clt_parameters,
				  threadsMax,  // maximal number of threads to launch
				  updateStatus,
				  debugLevel);
		  tp.clt_3d_passes.add(bgnd_data);
		  //    	  if (show_init_refine)
		  if ((debugLevel > -2) && clt_parameters.show_first_bg) {
			  tp.showScan(
					  tp.clt_3d_passes.get(0), // CLTPass3d   scan,
					  "bgnd_data-"+tp.clt_3d_passes.size());
		  }

		  //TODO: Move away from here?
		  boolean no_image_save = true;

    	  ImagePlus imp_bgnd = getBackgroundImage(
    			  no_image_save, // boolean no_image,
    			  clt_parameters,
    			  colorProcParameters,
    			  rgbParameters,
    			  name,               // .getTitle(), //String name=(String) imp_src.getProperty("name");
    			  ImageDtt.DISPARITY_INDEX_CM, // index of disparity value in disparity_map == 2 (0,2 or 4)
    			  threadsMax,  // maximal number of threads to launch
    			  updateStatus,
    			  debugLevel);
//		  if (debugLevel > -100) {
//			  return null;
//		  }


    	  bgnd_data.texture = (imp_bgnd == null)? null: ( imp_bgnd.getTitle()+ (clt_parameters.black_back? ".jpeg" : ".png"));

    	  // create x3d file
    	  X3dOutput x3dOutput = new X3dOutput(
    			  clt_parameters,
    			  correctionsParameters,
    			  geometryCorrection,
    			  tp.clt_3d_passes);

		  x3dOutput.generateBackground(clt_parameters.infinityDistance <= 0.0); // needs just first (background) scan

/*		  */


    	  // refine first measurement
    	  int bg_pass = tp.clt_3d_passes.size() - 1; // 0
    	  int refine_pass = tp.clt_3d_passes.size(); // 1

    	  //		  final boolean show_init_refine = true;
    	  //		  final boolean show_expand =      true;

    	  //    	  if (show_init_refine)
    	  if (debugLevel > -1) {
    		  tp.showScan(
    				  tp.clt_3d_passes.get(bg_pass), // CLTPass3d   scan,
    				  "after_bg-"+tp.clt_3d_passes.size());
    	  }

//    	  final double     weight_var = 1.0; // 1.0;   // weight of variance data (old, detects thin wires?)
 //   	  final double     weight_Y =   1.0;     // weight of average intensity
 //   	  final double     weight_RBmG = 5.0;  // weight of average color difference (0.5*(R+B)-G), shoukld be ~5*weight_Y

    	  // TODO: Make double pass - with only 	weight_var (thin wires) and weight_Y, weight_RBmG - larger objects
    	  // just use two instances of MacroCorrelation, run one after another (move code to MacroCorrelation class)
    	  // and then join

    	  MacroCorrelation mc = new MacroCorrelation(
    			  tp,
    			  clt_parameters.mc_disp8_trust,
    			  clt_parameters.mc_weight_var,   // final double     weight_var,   // weight of variance data (old, detects thin wires?)
    			  clt_parameters.mc_weight_Y,     // final double     weight_Y,     // weight of average intensity
    			  clt_parameters.mc_weight_RBmG   // final double     weight_RBmG,  // weight of average color difference (0.5*(R+B)-G), shoukld be ~5*weight_Y
    			  );


    	  double [][][] input_data = mc.CLTMacroSetData( // perform single pass according to prepared tiles operations and disparity
    			  bgnd_data);           // final CLTPass3d      src_scan, // results of the normal correlations (now expecting infinity)
    	  //    			  null);               // final double [][][]  other_channels, // other channels to correlate, such as average RGB (first index - subcamera, 2-nd - channel, 3-rd - pixel)


    	  TileProcessor mtp =  mc.CLTMacroScan( // perform single pass according to prepared tiles operations and disparity
    			  bgnd_data,          // final CLTPass3d                            src_scan, // results of the normal correlations (now expecting infinity)
    			  clt_parameters,     // EyesisCorrectionParameters.CLTParameters clt_parameters,
    			  geometryCorrection, // GeometryCorrection geometryCorrection,
    			  0.0,                // 	final double                             macro_disparity_low,
    			  clt_parameters.grow_disp_max / tp.getTileSize(), // final double                             macro_disparity_high,
    			  clt_parameters.mc_disp8_step, // final double                             macro_disparity_step,
    			  debugLevel); // 1); // 	final int                                debugLevel){

    	  CLTPass3d macro_combo = mtp.compositeScan(
    			  mtp.clt_3d_passes, // final ArrayList <CLTPass3d> passes,
    			  0, //  final int                   firstPass,
    			  mtp.clt_3d_passes.size(),                         //  final int                   lastPassPlus1,
    			  mtp.getTrustedCorrelation(),                      // 	final double                trustedCorrelation,
    			  mtp.getMaxOverexposure(),                         //  final double                max_overexposure,
    			  0.0, // clt_parameters.bgnd_range,                //	final double                disp_far,   // limit results to the disparity range
    			  clt_parameters.grow_disp_max / tp.getTileSize(),  //  final double                disp_near,
    			  clt_parameters.mc_strength,                       // final double                minStrength,
    			  Double.NaN, // clt_parameters.combine_min_hor,                   // final double                minStrengthHor,
    			  Double.NaN, // clt_parameters.combine_min_vert,                  // final double                minStrengthVert,
    			  // maybe temporarily? later keep weak?
    			  true,      // final boolean               no_weak,
    			  false, // final boolean               use_last,   //
    			  // TODO: when useCombo - pay attention to borders (disregard)
    			  false, // final boolean               usePoly)  // use polynomial method to find max), valid if useCombo == false
    			  true, // 	 final boolean               copyDebug)
    			  debugLevel);
    	  mtp.clt_3d_passes.add(macro_combo);
    	  if (clt_parameters.show_macro) {
    		  mtp.showScan(
    				  macro_combo, // CLTPass3d   scan,
    				  "macro_combo-"+mtp.clt_3d_passes.size());
    	  }

    	  for (int num_try = 0; num_try < 100; num_try++) {
    		  CLTPass3d refined_macro = mc.refineMacro(
    				  input_data,                                        // final double [][][]                      input_data,
    				  clt_parameters,                                    // EyesisCorrectionParameters.CLTParameters clt_parameters,
    				  geometryCorrection,                                // GeometryCorrection                       geometryCorrection,
    				  clt_parameters.mc_disp8_trust,                     // final double                             trustedCorrelation,
    				  0, //  final double                             disp_far,   // limit results to the disparity range
    				  clt_parameters.grow_disp_max / tp.getTileSize(), // final double                             disp_near,
    				  clt_parameters.mc_strength, // final double                             minStrength,
    				  clt_parameters.mc_unique_tol, // final double                             unique_tolerance,
    				  1); // final int                                debugLevel)
    		  if (refined_macro == null) break;
    		  mtp.clt_3d_passes.add(refined_macro);
    		  if (clt_parameters.show_macro) {
    			  mtp.showScan(
    					  refined_macro, // CLTPass3d   scan,
    					  "refined_macro-"+mtp.clt_3d_passes.size());
    		  }

    	  }
    	  CLTPass3d macro_combo1 = mtp.compositeScan(
    			  mtp.clt_3d_passes, // final ArrayList <CLTPass3d> passes,
    			  0, //  final int                   firstPass,
    			  mtp.clt_3d_passes.size(),                         //  final int                   lastPassPlus1,
    			  mtp.getTrustedCorrelation(),                       // 	 final double                trustedCorrelation,
    			  mtp.getMaxOverexposure(),                         //  final double                max_overexposure,
    			  0.0, // clt_parameters.bgnd_range,                //	 final double                disp_far,   // limit results to the disparity range
    			  clt_parameters.grow_disp_max / tp.getTileSize(),  //  final double                disp_near,
    			  clt_parameters.mc_strength,                       // final double                minStrength,
    			  Double.NaN, // clt_parameters.combine_min_hor,                   // final double                minStrengthHor,
    			  Double.NaN, // clt_parameters.combine_min_vert,                  // final double                minStrengthVert,
    			  // maybe temporarily? later keep weak?
    			  true,      // final boolean               no_weak,
    			  false, // final boolean               use_last,   //
    			  // TODO: when useCombo - pay attention to borders (disregard)
    			  false, // final boolean               usePoly)  // use polynomial method to find max), valid if useCombo == false
    			  true, // 	 final boolean               copyDebug)
    			  debugLevel);
    	  mtp.clt_3d_passes.add(macro_combo1);
    	  if (clt_parameters.show_macro) {
    		  mtp.showScan(
    			  macro_combo1, // CLTPass3d   scan,
    			  "macro_combo-"+mtp.clt_3d_passes.size());
    	  }


    	  ArrayList <CLTPass3d> new_meas = mc.prepareMeasurementsFromMacro(
    			  mtp.clt_3d_passes, // final ArrayList <CLTPass3d> macro_passes, // macro correlation measurements
    			  // in pixels
    			  3.0, // final double                disp_far,   // limit results to the disparity range
    			  clt_parameters.grow_disp_max, // final double                disp_near,
    			  clt_parameters.mc_strength, // final double                minStrength,
    			  clt_parameters.mc_strength, //final double                mc_trust_fin, //          =   0.3;   // When consolidating macro results, exclude high residual disparity
    			  clt_parameters.mc_strength, //final double                mc_trust_sigma, //        =   0.2;   // Gaussian sigma to reduce weight of large residual disparity
    			  clt_parameters.mc_strength, //final double                mc_ortho_weight, //       =   0.5;   // Weight from ortho neighbor supertiles
    			  clt_parameters.mc_strength, //final double                mc_diag_weight, //        =   0.25;  // Weight from diagonal neighbor supertiles
    			  clt_parameters.mc_strength, //final double                mc_gap, //                =   0.4;   // Do not remove measurements farther from the kept ones
    			  false, // final boolean               usePoly,  // use polynomial method to find max), valid if useCombo == false
    			  true, // final boolean               sort_disparity,  // sort results for increasing disparity (false - decreasing strength)
    			  clt_parameters.tileX, // final int                   dbg_x,
    			  clt_parameters.tileY, // final int                   dbg_y,
    			  debugLevel); // final int                   debugLevel);
    	  System.out.println("new_meas.size()="+new_meas.size());
    	  int indx = 0;
    	  if (clt_parameters.show_macro) {
    		  for (CLTPass3d pass: new_meas) {
    			  tp.showScan(
    					  pass, // CLTPass3d   scan,
    					  "meas-"+(indx++));

    		  }
    	  }

    	  int num_macro_refine = 3;
    	  for (CLTPass3d from_macro_pass: new_meas) {

    		  for (int nnn = 0; nnn < num_macro_refine; nnn ++){ //
    			  refine_pass = tp.clt_3d_passes.size(); // 1
    			  CLTPass3d refined = tp.refinePassSetup( // prepare tile tasks for the refine pass (re-measure disparities)
    					  //				  final double [][][]       image_data, // first index - number of image in a quad
    					  clt_parameters,
    					  clt_parameters.stUseRefine, // use supertiles
    					  bg_pass,
    					  // disparity range - differences from
    					  clt_parameters.bgnd_range, // double            disparity_far,
    					  clt_parameters.grow_disp_max, // other_range, //double            disparity_near,   //
    					  clt_parameters.ex_strength,  // double            this_sure,        // minimal strength to be considered definitely good
    					  clt_parameters.ex_nstrength, // double            ex_nstrength, // minimal 4-corr strength divided by channel diff for new (border) tiles
    					  clt_parameters.bgnd_maybe, // double            this_maybe,       // maximal strength to ignore as non-background
    					  clt_parameters.sure_smth, // sure_smth,        // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
    					  clt_parameters.pt_super_trust, //  final double      super_trust,      // If strength exceeds ex_strength * super_trust, do not apply ex_nstrength and plate_ds
    					  // using plates disparity/strength - averaged for small square sets of tiles. If null - just use raw tiles
    					  null, // 		final double [][] plate_ds,  // disparity/strength last time measured for the multi-tile squares. Strength =-1 - not measured. May be null
    					  true, // 		final boolean     keep_raw_fg,  // do not replace raw tiles by the plates, if raw is closer (like poles)
    					  0.0, // final double      scale_filtered_strength_pre, // scale plate_ds[1] before comparing to raw strength
    					  0.0, // final double      scale_filtered_strength_post,// scale plate_ds[1] when replacing raw (generally plate_ds is more reliable if it exists)

    					  ImageDtt.DISPARITY_INDEX_CM,  // index of disparity value in disparity_map == 2 (0,2 or 4)
    					  geometryCorrection,
    					  threadsMax,  // maximal number of threads to launch
    					  updateStatus,
    					  debugLevel); //2);
    			  tp.clt_3d_passes.add(refined);

    			  ///    		  if (debugLevel > 1)
//    			  if (debugLevel > 0)
       			  if ((debugLevel > -2) && clt_parameters.show_first_bg)
    				  tp.showScan(
    						  tp.clt_3d_passes.get(refine_pass), // CLTPass3d   scan,
    						  "before_makeUnique-"+refine_pass);
    			  int [] numLeftRemoved = tp.makeUnique(
    					  tp.clt_3d_passes,                  // final ArrayList <CLTPass3d> passes,
    					  0,                                 //  final int                   firstPass,
    					  refine_pass, // - 1,                   //  final int                   lastPassPlus1,
    					  tp.clt_3d_passes.get(refine_pass), //  final CLTPass3d             new_scan,
    					  clt_parameters.grow_disp_max,       // final double                grow_disp_max,
    					  clt_parameters.gr_unique_tol,   //  final double                unique_tolerance,
    					  clt_parameters.show_unique);      // final boolean               show_unique)
    			  if (debugLevel > -1){
    				  System.out.println("cycle makeUnique("+refine_pass+") -> left: "+numLeftRemoved[0]+", removed:" + numLeftRemoved[1]);
    			  }
    			  if (show_init_refine) tp.showScan(
    					  tp.clt_3d_passes.get(refine_pass), // CLTPass3d   scan,
    					  "after_refinePassSetup-"+tp.clt_3d_passes.size());

    			  CLTMeasure( // perform single pass according to prepared tiles operations and disparity
    					  image_data, // first index - number of image in a quad
    					  saturation_imp, //final boolean [][]  saturation_imp, // (near) saturated pixels or null
    					  clt_parameters,
    					  refine_pass,
    					  false, // true, // final boolean     save_textures,
    					  threadsMax,  // maximal number of threads to launch
    					  updateStatus,
    					  debugLevel);
    			  if (debugLevel > -1){
    				  System.out.println("CLTMeasure("+refine_pass+")");
    			  }
    			  if (show_init_refine) tp.showScan(
    					  tp.clt_3d_passes.get(refine_pass), // CLTPass3d   scan,
    					  "after_measure-"+tp.clt_3d_passes.size());
    			  if (nnn < (num_macro_refine-1)) {
    				  //        	  if (clt_parameters.combine_refine){
    				  CLTPass3d combo_pass = tp.compositeScan(
    						  tp.clt_3d_passes, // final ArrayList <CLTPass3d> passes,
    						  bg_pass, //  final int                   firstPass,
    						  tp.clt_3d_passes.size(),                          //  final int                   lastPassPlus1,
    						  //      				  tp.clt_3d_passes.get(bg_pass).getSelected(), // selected , // final boolean [] bg_tiles,          // get from selected in clt_3d_passes.get(0);
    						  //   				  clt_parameters.ex_min_over,// final double     ex_min_over,       // when expanding over previously detected (by error) background, disregard far tiles
    						  tp.getTrustedCorrelation(),                       // 	 final double                trustedCorrelation,
    		    			  tp.getMaxOverexposure(),                          //  final double                max_overexposure,
    						  0.0, // clt_parameters.bgnd_range,                //	 final double                disp_far,   // limit results to the disparity range
    						  clt_parameters.grow_disp_max,                     // final double                disp_near,
    						  clt_parameters.combine_min_strength,              // final double                minStrength,
    						  clt_parameters.combine_min_hor,                   // final double                minStrengthHor,
    						  clt_parameters.combine_min_vert,                  // final double                minStrengthVert,
    						  false,      // final boolean               no_weak,
    						  false, // final boolean               use_last,   //
    						  // TODO: when useCombo - pay attention to borders (disregard)
    						  false, // final boolean               usePoly)  // use polynomial method to find max), valid if useCombo == false
    						  true, // 	 final boolean               copyDebug)
    						  debugLevel);

    				  if (show_init_refine) tp.showScan(
    						  combo_pass, // CLTPass3d   scan,
    						  "after_compositeScan-"+tp.clt_3d_passes.size());

    				  tp.clt_3d_passes.add(combo_pass);

    			  }
    		  }
    		  // add new scan from macro
    		  tp.clt_3d_passes.add(from_macro_pass);
			  CLTMeasure( // perform single pass according to prepared tiles operations and disparity
					  image_data, // first index - number of image in a quad
					  saturation_imp, //final boolean [][]  saturation_imp, // (near) saturated pixels or null
					  clt_parameters,
					  tp.clt_3d_passes.size() -1, // refine_pass,
					  true, // final boolean     save_textures,
					  threadsMax,  // maximal number of threads to launch
					  updateStatus,
					  debugLevel);
			  if (debugLevel > -1){
				  System.out.println("CLTMeasure("+refine_pass+")");
			  }
			  if (show_init_refine) tp.showScan(
					  tp.clt_3d_passes.get(refine_pass), // CLTPass3d   scan,
					  "after_measure_macro-"+tp.clt_3d_passes.size());

    	  }

/// Refining after all added
		  if (debugLevel > -1){
			  System.out.println("---- Refining after all added ----");
		  }
    	  // first ("before_makeUnique-41-" was empty)
		  for (int nnn = 0; nnn < num_macro_refine; nnn ++){ //
			  refine_pass = tp.clt_3d_passes.size(); // 1
			  CLTPass3d refined = tp.refinePassSetup( // prepare tile tasks for the refine pass (re-measure disparities)
					  //				  final double [][][]       image_data, // first index - number of image in a quad
					  clt_parameters,
					  clt_parameters.stUseRefine, // use supertiles
					  bg_pass,
					  // disparity range - differences from
					  clt_parameters.bgnd_range, // double            disparity_far,
					  clt_parameters.grow_disp_max, // other_range, //double            disparity_near,   //
					  clt_parameters.ex_strength,  // double            this_sure,        // minimal strength to be considered definitely good
					  clt_parameters.ex_nstrength, // double            ex_nstrength, // minimal 4-corr strength divided by channel diff for new (border) tiles
					  clt_parameters.bgnd_maybe, // double            this_maybe,       // maximal strength to ignore as non-background
					  clt_parameters.sure_smth, // sure_smth,        // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
					  clt_parameters.pt_super_trust, //  final double      super_trust,      // If strength exceeds ex_strength * super_trust, do not apply ex_nstrength and plate_ds
					  // using plates disparity/strength - averaged for small square sets of tiles. If null - just use raw tiles
					  null, // 		final double [][] plate_ds,  // disparity/strength last time measured for the multi-tile squares. Strength =-1 - not measured. May be null
					  true, // 		final boolean     keep_raw_fg,  // do not replace raw tiles by the plates, if raw is closer (like poles)
					  0.0, // final double      scale_filtered_strength_pre, // scale plate_ds[1] before comparing to raw strength
					  0.0, // final double      scale_filtered_strength_post,// scale plate_ds[1] when replacing raw (generally plate_ds is more reliable if it exists)

					  ImageDtt.DISPARITY_INDEX_CM,  // index of disparity value in disparity_map == 2 (0,2 or 4)
					  geometryCorrection,
					  threadsMax,  // maximal number of threads to launch
					  updateStatus,
					  debugLevel); //2);
			  tp.clt_3d_passes.add(refined);

			  ///    		  if (debugLevel > 1)
			  if (debugLevel > 0)
				  tp.showScan(
						  tp.clt_3d_passes.get(refine_pass), // CLTPass3d   scan,
						  "before_makeUnique-"+refine_pass);
			  int [] numLeftRemoved = tp.makeUnique(
					  tp.clt_3d_passes,                  // final ArrayList <CLTPass3d> passes,
					  0,                                 //  final int                   firstPass,
					  refine_pass, // - 1,                   //  final int                   lastPassPlus1,
					  tp.clt_3d_passes.get(refine_pass), //  final CLTPass3d             new_scan,
					  clt_parameters.grow_disp_max,       // final double                grow_disp_max,
					  clt_parameters.gr_unique_tol,   //  final double                unique_tolerance,
					  clt_parameters.show_unique);      // final boolean               show_unique)
			  if (debugLevel > -1){
				  System.out.println("cycle makeUnique("+refine_pass+") -> left: "+numLeftRemoved[0]+", removed:" + numLeftRemoved[1]);
			  }
			  if (show_init_refine) tp.showScan(
					  tp.clt_3d_passes.get(refine_pass), // CLTPass3d   scan,
					  "after_refinePassSetup-"+tp.clt_3d_passes.size());

			  CLTMeasure( // perform single pass according to prepared tiles operations and disparity
					  image_data, // first index - number of image in a quad
					  saturation_imp, //final boolean [][]  saturation_imp, // (near) saturated pixels or null
					  clt_parameters,
					  refine_pass,
					  true, // final boolean     save_textures,
					  threadsMax,  // maximal number of threads to launch
					  updateStatus,
					  debugLevel);
//			  if (debugLevel > -1){
			  if (debugLevel > -2){
				  System.out.println("?.CLTMeasure("+refine_pass+")");
			  }
			  if (show_init_refine) tp.showScan(
					  tp.clt_3d_passes.get(refine_pass), // CLTPass3d   scan,
					  "after_measure-"+tp.clt_3d_passes.size());
//			  if (nnn < (num_macro_refine-1)) {
				  //        	  if (clt_parameters.combine_refine){
				  CLTPass3d combo_pass = tp.compositeScan(
						  tp.clt_3d_passes, // final ArrayList <CLTPass3d> passes,
						  bg_pass, //  final int                   firstPass,
						  tp.clt_3d_passes.size(),                          //  final int                   lastPassPlus1,
						  //      				  tp.clt_3d_passes.get(bg_pass).getSelected(), // selected , // final boolean [] bg_tiles,          // get from selected in clt_3d_passes.get(0);
						  //   				  clt_parameters.ex_min_over,// final double     ex_min_over,       // when expanding over previously detected (by error) background, disregard far tiles
						  tp.getTrustedCorrelation(),                       // 	 final double                trustedCorrelation,
		    			  tp.getMaxOverexposure(),                          //  final double                max_overexposure,
						  0.0, // clt_parameters.bgnd_range,                //	 final double                disp_far,   // limit results to the disparity range
						  clt_parameters.grow_disp_max,                     // final double                disp_near,
						  clt_parameters.combine_min_strength,              // final double                minStrength,
						  clt_parameters.combine_min_hor,                   // final double                minStrengthHor,
						  clt_parameters.combine_min_vert,                  // final double                minStrengthVert,
						  false,      // final boolean               no_weak,
						  false, // final boolean               use_last,   //
						  // TODO: when useCombo - pay attention to borders (disregard)
						  false, // final boolean               usePoly)  // use polynomial method to find max), valid if useCombo == false
						  true, // 	 final boolean               copyDebug)
						  debugLevel);

				  if (show_init_refine) tp.showScan(
						  combo_pass, // CLTPass3d   scan,
						  "after_compositeScan-"+tp.clt_3d_passes.size());

				  tp.clt_3d_passes.add(combo_pass);
		  }

 ///// Refining after all added   - end
		  Runtime.getRuntime().gc();
	      System.out.println("preExpandCLTQuad3d(): processing  finished at "+
			  IJ.d2s(0.000000001*(System.nanoTime()-this.startStepTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
		  return true;
	  }


//	  public ImagePlus expandCLTQuad3d(
	  public boolean extrinsicsCLT(
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  boolean 		   adjust_poly,
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel)
	  {
		  final boolean    batch_mode = clt_parameters.batch_run;
		  int debugLevelInner =  batch_mode ? -5: debugLevel;
		  boolean update_disp_from_latest = clt_parameters.lym_update_disp ; // true;
		  int max_tries =                   clt_parameters.lym_iter; // 25;
		  double min_sym_update =           clt_parameters.lym_change; //  4e-6; // stop iterations if no angle changes more than this
		  double min_poly_update =          clt_parameters.lym_poly_change; //  Parameter vector difference to exit from polynomial correction
		  int bg_scan = 0;
		  int combo_scan= tp.clt_3d_passes.size()-1;


		  if (!batch_mode && clt_parameters.show_extrinsic && (debugLevel >-1)) {
			  //		  if (!batch_mode && (debugLevel >-1)) {
			  tp.showScan(
					  tp.clt_3d_passes.get(bg_scan),   // CLTPass3d   scan,
					  "bg_scan"); //String title)
			  tp.showScan(
					  tp.clt_3d_passes.get(combo_scan),   // CLTPass3d   scan,
					  "combo_scan-"+combo_scan); //String title)
		  }

		  boolean [] bg_sel = null;
		  boolean [] bg_use = null;
		  double [] combo_disp = null;
		  double [] combo_str = null;
		  boolean [] combo_use = null;
		  double [] combo_overexp = null;
		  int num_combo = 0;
		  double [][] filtered_bgnd_disp_strength = tp.getFilteredDisparityStrength(
				  tp.clt_3d_passes, // final ArrayList <CLTPass3d> passes,// List, first, last - to search for the already tried disparity
				  bg_scan,          // final int        measured_scan_index, // will not look at higher scans
				  0,                  //	final int        start_scan_index,
				  null , // final boolean [] bg_tiles,          // get from selected in clt_3d_passes.get(0);
				  0.0,   // whatever as null above //  clt_parameters.ex_min_over,// final double     ex_min_over,       // when expanding over previously detected (by error) background, disregard far tiles

				  ImageDtt.DISPARITY_INDEX_CM,       // final int        disp_index,
				  ImageDtt.DISPARITY_STRENGTH_INDEX, // final int        str_index,
				  null,                    // final double []  tiltXY,    // null - free with limit on both absolute (2.0?) and relative (0.2) values
				  0.5, // clt_parameters.fcorr_inf_diff, // tp.getTrustedCorrelation(),//	final double     trustedCorrelation,
				  clt_parameters.fcorr_inf_strength,     //	final double     strength_floor,
				  clt_parameters.inf_str_pow,       // final double     strength_pow,
				  clt_parameters.ly_smpl_side,           // final int        smplSide, //        = 2;      // Sample size (side of a square)
				  clt_parameters.ly_smpl_num,            // final int        smplNum, //         = 3;      // Number after removing worst (should be >1)
				  clt_parameters.ly_smpl_rms,            //	final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
				  0.0,                                   // smplRelRms,         // final double     smplRelRms, //      = 0.005;  // Maximal RMS/disparity in addition to smplRms
				  clt_parameters.fds_smpl_wnd,            //	final boolean    smplWnd, //
				  clt_parameters.fds_abs_tilt,       //	final double     max_abs_tilt, //  = 2.0; // pix per tile
				  clt_parameters.fds_rel_tilt,       //	final double     max_rel_tilt, //  = 0.2; // (pix / disparity) per tile
				  clt_parameters.tileX, //   dbg_x,              //	final int        dbg_x,
				  clt_parameters.tileX, //  dbg_y,              // final int        dbg_y,
				  debugLevelInner);        //	final int        debugLevel)

		  // prepare re-measurements of background
		  bg_sel = tp.clt_3d_passes.get(bg_scan).getSelected();
		  bg_use = new boolean [bg_sel.length];
		  //		  double  [] bg_disp = tp.clt_3d_passes.get(bg_scan).getDisparity(0);
		  double [] bg_str =  tp.clt_3d_passes.get(bg_scan).getStrength();
		  double [] bg_overexp = tp.clt_3d_passes.get(bg_scan).getOverexposedFraction();
		  for (int  nTile = 0 ; nTile < bg_use.length; nTile++) {
			  if (bg_sel[nTile] &&
					  (filtered_bgnd_disp_strength[1][nTile] > 0.0) &&
					  (bg_str[nTile] > clt_parameters.fcorr_inf_strength) &&
					  ((bg_overexp == null) || (bg_overexp[nTile] < clt_parameters.lym_overexp))
					  ){
				  bg_use[nTile] = true;
			  }
		  }
		  int num_bg = tp.clt_3d_passes.get(bg_scan).setTileOpDisparity( // other minimal strength?
				  bg_use, // boolean [] selection,
				  null); // double []  disparity); // null for 0

		  // Prepare measurement of combo-scan - remove low strength and what was used for background
		  combo_disp = tp.clt_3d_passes.get(combo_scan).getDisparity(0);
		  combo_str =  tp.clt_3d_passes.get(combo_scan).getStrength();
		  combo_use = new boolean [bg_sel.length];
		  combo_overexp = tp.clt_3d_passes.get(combo_scan).getOverexposedFraction();
		  for (int  nTile = 0 ; nTile < bg_use.length; nTile++) {
			  if (!bg_use[nTile] &&
					  (combo_str[nTile] > clt_parameters.fcorr_inf_strength) &&
					  ((combo_overexp == null) || (combo_overexp[nTile] < clt_parameters.lym_overexp))
					  ){ // other minimal strength?
				  combo_use[nTile] = true;
			  }
		  }
		  num_combo = tp.clt_3d_passes.get(combo_scan).setTileOpDisparity(
				  combo_use, // boolean [] selection,
				  combo_disp); // double []  disparity);
		  if (debugLevel > -1) {
			  System.out.println("Number of background tiles = " + num_bg+", number of lazy eye tiles = " + num_combo);
		  }
		  // measure combo



		  CLTMeasure( // perform single pass according to prepared tiles operations and disparity
				  image_data, // first index - number of image in a quad
				  saturation_imp, // boolean [][] saturation_imp, // (near) saturated pixels or null
				  clt_parameters,
				  combo_scan,
				  false, // final boolean     save_textures,
				  true,  // final boolean     save_corr,
				  null,           // final double [][] mismatch,    // null or double [12][]
				  tp.threadsMax,  // maximal number of threads to launch
				  false, // updateStatus,
				  debugLevelInner - 1);
		  if (!batch_mode && clt_parameters.show_extrinsic && (debugLevel >-1)) {
			  tp.showScan(
					  tp.clt_3d_passes.get(bg_scan),   // CLTPass3d   scan,
					  "bg_scan_post"); //String title)
			  tp.showScan(
					  tp.clt_3d_passes.get(combo_scan),   // CLTPass3d   scan,
					  "combo_scan-"+combo_scan+"_post"); //String title)
		  }

		  double [][] filtered_combo_scand_isp_strength = tp.getFilteredDisparityStrength(
				  tp.clt_3d_passes, // final ArrayList <CLTPass3d> passes,// List, first, last - to search for the already tried disparity
				  combo_scan,          // final int        measured_scan_index, // will not look at higher scans
				  0,                  //	final int        start_scan_index,
				  null , // final boolean [] bg_tiles,          // get from selected in clt_3d_passes.get(0);
				  0.02,   // whatever as null above //  clt_parameters.ex_min_over,// final double     ex_min_over,       // when expanding over previously detected (by error) background, disregard far tiles

				  ImageDtt.DISPARITY_INDEX_CM,       // final int        disp_index,
				  ImageDtt.DISPARITY_STRENGTH_INDEX, // final int        str_index,
				  null,                    // final double []  tiltXY,    // null - free with limit on both absolute (2.0?) and relative (0.2) values
				  tp.getTrustedCorrelation(),//	final double     trustedCorrelation,
				  clt_parameters.fcorr_inf_strength,     //	final double     strength_floor,
				  clt_parameters.inf_str_pow,       // final double     strength_pow,
				  clt_parameters.ly_smpl_side,           // final int        smplSide, //        = 2;      // Sample size (side of a square)
				  clt_parameters.ly_smpl_num,            // final int        smplNum, //         = 3;      // Number after removing worst (should be >1)
				  clt_parameters.ly_smpl_rms,            //	final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
				  0.0,                                   // smplRelRms,         // final double     smplRelRms, //      = 0.005;  // Maximal RMS/disparity in addition to smplRms
				  clt_parameters.fds_smpl_wnd,            //	final boolean    smplWnd, //
				  clt_parameters.fds_abs_tilt,       //	final double     max_abs_tilt, //  = 2.0; // pix per tile
				  clt_parameters.fds_rel_tilt,       //	final double     max_rel_tilt, //  = 0.2; // (pix / disparity) per tile
				  clt_parameters.tileX, //   dbg_x,              //	final int        dbg_x,
				  clt_parameters.tileX, //  dbg_y,              // final int        dbg_y,
				  debugLevelInner);        //	final int        debugLevel)
		  // update selection after filtering

		  combo_disp = tp.clt_3d_passes.get(combo_scan).getDisparity(0);
		  combo_str =  tp.clt_3d_passes.get(combo_scan).getStrength();
		  combo_use = new boolean [bg_sel.length];
		  combo_overexp = tp.clt_3d_passes.get(combo_scan).getOverexposedFraction();
		  for (int  nTile = 0 ; nTile < bg_use.length; nTile++) {
			  if (!bg_use[nTile] &&
					  (combo_str[nTile] > clt_parameters.fcorr_inf_strength) &&
					  (filtered_combo_scand_isp_strength[1][nTile] > 0) &&
					  ((combo_overexp == null) || (combo_overexp[nTile] < clt_parameters.lym_overexp))
					  ){ // other minimal strength?
				  combo_use[nTile] = true;
			  }
		  }
		  int num_combo1 = tp.clt_3d_passes.get(combo_scan).setTileOpDisparity(
				  combo_use, // boolean [] selection,
				  combo_disp); // double []  disparity);
		  if (debugLevel > -1) {
			  System.out.println("Updated number of lazy eye tiles = " + num_combo1+" (was "+num_combo+")");
		  }

		  if (!batch_mode && clt_parameters.show_extrinsic && (debugLevel >-1)) {
			  String [] titles = {"bgnd_disp","bgnd_str","combo_disp","combo_str","bg_sel","bg_use","combo_use"};
			  double [] dbg_bg_sel = new double [bg_sel.length];
			  double [] dbg_bg_use =   new double [bg_sel.length];
			  double [] dbg_combo_use = new double [bg_sel.length];
			  for (int i= 0; i < bg_sel.length; i++) {
				  dbg_bg_sel[i] =    bg_sel[i]? 1.0:0.0; //only sky, no far mountains (too high disparity!)
				  dbg_bg_use[i] =    bg_use[i]? 1.0:0.0;
				  dbg_combo_use[i] = combo_use[i]? 1.0:0.0;
			  }
			  double [][]dbg_img = {
					  filtered_bgnd_disp_strength[0],
					  filtered_bgnd_disp_strength[1],
					  filtered_combo_scand_isp_strength[0],
					  filtered_combo_scand_isp_strength[1],
					  dbg_bg_sel,
					  dbg_bg_use, // too few
					  dbg_combo_use};
			  (new showDoubleFloatArrays()).showArrays(dbg_img,  tp.getTilesX(), tp.getTilesY(), true, "extrinsics_bgnd_combo",titles);
		  }
		  AlignmentCorrection ac = new AlignmentCorrection(this);
		  // iteration steps
		  double comp_diff = min_sym_update + 1; // (> min_sym_update)
		  for (int num_iter = 0; num_iter < max_tries; num_iter++){
			  if (update_disp_from_latest) {
				  tp.clt_3d_passes.get(combo_scan).updateDisparity();
			  }
			  double [][] bg_mismatch = new double[12][];
			  double [][] combo_mismatch = new double[12][];
			  CLTMeasure( // perform single pass according to prepared tiles operations and disparity
					  image_data, // first index - number of image in a quad
					  saturation_imp, // boolean [][] saturation_imp, // (near) saturated pixels or null
					  clt_parameters,
					  bg_scan,
					  false,             // final boolean     save_textures,
					  true,              // final boolean     save_corr,
					  bg_mismatch,    // final double [][] mismatch,    // null or double [12][]
					  tp.threadsMax,  // maximal number of threads to launch
					  false, // updateStatus,
					  debugLevelInner - 1);
			  CLTMeasure( // perform single pass according to prepared tiles operations and disparity
					  image_data, // first index - number of image in a quad
					  saturation_imp, // boolean [][] saturation_imp, // (near) saturated pixels or null
					  clt_parameters,
					  combo_scan,
					  false,             // final boolean     save_textures,
					  true,              // final boolean     save_corr,
					  combo_mismatch,    // final double [][] mismatch,    // null or double [12][]
					  tp.threadsMax,  // maximal number of threads to launch
					  false, // updateStatus,
					  debugLevelInner - 1);
			  double [][] scans14 = new double [28][];
			  scans14[14 * 0 + 0] =  tp.clt_3d_passes.get(bg_scan).disparity_map[ImageDtt.DISPARITY_INDEX_CM]; // .getDisparity(0);
			  scans14[14 * 0 + 1] =  tp.clt_3d_passes.get(bg_scan).getStrength();
			  scans14[14 * 1 + 0] =  tp.clt_3d_passes.get(combo_scan).disparity_map[ImageDtt.DISPARITY_INDEX_CM];
			  scans14[14 * 1 + 1] =  tp.clt_3d_passes.get(combo_scan).getStrength();
			  for (int i = 0; i < bg_mismatch.length; i++) {
				  scans14[14 * 0 + 2 + i] =  bg_mismatch[i];
				  scans14[14 * 1 + 2 + i] =  combo_mismatch[i];
			  }
			  if (debugLevelInner > 0) {
				  (new showDoubleFloatArrays()).showArrays(scans14, tp.getTilesX(), tp.getTilesY(), true, "scans_14"); //  , titles);
			  }

			  if (!batch_mode && clt_parameters.show_extrinsic && (debugLevel > 1)) {
				  tp.showScan(
						  tp.clt_3d_passes.get(bg_scan),   // CLTPass3d   scan,
						  "bg_scan_iter"); //String title)
				  tp.showScan(
						  tp.clt_3d_passes.get(combo_scan),   // CLTPass3d   scan,
						  "combo_scan-"+combo_scan+"_iter"); //String title)
			  }

			  double [][] target_disparity = {tp.clt_3d_passes.get(bg_scan).getDisparity(0), tp.clt_3d_passes.get(combo_scan).getDisparity(0)};





// TODO: fix above for using GT
// use 		lazyEyeCorrectionFromGT(..) when ground truth data is available
			  double [][][] new_corr = ac.lazyEyeCorrection(
					  adjust_poly,                       // final boolean use_poly,
					  true, // final boolean    restore_disp_inf, // Restore subtracted disparity for scan #0 (infinity)
					  clt_parameters.fcorr_radius,       // 	final double fcorr_radius,
					  clt_parameters.fcorr_inf_strength, //  final double min_strenth,
					  clt_parameters.fcorr_inf_diff,     // final double max_diff,
					  //								1.3, // final double comp_strength_var,
					  clt_parameters.inf_iters,          // 20, // 0, // final int max_iterations,
					  clt_parameters.inf_final_diff,     // 0.0001, // final double max_coeff_diff,
					  clt_parameters.inf_far_pull,       // 0.0, // 0.25, //   final double far_pull, //  = 0.2; // 1; //  0.5;
					  clt_parameters.inf_str_pow,        // 1.0, //   final double     strength_pow,
					  0.8*clt_parameters.disp_scan_step, // 1.5, // final double     lazyEyeCompDiff, // clt_parameters.fcorr_disp_diff
					  clt_parameters.ly_smpl_side,       // 3,   // final int        lazyEyeSmplSide, //        = 2;      // Sample size (side of a square)
					  clt_parameters.ly_smpl_num,        // 5,   // final int        lazyEyeSmplNum, //         = 3;      // Number after removing worst (should be >1)
					  clt_parameters.ly_smpl_rms,        // 0.1, // final double     lazyEyeSmplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
					  clt_parameters.ly_disp_var,        // 0.2, // final double     lazyEyeDispVariation, // 0.2, maximal full disparity difference between tgh tile and 8 neighborxs
					  clt_parameters.ly_disp_rvar,       // 0.2, // final double     lazyEyeDispRelVariation, // 0.02 Maximal relative full disparity difference to 8 neighbors
					  clt_parameters.ly_norm_disp,       //	final double     ly_norm_disp, //  =    5.0;     // Reduce weight of higher disparity tiles
					  clt_parameters.inf_smpl_side,      // 3, //   final int        smplSide, //        = 2;      // Sample size (side of a square)
					  clt_parameters.inf_smpl_num,       // 5, //   final int        smplNum,  //         = 3;      // Number after removing worst (should be >1)
					  clt_parameters.inf_smpl_rms,       // 0.1, // 0.05, //  final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
					  // histogram parameters
					  clt_parameters.ih_smpl_step,       // 8,    // final int        hist_smpl_side, // 8 x8 masked, 16x16 sampled
					  clt_parameters.ih_disp_min,        // -1.0, // final double     hist_disp_min,
					  clt_parameters.ih_disp_step,       // 0.05, // final double     hist_disp_step,
					  clt_parameters.ih_num_bins,        // 40,   // final int        hist_num_bins,
					  clt_parameters.ih_sigma,           // 0.1,  // final double     hist_sigma,
					  clt_parameters.ih_max_diff,        // 0.1,  // final double     hist_max_diff,
					  clt_parameters.ih_min_samples,     // 10,   // final int        hist_min_samples,
					  clt_parameters.ih_norm_center,     // true, // final boolean    hist_norm_center, // if there are more tiles that fit than min_samples, replace with
					  clt_parameters.ly_inf_frac,        // 0.5, // final double     inf_fraction,    // fraction of the weight for the infinity tiles
					  clt_parameters,  // EyesisCorrectionParameters.CLTParameters           clt_parameters,
					  scans14, // disp_strength, // scans,   // double [][] disp_strength,
					  target_disparity,          // double [][]      target_disparity, // null or programmed disparity (1 per each 14 entries of scans_14)
					  tp.getTilesX(), // int         tilesX,
					  clt_parameters.corr_magic_scale, // double      magic_coeff, // still not understood coefficent that reduces reported disparity value.  Seems to be around 8.5
					  debugLevelInner - 1); //  + (clt_parameters.fine_dbg ? 1:0)); // int debugLevel)
			  comp_diff = 0.0;
			  int num_pars = 0;
			  if (adjust_poly) {
				  apply_fine_corr(
						  new_corr,
						  debugLevelInner + 2);
				  for (int n = 0; n < new_corr.length; n++){
					  for (int d = 0; d < new_corr[n].length; d++){
						  for (int i = 0; i < new_corr[n][d].length; i++){
							  comp_diff += new_corr[n][d][i] * new_corr[n][d][i];
							  num_pars++;
						  }
					  }
				  }
				  comp_diff = Math.sqrt(comp_diff/num_pars);
				  if (debugLevel > -2) {
					  if ((debugLevel > -1) || (comp_diff < min_poly_update)) {
						  System.out.println("#### fine correction iteration step = "+(num_iter + 1) + " ( of "+max_tries+") change = "+
								  comp_diff + " ("+min_poly_update+")");
					  }
				  }

				  if (comp_diff < min_poly_update) { // add other parameter to exit from poly
					  break;
				  }
			  } else {
				  for (int i = 0; i < new_corr[0][0].length; i++){
					  comp_diff += new_corr[0][0][i] * new_corr[0][0][i];
				  }
				  comp_diff = Math.sqrt(comp_diff);
				  if (debugLevel > -10) { // should work even in batch mode
					  System.out.println("#### extrinsicsCLT(): iteration step = "+(num_iter + 1) + " ( of "+max_tries+") change = "+
							  comp_diff + " ("+min_sym_update+"), previous RMS = " + new_corr[0][1][0]+ " (debugLevel = "+debugLevel+")");
				  }

				  if (debugLevel > -2) {
					  if ((debugLevel > -1) || (comp_diff < min_sym_update)) {
//						  System.out.println("#### extrinsicsCLT(): iteration step = "+(num_iter + 1) + " ( of "+max_tries+") change = "+
//								  comp_diff + " ("+min_sym_update+"), previous RMS = " + new_corr[0][1][0]);
						  System.out.println("New extrinsic corrections:");
						  System.out.println(geometryCorrection.getCorrVector().toString());
					  }
				  }

				  if (comp_diff < min_sym_update) {
					  break;
				  }
			  }
		  }
		  return (comp_diff < (adjust_poly ? min_poly_update : min_sym_update));
	  }


	  public boolean extrinsicsCLTfromGT(
			  TwoQuadCLT       twoQuadCLT, //maybe null in no-rig mode, otherwise may contain rig measurements to be used as infinity ground truth
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  boolean 		   adjust_poly,
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel)
	  {

		  final boolean    filter_ds = clt_parameters.lyr_filter_ds;  //  = false; // true;
		  final boolean    filter_lyf = clt_parameters.lyr_filter_lyf; //  = false; // ~clt_parameters.lyf_filter, but may be different, now off for a single cameras

		  final boolean    batch_mode = clt_parameters.batch_run;
		  int debugLevelInner =  batch_mode ? -5: debugLevel;
		  int max_tries =                   clt_parameters.lym_iter; // 25;
		  double min_sym_update =           clt_parameters.lym_change; //  4e-6; // stop iterations if no angle changes more than this
		  double min_poly_update =          clt_parameters.lym_poly_change; //  Parameter vector difference to exit from polynomial correction
		  if ((twoQuadCLT == null) || (twoQuadCLT.getBiScan(0) == null)){
			  System.out.println("Rig data is not available, aborting");
			  return false;
		  }
		  BiScan scan = twoQuadCLT.getBiScan(0);
		  double [][] rig_disp_strength = 		scan.getDisparityStrength(
		    		true,   // final boolean only_strong,
		    		true,   // final boolean only_trusted,
		    		true) ; // final boolean only_enabled);

		  if (debugLevel > 20) {
			  boolean tmp_exit = true;
			  System.out.println("extrinsicsCLTfromGT()");
			  if (tmp_exit) {
				  System.out.println("will now exit. To continue - change variable tmp_exit in debugger" );
				  if (tmp_exit) {
					  return false;
				  }
			  }
		  }
		  GeometryCorrection geometryCorrection_main = null;
		  if (geometryCorrection.getRotMatrix(true) != null) {
			  geometryCorrection_main = twoQuadCLT.quadCLT_main.getGeometryCorrection();
			  double disparityScale =  geometryCorrection.getDisparityRadius()/geometryCorrection_main.getDisparityRadius();
			  for (int i = 0; i < rig_disp_strength[0].length; i++) {
				  rig_disp_strength[0][i] *= disparityScale;
			  }

			  if (debugLevel > -2) {
				  System.out.println("This is an AUX camera, using MAIN camera coordinates");
			  }


		  }

		  CLTPass3d comboScan = tp.compositeScan(
				  rig_disp_strength[0], // final double []             disparity,
				  rig_disp_strength[1], //  final double []             strength,
				  null, //  final boolean []            selected,
				  debugLevel); // final int                   debugLevel)
		  // comboScan will remain the same through iterations, no need to update disparity (maybe shrink selection?

		  AlignmentCorrection ac = new AlignmentCorrection(this);
		  // iteration steps
		  double comp_diff = min_sym_update + 1; // (> min_sym_update)
		  for (int num_iter = 0; num_iter < max_tries; num_iter++){

			  double [][] combo_mismatch = new double[12][];
			  CLTMeasure( // perform single pass according to prepared tiles operations and disparity
					  image_data,              // first index - number of image in a quad
					  saturation_imp,          // boolean [][] saturation_imp, // (near) saturated pixels or null
					  clt_parameters,
					  comboScan,               // final CLTPass3d     scan,
					  false,                   // final boolean     save_textures,
					  true,                    // final boolean     save_corr,
					  combo_mismatch,          // final double [][] mismatch,    // null or double [12][]
					  geometryCorrection_main, // final GeometryCorrection geometryCorrection_main, // If not null - covert to main camera coordinates
					  tp.threadsMax,           // maximal number of threads to launch
					  false,                   // updateStatus,
					  debugLevelInner - 1);

			  double [][] scans14 = new double [14][];
			  scans14[14 * 0 + 0] =  comboScan.disparity_map[ImageDtt.DISPARITY_INDEX_CM]; // .getDisparity(0);
			  scans14[14 * 0 + 1] =  comboScan.getStrength();
			  for (int i = 0; i < combo_mismatch.length; i++) {
				  scans14[14 * 0 + 2 + i] =  combo_mismatch[i];
			  }
			  if (debugLevelInner > 0) {
				  (new showDoubleFloatArrays()).showArrays(scans14, tp.getTilesX(), tp.getTilesY(), true, "scans_14"); //  , titles);
			  }

			  if (!batch_mode && clt_parameters.show_extrinsic && (debugLevel > 1)) {
				  tp.showScan(
						  comboScan,   // CLTPass3d   scan,
						  "combo_scan-"+num_iter+"_iter"); //String title)
			  }

			  double [][][] new_corr;
			  final double     inf_max_disparity = 2.0;
			  double [][][]    gt_disparity_strength = {rig_disp_strength};
			  new_corr = ac.lazyEyeCorrectionFromGT(
					  geometryCorrection_main, //final GeometryCorrection geometryCorrection_main, // if not null - this is an AUX camera of a rig
					  adjust_poly,                       // final boolean use_poly,
					  true, // final boolean    restore_disp_inf, // Restore subtracted disparity for scan #0 (infinity)
					  clt_parameters.fcorr_radius,       // 	final double fcorr_radius,
					  clt_parameters.fcorr_inf_strength, //  final double min_strenth,

					  clt_parameters.inf_str_pow,        // 1.0, //   final double     strength_pow,
					  0.8*clt_parameters.disp_scan_step, // 1.5, // final double     lazyEyeCompDiff, // clt_parameters.fcorr_disp_diff
					  clt_parameters.ly_smpl_side,       // 3,   // final int        lazyEyeSmplSide, //        = 2;      // Sample size (side of a square)
					  clt_parameters.ly_smpl_num,        // 5,   // final int        lazyEyeSmplNum, //         = 3;      // Number after removing worst (should be >1)
					  clt_parameters.ly_smpl_rms,        // 0.1, // final double     lazyEyeSmplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
					  clt_parameters.ly_disp_var,        // 0.2, // final double     lazyEyeDispVariation, // 0.2, maximal full disparity difference between tgh tile and 8 neighborxs
					  clt_parameters.ly_disp_rvar,       // 0.2, // final double     lazyEyeDispRelVariation, // 0.02 Maximal relative full disparity difference to 8 neighbors
					  clt_parameters.ly_norm_disp,       //	final double     ly_norm_disp, //  =    5.0;     // Reduce weight of higher disparity tiles
					  clt_parameters.inf_smpl_side,      // 3, //   final int        smplSide, //        = 2;      // Sample size (side of a square)
					  clt_parameters.inf_smpl_num,       // 5, //   final int        smplNum,  //         = 3;      // Number after removing worst (should be >1)
					  clt_parameters.inf_smpl_rms,       // 0.1, // 0.05, //  final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
					  // histogram parameters
					  clt_parameters.ih_smpl_step,       // 8,    // final int        hist_smpl_side, // 8 x8 masked, 16x16 sampled
					  clt_parameters.ih_disp_min,        // -1.0, // final double     hist_disp_min,
					  clt_parameters.ih_disp_step,       // 0.05, // final double     hist_disp_step,
					  clt_parameters.ih_num_bins,        // 40,   // final int        hist_num_bins,
					  clt_parameters.ih_sigma,           // 0.1,  // final double     hist_sigma,
					  clt_parameters.ih_max_diff,        // 0.1,  // final double     hist_max_diff,
					  clt_parameters.ih_min_samples,     // 10,   // final int        hist_min_samples,
					  clt_parameters.ih_norm_center,     // true, // final boolean    hist_norm_center, // if there are more tiles that fit than min_samples, replace with
					  clt_parameters.ly_inf_frac,        // 0.5, // final double     inf_fraction,    // fraction of the weight for the infinity tiles
					  inf_max_disparity, // final double     inf_max_disparity, // use all smaller disparities as inf_fraction
					  clt_parameters,  // EyesisCorrectionParameters.CLTParameters           clt_parameters,
					  scans14, // disp_strength, // scans,   // double [][] disp_strength,
					  gt_disparity_strength, // double [][][]    gt_disparity_strength, // 1 pair for each 14 entries of scans_14 (normally - just 1 scan
					  filter_ds,   // final boolean    filter_ds, //
					  filter_lyf,   // final boolean    filter_lyf, // ~clt_parameters.lyf_filter, but may be different, now off for a single cameras
					  tp.getTilesX(), // int         tilesX,
					  clt_parameters.corr_magic_scale, // double      magic_coeff, // still not understood coefficent that reduces reported disparity value.  Seems to be around 8.5
					  debugLevelInner - 1); //  + (clt_parameters.fine_dbg ? 1:0)); // int debugLevel)
			  comp_diff = 0.0;
			  int num_pars = 0;
			  if (adjust_poly) {
				  apply_fine_corr(
						  new_corr,
						  debugLevelInner + 2);
				  for (int n = 0; n < new_corr.length; n++){
					  for (int d = 0; d < new_corr[n].length; d++){
						  for (int i = 0; i < new_corr[n][d].length; i++){
							  comp_diff += new_corr[n][d][i] * new_corr[n][d][i];
							  num_pars++;
						  }
					  }
				  }
				  comp_diff = Math.sqrt(comp_diff/num_pars);
				  if (debugLevel > -2) {
					  if ((debugLevel > -1) || (comp_diff < min_poly_update)) {
						  System.out.println("#### fine correction iteration step = "+(num_iter + 1) + " ( of "+max_tries+") change = "+
								  comp_diff + " ("+min_poly_update+")");
					  }
				  }

				  if (comp_diff < min_poly_update) { // add other parameter to exit from poly
					  break;
				  }
			  } else {
				  for (int i = 0; i < new_corr[0][0].length; i++){
					  comp_diff += new_corr[0][0][i] * new_corr[0][0][i];
				  }
				  comp_diff = Math.sqrt(comp_diff);
				  if (debugLevel > -10) { // should work even in batch mode
					  System.out.println("#### extrinsicsCLTfromGT(): iteration step = "+(num_iter + 1) + " ( of "+max_tries+") change = "+
							  comp_diff + " ("+min_sym_update+"), previous RMS = " + new_corr[0][1][0]+ " (debugLevel = "+debugLevel+")");
				  }

				  if (debugLevel > -2) {
					  if ((debugLevel > -1) || (comp_diff < min_sym_update)) {
//						  System.out.println("#### extrinsicsCLT(): iteration step = "+(num_iter + 1) + " ( of "+max_tries+") change = "+
//								  comp_diff + " ("+min_sym_update+"), previous RMS = " + new_corr[0][1][0]);
						  System.out.println("New extrinsic corrections:");
						  System.out.println(geometryCorrection.getCorrVector().toString());
					  }
				  }
				  if (comp_diff < min_sym_update) {
					  break;
				  }
			  }
		  }
		  return (comp_diff < (adjust_poly ? min_poly_update : min_sym_update));
	  }






	  public boolean expandCLTQuad3d(
		  EyesisCorrectionParameters.CLTParameters           clt_parameters,
		  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
		  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
		  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
		  EyesisCorrectionParameters.RGBParameters             rgbParameters,
		  final int        threadsMax,  // maximal number of threads to launch
		  final boolean    updateStatus,
		  final int        debugLevel)
	  {
		  final boolean                       batch_mode = clt_parameters.batch_run; //disable any debug images
		  final int debugLevelInner =  batch_mode ? -3: debugLevel;
		  final double     trustedCorrelation = tp.getTrustedCorrelation();
		  final int        max_expand =  500; // 150; // 30;
  		  final boolean    show_retry_far =   clt_parameters.show_retry_far && false; // (max_expand <= 10);
		  boolean          show_expand =      false; //    clt_parameters.show_expand && (max_expand <= 10);

		  final int        disp_index =      ImageDtt.DISPARITY_INDEX_CM;
		  final int        str_index =       ImageDtt.DISPARITY_STRENGTH_INDEX;
		  final double     strength_floor =  clt_parameters.fds_str_floor;  // 0.6* clt_parameters.combine_min_strength;

		  // TODO: make parameters
		  final double     strength_pow    = clt_parameters.fds_str_pow ;   // 1.0;
		  final int        smplSide        = clt_parameters.fds_smpl_side ;       // 5; // 3;      // Sample size (side of a square)
		  final int        smplNum         = clt_parameters.fds_smpl_num ;        // 10; // 13; // 5;      // Number after removing worst (should be >1)
		  final double     smplRms         = clt_parameters.fds_smpl_rms ;        // 0.15; // Maximal RMS of the remaining tiles in a sample
		  final double     smplRelRms      = clt_parameters.fds_smpl_rel_rms ;     // 0.01; // 05;  // Maximal RMS/disparity in addition to smplRms
		  final boolean    smplWnd         = clt_parameters.fds_smpl_wnd ;        // true; //
		  final double     max_abs_tilt    = clt_parameters.fds_abs_tilt ;   // 2.0; // pix per tile
		  final double     max_rel_tilt    = clt_parameters.fds_rel_tilt ;   // 0.2; // (pix / disparity) per tile
		  final int        bg_pass = 0;

		  final int        dbg_x = 155;
		  final int        dbg_y = 207;

          int num_extended = -1;
// process once more to try combining of processed
          boolean last_pass = false;
//          for (int num_expand = 0; (num_expand < 4) && (num_extended != 0); num_expand++) {
          boolean over_infinity = false;
          int dbg_start_pass = 4; // 10; // 20;
          int dbg_end_pass =   -20; // 12; // -29;
//          final int        dbg_x0 = 73; //
//          final int        dbg_y0 = 195; //
  		  double [][] filtered_bgnd_disp_strength = tp.getFilteredDisparityStrength(
				  tp.clt_3d_passes, // final ArrayList <CLTPass3d> passes,// List, first, last - to search for the already tried disparity
				  bg_pass, // tp.clt_3d_passes.size() - 1,         // final int        measured_scan_index, // will not look at higher scans
  				  0,                  //	final int        start_scan_index,
  				  null , // final boolean [] bg_tiles,          // get from selected in clt_3d_passes.get(0);
				  0.0,   // whatever as null above //  clt_parameters.ex_min_over,// final double     ex_min_over,       // when expanding over previously detected (by error) background, disregard far tiles

  				  disp_index,         // final int        disp_index,
  				  str_index,          // final int        str_index,
  				  null,               // final double []  tiltXY,    // null - free with limit on both absolute (2.0?) and relative (0.2) values
  				  trustedCorrelation, //	final double     trustedCorrelation,
  				  strength_floor,     //	final double     strength_floor,
  				  strength_pow,       // final double     strength_pow,
  				  smplSide,           // final int        smplSide, //        = 2;      // Sample size (side of a square)
  				  smplNum,            // final int        smplNum, //         = 3;      // Number after removing worst (should be >1)
  				  smplRms,            //	final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
  				  smplRelRms,         // final double     smplRelRms, //      = 0.005;  // Maximal RMS/disparity in addition to smplRms

  				  smplWnd,            //	final boolean    smplWnd, //
  				  max_abs_tilt,       //	final double     max_abs_tilt, //  = 2.0; // pix per tile
  				  max_rel_tilt,       //	final double     max_rel_tilt, //  = 0.2; // (pix / disparity) per tile
  				  dbg_x,              //	final int        dbg_x,
  				  dbg_y,              // final int        dbg_y,
  				  debugLevelInner);   //	final int        debugLevel)



          for (int num_expand = 0; num_expand < max_expand; num_expand++) {
        	  boolean dbg_pass = (num_expand >= dbg_start_pass) && (num_expand <= dbg_end_pass);
        	  show_expand =      (clt_parameters.show_expand && (max_expand <= 10)) || dbg_pass;

        	  //        	  Runtime runtime = Runtime.getRuntime();
        	  //       	  runtime.gc();
        	  //       	  System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");

        	  num_extended = zMapExpansionStep(
        			  tp.clt_3d_passes,               // final ArrayList <CLTPass3d> passes,// List, first, last - to search for the already tried disparity
        			  clt_parameters,                 //final EyesisCorrectionParameters.CLTParameters           clt_parameters, // for refinePassSetup()
        			  0,                              // final int         firstPass,
        			  tp.clt_3d_passes.size(),        // final int         lastPassPlus1,
        			  bg_pass,                        // final int         bg_index,

        			  true,                           // final boolean     refine,         // now always should be true ?
        			  clt_parameters.gr_new_expand,   // final boolean     expand_neibs,   // expand from neighbors
        			  !clt_parameters.gr_new_expand,  // final boolean     expand_legacy,  // old mode of growing tiles (using max scanned). If both true, will try only if expand_neibs fails

        			  // Filtering by background data:
        			  filtered_bgnd_disp_strength,    // final double [][] filtered_bgnd_ds,       // if not null, will filter results not to have low disparity new tiles over supposed bgnd
        			  // for legacy expansion it is not used, but just checked for null, so double [0][] should work too
        			  clt_parameters.ex_min_over ,              // final double      ex_min_over,            // when expanding over previously detected (by error) background, disregard far tiles
        			  clt_parameters.gr_ovrbg_filtered ,  // final double      str_over_bg,            // Minimal filtered strength when replacing background data
        			  clt_parameters.gr_ovrbg_cmb ,     // final double      str_over_bg_combo,      // Minimal combined strength  when replacing background data
        			  clt_parameters.gr_ovrbg_cmb_hor , // final double      str_over_bg_combo_hor,  // Minimal combined strength  when replacing background data
        			  clt_parameters.gr_ovrbg_cmb_vert, // final double      str_over_bg_combo_vert, // Minimal combined strength  when replacing background data

        			  // Refine parameters use directly some clt_parameters, others below
        			  clt_parameters.stUseRefine,               // final boolean     stUseRefine, // use supertiles (now false)
        			  clt_parameters.bgnd_range,                //final double      bgnd_range,     // double            disparity_far,
        			  clt_parameters.ex_strength,               // final double      ex_strength,    // double            this_sure,        // minimal strength to be considered definitely good
        			  clt_parameters.ex_nstrength,              //final double      ex_nstrength,   // double            ex_nstrength, // minimal 4-corr strength divided by channel diff for new (border) tiles
        			  clt_parameters.bgnd_maybe,                // final double      bgnd_maybe,     // double            this_maybe,       // maximal strength to ignore as non-background
        			  clt_parameters.sure_smth,                 // final double      sure_smth,      // sure_smth,        // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
        			  clt_parameters.pt_super_trust,            //final double      pt_super_trust, //  final double      super_trust,      // If strength exceeds ex_strength * super_trust, do not apply ex_nstrength and plate_ds
        			  // using plates disparity/strength - averaged for small square sets of tiles. If null - just use raw tiles
        			  clt_parameters.pt_keep_raw_fg,            // final boolean     pt_keep_raw_fg, // 		final boolean     keep_raw_fg,  // do not replace raw tiles by the plates, if raw is closer (like poles)
        			  clt_parameters.pt_scale_pre,              // final double      pt_scale_pre,   // final double      scale_filtered_strength_pre, // scale plate_ds[1] before comparing to raw strength
        			  clt_parameters.pt_scale_post,             // final double      pt_scale_post,  // final double      scale_filtered_strength_post,// scale plate_ds[1] when replacing raw (generally plate_ds is more reliable if it exists)

        			  // Composite scan parameters
        			  clt_parameters.combine_min_strength,      // final double      combine_min_strength,              // final double                minStrength,
        			  clt_parameters.combine_min_hor,           // final double      combine_min_hor,                   // final double                minStrengthHor,
        			  clt_parameters.combine_min_vert,          // final double      combine_min_vert,                  // final double                minStrengthVert,

        			  //getFilteredDisparityStrength parameters

        			  // Consider separate parameters for FDS?
        			  clt_parameters.fds_str_floor,        // final double      fds_str_floor,
        			  clt_parameters.fds_str_pow,          // final double      fds_str_pow,
        			  clt_parameters.fds_smpl_side,             // final int         fds_smpl_size, // == 5
        			  clt_parameters.fds_smpl_num,              // final int         fds_smpl_points, // == 3
        			  clt_parameters.fds_smpl_rms,              // final double      fds_abs_rms,
        			  clt_parameters.fds_smpl_rel_rms,          // final double      fds_rel_rms,
        			  clt_parameters.fds_smpl_wnd,              // final boolean     fds_use_wnd,   // use window function fro the neighbors
        			  0.001,                                    // final double      fds_tilt_damp,
        			  clt_parameters.fds_abs_tilt,              // final double      fds_abs_tilt, //   = 2.0; // pix per tile
        			  clt_parameters.fds_rel_tilt,              // final double      fds_rel_tilt, //   = 0.2; // (pix / disparity) per tile

        			  // expansion from neighbors parameters:
        			  clt_parameters.gr_min_new ,               // final int         gr_min_new,         // discard variant if there are less new tiles
        			  clt_parameters.gr_steps_over,             // final int         gr_steps_over,  // how far to extend
        			  clt_parameters.gr_var_new_sngl ,          // final boolean     gr_var_new_sngl,
        			  clt_parameters.gr_var_new_fg ,            // final boolean     gr_var_new_fg,
        			  clt_parameters.gr_var_all_fg ,            // final boolean     gr_var_all_fg,
        			  clt_parameters.gr_var_new_bg ,            // final boolean     gr_var_new_bg,
        			  clt_parameters.gr_var_all_bg ,            // final boolean     gr_var_all_bg,
        			  clt_parameters.gr_var_next ,              // final boolean     gr_var_next,

        			  clt_parameters.gr_smpl_size , // final int         gr_smpl_size ,      // 5,      // 	final int         smpl_size, // == 5
        			  clt_parameters.gr_min_pnts , // final int         gr_min_pnts ,       // 5, // 4, // 3,      //    final int        min_points, // == 3
        			  clt_parameters.gr_use_wnd , // final boolean     gr_use_wnd ,        // true,   // 	final boolean     use_wnd,   // use window function fro the neighbors
        			  clt_parameters.gr_tilt_damp , // final double      gr_tilt_damp ,      // 0.001,  // 	final double      tilt_cost,
        			  clt_parameters.gr_split_rng , // final double      gr_split_rng ,      // 5.0,    // 	final double      split_threshold, // if full range of the values around the cell higher, need separate fg, bg
        			  clt_parameters.gr_same_rng , // final double      gr_same_rng ,       // 3.0,    // 	final double      same_range,      // modify
        			  clt_parameters.gr_diff_cont , // final double      gr_diff_cont ,      // 2.0,    // 	final double      diff_continue,   // maximal difference from the old value (for previously defined tiles
        			  clt_parameters.gr_abs_tilt , // final double      gr_abs_tilt ,       // 2.0,    //    final double     max_abs_tilt, //   = 2.0; // pix per tile
        			  clt_parameters.gr_rel_tilt , // final double      gr_rel_tilt ,       // 0.2,    //    final double     max_rel_tilt, //   = 0.2; // (pix / disparity) per tile
        			  clt_parameters.gr_smooth , // final int         gr_smooth ,         // 50,     // 	final int         max_tries, // maximal number of smoothing steps
        			  clt_parameters.gr_fin_diff , // final double      gr_fin_diff ,       // 0.01,   // 	final double      final_diff, // maximal change to finish iterations
        			  clt_parameters.gr_unique_pretol , // final double      gr_unique_pretol ,  // 1.0,    // 	final double      unique_pre_tolerance, // usually larger than clt_parameters.unique_tolerance

        			  // legacy expansion
        			  clt_parameters.grow_min_diff ,  // final double      grow_min_diff,   // =    0.5; // Grow more only if at least one channel has higher variance from others for the tile
        			  clt_parameters.grow_retry_far , // final boolean     grow_retry_far,  // final boolean     grow_retry_far,  // Retry tiles around known foreground that have low max_tried_disparity
        			  clt_parameters.grow_pedantic ,  // final boolean     grow_pedantic,   // final boolean     grow_pedantic,   // Scan full range between max_tried_disparity of the background and known foreground
        			  clt_parameters.grow_retry_inf , // final boolean     grow_retry_inf,  // final boolean     grow_retry_inf,  // Retry border tiles that were identified as infinity earlier

        			  // common expansion
        			  clt_parameters.gr_num_steps,    // final int         gr_num_steps,    // how far to extend
        			  clt_parameters.gr_unique_tol ,  // final double      gr_unique_tol,
        			  0.0, // clt_parameters.grow_disp_min ,  // final double      grow_disp_min,   // final double                disp_near,
        			  clt_parameters.grow_disp_max ,  // final double      grow_disp_max,   // final double                disp_near,
        			  clt_parameters.grow_disp_step , // final double      grow_disp_step,  //  =   6.0; // Increase disparity (from maximal tried) if nothing found in that tile // TODO: handle enclosed dips?

        			  // for debug level > 1, just corresponding clt_parameters.* will work, otherwise following 3 can override
        			  last_pass,                      // final boolean     last_pass, // just for more debug?
        			  true,                           // final boolean     remove_non_measurement, // save memory, true is OK
        			  show_expand,                    // final boolean     show_expand,
        			  clt_parameters. show_unique,    // final boolean     show_unique,
        			  show_retry_far,                 // final boolean     show_retry_far,

        			  dbg_x,
        			  dbg_y,             // final int        dbg_y,
        			  debugLevelInner + 2);        //	final int        debugLevel)

        	  if (last_pass) {
        		  break;
        	  } else if (num_extended == 0){
        		  System.out.println("**** processCLTQuad3d(): nothing to expand ***");
        		  System.out.println("!clt_parameters.ex_over_bgnd="+clt_parameters.ex_over_bgnd+" over_infinity="+over_infinity);
        		  if (!clt_parameters.ex_over_bgnd || over_infinity)  last_pass = true;
        		  else {
        			  over_infinity = true;
        			  if (debugLevel > -1){
        				  System.out.println("===== processCLTQuad3d(): trying to expand over previously identified background (may be by error)====");
        			  }
        		  }
        	  }

          }






          show_expand =      (clt_parameters.show_expand && (max_expand <= 10));
          if (debugLevelInner > -1) {
        	  tp.showScan(
        			  tp.clt_3d_passes.get(tp.clt_3d_passes.size()-2),   // CLTPass3d   scan,
        			  "after_pre_last_combo_pass-"+(tp.clt_3d_passes.size()-2)); //String title)
        	  tp.showScan(
        			  tp.clt_3d_passes.get(tp.clt_3d_passes.size()-1),   // CLTPass3d   scan,
        			  "after_last_combo_pass-"+(tp.clt_3d_passes.size()-1)); //String title)
          }

///          int refine_pass = tp.clt_3d_passes.size(); // refinePassSetup() will add one - not anymore!

///          int next_pass = tp.clt_3d_passes.size(); //
          tp.secondPassSetup( // prepare tile tasks for the second pass based on the previous one(s)
        		  clt_parameters,
        		  clt_parameters.stUsePass2, // use supertiles
        		  bg_pass,
				  // disparity range - differences from
				  clt_parameters.bgnd_range, // double            disparity_far,
//    				  -0.5, // 0.0, // clt_parameters.bgnd_range,                //	 final double                disp_far,   // limit results to the disparity range
				  clt_parameters.grow_disp_max, // other_range, //double            disparity_near,   //
				  clt_parameters.ex_strength,  // double            this_sure,        // minimal strength to be considered definitely good
				  clt_parameters.ex_nstrength, // double            ex_nstrength, // minimal 4-corr strength divided by channel diff for new (border) tiles
				  clt_parameters.bgnd_maybe, // double            this_maybe,       // maximal strength to ignore as non-background
				  clt_parameters.sure_smth, // sure_smth,        // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
				  clt_parameters.pt_super_trust, // super_trust,      // If strength exceeds ex_strength * super_trust, do not apply ex_nstrength and plate_ds

				  ImageDtt.DISPARITY_INDEX_CM,  // index of disparity value in disparity_map == 2 (0,2 or 4)
				  geometryCorrection,
				  threadsMax,  // maximal number of threads to launch
				  updateStatus,
				  debugLevelInner);

// Save tp.clt_3d_passes.size() to roll back without restarting the program
          tp.saveCLTPasses(false); // not rig, and reset rig data
          Runtime runtime = Runtime.getRuntime();
          runtime.gc();
          System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");
          return true; //  null;
	  }

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// called after composite scan is added to the list? composite scan is inside

	  public int zMapExpansionStep(
				final ArrayList <CLTPass3d> passes,// List, first, last - to search for the already tried disparity
				final EyesisCorrectionParameters.CLTParameters           clt_parameters, // for refinePassSetup()
				final int         firstPass,
				final int         lastPassPlus1,
				final int         bg_index,


				final boolean     refine,         // now always should be true ?
				final boolean     expand_neibs,   // expand from neighbors
				final boolean     expand_legacy,  // old mode of growing tiles (using max scanned). If both true, will try only if expand_neibs fails

				// Filtering by background data:
				final double [][] filtered_bgnd_ds,       // if not null, will filter results not to have low disparity new tiles over supposed bgnd
				                                          // for legacy expansion it is not used, but just checked for null, so double [0][] should work too
				final double      ex_min_over,            // when expanding over previously detected (by error) background, disregard far tiles
				final double      str_over_bg,            // Minimal filtered strength when replacing background data
				final double      str_over_bg_combo,      // Minimal combined strength  when replacing background data
				final double      str_over_bg_combo_hor,  // Minimal combined strength  when replacing background data
				final double      str_over_bg_combo_vert, // Minimal combined strength  when replacing background data

				// Refine parameters use directly some clt_parameters, others below
				final boolean     stUseRefine, // use supertiles
				final double      bgnd_range,     // double            disparity_far,
				final double      ex_strength,    // double            this_sure,        // minimal strength to be considered definitely good
				final double      ex_nstrength,   // double            ex_nstrength, // minimal 4-corr strength divided by channel diff for new (border) tiles
				final double      bgnd_maybe,     // double            this_maybe,       // maximal strength to ignore as non-background
				final double      sure_smth,      // sure_smth,        // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
				final double      pt_super_trust, //  final double      super_trust,      // If strength exceeds ex_strength * super_trust, do not apply ex_nstrength and plate_ds
				  // using plates disparity/strength - averaged for small square sets of tiles. If null - just use raw tiles
				final boolean     pt_keep_raw_fg, // 		final boolean     keep_raw_fg,  // do not replace raw tiles by the plates, if raw is closer (like poles)
				final double      pt_scale_pre,   // final double      scale_filtered_strength_pre, // scale plate_ds[1] before comparing to raw strength
				final double      pt_scale_post,  // final double      scale_filtered_strength_post,// scale plate_ds[1] when replacing raw (generally plate_ds is more reliable if it exists)

				// Composite scan parameters
				final double      combine_min_strength,              // final double                minStrength,
				final double      combine_min_hor,                   // final double                minStrengthHor,
				final double      combine_min_vert,                  // final double                minStrengthVert,

				//getFilteredDisparityStrength parameter
				final double      fds_str_floor,
				final double      fds_str_pow,
				final int         fds_smpl_size, // == 5
				final int         fds_smpl_points, // == 3
				final double      fds_abs_rms,
				final double      fds_rel_rms,
				final boolean     fds_use_wnd,   // use window function fro the neighbors
				final double      fds_tilt_damp,
				final double      fds_abs_tilt, //   = 2.0; // pix per tile
				final double      fds_rel_tilt, //   = 0.2; // (pix / disparity) per tile

				// expansion from neighbors parameters:
				final int         gr_min_new,         // discard variant if there are less new tiles
				final int         gr_steps_over,  // how far to extend
				final boolean     gr_var_new_sngl,
				final boolean     gr_var_new_fg,
				final boolean     gr_var_all_fg,
				final boolean     gr_var_new_bg,
				final boolean     gr_var_all_bg,
				final boolean     gr_var_next,

				final int         gr_smpl_size ,      // 5,      // 	final int         smpl_size, // == 5
				final int         gr_min_pnts ,       // 5, // 4, // 3,      //    final int        min_points, // == 3
				final boolean     gr_use_wnd ,        // true,   // 	final boolean     use_wnd,   // use window function fro the neighbors
				final double      gr_tilt_damp ,      // 0.001,  // 	final double      tilt_cost,
				final double      gr_split_rng ,      // 5.0,    // 	final double      split_threshold, // if full range of the values around the cell higher, need separate fg, bg
				final double      gr_same_rng ,       // 3.0,    // 	final double      same_range,      // modify
				final double      gr_diff_cont ,      // 2.0,    // 	final double      diff_continue,   // maximal difference from the old value (for previously defined tiles
				final double      gr_abs_tilt ,       // 2.0,    //    final double     max_abs_tilt, //   = 2.0; // pix per tile
				final double      gr_rel_tilt ,       // 0.2,    //    final double     max_rel_tilt, //   = 0.2; // (pix / disparity) per tile
				final int         gr_smooth ,         // 50,     // 	final int         max_tries, // maximal number of smoothing steps
				final double      gr_fin_diff ,       // 0.01,   // 	final double      final_diff, // maximal change to finish iterations
				final double      gr_unique_pretol ,  // 1.0,    // 	final double      unique_pre_tolerance, // usually larger than clt_parameters.unique_tolerance
				// legacy expansion

				final double      grow_min_diff,   // =    0.5; // Grow more only if at least one channel has higher variance from others for the tile
				final boolean     grow_retry_far,  // final boolean     grow_retry_far,  // Retry tiles around known foreground that have low max_tried_disparity
				final boolean     grow_pedantic,   // final boolean     grow_pedantic,   // Scan full range between max_tried_disparity of the background and known foreground
				final boolean     grow_retry_inf,  // final boolean     grow_retry_inf,  // Retry border tiles that were identified as infinity earlier

				// common expansion
				final int         gr_num_steps,    // how far to extend
				final double      gr_unique_tol,
				final double      grow_disp_min,   // final double                disp_near,
				final double      grow_disp_max,   // final double                disp_near,
				final double      grow_disp_step,  //  =   6.0; // Increase disparity (from maximal tried) if nothing found in that tile // TODO: handle enclosed dips?

				// for debug level > 1, just corresponding clt_parameters.* will work, otherwise following 3 can override
				final boolean     last_pass, // just for more debug?
				final boolean     remove_non_measurement, // save memory, true is OK
				final boolean     show_expand,
				final boolean     show_unique,
				final boolean     show_retry_far,

				final int         dbg_x,
				final int         dbg_y,
				final int         debugLevel)
	  {
    	  boolean dbg_pass = debugLevel > 1;
		  final CLTPass3d   bg_scan = passes.get(bg_index);         // background scan data, null - ignore background

    	  Runtime runtime = Runtime.getRuntime();
    	  runtime.gc();
    	  System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");


    	  // Get filtered (by flexible "plates" that can tilt to accommodate tiles disparity/strength map, that uses data from all previous
    	  // disparity "measurements". This data will be combined with individual tiles, quad, hor and vert correlation results

    	  double [][] filtered_disp_strength = tp.getFilteredDisparityStrength(
				  tp.clt_3d_passes, // final ArrayList <CLTPass3d> passes,// List, first, last - to search for the already tried disparity
    			  lastPassPlus1 - 1,                   // final int        measured_scan_index, // will not look at higher scans
    			  0,                                  // final int        start_scan_index,
    			  bg_scan.getSelected(),              // selected , // final boolean [] bg_tiles,          // get from selected in clt_3d_passes.get(0);
    			  ex_min_over,                        // final double     ex_min_over,       // when expanding over previously detected (by error) background, disregard far tiles
    			  ImageDtt.DISPARITY_INDEX_CM,        // final int        disp_index,
    			  ImageDtt.DISPARITY_STRENGTH_INDEX,  // final int        str_index,
    			  null,                               // final double []  tiltXY,    // null - free with limit on both absolute (2.0?) and relative (0.2) values
    			  tp.getTrustedCorrelation(),         // final double     trustedCorrelation,
    			  fds_str_floor,                 //	final double     strength_floor,
    			  fds_str_pow,                   // final double     strength_pow,
    			  fds_smpl_size,                      // final int        smplSide, //        = 2;      // Sample size (side of a square)
    			  fds_smpl_points,                    // final int        smplNum, //         = 3;      // Number after removing worst (should be >1)
    			  fds_abs_rms,                        //	final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
    			  fds_rel_rms,                        // final double     smplRelRms, //      = 0.005;  // Maximal RMS/disparity in addition to smplRms
    			  fds_use_wnd,                        //	final boolean    smplWnd, //
    			  fds_abs_tilt,                       //	final double     max_abs_tilt, //  = 2.0; // pix per tile
    			  fds_rel_tilt,                       //	final double     max_rel_tilt, //  = 0.2; // (pix / disparity) per tile
    			  dbg_x,                              //	final int        dbg_x,
    			  dbg_y,                              // final int        dbg_y,
    			  debugLevel+1);                      //	final int        debugLevel)

    	  // Optionally filter by background data to avoid far objects in the areas already identified as bgnd
    	  if (filtered_bgnd_ds != null) {
    		  tp.filterOverBackground(
    				  filtered_disp_strength,         // final double [][]           ds,
    				  filtered_bgnd_ds[1],            // final double []             bg_strength,
    				  bg_scan.getSelected(),          // final boolean []            bg_tiles,          // get from selected in clt_3d_passes.get(0);
    				  str_over_bg,                    // final double                minStrength,
    				  ex_min_over);                   // final double                ex_min_over        // when expanding over previously detected (by error) background, disregard far tiles
    	  }
    	  int refine_pass = passes.size() - 1; // last index plus 1

    	  // refine pass uses hor/vert thresholds inside
    	  // Now always start with refine - it will set disparity to improve previous results where residual disparity was significant
    	  if (refine) { // currently is broken if !refine
    		  CLTPass3d refined = tp.refinePassSetup( // prepare tile tasks for the refine pass (re-measure disparities), add it to the list
    				  clt_parameters,
    				  stUseRefine, // use supertiles
    				  bg_index, // bg_scan bg_pass,
    				  // disparity range - differences from
    				  bgnd_range,     // double            disparity_far,
    				  grow_disp_max,  // other_range, //double            disparity_near,   //
    				  ex_strength,    // double            this_sure,        // minimal strength to be considered definitely good
    				  ex_nstrength,   // double            ex_nstrength, // minimal 4-corr strength divided by channel diff for new (border) tiles
    				  bgnd_maybe,     // double            this_maybe,       // maximal strength to ignore as non-background
    				  sure_smth,      // sure_smth,        // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
    				  pt_super_trust, //  final double      super_trust,      // If strength exceeds ex_strength * super_trust, do not apply ex_nstrength and plate_ds
    				  // using plates disparity/strength - averaged for small square sets of tiles. If null - just use raw tiles
    				  filtered_disp_strength,        // 		final double [][] plate_ds,  // disparity/strength last time measured for the multi-tile squares. Strength =-1 - not measured. May be null
    				  pt_keep_raw_fg, // 		final boolean     keep_raw_fg,  // do not replace raw tiles by the plates, if raw is closer (like poles)
    				  pt_scale_pre,   // final double      scale_filtered_strength_pre, // scale plate_ds[1] before comparing to raw strength
    				  pt_scale_post,  // final double      scale_filtered_strength_post,// scale plate_ds[1] when replacing raw (generally plate_ds is more reliable if it exists)
    				  ImageDtt.DISPARITY_INDEX_CM,   // index of disparity value in disparity_map == 2 (0,2 or 4)
    				  geometryCorrection,
    				  tp.threadsMax,  // maximal number of threads to launch
    				  false, // updateStatus,
    				  debugLevel);
    		  passes.add(refined); // adding new scan, not yet measured
    		  refine_pass ++; // tp.refinePassSetup adds to the list
    		  if (show_expand || (clt_parameters.show_expand && dbg_pass)) {
    			  tp.showScan(
    					  passes.get(refine_pass), // CLTPass3d   scan,
    					  "after_refine-"+refine_pass);
    		  }
    	  }
    	  // Maybe it will not be used later, still calculate
    	  tp.calcMaxTried(
    			  passes, // final ArrayList <CLTPass3d> passes,
    			  bg_index, //  final int                   firstPass,
    			  refine_pass, // may add 1 to include current (for future?)     //  final int                   lastPassPlus1,
    			  passes.get(refine_pass)); //  final int                   lastPassPlus1,
    	  // repeated composite scan may be replaced by just a clone

    	  CLTPass3d extended_pass = tp.compositeScan(
    			  passes,             // final ArrayList <CLTPass3d> passes,
    			  bg_index,           //  final int                   firstPass,
    			  passes.size(),      //  final int                   lastPassPlus1,
    			  //      				  tp.clt_3d_passes.get(bg_pass).getSelected(), // selected , // final boolean [] bg_tiles,          // get from selected in clt_3d_passes.get(0);
    			  //    				  clt_parameters.ex_min_over,// final double     ex_min_over,       // when expanding over previously detected (by error) background, disregard far tiles
    			  tp.getTrustedCorrelation(),                       // 	 final double                trustedCorrelation,
    			  tp.getMaxOverexposure(),                          //  final double                max_overexposure,
    			  0.0, // clt_parameters.bgnd_range,                //	 final double                disp_far,   // limit results to the disparity range
    			  grow_disp_max,                     // final double                disp_near,
    			  combine_min_strength,              // final double                minStrength,
    			  combine_min_hor,                   // final double                minStrengthHor,
    			  combine_min_vert,                  // final double                minStrengthVert,
				  false,      // final boolean               no_weak,
    			  true, // false, // final boolean               use_last,   //
    			  // TODO: when useCombo - pay attention to borders (disregard)
    			  false, // final boolean               usePoly)  // use polynomial method to find max), valid if useCombo == false
    			  true, // 	 final boolean               copyDebug)
    			  debugLevel - 1);
    	  if (show_expand || (clt_parameters.show_expand && last_pass)) tp.showScan(
    			  passes.get(refine_pass), // CLTPass3d   scan,
    			  "after_refine-combine0-"+(passes.size() - 1));
    	  // TODO: Verify that IMG_DIFF0_INDEX has combined data
    	  if (filtered_bgnd_ds != null) {
    		  tp.filterOverBackground(
    				  extended_pass,           // final CLTPass3d            pass,
    				  str_over_bg_combo,       //	 final double             minStrength,
    				  str_over_bg_combo_hor,   //	 final double             minStrengthHor,
    				  str_over_bg_combo_vert,  //final double                minStrengthVert,
    				  bg_scan.getSelected(),   // selected , // final boolean [] bg_tiles,          // get from selected in clt_3d_passes.get(0);
    				  ex_min_over,             // final double     ex_min_over,       // when expanding over previously detected (by error) background, disregard far tiles
    				  filtered_disp_strength);
    	  }
    	  if (show_expand || (clt_parameters.show_expand && last_pass)) tp.showScan(
    			  passes.get(refine_pass), // CLTPass3d   scan,
    			  "after_refine-combine1-"+(passes.size() - 1));


    	  if (show_expand || (clt_parameters.show_expand && last_pass)) {
    		  final int tilesX = tp.getTilesX();
    		  final int tilesY = tp.getTilesY();
    		  String [] titles = {"disp","strength", "disp_combined","str_combined","max_tried","selected"};
    		  String title = "FDS_"+(passes.size() - 1);
    		  double [][] dbg_img = new double [titles.length][];
    		  dbg_img[0] = filtered_disp_strength[0];
    		  dbg_img[1] = filtered_disp_strength[1];
    		  dbg_img[2] = extended_pass.getDisparity();
    		  dbg_img[3] = extended_pass.getStrength();
    		  double [][] max_tried_disparity = extended_pass.getMaxTriedDisparity();
    		  if (max_tried_disparity != null){
    			  dbg_img[4] = new double [tilesX * tilesY];
    			  for (int ty = 0; ty < tilesY; ty++){
    				  for (int tx = 0; tx < tilesX; tx++){
    					  dbg_img[4][ty*tilesX + tx] = max_tried_disparity[ty][tx];
    				  }
    			  }
    		  }
    		  boolean [] s_selected =  extended_pass.selected;
    		  boolean [] s_border =    extended_pass.border_tiles;
    		  if ((s_selected != null) && (s_border != null)){
    			  dbg_img[5] = new double [tilesX * tilesY];
    			  for (int i = 0; i < dbg_img[5].length; i++){
    				  dbg_img[5][i] = 1.0 * ((s_selected[i]?1:0) + (s_border[i]?2:0));
    			  }
    		  }

    		  (new showDoubleFloatArrays()).showArrays(dbg_img,  tilesX, tilesY, true, title,titles);
    	  }
    	  int [] numLeftRemoved = {0, 0};
    	  int num_extended = 0;
    	  if (expand_neibs) { // show_expand){ // use new mode
    		  if (last_pass){
    			  System.out.println("+++++++++++ Last pass ++++++++++++");
    		  }

    		  //        			public boolean []
    		  tp.dbg_filtered_disp_strength = filtered_disp_strength; // to be shown inside
    		  boolean [] variants_flags = {
    				  gr_var_new_sngl,
    				  gr_var_new_fg,
    				  gr_var_all_fg,
    				  gr_var_new_bg,
    				  gr_var_all_bg,
    				  gr_var_next};
    		  numLeftRemoved = tp.prepareExpandVariant(
    				  extended_pass,                     // final CLTPass3d   scan,            // combined scan with max_tried_disparity, will be modified to re-scan
    				  passes.get(refine_pass), // final CLTPass3d   last_scan,       // last prepared tile - can use last_scan.disparity, .border_tiles and .selected
    				  /*null, // */ passes.get(bg_index),     // final CLTPass3d   bg_scan,         // background scan data
    				  passes,                  // final ArrayList <CLTPass3d> passes,// List, first, last - to search for the already tried disparity
    				  0,                                 // 	final int         firstPass,
    				  passes.size(),           //         					final int         lastPassPlus1,
    				  gr_min_new ,        // 20, // 	final int         min_new,         // discard variant if there are less new tiles
    				  variants_flags, // 0x1e, // 0x1f,                     // final int         variants_mask,
    				  gr_num_steps ,      // 8, // 	final int         num_steps,  // how far to extend
    				  gr_steps_over ,     // 8, // 4, // final int         num_steps_disc,  // how far to extend
    				  gr_smpl_size ,      // 5,      // 	final int         smpl_size, // == 5
    				  gr_min_pnts ,       // 5, // 4, // 3,      //    final int        min_points, // == 3
    				  gr_use_wnd ,        // true,   // 	final boolean     use_wnd,   // use window function fro the neighbors
    				  gr_tilt_damp ,      // 0.001,  // 	final double      tilt_cost,
    				  gr_split_rng ,      // 5.0,    // 	final double      split_threshold, // if full range of the values around the cell higher, need separate fg, bg
    				  gr_same_rng ,       // 3.0,    // 	final double      same_range,      // modify
    				  gr_diff_cont ,      // 2.0,    // 	final double      diff_continue,   // maximal difference from the old value (for previously defined tiles
    				  gr_abs_tilt ,       // 2.0,    //    final double     max_abs_tilt, //   = 2.0; // pix per tile
    				  gr_rel_tilt ,       // 0.2,    //    final double     max_rel_tilt, //   = 0.2; // (pix / disparity) per tile
    				  gr_smooth ,         // 50,     // 	final int         max_tries, // maximal number of smoothing steps
    				  gr_fin_diff ,       // 0.01,   // 	final double      final_diff, // maximal change to finish iterations
    				  grow_disp_min,      // final double                grow_disp_min,
    				  grow_disp_max,      // final double                grow_disp_max,
    				  grow_disp_step,      // final double                grow_disp_step,
    				  gr_unique_pretol ,  // 1.0,    // 	final double      unique_pre_tolerance, // usually larger than clt_parameters.unique_tolerance
//    				  last_pass?  gr_unique_tol : gr_unique_pretol, // 	final double      unique_tolerance,
    				  gr_unique_tol, // last_pass?  clt_parameters.gr_unique_tol : clt_parameters.gr_unique_pretol, // 	final double      unique_tolerance,
    				  show_expand,       // 	final boolean     show_expand,
    				  show_unique,       // 	final boolean     show_unique,
    				  dbg_x,
    				  dbg_y,
    				  debugLevel);
        	  num_extended = numLeftRemoved[0];
			  if (clt_parameters.show_expand || (clt_parameters.show_variant && (numLeftRemoved[1] > 1 ))) tp.showScan(
					  tp.clt_3d_passes.get(refine_pass), // CLTPass3d   scan,
					  "prepareExpandVariant-"+numLeftRemoved[1]+"-"+refine_pass); //String title)
    	  }
    	  if ((num_extended == 0) && expand_legacy) { // if both are on, will use legacy if neighbors failed

    		  boolean show_ex_debug = show_retry_far || (clt_parameters.show_retry_far && last_pass) || dbg_pass;
    		  num_extended = tp.setupExtendDisparity(
    				  extended_pass,                              // final CLTPass3d   scan,            // combined scan with max_tried_disparity, will be modified to re-scan
    				  passes.get(refine_pass), // final CLTPass3d   last_scan,       // last prepared tile - can use last_scan.disparity, .border_tiles and .selected
    				  (filtered_bgnd_ds != null)? null: passes.get(bg_index), // final CLTPass3d   bg_scan,         // background scan data
    						  gr_num_steps, // clt_parameters.grow_sweep,      // 8; // Try these number of tiles around known ones
    						  grow_disp_max,   //  =   50.0; // Maximal disparity to try
    						  0.5 * tp.getTrustedCorrelation(), // clt_parameters.grow_disp_trust, //  =  4.0; // Trust measured disparity within +/- this value
    						  grow_disp_step,  //  =   6.0; // Increase disparity (from maximal tried) if nothing found in that tile // TODO: handle enclosed dips?
    						  grow_min_diff,   // =    0.5; // Grow more only if at least one channel has higher variance from others for the tile
    						  grow_retry_far,  // final boolean     grow_retry_far,  // Retry tiles around known foreground that have low max_tried_disparity
    						  grow_pedantic,   // final boolean     grow_pedantic,   // Scan full range between max_tried_disparity of the background and known foreground
    						  grow_retry_inf,  // final boolean     grow_retry_inf,  // Retry border tiles that were identified as infinity earlier
    						  clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
    						  geometryCorrection, // GeometryCorrection geometryCorrection,
    						  show_ex_debug,     // true, // final boolean     show_debug,

    						  tp.threadsMax,  // maximal number of threads to launch
    						  false, // updateStatus,
    						  debugLevel);
    		  //TODO:  break if nothing wanted? - no, there are some left to be refined

    		  if (debugLevel > -1){
    			  System.out.println("=== setupExtendDisparity() added "+num_extended+" new tiles to scan");
    		  }

    		  numLeftRemoved = tp.makeUnique(
    				  passes,                  // final ArrayList <CLTPass3d> passes,
    				  0,                                 //  final int                   firstPass,
    				  refine_pass, // - 1,                   //  final int                   lastPassPlus1,
    				  extended_pass, // tp.clt_3d_passes.get(refine_pass), //  final CLTPass3d             new_scan,
    				  clt_parameters.grow_disp_max,       // final double                grow_disp_max,
    				  clt_parameters.gr_unique_tol,   //  final double                unique_tolerance,
    				  clt_parameters.show_unique);      // final boolean               show_unique)

        	  num_extended = numLeftRemoved[0];

    		  if (debugLevel > -1){
    			  System.out.println("last makeUnique("+refine_pass+") -> left: "+numLeftRemoved[0]+", removed:" + numLeftRemoved[1]);
    		  }



    		  //        		  num_extended = numLeftRemoved[0];
    		  //TODO:  break if nothing wanted? - here yes, will make sense



    	  }

		  refine_pass = passes.size(); // 1

		  passes.add(extended_pass);

		  if (show_expand) tp.showScan(
				  passes.get(refine_pass), // CLTPass3d   scan,
				  "before_measure-"+refine_pass); //String title)

//    	  num_extended = numLeftRemoved[0];

    	  CLTMeasure( // perform single pass according to prepared tiles operations and disparity
    			  image_data, // first index - number of image in a quad
				  saturation_imp, //final boolean [][]  saturation_imp, // (near) saturated pixels or null
    			  clt_parameters,
    			  refine_pass,
    			  false, // final boolean     save_textures,
    			  tp.threadsMax,  // maximal number of threads to launch
    			  false, // updateStatus,
    			  debugLevel - 1);
    	  if (remove_non_measurement) {
    		  if (debugLevel > 0){
    			  System.out.print("Removing non-measurement(derivative) scan data, current list size = "+passes.size());
    		  }
    		  tp.removeNonMeasurement();
    		  refine_pass = passes.size() - 1;
    		  if (debugLevel > 0){
    			  System.out.println(", new size = "+passes.size());
    		  }
    	  }


    	  if (show_expand || (clt_parameters.show_expand && last_pass)) {
    		  tp.showScan(
    				  passes.get(refine_pass), // CLTPass3d   scan,
    				  "after_measure-"+refine_pass); //String title)

    	  }
    	  if (debugLevel > -1){
    		  System.out.println("extending: CLTMeasure("+refine_pass+"),  num_extended = "+num_extended);
    	  }
    	  //num_expand == 9
    	  CLTPass3d combo_pass = tp.compositeScan(
    			  passes, // final ArrayList <CLTPass3d> passes,
    			  bg_index,                                         //  final int                   firstPass,
    			  passes.size(),                          //  final int                   lastPassPlus1,
    			  tp.getTrustedCorrelation(),                       // 	 final double                trustedCorrelation,
    			  tp.getMaxOverexposure(),                          //  final double                max_overexposure,
    			  -0.5, // 0.0, // clt_parameters.bgnd_range,                //	 final double                disp_far,   // limit results to the disparity range
    			  clt_parameters.grow_disp_max,                     // final double                disp_near,
    			  clt_parameters.combine_min_strength,              // final double                minStrength,
    			  clt_parameters.combine_min_hor,                   // final double                minStrengthHor,
    			  clt_parameters.combine_min_vert,                  // final double                minStrengthVert,
				  false,      // final boolean               no_weak,
    			  false, // final boolean               use_last,   //
    			  // TODO: when useCombo - pay attention to borders (disregard)
    			  false, // final boolean               usePoly)  // use polynomial method to find max), valid if useCombo == false
    			  true, // 	 final boolean               copyDebug)
    			  //    				  (num_expand == 9)? 2: debugLevel);
    			  debugLevel -1);
    	  passes.add(combo_pass);
    	  if (show_expand || (clt_parameters.show_expand && last_pass)) tp.showScan(
    			  passes.get(passes.size() -1), // refine_pass), // CLTPass3d   scan,
    			  "after_combo_pass-"+(passes.size() -1)); // (refine_pass)); //String title)
    	  return num_extended; // [0];
	  }

//*****************************************************************

//	  public ImagePlus output3d(
	  public boolean output3d(
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel)
	  {
		  final boolean    batch_mode = clt_parameters.batch_run;
		  this.startStepTime=System.nanoTime();
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  if (this.image_data == null){
			  return false;
		  }
		  double infinity_disparity = 	geometryCorrection.getDisparityFromZ(clt_parameters.infinityDistance);
		  X3dOutput x3dOutput = null;
		  WavefrontExport wfOutput = null;
		  if (clt_parameters.remove_scans){
			  System.out.println("Removing all scans but the first(background) and the last to save memory");
			  System.out.println("Will need to re-start the program to be able to output differently");
			  CLTPass3d   latest_scan = tp.clt_3d_passes.get(tp.clt_3d_passes.size() - 1);
			  tp.trimCLTPasses(1);
			  tp.clt_3d_passes.add(latest_scan); // put it back
		  }
		  int next_pass = tp.clt_3d_passes.size(); //
		  tp.thirdPassSetupSurf( // prepare tile tasks for the second pass based on the previous one(s) // needs last scan
				  clt_parameters,
				  //FIXME: make a special parameter?
				  infinity_disparity, //0.25 * clt_parameters.bgnd_range, // double            disparity_far,
				  clt_parameters.grow_disp_max, // other_range, //double            disparity_near,   //
				  geometryCorrection,
				  threadsMax,  // maximal number of threads to launch
				  updateStatus,
				  0); // final int         debugLevel)
		  if (!batch_mode && (debugLevel > -1)) {
			  tp.showScan(
					  tp.clt_3d_passes.get(next_pass-1),   // CLTPass3d   scan,
					  "after_pass3-"+(next_pass-1)); //String title)
		  }
 		  String x3d_path= correctionsParameters.selectX3dDirectory( // for x3d and obj
 				 correctionsParameters.getModelName(this.image_name), // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
 				  correctionsParameters.x3dModelVersion,
				  true,  // smart,
				  true);  //newAllowed, // save

 			// create x3d file
		  if (clt_parameters.output_x3d) {
			  x3dOutput = new X3dOutput(
					clt_parameters,
					correctionsParameters,
					geometryCorrection,
					tp.clt_3d_passes);
		  }
		  if (clt_parameters.output_obj && (x3d_path != null)) {
			  try {
				wfOutput = new WavefrontExport(
						  x3d_path,
						  correctionsParameters.getModelName(this.image_name),
						  clt_parameters,
						  correctionsParameters,
						  geometryCorrection,
						  tp.clt_3d_passes);
			} catch (IOException e) {
				System.out.println("Failed to open Wavefront files for writing");
				// TODO Auto-generated catch block
				e.printStackTrace();
				// do nothing, just keep
			}
		  }
		  if (x3dOutput != null) {
			  x3dOutput.generateBackground(clt_parameters.infinityDistance <= 0.0); // needs just first (background) scan
		  }

		  int bgndIndex = 0; // it already exists?
		  CLTPass3d bgndScan = tp.clt_3d_passes.get(bgndIndex);
//		  boolean [] bgnd_sel = bgndScan.getSelected().clone();
//		  int num_bgnd = 0;
//		  for (int i = 0; i < bgnd_sel.length; i++) if (bgnd_sel[i]) num_bgnd++;
//		  if (num_bgnd >= clt_parameters.min_bgnd_tiles) { // TODO: same for the backdrop too
//		  double infinity_disparity = 	geometryCorrection.getDisparityFromZ(clt_parameters.infinityDistance);
		  if (bgndScan.texture != null) { // TODO: same for the backdrop too
			  if (clt_parameters.infinityDistance > 0.0){ // generate background as a billboard
				  // grow selection, then grow once more and create border_tiles
				  // create/rstore, probably not needed
				  boolean [] bg_sel_backup = bgndScan.getSelected().clone();
				  boolean [] bg_border_backup = (bgndScan.getBorderTiles() == null) ? null: bgndScan.getBorderTiles().clone();
				  boolean [] bg_selected = bgndScan.getSelected();
				  // 			  tp.growTiles(
				  // 					  2,                   // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				  // 					 bg_selected,
				  // 					  null); // prohibit
				  boolean [] border_tiles = bg_selected.clone();
				  tp.growTiles(
						  2,                   // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
						  bg_selected,
						  null); // prohibit
				  // 			  for (int)

				  for (int i = 0; i < border_tiles.length; i++){
					  border_tiles[i] = !border_tiles[i] && bg_selected[i];
				  }
				  // update texture_tiles (remove what is known not to be background
				  for (int ty = 0; ty < tilesY; ty++){
					  for (int tx = 0; tx < tilesX; tx++){
						  if (!bg_selected[tx + tilesX*ty]){
							  bgndScan.texture_tiles[ty][tx] = null;
						  }
					  }
				  }

				  bgndScan.setBorderTiles(border_tiles);
				  // updates selection from non-null texture tiles
				  String texturePath = getPassImage( // get image from a single pass
						  clt_parameters,
						  colorProcParameters,
						  rgbParameters,
						  correctionsParameters.getModelName(this.image_name)+"-img_infinity", // +scanIndex,
						  bgndIndex,
						  threadsMax,  // maximal number of threads to launch
						  updateStatus,
						  batch_mode ? -5: debugLevel);
				  double [] scan_disparity = new double [tilesX * tilesY];
				  int indx = 0;
				  //		  boolean [] scan_selected = scan.getSelected();
				  for (int ty = 0; ty < tilesY; ty ++) for (int tx = 0; tx < tilesX; tx ++){
					  scan_disparity[indx++] = infinity_disparity;
				  }
				  //			  tp.showScan(
				  //    				  scan, // CLTPass3d   scan,
				  //    				  "infinityDistance");

				  boolean showTri = false; // ((scanIndex < next_pass + 1) && clt_parameters.show_triangles) ||((scanIndex - next_pass) == 73);
				  try {
					  generateClusterX3d(
							  x3dOutput,
							  wfOutput,  // output WSavefront if not null
							  texturePath,
							  "INFINITY", // id (scanIndex - next_pass), // id
							  "INFINITY", // class
							  bgndScan.getTextureBounds(),
							  bgndScan.selected,
							  scan_disparity, // scan.disparity_map[ImageDtt.DISPARITY_INDEX_CM],
							  clt_parameters.transform_size,
							  clt_parameters.correct_distortions, // requires backdrop image to be corrected also
							  showTri, // (scanIndex < next_pass + 1) && clt_parameters.show_triangles,
							  infinity_disparity,  // 0.3
							  clt_parameters.grow_disp_max, // other_range, // 2.0 'other_range - difference from the specified (*_CM)
							  clt_parameters.maxDispTriangle);
				  } catch (IOException e) {
					  // TODO Auto-generated catch block
					  e.printStackTrace();
					  return false;
				  }
				  // maybe not needed
				  bgndScan.setBorderTiles(bg_border_backup);
				  bgndScan.setSelected(bg_sel_backup);
			  }
		  }






		  for (int scanIndex = next_pass; scanIndex < tp.clt_3d_passes.size(); scanIndex++){
			  if (debugLevel > 0){
				  System.out.println("FPGA processing scan #"+scanIndex);
			  }
///			  CLTPass3d scan =
					  CLTMeasure( // perform single pass according to prepared tiles operations and disparity
					  image_data, // first index - number of image in a quad
//					  saturation_imp, //final boolean [][]  saturation_imp, // (near) saturated pixels or null
					  clt_parameters,
					  scanIndex,
    				  true,  // final boolean     save_textures,
    				  false, // final boolean     save_corr,
					  threadsMax,  // maximal number of threads to launch
					  updateStatus,
					  batch_mode ? -5: debugLevel);

		  }

		  // TEMPORARY EXIT

		  //      if (tp.clt_3d_passes.size() > 0) return null; // just to fool compiler


		  //	  int scan_limit = 10;
		  for (int scanIndex = next_pass; (scanIndex < tp.clt_3d_passes.size()) && (scanIndex < clt_parameters.max_clusters); scanIndex++){ // just temporary limiting
			  if (debugLevel > -1){
				  System.out.println("Generating cluster images (limit is set to "+clt_parameters.max_clusters+") largest, scan #"+scanIndex);
			  }
			  //		  ImagePlus cluster_image = getPassImage( // get image form a single pass
			  String texturePath = getPassImage( // get image from a single pass
					  clt_parameters,
					  colorProcParameters,
					  rgbParameters,
					  correctionsParameters.getModelName(this.image_name)+"-img"+scanIndex,
					  scanIndex,
					  threadsMax,  // maximal number of threads to launch
					  updateStatus,
					  batch_mode ? -5: debugLevel);

			  CLTPass3d scan = tp.clt_3d_passes.get(scanIndex);

/*
			  if ((scanIndex == 73) ) {
				  tp.showScan(
						  tp.clt_3d_passes.get(scanIndex), // CLTPass3d   scan,
						  "SELECTED-"+scanIndex);
			  }

*/

			  // TODO: use new updated disparity, for now just what was forced for the picture
			  double [] scan_disparity = new double [tilesX * tilesY];
			  int indx = 0;
			  //		  boolean [] scan_selected = scan.getSelected();
			  for (int ty = 0; ty < tilesY; ty ++) for (int tx = 0; tx < tilesX; tx ++){
				  //			  scan_selected[indx] = scan.tile_op[ty][tx] != 0;
				  scan_disparity[indx++] = scan.disparity[ty][tx];
			  }
			  if (clt_parameters.avg_cluster_disp){
				  double sw = 0.0, sdw = 0.0;
				  for (int i = 0; i< scan_disparity.length; i++){
					  if (scan.selected[i] && !scan.border_tiles[i]){
						  double w = scan.disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX][i];
						  sw +=w;
						  sdw += scan_disparity[i]*w;
					  }
				  }
				  sdw/=sw;
				  for (int i = 0; i< scan_disparity.length; i++){
					  scan_disparity[i] = sdw;
				  }
			  }
/*
			  if ((scanIndex == 73)) {
				  tp.showScan(
						  tp.clt_3d_passes.get(scanIndex), // CLTPass3d   scan,
						  "X3D-"+scanIndex);
			  }
*/
//			  boolean showTri = ((scanIndex < next_pass + 1) && clt_parameters.show_triangles) ||(scanIndex < 3);
			  boolean showTri = !batch_mode && (debugLevel > -1) && (((scanIndex < next_pass + 1) && clt_parameters.show_triangles) ||((scanIndex - next_pass) == 73));

//			  boolean showTri = ((scanIndex < next_pass + 1) && clt_parameters.show_triangles) ||(scanIndex == 49) || (scanIndex == 54);
			  try {
				generateClusterX3d(
						  x3dOutput,
						  wfOutput,  // output WSavefront if not null
						  texturePath,
						  "shape_id-"+(scanIndex - next_pass), // id
						  null, // class
						  scan.getTextureBounds(),
						  scan.selected,
						  scan_disparity, // scan.disparity_map[ImageDtt.DISPARITY_INDEX_CM],
						  clt_parameters.transform_size,
						  clt_parameters.correct_distortions, // requires backdrop image to be corrected also
						  showTri, // (scanIndex < next_pass + 1) && clt_parameters.show_triangles,
						  // FIXME: make a separate parameter:
						  infinity_disparity, //  0.25 * clt_parameters.bgnd_range,  // 0.3
						  clt_parameters.grow_disp_max, // other_range, // 2.0 'other_range - difference from the specified (*_CM)
						  clt_parameters.maxDispTriangle);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				return false;
			}
		  }

		  // now generate and save texture files (start with full, later use bounding rectangle?)


		  if ((x3d_path != null) && (x3dOutput != null)){
//			  x3d_path+=Prefs.getFileSeparator()+correctionsParameters.getModelName(this.image_name)+".x3d";
//			  x3dOutput.generateX3D(x3d_path);
			  x3dOutput.generateX3D(x3d_path+Prefs.getFileSeparator()+correctionsParameters.getModelName(this.image_name)+".x3d");
		  }
		  if (wfOutput != null){
			  wfOutput.close();
			  System.out.println("Wavefront object file saved to "+wfOutput.obj_path);
			  System.out.println("Wavefront material file saved to "+wfOutput.mtl_path);
		  }

		  // Save KML and ratings files if they do not exist (so not to overwrite edited ones), make them world-writable
		  writeKml        (debugLevel);
		  writeRatingFile (debugLevel);


		  Runtime.getRuntime().gc();
		  System.out.println("output3d(): generating 3d output files  finished at "+
				  IJ.d2s(0.000000001*(System.nanoTime()-this.startStepTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
		  return true;
	  }




	  public void generateClusterX3d(
			  X3dOutput       x3dOutput, // output x3d if not null
			  WavefrontExport wfOutput,  // output WSavefront if not null
			  String          texturePath,
			  String          id,
			  String          class_name,
			  Rectangle       bounds,
			  boolean []      selected,
			  double []       disparity, // if null, will use min_disparity
			  int             tile_size,
			  boolean         correctDistortions, // requires backdrop image to be corrected also
			  boolean         show_triangles,
			  double          min_disparity,
			  double          max_disparity,
			  double          maxDispTriangle
			  ) throws IOException
	  {
		  int [][] indices =  tp.getCoordIndices( // starting with 0, -1 - not selected
				  bounds,
				  selected);
		  double [][] texCoord = tp.getTexCoords( // get texture coordinates for indices
				  indices);
		  double [][] worldXYZ = tp.getCoords( // get world XYZ in meters for indices
				  disparity,
				  min_disparity,
				  max_disparity,
				  bounds,
				  indices,
				  tile_size,
				  correctDistortions, // requires backdrop image to be corrected also
				  this.geometryCorrection);

          double [] indexedDisparity = tp.getIndexedDisparities( // get disparity for each index
							disparity,
							min_disparity,
							max_disparity,
							bounds,
							indices,
							tile_size);

		  int [][] triangles = 	tp.triangulateCluster(
				  indices);


		  triangles = 	tp.filterTriangles(
					triangles,
					indexedDisparity, // disparities per vertex index
					maxDispTriangle); // maximal disparity difference in a triangle


		  if (show_triangles) {
			  double [] ddisp = (disparity == null)?(new double[1]):disparity;
			  if (disparity == null) {
				  ddisp[0] = min_disparity;
			  }
			  tp.testTriangles(
					  texturePath,
					  bounds,
					  selected,
					  ddisp, // disparity, // if disparity.length == 1 - use for all
					  tile_size,
					  indices,
					  triangles);
		  }
		  if (x3dOutput != null) {
		  x3dOutput.addCluster(
				  texturePath,
				  id,
				  class_name,
				  texCoord,
				  worldXYZ,
				  triangles);
		  }
		  if (wfOutput != null) {
			  wfOutput.addCluster(
				  texturePath,
				  id,
//				  class_name,
				  texCoord,
				  worldXYZ,
				  triangles);
		  }
	  }


	  public ImagePlus getBackgroundImage(
			  boolean    no_image_save,
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  String     name,
			  int        disparity_index, // index of disparity value in disparity_map == 2 (0,2 or 4)
//			  int        strength_index, // index of strength data in disparity map ==6
			  int        threadsMax,  // maximal number of threads to launch
			  boolean    updateStatus,
			  int        debugLevel
			  )
	  {
		  final boolean new_mode = false;



		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  showDoubleFloatArrays sdfa_instance = null;

		  if (clt_parameters.debug_filters && (debugLevel > -1))
//		  if ((debugLevel > -1))
			  sdfa_instance = new showDoubleFloatArrays(); // just for debugging?

		  CLTPass3d bgnd_data = tp.clt_3d_passes.get(0);
		  double [][][][] texture_tiles = bgnd_data.texture_tiles;

		  boolean [] bgnd_tiles =   tp.getBackgroundMask( // which tiles do belong to the background
				  clt_parameters.bgnd_range,     // disparity range to be considered background
				  clt_parameters.bgnd_sure,      // minimal strength to be considered definitely background
				  clt_parameters.bgnd_maybe,     // maximal strength to ignore as non-background
				  clt_parameters.sure_smth,      // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
				  clt_parameters.min_clstr_seed, // number of tiles in a cluster to seed (just background?)
				  clt_parameters.min_clstr_block,// number of tiles in a cluster to block (just non-background?)
				  disparity_index, // index of disparity value in disparity_map == 2 (0,2 or 4)
				  clt_parameters.show_bgnd_nonbgnd,
				  (clt_parameters.debug_filters ? (debugLevel) : -1));
		  boolean [] bgnd_tiles_new =   tp.getBackgroundMask_new( // which tiles do belong to the background
				  clt_parameters.bgnd_range,     // disparity range to be considered background
				  clt_parameters.bgnd_sure,      // minimal strength to be considered definitely background
				  clt_parameters.bgnd_maybe,     // maximal strength to ignore as non-background
				  clt_parameters.sure_smth,      // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
				  clt_parameters.min_clstr_seed, // number of tiles in a cluster to seed (just background?)
				  clt_parameters.min_clstr_block,// number of tiles in a cluster to block (just non-background?)
				  disparity_index, // index of disparity value in disparity_map == 2 (0,2 or 4)
				  clt_parameters.show_bgnd_nonbgnd,
				  (clt_parameters.debug_filters ? (debugLevel) : -1));
		  boolean [] bgnd_dbg =    bgnd_tiles.clone(); // only these have non 0 alpha
// TODO: fix mess - after modifying getBackgroundMask() to getBackgroundMask_new (before road with tractor was identified as a background because of a double
// tile glitch, the background was wrong. So temporarily both old/new are used and combined (new is grown twice)
// still does not work - using old variant for now

// background selections (slightly) influences the plane mertging / connections

		  for (int i = 0; i < bgnd_tiles.length; i++){
//			  bgnd_tiles[i] &= bgnd_tiles_new[i];
		  }
		  boolean [] bgnd_strict = bgnd_tiles.clone(); // only these have non 0 alpha
		  tp.growTiles(
				  clt_parameters.bgnd_grow,      // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				  bgnd_tiles,
				  null); // prohibit
		  boolean [] bgnd_tiles_grown = bgnd_tiles.clone(); // only these have non 0 alpha
		  tp.growTiles(
				  2,      // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				  bgnd_tiles,
				  null); // prohibit

// hacking - grow bgnd texture even more, without changing selection
		  boolean [] bgnd_tiles_grown2 = bgnd_tiles_grown.clone(); // only these have non 0 alpha

		  bgnd_data.selected = bgnd_tiles_grown; // selected for background w/o extra transparent layer
		  tp.growTiles(
				  2,      // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				  bgnd_tiles_grown2,
				  null); // prohibit
		  if (sdfa_instance!=null){
			  double [][] dbg_img = new double[5][tilesY * tilesX];
			  String [] titles = {"old","new","strict","grown","more_grown"};
			  for (int i = 0; i<dbg_img[0].length;i++){
				  //
				  dbg_img[0][i] =  bgnd_dbg[i]?1:0;
				  dbg_img[1][i] =  bgnd_tiles_new[i]?1:0;
				  dbg_img[2][i] =  bgnd_strict[i]?1:0;
				  dbg_img[3][i] =  bgnd_tiles_grown[i]?1:0;
				  dbg_img[4][i] =  bgnd_tiles[i]?1:0;
			  }
			sdfa_instance.showArrays(dbg_img,  tilesX, tilesY, true, "strict_grown",titles);
		  }

		  int num_bgnd = 0;
		  double [][][][] texture_tiles_bgnd = new double[tilesY][tilesX][][];
		  double [] alpha_zero = new double [4*clt_parameters.transform_size*clt_parameters.transform_size];
		  int alpha_index = 3;
		  for (int i = 0; i < alpha_zero.length; i++) alpha_zero[i]=0.0;
		  if (new_mode) {
			  for (int tileY = 0; tileY < tilesY; tileY++){
				  for (int tileX = 0; tileX < tilesX; tileX++){
					  texture_tiles_bgnd[tileY][tileX]= null;
					  if ((texture_tiles[tileY][tileX] != null) &&
							  bgnd_tiles_grown2[tileY * tilesX + tileX]) {
						  if (bgnd_tiles[tileY * tilesX + tileX]) {
							  texture_tiles_bgnd[tileY][tileX]= texture_tiles[tileY][tileX];
							  num_bgnd++;
						  }else{
							  texture_tiles_bgnd[tileY][tileX]= texture_tiles[tileY][tileX].clone();
							  texture_tiles_bgnd[tileY][tileX][alpha_index] = alpha_zero;
						  }
					  }
				  }
			  }
		  } else {
			  for (int tileY = 0; tileY < tilesY; tileY++){
				  for (int tileX = 0; tileX < tilesX; tileX++){
					  texture_tiles_bgnd[tileY][tileX]= null;
					  if ((texture_tiles[tileY][tileX] != null) &&
							  bgnd_tiles[tileY * tilesX + tileX]) {
						  if (bgnd_tiles_grown2[tileY * tilesX + tileX]) {
							  texture_tiles_bgnd[tileY][tileX]= texture_tiles[tileY][tileX];
							  num_bgnd++;
						  }else{
							  texture_tiles_bgnd[tileY][tileX]= texture_tiles[tileY][tileX].clone();
							  texture_tiles_bgnd[tileY][tileX][alpha_index] = alpha_zero;
						  }
					  }
				  }
			  }
		  }

		  if (num_bgnd < clt_parameters.min_bgnd_tiles){
			  return null; // no background to generate
		  }
		  ImageDtt image_dtt = new ImageDtt();
		  double [][] texture_overlap = image_dtt.combineRGBATiles(
				  texture_tiles_bgnd, // texture_tiles,               // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}
				  clt_parameters.transform_size,
				  true,                         // when false - output each tile as 16x16, true - overlap to make 8x8
				  clt_parameters.sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
				  threadsMax,                    // maximal number of threads to launch
				  debugLevel);
		  if (clt_parameters.alpha1 > 0){ // negative or 0 - keep alpha as it was
			  double scale = (clt_parameters.alpha1 > clt_parameters.alpha0) ? (1.0/(clt_parameters.alpha1 - clt_parameters.alpha0)) : 0.0;
			  for (int i = 0; i < texture_overlap[alpha_index].length; i++){
				  double d = texture_overlap[alpha_index][i];
				  if      (d >=clt_parameters.alpha1) d = 1.0;
				  else if (d <=clt_parameters.alpha0) d = 0.0;
				  else d = scale * (d- clt_parameters.alpha0);
				  texture_overlap[alpha_index][i] = d;
			  }
		  }
		  // for now - use just RGB. Later add oprion for RGBA
		  double [][] texture_rgb = {texture_overlap[0],texture_overlap[1],texture_overlap[2]};
		  double [][] texture_rgba = {texture_overlap[0],texture_overlap[1],texture_overlap[2],texture_overlap[3]};

		  //			  ImagePlus img_texture =
		  ImagePlus imp_texture_bgnd = linearStackToColor(
				  clt_parameters,
				  colorProcParameters,
				  rgbParameters,
				  name+"-texture-bgnd", // String name,
				  "", //String suffix, // such as disparity=...
				  true, // toRGB,
				  !this.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
				  true, // boolean saveShowIntermediate, // save/show if set globally
				  false, //true, // boolean saveShowFinal,        // save/show result (color image?)
				  ((clt_parameters.alpha1 > 0)? texture_rgba: texture_rgb),
				  tilesX *  clt_parameters.transform_size,
				  tilesY *  clt_parameters.transform_size,
				  1.0,         // double scaleExposure, // is it needed?
				  debugLevel);
		  // resize for backdrop here!
//	public double getFOVPix(){ // get ratio of 1 pixel X/Y to Z (distance to object)
		  ImagePlus imp_texture_bgnd_ext = resizeForBackdrop(
				  imp_texture_bgnd,
				  clt_parameters.black_back, //  boolean fillBlack,
				  clt_parameters.black_back, //  boolean noalpha,
				  debugLevel);
		  String path= correctionsParameters.selectX3dDirectory(
				  //TODO: Which one to use - name or this.image_name ?
 				  correctionsParameters.getModelName(this.image_name), // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
 				  correctionsParameters.x3dModelVersion,
// 				  name, // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
				  true,  // smart,
				  true);  //newAllowed, // save
		  // only show/save original size if debug or debug_filters)
		  int jpegQuality = clt_parameters.black_back ? 90: -1;
		  if (clt_parameters.debug_filters || (debugLevel > 0)) {
			  eyesisCorrections.saveAndShow(
					  imp_texture_bgnd,
					  path,
					  correctionsParameters.png && !clt_parameters.black_back,
					  clt_parameters.show_textures,
					  jpegQuality); // jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
		  }
		  if (!no_image_save) {
			  eyesisCorrections.saveAndShow(
					  imp_texture_bgnd_ext,
					  path,
					  correctionsParameters.png  && !clt_parameters.black_back,
					  clt_parameters.show_textures,
					  jpegQuality); // jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
		  }
		  return imp_texture_bgnd_ext;
	  }



	 public String getPassImage( // get image form a single pass, return relative path for x3d
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  String     name,
			  int        scanIndex,
			  int        threadsMax,  // maximal number of threads to launch
			  boolean    updateStatus,
			  int        debugLevel)
	  {
		 final int tilesX = tp.getTilesX();
		 final int tilesY = tp.getTilesY();

		  showDoubleFloatArrays sdfa_instance = null;
    	  if (debugLevel > -1) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		  CLTPass3d scan = tp.clt_3d_passes.get(scanIndex);
		  boolean [] borderTiles = scan.border_tiles;
		  double [][][][] texture_tiles = scan.texture_tiles;
		  scan.updateSelection(); // update .selected field (all selected, including border) and Rectangle bounds
		  double [][]alphaFade = tp.getAlphaFade(clt_parameters.transform_size);
		  if ((debugLevel > 0) && (scanIndex == 1)) {
			  String [] titles = new String[16];
			  for (int i = 0; i<titles.length;i++)  titles[i]=""+i;
			  sdfa_instance.showArrays(alphaFade, 2*clt_parameters.transform_size,2*clt_parameters.transform_size,true,"alphaFade",titles);
		  }
		  double [][][][] texture_tiles_cluster = new double[tilesY][tilesX][][];
		  double [] alpha_zero = new double [4*clt_parameters.transform_size*clt_parameters.transform_size];
		  int alpha_index = 3;
		  for (int i = 0; i < alpha_zero.length; i++) alpha_zero[i]=0.0;
		  for (int tileY = 0; tileY < tilesY; tileY++){
			  for (int tileX = 0; tileX < tilesX; tileX++){
				  texture_tiles_cluster[tileY][tileX]= null;
				  if (texture_tiles[tileY][tileX] != null) {
					  if (borderTiles[tileY * tilesX + tileX]) {
						  texture_tiles_cluster[tileY][tileX]= texture_tiles[tileY][tileX].clone();
						  if (clt_parameters.shAggrFade) {
							  texture_tiles_cluster[tileY][tileX][alpha_index] = alpha_zero;
						  } else {
							  if ((debugLevel > -1) && (scanIndex == 1)) {
								  System.out.println("getPassImage(): tileY="+tileY+", tileX = "+tileX+", tileY="+tileY);
							  }
							  int fade_mode=0;
							  if ((tileY > 0) &&              (texture_tiles[tileY - 1][tileX] != null) && !borderTiles[(tileY - 1) * tilesX + tileX]) fade_mode |= 1;
							  if ((tileX < (tilesX -1)) && (texture_tiles[tileY][tileX + 1] != null) && !borderTiles[tileY * tilesX + tileX + 1])   fade_mode |= 2;
							  if ((tileY < (tilesY -1)) && (texture_tiles[tileY + 1][tileX] != null) && !borderTiles[(tileY + 1) * tilesX + tileX]) fade_mode |= 4;
							  if ((tileX > 0) &&              (texture_tiles[tileY][tileX - 1] != null) && !borderTiles[tileY * tilesX + tileX - 1])   fade_mode |= 8;
							  texture_tiles_cluster[tileY][tileX][alpha_index] = alphaFade[fade_mode]; // alpha_zero;
						  }
					  }else{
						  texture_tiles_cluster[tileY][tileX]= texture_tiles[tileY][tileX];
					  }
				  }
			  }
		  }

		  ImageDtt image_dtt = new ImageDtt();
		  double [][] texture_overlap = image_dtt.combineRGBATiles(
				  texture_tiles_cluster, // texture_tiles,               // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}
				  clt_parameters.transform_size,
				  true,                         // when false - output each tile as 16x16, true - overlap to make 8x8
				  clt_parameters.sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
				  threadsMax,                    // maximal number of threads to launch
				  debugLevel);
		  if (clt_parameters.alpha1 > 0){ // negative or 0 - keep alpha as it was
			  double scale = (clt_parameters.alpha1 > clt_parameters.alpha0) ? (1.0/(clt_parameters.alpha1 - clt_parameters.alpha0)) : 0.0;
			  for (int i = 0; i < texture_overlap[alpha_index].length; i++){
				  double d = texture_overlap[alpha_index][i];
				  if      (d >=clt_parameters.alpha1) d = 1.0;
				  else if (d <=clt_parameters.alpha0) d = 0.0;
				  else d = scale * (d- clt_parameters.alpha0);
				  texture_overlap[alpha_index][i] = d;
			  }
		  }
		  // for now - use just RGB. Later add oprion for RGBA
		  double [][] texture_rgb = {texture_overlap[0],texture_overlap[1],texture_overlap[2]};
		  double [][] texture_rgba = {texture_overlap[0],texture_overlap[1],texture_overlap[2],texture_overlap[3]};
		  double [][] texture_rgbx = ((clt_parameters.alpha1 > 0)? texture_rgba: texture_rgb);

		  boolean resize = true;
		  if (resize) {
		  texture_rgbx = resizeGridTexture(
				  texture_rgbx,
				  clt_parameters.transform_size,
				  tilesX,
				  tilesY,
				  scan.getTextureBounds());
		  }

		  int width = resize ? (clt_parameters.transform_size * scan.getTextureBounds().width): (clt_parameters.transform_size * tilesX);
		  int height = resize ? (clt_parameters.transform_size * scan.getTextureBounds().height): (clt_parameters.transform_size * tilesY);

		  ImagePlus imp_texture_cluster = linearStackToColor(
				  clt_parameters,
				  colorProcParameters,
				  rgbParameters,
				  name+"-texture", // String name,
				  "", //String suffix, // such as disparity=...
				  true, // toRGB,
				  !this.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
				  true, // boolean saveShowIntermediate, // save/show if set globally
				  false, //true, // boolean saveShowFinal,        // save/show result (color image?)
				  texture_rgbx,
				  width, //tp.tilesX *  clt_parameters.transform_size,
				  height, //tp.tilesY *  clt_parameters.transform_size,
				  1.0,         // double scaleExposure, // is it needed?
				  debugLevel);


		  String path= correctionsParameters.selectX3dDirectory(
				  //TODO: Which one to use - name or this.image_name ?
 				  correctionsParameters.getModelName(this.image_name), // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
 				  correctionsParameters.x3dModelVersion,
				  true,  // smart,
				  true);  //newAllowed, // save
		  // only show/save original size if debug or debug_filters)
			  eyesisCorrections.saveAndShow(
					  imp_texture_cluster,
					  path,
					  correctionsParameters.png,
					  clt_parameters.show_textures,
					  -1); // jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
		  return imp_texture_cluster.getTitle()+".png"; // imp_texture_cluster;
	  }



	  public ImagePlus resizeForBackdrop(
			  ImagePlus imp,
			  boolean fillBlack,
			  boolean noalpha, // only with fillBlack, otherwize ignored
			  int debugLevel)
	  {
		  double backdropPixels = 2.0/geometryCorrection.getFOVPix();
		  if (debugLevel > -1) {
			  System.out.println("backdropPixels = "+backdropPixels);
		  }
		  // TODO: currently - just adding pixels, no rescaling (add later). Alternatively - just modify geometry earlier
		  int width = imp.getWidth();
		  int height = imp.getHeight();
		  int h_margin = (int) Math.round((backdropPixels - width)/2);
		  int v_margin = (int) Math.round((backdropPixels - height)/2);
		  int width2 =  width +  2 * h_margin;
		  int height2 = height + 2 * v_margin;
		  if (debugLevel > -1) {
			  System.out.println("backdropPixels = "+backdropPixels+" h_margin = "+h_margin+" v_margin = "+v_margin);
		  }
		  int [] src_pixels = (int []) imp.getProcessor().getPixels();
		  int [] pixels = new int [width2* height2];
		  int black = noalpha ? 0 :        0xff000000;
		  int mask =  noalpha ? 0xffffff : 0xffffffff;
		  if (fillBlack) {
			  for (int i = 0; i < pixels.length; i++){
				  pixels[i] = black;
			  }
		  }
		  int indx = 0;
		  int offset = v_margin *  width2 + h_margin;
		  if (fillBlack) {
			  for (int i = 0; i < height; i++){
				  for (int j = 0; j < width; j++){
					  int a =  (src_pixels[indx] >> 24) & 0xff;
					  if (a == 255) {
						  pixels[offset+ i * width2 + j] = src_pixels[indx] & mask;
					  } else if (a == 0) {
						  pixels[offset+ i * width2 + j] = black;
					  } else  {
						  int r = (src_pixels[indx] >> 16) & 0xff; //' maybe - swap with b?
						  int g = (src_pixels[indx] >>  8) & 0xff;
						  int b = (src_pixels[indx] >>  0) & 0xff;
						  r = (r * a) / 255;
						  g = (g * a) / 255;
						  b = (b * a) / 255;
						  pixels[offset+ i * width2 + j] = black | (r << 16) | (g << 8) | b;
					  }
					  indx++;
				  }
			  }
		  } else {
			  for (int i = 0; i < height; i++){
				  for (int j = 0; j < width; j++){
					  pixels[offset+ i * width2 + j] = src_pixels[indx++];
				  }
			  }
		  }
		  ColorProcessor cp=new ColorProcessor(width2,height2);
		  cp.setPixels(pixels);
		  ImagePlus imp_ext=new ImagePlus(imp.getTitle()+"-ext",cp);
		  return imp_ext;
	  }

//[tp.tilesY][tp.tilesX]["RGBA".length()][]
//linearStackToColor

	  public CLTPass3d CLTBackgroundMeas( // measure background
//			  final String        image_name,
			  final double [][][] image_data, // first index - number of image in a quad
			  final boolean [][]  saturation_imp, // (near) saturated pixels or null
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  final int           threadsMax,  // maximal number of threads to launch
			  final boolean       updateStatus,
			  final int           debugLevel)
	  {
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
//		  CLTPass3d scan_rslt = tp.new CLTPass3d(tp);
		  CLTPass3d scan_rslt = new CLTPass3d(tp);
		  int d = ImageDtt.setImgMask(0, 0xf);
		  d =     ImageDtt.setPairMask(d,0xf);
		  d =     ImageDtt.setForcedDisparity(d,true);
		  int [][]     tile_op =         tp.setSameTileOp(clt_parameters,  d, debugLevel);
		  double [][]  disparity_array = tp.setSameDisparity(0.0); // [tp.tilesY][tp.tilesX] - individual per-tile expected disparity
		  // undecided, so 2 modes of combining alpha - same as rgb, or use center tile only
		  double [][][][]     clt_corr_combo =    new double [ImageDtt.TCORR_TITLES.length][tilesY][tilesX][]; // will only be used inside?

//		  double min_corr_selected = clt_parameters.corr_normalize? clt_parameters.min_corr_normalized: clt_parameters.min_corr;
		  double min_corr_selected = clt_parameters.min_corr;

		  double [][] disparity_map = new double [ImageDtt.DISPARITY_TITLES.length][]; //[0] -residual disparity, [1] - orthogonal (just for debugging)

		  double [][] shiftXY = new double [4][2];
		  if (!clt_parameters.fine_corr_ignore) {
			  double [][] shiftXY0 = {
					  {clt_parameters.fine_corr_x_0,clt_parameters.fine_corr_y_0},
					  {clt_parameters.fine_corr_x_1,clt_parameters.fine_corr_y_1},
					  {clt_parameters.fine_corr_x_2,clt_parameters.fine_corr_y_2},
					  {clt_parameters.fine_corr_x_3,clt_parameters.fine_corr_y_3}};
			  shiftXY = shiftXY0;
		  }

		  double [][][][] texture_tiles =     new double [tilesY][tilesX][][]; // ["RGBA".length()][];
		  ImageDtt image_dtt = new ImageDtt();
//		  final double disparity_corr = (clt_parameters.z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/clt_parameters.z_correction);
		  double z_correction =  clt_parameters.z_correction;
		  if (clt_parameters.z_corr_map.containsKey(image_name)){
			  z_correction +=clt_parameters.z_corr_map.get(image_name);
		  }
		  final double disparity_corr = (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);

		  image_dtt.clt_aberrations_quad_corr(
				  clt_parameters.img_dtt,       // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
				  1,                            // final int  macro_scale, // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
				  tile_op,                      // per-tile operation bit codes
				  disparity_array,              // clt_parameters.disparity,     // final double            disparity,
				  image_data,                   // final double [][][]      imade_data, // first index - number of image in a quad
				  saturation_imp,               // boolean [][] saturation_imp, // (near) saturated pixels or null
				  // correlation results - final and partial
				  clt_corr_combo,               // [tp.tilesY][tp.tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				  null, // clt_corr_partial,    // [tp.tilesY][tp.tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				  null,    // [tp.tilesY][tp.tilesX][pair]{dx,dy,weight}[(2*transform_size-1)*(2*transform_size-1)] // transpose unapplied. null - do not calculate
				  //	Use it with disparity_maps[scan_step]?		  clt_mismatch,    // [tp.tilesY][tp.tilesX][pair]{dx,dy,weight}[(2*transform_size-1)*(2*transform_size-1)] // transpose unapplied. null - do not calculate
				  disparity_map,    // [12][tp.tilesY * tp.tilesX]
				  texture_tiles,        // [tp.tilesY][tp.tilesX]["RGBA".length()][];
				  tilesX * clt_parameters.transform_size, // imp_quad[0].getWidth(),       // final int width,
				  clt_parameters.fat_zero,      // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
				  clt_parameters.corr_sym,
				  clt_parameters.corr_offset,
				  clt_parameters.corr_red,
				  clt_parameters.corr_blue,
				  clt_parameters.corr_sigma,
				  clt_parameters.corr_normalize, // normalize correlation results by rms
				  min_corr_selected, // 0.0001; // minimal correlation value to consider valid
				  clt_parameters.max_corr_sigma,// 1.5;  // weights of points around global max to find fractional
				  clt_parameters.max_corr_radius,
//				  clt_parameters.enhortho_width,  // 2;    // reduce weight of center correlation pixels from center (0 - none, 1 - center, 2 +/-1 from center)
//				  clt_parameters.enhortho_scale,  // 0.2;  // multiply center correlation pixels (inside enhortho_width)
				  clt_parameters.max_corr_double, // Double pass when masking center of mass to reduce preference for integer values
				  clt_parameters.corr_mode,     // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
				  clt_parameters.min_shot,       // 10.0;  // Do not adjust for shot noise if lower than
				  clt_parameters.scale_shot,     // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
				  clt_parameters.diff_sigma,     // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_threshold, // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_gauss,     // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
				  clt_parameters.min_agree,      // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				  clt_parameters.dust_remove,    // Do not reduce average weight when only one image differes much from the average
				  clt_parameters.keep_weights,   // Add port weights to RGBA stack (debug feature)
				  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
				  null,                          // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
				  clt_kernels,                   // final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				  clt_parameters.kernel_step,
				  clt_parameters.transform_size,
				  clt_parameters.clt_window,
				  shiftXY, //
				  disparity_corr, // final double              disparity_corr, // disparity at infinity
				  (clt_parameters.fcorr_ignore? null: this.fine_corr),
				  clt_parameters.corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
				  clt_parameters.shift_x,       // final int               shiftX, // shift image horizontally (positive - right) - just for testing
				  clt_parameters.shift_y,       // final int               shiftY, // shift image vertically (positive - down)
				  clt_parameters.tileX,         // final int               debug_tileX,
				  clt_parameters.tileY,         // final int               debug_tileY,
				  (clt_parameters.dbg_mode & 64) != 0, // no fract shift
				  (clt_parameters.dbg_mode & 128) != 0, // no convolve
				  //				  (clt_parameters.dbg_mode & 256) != 0, // transpose convolve
				  threadsMax,
				  debugLevel);

		  scan_rslt.disparity = disparity_array;
		  scan_rslt.tile_op = tile_op;
		  scan_rslt.disparity_map = disparity_map;
		  scan_rslt.texture_tiles = texture_tiles;
		  scan_rslt.is_measured =   true;
		  scan_rslt.is_combo =      false;
		  scan_rslt.resetProcessed();

		  return scan_rslt;
	  }

	  public CLTPass3d  CLTMeasure( // perform single pass according to prepared tiles operations and disparity
			  final double [][][]       image_data, // first index - number of image in a quad
			  final boolean [][]  saturation_imp, // (near) saturated pixels or null
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  final int         scanIndex,
			  final boolean     save_textures,
			  final int         threadsMax,  // maximal number of threads to launch
			  final boolean     updateStatus,
			  final int         debugLevel)
	  {
		  return CLTMeasure( // perform single pass according to prepared tiles operations and disparity
				  image_data,     //  final double [][][]       image_data, // first index - number of image in a quad
				  saturation_imp, // boolean [][] saturation_imp, // (near) saturated pixels or null
				  clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
				  scanIndex,      // final int         scanIndex,
				  save_textures,  // final boolean     save_textures,
				  true,      // final boolean       save_corr,
				  null,           // final double [][] mismatch,    // null or double [12][]
				  threadsMax,     // final int         threadsMax,  // maximal number of threads to launch
				  updateStatus,   // final boolean     updateStatus,
				  debugLevel);    // final int         debugLevel);
	  }
	  public CLTPass3d  CLTMeasure( // perform single pass according to prepared tiles operations and disparity
			  final double [][][]       image_data, // first index - number of image in a quad
			  final boolean [][]  saturation_imp, // (near) saturated pixels or null
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  final int         scanIndex,
			  final boolean     save_textures,
			  final boolean     save_corr,
			  final int         threadsMax,  // maximal number of threads to launch
			  final boolean     updateStatus,
			  final int         debugLevel)
	  {
		  return CLTMeasure( // perform single pass according to prepared tiles operations and disparity
				  image_data,     //  final double [][][]       image_data, // first index - number of image in a quad
				  saturation_imp, // boolean [][] saturation_imp, // (near) saturated pixels or null
				  clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
				  scanIndex,      // final int         scanIndex,
				  save_textures,  // final boolean     save_textures,
				  save_corr,      // final boolean       save_corr,
				  null,           // final double [][] mismatch,    // null or double [12][]
				  threadsMax,     // final int         threadsMax,  // maximal number of threads to launch
				  updateStatus,   // final boolean     updateStatus,
				  debugLevel);    // final int         debugLevel);
	  }
	  public CLTPass3d  CLTMeasure( // perform single pass according to prepared tiles operations and disparity
			  final double [][][]       image_data, // first index - number of image in a quad
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  final int         scanIndex,
			  final boolean     save_textures,
			  final int         threadsMax,  // maximal number of threads to launch
			  final boolean     updateStatus,
			  final int         debugLevel)
	  {
		  return CLTMeasure( // perform single pass according to prepared tiles operations and disparity
				  image_data,     //  final double [][][]       image_data, // first index - number of image in a quad
				  null,           // boolean [][] saturation_imp, // (near) saturated pixels or null
				  clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
				  scanIndex,      // final int         scanIndex,
				  save_textures,  // final boolean     save_textures,
				  true,      // final boolean       save_corr,
				  null,           // final double [][] mismatch,    // null or double [12][]
				  threadsMax,     // final int         threadsMax,  // maximal number of threads to launch
				  updateStatus,   // final boolean     updateStatus,
				  debugLevel);    // final int         debugLevel);
	  }

	  public CLTPass3d  CLTMeasure( // perform single pass according to prepared tiles operations and disparity
			  final double [][][]       image_data, // first index - number of image in a quad
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  final int         scanIndex,
			  final boolean     save_textures,
			  final boolean     save_corr,
			  final int         threadsMax,  // maximal number of threads to launch
			  final boolean     updateStatus,
			  final int         debugLevel)
	  {
		  return CLTMeasure( // perform single pass according to prepared tiles operations and disparity
				  image_data,     //  final double [][][]       image_data, // first index - number of image in a quad
				  null,           // boolean [][] saturation_imp, // (near) saturated pixels or null
				  clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
				  scanIndex,      // final int         scanIndex,
				  save_textures,  // final boolean     save_textures,
				  save_corr,      // final boolean       save_corr,
				  null,           // final double [][] mismatch,    // null or double [12][]
				  threadsMax,     // final int         threadsMax,  // maximal number of threads to launch
				  updateStatus,   // final boolean     updateStatus,
				  debugLevel);    // final int         debugLevel);
	  }

	  public CLTPass3d  CLTMeasure( // perform single pass according to prepared tiles operations and disparity
//			  final String        image_name,
			  final double [][][] image_data, // first index - number of image in a quad
			  final boolean [][]  saturation_imp, // (near) saturated pixels or null
			  final EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  final int           scanIndex,
			  final boolean       save_textures,
			  final boolean       save_corr,
			  final double [][]   mismatch,    // null or double [12][]
			  final int           threadsMax,  // maximal number of threads to launch
			  final boolean       updateStatus,
			  final int           debugLevel)
	  {
		  final int dbg_x = -295-debugLevel;
		  final int dbg_y = -160-debugLevel;
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  CLTPass3d scan = tp.clt_3d_passes.get(scanIndex);
		  int [][]     tile_op =         scan.tile_op;


		  double [][]  disparity_array = scan.disparity;
		  // undecided, so 2 modes of combining alpha - same as rgb, or use center tile only
		  double [][][][]     clt_corr_combo =    new double [ImageDtt.TCORR_TITLES.length][tilesY][tilesX][]; // will only be used inside?
		  if (debugLevel > -1){
			  int numTiles = 0;
			  for (int ty = 0; ty < tile_op.length; ty ++) for (int tx = 0; tx < tile_op[ty].length; tx ++){
				  if (tile_op[ty][tx] != 0) numTiles ++;
			  }
			  System.out.println("CLTMeasure("+scanIndex+"): numTiles = "+numTiles);
			  if ((dbg_y >= 0) && (dbg_x >= 0) && (tile_op[dbg_y][dbg_x] != 0)){
				  System.out.println("CLTMeasure("+scanIndex+"): tile_op["+dbg_y+"]["+dbg_x+"] = "+tile_op[dbg_y][dbg_x]);
			  }
		  }
		  double min_corr_selected = clt_parameters.min_corr;

		  double [][] disparity_map = save_corr ? new double [ImageDtt.DISPARITY_TITLES.length][] : null; //[0] -residual disparity, [1] - orthogonal (just for debugging)

		  double [][] shiftXY = new double [4][2];
		  if (!clt_parameters.fine_corr_ignore) {
			  double [][] shiftXY0 = {
					  {clt_parameters.fine_corr_x_0,clt_parameters.fine_corr_y_0},
					  {clt_parameters.fine_corr_x_1,clt_parameters.fine_corr_y_1},
					  {clt_parameters.fine_corr_x_2,clt_parameters.fine_corr_y_2},
					  {clt_parameters.fine_corr_x_3,clt_parameters.fine_corr_y_3}};
			  shiftXY = shiftXY0;
		  }

		  double [][][][] texture_tiles =   save_textures ? new double [tilesY][tilesX][][] : null; // ["RGBA".length()][];
		  ImageDtt image_dtt = new ImageDtt();
//		  final double disparity_corr = (clt_parameters.z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/clt_parameters.z_correction);
		  double z_correction =  clt_parameters.z_correction;
		  if (clt_parameters.z_corr_map.containsKey(image_name)){
			  z_correction +=clt_parameters.z_corr_map.get(image_name);
		  }
		  final double disparity_corr = (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);

		  image_dtt.clt_aberrations_quad_corr(
				  clt_parameters.img_dtt,       // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
				  1,                            // final int  macro_scale, // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
				  tile_op,                      // per-tile operation bit codes
				  disparity_array,              // clt_parameters.disparity,     // final double            disparity,
				  image_data,                   // final double [][][]      imade_data, // first index - number of image in a quad
				  saturation_imp,               // boolean [][] saturation_imp, // (near) saturated pixels or null
				  // correlation results - final and partial
				  clt_corr_combo,               // [tp.tilesY][tp.tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				  null, // clt_corr_partial,    // [tp.tilesY][tp.tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				  mismatch, // null,    // [tp.tilesY][tp.tilesX][pair]{dx,dy,weight}[(2*transform_size-1)*(2*transform_size-1)] // transpose unapplied. null - do not calculate
				  //	Use it with disparity_maps[scan_step]?		  clt_mismatch,    // [tp.tilesY][tp.tilesX][pair]{dx,dy,weight}[(2*transform_size-1)*(2*transform_size-1)] // transpose unapplied. null - do not calculate
				  disparity_map,    // [12][tp.tilesY * tp.tilesX]
				  texture_tiles,        // [tp.tilesY][tp.tilesX]["RGBA".length()][];
				  tilesX * clt_parameters.transform_size, // imp_quad[0].getWidth(),       // final int width,
				  clt_parameters.fat_zero,      // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
				  clt_parameters.corr_sym,
				  clt_parameters.corr_offset,
				  clt_parameters.corr_red,
				  clt_parameters.corr_blue,
				  clt_parameters.corr_sigma,
				  clt_parameters.corr_normalize, // normalize correlation results by rms
				  min_corr_selected, // 0.0001; // minimal correlation value to consider valid
				  clt_parameters.max_corr_sigma,// 1.5;  // weights of points around global max to find fractional
				  clt_parameters.max_corr_radius,
//				  clt_parameters.enhortho_width,  // 2;    // reduce weight of center correlation pixels from center (0 - none, 1 - center, 2 +/-1 from center)
//				  clt_parameters.enhortho_scale,  // 0.2;  // multiply center correlation pixels (inside enhortho_width)
				  clt_parameters.max_corr_double, // Double pass when masking center of mass to reduce preference for integer values
				  clt_parameters.corr_mode,     // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
				  clt_parameters.min_shot,       // 10.0;  // Do not adjust for shot noise if lower than
				  clt_parameters.scale_shot,     // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
				  clt_parameters.diff_sigma,     // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_threshold, // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_gauss,     // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
				  clt_parameters.min_agree,      // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				  clt_parameters.dust_remove,    // Do not reduce average weight when only one image differes much from the average
				  clt_parameters.keep_weights,   // Add port weights to RGBA stack (debug feature)
				  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
				  null,                          // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
				  clt_kernels,                   // final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				  clt_parameters.kernel_step,
				  clt_parameters.transform_size,
				  clt_parameters.clt_window,
				  shiftXY, //
				  disparity_corr, // final double              disparity_corr, // disparity at infinity
				  (clt_parameters.fcorr_ignore? null: this.fine_corr),
				  clt_parameters.corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
				  clt_parameters.shift_x,       // final int               shiftX, // shift image horizontally (positive - right) - just for testing
				  clt_parameters.shift_y,       // final int               shiftY, // shift image vertically (positive - down)
				  clt_parameters.tileX,         // final int               debug_tileX,
				  clt_parameters.tileY,         // final int               debug_tileY,
				  (clt_parameters.dbg_mode & 64) != 0, // no fract shift
				  (clt_parameters.dbg_mode & 128) != 0, // no convolve
				  //				  (clt_parameters.dbg_mode & 256) != 0, // transpose convolve
				  threadsMax,
				  debugLevel);

		  scan.disparity_map = disparity_map;
		  scan.texture_tiles = texture_tiles;
		  scan.is_measured =   true;
		  scan.is_combo =      false;
		  scan.resetProcessed();
		  return scan;
	  }


	  public CLTPass3d  CLTMeasure( // perform single pass according to prepared tiles operations and disparity
			  final double [][][] image_data, // first index - number of image in a quad
			  final boolean [][]  saturation_imp, // (near) saturated pixels or null
			  final EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  final CLTPass3d     scan,
			  final boolean       save_textures,
			  final boolean       save_corr,
			  final double [][]   mismatch,    // null or double [12][]
			  final GeometryCorrection geometryCorrection_main, // If not null - covert to main camera coordinates
			  final int           threadsMax,  // maximal number of threads to launch
			  final boolean       updateStatus,
			  final int           debugLevel)
	  {
		  final int dbg_x = -295-debugLevel;
		  final int dbg_y = -160-debugLevel;
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  double [] disparity = scan.getDisparity();
		  double [] strength =  scan.getStrength();
		  boolean [] selection = scan.getSelected();
		  if (selection == null) {
			  selection = new boolean[tilesX*tilesY];
			  for (int nTile = 0; nTile < selection.length; nTile++) {
				  selection[nTile] = !Double.isNaN(disparity[nTile]) && (strength[nTile] > 0.0);
			  }
			  scan.setSelected(selection);
		  }
		  if ((scan.disparity == null) || (scan.tile_op == null)) {
			  scan.setTileOpDisparity(
					  scan.getSelected(), // boolean [] selection,
					  scan.getDisparity()); // double []  disparity)
		  }

		  int [][]     tile_op =         scan.tile_op;
		  double [][]  disparity_array = scan.disparity;
		  // undecided, so 2 modes of combining alpha - same as rgb, or use center tile only
		  double [][][][]     clt_corr_combo =    new double [ImageDtt.TCORR_TITLES.length][tilesY][tilesX][]; // will only be used inside?
		  if (debugLevel > -1){
			  int numTiles = 0;
			  for (int ty = 0; ty < tile_op.length; ty ++) for (int tx = 0; tx < tile_op[ty].length; tx ++){
				  if (tile_op[ty][tx] != 0) numTiles ++;
			  }
			  System.out.println("CLTMeasure(): numTiles = "+numTiles);
			  if ((dbg_y >= 0) && (dbg_x >= 0) && (tile_op[dbg_y][dbg_x] != 0)){
				  System.out.println("CLTMeasure(): tile_op["+dbg_y+"]["+dbg_x+"] = "+tile_op[dbg_y][dbg_x]);
			  }
		  }
		  double min_corr_selected = clt_parameters.min_corr;
		  double [][] disparity_map = save_corr ? new double [ImageDtt.DISPARITY_TITLES.length][] : null; //[0] -residual disparity, [1] - orthogonal (just for debugging)
		  double [][] shiftXY = new double [4][2];
		  if (!clt_parameters.fine_corr_ignore) {
			  double [][] shiftXY0 = {
					  {clt_parameters.fine_corr_x_0,clt_parameters.fine_corr_y_0},
					  {clt_parameters.fine_corr_x_1,clt_parameters.fine_corr_y_1},
					  {clt_parameters.fine_corr_x_2,clt_parameters.fine_corr_y_2},
					  {clt_parameters.fine_corr_x_3,clt_parameters.fine_corr_y_3}};
			  shiftXY = shiftXY0;
		  }

		  double [][][][] texture_tiles =   save_textures ? new double [tilesY][tilesX][][] : null; // ["RGBA".length()][];
		  ImageDtt image_dtt = new ImageDtt();
		  double z_correction =  clt_parameters.z_correction;
		  if (clt_parameters.z_corr_map.containsKey(image_name)){
			  z_correction +=clt_parameters.z_corr_map.get(image_name);
		  }
		  final double disparity_corr = (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);

		  image_dtt.clt_aberrations_quad_corr(
				  clt_parameters.img_dtt,       // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
				  1,                            // final int  macro_scale, // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
				  tile_op,                      // per-tile operation bit codes
				  disparity_array,              // clt_parameters.disparity,     // final double            disparity,
				  image_data,                   // final double [][][]      imade_data, // first index - number of image in a quad
				  saturation_imp,               // boolean [][] saturation_imp, // (near) saturated pixels or null
				  // correlation results - final and partial
				  clt_corr_combo,               // [tp.tilesY][tp.tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				  null, // clt_corr_partial,    // [tp.tilesY][tp.tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				  mismatch, // null,    // [tp.tilesY][tp.tilesX][pair]{dx,dy,weight}[(2*transform_size-1)*(2*transform_size-1)] // transpose unapplied. null - do not calculate
				  //	Use it with disparity_maps[scan_step]?		  clt_mismatch,    // [tp.tilesY][tp.tilesX][pair]{dx,dy,weight}[(2*transform_size-1)*(2*transform_size-1)] // transpose unapplied. null - do not calculate
				  disparity_map,    // [12][tp.tilesY * tp.tilesX]
				  texture_tiles,        // [tp.tilesY][tp.tilesX]["RGBA".length()][];
				  tilesX * clt_parameters.transform_size, // imp_quad[0].getWidth(),       // final int width,
				  clt_parameters.fat_zero,      // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
				  clt_parameters.corr_sym,
				  clt_parameters.corr_offset,
				  clt_parameters.corr_red,
				  clt_parameters.corr_blue,
				  clt_parameters.corr_sigma,
				  clt_parameters.corr_normalize, // normalize correlation results by rms
				  min_corr_selected, // 0.0001; // minimal correlation value to consider valid
				  clt_parameters.max_corr_sigma,// 1.5;  // weights of points around global max to find fractional
				  clt_parameters.max_corr_radius,
				  clt_parameters.max_corr_double, // Double pass when masking center of mass to reduce preference for integer values
				  clt_parameters.corr_mode,     // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
				  clt_parameters.min_shot,       // 10.0;  // Do not adjust for shot noise if lower than
				  clt_parameters.scale_shot,     // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
				  clt_parameters.diff_sigma,     // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_threshold, // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_gauss,     // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
				  clt_parameters.min_agree,      // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				  clt_parameters.dust_remove,    // Do not reduce average weight when only one image differes much from the average
				  clt_parameters.keep_weights,   // Add port weights to RGBA stack (debug feature)
				  geometryCorrection,           // final GeometryCorrection  geometryCorrection,
				  geometryCorrection_main, // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
				  clt_kernels,                  // final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				  clt_parameters.kernel_step,
				  clt_parameters.transform_size,
				  clt_parameters.clt_window,
				  shiftXY, //
				  disparity_corr, // final double              disparity_corr, // disparity at infinity
				  (clt_parameters.fcorr_ignore? null: this.fine_corr),
				  clt_parameters.corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
				  clt_parameters.shift_x,       // final int               shiftX, // shift image horizontally (positive - right) - just for testing
				  clt_parameters.shift_y,       // final int               shiftY, // shift image vertically (positive - down)
				  clt_parameters.tileX,         // final int               debug_tileX,
				  clt_parameters.tileY,         // final int               debug_tileY,
				  (clt_parameters.dbg_mode & 64) != 0, // no fract shift
				  (clt_parameters.dbg_mode & 128) != 0, // no convolve
				  //				  (clt_parameters.dbg_mode & 256) != 0, // transpose convolve
				  threadsMax,
				  debugLevel);

		  scan.disparity_map = disparity_map;
		  scan.texture_tiles = texture_tiles;
		  scan.is_measured =   true;
		  scan.is_combo =      false;
		  scan.resetProcessed();
		  return scan;
	  }




	  public ImagePlus [] conditionImageSetBatch( // used in batchCLT3d
			  final int                           nSet, // index of the 4-image set
			  final EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  final int [][]                      fileIndices, // =new int [numImagesToProcess][2]; // file index, channel number
			  final ArrayList<String>             setNames, //  = new ArrayList<String>();
			  final ArrayList<ArrayList<Integer>> setFiles, //  = new ArrayList<ArrayList<Integer>>();
			  final double []                     referenceExposures, // =eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
			  final double []                     scaleExposures, //  = new double[channelFiles.length]; //
			  final boolean [][]                  saturation_imp, //  = (clt_parameters.sat_level > 0.0)? new boolean[channelFiles.length][] : null;
			  final int                           debugLevel)
	  {
		  final boolean                       batch_mode = clt_parameters.batch_run; //disable any debug images
		  String [] sourceFiles=correctionsParameters.getSourcePaths();
		  int maxChn = 0;
		  for (int i = 0; i < setFiles.get(nSet).size(); i++){
			  int chn = fileIndices[setFiles.get(nSet).get(i)][1];
			  if (chn > maxChn) maxChn = chn;
		  }
		  int [] channelFiles = new int[maxChn+1];
		  for (int i =0; i < channelFiles.length; i++) channelFiles[i] = -1;
		  for (int i = 0; i < setFiles.get(nSet).size(); i++){
			  channelFiles[fileIndices[setFiles.get(nSet).get(i)][1]] = setFiles.get(nSet).get(i);
		  }

		  ImagePlus [] imp_srcs = new ImagePlus[channelFiles.length];
///		  boolean [][] saturation_imp = (clt_parameters.sat_level > 0.0)? new boolean[channelFiles.length][] : null;
///		  double [] scaleExposures = new double[channelFiles.length]; //
		  double [][] dbg_dpixels =  batch_mode? null : (new double [channelFiles.length][]);

		  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
			  int nFile=channelFiles[srcChannel];
			  imp_srcs[srcChannel]=null;
			  if (nFile >=0){
				  if (correctionsParameters.isJP4()){
					  int subchannel=eyesisCorrections.pixelMapping.getSubChannel(srcChannel);
					  if (this.correctionsParameters.swapSubchannels01) {
						  switch (subchannel){
						  case 0: subchannel=1; break;
						  case 1: subchannel=0; break;
						  }
					  }
					  if (debugLevel>0) System.out.println("Processing set " + setNames.get(nSet)+" channel "+srcChannel+" - subchannel "+subchannel+" of "+sourceFiles[nFile]);
					  ImagePlus imp_composite=eyesisCorrections.JP4_INSTANCE.open(
							  "", // path,
							  sourceFiles[nFile],
							  "",  //arg - not used in JP46 reader
							  true, // un-apply camera color gains
							  null, // new window
							  false); // do not show
					  imp_srcs[srcChannel]=eyesisCorrections.JP4_INSTANCE.demuxImage(imp_composite, subchannel);
					  if (imp_srcs[srcChannel] == null) imp_srcs[srcChannel] = imp_composite; // not a composite image
					  // do we need to add any properties?
				  } else {
					  imp_srcs[srcChannel]=new ImagePlus(sourceFiles[nFile]);
					  //					  (new JP46_Reader_camera(false)).decodeProperiesFromInfo(imp_src); // decode existent properties from info
					  eyesisCorrections.JP4_INSTANCE.decodeProperiesFromInfo(imp_srcs[srcChannel]); // decode existent properties from info
					  if (debugLevel>0) System.out.println("Processing "+sourceFiles[nFile]);
				  }
				  imp_srcs[srcChannel] =  padBayerToFullSize(
						  imp_srcs[srcChannel], // ImagePlus imp_src,
						  eyesisCorrections.pixelMapping.sensors[srcChannel].getSensorWH(),
						  true); // boolean replicate);
				  scaleExposures[srcChannel] = 1.0;
				  if (!Double.isNaN(referenceExposures[nFile]) && (imp_srcs[srcChannel].getProperty("EXPOSURE")!=null)){
					  scaleExposures[srcChannel] = referenceExposures[nFile]/Double.parseDouble((String) imp_srcs[srcChannel].getProperty("EXPOSURE"));
					  if (debugLevel > -1) {
						  System.out.println("Will scale intensity (to compensate for exposure) by  "+scaleExposures[srcChannel]+
								  ", EXPOSURE = "+imp_srcs[srcChannel].getProperty("EXPOSURE"));
					  }
				  }
				  imp_srcs[srcChannel].setProperty("name",    correctionsParameters.getNameFromSourceTiff(sourceFiles[nFile]));
				  imp_srcs[srcChannel].setProperty("channel", srcChannel); // it may already have channel
				  imp_srcs[srcChannel].setProperty("path",    sourceFiles[nFile]); // it may already have channel

				  if (this.correctionsParameters.pixelDefects && (eyesisCorrections.defectsXY!=null)&& (eyesisCorrections.defectsXY[srcChannel]!=null)){
					  // apply pixel correction
					  int numApplied=	eyesisCorrections.correctDefects(
							  imp_srcs[srcChannel],
							  srcChannel,
							  debugLevel);
					  if ((debugLevel>0) && (numApplied>0)) { // reduce verbosity after verified defect correction works
						  System.out.println("Corrected "+numApplied+" pixels in "+sourceFiles[nFile]);
					  }
				  }
				  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
				  int width =  imp_srcs[srcChannel].getWidth();
				  int height = imp_srcs[srcChannel].getHeight();
				  if (!batch_mode && (debugLevel > -1)) {
					  double [] max_pix= {0.0, 0.0, 0.0, 0.0};
//					  for (int y = 0; y < height-1; y+=2){
					  for (int y = 0; y < 499; y+=2){
//						  for (int x = 0; x < width-1; x+=2){
						  for (int x = width/2; x < width-1; x+=2){
							  if (pixels[y*width+x        ] > max_pix[0])  max_pix[0] = pixels[y*width+x        ];
							  if (pixels[y*width+x+      1] > max_pix[1])  max_pix[1] = pixels[y*width+x+      1];
							  if (pixels[y*width+x+width  ] > max_pix[2])  max_pix[2] = pixels[y*width+x+width  ];
							  if (pixels[y*width+x+width+1] > max_pix[3])  max_pix[3] = pixels[y*width+x+width+1];
						  }
					  }
					  System.out.println(String.format("channel %d max_pix[] = %6.2f %6.2f %6.2f %6.2f", srcChannel, max_pix[0], max_pix[1], max_pix[2], max_pix[3]));
					  dbg_dpixels[srcChannel] = new double [pixels.length];
					  for (int i = 0; i < pixels.length; i++) dbg_dpixels[srcChannel][i] = pixels[i];
					  //						  imp_srcs[srcChannel].show();
				  }
				  if (clt_parameters.sat_level > 0.0){
					  double [] saturations = {
							  Double.parseDouble((String) imp_srcs[srcChannel].getProperty("saturation_1")),
							  Double.parseDouble((String) imp_srcs[srcChannel].getProperty("saturation_0")),
							  Double.parseDouble((String) imp_srcs[srcChannel].getProperty("saturation_3")),
							  Double.parseDouble((String) imp_srcs[srcChannel].getProperty("saturation_2"))};
					  saturation_imp[srcChannel] = new boolean[width*height];
					  System.out.println(String.format("channel %d saturations = %6.2f %6.2f %6.2f %6.2f", srcChannel,
							  saturations[0],saturations[1],saturations[2],saturations[3]));
					  double [] scaled_saturations = new double [saturations.length];
					  for (int i = 0; i < scaled_saturations.length; i++){
						  scaled_saturations[i] = saturations[i] * clt_parameters.sat_level;
					  }
					  for (int y = 0; y < height-1; y+=2){
						  for (int x = 0; x < width-1; x+=2){
							  if (pixels[y*width+x        ] > scaled_saturations[0])  saturation_imp[srcChannel][y*width+x        ] = true;
							  if (pixels[y*width+x+      1] > scaled_saturations[1])  saturation_imp[srcChannel][y*width+x      +1] = true;
							  if (pixels[y*width+x+width  ] > scaled_saturations[2])  saturation_imp[srcChannel][y*width+x+width  ] = true;
							  if (pixels[y*width+x+width+1] > scaled_saturations[3])  saturation_imp[srcChannel][y*width+x+width+1] = true;
						  }
					  }
				  }


				  if (this.correctionsParameters.vignetting){
					  if ((eyesisCorrections.channelVignettingCorrection==null) || (srcChannel<0) || (srcChannel>=eyesisCorrections.channelVignettingCorrection.length) || (eyesisCorrections.channelVignettingCorrection[srcChannel]==null)){
						  System.out.println("No vignetting data for channel "+srcChannel);
						  return null;
					  }
					  ///						  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();


					  if (pixels.length!=eyesisCorrections.channelVignettingCorrection[srcChannel].length){
						  System.out.println("Vignetting data for channel "+srcChannel+" has "+eyesisCorrections.channelVignettingCorrection[srcChannel].length+" pixels, image "+sourceFiles[nFile]+" has "+pixels.length);
						  return null;
					  }
					  // TODO: Move to do it once:
					  double min_non_zero = 0.0;
					  for (int i=0;i<pixels.length;i++){
						  double d = eyesisCorrections.channelVignettingCorrection[srcChannel][i];
						  if ((d > 0.0) && ((min_non_zero == 0) || (min_non_zero > d))){
							  min_non_zero = d;
						  }
					  }
					  double max_vign_corr = clt_parameters.vignetting_range*min_non_zero;

					  System.out.println("Vignetting data: channel="+srcChannel+", min = "+min_non_zero);
					  for (int i=0;i<pixels.length;i++){
						  double d = eyesisCorrections.channelVignettingCorrection[srcChannel][i];
						  if (d > max_vign_corr) d = max_vign_corr;
						  pixels[i]*=d;
					  }
					  // Scale here, combine with vignetting later?
					  ///						  int width =  imp_srcs[srcChannel].getWidth();
					  ///						  int height = imp_srcs[srcChannel].getHeight();
					  for (int y = 0; y < height-1; y+=2){
						  for (int x = 0; x < width-1; x+=2){
							  pixels[y*width+x        ] *= clt_parameters.scale_g;
							  pixels[y*width+x+width+1] *= clt_parameters.scale_g;
							  pixels[y*width+x      +1] *= clt_parameters.scale_r;
							  pixels[y*width+x+width  ] *= clt_parameters.scale_b;
						  }
					  }

				  } else { // assuming GR/BG pattern
					  System.out.println("Applying fixed color gain correction parameters: Gr="+
							  clt_parameters.novignetting_r+", Gg="+clt_parameters.novignetting_g+", Gb="+clt_parameters.novignetting_b);
					  ///						  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
					  ///						  int width =  imp_srcs[srcChannel].getWidth();
					  ///						  int height = imp_srcs[srcChannel].getHeight();
					  double kr = clt_parameters.scale_r/clt_parameters.novignetting_r;
					  double kg = clt_parameters.scale_g/clt_parameters.novignetting_g;
					  double kb = clt_parameters.scale_b/clt_parameters.novignetting_b;
					  for (int y = 0; y < height-1; y+=2){
						  for (int x = 0; x < width-1; x+=2){
							  pixels[y*width+x        ] *= kg;
							  pixels[y*width+x+width+1] *= kg;
							  pixels[y*width+x      +1] *= kr;
							  pixels[y*width+x+width  ] *= kb;
						  }
					  }
				  }
			  }
		  }
		  // temporary applying scaleExposures[srcChannel] here, setting it to all 1.0
		  System.out.println("Temporarily applying scaleExposures[] here" );
		  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
			  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
			  for (int i = 0; i < pixels.length; i++){
				  pixels[i] *= scaleExposures[srcChannel];
			  }
			  scaleExposures[srcChannel] = 1.0;
		  }

		  // may need to equalize gains between channels
		  if (clt_parameters.gain_equalize || clt_parameters.colors_equalize){
			  channelGainsEqualize(
					  clt_parameters.gain_equalize,
					  clt_parameters.colors_equalize,
					  clt_parameters.nosat_equalize, // boolean nosat_equalize,
					  channelFiles,
					  imp_srcs,
					  saturation_imp, // boolean[][] saturated,
					  setNames.get(nSet), // just for debug messages == setNames.get(nSet)
					  debugLevel);
		  }


		  if (!batch_mode && (debugLevel > -1) && (saturation_imp != null)){
			  String [] titles = {"chn0","chn1","chn2","chn3"};
			  double [][] dbg_satur = new double [saturation_imp.length] [saturation_imp[0].length];
			  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
				  for (int i = 0; i < saturation_imp[srcChannel].length; i++){
					  dbg_satur[srcChannel][i] = saturation_imp[srcChannel][i]? 1.0 : 0.0;
				  }
			  }
			  int width =  imp_srcs[0].getWidth();
			  int height = imp_srcs[0].getHeight();
			  (new showDoubleFloatArrays()).showArrays(dbg_satur, width, height, true, "Saturated" , titles);

			  if (debugLevel > -1) { // 0){
				  double [][] dbg_dpixels_norm = new double [channelFiles.length][];
				  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
					  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
					  dbg_dpixels_norm[srcChannel] = new double[pixels.length];
					  for (int i = 0; i < pixels.length; i++){
						  dbg_dpixels_norm[srcChannel][i] = pixels[i];
					  }
				  }
				  (new showDoubleFloatArrays()).showArrays(dbg_dpixels, width, height, true, "dpixels" , titles);
				  (new showDoubleFloatArrays()).showArrays(dbg_dpixels_norm, width, height, true, "dpixels_norm" , titles);
				  double [][] dbg_dpixels_split = new double [4 * dbg_dpixels.length][dbg_dpixels[0].length / 4];
				  String [] dbg_titles = {"g1_0","r_0","b_0","g2_0","g1_2","r_1","b_1","g2_1","g1_2","r_2","b_2","g2_2","g1_3","r_3","b_3","g2_3"};
				  for (int srcChn = 0; srcChn < 4; srcChn++) {
					  for (int y = 0; y < height-1; y+=2){
						  for (int x = 0; x < width-1; x+=2){
							  dbg_dpixels_split[ 0 + 4 * srcChn][ y*width/4 +x/2 ] = dbg_dpixels_norm[srcChn][y * width + x             ];
							  dbg_dpixels_split[ 3 + 4 * srcChn][ y*width/4 +x/2 ] = dbg_dpixels_norm[srcChn][y * width + x  + width + 1];
							  dbg_dpixels_split[ 1 + 4 * srcChn][ y*width/4 +x/2 ] = dbg_dpixels_norm[srcChn][y * width + x  +         1];
							  dbg_dpixels_split[ 2 + 4 * srcChn][ y*width/4 +x/2 ] = dbg_dpixels_norm[srcChn][y * width + x  + width    ];
						  }
					  }
				  }
				  (new showDoubleFloatArrays()).showArrays(dbg_dpixels_split, width/2, height/2, true, "dpixels_split" , dbg_titles);
			  }
		  }

		  return imp_srcs;
	  }

	  public void batchCLT3d(
			  TwoQuadCLT       twoQuadCLT, //maybe null in no-rig mode, otherwise may contain rig measurements to be used as infinity ground truth
			  EyesisCorrectionParameters.CLTParameters          clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters      debayerParameters,
			  EyesisCorrectionParameters.ColorProcParameters    colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters          channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters          rgbParameters,
			  EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
//			  int          convolveFFTSize, // 128 - fft size, kernel size should be size/2
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel)
	  {
		  final int        debugLevelInner=clt_parameters.batch_run? -2: debugLevel;
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
					  // removeUnusedSensorData should be off!?
					  channels=this.eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (eyesisCorrections.isChannelEnabled(channels[i])){
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
		  double [] referenceExposures=eyesisCorrections.calcReferenceExposures(debugLevelInner); // multiply each image by this and divide by individual (if not NaN)
		  int [][] fileIndices=new int [numImagesToProcess][2]; // file index, channel number
		  int index=0;
		  for (int nFile=0;nFile<enabledFiles.length;nFile++){
			  if ((sourceFiles[nFile]!=null) && (sourceFiles[nFile].length()>1)) {
				  int [] channels={correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile])};
				  if (correctionsParameters.isJP4()){
					  int subCamera= channels[0]- correctionsParameters.firstSubCamera; // to match those in the sensor files
					  channels=eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (eyesisCorrections.isChannelEnabled(channels[i])){
						  fileIndices[index  ][0]=nFile;
						  fileIndices[index++][1]=channels[i];
					  }
				  }
			  }
		  }
		  ArrayList<String> setNames = new ArrayList<String>();
		  ArrayList<ArrayList<Integer>> setFiles = new ArrayList<ArrayList<Integer>>();

		  for (int iImage=0;iImage<fileIndices.length;iImage++){
			  int nFile=fileIndices[iImage][0];
			  String setName = correctionsParameters.getNameFromSourceTiff(sourceFiles[nFile]);
			  if (!setNames.contains(setName)) {
				  setNames.add(setName);
				  setFiles.add(new ArrayList<Integer>());
			  }
			  setFiles.get(setNames.indexOf(setName)).add(new Integer(nFile));
		  }

		  // enable debug for single-image when clt_batch_dbg1 is on
		  if (correctionsParameters.clt_batch_dbg1 && (setNames.size() < 2)) {
			  clt_parameters.batch_run = false; // disable batch_run for single image if  clt_batch_dbg1 is on
		  }

		  // Do per 4-image set processing
		  int nSet = 0;
		  for (nSet = 0; nSet < setNames.size(); nSet++){
			  if ((nSet > 0) &&(debugLevel > -2)) {
				  System.out.println("Processing set "+(nSet+0)+" (of "+setNames.size()+") finished at "+
				  IJ.d2s(0.000000001*(System.nanoTime()-this.startSetTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
			  }
			  this.startSetTime = System.nanoTime();
			  boolean [][] saturation_imp = (clt_parameters.sat_level > 0.0)? new boolean[QUAD][] : null;
			  double [] scaleExposures = new double[QUAD]; //
			  ImagePlus [] imp_srcs = conditionImageSetBatch(
					  nSet,               // final int                           nSet, // index of the 4-image set
					  clt_parameters,     // final EyesisCorrectionParameters.CLTParameters           clt_parameters,
					  fileIndices,        // final int [][]                      fileIndices, // =new int [numImagesToProcess][2]; // file index, channel number
					  setNames,           // final ArrayList<String>             setNames, //  = new ArrayList<String>();
					  setFiles,           // final ArrayList<ArrayList<Integer>> setFiles, //  = new ArrayList<ArrayList<Integer>>();
					  referenceExposures, //final double []                     referenceExposures, // =eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
					  scaleExposures,     // final double []                     scaleExposures, //  = new double[channelFiles.length]; //
					  saturation_imp,     // final boolean [][]                  saturation_imp, //  = (clt_parameters.sat_level > 0.0)? new boolean[channelFiles.length][] : null;
					  debugLevelInner);   // final int                           debugLevel)

			  // once per quad here
			  if (imp_srcs == null) continue;
			  // creating GeometryCorrection instance for  applyPixelShift()

			  if (correctionsParameters.clt_batch_apply_man) {
				  boolean fine_corr_set = !clt_parameters.fine_corr_ignore;
				  if (fine_corr_set) {
					  boolean nz = false;
					  double [][] shiftXY = {
							  {clt_parameters.fine_corr_x_0,clt_parameters.fine_corr_y_0},
							  {clt_parameters.fine_corr_x_1,clt_parameters.fine_corr_y_1},
							  {clt_parameters.fine_corr_x_2,clt_parameters.fine_corr_y_2},
							  {clt_parameters.fine_corr_x_3,clt_parameters.fine_corr_y_3}};

					  for (int i = 0; i < shiftXY.length; i++){
						  for (int j = 0; j < shiftXY[i].length; j++){
							  if (shiftXY[i][j] != 0.0) {
								  nz = true;
								  break;
							  }
						  }
					  }
					  if (nz) {
						  geometryCorrection.getCorrVector().applyPixelShift(
								  shiftXY); // double [][] pXY_shift)
						  clt_parameters.fine_corr_ignore = true;
						  System.out.println("Detected non-zero manual pixel correction, applying it to extrinsics (azimuth, tilt) and disabling");
					  }
				  }

			  } else if (!clt_parameters.fine_corr_ignore){ // temporary? Remove DC from the manual correction
				  double [][] shiftXY = {
						  {clt_parameters.fine_corr_x_0,clt_parameters.fine_corr_y_0},
						  {clt_parameters.fine_corr_x_1,clt_parameters.fine_corr_y_1},
						  {clt_parameters.fine_corr_x_2,clt_parameters.fine_corr_y_2},
						  {clt_parameters.fine_corr_x_3,clt_parameters.fine_corr_y_3}};
					double [] pXY_avg = {0.0,0.0};
					for (int i = 0; i < shiftXY.length; i++){
						for (int j = 0; j < 2; j++) {
							pXY_avg[j] += shiftXY[i][j]/shiftXY.length;
						}
					}
					for (int i = 0; i < shiftXY.length; i++){
						for (int j = 0; j < 2; j++) {
							shiftXY[i][j] -= pXY_avg[j];
						}
					}
					clt_parameters.fine_corr_x_0 = shiftXY[0][0];
					clt_parameters.fine_corr_y_0 = shiftXY[0][1];
					clt_parameters.fine_corr_x_1 = shiftXY[1][0];
					clt_parameters.fine_corr_y_1 = shiftXY[1][1];
					clt_parameters.fine_corr_x_2 = shiftXY[2][0];
					clt_parameters.fine_corr_y_2 = shiftXY[2][1];
					clt_parameters.fine_corr_x_3 = shiftXY[3][0];
					clt_parameters.fine_corr_y_3 = shiftXY[3][1];
			  }
			  if (correctionsParameters.clt_batch_extrinsic) {
				  if (tp != null) tp.resetCLTPasses();
				  boolean ok = preExpandCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
						  imp_srcs, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
						  saturation_imp, // boolean [][] saturation_imp, // (near) saturated pixels or null
						  clt_parameters,
						  debayerParameters,
//						  nonlinParameters,
						  colorProcParameters,
//						  channelGainParameters,
						  rgbParameters,
						  threadsMax,  // maximal number of threads to launch
						  updateStatus,
						  debugLevelInner);
				  if (ok) {
					  System.out.println("Adjusting extrinsics");
					  extrinsicsCLT(
//							  twoQuadCLT,   // TwoQuadCLT       twoQuadCLT, //maybe null in no-rig mode, otherwise may contain rig measurements to be used as infinity ground truth
							  clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
							  false,          // adjust_poly,
							  threadsMax,     //final int        threadsMax,  // maximal number of threads to launch
							  updateStatus,   // final boolean    updateStatus,
							  debugLevelInner);    // final int        debugLevel)
				  }
			  }
			  if (correctionsParameters.clt_batch_poly) {
				  if (tp != null) tp.resetCLTPasses();
				  boolean ok = preExpandCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
						  imp_srcs, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
						  saturation_imp, // boolean [][] saturation_imp, // (near) saturated pixels or null
						  clt_parameters,
						  debayerParameters,
//						  nonlinParameters,
						  colorProcParameters,
//						  channelGainParameters,
						  rgbParameters,
						  threadsMax,  // maximal number of threads to launch
						  updateStatus,
						  debugLevelInner);
				  if (ok) {
					  System.out.println("Adjusting polynomial fine crorection");
					  extrinsicsCLT(
//							  twoQuadCLT,   // TwoQuadCLT       twoQuadCLT, //maybe null in no-rig mode, otherwise may contain rig measurements to be used as infinity ground truth
							  clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
							  true,           // adjust_poly,
							  threadsMax,     //final int        threadsMax,  // maximal number of threads to launch
							  updateStatus,   // final boolean    updateStatus,
							  debugLevelInner);    // final int        debugLevel)
				  }

			  }
			  if (correctionsParameters.clt_batch_4img){
				  processCLTQuadCorr( // returns ImagePlus, but it already should be saved/shown
						  imp_srcs, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
						  saturation_imp, // boolean [][] saturation_imp, // (near) saturated pixels or null
						  clt_parameters,
						  debayerParameters,
//						  nonlinParameters,
						  colorProcParameters,
						  channelGainParameters,
						  rgbParameters,
//						  convolveFFTSize, // 128 - fft size, kernel size should be size/2
						  scaleExposures,
						  false, // apply_corr, // calculate and apply additional fine geometry correction
						  false, // infinity_corr, // calculate and apply geometry correction at infinity
						  threadsMax,  // maximal number of threads to launch
						  updateStatus,
						  debugLevelInner);
			  }
			  if (correctionsParameters.clt_batch_explore) {
				  if (tp != null) tp.resetCLTPasses();
				  boolean ok = preExpandCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
						  imp_srcs, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
						  saturation_imp, // boolean [][] saturation_imp, // (near) saturated pixels or null
						  clt_parameters,
						  debayerParameters,
//						  nonlinParameters,
						  colorProcParameters,
//						  channelGainParameters,
						  rgbParameters,
						  threadsMax,  // maximal number of threads to launch
						  updateStatus,
						  debugLevelInner);
				  if (ok) {
					  System.out.println("Explore 3d space");
					  expandCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
//							  imp_srcs, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
							  clt_parameters,
							  debayerParameters,
//							  nonlinParameters,
							  colorProcParameters,
							  channelGainParameters,
							  rgbParameters,
							  threadsMax,  // maximal number of threads to launch
							  updateStatus,
							  debugLevelInner);
				  } else continue;

			  } else continue; // if (correctionsParameters.clt_batch_explore)

			  if (correctionsParameters.clt_batch_surf) {
				  tp.showPlanes(
						  clt_parameters,
						  geometryCorrection,
						  threadsMax,
						  updateStatus,
						  debugLevelInner);

			  } else continue; // if (correctionsParameters.clt_batch_surf)

			  if (correctionsParameters.clt_batch_assign) {
				  double [][] assign_dbg = tp.assignTilesToSurfaces(
						  clt_parameters,
						  geometryCorrection,
						  threadsMax,
						  updateStatus,
						  debugLevelInner);
				  if (assign_dbg == null) continue;
			  } else continue; // if (correctionsParameters.clt_batch_assign)

			  if (correctionsParameters.clt_batch_gen3d) {
				  boolean ok = output3d(
						  clt_parameters,      // EyesisCorrectionParameters.CLTParameters           clt_parameters,
						  colorProcParameters, // EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
						  rgbParameters,       // EyesisCorrectionParameters.RGBParameters             rgbParameters,
						  threadsMax,          // final int        threadsMax,  // maximal number of threads to launch
						  updateStatus,        // final boolean    updateStatus,
//						  !clt_parameters.batch_run, // !batch_dbg,          //   final boolean    batch_mode,
						  debugLevelInner);         // final int        debugLevel)
				  if (!ok) continue;
			  } else continue; // if (correctionsParameters.clt_batch_gen3d)

			  Runtime.getRuntime().gc();
//			  if (debugLevel > -2) System.out.println("Processing set "+(nSet+1)+" (of "+fileIndices.length+") finished at "+
//					  IJ.d2s(0.000000001*(System.nanoTime()-this.startSetTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
			  if (eyesisCorrections.stopRequested.get()>0) {
				  System.out.println("User requested stop");
				  System.out.println("Processing "+(nSet + 1)+" file sets (of "+setNames.size()+") finished at "+
						  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
				  return;
			  }
		  }
		  if (debugLevel > -2) {
			  System.out.println("Processing set "+nSet+" (of "+setNames.size()+") finished at "+
			  IJ.d2s(0.000000001*(System.nanoTime()-this.startSetTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
		  }

		  System.out.println("Processing "+fileIndices.length+" files ("+setNames.size()+" file sets) finished in "+
				  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
	  }

	  public boolean setGpsLla(
			  String source_file)
	  {
		  ImagePlus imp=(new JP46_Reader_camera(false)).open(
				  "", // path,
				  source_file,
				  "",  //arg - not used in JP46 reader
				  true, // un-apply camera color gains
				  null, // new window
				  false); // do not show
		  if (imp.getProperty("LATITUDE") != null){
			  gps_lla = new double[3];
			  for (int i = 0; i < 3; i++) {
				  gps_lla[i] = Double.NaN;
			  }
				if (imp.getProperty("LATITUDE")  != null) gps_lla[0] =Double.parseDouble((String) imp.getProperty("LATITUDE"));
				if (imp.getProperty("LONGITUDE") != null) gps_lla[1] =Double.parseDouble((String) imp.getProperty("LONGITUDE"));
				if (imp.getProperty("ALTITUDE")  != null) gps_lla[2] =Double.parseDouble((String) imp.getProperty("ALTITUDE"));
				return true;
		  }
		  return false;
	  }


	  public boolean writeKml(
			  int debugLevel )
	  {
		  String [] sourceFiles_main=correctionsParameters.getSourcePaths();
		  SetChannels [] set_channels = setChannels(image_name,debugLevel); // only for specified image timestamp

		  ArrayList<String> path_list = new ArrayList<String>();
		  for (int i = 0; i < set_channels.length; i++) {
			  for (int fn:set_channels[i].file_number) {
				  path_list.add(sourceFiles_main[fn]);
			  }
		  }
		  for (String fname:path_list) {
			  if (setGpsLla(fname)) {
				  break;
			  }
		  }
		  if (gps_lla != null) {
			  String kml_copy_dir= correctionsParameters.selectX3dDirectory(
					  image_name, // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
					  null,
					  true,  // smart,
					  true);  //newAllowed, // save
			  double ts = Double.parseDouble(image_name.replace('_', '.'));
			  (new X3dOutput()).generateKML(
					  kml_copy_dir+ Prefs.getFileSeparator()+image_name+".kml", // String path,
					  false, // boolean overwrite,
					  "", // String icon_path, //<href>x3d/1487451413_967079.x3d</href> ?
					  ts, // double timestamp,
					  gps_lla); // double [] lla)
		  } else {
			  if (debugLevel > -2) {
				  System.out.println("GPS data not available, skipping KML file generation (TODO: maybe make some default LLA?)");
			  }
		  }
		  return true;
	  }

	  public boolean createThumbNailImage(
			  ImagePlus imp,
			  String dir,
			  String name,
			  int debugLevel)
	  {
		  String thumb_path = dir +  Prefs.getFileSeparator() + name+".jpeg";
			if (new File(thumb_path).exists()) {
				System.out.println("file "+thumb_path+" exists, skipping thumbnail generation");
				return false;
			}

		  int image_width = imp.getWidth();
		  int image_height = imp.getHeight();
		  double scale_h = 1.0 * (correctionsParameters.thumb_width + 1)/image_width;
		  double scale_v = 1.0 * (correctionsParameters.thumb_height + 1)/image_height;
		  double scale = ((scale_h > scale_v) ? scale_h : scale_v) / correctionsParameters.thumb_size;

		  ImageProcessor ip = imp.getProcessor().duplicate();

		  ip.setInterpolationMethod(ImageProcessor.BICUBIC);
		  ip.blurGaussian(2.0);
		  ip.scale(scale, scale);
		  int lm = (int) Math.round (((image_width*scale)-correctionsParameters.thumb_width)* correctionsParameters.thumb_h_center + (0.5*image_width*(1.0-scale)));
		  int tm = (int) Math.round (((image_height*scale)-correctionsParameters.thumb_height)* correctionsParameters.thumb_v_center + (0.5*image_height*(1.0-scale)));
		  Rectangle r = new Rectangle(lm,tm,correctionsParameters.thumb_width,correctionsParameters.thumb_height);
		  ip.setRoi(r);
		  ImageProcessor ip2 = ip.crop();
		  ImagePlus ip_thumb = new ImagePlus(name,ip2);
		  eyesisCorrections.saveAndShow(
				  ip_thumb,
				  dir,
				  false,
				  false,
				  correctionsParameters.JPEG_quality, // jpegQuality); // jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
				  (debugLevel > -2) ? debugLevel : 1); // int debugLevel (print what it saves)
		  return true;
	  }



	  public boolean writeRatingFile(
			  int debugLevel
			  )
	  {
		  String set_name = image_name;
		  if (set_name == null ) {
			  QuadCLT.SetChannels [] set_channels = setChannels(debugLevel);
			  set_name = set_channels[0].set_name;
		  }

		  String model_dir= correctionsParameters.selectX3dDirectory(
				  set_name, // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
				  null,
				  true,  // smart,
				  true);  //newAllowed, // save
		  String fname = model_dir+ Prefs.getFileSeparator()+"rating.txt";
		  File rating_file = new File(fname);
		  if (rating_file.exists()) {
			  if (debugLevel > -2){
				  System.out.println("file "+rating_file.getPath()+" exists, skipping overwrite");
			  }
			  return false;
		  }
		  List<String> lines = Arrays.asList(correctionsParameters.default_rating+"");
		  Path path = Paths.get(fname);
		  try {
			  Files.write(path,lines, Charset.forName("UTF-8"));
		  } catch (IOException e1) {
			  // TODO Auto-generated catch block
			  e1.printStackTrace();
			  return false;
		  }
		  try {
			  Path fpath = Paths.get(rating_file.getCanonicalPath());
			  Set<PosixFilePermission> perms =  Files.getPosixFilePermissions(fpath);
			  perms.add(PosixFilePermission.OTHERS_WRITE);
			  Files.setPosixFilePermissions(fpath, perms);
		  } catch (IOException e) {
			  // TODO Auto-generated catch block
			  e.printStackTrace();
			  return false;
		  }
		  return true;
	  }




}
