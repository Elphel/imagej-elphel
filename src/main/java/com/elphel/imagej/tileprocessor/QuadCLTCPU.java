package com.elphel.imagej.tileprocessor;
/**
 **
 ** QuadCLTCPU - Process images with CLT-based methods (code specific to ImageJ plugin)
 ** Using CPU
 **
 ** Copyright (C) 2017-2020 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  QuadCLTCPU.java is free software: you can redistribute it and/or modify
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
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Properties;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.calibration.PixelMapping;
import com.elphel.imagej.cameras.CLTParameters;
import com.elphel.imagej.cameras.ColorProcParameters;
import com.elphel.imagej.cameras.EyesisCorrectionParameters;
import com.elphel.imagej.cameras.ThermalColor;
import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.correction.CorrectionColorProc;
import com.elphel.imagej.correction.EyesisCorrections;
import com.elphel.imagej.jp4.JP46_Reader_camera;
import com.elphel.imagej.tileprocessor.GeometryCorrection.CorrVector;
import com.elphel.imagej.x3d.export.WavefrontExport;
import com.elphel.imagej.x3d.export.X3dOutput;

import ij.CompositeImage;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
//import ij.gui.Overlay;
import ij.io.FileSaver;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;


public class QuadCLTCPU {
	public static final String [] FGBG_TITLES_ADJ = {"disparity","strength"};
//	public static final String [] FGBG_TITLES = {"disparity","strength", "rms","rms-split","fg-disp","fg-str","bg-disp","bg-str"};
	public static final String [] FGBG_TITLES_AUX = {"disparity","strength", "rms","rms-split","fg-disp","fg-str","bg-disp","bg-str","aux-disp","aux-str"};
//	public static final enum      FGBG           {DISPARITY,  STRENGTH,   RMS,  RMS_SPLIT,  FG_DISP,  FG_STR,  BG_DISP,  BG_STR};
	public static final int       FGBG_DISPARITY = 0;
	public static final int       FGBG_STRENGTH =  1;
	public static final int       FGBG_RMS =       2;
	public static final int       FGBG_RMS_SPLIT = 3;
	public static final int       FGBG_FG_DISP =   4;
	public static final int       FGBG_FG_STR =    5;
	public static final int       FGBG_BG_DISP =   6;
	public static final int       FGBG_BG_STR =    7;
	public static final int       FGBG_AUX_DISP =  8; // AUX calculated disparity
	public static final int       FGBG_AUX_STR =   9; // AUX calculated strength

//	public GPUTileProcessor.GpuQuad gpuQuad =      null;

	static String []                                       fine_corr_coeff_names = {"A","B","C","D","E","F"};
	static String []                                       fine_corr_dir_names = {"X","Y"};
	public static String                                   PREFIX =     "EYESIS_DCT.";    // change later (first on save)
	public static String                                   PREFIX_AUX = "EYESIS_DCT_AUX."; // change later (first on save)
	static int                                             QUAD =  4; // number of cameras
	public Properties                                      properties = null;
//	public String                                          properties_prefix = "EYESIS_DCT.";
	public EyesisCorrections                               eyesisCorrections = null;
	public EyesisCorrectionParameters.CorrectionParameters correctionsParameters=null;
	double [][][][][][]                                    clt_kernels = null; // can be used to determine monochrome too?
	public GeometryCorrection                              geometryCorrection = null;
	double []                                              extrinsic_vect = new double [GeometryCorrection.CORR_NAMES.length]; // extrinsic corrections (needed from properties, before geometryCorrection
	public int                                             extra_items = 8; // number of extra items saved with kernels (center offset (partial, full, derivatives)
	public ImagePlus                                       eyesisKernelImage = null;
	public long                                            startTime;     // start of batch processing
	public long                                            startSetTime;  // start of set processing
	public long                                            startStepTime; // start of step processing

	public double [][][]                                   fine_corr  = new double [4][2][6]; // per port, per x/y, set of 6 coefficient for fine geometric corrections

	public TileProcessor                                   tp = null;

	public String                                          image_name = null;
	public String                                          image_path = null;
	double []                                              gps_lla =    null;
	public double [][][]                                   image_data = null;
	boolean new_image_data =                               false;
    boolean [][]                                           saturation_imp = null; // (near) saturated pixels or null
    boolean                                                is_aux = false;
    boolean                                                is_mono = false; // Use clt_kernels?
    double  []                                             lwir_offsets = null; // per image subtracted values
    double                                                 lwir_offset =  Double.NaN; // average of lwir_offsets[]
    // hot and cold are calculated during autoranging (when generating 4 images for restored (added lwir_offset)
    // absolute temperatures to be used instead of colorProcParameters lwir_low and lwir_high if autoranging
    // is enabled
    double []                                              lwir_cold_hot = null;
//    int []                                                 woi_tops; // used to calculate scanline timing
// just for debugging with the use of intermediate image
    public double [][]                                     ds_from_main = null;
    public boolean hasNewImageData() {
    	return new_image_data;
    }
    public double [][][]                                   getImageData(){
    	new_image_data = false;
    	return image_data;
    }

// magic scale should be set before using  TileProcessor (calculated disparities depend on it)
    public boolean isMonochrome() {return is_mono;}  // USED in lwir

    public boolean isAux()        {return is_aux;} // USED in lwir
    public String  sAux()         {return isAux()?"-AUX":"";} // USED in lwir
    public boolean isLwir()       {return !Double.isNaN(lwir_offset);} // clt_kernels // USED in lwir
    public double  getLwirOffset() {return lwir_offset;} // USED in lwir

    public double [] getColdHot() { // USED in lwir
    	return lwir_cold_hot;
    }
    public void setColdHot(double [] cold_hot) { // USED in lwir
    	lwir_cold_hot = cold_hot;
    }
    public void setColdHot(double cold, double hot) { // not used in lwir
    	lwir_cold_hot = new double[2];
    	lwir_cold_hot[0] = cold;
    	lwir_cold_hot[1] = hot;
    }

    public void    resetGroundTruthByRig() { // not used in lwir
    	tp.rig_disparity_strength = null;
    }
    public double [][] getGroundTruthByRig(){ // not used in lwir
    	if (tp == null) return null;
    	return tp.rig_disparity_strength;
    }
	public void setTiles (ImagePlus imp, // set tp.tilesX, tp.tilesY // USED in lwir
			CLTParameters    clt_parameters,
			int threadsMax
			){
		setTiles(clt_parameters,
				imp.getWidth()/clt_parameters.transform_size,
				imp.getHeight()/clt_parameters.transform_size,
				threadsMax);
	}

	public void setTiles ( // USED in lwir
			CLTParameters    clt_parameters,
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
					isMonochrome(),
					isLwir(),
					isAux(),
					clt_parameters.corr_magic_scale,
					clt_parameters.grow_disp_trust,
					clt_parameters.max_overexposure, // double maxOverexposure,
					threadsMax);
		}
	}

// used for aux camera only
	public boolean setupImageData( // not used in lwir Called for quadCLT_aux instance
			String image_name,
			String [] sourceFiles,
			CLTParameters       clt_parameters,
			ColorProcParameters                            colorProcParameters, //
			int threadsMax,
			int debugLevel) {
		  QuadCLTCPU.SetChannels [] set_channels_aux =  setChannels(image_name, debugLevel);
		  if ((set_channels_aux == null) || (set_channels_aux.length==0)) {
			  System.out.println("No files for the auxiliary camera match series "+image_name);
			  return false;
		  }
		  double [] referenceExposures_aux =  eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		  int [] channelFiles_aux =  set_channels_aux[0].fileNumber();
		  // make single
		  boolean [][] saturation_aux =  (clt_parameters.sat_level > 0.0)? new boolean[channelFiles_aux.length][] : null;
		  double [] scaleExposures_aux =  new double[channelFiles_aux.length];
//		  ImagePlus [] imp_srcs_aux = 
		  conditionImageSet(
				  clt_parameters,               // EyesisCorrectionParameters.CLTParameters  clt_parameters,
				  colorProcParameters,          //  ColorProcParameters                       colorProcParameters, //
				  sourceFiles,                  // String []                                 sourceFiles,
				  image_name,      // set_channels_aux[0].name(), // String                                    set_name,
				  referenceExposures_aux,       // double []                                 referenceExposures,
				  channelFiles_aux,             // int []                                    channelFiles,
				  scaleExposures_aux,           //output  // double [] scaleExposures
				  saturation_aux,               //output  // boolean [][]                              saturation_imp,
				  threadsMax,                 // int                                       threadsMax,
				  debugLevel); // int                                       debugLevel);

 /*		// 08/12/2020 Moved to condifuinImageSet
		  double [][][] double_stacks_aux = new double [imp_srcs_aux.length][][];
		  for (int i = 0; i < double_stacks_aux.length; i++){
			  double_stacks_aux[i] = eyesisCorrections.bayerToDoubleStack(
					  imp_srcs_aux[i], // source Bayer image, linearized, 32-bit (float))
					  null, // no margins, no oversample
					  this.is_mono);
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
		  new_image_data = true; // ?		  
		  saturation_imp =  saturation_aux;
		  tp.setTrustedCorrelation(clt_parameters.grow_disp_trust);
		  tp.resetCLTPasses();
*/		  
		  return true;

	}


	public QuadCLTCPU(
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
//		System.out.println("new QuadCLTCPU(), prefix = "+prefix);
		getProperties(prefix);
	}

	// TODO:Add saving just calibration

//	public void setProperties(){
//		setProperties(this.properties_prefix);
//	}
	public void setProperties(String prefix, Properties properties){ // save // USED in lwir
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
			gc = new GeometryCorrection(this.extrinsic_vect); // not used in lwir
		}
		for (int i = 0; i < GeometryCorrection.CORR_NAMES.length; i++){
			String name = prefix+"extrinsic_corr_"+GeometryCorrection.CORR_NAMES[i];
			properties.setProperty(name,  gc.getCorrVector().toArray()[i]+"");
		}
		if (is_aux && (gc.rigOffset != null)) {
			gc.rigOffset.setProperties(prefix,properties);
		}
	}


	public void copyPropertiesFrom(Properties other_properties, String other_prefix, String this_prefix){ // save // not used in lwir
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
			gc = new GeometryCorrection(this.extrinsic_vect);
		}
		for (int i = 0; i < GeometryCorrection.CORR_NAMES.length; i++){
			String other_name = other_prefix+"extrinsic_corr_"+GeometryCorrection.CORR_NAMES[i];
  			if (other_properties.getProperty(other_name)!=null) {
  				this.extrinsic_vect[i] = Double.parseDouble(other_properties.getProperty(other_name));
  				if (geometryCorrection != null){
  					geometryCorrection.getCorrVector().toArray()[i] = this.extrinsic_vect[i];
  				}
  			}
			String this_name =  this_prefix+"extrinsic_corr_"+GeometryCorrection.CORR_NAMES[i];
			properties.setProperty(this_name,  gc.getCorrVector().toArray()[i]+"");
//			System.out.println("copyPropertiesFrom():"+i+": setProperty("+this_name+","+gc.getCorrVector().toArray()[i]+"");
		}
//		System.out.println("Done copyPropertiesFrom");
	}

	public GeometryCorrection  getGeometryCorrection() { // USED in lwir
		return geometryCorrection;
	}
	public double [][][][][][] getCLTKernels(){ // USED in lwir
		return clt_kernels;
	}
	public void listGeometryCorrection(boolean full){ // not used in lwir
		GeometryCorrection gc = geometryCorrection;
		if (gc == null) { // if it was not yet created
			gc = new GeometryCorrection(this.extrinsic_vect);
		}
		gc.listGeometryCorrection(full);
	}

	public void getProperties(String prefix){ // restore // USED in lwir
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
  				if (this.extrinsic_vect == null) { // not used in lwir
  					// only create non-null array if there are saved values
  					this.extrinsic_vect = new double [GeometryCorrection.CORR_NAMES.length];
  				}
  				this.extrinsic_vect[i] = Double.parseDouble(properties.getProperty(name));
//  				System.out.println("getProperties():"+i+": getProperty("+name+") -> "+properties.getProperty(name)+"");

  				if (geometryCorrection != null){ // not used in lwir
//  					if (geometryCorrection.getCorrVector().toArray() == null) {
//  						geometryCorrection.resetCorrVector(); // make it array of zeros
//  					}
//  					geometryCorrection.getCorrVector().toArray()[i] = this.extrinsic_vect[i];
  					geometryCorrection.setCorrVector(i,this.extrinsic_vect[i]);
  				}
  			}
		}
//		if (is_aux && (geometryCorrection != null)) {
//			geometryCorrection.setRigOffsetFromProperies(prefix, properties);
//		}
		if (geometryCorrection == null) {
			double [] extrinsic_vect_saved = this.extrinsic_vect.clone();
			boolean OK = initGeometryCorrection(0); // int debugLevel);
			if (!OK) { // not used in lwir
				throw new IllegalArgumentException ("Failed to initialize geometry correction");
			}
			// Substitute vector generated in initGeometryCorrection with the saved from properties one:
			// it also replaces data inside geometryCorrection. TODO: redo to isolate this.extrinsic_vect from geometryCorrection
			this.extrinsic_vect = 	extrinsic_vect_saved;
			geometryCorrection.setCorrVector(this.extrinsic_vect);
//			geometryCorrection = new GeometryCorrection(this.extrinsic_vect);
		}

		if (is_aux) {
			geometryCorrection.setRigOffsetFromProperies(prefix, properties);
		}
	}

	public void setKernelImageFile(ImagePlus img_kernels){ // not used in lwir
		eyesisKernelImage = img_kernels;
	}

	public boolean kernelImageSet(){ // not used in lwir
		return eyesisKernelImage != null;
	}

	public boolean CLTKernelsAvailable(){ // USED in lwir
		return clt_kernels != null;
	}
	public boolean geometryCorrectionAvailable(){ // USED in lwir
		return (geometryCorrection != null) && geometryCorrection.isInitialized();
	}
	public void resetGeometryCorrection() {
		geometryCorrection = null;
//		extrinsic_vect = new double [GeometryCorrection.CORR_NAMES.length];
		extrinsic_vect = null;
	}
	public boolean initGeometryCorrection(int debugLevel){ // USED in lwir
		// keep rig offsets if edited
		if (geometryCorrection == null) {
			geometryCorrection = new GeometryCorrection(extrinsic_vect);
		}
		if (eyesisCorrections.pixelMapping == null) {
			// need to initialize sensor data
//			eyesisCorrections.initSensorFiles(.debugLevel..);
			eyesisCorrections.initPixelMapping(debugLevel);
		}
		PixelMapping.SensorData [] sensors =  eyesisCorrections.pixelMapping.sensors;
		// verify that all sensors have the same distortion parameters
		int numSensors = sensors.length;
		for (int i = 1; i < numSensors; i++){
			if (//	(sensors[0].focalLength !=           sensors[i].focalLength) || // null pointer
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
				return false; // not used in lwir
			}
		}


		// TODO: Verify correction sign!
		double f_avg = geometryCorrection.getCorrVector().setZoomsFromF(
				sensors[0].focalLength,
				sensors[1].focalLength,
				sensors[2].focalLength,
				sensors[3].focalLength);

		// following parameters are used for scaling extrinsic corrections
		geometryCorrection.focalLength = f_avg;
		geometryCorrection.pixelSize = sensors[0].pixelSize;
		geometryCorrection.distortionRadius = sensors[0].distortionRadius;

		for (int i = CorrVector.LENGTH_ANGLES; i < CorrVector.LENGTH; i++){
		}
		// set common distportion parameters
		geometryCorrection.setDistortion(
				f_avg, // sensors[0].focalLength,
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
		/*
		for (int i = 1; i < numSensors; i++){
			if (	(sensors[0].theta !=                 sensors[i].theta) || // elevation
					(sensors[0].heading !=               sensors[i].heading)){
				System.out.println("initGeometryCorrection(): All sensors have to have the same elevation and heading, but channels 0 and "+i+" mismatch");
				return false;
			}
		}
		*/
		double theta_avg = geometryCorrection.getCorrVector().setTiltsFromThetas(
				sensors[0].theta,
				sensors[1].theta,
				sensors[2].theta,
				sensors[3].theta);

		double heading_avg = geometryCorrection.getCorrVector().setAzimuthsFromHeadings(
				sensors[0].heading,
				sensors[1].heading,
				sensors[2].heading,
				sensors[3].heading);


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
				theta_avg, // sensors[0].theta,
				heading_avg, // sensors[0].heading,
				forward,
				right,
				height,
				roll,
				pXY0);
		geometryCorrection.planeProjectLenses(); // project all lenses to the common plane

		// calcualte reverse distortion as a table to be linear interpolated (now cubic!)
		geometryCorrection.calcReverseDistortionTable();

		if (numSensors == 4){
			geometryCorrection.adustSquare();
			System.out.println("Adjusted camera to orient X Y along the sides of a square");
		} else {
			System.out.println("============= Cannot adustSquare() as it requires exactly 4 sensors, "+numSensors+" provided ==========");
			return false; // not used in lwir
		}
		// Print parameters
		if (debugLevel > 0){
			geometryCorrection.listGeometryCorrection(debugLevel > 1);
			System.out.println("=== Extrinsic corrections ===");
			System.out.println(geometryCorrection.getCorrVector().toString());
		}

//listGeometryCorrection
		return true;
	}

//GeometryCorrection

	  public double [][][][][] calculateCLTKernel ( // per color/per tileY/ per tileX/per quadrant (plus offset as 5-th)/per pixel // not used in lwir
			  final PixelMapping.SensorData sensor, // to calculate extra shift
			  final ImageStack kernelStack,  // first stack with 3 colors/slices convolution kernels
			  final int          kernelSize, // 64
			  final CLTParameters clt_parameters,

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
					  ImageDtt image_dtt = new ImageDtt(
							  clt_parameters.transform_size,
							  isMonochrome(),
							  isLwir(),
							  clt_parameters.getScaleStrength(isAux()));
					  int chn,tileY,tileX;
					  DttRad2 dtt = new DttRad2(dtt_size);
					  ShowDoubleFloatArrays sdfa_instance = null;
					  if (globalDebugLevel > -1) sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?

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
								  kernelSize); // ,      // int src_size, // 64
///								  dtt_size);       // 8
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
///									  dtt_size,
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
								  clt_kernels[chn][tileY][tileX]); // , // 	double [][]   sym_kernels, // set of 4 SS, AS, SA, AA kdernels, each dtt_size * dtt_size (may have 5-th with center shift
///								  dtt_size); // 8
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
		  ImageDtt image_dtt = new ImageDtt(
				  clt_parameters.transform_size,
				  isMonochrome(),
				  isLwir(),
				  clt_parameters.getScaleStrength(isAux()));
		  image_dtt.clt_fill_coord_corr(
				  clt_parameters.kernel_step,  //  final int             kern_step, // distance between kernel centers, in pixels.
				  clt_kernels,                 // final double [][][][] clt_data,
				  threadsMax,                  // maximal number of threads to launch
				  globalDebugLevel);
		  return clt_kernels;
	  }

	  public double [][] flattenCLTKernels (      // per color, save 4 kernelas and displacement as (2*dtt_size+1)*(2*dtt_size) tiles in an image (last row - 4 values shift x,y) // not used in lwir
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

	  public void showCLTKernels( // not used in lwir
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

	  public void showCLTKernels( // not used in lwir
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
		  ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays();
		  sdfa_instance.showArrays(
				  flat_kernels,
				  clt_kernels[chn][0][0].length*(2*dtt_size),
				  clt_kernels[chn][0].length*(2*dtt_size+1),
				  true,
				  "clt_kernels-"+chn,
				  titles);
	  }

	  // USED in lwir
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


	  public boolean createCLTKernels( // not used in lwir
			  CLTParameters clt_parameters,
			  int          srcKernelSize,
			  int          threadsMax,  // maximal number of threads to launch
			  boolean      updateStatus,
			  int          debugLevel
			  ){
		  // get sensor geometry correction to apply to kernels as extra shifts
		  PixelMapping.SensorData [] sensors =  eyesisCorrections.pixelMapping.sensors;

		  String [] sharpKernelPaths= correctionsParameters.selectKernelChannelFiles(
				  0,  // 0 - sharp, 1 - smooth
				  correctionsParameters.firstSubCameraConfig,
				  correctionsParameters.numSubCameras,
//				  eyesisCorrections.usedChannels.length, // numChannels, // number of channels
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
		  ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays();

		  for (int chn=0;chn<eyesisCorrections.usedChannels.length;chn++){
			  if (eyesisCorrections.usedChannels[chn] && (sharpKernelPaths[chn]!=null) && (clt_kernels[chn]==null)){
				  ImagePlus imp_kernel_sharp=new ImagePlus(sharpKernelPaths[chn]);
				  if ((imp_kernel_sharp.getStackSize()<3) && (imp_kernel_sharp.getStackSize() != 1)) {
					  System.out.println("Need a 3-layer stack with Bayer or single for mono kernels");
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
				  if (flat_kernels.length ==1){
					  layerNames = new String[1];
					  layerNames[0] = "mono_clt_kernels";
				  }
				  ImageStack cltStack = sdfa_instance.makeStack(
						  flat_kernels,
						  width,
						  height,
						  layerNames);
				  String cltPath=correctionsParameters.cltKernelDirectory+
						  Prefs.getFileSeparator()+
						  correctionsParameters.cltKernelPrefix+
						  String.format("%02d",chn + correctionsParameters.firstSubCameraConfig)+
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
//				  fs.saveAsTiffStack(cltPath); // directory does not exist
				  fs.saveAsTiff(cltPath); // directory does not exist
			  }
		  }
		  return true;
	  }



	  public boolean readCLTKernels( // USED in lwir
			  CLTParameters clt_parameters,
			  int          threadsMax,  // maximal number of threads to launch
			  boolean      updateStatus,
			  int          debugLevel
			  ){
		  int dtt_size = clt_parameters.transform_size;
		  String [] cltKernelPaths = correctionsParameters.selectCLTChannelFiles(
				  correctionsParameters.firstSubCameraConfig,
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
		  ShowDoubleFloatArrays sdfa_instance = null;
		  if (debugLevel>0){
			  sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
		  }

		  for (int chn=0;chn<eyesisCorrections.usedChannels.length;chn++){
			  if (eyesisCorrections.usedChannels[chn] && (cltKernelPaths[chn]!=null)){
				  ImagePlus imp_kernel_clt=new ImagePlus(cltKernelPaths[chn]);
				  if ((imp_kernel_clt.getStackSize()< 3) && (imp_kernel_clt.getStackSize()!= 1))   {
					  System.out.println("Need a 3-layer stack or Bayer and 1-layer for mono with CLT kernels");
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
	  public double [] extractOneKernelFromStack( // not used in lwir
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

	  public double [] extractOneKernelFromStack( // not used in lwir
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
	  private void extractOneKernel(float [] pixels, //  array of combined square kernels, each // not used in lwir
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

	  public double [] reformatKernel( // not used in lwir
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
	  private void reformatKernel( // not used in lwir
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
	  public double []reformatKernel2( // averages by exactly 2 (decimate==2) // not used in lwir
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

	  private void reformatKernel2( // averages by exactly 2 (decimate==2) // not used in lwir
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

	  public void resetCLTKernels() // and geometry correction too // not used in lwir
	  {
		  clt_kernels = null;
		  geometryCorrection=null;

	  }

	  public ImageStack  YPrPbToRGB(double [][] yPrPb,  // USED in lwir
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

	  public double [][]  YPrPbToRBG(double [][] yPrPb, // not used in lwir
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

	  public void debayer_rbg( // not used in lwir
			  ImageStack stack_rbg){
		  debayer_rbg(stack_rbg, 1.0);
	  }

	  // Simple in-place debayer by (bi) linear approximation, assumes [0R/00], [00/B0], [G0/0G] slices
	  public void debayer_rbg( // not used in lwir
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

	  public void processCLTChannelImages( // not used in lwir
			  CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
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
				  int [] channels= fileChannelToSensorChannels(correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile]));
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
				  int [] channels= fileChannelToSensorChannels(correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile]));
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

			  imp_src = eyesisCorrections.getJp4Tiff(sourceFiles[nFile], this.geometryCorrection.woi_tops);

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

	  public ImagePlus processCLTChannelImage( // not used in lwir
			  ImagePlus imp_src, // should have properties "name"(base for saving results), "channel","path"
//			  EyesisCorrectionParameters.DCTParameters           dct_parameters,
			  CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
//			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  ColorProcParameters colorProcParameters,
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
		  ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?

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

//		  String title=name+"-"+String.format("%02d", channel);
		  String title=String.format("%s%s-%02d",name, sAux(), channel);
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
				  null, // no margins, no oversample
				  this.is_mono);

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
			  ImageDtt image_dtt = new ImageDtt(
					  clt_parameters.transform_size,
					  isMonochrome(),
					  isLwir(),
					  clt_parameters.getScaleStrength(isAux()));
			  for (int i =0 ; i < double_stack[0].length; i++){
				  double_stack[2][i]*=0.5; // Scale blue twice to compensate less pixels than green
			  }
			  double [][][][][] clt_data = image_dtt.clt_aberrations(
					  double_stack,                 // final double [][]       imade_data,
					  imp_src.getWidth(),           //	final int               width,
					  clt_kernels[channel],         // final double [][][][][] clt_kernels, // [color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
					  clt_parameters.kernel_step,
//					  clt_parameters.transform_size,
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
				  if (clt_parameters.getCorrSigma(image_dtt.isMonochrome()) > 0){ // no filter at all
					  for (int chn = 0; chn < clt_data.length; chn++) {
						  image_dtt.clt_lpf(
								  clt_parameters.getCorrSigma(image_dtt.isMonochrome()),
								  clt_data[chn],
///								  clt_parameters.transform_size,
								  threadsMax,
								  debugLevel);
					  }
				  }
/*
			  }
*/
			  int tilesY = imp_src.getHeight()/image_dtt.transform_size;
			  int tilesX = imp_src.getWidth()/image_dtt.transform_size;
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
			        			tilesX*image_dtt.transform_size,
			        			tilesY*image_dtt.transform_size,
			        			true,
			        			result.getTitle()+"-CLT");
			        }
			  }
			  double [][] iclt_data = new double [clt_data.length][];
			  for (int chn=0; chn<clt_data.length;chn++){
				  iclt_data[chn] = image_dtt.iclt_2d(
						  clt_data[chn],                  // scanline representation of dcd data, organized as dct_size x dct_size tiles
//						  image_dtt.transform_size,  // final int
						  clt_parameters.clt_window,      // window_type
						  15,                             // clt_parameters.iclt_mask,       //which of 4 to transform back
						  0,                              // clt_parameters.dbg_mode,        //which of 4 to transform back
						  threadsMax,
						  debugLevel);

			  }
					  if (debugLevel > -1) sdfa_instance.showArrays(
							  iclt_data,
							  (tilesX + 1) * image_dtt.transform_size,
							  (tilesY + 1) * image_dtt.transform_size,
							  true,
							  result.getTitle()+"-rbg_sigma");
				/*
				  }
			  }
			 */
			  if (debugLevel > 0) sdfa_instance.showArrays(iclt_data,
					  (tilesX + 1) * image_dtt.transform_size,
					  (tilesY + 1) * image_dtt.transform_size,
					  true,
					  result.getTitle()+"-ICLT-RGB");

			  // convert to ImageStack of 3 slices
			  String [] sliceNames = {"red", "blue", "green"};
			  stack = sdfa_instance.makeStack(
					  iclt_data,
					  (tilesX + 1) * image_dtt.transform_size,
					  (tilesY + 1) * image_dtt.transform_size,
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
	  public void processCLTSets( // not used in lwir
			  CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
//			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  ColorProcParameters colorProcParameters,
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
				  int [] channels= fileChannelToSensorChannels(correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile]));
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
				  int [] channels= fileChannelToSensorChannels(correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile]));
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
//			  setFiles.get(setNames.indexOf(setName)).add(new Integer(nFile));
			  setFiles.get(setNames.indexOf(setName)).add(nFile);
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
					  imp_srcs[srcChannel] = eyesisCorrections.getJp4Tiff(sourceFiles[nFile], this.geometryCorrection.woi_tops);

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

	  public ImagePlus processCLTSetImage( // not used in lwir
			  ImagePlus imp_src, // should have properties "name"(base for saving results), "channel","path"
			  CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
//			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  ColorProcParameters colorProcParameters,
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
		  ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?

		  // may use this.StartTime to report intermediate steps execution times
		  String name=(String) imp_src.getProperty("name");
		  //		int channel= Integer.parseInt((String) imp_src.getProperty("channel"));
		  int channel= (Integer) imp_src.getProperty("channel");
		  String path= (String) imp_src.getProperty("path");

//		  String title=name+"-"+String.format("%02d", channel);
		  String title=String.format("%s%s-%02d",name, sAux(), channel);
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
				  null, // no margins, no oversample
				  this.is_mono);

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
			  ImageDtt image_dtt = new ImageDtt(
					  clt_parameters.transform_size,
					  isMonochrome(),
					  isLwir(),
					  clt_parameters.getScaleStrength(isAux()));
			  for (int i =0 ; i < double_stack[0].length; i++){
				  double_stack[2][i]*=0.5; // Scale blue twice to compensate less pixels than green
			  }
			  double [][][][][] clt_data = image_dtt.clt_aberrations(
					  double_stack,                 // final double [][]       imade_data,
					  imp_src.getWidth(),           //	final int               width,
					  clt_kernels[channel],         // final double [][][][][] clt_kernels, // [color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
					  clt_parameters.kernel_step,
//					  image_dtt.transform_size,
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
				  if (clt_parameters.getCorrSigma(image_dtt.isMonochrome()) > 0){ // no filter at all
					  for (int chn = 0; chn < clt_data.length; chn++) {
						  image_dtt.clt_lpf(
								  clt_parameters.getCorrSigma(image_dtt.isMonochrome()),
								  clt_data[chn],
///								  image_dtt.transform_size,
								  threadsMax,
								  debugLevel);
					  }
				  }
/*
			  }
*/
			  int tilesY = imp_src.getHeight()/image_dtt.transform_size;
			  int tilesX = imp_src.getWidth()/image_dtt.transform_size;
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
			        			tilesX*image_dtt.transform_size,
			        			tilesY*image_dtt.transform_size,
			        			true,
			        			result.getTitle()+"-CLT");
			        }
			  }
			  double [][] iclt_data = new double [clt_data.length][];
			  for (int chn=0; chn<clt_data.length;chn++){
				  iclt_data[chn] = image_dtt.iclt_2d(
						  clt_data[chn],                  // scanline representation of dcd data, organized as dct_size x dct_size tiles
///						  image_dtt.transform_size,  // final int
						  clt_parameters.clt_window,      // window_type
						  15,                             // clt_parameters.iclt_mask,       //which of 4 to transform back
						  0,                              // clt_parameters.dbg_mode,        //which of 4 to transform back
						  threadsMax,
						  debugLevel);

			  }

//					  if (debugLevel > -1) System.out.println("Applyed LPF, sigma = "+dct_parameters.dbg_sigma);
					  if (debugLevel > 0) sdfa_instance.showArrays(
							  iclt_data,
							  (tilesX + 1) * image_dtt.transform_size,
							  (tilesY + 1) * image_dtt.transform_size,
							  true,
							  result.getTitle()+"-rbg_sigma");
				/*
				  }
			  }
			 */
			  if (debugLevel > 0) sdfa_instance.showArrays(iclt_data,
					  (tilesX + 0) * image_dtt.transform_size,
					  (tilesY + 0) * image_dtt.transform_size,
					  true,
					  result.getTitle()+"-ICLT-RGB");

			  // convert to ImageStack of 3 slices
			  String [] sliceNames = {"red", "blue", "green"};
			  stack = sdfa_instance.makeStack(
					  iclt_data,
					  (tilesX + 0) * image_dtt.transform_size,
					  (tilesY + 0) * image_dtt.transform_size,
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

	  public void processCLTQuads( // not used in lwir
			  CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
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
				  int [] channels= fileChannelToSensorChannels(correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile]));
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
				  int [] channels= fileChannelToSensorChannels(correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile]));
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
			  setFiles.get(setNames.indexOf(setName)).add(nFile); // add(new Integer(nFile))
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
			  this.geometryCorrection.woi_tops = new int [channelFiles.length];
			  double [] scaleExposures = new double[channelFiles.length];
			  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
				  int nFile=channelFiles[srcChannel];
				  imp_srcs[srcChannel]=null;
				  if (nFile >=0){
					  imp_srcs[srcChannel] = eyesisCorrections.getJp4Tiff(sourceFiles[nFile], this.geometryCorrection.woi_tops);

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

	  public ImagePlus [] processCLTQuad( // not used in lwir
			  ImagePlus [] imp_quad, // should have properties "name"(base for saving results), "channel","path"
			  CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
//			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  ColorProcParameters colorProcParameters,
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
		  ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?

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
					  null, // no margins, no oversample
					  this.is_mono);
		  }

//		  String [] rbg_titles = {"Red", "Blue", "Green"};
		  ImageStack stack;
		  // =================
		  ImageDtt image_dtt = new ImageDtt(
				  clt_parameters.transform_size,
				  isMonochrome(),
				  isLwir(),
				  clt_parameters.getScaleStrength(isAux()));
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
//				  image_dtt.transform_size,
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

//			  String title=name+"-"+String.format("%02d", iQuad);
			  String title=String.format("%s%s-%02d",name, sAux(), iQuad);
			  String titleFull=title+"-SPLIT";

			  if (clt_parameters.getCorrSigma(image_dtt.isMonochrome()) > 0){ // no filter at all
				  for (int chn = 0; chn < clt_data[iQuad].length; chn++) {
					  image_dtt.clt_lpf(
							  clt_parameters.getCorrSigma(image_dtt.isMonochrome()),
							  clt_data[iQuad][chn],
///							  image_dtt.transform_size,
							  threadsMax,
							  debugLevel);
				  }
			  }

			  int tilesY = imp_quad[iQuad].getHeight()/image_dtt.transform_size;
			  int tilesX = imp_quad[iQuad].getWidth()/image_dtt.transform_size;
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
							  tilesX*image_dtt.transform_size,
							  tilesY*image_dtt.transform_size,
							  true,
							  results[iQuad].getTitle()+"-CLT");
				  }
			  }
			  double [][] iclt_data = new double [clt_data[iQuad].length][];
			  for (int chn=0; chn<iclt_data.length;chn++){
				  iclt_data[chn] = image_dtt.iclt_2d(
						  clt_data[iQuad][chn],           // scanline representation of dcd data, organized as dct_size x dct_size tiles
///						  image_dtt.transform_size,  // final int
						  clt_parameters.clt_window,      // window_type
						  15,                             // clt_parameters.iclt_mask,       //which of 4 to transform back
						  0,                              // clt_parameters.dbg_mode,        //which of 4 to transform back
						  threadsMax,
						  debugLevel);

			  }
			  if (debugLevel > 0) sdfa_instance.showArrays(
					  iclt_data,
					  (tilesX + 0) * image_dtt.transform_size,
					  (tilesY + 0) * image_dtt.transform_size,
					  true,
					  results[iQuad].getTitle()+"-rbg_sigma");
			  if (debugLevel > 0) sdfa_instance.showArrays(iclt_data,
					  (tilesX + 0) * image_dtt.transform_size,
					  (tilesY + 0) * image_dtt.transform_size,
					  true,
					  results[iQuad].getTitle()+"-ICLT-RGB");

			  // convert to ImageStack of 3 slices
			  String [] sliceNames = {"red", "blue", "green"};
			  stack = sdfa_instance.makeStack(
					  iclt_data,
					  (tilesX + 0) * image_dtt.transform_size,
					  (tilesY + 0) * image_dtt.transform_size,
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

	  class SetChannels{ // USED in lwir
		  String set_name;    // set name (timestamp)
		  int [] file_number; // array of file numbers for channels
		  public SetChannels(String name, int[] fn){ // USED in lwir
			  set_name = name;
			  file_number = fn;
		  }
		  public String name() { // USED in lwir
			  return set_name;
		  }
		  public int [] fileNumber() { // USED in lwir
			  return file_number;
		  }
		  public int fileNumber(int i) { // not used in lwir
			  return file_number[i];
		  }
	  }

	  SetChannels [] setChannels( // USED in lwir
			  int debugLevel) {
		  return setChannels(null, debugLevel);
	  }

	  public int [] fileChannelToSensorChannels(int file_channel) { // USED in lwir
		  if (!eyesisCorrections.pixelMapping.subcamerasUsed()) { // not an Eyesis-type system
			  // Here use firstSubCameraConfig - subcamera, corresponding to sensors[0] of this PixelMapping instance (1 for Eyesis, 2 for Rig/LWIR)
			  return eyesisCorrections.pixelMapping.channelsForSubCamera(file_channel - correctionsParameters.firstSubCameraConfig);
		  } else 	if (correctionsParameters.isJP4()){ // not used in lwir
			  // Here use firstSubCamera - first filename index to be processed by this PixelMapping instance (1 for Eyesis, 2 for Rig/LWIR)
			  int subCamera= file_channel- correctionsParameters.firstSubCamera; // to match those in the sensor files
			  return  eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera);
		  } else {
			  int [] channels = {file_channel};
			  return channels;
		  }
	  }

	  SetChannels [] setChannels( // USED in lwir
			  String single_set_name, // process only files that contain specified series (timestamp) in the name
			  int debugLevel) {
		  String [] sourceFiles=correctionsParameters.getSourcePaths();
		  boolean [] enabledFiles=new boolean[sourceFiles.length];
		  for (int i=0;i<enabledFiles.length;i++) enabledFiles[i]=false;
		  int numFilesToProcess=0;
		  int numImagesToProcess=0;
		  for (int nFile=0;nFile<enabledFiles.length;nFile++){
			  if (    (sourceFiles[nFile]!=null) &&
					  (sourceFiles[nFile].length() > 1) &&
					  ((single_set_name == null) || (correctionsParameters.getNameFromTiff(sourceFiles[nFile]).contains(single_set_name)))) {
				  int [] channels= fileChannelToSensorChannels(correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile]));
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
			  return null; // not used in lwir
		  } else {
			  if (debugLevel>0) System.out.println(numFilesToProcess+ " files to process (of "+sourceFiles.length+"), "+numImagesToProcess+" images to process");
		  }
		  int [][] fileIndices=new int [numImagesToProcess][2]; // file index, channel number
		  int index=0;
		  for (int nFile=0;nFile<enabledFiles.length;nFile++){ // enabledFiles not used anymore?
			  if (    (sourceFiles[nFile]!=null) &&
					  (sourceFiles[nFile].length()>1) &&
					  ((single_set_name == null) || (correctionsParameters.getNameFromTiff(sourceFiles[nFile]).contains(single_set_name)))) { // not used in lwir
				  int [] channels= fileChannelToSensorChannels(correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile]));
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
			  String setName = correctionsParameters.getNameFromSourceTiff(sourceFiles[nFile]); // supports set directory name
			  if (!setNames.contains(setName)) {
				  setNames.add(setName);
				  setFiles.add(new ArrayList<Integer>());
			  }
			  setFiles.get(setNames.indexOf(setName)).add(iImage); // .add(new Integer(iImage));
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


	  int getTotalFiles(SetChannels [] sc) { // USED in lwir
		  int nf = 0;
		  for (int i = 0; i < sc.length; i++) nf+=sc[i].fileNumber().length;
		  return nf;
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
	   * @param threadsMax          maximal number of threads to use
	   * @param debugLevel          debug (verbosity) level
	   * @return array of per-channel ImagePlus objects to process (with saturation_imp)
	   */
	  public ImagePlus[] conditionImageSet( // USED in lwir
			  CLTParameters  clt_parameters,
			  ColorProcParameters                       colorProcParameters, //
			  String []                                 sourceFiles,
			  String                                    set_name,
			  double []                                 referenceExposures,
			  int []                                    channelFiles,
			  double []                                 scaleExposures,
			  boolean [][]                              saturation_imp,
			  int                                       threadsMax,
			  int                                       debugLevel)
	  {
		  this.image_name = set_name;
		  ImagePlus [] imp_srcs = new ImagePlus[channelFiles.length];
		  this.geometryCorrection.woi_tops = new int [channelFiles.length];
		  double [][] dbg_dpixels = new double [channelFiles.length][];
		  boolean is_lwir =            colorProcParameters.lwir_islwir;
		  boolean ignore_saturation =  is_lwir;
		  boolean lwir_subtract_dc =   colorProcParameters.lwir_subtract_dc;
		  boolean lwir_eq_chn =        colorProcParameters.lwir_eq_chn;
		  boolean correct_vignetting = colorProcParameters.correct_vignetting;
		  this.is_mono = is_lwir; // maybe add other monochrome?

		  for (int srcChannel=0; srcChannel < channelFiles.length; srcChannel++){
			  int nFile=channelFiles[srcChannel]; // channelFiles[srcChannel];
			  imp_srcs[srcChannel]=null;
			  if (nFile >=0){
				  imp_srcs[srcChannel] = eyesisCorrections.getJp4Tiff(sourceFiles[nFile], this.geometryCorrection.woi_tops);

				  scaleExposures[srcChannel] = 1.0;
				  if (!(referenceExposures == null) && !Double.isNaN(referenceExposures[nFile]) && (imp_srcs[srcChannel].getProperty("EXPOSURE")!=null)){
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
					  int numApplied=	eyesisCorrections.correctDefects( // not used in lwir
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
				  if ((debugLevel > -1) && (!isMonochrome())) {
					  double [] max_pix= {0.0, 0.0, 0.0, 0.0};
//					  for (int y = 0; y < height-1; y+=2){
					  for (int y = 0; (y < 499) && (y < height); y+=2){
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
					  saturation_imp[srcChannel] = new boolean[width*height];
					  if (!ignore_saturation) {
						  double [] saturations = {
								  Double.parseDouble((String) imp_srcs[srcChannel].getProperty("saturation_1")),
								  Double.parseDouble((String) imp_srcs[srcChannel].getProperty("saturation_0")),
								  Double.parseDouble((String) imp_srcs[srcChannel].getProperty("saturation_3")),
								  Double.parseDouble((String) imp_srcs[srcChannel].getProperty("saturation_2"))};
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
				  }

				  if (!is_lwir) { // no vigneting correction and no color scaling
					  if (this.correctionsParameters.vignetting && correct_vignetting){
						  if ((eyesisCorrections.channelVignettingCorrection==null) || (srcChannel<0) || (srcChannel>=eyesisCorrections.channelVignettingCorrection.length) || (eyesisCorrections.channelVignettingCorrection[srcChannel]==null)){
							  System.out.println("No vignetting data for channel "+srcChannel);
							  return null; // not used in lwir
						  }
						  ///						  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();

						  float [] vign_pixels = eyesisCorrections.channelVignettingCorrection[srcChannel];
						  if (pixels.length!=vign_pixels.length){
//							  System.out.println("Vignetting data for channel "+srcChannel+" has "+vign_pixels.length+" pixels, image "+sourceFiles[nFile]+" has "+pixels.length);
							  int woi_width =  Integer.parseInt((String) imp_srcs[srcChannel].getProperty("WOI_WIDTH"));
							  int woi_height = Integer.parseInt((String) imp_srcs[srcChannel].getProperty("WOI_HEIGHT"));
							  int woi_top =  Integer.parseInt((String) imp_srcs[srcChannel].getProperty("WOI_TOP"));
							  int woi_left =  Integer.parseInt((String) imp_srcs[srcChannel].getProperty("WOI_LEFT"));
							  int vign_width =  eyesisCorrections.pixelMapping.sensors[srcChannel].pixelCorrectionWidth;
							  int vign_height = eyesisCorrections.pixelMapping.sensors[srcChannel].pixelCorrectionHeight;

							  if (pixels.length != woi_width * woi_height){
								  System.out.println("Vignetting data for channel "+srcChannel+" has "+vign_pixels.length+" pixels, < "+
										  sourceFiles[nFile]+" has "+pixels.length);
								  woi_width = width;
								  woi_height = height;
							  }
							  if (vign_width < (woi_left + woi_width)) {
								  System.out.println("Vignetting data for channel "+srcChannel+
										  " has width + left ("+(woi_left+woi_width)+") > vign_width ("+vign_width+")");
								  return null;
							  }
							  if (vign_height < (woi_top + woi_height)) {
								  System.out.println("Vignetting data for channel "+srcChannel+
										  " has height + top ("+(woi_top+woi_height)+") > vign_height ("+vign_width+")");
								  return null;
							  }
							  if (pixels.length != woi_width * woi_height){
								  System.out.println("Vignetting data for channel "+srcChannel+" has "+vign_pixels.length+" pixels, < "+
										  sourceFiles[nFile]+" has "+pixels.length);
								  return null;
							  }
							  vign_pixels = new float[woi_width * woi_height];
							  for (int row = 0; row < woi_height; row++) {
								  System.arraycopy(
										  eyesisCorrections.channelVignettingCorrection[srcChannel], // src
										  (woi_top + row) * vign_width + woi_left, // srcPos,
										  vign_pixels,                             // dest,
										  row * woi_width,                         // destPos,
										  woi_width);                              // length);
							  }

						  }
						  // TODO: Move to do it once:
						  double min_non_zero = 0.0;
						  for (int i=0;i<pixels.length;i++){
							  double d = vign_pixels[i];
							  if ((d > 0.0) && ((min_non_zero == 0) || (min_non_zero > d))){
								  min_non_zero = d;
							  }
						  }
						  double max_vign_corr = clt_parameters.vignetting_range*min_non_zero;

						  System.out.println("Vignetting data: channel="+srcChannel+", min = "+min_non_zero);
						  for (int i=0;i<pixels.length;i++){
							  double d = vign_pixels[i];
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

					  } else { // assuming GR/BG pattern // not used in lwir
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
		  }
		  // temporary applying scaleExposures[srcChannel] here, setting it to all 1.0
		  System.out.println("Temporarily applying scaleExposures[] here - 1" );
		  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
			  if (!is_lwir) {
				  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
				  for (int i = 0; i < pixels.length; i++){
					  pixels[i] *= scaleExposures[srcChannel];
				  }
			  }
			  scaleExposures[srcChannel] = 1.0;

		  }

		  if ((debugLevel > -1) && (saturation_imp != null) && !is_lwir){
			  String [] titles = {"chn0","chn1","chn2","chn3"};
			  double [][] dbg_satur = new double [saturation_imp.length] [saturation_imp[0].length];
			  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
				  for (int i = 0; i < saturation_imp[srcChannel].length; i++){
					  dbg_satur[srcChannel][i] = saturation_imp[srcChannel][i]? 1.0 : 0.0;
				  }
			  }
			  int width =  imp_srcs[0].getWidth();
			  int height = imp_srcs[0].getHeight();
			  (new ShowDoubleFloatArrays()).showArrays(dbg_satur, width, height, true, "Saturated" , titles);

			  if ((debugLevel > -1) && !isMonochrome()) { // 0){
				  double [][] dbg_dpixels_norm = new double [channelFiles.length][];
				  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
					  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
					  dbg_dpixels_norm[srcChannel] = new double[pixels.length];
					  for (int i = 0; i < pixels.length; i++){
						  dbg_dpixels_norm[srcChannel][i] = pixels[i];
					  }
				  }
				  (new ShowDoubleFloatArrays()).showArrays(dbg_dpixels, width, height, true, "dpixels" , titles);
				  (new ShowDoubleFloatArrays()).showArrays(dbg_dpixels_norm, width, height, true, "dpixels_norm" , titles);
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
				  (new ShowDoubleFloatArrays()).showArrays(dbg_dpixels_split, width/2, height/2, true, "dpixels_split" , dbg_titles);
			  }

		  }
		  //			  Overlay ovl = imp_srcs[0].getOverlay();
		  // once per quad here
		  // may need to equalize gains between channels
		  if (!is_lwir && (clt_parameters.gain_equalize || clt_parameters.colors_equalize)){ // false, true
			  channelGainsEqualize(
					  clt_parameters.gain_equalize, //false
					  clt_parameters.colors_equalize, // true
					  clt_parameters.nosat_equalize, // boolean nosat_equalize, // true
					  channelFiles,
					  imp_srcs,
					  saturation_imp, // boolean[][] saturated,
					  set_name,       // setNames.get(nSet), // just for debug messages == setNames.get(nSet)
					  debugLevel);
		  }
		  if (is_lwir && (lwir_subtract_dc || lwir_eq_chn)) {
			   this.lwir_offsets = channelLwirEqualize(
						  channelFiles,
						  imp_srcs,
						  lwir_subtract_dc, // boolean      remove_dc,
						  set_name, // just for debug messages == setNames.get(nSet)
						  debugLevel);
			   int num_avg = 0;
			   this.lwir_offset = 0.0;
			   for (int srcChannel=0; srcChannel < channelFiles.length; srcChannel++){
				   int nFile=channelFiles[srcChannel];
				   if (nFile >=0){
					   this.lwir_offset += this.lwir_offsets[srcChannel];
					   num_avg++;
				   }
			   }
			   this.lwir_offset /= num_avg;
		  }
// 08/12/2020 common part moved here, from getRigImageStacks()
		  image_name = (String) imp_srcs[0].getProperty("name");
		  image_path=  (String) imp_srcs[0].getProperty("path");		  
		  this.saturation_imp = saturation_imp;
		  image_data =          new double [imp_srcs.length][][];
		  this.new_image_data = true;		  
		  for (int i = 0; i < image_data.length; i++){
			  image_data[i] = eyesisCorrections.bayerToDoubleStack(
					  imp_srcs[i], // source Bayer image, linearized, 32-bit (float))
					  null, // no margins, no oversample
					  is_mono);
			  
// TODO: Scale greens here ?
			  if (!is_mono && (image_data[i].length > 2)) {
				  for (int j =0 ; j < image_data[i][0].length; j++){
					  image_data[i][2][j]*=0.5; // Scale green 0.5 to compensate more pixels than R,B
				  }
			  }
		  }
		  setTiles (imp_srcs[0], // set global tp.tilesX, tp.tilesY
					clt_parameters,
					threadsMax); // where to get it? Use instance member
		  tp.setTrustedCorrelation(clt_parameters.grow_disp_trust);
		  tp.resetCLTPasses();
		  return imp_srcs;
	  }


	  public void processCLTQuadCorrs( // not used in lwir
			  CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  ColorProcParameters colorProcParameters,
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
		// multiply each image by this and divide by individual (if not NaN)
		  double [] referenceExposures = null;
		  if (!colorProcParameters.lwir_islwir) {
			  referenceExposures=eyesisCorrections.calcReferenceExposures(debugLevel);
		  }
		  for (int nSet = 0; nSet < set_channels.length; nSet++){
			  int [] channelFiles = set_channels[nSet].fileNumber();
			  boolean [][] saturation_imp = (clt_parameters.sat_level > 0.0)? new boolean[channelFiles.length][] : null;
			  double [] scaleExposures = new double[channelFiles.length];

			  ImagePlus [] imp_srcs = conditionImageSet(
					  clt_parameters,             // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					  colorProcParameters,
					  sourceFiles,                // String []                                 sourceFiles,
					  set_channels[nSet].name(),  // String                                    set_name,
					  referenceExposures,         // double []                                 referenceExposures,
					  channelFiles,               // int []                                    channelFiles,
					  scaleExposures,   //output  // double [] scaleExposures
					  saturation_imp,   //output  // boolean [][]                              saturation_imp,
					  threadsMax,                 // int                                       threadsMax,
					  debugLevel); // int                                       debugLevel);


			  // once per quad here
			  processCLTQuadCorrCPU( // returns ImagePlus, but it already should be saved/shown
					  imp_srcs, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
					  saturation_imp, // boolean [][] saturation_imp, // (near) saturated pixels or null
					  clt_parameters,
					  debayerParameters,
					  colorProcParameters,
					  channelGainParameters,
					  rgbParameters,
					  scaleExposures,
					  apply_corr, // calculate and apply additional fine geometry correction
					  infinity_corr, // calculate and apply geometry correction at infinity
					  threadsMax,  // maximal number of threads to launch
					  updateStatus,
					  debugLevel);
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



	  public void processCLTQuadCorrsTestERS(
			  CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  ColorProcParameters colorProcParameters,
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
		// multiply each image by this and divide by individual (if not NaN)
		  double [] referenceExposures = null;
		  if (!colorProcParameters.lwir_islwir) {
			  referenceExposures=eyesisCorrections.calcReferenceExposures(debugLevel);
		  }
		  for (int nSet = 0; nSet < set_channels.length; nSet++){
			  int [] channelFiles = set_channels[nSet].fileNumber();
			  boolean [][] saturation_imp = (clt_parameters.sat_level > 0.0)? new boolean[channelFiles.length][] : null;
			  double [] scaleExposures = new double[channelFiles.length];

			  ImagePlus [] imp_srcs = conditionImageSet(
					  clt_parameters,             // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					  colorProcParameters,
					  sourceFiles,                // String []                                 sourceFiles,
					  set_channels[nSet].name(),  // String                                    set_name,
					  referenceExposures,         // double []                                 referenceExposures,
					  channelFiles,               // int []                                    channelFiles,
					  scaleExposures,   //output  // double [] scaleExposures
					  saturation_imp,   //output  // boolean [][]                              saturation_imp,
					  threadsMax,                 // int                                       threadsMax,
					  debugLevel); // int                                       debugLevel);


			  // once per quad here
			  processCLTQuadCorrTestERS( // returns ImagePlus, but it already should be saved/shown
					  imp_srcs, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
					  saturation_imp, // boolean [][] saturation_imp, // (near) saturated pixels or null
					  clt_parameters,
					  debayerParameters,
					  colorProcParameters,
					  channelGainParameters,
					  rgbParameters,
					  scaleExposures,
					  apply_corr, // calculate and apply additional fine geometry correction
					  infinity_corr, // calculate and apply geometry correction at infinity
					  threadsMax,  // maximal number of threads to launch
					  updateStatus,
					  debugLevel);
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


	  public void processCLTQuadCorrsTest(
			  CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  ColorProcParameters colorProcParameters,
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
		// multiply each image by this and divide by individual (if not NaN)
		  double [] referenceExposures = null;
		  if (!colorProcParameters.lwir_islwir) {
			  referenceExposures=eyesisCorrections.calcReferenceExposures(debugLevel);
		  }
		  for (int nSet = 0; nSet < set_channels.length; nSet++){
			  int [] channelFiles = set_channels[nSet].fileNumber();
			  boolean [][] saturation_imp = (clt_parameters.sat_level > 0.0)? new boolean[channelFiles.length][] : null;
			  double [] scaleExposures = new double[channelFiles.length];

			  ImagePlus [] imp_srcs = conditionImageSet(
					  clt_parameters,             // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					  colorProcParameters,
					  sourceFiles,                // String []                                 sourceFiles,
					  set_channels[nSet].name(),  // String                                    set_name,
					  referenceExposures,         // double []                                 referenceExposures,
					  channelFiles,               // int []                                    channelFiles,
					  scaleExposures,   //output  // double [] scaleExposures
					  saturation_imp,   //output  // boolean [][]                              saturation_imp,
					  threadsMax,                 // int                                       threadsMax,
					  debugLevel); // int                                       debugLevel);


			  // once per quad here
			  processCLTQuadCorrTest( // returns ImagePlus, but it already should be saved/shown
					  imp_srcs, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
					  saturation_imp, // boolean [][] saturation_imp, // (near) saturated pixels or null
					  clt_parameters,
					  debayerParameters,
					  colorProcParameters,
					  channelGainParameters,
					  rgbParameters,
					  scaleExposures,
					  apply_corr, // calculate and apply additional fine geometry correction
					  infinity_corr, // calculate and apply geometry correction at infinity
					  threadsMax,  // maximal number of threads to launch
					  updateStatus,
					  debugLevel);
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



	  public void channelGainsEqualize( // USED in lwir
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
		  } else { // not used in lwir
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
	  public void channelGainsEqualize_old( // not used in lwir
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

	  public void channelGainsEqualize_new( // USED in lwir
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
					  if (gain_equalize){ // not used in lwir
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
	  public double []  channelLwirEqualize( // USED in lwir
			  int [] channelFiles,
			  ImagePlus [] imp_srcs,
			  boolean      remove_dc,
			  String setName, // just for debug messages == setNames.get(nSet)
			  int debugLevel){
		  double [] offsets = new double [channelFiles.length];
		  double [][] avr_pix = new double [channelFiles.length][2]; // val/weight
		  double [] wnd_x = {};
		  double [] wnd_y = {};
		  double total_s = 0.0, total_w = 0.0;
		  for (int srcChannel=0; srcChannel < channelFiles.length; srcChannel++){
			  int nFile=channelFiles[srcChannel];
			  if (nFile >=0){
				  avr_pix[srcChannel][0] = 0.0;
				  avr_pix[srcChannel][1] = 0.0;
				  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
				  int width =  imp_srcs[srcChannel].getWidth();
				  int height = imp_srcs[srcChannel].getHeight();
				  if (wnd_x.length != width) {
					  wnd_x = new double[width];
					  for (int i = 0; i < width; i++) {
						  wnd_x[i] = 0.5 - 0.5*Math.cos(2*Math.PI * (i+1) / (width + 1));
					  }
				  }
				  if (wnd_y.length != height) {
					  wnd_y = new double[height];
					  for (int i = 0; i < height; i++) {
						  wnd_y[i] = 0.5 - 0.5*Math.cos(2*Math.PI * (i+1) / (height + 1));
					  }
				  }
				  int indx = 0;
				  for (int y = 0; y < height; y++) {
					  for (int x = 0; x < width; x++) {
						  double w = wnd_y[y]*wnd_x[x];
						  avr_pix[srcChannel][0] += w * pixels[indx++];
						  avr_pix[srcChannel][1] += w;
					  }
				  }
				  total_s += avr_pix[srcChannel][0];
				  total_w += avr_pix[srcChannel][1];
				  avr_pix[srcChannel][0]/=avr_pix[srcChannel][1]; // weighted average
			  }
		  }
		  double avg = total_s/total_w;
		  if (!remove_dc) { // not used in lwir
			  for (int srcChannel=0; srcChannel < channelFiles.length; srcChannel++) if (channelFiles[srcChannel] >=0){
				  avr_pix[srcChannel][0] -= avg;
			  }

		  }
		  for (int srcChannel=0; srcChannel < channelFiles.length; srcChannel++){
			  int nFile=channelFiles[srcChannel];
			  if (nFile >=0) {
//				  offsets[srcChannel]= (avr_pix[srcChannel][0] - (remove_dc ? 0.0: avg));
				  offsets[srcChannel]= avr_pix[srcChannel][0];
				  float fd = (float)offsets[srcChannel];
				  float [] pixels = (float []) imp_srcs[srcChannel].getProcessor().getPixels();
				  for (int i = 0; i < pixels.length; i++) {
					  pixels[i] -= fd;
				  }
			  }
		  }

		  return offsets;
	  }

	  public ImagePlus [] processCLTQuadCorrCPU( // USED in lwir
			  ImagePlus [] imp_quad, // should have properties "name"(base for saving results), "channel","path"
			  boolean [][] saturation_imp, // (near) saturated pixels or null // Not needed use this.saturation_imp
			  CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  ColorProcParameters                              colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters         channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  double []	       scaleExposures, // probably not needed here
			  final boolean    apply_corr, // calculate and apply additional fine geometry correction
			  final boolean    infinity_corr, // calculate and apply geometry correction at infinity
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel){
		  final boolean      batch_mode = clt_parameters.batch_run; //disable any debug images
		  boolean advanced=  this.correctionsParameters.zcorrect || this.correctionsParameters.equirectangular;
		  boolean toRGB=     advanced? true: this.correctionsParameters.toRGB;
		  ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?

		  // may use this.StartTime to report intermediate steps execution times
///		  String name=this.correctionsParameters.getModelName((String) imp_quad[0].getProperty("name"));
///		  String path= (String) imp_quad[0].getProperty("path");

		  ImagePlus [] results = new ImagePlus[imp_quad.length];
		  for (int i = 0; i < results.length; i++) {
			  results[i] = imp_quad[i];
			  results[i].setTitle(results[i].getTitle()+"RAW");
		  }
		  if (debugLevel>1) System.out.println("processing: "+image_path);
/*		  // Moved to condifuinImageSet
		  double [][][] double_stacks = new double [imp_quad.length][][];
		  for (int i = 0; i < double_stacks.length; i++){
			  double_stacks[i] = eyesisCorrections.bayerToDoubleStack(
					  imp_quad[i], // source Bayer image, linearized, 32-bit (float))
					  null, // no margins, no oversample
					  this.is_mono);
		  }

		  for (int i = 0; i < double_stacks.length; i++){
			  if ( double_stacks[i].length > 2) {
				  for (int j =0 ; j < double_stacks[i][0].length; j++){
					  double_stacks[i][2][j]*=0.5; // Scale green 0.5 to compensate more pixels than R,B
				  }
			  } else {
				  for (int j =0 ; j < double_stacks[i][0].length; j++){
					  double_stacks[i][0][j]*=1.0; // Scale mono by 1/4 - to have the same overall "gain" as for bayer
				  }
			  }
		  }

		  setTiles (imp_quad[0], // set global tp.tilesX, tp.tilesY
				  clt_parameters,
				  threadsMax);
*/
		  ImageDtt image_dtt = new ImageDtt(
				  clt_parameters.transform_size,
				  isMonochrome(),
				  isLwir(),
				  clt_parameters.getScaleStrength(isAux()));
		  

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

		  if (clt_parameters.correlate){ // true
			  clt_corr_combo =    new double [ImageDtt.TCORR_TITLES.length][tilesY][tilesX][]; // will it work???
			  texture_tiles =     new double [tilesY][tilesX][][]; // ["RGBA".length()][];
			  for (int i = 0; i < tilesY; i++){
				  for (int j = 0; j < tilesX; j++){
					  for (int k = 0; k<clt_corr_combo.length; k++){
						  clt_corr_combo[k][i][j] = null;
					  }
					  texture_tiles[i][j] = null;
				  }
			  }
			  if (!infinity_corr && clt_parameters.corr_keep){ // true
				  clt_corr_partial = new double [tilesY][tilesX][][][];
				  for (int i = 0; i < tilesY; i++){
					  for (int j = 0; j < tilesX; j++){
						  clt_corr_partial[i][j] = null;
					  }
				  }
			  } // clt_parameters.corr_mismatch = false
			  if (clt_parameters.corr_mismatch || apply_corr || infinity_corr){ // added infinity_corr
				  clt_mismatch = new double [12][]; // What is 12?// not used in lwir
			  }
		  }
		  // Includes all 3 colors - will have zeros in unused
		  double [][] disparity_map = new double [ImageDtt.DISPARITY_TITLES.length][]; //[0] -residual disparity, [1] - orthogonal (just for debugging) last 4 - max pixel differences

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

		  double z_correction =  clt_parameters.z_correction;
		  if (clt_parameters.z_corr_map.containsKey(image_name)){
			  z_correction +=clt_parameters.z_corr_map.get(image_name);// not used in lwir
		  }
		  final double disparity_corr = (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);
		  double [][][][][][] clt_data = image_dtt.clt_aberrations_quad_corr(
				  clt_parameters.img_dtt,       // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
				  1,                            // final int  macro_scale, // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
				  tile_op,                      // per-tile operation bit codes
				  disparity_array,              // final double            disparity,
				  image_data, // double_stacks, // final double [][][]      imade_data, // first index - number of image in a quad
				  saturation_imp,               // boolean [][] saturation_imp, // (near) saturated pixels or null
				  // correlation results - final and partial
				  clt_corr_combo,               // [tp.tilesY][tp.tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				  clt_corr_partial,             // [tp.tilesY][tp.tilesX][pair][color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				  clt_mismatch,                 // [12][tp.tilesY * tp.tilesX] // transpose unapplied. null - do not calculate
				  disparity_map,                // [2][tp.tilesY * tp.tilesX]
				  texture_tiles,                // [tp.tilesY][tp.tilesX]["RGBA".length()][];
				  imp_quad[0].getWidth(),       // final int width,
				  clt_parameters.getFatZero(isMonochrome()),      // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
				  clt_parameters.corr_sym,
				  clt_parameters.corr_offset,
				  clt_parameters.corr_red,
				  clt_parameters.corr_blue,
				  clt_parameters.getCorrSigma(image_dtt.isMonochrome()),
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
///				  image_dtt.transform_size,
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
		  int first_color = isMonochrome()? ImageDtt.MONO_CHN:0;
		  if (debugLevel > -1){
			  System.out.println("clt_data.length="+clt_data.length+" clt_data[0].length="+clt_data[0].length
					  +" clt_data[0]["+first_color+"].length="+clt_data[0][first_color].length+" clt_data[0]["+first_color+"][0].length="+
					  clt_data[0][first_color][0].length);
		  }
		  // visualize texture tiles as RGBA slices
		  double [][] texture_nonoverlap = null;
		  double [][] texture_overlap = null;
		  String [] rbga_titles = {"red","blue","green","alpha"};
		  String [] rbga_weights_titles = {"red","blue","green","alpha","port0","port1","port2","port3","r-rms","b-rms","g-rms","w-rms"};
		  // In monochrome mode only G is used ImageDtt.MONO_CHN(==2), others are null
		  if (texture_tiles != null){
			  if (clt_parameters.show_nonoverlap){// not used in lwir
				  texture_nonoverlap = image_dtt.combineRBGATiles(
						  texture_tiles,                 // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}
///						  image_dtt.transform_size,
						  false,                         // when false - output each tile as 16x16, true - overlap to make 8x8
						  clt_parameters.sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
						  threadsMax,                    // maximal number of threads to launch
						  debugLevel);
				  sdfa_instance.showArrays(
						  texture_nonoverlap,
						  tilesX * (2 * image_dtt.transform_size),
						  tilesY * (2 * image_dtt.transform_size),
						  true,
						  image_name+sAux() + "-TXTNOL-D"+clt_parameters.disparity,
						  (clt_parameters.keep_weights?rbga_weights_titles:rbga_titles));

			  }
			  if (!infinity_corr && (clt_parameters.show_overlap || clt_parameters.show_rgba_color)){
				  int alpha_index = 3;
				  texture_overlap = image_dtt.combineRBGATiles(
						  texture_tiles,                 // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}
///						  image_dtt.transform_size,
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

				  if (!batch_mode && clt_parameters.show_overlap) {// not used in lwir
					  sdfa_instance.showArrays( // all but r-rms, b-rms
							  texture_overlap,
							  tilesX * image_dtt.transform_size,
							  tilesY * image_dtt.transform_size,
							  true,
							  image_name+sAux() + "-TXTOL-D"+clt_parameters.disparity,
							  (clt_parameters.keep_weights?rbga_weights_titles:rbga_titles));
				  }
				  if (!batch_mode && clt_parameters.show_rgba_color) {// not used in lwir
					  // for now - use just RGB. Later add option for RGBA
					  double [][] texture_rgb = {texture_overlap[0],texture_overlap[1],texture_overlap[2]};
					  double [][] texture_rgba = {texture_overlap[0],texture_overlap[1],texture_overlap[2],texture_overlap[3]};
//					  ImagePlus img_texture =
					  linearStackToColor(
							  clt_parameters,
							  colorProcParameters,
							  rgbParameters,
							  image_name+sAux()+"-texture", // String name,
							  "-D"+clt_parameters.disparity, //String suffix, // such as disparity=...
							  toRGB,
							  !this.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
							  true, // boolean saveShowIntermediate, // save/show if set globally
							  true, // boolean saveShowFinal,        // save/show result (color image?)
							  ((clt_parameters.alpha1 > 0)? texture_rgba: texture_rgb),
							  tilesX *  image_dtt.transform_size,
							  tilesY *  image_dtt.transform_size,
							  1.0,         // double scaleExposure, // is it needed?
							  debugLevel );
				  }
			  }
		  }
		  // visualize correlation results
		  // bo-b3 non-zero, r*, g* - zero
		  if (clt_corr_combo!=null){
			  if (disparity_map != null){
				  if (!batch_mode && clt_parameters.show_map &&  (debugLevel > -1)){
					  sdfa_instance.showArrays(
							  disparity_map,
							  tilesX,
							  tilesY,
							  true,
							  image_name+sAux()+"-DISP_MAP-D"+clt_parameters.disparity,
							  ImageDtt.DISPARITY_TITLES);
				  }
			  }

			  if (infinity_corr && (disparity_map != null)){// not used in lwir
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
								System.out.println("---- d{xy}[i] have disparity canceled, xy_mismatch will only reflect residual values---");
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
								  image_name+sAux() + "-inf_corr",
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
									  image_name+sAux() + "-inf_corr",
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


			  if (!batch_mode && !infinity_corr && clt_parameters.corr_show && (debugLevel > -1)){ // not used in lwir FALSE
				  double [][] corr_rslt = new double [clt_corr_combo.length][];
				  String [] titles = new String[clt_corr_combo.length]; // {"combo","sum"};
				  for (int i = 0; i< titles.length; i++) titles[i] = ImageDtt.TCORR_TITLES[i];
				  for (int i = 0; i<corr_rslt.length; i++) {
					  corr_rslt[i] = image_dtt.corr_dbg(
							  clt_corr_combo[i],
							  2*image_dtt.transform_size - 1,
							  clt_parameters.corr_border_contrast,
							  threadsMax,
							  debugLevel);
				  }
// all zeros
				  sdfa_instance.showArrays(
						  corr_rslt,
						  tilesX*(2*image_dtt.transform_size),
						  tilesY*(2*image_dtt.transform_size),
						  true,
						  image_name+sAux()+"-CORR-D"+clt_parameters.disparity,
						  titles );
			  }

			  if (!batch_mode && !infinity_corr && (clt_corr_partial!=null)){ // not used in lwir FALSE
				  if (debugLevel > -1){ // -1
					  String [] allColorNames = {"red","blue","green","combo"};
					  String [] titles = new String[clt_corr_partial.length];
					  for (int i = 0; i < titles.length; i++){
						  titles[i]=allColorNames[i % allColorNames.length]+"_"+(i / allColorNames.length);
					  }
					  double [][] corr_rslt_partial = image_dtt.corr_partial_dbg(
							  clt_corr_partial,
							  2*image_dtt.transform_size - 1,	//final int corr_size,
							  4,	// final int pairs,
							  4,    // final int colors,
							  clt_parameters.corr_border_contrast,
							  threadsMax,
							  debugLevel);
					  // titles.length = 15, corr_rslt_partial.length=16!
					  System.out.println("corr_rslt_partial.length = "+corr_rslt_partial.length+", titles.length = "+titles.length);
					  sdfa_instance.showArrays( // out of boundary 15
							  corr_rslt_partial,
							  tilesX*(2*image_dtt.transform_size),
							  tilesY*(2*image_dtt.transform_size),
							  true,
							  image_name+sAux()+"-PART_CORR-D"+clt_parameters.disparity);
//							  titles);
				  }
			  }
		  }
		  double [][][] iclt_data = new double [clt_data.length][][];
		  if (!infinity_corr && (clt_parameters.gen_chn_img || clt_parameters.gen_4_img || clt_parameters.gen_chn_stacks)) {
//			  ImagePlus [] imps_RGB = new ImagePlus[clt_data.length];
			  for (int iQuad = 0; iQuad < clt_data.length; iQuad++){

//				  String title=name+"-"+String.format("%02d", iQuad);
				  //				  String titleFull=title+"-SPLIT-D"+clt_parameters.disparity;

				  if (clt_parameters.getCorrSigma(image_dtt.isMonochrome()) > 0){ // no filter at all
					  for (int chn = 0; chn < clt_data[iQuad].length; chn++) if (clt_data[iQuad][chn] != null){
						  image_dtt.clt_lpf(
								  clt_parameters.getCorrSigma(image_dtt.isMonochrome()),
								  clt_data[iQuad][chn],
///								  image_dtt.transform_size,
								  threadsMax,
								  debugLevel);
					  }
				  }

				  if (debugLevel > 0){
					  System.out.println("--tilesX="+tilesX);
					  System.out.println("--tilesY="+tilesY);
				  }
				  if (!batch_mode && (debugLevel > 0)){ // FALSE
					  double [][] clt = new double [clt_data[iQuad].length*4][];
					  for (int chn = 0; chn < clt_data[iQuad].length; chn++) if (clt_data[iQuad][chn] != null){
						  double [][] clt_set = image_dtt.clt_dbg(
								  clt_data [iQuad][chn],
								  threadsMax,
								  debugLevel);
						  for (int ii = 0; ii < clt_set.length; ii++) clt[chn*4+ii] = clt_set[ii];
					  }

					  if (debugLevel > 0){
						  sdfa_instance.showArrays(clt,
								  tilesX*image_dtt.transform_size,
								  tilesY*image_dtt.transform_size,
								  true,
								  results[iQuad].getTitle()+"-CLT-D"+clt_parameters.disparity);
					  }
				  }
				  iclt_data[iQuad] = new double [clt_data[iQuad].length][];

				  for (int ncol=0; ncol<iclt_data[iQuad].length;ncol++) if (clt_data[iQuad][ncol] != null) {
					  iclt_data[iQuad][ncol] = image_dtt.iclt_2d(
							  clt_data[iQuad][ncol],           // scanline representation of dcd data, organized as dct_size x dct_size tiles
///							  image_dtt.transform_size,  // final int
							  clt_parameters.clt_window,      // window_type
							  15,                             // clt_parameters.iclt_mask,       //which of 4 to transform back
							  0,                              // clt_parameters.dbg_mode,        //which of 4 to transform back
							  threadsMax,
							  debugLevel);

				  }
				  if (clt_parameters.gen_chn_stacks) sdfa_instance.showArrays(
//				  if (clt_parameters.gen_chn_stacks || true) sdfa_instance.showArrays(
						  iclt_data[iQuad],
						  (tilesX + 0) * image_dtt.transform_size,
						  (tilesY + 0) * image_dtt.transform_size,
						  true,
						  results[iQuad].getTitle()+"-ICLT-RGB-D"+clt_parameters.disparity);
			  } // end of generating shifted channel images


			  // Use iclt_data here for LWIR autorange
			  if (colorProcParameters.isLwir() && colorProcParameters.lwir_autorange) {
				  double rel_low =  colorProcParameters.lwir_low;
				  double rel_high = colorProcParameters.lwir_high;
				  if (!Double.isNaN(getLwirOffset())) {
					  rel_low -=  getLwirOffset();
					  rel_high -= getLwirOffset();
				  }
				  double [] cold_hot =  autorange(
						  iclt_data, // double [][][] iclt_data, //  [iQuad][ncol][i] - normally only [][2][] is non-null
						  rel_low, // double hard_cold,// matches data, DC (this.lwir_offset)  subtracted
						  rel_high, // double hard_hot, // matches data, DC (this.lwir_offset)  subtracted
						  colorProcParameters.lwir_too_cold, // double too_cold, // pixels per image
						  colorProcParameters.lwir_too_hot, // double too_hot,  // pixels per image
						  1024); // int num_bins)
				  if (cold_hot != null) {
					  if (!Double.isNaN(getLwirOffset())) {
						  cold_hot[0] += getLwirOffset();
						  cold_hot[1] += getLwirOffset();
					  }
				  }
				  setColdHot(cold_hot); // will be used for shifted images and for texture tiles
			  }
			  ImagePlus [] imps_RGB = new ImagePlus[clt_data.length];
				  for (int iQuad = 0; iQuad < clt_data.length; iQuad++){
					  if (!clt_parameters.gen_chn_img) continue;
					  String title=String.format("%s%s-%02d",image_name, sAux(), iQuad);
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
							  iclt_data[iQuad],
							  tilesX *  image_dtt.transform_size,
							  tilesY *  image_dtt.transform_size,
							  scaleExposures[iQuad], // double scaleExposure, // is it needed?
							  debugLevel );
				  }


			  if (clt_parameters.gen_chn_img) {
				  // combine to a sliced color image
				  int [] slice_seq = {0,1,3,2}; //clockwise
				  int width = imps_RGB[0].getWidth();
				  int height = imps_RGB[0].getHeight();
				  ImageStack array_stack=new ImageStack(width,height);
				  for (int i = 0; i<slice_seq.length; i++){
					  if (imps_RGB[slice_seq[i]] != null) {
						  array_stack.addSlice("port_"+slice_seq[i], imps_RGB[slice_seq[i]].getProcessor().getPixels());
					  } else { // not used in lwir
						  array_stack.addSlice("port_"+slice_seq[i], results[slice_seq[i]].getProcessor().getPixels());
					  }
				  }
				  ImagePlus imp_stack = new ImagePlus(image_name+sAux()+"-SHIFTED-D"+clt_parameters.disparity, array_stack);
				  imp_stack.getProcessor().resetMinAndMax();
				  if (!batch_mode) {
					  imp_stack.updateAndDraw(); // not used in lwir
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
						  image_name, // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
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

				  String model_path= correctionsParameters.selectX3dDirectory(
						  image_name, // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
						  null,
						  true,  // smart,
						  true);  //newAllowed, // save

				  createThumbNailImage(
						  imps_RGB[0],
						  model_path,
						  "thumb"+sAux(),
						  debugLevel);


			  }
		  }
		  return results;
	  }

	  public ImagePlus [] processCLTQuadCorrTestERS(
			  ImagePlus [] imp_quad, // should have properties "name"(base for saving results), "channel","path"
			  boolean [][] saturation_imp, // (near) saturated pixels or null
			  CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  ColorProcParameters                              colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters         channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
//			  int              convolveFFTSize, // 128 - fft size, kernel size should be size/2
			  double []	       scaleExposures, // probably not needed here
			  final boolean    apply_corr, // calculate and apply additional fine geometry correction
			  final boolean    infinity_corr, // calculate and apply geometry correction at infinity
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  int              debugLevel){

		  if (debugLevel > -2) { // -1
			  debugLevel = clt_parameters.ly_debug_level;
		  }

		  String path= (String) imp_quad[0].getProperty("path");

		  ImagePlus [] results = new ImagePlus[imp_quad.length];
		  for (int i = 0; i < results.length; i++) {
			  results[i] = imp_quad[i];
			  results[i].setTitle(results[i].getTitle()+"RAW");
		  }
		  if (debugLevel>1) System.out.println("processing: "+path);
// 08/12/2020 Moved to condifuinImageSet - remove setTiles?
/*
   	     setTiles (imp_quad[0], // set global tp.tilesX, tp.tilesY
 				  clt_parameters,
				  threadsMax);
*/
		  final int tilesX=tp.getTilesX(); // imp_quad[0].getWidth()/clt_parameters.transform_size;
		  final int tilesY=tp.getTilesY(); // imp_quad[0].getHeight()/clt_parameters.transform_size;

		  final int clustersX= (tilesX + clt_parameters.tileStep - 1) / clt_parameters.tileStep;
		  final int clustersY= (tilesY + clt_parameters.tileStep - 1) / clt_parameters.tileStep;

		  double [][] dsxy = new double[clustersX * clustersY][];
		  ImagePlus imp_sel = WindowManager.getCurrentImage();
		  if ((imp_sel == null) || (imp_sel.getStackSize() != ExtrinsicAdjustment.INDX_LENGTH)) {
			  System.out.println("No image / wrong image selected, bailing out");
			  return null;
		  } else {
			  System.out.println("Image: "+imp_sel.getTitle());
				int width = imp_sel.getWidth();
				int height = imp_sel.getHeight();
				if ((width != clustersX) || (height != clustersY)) {
					  System.out.println(String.format("Image size mismatch: width=%d (%d), height=%d(%d)",
							  width, clustersX, height, clustersY));
					  return null;
				}
				ImageStack stack_sel = imp_sel.getStack();
				float [][] fpixels = new float [ExtrinsicAdjustment.INDX_LENGTH][];
				for (int i = 0; i < fpixels.length; i++) {
					fpixels[i] = (float[]) stack_sel.getPixels(i+1);
				}
				for (int tile = 0; tile < dsxy.length; tile ++) {
					if (fpixels[1][tile] > 0.0) {
						dsxy[tile] = new double[fpixels.length];
						for (int i = 0; i < fpixels.length; i++) {
							dsxy[tile][i] = fpixels[i][tile];
						}
					}
				}
		  }
		  ExtrinsicAdjustment ea = new ExtrinsicAdjustment(
				  geometryCorrection, // GeometryCorrection gc,
				  clt_parameters.tileStep,   // int         clusterSize,
				  clustersX, // 	int         clustersX,
				  clustersY); // int         clustersY);

		  double [] old_new_rms = new double[2];
		  boolean apply_extrinsic = (clt_parameters.ly_corr_scale != 0.0);

		  GeometryCorrection.CorrVector corr_vector =   ea.solveCorr (
				  clt_parameters.ly_marg_fract, 	// double      marg_fract,        // part of half-width, and half-height to reduce weights
				  clt_parameters.ly_inf_en,      // boolean     use_disparity,     // adjust disparity-related extrinsics
				  clt_parameters.ly_aztilt_en,   // boolean     use_aztilts,       // Adjust azimuths and tilts excluding disparity
				  clt_parameters.ly_diff_roll_en,//boolean     use_diff_rolls,    // Adjust differential rolls (3 of 4 angles)
//				  clt_parameters.ly_inf_force,   //	boolean     force_convergence, // if true try to adjust convergence (disparity, symmetrical parameter 0) even with no disparity
				  clt_parameters.ly_min_forced,  //	int         min_num_forced,    // minimal number of clusters with forced disparity to use it
				  // data, using just radial distortions
				  clt_parameters.ly_com_roll,    //boolean     common_roll,       // Enable common roll (valid for high disparity range only)
				  clt_parameters.ly_focalLength, //boolean     corr_focalLength,  // Correct scales (focal length temperature? variations)
				  clt_parameters.ly_ers_rot,     //	boolean     ers_rot,           // Enable ERS correction of the camera rotation
				  clt_parameters.ly_ers_forw,    //	boolean     ers_forw,          // Enable ERS correction of the camera linear movement in z direction
				  clt_parameters.ly_ers_side,    //	boolean     ers_side,          // Enable ERS correction of the camera linear movement in x direction
				  clt_parameters.ly_ers_vert,    //	boolean     ers_vert,          // Enable ERS correction of the camera linear movement in y direction
				  // add balancing-related here?
				  clt_parameters.ly_par_sel, // 	int         manual_par_sel,    // Manually select the parameter mask bit 0 - sym0, bit1 - sym1, ... (0 - use boolean flags, != 0 - ignore boolean flags)
				  1.0, // 	double      weight_disparity,
				  1.0, // 	double      weight_lazyeye,
				  dsxy, // double [][] measured_dsxy,
				  null, //	boolean [] force_disparity,    // boolean [] force_disparity,
				  false, // 	boolean     use_main, // corr_rots_aux != null;
				  geometryCorrection.getCorrVector(), // GeometryCorrection.CorrVector corr_vector,
				  old_new_rms, // double [] old_new_rms, // should be double[2]
				  debugLevel); //  + 5);// int debugLevel)

		  if (debugLevel > -2){
			  System.out.println("Old extrinsic corrections:");
			  System.out.println(geometryCorrection.getCorrVector().toString());
		  }
		  if (corr_vector != null) {
			  GeometryCorrection.CorrVector diff_corr = corr_vector.diffFromVector(geometryCorrection.getCorrVector());
			  if (debugLevel > -2){
					  System.out.println("New extrinsic corrections:");
					  System.out.println(corr_vector.toString());

					  System.out.println("Increment extrinsic corrections:");
					  System.out.println(diff_corr.toString());
					  // System.out.println("Correction scale = "+clt_parameters.ly_corr_scale);

			  }

			  if (apply_extrinsic){
				  geometryCorrection.setCorrVector(corr_vector) ;
/*
				  geometryCorrection.getCorrVector().incrementVector(diff_corr, clt_parameters.ly_corr_scale);
				  System.out.println("New (with correction scale applied) extrinsic corrections:");
				  System.out.println(geometryCorrection.getCorrVector().toString());
*/
				  System.out.println("Extrinsic correction updated (can be disabled by setting clt_parameters.ly_corr_scale = 0.0) ");

			  } else {
				  System.out.println("Correction is not applied according clt_parameters.ly_corr_scale == 0.0) ");
			  }
		  } else {
			  if (debugLevel > -3){
				  System.out.println("LMA failed");
			  }
		  }
		  double [][][] mismatch_corr_coefficients = new double [1][2][];
		  mismatch_corr_coefficients[0][0] = geometryCorrection.getCorrVector().toSymArray(null);
		  mismatch_corr_coefficients[0][1] = old_new_rms;
		  return null;
	  }

	  public ImagePlus [] processCLTQuadCorrTest(
			  ImagePlus [] imp_quad, // should have properties "name"(base for saving results), "channel","path"
			  boolean [][] saturation_imp, // (near) saturated pixels or null
			  CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  ColorProcParameters                              colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters         channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
//			  int              convolveFFTSize, // 128 - fft size, kernel size should be size/2
			  double []	       scaleExposures, // probably not needed here
			  final boolean    apply_corr, // calculate and apply additional fine geometry correction
			  final boolean    infinity_corr, // calculate and apply geometry correction at infinity
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel){
		  final boolean      batch_mode = clt_parameters.batch_run; //disable any debug images
		  String name=this.correctionsParameters.getModelName((String) imp_quad[0].getProperty("name"));
		  //		int channel= Integer.parseInt((String) imp_src.getProperty("channel"));
		  String path= (String) imp_quad[0].getProperty("path");

		  ImagePlus [] results = new ImagePlus[imp_quad.length];
		  for (int i = 0; i < results.length; i++) {
			  results[i] = imp_quad[i];
			  results[i].setTitle(results[i].getTitle()+"RAW");
		  }
		  if (debugLevel>1) System.out.println("processing: "+path);
/*		// 08/12/2020 Moved to condifuinImageSet				  

		  double [][][] double_stacks = new double [imp_quad.length][][];
		  for (int i = 0; i < double_stacks.length; i++){
			  double_stacks[i] = eyesisCorrections.bayerToDoubleStack(
					  imp_quad[i], // source Bayer image, linearized, 32-bit (float))
					  null, // no margins, no oversample
					  this.is_mono);
		  }

		  for (int i = 0; i < double_stacks.length; i++){
			  if ( double_stacks[i].length > 2) {
				  for (int j =0 ; j < double_stacks[i][0].length; j++){
					  double_stacks[i][2][j]*=0.5; // Scale green 0.5 to compensate more pixels than R,B
				  }
			  } else {
				  for (int j =0 ; j < double_stacks[i][0].length; j++){
					  double_stacks[i][0][j]*=1.0; // Scale mono by 1/4 - to have the same overall "gain" as for bayer
				  }
			  }
		  }

		  setTiles (imp_quad[0], // set global tp.tilesX, tp.tilesY
				  clt_parameters,
				  threadsMax);
*/		  
		  ImageDtt image_dtt = new ImageDtt(
				  clt_parameters.transform_size,
				  isMonochrome(),
				  isLwir(),
				  clt_parameters.getScaleStrength(isAux()));
		  
		  
		  double [][] disparity_array = tp.setSameDisparity(clt_parameters.disparity); // 0.0); // [tp.tilesY][tp.tilesX] - individual per-tile expected disparity
		  ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
		  ImagePlus imp_sel = WindowManager.getCurrentImage();
		  if ((imp_sel == null) || (imp_sel.getStackSize() != 3)) {
			  System.out.println("No image / wrong image selected, using infinity");
		  } else {
			  System.out.println("Image: "+imp_sel.getTitle());
				int width = imp_sel.getWidth();
				int height = imp_sel.getHeight();
				if ((width != tp.getTilesX()) || (height != tp.getTilesY())) {
					  System.out.println(String.format("Image size mismatch: width=%d (%d), height=%d(%d)",
							  width, tp.getTilesX(), height, tp.getTilesY()));
					  return null;
				}
				ImageStack stack_sel = imp_sel.getStack();
				float [] img_disparity = (float[]) stack_sel.getPixels(1); // first stack is disparity
				int indx = 0;
				for (int i = 0; i< disparity_array.length; i++){
					for (int j = 0; j< disparity_array[0].length; j++){
						double d = img_disparity[indx++];
						if (!Double.isNaN(d)) { // treat NaN as 0
							disparity_array[i][j] += d;
						}
					}
				}
		  }


		  // temporary setting up tile task file (one integer per tile, bitmask
		  // for testing defined for a window, later the tiles to process will be calculated based on previous passes results

		  int [][]    tile_op = tp.setSameTileOp(clt_parameters,  clt_parameters.tile_task_op, debugLevel);
//		  double [][] disparity_array = tp.setSameDisparity(clt_parameters.disparity); // 0.0); // [tp.tilesY][tp.tilesX] - individual per-tile expected disparity

		  //TODO: Add array of default disparity - use for combining images in force disparity mode (no correlation), when disparity is predicted from other tiles

		  double [][][][]     clt_corr_combo =   null;
		  double [][][][][]   clt_corr_partial = null; // [tp.tilesY][tp.tilesX][pair][color][(2*transform_size-1)*(2*transform_size-1)]
		  double [][]         clt_mismatch =     null; // [3*4][tp.tilesY * tp.tilesX] // transpose unapplied
		  double [][][][]     texture_tiles =    null; // [tp.tilesY][tp.tilesX]["RGBA".length()][]; // tiles will be 16x16, 2 visualization mode full 16 or overlapped
		  // undecided, so 2 modes of combining alpha - same as rgb, or use center tile only
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();

		  if (clt_parameters.correlate){ // true
			  clt_corr_combo =    new double [ImageDtt.TCORR_TITLES.length][tilesY][tilesX][];
			  texture_tiles =     new double [tilesY][tilesX][][]; // ["RGBA".length()][];
			  for (int i = 0; i < tilesY; i++){
				  for (int j = 0; j < tilesX; j++){
					  for (int k = 0; k<clt_corr_combo.length; k++){
						  clt_corr_combo[k][i][j] = null;
					  }
					  texture_tiles[i][j] = null;
				  }
			  }
			  if (!infinity_corr && clt_parameters.corr_keep){ // true
				  clt_corr_partial = new double [tilesY][tilesX][][][];
				  for (int i = 0; i < tilesY; i++){
					  for (int j = 0; j < tilesX; j++){
						  clt_corr_partial[i][j] = null;
					  }
				  }
			  } // clt_parameters.corr_mismatch = false
			  if (clt_parameters.corr_mismatch || apply_corr || infinity_corr){ // added infinity_corr
				  clt_mismatch = new double [12][]; // What is 12?// not used in lwir
			  }
		  }
		  // Includes all 3 colors - will have zeros in unused
		  double [][] disparity_map = new double [ImageDtt.DISPARITY_TITLES.length][]; //[0] -residual disparity, [1] - orthogonal (just for debugging) last 4 - max pixel differences

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

		  double z_correction =  clt_parameters.z_correction;
		  if (clt_parameters.z_corr_map.containsKey(name)){
			  z_correction +=clt_parameters.z_corr_map.get(name);// not used in lwir
		  }
		  final double disparity_corr = (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);
		  double [][] lazy_eye_data = image_dtt.cltMeasureLazyEye(
				  clt_parameters.img_dtt,       // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
				  tile_op,                      // per-tile operation bit codes
				  disparity_array,              // final double            disparity,
				  image_data, // double_stacks,                // final double [][][]      imade_data, // first index - number of image in a quad
				  saturation_imp,               // boolean [][] saturation_imp, // (near) saturated pixels or null
				  clt_mismatch,                 // [12][tp.tilesY * tp.tilesX] // transpose unapplied. null - do not calculate
				  disparity_map,                // [2][tp.tilesY * tp.tilesX]
				  imp_quad[0].getWidth(),       // final int width,
				  clt_parameters.getFatZero(isMonochrome()),      // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
				  clt_parameters.corr_red,
				  clt_parameters.corr_blue,
				  clt_parameters.getCorrSigma(image_dtt.isMonochrome()),
				  min_corr_selected, // 0.0001; // minimal correlation value to consider valid
				  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
				  null,                          // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
				  clt_kernels,                   // final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				  clt_parameters.kernel_step,
				  clt_parameters.clt_window,
				  shiftXY, //
				  disparity_corr, // final double              disparity_corr, // disparity at infinity
				  clt_parameters.shift_x,          // final int               shiftX, // shift image horizontally (positive - right) - just for testing
				  clt_parameters.shift_y,          // final int               shiftY, // shift image vertically (positive - down)
				  clt_parameters.tileStep,         // 	final int                 tileStep, // process tileStep x tileStep cluster of tiles when adjusting lazy eye parameters
				  clt_parameters.tileX, // -1234, // clt_parameters.tileX,         // final int               debug_tileX,
				  clt_parameters.tileY,         // final int               debug_tileY, -1234 will cause port coordinates debug images
				  threadsMax,
				  debugLevel);

		  if (lazy_eye_data != null) {
				  int clustersX= (tilesX + clt_parameters.tileStep - 1) / clt_parameters.tileStep;
				  int clustersY= (tilesY + clt_parameters.tileStep - 1) / clt_parameters.tileStep;
				  double [][] dbg_cluster = new double [ExtrinsicAdjustment.INDX_LENGTH][clustersY * clustersX];
				  for (int n = 0; n < lazy_eye_data.length; n++) {
					  if ((lazy_eye_data[n] != null) && (lazy_eye_data[n][ExtrinsicAdjustment.INDX_STRENGTH] >= clt_parameters.img_dtt.lma_diff_minw)) {
						  for (int i = 0; i < ExtrinsicAdjustment.INDX_LENGTH; i++ ) {
							  if (i == ExtrinsicAdjustment.INDX_STRENGTH) {
								  dbg_cluster[i][n] = lazy_eye_data[n][i] - clt_parameters.img_dtt.lma_diff_minw; // strength
							  } else {
								  dbg_cluster[i][n] = lazy_eye_data[n][i];
							  }
						  }
					  } else {
						  for (int i = 0; i < ExtrinsicAdjustment.INDX_LENGTH; i++ ) {
							  if (i == ExtrinsicAdjustment.INDX_STRENGTH) {
								  dbg_cluster[i][n] = 0.0; // strength
							  } else {
								  dbg_cluster[i][n] = Double.NaN;
							  }
						  }
					  }
				  }
				  //clt_parameters.img_dtt.lma_diff_minw
				  sdfa_instance.showArrays(
						  dbg_cluster,
						  clustersX,
						  clustersY,
						  true,
						  name+sAux()+"-CLT_MISMATCH-D"+clt_parameters.disparity+"_"+clt_parameters.tileStep+"x"+clt_parameters.tileStep,
						  ExtrinsicAdjustment.DATA_TITLES);

				  if (disparity_map != null){
					  int target_index = ImageDtt.DISPARITY_INDEX_INT;
					  int cm_index = ImageDtt.DISPARITY_INDEX_INT+1;
					  int lma_index = ImageDtt.DISPARITY_INDEX_VERT;
					  int strength_index = ImageDtt.DISPARITY_STRENGTH_INDEX;
					  double [][] scan_maps = new double[3][tilesX*tilesY];
					  for (int i = 0; i < scan_maps[0].length; i++) {
						  scan_maps[1][i] = disparity_map[lma_index][i];
						  scan_maps[2][i] = disparity_map[strength_index][i];
						  if (Double.isNaN(disparity_map[lma_index][i])) {
							  if (Double.isNaN(disparity_map[cm_index][i])) {
								  scan_maps[0][i] = disparity_map[target_index][i]; // TODO: add offset calculated from neighbours
							  } else {
								  scan_maps[0][i] = disparity_map[cm_index][i];
							  }

						  } else {
							  scan_maps[0][i] = disparity_map[lma_index][i];
						  }
					  }
					  String [] titles3 = {"combo", "lma", "strength"};
					  sdfa_instance.showArrays(
							  scan_maps,
							  tilesX,
							  tilesY,
							  true,
							  name+sAux()+"-DISP_MAP-D"+clt_parameters.disparity+"-CLT",
							  titles3);

				  }
		  }

/**/
		  if (disparity_map != null){
			  if (!batch_mode && clt_parameters.show_map &&  (debugLevel > -2)){
				  sdfa_instance.showArrays(
						  disparity_map,
						  tilesX,
						  tilesY,
						  true,
						  name+sAux()+"-DISP_MAP-D"+clt_parameters.disparity,
						  ImageDtt.DISPARITY_TITLES);
			  }
			  /*
			  if (clt_mismatch != null) {
				  sdfa_instance.showArrays(
						  clt_mismatch,
						  tilesX,
						  tilesY,
						  true,
						  name+sAux()+"-CLT_MISMATCH-D"+clt_parameters.disparity);
				  //					  ImageDtt.DISPARITY_TITLES);
				  double min_bw = 0.005;
//				  double [][] mismatch_w = new double [clt_mismatch.length][tilesX*tilesY];
				  double mismatch_sigma = clt_parameters.tileStep * clt_parameters.img_dtt.lma_diff_sigma;
				  DoubleGaussianBlur gb = new DoubleGaussianBlur();
				  for (int ntile = 0; ntile < clt_mismatch[0].length; ntile++) {
					  double w = disparity_map[ImageDtt.IMG_DIFF0_INDEX+0][ntile];
					  if ( !(w >=  clt_parameters.img_dtt.lma_diff_minw)) {
						  w = 0.0;
					  }
					  for (int n = 0; n < clt_mismatch.length/3; n++) {
						  if (w <= 0.0) {
							  clt_mismatch[3*n+0][ntile] = Double.NaN;
							  clt_mismatch[3*n+1][ntile] = Double.NaN;
							  clt_mismatch[3*n+2][ntile] = 0.0;
						  } else {
							  clt_mismatch[3*n+2][ntile]  = w;
						  }
					  }
				  }
				  sdfa_instance.showArrays(
						  clt_mismatch,
						  tilesX,
						  tilesY,
						  true,
						  name+sAux()+"-CLT_MISMATCH-FILTERED-D"+clt_parameters.disparity);


				  for (int ntile = 0; ntile < clt_mismatch[0].length; ntile++) {
					  double w = disparity_map[ImageDtt.IMG_DIFF0_INDEX+0][ntile];
					  if ( !(w >=  clt_parameters.img_dtt.lma_diff_minw)) {
						  w = 0.0;
					  }
					  for (int n = 0; n < clt_mismatch.length/3; n++) {
						  if (w <= 0.0) {
							  clt_mismatch[3*n+0][ntile] = 0.0;
							  clt_mismatch[3*n+1][ntile] = 0.0;
							  clt_mismatch[3*n+2][ntile] = 0.0;
						  } else {
							  clt_mismatch[3*n+0][ntile] *= w;
							  clt_mismatch[3*n+1][ntile] *= w;
							  clt_mismatch[3*n+2][ntile]  = w;
						  }
					  }
				  }
				  for (int n = 0; n < clt_mismatch.length/3; n++) {
					  gb.blurDouble(clt_mismatch[3*n+0] , tilesX, tilesY, mismatch_sigma, mismatch_sigma, 0.01);
					  gb.blurDouble(clt_mismatch[3*n+1] , tilesX, tilesY, mismatch_sigma, mismatch_sigma, 0.01);
					  gb.blurDouble(clt_mismatch[3*n+2] , tilesX, tilesY, mismatch_sigma, mismatch_sigma, 0.01);
					  for (int ntile = 0; ntile < clt_mismatch[0].length; ntile++) {
						  if (clt_mismatch[3*n+2][ntile] >= min_bw) {
							  clt_mismatch[3*n+0][ntile] /= clt_mismatch[3*n+2][ntile];
							  clt_mismatch[3*n+1][ntile] /= clt_mismatch[3*n+2][ntile];
							  if (n>0) {
								  double w = disparity_map[ImageDtt.IMG_DIFF0_INDEX+0][ntile] - clt_parameters.img_dtt.lma_diff_minw;
								  if ( w <  0.0) {
									  w = 0.0;
								  }
								  clt_mismatch[3*n+2][ntile]  = w;
							  }
						  } else {
							  clt_mismatch[3*n+0][ntile] = Double.NaN;
							  clt_mismatch[3*n+1][ntile] = Double.NaN;
							  clt_mismatch[3*n+2][ntile] = 0.0;
						  }
					  }
				  }
				  sdfa_instance.showArrays(
						  clt_mismatch,
						  tilesX,
						  tilesY,
						  true,
						  name+sAux()+"-CLT_MISMATCH-BLUR-D"+clt_parameters.disparity);
			  } */
		  }
/**/

		  return results;
	  }

	  double [][] resizeGridTexture( // USED in lwir
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
	  public int [] getLwirHistogram( // USED in lwir
			  double [] data,
			  double    hard_cold,
			  double    hard_hot,
			  int       num_bins) {
			  int [] hist = new int [num_bins];
			  double k = num_bins / (hard_hot - hard_cold);
			  for (double d:data) {
				  int bin = (int) ((d - hard_cold)*k);
				  if (bin < 0) bin = 0;
				  else if (bin >= num_bins) bin = (num_bins -1);
				  hist[bin]++;
			  }
			  return hist;
		  }
		  public int [] addHist( // USED in lwir
				  int [] this_hist,
				  int [] other_hist) {
			  for (int i = 0; i < this_hist.length; i++) {
				  this_hist[i] += other_hist[i];
			  }
			  return this_hist;
		  }
		  // get low/high (soft min/max) from the histogram
		  // returns value between 0.0 (low histogram limit and 1.0 - high histgram limit
		  public double getMarginFromHist( // USED in lwir
				  int [] hist, // histogram
				  double cumul_val, // cummulative number of items to be ignored
				  boolean high_marg) { // false - find low margin(output ~0.0) , true - find high margin (output ~1.0)
			  int n = 0;
			  int n_prev = 0;
			  int bin;
			  double s = 1.0 / hist.length;
			  double v;
			  if (high_marg) {
				  for (bin = hist.length -1; bin >= 0; bin--) {
					  n_prev = n;
					  n+= hist[bin];
					  if (n > cumul_val) break;
				  }
				  if (n <= cumul_val) { // not used in lwir
					  v =  0.0; // cumul_val > total number of samples
				  } else {
					  v = s* (bin + 1 - (cumul_val - n_prev)/(n - n_prev));
				  }

			  } else {
				  for (bin = 0; bin < hist.length; bin++) {
					  n_prev = n;
					  n+= hist[bin];
					  if (n > cumul_val) break;
				  }
				  if (n <= cumul_val) { // not used in lwir
					  v =  1.0; // cumul_val > total number of samples
				  } else {
					  v = s * (bin + (cumul_val - n_prev)/(n - n_prev));
				  }
			  }
			  return v;
		  }

		  public double [] autorange( // USED in lwir
				  double [][][] iclt_data, //  [iQuad][ncol][i] - normally only [][2][] is non-null
				  double hard_cold,// matches data, DC (this.lwir_offset)  subtracted
				  double hard_hot, // matches data, DC (this.lwir_offset)  subtracted
				  double too_cold, // pixels per image
				  double too_hot,  // pixels per image
				  int num_bins) {
			  int ncol;
			  for (ncol = 0; ncol < iclt_data[0].length; ncol++) {
				  if (iclt_data[0][ncol] != null) break;
			  }
			  too_cold *= iclt_data.length;
			  too_hot *= iclt_data.length;
			  int [] hist = null;
			  for (int iQuad = 0; iQuad < iclt_data.length; iQuad++) {
				  int [] this_hist = getLwirHistogram(
						  iclt_data[iQuad][ncol], // double [] data,
						  hard_cold,
						  hard_hot,
						  num_bins);
				  if (hist == null) {
					  hist = this_hist;
				  } else {
					  addHist(
							  hist,
							  this_hist);
				  }
			  }
			  double [] rel_lim = {
					  getMarginFromHist(
							  hist, // histogram
							  too_cold, // double cumul_val, // cummulative number of items to be ignored
							  false), // boolean high_marg)
					  getMarginFromHist(
							  hist, // histogram
							  too_hot, // double cumul_val, // cummulative number of items to be ignored
							  true)}; // boolean high_marg)
			  double [] abs_lim = {
					  rel_lim[0] * (hard_hot - hard_cold) + hard_cold,
					  rel_lim[1] * (hard_hot - hard_cold) + hard_cold,
			  };
			  return abs_lim;
		  }



	  // float
	  public ImagePlus linearStackToColor( // not used in lwir
			  CLTParameters         clt_parameters,
			  ColorProcParameters   colorProcParameters,
			  EyesisCorrectionParameters.RGBParameters         rgbParameters,
			  String name,
			  String suffix, // such as disparity=...
			  boolean toRGB,
			  boolean bpp16, // 16-bit per channel color mode for result
			  boolean saveShowIntermediate, // save/show if set globally
			  boolean saveShowFinal,        // save/show result (color image?)
			  float [][] iclt_data,
			  int width, // int tilesX,
			  int height, // int tilesY,
			  double scaleExposure,
			  int debugLevel
			  )
	  {
		  ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
		  // convert to ImageStack of 3 slices
		  String [] sliceNames = {"red", "blue", "green"};
		  int green_index = 2;
		  float [][] rbg_in = {iclt_data[0],iclt_data[1],iclt_data[2]};
		  float []   alpha = null; // (0..1.0)
//		  float [][] rgb_in = {iclt_data[0],iclt_data[1],iclt_data[2]};
		  if (iclt_data.length > 3) alpha = iclt_data[3];
		  if (isLwir()) {
			  String [] rgb_titles =  {"red","green","blue"};
			  String [] rgba_titles = {"red","green","blue","alpha"};
			  String [] titles = (alpha == null) ? rgb_titles : rgba_titles;
			  int num_slices = (alpha == null) ? 3 : 4;
			  double mn = colorProcParameters.lwir_low;
			  double mx = colorProcParameters.lwir_high;
			  double [] cold_hot = getColdHot();
			  if (cold_hot != null) {
				  mn = cold_hot[0];
				  mx = cold_hot[1];
			  }

			  double offset = getLwirOffset();
			  if (!Double.isNaN(offset)) {
				  mn -=  offset;
				  mx -=  offset;
			  }

			  ThermalColor tc = new ThermalColor(
					  colorProcParameters.lwir_palette,
					  mn,
					  mx,
					  255.0);
			  float [][] rgba = new float [num_slices][];
			  for (int i = 0; i < 3; i++) rgba[i] = new float [iclt_data[green_index].length];
			  for (int i = 0; i < rbg_in[green_index].length; i++) {
				  if (i == 700) {
					  System.out.println("linearStackToColor(): i="+i);
				  }
				  float [] rgb = tc.getRGB(iclt_data[green_index][i]);
				  rgba[0][i] = rgb[0]; // red
				  rgba[1][i] = rgb[1]; // green
				  rgba[2][i] = rgb[2]; // blue
			  }
			  if (alpha != null) {
				  rgba[3] = alpha; // 0..1
			  }
			  ImageStack stack = sdfa_instance.makeStack(
					  rgba,       // iclt_data,
					  width,      // (tilesX + 0) * clt_parameters.transform_size,
					  height,     // (tilesY + 0) * clt_parameters.transform_size,
					  titles,     // or use null to get chn-nn slice names
					  true);      // replace NaN with 0.0
			  ImagePlus imp_rgba =  EyesisCorrections.convertRGBAFloatToRGBA32(
					  stack,   // ImageStack stackFloat, //r,g,b,a
					  //						name+"ARGB"+suffix, // String title,
					  name+suffix, // String title,
					  0.0,   // double r_min,
					  255.0, // double r_max,
					  0.0,   // double g_min,
					  255.0, // double g_max,
					  0.0,   // double b_min,
					  255.0, // double b_max,
					  0.0,   // double alpha_min,
					  1.0);  // double alpha_max)
			  return imp_rgba;
		  }

		  ImageStack stack = sdfa_instance.makeStack(
//				  rgb_in, // iclt_data,
				  rbg_in, // iclt_data,
				  width,  // (tilesX + 0) * clt_parameters.transform_size,
				  height, // (tilesY + 0) * clt_parameters.transform_size,
				  sliceNames,  // or use null to get chn-nn slice names
				  true); // replace NaN with 0.0
		  return linearStackToColor(
				  clt_parameters, // EyesisCorrectionParameters.CLTParameters         clt_parameters,
				  colorProcParameters, // EyesisCorrectionParameters.ColorProcParameters   colorProcParameters,
				  rgbParameters, // EyesisCorrectionParameters.RGBParameters         rgbParameters,
				  name,   // String name,
				  suffix, // String suffix, // such as disparity=...
				  toRGB, // boolean toRGB,
				  bpp16, // boolean bpp16, // 16-bit per channel color mode for result
				  saveShowIntermediate, // boolean saveShowIntermediate, // save/show if set globally
				  saveShowFinal, // boolean saveShowFinal,        // save/show result (color image?)
				  stack, // ImageStack stack,
				  alpha, // float [] alpha_pixels,
				  width, // int width, // int tilesX,
				  height, // int height, // int tilesY,
				  scaleExposure, // double scaleExposure,
				  debugLevel); //int debugLevel
	  }
	  // double data


	  public ImagePlus linearStackToColor( // USED in lwir
			  CLTParameters         clt_parameters,
			  ColorProcParameters   colorProcParameters,
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
		  ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
		  // convert to ImageStack of 3 slices
		  String [] sliceNames = {"red", "blue", "green"};
		  int green_index = 2;
		  double []   alpha = null; // (0..1.0)
		  double [][] rbg_in = {iclt_data[0],iclt_data[1],iclt_data[2]};
		  float [] alpha_pixels = null;
		  if (iclt_data.length > 3) {
			  alpha = iclt_data[3];
			  if (alpha != null){
				  alpha_pixels = new float [alpha.length];
				  for (int i = 0; i <alpha.length; i++){
					  alpha_pixels[i] = (float) alpha[i];
				  }
			  }
		  }
		  if (isLwir()) {
			  String [] rgb_titles =  {"red","green","blue"};
			  String [] rgba_titles = {"red","green","blue","alpha"};
			  String [] titles = (alpha == null) ? rgb_titles : rgba_titles;
			  int num_slices = (alpha == null) ? 3 : 4;
			  double mn = colorProcParameters.lwir_low;
			  double mx = colorProcParameters.lwir_high;
			  double [] cold_hot = getColdHot();
			  if (cold_hot != null) {
				  mn = cold_hot[0];
				  mx = cold_hot[1];
			  }

			  double offset = getLwirOffset();
			  if (!Double.isNaN(offset)) {
				  mn -=  offset;
				  mx -=  offset;
			  }

			  ThermalColor tc = new ThermalColor(
					  colorProcParameters.lwir_palette,
					  mn,
					  mx,
					  255.0);
			  double [][] rgba = new double [num_slices][];
			  for (int i = 0; i < 3; i++) rgba[i] = new double [iclt_data[green_index].length];
			  for (int i = 0; i < rbg_in[green_index].length; i++) {
				  if (i == 700) {
					  System.out.println("linearStackToColor(): i="+i);
				  }
				  double [] rgb = tc.getRGB(iclt_data[green_index][i]);
				  rgba[0][i] = rgb[0]; // red
				  rgba[1][i] = rgb[1]; // green
				  rgba[2][i] = rgb[2]; // blue
			  }
			  if (alpha != null) {
				  rgba[3] = alpha; // 0..1
			  }
			  ImageStack stack = sdfa_instance.makeStack(
					  rgba,       // iclt_data,
					  width,      // (tilesX + 0) * clt_parameters.transform_size,
					  height,     // (tilesY + 0) * clt_parameters.transform_size,
					  titles,     // or use null to get chn-nn slice names
					  true);      // replace NaN with 0.0
			  ImagePlus imp_rgba =  EyesisCorrections.convertRGBAFloatToRGBA32(
					  stack,   // ImageStack stackFloat, //r,g,b,a
					  //						name+"ARGB"+suffix, // String title,
					  name+suffix, // String title,
					  0.0,   // double r_min,
					  255.0, // double r_max,
					  0.0,   // double g_min,
					  255.0, // double g_max,
					  0.0,   // double b_min,
					  255.0, // double b_max,
					  0.0,   // double alpha_min,
					  1.0);  // double alpha_max)
			  return imp_rgba;
		  }

		  ImageStack stack = sdfa_instance.makeStack(
				  rbg_in, // iclt_data,
				  width,  // (tilesX + 0) * clt_parameters.transform_size,
				  height, // (tilesY + 0) * clt_parameters.transform_size,
				  sliceNames,  // or use null to get chn-nn slice names
				  true); // replace NaN with 0.0

		  return linearStackToColor(
				  clt_parameters, // EyesisCorrectionParameters.CLTParameters         clt_parameters,
				  // gamma should be 1.0 for LWIR
				  colorProcParameters, // EyesisCorrectionParameters.ColorProcParameters   colorProcParameters,
				  rgbParameters, // EyesisCorrectionParameters.RGBParameters         rgbParameters,
				  name,   // String name,
				  suffix, // String suffix, // such as disparity=...
				  toRGB, // boolean toRGB,
				  bpp16, // boolean bpp16, // 16-bit per channel color mode for result
				  saveShowIntermediate, // boolean saveShowIntermediate, // save/show if set globally
				  saveShowFinal, // boolean saveShowFinal,        // save/show result (color image?)
				  stack, // ImageStack stack,
				  alpha_pixels, // float [] alpha_pixels,
				  width, // int width, // int tilesX,
				  height, // int height, // int tilesY,
				  scaleExposure, // double scaleExposure,
				  debugLevel); //int debugLevel
	  }




	  // Convert a single value pixels to color (r,b,g) values to be processed instead of the normal colors


	  public ImagePlus linearStackToColor(
			  CLTParameters         clt_parameters,
			  ColorProcParameters   colorProcParameters,
			  EyesisCorrectionParameters.RGBParameters         rgbParameters,
			  String name,
			  String suffix, // such as disparity=...
			  boolean toRGB,
			  boolean bpp16, // 16-bit per channel color mode for result
			  boolean saveShowIntermediate, // save/show if set globally
			  boolean saveShowFinal,        // save/show result (color image?)
			  ImageStack stack,
			  float [] alpha_pixels,
			  int width, // int tilesX,
			  int height, // int tilesY,
			  double scaleExposure,
			  int debugLevel
			  )
	  {
//		  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		  if (debugLevel > -1) { // 0){
			  double [] chn_avg = {0.0,0.0,0.0};
			  float [] pixels;
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
			  return null;// not used in lwir
		  }
		  if (debugLevel > 1) System.out.println("before colors.2");
//		  if (debugLevel > -1) System.out.println("before colors.2");
		  if (saveShowIntermediate && (debugLevel > 1)){
//		  if (saveShowIntermediate && (debugLevel > -1)){
			  ImagePlus imp_dbg=new ImagePlus(name+sAux()+"-preColors",stack);
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
			  titleFull=name+sAux()+"-RGB-float"+suffix;
			  //Trim stack to just first 3 slices
			  if (saveShowIntermediate && (debugLevel > 1)){ // 2){
				  ImagePlus imp_dbg=new ImagePlus("YPrPbToRGB",stack);
				  eyesisCorrections.saveAndShow(
						  imp_dbg,
						  this.correctionsParameters);
			  }
			  while (stack.getSize() > 3) stack.deleteLastSlice();
			  if (debugLevel > 1) System.out.println("Trimming color stack");
		  } else {// not used in lwir
			  titleFull=name+sAux()+"-YPrPb"+suffix;
			  if (debugLevel > 1) System.out.println("Using full stack, including YPbPr");
		  }

		  if (alpha_pixels != null){
			  stack.addSlice("alpha",alpha_pixels);
		  }

		  ImagePlus result= new ImagePlus(titleFull, stack);
		  if (debugLevel> 1){
			  result.show();
		  }

		  if (!toRGB && !this.correctionsParameters.jpeg){ // toRGB set for equirectangular// not used in lwir
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

		  titleFull=name+sAux()+"-RGB48"+suffix;
		  result= new ImagePlus(titleFull, stack);
		  result.updateAndDraw();
		  if (debugLevel > 1) {
//		  if (debugLevel > -1) {
			  System.out.println("result.updateAndDraw(), "+titleFull+"-RGB48");
			  result.show();
		  }

		  CompositeImage compositeImage=eyesisCorrections.convertToComposite(result);

		  if (!this.correctionsParameters.jpeg && bpp16){ // RGB48 was the end result // not used in lwir
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



	  public void apply_fine_corr( // not used in lwir
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
	  public void show_fine_corr() // not used in lwir
	  {
		  show_fine_corr("");
	  }
	  public void show_fine_corr(String prefix) // not used in lwir
	  {
		  show_fine_corr( this.fine_corr, prefix);
	  }

	  public void show_fine_corr( // not used in lwir
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


	  public void reset_fine_corr() // not used in lwir
	  {
		  this.fine_corr = new double [4][2][6]; // reset all coefficients to 0
	  }

	  public void showExtrinsicCorr(String name) // not used in lwir
	  {
		  System.out.println("Extrinsic corrections "+name);
		  if (geometryCorrection == null){
			  System.out.println("are not set, will be:");
			  System.out.println(new GeometryCorrection(this.extrinsic_vect).getCorrVector().toString());
		  } else {
			  System.out.println(geometryCorrection.getCorrVector().toString());
		  }
	  }

	  public boolean editExtrinsicCorr() // not used in lwir
	  {
		  if (geometryCorrection == null){
			  System.out.println("are not set, will be:");
			  return new GeometryCorrection(this.extrinsic_vect).getCorrVector().editVector(); // editIMU();
		  } else {
			  return geometryCorrection.getCorrVector().editVector(); // .editIMU();
		  }
	  }



	  public boolean editRig() // not used in lwir
	  {
		  if (!is_aux) {
			  System.out.println("Rig offsets can only be edited for the auxiliary camera, not for the main one");
			  return false;
		  }
//		  GeometryCorrection gc = this.geometryCorrection;
		  if (this.geometryCorrection == null){
			  System.out.println("geometryCorrection is not set, creating one");
			  this.geometryCorrection = new GeometryCorrection(this.extrinsic_vect);
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


	  public void resetExtrinsicCorr( // not used in lwir
			  CLTParameters           clt_parameters)
	  {
//		  this.extrinsic_vect = new double [GeometryCorrection.CORR_NAMES.length];
		  this.extrinsic_vect = null;
		  if (geometryCorrection != null){
			  geometryCorrection.resetCorrVector();
		  }
		  if (clt_parameters.fine_corr_apply){
			  clt_parameters.fine_corr_ignore = false;
		  }
	  }

	  public void cltDisparityScans( // not used in lwir
			  CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
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
				  int [] channels= fileChannelToSensorChannels(correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile]));
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
				  int [] channels= fileChannelToSensorChannels(correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile]));
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
			  setFiles.get(setNames.indexOf(setName)).add(nFile); // .add(new Integer(nFile));
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
			  this.geometryCorrection.woi_tops = new int [channelFiles.length];
			  boolean [][] saturation_imp = (clt_parameters.sat_level > 0.0)? new boolean[channelFiles.length][] : null;
			  double [] scaleExposures = new double[channelFiles.length];
			  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
				  int nFile=channelFiles[srcChannel];
				  imp_srcs[srcChannel]=null;
				  if (nFile >=0){
					  imp_srcs[srcChannel] = eyesisCorrections.getJp4Tiff(sourceFiles[nFile], this.geometryCorrection.woi_tops);

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

			  System.out.println("Temporarily applying scaleExposures[] here - 2" );
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
					  colorProcParameters,
					  channelGainParameters,
					  rgbParameters,
					  scaleExposures,
					  threadsMax,  // maximal number of threads to launch
					  updateStatus,
					  debugLevel);
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


	  public ImagePlus [] cltDisparityScan( // not used in lwir
			  ImagePlus [] imp_quad, // should have properties "name"(base for saving results), "channel","path"
			  boolean [][] saturation_imp, // (near) saturated pixels or null
			  CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
//			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
//			  int              convolveFFTSize, // 128 - fft size, kernel size should be size/2
			  double []	       scaleExposures, // probably not needed here
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel){
		  ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?

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
					  null, // no margins, no oversample
					  this.is_mono);
		  }
		  // =================
		  ImageDtt image_dtt = new ImageDtt(
				  clt_parameters.transform_size,
				  isMonochrome(),
				  isLwir(),
				  clt_parameters.getScaleStrength(isAux()));
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
					  clt_parameters.getFatZero(isMonochrome()),      // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
					  clt_parameters.corr_sym,
					  clt_parameters.corr_offset,
					  clt_parameters.corr_red,
					  clt_parameters.corr_blue,
					  clt_parameters.getCorrSigma(image_dtt.isMonochrome()),
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
///					  image_dtt.transform_size,
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
				  int num_tiles = disp_strength[0].length;
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
						  clt_parameters.getLyPerQuad(num_tiles), // final int        min_per_quadrant, // minimal tiles per quadrant (not counting the worst) tp proceed
						  clt_parameters.getLyInf(num_tiles),	  // final int        min_inf,          // minimal number of tiles at infinity to proceed
						  clt_parameters.getLyInfScale(num_tiles),// final int        min_inf_to_scale, // minimal number of tiles at infinity to apply weight scaling
						  clt_parameters.ly_right_left,      // false // equalize weights of right/left FoV (use with horizon in both halves and gross infinity correction)
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

			  ImagePlus imp_stack = new ImagePlus( name+sAux()+"-DISP_MAPS", array_stack);
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
			  ImagePlus imp_stack_trends = new ImagePlus( name+sAux()+"-DISP_TRENDS", trends_stack);
			  imp_stack_trends.getProcessor().resetMinAndMax();
			  imp_stack_trends.updateAndDraw();
			  eyesisCorrections.saveAndShow(imp_stack_trends, this.correctionsParameters);
		  }
		  return results;
	  }

	  public void process_infinity_corr( //from existing image // not used in lwir
			  CLTParameters clt_parameters,
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


	  public void processLazyEye( // not used in lwir
			  boolean dry_run,
			  CLTParameters clt_parameters,
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

			int num_tiles = disp_strength[0].length;

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
		    		clt_parameters.getLyPerQuad(num_tiles), // final int        min_per_quadrant, // minimal tiles per quadrant (not counting the worst) tp proceed
		    		clt_parameters.getLyInf(num_tiles),	    // final int        min_inf,          // minimal number of tiles at infinity to proceed
		    		clt_parameters.getLyInfScale(num_tiles),// final int        min_inf_to_scale, // minimal number of tiles at infinity to apply weight scaling
		    		clt_parameters.ly_right_left,      // false // equalize weights of right/left FoV (use with horizon in both halves and gross infinity correction)
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


	  public double [][] process_disparity_scan( // not used in lwir
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

	  public void showCLTPlanes( // not used in lwir
			  CLTParameters           clt_parameters,
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

	  public double [][]  assignCLTPlanes( // not used in lwir
			  CLTParameters           clt_parameters,
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
			setPassAvgRBGA(                      // get image from a single pass, return relative path for x3d // USED in lwir
					clt_parameters,                           // CLTParameters           clt_parameters,
					tp.clt_3d_passes.size() - 1, // int        scanIndex,
					threadsMax,                               // int        threadsMax,  // maximal number of threads to launch
					updateStatus,                             // boolean    updateStatus,
					debugLevel);                         // int        debugLevel)
		  	
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


	  public void out3d_old( // not used in lwir
			  CLTParameters           clt_parameters,
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



	  public void processCLTQuads3d( // not used in lwir
			  boolean                                              adjust_extrinsics,
			  boolean                                              adjust_poly,
			  TwoQuadCLT                                           twoQuadCLT, //maybe null in no-rig mode, otherwise may contain rig measurements to be used as infinity ground truth
			  CLTParameters                                        clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters         debayerParameters,
			  ColorProcParameters                                  colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters             channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
			  final int                                            threadsMax,  // maximal number of threads to launch
			  final boolean                                        updateStatus,
			  final int                                            debugLevel)
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
//			  boolean [][] 
					  saturation_imp = (clt_parameters.sat_level > 0.0)? new boolean[channelFiles.length][] : null;
			  double [] scaleExposures = new double[channelFiles.length];

			  ImagePlus [] imp_srcs = conditionImageSet(
					  clt_parameters,             // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					  colorProcParameters,        // ColorProcParameters                       colorProcParameters,
					  sourceFiles,                // String []                                 sourceFiles,
					  set_channels[nSet].name(),  // String                                    set_name,
					  referenceExposures,         // double []                                 referenceExposures,
					  channelFiles,               // int []                                    channelFiles,
					  scaleExposures,   //output  // double [] scaleExposures
					  saturation_imp,   //output  // boolean [][]                              saturation_imp,
					  threadsMax,                 // int                                       threadsMax,
					  debugLevel); // int                                       debugLevel);
			  boolean use_rig = (twoQuadCLT != null) && (twoQuadCLT.getBiScan(0) != null);
			  if (!adjust_extrinsics || !use_rig) {
				  // Difficult to fix: adjust extrinsics for aux - when it is adjusted alone, it will not match tiles to those of a rig!
				  // can use only far tiles with small gradients?

				  // once per quad here

// need to replace for low-res?
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
// Add here composite scans and show FG and BG images
				  // adjust extrinsics here


				  ArrayList<CLTPass3d> combo_pass_list = tp.compositeScan(
						  2, // just FG and BG
						  tp.clt_3d_passes, // final ArrayList <CLTPass3d> passes,
						  0, // bg_pass, //  final int                   firstPass,
						  tp.clt_3d_passes.size(),                          //  final int                   lastPassPlus1,
						  tp.getTrustedCorrelation(),                       // 	 final double                trustedCorrelation,
						  tp.getMaxOverexposure(),                          //  final double                max_overexposure,
						  0.0, // clt_parameters.bgnd_range,                //	 final double                disp_far,   // limit results to the disparity range
						  clt_parameters.grow_disp_max,                     // final double                disp_near,
						  clt_parameters.combine_min_strength,              // final double                minStrength,
						  clt_parameters.combine_min_hor,                   // final double                minStrengthHor,
						  clt_parameters.combine_min_vert,                  // final double                minStrengthVert,
						  false,                                            // final boolean               no_weak,
						  false,                                            // final boolean               use_last,   //
						  // TODO: when useCombo - pay attention to borders (disregard)
						  false,                                            // final boolean               usePoly)  // use polynomial method to find max), valid if useCombo == false
						  true,                                             // 	 final boolean               copyDebug)
						  debugLevel);

				  tp.ShowScansSFB(
						  combo_pass_list, // ArrayList<CLTPass3d> scans, // list of composite scans
						  this.image_name+sAux()+"-SFB0"); // String               title);

			  }
			  if (adjust_extrinsics) {
				  // temporarily
				  if (ds_from_main != null) {
					  System.out.println("Adjust AUX extrinsics using main camera measurements");
					  extrinsicsCLTfromGT(
//							  twoQuadCLT,   // TwoQuadCLT       twoQuadCLT, //maybe null in no-rig mode, otherwise may contain rig measurements to be used as infinity ground truth
							  null,
							  ds_from_main, // gt_disp_strength,
							  clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
							  adjust_poly,
							  threadsMax,  //final int        threadsMax,  // maximal number of threads to launch
							  updateStatus,// final boolean    updateStatus,
							  debugLevel + 2); // final int        debugLevel)

				  } else  if (use_rig) {
					  System.out.println("Adjust extrinsics using rig data here");
					  double [][] gt_disp_strength = getRigDSFromTwoQuadCL(
							  twoQuadCLT, //maybe null in no-rig mode, otherwise may contain rig measurements to be used as infinity ground truth
							  clt_parameters,
							  debugLevel + 2); // final int        debugLevel)
					  GeometryCorrection geometryCorrection_main = null;
					  if (geometryCorrection.getRotMatrix(true) != null) {
						  geometryCorrection_main = twoQuadCLT.quadCLT_main.getGeometryCorrection();
					  }
					  extrinsicsCLTfromGT(
//							  twoQuadCLT,   // TwoQuadCLT       twoQuadCLT, //maybe null in no-rig mode, otherwise may contain rig measurements to be used as infinity ground truth
							  geometryCorrection_main,
							  gt_disp_strength,
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

			  ArrayList<CLTPass3d> combo_pass_list = tp.compositeScan(
					  2, // just FG and BG
					  tp.clt_3d_passes, // final ArrayList <CLTPass3d> passes,
					  0, // bg_pass, //  final int                   firstPass,
					  tp.clt_3d_passes.size(),                          //  final int                   lastPassPlus1,
					  tp.getTrustedCorrelation(),                       // 	 final double                trustedCorrelation,
					  tp.getMaxOverexposure(),                          //  final double                max_overexposure,
					  0.0, // clt_parameters.bgnd_range,                //	 final double                disp_far,   // limit results to the disparity range
					  clt_parameters.grow_disp_max,                     // final double                disp_near,
					  clt_parameters.combine_min_strength,              // final double                minStrength,
					  clt_parameters.combine_min_hor,                   // final double                minStrengthHor,
					  clt_parameters.combine_min_vert,                  // final double                minStrengthVert,
					  false,                                            // final boolean               no_weak,
					  false,                                            // final boolean               use_last,   //
					  // TODO: when useCombo - pay attention to borders (disregard)
					  false,                                            // final boolean               usePoly)  // use polynomial method to find max), valid if useCombo == false
					  true,                                             // 	 final boolean               copyDebug)
					  debugLevel);

			  tp.ShowScansSFB(
					  combo_pass_list, // ArrayList<CLTPass3d> scans, // list of composite scans
					  this.image_name+sAux()+"-SFB1"); // String               title);


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

	public double [][] depthMapMainToAux(// USED in lwir
			double [][]        ds,
			GeometryCorrection geometryCorrection_main,
			GeometryCorrection geometryCorrection_aux,
			CLTParameters      clt_Parameters,
//			double             min_strength,
//			boolean            use_wnd,
			boolean            split_fg_bg,
//			double             split_fbg_rms,
			boolean            for_adjust, // for LY adjustment: only keep d,s and remove samples with high variations
			int                debug_level
			){
		class DS{// USED in lwir
			double disparity;  // gt disparity
			double strength;   // gt strength
			int tx;            // gt tile x
			int ty;            // gt tile x
			double fx;         // fractional aux tile X (0.0..1.0) for optional window
			double fy;         // fractional aux tile Y (0.0..1.0) for optional window
			DS (double disparity, double strength, int tx, int ty, double fx, double fy){
				this.disparity = disparity;
				this.strength =  strength;
				this.tx =        tx;
				this.ty =        ty;
				this.fx =        fx;
				this.fy =        fy;

			}
			@Override
			public String toString() { // not used in lwir
				return String.format("Disparity (str) = % 6f (%5f), tx=%d ty=%d fx=%5f fy=%5f\n", disparity, strength,tx,ty,fx,fy);
			}
		}
		int tile_size =  clt_Parameters.transform_size;
		int [] wh_main = geometryCorrection_main.getSensorWH();
		int [] wh_aux =  geometryCorrection_aux.getSensorWH();
		int tilesX_main = wh_main[0] / tile_size;
		int tilesY_main = wh_main[1] / tile_size;
		int tilesX_aux = wh_aux[0] / tile_size;
		int tilesY_aux = wh_aux[1] / tile_size;

		ArrayList<ArrayList<DS>> ds_list = new ArrayList<ArrayList<DS>>();
		for (int nt = 0; nt < tilesX_aux * tilesY_aux; nt++) {
			ds_list.add(new ArrayList<DS>());
		}
		for (int ty = 0; ty < tilesY_main; ty++) {
			double centerY = ty * tile_size + tile_size/2;
			for (int tx = 0; tx < tilesX_main; tx++) {
				int nt = ty*tilesX_main + tx;
				double centerX = tx * tile_size + tile_size/2;
				double disparity = ds[0][nt];
				double strength =  ds[1][nt];
				if ((strength >= clt_Parameters.ly_gt_strength) && !Double.isNaN(disparity)) {
					double [] dpxpy_aux =  geometryCorrection_aux.getFromOther(
							geometryCorrection_main, // GeometryCorrection other_gc,
							centerX,                 // double other_px,
							centerY,                 // double other_py,
							disparity);              // double other_disparity)
					double fx = dpxpy_aux[1]/tile_size;
					double fy = dpxpy_aux[2]/tile_size;
					int tx_aux = (int) Math.floor(fx);
					int ty_aux = (int) Math.floor(fy);
					fx -= tx_aux;
					fy -= ty_aux;
					if ((ty_aux >= 0) && (ty_aux < tilesY_aux) && (tx_aux >= 0) && (tx_aux < tilesX_aux)) {
						int nt_aux = ty_aux * tilesX_aux + tx_aux;
						ds_list.get(nt_aux).add(new DS(dpxpy_aux[0], strength, tx, ty, fx, fy));
					}
				}
			}
		}

		// simple average (ignoring below minimal)
		int num_slices = split_fg_bg? FGBG_TITLES_AUX.length:FGBG_TITLES_ADJ.length;
		double [][] ds_aux_avg = new double [num_slices][tilesX_aux * tilesY_aux];
		for (int ty = 0; ty < tilesY_aux; ty++) {
			for (int tx = 0; tx < tilesX_aux; tx++) {
//				if ((ty == 3) && (tx == 12)) {
//					System.out.println("tx = "+tx+", ty = "+ty);
//				}
				int nt = ty * tilesX_aux + tx;
				ds_aux_avg[FGBG_DISPARITY][nt] = Double.NaN;
				ds_aux_avg[FGBG_STRENGTH][nt] = 0.0;
				if(ds_list.get(nt).isEmpty()) continue;
	    		Collections.sort(ds_list.get(nt), new Comparator<DS>() {
	    		    @Override
	    		    public int compare(DS lhs, DS rhs) { // ascending
	    		        return rhs.disparity > lhs.disparity  ? -1 : (rhs.disparity  < lhs.disparity ) ? 1 : 0;
	    		    }
	    		});

				double sw = 0.0, swd = 0.0, swd2 = 0.0;
	    		for (DS dsi: ds_list.get(nt)) {
	    			double w = dsi.strength;
	    			if (clt_Parameters.ly_gt_use_wnd) {
	    				w *= Math.sin(Math.PI * dsi.fx) * Math.sin(Math.PI * dsi.fy);
	    			}
	    			sw +=  w;
	    			double wd = w * dsi.disparity;
	    			swd += wd;
	    			swd2 += wd * dsi.disparity;

	    		}
	    		ds_aux_avg[FGBG_DISPARITY][nt] = swd/sw;
	    		ds_aux_avg[FGBG_STRENGTH][nt] = sw/ds_list.get(nt).size();
	    		double rms = Math.sqrt( (swd2 * sw - swd * swd) / (sw * sw));
    			if (for_adjust && (rms >= clt_Parameters.ly_gt_rms)) { // remove ambiguous tiles
    	    		ds_aux_avg[FGBG_DISPARITY][nt] = Double.NaN;
    	    		ds_aux_avg[FGBG_STRENGTH][nt] = 0;

    			}
	    		if (split_fg_bg) {

	    			ds_aux_avg[FGBG_RMS ][nt] =      rms;
	    			ds_aux_avg[FGBG_RMS_SPLIT][nt] = ds_aux_avg[2][nt]; // rms
	    			ds_aux_avg[FGBG_FG_DISP][nt] =   ds_aux_avg[0][nt]; // fg disp
	    			ds_aux_avg[FGBG_FG_STR][nt] =    ds_aux_avg[1][nt]; // fg strength
	    			ds_aux_avg[FGBG_BG_DISP][nt] =   ds_aux_avg[0][nt]; // bg disp
	    			ds_aux_avg[FGBG_BG_STR][nt] =    ds_aux_avg[1][nt]; // bg strength
	    			if (rms >= clt_Parameters.ly_gt_rms) {
	    				// splitting while minimizing sum of 2 squared errors
	    	    		double [][] swfb =  new double [2][ds_list.get(nt).size() -1];
	    	    		double [][] swdfb = new double [2][ds_list.get(nt).size() -1];
	    	    		double []   s2fb =  new double [ds_list.get(nt).size() -1];
	    	    		for (int n = 0; n < s2fb.length; n++) { // split position
	    	    			double [] s2 = new double[2];
	    	    			for (int i = 0; i <= s2fb.length; i++) {
	    	    				int fg = (i > n)? 1 : 0; // 0 - bg, 1 - fg
	    	    				DS dsi = ds_list.get(nt).get(i);
	    		    			double w = dsi.strength;
	    		    			if (clt_Parameters.ly_gt_use_wnd) {
	    		    				w *= Math.sin(Math.PI * dsi.fx) * Math.sin(Math.PI * dsi.fy);
	    		    			}
	    		    			swfb[fg][n] +=  w;
	    		    			double wd =      w * dsi.disparity;
	    		    			swdfb[fg][n] += wd;
	    		    			s2[fg] +=        wd * dsi.disparity;
	    	    			}
	    	    			s2fb[n] =  ((s2[0] * swfb[0][n] - swdfb[0][n] * swdfb[0][n]) / swfb[0][n] +
	    	    					(s2[1] * swfb[1][n] - swdfb[1][n] * swdfb[1][n]) / swfb[1][n]) / (swfb[0][n] + swfb[1][n]);
	    	    		}
		    			// now find the n with lowest s2fb and use it to split fg/bg. Could be done in a single pass, but with saved arrays
		    			// it is easier to verify
		    			int nsplit = 0;
		    			for (int i = 1; i < s2fb.length; i++) if (s2fb[i] < s2fb[nsplit]) {
		    				nsplit = i;
		    			}
		    			ds_aux_avg[FGBG_RMS_SPLIT][nt] = s2fb[nsplit]; // rms split
		    			ds_aux_avg[FGBG_FG_DISP][nt] =   swdfb[1][nsplit] / swfb[1][nsplit] ;      // fg disp
		    			ds_aux_avg[FGBG_FG_STR][nt] =    swfb[1][nsplit]/ (s2fb.length - nsplit) ; // fg strength

		    			ds_aux_avg[FGBG_BG_DISP][nt] =   swdfb[0][nsplit] / swfb[0][nsplit] ;      // bg disp
		    			ds_aux_avg[FGBG_BG_STR][nt] =    swfb[0][nsplit]/ (nsplit + 1) ;           // bg strength
	    			}
	    		}
			}
		}
		return ds_aux_avg;

	}

	  public boolean preExpandCLTQuad3d( // USED in lwir
			  ImagePlus []                                     imp_quad, // should have properties "name"(base for saving results), "channel","path"
			  boolean [][]                                     saturation_imp,   // (near) saturated pixels or null
			  CLTParameters         clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  ColorProcParameters   colorProcParameters,
			  EyesisCorrectionParameters.RGBParameters         rgbParameters,
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel)
	  {
		  boolean no_macro = isLwir(); // make it a separate configurable parameter?

		  this.startStepTime=System.nanoTime();
		  final boolean    show_init_refine = clt_parameters.show_init_refine;

		  //max_expand
		  String name = (String) imp_quad[0].getProperty("name");
		  // should create data for the macro! (diff, rgb)
		  CLTPass3d bgnd_data = CLTBackgroundMeas( // measure background - both CPU and GPU (remove textures from GPU)
				  clt_parameters,
				  threadsMax,  // maximal number of threads to launch
				  updateStatus,
				  debugLevel);
		  tp.clt_3d_passes.add(bgnd_data);
		  //    	  if (show_init_refine)
		  if ((debugLevel > -2) && clt_parameters.show_first_bg) {
//		  if ((debugLevel > -3) && clt_parameters.show_first_bg) {
			  tp.showScan(
					  tp.clt_3d_passes.get(0), // CLTPass3d   scan,
					  "bgnd_data-"+tp.clt_3d_passes.size());
		  }

		  //TODO: Move away from here?
		  boolean no_image_save = false; // Restore? true;

		  boolean [] bgmask = getBackgroundImageMasks(
    			  clt_parameters,
    			  name,               // .getTitle(), //String name=(String) imp_src.getProperty("name");
    			  ImageDtt.DISPARITY_INDEX_CM, // index of disparity value in disparity_map == 2 (0,2 or 4)
    			  threadsMax,  // maximal number of threads to launch
    			  updateStatus,
    			  debugLevel);
		  
    	  ImagePlus imp_bgnd_int = getBackgroundImage(
    			  bgmask, // boolean []                                bgnd_tiles, 
    			  clt_parameters,
    			  colorProcParameters,
    			  rgbParameters,
    			  name,               // .getTitle(), //String name=(String) imp_src.getProperty("name");
    			  ImageDtt.DISPARITY_INDEX_CM, // index of disparity value in disparity_map == 2 (0,2 or 4)
    			  threadsMax,  // maximal number of threads to launch
    			  updateStatus,
    			  debugLevel);
		  
//    	  imp_bgnd_int.show();
		  
//		  if (debugLevel > -100) {
//			  return null;
//		  }
		  // resize for backdrop here!

    	  ImagePlus imp_bgnd = finalizeBackgroundImage( // USED in lwir
    			  imp_bgnd_int,   //  ImagePlus imp_texture_bgnd,
    			  no_image_save,  // boolean    no_image_save,
    			  clt_parameters, // CLTParameters           clt_parameters,
    			  name,           // String     name,
    			  debugLevel);    // int        debugLevel)

//    	  imp_bgnd.show();
    	  
    	  bgnd_data.texture = (imp_bgnd == null)? null: ( imp_bgnd.getTitle()+ (clt_parameters.black_back? ".jpeg" : ".png"));

    	  // create x3d file
    	  X3dOutput x3dOutput = new X3dOutput(
    			  clt_parameters,
    			  correctionsParameters,
    			  geometryCorrection,
    			  tp.clt_3d_passes);

		  x3dOutput.generateBackground(clt_parameters.infinityDistance <= 0.0); // needs just first (background) scan



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
    	  ArrayList <CLTPass3d> new_meas = null; // filled either from macro correlation, or just from plain disparity scan
    	  if (no_macro) {
    		  new_meas = prepareDisparityScan(
    				  clt_parameters.disp_scan_start, // double scan_start,
    				  clt_parameters.disp_scan_step, // 	double scan_step,
    				  clt_parameters.disp_scan_count); // 	int    scan_count)
    	  } else {
    		  MacroCorrelation mc = new MacroCorrelation(
    				  tp,
    				  clt_parameters.mc_disp8_trust,
    				  clt_parameters.mc_weight_var,   // final double     weight_var,   // weight of variance data (old, detects thin wires?)
    				  clt_parameters.mc_weight_Y,     // final double     weight_Y,     // weight of average intensity
    				  clt_parameters.mc_weight_RBmG   // final double     weight_RBmG,  // weight of average color difference (0.5*(R+B)-G), shoukld be ~5*weight_Y
    				  );


    		  double [][][] input_data = mc.CLTMacroSetData( // perform single pass according to prepared tiles operations and disparity
    				  bgnd_data);           // final CLTPass3d      src_scan, // results of the normal correlations (now expecting infinity)

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

    		  //    	  ArrayList <CLTPass3d>
    		  new_meas = mc.prepareMeasurementsFromMacro(
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
    	  }  // if (no_macro) {} else

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

    			  // Just debugging
    			  if (clt_parameters.gpu_debug_accum) {
    				  CLTMeasureCorrTesting( // perform single pass according to prepared tiles operations and disparity
    						  clt_parameters,
    						  0,
    						  false, // true, // final boolean     save_textures,
    						  threadsMax,  // maximal number of threads to launch
    						  updateStatus,
    						  debugLevel);

    				  CLTMeasureCorrTesting( // perform single pass according to prepared tiles operations and disparity
    						  clt_parameters,
    						  refine_pass,
    						  false, // true, // final boolean     save_textures,
    						  threadsMax,  // maximal number of threads to launch
    						  updateStatus,
    						  debugLevel);
    			  }
      			  // End of just debugging
    			  
    			  CLTMeasureCorr( // perform single pass according to prepared tiles operations and disparity
    					  clt_parameters,
    					  refine_pass,
    					  false, // true, // final boolean     save_textures,
    					  threadsMax,  // maximal number of threads to launch
    					  updateStatus,
    					  debugLevel);
    			  
    			  
    			  if (debugLevel > -3){
    				  System.out.println("CLTMeasure("+refine_pass+")-*");
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
    						  false, // final boolean               useP}oly)  // use polynomial method to find max), valid if useCombo == false
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
    		  CLTMeasureCorr( // perform single pass according to prepared tiles operations and disparity
//					  image_data, // first index - number of image in a quad
//					  saturation_imp, //final boolean [][]  saturation_imp, // (near) saturated pixels or null
					  clt_parameters,
					  tp.clt_3d_passes.size() -1, // refine_pass, - just added from macro?
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

			  CLTMeasureCorr( // perform single pass according to prepared tiles operations and disparity
//					  image_data, // first index - number of image in a quad
//					  saturation_imp, //final boolean [][]  saturation_imp, // (near) saturated pixels or null
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

			  CLTPass3d combo_pass = tp.compositeScan(
					  tp.clt_3d_passes, // final ArrayList <CLTPass3d> passes,
					  bg_pass, //  final int                   firstPass,
					  tp.clt_3d_passes.size(),                          //  final int                   lastPassPlus1,
					  tp.getTrustedCorrelation(),                       // 	 final double                trustedCorrelation,
					  tp.getMaxOverexposure(),                          //  final double                max_overexposure,
					  0.0, // clt_parameters.bgnd_range,                //	 final double                disp_far,   // limit results to the disparity range
					  clt_parameters.grow_disp_max,                     // final double                disp_near,
					  clt_parameters.combine_min_strength,              // final double                minStrength,
					  clt_parameters.combine_min_hor,                   // final double                minStrengthHor,
					  clt_parameters.combine_min_vert,                  // final double                minStrengthVert,
					  false,                                            // final boolean               no_weak,
					  false,                                            // final boolean               use_last,   //
					  // TODO: when useCombo - pay attention to borders (disregard)
					  false,                                            // final boolean               usePoly)  // use polynomial method to find max), valid if useCombo == false
					  true,                                             // 	 final boolean               copyDebug)
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

	  ArrayList <CLTPass3d> prepareDisparityScan( // USED in lwir
			  	double scan_start,
			  	double scan_step,
			  	int    scan_count){
		  ArrayList <CLTPass3d> measurements = new ArrayList <CLTPass3d>();
		  for (int si = 0; si < scan_count; si++ ) {
			  double disparity = scan_start + scan_step * si;
			  CLTPass3d pass = new CLTPass3d(tp, 0 );
			  int op = ImageDtt.setImgMask(0, 0xf);
			  op =     ImageDtt.setPairMask(op,0xf);
			  op =     ImageDtt.setForcedDisparity(op,true);
			  pass.disparity = tp.setSameDisparity(disparity);
			  pass.tile_op = tp.setSameTileOp(op);
			  measurements.add(pass);
		  }
		  return measurements;
	  }
	  
	  /**
	   *
	   * @param clt_parameters
	   * @param adjust_poly
	   * @param threadsMax
	   * @param updateStatus
	   * @param debugLevel
	   * @return true on success, false - on failure
	   */
	  public boolean extrinsicsCLT( // USED in lwir
			  CLTParameters           clt_parameters,
			  boolean 		   adjust_poly,
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel)
	  {
		  final boolean    batch_mode = clt_parameters.batch_run;
		  int debugLevelInner =  batch_mode ? -5: debugLevel;
		  boolean update_disp_from_latest = clt_parameters.lym_update_disp ; // true;
		  int max_tries =                   clt_parameters.lym_iter; // 25;
		  double min_sym_update =           clt_parameters.getLymChange(is_aux); //  4e-6; // stop iterations if no angle changes more than this
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
// combo_scan: normStrength - junk. Is it used?
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
					  (bg_str[nTile] > clt_parameters.fcorr_inf_strength) && // 0.13
					  ((bg_overexp == null) || (bg_overexp[nTile] < clt_parameters.lym_overexp)) //1e-4
					  ){
				  bg_use[nTile] = true;
			  }
		  }
		  if (true) {
			  String [] dbg_titles = {"fdisp", "fstr", "disp", "str", "overexp","sel","use"};
			  double [][] ddd = {filtered_bgnd_disp_strength[0],filtered_bgnd_disp_strength[1],null,bg_str,bg_overexp, null, null};
			  ddd[5] = new double [bg_sel.length];
			  ddd[6] = new double [bg_sel.length];
			  for (int  nTile = 0 ; nTile < bg_use.length; nTile++) {
				  ddd[5][nTile] = bg_sel[nTile]?1.0:0.0;
				  ddd[6][nTile] = bg_use[nTile]?1.0:0.0;
			  }
			  (new ShowDoubleFloatArrays()).showArrays(
				  ddd,
				  tp.getTilesX(),
				  tp.getTilesY(),
				  true,
				  "filtered_bgnd_disp_strength",dbg_titles);
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
/*
		  CLTMeasure( // perform single pass according to prepared tiles operations and disparity **** CPU?
				  clt_parameters,
				  combo_scan,
				  false, // final boolean     save_textures,
				  true,  // final boolean     save_corr,
				  null,           // final double [][] mismatch,    // null or double [12][]
				  tp.threadsMax,  // maximal number of threads to launch
				  false, // updateStatus,
				  debugLevelInner - 1);

 */
		  CLTMeasureCorr( // perform single pass according to prepared tiles operations and disparity **** CPU?
				  clt_parameters,
				  combo_scan,
				  false, // final boolean     save_textures,
				  tp.threadsMax,  // maximal number of threads to launch
				  false, // updateStatus,
				  debugLevelInner - 1);
		  if (!batch_mode && clt_parameters.show_extrinsic && (debugLevel >-3)) {
			  tp.showScan(
					  tp.clt_3d_passes.get(bg_scan),   // CLTPass3d   scan, badly filtered?
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
		  if (true) {
			  String [] dbg_titles = {"fdisp", "fstr", "disp", "str", "overexp","sel","use"};
			  double [][] ddd = {filtered_combo_scand_isp_strength[0],filtered_combo_scand_isp_strength[1],combo_disp,combo_str,combo_overexp, null, null};
			  ddd[5] = new double [combo_use.length];
			  ddd[6] = new double [combo_use.length];
			  for (int  nTile = 0 ; nTile < combo_use.length; nTile++) {
				  //ddd[5][nTile] = bg_sel[nTile]?1.0:0.0;
				  ddd[6][nTile] = combo_use[nTile]?1.0:0.0;
			  }
			  (new ShowDoubleFloatArrays()).showArrays(
					  ddd,
					  tp.getTilesX(),
					  tp.getTilesY(),
					  true,
					  "filtered_combo_scand_isp_strength",dbg_titles);
		  }
		  
		  
		  int num_combo1 = tp.clt_3d_passes.get(combo_scan).setTileOpDisparity( // GPU ==0 !
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
			  double [][]dbg_img = { // bg_use - all 0? (never assigned)?
					  filtered_bgnd_disp_strength[0],
					  filtered_bgnd_disp_strength[1],
					  filtered_combo_scand_isp_strength[0],
					  filtered_combo_scand_isp_strength[1],
					  dbg_bg_sel,
					  dbg_bg_use, // too few
					  dbg_combo_use};
			  (new ShowDoubleFloatArrays()).showArrays(dbg_img,  tp.getTilesX(), tp.getTilesY(), true, "extrinsics_bgnd_combo",titles);
		  }
		  AlignmentCorrection ac = null;
		  if (!clt_parameters.ly_lma_ers ) {
			  ac = new AlignmentCorrection(this);
		  }
		  // iteration steps
		  // if (!batch_mode && clt_parameters.show_extrinsic && (debugLevel >-1)) {
		  if (clt_parameters.show_extrinsic && (debugLevel > -1)) { // temporary
			  tp.showScan(
					  tp.clt_3d_passes.get(bg_scan),   // CLTPass3d   scan,
					  "bg_scan_post"); //String title)
			  tp.showScan(
					  tp.clt_3d_passes.get(combo_scan),   // CLTPass3d   scan,
					  "combo_scan-"+combo_scan+"_post"); //String title)
		  }


		  double comp_diff = min_sym_update + 1; // (> min_sym_update)

		  for (int num_iter = 0; num_iter < max_tries; num_iter++){
			  if (update_disp_from_latest) {
				  tp.clt_3d_passes.get(combo_scan).updateDisparity();
			  }
			  if (clt_parameters.ly_lma_ers) {
				  CLTMeasureLY( // perform single pass according to prepared tiles operations and disparity // USED in lwir
//						  image_data, // first index - number of image in a quad
//						  saturation_imp, // boolean [][] saturation_imp, // (near) saturated pixels or null
						  clt_parameters,
						  combo_scan,     // final int           scanIndex,
						  // only combine and calculate once, next passes keep
						  (num_iter >0)? -1: bg_scan,        // final int           bgIndex, // combine, if >=0
						  tp.threadsMax,  // maximal number of threads to launch
						  false, // updateStatus,
						  debugLevelInner -1); // - 1); // -5-1
				  if (debugLevel > -2) {
					  tp.showScan(
							  tp.clt_3d_passes.get(combo_scan),   // CLTPass3d   scan,
							  "LY_combo_scan-"+combo_scan+"_post"); //String title)
				  }

				  int tilesX = tp.getTilesX();
				  int tilesY = tp.getTilesY();
				  int cluster_size =clt_parameters.tileStep;
				  int clustersX= (tilesX + cluster_size - 1) / cluster_size;
				  int clustersY= (tilesY + cluster_size - 1) / cluster_size;


				  ExtrinsicAdjustment ea = new ExtrinsicAdjustment(
						  geometryCorrection, // GeometryCorrection gc,
						  clt_parameters.tileStep,   // int         clusterSize,
						  clustersX, // 	int         clustersX,
						  clustersY); // int         clustersY);

				  double [] old_new_rms = new double[2];
				  boolean apply_extrinsic = (clt_parameters.ly_corr_scale != 0.0);
				  CLTPass3d   scan = tp.clt_3d_passes.get(combo_scan);
				  GeometryCorrection.CorrVector corr_vector =   ea.solveCorr (
						  clt_parameters.ly_marg_fract, 	  // double      marg_fract,        // part of half-width, and half-height to reduce weights
						  clt_parameters.ly_inf_en,           // boolean     use_disparity,     // adjust disparity-related extrinsics
						  clt_parameters.ly_aztilt_en,        // boolean     use_aztilts,       // Adjust azimuths and tilts excluding disparity
						  clt_parameters.ly_diff_roll_en,     // boolean     use_diff_rolls,    // Adjust differential rolls (3 of 4 angles)
//						  clt_parameters.ly_inf_force,        // boolean     force_convergence, // if true try to adjust convergence (disparity, symmetrical parameter 0) even with no disparity
						  clt_parameters.ly_min_forced,       // int         min_num_forced,    // minimal number of clusters with forced disparity to use it
						  // data, using just radial distortions
						  clt_parameters.ly_com_roll,         // boolean     common_roll,       // Enable common roll (valid for high disparity range only)
						  clt_parameters.ly_focalLength ,     // boolean     corr_focalLength,  // Correct scales (focal length temperature? variations)
						  clt_parameters.ly_ers_rot,          // boolean     ers_rot,           // Enable ERS correction of the camera rotation
						  clt_parameters.ly_ers_forw,         // boolean     ers_forw,          // Enable ERS correction of the camera linear movement in z direction
						  clt_parameters.ly_ers_side,         // boolean     ers_side,          // Enable ERS correction of the camera linear movement in x direction
						  clt_parameters.ly_ers_vert,         // boolean     ers_vert,          // Enable ERS correction of the camera linear movement in y direction
						  // add balancing-related here?
						  clt_parameters.ly_par_sel,          // 	int         manual_par_sel,    // Manually select the parameter mask bit 0 - sym0, bit1 - sym1, ... (0 - use boolean flags, != 0 - ignore boolean flags)
						  1.0,                                // 	double      weight_disparity,
						  1.0,                                // 	double      weight_lazyeye,
						  scan.getLazyEyeData(),              // dsxy, // double [][] measured_dsxy,
						  scan.getLazyEyeForceDisparity(),    // null, //	boolean [] force_disparity,    // boolean [] force_disparity,
						  false, // 	boolean     use_main, // corr_rots_aux != null;
						  geometryCorrection.getCorrVector(), // GeometryCorrection.CorrVector corr_vector,
						  old_new_rms, // double [] old_new_rms, // should be double[2]
						  debugLevel); //  + 5);// int debugLevel)

				  if (debugLevel > -2){
					  System.out.println("Old extrinsic corrections:");
					  System.out.println(geometryCorrection.getCorrVector().toString());
				  }
				  if (corr_vector != null) {
					  GeometryCorrection.CorrVector diff_corr = corr_vector.diffFromVector(geometryCorrection.getCorrVector());
					  comp_diff = diff_corr.getNorm();

					  if (debugLevel > -2){
							  System.out.println("New extrinsic corrections:");
							  System.out.println(corr_vector.toString());
					  }
					  if (debugLevel > -3){
							  System.out.println("Increment extrinsic corrections:");
							  System.out.println(diff_corr.toString());
							  // System.out.println("Correction scale = "+clt_parameters.ly_corr_scale);

					  }

					  if (apply_extrinsic){
						  geometryCorrection.setCorrVector(corr_vector) ;
						  System.out.println("Extrinsic correction updated (can be disabled by setting clt_parameters.ly_corr_scale = 0.0) ");
					  } else {
						  System.out.println("Correction is not applied according clt_parameters.ly_corr_scale == 0.0) ");
					  }
				  } else {
					  if (debugLevel > -3){
						  System.out.println("LMA failed");
					  }
				  }

				  boolean done = (comp_diff < min_sym_update) || (num_iter == (max_tries - 1));
				  //				  System.out.println("done="+done);
				  if (debugLevel > -10) { // should work even in batch mode
					  System.out.println("#### extrinsicsCLT(): iteration step = "+(num_iter + 1) + " ( of "+max_tries+") change = "+
							  comp_diff + " ("+min_sym_update+"), previous RMS = " + old_new_rms[0]+
							  " final RMS = " + old_new_rms[1]+ " (debugLevel = "+debugLevel+")");
				  }
				  if (debugLevel > -10) {
					  if ((debugLevel > -3) || done) {
						  //						  System.out.println("#### extrinsicsCLT(): iteration step = "+(num_iter + 1) + " ( of "+max_tries+") change = "+
						  //								  comp_diff + " ("+min_sym_update+"), previous RMS = " + new_corr[0][1][0]);
						  System.out.println("New extrinsic corrections:");
						  System.out.println(geometryCorrection.getCorrVector().toString());
					  }
				  }

				  if (comp_diff < min_sym_update) {
					  break;
				  }
				  if (update_disp_from_latest) {
/*					  
					  CLTMeasure( // perform single pass according to prepared tiles operations and disparity
							  clt_parameters,
							  combo_scan,
							  false,             // final boolean     save_textures,
							  true,              // final boolean     save_corr,
							  null,              // combo_mismatch,    // final double [][] mismatch,    // null or double [12][]
							  tp.threadsMax,     // maximal number of threads to launch
							  false,             // updateStatus,
							  debugLevelInner - 1);
*/							  
					  CLTMeasureCorr( // perform single pass according to prepared tiles operations and disparity
							  clt_parameters,
							  combo_scan,
							  false,             // final boolean     save_textures,
							  tp.threadsMax,     // maximal number of threads to launch
							  false,             // updateStatus,
							  debugLevelInner - 1);
				  }

			  } else { // Old, no-GPU
				  double [][] bg_mismatch = new double[12][];
				  double [][] combo_mismatch = new double[12][];
				  CLTMeasure( // perform single pass according to prepared tiles operations and disparity
						  clt_parameters,
						  bg_scan,
						  false,             // final boolean     save_textures,
						  true,              // final boolean     save_corr,
						  bg_mismatch,    // final double [][] mismatch,    // null or double [12][]
						  tp.threadsMax,  // maximal number of threads to launch
						  false, // updateStatus,
						  debugLevelInner - 1);
				  CLTMeasure( // perform single pass according to prepared tiles operations and disparity
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
					  (new ShowDoubleFloatArrays()).showArrays(scans14, tp.getTilesX(), tp.getTilesY(), true, "scans_14"); //  , titles);
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

				  int num_tiles = tp.clt_3d_passes.get(combo_scan).getStrength().length;


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

						  clt_parameters.getLyPerQuad(num_tiles), // final int        min_per_quadrant, // minimal tiles per quadrant (not counting the worst) tp proceed
						  clt_parameters.getLyInf(num_tiles),	  // final int        min_inf,          // minimal number of tiles at infinity to proceed
						  clt_parameters.getLyInfScale(num_tiles),// final int        min_inf_to_scale, // minimal number of tiles at infinity to apply weight scaling

						  clt_parameters.ly_right_left,      // false // equalize weights of right/left FoV (use with horizon in both halves and gross infinity correction)
						  clt_parameters,  // EyesisCorrectionParameters.CLTParameters           clt_parameters,
						  scans14, // disp_strength, // scans,   // double [][] disp_strength,
						  target_disparity,          // double [][]      target_disparity, // null or programmed disparity (1 per each 14 entries of scans_14)
						  tp.getTilesX(), // int         tilesX,
						  clt_parameters.corr_magic_scale, // double      magic_coeff, // still not understood coefficent that reduces reported disparity value.  Seems to be around 8.5
						  debugLevelInner - 1); //  + (clt_parameters.fine_dbg ? 1:0)); // int debugLevel)
				  if (new_corr == null) {
					  return false;
				  }
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
					  boolean done = (comp_diff < min_sym_update) || (num_iter == (max_tries - 1));
					  //				  System.out.println("done="+done);
					  if (debugLevel > -10) { // should work even in batch mode
						  System.out.println("#### extrinsicsCLT(): iteration step = "+(num_iter + 1) + " ( of "+max_tries+") change = "+
								  comp_diff + " ("+min_sym_update+"), previous RMS = " + new_corr[0][1][0]+ " (debugLevel = "+debugLevel+")");
					  }
					  if (debugLevel > -10) {
						  if ((debugLevel > -1) || done) {
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
		  }
		  return true; // (comp_diff < (adjust_poly ? min_poly_update : min_sym_update));
	  }


	  public double [][] getRigDSFromTwoQuadCL( // not used in lwir
			  TwoQuadCLT       twoQuadCLT, //maybe null in no-rig mode, otherwise may contain rig measurements to be used as infinity ground truth
			  CLTParameters           clt_parameters,
			  final int        debugLevel) {
		  if ((twoQuadCLT == null) || (twoQuadCLT.getBiScan(0) == null)){
			  System.out.println("Rig data is not available, aborting");
			  return null;
		  }

		  BiScan scan = twoQuadCLT.getBiScan(0);
		  double [][] rig_disp_strength = 		scan.getDisparityStrength(
		    		true,   // final boolean only_strong,
		    		true,   // final boolean only_trusted,
		    		true) ; // final boolean only_enabled);
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
		  return rig_disp_strength;


	  }
/**
 *
 * @param geometryCorrection_main
 * @param rig_disp_strength
 * @param clt_parameters
 * @param adjust_poly
 * @param threadsMax
 * @param updateStatus
 * @param debugLevel
 * @return true on success, false on failure
 */
	  public boolean extrinsicsCLTfromGT( // USED in lwir
//			  TwoQuadCLT       twoQuadCLT, //maybe null in no-rig mode, otherwise may contain rig measurements to be used as infinity ground truth
			  GeometryCorrection geometryCorrection_main, // only used for aux camera if coordinates are for main (null for LWIR)
			  double [][] rig_disp_strength,
			  CLTParameters           clt_parameters,
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
		  double min_sym_update =           clt_parameters.getLymChange(is_aux); //  4e-6; // stop iterations if no angle changes more than this
		  double min_poly_update =          clt_parameters.lym_poly_change; //  Parameter vector difference to exit from polynomial correction
		  /*
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
*/
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
//					  image_data,              // first index - number of image in a quad
//					  saturation_imp,          // boolean [][] saturation_imp, // (near) saturated pixels or null
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
				  (new ShowDoubleFloatArrays()).showArrays(scans14, tp.getTilesX(), tp.getTilesY(), true, "scans_14"); //  , titles);
			  }

			  if (!batch_mode && clt_parameters.show_extrinsic && (debugLevel > 1)) {
				  tp.showScan(
						  comboScan,   // CLTPass3d   scan,
						  "combo_scan-"+num_iter+"_iter"); //String title)
			  }

			  double [][][] new_corr;
			  double [][][]    gt_disparity_strength = {rig_disp_strength};
			  int num_tiles = comboScan.getStrength().length;

			  new_corr = ac.lazyEyeCorrectionFromGT(
///					  geometryCorrection_main, //final GeometryCorrection geometryCorrection_main, // if not null - this is an AUX camera of a rig
					  adjust_poly,                       // final boolean use_poly,
					  true, // final boolean    restore_disp_inf, // Restore subtracted disparity for scan #0 (infinity)
					  clt_parameters.fcorr_radius,       // 	final double fcorr_radius,
					  clt_parameters.fcorr_inf_strength, //  final double min_strenth,

					  clt_parameters.inf_str_pow,        // 1.0, //   final double     strength_pow,
					  0.8*clt_parameters.disp_scan_step, // 1.5, // final double     lazyEyeCompDiff, // clt_parameters.fcorr_disp_diff
					  clt_parameters.ly_smpl_side,       // 3,   // final int        lazyEyeSmplSide, //        = 2;      // Sample size (side of a square)
					  clt_parameters.ly_smpl_num,        // 5,   // final int        lazyEyeSmplNum, //         = 3;      // Number after removing worst (should be >1)
					  clt_parameters.ly_smpl_rms,        // 0.1, // final double     lazyEyeSmplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
					  clt_parameters.ly_disp_var_gt,        // 0.2, // final double     lazyEyeDispVariation, // 0.2, maximal full disparity difference between tgh tile and 8 neighborxs
					  clt_parameters.ly_disp_rvar_gt,       // 0.2, // final double     lazyEyeDispRelVariation, // 0.02 Maximal relative full disparity difference to 8 neighbors
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
					  clt_parameters.getLyPerQuad(num_tiles), // final int        min_per_quadrant, // minimal tiles per quadrant (not counting the worst) tp proceed
					  0, // clt_parameters.getLyInf(num_tiles),	  // final int        min_inf,          // minimal number of tiles at infinity to proceed
					  clt_parameters.getLyInfScale(num_tiles),// final int        min_inf_to_scale, // minimal number of tiles at infinity to apply weight scaling

					  clt_parameters.ly_inf_max_disparity, // inf_max_disparity, // final double     inf_max_disparity, // use all smaller disparities as inf_fraction
					  clt_parameters,  // EyesisCorrectionParameters.CLTParameters           clt_parameters,
					  scans14, // disp_strength, // scans,   // double [][] disp_strength,
					  gt_disparity_strength, // double [][][]    gt_disparity_strength, // 1 pair for each 14 entries of scans_14 (normally - just 1 scan
					  filter_ds,   // final boolean    filter_ds, //
					  filter_lyf,   // final boolean    filter_lyf, // ~clt_parameters.lyf_filter, but may be different, now off for a single cameras
					  tp.getTilesX(), // int         tilesX,
					  clt_parameters.corr_magic_scale, // double      magic_coeff, // still not understood coefficent that reduces reported disparity value.  Seems to be around 8.5
					  debugLevelInner - 1); //  + (clt_parameters.fine_dbg ? 1:0)); // int debugLevel)
			  if (new_corr == null) {
				  return false; // not used in lwir
			  }
			  comp_diff = 0.0;
			  int num_pars = 0;
			  if (adjust_poly) { // not used in lwir
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
			  } else { // USED in lwir
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
		  return true; // (comp_diff < (adjust_poly ? min_poly_update : min_sym_update));
	  }






	  public boolean expandCLTQuad3d( // USED in lwir
		  CLTParameters           clt_parameters,
		  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
		  ColorProcParameters colorProcParameters,
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
        		  else { // not used in lwir
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
///          Runtime runtime = Runtime.getRuntime();
///          runtime.gc();
///          System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");
          return true; //  null;
	  }

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// called after composite scan is added to the list? composite scan is inside

	  public int zMapExpansionStep( // USED in lwir
				final ArrayList <CLTPass3d> passes,// List, first, last - to search for the already tried disparity
				final CLTParameters           clt_parameters, // for refinePassSetup()
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

    	  double [][] filtered_disp_strength = tp.getFilteredDisparityStrength( // disp all 0?, str -1/0
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
		  if (last_pass){
			  System.out.println("+++++++++++ Last pass - pre  compositeScan() ++++++++++++");
		  }

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

    		  (new ShowDoubleFloatArrays()).showArrays(dbg_img,  tilesX, tilesY, true, title,titles);
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
        	  num_extended = numLeftRemoved[0]; //0,0
			  if (clt_parameters.show_expand || (clt_parameters.show_variant && (numLeftRemoved[1] > 1 ))) tp.showScan( // not used in lwir
					  tp.clt_3d_passes.get(refine_pass), // CLTPass3d   scan,
					  "prepareExpandVariant-"+numLeftRemoved[1]+"-"+refine_pass); //String title)
    	  }
    	  if ((num_extended == 0) && expand_legacy) { // if both are on, will use legacy if neighbors failed // not used in lwir

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

		  CLTMeasureCorr( // perform single pass according to prepared tiles operations and disparity  BUG: gets with .disparity==null
//    			  image_data, // first index - number of image in a quad
//				  saturation_imp, //final boolean [][]  saturation_imp, // (near) saturated pixels or null
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
// Separate method to detect and remove periodic structures

	  public boolean showPeriodic( // not used in lwir
			  CLTParameters           clt_parameters,
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel) {
		  final boolean        usePoly =            false; // use polynomial method to find max), valid if useCombo == false

		  double [][] periodics = tp.getPeriodics(
				  tp.clt_3d_passes,                     // final ArrayList <CLTPass3d> passes,// List, first, last - to search for the already tried disparity
				  0,                                    // final int         firstPass,
				  tp.clt_3d_passes.size(),              // final int         lastPassPlus1,
				  clt_parameters.per_trustedCorrelation,// final double                trustedCorrelation,
				  clt_parameters.per_initial_diff,      // final double                initial_diff, // initial disparity difference to merge to maximum
				  clt_parameters.per_strength_floor,    // final double                strength_floor,
				  clt_parameters.per_strength_max_over, // final double                strength_max_over, // maximum should have strength by this more than the floor

				  clt_parameters.per_min_period,        // final double                min_period,
				  clt_parameters.per_min_num_periods,   // final int                   min_num_periods, // minimal number of periods
				  clt_parameters.per_disp_tolerance,    // final double                disp_tolerance, // maximal difference between the average of fundamental and 2-nd and first
				  // TODO: replace next parameter
				  clt_parameters.per_disp_tolerance,    // final double                disp_tol_center, // tolerance to match this (center) tile ds to that of the merged with neighbors - should be < min_period/2
				  clt_parameters.per_disp_match,        // final double                disp_match,     // disparity difference to match neighbors
				  clt_parameters.per_strong_match_inc,  // final double                strong_match_inc,   // extra strength to treat match as strong (for hysteresis)
				  usePoly,                              // final boolean               usePoly,  // use polynomial method to find max), valid if useCombo == false
				  clt_parameters.tileX,                 // final int               dbg_tileX,
				  clt_parameters.tileY,                 // final int               dbg_tileY,
				  threadsMax,                           // final int                   threadsMax,  // maximal number of threads to launch
				  updateStatus,                         // final boolean               updateStatus,
				  debugLevel+3);                        // final int                   debugLevel) // update status info

		  tp.periodics = periodics;

		  String [] dbg_titles= {"fundamental","period", "strength", "num_layers"};
			double [][] dbg_img = new double [4][tp.getTilesX()*tp.getTilesY()];
			for (int n = 0; n < periodics.length; n++) {
				for (int i = 0; i < dbg_img.length; i++) {
					dbg_img[i][n]= (periodics[n] == null) ? Double.NaN : periodics[n][i];
				}
			}

		  (new ShowDoubleFloatArrays()).showArrays(
				  dbg_img,
				  tp.getTilesX(),
				  tp.getTilesY(),
				  true,
				  image_name+sAux()+"-PERIODIC",
				  dbg_titles);
		  return true;
	  }






//*****************************************************************

//	  public ImagePlus output3d(
	  public boolean output3d( // USED in lwir
			  CLTParameters                            clt_parameters,
			  ColorProcParameters                      colorProcParameters,
			  EyesisCorrectionParameters.RGBParameters rgbParameters,
			  final int                                threadsMax,  // maximal number of threads to launch
			  final boolean                            updateStatus,
			  final int                                debugLevel)
	  {
		  final boolean    batch_mode = clt_parameters.batch_run;
		  this.startStepTime=System.nanoTime();
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  if (this.image_data == null){
			  return false; // not used in lwir
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
		  // Create tasks to scan, have tasks, disparity and border tiles in tp.clt_3d_passes
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
//TODO make it w/o need for  bgndScan.texture as GPU will calculate texture right before output
//use selection? or texture_selection instead?		  
		  if (bgndScan.texture != null) { // TODO: same for the backdrop too
			  if (clt_parameters.infinityDistance > 0.0){ // generate background as a billboard
				  // grow selection, then grow once more and create border_tiles
				  // create/rstore, probably not needed
				  boolean [] bg_sel_backup = bgndScan.getSelected().clone();
				  boolean [] bg_border_backup = (bgndScan.getBorderTiles() == null) ? null: bgndScan.getBorderTiles().clone();
				  boolean [] bg_selected = bgndScan.getSelected();
				  boolean [] border_tiles = bg_selected.clone();
				  tp.growTiles(
						  2,                   // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
						  bg_selected,
						  null); // prohibit

				  for (int i = 0; i < border_tiles.length; i++){
					  border_tiles[i] = !border_tiles[i] && bg_selected[i];
				  }
				  // update texture_tiles (remove what is known not to be background
				  if (bgndScan.texture_tiles != null) { // for CPU
					  for (int ty = 0; ty < tilesY; ty++){
						  for (int tx = 0; tx < tilesX; tx++){
							  if (!bg_selected[tx + tilesX*ty]){
								  bgndScan.texture_tiles[ty][tx] = null; //
							  }

						  }
					  }
				  } else {//  for GPU
					  for (int i = 0; i < bg_selected.length; i++) {
						  if (!bg_selected[i]) {
							  bgndScan.setTextureSelection(i,false);
						  }
					  }
				  }
				  
//TODO2020: set texture_selection
				  bgndScan.setBorderTiles(border_tiles);
				  // limit tile_op to selection?
				  // updates selection from non-null texture tiles
				  String texturePath = getPassImage( // get image from a single pass - both CPU and GPU
						  clt_parameters,
						  colorProcParameters,
						  rgbParameters,
						  correctionsParameters.getModelName(this.image_name)+"-img_infinity", // +scanIndex,
						  bgndIndex,
						  threadsMax,  // maximal number of threads to launch
						  updateStatus,
						  batch_mode ? -5: debugLevel);
				  if (texturePath != null) { // null if empty image
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
		  }

		  // With GPU - do nothing here or copy selected -> texture_selection?
		  for (int scanIndex = next_pass; scanIndex < tp.clt_3d_passes.size(); scanIndex++){
			  if (debugLevel > 0){
				  System.out.println("FPGA processing scan #"+scanIndex);
			  }
			  /*
			  CLTMeasure( // perform single pass according to prepared tiles operations and disparity
					  clt_parameters,
					  scanIndex,
					  true,  // final boolean     save_textures,
					  false, // final boolean     save_corr,
					  threadsMax,  // maximal number of threads to launch
					  updateStatus,
					  batch_mode ? -5: debugLevel);
		*/
			  CLTMeasureTextures( // has GPU version - will just copy selection
					  clt_parameters,
					  scanIndex,
					  threadsMax,  // maximal number of threads to launch
					  updateStatus,
					  batch_mode ? -5: debugLevel);

		  }

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
			  if (texturePath == null) { // not used in lwir
				  continue; // empty image
			  }
			  CLTPass3d scan = tp.clt_3d_passes.get(scanIndex);

			  // TODO: use new updated disparity, for now just what was forced for the picture
			  double [] scan_disparity = new double [tilesX * tilesY];
			  int indx = 0;
			  for (int ty = 0; ty < tilesY; ty ++) for (int tx = 0; tx < tilesX; tx ++){
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
			  boolean showTri = !batch_mode && (debugLevel > -1) && (((scanIndex < next_pass + 1) && clt_parameters.show_triangles) ||((scanIndex - next_pass) == 73));

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
				e.printStackTrace();
				return false;
			}
		  }

		  if ((x3d_path != null) && (x3dOutput != null)){
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




	  public void generateClusterX3d( // USED in lwir
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
		  if (bounds == null) {
			  return; // not used in lwir
		  }
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


//	  public ImagePlus getBackgroundImage( // USED in lwir
	  public boolean[] getBackgroundImageMasks( // USED in lwir
//			  boolean    no_image_save,
			  CLTParameters           clt_parameters,
//			  ColorProcParameters colorProcParameters,
//			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  String     name,
			  int        disparity_index, // index of disparity value in disparity_map == 2 (0,2 or 4)
			  int        threadsMax,  // maximal number of threads to launch
			  boolean    updateStatus,
			  int        debugLevel
			  )
	  {
//		  final boolean new_mode = false;
		  boolean dbg_gpu_transition = true;


		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  ShowDoubleFloatArrays sdfa_instance = null;

		  if ((clt_parameters.debug_filters && (debugLevel > -1)) || dbg_gpu_transition)
		  if ((debugLevel > -1))
			  sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?

		  CLTPass3d bgnd_data = tp.clt_3d_passes.get(0);
//		  double [][][][] texture_tiles = bgnd_data.texture_tiles;

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

		  bgnd_data.selected = bgnd_tiles_grown; // selected for background w/o extra transparent layer
		  boolean [] bgnd_tiles_grown2 = bgnd_tiles_grown.clone(); // only these have non 0 alpha
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
		  // not here - will be set/calculated for GPU only
///		  bgnd_data.setTextureSelection(bgnd_tiles_grown2); // selected for background w/o extra transparent layer
		  return bgnd_tiles ;
	  }
	  
	  
	  public ImagePlus getBackgroundImage( // USED in lwir
			  boolean []                                bgnd_tiles, 
			  CLTParameters                             clt_parameters,
			  ColorProcParameters                       colorProcParameters,
			  EyesisCorrectionParameters.RGBParameters  rgbParameters,
			  String     name,
			  int        disparity_index, // index of disparity value in disparity_map == 2 (0,2 or 4)
			  int        threadsMax,  // maximal number of threads to launch
			  boolean    updateStatus,
			  int        debugLevel
			  ) {
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  final boolean new_mode = false;
		  int num_bgnd = 0;

		  CLTPass3d bgnd_data = tp.clt_3d_passes.get(0);
		  double [][][][] texture_tiles = bgnd_data.texture_tiles;
		  double [][][][] texture_tiles_bgnd = new double[tilesY][tilesX][][];
		  double [] alpha_zero = new double [4*clt_parameters.transform_size*clt_parameters.transform_size];
		  int alpha_index = 3;
		  boolean [] bgnd_tiles_grown = bgnd_tiles.clone(); // only these have non 0 alpha
		  tp.growTiles(
				  2,      // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				  bgnd_tiles,
				  null); // prohibit
		  boolean [] bgnd_tiles_grown2 = bgnd_tiles_grown.clone(); // only these have non 0 alpha
		  tp.growTiles(
				  2,      // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				  bgnd_tiles_grown2,
				  null); // prohibit
		  
		  for (int i = 0; i < alpha_zero.length; i++) alpha_zero[i]=0.0;
		  // Seems to be wrong for the second?
		  
		  if (new_mode) { // not used
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
					  if ((texture_tiles[tileY][tileX] != null) && // null pointer
							  bgnd_tiles[tileY * tilesX + tileX]) {
						  if (bgnd_tiles_grown2[tileY * tilesX + tileX]) {
							  texture_tiles_bgnd[tileY][tileX]= texture_tiles[tileY][tileX];
							  num_bgnd++;
						  }else{ // not used in lwir
							  texture_tiles_bgnd[tileY][tileX]= texture_tiles[tileY][tileX].clone();
							  texture_tiles_bgnd[tileY][tileX][alpha_index] = alpha_zero;
						  }
					  }
				  }
			  }
		  }

		  if (num_bgnd < clt_parameters.min_bgnd_tiles){
			  return null; // no background to generate // not used in lwir
		  }
		  
		  
		  ImageDtt image_dtt = new ImageDtt(
				  clt_parameters.transform_size,
				  isMonochrome(),
				  isLwir(),
				  clt_parameters.getScaleStrength(isAux()));
		  double [][] texture_overlap = image_dtt.combineRBGATiles(
				  texture_tiles_bgnd, // texture_tiles,               // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}
///				  image_dtt.transform_size,
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
				  tilesX *  image_dtt.transform_size,
				  tilesY *  image_dtt.transform_size,
				  1.0,         // double scaleExposure, // is it needed?
				  debugLevel);
		  // resize for backdrop here!
//	public double getFOVPix(){ // get ratio of 1 pixel X/Y to Z (distance to object)
		  return imp_texture_bgnd;

	  }
	  
	  
	  public ImagePlus finalizeBackgroundImage( // USED in lwir
			  ImagePlus imp_texture_bgnd,
			  boolean    no_image_save,
			  CLTParameters           clt_parameters,
			  String     name,
			  int        debugLevel) {
		  
		  ImagePlus imp_texture_bgnd_ext = resizeForBackdrop(
				  imp_texture_bgnd,
				  clt_parameters.black_back, //  boolean fillBlack,
				  clt_parameters.black_back, //  boolean noalpha,
				  debugLevel);
		  String path= correctionsParameters.selectX3dDirectory(
				  //TODO: Which one to use - name or this.image_name ?
 				  correctionsParameters.getModelName(this.image_name), // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
 				  correctionsParameters.x3dModelVersion,
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



	 public String getPassImage( // get image from a single pass, return relative path for x3d // USED in lwir
			  CLTParameters           clt_parameters,
			  ColorProcParameters colorProcParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  String     name,
			  int        scanIndex,
			  int        threadsMax,  // maximal number of threads to launch
			  boolean    updateStatus,
			  int        debugLevel)
	  {
		 final int tilesX = tp.getTilesX();
		 final int tilesY = tp.getTilesY();

		  ShowDoubleFloatArrays sdfa_instance = null;
    	  if (debugLevel > -1) sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
		  CLTPass3d scan = tp.clt_3d_passes.get(scanIndex);
		  boolean [] borderTiles = scan.border_tiles;
		  double [][][][] texture_tiles = scan.texture_tiles;
		  // only place that uses updateSelection()
		  scan.updateSelection(); // update .selected field (all selected, including border) and Rectangle bounds
		  if (scan.getTextureBounds() == null) { // not used in lwir
			  System.out.println("getPassImage(): Empty image!");
			  return null;
		  }
		  ImageDtt image_dtt = new ImageDtt(
				  clt_parameters.transform_size,
				  isMonochrome(),
				  isLwir(),
				  clt_parameters.getScaleStrength(isAux()));
		  double [][]alphaFade = tp.getAlphaFade(image_dtt.transform_size);
		  if ((debugLevel > 0) && (scanIndex == 1)) { // not used in lwir
			  String [] titles = new String[16];
			  for (int i = 0; i<titles.length;i++)  titles[i]=""+i;
			  sdfa_instance.showArrays(alphaFade, 2*image_dtt.transform_size,2*image_dtt.transform_size,true,"alphaFade",titles);
		  }
		  double [][][][] texture_tiles_cluster = new double[tilesY][tilesX][][];
		  double [] alpha_zero = new double [4*image_dtt.transform_size*image_dtt.transform_size];
		  int alpha_index = 3;
		  for (int i = 0; i < alpha_zero.length; i++) alpha_zero[i]=0.0;
		  // border tiles are copied, alpha from alphaFade (not multiplied?)
		  for (int tileY = 0; tileY < tilesY; tileY++){
			  for (int tileX = 0; tileX < tilesX; tileX++){
				  texture_tiles_cluster[tileY][tileX]= null;
				  if (texture_tiles[tileY][tileX] != null) {
					  if (borderTiles[tileY * tilesX + tileX]) {
						  texture_tiles_cluster[tileY][tileX]= texture_tiles[tileY][tileX].clone();
						  if (clt_parameters.shAggrFade) { // not used in lwir
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

		  double [][] texture_overlap = image_dtt.combineRBGATiles(
				  texture_tiles_cluster, // texture_tiles,               // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}
///				  image_dtt.transform_size,
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
		  // for now - use just RGB. Later add option for RGBA (?)
		  double [][] texture_rgb = {texture_overlap[0],texture_overlap[1],texture_overlap[2]};
		  double [][] texture_rgba = {texture_overlap[0],texture_overlap[1],texture_overlap[2],texture_overlap[3]};
		  double [][] texture_rgbx = ((clt_parameters.alpha1 > 0)? texture_rgba: texture_rgb);

		  boolean resize = true;
		  if (resize) {
		  texture_rgbx = resizeGridTexture(
				  texture_rgbx,
				  image_dtt.transform_size,
				  tilesX,
				  tilesY,
				  scan.getTextureBounds());
		  }

		  int width = resize ? (image_dtt.transform_size * scan.getTextureBounds().width): (image_dtt.transform_size * tilesX);
		  int height = resize ? (image_dtt.transform_size * scan.getTextureBounds().height): (image_dtt.transform_size * tilesY);
		  if ((width <= 0) || (height <= 0)) {
			  System.out.println("***** BUG in getPassImage(): width="+width+", height="+height+", resize="+resize+" ****"); // not used in lwir
		  }

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
				  width, //tp.tilesX *  image_dtt.transform_size,
				  height, //tp.tilesY *  image_dtt.transform_size,
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



	  public ImagePlus resizeForBackdrop( // USED in lwir
			  ImagePlus imp,
			  boolean fillBlack,
			  boolean noalpha, // only with fillBlack, otherwise ignored
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
		  } else { // not used in lwir
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

	  
	  
	  public ImagePlus resizeToFull(
			  int       out_width,
			  int       out_height,
			  int       x0, // image offset-x pixels
			  int       y0, // image offset-y pixels
			  ImagePlus imp,
			  boolean   fillBlack,
			  boolean noalpha, // only with fillBlack, otherwise ignored
			  int debugLevel)
	  {
		  int width = imp.getWidth();
		  int height = imp.getHeight();
//		  int h_margin = (int) Math.round((backdropPixels - width)/2);
//		  int v_margin = (int) Math.round((backdropPixels - height)/2);
//		  int width2 =  width +  2 * h_margin;
//		  int height2 = height + 2 * v_margin;
		  int [] src_pixels = (int []) imp.getProcessor().getPixels();
		  int [] pixels = new int [out_width * out_height];
		  int black = noalpha ? 0 :        0xff000000;
		  int mask =  noalpha ? 0xffffff : 0xffffffff;
		  if (fillBlack) {
			  for (int i = 0; i < pixels.length; i++){
				  pixels[i] = black;
			  }
		  }
		  int indx = 0;
		  int offset = y0 * out_width + x0; // v_margin *  width2 + h_margin;
		  if (fillBlack) {
			  for (int i = 0; i < height; i++){
				  for (int j = 0; j < width; j++){
					  int a =  (src_pixels[indx] >> 24) & 0xff;
					  if (a == 255) {
						  pixels[offset+ i * out_width + j] = src_pixels[indx] & mask;
					  } else if (a == 0) {
						  pixels[offset+ i * out_width + j] = black;
					  } else  {
						  int r = (src_pixels[indx] >> 16) & 0xff; //' maybe - swap with b?
						  int g = (src_pixels[indx] >>  8) & 0xff;
						  int b = (src_pixels[indx] >>  0) & 0xff;
						  r = (r * a) / 255;
						  g = (g * a) / 255;
						  b = (b * a) / 255;
						  pixels[offset+ i * out_width + j] = black | (r << 16) | (g << 8) | b;
					  }
					  indx++;
				  }
			  }
		  } else { // not used in lwir
			  for (int i = 0; i < height; i++){
				  for (int j = 0; j < width; j++){
					  pixels[offset+ i * out_width + j] = src_pixels[indx++];
				  }
			  }
		  }
		  ColorProcessor cp=new ColorProcessor(out_width,out_height);
		  cp.setPixels(pixels);
		  ImagePlus imp_ext=new ImagePlus(imp.getTitle()+"-ext",cp);
		  return imp_ext;
	  }
	  
//[tp.tilesY][tp.tilesX]["RGBA".length()][]
//linearStackToColor

	  public CLTPass3d CLTBackgroundMeas( // measure background // USED in lwir
			  CLTParameters           clt_parameters,
			  final int           threadsMax,  // maximal number of threads to launch
			  final boolean       updateStatus,
			  final int           debugLevel)
	  {
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  CLTPass3d scan_rslt = new CLTPass3d(tp);
		  int d = ImageDtt.setImgMask(0, 0xf);
		  d =     ImageDtt.setPairMask(d,0xf);
		  d =     ImageDtt.setForcedDisparity(d,true);
		  int [][]     tile_op =         tp.setSameTileOp(clt_parameters,  d, debugLevel);
		  double [][]  disparity_array = tp.setSameDisparity(0.0); // [tp.tilesY][tp.tilesX] - individual per-tile expected disparity
		  // undecided, so 2 modes of combining alpha - same as rgb, or use center tile only
		  double [][][][]     clt_corr_combo =    new double [ImageDtt.TCORR_TITLES.length][tilesY][tilesX][]; // will only be used inside?

		  double min_corr_selected = clt_parameters.min_corr; // 0.02

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
		  ImageDtt image_dtt = new ImageDtt(
				  clt_parameters.transform_size,
				  isMonochrome(),
				  isLwir(),
				  clt_parameters.getScaleStrength(isAux()));
		  double z_correction =  clt_parameters.z_correction;
		  if (clt_parameters.z_corr_map.containsKey(image_name)){ // not used in lwir
			  z_correction +=clt_parameters.z_corr_map.get(image_name);
		  }
		  final double disparity_corr = (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);
		  
		  double [][][][][] clt_corr_partial = null;
		  if (clt_parameters.img_dtt.gpu_mode_debug) {
			  clt_corr_combo = null;
			  clt_corr_partial = new double [tilesY][tilesX][][][];
			  for (int i = 0; i < tilesY; i++){
				  for (int j = 0; j < tilesX; j++){
					  clt_corr_partial[i][j] = null;
				  }
			  }
		  }
		  clt_corr_combo = null; // Something is broke in old code, tries to use color==3 ==colors.weights.length
		  image_dtt.clt_aberrations_quad_corr(
				  clt_parameters.img_dtt,       // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
				  1,                            // final int  macro_scale, // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
				  tile_op,                      // per-tile operation bit codes
				  disparity_array,              // clt_parameters.disparity,     // final double            disparity,
				  image_data,                   // final double [][][]      imade_data, // first index - number of image in a quad
				  saturation_imp,               // boolean [][] saturation_imp, // (near) saturated pixels or null
				  // correlation results - final and partial
				  clt_corr_combo,               // [tp.tilesY][tp.tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				  clt_corr_partial, // null, // clt_corr_partial,    // [tp.tilesY][tp.tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				  null,    // [tp.tilesY][tp.tilesX][pair]{dx,dy,weight}[(2*transform_size-1)*(2*transform_size-1)] // transpose unapplied. null - do not calculate
				  //	Use it with disparity_maps[scan_step]?		  clt_mismatch,    // [tp.tilesY][tp.tilesX][pair]{dx,dy,weight}[(2*transform_size-1)*(2*transform_size-1)] // transpose unapplied. null - do not calculate
				  disparity_map,    // [12][tp.tilesY * tp.tilesX]
				  texture_tiles,        // [tp.tilesY][tp.tilesX]["RGBA".length()][];
				  tilesX * image_dtt.transform_size, // imp_quad[0].getWidth(),       // final int width,
				  clt_parameters.getFatZero(isMonochrome()),      // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
				  clt_parameters.corr_sym,
				  clt_parameters.corr_offset,
				  clt_parameters.corr_red,
				  clt_parameters.corr_blue,
				  clt_parameters.getCorrSigma(image_dtt.isMonochrome()),
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
				  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
				  null,                          // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
				  clt_kernels,                   // final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				  clt_parameters.kernel_step,
				  ///				  image_dtt.transform_size,
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
		  if (clt_corr_partial!=null){ // only to debug matching gpu/cpu
			  if (debugLevel > -3){ // -1
				  String [] allColorNames = {"red","blue","green","combo"};
				  String [] titles = new String[clt_corr_partial.length];
				  for (int i = 0; i < titles.length; i++){
					  titles[i]=allColorNames[i % allColorNames.length]+"_"+(i / allColorNames.length);
				  }
				  double [][] corr_rslt_partial = image_dtt.corr_partial_dbg(
						  clt_corr_partial,
						  2*image_dtt.transform_size - 1,	//final int corr_size,
						  4,	// final int pairs,
						  4,    // final int colors,
						  clt_parameters.corr_border_contrast,
						  threadsMax,
						  debugLevel);
				  // titles.length = 15, corr_rslt_partial.length=16!
				  System.out.println("corr_rslt_partial.length = "+corr_rslt_partial.length+", titles.length = "+titles.length);
				  (new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
						  corr_rslt_partial,
						  tilesX*(2*image_dtt.transform_size),
						  tilesY*(2*image_dtt.transform_size),
						  true,
						  image_name+sAux()+"-PART_CORR-CPU-D"+clt_parameters.disparity);
			  }
		  }

		  return scan_rslt;
	  }
	  
	  public CLTPass3d  CLTMeasureTextures( // perform single pass according to prepared tiles operations and disparity // not used in lwir
			  CLTParameters           clt_parameters,
			  final int         scanIndex,
			  final int         threadsMax,  // maximal number of threads to launch
			  final boolean     updateStatus,
			  final int         debugLevel)
	  {
		  return CLTMeasure( // measure background // USED in lwir
				  clt_parameters,
				  scanIndex,
				  true, // save_textures,
				  false, // save_corr,
				  threadsMax,  // maximal number of threads to launch
				  updateStatus,
				  debugLevel);
	  }
	  
	  public CLTPass3d  CLTMeasure( // perform single pass according to prepared tiles operations and disparity // not used in lwir
			  CLTParameters           clt_parameters,
			  final int         scanIndex,
			  final boolean     save_textures,
			  final boolean     save_corr,
			  final int         threadsMax,  // maximal number of threads to launch
			  final boolean     updateStatus,
			  final int         debugLevel)
	  {
		  return CLTMeasure( // perform single pass according to prepared tiles operations and disparity
				  clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
				  scanIndex,      // final int         scanIndex,
				  save_textures,  // final boolean     save_textures,
				  save_corr,      // final boolean       save_corr,
				  null,           // final double [][] mismatch,    // null or double [12][]
				  threadsMax,     // final int         threadsMax,  // maximal number of threads to launch
				  updateStatus,   // final boolean     updateStatus,
				  debugLevel);    // final int         debugLevel);
	  }
	  public CLTPass3d  CLTMeasureCorr( // perform single pass according to prepared tiles operations and disparity // not used in lwir
			  CLTParameters           clt_parameters,
			  final int         scanIndex,
			  final boolean     save_textures,
			  final int         threadsMax,  // maximal number of threads to launch
			  final boolean     updateStatus,
			  final int         debugLevel)
	  {
		  return CLTMeasure( // perform single pass according to prepared tiles operations and disparity
				  clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
				  scanIndex,      // final int         scanIndex,
				  save_textures,  // final boolean     save_textures,
				  true,           // final boolean       save_corr,
				  null,           // final double [][] mismatch,    // null or double [12][]
				  threadsMax,     // final int         threadsMax,  // maximal number of threads to launch
				  updateStatus,   // final boolean     updateStatus,
				  debugLevel);    // final int         debugLevel);
	  }

	  public CLTPass3d  CLTMeasure( // perform single pass according to prepared tiles operations and disparity // USED in lwir
			  final CLTParameters clt_parameters,
			  final int           scanIndex,
			  final boolean       save_textures,
			  final boolean       save_corr,
			  final double [][]   mismatch,    // null or double [12][] or [numClusters][] for new LMA
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
// Should not happen !
		  double [][]  disparity_array = scan.disparity;
		  if (scan.disparity == null) { // not used in lwir
			  System.out.println ("** BUG: should not happen - scan.disparity == null ! **");
			  System.out.println ("Trying to recover");
			  double [] backup_disparity = scan.getDisparity(0);
			  if (backup_disparity == null) {
				  System.out.println ("** BUG: no disparity at all !");
				  backup_disparity = new double[tilesX*tilesY];
			  }
			  scan.disparity = new double[tilesY][tilesX];
			  for (int ty = 0; ty < tilesY; ty++) {
				  for (int tx = 0; tx < tilesX; tx++) {
					  scan.disparity[ty][tx] = backup_disparity[ty*tilesX + tx];
					  if (Double.isNaN(scan.disparity[ty][tx])) {
						  scan.disparity[ty][tx] = 0;
						  tile_op[ty][tx] = 0;
					  }
				  }
			  }
			  disparity_array = scan.disparity;

		  }

		  // undecided, so 2 modes of combining alpha - same as rgb, or use center tile only
		  double [][][][]     clt_corr_combo = null; //    new double [ImageDtt.TCORR_TITLES.length][tilesY][tilesX][]; // will only be used inside?
		  // broken clt_aberrations_quad_corr_new
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
		  ImageDtt image_dtt = new ImageDtt(
				  clt_parameters.transform_size,
				  isMonochrome(),
				  isLwir(),
				  clt_parameters.getScaleStrength(isAux()));
		  double z_correction =  clt_parameters.z_correction;
		  if (clt_parameters.z_corr_map.containsKey(image_name)){ // not used in lwir
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
				  tilesX * image_dtt.transform_size, // imp_quad[0].getWidth(),       // final int width,
				  clt_parameters.getFatZero(isMonochrome()),      // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
				  clt_parameters.corr_sym,
				  clt_parameters.corr_offset,
				  clt_parameters.corr_red,
				  clt_parameters.corr_blue,
				  clt_parameters.getCorrSigma(image_dtt.isMonochrome()),
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
///				  image_dtt.transform_size,
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
		  scan.is_measured =   true; // but no disparity map/textures
		  scan.is_combo =      false;
		  scan.resetProcessed();
		  return scan;
	  }


	  public CLTPass3d  CLTMeasure( // perform single pass according to prepared tiles operations and disparity // USED in lwir
			  final CLTParameters           clt_parameters,
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
		  double [][][][]     clt_corr_combo =    null; // new double [ImageDtt.TCORR_TITLES.length][tilesY][tilesX][]; // will only be used inside?
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
		  ImageDtt image_dtt = new ImageDtt(
				  clt_parameters.transform_size,
				  isMonochrome(),
				  isLwir(),
				  clt_parameters.getScaleStrength(isAux()));
		  double z_correction =  clt_parameters.z_correction;
		  if (clt_parameters.z_corr_map.containsKey(image_name)){ // not used in lwir
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
				  tilesX * image_dtt.transform_size, // imp_quad[0].getWidth(),       // final int width,
				  clt_parameters.getFatZero(isMonochrome()),      // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
				  clt_parameters.corr_sym,
				  clt_parameters.corr_offset,
				  clt_parameters.corr_red,
				  clt_parameters.corr_blue,
				  clt_parameters.getCorrSigma(image_dtt.isMonochrome()),
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
///				  image_dtt.transform_size,
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


	  public CLTPass3d  CLTMeasureLY( // perform single pass according to prepared tiles operations and disparity // USED in lwir
			  final CLTParameters clt_parameters,
			  final int           scanIndex,
			  final int           bgIndex, // combine, if >=0
			  final int           threadsMax,  // maximal number of threads to launch
			  final boolean       updateStatus,
			  int           debugLevel)
	  {
		  final int dbg_x = -295-debugLevel;
		  final int dbg_y = -160-debugLevel;
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  final int cluster_size =clt_parameters.tileStep;
		  final int clustersX= (tilesX + cluster_size - 1) / cluster_size;
		  final int clustersY= (tilesY + cluster_size - 1) / cluster_size;

		  CLTPass3d scan = tp.clt_3d_passes.get(scanIndex);
		  scan.setLazyEyeClusterSize(cluster_size);
		  boolean [] force_disparity= new boolean[clustersX * clustersY];
		  //		  scan.setLazyEyeForceDisparity(force_disparity);
		  if (bgIndex >= 0) {
			  CLTPass3d    bg_scan = tp.clt_3d_passes.get(bgIndex);
			  // if at least one tile in a cluster is BG, use BG for the whole cluster and set lazy_eye_force_disparity
			  for (int cY = 0; cY < clustersY; cY ++) {
				  for (int cX = 0; cX < clustersX; cX ++) {
					  boolean has_bg = false;
					  for (int cty = 0; (cty < cluster_size) && !has_bg; cty++) {
						  int ty = cY * cluster_size + cty;
						  if (ty < tilesY) for (int ctx = 0; ctx < cluster_size; ctx++) {
							  int tx = cX * cluster_size + ctx;
							  if ((tx < tilesX ) && (bg_scan.tile_op[ty][tx] > 0)) {
								  has_bg = true;
								  break;
							  }
						  }
					  }
					  if (has_bg) {
						  for (int cty = 0; cty < cluster_size; cty++) {
							  int ty = cY * cluster_size + cty;
							  if (ty < tilesY) for (int ctx = 0; ctx < cluster_size; ctx++) {
								  int tx = cX * cluster_size + ctx;
								  if (tx < tilesX ) {
									  scan.tile_op[ty][tx] =   bg_scan.tile_op[ty][tx];
									  scan.disparity[ty][tx] = bg_scan.disparity[ty][tx];
								  }
							  }
						  }
						  force_disparity[cY * clustersX + cX] = true;
					  }
				  }
			  }
			  scan.setLazyEyeForceDisparity(force_disparity);
		  }
		  int [][]     tile_op =         scan.tile_op;
		  double [][]  disparity_array = scan.disparity;

		  // Should not happen !
		  if (scan.disparity == null) { // not used in lwir
			  System.out.println ("** BUG: should not happen - scan.disparity == null ! **");
			  System.out.println ("Trying to recover");
			  double [] backup_disparity = scan.getDisparity(0);
			  if (backup_disparity == null) {
				  System.out.println ("** BUG: no disparity at all !");
				  backup_disparity = new double[tilesX*tilesY];
			  }
			  scan.disparity = new double[tilesY][tilesX];
			  for (int ty = 0; ty < tilesY; ty++) {
				  for (int tx = 0; tx < tilesX; tx++) {
					  scan.disparity[ty][tx] = backup_disparity[ty*tilesX + tx];
					  if (Double.isNaN(scan.disparity[ty][tx])) {
						  scan.disparity[ty][tx] = 0;
						  tile_op[ty][tx] = 0;
					  }
				  }
			  }
			  disparity_array = scan.disparity;
		  }

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
		  double [][] shiftXY = new double [4][2];
		  if (!clt_parameters.fine_corr_ignore) {
			  double [][] shiftXY0 = {
					  {clt_parameters.fine_corr_x_0,clt_parameters.fine_corr_y_0},
					  {clt_parameters.fine_corr_x_1,clt_parameters.fine_corr_y_1},
					  {clt_parameters.fine_corr_x_2,clt_parameters.fine_corr_y_2},
					  {clt_parameters.fine_corr_x_3,clt_parameters.fine_corr_y_3}};
			  shiftXY = shiftXY0;
		  }

		  ImageDtt image_dtt = new ImageDtt(
				  clt_parameters.transform_size,
				  isMonochrome(),
				  isLwir(),
				  clt_parameters.getScaleStrength(isAux()));
		  double z_correction =  clt_parameters.z_correction;
		  if (clt_parameters.z_corr_map.containsKey(image_name)){ // not used in lwir
			  z_correction +=clt_parameters.z_corr_map.get(image_name);
		  }
		  final double disparity_corr = (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);

		  if (debugLevel > -5){
			  tp.showScan(
					  scan,   // CLTPass3d   scan,
					  "LY-combo_scan-"+scan+"_post"); //String title)
		  }
		  // use new, LMA-based mismatch calculation
		  double [][] lazy_eye_data = image_dtt.cltMeasureLazyEye ( // returns d,s lazy eye parameters
				  clt_parameters.img_dtt,       // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
				  tile_op,                      // per-tile operation bit codes
				  disparity_array,              // clt_parameters.disparity,     // final double            disparity,
				  image_data,                   // final double [][][]      imade_data, // first index - number of image in a quad
				  saturation_imp,               // boolean [][] saturation_imp, // (near) saturated pixels or null
				  null, // final double [][]         clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
				  // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
				  null, // disparity_map,    // [12][tp.tilesY * tp.tilesX]
				  tilesX * image_dtt.transform_size, // imp_quad[0].getWidth(),       // final int width,
				  clt_parameters.getFatZero(isMonochrome()),      // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
				  clt_parameters.corr_red,
				  clt_parameters.corr_blue,
				  clt_parameters.getCorrSigma(image_dtt.isMonochrome()),
				  min_corr_selected, // 0.0001; // minimal correlation value to consider valid
				  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
				  null,                          // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
				  clt_kernels,                   // final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				  clt_parameters.kernel_step,
//				  image_dtt.transform_size,
				  clt_parameters.clt_window,
				  shiftXY, //
				  disparity_corr, // final double              disparity_corr, // disparity at infinity
				  clt_parameters.shift_x,       // final int               shiftX, // shift image horizontally (positive - right) - just for testing
				  clt_parameters.shift_y,       // final int               shiftY, // shift image vertically (positive - down)
				  clt_parameters.tileStep,         // 	final int                 tileStep, // process tileStep x tileStep cluster of tiles when adjusting lazy eye parameters
				  clt_parameters.tileX,         // final int               debug_tileX,
				  clt_parameters.tileY,         // final int               debug_tileY,
				  threadsMax,
				  debugLevel - 2); // -0);
		  scan.setLazyEyeData(lazy_eye_data);
		  scan.is_measured =   true; // but no disparity map/textures
		  scan.is_combo =      false;
		  scan.resetProcessed();
		  return scan;
	  }






	  public ImagePlus [] conditionImageSetBatch( // used in batchCLT3d // not used in lwir
			  final int                           nSet, // index of the 4-image set
			  final CLTParameters           clt_parameters,
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
		  this.geometryCorrection.woi_tops = new int [channelFiles.length];
		  double [][] dbg_dpixels =  batch_mode? null : (new double [channelFiles.length][]);

		  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
			  int nFile=channelFiles[srcChannel];
			  imp_srcs[srcChannel]=null;
			  if (nFile >=0){
				  imp_srcs[srcChannel] = eyesisCorrections.getJp4Tiff(sourceFiles[nFile], this.geometryCorrection.woi_tops);

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
		  System.out.println("Temporarily applying scaleExposures[] here - 3" );
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
			  (new ShowDoubleFloatArrays()).showArrays(dbg_satur, width, height, true, "Saturated" , titles);

			  if (debugLevel > -1) { // 0){
				  double [][] dbg_dpixels_norm = new double [channelFiles.length][];
				  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
					  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
					  dbg_dpixels_norm[srcChannel] = new double[pixels.length];
					  for (int i = 0; i < pixels.length; i++){
						  dbg_dpixels_norm[srcChannel][i] = pixels[i];
					  }
				  }
				  (new ShowDoubleFloatArrays()).showArrays(dbg_dpixels, width, height, true, "dpixels" , titles);
				  (new ShowDoubleFloatArrays()).showArrays(dbg_dpixels_norm, width, height, true, "dpixels_norm" , titles);
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
				  (new ShowDoubleFloatArrays()).showArrays(dbg_dpixels_split, width/2, height/2, true, "dpixels_split" , dbg_titles);
			  }
		  }

		  return imp_srcs;
	  }

	  public void batchCLT3d( // Same can be ran for aux? // not used in lwir
			  TwoQuadCLT       twoQuadCLT, //maybe null in no-rig mode, otherwise may contain rig measurements to be used as infinity ground truth
			  CLTParameters          clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters      debayerParameters,
			  ColorProcParameters                               colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters          channelGainParameters, // also need aux!
			  EyesisCorrectionParameters.RGBParameters          rgbParameters,
			  EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
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
				  int [] channels= fileChannelToSensorChannels(correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile]));
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
				  int [] channels= fileChannelToSensorChannels(correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile]));
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
			  setFiles.get(setNames.indexOf(setName)).add(nFile); // .add(new Integer(nFile));
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
						  colorProcParameters,
						  rgbParameters,
						  threadsMax,  // maximal number of threads to launch
						  updateStatus,
						  debugLevelInner);
				  if (ok) {
					  System.out.println("Adjusting extrinsics");
					  extrinsicsCLT(
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
						  colorProcParameters,
						  rgbParameters,
						  threadsMax,  // maximal number of threads to launch
						  updateStatus,
						  debugLevelInner);
				  if (ok) {
					  System.out.println("Adjusting polynomial fine crorection");
					  extrinsicsCLT(
							  clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
							  true,           // adjust_poly,
							  threadsMax,     //final int        threadsMax,  // maximal number of threads to launch
							  updateStatus,   // final boolean    updateStatus,
							  debugLevelInner);    // final int        debugLevel)
				  }

			  }
			  if (correctionsParameters.clt_batch_4img){ // not used in lwir
				  processCLTQuadCorrCPU( // returns ImagePlus, but it already should be saved/shown
						  imp_srcs, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
						  saturation_imp, // boolean [][] saturation_imp, // (near) saturated pixels or null
						  clt_parameters,
						  debayerParameters,
						  colorProcParameters,
						  channelGainParameters,
						  rgbParameters,
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
						  colorProcParameters,
						  rgbParameters,
						  threadsMax,  // maximal number of threads to launch
						  updateStatus,
						  debugLevelInner);
				  if (ok) {
					  System.out.println("Explore 3d space");
					  expandCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
							  clt_parameters,
							  debayerParameters,
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
					// prepare average RGBA for the last scan
					setPassAvgRBGA(                      // get image from a single pass, return relative path for x3d // USED in lwir
							clt_parameters,                           // CLTParameters           clt_parameters,
							tp.clt_3d_passes.size() - 1, // int        scanIndex,
							threadsMax,                               // int        threadsMax,  // maximal number of threads to launch
							updateStatus,                             // boolean    updateStatus,
							debugLevelInner);                         // int        debugLevel)
				  
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
						  debugLevelInner);         // final int        debugLevel)
				  if (!ok) continue;
			  } else continue; // if (correctionsParameters.clt_batch_gen3d)

			  Runtime.getRuntime().gc();
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

	  public boolean setGpsLla( // USED in lwir
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
		  return false; // not used in lwir
	  }


	  public boolean writeKml( // USED in lwir
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

	  public boolean createThumbNailImage( // USED in lwir
			  ImagePlus imp,
			  String dir,
			  String name,
			  int debugLevel)
	  {
		  String thumb_path = dir +  Prefs.getFileSeparator() + name+".jpeg";
		  if (new File(thumb_path).exists() && !correctionsParameters.thumb_overwrite) {
			  System.out.println("file "+thumb_path+" exists, skipping thumbnail generation");
			  return false;
		  }

		  int image_width = imp.getWidth();
		  int image_height = imp.getHeight();
		  ImageProcessor ip = imp.getProcessor().duplicate();
		  if ((image_width >=  correctionsParameters.thumb_width) &&
				  (image_height >=  correctionsParameters.thumb_height)) {
			  double scale_h = 1.0 * (correctionsParameters.thumb_width + 1)/image_width;
			  double scale_v = 1.0 * (correctionsParameters.thumb_height + 1)/image_height;
			  double scale = ((scale_h > scale_v) ? scale_h : scale_v) / correctionsParameters.thumb_size;


			  ip.setInterpolationMethod(ImageProcessor.BICUBIC);
			  if (!isMonochrome()) {
				  ip.blurGaussian(2.0);
			  }
			  ip.scale(scale, scale);
			  int lm = (int) Math.round (((image_width*scale)-correctionsParameters.thumb_width)* correctionsParameters.thumb_h_center + (0.5*image_width*(1.0-scale)));
			  int tm = (int) Math.round (((image_height*scale)-correctionsParameters.thumb_height)* correctionsParameters.thumb_v_center + (0.5*image_height*(1.0-scale)));
			  Rectangle r = new Rectangle(lm,tm,correctionsParameters.thumb_width,correctionsParameters.thumb_height);
			  ip.setRoi(r);
			  ip = ip.crop();
		  } else {

		  }
		  ImagePlus ip_thumb = new ImagePlus(name,ip);
		  eyesisCorrections.saveAndShow(
				  ip_thumb,
				  dir,
				  false,
				  false,
				  correctionsParameters.JPEG_quality, // jpegQuality); // jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
				  (debugLevel > -2) ? debugLevel : 1); // int debugLevel (print what it saves)

		  return true;
	  }



	  public boolean writeRatingFile( // USED in lwir
			  int debugLevel
			  )
	  {
		  String set_name = image_name;
		  if (set_name == null ) {
			  QuadCLTCPU.SetChannels [] set_channels = setChannels(debugLevel);
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
	  
	  public void setPassAvgRBGA( // get image from a single pass, return relative path for x3d // USED in lwir
			  CLTParameters           clt_parameters,
			  int        scanIndex,
			  int        threadsMax,  // maximal number of threads to launch
			  boolean    updateStatus,
			  int        debugLevel)
	  {
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
//		  final int transform_size =clt_parameters.transform_size;
		  CLTPass3d scan = tp.clt_3d_passes.get(scanIndex);
		  double [][][][] texture_tiles = scan.texture_tiles; 
		  if (texture_tiles == null) return;
		  int num_layers = 0;
		  for (int ty = 0; ty < tilesY; ty++){
			  if (texture_tiles[ty] != null){
				  for (int tx = 0; tx < tilesX; tx++){
					  if (texture_tiles[ty][tx] != null){
						  num_layers = texture_tiles[ty][tx].length;
						  break;
					  }
				  }
				  if (num_layers > 0) break;
			  }
			  if (num_layers > 0) break;
		  }
		  int numTiles = tilesX * tilesY;
		  double [] scales = new double [num_layers];
		  for (int n = 0; n < num_layers; n++){
			  if       (n < 3)  scales[n] = 1.0/255.0; // R,B,G
			  else if  (n == 3) scales[n] = 1.0; //alpha
			  else if  (n < 8)  scales[n] = 1.0; // ports 0..3
			  else              scales[n] = 1.0/255.0; // RBG rms, in 1/255 units, but small
		  }
		  double [][] tileTones = new double [num_layers][numTiles];
		  for (int ty = 0; ty < tilesY; ty++ ) if (texture_tiles[ty] != null){
			  for (int tx = 0; tx < tilesX; tx++ ) if (texture_tiles[ty][tx] != null) {
				  int indx = ty * tilesX + tx;
				  for (int n = 0; n < num_layers; n++) if (texture_tiles[ty][tx][n] != null){
					  double s = 0.0;
					  for (int i = 0; i < texture_tiles[ty][tx][n].length; i++){
						  s += texture_tiles[ty][tx][n][i];
					  }
					  s /= (texture_tiles[ty][tx][n].length/4); // overlapping tiles
					  s *= scales[n];
					  tileTones[n][indx] = s;
				  }
			  }
		  }
		  scan.setTilesRBGA(tileTones); 
	  }
	  // non-trivial in QuadCLT (for the GPU)
	  public CLTPass3d  CLTMeasureCorrTesting( // perform single pass according to prepared tiles operations and disparity // not used in lwir
			  CLTParameters     clt_parameters,
			  final int         scanIndex,
			  final boolean     save_textures,
			  final int         threadsMax,  // maximal number of threads to launch
			  final boolean     updateStatus,
			  final int         debugLevel)
	  {
		 return CLTMeasureCorrTesting( // perform single pass according to prepared tiles operations and disparity // not used in lwir
				  clt_parameters,
				  scanIndex,
				  save_textures,
				  threadsMax,  // maximal number of threads to launch
				  updateStatus,
				  debugLevel); 
	  }
}
