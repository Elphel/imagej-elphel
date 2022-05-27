package com.elphel.imagej.cameras;
/**
** -----------------------------------------------------------------------------**
** EyesisCorrectionParameters.java
**
** Parameter classes for aberration correction for Eyesis4pi
**
**
** Copyright (C) 2012-2018 Elphel, Inc.
**
** -----------------------------------------------------------------------------**
**
**  EyesisCorrectionParameters.java is free software: you can redistribute it and/or modify
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

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Properties;

import javax.swing.JFileChooser;

import com.elphel.imagej.calibration.CalibrationFileManagement;
import com.elphel.imagej.calibration.DirectoryChoser;
import com.elphel.imagej.calibration.MultipleExtensionsFileFilter;
import com.elphel.imagej.common.GenericJTabbedDialog;
import com.elphel.imagej.common.WindowTools;

import ij.IJ;
import ij.Prefs;
import ij.gui.GenericDialog;


public class EyesisCorrectionParameters {
    public static class CorrectionParameters{
    	public static final String AUX_PREFIX = "AUX-";
    	public boolean swapSubchannels01=      true; // false; // (false: 0-1-2, true - 1-0-2)
  		public boolean split=                  true;
  		public boolean vignetting=             true;
  		public boolean pixelDefects=           true;
  		public double  pixelDefectsThreshold=  8.0; // normally none with less than 5.0 are stored?
  		public boolean debayer=                true;
  		public boolean showDebayerEnergy =     false;
  		public boolean saveDebayerEnergy =     true;
  		public boolean deconvolve =            true;
  		public boolean combine =               true;
  		public boolean showDenoiseMask =       false;
  		public boolean saveDenoiseMask =       true;
  		public boolean showChromaDenoiseMask = false;
  		public boolean saveChromaDenoiseMask = true;
  		public boolean showNoiseGains =        false;
  		public boolean saveNoiseGains =        false;
  		public boolean colorProc =             true;
  		public boolean blueProc =              true;
  		public boolean toRGB =                 true;
  		public boolean rotate =                true;
  		public boolean crop =                  true;  // crop to the sensor size
  		public int     equirectangularFormat=     0;  // 0 - 8 bit RGBA, 1 - 16 bit RGBA, 2 (32 int or 16 float!) ?, 3 - 32-bit FP RGBA. only 0, 1 and 3 currently supported
  		public double  outputRangeInt=          0.25;  // 1.0 intensity will be mapped to 65535*0.25
  		public double  outputRangeFP=          255.0; // 1.0 intensity will be saved as 255.0 (in float 32-bit mode)
  		public boolean imageJTags=             false; // encode ImageJ info data to the TIFF output header

  		public boolean jpeg =                  true;  // convert to RGB and save JPEG (if save is true)
  		public boolean png =                   true;  // use PNG instead of TIFF for 32-bit ARGB
  		public boolean save =                  true;
  		public boolean save16 =                false; // save 16-bit tiff also if the end result is 8 bit
  		public boolean save32 =                false; // save 32-bit tiff also if the end result is 8 or 16 bit
  		public boolean show =                  false ;
  		public int     JPEG_quality =          95;
  		public double  JPEG_scale =            0.5;
  		public boolean equirectangular=        true;
  		public boolean zcorrect=               true;
  		public boolean saveSettings =          true;
  		public String  tile_processor_gpu =    ""; // absolute path to tile_processor_gpu project or empty to use default GPU kernels

    	public String [] sourcePaths=          {};
//    	public String [] sourceSetPaths=       {}; // 2019 - directories with image sets
    	public boolean   use_set_dirs =        false; // each image set in a directory, use directory as a timestamp
    	public String sourceDirectory=         "";
    	public String sourcePrefix=            "";
    	public String sourceSuffix=            ".tiff"; //".jp4"
    	// first subcamera index as in properties of the sensor configuration and kernels and CLT kernels
    	public int    firstSubCameraConfig=    0; // channel index in config (sensor, clt) files
    	// first filename index to be processed by this PixelMapping instance
    	// only one instance for Eyesis, two for 8-rig and Talon (4EO+4LWIR)
    	public int    firstSubCamera=          1; // channel index in source file names
    	// number of subcameras in this PixelMapping instance
    	// number of source file indices and image channels (including mux) are found from the number
    	// of PixelMapping instance
    	public int    numSubCameras=           4; // channel index in source file names
    	public String sensorDirectory=         "";
    	public String sensorPrefix=            "sensor-";
    	public String sensorSuffix=            ".calib-tiff"; // fixed in PixelMapping

    	public String sharpKernelDirectory=    "";
    	public String sharpKernelPrefix=       "sharpKernel-";
    	public String sharpKernelSuffix=       ".kernel-tiff";
    	public String smoothKernelDirectory=   "";
    	public String smoothKernelPrefix=      "smoothKernel-";
    	public String smoothKernelSuffix=      ".kernel-tiff";
    	public String dctKernelDirectory=      "";
    	public String dctKernelPrefix=         "dct-";
    	public String dctSymSuffix=            ".sym-tiff";
    	public String dctAsymSuffix=           ".asym-tiff";
    	public String equirectangularDirectory="";
    	public String equirectangularPrefix=   "";
    	public String equirectangularSuffix=   ".eqr-tiff";
    	public boolean equirectangularCut=     true;
    	public String planeMapPrefix=          "";
    	public String planeMapSuffix=          ".plane-proj-tiff";
    	public boolean usePlaneProjection=     false;  //
    	public boolean planeAsJPEG=            true;   // save de-warped image as JPEG (only if equirectangularFormat==0)
    	public String resultsDirectory=        "";
    	public boolean removeUnusedSensorData= false;
    	public int exposureCorrectionMode=     2; // - 0 - none, 1 - absolute, 2 - relative
    	public double referenceExposure=       0.0003; // 3/10000 sec, used in absolute mode only
    	public double relativeExposure=        0.5; // 0.0 - use shortest (darken), 1.0 - use longest (brighten)

    	public String cltKernelDirectory=      "";
    	public String cltKernelPrefix=         "clt-";
    	public String cltSuffix=               ".clt-tiff";
  		public boolean use_x3d_subdirs =       true;
    	public String x3dSubdirPrefix=         "";
    	public String x3dSubdirSuffix=         "";

  		// CLT 3d batch parameters
    	public boolean process_main_sources = false; // not yet used
    	public boolean process_aux_sources =  true; // not yet used
    	public int     kml_sensors=                0xffffffff; // all sensors 
    	
    	public int     rig_batch_adjust_main = 0;
    	public int     rig_batch_adjust_aux =  0;
    	public int     rig_batch_adjust_rig =  0;

    	
    	public int     rig_batch_adjust_main_gt = 0; // adjust main camera using rig disparity as ground truth
    	public int     rig_batch_adjust_aux_gt =  0; // adjust aux camera using rig disparity as ground truth (TODO: finish geometry in derivatives)
    	public int     rig_batch_adjust_rig_gt =  0; // adjust rig after main and aux are adjusted with rig GT (late rig adjustment)

    	public boolean clt_batch_dsi1 =       true; // experimental for interscene
  		public boolean clt_batch_apply_man =  false;  // Apply (and disable) manual pixel shift
  		public boolean clt_batch_extrinsic =  false; // Calibrate extrinsic parameters for each set
  		public boolean clt_batch_poly =       false; // Calculate fine polynomial correction for each set
  		public boolean clt_batch_4img =       true;  // Create a set of 4 images, usually for disparity = 0
  		public boolean clt_batch_4img_aux =   true;  // Create a set of 4 images, usually for disparity = 0 for AUX camera
  		public boolean clt_batch_explore =    true;  // 1-st step of 3d reconstruction - explore disparities for each tile
  		public boolean clt_batch_surf =       true;  // Create super-tile 2.5d surfaces
  		public boolean clt_batch_assign =     true;  // Assign tiles to surfaces
  		public boolean clt_batch_gen3d =      true;  // Generate 3d output: x3d and/or obj+mtl
  		public boolean clt_batch_genMl =      true;  // Generate ML output
  		public boolean clt_batch_dbg1 =       true;  // Generate debug images if a single set is selected
  		public boolean clt_batch_dsi =        true;  // Create and save DSI combo image with the model
		public boolean clt_batch_dsi_aux =    false;  // Calculate and save aux camera DSI (currently it is offset from the main/rig data
		public boolean clt_batch_dsi_cm_strength = true;  // Use CM strength (no switch between LMA/no-LMA) for DSI export
		public boolean clt_batch_dsi_aux_full=false;  // more than just preExpandCLTQuad3d() (same as for Lazy Eye
  		public boolean clt_batch_save_extrinsics = true;  // Save cameras extrinsic parameters with the model
  		public boolean clt_batch_save_all =   true;  // Save all parameters with the model

  		public boolean clt_batch_skip_scenes =     false;  // Skip all per-scene processing, go directly to processing sequences
  		
  		public boolean clt_batch_pose_pairs_main = false;  // calculate pair-wise camera poses
  		public boolean clt_batch_pose_last_main =  false;  // calculate camera poses realtive to the last scene
  		public boolean clt_batch_pose_scene_main = false;  // calculate camera poses relative to all other ones
  		public int     clt_batch_offset_main =     0;      // when selecting multiple reference scene, offset from the last one
  		public int     clt_batch_step_main   =     3;      // step (decreasing timestamp) to select reference frames in a sequence
  		public boolean clt_batch_ml_last_main =    false;  // export ML files for the last (reference scene)
  		public boolean clt_batch_ml_all_main =     false;  // export ML files for all available reference scenes
  		
  		public boolean clt_batch_pose_pairs_aux =  false;  // calculate pair-wise camera poses
  		public boolean clt_batch_pose_last_aux =   false;  // calculate camera poses realtive to the last scene
  		public boolean clt_batch_pose_scene_aux =  false;  // calculate camera poses relative to all other ones
  		public int     clt_batch_offset_aux =      0;      // when selecting multiple reference scene, offset from the last one
  		public int     clt_batch_step_aux   =     10;      // step (decreasing timestamp) to select reference frames in a sequence
  		public boolean clt_batch_ml_last_aux =     false;  // export ML files for the last (reference scene)
  		public boolean clt_batch_ml_all_aux =      false;  // export ML files for all available reference scenes


    	public String x3dModelVersion="v01";
    	public String jp4SubDir="jp4"; // FIXME:

    	public String x3dDirectory="";

    	public String mlDirectory="ml";

    	public boolean  thumb_overwrite =     true;
    	public int      thumb_width =         200;
    	public int      thumb_height =        100;
    	public double   thumb_h_center =      0.5;
    	public double   thumb_v_center =      0.5;
    	public double   thumb_size =          0.75;
    	public int      default_rating =      5;



    	public CorrectionParameters getAux() {
    		return aux_camera;
    	}
    	public CorrectionParameters aux_camera = null; // auxiliary camera parameters
//  		public boolean use_aux =             true;  // Generate debug images if a single set is selected
		public void updateAuxFromMain() { // from master to aux
			if (aux_camera == null) {
				aux_camera = new CorrectionParameters();
				initAuxFromMain(aux_camera);
			} else {
				updateAuxFromMain(aux_camera);
			}
		}
		public void updateAuxFromMain(CorrectionParameters cp) { // from master to aux
  			cp.split =                  this.split;
  			cp.vignetting=              this.vignetting;
  			cp.pixelDefects=            this.pixelDefects;
  			cp.pixelDefectsThreshold=   this.pixelDefectsThreshold;
  			cp.debayer=                 this.debayer;
  			cp.showDebayerEnergy =      this.showDebayerEnergy;
  			cp.saveDebayerEnergy=  	    this.saveDebayerEnergy;
  			cp.deconvolve=              this.deconvolve;
  			cp.combine=                 this.combine;
  			cp.showDenoiseMask=  	    this.showDenoiseMask;
  			cp.saveDenoiseMask=         this.saveDenoiseMask;
  			cp.showChromaDenoiseMask=   this.showChromaDenoiseMask;
  			cp.saveChromaDenoiseMask=   this.saveChromaDenoiseMask;
  			cp.showNoiseGains=          this.showNoiseGains;
  			cp.saveNoiseGains=  	    this.saveNoiseGains;
  			cp.colorProc=  			    this.colorProc;
  			cp.blueProc=  			    this.blueProc;
  			cp.toRGB=  			        this.toRGB;
  			cp.rotate=  		        this.rotate;
  			cp.crop=  			        this.crop;
  			cp.equirectangularFormat=   this.equirectangularFormat;
  			cp.outputRangeInt=          this.outputRangeInt;
  			cp.outputRangeFP=  		    this.outputRangeFP;
  			cp.imageJTags=  		    this.imageJTags;
  			cp.jpeg=  			        this.jpeg;
  			cp.png=  			        this.png;
  			cp.save=  			        this.save;
  			cp.save16=  			    this.save16;
  			cp.save32=  			    this.save32;
  			cp.show=  			        this.show;
  			cp.JPEG_quality=            this.JPEG_quality;
  			cp.JPEG_scale=  		    this.JPEG_scale;
  			cp.equirectangular=  	    this.equirectangular;
  			cp.zcorrect=  			    this.zcorrect;
  			cp.saveSettings=  		    this.saveSettings;
  			cp.sourceDirectory=    	    this.sourceDirectory;
  			cp.tile_processor_gpu =     this.tile_processor_gpu;
  			cp.use_set_dirs =           this.use_set_dirs;
//  			cp.sourcePrefix=    	    this.sourcePrefix;
//  			cp.sourceSuffix=    	    this.sourceSuffix;
//  			cp.firstSubCamera=    	    this.firstSubCamera;
//  			cp.numSubCameras=    	    this.numSubCameras;
//  			cp.sensorDirectory=         this.sensorDirectory;
//  			cp.sensorPrefix=            this.sensorPrefix;
//  			cp.sensorSuffix=      	    this.sensorSuffix;
  			cp.sharpKernelDirectory=    this.sharpKernelDirectory;
  			cp.sharpKernelPrefix=       this.sharpKernelPrefix;
  			cp.sharpKernelSuffix=       this.sharpKernelSuffix;
  			cp.smoothKernelDirectory=   this.smoothKernelDirectory;
  			cp.smoothKernelPrefix=      this.smoothKernelPrefix;
  			cp.smoothKernelSuffix=      this.smoothKernelSuffix;
  			cp.dctKernelDirectory=      this.dctKernelDirectory;
  			cp.dctKernelPrefix=    	    this.dctKernelPrefix;
  			cp.dctSymSuffix=    	    this.dctSymSuffix;
  			cp.dctAsymSuffix=    	    this.dctAsymSuffix;
  			cp.equirectangularDirectory=this.equirectangularDirectory;
  			cp.equirectangularPrefix=   this.equirectangularPrefix;
  			cp.equirectangularSuffix=   this.equirectangularSuffix;
  			cp.equirectangularCut=      this.equirectangularCut;
  			cp.planeMapPrefix=    		this.planeMapPrefix;
  			cp.planeMapSuffix=          this.planeMapSuffix;
  			cp.usePlaneProjection=    	this.usePlaneProjection;
  			cp.planeAsJPEG=    		    this.planeAsJPEG;
//  			cp.resultsDirectory=    	this.resultsDirectory;
  			cp.removeUnusedSensorData=  this.removeUnusedSensorData;
    		if (this.sourcePaths!=null) {
    			cp.sourcePaths=new String[this.sourcePaths.length];
        		for (int i=0;i<this.sourcePaths.length;i++){
        			cp.sourcePaths[i] = this.sourcePaths[i];
        		}
    		}
  			cp.exposureCorrectionMode=  this.exposureCorrectionMode;
  			cp.referenceExposure=    	this.referenceExposure;
  			cp.relativeExposure=    	this.relativeExposure;
  			cp.swapSubchannels01=    	this.swapSubchannels01;
  			cp.x3dDirectory=    		this.x3dDirectory;
  			cp.mlDirectory=    	     	this.mlDirectory;
  			cp.use_x3d_subdirs=    		this.use_x3d_subdirs;
  			cp.x3dSubdirPrefix=    		this.x3dSubdirPrefix;
  			cp.x3dModelVersion=    		this.x3dModelVersion;
  			cp.jp4SubDir=    	     	this.jp4SubDir;


  			cp.process_main_sources=    this.process_main_sources;
  			cp.process_aux_sources=     this.process_aux_sources;
  			
  			cp.kml_sensors=             this.kml_sensors;
  			cp.rig_batch_adjust_main=   this.rig_batch_adjust_main;
  			cp.rig_batch_adjust_aux=	this.rig_batch_adjust_aux;
  			cp.rig_batch_adjust_rig=	this.rig_batch_adjust_rig;

  			cp.rig_batch_adjust_main_gt = this.rig_batch_adjust_main_gt;
  			cp.rig_batch_adjust_aux_gt =  this.rig_batch_adjust_aux_gt;
  			cp.rig_batch_adjust_rig_gt =  this.rig_batch_adjust_rig_gt;

  			cp.clt_batch_dsi1=		    this.clt_batch_dsi1;
  			cp.clt_batch_apply_man=		this.clt_batch_apply_man;
  			cp.clt_batch_extrinsic=		this.clt_batch_extrinsic;
  			cp.clt_batch_poly=    		this.clt_batch_poly;
  			cp.clt_batch_4img=    		this.clt_batch_4img;
  			cp.clt_batch_4img_aux=      this.clt_batch_4img_aux;
  			cp.clt_batch_explore=       this.clt_batch_explore;
  			cp.clt_batch_surf=    		this.clt_batch_surf;
  			cp.clt_batch_assign=    	this.clt_batch_assign;
  			cp.clt_batch_gen3d=    		this.clt_batch_gen3d;
  			cp.clt_batch_genMl=    		this.clt_batch_genMl;
  			cp.clt_batch_dbg1=    		this.clt_batch_dbg1;

  			cp.clt_batch_dsi=    		  this.clt_batch_dsi;
  			cp.clt_batch_dsi_aux=    	  this.clt_batch_dsi_aux;
  			cp.clt_batch_dsi_cm_strength= this.clt_batch_dsi_cm_strength;
  			cp.clt_batch_dsi_aux_full= 	  this.clt_batch_dsi_aux_full;
  			cp.clt_batch_save_extrinsics= this.clt_batch_save_extrinsics;
  			cp.clt_batch_save_all=        this.clt_batch_save_all;

  			cp.clt_batch_skip_scenes=        this.clt_batch_skip_scenes;

  			cp.clt_batch_pose_pairs_main=    this.clt_batch_pose_pairs_main;
  			cp.clt_batch_pose_last_main=     this.clt_batch_pose_last_main;
  			cp.clt_batch_pose_scene_main=    this.clt_batch_pose_scene_main;
  			
  			cp.clt_batch_offset_main=        this.clt_batch_offset_main;
  			cp.clt_batch_step_main=          this.clt_batch_step_main;
  			
  			cp.clt_batch_ml_last_main=       this.clt_batch_ml_last_main;
  			cp.clt_batch_ml_all_main=        this.clt_batch_ml_all_main;
  			
  			
  			cp.clt_batch_pose_pairs_aux=    this.clt_batch_pose_pairs_aux;
  			cp.clt_batch_pose_last_aux=     this.clt_batch_pose_last_aux;
  			cp.clt_batch_pose_scene_aux=    this.clt_batch_pose_scene_aux;
  			cp.clt_batch_offset_aux=        this.clt_batch_offset_aux;
  			cp.clt_batch_step_aux=          this.clt_batch_step_aux;
  			cp.clt_batch_ml_last_aux=       this.clt_batch_ml_last_aux;
  			cp.clt_batch_ml_all_aux=        this.clt_batch_ml_all_aux;

    		cp.thumb_overwrite =        this.thumb_overwrite;
    		cp.thumb_width =            this.thumb_width;
    		cp.thumb_height =           this.thumb_height;
    		cp.thumb_h_center =         this.thumb_h_center;
    		cp.thumb_v_center =         this.thumb_v_center;
    		cp.thumb_size =             this.thumb_size;
    		cp.default_rating =         this.default_rating;


		}


		public void initAuxFromMain(CorrectionParameters cp) { // from master to aux
			updateAuxFromMain(cp); // common parameters
			// empty to prevent accidental use of the wrong kernels/sesnor calibration files
  			cp.sensorDirectory=         ""; // this.sensorDirectory;
  			cp.cltKernelDirectory=    	""; // this.cltKernelDirectory;
  			cp.resultsDirectory=    	this.resultsDirectory+"/aux";
  			cp.firstSubCamera=    	    this.firstSubCamera + this.numSubCameras;
  			cp.firstSubCameraConfig=    this.firstSubCameraConfig; //  + this.numSubCameras; so old setups will have it zero
  			cp.numSubCameras=    	    this.numSubCameras;
  			cp.sensorPrefix=            ""; // this.sensorPrefix;
  			cp.sensorSuffix=      	    this.sensorSuffix;
  			cp.sourcePrefix=      	    this.sourcePrefix;
  			cp.sourceSuffix=      	    this.sourceSuffix;
  			cp.cltKernelPrefix=    		this.cltKernelPrefix;
  			cp.cltSuffix=               this.cltSuffix;
  			cp.x3dSubdirSuffix=    		this.x3dSubdirSuffix+"-aux";

    	}

		public void auxFromExternal(CorrectionParameters ecp) { // from master to aux
			this.aux_camera.sensorDirectory=      ecp.sensorDirectory;
			this.aux_camera.cltKernelDirectory=   ecp.cltKernelDirectory;
			this.aux_camera.resultsDirectory=     ecp.resultsDirectory+"/aux";
			this.aux_camera.firstSubCamera=       ecp.firstSubCamera;
			this.aux_camera.firstSubCameraConfig= ecp.firstSubCameraConfig;
			this.aux_camera.numSubCameras=    	  ecp.numSubCameras;
			this.aux_camera.sourcePrefix=         ecp.sourcePrefix;
			this.aux_camera.sourceSuffix=      	  ecp.sourceSuffix;
			this.aux_camera.sensorPrefix=         ecp.sensorPrefix;
			this.aux_camera.sensorSuffix=      	  ecp.sensorSuffix;
			this.aux_camera.cltKernelPrefix=      ecp.cltKernelPrefix;
			this.aux_camera.cltSuffix=            ecp.cltSuffix;
			this.aux_camera.x3dSubdirSuffix=      ecp.x3dSubdirSuffix + "-aux";
    	}


    	public void setProperties(String prefix,Properties properties){
  			properties.setProperty(prefix+"split",this.split+"");
  			properties.setProperty(prefix+"vignetting",this.vignetting+"");
  			properties.setProperty(prefix+"pixelDefects",this.pixelDefects+"");
  			properties.setProperty(prefix+"pixelDefectsThreshold",this.pixelDefectsThreshold+"");
  			properties.setProperty(prefix+"debayer",this.debayer+"");
  			properties.setProperty(prefix+"showDebayerEnergy",this.showDebayerEnergy+"");
  			properties.setProperty(prefix+"saveDebayerEnergy",this.saveDebayerEnergy+"");
  			properties.setProperty(prefix+"deconvolve",this.deconvolve+"");
  			properties.setProperty(prefix+"combine",this.combine+"");
  			properties.setProperty(prefix+"showDenoiseMask",this.showDenoiseMask+"");
  			properties.setProperty(prefix+"saveDenoiseMask",this.saveDenoiseMask+"");
  			properties.setProperty(prefix+"showChromaDenoiseMask",this.showChromaDenoiseMask+"");
  			properties.setProperty(prefix+"saveChromaDenoiseMask",this.saveChromaDenoiseMask+"");
  			properties.setProperty(prefix+"showNoiseGains",this.showNoiseGains+"");
  			properties.setProperty(prefix+"saveNoiseGains",this.saveNoiseGains+"");
  			properties.setProperty(prefix+"colorProc",this.colorProc+"");
  			properties.setProperty(prefix+"blueProc",this.blueProc+"");
  			properties.setProperty(prefix+"toRGB",this.toRGB+"");
  			properties.setProperty(prefix+"rotate",this.rotate+"");
  			properties.setProperty(prefix+"crop",this.crop+"");
  			properties.setProperty(prefix+"equirectangularFormat",this.equirectangularFormat+"");
  			properties.setProperty(prefix+"outputRangeInt",this.outputRangeInt+"");
  			properties.setProperty(prefix+"outputRangeFP",this.outputRangeFP+"");
  			properties.setProperty(prefix+"imageJTags",this.imageJTags+"");
  			properties.setProperty(prefix+"jpeg",this.jpeg+"");
  			properties.setProperty(prefix+"png",this.png+"");
  			properties.setProperty(prefix+"save",this.save+"");
  			properties.setProperty(prefix+"save16",this.save16+"");
  			properties.setProperty(prefix+"save32",this.save32+"");
  			properties.setProperty(prefix+"show",this.show+"");
  			properties.setProperty(prefix+"JPEG_quality",this.JPEG_quality+"");
  			properties.setProperty(prefix+"JPEG_scale",this.JPEG_scale+"");
  			properties.setProperty(prefix+"equirectangular",this.equirectangular+"");
  			properties.setProperty(prefix+"zcorrect",this.zcorrect+"");
  			properties.setProperty(prefix+"saveSettings",this.saveSettings+"");

    		properties.setProperty(prefix+"sourceDirectory",this.sourceDirectory);
    		properties.setProperty(prefix+"tile_processor_gpu",this.tile_processor_gpu);

    		properties.setProperty(prefix+"use_set_dirs",   this.use_set_dirs+"");

    		properties.setProperty(prefix+"sourcePrefix",this.sourcePrefix);
    		properties.setProperty(prefix+"sourceSuffix",this.sourceSuffix);
    		properties.setProperty(prefix+"firstSubCamera",this.firstSubCamera+"");
    		properties.setProperty(prefix+"firstSubCameraConfig",this.firstSubCameraConfig+"");
    		properties.setProperty(prefix+"numSubCameras", this.numSubCameras+"");

    		properties.setProperty(prefix+"sensorDirectory",this.sensorDirectory);
    		properties.setProperty(prefix+"sensorPrefix",this.sensorPrefix);
    		properties.setProperty(prefix+"sensorSuffix",this.sensorSuffix);


    		properties.setProperty(prefix+"sharpKernelDirectory",this.sharpKernelDirectory);
    		properties.setProperty(prefix+"sharpKernelPrefix",this.sharpKernelPrefix);
    		properties.setProperty(prefix+"sharpKernelSuffix",this.sharpKernelSuffix);
    		properties.setProperty(prefix+"smoothKernelDirectory",this.smoothKernelDirectory);
    		properties.setProperty(prefix+"smoothKernelPrefix",this.smoothKernelPrefix);
    		properties.setProperty(prefix+"smoothKernelSuffix",this.smoothKernelSuffix);

    		properties.setProperty(prefix+"dctKernelDirectory",this.dctKernelDirectory);
    		properties.setProperty(prefix+"dctKernelPrefix",this.dctKernelPrefix);
    		properties.setProperty(prefix+"dctSymSuffix",this.dctSymSuffix);
    		properties.setProperty(prefix+"dctAsymSuffix",this.dctAsymSuffix);

    		properties.setProperty(prefix+"equirectangularDirectory",this.equirectangularDirectory);
    		properties.setProperty(prefix+"equirectangularPrefix",this.equirectangularPrefix);
    		properties.setProperty(prefix+"equirectangularSuffix",this.equirectangularSuffix);
    		properties.setProperty(prefix+"equirectangularCut",this.equirectangularCut+"");

    		properties.setProperty(prefix+"planeMapPrefix",this.planeMapPrefix+"");
    		properties.setProperty(prefix+"planeMapSuffix",this.planeMapSuffix+"");
    		properties.setProperty(prefix+"usePlaneProjection",this.usePlaneProjection+"");
    		properties.setProperty(prefix+"planeAsJPEG",this.planeAsJPEG+"");

    		properties.setProperty(prefix+"resultsDirectory",this.resultsDirectory);
    		properties.setProperty(prefix+"removeUnusedSensorData",this.removeUnusedSensorData+"");
    		if (this.sourcePaths!=null) {
        		properties.setProperty(prefix+"sourcePaths",this.sourcePaths.length+"");
        		for (int i=0;i<this.sourcePaths.length;i++){
        			properties.setProperty(prefix+"sourcePath"+i,this.sourcePaths[i]);
        		}
    		}

    		properties.setProperty(prefix+"exposureCorrectionMode",this.exposureCorrectionMode+"");
    		properties.setProperty(prefix+"referenceExposure",     this.referenceExposure+"");
    		properties.setProperty(prefix+"relativeExposure",      this.relativeExposure+"");
    		properties.setProperty(prefix+"swapSubchannels01",     this.swapSubchannels01+"");

    		properties.setProperty(prefix+"cltKernelDirectory",    this.cltKernelDirectory);
    		properties.setProperty(prefix+"cltKernelPrefix",       this.cltKernelPrefix);
    		properties.setProperty(prefix+"cltSuffix",             this.cltSuffix);

    		properties.setProperty(prefix+"x3dDirectory",          this.x3dDirectory);
    		properties.setProperty(prefix+"use_x3d_subdirs",       this.use_x3d_subdirs+"");

    		properties.setProperty(prefix+"x3dSubdirPrefix",       this.x3dSubdirPrefix+"");
    		properties.setProperty(prefix+"x3dSubdirSuffix",       this.x3dSubdirSuffix+"");
    		properties.setProperty(prefix+"x3dModelVersion",       this.x3dModelVersion);
    		properties.setProperty(prefix+"jp4SubDir",             this.jp4SubDir);

    		properties.setProperty(prefix+"mlDirectory",           this.mlDirectory);

    		properties.setProperty(prefix+"process_main_sources",  this.process_main_sources+"");
    		properties.setProperty(prefix+"process_aux_sources",   this.process_aux_sources+"");

    		properties.setProperty(prefix+"kml_sensors",           this.kml_sensors+"");

    		properties.setProperty(prefix+"rig_batch_adjust_main", this.rig_batch_adjust_main+"");
    		properties.setProperty(prefix+"rig_batch_adjust_aux",  this.rig_batch_adjust_aux+"");
    		properties.setProperty(prefix+"rig_batch_adjust_rig",  this.rig_batch_adjust_rig+"");


    		properties.setProperty(prefix+"rig_batch_adjust_main_gt", this.rig_batch_adjust_main_gt+"");
    		properties.setProperty(prefix+"rig_batch_adjust_aux_gt",  this.rig_batch_adjust_aux_gt+"");
    		properties.setProperty(prefix+"rig_batch_adjust_rig_gt",  this.rig_batch_adjust_rig_gt+"");

    		properties.setProperty(prefix+"clt_batch_dsi1",        this.clt_batch_dsi1+"");
    		properties.setProperty(prefix+"clt_batch_apply_man",   this.clt_batch_apply_man+"");
    		properties.setProperty(prefix+"clt_batch_extrinsic",   this.clt_batch_extrinsic+"");
    		properties.setProperty(prefix+"clt_batch_poly",        this.clt_batch_poly+"");
    		properties.setProperty(prefix+"clt_batch_4img",        this.clt_batch_4img+"");
    		properties.setProperty(prefix+"clt_batch_4img_aux",    this.clt_batch_4img_aux+"");
    		properties.setProperty(prefix+"clt_batch_explore",     this.clt_batch_explore+"");
    		properties.setProperty(prefix+"clt_batch_surf",        this.clt_batch_surf+"");
    		properties.setProperty(prefix+"clt_batch_assign",      this.clt_batch_assign+"");
    		properties.setProperty(prefix+"clt_batch_gen3d",       this.clt_batch_gen3d+"");
    		properties.setProperty(prefix+"clt_batch_genMl",       this.clt_batch_genMl+"");

    		properties.setProperty(prefix+"clt_batch_dbg1",        this.clt_batch_dbg1+"");

    		properties.setProperty(prefix+"clt_batch_dsi",             this.clt_batch_dsi+"");
    		properties.setProperty(prefix+"clt_batch_dsi_aux",         this.clt_batch_dsi_aux+"");
    		properties.setProperty(prefix+"clt_batch_dsi_cm_strength", this.clt_batch_dsi_cm_strength+"");
    		properties.setProperty(prefix+"clt_batch_dsi_aux_full",    this.clt_batch_dsi_aux_full+"");
    		properties.setProperty(prefix+"clt_batch_save_extrinsics", this.clt_batch_save_extrinsics+"");
    		properties.setProperty(prefix+"clt_batch_save_all",        this.clt_batch_save_all+"");

    		properties.setProperty(prefix+"clt_batch_skip_scenes",          this.clt_batch_skip_scenes+"");

    		properties.setProperty(prefix+"clt_batch_pose_pairs_main",      this.clt_batch_pose_pairs_main+"");
    		properties.setProperty(prefix+"clt_batch_pose_last_main",       this.clt_batch_pose_last_main+"");
    		properties.setProperty(prefix+"clt_batch_pose_scene_main",      this.clt_batch_pose_scene_main+"");
    		
    		properties.setProperty(prefix+"clt_batch_offset_main",          this.clt_batch_offset_main+"");
    		properties.setProperty(prefix+"clt_batch_step_main",            this.clt_batch_step_main+"");
    		
    		properties.setProperty(prefix+"clt_batch_ml_last_main",         this.clt_batch_ml_last_main+"");
    		properties.setProperty(prefix+"clt_batch_ml_all_main",          this.clt_batch_ml_all_main+"");
    		
    		properties.setProperty(prefix+"clt_batch_pose_pairs_aux",      this.clt_batch_pose_pairs_aux+"");
    		properties.setProperty(prefix+"clt_batch_pose_last_aux",       this.clt_batch_pose_last_aux+"");
    		properties.setProperty(prefix+"clt_batch_pose_scene_aux",      this.clt_batch_pose_scene_aux+"");
    		properties.setProperty(prefix+"clt_batch_offset_aux",          this.clt_batch_offset_aux+"");
    		properties.setProperty(prefix+"clt_batch_step_aux",            this.clt_batch_step_aux+"");
    		properties.setProperty(prefix+"clt_batch_ml_last_aux",         this.clt_batch_ml_last_aux+"");
    		properties.setProperty(prefix+"clt_batch_ml_all_aux",          this.clt_batch_ml_all_aux+"");

    		properties.setProperty(prefix+"thumb_overwrite",       this.thumb_overwrite+"");
    		properties.setProperty(prefix+"thumb_width",           this.thumb_width+"");
    		properties.setProperty(prefix+"thumb_height",          this.thumb_height+"");
    		properties.setProperty(prefix+"thumb_h_center",        this.thumb_h_center+"");
    		properties.setProperty(prefix+"thumb_v_center",        this.thumb_v_center+"");
    		properties.setProperty(prefix+"thumb_size",            this.thumb_size+"");
    		properties.setProperty(prefix+"default_rating",        this.default_rating+"");


    		if (aux_camera != null) { // always
        		updateAuxFromMain();
    			String aux_prefix = prefix + AUX_PREFIX;
        		properties.setProperty(aux_prefix+"sensorDirectory",      this.aux_camera.sensorDirectory);
        		properties.setProperty(aux_prefix+"cltKernelDirectory",   this.aux_camera.cltKernelDirectory);
        		properties.setProperty(aux_prefix+"resultsDirectory",     this.aux_camera.resultsDirectory);
        		properties.setProperty(aux_prefix+"firstSubCamera",       this.aux_camera.firstSubCamera+"");
        		properties.setProperty(aux_prefix+"firstSubCameraConfig", this.aux_camera.firstSubCameraConfig+"");
        		properties.setProperty(aux_prefix+"numSubCameras",        this.aux_camera.numSubCameras+"");
        		properties.setProperty(aux_prefix+"sourcePrefix",         this.aux_camera.sourcePrefix);
        		properties.setProperty(aux_prefix+"sourceSuffix",         this.aux_camera.sourceSuffix);
        		properties.setProperty(aux_prefix+"sensorPrefix",         this.aux_camera.sensorPrefix);
        		properties.setProperty(aux_prefix+"sensorSuffix",         this.aux_camera.sensorSuffix);
        		properties.setProperty(aux_prefix+"cltKernelPrefix",      this.aux_camera.cltKernelPrefix);
        		properties.setProperty(aux_prefix+"cltSuffix",            this.aux_camera.cltSuffix);
        		properties.setProperty(aux_prefix+"x3dSubdirSuffix",      this.aux_camera.x3dSubdirSuffix);
    		}
    	}

    	public void getProperties(String prefix,Properties properties){
  		    if (properties.getProperty(prefix+"split")!=null) this.split=Boolean.parseBoolean(properties.getProperty(prefix+"split"));
  		    if (properties.getProperty(prefix+"vignetting")!=null) this.vignetting=Boolean.parseBoolean(properties.getProperty(prefix+"vignetting"));
  		    if (properties.getProperty(prefix+"pixelDefects")!=null) this.pixelDefects=Boolean.parseBoolean(properties.getProperty(prefix+"pixelDefects"));
  		    if (properties.getProperty(prefix+"pixelDefectsThreshold")!=null) this.pixelDefectsThreshold=Double.parseDouble(properties.getProperty(prefix+"pixelDefectsThreshold"));
  		    if (properties.getProperty(prefix+"debayer")!=null) this.debayer=Boolean.parseBoolean(properties.getProperty(prefix+"debayer"));
  		    if (properties.getProperty(prefix+"showDebayerEnergy")!=null) this.showDebayerEnergy=Boolean.parseBoolean(properties.getProperty(prefix+"showDebayerEnergy"));
  		    if (properties.getProperty(prefix+"saveDebayerEnergy")!=null) this.saveDebayerEnergy=Boolean.parseBoolean(properties.getProperty(prefix+"saveDebayerEnergy"));
  		    if (properties.getProperty(prefix+"deconvolve")!=null) this.deconvolve=Boolean.parseBoolean(properties.getProperty(prefix+"deconvolve"));
  		    if (properties.getProperty(prefix+"combine")!=null) this.combine=Boolean.parseBoolean(properties.getProperty(prefix+"combine"));
  		    if (properties.getProperty(prefix+"showDenoiseMask")!=null) this.showDenoiseMask=Boolean.parseBoolean(properties.getProperty(prefix+"showDenoiseMask"));
  		    if (properties.getProperty(prefix+"saveDenoiseMask")!=null) this.saveDenoiseMask=Boolean.parseBoolean(properties.getProperty(prefix+"saveDenoiseMask"));
  		    if (properties.getProperty(prefix+"showChromaDenoiseMask")!=null) this.showChromaDenoiseMask=Boolean.parseBoolean(properties.getProperty(prefix+"showChromaDenoiseMask"));
  		    if (properties.getProperty(prefix+"saveChromaDenoiseMask")!=null) this.saveChromaDenoiseMask=Boolean.parseBoolean(properties.getProperty(prefix+"saveChromaDenoiseMask"));
  		    if (properties.getProperty(prefix+"showNoiseGains")!=null) this.showNoiseGains=Boolean.parseBoolean(properties.getProperty(prefix+"showNoiseGains"));
  		    if (properties.getProperty(prefix+"saveNoiseGains")!=null) this.saveNoiseGains=Boolean.parseBoolean(properties.getProperty(prefix+"saveNoiseGains"));
  		    if (properties.getProperty(prefix+"colorProc")!=null) this.colorProc=Boolean.parseBoolean(properties.getProperty(prefix+"colorProc"));
  		    if (properties.getProperty(prefix+"blueProc")!=null) this.blueProc=Boolean.parseBoolean(properties.getProperty(prefix+"blueProc"));
  		    if (properties.getProperty(prefix+"toRGB")!=null) this.toRGB=Boolean.parseBoolean(properties.getProperty(prefix+"toRGB"));
  		    if (properties.getProperty(prefix+"rotate")!=null) this.rotate=Boolean.parseBoolean(properties.getProperty(prefix+"rotate"));
  		    if (properties.getProperty(prefix+"crop")!=null) this.crop=Boolean.parseBoolean(properties.getProperty(prefix+"crop"));   // crop to the sensor size
  		    if (properties.getProperty(prefix+"equirectangularFormat")!=null) this.equirectangularFormat=Integer.parseInt(properties.getProperty(prefix+"equirectangularFormat"));
  		    if (properties.getProperty(prefix+"outputRangeInt")!=null) this.outputRangeInt=Double.parseDouble(properties.getProperty(prefix+"outputRangeInt"));
  		    if (properties.getProperty(prefix+"outputRangeFP")!=null) this.outputRangeFP=Double.parseDouble(properties.getProperty(prefix+"outputRangeFP"));
  		    if (properties.getProperty(prefix+"imageJTags")!=null) this.imageJTags=Boolean.parseBoolean(properties.getProperty(prefix+"imageJTags"));
  		    if (properties.getProperty(prefix+"jpeg")!=null) this.jpeg=Boolean.parseBoolean(properties.getProperty(prefix+"jpeg"));   // convert to RGB and save jpeg (if save is true)
  		    if (properties.getProperty(prefix+"png")!=null) this.png=Boolean.parseBoolean(properties.getProperty(prefix+"png"));   // convert to RGB and save jpeg (if save is true)
  		    if (properties.getProperty(prefix+"save")!=null) this.save=Boolean.parseBoolean(properties.getProperty(prefix+"save"));
  		    if (properties.getProperty(prefix+"save16")!=null) this.save16=Boolean.parseBoolean(properties.getProperty(prefix+"save16")); // save 16-bit tiff also if the end result is 8 bit
  		    if (properties.getProperty(prefix+"save32")!=null) this.save32=Boolean.parseBoolean(properties.getProperty(prefix+"save32")); // save 32-bit tiff also if the end result is 8 or 16 bit
  		    if (properties.getProperty(prefix+"show")!=null) this.show=Boolean.parseBoolean(properties.getProperty(prefix+"show"));
  		    if (properties.getProperty(prefix+"JPEG_quality")!=null) this.JPEG_quality=Integer.parseInt(properties.getProperty(prefix+"JPEG_quality"));
  		    if (properties.getProperty(prefix+"JPEG_scale")!=null) this.JPEG_scale=Double.parseDouble(properties.getProperty(prefix+"JPEG_scale"));
  		    if (properties.getProperty(prefix+"equirectangular")!=null) this.equirectangular=Boolean.parseBoolean(properties.getProperty(prefix+"equirectangular"));
  		    if (properties.getProperty(prefix+"zcorrect")!=null) this.zcorrect=Boolean.parseBoolean(properties.getProperty(prefix+"zcorrect"));
  		    if (properties.getProperty(prefix+"saveSettings")!=null) this.saveSettings=Boolean.parseBoolean(properties.getProperty(prefix+"saveSettings"));
			if (properties.getProperty(prefix+"sourceDirectory")!=      null) this.sourceDirectory=properties.getProperty(prefix+"sourceDirectory");
			if (properties.getProperty(prefix+"tile_processor_gpu")!=      null) this.tile_processor_gpu=properties.getProperty(prefix+"tile_processor_gpu");
  		    if (properties.getProperty(prefix+"firstSubCamera")!=       null) this.firstSubCamera=Integer.parseInt(properties.getProperty(prefix+"firstSubCamera"));
  		    if (properties.getProperty(prefix+"firstSubCameraConfig")!= null) this.firstSubCameraConfig=Integer.parseInt(properties.getProperty(prefix+"firstSubCameraConfig"));
  		    if (properties.getProperty(prefix+"numSubCameras")!=        null) this.numSubCameras=Integer.parseInt(properties.getProperty(prefix+"numSubCameras"));
			if (properties.getProperty(prefix+"sensorDirectory")!=      null) this.sensorDirectory=properties.getProperty(prefix+"sensorDirectory");
			if (properties.getProperty(prefix+"use_set_dirs")!=         null) this.use_set_dirs=Boolean.parseBoolean(properties.getProperty(prefix+"use_set_dirs"));
			if (properties.getProperty(prefix+"sourcePrefix")!=         null) this.sourcePrefix=properties.getProperty(prefix+"sourcePrefix");
			if (properties.getProperty(prefix+"sourceSuffix")!=         null) this.sourceSuffix=properties.getProperty(prefix+"sourceSuffix");
			if (properties.getProperty(prefix+"sensorPrefix")!=         null) this.sensorPrefix=properties.getProperty(prefix+"sensorPrefix");
			if (properties.getProperty(prefix+"sensorSuffix")!=         null) this.sensorSuffix=properties.getProperty(prefix+"sensorSuffix");
			if (properties.getProperty(prefix+"sharpKernelDirectory")!= null) this.sharpKernelDirectory=properties.getProperty(prefix+"sharpKernelDirectory");
			if (properties.getProperty(prefix+"sharpKernelPrefix")!=    null) this.sharpKernelPrefix=properties.getProperty(prefix+"sharpKernelPrefix");
			if (properties.getProperty(prefix+"sharpKernelSuffix")!=    null) this.sharpKernelSuffix=properties.getProperty(prefix+"sharpKernelSuffix");
			if (properties.getProperty(prefix+"smoothKernelDirectory")!=null) this.smoothKernelDirectory=properties.getProperty(prefix+"smoothKernelDirectory");
			if (properties.getProperty(prefix+"smoothKernelPrefix")!=   null) this.smoothKernelPrefix=properties.getProperty(prefix+"smoothKernelPrefix");
			if (properties.getProperty(prefix+"smoothKernelSuffix")!=   null) this.smoothKernelSuffix=properties.getProperty(prefix+"smoothKernelSuffix");

			if (properties.getProperty(prefix+"dctKernelDirectory")!=   null) this.dctKernelDirectory=properties.getProperty(prefix+"dctKernelDirectory");
			if (properties.getProperty(prefix+"dctKernelPrefix")!=      null) this.dctKernelPrefix=properties.getProperty(prefix+"dctKernelPrefix");
			if (properties.getProperty(prefix+"dctSymSuffix")!=         null) this.dctSymSuffix=properties.getProperty(prefix+"dctSymSuffix");
			if (properties.getProperty(prefix+"dctAsymSuffix")!=        null) this.dctAsymSuffix=properties.getProperty(prefix+"dctAsymSuffix");

			if (properties.getProperty(prefix+"equirectangularDirectory")!=null) this.equirectangularDirectory=properties.getProperty(prefix+"equirectangularDirectory");
			if (properties.getProperty(prefix+"equirectangularPrefix")!=null) this.equirectangularPrefix=properties.getProperty(prefix+"equirectangularPrefix");
			if (properties.getProperty(prefix+"equirectangularSuffix")!=null) this.equirectangularSuffix=properties.getProperty(prefix+"equirectangularSuffix");

			if (properties.getProperty(prefix+"equirectangularCut")!=null)
				this.equirectangularCut=Boolean.parseBoolean(properties.getProperty(prefix+"equirectangularCut"));
//			if (properties.getProperty(prefix+"equirectangularSuffixA")!=null) this.equirectangularSuffixA=properties.getProperty(prefix+"equirectangularSuffixA");

			if (properties.getProperty(prefix+"planeMapPrefix")!=null) this.planeMapPrefix=properties.getProperty(prefix+"planeMapPrefix");
			if (properties.getProperty(prefix+"planeMapSuffix")!=null) this.planeMapSuffix=properties.getProperty(prefix+"planeMapSuffix");
			if (properties.getProperty(prefix+"usePlaneProjection")!=null)
				this.usePlaneProjection=Boolean.parseBoolean(properties.getProperty(prefix+"usePlaneProjection"));
			if (properties.getProperty(prefix+"planeAsJPEG")!=null)
				this.planeAsJPEG=Boolean.parseBoolean(properties.getProperty(prefix+"planeAsJPEG"));
			if (properties.getProperty(prefix+"resultsDirectory")!=     null) this.resultsDirectory=properties.getProperty(prefix+"resultsDirectory");
			if (properties.getProperty(prefix+"removeUnusedSensorData")!= null)
				this.removeUnusedSensorData=Boolean.parseBoolean(properties.getProperty(prefix+"removeUnusedSensorData"));
			if (properties.getProperty(prefix+"sourcePaths")!=   null){
				int numFiles=Integer.parseInt(properties.getProperty(prefix+"sourcePaths"));
				this.sourcePaths=new String[numFiles];
				for (int i=0;i<this.sourcePaths.length;i++){
					this.sourcePaths[i]=properties.getProperty(prefix+"sourcePath"+i);
        		}
			}
  		    if (properties.getProperty(prefix+"exposureCorrectionMode")!=null) this.exposureCorrectionMode=Integer.parseInt(properties.getProperty(prefix+"exposureCorrectionMode"));
  		    if (properties.getProperty(prefix+"referenceExposure")     !=null) this.referenceExposure=   Double.parseDouble(properties.getProperty(prefix+"referenceExposure"));
  		    if (properties.getProperty(prefix+"relativeExposure")      !=null) this.relativeExposure=    Double.parseDouble(properties.getProperty(prefix+"relativeExposure"));
  		    if (properties.getProperty(prefix+"swapSubchannels01")!=null) this.swapSubchannels01=Boolean.parseBoolean(properties.getProperty(prefix+"swapSubchannels01"));

			if (properties.getProperty(prefix+"cltKernelDirectory")!=   null) this.cltKernelDirectory=properties.getProperty(prefix+"cltKernelDirectory");
			if (properties.getProperty(prefix+"cltKernelPrefix")!=      null) this.cltKernelPrefix=properties.getProperty(prefix+"cltKernelPrefix");
			if (properties.getProperty(prefix+"cltSuffix")!=            null) this.cltSuffix=properties.getProperty(prefix+"cltSuffix");

			if (properties.getProperty(prefix+"x3dDirectory")!=         null) this.x3dDirectory=properties.getProperty(prefix+"x3dDirectory");

			if (properties.getProperty(prefix+"use_x3d_subdirs")!= null) this.use_x3d_subdirs=Boolean.parseBoolean(properties.getProperty(prefix+"use_x3d_subdirs"));

			if (properties.getProperty(prefix+"x3dSubdirPrefix")!=      null) this.x3dSubdirPrefix=properties.getProperty(prefix+"x3dSubdirPrefix");
			if (properties.getProperty(prefix+"x3dSubdirSuffix")!=      null) this.x3dSubdirSuffix=properties.getProperty(prefix+"x3dSubdirSuffix");

			if (properties.getProperty(prefix+"x3dModelVersion")!=      null) this.x3dModelVersion=properties.getProperty(prefix+"x3dModelVersion");
			if (properties.getProperty(prefix+"jp4SubDir")!=            null) this.jp4SubDir=properties.getProperty(prefix+"jp4SubDir");

			if (properties.getProperty(prefix+"mlDirectory")!=          null) this.mlDirectory=properties.getProperty(prefix+"mlDirectory");

			if (properties.getProperty(prefix+"process_main_sources")!= null) this.process_main_sources=Boolean.parseBoolean(properties.getProperty(prefix+"process_main_sources"));
			if (properties.getProperty(prefix+"process_aux_sources")!= null)  this.process_aux_sources=Boolean.parseBoolean(properties.getProperty(prefix+"process_aux_sources"));

  		    if (properties.getProperty(prefix+"kml_sensors")!=null) this.kml_sensors=Integer.parseInt(properties.getProperty(prefix+"kml_sensors"));
  		    
  		    if (properties.getProperty(prefix+"rig_batch_adjust_main")!=null) this.rig_batch_adjust_main=Integer.parseInt(properties.getProperty(prefix+"rig_batch_adjust_main"));
  		    if (properties.getProperty(prefix+"rig_batch_adjust_aux")!=null)  this.rig_batch_adjust_aux=Integer.parseInt(properties.getProperty(prefix+"rig_batch_adjust_aux"));
  		    if (properties.getProperty(prefix+"rig_batch_adjust_rig")!=null)  this.rig_batch_adjust_rig=Integer.parseInt(properties.getProperty(prefix+"rig_batch_adjust_rig"));
  		    if (properties.getProperty(prefix+"rig_batch_adjust_main_gt")!=null)  this.rig_batch_adjust_main_gt=Integer.parseInt(properties.getProperty(prefix+"rig_batch_adjust_main_gt"));
  		    if (properties.getProperty(prefix+"rig_batch_adjust_aux_gt")!=null)  this.rig_batch_adjust_aux_gt=Integer.parseInt(properties.getProperty(prefix+"rig_batch_adjust_aux_gt"));
  		    if (properties.getProperty(prefix+"rig_batch_adjust_rig_gt")!=null)  this.rig_batch_adjust_rig_gt=Integer.parseInt(properties.getProperty(prefix+"rig_batch_adjust_rig_gt"));

			if (properties.getProperty(prefix+"clt_batch_dsi1")!= null)      this.clt_batch_dsi1=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_dsi1"));
			if (properties.getProperty(prefix+"clt_batch_apply_man")!= null) this.clt_batch_apply_man=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_apply_man"));
			if (properties.getProperty(prefix+"clt_batch_extrinsic")!= null) this.clt_batch_extrinsic=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_extrinsic"));
			if (properties.getProperty(prefix+"clt_batch_poly")!= null)      this.clt_batch_poly=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_poly"));
			if (properties.getProperty(prefix+"clt_batch_4img")!= null)      this.clt_batch_4img=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_4img"));
			if (properties.getProperty(prefix+"clt_batch_4img_aux")!= null)  this.clt_batch_4img_aux=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_4img_aux"));
			if (properties.getProperty(prefix+"clt_batch_explore")!= null)   this.clt_batch_explore=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_explore"));
			if (properties.getProperty(prefix+"clt_batch_surf")!= null)      this.clt_batch_surf=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_surf"));
			if (properties.getProperty(prefix+"clt_batch_assign")!= null)    this.clt_batch_assign=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_assign"));
			if (properties.getProperty(prefix+"clt_batch_gen3d")!= null)     this.clt_batch_gen3d=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_gen3d"));

			if (properties.getProperty(prefix+"clt_batch_genMl")!= null)     this.clt_batch_genMl=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_genMl"));
			if (properties.getProperty(prefix+"clt_batch_dbg1")!= null)      this.clt_batch_dbg1=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_dbg1"));

			if (properties.getProperty(prefix+"clt_batch_dsi")!= null)             this.clt_batch_dsi=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_dsi"));
			if (properties.getProperty(prefix+"clt_batch_dsi_aux")!= null)         this.clt_batch_dsi_aux=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_dsi_aux"));
			if (properties.getProperty(prefix+"clt_batch_dsi_cm_strength")!= null) this.clt_batch_dsi_cm_strength=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_dsi_cm_strength"));
		
			if (properties.getProperty(prefix+"clt_batch_dsi_aux_full")!= null)    this.clt_batch_dsi_aux_full=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_dsi_aux_full"));
			if (properties.getProperty(prefix+"clt_batch_save_extrinsics")!= null) this.clt_batch_save_extrinsics=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_save_extrinsics"));
			if (properties.getProperty(prefix+"clt_batch_save_all")!= null)        this.clt_batch_save_all=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_save_all"));

			if (properties.getProperty(prefix+"clt_batch_skip_scenes")!= null)     this.clt_batch_skip_scenes=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_skip_scenes"));

			if (properties.getProperty(prefix+"clt_batch_pose_pairs_main")!= null) this.clt_batch_pose_pairs_main=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_pose_pairs_main"));
			if (properties.getProperty(prefix+"clt_batch_pose_last_main")!= null)  this.clt_batch_pose_last_main=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_pose_last_main"));
			if (properties.getProperty(prefix+"clt_batch_pose_scene_main")!= null) this.clt_batch_pose_scene_main=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_pose_scene_main"));
  		    if (properties.getProperty(prefix+"clt_batch_offset_main")!=null)      this.clt_batch_offset_main=Integer.parseInt(properties.getProperty(prefix+"clt_batch_offset_main"));
  		    if (properties.getProperty(prefix+"clt_batch_step_main")!=null)        this.clt_batch_step_main=Integer.parseInt(properties.getProperty(prefix+"clt_batch_step_main"));
			if (properties.getProperty(prefix+"clt_batch_ml_last_main")!= null)    this.clt_batch_ml_last_main=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_ml_last_main"));
			if (properties.getProperty(prefix+"clt_batch_ml_all_main")!= null)     this.clt_batch_ml_all_main=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_ml_all_main"));
			
			if (properties.getProperty(prefix+"clt_batch_pose_pairs_aux")!= null)  this.clt_batch_pose_pairs_aux=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_pose_pairs_aux"));
			if (properties.getProperty(prefix+"clt_batch_pose_last_aux")!= null)   this.clt_batch_pose_last_aux=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_pose_last_aux"));
			if (properties.getProperty(prefix+"clt_batch_pose_scene_aux")!= null)  this.clt_batch_pose_scene_aux=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_pose_scene_aux"));
  		    if (properties.getProperty(prefix+"clt_batch_offset_aux")!=null)       this.clt_batch_offset_aux=Integer.parseInt(properties.getProperty(prefix+"clt_batch_offset_aux"));
  		    if (properties.getProperty(prefix+"clt_batch_step_aux")!=null)         this.clt_batch_step_aux=Integer.parseInt(properties.getProperty(prefix+"clt_batch_step_aux"));
			if (properties.getProperty(prefix+"clt_batch_ml_last_aux")!= null)     this.clt_batch_ml_last_aux=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_ml_last_aux"));
			if (properties.getProperty(prefix+"clt_batch_ml_all_aux")!= null)      this.clt_batch_ml_all_aux=Boolean.parseBoolean(properties.getProperty(prefix+"clt_batch_ml_all_aux"));
			
			if (properties.getProperty(prefix+"thumb_overwrite")!= null)     this.thumb_overwrite=Boolean.parseBoolean(properties.getProperty(prefix+"thumb_overwrite"));
			if (properties.getProperty(prefix+"thumb_width")!=null)          this.thumb_width=Integer.parseInt(properties.getProperty(prefix+"thumb_width"));
  		    if (properties.getProperty(prefix+"thumb_height")!=null)         this.thumb_height=Integer.parseInt(properties.getProperty(prefix+"thumb_height"));
  		    if (properties.getProperty(prefix+"thumb_h_center")!=null)       this.thumb_h_center=   Double.parseDouble(properties.getProperty(prefix+"thumb_h_center"));
  		    if (properties.getProperty(prefix+"thumb_v_center")!=null)       this.thumb_v_center=   Double.parseDouble(properties.getProperty(prefix+"thumb_v_center"));
  		    if (properties.getProperty(prefix+"thumb_size")   !=null)        this.thumb_size=   Double.parseDouble(properties.getProperty(prefix+"thumb_size"));
  		    if (properties.getProperty(prefix+"default_rating") !=null)      this.default_rating=   Integer.parseInt(properties.getProperty(prefix+"default_rating"));

    		// copy common parameters to the auxiliary camera ones
    		updateAuxFromMain();
			String aux_prefix = prefix + AUX_PREFIX;
			if (properties.getProperty(aux_prefix+"sensorDirectory")!=      null) this.aux_camera.sensorDirectory=properties.getProperty(aux_prefix+"sensorDirectory");
			if (properties.getProperty(aux_prefix+"cltKernelDirectory")!=   null) this.aux_camera.cltKernelDirectory=properties.getProperty(aux_prefix+"cltKernelDirectory");
			if (properties.getProperty(aux_prefix+"resultsDirectory")!=     null) this.aux_camera.resultsDirectory=properties.getProperty(aux_prefix+"resultsDirectory");
  		    if (properties.getProperty(aux_prefix+"firstSubCamera")!=       null) this.aux_camera.firstSubCamera=Integer.parseInt(properties.getProperty(aux_prefix+"firstSubCamera"));
  		    if (properties.getProperty(aux_prefix+"firstSubCameraConfig")!= null) this.aux_camera.firstSubCameraConfig=Integer.parseInt(properties.getProperty(aux_prefix+"firstSubCameraConfig"));
  		    if (properties.getProperty(aux_prefix+"numSubCameras")!=        null) this.aux_camera.numSubCameras=Integer.parseInt(properties.getProperty(aux_prefix+"numSubCameras"));
  		    if (properties.getProperty(aux_prefix+"sourcePrefix")!=         null) this.aux_camera.sourcePrefix=properties.getProperty(aux_prefix+"sourcePrefix");
			if (properties.getProperty(aux_prefix+"sourceSuffix")!=         null) this.aux_camera.sourceSuffix=properties.getProperty(aux_prefix+"sourceSuffix");
			if (properties.getProperty(aux_prefix+"sensorPrefix")!=         null) this.aux_camera.sensorPrefix=properties.getProperty(aux_prefix+"sensorPrefix");
			if (properties.getProperty(aux_prefix+"sensorSuffix")!=         null) this.aux_camera.sensorSuffix=properties.getProperty(aux_prefix+"sensorSuffix");
			if (properties.getProperty(aux_prefix+"cltKernelPrefix")!=      null) this.aux_camera.cltKernelPrefix=properties.getProperty(aux_prefix+"cltKernelPrefix");
			if (properties.getProperty(aux_prefix+"cltSuffix")!=            null) this.aux_camera.cltSuffix=properties.getProperty(aux_prefix+"cltSuffix");
			if (properties.getProperty(aux_prefix+"x3dSubdirSuffix")!=      null) this.aux_camera.x3dSubdirSuffix=properties.getProperty(aux_prefix+"x3dSubdirSuffix");
		}

    	public boolean showJDialog(String title) {
//    		GenericDialog gd = new GenericDialog(title);
    		GenericJTabbedDialog gd = new GenericJTabbedDialog(title ,1000, 900);
            gd.addTab("Eyesis parameters","Eyesis camera parameters, most not all applicable to quad cameras");
    		gd.addCheckbox ("Splt into Bayer stack (if false will exit)",       this.split);
    		gd.addCheckbox ("Apply vignetting/color correction to source files",this.vignetting);
    		gd.addCheckbox ("Replace hot/warm/cold pixels with average of neighbors",this.pixelDefects);
    		gd.addNumericField("Pixel difference thershold to consider it \"bad\" on 255.0 scale (0 - use all)", this.pixelDefectsThreshold, 2,6,"8.0");
			String [] choices={"none","absolute","relative"};
			if (this.exposureCorrectionMode<0) this.exposureCorrectionMode=0;
			else if (this.exposureCorrectionMode>=choices.length) this.exposureCorrectionMode=choices.length-1;
			gd.addChoice      ("Exposure correction",choices, choices[this.exposureCorrectionMode]);
    		gd.addNumericField("Reference exposure (effective only in \"absolute\" mode)", 1000.0*this.referenceExposure, 2,6,"ms");
    		gd.addNumericField("Exposure scale (effective only in \"relative\" mode) 0 - darken, 1 - lighten", this.relativeExposure, 3,5,"");
    		gd.addCheckbox ("De-mosaic (if false will exit)",                   this.debayer);
    		gd.addCheckbox ("Show de-mosaic middle-frequency 'energy",          this.showDebayerEnergy);
    		gd.addCheckbox ("Save de-mosaic middle-frequency 'energy",          this.saveDebayerEnergy);
    		gd.addCheckbox ("Sharpen (convolve with calibration kernels)",      this.deconvolve);
    		gd.addCheckbox ("Denoise (convolve with Gaussian in smooth areas)", this.combine);
    		gd.addCheckbox ("Show denoise mask (white - use hi-res, black - low-res)", this.showDenoiseMask);
    		gd.addCheckbox ("Save denoise mask (white - use hi-res, black - low-res)", this.saveDenoiseMask);
    		gd.addCheckbox ("Show kernel noise gains",                          this.showNoiseGains);
    		gd.addCheckbox ("Save kernel noise gains",                          this.saveNoiseGains);
    		gd.addCheckbox ("Convert colors",                                   this.colorProc);
    		gd.addCheckbox ("Fix blue leak",                                    this.blueProc);
    		gd.addCheckbox ("Show chroma denoise mask (white - use hi-res, black - low-res)", this.showChromaDenoiseMask);
    		gd.addCheckbox ("Save chroma denoise mask (white - use hi-res, black - low-res)", this.saveChromaDenoiseMask);
    		gd.addCheckbox ("Rotate result image",                              this.rotate);
    		gd.addCheckbox ("Crop result image to the original size",           this.crop);
			String [] equirectangularFormatChoices={"RGBA 8-bit","RGBA 16-bit","RGBA 32-bit integer","RGBA 32-bit float","ImageJ stack"};
			int [] equirectangularFormats={0,1,2,3,4};
			int equirectangularFormatIndex=0;
			for ( int i=0;i<equirectangularFormats.length;i++) if (equirectangularFormats[i]==this.equirectangularFormat){
				equirectangularFormatIndex=i;
				break;
			}
			gd.addChoice   ("Equirectangular output format",equirectangularFormatChoices, equirectangularFormatChoices[equirectangularFormatIndex]);
    		gd.addNumericField("Map 1.0 intensity to this fraction of the full range 8/16/32-bit integer mode output", 100*this.outputRangeInt, 2,6,"%");
    		gd.addNumericField("Map 1.0 intensity to this value in 32-bit floating point output mode", this.outputRangeFP, 2,6,"");
    		gd.addCheckbox ("Encode ImageJ specific Info metadata to the output file TIFF header", this.imageJTags);

			gd.addCheckbox ("Convert to RGB48",                                 this.toRGB);
    		gd.addCheckbox ("Convert to 8 bit RGB (and save JPEG if save is enabled)", this.jpeg);
    		gd.addCheckbox ("Use PNG instead of TIFF for 32 bit (8 per color) RGBA", this.png);
    		gd.addCheckbox ("Save the result to file system",                   this.save);
    		gd.addCheckbox ("Save 16-bit tiff if the result is 8 bit",          this.save16);
    		gd.addCheckbox    ("Save 32-bit tiff if the result is 8 or 16 bit",    this.save32);
    		gd.addCheckbox    ("Show the result image",                            this.show);
    		gd.addNumericField("JPEG quality (%)",                                 this.JPEG_quality,0);
    		gd.addNumericField("JPEG scale   (%)",                            100* this.JPEG_scale,0);
    		gd.addCheckbox    ("Warp results to equirectangular",                  this.equirectangular);
    		gd.addCheckbox    ("Calculate distances in overlapping areas",         this.zcorrect);
    		gd.addCheckbox    ("Save current settings with results",               this.saveSettings);

            gd.addTab("Directories","Direcories paths");
    		gd.addStringField ("Source files directory",                           this.sourceDirectory, 60);
    		gd.addCheckbox    ("Select source directory",                          false);
    		gd.addCheckbox    ("Use individual subdirectory for image set",        this.use_set_dirs);

    		gd.addStringField ("Sensor calibration directory",                     this.sensorDirectory, 60);
    		gd.addCheckbox    ("Select sensor calibration directory",              false);

    		gd.addStringField ("Aberration kernels (sharp) directory",             this.sharpKernelDirectory, 60);
    		gd.addCheckbox    ("Select aberration kernels (sharp) directory",      false);
    		gd.addStringField ("Aberration kernels (smooth) directory",            this.smoothKernelDirectory, 60);
    		gd.addCheckbox    ("Select aberration kernels (smooth) directory",     false);

    		gd.addStringField ("Aberration kernels for DCT directory",             this.dctKernelDirectory, 60);
    		gd.addCheckbox    ("Select aberration kernels for DCT directory",      false);

    		gd.addStringField ("Aberration kernels for CLT directory",             this.cltKernelDirectory, 60);
    		gd.addCheckbox    ("Select aberration kernels for CLT directory",      false);

    		gd.addStringField ("x3d model version",                                this.x3dModelVersion, 20);    // 10a
    		gd.addStringField ("JP4 source image copy model subdirectory",         this.jp4SubDir, 20);    // 10b
    		gd.addStringField ("x3d output directory",                             this.x3dDirectory, 60);
    		gd.addCheckbox    ("Select x3d output directory",                      false);
    		gd.addCheckbox    ("Use individual subdirectory for each 3d model (timestamp as name)", this.use_x3d_subdirs);

    		gd.addStringField ("x3d subdirectory prefix",                          this.x3dSubdirPrefix, 10,
    				"When using timestamp as a subdirectory, add this prefix");
    		gd.addStringField ("x3d subdirectory suffix",                          this.x3dSubdirSuffix, 10,
    				"When using timestamp as a subdirectory, add this suffix");
    		gd.addStringField ("ML output directory",                              this.mlDirectory, 60,
    				"Non-empty directory with no \"/\" separator makes it a subdirectory of the model version directory");
    		gd.addCheckbox    ("Select ML output directory", false,"Erase text field or use \"/\" in it to enable absolute directory path selection");

    		gd.addStringField("Equirectangular maps directory (may be empty)",     this.equirectangularDirectory, 60);
    		gd.addCheckbox("Select equirectangular maps directory",                false);
    		gd.addStringField("Results directory",                                 this.resultsDirectory, 60);
    		gd.addCheckbox("Select results directory",                             false);

            gd.addTab("Prefix/suffix","Prefixes and suffixes for various file types");
    		gd.addStringField("Source files prefix",                               this.sourcePrefix, 60);
    		gd.addStringField("Source files suffix",                               this.sourceSuffix, 60);
    		gd.addNumericField("First subcamera (in the source filenames)",        this.firstSubCamera, 0);
    		gd.addNumericField("First subcamera (in config (clt, sensor) directories)",        this.firstSubCameraConfig, 0);
    		gd.addNumericField("Number of subcameras in this camera",              this.numSubCameras, 0);

    		gd.addStringField("Sensor files prefix",                               this.sensorPrefix, 40);
    		gd.addStringField("Sensor files suffix",                               this.sensorSuffix, 40);
    		gd.addStringField("Kernel files (sharp) prefix",                       this.sharpKernelPrefix, 40);
    		gd.addStringField("Kernel files (sharp) suffix",                       this.sharpKernelSuffix, 40);
    		gd.addStringField("Kernel files (smooth) prefix",                      this.smoothKernelPrefix, 40);
    		gd.addStringField("Kernel files (smooth) suffix",                      this.smoothKernelSuffix, 40);

    		gd.addStringField("DCT kernel files  prefix",                          this.dctKernelPrefix, 40);
    		gd.addStringField("DCT symmetical kernel files",                       this.dctSymSuffix, 40);
    		gd.addStringField("DCT asymmetrical kernel files suffix",              this.dctAsymSuffix, 40);
    		gd.addStringField("CLT kernel files prefix",                           this.cltKernelPrefix, 40);
    		gd.addStringField("CLT kernel files suffix",                           this.cltSuffix, 40);

    		gd.addStringField("Equirectangular maps prefix",     this.equirectangularPrefix, 40);
    		gd.addStringField("Equirectangular maps suffix",     this.equirectangularSuffix, 40);
    		gd.addCheckbox("Cut rolling-over equirectangular images in two", this.equirectangularCut);

    		gd.addStringField("Plane projection map prefix",     this.planeMapPrefix, 40);
    		gd.addStringField("Plane projection map suffix",     this.planeMapSuffix, 40);
    		gd.addCheckbox("Use projection to a common plane instead of the  equirectangular", this.usePlaneProjection);
    		gd.addCheckbox("Save de-warped images as JPEG instead of TIFF",  this.planeAsJPEG);


//    		gd.addStringField("Suffix for the second part of rolled-over equirectangular images",  this.equirectangularSuffixA, 40);

    		gd.addCheckbox   ("Remove unused sensor data",       this.removeUnusedSensorData);
    		gd.addCheckbox   ("Swap top and equator images",     this.swapSubchannels01);

    		gd.addTab("Thumbnails","Thumbnail image generation");
    		gd.addCheckbox("Overwrite existing thumbnail images",this.thumb_overwrite);
    		gd.addNumericField("Thumbnail image width",          this.thumb_width, 0,4,"pix",
    				"");
    		gd.addNumericField("Thumbnail image height",         this.thumb_height, 0,4,"pix",
    				"");
    		gd.addNumericField("Thumbnail image center horizontally", this.thumb_h_center, 2,6,"",
    				"0.0 - touch left margin, 1.0 - touch right margin");
    		gd.addNumericField("Thumbnail image center vertically", this.thumb_v_center, 2,6,"",
    				"0.0 - touch top margin, 1.0 - touch bottom margin");
    		gd.addNumericField("Thumbnail relative image size",  this.thumb_size, 2,6,"",
    				"1.0 - maximal to fit frame");
    		gd.addNumericField("Default scene rating",           this.default_rating, 0,2,"",
    				"Determins scene filtering");


//    		WindowTools.addScrollBars(gd);
    		gd.showDialog();
    		if (gd.wasCanceled()) return false;
    		this.split=                  gd.getNextBoolean();
    		this.vignetting=             gd.getNextBoolean();
    		this.pixelDefects=           gd.getNextBoolean();
    		this.pixelDefectsThreshold=  gd.getNextNumber();
			this.exposureCorrectionMode= gd.getNextChoiceIndex();
    		this.referenceExposure=0.001*gd.getNextNumber();
    		this.relativeExposure=       gd.getNextNumber();
    		this.debayer=           gd.getNextBoolean();
    		this.showDebayerEnergy= gd.getNextBoolean();
    		this.saveDebayerEnergy= gd.getNextBoolean();
    		this.deconvolve=        gd.getNextBoolean();
    		this.combine=           gd.getNextBoolean();
    		this.showDenoiseMask=   gd.getNextBoolean();
    		this.saveDenoiseMask=   gd.getNextBoolean();
    		this.showNoiseGains=    gd.getNextBoolean();
    		this.saveNoiseGains=    gd.getNextBoolean();
    		this.colorProc=         gd.getNextBoolean();
    		this.blueProc=         gd.getNextBoolean();
    		this.showChromaDenoiseMask=   gd.getNextBoolean();
    		this.saveChromaDenoiseMask=   gd.getNextBoolean();
    		this.rotate=            gd.getNextBoolean();
    		this.crop=              gd.getNextBoolean();
    		this.equirectangularFormat= equirectangularFormats[gd.getNextChoiceIndex()];
    		this.outputRangeInt=0.01*gd.getNextNumber();
    		this.outputRangeFP=     gd.getNextNumber();
    		this.imageJTags=        gd.getNextBoolean();
    		this.toRGB=             gd.getNextBoolean();
    		this.jpeg=              gd.getNextBoolean();
    		this.png=               gd.getNextBoolean();
    		this.save=              gd.getNextBoolean();
    		this.save16=            gd.getNextBoolean();
    		this.save32=            gd.getNextBoolean();
    		this.show=              gd.getNextBoolean();
    		this.JPEG_quality=(int) gd.getNextNumber();
    		this.JPEG_scale=   0.01*gd.getNextNumber();
    		this.equirectangular=   gd.getNextBoolean();
    		this.zcorrect=          gd.getNextBoolean();
    		this.saveSettings=      gd.getNextBoolean();

    		this.sourceDirectory=        gd.getNextString(); if (gd.getNextBoolean()) selectSourceDirectory(false, false);
    		this.use_set_dirs =          gd.getNextBoolean();
    		this.sensorDirectory=        gd.getNextString(); if (gd.getNextBoolean()) selectSensorDirectory(false, false);
    		this.sharpKernelDirectory=   gd.getNextString(); if (gd.getNextBoolean()) selectSharpKernelDirectory(false, false);
    		this.smoothKernelDirectory=  gd.getNextString(); if (gd.getNextBoolean()) selectSmoothKernelDirectory(false, true);
    		this.dctKernelDirectory=     gd.getNextString(); if (gd.getNextBoolean()) selectDCTKernelDirectory(false, true);
    		this.cltKernelDirectory=     gd.getNextString(); if (gd.getNextBoolean()) selectCLTKernelDirectory(false, true);
    		this.x3dModelVersion=        gd.getNextString(); // 10a
    		this.jp4SubDir=              gd.getNextString(); // 10b
    		this.x3dDirectory=           gd.getNextString(); if (gd.getNextBoolean()) selectX3dDirectory(false, true);
    		this.use_x3d_subdirs=        gd.getNextBoolean();

    		this.x3dSubdirPrefix=        gd.getNextString();
    		this.x3dSubdirSuffix=        gd.getNextString();

    		this.mlDirectory=            gd.getNextString(); if (gd.getNextBoolean()) selectMlDirectory(null, false, true);

    		this.equirectangularDirectory=  gd.getNextString(); if (gd.getNextBoolean()) selectEquirectangularDirectory(false, false);
    		this.resultsDirectory=       gd.getNextString(); if (gd.getNextBoolean()) selectResultsDirectory(false, true);
    		this.sourcePrefix=           gd.getNextString();
    		this.sourceSuffix=           gd.getNextString();
    		this.firstSubCamera=   (int) gd.getNextNumber();
    		this.firstSubCameraConfig=(int) gd.getNextNumber();
    		this.numSubCameras=    (int) gd.getNextNumber();
    		this.sensorPrefix=           gd.getNextString();
    		this.sensorSuffix=           gd.getNextString();
    		this.sharpKernelPrefix=      gd.getNextString();
    		this.sharpKernelSuffix=      gd.getNextString();
    		this.smoothKernelPrefix=     gd.getNextString();
    		this.smoothKernelSuffix=     gd.getNextString();
    		this.dctKernelPrefix=        gd.getNextString();
    		this.dctSymSuffix=           gd.getNextString();
    		this.dctAsymSuffix=          gd.getNextString();
    		this.cltKernelPrefix=        gd.getNextString();
    		this.cltSuffix=              gd.getNextString();
    		this.equirectangularPrefix=  gd.getNextString();
    		this.equirectangularSuffix=  gd.getNextString();
    		this.equirectangularCut=     gd.getNextBoolean();
    		this.planeMapPrefix=         gd.getNextString();
    		this.planeMapSuffix=         gd.getNextString();
    		this.usePlaneProjection=     gd.getNextBoolean();
    		this.planeAsJPEG=            gd.getNextBoolean();
//    		this.equirectangularSuffixA= gd.getNextString();
    		this.removeUnusedSensorData= gd.getNextBoolean();
    		this.swapSubchannels01=      gd.getNextBoolean();

    		this.thumb_overwrite =       gd.getNextBoolean();
    		this.thumb_width=      (int) gd.getNextNumber();
    		this.thumb_height=     (int) gd.getNextNumber();
    		this.thumb_h_center=         gd.getNextNumber();
    		this.thumb_v_center=         gd.getNextNumber();
    		this.thumb_size=             gd.getNextNumber();
    		this.default_rating=   (int) gd.getNextNumber();
    		return true;
    	}


    	public boolean showCLTBatchDialog(String title,
    			CLTParameters clt_parameters) {
    		GenericJTabbedDialog gd = new GenericJTabbedDialog(title,1000,1000);
    		updateAuxFromMain();


    		gd.addTab         ("File paths", "Select files and directories paths (common to main and optional auxiliary)");
			gd.addMessage     ("============ Common to the main and optional auxiliary camera============");
    		gd.addStringField ("GPU tile_processor_gpu project absolute path",     this.tile_processor_gpu, 60,
    				"Keep empty to use default GPU kernels");
    		gd.addCheckbox    ("Select GPU directory",                             false);

    		gd.addCheckbox    ("Save current settings with results",               this.saveSettings);           // 1
    		gd.addStringField ("Source files directory",                           this.sourceDirectory, 60);    // 2
    		gd.addCheckbox    ("Select source directory",                          false);                       // 3
    		gd.addCheckbox    ("Use individual subdirectory for each image set (timestamp as name)", this.use_set_dirs); //10

    		gd.addStringField ("x3d model version",                                this.x3dModelVersion, 60);    // 10a
    		gd.addStringField ("jp4 source copy subdirectory",                     this.jp4SubDir, 60);          // 10b
    		gd.addStringField ("x3d output directory",                             this.x3dDirectory, 60);       // 8
    		gd.addCheckbox    ("Select x3d output (top model) directory",          false);                       // 9

    		gd.addCheckbox    ("Use individual subdirectory for each 3d model (timestamp as name)", this.use_x3d_subdirs); //10

//    		gd.addStringField ("Source files prefix",                              this.sourcePrefix, 60);       // 13
//    		gd.addStringField ("Source files suffix",                              this.sourceSuffix, 60);       // 14

    		gd.addStringField ("x3d subdirectory prefix",                          this.x3dSubdirPrefix, 10,    // 14a
    				"When using timestamp as a subdirectory, add this prefix");

    		gd.addStringField ("ML output directory",                              this.mlDirectory, 60,
    				"Non-empty directory with no \"/\" separator makes it a subdirectory of the model version directory");
    		gd.addCheckbox    ("Select ML output directory", false,"Erase text field or use \"/\" in it to enable absolute directory path selection");

			gd.addMessage     ("============ Main camera============");

    		gd.addStringField ("Sensor calibration directory",                     this.sensorDirectory, 60);    // 4
    		gd.addCheckbox    ("Select sensor calibration directory",              false);                       // 5
    		gd.addStringField ("Aberration kernels for CLT directory",             this.cltKernelDirectory, 60); // 6
    		gd.addCheckbox    ("Select aberration kernels for CLT directory",      false);                       // 7
    		gd.addStringField ("Results directory",                                this.resultsDirectory, 60);   // 11
    		gd.addCheckbox    ("Select results directory",                         false);                       // 12
    		gd.addNumericField("First subcamera (in the source filename)",         this.firstSubCamera, 0);      // 15
    		gd.addNumericField("First subcamera (in config (clt, sensor) directories)", this.firstSubCameraConfig, 0);
    		gd.addNumericField("Number of subcameras in this camera ",             this.numSubCameras, 0); // 16
    		gd.addStringField ("Source files prefix",                              this.sourcePrefix, 60);       // 13
    		gd.addStringField ("Source files suffix",                              this.sourceSuffix, 60);       // 14
    		gd.addStringField ("Sensor files prefix",                              this.sensorPrefix, 40);       // 17
    		gd.addStringField ("Sensor files suffix",                              this.sensorSuffix, 40);       // 18

    		gd.addStringField ("CLT kernel files prefix",                          this.cltKernelPrefix, 40);    // 19
    		gd.addStringField ("CLT kernel files suffix",                          this.cltSuffix, 40);          // 20
    		gd.addStringField ("x3d subdirectory suffix",                          this.x3dSubdirSuffix, 10,     // 20a
    				"When using timestamp as a subdirectory, add this suffix");

			gd.addMessage     ("============ Auxiliary camera============");
    		gd.addStringField ("Aux sensor calibration directory",                     this.aux_camera.sensorDirectory, 60);    // 4b
    		gd.addCheckbox    ("Select aux sensor calibration directory",              false);                                  // 5b
    		gd.addStringField ("Aberration kernels for aux CLT directory",             this.aux_camera.cltKernelDirectory, 60); // 6b
    		gd.addCheckbox    ("Select aberration kernels for aux CLT directory",      false);                                  // 7b
    		gd.addStringField ("Aux results directory",                                this.aux_camera.resultsDirectory, 60);   // 11b
    		gd.addCheckbox    ("Select aux results directory",                         false);                                  // 12b
    		gd.addNumericField("First aux subcamera (in the source filename)",         this.aux_camera.firstSubCamera, 0);      // 15b
    		gd.addNumericField("First aux subcamera (in config (clt, sensor) directories)",this.aux_camera.firstSubCameraConfig, 0);
    		gd.addNumericField("Number of aux subcameras in this camera ",             this.aux_camera.numSubCameras, 0); // 16b
    		gd.addStringField ("Aux Source files prefix",                              this.aux_camera.sourcePrefix, 60);       // 13
    		gd.addStringField ("Aux Source files suffix",                              this.aux_camera.sourceSuffix, 60);       // 14
    		gd.addStringField ("Aux sensor files prefix",                              this.aux_camera.sensorPrefix, 40);       // 17b
    		gd.addStringField ("Aux sensor files suffix",                              this.aux_camera.sensorSuffix, 40);       // 18b
    		gd.addStringField ("Aux CLT kernel files prefix",                          this.aux_camera.cltKernelPrefix, 40);    // 19b
    		gd.addStringField ("Aux CLT kernel files suffix",                          this.aux_camera.cltSuffix, 40);          // 20b
    		gd.addStringField ("Aux x3d subdirectory suffix",                          this.aux_camera.x3dSubdirSuffix, 10,     // 20ba
    				"When using timestamp as a subdirectory, add this suffix");

  			gd.addTab         ("Batch", "Select Batch parameters");
    		gd.addCheckbox    ("Process main camera source images (false - ignore)",   this.process_main_sources); // 20c
    		gd.addCheckbox    ("Process AUX camera source images (false - ignore)",    this.process_aux_sources); // 20d
  			
			gd.addNumericField("Bitmask of channels were to look for GPS data",                            this.kml_sensors,  0);

			gd.addNumericField("Repeat main camera field adjustment (early, before rig)",                  this.rig_batch_adjust_main,  0);
			gd.addNumericField("Repeat aux camera field adjustment  (early, before rig)",                  this.rig_batch_adjust_aux,   0);
			gd.addNumericField("Repeat 2-quad camera rig field adjustment  (early, before late main/aux)", this.rig_batch_adjust_rig,   0);

			gd.addNumericField("Repeat main camera field adjustment (late, with GT disparity from rig)",   this.rig_batch_adjust_main_gt,  0);
			gd.addNumericField("Repeat aux camera field adjustment (late, with GT disparity from rig)",    this.rig_batch_adjust_aux_gt,   0);
			gd.addNumericField("Repeat 2-quad camera rig field adjustment (late, after all others)",       this.rig_batch_adjust_rig_gt,   0);

    		gd.addCheckbox    ("Experimental DSI",                                                   this.clt_batch_dsi1); // 21
    		gd.addCheckbox    ("Apply (and disable) manual pixel shift",                             this.clt_batch_apply_man); // 21
    		gd.addCheckbox    ("Calibrate extrinsic parameters for each set",                        this.clt_batch_extrinsic); // 22
    		gd.addCheckbox    ("Calculate fine polynomial correction for each set",                  this.clt_batch_poly);      // 23
    		gd.addCheckbox    ("Create a set of 4 images, usually for disparity = 0",                this.clt_batch_4img);      // 24
    		gd.addCheckbox    ("Create a set of 4 images for AUX (LWIR) camera",                     this.clt_batch_4img_aux);      // 24
    		gd.addCheckbox    ("1-st step of 3d reconstruction - explore disparities for each tile", this.clt_batch_explore);   // 25
    		gd.addCheckbox    ("Create super-tile 2.5d surfaces",                                    this.clt_batch_surf);      // 26
    		gd.addCheckbox    ("Assign tiles to surfaces",                                           this.clt_batch_assign);    // 27
    		gd.addCheckbox    ("Generate 3d output: x3d and/or obj+mtl",                             this.clt_batch_gen3d);     // 28
    		gd.addCheckbox    ("Generate ML output files",                                           this.clt_batch_genMl);     // 28

    		gd.addCheckbox    ("Generate debug images if a single set is selected",                  this.clt_batch_dbg1);      // 29

    		gd.addCheckbox    ("Create DSI combo image",                                             this.clt_batch_dsi,
    				"Save main camera, dual-quad rig and optionally aux camera combo DSI image with the model");
    		gd.addCheckbox    ("Include/genarate separate aux camera DSI data in the combo DSI",     this.clt_batch_dsi_aux,
    				"8-rig: DSI for the AUX camera is offset (by the rig baseline) from the main and rig DSI. Aux DSI requires extra processing time."+
    		"EO+LWIR - generate a separate GT+AUX file");
    		gd.addCheckbox    ("Use CM strength (no switch between LMA/no-LMA) for DSI",             this.clt_batch_dsi_cm_strength,
    				"Generate CM-only, single-tile strength for each tile keeping disparity and LMA-disparity from multi-tile"+
    		        "to use as a layer for interscene matching");
    		gd.addCheckbox    ("Additional steps to calculate Aux DSI (more than for LY adjustment)",   this.clt_batch_dsi_aux_full,
    				"(Not yet tested)");
    		gd.addCheckbox    ("Save field adjustment data with the model",                          this.clt_batch_save_extrinsics,
    				"This data can be used to restore specific filed-adjusted cameras extrinsics used when the model was generated");
    		gd.addCheckbox    ("Save all parameters with the model",                                 this.clt_batch_save_all,
    				"Save a copy of all parameters with the model");

			gd.addMessage     ("============ LWIR16 processing ============");
    		gd.addCheckbox    ("Skip scenes processing",                                               this.clt_batch_skip_scenes,
    				"Skip all per-scene processing, go directly to processing sequences");
    		
			gd.addMessage     ("============ RGB cameras ============");
    		gd.addCheckbox    ("RGB: Calculate pair-wise camera poses",                                 this.clt_batch_pose_pairs_main,
    				"RGB: Relative poses are calculated for pairs of consecututive scenes. Requires DSI for each scene");
    		gd.addCheckbox    ("RGB: Scene poses relative to the last",                                 this.clt_batch_pose_last_main,
    				"RGB: Relative camera poses to the reference (last) scene");
    		gd.addCheckbox    ("RGB: Scene poses relative to others",                                   this.clt_batch_pose_scene_main,
    				"RGB: Camera poses relative to all other scenes in the series, not just relative to the latest (not yet implemented)");

  			gd.addNumericField("RGB: Offset latest reference scene",                                    this.clt_batch_offset_main,   0, 3, "scenes",
  					"When selecting multiple reference scenes for ML files generation, offset from the last scene in the series");
  			gd.addNumericField("RGB: Step between reference scenes",                                    this.clt_batch_step_main,   0, 3, "scenes",
  					"When selecting multiple reference scenes for ML files generation, step between scenes");
    		
    		
    		gd.addCheckbox    ("RGB: Generate ML files for the last scene",                             this.clt_batch_ml_last_main,
    				"RGB: Generate ML output files for the last scene, requres 'Scene poses relative to the last'");
    		gd.addCheckbox    ("RGB: Generate ML files for each scene",                                 this.clt_batch_ml_all_main,
    				"RGB: Requires 'Scene poses relative to others', not yet implemented");
			
			gd.addMessage     ("============ LWIR cameras ============");
			gd.addCheckbox    ("LWIR: Calculate pair-wise camera poses",                                 this.clt_batch_pose_pairs_aux,
    				"LWIR: Relative poses are calculated for pairs of consecututive scenes. Requires DSI for each scene");
    		gd.addCheckbox    ("LWIR: Scene poses relative to the last",                                 this.clt_batch_pose_last_aux,
    				"LWIR: Relative camera poses to the reference (last) scene");
    		gd.addCheckbox    ("LWIR: Scene poses relative to others",                                   this.clt_batch_pose_scene_aux,
    				"LWIR: Camera poses relative to all other scenes in the series, not just relative to the latest (not yet implemented)");

  			gd.addNumericField("LWIR: Offset latest reference scene",                                    this.clt_batch_offset_aux,   0, 3, "scenes",
  					"When selecting multiple reference scenes for ML files generation, offset from the last scene in the series");
  			gd.addNumericField("LWIR: Step between reference scenes",                                    this.clt_batch_step_aux,   0, 3, "scenes",
  					"When selecting multiple reference scenes for ML files generation, step between scenes");
    		
    		gd.addCheckbox    ("LWIR: Generate ML files for the last scene",                             this.clt_batch_ml_last_aux,
    				"LWIR: Generate ML output files for the last scene, requres 'Scene poses relative to the last'");
    		gd.addCheckbox    ("LWIR: Generate ML files for each scene",                                 this.clt_batch_ml_all_aux,
    				"LWIR: Requires 'Scene poses relative to others', not yet implemented");

    		if (clt_parameters != null) {
//    			gd.addMessage     ("============ selected CLT parameters ============");
      			gd.addTab         ("CLT", "Modify selected CLT parameters");
    			gd.addNumericField("Maximal disparity to try",                                                            clt_parameters.grow_disp_max,  6);
      			gd.addCheckbox    ("Equalize green channel gain of the individual cnannels (bug fix for exposure)",       clt_parameters.gain_equalize);
      			gd.addNumericField("Inverse distance to infinity (misalignment correction)",                              clt_parameters.z_correction,  6);
      			gd.addNumericField("Number of clusters to keep",                                                          clt_parameters.tsNumClust,  0);
      			gd.addNumericField("Maximal number of output meshes to generate",                                         clt_parameters.max_clusters,   0);

    		}
//    		WindowTools.addScrollBars(gd);
    		gd.showDialog();
    		if (gd.wasCanceled()) return false;

    		this.tile_processor_gpu =    gd.getNextString(); if (gd.getNextBoolean()) selectGPUSourceDirectory(false, false);

    		this.saveSettings=           gd.getNextBoolean(); // 1

    		this.sourceDirectory=        gd.getNextString(); if (gd.getNextBoolean()) selectSourceDirectory(false, false);   // 3
    		this.use_set_dirs =          gd.getNextBoolean();
    		this.x3dModelVersion=        gd.getNextString(); //  10a
    		this.jp4SubDir=              gd.getNextString(); //  10b
    		this.x3dDirectory=           gd.getNextString(); if (gd.getNextBoolean()) selectX3dDirectory(false, true);       // 9
    		this.use_x3d_subdirs=        gd.getNextBoolean(); // 10
//    		this.sourcePrefix=           gd.getNextString();  // 13
//    		this.sourceSuffix=           gd.getNextString();  // 14
    		this.x3dSubdirPrefix=        gd.getNextString();  // 14a
    		this.mlDirectory=            gd.getNextString(); if (gd.getNextBoolean()) selectMlDirectory(null,false, true);       // 8d

// main camera
    		this.sensorDirectory=        gd.getNextString(); if (gd.getNextBoolean()) selectSensorDirectory(false, false);   // 5
    		this.cltKernelDirectory=     gd.getNextString(); if (gd.getNextBoolean()) selectCLTKernelDirectory(false, true); // 7
    		this.resultsDirectory=       gd.getNextString(); if (gd.getNextBoolean()) selectResultsDirectory(false, true);   // 12
    		this.firstSubCamera=   (int) gd.getNextNumber();  // 15
    		this.firstSubCameraConfig=(int) gd.getNextNumber();
    		this.numSubCameras=    (int) gd.getNextNumber();  // 16
    		this.sourcePrefix=           gd.getNextString();  // 13
    		this.sourceSuffix=           gd.getNextString();  // 14
    		this.sensorPrefix=           gd.getNextString();  // 17
    		this.sensorSuffix=           gd.getNextString();  // 18
    		this.cltKernelPrefix=        gd.getNextString();  // 19
    		this.cltSuffix=              gd.getNextString();  // 20
    		this.x3dSubdirSuffix=        gd.getNextString();  // 20a

// aux camera
    		this.aux_camera.sensorDirectory=        gd.getNextString(); if (gd.getNextBoolean()) aux_camera.selectSensorDirectory(false, false);   // 5b
    		this.aux_camera.cltKernelDirectory=     gd.getNextString(); if (gd.getNextBoolean()) aux_camera.selectCLTKernelDirectory(false, true); // 7b
    		this.aux_camera.resultsDirectory=       gd.getNextString(); if (gd.getNextBoolean()) aux_camera.selectResultsDirectory(false, true);   // 12b
    		this.aux_camera.firstSubCamera=   (int) gd.getNextNumber();  // 15b
    		this.aux_camera.firstSubCameraConfig=(int) gd.getNextNumber();
    		this.aux_camera.numSubCameras=    (int) gd.getNextNumber();  // 16b
    		this.aux_camera.sourcePrefix=           gd.getNextString();  // 13
    		this.aux_camera.sourceSuffix=           gd.getNextString();  // 14
    		this.aux_camera.sensorPrefix=           gd.getNextString();  // 17b
    		this.aux_camera.sensorSuffix=           gd.getNextString();  // 18b
    		this.aux_camera.cltKernelPrefix=        gd.getNextString();  // 19b
    		this.aux_camera.cltSuffix=              gd.getNextString();  // 20b
    		this.aux_camera.x3dSubdirSuffix=        gd.getNextString();  // 20ba

    		this.process_main_sources=              gd.getNextBoolean(); // 20c
    		this.process_aux_sources=               gd.getNextBoolean(); // 20d

    		this.kml_sensors =                (int) gd.getNextNumber();

    		this.rig_batch_adjust_main =      (int) gd.getNextNumber();
			this.rig_batch_adjust_aux =       (int) gd.getNextNumber();
			this.rig_batch_adjust_rig =       (int) gd.getNextNumber();

			this.rig_batch_adjust_main_gt =   (int) gd.getNextNumber();
			this.rig_batch_adjust_aux_gt =    (int) gd.getNextNumber();
			this.rig_batch_adjust_rig_gt =    (int) gd.getNextNumber();

    		this.clt_batch_dsi1=         gd.getNextBoolean(); // 21
    		this.clt_batch_apply_man=    gd.getNextBoolean(); // 21
    		this.clt_batch_extrinsic=    gd.getNextBoolean(); // 22
    		this.clt_batch_poly=         gd.getNextBoolean(); // 23
    		this.clt_batch_4img=         gd.getNextBoolean(); // 24
    		this.clt_batch_4img_aux=     gd.getNextBoolean(); // 24
    		this.clt_batch_explore=      gd.getNextBoolean(); // 25
    		this.clt_batch_surf=         gd.getNextBoolean(); // 26
    		this.clt_batch_assign=       gd.getNextBoolean(); // 27
    		this.clt_batch_gen3d=        gd.getNextBoolean(); // 28
    		this.clt_batch_genMl=        gd.getNextBoolean(); // 28
    		this.clt_batch_dbg1=         gd.getNextBoolean(); // 29
    		this.clt_batch_dsi=             gd.getNextBoolean();
    		this.clt_batch_dsi_aux=         gd.getNextBoolean();
    		this.clt_batch_dsi_cm_strength= gd.getNextBoolean();
    		this.clt_batch_dsi_aux_full=    gd.getNextBoolean();
    		this.clt_batch_save_extrinsics= gd.getNextBoolean();
    		this.clt_batch_save_all=        gd.getNextBoolean();
    		
    		this.clt_batch_skip_scenes=        gd.getNextBoolean();

    		this.clt_batch_pose_pairs_main=   gd.getNextBoolean();
    		this.clt_batch_pose_last_main=    gd.getNextBoolean();
    		this.clt_batch_pose_scene_main=   gd.getNextBoolean();
    		this.clt_batch_offset_main= (int) gd.getNextNumber();
    		this.clt_batch_step_main =  (int) gd.getNextNumber();
    		this.clt_batch_ml_last_main=      gd.getNextBoolean();
    		this.clt_batch_ml_all_main=       gd.getNextBoolean();

    		this.clt_batch_pose_pairs_aux=   gd.getNextBoolean();
    		this.clt_batch_pose_last_aux=    gd.getNextBoolean();
    		this.clt_batch_pose_scene_aux=   gd.getNextBoolean();
    		this.clt_batch_offset_aux= (int) gd.getNextNumber();
    		this.clt_batch_step_aux =  (int) gd.getNextNumber();
    		this.clt_batch_ml_last_aux=      gd.getNextBoolean();
    		this.clt_batch_ml_all_aux=       gd.getNextBoolean();
    		
    		if (clt_parameters != null) {
    			clt_parameters.grow_disp_max =      gd.getNextNumber();
    			clt_parameters.gain_equalize =      gd.getNextBoolean();
    			clt_parameters.z_correction =       gd.getNextNumber();
      			clt_parameters.tsNumClust =   (int) gd.getNextNumber();
      			clt_parameters.max_clusters = (int) gd.getNextNumber();
    		}
    		return true;
    	}



// TODO: extract timestamnp from JP4 or, at least combine movie timestamp+frame into a single filename string
    	public String [] getSourcePaths(){
    		String [] empty={};
    		return (this.sourcePaths!=null)?this.sourcePaths:empty;
    	}

    	public static int getChannelFromTiff(String path, String suffix){
    		int indexSuffix=path.length()-suffix.length();
    		int indexLastDash=indexSuffix-1; // in jp4 it will be underscore, not dash? Or should we change that?
    		while ((indexLastDash>0) &&
    				(indexLastDash>(indexSuffix-3)) &&
    				(path.charAt(indexLastDash)!='_') &&
    				(path.charAt(indexLastDash)!='-')) indexLastDash--;
    		return Integer.parseInt(path.substring(indexLastDash+1,indexSuffix));

    	}

    	public static int getChannelFromTiff(String path, String [] suffixes){
    		String suffix = null;
    		for (String s:suffixes) {
    			if (path.endsWith(s)) {
    				suffix = s;
    				break;
    			}
    		}
    		if (suffix == null) {
    			return -1;
    		}
    		int indexSuffix=path.length()-suffix.length();
    		int indexLastDash=indexSuffix-1; // in jp4 it will be underscore, not dash? Or should we change that?
    		while ((indexLastDash>0) &&
    				(indexLastDash>(indexSuffix-3)) &&
    				(path.charAt(indexLastDash)!='_') &&
    				(path.charAt(indexLastDash)!='-')) indexLastDash--;
    		return Integer.parseInt(path.substring(indexLastDash+1,indexSuffix));

    	}

    	// from source file name, if image set dirs are used, it's name will be used instead
    	public String getNameFromTiff(String path){
    		return getNameFromTiff(path, getSourceSuffixes());
    	}

    	public String getNameFromTiff(String path, String suffix){
    		if (use_set_dirs) { // use set name, not the file name
    			return getSetName(path);
    		}
    		int indexSuffix=path.length()-suffix.length();
    		int indexLastDash=indexSuffix-1; // in jp4 it will be underscore, not dash? Or should we change that?
    		while ((indexLastDash>0) &&
    				(indexLastDash>(indexSuffix-3)) &&
    				(path.charAt(indexLastDash)!='_') &&
    				(path.charAt(indexLastDash)!='-')) indexLastDash--;
    		int nameStart=path.lastIndexOf(Prefs.getFileSeparator())+1;
    		return path.substring(nameStart,indexLastDash);
    	}

    	public String getNameFromTiff(String path, String suffixes[]){
    		if (use_set_dirs) { // use set name, not the file name
    			return getSetName(path);
    		}
    		String suffix = null;
    		for (String s:suffixes) {
    			if (path.endsWith(s)) {
    				suffix = s;
    				break;
    			}
    		}
    		if (suffix == null) {
    			return null;
    		}
    		int indexSuffix=path.length()-suffix.length();
    		int indexLastDash=indexSuffix-1; // in jp4 it will be underscore, not dash? Or should we change that?
    		while ((indexLastDash>0) &&
    				(indexLastDash>(indexSuffix-3)) &&
    				(path.charAt(indexLastDash)!='_') &&
    				(path.charAt(indexLastDash)!='-')) indexLastDash--;
    		int nameStart=path.lastIndexOf(Prefs.getFileSeparator())+1;
    		return path.substring(nameStart,indexLastDash);
    	}


    	public boolean isJP4(){
			return this.sourceSuffix.equals(".jp4") || this.sourceSuffix.equals(".jp46");
    	}
    	public int getChannelFromSourceTiff(String path){
    		return getChannelFromTiff(path, getSourceSuffixes());
    	}
    	public String getNameFromSourceTiff(String path){
    		return getNameFromTiff(path, getSourceSuffixes());
    	}

    	public int getChannelFromKernelTiff(String path, int type){return getChannelFromTiff(path, (type==0)?this.sharpKernelSuffix:this.smoothKernelSuffix);}
    	public String getNameFromKernelTiff(String path, int type){return getNameFromTiff(path, (type==0)?this.sharpKernelSuffix:this.smoothKernelSuffix);}

    	public int getChannelFromDCTTiff(String path, int type){return getChannelFromTiff(path, (type==0)?this.dctSymSuffix:this.dctAsymSuffix);}
    	public String getNameFromDCTTiff(String path, int type){return getNameFromTiff(path, (type==0)?this.dctSymSuffix:this.dctAsymSuffix);}

    	public int getChannelFromCLTTiff(String path){return getChannelFromTiff(path, this.cltSuffix);}
    	public String getNameFromCLTTiff(String path){return getNameFromTiff(path, this.cltSuffix);}


    	public boolean selectSourceSets() {
    		return selectSourceSets(1);
    	}

    	public String getSetName(String filePath) {
//    		System.out.println(filePath);
    		if (filePath == null) {
    			return null;
    		}
    		return new File(filePath).getParentFile().getName();
    	}
    	public String getSetPath(String filePath) {
    		if (filePath == null) {
    			return null;
    		}
    		return (new File(filePath)).getParentFile().getPath();
    	}
    	public String [] getSetList(String [] filePaths) {
    		if ((filePaths == null) || (filePaths.length == 0)) {
    			return new String[0];
    		}
    		HashSet<String> sets= new HashSet<String>();
    		for (String fn:filePaths) {
    			String p =  getSetPath(fn);
    			if (p != null) {
    				sets.add(getSetPath(fn));
    			}
    		}
    		ArrayList<String> setList = new ArrayList<String>(sets);
    		Collections.sort(setList);
    		return setList.toArray(new String[0]);
    	}

    	public String [] getSourceSuffixes() {
    		String [] suffixes = null;//    		={this.sourceSuffix};
    		if (aux_camera != null) {
    			suffixes = new String[2];
    			suffixes[1] = aux_camera.sourceSuffix;
    		} else {
    			suffixes = new String[1];
    		}
    		suffixes[0] = sourceSuffix;
    		return suffixes;
    	}

    	public String [] getSourcePrefixes() {
    		String [] prefixes = null;//    		={this.sourceSuffix};
    		if (aux_camera != null) {
    			prefixes = new String[2];
    			prefixes[1] = aux_camera.sourcePrefix;
    		} else {
    			prefixes = new String[1];
    		}
    		prefixes[0] = sourcePrefix;
    		return prefixes;
    	}

    	public boolean selectSourceSets(int debugLevel) {
    		String [] defaultPaths = getSetList(this.sourcePaths); // returns non-null
    		File [] defaultFiles = new File[defaultPaths.length];
    		for (int i = 0; i < defaultPaths.length; i++) {
    			defaultFiles[i] = new File(defaultPaths[i]);
    		}

    		String [] extensions = getSourceSuffixes();//    		={this.sourceSuffix};
    		String [] prefixes = getSourcePrefixes();
    		int num_chn_main = numSubCameras;
    		int num_chn_aux =  ((aux_camera != null)?aux_camera.numSubCameras : 0);
    		int num_chn_files = num_chn_main + num_chn_aux;
    		extensions[0] = sourceSuffix;
    		prefixes[0] =   sourcePrefix;
			MultipleExtensionsFileFilter setFilter = new MultipleExtensionsFileFilter(prefixes,extensions,"Image sets");
			MultipleExtensionsFileFilter setFilterMain = new MultipleExtensionsFileFilter(
					new String[] {prefixes[0]},new String[] {extensions[0]},"Image sets main");
			MultipleExtensionsFileFilter setFilterAux = new MultipleExtensionsFileFilter(
					new String[] {prefixes[1]},new String[] {extensions[1]},"Image sets main");

	    	DirectoryChoser dc = new DirectoryChoser(
	    			setFilter,
	    			num_chn_files,
	    			0, // num_chn_files?
	    			null);
	    	dc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
	    	dc.setMultiSelectionEnabled(true);
	    	dc.setDialogTitle("Select Image sets (directories with simultaneous image files)");
	    	dc.setApproveButtonText("Select");
	    	File cur_dir = new File(this.sourceDirectory);
//	    	if (cur_dir != null) {
//	    		cur_dir = cur_dir.getParentFile();
//	    	}
	    	dc.setCurrentDirectory(cur_dir);
	    	dc.setSelectedFiles(defaultFiles);


	    	int returnVal = dc.showOpenDialog(IJ.getInstance());
	    	if (returnVal!=JFileChooser.APPROVE_OPTION)	return false;
	    	File [] files = dc.getSelectedFiles();
	    	if (files.length<1) return false;
	    	// Can not make it work correctly with multiple selection, giving up for now

	    	ArrayList<File>  setDirList =   new ArrayList<File>(); // list of set directories
	    	ArrayList<File>  setFilesList = new ArrayList<File>(); // list of set files
	    	for (int nFile=0;nFile<files.length;nFile++) {
//	    		String [] setChnFiles = files[nFile].list(setFilter);
	    		File [] setChnFiles =     files[nFile].listFiles(setFilter);
	    		File [] setMainChnFiles = files[nFile].listFiles(setFilterMain);
	    		File [] setAuxChnFiles =  files[nFile].listFiles(setFilterAux);
	    		int num_match = setChnFiles.length;
	    		if (num_match == num_chn_files) {
//	    			|| // all files for main and aux
//	    				(setMainChnFiles.length == num_chn_main) || // has all needed main camera files
//	    				(setAuxChnFiles.length == num_chn_aux))   //   has all needed camera files
//	    				{ // only use sets of exact number of files
	    			setDirList.add(files[nFile]);
	    			for (File f: setChnFiles) {
	    				setFilesList.add(f);
	    			}
	    		} else if ((setMainChnFiles.length == num_chn_main) || (setAuxChnFiles.length == num_chn_aux))  {
	    			setDirList.add(files[nFile]);
	    			if (setMainChnFiles.length == num_chn_main) {
		    			for (File f: setMainChnFiles) setFilesList.add(f);
	    			}
	    			if (setAuxChnFiles.length == num_chn_aux) {
		    			for (File f: setAuxChnFiles) setFilesList.add(f);
	    			}
	    		}
	    	}
	    	String [] sourceSetPaths = new String[setDirList.size()];
	    	for (int nFile = 0; nFile < sourceSetPaths.length; nFile++) {
	    		sourceSetPaths[nFile]= setDirList.get(nFile).getPath();
	    	}

	    	this.sourcePaths = new String[setFilesList.size()];
	    	for (int nFile = 0; nFile < sourcePaths.length; nFile++) {
	    		sourcePaths[nFile]= setFilesList.get(nFile).getPath();
	    	}

	    	if (setDirList.size() >1) {
	    		this.sourceDirectory = setDirList.get(0).getParentFile().getPath();
	    	}
    		return true;
    	}

// get list of source files in a directory
    	public String[]  selectSourceFileInSet(String setdir, int debugLevel) {
    		int num_chn_main = numSubCameras;
    		int num_chn_aux =  ((aux_camera != null)? aux_camera.numSubCameras : 0);
    		int num_chn_files = num_chn_main + num_chn_aux;
    		String [] extensions = getSourceSuffixes();//    		={this.sourceSuffix};
    		String [] prefixes = getSourcePrefixes();
    		extensions[0] = sourceSuffix;
    		prefixes[0] =   sourcePrefix;
			MultipleExtensionsFileFilter setFilter = new MultipleExtensionsFileFilter(prefixes,extensions,"Image sets");
	    	File fsetdir = new File(setdir);
	    	ArrayList<File>  setFilesList = new ArrayList<File>(); // list of set files
	    	File [] setChnFiles = fsetdir.listFiles(setFilter);
	    	int num_match = setChnFiles.length;
	    	if (    (num_match == num_chn_files) || // all files for main and aux
	    			(num_match == num_chn_main) || // only main camera files
	    			(num_match == num_chn_aux))   // only aux camera files
	    	{ // only use sets of exact number of files
	    		for (File f: setChnFiles) {
	    			setFilesList.add(f);
	    		}
	    	}

	    	String [] sourcePaths = new String[setFilesList.size()];
	    	for (int nFile = 0; nFile < sourcePaths.length; nFile++) {
	    		sourcePaths[nFile]= setFilesList.get(nFile).getPath();
	    	}
    		return sourcePaths;
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
    					"Select Source files, saved as "+ extensions[0],
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

    	public String [] selectSensorFiles(int debugLevel) { // will only open dialog if directory or files are not found
    		String []defaultPaths = new String[1];
    		if ((this.sensorDirectory==null) || (this.sensorDirectory.length()<=1)){ // empty or "/"
    			defaultPaths[0]="";
    		} else {
    			defaultPaths[0]=this.sensorDirectory+Prefs.getFileSeparator();
    		}
    		String [] extensions={this.sensorSuffix};
			MultipleExtensionsFileFilter sensorFilter =
				new MultipleExtensionsFileFilter(this.sensorPrefix,extensions,this.sensorPrefix+"*"+extensions[0]+" Sensor calibration files");
			if (debugLevel>1) System.out.println("selectSensorFiles("+debugLevel+"): defaultPaths[0]="+defaultPaths[0]+" "+this.sensorPrefix+"*"+this.sensorSuffix);

			String [] sensorFiles=null;
// try reading all matching files
			File dir= new File (this.sensorDirectory);
//			if (debugLevel>1) System.out.println("selectSensorFiles, dir="+this.sensorDirectory);
			File [] fileList=null;
			if (dir.exists()) {
				fileList=dir.listFiles(sensorFilter);
			}
			if ((fileList==null) || (fileList.length==0)){
				sensorFiles=CalibrationFileManagement.selectFiles(false,
    					"Select sensor calibration files",
    					"Select",
    					sensorFilter,
    					defaultPaths); // String [] defaultPaths); //this.sourceDirectory // null
    			if ((sensorFiles!=null) && (sensorFiles.length>0)){
    				this.sensorDirectory=sensorFiles[0].substring(0, sensorFiles[0].lastIndexOf(Prefs.getFileSeparator()));
    				dir= new File (this.sensorDirectory);
//    				if (debugLevel>1) System.out.println("selectSensorFiles, dir="+this.sensorDirectory);
    				fileList=dir.listFiles(sensorFilter);
    			}
			}
			if ((fileList==null) || (fileList.length==0)) return null;

			if (debugLevel>1) System.out.println("Sensor directory "+this.sensorDirectory+" has "+fileList.length+" matching sensor files.");
			sensorFiles = new String[fileList.length];
			for (int i=0;i<sensorFiles.length;i++) sensorFiles[i]=fileList[i].getPath();

			String directory=sensorFiles[0].substring(0, sensorFiles[0].lastIndexOf(Prefs.getFileSeparator()));
			String prefix=sensorFiles[0].substring(directory.length()+1, sensorFiles[0].length()-extensions[0].length()-2); // all but NN
			this.sensorDirectory=directory;
			this.sensorPrefix=prefix;

			return sensorFiles;
    	}

    	public String selectEquirectangularMapFile(
    			int channel,
    			int debugLevel) { // will only open dialog if directory or files are not found
    		String defaultPath="";
    		if ((this.equirectangularDirectory!=null) && (this.equirectangularDirectory.length()>1)){ // empty or "/"
    			defaultPath=this.equirectangularDirectory+Prefs.getFileSeparator();
    		}
    		String [] extensions={String.format("%02d",channel)+this.equirectangularSuffix}; // looking just for a single map
			MultipleExtensionsFileFilter equirectangularFilter =
				new MultipleExtensionsFileFilter(this.equirectangularPrefix,extensions,
						this.equirectangularPrefix+"*"+extensions[0]+" Equirectangular map for channel "+channel);
			if (debugLevel>1) System.out.println("selectEquirectangularMapFile("+debugLevel+"): defaultPath="+defaultPath+
					" "+this.equirectangularPrefix+"*"+this.equirectangularSuffix);

			String equirectangularFile=null;
// try reading all matching files
			File dir= new File (this.equirectangularDirectory);
//			if (debugLevel>1) System.out.println("selectSensorFiles, dir="+this.sensorDirectory);
			File [] fileList=null;
			if (dir.exists()) {
				fileList=dir.listFiles(equirectangularFilter);
			}
			if ((fileList==null) || (fileList.length==0)){
				equirectangularFile=CalibrationFileManagement.selectFile(false,
    					"Select Equirectangular map file for channel "+channel,
    					"Select",
    					equirectangularFilter,
    					defaultPath); // String [] defaultPaths); //this.sourceDirectory // null
    			if (equirectangularFile!=null) {
    				this.equirectangularDirectory=equirectangularFile.substring(0, equirectangularFile.lastIndexOf(Prefs.getFileSeparator()));
    				this.equirectangularPrefix=equirectangularFile.substring(this.equirectangularDirectory.length()+1, equirectangularFile.length()-extensions[0].length()-2);
    				return equirectangularFile;
    			} else return null;
			}
			if ((fileList==null) || (fileList.length==0)) return null;
			equirectangularFile=fileList[0].getPath();
			this.equirectangularDirectory=equirectangularFile.substring(0, equirectangularFile.lastIndexOf(Prefs.getFileSeparator()));
			this.equirectangularPrefix=equirectangularFile.substring(this.equirectangularDirectory.length()+1, equirectangularFile.length()-extensions[0].length()); // extensions[0] already includes channel
			if (fileList.length>1) {
				String msg = "Multiple files matched, prefix updated to match just the first one - "+equirectangularFile;
				System.out.println("Warning: "+msg);
				IJ.showMessage("Warning",msg);
			}
			if (debugLevel>1) System.out.println("selectEquirectangularMapFile() -> "+ equirectangularFile);
			return equirectangularFile;
    	}

    	public String selectPlaneMapFile(
//    			int channel,
    			int debugLevel) { // will only open dialog if directory or files are not found
    		String defaultPath="";
    		if ((this.equirectangularDirectory!=null) && (this.equirectangularDirectory.length()>1)){ // empty or "/"
    			defaultPath=this.equirectangularDirectory+Prefs.getFileSeparator();
    		}
    		String [] extensions={this.planeMapSuffix}; // looking just for a single map
			MultipleExtensionsFileFilter planeMapFilter =
				new MultipleExtensionsFileFilter(this.planeMapPrefix,extensions,
						this.planeMapPrefix+"*"+extensions[0]+" Plane map (all channels)");
			if (debugLevel>1) System.out.println("selectPlaneMapFile("+debugLevel+"): defaultPath="+defaultPath+
					" "+this.planeMapPrefix+"*"+this.planeMapSuffix);
			String planeMapFile=null;
// try reading all matching files
			File dir= new File (this.equirectangularDirectory);
//			if (debugLevel>1) System.out.println("selectSensorFiles, dir="+this.sensorDirectory);
			File [] fileList=null;
			if (dir.exists()) {
				fileList=dir.listFiles(planeMapFilter);
			}
			if ((fileList==null) || (fileList.length==0)){
				planeMapFile=CalibrationFileManagement.selectFile(false,
    					"SelectPlane map file for all channels",
    					"Select",
    					planeMapFilter,
    					defaultPath); // String [] defaultPaths); //this.sourceDirectory // null
    			if (planeMapFile!=null) {
    				this.equirectangularDirectory=planeMapFile.substring(0, planeMapFile.lastIndexOf(Prefs.getFileSeparator()));
    				this.planeMapPrefix=planeMapFile.substring(this.equirectangularDirectory.length()+1, planeMapFile.length()-extensions[0].length()-2);
    				return planeMapFile;
    			} else return null;
			}
			if ((fileList==null) || (fileList.length==0)) return null;
			planeMapFile=fileList[0].getPath();
			this.equirectangularDirectory=planeMapFile.substring(0, planeMapFile.lastIndexOf(Prefs.getFileSeparator()));
			this.planeMapPrefix=planeMapFile.substring(this.equirectangularDirectory.length()+1, planeMapFile.length()-extensions[0].length()); // extensions[0] already includes channel
			if (fileList.length>1) {
				String msg = "Multiple files matched, prefix updated to match just the first one - "+planeMapFile;
				System.out.println("Warning: "+msg);
				IJ.showMessage("Warning",msg);
			}
			if (debugLevel>1) System.out.println("selectPlaneMapFile() -> "+ planeMapFile);
			return planeMapFile;
    	}



    	public String [] selectKernelChannelFiles(
    			int type,  // 0 - sharp, 1 - smooth
    			int firstChannel,
    			int numChannels, // number of channels
    			int debugLevel) { // will only open dialog if directory or files are not found
    		String [] kernelFiles= selectKernelFiles(
        			type,  // 0 - sharp, 1 - smooth
        			debugLevel);
    		if (kernelFiles==null) return null;
    		String [] channelPaths=new String[numChannels];
    		for (int i=0;i<channelPaths.length;i++)channelPaths[i]=null;
    		for (int fileNum=0;fileNum<kernelFiles.length;fileNum++){
    			int chn=getChannelFromKernelTiff(kernelFiles[fileNum], type) - firstChannel;
    			if ((chn>=0) && (chn<numChannels)){
    				if (channelPaths[chn]==null){ // use first file for channel if there are multiple
    					channelPaths[chn]=kernelFiles[fileNum];
    				} else {
    					if (debugLevel>0) System.out.println("Multiple kernel files for channel "+
    							chn+": "+channelPaths[chn]+" and "+kernelFiles[fileNum]+". Using "+channelPaths[chn]);
    				}
    			}
    		}
    		return channelPaths;
    	}

    	public String [] selectKernelFiles(
    			int type,  // 0 - sharp, 1 - smooth
    			int debugLevel) { // will only open dialog if directory or files are not found
    		String []defaultPaths = new String[1];
    		String kernelDirectory=(type==0)?this.sharpKernelDirectory:this.smoothKernelDirectory;
    		if ((kernelDirectory==null) || (kernelDirectory.length()<=1)){ // empty or "/"
    			defaultPaths[0]="";
    		} else {
    			defaultPaths[0]=kernelDirectory+Prefs.getFileSeparator();
    		}
    		String [] extensions={(type==0)?this.sharpKernelSuffix:this.smoothKernelSuffix};
    		String  kernelPrefix= (type==0)?this.sharpKernelPrefix:this.smoothKernelPrefix;
			MultipleExtensionsFileFilter kernelFilter =
				new MultipleExtensionsFileFilter(kernelPrefix,extensions,kernelPrefix+
						"*"+extensions[0]+" "+((type==0)?"Sharp":"Smooth")+" kernel files");
			if (debugLevel>1) System.out.println("selectKernelFiles("+debugLevel+"): defaultPaths[0]="+defaultPaths[0]+" "+kernelPrefix+"*"+extensions[0]);

			String [] kernelFiles=null;
// try reading all matching files
			File dir= new File (kernelDirectory);
//			if (debugLevel>1) System.out.println("selectSensorFiles, dir="+this.sensorDirectory);
			File [] fileList=null;
			if (dir.exists()) {
				fileList=dir.listFiles(kernelFilter);
			}
			if ((fileList==null) || (fileList.length==0)){
				kernelFiles=CalibrationFileManagement.selectFiles(false,
    					"Select"+((type==0)?"sharp":"smooth")+" kernel files",
    					"Select",
    					kernelFilter,
    					defaultPaths); // String [] defaultPaths); //this.sourceDirectory // null
    			if ((kernelFiles!=null) && (kernelFiles.length>0)){
    				kernelDirectory=kernelFiles[0].substring(0, kernelFiles[0].lastIndexOf(Prefs.getFileSeparator()));
    				dir= new File (kernelDirectory);
//    				if (debugLevel>1) System.out.println("selectSensorFiles, dir="+this.sensorDirectory);
    				fileList=dir.listFiles(kernelFilter);
    				if (type==0) this.sharpKernelDirectory= kernelDirectory;
    				else         this.smoothKernelDirectory=kernelDirectory;
    			}
			}
			if ((fileList==null) || (fileList.length==0)) return null;
			if (debugLevel>1) System.out.println(((type==0)?"Sharp":"Smooth")+" kernel directory "+kernelDirectory+" has "+fileList.length+" matching files.");
			kernelFiles = new String[fileList.length];
			for (int i=0;i<kernelFiles.length;i++) kernelFiles[i]=fileList[i].getPath();
			String directory=kernelFiles[0].substring(0, kernelFiles[0].lastIndexOf(Prefs.getFileSeparator()));
			String prefix=kernelFiles[0].substring(directory.length()+1, kernelFiles[0].length()-extensions[0].length()-2); // all but NN
			if (type==0) this.sharpKernelDirectory=directory;
			else         this.smoothKernelDirectory=directory;
			if (type==0) this.sharpKernelPrefix=prefix;
			else         this.smoothKernelPrefix=prefix;
			return kernelFiles;
    	}

    	public String [] selectDCTChannelFiles(
    			int numChannels, // number of channels
    			int debugLevel) { // will only open dialog if directory or files are not found
    		String [] kernelFiles= selectDCTFiles(
        			debugLevel);
    		if (kernelFiles==null) return null;
    		String [] channelPaths=new String[numChannels];
    		for (int i=0;i<channelPaths.length;i++)channelPaths[i]=null;
    		for (int fileNum=0;fileNum<kernelFiles.length;fileNum++){
    			int chn=getChannelFromDCTTiff(kernelFiles[fileNum], 0); // 1 for asym files
    			if ((chn>=0) && (chn<numChannels)){
    				if (channelPaths[chn]==null){ // use first file for channel if there are multiple
    					channelPaths[chn]=kernelFiles[fileNum];
    				} else {
    					if (debugLevel>0) System.out.println("Multiple kernel files for channel "+
    							chn+": "+channelPaths[chn]+" and "+kernelFiles[fileNum]+". Using "+channelPaths[chn]);
    				}
    			}
    		}
    		return channelPaths;
    	}

    	public String [] selectDCTFiles(
    			int debugLevel) { // will only open dialog if directory or files are not found
    		String []defaultPaths = new String[1];
    		String kernelDirectory=this.dctKernelDirectory;
    		if ((kernelDirectory==null) || (kernelDirectory.length()<=1)){ // empty or "/"
    			defaultPaths[0]="";
    		} else {
    			defaultPaths[0]=kernelDirectory+Prefs.getFileSeparator();
    		}
    		String [] extensions={this.dctSymSuffix};
    		String  kernelPrefix= this.dctKernelPrefix;
			MultipleExtensionsFileFilter kernelFilter =
				new MultipleExtensionsFileFilter(kernelPrefix,extensions,kernelPrefix+
						"*"+extensions[0]+" DCT symmetrical kernel files");
			if (debugLevel>1) System.out.println("selectKernelFiles("+debugLevel+"): defaultPaths[0]="+defaultPaths[0]+" "+kernelPrefix+"*"+extensions[0]);

			String [] kernelFiles=null;
// try reading all matching files
			File dir= new File (kernelDirectory);
//			if (debugLevel>1) System.out.println("selectSensorFiles, dir="+this.sensorDirectory);
			File [] fileList=null;
			if (dir.exists()) {
				fileList=dir.listFiles(kernelFilter);
			}
			if ((fileList==null) || (fileList.length==0)){
				kernelFiles=CalibrationFileManagement.selectFiles(false,
    					"Select DCT symmetrical kernel files",
    					"Select",
    					kernelFilter,
    					defaultPaths); // String [] defaultPaths); //this.sourceDirectory // null
    			if ((kernelFiles!=null) && (kernelFiles.length>0)){
    				kernelDirectory=kernelFiles[0].substring(0, kernelFiles[0].lastIndexOf(Prefs.getFileSeparator()));
    				dir= new File (kernelDirectory);
//    				if (debugLevel>1) System.out.println("selectSensorFiles, dir="+this.sensorDirectory);
    				fileList=dir.listFiles(kernelFilter);
    				this.dctKernelDirectory= kernelDirectory;
    			}
			}
			if ((fileList==null) || (fileList.length==0)) return null;
			if (debugLevel>1) System.out.println("DCT kernel directory "+kernelDirectory+" has "+fileList.length+" matching files.");
			kernelFiles = new String[fileList.length];
			for (int i=0;i<kernelFiles.length;i++) kernelFiles[i]=fileList[i].getPath();
			String directory=kernelFiles[0].substring(0, kernelFiles[0].lastIndexOf(Prefs.getFileSeparator()));
			String prefix=kernelFiles[0].substring(directory.length()+1, kernelFiles[0].length()-extensions[0].length()-2); // all but NN
			this.dctKernelDirectory=directory;
			this.dctKernelPrefix=prefix;
			return kernelFiles;
    	}



    	public String [] selectCLTChannelFiles(
    			int firstChannel,
    			int numChannels, // number of channels
    			int debugLevel) { // will only open dialog if directory or files are not found
    		String [] kernelFiles= selectCLTFiles(
        			debugLevel);
    		if (kernelFiles==null) return null;
    		String [] channelPaths=new String[numChannels];
    		for (int i=0;i<channelPaths.length;i++)channelPaths[i]=null;
    		for (int fileNum=0;fileNum<kernelFiles.length;fileNum++){
    			int chn=getChannelFromCLTTiff(kernelFiles[fileNum]) - firstChannel;
    			if ((chn>=0) && (chn<numChannels)){
    				if (channelPaths[chn]==null){ // use first file for channel if there are multiple
    					channelPaths[chn]=kernelFiles[fileNum];
    				} else {
    					if (debugLevel>0) System.out.println("Multiple kernel files for channel "+
    							chn+": "+channelPaths[chn]+" and "+kernelFiles[fileNum]+". Using "+channelPaths[chn]);
    				}
    			}
    		}
    		return channelPaths;
    	}

    	public String [] selectCLTFiles(
    			int debugLevel) { // will only open dialog if directory or files are not found
    		String []defaultPaths = new String[1];
    		String kernelDirectory=this.cltKernelDirectory;
    		if ((kernelDirectory==null) || (kernelDirectory.length()<=1)){ // empty or "/"
    			defaultPaths[0]="";
    		} else {
    			defaultPaths[0]=kernelDirectory+Prefs.getFileSeparator();
    		}
    		String [] extensions={this.cltSuffix};
    		String  kernelPrefix= this.cltKernelPrefix;
			MultipleExtensionsFileFilter kernelFilter =
				new MultipleExtensionsFileFilter(kernelPrefix,extensions,kernelPrefix+
						"*"+extensions[0]+" CLT kernel files");
			if (debugLevel>1) System.out.println("selectKernelFiles("+debugLevel+"): defaultPaths[0]="+defaultPaths[0]+" "+kernelPrefix+"*"+extensions[0]);

			String [] kernelFiles=null;
// try reading all matching files
			File dir= new File (kernelDirectory);
			File [] fileList=null;
			if (dir.exists()) {
				fileList=dir.listFiles(kernelFilter);
			}
			if ((fileList==null) || (fileList.length==0)){
				kernelFiles=CalibrationFileManagement.selectFiles(false,
    					"Select CLT kernel files",
    					"Select",
    					kernelFilter,
    					defaultPaths); // String [] defaultPaths); //this.sourceDirectory // null
    			if ((kernelFiles!=null) && (kernelFiles.length>0)){
    				kernelDirectory=kernelFiles[0].substring(0, kernelFiles[0].lastIndexOf(Prefs.getFileSeparator()));
    				dir= new File (kernelDirectory);
    				fileList=dir.listFiles(kernelFilter);
    				this.cltKernelDirectory= kernelDirectory;
    			}
			}
			if ((fileList==null) || (fileList.length==0)) return null;
			if (debugLevel>1) System.out.println("CLT kernel directory "+kernelDirectory+" has "+fileList.length+" matching files.");
			kernelFiles = new String[fileList.length];
			for (int i=0;i<kernelFiles.length;i++) kernelFiles[i]=fileList[i].getPath();
			String directory=kernelFiles[0].substring(0, kernelFiles[0].lastIndexOf(Prefs.getFileSeparator()));
			String prefix=kernelFiles[0].substring(directory.length()+1, kernelFiles[0].length()-extensions[0].length()-2); // all but NN
			this.cltKernelDirectory=directory;
			this.cltKernelPrefix=prefix;
			return kernelFiles;
    	}

    	public String selectSourceDirectory(boolean smart, boolean newAllowed) { // normally newAllowed=false
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save
    				"Source (acquired from the camera) image directory", // title
    				"Select source directory", // button
    				null, // filter
    				this.sourceDirectory); // this.sourceDirectory);
    		if (dir!=null) this.sourceDirectory=dir;
    		return dir;
    	}
    	public String selectGPUSourceDirectory(boolean smart, boolean newAllowed) { // normally newAllowed=false
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save
    				"GPU kernel development project", // title
    				"Select GPU project directory", // button
    				null, // filter
    				this.tile_processor_gpu); // this.sourceDirectory);
    		if (dir!=null) this.tile_processor_gpu=dir;
    		return dir;
    	}
    	public String selectSensorDirectory(boolean smart,  boolean newAllowed) {
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save
    				"Sensor calibration directory", // title
    				"Select sensor calibration directory", // button
    				null, // filter
    				this.sensorDirectory); //this.sourceDirectory);
    		if (dir!=null) this.sensorDirectory=dir;
    		return dir;
    	}
    	public String selectSharpKernelDirectory(boolean smart, boolean newAllowed) {
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save
    				"Aberration kernels (sharp) directory", // title
    				"Select aberration kernels (sharp) directory", // button
    				null, // filter
    				this.sharpKernelDirectory); //this.sourceDirectory);
    		if (dir!=null) this.sharpKernelDirectory=dir;
    		return dir;
    	}

    	public String selectSmoothKernelDirectory(boolean smart, boolean newAllowed) {
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save
    				"Aberration kernels (smooth) directory", // title
    				"Select aberration kernels (smooth) directory", // button
    				null, // filter
    				this.smoothKernelDirectory); //this.sourceDirectory);
    		if (dir!=null) this.smoothKernelDirectory=dir;
    		return dir;
    	}

    	public String selectDCTKernelDirectory(boolean smart, boolean newAllowed) {
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save
    				"DCT aberration kernels directory (sym and asym files)", // title
    				"Select DCT aberration kernel sdirectory", // button
    				null, // filter
    				this.dctKernelDirectory); //this.sourceDirectory);
    		if (dir!=null) this.dctKernelDirectory=dir;
    		return dir;
    	}

    	public String selectCLTKernelDirectory(boolean smart, boolean newAllowed) {
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save
    				"CLT aberration kernels directory", // title
    				"Select CLT aberration kernels directory", // button
    				null, // filter
    				this.cltKernelDirectory); //this.sourceDirectory);
    		if (dir!=null) this.cltKernelDirectory=dir;
    		return dir;
    	}

    	public String selectX3dDirectory(boolean smart, boolean newAllowed) {
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save
    				"x3d output directory", // title
    				"Select x3d output directory", // button
    				null, // filter
    				this.x3dDirectory); //this.sourceDirectory);
    		if (dir!=null) this.x3dDirectory=dir;
    		return dir;
    	}

    	public String selectMlDirectory(String name, boolean smart, boolean newAllowed) {
    		if ((name != null) && (this.mlDirectory.length()>0) && (!this.mlDirectory.contains(Prefs.getFileSeparator()))) {
    			// relative to the X3D model version
    			String x3d_version_dir = selectX3dDirectory(name, this.x3dModelVersion, smart, newAllowed);
    			if (x3d_version_dir != null) {
    				return x3d_version_dir + Prefs.getFileSeparator() + this.mlDirectory;
    			}
    		}
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save
    				"ML output directory", // title
    				"Select ML output directory", // button
    				null, // filter
    				this.mlDirectory); //this.sourceDirectory);
    		if (dir!=null) this.mlDirectory=dir;
    		return dir;
    	}
    	// add prefix/suffix to the model name
    	public String getModelName(String name) {
    		return this.x3dSubdirPrefix + name + this.x3dSubdirSuffix;
    	}

    	// select qualified (by 'name' - quad timestamp) x3d subdirectory

    	public String selectX3dDirectory(String name, String version, boolean smart, boolean newAllowed) {

    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save
    				"x3d output directory", // title
    				"Select x3d output directory", // button
    				null, // filter
    				this.x3dDirectory); //this.sourceDirectory);
    		if (dir!=null) {
    			this.x3dDirectory=dir;
    			if (this.use_x3d_subdirs && (name != null) && !name.equals("")) {
//    				name = this.x3dDirectory + Prefs.getFileSeparator(); // +this.x3dSubdirPrefix + name + this.x3dSubdirSuffix;
    				name = this.x3dDirectory + Prefs.getFileSeparator()+ name ;
        			if ((version != null) && !version.equals("")) {
        				name = name + Prefs.getFileSeparator()+version;
        			}
    				dir= CalibrationFileManagement.selectDirectory(
    						smart,
    						newAllowed, // save
    						"x3d output sub-directory", // title
    						"Select x3d output sub-directory", // button
    						null, // filter
    						name); //this.x3dDirectory + Prefs.getFileSeparator()+name); //this.sourceDirectory);
    			}
    		}
    		return dir;
    	}

    	public String selectEquirectangularDirectory(boolean smart, boolean newAllowed) {
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save
    				"Equirectangular maps directory", // title
    				"Select equirectangular maps directory", // button
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
     }


    /* === Parameter classes === */
    public static class ProcessParameters {
  	    public int numEyesisChannels=3;
  	    public int numEyesisSubChannels=3;
  	    public boolean eyesisMode;
  		public boolean [][] frames=new boolean[3][3];
  		public boolean selectFile;
  		public boolean thisFileOnly;
  		public int     subChannelToProcess;
  		public boolean split;
  		public boolean debayer;
  		public boolean showDebayerEnergy;
  		public boolean saveDebayerEnergy;
  		public boolean deconvolve;
  		public boolean combine;
  		public boolean showDenoiseMask;
  		public boolean saveDenoiseMask;
  		public boolean showChromaDenoiseMask;
  		public boolean saveChromaDenoiseMask;
  		public boolean showNoiseGains;
  		public boolean saveNoiseGains;
  		public boolean colorProc;
  		public boolean blueProc;
  		public boolean toRGB;
  		public boolean rotate;
  		public boolean crop;   // crop to the sennor size
  		public boolean jpeg;   // convert to RGB and save jpeg (if save is true)
  		public boolean save;
  		public boolean save16; // save 16-bit tiff also if the end result is 8 bit
  		public boolean save32; // save 32-bit tiff also if the end result is 8 or 16 bit
  		public boolean show;
  		public int     JPEG_quality;
  		public double  JPEG_scale;
  		public boolean saveSettings;

  		public ProcessParameters(
  			boolean eyesisMode,
  			boolean frames_11,
  			boolean frames_12,
  			boolean frames_13,
  			boolean frames_21,
  			boolean frames_22,
  			boolean frames_23,
  			boolean frames_31,
  			boolean frames_32,
  			boolean frames_33,
  			boolean selectFile, // ask for file(s) to process
  			boolean thisFileOnly,
  			int     subChannelToProcess,
  			boolean split,
  			boolean debayer,
  			boolean showDebayerEnergy,
  			boolean saveDebayerEnergy,
  			boolean deconvolve,
  			boolean combine,
  			boolean showDenoiseMask,
  			boolean saveDenoiseMask,
  			boolean showChromaDenoiseMask,
  			boolean saveChromaDenoiseMask,
  			boolean showNoiseGains,
  			boolean saveNoiseGains,
  			boolean colorProc,
  			boolean blueProc,
  			boolean toRGB,
  			boolean rotate,
  			boolean crop,   // crop to the sennor size
  			boolean jpeg,   // convert to RGB and save jpeg (if save is true)
  			boolean save,
  			boolean save16, // save 16-bit tiff also if the end result is 8 bit
  			boolean save32, // save 32-bit tiff also if the end result is 8 or 16 bit
  			boolean show,
  			int     JPEG_quality,
  			double  JPEG_scale,
  			boolean saveSettings
  		) {
  			this.eyesisMode=eyesisMode;
  			this.frames[0][0]=frames_11;
  			this.frames[0][1]=frames_12;
  			this.frames[0][2]=frames_13;
  			this.frames[1][0]=frames_21;
  			this.frames[1][1]=frames_22;
  			this.frames[1][2]=frames_23;
  			this.frames[2][0]=frames_31;
  			this.frames[2][1]=frames_32;
  			this.frames[2][2]=frames_33;
  			this.selectFile=selectFile;
  			this.thisFileOnly=thisFileOnly;
  			this.subChannelToProcess=subChannelToProcess;
  			this.split=split;
  			this.debayer=debayer;
  			this.showDebayerEnergy=showDebayerEnergy;
  			this.saveDebayerEnergy=saveDebayerEnergy;
  			this.deconvolve=deconvolve;
  			this.combine=combine;
  			this.showDenoiseMask=showDenoiseMask;
  			this.saveDenoiseMask=saveDenoiseMask;
  			this.showNoiseGains=showNoiseGains;
  			this.saveNoiseGains=saveNoiseGains;
  			this.showChromaDenoiseMask=showChromaDenoiseMask;
  			this.saveChromaDenoiseMask=saveChromaDenoiseMask;
  			this.colorProc=colorProc;
  			this.blueProc=blueProc;
  			this.toRGB=toRGB;
  			this.rotate=rotate;
  			this.crop=crop;
  			this.jpeg=jpeg;
  			this.save=save;
  			this.save16=save16;
  			this.save32=save32;
  			this.show=show;
  			this.JPEG_quality=JPEG_quality;
  			this.JPEG_scale=  JPEG_scale;
  			this.saveSettings=saveSettings;
  		}
  		public void setProperties(String prefix,Properties properties){
  			int i,j;
  			properties.setProperty(prefix+"numEyesisChannels",this.numEyesisChannels+"");
  			properties.setProperty(prefix+"numEyesisSubChannels",this.numEyesisSubChannels+"");
  			properties.setProperty(prefix+"eyesisMode",this.eyesisMode+"");
  		    for (i=0;i<this.frames.length;i++) for (j=0;j<this.frames[0].length;j++)
  				properties.setProperty(prefix+"frames_"+i+"_"+j,this.frames[i][j]+"");
  			properties.setProperty(prefix+"selectFile",this.selectFile+"");
  			properties.setProperty(prefix+"thisFileOnly",this.thisFileOnly+"");
  			properties.setProperty(prefix+"subChannelToProcess",this.subChannelToProcess+"");
  			properties.setProperty(prefix+"split",this.split+"");
  			properties.setProperty(prefix+"debayer",this.debayer+"");
  			properties.setProperty(prefix+"showDebayerEnergy",this.showDebayerEnergy+"");
  			properties.setProperty(prefix+"saveDebayerEnergy",this.saveDebayerEnergy+"");
  			properties.setProperty(prefix+"deconvolve",this.deconvolve+"");
  			properties.setProperty(prefix+"combine",this.combine+"");
  			properties.setProperty(prefix+"showDenoiseMask",this.showDenoiseMask+"");
  			properties.setProperty(prefix+"saveDenoiseMask",this.saveDenoiseMask+"");
  			properties.setProperty(prefix+"showChromaDenoiseMask",this.showChromaDenoiseMask+"");
  			properties.setProperty(prefix+"saveChromaDenoiseMask",this.saveChromaDenoiseMask+"");
  			properties.setProperty(prefix+"showNoiseGains",this.showNoiseGains+"");
  			properties.setProperty(prefix+"saveNoiseGains",this.saveNoiseGains+"");
  			properties.setProperty(prefix+"colorProc",this.colorProc+"");
			properties.setProperty(prefix+"blueProc",this.blueProc+"");
  			properties.setProperty(prefix+"toRGB",this.toRGB+"");
  			properties.setProperty(prefix+"rotate",this.rotate+"");
  			properties.setProperty(prefix+"crop",this.crop+"");
  			properties.setProperty(prefix+"jpeg",this.jpeg+"");
  			properties.setProperty(prefix+"save",this.save+"");
  			properties.setProperty(prefix+"save16",this.save16+"");
  			properties.setProperty(prefix+"save32",this.save32+"");
  			properties.setProperty(prefix+"show",this.show+"");
  			properties.setProperty(prefix+"JPEG_quality",this.JPEG_quality+"");
  			properties.setProperty(prefix+"JPEG_scale",this.JPEG_scale+"");
  			properties.setProperty(prefix+"saveSettings",this.saveSettings+"");
  		}
  		public void getProperties(String prefix,Properties properties){
  			int i,j;
  			if (properties.getProperty(prefix+"numEyesisChannels")!=null) this.numEyesisChannels=Integer.parseInt(properties.getProperty(prefix+"numEyesisChannels"));
  			if (properties.getProperty(prefix+"numEyesisSubChannels")!=null) this.numEyesisSubChannels=Integer.parseInt(properties.getProperty(prefix+"numEyesisSubChannels"));
  			if (properties.getProperty(prefix+"eyesisMode")!=null) this.eyesisMode=Boolean.parseBoolean(properties.getProperty(prefix+"eyesisMode"));
  		    for (i=0;i<this.frames.length;i++) for (j=0;j<this.frames[0].length;j++)
  		    	if (properties.getProperty(prefix+"frames_"+i+"_"+j)!=null) this.frames[i][j]=Boolean.parseBoolean(properties.getProperty(prefix+"frames_"+i+"_"+j));
  		    if (properties.getProperty(prefix+"selectFile")!=null) this.selectFile=Boolean.parseBoolean(properties.getProperty(prefix+"selectFile"));
  		    if (properties.getProperty(prefix+"thisFileOnly")!=null) this.thisFileOnly=Boolean.parseBoolean(properties.getProperty(prefix+"thisFileOnly"));
  		    if (properties.getProperty(prefix+"subChannelToProcess")!=null) this.subChannelToProcess=Integer.parseInt(properties.getProperty(prefix+"subChannelToProcess"));
  		    if (properties.getProperty(prefix+"split")!=null) this.split=Boolean.parseBoolean(properties.getProperty(prefix+"split"));
  		    if (properties.getProperty(prefix+"debayer")!=null) this.debayer=Boolean.parseBoolean(properties.getProperty(prefix+"debayer"));
  		    if (properties.getProperty(prefix+"showDebayerEnergy")!=null) this.showDebayerEnergy=Boolean.parseBoolean(properties.getProperty(prefix+"showDebayerEnergy"));
  		    if (properties.getProperty(prefix+"saveDebayerEnergy")!=null) this.saveDebayerEnergy=Boolean.parseBoolean(properties.getProperty(prefix+"saveDebayerEnergy"));
  		    if (properties.getProperty(prefix+"deconvolve")!=null) this.deconvolve=Boolean.parseBoolean(properties.getProperty(prefix+"deconvolve"));
  		    if (properties.getProperty(prefix+"combine")!=null) this.combine=Boolean.parseBoolean(properties.getProperty(prefix+"combine"));
  		    if (properties.getProperty(prefix+"showDenoiseMask")!=null) this.showDenoiseMask=Boolean.parseBoolean(properties.getProperty(prefix+"showDenoiseMask"));
  		    if (properties.getProperty(prefix+"saveDenoiseMask")!=null) this.saveDenoiseMask=Boolean.parseBoolean(properties.getProperty(prefix+"saveDenoiseMask"));
  		    if (properties.getProperty(prefix+"showChromaDenoiseMask")!=null) this.showChromaDenoiseMask=Boolean.parseBoolean(properties.getProperty(prefix+"showChromaDenoiseMask"));
  		    if (properties.getProperty(prefix+"saveChromaDenoiseMask")!=null) this.saveChromaDenoiseMask=Boolean.parseBoolean(properties.getProperty(prefix+"saveChromaDenoiseMask"));
  		    if (properties.getProperty(prefix+"showNoiseGains")!=null) this.showNoiseGains=Boolean.parseBoolean(properties.getProperty(prefix+"showNoiseGains"));
  		    if (properties.getProperty(prefix+"saveNoiseGains")!=null) this.saveNoiseGains=Boolean.parseBoolean(properties.getProperty(prefix+"saveNoiseGains"));
  		    if (properties.getProperty(prefix+"colorProc")!=null) this.colorProc=Boolean.parseBoolean(properties.getProperty(prefix+"colorProc"));
  		    if (properties.getProperty(prefix+"blueProc")!=null) this.blueProc=Boolean.parseBoolean(properties.getProperty(prefix+"blueProc"));
  		    if (properties.getProperty(prefix+"toRGB")!=null) this.toRGB=Boolean.parseBoolean(properties.getProperty(prefix+"toRGB"));
  		    if (properties.getProperty(prefix+"rotate")!=null) this.rotate=Boolean.parseBoolean(properties.getProperty(prefix+"rotate"));
  		    if (properties.getProperty(prefix+"crop")!=null) this.crop=Boolean.parseBoolean(properties.getProperty(prefix+"crop"));   // crop to the sensor size
  		    if (properties.getProperty(prefix+"jpeg")!=null) this.jpeg=Boolean.parseBoolean(properties.getProperty(prefix+"jpeg"));   // convert to RGB and save jpeg (if save is true)
  		    if (properties.getProperty(prefix+"save")!=null) this.save=Boolean.parseBoolean(properties.getProperty(prefix+"save"));
  		    if (properties.getProperty(prefix+"save16")!=null) this.save16=Boolean.parseBoolean(properties.getProperty(prefix+"save16")); // save 16-bit tiff also if the end result is 8 bit
  		    if (properties.getProperty(prefix+"save32")!=null) this.save32=Boolean.parseBoolean(properties.getProperty(prefix+"save32")); // save 32-bit tiff also if the end result is 8 or 16 bit
  		    if (properties.getProperty(prefix+"show")!=null) this.show=Boolean.parseBoolean(properties.getProperty(prefix+"show"));
  		    if (properties.getProperty(prefix+"JPEG_quality")!=null) this.JPEG_quality=Integer.parseInt(properties.getProperty(prefix+"JPEG_quality"));
  		    if (properties.getProperty(prefix+"JPEG_scale")!=null) this.JPEG_scale=Double.parseDouble(properties.getProperty(prefix+"JPEG_scale"));
  		    if (properties.getProperty(prefix+"saveSettings")!=null) this.saveSettings=Boolean.parseBoolean(properties.getProperty(prefix+"saveSettings"));
  		}
  	}

  /* ======================================================================== */
    public static class FilesParameters {
  	  public String [][] rPSFNames=new String [3][3];
  	  public String [][] gaussianNames=new String [3][3];
  	  public String kernelDirectory;
  	  public String resultsDirectory;
  	  public String [] sourceFiles;
  	  public boolean useXML;
  	  public FilesParameters(
  			  String rPSFNames_11,
  			  String rPSFNames_12,
  			  String rPSFNames_13,
  			  String rPSFNames_21,
  			  String rPSFNames_22,
  			  String rPSFNames_23,
  			  String rPSFNames_31,
  			  String rPSFNames_32,
  			  String rPSFNames_33,
  			  String gaussianNames_11,
  			  String gaussianNames_12,
  			  String gaussianNames_13,
  			  String gaussianNames_21,
  			  String gaussianNames_22,
  			  String gaussianNames_23,
  			  String gaussianNames_31,
  			  String gaussianNames_32,
  			  String gaussianNames_33,
  			  String kernelDirectory,
    		      String resultsDirectory,
    		      boolean useXML
    		      ){
  		  this.rPSFNames[0][0]=rPSFNames_11;
  		  this.rPSFNames[0][1]=rPSFNames_12;
  		  this.rPSFNames[0][2]=rPSFNames_13;
  		  this.rPSFNames[1][0]=rPSFNames_21;
  		  this.rPSFNames[1][1]=rPSFNames_22;
  		  this.rPSFNames[1][2]=rPSFNames_23;
  		  this.rPSFNames[2][0]=rPSFNames_31;
  		  this.rPSFNames[2][1]=rPSFNames_32;
  		  this.rPSFNames[2][2]=rPSFNames_33;
  		  this.gaussianNames[0][0]=gaussianNames_11;
  		  this.gaussianNames[0][1]=gaussianNames_12;
  		  this.gaussianNames[0][2]=gaussianNames_13;
  		  this.gaussianNames[1][0]=gaussianNames_21;
  		  this.gaussianNames[1][1]=gaussianNames_22;
  		  this.gaussianNames[1][2]=gaussianNames_23;
  		  this.gaussianNames[2][0]=gaussianNames_31;
  		  this.gaussianNames[2][1]=gaussianNames_32;
  		  this.gaussianNames[2][2]=gaussianNames_33;
  		  this.kernelDirectory=    kernelDirectory;
  		  this.resultsDirectory=   resultsDirectory;
  		  this.useXML=useXML;
  	  }
  		public void setProperties(String prefix,Properties properties){
//  			properties.setProperty(prefix+"",this.+"");
  			int i,j;
  			for (i=0;i<this.rPSFNames.length;i++) for (j=0;j<this.rPSFNames[i].length;j++)
  				properties.setProperty(prefix+"rPSFNames_"+i+"_"+j,this.rPSFNames[i][j]);
  			for (i=0;i<this.gaussianNames.length;i++) for (j=0;j<this.gaussianNames[i].length;j++)
  				properties.setProperty(prefix+"gaussianNames_"+i+"_"+j,this.gaussianNames[i][j]);
  			properties.setProperty(prefix+"kernelDirectory",this.kernelDirectory);
  			properties.setProperty(prefix+"resultsDirectory",this.resultsDirectory);
  			properties.setProperty(prefix+"useXML",this.useXML+"");
  			j=(this.sourceFiles==null)?0:this.sourceFiles.length;
  			properties.setProperty(prefix+"sourceFiles_length",j+"");
  			for (i=0;i<j;i++)
  				properties.setProperty(prefix+"sourceFiles_"+i,this.sourceFiles[i]);
  		}
  		public void getProperties(String prefix,Properties properties){
  			int i,j;
  			for (i=0;i<this.rPSFNames.length;i++) for (j=0;j<this.rPSFNames[i].length;j++)
  				this.rPSFNames[i][j]=properties.getProperty(prefix+"rPSFNames_"+i+"_"+j);
  			for (i=0;i<this.gaussianNames.length;i++) for (j=0;j<this.gaussianNames[i].length;j++)
  				this.gaussianNames[i][j]=properties.getProperty(prefix+"gaussianNames_"+i+"_"+j);
  			this.kernelDirectory=properties.getProperty(prefix+"kernelDirectory");
  			this.resultsDirectory=properties.getProperty(prefix+"resultsDirectory");
  			this.useXML=Boolean.parseBoolean(properties.getProperty(prefix+"useXML"));
  			j=Integer.parseInt(properties.getProperty(prefix+"sourceFiles_length"));
  			this.sourceFiles=new String[j];
  			for (i=0;i<j;i++)
  				this.sourceFiles[i]=properties.getProperty(prefix+"sourceFiles_"+i);
  		}
    }

  /* ======================================================================== */
    public static class RGBParameters {
  		public double r_min = 0.075;
  		public double g_min = 0.075;
  		public double b_min = 0.075;
  		public double r_max = 1.0;
  		public double g_max = 1.0;
  		public double b_max = 1.0;
  		public double alpha_min = 0.0;
  		public double alpha_max = 1.0;

  		public RGBParameters(double r_min, double g_min, double b_min, double r_max, double g_max, double b_max, double alpha_min, double alpha_max) {
  			this.r_min = r_min;
  			this.g_min = g_min;
  			this.b_min = b_min;
  			this.r_max = r_max;
  			this.g_max = g_max;
  			this.b_max = b_max;
  			this.alpha_min = alpha_min;
  			this.alpha_max = alpha_max;
  		}
  		public void setProperties(String prefix,Properties properties){
  			properties.setProperty(prefix+"r_min",this.r_min+"");
  			properties.setProperty(prefix+"g_min",this.g_min+"");
  			properties.setProperty(prefix+"b_min",this.b_min+"");
  			properties.setProperty(prefix+"r_max",this.r_max+"");
  			properties.setProperty(prefix+"g_max",this.g_max+"");
  			properties.setProperty(prefix+"b_max",this.b_max+"");
  			properties.setProperty(prefix+"alpha_min",this.alpha_min+"");
  			properties.setProperty(prefix+"alpha_max",this.alpha_max+"");
  		}
  		public void getProperties(String prefix,Properties properties){
  			this.r_min=Double.parseDouble(properties.getProperty(prefix+"r_min"));
  			this.g_min=Double.parseDouble(properties.getProperty(prefix+"g_min"));
  			this.b_min=Double.parseDouble(properties.getProperty(prefix+"b_min"));
  			this.r_max=Double.parseDouble(properties.getProperty(prefix+"r_max"));
  			this.g_max=Double.parseDouble(properties.getProperty(prefix+"g_max"));
  			this.b_max=Double.parseDouble(properties.getProperty(prefix+"b_max"));
  			if (properties.getProperty(prefix+"alpha_min")!=null)  this.alpha_min=Double.parseDouble(properties.getProperty(prefix+"alpha_min"));
  			if (properties.getProperty(prefix+"alpha_max")!=null)  this.alpha_max=Double.parseDouble(properties.getProperty(prefix+"alpha_max"));
  		}

  	}
  /* ======================================================================== */

    /* ======================================================================== */
  // individual per-channel color balance and gain
    public static class ColorCalibParameters {
  		public double[][] gain=new double[3][3];
  		public double[][] balanceRed=new double[3][3];
  		public double[][] balanceBlue=new double[3][3];

  	  	public ColorCalibParameters(
  	  			double gain_11,
  	  			double gain_12,
  	  			double gain_13,
  	  			double gain_21,
  	  			double gain_22,
  	  			double gain_23,
  	  			double gain_31,
  	  			double gain_32,
  	  			double gain_33,
  	  			double balanceRed_11,
  	  			double balanceRed_12,
  	  			double balanceRed_13,
  	  			double balanceRed_21,
  	  			double balanceRed_22,
  	  			double balanceRed_23,
  	  			double balanceRed_31,
  	  			double balanceRed_32,
  	  			double balanceRed_33,
  	  			double balanceBlue_11,
  	  			double balanceBlue_12,
  	  			double balanceBlue_13,
  	  			double balanceBlue_21,
  	  			double balanceBlue_22,
  	  			double balanceBlue_23,
  	  			double balanceBlue_31,
  	  			double balanceBlue_32,
  	  			double balanceBlue_33){
  	  		this.gain[0][0]=gain_11;
  	  		this.gain[0][1]=gain_12;
  	  		this.gain[0][2]=gain_13;
  	  		this.gain[1][0]=gain_21;
  	  		this.gain[1][1]=gain_22;
  	  		this.gain[1][2]=gain_23;
  	  		this.gain[2][0]=gain_31;
  	  		this.gain[2][1]=gain_32;
  	  		this.gain[2][2]=gain_33;
  	  		this.balanceRed[0][0]=balanceRed_11;
  	  		this.balanceRed[0][1]=balanceRed_12;
  	  		this.balanceRed[0][2]=balanceRed_13;
  	  		this.balanceRed[1][0]=balanceRed_21;
  	  		this.balanceRed[1][1]=balanceRed_22;
  	  		this.balanceRed[1][2]=balanceRed_23;
  	  		this.balanceRed[2][0]=balanceRed_31;
  	  		this.balanceRed[2][1]=balanceRed_32;
  	  		this.balanceRed[2][2]=balanceRed_33;
  	  		this.balanceBlue[0][0]=balanceBlue_11;
  	  		this.balanceBlue[0][1]=balanceBlue_12;
  	  		this.balanceBlue[0][2]=balanceBlue_13;
  	  		this.balanceBlue[1][0]=balanceBlue_21;
  	  		this.balanceBlue[1][1]=balanceBlue_22;
  	  		this.balanceBlue[1][2]=balanceBlue_23;
  	  		this.balanceBlue[2][0]=balanceBlue_31;
  	  		this.balanceBlue[2][1]=balanceBlue_32;
  	  		this.balanceBlue[2][2]=balanceBlue_33;
  	    }
  		public void setProperties(String prefix,Properties properties){
  			int i,j;
  			for (i=0;i<this.gain.length;i++) for (j=0;j<this.gain[i].length;j++)
  			  properties.setProperty(prefix+"gain_"+i+"_"+j,this.gain[i][j]+"");
  			for (i=0;i<this.balanceRed.length;i++) for (j=0;j<this.balanceRed[i].length;j++)
  				  properties.setProperty(prefix+"balanceRed_"+i+"_"+j,this.balanceRed[i][j]+"");
  			for (i=0;i<this.balanceBlue.length;i++) for (j=0;j<this.balanceBlue[i].length;j++)
  				  properties.setProperty(prefix+"balanceBlue_"+i+"_"+j,this.balanceBlue[i][j]+"");
  		}
  		public void getProperties(String prefix,Properties properties){
//  			this.oversample=Integer.parseInt(properties.getProperty(prefix+"oversample"));
  			int i,j;
  			String s;
  			for (i=0;i<this.gain.length;i++) for (j=0;j<this.gain[i].length;j++) {
  				s=properties.getProperty(prefix+"gain_"+i+"_"+j);
  				if (s!=null) this.gain[i][j]=Double.parseDouble(s);
  			}
  			for (i=0;i<this.balanceRed.length;i++) for (j=0;j<this.balanceRed[i].length;j++) {
  				s=properties.getProperty(prefix+"balanceRed_"+i+"_"+j);
  				if (s!=null) this.balanceRed[i][j]=Double.parseDouble(s);
  			}
  			for (i=0;i<this.balanceBlue.length;i++) for (j=0;j<this.balanceBlue[i].length;j++) {
  				s=properties.getProperty(prefix+"balanceBlue_"+i+"_"+j);
  				if (s!=null) this.balanceBlue[i][j]=Double.parseDouble(s);
  			}
  		}

    }
    /* ======================================================================== */
    public static class NonlinParameters {
    	public boolean useRejectBlocksFilter;
    	public boolean combineBothModes;
   	public int    maskFFTSize; // 256
   	public int    blockPeriod; // 32
    	public double rejectFreqSigma; // 1.0, frequency domain
    	public double lowPassSigma;    // 5.0, spatial domain
    	public double filtMin;
    	public double filtMax;
  	public double [][] thresholdCorrection=new double[3][3]; // apply to filtMin and filtMax
  	public double [] thresholdCorr={
  			1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
  			1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
  			1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
  			1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0
  	};
    	public double threshold;
  	public boolean useDiffNoiseGains;
  	public double [] noiseGainWeights=new double[3];
  	public double blurSigma;     // blur sigma for mask calculation (blur convolution kernels for noise gain calculation
  	public double  noiseGainPower;
    	public boolean showMask;
  // ring filter
    	public boolean useRingFilter;    // filter out spots on denoise mask
      public double minMaxValue;       // minimal value (relative to filtMax) of the local maximum to be processed
      public double overRingThreshold; // ratio of local max. and maximal value in the surrounding ring to trigger filter
      public double overRingLimit;     // limit values in the center circle to scaled maximum in a ring
      public double ringIR;            // ring inner radius (center circle radius)
      public double ringOR;            // ring outer radius

    	public NonlinParameters(
    			boolean useRejectBlocksFilter,
    			boolean combineBothModes,
    			int maskFFTSize,
    			int blockPeriod,
    			double rejectFreqSigma,
    			double lowPassSigma,
    			double filtMin,
    			double filtMax,
    			double thresholdCorrection_11,
    			double thresholdCorrection_12,
    			double thresholdCorrection_13,
    			double thresholdCorrection_21,
    			double thresholdCorrection_22,
    			double thresholdCorrection_23,
    			double thresholdCorrection_31,
    			double thresholdCorrection_32,
    			double thresholdCorrection_33,
    			double threshold,
  			boolean useDiffNoiseGains,
  			double noiseGainWeights_0, // r
  			double noiseGainWeights_1, // b
  			double noiseGainWeights_2, // g
  			double blurSigma,     // blur sigma for mask calculation (blur convolution kernels for noise gain calculation
  			double noiseGainPower,
  			boolean useRingFilter,    // filter out spots on denoise mask
    		    double minMaxValue,       // minimal value (relative to filtMax) of the local maximum to be processed
    		    double overRingThreshold, // ratio of local max. and maximal value in the surrounding ring to trigger filter
    		    double overRingLimit,     // limit values in the center circle to scaled maximum in a ring
    		    double ringIR,            // ring inner radius (center circle radius)
    		    double ringOR             // ring outer radius

  ) {
    		this.useRejectBlocksFilter=useRejectBlocksFilter;
    		this.combineBothModes=combineBothModes;
    		this.maskFFTSize = maskFFTSize;
    		this.blockPeriod = blockPeriod;
    		this.rejectFreqSigma = rejectFreqSigma;
    		this.lowPassSigma = lowPassSigma;
    		this.filtMin = filtMin;
    		this.filtMax = filtMax;
    		this.thresholdCorrection[0][0]=thresholdCorrection_11;
    		this.thresholdCorrection[0][1]=thresholdCorrection_12;
    		this.thresholdCorrection[0][2]=thresholdCorrection_13;
    		this.thresholdCorrection[1][0]=thresholdCorrection_21;
    		this.thresholdCorrection[1][1]=thresholdCorrection_22;
    		this.thresholdCorrection[1][2]=thresholdCorrection_23;
    		this.thresholdCorrection[2][0]=thresholdCorrection_31;
    		this.thresholdCorrection[2][1]=thresholdCorrection_32;
    		this.thresholdCorrection[2][2]=thresholdCorrection_33;
    		this.threshold = threshold;
  		this.useDiffNoiseGains=useDiffNoiseGains;
  		this.noiseGainWeights[0]=noiseGainWeights_0;
  		this.noiseGainWeights[1]=noiseGainWeights_1;
  		this.noiseGainWeights[2]=noiseGainWeights_2;
  		this.blurSigma=blurSigma;
  		this.noiseGainPower=noiseGainPower;
  		this.useRingFilter=useRingFilter;
  		this.minMaxValue=minMaxValue;
  		this.overRingThreshold=overRingThreshold;
  		this.overRingLimit=overRingLimit;
  		this.ringIR=ringIR;
  		this.ringOR=ringOR;

    	}
    	public void modifyNumChannels(int numChannels){
    		if ((numChannels>0) && (numChannels!=this.thresholdCorr.length)){
    			double [] thresholdCorr1=this.thresholdCorr;
    			this.thresholdCorr=  new double[numChannels];
    			for (int i=0;i<numChannels;i++) {
    				int j=i;
    				if (j>=thresholdCorr1.length) j=thresholdCorr1.length-1;
    				this.thresholdCorr[i]=thresholdCorr1[j];
    			}
    		}
    	}

  	public void setProperties(String prefix,Properties properties){
//  		properties.setProperty(prefix+"oversample",this.oversample+"");
  		properties.setProperty(prefix+"useRejectBlocksFilter",this.useRejectBlocksFilter+"");
  		properties.setProperty(prefix+"combineBothModes",this.combineBothModes+"");
  		properties.setProperty(prefix+"maskFFTSize",this.maskFFTSize+"");
  		properties.setProperty(prefix+"blockPeriod",this.blockPeriod+"");
  		properties.setProperty(prefix+"rejectFreqSigma",this.rejectFreqSigma+"");
  		properties.setProperty(prefix+"lowPassSigma",this.lowPassSigma+"");
  		properties.setProperty(prefix+"filtMin",this.filtMin+"");
  		properties.setProperty(prefix+"filtMax",this.filtMax+"");
  		for (int i=0;i<this.thresholdCorrection.length;i++) for (int j=0;j<this.thresholdCorrection[i].length;j++)
  		  properties.setProperty(prefix+"thresholdCorrection_"+i+"_"+j,this.thresholdCorrection[i][j]+"");


  		properties.setProperty(prefix+"threshold",this.threshold+"");
  		properties.setProperty(prefix+"useDiffNoiseGains",this.useDiffNoiseGains+"");
  		for (int i=0;i<this.noiseGainWeights.length;i++)
  		   properties.setProperty(prefix+"noiseGainWeights_"+i,this.noiseGainWeights[i]+"");
  		properties.setProperty(prefix+"blurSigma",this.blurSigma+"");
  		properties.setProperty(prefix+"noiseGainPower",this.noiseGainPower+"");
  		properties.setProperty(prefix+"useRingFilter",this.useRingFilter+"");
  		properties.setProperty(prefix+"minMaxValue",this.minMaxValue+"");
  		properties.setProperty(prefix+"overRingThreshold",this.overRingThreshold+"");
  		properties.setProperty(prefix+"overRingLimit",this.overRingLimit+"");
  		properties.setProperty(prefix+"ringIR",this.ringIR+"");
  		properties.setProperty(prefix+"ringOR",this.ringOR+"");
  		properties.setProperty(prefix+"thresholdCorr",this.thresholdCorr.length+"");
  		for (int i =0;i<this.thresholdCorr.length;i++) properties.setProperty(prefix+"thresholdCorr_"+i,this.thresholdCorr[i]+"");
  	}
  	public void getProperties(String prefix,Properties properties){
  		//  		this.oversample=Integer.parseInt(properties.getProperty(prefix+"oversample"));
  		String s;
  		this.useRejectBlocksFilter=Boolean.parseBoolean(properties.getProperty(prefix+"useRejectBlocksFilter"));
  		this.combineBothModes=Boolean.parseBoolean(properties.getProperty(prefix+"combineBothModes"));
  		this.maskFFTSize=Integer.parseInt(properties.getProperty(prefix+"maskFFTSize"));
  		this.blockPeriod=Integer.parseInt(properties.getProperty(prefix+"blockPeriod"));
  		this.rejectFreqSigma=Double.parseDouble(properties.getProperty(prefix+"rejectFreqSigma"));
  		this.lowPassSigma=Double.parseDouble(properties.getProperty(prefix+"lowPassSigma"));
  		this.filtMin=Double.parseDouble(properties.getProperty(prefix+"filtMin"));
  		this.filtMax=Double.parseDouble(properties.getProperty(prefix+"filtMax"));
  		for (int i=0;i<this.thresholdCorrection.length;i++) for (int j=0;j<this.thresholdCorrection[i].length;j++)
  			this.thresholdCorrection[i][j]=Double.parseDouble(properties.getProperty(prefix+"thresholdCorrection_"+i+"_"+j));
  		this.threshold=Double.parseDouble(properties.getProperty(prefix+"threshold"));
  		this.useDiffNoiseGains=Boolean.parseBoolean(properties.getProperty(prefix+"useDiffNoiseGains"));
  		for (int i=0;i<this.noiseGainWeights.length;i++)
  			this.noiseGainWeights[i]=Double.parseDouble(properties.getProperty(prefix+"noiseGainWeights_"+i));
  		this.blurSigma=Double.parseDouble(properties.getProperty(prefix+"blurSigma"));
  		this.noiseGainPower=Double.parseDouble(properties.getProperty(prefix+"noiseGainPower"));
  		s=properties.getProperty(prefix+"useRingFilter");
  		if ((s==null) || (s=="")) return; // earlier revision
  		this.useRingFilter=Boolean.parseBoolean(properties.getProperty(prefix+"useRingFilter"));
  		this.minMaxValue=Double.parseDouble(properties.getProperty(prefix+"minMaxValue"));
  		this.overRingThreshold=Double.parseDouble(properties.getProperty(prefix+"overRingThreshold"));
  		this.overRingLimit=Double.parseDouble(properties.getProperty(prefix+"overRingLimit"));
  		this.ringIR=Double.parseDouble(properties.getProperty(prefix+"ringIR"));
  		this.ringOR=Double.parseDouble(properties.getProperty(prefix+"ringOR"));
  		if (properties.getProperty(prefix+"thresholdCorr")!=null){
  			this.thresholdCorr=new double [Integer.parseInt(properties.getProperty(prefix+"thresholdCorr"))];
  			for (int i=0;i<this.thresholdCorr.length;i++) this.thresholdCorr[i]=Double.parseDouble(properties.getProperty(prefix+"thresholdCorr_"+i));
  		}
  	}

    }
  /* ======================================================================== */
    public static class SplitParameters {
  		public int oversample;
  		public int addLeft;
  		public int addTop;
  		public int addRight;
  		public int addBottom;

  		public SplitParameters(int oversample, int addLeft, int addTop,
  				int addRight, int addBottom) {
  			this.oversample = oversample;
  			this.addLeft = addLeft;
  			this.addTop = addTop;
  			this.addRight = addRight;
  			this.addBottom = addBottom;
  		}
  		public void setProperties(String prefix,Properties properties){
  			properties.setProperty(prefix+"oversample",this.oversample+"");
  			properties.setProperty(prefix+"addLeft",   this.addLeft+"");
  			properties.setProperty(prefix+"addTop",    this.addTop+"");
  			properties.setProperty(prefix+"addRight",  this.addRight+"");
  			properties.setProperty(prefix+"addBottom", this.addBottom+"");
  		}
  		public void getProperties(String prefix,Properties properties){
  			this.oversample=Integer.parseInt(properties.getProperty(prefix+"oversample"));
  			this.addLeft=Integer.parseInt(properties.getProperty(prefix+"addLeft"));
  			this.addTop=Integer.parseInt(properties.getProperty(prefix+"addTop"));
  			this.addRight=Integer.parseInt(properties.getProperty(prefix+"addRight"));
  			this.addBottom=Integer.parseInt(properties.getProperty(prefix+"addBottom"));
  		}
  	}
    public static class DCTParameters {
  		public int dct_size =            8; //
  		public int asym_size =          15; //
  		public int asym_pixels =         4; // maximal number of non-zero pixels in direct convolution kernel
  		public int asym_distance =       1; // how far to try a new asym kernel pixel from existing ones
  		public int dct_window =          1; // currently only 3 types of windows - 0 (none), 1 and 2
  		public int LMA_steps =         100;
  		public double fact_precision=    0.003; // stop iterations if error rms less than this part of target kernel rms
  		public double compactness =      0.0;
  		public double sym_compactness =  0.5;
  		public double dc_weight =        1.0;  // importance of dc realtive to rms_pure
  		public int    asym_tax_free  =   0;    // "compactness" does not apply to pixels with |x|<=asym_tax_free  and |y| <= asym_tax_free
  		public int    seed_size =        1;    // number of initial cells in asym_kernel - should be 4*b + 1 (X around center cell) or 4*n + 0  (X around between cells)
  		public double asym_random  =    -1; // initialize asym_kernel with random numbers
  		public double dbg_x =            2.7;
  		public double dbg_y =            0.0;
  		public double dbg_x1 =          -1.3;
  		public double dbg_y1 =           2.0;
  		public double dbg_sigma =        0.8;
  		public double dbg_src_size =     8.0; // trying to slightly scale in dct space. == dct = 1:1, dct+1.0 - shrink dct(dct+1.0)
  		public double dbg_scale =        1.0; // Should ==DCT_PARAMETERS.dct_size / DCT_PARAMETERS.dbg_src_size
  		public double dbg_fold_scale =   1.0; // Modifies window during MDCT->DCT-IV folding
  		public String dbg_mask = ".........:::::::::.........:::::::::......*..:::::*:::.........:::::::::.........";
  		public int dbg_mode =            1; // 0 - old LMA, 1 - new LMA - *** not used anymore ***
  		public int dbg_window_mode =     1; // 0 - none, 1 - square, 2 - sin 3 - sin^2 Now _should_ be square !!!
  		public boolean centerWindowToTarget = true;
  		// parameters to extract a kernel from the kernel image file
  		public int    color_channel =    2; // green (<0 - use simulated kernel, also will use simulated if kernels are not set)
  		public int    decimation =       2; // decimate original kernel this much in each direction
  		public double decimateSigma =   -1.0; // special mode for 2:1 deciamtion
  		public int    tileX =            82;  // number of kernel tile (0..163)
  		public int    tileY =            62;  // number of kernel tile (0..122)
  		public int    kernel_chn =      -1; // camera channel calibration to use for aberration correction ( < 0 - no correction)
  		public boolean normalize =       true; // normalize both sym and asym kernels (asym to have sum==1, sym to have sum = dct_size
  		public boolean normalize_sym =   true; // normalize sym kernels separately
  		public boolean antiwindow =      false; // divide symmetrical kernel by a window function
  		public boolean skip_sym =        false; // do not apply symmetrical correction
  		public boolean convolve_direct = false; // do not apply symmetrical correction

  		// colors should be balanced before DCT color conversion!
  		public double novignetting_r    = 0.2644; // reg gain in the center of sensor calibration R (instead of vignetting)
  		public double novignetting_g    = 0.3733; // green gain in the center of sensor calibration G
  		public double novignetting_b    = 0.2034; // blue gain in the center of sensor calibration B

  		public double scale_r =           1.0; // extra gain correction after vignetting or nonvignetting, before other processing
  		public double scale_g =           1.0;
  		public double scale_b =           1.0;

  		public double vignetting_max    = 0.4; // value in vignetting data to correspond to 1x in the kernel
  		public double vignetting_range  = 5.0; // do not try to correct vignetting less than vignetting_max/vignetting_range

  		public boolean post_debayer     = false; // perform de-bayer after aberrations in pixel domain
  		public boolean color_DCT        = true; // false - use old color processing mode
  		public double  sigma_rb =         0.9; // additional (to G) blur for R and B
  		public double  sigma_y =          0.7; // blur for G contribution to Y
  		public double  sigma_color =      0.7; // blur for Pb and Pr in addition to that of Y
  		public double  line_thershold =   1.0; // line detection amplitude to apply line enhancement - not used?
  		public boolean nonlin =           true; // enable nonlinear processing (including denoise
  		public double  nonlin_max_y =     1.0; // maximal amount of nonlinear line/edge emphasis for Y component
  		public double  nonlin_max_c =     1.0; // maximal amount of nonlinear line/edge emphasis for C component
  		public double  nonlin_y =         0.1; // amount of nonlinear line/edge emphasis for Y component
  		public double  nonlin_c =         0.01; // amount of nonlinear line/edge emphasis for C component
  		public double  nonlin_corn =      0.5;  // relative weight for nonlinear corner elements
  		public boolean denoise =          true; // suppress noise during nonlinear processing
  		public double  denoise_y =        1.0;  // maximal total smoothing of the Y post-kernel (will compete with edge emphasis)
  		public double  denoise_c =        1.0;  // maximal total smoothing of the color differences post-kernel (will compete with edge emphasis)
  		public double  denoise_y_corn =   0.3;  // weight of the 4 corner pixels during denoise y (straight - 1-denoise_y_corn)
  		public double  denoise_c_corn =   0.3;  // weight of the 4 corner pixels during denoise y (straight - 1-denoise_c_corn)


  		public DCTParameters(){}

  		public DCTParameters(
  				int dct_size,
  				int asym_size,
  				int asym_pixels,
  				int asym_distance,
  				int dct_window,
  				double compactness,
  				int asym_tax_free,
  				int seed_size) {
  			this.dct_size =       dct_size;
  			this.asym_size =      asym_size;
  			this.asym_pixels =    asym_pixels;
  			this.asym_distance =  asym_distance;
  			this.dct_window =     dct_window;
  			this.compactness =    compactness;
  			this.asym_tax_free =  asym_tax_free;
  			this.seed_size =      seed_size;
  		}
  		public void setProperties(String prefix,Properties properties){
  			properties.setProperty(prefix+"dct_size",this.dct_size+"");
  			properties.setProperty(prefix+"asym_size",this.asym_size+"");
  			properties.setProperty(prefix+"asym_pixels",this.asym_pixels+"");
  			properties.setProperty(prefix+"asym_distance",this.asym_distance+"");
  			properties.setProperty(prefix+"dct_window",   this.dct_window+"");
  			properties.setProperty(prefix+"compactness",  this.compactness+"");
  			properties.setProperty(prefix+"sym_compactness",  this.sym_compactness+"");
  			properties.setProperty(prefix+"dc_weight",  this.dc_weight+"");
  			properties.setProperty(prefix+"fact_precision",  this.fact_precision+"");
  			properties.setProperty(prefix+"asym_tax_free",  this.asym_tax_free+"");
  			properties.setProperty(prefix+"seed_size",  this.seed_size+"");
  			properties.setProperty(prefix+"asym_random",  this.asym_random+"");
  			properties.setProperty(prefix+"LMA_steps",  this.LMA_steps+"");
  			properties.setProperty(prefix+"dbg_x",      this.dbg_x+"");
  			properties.setProperty(prefix+"dbg_y",      this.dbg_y+"");
  			properties.setProperty(prefix+"dbg_x1",     this.dbg_x1+"");
  			properties.setProperty(prefix+"dbg_y1",     this.dbg_y1+"");
  			properties.setProperty(prefix+"dbg_sigma",  this.dbg_sigma+"");
  			properties.setProperty(prefix+"dbg_src_size",this.dbg_src_size+"");
  			properties.setProperty(prefix+"dbg_scale",  this.dbg_scale+"");
  			properties.setProperty(prefix+"dbg_fold_scale",  this.dbg_fold_scale+"");
  			properties.setProperty(prefix+"dbg_mask",   this.dbg_mask+"");
  			properties.setProperty(prefix+"dbg_mode",   this.dbg_mode+"");
  			properties.setProperty(prefix+"dbg_window_mode",    this.dbg_window_mode+"");
  			properties.setProperty(prefix+"centerWindowToTarget",   this.centerWindowToTarget+"");
  			properties.setProperty(prefix+"color_channel",      this.color_channel+"");
  			properties.setProperty(prefix+"decimation",         this.decimation+"");
  			properties.setProperty(prefix+"decimateSigma",      this.decimateSigma+"");
  			properties.setProperty(prefix+"tileX",              this.tileX+"");
  			properties.setProperty(prefix+"tileY",              this.tileY+"");
  			properties.setProperty(prefix+"kernel_chn",         this.kernel_chn+"");
  			properties.setProperty(prefix+"normalize",          this.normalize+"");
  			properties.setProperty(prefix+"normalize_sym",      this.normalize_sym+"");
  			properties.setProperty(prefix+"antiwindow",         this.antiwindow+"");
  			properties.setProperty(prefix+"skip_sym",           this.skip_sym+"");
  			properties.setProperty(prefix+"convolve_direct",    this.convolve_direct+"");
  			properties.setProperty(prefix+"novignetting_r",     this.novignetting_r+"");
  			properties.setProperty(prefix+"novignetting_g",     this.novignetting_g+"");
  			properties.setProperty(prefix+"novignetting_b",     this.novignetting_b+"");
  			properties.setProperty(prefix+"scale_r",            this.scale_r+"");
  			properties.setProperty(prefix+"scale_g",            this.scale_g+"");
  			properties.setProperty(prefix+"scale_b",            this.scale_b+"");
  			properties.setProperty(prefix+"vignetting_max",     this.vignetting_max+"");
  			properties.setProperty(prefix+"vignetting_range",   this.vignetting_range+"");
  			properties.setProperty(prefix+"post_debayer",       this.post_debayer+"");
  			properties.setProperty(prefix+"color_DCT",          this.color_DCT+"");
  			properties.setProperty(prefix+"sigma_rb",           this.sigma_rb+"");
  			properties.setProperty(prefix+"sigma_y",            this.sigma_y+"");
  			properties.setProperty(prefix+"sigma_color",        this.sigma_color+"");
  			properties.setProperty(prefix+"line_thershold",     this.line_thershold+"");
  			properties.setProperty(prefix+"nonlin",             this.nonlin+"");
  			properties.setProperty(prefix+"nonlin_max_y",       this.nonlin_max_y+"");
  			properties.setProperty(prefix+"nonlin_max_c",       this.nonlin_max_c+"");
  			properties.setProperty(prefix+"nonlin_y",           this.nonlin_y+"");
  			properties.setProperty(prefix+"nonlin_c",           this.nonlin_c+"");
  			properties.setProperty(prefix+"nonlin_corn",        this.nonlin_corn+"");
  			properties.setProperty(prefix+"denoise",            this.denoise+"");
  			properties.setProperty(prefix+"denoise_y",          this.denoise_y+"");
  			properties.setProperty(prefix+"denoise_c",          this.denoise_c+"");
  			properties.setProperty(prefix+"denoise_y_corn",     this.denoise_y_corn+"");
  			properties.setProperty(prefix+"denoise_c_corn",     this.denoise_c_corn+"");



  		}
  		public void getProperties(String prefix,Properties properties){
  			if (properties.getProperty(prefix+"dct_size")!=null) this.dct_size=Integer.parseInt(properties.getProperty(prefix+"dct_size"));
  			if (properties.getProperty(prefix+"asym_size")!=null) this.asym_size=Integer.parseInt(properties.getProperty(prefix+"asym_size"));
  			if (properties.getProperty(prefix+"asym_pixels")!=null) this.asym_pixels=Integer.parseInt(properties.getProperty(prefix+"asym_pixels"));
  			if (properties.getProperty(prefix+"asym_distance")!=null) this.asym_distance=Integer.parseInt(properties.getProperty(prefix+"asym_distance"));
  			if (properties.getProperty(prefix+"dct_window")!=null) this.dct_window=Integer.parseInt(properties.getProperty(prefix+"dct_window"));
  			if (properties.getProperty(prefix+"compactness")!=null) this.compactness=Double.parseDouble(properties.getProperty(prefix+"compactness"));
  			if (properties.getProperty(prefix+"sym_compactness")!=null) this.sym_compactness=Double.parseDouble(properties.getProperty(prefix+"sym_compactness"));
  			if (properties.getProperty(prefix+"dc_weight")!=null) this.dc_weight=Double.parseDouble(properties.getProperty(prefix+"dc_weight"));
  			if (properties.getProperty(prefix+"fact_precision")!=null) this.fact_precision=Double.parseDouble(properties.getProperty(prefix+"fact_precision"));
  			if (properties.getProperty(prefix+"asym_tax_free")!=null) this.asym_tax_free=Integer.parseInt(properties.getProperty(prefix+"asym_tax_free"));
  			if (properties.getProperty(prefix+"seed_size")!=null) this.seed_size=Integer.parseInt(properties.getProperty(prefix+"seed_size"));
  			if (properties.getProperty(prefix+"asym_random")!=null) this.asym_random=Double.parseDouble(properties.getProperty(prefix+"asym_random"));
  			if (properties.getProperty(prefix+"LMA_steps")!=null) this.LMA_steps=Integer.parseInt(properties.getProperty(prefix+"LMA_steps"));
  			if (properties.getProperty(prefix+"dbg_x")!=null) this.dbg_x=Double.parseDouble(properties.getProperty(prefix+"dbg_x"));
  			if (properties.getProperty(prefix+"dbg_y")!=null) this.dbg_y=Double.parseDouble(properties.getProperty(prefix+"dbg_y"));
  			if (properties.getProperty(prefix+"dbg_x1")!=null) this.dbg_x1=Double.parseDouble(properties.getProperty(prefix+"dbg_x1"));
  			if (properties.getProperty(prefix+"dbg_y1")!=null) this.dbg_y1=Double.parseDouble(properties.getProperty(prefix+"dbg_y1"));
  			if (properties.getProperty(prefix+"dbg_sigma")!=null) this.dbg_sigma=Double.parseDouble(properties.getProperty(prefix+"dbg_sigma"));
  			if (properties.getProperty(prefix+"dbg_src_size")!=null) this.dbg_src_size=Double.parseDouble(properties.getProperty(prefix+"dbg_src_size"));
  			if (properties.getProperty(prefix+"dbg_scale")!=null) this.dbg_scale=Double.parseDouble(properties.getProperty(prefix+"dbg_scale"));
  			if (properties.getProperty(prefix+"dbg_fold_scale")!=null) this.dbg_fold_scale=Double.parseDouble(properties.getProperty(prefix+"dbg_fold_scale"));
  			if (properties.getProperty(prefix+"dbg_mask")!=null) this.dbg_mask=properties.getProperty(prefix+"dbg_mask");
  			if (properties.getProperty(prefix+"dbg_mode")!=null) this.dbg_mode=Integer.parseInt(properties.getProperty(prefix+"dbg_mode"));
  			if (properties.getProperty(prefix+"centerWindowToTarget")!=null) this.centerWindowToTarget=Boolean.parseBoolean(properties.getProperty(prefix+"centerWindowToTarget"));
  			if (properties.getProperty(prefix+"color_channel")!=null) this.color_channel=Integer.parseInt(properties.getProperty(prefix+"color_channel"));
  			if (properties.getProperty(prefix+"decimation")!=null) this.decimation=Integer.parseInt(properties.getProperty(prefix+"decimation"));
  			if (properties.getProperty(prefix+"decimateSigma")!=null) this.decimateSigma=Double.parseDouble(properties.getProperty(prefix+"decimateSigma"));
  			if (properties.getProperty(prefix+"tileX")!=null) this.tileX=Integer.parseInt(properties.getProperty(prefix+"tileX"));
  			if (properties.getProperty(prefix+"tileY")!=null) this.tileY=Integer.parseInt(properties.getProperty(prefix+"tileY"));
  			if (properties.getProperty(prefix+"dbg_window_mode")!=null) this.dbg_window_mode=Integer.parseInt(properties.getProperty(prefix+"dbg_window_mode"));
  			if (properties.getProperty(prefix+"kernel_chn")!=null) this.kernel_chn=Integer.parseInt(properties.getProperty(prefix+"kernel_chn"));
  			if (properties.getProperty(prefix+"normalize")!=null) this.normalize=Boolean.parseBoolean(properties.getProperty(prefix+"normalize"));
  			if (properties.getProperty(prefix+"normalize_sym")!=null) this.normalize_sym=Boolean.parseBoolean(properties.getProperty(prefix+"normalize_sym"));
  			if (properties.getProperty(prefix+"antiwindow")!=null) this.antiwindow=Boolean.parseBoolean(properties.getProperty(prefix+"antiwindow"));
  			if (properties.getProperty(prefix+"skip_sym")!=null) this.skip_sym=Boolean.parseBoolean(properties.getProperty(prefix+"skip_sym"));
  			if (properties.getProperty(prefix+"convolve_direct")!=null) this.convolve_direct=Boolean.parseBoolean(properties.getProperty(prefix+"convolve_direct"));
  			if (properties.getProperty(prefix+"novignetting_r")!=null) this.novignetting_r=Double.parseDouble(properties.getProperty(prefix+"novignetting_r"));
  			if (properties.getProperty(prefix+"novignetting_g")!=null) this.novignetting_g=Double.parseDouble(properties.getProperty(prefix+"novignetting_g"));
  			if (properties.getProperty(prefix+"novignetting_b")!=null) this.novignetting_b=Double.parseDouble(properties.getProperty(prefix+"novignetting_b"));
  			if (properties.getProperty(prefix+"scale_r")!=null)        this.scale_r=Double.parseDouble(properties.getProperty(prefix+"scale_r"));
  			if (properties.getProperty(prefix+"scale_g")!=null)        this.scale_g=Double.parseDouble(properties.getProperty(prefix+"scale_g"));
  			if (properties.getProperty(prefix+"scale_b")!=null)        this.scale_b=Double.parseDouble(properties.getProperty(prefix+"scale_b"));
  			if (properties.getProperty(prefix+"vignetting_max")!=null) this.vignetting_max=Double.parseDouble(properties.getProperty(prefix+"vignetting_max"));
  			if (properties.getProperty(prefix+"vignetting_range")!=null) this.vignetting_range=Double.parseDouble(properties.getProperty(prefix+"vignetting_range"));
  			if (properties.getProperty(prefix+"post_debayer")!=null)   this.post_debayer=Boolean.parseBoolean(properties.getProperty(prefix+"post_debayer"));
  			if (properties.getProperty(prefix+"color_DCT")!=null)      this.color_DCT=Boolean.parseBoolean(properties.getProperty(prefix+"color_DCT"));
  			if (properties.getProperty(prefix+"sigma_rb")!=null)       this.sigma_rb=Double.parseDouble(properties.getProperty(prefix+"sigma_rb"));
  			if (properties.getProperty(prefix+"sigma_y")!=null)        this.sigma_y=Double.parseDouble(properties.getProperty(prefix+"sigma_y"));
  			if (properties.getProperty(prefix+"sigma_color")!=null)    this.sigma_color=Double.parseDouble(properties.getProperty(prefix+"sigma_color"));
  			if (properties.getProperty(prefix+"line_thershold")!=null) this.line_thershold=Double.parseDouble(properties.getProperty(prefix+"line_thershold"));
  			if (properties.getProperty(prefix+"nonlin")!=null)         this.nonlin=Boolean.parseBoolean(properties.getProperty(prefix+"nonlin"));
  			if (properties.getProperty(prefix+"nonlin_max_y")!=null)   this.nonlin_max_y=Double.parseDouble(properties.getProperty(prefix+"nonlin_max_y"));
  			if (properties.getProperty(prefix+"nonlin_max_c")!=null)   this.nonlin_max_c=Double.parseDouble(properties.getProperty(prefix+"nonlin_max_c"));
  			if (properties.getProperty(prefix+"nonlin_y")!=null)       this.nonlin_y=Double.parseDouble(properties.getProperty(prefix+"nonlin_y"));
  			if (properties.getProperty(prefix+"nonlin_c")!=null)       this.nonlin_c=Double.parseDouble(properties.getProperty(prefix+"nonlin_c"));
  			if (properties.getProperty(prefix+"nonlin_corn")!=null)    this.nonlin_corn=Double.parseDouble(properties.getProperty(prefix+"nonlin_corn"));
  			if (properties.getProperty(prefix+"denoise")!=null)        this.denoise=Boolean.parseBoolean(properties.getProperty(prefix+"denoise"));
  			if (properties.getProperty(prefix+"denoise_y")!=null)      this.denoise_y=Double.parseDouble(properties.getProperty(prefix+"denoise_y"));
  			if (properties.getProperty(prefix+"denoise_c")!=null)      this.denoise_c=Double.parseDouble(properties.getProperty(prefix+"denoise_c"));
  			if (properties.getProperty(prefix+"denoise_y_corn")!=null) this.denoise_y_corn=Double.parseDouble(properties.getProperty(prefix+"denoise_y_corn"));
  			if (properties.getProperty(prefix+"denoise_c_corn")!=null) this.denoise_c_corn=Double.parseDouble(properties.getProperty(prefix+"denoise_c_corn"));

  		}
  		public boolean showDialog() {
  			GenericDialog gd = new GenericDialog("Set DCT parameters");
  			gd.addNumericField("DCT size",                                                       this.dct_size,            0);
  			gd.addNumericField("Size of asymmetrical (non-DCT) kernel",                          this.asym_size,           0);
  			gd.addNumericField("Maximal number of non-zero pixels in direct convolution kernel", this.asym_pixels,         0);
  			gd.addNumericField("How far to try a new asym kernel pixel from existing ones",      this.asym_distance,       0);
  			gd.addNumericField("MDCT window type (0,1,2)",                                       this.dct_window,          0);
  			gd.addNumericField("LMA_steps",                                                      this.LMA_steps,           0);
  			gd.addNumericField("Compactness (punish off-center asym_kernel pixels (proportional to r^2)", this.compactness,2);
  			gd.addNumericField("Symmetrical kernel compactness (proportional to r^2)",           this.sym_compactness,     2);
  			gd.addNumericField("Relative importance of DC error to RMS",                         this.dc_weight,           2);
  			gd.addNumericField("Factorization target precision (stop if achieved)",              this.fact_precision,      4);
  			gd.addNumericField("Do not punish pixels in the square around center",               this.asym_tax_free,       0);
  			gd.addNumericField("Start asym_kernel with this number of pixels (0 - single, 4n+0 (X between cells), 4*n+1 - x around center cell",               this.seed_size,     0); //0..2
  			gd.addNumericField("Initialize asym_kernel with random numbers (amplitude)",         this.asym_random,         2);
  			gd.addNumericField("dbg_x",                                                          this.dbg_x,               2);
  			gd.addNumericField("dbg_y",                                                          this.dbg_y,               2);
  			gd.addNumericField("dbg_x1",                                                         this.dbg_x1,              2);
  			gd.addNumericField("dbg_y1",                                                         this.dbg_y1,              2);
  			gd.addNumericField("dbg_sigma",                                                      this.dbg_sigma,           3);
  			gd.addNumericField("== dct_size = 1:1, dct+1.0 - shrink dct(dct+1.0)",               this.dbg_src_size,        3);
  			gd.addNumericField("Should ==DCT_PARAMETERS.dct_size / DCT_PARAMETERS.dbg_src_size", this.dbg_scale,           3);
  			gd.addNumericField("Modifies window during MDCT->DCT-IV folding",                    this.dbg_fold_scale,      3);
  			gd.addStringField ("Debug mask (anything but * is false)",                           this.dbg_mask,          100);
  			gd.addNumericField("LMA implementation: 0 - old, 1 - new",                           this.dbg_mode,            0);
  			gd.addNumericField("Convolution window: 0 - none, [1 - square], 2 - sin, 3 - sin^2", this.dbg_window_mode,     0);
  			gd.addCheckbox    ("Center convolution window around target kernel center",          this.centerWindowToTarget);
  			gd.addNumericField("Color channel to extract kernel (<0 - use synthetic)",           this.color_channel,       0);
  			gd.addNumericField("Convolution kernel decimation (original is normally 2x)",        this.decimation,          0);
  			gd.addNumericField("Smooth convolution kernel before decimation",                    this.decimateSigma,       3);
  			gd.addNumericField("Tile X to extract (0..163)",                                     this.tileX,               0);
  			gd.addNumericField("Tile Y to extract (0..122)",                                     this.tileY,               0);
  			gd.addNumericField("Calibration channel to use for aberration ( <0 - no correction)",this.kernel_chn,          0);
  			gd.addCheckbox    ("Normalize both sym and asym kernels ",                           this.normalize);
  			gd.addCheckbox    ("Normalize sym kernels separately",                               this.normalize_sym);
  			gd.addCheckbox    ("Divide symmetrical kernel by a window function",                 this.antiwindow);
  			gd.addCheckbox    ("Do not apply symmetrical (DCT) correction ",                     this.skip_sym);
  			gd.addCheckbox    ("Convolve directly with symmetrical kernel (debug feature) ",     this.convolve_direct);
  			gd.addNumericField("Reg gain in the center of sensor calibration R (instead of vignetting)",this.novignetting_r,   4);
  			gd.addNumericField("Green gain in the center of sensor calibration G (instead of vignetting)",this.novignetting_g, 4);
  			gd.addNumericField("Blue gain in the center of sensor calibration B (instead of vignetting)",this.novignetting_b,  4);
  			gd.addNumericField("Extra red correction to compensate for light temperature",               this.scale_r,  4);
  			gd.addNumericField("Extra green correction to compensate for light temperature",             this.scale_g,  4);
  			gd.addNumericField("Extra blue correction to compensate for light temperature",              this.scale_b,  4);
  			gd.addNumericField("Value (max) in vignetting data to correspond to 1x in the kernel",    this.vignetting_max,      3);
  			gd.addNumericField("Do not try to correct vignetting smaller than this fraction of max",  this.vignetting_range,  3);
  			gd.addCheckbox    ("Perform de-bayer after aberrations in pixel domain",              this.post_debayer             );
  			gd.addCheckbox    ("Use DCT-based color conversion (false - just LPF RGB with dbg_sigma)",this.color_DCT             );
  			gd.addNumericField("Gaussian sigma to apply to R and B (in addition to G), pix",      this.sigma_rb,            3);
  			gd.addNumericField("Gaussian sigma to apply to Y in the MDCT domain, pix",            this.sigma_y,             3);
  			gd.addNumericField("Gaussian sigma to apply to Pr and Pb in the MDCT domain, pix",    this.sigma_color,         3);
  			gd.addNumericField("Threshold for line detection (not yet used)",                     this.line_thershold,      3);

  			gd.addCheckbox    ("Use non-linear line emphasis and denoise",                        this.nonlin             );
  			gd.addNumericField("Maximal amount of non-linear emphasis for linear edges for Y component",  this.nonlin_max_y,3);
  			gd.addNumericField("Maximal amount of non-linear emphasis for linear edges for color diffs.", this.nonlin_max_c,3);
  			gd.addNumericField("Sensitivity of non-linear emphasis for linear edges for Y component",  this.nonlin_y,       3);
  			gd.addNumericField("Sensitivity of non-linear emphasis for linear edges for color diffs.", this.nonlin_c,       3);
  			gd.addNumericField("Corretion for diagonal/corner emphasis elements",                 this.nonlin_corn,         3);
  			gd.addCheckbox    ("Suppress noise during nonlinear processing",                      this.denoise             );
  			gd.addNumericField("Maximal total smoothing of the Y post-kernel (will compete with edge emphasis)",  this.denoise_y,3);
  			gd.addNumericField("Maximal total smoothing of the color differences post-kernel (will compete with edge emphasis)",  this.denoise_c,3);
  			gd.addNumericField("Weight of the 4 corner pixels during denoising y (straight - 1.0-denoise_y_corn)",  this.denoise_y_corn,3);
  			gd.addNumericField("Weight of the 4 corner pixels during denoising color ((straight - 1.0-denoise_c_corn))",  this.denoise_c_corn,3);

  			WindowTools.addScrollBars(gd);
  			gd.showDialog();

  			if (gd.wasCanceled()) return false;
  			this.dct_size=        (int) gd.getNextNumber();
  			this.asym_size=       (int) gd.getNextNumber();
  			this.asym_pixels=     (int) gd.getNextNumber();
  			this.asym_distance=   (int) gd.getNextNumber();
  			this.dct_window=      (int) gd.getNextNumber();
  			this.LMA_steps =      (int) gd.getNextNumber();
  			this.compactness =          gd.getNextNumber();
  			this.sym_compactness =      gd.getNextNumber();
  			this.dc_weight =            gd.getNextNumber();
  			this.fact_precision =       gd.getNextNumber();
  			this.asym_tax_free =  (int) gd.getNextNumber();
  			this.seed_size =      (int) gd.getNextNumber();
  			this.asym_random=           gd.getNextNumber();
  			this.dbg_x=                 gd.getNextNumber();
  			this.dbg_y=                 gd.getNextNumber();
  			this.dbg_x1=                gd.getNextNumber();
  			this.dbg_y1=                gd.getNextNumber();
  			this.dbg_sigma=             gd.getNextNumber();
  			this.dbg_src_size=          gd.getNextNumber();
  			this.dbg_scale=             gd.getNextNumber();
  			this.dbg_fold_scale=        gd.getNextNumber();
  			this.dbg_mask=              gd.getNextString();
  			this.dbg_mode=        (int) gd.getNextNumber();
  			this.dbg_window_mode= (int) gd.getNextNumber();
  			this.centerWindowToTarget=  gd.getNextBoolean();
  			this.color_channel=   (int) gd.getNextNumber();
  			this.decimation=      (int) gd.getNextNumber();
  			this.decimateSigma=         gd.getNextNumber();
  			this.tileX=           (int) gd.getNextNumber();
  			this.tileY=           (int) gd.getNextNumber();
  			this.kernel_chn=      (int) gd.getNextNumber();
  			this.normalize=             gd.getNextBoolean();
  			this.normalize_sym=         gd.getNextBoolean();
  			this.antiwindow=            gd.getNextBoolean();
  			this.skip_sym=              gd.getNextBoolean();
  			this.convolve_direct=       gd.getNextBoolean();

  			this.novignetting_r=        gd.getNextNumber();
  			this.novignetting_g=        gd.getNextNumber();
  			this.novignetting_b=        gd.getNextNumber();
  			this.scale_r=               gd.getNextNumber();
  			this.scale_g=               gd.getNextNumber();
  			this.scale_b=               gd.getNextNumber();
  			this.vignetting_max=        gd.getNextNumber();
  			this.vignetting_range=      gd.getNextNumber();

  			this.post_debayer=          gd.getNextBoolean();
  			this.color_DCT=             gd.getNextBoolean();
  			this.sigma_rb=              gd.getNextNumber();
  			this.sigma_y=               gd.getNextNumber();
  			this.sigma_color=           gd.getNextNumber();
  			this.line_thershold=        gd.getNextNumber();

  			this.nonlin=             gd.getNextBoolean();
  			this.nonlin_max_y=          gd.getNextNumber();
  			this.nonlin_max_c=          gd.getNextNumber();
  			this.nonlin_y=              gd.getNextNumber();
  			this.nonlin_c=              gd.getNextNumber();
  			this.nonlin_corn=           gd.getNextNumber();

  			this.denoise=             gd.getNextBoolean();
  			this.denoise_y=          gd.getNextNumber();
  			this.denoise_c=          gd.getNextNumber();
  			this.denoise_y_corn=          gd.getNextNumber();
  			this.denoise_c_corn=          gd.getNextNumber();

  			//  	    MASTER_DEBUG_LEVEL= (int) gd.getNextNumber();
  			return true;
  		}
  	}


    /* ======================================================================== */
    public static class DebayerParameters {
  		public int size;
  		public double polarStep;
  		public double debayerThreshold;
  		public double debayerRelativeWidthGreen;
  		public double debayerRelativeWidthRedblue;
  		public double debayerRelativeWidthRedblueMain;
  		public double debayerRelativeWidthRedblueClones;
  		public double debayerGamma;
  		public double debayerBonus;
  		public double mainToAlias;
  		public double debayerMaskBlur;
  		public boolean debayerUseScissors;
  		public boolean debug;
  		public int xDebug;
  		public int yDebug;
  		public boolean debayerStacks;
  		public DebayerParameters(int size, double polarStep,
  				double debayerThreshold, double debayerRelativeWidthGreen,
  				double debayerRelativeWidthRedblue,
  				double debayerRelativeWidthRedblueMain,
  				double debayerRelativeWidthRedblueClones, double debayerGamma,
  				double debayerBonus, double mainToAlias, double debayerMaskBlur,
  				boolean debayerUseScissors,
  				boolean debug, int xDebug, int yDebug,
  				boolean debayerStacks) {
  			this.size = size;
  			this.polarStep = polarStep;
  			this.debayerThreshold = debayerThreshold;
  			this.debayerRelativeWidthGreen = debayerRelativeWidthGreen;
  			this.debayerRelativeWidthRedblue = debayerRelativeWidthRedblue;
  			this.debayerRelativeWidthRedblueMain = debayerRelativeWidthRedblueMain;
  			this.debayerRelativeWidthRedblueClones = debayerRelativeWidthRedblueClones;
  			this.debayerGamma = debayerGamma;
  			this.debayerBonus = debayerBonus;
  			this.mainToAlias = mainToAlias;
  			this.debayerMaskBlur = debayerMaskBlur;
  			this.debayerUseScissors = debayerUseScissors;
  			this.debug = debug;
  			this.xDebug = xDebug;
  			this.yDebug = yDebug;
  			this.debayerStacks = debayerStacks;
  		}
  		public void setProperties(String prefix,Properties properties){
//  			properties.setProperty(prefix+"oversample",this.oversample+"");
  			properties.setProperty(prefix+"size",this.size+"");
  			properties.setProperty(prefix+"polarStep",this.polarStep+"");
  			properties.setProperty(prefix+"debayerThreshold",this.debayerThreshold+"");
  			properties.setProperty(prefix+"debayerRelativeWidthGreen",this.debayerRelativeWidthGreen+"");
  			properties.setProperty(prefix+"debayerRelativeWidthRedblue",this.debayerRelativeWidthRedblue+"");
  			properties.setProperty(prefix+"debayerRelativeWidthRedblueMain",this.debayerRelativeWidthRedblueMain+"");
  			properties.setProperty(prefix+"debayerRelativeWidthRedblueClones",this.debayerRelativeWidthRedblueClones+"");
  			properties.setProperty(prefix+"debayerGamma",this.debayerGamma+"");
  			properties.setProperty(prefix+"debayerBonus",this.debayerBonus+"");
  			properties.setProperty(prefix+"mainToAlias",this.mainToAlias+"");
  			properties.setProperty(prefix+"debayerMaskBlur",this.debayerMaskBlur+"");
  			properties.setProperty(prefix+"debayerUseScissors",this.debayerUseScissors+"");
  			properties.setProperty(prefix+"debug",this.debug+"");
  			properties.setProperty(prefix+"xDebug",this.xDebug+"");
  			properties.setProperty(prefix+"yDebug",this.yDebug+"");
  			properties.setProperty(prefix+"debayerStacks",this.debayerStacks+"");
  		}
  		public void getProperties(String prefix,Properties properties){
  			this.size=                             Integer.parseInt(properties.getProperty(prefix+"size"));
  			this.polarStep=                        Double.parseDouble(properties.getProperty(prefix+"polarStep"));
  			this.debayerThreshold=                 Double.parseDouble(properties.getProperty(prefix+"debayerThreshold"));
  			this.debayerRelativeWidthGreen=        Double.parseDouble(properties.getProperty(prefix+"debayerRelativeWidthGreen"));
  			this.debayerRelativeWidthRedblue=      Double.parseDouble(properties.getProperty(prefix+"debayerRelativeWidthRedblue"));
  			this.debayerRelativeWidthRedblueMain=  Double.parseDouble(properties.getProperty(prefix+"debayerRelativeWidthRedblueMain"));
  			this.debayerRelativeWidthRedblueClones=Double.parseDouble(properties.getProperty(prefix+"debayerRelativeWidthRedblueClones"));
  			this.debayerGamma=                     Double.parseDouble(properties.getProperty(prefix+"debayerGamma"));
  			this.debayerBonus=                     Double.parseDouble(properties.getProperty(prefix+"debayerBonus"));
  			this.mainToAlias=                      Double.parseDouble(properties.getProperty(prefix+"mainToAlias"));
  			this.debayerMaskBlur=                  Double.parseDouble(properties.getProperty(prefix+"debayerMaskBlur"));
  			this.debayerUseScissors=               Boolean.parseBoolean(properties.getProperty(prefix+"debayerUseScissors"));
  			this.debug=                            Boolean.parseBoolean(properties.getProperty(prefix+"debug"));
  			this.xDebug=                           Integer.parseInt(properties.getProperty(prefix+"xDebug"));
  			this.yDebug=                           Integer.parseInt(properties.getProperty(prefix+"yDebug"));
  			this.debayerStacks=                    Boolean.parseBoolean(properties.getProperty(prefix+"debayerStacks"));
  		}
  	}

    public static class EquirectangularParameters {

		public double longitudeLeft=    -180.0;
		public double longitudeRight=    180.0;
		public double latitudeTop=        90.0;
		public double latitudeBottom=    -90.0;
		public int pixelsHorizontal=   14268;
		public int imageWidth=          2592;
		public int imageHeight=         1936;
		public double resolutionScale=     1.0;
		public double x0=                  0.0;
		public double y0=                  0.0;
		public int longitudeWidth=      3000; //pix
		public boolean clearFullMap=      true;
		public boolean clearAllMaps=      true;

		public boolean needRebuild=      false;
// common plane parameters (dual camera, triclope camera)
		public boolean generateCommonPlane =  false;
		public double projectionElevation=     0.0;
		public double projectionYaw=           0.0;
		public double projectionRoll=          0.0;

		public boolean matchPixelSize=    true; // disregard next value, calculate projectionPixelSize from the equirectangular map
		public double projectionPixelSize=0.00044036902;
		public int    projectionWidth=   2920;
		public int     projectionHeight= 2220;
		public double projectionCenterX= 0.5*this.projectionWidth;
		public double projectionCenterY= 0.5*this.projectionHeight;
		public double nominalHorizontalDisparity=60.0; // nominal distance between horizontal cameras, mm
		public boolean [] channelSelection=null;


    	public EquirectangularParameters(){
    	}
    	public boolean isNeedRebuildSet(){
    		boolean result=this.needRebuild;
    		this.needRebuild=false;
    		return result;
    	}

    	public boolean [] getChannelsToProcess(){ return this.channelSelection;}

    	public void setProperties(String prefix,Properties properties){
    		properties.setProperty(prefix+"longitudeLeft",this.longitudeLeft+"");
    		properties.setProperty(prefix+"longitudeRight",this.longitudeRight+"");
    		properties.setProperty(prefix+"latitudeTop",this.latitudeTop+"");
    		properties.setProperty(prefix+"latitudeBottom",this.latitudeBottom+"");

    		properties.setProperty(prefix+"pixelsHorizontal",this.pixelsHorizontal+"");
    		properties.setProperty(prefix+"imageWidth",this.imageWidth+"");
    		properties.setProperty(prefix+"imageHeight",this.imageHeight+"");

    		properties.setProperty(prefix+"resolutionScale",this.resolutionScale+"");
    		properties.setProperty(prefix+"x0",this.x0+"");
    		properties.setProperty(prefix+"y0",this.y0+"");

    		properties.setProperty(prefix+"longitudeWidth",this.longitudeWidth+"");

    		properties.setProperty(prefix+"clearFullMap",this.clearFullMap+"");
    		properties.setProperty(prefix+"clearAllMaps",this.clearAllMaps+"");
    		properties.setProperty(prefix+"generateCommonPlane",this.generateCommonPlane+"");
    		properties.setProperty(prefix+"projectionElevation",this.projectionElevation+"");
    		properties.setProperty(prefix+"projectionYaw",this.projectionYaw+"");
    		properties.setProperty(prefix+"projectionRoll",this.projectionRoll+"");
    		properties.setProperty(prefix+"matchPixelSize",this.matchPixelSize+"");
    		properties.setProperty(prefix+"projectionPixelSize",this.projectionPixelSize+"");
    		properties.setProperty(prefix+"projectionWidth",this.projectionWidth+"");
    		properties.setProperty(prefix+"projectionHeight",this.projectionHeight+"");
    		properties.setProperty(prefix+"projectionCenterX",this.projectionCenterX+"");
    		properties.setProperty(prefix+"projectionCenterY",this.projectionCenterY+"");
    		properties.setProperty(prefix+"nominalHorizontalDisparity",this.nominalHorizontalDisparity+"");
    		if (this.channelSelection!=null){
    			properties.setProperty(prefix+"numberProjectedChannels",this.channelSelection.length+"");
    			for (int i=0;i<this.channelSelection.length;i++){
    				properties.setProperty(prefix+"projectedChannel"+i,this.channelSelection[i]+"");
    			}
    		}

    	}
    	public void getProperties(String prefix,Properties properties){

    		if (properties.getProperty(prefix+"longitudeLeft")!=null) this.longitudeLeft=      Double.parseDouble(properties.getProperty(prefix+"longitudeLeft"));
    		if (properties.getProperty(prefix+"longitudeRight")!=null) this.longitudeRight=      Double.parseDouble(properties.getProperty(prefix+"longitudeRight"));
    		if (properties.getProperty(prefix+"latitudeTop")!=null)this.latitudeTop=       Double.parseDouble(properties.getProperty(prefix+"latitudeTop"));
    		if (properties.getProperty(prefix+"latitudeBottom")!=null)this.latitudeBottom=       Double.parseDouble(properties.getProperty(prefix+"latitudeBottom"));

    		if (properties.getProperty(prefix+"pixelsHorizontal")!=null)this.pixelsHorizontal=          Integer.parseInt(properties.getProperty(prefix+"pixelsHorizontal"));
    		if (properties.getProperty(prefix+"imageWidth")!=null)this.imageWidth=          Integer.parseInt(properties.getProperty(prefix+"imageWidth"));
    		if (properties.getProperty(prefix+"imageHeight")!=null)this.imageHeight=          Integer.parseInt(properties.getProperty(prefix+"imageHeight"));

    		if (properties.getProperty(prefix+"resolutionScale")!=null)this.resolutionScale=       Double.parseDouble(properties.getProperty(prefix+"resolutionScale"));
    		if (properties.getProperty(prefix+"x0")!=null)this.x0=       Double.parseDouble(properties.getProperty(prefix+"x0"));
    		if (properties.getProperty(prefix+"y0")!=null)this.y0=       Double.parseDouble(properties.getProperty(prefix+"y0"));

    		if (properties.getProperty(prefix+"longitudeWidth")!=null)this.longitudeWidth=          Integer.parseInt(properties.getProperty(prefix+"longitudeWidth"));

    		if (properties.getProperty(prefix+"clearFullMap")!=null)this.clearFullMap=   Boolean.parseBoolean(properties.getProperty(prefix+"clearFullMap"));
    		if (properties.getProperty(prefix+"clearAllMaps")!=null)this.clearAllMaps=   Boolean.parseBoolean(properties.getProperty(prefix+"clearAllMaps"));

    		if (properties.getProperty(prefix+"generateCommonPlane")!=null)this.generateCommonPlane=       Boolean.parseBoolean(properties.getProperty(prefix+"generateCommonPlane"));
    		if (properties.getProperty(prefix+"projectionElevation")!=null)this.projectionElevation=       Double.parseDouble(properties.getProperty(prefix+"projectionElevation"));
    		if (properties.getProperty(prefix+"projectionYaw")!=null)this.projectionYaw=       Double.parseDouble(properties.getProperty(prefix+"projectionYaw"));
    		if (properties.getProperty(prefix+"projectionRoll")!=null)this.projectionRoll=       Double.parseDouble(properties.getProperty(prefix+"projectionRoll"));
    		if (properties.getProperty(prefix+"matchPixelSize")!=null)this.matchPixelSize=   Boolean.parseBoolean(properties.getProperty(prefix+"matchPixelSize"));
    		if (properties.getProperty(prefix+"projectionPixelSize")!=null)this.projectionPixelSize=       Double.parseDouble(properties.getProperty(prefix+"projectionPixelSize"));

    		if (properties.getProperty(prefix+"projectionWidth")!=null)this.projectionWidth=          Integer.parseInt(properties.getProperty(prefix+"projectionWidth"));
    		if (properties.getProperty(prefix+"projectionHeight")!=null)this.projectionHeight=          Integer.parseInt(properties.getProperty(prefix+"projectionHeight"));

    		if (properties.getProperty(prefix+"projectionCenterX")!=null)this.projectionCenterX=       Double.parseDouble(properties.getProperty(prefix+"projectionCenterX"));
    		if (properties.getProperty(prefix+"projectionCenterY")!=null)this.projectionCenterY=       Double.parseDouble(properties.getProperty(prefix+"projectionCenterY"));
    		if (properties.getProperty(prefix+"nominalHorizontalDisparity")!=null)this.nominalHorizontalDisparity=       Double.parseDouble(properties.getProperty(prefix+"nominalHorizontalDisparity"));
    		if (properties.getProperty(prefix+"numberProjectedChannels")!=null){
    			int numberProjectedChannels=Integer.parseInt(properties.getProperty(prefix+"numberProjectedChannels"));
    			this.channelSelection=new boolean[numberProjectedChannels];
    			for (int i=0;i<this.channelSelection.length;i++){
    	    		if (properties.getProperty(prefix+"projectedChannel"+i)!=null)
    	    			this.channelSelection[i]=   Boolean.parseBoolean(properties.getProperty(prefix+"projectedChannel"+i));
    			}
    		}

    	}


    	public boolean showDialog() {
    		needRebuild=false;
			GenericDialog gd=new GenericDialog("Select parameters for equirectangular->sensor pixel mapping");
			gd.addMessage("Equirectangular area");
			gd.addNumericField("Longitude left", this.longitudeLeft, 1,6,"degrees" );
			gd.addNumericField("Longitude right", this.longitudeRight, 1,6,"degrees" );
			gd.addNumericField("Latitude top", this.latitudeTop, 1,6,"degrees" );
			gd.addNumericField("Latitude bottom", this.latitudeBottom, 1,6,"degrees" );
			gd.addNumericField("Pixels horizontal ", this.pixelsHorizontal,0,5,"image pix");
			gd.addMessage("Source image parameters");
			gd.addNumericField("Input image width", this.imageWidth,0,4,"image pix");
			gd.addNumericField("Input image height", this.imageHeight,0,4,"image pix");
			gd.addNumericField("Input image resolution scale (2.0 - twice resolution, 0.5 - half)", this.resolutionScale, 4,6,"x" );
			gd.addNumericField("Input image left margin", this.x0, 1,6,"sensor pix" );
			gd.addNumericField("Input image top margin", this.y0, 1,6,"sensor pix" );
			gd.addMessage("Reduction of memory usage");
			gd.addNumericField("Crop files horizontally to ", this.longitudeWidth,0,4,"longitude pix");
			gd.addCheckbox    ("Clear full map",    this.clearFullMap);
			gd.addCheckbox    ("Clear all data",    this.clearFullMap);
			gd.addMessage("Parameters for the common projection plane (binocular/trinocular cameras)");
			gd.addCheckbox    ("Generate common projection plane",    this.generateCommonPlane);
			gd.addNumericField("View axis elevation (orthogonal to projection plane)",this.projectionElevation,2,6,"degrees");
			gd.addNumericField("View axis heading   (orthogonal to projection plane)",this.projectionYaw,2,6,"degrees");
			gd.addNumericField("View plane rotation (roll) around the view axis",     this.projectionRoll,2,6,"degrees");
			gd.addCheckbox    ("Match projection pixel size to that of the equirectangular map",    this.matchPixelSize);
			gd.addNumericField("Projection pixel size (relative)     ",1000*this.projectionPixelSize,4,8,"x1/1000");
			gd.addNumericField("Projection plane width", this.projectionWidth,0,5,"pix");
			gd.addNumericField("Projection plane height",this.projectionHeight,0,5,"pix");
			gd.addNumericField("Projection plane Center X (point orthogonal to the view axis), right",   this.projectionCenterX,2,8,"pix");
			gd.addNumericField("Projection plane Center Y (point orthogonal to the view axis), down",    this.projectionCenterY,2,8,"pix");
			gd.addNumericField("Nominal distance between the 2 horizontal cameras",    this.nominalHorizontalDisparity,2,8,"mm");
			gd.enableYesNoCancel("OK", "Rebuild map files");
    		WindowTools.addScrollBars(gd);
			gd.showDialog();
			if (gd.wasCanceled()) return false;
			this.longitudeLeft=          gd.getNextNumber();
			this.longitudeRight=         gd.getNextNumber();
			this.latitudeTop=            gd.getNextNumber();
			this.latitudeBottom=         gd.getNextNumber();
			this.pixelsHorizontal= (int) gd.getNextNumber();
			this.imageWidth=      (int)  gd.getNextNumber();
			this.imageHeight=     (int)  gd.getNextNumber();
			this.resolutionScale=        gd.getNextNumber();
			this.x0=                     gd.getNextNumber();
			this.y0=                     gd.getNextNumber();
			this.longitudeWidth=   (int) gd.getNextNumber();
			this.clearFullMap=           gd.getNextBoolean();
			this.clearAllMaps=           gd.getNextBoolean();
    		this.generateCommonPlane =   gd.getNextBoolean();
			this.projectionElevation=    gd.getNextNumber();
			this.projectionYaw=          gd.getNextNumber();
			this.projectionRoll=         gd.getNextNumber();
			this.matchPixelSize =        gd.getNextBoolean();
			this.projectionPixelSize=    0.001*gd.getNextNumber();
			this.projectionWidth=  (int) gd.getNextNumber();
			this.projectionHeight= (int) gd.getNextNumber();
			this.projectionCenterX=      gd.getNextNumber();
			this.projectionCenterY=      gd.getNextNumber();
			this.nominalHorizontalDisparity=      gd.getNextNumber();

			if (!gd.wasOKed()) needRebuild=true;
    		return true;
    	}
    	public boolean selectChannelsToProcess(String title, int numChannels) {
    		if (numChannels<=0){
    			this.channelSelection=null;
    			return true;
    		}
    		boolean [] newSelecttion=new boolean [numChannels];
    		boolean lastChoice=true; // selected
    		for (int i=0;i<numChannels;i++){
    			if ((this.channelSelection!=null) && (i<this.channelSelection.length)){
    				newSelecttion[i]=this.channelSelection[i];
    				lastChoice=newSelecttion[i];
    			} else newSelecttion[i]=lastChoice;
    		}
    		while (true) {
    			GenericDialog gd = new GenericDialog(title);
    			for (int i=0;i<numChannels;i++) gd.addCheckbox("channel "+i, newSelecttion[i]);
    			gd.enableYesNoCancel("OK", "All like channel 0");
    			WindowTools.addScrollBars(gd);
    			gd.showDialog();
    			if (gd.wasCanceled()) return false; // but do not modify this.channelSelection
    			for (int i=0;i<numChannels;i++) newSelecttion[i]=gd.getNextBoolean();
    			if (gd.wasOKed()){
    				this.channelSelection=newSelecttion;
    				return true;
    			} else {
    				for (int i=1;i<numChannels;i++) newSelecttion[i]=newSelecttion[0];
    			}
    		}
    	}


    }



  /* ======================================================================== */

}
