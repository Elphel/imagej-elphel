package com.elphel.imagej.tileprocessor;
import java.awt.Rectangle;
/**
 ** TwoQuadCLT - Process images from a pair of Quad/Octal cameras
 **
 ** Copyright (C) 2018 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  TwoQuadCLT.java is free software: you can redistribute it and/or modify
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
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;
import java.nio.channels.Channels;
import java.nio.channels.FileChannel;
import java.nio.channels.WritableByteChannel;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Properties;
import java.util.Random;

import com.elphel.imagej.cameras.CLTParameters;
import com.elphel.imagej.cameras.ColorProcParameters;
import com.elphel.imagej.cameras.EyesisCorrectionParameters;
import com.elphel.imagej.common.GenericJTabbedDialog;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.correction.CorrectionColorProc;
import com.elphel.imagej.gpu.GPUTileProcessor;
import com.elphel.imagej.jp4.JP46_Reader_camera;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.io.FileSaver;

public class TwoQuadCLT {
	public static  int DSI_DISPARITY_MAIN = 0;
	public static  int DSI_DISPARITY_AUX =  1;
	public static  int DSI_DISPARITY_RIG =  2;
	public static  int DSI_DISPARITY_X3D =  3;
	public static  int DSI_STRENGTH_MAIN =  4;
	public static  int DSI_STRENGTH_AUX =   5;
	public static  int DSI_STRENGTH_RIG =   6;

	public static String DSI_COMBO_SUFFIX = "-DSI_COMBO";
	public static String DSI_MAIN_SUFFIX =  "-DSI_MAIN";

	public static String [] DSI_SLICES =
		{       "disparity_main",
				"disparity_aux",
				"disparity_rig",
				"disparity_x3d",
			    "strength_main",
			    "strength_aux",
			    "strength_rig"};

	public long                                            startTime;     // start of batch processing
	public long                                            startSetTime;  // start of set processing
	public long                                            startStepTime; // start of step processing
	public BiCamDSI                                        biCamDSI_persistent; // may be removed later to save memory, now to be able to continue
	PoleProcessor                                          poleProcessor_persistent;
	public QuadCLT quadCLT_main =  null;
	public QuadCLT quadCLT_aux =   null;
	public double [][] dsi = new double [DSI_SLICES.length][];
	public double [][] dsi_aux_from_main; // Main camera DSI converted into the coordinates of the AUX one (see QuadCLT.FGBG_TITLES)-added two slises from aux ds

	public TwoQuadCLT(
			QuadCLT quadCLT_main,
			QuadCLT quadCLT_aux)
	{
		this.quadCLT_main = quadCLT_main;
		this.quadCLT_aux =  quadCLT_aux;
	}
	public QuadCLT getMain()
	{
		return quadCLT_main;
	}
	public QuadCLT getAux()
	{
		return quadCLT_aux;
	}
	public void resetRig(boolean all_cams)
	{
		biCamDSI_persistent = null;
		poleProcessor_persistent = null;
		if (all_cams) {
			if ((quadCLT_main != null) && (quadCLT_main.tp != null)) {
				quadCLT_main.tp.resetCLTPasses();
			}
			if ((quadCLT_aux != null) && (quadCLT_aux.tp != null)) {
				quadCLT_aux.tp.resetCLTPasses();
			}

		}
	}

	public BiScan getBiScan(int indx) {
		if ((biCamDSI_persistent == null) ||
				(biCamDSI_persistent.biScans == null) ||
				(biCamDSI_persistent.biScans.size() <= indx)) {
			return null;
		}
		return biCamDSI_persistent.biScans.get(indx);
	}
/*
	public void processCLTTwoQuadCorrs( // not referenced?
			QuadCLT quadCLT_main,
			QuadCLT quadCLT_aux,
			CLTParameters       clt_parameters,
			EyesisCorrectionParameters.DebayerParameters   debayerParameters,
			ColorProcParameters                            colorProcParameters,
			ColorProcParameters                            colorProcParameters_aux,
			CorrectionColorProc.ColorGainsParameters       channelGainParameters_main,
			CorrectionColorProc.ColorGainsParameters       channelGainParameters_aux,
			EyesisCorrectionParameters.RGBParameters       rgbParameters,
			final int                                      threadsMax,  // maximal number of threads to launch
			final boolean                                  updateStatus,
			final int                                      debugLevel) throws Exception
	{

		this.startTime=System.nanoTime();
		String [] sourceFiles=quadCLT_main.correctionsParameters.getSourcePaths();
		QuadCLT.SetChannels [] set_channels_main = quadCLT_main.setChannels(debugLevel);
		QuadCLT.SetChannels [] set_channels_aux =  quadCLT_aux.setChannels(debugLevel);
		if ((set_channels_main == null) || (set_channels_main.length==0) || (set_channels_aux == null) || (set_channels_aux.length==0)) {
			System.out.println("No files to process (of "+sourceFiles.length+")");
			return;
		}
		double [] referenceExposures_main = quadCLT_main.eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		double [] referenceExposures_aux =  quadCLT_aux.eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		for (int nSet = 0; nSet < set_channels_main.length; nSet++){
			// check it is the same set for both cameras
			if (set_channels_aux.length <= nSet ) {
				throw new Exception ("Set names for cameras do not match: main camera: '"+set_channels_main[nSet].name()+"', aux. camera: nothing");
			}
			if (!set_channels_main[nSet].name().equals(set_channels_aux[nSet].name())) {
				throw new Exception ("Set names for cameras do not match: main camera: '"+set_channels_main[nSet].name()+"', aux. camera: '"+set_channels_main[nSet].name()+"'");
			}

			int [] channelFiles_main = set_channels_main[nSet].fileNumber();
			int [] channelFiles_aux =  set_channels_aux[nSet].fileNumber();
			boolean [][] saturation_imp_main = (clt_parameters.sat_level > 0.0)? new boolean[channelFiles_main.length][] : null;
			boolean [][] saturation_imp_aux =  (clt_parameters.sat_level > 0.0)? new boolean[channelFiles_main.length][] : null;
			double [] scaleExposures_main = new double[channelFiles_main.length];
			double [] scaleExposures_aux =  new double[channelFiles_main.length];

			ImagePlus [] imp_srcs_main = quadCLT_main.conditionImageSet(
					clt_parameters,                 // EyesisCorrectionParameters.CLTParameters  clt_parameters,
				    colorProcParameters,            //  ColorProcParameters                       colorProcParameters, //

					sourceFiles,                    // String []                                 sourceFiles,
					set_channels_main[nSet].name(), // String                                    set_name,
					referenceExposures_main,        // double []                                 referenceExposures,
					channelFiles_main,              // int []                                    channelFiles,
					scaleExposures_main,            //output  // double [] scaleExposures
					saturation_imp_main,            //output  // boolean [][]                              saturation_imp,
					debugLevel); // int                                       debugLevel);

			ImagePlus [] imp_srcs_aux = quadCLT_aux.conditionImageSet(
					clt_parameters,                 // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					colorProcParameters_aux,        //  ColorProcParameters                       colorProcParameters, //
					sourceFiles,                    // String []                                 sourceFiles,
					set_channels_aux[nSet].name(), // String                                    set_name,
					referenceExposures_aux,        // double []                                 referenceExposures,
					channelFiles_aux,              // int []                                    channelFiles,
					scaleExposures_aux,            //output  // double [] scaleExposures
					saturation_imp_aux,            //output  // boolean [][]                              saturation_imp,
					debugLevel); // int                                       debugLevel);

			// Temporarily processing individaully with the old code
			quadCLT_main.processCLTQuadCorr( // returns ImagePlus, but it already should be saved/shown
					imp_srcs_main, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
					saturation_imp_main, // boolean [][] saturation_imp, // (near) saturated pixels or null
					clt_parameters,
					debayerParameters,
					colorProcParameters,
					channelGainParameters_main,
					rgbParameters,
					scaleExposures_main,
					false, // apply_corr, // calculate and apply additional fine geometry correction
					false, // infinity_corr, // calculate and apply geometry correction at infinity
					threadsMax,  // maximal number of threads to launch
					updateStatus,
					debugLevel);

			Runtime.getRuntime().gc();
			if (debugLevel >-1) System.out.println("Processing set "+(nSet+1)+" (of "+set_channels_main.length+") finished at "+
					IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");

			quadCLT_aux.processCLTQuadCorr( // returns ImagePlus, but it already should be saved/shown
					imp_srcs_aux, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
					saturation_imp_aux, // boolean [][] saturation_imp, // (near) saturated pixels or null
					clt_parameters,
					debayerParameters,
					colorProcParameters_aux,
					channelGainParameters_aux,
					rgbParameters,
					scaleExposures_aux,
					false, // apply_corr, // calculate and apply additional fine geometry correction
					false, // infinity_corr, // calculate and apply geometry correction at infinity
					threadsMax,  // maximal number of threads to launch
					updateStatus,
					debugLevel);

			Runtime.getRuntime().gc();
			if (debugLevel >-1) System.out.println("Processing set "+(nSet+1)+" (of "+set_channels_aux.length+") finished at "+
					IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");

			if (quadCLT_aux.eyesisCorrections.stopRequested.get()>0) {
				System.out.println("User requested stop");
				System.out.println("Processing "+(nSet + 1)+" file sets (of "+set_channels_main.length+") finished at "+
						IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
				return;
			}
		}
		System.out.println("processCLTTwoQuadCorrs(): processing "+(quadCLT_main.getTotalFiles(set_channels_main)+quadCLT_aux.getTotalFiles(set_channels_aux))+" files ("+set_channels_main.length+" file sets) finished at "+
				IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");

	}

*/



	public void processCLTQuadCorrPairs(
			QuadCLT                                        quadCLT_main,
			QuadCLT                                        quadCLT_aux,
			CLTParameters       clt_parameters,
			EyesisCorrectionParameters.DebayerParameters   debayerParameters,
			ColorProcParameters                            colorProcParameters,
			ColorProcParameters                            colorProcParameters_aux,
			EyesisCorrectionParameters.RGBParameters       rgbParameters,
			final int                                      threadsMax,  // maximal number of threads to launch
			final boolean                                  updateStatus,
			final int                                      debugLevel) throws Exception
	{

		this.startTime=System.nanoTime();
		String [] sourceFiles=quadCLT_main.correctionsParameters.getSourcePaths();
		QuadCLT.SetChannels [] set_channels_main = quadCLT_main.setChannels(debugLevel);
		QuadCLT.SetChannels [] set_channels_aux =  quadCLT_aux.setChannels(debugLevel);
		if ((set_channels_main == null) || (set_channels_main.length==0) || (set_channels_aux == null) || (set_channels_aux.length==0)) {
			System.out.println("No files to process (of "+sourceFiles.length+")");
			return;
		}
		double [] referenceExposures_main = quadCLT_main.eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		double [] referenceExposures_aux =  quadCLT_aux.eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		for (int nSet = 0; nSet < set_channels_main.length; nSet++){
			// check it is the same set for both cameras
			if (set_channels_aux.length <= nSet ) {
				throw new Exception ("Set names for cameras do not match: main camera: '"+set_channels_main[nSet].name()+"', aux. camera: nothing");
			}
			if (!set_channels_main[nSet].name().equals(set_channels_aux[nSet].name())) {
				throw new Exception ("Set names for cameras do not match: main camera: '"+set_channels_main[nSet].name()+"', aux. camera: '"+set_channels_main[nSet].name()+"'");
			}

			int [] channelFiles_main = set_channels_main[nSet].fileNumber();
			int [] channelFiles_aux =  set_channels_aux[nSet].fileNumber();
			boolean [][] saturation_imp_main = (clt_parameters.sat_level > 0.0)? new boolean[channelFiles_main.length][] : null;
			boolean [][] saturation_imp_aux =  (clt_parameters.sat_level > 0.0)? new boolean[channelFiles_main.length][] : null;
			double [] scaleExposures_main = new double[channelFiles_main.length];
			double [] scaleExposures_aux =  new double[channelFiles_main.length];

			ImagePlus [] imp_srcs_main = quadCLT_main.conditionImageSet(
					clt_parameters,                 // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					colorProcParameters,            //  ColorProcParameters                       colorProcParameters, //
					sourceFiles,                    // String []                                 sourceFiles,
					set_channels_main[nSet].name(), // String                                    set_name,
					referenceExposures_main,        // double []                                 referenceExposures,
					channelFiles_main,              // int []                                    channelFiles,
					scaleExposures_main,            //output  // double [] scaleExposures
					saturation_imp_main,            //output  // boolean [][]                              saturation_imp,
					threadsMax,                 // int                                       threadsMax,
					debugLevel); // int                                       debugLevel);

			ImagePlus [] imp_srcs_aux = quadCLT_aux.conditionImageSet(
					clt_parameters,                 // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					colorProcParameters_aux,        //  ColorProcParameters                       colorProcParameters, //
					sourceFiles,                    // String []                                 sourceFiles,
					set_channels_aux[nSet].name(), // String                                    set_name,
					referenceExposures_aux,        // double []                                 referenceExposures,
					channelFiles_aux,              // int []                                    channelFiles,
					scaleExposures_aux,            //output  // double [] scaleExposures
					saturation_imp_aux,            //output  // boolean [][]                              saturation_imp,
					threadsMax,                 // int                                       threadsMax,
					debugLevel); // int                                       debugLevel);

			// Tempporarily processing individaully with the old code
			processCLTQuadCorrPair(
					quadCLT_main,               // QuadCLT                                        quadCLT_main,
					quadCLT_aux,                // QuadCLT                                        quadCLT_aux,
					imp_srcs_main,              // ImagePlus []                                   imp_quad_main,
					imp_srcs_aux,               // ImagePlus []                                   imp_quad_aux,
					saturation_imp_main,        // boolean [][]                                   saturation_main, // (near) saturated pixels or null
					saturation_imp_aux,         // boolean [][]                                   saturation_aux, // (near) saturated pixels or null
					clt_parameters,             // EyesisCorrectionParameters.CLTParameters       clt_parameters,
					debayerParameters,          // EyesisCorrectionParameters.DebayerParameters   debayerParameters,
					colorProcParameters,        // EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
					colorProcParameters_aux,
					rgbParameters,              // EyesisCorrectionParameters.RGBParameters       rgbParameters,
					scaleExposures_main,        // double []	                                     scaleExposures_main, // probably not needed here - restores brightness of the final image
					scaleExposures_aux,         // double []	                                     scaleExposures_aux, // probably not needed here - restores brightness of the final image
					false,                      //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
					// averages measurements
					clt_parameters.rig.lt_avg_radius,// final int                                      lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using

					//			  final boolean    apply_corr, // calculate and apply additional fine geometry correction
					//			  final boolean    infinity_corr, // calculate and apply geometry correction at infinity
					threadsMax,                 // final int        threadsMax,  // maximal number of threads to launch
					updateStatus,               // final boolean    updateStatus,
					debugLevel);                // final int        debugLevel);

			Runtime.getRuntime().gc();
			if (debugLevel >-1) System.out.println("Processing set "+(nSet+1)+" (of "+set_channels_aux.length+") finished at "+
					IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");

			if (quadCLT_aux.eyesisCorrections.stopRequested.get()>0) {
				System.out.println("User requested stop");
				System.out.println("Processing "+(nSet + 1)+" file sets (of "+set_channels_main.length+") finished at "+
						IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
				return;
			}
		}
		System.out.println("processCLTQuadCorrPairs(): processing "+(quadCLT_main.getTotalFiles(set_channels_main)+quadCLT_aux.getTotalFiles(set_channels_aux))+" files ("+set_channels_main.length+" file sets) finished at "+
				IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");

	}

	public void prepareFilesForGPUDebug(
			String                                         save_prefix, // absolute path to the cuda project root
			QuadCLT                                        quadCLT_main,
			QuadCLT                                        quadCLT_aux,
			CLTParameters       clt_parameters,
			EyesisCorrectionParameters.DebayerParameters   debayerParameters,
			ColorProcParameters                            colorProcParameters,
			ColorProcParameters                            colorProcParameters_aux,
			EyesisCorrectionParameters.RGBParameters       rgbParameters,
			final int                                      threadsMax,  // maximal number of threads to launch
			final boolean                                  updateStatus,
			final int                                      debugLevel) throws Exception
	{

		this.startTime=System.nanoTime();
		String [] sourceFiles=quadCLT_main.correctionsParameters.getSourcePaths();
		QuadCLT.SetChannels [] set_channels_main = quadCLT_main.setChannels(debugLevel);
		QuadCLT.SetChannels [] set_channels_aux =  quadCLT_aux.setChannels(debugLevel);
		if ((set_channels_main == null) || (set_channels_main.length==0) || (set_channels_aux == null) || (set_channels_aux.length==0)) {
			System.out.println("No files to process (of "+sourceFiles.length+")");
			return;
		}
		double [] referenceExposures_main = quadCLT_main.eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		double [] referenceExposures_aux =  quadCLT_aux.eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		for (int nSet = 0; nSet < set_channels_main.length; nSet++){
			// check it is the same set for both cameras
			if (set_channels_aux.length <= nSet ) {
				throw new Exception ("Set names for cameras do not match: main camera: '"+set_channels_main[nSet].name()+"', aux. camera: nothing");
			}
			if (!set_channels_main[nSet].name().equals(set_channels_aux[nSet].name())) {
				throw new Exception ("Set names for cameras do not match: main camera: '"+set_channels_main[nSet].name()+"', aux. camera: '"+set_channels_main[nSet].name()+"'");
			}

			int [] channelFiles_main = set_channels_main[nSet].fileNumber();
			int [] channelFiles_aux =  set_channels_aux[nSet].fileNumber();
			boolean [][] saturation_imp_main = (clt_parameters.sat_level > 0.0)? new boolean[channelFiles_main.length][] : null;
			boolean [][] saturation_imp_aux =  (clt_parameters.sat_level > 0.0)? new boolean[channelFiles_main.length][] : null;
			double [] scaleExposures_main = new double[channelFiles_main.length];
			double [] scaleExposures_aux =  new double[channelFiles_main.length];

			ImagePlus [] imp_srcs_main = quadCLT_main.conditionImageSet(
					clt_parameters,                 // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					colorProcParameters,            //  ColorProcParameters                       colorProcParameters, //
					sourceFiles,                    // String []                                 sourceFiles,
					set_channels_main[nSet].name(), // String                                    set_name,
					referenceExposures_main,        // double []                                 referenceExposures,
					channelFiles_main,              // int []                                    channelFiles,
					scaleExposures_main,            //output  // double [] scaleExposures
					saturation_imp_main,            //output  // boolean [][]                              saturation_imp,
					threadsMax,                 // int                                       threadsMax,
					debugLevel); // int                                       debugLevel);

			ImagePlus [] imp_srcs_aux = quadCLT_aux.conditionImageSet(
					clt_parameters,                 // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					colorProcParameters_aux,        //  ColorProcParameters                       colorProcParameters, //
					sourceFiles,                    // String []                                 sourceFiles,
					set_channels_aux[nSet].name(), // String                                    set_name,
					referenceExposures_aux,        // double []                                 referenceExposures,
					channelFiles_aux,              // int []                                    channelFiles,
					scaleExposures_aux,            //output  // double [] scaleExposures
					saturation_imp_aux,            //output  // boolean [][]                              saturation_imp,
					threadsMax,                 // int                                       threadsMax,
					debugLevel); // int                                       debugLevel);

			// Tempporarily processing individually with the old code
			processCLTQuadCorrPairForGPU(
					save_prefix,                // String save_prefix,
					quadCLT_main,               // QuadCLT                                        quadCLT_main,
					quadCLT_aux,                // QuadCLT                                        quadCLT_aux,
					imp_srcs_main,              // ImagePlus []                                   imp_quad_main,
					imp_srcs_aux,               // ImagePlus []                                   imp_quad_aux,
					saturation_imp_main,        // boolean [][]                                   saturation_main, // (near) saturated pixels or null
					saturation_imp_aux,         // boolean [][]                                   saturation_aux, // (near) saturated pixels or null
					clt_parameters,             // EyesisCorrectionParameters.CLTParameters       clt_parameters,
					debayerParameters,          // EyesisCorrectionParameters.DebayerParameters   debayerParameters,
					colorProcParameters,        // EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
					colorProcParameters_aux,        // EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
					rgbParameters,              // EyesisCorrectionParameters.RGBParameters       rgbParameters,
					scaleExposures_main,        // double []	                                     scaleExposures_main, // probably not needed here - restores brightness of the final image
					scaleExposures_aux,         // double []	                                     scaleExposures_aux, // probably not needed here - restores brightness of the final image
					false,                      //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
					// averages measurements
					clt_parameters.rig.lt_avg_radius,// final int                                      lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using

					//			  final boolean    apply_corr, // calculate and apply additional fine geometry correction
					//			  final boolean    infinity_corr, // calculate and apply geometry correction at infinity
					threadsMax,                 // final int        threadsMax,  // maximal number of threads to launch
					updateStatus,               // final boolean    updateStatus,
					debugLevel);                // final int        debugLevel);

			Runtime.getRuntime().gc();
			if (debugLevel >-1) System.out.println("Processing set "+(nSet+1)+" (of "+set_channels_aux.length+") finished at "+
					IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");

			if (quadCLT_aux.eyesisCorrections.stopRequested.get()>0) {
				System.out.println("User requested stop");
				System.out.println("Processing "+(nSet + 1)+" file sets (of "+set_channels_main.length+") finished at "+
						IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
				return;
			}
		}
		System.out.println("prepareFilesForGPUDebug(): processing "+(quadCLT_main.getTotalFiles(set_channels_main)+quadCLT_aux.getTotalFiles(set_channels_aux))+" files ("+set_channels_main.length+" file sets) finished at "+
				IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");

	}

	public void processCLTQuadCorrPairsGpu(
			GPUTileProcessor.GpuQuad                        gpuQuad_main,
			GPUTileProcessor.GpuQuad                        gpuQuad_aux,
			QuadCLT                                         quadCLT_main,
			QuadCLT                                         quadCLT_aux,
			CLTParameters       clt_parameters,
			EyesisCorrectionParameters.CorrectionParameters ecp,
			EyesisCorrectionParameters.DebayerParameters    debayerParameters,
			ColorProcParameters                             colorProcParameters,
			ColorProcParameters                             colorProcParameters_aux,
			EyesisCorrectionParameters.RGBParameters        rgbParameters,
			final int                                       threadsMax,  // maximal number of threads to launch
			final boolean                                   updateStatus,
			final int                                       debugLevel) throws Exception
	{
//		TwoQuadGPU twoQuadGPU = new TwoQuadGPU(this);
		this.startTime=System.nanoTime();
		String [] sourceFiles=quadCLT_main.correctionsParameters.getSourcePaths();
		QuadCLT.SetChannels [] set_channels_main = quadCLT_main.setChannels(debugLevel);
		QuadCLT.SetChannels [] set_channels_aux =  quadCLT_aux.setChannels(debugLevel);
		if ((set_channels_main == null) || (set_channels_main.length==0) || (set_channels_aux == null) || (set_channels_aux.length==0)) {
			System.out.println("No files to process (of "+sourceFiles.length+")");
			return;
		}
		double [] referenceExposures_main = quadCLT_main.eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		double [] referenceExposures_aux =  quadCLT_aux.eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		for (int nSet = 0; nSet < set_channels_main.length; nSet++){
			// check it is the same set for both cameras
			if (set_channels_aux.length <= nSet ) {
				throw new Exception ("Set names for cameras do not match: main camera: '"+set_channels_main[nSet].name()+"', aux. camera: nothing");
			}
			if (!set_channels_main[nSet].name().equals(set_channels_aux[nSet].name())) {
				throw new Exception ("Set names for cameras do not match: main camera: '"+set_channels_main[nSet].name()+"', aux. camera: '"+set_channels_main[nSet].name()+"'");
			}

			int [] channelFiles_main = set_channels_main[nSet].fileNumber();
			int [] channelFiles_aux =  set_channels_aux[nSet].fileNumber();
			boolean [][] saturation_imp_main = (clt_parameters.sat_level > 0.0)? new boolean[channelFiles_main.length][] : null;
			boolean [][] saturation_imp_aux =  (clt_parameters.sat_level > 0.0)? new boolean[channelFiles_main.length][] : null;
			double [] scaleExposures_main = new double[channelFiles_main.length];
			double [] scaleExposures_aux =  new double[channelFiles_main.length];

			ImagePlus [] imp_srcs_main = quadCLT_main.conditionImageSet(
					clt_parameters,                 // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					colorProcParameters,            //  ColorProcParameters                       colorProcParameters, //
					sourceFiles,                    // String []                                 sourceFiles,
					set_channels_main[nSet].name(), // String                                    set_name,
					referenceExposures_main,        // double []                                 referenceExposures,
					channelFiles_main,              // int []                                    channelFiles,
					scaleExposures_main,            //output  // double [] scaleExposures
					saturation_imp_main,            //output  // boolean [][]                              saturation_imp,
					threadsMax,                 // int                                       threadsMax,
					debugLevel); // int                                       debugLevel);

			ImagePlus [] imp_srcs_aux = quadCLT_aux.conditionImageSet(
					clt_parameters,                 // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					colorProcParameters_aux,        //  ColorProcParameters                       colorProcParameters, //
					sourceFiles,                    // String []                                 sourceFiles,
					set_channels_aux[nSet].name(), // String                                    set_name,
					referenceExposures_aux,        // double []                                 referenceExposures,
					channelFiles_aux,              // int []                                    channelFiles,
					scaleExposures_aux,            //output  // double [] scaleExposures
					saturation_imp_aux,            //output  // boolean [][]                              saturation_imp,
					threadsMax,                 // int                                       threadsMax,
					debugLevel); // int                                       debugLevel);

			// Tempporarily processing individaully with the old code
			QuadCLT.processCLTQuadCorrPairGpu(
//					gpuQuad_main,               // GPUTileProcessor.GpuQuad                       gpuQuad_main,
//					gpuQuad_aux,                // GPUTileProcessor.GpuQuad                       gpuQuad_aux,
					quadCLT_main,               // QuadCLT                                        quadCLT_main,
					quadCLT_aux,                // QuadCLT                                        quadCLT_aux,
					imp_srcs_main,              // ImagePlus []                                   imp_quad_main,
					imp_srcs_aux,               // ImagePlus []                                   imp_quad_aux,
					saturation_imp_main,        // boolean [][]                                   saturation_main, // (near) saturated pixels or null
					saturation_imp_aux,         // boolean [][]                                   saturation_aux, // (near) saturated pixels or null
					clt_parameters,             // EyesisCorrectionParameters.CLTParameters       clt_parameters,
					ecp,                        // EyesisCorrectionParameters.CorrectionParameters ecp,
					debayerParameters,          // EyesisCorrectionParameters.DebayerParameters   debayerParameters,
					colorProcParameters,        // EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
					colorProcParameters_aux,    // EyesisCorrectionParameters.ColorProcParameters colorProcParameters_aux,
					rgbParameters,              // EyesisCorrectionParameters.RGBParameters       rgbParameters,
					scaleExposures_main,        // double []	                                     scaleExposures_main, // probably not needed here - restores brightness of the final image
					scaleExposures_aux,         // double []	                                     scaleExposures_aux, // probably not needed here - restores brightness of the final image
					false,                      //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
					// averages measurements
					clt_parameters.rig.lt_avg_radius,// final int                                      lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using

					//			  final boolean    apply_corr, // calculate and apply additional fine geometry correction
					//			  final boolean    infinity_corr, // calculate and apply geometry correction at infinity
					threadsMax,                 // final int        threadsMax,  // maximal number of threads to launch
					updateStatus,               // final boolean    updateStatus,
					debugLevel);                // final int        debugLevel);

			Runtime.getRuntime().gc();
			if (debugLevel >-1) System.out.println("Processing set "+(nSet+1)+" (of "+set_channels_aux.length+") finished at "+
					IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");

			if (quadCLT_aux.eyesisCorrections.stopRequested.get()>0) {
				System.out.println("User requested stop");
				System.out.println("Processing "+(nSet + 1)+" file sets (of "+set_channels_main.length+") finished at "+
						IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
				return;
			}
		}
		System.out.println("processCLTQuadCorrPairsGpu(): processing "+(quadCLT_main.getTotalFiles(set_channels_main)+quadCLT_aux.getTotalFiles(set_channels_aux))+" files ("+set_channels_main.length+" file sets) finished at "+
				IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");

	}
	
	public ImagePlus [] processCLTQuadCorrPair(
			QuadCLT                                        quadCLT_main,
			QuadCLT                                        quadCLT_aux,
			ImagePlus []                                   imp_quad_main,
			ImagePlus []                                   imp_quad_aux,
			boolean [][]                                   saturation_main, // (near) saturated pixels or null
			boolean [][]                                   saturation_aux, // (near) saturated pixels or null
			CLTParameters       clt_parameters,
			EyesisCorrectionParameters.DebayerParameters   debayerParameters,
			ColorProcParameters                            colorProcParameters,
			ColorProcParameters                            colorProcParameters_aux,
			EyesisCorrectionParameters.RGBParameters       rgbParameters,
			double []	                                     scaleExposures_main, // probably not needed here - restores brightness of the final image
			double []	                                     scaleExposures_aux, // probably not needed here - restores brightness of the final image
			boolean                                        notch_mode, // use pole-detection mode for inter-camera correlation
			final int                                      lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
			final int        threadsMax,  // maximal number of threads to launch
			final boolean    updateStatus,
			final int        debugLevel){
		final boolean      batch_mode = clt_parameters.batch_run; //disable any debug images
		final boolean get_ers = !batch_mode;
		//		  boolean batch_mode = false;
		boolean infinity_corr = false;
		double [][] scaleExposures= {scaleExposures_main, scaleExposures_aux};
		boolean toRGB=     quadCLT_main.correctionsParameters.toRGB;
		ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging? - TODO - move where it belongs
		// may use this.StartTime to report intermediate steps execution times
		String name=quadCLT_main.correctionsParameters.getModelName((String) imp_quad_main[0].getProperty("name"));
		String path= (String) imp_quad_main[0].getProperty("path"); // Only for debug output
		ImagePlus [] results = new ImagePlus[imp_quad_main.length + imp_quad_main.length];
		for (int i = 0; i < results.length; i++) {
			if (i< imp_quad_main.length) {
				results[i] = imp_quad_main[i];
			} else {
				results[i] = imp_quad_main[i-imp_quad_main.length];
			}
			results[i].setTitle(results[i].getTitle()+"RAW");
		}
		if (debugLevel>1) System.out.println("processing: "+path);
/*		// 08/12/2020 Moved to conditionImageSet		
		getRigImageStacks(
				clt_parameters,  // EyesisCorrectionParameters.CLTParameters       clt_parameters,
				quadCLT_main,    // QuadCLT                                         quadCLT_main,
				quadCLT_aux,     // QuadCLT                                          quadCLT_aux,
				imp_quad_main,   // ImagePlus []                                   imp_quad_main,
				imp_quad_aux,    // ImagePlus []                                    imp_quad_aux,
				saturation_main, // boolean [][]        saturation_main, // (near) saturated pixels or null
				saturation_aux, // boolean [][]        saturation_aux, // (near) saturated pixels or null
				threadsMax,      // maximal number of threads to launch
				debugLevel);     // final int        debugLevel);
*/
		// temporary setting up tile task file (one integer per tile, bitmask
		// for testing defined for a window, later the tiles to process will be calculated based on previous passes results

		int [][]    tile_op_main = quadCLT_main.tp.setSameTileOp(clt_parameters,  clt_parameters.tile_task_op, debugLevel);
		//		  int [][]    tile_op_aux =  quadCLT_aux.tp.setSameTileOp (clt_parameters,  clt_parameters.tile_task_op, debugLevel);

		double [][] disparity_array_main = quadCLT_main.tp.setSameDisparity(clt_parameters.disparity); // [tp.tilesY][tp.tilesX] - individual per-tile expected disparity
		//TODO: Add array of default disparity - use for combining images in force disparity mode (no correlation), when disparity is predicted from other tiles
		// start with all old arrays, remove some later
		double [][][][]     clt_corr_combo =        null;
		double [][][][]     texture_tiles_main =    null; // [tp.tilesY][tp.tilesX]["RGBA".length()][]; // tiles will be 16x16, 2 visualization mode full 16 or overlapped
		double [][][][]     texture_tiles_aux =     null; // [tp.tilesY][tp.tilesX]["RGBA".length()][]; // tiles will be 16x16, 2 visualization mode full 16 or overlapped
		// undecided, so 2 modes of combining alpha - same as rgb, or use center tile only
		// Assuming same size of both image sets
		final int tilesX = quadCLT_main.tp.getTilesX();
		final int tilesY = quadCLT_main.tp.getTilesY();

		final boolean keep_lt_corr_combo = clt_parameters.corr_keep;
		// FIXME: Remove those huge arrays when not needed!
		if (clt_parameters.correlate){
			if (keep_lt_corr_combo) clt_corr_combo =    new double [ImageDtt.TCORR_TITLES.length][tilesY][tilesX][];
			texture_tiles_main =     new double [tilesY][tilesX][][]; // ["RGBA".length()][];
			texture_tiles_aux =     new double [tilesY][tilesX][][]; // ["RGBA".length()][];
			for (int i = 0; i < tilesY; i++){
				for (int j = 0; j < tilesX; j++){
					if (keep_lt_corr_combo) for (int k = 0; k<clt_corr_combo.length; k++) clt_corr_combo[k][i][j] = null;
					texture_tiles_main[i][j] = null;
					texture_tiles_aux[i][j] = null;
				}
			}
		}

		double [][] disparity_bimap  = new double [ImageDtt.BIDISPARITY_TITLES.length][]; //[0] -residual disparity, [1] - orthogonal (just for debugging) last 4 - max pixel differences

		ImageDtt image_dtt = new ImageDtt(
				clt_parameters.transform_size,
				quadCLT_main.isMonochrome(),
				quadCLT_main.isLwir(),
				clt_parameters.getScaleStrength(false));

		double [][] ml_data = null;
//		int [][] woi_tops = {quadCLT_main.woi_tops,quadCLT_aux.woi_tops};
		final double [][][]       ers_delay = get_ers?(new double [2][][]):null;


		final double [][][][][][][] clt_bidata = // new double[2][quad][nChn][tilesY][tilesX][][]; // first index - main/aux
				image_dtt.clt_bi_quad (
						clt_parameters,                       // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
						clt_parameters.getFatZero(image_dtt.isMonochrome()),  // final double              fatzero,         // May use correlation fat zero from 2 different parameters - fat_zero and rig.ml_fatzero
						notch_mode,                           //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
						lt_rad,                               // final int                                      lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
						// first measurement - use default setting
						clt_parameters.rig.no_int_x0, // boolean   no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
						tile_op_main,                         // final int [][]            tile_op_main,    // [tilesY][tilesX] - what to do - 0 - nothing for this tile
						disparity_array_main,                 // final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
						quadCLT_main.image_data,              // final double [][][]       image_data_main, // first index - number of image in a quad
						quadCLT_aux.image_data,               // final double [][][]       image_data_aux,  // first index - number of image in a quad
						quadCLT_main.saturation_imp,          // final boolean [][]        saturation_main, // (near) saturated pixels or null
						quadCLT_aux.saturation_imp,           // final boolean [][]        saturation_aux,  // (near) saturated pixels or null
						// correlation results - combo will be for the correation between two quad cameras
						clt_corr_combo,                       // final double [][][][]     clt_corr_combo,  // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
						// [type][tilesY][tilesX] should be set by caller
						// types: 0 - selected correlation (product+offset), 1 - sum
						disparity_bimap,                      // final double [][]    disparity_bimap, // [23][tilesY][tilesX]
						ml_data,                              // 	final double [][]         ml_data,         // data for ML - 10 layers - 4 center areas (3x3, 5x5,..) per camera-per direction, 1 - composite, and 1 with just 1 data (target disparity)
						texture_tiles_main,                   // final double [][][][]     texture_tiles_main, // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
						texture_tiles_aux,                    // final double [][][][]     texture_tiles_aux,  // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
						imp_quad_main[0].getWidth(),          // final int                 width,
						quadCLT_main.getGeometryCorrection(), // final GeometryCorrection  geometryCorrection_main,
						quadCLT_aux.getGeometryCorrection(),  // final GeometryCorrection  geometryCorrection_aux,
						quadCLT_main.getCLTKernels(),         // final double [][][][][][] clt_kernels_main, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
						quadCLT_aux.getCLTKernels(),          // final double [][][][][][] clt_kernels_aux,  // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
						clt_parameters.corr_magic_scale,      // final double              corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
						true,                                 // 	final boolean             keep_clt_data,
//						woi_tops,                             // final int [][]            woi_tops,
						ers_delay,                            // final double [][][]       ers_delay,        // if not null - fill with tile center acquisition delay
						threadsMax,                           // final int                 threadsMax,  // maximal number of threads to launch
						debugLevel);                          // final int                 globalDebugLevel);


		if (ers_delay !=null) {
			showERSDelay(ers_delay);
		}

		double [][] texture_nonoverlap_main = null;
		double [][] texture_nonoverlap_aux = null;
		double [][] texture_overlap_main = null;
		double [][] texture_overlap_aux = null;
		String [] rgba_titles = {"red","blue","green","alpha"};
		String [] rgba_weights_titles = {"red","blue","green","alpha","port0","port1","port2","port3","r-rms","b-rms","g-rms","w-rms"};
		if ((texture_tiles_main != null) && (texture_tiles_aux != null)){
			if (clt_parameters.show_nonoverlap){
				texture_nonoverlap_main = image_dtt.combineRBGATiles(
						texture_tiles_main,                 // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}
//						image_dtt.transform_size,
						false,                         // when false - output each tile as 16x16, true - overlap to make 8x8
						clt_parameters.sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
						threadsMax,                    // maximal number of threads to launch
						debugLevel);
				sdfa_instance.showArrays(
						texture_nonoverlap_main,
						tilesX * (2 * image_dtt.transform_size),
						tilesY * (2 * image_dtt.transform_size),
						true,
						name + "-TXTNOL-D"+clt_parameters.disparity+"-MAIN",
						(clt_parameters.keep_weights?rgba_weights_titles:rgba_titles));

				texture_nonoverlap_aux = image_dtt.combineRBGATiles(
						texture_tiles_aux,                 // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}
//						image_dtt.transform_size,
						false,                         // when false - output each tile as 16x16, true - overlap to make 8x8
						clt_parameters.sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
						threadsMax,                    // maximal number of threads to launch
						debugLevel);
				sdfa_instance.showArrays(
						texture_nonoverlap_aux,
						tilesX * (2 * image_dtt.transform_size),
						tilesY * (2 * image_dtt.transform_size),
						true,
						name + "-TXTNOL-D"+clt_parameters.disparity+"-AUX",
						(clt_parameters.keep_weights?rgba_weights_titles:rgba_titles));
			}
			if (!infinity_corr && (clt_parameters.show_overlap || clt_parameters.show_rgba_color)){
				int alpha_index = 3;
				texture_overlap_main = image_dtt.combineRBGATiles(
						texture_tiles_main,                 // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}
//						image_dtt.transform_size,
						true,                         // when false - output each tile as 16x16, true - overlap to make 8x8
						clt_parameters.sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
						threadsMax,                    // maximal number of threads to launch
						debugLevel);
				if (clt_parameters.alpha1 > 0){ // negative or 0 - keep alpha as it was
					double scale = (clt_parameters.alpha1 > clt_parameters.alpha0) ? (1.0/(clt_parameters.alpha1 - clt_parameters.alpha0)) : 0.0;
					for (int i = 0; i < texture_overlap_main[alpha_index].length; i++){
						double d = texture_overlap_main[alpha_index][i];
						if      (d >=clt_parameters.alpha1) d = 1.0;
						else if (d <=clt_parameters.alpha0) d = 0.0;
						else d = scale * (d- clt_parameters.alpha0);
						texture_overlap_main[alpha_index][i] = d;
					}
				}

				texture_overlap_aux = image_dtt.combineRBGATiles(
						texture_tiles_aux,                 // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}
//						image_dtt.transform_size,
						true,                         // when false - output each tile as 16x16, true - overlap to make 8x8
						clt_parameters.sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
						threadsMax,                    // maximal number of threads to launch
						debugLevel);
				if (clt_parameters.alpha1 > 0){ // negative or 0 - keep alpha as it was
					double scale = (clt_parameters.alpha1 > clt_parameters.alpha0) ? (1.0/(clt_parameters.alpha1 - clt_parameters.alpha0)) : 0.0;
					for (int i = 0; i < texture_overlap_aux[alpha_index].length; i++){
						double d = texture_overlap_aux[alpha_index][i];
						if      (d >=clt_parameters.alpha1) d = 1.0;
						else if (d <=clt_parameters.alpha0) d = 0.0;
						else d = scale * (d- clt_parameters.alpha0);
						texture_overlap_aux[alpha_index][i] = d;
					}
				}

				if (!batch_mode && clt_parameters.show_overlap) {
					sdfa_instance.showArrays(
							texture_overlap_main,
							tilesX * image_dtt.transform_size,
							tilesY * image_dtt.transform_size,
							true,
							name + "-TXTOL-D"+clt_parameters.disparity+"-MAIN",
							(clt_parameters.keep_weights?rgba_weights_titles:rgba_titles));
				}
				if (!batch_mode && clt_parameters.show_overlap) {
					sdfa_instance.showArrays(
							texture_overlap_aux,
							tilesX * image_dtt.transform_size,
							tilesY * image_dtt.transform_size,
							true,
							name + "-TXTOL-D"+clt_parameters.disparity+"-AUX",
							(clt_parameters.keep_weights?rgba_weights_titles:rgba_titles));
				}


				if (!batch_mode && clt_parameters.show_rgba_color) {
					// for now - use just RGB. Later add option for RGBA
					double [][] texture_rgb_main = {texture_overlap_main[0],texture_overlap_main[1],texture_overlap_main[2]};
					double [][] texture_rgba_main = {texture_overlap_main[0],texture_overlap_main[1],texture_overlap_main[2],texture_overlap_main[3]};
					double [][] texture_rgb_aux = {texture_overlap_aux[0],texture_overlap_aux[1],texture_overlap_aux[2]};
					double [][] texture_rgba_aux = {texture_overlap_aux[0],texture_overlap_aux[1],texture_overlap_aux[2],texture_overlap_aux[3]};
					ImagePlus imp_texture_main = quadCLT_main.linearStackToColor(
							clt_parameters,
							colorProcParameters,
							rgbParameters,
							name+"-texture", // String name,
							"-D"+clt_parameters.disparity+"-MAIN", //String suffix, // such as disparity=...
							toRGB,
							!quadCLT_main.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
							false, // true, // boolean saveShowIntermediate, // save/show if set globally
							false, // true, // boolean saveShowFinal,        // save/show result (color image?)
							((clt_parameters.alpha1 > 0)? texture_rgba_main: texture_rgb_main),
							tilesX *  image_dtt.transform_size,
							tilesY *  image_dtt.transform_size,
							1.0,         // double scaleExposure, // is it needed?
							debugLevel );
					ImagePlus imp_texture_aux = quadCLT_aux.linearStackToColor(
							clt_parameters,
							colorProcParameters_aux,
							rgbParameters,
							name+"-texture", // String name,
							"-D"+clt_parameters.disparity+"-AUX", //String suffix, // such as disparity=...
							toRGB,
							!quadCLT_aux.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
							false, // true, // boolean saveShowIntermediate, // save/show if set globally
							false, // true, // boolean saveShowFinal,        // save/show result (color image?)
							((clt_parameters.alpha1 > 0)? texture_rgba_aux: texture_rgb_aux),
							tilesX *  image_dtt.transform_size,
							tilesY *  image_dtt.transform_size,
							1.0,         // double scaleExposure, // is it needed?
							debugLevel );
					int width = imp_texture_main.getWidth();
					int height =imp_texture_main.getHeight();
					ImageStack texture_stack=new ImageStack(width,height);
					texture_stack.addSlice("main",      imp_texture_main.getProcessor().getPixels());
					texture_stack.addSlice("auxiliary", imp_texture_aux. getProcessor().getPixels());
					ImagePlus imp_texture_stack = new ImagePlus(name+"-TEXTURES-D"+clt_parameters.disparity, texture_stack);
					imp_texture_stack.getProcessor().resetMinAndMax();
					//					  imp_texture_stack.updateAndDraw();
					imp_texture_stack.show();
				}
			}
		}
		// visualize correlation results (will be for inter-pair correlation)
		if (clt_corr_combo!=null){
			if (disparity_bimap != null){
				if (!batch_mode && clt_parameters.show_map &&  (debugLevel > -2)){
					sdfa_instance.showArrays(
							disparity_bimap,
							tilesX,
							tilesY,
							true,
							name+"-DISP_MAP-D"+clt_parameters.disparity,
							ImageDtt.BIDISPARITY_TITLES);
					boolean [] trusted = 	  getTrustedDisparity(
							quadCLT_main,  // QuadCLT            quadCLT_main,  // tiles should be set
							quadCLT_aux,   // QuadCLT            quadCLT_aux,
							true,          // boolean            use_individual,
							0.14,          // double             min_combo_strength,    // check correlation strength combined for all 3 correlations
							clt_parameters.grow_disp_trust,       // double             max_trusted_disparity, // 4.0 -> change to rig_trust
							1.0,           // double             trusted_tolerance,
							null,          // boolean []         was_trusted,
							disparity_bimap ); // double [][]        bimap // current state of measurements
					for (int layer = 0; layer < disparity_bimap.length; layer ++) if (disparity_bimap[layer] != null){
						for (int nTile = 0; nTile < disparity_bimap[layer].length; nTile++) {
							if (!trusted[nTile]) disparity_bimap[layer][nTile] = Double.NaN;
						}
					}
					sdfa_instance.showArrays(
							disparity_bimap,
							tilesX,
							tilesY,
							true,
							name+"-DISP_MAP_TRUSTED-D"+clt_parameters.disparity,
							ImageDtt.BIDISPARITY_TITLES);
				}
			}

			if (!batch_mode && !infinity_corr && clt_parameters.corr_show && (debugLevel > -2)){
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

				sdfa_instance.showArrays(
						corr_rslt,
						tilesX*(2*image_dtt.transform_size),
						tilesY*(2*image_dtt.transform_size),
						true,
						name + "-CORR-D"+clt_parameters.disparity,
						titles );
			}
		}
		//		  final double [][][][][][][] clt_bidata = // new double[2][quad][nChn][tilesY][tilesX][][]; // first index - main/aux
		int quad_main = clt_bidata[0].length;
		int quad_aux =  clt_bidata[1].length;

		if (!infinity_corr && (clt_parameters.gen_chn_img || clt_parameters.gen_4_img || clt_parameters.gen_chn_stacks)) {
			ImagePlus [] imps_RGB = new ImagePlus[quad_main+quad_aux];
			for (int iQuadComb = 0; iQuadComb < imps_RGB.length; iQuadComb++){
				int iAux = (iQuadComb >= quad_main) ? 1 : 0;
				int iSubCam= iQuadComb - iAux * quad_main;

				// Uncomment to have master/aux names
				//				  String title=name+"-"+String.format("%s%02d", ((iAux>0)?"A":"M"),iSubCam);
				String title=name+"-"+String.format("%02d", iQuadComb);

				if (clt_parameters.getCorrSigma(image_dtt.isMonochrome()) > 0){ // no filter at all
					for (int chn = 0; chn < clt_bidata[iAux][iSubCam].length; chn++) {
						image_dtt.clt_lpf(
								clt_parameters.getCorrSigma(image_dtt.isMonochrome()),
								clt_bidata[iAux][iSubCam][chn],
//								image_dtt.transform_size,
								threadsMax,
								debugLevel);
					}
				}

				if (debugLevel > 0){
					System.out.println("--tilesX="+tilesX);
					System.out.println("--tilesY="+tilesY);
				}
				if (!batch_mode && (debugLevel > 0)){
					double [][] clt = new double [clt_bidata[iAux][iSubCam].length*4][];
					for (int chn = 0; chn < clt_bidata[iAux][iSubCam].length; chn++) {
						double [][] clt_set = image_dtt.clt_dbg(
								clt_bidata[iAux][iSubCam][chn],
								threadsMax,
								debugLevel);
						for (int ii = 0; ii < clt_set.length; ii++) clt[chn*4+ii] = clt_set[ii];
					}

					if (debugLevel > 0){
						sdfa_instance.showArrays(clt,
								tilesX*image_dtt.transform_size,
								tilesY*image_dtt.transform_size,
								true,
								results[iQuadComb].getTitle()+"-CLT-D"+clt_parameters.disparity);
					}
				}
				double [][] iclt_data = new double [clt_bidata[iAux][iSubCam].length][];
				for (int chn=0; chn<iclt_data.length;chn++){
					iclt_data[chn] = image_dtt.iclt_2d(
							clt_bidata[iAux][iSubCam][chn], // scanline representation of dcd data, organized as dct_size x dct_size tiles
//							image_dtt.transform_size,  // final int
							clt_parameters.clt_window,      // window_type
							15,                             // clt_parameters.iclt_mask,       //which of 4 to transform back
							0,                              // clt_parameters.dbg_mode,        //which of 4 to transform back
							threadsMax,
							debugLevel);

				}

				if (clt_parameters.gen_chn_stacks) sdfa_instance.showArrays(iclt_data,
						(tilesX + 0) * image_dtt.transform_size,
						(tilesY + 0) * image_dtt.transform_size,
						true,
						results[iQuadComb].getTitle()+"-ICLT-RGB-D"+clt_parameters.disparity);
				if (!clt_parameters.gen_chn_img) continue;

				imps_RGB[iQuadComb] = quadCLT_main.linearStackToColor( // probably no need to separate and process the second half with quadCLT_aux
						clt_parameters,
						colorProcParameters,
						rgbParameters,
						title, // String name,
						"-D"+clt_parameters.disparity, //String suffix, // such as disparity=...
						toRGB,
						!quadCLT_main.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
						!batch_mode, // true, // boolean saveShowIntermediate, // save/show if set globally
						false, // boolean saveShowFinal,        // save/show result (color image?)
						iclt_data,
						tilesX *  image_dtt.transform_size,
						tilesY *  image_dtt.transform_size,
						scaleExposures[iAux][iSubCam], // double scaleExposure, // is it needed?
						debugLevel );
			} // end of generating shifted channel images



			if (clt_parameters.gen_chn_img) {
				// combine to a sliced color image
				// assuming total number of images to be multiple of 4
				//			  int [] slice_seq = {0,1,3,2}; //clockwise
				int [] slice_seq = new int[results.length];
				for (int i = 0; i < slice_seq.length; i++) {
					slice_seq[i] = i ^ ((i >> 1) & 1); // 0,1,3,2,4,5,7,6, ...
				}
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
				quadCLT_main.eyesisCorrections.saveAndShowEnable(
						imp_stack,  // ImagePlus             imp,
						quadCLT_main.correctionsParameters, // EyesisCorrectionParameters.CorrectionParameters  correctionsParameters,
						true, // boolean               enableSave,
						!batch_mode) ;// boolean               enableShow);
			}
			if (clt_parameters.gen_4_img) {
				// Save as individual JPEG images in the model directory
				String x3d_path= quadCLT_main.correctionsParameters.selectX3dDirectory(
						name, // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
						quadCLT_main.correctionsParameters.x3dModelVersion,
						true,  // smart,
						true);  //newAllowed, // save
				for (int sub_img = 0; sub_img < imps_RGB.length; sub_img++){
					quadCLT_main.eyesisCorrections.saveAndShow(
							imps_RGB[sub_img],
							x3d_path,
							quadCLT_main.correctionsParameters.png && !clt_parameters.black_back,
							!batch_mode && clt_parameters.show_textures,
							quadCLT_main.correctionsParameters.JPEG_quality, // jpegQuality); // jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
							(debugLevel > 0) ? debugLevel : 1); // int debugLevel (print what it saves)
				}
				String model_path= quadCLT_main.correctionsParameters.selectX3dDirectory(
						name, // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
						null,
						true,  // smart,
						true);  //newAllowed, // save

				quadCLT_main.createThumbNailImage(
						imps_RGB[0],
						model_path,
						"thumb",
						debugLevel);

			}
		}

		return results;
	}

	public static void showImageFromGPU() {
		int width =  2592+8;
		int height = 1936+8;
		int l = width*height;
		String path = "/home/eyesis/workspace-python3/nvidia_dct8x8/clt/main_chn0.rbg";
		String [] titles= {"R","B","G"};
		float [] img_rbg = getFloatsFromFile(path);
		float [][] img = new float [3][l];
		for(int nc = 0; nc < 3; nc++) {
			System.arraycopy(img_rbg, l * nc, img[nc], 0 ,l);
		}

		(new ShowDoubleFloatArrays()).showArrays(
				img,
				width,
				height,
				true,
				"RBG",
				titles);

	}

	public static float [] getFloatsFromFile(String filepath) {
		float [] fdata  = null;
		try {
			FileInputStream inFile = new FileInputStream(filepath);
//			DataInputStream din = new DataInputStream(inFile);
			FileChannel inChannel = inFile.getChannel();
			int cl = (int) inChannel.size();
			ByteBuffer buffer = ByteBuffer.allocateDirect(cl); //1024*1024*60);
			buffer.order(ByteOrder.LITTLE_ENDIAN);
			buffer.clear();
			inChannel.read(buffer);
			buffer.flip();
			FloatBuffer fb = buffer.asFloatBuffer();
			fdata = new float[fb.limit()];
			fb.get(fdata);
//			fdata = fb.array();
			inFile.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return fdata;

	}

	public void saveFloatKernels(String file_prefix,
			                     double [][][][][][] clt_kernels,
			                     double [][][]       image_data,
			                     double [][][]       port_xy,
			                     boolean transpose) throws IOException {
		if (clt_kernels != null) {
			for (int chn = 0; chn < clt_kernels.length; chn++) {
				String kern_path = file_prefix+"_chn"+chn+(transpose?"_transposed":"")+".kernel";
				String offs_path = file_prefix+"_chn"+chn+(transpose?"_transposed":"")+".kernel_offsets";
				FileOutputStream fos = new FileOutputStream(kern_path);
				DataOutputStream dos = new DataOutputStream(fos);
				WritableByteChannel channel = Channels.newChannel(dos);
				int float_buffer_size = clt_kernels[chn].length * clt_kernels[chn][0].length* clt_kernels[chn][0][0].length * 4 * 64;
				ByteBuffer bb = ByteBuffer.allocate(float_buffer_size * 4);
				bb.order(ByteOrder.LITTLE_ENDIAN);
				bb.clear();
				for (int ty = 0; ty <  clt_kernels[chn][0].length; ty++) {
					for (int tx = 0; tx <  clt_kernels[chn][0][ty].length; tx++) {
						for (int col = 0; col <  clt_kernels[chn].length; col++) {
							for (int p = 0; p < 4; p++) {
								double [] pa = clt_kernels[chn][col][ty][tx][p];
								for (int i0 = 0; i0 < 64; i0++) {
									int i;
									if (transpose) {
										i = ((i0 & 7) << 3) + ((i0 >>3) & 7);
									} else {
										i = i0;
									}
//									dos.writeFloat((float)pa[i]);
									bb.putFloat((float)pa[i]);
								}
							}
						}
					}
				}
				bb.flip();
				channel.write(bb);
				dos.close();

				fos = new FileOutputStream(offs_path);
				dos = new DataOutputStream(fos);
				channel = Channels.newChannel(dos);
				float_buffer_size = clt_kernels[chn][0].length * clt_kernels[chn][0].length* clt_kernels[chn][0][0].length * 4 * clt_kernels[chn][0][0][0][4].length;
				bb = ByteBuffer.allocate(float_buffer_size * 4);
				bb.order(ByteOrder.LITTLE_ENDIAN);
				bb.clear();
				for (int ty = 0; ty <  clt_kernels[chn][0].length; ty++) {
					for (int tx = 0; tx <  clt_kernels[chn][0][ty].length; tx++) {
						for (int col = 0; col <  clt_kernels[chn].length; col++) {
							double [] pa = clt_kernels[chn][col][ty][tx][4];
							for (int i = 0; i < pa.length; i++) {
//								dos.writeFloat((float)pa[i]);
								bb.putFloat((float)pa[i]);
							}
						}
					}
				}
				bb.flip();
				channel.write(bb);
				dos.close();
			}
		}

		if (image_data != null) {
			for (int chn = 0; chn < image_data.length; chn++) {
				String img_path =  file_prefix+"_chn"+chn+".bayer";
				FileOutputStream fos = new FileOutputStream(img_path);
				DataOutputStream dos = new DataOutputStream(fos);
				WritableByteChannel channel = Channels.newChannel(dos);
				ByteBuffer bb = ByteBuffer.allocate(image_data[chn][0].length * 4);
				bb.order(ByteOrder.LITTLE_ENDIAN);
				bb.clear();
				for (int i = 0; i <  image_data[chn][0].length; i++) {
//					dos.writeFloat((float) (image_data[chn][0][i] + image_data[chn][1][i] + image_data[chn][2][i]));
					double d = 0;
					for (int c = 0; c < image_data[chn].length; c++) {
						d += image_data[chn][c][i];
					}
//					bb.putFloat((float) (image_data[chn][0][i] + image_data[chn][1][i] + image_data[chn][2][i]));
					bb.putFloat((float) d);
				}
				bb.flip();
				channel.write(bb);
				dos.close();
			}
		}
		if (port_xy != null) {
			for (int chn = 0; chn < port_xy[0].length; chn++) {
				String img_path =  file_prefix+"_chn"+chn+".portsxy";
				FileOutputStream fos = new FileOutputStream(img_path);
				DataOutputStream dos = new DataOutputStream(fos);
				WritableByteChannel channel = Channels.newChannel(dos);
				ByteBuffer bb = ByteBuffer.allocate(port_xy.length * 2 * 4);
				bb.order(ByteOrder.LITTLE_ENDIAN);
				bb.clear();
				for (int i = 0; i <  port_xy.length; i++) {
					bb.putFloat((float) (port_xy[i][chn][0])); // x-offset
					bb.putFloat((float) (port_xy[i][chn][1])); // y-offset
				}
				bb.flip();
				channel.write(bb);
				dos.close();
			}
		}
	}



	public void saveFloatKernelsBigEndian(String file_prefix,
			double [][][][][][] clt_kernels,
			double [][][]       image_data,
			double [][][]       port_xy,
			boolean transpose) throws IOException {
		if (clt_kernels != null) {
			for (int chn = 0; chn < clt_kernels.length; chn++) {
				String kern_path = file_prefix+"_chn"+chn+(transpose?"_transposed":"")+".kernel";
				String offs_path = file_prefix+"_chn"+chn+(transpose?"_transposed":"")+".kernel_offsets";
				FileOutputStream fos = new FileOutputStream(kern_path);
				DataOutputStream dos = new DataOutputStream(fos);
				for (int ty = 0; ty <  clt_kernels[chn][0].length; ty++) {
					for (int tx = 0; tx <  clt_kernels[chn][0][ty].length; tx++) {
						for (int col = 0; col <  clt_kernels[chn].length; col++) {
							for (int p = 0; p < 4; p++) {
								double [] pa = clt_kernels[chn][col][ty][tx][p];
								for (int i0 = 0; i0 < 64; i0++) {
									int i;
									if (transpose) {
										i = ((i0 & 7) << 3) + ((i0 >>3) & 7);
									} else {
										i = i0;
									}
									dos.writeFloat((float)pa[i]);
								}
							}
						}
					}
				}
				dos.close();
				fos = new FileOutputStream(offs_path);
				dos = new DataOutputStream(fos);

				for (int ty = 0; ty <  clt_kernels[chn][0].length; ty++) {
					for (int tx = 0; tx <  clt_kernels[chn][0][ty].length; tx++) {
						for (int col = 0; col <  clt_kernels[chn].length; col++) {
							double [] pa = clt_kernels[chn][col][ty][tx][4];
							for (int i = 0; i < pa.length; i++) {
								dos.writeFloat((float)pa[i]);
							}
						}
					}
				}

				dos.close();
			}
		}

		if (image_data != null) {
			for (int chn = 0; chn < image_data.length; chn++) {
				String img_path =  file_prefix+"_chn"+chn+".bayer";
				FileOutputStream fos = new FileOutputStream(img_path);
				DataOutputStream dos = new DataOutputStream(fos);
				for (int i = 0; i <  image_data[chn][0].length; i++) {
					dos.writeFloat((float) (image_data[chn][0][i] + image_data[chn][1][i] + image_data[chn][2][i]));
				}
				dos.close();
			}
		}
		if (port_xy != null) {
			for (int chn = 0; chn < port_xy[0].length; chn++) {
				String img_path =  file_prefix+"_chn"+chn+".portsxy";
				FileOutputStream fos = new FileOutputStream(img_path);
				DataOutputStream dos = new DataOutputStream(fos);
				for (int i = 0; i <  port_xy.length; i++) {
					dos.writeFloat((float) (port_xy[i][chn][0])); // x-offset
					dos.writeFloat((float) (port_xy[i][chn][1])); // y-offset
				}
				dos.close();
			}

		}

	}


	public ImagePlus [] processCLTQuadCorrPairForGPU(
			String                                         save_prefix,
			QuadCLT                                        quadCLT_main,
			QuadCLT                                        quadCLT_aux,
			ImagePlus []                                   imp_quad_main,
			ImagePlus []                                   imp_quad_aux,
			boolean [][]                                   saturation_main, // (near) saturated pixels or null
			boolean [][]                                   saturation_aux, // (near) saturated pixels or null
			CLTParameters       clt_parameters,
			EyesisCorrectionParameters.DebayerParameters   debayerParameters,
			ColorProcParameters                            colorProcParameters,
			ColorProcParameters                            colorProcParameters_aux,
			EyesisCorrectionParameters.RGBParameters       rgbParameters,
			double []	                                     scaleExposures_main, // probably not needed here - restores brightness of the final image
			double []	                                     scaleExposures_aux, // probably not needed here - restores brightness of the final image
			boolean                                        notch_mode, // use pole-detection mode for inter-camera correlation
			final int                                      lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
			final int        threadsMax,  // maximal number of threads to launch
			final boolean    updateStatus,
			final int        debugLevel){
		final boolean      batch_mode = clt_parameters.batch_run; //disable any debug images
		final boolean get_ers = !batch_mode;
		//		  boolean batch_mode = false;
		boolean infinity_corr = false;
		double [][] scaleExposures= {scaleExposures_main, scaleExposures_aux};
		boolean toRGB=     quadCLT_main.correctionsParameters.toRGB;
		ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging? - TODO - move where it belongs
		// may use this.StartTime to report intermediate steps execution times
		String name=quadCLT_main.correctionsParameters.getModelName((String) imp_quad_main[0].getProperty("name"));
		String path= (String) imp_quad_main[0].getProperty("path"); // Only for debug output
		ImagePlus [] results = new ImagePlus[imp_quad_main.length + imp_quad_main.length];
		for (int i = 0; i < results.length; i++) {
			if (i< imp_quad_main.length) {
				results[i] = imp_quad_main[i];
			} else {
				results[i] = imp_quad_main[i-imp_quad_main.length];
			}
			results[i].setTitle(results[i].getTitle()+"RAW");
		}
		if (debugLevel>1) System.out.println("processing: "+path);
/*		// 08/12/2020 Moved to conditionImageSet		
		getRigImageStacks(
				clt_parameters,  // EyesisCorrectionParameters.CLTParameters       clt_parameters,
				quadCLT_main,    // QuadCLT                                         quadCLT_main,
				quadCLT_aux,     // QuadCLT                                          quadCLT_aux,
				imp_quad_main,   // ImagePlus []                                   imp_quad_main,
				imp_quad_aux,    // ImagePlus []                                    imp_quad_aux,
				saturation_main, // boolean [][]        saturation_main, // (near) saturated pixels or null
				saturation_aux, // boolean [][]        saturation_aux, // (near) saturated pixels or null
				threadsMax,      // maximal number of threads to launch
				debugLevel);     // final int        debugLevel);
*/
		// temporary setting up tile task file (one integer per tile, bitmask
		// for testing defined for a window, later the tiles to process will be calculated based on previous passes results

		int [][]    tile_op_main = quadCLT_main.tp.setSameTileOp(clt_parameters,  clt_parameters.tile_task_op, debugLevel);
		//		  int [][]    tile_op_aux =  quadCLT_aux.tp.setSameTileOp (clt_parameters,  clt_parameters.tile_task_op, debugLevel);

		double [][] disparity_array_main = quadCLT_main.tp.setSameDisparity(clt_parameters.disparity); // [tp.tilesY][tp.tilesX] - individual per-tile expected disparity

		//TODO: Add array of default disparity - use for combining images in force disparity mode (no correlation), when disparity is predicted from other tiles

		// start with all old arrays, remove some later
		double [][][][]     clt_corr_combo =        null;
		double [][][][]     texture_tiles_main =    null; // [tp.tilesY][tp.tilesX]["RGBA".length()][]; // tiles will be 16x16, 2 visualization mode full 16 or overlapped
		double [][][][]     texture_tiles_aux =     null; // [tp.tilesY][tp.tilesX]["RGBA".length()][]; // tiles will be 16x16, 2 visualization mode full 16 or overlapped
		// undecided, so 2 modes of combining alpha - same as rgb, or use center tile only
		// Assuming same size of both image sets
		final int tilesX = quadCLT_main.tp.getTilesX();
		final int tilesY = quadCLT_main.tp.getTilesY();

		final boolean keep_lt_corr_combo = clt_parameters.corr_keep;
		// FIXME: Remove those huge arrays when not needed!
		if (clt_parameters.correlate){
			if (keep_lt_corr_combo) clt_corr_combo =    new double [ImageDtt.TCORR_TITLES.length][tilesY][tilesX][];
			texture_tiles_main =     new double [tilesY][tilesX][][]; // ["RGBA".length()][];
			texture_tiles_aux =     new double [tilesY][tilesX][][]; // ["RGBA".length()][];
			for (int i = 0; i < tilesY; i++){
				for (int j = 0; j < tilesX; j++){
					if (keep_lt_corr_combo) for (int k = 0; k<clt_corr_combo.length; k++) clt_corr_combo[k][i][j] = null;
					texture_tiles_main[i][j] = null;
					texture_tiles_aux[i][j] = null;
				}
			}

		}
		double [][] disparity_bimap  = new double [ImageDtt.BIDISPARITY_TITLES.length][]; //[0] -residual disparity, [1] - orthogonal (just for debugging) last 4 - max pixel differences
		ImageDtt image_dtt = new ImageDtt(
				clt_parameters.transform_size,
				quadCLT_main.isMonochrome(),
				quadCLT_main.isLwir(),
				clt_parameters.getScaleStrength(false));
		double [][] ml_data = null;
		double [][][]       ers_delay = get_ers?(new double [2][][]):null;
		// here all data is ready (images, kernels) to try GPU code
		float [][] main_bayer = new float [quadCLT_main.image_data.length][quadCLT_main.image_data[0][0].length];
		float [][] dst_bayer =  new float [quadCLT_main.image_data.length][quadCLT_main.image_data[0][0].length];
		for (int nc = 0; nc < main_bayer.length; nc++) {
			int nc1 = (nc +1) % 4;
			for (int i = 0; i < main_bayer[nc].length; i++) {
				main_bayer[nc1][i] = (float) (quadCLT_main.image_data[nc][0][i] + quadCLT_main.image_data[nc][1][i] + quadCLT_main.image_data[nc][2][i]);
				dst_bayer[nc][i]= nc*main_bayer[nc].length + i;
			}
		}

		double [][][]       port_xy_main_dbg = new double [tilesX*tilesY][][];
		double [][][]       port_xy_aux_dbg =  new double [tilesX*tilesY][][];
//		double [][][]       corr2ddata =       new double [1][][];
		double [][] disparity_map = new double [ImageDtt.DISPARITY_TITLES.length][];

		final double [][][][][][][] clt_bidata = // new double[2][quad][nChn][tilesY][tilesX][][]; // first index - main/aux
				image_dtt.clt_bi_quad_dbg (
						clt_parameters,                       // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
						clt_parameters.getFatZero(image_dtt.isMonochrome()),              // final double              fatzero,         // May use correlation fat zero from 2 different parameters - fat_zero and rig.ml_fatzero
						notch_mode,                           //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
						lt_rad,                               // final int                                      lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
						// first measurement - use default setting
						clt_parameters.rig.no_int_x0, // boolean   no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
						tile_op_main,                         // final int [][]            tile_op_main,    // [tilesY][tilesX] - what to do - 0 - nothing for this tile
						disparity_array_main,                 // final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
						quadCLT_main.image_data,              // final double [][][]       image_data_main, // first index - number of image in a quad
						quadCLT_aux.image_data,               // final double [][][]       image_data_aux,  // first index - number of image in a quad
						quadCLT_main.saturation_imp,          // final boolean [][]        saturation_main, // (near) saturated pixels or null
						quadCLT_aux.saturation_imp,           // final boolean [][]        saturation_aux,  // (near) saturated pixels or null
						// correlation results - combo will be for the correation between two quad cameras
						clt_corr_combo,                       // final double [][][][]     clt_corr_combo,  // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
						// [type][tilesY][tilesX] should be set by caller
						// types: 0 - selected correlation (product+offset), 1 - sum
						disparity_bimap,                      // final double [][]    disparity_bimap, // [23][tilesY][tilesX]
						ml_data,                              // 	final double [][]         ml_data,         // data for ML - 10 layers - 4 center areas (3x3, 5x5,..) per camera-per direction, 1 - composite, and 1 with just 1 data (target disparity)
						disparity_map,                        // final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
						texture_tiles_main,                   // final double [][][][]     texture_tiles_main, // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
						texture_tiles_aux,                    // final double [][][][]     texture_tiles_aux,  // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
						imp_quad_main[0].getWidth(),          // final int                 width,
						quadCLT_main.getGeometryCorrection(), // final GeometryCorrection  geometryCorrection_main,
						quadCLT_aux.getGeometryCorrection(),  // final GeometryCorrection  geometryCorrection_aux,
						quadCLT_main.getCLTKernels(),         // final double [][][][][][] clt_kernels_main, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
						quadCLT_aux.getCLTKernels(),          // final double [][][][][][] clt_kernels_aux,  // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
						clt_parameters.corr_magic_scale,      // final double              corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
						true,                                 // 	final boolean             keep_clt_data,
						ers_delay,                            // final double [][][]       ers_delay,        // if not null - fill with tile center acquisition delay
						threadsMax,                           // final int                 threadsMax,  // maximal number of threads to launch
						debugLevel,                           // final int                 globalDebugLevel);
						port_xy_main_dbg,                     // final double [][][]       port_xy_main_dbg, // for each tile/port save x,y pixel coordinates (gpu code development)
						port_xy_aux_dbg);                      // final double [][][]       port_xy_aux_dbg) // for each tile/port save x,y pixel coordinates (gpu code development)

/////		  double [][] disparity_map = new double [ImageDtt.DISPARITY_TITLES.length][]; //[0] -residual disparity, [1] - orthogonal (just for debugging)

		String [] sub_titles = new String [GPUTileProcessor.NUM_CAMS * (GPUTileProcessor.NUM_COLORS+1)];
		double [][] sub_disparity_map = new double [sub_titles.length][];
		for (int ncam = 0; ncam < GPUTileProcessor.NUM_CAMS; ncam++) {
			sub_disparity_map[ncam] = disparity_map[ncam + ImageDtt.IMG_DIFF0_INDEX];
			sub_titles[ncam] = ImageDtt.DISPARITY_TITLES[ncam + ImageDtt.IMG_DIFF0_INDEX];
			for (int ncol = 0; ncol < GPUTileProcessor.NUM_COLORS; ncol++) {
				sub_disparity_map[ncam + (ncol + 1)* GPUTileProcessor.NUM_CAMS] =
						disparity_map[ncam +ncol* GPUTileProcessor.NUM_CAMS+ ImageDtt.IMG_TONE_RGB];
				sub_titles[ncam + (ncol + 1)* GPUTileProcessor.NUM_CAMS] =
						ImageDtt.DISPARITY_TITLES[ncam +ncol* GPUTileProcessor.NUM_CAMS+ ImageDtt.IMG_TONE_RGB];
			}
		}
//		String [] sub_titles = {ImageDtt.DISPARITY_TITLES[ImageDtt.IMG_DIFF0_INDEX]

		(new ShowDoubleFloatArrays()).showArrays(
				sub_disparity_map,
	    		tilesX,
	    		tilesY,
				true,
				name + "-CPU-EXTRA-D"+clt_parameters.disparity,
				sub_titles);
		// Create list of all correlation pairs
		double [][][][][][] clt_data = clt_bidata[0];
		int numTiles = tilesX * tilesY;
		int numPairs = GPUTileProcessor.NUM_PAIRS;
		int [] corr_indices = new int [numTiles * numPairs];
		int indx=0;
		for (int i = 0; i < numTiles; i++) {
			for (int j = 0; j < numPairs; j++) {
				corr_indices[indx++] = (i << GPUTileProcessor.CORR_NTILE_SHIFT) + j;
			}
		}
		double [][] corrs2d = image_dtt.get2DCorrs(
				clt_parameters,  // final CLTParameters       clt_parameters,
				clt_data,        // final double [][][][][][] clt_data, // [channel_in_quad][color][tileY][tileX][band][pixel];
				corr_indices,      // final int    []           pairs_list,
				threadsMax,      // final int                 threadsMax,  // maximal number of threads to launch
				debugLevel);     // final int                 debugLevel
		float [][] fcorrs2d = new float [corrs2d.length][corrs2d[0].length];
		// for compatibility with the actual GPUI output
		for (int n = 0; n < corrs2d.length; n++) {
			for (int i = 0; i < corrs2d[0].length; i++) {
				fcorrs2d[n][i] = (float) corrs2d[n][i];
			}
		}
		int [] wh = new int[2];
		double [][] dbg_corr = GPUTileProcessor.getCorr2DView(
	    		tilesX,
	    		tilesY,
	    		corr_indices,
	    		fcorrs2d,
	    		wh);
		(new ShowDoubleFloatArrays()).showArrays(
				dbg_corr,
				wh[0],
				wh[1],
				true,
				name + "-CPU-CORR2D-D"+clt_parameters.disparity,
				GPUTileProcessor.getCorrTitles());

		if ((save_prefix != null) && (save_prefix != "")) {

			if (debugLevel < -1000) {
				return null;
			}

			String kernel_dir = save_prefix+"clt/";
			File kdir = new File(kernel_dir);
			kdir.mkdir();
			//		boolean [][] what_to_save = {{false,false,true}, {false,false,true}};
			boolean [][] what_to_save = {{true,true,true}, {true,true,true}};
			try {
				saveFloatKernels(
						kernel_dir +"main", // String file_prefix,
//						(what_to_save[0][0]?clt_kernels_main:null), // double [][][][][][] clt_kernels, // null
						(what_to_save[0][0]?quadCLT_main.getCLTKernels():null), // double [][][][][][] clt_kernels, // null
						(what_to_save[0][1]?quadCLT_main.image_data:null),
						(what_to_save[0][2]?port_xy_main_dbg:null), // double [][][]       port_xy,
						true);
			} catch (IOException e) {
				System.out.println("Failed to save flattened kernels tp "+kernel_dir);
				// TODO Auto-generated catch block
				e.printStackTrace();
			} // boolean transpose);

			try {
				saveFloatKernels(
						kernel_dir +"aux", // String file_prefix,
//						(what_to_save[1][0]?clt_kernels_aux:null), // double [][][][][][] clt_kernels, // null
						(what_to_save[1][0]?quadCLT_aux.getCLTKernels():null), // double [][][][][][] clt_kernels, // null
						(what_to_save[1][1]?quadCLT_aux.image_data:null),
						(what_to_save[1][2]?port_xy_aux_dbg:null), // double [][][]       port_xy,
						true);
			} catch (IOException e) {
				System.out.println("Failed to save flattened kernels tp "+kernel_dir);
				e.printStackTrace();
			} // boolean transpose);

			try {
				quadCLT_main.getGeometryCorrection().saveFloatsGPU(kernel_dir +"main");
			} catch (IOException e) {
				System.out.println("Failed to save geometry correction data (float) to "+kernel_dir);
				e.printStackTrace();
			}

			try {
				quadCLT_main.getGeometryCorrection().saveDoublesGPU(kernel_dir +"main");
			} catch (IOException e) {
				System.out.println("Failed to save geometry correction data (double) to "+kernel_dir);
				e.printStackTrace();
			}

			quadCLT_main.getGeometryCorrection().getCorrVector().getRotMatricesDbg();
			quadCLT_main.getGeometryCorrection().getCorrVector().getRotDeriveMatricesDbg();

			if (debugLevel < -1000) {
				return null;
			}
			if (ers_delay !=null) {
				showERSDelay(ers_delay);
			}

		}
		double [][] texture_nonoverlap_main = null;
		double [][] texture_nonoverlap_aux = null;
		double [][] texture_overlap_main = null;
		double [][] texture_overlap_aux = null;
		String [] rgba_titles = {"red","blue","green","alpha"};
		String [] rgba_weights_titles = {"red","blue","green","alpha","port0","port1","port2","port3","r-rms","b-rms","g-rms","w-rms"};
		if ((texture_tiles_main != null) && (texture_tiles_aux != null)){
			if ((debugLevel > -1) && (clt_parameters.tileX >= 0) && (clt_parameters.tileY >= 0) && (clt_parameters.tileX < tilesX) && (clt_parameters.tileY < tilesY)) {
				double [][] texture_tile = texture_tiles_main[clt_parameters.tileY][clt_parameters.tileX];
				int tile = +clt_parameters.tileY * tilesX  +clt_parameters.tileX;
    			System.out.println("=== tileX= "+clt_parameters.tileX+" tileY= "+clt_parameters.tileY+" tile="+tile+" ===");

    			for (int slice =0; slice < texture_tile.length; slice++) {
    				System.out.println("\n=== Slice="+slice+" ===");
    				for (int i = 0; i < 2 * GPUTileProcessor.DTT_SIZE; i++) {
    					for (int j = 0; j < 2 * GPUTileProcessor.DTT_SIZE; j++) {
    						System.out.print(String.format("%10.4f ",
    								texture_tile[slice][2 * GPUTileProcessor.DTT_SIZE * i + j]));
    					}
    					System.out.println();
    				}
    			}
    			(new ShowDoubleFloatArrays()).showArrays(
						texture_tile,
						2 * image_dtt.transform_size,
						2 * image_dtt.transform_size,
						true,
						name + "-TXTNOL-CPU-D"+clt_parameters.disparity+"-X"+clt_parameters.tileX+"-Y"+clt_parameters.tileY,
						(clt_parameters.keep_weights?rgba_weights_titles:rgba_titles));

			}

			if (clt_parameters.show_nonoverlap){
				texture_nonoverlap_main = image_dtt.combineRBGATiles(
						texture_tiles_main,                 // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}
//						image_dtt.transform_size,
						false,                         // when false - output each tile as 16x16, true - overlap to make 8x8
						clt_parameters.sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
						threadsMax,                    // maximal number of threads to launch
						debugLevel);
				sdfa_instance.showArrays(
						texture_nonoverlap_main,
						tilesX * (2 * image_dtt.transform_size),
						tilesY * (2 * image_dtt.transform_size),
						true,
						name + "-TXTNOL-D"+clt_parameters.disparity+"-MAIN",
						(clt_parameters.keep_weights?rgba_weights_titles:rgba_titles));

				texture_nonoverlap_aux = image_dtt.combineRBGATiles(
						texture_tiles_aux,                 // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}
//						image_dtt.transform_size,
						false,                         // when false - output each tile as 16x16, true - overlap to make 8x8
						clt_parameters.sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
						threadsMax,                    // maximal number of threads to launch
						debugLevel);
				sdfa_instance.showArrays(
						texture_nonoverlap_aux,
						tilesX * (2 * image_dtt.transform_size),
						tilesY * (2 * image_dtt.transform_size),
						true,
						name + "-TXTNOL-D"+clt_parameters.disparity+"-AUX",
						(clt_parameters.keep_weights?rgba_weights_titles:rgba_titles));
			}
			if (!infinity_corr && (clt_parameters.show_overlap || clt_parameters.show_rgba_color)){
				int alpha_index = 3;
				texture_overlap_main = image_dtt.combineRBGATiles(
						texture_tiles_main,                 // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}
//						image_dtt.transform_size,
						true,                         // when false - output each tile as 16x16, true - overlap to make 8x8
						clt_parameters.sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
						threadsMax,                    // maximal number of threads to launch
						debugLevel);
				if (clt_parameters.alpha1 > 0){ // negative or 0 - keep alpha as it was
					double scale = (clt_parameters.alpha1 > clt_parameters.alpha0) ? (1.0/(clt_parameters.alpha1 - clt_parameters.alpha0)) : 0.0;
					for (int i = 0; i < texture_overlap_main[alpha_index].length; i++){
						double d = texture_overlap_main[alpha_index][i];
						if      (d >=clt_parameters.alpha1) d = 1.0;
						else if (d <=clt_parameters.alpha0) d = 0.0;
						else d = scale * (d- clt_parameters.alpha0);
						texture_overlap_main[alpha_index][i] = d;
					}
				}

				texture_overlap_aux = image_dtt.combineRBGATiles(
						texture_tiles_aux,                 // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}
//						image_dtt.transform_size,
						true,                         // when false - output each tile as 16x16, true - overlap to make 8x8
						clt_parameters.sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
						threadsMax,                    // maximal number of threads to launch
						debugLevel);
				if (clt_parameters.alpha1 > 0){ // negative or 0 - keep alpha as it was
					double scale = (clt_parameters.alpha1 > clt_parameters.alpha0) ? (1.0/(clt_parameters.alpha1 - clt_parameters.alpha0)) : 0.0;
					for (int i = 0; i < texture_overlap_aux[alpha_index].length; i++){
						double d = texture_overlap_aux[alpha_index][i];
						if      (d >=clt_parameters.alpha1) d = 1.0;
						else if (d <=clt_parameters.alpha0) d = 0.0;
						else d = scale * (d- clt_parameters.alpha0);
						texture_overlap_aux[alpha_index][i] = d;
					}
				}

				if (!batch_mode && clt_parameters.show_overlap) {
					sdfa_instance.showArrays(
							texture_overlap_main,
							tilesX * image_dtt.transform_size,
							tilesY * image_dtt.transform_size,
							true,
							name + "-TXTOL-D"+clt_parameters.disparity+"-MAIN",
							(clt_parameters.keep_weights?rgba_weights_titles:rgba_titles));
				}
				if (!batch_mode && clt_parameters.show_overlap) {
					sdfa_instance.showArrays(
							texture_overlap_aux,
							tilesX * image_dtt.transform_size,
							tilesY * image_dtt.transform_size,
							true,
							name + "-TXTOL-D"+clt_parameters.disparity+"-AUX",
							(clt_parameters.keep_weights?rgba_weights_titles:rgba_titles));
				}


				if (!batch_mode && clt_parameters.show_rgba_color) {
					// for now - use just RGB. Later add option for RGBA
					double [][] texture_rgb_main = {texture_overlap_main[0],texture_overlap_main[1],texture_overlap_main[2]};
					double [][] texture_rgba_main = {texture_overlap_main[0],texture_overlap_main[1],texture_overlap_main[2],texture_overlap_main[3]};
					double [][] texture_rgb_aux = {texture_overlap_aux[0],texture_overlap_aux[1],texture_overlap_aux[2]};
					double [][] texture_rgba_aux = {texture_overlap_aux[0],texture_overlap_aux[1],texture_overlap_aux[2],texture_overlap_aux[3]};
					ImagePlus imp_texture_main = quadCLT_main.linearStackToColor(
							clt_parameters,
							colorProcParameters,
							rgbParameters,
							name+"-texture", // String name,
							"-D"+clt_parameters.disparity+"-MAINCPU", //String suffix, // such as disparity=...
							toRGB,
							!quadCLT_main.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
							false, // true, // boolean saveShowIntermediate, // save/show if set globally
							false, // true, // boolean saveShowFinal,        // save/show result (color image?)
							((clt_parameters.alpha1 > 0)? texture_rgba_main: texture_rgb_main),
							tilesX *  image_dtt.transform_size,
							tilesY *  image_dtt.transform_size,
							1.0,         // double scaleExposure, // is it needed?
							debugLevel );
					ImagePlus imp_texture_aux = quadCLT_aux.linearStackToColor(
							clt_parameters,
							colorProcParameters_aux,
							rgbParameters,
							name+"-texture", // String name,
							"-D"+clt_parameters.disparity+"-AUX", //String suffix, // such as disparity=...
							toRGB,
							!quadCLT_aux.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
							false, // true, // boolean saveShowIntermediate, // save/show if set globally
							false, // true, // boolean saveShowFinal,        // save/show result (color image?)
							((clt_parameters.alpha1 > 0)? texture_rgba_aux: texture_rgb_aux),
							tilesX *  image_dtt.transform_size,
							tilesY *  image_dtt.transform_size,
							1.0,         // double scaleExposure, // is it needed?
							debugLevel );
					int width = imp_texture_main.getWidth();
					int height =imp_texture_main.getHeight();
					ImageStack texture_stack=new ImageStack(width,height);
					texture_stack.addSlice("main",      imp_texture_main.getProcessor().getPixels());
					texture_stack.addSlice("auxiliary", imp_texture_aux. getProcessor().getPixels());
					ImagePlus imp_texture_stack = new ImagePlus(name+"-TEXTURES-D"+clt_parameters.disparity, texture_stack);
					imp_texture_stack.getProcessor().resetMinAndMax();
					//					  imp_texture_stack.updateAndDraw();
					imp_texture_stack.show();
				}
			}
		}
		// visualize correlation results (will be for inter-pair correlation)
		if (clt_corr_combo!=null){
			if (disparity_bimap != null){
				if (!batch_mode && clt_parameters.show_map &&  (debugLevel > -2)){
					sdfa_instance.showArrays(
							disparity_bimap,
							tilesX,
							tilesY,
							true,
							name+"-DISP_MAP-D"+clt_parameters.disparity,
							ImageDtt.BIDISPARITY_TITLES);
					boolean [] trusted = 	  getTrustedDisparity(
							quadCLT_main,  // QuadCLT            quadCLT_main,  // tiles should be set
							quadCLT_aux,   // QuadCLT            quadCLT_aux,
							true,          // boolean            use_individual,
							0.14,          // double             min_combo_strength,    // check correlation strength combined for all 3 correlations
							clt_parameters.grow_disp_trust,       // double             max_trusted_disparity, // 4.0 -> change to rig_trust
							1.0,           // double             trusted_tolerance,
							null,          // boolean []         was_trusted,
							disparity_bimap ); // double [][]        bimap // current state of measurements
					for (int layer = 0; layer < disparity_bimap.length; layer ++) if (disparity_bimap[layer] != null){
						for (int nTile = 0; nTile < disparity_bimap[layer].length; nTile++) {
							if (!trusted[nTile]) disparity_bimap[layer][nTile] = Double.NaN;
						}
					}
					sdfa_instance.showArrays(
							disparity_bimap,
							tilesX,
							tilesY,
							true,
							name+"-DISP_MAP_TRUSTED-D"+clt_parameters.disparity,
							ImageDtt.BIDISPARITY_TITLES);
				}
			}

			if (!batch_mode && !infinity_corr && clt_parameters.corr_show && (debugLevel > -2)){
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

				sdfa_instance.showArrays(
						corr_rslt,
						tilesX*(2*image_dtt.transform_size),
						tilesY*(2*image_dtt.transform_size),
						true,
						name + "-CORR-D"+clt_parameters.disparity,
						titles );
			}
		}
		//		  final double [][][][][][][] clt_bidata = // new double[2][quad][nChn][tilesY][tilesX][][]; // first index - main/aux
		int quad_main = clt_bidata[0].length;
		int quad_aux =  clt_bidata[1].length;

		if (!infinity_corr && (clt_parameters.gen_chn_img || clt_parameters.gen_4_img || clt_parameters.gen_chn_stacks)) {
			ImagePlus [] imps_RGB = new ImagePlus[quad_main+quad_aux];
			for (int iQuadComb = 0; iQuadComb < imps_RGB.length; iQuadComb++){
				int iAux = (iQuadComb >= quad_main) ? 1 : 0;
				int iSubCam= iQuadComb - iAux * quad_main;

				// Uncomment to have master/aux names
				//				  String title=name+"-"+String.format("%s%02d", ((iAux>0)?"A":"M"),iSubCam);
				String title=name+"-"+String.format("%02d", iQuadComb);
				if (clt_parameters.getCorrSigma(image_dtt.isMonochrome()) > 0){ // no filter at all
					for (int chn = 0; chn < clt_bidata[iAux][iSubCam].length; chn++) {
		                int debug_lpf = ((iQuadComb ==0) && (chn==0))?3: debugLevel;
						image_dtt.clt_lpf(
								clt_parameters.getCorrSigma(image_dtt.isMonochrome()),
								clt_bidata[iAux][iSubCam][chn],
//								image_dtt.transform_size,
								threadsMax,
								debug_lpf);
					}
				}

				if (debugLevel > 0){
					System.out.println("--tilesX="+tilesX);
					System.out.println("--tilesY="+tilesY);
				}
				if (!batch_mode && (debugLevel > 0)){
					double [][] clt = new double [clt_bidata[iAux][iSubCam].length*4][];
					for (int chn = 0; chn < clt_bidata[iAux][iSubCam].length; chn++) {
						double [][] clt_set = image_dtt.clt_dbg(
								clt_bidata[iAux][iSubCam][chn],
								threadsMax,
								debugLevel);
						for (int ii = 0; ii < clt_set.length; ii++) clt[chn*4+ii] = clt_set[ii];
					}

					if (debugLevel > 0){
						sdfa_instance.showArrays(clt,
								tilesX*image_dtt.transform_size,
								tilesY*image_dtt.transform_size,
								true,
								results[iQuadComb].getTitle()+"-CLT-D"+clt_parameters.disparity);
					}
				}
				double [][] iclt_data = new double [clt_bidata[iAux][iSubCam].length][];
				for (int chn=0; chn<iclt_data.length;chn++){
					iclt_data[chn] = image_dtt.iclt_2d_debug_gpu(
							clt_bidata[iAux][iSubCam][chn], // scanline representation of dcd data, organized as dct_size x dct_size tiles
//							image_dtt.transform_size,  // final int
							clt_parameters.clt_window,      // window_type
							15,                             // clt_parameters.iclt_mask,       //which of 4 to transform back
							0,                              // clt_parameters.dbg_mode,        //which of 4 to transform back
							threadsMax,
							debugLevel,
							clt_parameters.tileX,           // final int                 debug_tileX
							clt_parameters.tileY);          // final int                 debug_tileY

				}

				if (clt_parameters.gen_chn_stacks) sdfa_instance.showArrays(iclt_data,
						(tilesX + 0) * image_dtt.transform_size,
						(tilesY + 0) * image_dtt.transform_size,
						true,
						results[iQuadComb].getTitle()+"-ICLT-RGB-D"+clt_parameters.disparity);
				if (!clt_parameters.gen_chn_img) continue;

				imps_RGB[iQuadComb] = quadCLT_main.linearStackToColor( // probably no need to separate and process the second half with quadCLT_aux
						clt_parameters,
						colorProcParameters,
						rgbParameters,
						title, // String name,
						"-D"+clt_parameters.disparity, //String suffix, // such as disparity=...
						toRGB,
						!quadCLT_main.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
						!batch_mode, // true, // boolean saveShowIntermediate, // save/show if set globally
						false, // boolean saveShowFinal,        // save/show result (color image?)
						iclt_data,
						tilesX *  image_dtt.transform_size,
						tilesY *  image_dtt.transform_size,
						scaleExposures[iAux][iSubCam], // double scaleExposure, // is it needed?
						debugLevel );
			} // end of generating shifted channel images



			if (clt_parameters.gen_chn_img) {
				// combine to a sliced color image
				// assuming total number of images to be multiple of 4
				//			  int [] slice_seq = {0,1,3,2}; //clockwise
				int [] slice_seq = new int[results.length];
				for (int i = 0; i < slice_seq.length; i++) {
					slice_seq[i] = i ^ ((i >> 1) & 1); // 0,1,3,2,4,5,7,6, ...
				}
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
				quadCLT_main.eyesisCorrections.saveAndShowEnable(
						imp_stack,  // ImagePlus             imp,
						quadCLT_main.correctionsParameters, // EyesisCorrectionParameters.CorrectionParameters  correctionsParameters,
						true, // boolean               enableSave,
						!batch_mode) ;// boolean               enableShow);
			}
			if (clt_parameters.gen_4_img) {
				// Save as individual JPEG images in the model directory
				String x3d_path= quadCLT_main.correctionsParameters.selectX3dDirectory(
						name, // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
						quadCLT_main.correctionsParameters.x3dModelVersion,
						true,  // smart,
						true);  //newAllowed, // save
				for (int sub_img = 0; sub_img < imps_RGB.length; sub_img++){
					quadCLT_main.eyesisCorrections.saveAndShow(
							imps_RGB[sub_img],
							x3d_path,
							quadCLT_main.correctionsParameters.png && !clt_parameters.black_back,
							!batch_mode && clt_parameters.show_textures,
							quadCLT_main.correctionsParameters.JPEG_quality, // jpegQuality); // jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
							(debugLevel > 0) ? debugLevel : 1); // int debugLevel (print what it saves)
				}
				String model_path= quadCLT_main.correctionsParameters.selectX3dDirectory(
						name, // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
						null,
						true,  // smart,
						true);  //newAllowed, // save

				quadCLT_main.createThumbNailImage(
						imps_RGB[0],
						model_path,
						"thumb",
						debugLevel);

			}
		}

		return results;
	}
	
	public void showERSDelay(double [][][] ers_delay)
	{
		int tilesX = quadCLT_main.tp.getTilesX();
		int tilesY = quadCLT_main.tp.getTilesY();
		int nTiles = tilesX * tilesY;
		int nchn_main = ers_delay[0].length;
		int nchn_aux =  ers_delay[1].length;
		double [][] dbg_img = new double[nchn_main + nchn_aux][tilesX*tilesY];
		for (int nTile = 0; nTile< nTiles; nTile++) {
			double avg_dly = 0;
			for (int i = 0; i < nchn_main; i++) {
				avg_dly += ers_delay[0][i][nTile];
			}
			avg_dly /= nchn_main;
			for (int i = 0; i < nchn_main; i++) {
				dbg_img[i][nTile] = ers_delay[0][i][nTile] - avg_dly;
			}
			for (int i = 0; i < nchn_aux; i++) {
				dbg_img[i + nchn_main][nTile] = ers_delay[1][i][nTile] - avg_dly;
			}
		}
		(new ShowDoubleFloatArrays()).showArrays(
				dbg_img,
				tilesX,
				tilesY,
				true,
				"ERS_DELAYS");
	}


	//	  public double [][][][] getRigImageStacks(
// 08/12/2020 Moved to conditionImageSet
/*	
	public void getRigImageStacks( // Removed all references
			CLTParameters       clt_parameters,
			QuadCLT                                        quadCLT_main,
			QuadCLT                                        quadCLT_aux,
			ImagePlus []                                   imp_quad_main,
			ImagePlus []                                   imp_quad_aux,
			boolean [][]                                   saturation_main, // (near) saturated pixels or null
			boolean [][]                                   saturation_aux, // (near) saturated pixels or null
			int                                            threadsMax,  // maximal number of threads to launch
			int                                            debugLevel){
		double [][][] double_stacks_main = new double [imp_quad_main.length][][];
		for (int i = 0; i < double_stacks_main.length; i++){
			double_stacks_main[i] = quadCLT_main.eyesisCorrections.bayerToDoubleStack(
					imp_quad_main[i], // source Bayer image, linearized, 32-bit (float))
					  null, // no margins, no oversample
					  quadCLT_main.is_mono);
		}

		double [][][] double_stacks_aux = new double [imp_quad_aux.length][][];
		for (int i = 0; i < double_stacks_aux.length; i++){
			double_stacks_aux[i] = quadCLT_aux.eyesisCorrections.bayerToDoubleStack(
					imp_quad_aux[i], // source Bayer image, linearized, 32-bit (float))
					  null, // no margins, no oversample
					  quadCLT_aux.is_mono);
		}

		for (int i = 0; i < double_stacks_main.length; i++){
			for (int j =0 ; j < double_stacks_main[i][0].length; j++){
				double_stacks_main[i][2][j]*=0.5; // Scale green 0.5 to compensate more pixels than R,B
			}
		}

		for (int i = 0; i < double_stacks_aux.length; i++){
			if (double_stacks_aux[i].length > 2) { // skip for monochrome, only if color
				for (int j =0 ; j < double_stacks_aux[i][0].length; j++){
					double_stacks_aux[i][2][j]*=0.5; // Scale green 0.5 to compensate more pixels than R,B
				}
			}
		}
		quadCLT_main.setTiles (imp_quad_main[0], // set global tp.tilesX, tp.tilesY
				clt_parameters,
				threadsMax);
		quadCLT_aux.setTiles (imp_quad_aux[0], // set global tp.tilesX, tp.tilesY
				clt_parameters,
				threadsMax);
		String name_main = (String) imp_quad_main[0].getProperty("name");
		String name_aux =  (String) imp_quad_main[0].getProperty("name");
		quadCLT_main.image_name =     name_main;
		quadCLT_aux.image_name =      name_aux;
		quadCLT_main.image_data =     double_stacks_main;
		quadCLT_aux.image_data =      double_stacks_aux;
		quadCLT_main.saturation_imp = saturation_main;
		quadCLT_aux.saturation_imp =  saturation_aux;
		//		  quadCLT_main.tp.resetCLTPasses();
		quadCLT_main.tp.setTrustedCorrelation(clt_parameters.grow_disp_trust);
		quadCLT_aux.tp.setTrustedCorrelation(clt_parameters.grow_disp_trust);
	}
*/

	public void processInfinityRigs( // actually there is no sense to process multiple image sets. Combine with other processing?
			QuadCLT quadCLT_main,
			QuadCLT quadCLT_aux,
			CLTParameters       clt_parameters,
			ColorProcParameters                            colorProcParameters, //
			ColorProcParameters                            colorProcParameters_aux, //
			final int                                      threadsMax,  // maximal number of threads to launch
			final boolean                                  updateStatus,
			final int                                      debugLevel) throws Exception
	{

		this.startTime=System.nanoTime();
		String [] sourceFiles=quadCLT_main.correctionsParameters.getSourcePaths();
		QuadCLT.SetChannels [] set_channels_main = quadCLT_main.setChannels(debugLevel);
		QuadCLT.SetChannels [] set_channels_aux =  quadCLT_aux.setChannels(debugLevel);
		if ((set_channels_main == null) || (set_channels_main.length==0) || (set_channels_aux == null) || (set_channels_aux.length==0)) {
			System.out.println("No files to process (of "+sourceFiles.length+")");
			return;
		}
		double [] referenceExposures_main = quadCLT_main.eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		double [] referenceExposures_aux =  quadCLT_aux.eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		for (int nSet = 0; nSet < set_channels_main.length; nSet++){
			// check it is the same set for both cameras
			if (set_channels_aux.length <= nSet ) {
				throw new Exception ("Set names for cameras do not match: main camera: '"+set_channels_main[nSet].name()+"', aux. camera: nothing");
			}
			if (!set_channels_main[nSet].name().equals(set_channels_aux[nSet].name())) {
				throw new Exception ("Set names for cameras do not match: main camera: '"+set_channels_main[nSet].name()+"', aux. camera: '"+set_channels_main[nSet].name()+"'");
			}

			int [] channelFiles_main = set_channels_main[nSet].fileNumber();
			int [] channelFiles_aux =  set_channels_aux[nSet].fileNumber();
			boolean [][] saturation_imp_main = (clt_parameters.sat_level > 0.0)? new boolean[channelFiles_main.length][] : null;
			boolean [][] saturation_imp_aux =  (clt_parameters.sat_level > 0.0)? new boolean[channelFiles_main.length][] : null;
			double [] scaleExposures_main = new double[channelFiles_main.length];
			double [] scaleExposures_aux =  new double[channelFiles_main.length];

			ImagePlus [] imp_srcs_main = quadCLT_main.conditionImageSet(
					clt_parameters,                 // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					colorProcParameters,            //  ColorProcParameters                       colorProcParameters, //
					sourceFiles,                    // String []                                 sourceFiles,
					set_channels_main[nSet].name(), // String                                    set_name,
					referenceExposures_main,        // double []                                 referenceExposures,
					channelFiles_main,              // int []                                    channelFiles,
					scaleExposures_main,            //output  // double [] scaleExposures
					saturation_imp_main,            //output  // boolean [][]                              saturation_imp,
					threadsMax,                 // int                                       threadsMax,
					debugLevel); // int                                       debugLevel);

			ImagePlus [] imp_srcs_aux = quadCLT_aux.conditionImageSet(
					clt_parameters,                 // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					colorProcParameters_aux,        //  ColorProcParameters                       colorProcParameters, //
					sourceFiles,                    // String []                                 sourceFiles,
					set_channels_aux[nSet].name(), // String                                    set_name,
					referenceExposures_aux,        // double []                                 referenceExposures,
					channelFiles_aux,              // int []                                    channelFiles,
					scaleExposures_aux,            //output  // double [] scaleExposures
					saturation_imp_aux,            //output  // boolean [][]                              saturation_imp,
					threadsMax,                 // int                                       threadsMax,
					debugLevel); // int                                       debugLevel);

			// Temporary processing individually with the old code
			processInfinityRig(
					quadCLT_main,               // QuadCLT                                        quadCLT_main,
					quadCLT_aux,                // QuadCLT                                        quadCLT_aux,
					imp_srcs_main,              // ImagePlus []                                   imp_quad_main,
					imp_srcs_aux,               // ImagePlus []                                   imp_quad_aux,
					saturation_imp_main,        // boolean [][]                                   saturation_main, // (near) saturated pixels or null
					saturation_imp_aux,         // boolean [][]                                   saturation_aux, // (near) saturated pixels or null
					clt_parameters,             // EyesisCorrectionParameters.CLTParameters       clt_parameters,
					threadsMax,                 // final int        threadsMax,  // maximal number of threads to launch
					updateStatus,               // final boolean    updateStatus,
					debugLevel);                // final int        debugLevel);


			Runtime.getRuntime().gc();
			if (debugLevel >-1) System.out.println("Processing set "+(nSet+1)+" (of "+set_channels_aux.length+") finished at "+
					IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");

			if (quadCLT_aux.eyesisCorrections.stopRequested.get()>0) {
				System.out.println("User requested stop");
				System.out.println("Processing "+(nSet + 1)+" file sets (of "+set_channels_main.length+") finished at "+
						IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
				return;
			}
		}
		System.out.println("processInfinityRigs(): processing "+(quadCLT_main.getTotalFiles(set_channels_main)+quadCLT_aux.getTotalFiles(set_channels_aux))+" files ("+set_channels_main.length+" file sets) finished at "+
				IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");

	}


	// use macro mode (strengths with limited disparity) to align at infinity if it is way off.
	public boolean prealignInfinityRig(
			QuadCLT                                        quadCLT_main,
			QuadCLT                                        quadCLT_aux,
			CLTParameters       clt_parameters,
			final int        threadsMax,  // maximal number of threads to launch
			final boolean    updateStatus,
			final int        debugLevel){
		double [][] disparity_bimap = measureNewRigDisparity(
				quadCLT_main,      // QuadCLT             quadCLT_main,    // tiles should be set
				quadCLT_aux,       // QuadCLT             quadCLT_aux,
				null,              // double [][]         src_bimap,       // current state of measurements (or null for new measurement)
				0.0,               // double              disparity,
				null,              // ArrayList<Integer>  tile_list,       // or null. If non-null - do not remeasure members of the list
				clt_parameters,    // EyesisCorrectionParameters.CLTParameters       clt_parameters,
				false,             //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
				0,                              // final int  // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
				// first measurement - use default setting
				clt_parameters.rig.no_int_x0, // boolean   no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
				threadsMax,        // final int           threadsMax,      // maximal number of threads to launch
				updateStatus,      // final boolean       updateStatus,
				debugLevel);       // final int           debugLevel);
		// create two arrays of main and aux strengths and small disparity
		//	public double  inf_max_disp_main =         0.15;
		//	public double  inf_max_disp_aux =          0.15;
		final int tilesX = quadCLT_main.tp.getTilesX();
		final int tilesY = quadCLT_main.tp.getTilesY();
		final int tile_op_all = clt_parameters.tile_task_op; //FIXME Use some constant?
		double max_main_disparity = 3* clt_parameters.rig.inf_max_disp_main;
		double max_aux_disparity = 3* clt_parameters.rig.inf_max_disp_aux;
		double min_strength_main =0.2;
		double min_strength_aux =0.2;

		int numTiles = disparity_bimap[ImageDtt.BI_STR_FULL_INDEX].length;
		double [][][] macro_pair = new double[2][1][numTiles]; // second index for compatibility with multiple colors
		for (int nTile = 0; nTile < numTiles; nTile++) {
			double d_main= disparity_bimap[ImageDtt.BI_DISP_FULL_INDEX ][nTile];
			double d_aux=  disparity_bimap[ImageDtt.BI_ADISP_FULL_INDEX][nTile];
			double s_main= disparity_bimap[ImageDtt.BI_STR_FULL_INDEX ][nTile];
			double s_aux=  disparity_bimap[ImageDtt.BI_ASTR_FULL_INDEX][nTile];
			if ((Math.abs(d_main) <= max_main_disparity) && (s_main > min_strength_main)) macro_pair[0][0][nTile] = s_main-min_strength_main;
			if ((Math.abs(d_aux)  <= max_aux_disparity)  && (s_aux  > min_strength_aux))  macro_pair[1][0][nTile] = s_aux- min_strength_aux;
			//			  macro_pair[0][0][nTile] = s_main;
			//			  macro_pair[1][0][nTile] = s_aux;
		}
		if (debugLevel >-2) {
			double [][] dbg_macro = {macro_pair[0][0],macro_pair[1][0]};
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_macro,
					tilesX,
					tilesY,
					true,
					"MACRO-INPUT");
		}

		ImageDtt image_dtt = new ImageDtt(
				clt_parameters.transform_size,
				quadCLT_main.isMonochrome(),
				quadCLT_main.isLwir(),
				clt_parameters.getScaleStrength(false));

		int macro_scale = image_dtt.transform_size;
		int mTilesX = tilesX/macro_scale;
		int mTilesY = tilesY/macro_scale;
		int [][] mtile_op = new int [mTilesY][mTilesX];
		for (int mTileY=0; mTileY < mTilesY; mTileY++) {
			for (int mTileX=0; mTileX < mTilesX; mTileX++) {
				mtile_op[mTileY][mTileX] = tile_op_all;
			}
		}
		double [][] mdisparity_array = new double [mTilesY][mTilesX]; // keep all zeros
		double [][] mdisparity_bimap = new double [ImageDtt.BIDISPARITY_TITLES.length][];
		image_dtt.clt_bi_macro(
				clt_parameters,                       // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
				clt_parameters.getFatZero(image_dtt.isMonochrome()),              // final double              fatzero,         // May use correlation fat zero from 2 different parameters - fat_zero and rig.ml_fatzero
				macro_scale, // final int                 macro_scale,
				mtile_op, // final int [][]            tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
				mdisparity_array, // final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
				macro_pair, // 	final double [][][]       image_rig_data,  // [2][1][pixels] (single color channel)
				mdisparity_bimap, // final double [][]         disparity_bimap, // [23][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
				tilesX, // 	final int                 width,
				quadCLT_main.getGeometryCorrection(), // final GeometryCorrection  geometryCorrection_main,
				quadCLT_aux.getGeometryCorrection(),  // final GeometryCorrection  geometryCorrection_aux,
				clt_parameters.corr_magic_scale,      //  final double              corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
				threadsMax,        // final int           threadsMax,      // maximal number of threads to launch
				debugLevel + 1);       // final int           debugLevel);
		// Display macro correlation results (later convert to a full size disparity_bimap for compatibility with infinity alignment
		if (debugLevel > -2) {
			(new ShowDoubleFloatArrays()).showArrays(
					mdisparity_bimap,
					mTilesX,
					mTilesY,
					true,
					"MACRO-DISP_MAP",
					ImageDtt.BIDISPARITY_TITLES);
		}

		/*
		 * ImageDtt.BIDISPARITY_TITLES.length
	public void  clt_bi_macro(
			final EyesisCorrectionParameters.CLTParameters       clt_parameters,
			final int                 macro_scale,
			final int [][]            tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
			final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
			final double [][][]       image_rig_data,  // [2][1][pixels] (single color channel)
			final double [][]         disparity_bimap, // [23][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
			final int                 width,
			final GeometryCorrection  geometryCorrection_main,
			final GeometryCorrection  geometryCorrection_aux, // it has rig offset)
			final double              corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
			final int                 threadsMax,  // maximal number of threads to launch
			final int                 debugLevel)

		 */
		return true;
	}

	public boolean processInfinityRig(
			QuadCLT                                        quadCLT_main,
			QuadCLT                                        quadCLT_aux,
			ImagePlus []                                   imp_quad_main,
			ImagePlus []                                   imp_quad_aux,
			boolean [][]                                   saturation_main, // (near) saturated pixels or null
			boolean [][]                                   saturation_aux, // (near) saturated pixels or null
			CLTParameters       clt_parameters,
			final int        threadsMax,  // maximal number of threads to launch
			final boolean    updateStatus,
			final int        debugLevel0){
		final int        debugLevel = debugLevel0 + (clt_parameters.rig.rig_mode_debug?2:0);
		//		  double [][][][] double_stacks =
/*		// 08/12/2020 Moved to condifuinImageSet				  		
		getRigImageStacks(
				clt_parameters,  // EyesisCorrectionParameters.CLTParameters       clt_parameters,
				quadCLT_main,    // QuadCLT                                         quadCLT_main,
				quadCLT_aux,     // QuadCLT                                          quadCLT_aux,
				imp_quad_main,   // ImagePlus []                                   imp_quad_main,
				imp_quad_aux,    // ImagePlus []                                    imp_quad_aux,
				saturation_main, // boolean [][]        saturation_main, // (near) saturated pixels or null
				saturation_aux, // boolean [][]        saturation_aux, // (near) saturated pixels or null
				threadsMax,      // maximal number of threads to launch
				debugLevel);     // final int        debugLevel);

		quadCLT_main.tp.resetCLTPasses();
		quadCLT_aux.tp.resetCLTPasses();
*/		
		/*
FIXME - make it work
 		  prealignInfinityRig(
 				  quadCLT_main, // QuadCLT                                        quadCLT_main,
 				  quadCLT_aux, // QuadCLT                                        quadCLT_aux,
 				  clt_parameters,    // EyesisCorrectionParameters.CLTParameters       clt_parameters,
 				  threadsMax,      // maximal number of threads to launch
 				  updateStatus,    // final boolean    updateStatus,
 				  debugLevel);     // final int        debugLevel);
if (debugLevel > -100) return true; // temporarily !
		 */
		double disparity_correction_main = 0.0; // actual disparity in main camera pixels for "infinity" objects
		if (clt_parameters.infinity_distace_map.containsKey(quadCLT_main.image_name)){
			double infinity_distance = clt_parameters.infinity_distace_map.get(quadCLT_main.image_name);
			disparity_correction_main = quadCLT_main.geometryCorrection.getDisparityFromZ(infinity_distance);
			if (debugLevel > -5) {
				System.out.println("Found infinity distance record for "+quadCLT_main.image_name+", it is "+infinity_distance+
						" m, corresponding to "+disparity_correction_main+" pixels of the main camera.");
			}
		} else {
			if (debugLevel > -5) {
				System.out.println("No infinity distance record for "+quadCLT_main.image_name+" found, considering far objects to be true infinity");
			}
		}
		if (debugLevel > -5) {
			System.out.println("disparity_correction_main= "+disparity_correction_main);
		}

		final int tilesX = quadCLT_main.tp.getTilesX();

		// perform full re-measure cycles
		double [][] disparity_bimap = null;
		int [] num_new = new int[1];
		for (int num_full_cycle = 0; num_full_cycle < clt_parameters.rig.rig_adjust_full_cycles;num_full_cycle++) {
			if (debugLevel > -5) {
				System.out.println("processInfinityRig(), reselecting and re-measuring tiles (may increase RMS), full_cycle "+
						(num_full_cycle + 1)+" of "+clt_parameters.rig.rig_adjust_full_cycles);
			}
			disparity_bimap = null;
			ArrayList<Integer> tile_list = new ArrayList<Integer>();
			// measure and refine
			for (int disp_step = 0; disp_step < clt_parameters.rig.rig_num_disp_steps; disp_step++) {
				double disparity = disparity_correction_main;
				if (disp_step > 0) {
					disparity += disp_step * clt_parameters.rig.rig_disp_range/(clt_parameters.rig.rig_num_disp_steps -1);
				}
				disparity_bimap = measureNewRigDisparity(
						quadCLT_main,      // QuadCLT             quadCLT_main,    // tiles should be set
						quadCLT_aux,       // QuadCLT             quadCLT_aux,
						disparity_bimap,   // double [][]         src_bimap,       // current state of measurements (or null for new measurement)
						disparity,         // double              disparity,
						tile_list,         // ArrayList<Integer>  tile_list,       // or null. If non-null - do not remeasure members of the list
						clt_parameters,    // EyesisCorrectionParameters.CLTParameters       clt_parameters,
						false,             //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
						0,                 // final int  // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
						// first measurement - use default setting
						clt_parameters.rig.no_int_x0, // boolean   no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
						threadsMax,        // final int           threadsMax,      // maximal number of threads to launch
						updateStatus,      // final boolean       updateStatus,
						debugLevel);       // final int           debugLevel);

				if (disparity_bimap != null){
					if (clt_parameters.show_map &&  (debugLevel > -2) && clt_parameters.rig.rig_mode_debug){
						(new ShowDoubleFloatArrays()).showArrays(
								disparity_bimap,
								tilesX,
								disparity_bimap[0].length/tilesX,
								true,
								"DISP_MAP-D"+disparity,
								ImageDtt.BIDISPARITY_TITLES);
					}
				}
				if (disp_step > 0) { // refine non-infinity passes
					// first - refine once for the main camera
					double [][] prev_bimap = null;
					double [] scale_bad = null;
					for (int nrefine = 0; nrefine < clt_parameters.rig.num_refine_master; nrefine++) {
						double [][] disparity_bimap_new =  refineRig(
								quadCLT_main,    // QuadCLT                        quadCLT_main,    // tiles should be set
								quadCLT_aux,     // QuadCLT                        quadCLT_aux,
								disparity_bimap, // double [][]                    src_bimap,       // current state of measurements (or null for new measurement)
								prev_bimap,      // double [][]                    prev_bimap, // previous state of measurements or null
								null, // scale_bad,       // double []                      scale_bad,
								2, //0,               // int                            refine_mode,     // 0 - by main, 1 - by aux, 2 - by inter
								true,            // boolean                        keep_inf,        // keep expected disparity 0.0 if it was so
								disparity_correction_main, // double                                   inf_disparity,
								clt_parameters.rig.refine_min_strength , // double refine_min_strength, // do not refine weaker tiles
								clt_parameters.rig.refine_tolerance ,    // double refine_tolerance,    // do not refine if absolute disparity below
								null, // tile_list,       // ArrayList<Integer>             tile_list,       // or null
								num_new,         // int     []                                     num_new,
								clt_parameters,  // EyesisCorrectionParameters.CLTParameters clt_parameters,
								false,           //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
								0,                              // final int  // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
								// refine mode - no window offset
								true,            // boolean                                  no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
								threadsMax,      // final int                      threadsMax,      // maximal number of threads to launch
								updateStatus,    // final boolean                  updateStatus,
								debugLevel);     // final int                      debugLevel);
						prev_bimap = disparity_bimap;
						disparity_bimap = disparity_bimap_new;
						if (disparity_bimap != null){
							if (clt_parameters.show_map &&  (debugLevel > 2) && clt_parameters.rig.rig_mode_debug){
								(new ShowDoubleFloatArrays()).showArrays(
										disparity_bimap,
										tilesX,
										disparity_bimap[0].length/tilesX,
										true,
										"DISP_MAP-refine_master-"+nrefine,
										ImageDtt.BIDISPARITY_TITLES);
							}
						}

					}
					// now refine for the inter-camera correlation
					prev_bimap = null;
					for (int nrefine = 0; nrefine < clt_parameters.rig.num_refine_inter; nrefine++) {
						double [][] disparity_bimap_new =  refineRig(
								quadCLT_main,    // QuadCLT                        quadCLT_main,    // tiles should be set
								quadCLT_aux,     // QuadCLT                        quadCLT_aux,
								disparity_bimap, // double [][]                    src_bimap,       // current state of measurements (or null for new measurement)
								prev_bimap,      // double [][]                    prev_bimap, // previous state of measurements or null
								scale_bad,       // double []                      scale_bad,
								2,               // int                            refine_mode,     // 0 - by main, 1 - by aux, 2 - by inter
								true,            // boolean                        keep_inf,        // keep expected disparity 0.0 if it was so
								disparity_correction_main, // double                                   inf_disparity,
								clt_parameters.rig.refine_min_strength , // double refine_min_strength, // do not refine weaker tiles
								clt_parameters.rig.refine_tolerance ,    // double refine_tolerance,    // do not refine if absolute disparity below
								null, // tile_list,       // ArrayList<Integer>             tile_list,       // or null
								num_new,         // int     []                                     num_new,
								clt_parameters,  // EyesisCorrectionParameters.CLTParameters clt_parameters,
								false,             //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
								0,                              // final int  // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
								// refine mode - no window offset
								true,            // boolean                                  no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
								threadsMax,      // final int                      threadsMax,      // maximal number of threads to launch
								updateStatus,    // final boolean                  updateStatus,
								debugLevel);     // final int                      debugLevel);
						prev_bimap = disparity_bimap;
						disparity_bimap = disparity_bimap_new;
						if (disparity_bimap != null){
							if (clt_parameters.show_map &&  (debugLevel > 2) && clt_parameters.rig.rig_mode_debug){
								(new ShowDoubleFloatArrays()).showArrays(
										disparity_bimap,
										tilesX,
										disparity_bimap[0].length/tilesX,
										true,
										"DISP_MAP-refine_inter-"+nrefine,
										ImageDtt.BIDISPARITY_TITLES);
							}
						}

					}
				} // if (disparity > 0.0) { // refine non-infinity passes
				// rebuild tiles list that match requirements, debug-show results
				tile_list = selectRigTiles(
						clt_parameters, // EyesisCorrectionParameters.CLTParameters       clt_parameters,
						true, // boolean select_infinity,
						true, // boolean select_noninfinity,
						disparity_correction_main, // double                                   inf_disparity,
						disparity_bimap, // double[][] disparity_bimap,
						tilesX);  // int tilesX,
				if (debugLevel > -1) showListedRigTiles(
						"rig_tiles_D"+disparity, // title
						tile_list, //ArrayList<Integer> tile_list,
						disparity_bimap, // double[][] disparity_bimap,
						tilesX); // int tilesX

			} //  for (int disp_step = 0; disp_step < clt_parameters.rig.rig_num_disp_steps; disp_step++)

			// short cycle (remeasure only for the list). May also break from it if RMS is not improving
			for (int num_short_cycle = 0; num_short_cycle < clt_parameters.rig.rig_adjust_short_cycles;num_short_cycle++) {
				// refine for the existing list - all listed tiles, no thresholds
				disparity_bimap =  refineRig(
						quadCLT_main,    // QuadCLT                        quadCLT_main,    // tiles should be set
						quadCLT_aux,     // QuadCLT                        quadCLT_aux,
						disparity_bimap, // double [][]                    src_bimap,       // current state of measurements (or null for new measurement)
						null,            // double [][]                                    prev_bimap, // previous state of measurements or null
						null,       // double []                           scale_bad,
						2,               // int                            refine_mode,     // 0 - by main, 1 - by aux, 2 - by inter
						// will still re-measure infinity if refine_min_strength == 0.0
						true,            // boolean                        keep_inf,        // keep expected disparity 0.0 if it was so
						disparity_correction_main, // double                                   inf_disparity,
						0.0,             // double refine_min_strength, // do not refine weaker tiles
						0.0,             // double refine_tolerance,    // do not refine if absolute disparity below
						tile_list,       // ArrayList<Integer>             tile_list,       // or null
						num_new,         // int     []                                     num_new,
						clt_parameters,  // EyesisCorrectionParameters.CLTParameters clt_parameters,
						false,             //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
						0,                              // final int  // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
						// refine mode - disable window offset
						true,            // boolean                                  no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
						threadsMax,      // final int                      threadsMax,      // maximal number of threads to launch
						updateStatus,    // final boolean                  updateStatus,
						debugLevel);     // final int                      debugLevel);
				// show updated results for the list
				if (debugLevel > -2) showListedRigTiles(
						"selected_rig_tiles", // title
						tile_list, //ArrayList<Integer> tile_list,
						disparity_bimap, // double[][] disparity_bimap,
						tilesX); // int tilesX

				// do actual adjustment step, update rig parameters
				quadCLT_aux.geometryCorrection.getRigCorrection(
						clt_parameters.rig.inf_weight ,     // double             infinity_importance, // of all measurements
						clt_parameters.rig.inf_weight_disp, // double             dx_max, //  = 0.3;
						clt_parameters.rig.inf_weight_disp_pow,                                // double             dx_pow, //  = 1.0;
						clt_parameters.rig.rig_adjust_orientation,        // boolean            adjust_orientation,
						clt_parameters.rig.rig_adjust_roll,               // boolean            adjust_roll,
						clt_parameters.rig.rig_adjust_zoom,               // boolean            adjust_zoom,
						clt_parameters.rig.rig_adjust_angle,              // boolean            adjust_angle,
						clt_parameters.rig.rig_adjust_distance,           // boolean            adjust_distance,
						clt_parameters.rig.rig_adjust_forward,            // boolean            adjust_forward, // not used
						clt_parameters.rig.rig_correction_scale,          // double             scale_correction,
						disparity_correction_main, // double             infinity_disparity,
						tile_list,                                        // ArrayList<Integer> tile_list,
						quadCLT_main,                                     // QuadCLT            qc_main,
						disparity_bimap[ImageDtt.BI_STR_CROSS_INDEX],     // double []          strength,
						disparity_bimap[ImageDtt.BI_DISP_CROSS_DX_INDEX], // double []          diff_x, // used only with target_disparity == 0
						disparity_bimap[ImageDtt.BI_DISP_CROSS_DY_INDEX], // double []          diff_y,
						disparity_bimap[ImageDtt.BI_TARGET_INDEX],        // double []          target_disparity,
						debugLevel+1);

			} // end of for (int num_short_cycle = 0; num_short_cycle < clt_parameters.rig.rig_adjust_short_cycles;num_short_cycle++) {


		} // end of for (int num_full_cycle = 0; num_full_cycle < clt_parameters.rig.rig_adjust_full_cycles;num_full_cycle++) {
		if (disparity_bimap != null){
			if (clt_parameters.show_map &&  (debugLevel > 0) && clt_parameters.rig.rig_mode_debug){
				(new ShowDoubleFloatArrays()).showArrays(
						disparity_bimap,
						tilesX,
						disparity_bimap[0].length/tilesX,
						true,
						"DISP_MAP",
						ImageDtt.BIDISPARITY_TITLES);
			}
		}

		// loop here with the same tileList

		return true;
	}



	public boolean processPoles(
			QuadCLT            quadCLT_main,  // tiles should be set
			QuadCLT            quadCLT_aux,
			BiCamDSI                                       biCamDSI,
			CLTParameters       clt_parameters,
			final int                                      threadsMax,  // maximal number of threads to launch
			final boolean                                  updateStatus,
			final int                                      globalDebugLevel)
	{
		// -1 and above - follow clt_parameters.poles.poles_debug_level, -2 and below - follow global debug level
		final int debugLevel = (globalDebugLevel > -2) ? clt_parameters.poles.poles_debug_level:globalDebugLevel;
		if ((quadCLT_main == null) ||
				(quadCLT_main.tp == null) ||
				(quadCLT_main.tp.clt_3d_passes == null) ||
				(biCamDSI_persistent== null) ||
				(biCamDSI_persistent.biScans== null) ||
				(quadCLT_main.tp.rig_pre_poles_ds == null)
				) {
			String msg = "Data is not available. Please run \"Ground truth\" first";
			IJ.showMessage("Error",msg);
			System.out.println(msg);
			return false;
		}
		int scan_index = biCamDSI_persistent.biScans.size()-1;
		BiScan biScan =  biCamDSI_persistent.biScans.get(scan_index);
		if (debugLevel> -1) {
			biScan.showScan(quadCLT_main.image_name+"LastBiScan-"+scan_index, null);
		}
		//		  if (poleProcessor_persistent == null) {
		poleProcessor_persistent = new PoleProcessor(
				biCamDSI_persistent,
				this, //	TwoQuadCLT twoQuadCLT,
				quadCLT_main.tp.getTilesX(),
				quadCLT_main.tp.getTilesY());
		//		  }
		PoleProcessor pp = poleProcessor_persistent;
		boolean [] selection =  quadCLT_main.tp.rig_pre_poles_sel.clone();
		int num_bugs_corrected= pp.zero_tiles_check(
				"pre-poles",
				quadCLT_main.tp.rig_pre_poles_ds, // double [][] ds,
				true); // boolean debug)
		if (num_bugs_corrected > 0) {
			System.out.println("Corrected "+num_bugs_corrected+" bugs, where disparity is <= 0 with non-zero strength before poles");
		}

		double [][] poles_ds = pp.processPoles(
				quadCLT_main,  // QuadCLT            quadCLT_main,  // tiles should be set
				quadCLT_aux,   //QuadCLT            quadCLT_aux,
				quadCLT_main.tp.rig_pre_poles_ds, //  double [][]                                    src_ds, // source disparity, strength pair
				selection, // boolean []                                     selection, // source tile selection, will be modified
				clt_parameters, // EyesisCorrectionParameters.CLTParameters       clt_parameters,
				threadsMax,  // final int                                      threadsMax,  // maximal number of threads to launch
				updateStatus, // final boolean                                  updateStatus,
				debugLevel); // final int                                      globalDebugLevel)

		num_bugs_corrected= pp.zero_tiles_check(
				"post-poles",
				poles_ds, // double [][] ds,
				true); // boolean debug)
		if (num_bugs_corrected > 0) {
			System.out.println("Corrected "+num_bugs_corrected+" bugs, where disparity is <= 0 with non-zero strength after poles");
		}
		quadCLT_main.tp.rig_post_poles_ds =        poles_ds; // Rig disparity and strength before processing poles
		quadCLT_main.tp.rig_post_poles_sel =       selection; // Rig tile selection before processing poles

		CLTPass3d scan_last = quadCLT_main.tp.clt_3d_passes.get( quadCLT_main.tp.clt_3d_passes.size() -1); // get really last one

		quadCLT_main.tp.trimCLTPasses(false); // remove rig composite scan if any
		CLTPass3d rig_scan = quadCLT_main.tp.compositeScan(
				poles_ds[0],  // final double []             disparity,
				poles_ds[1],  // final double []             strength,
				selection,                  // final boolean []            selected,
				debugLevel);                // final int                   debugLevel)
		rig_scan.texture_tiles = scan_last.texture_tiles;
		// scan_last
		quadCLT_main.tp.clt_3d_passes.add(rig_scan);
		quadCLT_main.tp.saveCLTPasses(true);       // rig pass
		return true;
	}

	//  improve DSI acquired for a single camera by use of a pair
	// Run this after "CLT 3D"

	public void groundTruth(
			QuadCLT            quadCLT_main,  // tiles should be set
			QuadCLT            quadCLT_aux,
			CLTParameters       clt_parameters,
			ColorProcParameters                            colorProcParameters, //
			ColorProcParameters                            colorProcParameters_aux, //
			final int                                      threadsMax,  // maximal number of threads to launch
			final boolean                                  updateStatus,
			final int                                      debugLevel)// throws Exception
	{
		if ((quadCLT_main.tp == null) || (quadCLT_main.tp.clt_3d_passes == null)) {
			String msg = "DSI data not available. Please run\"CLT 3D\" first";
			IJ.showMessage("ERROR",msg);
			System.out.println(msg);
			return;
		}

		double [][] rig_disparity_strength =  quadCLT_main.getGroundTruthByRig(); // saved pair

		if (rig_disparity_strength == null) {
			rig_disparity_strength = groundTruthByRigPlanes(
					quadCLT_main,   // QuadCLT            quadCLT_main,  // tiles should be set
					quadCLT_aux,    // QuadCLT            quadCLT_aux,
					clt_parameters, // EyesisCorrectionParameters.CLTParameters       clt_parameters,
					colorProcParameters,          //  ColorProcParameters                       colorProcParameters, //
					colorProcParameters_aux,          //  ColorProcParameters                       colorProcParameters, //
					threadsMax,     // final int                                      threadsMax,  // maximal number of threads to launch
					updateStatus,   // final boolean                                  updateStatus,
					debugLevel);    // final int                                      debugLevel);
			if (rig_disparity_strength != null) {
				quadCLT_main.tp.rig_disparity_strength = rig_disparity_strength;
			}
		}
		CLTPass3d scan_bg =   quadCLT_main.tp.clt_3d_passes.get( 0); // get bg scan
		//		  quadCLT_main.tp.trimCLTPasses(false); // remove rig composite scan if any
		// last but not including any rid data
		CLTPass3d scan_last = quadCLT_main.tp.clt_3d_passes.get( quadCLT_main.tp.clt_3d_passes_size -1); // get last one

		// TODO: combine in a single function to always call after groundTruthByRig. Or before use?
		boolean [] infinity_select = scan_bg.getSelected(); // null;
		boolean [] was_select = scan_last.getSelected(); // null;
		boolean[] selection = null;
		final int tilesX = quadCLT_main.tp.getTilesX();
		final int tilesY = quadCLT_main.tp.getTilesY();
		BiCamDSI biCamDSI = new BiCamDSI( tilesX, tilesY, threadsMax);
		//          boolean [][] dbg_sel = (debugLevel > -4)? new boolean [8][]:null;// was 4
		boolean [][] dbg_sel = (debugLevel > -2)? new boolean [8][]:null;// was 4
		if (dbg_sel!=null) dbg_sel[0] = infinity_select;
		if (dbg_sel!=null) dbg_sel[1] = was_select;
		if (clt_parameters.rig.ltfar_en) {
			double strength_floor = (clt_parameters.rig.ltfar_auto_floor)? Double.NaN: (clt_parameters.rig.lt_trusted_strength*clt_parameters.rig.lt_strength_rfloor);
			double ltfar_min_strength = clt_parameters.rig.ltfar_min_rstrength * clt_parameters.rig.lt_trusted_strength;

			selection = biCamDSI.selectFarObjects(
					strength_floor, //double     strength_floor,
					clt_parameters.rig.ltfar_min_disparity , //	double     min_disparity,
					clt_parameters.rig.ltfar_min_mean , //	double     min_mean,
					clt_parameters.rig.ltfar_max_disparity , //double     max_disparity,
					clt_parameters.rig.ltfar_max_mean , //double     max_mean,
					clt_parameters.rig.ltfar_min_disp_to_rms , //double     min_disp_to_rms,
					ltfar_min_strength , //double     min_strength,
					clt_parameters.rig.ltfar_neib_dist , //int        neib_dist, // >=1
					clt_parameters.rig.ltfar_rsigma , //double     rsigma,
					clt_parameters.rig.ltfar_frac, // double     tile_frac,
					//ltfar_frac
					infinity_select, // null, // TODO: add real	after debug boolean [] infinity_select,
					rig_disparity_strength[0], // 	double []  disparity,
					rig_disparity_strength[1], // 	double []  strength,
					tilesX, // 	int        tilesX, debug only
					debugLevel); // int        debugLevel);
			if (dbg_sel!=null) dbg_sel[2] = selection.clone();
		}

		selection=biCamDSI.selectNearObjects(
				0.0, // double     min_strength,
				infinity_select, // boolean [] infinity_select,
				selection, // boolean [] selection,
				rig_disparity_strength[0], // 	double []  disparity,
				rig_disparity_strength[1]);  // 	double []  strength,



		if (dbg_sel!=null) dbg_sel[3] = selection.clone();

		if (clt_parameters.rig.rf_master_infinity) {
			selection = biCamDSI.combineSelections(
					selection, // boolean [] selection1,
					was_select, // boolean [] selection2,
					infinity_select, // boolean [] enable,
					false); // boolean    invert_enable)
		}
		if (clt_parameters.rig.rf_master_near) {
			selection = biCamDSI.combineSelections(
					selection, // boolean [] selection1,
					was_select, // boolean [] selection2,
					infinity_select, // boolean [] enable,
					true); // boolean    invert_enable)
		}
		if (dbg_sel!=null) dbg_sel[4] = selection.clone(); // far+near+old

		boolean [] selection_lone = biCamDSI.selectLoneFar(
				clt_parameters.rig.ltfar_trusted_s, //  double     min_far_strength,
				clt_parameters.rig.ltfar_trusted_d, // 	double     min_far_disparity,
				infinity_select,                    // boolean [] infinity_select,
				null,                               // boolean [] selection,
				rig_disparity_strength[0],          // double []  disparity,
				rig_disparity_strength[1]);         // double []  strength)

		selection_lone = biCamDSI.selectNearOverInfinity(
				clt_parameters.rig.ltgrp_min_strength,       // double     min_strength,
				clt_parameters.rig.ltgrp_min_neibs,          // int        min_neibs,
				clt_parameters.rig.ltgrp_gr_min_disparity,   // double     min_disparity,
				clt_parameters.rig.ltgrp_gr_disp_atolerance, // double     disp_atolerance,
				clt_parameters.rig.ltgrp_gr_disp_rtolerance, // double     disp_rtolerance,
				infinity_select,                             // boolean [] infinity_select,
				selection,                                   // boolean [] selection,
				rig_disparity_strength[0],                   // double []  disparity,
				rig_disparity_strength[1]);                  // double []  strength,

		// add this selection, but keep it too till after filtering
		selection = biCamDSI.combineSelections(
				selection,                          // boolean [] selection1,
				selection_lone,                     // boolean [] selection2,
				null,                               // boolean [] enable,
				false);                             // boolean    invert_enable)
		if (dbg_sel!=null) dbg_sel[5] = selection.clone(); // far+near+old+lone

		// filter only near objects
		biCamDSI.expandShrinkExpandFilter(
				clt_parameters.rig.rf_pre_expand_near,   // int        pre_expand,
				clt_parameters.rig.rf_shrink_near,       // 	int        shrink,
				clt_parameters.rig.rf_post_expand_near,  // int        post_expand,
				selection,                          // 	boolean [] selection,
				infinity_select); // 	boolean [] prohibit)
		// filter all - infinity and near (less aggressive)

		if (dbg_sel!=null) dbg_sel[6] = selection.clone(); // far+near+old+lone- filtered near

		biCamDSI.expandShrinkExpandFilter(
				clt_parameters.rig.rf_pre_expand,   // int        pre_expand,
				clt_parameters.rig.rf_shrink,       // 	int        shrink,
				clt_parameters.rig.rf_post_expand,  // int        post_expand,
				selection,                          // 	boolean [] selection,
				null); // 	boolean [] prohibit)

		if (dbg_sel!=null) dbg_sel[7] = selection.clone(); // far+near+old+lone- filtered near-filtered far

		// restore lone tile selections (may be removed by filter)
		selection = biCamDSI.combineSelections(
				selection,                          // boolean [] selection1,
				selection_lone,                     // boolean [] selection2,
				null,                               // boolean [] enable,
				false);                             // boolean    invert_enable)
		// restore master camera selections (poles may be removed by a filter
		if (clt_parameters.rig.rf_master_infinity) {
			selection = biCamDSI.combineSelections(
					selection, // boolean [] selection1,
					was_select, // boolean [] selection2,
					infinity_select, // boolean [] enable,
					false); // boolean    invert_enable)
		}
		if (clt_parameters.rig.rf_master_near) {
			selection = biCamDSI.combineSelections(
					selection, // boolean [] selection1,
					was_select, // boolean [] selection2,
					infinity_select, // boolean [] enable,
					true); // boolean    invert_enable)
		}

		if (dbg_sel != null) {
			double [][] dbg_img = new double[7][selection.length];
			for (int nTile = 0; nTile < selection.length;nTile++) {
				dbg_img[0][nTile] = (dbg_sel[0][nTile]? 1.0:0.0) + (dbg_sel[1][nTile]?2.0:0.0);
				dbg_img[1][nTile] = (dbg_sel[2][nTile]? 1.0:0.0) + (dbg_sel[3][nTile]?2.0:0.0);
				dbg_img[2][nTile] = (dbg_sel[4][nTile]? 1.0:0.0) + (dbg_sel[5][nTile]?2.0:0.0);
				dbg_img[3][nTile] = (dbg_sel[6][nTile]? 1.0:0.0);
				dbg_img[4][nTile] = (dbg_sel[7][nTile]? 1.0:0.0);
				dbg_img[5][nTile] = (selection[nTile]? 1.0:0.0);
				dbg_img[6][nTile] = (selection_lone[nTile]? 1.0:0.0);
			}
			String [] titles = {"old_sel","new_sel","combo-lone", "f-near", "f-far","final","lone"};
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					tilesX,
					dbg_img[0].length/tilesX,
					true,
					"Selections",
					titles);
		}
		// set composite scan

		boolean [] cond_sel = (clt_parameters.rig.rf_remove_unselected)? selection:null;
		double [][] f_rig_disparity_strength =  biCamDSI.filterDisparityStrength(
				clt_parameters.rig.rf_min_disp,     // double  min_disparity,
				rig_disparity_strength[0],          // double []  disparity,
				rig_disparity_strength[1],          // double []  strength)
				cond_sel);                          // 	boolean [] selected)
		// CLT ASSIGN needs best texture for each tile. Initially will just copy from the previous master
		// composite scan, later - fill disparity gaps and re-measure

		/*
		 *  TODO: interpolate disparities before measuring to fill gaps?
	  public double [][][][][] getRigTextures(
			  boolean                                  need_master,
			  boolean                                  need_aux,
			  double []                                disparity, // non-nan - measure
			  QuadCLT                                  quadCLT_main,  // tiles should be set
			  QuadCLT                                  quadCLT_aux,
			  EyesisCorrectionParameters.CLTParameters clt_parameters,
			  final int                                threadsMax,  // maximal number of threads to launch
			  final boolean                            updateStatus,
			  final int                                debugLevel) // throws Exception

		 */

		quadCLT_main.tp.trimCLTPasses(false); // remove rig composite scan if any
		//		 quadCLT_main.tp.rig_pre_poles_ds =        f_rig_disparity_strength; // Rig disparity and strength before processing poles
		//		 quadCLT_main.tp.rig_pre_poles_sel =        selection; // Rig tile selection before processing poles
		quadCLT_main.tp.rig_pre_poles_ds =        rig_disparity_strength; // Rig disparity and strength before processing poles
		quadCLT_main.tp.rig_pre_poles_sel =        selection; // Rig tile selection before processing poles
		quadCLT_main.tp.rig_post_poles_ds =        null;
		quadCLT_main.tp.rig_post_poles_sel =       null;

		if (debugLevel > -4) {
			System.out.println("Saved quadCLT_main.tp.rig_pre_poles_ds and quadCLT_main.tp.rig_pre_poles_sel");
		}
		// process poles if enabled
		if (clt_parameters.poles.poles_en) {
			if (debugLevel > -4) {
				System.out.println("Processing poles");
			}

			//			  if (poleProcessor_persistent == null) {
			poleProcessor_persistent = new PoleProcessor(
					biCamDSI_persistent,
					this, //	TwoQuadCLT twoQuadCLT,
					quadCLT_main.tp.getTilesX(),
					quadCLT_main.tp.getTilesY());
			//			  }
			PoleProcessor pp = poleProcessor_persistent;
			selection =  quadCLT_main.tp.rig_pre_poles_sel.clone();

			int num_bugs_corrected= pp.zero_tiles_check(
					"pre-poles",
					quadCLT_main.tp.rig_pre_poles_ds, // double [][] ds,
					true); // boolean debug)
			if (num_bugs_corrected > 0) {
				System.out.println("Corrected "+num_bugs_corrected+" bugs, where disparity is <= 0 with non-zero strength before poles");
			}


			double [][] poles_ds = pp.processPoles(
					quadCLT_main,  // QuadCLT            quadCLT_main,  // tiles should be set
					quadCLT_aux,   //QuadCLT            quadCLT_aux,
					quadCLT_main.tp.rig_pre_poles_ds, //  double [][]                                    src_ds, // source disparity, strength pair
					selection, // boolean []                                     selection, // source tile selection, will be modified
					clt_parameters, // EyesisCorrectionParameters.CLTParameters       clt_parameters,
					threadsMax,  // final int                                      threadsMax,  // maximal number of threads to launch
					updateStatus, // final boolean                                  updateStatus,

					(debugLevel > -2) ? clt_parameters.poles.poles_debug_level:debugLevel);//
			//					 (debugLevel > -4) ? clt_parameters.poles.poles_debug_level:debugLevel);//
			//					 debugLevel); // final int                                      globalDebugLevel)

			num_bugs_corrected= pp.zero_tiles_check(
					"post-poles",
					poles_ds, // double [][] ds,
					true); // boolean debug)
			if (num_bugs_corrected > 0) {
				System.out.println("Corrected "+num_bugs_corrected+" bugs, where disparity is <= 0 with non-zero strength after poles");
			}

			quadCLT_main.tp.rig_post_poles_ds =        poles_ds; // Rig disparity and strength before processing poles
			quadCLT_main.tp.rig_post_poles_sel =       selection; // Rig tile selection before processing poles
			if (debugLevel > -4) {
				System.out.println("Saved quadCLT_main.tp.rig_post_poles_ds and quadCLT_main.tp.rig_post_poles_sel");
			}

			f_rig_disparity_strength = poles_ds;
		} else {
			if (debugLevel > -4) {
				System.out.println("Skipped processing poles");
			}
		}

		// Add areas from the main camera that correspond do occlusions of the aux (and so missing rig data). No real analysis of the occlusions here,
		// just filling gaps where rig data is 0.0 strengtth

		if (clt_parameters.rig.oc_fill_aux_occl) {
			if (debugLevel > -4) {
				System.out.println("Trying to fill rig gaps with main camera data (assuming aux cam occlusions)");
			}
		//scan_last
			double [][] main_last_scan = quadCLT_main.tp.getShowDS(scan_last,
					false); // boolean force_final);
			double disp_scale_inter = quadCLT_main.geometryCorrection.getDisparityRadius()/quadCLT_aux.geometryCorrection.getBaseline(); // 4.8
			// TODO: Move to parameters
			double      disp_atol = 0.2;
			double      disp_rtol = 0.1;
			double      strength_boost = 1.5;
//			boolean

			boolean [] sel_occl =  selectAuxOcclusions(
					main_last_scan,                // double [][] ds_main,
					f_rig_disparity_strength,      // double [][] ds_rig,
					disp_scale_inter,              // double      disp_rig_scale,
					clt_parameters.rig.oc_min_strength, // double      min_strength,
					clt_parameters.grow_disp_max,                 // double      max_disparity,
					disp_atol,                     // double      disp_atol,
					disp_rtol,                     // double      disp_rtol,
					clt_parameters.transform_size, // int         tileSize,     // 8
					1);                            // int         debugLevel);
			if (sel_occl != null) {
				if (debugLevel > -4) {
					double [][] dbg_img = new double [5][];
					String [] titles = {"disp_main","disp_rig","disp_occl","str_main","str_rig"};
					dbg_img[0] = main_last_scan[0];
					dbg_img[1] = f_rig_disparity_strength[0];
					dbg_img[2] = main_last_scan[0].clone();
					dbg_img[3] = main_last_scan[1];
					dbg_img[4] = f_rig_disparity_strength[1];
					for (int nTile = 0; nTile < sel_occl.length; nTile++) {
						if (!sel_occl[nTile]) {
							dbg_img[2][nTile] = Double.NaN;
						}
					}
					(new ShowDoubleFloatArrays()).showArrays(
							dbg_img,
							tilesX,
							dbg_img[0].length/tilesX,
							true,
							quadCLT_main.image_name+"Aux_occlusions",
							titles);
				}
				int num_replaced_occl = 0;
				for (int nTile = 0; nTile <  f_rig_disparity_strength[1].length; nTile++) { // only missing tiles
//					if (nTile== 50582) {
//						System.out.println("Replacing occlusions, nTile="+nTile);
//					}
					if (sel_occl[nTile]) {
						f_rig_disparity_strength[0][nTile] = main_last_scan[0][nTile];
						f_rig_disparity_strength[1][nTile] = strength_boost*  main_last_scan[1][nTile];
						selection[nTile] = true;
						num_replaced_occl ++;
					}
				}
				if (debugLevel > -4) {
					System.out.println("Replaced "+num_replaced_occl+" empty tiles in rig DSI with main camera data");
				}
			} else {
				if (debugLevel > -4) {
					System.out.println("No occluded tiles replaced");
				}

			}

/*
			int num_replaced_occl = 0;
			for (int nTile = 0; nTile <  f_rig_disparity_strength[1].length; nTile++) { // only missing tiles
				if (nTile== 50582) {
					System.out.println("Replacing occlusions, nTile="+nTile);
				}
				if (f_rig_disparity_strength[1][nTile] <= 0.0) { // only missing tiles
					if (    (main_last_scan[0][nTile] >= clt_parameters.rig.oc_min_disparity) && // only near objects
							(main_last_scan[1][nTile] >= clt_parameters.rig.oc_min_strength)) { // only strong tiles
						f_rig_disparity_strength[0][nTile] = main_last_scan[0][nTile];
						f_rig_disparity_strength[1][nTile] = main_last_scan[1][nTile];
						selection[nTile] = true;
						num_replaced_occl ++;
					}
				}
			}
*/
		}

		CLTPass3d rig_scan = quadCLT_main.tp.compositeScan(
				f_rig_disparity_strength[0],  // final double []             disparity,
				f_rig_disparity_strength[1],  // final double []             strength,
				selection,                  // final boolean []            selected,
				debugLevel);                // final int                   debugLevel)
		rig_scan.texture_tiles = scan_last.texture_tiles;

		// scan_last

		quadCLT_main.tp.clt_3d_passes.add(rig_scan);
		quadCLT_main.tp.saveCLTPasses(true);       // rig pass
	}

	public boolean [] selectAuxOcclusions(
			double [][] ds_main,
			double [][] ds_rig,
			double      disp_rig_scale,
			double      min_strength,
			double      max_disparity,
			double      disp_atol,
			double      disp_rtol,
			int         tileSize,     // 8
			int         debugLevel)
	{
		TileNeibs         tnImage =  biCamDSI_persistent.tnImage;
		int tilesX=  tnImage.getSizeX();
		int tilesY=  tnImage.getSizeY();
		int numTiles = tilesY*tilesX;
		boolean [] occl_sel = new boolean [numTiles];
		int num_replaced_occl = 0;
		for (int nTile = 0; nTile < numTiles; nTile++) { // only missing tiles
			if (nTile== 47120) {
				System.out.println("Replacing occlusions, nTile="+nTile);
			}
			if ((ds_rig[1][nTile] <= 0.0) && (ds_main[1][nTile] >= min_strength)) { // only missing tiles
				// See if aux camera may be potentially occluded for this tile
				int tileX = nTile % tilesX;
				int tileY = nTile / tilesX;
				double disp = ds_main[0][nTile];
				for (int tx = tileX+1; tx < tilesX; tx++) {
					disp += tileSize*disp_rig_scale; // disparity of the next tile right (here we assume aux is to the right of the main!)
					if (disp > max_disparity) {
						break;
					}
					int nTile1 = tnImage.getIndex(tx, tileY);
					if (nTile1>= 0) {
						double mn = (ds_rig[1][nTile1] > 0.0) ? ds_rig[0][nTile1]: Double.NaN;
						double mx = mn;
						for (int dir = 0; dir < 8; dir++) {
							int nTile2 = tnImage.getNeibIndex(nTile1, dir);
							if (nTile2 >= 0) if ((ds_rig[1][nTile2] > 0.0) && !Double.isNaN(ds_rig[0][nTile2])){
								if (!(ds_rig[0][nTile2] > mn)) mn = ds_rig[0][nTile2]; // mn==NaN - assign
								if (!(ds_rig[0][nTile2] < mx)) mx = ds_rig[0][nTile2]; // mn==NaN - assign
							}
						}
						if (!Double.isNaN(mn)) {
							double disp_avg = 0.5*(mn+mx);
							double disp_tol = disp_atol + disp_avg*disp_rtol;
							if (mx < (disp_avg + disp_tol)) mx =  disp_avg + disp_tol;
							if (mn < (disp_avg - disp_tol)) mn =  disp_avg - disp_tol;
							if ((disp >= mn) && (disp <=mx)) { // that may be occlusion
								occl_sel[nTile] = true;
								num_replaced_occl ++;
								if (debugLevel > 0) {
									System.out.println("Tile ("+tileX+", "+tileY+") may have an occlusion for aux camera by ("+
											(nTile1 % tilesX)+", "+(nTile1 / tilesX)+")");
								}
								break; // for (int nTile = 0; nTile < numTiles; nTile++)
							}
						}
					}
				}
			}
		}
		if (debugLevel > 0) {
			System.out.println("Found "+num_replaced_occl+" potential aux camera occluded tiles");
		}
		if (num_replaced_occl == 0) {
			occl_sel = null; // not caller to bother
		}
		return occl_sel;
	}

	public double [][][][][] getRigTextures( // never used
			boolean                                  need_master,
			boolean                                  need_aux,
			double []                                disparity, // non-nan - measure
			QuadCLT                                  quadCLT_main,  // tiles should be set
			QuadCLT                                  quadCLT_aux,
			CLTParameters clt_parameters,
			final int                                threadsMax,  // maximal number of threads to launch
			final boolean                            updateStatus,
			final int                                debugLevel) // throws Exception
	{
		final int tile_op_all = clt_parameters.tile_task_op; //FIXME Use some constant?
		final int tilesX = quadCLT_main.tp.getTilesX();
		final int tilesY = quadCLT_main.tp.getTilesY();
		final int [][]     tile_op = new int [tilesY][tilesX];
		final double [][]  disparity_array = new double [tilesY][tilesX];
		boolean [] need_textures = {need_master, need_aux};
		final double [][][][][] texture_tiles = new double [need_textures.length][][][][];
		for (int ncam = 0; ncam < need_textures.length; ncam++) if (need_textures[ncam]) {
			texture_tiles[ncam] = new double [tilesY][tilesX][][];
			for (int tileY = 0; tileY < tilesY; tileY++){
				for (int tileX = 0; tileX < tilesX; tileX++){
					texture_tiles[ncam][tileY][tileX] = null;
					int nTile = tileY*tilesX + tileX;
					if (!Double.isNaN(disparity[nTile])) {
						tile_op[tileY][tileX] = tile_op_all;
						disparity_array[tileY][tileX] = disparity[nTile];
					}
				}
			}
		}
		ImageDtt image_dtt = new ImageDtt(
				clt_parameters.transform_size,
				quadCLT_main.isMonochrome(),
				quadCLT_main.isLwir(),
				clt_parameters.getScaleStrength(false));
		image_dtt.clt_bi_quad (
				clt_parameters,                       // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
				clt_parameters.getFatZero(image_dtt.isMonochrome()),              // final double              fatzero,         // May use correlation fat zero from 2 different parameters - fat_zero and rig.ml_fatzero
				false,                                //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
				0,                                    // final int  // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
				// does not matter here
				true,                                 // boolean                                  no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
				tile_op,                              // final int [][]            tile_op_main,    // [tilesY][tilesX] - what to do - 0 - nothing for this tile
				disparity_array,                      // final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
				quadCLT_main.image_data,              // final double [][][]       image_data_main, // first index - number of image in a quad
				quadCLT_aux.image_data,               // final double [][][]       image_data_aux,  // first index - number of image in a quad
				quadCLT_main.saturation_imp,          // final boolean [][]        saturation_main, // (near) saturated pixels or null
				quadCLT_aux.saturation_imp,           // final boolean [][]        saturation_aux,  // (near) saturated pixels or null
				// correlation results - combo will be for the correation between two quad cameras
				null,                                 // final double [][][][]     clt_corr_combo,  // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				null, // disparity_bimap,             // final double [][]    disparity_bimap, // [23][tilesY][tilesX]
				null, // ml_data,                     // 	final double [][]         ml_data,         // data for ML - 10 layers - 4 center areas (3x3, 5x5,..) per camera-per direction, 1 - composite, and 1 with just 1 data (target disparity)
				texture_tiles[0],                     // final double [][][][]     texture_tiles_main, // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
				texture_tiles[1],                     // final double [][][][]     texture_tiles_aux,  // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
				quadCLT_main.tp.getTilesX()*image_dtt.transform_size, // final int                 width,

				quadCLT_main.getGeometryCorrection(), // final GeometryCorrection  geometryCorrection_main,
				quadCLT_aux.getGeometryCorrection(),  // final GeometryCorrection  geometryCorrection_aux,
				quadCLT_main.getCLTKernels(),         // final double [][][][][][] clt_kernels_main, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				quadCLT_aux.getCLTKernels(),          // final double [][][][][][] clt_kernels_aux,  // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				clt_parameters.corr_magic_scale,      // final double              corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
				false, // true,                                 // 	final boolean             keep_clt_data,
//				woi_tops,                             // final int [][]            woi_tops,
				null,                                 // final double [][][]       ers_delay,        // if not null - fill with tile center acquisition delay
				threadsMax,                           // final int                 threadsMax,  // maximal number of threads to launch
				debugLevel-2);                        // final int                 globalDebugLevel);

		return texture_tiles;

	}



	public void copyJP4src(
			QuadCLT                                  quadCLT_main,  // tiles should be set
			QuadCLT                                  quadCLT_aux,
			CLTParameters clt_parameters,
			final int                                debugLevel) // throws Exception
	{
		//		  quadCLT_main.writeKml(debugLevel);
		//		  quadCLT_main.writeRatingFile(debugLevel);


		String [] sourceFiles_main=quadCLT_main.correctionsParameters.getSourcePaths();
		//



		//		  String [] sourceFiles_aux=quadCLT_main.correctionsParameters.getSourcePaths();
		String set_name = quadCLT_main.image_name;
		if (set_name == null ) {
			QuadCLT.SetChannels [] set_channels=quadCLT_main.setChannels(debugLevel);
			set_name = set_channels[0].set_name;
		}

		QuadCLT.SetChannels [] set_channels_main = quadCLT_main.setChannels(set_name,debugLevel); // only for specified image timestamp
		QuadCLT.SetChannels [] set_channels_aux =  quadCLT_aux.setChannels(set_name,debugLevel);

		ArrayList<String> path_list = new ArrayList<String>();
		for (int i = 0; i < set_channels_main.length; i++) {
			for (int fn:set_channels_main[i].file_number) {
				path_list.add(sourceFiles_main[fn]);
			}
		}
		for (int i = 0; i < set_channels_aux.length; i++) {
			for (int fn:set_channels_aux[i].file_number) {
				path_list.add(sourceFiles_main[fn]);
			}
		}

		quadCLT_main.writeKml(debugLevel ); // also generated with x3d model

		String jp4_copy_path= quadCLT_main.correctionsParameters.selectX3dDirectory(
				set_name, // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
				quadCLT_main.correctionsParameters.jp4SubDir,
				true,  // smart,
				true);  //newAllowed, // save
		File dir = (new File(jp4_copy_path)); // .getParentFile();
		if (!dir.exists()){
			dir.mkdirs();
			System.out.println("Created "+dir);
		}
		for (String fname:path_list) {
			if (fname.contains(quadCLT_main.image_name)) { // only files containing set name // TODO:improve
				File file = new File(fname);
				try {
					Files.copy(
							(file).toPath(),
							(new File(jp4_copy_path + Prefs.getFileSeparator()+file.getName())).toPath(),
							StandardCopyOption.REPLACE_EXISTING);
					System.out.println("Copied "+fname+" -> "+dir);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					System.out.println("Failed to copy "+fname+" -> "+dir);
				}
			}
		}
		System.out.println("jp4_copy_path = "+jp4_copy_path);
		//		  System.out.println("Do something useful here");
	}



	public void outputMLData(
			QuadCLT                                  quadCLT_main,  // tiles should be set
			QuadCLT                                  quadCLT_aux,
			CLTParameters clt_parameters,
			String                                   ml_directory,       // full path or null (will use config one)
			final int                                threadsMax,  // maximal number of threads to launch
			final boolean                            updateStatus,
			final int                                debugLevel) // throws Exception
	{
		if ((quadCLT_main.tp == null) || (quadCLT_main.tp.rig_pre_poles_ds == null)) {
			String msg = "DSI data not available. Please run \"Ground truth\" first";
			IJ.showMessage("Error",msg);
			System.out.println(msg);
			return;
		}
		boolean post_poles =  clt_parameters.rig.ml_poles;
		double [][] rig_disparity_strength = clt_parameters.rig.ml_poles?quadCLT_main.tp.rig_post_poles_ds : quadCLT_main.tp.rig_pre_poles_ds;
		if ((quadCLT_main.tp.main_ds_ml == null) && (this.dsi != null) &&  (this.dsi[DSI_DISPARITY_MAIN] != null) &&  (this.dsi[DSI_STRENGTH_MAIN] != null)){
			quadCLT_main.tp.main_ds_ml = new double[2][];
			quadCLT_main.tp.main_ds_ml[0] = this.dsi[DSI_DISPARITY_MAIN];
			quadCLT_main.tp.main_ds_ml[1] = this.dsi[DSI_STRENGTH_MAIN];
		}
		double [][] main_disparity_strength = quadCLT_main.tp.main_ds_ml;
		if (rig_disparity_strength == null) {
			System.out.println("DSI data for the scene after poles extraction is not available. You may enable it and re-run \"Ground truth\" command or run \"Poles GT\"");
			rig_disparity_strength = quadCLT_main.tp.rig_pre_poles_ds;
			System.out.println("Using pre-poles data for ML output");
		}
		if (debugLevel > -6) {
			if (post_poles) {
				System.out.println("==== Generating ML data for the DSI that includes extracted vertical poles ====");
			} else {
				System.out.println("==== Generating ML data for the DSI that DOES NOT include extracted vertical poles ====");
			}
		}
		// Create filtered/expanded min and rig DSI
		// FIXME: Make quadCLT_main.tp.main_ds_ml dusing COMBO_DSI write (too late) or just write earlier?
		double [][] main_dsi = null;
		if (main_disparity_strength != null) {
			main_dsi = enhanceMainDSI(
				main_disparity_strength,              // double [][] main_dsi,
				rig_disparity_strength,               // double [][] rig_dsi,
				clt_parameters.rig.ml_rig_tolerance,  // double      rig_tolerance,
				clt_parameters.rig.ml_rnd_offset,     // double      rnd_offset,
				clt_parameters.rig.ml_main_tolerance, // double      main_tolerance,
				clt_parameters.rig.ml_grow_steps,     // int         grow_steps,
				clt_parameters.rig.ml_grow_mode,      // int         grow_mode,
				clt_parameters.rig.ml_new_strength,   // double      new_strength)
//				debugLevel + 3);                        // int         debugLevel);
				debugLevel + 0);                        // int         debugLevel);
		}
		double [][] rig_dsi = enhanceMainDSI(
				rig_disparity_strength,              // double [][] main_dsi,
				rig_disparity_strength,               // double [][] rig_dsi,
				clt_parameters.rig.ml_rig_tolerance,  // double      rig_tolerance,
				clt_parameters.rig.ml_rnd_offset,     // double      rnd_offset,
				clt_parameters.rig.ml_main_tolerance, // double      main_tolerance,
				clt_parameters.rig.ml_grow_steps,     // int         grow_steps,
				clt_parameters.rig.ml_grow_mode,      // int         grow_mode,
				clt_parameters.rig.ml_new_strength,   // double      new_strength)
//				debugLevel + 3);                        // int         debugLevel);
				debugLevel + 0);                        // int         debugLevel);

		if (ml_directory == null) ml_directory= quadCLT_main.correctionsParameters.selectMlDirectory(
				quadCLT_main.image_name,
				true,  // smart,
				true);  //newAllowed, // save
		Correlation2d corr2d = new Correlation2d(
				clt_parameters.img_dtt,              // ImageDttParameters  imgdtt_params,
				clt_parameters.transform_size,             // int transform_size,
				2.0,                        //  double wndx_scale, // (wndy scale is always 1.0)
				quadCLT_main.isMonochrome(),
				(debugLevel > -1));   //   boolean debug)

// Create test data that does not rely on the rig measurements
		String img_name_main =quadCLT_main.image_name+"-ML_DATA-";
		// zero random offset:
		if ( clt_parameters.rig.ml_main && (main_dsi != null)) {
			double [][] ml_data_main = remeasureRigML(
					0.0,                     // double               disparity_offset_low,
					Double.NaN, // 0.0,                    // double               disparity_offset_high,
					quadCLT_main,                             // QuadCLT              quadCLT_main,    // tiles should be set
					quadCLT_aux,                              // QuadCLT              quadCLT_aux,
					//	    			  disparity_bimap,                          // double [][]          src_bimap,
					main_dsi[0],               // double []            disparity_main, // main camera disparity to use - if null, calculate from the rig one
					rig_disparity_strength[0],                // double []            disparity,
					rig_disparity_strength[1],                // double []            strength,

					clt_parameters,                           // EyesisCorrectionParameters.CLTParameters clt_parameters,
					clt_parameters.rig.ml_hwidth,             // int ml_hwidth
					clt_parameters.rig.ml_fatzero,            // double               fatzero,
					//change if needed?
					0, //  int              lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using (2*notch_mode+1)^2 square
					threadsMax,                               // final int            threadsMax,  // maximal number of threads to launch
					updateStatus,                             // final boolean        updateStatus,
					debugLevel);                              // final int            debugLevel);
			saveMlFile(
					img_name_main,                                 // String               ml_title,
					ml_directory,                             // String               ml_directory,
					Double.NaN,                               // double               disp_offset_low,
					Double.NaN,                               // double               disp_offset_high,
					quadCLT_main,                             // QuadCLT              quadCLT_main,
					quadCLT_aux,                              // QuadCLT              quadCLT_aux,
					corr2d,                                   //Correlation2d        corr2d, // to access "other" layer
					clt_parameters.rig.ml_8bit,               // boolean              use8bpp,
					clt_parameters.rig.ml_limit_extrim,       // double               limit_extrim,
					clt_parameters.rig.ml_keep_aux,           // boolean              keep_aux,
					clt_parameters.rig.ml_keep_inter,         // boolean              keep_inter,
					clt_parameters.rig.ml_keep_hor_vert,      // boolean              keep_hor_vert,
					clt_parameters.rig.ml_keep_tbrl,          // boolean              ml_keep_tbrl,
					clt_parameters.rig.ml_keep_debug,         // boolean              keep_debug,
					clt_parameters.rig.ml_fatzero,            // double               ml_fatzero,
					clt_parameters.rig.ml_hwidth,             // int                  ml_hwidth,
					ml_data_main,                             // double [][]          ml_data,
					clt_parameters.rig.ml_show_ml,            // boolean              show,
					debugLevel);                              // int                  debugLevel

		}
		// Create images from main cameras, but adding random disparity offset
		if ( clt_parameters.rig.ml_main_rnd && (main_dsi != null)) {
			double [][] ml_data_main = remeasureRigML(
					0.0,                     // double               disparity_offset_low,
					clt_parameters.rig.ml_disparity_sweep, // 0.0,                    // double               disparity_offset_high,
					quadCLT_main,                             // QuadCLT              quadCLT_main,    // tiles should be set
					quadCLT_aux,                              // QuadCLT              quadCLT_aux,
					//	    			  disparity_bimap,                          // double [][]          src_bimap,
					main_dsi[0],               // double []            disparity_main, // main camera disparity to use - if null, calculate from the rig one
					rig_disparity_strength[0],                // double []            disparity,
					rig_disparity_strength[1],                // double []            strength,

					clt_parameters,                           // EyesisCorrectionParameters.CLTParameters clt_parameters,
					clt_parameters.rig.ml_hwidth,             // int ml_hwidth
					clt_parameters.rig.ml_fatzero,            // double               fatzero,
					//change if needed?
					0, //  int              lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using (2*notch_mode+1)^2 square
					threadsMax,                               // final int            threadsMax,  // maximal number of threads to launch
					updateStatus,                             // final boolean        updateStatus,
					debugLevel);                              // final int            debugLevel);
			saveMlFile(
					img_name_main,                            // String               ml_title,
					ml_directory,                             // String               ml_directory,
					Double.NaN,                               // double               disp_offset_low,
					clt_parameters.rig.ml_disparity_sweep,    // double               disp_offset_high,
					quadCLT_main,                             // QuadCLT              quadCLT_main,
					quadCLT_aux,                              // QuadCLT              quadCLT_aux,
					corr2d,                                   //Correlation2d        corr2d, // to access "other" layer
					clt_parameters.rig.ml_8bit,               // boolean              use8bpp,
					clt_parameters.rig.ml_limit_extrim,       // double               limit_extrim,
					clt_parameters.rig.ml_keep_aux,           // boolean              keep_aux,
					clt_parameters.rig.ml_keep_inter,         // boolean              keep_inter,
					clt_parameters.rig.ml_keep_hor_vert,      // boolean              keep_hor_vert,
					clt_parameters.rig.ml_keep_tbrl,          // boolean              ml_keep_tbrl,
					clt_parameters.rig.ml_keep_debug,         // boolean              keep_debug,
					clt_parameters.rig.ml_fatzero,            // double               ml_fatzero,
					clt_parameters.rig.ml_hwidth,             // int                  ml_hwidth,
					ml_data_main,                             // double [][]          ml_data,
					clt_parameters.rig.ml_show_ml,            // boolean              show,
					debugLevel);                              // int                  debugLevel

		}

		// Create images from main cameras, but adding random disparity offset
		if ( clt_parameters.rig.ml_rig_rnd && (rig_dsi != null)) {
			double [][] ml_data_main = remeasureRigML(
					0.0,                     // double               disparity_offset_low,
					clt_parameters.rig.ml_disparity_sweep, // 0.0,                    // double               disparity_offset_high,
					quadCLT_main,                             // QuadCLT              quadCLT_main,    // tiles should be set
					quadCLT_aux,                              // QuadCLT              quadCLT_aux,
					rig_dsi[0],                               // double []            disparity_main, // main camera disparity to use - if null, calculate from the rig one
					rig_disparity_strength[0],                // double []            disparity,
					rig_disparity_strength[1],                // double []            strength,

					clt_parameters,                           // EyesisCorrectionParameters.CLTParameters clt_parameters,
					clt_parameters.rig.ml_hwidth,             // int ml_hwidth
					clt_parameters.rig.ml_fatzero,            // double               fatzero,
					//change if needed?
					0, //  int              lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using (2*notch_mode+1)^2 square
					threadsMax,                               // final int            threadsMax,  // maximal number of threads to launch
					updateStatus,                             // final boolean        updateStatus,
					debugLevel);                              // final int            debugLevel);
			saveMlFile(
					img_name_main,                            // String               ml_title,
					ml_directory,                             // String               ml_directory,
					Double.NaN,                               // double               disp_offset_low,
					-clt_parameters.rig.ml_disparity_sweep,   // double               disp_offset_high,
					quadCLT_main,                             // QuadCLT              quadCLT_main,
					quadCLT_aux,                              // QuadCLT              quadCLT_aux,
					corr2d,                                   //Correlation2d        corr2d, // to access "other" layer
					clt_parameters.rig.ml_8bit,               // boolean              use8bpp,
					clt_parameters.rig.ml_limit_extrim,       // double               limit_extrim,
					clt_parameters.rig.ml_keep_aux,           // boolean              keep_aux,
					clt_parameters.rig.ml_keep_inter,         // boolean              keep_inter,
					clt_parameters.rig.ml_keep_hor_vert,      // boolean              keep_hor_vert,
					clt_parameters.rig.ml_keep_tbrl,          // boolean              ml_keep_tbrl,
					clt_parameters.rig.ml_keep_debug,         // boolean              keep_debug,
					clt_parameters.rig.ml_fatzero,            // double               ml_fatzero,
					clt_parameters.rig.ml_hwidth,             // int                  ml_hwidth,
					ml_data_main,                             // double [][]          ml_data,
					clt_parameters.rig.ml_show_ml,            // boolean              show,
					debugLevel);                              // int                  debugLevel

		}
		// Old code:
		// Create images from tig with random offset
		// TODO: add extrapolation, similar to main camera?
		for (int sweep_step = 0; sweep_step < clt_parameters.rig.ml_sweep_steps; sweep_step++){
			double disparity_offset_low = 0;
			double disparity_offset_high = Double.NaN;
			if (clt_parameters.rig.ml_sweep_steps > 1) {
				disparity_offset_low = clt_parameters.rig.ml_disparity_sweep * (2.0 * sweep_step/(clt_parameters.rig.ml_sweep_steps - 1.0) -1.0);
				if (clt_parameters.rig.ml_randomize) {
					disparity_offset_high = clt_parameters.rig.ml_disparity_sweep * (2.0 * (sweep_step + 1)/(clt_parameters.rig.ml_sweep_steps - 1.0) -1.0);
				}
			}

			String img_name = clt_parameters.rig.ml_randomize ? quadCLT_main.image_name+"-ML_DATARND-" : quadCLT_main.image_name+"-ML_DATA-";
			double [][] ml_data = remeasureRigML(
					disparity_offset_low,                     // double               disparity_offset_low,
					disparity_offset_high,                    // double               disparity_offset_high,
					quadCLT_main,                             // QuadCLT              quadCLT_main,    // tiles should be set
					quadCLT_aux,                              // QuadCLT              quadCLT_aux,
					//	    			  disparity_bimap,                          // double [][]          src_bimap,
					null, // double []                                      disparity_main, // main camera disparity to use - if null, calculate from the rig one
					rig_disparity_strength[0],                // double []            disparity,
					rig_disparity_strength[1],                // double []            strength,

					clt_parameters,                           // EyesisCorrectionParameters.CLTParameters clt_parameters,
					clt_parameters.rig.ml_hwidth,             // int ml_hwidth
					clt_parameters.rig.ml_fatzero,            // double               fatzero,
					//change if needed?
					0, //  int              lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using (2*notch_mode+1)^2 square
					threadsMax,                               // final int            threadsMax,  // maximal number of threads to launch
					updateStatus,                             // final boolean        updateStatus,
					debugLevel);                              // final int            debugLevel);
			saveMlFile(
					img_name,                                 // String               ml_title,
					ml_directory,                             // String               ml_directory,
					disparity_offset_low,                     // double               disp_offset_low,
					disparity_offset_high,                    // double               disp_offset_high,
					quadCLT_main,                             // QuadCLT              quadCLT_main,
					quadCLT_aux,                              // QuadCLT              quadCLT_aux,
					corr2d,                                   //Correlation2d        corr2d, // to access "other" layer
					clt_parameters.rig.ml_8bit,               // boolean              use8bpp,
					clt_parameters.rig.ml_limit_extrim,       // double               limit_extrim,
					clt_parameters.rig.ml_keep_aux,           // boolean              keep_aux,
					clt_parameters.rig.ml_keep_inter,         // boolean              keep_inter,
					clt_parameters.rig.ml_keep_hor_vert,      // boolean              keep_hor_vert,
					clt_parameters.rig.ml_keep_tbrl,          // boolean              ml_keep_tbrl,
					clt_parameters.rig.ml_keep_debug,         // boolean              keep_debug,
					clt_parameters.rig.ml_fatzero,            // double               ml_fatzero,
					clt_parameters.rig.ml_hwidth,             // int                  ml_hwidth,
					ml_data,                                  // double [][]          ml_data,
					clt_parameters.rig.ml_show_ml,            // boolean              show,
					debugLevel);                              // int                  debugLevel
			Runtime.getRuntime().gc();
			if (clt_parameters.rig.ml_randomize) {
				System.out.println("Generated ML data, randomized offset = "+String.format("%8.5f...%8.5f",disparity_offset_low,disparity_offset_high)+", --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
			} else {
				System.out.println("Generated ML data, offset = "+String.format("%8.5f",disparity_offset_low)+", --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
			}
			if (clt_parameters.rig.ml_randomize && (sweep_step == (clt_parameters.rig.ml_sweep_steps-2))) { // reduce by 1 for random
				break;
			}
		}
	}



	public void outputMLDataLwir(
			QuadCLT                                  quadCLT_aux,
			CLTParameters clt_parameters,
			String                                   ml_directory,       // full path or null (will use config one)
			final int                                threadsMax,  // maximal number of threads to launch
			final boolean                            updateStatus,
			final int                                debugLevel) // throws Exception
	{

		if ((quadCLT_aux == null) || (quadCLT_aux.tp == null) || (this.dsi_aux_from_main  == null)) {
			String msg = "DSI data for AUX camera is not available. Please generateload it first";
			IJ.showMessage("Error",msg);
			System.out.println(msg);
			return;
		}
		String img_name =quadCLT_main.image_name+"-ML_DATA-";
		if (ml_directory == null) ml_directory= quadCLT_main.correctionsParameters.selectMlDirectory(
				quadCLT_main.image_name,
				true,  // smart,
				true);  //newAllowed, // save
		Correlation2d corr2d = new Correlation2d(
				clt_parameters.img_dtt,         // ImageDttParameters  imgdtt_params,
				clt_parameters.transform_size,  // int transform_size,
				2.0,                            // double wndx_scale, // (wndy scale is always 1.0)
				quadCLT_aux.isMonochrome(),     // boolean monochrome,
				(debugLevel > -1));             // boolean debug)

		if (clt_parameters.rig.ml_aux_ag) {
			double [][] ml_data =  remeasureAuxML(
					quadCLT_aux,                                    // QuadCLT                                        quadCLT_aux,
					this.dsi_aux_from_main[QuadCLT.FGBG_DISPARITY], // double []                                      disparity,
					this.dsi_aux_from_main,                         // double [][]                                    dsi_aux_from_main, // Main camera DSI converted into the coordinates of the AUX one (see QuadCLT.FGBG_TITLES)-
					clt_parameters,                                 // CLTParameters       clt_parameters,
					clt_parameters.rig.ml_hwidth,                   // int                                            ml_hwidth,
					clt_parameters.rig.ml_fatzero,                  // double                                         fatzero,
					threadsMax,                                     // final int            threadsMax,  // maximal number of threads to launch
					updateStatus,                                   // final boolean        updateStatus,
					debugLevel);                                    // final int            debugLevel);

			saveMlFile(
					img_name,                                 // String               ml_title,
					ml_directory,                             // String               ml_directory,
					Double.NaN,                               // double               disp_offset_low,
					0.0,   // AG                              // double               disp_offset_high,
					null, // quadCLT_main,                    // QuadCLT              quadCLT_main,
					quadCLT_aux,                              // QuadCLT              quadCLT_aux,
					corr2d,                                   // Correlation2d        corr2d, // to access "other" layer
					clt_parameters.rig.ml_8bit,               // boolean              use8bpp, // true
					clt_parameters.rig.ml_limit_extrim,       // double               limit_extrim, 1E-5
					clt_parameters.rig.ml_keep_aux,           // boolean              keep_aux, false
					clt_parameters.rig.ml_keep_inter,         // boolean              keep_inter, // false
					clt_parameters.rig.ml_keep_hor_vert,      // boolean              keep_hor_vert, // true
					clt_parameters.rig.ml_keep_tbrl,          // boolean              ml_keep_tbrl,  // false
					clt_parameters.rig.ml_keep_debug,         // boolean              keep_debug,    // false
					clt_parameters.rig.ml_fatzero,            // double               ml_fatzero,    // 0.05
					clt_parameters.rig.ml_hwidth,             // int                  ml_hwidth,     // 4
					ml_data,                                  // double [][]          ml_data,       //false
					clt_parameters.rig.ml_show_ml,            // boolean              show,
					debugLevel);                              // int                  debugLevel
		}
		if (clt_parameters.rig.ml_aux_fg) {
			double [][] ml_data =  remeasureAuxML(
					quadCLT_aux,                                  // QuadCLT                                        quadCLT_aux,
					this.dsi_aux_from_main[QuadCLT.FGBG_FG_DISP], // double []                                      disparity,
					this.dsi_aux_from_main,                       // double [][]                                    dsi_aux_from_main, // Main camera DSI converted into the coordinates of the AUX one (see QuadCLT.FGBG_TITLES)-
					clt_parameters,                               // CLTParameters       clt_parameters,
					clt_parameters.rig.ml_hwidth,                 // int                                            ml_hwidth,
					clt_parameters.rig.ml_fatzero,                // double                                         fatzero,
					threadsMax,                                   // final int            threadsMax,  // maximal number of threads to launch
					updateStatus,                                 // final boolean        updateStatus,
					debugLevel);                                  // final int            debugLevel);

			saveMlFile(
					img_name,                                 // String               ml_title,
					ml_directory,                             // String               ml_directory,
					Double.NaN,                               // double               disp_offset_low,
					1.0,  // FG                               // double               disp_offset_high,
					null, // quadCLT_main,                    // QuadCLT              quadCLT_main,
					quadCLT_aux,                              // QuadCLT              quadCLT_aux,
					corr2d,                                   // Correlation2d        corr2d, // to access "other" layer
					clt_parameters.rig.ml_8bit,               // boolean              use8bpp, // true
					clt_parameters.rig.ml_limit_extrim,       // double               limit_extrim, 1E-5
					clt_parameters.rig.ml_keep_aux,           // boolean              keep_aux, false
					clt_parameters.rig.ml_keep_inter,         // boolean              keep_inter, // false
					clt_parameters.rig.ml_keep_hor_vert,      // boolean              keep_hor_vert, // true
					clt_parameters.rig.ml_keep_tbrl,          // boolean              ml_keep_tbrl,  // false
					clt_parameters.rig.ml_keep_debug,         // boolean              keep_debug,    // false
					clt_parameters.rig.ml_fatzero,            // double               ml_fatzero,    // 0.05
					clt_parameters.rig.ml_hwidth,             // int                  ml_hwidth,     // 4
					ml_data,                                  // double [][]          ml_data,       //false
					clt_parameters.rig.ml_show_ml,            // boolean              show,
					debugLevel);                              // int                  debugLevel
		}

		if (clt_parameters.rig.ml_aux_bg) {
			double [][] ml_data =  remeasureAuxML(
					quadCLT_aux,                                  // QuadCLT                                        quadCLT_aux,
					this.dsi_aux_from_main[QuadCLT.FGBG_BG_DISP], // double []                                      disparity,
					this.dsi_aux_from_main,                       // double [][]                                    dsi_aux_from_main, // Main camera DSI converted into the coordinates of the AUX one (see QuadCLT.FGBG_TITLES)-
					clt_parameters,                               // CLTParameters       clt_parameters,
					clt_parameters.rig.ml_hwidth,                 // int                                            ml_hwidth,
					clt_parameters.rig.ml_fatzero,                // double                                         fatzero,
					threadsMax,                                   // final int            threadsMax,  // maximal number of threads to launch
					updateStatus,                                 // final boolean        updateStatus,
					debugLevel);                                  // final int            debugLevel);

			saveMlFile(
					img_name,                                 // String               ml_title,
					ml_directory,                             // String               ml_directory,
					Double.NaN,                               // double               disp_offset_low,
					-1.0,  // BG                              // double               disp_offset_high,
					null, // quadCLT_main,                    // QuadCLT              quadCLT_main,
					quadCLT_aux,                              // QuadCLT              quadCLT_aux,
					corr2d,                                   // Correlation2d        corr2d, // to access "other" layer
					clt_parameters.rig.ml_8bit,               // boolean              use8bpp, // true
					clt_parameters.rig.ml_limit_extrim,       // double               limit_extrim, 1E-5
					clt_parameters.rig.ml_keep_aux,           // boolean              keep_aux, false
					clt_parameters.rig.ml_keep_inter,         // boolean              keep_inter, // false
					clt_parameters.rig.ml_keep_hor_vert,      // boolean              keep_hor_vert, // true
					clt_parameters.rig.ml_keep_tbrl,          // boolean              ml_keep_tbrl,  // false
					clt_parameters.rig.ml_keep_debug,         // boolean              keep_debug,    // false
					clt_parameters.rig.ml_fatzero,            // double               ml_fatzero,    // 0.05
					clt_parameters.rig.ml_hwidth,             // int                  ml_hwidth,     // 4
					ml_data,                                  // double [][]          ml_data,       //false
					clt_parameters.rig.ml_show_ml,            // boolean              show,
					debugLevel);                              // int                  debugLevel
		}
		if (clt_parameters.rig.ml_aux_step > 0.0) {
			double [] target_disparity = new double [this.dsi_aux_from_main[QuadCLT.FGBG_DISPARITY].length];
			for (double disparity = clt_parameters.rig.ml_aux_low; disparity <= clt_parameters.rig.ml_aux_high; disparity += clt_parameters.rig.ml_aux_step) {
				for (int nt = 0; nt < target_disparity.length; nt++) {
					target_disparity[nt] = disparity;
				}
				double [][] ml_data =  remeasureAuxML(
						quadCLT_aux,                          // QuadCLT                                        quadCLT_aux,
						target_disparity,                     // double []                                      disparity,
						this.dsi_aux_from_main,               // double [][]                                    dsi_aux_from_main, // Main camera DSI converted into the coordinates of the AUX one (see QuadCLT.FGBG_TITLES)-
						clt_parameters,                       // CLTParameters       clt_parameters,
						clt_parameters.rig.ml_hwidth,         // int                                            ml_hwidth,
						clt_parameters.rig.ml_fatzero,        // double                                         fatzero,
						threadsMax,                           // final int            threadsMax,  // maximal number of threads to launch
						updateStatus,                         // final boolean        updateStatus,
						debugLevel);                          // final int            debugLevel);

				saveMlFile(
						img_name,                                 // String               ml_title,
						ml_directory,                             // String               ml_directory,
						disparity, // absolute disparity          // double               disp_offset_low,
						0.0,  // not used                         // double               disp_offset_high,
						null, // -> aux_mode (LWIR),              // QuadCLT              quadCLT_main,
						quadCLT_aux,                              // QuadCLT              quadCLT_aux,
						corr2d,                                   // Correlation2d        corr2d, // to access "other" layer
						clt_parameters.rig.ml_8bit,               // boolean              use8bpp, // true
						clt_parameters.rig.ml_limit_extrim,       // double               limit_extrim, 1E-5
						clt_parameters.rig.ml_keep_aux,           // boolean              keep_aux, false
						clt_parameters.rig.ml_keep_inter,         // boolean              keep_inter, // false
						clt_parameters.rig.ml_keep_hor_vert,      // boolean              keep_hor_vert, // true
						clt_parameters.rig.ml_keep_tbrl,          // boolean              ml_keep_tbrl,  // false
						clt_parameters.rig.ml_keep_debug,         // boolean              keep_debug,    // false
						clt_parameters.rig.ml_fatzero,            // double               ml_fatzero,    // 0.05
						clt_parameters.rig.ml_hwidth,             // int                  ml_hwidth,     // 4
						ml_data,                                  // double [][]          ml_data,       //false
						clt_parameters.rig.ml_show_ml,            // boolean              show,
						debugLevel);                              // int                  debugLevel
			}
		}
		Runtime.getRuntime().gc();
    	System.out.println("Generated ML data, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
	}


	public double [][] prepareRefineExistingDSI(
			QuadCLT            quadCLT_main,  // tiles should be set
			QuadCLT            quadCLT_aux,
			CLTParameters       clt_parameters,
			ColorProcParameters                            colorProcParameters, //
			ColorProcParameters                            colorProcParameters_aux, //
			final int                                      threadsMax,  // maximal number of threads to launch
			final boolean                                  updateStatus,
			final int                                      debugLevel) //  throws Exception
	{
		final int refine_inter = 2; // 3; // 3 - dx, 2 - disparity
		final int tilesX = quadCLT_main.tp.getTilesX();
		//		  final int tilesY = quadCLT_main.tp.getTilesY();
		if ((quadCLT_main == null) || (quadCLT_aux == null)) {
			System.out.println("QuadCLT instances are not initilaized");
			return null;
		}
		// verify main camera has measured data
		// Measure with target disparity == 0
		if ((quadCLT_main.tp == null) || (quadCLT_main.tp.clt_3d_passes == null) || (quadCLT_main.tp.clt_3d_passes.size() ==0)){
			System.out.println("No DSI data for the main camera is available. Please run \"CLT 3D\" command");
			return null;
		}
		// See if auxiliary camera has images configured, if not - do it now.
		if (quadCLT_aux.image_name ==null) {
			boolean aux_OK = quadCLT_aux.setupImageData(
					quadCLT_main.image_name,                             // String image_name,
					quadCLT_main.correctionsParameters.getSourcePaths(), // String [] sourceFiles,
					clt_parameters,                                      // EyesisCorrectionParameters.CLTParameters       clt_parameters,
					colorProcParameters,                                 //  ColorProcParameters                       colorProcParameters, //
					threadsMax,                                          // int threadsMax,
					debugLevel);                                         // int debugLevel);
			if (!aux_OK) {
				return null;
			}
		}
		// Re-measure background
		double [][] disparity_bimap_infinity = measureNewRigDisparity(
				quadCLT_main,      // QuadCLT             quadCLT_main,    // tiles should be set
				quadCLT_aux,       // QuadCLT             quadCLT_aux,
				null,              // double [][]         src_bimap,       // current state of measurements (or null for new measurement)
				0,                 // double              disparity,
				null,              // ArrayList<Integer>  tile_list,       // or null. If non-null - do not remeasure members of the list
				clt_parameters,    // EyesisCorrectionParameters.CLTParameters       clt_parameters,
				false,             //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
				0,                 // final int  // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
				// first measurement - use as set in parameters
				clt_parameters.rig.no_int_x0,           // boolean                                  no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
				threadsMax,        // final int           threadsMax,      // maximal number of threads to launch
				updateStatus,      // final boolean       updateStatus,
				debugLevel);       // final int           debugLevel);


		int [] num_new = new int[1];
		boolean [] trusted_infinity = 	  getTrustedDisparity(
				quadCLT_main,                            // QuadCLT            quadCLT_main,  // tiles should be set
				quadCLT_aux,                             // QuadCLT            quadCLT_aux,
				true,                                    // boolean            use_individual,
				clt_parameters.rig.min_trusted_strength, // double             min_combo_strength,    // check correlation strength combined for all 3 correlations
				clt_parameters.grow_disp_trust,          // double             max_trusted_disparity, // 4.0 -> change to rig_trust
				clt_parameters.rig.trusted_tolerance,    // double             trusted_tolerance,
				null,                                    // boolean []         was_trusted,
				disparity_bimap_infinity );              // double [][]        bimap // current state of measurements

		if (clt_parameters.show_map &&  (debugLevel > 0) && clt_parameters.rig.rig_mode_debug){
			(new ShowDoubleFloatArrays()).showArrays(
					disparity_bimap_infinity,
					tilesX,
					disparity_bimap_infinity[0].length/tilesX,
					true,
					quadCLT_main.image_name+"DISP_MAP-INFINITY",
					ImageDtt.BIDISPARITY_TITLES);

		}
		double [][] prev_bimap = null;
		double [] scale_bad = new double [trusted_infinity.length];
		for (int i = 0; i < scale_bad.length; i++) scale_bad[i] = 1.0;
		for (int nref = 0; nref < clt_parameters.rig.num_inf_refine; nref++) {
			// refine infinity using inter correlation
			double refine_tolerance = (nref > 0)?clt_parameters.rig.refine_tolerance:0.0;
			// only for infinity to prevent remaining 0.0?
			double [][] disparity_bimap_new =  refineRigSel(
					quadCLT_main,    // QuadCLT                        quadCLT_main,    // tiles should be set
					quadCLT_aux,     // QuadCLT                        quadCLT_aux,
					disparity_bimap_infinity, // double [][]                    src_bimap,       // current state of measurements (or null for new measurement)
					prev_bimap,      // double [][]                    prev_bimap, // previous state of measurements or null
					scale_bad,       // double []                      scale_bad,
					refine_inter,    // int                            refine_mode,     // 0 - by main, 1 - by aux, 2 - by inter
					false,           // boolean                        keep_inf,        // keep expected disparity 0.0 if it was so
					0.0,             // double                                         inf_disparity,
					0.0, // clt_parameters.rig.refine_min_strength , // double refine_min_strength, // do not refine weaker tiles
					refine_tolerance ,    // double refine_tolerance,    // do not refine if absolute disparity below
					trusted_infinity, // tile_list,       // ArrayList<Integer>             tile_list,       // or null
					num_new,         // int     []                                     num_new,
					clt_parameters,  // EyesisCorrectionParameters.CLTParameters clt_parameters,
					false,             //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
					0,                 // final int  // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
					// refine - no window offset
					true,           // boolean                                  no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
					threadsMax,      // final int                      threadsMax,      // maximal number of threads to launch
					updateStatus,    // final boolean                  updateStatus,
					debugLevel);     // final int                      debugLevel);
			prev_bimap = disparity_bimap_infinity;
			disparity_bimap_infinity = disparity_bimap_new;
			trusted_infinity = 	  getTrustedDisparity(
					quadCLT_main,                             // QuadCLT            quadCLT_main,  // tiles should be set
					quadCLT_aux,                              // QuadCLT            quadCLT_aux,
					true,                                     // boolean            use_individual,
					clt_parameters.rig.min_trusted_strength,  // double             min_combo_strength,    // check correlation strength combined for all 3 correlations
					clt_parameters.grow_disp_trust,           // double             max_trusted_disparity, // 4.0 -> change to rig_trust
					clt_parameters.rig.trusted_tolerance,     // double             trusted_tolerance,
					trusted_infinity, // null,                // boolean []         was_trusted,
					disparity_bimap_infinity );               // double [][]        bimap // current state of measurements
			if (debugLevel > -2) {
				System.out.println("enhanceByRig(): refined (infinity) "+num_new[0]+" tiles");
			}

			if (num_new[0] < clt_parameters.rig.min_new) break;
		}



		if (clt_parameters.show_map &&  (debugLevel > 0) && clt_parameters.rig.rig_mode_debug){

			for (int layer = 0; layer < disparity_bimap_infinity.length; layer ++) if (disparity_bimap_infinity[layer] != null){
				for (int nTile = 0; nTile < disparity_bimap_infinity[layer].length; nTile++) {
					if (!trusted_infinity[nTile]) disparity_bimap_infinity[layer][nTile] = Double.NaN;
				}
			}
			if (scale_bad!= null) {
				int num_bad = 0, num_trusted = 0;
				for (int nTile = 0; nTile < scale_bad.length; nTile++) {
					if (!trusted_infinity[nTile]) scale_bad[nTile] = Double.NaN;
					else {
						if (scale_bad[nTile] < 1.0) num_bad++;
						scale_bad[nTile] =  -scale_bad[nTile];
						num_trusted ++;

					}
				}
				System.out.println("num_trusted = "+num_trusted+", num_bad = "+num_bad);
				(new ShowDoubleFloatArrays()).showArrays(
						scale_bad,
						tilesX,
						disparity_bimap_infinity[0].length/tilesX,
						quadCLT_main.image_name+"-INFINITY-SCALE_BAD"+clt_parameters.disparity);

			}
			(new ShowDoubleFloatArrays()).showArrays(
					disparity_bimap_infinity,
					tilesX,
					disparity_bimap_infinity[0].length/tilesX,
					true,
					quadCLT_main.image_name+"-INFINITY-REFINED-TRUSTED"+clt_parameters.disparity,
					ImageDtt.BIDISPARITY_TITLES);
		}

		// Get DSI from the main camera
		quadCLT_main.tp.trimCLTPasses(false); // remove rig composite scan if any
		// last but not including any rig data
		CLTPass3d scan_last = quadCLT_main.tp.clt_3d_passes.get( quadCLT_main.tp.clt_3d_passes_size -1); // get last one
		double [][] disparity_bimap = setBimapFromCLTPass3d(
				scan_last,      // CLTPass3d                                scan,
				quadCLT_main,   // QuadCLT                                  quadCLT_main,  // tiles should be set
				quadCLT_aux,    // QuadCLT                                  quadCLT_aux,
				clt_parameters, // EyesisCorrectionParameters.CLTParameters clt_parameters,
				threadsMax,     // final int        threadsMax,  // maximal number of threads to launch
				updateStatus,   // final boolean    updateStatus,
				debugLevel);    // final int        debugLevel);

		//		  int [] num_new = new int[1];
		boolean [] trusted_near = 	  getTrustedDisparity(
				quadCLT_main,                                      // QuadCLT            quadCLT_main,  // tiles should be set
				quadCLT_aux,                                       // QuadCLT            quadCLT_aux,
				true,                                              // boolean            use_individual,
				0.5*clt_parameters.rig.min_trusted_strength,       // double             min_combo_strength,    // check correlation strength combined for all 3 correlations
				clt_parameters.grow_disp_trust,                    // double             max_trusted_disparity, // 4.0 -> change to rig_trust
				clt_parameters.rig.trusted_tolerance,              // double             trusted_tolerance,
				null,                                              // boolean []         was_trusted,
				disparity_bimap); // double [][]        bimap // current state of measurements
		if (clt_parameters.show_map &&  (debugLevel > 0) && clt_parameters.rig.rig_mode_debug){
			(new ShowDoubleFloatArrays()).showArrays(
					disparity_bimap,
					tilesX,
					disparity_bimap[0].length/tilesX,
					true,
					quadCLT_main.image_name+"DISP_MAP-NONINFINITY",
					ImageDtt.BIDISPARITY_TITLES);

		}

		for (int i = 0; i < scale_bad.length; i++) scale_bad[i] = 1.0;
		prev_bimap = null;
		for (int nref = 0; nref < clt_parameters.rig.num_near_refine; nref++) {
			// refine infinity using inter correlation
			double [][] disparity_bimap_new =  refineRigSel(
					quadCLT_main,    // QuadCLT                        quadCLT_main,    // tiles should be set
					quadCLT_aux,     // QuadCLT                        quadCLT_aux,
					disparity_bimap, // double [][]                    src_bimap,       // current state of measurements (or null for new measurement)
					prev_bimap,      // double [][]                    prev_bimap, // previous state of measurements or null
					scale_bad,       // double []                      scale_bad,
					refine_inter,    // int                            refine_mode,     // 0 - by main, 1 - by aux, 2 - by inter
					false,           // boolean                        keep_inf,        // keep expected disparity 0.0 if it was so
					0.0,             // double                                         inf_disparity,
					0.0, // clt_parameters.rig.refine_min_strength , // double refine_min_strength, // do not refine weaker tiles
					clt_parameters.rig.refine_tolerance ,    // double refine_tolerance,    // do not refine if absolute disparity below
					trusted_near,    // tile_list,       // ArrayList<Integer>             tile_list,       // or null
					num_new,         // int     []                                     num_new,
					clt_parameters,  // EyesisCorrectionParameters.CLTParameters clt_parameters,
					false,           //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
					0,                 // final int  // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
					// refine - no window offset
					true,           // boolean                                  no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
					threadsMax,      // final int                      threadsMax,      // maximal number of threads to launch
					updateStatus,    // final boolean                  updateStatus,
					debugLevel);     // final int                      debugLevel);
			prev_bimap = disparity_bimap;
			disparity_bimap = disparity_bimap_new;
			trusted_near = 	  getTrustedDisparity(
					quadCLT_main,                            // QuadCLT            quadCLT_main,  // tiles should be set
					quadCLT_aux,                             // QuadCLT            quadCLT_aux,
					true,                                     // boolean            use_individual,
					clt_parameters.rig.min_trusted_strength, // double             min_combo_strength,    // check correlation strength combined for all 3 correlations
					clt_parameters.grow_disp_trust,          // double             max_trusted_disparity, // 4.0 -> change to rig_trust
					clt_parameters.rig.trusted_tolerance,    // double             trusted_tolerance,
					trusted_near, // null,                   // boolean []         was_trusted,
					disparity_bimap );                       // double [][]        bimap // current state of measurements
			if (debugLevel > -2) {
				System.out.println("enhanceByRig(): refined (near) "+num_new[0]+" tiles");
			}
			if (num_new[0] < clt_parameters.rig.min_new) break;
		}

		if (clt_parameters.show_map &&  (debugLevel > 0) && clt_parameters.rig.rig_mode_debug){

			for (int layer = 0; layer < disparity_bimap.length; layer ++) if (disparity_bimap[layer] != null){
				for (int nTile = 0; nTile < disparity_bimap[layer].length; nTile++) {
					if (!trusted_near[nTile]) disparity_bimap[layer][nTile] = Double.NaN;
				}
			}
			if (scale_bad!= null) {
				int num_bad = 0, num_trusted = 0;
				for (int nTile = 0; nTile < scale_bad.length; nTile++) {
					if (!trusted_near[nTile]) scale_bad[nTile] = Double.NaN;
					else {
						if (scale_bad[nTile] < 1.0) num_bad++;
						scale_bad[nTile] =  -scale_bad[nTile];
						num_trusted ++;

					}
				}
				System.out.println("num_trusted = "+num_trusted+", num_bad = "+num_bad);
				(new ShowDoubleFloatArrays()).showArrays(
						scale_bad,
						tilesX,
						disparity_bimap[0].length/tilesX,
						quadCLT_main.image_name+"-NEAR-SCALE_BAD"+clt_parameters.disparity);

			}
			(new ShowDoubleFloatArrays()).showArrays(
					disparity_bimap,
					tilesX,
					disparity_bimap[0].length/tilesX,
					true,
					quadCLT_main.image_name+"-NEAR-REFINED-TRUSTED"+clt_parameters.disparity,
					ImageDtt.BIDISPARITY_TITLES);
		}
		// Combine infinity and non-infinity
		for (int nTile = 0; nTile < disparity_bimap[0].length; nTile++) {
			if (trusted_infinity[nTile] &&
					(!trusted_near[nTile] ||
							(disparity_bimap_infinity[ImageDtt.BI_STR_ALL_INDEX][nTile] > disparity_bimap[ImageDtt.BI_STR_ALL_INDEX][nTile]))) {
				for (int i = 0; i < disparity_bimap.length; i++) if (disparity_bimap_infinity[i]!=null) {
					disparity_bimap[i][nTile] = disparity_bimap_infinity[i][nTile];
				}
				trusted_near[nTile] = true;
			}
		}

		if (clt_parameters.show_map &&  (debugLevel > -3) && clt_parameters.rig.rig_mode_debug && !clt_parameters.batch_run){
			(new ShowDoubleFloatArrays()).showArrays(
					disparity_bimap,
					tilesX,
					disparity_bimap[0].length/tilesX,
					true,
					quadCLT_main.image_name+"DSI_ALL",
					ImageDtt.BIDISPARITY_TITLES);
			for (int layer = 0; layer < disparity_bimap.length; layer ++) if (disparity_bimap[layer] != null){
				for (int nTile = 0; nTile < disparity_bimap[layer].length; nTile++) {
					if (!trusted_near[nTile]) disparity_bimap[layer][nTile] = Double.NaN;
				}
			}
			for (int layer:ImageDtt.BIDISPARITY_STRENGTHS) if (disparity_bimap[layer] != null){
				for (int nTile = 0; nTile < disparity_bimap[layer].length; nTile++) {
					if (!trusted_near[nTile]) disparity_bimap[layer][nTile] = 0.0;
				}
			}
			(new ShowDoubleFloatArrays()).showArrays(
					disparity_bimap,
					tilesX,
					disparity_bimap[0].length/tilesX,
					true,
					quadCLT_main.image_name+"-DSI-ALL-TRUSTED",
					ImageDtt.BIDISPARITY_TITLES);
		}

		// just testing - copy and run in pole detection mode (TODO: Remove)
		if (debugLevel > -100) return disparity_bimap;

		// need to re-measure
		double [][] disparity_bimap_poles =  measureNewRigDisparity(
				quadCLT_main,    // QuadCLT                        quadCLT_main,    // tiles should be set
				quadCLT_aux,     // QuadCLT                        quadCLT_aux,
				disparity_bimap[ImageDtt.BI_TARGET_INDEX], // double []                                      disparity, // Double.NaN - skip, ohers - measure
				clt_parameters, // EyesisCorrectionParameters.CLTParameters       clt_parameters,
				true, // boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
				0,                 // final int  // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
				// any - in notch mode it is disabled
				clt_parameters.rig.no_int_x0,           // boolean                                  no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
				threadsMax,      // final int                      threadsMax,      // maximal number of threads to launch
				updateStatus,    // final boolean                  updateStatus,
				debugLevel); // +4);     // final int                      debugLevel);

		(new ShowDoubleFloatArrays()).showArrays(
				disparity_bimap_poles,
				tilesX,
				disparity_bimap[0].length/tilesX,
				true,
				quadCLT_main.image_name+"-INITIAL-POLES",
				ImageDtt.BIDISPARITY_TITLES);

		for (int i = 0; i < scale_bad.length; i++) scale_bad[i] = 1.0;
		prev_bimap = null;
		for (int nref = 0; nref < clt_parameters.rig.num_near_refine; nref++) {
			// refine infinity using inter correlation
			double [][] disparity_bimap_new =  refineRigSel(
					quadCLT_main,    // QuadCLT                        quadCLT_main,    // tiles should be set
					quadCLT_aux,     // QuadCLT                        quadCLT_aux,
					disparity_bimap_poles, // double [][]                    src_bimap,       // current state of measurements (or null for new measurement)
					prev_bimap,      // double [][]                    prev_bimap, // previous state of measurements or null
					scale_bad,       // double []                      scale_bad,
					refine_inter,    // int                            refine_mode,     // 0 - by main, 1 - by aux, 2 - by inter
					false,           // boolean                        keep_inf,        // keep expected disparity 0.0 if it was so
					0.0,             // double                                         inf_disparity,
					0.0, // clt_parameters.rig.refine_min_strength , // double refine_min_strength, // do not refine weaker tiles
					clt_parameters.rig.refine_tolerance ,    // double refine_tolerance,    // do not refine if absolute disparity below
					trusted_near,    // tile_list,       // ArrayList<Integer>             tile_list,       // or null
					num_new,         // int     []                                     num_new,
					clt_parameters,  // EyesisCorrectionParameters.CLTParameters clt_parameters,
					true,            //  final boolean                  notch_mode,      // use notch filter for inter-camera correlation to detect poles
					0,               // final int  // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
					true,            // boolean                                  no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
					threadsMax,      // final int                      threadsMax,      // maximal number of threads to launch
					updateStatus,    // final boolean                  updateStatus,
					debugLevel);     // final int                      debugLevel);
			prev_bimap = disparity_bimap_poles;
			disparity_bimap_poles = disparity_bimap_new;
			trusted_near = 	  getTrustedDisparity(
					quadCLT_main,                            // QuadCLT            quadCLT_main,  // tiles should be set
					quadCLT_aux,                             // QuadCLT            quadCLT_aux,
					true,                                    // boolean            use_individual,
					clt_parameters.rig.min_trusted_strength, // double             min_combo_strength,    // check correlation strength combined for all 3 correlations
					clt_parameters.grow_disp_trust,          // double             max_trusted_disparity, // 4.0 -> change to rig_trust
					clt_parameters.rig.trusted_tolerance,    // double             trusted_tolerance,
					trusted_near, // null,                   // boolean []         was_trusted,
					disparity_bimap_poles );                       // double [][]        bimap // current state of measurements
			if (debugLevel > -2) {
				System.out.println("enhanceByRig(): refined (poles) "+num_new[0]+" tiles");
			}
			if (num_new[0] < clt_parameters.rig.min_new) break;
		}

		if (clt_parameters.show_map &&  (debugLevel > -2) && clt_parameters.rig.rig_mode_debug){

			if (scale_bad!= null) {
				int num_bad = 0, num_trusted = 0;
				for (int nTile = 0; nTile < scale_bad.length; nTile++) {
					if (!trusted_near[nTile]) scale_bad[nTile] = Double.NaN;
					else {
						if (scale_bad[nTile] < 1.0) num_bad++;
						scale_bad[nTile] =  -scale_bad[nTile];
						num_trusted ++;

					}
				}
				System.out.println("num_trusted = "+num_trusted+", num_bad = "+num_bad);
				(new ShowDoubleFloatArrays()).showArrays(
						scale_bad,
						tilesX,
						disparity_bimap_poles[0].length/tilesX,
						quadCLT_main.image_name+"-NEAR-SCALE_BAD_POLES"+clt_parameters.disparity);

			}
			(new ShowDoubleFloatArrays()).showArrays(
					disparity_bimap_poles,
					tilesX,
					disparity_bimap_poles[0].length/tilesX,
					true,
					quadCLT_main.image_name+"-POLES-REFINED-TRUSTED",
					ImageDtt.BIDISPARITY_TITLES);
		}
		// END OF just testing - copy and run in pole detection mode (TODO: Remove)
		return disparity_bimap;

	}

	public BiScan rigInitialScan(
			QuadCLT            quadCLT_main,  // tiles should be set
			QuadCLT            quadCLT_aux,
			CLTParameters       clt_parameters,
			ColorProcParameters                            colorProcParameters, //
			ColorProcParameters                            colorProcParameters_aux, //
			final int                                      threadsMax,  // maximal number of threads to launch
			final boolean                                  updateStatus,
			final int                                      debugLevel) //  throws Exception
	{
		if (getBiScan(0) != null) {
			return getBiScan(0);
		};

		final int tilesX = quadCLT_main.tp.getTilesX();
		final int tilesY = quadCLT_main.tp.getTilesY();
		BiCamDSI biCamDSI = new BiCamDSI( tilesX, tilesY,threadsMax);
		System.out.println("rigInitialScan()");
		double [][] disparity_bimap = prepareRefineExistingDSI(
				quadCLT_main,  // QuadCLT            quadCLT_main,  // tiles should be set
				quadCLT_aux,   // QuadCLT            quadCLT_aux,
				clt_parameters, // EyesisCorrectionParameters.CLTParameters       clt_parameters,
				colorProcParameters,          //  ColorProcParameters                       colorProcParameters, //
				colorProcParameters_aux,          //  ColorProcParameters                       colorProcParameters, //
				threadsMax,        // final int           threadsMax,      // maximal number of threads to launch
				updateStatus,      // final boolean       updateStatus,
				debugLevel);       // final int           debugLevel);
		if (disparity_bimap == null) {
			String msg = "Failed to get (and refine) initial rig DSI from the existing data";
			System.out.println(msg);
			IJ.showMessage("ERROR",msg);
			return null;
		}

		biCamDSI.addBiScan(disparity_bimap, BiScan.BISCAN_SINGLECORR);
		biCamDSI.getBiScan(0).calcTrusted(   // finds strong trusted and validates week ones if they fit planes
				clt_parameters.rig.pf_trusted_strength, // final double     trusted_strength, // trusted correlation strength
				clt_parameters.rig.pf_strength_rfloor,  // final double     strength_rfloor,   // strength floor - relative to trusted
				clt_parameters.rig.pf_cond_rtrusted,    // final double     cond_rtrusted,     // minimal strength to consider - fraction of trusted
				clt_parameters.rig.pf_strength_pow,     // final double     strength_pow,      // raise strength-floor to this power
				clt_parameters.rig.pf_smpl_radius,      // final int        smpl_radius,
				clt_parameters.rig.pf_smpl_num,         // final int        smpl_num,   //         = 3;      // Number after removing worst (should be >1)
				clt_parameters.rig.pf_smpl_fract,       // final double     smpl_fract, // Number of friends among all neighbors
				clt_parameters.rig.pf_max_adiff,        // final double     max_adiff,  // Maximal absolute difference betweenthe center tile and friends
				clt_parameters.rig.pf_max_rdiff,        // final double     max_rdiff, //  Maximal relative difference between the center tile and friends
				clt_parameters.rig.pf_max_atilt,        // final double     max_atilt, //  = 2.0; // pix per tile
				clt_parameters.rig.pf_max_rtilt,        // final double     max_rtilt, //  = 0.2; // (pix / disparity) per tile
				clt_parameters.rig.pf_smpl_arms,        // final double     smpl_arms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
				clt_parameters.rig.pf_smpl_rrms,        // final double     smpl_rrms,        //      = 0.005;  // Maximal RMS/disparity in addition to smplRms
				clt_parameters.rig.pf_damp_tilt,        // final double     damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
				clt_parameters.rig.pf_rwsigma,          // 						final double     rwsigma,           //  = 0.7; // influence of far neighbors diminish as a Gaussian with this sigma
				clt_parameters.tileX,                   // final int        dbg_x,
				clt_parameters.tileY,                   // final int        dbg_y,
				debugLevel);                            // final int        debugLevel);


		if (clt_parameters.show_map &&  (debugLevel > 0) && clt_parameters.rig.rig_mode_debug){
			biCamDSI.getLastBiScan(BiScan.BISCAN_SINGLECORR).showScan(
					quadCLT_main.image_name+"-BISCAN_initial");
		}
		biCamDSI_persistent = biCamDSI;
		return getBiScan(0);
	}



	public double [][] groundTruthByRigPlanes(
			QuadCLT            quadCLT_main,  // tiles should be set
			QuadCLT            quadCLT_aux,
			CLTParameters       clt_parameters,
			ColorProcParameters                            colorProcParameters, //
			ColorProcParameters                            colorProcParameters_aux, //
			final int                                      threadsMax,  // maximal number of threads to launch
			final boolean                                  updateStatus,
			final int                                      debugLevel) //  throws Exception
	{
		if (getBiScan(0) != null) {
			if (getBiScan(1) != null) {
				System.out.println("Expected just a single BiScan here, got "+(biCamDSI_persistent.biScans.size())+". Trimming");
				while (biCamDSI_persistent.biScans.size()>1) {
					biCamDSI_persistent.biScans.remove(1);
				}
				Runtime.getRuntime().gc();
				System.out.println("--- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");


			}
		} else { // create a new one
			if (rigInitialScan(				quadCLT_main,  // QuadCLT            quadCLT_main,  // tiles should be set
					quadCLT_aux,   // QuadCLT            quadCLT_aux,
					clt_parameters, // EyesisCorrectionParameters.CLTParameters       clt_parameters,
					colorProcParameters,          //  ColorProcParameters                       colorProcParameters, //
					colorProcParameters_aux,          //  ColorProcParameters                       colorProcParameters, //

					threadsMax,        // final int           threadsMax,      // maximal number of threads to launch
					updateStatus,      // final boolean       updateStatus,
					debugLevel)== null)       // final int           debugLevel);
			{
				return null; // Failed to get
			}
		}

		final int num_tries_strongest_by_fittest = 5;
		final int num_full_cycles =       clt_parameters.rig.pf_en_trim_fg? 3 : 1;  // Number of full low-texture cycles that include growing flat LT and trimmin weak FG over BG
		//		  final int num_cross_gaps_cycles = 12+ (num_simple_expand_cysles * dxy.length); // maximal number of adding new tiles cycles while "crossing the gaps)
		final int min_cross_gaps_new =    20; // minimal number of the new added tiles
		final int refine_inter = 2; // 3; // 3 - dx, 2 - disparity
		final int tilesX = quadCLT_main.tp.getTilesX();
		final int tilesY = quadCLT_main.tp.getTilesY();


		/*
		BiCamDSI biCamDSI = new BiCamDSI( tilesX, tilesY,threadsMax);
		System.out.println("groundTruthByRigPlanes()");


		double [][] disparity_bimap = prepareRefineExistingDSI(
				quadCLT_main,  // QuadCLT            quadCLT_main,  // tiles should be set
				quadCLT_aux,   // QuadCLT            quadCLT_aux,
				clt_parameters, // EyesisCorrectionParameters.CLTParameters       clt_parameters,
				threadsMax,        // final int           threadsMax,      // maximal number of threads to launch
				updateStatus,      // final boolean       updateStatus,
				debugLevel);       // final int           debugLevel);
		if (disparity_bimap == null) {
			String msg = "Failed to get (and refine) initial rig DSI from the existing data";
			System.out.println(msg);
			IJ.showMessage("ERROR",msg);
			return null;
		}

		biCamDSI.addBiScan(disparity_bimap, BiScan.BISCAN_SINGLECORR);
		if (clt_parameters.show_map &&  (debugLevel > 0) && clt_parameters.rig.rig_mode_debug){
			biCamDSI.getLastBiScan(BiScan.BISCAN_SINGLECORR).showScan(
					quadCLT_main.image_name+"-BISCAN_initial");
		}
 		*/
		BiCamDSI biCamDSI = biCamDSI_persistent;



		int [][] dxy = {{0, -1},{0,1},{-1,0},{1,0}};
		int num_simple_expand_cysles = 8;
		final int num_cross_gaps_cycles = 12 + (num_simple_expand_cysles * dxy.length); // maximal number of adding new tiles cycles while "crossing the gaps)
		for (int num_fcycle = 0; num_fcycle < num_full_cycles; num_fcycle++) {
			// Grow tiles, cross gaps (do not trim yet
			//
			int [] num_simple_added = new int[dxy.length];
			for (int i = 0; i < num_simple_added.length; i++) {
				num_simple_added[i] = -1;
			}
			for (int num_cycle = 0; num_cycle < num_cross_gaps_cycles; num_cycle++) {
				BiScan last_scan = biCamDSI.getLastBiScan(BiScan.BISCAN_SINGLECORR);
				if (clt_parameters.rig.lonefg_disp_incr > 0.0) {
					int removed_lone = last_scan.trimWeakLoneFG(
							clt_parameters.rig.pf_trusted_strength, // final double    trusted_strength, // trusted correlation strength
							clt_parameters.rig.lonefg_rstrength,  // final double     min_rstrength,   // strength floor - relative to trusted
							clt_parameters.rig.lonefg_disp_incr ,//final double     max_disp_inc,
							clt_parameters.tileX,                   // final int        dbg_x,
							clt_parameters.tileY,                   // final int        dbg_y,
							debugLevel+2);                            // final int        debugLevel);
					System.out.println("trimWeakLoneFG() -> "+removed_lone);
				}
				int [] trusted_stats = last_scan.calcTrusted(   // finds strong trusted and validates week ones if they fit planes
						clt_parameters.rig.pf_trusted_strength, // final double     trusted_strength, // trusted correlation strength
						clt_parameters.rig.pf_strength_rfloor,  // final double     strength_rfloor,   // strength floor - relative to trusted
						clt_parameters.rig.pf_cond_rtrusted,    // final double     cond_rtrusted,     // minimal strength to consider - fraction of trusted
						clt_parameters.rig.pf_strength_pow,     // final double     strength_pow,      // raise strength-floor to this power
						clt_parameters.rig.pf_smpl_radius,      // final int        smpl_radius,
						clt_parameters.rig.pf_smpl_num,         // final int        smpl_num,   //         = 3;      // Number after removing worst (should be >1)
						clt_parameters.rig.pf_smpl_fract,       // final double     smpl_fract, // Number of friends among all neighbors
						clt_parameters.rig.pf_max_adiff,        // final double     max_adiff,  // Maximal absolute difference betweenthe center tile and friends
						clt_parameters.rig.pf_max_rdiff,        // final double     max_rdiff, //  Maximal relative difference between the center tile and friends
						clt_parameters.rig.pf_max_atilt,        // final double     max_atilt, //  = 2.0; // pix per tile
						clt_parameters.rig.pf_max_rtilt,        // final double     max_rtilt, //  = 0.2; // (pix / disparity) per tile
						clt_parameters.rig.pf_smpl_arms,        // final double     smpl_arms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
						clt_parameters.rig.pf_smpl_rrms,        // final double     smpl_rrms,        //      = 0.005;  // Maximal RMS/disparity in addition to smplRms
						clt_parameters.rig.pf_damp_tilt,        // final double     damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
						clt_parameters.rig.pf_rwsigma,          // 						final double     rwsigma,           //  = 0.7; // influence of far neighbors diminish as a Gaussian with this sigma
						clt_parameters.tileX,                   // final int        dbg_x,
						clt_parameters.tileY,                   // final int        dbg_y,
						debugLevel);                            // final int        debugLevel);
				if (clt_parameters.show_map &&  (debugLevel > 0) && clt_parameters.rig.rig_mode_debug){
					biCamDSI.getLastBiScan(BiScan.BISCAN_SINGLECORR).showScan(
							quadCLT_main.image_name+"-BISCAN_TRUSTED"+num_fcycle+"-"+num_cycle);
				}

				/*
				 * @return array of 3 numbers: number of trusted strong tiles, number of additional trusted by plane fitting, and number of all
				 * somewhat strong tiles
				 */
				if (debugLevel > -4) {
					System.out.println("groundTruthByRigPlanes() grow pass "+num_cycle+" of "+ num_cross_gaps_cycles+
							" strong trusted: "+trusted_stats[0]+ " neib trusted: "+trusted_stats[1]+" weak trusted: " + trusted_stats[2]);
				}
				int num_added_tiles =0;
				if ((num_cycle < num_simple_expand_cysles * dxy.length) && (num_cycle >= dxy.length)) {
					boolean all_last_zeros = true;
					for (int i =0 ; i <=dxy.length; i++) {
						if (num_simple_added[(num_cycle-1-i) % dxy.length] != 0) {
							all_last_zeros = false;
							break;
						}
					}
					if (all_last_zeros) {
						if (debugLevel > -4) {
							System.out.println("==== groundTruthByRigPlanes() grow pass "+num_cycle+" of "+ num_cross_gaps_cycles+
									" all "+(num_simple_added.length)+" direction expansions added 0 tiles, proceeding to low texture expansion ====");
						}
						num_cycle = num_simple_expand_cysles * dxy.length;
					}
				}

				if (num_cycle < num_simple_expand_cysles * dxy.length) {
					// simple duplicating one step in 4 directions
					num_added_tiles = last_scan.suggestNewScan(
							dxy[num_cycle % dxy.length], // final int []     dxy,               //up,down,right,left
							clt_parameters.rig.pf_discard_cond,     // final boolean    discard_cond,      // consider conditionally trusted tiles (not promoted to trusted) as empty
							clt_parameters.rig.pf_discard_weak,     // final boolean    discard_weak,      // consider conditionally trusted tiles (not promoted to trusted) as empty
							clt_parameters.rig.pf_discard_strong,   // final boolean    discard_strong,      // consider conditionally trusted tiles (not promoted to trusted) as empty
							clt_parameters.rig.pf_new_diff,         // final double     new_diff,            // minimal difference between the new suggested and the already tried/measured one
							true,                                   // final boolean    remove_all_tried,  // remove from suggested - not only disabled, but all tried
							clt_parameters.tileX,                   // final int        dbg_x,
							clt_parameters.tileY,                   // final int        dbg_y,
							debugLevel);                            // final int        debugLevel);
					//TODO: add expanding FG over existing BG. Use "Strong enough" for FG to beat BG. Maybe expand by multiple steps?
					num_simple_added[num_cycle % dxy.length] = num_added_tiles;
				} else {

					// suggest new disparities, using plane surfaces (extending around that may cause false surfaces)
					/*
					 * Suggest disparities to try for the tiles in poorly textured areas by fitting planes in DSI
					 * calcTrusted should be called before to set up trusted/cond_trusted tiles
					 * suggested tiles will be compared against and made sure they differ by more than a specified margin
					 * 1) current measured (refined) disparity value
					 * 2) target disparity that lead to the current measurement after refinement
					 * 3) any other disable measurement
					 * 4) any target disparity that lead to the disabled measurement
					 * @return number of new tiles to measure in the  array of suggested disparities - Double.NaN - nothing suggested
					 *  for the tile. May need additional filtering to avoid suggested already tried disparities
					 */

					num_added_tiles = last_scan.suggestNewScan(
							null,                                   // final boolean [] area_of_interest,
							null,                                   // final double [][] disparityStrength,
							clt_parameters.rig.pf_trusted_strength, // final double     trusted_strength, // trusted correlation strength
							clt_parameters.rig.pf_strength_rfloor,  // final double     strength_rfloor,   // strength floor - relative to trusted
							clt_parameters.rig.pf_discard_cond,     // final boolean    discard_cond,      // consider conditionally trusted tiles (not promoted to trusted) as empty
							clt_parameters.rig.pf_discard_weak,     // final boolean    discard_weak,      // consider conditionally trusted tiles (not promoted to trusted) as empty
							clt_parameters.rig.pf_discard_strong,   // final boolean    discard_strong,      // consider conditionally trusted tiles (not promoted to trusted) as empty
							clt_parameters.rig.pf_strength_pow,     // final double     strength_pow,      // raise strength-floor to this power
							clt_parameters.rig.pf_smpl_radius,      // final int        smpl_radius,
							clt_parameters.rig.pf_smpl_num,         // final int        smpl_num,   //         = 3;      // Number after removing worst (should be >1)
							clt_parameters.rig.pf_smpl_fract,       // final double     smpl_fract, // Number of friends among all neighbors
							clt_parameters.rig.pf_smpl_num_narrow,  // final int        smpl_num_narrow,   //         = 3;      // Number after removing worst (should be >1)
							clt_parameters.rig.pf_max_adiff,        // final double     max_adiff,  // Maximal absolute difference betweenthe center tile and friends
							clt_parameters.rig.pf_max_rdiff,        // final double     max_rdiff, //  Maximal relative difference between the center tile and friends
							clt_parameters.rig.pf_max_atilt,        // final double     max_atilt, //  = 2.0; // pix per tile
							clt_parameters.rig.pf_max_rtilt,        // final double     max_rtilt, //  = 0.2; // (pix / disparity) per tile
							clt_parameters.rig.pf_smpl_arms,        // final double     smpl_arms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
							clt_parameters.rig.pf_smpl_rrms,        // final double     smpl_rrms,        //      = 0.005;  // Maximal RMS/disparity in addition to smplRms
							clt_parameters.rig.pf_damp_tilt,        // final double     damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
							clt_parameters.rig.pf_rwsigma,          // final double     rwsigma,           //  = 0.7; // influence of far neighbors diminish as a Gaussian with this sigma
							clt_parameters.rig.pf_rwsigma_narrow,   // final double     rwsigma_narrow,    //  = used to determine initial tilt
							clt_parameters.rig.pf_new_diff,         // final double     new_diff,            // minimal difference between the new suggested and the already tried/measured one
							true,                                   // final boolean    remove_all_tried,  // remove from suggested - not only disabled, but all tried
							0.0,                                    // final double     center_weight,     // use center tile too (0.0 - do not use)
							clt_parameters.rig.pf_use_alt,          // final boolean    use_alt,           // use tiles from other scans if they fit better
							clt_parameters.rig.pf_goal_fraction_rms,// final double     goal_fraction_rms, // Try to make rms to be this fraction of maximal acceptable by removing outliers
							clt_parameters.rig.pf_boost_low_density,// NOT USED HERE, MAY BE 0,                                    // final double     boost_low_density, // 0 - strength is proportional to 1/density, 1.0 - same as remaining tiles
							null,                                   // final double []  smooth_strength,   // optionally fill strength array when used for smoothing DSI
							0,                                      // final int        fourq_min,         // each of the 4 corners should have at least this number of tiles.
							0,                                      // final int        fourq_gap,         // symmetrical vertical and horizontal center areas that do not belong to any corner
							clt_parameters.tileX,                   // final int        dbg_x,
							clt_parameters.tileY,                   // final int        dbg_y,
							debugLevel);                            // final int        debugLevel);
				}
				if (debugLevel > -4) {
					System.out.println("groundTruthByRigPlanes() full cycle = "+num_fcycle+", grow pass "+num_cycle+" of "+ num_cross_gaps_cycles+
							" suggestNewScan() -> "+num_added_tiles);
				}
				//num_cycle < num_cross_gaps_cycles;
				//				  boolean last_cycle = (num_added_tiles < min_cross_gaps_new) || (num_cycle >= (num_cross_gaps_cycles-1));
				boolean last_cycle = ((num_added_tiles < min_cross_gaps_new) && (num_cycle >= num_simple_expand_cysles * dxy.length)) || (num_cycle >= (num_cross_gaps_cycles-1));
				if (clt_parameters.show_map &&  (debugLevel > -2) && clt_parameters.rig.rig_mode_debug){
					if (last_cycle) //  || (num_cycle < 2 * dxy.length))
						biCamDSI.getLastBiScan(BiScan.BISCAN_SINGLECORR).showScan(
								quadCLT_main.image_name+"-BISCAN_SUGGESTED"+num_fcycle+"-"+num_cycle);
				}

				if (last_cycle && clt_parameters.rig.pf_en_trim_fg) { // last cycle and trimming enabled
					if (debugLevel > -4) {
						//						  System.out.println("groundTruthByRigPlanes(): num_added_tiles= "+num_added_tiles+" > "+min_cross_gaps_new+", done growing over gaps");
						System.out.println("groundTruthByRigPlanes(): that was the last growing over gaps cycle, performing trimming hanging weak FG over BG");
					}

					/*
					 * Disable low-textured tiles are not between strong tiles, but on one side of it.
					 * This method relies on the assumption that FG edge should have strong correlation, so it tries multiple directions
					 * from the weak (not trusted strong) tiles and trims tiles that either do not have anything in that direction or have
					 * farther tiles.
					 * Trimming(disabling) weak (trusted but not strong_trusted) tiles if on any one side:
					 *   a) there are no same plane or closer tiles
					 *   b) there are no in-plane or closer strong tiles, but there are some (strong or any?) farther tiles
					 * repeat while more are trimmed
					 * maybe, if there are both strong in-plane and far - see which are closer
					 */
					int num_trimmed = last_scan.trimWeakFG(
							clt_parameters.rig.pf_trusted_strength, // final double     trusted_strength, // trusted correlation strength
							clt_parameters.rig.pf_strength_rfloor,  // final double     strength_rfloor,   // strength floor - relative to trusted
							clt_parameters.rig.pf_cond_rtrusted,    // final double     cond_rtrusted,     // minimal strength to consider - fraction of trusted
							clt_parameters.rig.pf_strength_pow,     // final double     strength_pow,      // raise strength-floor to this power
							clt_parameters.rig.pf_smpl_radius,      // final int        smpl_radius,
							clt_parameters.rig.pf_smpl_num,         // final int        smpl_num,   //         = 3;      // Number after removing worst (should be >1)
							clt_parameters.rig.pf_smpl_fract,       // final double     smpl_fract, // Number of friends among all neighbors
							clt_parameters.rig.pf_max_adiff,        // final double     max_adiff,  // Maximal absolute difference betweenthe center tile and friends
							clt_parameters.rig.pf_max_rdiff,        // final double     max_rdiff, //  Maximal relative difference between the center tile and friends
							clt_parameters.rig.pf_max_atilt,        // final double     max_atilt, //  = 2.0; // pix per tile
							clt_parameters.rig.pf_max_rtilt,        // final double     max_rtilt, //  = 0.2; // (pix / disparity) per tile
							clt_parameters.rig.pf_smpl_arms,        // final double     smpl_arms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
							clt_parameters.rig.pf_smpl_rrms,        // final double     smpl_rrms,        //      = 0.005;  // Maximal RMS/disparity in addition to smplRms
							clt_parameters.rig.pf_damp_tilt,        // final double     damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
							clt_parameters.rig.pf_rwsigma,          // final double     rwsigma,           //  = 0.7; // influence of far neighbors diminish as a Gaussian with this sigma

							clt_parameters.rig.pf_atolerance,       // final double     atolerance,  // When deciding closer/farther
							clt_parameters.rig.pf_rtolerance,       // final double     rtolerance,  // same, scaled with disparity
							clt_parameters.rig.pf_num_dirs,         // final int        num_dirs,    // number of directions to try
							clt_parameters.rig.pf_blind_dist,       // final double     blind_dist,  // analyze only tiles farther than this in the selected direction
							clt_parameters.rig.pf_strong_only_far,  // final boolean    strong_only_far, // in variant b) only compare with strong far
							clt_parameters.rig.pf_num_strong_far,   // final int        num_strong_far,    // number of directions to try
							clt_parameters.rig.pf_num_weak_far,     // final int        num_weak_far,     // number of directions to try
							clt_parameters.tileX,                   // final int        dbg_x,
							clt_parameters.tileY,                   // final int        dbg_y,
							debugLevel);                            // final int        debugLevel);
					if (debugLevel > -4) {
						System.out.println("groundTruthByRigPlanes(): full cycle="+num_fcycle+" num_trimmed= "+num_trimmed+" tiles");
					}
					if (clt_parameters.show_map &&  (debugLevel > 0) && clt_parameters.rig.rig_mode_debug){
						biCamDSI.getLastBiScan(BiScan.BISCAN_SINGLECORR).showScan(
								quadCLT_main.image_name+"-BISCAN_TRIMMED"+num_fcycle+"-"+num_cycle);
					}

					// suggest again, after trimming
					int num_added_tiles_trimmed = last_scan.suggestNewScan(
							null,                                   // final boolean [] area_of_interest,
							null,                                   // final double [][] disparityStrength,
							clt_parameters.rig.pf_trusted_strength, // final double     trusted_strength, // trusted correlation strength
							clt_parameters.rig.pf_strength_rfloor,  // final double     strength_rfloor,   // strength floor - relative to trusted
							clt_parameters.rig.pf_discard_cond,     // final boolean    discard_cond,      // consider conditionally trusted tiles (not promoted to trusted) as empty
							clt_parameters.rig.pf_discard_weak,     // final boolean    discard_weak,      // consider conditionally trusted tiles (not promoted to trusted) as empty
							clt_parameters.rig.pf_discard_strong,   // final boolean    discard_strong,      // consider conditionally trusted tiles (not promoted to trusted) as empty
							clt_parameters.rig.pf_strength_pow,     // final double     strength_pow,      // raise strength-floor to this power
							clt_parameters.rig.pf_smpl_radius,      // final int        smpl_radius,
							clt_parameters.rig.pf_smpl_num,         // final int        smpl_num,   //         = 3;      // Number after removing worst (should be >1)
							clt_parameters.rig.pf_smpl_fract,       // final double     smpl_fract, // Number of friends among all neighbors
							clt_parameters.rig.pf_smpl_num_narrow,  // final int        smpl_num_narrow,   //         = 3;      // Number after removing worst (should be >1)
							clt_parameters.rig.pf_max_adiff,        // final double     max_adiff,  // Maximal absolute difference betweenthe center tile and friends
							clt_parameters.rig.pf_max_rdiff,        // final double     max_rdiff, //  Maximal relative difference between the center tile and friends
							clt_parameters.rig.pf_max_atilt,        // final double     max_atilt, //  = 2.0; // pix per tile
							clt_parameters.rig.pf_max_rtilt,        // final double     max_rtilt, //  = 0.2; // (pix / disparity) per tile
							clt_parameters.rig.pf_smpl_arms,        // final double     smpl_arms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
							clt_parameters.rig.pf_smpl_rrms,        // final double     smpl_rrms,        //      = 0.005;  // Maximal RMS/disparity in addition to smplRms
							clt_parameters.rig.pf_damp_tilt,        // final double     damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
							clt_parameters.rig.pf_rwsigma,          // final double     rwsigma,           //  = 0.7; // influence of far neighbors diminish as a Gaussian with this sigma
							clt_parameters.rig.pf_rwsigma_narrow,   // final double     rwsigma_narrow,    //  = used to determine initial tilt
							clt_parameters.rig.pf_new_diff,         // final double     new_diff,            // minimal difference between the new suggested and the already tried/measured one
							true,                                   // final boolean    remove_all_tried,  // remove from suggested - not only disabled, but all tried
							0.0,                                    // final double     center_weight,     // use center tile too (0.0 - do not use)
							clt_parameters.rig.pf_use_alt,          // final boolean    use_alt,           // use tiles from other scans if they fit better
							clt_parameters.rig.pf_goal_fraction_rms,// final double     goal_fraction_rms, // Try to make rms to be this fraction of maximal acceptable by removing outliers
							clt_parameters.rig.pf_boost_low_density,// NOT USED HERE, MAY BE 0,                                    // final double     boost_low_density, // 0 - strength is proportional to 1/density, 1.0 - same as remaining tiles
							null,                                   // final double []  smooth_strength,   // optionally fill strength array when used for smoothing DSI
							0,                                      // final int        fourq_min,         // each of the 4 corners should have at least this number of tiles.
							0,                                      // final int        fourq_gap,         // symmetrical vertical and horizontal center areas that do not belong to any corner
							clt_parameters.tileX,                   // final int        dbg_x,
							clt_parameters.tileY,                   // final int        dbg_y,
							debugLevel);                            // final int        debugLevel);
					if (debugLevel > -4) {
						System.out.println("groundTruthByRigPlanes() full cycle = "+num_fcycle+", grow pass "+num_cycle+" of "+ num_cross_gaps_cycles+
								" suggestNewScan() -> "+num_added_tiles_trimmed+"( after trimming)");
					}
					if (clt_parameters.show_map &&  (debugLevel > -2) && clt_parameters.rig.rig_mode_debug){
						biCamDSI.getLastBiScan(BiScan.BISCAN_SINGLECORR).showScan(
								quadCLT_main.image_name+"-BISCAN_TRIMMED_SUGGESTED"+num_fcycle+"-"+num_cycle);
					}

					//					  break; // too few added before trimmimng or number of steps exceeded limit
				} // if (last_cycle && clt_parameters.rig.pf_en_trim_fg) { // last cycle and trimming enabled
				// measure and refine
				double [] target_disparity = biCamDSI.getTargetDisparity(-1); // get last

				// Measure provided tiles (break after, if it was the last cycle)

				double [][] disparity_bimap = measureNewRigDisparity(
						quadCLT_main,      // QuadCLT             quadCLT_main,    // tiles should be set
						quadCLT_aux,       // QuadCLT             quadCLT_aux,
						target_disparity,  // double []                                      disparity, // Double.NaN - skip, ohers - measure
						clt_parameters,    // EyesisCorrectionParameters.CLTParameters       clt_parameters,
						false,             //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
						0,                 // final int  // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
						// first measurement - use default value:
						clt_parameters.rig.no_int_x0,           // boolean                                  no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
						threadsMax,        // final int           threadsMax,      // maximal number of threads to launch
						updateStatus,      // final boolean       updateStatus,
						debugLevel);       // final int           debugLevel);

				// refine measurements
				int [] num_new = new int[1];
				// at least for small separation FG/BG individual cameras may not provide trusted results - ignore them only use rig
				boolean [] trusted_measurements = 	  getTrustedDisparity(
						quadCLT_main,                            // QuadCLT            quadCLT_main,  // tiles should be set
						quadCLT_aux,                             // QuadCLT            quadCLT_aux,
						//						  false,                                   // boolean            use_individual,
						true,                                   // boolean            use_individual,
						0.8*clt_parameters.rig.min_trusted_strength, // double             min_combo_strength,    // check correlation strength combined for all 3 correlations
						clt_parameters.grow_disp_trust,          // double             max_trusted_disparity, // 4.0 -> change to rig_trust
						clt_parameters.rig.trusted_tolerance,    // double             trusted_tolerance,
						null,                                    // boolean []         was_trusted,
						disparity_bimap);              // double [][]        bimap // current state of measurements

				if (clt_parameters.show_map &&  (debugLevel > 0) && clt_parameters.rig.rig_mode_debug){
					(new ShowDoubleFloatArrays()).showArrays(
							disparity_bimap,
							tilesX,
							tilesY,
							true,
							quadCLT_main.image_name+"-gaps_cycle"+num_cycle,
							ImageDtt.BIDISPARITY_TITLES);

				}
				double [][] prev_bimap = null;
				double [] scale_bad = new double [trusted_measurements.length];
				for (int i = 0; i < scale_bad.length; i++) scale_bad[i] = 1.0;
				for (int nref = 0; nref < clt_parameters.rig.num_inf_refine; nref++) {
					// refine infinity using inter correlation
					double [][] disparity_bimap_new =  refineRigSel(
							quadCLT_main,    // QuadCLT                        quadCLT_main,    // tiles should be set
							quadCLT_aux,     // QuadCLT                        quadCLT_aux,
							disparity_bimap, // double [][]                    src_bimap,       // current state of measurements (or null for new measurement)
							prev_bimap,      // double [][]                    prev_bimap, // previous state of measurements or null
							scale_bad,       // double []                      scale_bad,
							refine_inter,    // int                            refine_mode,     // 0 - by main, 1 - by aux, 2 - by inter
							false,           // boolean                        keep_inf,        // keep expected disparity 0.0 if it was so
							0.0,             // double                                         inf_disparity,
							0.0, // clt_parameters.rig.refine_min_strength , // double refine_min_strength, // do not refine weaker tiles
							clt_parameters.rig.refine_tolerance ,    // double refine_tolerance,    // do not refine if absolute disparity below
							trusted_measurements, // tile_list,       // ArrayList<Integer>             tile_list,       // or null
							num_new,         // int     []                                     num_new,
							clt_parameters,  // EyesisCorrectionParameters.CLTParameters clt_parameters,
							false,           //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
							0,               // final int  // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
							// disable window preset in refine mode
							true,            // boolean                                  no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
							threadsMax,      // final int                      threadsMax,      // maximal number of threads to launch
							updateStatus,    // final boolean                  updateStatus,
							debugLevel);     // final int                      debugLevel);
					prev_bimap = disparity_bimap;
					disparity_bimap = disparity_bimap_new;
					trusted_measurements = 	  getTrustedDisparityInter(
							0.0, // clt_parameters.rig.lt_trusted_strength*clt_parameters.rig.lt_need_friends, // double             min_inter_strength,    // check correlation strength combined for all 3 correlations
							clt_parameters.grow_disp_trust,           // double             max_trusted_disparity,
							trusted_measurements,                     // boolean []         was_trusted,
							disparity_bimap );                        // double [][]        bimap // current state of measurements

					if (debugLevel > -2) {
						System.out.println("groundTruthByRigPlanes(): cycle="+num_cycle+", refinement step="+nref+" num_new= "+num_new[0]+" tiles");
					}
					if (num_new[0] < clt_parameters.rig.pf_min_new) {
						break;
					}
				}


				//FIXME: 	show only for the last_cycle
				if (clt_parameters.show_map &&  (debugLevel > 0) && clt_parameters.rig.rig_mode_debug){
					double [][] dbg_img = new double[disparity_bimap.length][];
					for (int layer = 0; layer < disparity_bimap.length; layer ++) if (disparity_bimap[layer] != null){
						dbg_img [layer]= disparity_bimap[layer].clone();
						for (int nTile = 0; nTile < disparity_bimap[layer].length; nTile++) {
							if (!trusted_measurements[nTile]) dbg_img[layer][nTile] = Double.NaN;
						}
					}
					if ((scale_bad!= null) && (debugLevel > 0)){
						int num_bad = 0, num_trusted = 0;
						for (int nTile = 0; nTile < scale_bad.length; nTile++) {
							if (!trusted_measurements[nTile]) scale_bad[nTile] = Double.NaN;
							else {
								if (scale_bad[nTile] < 1.0) num_bad++;
								scale_bad[nTile] =  -scale_bad[nTile];
								num_trusted ++;

							}
						}
						System.out.println("num_trusted = "+num_trusted+", num_bad = "+num_bad);
						(new ShowDoubleFloatArrays()).showArrays(
								scale_bad,
								tilesX,
								tilesY,
								quadCLT_main.image_name+"-gaps_cycle"+num_cycle+"-scale_bad");

					}
					(new ShowDoubleFloatArrays()).showArrays(
							dbg_img,
							tilesX,
							tilesY,
							true,
							quadCLT_main.image_name+"-gaps_cycle"+num_cycle+"-REFINED_TRUSTED",
							ImageDtt.BIDISPARITY_TITLES);
				}

				// add refined data
				biCamDSI.addBiScan(disparity_bimap, BiScan.BISCAN_SINGLECORR);
				// find strongest if FG over BG - not the strongest, but either strongest or just strong and nearer
				biCamDSI.getLastBiScan(BiScan.BISCAN_SINGLECORR).copyLastStrongestEnabled(
						clt_parameters.rig.pf_last_priority); // final boolean last_priority)

				double afloor = clt_parameters.rig.pf_trusted_strength * clt_parameters.rig.pf_strength_rfloor;

				// replace strongest by fittest
				for (int nfit = 0; nfit < num_tries_strongest_by_fittest; nfit++) {
					int num_replaced = biCamDSI.getLastBiScan(BiScan.BISCAN_SINGLECORR).copyFittestEnabled(
							afloor,                             // final double  str_floor,      // absolute strength floor
							clt_parameters.rig.pf_disp_afloor,  // final double  pf_disp_afloor, // =            0.1;    // When selecting the best fit from the alternative disparities, divide by difference increased by this
							clt_parameters.rig.pf_disp_rfloor); // 	final double  pf_disp_rfloor) //  =            0.02;   // Increase pf_disp_afloor for large disparities
					if ((debugLevel > -2) && clt_parameters.rig.rig_mode_debug){
						System.out.println("groundTruthByRigPlanes(): Replacing strongest by fittest: ntry = "+nfit+", replaced "+num_replaced+" tiles");
					}
					if (num_replaced == 0) {
						break;
					}
				}
				double  fg_str_good_enough =  clt_parameters.rig.pf_trusted_strength  * 0.6; // 4; // absolute strength floor for good enough
				double  fg_min_FGtoBG =      1.0;    // minimal disparity difference over
				double  fg_disp_atolerance = 0.1;    // Maximal absolute disparity difference to qualifying neighbor
				double  fg_disp_rtolerance = 0.02;   // Maximal relative (to absolute disparity) disparity difference to qualifying neighbor
				int     fg_min_neib =        2;      // minimal number of qualifying neighbors to promote FG tile

				// promote thin FG objects over even stronger BG ones (as thin stick in FG over textured BG)

				int num_replaced_fg = biCamDSI.getLastBiScan(BiScan.BISCAN_SINGLECORR).copyStrongFGEnabled(
						fg_str_good_enough, // final double  str_good_enough, // absolute strength floor for good enough
						fg_min_FGtoBG,      // final double  min_FGtoBG,      // minimal disparity difference over
						fg_disp_atolerance, // final double  disp_atolerance, // =  0.1;    // Maximal absolute disparity difference to qualifying neighbor
						fg_disp_rtolerance, // final double  disp_rtolerance, // =  0.02;   // Maximal relative (to absolute disparity) disparity difference to qualifying neighbor
						fg_min_neib);       // final int     min_neib)        // minimal number of qualifying neighbors to promote FG tile
				//				  if ((debugLevel > -2) && clt_parameters.rig.rig_mode_debug){
				if ((debugLevel > -4)){
					System.out.println("groundTruthByRigPlanes(): Replacing BG with FG tiles,  replaced "+num_replaced_fg+" tiles");
				}

				if (clt_parameters.show_map &&  (debugLevel > 0) && clt_parameters.rig.rig_mode_debug){
					biCamDSI.getLastBiScan(BiScan.BISCAN_SINGLECORR).showScan(
							quadCLT_main.image_name+"-BISCAN_"+num_fcycle+"-"+num_cycle);
				}
				if (last_cycle) { // last cycle
					if ((debugLevel > -2) && clt_parameters.rig.rig_mode_debug){
						System.out.println("groundTruthByRigPlanes(): that was refinement measurements after the trimming of hanging weak FG over BG");
					}
					break;
				}
				//	public void showScan(String title) {
				if (debugLevel > -2){
					System.out.println("groundTruthByRigPlanes(): num_cycle="+num_cycle);
				}
			} // for (int num_cycle = 0; num_cycle < num_cross_gaps_cycles; num_cycle++) {
			if (debugLevel > -2){
				System.out.println("groundTruthByRigPlanes(): num_fcycle="+num_fcycle);
			}
		}// 		  for (int num_fcycle = 0; num_fcycle < num_full_cycles; num_fcycle++) {



		// Fill in low-textured areas using averaged correlation
		double [][] rig_disparity_strength = null;
		biCamDSI_persistent = biCamDSI; // save for pole detection
		if (clt_parameters.rig.ltavg_en) {
			if (debugLevel > -2) {
				System.out.println("groundTruthByRigPlanes(): Processing low-textured areas with multi-tile correlation averaging");
			}
			// next method adds to the list of BiScans
			//			  double [][] ds_avg =
			//requires biCamDSI_persistent
			measureLowTextureAreas(
					quadCLT_main,   // QuadCLT            quadCLT_main,  // tiles should be set
					quadCLT_aux,    // QuadCLT            quadCLT_aux,  // tiles should be set
					clt_parameters, // EyesisCorrectionParameters.CLTParameters       clt_parameters,
					threadsMax,     // final int     threadsMax,  // maximal number of threads to launch
					updateStatus,   // final boolean updateStatus,
					debugLevel-2);    // final int    debugLevel) //
			// get last that was just added
			rig_disparity_strength = biCamDSI.getLastBiScan(BiScan.BISCAN_ANY).getDisparityStrength(
					false, // boolean only_strong,
					false, // boolean only_trusted,
					true); // boolean only_enabled,

		} else {
			// ignore any added low-texture areas
			rig_disparity_strength = biCamDSI.getLastBiScan(BiScan.BISCAN_SINGLECORR).getDisparityStrength(
					false, // boolean only_strong,
					false, // boolean only_trusted,
					true); // boolean only_enabled,
		}
		return rig_disparity_strength;
	}

	public boolean showBiScan(
			QuadCLT            quadCLT_main,  // tiles should be set
			QuadCLT            quadCLT_aux,  // tiles should be set
			CLTParameters       clt_parameters,
			final int                                      threadsMax,  // maximal number of threads to launch
			final boolean                                  updateStatus,
			final int                                      debugLevel) //  throws Exception
	{

		if ((quadCLT_main == null) ||
				(quadCLT_main.tp == null) ||
				(quadCLT_main.tp.clt_3d_passes == null) ||
				(biCamDSI_persistent== null) ||
				(biCamDSI_persistent.biScans== null)) {
			String msg = "Data is not available. Please run \"Ground truth\" first";
			IJ.showMessage("Error",msg);
			System.out.println(msg);
			return false;
		}
		int scan_index = biCamDSI_persistent.biScans.size()-1;
		boolean show_smooth =     false;
		boolean keep_unreliable = false;
		boolean keep_weak =       false;
		boolean keep_strong =     false;
		double  center_weight =   1.0;
		// from clt_parameters.rig
		double trusted_strength = clt_parameters.rig.pf_trusted_strength;
		double cond_rtrusted =    clt_parameters.rig.pf_cond_rtrusted;
		double strength_rfloor =  clt_parameters.rig.pf_strength_rfloor;
		double strength_pow =     clt_parameters.rig.pf_strength_pow;
		int smpl_radius =         clt_parameters.rig.pf_fourq_radius; // clt_parameters.rig.pf_smpl_radius;
		int smpl_num =            clt_parameters.rig.pf_smpl_num;
		int smpl_num_narrow =     clt_parameters.rig.pf_smpl_num_narrow;
		double smpl_fract =       clt_parameters.rig.pf_smpl_fract;
		double max_adiff =        clt_parameters.rig.pf_max_adiff;
		double max_rdiff =        clt_parameters.rig.pf_max_rdiff;
		double max_atilt =        clt_parameters.rig.pf_max_atilt;
		double max_rtilt =        clt_parameters.rig.pf_max_rtilt;
		double smpl_arms =        clt_parameters.rig.pf_smpl_arms;
		double smpl_rrms =        clt_parameters.rig.pf_smpl_rrms;
		double damp_tilt =        clt_parameters.rig.pf_damp_tilt;
		double rwsigma =          clt_parameters.rig.pf_rwsigma;
		double rwsigma_narrow =   clt_parameters.rig.pf_rwsigma_narrow;

		boolean use_alt =         clt_parameters.rig.pf_use_alt;
		double goal_fraction_rms= clt_parameters.rig.pf_goal_fraction_rms;    // Try to make rms to be this fraction of maximal acceptable by removing outliers
		double boost_low_density= clt_parameters.rig.pf_boost_low_density;    // Strength assigned to fake tiles from neighbors (the lower - the higher)

		int        fourq_min =    clt_parameters.rig.pf_fourq_min;
		int        fourq_gap =    clt_parameters.rig.pf_fourq_gap;

		boolean  run_avg =         true; // false; //ltavg_en
		int      lt_radius =       clt_parameters.rig.ltavg_radius;
		boolean  strong_only =     clt_parameters.rig.ltavg_dens_strong;
		int      need_tiles =      clt_parameters.rig.ltavg_dens_tiles;
		int      max_radius =      clt_parameters.rig.ltavg_dens_radius;

		double   min_disparity =   clt_parameters.rig.ltavg_min_disparity;

		double   max_density =     clt_parameters.rig.ltavg_max_density;
		int      gap_hwidth =      clt_parameters.rig.ltavg_gap_hwidth;
		int      clust_hwidth =    clt_parameters.rig.ltavg_clust_hwidth;
		int      extra_grow =      clt_parameters.rig.ltavg_extra_grow;
		// smoothing parameters
		boolean  smooth_strength = clt_parameters.rig.ltavg_smooth_strength;  // provide tile strength when smoothing target disparity
		double   neib_pull =       clt_parameters.rig.ltavg_neib_pull;    // pull to weighted average relative to pull to the original disparity value. If 0.0 - will only update former NaN-s
		int      max_iter =        clt_parameters.rig.ltavg_max_iter;      //
		double   min_change =      clt_parameters.rig.ltavg_min_change;   //


		int ref_smpl_radius =  clt_parameters.rig.ltavg_ref_smpl_radius; // 11;          // final int        smpl_radius,
		int ref_smpl_num =     clt_parameters.rig.ltavg_ref_smpl_num; // 40;          // final int        smpl_num_narrow,   //         = 3;      // Number after removing worst (should be >1)
		double ref_max_adiff = clt_parameters.rig.ltavg_ref_max_adiff; //0.15;        // final double     max_adiff,  // Maximal absolute difference between the center tile and friends
		double ref_max_rdiff = clt_parameters.rig.ltavg_ref_max_rdiff; //0.04;        // final double     max_rdiff, //  Maximal relative difference between the center tile and friends
		double ref_smpl_arms = clt_parameters.rig.ltavg_ref_smpl_arms; //0.1;         // final double     smpl_arms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
		double ref_smpl_rrms = clt_parameters.rig.ltavg_ref_smpl_rrms; //0.01;        // final double     smpl_rrms,        //      = 0.005;  // Maximal RMS/disparity in addition to smplRms
		int num_lt_refine =    clt_parameters.rig.ltavg_num_lt_refine; // 20; // make a parameter
		double strong_tol =    clt_parameters.rig.ltavg_strong_tol ; //
		double weak_tol =      clt_parameters.rig.ltavg_weak_tol ; //

		GenericJTabbedDialog gd = new GenericJTabbedDialog("Set CLT parameters",900,1100);
		gd.addTab("Genearl","Select bi-scan to show and process");

		gd.addNumericField("Scan index (0..."+(biCamDSI_persistent.biScans.size()-1),  scan_index, 0, 2, "",  "Display scan by index");

		gd.addCheckbox    ("Show smooth disparity/strength for the selected scan",  show_smooth, 		"Unchecked - just as is");
		gd.addCheckbox    ("Keep unreliable tiles",                                 keep_unreliable, 	"Unchecked - overwrite with smooth data");
		gd.addCheckbox    ("Keep weak (but trusted) tiles",                         keep_weak,       	"Unchecked - overwrite with smooth data");
		gd.addCheckbox    ("Keep strng trusted tiles",                              keep_strong, 	    "Unchecked - overwrite with smooth data");
		gd.addNumericField("Center weight - relative weight of the existing tile ", center_weight,  4,6,"",
				"0.0 - suggest new disparity over existing tiles without ant regard to the original value, 1.0 - same influence as any other tile");

		gd.addMessage     ("Parameters that are copied from the CLT parameters");
		gd.addNumericField("Strength sufficient without neighbors",                                               trusted_strength,  4,6,"",
				"Unconditionally trusted tile. Other stength values are referenceds as fraction of this value (after strength floor subtraction)");
		gd.addNumericField("Strength sufficient with neighbors support, fraction of the trusted strength",      cond_rtrusted,  4,6,"",
				"Strength that may be valid for the tile if there are neighbors in the same possibly tilted plane of the DSI (floor corrected)");

		gd.addNumericField("Fraction of trusted strength to subtract",                                            strength_rfloor,  4,6,"",
				"Strength floor to subtract from all strength values");
		gd.addNumericField("Raise strength-floor to this power",                                                  strength_pow,  4,6,"",
				"Currently just 1.0 - lenear");
		gd.addNumericField("How far to extend around known tiles (probably should increase this value up to?",    smpl_radius,  0,3,"tiles",
				"Process a aquare centered at the current tile withthe side of twice this value plus 1 (2*pf_smpl_radius + 1)");
		gd.addNumericField("Number after remaining in the sample square after removing worst fitting tiles",      smpl_num,  0,3,"",
				"When fitting planes the outliers are removed until the number of remaining tiles equals this value");
		gd.addNumericField("Number of remaining tiles when using narrow selection",                               smpl_num_narrow,  0,3,"",
				"Number of remaining tiles during initial palne fitting to the center pixels (it is later extended to include farther tiles)");
		gd.addNumericField("Fraction of the reamining tiles of all non-zero tiles?",                              smpl_fract,  4,6,"",
				"This value is combined to the previous one (absilute). Maximal of absolute and relative times number of all non-empty tiles is used");
		gd.addNumericField("Maximal absolute disparity difference between the plane and tiles that fit",          max_adiff,  4,6,"pix",
				"Maximal absolute disparity difference for fitting. Combined with the next one (relative) ");
		gd.addNumericField("Maximal relative (to center disparity) difference between the plane and tiles that fit",max_rdiff,  4,6,"pix/pix",
				"This value is multipled by the tile disparity and added to the maximal absolute difference");
		gd.addNumericField("Maximal absolute tile tilt in DSI space",                                             max_atilt,  4,6,"pix/tile",
				"Maximal disparity difference betweeing neighbor tiles for the tilted plane. Combined with the relative one (next), min of both limits applies");
		gd.addNumericField("Maximal relative (per pixel of disparity) tile tilt in DSI space",                    max_rtilt,  4,6,"1/tile",
				"Maximal relative (to center disparity) tilt. Near tiles (larger disparity may have larger differnce.");
		gd.addNumericField("Maximal absolute RMS of the remaining tiles in a sample",                             smpl_arms,  4,6,"pix",
				"After removing outliers RMS of the remaining tiles must be less than this value");
		gd.addNumericField("Maximal relative (to center disparity) RMS of the remaining tiles in a sample",       smpl_rrms,  4,6,"pix/pix",
				"Relative RMS times disparity is added to the absolute one");
		gd.addNumericField("Tilt cost for damping insufficient plane data",                                       damp_tilt,  4,6,"",
				"Regularisation to handle co-linear and even single-point planes, forcing fronto-parallel for single point, and minimal tilt for co-linear set");
		gd.addNumericField("Influence of far neighbors is reduced as a Gaussian with this sigma",                 rwsigma,  4,6,"",
				"Sigma is relative to selection radius (square half-side)");
		gd.addNumericField("Weight function Gaussian sigma (relative to radius) for initial plane fitting",       rwsigma_narrow,  4,6,"",
				"Weight function Gaussian sigma (relative to selection radius) for initial plane fitting. May be ~=1/radius");
		gd.addCheckbox    ("When fitting planes, look for alternative measure tiles",use_alt, 	        "Unchecked - only use the latest (current) tile");
		gd.addNumericField("Try to make rms to be this fraction of maximal acceptable by removing outliers",      goal_fraction_rms,  4,6,"pix",
				"When removing outliers to fit planes, stop removing when the RMS of the remaining drops below this fraction of the maxium allowed RMS (should be < 1.0");
		gd.addNumericField("Strength assigned to fake tiles from neighbors (the lower - the higher)",             boost_low_density,  4,6,"pix",
				"Returned strength assigned to the tiles increases with this value - seems to be a bug");

		gd.addNumericField("Each of the 4 corners should have at least this number of tiles",                     fourq_min,  0,3,"",
				"Apply (>0) only for filling gaps, not during expansion. It requires that every of the 4 corners of the sample square has this number of tiles for a plane");
		gd.addNumericField("Four corners center gap half-width (1 - 1 tile, 2 - 3 tiles, 3 - 5 tiles, ...",       fourq_gap,  0,3,"",
				"Specifies corners of the sample square that should have tiles remain, after removing centre columns and center rows");

		gd.addTab("LT Avg","Low texture correlatinaveraging");
		gd.addCheckbox    ("Measure with tile averaging",                                                         run_avg, 	        "");
		gd.addNumericField("Averaging radius (1 - 3x3 square, 2  - 5x5, ...",                                     lt_radius,  0,3,"",  "");
		gd.addCheckbox    ("Calculate density of strong trusted only (false include weak trusted)",               strong_only, 	    "");
		gd.addNumericField("Minimal number tiles to calculate density)",                                          need_tiles,  0,3,"",  "");
		gd.addNumericField("Maximal radius for measuruing density)",                                              max_radius,  0,3,"",  "");

		gd.addNumericField("Minimal disparity to apply filter",                                                   min_disparity,  4,6,"pix",
				"Farther objects will not be filtered");

		gd.addNumericField("Maximal density to consider it to be low textured area",                              max_density,  4,6,"",
				"Select areas with lower density");
		gd.addNumericField("Maximal radius of a void in low-texture selection to fill"  ,                         gap_hwidth,  0,3,"",
				"Low textured selection may have gaps that will be filled");
		gd.addNumericField("Minimal radius of a low-textured cluster to process",                                 clust_hwidth,  0,3,"",
				"Remove low-textured areas smaller that twice this size in each orthogonal directions");
		gd.addNumericField("Additionally grow low-textured areas selections",                                     extra_grow,  0,3,"",
				"Low textured areas will be grown by the radius of correlation averaging plus this value");




		gd.addCheckbox    ("Use tile strengths when filling gaps/smoothing",                                      smooth_strength,
				"Unchecked - consider all tiles to have the same strength");
		gd.addNumericField("Relative pull of the nieghbor tiles compared to the original disparity" ,             neib_pull,   4,6,"",
				"If set to 0.0 - only gaps will be filled, defined disparities will not be modified");
		gd.addNumericField("Maximal number of smoothing / gap filling iterations to perform",                     max_iter,    0,3,"",
				"Safety limit for smoothing iterations ");
		gd.addNumericField("Minimal disparity change to continue smoothing",                                      min_change,  4,6,"pix","");

		gd.addNumericField("How far to extend around a tile when refining averaging correlation measuremnts by planes ",    ref_smpl_radius,  0,3,"tiles",
				"Process a aquare centered at the current tile withthe side of twice this value plus 1 (2*pf_smpl_radius + 1)");
		gd.addNumericField("Number after remaining in the sample square after removing worst fitting tiles",      ref_smpl_num,  0,3,"",
				"When fitting planes the outliers are removed until the number of remaining tiles equals this value");
		gd.addNumericField("Maximal absolute disparity difference between the plane and tiles that fit",          ref_max_adiff,  4,6,"pix",
				"Maximal absolute disparity difference for fitting. Combined with the next one (relative) ");
		gd.addNumericField("Maximal relative (to center disparity) difference between the plane and tiles that fit",ref_max_rdiff,  4,6,"pix/pix",
				"This value is multipled by the tile disparity and added to the maximal absolute difference");
		gd.addNumericField("Maximal absolute RMS of the remaining tiles in a sample",                             ref_smpl_arms,  4,6,"pix",
				"After removing outliers RMS of the remaining tiles must be less than this value");
		gd.addNumericField("Maximal relative (to center disparity) RMS of the remaining tiles in a sample",       ref_smpl_rrms,  4,6,"pix/pix",
				"Relative RMS times disparity is added to the absolute one");


		gd.addNumericField("Maximal number of low texture/averaging correlation passes",                          num_lt_refine,  0,3,"",
				"Will also exit when maximal tile disparity change falls below ltavg_min_change (..continue smoothing above)");
		gd.addNumericField("Strong tile difference to averaged to be accepted",                                   strong_tol,  4,6,"pix",
				"When combining normal measurements with low texture/correlation averaging use strong normal if they are close to averaged");
		gd.addNumericField("Weak trusted tile difference to averaged to be accepted",                             weak_tol,  4,6,"pix",
				"When combining normal measurements with low texture/correlation averaging use weak normal if they are close to averaged");



		//		  boolean run_avg = false;
		//		  int     lt_radius = 1;


		gd.showDialog();
		if (gd.wasCanceled()) return false;
		scan_index =      (int) gd.getNextNumber();
		show_smooth =           gd.getNextBoolean();
		keep_unreliable =       gd.getNextBoolean();
		keep_weak =             gd.getNextBoolean();
		keep_strong  =          gd.getNextBoolean();
		center_weight =         gd.getNextNumber();
		// from clt_parameters.rig
		trusted_strength =      gd.getNextNumber();
		cond_rtrusted =         gd.getNextNumber();
		strength_rfloor =       gd.getNextNumber();
		strength_pow =          gd.getNextNumber();
		smpl_radius  =    (int) gd.getNextNumber();
		smpl_num =        (int) gd.getNextNumber();
		smpl_num_narrow = (int) gd.getNextNumber();
		smpl_fract =            gd.getNextNumber();
		max_adiff =             gd.getNextNumber();
		max_rdiff =             gd.getNextNumber();
		max_atilt =             gd.getNextNumber();
		max_rtilt =             gd.getNextNumber();
		smpl_arms =             gd.getNextNumber();
		smpl_rrms =             gd.getNextNumber();
		damp_tilt =             gd.getNextNumber();
		rwsigma =               gd.getNextNumber();
		rwsigma_narrow =        gd.getNextNumber();

		use_alt =               gd.getNextBoolean();
		goal_fraction_rms=      gd.getNextNumber();
		boost_low_density=      gd.getNextNumber();

		fourq_min=        (int) gd.getNextNumber();
		fourq_gap=        (int) gd.getNextNumber();

		//		  gd.addTab("LT Avg","Low texture correlatinaveraging");

		run_avg  =              gd.getNextBoolean();
		lt_radius=        (int) gd.getNextNumber();
		strong_only =           gd.getNextBoolean();
		need_tiles=       (int) gd.getNextNumber();
		max_radius=       (int) gd.getNextNumber();

		min_disparity =         gd.getNextNumber();

		max_density =           gd.getNextNumber();
		gap_hwidth=             (int) gd.getNextNumber();
		clust_hwidth=           (int) gd.getNextNumber();
		extra_grow=             (int) gd.getNextNumber();

		smooth_strength =       gd.getNextBoolean();
		neib_pull =             gd.getNextNumber();
		max_iter=         (int) gd.getNextNumber();
		min_change =            gd.getNextNumber();

		ref_smpl_radius=  (int) gd.getNextNumber();
		ref_smpl_num=     (int) gd.getNextNumber();
		ref_max_adiff=          gd.getNextNumber();
		ref_max_rdiff=          gd.getNextNumber();
		ref_smpl_arms=          gd.getNextNumber();
		ref_smpl_rrms=          gd.getNextNumber();

		num_lt_refine=    (int) gd.getNextNumber();
		strong_tol=             gd.getNextNumber();
		weak_tol=               gd.getNextNumber();

		BiScan biScan =  biCamDSI_persistent.biScans.get(scan_index);
		double [][] ds = null;
		if (show_smooth) {
			@SuppressWarnings("unused")
			double [][] final_ds=  measureLowTextureAreas(
					quadCLT_main,      // QuadCLT            quadCLT_main,  // tiles should be set
					quadCLT_aux,       // QuadCLT            quadCLT_aux,  // tiles should be set
					clt_parameters,    // EyesisCorrectionParameters.CLTParameters       clt_parameters,
					// from clt_parameters.rig
					trusted_strength,  // double      trusted_strength,
					cond_rtrusted,     // double      cond_rtrusted,
					strength_rfloor,   // double      strength_rfloor,
					strength_pow,      // double      strength_pow,
					smpl_radius,       // int         smpl_radius,
					smpl_num,          // int         smpl_num,
					smpl_num_narrow,   // int         smpl_num_narrow,
					smpl_fract,        // double      smpl_fract,
					max_adiff,         // double      max_adiff,
					max_rdiff,         // double      max_rdiff,
					max_atilt,         // double      max_atilt,
					max_rtilt,         // double      max_rtilt,
					smpl_arms,         // double      smpl_arms,
					smpl_rrms,         // double      smpl_rrms,
					damp_tilt,         // double      damp_tilt,
					rwsigma,           // double      rwsigma,
					rwsigma_narrow,    // double      rwsigma_narrow,
					use_alt,           // boolean     use_alt,
					goal_fraction_rms, // double      goal_fraction_rms,
					boost_low_density, // double      boost_low_density,
					fourq_min,         // int         fourq_min,
					fourq_gap,         // int         fourq_gap,
					lt_radius,         // int         lt_radius,
					strong_only,       // boolean     strong_only,
					need_tiles,        // int         need_tiles,
					max_radius,        // int         max_radius,
					min_disparity,     // double      min_disparity,
					max_density,       // double      max_density,
					gap_hwidth,        // int         gap_hwidth,
					clust_hwidth,      // int         clust_hwidth,
					extra_grow,        // int         extra_grow,
					// smoothing parameters
					smooth_strength,   // boolean     smooth_strength,
					neib_pull,         // double      neib_pull,
					max_iter,          // int         max_iter,
					min_change,        // double      min_change,
					ref_smpl_radius,   // int         ref_smpl_radius,
					ref_smpl_num,      // int         ref_smpl_num,
					ref_max_adiff,     // double      ref_max_adiff,
					ref_max_rdiff,     // double      ref_max_rdiff,
					ref_smpl_arms,     // double      ref_smpl_arms,
					ref_smpl_rrms,     // double      ref_smpl_rrms,
					num_lt_refine,     // int         num_lt_refine,
					strong_tol,        // double      strong_tol,
					weak_tol,          // double      weak_tol,
					clt_parameters.rig.ltavg_expand_lt,      // boolean    expand_lt, //  =           true;
					clt_parameters.rig.ltavg_expand_dist,    // int        expand_dist, //  =         4;
					clt_parameters.rig.ltavg_expand_tol,     // double     expand_tol, //  =          0.15;   // expand LT right and left if it ends with same or nearer tile
					clt_parameters.rig.ltavg_expand_floor,   // double      expand_floor, //  =      0.5;    // multiply single-tile strength floor for correlation-average
					clt_parameters.rig.ltavg_expand_sample_num,// int         expand_sample_num, //  =   5;      // minimal number of samples in expansion mode
					threadsMax,        // maximal number of threads to launch
					updateStatus,
					debugLevel);
		} else {
			biScan.showScan(quadCLT_main.image_name+"BiScan-"+scan_index,ds);
		}

		return true;
	}

	public double [][] measureLowTextureAreas(
			QuadCLT            quadCLT_main,  // tiles should be set
			QuadCLT            quadCLT_aux,  // tiles should be set
			CLTParameters       clt_parameters,
			final int     threadsMax,  // maximal number of threads to launch
			final boolean updateStatus,
			final int    debugLevel) //  throws Exception
	{
		return  measureLowTextureAreas(
				quadCLT_main,      // QuadCLT            quadCLT_main,  // tiles should be set
				quadCLT_aux,       // QuadCLT            quadCLT_aux,  // tiles should be set
				clt_parameters,    // EyesisCorrectionParameters.CLTParameters       clt_parameters,
				// from clt_parameters.rig
				clt_parameters.rig.pf_trusted_strength,  // double      trusted_strength,
				clt_parameters.rig.pf_cond_rtrusted,     // double      cond_rtrusted,
				clt_parameters.rig.pf_strength_rfloor,   // double      strength_rfloor,
				clt_parameters.rig.pf_strength_pow,      // double      strength_pow,
				clt_parameters.rig.pf_fourq_radius,      // clt_parameters.rig.pf_smpl_radius, // smpl_radius,       // int         smpl_radius,
				clt_parameters.rig.pf_smpl_num,          // int         smpl_num,
				clt_parameters.rig.pf_smpl_num_narrow,   // int         smpl_num_narrow,
				clt_parameters.rig.pf_smpl_fract,        // double      smpl_fract,
				clt_parameters.rig.pf_max_adiff,         // double      max_adiff,
				clt_parameters.rig.pf_max_rdiff,         // double      max_rdiff,
				clt_parameters.rig.pf_max_atilt,         // double      max_atilt,
				clt_parameters.rig.pf_max_rtilt,         // double      max_rtilt,
				clt_parameters.rig.pf_smpl_arms,         // double      smpl_arms,
				clt_parameters.rig.pf_smpl_rrms,         // double      smpl_rrms,
				clt_parameters.rig.pf_damp_tilt,         // double      damp_tilt,
				clt_parameters.rig.pf_rwsigma,           // double      rwsigma,
				clt_parameters.rig.pf_rwsigma_narrow,    // double      rwsigma_narrow,
				clt_parameters.rig.pf_use_alt,           // boolean     use_alt,
				clt_parameters.rig.pf_goal_fraction_rms, // double      goal_fraction_rms, // Try to make rms to be this fraction of maximal acceptable by removing outliers
				clt_parameters.rig.pf_boost_low_density, // double      boost_low_density,// Strength assigned to fake tiles from neighbors (the lower - the higher)
				clt_parameters.rig.pf_fourq_min,         // int         fourq_min,
				clt_parameters.rig.pf_fourq_gap,         // int         fourq_gap,
				clt_parameters.rig.ltavg_radius,         // int         lt_radius,
				clt_parameters.rig.ltavg_dens_strong,    // boolean     strong_only,
				clt_parameters.rig.ltavg_dens_tiles,     // int         need_tiles,
				clt_parameters.rig.ltavg_dens_radius,    // int         max_radius,
				clt_parameters.rig.ltavg_min_disparity,  // double      min_disparity,
				clt_parameters.rig.ltavg_max_density,    // double      max_density,
				clt_parameters.rig.ltavg_gap_hwidth,     // int         gap_hwidth,
				clt_parameters.rig.ltavg_clust_hwidth,   // int         clust_hwidth,
				clt_parameters.rig.ltavg_extra_grow,     // int         extra_grow,
				// smoothing parameters
				clt_parameters.rig.ltavg_smooth_strength,// boolean     smooth_strength,
				clt_parameters.rig.ltavg_neib_pull,      // double      neib_pull,
				clt_parameters.rig.ltavg_max_iter,       // int         max_iter,
				clt_parameters.rig.ltavg_min_change,     // double      min_change,
				clt_parameters.rig.ltavg_ref_smpl_radius,// int         ref_smpl_radius,
				clt_parameters.rig.ltavg_ref_smpl_num,   // int         ref_smpl_num,
				clt_parameters.rig.ltavg_ref_max_adiff,  // double      ref_max_adiff,
				clt_parameters.rig.ltavg_ref_max_rdiff,  // double      ref_max_rdiff,
				clt_parameters.rig.ltavg_ref_smpl_arms,  // double      ref_smpl_arms,
				clt_parameters.rig.ltavg_ref_smpl_rrms,  // double      ref_smpl_rrms,
				clt_parameters.rig.ltavg_num_lt_refine,  // int         num_lt_refine,
				clt_parameters.rig.ltavg_strong_tol,     // double      strong_tol,
				clt_parameters.rig.ltavg_weak_tol,       // double      weak_tol,
				clt_parameters.rig.ltavg_expand_lt,      // boolean    expand_lt, //  =           true;
				clt_parameters.rig.ltavg_expand_dist,    // int        expand_dist, //  =         4;
				clt_parameters.rig.ltavg_expand_tol,     // double     expand_tol, //  =          0.15;   // expand LT right and left if it ends with same or nearer tile
				clt_parameters.rig.ltavg_expand_floor,   // double      expand_floor, //  =      0.5;    // multiply single-tile strength floor for correlation-average
				clt_parameters.rig.ltavg_expand_sample_num,// int         expand_sample_num, //  =   5;      // minimal number of samples in expansion mode
				threadsMax,        // maximal number of threads to launch
				updateStatus,
				debugLevel);
	}


	public double [][] measureLowTextureAreas(
			QuadCLT            quadCLT_main,  // tiles should be set
			QuadCLT            quadCLT_aux,  // tiles should be set
			CLTParameters       clt_parameters,
			// from clt_parameters.rig
			double      trusted_strength,
			double      cond_rtrusted,
			double      strength_rfloor,
			double      strength_pow,
			int         smpl_radius,
			int         smpl_num,
			int         smpl_num_narrow,
			double      smpl_fract,
			double      max_adiff,
			double      max_rdiff,
			double      max_atilt,
			double      max_rtilt,
			double      smpl_arms,
			double      smpl_rrms,
			double      damp_tilt,
			double      rwsigma,
			double      rwsigma_narrow,
			boolean     use_alt,
			double      goal_fraction_rms,
			double      boost_low_density,
			int         fourq_min,
			int         fourq_gap,
			//			  boolean     run_avg,
			int         lt_radius,
			boolean     strong_only,
			int         need_tiles,
			int         max_radius,
			double      min_disparity,
			double      max_density,
			int         gap_hwidth,
			int         clust_hwidth,
			int         extra_grow,
			// smoothing parameters
			boolean     smooth_strength,
			double      neib_pull,
			int         max_iter,
			double      min_change,
			int         ref_smpl_radius,
			int         ref_smpl_num,
			double      ref_max_adiff,
			double      ref_max_rdiff,
			double      ref_smpl_arms,
			double      ref_smpl_rrms,
			int         num_lt_refine,
			double      strong_tol,
			double      weak_tol,
			boolean     expand_lt, //  =           true;
			int         expand_dist, //  =         4;
			double      expand_tol, //  =          0.15;   // expand LT right and left if it ends with same or nearer tile
			double      expand_floor, //  =      0.5;    // multiply single-tile strength floor for correlation-average
			int         expand_sample_num, //  =   5;      // minimal number of samples in expansion mode
			final int     threadsMax,  // maximal number of threads to launch
			final boolean updateStatus,
			final int    debugLevel) //  throws Exception
	{
		boolean show_smooth =     true;
		boolean keep_unreliable = false;
		boolean keep_weak =       false;
		boolean keep_strong =     false;
		double  center_weight =   1.0;

		if (debugLevel > -2){
			System.out.println(" === showBiScan( parameters : =====");
			//			  System.out.println("       scan_index= "+scan_index);
			System.out.println("      show_smooth= "+show_smooth);
			System.out.println("  keep_unreliable= "+keep_unreliable);
			System.out.println("        keep_weak= "+keep_weak);
			System.out.println("      keep_strong= "+keep_strong);
			System.out.println("    center_weight= "+center_weight);

			System.out.println(" trusted_strength= "+trusted_strength);
			System.out.println("    cond_rtrusted= "+cond_rtrusted);
			System.out.println("  strength_rfloor= "+strength_rfloor);
			System.out.println("     strength_pow= "+strength_pow);
			System.out.println("      smpl_radius= "+smpl_radius);
			System.out.println("         smpl_num= "+smpl_num);
			System.out.println("  smpl_num_narrow= "+smpl_num_narrow);
			System.out.println("       smpl_fract= "+smpl_fract);
			System.out.println("        max_adiff= "+max_adiff);
			System.out.println("        max_rdiff= "+max_rdiff);
			System.out.println("        max_atilt= "+max_atilt);
			System.out.println("        max_rtilt= "+max_rtilt);
			System.out.println("        smpl_arms= "+smpl_arms);
			System.out.println("        smpl_rrms= "+smpl_rrms);
			System.out.println("        damp_tilt= "+damp_tilt);
			System.out.println("          rwsigma= "+rwsigma);
			System.out.println("   rwsigma_narrow= "+rwsigma_narrow);

			System.out.println("          use_alt= "+use_alt);
			System.out.println("goal_fraction_rms= "+goal_fraction_rms);
			System.out.println("boost_low_density= "+boost_low_density);

			System.out.println("        fourq_min= "+fourq_min);
			System.out.println("        fourq_gap= "+fourq_gap);

			//			  System.out.println("          run_avg= "+run_avg);
			System.out.println("        lt_radius= "+lt_radius);
			System.out.println("      strong_only= "+strong_only);
			System.out.println("       need_tiles= "+need_tiles);
			System.out.println("       max_radius= "+max_radius);

			System.out.println("    min_disparity= "+min_disparity);
			System.out.println("      max_density= "+max_density);
			System.out.println("       gap_hwidth= "+gap_hwidth);
			System.out.println("       extra_grow= "+extra_grow);

			System.out.println("     clust_hwidth= "+clust_hwidth);
			System.out.println("  smooth_strength= "+smooth_strength);
			System.out.println("        neib_pull= "+neib_pull);
			System.out.println("         max_iter= "+max_iter);
			System.out.println("       min_change= "+min_change);

			System.out.println("  ref_smpl_radius= "+ref_smpl_radius);
			System.out.println("     ref_smpl_num= "+ref_smpl_num);
			System.out.println("    ref_max_adiff= "+ref_max_adiff);
			System.out.println("    ref_max_rdiff= "+ref_max_rdiff);
			System.out.println("    ref_smpl_arms= "+ref_smpl_arms);
			System.out.println("    ref_smpl_rrms= "+ref_smpl_rrms);

			System.out.println("    num_lt_refine= "+num_lt_refine);
			System.out.println("       strong_tol= "+strong_tol);
			System.out.println("         weak_tol= "+weak_tol);

			System.out.println("        expand_lt= "+expand_lt);
			System.out.println("      expand_dist= "+expand_dist);
			System.out.println("       expand_tol= "+expand_tol);

			System.out.println("     expand_floor= "+expand_floor);
			System.out.println("expand_sample_num= "+expand_sample_num);
		}
		BiScan biScan =  biCamDSI_persistent.getLastBiScan(BiScan.BISCAN_SINGLECORR); // biScans.get(scan_index);
		biCamDSI_persistent.getLastBiScan(BiScan.BISCAN_SINGLECORR).copyLastStrongestEnabled(
				clt_parameters.rig.pf_last_priority); // final boolean last_priority)
		double afloor = trusted_strength * strength_rfloor;
		int num_tries_strongest_by_fittest = 5;
		// replace strongest by fittest
		for (int nfit = 0; nfit < num_tries_strongest_by_fittest; nfit++) {
			int num_replaced = biCamDSI_persistent.getLastBiScan(BiScan.BISCAN_SINGLECORR).copyFittestEnabled(
					afloor,                             // final double  str_floor,      // absolute strength floor
					clt_parameters.rig.pf_disp_afloor,  // final double  pf_disp_afloor, // =            0.1;    // When selecting the best fit from the alternative disparities, divide by difference increased by this
					clt_parameters.rig.pf_disp_rfloor); // 	final double  pf_disp_rfloor) //  =            0.02;   // Increase pf_disp_afloor for large disparities
			if ((debugLevel > -2) && clt_parameters.rig.rig_mode_debug){
				System.out.println("measureLowTextureAreas(): Replacing strongest by fittest: ntry = "+nfit+", replaced "+num_replaced+" tiles");
			}
			if (num_replaced == 0) {
				break;
			}
		}

		double  fg_str_good_enough = trusted_strength * 0.4; // absolute strength floor for good enough
		double  fg_min_FGtoBG =      1.0;    // minimal disparity difference over
		double  fg_disp_atolerance = 0.1;    // Maximal absolute disparity difference to qualifying neighbor
		double  fg_disp_rtolerance = 0.02;   // Maximal relative (to absolute disparity) disparity difference to qualifying neighbor
		int     fg_min_neib =        2;      // minimal number of qualifying neighbors to promote FG tile

		// promote thin FG objects over even stronger BG ones (as thin stick in FG over textured BG)

		int num_replaced_fg = biCamDSI_persistent.getLastBiScan(BiScan.BISCAN_SINGLECORR).copyStrongFGEnabled(
				fg_str_good_enough, // final double  str_good_enough, // absolute strength floor for good enough
				fg_min_FGtoBG,      // final double  min_FGtoBG,      // minimal disparity difference over
				fg_disp_atolerance, // final double  disp_atolerance, // =  0.1;    // Maximal absolute disparity difference to qualifying neighbor
				fg_disp_rtolerance, // final double  disp_rtolerance, // =  0.02;   // Maximal relative (to absolute disparity) disparity difference to qualifying neighbor
				fg_min_neib);       // final int     min_neib)        // minimal number of qualifying neighbors to promote FG tile
		//		  if ((debugLevel > -2) && clt_parameters.rig.rig_mode_debug){
		if ((debugLevel > -2)){
			System.out.println("measureLowTextureAreas(): Replacing BG with FG tiles,  replaced "+num_replaced_fg+" tiles");
		}



		biScan.calcTrusted(   // finds strong trusted and validates week ones if they fit planes
				trusted_strength, // final double     trusted_strength, // trusted correlation strength
				strength_rfloor,  // final double     strength_rfloor,   // strength floor - relative to trusted
				cond_rtrusted,    // final double     cond_rtrusted,     // minimal strength to consider - fraction of trusted
				strength_pow,     // final double     strength_pow,      // raise strength-floor to this power
				smpl_radius,      // final int        smpl_radius,
				smpl_num,         // final int        smpl_num,   //         = 3;      // Number after removing worst (should be >1)
				smpl_fract,       // final double     smpl_fract, // Number of friends among all neighbors
				max_adiff,        // final double     max_adiff,  // Maximal absolute difference betweenthe center tile and friends
				max_rdiff,        // final double     max_rdiff, //  Maximal relative difference between the center tile and friends
				max_atilt,        // final double     max_atilt, //  = 2.0; // pix per tile
				max_rtilt,        // final double     max_rtilt, //  = 0.2; // (pix / disparity) per tile
				smpl_arms,        // final double     smpl_arms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
				smpl_rrms,        // final double     smpl_rrms,        //      = 0.005;  // Maximal RMS/disparity in addition to smplRms
				damp_tilt,        // final double     damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
				rwsigma,          // 						final double     rwsigma,           //  = 0.7; // influence of far neighbors diminish as a Gaussian with this sigma
				clt_parameters.tileX,                   // final int        dbg_x,
				clt_parameters.tileY,                   // final int        dbg_y,
				debugLevel);                            // final int        debugLevel);

		double [] density = 		  biScan.getDensity(
				strong_only, // final boolean strong_only,
				need_tiles, // 20, // 10,    // final int need_tiles,
				max_radius, // 20, // 15,    // final int max_radius,
				clt_parameters.tileX, // final int        dbg_x,
				clt_parameters.tileY, // final int        dbg_y,
				debugLevel+2);        // final int        debugLevel
		boolean [] pre_select = biScan.selectLowTextures(
				min_disparity,        // double    min_disparity,
				max_density,          // double    max_density,
				lt_radius+extra_grow, // int       grow,
				gap_hwidth,           // int       max_gap_radius,
				clust_hwidth,         // int       min_clust_radius,
				density,              // double [] density,
				null);                // double [] src_disparity);
		double []   dbg_presel = new double [pre_select.length];
		for (int i = 0; i < pre_select.length; i++) dbg_presel[i] = pre_select[i]? 1.0:0.0;
		double [][] dbg_dens_str = {density, dbg_presel};
		//		  biScan.showScan(quadCLT_main.image_name+"-density-"+scan_index,dbg_dens_str); //list_index
		if (debugLevel > 0) {
			biScan.showScan(quadCLT_main.image_name+"-density-"+biScan.list_index,dbg_dens_str); //list_index
		}

		double [][] ds = biScan.getFilteredDisparityStrength(
				pre_select,           // final boolean [] area_of_interest,
				null,                 // final double [][] disparityStrength,
				min_disparity,        // final double     min_disparity,    // keep original disparity far tiles
				trusted_strength,     // final double     trusted_strength, // trusted correlation strength
				strength_rfloor,      // final double     strength_rfloor,   // strength floor - relative to trusted
				!keep_unreliable,     // final boolean    discard_unreliable,// replace v
				!keep_weak,           // final boolean    discard_weak,      // consider weak trusted tiles (not promoted to trusted) as empty
				!keep_strong,         // final boolean    discard_strong,    // suggest new disparities even for strong tiles
				strength_pow,         // final double     strength_pow,      // raise strength-floor to this power
				null,                 // final double []  smpl_radius_array, // space-variant radius
				smpl_radius,          // final int        smpl_radius,
				smpl_num,             // final int        smpl_num,   //         = 3;      // Number after removing worst (should be >1)
				smpl_fract,           // final double     smpl_fract, // Number of friends among all neighbors
				smpl_num_narrow,      // final int        smpl_num_narrow,   //         = 3;      // Number after removing worst (should be >1)
				max_adiff,            // final double     max_adiff,  // Maximal absolute difference between the center tile and friends
				max_rdiff,            // final double     max_rdiff, //  Maximal relative difference between the center tile and friends
				max_atilt,            // final double     max_atilt, //  = 2.0; // pix per tile
				max_rtilt,            // final double     max_rtilt, //  = 0.2; // (pix / disparity) per tile
				smpl_arms,            // final double     smpl_arms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
				smpl_rrms,            // final double     smpl_rrms,        //      = 0.005;  // Maximal RMS/disparity in addition to smplRms
				damp_tilt,            // final double     damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
				rwsigma,              // final double     rwsigma,           //  = 0.7; // influence of far neighbors diminish as a Gaussian with this sigma
				rwsigma_narrow,       // final double     rwsigma_narrow,    //  = used to determine initial tilt
				center_weight,        // final double     center_weight,     // use center tile too (0.0 - do not use)
				use_alt,              // final boolean    use_alt,           // use tiles from other scans if they fit better
				goal_fraction_rms,    // final double     goal_fraction_rms, // Try to make rms to be this fraction of maximal acceptable by removing outliers
				boost_low_density,    //final double     boost_low_density, // 0 - strength is proportional to 1/density, 1.0 - same as remaining tiles
				fourq_min,            // final int        fourq_min,         // each of the 4 corners should have at least this number of tiles.
				fourq_gap,            // final int        fourq_gap,         // symmetrical vertical and horizontal center areas that do not belong to any corner
				clt_parameters.tileX, // final int        dbg_x,
				clt_parameters.tileY, // final int        dbg_y,
				debugLevel+0);          // final int        debugLevel
		if (debugLevel > 0) {
			biScan.showScan(quadCLT_main.image_name+"-BiScan-"+biScan.list_index,ds);
		}

		boolean [] lt_select = biScan.selectLowTextures(
				min_disparity,        // double    min_disparity,
				max_density,          // double    max_density,
				lt_radius+extra_grow, // 	int       grow,
				gap_hwidth,           // int       max_gap_radius,
				clust_hwidth,         // int       min_clust_radius,
				density,              // double [] density,
				ds[0]);               // double [] src_disparity);

		// compare iterations at lt_compare;
		boolean [] lt_compare = lt_select.clone();
		biCamDSI_persistent.tnImage.shrinkSelection(
				2*(lt_radius), // int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				lt_compare, // boolean [] tiles,
				null); // boolean [] prohibit)


		double [][] ds1 = {ds[0].clone(), ds[1].clone()} ;
		for (int i = 0; i < lt_select.length; i++) if (!lt_select[i]) {
			ds1[0][i] = Double.NaN;
			ds1[1][i] = 0.0;
		}
		if (debugLevel > 0) {
			biScan.showScan(quadCLT_main.image_name+"-selection-"+biScan.list_index,ds1);
		}
		double [] lt_strength = smooth_strength? ds1[1]:null;

		double [][] ds2 = biScan.fillAndSmooth(
				ds1[0],               // final double [] src_disparity,
				lt_strength,          // final double [] src_strength, // if not null will be used for weighted pull
				lt_select,            // final boolean [] selection,
				neib_pull,            // final double     neib_pull, // pull to weighted average relative to pull to the original disparity value. If 0.0 - will only update former NaN-s
				max_iter,             // final int max_iterations,
				min_change,           // final double min_change,
				clt_parameters.tileX, // final int        dbg_x,
				clt_parameters.tileY, // final int        dbg_y,
				debugLevel+0);        // final int        debugLevel
		if (debugLevel > 0) {
			biScan.showScan(quadCLT_main.image_name+"-smooth-"+biScan.list_index,ds2);
		}

		//		  if (run_avg) {
		int tilesX = quadCLT_main.tp.getTilesX();
		double [][] disparity_bimap = measureNewRigDisparity(
				quadCLT_main,   // QuadCLT                                  quadCLT_main,  // tiles should be set
				quadCLT_aux,    // QuadCLT                                  quadCLT_aux,
				ds2[0],          // double []                                      disparity, // Double.NaN - skip, ohers - measure
				clt_parameters, // EyesisCorrectionParameters.CLTParameters clt_parameters,
				false, // boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
				lt_radius,      // int                 lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using (2*notch_mode+1)^2 square
				// use set from parameters
				clt_parameters.rig.no_int_x0,    // boolean                                  no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
				threadsMax,     // final int        threadsMax,  // maximal number of threads to launch
				updateStatus,  // updateStatus,   // final boolean    updateStatus,
				debugLevel);    // final int        debugLevel)
		if (debugLevel > 0) {
			(new ShowDoubleFloatArrays()).showArrays(
					disparity_bimap,
					tilesX,
					disparity_bimap[0].length/tilesX,
					true,
					quadCLT_main.image_name+"LPF"+lt_radius,
					ImageDtt.BIDISPARITY_TITLES);
		}

		//try to refine
		int [] num_new = new int[1];
		/*
			  boolean [] trusted_measurements = 	  getTrustedDisparity(
					  quadCLT_main,                            // QuadCLT            quadCLT_main,  // tiles should be set
					  quadCLT_aux,                             // QuadCLT            quadCLT_aux,
					  clt_parameters.rig.min_trusted_strength, // double             min_combo_strength,    // check correlation strength combined for all 3 correlations
					  clt_parameters.grow_disp_trust,          // double             max_trusted_disparity, // 4.0 -> change to rig_trust
					  clt_parameters.rig.trusted_tolerance,    // double             trusted_tolerance,
					  null,                                    // boolean []         was_trusted,
					  disparity_bimap);              // double [][]        bimap // current state of measurements
		 */
		double [][] prev_bimap = null;
		double [] scale_bad = new double [ds2[0].length];
		for (int i = 0; i < scale_bad.length; i++) scale_bad[i] = 1.0;
		for (int nref = 0; nref < num_lt_refine; nref++) { // clt_parameters.rig.num_inf_refine; nref++) {
			//			  for (int nref = 0; nref < clt_parameters.rig.num_inf_refine; nref++) {

			double [][] disparity_bimap_new =  refineRigAvg(
					quadCLT_main,    // QuadCLT                        quadCLT_main,    // tiles should be set
					quadCLT_aux,     // QuadCLT                        quadCLT_aux,
					disparity_bimap, // double [][]                    src_bimap,       // current state of measurements (or null for new measurement)
					prev_bimap,      // double [][]                    prev_bimap, // previous state of measurements or null
					scale_bad,       // double []                      scale_bad,
					lt_select,       // final boolean []  area_of_interest,
					num_new,         // int     []                               num_new,
					clt_parameters,  // EyesisCorrectionParameters.CLTParameters clt_parameters,
					lt_radius,          // final int                                lt_radius,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using (2*notch_mode+1)^2 square
					biScan,          // final BiScan                             biScan,
					min_disparity, // final double     min_disparity,    // keep original disparity far tiles
					trusted_strength, // final double     trusted_strength, // trusted correlation strength
					expand_floor * strength_rfloor,// final double     strength_rfloor,   // strength floor - relative to trusted
					strength_pow,      //final double     strength_pow,      // raise strength-floor to this power
					ref_smpl_radius,   // final int        smpl_radius,
					smpl_fract,     // final double     smpl_fract, // Number of friends among all neighbors
					ref_smpl_num,   // final int        ref_smpl_num,   //         = 3;      // Number after removing worst (should be >1)
					ref_max_adiff,      // final double     max_adiff,  // Maximal absolute difference betweenthe center tile and friends
					ref_max_rdiff,      // final double     max_rdiff, //  Maximal relative difference between the center tile and friends
					max_atilt,      // final double     max_atilt, //  = 2.0; // pix per tile
					max_rtilt,      // final double     max_rtilt, //  = 0.2; // (pix / disparity) per tile
					ref_smpl_arms,  // final double     smpl_arms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
					ref_smpl_rrms,  // 							  final double     smpl_rrms,        //      = 0.005;  // Maximal RMS/disparity in addition to smplRms
					damp_tilt,      // final double     damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
					rwsigma,        // final double     rwsigma,           //  = 0.7; // influence of far neighbors diminish as a Gaussian with this sigma
					goal_fraction_rms, //final double     goal_fraction_rms, // Try to make rms to be this fraction of maximal acceptable by removing outliers
					max_iter,       // final int        max_iterations,
					min_change, // final double     min_change,
					clt_parameters.tileX, // final int        dbg_x,
					clt_parameters.tileY, // final int        dbg_y,
					threadsMax,      // final int                      threadsMax,      // maximal number of threads to launch
					updateStatus,    // final boolean                  updateStatus,
					debugLevel);     // final int                      debugLevel);

			prev_bimap = disparity_bimap;
			disparity_bimap = disparity_bimap_new;

			double max_diff =0.0, sw = 0.0, swd = 0.0, swd2 = 0.0;
			for (int nTile = 0; nTile < lt_compare.length; nTile++) if (lt_compare[nTile]){
				double w = disparity_bimap[ImageDtt.BI_STR_CROSS_INDEX][nTile]; // subtract floor?
				double d = disparity_bimap[ImageDtt.BI_TARGET_INDEX][nTile]-prev_bimap[ImageDtt.BI_TARGET_INDEX][nTile];
				if (Double.isNaN(d)) {
					System.out.println("showBiScan(): got NaN: disparity_bimap[ImageDtt.BI_TARGET_INDEX]["+nTile+"]="+disparity_bimap[ImageDtt.BI_TARGET_INDEX][nTile]+
							", prev_bimap[ImageDtt.BI_TARGET_INDEX]["+nTile+"]="+prev_bimap[ImageDtt.BI_TARGET_INDEX][nTile]);
				} else {
					sw += w;
					swd += w*d;
					swd2 +=w*d*d;
					max_diff = Math.max(max_diff, Math.abs(d));
				}
			}
			double mean = swd / sw;
			double rms = Math.sqrt(swd2 / sw);
			if (debugLevel > 0) {
				System.out.println("showBiScan() iteration "+nref+": mean ="+mean+", rms = "+ rms+", max diff. = "+max_diff );
			}

			if (debugLevel > 10) {(new ShowDoubleFloatArrays()).showArrays(
					disparity_bimap,
					tilesX,
					disparity_bimap[0].length/tilesX,
					true,
					quadCLT_main.image_name+"RE-MEASURED_R"+lt_radius+"-N"+nref,
					ImageDtt.BIDISPARITY_TITLES);
			}



			if (debugLevel > 0) {
				System.out.println("measureLowTextureAreas():  refinement step="+nref+" num_new= "+num_new[0]+" tiles");
			}
			if (num_new[0] < clt_parameters.rig.pf_min_new) break; // currently will never happen
			if (( max_diff < min_change) || (nref == (num_lt_refine - 1))) {
				if (debugLevel > -2) {
					System.out.println("showBiScan() final iteration "+nref+": mean ="+mean+", rms = "+ rms+", max diff. = "+max_diff );
				}
				break;
			}

		}
		if (debugLevel > 0) {
			(new ShowDoubleFloatArrays()).showArrays(
					disparity_bimap,
					tilesX,
					disparity_bimap[0].length/tilesX,
					true,
					quadCLT_main.image_name+"CORR-AVG"+lt_radius,
					ImageDtt.BIDISPARITY_TITLES);
		}

		double [][] avg_ds = {disparity_bimap[ImageDtt.BI_TARGET_INDEX],disparity_bimap[ImageDtt.BI_STR_CROSS_INDEX]};
		// maybe trim all previously added to the last BiScan.BISCAN_SINGLECORR?
		// so far just add
		for (int nTile = 0; nTile < lt_select.length; nTile++) if (!lt_select[nTile]){ // keep border tiles
			avg_ds[0][nTile] = Double.NaN;
			avg_ds[1][nTile] = 0.0;
		}

		boolean [] strong = biScan.strong_trusted;
		boolean [] weak =   biScan.trusted;
		double [][] ds_single = biScan.getDisparityStrength(
				false, // only_strong,
				false, // only_trusted,
				true); // only_enabled);
		if (expand_lt) {
			boolean [] expanded_lt = lt_select.clone();
			TileNeibs         tnImage =  biCamDSI_persistent.tnImage;
			tnImage.growSelection(
					2* expand_dist,   // int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
					expanded_lt,      // boolean [] tiles,
					null);            // boolean [] prohibit)
			for (int nTile = 0; nTile < expanded_lt.length; nTile++) {
				expanded_lt[nTile] &= ! lt_select[nTile] && !weak[nTile]; // Or just !Double.isNaN(ds_single[0][nTile]) ???
			}
			if (debugLevel > 0) {
				double [][] dbg_sel = new double [2][expanded_lt.length];
				for (int i = 0; i < expanded_lt.length; i++) {
					dbg_sel[0][i] = (lt_select[i]? 1:0) + (expanded_lt[i]? 2:0);
					dbg_sel[1][i] = (lt_select[i]? 1:0);
				}
				biScan.showScan(quadCLT_main.image_name+"-lt-selections"+biScan.list_index, dbg_sel);
				biScan.showScan(quadCLT_main.image_name+"-avg_ds"+biScan.list_index, avg_ds);

			}

			int expand_radius = 2*expand_dist; // where it looks for valid tiles
			double [][] ds_preexpanded =	biScan.getFilteredDisparityStrength(
					expanded_lt,          // final boolean [] area_of_interest,
					avg_ds,               // final double [][] disparityStrength,
					min_disparity,        // final double     min_disparity,    // keep original disparity far tiles
					trusted_strength,     // final double     trusted_strength, // trusted correlation strength
					expand_floor * strength_rfloor,// final double     strength_rfloor,   // strength floor - relative to trusted
					true,     // final boolean    discard_unreliable,// replace v
					true,           // final boolean    discard_weak,      // consider weak trusted tiles (not promoted to trusted) as empty
					true,         // final boolean    discard_strong,    // suggest new disparities even for strong tiles
					strength_pow,         // final double     strength_pow,      // raise strength-floor to this power
					null,                 // final double []  smpl_radius_array, // space-variant radius
					expand_radius,        // final int        smpl_radius,
					0,                    // final int        smpl_num,   //         = 3;      // Number after removing worst (should be >1)
					//					  0.5 * smpl_fract,        // final double     smpl_fract, // Number of friends among all neighbors
					smpl_fract,           // final double     smpl_fract, // Number of friends among all neighbors
					expand_sample_num,    // ref_smpl_num,         // final int        smpl_num_narrow,   //         = 3;      // Number after removing worst (should be >1)
					ref_max_adiff,        // final double     max_adiff,  // Maximal absolute difference between the center tile and friends
					ref_max_rdiff,        // final double     max_rdiff, //  Maximal relative difference between the center tile and friends
					max_atilt,            // final double     max_atilt, //  = 2.0; // pix per tile
					max_rtilt,            // final double     max_rtilt, //  = 0.2; // (pix / disparity) per tile
					ref_smpl_arms,        // final double     smpl_arms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
					ref_smpl_rrms,        // final double     smpl_rrms,        //      = 0.005;  // Maximal RMS/disparity in addition to smplRms
					damp_tilt,            // final double     damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
					0.0,                  // final double     rwsigma,           //  = 0.7; // influence of far neighbors diminish as a Gaussian with this sigma
					rwsigma,              // final double     rwsigma_narrow,    //  = used to determine initial tilt
					1.0,                  // final double     center_weight,     // use center tile too (0.0 - do not use)
					false,                // final boolean    use_alt,           // use tiles from other scans if they fit better
					goal_fraction_rms,    // final double     goal_fraction_rms, // Try to make rms to be this fraction of maximal acceptable by removing outliers
					0.8,                  // boost_low_density,    //final double     boost_low_density, // 0 - strength is proportional to 1/density, 1.0 - same as remaining tiles
					0,                    // final int        fourq_min,         // each of the 4 corners should have at least this number of tiles.
					0,                    // final int        fourq_gap,         // symmetrical vertical and horizontal center areas that do not belong to any corner
					clt_parameters.tileX, // final int        dbg_x,
					clt_parameters.tileY, // final int        dbg_y,
					debugLevel+0);          // final int        debugLevel
			if (debugLevel > -4) {
				biScan.showScan(quadCLT_main.image_name+"-preexpand"+biScan.list_index, ds_preexpanded);
			}

			for (int nTile = 0; nTile < expanded_lt.length; nTile++) {
				if (Double.isNaN(ds_preexpanded[0][nTile]) && !Double.isNaN(avg_ds[0][nTile])){
					ds_preexpanded[0][nTile] = avg_ds[0][nTile];
					ds_preexpanded[1][nTile] = avg_ds[1][nTile];
				}

			}
			if (debugLevel > -4) {// 0) {
				biScan.showScan(quadCLT_main.image_name+"-combo-expand"+biScan.list_index, ds_preexpanded); // already wrong strength
			}

			//			  double [][] ds_expanded =
			avg_ds = biScan.getLTExpanded(
					expand_tol,       // final double      tolerance, // should be not NaN over lt
					ds_preexpanded,   // final double [][] ds_lt, // should be not NaN over lt
					ds_single[0],     // final double []   d_single,
					lt_select,        // final boolean []  lt_sel,
					expanded_lt,      // final boolean []  exp_sel,
					weak);            // final boolean []  trusted);
			// TODO use fillAndSmooth() to calculate new weights (keeping disparity as it was)
			/*
			  double [][] avg_ds_strength = biScan.fillAndSmooth(
					  ds_preexpanded[0],    // final double [] src_disparity,
					  ds_preexpanded[1],    // final double [] src_strength, // if not null will be used for weighted pull
					  lt_select,            // final boolean [] selection,
					  0.0,                  // only gaps neib_pull, // final double     neib_pull, // pull to weighted average relative to pull to the original disparity value. If 0.0 - will only update former NaN-s
					  max_iter,             // final int max_iterations,
					  min_change,           // final double min_change,
					  clt_parameters.tileX, // final int        dbg_x,
					  clt_parameters.tileY, // final int        dbg_y,
					  debugLevel+0);        // final int        debugLevel

			  avg_ds[1] = avg_ds_strength[1];
			 */
			if (debugLevel > -4) { //-2) {
				biScan.showScan(quadCLT_main.image_name+"-lt-expanded"+biScan.list_index, avg_ds);
			}

		}

		int new_index = biCamDSI_persistent.addBiScan(
				avg_ds[0], // double [] disparity, // this will be "measured"
				avg_ds[1], // double [] strength,
				null, // boolean [] trusted,
				null, // boolean [] disabled,
				BiScan.BISCAN_AVGCORR); //int        scan_type)
		BiScan newScan = biCamDSI_persistent.getBiScan(new_index);
		if (debugLevel > -2) {
			System.out.println("Added scan #new_index");
		}

		for (int nTile = 0; nTile < lt_compare.length; nTile++){
			if (!lt_compare[nTile]) {
				if (!Double.isNaN(ds_single[0][nTile])) { // keep border from low texture if there are no normal measurements for this tile
					newScan.src_index[nTile] = biScan.src_index[nTile]; // will use same disparity/strength
				}
			} else {
				if (
						(strong[nTile] && (Math.abs(ds_single[0][nTile] - avg_ds[0][nTile]) <= strong_tol)) ||
						(weak[nTile] &&   (Math.abs(ds_single[0][nTile] - avg_ds[0][nTile]) <= weak_tol))) {
					newScan.src_index[nTile] = biScan.src_index[nTile]; // will use same disparity/strength
				}
			}
		}

		double [][] merged_ds = newScan.getDisparityStrength(
				false, // only_strong,
				false, // only_trusted,
				true); // only_enabled);



		// TODO: strength floor with averaged
		double avg_rfloor =        0.8 * strength_rfloor;
		double avg_cond_rtrusted = 0.8 * cond_rtrusted;

		int [] trusted_stats = newScan.calcTrusted(   // finds strong trusted and validates week ones if they fit planes
				trusted_strength,       // final double     trusted_strength, // trusted correlation strength
				avg_rfloor,             // final double     strength_rfloor,   // strength floor - relative to trusted
				avg_cond_rtrusted,      // final double     cond_rtrusted,     // minimal strength to consider - fraction of trusted
				strength_pow,           // final double     strength_pow,      // raise strength-floor to this power
				smpl_radius,            // final int        smpl_radius,
				smpl_num,               // final int        smpl_num,   //         = 3;      // Number after removing worst (should be >1)
				smpl_fract,             // final double     smpl_fract, // Number of friends among all neighbors
				max_adiff,              // final double     max_adiff,  // Maximal absolute difference betweenthe center tile and friends
				max_rdiff,              // final double     max_rdiff, //  Maximal relative difference between the center tile and friends
				max_atilt,              // final double     max_atilt, //  = 2.0; // pix per tile
				max_rtilt,              // final double     max_rtilt, //  = 0.2; // (pix / disparity) per tile
				smpl_arms,              // final double     smpl_arms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
				smpl_rrms,              // final double     smpl_rrms,        //      = 0.005;  // Maximal RMS/disparity in addition to smplRms
				damp_tilt,              // final double     damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
				rwsigma,                // 						final double     rwsigma,           //  = 0.7; // influence of far neighbors diminish as a Gaussian with this sigma
				clt_parameters.tileX,                   // final int        dbg_x,
				clt_parameters.tileY,                   // final int        dbg_y,
				debugLevel);                            // final int        debugLevel);


		if (debugLevel > -2) {
			System.out.println("measureLowTextureAreas()  strong trusted: "+trusted_stats[0]+
					" neib trusted: "+trusted_stats[1]+" weak trusted: " + trusted_stats[2]);
		}

		if (debugLevel > -2) {
			newScan.showScan(quadCLT_main.image_name+"-smooth-"+newScan.list_index, avg_ds);
		}
		return merged_ds;
	}




	public double [][] refineRigAvg(
			QuadCLT                                  quadCLT_main,  // tiles should be set
			QuadCLT                                  quadCLT_aux,
			double [][]                              src_bimap, // current state of measurements
			double [][]                              prev_bimap, // previous state of measurements or null
			double []                                scale_bad,
			final boolean []  area_of_interest,
			int     []                               num_new,
			CLTParameters clt_parameters,
			final int                                lt_radius,      // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using (2*notch_mode+1)^2 square
			final BiScan                             biScan,
			//			  final double [][] disparityStrength,
			final double     min_disparity,    // keep original disparity far tiles
			final double     trusted_strength, // trusted correlation strength
			final double     avg_strength_rfloor,   // strength floor - relative to trusted
			final double     strength_pow,      // raise strength-floor to this power
			final int        smpl_radius,
			final double     smpl_fract, // Number of friends among all neighbors
			final int        ref_smpl_num,   //         = 3;      // Number after removing worst (should be >1)
			final double     max_adiff,  // Maximal absolute difference betweenthe center tile and friends
			final double     max_rdiff, //  Maximal relative difference between the center tile and friends
			final double     max_atilt, //  = 2.0; // pix per tile
			final double     max_rtilt, //  = 0.2; // (pix / disparity) per tile
			final double     smpl_arms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
			final double     smpl_rrms,        //      = 0.005;  // Maximal RMS/disparity in addition to smplRms
			final double     damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
			final double     rwsigma,           //  = 0.7; // influence of far neighbors diminish as a Gaussian with this sigma
			final double     goal_fraction_rms, // Try to make rms to be this fraction of maximal acceptable by removing outliers
			final int        max_iterations,
			final double     min_change,
			final int        dbg_x,
			final int        dbg_y,
			final int        threadsMax,  // maximal number of threads to launch
			final boolean    updateStatus,
			final int        debugLevel)
	{
		int tilesX =quadCLT_main.tp.getTilesX();
		int tilesY =quadCLT_main.tp.getTilesY();
		int [][] tile_op = new int [tilesY][tilesX];
		double [][] disparity_array = new double [tilesY][tilesX];
		double disp_scale_main =  1.0/clt_parameters.corr_magic_scale; // Is it needed?
		double disp_scale_aux =   disp_scale_main * quadCLT_main.geometryCorrection.getDisparityRadius()/quadCLT_aux.geometryCorrection.getDisparityRadius();
		double disp_scale_inter = disp_scale_main * quadCLT_main.geometryCorrection.getDisparityRadius()/quadCLT_aux.geometryCorrection.getBaseline();
		int tile_op_all = clt_parameters.tile_task_op; //FIXME Use some constant?
		int numMeas = 0;
		double [][] ds_ref = new double[2][];
		ds_ref[0] = new double [tilesX*tilesY];

		for (int tileY = 0; tileY<tilesY;tileY++) {
			for (int tileX = 0; tileX<tilesX;tileX++) {
				int nTile = tileY * tilesX + tileX;
				if (((area_of_interest == null) || area_of_interest[nTile]) && !Double.isNaN(src_bimap[ImageDtt.BI_TARGET_INDEX][nTile])) {
					if (prepRefineTile(
							(lt_radius > 0),
							clt_parameters, // EyesisCorrectionParameters.CLTParameters       clt_parameters,
							tile_op_all,    // int                                            tile_op_all,
							src_bimap, // double [][]                                     src_bimap, // current state of measurements
							prev_bimap, // double [][]                                    prev_bimap, // previous state of measurements or null
							scale_bad,  // double []                                      scale_bad,
							tile_op, // int [][]                                          tile_op, // common for both amin and aux
							disparity_array, // double [][]                                    disparity_array,
							2, // int                                            refine_mode, // 0 - by main, 1 - by aux, 2 - by inter
							false,    // boolean                                        keep_inf,    // keep expected disparity 0.0 if it was so
							0.0,             // double                                         inf_disparity,
							0.0, // double                                         refine_min_strength, // do not refine weaker tiles
							0.0,    // double                                         refine_tolerance,    // do not refine if absolute disparity below
							disp_scale_main,  // double                                         disp_scale_main,  // 1.0
							disp_scale_aux,   //double                                         disp_scale_aux,   // ~0.58
							disp_scale_inter, //double                                         disp_scale_inter, // ~4.86
							//								  scale_step,       // double                                         scale_step,  // scale for "unstable tiles"
							tileX, // int                                            tileX,
							tileY, // int                                            tileY,
							nTile )) {
						numMeas++; //int                                            nTile
						ds_ref[0][nTile] = disparity_array[tileY][tileX];
					}
				}
			}
		}
		ds_ref[1] = src_bimap[ImageDtt.BI_STR_CROSS_INDEX];
		if (debugLevel > 0) {
			System.out.println("refineRigAvg(): prepared "+numMeas+" to measure");
		}

		double [][] ds_planes = biScan.getFilteredDisparityStrength(
				area_of_interest,     // final boolean [] area_of_interest,
				ds_ref,                 // final double [][] disparityStrength,
				min_disparity,        // final double     min_disparity,    // keep original disparity far tiles
				trusted_strength,     // final double     trusted_strength, // trusted correlation strength
				avg_strength_rfloor,// final double     strength_rfloor,   // strength floor - relative to trusted
				true,     // final boolean    discard_unreliable,// replace v
				true,           // final boolean    discard_weak,      // consider weak trusted tiles (not promoted to trusted) as empty
				true,         // final boolean    discard_strong,    // suggest new disparities even for strong tiles
				strength_pow,         // final double     strength_pow,      // raise strength-floor to this power
				null,                 // final double []  smpl_radius_array, // space-variant radius
				smpl_radius,          // final int        smpl_radius,
				0,                    // final int        smpl_num,   //         = 3;      // Number after removing worst (should be >1)
				smpl_fract,           // final double     smpl_fract, // Number of friends among all neighbors
				ref_smpl_num,         // final int        smpl_num_narrow,   //         = 3;      // Number after removing worst (should be >1)
				max_adiff,            // final double     max_adiff,  // Maximal absolute difference between the center tile and friends
				max_rdiff,            // final double     max_rdiff, //  Maximal relative difference between the center tile and friends
				max_atilt,            // final double     max_atilt, //  = 2.0; // pix per tile
				max_rtilt,            // final double     max_rtilt, //  = 0.2; // (pix / disparity) per tile
				smpl_arms,            // final double     smpl_arms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
				smpl_rrms,            // final double     smpl_rrms,        //      = 0.005;  // Maximal RMS/disparity in addition to smplRms
				damp_tilt,            // final double     damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
				0.0,                  // final double     rwsigma,           //  = 0.7; // influence of far neighbors diminish as a Gaussian with this sigma
				rwsigma,              // final double     rwsigma_narrow,    //  = used to determine initial tilt
				1.0,                  // final double     center_weight,     // use center tile too (0.0 - do not use)
				false,                // final boolean    use_alt,           // use tiles from other scans if they fit better
				goal_fraction_rms,    // final double     goal_fraction_rms, // Try to make rms to be this fraction of maximal acceptable by removing outliers
				0.8,                  // boost_low_density,    //final double     boost_low_density, // 0 - strength is proportional to 1/density, 1.0 - same as remaining tiles
				0,                    // final int        fourq_min,         // each of the 4 corners should have at least this number of tiles.
				0,                    // final int        fourq_gap,         // symmetrical vertical and horizontal center areas that do not belong to any corner
				clt_parameters.tileX, // final int        dbg_x,
				clt_parameters.tileY, // final int        dbg_y,
				debugLevel+0);          // final int        debugLevel

		// fill NaN gaps:
		double [][] ds_no_gaps = biScan.fillAndSmooth(
				ds_planes[0],         // disparity_bimap[ImageDtt.BI_TARGET_INDEX], // ds1[0], // final double [] src_disparity,
				null,                 // lt_strength, // final double [] src_strength, // if not null will be used for weighted pull
				area_of_interest,     // final boolean [] selection,
				0.0,                  // only gaps neib_pull, // final double     neib_pull, // pull to weighted average relative to pull to the original disparity value. If 0.0 - will only update former NaN-s
				max_iterations,       // final int max_iterations,
				min_change,           // final double min_change,
				clt_parameters.tileX, // final int        dbg_x,
				clt_parameters.tileY, // final int        dbg_y,
				debugLevel+0);        // final int        debugLevel

		// reformat back to disparity_array
		for (int tileY = 0; tileY<tilesY;tileY++) {
			for (int tileX = 0; tileX<tilesX;tileX++) {
				int nTile = tileY * tilesX + tileX;
				disparity_array[tileY][tileX] = ds_no_gaps[0][nTile];
			}
		}

		double [][] disparity_bimap  =  measureRig(
				quadCLT_main,        // QuadCLT                                        quadCLT_main,  // tiles should be set
				quadCLT_aux,         // QuadCLT                                        quadCLT_aux,
				tile_op,             // int [][]                                       tile_op, // common for both amin and aux
				disparity_array,     // double [][]                                    disparity_array,
				null, // double [][]                                    ml_data,         // data for ML - 10 layers - 4 center areas (3x3, 5x5,..) per camera-per direction, 1 - composite, and 1 with just 1 data (target disparity)
				clt_parameters,      // EyesisCorrectionParameters.CLTParameters       clt_parameters,
				clt_parameters.getFatZero(quadCLT_main.isMonochrome()), // double                                         fatzero,
				false,                          //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
				lt_radius,                       // final int  // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
				true,                            // final boolean             no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
				threadsMax,          //final int        threadsMax,  // maximal number of threads to launch
				updateStatus,        // final boolean    updateStatus,
				debugLevel);          // final int        debugLevel)

		// combine with old results for tiles that were not re-measured
		// not needed here as so far everything selected is re-measured in the averaging mode)
		/*
		  for (int tileY = 0; tileY<tilesY;tileY++) {
			  for (int tileX = 0; tileX<tilesX;tileX++) {
				  int nTile = tileY * tilesX + tileX;
				  if ((selection == null) || selection[nTile]) {
					  if (Double.isNaN(disparity_bimap[ImageDtt.BI_TARGET_INDEX][nTile]) && !Double.isNaN(src_bimap[ImageDtt.BI_TARGET_INDEX][nTile])) {
						  for (int i = 0; i < disparity_bimap.length; i++) {
							  disparity_bimap[i][nTile] = src_bimap[i][nTile];
						  }
					  }
				  }
			  }
		  }
		 */
		if (num_new != null) {
			num_new[0] = numMeas;
		}
		return disparity_bimap;
	}




	public double [][] fillPoorTextureByInter(
			QuadCLT            quadCLT_main,  // tiles should be set
			QuadCLT            quadCLT_aux,
			CLTParameters       clt_parameters,
			double [][]                                    disparity_bimap,
			BiCamDSI                                       biCamDSI,
			final int                                      threadsMax,  // maximal number of threads to launch
			final boolean                                  updateStatus,
			final int                                      debugLevel)// throws Exception
	{
		final int refine_inter = 2; // use inter-cam disparity for refinement
		final int tilesX = quadCLT_main.tp.getTilesX();
		//		  final int tilesY = quadCLT_main.tp.getTilesY();
		double [] suggestedLTMeasurements =  biCamDSI.suggestLTTiles(
				disparity_bimap,                           // double [][] disparity_bimap,
				null, // boolean []  trusted,       // may be null if disparity is alreasdy NaN-ed
				clt_parameters.rig.lt_min_disparity,       // double      min_disparity, //  =         0.0;    // apply low texture to near objects
				clt_parameters.rig.lt_trusted_strength,    // double      trusted_strength, // =       0.2;    // strength sufficient without neighbors
				clt_parameters.rig.lt_strength_rfloor,     // double      strength_rfloor,  // =       0.28;   // strength floor relative to trusted_strength
				clt_parameters.rig.lt_need_friends,        // double      need_friends, // =           0.4;    // strength sufficient with neighbors support, fraction of lt_trusted_strength
				clt_parameters.rig.lt_extend_dist,         // int         extend_dist, // =            3;      // how far to extend around known tiles (probably should increase this value up to?
				// dealing with neighbors variance
				clt_parameters.rig.lt_wsigma,              // double      wsigma,     //  = 1.0; // influence of far neighbors diminish as a Gaussian with this sigma
				clt_parameters.rig.lt_max_asigma,          // double      max_asigma, // =             .15;     // Maximal acceptable standard deviation of the neighbors (remove, then add)
				clt_parameters.rig.lt_max_rsigma,          // double      max_rsigma, // =             .05;     // Maximal acceptable standard deviation of the neighbors (remove, then add)
				debugLevel);                               // int         debugLevel


		double [][] disparity_bimap_lt = setBimapFromDisparityNaN(
				suggestedLTMeasurements,      // double []                                      disparity,
				quadCLT_main,         // QuadCLT                                  quadCLT_main,  // tiles should be set
				quadCLT_aux,          // QuadCLT                                  quadCLT_aux,
				clt_parameters,       // EyesisCorrectionParameters.CLTParameters clt_parameters,
				false,                //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
				0,                 // final int  // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
				// first measurement - use default value:
				clt_parameters.rig.no_int_x0,           // boolean                                  no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
				threadsMax,           // final int        threadsMax,  // maximal number of threads to launch
				updateStatus,         // final boolean    updateStatus,
				debugLevel);          // final int        debugLevel);

		if (clt_parameters.show_map &&  (debugLevel > -2) && clt_parameters.rig.rig_mode_debug){
			(new ShowDoubleFloatArrays()).showArrays(
					disparity_bimap_lt,
					tilesX,
					disparity_bimap[0].length/tilesX,
					true,
					quadCLT_main.image_name+"NEW_LT_MEASURED",
					ImageDtt.BIDISPARITY_TITLES);
		}

		// refine just the new suggested measurements
		double [][] prev_bimap = null;
		int [] num_new =        new int[1];
		double [] scale_bad = new double [suggestedLTMeasurements.length];
		for (int i = 0; i < scale_bad.length; i++) scale_bad[i] = 1.0;
		prev_bimap = null;
		boolean [] trusted_lt = null;
		for (int nref = 0; nref < clt_parameters.rig.num_near_refine; nref++) {
			// refine infinity using inter correlation
			double [][] disparity_bimap_new =  refineRigSel(
					quadCLT_main,    // QuadCLT                        quadCLT_main,    // tiles should be set
					quadCLT_aux,     // QuadCLT                        quadCLT_aux,
					disparity_bimap_lt, // double [][]                    src_bimap,       // current state of measurements (or null for new measurement)
					prev_bimap,      // double [][]                    prev_bimap, // previous state of measurements or null
					scale_bad,       // double []                      scale_bad,
					refine_inter,    // int                            refine_mode,     // 0 - by main, 1 - by aux, 2 - by inter
					false,           // boolean                        keep_inf,        // keep expected disparity 0.0 if it was so
					0.0,             // double                                         inf_disparity,
					0.0, // clt_parameters.rig.refine_min_strength , // double refine_min_strength, // do not refine weaker tiles
					clt_parameters.rig.refine_tolerance ,    // double refine_tolerance,    // do not refine if absolute disparity below
					trusted_lt,      // null, // trusted_lt,    // tile_list,       // ArrayList<Integer>             tile_list,       // or null
					num_new,         // int     []                                     num_new,
					clt_parameters,  // EyesisCorrectionParameters.CLTParameters clt_parameters,
					false,           //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
					0,                 // final int  // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
					// in refine mode disable int preset of the window
					true,             // boolean                                  no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide

					threadsMax,      // final int                      threadsMax,      // maximal number of threads to launch
					updateStatus,    // final boolean                  updateStatus,
					debugLevel);     // final int                      debugLevel);
			prev_bimap = disparity_bimap_lt;
			disparity_bimap_lt = disparity_bimap_new;
			//low texture may have very poor individual correlations, do not check them at all
			trusted_lt = 	getTrustedDisparityInter(
					0.0, // clt_parameters.rig.lt_trusted_strength*clt_parameters.rig.lt_need_friends, // double             min_inter_strength,    // check correlation strength combined for all 3 correlations
					clt_parameters.grow_disp_trust,          // double             max_trusted_disparity,
					trusted_lt, // boolean []         was_trusted,
					disparity_bimap_lt );                       // double [][]        bimap // current state of measurements

			if (debugLevel > -2) {
				System.out.println("enhanceByRig(): refined (lt) "+num_new[0]+" tiles");
			}
			if (num_new[0] < clt_parameters.rig.min_new) break;
		}
		trusted_lt = 	getTrustedDisparityInter(
				clt_parameters.rig.lt_trusted_strength*clt_parameters.rig.lt_need_friends, // double             min_inter_strength,    // check correlation strength combined for all 3 correlations
				clt_parameters.grow_disp_trust,          // double             max_trusted_disparity,
				trusted_lt, // boolean []         was_trusted,
				disparity_bimap_lt );                       // double [][]        bimap // current state of measurements


		if (clt_parameters.show_map &&  (debugLevel > -2) && clt_parameters.rig.rig_mode_debug){
			(new ShowDoubleFloatArrays()).showArrays(
					disparity_bimap_lt,
					tilesX,
					disparity_bimap[0].length/tilesX,
					true,
					quadCLT_main.image_name+"NEW_LT_REFINED",
					ImageDtt.BIDISPARITY_TITLES);
		}

		// combine new measured results with the previously known (new overwrites old
		for (int nTile = 0; nTile < disparity_bimap_lt[0].length; nTile++) {
			//    		  if (trusted_lt[nTile] &&
			//   				  (!trusted_near[nTile] ||
			//    						  (disparity_bimap_infinity[ImageDtt.BI_STR_ALL_INDEX][nTile] > disparity_bimap[ImageDtt.BI_STR_ALL_INDEX][nTile]))) {
			// start with unconditional use of new
			if (trusted_lt[nTile]) {
				for (int i = 0; i < disparity_bimap.length; i++) if (disparity_bimap!=null) {
					disparity_bimap[i][nTile] = disparity_bimap_lt[i][nTile];
				}
				//    			  trusted_near[nTile] = true;
			}
		}
		if (clt_parameters.show_map &&  (debugLevel > -2) && clt_parameters.rig.rig_mode_debug){
			(new ShowDoubleFloatArrays()).showArrays(
					disparity_bimap_lt,
					tilesX,
					disparity_bimap[0].length/tilesX,
					true,
					quadCLT_main.image_name+"DSI_ADDED",
					ImageDtt.BIDISPARITY_TITLES);
		}
		return disparity_bimap;
	}




	public void saveMlFile(
			String               ml_title,
			String               ml_directory,
			// new meanings in aux mode:
			// disp_offset_low NaN:
			//   disp_offset_high = 0.0 use average GT disparity
			//   disp_offset_high > 0.0 use FG GT disparity
			//   disp_offset_high < 0.0 use BG GT disparity
			// disp_offset_low != NaN: it is absolute disparity

			double               disp_offset_low,  // NaN - Main camera is used
			double               disp_offset_high, // !NaN and  isNaN(disp_offset_low) - random amplitude: positive - from main, negative - from rig
			QuadCLT              quadCLT_main, // use null for aux mode
			QuadCLT              quadCLT_aux,
			Correlation2d        corr2d, // to access "other" layer
			boolean              use8bpp,
			double               limit_extrim,
			boolean              keep_aux,
			boolean              keep_inter,
			boolean              keep_hor_vert,
			boolean              ml_keep_tbrl,
			boolean              keep_debug,
			double               ml_fatzero,
			int                  ml_hwidth,
			double [][]          ml_data,
			boolean              show,
			int                  debugLevel
			) {
		final boolean has_main = (quadCLT_main != null);
		final boolean aux_mode = !has_main;
		final int tilesX = (has_main) ? quadCLT_main.tp.getTilesX() : quadCLT_aux.tp.getTilesX();
		final int tilesY = (has_main) ? quadCLT_main.tp.getTilesY() : quadCLT_aux.tp.getTilesY();
		int ml_width = 2 * ml_hwidth + 1;
		int width =  tilesX * ml_width;
		int height = tilesY * ml_width;
		String title = ml_title+ (use8bpp?"08":"32")+"B-"+(keep_aux?"A":"")+(keep_inter?"I":"")+(keep_hor_vert?"O":"")+(ml_keep_tbrl?"T":"")+
				(keep_debug?"D":"")+"-FZ"+ml_fatzero;
		if (aux_mode) {
			if (!Double.isNaN(disp_offset_low)) {
				title += "-D" + String.format("%08.5f",disp_offset_low).trim(); // absolute disparity
			} else if (disp_offset_high > 0.0){ // foreground ground truth disparity
				title += "-FG";
			} else if (disp_offset_high < 0.0){ // background ground truth disparity
				title += "-BG";
			} else { // average ground truth disparity
				title += "-AG";
			}
		} else {
			if (!Double.isNaN(disp_offset_low)) {
				title += "-OFFS" + String.format("%08.5f",disp_offset_low).trim();
				if (disp_offset_high > disp_offset_low) {
					title+="_";
					title+=String.format("%08.5f",disp_offset_high).trim();
				}
			} else {
				if (Double.isNaN(disp_offset_high)) {
					title += "-MAIN";
				} else {
					if (disp_offset_high > 0) {
						title += "-MAIN_RND";
						title+=String.format("%08.5f",disp_offset_high).trim();
					} else {
						disp_offset_high = -disp_offset_high;
						title+="-RIG_RND";
						title+=String.format("%08.5f",disp_offset_high).trim();
					}
				}
			}
		}
		int [] main_indices = {
				ImageDtt.ML_TOP_INDEX,    // 8 - top pair 2d correlation center area (auxiliary camera)
				ImageDtt.ML_BOTTOM_INDEX, // 9 - bottom pair 2d correlation center area (auxiliary camera)
				ImageDtt.ML_LEFT_INDEX,   //10 - left pair 2d correlation center area (auxiliary camera)
				ImageDtt.ML_RIGHT_INDEX,  //11 - right pair 2d correlation center area (auxiliary camera)
				ImageDtt.ML_DIAGM_INDEX,  //12 - main diagonal (top-left to bottom-right) pair 2d correlation center area (auxiliary camera)
				ImageDtt.ML_DIAGO_INDEX,  //13 - other diagonal (bottom-left to top-right) pair 2d correlation center area (auxiliary camera)
				ImageDtt.ML_HOR_INDEX,    //14 - horizontal pairs combined 2d correlation center area (auxiliary camera)
				ImageDtt.ML_VERT_INDEX    //15 - vertical pairs combined 2d correlation center area (auxiliary camera)
		};
		int [] aux_indices = {
				ImageDtt.ML_TOP_AUX_INDEX,    // 8 - top pair 2d correlation center area (auxiliary camera)
				ImageDtt.ML_BOTTOM_AUX_INDEX, // 9 - bottom pair 2d correlation center area (auxiliary camera)
				ImageDtt.ML_LEFT_AUX_INDEX,   //10 - left pair 2d correlation center area (auxiliary camera)
				ImageDtt.ML_RIGHT_AUX_INDEX,  //11 - right pair 2d correlation center area (auxiliary camera)
				ImageDtt.ML_DIAGM_AUX_INDEX,  //12 - main diagonal (top-left to bottom-right) pair 2d correlation center area (auxiliary camera)
				ImageDtt.ML_DIAGO_AUX_INDEX,  //13 - other diagonal (bottom-left to top-right) pair 2d correlation center area (auxiliary camera)
				ImageDtt.ML_HOR_AUX_INDEX,    //14 - horizontal pairs combined 2d correlation center area (auxiliary camera)
				ImageDtt.ML_VERT_AUX_INDEX    //15 - vertical pairs combined 2d correlation center area (auxiliary camera)
		};
		int [] inter_indices = {
				ImageDtt.ML_INTER_INDEX       //16 - inter-camera (between two quad ones) correlation center area
		};
		int [] hor_vert_indices = {
				ImageDtt.ML_HOR_INDEX,        // 6 - horizontal pairs combined 2d correlation center area
				ImageDtt.ML_VERT_INDEX,       // 7 - vertical pairs combined 2d correlation center area
				ImageDtt.ML_HOR_AUX_INDEX,    //14 - horizontal pairs combined 2d correlation center area (auxiliary camera)
				ImageDtt.ML_VERT_AUX_INDEX    //15 - vertical pairs combined 2d correlation center area (auxiliary camera)
		};

		int [] tbrl_indices = {
				ImageDtt.ML_TOP_INDEX,        // 0 - top pair 2d correlation center area
				ImageDtt.ML_BOTTOM_INDEX,     // 1 - bottom pair 2d correlation center area
				ImageDtt.ML_LEFT_INDEX,       // 2 - left pair 2d correlation center area
				ImageDtt.ML_RIGHT_INDEX,      // 3 - right pair 2d correlation center area
				ImageDtt.ML_TOP_AUX_INDEX,    // 8 - top pair 2d correlation center area (auxiliary camera)
				ImageDtt.ML_BOTTOM_AUX_INDEX, // 9 - bottom pair 2d correlation center area (auxiliary camera)
				ImageDtt.ML_LEFT_AUX_INDEX,   //10 - left pair 2d correlation center area (auxiliary camera)
				ImageDtt.ML_RIGHT_AUX_INDEX   //11 - right pair 2d correlation center area (auxiliary camera)
		};
		int [] dbg_indices = {
				ImageDtt.ML_DBG1_INDEX        //18 - just debug data (first - auto phase correlation)
		};
		int [] non_corr_indices = {
				ImageDtt.ML_OTHER_INDEX,      //17 - other data: 0 (top left tile corner) - preset disparity of the tile, 1: (next element) - ground trouth data, 2:
				ImageDtt.ML_DBG1_INDEX        //18 - just debug data (first - auto phase correlation)
		};

		boolean [] skip_layers = new boolean [ImageDtt.ML_TITLES.length];
		if (!has_main)       for (int nl:main_indices)     skip_layers[nl] = true;
		if (!keep_aux)       for (int nl:aux_indices)      skip_layers[nl] = true;
		if (!keep_inter)     for (int nl:inter_indices)    skip_layers[nl] = true;
		if (!keep_hor_vert)  for (int nl:hor_vert_indices) skip_layers[nl] = true;
		if (!ml_keep_tbrl)   for (int nl:tbrl_indices)     skip_layers[nl] = true;

		if (!keep_debug)     for (int nl:dbg_indices)      skip_layers[nl] = true;

		ImageStack array_stack=new ImageStack(width,height);
		double soft_mn = Double.NaN,soft_mx = Double.NaN;
		if (use8bpp) {
			int num_bins = 256;
			boolean [] skip_histogram =   skip_layers.clone();
			for (int nl:non_corr_indices) skip_histogram[nl] = true;
			double mn = 0.0, mx = 0.0; // data has both positive and negative values
			for (int nl = 0; nl < ml_data.length; nl++) if (!skip_histogram[nl]) {
				for (int i = 0; i < ml_data[nl].length; i++) if (!Double.isNaN(ml_data[nl][i])){
					if      (ml_data[nl][i] > mx) mx = ml_data[nl][i];
					else if (ml_data[nl][i] < mn) mn = ml_data[nl][i];
				}
			}
			if (debugLevel > -2) {
				System.out.println("saveMlFile(): min="+mn+", max="+mx);
			}
			int [] histogram = new int [num_bins];
			int num_values = 0;
			for (int nl = 0; nl < ml_data.length; nl++) if (!skip_histogram[nl]) {
				for (int i = 0; i < ml_data[nl].length; i++) if (!Double.isNaN(ml_data[nl][i])){
					int bin = (int) Math.round(num_bins*(ml_data[nl][i] - mn)/(mx-mn));
					// rounding errors?
					if (bin < 0) bin = 0;
					else if (bin >= num_bins) bin = num_bins-1;
					histogram[bin]++;
					num_values++;
				}
			}
			double ignore_vals = limit_extrim*num_values;
			soft_mn = mn;
			soft_mx = mx;
			{
				double sl = 0.0;
				int i = 0;
				while (sl < ignore_vals) {
					i++;
					sl+= histogram[i];
					soft_mn += (mx-mn)/num_bins;
				}
				double f = (sl - ignore_vals)/histogram[i];
				soft_mn -= (mx-mn)/num_bins*(1.0 - f);

				sl = 0.0;
				i = num_bins-1;
				while (sl < ignore_vals) {
					i--;
					sl+= histogram[i];
					soft_mx -= (mx-mn)/num_bins;
				}
				f = (sl - ignore_vals)/histogram[i];
				soft_mn += (mx-mn)/num_bins*(1.0 - f);
				if (debugLevel > -2) {
					System.out.println("saveMlFile(): soft min="+soft_mn+", soft max="+soft_mx);
				}
			}
			// convert double data  to byte, so v<=soft_mn -> 1; v>= soft_mx -> 255, NaN -> 0
			int [] signed_data = {
					ImageDtt.ML_OTHER_TARGET,
					ImageDtt.ML_OTHER_GTRUTH,
					ImageDtt.ML_OTHER_GTRUTH_FG_DISP,
					ImageDtt.ML_OTHER_GTRUTH_BG_DISP,
					ImageDtt.ML_OTHER_AUX_DISP};
			int [] unsigned_data = {
					ImageDtt.ML_OTHER_GTRUTH_STRENGTH,
					ImageDtt.ML_OTHER_GTRUTH_RMS,
					ImageDtt.ML_OTHER_GTRUTH_RMS_SPLIT,
					ImageDtt.ML_OTHER_GTRUTH_FG_STR,
					ImageDtt.ML_OTHER_GTRUTH_BG_STR,
					ImageDtt.ML_OTHER_AUX_STR};
			byte [][] iml_data = new byte [ml_data.length][];
			for (int nl = 0; nl < ml_data.length; nl++) if (!skip_layers[nl]) {
				iml_data[nl] = new byte [ml_data[nl].length];
				if (nl == ImageDtt.ML_OTHER_INDEX) {
					// special treatment - make 2 bytes of one disparity value
					for (int tileY = 0; tileY < tilesY; tileY++) {
						for (int tileX = 0; tileX < tilesX; tileX++) {
							/*
							double target_disparity = corr2d.restoreMlTilePixel(
									tileX,                            // int         tileX,
									tileY,                            // int         tileY,
									ml_hwidth,                        // int         ml_hwidth,
									ml_data,                          // double [][] ml_data,
									ImageDtt.ML_OTHER_INDEX,          // int         ml_layer,
									ImageDtt.ML_OTHER_TARGET ,        // int         ml_index,
									tilesX);                          // int         tilesX);
							double gtruth_disparity = corr2d.restoreMlTilePixel(
									tileX,                            // int         tileX,
									tileY,                            // int         tileY,
									ml_hwidth,                        // int         ml_hwidth,
									ml_data,                          // double [][] ml_data,
									ImageDtt.ML_OTHER_INDEX,          // int         ml_layer,
									ImageDtt.ML_OTHER_GTRUTH ,        // int         ml_index,
									tilesX);                          // int         tilesX);
							double gtruth_strength = corr2d.restoreMlTilePixel(
									tileX,                            // int         tileX,
									tileY,                            // int         tileY,
									ml_hwidth,                        // int         ml_hwidth,
									ml_data,                          // double [][] ml_data,
									ImageDtt.ML_OTHER_INDEX,          // int         ml_layer,
									ImageDtt.ML_OTHER_GTRUTH_STRENGTH ,     // int         ml_index,
									tilesX);                          // int         tilesX);
							// converting disparity to 9.7 ( 1/128 pixel step, +/-256 pixels disparity range), 0x8000 - zero disparity
							// converting strength to 2 bytes 0.16 fixed point
							int itd = (int) Math.round(128 * target_disparity) + 0x8000;
							int [] itarget_disparity = {itd >> 8, itd & 0xff};
							int igt = (int) Math.round(128 * gtruth_disparity) + 0x8000;
							int [] igtruth_disparity = {igt >> 8, igt & 0xff};
							int igs = (int) Math.round(0x10000 * gtruth_strength);
							int [] igtruth_strength =  {igs >> 8, igs & 0xff};
							for (int nb = 0; nb<2; nb++) {
								if (!Double.isNaN(target_disparity)) {
									int indx =  corr2d.getMlTilePixelIndex(tileX,tileY, ml_hwidth, ImageDtt.ML_OTHER_TARGET + nb, tilesX);
									iml_data[nl][indx] = (byte) itarget_disparity[nb];
								}
								if (!Double.isNaN(gtruth_disparity)) {
									int indx =  corr2d.getMlTilePixelIndex(tileX,tileY, ml_hwidth, ImageDtt.ML_OTHER_GTRUTH + nb, tilesX);
									iml_data[nl][indx] = (byte) igtruth_disparity[nb];
								}
								if (gtruth_strength > 0.0) {
									int indx =  corr2d.getMlTilePixelIndex(tileX,tileY, ml_hwidth, ImageDtt.ML_OTHER_GTRUTH_STRENGTH + nb, tilesX);
									iml_data[nl][indx] = (byte) igtruth_strength[nb];
								}
							}
							*/
							for (int data:signed_data) if (aux_mode || (data <= ImageDtt.ML_OTHER_GTRUTH_STRENGTH)) {
								double d =  corr2d.restoreMlTilePixel(
										tileX,                            // int         tileX,
										tileY,                            // int         tileY,
										ml_hwidth,                        // int         ml_hwidth,
										ml_data,                          // double [][] ml_data,
										ImageDtt.ML_OTHER_INDEX,          // int         ml_layer,
										data ,                           // int         ml_index,
										tilesX);                          // int         tilesX);
								if (!Double.isNaN(d)) {
									int id = (int) Math.round(128 * d) + 0x8000;
									int [] ida = {id >> 8, id & 0xff};
									for (int nb = 0; nb<2; nb++) {
										int indx =  corr2d.getMlTilePixelIndex(tileX,tileY, ml_hwidth, data + nb, tilesX);
										iml_data[nl][indx] = (byte) ida[nb];
									}
								}
							}
							for (int data:unsigned_data) if (aux_mode || (data <= ImageDtt.ML_OTHER_GTRUTH_STRENGTH)) {
								double d =  corr2d.restoreMlTilePixel(
										tileX,                            // int         tileX,
										tileY,                            // int         tileY,
										ml_hwidth,                        // int         ml_hwidth,
										ml_data,                          // double [][] ml_data,
										ImageDtt.ML_OTHER_INDEX,          // int         ml_layer,
										data ,                           // int         ml_index,
										tilesX);                          // int         tilesX);
								if (!Double.isNaN(d) && (d > 0.0)) {
									int iscale = 0x10000;
									if ((data == ImageDtt.ML_OTHER_GTRUTH_RMS) || (data == ImageDtt.ML_OTHER_GTRUTH_RMS_SPLIT)) {
										iscale = 0x1000;
									}
									int id = (int) Math.round(iscale * d);
									int [] ida = {id >> 8, id & 0xff};
									for (int nb = 0; nb<2; nb++) {
										int indx =  corr2d.getMlTilePixelIndex(tileX,tileY, ml_hwidth, data + nb, tilesX);
										iml_data[nl][indx] = (byte) ida[nb];
									}
								}
							}
						}
					}
				} else {
					double k = 254.0/(soft_mx-soft_mn);
					for (int i = 0; i < ml_data[nl].length;i++) {
						if (Double.isNaN(ml_data[nl][i])){
							iml_data[nl][i] = 0; // -128;
						} else {
							int iv = (int) Math.round(k*(ml_data[nl][i]-soft_mn));
							if      (iv < 0) iv = 0;
							else if (iv > 254) iv = 254;
							iml_data[nl][i] = (byte) (iv + 1); //  (iv - 127); // NaN will stay 0;
						}
					}
				}
			}
			for (int nl = 0; nl< ml_data.length; nl++) if (!skip_layers[nl]) {
				array_stack.addSlice(ImageDtt.ML_TITLES[nl], iml_data[nl]);
			}
		} else {
			float [] fpixels;
			for (int nl = 0; nl< ml_data.length; nl++) if (!skip_layers[nl]) {
				fpixels=new float[ml_data[nl].length];
				for (int j=0;j<fpixels.length;j++) fpixels[j]=(float) ml_data[nl][j];
				array_stack.addSlice(ImageDtt.ML_TITLES[nl],    fpixels);
			}
		}

		double disparityRadiusMain =  has_main ? quadCLT_main.geometryCorrection.getDisparityRadius():Double.NaN;
		double disparityRadiusAux =   quadCLT_aux.geometryCorrection.getDisparityRadius();
		double intercameraBaseline =  has_main ? quadCLT_aux.geometryCorrection.getBaseline():Double.NaN;

		ImagePlus imp_ml = new ImagePlus(title, array_stack);
		imp_ml.setProperty("VERSION",  "1.2");
		imp_ml.setProperty("tileWidth",   ""+ml_width);
		imp_ml.setProperty("dispOffset",  ""+disp_offset_low);
		if (disp_offset_high>disp_offset_low) {
			imp_ml.setProperty("dispOffsetLow",  ""+disp_offset_low);
			imp_ml.setProperty("dispOffsetHigh",  ""+disp_offset_high);
		}
		imp_ml.setProperty("ML_OTHER_TARGET",           ""+ImageDtt.ML_OTHER_TARGET);
		imp_ml.setProperty("ML_OTHER_GTRUTH",           ""+ImageDtt.ML_OTHER_GTRUTH);
		imp_ml.setProperty("ML_OTHER_GTRUTH_STRENGTH",  ""+ImageDtt.ML_OTHER_GTRUTH_STRENGTH);
		if (aux_mode) {
			imp_ml.setProperty("ML_OTHER_GTRUTH_RMS",       ImageDtt.ML_OTHER_GTRUTH_RMS);
			imp_ml.setProperty("ML_OTHER_GTRUTH_RMS_SPLIT", ImageDtt.ML_OTHER_GTRUTH_RMS_SPLIT);
			imp_ml.setProperty("ML_OTHER_GTRUTH_FG_DISP",   ImageDtt.ML_OTHER_GTRUTH_FG_DISP);
			imp_ml.setProperty("ML_OTHER_GTRUTH_FG_STR",    ImageDtt.ML_OTHER_GTRUTH_FG_STR);
			imp_ml.setProperty("ML_OTHER_GTRUTH_BG_DISP",   ImageDtt.ML_OTHER_GTRUTH_BG_DISP);
			imp_ml.setProperty("ML_OTHER_GTRUTH_BG_STR",    ImageDtt.ML_OTHER_GTRUTH_BG_STR);
			imp_ml.setProperty("ML_OTHER_AUX_DISP",         ImageDtt.ML_OTHER_AUX_DISP);
			imp_ml.setProperty("ML_OTHER_AUX_STR",          ImageDtt.ML_OTHER_AUX_STR);
		}
		if (has_main) {
			imp_ml.setProperty("disparityRadiusMain",  ""+disparityRadiusMain);
			imp_ml.setProperty("intercameraBaseline",  ""+intercameraBaseline);
		}
		imp_ml.setProperty("disparityRadiusAux",  ""+disparityRadiusAux);
		if (use8bpp) {
			imp_ml.setProperty("data_min",  ""+soft_mn);
			imp_ml.setProperty("data_max",  ""+soft_mx);
		}

		imp_ml.setProperty("comment_tileWidth",   "Square tile size for each 2d correlation, always odd");
		imp_ml.setProperty("comment_dispOffset",  "Tile target disparity minus ground truth disparity");
		if (disp_offset_high>disp_offset_low) {
			imp_ml.setProperty("comment_dispOffsetLow",  "Tile target disparity minus ground truth disparity, low margin for random distribution");
			imp_ml.setProperty("comment_dispOffsetHigh", "Tile target disparity minus ground truth disparity, high margin for random distribution");
		}
		imp_ml.setProperty("comment_ML_OTHER_TARGET",  "Offset of the target disparity in the \"other\" layer tile");
		imp_ml.setProperty("comment_ML_OTHER_GTRUTH",  "Offset of the ground truth disparity in the \"other\" layer tile");
		imp_ml.setProperty("comment_ML_OTHER_GTRUTH_STRENGTH",  "Offset of the ground truth strength in the \"other\" layer tile");
		if (aux_mode) {
			imp_ml.setProperty("comment_ML_OTHER_GTRUTH_RMS","Offset of the GT disparity RMS");
			imp_ml.setProperty("comment_ML_OTHER_GTRUTH_RMS_SPLIT", "Offset of the GT disparity RMS combined from FG and BG");
			imp_ml.setProperty("comment_ML_OTHER_GTRUTH_FG_DISP", "Offset of the GT FG disparity");
			imp_ml.setProperty("comment_ML_OTHER_GTRUTH_FG_STR", "Offset of the GT FG strength");
			imp_ml.setProperty("comment_ML_OTHER_GTRUTH_BG_DISP", "Offset of the GT BG disparity");
			imp_ml.setProperty("comment_ML_OTHER_GTRUTH_BG_STR", "Offset of the GT BG strength");
		}
		if (use8bpp) {
			imp_ml.setProperty("comment_data_min",  "Defined only for 8bpp mode - value, corresponding to -127 (-128 is NaN)");
			imp_ml.setProperty("comment_data_max",  "Defined only for 8bpp mode - value, corresponding to +127 (-128 is NaN)");
		}
		if (has_main) {
			imp_ml.setProperty("comment_disparityRadiusMain",  "Side of the square where 4 main camera subcameras are located (mm)");
			imp_ml.setProperty("comment_intercameraBaseline",  "Horizontal distance between the main and the auxiliary camera centers (mm). Disparity is specified for the main camera");
		}
		imp_ml.setProperty("comment_disparityRadiusAux",  "Side of the square where 4 main camera subcameras are located (mm). Disparity is specified for the main camera");

		(new JP46_Reader_camera(false)).encodeProperiesToInfo(imp_ml);
		imp_ml.getProcessor().resetMinAndMax();
		if (show ) {
			imp_ml.show();
		}
		File dir = new File(ml_directory);
		if (!dir.exists()){
			dir.mkdirs();
			System.out.println("Created "+dir);
		}

		String path = ml_directory+=Prefs.getFileSeparator()+imp_ml.getTitle();
		FileSaver fs=new FileSaver(imp_ml);
		fs.saveAsTiff(path+".tiff");
		if (debugLevel > -4) {
			System.out.println("Saved ML data to "+path+".tiff");
		}
	}



	/**
	 * Get disparity from the DSI data of the main camera and re-measure with the rig
	 * @param scan - main camera DSI scan to get initial disparity and selection from
	 * @param quadCLT_main main camera QuadCLT instance (should have tp initialized)
	 * @param quadCLT_aux auxiliary camera QuadCLT instance (should have tp initialized)
	 * @param clt_parameters various configuration parameters
	 * @param threadsMax maximal number of threads to use
	 * @param updateStatus update IJ status bar
	 * @param debugLevel debug level
	 * @return rig measurement results
	 */

	public double [][] setBimapFromCLTPass3d(
			CLTPass3d                                      scan,
			QuadCLT                                        quadCLT_main,  // tiles should be set
			QuadCLT                                        quadCLT_aux,
			CLTParameters       clt_parameters,
			final int                                      threadsMax,  // maximal number of threads to launch
			final boolean                                  updateStatus,
			final int                                      debugLevel){
		int tile_op_all = clt_parameters.tile_task_op; //FIXME Use some constant?
		final int tilesX = quadCLT_main.tp.getTilesX();
		final int tilesY = quadCLT_main.tp.getTilesY();
		int [][]                                           tile_op = new int[tilesY][tilesX]; // common for both amin and aux
		double [][]                                        disparity_array = new double[tilesY][tilesX];

		double [] disparity =    scan.getDisparity();
		double [] strength =     scan.getOriginalStrength();
		//		  boolean [] selection =   scan.getSelected();
		boolean [] selection =   new boolean [strength.length];
		for (int nTile = 0; nTile < selection.length; nTile++) {
			selection[nTile] = strength[nTile] > 0.0;
		}

		for (int nTile = 0; nTile < disparity.length; nTile++) {
			if (((selection == null) || (selection[nTile]) && !Double.isNaN(disparity[nTile]))) {
				int tileY = nTile / tilesX;
				int tileX = nTile % tilesX;
				tile_op[tileY][tileX] = tile_op_all;
				disparity_array[tileY][tileX] = disparity[nTile];
			}
		}
		double [][] disparity_bimap = measureRig(
				quadCLT_main,    // QuadCLT                                  quadCLT_main,  // tiles should be set
				quadCLT_aux,     //QuadCLT                                   quadCLT_aux,
				tile_op,         // int [][]                                 tile_op, // common for both amin and aux
				disparity_array, // double [][]                              disparity_array,
				null, // double [][]                                    ml_data,         // data for ML - 10 layers - 4 center areas (3x3, 5x5,..) per camera-per direction, 1 - composite, and 1 with just 1 data (target disparity)
				clt_parameters,  // EyesisCorrectionParameters.CLTParameters clt_parameters,
				clt_parameters.getFatZero(quadCLT_main.isMonochrome()), // double                                         fatzero,
				false,                   //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
				0,                          // final int                                      lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
				clt_parameters.rig.no_int_x0, // final boolean             no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very
				threadsMax,      // final int        threadsMax,  // maximal number of threads to launch
				updateStatus,    // final boolean    updateStatus,
				debugLevel);     // final int        debugLevel);

		return disparity_bimap;
	}

	/**
	 * Perform rig measurement from the new predicted disparity values
	 * @param disparity array of predicted disparities, NaN - do not measure
	 * @param quadCLT_main main camera QuadCLT instance (should have tp initialized)
	 * @param quadCLT_aux auxiliary camera QuadCLT instance (should have tp initialized)
	 * @param clt_parameters various configuration parameters
	 * @param threadsMax maximal number of threads to use
	 * @param updateStatus update IJ status bar
	 * @param debugLevel debug level
	 * @return rig measurement results
	 */

	public double [][] setBimapFromDisparityNaN(
			double []                                      disparity,
			QuadCLT                                        quadCLT_main,  // tiles should be set
			QuadCLT                                        quadCLT_aux,
			CLTParameters       clt_parameters,
			boolean                                        notch_mode,      // use notch filter for inter-camera correlation to detect poles
			int                                            lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using (2*notch_mode+1)^2 square
			boolean                                        no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
			final int                                      threadsMax,  // maximal number of threads to launch
			final boolean                                  updateStatus,
			final int                                      debugLevel){
		int tile_op_all = clt_parameters.tile_task_op; //FIXME Use some constant?
		final int tilesX = quadCLT_main.tp.getTilesX();
		final int tilesY = quadCLT_main.tp.getTilesY();
		int [][]                                           tile_op = new int[tilesY][tilesX]; // common for both amin and aux
		double [][]                                        disparity_array = new double[tilesY][tilesX];


		for (int nTile = 0; nTile < disparity.length; nTile++) {
			if (!Double.isNaN(disparity[nTile])) {
				int tileY = nTile / tilesX;
				int tileX = nTile % tilesX;
				tile_op[tileY][tileX] = tile_op_all;
				disparity_array[tileY][tileX] = disparity[nTile];
			}
		}
		double [][] disparity_bimap = measureRig(
				quadCLT_main,    // QuadCLT                                  quadCLT_main,  // tiles should be set
				quadCLT_aux,     //QuadCLT                                   quadCLT_aux,
				tile_op,         // int [][]                                 tile_op, // common for both amin and aux
				disparity_array, // double [][]                              disparity_array,
				null, // double [][]                                    ml_data,         // data for ML - 10 layers - 4 center areas (3x3, 5x5,..) per camera-per direction, 1 - composite, and 1 with just 1 data (target disparity)
				clt_parameters,  // EyesisCorrectionParameters.CLTParameters clt_parameters,
				clt_parameters.getFatZero(quadCLT_main.isMonochrome()), // double                                         fatzero,
				notch_mode,                          //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
				lt_rad,                               // final int                                      lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
				no_int_x0,       // boolean                                  no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
				threadsMax,      // final int        threadsMax,  // maximal number of threads to launch
				updateStatus,    // final boolean    updateStatus,
				debugLevel);     // final int        debugLevel);

		return disparity_bimap;
	}

	/**
	 * Select tiles that have small enough residual disparity to trust on each of the two cameras and the whole rig.
	 * @param quadCLT_main main camera QuadCLT instance (should have tp initialized)
	 * @param quadCLT_aux auxiliary camera QuadCLT instance (should have tp initialized)
	 * @param min_combo_strength minimal combined correlation strength for each camera and the rig
	 * @param max_trusted_disparity maximal disparity on a rig to trust (currently 4.0 pixels from zero each way)
	 * @param trusted_tolerance allow other cameras with scaled disparity if their disparity is under this value (for small baselines)
	 * @param bimap measured data
	 * @return per-tile array of trusted tiles
	 */
	boolean [] getTrustedDisparity(
			QuadCLT            quadCLT_main,  // tiles should be set
			QuadCLT            quadCLT_aux,
			boolean            use_individual,
			double             min_combo_strength,    // check correlation strength combined for all 3 correlations
			double             max_trusted_disparity,
			double             trusted_tolerance,
			boolean []         was_trusted,
			double [][]        bimap // current state of measurements
			) {
		double trusted_inter =    max_trusted_disparity;
		double trusted_main =     trusted_inter * quadCLT_main.geometryCorrection.getDisparityRadius()/quadCLT_aux.geometryCorrection.getBaseline();
		double trusted_aux =      trusted_inter * quadCLT_aux.geometryCorrection.getDisparityRadius()/quadCLT_aux.geometryCorrection.getBaseline();
		if (trusted_main < trusted_tolerance) trusted_main = trusted_tolerance;
		if (trusted_aux  < trusted_tolerance) trusted_aux =  trusted_tolerance;
		boolean [] trusted = new boolean [bimap[ImageDtt.BI_DISP_CROSS_INDEX].length];
		if (use_individual) {
			for (int i = 0; i < trusted.length; i++) {
				trusted[i] = (Math.abs(bimap[ImageDtt.BI_DISP_CROSS_INDEX][i]) <= trusted_inter) &&
						(Math.abs(bimap[ImageDtt.BI_DISP_FULL_INDEX][i])  <= trusted_main) &&
						(Math.abs(bimap[ImageDtt.BI_ADISP_FULL_INDEX][i]) <= trusted_aux) &&
						(bimap[ImageDtt.BI_STR_ALL_INDEX][i]              >= min_combo_strength) &&
						((was_trusted == null) || was_trusted[i]);
			}
		} else {
			for (int i = 0; i < trusted.length; i++) {
				trusted[i] = (Math.abs(bimap[ImageDtt.BI_DISP_CROSS_INDEX][i]) <= trusted_inter) &&
						(bimap[ImageDtt.BI_STR_CROSS_INDEX][i]                 >= min_combo_strength) &&
						((was_trusted == null) || was_trusted[i]);
			}
		}
		return trusted;
	}
	/**
	 * Select tiles that have small enough residual disparity and sufficient strength to trust for the whole rig only
	 * @param min_inter_strength minimal inter-camera correlation strength
	 * @param max_trusted_disparity maximal disparity on a rig to trust (currently 4.0 pixels from zero each way)
	 * @param bimap measured data
	 * @return per-tile array of trusted tiles
	 */
	boolean [] getTrustedDisparityInter(
			double             min_inter_strength,    // check correlation strength combined for all 3 correlations
			double             max_trusted_disparity,
			boolean []         was_trusted,
			double [][]        bimap // current state of measurements
			) {
		double trusted_inter =    max_trusted_disparity;
		boolean [] trusted = new boolean [bimap[ImageDtt.BI_DISP_CROSS_INDEX].length];
		for (int i = 0; i < trusted.length; i++) {
			trusted[i] = (Math.abs(bimap[ImageDtt.BI_DISP_CROSS_INDEX][i]) <= trusted_inter) &&
					(bimap[ImageDtt.BI_STR_CROSS_INDEX][i]              >= min_inter_strength) &&
					((was_trusted == null) || was_trusted[i]);
		}
		return trusted;
	}

	/**
	 * Refine (re-measure with updated expected disparity) tiles. If refine_min_strength and refine_tolerance are both
	 * set to 0.0, all (or listed) tiles will be re-measured, use camera after extrinsics are changed
	 * With refine_min_strength == 0, will re-measure infinity (have keep_inf == true)
	 * @param quadCLT_main main camera QuadCLT instance (should have tp initialized)
	 * @param quadCLT_aux auxiliary camera QuadCLT instance (should have tp initialized)
	 * @param src_bimap results of the older measurements (now includes expected disparity)
	 * @param prev_bimap results of the even older measurements to interpolate if there was an overshoot
	 * @param refine_mode reference camera data: 0 - main camera, 1 - aux camera, 2 - cross-camera
	 * @param keep_inf do not refine expected disparity for infinity, unless  refine_min_strength == 0
	 * @param inf_sel if not null - selects "infinity" tiles that shoild not be remeasured
	 * @param refine_min_strength do not refine weaker tiles
	 * @param refine_tolerance do not refine if residual disparity (after FD pre-shift by expected disparity) less than this
	 * @param tile_list list of selected tiles or null. If null - try to refine all tiles, otherwise - only listed tiles
	 * @param num_new - otional int [1] to return number of new tiles
	 * @param clt_parameters various configuration parameters
	 * @param threadsMax maximal number of threads to use
	 * @param updateStatus update IJ status bar
	 * @param debugLevel debug level
	 * @return results of the new measurements combined with the old results
	 */

	public double [][] refineRig(
			QuadCLT                                  quadCLT_main,  // tiles should be set
			QuadCLT                                  quadCLT_aux,
			double [][]                              src_bimap, // current state of measurements
			double [][]                              prev_bimap, // previous state of measurements or null
			double []                                scale_bad,
			int                                      refine_mode, // 0 - by main, 1 - by aux, 2 - by inter
			boolean                                  keep_inf,    // keep expected disparity 0.0 if it was so
			double                                   inf_disparity,
			double                                   refine_min_strength, // do not refine weaker tiles
			double                                   refine_tolerance,    // do not refine if absolute disparity below
			ArrayList<Integer>                       tile_list, // or null
			int     []                               num_new,
			CLTParameters clt_parameters,
			boolean                                  notch_mode,      // use notch filter for inter-camera correlation to detect poles
			int                                      lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using (2*notch_mode+1)^2 square
			final boolean                            no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
			final int                                threadsMax,  // maximal number of threads to launch
			final boolean                            updateStatus,
			final int                                debugLevel){
		boolean [] selection = null;
		if (tile_list != null) {
			selection = new boolean [quadCLT_main.tp.getTilesX() * quadCLT_main.tp.getTilesY()];
			for (int nTile:tile_list) selection[nTile] = true;
		}
		return refineRigSel(
				quadCLT_main,  // tiles should be set
				quadCLT_aux,
				src_bimap, // current state of measurements
				prev_bimap, // previous state of measurements or null
				scale_bad,  //
				refine_mode, // 0 - by main, 1 - by aux, 2 - by inter
				keep_inf,    // keep expected disparity 0.0 if it was so
				inf_disparity, // double                                   inf_disparity,
				refine_min_strength, // do not refine weaker tiles
				refine_tolerance,    // do not refine if absolute disparity below
				selection,
				num_new,
				clt_parameters,
				notch_mode,                          //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
				lt_rad,                              // final int  // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
				no_int_x0,           // boolean                                  no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
				threadsMax,  // maximal number of threads to launch
				updateStatus,
				debugLevel);
	}

	public double [][] refineRigSel(
			QuadCLT                                        quadCLT_main,  // tiles should be set
			QuadCLT                                        quadCLT_aux,
			double [][]                                    src_bimap, // current state of measurements
			double [][]                                    prev_bimap, // previous state of measurements or null
			double []                                      scale_bad,
			int                                            refine_mode, // 0 - by main, 1 - by aux, 2 - by inter
			boolean                                        keep_inf,    // keep expected disparity 0.0 if it was so
			double                                         inf_disparity,
			double                                         refine_min_strength, // do not refine weaker tiles
			double                                         refine_tolerance,    // do not refine if absolute disparity below
			boolean []                                     selection,
			int     []                                     num_new,
			CLTParameters       clt_parameters,
			final boolean                                  notch_mode,      // use notch filter for inter-camera correlation to detect poles
			final int                                      lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using (2*notch_mode+1)^2 square
			final boolean                                  no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
			final int        threadsMax,  // maximal number of threads to launch
			final boolean    updateStatus,
			final int        debugLevel){
		int tilesX =quadCLT_main.tp.getTilesX();
		int tilesY =quadCLT_main.tp.getTilesY();
		int [][] tile_op = new int [tilesY][tilesX];
		double [][] disparity_array = new double [tilesY][tilesX];
		double disp_scale_main =  1.0/clt_parameters.corr_magic_scale; // Is it needed?
		double disp_scale_aux =   disp_scale_main * quadCLT_main.geometryCorrection.getDisparityRadius()/quadCLT_aux.geometryCorrection.getDisparityRadius();
		double disp_scale_inter = disp_scale_main * quadCLT_main.geometryCorrection.getDisparityRadius()/quadCLT_aux.geometryCorrection.getBaseline();
		int tile_op_all = clt_parameters.tile_task_op; //FIXME Use some constant?
		int numMeas = 0;
		for (int tileY = 0; tileY<tilesY;tileY++) {
			for (int tileX = 0; tileX<tilesX;tileX++) {
				int nTile = tileY * tilesX + tileX;
				if (((selection == null) || selection[nTile]) && !Double.isNaN(src_bimap[ImageDtt.BI_TARGET_INDEX][nTile])) {
					if (prepRefineTile(
							(lt_rad > 0),
							clt_parameters, // EyesisCorrectionParameters.CLTParameters       clt_parameters,
							tile_op_all,    // int                                            tile_op_all,
							src_bimap, // double [][]                                     src_bimap, // current state of measurements
							prev_bimap, // double [][]                                    prev_bimap, // previous state of measurements or null
							scale_bad,  // double []                                      scale_bad,
							tile_op, // int [][]                                          tile_op, // common for both amin and aux
							disparity_array, // double [][]                                    disparity_array,
							refine_mode, // int                                            refine_mode, // 0 - by main, 1 - by aux, 2 - by inter
							keep_inf,    // boolean                                        keep_inf,    // keep expected disparity 0.0 if it was so
							inf_disparity, // double                                         inf_disparity,
							refine_min_strength, // double                                         refine_min_strength, // do not refine weaker tiles
							refine_tolerance,    // double                                         refine_tolerance,    // do not refine if absolute disparity below
							disp_scale_main,  // double                                         disp_scale_main,  // 1.0
							disp_scale_aux,   //double                                         disp_scale_aux,   // ~0.58
							disp_scale_inter, //double                                         disp_scale_inter, // ~4.86
							//								  scale_step,       // double                                         scale_step,  // scale for "unstable tiles"
							tileX, // int                                            tileX,
							tileY, // int                                            tileY,
							nTile )) numMeas++; //int                                            nTile
				}
			}
		}
		if (debugLevel >0) {
			System.out.println("refineRig() mode="+refine_mode+": Prepared "+numMeas+" to measure");
		}
		double [][] disparity_bimap  =  measureRig(
				quadCLT_main,        // QuadCLT                                        quadCLT_main,  // tiles should be set
				quadCLT_aux,         // QuadCLT                                        quadCLT_aux,
				tile_op,             // int [][]                                       tile_op, // common for both amin and aux
				disparity_array,     // double [][]                                    disparity_array,
				null, // double [][]                                    ml_data,         // data for ML - 10 layers - 4 center areas (3x3, 5x5,..) per camera-per direction, 1 - composite, and 1 with just 1 data (target disparity)
				clt_parameters,      // EyesisCorrectionParameters.CLTParameters       clt_parameters,
				clt_parameters.getFatZero(quadCLT_main.isMonochrome()), // double                                         fatzero,
				notch_mode,                          //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
				lt_rad,                              // final int  // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
				no_int_x0,                            // final boolean             no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
				threadsMax,          //final int        threadsMax,  // maximal number of threads to launch
				updateStatus,        // final boolean    updateStatus,
				debugLevel);          // final int        debugLevel)
		// combine with old results for tiles that were not re-measured

		for (int tileY = 0; tileY<tilesY;tileY++) {
			for (int tileX = 0; tileX<tilesX;tileX++) {
				int nTile = tileY * tilesX + tileX;
				if ((selection == null) || selection[nTile]) {
					if (Double.isNaN(disparity_bimap[ImageDtt.BI_TARGET_INDEX][nTile]) && !Double.isNaN(src_bimap[ImageDtt.BI_TARGET_INDEX][nTile])) {
						for (int i = 0; i < disparity_bimap.length; i++) {
							disparity_bimap[i][nTile] = src_bimap[i][nTile];
						}
					}
				}
			}
		}
		if (num_new != null) {
			num_new[0] = numMeas;
		}
		return disparity_bimap;
	}

	public double [][] refineRig_new(
			QuadCLT                                  quadCLT_main,  // tiles should be set
			QuadCLT                                  quadCLT_aux,
			double [][]                              src_bimap, // current state of measurements
			double [][]                              prev_bimap, // previous state of measurements or null
			double []                                scale_bad,
			int                                      refine_mode, // 0 - by main, 1 - by aux, 2 - by inter
			boolean                                  keep_inf,    // keep expected disparity 0.0 if it was so
			boolean []                               inf_sel,     // infinity selection
			double                                   refine_min_strength, // do not refine weaker tiles
			double                                   refine_tolerance,    // do not refine if absolute disparity below
			ArrayList<Integer>                       tile_list, // or null
			int     []                               num_new,
			CLTParameters clt_parameters,
			boolean                                  notch_mode,      // use notch filter for inter-camera correlation to detect poles
			int                                      lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using (2*notch_mode+1)^2 square
			final boolean                            no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
			final int                                threadsMax,  // maximal number of threads to launch
			final boolean                            updateStatus,
			final int                                debugLevel){
		boolean [] selection = null;
		if (tile_list != null) {
			selection = new boolean [quadCLT_main.tp.getTilesX() * quadCLT_main.tp.getTilesY()];
			if (inf_sel == null) {
				for (int nTile:tile_list) selection[nTile] = true;
			} else {
				for (int nTile:tile_list) selection[nTile] = !inf_sel[nTile];
			}
		} else if (inf_sel != null) {
			selection = new boolean [quadCLT_main.tp.getTilesX() * quadCLT_main.tp.getTilesY()];
			for (int nTile = 0; nTile < selection.length; nTile++) {
				selection[nTile] = !inf_sel[nTile];
			}

		}
		return refineRigSel_new(
				quadCLT_main,  // tiles should be set
				quadCLT_aux,
				src_bimap, // current state of measurements
				prev_bimap, // previous state of measurements or null
				scale_bad,  //
				refine_mode, // 0 - by main, 1 - by aux, 2 - by inter
				keep_inf,    // keep expected disparity 0.0 if it was so
				refine_min_strength, // do not refine weaker tiles
				refine_tolerance,    // do not refine if absolute disparity below
				selection,
				num_new,
				clt_parameters,
				notch_mode,                          //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
				lt_rad,                              // final int  // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
				no_int_x0,           // boolean                                  no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
				threadsMax,  // maximal number of threads to launch
				updateStatus,
				debugLevel);
	}

	public double [][] refineRigSel_new(
			QuadCLT                                        quadCLT_main,  // tiles should be set
			QuadCLT                                        quadCLT_aux,
			double [][]                                    src_bimap, // current state of measurements
			double [][]                                    prev_bimap, // previous state of measurements or null
			double []                                      scale_bad,
			int                                            refine_mode, // 0 - by main, 1 - by aux, 2 - by inter
			boolean                                        keep_inf,    // keep expected disparity 0.0 if it was so
			double                                         refine_min_strength, // do not refine weaker tiles
			double                                         refine_tolerance,    // do not refine if absolute disparity below
			boolean []                                     selection,
			int     []                                     num_new,
			CLTParameters       clt_parameters,
			final boolean                                  notch_mode,      // use notch filter for inter-camera correlation to detect poles
			final int                                      lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using (2*notch_mode+1)^2 square
			final boolean                                  no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
			final int        threadsMax,  // maximal number of threads to launch
			final boolean    updateStatus,
			final int        debugLevel){
		int tilesX =quadCLT_main.tp.getTilesX();
		int tilesY =quadCLT_main.tp.getTilesY();
		int [][] tile_op = new int [tilesY][tilesX];
		double [][] disparity_array = new double [tilesY][tilesX];
		double disp_scale_main =  1.0/clt_parameters.corr_magic_scale; // Is it needed?
		double disp_scale_aux =   disp_scale_main * quadCLT_main.geometryCorrection.getDisparityRadius()/quadCLT_aux.geometryCorrection.getDisparityRadius();
		double disp_scale_inter = disp_scale_main * quadCLT_main.geometryCorrection.getDisparityRadius()/quadCLT_aux.geometryCorrection.getBaseline();
		int tile_op_all = clt_parameters.tile_task_op; //FIXME Use some constant?
		int numMeas = 0;
		for (int tileY = 0; tileY<tilesY;tileY++) {
			for (int tileX = 0; tileX<tilesX;tileX++) {
				int nTile = tileY * tilesX + tileX;
				if (((selection == null) || selection[nTile]) && !Double.isNaN(src_bimap[ImageDtt.BI_TARGET_INDEX][nTile])) {
					if (prepRefineTile(
							(lt_rad > 0),
							clt_parameters, // EyesisCorrectionParameters.CLTParameters       clt_parameters,
							tile_op_all,    // int                                            tile_op_all,
							src_bimap, // double [][]                                     src_bimap, // current state of measurements
							prev_bimap, // double [][]                                    prev_bimap, // previous state of measurements or null
							scale_bad,  // double []                                      scale_bad,
							tile_op, // int [][]                                          tile_op, // common for both amin and aux
							disparity_array, // double [][]                                    disparity_array,
							refine_mode, // int                                            refine_mode, // 0 - by main, 1 - by aux, 2 - by inter
							keep_inf,    // boolean                                        keep_inf,    // keep expected disparity 0.0 if it was so
							0.0,             // double                                         inf_disparity,
							refine_min_strength, // double                                         refine_min_strength, // do not refine weaker tiles
							refine_tolerance,    // double                                         refine_tolerance,    // do not refine if absolute disparity below
							disp_scale_main,  // double                                         disp_scale_main,  // 1.0
							disp_scale_aux,   //double                                         disp_scale_aux,   // ~0.58
							disp_scale_inter, //double                                         disp_scale_inter, // ~4.86
							//								  scale_step,       // double                                         scale_step,  // scale for "unstable tiles"
							tileX, // int                                            tileX,
							tileY, // int                                            tileY,
							nTile )) numMeas++; //int                                            nTile
				}
			}
		}
		if (debugLevel >0) {
			System.out.println("refineRigSel() mode="+refine_mode+": Prepared "+numMeas+" to measure");
		}
		double [][] disparity_bimap  =  measureRig(
				quadCLT_main,        // QuadCLT                                        quadCLT_main,  // tiles should be set
				quadCLT_aux,         // QuadCLT                                        quadCLT_aux,
				tile_op,             // int [][]                                       tile_op, // common for both amin and aux
				disparity_array,     // double [][]                                    disparity_array,
				null, // double [][]                                    ml_data,         // data for ML - 10 layers - 4 center areas (3x3, 5x5,..) per camera-per direction, 1 - composite, and 1 with just 1 data (target disparity)
				clt_parameters,      // EyesisCorrectionParameters.CLTParameters       clt_parameters,
				clt_parameters.getFatZero(quadCLT_main.isMonochrome()), // double                                         fatzero,
				notch_mode,                          //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
				lt_rad,                              // final int  // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
				no_int_x0,                            // final boolean             no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
				threadsMax,          //final int        threadsMax,  // maximal number of threads to launch
				updateStatus,        // final boolean    updateStatus,
				debugLevel);          // final int        debugLevel)
		// combine with old results for tiles that were not re-measured

		for (int tileY = 0; tileY<tilesY;tileY++) {
			for (int tileX = 0; tileX<tilesX;tileX++) {
				int nTile = tileY * tilesX + tileX;
// FIXME: see if the parameter is needed
//				if ((selection == null) || selection[nTile]) {
					if (Double.isNaN(disparity_bimap[ImageDtt.BI_TARGET_INDEX][nTile]) && !Double.isNaN(src_bimap[ImageDtt.BI_TARGET_INDEX][nTile])) {
						if (nTile == 35509) {
							System.out.println("refineRigSel(): nTile="+nTile);
						}

						for (int i = 0; i < disparity_bimap.length; i++) {
							disparity_bimap[i][nTile] = src_bimap[i][nTile];
						}
					}
//				}
			}
		}
		if (num_new != null) {
			num_new[0] = numMeas;
		}
		return disparity_bimap;
	}

	/**
	 * Add measurements with new specified disparity of the main camera
	 * @param quadCLT_main main camera QuadCLT instance (should have tp initialized)
	 * @param quadCLT_aux auxiliary camera QuadCLT instance (should have tp initialized)
	 * @param src_bimap results of the older measurements (now includes expected disparity) or null (no old results available)
	 * @param disparity new expected disparity value to try
	 * @param tile_list list of selected tiles or null. If not null, will not re-measure listed tiles
	 * @param clt_parameters various configuration parameters
	 * @param threadsMax maximal number of threads to use
	 * @param updateStatus update IJ status bar
	 * @param debugLevel debug level
	 * @return results of the new measurements combined with the old results (if available)
	 */

	public double [][] measureNewRigDisparity(
			QuadCLT                                        quadCLT_main,  // tiles should be set
			QuadCLT                                        quadCLT_aux,
			double [][]                                    src_bimap, // current state of measurements (or null for new measurement)
			double                                         disparity,
			ArrayList<Integer>                             tile_list, // or null. If non-null - do not remeasure members of the list
			CLTParameters       clt_parameters,
			boolean                                        notch_mode,      // use notch filter for inter-camera correlation to detect poles
			int                                            lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using (2*notch_mode+1)^2 square
			boolean                                        no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
			final int        threadsMax,  // maximal number of threads to launch
			final boolean    updateStatus,
			final int        debugLevel){
		if ((debugLevel > 0) && (tile_list != null)) {
			System.out.println("measureNewRigDisparity(), disparity = "+disparity+", tile_list.size()="+tile_list.size());
		}
		int tile_op_all = clt_parameters.tile_task_op; //FIXME Use some constant?
		int tilesX =quadCLT_main.tp.getTilesX();
		int tilesY =quadCLT_main.tp.getTilesY();
		int [][] tile_op = new int [tilesY][tilesX];
		double [][] disparity_array = new double [tilesY][tilesX];

		boolean [] selected = new boolean [tilesX * tilesY];
		if (tile_list != null) {
			for (int nTile:tile_list) {
				selected[nTile] = true;
			}
		}
		for (int tileY = 0; tileY<tilesY;tileY++) {
			for (int tileX = 0; tileX<tilesX;tileX++) {
				int nTile = tileY * tilesX + tileX;
				if ((src_bimap == null) || Double.isNaN(src_bimap[ImageDtt.BI_TARGET_INDEX][nTile]) || !selected[nTile]) {
					tile_op[tileY][tileX] = tile_op_all;
					disparity_array[tileY][tileX] = disparity;
				}
			}
		}
		double [][] disparity_bimap  =  measureRig(
				quadCLT_main,        // QuadCLT                                        quadCLT_main,  // tiles should be set
				quadCLT_aux,         // QuadCLT                                        quadCLT_aux,
				tile_op,             // int [][]                                       tile_op, // common for both amin and aux
				disparity_array,     // double [][]                                    disparity_array,
				null, // double [][]                                    ml_data,         // data for ML - 10 layers - 4 center areas (3x3, 5x5,..) per camera-per direction, 1 - composite, and 1 with just 1 data (target disparity)
				clt_parameters,      // EyesisCorrectionParameters.CLTParameters       clt_parameters,
				clt_parameters.getFatZero(quadCLT_main.isMonochrome()), // double                                         fatzero,
				notch_mode,                          //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
				lt_rad,                              // final int  // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
				no_int_x0,           // boolean                                  no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
				threadsMax,          //final int        threadsMax,  // maximal number of threads to launch
				updateStatus,        // final boolean    updateStatus,
				debugLevel);          // final int        debugLevel)

		// combine with old results (if available) for the tiles that were not re-measured
		if (src_bimap != null) {
			for (int tileY = 0; tileY<tilesY;tileY++) {
				for (int tileX = 0; tileX<tilesX;tileX++) {
					int nTile = tileY * tilesX + tileX;
					if (nTile == 35509) {
						System.out.println("measureNewRigDisparity(): nTile="+nTile);
					}
					boolean use_old = false;
					//				  if (Double.isNaN(disparity_bimap[ImageDtt.BI_TARGET_INDEX][nTile]) && !Double.isNaN(src_bimap[ImageDtt.BI_TARGET_INDEX][nTile])) {
					if (Double.isNaN(disparity_bimap[ImageDtt.BI_TARGET_INDEX][nTile])) {
						use_old = true;
					} else if (!Double.isNaN(src_bimap[ImageDtt.BI_TARGET_INDEX][nTile])) {
						double comp_strength_old = src_bimap[ImageDtt.BI_STR_FULL_INDEX][nTile] * src_bimap[ImageDtt.BI_ASTR_FULL_INDEX][nTile];
						double comp_strength_new = disparity_bimap[ImageDtt.BI_STR_FULL_INDEX][nTile]* disparity_bimap[ImageDtt.BI_ASTR_FULL_INDEX][nTile];
						if (comp_strength_old > comp_strength_new) {
							use_old = true;
						} else if (!(comp_strength_old < comp_strength_new)) { //one or both are NaN
							// for now - enough,
							if (Double.isNaN (comp_strength_new)) {
								use_old = true;
							}
						}
					}
					if (use_old) {
						for (int i = 0; i < disparity_bimap.length; i++) {
							disparity_bimap[i][nTile] = src_bimap[i][nTile];
						}
					}
				}
			}
		}
		return disparity_bimap;
	}

	/**
	 * Measure with specified disparity array, skip Double.NaN tiles
	 * @param quadCLT_main main camera QuadCLT instance (should have tp initialized)
	 * @param quadCLT_aux auxiliary camera QuadCLT instance (should have tp initialized)
	 * @param disparity expected per-tile disparities. NaN - do not measure
	 * @param clt_parameters various configuration parameters
	 * @param threadsMax maximal number of threads to use
	 * @param updateStatus update IJ status bar
	 * @param debugLevel debug level
	 * @return results of the new measurements combined with the old results (if available)
	 */
	public double [][] measureNewRigDisparity(
			QuadCLT                                        quadCLT_main,  // tiles should be set
			QuadCLT                                        quadCLT_aux,
			double []                                      disparity, // Double.NaN - skip, ohers - measure
			CLTParameters       clt_parameters,
			boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
			int                 lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using (2*notch_mode+1)^2 square
			boolean                                  no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
			final int        threadsMax,  // maximal number of threads to launch
			final boolean    updateStatus,
			final int        debugLevel){
		int tile_op_all = clt_parameters.tile_task_op; //FIXME Use some constant?
		int tilesX =quadCLT_main.tp.getTilesX();
		int tilesY =quadCLT_main.tp.getTilesY();
		int [][] tile_op = new int [tilesY][tilesX];
		double [][] disparity_array = new double [tilesY][tilesX];
		for (int tileY = 0; tileY<tilesY;tileY++) {
			for (int tileX = 0; tileX<tilesX;tileX++) {
				int nTile = tileY * tilesX + tileX;
				disparity_array[tileY][tileX] = disparity[nTile];
				if (!Double.isNaN(disparity[nTile])) {
					tile_op[tileY][tileX] = tile_op_all;
					//disparity_array[tileY][tileX] = disparity[nTile];
				}
			}
		}
		double [][] disparity_bimap  =  measureRig(
				quadCLT_main,        // QuadCLT                                        quadCLT_main,  // tiles should be set
				quadCLT_aux,         // QuadCLT                                        quadCLT_aux,
				tile_op,             // int [][]                                       tile_op, // common for both amin and aux
				disparity_array,     // double [][]                                    disparity_array,
				null, // double [][]                                    ml_data,         // data for ML - 10 layers - 4 center areas (3x3, 5x5,..) per camera-per direction, 1 - composite, and 1 with just 1 data (target disparity)
				clt_parameters,      // EyesisCorrectionParameters.CLTParameters       clt_parameters,
				clt_parameters.getFatZero(quadCLT_main.isMonochrome()), // double                                         fatzero,
				notch_mode,          //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
				lt_rad,                              // final int  // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
				no_int_x0,           // boolean                                  no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
				threadsMax,          //final int        threadsMax,  // maximal number of threads to launch
				updateStatus,        // final boolean    updateStatus,
				debugLevel);         // final int        debugLevel)
		return disparity_bimap;
	}


	private boolean prepRefineTile(
			boolean                                        remeasure_nan, // use old suggestion if result of the measurement was NaN
			CLTParameters       clt_parameters,
			int                                            tile_op_all,
			double [][]                                    src_bimap, // current state of measurements
			double [][]                                    prev_bimap, // previous state of measurements or null
			double []                                      scale_bad,
			int [][]                                       tile_op, // common for both amin and aux
			double [][]                                    disparity_array,
			int                                            refine_mode, // 0 - by main, 1 - by aux, 2 - by inter, 3 - inter-dx
			boolean                                        keep_inf,    // keep expected disparity 0.0 if it was so
			double                                         inf_disparity,
			double                                         refine_min_strength, // do not refine weaker tiles
			double                                         refine_tolerance,    // do not refine if absolute disparity below
			double                                         disp_scale_main,  // 1.0
			double                                         disp_scale_aux,   // ~0.58
			double                                         disp_scale_inter, // ~4.86
			//			  double                                         scale_step,  // scale for "unstable tiles"
			int                                            tileX,
			int                                            tileY,
			int                                            nTile) {

		boolean debug_this = false; // nTile==40661; // 61924;
		// check if it was measured (skip NAN)
		if (Double.isNaN(src_bimap[ImageDtt.BI_TARGET_INDEX][nTile])) return false;
		// check if it is infinity and change is prohibited
		if (keep_inf && (src_bimap[ImageDtt.BI_TARGET_INDEX][nTile] == inf_disparity)) {
			if ((refine_min_strength == 0.0) || (refine_tolerance == 0.0)) {
				tile_op[tileY][tileX] = tile_op_all;
				disparity_array[tileY][tileX] = inf_disparity;
				return true;
			}
			return false;
		}
		double diff_disp, strength, disp_scale, diff_prev;
		switch (refine_mode) {
		case 0:
			diff_disp = src_bimap[ImageDtt.BI_DISP_FULL_INDEX][nTile];
			diff_prev= (prev_bimap == null)? Double.NaN:prev_bimap[ImageDtt.BI_DISP_FULL_INDEX][nTile];
			strength = src_bimap[ImageDtt.BI_STR_FULL_INDEX][nTile];
			disp_scale = disp_scale_main;
			break;
		case 1:
			diff_disp = src_bimap[ImageDtt.BI_ADISP_FULL_INDEX][nTile];
			diff_prev= (prev_bimap == null)? Double.NaN:prev_bimap[ImageDtt.BI_ADISP_FULL_INDEX][nTile];
			strength = src_bimap[ImageDtt.BI_ASTR_FULL_INDEX][nTile];
			disp_scale = disp_scale_aux;
			break;
		case 2:
			diff_disp = src_bimap[ImageDtt.BI_DISP_CROSS_INDEX][nTile];
			diff_prev= (prev_bimap == null)? Double.NaN:prev_bimap[ImageDtt.BI_DISP_CROSS_INDEX][nTile];
			strength = src_bimap[ImageDtt.BI_STR_CROSS_INDEX][nTile];
			disp_scale = disp_scale_inter;
			break;
		default: // case 3
		diff_disp = src_bimap[ImageDtt.BI_DISP_CROSS_DX_INDEX][nTile];
		diff_prev= (prev_bimap == null)? Double.NaN:prev_bimap[ImageDtt.BI_DISP_CROSS_DX_INDEX][nTile];
		strength = src_bimap[ImageDtt.BI_STR_CROSS_INDEX][nTile];
		disp_scale = disp_scale_inter;
		}
		if (Double.isNaN(diff_disp)) {
			if (remeasure_nan) {
				tile_op[tileY][tileX] = tile_op_all;
				disparity_array[tileY][tileX] = src_bimap[ImageDtt.BI_TARGET_INDEX][nTile]; // repeat previous - only makes sense with averaging
				return true;

			} else {
				return false;
			}
		}

		// strong enough?
		if (strength < refine_min_strength) return false;
		// residual disparity large enough to bother
		if (Math.abs(diff_disp) < refine_tolerance) return false;
		// or use extrapolate too?
		if (debug_this) {
			//			  System.out.println("disp_scale="+disp_scale);
			if (prev_bimap != null) {
				//				  System.out.println(
				///						  "prepRefineTile(): prev_bimap[ImageDtt.BI_TARGET_INDEX][nTile] = "+prev_bimap[ImageDtt.BI_TARGET_INDEX][nTile]+
				//						  ", src_bimap[ImageDtt.BI_TARGET_INDEX][nTile]="+src_bimap[ImageDtt.BI_TARGET_INDEX][nTile]+", diff_prev="+diff_prev+", diff_disp="+diff_disp);
			} else {
				//				  System.out.println("prepRefineTile():  src_bimap[ImageDtt.BI_TARGET_INDEX][nTile]="+src_bimap[ImageDtt.BI_TARGET_INDEX][nTile]+", diff_prev="+diff_prev+", diff_disp="+diff_disp);
			}
		}
		boolean from_prev = false;
		if ((scale_bad != null) && (Math.abs(diff_disp) > Math.abs(diff_prev))){
			scale_bad[nTile] *= 0.9; // reduce step
			from_prev = true;
		}
		double new_disp;
		if (Double.isNaN(diff_prev) || (diff_prev * diff_disp > 0)) {
			new_disp = src_bimap[ImageDtt.BI_TARGET_INDEX][nTile] + diff_disp*disp_scale;
			//			  if (debug_this) System.out.println(" >> 1 => disparity_array["+tileY+"]["+tileX+"]="+new_disp);
		}  else { // interpolate
			new_disp = (src_bimap[ImageDtt.BI_TARGET_INDEX][nTile]*diff_prev - prev_bimap[ImageDtt.BI_TARGET_INDEX][nTile]*diff_disp)/(diff_prev-diff_disp);
			//			  if (debug_this) System.out.println(" >> 2  => disparity_array["+tileY+"]["+tileX+"]="+new_disp);
		}
		double ref_target = from_prev ? prev_bimap[ImageDtt.BI_TARGET_INDEX][nTile] : src_bimap[ImageDtt.BI_TARGET_INDEX][nTile];
		if ((scale_bad != null) && (scale_bad[nTile] < 1.0)) {
			new_disp = ref_target * (1.0 - scale_bad[nTile]) +  new_disp * scale_bad[nTile];
			//			  if (debug_this) System.out.println(" scale_bad["+nTile+"]= "+scale_bad[nTile]+" , corrected new_disp="+new_disp+" (was "+ src_bimap[ImageDtt.BI_TARGET_INDEX][nTile]+
			//					  ", prev "+prev_bimap[ImageDtt.BI_TARGET_INDEX][nTile]+")");
		} else if (prev_bimap!=null){
			//			  if (debug_this) System.out.println("new_disp="+new_disp+" (was "+ src_bimap[ImageDtt.BI_TARGET_INDEX][nTile]+
			//					  ", prev "+prev_bimap[ImageDtt.BI_TARGET_INDEX][nTile]+")");
		}
		if (debug_this) System.out.println("prepRefineTile():target_diff "+ src_bimap[ImageDtt.BI_TARGET_INDEX][nTile]+","+diff_disp+","+scale_bad[nTile]);

		//		  if (Math.abs((new_disp - ref_target)/new_disp) < refine_tolerance) return false;
		if (Math.abs(new_disp - ref_target) < refine_tolerance) return false;
		disparity_array[tileY][tileX] = new_disp;
		tile_op[tileY][tileX] = tile_op_all;
		return true;
	}

	double [][] measureRig(
			QuadCLT                                  quadCLT_main,  // tiles should be set
			QuadCLT                                  quadCLT_aux,
			int [][]                                 tile_op, // common for both amin and aux
			double [][]                              disparity_array,
			double [][]                              ml_data,         // data for ML - 10 layers - 4 center areas (3x3, 5x5,..) per camera-per direction, 1 - composite, and 1 with just 1 data (target disparity)
			CLTParameters clt_parameters,
			double                                   fatzero,
			boolean                                  notch_mode, // use pole-detection mode for inter-camera correlation
			int                                      lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using (2*notch_mode+1)^2 square
			boolean                                  no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
			final int                                threadsMax,  // maximal number of threads to launch
			final boolean                            updateStatus,
			final int                                debugLevel){
		ImageDtt image_dtt = new ImageDtt(
				clt_parameters.transform_size,
				quadCLT_main.isMonochrome(),
				quadCLT_main.isLwir(),
				clt_parameters.getScaleStrength(false));

		double [][] disparity_bimap  = new double [ImageDtt.BIDISPARITY_TITLES.length][]; //[0] -residual disparity, [1] - orthogonal (just for debugging) last 4 - max pixel differences

//		int [][] woi_tops = {quadCLT_main.woi_tops,quadCLT_aux.woi_tops};
		image_dtt.clt_bi_quad (
				clt_parameters,                       // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
				fatzero,                              // final double              fatzero,         // May use correlation fat zero from 2 different parameters - fat_zero and rig.ml_fatzero
				notch_mode,                           //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
				lt_rad,                               // final int                                   lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
				no_int_x0,                            // final boolean             no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
				tile_op,                              // final int [][]            tile_op_main,    // [tilesY][tilesX] - what to do - 0 - nothing for this tile
				disparity_array,                      // final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
				quadCLT_main.image_data,              // final double [][][]       image_data_main, // first index - number of image in a quad
				quadCLT_aux.image_data,               // final double [][][]       image_data_aux,  // first index - number of image in a quad
				quadCLT_main.saturation_imp,          // final boolean [][]        saturation_main, // (near) saturated pixels or null
				quadCLT_aux.saturation_imp,           // final boolean [][]        saturation_aux,  // (near) saturated pixels or null
				// correlation results - combo will be for the correation between two quad cameras
				null,                                 // final double [][][][]     clt_corr_combo,  // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				disparity_bimap,                      // final double [][]    disparity_bimap, // [23][tilesY][tilesX]
				ml_data,                              // 	final double [][]         ml_data,         // data for ML - 10 layers - 4 center areas (3x3, 5x5,..) per camera-per direction, 1 - composite, and 1 with just 1 data (target disparity)
				null,                                 // final double [][][][]     texture_tiles_main, // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
				null,                                 // final double [][][][]     texture_tiles_aux,  // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
				quadCLT_main.tp.getTilesX()*image_dtt.transform_size, // final int                 width,

				quadCLT_main.getGeometryCorrection(), // final GeometryCorrection  geometryCorrection_main,
				quadCLT_aux.getGeometryCorrection(),  // final GeometryCorrection  geometryCorrection_aux,
				quadCLT_main.getCLTKernels(),         // final double [][][][][][] clt_kernels_main, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				quadCLT_aux.getCLTKernels(),          // final double [][][][][][] clt_kernels_aux,  // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				clt_parameters.corr_magic_scale,      // final double              corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
				false, // true,                                 // 	final boolean             keep_clt_data,
//				woi_tops,                             // final int [][]            woi_tops,
				null,                                 // final double [][][]       ers_delay,        // if not null - fill with tile center acquisition delay
				threadsMax,                           // final int                 threadsMax,  // maximal number of threads to launch
				debugLevel-2);                        // final int                 globalDebugLevel);
		return disparity_bimap;
	}



	public double [][] remeasureRigML(
			double                                         disparity_offset_low,
			double                                         disparity_offset_high,
			QuadCLT                                        quadCLT_main,  // tiles should be set
			QuadCLT                                        quadCLT_aux,
			double []                                      disparity_main, // main camera disparity to use - if null, calculate from the rig one
			double []                                      disparity,
			double []                                      strength,
			CLTParameters       clt_parameters,
			int                                            ml_hwidth,
			double                                         fatzero,
			int              lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using (2*notch_mode+1)^2 square
			final int        threadsMax,  // maximal number of threads to launch
			final boolean    updateStatus,
			final int        debugLevel){
		int tile_op_all = clt_parameters.tile_task_op; //FIXME Use some constant?
		final int tilesX = quadCLT_main.tp.getTilesX();
		final int tilesY = quadCLT_main.tp.getTilesY();
		//		  int ml_hwidth =      2;
		int ml_width = 2 * ml_hwidth + 1;

		double [][] ml_data  = new double [ImageDtt.ML_TITLES.length][tilesX * tilesY * ml_width * ml_width];
		for (int nLayer=0; nLayer < ml_data.length; nLayer++) {
			for (int nTile=0; nTile <  ml_data[nLayer].length; nTile++) {
				ml_data[nLayer][nTile] = Double.NaN;
			}
		}

		int [][]           tile_op = new int[tilesY][tilesX]; // common for both main and aux
		double [][]        disparity_array = new double[tilesY][tilesX];
		boolean [] selection =   new boolean [strength.length];
		if (disparity_main != null) {
			for (int nTile = 0; nTile < selection.length; nTile++) {
				//			selection[nTile] = strength[nTile] > 0.0;
				selection[nTile] = (strength[nTile] > 0.0) || !Double.isNaN(disparity_main[nTile]); // measuring correlation for clusters - around defined tiles
				if (selection[nTile] && Double.isNaN(disparity[nTile])) {
					disparity[nTile] = 0.0; // it will not be used
					strength[nTile]  = 0.0; // should already be set
				}
			}
		}
		if (!(disparity_offset_high > disparity_offset_low)) { // otherwise offset to the rig data is used
			for (int nTile = 0; nTile < selection.length; nTile++) {
				if (selection[nTile] && Double.isNaN(disparity[nTile])) {
					disparity[nTile] = 0.0; // it will not be used
					strength[nTile]  = 0.0; // should already be set
				}
			}
		}
		Random rnd = new Random(System.nanoTime());
		Correlation2d corr2d = new Correlation2d(
				clt_parameters.img_dtt,              // ImageDttParameters  imgdtt_params,
				clt_parameters.transform_size,             // int transform_size,
				2.0,                        //  double wndx_scale, // (wndy scale is always 1.0)
				quadCLT_main.isMonochrome(),
				(debugLevel > -1));   //   boolean debug)

		for (int nTile = 0; nTile < disparity.length; nTile++) {
			if ((selection == null) || selection[nTile]) {
				int tileY = nTile / tilesX;
				int tileX = nTile % tilesX;
				tile_op[tileY][tileX] = tile_op_all;
				double disparity_offset = disparity_offset_low;
				if (disparity_offset_high > disparity_offset) { // will not happen if disparity_offset_high is NaN
					disparity_offset = disparity_offset_low + (disparity_offset_high-disparity_offset_low)*rnd.nextDouble();
				}
				if ((disparity_main != null) && !Double.isNaN(disparity_main[nTile])) {
					// Use actuall disparity from the main camera
					disparity_offset = disparity_main[nTile] - disparity[nTile];
					if (!Double.isNaN(disparity_offset_high) && (disparity_offset_high > 0.0)) {
						disparity_offset += disparity_offset_high* (2 * rnd.nextDouble() - 1.0);
					}
				}
				disparity_array[tileY][tileX] = disparity[nTile]+disparity_offset;
				corr2d.saveMlTilePixel(
						tileX,                         // int         tileX,
						tileY,                         // int         tileY,
						ml_hwidth,                     // int         ml_hwidth,
						ml_data,                       // double [][] ml_data,
						ImageDtt.ML_OTHER_INDEX,       // int         ml_layer,
						ImageDtt.ML_OTHER_GTRUTH ,     // int         ml_index,
						disparity[nTile],              // target disparitydouble      ml_value,
						tilesX);                       // int         tilesX);
				corr2d.saveMlTilePixel(
						tileX,                         // int         tileX,
						tileY,                         // int         tileY,
						ml_hwidth,                     // int         ml_hwidth,
						ml_data,                       // double [][] ml_data,
						ImageDtt.ML_OTHER_INDEX,       // int         ml_layer,
						ImageDtt.ML_OTHER_GTRUTH_STRENGTH, // int         ml_index,
						strength[nTile],               // target disparitydouble      ml_value,
						tilesX);                       // int         tilesX);
			}
		}
		measureRig(
				quadCLT_main, // QuadCLT                                        quadCLT_main,  // tiles should be set
				quadCLT_aux, // QuadCLT                                        quadCLT_aux,
				tile_op, // int [][]                                       tile_op, // common for both amin and aux
				disparity_array, // double [][]                                    disparity_array,
				ml_data, // double [][]                                    ml_data,         // data for ML - 10 layers - 4 center areas (3x3, 5x5,..) per camera-per direction, 1 - composite, and 1 with just 1 data (target disparity)
				clt_parameters,   // EyesisCorrectionParameters.CLTParameters       clt_parameters,
				fatzero,          // double                                         fatzero,
				false,            //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
				lt_rad,                              // final int  // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
				true,             // whatever here. final boolean             no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very
				threadsMax,       // maximal number of threads to launch // final int        threadsMax,  // maximal number of threads to launch
				updateStatus,     // final boolean    updateStatus,
				debugLevel);      // final int        debugLevel)
		return ml_data;
	}

	double [][] measureAux(
//			QuadCLT                                  quadCLT_main,  // tiles should be set
			QuadCLT                                  quadCLT_aux,
			int [][]                                 tile_op, // common for both amin and aux
			double [][]                              disparity_array,
			double [][]                              ml_data,         // data for ML - 10 layers - 4 center areas (3x3, 5x5,..) per camera-per direction, 1 - composite, and 1 with just 1 data (target disparity)
			CLTParameters clt_parameters,
			double                                   fatzero,
			boolean                                  notch_mode, // use pole-detection mode for inter-camera correlation
			int                                      lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using (2*notch_mode+1)^2 square
			boolean                                  no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
			final int                                threadsMax,  // maximal number of threads to launch
			final boolean                            updateStatus,
			final int                                debugLevel){
		ImageDtt image_dtt = new ImageDtt(
				clt_parameters.transform_size,
				quadCLT_aux.isMonochrome(),
				quadCLT_aux.isLwir(),
				clt_parameters.getScaleStrength(true));
		double [][] disparity_bimap  = new double [ImageDtt.BIDISPARITY_TITLES.length][]; //[0] -residual disparity, [1] - orthogonal (just for debugging) last 4 - max pixel differences
		image_dtt.clt_bi_quad (
				clt_parameters,                       // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
				fatzero,                              // final double              fatzero,         // May use correlation fat zero from 2 different parameters - fat_zero and rig.ml_fatzero
				notch_mode,                           //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
				lt_rad,                               // final int                                   lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
				no_int_x0,                            // final boolean             no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
				tile_op,                              // final int [][]            tile_op_main,    // [tilesY][tilesX] - what to do - 0 - nothing for this tile
				disparity_array,                      // final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
				null, // quadCLT_main.image_data,              // final double [][][]       image_data_main, // first index - number of image in a quad
				quadCLT_aux.image_data,               // final double [][][]       image_data_aux,  // first index - number of image in a quad
				null, // quadCLT_main.saturation_imp,          // final boolean [][]        saturation_main, // (near) saturated pixels or null
				quadCLT_aux.saturation_imp,           // final boolean [][]        saturation_aux,  // (near) saturated pixels or null
				// correlation results - combo will be for the correation between two quad cameras
				null,                                 // final double [][][][]     clt_corr_combo,  // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				disparity_bimap,                      // final double [][]    disparity_bimap, // [23][tilesY][tilesX]
				ml_data,                              // 	final double [][]         ml_data,         // data for ML - 10 layers - 4 center areas (3x3, 5x5,..) per camera-per direction, 1 - composite, and 1 with just 1 data (target disparity)
				null,                                 // final double [][][][]     texture_tiles_main, // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
				null,                                 // final double [][][][]     texture_tiles_aux,  // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
//				quadCLT_main.tp.getTilesX()*image_dtt.transform_size, // final int                 width,
				quadCLT_aux.tp.getTilesX()*image_dtt.transform_size, // final int                 width,

				null, // quadCLT_main.getGeometryCorrection(), // final GeometryCorrection  geometryCorrection_main,
				quadCLT_aux.getGeometryCorrection(),  // final GeometryCorrection  geometryCorrection_aux,
				null, // quadCLT_main.getCLTKernels(),         // final double [][][][][][] clt_kernels_main, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				quadCLT_aux.getCLTKernels(),          // final double [][][][][][] clt_kernels_aux,  // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				clt_parameters.corr_magic_scale,      // final double              corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
				false, // true,                                 // 	final boolean             keep_clt_data,
				null,                                 // final double [][][]       ers_delay,        // if not null - fill with tile center acquisition delay
				threadsMax,                           // final int                 threadsMax,  // maximal number of threads to launch
				debugLevel-2);                        // final int                 globalDebugLevel);
		return disparity_bimap;
	}

	public double [][] remeasureAuxML(
			QuadCLT                                        quadCLT_aux,
			double []                                      disparity, // Double.NaN - do not measure
			double [][]                                    dsi_aux_from_main, // Main camera DSI converted into the coordinates of the AUX one (see QuadCLT.FGBG_TITLES)-
			CLTParameters       clt_parameters,
			int                                            ml_hwidth,
			double                                         fatzero,
			final int        threadsMax,  // maximal number of threads to launch
			final boolean    updateStatus,
			final int        debugLevel){
		int tile_op_all = clt_parameters.tile_task_op; //FIXME Use some constant?
		final int tilesX = quadCLT_aux.tp.getTilesX();
		final int tilesY = quadCLT_aux.tp.getTilesY();
		//		  int ml_hwidth =      2;
		int ml_width = 2 * ml_hwidth + 1;

		double [][] ml_data  = new double [ImageDtt.ML_TITLES.length][tilesX * tilesY * ml_width * ml_width];
		for (int nLayer=0; nLayer < ml_data.length; nLayer++) {
			for (int nTile=0; nTile <  ml_data[nLayer].length; nTile++) {
				ml_data[nLayer][nTile] = Double.NaN;
			}
		}

		int [][]           tile_op = new int[tilesY][tilesX]; // common for both main and aux
		double [][]        disparity_array = new double[tilesY][tilesX];
		Correlation2d corr2d = new Correlation2d(
				clt_parameters.img_dtt,              // ImageDttParameters  imgdtt_params,
				clt_parameters.transform_size,             // int transform_size,
				2.0,                        //  double wndx_scale, // (wndy scale is always 1.0)
				quadCLT_aux.isMonochrome(),
				(debugLevel > -1));   //   boolean debug)

		for (int nTile = 0; nTile < disparity.length; nTile++) if (!Double.isNaN(disparity[nTile])) {
				int tileY = nTile / tilesX;
				int tileX = nTile % tilesX;
				tile_op[tileY][tileX] = tile_op_all;
				disparity_array[tileY][tileX] = disparity[nTile];
				corr2d.saveMlTilePixel(
						tileX,                                            // int         tileX,
						tileY,                                            // int         tileY,
						ml_hwidth,                                        // int         ml_hwidth,
						ml_data,                                          // double [][] ml_data,
						ImageDtt.ML_OTHER_INDEX,                          // int         ml_layer,
						ImageDtt.ML_OTHER_GTRUTH ,                        // int         ml_index,
						dsi_aux_from_main[QuadCLT.FGBG_DISPARITY][nTile], // double      ml_value,
						tilesX);                                          // int         tilesX);
				corr2d.saveMlTilePixel(
						tileX,                                            // int         tileX,
						tileY,                                            // int         tileY,
						ml_hwidth,                                        // int         ml_hwidth,
						ml_data,                                          // double [][] ml_data,
						ImageDtt.ML_OTHER_INDEX,                          // int         ml_layer,
						ImageDtt.ML_OTHER_GTRUTH_STRENGTH,                // int         ml_index,
						dsi_aux_from_main[QuadCLT.FGBG_STRENGTH][nTile],  // double      ml_value,
						tilesX);                                          // int         tilesX);
				if (dsi_aux_from_main.length < (QuadCLT.FGBG_RMS -1)) continue;
				corr2d.saveMlTilePixel(
						tileX,                                            // int         tileX,
						tileY,                                            // int         tileY,
						ml_hwidth,                                        // int         ml_hwidth,
						ml_data,                                          // double [][] ml_data,
						ImageDtt.ML_OTHER_INDEX,                          // int         ml_layer,
						ImageDtt.ML_OTHER_GTRUTH_RMS,                     // int         ml_index,
						dsi_aux_from_main[QuadCLT.FGBG_RMS][nTile],       // double      ml_value,
						tilesX);                                          // int         tilesX);
				if (dsi_aux_from_main.length < (QuadCLT.FGBG_RMS_SPLIT -1)) continue;
				corr2d.saveMlTilePixel(
						tileX,                                            // int         tileX,
						tileY,                                            // int         tileY,
						ml_hwidth,                                        // int         ml_hwidth,
						ml_data,                                          // double [][] ml_data,
						ImageDtt.ML_OTHER_INDEX,                          // int         ml_layer,
						ImageDtt.ML_OTHER_GTRUTH_RMS_SPLIT,               // int         ml_index,
						dsi_aux_from_main[QuadCLT.FGBG_RMS_SPLIT][nTile], // double      ml_value,
						tilesX);                                          // int         tilesX);
				if (dsi_aux_from_main.length < (QuadCLT.FGBG_FG_DISP -1)) continue;
				corr2d.saveMlTilePixel(
						tileX,                                            // int         tileX,
						tileY,                                            // int         tileY,
						ml_hwidth,                                        // int         ml_hwidth,
						ml_data,                                          // double [][] ml_data,
						ImageDtt.ML_OTHER_INDEX,                          // int         ml_layer,
						ImageDtt.ML_OTHER_GTRUTH_FG_DISP,                 // int         ml_index,
						dsi_aux_from_main[QuadCLT.FGBG_FG_DISP][nTile],   // double      ml_value,
						tilesX);                                          // int         tilesX);
				if (dsi_aux_from_main.length < (QuadCLT.FGBG_FG_STR -1)) continue;
				corr2d.saveMlTilePixel(
						tileX,                                            // int         tileX,
						tileY,                                            // int         tileY,
						ml_hwidth,                                        // int         ml_hwidth,
						ml_data,                                          // double [][] ml_data,
						ImageDtt.ML_OTHER_INDEX,                          // int         ml_layer,
						ImageDtt.ML_OTHER_GTRUTH_FG_STR,                  // int         ml_index,
						dsi_aux_from_main[QuadCLT.FGBG_FG_STR][nTile],    // double      ml_value,
						tilesX);                                          // int         tilesX);
				if (dsi_aux_from_main.length < (QuadCLT.FGBG_BG_DISP -1)) continue;
				corr2d.saveMlTilePixel(
						tileX,                                            // int         tileX,
						tileY,                                            // int         tileY,
						ml_hwidth,                                        // int         ml_hwidth,
						ml_data,                                          // double [][] ml_data,
						ImageDtt.ML_OTHER_INDEX,                          // int         ml_layer,
						ImageDtt.ML_OTHER_GTRUTH_BG_DISP,                 // int         ml_index,
						dsi_aux_from_main[QuadCLT.FGBG_BG_DISP][nTile],   // double      ml_value,
						tilesX);                                          // int         tilesX);
				if (dsi_aux_from_main.length < (QuadCLT.FGBG_BG_STR -1)) continue;
				corr2d.saveMlTilePixel(
						tileX,                                            // int         tileX,
						tileY,                                            // int         tileY,
						ml_hwidth,                                        // int         ml_hwidth,
						ml_data,                                          // double [][] ml_data,
						ImageDtt.ML_OTHER_INDEX,                          // int         ml_layer,
						ImageDtt.ML_OTHER_GTRUTH_BG_STR,                  // int         ml_index,
						dsi_aux_from_main[QuadCLT.FGBG_BG_STR][nTile],    // double      ml_value,
						tilesX);                                          // int         tilesX);

				if (dsi_aux_from_main.length < (QuadCLT.FGBG_AUX_DISP -1)) continue;
				corr2d.saveMlTilePixel(
						tileX,                                            // int         tileX,
						tileY,                                            // int         tileY,
						ml_hwidth,                                        // int         ml_hwidth,
						ml_data,                                          // double [][] ml_data,
						ImageDtt.ML_OTHER_INDEX,                          // int         ml_layer,
						ImageDtt.ML_OTHER_AUX_DISP,                 // int         ml_index,
						dsi_aux_from_main[QuadCLT.FGBG_AUX_DISP][nTile],   // double      ml_value,
						tilesX);                                          // int         tilesX);
				if (dsi_aux_from_main.length < (QuadCLT.FGBG_AUX_STR -1)) continue;
				corr2d.saveMlTilePixel(
						tileX,                                            // int         tileX,
						tileY,                                            // int         tileY,
						ml_hwidth,                                        // int         ml_hwidth,
						ml_data,                                          // double [][] ml_data,
						ImageDtt.ML_OTHER_INDEX,                          // int         ml_layer,
						ImageDtt.ML_OTHER_AUX_STR,                  // int         ml_index,
						dsi_aux_from_main[QuadCLT.FGBG_AUX_STR][nTile],    // double      ml_value,
						tilesX);                                          // int         tilesX);
		}

		measureAux(
				quadCLT_aux, // QuadCLT                                        quadCLT_aux,
				tile_op, // int [][]                                       tile_op, // common for both amin and aux
				disparity_array, // double [][]                                    disparity_array,
				ml_data, // double [][]                                    ml_data,         // data for ML - 10 layers - 4 center areas (3x3, 5x5,..) per camera-per direction, 1 - composite, and 1 with just 1 data (target disparity)
				clt_parameters,   // EyesisCorrectionParameters.CLTParameters       clt_parameters,
				fatzero,          // double                                         fatzero,
				false,            //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
				0,                // lt_rad,    // final int  // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
				true,             // whatever here. final boolean             no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very
				threadsMax,       // maximal number of threads to launch // final int        threadsMax,  // maximal number of threads to launch
				updateStatus,     // final boolean    updateStatus,
				debugLevel);      // final int        debugLevel)
		return ml_data;
	}





	ArrayList<Integer> selectRigTiles(
			CLTParameters       clt_parameters,
			boolean select_infinity,
			boolean select_noninfinity,
			double                                   inf_disparity,
			double[][] disparity_bimap,
			int tilesX)
	{
		ArrayList<Integer> tilesList = new ArrayList<Integer>();
		int numTiles = disparity_bimap[ImageDtt.BI_STR_FULL_INDEX].length;
		if (select_infinity) {
			for (int nTile = 0; nTile < numTiles; nTile++) {
				if (    (disparity_bimap[ImageDtt.BI_TARGET_INDEX][nTile] == inf_disparity) && // expected disparity was 0.0 (infinity)
						(disparity_bimap[ImageDtt.BI_STR_FULL_INDEX][nTile] >= clt_parameters.rig.inf_min_strength_main) &&
						(disparity_bimap[ImageDtt.BI_ASTR_FULL_INDEX][nTile] >= clt_parameters.rig.inf_min_strength_aux) &&
						(disparity_bimap[ImageDtt.BI_STR_CROSS_INDEX][nTile] >= clt_parameters.rig.inf_min_strength_rig) &&

						(disparity_bimap[ImageDtt.BI_DISP_FULL_INDEX][nTile] <= clt_parameters.rig.inf_max_disp_main) &&
						(disparity_bimap[ImageDtt.BI_ADISP_FULL_INDEX][nTile] <= clt_parameters.rig.inf_max_disp_aux) &&
						(disparity_bimap[ImageDtt.BI_DISP_CROSS_INDEX][nTile] <= clt_parameters.rig.inf_max_disp_rig) &&

						(disparity_bimap[ImageDtt.BI_DISP_FULL_INDEX][nTile] >= -clt_parameters.rig.inf_max_disp_main) &&
						(disparity_bimap[ImageDtt.BI_ADISP_FULL_INDEX][nTile] >= -clt_parameters.rig.inf_max_disp_aux) &&
						(disparity_bimap[ImageDtt.BI_DISP_CROSS_INDEX][nTile] >= -clt_parameters.rig.inf_max_disp_rig * clt_parameters.rig.inf_neg_tolerance)) {
					tilesList.add(nTile);
				}
			}
		}
		if (select_noninfinity) {
			for (int nTile = 0; nTile < numTiles; nTile++) {
				if (    (disparity_bimap[ImageDtt.BI_TARGET_INDEX][nTile] > inf_disparity ) && // expected disparity was > 0.0 (not infinity)
						(disparity_bimap[ImageDtt.BI_STR_FULL_INDEX][nTile] >= clt_parameters.rig.inf_min_strength_main)&&
						(disparity_bimap[ImageDtt.BI_ASTR_FULL_INDEX][nTile] >= clt_parameters.rig.inf_min_strength_aux)&&
						(disparity_bimap[ImageDtt.BI_STR_CROSS_INDEX][nTile] >= clt_parameters.rig.inf_min_strength_rig)&&
						(Math.abs(disparity_bimap[ImageDtt.BI_DISP_FULL_INDEX][nTile]) <= clt_parameters.rig.near_max_disp_main)&&
						(Math.abs(disparity_bimap[ImageDtt.BI_ADISP_FULL_INDEX][nTile]) <= clt_parameters.rig.near_max_disp_aux)&&
						(Math.abs(disparity_bimap[ImageDtt.BI_DISP_CROSS_INDEX][nTile]) <= clt_parameters.rig.near_max_disp_rig)) {
					tilesList.add(nTile);
				}
			}
		}
		return tilesList;
	}

	ArrayList<Integer> selectRigTiles_new(
			CLTParameters       clt_parameters,
			boolean ignore_target,    // to select non-zero "infinity"
			boolean []  disabled,
			boolean select_infinity,
			boolean select_noninfinity,
			double[][] disparity_bimap,
			int tilesX)
	{
		ArrayList<Integer> tilesList = new ArrayList<Integer>();
		int numTiles = disparity_bimap[ImageDtt.BI_STR_FULL_INDEX].length;
		if (select_infinity) {
			for (int nTile = 0; nTile < numTiles; nTile++) {
				if (    ((disabled == null)  || !disabled[nTile]) &&
						!Double.isNaN(disparity_bimap[ImageDtt.BI_DISP_CROSS_DX_INDEX][nTile]) &&
						!Double.isNaN(disparity_bimap[ImageDtt.BI_DISP_CROSS_DY_INDEX][nTile]) &&
						((disparity_bimap[ImageDtt.BI_TARGET_INDEX][nTile] == 0.0 ) || ignore_target)&& // expected disparity was 0.0 (infinity)
						(disparity_bimap[ImageDtt.BI_STR_FULL_INDEX][nTile] >= clt_parameters.rig.inf_min_strength_main) &&
						(disparity_bimap[ImageDtt.BI_ASTR_FULL_INDEX][nTile] >= clt_parameters.rig.inf_min_strength_aux) &&
						(disparity_bimap[ImageDtt.BI_STR_CROSS_INDEX][nTile] >= clt_parameters.rig.inf_min_strength_rig) &&

						(disparity_bimap[ImageDtt.BI_DISP_FULL_INDEX][nTile] <= clt_parameters.rig.inf_max_disp_main) &&
						(disparity_bimap[ImageDtt.BI_ADISP_FULL_INDEX][nTile] <= clt_parameters.rig.inf_max_disp_aux) &&
						(disparity_bimap[ImageDtt.BI_DISP_CROSS_INDEX][nTile] <= clt_parameters.rig.inf_max_disp_rig) &&

						(disparity_bimap[ImageDtt.BI_DISP_FULL_INDEX][nTile] >= -clt_parameters.rig.inf_max_disp_main) &&
						(disparity_bimap[ImageDtt.BI_ADISP_FULL_INDEX][nTile] >= -clt_parameters.rig.inf_max_disp_aux) &&
						(disparity_bimap[ImageDtt.BI_DISP_CROSS_INDEX][nTile] >= -clt_parameters.rig.inf_max_disp_rig * clt_parameters.rig.inf_neg_tolerance)) {
					tilesList.add(nTile);
				}
			}
		}
		if (select_noninfinity) {
			for (int nTile = 0; nTile < numTiles; nTile++) {
				if (    ((disabled == null)  || !disabled[nTile]) &&
						!Double.isNaN(disparity_bimap[ImageDtt.BI_DISP_CROSS_DX_INDEX][nTile]) &&
						!Double.isNaN(disparity_bimap[ImageDtt.BI_DISP_CROSS_DY_INDEX][nTile]) &&
						(disparity_bimap[ImageDtt.BI_TARGET_INDEX][nTile] > 0.0 ) && // expected disparity was > 0.0 (not infinity)
						(disparity_bimap[ImageDtt.BI_STR_FULL_INDEX][nTile] >= clt_parameters.rig.inf_min_strength_main)&&
						(disparity_bimap[ImageDtt.BI_ASTR_FULL_INDEX][nTile] >= clt_parameters.rig.inf_min_strength_aux)&&
						(disparity_bimap[ImageDtt.BI_STR_CROSS_INDEX][nTile] >= clt_parameters.rig.inf_min_strength_rig)&&
						(Math.abs(disparity_bimap[ImageDtt.BI_DISP_FULL_INDEX][nTile]) <= clt_parameters.rig.near_max_disp_main)&&
						(Math.abs(disparity_bimap[ImageDtt.BI_ADISP_FULL_INDEX][nTile]) <= clt_parameters.rig.near_max_disp_aux)&&
						(Math.abs(disparity_bimap[ImageDtt.BI_DISP_CROSS_INDEX][nTile]) <= clt_parameters.rig.near_max_disp_rig)) {
					tilesList.add(nTile);
				}
			}
		}
		return tilesList;
	}



	public void showListedRigTiles(
			String title,
			ArrayList<Integer> tile_list,
			double[][] disparity_bimap,
			int tilesX ){
		int numTiles = disparity_bimap[ImageDtt.BI_STR_FULL_INDEX].length;
		String [] titles = {"disparity", "dx","dy","strength","target"};
		double [][] dbg_inf = new double [5][numTiles];
		for (int nTile = 0; nTile < numTiles; nTile++){
			dbg_inf[0][nTile] = Double.NaN;
			dbg_inf[1][nTile] = Double.NaN;
			dbg_inf[2][nTile] = Double.NaN;
			dbg_inf[3][nTile] = Double.NaN; // 0.0;
			dbg_inf[4][nTile] = Double.NaN; // 0.0;
		}
		for (int nTile:tile_list) {
			dbg_inf[0][nTile] = disparity_bimap[ImageDtt.BI_DISP_CROSS_INDEX][nTile];
			dbg_inf[1][nTile] = disparity_bimap[ImageDtt.BI_DISP_CROSS_DX_INDEX][nTile];
			dbg_inf[2][nTile] = disparity_bimap[ImageDtt.BI_DISP_CROSS_DY_INDEX][nTile];
			dbg_inf[3][nTile] = disparity_bimap[ImageDtt.BI_STR_CROSS_INDEX][nTile];
			dbg_inf[4][nTile] = disparity_bimap[ImageDtt.BI_TARGET_INDEX][nTile];
		}
		(new ShowDoubleFloatArrays()).showArrays(
				dbg_inf,
				tilesX,
				numTiles/tilesX,
				true,
				title,
				titles );
	}

	public void batchRig(
			QuadCLT                                              quadCLT_main,  // tiles should be set
			QuadCLT                                              quadCLT_aux,
			CLTParameters             clt_parameters,
			EyesisCorrectionParameters.DebayerParameters         debayerParameters,
			ColorProcParameters                                  colorProcParameters,
			ColorProcParameters                                  colorProcParameters_aux,
			CorrectionColorProc.ColorGainsParameters             channelGainParameters,
			EyesisCorrectionParameters.RGBParameters             rgbParameters,
			EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
			Properties                                           properties,
			final int        threadsMax,  // maximal number of threads to launch
			final boolean    updateStatus,
			final int        debugLevel)  throws Exception
	{
		//		  final boolean    batch_mode = clt_parameters.batch_run;
		final int        debugLevelInner=clt_parameters.batch_run? -2: debugLevel;
		this.startTime=System.nanoTime();
		String [] sourceFiles=quadCLT_main.correctionsParameters.getSourcePaths();
		QuadCLT.SetChannels [] set_channels_main = quadCLT_main.setChannels(debugLevel);
		QuadCLT.SetChannels [] set_channels_aux =  quadCLT_aux.setChannels(debugLevel);
		if ((set_channels_main == null) || (set_channels_main.length==0) || (set_channels_aux == null) || (set_channels_aux.length==0)) {
			System.out.println("No files to process (of "+sourceFiles.length+")");
			return;
		}
		double [] referenceExposures_main = quadCLT_main.eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		double [] referenceExposures_aux =  quadCLT_aux.eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		for (int nSet = 0; nSet < set_channels_main.length; nSet++){
			// check it is the same set for both cameras
			if (set_channels_aux.length <= nSet ) {
				throw new Exception ("Set names for cameras do not match: main camera: '"+set_channels_main[nSet].name()+"', aux. camera: nothing");
			}
			if (!set_channels_main[nSet].name().equals(set_channels_aux[nSet].name())) {
				throw new Exception ("Set names for cameras do not match: main camera: '"+set_channels_main[nSet].name()+"', aux. camera: '"+set_channels_main[nSet].name()+"'");
			}

			int [] channelFiles_main = set_channels_main[nSet].fileNumber();
			int [] channelFiles_aux =  set_channels_aux[nSet].fileNumber();
			boolean [][] saturation_imp_main = (clt_parameters.sat_level > 0.0)? new boolean[channelFiles_main.length][] : null;
			boolean [][] saturation_imp_aux =  (clt_parameters.sat_level > 0.0)? new boolean[channelFiles_main.length][] : null;
			double [] scaleExposures_main = new double[channelFiles_main.length];
			double [] scaleExposures_aux =  new double[channelFiles_main.length];
			if (updateStatus) IJ.showStatus("Conditioning main camera image set for "+quadCLT_main.image_name);
			ImagePlus [] imp_srcs_main = quadCLT_main.conditionImageSet(
					clt_parameters,                 // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					colorProcParameters,            //  ColorProcParameters                       colorProcParameters, //
					sourceFiles,                    // String []                                 sourceFiles,
					set_channels_main[nSet].name(), // String                                    set_name,
					referenceExposures_main,        // double []                                 referenceExposures,
					channelFiles_main,              // int []                                    channelFiles,
					scaleExposures_main,            //output  // double [] scaleExposures
					saturation_imp_main,            //output  // boolean [][]                              saturation_imp,
					threadsMax,                 // int                                       threadsMax,
					debugLevelInner); // int                                       debugLevel);
			if (updateStatus) IJ.showStatus("Conditioning aux camera image set for "+quadCLT_main.image_name);
			ImagePlus [] imp_srcs_aux = quadCLT_aux.conditionImageSet(
					clt_parameters,                 // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					colorProcParameters_aux,        //  ColorProcParameters                       colorProcParameters, //
					sourceFiles,                    // String []                                 sourceFiles,
					set_channels_aux[nSet].name(), // String                                    set_name,
					referenceExposures_aux,        // double []                                 referenceExposures,
					channelFiles_aux,              // int []                                    channelFiles,
					scaleExposures_aux,            //output  // double [] scaleExposures
					saturation_imp_aux,            //output  // boolean [][]                              saturation_imp,
					threadsMax,                 // int                                       threadsMax,
					debugLevelInner); // int                                       debugLevel);

			// optionally adjust main, aux, rig here
			// Early main camera adjustment, rig data is not available
			for (int num_adjust_main = 0; num_adjust_main < quadCLT_main.correctionsParameters.rig_batch_adjust_main; num_adjust_main++) {
				if (updateStatus) IJ.showStatus("Building basic  DSI for the main camera image set "+quadCLT_main.image_name+
						" (w/o rig), pass "+(num_adjust_main+1)+" of "+quadCLT_main.correctionsParameters.rig_batch_adjust_main);
				if (debugLevel > -5) {
					System.out.println("Building basic  DSI for the main camera image set "+quadCLT_main.image_name+
							" (w/o rig), pass "+(num_adjust_main+1)+" of "+quadCLT_main.correctionsParameters.rig_batch_adjust_main);
				}

				quadCLT_main.preExpandCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
						imp_srcs_main, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
						saturation_imp_main, // boolean [][] saturation_imp, // (near) saturated pixels or null
						clt_parameters,
						debayerParameters,
						colorProcParameters,
						rgbParameters,
						threadsMax,  // maximal number of threads to launch
						updateStatus,
						debugLevelInner);
				// adjust extrinsics here
				System.out.println("Adjust main extrinsics here");
				if (updateStatus) IJ.showStatus("Adjusting main camera image set for "+quadCLT_main.image_name+
						" (w/o rig), pass "+(num_adjust_main+1)+" of "+quadCLT_main.correctionsParameters.rig_batch_adjust_main);
				if (debugLevel > -5) {
					System.out.println("Adjusting main camera image set for "+quadCLT_main.image_name+
							" (w/o rig), pass "+(num_adjust_main+1)+" of "+quadCLT_main.correctionsParameters.rig_batch_adjust_main);
				}
				boolean ok = quadCLT_main.extrinsicsCLT(
						clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
						false, // adjust_poly,
						  -1.0, // double inf_min,
						  1.0,  // double inf_max,
						threadsMax,  //final int        threadsMax,  // maximal number of threads to launch
						updateStatus,// final boolean    updateStatus,
						debugLevelInner); // final int        debugLevel)
// clear memory for main
				quadCLT_main.tp.resetCLTPasses();
				if (!ok) break;

			}
			// Early aux camera adjustment, rig data is not available
			for (int num_adjust_aux = 0; num_adjust_aux < quadCLT_main.correctionsParameters.rig_batch_adjust_aux; num_adjust_aux++) {
				if (updateStatus) IJ.showStatus("Building basic DSI for the aux camera image set "+quadCLT_main.image_name+
						" (w/o rig), pass "+(num_adjust_aux+1)+" of "+quadCLT_main.correctionsParameters.rig_batch_adjust_aux);
				if (debugLevel > -5) {
					System.out.println("Building basic DSI for the aux camera image set "+quadCLT_main.image_name+
							" (w/o rig), pass "+(num_adjust_aux+1)+" of "+quadCLT_main.correctionsParameters.rig_batch_adjust_aux);
				}
				quadCLT_aux.preExpandCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
						imp_srcs_aux, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
						saturation_imp_aux, // boolean [][] saturation_imp, // (near) saturated pixels or null
						clt_parameters,
						debayerParameters,
						colorProcParameters_aux,
						rgbParameters,
						threadsMax,  // maximal number of threads to launch
						updateStatus,
						debugLevelInner);
				// adjust extrinsics here
				if (updateStatus) IJ.showStatus("Adjusting aux camera image set for "+quadCLT_main.image_name+
						" (w/o rig), pass "+(num_adjust_aux+1)+" of "+quadCLT_main.correctionsParameters.rig_batch_adjust_aux);
				if (debugLevel > -5) {
					System.out.println("Adjusting aux camera image set for "+quadCLT_main.image_name+
							" (w/o rig), pass "+(num_adjust_aux+1)+" of "+quadCLT_main.correctionsParameters.rig_batch_adjust_aux);
				}
			    double inf_min = -1.0;
			    double inf_max =  1.0;
			    if (num_adjust_aux >= (quadCLT_main.correctionsParameters.rig_batch_adjust_aux/2)) {
			        inf_min = -0.2;
			        inf_max = 0.05;
			    }
				
				boolean ok = quadCLT_aux.extrinsicsCLT(
						clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
						false,          // adjust_poly,
			            inf_min, // double inf_min,
			            inf_max,  // double inf_max,
						threadsMax,     //final int        threadsMax,  // maximal number of threads to launch
						updateStatus,   // final boolean    updateStatus,
						debugLevelInner);    // final int        debugLevel)
				// clear memory for aux
				quadCLT_aux.tp.resetCLTPasses();
				if (!ok) break;
			}
			// Early rig adjustment, main/aux are not adjusted with rig data
			for (int num_adjust_rig = 0; num_adjust_rig < quadCLT_main.correctionsParameters.rig_batch_adjust_rig; num_adjust_rig++) {
				if (updateStatus) IJ.showStatus("Adjusting dual camera rig infinity (early) for "+quadCLT_main.image_name+
						", pass "+(num_adjust_rig+1)+" of "+quadCLT_main.correctionsParameters.rig_batch_adjust_rig);
				if (debugLevel > -5) {
					System.out.println("Adjusting dual camera rig infinity (early) for "+quadCLT_main.image_name+
							", pass "+(num_adjust_rig+1)+" of "+quadCLT_main.correctionsParameters.rig_batch_adjust_rig);
				}

				processInfinityRig(
						quadCLT_main,               // QuadCLT                                        quadCLT_main,
						quadCLT_aux,                // QuadCLT                                        quadCLT_aux,
						imp_srcs_main,              // ImagePlus []                                   imp_quad_main,
						imp_srcs_aux,               // ImagePlus []                                   imp_quad_aux,
						saturation_imp_main,        // boolean [][]                                   saturation_main, // (near) saturated pixels or null
						saturation_imp_aux,         // boolean [][]                                   saturation_aux, // (near) saturated pixels or null
						clt_parameters,             // EyesisCorrectionParameters.CLTParameters       clt_parameters,
						threadsMax,                 // final int        threadsMax,  // maximal number of threads to launch
						updateStatus,               // final boolean    updateStatus,
						debugLevelInner-2);                // final int        debugLevel);
			}
			// Delete previous results if any camera adjustments were made
			if (    (quadCLT_main.correctionsParameters.rig_batch_adjust_main > 0) ||
					(quadCLT_main.correctionsParameters.rig_batch_adjust_aux > 0) ||
					(quadCLT_main.correctionsParameters.rig_batch_adjust_rig > 0)) {

				if (debugLevel > -5) {
					System.out.println("\n---- Field adjustments after the first round ----");
					quadCLT_main.showExtrinsicCorr("main");  // show_fine_corr("main");
					quadCLT_aux.showExtrinsicCorr("aux");    // show_fine_corr("aux");
					quadCLT_aux.geometryCorrection.showRig();// show_fine_corr("aux");
				}
				biCamDSI_persistent = null;
				quadCLT_main.tp.resetCLTPasses();
				quadCLT_aux.tp.resetCLTPasses();

			}

			if ((quadCLT_main.correctionsParameters.rig_batch_adjust_main_gt > 0) ||
					(quadCLT_main.correctionsParameters.rig_batch_adjust_main_gt > 0)){
				{
					if (updateStatus) IJ.showStatus("Building basic DSI for the main camera image set "+quadCLT_main.image_name+ " (post rig adjustment)");
					quadCLT_main.preExpandCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
							imp_srcs_main, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
							saturation_imp_main, // boolean [][] saturation_imp, // (near) saturated pixels or null
							clt_parameters,
							debayerParameters,
							colorProcParameters,
							rgbParameters,
							threadsMax,  // maximal number of threads to launch
							updateStatus,
							debugLevelInner);
					if (updateStatus) IJ.showStatus("Expanding DSI for the main camera image set "+quadCLT_main.image_name+ " (post rig adjustment)");
					quadCLT_main.expandCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
							clt_parameters,
							debayerParameters,
							colorProcParameters,
							channelGainParameters,
							rgbParameters,
							threadsMax,  // maximal number of threads to launch
							updateStatus,
							debugLevel);
					if (updateStatus) IJ.showStatus("Refining main camera DSI with the rig for "+quadCLT_main.image_name);
					// Prepare DSI from the main camera, then refine it with the rig to get the "GT" data for the main/aux
					BiScan scan = rigInitialScan(
							quadCLT_main,  // QuadCLT            quadCLT_main,  // tiles should be set
							quadCLT_aux,  // QuadCLT            quadCLT_aux,
							clt_parameters, // EyesisCorrectionParameters.CLTParameters       clt_parameters,
							colorProcParameters,           //  ColorProcParameters                       colorProcParameters, //
							colorProcParameters_aux,          //  ColorProcParameters                       colorProcParameters, //
							threadsMax,  // final int                                      threadsMax,  // maximal number of threads to launch
							updateStatus, // final boolean                                  updateStatus,
							debugLevel); // final int                                      debugLevel)
					if (scan == null) {
						System.out.println("***** Failed to get main camera DSI and refine it with the rig, aborting this image set *****");
						continue;
					}

				}
			}

			for (int num_adjust_main = 0; num_adjust_main < quadCLT_main.correctionsParameters.rig_batch_adjust_main_gt; num_adjust_main++) {
				if (updateStatus) IJ.showStatus("Adjusting main camera image set for "+quadCLT_main.image_name+
						" (with rig DSI), pass "+(num_adjust_main+1)+" of "+quadCLT_main.correctionsParameters.rig_batch_adjust_main_gt);
				if (debugLevel > -5) {
					System.out.println("Adjusting main camera image set for "+quadCLT_main.image_name+
							" (with rig DSI), pass "+(num_adjust_main+1)+" of "+quadCLT_main.correctionsParameters.rig_batch_adjust_main_gt);
				}
				double [][] gt_disp_strength = quadCLT_main.getRigDSFromTwoQuadCL(
						this, //maybe null in no-rig mode, otherwise may contain rig measurements to be used as infinity ground truth
						clt_parameters,
						debugLevelInner); // final int        debugLevel)

//				  GeometryCorrection geometryCorrection_main = null;
//				  if (geometryCorrection.getRotMatrix(true) != null) {
//					  geometryCorrection_main = twoQuadCLT.quadCLT_main.getGeometryCorrection();
//				  }

				boolean ok = quadCLT_main.extrinsicsCLTfromGT(
						//						  this,   // TwoQuadCLT       twoQuadCLT, //maybe null in no-rig mode, otherwise may contain rig measurements to be used as infinity ground truth
						null,
						gt_disp_strength,
						clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
						false,
						threadsMax,  //final int        threadsMax,  // maximal number of threads to launch
						updateStatus,// final boolean    updateStatus,
						debugLevelInner); // final int        debugLevel)
				if (!ok) break;
			}

			for (int num_adjust_aux = 0; num_adjust_aux < quadCLT_main.correctionsParameters.rig_batch_adjust_aux_gt; num_adjust_aux++) {
				if (updateStatus) IJ.showStatus("Adjusting aux camera image set for "+quadCLT_main.image_name+
						" (with rig DSI), pass "+(num_adjust_aux+1)+" of "+quadCLT_main.correctionsParameters.rig_batch_adjust_aux_gt);
				if (debugLevel > -5) {
					System.out.println("Adjusting aux camera image set for "+quadCLT_main.image_name+
							" (with rig DSI), pass "+(num_adjust_aux+1)+" of "+quadCLT_main.correctionsParameters.rig_batch_adjust_aux_gt);
				}
				double [][] gt_disp_strength = quadCLT_aux.getRigDSFromTwoQuadCL(
						this, //maybe null in no-rig mode, otherwise may contain rig measurements to be used as infinity ground truth
						clt_parameters,
						debugLevelInner); // final int        debugLevel)

				boolean ok = quadCLT_aux.extrinsicsCLTfromGT(
//						  this,   // TwoQuadCLT       twoQuadCLT, //maybe null in no-rig mode, otherwise may contain rig measurements to be used as infinity ground truth
						  quadCLT_main.getGeometryCorrection(),
						  gt_disp_strength,
						  clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
						  false,
						  threadsMax,  //final int        threadsMax,  // maximal number of threads to launch
						  updateStatus,// final boolean    updateStatus,
						  debugLevelInner); // final int        debugLevel)
				if (!ok) break;
			}
			// Late rig adjustment, after main/aux are adjusted with rig data as ground truth
			// keeping the same DSI, required measurements will be performed anyway
			for (int num_adjust_rig = 0; num_adjust_rig < quadCLT_main.correctionsParameters.rig_batch_adjust_rig_gt; num_adjust_rig++) {
				if (updateStatus) IJ.showStatus("Adjusting dual camera rig infinity (late, after re-adjustment of the main/aux) for "+quadCLT_main.image_name+
						", pass "+(num_adjust_rig+1)+" of "+quadCLT_main.correctionsParameters.rig_batch_adjust_rig_gt);
				if (debugLevel > -5) {
					System.out.println("Adjusting dual camera rig infinity (late, after re-adjustment of the main/aux) for "+quadCLT_main.image_name+
							", pass "+(num_adjust_rig+1)+" of "+quadCLT_main.correctionsParameters.rig_batch_adjust_rig_gt);
				}

				processInfinityRig(
						quadCLT_main,               // QuadCLT                                        quadCLT_main,
						quadCLT_aux,                // QuadCLT                                        quadCLT_aux,
						imp_srcs_main,              // ImagePlus []                                   imp_quad_main,
						imp_srcs_aux,               // ImagePlus []                                   imp_quad_aux,
						saturation_imp_main,        // boolean [][]                                   saturation_main, // (near) saturated pixels or null
						saturation_imp_aux,         // boolean [][]                                   saturation_aux, // (near) saturated pixels or null
						clt_parameters,             // EyesisCorrectionParameters.CLTParameters       clt_parameters,
						threadsMax,                 // final int        threadsMax,  // maximal number of threads to launch
						updateStatus,               // final boolean    updateStatus,
						debugLevelInner-2);                // final int        debugLevel);
			}
			// Delete previous results if any of the late adjustments were performed
			if (    (quadCLT_main.correctionsParameters.rig_batch_adjust_main_gt > 0) ||
					(quadCLT_main.correctionsParameters.rig_batch_adjust_aux_gt > 0) ||
					(quadCLT_main.correctionsParameters.rig_batch_adjust_rig_gt > 0)) {

				if (debugLevel > -5) {
					System.out.println("\n---- Field adjustments after the second round ----");
					quadCLT_main.showExtrinsicCorr("main");  // show_fine_corr("main");
					quadCLT_aux.showExtrinsicCorr("aux");    // show_fine_corr("aux");
					quadCLT_aux.geometryCorrection.showRig();// show_fine_corr("aux");
				}
				biCamDSI_persistent = null;
				quadCLT_main.tp.resetCLTPasses();
				quadCLT_aux.tp.resetCLTPasses();

			}

			// Generate 8 images and thumbnails
			if (quadCLT_main.correctionsParameters.clt_batch_4img){
				if (updateStatus) IJ.showStatus("Rendering 8 image set (disparity = 0) for "+quadCLT_main.image_name);
				processCLTQuadCorrPair(
						quadCLT_main,               // QuadCLT                                        quadCLT_main,
						quadCLT_aux,                // QuadCLT                                        quadCLT_aux,
						imp_srcs_main,              // ImagePlus []                                   imp_quad_main,
						imp_srcs_aux,               // ImagePlus []                                   imp_quad_aux,
						saturation_imp_main,        // boolean [][]                                   saturation_main, // (near) saturated pixels or null
						saturation_imp_aux,         // boolean [][]                                   saturation_aux, // (near) saturated pixels or null
						clt_parameters,             // EyesisCorrectionParameters.CLTParameters       clt_parameters,
						debayerParameters,          // EyesisCorrectionParameters.DebayerParameters   debayerParameters,
						colorProcParameters,        // EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
						colorProcParameters_aux,    // EyesisCorrectionParameters.ColorProcParameters colorProcParameters_aux,
						rgbParameters,              // EyesisCorrectionParameters.RGBParameters       rgbParameters,
						scaleExposures_main,        // double []	                                     scaleExposures_main, // probably not needed here - restores brightness of the final image
						scaleExposures_aux,         // double []	                                     scaleExposures_aux, // probably not needed here - restores brightness of the final image
						false,                      //  final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
						// averages measurements
						clt_parameters.rig.lt_avg_radius,// final int                                      lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
						threadsMax,                 // final int        threadsMax,  // maximal number of threads to launch
						updateStatus,               // final boolean    updateStatus,
						debugLevelInner);                // final int        debugLevel);
			}
			// ground truth and then all the rest

			quadCLT_main.tp.resetCLTPasses();
			quadCLT_aux.tp.resetCLTPasses();
			if (quadCLT_main.correctionsParameters.clt_batch_dsi_aux) {
				if (updateStatus) IJ.showStatus("Building basic DSI for the aux camera image set "+quadCLT_main.image_name+" (for DSI export)");
				quadCLT_aux.preExpandCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
						imp_srcs_aux, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
						saturation_imp_aux, // boolean [][] saturation_imp, // (near) saturated pixels or null
						clt_parameters,
						debayerParameters,
						colorProcParameters_aux,
						rgbParameters,
						threadsMax,  // maximal number of threads to launch
						updateStatus,
						debugLevelInner);
				if (updateStatus) IJ.showStatus("Expanding DSI for the aux camera image set "+quadCLT_main.image_name+" (for DSI export)");
				quadCLT_aux.expandCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
						clt_parameters,
						debayerParameters,
						colorProcParameters_aux,
						channelGainParameters,
						rgbParameters,
						threadsMax,  // maximal number of threads to launch
						updateStatus,
						debugLevel);
				double [][] aux_last_scan = quadCLT_aux.tp.getShowDS(
						quadCLT_aux.tp.clt_3d_passes.get( quadCLT_aux.tp.clt_3d_passes.size() -1),
						false); // boolean force_final);

				dsi[DSI_DISPARITY_AUX] = aux_last_scan[0];
				dsi[DSI_STRENGTH_AUX] =  aux_last_scan[1];
				quadCLT_aux.tp.resetCLTPasses();
			}

			if (quadCLT_main.correctionsParameters.clt_batch_explore) {
				if (updateStatus) IJ.showStatus("Building basic DSI for the main camera image set "+quadCLT_main.image_name+" (after all adjustments)");
				quadCLT_main.preExpandCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
						imp_srcs_main, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
						saturation_imp_main, // boolean [][] saturation_imp, // (near) saturated pixels or null
						clt_parameters,
						debayerParameters,
						colorProcParameters,
						rgbParameters,
						threadsMax,  // maximal number of threads to launch
						updateStatus,
						debugLevelInner);
				if (updateStatus) IJ.showStatus("Expanding DSI for the main camera image set "+quadCLT_main.image_name+" (after all adjustments)");
				quadCLT_main.expandCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
						clt_parameters,
						debayerParameters,
						colorProcParameters,
						channelGainParameters,
						rgbParameters,
						threadsMax,  // maximal number of threads to launch
						updateStatus,
						debugLevel);
				double [][] main_last_scan = quadCLT_main.tp.getShowDS(
						quadCLT_main.tp.clt_3d_passes.get( quadCLT_main.tp.clt_3d_passes.size() -1),
						false); // boolean force_final);

				dsi[DSI_DISPARITY_MAIN] = main_last_scan[0];
				dsi[DSI_STRENGTH_MAIN] =  main_last_scan[1];

//				CLTPass3d scan_last = quadCLT_main.tp.clt_3d_passes.get( quadCLT_main.tp.clt_3d_passes.size() -1);

				if (updateStatus) IJ.showStatus("Creating \"ground truth\" DSI using dual-camera rig "+quadCLT_main.image_name);
				groundTruth(                        // actually there is no sense to process multiple image sets. Combine with other processing?
						quadCLT_main,               // QuadCLT quadCLT_main,
						quadCLT_aux,                // QuadCLT quadCLT_aux,
						clt_parameters,             // EyesisCorrectionParameters.CLTParameters       clt_parameters,
						colorProcParameters,          //  ColorProcParameters                       colorProcParameters, //
						colorProcParameters_aux,          //  ColorProcParameters                       colorProcParameters, //
						threadsMax,                 // final int        threadsMax,  // maximal number of threads to launch
						updateStatus,               // final boolean    updateStatus,
						debugLevelInner-2);                // final int        debugLevel);
				double [][] rig_last_scan = quadCLT_main.tp.getShowDS(
						quadCLT_main.tp.clt_3d_passes.get( quadCLT_main.tp.clt_3d_passes.size() -1),
						false); // boolean force_final);

				dsi[DSI_DISPARITY_RIG] = rig_last_scan[0];
				dsi[DSI_STRENGTH_RIG] =  rig_last_scan[1];

			} else {
				continue;
			}
			resetRig   (false); // includes biCamDSI_persistent = null
			Runtime.getRuntime().gc();
			System.out.println("--- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");

			if (quadCLT_main.correctionsParameters.clt_batch_surf) {
				if (updateStatus) IJ.showStatus("Creating and filtering supertile plane surfaces from the DSI "+quadCLT_main.image_name);

				quadCLT_main.tp.showPlanes(
						clt_parameters,
						quadCLT_main.geometryCorrection,
						threadsMax,
						updateStatus,
						debugLevelInner);

			} else continue; // if (correctionsParameters.clt_batch_surf)

			if (quadCLT_main.correctionsParameters.clt_batch_assign) {
				if (updateStatus) IJ.showStatus("Assigning tiles to candidate surfaces "+quadCLT_main.image_name);
				// prepare average RGBA for the last scan
				quadCLT_main.setPassAvgRBGA(                      // get image from a single pass, return relative path for x3d // USED in lwir
						clt_parameters,                           // CLTParameters           clt_parameters,
						quadCLT_main.tp.clt_3d_passes.size() - 1, // int        scanIndex,
						threadsMax,                               // int        threadsMax,  // maximal number of threads to launch
						updateStatus,                             // boolean    updateStatus,
						debugLevelInner);                         // int        debugLevel)
				double [][] assignments_dbg = quadCLT_main.tp.assignTilesToSurfaces(
						clt_parameters,
						quadCLT_main.geometryCorrection,
						threadsMax,
						updateStatus,
						debugLevelInner);
				if (assignments_dbg == null) continue;
				dsi[DSI_DISPARITY_X3D] = assignments_dbg[TileSurface.ASGN_A_DISP];

// TODO use assignments_dbg

			} else continue; // if (correctionsParameters.clt_batch_assign)
			// generate ML data if enabled
			if (quadCLT_main.correctionsParameters.clt_batch_genMl) { // rig.ml_generate) { //clt_batch_genMl
				outputMLData(
						quadCLT_main,   // QuadCLT                                  quadCLT_main,  // tiles should be set
						quadCLT_aux,    // QuadCLT                                  quadCLT_aux,
						clt_parameters, // EyesisCorrectionParameters.CLTParameters clt_parameters,
	    				null,           //String                                   ml_directory,       // full path or null (will use config one)
						threadsMax,     // final int                                threadsMax,  // maximal number of threads to launch
						updateStatus,   // final boolean                            updateStatus,
						debugLevel);    // final int                                debugLevel)
				/*
				if (clt_parameters.rig.ml_copyJP4) {
					copyJP4src(
							quadCLT_main,   // QuadCLT                                  quadCLT_main,  // tiles should be set
							quadCLT_aux,    // QuadCLT                                  quadCLT_aux,
							clt_parameters, // EyesisCorrectionParameters.CLTParameters clt_parameters,
							debugLevel);    // final int                                debugLevel)
				}*/
			}
			// copy regardless of ML generation
			if (clt_parameters.rig.ml_copyJP4) {
				copyJP4src(
						quadCLT_main,   // QuadCLT                                  quadCLT_main,  // tiles should be set
						quadCLT_aux,    // QuadCLT                                  quadCLT_aux,
						clt_parameters, // EyesisCorrectionParameters.CLTParameters clt_parameters,
						debugLevel);    // final int                                debugLevel)
			}

			if (quadCLT_main.correctionsParameters.clt_batch_gen3d) {
				if (updateStatus) IJ.showStatus("Generating and exporting 3D scene model "+quadCLT_main.image_name);
				boolean ok = quadCLT_main.output3d(
						clt_parameters,      // EyesisCorrectionParameters.CLTParameters           clt_parameters,
						colorProcParameters, // EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
						rgbParameters,       // EyesisCorrectionParameters.RGBParameters             rgbParameters,
						threadsMax,          // final int        threadsMax,  // maximal number of threads to launch
						updateStatus,        // final boolean    updateStatus,
						debugLevelInner);         // final int        debugLevel)
				if (!ok) continue;
			} else continue; // if (correctionsParameters.clt_batch_gen3d)


			if (quadCLT_main.correctionsParameters.clt_batch_dsi) {
				quadCLT_main.saveDSI (
						dsi);
			}

			if (quadCLT_main.correctionsParameters.clt_batch_save_extrinsics) {
				 saveProperties(
							null,        // String path,                // full name with extension or w/o path to use x3d directory
							null,        // Properties properties,    // if null - will only save extrinsics)
							debugLevel);
			}
			if (quadCLT_main.correctionsParameters.clt_batch_save_all) {
				 saveProperties(
							null,        // String path,                // full name with extension or w/o path to use x3d directory
							properties,  // Properties properties,    // if null - will only save extrinsics)
							debugLevel);
			}
			Runtime.getRuntime().gc();
			if (debugLevel >-1) System.out.println("Processing set "+(nSet+1)+" (of "+set_channels_aux.length+") finished at "+
					IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");

			if (quadCLT_aux.eyesisCorrections.stopRequested.get()>0) {
				System.out.println("User requested stop");
				System.out.println("Processing "+(nSet + 1)+" file sets (of "+set_channels_main.length+") finished at "+
						IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
				return;
			}
		}
		System.out.println("batchRig(): processing "+(quadCLT_main.getTotalFiles(set_channels_main)+quadCLT_aux.getTotalFiles(set_channels_aux))+" files ("+set_channels_main.length+" file sets) finished at "+
				IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");

	}

	public void TestInterScene(
			QuadCLT                                              quadCLT_main, // tiles should be set
//			QuadCLT                                              quadCLT_aux,
			CLTParameters             clt_parameters,
			EyesisCorrectionParameters.DebayerParameters         debayerParameters,
			ColorProcParameters                                  colorProcParameters,
			ColorProcParameters                                  colorProcParameters_aux,
			CorrectionColorProc.ColorGainsParameters             channelGainParameters,
			EyesisCorrectionParameters.RGBParameters             rgbParameters,
			EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
			Properties                                           properties,
			final int        threadsMax,  // maximal number of threads to launch
			final boolean    updateStatus,
			final int        debugLevel)  throws Exception
	{
		if ((quadCLT_main != null) && (quadCLT_main.getGPU() != null)) {
			quadCLT_main.getGPU().resetGeometryCorrection();
			quadCLT_main.gpuResetCorrVector(); // .getGPU().resetGeometryCorrectionVector();
		}
		// final boolean    batch_mode = clt_parameters.batch_run;
		this.startTime=System.nanoTime();
		String [] sourceFiles0=quadCLT_main.correctionsParameters.getSourcePaths();
		QuadCLT.SetChannels [] set_channels_main = quadCLT_main.setChannels(debugLevel);
		if ((set_channels_main == null) || (set_channels_main.length==0)) {
			System.out.println("No files to process (of "+sourceFiles0.length+")");
			return;
		}
		QuadCLT.SetChannels [] set_channels=quadCLT_main.setChannels(debugLevel);
//		String set_name = set_channels[0].set_name;
		
		QuadCLT [] quadCLTs = new QuadCLT [set_channels.length]; 
		for (int i = 0; i < quadCLTs.length; i++) {
			quadCLTs[i] = quadCLT_main.spawnQuadCLT(
					set_channels[i].set_name,
					clt_parameters,
					colorProcParameters, //
					threadsMax,
					debugLevel);
			// temporarily fix wrong sign:
			ErsCorrection ers = (ErsCorrection) (quadCLTs[i].getGeometryCorrection());
			ers.setupERSfromExtrinsics();			
			quadCLTs[i].setDSRBG(
					clt_parameters, // CLTParameters  clt_parameters,
					threadsMax,     // int            threadsMax,  // maximal number of threads to launch
					updateStatus,   // boolean        updateStatus,
					debugLevel);    // int            debugLevel)
///			quadCLTs[i].showDSIMain();
		}
		
		
		double k_prev = 0.75;
		
		double corr_scale = 0.75;
///		for (int i = 0; i < quadCLTs.length; i++) {
		OpticalFlow opticalFlow = new OpticalFlow(
				threadsMax, // int            threadsMax,  // maximal number of threads to launch
				updateStatus); // boolean        updateStatus);
//		int         margin = 4; //  extra margins over 16x16 tiles to accommodate distorted destination tiles
		/*
		
		double      tolerance_absolute = 0.25; // absolute disparity half-range in each tile
		double      tolerance_relative = 0.2; // relative disparity half-range in each tile
		double      center_occupancy =   0.25;   // fraction of remaining  tiles in the center 8x8 area (<1.0)

		int         num_passes = 100;
		double      max_change = 0.005 ;
		*/
		
		for (int i = 1; i < quadCLTs.length; i++) {
			QuadCLT qPrev = (i > 0) ? quadCLTs[i - 1] : null;
/*			
			double [][][] source_tiles = opticalFlow.prepareReferenceTiles(
					quadCLTs[i],        // final QuadCLT     qthis,
//					margin,             // final int         margin, // extra margins over 16x16 tiles to accommodate distorted destination tiles
					tolerance_absolute, // final double      tolerance_absolute, // absolute disparity half-range in each tile
					tolerance_relative, // final double      tolerance_relative, // relative disparity half-range in each tile
					center_occupancy,   // final double      center_occupancy,   // fraction of remaining  tiles in the center 8x8 area (<1.0)
					2); // final int         debug_level)
			
			opticalFlow.fillTilesNans(
					source_tiles,          // final double [][][] nan_tiles,
					quadCLTs[i],           // final QuadCLT     qthis,
//					margin,                // final int         margin, // extra margins over 16x16 tiles to accommodate distorted destination tiles
//					0.5 * Math.sqrt(2.0),  // double    diagonal_weight, // relative to ortho
					num_passes,            // final int         num_passes,
					max_change,            // final double      max_change,
					2);                    // final int         debug_level)
			*/
			

			double [][][] pair_sets =
					opticalFlow.get_pair(
					k_prev,
					quadCLTs[i],
					qPrev,
					corr_scale,
					1); // -1); // int debug_level);
		}

		/*
		for (int i = 0; i < quadCLTs.length; i++) {
			quadCLTs[i].saveInterProperties(// save properties for interscene processing (extrinsics, ers, ...)
					null, // String path,             // full name with extension or w/o path to use x3d directory
					debugLevel);
		}		
		 */

		System.out.println("End of test");


	}
/*
 * 		    public double [][]  getDSRBG (){
		    	return dsrbg;
		    }


 */
	
	public void batchLwirRig(
			QuadCLT                                              quadCLT_main, // tiles should be set
			QuadCLT                                              quadCLT_aux,
			CLTParameters             clt_parameters,
			EyesisCorrectionParameters.DebayerParameters         debayerParameters,
			ColorProcParameters                                  colorProcParameters,
			ColorProcParameters                                  colorProcParameters_aux,
			CorrectionColorProc.ColorGainsParameters             channelGainParameters,
			EyesisCorrectionParameters.RGBParameters             rgbParameters,
			EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
			Properties                                           properties,
			final int        threadsMax,  // maximal number of threads to launch
			final boolean    updateStatus,
			final int        debugLevel)  throws Exception
	{
		if ((quadCLT_main != null) && (quadCLT_main.getGPU() != null)) {
			quadCLT_main.getGPU().resetGeometryCorrection();
			quadCLT_main.gpuResetCorrVector(); // .getGPU().resetGeometryCorrectionVector();
		}
		if ((quadCLT_aux != null) && (quadCLT_aux.getGPU() != null)) {
			quadCLT_aux.getGPU().resetGeometryCorrection();
			quadCLT_aux.gpuResetCorrVector(); // .getGPU().resetGeometryCorrectionVector();
		}

		// final boolean    batch_mode = clt_parameters.batch_run;
		// Reset dsi data (only 2 slices will be used)
		this.dsi =           new double [DSI_SLICES.length][];
		this.dsi_aux_from_main =   null; // full data, including rms, fg and bg data
		quadCLT_aux.ds_from_main = null; // this is for adjust6ment only (short)

		final int        debugLevelInner=clt_parameters.batch_run? -2: debugLevel;
		this.startTime=System.nanoTime();
		String [] sourceFiles=quadCLT_main.correctionsParameters.getSourcePaths();
		QuadCLT.SetChannels [] set_channels_main = quadCLT_main.setChannels(debugLevel);
		QuadCLT.SetChannels [] set_channels_aux =  quadCLT_aux.setChannels(debugLevel);
		if ((set_channels_main == null) || (set_channels_main.length==0) || (set_channels_aux == null) || (set_channels_aux.length==0)) {
			System.out.println("No files to process (of "+sourceFiles.length+")");
			return;
		}
		double [] referenceExposures_main = null;
		double [] referenceExposures_aux =  null;
		if (!colorProcParameters.lwir_islwir)      referenceExposures_main = quadCLT_main.eyesisCorrections.calcReferenceExposures(debugLevel);
		if (!colorProcParameters_aux.lwir_islwir)  referenceExposures_aux =  quadCLT_aux.eyesisCorrections.calcReferenceExposures(debugLevel);
		for (int nSet = 0; nSet < set_channels_main.length; nSet++){
			// check it is the same set for both cameras
			if (set_channels_aux.length <= nSet ) {
				throw new Exception ("Set names for cameras do not match: main camera: '"+set_channels_main[nSet].name()+"', aux. camera: nothing");
			}
			if (!set_channels_main[nSet].name().equals(set_channels_aux[nSet].name())) {
				throw new Exception ("Set names for cameras do not match: main camera: '"+set_channels_main[nSet].name()+"', aux. camera: '"+set_channels_main[nSet].name()+"'");
			}

			int [] channelFiles_main = set_channels_main[nSet].fileNumber();
			int [] channelFiles_aux =  set_channels_aux[nSet].fileNumber();
			boolean [][] saturation_imp_main = (clt_parameters.sat_level > 0.0)? new boolean[channelFiles_main.length][] : null;
			boolean [][] saturation_imp_aux =  (clt_parameters.sat_level > 0.0)? new boolean[channelFiles_main.length][] : null;
			double [] scaleExposures_main = new double[channelFiles_main.length];
			double [] scaleExposures_aux =  new double[channelFiles_main.length];
			if (updateStatus) IJ.showStatus("Conditioning main camera image set for "+quadCLT_main.image_name);
			ImagePlus [] imp_srcs_main = quadCLT_main.conditionImageSet(
					clt_parameters,                 // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					colorProcParameters,            //  ColorProcParameters                       colorProcParameters, //
					sourceFiles,                    // String []                                 sourceFiles,
					set_channels_main[nSet].name(), // String                                    set_name,
					referenceExposures_main,        // double []                                 referenceExposures,
					channelFiles_main,              // int []                                    channelFiles,
					scaleExposures_main,            //output  // double [] scaleExposures
					saturation_imp_main,            //output  // boolean [][]                              saturation_imp,
					threadsMax,                 // int                                       threadsMax,
					debugLevelInner); // int                                       debugLevel);
			if (updateStatus) IJ.showStatus("Conditioning aux camera image set for "+quadCLT_main.image_name);

			// optionally adjust main, aux (aux always will use main - calculate if needed
			// Early main camera adjustment, rig data is not available
			// with LWIR only 1 type of adjustments is possibkle - pre for main, post for aux. Combine configuration fields made for the 8-rig
			int adjust_main = (quadCLT_main.correctionsParameters.rig_batch_adjust_main > quadCLT_main.correctionsParameters.rig_batch_adjust_main_gt)?
					quadCLT_main.correctionsParameters.rig_batch_adjust_main : quadCLT_main.correctionsParameters.rig_batch_adjust_main_gt;
			int adjust_aux = (quadCLT_main.correctionsParameters.rig_batch_adjust_aux > quadCLT_main.correctionsParameters.rig_batch_adjust_aux_gt)?
					quadCLT_main.correctionsParameters.rig_batch_adjust_aux : quadCLT_main.correctionsParameters.rig_batch_adjust_aux_gt;

			for (int num_adjust_main = 0; num_adjust_main < adjust_main; num_adjust_main++) {
				if (updateStatus) IJ.showStatus("Building basic  DSI for the main camera image set "+quadCLT_main.image_name+
						", pass "+(num_adjust_main+1)+" of "+adjust_main);
				if (debugLevel > -5) {
					System.out.println("Building basic  DSI for the main camera image set "+quadCLT_main.image_name+
							", pass "+(num_adjust_main+1)+" of "+adjust_main);
				}
//Generates background image in model tree - should be done later, after adjustment (It is overwritten later, so OK)
				quadCLT_main.preExpandCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
						imp_srcs_main, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
						saturation_imp_main, // boolean [][] saturation_imp, // (near) saturated pixels or null
						clt_parameters,
						debayerParameters,
						colorProcParameters,
						rgbParameters,
						threadsMax,  // maximal number of threads to launch
						updateStatus,
						debugLevelInner);
				// adjust extrinsics here
				System.out.println("Adjust main extrinsics here");
				if (updateStatus) IJ.showStatus("Adjusting main camera image set for "+quadCLT_main.image_name+
						", pass "+(num_adjust_main+1)+" of "+adjust_main);
				if (debugLevel > -5) {
					System.out.println("Adjusting main camera image set for "+quadCLT_main.image_name+
							", pass "+(num_adjust_main+1)+" of "+adjust_main);
				}
				if (debugLevel > -1){
					int scan_index =  quadCLT_main.tp.clt_3d_passes.size() -1;
					quadCLT_main.tp.showScan(
							quadCLT_main.tp.clt_3d_passes.get(scan_index),   // CLTPass3d   scan,
							"pre-adjust-extrinsic-scan-"+scan_index); //String title)
					for (int s = 0; (s < 5) && (s < scan_index); s++) {
						quadCLT_main.tp.showScan(
								quadCLT_main.tp.clt_3d_passes.get(s),   // CLTPass3d   scan,
								"pre-adjust-extrinsic-scan-"+s); //String title)
					}
				}
			    double inf_min = clt_parameters.ly_inf_min_broad; // -0.5;
			    double inf_max = clt_parameters.ly_inf_max_broad; // 0.5;
			    if (clt_parameters.ly_inf_force_fine || (num_adjust_main >= (adjust_main/2))) {
			        inf_min = clt_parameters.ly_inf_min_narrow; // -0.2;
			        inf_max = clt_parameters.ly_inf_max_narrow; // 0.05;
			        System.out.println("Late adjustment, using narrow band infinity detection, inf_min="+inf_min+", inf_max="+inf_max);
			    } else {
			        System.out.println("Early adjustment, using wide band infinity detection, inf_min="+inf_min+", inf_max="+inf_max);
			    }
				boolean ok = quadCLT_main.extrinsicsCLT(
						clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
						false, // adjust_poly,
						inf_min, // double inf_min,
						inf_max,  // double inf_max,
						threadsMax,  //final int        threadsMax,  // maximal number of threads to launch
						updateStatus,// final boolean    updateStatus,
						debugLevelInner); // final int        debugLevel)
// clear memory for main
				quadCLT_main.tp.resetCLTPasses();
				if (!ok) break;
			}
			if (quadCLT_main.correctionsParameters.clt_batch_dsi1){
		        System.out.println("Trying experimental features DSI/ERS");
				quadCLT_main.preExpandCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
						imp_srcs_main, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
						saturation_imp_main, // boolean [][] saturation_imp, // (near) saturated pixels or null
						clt_parameters,
						debayerParameters,
						colorProcParameters,
						rgbParameters,
						threadsMax,  // maximal number of threads to launch
						updateStatus,
						debugLevelInner);
				
				double [][] dsi_ly = quadCLT_main.filterByLY(
						clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
				        clt_parameters.ly_inf_min_narrow, // double inf_min,
				        clt_parameters.ly_inf_max_narrow,  // double inf_max,
						threadsMax,  //final int        threadsMax,  // maximal number of threads to launch
						updateStatus,// final boolean    updateStatus,
						debugLevelInner); // final int        debugLevel)				
				dsi[DSI_DISPARITY_MAIN] = dsi_ly[0];
				dsi[DSI_STRENGTH_MAIN] =  dsi_ly[1];
//				if (quadCLT_main.correctionsParameters.clt_batch_dsi) { // Should be always enabled ?
				quadCLT_main.saveDSIMain (dsi);
//				}
				// clear memory for main
				quadCLT_main.tp.resetCLTPasses();
				// copy regardless of ML generation
				// See if it will copy all files, not just the main camera ones

				if (clt_parameters.rig.ml_copyJP4) {
					copyJP4src(
							quadCLT_main,   // QuadCLT                                  quadCLT_main,  // tiles should be set
							quadCLT_aux,    // QuadCLT                                  quadCLT_aux,
							clt_parameters, // EyesisCorrectionParameters.CLTParameters clt_parameters,
							debugLevel);    // final int                                debugLevel)
				}
			}			
			
			
			// Generate 4 main camera images and thumbnail
			if (quadCLT_main.correctionsParameters.clt_batch_4img){
				if (clt_parameters.gpu_use_main) {
					if (updateStatus) IJ.showStatus("GPU: Rendering 4 image set (disparity = 0) for "+quadCLT_main.image_name+ "and a thumb nail");
					quadCLT_main.processCLTQuadCorrGPU(
							imp_srcs_main,       // ImagePlus []                                    imp_quad,
							saturation_imp_main, // boolean [][] saturation_imp, // (near) saturated pixels or null
							clt_parameters,      // CLTParameters                                   clt_parameters,
							debayerParameters,   // EyesisCorrectionParameters.DebayerParameters    debayerParameters,
							colorProcParameters, // ColorProcParameters                             colorProcParameters,
							channelGainParameters,
							rgbParameters,       // EyesisCorrectionParameters.RGBParameters        rgbParameters,
							scaleExposures_main, // double []	                                    scaleExposures, // probably not needed here - restores brightness of the final image
							threadsMax,          // final int        threadsMax,  // maximal number of threads to launch
							updateStatus,        // final boolean    updateStatus,
							debugLevel);        // final int        debugLevel);
				} else {
					if (updateStatus) IJ.showStatus("CPU: Rendering 4 image set (disparity = 0) for "+quadCLT_main.image_name+ "and a thumb nail");

					quadCLT_main.processCLTQuadCorrCPU( // returns ImagePlus, but it already should be saved/shown
							imp_srcs_main, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
							saturation_imp_main, // boolean [][] saturation_imp, // (near) saturated pixels or null
							clt_parameters,
							debayerParameters,
							colorProcParameters,
							channelGainParameters,
							rgbParameters,
							scaleExposures_main,
							false, // calculate and apply additional fine geometry correction
							false, // calculate and apply geometry correction at infinity
							threadsMax,  // maximal number of threads to launch
							updateStatus,
							debugLevel);
				}
				quadCLT_main.tp.resetCLTPasses();
			}


			if (quadCLT_main.correctionsParameters.clt_batch_explore) {
				if (updateStatus) IJ.showStatus("Building basic DSI for the main camera image set "+quadCLT_main.image_name+" (after all adjustments)");
				quadCLT_main.preExpandCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
						imp_srcs_main, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
						saturation_imp_main, // boolean [][] saturation_imp, // (near) saturated pixels or null
						clt_parameters,
						debayerParameters,
						colorProcParameters,
						rgbParameters,
						threadsMax,  // maximal number of threads to launch
						updateStatus,
						debugLevelInner);
				if (updateStatus) IJ.showStatus("Expanding DSI for the main camera image set "+quadCLT_main.image_name+" (after all adjustments)");
				quadCLT_main.expandCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
						clt_parameters,
						debayerParameters,
						colorProcParameters,
						channelGainParameters,
						rgbParameters,
						threadsMax,  // maximal number of threads to launch
						updateStatus,
						debugLevel);
				double [][] main_last_scan = quadCLT_main.tp.getShowDS(
						quadCLT_main.tp.clt_3d_passes.get( quadCLT_main.tp.clt_3d_passes.size() -1),
						false); // boolean force_final);

				  if (debugLevel > -1) { //-5){
					  int scan_index =  quadCLT_main.tp.clt_3d_passes.size() -1;
					  quadCLT_main.tp.showScan(
							  quadCLT_main.tp.clt_3d_passes.get(scan_index),   // CLTPass3d   scan,
							  "test_pre-after-"+scan_index); //String title)
				  }

				
				
				dsi[DSI_DISPARITY_MAIN] = main_last_scan[0];
				dsi[DSI_STRENGTH_MAIN] =  main_last_scan[1];
				if (quadCLT_main.correctionsParameters.clt_batch_dsi) { // Should be always enabled ?
					quadCLT_main.saveDSIMain (
							dsi);
				}


				Runtime.getRuntime().gc();
				System.out.println("--- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");

				if (quadCLT_main.correctionsParameters.clt_batch_surf) {
					if (updateStatus) IJ.showStatus("Creating and filtering supertile plane surfaces from the DSI "+quadCLT_main.image_name);
					quadCLT_main.tp.showPlanes(
							clt_parameters,
							quadCLT_main.geometryCorrection,
							threadsMax,
							updateStatus,
							debugLevelInner);

					if (quadCLT_main.correctionsParameters.clt_batch_assign) {
						if (updateStatus) IJ.showStatus("Assigning tiles to candidate surfaces "+quadCLT_main.image_name);
						// prepare average RGBA for the last scan
						quadCLT_main.setPassAvgRBGA(                      // get image from a single pass, return relative path for x3d // USED in lwir
								clt_parameters,                           // CLTParameters           clt_parameters,
								quadCLT_main.tp.clt_3d_passes.size() - 1, // int        scanIndex,
								threadsMax,                               // int        threadsMax,  // maximal number of threads to launch
								updateStatus,                             // boolean    updateStatus,
								debugLevelInner);                         // int        debugLevel)
						double [][] assignments_dbg = quadCLT_main.tp.assignTilesToSurfaces(
								clt_parameters,
								quadCLT_main.geometryCorrection,
								threadsMax,
								updateStatus,
								debugLevelInner);
						if (assignments_dbg == null) continue;
						dsi[DSI_DISPARITY_X3D] = assignments_dbg[TileSurface.ASGN_A_DISP];

						// TODO use assignments_dbg

						// generate ML data if enabled
						/*
						if (quadCLT_main.correctionsParameters.clt_batch_genMl) { // rig.ml_generate) { //clt_batch_genMl
							outputMLData(
									quadCLT_main,   // QuadCLT                                  quadCLT_main,  // tiles should be set
									quadCLT_aux,    // QuadCLT                                  quadCLT_aux,
									clt_parameters, // EyesisCorrectionParameters.CLTParameters clt_parameters,
									null,           //String                                   ml_directory,       // full path or null (will use config one)
									threadsMax,     // final int                                threadsMax,  // maximal number of threads to launch
									updateStatus,   // final boolean                            updateStatus,
									debugLevel);    // final int                                debugLevel)
						}
						*/
						// copy regardless of ML generation
						// See if it will copy all files, not just the main camera ones

						if (clt_parameters.rig.ml_copyJP4) {
							copyJP4src(
									quadCLT_main,   // QuadCLT                                  quadCLT_main,  // tiles should be set
									quadCLT_aux,    // QuadCLT                                  quadCLT_aux,
									clt_parameters, // EyesisCorrectionParameters.CLTParameters clt_parameters,
									debugLevel);    // final int                                debugLevel)
						}

						if (quadCLT_main.correctionsParameters.clt_batch_gen3d) {
							if (updateStatus) IJ.showStatus("Generating and exporting 3D scene model "+quadCLT_main.image_name);
							boolean ok = quadCLT_main.output3d(
									clt_parameters,      // EyesisCorrectionParameters.CLTParameters           clt_parameters,
									colorProcParameters, // EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
									rgbParameters,       // EyesisCorrectionParameters.RGBParameters             rgbParameters,
									threadsMax,          // final int        threadsMax,  // maximal number of threads to launch
									updateStatus,        // final boolean    updateStatus,
									debugLevelInner);         // final int        debugLevel)
							if (!ok) continue;
						}
						// save assigned disparity also? - with "-DSI_COMBO" suffix
						if (quadCLT_main.correctionsParameters.clt_batch_dsi) {
							quadCLT_main.saveDSI (
									dsi);
						}

					}
				}
			} else { // if (quadCLT_main.correctionsParameters.clt_batch_explore) {
				int num_restored =  0;
				try {
					num_restored =  quadCLT_main.restoreDSI(DSI_MAIN_SUFFIX, // "-DSI_COMBO", "-DSI_MAIN"
							dsi);
							
				} catch (Exception e) {

				}
				if (num_restored < 2) {
					System.out.println("No DSI from the main camera is available. Please re-run with 'clt_batch_explore' enabled to generate it");
					if (quadCLT_main.correctionsParameters.clt_batch_save_extrinsics) {
						saveProperties(
								null,        // String path,                // full name with extension or w/o path to use x3d directory
								null,        // Properties properties,      // if null - will only save extrinsics)
								debugLevel);
						
						quadCLT_main.saveInterProperties( // save properties for interscene processing (extrinsics, ers, ...)
								null, // String path,             // full name with extension or w/o path to use x3d directory
//								null, // Properties properties,   // if null - will only save extrinsics)
								debugLevel);
					}
					
					
					
					if (quadCLT_main.correctionsParameters.clt_batch_save_all) {
						saveProperties(
								null,        // String path,                // full name with extension or w/o path to use x3d directory
								properties,  // Properties properties,    // if null - will only save extrinsics)
								debugLevel);
					}
					continue; // skipping to the next file
				}
			}

			// Process AUX (LWIR) camera data
			// 1) Prepare DS for adjustments (just d/s, with ambiguous disparity tiles removed)
			// 2) Prepare full D/S and FG/BG data to be embedded within the ML files
			double [][] main_ds = {dsi[DSI_DISPARITY_MAIN], dsi[DSI_STRENGTH_MAIN]};

			if ((adjust_aux == 0) &&
					!quadCLT_main.correctionsParameters.clt_batch_4img_aux &&
					!quadCLT_main.correctionsParameters.clt_batch_dsi_aux &&
					!quadCLT_main.correctionsParameters.clt_batch_genMl &&
					!quadCLT_main.correctionsParameters.clt_batch_save_extrinsics &&
					!quadCLT_main.correctionsParameters.clt_batch_save_all) {
				continue; 
			}
			
			
			quadCLT_aux.ds_from_main = quadCLT_aux.depthMapMainToAux( // only 2 layers for adjustments
					main_ds, // double [][] ds,
					quadCLT_main.getGeometryCorrection(), //  GeometryCorrection geometryCorrection_main,
					quadCLT_aux.getGeometryCorrection(), //  GeometryCorrection geometryCorrection_aux,
					clt_parameters,
					false,             // split_fg_bg,
					true,              // for_adjust,
					debugLevel);       // DEBUG_LEVEL); // int debug_level

			this.dsi_aux_from_main = quadCLT_aux.depthMapMainToAux( // 8 layers for ML generation/exporting + 2 zero layers
					main_ds, // double [][] ds,
					quadCLT_main.getGeometryCorrection(), //  GeometryCorrection geometryCorrection_main,
					quadCLT_aux.getGeometryCorrection(),  //  GeometryCorrection geometryCorrection_aux,
					clt_parameters,
					true,                                 // split_fg_bg,
					false,                                // for_adjust,
					debugLevel);                          // int debug_level

			ImagePlus [] imp_srcs_aux = quadCLT_aux.conditionImageSet(
					clt_parameters,                 // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					colorProcParameters_aux,        //  ColorProcParameters                       colorProcParameters, //
					sourceFiles,                    // String []                                 sourceFiles,
					set_channels_aux[nSet].name(), // String                                    set_name,
					referenceExposures_aux,        // double []                                 referenceExposures,
					channelFiles_aux,              // int []                                    channelFiles,
					scaleExposures_aux,            //output  // double [] scaleExposures
					saturation_imp_aux,            //output  // boolean [][]                              saturation_imp,
					threadsMax,                 // int                                       threadsMax,
					debugLevelInner); // int                                       debugLevel);

			// optionally adjust AUX extrinsics (using quadCLT_aux.ds_from_main )
			for (int num_adjust_aux = 0; num_adjust_aux < adjust_aux; num_adjust_aux++) {
				if (updateStatus) IJ.showStatus("Building basic  DSI for the AUX camera image set "+quadCLT_main.image_name+
						" using main camera DSI, pass "+(num_adjust_aux+1)+" of "+num_adjust_aux);
				if (debugLevel > -5) {
					System.out.println("Building basic  DSI for the AUX camera image set "+quadCLT_main.image_name+
							" using main camera DSI, pass "+(num_adjust_aux+1)+" of "+num_adjust_aux);
				}
				quadCLT_aux.preExpandCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
						imp_srcs_aux, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
						saturation_imp_aux, // boolean [][] saturation_imp, // (near) saturated pixels or null
						clt_parameters,
						debayerParameters,
						colorProcParameters_aux,
						rgbParameters,
						threadsMax,  // maximal number of threads to launch
						updateStatus,
						debugLevelInner);
				// adjust extrinsics here
				System.out.println("Adjust AUX extrinsics here");
				if (updateStatus) IJ.showStatus("Adjusting AUX camera image set for "+quadCLT_aux.image_name+
						", pass "+(num_adjust_aux+1)+" of "+adjust_aux);
				if (debugLevel > -5) {
					System.out.println("Adjusting AUX camera image set for "+quadCLT_aux.image_name+
							", pass "+(num_adjust_aux+1)+" of "+adjust_aux);
				}
				if (quadCLT_aux.ds_from_main == null) {
					System.out.println("BUG: quadCLT_aux.ds_from_main should be not null here!");
				    double inf_min = -1.0;
				    double inf_max =  1.0;
				    if (num_adjust_aux >= (adjust_aux/2)) {
				        inf_min = -0.2;
				        inf_max = 0.05;
				    }
					// adjust w/o main camera - maybe will be used in the future
					boolean ok = quadCLT_aux.extrinsicsCLT(
							clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
							false, // adjust_poly,
				            inf_min, // double inf_min,
				            inf_max,  // double inf_max,
							threadsMax,  //final int        threadsMax,  // maximal number of threads to launch
							updateStatus,// final boolean    updateStatus,
							debugLevelInner); // final int        debugLevel)
					if (!ok) break;
				}
				boolean ok = quadCLT_aux.extrinsicsCLTfromGT(
						  null,
						  quadCLT_aux.ds_from_main, // gt_disp_strength,
						  clt_parameters,           // EyesisCorrectionParameters.CLTParameters           clt_parameters,
						  false,                    // adjust_poly,
						  threadsMax,               // final int        threadsMax,  // maximal number of threads to launch
						  updateStatus,             // final boolean    updateStatus,
						  debugLevel + 2);          // final int        debugLevel)
// clear memory for AUX
				quadCLT_aux.tp.resetCLTPasses();
				if (!ok) break;
			}
			// Generate 4 AUX camera images and thumbnail
			if (quadCLT_main.correctionsParameters.clt_batch_4img_aux){
				if (updateStatus) IJ.showStatus("Rendering 4 AUX image set (disparity = 0) for "+quadCLT_aux.image_name);

				quadCLT_aux.processCLTQuadCorrCPU( // returns ImagePlus, but it already should be saved/shown
						imp_srcs_aux, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
						saturation_imp_aux, // boolean [][] saturation_imp, // (near) saturated pixels or null
						  clt_parameters,
						  debayerParameters,
						  colorProcParameters_aux,
						  channelGainParameters,
						  rgbParameters,
						  scaleExposures_aux,
						  false, // calculate and apply additional fine geometry correction
						  false, // calculate and apply geometry correction at infinity
						  threadsMax,  // maximal number of threads to launch
						  updateStatus,
						  debugLevel);
				quadCLT_aux.tp.resetCLTPasses();
			}
			// Currently - no LWIR 3D model generation, maybe it will be added later
			// Generate AUX DS
			if (quadCLT_main.correctionsParameters.clt_batch_dsi_aux) {
				if (updateStatus) IJ.showStatus("Building basic DSI for the aux camera image set "+quadCLT_main.image_name+" (for DSI export)");
				quadCLT_aux.preExpandCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
						imp_srcs_aux, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
						saturation_imp_aux, // boolean [][] saturation_imp, // (near) saturated pixels or null
						clt_parameters,
						debayerParameters,
						colorProcParameters_aux,
						rgbParameters,
						threadsMax,  // maximal number of threads to launch
						updateStatus,
						debugLevelInner);
				if (updateStatus) IJ.showStatus("Expanding DSI for the aux camera image set "+quadCLT_main.image_name+" (for DSI export)");
				quadCLT_aux.expandCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
						clt_parameters,
						debayerParameters,
						colorProcParameters_aux,
						channelGainParameters,
						rgbParameters,
						threadsMax,  // maximal number of threads to launch
						updateStatus,
						debugLevel);
				double [][] aux_last_scan = quadCLT_aux.tp.getShowDS(
						quadCLT_aux.tp.clt_3d_passes.get( quadCLT_aux.tp.clt_3d_passes.size() -1),
						false); // boolean force_final);

//				dsi[DSI_DISPARITY_AUX] = aux_last_scan[0]; // incompatible dimensions
//				dsi[DSI_STRENGTH_AUX] =  aux_last_scan[1]; // incompatible dimensions
				dsi_aux_from_main[QuadCLT.FGBG_AUX_DISP] = aux_last_scan[0];
				dsi_aux_from_main[QuadCLT.FGBG_AUX_STR] =  aux_last_scan[1];
				quadCLT_main.saveDSIGTAux( // GT from main and AUX DS
						quadCLT_aux,
						dsi_aux_from_main);
				quadCLT_aux.tp.resetCLTPasses();
			}
			//

			// TODO: Add new ML generation here
			if (quadCLT_main.correctionsParameters.clt_batch_genMl) { // rig.ml_generate) { //clt_batch_genMl
				outputMLDataLwir(
						quadCLT_aux,    // QuadCLT                                  quadCLT_aux,
						clt_parameters, // EyesisCorrectionParameters.CLTParameters clt_parameters,
						null,           //String                                   ml_directory,       // full path or null (will use config one)
						threadsMax,     // final int                                threadsMax,  // maximal number of threads to launch
						updateStatus,   // final boolean                            updateStatus,
						debugLevel);    // final int                                debugLevel)
			}
			
			

			if (quadCLT_main.correctionsParameters.clt_batch_save_extrinsics) {
				saveProperties( // uses global quadCLT_main
						null,        // String path,                // full name with extension or w/o path to use x3d directory
						null,        // Properties properties,      // if null - will only save extrinsics)
						debugLevel);
				quadCLT_main.saveInterProperties( // save properties for interscene processing (extrinsics, ers, ...)
						null, // String path,             // full name with extension or w/o path to use x3d directory
//						null, // Properties properties,   // if null - will only save extrinsics)
						debugLevel);
				
			}
			if (quadCLT_main.correctionsParameters.clt_batch_save_all) {
				saveProperties(
						null,        // String path,                // full name with extension or w/o path to use x3d directory
						properties,  // Properties properties,    // if null - will only save extrinsics)
						debugLevel);
			}

			Runtime.getRuntime().gc();
			if (debugLevel >-1) System.out.println("Processing set "+(nSet+1)+" (of "+set_channels_aux.length+") finished at "+
					IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");

			if (quadCLT_aux.eyesisCorrections.stopRequested.get()>0) {
				System.out.println("User requested stop");
				System.out.println("Processing "+(nSet + 1)+" file sets (of "+set_channels_main.length+") finished at "+
						IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
				return;
			}
		}
		System.out.println("batchLwirRig(): processing "+(quadCLT_main.getTotalFiles(set_channels_main)+quadCLT_aux.getTotalFiles(set_channels_aux))+" files ("+set_channels_main.length+" file sets) finished at "+
				IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
	}

	public void showDSI()
	{
		  quadCLT_main.showDSI(dsi);
	}

/*
	public void saveDSIMain(
			QuadCLT quadCLT,
			double [][] dsi) // DSI_SLICES.length
	{
		String x3d_path= quadCLT.correctionsParameters.selectX3dDirectory( // for x3d and obj
				quadCLT.correctionsParameters.getModelName(quadCLT.image_name), // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
				quadCLT.correctionsParameters.x3dModelVersion,
				true,  // smart,
				true);  //newAllowed, // save
		String title = quadCLT.image_name+"-DSI_MAIN";
		String []   titles =   {DSI_SLICES[DSI_DISPARITY_MAIN], DSI_SLICES[DSI_STRENGTH_MAIN]};
		double [][] dsi_main = {dsi[DSI_DISPARITY_MAIN],        dsi[DSI_STRENGTH_MAIN]};

		ImagePlus imp = (new ShowDoubleFloatArrays()).makeArrays(dsi_main,quadCLT.tp.getTilesX(), quadCLT.tp.getTilesY(),  title, titles);
		quadCLT.eyesisCorrections.saveAndShow(
				imp,      // ImagePlus             imp,
				x3d_path, // String                path,
				false,    // boolean               png,
				false,    // boolean               show,
				0);       // int                   jpegQuality)
	}


	// Save GT from main and AUX calculated DS
	public void saveDSIGTAux(
			QuadCLT quadCLT_main,
			QuadCLT quadCLT_aux,
			double [][] dsi_aux_from_main)
	{
		String x3d_path= quadCLT_main.correctionsParameters.selectX3dDirectory( // for x3d and obj
				quadCLT_main.correctionsParameters.getModelName(quadCLT_main.image_name), // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
				quadCLT_main.correctionsParameters.x3dModelVersion,
				true,  // smart,
				true);  //newAllowed, // save
		String title = quadCLT_aux.image_name+"-DSI_GT-AUX";
//		String []   titles =   {DSI_SLICES[DSI_DISPARITY_MAIN], DSI_SLICES[DSI_STRENGTH_MAIN]};
//		double [][] dsi_main = {dsi[DSI_DISPARITY_MAIN],        dsi[DSI_STRENGTH_MAIN]};

		ImagePlus imp = (new ShowDoubleFloatArrays()).makeArrays(
				dsi_aux_from_main, // dsi_main,
				quadCLT_aux.tp.getTilesX(),
				quadCLT_aux.tp.getTilesY(),
				title,
				QuadCLT.FGBG_TITLES_AUX); // titles);
		quadCLT_main.eyesisCorrections.saveAndShow(
				imp,      // ImagePlus             imp,
				x3d_path, // String                path,
				false,    // boolean               png,
				false,    // boolean               show,
				0);       // int                   jpegQuality)
	}

	public void showDSIMain(
			QuadCLT quadCLT_main,
			double [][] dsi)
	{
		  String title = quadCLT_main.image_name+"-DSI_MAIN";
		  String []   titles =   {DSI_SLICES[DSI_DISPARITY_MAIN], DSI_SLICES[DSI_STRENGTH_MAIN]};
		  double [][] dsi_main = {dsi[DSI_DISPARITY_MAIN],        dsi[DSI_STRENGTH_MAIN]};

		  (new ShowDoubleFloatArrays()).showArrays(dsi_main,quadCLT_main.tp.getTilesX(), quadCLT_main.tp.getTilesY(), true, title, titles);
	}
*/


	public double [][] getRigDSI(
			String                                         path_DSI,   // Combo DSI path
			boolean                                        main) // false - rig
	{
		int [] slices_rig =  {TwoQuadCLT.DSI_DISPARITY_RIG,  TwoQuadCLT.DSI_STRENGTH_RIG};
		int [] slices_main = {TwoQuadCLT.DSI_DISPARITY_MAIN, TwoQuadCLT.DSI_STRENGTH_MAIN};
		int [] slices = main? slices_main:slices_rig;
		double[][] dsi = new double [slices.length][];
		ImagePlus imp_dsi=new ImagePlus(path_DSI);
		ImageStack dsi_stack=  imp_dsi.getStack();
		int nLayers = dsi_stack.getSize();
		for (int nl = 0; nl < nLayers; nl++){
			for (int i = 0; i < slices.length; i++) {
				if (TwoQuadCLT.DSI_SLICES[slices[i]].equals(dsi_stack.getSliceLabel(nl + 1))) {
					dsi[i] = new double [((float[]) dsi_stack.getPixels(nl + 1)).length];
					for (int j = 0; j < dsi[i].length; j++) {
						dsi[i][j] = ((float[]) dsi_stack.getPixels(nl + 1))[j];
					}
					break;
				}
			}
		}
		return dsi;
	}
	/**
	 * Enhance main DSI data to prepare ML files:
	 * 1. remove tiles that are too different from the rig (> rig_tolerance) - that may me different object
	 * 2. replace undefined tiles that have rig data with average of 8 neighbors within
	 *    rig_tolerance from rig disparity. If there are no suitable neighbors, use
	 *    random offset from rig disparity within +/- rnd_offset
	 * 3. grow defined data, each layer using min/max/average from the known neighbors
	 * 4. assign new_strength to previously undefined tiles
	 * @param main_dsi
	 * @param rig_dsi
	 * @param rig_tolerance
	 * @param rnd_offset
	 * @param main_tolerance
	 * @param grow_steps
	 * @param grow_mode : -1 - min, 0 - average, +1 - max
	 * @param new_strength
	 * @return 'enhanced' disparity/strength for the main camera
	 */

	public double [][] enhanceMainDSI(
			double [][] main_dsi,
			double [][] rig_dsi,
			double      rig_tolerance,
			double      rnd_offset,
			double      main_tolerance,
			int         grow_steps,
			int         grow_mode,
			double      inher_strength,
			int         debugLevel)
	{
		double [][] dbg_img = null;
		int num_dbg = 0;
		String [] dbg_titles = null;
		if (debugLevel > 1) {
			num_dbg = grow_steps + 4;
			dbg_img = new double [2 * num_dbg][];
			dbg_titles = new String [2 * num_dbg];
			dbg_img[0]           = main_dsi[0].clone();
			dbg_img[0 + num_dbg] = main_dsi[1].clone();
			dbg_titles[0] =           "main_d";
			dbg_titles[0 + num_dbg] = "main_s";
			dbg_img[1]           = rig_dsi [0].clone();
			dbg_img[1 + num_dbg] = rig_dsi [1].clone();
			dbg_titles[1] =           "rig_d";
			dbg_titles[1 + num_dbg] = "rig_s";
		}
		final int tilesX = this.quadCLT_main.tp.getTilesX();
		final int tilesY = this.quadCLT_main.tp.getTilesY();
		TileNeibs         tnImage =  new TileNeibs(tilesX, tilesY); // biCamDSI_persistent.tnImage;

		double [][]       enh_dsi = {main_dsi[0].clone(),main_dsi[1].clone()};
		int num_too_far =      0;
		int num_by_neib =      0;
		int num_rnd =          0;
		int num_extrapolated = 0;

		for (int nTile = 0; nTile < enh_dsi[0].length;nTile++) {
			// Remove from main camera measurements those that do not have rig data or differ from rig by more than rig_tolerance
			if (Double.isNaN(rig_dsi[0][nTile]) || (Math.abs(rig_dsi[0][nTile] - enh_dsi[0][nTile]) > rig_tolerance)) {
				enh_dsi[1][nTile] = 0.0;
				enh_dsi[0][nTile] = Double.NaN;
				if (!Double.isNaN(rig_dsi[0][nTile])) {
					num_too_far++; // count only too far, not NaN-s
				}
			}
		}
		double [] new_disp = enh_dsi[0].clone();
		double [] new_str =  enh_dsi[1].clone();
		Random rnd = new Random(System.nanoTime());
		for (int nTile = 0; nTile < enh_dsi[0].length;nTile++) if (Double.isNaN(enh_dsi[0][nTile]) && !Double.isNaN(rig_dsi[0][nTile])){
			double sw = 0.0;
			double sdw = 0.0;
			int nneibs = 0;
			for (int dir = 0; dir < 8; dir++) {
				int nTile1 = tnImage.getNeibIndex(nTile, dir);
				if (    (nTile1 >= 0) &&
						!Double.isNaN(enh_dsi[0][nTile1]) &&
						(Math.abs(enh_dsi[0][nTile1] - rig_dsi[0][nTile]) <= main_tolerance)) {
					double w = enh_dsi[1][nTile1];
					sw +=  w;
					sdw += w * enh_dsi[0][nTile1];
					nneibs++;
				}
			}
			if (sw > 0) {
				new_disp[nTile] = sdw/sw;
				new_str[nTile] = inher_strength * sw/nneibs;
				num_by_neib++;
			} else {
				new_disp[nTile] = rig_dsi[0][nTile] + rnd_offset*(2 * rnd.nextDouble() - 1.0);
				new_str[nTile] = inher_strength * rig_dsi[1][nTile];
				num_rnd++;
			}

		}
		for (int nTile = 0; nTile < enh_dsi[0].length;nTile++) if (Double.isNaN(enh_dsi[0][nTile]) && !Double.isNaN(rig_dsi[0][nTile])){
			enh_dsi[0][nTile] = new_disp[nTile];
			enh_dsi[1][nTile] = new_str[nTile];
		}
		if (dbg_img != null) {
			dbg_img[2]           = enh_dsi[0].clone();
			dbg_img[2 + num_dbg] = enh_dsi[1].clone();
			dbg_titles[2] =           "enh_preexp_d";
			dbg_titles[2 + num_dbg] = "enh_preexp_s";
		}
		//		double inher_strength = 0.5;
		boolean [] exp_full = null;
		if (grow_steps > 0) for (int n_expand = 0; n_expand < (grow_steps+1); n_expand++) {
			boolean [] selection = new boolean[enh_dsi[0].length];
			for (int nTile = 0; nTile < enh_dsi[0].length;nTile++) {
				selection[nTile] = !Double.isNaN(enh_dsi[0][nTile]);
			}
			if (n_expand == (grow_steps-1)) {
				exp_full = selection.clone();
				tnImage.growSelection(2, exp_full, null); // hor, vert and diagonal
				tnImage.growSelection(1, selection, null); // only hor/vert
			} else if ((n_expand == grow_steps) && (exp_full != null)) {
				selection = exp_full;
			} else {
				tnImage.growSelection(1, selection, null); // only hor/vert
			}

//			tnImage.growSelection(2, selection, null);
			for (int nTile = 0; nTile < enh_dsi[0].length;nTile++) if (selection[nTile] && Double.isNaN((enh_dsi[0][nTile]))){
				double sw = 0.0;
				double sdw = 0.0;
				int ntiles = 0;
				double d_min = Double.POSITIVE_INFINITY;
				double d_max = Double.NEGATIVE_INFINITY;
				for (int dir = 0; dir < 8; dir++) {
					int nTile1 = tnImage.getNeibIndex(nTile, dir);
					if ((nTile1 >= 0) && !Double.isNaN(enh_dsi[0][nTile1])){
						double d = enh_dsi[0][nTile1];
						if (d < d_min) {
							d_min = d;
						}
						if (d > d_max) {
							d_max = d;
						}
						double w = enh_dsi[1][nTile1];
						sw +=  w;
						sdw += w * enh_dsi[0][nTile1];
						ntiles++;
					}
				}
				if (sw > 0.0) {
					double d_avg = sdw/sw;
					switch (grow_mode) {
					case -1: new_disp[nTile] = d_min; break;
					case  1: new_disp[nTile] = d_max; break;
					case  0: new_disp[nTile] = d_avg;   break;
					default: {
						// find - which of the min,max, avg is closer to rig (if available, if not - use avg)
						if (!Double.isNaN(rig_dsi[0][nTile])) {
							double [] diffs = {Math.abs(d_min - rig_dsi[0][nTile]), Math.abs(d_max - rig_dsi[0][nTile]),Math.abs(d_avg - rig_dsi[0][nTile])};
							if ((diffs[0] < diffs[1]) && (diffs[0] < diffs[2])){
								new_disp[nTile] = d_min;
							} else if (diffs[1] < diffs[2]) {
								new_disp[nTile] = d_max;
							} else {
								new_disp[nTile] = d_avg;
							}
						} else {
							new_disp[nTile] = d_avg;
						}
					}
					}
					new_str[nTile] = inher_strength * sw/ntiles;
					num_extrapolated++;
				} else {
					new_disp[nTile] = Double.NaN;
					new_str[nTile] = 0.0;
				}

			}
			for (int nTile = 0; nTile < enh_dsi[0].length;nTile++) if (selection[nTile] && Double.isNaN(enh_dsi[0][nTile]) && !Double.isNaN(new_disp[nTile]) ){
				enh_dsi[0][nTile] = new_disp[nTile];
				enh_dsi[1][nTile] = new_str[nTile]; // new_strength;
			}
			if (dbg_img != null) {
				dbg_img[3+n_expand]           = enh_dsi[0].clone();
				dbg_img[3+n_expand + num_dbg] = enh_dsi[1].clone();
				dbg_titles[3+n_expand] =           "enh_exp"+n_expand+"_d";
				dbg_titles[3+n_expand + num_dbg] = "enh_exp"+n_expand+"_s";
			}

		}
		if (debugLevel > 0) {
			System.out.println("enhanceMainDSI(): num_too_far="+num_too_far+", num_by_neib="+num_by_neib+", num_rnd="+num_rnd+", num_extrapolated="+num_extrapolated);
		}
		if (dbg_img != null) {
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					tilesX,
					tilesY,
					true,
					"enhanced_dsi_dbg",
					dbg_titles);

		}
		return enh_dsi;
	}


	public void regenerateML(
			String                                         path_DSI,   // Combo DSI path
			String                                         model_dir,  // model/version directory
			String                                         ml_subdir,  // new ML subdir (or null to use configuartion one)
			QuadCLT                                        quadCLT_main,
			QuadCLT                                        quadCLT_aux,
			CLTParameters       clt_parameters,
			EyesisCorrectionParameters.DebayerParameters   debayerParameters,
			ColorProcParameters                            colorProcParameters,
			ColorProcParameters                            colorProcParameters_aux,
			EyesisCorrectionParameters.RGBParameters       rgbParameters,
			final int                                      threadsMax,  // maximal number of threads to launch
			final boolean                                  updateStatus,
			final int                                      debugLevel) throws Exception
	{
		this.startTime=System.nanoTime();
		boolean new_dir_only = true;
		String ml_dir = (ml_subdir == null)? null: (model_dir+Prefs.getFileSeparator()+ml_subdir);
		if (new_dir_only) {
			File dir = new File(ml_dir);
			if (dir.exists()){
				System.out.println("Directory already exists: "+dir+", and new_dir_only=true -> skipping generation");
				return;
			}
		}

		double [][] rig_dsi =  getRigDSI(path_DSI, false);
		double [][] main_dsi = getRigDSI(path_DSI, true);

		String [] sourceFiles=quadCLT_main.correctionsParameters.getSourcePaths();
		QuadCLT.SetChannels [] set_channels_main = quadCLT_main.setChannels(debugLevel);
		QuadCLT.SetChannels [] set_channels_aux =  quadCLT_aux.setChannels(debugLevel);
		if ((set_channels_main == null) || (set_channels_main.length==0) || (set_channels_aux == null) || (set_channels_aux.length==0)) {
			System.out.println("No files to process (of "+sourceFiles.length+")");
			return;
		}
		double [] referenceExposures_main = quadCLT_main.eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		double [] referenceExposures_aux =  quadCLT_aux.eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		if ((set_channels_main.length != 1) || (set_channels_aux.length != 1)) {
			System.out.println("Wrong number of source file sets, expecting 1, got "+set_channels_main.length+", "+set_channels_aux.length);
		}
		int nSet = 0;
		if (set_channels_aux.length <= nSet ) {
			throw new Exception ("Set names for cameras do not match: main camera: '"+set_channels_main[nSet].name()+"', aux. camera: nothing");
		}
		if (!set_channels_main[nSet].name().equals(set_channels_aux[nSet].name())) {
			throw new Exception ("Set names for cameras do not match: main camera: '"+set_channels_main[nSet].name()+"', aux. camera: '"+set_channels_main[nSet].name()+"'");
		}

		int [] channelFiles_main = set_channels_main[nSet].fileNumber();
		int [] channelFiles_aux =  set_channels_aux[nSet].fileNumber();
		boolean [][] saturation_main = (clt_parameters.sat_level > 0.0)? new boolean[channelFiles_main.length][] : null;
		boolean [][] saturation_aux =  (clt_parameters.sat_level > 0.0)? new boolean[channelFiles_main.length][] : null;
		double [] scaleExposures_main = new double[channelFiles_main.length];
		double [] scaleExposures_aux =  new double[channelFiles_main.length];

		ImagePlus [] imp_quad_main = quadCLT_main.conditionImageSet(
				clt_parameters,                 // EyesisCorrectionParameters.CLTParameters  clt_parameters,
				colorProcParameters,            //  ColorProcParameters                       colorProcParameters, //
				sourceFiles,                    // String []                                 sourceFiles,
				set_channels_main[nSet].name(), // String                                    set_name,
				referenceExposures_main,        // double []                                 referenceExposures,
				channelFiles_main,              // int []                                    channelFiles,
				scaleExposures_main,            //output  // double [] scaleExposures
				saturation_main,                //output  // boolean [][]                              saturation_imp,
				threadsMax,                 // int                                       threadsMax,
				debugLevel);                   // int                                       debugLevel);

		ImagePlus [] imp_quad_aux = quadCLT_aux.conditionImageSet(
				clt_parameters,                 // EyesisCorrectionParameters.CLTParameters  clt_parameters,
				colorProcParameters_aux,        //  ColorProcParameters                       colorProcParameters, //
				sourceFiles,                    // String []                                 sourceFiles,
				set_channels_aux[nSet].name(),  // String                                    set_name,
				referenceExposures_aux,         // double []                                 referenceExposures,
				channelFiles_aux,               // int []                                    channelFiles,
				scaleExposures_aux,             //output  // double [] scaleExposures
				saturation_aux,                 //output  // boolean [][]                              saturation_imp,
				threadsMax,                 // int                                       threadsMax,
				debugLevel);                    // int                                       debugLevel);
/*		// 08/12/2020 Moved to conditionImageSet
		getRigImageStacks(
				clt_parameters,  // EyesisCorrectionParameters.CLTParameters       clt_parameters,
				quadCLT_main,    // QuadCLT                                         quadCLT_main,
				quadCLT_aux,     // QuadCLT                                          quadCLT_aux,
				imp_quad_main,   // ImagePlus []                                   imp_quad_main,
				imp_quad_aux,    // ImagePlus []                                    imp_quad_aux,
				saturation_main, // boolean [][]        saturation_main, // (near) saturated pixels or null
				saturation_aux, // boolean [][]        saturation_aux, // (near) saturated pixels or null
				threadsMax,      // maximal number of threads to launch
				debugLevel);     // final int        debugLevel);
		// now tp is defined
*/ 

/*
		main_dsi = enhanceMainDSI(
				main_dsi,                             // double [][] main_dsi,
				rig_dsi,                              // double [][] rig_dsi,
				clt_parameters.rig.ml_rig_tolerance,  // double      rig_tolerance,
				clt_parameters.rig.ml_rnd_offset,     // double      rnd_offset,
				clt_parameters.rig.ml_main_tolerance, // double      main_tolerance,
				clt_parameters.rig.ml_grow_steps,     // int         grow_steps,
				clt_parameters.rig.ml_grow_mode,      // int         grow_mode,
				clt_parameters.rig.ml_new_strength,   // double      new_strength)
//				debugLevel + 3);                        // int         debugLevel);
				debugLevel + 0);                        // int         debugLevel);
*/
		quadCLT_main.tp.rig_pre_poles_ds = rig_dsi;  // use rig data from the COMBO-DSI file
		quadCLT_main.tp.main_ds_ml =       main_dsi; // use rig data from the COMBO-DSI file

//		quadCLT_main.tp.resetCLTPasses();
//		quadCLT_aux.tp.resetCLTPasses();
//		String ml_dir = (ml_subdir == null)? null: (model_dir+Prefs.getFileSeparator()+ml_subdir);
		outputMLData(
				quadCLT_main,   // QuadCLT                                  quadCLT_main,  // tiles should be set
				quadCLT_aux,    // QuadCLT                                  quadCLT_aux,
				clt_parameters, // EyesisCorrectionParameters.CLTParameters clt_parameters,
				ml_dir,         // String                                   ml_directory,       // full path or null (will use config one)
				threadsMax,     // final int                                threadsMax,  // maximal number of threads to launch
				updateStatus,   // final boolean                            updateStatus,
				debugLevel);    // final int                                debugLevel)
		System.out.println();
	}


	public void saveProperties(
			String path,                // full name with extension or w/o path to use x3d directory
			Properties properties,    // if null - will only save extrinsics)
			int debugLevel)
	{
		// update properties from potentially modified parameters (others should be updated
		if (path == null) {
			path = quadCLT_main.image_name + ((properties == null) ? "-EXTRINSICS":"")+".corr-xml";

		}
		if (!path.contains(Prefs.getFileSeparator())) {
			  String x3d_path= quadCLT_main.correctionsParameters.selectX3dDirectory( // for x3d and obj
					  quadCLT_main.correctionsParameters.getModelName(quadCLT_main.image_name), // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
					  quadCLT_main.correctionsParameters.x3dModelVersion,
						  true,  // smart,
						  true);  //newAllowed, // save
			  path = x3d_path+Prefs.getFileSeparator()+path;
		}
		if (properties == null) {
			properties = new Properties();
		}
		quadCLT_main.setProperties(QuadCLT.PREFIX,properties);
		quadCLT_aux.setProperties(QuadCLT.PREFIX_AUX,properties);
		OutputStream os;
		try {
			os = new FileOutputStream(path);
		} catch (FileNotFoundException e1) {
			// missing config directory
			File dir = (new File(path)).getParentFile();
			if (!dir.exists()){
				dir.mkdirs();
				try {
					os = new FileOutputStream(path);
				} catch (FileNotFoundException e2) {
					IJ.showMessage("Error","Failed to create directory "+dir.getName()+" to save configuration file: "+path);
					return;
				}
			} else {
				IJ.showMessage("Error","Failed to open configuration file: "+path);
				return;
			}
		}
		try {
			properties.storeToXML(os,
					"last updated " + new java.util.Date(), "UTF8");

		} catch (IOException e) {
			IJ.showMessage("Error","Failed to write XML configuration file: "+path);
			return;
		}
		try {
			os.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		if (debugLevel> -3) {
			System.out.println("Configuration parameters are saved to "+path);
		}
	}
}
