/**
 **
 ** LwirWorld.java - Incrementally build LWIR 3D world(s)
 **
 ** Copyright (C) 2022 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  LwirWorld.java is free software: you can redistribute it and/or modify
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

package com.elphel.imagej.tileprocessor.lwoc;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Properties;

import com.elphel.imagej.cameras.CLTParameters;
import com.elphel.imagej.cameras.ColorProcParameters;
import com.elphel.imagej.cameras.EyesisCorrectionParameters;
import com.elphel.imagej.common.MultiThreading;
import com.elphel.imagej.correction.CorrectionColorProc;
import com.elphel.imagej.tileprocessor.CLTPass3d;
import com.elphel.imagej.tileprocessor.OpticalFlow;
import com.elphel.imagej.tileprocessor.QuadCLT;
import com.elphel.imagej.tileprocessor.TileProcessor;
import com.elphel.imagej.tileprocessor.TwoQuadCLT;

import ij.IJ;

public class LwirWorld {
	public static double [] ZERO3 = new double[3];
    public static void buildWorld(
			QuadCLT                                              quadCLT_main, // tiles should be set
			int                                                  ref_index_unused, // -1 - last
			int                                                  ref_step_unused, // not used here
			CLTParameters             clt_parameters,
			EyesisCorrectionParameters.DebayerParameters         debayerParameters,
			ColorProcParameters                                  colorProcParameters,
			CorrectionColorProc.ColorGainsParameters             channelGainParameters,
			EyesisCorrectionParameters.RGBParameters             rgbParameters,
			EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
			Properties                                           properties,
			boolean                                              reset_from_extrinsics,
			final int        threadsMax,  // maximal number of threads to launch
			final boolean    updateStatus,
			final int        debugLevel)  throws Exception
    {
    	System.out.println("buildWorld");
		double    world_hsize =      clt_parameters.lwp.world_hsize;			
		double    min_hsize =        clt_parameters.lwp.min_hsize;
		double    max_hsize =        clt_parameters.lwp.max_hsize;			
		int       max_mesh_centers = clt_parameters.lwp.max_mesh_centers;
		int       max_cameras =      clt_parameters.lwp.max_cameras;    	
    	int     min_num_scenes =     clt_parameters.imp.min_num_scenes; // abandon series if there are less than this number of scenes in it 
    	if (min_num_scenes < 1) {
    		min_num_scenes = 1;
    	}
    	long start_time_all = System.nanoTime();
    	
		EyesisCorrectionParameters.CorrectionParameters.PathFirstLast[] pathFirstLast = null;
		int num_seq = 1;
		if (quadCLT_main.correctionsParameters.useSourceList) {
			pathFirstLast = quadCLT_main.correctionsParameters.getSourceSets(
					quadCLT_main.correctionsParameters.sourceSequencesList);
			if (pathFirstLast != null) {
				num_seq = pathFirstLast.length; 
			}
		}
		LwocWorld.initWorlds();
		LwocWorld.newWorld(
				ZERO3,            // double [] center,
				world_hsize,      // double    world_hsize,			
				min_hsize,        // double    min_hsize,
				max_hsize,        // double    max_hsize,			
				max_mesh_centers, // int       max_mesh_centers,
				max_cameras,      // int       max_cameras,
				ZERO3,            // double [] atr,
				ZERO3);           // double [] xyz) 
		// LwocWorld.getWorld()
		for (int nseq = 0; nseq < num_seq; nseq++) { // 99-scene series
			long start_time_seq = System.nanoTime();
			System.out.println("\nSTARTED PROCESSING SCENE SEQUENCE "+nseq+" (last is "+(num_seq-1)+")");
			if (pathFirstLast != null) {
				File [] scene_dirs = (new File(pathFirstLast[nseq].path)).listFiles(); // may contain non-directories, will be filtered by filterScenes
				quadCLT_main.correctionsParameters.filterScenes(
						scene_dirs, // File [] scene_dirs,
						pathFirstLast[nseq].first, // int scene_first, // first scene to process
						pathFirstLast[nseq].last); // int scene_last);  // last scene to process (negative - add length
				
				// Now quadCLT_main.correctionsParameters.sourcePaths[] contain 16*99 full paths
				
				if (pathFirstLast[nseq].movement_size < 0) {
					clt_parameters.imp.mov_en = false;
					if (debugLevel > -4) {
						System.out.println("Disabling movement detection for this scene.");
					}
				} else {
					clt_parameters.imp.mov_en = true;
					clt_parameters.imp.mov_max_len = pathFirstLast[nseq].movement_size;
					if (debugLevel > -4) {
						System.out.println("Enabling movement detection for this scene with maximum cluster linear size of "+
								clt_parameters.imp.mov_max_len+" tiles.");
					}
				}
			}
			int ref_index = -1; // -1 - last
			int [] start_ref_pointers = new int[2];
			while ((ref_index < 0) || ((ref_index + 1) >= min_num_scenes)) {
				String model_directory = buildSeries(
						(pathFirstLast != null),   //boolean                                              batch_mode,
						quadCLT_main,              // QuadCLT                                              quadCLT_main, // tiles should be set
						ref_index,                 // int                                                  ref_index, // -1 - last
						clt_parameters,            // CLTParameters             clt_parameters,
						debayerParameters,         // EyesisCorrectionParameters.DebayerParameters         debayerParameters,
						colorProcParameters,       // ColorProcParameters                                  colorProcParameters,
						channelGainParameters,     // CorrectionColorProc.ColorGainsParameters             channelGainParameters,
						rgbParameters,             // EyesisCorrectionParameters.RGBParameters             rgbParameters,
						equirectangularParameters, // EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
						properties,                // Properties                                           properties,
						reset_from_extrinsics,     // boolean                                              reset_from_extrinsics,
						start_ref_pointers,        // int []                                               start_ref_pointers,
						updateStatus, // final boolean    updateStatus,
						debugLevel+2); // final int        debugLevel)
				if (model_directory == null) {
					System.out.println("Failed to build sequence for series "+ref_index);
					break; // and go to the to next scene sequence from the list
				}
				String series_action = (start_ref_pointers[0] < (min_num_scenes-1))?"is FINISHED ":("will continue down from scene "+(start_ref_pointers[0]));
				System.out.println("PROCESSING SCENE SEQUENCE "+nseq+" (last is "+(num_seq-1)+") "+series_action+" in "+
						IJ.d2s(0.000000001*(System.nanoTime()-start_time_seq),3)+" sec ("+
						IJ.d2s(0.000000001*(System.nanoTime()-start_time_all),3)+" sec from the overall start");
				// will open dialog if does not exist
				String linkedModelsDirectory = quadCLT_main.correctionsParameters.selectLinkedModelsDirectory(true,true);
				if ((linkedModelsDirectory != null) && (linkedModelsDirectory.length() > 0)) {
					Path pathAbsolute = Paths.get(model_directory);
					Path pathBase = Paths.get(linkedModelsDirectory);
					Path pathRelative = pathBase.relativize(pathAbsolute);
					File linkDir = new File(linkedModelsDirectory);
					linkDir.mkdirs();
					File link = new File(linkDir, pathAbsolute.getFileName().toString());
					if (link.exists()) {
						link.delete();
					}
					Files.createSymbolicLink(link.toPath(), pathRelative);
				}
				if (start_ref_pointers[0] < (min_num_scenes-1)) {
					break;
				}
				ref_index = start_ref_pointers[0]; // continue from the same attached to the previous reference
			}			
		}
    }	

    public static String buildSeries(
    		boolean                                              batch_mode,
			QuadCLT                                              quadCLT_main, // tiles should be set
			int                                                  ref_index, // -1 - last
			CLTParameters                                        clt_parameters,
			EyesisCorrectionParameters.DebayerParameters         debayerParameters,
			ColorProcParameters                                  colorProcParameters,
			CorrectionColorProc.ColorGainsParameters             channelGainParameters,
			EyesisCorrectionParameters.RGBParameters             rgbParameters,
			EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
			Properties                                           properties,
			boolean                                              reset_from_extrinsics,
			int []                                               start_ref_pointers, // [0] - earliest valid scene, [1] ref_index
			                                                     // each element is 0 for non-stereo and full width for stereo 
			final boolean    updateStatus,
			final int        debugLevel)  throws Exception
    {
    	boolean   build_ref_dsi =              clt_parameters.imp.force_ref_dsi;

    	final int debugLevelInner=             clt_parameters.batch_run? -2: debugLevel; // copied from TQ
//    	boolean  photo_en =                  clt_parameters.photo_en; //          false; // perform photogrammetric calibration to equalize pixel values
    	boolean  photo_each =                clt_parameters.photo_each; //        true;  // perform photogrammetric calibration to equalize pixel values
		boolean  photo_to_main=              clt_parameters.photo_to_main; // maybe it will not be needed, it will apply this calibration to the next scene sequence

    	int      photo_num_full =            clt_parameters.photo_num_full; //     1;    // Number of full recalibrations with re-processing of the images  
    	int      photo_num_refines =         clt_parameters.photo_num_refines; //  3;    // Calibrate, remove outliers, recalibrate, ... 
    	int      photo_min_good =            clt_parameters.photo_min_good;    //  1000;    // Minimal number of good pixels for photometric calibration 
    	double   photo_min_strength =        clt_parameters.photo_min_strength; // 0.0;  // maybe add to filter out weak tiles
    	double   photo_max_diff =            clt_parameters.photo_max_diff; //    40.0;  // To filter mismatches. Normal (adjusted) have RMSE ~9
    	int      photo_order =               clt_parameters.photo_order;
    	double   photo_std_1 =               clt_parameters.photo_std_1; //  50.0;   // Minimal standard deviation of the filtered values for poly order 1
    	double   photo_std_2 =               clt_parameters.photo_std_2; // 200.0;   // Minimal standard deviation of the filtered values for poly order 2
    	int      photo_offs_set =            clt_parameters.photo_offs_set; //  0;    // 0 - keep weighted offset average, 1 - balance result image, 2 - set weighted average to specific value
    	double   photo_offs =                clt_parameters.photo_offs;     // 21946; // weighted average offset target value, if photo_offs_set (and not photo_offs_balance)
//    	boolean  photo_debug =               clt_parameters.photo_debug; //       false; // Generate images and text

 		double sky_highest_min =     clt_parameters.imp.sky_highest_min;
 		double cold_frac =           clt_parameters.imp.cold_frac;
 		double hot_frac =            clt_parameters.imp.hot_frac;
 		double cold_scale =          clt_parameters.imp.cold_scale;
 		double sky_seed =            clt_parameters.imp.sky_seed;
 		double lma_seed =            clt_parameters.imp.lma_seed;
 		int    sky_shrink =          clt_parameters.imp.sky_shrink;
 		int    sky_bottleneck =      clt_parameters.imp.sky_bottleneck;
 		int    seed_rows =           clt_parameters.imp.seed_rows;
 		double sky_lim =             clt_parameters.imp.sky_lim;
 		int    sky_expand_extra =    clt_parameters.imp.sky_expand_extra;
 		double min_strength =        clt_parameters.imp.min_strength;
 		int    lowest_sky_row =      clt_parameters.imp.lowest_sky_row;
 		double sky_bottom_override = clt_parameters.imp.sky_bottom_override;
 		int    sky_override_shrink = clt_parameters.imp.sky_override_shrink;
		double disp_boost_min =      clt_parameters.imp.disp_boost_min; //  0.5;
		double disp_boost_diff  =    clt_parameters.imp.disp_boost_diff; // 0.35;
		int    disp_boost_neibs  =   clt_parameters.imp.disp_boost_neibs; // 2;
		double disp_boost_amount  =  clt_parameters.imp.disp_boost_amount; // 2.0;
    	
		boolean [] ref_blue_sky = null; // turn off "lma" in the ML output  
    	int earliest_scene =        0;      // increase on failure
    	
    	
		if ((quadCLT_main != null) && (quadCLT_main.getGPU() != null)) {
			quadCLT_main.getGPU().resetGeometryCorrection();
			quadCLT_main.gpuResetCorrVector(); // .getGPU().resetGeometryCorrectionVector();
		}
		// final boolean    batch_mode = clt_parameters.batch_run;
//		this.startTime=System.nanoTime();
		String [] sourceFiles0=quadCLT_main.correctionsParameters.getSourcePaths();
		QuadCLT.SetChannels [] set_channels_main = quadCLT_main.setChannels(debugLevel);
		if ((set_channels_main == null) || (set_channels_main.length==0)) {
			System.out.println("buildSeriesTQ(): No files to process (of "+sourceFiles0.length+")");
			return null;
		}
		// set_channels will include all 99 scenes even as quadCLTs.length matches ref_index
		QuadCLT.SetChannels [] set_channels=quadCLT_main.setChannels(debugLevel);
		
		// Each of (usually 99) set_channels entries has set name (timestamp) and an array of corresponding (16) file indices
		
		if (ref_index < 0) {
			ref_index += set_channels.length; 
		}
		if (start_ref_pointers != null) {
			start_ref_pointers[0] = 0;
			start_ref_pointers[1] = ref_index;
		}
		
		QuadCLT [] quadCLTs = new QuadCLT [ref_index+1]; //  [set_channels.length];
		
		//start_index
		
		double [][][] scenes_xyzatr =      new double [quadCLTs.length][][]; // previous scene relative to the next one
		scenes_xyzatr[ref_index] = new double[2][3]; // all zeros
		// See if build_ref_dsi is needed
		if (!build_ref_dsi) {
			quadCLTs[ref_index] = (QuadCLT) quadCLT_main.spawnNoModelQuadCLT( // will conditionImageSet
					set_channels[ref_index].set_name,
					clt_parameters,
					colorProcParameters, //
					MultiThreading.THREADS_MAX, // threadsMax,
					debugLevel-2);
			if ((quadCLTs[ref_index] == null) || !quadCLTs[ref_index].dsiExists()) {
				if (debugLevel >-2) {
					System.out.println("DSI data for scene "+set_channels[ref_index].set_name+" does not exist, forcing its calculation.");
				}
				build_ref_dsi = true;
			}
		}
		if (!build_ref_dsi && (quadCLTs[ref_index] != null)) {
			quadCLTs[ref_index].restoreInterProperties(null, false, debugLevel); //null
		}
		// 1. Reference scene DSI
		
		
/* */
		while ((quadCLTs[ref_index] == null) || (quadCLTs[ref_index].getBlueSky() == null)) { // null
			if (build_ref_dsi) {
				TwoQuadCLT.copyJP4src( // actually there is no sense to process multiple image sets. Combine with other
						// processing?
						set_channels[ref_index].set_name, // String set_name
						quadCLT_main, // QuadCLT quadCLT_main,
						null, // QuadCLT quadCLT_aux,
						clt_parameters, // EyesisCorrectionParameters.DCTParameters dct_parameters,
						true, // boolean                                  skip_existing,
						true, // false, // boolean                                  search_KML,
						debugLevel);
				quadCLTs[ref_index] = (QuadCLT) quadCLT_main.spawnNoModelQuadCLT( // will conditionImageSet
						set_channels[ref_index].set_name,
						clt_parameters,
						colorProcParameters, //
						MultiThreading.THREADS_MAX, // threadsMax,
						debugLevel-2);
				quadCLTs[ref_index].saveQuadClt(); // to re-load new set of Bayer images to the GPU (do nothing for CPU) and Geometry
				quadCLTs[ref_index].preExpandCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
						clt_parameters,
						debayerParameters,
						colorProcParameters,
						rgbParameters,
						MultiThreading.THREADS_MAX, // threadsMax,
						updateStatus,
						debugLevelInner);
				double [][] aux_last_scan = TileProcessor.getDSLMA(
						quadCLTs[ref_index].tp.clt_3d_passes.get( quadCLTs[ref_index].tp.clt_3d_passes.size() -1),
						false); // boolean force_final);
				double [][] dsi = new double [TwoQuadCLT.DSI_SLICES.length][];
				dsi[TwoQuadCLT.DSI_DISPARITY_AUX] =     aux_last_scan[0];
				dsi[TwoQuadCLT.DSI_STRENGTH_AUX] =      aux_last_scan[1];
				dsi[TwoQuadCLT.DSI_DISPARITY_AUX_LMA] = aux_last_scan[2];
				dsi[TwoQuadCLT.DSI_SPREAD_AUX] =        aux_last_scan[3];
				dsi[TwoQuadCLT.DSI_AVGVAL_AUX] =        aux_last_scan[4];

				if (quadCLT_main.correctionsParameters.clt_batch_dsi_cm_strength) {
					CLTPass3d scan = new CLTPass3d(quadCLTs[ref_index].tp);
					scan.setTileOpDisparity(aux_last_scan[0]); // measure w/o LMA, use just strength
					quadCLTs[ref_index].CLTMeas( // perform single pass according to prepared tiles operations and disparity
							clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
							scan,           // final CLTPass3d   scan,
							false,          // final boolean     save_textures,
							false,          // final boolean       need_diffs,     // calculate diffs even if textures are not needed 
							0,              // final int         clust_radius,
							true,           // final boolean     save_corr,
							false,          // final boolean     run_lma, // =    true;
							0.0,            // final double        max_chn_diff, // filter correlation results by maximum difference between channels
							-1.0,           // final double        mismatch_override, // keep tile with large mismatch if there is LMA with really strong correlation
							MultiThreading.THREADS_MAX, // threadsMax,
							updateStatus,   // final boolean     updateStatus,
							debugLevel-2);  // final int         debugLevel);
					dsi[TwoQuadCLT.DSI_STRENGTH_AUX] =      scan.getStrength();
					if (debugLevel > 1) {
						quadCLTs[ref_index].tp.showScan(
								scan, // CLTPass3d   scan,
								"test-strength");
					}
				} else {
					dsi[TwoQuadCLT.DSI_STRENGTH_AUX] =      aux_last_scan[1];
				}
				// perform photometric here, after first DSI
				
				boolean ran_photo_each = false;				
				if (photo_each && (!quadCLTs[ref_index].isPhotometricThis() || !batch_mode)) {
//					if (debugLevel > -3) {
//						System.out.println("**** Running photometric equalization for "+quadCLT_ref.getImageName()+
//								", current was from scene "+quadCLT_ref.getPhotometricScene()+" ****");
//					}
					OpticalFlow.photoEach(
							clt_parameters,      // CLTParameters                                 clt_parameters,
							colorProcParameters, // ColorProcParameters                           colorProcParameters,
							quadCLT_main,        // QuadCLT                                       quadCLT_main, // tiles should be set
							quadCLTs[ref_index],         // QuadCLT                                       quadCLT_ref, // tiles should be set
							dsi,                 // final double [][]                             dsi,
					    	clt_parameters.photo_num_full, // int      photo_num_full =              
							batch_mode,          // final boolean                                 batch_mode,
							MultiThreading.THREADS_MAX, // threadsMax,	
							updateStatus,        // final boolean                                 updateStatus,
							debugLevel);         // final int                                     debugLevel)			
					ran_photo_each = true; // will need to re-run after blue sky detection
					
				} else {
					System.out.println("(Re)using photometric calibration from this sequence reference "+quadCLTs[ref_index].getPhotometricScene());
					quadCLTs[ref_index].setQuadClt(); // just in case ?
				}
				quadCLTs[ref_index].setBlueSky  (
						sky_seed,           // double sky_seed, //  =       7.0;  // start with product of strength by diff_second below this
						lma_seed,
						sky_lim,            // double sky_lim, //   =      15.0; // then expand to product of strength by diff_second below this
						sky_shrink,         // int    sky_shrink, //  =       4;
						sky_expand_extra,   // int    sky_expand_extra, //  = 100; // 1?
						sky_bottleneck,     //int    sky_bottleneck, // 
						cold_scale, // =       0.2;  // <=1.0. 1.0 - disables temperature dependence
						cold_frac, // =        0.005; // this and lower will scale fom by  cold_scale
						hot_frac, // =         0.9;    // this and above will scale fom by 1.0
						min_strength, // =     0.08;
						seed_rows, // =        5; // sky should appear in this top rows 
						lowest_sky_row,        //  =   50;// appears that low - invalid, remove completely
						sky_bottom_override,   // double sky_temp_override,     // really cold average seed - ignore lowest_sky_row filter
						sky_override_shrink,   // int    shrink_for_temp,       // shrink before finding hottest sky
						sky_highest_min,       //  =   100; // lowest absolute value should not be higher (requires photometric)
						disp_boost_min,        // double disp_boost_min,    //  = 0.5;
						disp_boost_diff,       //double disp_boost_diff,   //  = 0.35;
						disp_boost_neibs,   //int    disp_boost_neibs,  //  = 2;
						disp_boost_amount,  //double disp_boost_amount, //  = 2.0;
						dsi[TwoQuadCLT.DSI_STRENGTH_AUX], // double [] strength,
						dsi[TwoQuadCLT.DSI_SPREAD_AUX], // double [] spread,
						dsi[TwoQuadCLT.DSI_DISPARITY_AUX_LMA], // double [] spread,
						dsi[TwoQuadCLT.DSI_AVGVAL_AUX],//	double [] avg_val,
						batch_mode? -1: 1); /// debugLevel);        // int debugLevel)
				if (ran_photo_each) {
					quadCLTs[ref_index].setBlueSky(null); // Reset blue sky - is it needed?
					// see if blue sky was detected - rerun photoEach
					boolean [] blue_sky = quadCLTs[ref_index].getBlueSky();
					boolean has_blue_sky = false;
					for (int i = 0; i < blue_sky.length; i++) if (blue_sky[i]) {
						has_blue_sky = true;
						break;
					}
					if (has_blue_sky) {
						if (debugLevel > -3) {
							System.out.println("Detected non-empty Blue Sky after initial DSI in "+quadCLTs[ref_index].getImageName()+
									", re-running photoEach()");
						}
						OpticalFlow.photoEach(
								clt_parameters,      // CLTParameters                                 clt_parameters,
								colorProcParameters, // ColorProcParameters                           colorProcParameters,
								quadCLT_main,        // QuadCLT                                       quadCLT_main, // tiles should be set
								quadCLTs[ref_index],         // QuadCLT                                       quadCLT_ref, // tiles should be set
								dsi,                 // final double [][]                             dsi,
								// just once?
						    	clt_parameters.photo_num_full, // int      photo_num_full =              
								batch_mode,          // final boolean                                 batch_mode,
								MultiThreading.THREADS_MAX, // threadsMax,	
								updateStatus,        // final boolean                                 updateStatus,
								debugLevel);         // final int                                     debugLevel)			

					}
				}
				quadCLTs[ref_index].saveDSIAll (
						"-DSI_MAIN", // String suffix, // "-DSI_MAIN"
						dsi);
				quadCLTs[ref_index].saveInterProperties( // save properties for interscene processing (extrinsics, ers, ...)
						null, // String path,             // full name with extension or w/o path to use x3d directory
						//							 	null, // Properties properties,   // if null - will only save extrinsics)
						debugLevel);
				// Copy source files to the model directory - they are needed for poses and interscene
				for (int scene_index =  ref_index - 1; scene_index >= 0 ; scene_index--) {
					TwoQuadCLT.copyJP4src( // actually there is no sense to process multiple image sets. Combine with other
							// processing?
							set_channels[scene_index].set_name, // String set_name
							quadCLT_main, // QuadCLT quadCLT_main,
							null, // QuadCLT quadCLT_aux,
							clt_parameters, // EyesisCorrectionParameters.DCTParameters dct_parameters,
							true, // boolean                                  skip_existing,					
							true, // false, // boolean                                  search_KML,
							debugLevel);
				} // split cycles to remove output clutter

				quadCLTs[ref_index].set_orient(0); // reset orientations
				quadCLTs[ref_index].set_accum(0);  // reset accumulations ("build_interscene") number
				earliest_scene = 0;  // reset failures, try to use again all scenes
			} else {// if (build_ref_dsi) {
				// need to read photometric from reference scene -INTERFRAME.corr-xml, and if it exists - set and propagate to main?
				if (photo_each  && !quadCLTs[ref_index].isPhotometricThis()) {
					if (debugLevel > -3) {
						System.out.println("Per-sequence photogrammetric calibration is required, but it does not exist");
					}
					if (photo_to_main) {
						quadCLT_main.setLwirOffsets(quadCLTs[ref_index].getLwirOffsets());
						quadCLT_main.setLwirScales (quadCLTs[ref_index].getLwirScales ());
						quadCLT_main.setLwirScales2(quadCLTs[ref_index].getLwirScales2());
						quadCLT_main.setPhotometricScene(quadCLTs[ref_index].getPhotometricScene());
						if (debugLevel > -3) {
							System.out.println("Propagating per-sequence photogrammetric calibration to main instance, will apply to next sequence and config file");
						}
					}
				} else {
					System.out.println("Without building reference DSI: (re)using photometric calibration from this sequence reference "+quadCLTs[ref_index].getPhotometricScene());
				}
				// read DSI_MAIN
				double [][] dsi = quadCLTs[ref_index].readDsiMain();
				if (dsi[TwoQuadCLT.DSI_SPREAD_AUX] == null) {
					System.out.println("DSI_MAIN file has old format and does not have spread data, will recalculate.");
				} else {
					quadCLTs[ref_index].setBlueSky  (
							sky_seed,                              // double sky_seed, //  =       7.0;  // start with product of strength by diff_second below this
							lma_seed,                              //          2.0;  // seed - disparity_lma limit
							sky_lim,                               // double sky_lim, //   =      15.0; // then expand to product of strength by diff_second below this
							sky_shrink,                            // int    sky_shrink, //  =       4;
							sky_expand_extra,                      // int    sky_expand_extra, //  = 100; // 1?
							sky_bottleneck,       //int    sky_bottleneck, // 
							cold_scale, // =       0.2;  // <=1.0. 1.0 - disables temperature dependence
							cold_frac, // =        0.005; // this and lower will scale fom by  cold_scale
							hot_frac, // =         0.9;    // this and above will scale fom by 1.0
							min_strength, // =     0.08;
							seed_rows, // =        5; // sky should appear in this top rows
							lowest_sky_row,        //  =     50;// appears that low - invalid, remove completely
							sky_bottom_override,   // double sky_temp_override,     // really cold average seed - ignore lowest_sky_row filter
							sky_override_shrink,   // int    shrink_for_temp,       // shrink before finding hottest sky
							sky_highest_min,       //  =   100; // lowest absolute value should not be higher (requires photometric)
							disp_boost_min,        // double disp_boost_min,    //  = 0.5;
							disp_boost_diff,       //double disp_boost_diff,   //  = 0.35;
							disp_boost_neibs,      //int    disp_boost_neibs,  //  = 2;
							disp_boost_amount,     //double disp_boost_amount, //  = 2.0;
							dsi[TwoQuadCLT.DSI_STRENGTH_AUX],      // double [] strength,
							dsi[TwoQuadCLT.DSI_SPREAD_AUX],        // double [] spread,
							dsi[TwoQuadCLT.DSI_DISPARITY_AUX_LMA], //double [] disp_lma,
							dsi[TwoQuadCLT.DSI_AVGVAL_AUX],//	double [] avg_val,
							debugLevel);        // int debugLevel)
				}
			}
		} // while (blue_sky == null) 
 /* */
    	
    	return null;
    }
}
