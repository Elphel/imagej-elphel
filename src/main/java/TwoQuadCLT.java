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
import java.util.ArrayList;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.io.FileSaver;

public class TwoQuadCLT {
	public long                                            startTime;     // start of batch processing
	public long                                            startSetTime;  // start of set processing
	public long                                            startStepTime; // start of step processing


	  public void processCLTQuadCorrs(
			  QuadCLT quadCLT_main,
			  QuadCLT quadCLT_aux,
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters   debayerParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
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
					  sourceFiles,                    // String []                                 sourceFiles,
					  set_channels_main[nSet].name(), // String                                    set_name,
					  referenceExposures_main,        // double []                                 referenceExposures,
					  channelFiles_main,              // int []                                    channelFiles,
					  scaleExposures_main,            //output  // double [] scaleExposures
					  saturation_imp_main,            //output  // boolean [][]                              saturation_imp,
					  debugLevel); // int                                       debugLevel);

			  ImagePlus [] imp_srcs_aux = quadCLT_aux.conditionImageSet(
					  clt_parameters,                 // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					  sourceFiles,                    // String []                                 sourceFiles,
					  set_channels_aux[nSet].name(), // String                                    set_name,
					  referenceExposures_aux,        // double []                                 referenceExposures,
					  channelFiles_aux,              // int []                                    channelFiles,
					  scaleExposures_aux,            //output  // double [] scaleExposures
					  saturation_imp_aux,            //output  // boolean [][]                              saturation_imp,
					  debugLevel); // int                                       debugLevel);

			  // Tempporarily processing individaully with the old code
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
					  colorProcParameters,
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
		  System.out.println("processCLTQuadCorrs(): processing "+(quadCLT_main.getTotalFiles(set_channels_main)+quadCLT_aux.getTotalFiles(set_channels_aux))+" files ("+set_channels_main.length+" file sets) finished at "+
				  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");

	  }


	  public void processCLTQuadCorrPairs(
			  QuadCLT quadCLT_main,
			  QuadCLT quadCLT_aux,
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters   debayerParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
//			  CorrectionColorProc.ColorGainsParameters       channelGainParameters_main,
//			  CorrectionColorProc.ColorGainsParameters       channelGainParameters_aux,
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
					  sourceFiles,                    // String []                                 sourceFiles,
					  set_channels_main[nSet].name(), // String                                    set_name,
					  referenceExposures_main,        // double []                                 referenceExposures,
					  channelFiles_main,              // int []                                    channelFiles,
					  scaleExposures_main,            //output  // double [] scaleExposures
					  saturation_imp_main,            //output  // boolean [][]                              saturation_imp,
					  debugLevel); // int                                       debugLevel);

			  ImagePlus [] imp_srcs_aux = quadCLT_aux.conditionImageSet(
					  clt_parameters,                 // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					  sourceFiles,                    // String []                                 sourceFiles,
					  set_channels_aux[nSet].name(), // String                                    set_name,
					  referenceExposures_aux,        // double []                                 referenceExposures,
					  channelFiles_aux,              // int []                                    channelFiles,
					  scaleExposures_aux,            //output  // double [] scaleExposures
					  saturation_imp_aux,            //output  // boolean [][]                              saturation_imp,
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
//					  channelGainParameters_main, // CorrectionColorProc.ColorGainsParameters       channelGainParameters_main,
//					  channelGainParameters_aux,  // CorrectionColorProc.ColorGainsParameters       channelGainParameters_aux,
					  rgbParameters,              // EyesisCorrectionParameters.RGBParameters       rgbParameters,
					  scaleExposures_main,        // double []	                                     scaleExposures_main, // probably not needed here - restores brightness of the final image
					  scaleExposures_aux,         // double []	                                     scaleExposures_aux, // probably not needed here - restores brightness of the final image
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
		  System.out.println("processCLTQuadCorrs(): processing "+(quadCLT_main.getTotalFiles(set_channels_main)+quadCLT_aux.getTotalFiles(set_channels_aux))+" files ("+set_channels_main.length+" file sets) finished at "+
				  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");

	  }




	  public ImagePlus [] processCLTQuadCorrPair(
			  QuadCLT                                        quadCLT_main,
			  QuadCLT                                        quadCLT_aux,
			  ImagePlus []                                   imp_quad_main,
			  ImagePlus []                                   imp_quad_aux,
			  boolean [][]                                   saturation_main, // (near) saturated pixels or null
			  boolean [][]                                   saturation_aux, // (near) saturated pixels or null
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters   debayerParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  EyesisCorrectionParameters.RGBParameters       rgbParameters,
			  double []	                                     scaleExposures_main, // probably not needed here - restores brightness of the final image
			  double []	                                     scaleExposures_aux, // probably not needed here - restores brightness of the final image
			  //			  final boolean    apply_corr, // calculate and apply additional fine geometry correction
			  //			  final boolean    infinity_corr, // calculate and apply geometry correction at infinity
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel){
		  //		  final boolean      batch_mode = clt_parameters.batch_run; //disable any debug images
		  boolean batch_mode = false;
		  boolean infinity_corr = false;
		  double [][] scaleExposures= {scaleExposures_main, scaleExposures_aux};
		  boolean toRGB=     quadCLT_main.correctionsParameters.toRGB;
		  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging? - TODO - move where it belongs
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

		  ImageDtt image_dtt = new ImageDtt();
		  double [][] ml_data = null;

		  final double [][][][][][][] clt_bidata = // new double[2][quad][nChn][tilesY][tilesX][][]; // first index - main/aux
				  image_dtt.clt_bi_quad (
						  clt_parameters,                       // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
						  clt_parameters.fat_zero,              // final double              fatzero,         // May use correlation fat zero from 2 different parameters - fat_zero and rig.ml_fatzero
						  tile_op_main,                         // final int [][]            tile_op_main,    // [tilesY][tilesX] - what to do - 0 - nothing for this tile
//						  tile_op_aux,                          // final int [][]            tile_op_aux,     // [tilesY][tilesX] - what to do - 0 - nothing for this tile
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
						  threadsMax,                           // final int                 threadsMax,  // maximal number of threads to launch
						  debugLevel);                          // final int                 globalDebugLevel);



		  double [][] texture_nonoverlap_main = null;
		  double [][] texture_nonoverlap_aux = null;
		  double [][] texture_overlap_main = null;
		  double [][] texture_overlap_aux = null;
		  String [] rgba_titles = {"red","blue","green","alpha"};
		  String [] rgba_weights_titles = {"red","blue","green","alpha","port0","port1","port2","port3","r-rms","b-rms","g-rms","w-rms"};
		  if ((texture_tiles_main != null) && (texture_tiles_aux != null)){
			  if (clt_parameters.show_nonoverlap){
				  texture_nonoverlap_main = image_dtt.combineRGBATiles(
						  texture_tiles_main,                 // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}
						  clt_parameters.transform_size,
						  false,                         // when false - output each tile as 16x16, true - overlap to make 8x8
						  clt_parameters.sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
						  threadsMax,                    // maximal number of threads to launch
						  debugLevel);
				  sdfa_instance.showArrays(
						  texture_nonoverlap_main,
						  tilesX * (2 * clt_parameters.transform_size),
						  tilesY * (2 * clt_parameters.transform_size),
						  true,
						  name + "-TXTNOL-D"+clt_parameters.disparity+"-MAIN",
						  (clt_parameters.keep_weights?rgba_weights_titles:rgba_titles));

				  texture_nonoverlap_aux = image_dtt.combineRGBATiles(
						  texture_tiles_aux,                 // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}
						  clt_parameters.transform_size,
						  false,                         // when false - output each tile as 16x16, true - overlap to make 8x8
						  clt_parameters.sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
						  threadsMax,                    // maximal number of threads to launch
						  debugLevel);
				  sdfa_instance.showArrays(
						  texture_nonoverlap_aux,
						  tilesX * (2 * clt_parameters.transform_size),
						  tilesY * (2 * clt_parameters.transform_size),
						  true,
						  name + "-TXTNOL-D"+clt_parameters.disparity+"-AUX",
						  (clt_parameters.keep_weights?rgba_weights_titles:rgba_titles));
			  }
			  if (!infinity_corr && (clt_parameters.show_overlap || clt_parameters.show_rgba_color)){
				  int alpha_index = 3;
				  texture_overlap_main = image_dtt.combineRGBATiles(
						  texture_tiles_main,                 // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}
						  clt_parameters.transform_size,
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

				  texture_overlap_aux = image_dtt.combineRGBATiles(
						  texture_tiles_aux,                 // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}
						  clt_parameters.transform_size,
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
							  tilesX * clt_parameters.transform_size,
							  tilesY * clt_parameters.transform_size,
							  true,
							  name + "-TXTOL-D"+clt_parameters.disparity+"-MAIN",
							  (clt_parameters.keep_weights?rgba_weights_titles:rgba_titles));
				  }
				  if (!batch_mode && clt_parameters.show_overlap) {
					  sdfa_instance.showArrays(
							  texture_overlap_aux,
							  tilesX * clt_parameters.transform_size,
							  tilesY * clt_parameters.transform_size,
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
							  tilesX *  clt_parameters.transform_size,
							  tilesY *  clt_parameters.transform_size,
							  1.0,         // double scaleExposure, // is it needed?
							  debugLevel );
					  ImagePlus imp_texture_aux = quadCLT_aux.linearStackToColor(
							  clt_parameters,
							  colorProcParameters,
							  rgbParameters,
							  name+"-texture", // String name,
							  "-D"+clt_parameters.disparity+"-AUX", //String suffix, // such as disparity=...
							  toRGB,
							  !quadCLT_aux.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
							  false, // true, // boolean saveShowIntermediate, // save/show if set globally
							  false, // true, // boolean saveShowFinal,        // save/show result (color image?)
							  ((clt_parameters.alpha1 > 0)? texture_rgba_aux: texture_rgb_aux),
							  tilesX *  clt_parameters.transform_size,
							  tilesY *  clt_parameters.transform_size,
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
							  quadCLT_aux,   //QuadCLT            quadCLT_aux,
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
		  }
		  //		  final double [][][][][][][] clt_bidata = // new double[2][quad][nChn][tilesY][tilesX][][]; // first index - main/aux
		  int quad_main = clt_bidata[0].length;
		  int quad_aux =  clt_bidata[1].length;

		  if (!infinity_corr && (clt_parameters.gen_chn_img || clt_parameters.gen_4_img || clt_parameters.gen_chn_stacks)) {
			  ImagePlus [] imps_RGB = new ImagePlus[quad_main+quad_aux];
			  for (int iQuadComb = 0; iQuadComb < imps_RGB.length; iQuadComb++){
				  int iAux = (iQuadComb >= quad_main) ? 1 : 0;
				  int iSubCam= iQuadComb - iAux * quad_main;

				  String title=name+"-"+String.format("%s%02d", ((iAux>0)?"A":"M"),iSubCam);

				  if (clt_parameters.corr_sigma > 0){ // no filter at all
					  for (int chn = 0; chn < clt_bidata[iAux][iSubCam].length; chn++) {
						  image_dtt.clt_lpf(
								  clt_parameters.corr_sigma,
								  clt_bidata[iAux][iSubCam][chn],
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
								  tilesX*clt_parameters.transform_size,
								  tilesY*clt_parameters.transform_size,
								  true,
								  results[iQuadComb].getTitle()+"-CLT-D"+clt_parameters.disparity);
					  }
				  }
				  double [][] iclt_data = new double [clt_bidata[iAux][iSubCam].length][];
				  for (int chn=0; chn<iclt_data.length;chn++){
					  iclt_data[chn] = image_dtt.iclt_2d(
							  clt_bidata[iAux][iSubCam][chn], // scanline representation of dcd data, organized as dct_size x dct_size tiles
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
						  tilesX *  clt_parameters.transform_size,
						  tilesY *  clt_parameters.transform_size,
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
			  }
		  }

		  return results;
	  }




//	  public double [][][][] getRigImageStacks(
	  public void getRigImageStacks(
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
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
					  null); // no margins, no oversample
		  }

		  double [][][] double_stacks_aux = new double [imp_quad_aux.length][][];
		  for (int i = 0; i < double_stacks_aux.length; i++){
			  double_stacks_aux[i] = quadCLT_aux.eyesisCorrections.bayerToDoubleStack(
					  imp_quad_aux[i], // source Bayer image, linearized, 32-bit (float))
					  null); // no margins, no oversample
		  }

		  for (int i = 0; i < double_stacks_main.length; i++){
			  for (int j =0 ; j < double_stacks_main[i][0].length; j++){
				  double_stacks_main[i][2][j]*=0.5; // Scale green 0.5 to compensate more pixels than R,B
			  }
		  }
		  for (int i = 0; i < double_stacks_aux.length; i++){
			  for (int j =0 ; j < double_stacks_aux[i][0].length; j++){
				  double_stacks_aux[i][2][j]*=0.5; // Scale green 0.5 to compensate more pixels than R,B
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


	  public void processInfinityRigs( // actually there is no sense to process multiple image sets. Combine with other processing?
			  QuadCLT quadCLT_main,
			  QuadCLT quadCLT_aux,
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
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
					  sourceFiles,                    // String []                                 sourceFiles,
					  set_channels_main[nSet].name(), // String                                    set_name,
					  referenceExposures_main,        // double []                                 referenceExposures,
					  channelFiles_main,              // int []                                    channelFiles,
					  scaleExposures_main,            //output  // double [] scaleExposures
					  saturation_imp_main,            //output  // boolean [][]                              saturation_imp,
					  debugLevel); // int                                       debugLevel);

			  ImagePlus [] imp_srcs_aux = quadCLT_aux.conditionImageSet(
					  clt_parameters,                 // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					  sourceFiles,                    // String []                                 sourceFiles,
					  set_channels_aux[nSet].name(), // String                                    set_name,
					  referenceExposures_aux,        // double []                                 referenceExposures,
					  channelFiles_aux,              // int []                                    channelFiles,
					  scaleExposures_aux,            //output  // double [] scaleExposures
					  saturation_imp_aux,            //output  // boolean [][]                              saturation_imp,
					  debugLevel); // int                                       debugLevel);

			  // Tempporarily processing individaully with the old code
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
		  System.out.println("processCLTQuadCorrs(): processing "+(quadCLT_main.getTotalFiles(set_channels_main)+quadCLT_aux.getTotalFiles(set_channels_aux))+" files ("+set_channels_main.length+" file sets) finished at "+
				  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");

	  }


	// use macro mode (strengths with limited disparity) to align at infinity if it is way off.
	  public boolean prealignInfinityRig(
			  QuadCLT                                        quadCLT_main,
			  QuadCLT                                        quadCLT_aux,
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
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
			  (new showDoubleFloatArrays()).showArrays(
					  dbg_macro,
					  tilesX,
					  tilesY,
					  true,
					  "MACRO-INPUT");
		  }

		  int macro_scale = clt_parameters.transform_size;
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
		  ImageDtt image_dtt = new ImageDtt();
		  image_dtt.clt_bi_macro(
				  clt_parameters,                       // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
				  clt_parameters.fat_zero,              // final double              fatzero,         // May use correlation fat zero from 2 different parameters - fat_zero and rig.ml_fatzero
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
		  (new showDoubleFloatArrays()).showArrays(
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
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel0){
		  final int        debugLevel = debugLevel0 + (clt_parameters.rig.rig_mode_debug?2:0);
//		  double [][][][] double_stacks =
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
		  final int tilesX = quadCLT_main.tp.getTilesX();

// perform full re-measure cycles
		  double [][] disparity_bimap = null;
          int [] num_new = new int[1];
		  for (int num_full_cycle = 0; num_full_cycle < clt_parameters.rig.rig_adjust_full_cycles;num_full_cycle++) {
			  disparity_bimap = null;
			  ArrayList<Integer> tile_list = new ArrayList<Integer>();
			  // measure and refine
			  for (int disp_step = 0; disp_step < clt_parameters.rig.rig_num_disp_steps; disp_step++) {
				  double disparity = 0.0;
				  if (clt_parameters.rig.rig_num_disp_steps > 1) {
					  disparity = disp_step * clt_parameters.rig.rig_disp_range/(clt_parameters.rig.rig_num_disp_steps -1);
				  }
				  disparity_bimap = measureNewRigDisparity(
						  quadCLT_main,      // QuadCLT             quadCLT_main,    // tiles should be set
						  quadCLT_aux,       // QuadCLT             quadCLT_aux,
						  disparity_bimap,   // double [][]         src_bimap,       // current state of measurements (or null for new measurement)
						  disparity,         // double              disparity,
						  tile_list,         // ArrayList<Integer>  tile_list,       // or null. If non-null - do not remeasure members of the list
						  clt_parameters,    // EyesisCorrectionParameters.CLTParameters       clt_parameters,
						  threadsMax,        // final int           threadsMax,      // maximal number of threads to launch
						  updateStatus,      // final boolean       updateStatus,
						  debugLevel);       // final int           debugLevel);

				  if (disparity_bimap != null){
					  if (clt_parameters.show_map &&  (debugLevel > 2) && clt_parameters.rig.rig_mode_debug){
						  (new showDoubleFloatArrays()).showArrays(
								  disparity_bimap,
								  tilesX,
								  disparity_bimap[0].length/tilesX,
								  true,
								  "DISP_MAP-D"+disparity,
								  ImageDtt.BIDISPARITY_TITLES);
					  }
				  }
				  if (disparity > 0.0) { // refine non-infinity passes
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
								  clt_parameters.rig.refine_min_strength , // double refine_min_strength, // do not refine weaker tiles
								  clt_parameters.rig.refine_tolerance ,    // double refine_tolerance,    // do not refine if absolute disparity below
								  null, // tile_list,       // ArrayList<Integer>             tile_list,       // or null
								  num_new,         // int     []                                     num_new,
								  clt_parameters,  // EyesisCorrectionParameters.CLTParameters clt_parameters,
								  threadsMax,      // final int                      threadsMax,      // maximal number of threads to launch
								  updateStatus,    // final boolean                  updateStatus,
								  debugLevel);     // final int                      debugLevel);
						  prev_bimap = disparity_bimap;
						  disparity_bimap = disparity_bimap_new;
						  if (disparity_bimap != null){
							  if (clt_parameters.show_map &&  (debugLevel > 2) && clt_parameters.rig.rig_mode_debug){
								  (new showDoubleFloatArrays()).showArrays(
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
								  clt_parameters.rig.refine_min_strength , // double refine_min_strength, // do not refine weaker tiles
								  clt_parameters.rig.refine_tolerance ,    // double refine_tolerance,    // do not refine if absolute disparity below
								  null, // tile_list,       // ArrayList<Integer>             tile_list,       // or null
								  num_new,         // int     []                                     num_new,
								  clt_parameters,  // EyesisCorrectionParameters.CLTParameters clt_parameters,
								  threadsMax,      // final int                      threadsMax,      // maximal number of threads to launch
								  updateStatus,    // final boolean                  updateStatus,
								  debugLevel);     // final int                      debugLevel);
						  prev_bimap = disparity_bimap;
						  disparity_bimap = disparity_bimap_new;
						  if (disparity_bimap != null){
							  if (clt_parameters.show_map &&  (debugLevel > 2) && clt_parameters.rig.rig_mode_debug){
								  (new showDoubleFloatArrays()).showArrays(
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
				  // refine for the existing list - all listed tiles, no thersholds
				  disparity_bimap =  refineRig(
						  quadCLT_main,    // QuadCLT                        quadCLT_main,    // tiles should be set
						  quadCLT_aux,     // QuadCLT                        quadCLT_aux,
						  disparity_bimap, // double [][]                    src_bimap,       // current state of measurements (or null for new measurement)
						  null,            // double [][]                                    prev_bimap, // previous state of measurements or null
						  null,       // double []                           scale_bad,
						  2,               // int                            refine_mode,     // 0 - by main, 1 - by aux, 2 - by inter
						  // will still re-measure infinity if refine_min_strength == 0.0
						  true,            // boolean                        keep_inf,        // keep expected disparity 0.0 if it was so
						  0.0,             // double refine_min_strength, // do not refine weaker tiles
						  0.0,             // double refine_tolerance,    // do not refine if absolute disparity below
						  tile_list,       // ArrayList<Integer>             tile_list,       // or null
						  num_new,         // int     []                                     num_new,
						  clt_parameters,  // EyesisCorrectionParameters.CLTParameters clt_parameters,
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
						  clt_parameters.rig.inf_weight , // double             infinity_importance, // of all measurements
						  clt_parameters.rig.rig_adjust_orientation,        // boolean            adjust_orientation,
						  clt_parameters.rig.rig_adjust_zoom,               // boolean            adjust_zoom,
						  clt_parameters.rig.rig_adjust_angle,              // boolean            adjust_angle,
						  clt_parameters.rig.rig_adjust_distance,           // boolean            adjust_distance,
						  clt_parameters.rig.rig_adjust_forward,            // boolean            adjust_forward, // not used
						  clt_parameters.rig.rig_correction_scale,          // double             scale_correction,
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
				  (new showDoubleFloatArrays()).showArrays(
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


	  //  improve DSI acquired for a single camera by use of a pair
	  // Run this after "CLT 3D"

	  public void enhanceByRig(
			  QuadCLT            quadCLT_main,  // tiles should be set
			  QuadCLT            quadCLT_aux,
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
			  final boolean                                  use_planes,
			  final int                                      threadsMax,  // maximal number of threads to launch
			  final boolean                                  updateStatus,
			  final int                                      debugLevel)// throws Exception
	  {
//		  boolean combine_oldsel_far = clt_parameters.rig.rf_master_infinity; // = true;
//		  boolean combine_oldsel_near = clt_parameters.rig.rf_master_near;    // = false; //
		  if ((quadCLT_main.tp == null) || (quadCLT_main.tp.clt_3d_passes == null)) {
			  String msg = "DSI data not available. Please run\"CLT 3D\" first";
			  IJ.showMessage("ERROR",msg);
			  System.out.println(msg);
			  return;
		  }

		  double [][] rig_disparity_strength =  quadCLT_main.getGroundTruthByRig();
		  if (rig_disparity_strength == null) {
			  if (use_planes) {
				  rig_disparity_strength = groundTruthByRigPlanes(
						  quadCLT_main,   // QuadCLT            quadCLT_main,  // tiles should be set
						  quadCLT_aux,    // QuadCLT            quadCLT_aux,
						  clt_parameters, // EyesisCorrectionParameters.CLTParameters       clt_parameters,
						  threadsMax,     // final int                                      threadsMax,  // maximal number of threads to launch
						  updateStatus,   // final boolean                                  updateStatus,
						  debugLevel);    // final int                                      debugLevel);
			  } else {
				  rig_disparity_strength = groundTruthByRig(
						  quadCLT_main,   // QuadCLT            quadCLT_main,  // tiles should be set
						  quadCLT_aux,    // QuadCLT            quadCLT_aux,
						  clt_parameters, // EyesisCorrectionParameters.CLTParameters       clt_parameters,
						  threadsMax,     // final int                                      threadsMax,  // maximal number of threads to launch
						  updateStatus,   // final boolean                                  updateStatus,
						  debugLevel);    // final int                                      debugLevel);
			  }
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

//		 if ((combine_oldsel_far || combine_oldsel_near) && (was_select != null) ) {
//			 for (int nTile=0; nTile < selection.length; nTile++) {
//				 selection[nTile] |= was_select[nTile] && (infinity_select[nTile]? combine_oldsel_far : combine_oldsel_near) ;
//			 }
//		 }
		 boolean [] selection_lone = biCamDSI.selectLoneFar(
				 clt_parameters.rig.ltfar_trusted_s, //  double     min_far_strength,
				 clt_parameters.rig.ltfar_trusted_d, // 	double     min_far_disparity,
				 infinity_select,                    // boolean [] infinity_select,
				 null,                               // boolean [] selection,
				 rig_disparity_strength[0],          // double []  disparity,
				 rig_disparity_strength[1]);         // double []  strength)
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
		 // restore master camera selections (poles may be remode by a filter
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
				 dbg_img[1][nTile] = (dbg_sel[2][nTile]? 1.0:0.0) + (dbg_sel[3][nTile]?1.0:0.0);
				 dbg_img[2][nTile] = (dbg_sel[4][nTile]? 1.0:0.0) + (dbg_sel[5][nTile]?1.0:0.0);
				 dbg_img[3][nTile] = (dbg_sel[6][nTile]? 1.0:0.0);
				 dbg_img[4][nTile] = (dbg_sel[7][nTile]? 1.0:0.0);
				 dbg_img[5][nTile] = (selection[nTile]? 1.0:0.0);
				 dbg_img[6][nTile] = (selection_lone[nTile]? 1.0:0.0);
			 }
			 String [] titles = {"old_sel","new_sel","combo-lone", "f-near", "f-far","final","lone"};
			  (new showDoubleFloatArrays()).showArrays(
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
// CLT ASSIGN needs best texture for each tile. Initially will just copy from teh previous master
// composite scan, later - fill disparity gaps and re-measure


/*
 *  TODO: interpolate disaprities before measuring to fill gaps?
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
		  ImageDtt image_dtt = new ImageDtt();
		  image_dtt.clt_bi_quad (
				  clt_parameters,                       // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
				  clt_parameters.fat_zero,              // final double              fatzero,         // May use correlation fat zero from 2 different parameters - fat_zero and rig.ml_fatzero
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
				  quadCLT_main.tp.getTilesX()*clt_parameters.transform_size, // final int                 width,

				  quadCLT_main.getGeometryCorrection(), // final GeometryCorrection  geometryCorrection_main,
				  quadCLT_aux.getGeometryCorrection(),  // final GeometryCorrection  geometryCorrection_aux,
				  quadCLT_main.getCLTKernels(),         // final double [][][][][][] clt_kernels_main, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				  quadCLT_aux.getCLTKernels(),          // final double [][][][][][] clt_kernels_aux,  // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				  clt_parameters.corr_magic_scale,      // final double              corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
				  false, // true,                                 // 	final boolean             keep_clt_data,
				  threadsMax,                           // final int                 threadsMax,  // maximal number of threads to launch
				  debugLevel-2);                        // final int                 globalDebugLevel);

		  return texture_tiles;

	  }

	  public void outputMLData(
			  QuadCLT                                  quadCLT_main,  // tiles should be set
			  QuadCLT                                  quadCLT_aux,
			  EyesisCorrectionParameters.CLTParameters clt_parameters,
			  final int                                threadsMax,  // maximal number of threads to launch
			  final boolean                            updateStatus,
			  final int                                debugLevel) // throws Exception
	  {
		  if ((quadCLT_main.tp == null) || (quadCLT_main.tp.clt_3d_passes == null)) {
			  String msg = "DSI data not available. Please run \"CLT 3D\" first";
			  IJ.showMessage("Error",msg);
			  System.out.println(msg);
			  return;
		  }
		  double [][] rig_disparity_strength =  quadCLT_main.getGroundTruthByRig();
		  if (rig_disparity_strength == null) {
			  rig_disparity_strength = groundTruthByRig(
					  quadCLT_main,   // QuadCLT            quadCLT_main,  // tiles should be set
					  quadCLT_aux,    // QuadCLT            quadCLT_aux,
					  clt_parameters, // EyesisCorrectionParameters.CLTParameters       clt_parameters,
					  threadsMax,     // final int                                      threadsMax,  // maximal number of threads to launch
					  updateStatus,   // final boolean                                  updateStatus,
					  debugLevel-2);  // final int                                      debugLevel);
			  if (rig_disparity_strength != null) {
				  quadCLT_main.tp.rig_disparity_strength = rig_disparity_strength;
			  } else {
				  System.out.println("outputMLData(): failed to get ground truth data, aborting");
				  return;
			  }
		  }

		  String ml_directory= quadCLT_main.correctionsParameters.selectMlDirectory(
		  true,  // smart,
		  true);  //newAllowed, // save
			Correlation2d corr2d = new Correlation2d(
					clt_parameters.img_dtt,              // ImageDttParameters  imgdtt_params,
					clt_parameters.transform_size,             // int transform_size,
					2.0,                        //  double wndx_scale, // (wndy scale is always 1.0)
					(debugLevel > -1));   //   boolean debug)

		  for (int sweep_step = 0; sweep_step < clt_parameters.rig.ml_sweep_steps; sweep_step++){
			  double disparity_offset = 0; // clt_parameters.rig.ml_disparity_sweep * (2.0 * sweep_step/(clt_parameters.rig.ml_sweep_steps - 1.0) -1.0);
			  if (clt_parameters.rig.ml_sweep_steps > 1) {
				  disparity_offset = clt_parameters.rig.ml_disparity_sweep * (2.0 * sweep_step/(clt_parameters.rig.ml_sweep_steps - 1.0) -1.0);
			  }
	    	  double [][] ml_data = remeasureRigML(
	    			  disparity_offset,                         // double               disparity_offset,
	    			  quadCLT_main,                             // QuadCLT              quadCLT_main,    // tiles should be set
	    			  quadCLT_aux,                              // QuadCLT              quadCLT_aux,
//	    			  disparity_bimap,                          // double [][]          src_bimap,
	    			  rig_disparity_strength[0],                // double []            disparity,
	    			  rig_disparity_strength[1],                // double []            strength,

	    			  clt_parameters,                           // EyesisCorrectionParameters.CLTParameters clt_parameters,
	    			  clt_parameters.rig.ml_hwidth,             // int ml_hwidth
	    			  clt_parameters.rig.ml_fatzero,            // double               fatzero,
	    			  threadsMax,                               // final int            threadsMax,  // maximal number of threads to launch
	    			  updateStatus,                             // final boolean        updateStatus,
	    			  debugLevel);                              // final int            debugLevel);
	    	  saveMlFile(
	    			  quadCLT_main.image_name+"-ML_DATA-", // String               ml_title,
	    			  ml_directory,                             // String               ml_directory,
	    			  disparity_offset,                         // double               disp_offset,
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
		  }
	  }

	  public double [][] prepareRefineExistingDSI(
			  QuadCLT            quadCLT_main,  // tiles should be set
			  QuadCLT            quadCLT_aux,
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
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
				  threadsMax,        // final int           threadsMax,      // maximal number of threads to launch
				  updateStatus,      // final boolean       updateStatus,
				  debugLevel);       // final int           debugLevel);


		  int [] num_new = new int[1];
          boolean [] trusted_infinity = 	  getTrustedDisparity(
        		  quadCLT_main,                            // QuadCLT            quadCLT_main,  // tiles should be set
        		  quadCLT_aux,                             // QuadCLT            quadCLT_aux,
        		  clt_parameters.rig.min_trusted_strength, // double             min_combo_strength,    // check correlation strength combined for all 3 correlations
        		  clt_parameters.grow_disp_trust,          // double             max_trusted_disparity, // 4.0 -> change to rig_trust
        		  clt_parameters.rig.trusted_tolerance,    // double             trusted_tolerance,
    			  null,                                    // boolean []         was_trusted,
        		  disparity_bimap_infinity );              // double [][]        bimap // current state of measurements

          if (clt_parameters.show_map &&  (debugLevel > 0) && clt_parameters.rig.rig_mode_debug){
        	  (new showDoubleFloatArrays()).showArrays(
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
        	  double [][] disparity_bimap_new =  refineRigSel(
        			  quadCLT_main,    // QuadCLT                        quadCLT_main,    // tiles should be set
        			  quadCLT_aux,     // QuadCLT                        quadCLT_aux,
        			  disparity_bimap_infinity, // double [][]                    src_bimap,       // current state of measurements (or null for new measurement)
        			  prev_bimap,      // double [][]                    prev_bimap, // previous state of measurements or null
					  scale_bad,       // double []                      scale_bad,
					  refine_inter,    // int                            refine_mode,     // 0 - by main, 1 - by aux, 2 - by inter
        			  false,           // boolean                        keep_inf,        // keep expected disparity 0.0 if it was so
        			  0.0, // clt_parameters.rig.refine_min_strength , // double refine_min_strength, // do not refine weaker tiles
        			  clt_parameters.rig.refine_tolerance ,    // double refine_tolerance,    // do not refine if absolute disparity below
        			  trusted_infinity, // tile_list,       // ArrayList<Integer>             tile_list,       // or null
        			  num_new,         // int     []                                     num_new,
        			  clt_parameters,  // EyesisCorrectionParameters.CLTParameters clt_parameters,
        			  threadsMax,      // final int                      threadsMax,      // maximal number of threads to launch
        			  updateStatus,    // final boolean                  updateStatus,
        			  debugLevel);     // final int                      debugLevel);
        	  prev_bimap = disparity_bimap_infinity;
        	  disparity_bimap_infinity = disparity_bimap_new;
              trusted_infinity = 	  getTrustedDisparity(
            		  quadCLT_main,                             // QuadCLT            quadCLT_main,  // tiles should be set
            		  quadCLT_aux,                              // QuadCLT            quadCLT_aux,
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
            	  (new showDoubleFloatArrays()).showArrays(
            			  scale_bad,
            			  tilesX,
            			  disparity_bimap_infinity[0].length/tilesX,
            			  quadCLT_main.image_name+"-INFINITY-SCALE_BAD"+clt_parameters.disparity);

        	  }
        	  (new showDoubleFloatArrays()).showArrays(
        			  disparity_bimap_infinity,
        			  tilesX,
        			  disparity_bimap_infinity[0].length/tilesX,
        			  true,
        			  quadCLT_main.image_name+"-INFINITY-REFINED-TRUSTED"+clt_parameters.disparity,
        			  ImageDtt.BIDISPARITY_TITLES);
          }

          // Get DSI from the main camera
		  quadCLT_main.tp.trimCLTPasses(false); // remove rig composite scan if any
		  // last but not including any rid data
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
        		  0.5*clt_parameters.rig.min_trusted_strength,       // double             min_combo_strength,    // check correlation strength combined for all 3 correlations
        		  clt_parameters.grow_disp_trust,                    // double             max_trusted_disparity, // 4.0 -> change to rig_trust
        		  clt_parameters.rig.trusted_tolerance,              // double             trusted_tolerance,
    			  null,                                              // boolean []         was_trusted,
        		  disparity_bimap); // double [][]        bimap // current state of measurements
          if (clt_parameters.show_map &&  (debugLevel > 0) && clt_parameters.rig.rig_mode_debug){
        	  (new showDoubleFloatArrays()).showArrays(
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
        			  0.0, // clt_parameters.rig.refine_min_strength , // double refine_min_strength, // do not refine weaker tiles
        			  clt_parameters.rig.refine_tolerance ,    // double refine_tolerance,    // do not refine if absolute disparity below
        			  trusted_near,    // tile_list,       // ArrayList<Integer>             tile_list,       // or null
        			  num_new,         // int     []                                     num_new,
        			  clt_parameters,  // EyesisCorrectionParameters.CLTParameters clt_parameters,
        			  threadsMax,      // final int                      threadsMax,      // maximal number of threads to launch
        			  updateStatus,    // final boolean                  updateStatus,
        			  debugLevel);     // final int                      debugLevel);
        	  prev_bimap = disparity_bimap;
        	  disparity_bimap = disparity_bimap_new;
              trusted_near = 	  getTrustedDisparity(
            		  quadCLT_main,                            // QuadCLT            quadCLT_main,  // tiles should be set
            		  quadCLT_aux,                             // QuadCLT            quadCLT_aux,
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
            	  (new showDoubleFloatArrays()).showArrays(
            			  scale_bad,
            			  tilesX,
            			  disparity_bimap[0].length/tilesX,
            			  quadCLT_main.image_name+"-NEAR-SCALE_BAD"+clt_parameters.disparity);

        	  }
        	  (new showDoubleFloatArrays()).showArrays(
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

          if (clt_parameters.show_map &&  (debugLevel > 0) && clt_parameters.rig.rig_mode_debug){
        	  (new showDoubleFloatArrays()).showArrays(
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
				  (new showDoubleFloatArrays()).showArrays(
						  disparity_bimap,
						  tilesX,
						  disparity_bimap[0].length/tilesX,
						  true,
						  quadCLT_main.image_name+"-DSI-ALL-TRUSTED"+clt_parameters.disparity,
						  ImageDtt.BIDISPARITY_TITLES);
          }

          return disparity_bimap;

	  }

	  public double [][] groundTruthByRigPlanes(
			  QuadCLT            quadCLT_main,  // tiles should be set
			  QuadCLT            quadCLT_aux,
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
			  final int                                      threadsMax,  // maximal number of threads to launch
			  final boolean                                  updateStatus,
			  final int                                      debugLevel) //  throws Exception
	  {
		  final int num_full_cycles =       3;  // Number of full low-texture cycles that include growing flat LT and trimmin weak FG over BG
		  final int num_cross_gaps_cycles = 20; // maximalnumger of adding new tiles cycles while "crossing the gaps)
		  final int min_cross_gaps_new =    20; // minimal number of the new added tiles
		  final int refine_inter = 2; // 3; // 3 - dx, 2 - disparity
		  final int tilesX = quadCLT_main.tp.getTilesX();
		  final int tilesY = quadCLT_main.tp.getTilesY();
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

		  biCamDSI.addBiScan(disparity_bimap);
		  if (clt_parameters.show_map &&  (debugLevel > 0) && clt_parameters.rig.rig_mode_debug){
			  biCamDSI.getLastBiScan().showScan(
					  quadCLT_main.image_name+"-BISCAN_initial");
		  }
		  int [][] dxy = {{0, -1},{0,1},{-1,0},{1,0}};
		  for (int num_fcycle = 0; num_fcycle < num_full_cycles; num_fcycle++) {
			  // Grow tiles, cross gaps (do not trim yet
			  //
			  for (int num_cycle = 0; num_cycle < num_cross_gaps_cycles; num_cycle++) {
				  BiScan last_scan = biCamDSI.getLastBiScan();
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
					  biCamDSI.getLastBiScan().showScan(
							  quadCLT_main.image_name+"-BISCAN_TRUSTED"+num_fcycle+"-"+num_cycle);
				  }

				  /*
				   * @return array of 3 numbers: number of trusted strong tiles, number of additional trusted by plane fitting, and number of all
				   * somewhat strong tiles
				   */
				  if (debugLevel > -2) {
					  System.out.println("groundTruthByRigPlanes() grow pass "+num_cycle+" of "+ num_cross_gaps_cycles+
							  " strong trusted: "+trusted_stats[0]+ " neib trusted: "+trusted_stats[1]+" weak trusted: " + trusted_stats[2]);
				  }
				  int num_added_tiles =0;
				  if (num_cycle < 2*dxy.length) {
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
				  } else {

					  // suggest new disparities, using plane surfaces (extending around that may cause false surfaces)
					  /*
					   * Suggest disparities to try for the tiles in poorly textured areas by fitting planes in DSI
					   * calcTrusted should be called before to set up trusted/cond_trusted tiles
					   * suggested tiles will be compared against and made sure they differ by more than a specified margin
					   * 1) current measured (refined) disparity value
					   * 2) target disaprity that lead to the current measurement after refinement
					   * 3) any other disable measurement
					   * 4) any target disparity that lead to the disabled measurement
					   * @return number of new tiles to measure in the  array of suggested disparities - Double.NaN - nothing suggested
					   *  for the tile. May need additional filtering to avoid suggested already tried disparities
					   */

					  num_added_tiles = last_scan.suggestNewScan(
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
							  clt_parameters.tileX,                   // final int        dbg_x,
							  clt_parameters.tileY,                   // final int        dbg_y,
							  debugLevel);                            // final int        debugLevel);
				  }
				  if (debugLevel > -2) {
					  System.out.println("groundTruthByRigPlanes() full cycle = "+num_fcycle+", grow pass "+num_cycle+" of "+ num_cross_gaps_cycles+
							  " suggestNewScan() -> "+num_added_tiles);
				  }
				  //num_cycle < num_cross_gaps_cycles;
				  boolean last_cycle = (num_added_tiles < min_cross_gaps_new) || (num_cycle >= (num_cross_gaps_cycles-1));
				  if (clt_parameters.show_map &&  (debugLevel > -2) && clt_parameters.rig.rig_mode_debug){
					 if (last_cycle || (num_cycle < 2 * dxy.length))
						 biCamDSI.getLastBiScan().showScan(
							  quadCLT_main.image_name+"-BISCAN_SUGGESTED"+num_fcycle+"-"+num_cycle);
				  }

				  if (last_cycle && clt_parameters.rig.pf_en_trim_fg) { // last cycle and trimming enabled
					  if (debugLevel > -2) {
						  //						  System.out.println("groundTruthByRigPlanes(): num_added_tiles= "+num_added_tiles+" > "+min_cross_gaps_new+", done growing over gaps");
						  System.out.println("groundTruthByRigPlanes(): that was the last growing over gaps cycle, performing trimming hanging weak FG over BG");
					  }

					  /*
					   * Disable low-textured tiles are not between strong tiles, but on one side of it.
					   * This method relies on the assumption that FG edge should bave strong correlation, so it tries multiple directions
					   * from the weak (not trusted strong) tiles and trims tiles that eithre do not have anything in that direction or have
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
					  if (debugLevel > -2) {
						  System.out.println("groundTruthByRigPlanes(): full cycle="+num_fcycle+" num_trimmed= "+num_trimmed+" tiles");
					  }
					  if (clt_parameters.show_map &&  (debugLevel > 0) && clt_parameters.rig.rig_mode_debug){
						  biCamDSI.getLastBiScan().showScan(
								  quadCLT_main.image_name+"-BISCAN_TRIMMED"+num_fcycle+"-"+num_cycle);
					  }

					  // suggest again, after trimming
					  int num_added_tiles_trimmed = last_scan.suggestNewScan(
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
							  clt_parameters.tileX,                   // final int        dbg_x,
							  clt_parameters.tileY,                   // final int        dbg_y,
							  debugLevel);                            // final int        debugLevel);
					  if (debugLevel > -2) {
						  System.out.println("groundTruthByRigPlanes() full cycle = "+num_fcycle+", grow pass "+num_cycle+" of "+ num_cross_gaps_cycles+
								  " suggestNewScan() -> "+num_added_tiles_trimmed+"( after trimming)");
					  }
					  if (clt_parameters.show_map &&  (debugLevel > -2) && clt_parameters.rig.rig_mode_debug){
						  biCamDSI.getLastBiScan().showScan(
								  quadCLT_main.image_name+"-BISCAN_TRIMMED_SUGGESTED"+num_fcycle+"-"+num_cycle);
					  }

					  //					  break; // too few added before trimmimng or number of steps exceeded limit
				  } // if (last_cycle && clt_parameters.rig.pf_en_trim_fg) { // last cycle and trimming enabled
				  // measure and refine
				  double [] target_disparity = biCamDSI.getTargetDisparity(-1); // get last

				  // Measure provided tiles (brdeak after, if it was the last cycle)

				  disparity_bimap = measureNewRigDisparity(
						  quadCLT_main,      // QuadCLT             quadCLT_main,    // tiles should be set
						  quadCLT_aux,       // QuadCLT             quadCLT_aux,
						  target_disparity,  // double []                                      disparity, // Double.NaN - skip, ohers - measure
						  clt_parameters,    // EyesisCorrectionParameters.CLTParameters       clt_parameters,
						  threadsMax,        // final int           threadsMax,      // maximal number of threads to launch
						  updateStatus,      // final boolean       updateStatus,
						  debugLevel);       // final int           debugLevel);

				  // refine measurements
				  int [] num_new = new int[1];
				  boolean [] trusted_measurements = 	  getTrustedDisparity(
						  quadCLT_main,                            // QuadCLT            quadCLT_main,  // tiles should be set
						  quadCLT_aux,                             // QuadCLT            quadCLT_aux,
						  clt_parameters.rig.min_trusted_strength, // double             min_combo_strength,    // check correlation strength combined for all 3 correlations
						  clt_parameters.grow_disp_trust,          // double             max_trusted_disparity, // 4.0 -> change to rig_trust
						  clt_parameters.rig.trusted_tolerance,    // double             trusted_tolerance,
						  null,                                    // boolean []         was_trusted,
						  disparity_bimap);              // double [][]        bimap // current state of measurements

				  if (clt_parameters.show_map &&  (debugLevel > 0) && clt_parameters.rig.rig_mode_debug){
					  (new showDoubleFloatArrays()).showArrays(
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
							  0.0, // clt_parameters.rig.refine_min_strength , // double refine_min_strength, // do not refine weaker tiles
							  clt_parameters.rig.refine_tolerance ,    // double refine_tolerance,    // do not refine if absolute disparity below
							  trusted_measurements, // tile_list,       // ArrayList<Integer>             tile_list,       // or null
							  num_new,         // int     []                                     num_new,
							  clt_parameters,  // EyesisCorrectionParameters.CLTParameters clt_parameters,
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
					  if (num_new[0] < clt_parameters.rig.pf_min_new) break;
				  }


				  //FIXME: 	show only for 	last_cycle
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
						  (new showDoubleFloatArrays()).showArrays(
								  scale_bad,
								  tilesX,
								  tilesY,
								  quadCLT_main.image_name+"-gaps_cycle"+num_cycle+"-scale_bad");

					  }
					  (new showDoubleFloatArrays()).showArrays(
							  dbg_img,
							  tilesX,
							  tilesY,
							  true,
							  quadCLT_main.image_name+"-gaps_cycle"+num_cycle+"-REFINED_TRUSTED",
							  ImageDtt.BIDISPARITY_TITLES);
				  }

				  // add refined data
				  biCamDSI.addBiScan(disparity_bimap);
				  // find strongest
				  biCamDSI.getLastBiScan().copyLastStrongestEnabled(
						  clt_parameters.rig.pf_last_priority); // final boolean last_priority)

				  double afloor = clt_parameters.rig.pf_trusted_strength * clt_parameters.rig.pf_strength_rfloor;

				  // replace strongest by fittest
				  for (int nfit = 0; nfit < 5; nfit++) {
					  int num_replaced = biCamDSI.getLastBiScan().copyFittestEnabled(
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

				  if (clt_parameters.show_map &&  (debugLevel > 0) && clt_parameters.rig.rig_mode_debug){
					  biCamDSI.getLastBiScan().showScan(
							  quadCLT_main.image_name+"-BISCAN_"+num_fcycle+"-"+num_cycle);
				  }
				  if (last_cycle) { // last cycle
					  if ((debugLevel > -2) && clt_parameters.rig.rig_mode_debug){
						  System.out.println("groundTruthByRigPlanes(): that was refinement measurements after the trimming of hanging weak FG over BG");
					  }
					  break;
				  }
//	public void showScan(String title) {
			  } // for (int num_cycle = 0; num_cycle < num_cross_gaps_cycles; num_cycle++) {
		  }// 		  for (int num_fcycle = 0; num_fcycle < num_full_cycles; num_fcycle++) {


		  double [][] rig_disparity_strength = biCamDSI.getLastBiScan().getDisparityStrength(
				    		false, // boolean only_strong,
				    		false, // boolean only_trusted,
				    		true); // boolean only_enabled,

		  return rig_disparity_strength;
	  }




	  public double [][] groundTruthByRig(
			  QuadCLT            quadCLT_main,  // tiles should be set
			  QuadCLT            quadCLT_aux,
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
			  final int                                      threadsMax,  // maximal number of threads to launch
			  final boolean                                  updateStatus,
			  final int                                      debugLevel) //  throws Exception
	  {
		  int refine_inter = 2; // 3; // 3 - dx, 2 - disparity
		  System.out.println("enhanceByRig()");
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
		  // See if auxiliary camera has images configured, if not - do it now. First need to get image name
		  if (quadCLT_aux.image_name ==null) {
//			  String image_name = quadCLT_main.image_name;
			  String [] sourceFiles=quadCLT_main.correctionsParameters.getSourcePaths();
			  QuadCLT.SetChannels [] set_channels_aux =  quadCLT_aux.setChannels(quadCLT_main.image_name, debugLevel);
			  if ((set_channels_aux == null) || (set_channels_aux.length==0)) {
				  System.out.println("No files for the auxiliary camera match series "+quadCLT_main.image_name);
				  return null;
			  }
			  double [] referenceExposures_aux =  quadCLT_aux.eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
			  int [] channelFiles_aux =  set_channels_aux[0].fileNumber();
			  // make single
			  boolean [][] saturation_aux =  (clt_parameters.sat_level > 0.0)? new boolean[channelFiles_aux.length][] : null;
			  double [] scaleExposures_aux =  new double[channelFiles_aux.length];
			  ImagePlus [] imp_srcs_aux = quadCLT_aux.conditionImageSet(
					  clt_parameters,               // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					  sourceFiles,                  // String []                                 sourceFiles,
					  quadCLT_main.image_name,      // set_channels_aux[0].name(), // String                                    set_name,
					  referenceExposures_aux,       // double []                                 referenceExposures,
					  channelFiles_aux,             // int []                                    channelFiles,
					  scaleExposures_aux,           //output  // double [] scaleExposures
					  saturation_aux,               //output  // boolean [][]                              saturation_imp,
					  debugLevel); // int                                       debugLevel);

			  double [][][] double_stacks_aux = new double [imp_srcs_aux.length][][];
			  for (int i = 0; i < double_stacks_aux.length; i++){
				  double_stacks_aux[i] = quadCLT_aux.eyesisCorrections.bayerToDoubleStack(
						  imp_srcs_aux[i], // source Bayer image, linearized, 32-bit (float))
						  null); // no margins, no oversample
			  }

			  for (int i = 0; i < double_stacks_aux.length; i++){
				  for (int j =0 ; j < double_stacks_aux[i][0].length; j++){
					  double_stacks_aux[i][2][j]*=0.5; // Scale green 0.5 to compensate more pixels than R,B
				  }
			  }
			  quadCLT_aux.setTiles (imp_srcs_aux[0], // set global tp.tilesX, tp.tilesY
					  clt_parameters,
					  threadsMax);
			  quadCLT_aux.image_name =      quadCLT_main.image_name;
			  quadCLT_aux.image_data =      double_stacks_aux;
			  quadCLT_aux.saturation_imp =  saturation_aux;
			  quadCLT_aux.tp.setTrustedCorrelation(clt_parameters.grow_disp_trust);
			  quadCLT_aux.tp.resetCLTPasses();
		  }
		  //Re-measure background
		  final int tilesX = quadCLT_main.tp.getTilesX();
		  final int tilesY = quadCLT_main.tp.getTilesY();
		  double [][] disparity_bimap_infinity = measureNewRigDisparity(
				  quadCLT_main,      // QuadCLT             quadCLT_main,    // tiles should be set
				  quadCLT_aux,       // QuadCLT             quadCLT_aux,
				  null,              // double [][]         src_bimap,       // current state of measurements (or null for new measurement)
				  0,                 // double              disparity,
				  null,              // ArrayList<Integer>  tile_list,       // or null. If non-null - do not remeasure members of the list
				  clt_parameters,    // EyesisCorrectionParameters.CLTParameters       clt_parameters,
				  threadsMax,        // final int           threadsMax,      // maximal number of threads to launch
				  updateStatus,      // final boolean       updateStatus,
				  debugLevel);       // final int           debugLevel);


		  int [] num_new = new int[1];
          boolean [] trusted_infinity = 	  getTrustedDisparity(
        		  quadCLT_main,                            // QuadCLT            quadCLT_main,  // tiles should be set
        		  quadCLT_aux,                             // QuadCLT            quadCLT_aux,
        		  clt_parameters.rig.min_trusted_strength, // double             min_combo_strength,    // check correlation strength combined for all 3 correlations
        		  clt_parameters.grow_disp_trust,          // double             max_trusted_disparity, // 4.0 -> change to rig_trust
        		  clt_parameters.rig.trusted_tolerance,    // double             trusted_tolerance,
    			  null,                                    // boolean []         was_trusted,
        		  disparity_bimap_infinity );              // double [][]        bimap // current state of measurements

          if (clt_parameters.show_map &&  (debugLevel > 0) && clt_parameters.rig.rig_mode_debug){
        	  (new showDoubleFloatArrays()).showArrays(
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
        	  double [][] disparity_bimap_new =  refineRigSel(
        			  quadCLT_main,    // QuadCLT                        quadCLT_main,    // tiles should be set
        			  quadCLT_aux,     // QuadCLT                        quadCLT_aux,
        			  disparity_bimap_infinity, // double [][]                    src_bimap,       // current state of measurements (or null for new measurement)
        			  prev_bimap,      // double [][]                    prev_bimap, // previous state of measurements or null
					  scale_bad,       // double []                      scale_bad,
					  refine_inter,    // int                            refine_mode,     // 0 - by main, 1 - by aux, 2 - by inter
        			  false,           // boolean                        keep_inf,        // keep expected disparity 0.0 if it was so
        			  0.0, // clt_parameters.rig.refine_min_strength , // double refine_min_strength, // do not refine weaker tiles
        			  clt_parameters.rig.refine_tolerance ,    // double refine_tolerance,    // do not refine if absolute disparity below
        			  trusted_infinity, // tile_list,       // ArrayList<Integer>             tile_list,       // or null
        			  num_new,         // int     []                                     num_new,
        			  clt_parameters,  // EyesisCorrectionParameters.CLTParameters clt_parameters,
        			  threadsMax,      // final int                      threadsMax,      // maximal number of threads to launch
        			  updateStatus,    // final boolean                  updateStatus,
        			  debugLevel);     // final int                      debugLevel);
        	  prev_bimap = disparity_bimap_infinity;
        	  disparity_bimap_infinity = disparity_bimap_new;
              trusted_infinity = 	  getTrustedDisparity(
            		  quadCLT_main,                             // QuadCLT            quadCLT_main,  // tiles should be set
            		  quadCLT_aux,                              // QuadCLT            quadCLT_aux,
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
            	  (new showDoubleFloatArrays()).showArrays(
            			  scale_bad,
            			  tilesX,
            			  disparity_bimap_infinity[0].length/tilesX,
            			  quadCLT_main.image_name+"-INFINITY-SCALE_BAD"+clt_parameters.disparity);

        	  }
        	  (new showDoubleFloatArrays()).showArrays(
        			  disparity_bimap_infinity,
        			  tilesX,
        			  disparity_bimap_infinity[0].length/tilesX,
        			  true,
        			  quadCLT_main.image_name+"-INFINITY-REFINED-TRUSTED"+clt_parameters.disparity,
        			  ImageDtt.BIDISPARITY_TITLES);
          }

          // Get DSI from the main camera
		  quadCLT_main.tp.trimCLTPasses(false); // remove rig composite scan if any
		  // last but not including any rid data
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
        		  0.5*clt_parameters.rig.min_trusted_strength,       // double             min_combo_strength,    // check correlation strength combined for all 3 correlations
        		  clt_parameters.grow_disp_trust,                    // double             max_trusted_disparity, // 4.0 -> change to rig_trust
        		  clt_parameters.rig.trusted_tolerance,              // double             trusted_tolerance,
    			  null,                                              // boolean []         was_trusted,
        		  disparity_bimap); // double [][]        bimap // current state of measurements
          if (clt_parameters.show_map &&  (debugLevel > 0) && clt_parameters.rig.rig_mode_debug){
        	  (new showDoubleFloatArrays()).showArrays(
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
        			  0.0, // clt_parameters.rig.refine_min_strength , // double refine_min_strength, // do not refine weaker tiles
        			  clt_parameters.rig.refine_tolerance ,    // double refine_tolerance,    // do not refine if absolute disparity below
        			  trusted_near,    // tile_list,       // ArrayList<Integer>             tile_list,       // or null
        			  num_new,         // int     []                                     num_new,
        			  clt_parameters,  // EyesisCorrectionParameters.CLTParameters clt_parameters,
        			  threadsMax,      // final int                      threadsMax,      // maximal number of threads to launch
        			  updateStatus,    // final boolean                  updateStatus,
        			  debugLevel);     // final int                      debugLevel);
        	  prev_bimap = disparity_bimap;
        	  disparity_bimap = disparity_bimap_new;
              trusted_near = 	  getTrustedDisparity(
            		  quadCLT_main,                            // QuadCLT            quadCLT_main,  // tiles should be set
            		  quadCLT_aux,                             // QuadCLT            quadCLT_aux,
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
            	  (new showDoubleFloatArrays()).showArrays(
            			  scale_bad,
            			  tilesX,
            			  disparity_bimap[0].length/tilesX,
            			  quadCLT_main.image_name+"-NEAR-SCALE_BAD"+clt_parameters.disparity);

        	  }
        	  (new showDoubleFloatArrays()).showArrays(
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

          if (clt_parameters.show_map &&  (debugLevel > 0) && clt_parameters.rig.rig_mode_debug){
        	  (new showDoubleFloatArrays()).showArrays(
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
				  (new showDoubleFloatArrays()).showArrays(
						  disparity_bimap,
						  tilesX,
						  disparity_bimap[0].length/tilesX,
						  true,
						  quadCLT_main.image_name+"-DSI-ALL-TRUSTED"+clt_parameters.disparity,
						  ImageDtt.BIDISPARITY_TITLES);
          }


// repeat multiple times (need to compare added tiles (after cleanup? - return number of tiles after cleanu-up and compare improvements)
          int num_trusted = 0;
//		  final int tilesX = quadCLT_main.tp.getTilesX();
//		  final int tilesY = quadCLT_main.tp.getTilesY();
          // grow around using all camera and inter-camera correlations (try to get low-textured,like our street pavement)
          BiCamDSI biCamDSI = new BiCamDSI( tilesX, tilesY,threadsMax);
          int min_added_tiles = clt_parameters.rig.lt_min_new;
          for (int num_fill = 0; num_fill < clt_parameters.rig.lt_repeat; num_fill++) {
        	  int num_new_trusted = biCamDSI.removeLTUntrusted(
        			  disparity_bimap,                           // double [][] disparity_bimap,
        			  clt_parameters.rig.lt_min_disparity,       // double      min_disparity, //  =          0.0;    // apply low texture to near objects
        			  clt_parameters.rig.lt_trusted_strength,    // double      trusted_strength, // =       0.2;    // strength sufficient without neighbors
        			  clt_parameters.rig.lt_need_friends,        // double      need_friends, // =           0.4;    // strength sufficient with neighbors support, fraction of lt_trusted_strength
        			  clt_parameters.rig.lt_friends_diff,        // double      friends_diff, // =           0.2;    // pix difference to neighbors to be considered a match (TODO: use tilted)
        			  clt_parameters.rig.lt_friends_rdiff,       // double      friends_rdiff, // =          0.04;   // additional relative pix per pixel of disparity
        			  clt_parameters.rig.lt_min_friends_any,     // int         min_friends_any, // =        2;      // minimal number of even weak friends
        			  clt_parameters.rig.lt_min_friends_trusted, // int         min_friends_trusted, // =    2;      // minimal number of trusted (strong or already confirmed)
        			  clt_parameters.rig.lt_friends_dist,        // int         friends_dist, // =           3;      // how far to look for friends
        			  debugLevel); // int         debugLevel
        	  if ((num_new_trusted - num_trusted) < min_added_tiles) {
        		  if (debugLevel > -2) {
        			  System.out.println("enhanceByRig(): pass="+num_fill+", number of added tiles = "+(num_new_trusted - num_trusted)+" < " +min_added_tiles+", done adding");
        			  break;
        		  }
        	  } else {
        		  if (debugLevel > -2) {
        			  System.out.println("enhanceByRig(): pass="+num_fill+", number of added tiles = "+(num_new_trusted - num_trusted));
        		  }
        	  }
        	  num_trusted = num_new_trusted;
        	  disparity_bimap = fillPoorTextureByInter(
        			  quadCLT_main,     // QuadCLT            quadCLT_main,  // tiles should be set
        			  quadCLT_aux,      // QuadCLT            quadCLT_aux,
        			  clt_parameters,   // EyesisCorrectionParameters.CLTParameters       clt_parameters,
        			  disparity_bimap,  //double [][]                                    disparity_bimap,
        			  biCamDSI,          // BiCamDSI                                       biCamDSI,
        			  threadsMax,       // final int                                      threadsMax,  // maximal number of threads to launch
        			  updateStatus,     // final boolean                                  updateStatus,
        			  debugLevel-2);      // final int                                      debugLevel)// throws Exception


        	  if (clt_parameters.show_map &&  (debugLevel > 0) && clt_parameters.rig.rig_mode_debug){ //OK
        		  (new showDoubleFloatArrays()).showArrays(
        				  disparity_bimap,
        				  tilesX,
        				  disparity_bimap[0].length/tilesX,
        				  true,
        				  quadCLT_main.image_name+"DSI_LT-N"+num_fill,
        				  ImageDtt.BIDISPARITY_TITLES);
        	  }
          }
    	  int num_new_trusted = biCamDSI.removeLTUntrusted(
    			  disparity_bimap,                           // double [][] disparity_bimap,
    			  clt_parameters.rig.lt_min_disparity,       // double      min_disparity, //  =          0.0;    // apply low texture to near objects
    			  clt_parameters.rig.lt_trusted_strength,    // double      trusted_strength, // =       0.2;    // strength sufficient without neighbors
    			  clt_parameters.rig.lt_need_friends,        // double      need_friends, // =           0.4;    // strength sufficient with neighbors support, fraction of lt_trusted_strength
    			  clt_parameters.rig.lt_friends_diff,        // double      friends_diff, // =           0.2;    // pix difference to neighbors to be considered a match (TODO: use tilted)
    			  clt_parameters.rig.lt_friends_rdiff,       // double      friends_rdiff, // =          0.04;   // additional relative pix per pixel of disparity
    			  clt_parameters.rig.lt_min_friends_any,     // int         min_friends_any, // =        2;      // minimal number of even weak friends
    			  clt_parameters.rig.lt_min_friends_trusted, // int         min_friends_trusted, // =    2;      // minimal number of trusted (strong or already confirmed)
    			  clt_parameters.rig.lt_friends_dist,        // int         friends_dist, // =           3;      // how far to look for friends
    			  debugLevel); // int         debugLevel

    	  if (clt_parameters.show_map &&  (debugLevel > -2) && clt_parameters.rig.rig_mode_debug){ //OK
    		  System.out.println("There are total "+num_new_trusted+" trusted tiles in ground truth data");
    		  (new showDoubleFloatArrays()).showArrays(
    				  disparity_bimap,
    				  tilesX,
    				  disparity_bimap[0].length/tilesX,
    				  true,
    				  quadCLT_main.image_name+"DSI_LT-FINAL",
    				  ImageDtt.BIDISPARITY_TITLES);
    	  }
    	  if (clt_parameters.rig.lt_remove_stray) {
    		  int num_stray_removed = biCamDSI.removeFalseMatches(
    				  disparity_bimap,                        // double [][] disparity_bimap,
    				  clt_parameters.rig.lt_trusted_strength, //double      trusted_strength, // =       0.2;    // strength sufficient without neighbors
    				  clt_parameters.rig.lt_stray_rstrength,  // double      stray_rstrength, // = 2.0;    // Relative to trusted strength - trust above that
    				  clt_parameters.rig.lt_strength_rfloor,  // double      strength_rfloor, // = 0.28;    // Relative to trusted strength
    				  clt_parameters.rig.lt_stray_dist ,      //	int         stray_dist,      // = 2;      // How far to look (and erase) around a potentially falsely matched tile
    				  clt_parameters.rig.lt_stray_over ,      //	double      stray_over,      // = 2.0;    // Stray tile should be this stronger than the strongest neighbor to be recognized
    				  debugLevel);                            // int         debugLevel
			  if (debugLevel > -2) {
				  System.out.println("removeFalseMatches() found ="+num_stray_removed+" false matches and erased those tiles and tiles around");
			  }
    		  if (num_stray_removed > 0) {
    			  // refine again:
    			  num_trusted = 0;
    			  for (int num_fill = 0; num_fill < clt_parameters.rig.lt_repeat; num_fill++) {
    				  num_new_trusted = biCamDSI.removeLTUntrusted(
    						  disparity_bimap,                           // double [][] disparity_bimap,
    						  clt_parameters.rig.lt_min_disparity,       // double      min_disparity, //  =          0.0;    // apply low texture to near objects
    						  clt_parameters.rig.lt_trusted_strength,    // double      trusted_strength, // =       0.2;    // strength sufficient without neighbors
    						  clt_parameters.rig.lt_need_friends,        // double      need_friends, // =           0.4;    // strength sufficient with neighbors support, fraction of lt_trusted_strength
    						  clt_parameters.rig.lt_friends_diff,        // double      friends_diff, // =           0.2;    // pix difference to neighbors to be considered a match (TODO: use tilted)
    						  clt_parameters.rig.lt_friends_rdiff,       // double      friends_rdiff, // =          0.04;   // additional relative pix per pixel of disparity
    						  clt_parameters.rig.lt_min_friends_any,     // int         min_friends_any, // =        2;      // minimal number of even weak friends
    						  clt_parameters.rig.lt_min_friends_trusted, // int         min_friends_trusted, // =    2;      // minimal number of trusted (strong or already confirmed)
    						  clt_parameters.rig.lt_friends_dist,        // int         friends_dist, // =           3;      // how far to look for friends
    						  debugLevel); // int         debugLevel
    				  if ((num_new_trusted - num_trusted) < min_added_tiles) {
    					  if (debugLevel > -2) {
    						  System.out.println("enhanceByRig(): pass="+num_fill+", number of added tiles = "+(num_new_trusted - num_trusted)+" < " +min_added_tiles+", done adding");
    						  break;
    					  }
    				  } else {
    					  if (debugLevel > -2) {
    						  System.out.println("enhanceByRig(): pass="+num_fill+", number of added tiles = "+(num_new_trusted - num_trusted));
    					  }
    				  }
    				  num_trusted = num_new_trusted;
    				  disparity_bimap = fillPoorTextureByInter(
    						  quadCLT_main,     // QuadCLT            quadCLT_main,  // tiles should be set
    						  quadCLT_aux,      // QuadCLT            quadCLT_aux,
    						  clt_parameters,   // EyesisCorrectionParameters.CLTParameters       clt_parameters,
    						  disparity_bimap,  //double [][]                                    disparity_bimap,
    						  biCamDSI,          // BiCamDSI                                       biCamDSI,
    						  threadsMax,       // final int                                      threadsMax,  // maximal number of threads to launch
    						  updateStatus,     // final boolean                                  updateStatus,
    						  debugLevel-2);      // final int                                      debugLevel)// throws Exception


    				  if (clt_parameters.show_map &&  (debugLevel > 0) && clt_parameters.rig.rig_mode_debug){ //OK
    					  (new showDoubleFloatArrays()).showArrays(
    							  disparity_bimap,
    							  tilesX,
    							  disparity_bimap[0].length/tilesX,
    							  true,
    							  quadCLT_main.image_name+"DSI_LT-N"+num_fill,
    							  ImageDtt.BIDISPARITY_TITLES);
    				  }
    			  }
    		  }
    	  }

    	  double [][] rig_disparity_strength = {disparity_bimap[ImageDtt.BI_TARGET_INDEX],disparity_bimap[ImageDtt.BI_STR_CROSS_INDEX]};

    	  return rig_disparity_strength;
	  }




	  public double [][] fillPoorTextureByInter(
			  QuadCLT            quadCLT_main,  // tiles should be set
			  QuadCLT            quadCLT_aux,
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
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
    			  threadsMax,           // final int        threadsMax,  // maximal number of threads to launch
    			  updateStatus,         // final boolean    updateStatus,
    			  debugLevel);          // final int        debugLevel);

          if (clt_parameters.show_map &&  (debugLevel > -2) && clt_parameters.rig.rig_mode_debug){
        	  (new showDoubleFloatArrays()).showArrays(
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
        			  0.0, // clt_parameters.rig.refine_min_strength , // double refine_min_strength, // do not refine weaker tiles
        			  clt_parameters.rig.refine_tolerance ,    // double refine_tolerance,    // do not refine if absolute disparity below
        			  trusted_lt,      // null, // trusted_lt,    // tile_list,       // ArrayList<Integer>             tile_list,       // or null
        			  num_new,         // int     []                                     num_new,
        			  clt_parameters,  // EyesisCorrectionParameters.CLTParameters clt_parameters,
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
        	  (new showDoubleFloatArrays()).showArrays(
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
        	  (new showDoubleFloatArrays()).showArrays(
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
			  double               disp_offset,
			  QuadCLT              quadCLT_main,
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
		  final int tilesX = quadCLT_main.tp.getTilesX();
		  final int tilesY = quadCLT_main.tp.getTilesY();
		  int ml_width = 2 * ml_hwidth + 1;
		  int width =  tilesX * ml_width;
		  int height = tilesY * ml_width;
		  String title = ml_title+ (use8bpp?"08":"32")+"B-"+(keep_aux?"A":"")+(keep_inter?"I":"")+(keep_hor_vert?"O":"")+(ml_keep_tbrl?"T":"")+
				  (keep_debug?"D":"")+"-FZ"+ml_fatzero+"-OFFS"+disp_offset;
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
			  byte [][] iml_data = new byte [ml_data.length][];
			  for (int nl = 0; nl < ml_data.length; nl++) if (!skip_layers[nl]) {
				  iml_data[nl] = new byte [ml_data[nl].length];
				  if (nl == ImageDtt.ML_OTHER_INDEX) {
					  // special treatment - make 2 bytes of one disparity value
					  for (int tileY = 0; tileY < tilesY; tileY++) {
						  for (int tileX = 0; tileX < tilesX; tileX++) {
							  int nTile = tileY * tilesX + tileX;
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
							  iml_data[nl][i] = (byte) iv; //  (iv - 127); // NaN will stay 0;
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

	      double disparityRadiusMain =  quadCLT_main.geometryCorrection.getDisparityRadius();
	      double disparityRadiusAux =   quadCLT_aux.geometryCorrection.getDisparityRadius();
	      double intercameraBaseline =  quadCLT_aux.geometryCorrection.getBaseline();

	      ImagePlus imp_ml = new ImagePlus(title, array_stack);
	      imp_ml.setProperty("VERSION",  "1.0");
	      imp_ml.setProperty("tileWidth",   ""+ml_width);
	      imp_ml.setProperty("dispOffset",  ""+disp_offset);
	      imp_ml.setProperty("ML_OTHER_TARGET",  ""+ImageDtt.ML_OTHER_TARGET);
	      imp_ml.setProperty("ML_OTHER_GTRUTH",  ""+ImageDtt.ML_OTHER_GTRUTH);
	      imp_ml.setProperty("ML_OTHER_GTRUTH_STRENGTH",  ""+ImageDtt.ML_OTHER_GTRUTH_STRENGTH);
	      imp_ml.setProperty("disparityRadiusMain",  ""+disparityRadiusMain);
	      imp_ml.setProperty("disparityRadiusAux",  ""+disparityRadiusAux);
	      imp_ml.setProperty("intercameraBaseline",  ""+intercameraBaseline);
	      imp_ml.setProperty("data_min",  ""+soft_mn);
	      imp_ml.setProperty("data_max",  ""+soft_mx);

	      imp_ml.setProperty("comment_tileWidth",   "Square tile size for each 2d correlation, always odd");
	      imp_ml.setProperty("comment_dispOffset",  "Tile target disparity minum ground truth disparity");
	      imp_ml.setProperty("comment_ML_OTHER_TARGET",  "Offset of the target disparity in the \"other\" layer tile");
	      imp_ml.setProperty("comment_ML_OTHER_GTRUTH",  "Offset of the ground truth disparity in the \"other\" layer tile");
	      imp_ml.setProperty("comment_ML_OTHER_GTRUTH_STRENGTH",  "Offset of the ground truth strength in the \"other\" layer tile");
	      imp_ml.setProperty("comment_data_min",  "Defined only for 8bpp mode - value, corresponding to -127 (-128 is NaN)");
	      imp_ml.setProperty("comment_data_max",  "Defined only for 8bpp mode - value, corresponding to +127 (-128 is NaN)");

	      imp_ml.setProperty("comment_disparityRadiusMain",  "Side of the square where 4 main camera subcameras are located (mm)");
	      imp_ml.setProperty("comment_disparityRadiusAux",  "Side of the square where 4 main camera subcameras are located (mm). Disparity is specified for the main camera");
	      imp_ml.setProperty("comment_intercameraBaseline",  "Horizontal distance between the main and the auxiliary camera centers (mm). Disparity is specified for the main camera");
	      (new JP46_Reader_camera(false)).encodeProperiesToInfo(imp_ml);
	      imp_ml.getProcessor().resetMinAndMax();
	      if (show ) {
			   imp_ml.show();
 	      }
	      String path = ml_directory+=Prefs.getFileSeparator()+imp_ml.getTitle();
	      FileSaver fs=new FileSaver(imp_ml);
	      fs.saveAsTiff(path+".tiff");
	      if (debugLevel > -2) {
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
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
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
				  clt_parameters.fat_zero, // double                                         fatzero,
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
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
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
				  clt_parameters.fat_zero, // double                                         fatzero,
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
		  for (int i = 0; i < trusted.length; i++) {
			  trusted[i] = (Math.abs(bimap[ImageDtt.BI_DISP_CROSS_INDEX][i]) <= trusted_inter) &&
					  (Math.abs(bimap[ImageDtt.BI_DISP_FULL_INDEX][i])  <= trusted_main) &&
					  (Math.abs(bimap[ImageDtt.BI_ADISP_FULL_INDEX][i]) <= trusted_aux) &&
					  (bimap[ImageDtt.BI_STR_ALL_INDEX][i]              >= min_combo_strength) &&
					  ((was_trusted == null) || was_trusted[i]);
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
			  QuadCLT                                        quadCLT_main,  // tiles should be set
			  QuadCLT                                        quadCLT_aux,
			  double [][]                                    src_bimap, // current state of measurements
			  double [][]                                    prev_bimap, // previous state of measurements or null
			  double []                                      scale_bad,
			  int                                            refine_mode, // 0 - by main, 1 - by aux, 2 - by inter
			  boolean                                        keep_inf,    // keep expected disparity 0.0 if it was so
			  double                                         refine_min_strength, // do not refine weaker tiles
			  double                                         refine_tolerance,    // do not refine if absolute disparity below
			  ArrayList<Integer>                             tile_list, // or null
			  int     []                                     num_new,
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel){
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
				  refine_min_strength, // do not refine weaker tiles
				  refine_tolerance,    // do not refine if absolute disparity below
				  selection,
				  num_new,
				  clt_parameters,
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
			  double                                         refine_min_strength, // do not refine weaker tiles
			  double                                         refine_tolerance,    // do not refine if absolute disparity below
			  boolean []                                     selection,
			  int     []                                     num_new,
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
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
							  clt_parameters, // EyesisCorrectionParameters.CLTParameters       clt_parameters,
							  tile_op_all,    // int                                            tile_op_all,
							  src_bimap, // double [][]                                     src_bimap, // current state of measurements
							  prev_bimap, // double [][]                                    prev_bimap, // previous state of measurements or null
							  scale_bad,  // double []                                      scale_bad,
							  tile_op, // int [][]                                          tile_op, // common for both amin and aux
							  disparity_array, // double [][]                                    disparity_array,
							  refine_mode, // int                                            refine_mode, // 0 - by main, 1 - by aux, 2 - by inter
							  keep_inf,    // boolean                                        keep_inf,    // keep expected disparity 0.0 if it was so
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
				  clt_parameters.fat_zero, // double                                         fatzero,
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

	  /**
	   * Add measurements with new specified disparity of the main camera
	   * @param quadCLT_main main camera QuadCLT instance (should have tp initialized)
	   * @param quadCLT_aux auxiliary camera QuadCLT instance (should have tp initialized)
	   * @param src_bimap results of the older measurements (now includes expected disparity) or null (no old results available)
	   * @param disparity new expected disparity value to try
	   * @param tile_list list of selected tiles or null. If not null, will not re-measure listed tiles
	   * @param saturation_main saturated pixels bitmaps for the main camera
	   * @param saturation_aux saturated pixels bitmaps for the auxiliary camera
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
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
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
				  clt_parameters.fat_zero, // double                                         fatzero,
				  threadsMax,          //final int        threadsMax,  // maximal number of threads to launch
				  updateStatus,        // final boolean    updateStatus,
				  debugLevel);          // final int        debugLevel)

		  // combine with old results (if available) for the tiles that were not re-measured
		  if (src_bimap != null) {
			  for (int tileY = 0; tileY<tilesY;tileY++) {
				  for (int tileX = 0; tileX<tilesX;tileX++) {
					  int nTile = tileY * tilesX + tileX;
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
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
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
				  clt_parameters.fat_zero, // double                                         fatzero,
				  threadsMax,          //final int        threadsMax,  // maximal number of threads to launch
				  updateStatus,        // final boolean    updateStatus,
				  debugLevel);          // final int        debugLevel)
		  return disparity_bimap;
	  }




	  private boolean prepRefineTile(
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
			  int                                            tile_op_all,
			  double [][]                                    src_bimap, // current state of measurements
			  double [][]                                    prev_bimap, // previous state of measurements or null
			  double []                                      scale_bad,
			  int [][]                                       tile_op, // common for both amin and aux
			  double [][]                                    disparity_array,
			  int                                            refine_mode, // 0 - by main, 1 - by aux, 2 - by inter, 3 - inter-dx
			  boolean                                        keep_inf,    // keep expected disparity 0.0 if it was so
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
		  if (keep_inf && (src_bimap[ImageDtt.BI_TARGET_INDEX][nTile] == 0.0)) {
			  if ((refine_min_strength == 0.0) || (refine_tolerance == 0.0)) {
				  tile_op[tileY][tileX] = tile_op_all;
				  disparity_array[tileY][tileX] = 0.0;
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


	  private double [][] measureRig(
			  QuadCLT                                        quadCLT_main,  // tiles should be set
			  QuadCLT                                        quadCLT_aux,
			  int [][]                                       tile_op, // common for both amin and aux
			  double [][]                                    disparity_array,
			  double [][]                                    ml_data,         // data for ML - 10 layers - 4 center areas (3x3, 5x5,..) per camera-per direction, 1 - composite, and 1 with just 1 data (target disparity)
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
			  double                                         fatzero,
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel){
		  ImageDtt image_dtt = new ImageDtt();

		  double [][] disparity_bimap  = new double [ImageDtt.BIDISPARITY_TITLES.length][]; //[0] -residual disparity, [1] - orthogonal (just for debugging) last 4 - max pixel differences

		  image_dtt.clt_bi_quad (
				  clt_parameters,                       // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
				  fatzero,                              // final double              fatzero,         // May use correlation fat zero from 2 different parameters - fat_zero and rig.ml_fatzero
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
				  quadCLT_main.tp.getTilesX()*clt_parameters.transform_size, // final int                 width,

				  quadCLT_main.getGeometryCorrection(), // final GeometryCorrection  geometryCorrection_main,
				  quadCLT_aux.getGeometryCorrection(),  // final GeometryCorrection  geometryCorrection_aux,
				  quadCLT_main.getCLTKernels(),         // final double [][][][][][] clt_kernels_main, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				  quadCLT_aux.getCLTKernels(),          // final double [][][][][][] clt_kernels_aux,  // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				  clt_parameters.corr_magic_scale,      // final double              corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
				  false, // true,                                 // 	final boolean             keep_clt_data,
				  threadsMax,                           // final int                 threadsMax,  // maximal number of threads to launch
				  debugLevel-2);                        // final int                 globalDebugLevel);
		  return disparity_bimap;
	  }

	  public double [][] remeasureRigML(
			  double                                         disparity_offset,
			  QuadCLT                                        quadCLT_main,  // tiles should be set
			  QuadCLT                                        quadCLT_aux,
//			  double [][]                                    src_bimap,
			  double []                                      disparity,
			  double []                                      strength,
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
			  int                                            ml_hwidth,
			  double                                         fatzero,
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

//		  double [] disparity =    src_bimap[ImageDtt.BI_TARGET_INDEX];
//		  double [] strength =     src_bimap[ImageDtt.BI_STR_CROSS_INDEX];
		  boolean [] selection =   new boolean [strength.length];
		  for (int nTile = 0; nTile < selection.length; nTile++) {
			  selection[nTile] = strength[nTile] > 0.0;
		  }

			Correlation2d corr2d = new Correlation2d(
					clt_parameters.img_dtt,              // ImageDttParameters  imgdtt_params,
					clt_parameters.transform_size,             // int transform_size,
					2.0,                        //  double wndx_scale, // (wndy scale is always 1.0)
					(debugLevel > -1));   //   boolean debug)

		  for (int nTile = 0; nTile < disparity.length; nTile++) {
			  if (((selection == null) || (selection[nTile]) && !Double.isNaN(disparity[nTile]))) {
				  int tileY = nTile / tilesX;
				  int tileX = nTile % tilesX;
				  tile_op[tileY][tileX] = tile_op_all;
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
				  threadsMax,       // maximal number of threads to launch // final int        threadsMax,  // maximal number of threads to launch
				  updateStatus,     // final boolean    updateStatus,
				  debugLevel);      // final int        debugLevel)



		  return ml_data;
	  }


	  ArrayList<Integer> selectRigTiles(
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
			  boolean select_infinity,
			  boolean select_noninfinity,
			  double[][] disparity_bimap,
			  int tilesX)
	  {
		  ArrayList<Integer> tilesList = new ArrayList<Integer>();
		  int numTiles = disparity_bimap[ImageDtt.BI_STR_FULL_INDEX].length;
		  if (select_infinity) {
			  for (int nTile = 0; nTile < numTiles; nTile++) {
				  if (    (disparity_bimap[ImageDtt.BI_TARGET_INDEX][nTile] == 0.0 ) && // expected disparity was 0.0 (infinity)
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
				  if (    (disparity_bimap[ImageDtt.BI_TARGET_INDEX][nTile] > 0.0 ) && // expected disparity was > 0.0 (not infinity)
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
			  (new showDoubleFloatArrays()).showArrays(
					  dbg_inf,
					  tilesX,
					  numTiles/tilesX,
					  true,
					  title,
					  titles );
	  }



}
