import java.util.ArrayList;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;

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
				  throw new Exception ("Set naims for cameras do not match: main camera: '"+set_channels_main[nSet].name()+"', aux. camera: nothing");
			  }
			  if (!set_channels_main[nSet].name().equals(set_channels_aux[nSet].name())) {
				  throw new Exception ("Set naims for cameras do not match: main camera: '"+set_channels_main[nSet].name()+"', aux. camera: '"+set_channels_main[nSet].name()+"'");
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
				  throw new Exception ("Set naims for cameras do not match: main camera: '"+set_channels_main[nSet].name()+"', aux. camera: nothing");
			  }
			  if (!set_channels_main[nSet].name().equals(set_channels_aux[nSet].name())) {
				  throw new Exception ("Set naims for cameras do not match: main camera: '"+set_channels_main[nSet].name()+"', aux. camera: '"+set_channels_main[nSet].name()+"'");
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
		  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging? - TODO - move whete it belongs
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

		  double [][][][] double_stacks = 	  getRigImageStacks(
				  clt_parameters,  // EyesisCorrectionParameters.CLTParameters       clt_parameters,
				  quadCLT_main,    // QuadCLT                                         quadCLT_main,
				  quadCLT_aux,     // QuadCLT                                          quadCLT_aux,
				  imp_quad_main,   // ImagePlus []                                   imp_quad_main,
				  imp_quad_aux,    // ImagePlus []                                    imp_quad_aux,
				  threadsMax,      // maximal number of threads to launch
				  debugLevel);     // final int        debugLevel);

		  // temporary setting up tile task file (one integer per tile, bitmask
		  // for testing defined for a window, later the tiles to process will be calculated based on previous passes results

		  int [][]    tile_op_main = quadCLT_main.tp.setSameTileOp(clt_parameters,  clt_parameters.tile_task_op, debugLevel);
		  int [][]    tile_op_aux =  quadCLT_aux.tp.setSameTileOp (clt_parameters,  clt_parameters.tile_task_op, debugLevel);

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

		  double min_corr_selected = clt_parameters.min_corr;

		  ImageDtt image_dtt = new ImageDtt();

		  final double [][][][][][][] clt_bidata = // new double[2][quad][nChn][tilesY][tilesX][][]; // first index - main/aux
				  image_dtt.clt_bi_quad (
						  clt_parameters,                       // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
						  tile_op_main,                         // final int [][]            tile_op_main,    // [tilesY][tilesX] - what to do - 0 - nothing for this tile
						  tile_op_aux,                          // final int [][]            tile_op_aux,     // [tilesY][tilesX] - what to do - 0 - nothing for this tile
						  disparity_array_main,                 // final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
						  double_stacks[0],                     // final double [][][]       image_data_main, // first index - number of image in a quad
						  double_stacks[1],                     // final double [][][]       image_data_aux,  // first index - number of image in a quad
						  saturation_main,                      // final boolean [][]        saturation_main, // (near) saturated pixels or null
						  saturation_aux,                       // final boolean [][]        saturation_aux,  // (near) saturated pixels or null
						  // correlation results - combo will be for the correation between two quad cameras
						  clt_corr_combo,                       // final double [][][][]     clt_corr_combo,  // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
						  // [type][tilesY][tilesX] should be set by caller
						  // types: 0 - selected correlation (product+offset), 1 - sum
						  disparity_bimap,                      // final double [][]    disparity_bimap, // [23][tilesY][tilesX]
						  texture_tiles_main,                   // final double [][][][]     texture_tiles_main, // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
						  texture_tiles_aux,                    // final double [][][][]     texture_tiles_aux,  // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
						  imp_quad_main[0].getWidth(),          // final int                 width,

						  min_corr_selected,                    // final double              min_corr,        // 0.02; // minimal correlation value to consider valid
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
					  quadCLT_main.linearStackToColor(
							  clt_parameters,
							  colorProcParameters,
							  rgbParameters,
							  name+"-texture", // String name,
							  "-D"+clt_parameters.disparity+"-MAIN", //String suffix, // such as disparity=...
							  toRGB,
							  !quadCLT_main.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
							  true, // boolean saveShowIntermediate, // save/show if set globally
							  true, // boolean saveShowFinal,        // save/show result (color image?)
							  ((clt_parameters.alpha1 > 0)? texture_rgba_main: texture_rgb_main),
							  tilesX *  clt_parameters.transform_size,
							  tilesY *  clt_parameters.transform_size,
							  1.0,         // double scaleExposure, // is it needed?
							  debugLevel );
					  quadCLT_aux.linearStackToColor(
							  clt_parameters,
							  colorProcParameters,
							  rgbParameters,
							  name+"-texture", // String name,
							  "-D"+clt_parameters.disparity+"-AUX", //String suffix, // such as disparity=...
							  toRGB,
							  !quadCLT_aux.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
							  true, // boolean saveShowIntermediate, // save/show if set globally
							  true, // boolean saveShowFinal,        // save/show result (color image?)
							  ((clt_parameters.alpha1 > 0)? texture_rgba_aux: texture_rgb_aux),
							  tilesX *  clt_parameters.transform_size,
							  tilesY *  clt_parameters.transform_size,
							  1.0,         // double scaleExposure, // is it needed?
							  debugLevel );
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





	  public double [][][][] getRigImageStacks(
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
			  QuadCLT                                        quadCLT_main,
			  QuadCLT                                        quadCLT_aux,
			  ImagePlus []                                   imp_quad_main,
			  ImagePlus []                                   imp_quad_aux,
			  int                                            threadsMax,  // maximal number of threads to launch
			  int                                            debugLevel){
		  double [][][] double_stacks_main = new double [imp_quad_main.length][][];
		  for (int i = 0; i < double_stacks_main.length; i++){
			  double_stacks_main[i] = quadCLT_main.eyesisCorrections.bayerToDoubleStack(
					  imp_quad_main[i], // source Bayer image, linearized, 32-bit (float))
					  null); // no margins, no oversample
		  }

		  double [][][] double_stacks_aux = new double [imp_quad_main.length][][];
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
		  double [][][][] rigStacks = {double_stacks_main,double_stacks_aux};
		  quadCLT_main.setTiles (imp_quad_main[0], // set global tp.tilesX, tp.tilesY
				  clt_parameters,
				  threadsMax);
		  quadCLT_aux.setTiles (imp_quad_aux[0], // set global tp.tilesX, tp.tilesY
				  clt_parameters,
				  threadsMax);

		  return rigStacks;
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
				  throw new Exception ("Set naims for cameras do not match: main camera: '"+set_channels_main[nSet].name()+"', aux. camera: nothing");
			  }
			  if (!set_channels_main[nSet].name().equals(set_channels_aux[nSet].name())) {
				  throw new Exception ("Set naims for cameras do not match: main camera: '"+set_channels_main[nSet].name()+"', aux. camera: '"+set_channels_main[nSet].name()+"'");
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
			  final int        debugLevel){

		  double [][][][] double_stacks = 	  getRigImageStacks(
				  clt_parameters,  // EyesisCorrectionParameters.CLTParameters       clt_parameters,
				  quadCLT_main,    // QuadCLT                                         quadCLT_main,
				  quadCLT_aux,     // QuadCLT                                          quadCLT_aux,
				  imp_quad_main,   // ImagePlus []                                   imp_quad_main,
				  imp_quad_aux,    // ImagePlus []                                    imp_quad_aux,
				  threadsMax,      // maximal number of threads to launch
				  debugLevel);     // final int        debugLevel);

		  final int tilesX = quadCLT_main.tp.getTilesX();

		  double [][] disparity_bimap = measureInfinityRig(
				  quadCLT_main,  //QuadCLT                                        quadCLT_main,  // tiles should be set
				  quadCLT_aux,  //QuadCLT                                        quadCLT_aux,
				  double_stacks, // double [][][][]                                double_stacks,
				  null, // ArrayList<Integer>                             tileList,      // or null
				  saturation_main, // boolean [][]                                   saturation_main, // (near) saturated pixels or null
				  saturation_aux, //boolean [][]                                   saturation_aux, // (near) saturated pixels or null
			       clt_parameters, // EyesisCorrectionParameters.CLTParameters       clt_parameters,
				  threadsMax,  // final int        threadsMax,  // maximal number of threads to launch
				  updateStatus, // final boolean    updateStatus,
				  debugLevel); // final int        debugLevel);

		  if (disparity_bimap != null) {
			  ArrayList<Integer> tileList =  selectInfinityTiles(
					  clt_parameters,  // EyesisCorrectionParameters.CLTParameters       clt_parameters,
					  disparity_bimap, // double[][] disparity_bimap,
					  tilesX,          // int tilesX,
					  clt_parameters.rig.rig_mode_debug); //  boolean debug
			  System.out.println("Selected "+tileList.size()+" tiles for infinity correction");
		  }
// loop here with the same tileList

		  return true;
	  }



	  public double [][] measureInfinityRig(
			  QuadCLT                                        quadCLT_main,  // tiles should be set
			  QuadCLT                                        quadCLT_aux,
			  double [][][][]                                double_stacks,
			  ArrayList<Integer>                             tileList,      // or null
			  boolean [][]                                   saturation_main, // (near) saturated pixels or null
			  boolean [][]                                   saturation_aux, // (near) saturated pixels or null
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel){

		  int tilesX = quadCLT_main.tp.getTilesX();
		  int tilesY = quadCLT_main.tp.getTilesX();
		  int [][]                                       tile_op;
		  if (tileList == null ) {
			  tile_op = quadCLT_main.tp.setSameTileOp(clt_parameters,  clt_parameters.tile_task_op, debugLevel);
		  } else {
			  tile_op = new int [tilesY][tilesX];
			  for (int nTile:tileList) {
				  tile_op[nTile/tilesX][nTile%tilesX] = clt_parameters.tile_task_op;
			  }
		  }
		  // for infinity use disparity of all 0-s;
		  double [][]                                    disparity_array = new double [quadCLT_main.tp.getTilesY()][quadCLT_main.tp.getTilesX()];
		  return measureRig(
				  quadCLT_main, // QuadCLT                                        quadCLT_main,  // tiles should be set
				  quadCLT_aux,  // QuadCLT                                        quadCLT_aux,
				  double_stacks, // double [][][][]                                double_stacks,
				  tile_op,  // int [][]                                       tile_op, // common for both amin and aux
				  disparity_array, // double [][]                                    disparity_array,
				  saturation_main, // boolean [][]                                   saturation_main, // (near) saturated pixels or null
				  saturation_aux, // boolean [][]                                   saturation_aux, // (near) saturated pixels or null
				  clt_parameters, // EyesisCorrectionParameters.CLTParameters       clt_parameters,
				  threadsMax,  //final int        threadsMax,  // maximal number of threads to launch
				  updateStatus, // final boolean    updateStatus,
				  debugLevel); // final int        debugLevel)

	  }




	  public double [][] measureRig(
			  QuadCLT                                        quadCLT_main,  // tiles should be set
			  QuadCLT                                        quadCLT_aux,
			  double [][][][]                                double_stacks,
			  int [][]                                       tile_op, // common for both amin and aux
			  double [][]                                    disparity_array,
			  boolean [][]                                   saturation_main, // (near) saturated pixels or null
			  boolean [][]                                   saturation_aux, // (near) saturated pixels or null
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel){
		  ImageDtt image_dtt = new ImageDtt();

		  double [][] disparity_bimap  = new double [ImageDtt.BIDISPARITY_TITLES.length][]; //[0] -residual disparity, [1] - orthogonal (just for debugging) last 4 - max pixel differences

		  double min_corr_selected = clt_parameters.min_corr;
		  image_dtt.clt_bi_quad (
				  clt_parameters,                       // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
				  tile_op,                              // final int [][]            tile_op_main,    // [tilesY][tilesX] - what to do - 0 - nothing for this tile
				  tile_op,                              // final int [][]            tile_op_aux,     // [tilesY][tilesX] - what to do - 0 - nothing for this tile
				  disparity_array,                      // final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
				  double_stacks[0],                     // final double [][][]       image_data_main, // first index - number of image in a quad
				  double_stacks[1],                     // final double [][][]       image_data_aux,  // first index - number of image in a quad
				  saturation_main,                      // final boolean [][]        saturation_main, // (near) saturated pixels or null
				  saturation_aux,                       // final boolean [][]        saturation_aux,  // (near) saturated pixels or null
				  // correlation results - combo will be for the correation between two quad cameras
				  null,                                 // final double [][][][]     clt_corr_combo,  // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				  disparity_bimap,                      // final double [][]    disparity_bimap, // [23][tilesY][tilesX]
				  null,                                 // final double [][][][]     texture_tiles_main, // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
				  null,                                 // final double [][][][]     texture_tiles_aux,  // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
				  quadCLT_main.tp.getTilesX()*clt_parameters.transform_size, // final int                 width,

				  min_corr_selected,                    // final double              min_corr,        // 0.02; // minimal correlation value to consider valid
				  quadCLT_main.getGeometryCorrection(), // final GeometryCorrection  geometryCorrection_main,
				  quadCLT_aux.getGeometryCorrection(),  // final GeometryCorrection  geometryCorrection_aux,
				  quadCLT_main.getCLTKernels(),         // final double [][][][][][] clt_kernels_main, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				  quadCLT_aux.getCLTKernels(),          // final double [][][][][][] clt_kernels_aux,  // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				  clt_parameters.corr_magic_scale,      // final double              corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
				  false, // true,                                 // 	final boolean             keep_clt_data,
				  threadsMax,                           // final int                 threadsMax,  // maximal number of threads to launch
				  debugLevel);                          // final int                 globalDebugLevel);
		  return disparity_bimap;
	  }





	  ArrayList<Integer> selectInfinityTiles(
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
			  double[][] disparity_bimap,
			  int tilesX,
			  boolean debug
			  ){
		  ArrayList<Integer> tilesList = new ArrayList<Integer>();
		  int numTiles = disparity_bimap[ImageDtt.BI_STR_FULL_INDEX].length;
		  for (int nTile = 0; nTile < numTiles; nTile++) {
			  if ((disparity_bimap[ImageDtt.BI_STR_FULL_INDEX][nTile] >= clt_parameters.rig.inf_min_strength_main)&&
					  (disparity_bimap[ImageDtt.BI_ASTR_FULL_INDEX][nTile] >= clt_parameters.rig.inf_min_strength_aux)&&
					  (disparity_bimap[ImageDtt.BI_STR_CROSS_INDEX][nTile] >= clt_parameters.rig.inf_min_strength_rig)&&
					  (Math.abs(disparity_bimap[ImageDtt.BI_DISP_FULL_INDEX][nTile]) <= clt_parameters.rig.inf_max_disp_main)&&
					  (Math.abs(disparity_bimap[ImageDtt.BI_ADISP_FULL_INDEX][nTile]) <= clt_parameters.rig.inf_max_disp_aux)&&
					  (Math.abs(disparity_bimap[ImageDtt.BI_DISP_CROSS_INDEX][nTile]) <= clt_parameters.rig.inf_max_disp_rig)) {
				  tilesList.add(nTile);
			  }
		  }
		  if (debug) {
			  String [] titles = {"disparity", "dx","dy","strength"};
			  double [][] dbg_inf = new double [4][numTiles];
			  for (int nTile = 0; nTile < numTiles; nTile++){
				  dbg_inf[0][nTile] = Double.NaN;
				  dbg_inf[1][nTile] = Double.NaN;
				  dbg_inf[2][nTile] = Double.NaN;
				  dbg_inf[3][nTile] = Double.NaN; // 0.0;
			  }
			  for (int nTile:tilesList) {
				  dbg_inf[0][nTile] = disparity_bimap[ImageDtt.BI_DISP_CROSS_INDEX][nTile];
				  dbg_inf[1][nTile] = disparity_bimap[ImageDtt.BI_DISP_CROSS_DX_INDEX][nTile];
				  dbg_inf[2][nTile] = disparity_bimap[ImageDtt.BI_DISP_CROSS_DY_INDEX][nTile];
				  dbg_inf[3][nTile] = disparity_bimap[ImageDtt.BI_STR_CROSS_INDEX][nTile];
			  }
			  (new showDoubleFloatArrays()).showArrays(
					  dbg_inf,
					  tilesX,
					  numTiles/tilesX,
					  true,
					  "Infinity_selected_data",
					  titles );
		  }
//new showDoubleFloatArrays()
		  return tilesList;
	  }





}
