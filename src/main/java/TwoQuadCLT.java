import java.util.ArrayList;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.io.FileSaver;

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

		  double min_corr_selected = clt_parameters.min_corr;

		  ImageDtt image_dtt = new ImageDtt();
		  double [][] ml_data = null;

		  final double [][][][][][][] clt_bidata = // new double[2][quad][nChn][tilesY][tilesX][][]; // first index - main/aux
				  image_dtt.clt_bi_quad (
						  clt_parameters,                       // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
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



		  final int tilesX = quadCLT_main.tp.getTilesX();

// perform full re-measure cycles
		  double [][] disparity_bimap = null;
          int [] num_new = new int[1];
		  for (int num_full_cycle = 0; num_full_cycle < clt_parameters.rig.rig_adjust_full_cycles;num_full_cycle++) {
			  disparity_bimap = null;
			  ArrayList<Integer> tile_list = new ArrayList<Integer>();
			  // measuere and refine
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
						  tile_list,         // ArrayList<Integer>  tile_list,       // or null. If non-null - do not remeasure members of teh list
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

			  // short cycle (remeasure only for teh list). May also break from it if RMS is not improving
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
	  public double [][] enhanceByRig(
			  QuadCLT            quadCLT_main,  // tiles should be set
			  QuadCLT            quadCLT_aux,
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
			  final int                                      threadsMax,  // maximal number of threads to launch
			  final boolean                                  updateStatus,
			  final int                                      debugLevel) throws Exception
	  {
//		  int min_new = 100; // make a parameter
//		  int num_inf_refine =  20;
//		  int num_near_refine = 20;
//		  double min_trusted_strength = 0.1; // 14;
//		  double trusted_tolerance  =   1.0;
		  System.out.println("enhanceByRig()");
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
		  double [][] disparity_bimap_infinity = measureNewRigDisparity(
				  quadCLT_main,      // QuadCLT             quadCLT_main,    // tiles should be set
				  quadCLT_aux,       // QuadCLT             quadCLT_aux,
				  null,              // double [][]         src_bimap,       // current state of measurements (or null for new measurement)
				  0,                 // double              disparity,
				  null,              // ArrayList<Integer>  tile_list,       // or null. If non-null - do not remeasure members of teh list
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

          if (clt_parameters.show_map &&  (debugLevel > -2) && clt_parameters.rig.rig_mode_debug){
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
        			  2,               // int                            refine_mode,     // 0 - by main, 1 - by aux, 2 - by inter
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

          if (clt_parameters.show_map &&  (debugLevel > -2) && clt_parameters.rig.rig_mode_debug){

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
          CLTPass3d scan_last = quadCLT_main.tp.clt_3d_passes.get( quadCLT_main.tp.clt_3d_passes.size() -1); // get last one
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
          if (clt_parameters.show_map &&  (debugLevel > -2) && clt_parameters.rig.rig_mode_debug){
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
        			  2,               // int                            refine_mode,     // 0 - by main, 1 - by aux, 2 - by inter
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

          if (clt_parameters.show_map &&  (debugLevel > -2) && clt_parameters.rig.rig_mode_debug){

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

          if (clt_parameters.show_map &&  (debugLevel > -2) && clt_parameters.rig.rig_mode_debug){
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

          // grow around using all camera and inter-camera correlations (try to get low-textured,like our street)

          // re-measure with ML data output
//		  int ml_hwidth =      2; // move to clt_parameters
//		  int ml_width = 2 * clt_parameters.rig.ml_hwidth + 1;
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
	    			  disparity_bimap,                          // double [][]                                    src_bimap,
	    			  clt_parameters,                           // EyesisCorrectionParameters.CLTParameters clt_parameters,
	    			  clt_parameters.rig.ml_hwidth,             // int ml_hwidth
	    			  threadsMax,                               // final int            threadsMax,  // maximal number of threads to launch
	    			  updateStatus,                             // final boolean        updateStatus,
	    			  debugLevel);                              // final int            debugLevel);
	    	  saveMlFile(
	    			  quadCLT_main.image_name+"-ML_DATA-OFFS_", // String               ml_title,
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
	    			  clt_parameters.rig.ml_keep_debug,         // boolean              keep_debug,
	    			  clt_parameters.rig.ml_hwidth,             // int                  ml_hwidth,
	    			  ml_data,                                  // double [][]          ml_data,
	    			  clt_parameters.rig.ml_show_ml,            // boolean              show,
	    			  debugLevel);                              // int                  debugLevel
		  }
          return null;
          //clt_3d_passes
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
			  boolean              keep_debug,
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
		  String title = ml_title+(keep_aux?"A":"")+(keep_inter?"I":"")+(keep_hor_vert?"O":"")+(keep_debug?"D":"")+disp_offset;
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
    			  threadsMax,      // final int        threadsMax,  // maximal number of threads to launch
    			  updateStatus,    // final boolean    updateStatus,
    			  debugLevel);     // final int        debugLevel);

		  return disparity_bimap;
	  }
	  /**
	   * Select tiles that have small enough residual disparity to trust on each of the two cameras and the whole rig.
	   * @param quadCLT_main main camera QuadCLT instance (should have tp initialized)
	   * @param quadCLT_aux auxiliary camera QuadCLT instance (should have tp initialized)
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
			  ArrayList<Integer>                             tile_list, // or null. If non-null - do not remeasure members of teh list
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel){
		  if (debugLevel > 0) {
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

	  private boolean prepRefineTile(
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
			  int                                            tile_op_all,
			  double [][]                                    src_bimap, // current state of measurements
			  double [][]                                    prev_bimap, // previous state of measurements or null
			  double []                                      scale_bad,
			  int [][]                                       tile_op, // common for both amin and aux
			  double [][]                                    disparity_array,
			  int                                            refine_mode, // 0 - by main, 1 - by aux, 2 - by inter
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

		  boolean debug_this = nTile==40661; // 61924;
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
		  default:
			  diff_disp = src_bimap[ImageDtt.BI_DISP_CROSS_INDEX][nTile];
			  diff_prev= (prev_bimap == null)? Double.NaN:prev_bimap[ImageDtt.BI_DISP_CROSS_INDEX][nTile];
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

		  if (Math.abs((new_disp - ref_target)/new_disp) < refine_tolerance) return false;

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
			  final int        threadsMax,  // maximal number of threads to launch
			  final boolean    updateStatus,
			  final int        debugLevel){
		  ImageDtt image_dtt = new ImageDtt();

		  double [][] disparity_bimap  = new double [ImageDtt.BIDISPARITY_TITLES.length][]; //[0] -residual disparity, [1] - orthogonal (just for debugging) last 4 - max pixel differences

		  double min_corr_selected = clt_parameters.min_corr;
		  image_dtt.clt_bi_quad (
				  clt_parameters,                       // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
				  tile_op,                              // final int [][]            tile_op_main,    // [tilesY][tilesX] - what to do - 0 - nothing for this tile
//				  tile_op,                              // final int [][]            tile_op_aux,     // [tilesY][tilesX] - what to do - 0 - nothing for this tile
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

				  min_corr_selected,                    // final double              min_corr,        // 0.02; // minimal correlation value to consider valid
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
			  double [][]                                    src_bimap,
			  EyesisCorrectionParameters.CLTParameters       clt_parameters,
			  int                                            ml_hwidth,
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

		  double [] disparity =    src_bimap[ImageDtt.BI_TARGET_INDEX];
		  double [] strength =     src_bimap[ImageDtt.BI_STR_ALL_INDEX];
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
