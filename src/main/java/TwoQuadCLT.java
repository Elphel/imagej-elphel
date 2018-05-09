import ij.IJ;
import ij.ImagePlus;

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

/*
	public QuadCLT quadCLT_main;
	public QuadCLT quadCLT_aux;
	public TwoQuadCLT (
			QuadCLT quadCLT_main,
			QuadCLT quadCLT_aux
			) {
		this.quadCLT_main = quadCLT_main;
		this.quadCLT_aux =  quadCLT_aux;
	}
*/

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

}
