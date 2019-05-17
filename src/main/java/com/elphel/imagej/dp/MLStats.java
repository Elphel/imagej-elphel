package com.elphel.imagej.dp;
/**
 ** MLStats - Generate reports over multiple scenes
 **
 ** Copyright (C) 2018 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  MLStats.java is free software: you can redistribute it and/or modify
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

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.FileVisitOption;
import java.nio.file.FileVisitResult;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.SimpleFileVisitor;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.EnumSet;
import java.util.List;
import java.util.Properties;

import com.elphel.imagej.cameras.EyesisCorrectionParameters;
import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.jp4.JP46_Reader_camera;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.text.TextWindow;

public class MLStats {
//	TwoQuadCLT twoQuadCLT;
//	public MLStats(TwoQuadCLT twoQuadCLT) {
//		this.twoQuadCLT = twoQuadCLT;
//	}
	public static boolean dsiHistogram(String dir) {
		Path path= Paths.get(dir);
		int disparity_bins =        1000;
		int strength_bins =          100;

		double disparity_min_drop =  -0.1;
		double disparity_min_clip =  -0.1;
		double disparity_max_drop = 100.0; //
		double disparity_max_clip = 100.0; //
		double strength_min_drop =    0.1;
		double strength_min_clip =    0.1;
		double strength_max_drop =    1.0; //
		double strength_max_clip =    0.9; //
		boolean normalize =           true;
		double master_weight_power =  1.0;
		double master_weight_floor =  0.08;
		double disparity_outlier =    1.0;
		double pre_log_offs =         0.01; // add before log to avoid -infinity
		double log_sigma =            2.00; // blur logarithm of the histogram (in bins)
		double mask_threshold =       0.25; // relative tile population
		double log_sigma_d =          0.5;  // blur logarithm of the histogram in disparity pixels
		double log_sigma_s =          0.01; // blur logarithm of the histogram in strength units
		int    result_disparity_step =  10; // bins


		String mask = ".*-DSI_COMBO\\.tiff";

		GenericDialog gd = new GenericDialog("Select file mask and histogram parameters");
		gd.addStringField ("Combined DSI file mask: ",                       mask,                        40);
		gd.addNumericField("Number of disparity bins",                       disparity_bins,               0);
		gd.addNumericField("Number of strength bins",                        strength_bins,                0);
		gd.addNumericField("Drop tiles with disparities below",              disparity_min_drop,           3);
		gd.addNumericField("Clip low disparities with",                      disparity_min_clip,           3);
		gd.addNumericField("Drop tiles with disparities above",              disparity_max_drop,           3);
		gd.addNumericField("Clip high disparities with",                     disparity_max_clip,           3);
		gd.addNumericField("Drop tiles with strength below",                 strength_min_drop,            3);
		gd.addNumericField("Clip low strength with",                         strength_min_clip,            3);
		gd.addNumericField("Drop tiles with strength above",                 strength_max_drop,            3);
		gd.addNumericField("Clip high strength with",                        strength_max_clip,            3);
		gd.addCheckbox("Normalize histogram to average 1.0",                 normalize);
		gd.addNumericField("Master weight power (after floor)",              master_weight_power,          3);
		gd.addNumericField("Master weight floor",                            master_weight_floor,          3);
		gd.addNumericField("Ignore tiles with disparity difference higher",  disparity_outlier,            3);

		gd.addNumericField("Mask: add before log to avoid -infinity",        pre_log_offs,                 3);
		gd.addNumericField("Mask: blur logarithm of the histogram",          log_sigma,                    3);
		gd.addNumericField("Mask: threshold (relative tile population)",     mask_threshold,               3);

		gd.addNumericField("Blur logarithm of the histogram in disparity pixels", log_sigma_d,             3);
		gd.addNumericField("Blur logarithm of the histogram in strength units",   log_sigma_s,             3);
		gd.addNumericField("Report RMS for each disaprity bins",             result_disparity_step,        0);

		gd.showDialog ();
		if (gd.wasCanceled()) return false;
		mask =             gd.getNextString();
		disparity_bins = (int) gd.getNextNumber();
		strength_bins =  (int) gd.getNextNumber();
		disparity_min_drop =   gd.getNextNumber();
		disparity_min_clip =   gd.getNextNumber();
		disparity_max_drop =   gd.getNextNumber();
		disparity_max_clip =   gd.getNextNumber();
		strength_min_drop =    gd.getNextNumber();
		strength_min_clip =    gd.getNextNumber();
		strength_max_drop =    gd.getNextNumber();
		strength_max_clip =    gd.getNextNumber();
		normalize =            gd.getNextBoolean();
		master_weight_power =  gd.getNextNumber();
		master_weight_floor =  gd.getNextNumber();
		disparity_outlier =    gd.getNextNumber();

		pre_log_offs =         gd.getNextNumber();
		log_sigma =            gd.getNextNumber();
		mask_threshold =       gd.getNextNumber();
		log_sigma_d =          gd.getNextNumber();
		log_sigma_s =          gd.getNextNumber();
		result_disparity_step = (int) gd.getNextNumber();


		// get list of all files:
		System.out.println("File mask = "+mask);

		final List<Path> files=new ArrayList<>();
		final String fMask = mask;
		try {
			Files.walkFileTree(path, EnumSet.of(FileVisitOption.FOLLOW_LINKS), Integer.MAX_VALUE, new SimpleFileVisitor<Path>(){
				@Override
				public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
					if(!attrs.isDirectory()){
						if (file.toString().matches(fMask)) {
							files.add(file);
						}
					}
					return FileVisitResult.CONTINUE;
				}
			});
		} catch (IOException e) {
			e.printStackTrace();
		}
		int [][] hist = new int [disparity_bins][strength_bins];
		int [] slices = {TwoQuadCLT.DSI_DISPARITY_RIG,TwoQuadCLT.DSI_STRENGTH_RIG, TwoQuadCLT.DSI_DISPARITY_MAIN,TwoQuadCLT.DSI_STRENGTH_MAIN};
		double [][][] ds_error = new double [disparity_bins][strength_bins][4]; // adding averaged with neighbors
		double [] dir_weights = {0.7,0.5}; // center weight = 1.0
		double disparity_outlier2 = disparity_outlier*disparity_outlier;
		double disparity_step = (disparity_max_clip - disparity_min_clip) / disparity_bins;
		double strength_step =  (strength_max_clip -  strength_min_clip)  / strength_bins;
		double disparity_offs = disparity_min_clip - disparity_step/2; // last and first bin that include clip will be 0.5 width
		double strength_offs =  strength_min_clip -  strength_step/2; // last and first bin that include clip will be 0.5 width
		int total_tiles_used = 0;
		int nfile = 0;
		for (Path p:files) {
			ImagePlus imp_dsi=new ImagePlus(p.normalize().toString());
			  ImageStack dsi_stack=  imp_dsi.getStack();
			  int width = dsi_stack.getWidth();
			  int height = dsi_stack.getHeight();
			  TileNeibs         tnImage = new TileNeibs(width,height);
			  float [][] dsi_float = new float [slices.length][];
			  int nLayers = dsi_stack.getSize();
			  for (int nl = 0; nl < nLayers; nl++){
				  for (int i = 0; i < slices.length; i++) {
					  if (TwoQuadCLT.DSI_SLICES[slices[i]].equals(dsi_stack.getSliceLabel(nl + 1))) {
						  dsi_float[i]= (float[]) dsi_stack.getPixels(nl + 1);
						  break;
					  }
				  }
			  }
			  int nut = 0;
			  for (int nTile = 0; nTile < dsi_float[0].length; nTile++) if (!Float.isNaN( dsi_float[0][nTile])){
				  double d = dsi_float[0][nTile];
				  double s = dsi_float[1][nTile];
				  if (
						  (d >= disparity_min_drop) &&
						  (d <= disparity_max_drop) &&
						  (s >= strength_min_drop) &&
						  (s <= strength_max_drop)) {
					  if (d < disparity_min_clip) d = disparity_min_clip;
					  if (d > disparity_max_clip) d = disparity_max_clip;
					  if (s < strength_min_clip)  s = strength_min_clip;
					  if (s > strength_max_clip)  s = strength_max_clip;
					  int dbin = (int) ((d - disparity_offs)/disparity_step);
					  if (dbin >= disparity_bins) dbin = disparity_bins - 1;
					  int sbin = (int) ((s - strength_offs)/strength_step);
					  if (sbin >= strength_bins) sbin = strength_bins - 1;
//					int [][] hist = new int [disparity_bins][strength_bins];
					  hist[dbin][sbin]++;
					  nut++;
					  double dm = dsi_float[2][nTile];
					  double sm = dsi_float[3][nTile] - master_weight_floor;
					  double de2 = (dm - d) * (dm - d);
					  if ((de2 <= disparity_outlier2) && (sm > 0.0)) {
						  double w = 1.0;
						  if (master_weight_power > 0.0) {
							  w = sm;
							  if (master_weight_power != 1.0) {
								w = Math.pow(w, master_weight_power);
							  }
						  }
						  ds_error[dbin][sbin][0] += w* de2;
						  ds_error[dbin][sbin][1] += w;
						  // combine with neighbors
						  double sw =  w;
						  double sew = w * (dm - d);
						  for (int direction = 0; direction < 8; direction++) {
							  int nTile1 = tnImage.getNeibIndex(nTile, direction);
							  if (nTile1 >= 0) {
								  d =  dsi_float[0][nTile1];
								  s =  dsi_float[1][nTile1] - strength_min_drop;
								  if ((s > 0.0) && (s <= strength_max_drop)) { // s was infinity here!
									  sm = dsi_float[3][nTile1] - master_weight_floor;
									  if (sm > 0.0) {
										  double de = dsi_float[2][nTile1] - d;
										  de2 = de*de;
										  if (de2 < disparity_outlier2) {
											  w = 1.0;
											  if (master_weight_power > 0.0) {
												  w = sm;
												  if (master_weight_power != 1.0) {
													w = Math.pow(w, master_weight_power);
												  }
											  }
											  if (! ((w >0.0) && (w < 1.0))){
												  System.out.println("strange w="+w);
											  }
											  w *= dir_weights[direction & 1] * s; //
											  sw  += w;
											  sew += w * de;
											  if (Double.isInfinite(w)){
												  System.out.println(Double.isInfinite(w));
											  }

										  }
									  }
								  }
							  }
						  }
//						  sew /= sw;
						  if (sw > 0.0) {
							  ds_error[dbin][sbin][2] += sew * sew / sw;
							  ds_error[dbin][sbin][3] += sw;
							  if (Double.isNaN(ds_error[dbin][sbin][2])) {
								  System.out.println("Double.isNaN(ds_error[dbin][sbin][2])");
							  }
						  }
					  }
				  }
			  }
			  System.out.println((++nfile)+": "+p.getFileName()+": "+nut+" useful tiles counted");
			  total_tiles_used += nut;
		}
		System.out.println("Total number of useful tiles: "+total_tiles_used+ " of "+nfile+" files");
		String [] titles = {"histogram", "histogram_masked", "histogram_ideal", "disp_err","disp_err9", "masked_err","masked_err9"};

		double [][] hist_double = new double [titles.length][disparity_bins*strength_bins];

		double scale = 1.0;
		if (normalize) {
			scale *= (1.0* disparity_bins * strength_bins) / total_tiles_used;
		}

		for (int nTile = 0; nTile < hist_double[0].length; nTile++) {
			int dbin = nTile % disparity_bins;
			int sbin = nTile / disparity_bins;
			hist_double[0][nTile] = scale * hist[dbin][sbin];
		}
		// create a mask of relatively frequent d/s cells
//		double pre_log_offs =         0.01; // add before log to avoid -infinity
//		double log_sigma =            2.00; // blur logarithm of the histogram
//		double mask_threshold =       0.25; // relative tile population
		double [] mask_calc = new double [disparity_bins*strength_bins]; // hist_double[0].clone();
		for (int nTile = 0; nTile < mask_calc.length; nTile++ ) {
			mask_calc[nTile] = Math.log( hist_double[0][nTile] + pre_log_offs);
		}
		if (log_sigma > 0.0) {
			(new DoubleGaussianBlur()).blurDouble(
					mask_calc,
					disparity_bins, // int width,
					strength_bins, // int height,
					log_sigma, // double sigmaX,
					log_sigma, // double sigmaY,
					0.01); // double accuracy)
		}
		double log_threshold = Math.log(mask_threshold+pre_log_offs);
		boolean [] ds_mask = new boolean [disparity_bins*strength_bins];
		for (int nTile = 0; nTile < mask_calc.length; nTile++ ) {
			ds_mask[nTile] = mask_calc[nTile] >= log_threshold;
		}
//		log_sigma_d =          gd.getNextNumber();
//		log_sigma_s =          gd.getNextNumber();
//		double disparity_step = (disparity_max_clip - disparity_min_clip) / disparity_bins;
//		double strength_step =  (strength_max_clip -  strength_min_clip)  / strength_bins;
		(new DoubleGaussianBlur()).blurDouble(
				mask_calc,
				disparity_bins,             // int width,
				strength_bins,              // int height,
				log_sigma_d/disparity_step, // double sigmaX,
				log_sigma_s/disparity_step, // double sigmaY,
				0.01);                      // double accuracy)
		for (int nTile = 0; nTile < mask_calc.length; nTile++ ) {
			if (!ds_mask[nTile]) {
				mask_calc[nTile] = 0.0;
			} else {
				mask_calc[nTile] = Math.exp(mask_calc[nTile]) - log_threshold;
			}
		}


		for (int nTile = 0; nTile < hist_double[0].length; nTile++) {
			int dbin = nTile % disparity_bins;
			int sbin = nTile / disparity_bins;

			hist_double[1][nTile] =  ds_mask[nTile] ? hist_double[0][nTile]: 0.0;

			hist_double[2][nTile] =  mask_calc[nTile];

			if (ds_error[dbin][sbin][2] > 0.0) {
				hist_double[3][nTile] = Math.sqrt(ds_error[dbin][sbin][0]/ds_error[dbin][sbin][1]);
			} else {
				hist_double[3][nTile] = Double.NaN;
			}


			if (ds_error[dbin][sbin][3] > 0.0) {
				hist_double[4][nTile] = Math.sqrt(ds_error[dbin][sbin][2]/ds_error[dbin][sbin][3]);
			} else {
				hist_double[4][nTile] = Double.NaN;
			}


			if (ds_mask[nTile] && (ds_error[dbin][sbin][2] > 0.0)) {
				hist_double[5][nTile] = Math.sqrt(ds_error[dbin][sbin][0]/ds_error[dbin][sbin][1]);
			} else {
				hist_double[5][nTile] = Double.NaN;
			}

			if (ds_mask[nTile] && (ds_error[dbin][sbin][3] > 0.0)) {
				hist_double[6][nTile] = Math.sqrt(ds_error[dbin][sbin][2]/ds_error[dbin][sbin][3]);
			} else {
				hist_double[6][nTile] = Double.NaN;
			}
		}

//		result_disparity_step
		// Calculate weighted by both master and strength bin (rig strength), report for each  result_disparity_step disparity bins (running totals)
		double sew=0.0, sw=0.0, sew9 = 0.0, sw9 = 0.0;
		for (int disp0 = 0; disp0 < disparity_bins; disp0 += result_disparity_step) {
			int dbin = disp0;
			for (; (dbin < (disp0 + result_disparity_step)) && (dbin < disparity_bins); dbin ++) {
				for (int sbin = 0; sbin < strength_bins; sbin++) {
					int nTile = dbin+sbin*disparity_bins;
					if (ds_mask[nTile]) {
						sew +=  ds_error[dbin][sbin][0] * sbin; //  * strength_step);
						sw +=   ds_error[dbin][sbin][1] * sbin;
						sew9 += ds_error[dbin][sbin][2] * sbin;
						sw9 +=  ds_error[dbin][sbin][3] * sbin;
						if (Double.isNaN(sew9)) {
							System.out.println("Double.isNaN(sew9)");
						}
					}
				}
			}
			double run_disp = disparity_offs + disparity_step * dbin;
			double rms = Math.sqrt(sew/sw);
			double rms9 = Math.sqrt(sew9/sw9);
			System.out.println(String.format("disparity: %7.3f pix rms= %5.3f  rms9= %5.3f", run_disp, rms, rms9));
		}

		ImagePlus imp = (new ShowDoubleFloatArrays()).makeArrays(
				hist_double,
				disparity_bins,
				strength_bins,
				"DSI_metrics",
				titles
				);
		imp.setProperty("disparity_bins",  disparity_bins+"");
		imp.setProperty("comment_disparity_bins",  "Number of disparity bins");
		imp.setProperty("strength_bins",  strength_bins+"");
		imp.setProperty("comment_strength_bins",  "Number of strength (confidence) bins");

		imp.setProperty("disparity_min_drop",  disparity_min_drop+"");
		imp.setProperty("comment_disparity_min_drop",  "Drop tiles with disparities below");
		imp.setProperty("disparity_min_clip",  disparity_min_clip+"");
		imp.setProperty("comment_disparity_min_clip",  "Clip low disparities with");
		imp.setProperty("disparity_max_drop",  disparity_max_drop+"");
		imp.setProperty("comment_disparity_max_drop",  "Drop tiles with disparities above");
		imp.setProperty("disparity_max_clip",  disparity_max_clip+"");
		imp.setProperty("comment_disparity_max_clip",  "Clip high disparities with");

		imp.setProperty("strength_min_drop",  strength_min_drop+"");
		imp.setProperty("comment_strength_min_drop",  "Drop tiles with strength below");
		imp.setProperty("strength_min_clip",  strength_min_clip+"");
		imp.setProperty("comment_strength_min_clip",  "Clip low strength with");
		imp.setProperty("strength_max_drop",  strength_max_drop+"");
		imp.setProperty("comment_strength_max_drop",  "Drop tiles with strength above");
		imp.setProperty("strength_max_clip",  strength_max_clip+"");
		imp.setProperty("comment_strength_max_clip",  "Clip high strength with");

		imp.setProperty("total_tiles_used",  total_tiles_used+"");
		imp.setProperty("comment_total_tiles_used",  "Total number of tiles used");

		imp.setProperty("master_weight_power",  master_weight_power+"");
		imp.setProperty("comment_master_weight_power",  "Master weight power (after floor)");

		imp.setProperty("master_weight_floor",  master_weight_floor+"");
		imp.setProperty("comment_master_weight_floor",  "Master weight floor");

		imp.setProperty("disparity_outlier",  disparity_outlier+"");
		imp.setProperty("comment_disparity_outlier",  "Ignore tiles with disparity difference higher");

		(new JP46_Reader_camera(false)).encodeProperiesToInfo(imp);
		imp.show();
		return true;

	}

	public static boolean mlRecalc(String dir,
			EyesisCorrectionParameters.CLTParameters       clt_parameters,
			EyesisCorrectionParameters.DebayerParameters   debayerParameters,
			EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			EyesisCorrectionParameters.RGBParameters       rgbParameters,
			final int                                      threadsMax,  // maximal number of threads to launch
			final boolean                                  updateStatus,
			final int                                      debugLevel) throws Exception

	{
		String mask = ".*EXTRINSICS\\.corr-xml";
		String full_conf_suffix = ".corr-xml";
		String dsi_suffix = "-DSI_COMBO.tiff";
		String correction_parameters_prefix = "CORRECTION_PARAMETERS.";

		System.out.println("File mask = "+mask);

		final List<Path> files=new ArrayList<>();
		final String fMask = mask;
		Path path= Paths.get(dir);
		try {
			Files.walkFileTree(path, EnumSet.of(FileVisitOption.FOLLOW_LINKS), Integer.MAX_VALUE, new SimpleFileVisitor<Path>(){
				@Override
				public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
					if(!attrs.isDirectory()){
						if (file.toString().matches(fMask)) {
							files.add(file);
						}
					}
					return FileVisitResult.CONTINUE;
				}
			});
		} catch (IOException e) {
			e.printStackTrace();
			return false;
		}
		int indx = 0;
		for (Path p:files) {

			int count =     p.getNameCount();
			if (count >=3) {
				String model =   p.getName(count-3).toString();
				String version = p.getName(count-2).toString();
				String name =    p.getName(count-1).toString();
//				String root = p.getRoot().toString();
				Path parent = p.getParent();
//				System.out.println(indx+": model:"+model+", version:"+version+", name: "+name);
				Path full_conf_path = p.resolveSibling(model+full_conf_suffix);
				Path dsi_combo_path = p.resolveSibling(model+dsi_suffix);
				System.out.println(indx+": model:"+model+", version:"+version+", name: "+name+ ", parent="+parent+", full_conf_path="+full_conf_path);
				// See if there is a full configuration file and DSI combo
				if (!full_conf_path.toFile().exists()) {
					System.out.println("Full configuration file for model:"+model+", version:"+version+": "+full_conf_path+" does not exist");
					continue;
				}
				if (!dsi_combo_path.toFile().exists()) {
					System.out.println("DSI combo file (GT data) for model:"+model+", version:"+version+": "+dsi_combo_path+" does not exist");
					continue;
				}
				Properties full_properties = loadProperties(
						full_conf_path.toString(), // String path,
						null); // Properties properties)
//				String [] sourcePaths = null;
				ArrayList<String> source_list = new ArrayList<String>();
				if (full_properties.getProperty(correction_parameters_prefix+"sourcePaths")!=   null){
					int numFiles=Integer.parseInt(full_properties.getProperty(correction_parameters_prefix+"sourcePaths"));
//					sourcePaths=new String[numFiles];
					for (int i = 0; i < numFiles; i++){
						String src=full_properties.getProperty(correction_parameters_prefix+"sourcePath"+i);
						Path src_path=Paths.get(src);
						if (src_path.getName(src_path.getNameCount()-1).toString().contains(model)){
							source_list.add(src);
						}
	        		}
				} else {
					System.out.println("Source files are not povided for model:"+model+", version:"+version+" in "+full_conf_path);
					continue;
				}
				System.out.println("Got "+ source_list.size() +" source files");

				if (indx > -1) {
					return true; // temporarily
				}

				/*
				 *
				this.sourcePaths=new String[numFiles];
				for (int i=0;i<this.sourcePaths.length;i++){
					this.sourcePaths[i]=properties.getProperty(prefix+"sourcePath"+i);
        		}

				 *
				 *
				Properties properties = loadProperties(
						p.toString(), // String path,
						null); // Properties properties)
				QuadCLT qcm = new QuadCLT(
						QuadCLT.PREFIX, // String                                          prefix,
						properties, // Properties                                      properties,
						null, // EyesisCorrections                               eyesisCorrections,
						null // EyesisCorrectionParameters.CorrectionParameters correctionsParameters
						);
				QuadCLT qca = new QuadCLT(
						QuadCLT.PREFIX_AUX, // String                                          prefix,
						properties, // Properties                                      properties,
						null, // EyesisCorrections                               eyesisCorrections,
						null // EyesisCorrectionParameters.CorrectionParameters correctionsParameters
						);
				System.out.println(indx+": model:"+model+", version:"+version+", name: "+name);
				GeometryCorrection.CorrVector cvm = qcm.geometryCorrection.getCorrVector();
				GeometryCorrection.CorrVector cva = qca.geometryCorrection.getCorrVector();
				*/
				indx++;
			}
		}



		return true;
	}


	public static boolean listExtrinsics(String dir) // , String mask)
	{
		Path path= Paths.get(dir);
		boolean inPixels =        true;
		boolean inMrad =          false;
		boolean showATR =         true;
		boolean showZooms =       true;
		boolean showSym =         true;
		boolean showRigATR =      true;
		boolean showRigZoom =     true;
		boolean showRigAngle =    true;
		boolean showRigBaseline = false;
		String mask = ".*EXTRINSICS\\.corr-xml";

		GenericDialog gd = new GenericDialog("Select file mask and output format");
		gd.addStringField ("Extrinsics (or general configuration) file mask: ", mask, 40);
		gd.addCheckbox("Show results in pixels (false in angular units)",inPixels);
		gd.addCheckbox("Show results in mrad (false in arcseconds)",     inMrad);
		gd.addCheckbox("Show azimuths, tilts, rolls",                    showATR);
		gd.addCheckbox("Show zooms",                                     showZooms);
		gd.addCheckbox("Show symmetric angles",                          showSym);
		gd.addCheckbox("Show rig azimuth, tilt, roll",                   showRigATR);
		gd.addCheckbox("Show rig zoom",                                  showRigZoom);
		gd.addCheckbox("Show rig angle",                                 showRigAngle);
		gd.addCheckbox("Show rig baseline",                              showRigBaseline);
		gd.showDialog ();
		if (gd.wasCanceled()) return false;
		mask =             gd.getNextString();
		inPixels =         gd.getNextBoolean();
		inMrad =           gd.getNextBoolean();
		showATR =          gd.getNextBoolean();
		showZooms =        gd.getNextBoolean();
		showSym =          gd.getNextBoolean();
		showRigATR =       gd.getNextBoolean();
		showRigZoom =      gd.getNextBoolean();
		showRigAngle =     gd.getNextBoolean();
		showRigBaseline =  gd.getNextBoolean();

		String units =  inPixels ? "pix":(inMrad?"mil":"\"");
		String zunits = inPixels ? "pix":(inMrad?"mil":"\"");
		double scale = inPixels ? 1.0 : (inMrad?1000.0:(180.0/Math.PI*60*60)); //leave pixels as is, convert radians to arc-sec
		String fmt =  "\t"+(inPixels ? "%8.4f":(inMrad?"%8.4f":"%8.2f"));
		String fmt_angle =  "\t%8.3f";
		String fmt_len =    "\t%8.3f";
		int num_sym = showZooms? GeometryCorrection.CorrVector.LENGTH:GeometryCorrection.CorrVector.ZOOM_INDEX;
		class ModVerString{
			String model;
			String version;
			String txt;
			ModVerString (String model, String version, String txt){
				this.model = model;
				this.version = version;
				this.txt = txt;
			}
			@Override
			public String toString() {
				return model+"\t"+version+txt;
			}

		}

		ArrayList<ModVerString> line_list = new ArrayList<ModVerString>();
		String title = "Extrinsic_parameters_variations_in_"+(inPixels ? "pixels":(inMrad?"mrad":"arcseconds"));
		String header="#\tModel\tVersion";
		int num_col = 0;
		if (showATR) {
			header+=String.format("\taz m0 (%s)\taz m1 (%s)\taz m2 (%s)\taz m3 (%s)"+
					"\ttl m0 (%s)\ttl  m1 (%s) \ttl m2 (%s)\ttl m3 (%s)"+
					"\trl m0 (%s)\trl  m1 (%s) \trl m2 (%s)\trl m3 (%s)",
					units,units,units,units,units,units,units,units,units,units,units,units);
			num_col+=12;
		}
		if (showZooms) {
			header+=String.format("\tzm m0 (%s)\tzm m1 (%s)\tzm m2 (%s)\tzm m3 (%s)", zunits, zunits, zunits, zunits);
			num_col+=4;
		}

		if (showSym) {
			for (int i = 0; i < num_sym; i++) {
				header+=String.format("\tsym%02d-m", i);
			}
			num_col+=num_sym;

		}

		if (showATR) {
			header+=String.format("\taz a0 (%s)\taz a1 (%s)\taz a2 (%s)\taz a3 (%s)"+
					"\ttl a0 (%s)\ttl  a1 (%s) \ttl a2 (%s)\ttl a3 (%s)"+
					"\trl a0 (%s)\trl  a1 (%s) \trl a2 (%s)\trl a3 (%s)",
					units,units,units,units,units,units,units,units,units,units,units,units);
			num_col+=12;
		}
		if (showZooms) {
			header+=String.format("\tzm a0 (%s)\tzm a1 (%s)\tzm a2 (%s)\tzm a3 (%s)", zunits, zunits, zunits, zunits);
			num_col+=4;
		}
		if (showSym) {
			for (int i = 0; i < num_sym; i++) {
				header+=String.format("\tsym%02d-a", i);
			}
			num_col+=num_sym;
		}
		if (showRigATR) {
			header+=String.format("\trig azmth (%s)\trig tilt(%s)\trig roll (%s)", units, units, units);
			num_col+=3;
		}
		if (showRigZoom) {
			header+=String.format("\trig zoom (%s)", zunits);
			num_col+=1;
		}
		//
		if (showRigAngle) {
			header+="\trig angle (°)";
			num_col+=1;
		}
		if (showRigBaseline) {
			header+="\trig baseline (mm)";
			num_col+=1;
		}


		System.out.println("File mask = "+mask);

		final List<Path> files=new ArrayList<>();
		final String fMask = mask;
		try {
			Files.walkFileTree(path, EnumSet.of(FileVisitOption.FOLLOW_LINKS), Integer.MAX_VALUE, new SimpleFileVisitor<Path>(){
				@Override
				public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
					if(!attrs.isDirectory()){
						if (file.toString().matches(fMask)) {
							files.add(file);
						}
					}
					return FileVisitResult.CONTINUE;
				}
			});
		} catch (IOException e) {
			e.printStackTrace();
		}
		int indx=1;
		double [][] stats = new double [num_col][2];
		String [] fmts = new String[num_col];
		for (Path p:files) {

			int count =     p.getNameCount();
			if (count >=3) {
				String model =   p.getName(count-3).toString();
				String version = p.getName(count-2).toString();
				String name =    p.getName(count-1).toString();
				System.out.println(indx+": model:"+model+", version:"+version+", name: "+name);
				Properties properties = loadProperties(
						p.toString(), // String path,
						null); // Properties properties)
				QuadCLT qcm = new QuadCLT(
						QuadCLT.PREFIX, // String                                          prefix,
						properties, // Properties                                      properties,
						null, // EyesisCorrections                               eyesisCorrections,
						null // EyesisCorrectionParameters.CorrectionParameters correctionsParameters
						);
				QuadCLT qca = new QuadCLT(
						QuadCLT.PREFIX_AUX, // String                                          prefix,
						properties, // Properties                                      properties,
						null, // EyesisCorrections                               eyesisCorrections,
						null // EyesisCorrectionParameters.CorrectionParameters correctionsParameters
						);
				System.out.println(indx+": model:"+model+", version:"+version+", name: "+name);
				GeometryCorrection.CorrVector cvm = qcm.geometryCorrection.getCorrVector();
				GeometryCorrection.CorrVector cva = qca.geometryCorrection.getCorrVector();
//				double [] vect_main = qcm.geometryCorrection.getCorrVector().toArray();
				StringBuffer sb = new StringBuffer();
				int ncol = 0;

//				sb.append(indx+"\t"+model+"\t"+version);
				if (showATR) { // main camera
					double [] v = new double[4];
					for (int i = 0; i <4; i++) {
						if (i < 3) {
							v[i] =  scale * cvm.getExtrinsicParameterValue(i+GeometryCorrection.CorrVector.AZIMUTH_INDEX, inPixels);
							v[3] -= v[i];
						}
						sb.append(String.format(fmt,v[i])); // azimuths
						fmts [ncol] = fmt;
						stats[ncol  ][0]+=v[i];
						stats[ncol++][1]+=v[i]*v[i];
					}
					v = new double[4];
					for (int i = 0; i <4; i++) {
						if (i < 3) {
							v[i] =  scale * cvm.getExtrinsicParameterValue(i+GeometryCorrection.CorrVector.TILT_INDEX, inPixels);
							v[3] -= v[i];
						}
						sb.append(String.format(fmt,v[i])); // tilts
						fmts [ncol] = fmt;
						stats[ncol  ][0]+=v[i];
						stats[ncol++][1]+=v[i]*v[i];
					}
					v = new double[4];
					for (int i = 0; i <4; i++) {
						if (i < 3) {
							v[i] =  scale * cvm.getExtrinsicParameterValue(i+GeometryCorrection.CorrVector.ROLL_INDEX, inPixels);
							v[3] -= v[i];
						}
						sb.append(String.format(fmt,v[i])); // rolls
						fmts [ncol] = fmt;
						stats[ncol  ][0]+=v[i];
						stats[ncol++][1]+=v[i]*v[i];
					}
				}
				if (showZooms) { // main camera
					double [] v = new double[4];
					for (int i = 0; i <4; i++) {
						if (i < 3) {
							v[i] =  scale * cvm.getExtrinsicParameterValue(i+GeometryCorrection.CorrVector.ZOOM_INDEX, inPixels);
							v[3] -= v[i];
						}
						sb.append(String.format(fmt,v[i])); // zooms
						fmts [ncol] = fmt;
						stats[ncol  ][0]+=v[i];
						stats[ncol++][1]+=v[i]*v[i];
					}
				}

				if (showSym) {
					for (int i = 0; i < num_sym; i++) {
						double v =  scale * cvm.getExtrinsicSymParameterValue(i, inPixels);
						sb.append(String.format(fmt,v)); // sym parameters
						fmts [ncol] = fmt;
						stats[ncol  ][0]+=v;
						stats[ncol++][1]+=v*v;
					}
				}


				if (showATR) { // aux camera
					double [] v = new double[4];
					for (int i = 0; i <4; i++) {
						if (i < 3) {
							v[i] =  scale * cva.getExtrinsicParameterValue(i+GeometryCorrection.CorrVector.AZIMUTH_INDEX, inPixels);
							v[3] -= v[i];
						}
						sb.append(String.format(fmt,v[i])); // azimuths
						fmts [ncol] = fmt;
						stats[ncol  ][0]+=v[i];
						stats[ncol++][1]+=v[i]*v[i];
					}
					v = new double[4];
					for (int i = 0; i <4; i++) {
						if (i < 3) {
							v[i] =  scale * cva.getExtrinsicParameterValue(i+GeometryCorrection.CorrVector.TILT_INDEX, inPixels);
							v[3] -= v[i];
						}
						sb.append(String.format(fmt,v[i])); // tilts
						fmts [ncol] = fmt;
						stats[ncol  ][0]+=v[i];
						stats[ncol++][1]+=v[i]*v[i];
					}
					v = new double[4];
					for (int i = 0; i <4; i++) {
						if (i < 3) {
							v[i] =  scale * cva.getExtrinsicParameterValue(i+GeometryCorrection.CorrVector.ROLL_INDEX, inPixels);
							v[3] -= v[i];
						}
						sb.append(String.format(fmt,v[i])); // rolls
						fmts [ncol] = fmt;
						stats[ncol  ][0]+=v[i];
						stats[ncol++][1]+=v[i]*v[i];
					}
				}
				if (showZooms) { // aux camera
					double [] v = new double[4];
					for (int i = 0; i <4; i++) {
						if (i < 3) {
							v[i] =  scale * cva.getExtrinsicParameterValue(i+GeometryCorrection.CorrVector.ZOOM_INDEX, inPixels);
							v[3] -= v[i];
						}
						sb.append(String.format(fmt,v[i])); // zooms
						fmts [ncol] = fmt;
						stats[ncol  ][0]+=v[i];
						stats[ncol++][1]+=v[i]*v[i];
					}
				}
				if (showSym) {
					for (int i = 0; i < num_sym; i++) {
						double v =  scale * cva.getExtrinsicSymParameterValue(i, inPixels);
						sb.append(String.format(fmt,v)); // sym parameters
						fmts [ncol] = fmt;
						stats[ncol  ][0]+=v;
						stats[ncol++][1]+=v*v;
					}
				}
				if (showRigATR) {
					int [] indices = {
							GeometryCorrection.RigOffset.AUX_AZIMUTH_INDEX,
							GeometryCorrection.RigOffset.AUX_TILT_INDEX,
							GeometryCorrection.RigOffset.AUX_ROLL_INDEX};
					for (int i = 0; i < indices.length; i++) {
						double v = scale * qca.geometryCorrection.getRigOffsetParameter(indices[i], inPixels);
						sb.append(String.format(fmt,v)); // tig atz parameters
						fmts [ncol] = fmt;
						stats[ncol  ][0]+=v;
						stats[ncol++][1]+=v*v;
					}
				}
				if (showRigZoom) {
					double v = scale * qca.geometryCorrection.getRigOffsetParameter(GeometryCorrection.RigOffset.AUX_ZOOM_INDEX, inPixels);
					sb.append(String.format(fmt,v)); // tig atz parameters
					fmts [ncol] = fmt;
					stats[ncol  ][0]+=v;
					stats[ncol++][1]+=v*v;
				}

				if (showRigAngle) {
//					sb.append(String.format(fmt_angle, 180/Math.PI * qca.geometryCorrection.getRigOffsetParameter(GeometryCorrection.RigOffset.AUX_ANGLE_INDEX,    false)));
					double v = 180/Math.PI * qca.geometryCorrection.getRigOffsetParameter(GeometryCorrection.RigOffset.AUX_ANGLE_INDEX, inPixels);
					sb.append(String.format(fmt_angle,v)); // tig atz parameters
					fmts [ncol] = fmt_angle;
					stats[ncol  ][0]+=v;
					stats[ncol++][1]+=v*v;

				}
				if (showRigBaseline) {
//					sb.append(String.format(fmt_len,                 qca.geometryCorrection.getRigOffsetParameter(GeometryCorrection.RigOffset.AUX_BASELINE_INDEX,    false)));
					double v = qca.geometryCorrection.getRigOffsetParameter(GeometryCorrection.RigOffset.AUX_ANGLE_INDEX, inPixels);
					sb.append(String.format(fmt_len,v)); // tig atz parameters
					fmts [ncol] = fmt_len;
					stats[ncol  ][0]+=v;
					stats[ncol++][1]+=v*v;
				}
				sb.append("\n");
				line_list.add( new ModVerString(model,version,sb.toString()));
				indx++;
			}
		}
		Collections.sort(line_list, new Comparator<ModVerString>() {
			@Override
			public int compare(ModVerString lhs, ModVerString rhs) {
				// -1 - less than, 1 - greater than, 0 - equal, not inverted for ascending disparity
				int rslt = lhs.model.compareTo(rhs.model);
				if (rslt == 0) rslt = lhs.version.compareTo(rhs.version);
				return rslt;
			}
		});

		StringBuffer sb_avg = new StringBuffer();
		StringBuffer sb_rms = new StringBuffer();
		sb_avg.append("--\tAverage\t");
		sb_rms.append("--\tStandard deviation\t");
		int nrows = indx -1;
		if (indx > 1) {
			for (int ncol = 0; ncol < num_col; ncol++) {
				stats[ncol][0] /= nrows;
				stats[ncol][1] /= nrows;
				stats[ncol][1] = Math.sqrt(stats[ncol][1] - stats[ncol][0]*stats[ncol][0]);
				if (Double.isNaN(stats[ncol][1])) { // small negative by rounding error?
					stats[ncol][1] = 0.0;
				}
				sb_avg.append(String.format(fmts[ncol],stats[ncol][0]));
				sb_rms.append(String.format(fmts[ncol],stats[ncol][1]));
			}
		}
		// ~1280

		//		 int a = GeometryCorrection.RigOffset.VECTOR_LENGTH;
		StringBuffer sb = new StringBuffer();
		indx = 1;
		for (ModVerString mvs:line_list) {
			sb.append((indx++)+"\t"+mvs.toString());
		}
		sb.append(sb_avg.toString()+"\n");
		sb.append(sb_rms.toString()+"\n");
	    new TextWindow (title, header, sb.toString(), 1200,800);
		return true;
	}
	public static Properties loadProperties(
			String path,
			Properties properties){
		if (properties == null) {
			properties = new Properties();
		}
		InputStream is;
		try {
			is = new FileInputStream(path);
		} catch (FileNotFoundException e) {
			IJ.showMessage("Error","Failed to open configuration file: "+path);
			return null;
		}
		try {
			properties.loadFromXML(is);

		} catch (IOException e) {
			IJ.showMessage("Error","Failed to read XML configuration file: "+path);
			return null;
		}
		try {
			is.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return properties;
		//	     getAllProperties(properties);
		//		 if (DEBUG_LEVEL>0) System.out.println("Configuration parameters are restored from "+path);
	}


}
