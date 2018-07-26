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
		int disparity_bins =         400;
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

		String mask = ".*-DSI_COMBO\\.tiff";

		GenericDialog gd = new GenericDialog("Select file mask and histogram parameters");
		gd.addStringField ("Combined DSI file mask: ",           mask, 40);
		gd.addNumericField("Number of disparity bins",           disparity_bins,               0);
		gd.addNumericField("Number of strength bins",            strength_bins,                0);
		gd.addNumericField("Drop tiles with disparities below",  disparity_min_drop,           3);
		gd.addNumericField("Clip low disparities with",          disparity_min_clip,           3);
		gd.addNumericField("Drop tiles with disparities above",  disparity_max_drop,           3);
		gd.addNumericField("Clip high disparities with",         disparity_max_clip,           3);
		gd.addNumericField("Drop tiles with strength below",     strength_min_drop,            3);
		gd.addNumericField("Clip low strength with",             strength_min_clip,            3);
		gd.addNumericField("Drop tiles with strength above",     strength_max_drop,            3);
		gd.addNumericField("Clip high strength with",            strength_max_clip,            3);
		gd.addCheckbox("Normalize histogram to average 1.0",    normalize);

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
		int [] slices = {TwoQuadCLT.DSI_DISPARITY_RIG,TwoQuadCLT.DSI_STRENGTH_RIG};
		double disparity_step = (disparity_max_clip - disparity_min_clip) / disparity_bins;
		double strength_step =  (strength_max_clip -  strength_min_clip)  / strength_bins;
		double disparity_offs = disparity_min_clip - disparity_step/2; // last and first bin that include clip will be 0.5 width
		double strength_offs =  strength_min_clip -  strength_step/2; // last and first bin that include clip will be 0.5 width
		int total_tiles_used = 0;
		for (Path p:files) {
			ImagePlus imp_dsi=new ImagePlus(p.normalize().toString());
			  ImageStack dsi_stack=  imp_dsi.getStack();
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

				  }
			  }
			  System.out.println(p.getFileName()+": "+nut+" useful tiles counted");
			  total_tiles_used += nut;
		}
	   System.out.println("Total number of useful tiles: "+total_tiles_used);
		double [] hist_double = new double [disparity_bins*strength_bins];
		double scale = 1.0;
		if (normalize) {
			scale *= (1.0* disparity_bins * strength_bins) / total_tiles_used;
		}
		for (int nTile = 0; nTile < hist_double.length; nTile++) {
			int dbin = nTile % disparity_bins;
			int sbin = nTile / disparity_bins;
			hist_double[nTile] = scale * hist[dbin][sbin];
		}
//		  ImagePlus imp= makeArrays(pixels, width, height, title);
//		  if (imp!=null) imp.show();

		ImagePlus imp = (new showDoubleFloatArrays()).makeArrays(
				hist_double,
				disparity_bins,
				strength_bins,
				"DSI_histogram");
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


		(new JP46_Reader_camera(false)).encodeProperiesToInfo(imp);
		imp.show();
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
