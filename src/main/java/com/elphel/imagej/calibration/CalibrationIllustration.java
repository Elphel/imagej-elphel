package com.elphel.imagej.calibration;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Rectangle;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.calibration.DistortionCalibrationData.GridImageParameters;
import com.elphel.imagej.cameras.ThermalColor;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.correction.EyesisCorrections;
import com.elphel.imagej.lwir.LwirReaderParameters;
import com.elphel.imagej.readers.ImagejJp4Tiff;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.Line;
import ij.process.ImageProcessor;
import loci.formats.FormatException;

public class CalibrationIllustration {
	public static final int MAX_THREADS = 100; // combine from all classes?
	public static final int CONTRAST_INDEX = 2;
	public static final String FOOTAGE_DIR = "footage";

	EyesisAberrations              eyesisAberrations;
	Distortions                    distortions;
	CalibrationIllustrationParameters illustrationParameters;	
	LwirReaderParameters           lwirReaderParameters;
	AtomicInteger                  stopRequested;
	int                            debug_level;
	int                            numStations;
	int                            numSubCameras;
	String [][]                    gridFileDirs; // =       new String [numStations][];
	boolean [][][]                 gridUseChn; //  =     new boolean [numStations][][]; // channels to use in each scene  
	String []                      sourceStationDirs; // = new String [numStations];      // directories of the source files per station
	String []                      grid_extensions={".tiff"};
	String []                      src_extensions={".tiff"};
	double [][]                    offs_scale;// [chn][0] - offset to subtract to normalize, [chn][1] - scale to divide by to normalize     
	
	public CalibrationIllustration (
			LwirReaderParameters           lwirReaderParameters,
			CalibrationIllustrationParameters illustrationParameters,			
			EyesisAberrations              eyesisAberrations,
			Distortions                    distortions,
			AtomicInteger                  stopRequested,
			int                            debug_level) {
		this.lwirReaderParameters =        lwirReaderParameters;
		this.illustrationParameters =      illustrationParameters;
		this.eyesisAberrations =           eyesisAberrations;
		this.distortions =                 distortions;
		this.stopRequested =               stopRequested;
		this.debug_level =                 debug_level;
		offs_scale = new double[lwirReaderParameters.getTypeMap().length][2];
		for (int i = 0; i < offs_scale.length; i++) {
			offs_scale[i][0] = 0.0;
			offs_scale[i][1] = 1.0;
		}
		
	}			

	public void plotGrid(
			int       numImg,
			int       line_width,
			ImagePlus imp,
			Color     color_grid,
			Color     color_grid_weak,
			Color     color_grid_extra, // at least one end points to extra (unreliable) nodes (may be null)
			double    weak_threshold
			) {
		if (numImg == 2276) {
			System.out.println(">>>numImg="+numImg);
		}
		
//		int line_width =  3; //https://imagej.nih.gov/ij/developer/api/ij/ij/gui/Line.html#drawPixels(ij.process.ImageProcessor)
		DistortionCalibrationData dcd = distortions.fittingStrategy.distortionCalibrationData;
		GridImageParameters gip = dcd.gIP[numImg];
		int [][] pUV =          gip.pixelsUV;
		int [][] pUV_extra =    gip.pixelsUV_extra;
		double [][] pXY =       gip.pixelsXY;
		double [][] pXY_extra = gip.pixelsXY_extra;
		if ((pUV == null) || (pUV.length==0)) {
			return;
		}
		int minU=pUV[0][0],maxU=minU,minV=pUV[0][1],maxV=minV;
		for (int i = 0; i < pUV.length; i++) {
			if      (pUV[i][0] < minU) minU = pUV[i][0];
			else if (pUV[i][0] > maxU) maxU = pUV[i][0];
			if      (pUV[i][1] < minV) minV = pUV[i][1];
			else if (pUV[i][1] > maxV) maxV = pUV[i][1];
		}
		if (color_grid_extra != null) {
			for (int i = 0; i < pUV_extra.length; i++) {
				if      (pUV_extra[i][0] < minU) minU = pUV_extra[i][0];
				else if (pUV_extra[i][0] > maxU) maxU = pUV_extra[i][0];
				if      (pUV_extra[i][1] < minV) minV = pUV_extra[i][1];
				else if (pUV_extra[i][1] > maxV) maxV = pUV_extra[i][1];
			}
		}
		
		Rectangle dimsUV= new Rectangle(minU-1, minV-1, maxU-minU+2, maxV-minV+2);
		int [] grid = new int [dimsUV.height * dimsUV.width];
		Arrays.fill(grid,  -1);
		int main_len = pUV.length;
		for (int i = 0; i < main_len; i++) {
			int indx = (pUV[i][0]-dimsUV.x) + (pUV[i][1]-dimsUV.y)* dimsUV.width;
			grid[indx] = i;
		}
		if (color_grid_extra != null) {
			for (int i = 0; i < pUV_extra.length; i++) {
				int indx = (pUV_extra[i][0]-dimsUV.x) + (pUV_extra[i][1]-dimsUV.y)* dimsUV.width;
				grid[indx] = i + main_len;
			}
		}
		ImageProcessor ip = imp.getProcessor();
		ip.setLineWidth(line_width);
//		only use stroke for line_width>1
		/*
		BasicStroke stroke = null;
		if (line_width>1) {
			stroke = new BasicStroke(line_width, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
		}
		*/
//		BasicStroke stroke = new BasicStroke(line_width, BasicStroke.CAP_BUTT, BasicStroke.JOIN_ROUND);
//		BasicStroke stroke = new BasicStroke(line_width, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
//		roi.setStrokeColor(new Color(red, green, blue, alpha));
//		roi.setStroke(stroke);		
//		int upper = dimsUV.width * (dimsUV.height-1);
		int upper = dimsUV.width * dimsUV.height;
		for (int indx = dimsUV.width + 1; indx < upper; indx++) if (grid[indx] >= 0){
			// try to the right and down
			int g0 = grid[indx];
			if ((numImg == 2276) && (g0 == 78)) {
				System.out.println(">>>g0="+g0);
			}
			for (int dir = 0; dir <2; dir++) { // 90 - right, 1 - down
				int indx1 = indx + ((dir > 0) ? dimsUV.width : 1);
				if (indx1 < grid.length) {
					int g1 = grid[indx1];
					if (g1 >= 0) {
						double [] pXY0 =  (g0 < main_len) ? pXY[g0] : pXY_extra[g0 - main_len];
						double [] pXY1 =  (g1 < main_len) ? pXY[g1] : pXY_extra[g1 - main_len];
//						Line line = new Line(pXY0[0], pXY0[1], pXY1[0],pXY1[1]);
						if ((g0 < main_len) && (g1 < main_len)) {
							if ((pXY0[CONTRAST_INDEX] >= weak_threshold) && (pXY1[CONTRAST_INDEX] >= weak_threshold)) {
								if (color_grid == null) {
									continue;
								}
								ip.setColor(color_grid);
							} else {
								if (color_grid_weak == null) {
									continue;
								}
								if ((pXY0[CONTRAST_INDEX] == 0.0) && (pXY1[CONTRAST_INDEX] == 0.0)) { // treat strength==0 as extra
									if (color_grid_extra == null) {
										continue;
									}
									ip.setColor(color_grid_extra); // should not get here with color_grid_extra==null
								} else {
									ip.setColor(color_grid_weak);
								}
							}
						} else {
							if (color_grid_extra == null) {
								continue;
							}
							ip.setColor(color_grid_extra); // should not get here with color_grid_extra==null
						}
						/*
						if (stroke != null) {
							if ((g0 < main_len) && (g1 < main_len)) {
								line.setStrokeColor(color_grid);
							} else {
								line.setStrokeColor(color_grid_extra); // should not get here with color_grid_extra==null
							}
							line.setStroke(stroke);
						}
						line.setStrokeWidth(1.0);
						*/
						for (int dy = -line_width+1; dy < line_width; dy+=2) {
							for (int dx = -line_width+1; dx < line_width; dx+=2) {
								Line line = new Line(pXY0[0] + 0.5*dx, pXY0[1] + 0.5*dy, pXY1[0] + 0.5*dx,pXY1[1] + 0.5*dy);
								line.drawPixels(ip);
							}
						}
//						Line line = new Line(pXY0[0], pXY0[1], pXY1[0],pXY1[1]);
//						line.drawPixels(ip);
						/*
						for (double lw = 0.5; lw <= line_width; lw+=0.5) {
							line.setStrokeWidth(lw);
							line.drawPixels(ip);
						}
						*/
					}
				}
			}
		}
	}
	
	public boolean convertKernels(
			int dsize, // direct kernel size
			int rsize  // inverse kernel size
			) {
		long startTime=System.nanoTime(); // restart timer after possible interactive dialogs
		boolean same_size_inv = true;
//		final boolean [] selectedChannels = eyesisAberrations.aberrationParameters.getChannelSelection(distortions);
		final boolean [] selectedChannels =  illustrationParameters.getSelectedChannels();
//		final DistortionCalibrationData dcd = distortions.fittingStrategy.distortionCalibrationData;
//		final MultipleExtensionsFileFilter sourceFilter =
//				new MultipleExtensionsFileFilter("",src_extensions,"Source calibration images");
//		final int lwir0 =          illustrationParameters.getLwirReaderParameters().getLwirChn0();
//		final int eo0 =            illustrationParameters.getLwirReaderParameters().getEoChn0();
		int [][] whall = new int [selectedChannels.length][2];
		String illustrations_dir = eyesisAberrations.aberrationParameters.illustrationsDirectory;
		String partial_dir =       eyesisAberrations.aberrationParameters.partialKernelDirectory;
		String direct_dir =        eyesisAberrations.aberrationParameters.psfKernelDirectory;
		String inverted_dir =      eyesisAberrations.aberrationParameters.aberrationsKernelDirectory;
		String partial_dest = illustrations_dir+ (illustrations_dir.endsWith(Prefs.getFileSeparator())?"":Prefs.getFileSeparator())+"partial";
		String direct_dest =  illustrations_dir+ (illustrations_dir.endsWith(Prefs.getFileSeparator())?"":Prefs.getFileSeparator())+"direct";
		String inverted_dest =  illustrations_dir+ (illustrations_dir.endsWith(Prefs.getFileSeparator())?"":Prefs.getFileSeparator())+"inverted";

		if (illustrationParameters.kernel_process_partial) {
			convertPartial(
					selectedChannels,
					dsize, 
					partial_dir, // String source_dir,
					partial_dest);
		}
		if (illustrationParameters.kernel_process_direct) {
			convertDirect(
					selectedChannels,
					dsize, 
					direct_dir, // String source_dir,
					direct_dest,
					whall);
		}
		if (illustrationParameters.kernel_process_inverse) {
			if (!same_size_inv) {
				whall = null;
			}
			convertInverted(
					selectedChannels,
					rsize, 
					inverted_dir, // String source_dir,
					inverted_dest,
					whall);
		}
		
		System.out.println("All done in "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+" sec.");
		return true;
	}
	
	public boolean convertPartial(
			final boolean [] selectedChannels,
			int src_step,
			String source_dir,
			String dest_dir) {
		String partialSuffix= eyesisAberrations.aberrationParameters.partialSuffix; // .".ppsf-tiff";
		File destDir= new File (dest_dir);
		if (!destDir.exists()){
			if (!destDir.mkdirs()) {
				IJ.showMessage("Error","Failed to create results directory "+dest_dir);
				return false;
			}
		}
		String [][] src_files = new String [selectedChannels.length][];
		for (int nChn = 0; nChn < src_files.length; nChn++) if (selectedChannels[nChn]) {
			MultipleExtensionsFileFilter sourceFilter =
					new MultipleExtensionsFileFilter("",new String[] {String.format("%02d",  nChn)+partialSuffix},"partial kernels");
			src_files[nChn] = (new File(source_dir)).list(sourceFilter); // are these full files?
		}
		// create list of partial kernel files
		// partial-1623732618_112155-00.ppsf-tiff
   		final Thread[] threads = newThreadArray(MAX_THREADS);
   		final AtomicInteger indxAtomic = new AtomicInteger(0);
		
		for (int iChn = 0; iChn < src_files.length; iChn++) if (src_files[iChn]!=null) {
			final int nChn=iChn;
			indxAtomic.set(0);
			
			int sensor_type = distortions.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getSensorType(nChn);
			int decimate =    (sensor_type==1) ? illustrationParameters.kernel_decimate_direct_lwir : illustrationParameters.kernel_decimate_direct_eo;
			int kernel_widh = (sensor_type==1) ? illustrationParameters.kernel_dia_direct_lwir : illustrationParameters.kernel_dia_direct_eo;

			
//			String illustrations_dir = eyesisAberrations.aberrationParameters.illustrationsDirectory;
			final String chn_ill_dir = dest_dir+Prefs.getFileSeparator()+ illustrationParameters.getChannelPrefix()+String.format("%02d", nChn);
			// create directory if it does not exist
			final File destChnDir= new File (chn_ill_dir);
			if (!destChnDir.exists()){
				if (!destChnDir.mkdirs()) {
					IJ.showMessage("Error","Failed to create results directory "+chn_ill_dir);
					continue;
				}
			}
			
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						//   					for (int nChn = indxAtomic.getAndIncrement(); nChn < selectedChannels.length; nChn = indxAtomic.getAndIncrement()) if (selectedChannels[nChn]) {
						// iterate through all image set (some grids may be missing)
						for (int nImg = indxAtomic.getAndIncrement();  nImg < src_files[nChn].length; nImg = indxAtomic.getAndIncrement()) {
							//			for (int nImg=0; nImg < src_files[nChn].length; nImg++) {
							ImagePlus imp = new ImagePlus(source_dir+Prefs.getFileSeparator()+src_files[nChn][nImg]);
							String title = imp.getTitle();
							if (title.endsWith(partialSuffix)) {
								title = title.substring(0, title.lastIndexOf("."));
							}
							int [] wh = new int [2];
							double [][] dpixels = getDecimatedStack(
									imp, // ImagePlus imp,
									wh,  // int []    wh, // will return width, height
									src_step, // int       src_step, // margins == src_step/2 
									decimate, // 1 - no decimation
									kernel_widh,
									0.0001); // if none kernel pixel above - return null
							if (dpixels == null) {
								continue; // no kernels
							}
							ImagePlus imp_ill;
							if (sensor_type==1) {
								imp_ill = normalizeColorizeLWIRKernels(
										title+"-ill", // String    title,
										dpixels[0],      // double [] dpixels,
										wh[0],        // int       width,
										illustrationParameters.kernel_direct_lwir_min, // double    mn,
										illustrationParameters.kernel_direct_lwir_max, // double    mx,
										illustrationParameters.kernel_lwir_palette); // int palette);
							} else {
								imp_ill = normalizeColorizeEOKernels(
										title+"-ill",                                   // String    title,
										dpixels,                                        // double [] dpixels,
										wh[0],                                          // int       width,
										illustrationParameters.kernel_direct_red_min,   // double    mnR,
										illustrationParameters.kernel_direct_red_max,   // double    mxR,
										illustrationParameters.kernel_direct_green_min, // double    mnG,
										illustrationParameters.kernel_direct_green_max, // double    mxG,
										illustrationParameters.kernel_direct_blue_min,  // double    mnB,
										illustrationParameters.kernel_direct_blue_max); // double    mnB, 
							}
							// save result
							EyesisCorrections.saveAndShow(
									imp_ill,
									chn_ill_dir,
									illustrationParameters.usePNG(),
									false, // show
									illustrationParameters.JPEG_quality, //  <0 - keep current, 0 - force Tiff, >0 use for JPEG
									0); // debug_level); 

							//	            return true; // temporarily
						}
					}
				};
			}
			startAndJoin(threads);
		}
		///imp_psf.getStack()	
		//int sensor_type = dcd.eyesisCameraParameters.getSensorType(nChn);		
		return true;
	}
	
	public boolean convertDirect(
			final boolean [] selectedChannels,
			final int src_step,
			final String source_dir,
			final String dest_dir,
			final int [][] whall) { // null or per-channel
		String directSuffix= eyesisAberrations.aberrationParameters.psfSuffix; // .".psf-tiff";
		File destDir= new File (dest_dir);
		if (!destDir.exists()){
			if (!destDir.mkdirs()) {
				IJ.showMessage("Error","Failed to create results directory "+dest_dir);
				return false;
			}
		}
		String [][] src_files = new String [selectedChannels.length][]; // should be one file per channel
		for (int nChn = 0; nChn < src_files.length; nChn++) if (selectedChannels[nChn]) {
			MultipleExtensionsFileFilter sourceFilter =
					new MultipleExtensionsFileFilter("",new String[] {String.format("%02d",  nChn)+directSuffix},"direct kernels");
			src_files[nChn] = (new File(source_dir)).list(sourceFilter); // are these full files?
		}
		// create list of partial kernel files
		// partial-1623732618_112155-00.ppsf-tiff
		
		for (int iChn = 0; iChn < src_files.length; iChn++) if (src_files[iChn]!=null) {
			final int nChn=iChn;
			
			int sensor_type = distortions.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getSensorType(nChn);
			int decimate =    (sensor_type==1) ? illustrationParameters.kernel_decimate_direct_lwir : illustrationParameters.kernel_decimate_direct_eo;
			int kernel_widh = (sensor_type==1) ? illustrationParameters.kernel_dia_direct_lwir : illustrationParameters.kernel_dia_direct_eo;

			/*
//			String illustrations_dir = eyesisAberrations.aberrationParameters.illustrationsDirectory;
			final String chn_ill_dir = dest_dir+Prefs.getFileSeparator()+ illustrationParameters.getChannelPrefix()+String.format("%02d", nChn);
			// create directory if it does not exist
			final File destChnDir= new File (chn_ill_dir);
			if (!destChnDir.exists()){
				if (!destChnDir.mkdirs()) {
					IJ.showMessage("Error","Failed to create results directory "+chn_ill_dir);
					continue;
				}
			}
			*/
			//			for (int nImg=0; nImg < src_files[nChn].length; nImg++) {
			ImagePlus imp = new ImagePlus(source_dir+Prefs.getFileSeparator()+src_files[nChn][0]);
			String title = imp.getTitle();
			if (title.endsWith(directSuffix)) {
				title = title.substring(0, title.lastIndexOf("."));
			}
			int [] wh = (whall == null) ? new int [2]: whall[nChn];
			double [][] dpixels = getDecimatedStack(
					imp, // ImagePlus imp,
					wh,  // int []    wh, // will return width, height
					src_step, // int       src_step, // margins == src_step/2 
					decimate, // 1 - no decimation
					kernel_widh,
					0.0); // if none kernel pixel above - return null (here it is not needed) 
			ImagePlus imp_ill;
			if (sensor_type==1) {
				imp_ill = normalizeColorizeLWIRKernels(
						title+"-ill", // String    title,
						dpixels[0],      // double [] dpixels,
						wh[0],        // int       width,
						illustrationParameters.kernel_direct_lwir_min, // double    mn,
						illustrationParameters.kernel_direct_lwir_max, // double    mx,
						illustrationParameters.kernel_lwir_palette); // int palette); 
			} else {
				imp_ill = normalizeColorizeEOKernels(
						title+"-ill",                                   // String    title,
						dpixels,                                        // double [] dpixels,
						wh[0],                                          // int       width,
						illustrationParameters.kernel_direct_red_min,   // double    mnR,
						illustrationParameters.kernel_direct_red_max,   // double    mxR,
						illustrationParameters.kernel_direct_green_min, // double    mnG,
						illustrationParameters.kernel_direct_green_max, // double    mxG,
						illustrationParameters.kernel_direct_blue_min,  // double    mnB,
						illustrationParameters.kernel_direct_blue_max); // double    mnB, 
			}
			// save result
			EyesisCorrections.saveAndShow(
					imp_ill,
					dest_dir, // chn_ill_dir,
					illustrationParameters.usePNG(),
					false, // show
					illustrationParameters.JPEG_quality, //  <0 - keep current, 0 - force Tiff, >0 use for JPEG
					0); // debug_level); 

			//	            return true; // temporarily
		}
		//int sensor_type = dcd.eyesisCameraParameters.getSensorType(nChn);		
		return true;
	}

	public boolean convertInverted(
			final boolean [] selectedChannels,
			final int src_step,
			final String source_dir,
			final String dest_dir,
			final int [][] whall) { // null or per-channel
		String invertedSuffix= eyesisAberrations.aberrationParameters.aberrationsSuffix; // .".psf-tiff";
		File destDir= new File (dest_dir);
		if (!destDir.exists()){
			if (!destDir.mkdirs()) {
				IJ.showMessage("Error","Failed to create results directory "+dest_dir);
				return false;
			}
		}
		String [][] src_files = new String [selectedChannels.length][]; // should be one file per channel
		for (int nChn = 0; nChn < src_files.length; nChn++) if (selectedChannels[nChn]) {
			MultipleExtensionsFileFilter sourceFilter =
					new MultipleExtensionsFileFilter("",new String[] {String.format("%02d",  nChn)+invertedSuffix},"inverted kernels");
			src_files[nChn] = (new File(source_dir)).list(sourceFilter); // are these full files?
		}
		// create list of partial kernel files
		// partial-1623732618_112155-00.ppsf-tiff
		
		for (int iChn = 0; iChn < src_files.length; iChn++) if (src_files[iChn]!=null) {
			final int nChn=iChn;
			
			int sensor_type = distortions.fittingStrategy.distortionCalibrationData.eyesisCameraParameters.getSensorType(nChn);
			int decimate =    (sensor_type==1) ? illustrationParameters.kernel_decimate_inverse_lwir : illustrationParameters.kernel_decimate_inverse_eo;
			int kernel_widh = (sensor_type==1) ? illustrationParameters.kernel_dia_inverse_lwir : illustrationParameters.kernel_dia_inverse_eo;

			ImagePlus imp = new ImagePlus(source_dir+Prefs.getFileSeparator()+src_files[nChn][0]);
			String title = imp.getTitle();
			if (title.endsWith(invertedSuffix)) {
				title = title.substring(0, title.lastIndexOf("."));
			}
//			int [] wh = new int [2];
			int [] wh = (whall == null) ? new int [2]: whall[nChn];
			double [][] dpixels = getDecimatedStack(
					imp, // ImagePlus imp,
					wh,  // int []    wh, // will return width, height
					src_step, // int       src_step, // margins == src_step/2 
					decimate, // 1 - no decimation
					kernel_widh,
					0.0); // if none kernel pixel above - return null (here it is not needed) 
			ImagePlus imp_ill;
			if (sensor_type==1) {
				imp_ill = normalizeColorizeLWIRKernels(
						title+"-ill", // String    title,
						dpixels[0],      // double [] dpixels,
						wh[0],        // int       width,
						illustrationParameters.kernel_inverse_lwir_min, // double    mn,
						illustrationParameters.kernel_inverse_lwir_max, // double    mx,
						illustrationParameters.kernel_lwir_palette); // int palette); 
			} else {
				imp_ill = normalizeColorizeEOKernels(
						title+"-ill",                                   // String    title,
						dpixels,                                        // double [] dpixels,
						wh[0],                                          // int       width,
						illustrationParameters.kernel_inverse_red_min,   // double    mnR,
						illustrationParameters.kernel_inverse_red_max,   // double    mxR,
						illustrationParameters.kernel_inverse_green_min, // double    mnG,
						illustrationParameters.kernel_inverse_green_max, // double    mxG,
						illustrationParameters.kernel_inverse_blue_min,  // double    mnB,
						illustrationParameters.kernel_inverse_blue_max); // double    mnB, 
			}
			// save result
			EyesisCorrections.saveAndShow(
					imp_ill,
					dest_dir, // chn_ill_dir,
					illustrationParameters.usePNG(),
					false, // show
					illustrationParameters.JPEG_quality, //  <0 - keep current, 0 - force Tiff, >0 use for JPEG
					0); // debug_level); 

			//	            return true; // temporarily
		}
		//int sensor_type = dcd.eyesisCameraParameters.getSensorType(nChn);		
		return true;
	}

	
	public double [][] getDecimatedStack(
		ImagePlus imp,
		int []    wh, // will return width, height. If wh[0]!=0, will adjust to it
		int       src_step, // margins == src_step/2 
		int       decimate, // 1 - no decimation
		int       out_kernel_widh,
		double    mingood){ // minimal value to ne considered non-empty kernel
        int width = imp.getWidth();
        int height = imp.getHeight();
		ImageStack stack=imp.getStack();
		float [][] pixels=new float[stack.getSize()][]; // now - 8 (x,y,u,v,contrast, vignR,vignG,vignB
		for (int i=0;i<pixels.length;i++) pixels[i]= (float[]) stack.getPixels(i+1); // pixel X : negative - no grid here
		int khor =   width / src_step; // number of available kernels in a row
		int kvert =  height / src_step; // number of available kernels in a column
		
        int thor =  khor / decimate;
        int tvert = kvert / decimate;
        
//		int margin_left = (src_step + (thor % decimate)) / 2; // with decimation will be increased to make more center-symmetrical
//		int margin_top =  (src_step +(tvert % decimate)) / 2;

//		int margin_left = src_step/2  + src_step*((thor %  decimate) / 2); // with decimation will be increased to make more center-symmetrical
//		int margin_top =  src_step/2  + src_step*((tvert % decimate) / 2);
		
		if ((wh[0] > 0) && (wh[1]>0)) {
			thor = wh[0] / out_kernel_widh;
			tvert = wh[1] / out_kernel_widh;
		}
		int margin_left = src_step/2 + src_step * ((khor -  ((thor-1) *   decimate + 1))/2);
		int margin_top =  src_step/2 + src_step * ((kvert - ((tvert -1) * decimate + 1))/2);
		
        
        wh[0] = thor *   out_kernel_widh;
        wh[1] = tvert *  out_kernel_widh;
        double [][] rslt = new double [pixels.length][wh[0]*wh[1]];
//        for (int i = 0; i < pixels.length; i++) {
//        	Arrays.fill(rslt[i], Double.NaN);
//      }
        boolean nonzero = mingood <= 0; // 0 - do not check
        for (int slice=0; slice < pixels.length; slice++) {
        	for (int tileY=0; tileY < tvert; tileY ++) {
        		int sbase = width* (margin_top + tileY*src_step*decimate - out_kernel_widh / 2);
        		int dbase = wh[0]* (tileY*out_kernel_widh);
        		for (int tileX=0; tileX < thor; tileX ++) {
            		int sindex = sbase + (margin_left + tileX*src_step*decimate - out_kernel_widh / 2);
            		int dindex = dbase + (tileX * out_kernel_widh);
            		
            		for (int i = 0; i < out_kernel_widh; i++) {
            			if (!nonzero) {
            				for (int j = 0; j < out_kernel_widh; j++) {
            					double d = pixels[slice][sindex++];
            					rslt[slice][dindex++] = d;
            					nonzero |= (d > mingood);
            				}
            			} else { // no need to check, faster
            				for (int j = 0; j < out_kernel_widh; j++) {
            					rslt[slice][dindex++] = pixels[slice][sindex++];
            				}

            			}
                		sindex+=(width-out_kernel_widh);
                		dindex+=(wh[0]-out_kernel_widh);
            		}
            		
        		}
        	}
        }
//        (new ShowDoubleFloatArrays()).showArrays(rslt, wh[0],wh[1], true,"kernels_compressed");
		return nonzero ? rslt : null;
	}
	
	public ImagePlus normalizeColorizeLWIRKernels(
			String    title,
			double [] dpixels,
			int       width,
			double    mn,
			double    mx,
			int palette) {
		double [][] pseudo_pixels = new double [4] [dpixels.length];
		ThermalColor tc = new ThermalColor(
				illustrationParameters.getPalette(), // 	public int     lwir_palette =           0; // 0 - white - hot, 1 - black - hot, 2+ - colored
				mn,
				mx,
				255.0);
		for (int i = 0; i < dpixels.length; i++) {
			double [] rgb = tc.getRGB((double) dpixels[i]);
			pseudo_pixels[0][i] = rgb[0]; // red
			pseudo_pixels[1][i] = rgb[1]; // green
			pseudo_pixels[2][i] = rgb[2]; // blue
			pseudo_pixels[3][i] = 1.0; // alpha
		}
		String [] rgb_titles =  {"red","green","blue","alpha"};
		ImageStack stack = (new  ShowDoubleFloatArrays()).makeStack(
				pseudo_pixels, 
				width,         
				dpixels.length/width,
				rgb_titles,
				true);         // replace NaN with 0.0
		ImagePlus imp =  EyesisCorrections.convertRGBAFloatToRGBA32(
				stack,   // ImageStack stackFloat, //r,g,b,a
				//						name+"ARGB"+suffix, // String title,
				title, // String title,
				0.0,   // double r_min,
				255.0, // double r_max,
				0.0,   // double g_min,
				255.0, // double g_max,
				0.0,   // double b_min,
				255.0, // double b_max,
				0.0,   // double alpha_min,
				1.0);  // double alpha_max)
		return imp;
	}

	public ImagePlus normalizeColorizeEOKernels(
			String      title,
			double [][] dpixels, // slices: R,B,G
			int         width,
			double      mnR,
			double      mxR,
			double      mnG,
			double      mxG,
			double      mnB,
			double      mxB) {
		double [][] pseudo_pixels = new double [4] [dpixels[0].length];
		double kR = 255.0/(mxR-mnR);
		double kG = 255.0/(mxG-mnG);
		double kB = 255.0/(mxB-mnB);
		for (int i = 0; i < dpixels[0].length; i++) {
			pseudo_pixels[0][i] = kR * (dpixels[0][i] - mnR); // red
			pseudo_pixels[1][i] = kG * (dpixels[2][i] - mnG); // green
			pseudo_pixels[2][i] = kB * (dpixels[1][i] - mnB); // blue
			pseudo_pixels[3][i] = 1.0; // alpha
		}
		String [] rgb_titles =  {"red","green","blue","alpha"};
		ImageStack stack = (new  ShowDoubleFloatArrays()).makeStack(
				pseudo_pixels, 
				width,         
				dpixels[0].length/width,
				rgb_titles,
				true);         // replace NaN with 0.0
		ImagePlus imp =  EyesisCorrections.convertRGBAFloatToRGBA32(
				stack,   // ImageStack stackFloat, //r,g,b,a
				//						name+"ARGB"+suffix, // String title,
				title, // String title,
				0.0,   // double r_min,
				255.0, // double r_max,
				0.0,   // double g_min,
				255.0, // double g_max,
				0.0,   // double b_min,
				255.0, // double b_max,
				0.0,   // double alpha_min,
				1.0);  // double alpha_max)
		return imp;
	}
	
	public boolean addMissingAsLinks() {
		String footage_dir = eyesisAberrations.aberrationParameters.illustrationsDirectory+Prefs.getFileSeparator()+FOOTAGE_DIR;
		File f_footage_dir=new File(footage_dir);
		if (!f_footage_dir.isDirectory()) {
			String msg="Directory "+footage_dir+" does not exist, aborting";
			IJ.showMessage("Warning",msg);
    		System.out.println("Warning: "+msg);
    		return false;
		}
		int ref_chn = 0; // assume 
// 1627878701_134650_14-footage.png
// 1627878701_134650_16-footage.jpeg
//			String chn_foot_dir = footage_dir+Prefs.getFileSeparator()+ illustrationParameters.getChannelPrefix()+String.format("%02d", nChn);
//file:///home/eyesis/lwir16-proc/captures/floras_lake/footage-illustrations/footage/chn_00/1627878701_134650_0-footage.png		
		String ref_dir = footage_dir+Prefs.getFileSeparator()+ illustrationParameters.getChannelPrefix()+String.format("%02d", ref_chn);
		File f_ref_dir=new File(ref_dir);
		File[] ref_files=f_ref_dir.listFiles(); // all files
//		String[] ref_ts = new String [ref_files.length];
		ArrayList<String> ref_ts_list = new ArrayList<String>();

		for (int i = 0; i < ref_files.length; i++) {
			String p = ref_files[i].getPath();
			int [] start_end = getStartEndTS(p);
			if (start_end != null) {
				ref_ts_list.add(p.substring(start_end[0],start_end[1]));
			}
		}
		Collections.sort(ref_ts_list);
		String [] ref_ts = ref_ts_list.toArray(new String[0]);
		boolean [] selectedChannels =  illustrationParameters.getSelectedChannels();
		for (int nChn = 0; nChn < selectedChannels.length; nChn++) if (selectedChannels[nChn] && (lwirReaderParameters.getTypeMap()[nChn]==LwirReaderParameters.TYPE_EO)) {
			String eo_dir = footage_dir+Prefs.getFileSeparator()+ illustrationParameters.getChannelPrefix()+String.format("%02d", nChn);
			File f_eo_dir=new File(eo_dir);
			File[] eo_files=f_eo_dir.listFiles(); // all files
			String [] eo_paths = new String [eo_files.length];
			for (int i = 0; i < eo_files.length; i++) {
				eo_paths[i] = eo_files[i].getPath();
			}
			Arrays.sort(eo_paths); // assuming single directory
			String last_path = eo_paths[0];
			int chn_indx = 0;
			for (int ref_indx = 0; ref_indx < ref_ts.length; ref_indx++) {
				int [] start_end = null;
				String chn_ts = ""; //eo_paths[chn_indx];
				if (chn_indx < eo_paths.length) {
					start_end = getStartEndTS(eo_paths[chn_indx]);
				}
				if (start_end != null) {
					chn_ts =eo_paths[chn_indx].substring(start_end[0],start_end[1]);
				}
				if (ref_ts[ref_indx].equals(chn_ts)) {
					if (! Files.isSymbolicLink(Paths.get(eo_paths[chn_indx]))) {
						last_path = eo_paths[chn_indx];
					}
					chn_indx++;
				} else {
					start_end = getStartEndTS(last_path);
					String linked_path = last_path.substring(0,start_end[0])+ref_ts[ref_indx]+last_path.substring(start_end[1]);
					Path newLink = Paths.get(linked_path);
					Path target =  Paths.get(last_path);
					try {
					    Files.createSymbolicLink(newLink, target);
					} catch (IOException x) {
					    System.err.println(x);
					} catch (UnsupportedOperationException x) {
					    // Some file systems do not support symbolic links.
					    System.err.println(x);
					}
				}
				
			}
		}
		return true;
	}
	
	public int [] getStartEndTS(String p) {
		int istart = p.lastIndexOf(Prefs.getFileSeparator());
		if (istart < 0) {
			istart = 0;
		} else {
			istart++;
		}
		int iend = p.indexOf('_', istart);
		if (iend > 0) {
			iend = p.indexOf('_', iend+1); // second "_"
			if (iend > 0) {
				return new int[] {istart,iend};
			}
		}
		return null;
	}
	
	public boolean convertCapturedLwirFiles() {
		long startTime=System.nanoTime(); // restart timer after possible interactive dialogs
		final boolean [] selectedChannels =  illustrationParameters.getSelectedChannels();
		boolean has_lwir = false;
		for (int iChn = 0; iChn < selectedChannels.length; iChn++) {
			if (selectedChannels[iChn] && (lwirReaderParameters.getTypeMap()[iChn]==LwirReaderParameters.TYPE_LWIR)) {
				has_lwir = true;
				break;
			}
		}
		if (!has_lwir) {
			System.out.println("No LWIR channels to process");
			return false;
		}


		
		final CapturedScene [] captured_scenes = listCapturedScenes(
				eyesisAberrations.aberrationParameters.capturedDirectory, // String   captured_path,
				illustrationParameters.min_ts,// double   min_ts,
				illustrationParameters.max_ts,// double   max_ts,
				illustrationParameters.captures_all_lwir,
				illustrationParameters.captures_all_eo,
				illustrationParameters.captures_all);
		
		// optionally perform balancing
		double [][] windows = null;
		if (illustrationParameters.calib_offs_gain) {
			CapturedScene [] balancing_scenes = captured_scenes;
			if (    (illustrationParameters.calib_offs_gain_ts < illustrationParameters.min_ts) ||
					((illustrationParameters.calib_offs_gain_ts + illustrationParameters.calib_offs_gain_dur) > illustrationParameters.max_ts)) {
				balancing_scenes = listCapturedScenes(
						eyesisAberrations.aberrationParameters.capturedDirectory, // String   captured_path,
						illustrationParameters.calib_offs_gain_ts,// double   min_ts,
						illustrationParameters.calib_offs_gain_ts + illustrationParameters.calib_offs_gain_dur,// double   max_ts,
						illustrationParameters.captures_all_lwir,
						illustrationParameters.captures_all_eo,
						illustrationParameters.captures_all);


			}
			// read once scene to find width/height
			ImagePlus[] imps = getImagesMultithreaded(
					balancing_scenes[0].images,             // final String []   image_paths,
					(1 << LwirReaderParameters.TYPE_LWIR),  // final int         types_mask, // +1 - TYPE_EO, +2 - TYPE_LWIR
					null); // final double [][] pixels) // if not null - will fill
			windows = getWindows (imps,illustrationParameters.auto_range_wnd_type);
//		(new ShowDoubleFloatArrays()).showArrays(windows, width, height, true,  "windows");
			double [][] offs_gains = balanceOffsGains(
					balancing_scenes, // CapturedScene [] scenes,
					windows,          // double [][]      window,
					illustrationParameters.calib_offs_gain_ts,//  double   min_ts,
					illustrationParameters.calib_offs_gain_ts + illustrationParameters.calib_offs_gain_dur,// double   max_ts,
					illustrationParameters.min_sigma,         //  double  min_sigma, // do not use images with small sigma (noise will be a dominant factor 
					illustrationParameters.noise_sigma);      //  noise_sigma); // estimate of sigma in noise-only images (OK to use 0)
			
			illustrationParameters.setLWIROffsetGains(offs_gains); // will not overwrite unused channels
//			return false; // temporarily
		}
		double [][] offs_gains = illustrationParameters.getLWIROffsetGains();
		// create directories before threads
		for (int nChn = 0; nChn < selectedChannels.length; nChn++) if (selectedChannels[nChn] && (lwirReaderParameters.getTypeMap()[nChn]==LwirReaderParameters.TYPE_LWIR)) {
//		for (int nChn = 0; nChn < selectedChannels.length; nChn++) if (selectedChannels[nChn]) {
			String footage_dir = eyesisAberrations.aberrationParameters.illustrationsDirectory+Prefs.getFileSeparator()+FOOTAGE_DIR;
			String chn_foot_dir = footage_dir+Prefs.getFileSeparator()+ illustrationParameters.getChannelPrefix()+String.format("%02d", nChn);
			// create directory if it does not exist
			File destDir= new File (chn_foot_dir);
			if (!destDir.exists()){
				if (!destDir.mkdirs()) {
					IJ.showMessage("Error","Failed to create results directory "+chn_foot_dir);
					continue;
				}
			}
		}
		final double [][] dpixels = new double [selectedChannels.length][];
		double [][] min_max = new double [selectedChannels.length][];
		
		double    percentile_min = illustrationParameters.percentile_min;            
		double    percentile_max = illustrationParameters.percentile_max;
		double    max_range =      illustrationParameters.max_range;
		double    hot_importance = illustrationParameters.hot_importance;
		boolean   auto_range =     illustrationParameters.auto_range;
		boolean   auto_lim_range = illustrationParameters.auto_lim_range;
		double    autorange_offs_up =     illustrationParameters.autorange_offs_up; //        0.2;  // shift range up (fraction of new value)
		double    autorange_offs_down =   illustrationParameters.autorange_offs_down; //        0.05; // shift range down (fraction of new value)
		double    autorange_range_up =    illustrationParameters.autorange_range_up; //        0.1;  // increase range (fraction of new value)
		double    autorange_range_down =  illustrationParameters.autorange_range_down; //        0.02; // decrease range (fraction of new value)
		final int captures_palette = illustrationParameters.captures_palette;
		
		
		double    prev_lim_high =  0;
		double    prev_lim_range = 0;
		for (int iScene = 0; iScene < captured_scenes.length; iScene++) {
			final int nScene = iScene;
			final ImagePlus[] imps = getImagesMultithreaded(
					captured_scenes[nScene].images,         // final String []   image_paths,
					(1 << LwirReaderParameters.TYPE_LWIR),  // final int         types_mask, // +1 - TYPE_EO, +2 - TYPE_LWIR
					dpixels); // final double [][] pixels)  // if not null - will fill
			
			// TODO: Make sure all images are read in, continue if any failed !!!\
			boolean has_bad_images = false;
			int sensor_types[] = lwirReaderParameters.getTypeMap(); // eo - 0, lwir - 1
			for (int nChn = 0; nChn < captured_scenes[nScene].images.length; nChn++) {
				if ((captured_scenes[nScene].images[nChn] != null) && ((sensor_types[nChn] & (1 << LwirReaderParameters.TYPE_LWIR)) != 0)  && (imps[nChn] == null)) {
					has_bad_images = true;
					System.out.println("***** Failed to read image "+captured_scenes[nScene].images[nChn]+", will skip this scene *****");
				}
			}
			if (has_bad_images) {
				continue;
			}
//			double lim_low = illustrationParameters.lwir_ranges[nChn][0];
//			double lim_high = illustrationParameters.lwir_ranges[nChn][1];
			final double [][] lim_low_high = new double [illustrationParameters.lwir_ranges.length][2];
			for (int i = 0; i < lim_low_high.length; i++) {
				lim_low_high[i][0] = illustrationParameters.lwir_ranges[i][0]; // per channel min/max limits
				lim_low_high[i][1] = illustrationParameters.lwir_ranges[i][1];
			}
			if (auto_range) {
				if (windows == null) {
					windows = getWindows (imps,illustrationParameters.auto_range_wnd_type);
				}
				double [][] histograms = getMultiHistograms( // calculate min, max and normalized histogram (assuming input data are actually integer
						dpixels, // final double [][] dpixels,
						windows, // final double [][] windows,
						min_max); // final double [][] min_max) // histogram bin [0] corresponds to [(int) min, (int) min +1) 
				// Normalize and combine min_max
				double mn=Double.NaN, mx= Double.NaN;
				for (int nChn = 0; nChn < min_max.length; nChn++) if (min_max[nChn] != null) {
					double scaled_mn = min_max[nChn][0]/offs_gains[nChn][1] - offs_gains[nChn][0];
					double scaled_mx = min_max[nChn][1]/offs_gains[nChn][1] - offs_gains[nChn][0];
					if (!(mn <= scaled_mn)) mn = scaled_mn;
					if (!(mx >= scaled_mx)) mx = scaled_mx;
				}
				// now got common min and max

				int ioffs = (int) mn - 1; // extra 1 to avoid checking limits
				int hist_len = (int) Math.ceil(mx) - ioffs +1; //range is +2 to avoid checkin limits 
				double [] histogram = new double[hist_len];
				for (int nChn = 0; nChn < min_max.length; nChn++) if (min_max[nChn] != null) {
					double scale = 1.0/offs_gains[nChn][1];
					//				double offs_add = -offs_gains[nChn][0] - ioffs + 1; // to get common bin >
					double offs_add = (((int) min_max[nChn][0])*scale - offs_gains[nChn][0] - (ioffs+1)); 

					for (int i = 0; i < histograms[nChn].length; i++) {
						int b = (int) (offs_add  + i* scale);
						histogram[b]+=histograms[nChn][i];
					}
				}
				double sum_hist = 0.0;
				for (int i = 0; i < histogram.length; i++) {
					sum_hist += histogram[i];
				}
				// normalize, sum(hist) == 1.0;
				double k = 1.0/sum_hist;
				for (int b = 0; b < histogram.length; b++) {
					histogram[b] *= k; 
				}
				// Find desired limits (or low and range)
				double lim_low = ioffs; // first bin
				double lim_high = ioffs + histogram.length -1; //  last bin
				if ((percentile_min + percentile_max) >= 1.0) {
					percentile_min /= (percentile_min + percentile_max);
					percentile_max = 1.0 - percentile_min; 
				}
				if (percentile_min > 0) {
					double pix_low = percentile_min; // * ipixels.length * 
					double cumul = 0.0;
					int i;
					for (i = 0; i < histogram.length; i++) {
						cumul += histogram[i];
						if (cumul > pix_low) {
							break;
						}
					}
					if (i >=histogram.length) { // should never happen
						lim_low = ioffs + histogram.length -1;
					} else {
						double lo = cumul - histogram[i];
						double kk = (pix_low - lo) / histogram[i];
						lim_low = ioffs + i - 0.5 + kk;
					}
				}
				if (percentile_max > 0) {
					double pix_hi = percentile_max; // *  ipixels.length
					double cumul = 0.0;
					int i;
					for (i = histogram.length - 1; i >= 0; i--) {
						cumul += histogram[i];
						if (cumul > pix_hi) {
							break;
						}
					}
					if (i < 0) { // should never happen
						lim_high = ioffs;
					} else {
						double lo = cumul - histogram[i];
						double kk = (pix_hi - lo) / histogram[i];
						lim_high = ioffs + i + 0.5 - kk;
					}
				}
				if (auto_lim_range && ((lim_high - lim_low)  > max_range)) {
					double [] weighted_cumul = new double [histogram.length];
					double wc = 0.0;
					for (int i = 0; i < histogram.length; i++) {
						double w = i + (1.0 - hot_importance) * lim_low;
						wc+= w * histogram[i];
						weighted_cumul[i] = wc;
					}
					int irange = (int) max_range;
					int ibest = 0;
					double best = weighted_cumul[irange];
					for (int i = 1; i < (histogram.length - irange); i++) {
						double w = weighted_cumul[irange + i] - weighted_cumul[i];
						if (w > best) {
							best = w;
							ibest = i;
						}
					}
					lim_low = ioffs + ibest;
					lim_high = lim_low + irange;
				}
				double lim_range = lim_high - lim_low;
				if (nScene > 0) { // lpf
					if (lim_high > prev_lim_high) {
						lim_high = lim_high * autorange_offs_up +   prev_lim_high * (1.0 - autorange_offs_up);
					} else {
						lim_high = lim_high * autorange_offs_down + prev_lim_high * (1.0 - autorange_offs_down);
					}
					if (lim_range > prev_lim_range) {
						lim_range = lim_range * autorange_range_up +   prev_lim_range * (1.0 - autorange_range_up);
					} else {
						lim_range = lim_range * autorange_range_down + prev_lim_range * (1.0 - autorange_range_down);
					}
					lim_low = lim_high - lim_range;
				}
				prev_lim_high =  lim_high;
				prev_lim_range = lim_range;
				// Now apply shifts/scales to find individual limits
				for (int nChn = 0; nChn < lim_low_high.length; nChn++) {
					lim_low_high[nChn][0] = (lim_low + offs_gains[nChn][0]) * offs_gains[nChn][1]; 
					lim_low_high[nChn][1] = (lim_high + offs_gains[nChn][0]) * offs_gains[nChn][1];
				}
				
				System.out.println(String.format("---- Scene: %s, lim_low=%8.2f, lim_high=%8.2f, lim_range=%8.2f",
						captured_scenes[nScene].name,lim_low,lim_high,lim_range));
				
			} // if (auto_range)
			
	   		final Thread[] threads = newThreadArray(MAX_THREADS);
	   		final AtomicInteger indxAtomic = new AtomicInteger(0);
//	   		final ImagePlus[] imp_array = new ImagePlus[image_paths.length];
	        for (int ithread = 0; ithread < threads.length; ithread++) {
	            threads[ithread] = new Thread() {
	                @Override
	                public void run() {
	                	for (int nChn = indxAtomic.getAndIncrement();  nChn < imps.length; nChn = indxAtomic.getAndIncrement()) if (imps[nChn] != null) {
//	                		imp_array[nChn] = null; // important to set for unused sensor types
	// dpixels
	                		int width = imps[nChn].getWidth();
	                		int height = imps[nChn].getHeight();
	                		String title = imps[nChn].getTitle();
	                		if (title.lastIndexOf(".") > 0) {
	                			title = title.substring(0, title.lastIndexOf("."));
	                		}
	                		if (title.lastIndexOf(Prefs.getFileSeparator()) > 0) {
	                			title = title.substring(title.lastIndexOf(Prefs.getFileSeparator())+1);
	                		}
	                		String title_footage = title+"-footage";
	                		
	                		ImageStack stack = null;
	                		double [][]  pseudo_pixels;
//	                		int line_width = 1; // illustrationParameters.getLineWidthLwir();;
	                		pseudo_pixels = new double [4] [dpixels[nChn].length];
	                		ThermalColor tc = new ThermalColor(
	                				captures_palette, // 	public int     lwir_palette =           0; // 0 - white - hot, 1 - black - hot, 2+ - colored
	                				lim_low_high[nChn][0],
	                				lim_low_high[nChn][1],
	                				255.0);
	                		for (int i = 0; i < dpixels[nChn].length; i++) {
	                			double [] rgb = tc.getRGB(dpixels[nChn][i]);
	                			pseudo_pixels[0][i] = rgb[0]; // red
	                			pseudo_pixels[1][i] = rgb[1]; // green
	                			pseudo_pixels[2][i] = rgb[2]; // blue
	                			pseudo_pixels[3][i] = 1.0; // alpha
	                		}

	                		String [] rgb_titles =  {"red","green","blue","alpha"};
	                		stack = (new  ShowDoubleFloatArrays()).makeStack(
	                				pseudo_pixels, // iclt_data,
	                				width,         // (tilesX + 0) * clt_parameters.transform_size,
	                				height,        // (tilesY + 0) * clt_parameters.transform_size,
	                				rgb_titles,    // or use null to get chn-nn slice names
	                				true);         // replace NaN with 0.0
	                		ImagePlus imp_out =  EyesisCorrections.convertRGBAFloatToRGBA32(
	                				stack,   // ImageStack stackFloat, //r,g,b,a
	                				//						name+"ARGB"+suffix, // String title,
	                				title_footage, // String title,
	                				0.0,   // double r_min,
	                				255.0, // double r_max,
	                				0.0,   // double g_min,
	                				255.0, // double g_max,
	                				0.0,   // double b_min,
	                				255.0, // double b_max,
	                				0.0,   // double alpha_min,
	                				1.0);  // double alpha_max)
	                		
	                		
							if (imp_out == null) {
								continue;
							}
							if (illustrationParameters.captures_annotate) {
								String scene_title = captured_scenes[nScene].name;
								if (scene_title.lastIndexOf(Prefs.getFileSeparator()) > 0) {
									scene_title = scene_title.substring(scene_title.lastIndexOf(Prefs.getFileSeparator())+1);
								}
								ImageProcessor ip = imp_out.getProcessor();
								int posX= width - 119; // 521;
								int posY= height + 1;  // 513;
//								int posX=521;
//								int posY=513;
								Font font =    new Font("Monospaced", Font.PLAIN, 12);
								ip.setColor(illustrationParameters.color_annotate); // Color.BLUE);
								ip.setFont(font);
								ip.drawString(scene_title, posX, posY,Color.BLACK); 
							}
	                		
							String footage_dir = eyesisAberrations.aberrationParameters.illustrationsDirectory+Prefs.getFileSeparator()+FOOTAGE_DIR;
							String chn_foot_dir = footage_dir+Prefs.getFileSeparator()+ illustrationParameters.getChannelPrefix()+String.format("%02d", nChn);
							// create directory if it does not exist
							File destDir= new File (chn_foot_dir);
							if (!destDir.exists()){
								if (!destDir.mkdirs()) {
									IJ.showMessage("Error","Failed to create results directory "+chn_foot_dir);
									continue;
								}
							}
							EyesisCorrections.saveAndShow(
									imp_out,
									chn_foot_dir,
									illustrationParameters.usePNG(),
									false, // show
									illustrationParameters.JPEG_quality, //  <0 - keep current, 0 - force Tiff, >0 use for JPEG
									((nChn==0)?1:0)); // print only for channel 0 
	                		
	                	}
	                }
	            };
	        }
	        startAndJoin(threads);
			
		}	// for (int nScene = 0; nScene < captured_scenes.length; nScene++) {	
	
		/*
		
   		final Thread[] threads = newThreadArray(MAX_THREADS);
   		final AtomicInteger indxAtomic = new AtomicInteger(0);
		
		for (int iChn = 0; iChn < selectedChannels.length; iChn++) if (selectedChannels[iChn] && (lwirReaderParameters.getTypeMap()[iChn]==LwirReaderParameters.TYPE_LWIR)) {
			final int nChn=iChn;
			indxAtomic.set(0);
			// Create directory before threads
			String footage_dir = eyesisAberrations.aberrationParameters.illustrationsDirectory+Prefs.getFileSeparator()+FOOTAGE_DIR;
			String chn_foot_dir = footage_dir+Prefs.getFileSeparator()+ illustrationParameters.getChannelPrefix()+String.format("%02d", nChn);
			// create directory if it does not exist
			File destDir= new File (chn_foot_dir);
			if (!destDir.exists()){
				if (!destDir.mkdirs()) {
					IJ.showMessage("Error","Failed to create results directory "+chn_foot_dir);
					continue;
				}
			}
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						ImagejJp4Tiff imagejJp4Tiff = new ImagejJp4Tiff();
						for (int nScene = indxAtomic.getAndIncrement(); nScene < captured_scenes.length; nScene = indxAtomic.getAndIncrement()) {
							ImagePlus imp_out =  convertCapturedLWIR(
									imagejJp4Tiff,                               // ImagejJp4Tiff imagejJp4Tiff,,
									captured_scenes[nScene].images[nChn],        // String    src_path,
									illustrationParameters.lwir_ranges[nChn][0], // double    abs_min,
									illustrationParameters.lwir_ranges[nChn][1], // double    abs_max,
									illustrationParameters.percentile_min,       // double    percentile_min,            
									illustrationParameters.percentile_max,       // double    percentile_max,
									illustrationParameters.max_range,            // double    max_range,
									illustrationParameters.hot_importance,       // double    hot_importance,
									illustrationParameters.auto_range,           // boolean   auto_range,
									illustrationParameters.auto_lim_range,       // boolean   auto_lim_range, 
									illustrationParameters.captures_annotate,    // boolean   captures_annotate, 
									illustrationParameters.color_annotate,       // Color     color_annotate,
									illustrationParameters.captures_palette);    // int       captures_palette)
							if (imp_out == null) {
								continue;
							}
							if (illustrationParameters.captures_annotate) {
								String scene_title = captured_scenes[nScene].name;
								if (scene_title.lastIndexOf(Prefs.getFileSeparator()) > 0) {
									scene_title = scene_title.substring(scene_title.lastIndexOf(Prefs.getFileSeparator())+1);
								}
								ImageProcessor ip = imp_out.getProcessor();
//								int posX=521;
//								int posY=513;
								int width =  imp_out.getWidth();
								int height = imp_out.getHeight();
								int posX= width - 119; // 521;
								int posY= height + 1;  // 513;
 
								Font font =    new Font("Monospaced", Font.PLAIN, 12);
								ip.setColor(illustrationParameters.color_annotate); // Color.BLUE);
								ip.setFont(font);
								ip.drawString(scene_title, posX, posY,Color.BLACK); 
							}
							String footage_dir = eyesisAberrations.aberrationParameters.illustrationsDirectory+Prefs.getFileSeparator()+FOOTAGE_DIR;
							String chn_foot_dir = footage_dir+Prefs.getFileSeparator()+ illustrationParameters.getChannelPrefix()+String.format("%02d", nChn);
							// create directory if it does not exist
							File destDir= new File (chn_foot_dir);
							if (!destDir.exists()){
								if (!destDir.mkdirs()) {
									IJ.showMessage("Error","Failed to create results directory "+chn_foot_dir);
									continue;
								}
							}
							EyesisCorrections.saveAndShow(
									imp_out,
									chn_foot_dir,
									illustrationParameters.usePNG(),
									false, // show
									illustrationParameters.JPEG_quality, //  <0 - keep current, 0 - force Tiff, >0 use for JPEG
									0); // debug_level); 
						}
					}
				};
			}
			startAndJoin(threads);
		}
		*/
		System.out.println("All done in "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+" sec.");
		return true;
	}
	

	public boolean convertCapturedEoFiles() {
		long startTime=System.nanoTime(); // restart timer after possible interactive dialogs
		final boolean [] selectedChannels =  illustrationParameters.getSelectedChannels();
		boolean has_eo = false;
		for (int iChn = 0; iChn < selectedChannels.length; iChn++) {
			if (selectedChannels[iChn] && (lwirReaderParameters.getTypeMap()[iChn]==LwirReaderParameters.TYPE_EO)) {
				has_eo = true;
				break;
			}
		}
		if (!has_eo) {
			System.out.println("No EO channels to process");
			return false;
		}
		final CapturedScene [] captured_scenes = listCapturedScenes( // will return only scenes that have all 4 EO channels
				eyesisAberrations.aberrationParameters.capturedDirectory, // String   captured_path,
				illustrationParameters.min_ts,// double   min_ts,
				illustrationParameters.max_ts,// double   max_ts,
				illustrationParameters.captures_all_lwir,
				illustrationParameters.captures_all_eo,
				illustrationParameters.captures_all); // true); // illustrationParameters.captures_all_lwir); // illustrationParameters.captures_all);
			
		
   		final Thread[] threads = newThreadArray(MAX_THREADS);
   		final AtomicInteger indxAtomic = new AtomicInteger(0);
		final int eo0 =   illustrationParameters.getLwirReaderParameters().getEoChn0();
		
		for (int iChn = 0; iChn < selectedChannels.length; iChn++) if (selectedChannels[iChn] && (lwirReaderParameters.getTypeMap()[iChn]==LwirReaderParameters.TYPE_EO)) {
			final int nChn=iChn;
			indxAtomic.set(0);
			// Create directory before threads
			String footage_dir = eyesisAberrations.aberrationParameters.illustrationsDirectory+Prefs.getFileSeparator()+FOOTAGE_DIR;
			String chn_foot_dir = footage_dir+Prefs.getFileSeparator()+ illustrationParameters.getChannelPrefix()+String.format("%02d", nChn);
			// create directory if it does not exist
			File destDir= new File (chn_foot_dir);
			if (!destDir.exists()){
				if (!destDir.mkdirs()) {
					IJ.showMessage("Error","Failed to create results directory "+chn_foot_dir);
					continue;
				}
			}
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						ImagejJp4Tiff imagejJp4Tiff = new ImagejJp4Tiff();
						for (int nScene = indxAtomic.getAndIncrement(); nScene < captured_scenes.length; nScene = indxAtomic.getAndIncrement()) {
							ImagePlus imp_out =  convertCapturedEO(
									imagejJp4Tiff,                               // ImagejJp4Tiff imagejJp4Tiff,
									captured_scenes[nScene].images[nChn],        // String      src_path,
									nChn - eo0,                                  // int         eo_chn,
									illustrationParameters.eo_rb2g_hi,           // double [][] eo_rb2g_hi,
									illustrationParameters.getSaturation(),      // double      saturation,
									illustrationParameters.getGamma(),           // double      gamma,
									illustrationParameters.getMinLin());          // double      minlin_gamma, // do not apply gamma to lower values
							if (imp_out == null) {
								continue;
							}
							if (illustrationParameters.captures_annotate) {
								String scene_title = captured_scenes[nScene].name;
								if (scene_title.lastIndexOf(Prefs.getFileSeparator()) > 0) {
									scene_title = scene_title.substring(scene_title.lastIndexOf(Prefs.getFileSeparator())+1);
								}
								ImageProcessor ip = imp_out.getProcessor();
								int width =  imp_out.getWidth();
								int height = imp_out.getHeight();
								int posX= width - 119; // 521;
								int posY= height + 1;  // 513;
								Font font =    new Font("Monospaced", Font.PLAIN, 12);
								ip.setColor(illustrationParameters.color_annotate); // Color.BLUE);
								ip.setFont(font);
								ip.drawString(scene_title, posX, posY,Color.BLACK); 
							}
							String footage_dir = eyesisAberrations.aberrationParameters.illustrationsDirectory+Prefs.getFileSeparator()+FOOTAGE_DIR;
							String chn_foot_dir = footage_dir+Prefs.getFileSeparator()+ illustrationParameters.getChannelPrefix()+String.format("%02d", nChn);
							// create directory if it does not exist
							File destDir= new File (chn_foot_dir);
							if (!destDir.exists()){
								if (!destDir.mkdirs()) {
									IJ.showMessage("Error","Failed to create results directory "+chn_foot_dir);
									continue;
								}
							}
							EyesisCorrections.saveAndShow(
									imp_out,
									chn_foot_dir,
									illustrationParameters.usePNG_EO(),
									false, // show
									illustrationParameters.JPEG_quality, //  <0 - keep current, 0 - force Tiff, >0 use for JPEG
									0); // debug_level); 
						}
					}
				};
			}
			startAndJoin(threads);
		}
		System.out.println("All done in "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+" sec.");
		return true;
	}
	
	
	
	
	ImagePlus convertCapturedLWIR( // not used
			ImagejJp4Tiff imagejJp4Tiff,
			String    src_path,
			double    abs_min,
			double    abs_max,
			double    percentile_min,            
			double    percentile_max,
			double    max_range,
			double    hot_importance,
			boolean   auto_range,
			boolean   auto_lim_range, 
			boolean   captures_annotate, 
		    Color     color_annotate,
			int       captures_palette) {
// read source image		
		ImagePlus imp_src = null;
		try {
			imp_src= imagejJp4Tiff.readTiffJp4(src_path);
		} catch (IOException e) {
			System.out.println("convertCapturedLWIR IOException " + src_path);
		} catch (FormatException e) {
			System.out.println("convertCapturedLWIR FormatException " + src_path);
		}
		if (imp_src == null) {
			return null;
		}
		int width = imp_src.getWidth();
		int height = imp_src.getHeight();
		String title = imp_src.getTitle();
		if (title.lastIndexOf(".") > 0) {
			title = title.substring(0, title.lastIndexOf("."));
		}
		if (title.lastIndexOf(Prefs.getFileSeparator()) > 0) {
			title = title.substring(title.lastIndexOf(Prefs.getFileSeparator())+1);
		}
		String title_footage = title+"-footage";
		
		double lim_low = abs_min;
		double lim_high = abs_max;
		float [] fpixels = (float[]) imp_src.getProcessor().getPixels();
		if (auto_range) {
			// got float processor, but actual values arde int16
			int imin = 0xffff;
			int imax = 0;
			int [] ipixels = new int [fpixels.length];
			for (int i = 0; i < ipixels.length; i++) {
				int d = (int) fpixels[i];
				if (d < imin) imin=d;
				else if (d > imax) imax=d;
				ipixels[i] = d;
			}
			int [] hist = new int [imax -imin +1];
			Arrays.fill(hist,  0);
			for (int i = 0; i < ipixels.length; i++) {
				hist[ipixels[i]-imin]++;
			}
			lim_low = imin;
			lim_high = imax;
			if ((percentile_min + percentile_max) >= 1.0) {
				percentile_min /= (percentile_min + percentile_max);
				percentile_max = 1.0 - percentile_min; 
			}
			if (percentile_min > 0) {
				double pix_low = ipixels.length * percentile_min;
				double cumul = 0.0;
				int i;
				for (i = 0; i < hist.length; i++) {
					cumul += hist[i];
					if (cumul > pix_low) {
						break;
					}
				}
				if (i >=hist.length) { // should never happen
					lim_low = imax;
				} else {
					double lo = cumul - hist[i];
					double k = (pix_low - lo) / hist[i];
					lim_low = imin + i - 0.5 + k;
				}
			}
			if (percentile_max > 0) {
				double pix_hi = ipixels.length * percentile_max;
				double cumul = 0.0;
				int i;
				for (i = hist.length - 1; i >= 0; i--) {
					cumul += hist[i];
					if (cumul > pix_hi) {
						break;
					}
				}
				if (i < 0) { // should never happen
					lim_high = imin;
				} else {
					double lo = cumul - hist[i];
					double k = (pix_hi - lo) / hist[i];
					lim_high = imin + i + 0.5 - k;
				}
			}
			if (auto_lim_range && ((lim_high - lim_low)  > max_range)) {
				double [] weighted_cumul = new double [hist.length];
				double wc = 0.0;
				for (int i = 0; i < hist.length; i++) {
					double w = i + (1.0 - hot_importance) * lim_low;
					wc+= w * hist[i];
					weighted_cumul[i] = wc;
				}
				int irange = (int) max_range;
				int ibest = 0;
				double best = weighted_cumul[irange];
				for (int i = 1; i < (hist.length - irange); i++) {
					double w = weighted_cumul[irange + i] - weighted_cumul[i];
					if (w > best) {
						best = w;
						ibest = i;
					}
				}
				lim_low = imin + ibest;
				lim_high = lim_low + irange;
			}
		}
		
		ImageStack stack = null;
		double [][]  pseudo_pixels;
//		int line_width = 1; // illustrationParameters.getLineWidthLwir();;
		pseudo_pixels = new double [4] [fpixels.length];
		ThermalColor tc = new ThermalColor(
				captures_palette, // 	public int     lwir_palette =           0; // 0 - white - hot, 1 - black - hot, 2+ - colored
				lim_low,
				lim_high,
				255.0);
		for (int i = 0; i < fpixels.length; i++) {
			double [] rgb = tc.getRGB((double) fpixels[i]);
			pseudo_pixels[0][i] = rgb[0]; // red
			pseudo_pixels[1][i] = rgb[1]; // green
			pseudo_pixels[2][i] = rgb[2]; // blue
			pseudo_pixels[3][i] = 1.0; // alpha
		}

		String [] rgb_titles =  {"red","green","blue","alpha"};
		stack = (new  ShowDoubleFloatArrays()).makeStack(
				pseudo_pixels, // iclt_data,
				width,         // (tilesX + 0) * clt_parameters.transform_size,
				height,        // (tilesY + 0) * clt_parameters.transform_size,
				rgb_titles,    // or use null to get chn-nn slice names
				true);         // replace NaN with 0.0
		ImagePlus imp_footage =  EyesisCorrections.convertRGBAFloatToRGBA32(
				stack,   // ImageStack stackFloat, //r,g,b,a
				//						name+"ARGB"+suffix, // String title,
				title_footage, // String title,
				0.0,   // double r_min,
				255.0, // double r_max,
				0.0,   // double g_min,
				255.0, // double g_max,
				0.0,   // double b_min,
				255.0, // double b_max,
				0.0,   // double alpha_min,
				1.0);  // double alpha_max)
		return imp_footage;
	}
	

	ImagePlus convertCapturedEO(
			ImagejJp4Tiff imagejJp4Tiff,
			String      src_path,
			int         eo_chn,
			double [][] eo_rb2g_hi,
			double      saturation,
			double      gamma,
			double      minlin_gamma) { // do not apply gamma to lower values
//			boolean     captures_annotate, 
//		    Color       color_annotate) {
// read source image		
		ImagePlus imp_src = null;
		try {
			imp_src= imagejJp4Tiff.readTiffJp4(src_path);
		} catch (IOException e) {
			System.out.println("convertCapturedEO IOException " + src_path);
		} catch (FormatException e) {
			System.out.println("convertCapturedEO FormatException " + src_path);
		}
		if (imp_src == null) {
			return null;
		}
		int width = imp_src.getWidth();
		int height = imp_src.getHeight();
		String title = imp_src.getTitle();
		if (title.lastIndexOf(".") > 0) {
			title = title.substring(0, title.lastIndexOf("."));
		}
		if (title.lastIndexOf(Prefs.getFileSeparator()) > 0) {
			title = title.substring(title.lastIndexOf(Prefs.getFileSeparator())+1);
		}
		String title_footage = title+"-footage";
		
		ImageStack stack = null;
		double [][]  pseudo_pixels;
//		int line_width = 1; // illustrationParameters.getLineWidthLwir();
		double [][] drgb = MatchSimulatedPattern.simpleDemosaic(
				imp_src,
				eo_rb2g_hi[eo_chn][0],  // r2g,
				eo_rb2g_hi[eo_chn][1],  // b2g,
				saturation,             // saturation,
				gamma,                  // gamma,
				minlin_gamma,           //minlin_gamma, // do not apply gamma to lower values
				eo_rb2g_hi[eo_chn][2]); // ,rgb_hi);        // map to 255, gamma will preserve
		pseudo_pixels = new double [4][];
		for (int i = 0; i < drgb.length; i++) {
			pseudo_pixels[i] = drgb[i];
		}
		pseudo_pixels[3] = new double [pseudo_pixels[0].length];
		Arrays.fill(pseudo_pixels[3], 1.0);
		

		String [] rgb_titles =  {"red","green","blue","alpha"};
		stack = (new  ShowDoubleFloatArrays()).makeStack(
				pseudo_pixels, // iclt_data,
				width,         // (tilesX + 0) * clt_parameters.transform_size,
				height,        // (tilesY + 0) * clt_parameters.transform_size,
				rgb_titles,    // or use null to get chn-nn slice names
				true);         // replace NaN with 0.0
		ImagePlus imp_footage =  EyesisCorrections.convertRGBAFloatToRGBA32(
				stack,   // ImageStack stackFloat, //r,g,b,a
				//						name+"ARGB"+suffix, // String title,
				title_footage, // String title,
				0.0,   // double r_min,
				255.0, // double r_max,
				0.0,   // double g_min,
				255.0, // double g_max,
				0.0,   // double b_min,
				255.0, // double b_max,
				0.0,   // double alpha_min,
				1.0);  // double alpha_max)
		return imp_footage;
	}
	
	
	
	
	
	class CapturedScene{
		String   name;
		String[] images; // indexed by channel
		double   ts;
		CapturedScene (String name, String [] images, double ts){
			this.name = name;
			this.images = images;
			this.ts = ts;
		}
	}
	
	
	CapturedScene [] listCapturedScenes(
			String   captured_path,
			double   min_ts,
			double   max_ts,
			boolean  captures_all_lwir,
			boolean  captures_all_eo,
			boolean  captures_all) {
		int sensor_types[] = lwirReaderParameters.getTypeMap(); // eo - 0, lwir - 1
		int num_sensors[] = {0,0};
		int num_all_sensors = 0; 
		
		for (int st:sensor_types) {
			num_sensors[st]++;
			num_all_sensors++;
		}
		File dFile=new File(captured_path);
		File[] scenesFiles=dFile.listFiles(); // all files
		String[] scenePaths = new String [scenesFiles.length];
		for (int i = 0; i < scenesFiles.length; i++) {
			scenePaths[i] = scenesFiles[i].getPath();
		}
		Arrays.sort(scenePaths);
		// Filter by number of files
		ArrayList<CapturedScene> filteredScenesList = new ArrayList<CapturedScene>();
		
		
//		for (File sceneDir: scenesFiles) {
		for (String scenePath: scenePaths) {
			int basename_start = scenePath.lastIndexOf(Prefs.getFileSeparator());
			if (basename_start >= 0) {
				basename_start++;
			} else {
				basename_start = 0;
			}
			//			double ts = Double.parseDouble(sceneDir.getPath().replace('_','.'));
			double ts = Double.parseDouble(scenePath.substring(basename_start).replace('_','.'));
			if (ts < min_ts) {
				continue;
			}
			if (ts > max_ts) {
				continue;
			}
			int ns[] = {0,0};
			int nsa = 0;
			File sceneDir = new File(scenePath);
			File [] sFiles = sceneDir.listFiles();
			for (File ifile: sFiles) {
				int chn = DistortionCalibrationData.getChannelFromPath(ifile.getPath());
				ns[sensor_types[chn]]++;
				nsa++;
			}
			if (captures_all && (nsa < num_all_sensors)) {
				continue;
			}
			if (captures_all_lwir &&
//					(ns[LwirReaderParameters.TYPE_LWIR] > 0 ) && 
					(ns[LwirReaderParameters.TYPE_LWIR] < num_sensors[LwirReaderParameters.TYPE_LWIR])) {
				continue;
			}
			if (captures_all_eo &&
					((ns[LwirReaderParameters.TYPE_EO] > 0 ) || captures_all) && // captures_all - zero is not an option
					(ns[LwirReaderParameters.TYPE_EO] < num_sensors[LwirReaderParameters.TYPE_EO])) {
				continue;
			}
			String [] images = new String [num_all_sensors];
			for (File ifile: sFiles) {
				int chn = DistortionCalibrationData.getChannelFromPath(ifile.getPath());
				images[chn] = ifile.getPath();
				ns[sensor_types[chn]]++;
			}
			filteredScenesList.add(new CapturedScene(sceneDir.getPath(), images, ts));
			
		}
		System.out.println("Selected "+(filteredScenesList.size())+" scenes");
		return filteredScenesList.toArray(new CapturedScene[0]);
	}

	double [][] getNewMeanSigma(
		final double [][] pixels,
		final double [][] wnd,
		final double      noise_sigma) // estimation of sigma with noise only
	{ 
		final double [][] new_mean_sigma = new double [pixels.length][];
		
   		final Thread[] threads = newThreadArray(MAX_THREADS);
   		final AtomicInteger indxAtomic = new AtomicInteger(0);
		
			
        for (int ithread = 0; ithread < threads.length; ithread++) {
            threads[ithread] = new Thread() {
                @Override
                public void run() {
                    for (int nChn = indxAtomic.getAndIncrement();  nChn < new_mean_sigma.length; nChn = indxAtomic.getAndIncrement()) {
                    	double [] pix = pixels[nChn];
                    	if (pix != null) {
                    		double s0=0.0, sx=0.0, sx2=0.0;
                    		for (int i = 0; i < pix.length; i++) {
                    			double w = (wnd == null)? 1.0:wnd[nChn][i];
                    			s0+= w;
                    			double d = pix[i];
                    			double wx = w * d;
                    			sx += wx;
                    			sx2+= wx * pix[i];
                    		}
                    		new_mean_sigma[nChn] = new double[2];
                    		new_mean_sigma[nChn][0] = sx/s0;                          // mean
                    		double sigma2 = (s0*sx2 - sx*sx) /s0/s0; //  - (noise_sigma * noise_sigma);
                    		if (noise_sigma > 0.0) {
                    			sigma2 -= noise_sigma * noise_sigma;
                    			sigma2 = Math.max(sigma2, 0.0);
                    		}
                    		new_mean_sigma[nChn][1] = Math.sqrt(sigma2); // sigma
                    	}
                    }
                }
            };
        }
        startAndJoin(threads);
		return new_mean_sigma;
	}
	
	ImagePlus[] getImagesMultithreaded(
			final String []   image_paths,
			final int         types_mask, // +1 - TYPE_EO, +2 - TYPE_LWIR
			final double [][] pixels) // if not null - will fill
		{ 
			int sensor_types[] = lwirReaderParameters.getTypeMap(); // eo - 0, lwir - 1
			
	   		final Thread[] threads = newThreadArray(MAX_THREADS);
	   		final AtomicInteger indxAtomic = new AtomicInteger(0);
	   		final ImagePlus[] imp_array = new ImagePlus[image_paths.length];
	        for (int ithread = 0; ithread < threads.length; ithread++) {
	            threads[ithread] = new Thread() {
	                @Override
	                public void run() {
	                	ImagejJp4Tiff imagejJp4Tiff = new ImagejJp4Tiff();
	                	for (int nChn = indxAtomic.getAndIncrement();  nChn < image_paths.length; nChn = indxAtomic.getAndIncrement()) {
	                		imp_array[nChn] = null; // important to set for unused sensor types
	                		if (((1 << sensor_types[nChn]) & types_mask) != 0) {
	                			try {
	                				imp_array[nChn]= imagejJp4Tiff.readTiffJp4(image_paths[nChn]);
	                			} catch (IOException e) {
	                				System.out.println("getImagesMultithreaded IOException " + image_paths[nChn]);
	                			} catch (FormatException e) {
	                				System.out.println("getImagesMultithreaded FormatException " + image_paths[nChn]);
	                			}
	                			if (imp_array[nChn] != null) {
	                				if (pixels != null) {
	                					float [] fpixels = (float[]) imp_array[nChn].getProcessor().getPixels();
	                					pixels[nChn] = new double[fpixels.length];
	                					for (int i = 0; i < fpixels.length; i++) {
	                						pixels[nChn][i] = fpixels[i]; 
	                					}
	                				}
	                			}
	                		}
	                	}
	                }
	            };
	        }
	        startAndJoin(threads);
			return imp_array;
		}

	double [][] getMultiHistograms( // calculate min, max and normalized histogram (assuming input data are actually integer
			final double [][] dpixels,
			final double [][] windows,
			final double [][] min_max) // histogram bin [0] corresponds to [(int) min, (int) min +1) 
		{ 
			final double [][] histograms = new double[dpixels.length][];
	   		final Thread[] threads = newThreadArray(MAX_THREADS);
	   		final AtomicInteger indxAtomic = new AtomicInteger(0);
	        for (int ithread = 0; ithread < threads.length; ithread++) {
	            threads[ithread] = new Thread() {
	                @Override
	                public void run() {
	                	ImagejJp4Tiff imagejJp4Tiff = new ImagejJp4Tiff();
	                	for (int nChn = indxAtomic.getAndIncrement();  nChn < dpixels.length; nChn = indxAtomic.getAndIncrement()) {
	                		double [] pixels = dpixels[nChn];
	                		if (pixels == null) {
		                		histograms[nChn] = null;
		                		min_max[nChn] =    null;
	                		} else  {
	                			double [] window = windows[nChn];
		                		double mn = pixels[0], mx = pixels[0];
		                		for (int i = 1; i < pixels.length; i++) {
		                			if      (mn > pixels[i]) mn = pixels[i];
		                			else if (mx < pixels[i]) mx = pixels[i];
		                		}
		                		int offs = (int) mn;
		                		int hist_len = (int) Math.ceil(mx) - offs + 1;
		                		double sum_hist = 0.0;
		                		double [] histogram = new double[hist_len];
		                		for (int i = 0; i < pixels.length; i++) {
		                			int b = (int) pixels[i] - offs;
		                			double w = window[i];
		                			histogram[b] += w;
		                			sum_hist += w;
		                		}
		                		// noprmailize, sum(hist) == 1.0;
		                		double k = 1.0/sum_hist;
		                		for (int b = 0; b < histogram.length; b++) {
		                			histogram[b] *= k; 
		                		}
		                		histograms[nChn] = histogram;
		                		min_max[nChn] = new double[] {mn,mx};
	                		}
	                	}
	                }
	            };
	        }
	        startAndJoin(threads);
			return histograms;
		}
	
	
	
	double [][] balanceOffsGains(
			CapturedScene [] scenes,
			double [][]      window,
			double           min_ts, // if not overlapping, will need separate listCapturedScenes() run
			double           max_ts,
			double           min_sigma, // do not use images with small sigma (noise will be a dominant factor
		                                 // so use large thresholds/good contrast images
			double      noise_sigma) { // estimate of sigma in noise-only images (OK to use 0)
//		double [][]                    offs_scale;// [chn][0] - offset to subtract to normalize, [chn][1] - scale to divide by to normalize
		int sensor_types[] = lwirReaderParameters.getTypeMap(); // eo - 0, lwir - 1
		double [][] offs_scale = new double [sensor_types.length][2];
		int num_imgs[] = new int[sensor_types.length];
		double [][] dpixels = new double [sensor_types.length][];
		for (CapturedScene scene: scenes) {
			if ((scene.ts >= min_ts) && (scene.ts <= max_ts)) {
				ImagePlus [] imps = getImagesMultithreaded(
						scene.images, // final String []   image_paths,
						(1 << LwirReaderParameters.TYPE_LWIR), // final int         types_mask, // +1 - TYPE_EO, +2 - TYPE_LWIR
						dpixels); // final double [][] pixels) // if not null - will fill
				
				double [][] mean_sigmas = getNewMeanSigma(
						dpixels, // final double [][] pixels,
						window, // final double [] wnd);
						noise_sigma); // final double      noise_sigma);
				for (int i = 0; i < mean_sigmas.length; i++) {
					if (mean_sigmas[i] != null) {
						if (mean_sigmas[i][1] >= min_sigma) {
							num_imgs[i]++;
							offs_scale[i][0] += mean_sigmas[i][0];
							offs_scale[i][1] += mean_sigmas[i][1];
						}
					}
				}
			}
		}
		double mean_sigma_log = 0.0;
		int num_avg = 0;
		for (int i = 0; i < num_imgs.length; i++) {
			if (num_imgs[i] > 0) {
				offs_scale[i][0] /= num_imgs[i];
				offs_scale[i][1] /= num_imgs[i];
				num_avg++;
				mean_sigma_log += Math.log(offs_scale[i][1]); 
			}
		}		
		mean_sigma_log /= num_avg;
		double mean_offs = 0;
		double mean_sigma = Math.exp(mean_sigma_log); // replace with logs!
		for (int i = 0; i < num_imgs.length; i++) {
			if (num_imgs[i] > 0) {
				offs_scale[i][1] /= mean_sigma;
				offs_scale[i][0] /= offs_scale[i][1]; // normalize gain
				mean_offs += offs_scale[i][0]; 
			}
		}
		mean_offs /= num_avg;
		for (int i = 0; i < num_imgs.length; i++) {
			if (num_imgs[i] > 0) {
				 offs_scale[i][0] -= mean_offs;
			} else {
				offs_scale[i] = null;
			}
		}
		return offs_scale; // first scale, then shift
	}

/*
		double mean_offs = 0;
		double mean_sigma_log = 0.0;
		int num_avg = 0;
		for (int i = 0; i < num_imgs.length; i++) {
			if (num_imgs[i] > 0) {
				offs_scale[i][0] /= num_imgs[i];
				offs_scale[i][1] /= num_imgs[i];
				num_avg++;
				mean_offs += offs_scale[i][0]; 
				mean_sigma_log += Math.log(offs_scale[i][1]); 
			}
		}		
		mean_offs /= num_avg;
		mean_sigma_log /= num_avg;
		double mean_sigma = Math.exp(mean_sigma_log); // replace with logs!
		for (int i = 0; i < num_imgs.length; i++) {
			if (num_imgs[i] > 0) {
				 offs_scale[i][0] -= mean_offs;
				 offs_scale[i][1] /= mean_sigma;
			} else {
				offs_scale[i] = null;
			}
		}
	
 */
	double [][] getWindows (ImagePlus [] imps, int wnd_type) {
		double [][] windows = new double [imps.length][];
		int i0 = -1;
		int width = 0;
		int height = 0;
		for (int i = 0; i < imps.length; i++) {
			if (imps[i] != null) {
				width = imps[i].getWidth();
				height = imps[i].getWidth();
				if ((i0 < 0) || (width*height != windows[i0].length)) { // Does not handle different image sizes - possible if needed
					i0 = i;
					windows[i] = getWindow (
							width,
							height,
							wnd_type); //  3; // 0 - piramid, 1 half-sin, 2-piramid squared, 3 - sin^2
				} else {
					windows[i] = windows[i0]; // by reference
				}
			}
		}
		return windows;
	}

	
	double [] getWindow (int width, int height, int type) {
		double [] wnd = new double [width * height];
		double sy, sx;
		double kyl = 2.0/height;
		double kxl = 2.0/width;
		boolean lin = (type == 0) || (type == 2);
		for (int y = 0; y < height; y++) {
			if (lin) {
				sy = kyl * ((y >= height / 2) ? (height - y) : (y+1));
			} else {
				sy = Math.sin(Math.PI * y / height);

			}
			for (int x = 0; x < width; x++) {
				if (lin) {
					sx = kxl * ((x >= width / 2) ? (width - y) : (y+1));
				} else {
					sx = Math.sin(Math.PI * x / width);
				}
				wnd[width*y + x] = sy * sx;
			}
		}
		if ((type == 1) || (type == 3)){ // squared
			for (int i = 0; i < wnd.length; i++) {
				wnd[i] *= wnd[i];
			}
		}

		return wnd;
	}
	
	public boolean convertSourceFiles() {
		long startTime=System.nanoTime(); // restart timer after possible interactive dialogs
//		final boolean [] selectedChannels = eyesisAberrations.aberrationParameters.getChannelSelection(distortions);
		final boolean [] selectedChannels =  illustrationParameters.getSelectedChannels();
		final DistortionCalibrationData dcd = distortions.fittingStrategy.distortionCalibrationData;
		final MultipleExtensionsFileFilter sourceFilter =
				new MultipleExtensionsFileFilter("",src_extensions,"Source calibration images");
		final int lwir0 = illustrationParameters.getLwirReaderParameters().getLwirChn0();
		final int eo0 =   illustrationParameters.getLwirReaderParameters().getEoChn0();
   		final Thread[] threads = newThreadArray(MAX_THREADS);
   		final AtomicInteger indxAtomic = new AtomicInteger(0);
		for (int iChn = 0; iChn < selectedChannels.length; iChn++) if (selectedChannels[iChn]) {
			final int nChn=iChn;
			indxAtomic.set(0);
			// Create directory before threads
			String illustrations_dir = eyesisAberrations.aberrationParameters.illustrationsDirectory;
			String chn_ill_dir = illustrations_dir+Prefs.getFileSeparator()+ illustrationParameters.getChannelPrefix()+String.format("%02d", nChn);
			// create directory if it does not exist
			File destDir= new File (chn_ill_dir);
			if (!destDir.exists()){
				if (!destDir.mkdirs()) {
					IJ.showMessage("Error","Failed to create results directory "+chn_ill_dir);
					continue;
				}
			}
			
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						//   					for (int nChn = indxAtomic.getAndIncrement(); nChn < selectedChannels.length; nChn = indxAtomic.getAndIncrement()) if (selectedChannels[nChn]) {
						// iterate through all image set (some grids may be missing)
						for (int nSet = indxAtomic.getAndIncrement(); nSet < dcd.gIS.length; nSet = indxAtomic.getAndIncrement()) {
							if (nSet == 679) {
								System.out.println("Debug Set="+nSet);
							}
							int station = dcd.gIS[nSet].getStationNumber();
							if (illustrationParameters.useStation(station)) { // some stations only
								// construct source file name
								String srcPath = null;
								int numImg = -1;
								if (dcd.gIS[nSet].imageSet[nChn] != null) {
									srcPath = dcd.gIS[nSet].imageSet[nChn].source_path;
									numImg = dcd.gIS[nSet].imageSet[nChn].getImageNumber();
								} else {
									// find other non-null
									for (int i = 0; i < dcd.gIS[nSet].imageSet.length; i++) {
										if (dcd.gIS[nSet].imageSet[i] != null) {
											String other_path = dcd.gIS[nSet].imageSet[i].source_path;
											String set_path = other_path.substring(0,other_path.lastIndexOf(Prefs.getFileSeparator()));
											File set_dir = new File(set_path);
											String [] sfiles = set_dir.list(sourceFilter);
											for (String spath:sfiles) {
												int chn = DistortionCalibrationData.pathToChannel(spath);
												if (chn == nChn) {
													srcPath = (new File(set_dir,spath)).getPath();
													break;
												}
											}
											break;
										}
									}
								}
								// open 32-bit image
								if (srcPath == null) {
									System.out.println ("Source image for set "+nSet+", channel "+nChn+" does not exist");
									continue;
								}
								ImagePlus imp = new ImagePlus(srcPath);
								// convert to 8-bit color?
								int width = imp.getWidth();
								int height = imp.getHeight();
								float [] pixels = (float[]) imp.getProcessor().getPixels();
								String title = imp.getTitle();
								if (title.lastIndexOf(".") > 0) {
									title = title.substring(0, title.lastIndexOf("."));
								}
								String title_annot = title+"-annot";
								if (illustrationParameters.useStationInFilenames()) {
									title_annot = "station-"+station+"_"+title_annot;
								}
								if (numImg >=0) { // check if it is a bad image
									int num_above = 0; 
									GridImageParameters gip = dcd.gIP[numImg];
									double [][] pXY =       gip.pixelsXY;
									for (int i = 0; i < pXY.length; i++ ) {
										if (pXY[i][CONTRAST_INDEX] > illustrationParameters.getThresholdContrast()) num_above++;
									}
									if (num_above < illustrationParameters.getThresholdNumber()) {
										title_annot += "-BAD";
									}
								} else {
									title_annot += "-EMPTY";
								}
								int sensor_type = dcd.eyesisCameraParameters.getSensorType(nChn);
								ImageStack stack = null;
								double [][]  pseudo_pixels;
								int line_width = 1;
								if (sensor_type == 1) {
									pseudo_pixels = new double [4] [pixels.length];
									ThermalColor tc = new ThermalColor(
											illustrationParameters.getPalette(), // 	public int     lwir_palette =           0; // 0 - white - hot, 1 - black - hot, 2+ - colored
											illustrationParameters.getLwirRange(nChn- lwir0)[0],
											illustrationParameters.getLwirRange(nChn- lwir0)[1],
											255.0);
									for (int i = 0; i < pixels.length; i++) {
										double [] rgb = tc.getRGB((double) pixels[i]);
										pseudo_pixels[0][i] = rgb[0]; // red
										pseudo_pixels[1][i] = rgb[1]; // green
										pseudo_pixels[2][i] = rgb[2]; // blue
										pseudo_pixels[3][i] = 1.0; // alpha
									}
									line_width = illustrationParameters.getLineWidthLwir();
								} else { // eo
									double [][] drgb = MatchSimulatedPattern.simpleDemosaic(
											imp,
											illustrationParameters.eo_rb2g_hi[nChn-eo0][0], // r2g,
											illustrationParameters.eo_rb2g_hi[nChn-eo0][1], // b2g,
											illustrationParameters.getSaturation(),         // saturation,
											illustrationParameters.getGamma(),              // gamma,
											illustrationParameters.getMinLin(),             //minlin_gamma, // do not apply gamma to lower values
											illustrationParameters.eo_rb2g_hi[nChn-eo0][2]); // ,rgb_hi);        // map to 255, gamma will preserve
									pseudo_pixels = new double [4][];
									for (int i = 0; i < drgb.length; i++) {
										pseudo_pixels[i] = drgb[i];
									}
									pseudo_pixels[3] = new double [pseudo_pixels[0].length];
									Arrays.fill(pseudo_pixels[3], 1.0);
									line_width = illustrationParameters.getLineWidthEo();
								}

								String [] rgb_titles =  {"red","green","blue","alpha"};
								stack = (new  ShowDoubleFloatArrays()).makeStack(
										pseudo_pixels, // iclt_data,
										width,         // (tilesX + 0) * clt_parameters.transform_size,
										height,        // (tilesY + 0) * clt_parameters.transform_size,
										rgb_titles,    // or use null to get chn-nn slice names
										true);         // replace NaN with 0.0
								ImagePlus imp_annot =  EyesisCorrections.convertRGBAFloatToRGBA32(
										stack,   // ImageStack stackFloat, //r,g,b,a
										//						name+"ARGB"+suffix, // String title,
										title_annot, // String title,
										0.0,   // double r_min,
										255.0, // double r_max,
										0.0,   // double g_min,
										255.0, // double g_max,
										0.0,   // double b_min,
										255.0, // double b_max,
										0.0,   // double alpha_min,
										1.0);  // double alpha_max)
								if (numImg >=0) {
										if (numImg == 2276) {
											System.out.println(">>>numImg="+numImg);
										}
									plotGrid(numImg,
											line_width,
											imp_annot,
											illustrationParameters.getGridColor(), // new Color(250, 0,  0), //  color_grid,
											illustrationParameters.getGridWeakColor(), // new Color(250, 0,  0), //  color_grid,
											illustrationParameters.getGridExtraColor(), // new Color(200, 200,0) // null // //at least one end points to extra (unreliable) nodes
											illustrationParameters.getThresholdContrast()
											);
								}

								String illustrations_dir = eyesisAberrations.aberrationParameters.illustrationsDirectory;
								String chn_ill_dir = illustrations_dir+Prefs.getFileSeparator()+ illustrationParameters.getChannelPrefix()+String.format("%02d", nChn);
								// create directory if it does not exist
								File destDir= new File (chn_ill_dir);
								if (!destDir.exists()){ // Should be created before threads !
									if (!destDir.mkdirs()) {
										IJ.showMessage("Error","Failed to create results directory "+chn_ill_dir+". *** It should be already created! ***");
										continue;
									}
								}
								EyesisCorrections.saveAndShow(
										imp_annot,
										chn_ill_dir,
										illustrationParameters.usePNG(),
										false, // show
										illustrationParameters.JPEG_quality, //  <0 - keep current, 0 - force Tiff, >0 use for JPEG
										0); // debug_level); 
							}
						}
					}
				};
			}
			startAndJoin(threads);
		}
		System.out.println("All done in "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+" sec.");
		return true;
	}

	
	public boolean removeBadGrids() {
		long startTime=System.nanoTime(); // restart timer after possible interactive dialogs
//		final boolean [] selectedChannels = eyesisAberrations.aberrationParameters.getChannelSelection(distortions);
		final boolean [] selectedChannels =  illustrationParameters.getSelectedChannels();
		final DistortionCalibrationData dcd = distortions.fittingStrategy.distortionCalibrationData;
   		final Thread[] threads = newThreadArray(MAX_THREADS);
   		final AtomicInteger indxAtomic = new AtomicInteger(0);
   		final AtomicInteger numRemoved = new AtomicInteger(0);
		for (int iChn = 0; iChn < selectedChannels.length; iChn++) if (selectedChannels[iChn]) {
			final int nChn=iChn;
			indxAtomic.set(0);
			// Create directory before threads
			String illustrations_dir = eyesisAberrations.aberrationParameters.illustrationsDirectory;
			String chn_ill_dir = illustrations_dir+Prefs.getFileSeparator()+ illustrationParameters.getChannelPrefix()+String.format("%02d", nChn);
			// create directory if it does not exist
			File destDir= new File (chn_ill_dir);
			if (!destDir.exists()){
				if (!destDir.mkdirs()) {
					IJ.showMessage("Error","Failed to create results directory "+chn_ill_dir);
					continue;
				}
			}
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						//   					for (int nChn = indxAtomic.getAndIncrement(); nChn < selectedChannels.length; nChn = indxAtomic.getAndIncrement()) if (selectedChannels[nChn]) {
						// iterate through all image set (some grids may be missing)
						for (int nSet = indxAtomic.getAndIncrement(); nSet < dcd.gIS.length; nSet = indxAtomic.getAndIncrement()) {
							int station = dcd.gIS[nSet].getStationNumber();
							if (illustrationParameters.useStation(station)) { // some stations only
								// construct source file name
								if (dcd.gIS[nSet].imageSet[nChn] != null) {
									int numImg = dcd.gIS[nSet].imageSet[nChn].getImageNumber();
									int num_above = 0; 
									GridImageParameters gip = dcd.gIP[numImg];
									double [][] pXY =       gip.pixelsXY;
									for (int i = 0; i < pXY.length; i++ ) {
										if (pXY[i][CONTRAST_INDEX] > illustrationParameters.getThresholdContrast()) num_above++;
									}
									if (num_above < illustrationParameters.getThresholdNumber()) {
										String grid_path = dcd.getImagePath(numImg);
										System.out.println("Removing bad grid file: "+grid_path);
										new File(dcd.getImagePath(numImg)).delete();
										numRemoved.getAndIncrement();
									}
								}
							}
						}
					}
				};
			}
			startAndJoin(threads);
		}
		System.out.println("Removed "+(numRemoved.get())+" bad grid files in "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+" sec.");
		return true;
	}

	
	
	
	/* Create a Thread[] array as large as the number of processors available.
	 * From Stephan Preibisch's Multithreading.java class. See:
	 * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD
	 */
	private static Thread[] newThreadArray(int maxCPUs) {
		int n_cpus = Runtime.getRuntime().availableProcessors();
		if (n_cpus>maxCPUs)n_cpus=maxCPUs;
		return new Thread[n_cpus];
	}
/* Start all given threads and wait on each of them until all are done.
	 * From Stephan Preibisch's Multithreading.java class. See:
	 * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD
	 */
	private static void startAndJoin(Thread[] threads)
	{
		for (int ithread = 0; ithread < threads.length; ++ithread)
		{
			threads[ithread].setPriority(Thread.NORM_PRIORITY);
			threads[ithread].start();
		}

		try
		{
			for (int ithread = 0; ithread < threads.length; ++ithread)
				threads[ithread].join();
		} catch (InterruptedException ie)
		{
			throw new RuntimeException(ie);
		}
	}
	
}
