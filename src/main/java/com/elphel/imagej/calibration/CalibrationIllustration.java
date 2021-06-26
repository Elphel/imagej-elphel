package com.elphel.imagej.calibration;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Rectangle;
import java.io.File;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.calibration.DistortionCalibrationData.GridImageParameters;
import com.elphel.imagej.cameras.ThermalColor;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.correction.EyesisCorrections;
import com.elphel.imagej.lwir.LwirReaderParameters;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.Line;
import ij.process.ImageProcessor;

public class CalibrationIllustration {
	public static final int MAX_THREADS = 100; // combine from all classes?
	public static final int CONTRAST_INDEX = 2;

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
	
	public CalibrationIllustration (
			CalibrationIllustrationParameters illustrationParameters,			
			EyesisAberrations              eyesisAberrations,
			Distortions                    distortions,
			AtomicInteger                  stopRequested,
			int                            debug_level) {
		this.illustrationParameters =      illustrationParameters;
		this.eyesisAberrations =           eyesisAberrations;
		this.distortions =                 distortions;
		this.stopRequested =               stopRequested;
		this.debug_level =                 debug_level;
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
								ip.setColor(color_grid_weak);
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
	
	public boolean convertSourceFiles() {
		long startTime=System.nanoTime(); // restart timer after possible interactive dialogs
		final boolean [] selectedChannels = eyesisAberrations.aberrationParameters.getChannelSelection(distortions);
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
										illustrationParameters.save_png,
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
		final boolean [] selectedChannels = eyesisAberrations.aberrationParameters.getChannelSelection(distortions);
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
		System.out.println("Remoded "+(numRemoved.get())+" bad grid files in "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+" sec.");
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
