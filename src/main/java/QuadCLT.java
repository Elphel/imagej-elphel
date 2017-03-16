/**
 **
 ** QuadCLT - Process images with CLT-based methods (code specific to ImageJ plugin)
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  QuadCLT.java is free software: you can redistribute it and/or modify
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
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Properties;
import java.util.concurrent.atomic.AtomicInteger;

import ij.CompositeImage;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.io.FileSaver;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;


public class QuadCLT {
	static String []                                       fine_corr_coeff_names = {"A","B","C","D","E","F"};
	static String []                                       fine_corr_dir_names = {"X","Y"};
	static String                                          prefix = "EYESIS_DCT."; // change later (first on save)
	public Properties                                      properties = null;
	public EyesisCorrections                               eyesisCorrections = null;
	public EyesisCorrectionParameters.CorrectionParameters correctionsParameters=null;
	double [][][][][][]                                    clt_kernels = null;
	GeometryCorrection                                     geometryCorrection = null;
	public int                                             extra_items = 8; // number of extra items saved with kernels (center offset (partial, full, derivatives)
	public ImagePlus                                       eyesisKernelImage = null;
	public long                                            startTime;
	public double [][][]                                   fine_corr  = new double [4][2][6]; // per port, per x/y, set of 6 coefficient for fine geometric corrections 
	
	TileProcessor                                          tp = null;
	
// magic scale should be set before using  TileProcessor (calculated disparities depend on it)

	public void setTiles (ImagePlus imp, // set tp.tilesX, tp.tilesY
			EyesisCorrectionParameters.CLTParameters    clt_parameters,
			int threadsMax
			){
		if (tp == null){
			tp = new TileProcessor(imp.getWidth()/clt_parameters.transform_size,
					imp.getHeight()/clt_parameters.transform_size,
					clt_parameters.transform_size,
					clt_parameters.stSize,
					clt_parameters.corr_magic_scale,
					threadsMax);
		}  
	}
	
	public QuadCLT(
			Properties                                      properties,
			EyesisCorrections                               eyesisCorrections,
			EyesisCorrectionParameters.CorrectionParameters correctionsParameters
			){
		this.eyesisCorrections=      eyesisCorrections;
		this.correctionsParameters = correctionsParameters;
		this.properties =            properties;
		getProperties();
	}

	// TODO:Add saving just calibration
	public void setProperties(){
		for (int n = 0; n < fine_corr.length; n++){
			for (int d = 0; d < fine_corr[n].length; d++){
				for (int i = 0; i < fine_corr[n][d].length; i++){
					String name = prefix+"fine_corr_"+n+fine_corr_dir_names[d]+fine_corr_coeff_names[i];
					properties.setProperty(name,   this.fine_corr[n][d][i]+"");
				}
			}
		}
	}	

	public void getProperties(){
		for (int n = 0; n < fine_corr.length; n++){
			for (int d = 0; d < fine_corr[n].length; d++){
				for (int i = 0; i < fine_corr[n][d].length; i++){
					String name = prefix+"fine_corr_"+n+fine_corr_dir_names[d]+fine_corr_coeff_names[i];
		  			if (properties.getProperty(name)!=null) this.fine_corr[n][d][i]=Double.parseDouble(properties.getProperty(name));
				}
			}
		}
	}	


	public void setKernelImageFile(ImagePlus img_kernels){
		eyesisKernelImage = img_kernels;
	}
	
	public boolean kernelImageSet(){
		return eyesisKernelImage != null;
	}

	public boolean CLTKernelsAvailable(){
		return clt_kernels != null;
	}
	public boolean geometryCorrectionAvailable(){
		return geometryCorrection != null;
	}
	public boolean initGeometryCorrection(int debugLevel){
		geometryCorrection = new GeometryCorrection();
		PixelMapping.SensorData [] sensors =  eyesisCorrections.pixelMapping.sensors;
		// verify that all sensors have the same distortion parameters
		int numSensors = sensors.length;
		for (int i = 1; i < numSensors; i++){
			if (	(sensors[0].focalLength !=           sensors[i].focalLength) ||
					(sensors[0].distortionC !=           sensors[i].distortionC) ||
					(sensors[0].distortionB !=           sensors[i].distortionB) ||
					(sensors[0].distortionA !=           sensors[i].distortionA) ||
					(sensors[0].distortionA5 !=          sensors[i].distortionA5) ||
					(sensors[0].distortionA6 !=          sensors[i].distortionA6) ||
					(sensors[0].distortionA7 !=          sensors[i].distortionA7) ||
					(sensors[0].distortionA8 !=          sensors[i].distortionA8) ||
					(sensors[0].distortionRadius !=      sensors[i].distortionRadius) ||
					(sensors[0].pixelCorrectionWidth !=  sensors[i].pixelCorrectionWidth) ||
					(sensors[0].pixelCorrectionHeight != sensors[i].pixelCorrectionHeight) ||
					(sensors[0].pixelSize !=             sensors[i].pixelSize)){
				System.out.println("initGeometryCorrection(): All sensors have to have the same distortion model, but channels 0 and "+i+" mismatch");
				return false;
			}
		}
		// set common distportion parameters
		geometryCorrection.setDistortion(
				sensors[0].focalLength,
				sensors[0].distortionC,
				sensors[0].distortionB,
				sensors[0].distortionA,
				sensors[0].distortionA5,
				sensors[0].distortionA6,
				sensors[0].distortionA7,
				sensors[0].distortionA8,
				sensors[0].distortionRadius,
				sensors[0].pixelCorrectionWidth,   // virtual camera center is at (pixelCorrectionWidth/2, pixelCorrectionHeight/2)
				sensors[0].pixelCorrectionHeight,
				sensors[0].pixelSize);
		// set other/individual sensor parameters
		for (int i = 1; i < numSensors; i++){
			if (	(sensors[0].theta !=                 sensors[i].theta) || // elevation
					(sensors[0].heading !=               sensors[i].heading)){
				System.out.println("initGeometryCorrection(): All sensors have to have the same elevation and heading, but channels 0 and "+i+" mismatch");
				return false;
			}
		}
		double []   forward = new double[numSensors];
		double []   right =   new double[numSensors];
		double []   height =  new double[numSensors];
		double []   roll =    new double[numSensors];
		double [][] pXY0 =    new double[numSensors][2];
		for (int i = 0; i < numSensors; i++){
			forward[i] = sensors[i].forward;
			right[i] =   sensors[i].right;
			height[i] =  sensors[i].height;
			roll[i] =    sensors[i].psi;
			pXY0[i][0] = sensors[i].px0;
			pXY0[i][1] = sensors[i].py0;
		}
		geometryCorrection.setSensors(
				numSensors,
				sensors[0].theta,
				sensors[0].heading,
				forward,
				right,
				height,
				roll, 
				pXY0);
		geometryCorrection.planeProjectLenses(); // project all lenses to the common plane 
		
		// calcualte reverse distortion as a table to be linear intr4epolated
		geometryCorrection.calcReverseDistortionTable();
		
		if (numSensors == 4){
			geometryCorrection.adustSquare();
			System.out.println("Adjusted camera to orient X Y along the sides of a square");
		} else {
			System.out.println("============= Cannot adustSquare() as it requires exactly 4 sensors, "+numSensors+" provided ==========");
			return false;
		}
		// Print parameters
		if (debugLevel > 0){
			geometryCorrection.listGeometryCorrection(debugLevel > 1);
		}
		
//listGeometryCorrection		
		return true;
	}
	
//GeometryCorrection

	  public double [][][][][] calculateCLTKernel ( // per color/per tileY/ per tileX/per quadrant (plus offset as 5-th)/per pixel
			  final ImageStack kernelStack,  // first stack with 3 colors/slices convolution kernels
			  final int          kernelSize, // 64
			  final EyesisCorrectionParameters.CLTParameters clt_parameters,

			  final int          threadsMax,  // maximal number of threads to launch                         
			  final boolean      updateStatus,
			  final int          globalDebugLevel) // update status info
	  {
		  if (kernelStack==null) return null;
		  final int kernelWidth=kernelStack.getWidth();
		  final int kernelNumHor=kernelWidth/kernelSize;
		  final int kernelNumVert=kernelStack.getHeight()/kernelSize;
		  final int nChn=kernelStack.getSize();
		  final int dtt_size =      clt_parameters.transform_size; 
		  final int dtt_len = dtt_size* dtt_size;  
		  final double [][][][][] clt_kernels = new double [nChn][kernelNumVert][kernelNumHor][5][];
		  for (int chn = 0; chn < nChn; chn++){
			  for (int tileY = 0; tileY < kernelNumVert ; tileY++){
				  for (int tileX = 0; tileX < kernelNumHor ; tileX++){
					  for (int n = 0; n<4; n++){
						  clt_kernels[chn][tileY][tileX][n] = new double [dtt_len];
					  }
					  clt_kernels[chn][tileY][tileX][4] = new double [extra_items];
				  }
			  }
		  }
		  // currently each 64x64 kernel corresponds to 16x16 original pixels tile, 2 tiles margin each side
		  final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		  final AtomicInteger ai = new AtomicInteger(0);
		  final int numberOfKernels=     kernelNumHor*kernelNumVert*nChn;
		  final int numberOfKernelsInChn=kernelNumHor*kernelNumVert;
		  final double [] norm_sym_weights = clt_parameters.norm_kern ? new double [dtt_size*dtt_size]:null;
		  if (norm_sym_weights != null) {
			  for (int i = 0; i < dtt_size; i++){
				  for (int j = 0; j < dtt_size; j++){
					  norm_sym_weights[i*dtt_size+j] = Math.cos(Math.PI*i/(2*dtt_size))*Math.cos(Math.PI*j/(2*dtt_size));
				  }
			  }
		  }
		  
		  final long startTime = System.nanoTime();
		  System.out.println("calculateCLTKernel():numberOfKernels="+numberOfKernels);
		  for (int ithread = 0; ithread < threads.length; ithread++) {
			  threads[ithread] = new Thread() {
				  public void run() {
					  float [] kernelPixels= null; // will be initialized at first use NOT yet?
					  double [] kernel=      new double[kernelSize*kernelSize];
					  int centered_len = (2*dtt_size-1) * (2*dtt_size-1);
					  double [] kernel_centered = new double [centered_len + extra_items];
					  ImageDtt image_dtt = new ImageDtt();
					  int chn,tileY,tileX;
					  DttRad2 dtt = new DttRad2(dtt_size);
					  showDoubleFloatArrays sdfa_instance = null;
					  if (globalDebugLevel > -1) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?

					  for (int nTile = ai.getAndIncrement(); nTile < numberOfKernels; nTile = ai.getAndIncrement()) {
						  chn=nTile/numberOfKernelsInChn;
						  tileY =(nTile % numberOfKernelsInChn)/kernelNumHor;
						  tileX = nTile % kernelNumHor;
						  if (tileX==0) {
							  if (updateStatus) IJ.showStatus("Processing kernels, channel "+(chn+1)+" of "+nChn+", row "+(tileY+1)+" of "+kernelNumVert);
							  if (globalDebugLevel>2) System.out.println("Processing kernels, channel "+(chn+1)+" of "+nChn+", row "+(tileY+1)+" of "+kernelNumVert+" : "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
						  }
						  kernelPixels=(float[]) kernelStack.getPixels(chn+1);
						  
						  /* read convolution kernel */
						  extractOneKernel(
								  kernelPixels, //  array of combined square kernels, each 
								  kernel, // will be filled, should have correct size before call
								  kernelNumHor, // number of kernels in a row
								  tileX, // horizontal number of kernel to extract
								  tileY); // vertical number of kernel to extract
						  if ((globalDebugLevel > 0) && (tileY == clt_parameters.tileY/2)  && (tileX == clt_parameters.tileX/2)) {
							  int length=kernel.length;
							  int size=(int) Math.sqrt(length);
							  double s =0.0;
							  for (int i=0;i<kernel.length;i++) s+=kernel[i];
							  System.out.println("calculateCLTKernel(): sum(kernel_raw)="+s);
							  sdfa_instance.showArrays(
									  kernel,
									  size,
									  size,
									  "raw_kernel-"+chn+"-X"+(clt_parameters.tileX/2)+"-Y"+(clt_parameters.tileY/2));
							  
						  }
						  
						  // now has 64x64
						  image_dtt.clt_convert_double_kernel( // converts double resolution kernel
								  kernel,          // double []   src_kernel, //
								  kernel_centered, // double []   dst_kernel, // should be (2*dtt_size-1) * (2*dtt_size-1) + extra_items size - kernel and dx, dy to the nearest 1/2 pixels
								                   // also actual full center shifts in sensor pixels
								  kernelSize,      // int src_size, // 64
								  dtt_size);       // 8
						  if ((globalDebugLevel > 0) && (tileY == clt_parameters.tileY/2)  && (tileX == clt_parameters.tileX/2)) {
							  int length=kernel_centered.length;
							  int size=(int) Math.sqrt(length);
							  double s =0.0;
							  int klen = (2*dtt_size-1) * (2*dtt_size-1);
							  for (int i = 0; i < klen; i++) s += kernel_centered[i];
							  System.out.println("calculateCLTKernel(): sum(kernel_centered)="+s);
							  sdfa_instance.showArrays(
									  kernel_centered,
									  size,
									  size,
									  "kernel_centered-"+chn+"-X"+(clt_parameters.tileX/2)+"-Y"+(clt_parameters.tileY/2));
						  }
						  
						  if (norm_sym_weights != null) {
							  image_dtt.clt_normalize_kernel( // 
									  kernel_centered, // double []   kernel, // should be (2*dtt_size-1) * (2*dtt_size-1) + extra_items size (last (2*dtt_size-1) are not modified)
									  norm_sym_weights, // double []   window, // normalizes result kernel * window to have sum of elements == 1.0 
									  dtt_size,
									  (globalDebugLevel > 0) && (tileY == clt_parameters.tileY/2)  && (tileX == clt_parameters.tileX/2)); // 8
							  if ((globalDebugLevel > 0) && (tileY == clt_parameters.tileY/2)  && (tileX == clt_parameters.tileX/2)) {
								  int length=kernel_centered.length;
								  int size=(int) Math.sqrt(length);
								  double s =0.0;
								  int klen = (2*dtt_size-1) * (2*dtt_size-1);
								  for (int i = 0; i < klen; i++) s += kernel_centered[i];
								  System.out.println("calculateCLTKernel(): sum(kernel_normalized)="+s);
								  sdfa_instance.showArrays(
										  kernel_centered,
										  size,
										  size,
										  "kernel_normalized-"+chn+"-X"+(clt_parameters.tileX/2)+"-Y"+(clt_parameters.tileY/2));

							  }
						  }
						  image_dtt.clt_symmetrize_kernel( // 
								  kernel_centered, // double []     kernel,      // should be (2*dtt_size-1) * (2*dtt_size-1) +4 size (last 4 are not modified)
								  clt_kernels[chn][tileY][tileX], // 	double [][]   sym_kernels, // set of 4 SS, AS, SA, AA kdernels, each dtt_size * dtt_size (may have 5-th with center shift  
								  dtt_size); // 8
						  for (int i = 0; i < extra_items; i++){
							  clt_kernels[chn][tileY][tileX][4][i] = kernel_centered [centered_len + i];
						  }
						  if ((globalDebugLevel > 0) && (tileY == clt_parameters.tileY/2)  && (tileX == clt_parameters.tileX/2)) {
							  double [][] dbg_clt = {
									  clt_kernels[chn][tileY][tileX][0],
									  clt_kernels[chn][tileY][tileX][1],
									  clt_kernels[chn][tileY][tileX][2],
									  clt_kernels[chn][tileY][tileX][3]};
							  String [] titles = {"CC", "SC", "CS", "SS"};
							  int length=dbg_clt[0].length;
							  int size=(int) Math.sqrt(length);
							  sdfa_instance.showArrays(
									  dbg_clt,
									  size,
									  size,
									  true,
									  "pre_clt_kernels-"+chn,
									  titles);
						  }						  
						  image_dtt.clt_dtt3_kernel( // 
								  clt_kernels[chn][tileY][tileX], // double [][]   kernels, // set of 4 SS, AS, SA, AA kdernels, each dtt_size * dtt_size (may have 5-th with center shift  
								  dtt_size, // 8
								  dtt);
						  if ((globalDebugLevel > 0) && (tileY == clt_parameters.tileY/2)  && (tileX == clt_parameters.tileX/2)) {
							  double [][] dbg_clt = {
									  clt_kernels[chn][tileY][tileX][0],
									  clt_kernels[chn][tileY][tileX][1],
									  clt_kernels[chn][tileY][tileX][2],
									  clt_kernels[chn][tileY][tileX][3]};
							  String [] titles = {"CC", "SC", "CS", "SS"};
							  int length=dbg_clt[0].length;
							  int size=(int) Math.sqrt(length);
							  sdfa_instance.showArrays(
									  dbg_clt,
									  size,
									  size,
									  true,
									  "dbg_clt_kernels-"+chn,
									  titles);
							  System.out.println("calculateCLTKernel() chn="+chn+" "+
									  "tileX = "+clt_parameters.tileX+" ("+(clt_parameters.tileX/2)+") "+
									  "tileY = "+clt_parameters.tileY+" ("+(clt_parameters.tileY/2)+") "+
									  "center_x = "+clt_kernels[chn][tileY][tileX][4][0]+", "+
									  "center_y = "+clt_kernels[chn][tileY][tileX][4][1]+", "+
									  "full_dx =  "+clt_kernels[chn][tileY][tileX][4][2]+", "+
									  "full_dy =  "+clt_kernels[chn][tileY][tileX][4][3]);
						  }						  
					  }
				  }
			  };
		  }		      
		  ImageDtt.startAndJoin(threads);
		  if (globalDebugLevel > 1) System.out.println("Threads done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		  System.out.println("1.Threads done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		  // Calculate differential offsets to interpolate for tiles between kernel centers
		  ImageDtt image_dtt = new ImageDtt();
		  image_dtt.clt_fill_coord_corr(
				  clt_parameters.kernel_step,  //  final int             kern_step, // distance between kernel centers, in pixels.
				  clt_kernels,                 // final double [][][][] clt_data,
				  threadsMax,                  // maximal number of threads to launch                         
				  globalDebugLevel);
		  return clt_kernels;
	  }

	  public double [][] flattenCLTKernels (      // per color, save 4 kernelas and displacement as (2*dtt_size+1)*(2*dtt_size) tiles in an image (last row - 4 values shift x,y)
			  final double [][][][][] kernels,    // per color/per tileY/ per tileX/per quadrant (plus offset as 5-th)/per pixel
			  final int          threadsMax,      // maximal number of threads to launch                         
			  final boolean      updateStatus,
			  final int          globalDebugLevel) // update status info
	  {
		  if (kernels==null) return null;
		  final int nChn = kernels.length;
		  final int kernelNumVert=kernels[0].length;
		  final int kernelNumHor=kernels[0][0].length;
		  final int dtt_len = kernels[0][0][0][0].length;  
		  final int dtt_size =      (int) Math.sqrt(dtt_len);
		  final int tileWidth =  2 * dtt_size;
		  final int tileHeight = 2 * dtt_size + 1; // last row - shift with 0.5 pix steps
		  final int width =  tileWidth *  kernelNumHor;
		  final int height = tileHeight * kernelNumVert;
		  final double [][] clt_flat = new double [nChn][width * height];
		  // currently each 64x64 kernel corresponds to 16x16 original pixels tile, 2 tiles margin each side
		  final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		  final AtomicInteger ai = new AtomicInteger(0);
		  final int numberOfKernels=     kernelNumHor*kernelNumVert*nChn;
		  final int numberOfKernelsInChn=kernelNumHor*kernelNumVert;
		  
		  final long startTime = System.nanoTime();
		  System.out.println("flattenCLTKernels():numberOfKernels="+numberOfKernels);
		  for (int ithread = 0; ithread < threads.length; ithread++) {
			  threads[ithread] = new Thread() {
				  public void run() {
					  int chn,tileY,tileX;
					  for (int nTile = ai.getAndIncrement(); nTile < numberOfKernels; nTile = ai.getAndIncrement()) {
						  chn=nTile/numberOfKernelsInChn;
						  tileY =(nTile % numberOfKernelsInChn)/kernelNumHor;
						  tileX = nTile % kernelNumHor;
						  for (int i = 0; i < dtt_size; i++){
							System.arraycopy(
									kernels[chn][tileY][tileX][0],
									i * dtt_size,
									clt_flat[chn],
									(tileY*tileHeight + i)            * width + (tileX * tileWidth),
									dtt_size);
							System.arraycopy(
									kernels[chn][tileY][tileX][1],
									i * dtt_size,
									clt_flat[chn],
									(tileY*tileHeight + i)            * width + (tileX * tileWidth) + dtt_size,
									dtt_size);

							System.arraycopy(
									kernels[chn][tileY][tileX][2],
									i * dtt_size,
									clt_flat[chn],
									(tileY*tileHeight + i + dtt_size) * width + (tileX * tileWidth),
									dtt_size);
							System.arraycopy(
									kernels[chn][tileY][tileX][3],
									i * dtt_size,
									clt_flat[chn],
									(tileY*tileHeight + i + dtt_size) * width + (tileX * tileWidth) + 1 * dtt_size,
									dtt_size);
						  }
						  System.arraycopy(
								  kernels[chn][tileY][tileX][4], // just 2 values
								  0,
								  clt_flat[chn],
								  (tileY*tileHeight + 2 * dtt_size) * width + (tileX * tileWidth),
								  extra_items);
					  }
				  }
			  };
		  }		      
		  ImageDtt.startAndJoin(threads);
		  if (globalDebugLevel > 1) System.out.println("Threads done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		  System.out.println("1.Threads done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		  /* prepare result stack to return */
		  return clt_flat;
	  }

	  public void showCLTKernels(
			  final int          threadsMax,      // maximal number of threads to launch                         
			  final boolean      updateStatus,
			  final int          globalDebugLevel) // update status info
	  {
		  for (int chn=0;chn < clt_kernels.length; chn++){
			  if (clt_kernels[chn]!=null){
				  //					System.out.println("showKernels("+chn+")");
				  showCLTKernels(
						  chn,
						  threadsMax,
						  updateStatus,
						  globalDebugLevel);
			  }
		  }
	  }

	  public void showCLTKernels(
			  int chn,
			  final int          threadsMax,      // maximal number of threads to launch                         
			  final boolean      updateStatus,
			  final int          globalDebugLevel) // update status info
	  {
		  double [][] flat_kernels = flattenCLTKernels (      // per color, save 4 kernelas and displacement as (2*dtt_size+1)*(2*dtt_size) tiles in an image (last row - shift x,y)
				  clt_kernels[chn],     // per color/per tileY/ per tileX/per quadrant (plus offset as 5-th)/per pixel
				  threadsMax,      // maximal number of threads to launch                         
				  updateStatus,
				  globalDebugLevel); // update status info
		  int dtt_len = clt_kernels[chn][0][0][0][0].length;
		  int dtt_size= (int)Math.sqrt(dtt_len);
		  String [] titles = {"red", "blue", "green"}; 
		  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
		  sdfa_instance.showArrays(
				  flat_kernels,
				  clt_kernels[chn][0][0].length*(2*dtt_size),
				  clt_kernels[chn][0].length*(2*dtt_size+1),
				  true,
				  "clt_kernels-"+chn,
				  titles);
	  }
	  
 	  
	  public double [][][][][] extractCLTKernels (      // per color, save 4 kernelas and displacement as (2*dtt_size+1)*(2*dtt_size) tiles in an image (last row - shift x,y)
			  final float [][]   flat_kernels,    // per color/per tileY/ per tileX/per quadrant (plus offset as 5-th)/per pixel
			  final int          width,
			  final int          dtt_size,
			  final int          threadsMax,      // maximal number of threads to launch                         
			  final boolean      updateStatus,
			  final int          globalDebugLevel) // update status info
	  {
		  if (flat_kernels==null) return null;
		  final int nChn =       flat_kernels.length;
		  final int height =     flat_kernels[0].length/width;
		  final int tileWidth =  2 * dtt_size;
		  final int tileHeight = 2 * dtt_size + 1; // last row - shift with 0.5 pix steps
		  final int kernelNumHor =  width / tileWidth;
		  final int kernelNumVert = height / tileHeight;
		  final int dtt_len = dtt_size*dtt_size;  
		  final double [][][][][] clt_kernels = new double [nChn][kernelNumVert][kernelNumHor][5][];
		  for (int chn = 0; chn < nChn; chn++){
			  for (int tileY = 0; tileY < kernelNumVert ; tileY++){
				  for (int tileX = 0; tileX < kernelNumHor ; tileX++){
					  for (int n = 0; n<4; n++){
						  clt_kernels[chn][tileY][tileX][n] = new double [dtt_len];
					  }
					  clt_kernels[chn][tileY][tileX][4] = new double [extra_items];
				  }
			  }
		  }
		  
		  // currently each 64x64 kernel corresponds to 16x16 original pixels tile, 2 tiles margin each side
		  final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		  final AtomicInteger ai = new AtomicInteger(0);
		  final int numberOfKernels=     kernelNumHor*kernelNumVert*nChn;
		  final int numberOfKernelsInChn=kernelNumHor*kernelNumVert;
		  
		  final long startTime = System.nanoTime();
		  System.out.println("flattenCLTKernels():numberOfKernels="+numberOfKernels);
		  for (int ithread = 0; ithread < threads.length; ithread++) {
			  threads[ithread] = new Thread() {
				  public void run() {
					  int chn,tileY,tileX;
					  for (int nTile = ai.getAndIncrement(); nTile < numberOfKernels; nTile = ai.getAndIncrement()) {
						  chn=nTile/numberOfKernelsInChn;
						  tileY =(nTile % numberOfKernelsInChn)/kernelNumHor;
						  tileX = nTile % kernelNumHor;
						  for (int i = 0; i < dtt_size; i++){
							  for (int j = 0; j<dtt_size; j++){
								  int indx = i*dtt_size+j;
								  int baddr = (tileY*tileHeight + i) * width + (tileX * tileWidth) + j;
								  clt_kernels[chn][tileY][tileX][0][indx] = flat_kernels[chn][baddr];
								  clt_kernels[chn][tileY][tileX][1][indx] = flat_kernels[chn][baddr + dtt_size];
								  clt_kernels[chn][tileY][tileX][2][indx] = flat_kernels[chn][baddr + dtt_size * width];
								  clt_kernels[chn][tileY][tileX][3][indx] = flat_kernels[chn][baddr + dtt_size * width + dtt_size];
							  }
						  }
						  for (int i = 0; i < extra_items; i ++) {
							  clt_kernels[chn][tileY][tileX][4][i] = flat_kernels[chn][(tileY*tileHeight + 2 * dtt_size) * width + (tileX * tileWidth) + i];
						  }
					  }
				  }
			  };
		  }		      
		  ImageDtt.startAndJoin(threads);
		  if (globalDebugLevel > 1) System.out.println("Threads done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		  System.out.println("1.Threads done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		  /* prepare result stack to return */
		  return clt_kernels;
	  }
	  
	  
	  public boolean createCLTKernels(
			  EyesisCorrectionParameters.CLTParameters clt_parameters,				
			  int          srcKernelSize,
			  int          threadsMax,  // maximal number of threads to launch                         
			  boolean      updateStatus,
			  int          debugLevel
			  ){
		  String [] sharpKernelPaths= correctionsParameters.selectKernelChannelFiles(
				  0,  // 0 - sharp, 1 - smooth
				  eyesisCorrections.usedChannels.length, // numChannels, // number of channels
				  eyesisCorrections.debugLevel);
		  if (sharpKernelPaths==null) return false;
		  for (int i=0;i<sharpKernelPaths.length;i++){
			  System.out.println(i+":"+sharpKernelPaths[i]);
		  }
			if (clt_kernels == null){
				clt_kernels = new double[eyesisCorrections.usedChannels.length][][][][][];
				for (int chn=0;chn<eyesisCorrections.usedChannels.length;chn++){
					clt_kernels[chn] = null;
				}
			}
		  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();

		  for (int chn=0;chn<eyesisCorrections.usedChannels.length;chn++){
			  if (eyesisCorrections.usedChannels[chn] && (sharpKernelPaths[chn]!=null) && (clt_kernels[chn]==null)){
				  ImagePlus imp_kernel_sharp=new ImagePlus(sharpKernelPaths[chn]);
				  if (imp_kernel_sharp.getStackSize()<3) {
					  System.out.println("Need a 3-layer stack with kernels");
					  sharpKernelPaths[chn]=null;
					  continue;
				  }
				  ImageStack kernel_sharp_stack= imp_kernel_sharp.getStack();
				  System.out.println("debugLevel = "+debugLevel+" kernel_sharp_stack.getWidth() = "+kernel_sharp_stack.getWidth()+
						  " kernel_sharp_stack.getHeight() = "+kernel_sharp_stack.getHeight());

				  double [][][][][] kernels = calculateCLTKernel ( // per color/per tileY/ per tileX/per quadrant (plus offset as 5-th)/per pixel
						  kernel_sharp_stack,                      // final ImageStack kernelStack,  // first stack with 3 colors/slices convolution kernels
						  srcKernelSize,                           // final int          kernelSize, // 64
						  clt_parameters,                          // final EyesisCorrectionParameters.CLTParameters clt_parameters,
						  threadsMax,  // maximal number of threads to launch                         
						  updateStatus,
						  debugLevel); // update status info
				  
				  double [][] flat_kernels = flattenCLTKernels (      // per color, save 4 kernelas and displacement as (2*dtt_size+1)*(2*dtt_size) tiles in an image (last row - shift x,y)
						  kernels,    // per color/per tileY/ per tileX/per quadrant (plus offset as 5-th)/per pixel
						  threadsMax,  // maximal number of threads to launch                         
						  updateStatus,
						  debugLevel); // update status info
				  int kernelNumHor=kernels[0][0].length;
				  int dtt_len = kernels[0][0][0][0].length;  
				  int dtt_size =      (int) Math.sqrt(dtt_len);
				  int tileWidth =  2 * dtt_size;
				  int width =  tileWidth *  kernelNumHor;
				  int height = flat_kernels[0].length/width;

				  String [] layerNames = {"red_clt_kernels","blue_clt_kernels","green_clt_kernels"};
				  ImageStack cltStack = sdfa_instance.makeStack(
						  flat_kernels,
						  width,
						  height,
						  layerNames);

				  String cltPath=correctionsParameters.cltKernelDirectory+
						  Prefs.getFileSeparator()+
						  correctionsParameters.cltKernelPrefix+
						  String.format("%02d",chn)+
						  correctionsParameters.cltSuffix;
				  String msg="Saving CLT convolution kernels to "+cltPath;
				  IJ.showStatus(msg);
				  if (debugLevel>0) System.out.println(msg);
				  ImagePlus imp_clt=new ImagePlus(imp_kernel_sharp.getTitle()+"-clt",cltStack);
				  if (debugLevel > 0) {
					  imp_clt.getProcessor().resetMinAndMax();
					  imp_clt.show();
				  }
				  FileSaver fs=new FileSaver(imp_clt);
				  fs.saveAsTiffStack(cltPath);
			  }
		  }
		  return true;
	  }


	  
	  public boolean readCLTKernels(
			  EyesisCorrectionParameters.CLTParameters clt_parameters,
			  int          threadsMax,  // maximal number of threads to launch                         
			  boolean      updateStatus,
			  int          debugLevel
			  ){
		  int dtt_size = clt_parameters.transform_size;
		  String [] cltKernelPaths = correctionsParameters.selectCLTChannelFiles(
				  //					0,  // 0 - sharp, 1 - smooth
				  eyesisCorrections.usedChannels.length, // numChannels, // number of channels
				  eyesisCorrections.debugLevel);
		  if (cltKernelPaths==null) return false;

		  for (int i=0;i<cltKernelPaths.length;i++){
			  System.out.println(i+":"+cltKernelPaths[i]); // some may be null!
		  }
		  
		  if (clt_kernels == null){
			  clt_kernels = new double[eyesisCorrections.usedChannels.length][][][][][];
			  for (int chn=0;chn<eyesisCorrections.usedChannels.length;chn++){
				  clt_kernels[chn] = null;
			  }
		  }
		  showDoubleFloatArrays sdfa_instance = null;
		  if (debugLevel>0){
			  sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		  }

		  for (int chn=0;chn<eyesisCorrections.usedChannels.length;chn++){
			  if (eyesisCorrections.usedChannels[chn] && (cltKernelPaths[chn]!=null)){
				  ImagePlus imp_kernel_clt=new ImagePlus(cltKernelPaths[chn]);
				  if (imp_kernel_clt.getStackSize()<3) {
					  System.out.println("Need a 3-layer stack with symmetrical DCT kernels");
					  cltKernelPaths[chn]=null;
					  continue;
				  }

				  ImageStack kernel_clt_stack=  imp_kernel_clt.getStack();
				  if (debugLevel>0){
					  System.out.println(" kernel_clt_stack.getWidth() = "+kernel_clt_stack.getWidth()+
							  " kernel_clt_stack.getHeight() = "+kernel_clt_stack.getHeight());
				  }
				  int nColors = kernel_clt_stack.getSize();
				  float [][] flat_kernels = new float [nColors][];
				  for (int nc = 0; nc < nColors; nc++){
					  flat_kernels[nc]= (float[]) kernel_clt_stack.getPixels(nc + 1);
				  }
				  clt_kernels[chn] =  extractCLTKernels ( // per color, save 4 kernelas and displacement as (2*dtt_size+1)*(2*dtt_size) tiles in an image (last row - shift x,y)
						  flat_kernels,                   // per color/per tileY/ per tileX/per quadrant (plus offset as 5-th)/per pixel
						  kernel_clt_stack.getWidth(),    // final int          width,
						  dtt_size,
						  threadsMax,                     // maximal number of threads to launch                         
						  updateStatus,
						  debugLevel);                    // update status info
				  
				  
				  if (sdfa_instance != null){
					  for (int nc = 0; nc < clt_kernels[chn].length; nc++){
					  double [][] dbg_clt = {
							  clt_kernels[chn][nc][clt_parameters.tileY/2][clt_parameters.tileX/2][0],
							  clt_kernels[chn][nc][clt_parameters.tileY/2][clt_parameters.tileX/2][1],
							  clt_kernels[chn][nc][clt_parameters.tileY/2][clt_parameters.tileX/2][2],
							  clt_kernels[chn][nc][clt_parameters.tileY/2][clt_parameters.tileX/2][3]};
					  String [] titles = {"CC", "SC", "CS", "SS"};
					  int length=dbg_clt[0].length;
					  int size=(int) Math.sqrt(length);
					  sdfa_instance.showArrays(
							  dbg_clt,
							  size,
							  size,
							  true,
							  "dbg_clt-"+nc,
							  titles);
					  System.out.println("readCLTKernels() chn="+chn+", color="+nc+" "+
							  "tileX = "+clt_parameters.tileX+" ("+(clt_parameters.tileX/2)+") "+
							  "tileY = "+clt_parameters.tileY+" ("+(clt_parameters.tileY/2)+") "+
							  "center_x = "+clt_kernels[chn][nc][clt_parameters.tileY/2][clt_parameters.tileX/2][4][0]+", "+
							  "center_y = "+clt_kernels[chn][nc][clt_parameters.tileY/2][clt_parameters.tileX/2][4][1]);
					  }
				  }						  
			  }
		  }
		  return true;
	  }


// mostly for testing
//eyesisKernelImage	  
	  public double [] extractOneKernelFromStack(
			  final int          kernelSize, // 64
			  final int chn,
			  final int xTile, // horizontal number of kernel to extract
			  final int yTile)  // vertical number of kernel to extract
	  {
		  if (eyesisKernelImage == null) return null;
		  final ImageStack kernelStack = eyesisKernelImage.getStack();
		  return extractOneKernelFromStack(
				  kernelStack,  // first stack with 3 colors/slices convolution kernels
				  kernelSize, // 64
				  chn,
				  xTile, // horizontal number of kernel to extract
				  yTile);  // vertical number of kernel to extract
	  }
	  
	  public double [] extractOneKernelFromStack(
			  final ImageStack kernelStack,  // first stack with 3 colors/slices convolution kernels
			  final int          kernelSize, // 64
			  final int chn,
			  final int xTile, // horizontal number of kernel to extract
			  final int yTile)  // vertical number of kernel to extract
	  {
	  
		  final int kernelWidth=kernelStack.getWidth();
		  final int kernelNumHor=kernelWidth/kernelSize;
		  double [] kernel = new double [kernelSize*kernelSize];
		  extractOneKernel((float[]) kernelStack.getPixels(chn+1), //  array of combined square kernels, each 
				  kernel,        // will be filled, should have correct size before call
				  kernelNumHor,  // number of kernels in a row
				  xTile,         // horizontal number of kernel to extract
				  yTile);
		  return kernel;
	  }
	  // to be used in threaded method
	  private void extractOneKernel(float [] pixels, //  array of combined square kernels, each 
			  double [] kernel, // will be filled, should have correct size before call
			  int numHor, // number of kernels in a row
			  int xTile, // horizontal number of kernel to extract
			  int yTile) { // vertical number of kernel to extract
		  int length=kernel.length;
		  int size=(int) Math.sqrt(length);
		  int i,j;
		  int pixelsWidth=numHor*size;
		  int pixelsHeight=pixels.length/pixelsWidth;
		  int numVert=pixelsHeight/size;
		  /* limit tile numbers - effectively add margins around the known kernels */
		  if (xTile<0) xTile=0;
		  else if (xTile>=numHor) xTile=numHor-1;
		  if (yTile<0) yTile=0;
		  else if (yTile>=numVert) yTile=numVert-1;
		  int base=(yTile*pixelsWidth+xTile)*size;
		  for (i=0;i<size;i++) for (j=0;j<size;j++) kernel [i*size+j]=pixels[base+i*pixelsWidth+j];
	  }

	  public double [] reformatKernel(
			  double [] src_kernel,// will be blured in-place
			  int       src_size,  // typical 64
			  int       dst_size,  // typical 15 // destination size
			  int       decimation,// typical 2
			  double    sigma)
	  {
		  double [] dst_kernel = new double [dst_size*dst_size];
		  DoubleGaussianBlur gb = null;
		  if (sigma > 0) gb = new DoubleGaussianBlur();
		  reformatKernel(
				  src_kernel,
				  dst_kernel,
				  src_size,
				  dst_size,
				  decimation,
				  sigma,
				  gb);
		  return dst_kernel;

	  }
	  // to be used in threaded method
	  private void reformatKernel(
			  double [] src_kernel, // will be blured in-place
			  double [] dst_kernel,
			  int src_size,
			  int dst_size,
			  int decimation,
			  double sigma,
			  DoubleGaussianBlur gb)
	  {
		  if (gb != null) gb.blurDouble(src_kernel, src_size, src_size, sigma, sigma, 0.01);
		  int src_center = src_size / 2; // 32
		  int dst_center = dst_size / 2; // 7
		  for (int i = 0; i< dst_size; i++){
			  int src_i = (i - dst_center)*decimation + src_center;
			  if ((src_i >= 0) && (src_i < src_size)) {
				  for (int j = 0; j< dst_size; j++) {
					  int src_j = (j - dst_center)*decimation + src_center;
					  if ((src_j >= 0) && (src_j < src_size)) {
						  dst_kernel[i*dst_size + j] = src_kernel[src_i*src_size + src_j];
					  } else {
						  dst_kernel[i*dst_size + j] = 0;
					  }
				  }
			  } else {
				  for (int j = 0; j< dst_size; j++) dst_kernel[i*dst_size + j] = 0;
			  }
		  }
	  }
	  public double []reformatKernel2( // averages by exactly 2 (decimate==2)
			  double [] src_kernel, //
			  int src_size,
			  int dst_size){
		  double [] dst_kernel = new double [dst_size*dst_size];
		  reformatKernel2(
				  src_kernel,
				  dst_kernel,
				  src_size,
				  dst_size);
		  return dst_kernel;
	  }

	  private void reformatKernel2( // averages by exactly 2 (decimate==2)
			  double [] src_kernel, //
			  double [] dst_kernel,
			  int src_size,
			  int dst_size)
	  {
		  int decimation = 2;
		  int [] indices = {0,-src_size,-1,1,src_size,-src_size-1,-src_size+1,src_size-1,src_size+1};
		  double [] weights = {0.25,0.125,0.125,0.125,0.125,0.0625,0.0625,0.0625,0.0625};
		  int src_center = src_size / 2; // 32
		  int dst_center = dst_size / 2; // 7
		  int src_len = src_size*src_size;  
		  for (int i = 0; i< dst_size; i++){
			  int src_i = (i - dst_center)*decimation + src_center;
			  if ((src_i >= 0) && (src_i < src_size)) {
				  for (int j = 0; j< dst_size; j++) {
					  int src_j = (j - dst_center)*decimation + src_center;
					  int dst_index = i*dst_size + j;
					  dst_kernel[dst_index] = 0.0;
					  if ((src_j >= 0) && (src_j < src_size)) {
						  int src_index = src_i*src_size + src_j;
						  for (int k = 0; k < indices.length; k++){
							  int indx = src_index + indices[k]; // normally source kernel should be larger, these lines just to save from "out of bounds"
							  if      (indx < 0)        indx += src_len;
							  else if (indx >= src_len) indx -= src_len;
							  dst_kernel[dst_index] += weights[k]*src_kernel[indx];
						  }
					  }
				  }
			  } else {
				  for (int j = 0; j< dst_size; j++) dst_kernel[i*dst_size + j] = 0;
			  }
		  }
	  }
	  
	  public void resetCLTKernels() // and geometry corection too
	  {
		  clt_kernels = null;
		  geometryCorrection=null;

	  }

	  public ImageStack  YPrPbToRGB(double [][] yPrPb,
			  double Kr,        // 0.299;
			  double Kb,        // 0.114;
			  int width
			  ) {
		  int length = yPrPb[0].length;
		  int height = length/width;

		  float [] fpixels_r= new float [length];
		  float [] fpixels_g= new float [length];
		  float [] fpixels_b= new float [length];
		  double Kg=1.0-Kr-Kb;
		  int i;
		  /**
	R= Y+ Pr*2.0*(1-Kr)
	B= Y+ Pb*2.0*(1-Kb)
	G= Y  +Pr*(- 2*Kr*(1-Kr))/Kg + Pb*(-2*Kb*(1-Kb))/Kg

		   */
		  double KPrR=  2.0*(1-Kr);
		  double KPbB=  2.0*(1-Kb);
		  double KPrG= -2.0*Kr*(1-Kr)/Kg;
		  double KPbG= -2.0*Kb*(1-Kb)/Kg;
		  double Y,Pr,Pb;
		  for (i=0;i<length;i++) {
			  Pb=yPrPb[2][i];
			  Pr=yPrPb[1][i];
			  Y =yPrPb[0][i];
			  fpixels_r[i]=(float) (Y+ Pr*KPrR);
			  fpixels_b[i]=(float) (Y+ Pb*KPbB);
			  fpixels_g[i]=(float) (Y+ Pr*KPrG + Pb*KPbG);
		  }
		  ImageStack stack=new ImageStack(width,height);
		  stack.addSlice("red",   fpixels_r);
		  stack.addSlice("green", fpixels_g);
		  stack.addSlice("blue",  fpixels_b);
		  return stack;

	  }

	  public double [][]  YPrPbToRBG(double [][] yPrPb,
			  double Kr,        // 0.299;
			  double Kb,        // 0.114;
			  int width
			  ) {
		  int length = yPrPb[0].length;
		  //			int height = length/width;

		  double [][]rbg = new double[3][length];
		  double Kg=1.0-Kr-Kb;
		  int i;
		  /**
	R= Y+ Pr*2.0*(1-Kr)
	B= Y+ Pb*2.0*(1-Kb)
	G= Y  +Pr*(- 2*Kr*(1-Kr))/Kg + Pb*(-2*Kb*(1-Kb))/Kg

		   */
		  double KPrR=  2.0*(1-Kr);
		  double KPbB=  2.0*(1-Kb);
		  double KPrG= -2.0*Kr*(1-Kr)/Kg;
		  double KPbG= -2.0*Kb*(1-Kb)/Kg;
		  double Y,Pr,Pb;
		  for (i=0;i<length;i++) {
			  Pb=yPrPb[2][i];
			  Pr=yPrPb[1][i];
			  Y =yPrPb[0][i];
			  rbg[0][i]=(float) (Y+ Pr*KPrR);
			  rbg[1][i]=(float) (Y+ Pb*KPbB);
			  rbg[2][i]=(float) (Y+ Pr*KPrG + Pb*KPbG);
		  }
		  return rbg;
	  }
	  
	  public void debayer_rbg(
			  ImageStack stack_rbg){
		  debayer_rbg(stack_rbg, 1.0);
	  }

	  // Simple in-place debayer by (bi) linear approximation, assumes [0R/00], [00/B0], [G0/0G] slices
	  public void debayer_rbg(
			  ImageStack stack_rbg,
			  double scale)
	  {
		  int width = stack_rbg.getWidth();
		  int height = stack_rbg.getHeight();
		  float [] fpixels_r = (float[]) stack_rbg.getPixels(1);
		  float [] fpixels_b = (float[]) stack_rbg.getPixels(2);
		  float [] fpixels_g = (float[]) stack_rbg.getPixels(3);
		  int [][][] av_row = {
				  {{1,1},{-1,1},{-1,-1}},
				  {{1,1},{-1,1},{-1,-1}},
				  {{1,1},{-1,1},{-1,-1}}};
		  int [][][] av_col =  {
				  {{ width, width},{ width, width},{ width, width}},
				  {{-width, width},{-width, width},{-width, width}},
				  {{-width,-width},{-width,-width},{-width,-width}}};

		  int [][][] av_xcross =  {
				  {{ width+1, width+1, width+1, width+1}, { width-1, width+1, width-1, width+1}, { width-1, width-1, width-1, width-1}},
				  {{-width+1, width+1,-width+1, width+1}, {-width-1,-width+1, width-1, width+1}, {-width-1,-width-1, width-1, width-1}},
				  {{-width+1,-width+1,-width+1,-width+1}, {-width-1,-width+1,-width-1,-width+1}, {-width-1,-width-1,-width-1,-width-1}}};
		  int [][][] av_plus =  {
				  {{ width,         1,       1, width},   { width,        -1,       1, width  }, { width,        -1,      -1, width  }},
				  {{-width,         1,       1, width},   {-width,        -1,       1, width  }, {-width,        -1,      -1, width  }},
				  {{-width,         1,       1,-width},   {-width,        -1,       1,-width  }, {-width,        -1,      -1,-width  }}};
		  for (int y = 0; y < height; y++){
			  boolean odd_row = (y & 1) != 0;
			  int row_type = (y==0)? 0: ((y==(height-1))?2:1);
			  for (int x = 0; x < width; x++){
				  int indx = y*width+x;
				  boolean odd_col =   (x & 1) != 0;
				  int col_type = (x==0)? 0: ((x==(width-1))?2:1);
				  if (odd_row){
					  if (odd_col){ // GB site
						  fpixels_r[indx] = 0.5f*(
								  fpixels_r[indx+av_col[row_type][col_type][0]]+
								  fpixels_r[indx+av_col[row_type][col_type][1]]); 
						  fpixels_b[indx] = 0.5f*(
								  fpixels_b[indx+av_row[row_type][col_type][0]]+
								  fpixels_b[indx+av_row[row_type][col_type][1]]); 
					  } else { // !odd col  // B site
						  fpixels_r[indx] = 0.25f*(
								  fpixels_r[indx+av_xcross[row_type][col_type][0]]+
								  fpixels_r[indx+av_xcross[row_type][col_type][1]]+
								  fpixels_r[indx+av_xcross[row_type][col_type][2]]+
								  fpixels_r[indx+av_xcross[row_type][col_type][3]]);
						  fpixels_g[indx] = 0.25f*(
								  fpixels_g[indx+av_plus[row_type][col_type][0]]+
								  fpixels_g[indx+av_plus[row_type][col_type][1]]+
								  fpixels_g[indx+av_plus[row_type][col_type][2]]+
								  fpixels_g[indx+av_plus[row_type][col_type][3]]); 
					  }

				  } else { // !odd_row 
					  if (odd_col){  // R site
						  fpixels_b[indx] = 0.25f*(
								  fpixels_b[indx+av_xcross[row_type][col_type][0]]+
								  fpixels_b[indx+av_xcross[row_type][col_type][1]]+
								  fpixels_b[indx+av_xcross[row_type][col_type][2]]+
								  fpixels_b[indx+av_xcross[row_type][col_type][3]]);
						  fpixels_g[indx] = 0.25f*(
								  fpixels_g[indx+av_plus[row_type][col_type][0]]+
								  fpixels_g[indx+av_plus[row_type][col_type][1]]+
								  fpixels_g[indx+av_plus[row_type][col_type][2]]+
								  fpixels_g[indx+av_plus[row_type][col_type][3]]); 
					  } else { // !odd col  // G site
						  fpixels_r[indx] = 0.5f*(
								  fpixels_r[indx+av_row[row_type][col_type][0]]+
								  fpixels_r[indx+av_row[row_type][col_type][1]]); 
						  fpixels_b[indx] = 0.5f*(
								  fpixels_b[indx+av_col[row_type][col_type][0]]+
								  fpixels_b[indx+av_col[row_type][col_type][1]]); 
					  }
				  }
			  }				
		  }
		  if (scale !=1.0){
			  for (int i = 0; i< fpixels_r.length; i++){
				  fpixels_r[i] *= scale;
				  fpixels_b[i] *= scale;
				  fpixels_g[i] *= scale;
			  }
		  }
	  }

	  public void processCLTChannelImages(
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
			  int          convolveFFTSize, // 128 - fft size, kernel size should be size/2
			  final int          threadsMax,  // maximal number of threads to launch                         
			  final boolean    updateStatus,
			  final int        debugLevel)
	  {
		  this.startTime=System.nanoTime();
		  String [] sourceFiles=correctionsParameters.getSourcePaths();
		  boolean [] enabledFiles=new boolean[sourceFiles.length];
		  for (int i=0;i<enabledFiles.length;i++) enabledFiles[i]=false;
		  int numFilesToProcess=0;
		  int numImagesToProcess=0;
		  for (int nFile=0;nFile<enabledFiles.length;nFile++){
			  if ((sourceFiles[nFile]!=null) && (sourceFiles[nFile].length()>1)) {
				  int [] channels={correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile])};
				  if (correctionsParameters.isJP4()){
					  int subCamera= channels[0]- correctionsParameters.firstSubCamera; // to match those in the sensor files
					  // removeUnusedSensorData should be off!?
					  channels=this.eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (eyesisCorrections.isChannelEnabled(channels[i])){
						  if (!enabledFiles[nFile]) numFilesToProcess++;
						  enabledFiles[nFile]=true;
						  numImagesToProcess++;
					  }
				  }
			  }
		  }
		  if (numFilesToProcess==0){
			  System.out.println("No files to process (of "+sourceFiles.length+")");
			  return;
		  } else {
			  if (debugLevel>0) System.out.println(numFilesToProcess+ " files to process (of "+sourceFiles.length+"), "+numImagesToProcess+" images to process");
		  }
		  double [] referenceExposures=eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		  int [][] fileIndices=new int [numImagesToProcess][2]; // file index, channel number
		  int index=0;
		  for (int nFile=0;nFile<enabledFiles.length;nFile++){
			  if ((sourceFiles[nFile]!=null) && (sourceFiles[nFile].length()>1)) {
				  int [] channels={correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile])};
				  if (correctionsParameters.isJP4()){
					  int subCamera= channels[0]- correctionsParameters.firstSubCamera; // to match those in the sensor files
					  channels=eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (eyesisCorrections.isChannelEnabled(channels[i])){
						  fileIndices[index  ][0]=nFile;
						  fileIndices[index++][1]=channels[i];
					  }
				  }
			  }
		  }
		  for (int iImage=0;iImage<fileIndices.length;iImage++){
			  int nFile=fileIndices[iImage][0];
			  ImagePlus imp_src=null;
			  //				  int srcChannel=correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile]);
			  int srcChannel=fileIndices[iImage][1];
			  if (correctionsParameters.isJP4()){
				  int subchannel=eyesisCorrections.pixelMapping.getSubChannel(srcChannel);
				  if (this.correctionsParameters.swapSubchannels01) {
					  switch (subchannel){
					  case 0: subchannel=1; break;
					  case 1: subchannel=0; break;
					  }
				  }
				  if (debugLevel>0) System.out.println("Processing channel "+fileIndices[iImage][1]+" - subchannel "+subchannel+" of "+sourceFiles[nFile]);
				  ImagePlus imp_composite=eyesisCorrections.JP4_INSTANCE.open(
						  "", // path,
						  sourceFiles[nFile],
						  "",  //arg - not used in JP46 reader
						  true, // un-apply camera color gains
						  null, // new window
						  false); // do not show
				  imp_src=eyesisCorrections.JP4_INSTANCE.demuxImage(imp_composite, subchannel);
				  if (imp_src==null) imp_src=imp_composite; // not a composite image

				  // do we need to add any properties?				  
			  } else { 
				  imp_src=new ImagePlus(sourceFiles[nFile]);
				  //					  (new JP46_Reader_camera(false)).decodeProperiesFromInfo(imp_src); // decode existent properties from info
				  eyesisCorrections.JP4_INSTANCE.decodeProperiesFromInfo(imp_src); // decode existent properties from info
				  if (debugLevel>0) System.out.println("Processing "+sourceFiles[nFile]);
			  }
			  double scaleExposure=1.0;
			  if (!Double.isNaN(referenceExposures[nFile]) && (imp_src.getProperty("EXPOSURE")!=null)){
				  scaleExposure=referenceExposures[nFile]/Double.parseDouble((String) imp_src.getProperty("EXPOSURE"));
				  //					  imp_src.setProperty("scaleExposure", scaleExposure); // it may already have channel
				  if (debugLevel>0) System.out.println("Will scale intensity (to compensate for exposure) by  "+scaleExposure);
			  }
			  imp_src.setProperty("name",    correctionsParameters.getNameFromSourceTiff(sourceFiles[nFile]));
			  imp_src.setProperty("channel", srcChannel); // it may already have channel
			  imp_src.setProperty("path",    sourceFiles[nFile]); // it may already have channel
			  //				  ImagePlus result=processChannelImage( // returns ImagePlus, but it already should be saved/shown
			  processCLTChannelImage( // returns ImagePlus, but it already should be saved/shown
					  imp_src, // should have properties "name"(base for saving results), "channel","path"
					  clt_parameters,
					  debayerParameters,
					  nonlinParameters,
					  colorProcParameters,
					  channelGainParameters,
					  rgbParameters,
					  convolveFFTSize, // 128 - fft size, kernel size should be size/2
					  scaleExposure,
					  threadsMax,  // maximal number of threads to launch                         
					  updateStatus,
					  debugLevel);
			  // warp result (add support for different color modes)
			  if (this.correctionsParameters.equirectangular){
				  if (equirectangularParameters.clearFullMap) eyesisCorrections.pixelMapping.deleteEquirectangularMapFull(srcChannel); // save memory? //removeUnusedSensorData - no, use equirectangular specific settings
				  if (equirectangularParameters.clearAllMaps) eyesisCorrections.pixelMapping.deleteEquirectangularMapAll(srcChannel); // save memory? //removeUnusedSensorData - no, use equirectangular specific settings
			  }
			  //pixelMapping
			  Runtime.getRuntime().gc();
			  if (debugLevel >-1) System.out.println("Processing image "+(iImage+1)+" (of "+fileIndices.length+") finished at "+
					  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
			  if (eyesisCorrections.stopRequested.get()>0) {
				  System.out.println("User requested stop");
				  return;
			  }
		  }
		  System.out.println("Processing "+fileIndices.length+" files finished at "+
				  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
	  }		

	  public ImagePlus processCLTChannelImage(
			  ImagePlus imp_src, // should have properties "name"(base for saving results), "channel","path"
//			  EyesisCorrectionParameters.DCTParameters           dct_parameters,
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  int          convolveFFTSize, // 128 - fft size, kernel size should be size/2
			  double 		     scaleExposure,
			  final int        threadsMax,  // maximal number of threads to launch                         
			  final boolean    updateStatus,
			  final int        debugLevel){
		  boolean advanced=this.correctionsParameters.zcorrect || this.correctionsParameters.equirectangular;
		  boolean crop=      advanced? true: this.correctionsParameters.crop; 
		  boolean rotate=    advanced? false: this.correctionsParameters.rotate; 
		  double JPEG_scale= advanced? 1.0: this.correctionsParameters.JPEG_scale;
		  boolean toRGB=     advanced? true: this.correctionsParameters.toRGB; 
		  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?

		  // may use this.StartTime to report intermediate steps execution times
		  String name=(String) imp_src.getProperty("name");
		  //		int channel= Integer.parseInt((String) imp_src.getProperty("channel"));
		  int channel= (Integer) imp_src.getProperty("channel");
		  String path= (String) imp_src.getProperty("path");
		  if (this.correctionsParameters.pixelDefects && (eyesisCorrections.defectsXY!=null)&& (eyesisCorrections.defectsXY[channel]!=null)){
			  // apply pixel correction
			  int numApplied=	eyesisCorrections.correctDefects(
					  imp_src,
					  channel,
					  debugLevel);
			  if ((debugLevel>0) && (numApplied>0)) { // reduce verbosity after verified defect correction works
				  System.out.println("Corrected "+numApplied+" pixels in "+path);
			  }
		  }
		  if (this.correctionsParameters.vignetting){
			  if ((eyesisCorrections.channelVignettingCorrection==null) || (channel<0) || (channel>=eyesisCorrections.channelVignettingCorrection.length) || (eyesisCorrections.channelVignettingCorrection[channel]==null)){
				  System.out.println("No vignetting data for channel "+channel);
				  return null;
			  }
			  float [] pixels=(float []) imp_src.getProcessor().getPixels();
			  if (pixels.length!=eyesisCorrections.channelVignettingCorrection[channel].length){
				  System.out.println("Vignetting data for channel "+channel+" has "+eyesisCorrections.channelVignettingCorrection[channel].length+" pixels, image "+path+" has "+pixels.length);
				  return null;
			  }
			  // TODO: Move to do it once:
			  double min_non_zero = 0.0;
			  for (int i=0;i<pixels.length;i++){
				  double d = eyesisCorrections.channelVignettingCorrection[channel][i];
				  if ((d > 0.0) && ((min_non_zero == 0) || (min_non_zero > d))){
					  min_non_zero = d;
				  }
			  }
			  double max_vign_corr = clt_parameters.vignetting_range*min_non_zero;

			  System.out.println("Vignetting data: channel="+channel+", min = "+min_non_zero);
			  for (int i=0;i<pixels.length;i++){
				  double d = eyesisCorrections.channelVignettingCorrection[channel][i];
				  if (d > max_vign_corr) d = max_vign_corr;
				  pixels[i]*=d;
			  }
			  // Scale here, combine with vignetting later?
			  int width =  imp_src.getWidth();
			  int height = imp_src.getHeight();
			  for (int y = 0; y < height-1; y+=2){
				  for (int x = 0; x < width-1; x+=2){
					  pixels[y*width+x        ] *= clt_parameters.scale_g;
					  pixels[y*width+x+width+1] *= clt_parameters.scale_g;
					  pixels[y*width+x      +1] *= clt_parameters.scale_r;
					  pixels[y*width+x+width  ] *= clt_parameters.scale_b;
				  }
			  }
			  
		  } else { // assuming GR/BG pattern
			  System.out.println("Applying fixed color gain correction parameters: Gr="+
					  clt_parameters.novignetting_r+", Gg="+clt_parameters.novignetting_g+", Gb="+clt_parameters.novignetting_b);
			  float [] pixels=(float []) imp_src.getProcessor().getPixels();
			  int width =  imp_src.getWidth();
			  int height = imp_src.getHeight();
			  double kr = clt_parameters.scale_r/clt_parameters.novignetting_r;
			  double kg = clt_parameters.scale_g/clt_parameters.novignetting_g;
			  double kb = clt_parameters.scale_b/clt_parameters.novignetting_b;
			  for (int y = 0; y < height-1; y+=2){
				  for (int x = 0; x < width-1; x+=2){
					  pixels[y*width+x        ] *= kg;
					  pixels[y*width+x+width+1] *= kg;
					  pixels[y*width+x      +1] *= kr;
					  pixels[y*width+x+width  ] *= kb;
				  }
			  }
		  }
		  if (clt_parameters.gain_equalize){
			  
		  }
		  
		  String title=name+"-"+String.format("%02d", channel);
		  ImagePlus result=imp_src;
		  if (debugLevel>1) System.out.println("processing: "+path);
		  result.setTitle(title+"RAW");
		  if (!this.correctionsParameters.split){
			  eyesisCorrections.saveAndShow(result, this.correctionsParameters);
			  return result;
		  }
		  // Generate split parameters for DCT processing mode
		  EyesisCorrectionParameters.SplitParameters splitParameters = new EyesisCorrectionParameters.SplitParameters(
				                 1,  // oversample; // currently source kernels are oversampled
						  clt_parameters.transform_size/2, // addLeft
						  clt_parameters.transform_size/2, // addTop
						  clt_parameters.transform_size/2, // addRight
						  clt_parameters.transform_size/2  // addBottom
				  );		   

		  // Split into Bayer components, oversample, increase canvas
		  
		  double [][] double_stack = eyesisCorrections.bayerToDoubleStack(
				  result, // source Bayer image, linearized, 32-bit (float))
				  null); // no margins, no oversample
		  
//		  ImageStack stack= eyesisCorrections.bayerToStack(
//				  result, // source Bayer image, linearized, 32-bit (float))
//				  splitParameters);
		  String titleFull=title+"-SPLIT";
		  if (debugLevel > -1){
			  double [] chn_avg = {0.0,0.0,0.0};
			  int width =  imp_src.getWidth();
			  int height = imp_src.getHeight();
			  for (int c = 0; c < 3; c++){
				  for (int i = 0; i<double_stack[c].length; i++){
					  chn_avg[c] += double_stack[c][i];
				  }
			  }
			  chn_avg[0] /= width*height/4;
			  chn_avg[1] /= width*height/4;
			  chn_avg[2] /= width*height/2;
			  System.out.println("Split channels averages: R="+chn_avg[0]+", G="+chn_avg[2]+", B="+chn_avg[1]); 
		  }
		  String [] rbg_titles = {"Red", "Blue", "Green"};
		  ImageStack stack;
		  if (!this.correctionsParameters.debayer) {
//			  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
//			  ImageStack 
			  stack = sdfa_instance.makeStack(double_stack, imp_src.getWidth(), imp_src.getHeight(), rbg_titles);
			  result= new ImagePlus(titleFull, stack);    			  
			  eyesisCorrections.saveAndShow(result, this.correctionsParameters);
			  return result;
		  }
		  // =================
		  if (debugLevel > 0) {
			  System.out.println("Showing image BEFORE_CLT_PROC");
			  sdfa_instance.showArrays(double_stack,  imp_src.getWidth(), imp_src.getHeight(), true, "BEFORE_CLT_PROC", rbg_titles);
		  }
		  if (this.correctionsParameters.deconvolve) { // process with DCT, otherwise use simple debayer
			  ImageDtt image_dtt = new ImageDtt();
/*			  
			  double [][][][][] clt_data = image_dtt.cltStack(
					  stack,
					  channel,
					  clt_parameters,
					  clt_parameters.ishift_x, //final int shiftX, // shift image horizontally (positive - right)
					  clt_parameters.ishift_y, //final int shiftY, // shift image vertically (positive - down)
					  threadsMax,
					  debugLevel,
					  updateStatus);
*/
			  for (int i =0 ; i < double_stack[0].length; i++){
//				  double_stack[0][i]*=2.0; // Scale red twice to compensate less pixels than green
//				  double_stack[1][i]*=2.0; // Scale blue twice to compensate less pixels than green
				  double_stack[2][i]*=0.5; // Scale blue twice to compensate less pixels than green
			  }
			  double [][][][][] clt_data = image_dtt.clt_aberrations( 
					  double_stack,                 // final double [][]       imade_data,
					  imp_src.getWidth(),           //	final int               width,
					  clt_kernels[channel],         // final double [][][][][] clt_kernels, // [color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
					  clt_parameters.kernel_step,
					  clt_parameters.transform_size,
					  clt_parameters.clt_window,
					  clt_parameters.shift_x,       // final int               shiftX, // shift image horizontally (positive - right) - just for testing
					  clt_parameters.shift_y,       // final int               shiftY, // shift image vertically (positive - down)
					  clt_parameters.tileX,         // final int               debug_tileX,
					  clt_parameters.tileY,         // final int               debug_tileY,
					  (clt_parameters.dbg_mode & 64) != 0, // no fract shift
					  (clt_parameters.dbg_mode & 128) != 0, // no convolve
					  (clt_parameters.dbg_mode & 256) != 0, // transpose convolve
					  threadsMax,
					  debugLevel);
//					  updateStatus);

			  
			  System.out.println("clt_data.length="+clt_data.length+" clt_data[0].length="+clt_data[0].length
					  +" clt_data[0][0].length="+clt_data[0][0].length+" clt_data[0][0][0].length="+
					  clt_data[0][0][0].length+" clt_data[0][0][0][0].length="+clt_data[0][0][0][0].length);
			  /*
			  if (dct_parameters.color_DCT){ // convert RBG -> YPrPb
				  dct_data = image_dtt.dct_color_convert(
						  dct_data,
						  colorProcParameters.kr,
						  colorProcParameters.kb,
						  dct_parameters.sigma_rb,        // blur of channels 0,1 (r,b) in addition to 2 (g)
						  dct_parameters.sigma_y,         // blur of Y from G
						  dct_parameters.sigma_color,     // blur of Pr, Pb in addition to Y
						  threadsMax,
						  debugLevel);
			  } else { // just LPF RGB
			  */
				  if (clt_parameters.corr_sigma > 0){ // no filter at all
					  for (int chn = 0; chn < clt_data.length; chn++) {
						  image_dtt.clt_lpf(
								  clt_parameters.corr_sigma,
								  clt_data[chn],
								  clt_parameters.transform_size,
								  threadsMax,
								  debugLevel);
					  }
				  }
/*
			  }
*/
			  int tilesY = imp_src.getHeight()/clt_parameters.transform_size;
			  int tilesX = imp_src.getWidth()/clt_parameters.transform_size;
			  if (debugLevel > 0){
				  System.out.println("--tp.tilesX="+tilesX);
				  System.out.println("--tp.tilesY="+tilesY);
			  }
			  if (debugLevel > 1){
			        double [][] clt = new double [clt_data.length*4][];
			        for (int chn = 0; chn < clt_data.length; chn++) {
			        	double [][] clt_set = image_dtt.clt_dbg(
			        			clt_data [chn],
								  threadsMax,
								  debugLevel);
			        	for (int ii = 0; ii < clt_set.length; ii++) clt[chn*4+ii] = clt_set[ii];
			        }

			        if (debugLevel > 0){
			        	sdfa_instance.showArrays(clt,
			        			tilesX*clt_parameters.transform_size,
			        			tilesY*clt_parameters.transform_size,
			        			true,
			        			result.getTitle()+"-CLT");  
			        }
			  }
			  double [][] iclt_data = new double [clt_data.length][];
			  for (int chn=0; chn<clt_data.length;chn++){
				  iclt_data[chn] = image_dtt.iclt_2d(
						  clt_data[chn],                  // scanline representation of dcd data, organized as dct_size x dct_size tiles  
						  clt_parameters.transform_size,  // final int
						  clt_parameters.clt_window,      // window_type
						  15,                             // clt_parameters.iclt_mask,       //which of 4 to transform back
						  0,                              // clt_parameters.dbg_mode,        //which of 4 to transform back
						  threadsMax,
						  debugLevel);
				  
			  }
					  if (debugLevel > -1) sdfa_instance.showArrays(
							  iclt_data,
							  (tilesX + 1) * clt_parameters.transform_size,
							  (tilesY + 1) * clt_parameters.transform_size,
							  true,
							  result.getTitle()+"-rbg_sigma");
				/*	  
				  }
			  }
			 */
			  if (debugLevel > 0) sdfa_instance.showArrays(iclt_data,
					  (tilesX + 1) * clt_parameters.transform_size,
					  (tilesY + 1) * clt_parameters.transform_size,
					  true,
					  result.getTitle()+"-ICLT-RGB");

			  // convert to ImageStack of 3 slices
			  String [] sliceNames = {"red", "blue", "green"};
			  stack = sdfa_instance.makeStack(
					  iclt_data,
					  (tilesX + 1) * clt_parameters.transform_size,
					  (tilesY + 1) * clt_parameters.transform_size,
					  sliceNames); // or use null to get chn-nn slice names


		  } else { // if (this.correctionsParameters.deconvolve) - here use a simple debayer
			  System.out.println("Bypassing CLT-based aberration correction");
			  stack = sdfa_instance.makeStack(double_stack, imp_src.getWidth(), imp_src.getHeight(), rbg_titles);
			  debayer_rbg(stack, 0.25); // simple standard 3x3 kernel debayer
		  }
		  if (debugLevel > -1){
			  double [] chn_avg = {0.0,0.0,0.0};
			  float [] pixels;
			  int width =  stack.getWidth();
			  int height = stack.getHeight();
			  
			  for (int c = 0; c <3; c++){
				  pixels = (float[]) stack.getPixels(c+1);
				  for (int i = 0; i<pixels.length; i++){
					  chn_avg[c] += pixels[i];
				  }
			  }
			  chn_avg[0] /= width*height;
			  chn_avg[1] /= width*height;
			  chn_avg[2] /= width*height;
			  System.out.println("Processed channels averages: R="+chn_avg[0]+", G="+chn_avg[2]+", B="+chn_avg[1]); 
		  }
		  
		  if (!this.correctionsParameters.colorProc){
			  result= new ImagePlus(titleFull, stack);    			  
			  eyesisCorrections.saveAndShow(
					  result,
					  this.correctionsParameters);
			  return result;
		  }
		  if (debugLevel > 1) System.out.println("before colors.1");
		  //Processing colors - changing stack sequence to r-g-b (was r-b-g)
		  if (!eyesisCorrections.fixSliceSequence(
				  stack,
				  debugLevel)){
			  if (debugLevel > -1) System.out.println("fixSliceSequence() returned false");
			  return null;
		  }
		  if (debugLevel > 1) System.out.println("before colors.2");
		  if (debugLevel > 1){
			  ImagePlus imp_dbg=new ImagePlus(imp_src.getTitle()+"-"+channel+"-preColors",stack);
			  eyesisCorrections.saveAndShow(
					  imp_dbg,
					  this.correctionsParameters);
		  }
		  if (debugLevel > 1) System.out.println("before colors.3, scaleExposure="+scaleExposure+" scale = "+(255.0/eyesisCorrections.psfSubpixelShouldBe4/eyesisCorrections.psfSubpixelShouldBe4/scaleExposure));
		  CorrectionColorProc correctionColorProc=new CorrectionColorProc(eyesisCorrections.stopRequested);
		  double [][] yPrPb=new double [3][];
		  //			if (dct_parameters.color_DCT){
		  // need to get YPbPr - not RGB here				
		  //			} else {
		  correctionColorProc.processColorsWeights(stack, // just gamma convert? TODO: Cleanup? Convert directly form the linear YPrPb
				  //					  255.0/this.psfSubpixelShouldBe4/this.psfSubpixelShouldBe4, //  double scale,     // initial maximal pixel value (16))
				  //					  255.0/eyesisCorrections.psfSubpixelShouldBe4/eyesisCorrections.psfSubpixelShouldBe4/scaleExposure, //  double scale,     // initial maximal pixel value (16))
				  //					  255.0/2/2/scaleExposure, //  double scale,     // initial maximal pixel value (16))
				  255.0/scaleExposure, //  double scale,     // initial maximal pixel value (16))
				  colorProcParameters,
				  channelGainParameters,
				  channel,
				  null, //correctionDenoise.getDenoiseMask(),
				  this.correctionsParameters.blueProc,
				  debugLevel);
		  if (debugLevel > 1) System.out.println("Processed colors to YPbPr, total number of slices="+stack.getSize());
		  if (debugLevel > 1) {
			  ImagePlus imp_dbg=new ImagePlus("procColors",stack);
			  eyesisCorrections.saveAndShow(
					  imp_dbg,
					  this.correctionsParameters);
		  }
		  float [] fpixels;
		  int [] slices_YPrPb = {8,6,7};
		  yPrPb=new double [3][];
		  for (int n = 0; n < slices_YPrPb.length; n++){
			  fpixels = (float[]) stack.getPixels(slices_YPrPb[n]);
			  yPrPb[n] = new double [fpixels.length];
			  for (int i = 0; i < fpixels.length; i++) yPrPb[n][i] = fpixels[i];
		  }

		  if (toRGB) {
			  System.out.println("correctionColorProc.YPrPbToRGB");
			  stack =  YPrPbToRGB(yPrPb,
					  colorProcParameters.kr,        // 0.299;
					  colorProcParameters.kb,        // 0.114;
					  stack.getWidth());

			  title=titleFull; // including "-DECONV" or "-COMBO"
			  titleFull=title+"-RGB-float";
			  //Trim stack to just first 3 slices
			  if (debugLevel > 1){ // 2){
				  ImagePlus imp_dbg=new ImagePlus("YPrPbToRGB",stack);
				  eyesisCorrections.saveAndShow(
						  imp_dbg,
						  this.correctionsParameters);
			  }
			  while (stack.getSize() > 3) stack.deleteLastSlice();
			  if (debugLevel > 1) System.out.println("Trimming color stack");
		  } else {
			  title=titleFull; // including "-DECONV" or "-COMBO"
			  titleFull=title+"-YPrPb"; // including "-DECONV" or "-COMBO"
			  if (debugLevel > 1) System.out.println("Using full stack, including YPbPr");
		  }

		  result= new ImagePlus(titleFull, stack);    			  
		  // Crop image to match original one (scaled to oversampling)
		  if (crop){ // always crop if equirectangular
			  if (debugLevel > 1) System.out.println("cropping");
			  stack = eyesisCorrections.cropStack32(stack,splitParameters);
			  if (debugLevel > 2) { // 2){
				  ImagePlus imp_dbg=new ImagePlus("cropped",stack);
				  eyesisCorrections.saveAndShow(
						  imp_dbg,
						  this.correctionsParameters);
			  }
		  }
		  // rotate the result			  
		  if (rotate){ // never rotate for equirectangular
			  stack=eyesisCorrections.rotateStack32CW(stack);
		  }
		  if (!toRGB && !this.correctionsParameters.jpeg){ // toRGB set for equirectangular
			  if (debugLevel > 1) System.out.println("!toRGB && !this.correctionsParameters.jpeg");
			  eyesisCorrections.saveAndShow(result, this.correctionsParameters);
			  return result;
		  } else { // that's not the end result, save if required
			  if (debugLevel > 1) System.out.println("!toRGB && !this.correctionsParameters.jpeg - else");
			  eyesisCorrections.saveAndShow(result,
					  eyesisCorrections.correctionsParameters,
					  eyesisCorrections.correctionsParameters.save32,
					  false,
					  eyesisCorrections.correctionsParameters.JPEG_quality); // save, no show
		  }
		  // convert to RGB48 (16 bits per color component)
		  ImagePlus imp_RGB;
		  stack=eyesisCorrections.convertRGB32toRGB16Stack(
				  stack,
				  rgbParameters); 

		  titleFull=title+"-RGB48";
		  result= new ImagePlus(titleFull, stack);
		  //			  ImagePlus imp_RGB24;
		  result.updateAndDraw();
		  if (debugLevel > 1) System.out.println("result.updateAndDraw(), "+titleFull+"-RGB48");

		  CompositeImage compositeImage=eyesisCorrections.convertToComposite(result);

		  if (!this.correctionsParameters.jpeg && !advanced){ // RGB48 was the end result
			  if (debugLevel > 1) System.out.println("if (!this.correctionsParameters.jpeg && !advanced)");
			  eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters);
			  return result;
		  } else { // that's not the end result, save if required
			  if (debugLevel > 1) System.out.println("if (!this.correctionsParameters.jpeg && !advanced) - else");
			  eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters, this.correctionsParameters.save16, false); // save, no show
			  //				  eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters, this.correctionsParameters.save16, true); // save, no show
		  }

		  imp_RGB=eyesisCorrections.convertRGB48toRGB24(
				  stack,
				  title+"-RGB24",
				  0, 65536, // r range 0->0, 65536->256
				  0, 65536, // g range
				  0, 65536,// b range
				  0, 65536);// alpha range
		  if (JPEG_scale!=1.0){
			  ImageProcessor ip=imp_RGB.getProcessor();
			  ip.setInterpolationMethod(ImageProcessor.BICUBIC);
			  ip=ip.resize((int)(ip.getWidth()*JPEG_scale),(int) (ip.getHeight()*JPEG_scale));
			  imp_RGB= new ImagePlus(imp_RGB.getTitle(),ip);
			  imp_RGB.updateAndDraw();
		  }
		  eyesisCorrections.saveAndShow(imp_RGB, this.correctionsParameters);

		  return result;
	  }
	  
// Processing sets of 4 images together	  
	  public void processCLTSets(
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
			  int          convolveFFTSize, // 128 - fft size, kernel size should be size/2
			  final int          threadsMax,  // maximal number of threads to launch                         
			  final boolean    updateStatus,
			  final int        debugLevel)
	  {
		  this.startTime=System.nanoTime();
		  String [] sourceFiles=correctionsParameters.getSourcePaths();
		  boolean [] enabledFiles=new boolean[sourceFiles.length];
		  for (int i=0;i<enabledFiles.length;i++) enabledFiles[i]=false;
		  int numFilesToProcess=0;
		  int numImagesToProcess=0;
		  for (int nFile=0;nFile<enabledFiles.length;nFile++){
			  if ((sourceFiles[nFile]!=null) && (sourceFiles[nFile].length()>1)) {
				  int [] channels={correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile])};
				  if (correctionsParameters.isJP4()){
					  int subCamera= channels[0]- correctionsParameters.firstSubCamera; // to match those in the sensor files
					  // removeUnusedSensorData should be off!?
					  channels=this.eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (eyesisCorrections.isChannelEnabled(channels[i])){
						  if (!enabledFiles[nFile]) numFilesToProcess++;
						  enabledFiles[nFile]=true;
						  numImagesToProcess++;
					  }
				  }
			  }
		  }
		  if (numFilesToProcess==0){
			  System.out.println("No files to process (of "+sourceFiles.length+")");
			  return;
		  } else {
			  if (debugLevel>0) System.out.println(numFilesToProcess+ " files to process (of "+sourceFiles.length+"), "+numImagesToProcess+" images to process");
		  }
		  double [] referenceExposures=eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		  int [][] fileIndices=new int [numImagesToProcess][2]; // file index, channel number
		  int index=0;
		  for (int nFile=0;nFile<enabledFiles.length;nFile++){
			  if ((sourceFiles[nFile]!=null) && (sourceFiles[nFile].length()>1)) {
				  int [] channels={correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile])};
				  if (correctionsParameters.isJP4()){
					  int subCamera= channels[0]- correctionsParameters.firstSubCamera; // to match those in the sensor files
					  channels=eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (eyesisCorrections.isChannelEnabled(channels[i])){
						  fileIndices[index  ][0]=nFile;
						  fileIndices[index++][1]=channels[i];
					  }
				  }
			  }
		  }
		  ArrayList<String> setNames = new ArrayList<String>();
		  ArrayList<ArrayList<Integer>> setFiles = new ArrayList<ArrayList<Integer>>();

		  for (int iImage=0;iImage<fileIndices.length;iImage++){
			  int nFile=fileIndices[iImage][0];
			  String setName = correctionsParameters.getNameFromSourceTiff(sourceFiles[nFile]);
			  if (!setNames.contains(setName)) {
				  setNames.add(setName);
				  setFiles.add(new ArrayList<Integer>());
			  }
			  setFiles.get(setNames.indexOf(setName)).add(new Integer(nFile));
		  }		  
		  int iImage = 0;
		  for (int nSet = 0; nSet < setNames.size(); nSet++){
			  int maxChn = 0;
			  for (int i = 0; i < setFiles.get(nSet).size(); i++){
				  int chn = fileIndices[setFiles.get(nSet).get(i)][1];
				  if (chn > maxChn) maxChn = chn;
			  }
			  int [] channelFiles = new int[maxChn+1];
			  for (int i =0; i < channelFiles.length; i++) channelFiles[i] = -1;
			  for (int i = 0; i < setFiles.get(nSet).size(); i++){
				  channelFiles[fileIndices[setFiles.get(nSet).get(i)][1]] = setFiles.get(nSet).get(i);
			  }
			  
			  ImagePlus [] imp_srcs = new ImagePlus[channelFiles.length];
			  double [] scaleExposure = new double[channelFiles.length];
			  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
				  int nFile=channelFiles[srcChannel];
				  imp_srcs[srcChannel]=null;
				  if (nFile >=0){
					  if (correctionsParameters.isJP4()){
						  int subchannel=eyesisCorrections.pixelMapping.getSubChannel(srcChannel);
						  if (this.correctionsParameters.swapSubchannels01) {
							  switch (subchannel){
							  case 0: subchannel=1; break;
							  case 1: subchannel=0; break;
							  }
						  }
						  if (debugLevel>0) System.out.println("Processing set " + setNames.get(nSet)+" channel "+srcChannel+" - subchannel "+subchannel+" of "+sourceFiles[nFile]);
						  ImagePlus imp_composite=eyesisCorrections.JP4_INSTANCE.open(
								  "", // path,
								  sourceFiles[nFile],
								  "",  //arg - not used in JP46 reader
								  true, // un-apply camera color gains
								  null, // new window
								  false); // do not show
						  imp_srcs[srcChannel]=eyesisCorrections.JP4_INSTANCE.demuxImage(imp_composite, subchannel);
						  if (imp_srcs[srcChannel] == null) imp_srcs[srcChannel] = imp_composite; // not a composite image
						  // do we need to add any properties?				  
					  } else { 
						  imp_srcs[srcChannel]=new ImagePlus(sourceFiles[nFile]);
						  //					  (new JP46_Reader_camera(false)).decodeProperiesFromInfo(imp_src); // decode existent properties from info
						  eyesisCorrections.JP4_INSTANCE.decodeProperiesFromInfo(imp_srcs[srcChannel]); // decode existent properties from info
						  if (debugLevel>0) System.out.println("Processing "+sourceFiles[nFile]);
					  }
					  scaleExposure[srcChannel] = 1.0;
					  if (!Double.isNaN(referenceExposures[nFile]) && (imp_srcs[srcChannel].getProperty("EXPOSURE")!=null)){
						  scaleExposure[srcChannel] = referenceExposures[nFile]/Double.parseDouble((String) imp_srcs[srcChannel].getProperty("EXPOSURE"));
						  if (debugLevel>0) System.out.println("Will scale intensity (to compensate for exposure) by  "+scaleExposure);
					  }
					  imp_srcs[srcChannel].setProperty("name",    correctionsParameters.getNameFromSourceTiff(sourceFiles[nFile]));
					  imp_srcs[srcChannel].setProperty("channel", srcChannel); // it may already have channel
					  imp_srcs[srcChannel].setProperty("path",    sourceFiles[nFile]); // it may already have channel

					  if (this.correctionsParameters.pixelDefects && (eyesisCorrections.defectsXY!=null)&& (eyesisCorrections.defectsXY[srcChannel]!=null)){
						  // apply pixel correction
						  int numApplied=	eyesisCorrections.correctDefects(
								  imp_srcs[srcChannel],
								  srcChannel,
								  debugLevel);
						  if ((debugLevel>0) && (numApplied>0)) { // reduce verbosity after verified defect correction works
							  System.out.println("Corrected "+numApplied+" pixels in "+sourceFiles[nFile]);
						  }
					  }

					  if (this.correctionsParameters.vignetting){
						  if ((eyesisCorrections.channelVignettingCorrection==null) || (srcChannel<0) || (srcChannel>=eyesisCorrections.channelVignettingCorrection.length) || (eyesisCorrections.channelVignettingCorrection[srcChannel]==null)){
							  System.out.println("No vignetting data for channel "+srcChannel);
							  return;
						  }
						  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
						  if (pixels.length!=eyesisCorrections.channelVignettingCorrection[srcChannel].length){
							  System.out.println("Vignetting data for channel "+srcChannel+" has "+eyesisCorrections.channelVignettingCorrection[srcChannel].length+" pixels, image "+sourceFiles[nFile]+" has "+pixels.length);
							  return;
						  }
						  // TODO: Move to do it once:
						  double min_non_zero = 0.0;
						  for (int i=0;i<pixels.length;i++){
							  double d = eyesisCorrections.channelVignettingCorrection[srcChannel][i];
							  if ((d > 0.0) && ((min_non_zero == 0) || (min_non_zero > d))){
								  min_non_zero = d;
							  }
						  }
						  double max_vign_corr = clt_parameters.vignetting_range*min_non_zero;

						  System.out.println("Vignetting data: channel="+srcChannel+", min = "+min_non_zero);
						  for (int i=0;i<pixels.length;i++){
							  double d = eyesisCorrections.channelVignettingCorrection[srcChannel][i];
							  if (d > max_vign_corr) d = max_vign_corr;
							  pixels[i]*=d;
						  }
						  // Scale here, combine with vignetting later?
						  int width =  imp_srcs[srcChannel].getWidth();
						  int height = imp_srcs[srcChannel].getHeight();
						  for (int y = 0; y < height-1; y+=2){
							  for (int x = 0; x < width-1; x+=2){
								  pixels[y*width+x        ] *= clt_parameters.scale_g;
								  pixels[y*width+x+width+1] *= clt_parameters.scale_g;
								  pixels[y*width+x      +1] *= clt_parameters.scale_r;
								  pixels[y*width+x+width  ] *= clt_parameters.scale_b;
							  }
						  }

					  } else { // assuming GR/BG pattern
						  System.out.println("Applying fixed color gain correction parameters: Gr="+
								  clt_parameters.novignetting_r+", Gg="+clt_parameters.novignetting_g+", Gb="+clt_parameters.novignetting_b);
						  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
						  int width =  imp_srcs[srcChannel].getWidth();
						  int height = imp_srcs[srcChannel].getHeight();
						  double kr = clt_parameters.scale_r/clt_parameters.novignetting_r;
						  double kg = clt_parameters.scale_g/clt_parameters.novignetting_g;
						  double kb = clt_parameters.scale_b/clt_parameters.novignetting_b;
						  for (int y = 0; y < height-1; y+=2){
							  for (int x = 0; x < width-1; x+=2){
								  pixels[y*width+x        ] *= kg;
								  pixels[y*width+x+width+1] *= kg;
								  pixels[y*width+x      +1] *= kr;
								  pixels[y*width+x+width  ] *= kb;
							  }
						  }
					  }
				  }
			  }
			  // may need to equalize gains between channels
			  // may need to equalize gains between channels
			  if (clt_parameters.gain_equalize || clt_parameters.colors_equalize){
				  channelGainsEqualize(
						  clt_parameters.gain_equalize,
						  clt_parameters.colors_equalize,
						  channelFiles,
						  imp_srcs,
						  setNames.get(nSet), // just for debug messeges == setNames.get(nSet)
						  debugLevel);
			  }
			  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
				  int nFile=channelFiles[srcChannel];
				  if (nFile >=0){
					  processCLTSetImage( // returns ImagePlus, but it already should be saved/shown
							  imp_srcs[srcChannel], // should have properties "name"(base for saving results), "channel","path"
							  clt_parameters,
							  debayerParameters,
							  nonlinParameters,
							  colorProcParameters,
							  channelGainParameters,
							  rgbParameters,
							  convolveFFTSize, // 128 - fft size, kernel size should be size/2
							  scaleExposure[srcChannel],
							  threadsMax,  // maximal number of threads to launch                         
							  updateStatus,
							  debugLevel);
					  // warp result (add support for different color modes)
					  if (this.correctionsParameters.equirectangular){
						  if (equirectangularParameters.clearFullMap) eyesisCorrections.pixelMapping.deleteEquirectangularMapFull(srcChannel); // save memory? //removeUnusedSensorData - no, use equirectangular specific settings
						  if (equirectangularParameters.clearAllMaps) eyesisCorrections.pixelMapping.deleteEquirectangularMapAll(srcChannel); // save memory? //removeUnusedSensorData - no, use equirectangular specific settings
					  }
					  //pixelMapping
					  Runtime.getRuntime().gc();
					  if (debugLevel >-1) System.out.println("Processing image "+(iImage+1)+" (of "+fileIndices.length+") finished at "+
							  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
					  if (eyesisCorrections.stopRequested.get()>0) {
						  System.out.println("User requested stop");
						  return;
					  }
					  iImage++;
				  }
			  }
		  }
		  System.out.println("Processing "+fileIndices.length+" files finished at "+
				  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
	  }		
	  
	  public ImagePlus processCLTSetImage(
			  ImagePlus imp_src, // should have properties "name"(base for saving results), "channel","path"
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  int          convolveFFTSize, // 128 - fft size, kernel size should be size/2
			  double 		     scaleExposure,
			  final int        threadsMax,  // maximal number of threads to launch                         
			  final boolean    updateStatus,
			  final int        debugLevel){
		  boolean advanced=this.correctionsParameters.zcorrect || this.correctionsParameters.equirectangular;
//		  boolean crop=      advanced? true: this.correctionsParameters.crop; 
		  boolean rotate=    advanced? false: this.correctionsParameters.rotate; 
		  double JPEG_scale= advanced? 1.0: this.correctionsParameters.JPEG_scale;
		  boolean toRGB=     advanced? true: this.correctionsParameters.toRGB; 
		  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?

		  // may use this.StartTime to report intermediate steps execution times
		  String name=(String) imp_src.getProperty("name");
		  //		int channel= Integer.parseInt((String) imp_src.getProperty("channel"));
		  int channel= (Integer) imp_src.getProperty("channel");
		  String path= (String) imp_src.getProperty("path");
		  		  
		  String title=name+"-"+String.format("%02d", channel);
		  ImagePlus result=imp_src;
		  if (debugLevel>1) System.out.println("processing: "+path);
		  result.setTitle(title+"RAW");
		  if (!this.correctionsParameters.split){
			  eyesisCorrections.saveAndShow(result, this.correctionsParameters);
			  return result;
		  }
		  // Generate split parameters for DCT processing mode
		  // Split into Bayer components, oversample, increase canvas
		  
		  double [][] double_stack = eyesisCorrections.bayerToDoubleStack(
				  result, // source Bayer image, linearized, 32-bit (float))
				  null); // no margins, no oversample
		  
		  String titleFull=title+"-SPLIT";
		  if (debugLevel > -1){
			  double [] chn_avg = {0.0,0.0,0.0};
			  int width =  imp_src.getWidth();
			  int height = imp_src.getHeight();
			  for (int c = 0; c < 3; c++){
				  for (int i = 0; i<double_stack[c].length; i++){
					  chn_avg[c] += double_stack[c][i];
				  }
			  }
			  chn_avg[0] /= width*height/4;
			  chn_avg[1] /= width*height/4;
			  chn_avg[2] /= width*height/2;
			  System.out.println("Split channels averages: R="+chn_avg[0]+", G="+chn_avg[2]+", B="+chn_avg[1]); 
		  }
		  String [] rbg_titles = {"Red", "Blue", "Green"};
		  ImageStack stack;
		  if (!this.correctionsParameters.debayer) {
//			  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
//			  ImageStack 
			  stack = sdfa_instance.makeStack(double_stack, imp_src.getWidth(), imp_src.getHeight(), rbg_titles);
			  result= new ImagePlus(titleFull, stack);    			  
			  eyesisCorrections.saveAndShow(result, this.correctionsParameters);
			  return result;
		  }
		  // =================
		  if (debugLevel > 0) {
			  System.out.println("Showing image BEFORE_CLT_PROC");
			  sdfa_instance.showArrays(double_stack,  imp_src.getWidth(), imp_src.getHeight(), true, "BEFORE_CLT_PROC", rbg_titles);
		  }
		  if (this.correctionsParameters.deconvolve) { // process with DCT, otherwise use simple debayer
			  ImageDtt image_dtt = new ImageDtt();
/*			  
			  double [][][][][] clt_data = image_dtt.cltStack(
					  stack,
					  channel,
					  clt_parameters,
					  clt_parameters.ishift_x, //final int shiftX, // shift image horizontally (positive - right)
					  clt_parameters.ishift_y, //final int shiftY, // shift image vertically (positive - down)
					  threadsMax,
					  debugLevel,
					  updateStatus);
*/
			  for (int i =0 ; i < double_stack[0].length; i++){
//				  double_stack[0][i]*=2.0; // Scale red twice to compensate less pixels than green
//				  double_stack[1][i]*=2.0; // Scale blue twice to compensate less pixels than green
				  double_stack[2][i]*=0.5; // Scale blue twice to compensate less pixels than green
			  }
			  double [][][][][] clt_data = image_dtt.clt_aberrations( 
					  double_stack,                 // final double [][]       imade_data,
					  imp_src.getWidth(),           //	final int               width,
					  clt_kernels[channel],         // final double [][][][][] clt_kernels, // [color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
					  clt_parameters.kernel_step,
					  clt_parameters.transform_size,
					  clt_parameters.clt_window,
					  clt_parameters.shift_x,       // final int               shiftX, // shift image horizontally (positive - right) - just for testing
					  clt_parameters.shift_y,       // final int               shiftY, // shift image vertically (positive - down)
					  clt_parameters.tileX,         // final int               debug_tileX,
					  clt_parameters.tileY,         // final int               debug_tileY,
					  (clt_parameters.dbg_mode & 64) != 0, // no fract shift
					  (clt_parameters.dbg_mode & 128) != 0, // no convolve
					  (clt_parameters.dbg_mode & 256) != 0, // transpose convolve
					  threadsMax,
					  debugLevel);
//					  updateStatus);

			  
			  System.out.println("clt_data.length="+clt_data.length+" clt_data[0].length="+clt_data[0].length
					  +" clt_data[0][0].length="+clt_data[0][0].length+" clt_data[0][0][0].length="+
					  clt_data[0][0][0].length+" clt_data[0][0][0][0].length="+clt_data[0][0][0][0].length);
			  /*
			  if (dct_parameters.color_DCT){ // convert RBG -> YPrPb
				  dct_data = image_dtt.dct_color_convert(
						  dct_data,
						  colorProcParameters.kr,
						  colorProcParameters.kb,
						  dct_parameters.sigma_rb,        // blur of channels 0,1 (r,b) in addition to 2 (g)
						  dct_parameters.sigma_y,         // blur of Y from G
						  dct_parameters.sigma_color,     // blur of Pr, Pb in addition to Y
						  threadsMax,
						  debugLevel);
			  } else { // just LPF RGB
			  */
				  if (clt_parameters.corr_sigma > 0){ // no filter at all
					  for (int chn = 0; chn < clt_data.length; chn++) {
						  image_dtt.clt_lpf(
								  clt_parameters.corr_sigma,
								  clt_data[chn],
								  clt_parameters.transform_size,
								  threadsMax,
								  debugLevel);
					  }
				  }
/*
			  }
*/
			  int tilesY = imp_src.getHeight()/clt_parameters.transform_size;
			  int tilesX = imp_src.getWidth()/clt_parameters.transform_size;
			  if (debugLevel > 0){
				  System.out.println("--tilesX="+tilesX);
				  System.out.println("--tilesY="+tilesY);
			  }
			  if (debugLevel > 1){
			        double [][] clt = new double [clt_data.length*4][];
			        for (int chn = 0; chn < clt_data.length; chn++) {
			        	double [][] clt_set = image_dtt.clt_dbg(
			        			clt_data [chn],
								  threadsMax,
								  debugLevel);
			        	for (int ii = 0; ii < clt_set.length; ii++) clt[chn*4+ii] = clt_set[ii];
			        }

			        if (debugLevel > 0){
			        	sdfa_instance.showArrays(clt,
			        			tilesX*clt_parameters.transform_size,
			        			tilesY*clt_parameters.transform_size,
			        			true,
			        			result.getTitle()+"-CLT");  
			        }
			  }
			  double [][] iclt_data = new double [clt_data.length][];
			  for (int chn=0; chn<clt_data.length;chn++){
				  iclt_data[chn] = image_dtt.iclt_2d(
						  clt_data[chn],                  // scanline representation of dcd data, organized as dct_size x dct_size tiles  
						  clt_parameters.transform_size,  // final int
						  clt_parameters.clt_window,      // window_type
						  15,                             // clt_parameters.iclt_mask,       //which of 4 to transform back
						  0,                              // clt_parameters.dbg_mode,        //which of 4 to transform back
						  threadsMax,
						  debugLevel);
				  
			  }
			  
//					  if (debugLevel > -1) System.out.println("Applyed LPF, sigma = "+dct_parameters.dbg_sigma);
					  if (debugLevel > 0) sdfa_instance.showArrays(
							  iclt_data,
							  (tilesX + 1) * clt_parameters.transform_size,
							  (tilesY + 1) * clt_parameters.transform_size,
							  true,
							  result.getTitle()+"-rbg_sigma");
				/*	  
				  }
			  }
			 */
			  if (debugLevel > 0) sdfa_instance.showArrays(iclt_data,
					  (tilesX + 0) * clt_parameters.transform_size,
					  (tilesY + 0) * clt_parameters.transform_size,
					  true,
					  result.getTitle()+"-ICLT-RGB");

			  // convert to ImageStack of 3 slices
			  String [] sliceNames = {"red", "blue", "green"};
			  stack = sdfa_instance.makeStack(
					  iclt_data,
					  (tilesX + 0) * clt_parameters.transform_size,
					  (tilesY + 0) * clt_parameters.transform_size,
					  sliceNames); // or use null to get chn-nn slice names


		  } else { // if (this.correctionsParameters.deconvolve) - here use a simple debayer
			  System.out.println("Bypassing CLT-based aberration correction");
			  stack = sdfa_instance.makeStack(double_stack, imp_src.getWidth(), imp_src.getHeight(), rbg_titles);
			  debayer_rbg(stack, 0.25); // simple standard 3x3 kernel debayer
		  }
		  if (debugLevel > -1){
			  double [] chn_avg = {0.0,0.0,0.0};
			  float [] pixels;
			  int width =  stack.getWidth();
			  int height = stack.getHeight();
			  
			  for (int c = 0; c <3; c++){
				  pixels = (float[]) stack.getPixels(c+1);
				  for (int i = 0; i<pixels.length; i++){
					  chn_avg[c] += pixels[i];
				  }
			  }
			  chn_avg[0] /= width*height;
			  chn_avg[1] /= width*height;
			  chn_avg[2] /= width*height;
			  System.out.println("Processed channels averages: R="+chn_avg[0]+", G="+chn_avg[2]+", B="+chn_avg[1]); 
		  }
		  
		  if (!this.correctionsParameters.colorProc){
			  result= new ImagePlus(titleFull, stack);    			  
			  eyesisCorrections.saveAndShow(
					  result,
					  this.correctionsParameters);
			  return result;
		  }
		  if (debugLevel > 1) System.out.println("before colors.1");
		  //Processing colors - changing stack sequence to r-g-b (was r-b-g)
		  if (!eyesisCorrections.fixSliceSequence(
				  stack,
				  debugLevel)){
			  if (debugLevel > -1) System.out.println("fixSliceSequence() returned false");
			  return null;
		  }
		  if (debugLevel > 1) System.out.println("before colors.2");
		  if (debugLevel > 1){
			  ImagePlus imp_dbg=new ImagePlus(imp_src.getTitle()+"-"+channel+"-preColors",stack);
			  eyesisCorrections.saveAndShow(
					  imp_dbg,
					  this.correctionsParameters);
		  }
		  if (debugLevel > 1) System.out.println("before colors.3, scaleExposure="+scaleExposure+" scale = "+(255.0/eyesisCorrections.psfSubpixelShouldBe4/eyesisCorrections.psfSubpixelShouldBe4/scaleExposure));
		  CorrectionColorProc correctionColorProc=new CorrectionColorProc(eyesisCorrections.stopRequested);
		  double [][] yPrPb=new double [3][];
		  //			if (dct_parameters.color_DCT){
		  // need to get YPbPr - not RGB here				
		  //			} else {
		  correctionColorProc.processColorsWeights(stack, // just gamma convert? TODO: Cleanup? Convert directly form the linear YPrPb
				  //					  255.0/this.psfSubpixelShouldBe4/this.psfSubpixelShouldBe4, //  double scale,     // initial maximal pixel value (16))
				  //					  255.0/eyesisCorrections.psfSubpixelShouldBe4/eyesisCorrections.psfSubpixelShouldBe4/scaleExposure, //  double scale,     // initial maximal pixel value (16))
				  //					  255.0/2/2/scaleExposure, //  double scale,     // initial maximal pixel value (16))
				  255.0/scaleExposure, //  double scale,     // initial maximal pixel value (16))
				  colorProcParameters,
				  channelGainParameters,
				  channel,
				  null, //correctionDenoise.getDenoiseMask(),
				  this.correctionsParameters.blueProc,
				  debugLevel);
		  if (debugLevel > 1) System.out.println("Processed colors to YPbPr, total number of slices="+stack.getSize());
		  if (debugLevel > 1) {
			  ImagePlus imp_dbg=new ImagePlus("procColors",stack);
			  eyesisCorrections.saveAndShow(
					  imp_dbg,
					  this.correctionsParameters);
		  }
		  float [] fpixels;
		  int [] slices_YPrPb = {8,6,7};
		  yPrPb=new double [3][];
		  for (int n = 0; n < slices_YPrPb.length; n++){
			  fpixels = (float[]) stack.getPixels(slices_YPrPb[n]);
			  yPrPb[n] = new double [fpixels.length];
			  for (int i = 0; i < fpixels.length; i++) yPrPb[n][i] = fpixels[i];
		  }

		  if (toRGB) {
			  System.out.println("correctionColorProc.YPrPbToRGB");
			  stack =  YPrPbToRGB(yPrPb,
					  colorProcParameters.kr,        // 0.299;
					  colorProcParameters.kb,        // 0.114;
					  stack.getWidth());

			  title=titleFull; // including "-DECONV" or "-COMBO"
			  titleFull=title+"-RGB-float";
			  //Trim stack to just first 3 slices
			  if (debugLevel > 1){ // 2){
				  ImagePlus imp_dbg=new ImagePlus("YPrPbToRGB",stack);
				  eyesisCorrections.saveAndShow(
						  imp_dbg,
						  this.correctionsParameters);
			  }
			  while (stack.getSize() > 3) stack.deleteLastSlice();
			  if (debugLevel > 1) System.out.println("Trimming color stack");
		  } else {
			  title=titleFull; // including "-DECONV" or "-COMBO"
			  titleFull=title+"-YPrPb"; // including "-DECONV" or "-COMBO"
			  if (debugLevel > 1) System.out.println("Using full stack, including YPbPr");
		  }

		  result= new ImagePlus(titleFull, stack);    			  
		  // Crop image to match original one (scaled to oversampling)
/*		  
		  if (crop){ // always crop if equirectangular
			  if (debugLevel > 1) System.out.println("cropping");
			  stack = eyesisCorrections.cropStack32(stack,splitParameters);
			  if (debugLevel > 2) { // 2){
				  ImagePlus imp_dbg=new ImagePlus("cropped",stack);
				  eyesisCorrections.saveAndShow(
						  imp_dbg,
						  this.correctionsParameters);
			  }
		  }
*/		  
		  // rotate the result			  
		  if (rotate){ // never rotate for equirectangular
			  stack=eyesisCorrections.rotateStack32CW(stack);
		  }
		  if (!toRGB && !this.correctionsParameters.jpeg){ // toRGB set for equirectangular
			  if (debugLevel > 1) System.out.println("!toRGB && !this.correctionsParameters.jpeg");
			  eyesisCorrections.saveAndShow(result, this.correctionsParameters);
			  return result;
		  } else { // that's not the end result, save if required
			  if (debugLevel > 1) System.out.println("!toRGB && !this.correctionsParameters.jpeg - else");
			  eyesisCorrections.saveAndShow(result,
					  eyesisCorrections.correctionsParameters,
					  eyesisCorrections.correctionsParameters.save32,
					  false,
					  eyesisCorrections.correctionsParameters.JPEG_quality); // save, no show
		  }
		  // convert to RGB48 (16 bits per color component)
		  ImagePlus imp_RGB;
		  stack=eyesisCorrections.convertRGB32toRGB16Stack(
				  stack,
				  rgbParameters); 

		  titleFull=title+"-RGB48";
		  result= new ImagePlus(titleFull, stack);
		  //			  ImagePlus imp_RGB24;
		  result.updateAndDraw();
		  if (debugLevel > 1) System.out.println("result.updateAndDraw(), "+titleFull+"-RGB48");

		  CompositeImage compositeImage=eyesisCorrections.convertToComposite(result);

		  if (!this.correctionsParameters.jpeg && !advanced){ // RGB48 was the end result
			  if (debugLevel > 1) System.out.println("if (!this.correctionsParameters.jpeg && !advanced)");
			  eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters);
			  return result;
		  } else { // that's not the end result, save if required
			  if (debugLevel > 1) System.out.println("if (!this.correctionsParameters.jpeg && !advanced) - else");
			  eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters, this.correctionsParameters.save16, false); // save, no show
			  //				  eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters, this.correctionsParameters.save16, true); // save, no show
		  }

		  imp_RGB=eyesisCorrections.convertRGB48toRGB24(
				  stack,
				  title+"-RGB24",
				  0, 65536, // r range 0->0, 65536->256
				  0, 65536, // g range
				  0, 65536,// b range
				  0, 65536);// alpha range
		  if (JPEG_scale!=1.0){
			  ImageProcessor ip=imp_RGB.getProcessor();
			  ip.setInterpolationMethod(ImageProcessor.BICUBIC);
			  ip=ip.resize((int)(ip.getWidth()*JPEG_scale),(int) (ip.getHeight()*JPEG_scale));
			  imp_RGB= new ImagePlus(imp_RGB.getTitle(),ip);
			  imp_RGB.updateAndDraw();
		  }
		  eyesisCorrections.saveAndShow(imp_RGB, this.correctionsParameters);

		  return result;
	  }
	  
	  public void processCLTQuads(
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
			  int          convolveFFTSize, // 128 - fft size, kernel size should be size/2
			  final int          threadsMax,  // maximal number of threads to launch                         
			  final boolean    updateStatus,
			  final int        debugLevel)
	  {
		  this.startTime=System.nanoTime();
		  String [] sourceFiles=correctionsParameters.getSourcePaths();
		  boolean [] enabledFiles=new boolean[sourceFiles.length];
		  for (int i=0;i<enabledFiles.length;i++) enabledFiles[i]=false;
		  int numFilesToProcess=0;
		  int numImagesToProcess=0;
		  for (int nFile=0;nFile<enabledFiles.length;nFile++){
			  if ((sourceFiles[nFile]!=null) && (sourceFiles[nFile].length()>1)) {
				  int [] channels={correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile])};
				  if (correctionsParameters.isJP4()){
					  int subCamera= channels[0]- correctionsParameters.firstSubCamera; // to match those in the sensor files
					  // removeUnusedSensorData should be off!?
					  channels=this.eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (eyesisCorrections.isChannelEnabled(channels[i])){
						  if (!enabledFiles[nFile]) numFilesToProcess++;
						  enabledFiles[nFile]=true;
						  numImagesToProcess++;
					  }
				  }
			  }
		  }
		  if (numFilesToProcess==0){
			  System.out.println("No files to process (of "+sourceFiles.length+")");
			  return;
		  } else {
			  if (debugLevel>0) System.out.println(numFilesToProcess+ " files to process (of "+sourceFiles.length+"), "+numImagesToProcess+" images to process");
		  }
		  double [] referenceExposures=eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		  int [][] fileIndices=new int [numImagesToProcess][2]; // file index, channel number
		  int index=0;
		  for (int nFile=0;nFile<enabledFiles.length;nFile++){
			  if ((sourceFiles[nFile]!=null) && (sourceFiles[nFile].length()>1)) {
				  int [] channels={correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile])};
				  if (correctionsParameters.isJP4()){
					  int subCamera= channels[0]- correctionsParameters.firstSubCamera; // to match those in the sensor files
					  channels=eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (eyesisCorrections.isChannelEnabled(channels[i])){
						  fileIndices[index  ][0]=nFile;
						  fileIndices[index++][1]=channels[i];
					  }
				  }
			  }
		  }
		  ArrayList<String> setNames = new ArrayList<String>();
		  ArrayList<ArrayList<Integer>> setFiles = new ArrayList<ArrayList<Integer>>();

		  for (int iImage=0;iImage<fileIndices.length;iImage++){
			  int nFile=fileIndices[iImage][0];
			  String setName = correctionsParameters.getNameFromSourceTiff(sourceFiles[nFile]);
			  if (!setNames.contains(setName)) {
				  setNames.add(setName);
				  setFiles.add(new ArrayList<Integer>());
			  }
			  setFiles.get(setNames.indexOf(setName)).add(new Integer(nFile));
		  }		  
		  for (int nSet = 0; nSet < setNames.size(); nSet++){
			  int maxChn = 0;
			  for (int i = 0; i < setFiles.get(nSet).size(); i++){
				  int chn = fileIndices[setFiles.get(nSet).get(i)][1];
				  if (chn > maxChn) maxChn = chn;
			  }
			  int [] channelFiles = new int[maxChn+1];
			  for (int i =0; i < channelFiles.length; i++) channelFiles[i] = -1;
			  for (int i = 0; i < setFiles.get(nSet).size(); i++){
				  channelFiles[fileIndices[setFiles.get(nSet).get(i)][1]] = setFiles.get(nSet).get(i);
			  }
			  
			  ImagePlus [] imp_srcs = new ImagePlus[channelFiles.length];
			  double [] scaleExposures = new double[channelFiles.length];
			  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
				  int nFile=channelFiles[srcChannel];
				  imp_srcs[srcChannel]=null;
				  if (nFile >=0){
					  if (correctionsParameters.isJP4()){
						  int subchannel=eyesisCorrections.pixelMapping.getSubChannel(srcChannel);
						  if (this.correctionsParameters.swapSubchannels01) {
							  switch (subchannel){
							  case 0: subchannel=1; break;
							  case 1: subchannel=0; break;
							  }
						  }
						  if (debugLevel>0) System.out.println("Processing set " + setNames.get(nSet)+" channel "+srcChannel+" - subchannel "+subchannel+" of "+sourceFiles[nFile]);
						  ImagePlus imp_composite=eyesisCorrections.JP4_INSTANCE.open(
								  "", // path,
								  sourceFiles[nFile],
								  "",  //arg - not used in JP46 reader
								  true, // un-apply camera color gains
								  null, // new window
								  false); // do not show
						  imp_srcs[srcChannel]=eyesisCorrections.JP4_INSTANCE.demuxImage(imp_composite, subchannel);
						  if (imp_srcs[srcChannel] == null) imp_srcs[srcChannel] = imp_composite; // not a composite image
						  // do we need to add any properties?				  
					  } else { 
						  imp_srcs[srcChannel]=new ImagePlus(sourceFiles[nFile]);
						  //					  (new JP46_Reader_camera(false)).decodeProperiesFromInfo(imp_src); // decode existent properties from info
						  eyesisCorrections.JP4_INSTANCE.decodeProperiesFromInfo(imp_srcs[srcChannel]); // decode existent properties from info
						  if (debugLevel>0) System.out.println("Processing "+sourceFiles[nFile]);
					  }
					  scaleExposures[srcChannel] = 1.0;
					  if (!Double.isNaN(referenceExposures[nFile]) && (imp_srcs[srcChannel].getProperty("EXPOSURE")!=null)){
						  scaleExposures[srcChannel] = referenceExposures[nFile]/Double.parseDouble((String) imp_srcs[srcChannel].getProperty("EXPOSURE"));
						  if (debugLevel>0) System.out.println("Will scale intensity (to compensate for exposure) by  "+scaleExposures[srcChannel]);
					  }
					  imp_srcs[srcChannel].setProperty("name",    correctionsParameters.getNameFromSourceTiff(sourceFiles[nFile]));
					  imp_srcs[srcChannel].setProperty("channel", srcChannel); // it may already have channel
					  imp_srcs[srcChannel].setProperty("path",    sourceFiles[nFile]); // it may already have channel

					  if (this.correctionsParameters.pixelDefects && (eyesisCorrections.defectsXY!=null)&& (eyesisCorrections.defectsXY[srcChannel]!=null)){
						  // apply pixel correction
						  int numApplied=	eyesisCorrections.correctDefects(
								  imp_srcs[srcChannel],
								  srcChannel,
								  debugLevel);
						  if ((debugLevel>0) && (numApplied>0)) { // reduce verbosity after verified defect correction works
							  System.out.println("Corrected "+numApplied+" pixels in "+sourceFiles[nFile]);
						  }
					  }

					  if (this.correctionsParameters.vignetting){
						  if ((eyesisCorrections.channelVignettingCorrection==null) || (srcChannel<0) || (srcChannel>=eyesisCorrections.channelVignettingCorrection.length) || (eyesisCorrections.channelVignettingCorrection[srcChannel]==null)){
							  System.out.println("No vignetting data for channel "+srcChannel);
							  return;
						  }
						  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
						  if (pixels.length!=eyesisCorrections.channelVignettingCorrection[srcChannel].length){
							  System.out.println("Vignetting data for channel "+srcChannel+" has "+eyesisCorrections.channelVignettingCorrection[srcChannel].length+" pixels, image "+sourceFiles[nFile]+" has "+pixels.length);
							  return;
						  }
						  // TODO: Move to do it once:
						  double min_non_zero = 0.0;
						  for (int i=0;i<pixels.length;i++){
							  double d = eyesisCorrections.channelVignettingCorrection[srcChannel][i];
							  if ((d > 0.0) && ((min_non_zero == 0) || (min_non_zero > d))){
								  min_non_zero = d;
							  }
						  }
						  double max_vign_corr = clt_parameters.vignetting_range*min_non_zero;

						  System.out.println("Vignetting data: channel="+srcChannel+", min = "+min_non_zero);
						  for (int i=0;i<pixels.length;i++){
							  double d = eyesisCorrections.channelVignettingCorrection[srcChannel][i];
							  if (d > max_vign_corr) d = max_vign_corr;
							  pixels[i]*=d;
						  }
						  // Scale here, combine with vignetting later?
						  int width =  imp_srcs[srcChannel].getWidth();
						  int height = imp_srcs[srcChannel].getHeight();
						  for (int y = 0; y < height-1; y+=2){
							  for (int x = 0; x < width-1; x+=2){
								  pixels[y*width+x        ] *= clt_parameters.scale_g;
								  pixels[y*width+x+width+1] *= clt_parameters.scale_g;
								  pixels[y*width+x      +1] *= clt_parameters.scale_r;
								  pixels[y*width+x+width  ] *= clt_parameters.scale_b;
							  }
						  }

					  } else { // assuming GR/BG pattern
						  System.out.println("Applying fixed color gain correction parameters: Gr="+
								  clt_parameters.novignetting_r+", Gg="+clt_parameters.novignetting_g+", Gb="+clt_parameters.novignetting_b);
						  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
						  int width =  imp_srcs[srcChannel].getWidth();
						  int height = imp_srcs[srcChannel].getHeight();
						  double kr = clt_parameters.scale_r/clt_parameters.novignetting_r;
						  double kg = clt_parameters.scale_g/clt_parameters.novignetting_g;
						  double kb = clt_parameters.scale_b/clt_parameters.novignetting_b;
						  for (int y = 0; y < height-1; y+=2){
							  for (int x = 0; x < width-1; x+=2){
								  pixels[y*width+x        ] *= kg;
								  pixels[y*width+x+width+1] *= kg;
								  pixels[y*width+x      +1] *= kr;
								  pixels[y*width+x+width  ] *= kb;
							  }
						  }
					  }
				  }
			  }
			  // once per quad here
			  // may need to equalize gains between channels
			  // may need to equalize gains between channels
			  if (clt_parameters.gain_equalize || clt_parameters.colors_equalize){
				  channelGainsEqualize(
						  clt_parameters.gain_equalize,
						  clt_parameters.colors_equalize,
						  channelFiles,
						  imp_srcs,
						  setNames.get(nSet), // just for debug messeges == setNames.get(nSet)
						  debugLevel);
			  }
			  // once per quad here
			  processCLTQuad( // returns ImagePlus, but it already should be saved/shown
					  imp_srcs, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
					  clt_parameters,
					  debayerParameters,
					  nonlinParameters,
					  colorProcParameters,
					  channelGainParameters,
					  rgbParameters,
					  convolveFFTSize, // 128 - fft size, kernel size should be size/2
					  scaleExposures,
					  threadsMax,  // maximal number of threads to launch                         
					  updateStatus,
					  debugLevel);
			  Runtime.getRuntime().gc();
			  if (debugLevel >-1) System.out.println("Processing set "+(nSet+1)+" (of "+fileIndices.length+") finished at "+
					  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
			  if (eyesisCorrections.stopRequested.get()>0) {
				  System.out.println("User requested stop");
				  return;
			  }
		  }
		  System.out.println("Processing "+fileIndices.length+" files finished at "+
				  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
	  }		
	  
	  public ImagePlus [] processCLTQuad(
			  ImagePlus [] imp_quad, // should have properties "name"(base for saving results), "channel","path"
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  int          convolveFFTSize, // 128 - fft size, kernel size should be size/2
			  double []	       scaleExposures, // probably not needed here
			  final int        threadsMax,  // maximal number of threads to launch                         
			  final boolean    updateStatus,
			  final int        debugLevel){
		  boolean advanced=this.correctionsParameters.zcorrect || this.correctionsParameters.equirectangular;
//		  boolean crop=      advanced? true: this.correctionsParameters.crop; 
		  boolean rotate=    advanced? false: this.correctionsParameters.rotate; 
		  double JPEG_scale= advanced? 1.0: this.correctionsParameters.JPEG_scale;
		  boolean toRGB=     advanced? true: this.correctionsParameters.toRGB; 
		  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?

		  // may use this.StartTime to report intermediate steps execution times
		  String name=(String) imp_quad[0].getProperty("name");
		  //		int channel= Integer.parseInt((String) imp_src.getProperty("channel"));
		  int channel= (Integer) imp_quad[0].getProperty("channel");
		  String path= (String) imp_quad[0].getProperty("path");
		  		  
		  ImagePlus [] results = new ImagePlus[imp_quad.length];
		  for (int i = 0; i < results.length; i++) {
			  results[i] = imp_quad[i];
			  results[i].setTitle(results[i].getTitle()+"RAW");			  
		  }
		  if (debugLevel>1) System.out.println("processing: "+path);
		  double [][][] double_stacks = new double [imp_quad.length][][];
		  for (int i = 0; i < double_stacks.length; i++){
			  double_stacks[i] = eyesisCorrections.bayerToDoubleStack(
					  imp_quad[i], // source Bayer image, linearized, 32-bit (float))
					  null); // no margins, no oversample
		  }
		  
//		  String [] rbg_titles = {"Red", "Blue", "Green"};
		  ImageStack stack;
		  // =================
		  ImageDtt image_dtt = new ImageDtt();
		  for (int i = 0; i < double_stacks.length; i++){
			  for (int j =0 ; j < double_stacks[i][0].length; j++){
				  double_stacks[i][2][j]*=0.5; // Scale green 0.5 to compensate more pixels than R,B
			  }
		  }
		  double [][][][][][] clt_data = image_dtt.clt_aberrations_quad( 
				  clt_parameters.disparity,     // final double            disparity,
				  double_stacks,                // final double [][][]      imade_data, // first index - number of image in a quad
				  imp_quad[0].getWidth(),       //	final int               width,
				  geometryCorrection,           // final GeometryCorrection  geometryCorrection,
				  clt_kernels,                  // final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				  clt_parameters.kernel_step,
				  clt_parameters.transform_size,
				  clt_parameters.clt_window,
				  clt_parameters.shift_x,       // final int               shiftX, // shift image horizontally (positive - right) - just for testing
				  clt_parameters.shift_y,       // final int               shiftY, // shift image vertically (positive - down)
				  clt_parameters.tileX,         // final int               debug_tileX,
				  clt_parameters.tileY,         // final int               debug_tileY,
				  (clt_parameters.dbg_mode & 64) != 0, // no fract shift
				  (clt_parameters.dbg_mode & 128) != 0, // no convolve
				  (clt_parameters.dbg_mode & 256) != 0, // transpose convolve
				  threadsMax,
				  debugLevel);

		  System.out.println("clt_data.length="+clt_data.length+" clt_data[0].length="+clt_data[0].length
				  +" clt_data[0][0].length="+clt_data[0][0].length+" clt_data[0][0][0].length="+
				  clt_data[0][0][0].length+" clt_data[0][0][0][0].length="+clt_data[0][0][0][0].length+
				  " clt_data[0][0][0][0][0].length="+clt_data[0][0][0][0][0].length);
		  
		  for (int iQuad = 0; iQuad <clt_data.length; iQuad++){
			  
			  String title=name+"-"+String.format("%02d", iQuad);
			  String titleFull=title+"-SPLIT";
			  
			  if (clt_parameters.corr_sigma > 0){ // no filter at all
				  for (int chn = 0; chn < clt_data[iQuad].length; chn++) {
					  image_dtt.clt_lpf(
							  clt_parameters.corr_sigma,
							  clt_data[iQuad][chn],
							  clt_parameters.transform_size,
							  threadsMax,
							  debugLevel);
				  }
			  }

			  int tilesY = imp_quad[iQuad].getHeight()/clt_parameters.transform_size;
			  int tilesX = imp_quad[iQuad].getWidth()/clt_parameters.transform_size;
			  if (debugLevel > 0){
				  System.out.println("--tp.tilesX="+tilesX);
				  System.out.println("--tp.tilesY="+tilesY);
			  }
			  if (debugLevel > 1){
				  double [][] clt = new double [clt_data[iQuad].length*4][];
				  for (int chn = 0; chn < clt_data[iQuad].length; chn++) {
					  double [][] clt_set = image_dtt.clt_dbg(
							  clt_data [iQuad][chn],
							  threadsMax,
							  debugLevel);
					  for (int ii = 0; ii < clt_set.length; ii++) clt[chn*4+ii] = clt_set[ii];
				  }

				  if (debugLevel > 0){
					  sdfa_instance.showArrays(clt,
							  tilesX*clt_parameters.transform_size,
							  tilesY*clt_parameters.transform_size,
							  true,
							  results[iQuad].getTitle()+"-CLT");  
				  }
			  }
			  double [][] iclt_data = new double [clt_data[iQuad].length][];
			  for (int chn=0; chn<iclt_data.length;chn++){
				  iclt_data[chn] = image_dtt.iclt_2d(
						  clt_data[iQuad][chn],           // scanline representation of dcd data, organized as dct_size x dct_size tiles  
						  clt_parameters.transform_size,  // final int
						  clt_parameters.clt_window,      // window_type
						  15,                             // clt_parameters.iclt_mask,       //which of 4 to transform back
						  0,                              // clt_parameters.dbg_mode,        //which of 4 to transform back
						  threadsMax,
						  debugLevel);

			  }
			  if (debugLevel > 0) sdfa_instance.showArrays(
					  iclt_data,
					  (tilesX + 0) * clt_parameters.transform_size,
					  (tilesY + 0) * clt_parameters.transform_size,
					  true,
					  results[iQuad].getTitle()+"-rbg_sigma");
			  if (debugLevel > 0) sdfa_instance.showArrays(iclt_data,
					  (tilesX + 0) * clt_parameters.transform_size,
					  (tilesY + 0) * clt_parameters.transform_size,
					  true,
					  results[iQuad].getTitle()+"-ICLT-RGB");

			  // convert to ImageStack of 3 slices
			  String [] sliceNames = {"red", "blue", "green"};
			  stack = sdfa_instance.makeStack(
					  iclt_data,
					  (tilesX + 0) * clt_parameters.transform_size,
					  (tilesY + 0) * clt_parameters.transform_size,
					  sliceNames); // or use null to get chn-nn slice names

			  if (debugLevel > -1){
				  double [] chn_avg = {0.0,0.0,0.0};
				  float [] pixels;
				  int width =  stack.getWidth();
				  int height = stack.getHeight();

				  for (int c = 0; c <3; c++){
					  pixels = (float[]) stack.getPixels(c+1);
					  for (int i = 0; i<pixels.length; i++){
						  chn_avg[c] += pixels[i];
					  }
				  }
				  chn_avg[0] /= width*height;
				  chn_avg[1] /= width*height;
				  chn_avg[2] /= width*height;
				  System.out.println("Processed channels averages: R="+chn_avg[0]+", G="+chn_avg[2]+", B="+chn_avg[1]); 
			  }

			  if (!this.correctionsParameters.colorProc){
				  results[iQuad]= new ImagePlus(titleFull, stack);    			  
				  eyesisCorrections.saveAndShow(
						  results[iQuad],
						  this.correctionsParameters);
				  continue; // return results;
			  }
			  if (debugLevel > 1) System.out.println("before colors.1");
			  //Processing colors - changing stack sequence to r-g-b (was r-b-g)
			  if (!eyesisCorrections.fixSliceSequence(
					  stack,
					  debugLevel)){
				  if (debugLevel > -1) System.out.println("fixSliceSequence() returned false");
				  return null;
			  }
			  if (debugLevel > 1) System.out.println("before colors.2");
			  if (debugLevel > 1){
				  ImagePlus imp_dbg=new ImagePlus(imp_quad[iQuad].getTitle()+"-"+channel+"-preColors",stack);
				  eyesisCorrections.saveAndShow(
						  imp_dbg,
						  this.correctionsParameters);
			  }
			  if (debugLevel > 1) System.out.println("before colors.3, scaleExposure="+scaleExposures[iQuad]+" scale = "+(255.0/eyesisCorrections.psfSubpixelShouldBe4/eyesisCorrections.psfSubpixelShouldBe4/scaleExposures[iQuad]));
			  CorrectionColorProc correctionColorProc=new CorrectionColorProc(eyesisCorrections.stopRequested);
			  double [][] yPrPb=new double [3][];
			  //			if (dct_parameters.color_DCT){
			  // need to get YPbPr - not RGB here				
			  //			} else {
			  correctionColorProc.processColorsWeights(stack, // just gamma convert? TODO: Cleanup? Convert directly form the linear YPrPb
					  //					  255.0/this.psfSubpixelShouldBe4/this.psfSubpixelShouldBe4, //  double scale,     // initial maximal pixel value (16))
					  //					  255.0/eyesisCorrections.psfSubpixelShouldBe4/eyesisCorrections.psfSubpixelShouldBe4/scaleExposure, //  double scale,     // initial maximal pixel value (16))
					  //					  255.0/2/2/scaleExposure, //  double scale,     // initial maximal pixel value (16))
					  255.0/scaleExposures[iQuad], //  double scale,     // initial maximal pixel value (16))
					  colorProcParameters,
					  channelGainParameters,
					  channel,
					  null, //correctionDenoise.getDenoiseMask(),
					  this.correctionsParameters.blueProc,
					  debugLevel);
			  if (debugLevel > 1) System.out.println("Processed colors to YPbPr, total number of slices="+stack.getSize());
			  if (debugLevel > 1) {
				  ImagePlus imp_dbg=new ImagePlus("procColors",stack);
				  eyesisCorrections.saveAndShow(
						  imp_dbg,
						  this.correctionsParameters);
			  }
			  float [] fpixels;
			  int [] slices_YPrPb = {8,6,7};
			  yPrPb=new double [3][];
			  for (int n = 0; n < slices_YPrPb.length; n++){
				  fpixels = (float[]) stack.getPixels(slices_YPrPb[n]);
				  yPrPb[n] = new double [fpixels.length];
				  for (int i = 0; i < fpixels.length; i++) yPrPb[n][i] = fpixels[i];
			  }

			  if (toRGB) {
				  System.out.println("correctionColorProc.YPrPbToRGB");
				  stack =  YPrPbToRGB(yPrPb,
						  colorProcParameters.kr,        // 0.299;
						  colorProcParameters.kb,        // 0.114;
						  stack.getWidth());

				  title=titleFull; // including "-DECONV" or "-COMBO"
				  titleFull=title+"-RGB-float";
				  //Trim stack to just first 3 slices
				  if (debugLevel > 1){ // 2){
					  ImagePlus imp_dbg=new ImagePlus("YPrPbToRGB",stack);
					  eyesisCorrections.saveAndShow(
							  imp_dbg,
							  this.correctionsParameters);
				  }
				  while (stack.getSize() > 3) stack.deleteLastSlice();
				  if (debugLevel > 1) System.out.println("Trimming color stack");
			  } else {
				  title=titleFull; // including "-DECONV" or "-COMBO"
				  titleFull=title+"-YPrPb"; // including "-DECONV" or "-COMBO"
				  if (debugLevel > 1) System.out.println("Using full stack, including YPbPr");
			  }

			  results[iQuad]= new ImagePlus(titleFull, stack);    			  
			  // rotate the result			  
			  if (rotate){ // never rotate for equirectangular
				  stack=eyesisCorrections.rotateStack32CW(stack);
			  }
			  if (!toRGB && !this.correctionsParameters.jpeg){ // toRGB set for equirectangular
				  if (debugLevel > 1) System.out.println("!toRGB && !this.correctionsParameters.jpeg");
				  eyesisCorrections.saveAndShow(results[iQuad], this.correctionsParameters);
				  continue; // return result;
			  } else { // that's not the end result, save if required
				  if (debugLevel > 1) System.out.println("!toRGB && !this.correctionsParameters.jpeg - else");
				  eyesisCorrections.saveAndShow(results[iQuad],
						  eyesisCorrections.correctionsParameters,
						  eyesisCorrections.correctionsParameters.save32,
						  false,
						  eyesisCorrections.correctionsParameters.JPEG_quality); // save, no show
			  }
			  // convert to RGB48 (16 bits per color component)
			  ImagePlus imp_RGB;
			  stack=eyesisCorrections.convertRGB32toRGB16Stack(
					  stack,
					  rgbParameters); 

			  titleFull=title+"-RGB48";
			  results[iQuad]= new ImagePlus(titleFull, stack);    			  
			  
			  //			  ImagePlus imp_RGB24;
			  results[iQuad].updateAndDraw();
			  if (debugLevel > 1) System.out.println("result.updateAndDraw(), "+titleFull+"-RGB48");

			  CompositeImage compositeImage=eyesisCorrections.convertToComposite(results[iQuad]);

			  if (!this.correctionsParameters.jpeg && !advanced){ // RGB48 was the end result
				  if (debugLevel > 1) System.out.println("if (!this.correctionsParameters.jpeg && !advanced)");
				  eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters);
				  continue; // return result;
			  } else { // that's not the end result, save if required
				  if (debugLevel > 1) System.out.println("if (!this.correctionsParameters.jpeg && !advanced) - else");
				  eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters, this.correctionsParameters.save16, false); // save, no show
				  //				  eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters, this.correctionsParameters.save16, true); // save, no show
			  }

			  imp_RGB=eyesisCorrections.convertRGB48toRGB24(
					  stack,
					  title+"-RGB24",
					  0, 65536, // r range 0->0, 65536->256
					  0, 65536, // g range
					  0, 65536,// b range
					  0, 65536);// alpha range
			  if (JPEG_scale!=1.0){
				  ImageProcessor ip=imp_RGB.getProcessor();
				  ip.setInterpolationMethod(ImageProcessor.BICUBIC);
				  ip=ip.resize((int)(ip.getWidth()*JPEG_scale),(int) (ip.getHeight()*JPEG_scale));
				  imp_RGB= new ImagePlus(imp_RGB.getTitle(),ip);
				  imp_RGB.updateAndDraw();
			  }
			  eyesisCorrections.saveAndShow(imp_RGB, this.correctionsParameters);
		  }
		  return results;
	  }
	  	  
	  public void processCLTQuadCorrs(
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
			  int          convolveFFTSize, // 128 - fft size, kernel size should be size/2
			  final boolean    apply_corr, // calculate and apply additional fine geometry correction 
			  final int          threadsMax,  // maximal number of threads to launch                         
			  final boolean    updateStatus,
			  final int        debugLevel)
	  {
		  this.startTime=System.nanoTime();
		  String [] sourceFiles=correctionsParameters.getSourcePaths();
		  boolean [] enabledFiles=new boolean[sourceFiles.length];
		  for (int i=0;i<enabledFiles.length;i++) enabledFiles[i]=false;
		  int numFilesToProcess=0;
		  int numImagesToProcess=0;
		  for (int nFile=0;nFile<enabledFiles.length;nFile++){
			  if ((sourceFiles[nFile]!=null) && (sourceFiles[nFile].length()>1)) {
				  int [] channels={correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile])};
				  if (correctionsParameters.isJP4()){
					  int subCamera= channels[0]- correctionsParameters.firstSubCamera; // to match those in the sensor files
					  // removeUnusedSensorData should be off!?
					  channels=this.eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (eyesisCorrections.isChannelEnabled(channels[i])){
						  if (!enabledFiles[nFile]) numFilesToProcess++;
						  enabledFiles[nFile]=true;
						  numImagesToProcess++;
					  }
				  }
			  }
		  }
		  if (numFilesToProcess==0){
			  System.out.println("No files to process (of "+sourceFiles.length+")");
			  return;
		  } else {
			  if (debugLevel>0) System.out.println(numFilesToProcess+ " files to process (of "+sourceFiles.length+"), "+numImagesToProcess+" images to process");
		  }
		  double [] referenceExposures=eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		  int [][] fileIndices=new int [numImagesToProcess][2]; // file index, channel number
		  int index=0;
		  for (int nFile=0;nFile<enabledFiles.length;nFile++){
			  if ((sourceFiles[nFile]!=null) && (sourceFiles[nFile].length()>1)) {
				  int [] channels={correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile])};
				  if (correctionsParameters.isJP4()){
					  int subCamera= channels[0]- correctionsParameters.firstSubCamera; // to match those in the sensor files
					  channels=eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (eyesisCorrections.isChannelEnabled(channels[i])){
						  fileIndices[index  ][0]=nFile;
						  fileIndices[index++][1]=channels[i];
					  }
				  }
			  }
		  }
		  ArrayList<String> setNames = new ArrayList<String>();
		  ArrayList<ArrayList<Integer>> setFiles = new ArrayList<ArrayList<Integer>>();

		  for (int iImage=0;iImage<fileIndices.length;iImage++){
			  int nFile=fileIndices[iImage][0];
			  String setName = correctionsParameters.getNameFromSourceTiff(sourceFiles[nFile]);
			  if (!setNames.contains(setName)) {
				  setNames.add(setName);
				  setFiles.add(new ArrayList<Integer>());
			  }
			  setFiles.get(setNames.indexOf(setName)).add(new Integer(nFile));
		  }		  
		  for (int nSet = 0; nSet < setNames.size(); nSet++){
			  int maxChn = 0;
			  for (int i = 0; i < setFiles.get(nSet).size(); i++){
				  int chn = fileIndices[setFiles.get(nSet).get(i)][1];
				  if (chn > maxChn) maxChn = chn;
			  }
			  int [] channelFiles = new int[maxChn+1];
			  for (int i =0; i < channelFiles.length; i++) channelFiles[i] = -1;
			  for (int i = 0; i < setFiles.get(nSet).size(); i++){
				  channelFiles[fileIndices[setFiles.get(nSet).get(i)][1]] = setFiles.get(nSet).get(i);
			  }

			  ImagePlus [] imp_srcs = new ImagePlus[channelFiles.length];
			  double [] scaleExposures = new double[channelFiles.length]; //
			  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
				  int nFile=channelFiles[srcChannel];
				  imp_srcs[srcChannel]=null;
				  if (nFile >=0){
					  if (correctionsParameters.isJP4()){
						  int subchannel=eyesisCorrections.pixelMapping.getSubChannel(srcChannel);
						  if (this.correctionsParameters.swapSubchannels01) {
							  switch (subchannel){
							  case 0: subchannel=1; break;
							  case 1: subchannel=0; break;
							  }
						  }
						  if (debugLevel>0) System.out.println("Processing set " + setNames.get(nSet)+" channel "+srcChannel+" - subchannel "+subchannel+" of "+sourceFiles[nFile]);
						  ImagePlus imp_composite=eyesisCorrections.JP4_INSTANCE.open(
								  "", // path,
								  sourceFiles[nFile],
								  "",  //arg - not used in JP46 reader
								  true, // un-apply camera color gains
								  null, // new window
								  false); // do not show
						  imp_srcs[srcChannel]=eyesisCorrections.JP4_INSTANCE.demuxImage(imp_composite, subchannel);
						  if (imp_srcs[srcChannel] == null) imp_srcs[srcChannel] = imp_composite; // not a composite image
						  // do we need to add any properties?				  
					  } else { 
						  imp_srcs[srcChannel]=new ImagePlus(sourceFiles[nFile]);
						  //					  (new JP46_Reader_camera(false)).decodeProperiesFromInfo(imp_src); // decode existent properties from info
						  eyesisCorrections.JP4_INSTANCE.decodeProperiesFromInfo(imp_srcs[srcChannel]); // decode existent properties from info
						  if (debugLevel>0) System.out.println("Processing "+sourceFiles[nFile]);
					  }
					  scaleExposures[srcChannel] = 1.0;
					  if (!Double.isNaN(referenceExposures[nFile]) && (imp_srcs[srcChannel].getProperty("EXPOSURE")!=null)){
						  scaleExposures[srcChannel] = referenceExposures[nFile]/Double.parseDouble((String) imp_srcs[srcChannel].getProperty("EXPOSURE"));
						  if (debugLevel>0) System.out.println("Will scale intensity (to compensate for exposure) by  "+scaleExposures[srcChannel]);
					  }
					  imp_srcs[srcChannel].setProperty("name",    correctionsParameters.getNameFromSourceTiff(sourceFiles[nFile]));
					  imp_srcs[srcChannel].setProperty("channel", srcChannel); // it may already have channel
					  imp_srcs[srcChannel].setProperty("path",    sourceFiles[nFile]); // it may already have channel

					  if (this.correctionsParameters.pixelDefects && (eyesisCorrections.defectsXY!=null)&& (eyesisCorrections.defectsXY[srcChannel]!=null)){
						  // apply pixel correction
						  int numApplied=	eyesisCorrections.correctDefects(
								  imp_srcs[srcChannel],
								  srcChannel,
								  debugLevel);
						  if ((debugLevel>0) && (numApplied>0)) { // reduce verbosity after verified defect correction works
							  System.out.println("Corrected "+numApplied+" pixels in "+sourceFiles[nFile]);
						  }
					  }

					  if (this.correctionsParameters.vignetting){
						  if ((eyesisCorrections.channelVignettingCorrection==null) || (srcChannel<0) || (srcChannel>=eyesisCorrections.channelVignettingCorrection.length) || (eyesisCorrections.channelVignettingCorrection[srcChannel]==null)){
							  System.out.println("No vignetting data for channel "+srcChannel);
							  return;
						  }
						  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
						  if (pixels.length!=eyesisCorrections.channelVignettingCorrection[srcChannel].length){
							  System.out.println("Vignetting data for channel "+srcChannel+" has "+eyesisCorrections.channelVignettingCorrection[srcChannel].length+" pixels, image "+sourceFiles[nFile]+" has "+pixels.length);
							  return;
						  }
						  // TODO: Move to do it once:
						  double min_non_zero = 0.0;
						  for (int i=0;i<pixels.length;i++){
							  double d = eyesisCorrections.channelVignettingCorrection[srcChannel][i];
							  if ((d > 0.0) && ((min_non_zero == 0) || (min_non_zero > d))){
								  min_non_zero = d;
							  }
						  }
						  double max_vign_corr = clt_parameters.vignetting_range*min_non_zero;

						  System.out.println("Vignetting data: channel="+srcChannel+", min = "+min_non_zero);
						  for (int i=0;i<pixels.length;i++){
							  double d = eyesisCorrections.channelVignettingCorrection[srcChannel][i];
							  if (d > max_vign_corr) d = max_vign_corr;
							  pixels[i]*=d;
						  }
						  // Scale here, combine with vignetting later?
						  int width =  imp_srcs[srcChannel].getWidth();
						  int height = imp_srcs[srcChannel].getHeight();
						  for (int y = 0; y < height-1; y+=2){
							  for (int x = 0; x < width-1; x+=2){
								  pixels[y*width+x        ] *= clt_parameters.scale_g;
								  pixels[y*width+x+width+1] *= clt_parameters.scale_g;
								  pixels[y*width+x      +1] *= clt_parameters.scale_r;
								  pixels[y*width+x+width  ] *= clt_parameters.scale_b;
							  }
						  }

					  } else { // assuming GR/BG pattern
						  System.out.println("Applying fixed color gain correction parameters: Gr="+
								  clt_parameters.novignetting_r+", Gg="+clt_parameters.novignetting_g+", Gb="+clt_parameters.novignetting_b);
						  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
						  int width =  imp_srcs[srcChannel].getWidth();
						  int height = imp_srcs[srcChannel].getHeight();
						  double kr = clt_parameters.scale_r/clt_parameters.novignetting_r;
						  double kg = clt_parameters.scale_g/clt_parameters.novignetting_g;
						  double kb = clt_parameters.scale_b/clt_parameters.novignetting_b;
						  for (int y = 0; y < height-1; y+=2){
							  for (int x = 0; x < width-1; x+=2){
								  pixels[y*width+x        ] *= kg;
								  pixels[y*width+x+width+1] *= kg;
								  pixels[y*width+x      +1] *= kr;
								  pixels[y*width+x+width  ] *= kb;
							  }
						  }
					  }
				  }
			  }
			  // once per quad here
			  // may need to equalize gains between channels
			  if (clt_parameters.gain_equalize || clt_parameters.colors_equalize){
				  channelGainsEqualize(
						  clt_parameters.gain_equalize,
						  clt_parameters.colors_equalize,
						  channelFiles,
						  imp_srcs,
						  setNames.get(nSet), // just for debug messeges == setNames.get(nSet)
						  debugLevel);
			  }
			  // once per quad here
			  processCLTQuadCorr( // returns ImagePlus, but it already should be saved/shown
					  imp_srcs, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
					  clt_parameters,
					  debayerParameters,
					  nonlinParameters,
					  colorProcParameters,
					  channelGainParameters,
					  rgbParameters,
					  convolveFFTSize, // 128 - fft size, kernel size should be size/2
					  scaleExposures,
					  apply_corr, // calculate and apply additional fine geometry correction 
					  threadsMax,  // maximal number of threads to launch                         
					  updateStatus,
					  debugLevel);
			  Runtime.getRuntime().gc();
			  if (debugLevel >-1) System.out.println("Processing set "+(nSet+1)+" (of "+fileIndices.length+") finished at "+
					  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
			  if (eyesisCorrections.stopRequested.get()>0) {
				  System.out.println("User requested stop");
				  return;
			  }
		  }
		  System.out.println("Processing "+fileIndices.length+" files finished at "+
				  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
	  }		

	  public void channelGainsEqualize(
			  boolean gain_equalize,
			  boolean colors_equalize,
			  int [] channelFiles,
			  ImagePlus [] imp_srcs,
			  String setName, // just for debug messeges == setNames.get(nSet)
			  int debugLevel){
		  double [][] avr_pix = new double [channelFiles.length][3];
		  double [] avr_RGB = {0.0,0.0,0.0};
		  int numChn = 0;
		  for (int srcChannel=0; srcChannel < channelFiles.length; srcChannel++){
			  int nFile=channelFiles[srcChannel];
			  if (nFile >=0){
				  for (int i = 0; i < avr_pix[srcChannel].length; i++) avr_pix[srcChannel][i] = 0;
				  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
				  int width =  imp_srcs[srcChannel].getWidth();
				  int height = imp_srcs[srcChannel].getHeight();
				  for (int y = 0; y < height-1; y+=2){
					  for (int x = 0; x < width-1; x+=2){
						  avr_pix[srcChannel][0] += pixels[y*width+x      +1]; 
						  avr_pix[srcChannel][2] += pixels[y*width+x+width  ];
						  avr_pix[srcChannel][1] += pixels[y*width+x        ];
						  avr_pix[srcChannel][1] += pixels[y*width+x+width+1];
					  }
				  }
				  avr_pix[srcChannel][0] /= 0.25*width*height;
				  avr_pix[srcChannel][1] /= 0.5 *width*height;
				  avr_pix[srcChannel][2] /= 0.25*width*height;
				  for (int j=0; j < avr_RGB.length; j++) avr_RGB[j] += avr_pix[srcChannel][j];
				  numChn++;
				  if (debugLevel>-1) {
					  System.out.println("processCLTSets(): set "+ setName + " channel "+srcChannel+
							  " R"+avr_pix[srcChannel][0]+" G"+avr_pix[srcChannel][1]+" B"+avr_pix[srcChannel][2]);
				  }

			  }
		  }
		  for (int j=0; j < avr_RGB.length; j++) avr_RGB[j] /= numChn;
		  if (debugLevel>-1) {
			  System.out.println("processCLTSets(): set "+ setName + "average color values: "+
					  " R="+avr_RGB[0]+" G=" + avr_RGB[1]+" B=" + avr_RGB[2]);
		  }
		  for (int srcChannel=0; srcChannel < channelFiles.length; srcChannel++){
			  int nFile=channelFiles[srcChannel];
			  if (nFile >=0){
				  double [] scales = new double [avr_RGB.length];
				  for (int j=0;j < scales.length; j++){
					  scales[j] = 1.0;
					  if (gain_equalize){
						  scales[j] *=  avr_RGB[1]/avr_pix[srcChannel][1]; // 1 - index of green color
					  }
					  if (colors_equalize){
						  scales[j] *=  avr_RGB[j]/avr_pix[srcChannel][j] / (avr_RGB[1]/avr_pix[srcChannel][1]);
					  }
				  }
				  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
				  int width =  imp_srcs[srcChannel].getWidth();
				  int height = imp_srcs[srcChannel].getHeight();
				  for (int y = 0; y < height-1; y+=2){
					  for (int x = 0; x < width-1; x+=2){
						  pixels[y*width+x        ] *= scales[1];
						  pixels[y*width+x+width+1] *= scales[1];
						  pixels[y*width+x      +1] *= scales[0];
						  pixels[y*width+x+width  ] *= scales[2];
					  }
				  }
			  }
		  }
	  }
	  
	  public ImagePlus [] processCLTQuadCorr(
			  ImagePlus [] imp_quad, // should have properties "name"(base for saving results), "channel","path"
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  int              convolveFFTSize, // 128 - fft size, kernel size should be size/2
			  double []	       scaleExposures, // probably not needed here
			  final boolean    apply_corr, // calculate and apply additional fine geometry correction 
			  final int        threadsMax,  // maximal number of threads to launch                         
			  final boolean    updateStatus,
			  final int        debugLevel){
		  boolean advanced=this.correctionsParameters.zcorrect || this.correctionsParameters.equirectangular;
		  //		  boolean crop=      advanced? true: this.correctionsParameters.crop; 
//		  boolean rotate=    advanced? false: this.correctionsParameters.rotate; 
//		  double JPEG_scale= advanced? 1.0: this.correctionsParameters.JPEG_scale;
		  boolean toRGB=     advanced? true: this.correctionsParameters.toRGB; 
		  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?

		  // may use this.StartTime to report intermediate steps execution times
		  String name=(String) imp_quad[0].getProperty("name");
		  //		int channel= Integer.parseInt((String) imp_src.getProperty("channel"));
//		  int channel= (Integer) imp_quad[0].getProperty("channel");
		  String path= (String) imp_quad[0].getProperty("path");

		  ImagePlus [] results = new ImagePlus[imp_quad.length];
		  for (int i = 0; i < results.length; i++) {
			  results[i] = imp_quad[i];
			  results[i].setTitle(results[i].getTitle()+"RAW");			  
		  }
		  if (debugLevel>1) System.out.println("processing: "+path);
		  double [][][] double_stacks = new double [imp_quad.length][][];
		  for (int i = 0; i < double_stacks.length; i++){
			  double_stacks[i] = eyesisCorrections.bayerToDoubleStack(
					  imp_quad[i], // source Bayer image, linearized, 32-bit (float))
					  null); // no margins, no oversample
		  }

		  //		  String [] rbg_titles = {"Red", "Blue", "Green"};
//		  ImageStack stack;
		  // =================
		  ImageDtt image_dtt = new ImageDtt();
		  for (int i = 0; i < double_stacks.length; i++){
			  for (int j =0 ; j < double_stacks[i][0].length; j++){
				  double_stacks[i][2][j]*=0.5; // Scale green 0.5 to compensate more pixels than R,B
			  }
		  }

		  setTiles (imp_quad[0], // set global tp.tilesX, tp.tilesY
				  clt_parameters,
				  threadsMax);
		  
		  
		  // temporary setting up tile task file (one integer per tile, bitmask
		  // for testing defined for a window, later the tiles to process will be calculated based on previous passes results

		  int [][]    tile_op = tp.setSameTileOp(clt_parameters,  clt_parameters.tile_task_op, debugLevel);
		  double [][] disparity_array = tp.setSameDisparity(clt_parameters.disparity); // 0.0); // [tp.tilesY][tp.tilesX] - individual per-tile expected disparity
		  
		  //TODO: Add array of default disparity - use for combining images in force disparity mode (no correlation), when disparity is predicted from other tiles

		  double [][][][]     clt_corr_combo =   null;
		  double [][][][][]   clt_corr_partial = null; // [tp.tilesY][tp.tilesX][pair][color][(2*transform_size-1)*(2*transform_size-1)]
		  double [][]         clt_mismatch =     null; // [3*4][tp.tilesY * tp.tilesX] // transpose unapplied
		  double [][][][]     texture_tiles =    null; // [tp.tilesY][tp.tilesX]["RGBA".length()][]; // tiles will be 16x16, 2 visualizaion mode full 16 or overlapped
		  // undecided, so 2 modes of combining alpha - same as rgb, or use center tile only
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  
		  if (clt_parameters.correlate){
			  //			  clt_corr_combo =    new double [2][tp.tilesY][tp.tilesX][];
			  clt_corr_combo =    new double [ImageDtt.TCORR_TITLES.length][tilesY][tilesX][];
			  texture_tiles =     new double [tilesY][tilesX][][]; // ["RGBA".length()][];
			  for (int i = 0; i < tilesY; i++){
				  for (int j = 0; j < tilesX; j++){
					  for (int k = 0; k<clt_corr_combo.length; k++){
						  clt_corr_combo[k][i][j] = null;
					  }
					  //					  clt_corr_combo[1][i][j] = null;
					  texture_tiles[i][j] = null;
				  }
			  }
			  if (clt_parameters.corr_keep){
				  clt_corr_partial = new double [tilesY][tilesX][][][];
				  for (int i = 0; i < tilesY; i++){
					  for (int j = 0; j < tilesX; j++){
						  clt_corr_partial[i][j] = null;
					  }
				  }
			  }
			  if (clt_parameters.corr_mismatch || apply_corr){
				  clt_mismatch = new double [12][];
			  }
		  }
		  double [][] disparity_map = new double [ImageDtt.DISPARITY_TITLES.length][]; //[0] -residual disparity, [1] - orthogonal (just for debugging) last 4 - max pixel differences
		  
//		  double min_corr_selected = clt_parameters.corr_normalize? clt_parameters.min_corr_normalized: clt_parameters.min_corr;
		  double min_corr_selected = clt_parameters.min_corr;
		  double [][] shiftXY = {
				  {clt_parameters.fine_corr_x_0,clt_parameters.fine_corr_y_0},
				  {clt_parameters.fine_corr_x_1,clt_parameters.fine_corr_y_1},
				  {clt_parameters.fine_corr_x_2,clt_parameters.fine_corr_y_2},
				  {clt_parameters.fine_corr_x_3,clt_parameters.fine_corr_y_3}};
		  
		  double [][][][][][] clt_data = image_dtt.clt_aberrations_quad_corr(
				  tile_op,                      // per-tile operation bit codes
				  disparity_array,              // final double            disparity,
				  double_stacks,                // final double [][][]      imade_data, // first index - number of image in a quad
				  // correlation results - final and partial          
				  clt_corr_combo,               // [tp.tilesY][tp.tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				  clt_corr_partial,             // [tp.tilesY][tp.tilesX][pair][color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				  clt_mismatch,                 // [12][tp.tilesY * tp.tilesX] // transpose unapplied. null - do not calculate
				  disparity_map,                // [2][tp.tilesY * tp.tilesX]
				  texture_tiles,                // [tp.tilesY][tp.tilesX]["RGBA".length()][]; 			  
				  imp_quad[0].getWidth(),       // final int width,
				  clt_parameters.fat_zero,      // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
				  clt_parameters.corr_sym,
				  clt_parameters.corr_offset,
				  clt_parameters.corr_red,
				  clt_parameters.corr_blue,
				  clt_parameters.corr_sigma,
//				  clt_parameters.corr_mask,
				  clt_parameters.corr_normalize, // normalize correlation results by rms 
				  min_corr_selected, // 0.0001; // minimal correlation value to consider valid 
				  clt_parameters.max_corr_sigma,// 1.5;  // weights of points around global max to find fractional
				  clt_parameters.max_corr_radius,
				  clt_parameters.enhortho_width,  // 2;    // reduce weight of center correlation pixels from center (0 - none, 1 - center, 2 +/-1 from center)
				  clt_parameters.enhortho_scale,  // 0.2;  // multiply center correlation pixels (inside enhortho_width)
				  clt_parameters.max_corr_double, // Double pass when masking center of mass to reduce preference for integer values
				  clt_parameters.corr_mode,     // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
				  clt_parameters.min_shot,       // 10.0;  // Do not adjust for shot noise if lower than
				  clt_parameters.scale_shot,     // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
				  clt_parameters.diff_sigma,     // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_threshold, // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_gauss,     // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
				  clt_parameters.min_agree,      // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				  clt_parameters.dust_remove,    // Do not reduce average weight when only one image differes much from the average
				  clt_parameters.keep_weights,   // Add port weights to RGBA stack (debug feature)
				  geometryCorrection,           // final GeometryCorrection  geometryCorrection,
				  clt_kernels,                  // final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				  clt_parameters.kernel_step,
				  clt_parameters.transform_size,
				  clt_parameters.clt_window,
				  shiftXY, // 
				  (clt_parameters.fcorr_ignore? null: this.fine_corr),
				  clt_parameters.corr_magic_scale, // stil not understood coefficent that reduces reported disparity value.  Seems to be around 8.5 
				  
				  clt_parameters.shift_x,       // final int               shiftX, // shift image horizontally (positive - right) - just for testing
				  clt_parameters.shift_y,       // final int               shiftY, // shift image vertically (positive - down)
				  clt_parameters.tileX,         // final int               debug_tileX,
				  clt_parameters.tileY,         // final int               debug_tileY,
				  (clt_parameters.dbg_mode & 64) != 0, // no fract shift
				  (clt_parameters.dbg_mode & 128) != 0, // no convolve
				  //				  (clt_parameters.dbg_mode & 256) != 0, // transpose convolve
				  threadsMax,
				  debugLevel);
		  if (debugLevel > -1){
			  System.out.println("clt_data.length="+clt_data.length+" clt_data[0].length="+clt_data[0].length
					  +" clt_data[0][0].length="+clt_data[0][0].length+" clt_data[0][0][0].length="+
					  clt_data[0][0][0].length);
		  }
//		  +" clt_data[0][0][0][0].length="+clt_data[0][0][0][0].length+
//				  " clt_data[0][0][0][0][0].length="+clt_data[0][0][0][0][0].length);
		  // visualize texture tiles as RGBA slices
		  double [][] texture_nonoverlap = null;
		  double [][] texture_overlap = null;
		  String [] rgba_titles = {"red","blue","green","alpha"};
		  String [] rgba_weights_titles = {"red","blue","green","alpha","port0","port1","port2","port3","r-rms","b-rms","g-rms","w-rms"};
		  if (texture_tiles != null){
			  if (clt_parameters.show_nonoverlap){
				  texture_nonoverlap = image_dtt.combineRGBATiles(
						  texture_tiles,                 // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}   
						  clt_parameters.transform_size,
						  false,                         // when false - output each tile as 16x16, true - overlap to make 8x8
						  clt_parameters.sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only 
						  threadsMax,                    // maximal number of threads to launch                         
						  debugLevel);
				  sdfa_instance.showArrays(
						  texture_nonoverlap,
						  tilesX * (2 * clt_parameters.transform_size),
						  tilesY * (2 * clt_parameters.transform_size),
						  true,
						  name + "-TXTNOL-D"+clt_parameters.disparity,
						  (clt_parameters.keep_weights?rgba_weights_titles:rgba_titles));

			  }
			  if (clt_parameters.show_overlap || clt_parameters.show_rgba_color){
				  int alpha_index = 3;
				  texture_overlap = image_dtt.combineRGBATiles(
						  texture_tiles,                 // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}   
						  clt_parameters.transform_size,
						  true,                         // when false - output each tile as 16x16, true - overlap to make 8x8
						  clt_parameters.sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only 
						  threadsMax,                    // maximal number of threads to launch                         
						  debugLevel);
				  if (clt_parameters.alpha1 > 0){ // negative or 0 - keep alpha as it was
					  double scale = (clt_parameters.alpha1 > clt_parameters.alpha0) ? (1.0/(clt_parameters.alpha1 - clt_parameters.alpha0)) : 0.0; 
					  for (int i = 0; i < texture_overlap[alpha_index].length; i++){
						  double d = texture_overlap[alpha_index][i];
						  if      (d >=clt_parameters.alpha1) d = 1.0;
						  else if (d <=clt_parameters.alpha0) d = 0.0;
						  else d = scale * (d- clt_parameters.alpha0);
						  texture_overlap[alpha_index][i] = d;
					  }
				  }
				  
				  if (clt_parameters.show_overlap) {
					  sdfa_instance.showArrays(
							  texture_overlap,
							  tilesX * clt_parameters.transform_size,
							  tilesY * clt_parameters.transform_size,
							  true,
							  name + "-TXTOL-D"+clt_parameters.disparity,
							  (clt_parameters.keep_weights?rgba_weights_titles:rgba_titles));
				  }
				  if (clt_parameters.show_rgba_color) {
					  // for now - use just RGB. Later add oprion for RGBA
					  double [][] texture_rgb = {texture_overlap[0],texture_overlap[1],texture_overlap[2]};
					  double [][] texture_rgba = {texture_overlap[0],texture_overlap[1],texture_overlap[2],texture_overlap[3]};
//					  ImagePlus img_texture = 
					  linearStackToColor(
							  clt_parameters,
							  colorProcParameters,
							  rgbParameters,
							  name+"-texture", // String name,
							  "-D"+clt_parameters.disparity, //String suffix, // such as disparity=...
							  toRGB,
							  !this.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
							  true, // boolean saveShowIntermediate, // save/show if set globally
							  true, // boolean saveShowFinal,        // save/show result (color image?)
							  ((clt_parameters.alpha1 > 0)? texture_rgba: texture_rgb),
							  tilesX *  clt_parameters.transform_size,
							  tilesY *  clt_parameters.transform_size,
							  1.0,         // double scaleExposure, // is it needed?
							  debugLevel );
				  }
			  }
		  }
		  // visualize correlation results
		  if (clt_corr_combo!=null){
			  if (disparity_map != null){
				  if (clt_parameters.show_map){
					  sdfa_instance.showArrays(
							  disparity_map,
							  tilesX,
							  tilesY,
							  true,
							  name+"-DISP_MAP-D"+clt_parameters.disparity,
							  ImageDtt.DISPARITY_TITLES);
				  }			  
			  }

			  if (clt_mismatch != null){
				  if (debugLevel > -1){
					  String [] disparity_titles = {"dx0", "dy0","strength0","dx1", "dy1","strength1","dx2", "dy2","strength2","dx3", "dy3","strength3"};
					  sdfa_instance.showArrays(
							  clt_mismatch,
							  tilesX,
							  tilesY,
							  true,
							  name+"-MISMATCH_XYW-D"+clt_parameters.disparity,
							  disparity_titles);
				  }
				  if (apply_corr && (disparity_map != null)){
					  System.out.println("=== applying geometry correction coefficients ===");
					  double [][][] new_corr = fine_geometry_correction(
							  clt_parameters,
							  disparity_map,
							  clt_mismatch,
							  tilesX,
							  clt_parameters.corr_magic_scale, // stil not understood coefficent that reduces reported disparity value.  Seems to be around 8.5  
							  debugLevel); // int debugLevel)
					  apply_fine_corr(
							  new_corr,
							  debugLevel + 1);
				  }
			  }


			  if (clt_parameters.corr_show){
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

			  if (clt_corr_partial!=null){
				  if (debugLevel > -1){ // -1
					  String [] allColorNames = {"red","blue","green","combo"};
					  String [] titles = new String[clt_corr_partial.length];
					  for (int i = 0; i < titles.length; i++){
						  titles[i]=allColorNames[i % allColorNames.length]+"_"+(i / allColorNames.length);
					  }
					  double [][] corr_rslt_partial = image_dtt.corr_partial_dbg(
							  clt_corr_partial,
							  2*clt_parameters.transform_size - 1,	//final int corr_size,
							  4,	// final int pairs,
							  4,    // final int colors,
							  clt_parameters.corr_border_contrast,
							  threadsMax,
							  debugLevel);
					  sdfa_instance.showArrays(
							  corr_rslt_partial,
							  tilesX*(2*clt_parameters.transform_size),
							  tilesY*(2*clt_parameters.transform_size),
							  true,
							  name+"-PART_CORR-D"+clt_parameters.disparity,
							  titles);
				  }
			  }
		  }

		  if (clt_parameters.gen_chn_img || clt_parameters.gen_chn_stacks) {
			  ImagePlus [] imps_RGB = new ImagePlus[clt_data.length];
			  for (int iQuad = 0; iQuad < clt_data.length; iQuad++){

				  String title=name+"-"+String.format("%02d", iQuad);
//				  String titleFull=title+"-SPLIT-D"+clt_parameters.disparity;

				  if (clt_parameters.corr_sigma > 0){ // no filter at all
					  for (int chn = 0; chn < clt_data[iQuad].length; chn++) {
						  image_dtt.clt_lpf(
								  clt_parameters.corr_sigma,
								  clt_data[iQuad][chn],
								  clt_parameters.transform_size,
								  threadsMax,
								  debugLevel);
					  }
				  }

				  //			  int tp.tilesY = imp_quad[iQuad].getHeight()/clt_parameters.transform_size;
				  //			  int tp.tilesX = imp_quad[iQuad].getWidth()/clt_parameters.transform_size;
				  if (debugLevel > 0){
					  System.out.println("--tilesX="+tilesX);
					  System.out.println("--tilesY="+tilesY);
				  }
				  if (debugLevel > 0){
					  double [][] clt = new double [clt_data[iQuad].length*4][];
					  for (int chn = 0; chn < clt_data[iQuad].length; chn++) {
						  double [][] clt_set = image_dtt.clt_dbg(
								  clt_data [iQuad][chn],
								  threadsMax,
								  debugLevel);
						  for (int ii = 0; ii < clt_set.length; ii++) clt[chn*4+ii] = clt_set[ii];
					  }

					  if (debugLevel > 0){
						  sdfa_instance.showArrays(clt,
								  tilesX*clt_parameters.transform_size,
								  tilesY*clt_parameters.transform_size,
								  true,
								  results[iQuad].getTitle()+"-CLT-D"+clt_parameters.disparity);  
					  }
				  }
				  double [][] iclt_data = new double [clt_data[iQuad].length][];
				  for (int chn=0; chn<iclt_data.length;chn++){
					  iclt_data[chn] = image_dtt.iclt_2d(
							  clt_data[iQuad][chn],           // scanline representation of dcd data, organized as dct_size x dct_size tiles  
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
						  results[iQuad].getTitle()+"-ICLT-RGB-D"+clt_parameters.disparity);
				  if (!clt_parameters.gen_chn_img) continue;
				  
				  imps_RGB[iQuad] = linearStackToColor(
						  clt_parameters,
						  colorProcParameters,
						  rgbParameters,
						  title, // String name,
						  "-D"+clt_parameters.disparity, //String suffix, // such as disparity=...
						  toRGB,
						  !this.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
						  true, // boolean saveShowIntermediate, // save/show if set globally
						  false, // boolean saveShowFinal,        // save/show result (color image?)
						  iclt_data,
						  tilesX *  clt_parameters.transform_size,
						  tilesY *  clt_parameters.transform_size,
						  scaleExposures[iQuad], // double scaleExposure, // is it needed?
						  debugLevel );
				  
			  } // end of generating shifted channel images
			  if (clt_parameters.gen_chn_img) {
				  // combine to a sliced color image
				  int [] slice_seq = {0,1,3,2}; //clockwise
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
				  imp_stack.updateAndDraw();
				  //imp_stack.getProcessor().resetMinAndMax();
				  //imp_stack.show();
				  eyesisCorrections.saveAndShow(imp_stack, this.correctionsParameters);
			  }
		  }
		  return results;
	  }

	  double [][] resizeGridTexture(
			  double [][] imgData,
			  int tileSize,
			  int tilesX,
			  int tilesY,
			  Rectangle bounds){
//		  int width =   tileSize*bounds.width + 1; 
//		  int height =  tileSize*bounds.height + 1;
		  int width =   tileSize*bounds.width; 
		  int height =  tileSize*bounds.height;
		  // Adding row/column of all 0, assuming java zeroes arrays
		  int numSlices = imgData.length;
		  double [][] rslt = new double [numSlices][width*height];
		  int offset = (tileSize*bounds.y)*tileSize*tilesX + (tileSize*bounds.x); 
//		  System.out.println("bounds = {.x:"+bounds.x + ", .y=" + bounds.y + ", .width=" + bounds.width + ", .height=" + bounds.height+ ", numSlices="+numSlices);
//		  System.out.println("offset = " + offset+ " tileSize = "+tileSize);
//		  for (int y = 0; y < (height - 1); y ++){
//			  for (int x = 0; x < width-1; x ++){
		  for (int y = 0; y < height; y ++){
			  for (int x = 0; x < width; x ++){
				  int indx = width * y + x; 
				  int indx_in = indx + offset;
				  for (int i = 0; i < numSlices; i++) {
					  rslt[i][indx] = Double.isNaN(imgData[i][indx_in])? 0.0:imgData[i][indx_in];
//					  if ((imgData[i][indx_in] !=0.0) && (i ==3))
//					  {
//						  int a = 0;
//						  int b = a;
//					  }
				  }
			  }
			  offset += tilesX * tileSize - width;
			  
		  }
		  return rslt;
	  }
	  
	  public ImagePlus linearStackToColor(
			  EyesisCorrectionParameters.CLTParameters         clt_parameters,
			  EyesisCorrectionParameters.ColorProcParameters   colorProcParameters,
			  EyesisCorrectionParameters.RGBParameters         rgbParameters,
			  String name,
			  String suffix, // such as disparity=...
			  boolean toRGB,
			  boolean bpp16, // 16-bit per channel color mode for result
			  boolean saveShowIntermediate, // save/show if set globally
			  boolean saveShowFinal,        // save/show result (color image?)
			  double [][] iclt_data,
			  int width, // int tilesX,
			  int height, // int tilesY,
			  double scaleExposure,
			  int debugLevel
			  )
	  {
		  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		  // convert to ImageStack of 3 slices
		  String [] sliceNames = {"red", "blue", "green"};
		  double []   alpha = null; // (0..1.0)
		  double [][] rgb_in = {iclt_data[0],iclt_data[1],iclt_data[2]};
		  if (iclt_data.length > 3) alpha = iclt_data[3]; 
		  ImageStack stack = sdfa_instance.makeStack(
				  rgb_in, // iclt_data,
				  width,  // (tilesX + 0) * clt_parameters.transform_size,
				  height, // (tilesY + 0) * clt_parameters.transform_size,
				  sliceNames,  // or use null to get chn-nn slice names
				  true); // replace NaN with 0.0
		  if (debugLevel > 0){
			  double [] chn_avg = {0.0,0.0,0.0};
			  float [] pixels;
//			  int width =  stack.getWidth();
//			  int height = stack.getHeight();

			  for (int c = 0; c <3; c++){
				  pixels = (float[]) stack.getPixels(c+1);
				  for (int i = 0; i<pixels.length; i++){
					  chn_avg[c] += pixels[i];
				  }
			  }
			  chn_avg[0] /= width*height;
			  chn_avg[1] /= width*height;
			  chn_avg[2] /= width*height;
			  System.out.println("Processed channels averages: R="+chn_avg[0]+", G="+chn_avg[2]+", B="+chn_avg[1]); 
		  }

		  if (debugLevel > 1) System.out.println("before colors.1");
		  //Processing colors - changing stack sequence to r-g-b (was r-b-g)
		  if (!eyesisCorrections.fixSliceSequence(
				  stack,
				  debugLevel)){
			  if (debugLevel > -1) System.out.println("fixSliceSequence() returned false");
			  return null;
		  }
		  if (debugLevel > 1) System.out.println("before colors.2");
		  if (saveShowIntermediate && (debugLevel > 1)){
			  ImagePlus imp_dbg=new ImagePlus(name+"-preColors",stack);
			  eyesisCorrections.saveAndShow(
					  imp_dbg,
					  this.correctionsParameters);
		  }
		  if (debugLevel > 1) System.out.println("before colors.3, scaleExposure="+scaleExposure+" scale = "+(255.0/eyesisCorrections.psfSubpixelShouldBe4/eyesisCorrections.psfSubpixelShouldBe4/scaleExposure));
		  CorrectionColorProc correctionColorProc=new CorrectionColorProc(eyesisCorrections.stopRequested);
		  double [][] yPrPb=new double [3][];
		  //			if (dct_parameters.color_DCT){
		  // need to get YPbPr - not RGB here				
		  //			} else {
		  correctionColorProc.processColorsWeights(stack, // just gamma convert? TODO: Cleanup? Convert directly form the linear YPrPb
				  //					  255.0/this.psfSubpixelShouldBe4/this.psfSubpixelShouldBe4, //  double scale,     // initial maximal pixel value (16))
				  //					  255.0/eyesisCorrections.psfSubpixelShouldBe4/eyesisCorrections.psfSubpixelShouldBe4/scaleExposure, //  double scale,     // initial maximal pixel value (16))
				  //					  255.0/2/2/scaleExposure, //  double scale,     // initial maximal pixel value (16))
				  255.0/scaleExposure, //  double scale,     // initial maximal pixel value (16))
				  colorProcParameters,
				  null, // channelGainParameters,
				  -1, // channel,
				  null, //correctionDenoise.getDenoiseMask(),
				  this.correctionsParameters.blueProc,
				  debugLevel);
		  if (debugLevel > 1) System.out.println("Processed colors to YPbPr, total number of slices="+stack.getSize());
		  if (saveShowIntermediate && (debugLevel > 1)) {
			  ImagePlus imp_dbg=new ImagePlus("procColors",stack);
			  eyesisCorrections.saveAndShow(
					  imp_dbg,
					  this.correctionsParameters);
		  }
		  float [] fpixels;
		  int [] slices_YPrPb = {8,6,7};
		  yPrPb=new double [3][];
		  for (int n = 0; n < slices_YPrPb.length; n++){
			  fpixels = (float[]) stack.getPixels(slices_YPrPb[n]);
			  yPrPb[n] = new double [fpixels.length];
			  for (int i = 0; i < fpixels.length; i++) yPrPb[n][i] = fpixels[i];
		  }
		  String titleFull = "";
		  if (toRGB) {
			  System.out.println("correctionColorProc.YPrPbToRGB");
			  stack =  YPrPbToRGB(yPrPb,
					  colorProcParameters.kr,        // 0.299;
					  colorProcParameters.kb,        // 0.114;
					  stack.getWidth());
			  titleFull=name+"-RGB-float"+suffix;
			  //Trim stack to just first 3 slices
			  if (saveShowIntermediate && (debugLevel > 1)){ // 2){
				  ImagePlus imp_dbg=new ImagePlus("YPrPbToRGB",stack);
				  eyesisCorrections.saveAndShow(
						  imp_dbg,
						  this.correctionsParameters);
			  }
			  while (stack.getSize() > 3) stack.deleteLastSlice();
			  if (debugLevel > 1) System.out.println("Trimming color stack");
		  } else {
			  titleFull=name+"-YPrPb"+suffix;
			  if (debugLevel > 1) System.out.println("Using full stack, including YPbPr");
		  }
		  if (alpha != null){
			  float [] alpha_pixels = new float [alpha.length];
			  for (int i = 0; i <alpha.length; i++){
				  alpha_pixels[i] = (float) alpha[i];
			  }
			  stack.addSlice("alpha",alpha_pixels);
		  }
		  
		  ImagePlus result= new ImagePlus(titleFull, stack);    			  

		  if (!toRGB && !this.correctionsParameters.jpeg){ // toRGB set for equirectangular
			  if (debugLevel > 1) System.out.println("!toRGB && !this.correctionsParameters.jpeg");
			  if (saveShowIntermediate) eyesisCorrections.saveAndShow(result, this.correctionsParameters);
			  return result;
		  } else { // that's not the end result, save if required
			  if (debugLevel > 1) System.out.println("!toRGB && !this.correctionsParameters.jpeg - else");
			  
			  
			  if (saveShowIntermediate) eyesisCorrections.saveAndShow(result, // saves OK with alpha - 32-bit float
					  eyesisCorrections.correctionsParameters,
					  eyesisCorrections.correctionsParameters.save32,
					  false, // true, // false,
					  eyesisCorrections.correctionsParameters.JPEG_quality); // save, no show
		  }
		  stack=eyesisCorrections.convertRGB32toRGB16Stack(
				  stack,
				  rgbParameters); 

		  titleFull=name+"-RGB48"+suffix;
		  result= new ImagePlus(titleFull, stack);    			  
		  result.updateAndDraw();
		  if (debugLevel > 1) System.out.println("result.updateAndDraw(), "+titleFull+"-RGB48");

		  CompositeImage compositeImage=eyesisCorrections.convertToComposite(result);

		  if (!this.correctionsParameters.jpeg && bpp16){ // RGB48 was the end result
			  if (debugLevel > 1) System.out.println("if (!this.correctionsParameters.jpeg && !advanced)");
			  if (saveShowIntermediate) eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters);
			  return compositeImage; // return result;
		  } else { // that's not the end result, save if required
			  if (debugLevel > 1) System.out.println("if (!this.correctionsParameters.jpeg && !advanced) - else");
			  if (saveShowIntermediate) eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters, this.correctionsParameters.save16, false); // save, no show
		  }
		  result = eyesisCorrections.convertRGB48toRGB24(
				  stack,
//				  name+"-RGB24"+suffix,
				  name+suffix,
				  0, 65536, // r range 0->0, 65536->256
				  0, 65536, // g range
				  0, 65536,// b range
				  0, 65536);// alpha range
		  // next will save either JPEG (if no alpha) or RGBA tiff (if alpha is present). ImageJ shows just RGB (no alpha)
		  if (saveShowFinal) eyesisCorrections.saveAndShow(result, this.correctionsParameters);
		  return result;

	  }

	  
	  public double [][][] fine_geometry_correction(
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  double [][] disparity_map,
			  double [][] clt_mismatch,
			  int         tilesX,
			  double      magic_coeff, // stil not understood coefficent that reduces reported disparity value.  Seems to be around 8.5  
			  int debugLevel)
	  {
		  // TODO: update the following indices when format of the arrays will change (removed unneeded data)
		  final int index_disparity =     ImageDtt.DISPARITY_INDEX_CM; // 2; // "cm_disparity" 
		  final int index_strength =      ImageDtt.DISPARITY_STRENGTH_INDEX; //6; // "strength"
		  final int [] indices_mismatch = {1,4,6,9}; //   "dy0", "dy1", "dx2", "dx3"
		  final int numTiles =            disparity_map[0].length;
		  final int tilesY =              numTiles/tilesX;

		  int numUsed = 0;
		  for (int i = 0; i < numTiles; i++) {
			  if ((disparity_map[index_strength][i] >=  clt_parameters.fcorr_min_stength) &&
					  (Math.abs(disparity_map[index_disparity][i]) <= clt_parameters.fcorr_disp_diff)) numUsed ++;
		  }
		  double [][][] mdata = new double[numUsed][3][];
		  int indx = 0;
		  for (int i = 0; i < numTiles; i++) {
			  if ((disparity_map[index_strength][i] >=  clt_parameters.fcorr_min_stength) &&
					  (Math.abs(disparity_map[index_disparity][i]) <= clt_parameters.fcorr_disp_diff)) {
				  int tileX = i % tilesX;
				  int tileY = i / tilesX;
				  mdata[indx][0] = new double [2];
				  mdata[indx][0][0] =  (2.0 * tileX)/tilesX - 1.0; // -1.0 to +1.0;
				  mdata[indx][0][1] =  (2.0 * tileY)/tilesY - 1.0; // -1.0 to +1.0
				  mdata[indx][1] = new double [indices_mismatch.length]; // 4
				  for (int ip = 0; ip < indices_mismatch.length; ip++) {
					  mdata[indx][1][ip] =  clt_mismatch[indices_mismatch[ip]][i];
				  }
				  mdata[indx][2] = new double [1];
				  mdata[indx][2][0] =  disparity_map[index_strength][i];
				  indx ++;
			  }
		  }
		  PolynomialApproximation pa = new PolynomialApproximation();
		  double thresholdLin = 1.0E-20;  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail)
		  double thresholdQuad = 1.0E-30; // threshold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)

		  /*
		   * returns array of vectors or null
		   * each vector (one per each z component) is either 6-element-  (A,B,C,D,E,F) if quadratic is possible and enabled
		   * or 3-element - (D,E,F) if linear is possible and quadratic is not possible or disabled
		   * returns null if not enough data even for the linear approximation
		   */
		  
		  double [][] coeffs = pa.quadraticApproximation(
				  mdata,
				  !clt_parameters.fcorr_quadratic, // boolean forceLinear,  // use linear approximation
				  thresholdLin,  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail)
				  thresholdQuad,  // threshold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
				  debugLevel);
		  for (int i = 0; i < coeffs.length; i++){
			  if (coeffs[i] == null){
				  coeffs[i] = new double [6];
				  for (int j = 0; j < coeffs[i].length; i++) coeffs[i][j] = 0.0;
			  } else if (coeffs[i].length < 6){
				  double [] short_coeffs = coeffs[i];
				  coeffs[i] = new double [6];
				  for (int j = 0; j < coeffs[i].length; j++) {
					  if (j < (coeffs[i].length - short_coeffs.length)) coeffs[i][j] = 0.0;
					  else coeffs[i][j] = short_coeffs[j - (coeffs[i].length - short_coeffs.length)]; 
				  }
			  }
		  }
		  // convert to 8 sets of coefficient for px0, py0, px1, py1, ... py3.  
		  double [][][] coeff_full = new double [4][2][6];
		  double scale = 0.5/magic_coeff;
		  for (int j = 0; j<6; j++){
			  coeff_full[0][0][j] = -coeffs[2][j] * scale;
			  coeff_full[0][1][j] = -coeffs[0][j] * scale;
			  coeff_full[1][0][j] = -coeffs[3][j] * scale;
			  coeff_full[1][1][j] =  coeffs[0][j] * scale;
			  coeff_full[2][0][j] =  coeffs[2][j] * scale;
			  coeff_full[2][1][j] = -coeffs[1][j] * scale;
			  coeff_full[3][0][j] =  coeffs[3][j] * scale;
			  coeff_full[3][1][j] =  coeffs[1][j] * scale;
		  }
		  if (debugLevel > -1){
			  int tilesX8 = (tilesX - 1) / 8 + 1;
			  int tilesY8 = (tilesY - 1) / 8 + 1;
			  int len8 = tilesX8*tilesY8;
			  double [][] mismatch_8 = new double[11][len8];
			  for (int i = 0; i < numTiles; i++) {
				  if ((disparity_map[index_strength][i] >=  clt_parameters.fcorr_min_stength) &&
						  (Math.abs(disparity_map[index_disparity][i]) <= clt_parameters.fcorr_disp_diff)) {
					  int tileX = i % tilesX;
					  int tileY = i / tilesX;
					  int tX8 = tileX/8;
					  int tY8 = tileY/8;
					  int indx8 = tY8*tilesX8+tX8;
					  double tX =  (2.0 * tileX)/tilesX - 1.0; // -1.0 to +1.0;
					  double tY =  (2.0 * tileY)/tilesY - 1.0; // -1.0 to +1.0
					  double w = disparity_map[index_strength][i];
					  for (int ip = 0; ip < 4; ip++) {
						  mismatch_8[ip][indx8] += w * clt_mismatch[indices_mismatch[ip]][i];
					  }
					  mismatch_8[8][indx8] += w * tX;
					  mismatch_8[9][indx8] += w * tY;
					  mismatch_8[10][indx8] += w;
				  }
			  }
			  for (int i = 0; i < len8; i++) {
				  if (mismatch_8[10][i] > 0){
					  for (int n = 0; n<4; n++){
						  mismatch_8[n][i] /= mismatch_8[10][i]; 
					  }
					  mismatch_8[8][i] /= mismatch_8[10][i]; 
					  mismatch_8[9][i] /= mismatch_8[10][i]; 
				  } else {
					  for (int n = 0; n<4; n++){
						  mismatch_8[n][i] = Double.NaN; 
					  }
					  mismatch_8[8][i] = Double.NaN; 
					  mismatch_8[9][i] = Double.NaN; 
					  
				  }
				  for (int n = 0; n<4; n++){
					  double tX = mismatch_8[8][i];
					  double tY = mismatch_8[9][i];
					  //f(x,y)=A*x^2+B*y^2+C*x*y+D*x+E*y+F
					  mismatch_8[4+n][i] = (
							  coeffs[n][0]*tX*tX+
							  coeffs[n][1]*tY*tY+
							  coeffs[n][2]*tX*tY+
							  coeffs[n][3]*tX+
							  coeffs[n][4]*tY+
							  coeffs[n][5]);
				  }
			  }
			  String [] titles = {"dy0","dy1","dx2","dx3","cy0","cy1","cx2","cx3","tx","ty","w"};
			  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
			  sdfa_instance.showArrays(
					  mismatch_8,
					  tilesX8,
					  tilesY8,
					  true,
					  "mismatch_8",
					  titles);
		  }
		  
		  
		  return coeff_full;
	  }
	  
	  public void apply_fine_corr(
			  double [][][] corr,
			  int debugLevel)
	  {
		  if (debugLevel > 1){
			  if (debugLevel > 1){
				  show_fine_corr( this.fine_corr, "  was");
			  }
		  }
		  if (debugLevel > 1){
			  if (debugLevel > 1){
				  show_fine_corr(corr, "added");
			  }
		  }
		  
		  for (int n = 0; n< corr.length; n++){
			  for (int i = 0; i< corr[n].length; i++){
				  for (int j = 0; j< corr[n][i].length; j++){
					  this.fine_corr[n][i][j]+=corr[n][i][j];
				  }
			  }
		  }
		  if (debugLevel > 0){
			  show_fine_corr( this.fine_corr, "");
		  }
	  }
	  public void show_fine_corr()
	  {
		  show_fine_corr( this.fine_corr, "");
	  }
	  
	  public void show_fine_corr(
			  double [][][] corr,
			  String prefix)
	  {
		  String sadd = (prefix.length() > 0)?(prefix+" "):""; 

		  for (int n = 0; n< corr.length; n++){
			  for (int i = 0; i< corr[n].length; i++){
				  System.out.print(sadd+"port"+n+": "+fine_corr_dir_names[i]+": ");
				  for (int j = 0; j< corr[n][i].length; j++){
					  System.out.print(fine_corr_coeff_names[j]+"="+corr[n][i][j]+" ");
				  }
				  System.out.println();
			  }
		  }
	  }
	  
	  public void reset_fine_corr()
	  {
		  this.fine_corr = new double [4][2][6]; // reset all coefficients to 0  
	  }
	  
	  public void cltDisparityScans(
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
			  int          convolveFFTSize, // 128 - fft size, kernel size should be size/2
			  final int          threadsMax,  // maximal number of threads to launch                         
			  final boolean    updateStatus,
			  final int        debugLevel)
	  {
		  this.startTime=System.nanoTime();
		  String [] sourceFiles=correctionsParameters.getSourcePaths();
		  boolean [] enabledFiles=new boolean[sourceFiles.length];
		  for (int i=0;i<enabledFiles.length;i++) enabledFiles[i]=false;
		  int numFilesToProcess=0;
		  int numImagesToProcess=0;
		  for (int nFile=0;nFile<enabledFiles.length;nFile++){
			  if ((sourceFiles[nFile]!=null) && (sourceFiles[nFile].length()>1)) {
				  int [] channels={correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile])};
				  if (correctionsParameters.isJP4()){
					  int subCamera= channels[0]- correctionsParameters.firstSubCamera; // to match those in the sensor files
					  // removeUnusedSensorData should be off!?
					  channels=this.eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (eyesisCorrections.isChannelEnabled(channels[i])){
						  if (!enabledFiles[nFile]) numFilesToProcess++;
						  enabledFiles[nFile]=true;
						  numImagesToProcess++;
					  }
				  }
			  }
		  }
		  if (numFilesToProcess==0){
			  System.out.println("No files to process (of "+sourceFiles.length+")");
			  return;
		  } else {
			  if (debugLevel>0) System.out.println(numFilesToProcess+ " files to process (of "+sourceFiles.length+"), "+numImagesToProcess+" images to process");
		  }
		  double [] referenceExposures=eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		  int [][] fileIndices=new int [numImagesToProcess][2]; // file index, channel number
		  int index=0;
		  for (int nFile=0;nFile<enabledFiles.length;nFile++){
			  if ((sourceFiles[nFile]!=null) && (sourceFiles[nFile].length()>1)) {
				  int [] channels={correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile])};
				  if (correctionsParameters.isJP4()){
					  int subCamera= channels[0]- correctionsParameters.firstSubCamera; // to match those in the sensor files
					  channels=eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (eyesisCorrections.isChannelEnabled(channels[i])){
						  fileIndices[index  ][0]=nFile;
						  fileIndices[index++][1]=channels[i];
					  }
				  }
			  }
		  }
		  ArrayList<String> setNames = new ArrayList<String>();
		  ArrayList<ArrayList<Integer>> setFiles = new ArrayList<ArrayList<Integer>>();

		  for (int iImage=0;iImage<fileIndices.length;iImage++){
			  int nFile=fileIndices[iImage][0];
			  String setName = correctionsParameters.getNameFromSourceTiff(sourceFiles[nFile]);
			  if (!setNames.contains(setName)) {
				  setNames.add(setName);
				  setFiles.add(new ArrayList<Integer>());
			  }
			  setFiles.get(setNames.indexOf(setName)).add(new Integer(nFile));
		  }		  
		  for (int nSet = 0; nSet < setNames.size(); nSet++){
			  int maxChn = 0;
			  for (int i = 0; i < setFiles.get(nSet).size(); i++){
				  int chn = fileIndices[setFiles.get(nSet).get(i)][1];
				  if (chn > maxChn) maxChn = chn;
			  }
			  int [] channelFiles = new int[maxChn+1];
			  for (int i =0; i < channelFiles.length; i++) channelFiles[i] = -1;
			  for (int i = 0; i < setFiles.get(nSet).size(); i++){
				  channelFiles[fileIndices[setFiles.get(nSet).get(i)][1]] = setFiles.get(nSet).get(i);
			  }

			  ImagePlus [] imp_srcs = new ImagePlus[channelFiles.length];
			  double [] scaleExposures = new double[channelFiles.length];
			  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
				  int nFile=channelFiles[srcChannel];
				  imp_srcs[srcChannel]=null;
				  if (nFile >=0){
					  if (correctionsParameters.isJP4()){
						  int subchannel=eyesisCorrections.pixelMapping.getSubChannel(srcChannel);
						  if (this.correctionsParameters.swapSubchannels01) {
							  switch (subchannel){
							  case 0: subchannel=1; break;
							  case 1: subchannel=0; break;
							  }
						  }
						  if (debugLevel>0) System.out.println("Processing set " + setNames.get(nSet)+" channel "+srcChannel+" - subchannel "+subchannel+" of "+sourceFiles[nFile]);
						  ImagePlus imp_composite=eyesisCorrections.JP4_INSTANCE.open(
								  "", // path,
								  sourceFiles[nFile],
								  "",  //arg - not used in JP46 reader
								  true, // un-apply camera color gains
								  null, // new window
								  false); // do not show
						  imp_srcs[srcChannel]=eyesisCorrections.JP4_INSTANCE.demuxImage(imp_composite, subchannel);
						  if (imp_srcs[srcChannel] == null) imp_srcs[srcChannel] = imp_composite; // not a composite image
						  // do we need to add any properties?				  
					  } else { 
						  imp_srcs[srcChannel]=new ImagePlus(sourceFiles[nFile]);
						  //					  (new JP46_Reader_camera(false)).decodeProperiesFromInfo(imp_src); // decode existent properties from info
						  eyesisCorrections.JP4_INSTANCE.decodeProperiesFromInfo(imp_srcs[srcChannel]); // decode existent properties from info
						  if (debugLevel>0) System.out.println("Processing "+sourceFiles[nFile]);
					  }
					  scaleExposures[srcChannel] = 1.0;
					  if (!Double.isNaN(referenceExposures[nFile]) && (imp_srcs[srcChannel].getProperty("EXPOSURE")!=null)){
						  scaleExposures[srcChannel] = referenceExposures[nFile]/Double.parseDouble((String) imp_srcs[srcChannel].getProperty("EXPOSURE"));
						  if (debugLevel>0) System.out.println("Will scale intensity (to compensate for exposure) by  "+scaleExposures[srcChannel]);
					  }
					  imp_srcs[srcChannel].setProperty("name",    correctionsParameters.getNameFromSourceTiff(sourceFiles[nFile]));
					  imp_srcs[srcChannel].setProperty("channel", srcChannel); // it may already have channel
					  imp_srcs[srcChannel].setProperty("path",    sourceFiles[nFile]); // it may already have channel

					  if (this.correctionsParameters.pixelDefects && (eyesisCorrections.defectsXY!=null)&& (eyesisCorrections.defectsXY[srcChannel]!=null)){
						  // apply pixel correction
						  int numApplied=	eyesisCorrections.correctDefects(
								  imp_srcs[srcChannel],
								  srcChannel,
								  debugLevel);
						  if ((debugLevel>0) && (numApplied>0)) { // reduce verbosity after verified defect correction works
							  System.out.println("Corrected "+numApplied+" pixels in "+sourceFiles[nFile]);
						  }
					  }

					  if (this.correctionsParameters.vignetting){
						  if ((eyesisCorrections.channelVignettingCorrection==null) || (srcChannel<0) || (srcChannel>=eyesisCorrections.channelVignettingCorrection.length) || (eyesisCorrections.channelVignettingCorrection[srcChannel]==null)){
							  System.out.println("No vignetting data for channel "+srcChannel);
							  return;
						  }
						  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
						  if (pixels.length!=eyesisCorrections.channelVignettingCorrection[srcChannel].length){
							  System.out.println("Vignetting data for channel "+srcChannel+" has "+eyesisCorrections.channelVignettingCorrection[srcChannel].length+" pixels, image "+sourceFiles[nFile]+" has "+pixels.length);
							  return;
						  }
						  // TODO: Move to do it once:
						  double min_non_zero = 0.0;
						  for (int i=0;i<pixels.length;i++){
							  double d = eyesisCorrections.channelVignettingCorrection[srcChannel][i];
							  if ((d > 0.0) && ((min_non_zero == 0) || (min_non_zero > d))){
								  min_non_zero = d;
							  }
						  }
						  double max_vign_corr = clt_parameters.vignetting_range*min_non_zero;

						  System.out.println("Vignetting data: channel="+srcChannel+", min = "+min_non_zero);
						  for (int i=0;i<pixels.length;i++){
							  double d = eyesisCorrections.channelVignettingCorrection[srcChannel][i];
							  if (d > max_vign_corr) d = max_vign_corr;
							  pixels[i]*=d;
						  }
						  // Scale here, combine with vignetting later?
						  int width =  imp_srcs[srcChannel].getWidth();
						  int height = imp_srcs[srcChannel].getHeight();
						  for (int y = 0; y < height-1; y+=2){
							  for (int x = 0; x < width-1; x+=2){
								  pixels[y*width+x        ] *= clt_parameters.scale_g;
								  pixels[y*width+x+width+1] *= clt_parameters.scale_g;
								  pixels[y*width+x      +1] *= clt_parameters.scale_r;
								  pixels[y*width+x+width  ] *= clt_parameters.scale_b;
							  }
						  }

					  } else { // assuming GR/BG pattern
						  System.out.println("Applying fixed color gain correction parameters: Gr="+
								  clt_parameters.novignetting_r+", Gg="+clt_parameters.novignetting_g+", Gb="+clt_parameters.novignetting_b);
						  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
						  int width =  imp_srcs[srcChannel].getWidth();
						  int height = imp_srcs[srcChannel].getHeight();
						  double kr = clt_parameters.scale_r/clt_parameters.novignetting_r;
						  double kg = clt_parameters.scale_g/clt_parameters.novignetting_g;
						  double kb = clt_parameters.scale_b/clt_parameters.novignetting_b;
						  for (int y = 0; y < height-1; y+=2){
							  for (int x = 0; x < width-1; x+=2){
								  pixels[y*width+x        ] *= kg;
								  pixels[y*width+x+width+1] *= kg;
								  pixels[y*width+x      +1] *= kr;
								  pixels[y*width+x+width  ] *= kb;
							  }
						  }
					  }
				  }
			  }
			  // once per quad here
			  // may need to equalize gains between channels
			  if (clt_parameters.gain_equalize || clt_parameters.colors_equalize){
				  channelGainsEqualize(
						  clt_parameters.gain_equalize,
						  clt_parameters.colors_equalize,
						  channelFiles,
						  imp_srcs,
						  setNames.get(nSet), // just for debug messeges == setNames.get(nSet)
						  debugLevel);
			  }
			  // once per quad here
			  cltDisparityScan( // returns ImagePlus, but it already should be saved/shown
					  imp_srcs, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
					  clt_parameters,
					  debayerParameters,
					  nonlinParameters,
					  colorProcParameters,
					  channelGainParameters,
					  rgbParameters,
					  convolveFFTSize, // 128 - fft size, kernel size should be size/2
					  scaleExposures,
					  threadsMax,  // maximal number of threads to launch                         
					  updateStatus,
					  debugLevel);
			  Runtime.getRuntime().gc();
			  if (debugLevel >-1) System.out.println("Processing set "+(nSet+1)+" (of "+fileIndices.length+") finished at "+
					  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
			  if (eyesisCorrections.stopRequested.get()>0) {
				  System.out.println("User requested stop");
				  return;
			  }
		  }
		  System.out.println("Processing "+fileIndices.length+" files finished at "+
				  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
	  }		


	  public ImagePlus [] cltDisparityScan(
			  ImagePlus [] imp_quad, // should have properties "name"(base for saving results), "channel","path"
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  int              convolveFFTSize, // 128 - fft size, kernel size should be size/2
			  double []	       scaleExposures, // probably not needed here
			  final int        threadsMax,  // maximal number of threads to launch                         
			  final boolean    updateStatus,
			  final int        debugLevel){
		  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?

		  // may use this.StartTime to report intermediate steps execution times
		  String name=(String) imp_quad[0].getProperty("name");
		  //		int channel= Integer.parseInt((String) imp_src.getProperty("channel"));
		  String path= (String) imp_quad[0].getProperty("path");

		  ImagePlus [] results = new ImagePlus[imp_quad.length];
		  for (int i = 0; i < results.length; i++) {
			  results[i] = imp_quad[i];
			  results[i].setTitle(results[i].getTitle()+"RAW");			  
		  }
		  if (debugLevel>1) System.out.println("processing: "+path);
		  double [][][] double_stacks = new double [imp_quad.length][][];
		  for (int i = 0; i < double_stacks.length; i++){
			  double_stacks[i] = eyesisCorrections.bayerToDoubleStack(
					  imp_quad[i], // source Bayer image, linearized, 32-bit (float))
					  null); // no margins, no oversample
		  }
		  // =================
		  ImageDtt image_dtt = new ImageDtt();
		  for (int i = 0; i < double_stacks.length; i++){
			  for (int j =0 ; j < double_stacks[i][0].length; j++){
				  double_stacks[i][2][j]*=0.5; // Scale green 0.5 to compensate more pixels than R,B
			  }
		  }
		  
		  setTiles (imp_quad[0], // set global tp.tilesX, tp.tilesY
				  clt_parameters,
				  threadsMax);
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  
		  // temporary setting up tile task file (one integer per tile, bitmask
		  // for testing defined for a window, later the tiles to process will be calculated based on previous passes results
		  int [][]    tile_op = tp.setSameTileOp(clt_parameters,  clt_parameters.tile_task_op, debugLevel);
		  
		  //TODO: Add array of default disparity - use for combining images in force disparity mode (no correlation), when disparity is predicted from other tiles

		  // undecided, so 2 modes of combining alpha - same as rgb, or use center tile only
		  double [][][][]     clt_corr_combo =    new double [ImageDtt.TCORR_TITLES.length][tilesY][tilesX][]; // will only be used inside?
		  double min_corr_selected = clt_parameters.min_corr;
		  
		  double [][][] disparity_maps = new double [clt_parameters.disp_scan_count][ImageDtt.DISPARITY_TITLES.length][]; //[0] -residual disparity, [1] - orthogonal (just for debugging)

		  for (int scan_step = 0; scan_step < clt_parameters.disp_scan_count; scan_step++) {
			  double disparity = clt_parameters.disp_scan_start + scan_step * clt_parameters.disp_scan_step;
			  double [][] disparity_array = tp.setSameDisparity(disparity); // [tp.tilesY][tp.tilesX] - individual per-tile expected disparity
			  
			  double [][] shiftXY = {
					  {clt_parameters.fine_corr_x_0,clt_parameters.fine_corr_y_0},
					  {clt_parameters.fine_corr_x_1,clt_parameters.fine_corr_y_1},
					  {clt_parameters.fine_corr_x_2,clt_parameters.fine_corr_y_2},
					  {clt_parameters.fine_corr_x_3,clt_parameters.fine_corr_y_3}};

			  image_dtt.clt_aberrations_quad_corr(
					  tile_op,                      // per-tile operation bit codes
					  disparity_array,              // clt_parameters.disparity,     // final double            disparity,
					  double_stacks,                // final double [][][]      imade_data, // first index - number of image in a quad
					  // correlation results - final and partial          
					  clt_corr_combo,               // [tp.tilesY][tp.tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					  null, // clt_corr_partial,    // [tp.tilesY][tp.tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					  null,    // [tp.tilesY][tp.tilesX][pair]{dx,dy,weight}[(2*transform_size-1)*(2*transform_size-1)] // transpose unapplied. null - do not calculate
//	Use it with disparity_maps[scan_step]?		  clt_mismatch,    // [tp.tilesY][tp.tilesX][pair]{dx,dy,weight}[(2*transform_size-1)*(2*transform_size-1)] // transpose unapplied. null - do not calculate
					  disparity_maps[scan_step],    // [2][tp.tilesY * tp.tilesX]
					  null, //texture_tiles,        // [tp.tilesY][tp.tilesX]["RGBA".length()][]; 			  
					  imp_quad[0].getWidth(),       // final int width,
					  clt_parameters.fat_zero,      // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
					  clt_parameters.corr_sym,
					  clt_parameters.corr_offset,
					  clt_parameters.corr_red,
					  clt_parameters.corr_blue,
					  clt_parameters.corr_sigma,
					  clt_parameters.corr_normalize, // normalize correlation results by rms 
					  min_corr_selected, // 0.0001; // minimal correlation value to consider valid 
					  clt_parameters.max_corr_sigma,// 1.5;  // weights of points around global max to find fractional
					  clt_parameters.max_corr_radius,
					  clt_parameters.enhortho_width,  // 2;    // reduce weight of center correlation pixels from center (0 - none, 1 - center, 2 +/-1 from center)
					  clt_parameters.enhortho_scale,  // 0.2;  // multiply center correlation pixels (inside enhortho_width)
					  clt_parameters.max_corr_double, // Double pass when masking center of mass to reduce preference for integer values
					  clt_parameters.corr_mode,     // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
					  clt_parameters.min_shot,       // 10.0;  // Do not adjust for shot noise if lower than
					  clt_parameters.scale_shot,     // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
					  clt_parameters.diff_sigma,     // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
					  clt_parameters.diff_threshold, // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
					  clt_parameters.diff_gauss,     // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
					  clt_parameters.min_agree,      // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
					  clt_parameters.dust_remove,    // Do not reduce average weight when only one image differes much from the average
					  clt_parameters.keep_weights,   // Add port weights to RGBA stack (debug feature)
					  geometryCorrection,           // final GeometryCorrection  geometryCorrection,
					  clt_kernels,                  // final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
					  clt_parameters.kernel_step,
					  clt_parameters.transform_size,
					  clt_parameters.clt_window,
					  shiftXY, // 
					  (clt_parameters.fcorr_ignore? null: this.fine_corr),
					  clt_parameters.corr_magic_scale, // stil not understood coefficent that reduces reported disparity value.  Seems to be around 8.5 
					  clt_parameters.shift_x,       // final int               shiftX, // shift image horizontally (positive - right) - just for testing
					  clt_parameters.shift_y,       // final int               shiftY, // shift image vertically (positive - down)
					  clt_parameters.tileX,         // final int               debug_tileX,
					  clt_parameters.tileY,         // final int               debug_tileY,
					  (clt_parameters.dbg_mode & 64) != 0, // no fract shift
					  (clt_parameters.dbg_mode & 128) != 0, // no convolve
					  //				  (clt_parameters.dbg_mode & 256) != 0, // transpose convolve
					  threadsMax,
					  debugLevel);
		  }
		  int [] disp_indices = {
				  ImageDtt.DISPARITY_INDEX_INT,
				  ImageDtt.DISPARITY_INDEX_CM,
				  ImageDtt.DISPARITY_INDEX_POLY,
				  ImageDtt.DISPARITY_STRENGTH_INDEX,
				  ImageDtt.DISPARITY_VARIATIONS_INDEX};
		  String [] disparity_titles = new String [disp_indices.length];
		  for (int i = 0; i < disparity_titles.length; i++ ) disparity_titles[i] = ImageDtt.DISPARITY_TITLES[i];
		  
//				  2,4,6,7};
		  String [] disparities_titles = new String [disparity_titles.length * clt_parameters.disp_scan_count];
		  double [][] disparities_maps = new double [disparity_titles.length * clt_parameters.disp_scan_count][];
		  int indx = 0;

		  for (int scan_step = 0; scan_step < clt_parameters.disp_scan_count; scan_step++) {
			  double disparity = clt_parameters.disp_scan_start + scan_step * clt_parameters.disp_scan_step;
			  for (int i = 0; i < disparity_titles.length; i++){
				  disparities_titles[indx] = disparity_titles[i]+"_"+disparity;
				  disparities_maps[indx++] = disparity_maps[scan_step][disp_indices[i]]; 
			  }
		  }
		  
		  ImageStack array_stack = sdfa_instance.makeStack(
				  disparities_maps,
				  tilesX,
				  tilesY,
				  disparities_titles);
			  		  
		  ImagePlus imp_stack = new ImagePlus( name+"-DISP_MAPS", array_stack);
		  imp_stack.getProcessor().resetMinAndMax();
		  imp_stack.updateAndDraw();
		  //imp_stack.getProcessor().resetMinAndMax();
		  //imp_stack.show();
		  eyesisCorrections.saveAndShow(imp_stack, this.correctionsParameters);
// process scan results
		  double [][] scan_trends = process_disparity_scan(
				  disparities_maps,
				  clt_parameters.disp_scan_step,
				  clt_parameters.disp_scan_start,
				  0.0); // min_corr_selected); // all what is remaining 
		  String [] trend_titles={"b_cm","b_poly","a_cm","a_poly","rms_cm","rms_poly","strength","samples"};

		  ImageStack trends_stack = sdfa_instance.makeStack(
				  scan_trends,
				  tilesX,
				  tilesY,
				  trend_titles);
		  ImagePlus imp_stack_trends = new ImagePlus( name+"-DISP_TRENDS", trends_stack);
		  imp_stack_trends.getProcessor().resetMinAndMax();
		  imp_stack_trends.updateAndDraw();
		  eyesisCorrections.saveAndShow(imp_stack_trends, this.correctionsParameters);
		  
		  return results;
	  }
	  
	  public double [][] process_disparity_scan(
			  double [][] disparities_maps,
			  double disp_step,
			  double disp_start,
			  double min_strength)
	  {
		  final int num_items = 5;
		  final int index_strength = 3;
		  final int index_cm = 1;
		  final int index_poly = 1;
		  
		  final int ind_b_cm =   0;
		  final int ind_b_poly = 1;
		  final int ind_a_cm =   2;
		  final int ind_a_poly = 3;
		  final int ind_rms_cm =   4;
		  final int ind_rms_poly = 5;
		  final int ind_strength =   6;
		  final int ind_samples = 7;
		  
		  final int num_steps = disparities_maps.length / num_items; // int, cm, poly, strength, variety
		  final int disp_len =  disparities_maps[0].length;
		  double [][] rslt = new double [8][disp_len];
		  for (int i = 0; i < disp_len; i++){
			  double s0 = 0.0, sx = 0.0, sx2 = 0.0, sy_cm = 0.0, sxy_cm = 0.0, sy_poly = 0.0, sxy_poly = 0.0;
			  int samples = 0;
			  for (int step = 0; step < num_steps; step++ ){
				  double wi = disparities_maps[num_items*step + index_strength][i];
				  if (wi > min_strength) {
					  double xi =      disp_start + step * disp_step;
					  double yi_cm =   disparities_maps[num_items*step + index_cm][i];
					  double yi_poly = disparities_maps[num_items*step + index_poly][i];
					  if (Double.isNaN(yi_cm) || Double.isNaN(yi_poly)) continue;
					  s0 +=     wi;
					  sx +=     wi*xi;
					  sx2 +=    wi*xi*xi;
					  sy_cm +=  wi*yi_cm;
					  sxy_cm += wi*xi*yi_cm;
					  sy_poly +=  wi*yi_poly;
					  sxy_poly += wi*xi*yi_poly;
					  samples++;
				  }
			  }
			  double denom = (s0*sx2 - sx*sx);
			  rslt[ind_strength][i] = s0; 
			  rslt[ind_samples][i] =  samples;
			  if (denom != 0.0) {
				  rslt[ind_a_cm][i] =     (s0*sxy_cm - sx*sy_cm) /  denom;
				  rslt[ind_b_cm][i] =     (sy_cm*sx2 - sx*sxy_cm) / denom;
				  rslt[ind_a_poly][i] =   (s0*sxy_poly - sx*sy_poly) /  denom;
				  rslt[ind_b_poly][i] =   (sy_poly*sx2 - sx*sxy_poly) / denom;
				  rslt[ind_rms_cm][i] =   0.0;
				  rslt[ind_rms_poly][i] = 0.0;
				  for (int step = 0; step < num_steps; step++ ){
					  double wi = disparities_maps[num_items*step + index_strength][i];
					  if (wi > min_strength) {
						  double xi =      disp_start + step * disp_step;
						  double yi_cm =   disparities_maps[num_items*step + index_cm][i];
						  double yi_poly = disparities_maps[num_items*step + index_poly][i];
						  if (Double.isNaN(yi_cm) || Double.isNaN(yi_poly)) continue;
						  double d_cm =   yi_cm -   (rslt[ind_a_cm][i]*xi +rslt[ind_b_cm][i]); 
						  double d_poly = yi_poly - (rslt[ind_a_poly][i]*xi +rslt[ind_b_poly][i]);
						  rslt[ind_rms_cm][i] +=   wi*d_cm*d_cm;
						  rslt[ind_rms_poly][i] += wi*d_poly*d_poly;
					  }
				  }
				  rslt[ind_rms_cm][i] =   Math.sqrt(rslt[ind_rms_cm][i]/s0);
				  rslt[ind_rms_poly][i] = Math.sqrt(rslt[ind_rms_cm][i]/s0);
			  } else {
				  rslt[ind_a_cm][i] =     Double.NaN;
				  rslt[ind_b_cm][i] =     Double.NaN;
				  rslt[ind_a_poly][i] =   Double.NaN;
				  rslt[ind_b_poly][i] =   Double.NaN;
				  rslt[ind_rms_cm][i] =   Double.NaN;
				  rslt[ind_rms_poly][i] = Double.NaN;
			  }
			  
		  }
		  return rslt;
	  }
//	  public ImagePlus [] cltDisparityScan(
	  
	  public void processCLTQuads3d(
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
			  final int          threadsMax,  // maximal number of threads to launch                         
			  final boolean    updateStatus,
			  final int        debugLevel)
	  {
		  this.startTime=System.nanoTime();
		  String [] sourceFiles=correctionsParameters.getSourcePaths();
		  boolean [] enabledFiles=new boolean[sourceFiles.length];
		  for (int i=0;i<enabledFiles.length;i++) enabledFiles[i]=false;
		  int numFilesToProcess=0;
		  int numImagesToProcess=0;
		  for (int nFile=0;nFile<enabledFiles.length;nFile++){
			  if ((sourceFiles[nFile]!=null) && (sourceFiles[nFile].length()>1)) {
				  int [] channels={correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile])};
				  if (correctionsParameters.isJP4()){
					  int subCamera= channels[0]- correctionsParameters.firstSubCamera; // to match those in the sensor files
					  // removeUnusedSensorData should be off!?
					  channels=this.eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (eyesisCorrections.isChannelEnabled(channels[i])){
						  if (!enabledFiles[nFile]) numFilesToProcess++;
						  enabledFiles[nFile]=true;
						  numImagesToProcess++;
					  }
				  }
			  }
		  }
		  if (numFilesToProcess==0){
			  System.out.println("No files to process (of "+sourceFiles.length+")");
			  return;
		  } else {
			  if (debugLevel>0) System.out.println(numFilesToProcess+ " files to process (of "+sourceFiles.length+"), "+numImagesToProcess+" images to process");
		  }
		  double [] referenceExposures=eyesisCorrections.calcReferenceExposures(debugLevel); // multiply each image by this and divide by individual (if not NaN)
		  int [][] fileIndices=new int [numImagesToProcess][2]; // file index, channel number
		  int index=0;
		  for (int nFile=0;nFile<enabledFiles.length;nFile++){
			  if ((sourceFiles[nFile]!=null) && (sourceFiles[nFile].length()>1)) {
				  int [] channels={correctionsParameters.getChannelFromSourceTiff(sourceFiles[nFile])};
				  if (correctionsParameters.isJP4()){
					  int subCamera= channels[0]- correctionsParameters.firstSubCamera; // to match those in the sensor files
					  channels=eyesisCorrections.pixelMapping.channelsForSubCamera(subCamera);
				  }
				  if (channels!=null){
					  for (int i=0;i<channels.length;i++) if (eyesisCorrections.isChannelEnabled(channels[i])){
						  fileIndices[index  ][0]=nFile;
						  fileIndices[index++][1]=channels[i];
					  }
				  }
			  }
		  }
		  ArrayList<String> setNames = new ArrayList<String>();
		  ArrayList<ArrayList<Integer>> setFiles = new ArrayList<ArrayList<Integer>>();

		  for (int iImage=0;iImage<fileIndices.length;iImage++){
			  int nFile=fileIndices[iImage][0];
			  String setName = correctionsParameters.getNameFromSourceTiff(sourceFiles[nFile]);
			  if (!setNames.contains(setName)) {
				  setNames.add(setName);
				  setFiles.add(new ArrayList<Integer>());
			  }
			  setFiles.get(setNames.indexOf(setName)).add(new Integer(nFile));
		  }		  
		  for (int nSet = 0; nSet < setNames.size(); nSet++){
			  int maxChn = 0;
			  for (int i = 0; i < setFiles.get(nSet).size(); i++){
				  int chn = fileIndices[setFiles.get(nSet).get(i)][1];
				  if (chn > maxChn) maxChn = chn;
			  }
			  int [] channelFiles = new int[maxChn+1];
			  for (int i =0; i < channelFiles.length; i++) channelFiles[i] = -1;
			  for (int i = 0; i < setFiles.get(nSet).size(); i++){
				  channelFiles[fileIndices[setFiles.get(nSet).get(i)][1]] = setFiles.get(nSet).get(i);
			  }

			  ImagePlus [] imp_srcs = new ImagePlus[channelFiles.length];
			  double [] scaleExposures = new double[channelFiles.length];
			  for (int srcChannel=0; srcChannel<channelFiles.length; srcChannel++){
				  int nFile=channelFiles[srcChannel];
				  imp_srcs[srcChannel]=null;
				  if (nFile >=0){
					  if (correctionsParameters.isJP4()){
						  int subchannel=eyesisCorrections.pixelMapping.getSubChannel(srcChannel);
						  if (this.correctionsParameters.swapSubchannels01) {
							  switch (subchannel){
							  case 0: subchannel=1; break;
							  case 1: subchannel=0; break;
							  }
						  }
						  if (debugLevel>0) System.out.println("Processing set " + setNames.get(nSet)+" channel "+srcChannel+" - subchannel "+subchannel+" of "+sourceFiles[nFile]);
						  ImagePlus imp_composite=eyesisCorrections.JP4_INSTANCE.open(
								  "", // path,
								  sourceFiles[nFile],
								  "",  //arg - not used in JP46 reader
								  true, // un-apply camera color gains
								  null, // new window
								  false); // do not show
						  imp_srcs[srcChannel]=eyesisCorrections.JP4_INSTANCE.demuxImage(imp_composite, subchannel);
						  if (imp_srcs[srcChannel] == null) imp_srcs[srcChannel] = imp_composite; // not a composite image
						  // do we need to add any properties?				  
					  } else { 
						  imp_srcs[srcChannel]=new ImagePlus(sourceFiles[nFile]);
						  //					  (new JP46_Reader_camera(false)).decodeProperiesFromInfo(imp_src); // decode existent properties from info
						  eyesisCorrections.JP4_INSTANCE.decodeProperiesFromInfo(imp_srcs[srcChannel]); // decode existent properties from info
						  if (debugLevel>0) System.out.println("Processing "+sourceFiles[nFile]);
					  }
					  scaleExposures[srcChannel] = 1.0;
					  if (!Double.isNaN(referenceExposures[nFile]) && (imp_srcs[srcChannel].getProperty("EXPOSURE")!=null)){
						  scaleExposures[srcChannel] = referenceExposures[nFile]/Double.parseDouble((String) imp_srcs[srcChannel].getProperty("EXPOSURE"));
						  if (debugLevel>0) System.out.println("Will scale intensity (to compensate for exposure) by  "+scaleExposures[srcChannel]);
					  }
					  imp_srcs[srcChannel].setProperty("name",    correctionsParameters.getNameFromSourceTiff(sourceFiles[nFile]));
					  imp_srcs[srcChannel].setProperty("channel", srcChannel); // it may already have channel
					  imp_srcs[srcChannel].setProperty("path",    sourceFiles[nFile]); // it may already have channel

					  if (this.correctionsParameters.pixelDefects && (eyesisCorrections.defectsXY!=null)&& (eyesisCorrections.defectsXY[srcChannel]!=null)){
						  // apply pixel correction
						  int numApplied=	eyesisCorrections.correctDefects(
								  imp_srcs[srcChannel],
								  srcChannel,
								  debugLevel);
						  if ((debugLevel>0) && (numApplied>0)) { // reduce verbosity after verified defect correction works
							  System.out.println("Corrected "+numApplied+" pixels in "+sourceFiles[nFile]);
						  }
					  }

					  if (this.correctionsParameters.vignetting){
						  if ((eyesisCorrections.channelVignettingCorrection==null) || (srcChannel<0) || (srcChannel>=eyesisCorrections.channelVignettingCorrection.length) || (eyesisCorrections.channelVignettingCorrection[srcChannel]==null)){
							  System.out.println("No vignetting data for channel "+srcChannel);
							  return;
						  }
						  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
						  if (pixels.length!=eyesisCorrections.channelVignettingCorrection[srcChannel].length){
							  System.out.println("Vignetting data for channel "+srcChannel+" has "+eyesisCorrections.channelVignettingCorrection[srcChannel].length+" pixels, image "+sourceFiles[nFile]+" has "+pixels.length);
							  return;
						  }
						  // TODO: Move to do it once:
						  double min_non_zero = 0.0;
						  for (int i=0;i<pixels.length;i++){
							  double d = eyesisCorrections.channelVignettingCorrection[srcChannel][i];
							  if ((d > 0.0) && ((min_non_zero == 0) || (min_non_zero > d))){
								  min_non_zero = d;
							  }
						  }
						  double max_vign_corr = clt_parameters.vignetting_range*min_non_zero;

						  System.out.println("Vignetting data: channel="+srcChannel+", min = "+min_non_zero);
						  for (int i=0;i<pixels.length;i++){
							  double d = eyesisCorrections.channelVignettingCorrection[srcChannel][i];
							  if (d > max_vign_corr) d = max_vign_corr;
							  pixels[i]*=d;
						  }
						  // Scale here, combine with vignetting later?
						  int width =  imp_srcs[srcChannel].getWidth();
						  int height = imp_srcs[srcChannel].getHeight();
						  for (int y = 0; y < height-1; y+=2){
							  for (int x = 0; x < width-1; x+=2){
								  pixels[y*width+x        ] *= clt_parameters.scale_g;
								  pixels[y*width+x+width+1] *= clt_parameters.scale_g;
								  pixels[y*width+x      +1] *= clt_parameters.scale_r;
								  pixels[y*width+x+width  ] *= clt_parameters.scale_b;
							  }
						  }

					  } else { // assuming GR/BG pattern
						  System.out.println("Applying fixed color gain correction parameters: Gr="+
								  clt_parameters.novignetting_r+", Gg="+clt_parameters.novignetting_g+", Gb="+clt_parameters.novignetting_b);
						  float [] pixels=(float []) imp_srcs[srcChannel].getProcessor().getPixels();
						  int width =  imp_srcs[srcChannel].getWidth();
						  int height = imp_srcs[srcChannel].getHeight();
						  double kr = clt_parameters.scale_r/clt_parameters.novignetting_r;
						  double kg = clt_parameters.scale_g/clt_parameters.novignetting_g;
						  double kb = clt_parameters.scale_b/clt_parameters.novignetting_b;
						  for (int y = 0; y < height-1; y+=2){
							  for (int x = 0; x < width-1; x+=2){
								  pixels[y*width+x        ] *= kg;
								  pixels[y*width+x+width+1] *= kg;
								  pixels[y*width+x      +1] *= kr;
								  pixels[y*width+x+width  ] *= kb;
							  }
						  }
					  }
				  }
			  }
			  // once per quad here
			  // may need to equalize gains between channels
			  if (clt_parameters.gain_equalize || clt_parameters.colors_equalize){
				  channelGainsEqualize(
						  clt_parameters.gain_equalize,
						  clt_parameters.colors_equalize,
						  channelFiles,
						  imp_srcs,
						  setNames.get(nSet), // just for debug messeges == setNames.get(nSet)
						  debugLevel);
			  }
			  // once per quad here
			  processCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
					  imp_srcs, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
					  clt_parameters,
					  debayerParameters,
					  nonlinParameters,
					  colorProcParameters,
					  channelGainParameters,
					  rgbParameters,
					  threadsMax,  // maximal number of threads to launch                         
					  updateStatus,
					  debugLevel);
			  Runtime.getRuntime().gc();
			  if (debugLevel >-1) System.out.println("Processing set "+(nSet+1)+" (of "+fileIndices.length+") finished at "+
					  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
			  if (eyesisCorrections.stopRequested.get()>0) {
				  System.out.println("User requested stop");
				  return;
			  }
		  }
		  System.out.println("Processing "+fileIndices.length+" files finished at "+
				  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
	  }		

	  public ImagePlus processCLTQuad3d(
			  ImagePlus [] imp_quad, // should have properties "name"(base for saving results), "channel","path"
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			  EyesisCorrectionParameters.NonlinParameters       nonlinParameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  final int        threadsMax,  // maximal number of threads to launch                         
			  final boolean    updateStatus,
			  final int        debugLevel)
	  {
		  
		  String name = (String) imp_quad[0].getProperty("name");
		  double [][][] image_data = new double [imp_quad.length][][];
		  for (int i = 0; i < image_data.length; i++){
			  image_data[i] = eyesisCorrections.bayerToDoubleStack(
					  imp_quad[i], // source Bayer image, linearized, 32-bit (float))
					  null); // no margins, no oversample
		  }
		  for (int i = 0; i < image_data.length; i++){
			  for (int j =0 ; j < image_data[i][0].length; j++){
				  image_data[i][2][j]*=0.5; // Scale green 0.5 to compensate more pixels than R,B
			  }
		  }
		  setTiles (imp_quad[0], // set global tp.tilesX, tp.tilesY
				  clt_parameters,
				  threadsMax);
		  tp.resetCLTPasses();
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  
		  
		  TileProcessor.CLTPass3d bgnd_data = CLTBackgroundMeas( // measure background
				  image_data, // 
				  clt_parameters,
				  threadsMax,  // maximal number of threads to launch                         
				  updateStatus,
				  debugLevel);
		  tp.clt_3d_passes.add(bgnd_data);
		  
		  ImagePlus imp_bgnd = getBackgroundImage(
				  clt_parameters,
				  colorProcParameters,
				  rgbParameters,
				  name,               // .getTitle(), //String name=(String) imp_src.getProperty("name");
				  ImageDtt.DISPARITY_INDEX_CM, // index of disparity value in disparity_map == 2 (0,2 or 4)
				  threadsMax,  // maximal number of threads to launch                         
				  updateStatus,
				  debugLevel);
		  bgnd_data.texture = imp_bgnd.getTitle()+".png";
		  
// create x3d file
		  X3dOutput x3dOutput = new X3dOutput(
					clt_parameters,
					correctionsParameters,
					geometryCorrection,
					tp.clt_3d_passes);
		  
		  x3dOutput.generateBackground();
 		  String x3d_path= correctionsParameters.selectX3dDirectory(
				  true,  // smart,
				  true);  //newAllowed, // save
 		  
		  // testing 2-nd pass
		  tp.secondPassSetup( // prepare tile tasks for the second pass based on the previous one(s)
//				  final double [][][]       image_data, // first index - number of image in a quad
				  clt_parameters,
				  // disparity range - differences from 
				  clt_parameters.bgnd_range, // double            disparity_far,  
				  clt_parameters.other_range, //double            disparity_near,   //
				  clt_parameters.bgnd_sure, // double            this_sure,        // minimal strength to be considered definitely background
				  clt_parameters.bgnd_maybe, // double            this_maybe,       // maximal strength to ignore as non-background
				  clt_parameters.sure_smth, // sure_smth,        // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
				  ImageDtt.DISPARITY_INDEX_CM,  // index of disparity value in disparity_map == 2 (0,2 or 4)
				  geometryCorrection,
				  threadsMax,  // maximal number of threads to launch                         
				  updateStatus,
				  debugLevel);
		  
		  // get images for predefined regions and dispariteies. First - with just fixed scans 1 .. list.size()
		  for (int scanIndex = 1; scanIndex < tp.clt_3d_passes.size(); scanIndex++){
			  if (debugLevel > 0){
				  System.out.println("FPGA processing scan #"+scanIndex);
			  }
				  TileProcessor.CLTPass3d scan = CLTMeasure( // perform single pass according to prepared tiles operations and disparity
						  image_data, // first index - number of image in a quad
						  clt_parameters,
						  scanIndex,
						  threadsMax,  // maximal number of threads to launch                         
						  updateStatus,
						  debugLevel);

			  
		  }
//		  int scan_limit = 10; 
		  for (int scanIndex = 1; (scanIndex < tp.clt_3d_passes.size()) && (scanIndex < clt_parameters.max_clusters); scanIndex++){ // just temporary limiting
			  if (debugLevel > -1){
				  System.out.println("Generating cluster images (limit is set to "+clt_parameters.max_clusters+") largest, scan #"+scanIndex);
			  }
//			  ImagePlus cluster_image = getPassImage( // get image form a single pass
			  String texturePath = getPassImage( // get image form a single pass
					  clt_parameters,
					  colorProcParameters,
					  rgbParameters,
					  name+"-img"+scanIndex,
					  scanIndex, 
					  threadsMax,  // maximal number of threads to launch                         
					  updateStatus,
					  debugLevel);
			  
			  TileProcessor.CLTPass3d scan = tp.clt_3d_passes.get(scanIndex);
			  
			  // TODO: use new updated disparity, for now just what was forced for the picture
			  double [] scan_disparity = new double [tilesX * tilesY];
			  int indx = 0;
			  for (int ty = 0; ty < tilesY; ty ++) for (int tx = 0; tx < tilesX; tx ++){
				  scan_disparity[indx++] = scan.disparity[ty][tx];
			  }
			  if (clt_parameters.avg_cluster_disp){
				  double sw = 0.0, sdw = 0.0;
				  for (int i = 0; i< scan_disparity.length; i++){
					  if (scan.selected[i] && !scan.border_tiles[i]){
						  double w = scan.disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX][i];
						  sw +=w;
						  sdw += scan_disparity[i]*w;
					  }
				  }
				  sdw/=sw;
				  for (int i = 0; i< scan_disparity.length; i++){
					  scan_disparity[i] = sdw;
				  }
				  
				  
			  }
			  generateClusterX3d(
					  x3dOutput,
					  texturePath,
					  scan.bounds,
					  scan.selected,
					  scan_disparity, // scan.disparity_map[ImageDtt.DISPARITY_INDEX_CM],
					  clt_parameters.transform_size,
					  clt_parameters.correct_distortions, // requires backdrop image to be corrected also
					  (scanIndex <2) && clt_parameters.show_triangles,
					  clt_parameters.bgnd_range,
					  clt_parameters.other_range,
					  clt_parameters.maxDispTriangle); 
			  
		  }

		  // now generate and save texture files (start with full, later use bounding rectangle?)
		  
		  
		  
		  
		  
 		  if (x3d_path != null){
 			 x3d_path+=Prefs.getFileSeparator()+name+".x3d";
 			  x3dOutput.generateX3D(x3d_path);
 		  }
 		  
		  return imp_bgnd; // relative (to x3d directory) path - (String) imp_bgnd.getProperty("name");
	  }

	  public void generateClusterX3d(
			  X3dOutput  x3dOutput,
			  String     texturePath,
			  Rectangle  bounds,
			  boolean [] selected,
			  double []  disparity,
			  int        tile_size,
			  boolean    correctDistortions, // requires backdrop image to be corrected also
			  boolean    show_triangles,
			  double     min_disparity,
			  double     max_disparity,
			  double     maxDispTriangle
			  ) 
	  {
		  int [][] indices =  tp.getCoordIndices( // starting with 0, -1 - not selected
				  bounds,
				  selected);
		  double [][] texCoord = tp.getTexCoords( // get texture coordinates for indices
				  indices);
		  double [][] worldXYZ = tp.getCoords( // get world XYZ in meters for indices
				  disparity,
				  min_disparity,
				  max_disparity,
				  bounds,
				  indices,
				  tile_size,
				  correctDistortions, // requires backdrop image to be corrected also
				  this.geometryCorrection);
				  
          double [] indexedDisparity = tp.getIndexedDisparities( // get disparity for each index
							disparity,
							min_disparity,
							max_disparity,
							bounds,
							indices,
							tile_size);

		  int [][] triangles = 	tp.triangulateCluster(
				  indices);
		  
		  
		  triangles = 	tp.filterTriangles(
					triangles,
					indexedDisparity, // disparities per vertex index
					maxDispTriangle); // maximal disparity difference in a triangle
		  
		  
		  if (show_triangles) {
			  tp.testTriangles(
					  texturePath,
					  bounds,
					  selected,
					  disparity,
					  tile_size,
					  indices,
					  triangles);
		  }

		  x3dOutput.addCluster(
				  texturePath,
				  texCoord,
				  worldXYZ,
				  triangles);
	  }
	  
	  
	  public ImagePlus getBackgroundImage(
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  String     name,
			  int        disparity_index, // index of disparity value in disparity_map == 2 (0,2 or 4)
//			  int        strength_index, // index of strength data in disparity map ==6
			  int        threadsMax,  // maximal number of threads to launch                         
			  boolean    updateStatus,
			  int        debugLevel
			  )
	  {
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  showDoubleFloatArrays sdfa_instance = null;
		  
		  if (clt_parameters.debug_filters && (debugLevel > -1)) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		  
		  TileProcessor.CLTPass3d bgnd_data = tp.clt_3d_passes.get(0); 
		  double [][][][] texture_tiles = bgnd_data.texture_tiles;
		  
		  boolean [] bgnd_tiles =   tp.getBackgroundMask( // which tiles do belong to the background
				  clt_parameters.bgnd_range,     // disparity range to be considered background
				  clt_parameters.bgnd_sure,      // minimal strength to be considered definitely background
				  clt_parameters.bgnd_maybe,     // maximal strength to ignore as non-background
				  clt_parameters.sure_smth,      // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
				  clt_parameters.min_clstr_seed, // number of tiles in a cluster to seed (just background?)
				  clt_parameters.min_clstr_block,// number of tiles in a cluster to block (just non-background?)
				  disparity_index, // index of disparity value in disparity_map == 2 (0,2 or 4)
				  (clt_parameters.debug_filters ? debugLevel : -1));
		  boolean [] bgnd_strict = bgnd_tiles.clone(); // only these have non 0 alpha
		  tp.growTiles(
				  clt_parameters.bgnd_grow,      // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				  bgnd_tiles,
				  null); // prohibit
		  boolean [] bgnd_tiles_grown = bgnd_tiles.clone(); // only these have non 0 alpha
		  tp.growTiles(
				  2,      // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				  bgnd_tiles,
				  null); // prohibit

		  bgnd_data.selected = bgnd_tiles_grown; // selected for background w/o extra transparent layer
		  
		  if (sdfa_instance!=null){
			  double [][] dbg_img = new double[3][tilesY * tilesX];
			  String [] titles = {"strict","grown","more_grown"};
			  for (int i = 0; i<dbg_img[0].length;i++){
				  dbg_img[0][i] = bgnd_strict[i]?1:0;
				  dbg_img[1][i] =  bgnd_tiles_grown[i]?1:0;
				  dbg_img[2][i] =  bgnd_tiles[i]?1:0;
			  }
			sdfa_instance.showArrays(dbg_img,  tilesX, tilesY, true, "strict_grown",titles);
		  }
		  
		  
		  
		  double [][][][] texture_tiles_bgnd = new double[tilesY][tilesX][][];
		  double [] alpha_zero = new double [4*clt_parameters.transform_size*clt_parameters.transform_size];
		  int alpha_index = 3;
		  for (int i = 0; i < alpha_zero.length; i++) alpha_zero[i]=0.0;
		  for (int tileY = 0; tileY < tilesY; tileY++){
			  for (int tileX = 0; tileX < tilesX; tileX++){
				  texture_tiles_bgnd[tileY][tileX]= null;
				  if ((texture_tiles[tileY][tileX] != null) &&
						  bgnd_tiles[tileY * tilesX + tileX]) {
					  if (bgnd_tiles_grown[tileY * tilesX + tileX]) {
						  texture_tiles_bgnd[tileY][tileX]= texture_tiles[tileY][tileX];
					  }else{
						  texture_tiles_bgnd[tileY][tileX]= texture_tiles[tileY][tileX].clone();
						  texture_tiles_bgnd[tileY][tileX][alpha_index] = alpha_zero;
					  }
					  
					  
				  }
			  }
		  }
		  
		  ImageDtt image_dtt = new ImageDtt();
		  double [][] texture_overlap = image_dtt.combineRGBATiles(
				  texture_tiles_bgnd, // texture_tiles,               // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}   
				  clt_parameters.transform_size,
				  true,                         // when false - output each tile as 16x16, true - overlap to make 8x8
				  clt_parameters.sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only 
				  threadsMax,                    // maximal number of threads to launch                         
				  debugLevel);
		  if (clt_parameters.alpha1 > 0){ // negative or 0 - keep alpha as it was
			  double scale = (clt_parameters.alpha1 > clt_parameters.alpha0) ? (1.0/(clt_parameters.alpha1 - clt_parameters.alpha0)) : 0.0; 
			  for (int i = 0; i < texture_overlap[alpha_index].length; i++){
				  double d = texture_overlap[alpha_index][i];
				  if      (d >=clt_parameters.alpha1) d = 1.0;
				  else if (d <=clt_parameters.alpha0) d = 0.0;
				  else d = scale * (d- clt_parameters.alpha0);
				  texture_overlap[alpha_index][i] = d;
			  }
		  }
		  // for now - use just RGB. Later add oprion for RGBA
		  double [][] texture_rgb = {texture_overlap[0],texture_overlap[1],texture_overlap[2]};
		  double [][] texture_rgba = {texture_overlap[0],texture_overlap[1],texture_overlap[2],texture_overlap[3]};
		  
		  //			  ImagePlus img_texture = 
		  ImagePlus imp_texture_bgnd = linearStackToColor(
				  clt_parameters,
				  colorProcParameters,
				  rgbParameters,
				  name+"-texture-bgnd", // String name,
				  "", //String suffix, // such as disparity=...
				  true, // toRGB,
				  !this.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
				  true, // boolean saveShowIntermediate, // save/show if set globally
				  false, //true, // boolean saveShowFinal,        // save/show result (color image?)
				  ((clt_parameters.alpha1 > 0)? texture_rgba: texture_rgb),
				  tilesX *  clt_parameters.transform_size,
				  tilesY *  clt_parameters.transform_size,
				  1.0,         // double scaleExposure, // is it needed?
				  debugLevel);
		  // resize for backdrop here!
//	public double getFOVPix(){ // get ratio of 1 pixel X/Y to Z (distance to object)
		  ImagePlus imp_texture_bgnd_ext = resizeForBackdrop(
				  imp_texture_bgnd,
				  debugLevel); 
		  
		  
		  String path= correctionsParameters.selectX3dDirectory(
				  true,  // smart,
				  true);  //newAllowed, // save
		  // only show/save original size if debug or debug_filters)  
		  if (clt_parameters.debug_filters || (debugLevel > 0)) {
			  eyesisCorrections.saveAndShow(
					  imp_texture_bgnd,
					  path,
					  correctionsParameters.png,
					  clt_parameters.show_textures,
					  -1); // jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
		  }
		  eyesisCorrections.saveAndShow(
				  imp_texture_bgnd_ext,
				  path,
				  correctionsParameters.png,
				  clt_parameters.show_textures,
				  -1); // jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
		  return imp_texture_bgnd_ext;
	  }
	  
	  
	  
//	  public ImagePlus getPassImage( // get image form a single pass
	 public String getPassImage( // get image form a single pass, return relative path for x3d
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  EyesisCorrectionParameters.ColorProcParameters colorProcParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  String     name,
			  int        scanIndex, 
			  int        threadsMax,  // maximal number of threads to launch                         
			  boolean    updateStatus,
			  int        debugLevel)
	  {
		 final int tilesX = tp.getTilesX();
		 final int tilesY = tp.getTilesY();
		 
		  showDoubleFloatArrays sdfa_instance = null;
//    	  if (clt_parameters.debug_filters && (debugLevel > -1)) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?		  
    	  if (debugLevel > -1) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?		  
		  TileProcessor.CLTPass3d scan = tp.clt_3d_passes.get(scanIndex);
		  boolean [] borderTiles = scan.border_tiles;
		  double [][][][] texture_tiles = scan.texture_tiles;
		  scan.updateSelection(); // update .selected field (all selected, including border) and Rectangle bounds
		  double [][]alphaFade = tp.getAlphaFade(clt_parameters.transform_size);
		  if ((debugLevel > 0) && (scanIndex == 1)) {
			  String [] titles = new String[16];
			  for (int i = 0; i<titles.length;i++)  titles[i]=""+i;
			  sdfa_instance.showArrays(alphaFade, 2*clt_parameters.transform_size,2*clt_parameters.transform_size,true,"alphaFade",titles);
		  }
		  double [][][][] texture_tiles_cluster = new double[tilesY][tilesX][][];
		  double [] alpha_zero = new double [4*clt_parameters.transform_size*clt_parameters.transform_size];
		  int alpha_index = 3;
		  for (int i = 0; i < alpha_zero.length; i++) alpha_zero[i]=0.0;
		  for (int tileY = 0; tileY < tilesY; tileY++){
			  for (int tileX = 0; tileX < tilesX; tileX++){
				  texture_tiles_cluster[tileY][tileX]= null;
				  if (texture_tiles[tileY][tileX] != null) {
					  if (borderTiles[tileY * tilesX + tileX]) {
						  texture_tiles_cluster[tileY][tileX]= texture_tiles[tileY][tileX].clone();
						  if (clt_parameters.shAggrFade) {
							  texture_tiles_cluster[tileY][tileX][alpha_index] = alpha_zero;
						  } else {
							  if ((debugLevel > -1) && (scanIndex == 1)) {
								  System.out.println("getPassImage(): tileY="+tileY+", tileX = "+tileX+", tileY="+tileY);
							  }
							  int fade_mode=0;
							  if ((tileY > 0) &&              (texture_tiles[tileY - 1][tileX] != null) && !borderTiles[(tileY - 1) * tilesX + tileX]) fade_mode |= 1;
							  if ((tileX < (tilesX -1)) && (texture_tiles[tileY][tileX + 1] != null) && !borderTiles[tileY * tilesX + tileX + 1])   fade_mode |= 2;
							  if ((tileY < (tilesY -1)) && (texture_tiles[tileY + 1][tileX] != null) && !borderTiles[(tileY + 1) * tilesX + tileX]) fade_mode |= 4;
							  if ((tileX > 0) &&              (texture_tiles[tileY][tileX - 1] != null) && !borderTiles[tileY * tilesX + tileX - 1])   fade_mode |= 8;
							  texture_tiles_cluster[tileY][tileX][alpha_index] = alphaFade[fade_mode]; // alpha_zero;
						  }
					  }else{
						  texture_tiles_cluster[tileY][tileX]= texture_tiles[tileY][tileX];
					  }
				  }
			  }
		  }
		  
		  ImageDtt image_dtt = new ImageDtt();
		  double [][] texture_overlap = image_dtt.combineRGBATiles(
				  texture_tiles_cluster, // texture_tiles,               // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}   
				  clt_parameters.transform_size,
				  true,                         // when false - output each tile as 16x16, true - overlap to make 8x8
				  clt_parameters.sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only 
				  threadsMax,                    // maximal number of threads to launch                         
				  debugLevel);
		  if (clt_parameters.alpha1 > 0){ // negative or 0 - keep alpha as it was
			  double scale = (clt_parameters.alpha1 > clt_parameters.alpha0) ? (1.0/(clt_parameters.alpha1 - clt_parameters.alpha0)) : 0.0; 
			  for (int i = 0; i < texture_overlap[alpha_index].length; i++){
				  double d = texture_overlap[alpha_index][i];
				  if      (d >=clt_parameters.alpha1) d = 1.0;
				  else if (d <=clt_parameters.alpha0) d = 0.0;
				  else d = scale * (d- clt_parameters.alpha0);
				  texture_overlap[alpha_index][i] = d;
			  }
		  }
		  // for now - use just RGB. Later add oprion for RGBA
		  double [][] texture_rgb = {texture_overlap[0],texture_overlap[1],texture_overlap[2]};
		  double [][] texture_rgba = {texture_overlap[0],texture_overlap[1],texture_overlap[2],texture_overlap[3]};
		  double [][] texture_rgbx = ((clt_parameters.alpha1 > 0)? texture_rgba: texture_rgb);

//    	  sdfa_instance = new showDoubleFloatArrays(); // just for debugging?		  
//     	  sdfa_instance.showArrays(texture_rgbx, clt_parameters.transform_size * tp.tilesX, clt_parameters.transform_size * tp.tilesY, true, "texture_rgbx_full");
		  
		  boolean resize = true;
		  if (resize) {
		  texture_rgbx = resizeGridTexture(
				  texture_rgbx,
				  clt_parameters.transform_size,
				  tilesX,
				  tilesY,
				  scan.bounds);
		  }
		  
		  
		  
//		  int width = resize ? (clt_parameters.transform_size * scan.bounds.width + 1): (clt_parameters.transform_size * tp.tilesX);
//		  int height = resize ? (clt_parameters.transform_size * scan.bounds.height + 1): (clt_parameters.transform_size * tp.tilesY);
		  int width = resize ? (clt_parameters.transform_size * scan.bounds.width): (clt_parameters.transform_size * tilesX);
		  int height = resize ? (clt_parameters.transform_size * scan.bounds.height): (clt_parameters.transform_size * tilesY);

//    	  sdfa_instance = new showDoubleFloatArrays(); // just for debugging?		  
//     	  sdfa_instance.showArrays(texture_rgbx, width, height, true, "texture_rgbx");
		  
		  ImagePlus imp_texture_cluster = linearStackToColor(
				  clt_parameters,
				  colorProcParameters,
				  rgbParameters,
				  name+"-texture", // String name,
				  "", //String suffix, // such as disparity=...
				  true, // toRGB,
				  !this.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
				  true, // boolean saveShowIntermediate, // save/show if set globally
				  false, //true, // boolean saveShowFinal,        // save/show result (color image?)
				  texture_rgbx, 
				  width, //tp.tilesX *  clt_parameters.transform_size,
				  height, //tp.tilesY *  clt_parameters.transform_size,
				  1.0,         // double scaleExposure, // is it needed?
				  debugLevel);
		  
		  
		  String path= correctionsParameters.selectX3dDirectory(
				  true,  // smart,
				  true);  //newAllowed, // save
		  // only show/save original size if debug or debug_filters)  
			  eyesisCorrections.saveAndShow(
					  imp_texture_cluster,
					  path,
					  correctionsParameters.png,
					  clt_parameters.show_textures,
					  -1); // jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
		  return imp_texture_cluster.getTitle()+".png"; // imp_texture_cluster;
	  }

	  
	  
	  public ImagePlus resizeForBackdrop(ImagePlus imp, int debugLevel){
		  double backdropPixels = 2.0/geometryCorrection.getFOVPix();
		  if (debugLevel > -1) {
			  System.out.println("backdropPixels = "+backdropPixels);
		  }
		  // TODO: currently - just adding pixels, no rescaling (add later). Alternatively - just modify geometry earlier
		  int width = imp.getWidth();
		  int height = imp.getHeight(); 
		  int h_margin = (int) Math.round((backdropPixels - width)/2);
		  int v_margin = (int) Math.round((backdropPixels - height)/2);
		  int width2 =  width +  2 * h_margin;
		  int height2 = height + 2 * v_margin;
		  if (debugLevel > -1) {
			  System.out.println("backdropPixels = "+backdropPixels+" h_margin = "+h_margin+" v_margin = "+v_margin);
		  }
		  int [] src_pixels = (int []) imp.getProcessor().getPixels();
		  int [] pixels = new int [width2* height2];
		  int indx = 0;
		  int offset = v_margin *  width2 + h_margin;
		  for (int i = 0; i < height; i++){
			  for (int j = 0; j < width; j++){
				  pixels[offset+ i * width2 + j] = src_pixels[indx++]; 
			  }
		  }
		  ColorProcessor cp=new ColorProcessor(width2,height2);
		  cp.setPixels(pixels);
		  ImagePlus imp_ext=new ImagePlus(imp.getTitle()+"-ext",cp);
		  return imp_ext;
	  }
	  
//[tp.tilesY][tp.tilesX]["RGBA".length()][]	  
//linearStackToColor

	  public TileProcessor.CLTPass3d CLTBackgroundMeas( // measure background
			  final double [][][]       image_data, // first index - number of image in a quad
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  final int         threadsMax,  // maximal number of threads to launch                         
			  final boolean     updateStatus,
			  final int         debugLevel)
	  {
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  TileProcessor.CLTPass3d scan_rslt = tp.new CLTPass3d();
		  int d = ImageDtt.setImgMask(0, 0xf);
		  d =     ImageDtt.setPairMask(d,0xf);
		  d =     ImageDtt.setForcedDisparity(d,true);
		  int [][]     tile_op =         tp.setSameTileOp(clt_parameters,  d, debugLevel);
		  double [][]  disparity_array = tp.setSameDisparity(0.0); // [tp.tilesY][tp.tilesX] - individual per-tile expected disparity
		  // undecided, so 2 modes of combining alpha - same as rgb, or use center tile only
		  double [][][][]     clt_corr_combo =    new double [ImageDtt.TCORR_TITLES.length][tilesY][tilesX][]; // will only be used inside?
		  
//		  double min_corr_selected = clt_parameters.corr_normalize? clt_parameters.min_corr_normalized: clt_parameters.min_corr;
		  double min_corr_selected = clt_parameters.min_corr;
		  
		  double [][] disparity_map = new double [ImageDtt.DISPARITY_TITLES.length][]; //[0] -residual disparity, [1] - orthogonal (just for debugging)

		  double [][] shiftXY = {
				  {clt_parameters.fine_corr_x_0,clt_parameters.fine_corr_y_0},
				  {clt_parameters.fine_corr_x_1,clt_parameters.fine_corr_y_1},
				  {clt_parameters.fine_corr_x_2,clt_parameters.fine_corr_y_2},
				  {clt_parameters.fine_corr_x_3,clt_parameters.fine_corr_y_3}};
		  
		  double [][][][] texture_tiles =     new double [tilesY][tilesX][][]; // ["RGBA".length()][];
		  ImageDtt image_dtt = new ImageDtt();
		  image_dtt.clt_aberrations_quad_corr(
				  tile_op,                      // per-tile operation bit codes
				  disparity_array,              // clt_parameters.disparity,     // final double            disparity,
				  image_data,                   // final double [][][]      imade_data, // first index - number of image in a quad
				  // correlation results - final and partial          
				  clt_corr_combo,               // [tp.tilesY][tp.tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				  null, // clt_corr_partial,    // [tp.tilesY][tp.tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				  null,    // [tp.tilesY][tp.tilesX][pair]{dx,dy,weight}[(2*transform_size-1)*(2*transform_size-1)] // transpose unapplied. null - do not calculate
				  //	Use it with disparity_maps[scan_step]?		  clt_mismatch,    // [tp.tilesY][tp.tilesX][pair]{dx,dy,weight}[(2*transform_size-1)*(2*transform_size-1)] // transpose unapplied. null - do not calculate
				  disparity_map,    // [12][tp.tilesY * tp.tilesX]
				  texture_tiles,        // [tp.tilesY][tp.tilesX]["RGBA".length()][]; 			  
				  tilesX * clt_parameters.transform_size, // imp_quad[0].getWidth(),       // final int width,
				  clt_parameters.fat_zero,      // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
				  clt_parameters.corr_sym,
				  clt_parameters.corr_offset,
				  clt_parameters.corr_red,
				  clt_parameters.corr_blue,
				  clt_parameters.corr_sigma,
				  clt_parameters.corr_normalize, // normalize correlation results by rms 
				  min_corr_selected, // 0.0001; // minimal correlation value to consider valid 
				  clt_parameters.max_corr_sigma,// 1.5;  // weights of points around global max to find fractional
				  clt_parameters.max_corr_radius,
				  clt_parameters.enhortho_width,  // 2;    // reduce weight of center correlation pixels from center (0 - none, 1 - center, 2 +/-1 from center)
				  clt_parameters.enhortho_scale,  // 0.2;  // multiply center correlation pixels (inside enhortho_width)
				  clt_parameters.max_corr_double, // Double pass when masking center of mass to reduce preference for integer values
				  clt_parameters.corr_mode,     // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
				  clt_parameters.min_shot,       // 10.0;  // Do not adjust for shot noise if lower than
				  clt_parameters.scale_shot,     // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
				  clt_parameters.diff_sigma,     // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_threshold, // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_gauss,     // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
				  clt_parameters.min_agree,      // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				  clt_parameters.dust_remove,    // Do not reduce average weight when only one image differes much from the average
				  clt_parameters.keep_weights,   // Add port weights to RGBA stack (debug feature)
				  geometryCorrection,           // final GeometryCorrection  geometryCorrection,
				  clt_kernels,                  // final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				  clt_parameters.kernel_step,
				  clt_parameters.transform_size,
				  clt_parameters.clt_window,
				  shiftXY, // 
				  (clt_parameters.fcorr_ignore? null: this.fine_corr),
				  clt_parameters.corr_magic_scale, // stil not understood coefficent that reduces reported disparity value.  Seems to be around 8.5 
				  clt_parameters.shift_x,       // final int               shiftX, // shift image horizontally (positive - right) - just for testing
				  clt_parameters.shift_y,       // final int               shiftY, // shift image vertically (positive - down)
				  clt_parameters.tileX,         // final int               debug_tileX,
				  clt_parameters.tileY,         // final int               debug_tileY,
				  (clt_parameters.dbg_mode & 64) != 0, // no fract shift
				  (clt_parameters.dbg_mode & 128) != 0, // no convolve
				  //				  (clt_parameters.dbg_mode & 256) != 0, // transpose convolve
				  threadsMax,
				  debugLevel);
		  
		  scan_rslt.disparity = disparity_array;
		  scan_rslt.tile_op = tile_op;
		  scan_rslt.disparity_map = disparity_map;
		  scan_rslt.texture_tiles = texture_tiles;
		  return scan_rslt;
	  }

	  public TileProcessor.CLTPass3d  CLTMeasure( // perform single pass according to prepared tiles operations and disparity
			  final double [][][]       image_data, // first index - number of image in a quad
			  EyesisCorrectionParameters.CLTParameters           clt_parameters,
			  final int         scanIndex,
			  final int         threadsMax,  // maximal number of threads to launch                         
			  final boolean     updateStatus,
			  final int         debugLevel)
	  {
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  TileProcessor.CLTPass3d scan = tp.clt_3d_passes.get(scanIndex);
		  int [][]     tile_op =         scan.tile_op;
		  double [][]  disparity_array = scan.disparity;
		  // undecided, so 2 modes of combining alpha - same as rgb, or use center tile only
		  double [][][][]     clt_corr_combo =    new double [ImageDtt.TCORR_TITLES.length][tilesY][tilesX][]; // will only be used inside?
		  
		  double min_corr_selected = clt_parameters.min_corr;
		  
		  double [][] disparity_map = new double [ImageDtt.DISPARITY_TITLES.length][]; //[0] -residual disparity, [1] - orthogonal (just for debugging)

		  double [][] shiftXY = {
				  {clt_parameters.fine_corr_x_0,clt_parameters.fine_corr_y_0},
				  {clt_parameters.fine_corr_x_1,clt_parameters.fine_corr_y_1},
				  {clt_parameters.fine_corr_x_2,clt_parameters.fine_corr_y_2},
				  {clt_parameters.fine_corr_x_3,clt_parameters.fine_corr_y_3}};
		  
		  double [][][][] texture_tiles =     new double [tilesY][tilesX][][]; // ["RGBA".length()][];
		  ImageDtt image_dtt = new ImageDtt();
		  image_dtt.clt_aberrations_quad_corr(
				  tile_op,                      // per-tile operation bit codes
				  disparity_array,              // clt_parameters.disparity,     // final double            disparity,
				  image_data,                   // final double [][][]      imade_data, // first index - number of image in a quad
				  // correlation results - final and partial          
				  clt_corr_combo,               // [tp.tilesY][tp.tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				  null, // clt_corr_partial,    // [tp.tilesY][tp.tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				  null,    // [tp.tilesY][tp.tilesX][pair]{dx,dy,weight}[(2*transform_size-1)*(2*transform_size-1)] // transpose unapplied. null - do not calculate
				  //	Use it with disparity_maps[scan_step]?		  clt_mismatch,    // [tp.tilesY][tp.tilesX][pair]{dx,dy,weight}[(2*transform_size-1)*(2*transform_size-1)] // transpose unapplied. null - do not calculate
				  disparity_map,    // [12][tp.tilesY * tp.tilesX]
				  texture_tiles,        // [tp.tilesY][tp.tilesX]["RGBA".length()][]; 			  
				  tilesX * clt_parameters.transform_size, // imp_quad[0].getWidth(),       // final int width,
				  clt_parameters.fat_zero,      // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
				  clt_parameters.corr_sym,
				  clt_parameters.corr_offset,
				  clt_parameters.corr_red,
				  clt_parameters.corr_blue,
				  clt_parameters.corr_sigma,
				  clt_parameters.corr_normalize, // normalize correlation results by rms 
				  min_corr_selected, // 0.0001; // minimal correlation value to consider valid 
				  clt_parameters.max_corr_sigma,// 1.5;  // weights of points around global max to find fractional
				  clt_parameters.max_corr_radius,
				  clt_parameters.enhortho_width,  // 2;    // reduce weight of center correlation pixels from center (0 - none, 1 - center, 2 +/-1 from center)
				  clt_parameters.enhortho_scale,  // 0.2;  // multiply center correlation pixels (inside enhortho_width)
				  clt_parameters.max_corr_double, // Double pass when masking center of mass to reduce preference for integer values
				  clt_parameters.corr_mode,     // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
				  clt_parameters.min_shot,       // 10.0;  // Do not adjust for shot noise if lower than
				  clt_parameters.scale_shot,     // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
				  clt_parameters.diff_sigma,     // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_threshold, // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_gauss,     // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
				  clt_parameters.min_agree,      // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				  clt_parameters.dust_remove,    // Do not reduce average weight when only one image differes much from the average
				  clt_parameters.keep_weights,   // Add port weights to RGBA stack (debug feature)
				  geometryCorrection,           // final GeometryCorrection  geometryCorrection,
				  clt_kernels,                  // final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				  clt_parameters.kernel_step,
				  clt_parameters.transform_size,
				  clt_parameters.clt_window,
				  shiftXY, // 
				  (clt_parameters.fcorr_ignore? null: this.fine_corr),
				  clt_parameters.corr_magic_scale, // stil not understood coefficent that reduces reported disparity value.  Seems to be around 8.5 
				  clt_parameters.shift_x,       // final int               shiftX, // shift image horizontally (positive - right) - just for testing
				  clt_parameters.shift_y,       // final int               shiftY, // shift image vertically (positive - down)
				  clt_parameters.tileX,         // final int               debug_tileX,
				  clt_parameters.tileY,         // final int               debug_tileY,
				  (clt_parameters.dbg_mode & 64) != 0, // no fract shift
				  (clt_parameters.dbg_mode & 128) != 0, // no convolve
				  //				  (clt_parameters.dbg_mode & 256) != 0, // transpose convolve
				  threadsMax,
				  debugLevel);
		  
		  scan.disparity_map = disparity_map;
		  scan.texture_tiles = texture_tiles;
		  return scan;
	  }
	  
	  
	  
}
