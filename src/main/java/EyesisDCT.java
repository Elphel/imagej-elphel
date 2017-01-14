/**
 **
 ** EyesisDCT - Process images with DTT-based methods (code specific to ImageJ plugin)
 **
 ** Copyright (C) 2016 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  EyesisDCT.java is free software: you can redistribute it and/or modify
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

import java.util.concurrent.atomic.AtomicInteger;

import ij.CompositeImage;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.io.FileSaver;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;


public class EyesisDCT {
	public EyesisCorrections eyesisCorrections = null;
	public EyesisCorrectionParameters.CorrectionParameters correctionsParameters=null;
	public EyesisCorrectionParameters.DCTParameters dctParameters = null;
	public DCTKernels [] kernels = null;
	public ImagePlus eyesisKernelImage = null;
	public long startTime;
	
	public EyesisDCT(
			EyesisCorrections eyesisCorrections,
			EyesisCorrectionParameters.CorrectionParameters correctionsParameters,
			EyesisCorrectionParameters.DCTParameters dctParameters
			){
		this.eyesisCorrections= eyesisCorrections;
		this.correctionsParameters = correctionsParameters;
		this.dctParameters= dctParameters;
	}
	public class DCTKernels{
		// -- old --
		public int           size = 32;     // kernel (DCT) size
		public int           img_step = 32; // pixel step in the image for each kernel
		public double [][][] offsets = null; // per color, per kernel,per coord
		public double [][]   kern = null;  // kernel image in linescan order
		// -- new --
		public int           numHor =        164; // number of kernel tiles in a row
		public int           dct_size =        8;  // DCT-II size, sym. kernel square side is 2*dct_size-1 
		public int           asym_size =      15; // asymmetrical convolution limits, odd
		public int           asym_nonzero =   4; // maximal number of non-zero elements in the asymmetrical kernels 
		public double [][]   sym_kernels =  null; // per-color channel, DCT kernels in linescan order
		public double [][]   sym_direct =   null; // per-color channel, DCT kernels in linescan order (direct, not dct-iii transformed) - debug feature
		public double [][]   asym_kernels = null; // per-color channel, asymmetrical kernels (all but asym_nonzero elements are strictly 0)
		public double [][][][] st_kernels = null; // [color][tileY][tileX][pixel]
		public double [][][][] st_direct =  null; // [color][tileY][tileX][pixel] - direct, not converted with DCT-III - debug feature
		public double [][][][] asym_val =   null; // [color][tileY][tileX][index] // value - asym kernel for elements
		public int    [][][][] asym_indx =  null; // [color][tileY][tileX][index] // value - index of non-zero elements in the list 
	}

	public void setKernelImageFile(ImagePlus img_kernels){
		eyesisKernelImage = img_kernels;
	}
	
	public boolean kernelImageSet(){
		return eyesisKernelImage != null;
	}
	
	public boolean DCTKernelsAvailable(){
		return kernels != null;
	}
	
	public boolean createDCTKernels(
			EyesisCorrectionParameters.DCTParameters dct_parameters,
/*			
			PixelMapping pixelMapping,
*/			
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
		if (kernels == null){
			kernels = new DCTKernels[eyesisCorrections.usedChannels.length];
			for (int chn=0;chn<eyesisCorrections.usedChannels.length;chn++){
				kernels[chn] = null;
			}
		}
	    showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
	    

		for (int chn=0;chn<eyesisCorrections.usedChannels.length;chn++){
			if (eyesisCorrections.usedChannels[chn] && (sharpKernelPaths[chn]!=null) && (kernels[chn]==null)){
				ImagePlus imp_kernel_sharp=new ImagePlus(sharpKernelPaths[chn]);
				if (imp_kernel_sharp.getStackSize()<3) {
					System.out.println("Need a 3-layer stack with kernels");
					sharpKernelPaths[chn]=null;
					continue;
				}
				ImageStack kernel_sharp_stack= imp_kernel_sharp.getStack();
				System.out.println("debugLevel = "+debugLevel+" kernel_sharp_stack.getWidth() = "+kernel_sharp_stack.getWidth()+
						" kernel_sharp_stack.getHeight() = "+kernel_sharp_stack.getHeight());
				DCTKernels kernels = calculateDCTKernel (
						kernel_sharp_stack,            // final ImageStack kernelStack,  // first stack with 3 colors/slices convolution kernels
						srcKernelSize,                 // final int          kernelSize, // 64
						dct_parameters,                // final double       blurSigma,
						threadsMax,  // maximal number of threads to launch                         
						updateStatus,
						debugLevel); // update status info
				int sym_width =  kernels.numHor * kernels.dct_size;
				int sym_height = kernels.sym_kernels[0].length /sym_width;
// save files
				String [] symNames = {"red_sym","blue_sym","green_sym"};
				String [] asymNames = {"red_asym","blue_asym","green_asym"};
				ImageStack symStack = sdfa_instance.makeStack(
						kernels.sym_kernels,
						sym_width,
						sym_height,
						symNames);
	    		String symPath=correctionsParameters.dctKernelDirectory+
	    				           Prefs.getFileSeparator()+
	    				           correctionsParameters.dctKernelPrefix+
	    				           String.format("%02d",chn)+
	    				           correctionsParameters.dctSymSuffix;
	    		String msg="Saving symmetrical convolution kernels to "+symPath;
	    		IJ.showStatus(msg);
	    		if (debugLevel>0) System.out.println(msg);
				ImagePlus imp_sym=new ImagePlus(imp_kernel_sharp.getTitle()+"-sym",symStack);
	    		if (debugLevel > 1) {
	    			imp_sym.getProcessor().resetMinAndMax();
	    			imp_sym.show();
	    		}
        		FileSaver fs=new FileSaver(imp_sym);
        		fs.saveAsTiffStack(symPath);
				
//				sdfa_instance.showArrays(kernels.sym_kernels,  sym_width, sym_height, true, imp_kernel_sharp.getTitle()+"-sym");

				int asym_width =  kernels.numHor * kernels.asym_size;
				int asym_height = kernels.asym_kernels[0].length /asym_width;
				ImageStack asymStack = sdfa_instance.makeStack(
						kernels.asym_kernels,
						asym_width,
						asym_height,
						asymNames);
	    		String asymPath=correctionsParameters.dctKernelDirectory+
	    				           Prefs.getFileSeparator()+
	    				           correctionsParameters.dctKernelPrefix+
	    				           String.format("%02d",chn)+
	    				           correctionsParameters.dctAsymSuffix;
	    		msg="Saving asymmetrical convolution kernels "+asymPath;
	    		IJ.showStatus(msg);
	    		if (debugLevel>0) System.out.println(msg);
				ImagePlus imp_asym=new ImagePlus(imp_kernel_sharp.getTitle()+"-asym",asymStack);
	    		if (debugLevel > 1) {
	    			imp_asym.getProcessor().resetMinAndMax();
	    			imp_asym.show();
	    		}
        		fs=new FileSaver(imp_asym);
        		fs.saveAsTiffStack(asymPath);
//				sdfa_instance.showArrays(kernels.asym_kernels,  asym_width, asym_height, true, imp_kernel_sharp.getTitle()+"-asym");
			}
		}
		return true;
	}
	

	  public DCTKernels calculateDCTKernel (
			  final ImageStack kernelStack,  // first stack with 3 colors/slices convolution kernels
			  final int          kernelSize, // 64
			  final EyesisCorrectionParameters.DCTParameters dct_parameters,
/*			  
			  final double [][]  vignetting,
			  int                vign_width,
			  int                vign_height,
			  int                vign_decimation,
*/			  
			  final int          threadsMax,  // maximal number of threads to launch                         
			  final boolean      updateStatus,
			  final int          globalDebugLevel) // update status info
	  {
		  if (kernelStack==null) return null;
		  final int kernelWidth=kernelStack.getWidth();
		  final int kernelNumHor=kernelWidth/kernelSize;
		  final int kernelNumVert=kernelStack.getHeight()/kernelSize;
		  final int nChn=kernelStack.getSize();
//		  final int length=kernelNumHor*kernelNumVert* dct_parameters.dct_size * dct_parameters.dct_size;// size of kernel data 
		  final DCTKernels dct_kernel = new DCTKernels();
		  dct_kernel.size = dct_parameters.dct_size;
		  dct_kernel.img_step = kernelSize/2/dct_parameters.decimation ; // May be wrong
		  dct_kernel.sym_kernels =  new double [nChn][kernelNumHor*kernelNumVert*dct_parameters.dct_size * dct_parameters.dct_size];
		  dct_kernel.asym_kernels = new double [nChn][kernelNumHor*kernelNumVert*dct_parameters.asym_size * dct_parameters.asym_size];
		  dct_kernel.asym_nonzero = dct_parameters.asym_pixels;
		  // currently each 64x64 kernel corresponds to 16x16 original pixels tile, 2 tiles margin each side
		  final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		  final AtomicInteger ai = new AtomicInteger(0);
		  final int numberOfKernels=     kernelNumHor*kernelNumVert*nChn;
		  final int numberOfKernelsInChn=kernelNumHor*kernelNumVert;
		  final int dct_size =      dct_parameters.dct_size; 
		  final int preTargetSize = 4 * dct_size; 
		  final int targetSize =    2 * dct_size; // normally 16
		  final double [] anitperiodic_window = createAntiperiodicWindow(dct_size);
/*		  
		  final int vgn_step =   dct_kernel.img_step/vign_decimation;
		  final int vgn_width =  vign_width/vign_decimation;
		  final int vgn_height = vign_height/vign_decimation;
		  final int vgn_left =   (vign_width  - (dct_kernel.img_step * (kernelNumHor-1)))/(2* vign_decimation); // in decimated pixels
		  final int vgn_top =    (vign_height  -(dct_kernel.img_step * (kernelNumVert-1)))/(2* vign_decimation); // in decimated pixels
		  final double vgn_max = dct_parameters.vignetting_max;
		  final double vgn_min = vgn_max/dct_parameters.vignetting_range;
*/		  
		  final int chn_green = 2; // all others multiply by 4 as there 1 in 4 Bayer for those, green - by 2
		  final long startTime = System.nanoTime();
		  System.out.println("calculateDCTKernel():numberOfKernels="+numberOfKernels);
		  for (int ithread = 0; ithread < threads.length; ithread++) {
			  threads[ithread] = new Thread() {
				  public void run() {
					  DoubleGaussianBlur gb=null;
					  if (dct_parameters.decimateSigma > 0)	 gb=new DoubleGaussianBlur();
					  float [] kernelPixels= null; // will be initialized at first use NOT yet?
					  double [] kernel=      new double[kernelSize*kernelSize];
					  
//					  int targetSize = dct_parameters.asym_size + 2 * dct_parameters.dct_size - 2;
					  
					  double [] pre_target_kernel= new double [preTargetSize * preTargetSize]; // before made antiperiodic 
					  double [] target_kernel = new double [targetSize * targetSize]; // strictly antiperiodic in both x and y directions
					  
					  FactorConvKernel factorConvKernel = new FactorConvKernel();
					  factorConvKernel.setDebugLevel       (0); // globalDebugLevel);
//					  factorConvKernel.setTargetWindowMode (dct_parameters.dbg_window_mode, dct_parameters.centerWindowToTarget);
					  factorConvKernel.setTargetWindowMode (dct_parameters.centerWindowToTarget);
					  factorConvKernel.numIterations =     dct_parameters.LMA_steps;
					  factorConvKernel.setAsymCompactness  (dct_parameters.compactness,	dct_parameters.asym_tax_free);
					  factorConvKernel.setSymCompactness   (dct_parameters.sym_compactness);
					  factorConvKernel.setDCWeight         (dct_parameters.dc_weight);
					  
					  int chn,tileY,tileX;
//					  int chn0=-1;
//					  int i;
//					  double sum;
					  for (int nTile = ai.getAndIncrement(); nTile < numberOfKernels; nTile = ai.getAndIncrement()) {
						  chn=nTile/numberOfKernelsInChn;
						  tileY =(nTile % numberOfKernelsInChn)/kernelNumHor;
						  tileX = nTile % kernelNumHor;
/*						  
						  if (vignetting != null){
							  int vh = vgn_left + vgn_step * tileX;
							  if (vh < 0)           vh = 0;
							  if (vh >= vgn_width)  vh = vgn_width -1;
							  int vv = vgn_top + vgn_step * tileY;
							  if (vv < 0)           vv = 0;
							  if (vv >= vgn_height) vh = vgn_height -1;
							  kscale = vignetting[chn][vv*vgn_width+vh];
							  if (kscale < vgn_min) kscale = vgn_min; 
						  }
						  kscale = vgn_max/kscale * ((chn == chn_green)? 2:4);
*/						  
////						  kscale = (chn == chn_green)? 2:4;
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
						  
						  if ((dct_parameters.decimation == 2) && (dct_parameters.decimateSigma<0)) {
							  reformatKernel2( // averages by exactly 2 (decimate==2)
									  kernel,
									  pre_target_kernel, // expand/crop, blur/decimate result (32x32)
									  kernelSize,
									  preTargetSize); // 32
						  } else {
							  reformatKernel(
									  kernel,        // will be blurred in-place
									  pre_target_kernel, // expand/crop, blur/decimate result (32x32)
									  kernelSize,
									  preTargetSize, // 32
									  dct_parameters.decimation,
									  dct_parameters.decimateSigma,
									  gb);
						  }
						  if (dct_parameters.normalize) { // or should it be normalized after antiperiodic?
							  double s =0.0;
							  for (int i = 0; i < pre_target_kernel.length; i++){
								  s+=pre_target_kernel[i];
							  }
							  s =  1.0 / s;
							  for (int i = 0; i < pre_target_kernel.length; i++){
								  pre_target_kernel[i] *= s;
							  }
							  if (globalDebugLevel > 1){ // was already close to 1.0
								  System.out.println(tileX+"/"+tileY+ " s="+s);
							  }
						  }
						  // make exactly anitperiodic
						  makeAntiperiodic(
								  dct_size,
								  pre_target_kernel,   // 16*dct_zize*dct_zize
								  anitperiodic_window, // 16*dct_zize*dct_zize
								  target_kernel);      //  4*dct_zize*dct_zize
						  
						  
						  factorConvKernel.calcKernels(
								  target_kernel,
								  dct_parameters.asym_size,
								  dct_parameters.dct_size,
								  dct_parameters.fact_precision,
								  dct_parameters.asym_pixels,         // maximal number of non-zero pixels in asymmmetrical kernel
								  dct_parameters.asym_distance,       // how far to seed a new pixel
								  dct_parameters.seed_size);
						  double [] sym_kernel =  factorConvKernel.getSymKernel(); 
						  double [] asym_kernel = factorConvKernel.getAsymKernel();
						  int sym_kernel_inc_index =   kernelNumHor * dct_parameters.dct_size;
						  int sym_kernel_start_index = (sym_kernel_inc_index * tileY + tileX) * dct_parameters.dct_size;
//						  int indx = 0;
						  for (int i = 0; i<dct_parameters.dct_size;i++){
							  
								System.arraycopy(
										sym_kernel,
										i * dct_parameters.dct_size,
										dct_kernel.sym_kernels[chn],
										sym_kernel_start_index + i * sym_kernel_inc_index,
										dct_parameters.dct_size);

							  /*
							  int dst_start = sym_kernel_start_index + i * sym_kernel_inc_index;
							  for (int j = 0; j < dct_parameters.dct_size; j++){
								  dct_kernel.sym_kernels[chn][dst_start++] = sym_kernel[indx++];
							  }
							  */
						  }
						  int asym_kernel_inc_index =   kernelNumHor * dct_parameters.asym_size;
						  int asym_kernel_start_index = (asym_kernel_inc_index * tileY + tileX)* dct_parameters.asym_size;
						  for (int i = 0; i<dct_parameters.asym_size;i++){
								System.arraycopy(
										asym_kernel,
										i * dct_parameters.asym_size,
										dct_kernel.asym_kernels[chn],
										asym_kernel_start_index + i * asym_kernel_inc_index,
										dct_parameters.asym_size);
						  }
					  }
				  }
			  };
		  }		      
		  ImageDtt.startAndJoin(threads);
		  if (globalDebugLevel > 1) System.out.println("Threads done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		  System.out.println("1.Threads done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		  /* prepare result stack to return */
//		  ImageStack outStack=new ImageStack(kernelNumHor,kernelNumVert);
		  return dct_kernel;
	  }
//processChannelImage	
//convolveStackWithKernelStack	

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
//		  final int kernelNumVert=kernelStack.getHeight()/kernelSize;
//		  final int nChn=kernelStack.getSize();
		  double [] kernel = new double [kernelSize*kernelSize];
		  extractOneKernel((float[]) kernelStack.getPixels(chn+1), //  array of combined square kernels, each 
				  kernel,        // will be filled, should have correct size before call
				  kernelNumHor,  // number of kernels in a row
				  xTile,         // horizontal number of kernel to extract
				  yTile);
		  return kernel;
	  }
	  
	  //imp2.getStack()
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
	  
	  private double [] createAntiperiodicWindow(
			  int dct_size)
	  {
		  double [] wnd = new double[4*dct_size];
		  for (int i =0; i<dct_size; i++){
			  wnd[i]=                 0.5 * (1.0-Math.cos(Math.PI*i/dct_size));
			  wnd[i + 1 * dct_size] = 1.0;
			  wnd[i + 2 * dct_size] = 1.0;
			  wnd[i + 3 * dct_size] = 1.0 - wnd[i];
		  }
		  int n4 = dct_size*4;
		  double [] window = new double [n4 * n4];
		  for (int i =0; i < n4; i++){
			  for (int j =0; j < n4; j++){
				  window[i * n4 + j] = wnd[i]*wnd[j]; 
			  }
		  }
		  return window;
	  }
	  
	  public double []makeAntiperiodic(
			  int       dct_size,
			  double [] src_kernel) // 16*dct_zize*dct_zize
	  {
		 double [] window =  createAntiperiodicWindow(dct_size);
		 double [] antiperiodic_kernel = new double [4*dct_size*dct_size];
		 makeAntiperiodic(
				  dct_size,
				  src_kernel, // 16*dct_zize*dct_zize
				  window,     // 16*dct_zize*dct_zize
				  antiperiodic_kernel); //  4*dct_zize*dct_zize
		 
		 return antiperiodic_kernel;
	  }

	  
	  private void makeAntiperiodic(
			  int       dct_size,
			  double [] src_kernel,          // 16*dct_zize*dct_zize
			  double [] window,              // 16*dct_zize*dct_zize
			  double [] antiperiodic_kernel) //  4*dct_zize*dct_zize
	  {
		  int n2 = dct_size * 2;
		  int n4 = dct_size * 4;
		  for (int i = 0; i < n2; i++){
			  for (int j = 0; j < n2; j++){
				  int dst_index = i*n2+j;
				  int isrc = i+dct_size;
				  int jsrc =  j+dct_size;
				  int isrcp = (isrc + n2) % n4;
				  int jsrcp = (jsrc + n2) % n4;
				  int src_index0= (isrc)  * n4 + jsrc;
				  int src_index1= (isrcp) * n4 + jsrc;
				  int src_index2= (isrc) *  n4 + jsrcp;
				  int src_index3= (isrcp) * n4 + jsrcp;
				  
				  antiperiodic_kernel[dst_index] =
						    src_kernel[src_index0] * window[src_index0]
						  - src_kernel[src_index1] * window[src_index1]
						  - src_kernel[src_index2] * window[src_index2]
						  + src_kernel[src_index3] * window[src_index3];				  
			  }
		  }
	  }
	  
		public boolean readDCTKernels(
				EyesisCorrectionParameters.DCTParameters dct_parameters,
				int          srcKernelSize,
				int          threadsMax,  // maximal number of threads to launch                         
				boolean      updateStatus,
				int          debugLevel
				){
			String [] symKernelPaths = correctionsParameters.selectDCTChannelFiles(
//					0,  // 0 - sharp, 1 - smooth
					eyesisCorrections.usedChannels.length, // numChannels, // number of channels
					eyesisCorrections.debugLevel);
			if (symKernelPaths==null) return false;
			
			String [] asymKernelPaths = new String[symKernelPaths.length];
			for (int chn = 0; chn <symKernelPaths.length; chn++ ) if (symKernelPaths[chn] != null){
				int indexPeriod=symKernelPaths[chn].indexOf('.',symKernelPaths[chn].lastIndexOf(Prefs.getFileSeparator()));
				asymKernelPaths[chn] = symKernelPaths[chn].substring(0, indexPeriod) + correctionsParameters.dctAsymSuffix;
			}
			
			for (int i=0;i<symKernelPaths.length;i++){
				System.out.println(i+":"+symKernelPaths[i]+", "+asymKernelPaths[i]); // some may be null!
			}
			if (kernels == null){
				kernels = new DCTKernels[eyesisCorrections.usedChannels.length];
				for (int chn=0;chn<eyesisCorrections.usedChannels.length;chn++){
					kernels[chn] = null;
				}
			}
		    showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
			DttRad2 dtt = new DttRad2(dct_parameters.dct_size);

			for (int chn=0;chn<eyesisCorrections.usedChannels.length;chn++){
//				if (eyesisCorrections.usedChannels[chn] && (symKernelPaths[chn]!=null) && (kernels[chn]==null)){
				if (eyesisCorrections.usedChannels[chn] && (symKernelPaths[chn]!=null)){
					ImagePlus imp_kernel_sym=new ImagePlus(symKernelPaths[chn]);
					if (imp_kernel_sym.getStackSize()<3) {
						System.out.println("Need a 3-layer stack with symmetrical DCT kernels");
						symKernelPaths[chn]=null;
						continue;
					}
					ImagePlus imp_kernel_asym=new ImagePlus(asymKernelPaths[chn]);
					if (imp_kernel_sym.getStackSize()<3) {
						System.out.println("Need a 3-layer stack with asymmetrical kernels");
						asymKernelPaths[chn]=null;
						continue;
					}
					
					ImageStack kernel_sym_stack=  imp_kernel_sym.getStack();
					ImageStack kernel_asym_stack= imp_kernel_asym.getStack();
					if (debugLevel>0){
						System.out.println(" kernel_asym_stack.getWidth() = "+kernel_asym_stack.getWidth()+
								" kernel_asym_stack.getHeight() = "+kernel_asym_stack.getHeight());
					}
					int nColors = kernel_sym_stack.getSize();
					kernels[chn]=  new DCTKernels();
					kernels[chn].size = dct_parameters.dct_size;
					kernels[chn].img_step = srcKernelSize/2/dct_parameters.decimation ; // May be wrong
					kernels[chn].asym_nonzero = dct_parameters.asym_pixels;
					
					kernels[chn].sym_kernels =  new double [nColors][];
					kernels[chn].asym_kernels = new double [nColors][];
					
					for (int nc = 0; nc < nColors; nc++){
						float [] pixels = (float[]) kernel_sym_stack.getPixels(nc + 1);
						kernels[chn].sym_kernels[nc]= new double[pixels.length]; 
						for (int i = 0; i<pixels.length; i++){
							kernels[chn].sym_kernels[nc][i] = pixels[i];
						}
						pixels = (float[]) kernel_asym_stack.getPixels(nc + 1);
						kernels[chn].asym_kernels[nc]= new double[pixels.length]; 
						for (int i = 0; i<pixels.length; i++){
							kernels[chn].asym_kernels[nc][i] = pixels[i];
						}
					}
					int dct_size = kernels[chn].dct_size;
					int asym_size= kernels[chn].asym_size;
					int sym_width =  kernels[chn].numHor * dct_size;
					int sym_height = kernels[chn].sym_kernels[0].length /sym_width;
					int asym_nonzero =kernels[chn].asym_nonzero;
//					sdfa_instance.showArrays(kernels[chn].sym_kernels,  sym_width, sym_height, true, symKernelPaths[chn]);

					int asym_width =  kernels[chn].numHor * kernels[chn].asym_size;
					int asym_height = kernels[chn].asym_kernels[0].length /asym_width;
//					sdfa_instance.showArrays(kernels[chn].asym_kernels,  asym_width, asym_height, true, asymKernelPaths[chn]);
					int numHor =   kernels[chn].numHor;
					int numVert =  kernels[chn].sym_kernels[0].length / (dct_size * dct_size * numHor);
					kernels[chn].st_kernels = new double [nColors][numVert][numHor][dct_size * dct_size];
					kernels[chn].st_direct =  new double [nColors][numVert][numHor][dct_size * dct_size];
					kernels[chn].asym_val =   new double [nColors][numVert][numHor][asym_nonzero];
					kernels[chn].asym_indx =  new int [nColors][numVert][numHor][asym_nonzero];
					int sym_kernel_inc_index =   numHor * dct_size;
					int asym_kernel_inc_index =   numHor * asym_size;
					if (debugLevel>0) {
						System.out.println("readDCTKernels() debugLevel = "+debugLevel+
								" kernels["+chn+"].size = "+kernels[chn].size+
								" kernels["+chn+"].img_step = "+kernels[chn].img_step+
								" kernels["+chn+"].asym_nonzero = "+kernels[chn].asym_nonzero+
								" nColors = "+ nColors+
								" numVert = "+ numVert+
								" numHor =  "+ numHor);
					}
					for (int nc = 0; nc < nColors; nc++){
						for (int tileY = 0; tileY < numVert; tileY++){
							for (int tileX = 0; tileX < numHor; tileX++){
								// extract asymmetrical kernel and convert it to list of values and indices (as arrays as the length is known)
								int asym_kernel_start_index = (asym_kernel_inc_index * tileY + tileX)* asym_size;
								int indx = 0;
								for (int i = 0; (i < dct_parameters.asym_size) && (indx < asym_nonzero); i++){
									for (int j = 0; (j < dct_parameters.asym_size) && (indx < asym_nonzero); j++){
										double v = kernels[chn].asym_kernels[nc][asym_kernel_start_index + i * asym_kernel_inc_index +j];
										if (v!=0.0){
											if ((debugLevel>0) && (tileY==67) && (tileX==125)) {
												System.out.println("i="+i+" j="+j+" v="+v+" indx="+indx+" i * asym_size + j="+(i * asym_size + j));
												System.out.println("asym_kernel_start_index + i * asym_kernel_inc_index +j="+(asym_kernel_start_index + i * asym_kernel_inc_index +j));
											}
										
											kernels[chn].asym_val[nc][tileY][tileX][indx] =    v;
											kernels[chn].asym_indx[nc][tileY][tileX][indx++] = i * asym_size + j;
										}
									}
								}
								if ((debugLevel>0) && (tileY==67) && (tileX==125)) {
									for (int i=0; i<kernels[chn].asym_indx[nc][tileY][tileX].length; i++){
										System.out.println("kernels["+chn+"].asym_val["+nc+"]["+tileY+"]["+tileX+"]["+i+"]="+kernels[chn].asym_val[nc][tileY][tileX][i]);
										System.out.println("kernels["+chn+"].asym_indx["+nc+"]["+tileY+"]["+tileX+"]["+i+"]="+kernels[chn].asym_indx[nc][tileY][tileX][i]);
									}
								}
								
								double scale_asym = 0.0;
								if (dct_parameters.normalize) {
									for (int i = 0; i < kernels[chn].asym_val[nc][tileY][tileX].length;i++){
										scale_asym += kernels[chn].asym_val[nc][tileY][tileX][i];
										if ((debugLevel>0) && (tileY==67) && (tileX==125)) {
											System.out.println("i="+i+", sum="+scale_asym);
										}
									}
									double k = ((nc == 2)? 1.0:2.0) / scale_asym;
									for (int i = 0; i < kernels[chn].asym_val[nc][tileY][tileX].length;i++){
//										kernels[chn].asym_val[nc][tileY][tileX][i] /= scale_asym;
										kernels[chn].asym_val[nc][tileY][tileX][i] *= k; // includes correction for different number of pixels in r,b(1/4) and G (2/4)
									}
									if ((debugLevel>0) && (tileY==67) && (tileX==125)) {
										System.out.println("nc="+nc+" sum="+scale_asym+", normalized:");

										for (int i=0; i<kernels[chn].asym_indx[nc][tileY][tileX].length; i++){
											System.out.println("kernels["+chn+"].asym_val["+nc+"]["+tileY+"]["+tileX+"]["+i+"]="+kernels[chn].asym_val[nc][tileY][tileX][i]);
											System.out.println("kernels["+chn+"].asym_indx["+nc+"]["+tileY+"]["+tileX+"]["+i+"]="+kernels[chn].asym_indx[nc][tileY][tileX][i]);
										}
									}
								} else {
									scale_asym = 1.0;
								}

								// extract DCT (symmetrical) kernels
								int sym_kernel_start_index = (sym_kernel_inc_index * tileY + tileX) * dct_size;
								for (int i = 0; i < dct_size;i++){
									System.arraycopy( // copy one kernel line
											kernels[chn].sym_kernels[nc],
											sym_kernel_start_index + i * sym_kernel_inc_index,
											kernels[chn].st_kernels[nc][tileY][tileX],
											i * dct_size,
											dct_size);
								}
								if (scale_asym != 1.0){
									for (int i=0; i < kernels[chn].st_kernels[nc][tileY][tileX].length;i++) {
										kernels[chn].st_kernels[nc][tileY][tileX][i] *= scale_asym;  
									}
								}
								if (dct_parameters.dbg_mode == 0){ // normalize sym kernel regardless of asym:
									double scale_sym = 0.0;
									for (int i = 0; i< dct_size; i++){
										for (int j = 0; j< dct_size; j++){
											double d = kernels[chn].st_kernels[nc][tileY][tileX][i*dct_size+j];
											if (i > 0) d*=2;
											if (j > 0) d*=2;
											scale_sym +=d;
										}
									}
									for (int i=0; i < kernels[chn].st_kernels[nc][tileY][tileX].length;i++) {
										kernels[chn].st_kernels[nc][tileY][tileX][i] /= scale_sym;  
									}
								}
								// Make a copy of direct kernels (debug feature, may be removed later)
								for (int i = 0; i < dct_size;i++){
									System.arraycopy( // copy one kernel line
											kernels[chn].st_kernels[nc][tileY][tileX],
											i * dct_size,
											kernels[chn].st_direct[nc][tileY][tileX],
											i * dct_size,
											dct_size);
								}
								// scale so multiplication will not change normalization
								for (int i=0; i < kernels[chn].st_kernels[nc][tileY][tileX].length;i++) {
									kernels[chn].st_kernels[nc][tileY][tileX][i] *= dct_size;  
								}

//								kernels[chn].st_kernels[nc][tileY][tileX]= dtt.dttt_iii(kernels[chn].st_kernels[nc][tileY][tileX]);
								kernels[chn].st_kernels[nc][tileY][tileX]= dtt.dttt_iiie(kernels[chn].st_kernels[nc][tileY][tileX]);
							}
//							System.out.println("tileY="+tileY);
						}
					}
					// Debug will be removed later, the
					if (debugLevel > 0) {
						sdfa_instance.showArrays(kernels[chn].sym_kernels,  sym_width, sym_height, true, symKernelPaths[chn]);
						sdfa_instance.showArrays(kernels[chn].asym_kernels,  asym_width, asym_height, true, asymKernelPaths[chn]);
					}
					kernels[chn].sym_kernels = null;  // not needed anymore
					kernels[chn].asym_kernels = null; // not needed anymore
				}
			}
			return true;
		}
		
		public void showKernels(){
//			System.out.println("showKernels(): kernels.length="+kernels.length);
			for (int chn=0;chn < kernels.length; chn++){
				if (kernels[chn]!=null){
//					System.out.println("showKernels("+chn+")");
					 showKernels(chn);
				}
			}
		}
		
		public void showKernels(int chn){
			int nColors = kernels[chn].st_kernels.length;
			int numVert = kernels[chn].st_kernels[0].length;
			int numHor =  kernels[chn].st_kernels[0][0].length;
			int dct_size = kernels[chn].dct_size;
			int asym_size= kernels[chn].asym_size;
			int sym_width =  numHor * dct_size;
			int sym_height = numVert * dct_size;
//			int asym_nonzero =kernels[chn].asym_nonzero;
			int asym_width =  numHor * kernels[chn].asym_size;
			int asym_height = numVert * kernels[chn].asym_size;
			kernels[chn].sym_kernels =  new double [nColors][sym_width*sym_height];
			if (kernels[chn].st_direct != null) {
				kernels[chn].sym_direct =   new double [nColors][sym_width*sym_height];
			}
			kernels[chn].asym_kernels = new double [nColors][asym_width*asym_height];
			int sym_kernel_inc_index =   numHor * dct_size;
			int asym_kernel_inc_index =   numHor * asym_size;
		    showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
			for (int nc = 0; nc < nColors; nc++){
				for (int tileY = 0; tileY < numVert; tileY++){
					for (int tileX = 0; tileX < numHor; tileX++){
						// set DCT (symmetrical) kernels
						int sym_kernel_start_index = (sym_kernel_inc_index * tileY + tileX) * dct_size;
						for (int i = 0; i < dct_size;i++){
							System.arraycopy( // copy one kernel line
									kernels[chn].st_kernels[nc][tileY][tileX],
									i * dct_size,
									kernels[chn].sym_kernels[nc],
									sym_kernel_start_index + i * sym_kernel_inc_index,
									dct_size);
						}
						if (kernels[chn].st_direct != null){
							for (int i = 0; i < dct_size;i++){
								System.arraycopy( // copy one kernel line
										kernels[chn].st_direct[nc][tileY][tileX],
										i * dct_size,
										kernels[chn].sym_direct[nc],
										sym_kernel_start_index + i * sym_kernel_inc_index,
										dct_size);
							}
						}
						
						// set asymmetrical kernel from the list of values and indices
						double [] asym_kernel = new double[asym_size*asym_size];
						for (int i=0;i<asym_kernel.length; i++) asym_kernel[i] = 0.0;
						for (int i = 0; i<kernels[chn].asym_indx[nc][tileY][tileX].length; i++){
							asym_kernel[kernels[chn].asym_indx[nc][tileY][tileX][i]] = kernels[chn].asym_val[nc][tileY][tileX][i];
							/*
							if ((tileY==67) && (tileX==125)) {
								System.out.println("kernels["+chn+"].asym_val["+nc+"]["+tileY+"]["+tileX+"]["+i+"]="+kernels[chn].asym_val[nc][tileY][tileX][i]);
								System.out.println("kernels["+chn+"].asym_indx["+nc+"]["+tileY+"]["+tileX+"]["+i+"]="+kernels[chn].asym_indx[nc][tileY][tileX][i]);
								System.out.println("asym_kernel["+(kernels[chn].asym_indx[nc][tileY][tileX][i])+"]="+asym_kernel[kernels[chn].asym_indx[nc][tileY][tileX][i]]);
							}
							*/
						}
						int asym_kernel_start_index = (asym_kernel_inc_index * tileY + tileX)* asym_size;
						for (int i = 0; i < asym_size;i++){
							System.arraycopy( // copy one kernel line
									asym_kernel,
									i * asym_size,
									kernels[chn].asym_kernels[nc],
									asym_kernel_start_index + i * asym_kernel_inc_index,
									asym_size);
						}
					}
				}
			}
			System.out.println("sym_width="+sym_width+" sym_height="+sym_height);
			System.out.println("kernels["+chn+"].sym_kernels.length="+kernels[chn].sym_kernels.length);
			System.out.println("kernels["+chn+"][0].sym_kernels.length="+kernels[chn].sym_kernels[0].length);
			System.out.println("asym_width="+asym_width+" asym_height="+asym_height);
			System.out.println("kernels["+chn+"].asym_kernels.length="+kernels[chn].asym_kernels.length);
			System.out.println("kernels["+chn+"][0].asym_kernels.length="+kernels[chn].asym_kernels[0].length);
			sdfa_instance.showArrays(kernels[chn].sym_kernels,  sym_width, sym_height, true, "restored-sym-"+chn);
			sdfa_instance.showArrays(kernels[chn].asym_kernels,  asym_width, asym_height, true, "restored-asym-"+chn);
			if (kernels[chn].st_direct != null){
				sdfa_instance.showArrays(kernels[chn].sym_direct,  sym_width, sym_height, true, "restored-direct-"+chn);
			}
			kernels[chn].sym_kernels = null;  // not needed anymore
			kernels[chn].asym_kernels = null; // not needed anymore
			kernels[chn].sym_direct = null;  // not needed anymore
		}

//		public boolean isChannelEnabled(int channel){
//			return ((channel>=0) && (channel<this.usedChannels.length) && this.usedChannels[channel]);  
//		}
		
		
		public void processDCTChannelImages(
//				EyesisCorrectionParameters.SplitParameters         splitParameters,
				EyesisCorrectionParameters.DCTParameters           dct_parameters,
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
				  processDCTChannelImage( // returns ImagePlus, but it already should be saved/shown
						  imp_src, // should have properties "name"(base for saving results), "channel","path"
//						  splitParameters,
						  dct_parameters,
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
				  if (debugLevel>0) System.out.println("Processing image "+(iImage+1)+" (of "+fileIndices.length+") finished at "+
						  IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
				  if (eyesisCorrections.stopRequested.get()>0) {
					  System.out.println("User requested stop");
					  return;
				  }
			  }
			
		}		

		public ImagePlus processDCTChannelImage(
				ImagePlus imp_src, // should have properties "name"(base for saving results), "channel","path"
//				EyesisCorrectionParameters.SplitParameters         splitParameters, // will not be used !
				EyesisCorrectionParameters.DCTParameters           dct_parameters,
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
				double max_vign_corr = dct_parameters.vignetting_range*min_non_zero;
				
				System.out.println("Vignetting data: channel="+channel+", min = "+min_non_zero);
				for (int i=0;i<pixels.length;i++){
					double d = eyesisCorrections.channelVignettingCorrection[channel][i];
					if (d > max_vign_corr) d = max_vign_corr;
					pixels[i]*=d;
				}

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
    				(dct_parameters.decimation == 2)?1:2,  // oversample; // currently source kernels are oversampled
    				dct_parameters.dct_size/2, // addLeft
    				dct_parameters.dct_size/2, // addTop
    				dct_parameters.dct_size/2, // addRight
    				dct_parameters.dct_size/2  // addBottom
    				);		   
			
			// Split into Bayer components, oversample, increase canvas    		  
			ImageStack stack= eyesisCorrections.bayerToStack(
					result, // source Bayer image, linearized, 32-bit (float))
					splitParameters);
			String titleFull=title+"-SPLIT";
			if (!this.correctionsParameters.debayer) {
				result= new ImagePlus(titleFull, stack);    			  
				eyesisCorrections.saveAndShow(result, this.correctionsParameters);
				return result;
			}
// =================
			showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
			if (this.correctionsParameters.deconvolve) { // process with DCT, otherwise use simple debayer
				ImageDtt image_dtt = new ImageDtt();
				double [][][][] dctdc_data = image_dtt.mdctDcStack(
						stack,
						dct_parameters,
						this,
						threadsMax,
						debugLevel,
						updateStatus);
				System.out.println("dctdc_data.length="+dctdc_data.length+" dctdc_data[0].length="+dctdc_data[0].length
						+" dctdc_data[0][0].length="+dctdc_data[0][0].length+" dctdc_data[0][0][0].length="+dctdc_data[0][0][0].length);
				if (dct_parameters.color_DCT){ // convert RBG -> YPrPb
					dctdc_data = image_dtt.dct_color_convert(
							dctdc_data,
							colorProcParameters.kr,
							colorProcParameters.kb,
							dct_parameters.sigma_rb,        // blur of channels 0,1 (r,b) in addition to 2 (g)
							dct_parameters.sigma_y,         // blur of Y from G
							dct_parameters.sigma_color,     // blur of Pr, Pb in addition to Y
							threadsMax,
							debugLevel);
				} else { // just LPF RGB
					for (int chn = 0; chn < dctdc_data.length; chn++) {
						image_dtt.dct_lpf(
								dct_parameters.dbg_sigma,
								dctdc_data[chn],
								threadsMax,
								debugLevel);
					}
				}

				int tilesY = stack.getHeight()/dct_parameters.dct_size - 1;
				int tilesX = stack.getWidth()/dct_parameters.dct_size - 1;
				if (debugLevel>0){
					System.out.println("--tilesX="+tilesX);
					System.out.println("--tilesY="+tilesY);
				}
				double [][] dct_dc = new double [dctdc_data.length][];
				double [][] dct_ac = new double [dctdc_data.length][];
				for (int chn = 0; chn < dct_dc.length; chn++) {
					if (!dct_parameters.color_DCT){ // convert RBG -> YPrPb
						dct_dc[chn] = image_dtt.lapped_dct_dcac(
								false, //       out_ac, // false - output DC, true - output AC
								dctdc_data [chn],
								threadsMax,
								debugLevel);
					}
					dct_ac[chn] = image_dtt.lapped_dct_dcac(
							true, //       out_ac, // false - output DC, true - output AC
							dctdc_data [chn],
							threadsMax,
							debugLevel);
				}
				//	        System.out.println("dct_dc.length="+dct_dc.length+" dct_ac.length="+dct_ac.length);
				if (debugLevel > 0){
					sdfa_instance.showArrays(dct_ac,
							tilesX*dct_parameters.dct_size,
							tilesY*dct_parameters.dct_size,
							true,
							result.getTitle()+"-DCT_AC");
					if (!dct_parameters.color_DCT){ // convert RBG -> YPrPb
						sdfa_instance.showArrays(dct_dc,
								tilesX,
								tilesY,
								true,
								result.getTitle()+"-DCT_DC");
					}
				}
				double [][] idct_data = new double [dctdc_data.length][];
				for (int chn=0; chn<idct_data.length;chn++){
					idct_data[chn] = image_dtt.lapped_idctdc(
							dctdc_data[chn],                  // scanline representation of dcd data, organized as dct_size x dct_size tiles  
							dct_parameters.dct_size,        // final int
							dct_parameters.dct_window,      //window_type
							threadsMax,
							debugLevel);
				}
				if (dct_parameters.color_DCT){ // convert RBG -> YPrPb
					sdfa_instance.showArrays(
							idct_data,
							(tilesX + 1) * dct_parameters.dct_size,
							(tilesY + 1) * dct_parameters.dct_size,
							true,
							result.getTitle()+"-IDCTDC-YPrPb");
					if (dct_parameters.nonlin && ((dct_parameters.nonlin_y != 0.0) || (dct_parameters.nonlin_c != 0.0))) {
						System.out.println("Applying edge emphasis, nonlin_y="+dct_parameters.nonlin_y+
								", nonlin_c="+dct_parameters.nonlin_c+", nonlin_corn="+dct_parameters.nonlin_corn);
						idct_data = edge_emphasis(
								idct_data,                              // final double [][] yPrPb,
								(tilesX + 1) * dct_parameters.dct_size, // final int         width,
								dct_parameters.dct_size,                // final int         step,  //(does not need to be this)            // just for multi-threading efficiency?
								dct_parameters.nonlin_max_y,            // final double      nonlin_max_y =     1.0; // maximal amount of nonlinear line/edge emphasis for Y component
								dct_parameters.nonlin_max_c,            // final double      nonlin_max_c =     1.0; // maximal amount of nonlinear line/edge emphasis for C component
								dct_parameters.nonlin_y,                // final double      nonlin_y,         //  =        0.01; // amount of nonlinear line/edge emphasis for Y component
								dct_parameters.nonlin_c,                // final double      nonlin_c,         //  =        0.01; // amount of nonlinear line/edge emphasis for C component
								dct_parameters.nonlin_corn,             // final double      nonlin_corn,      // =     0.5;  // relative weight for nonlinear corner elements
								(dct_parameters.denoise? dct_parameters.denoise_y:0.0), // final double      denoise_y,        // =        1.0;  // maximal total smoothing of the Y post-kernel (will compete with edge emphasis)
								(dct_parameters.denoise? dct_parameters.denoise_c:0.0), // final double      denoise_c,        //  =        1.0;  // maximal total smoothing of the color differences post-kernel (will compete with edge emphasis)
								dct_parameters.denoise_y_corn,          // final double      denoise_y_corn,   // =   0.5;  // weight of the 4 corner pixels during denoise y (relative to 4 straight)
								dct_parameters.denoise_c_corn,          // final double      denoise_c_corn,   // =   0.5;  // weight of the 4 corner pixels during denoise y (relative to 4 straight)
								threadsMax,                             // final int         threadsMax,       // maximal number of threads to launch                         
								debugLevel);                            // final int         globalDebugLevel)
						sdfa_instance.showArrays(
								idct_data,
								(tilesX + 1) * dct_parameters.dct_size,
								(tilesY + 1) * dct_parameters.dct_size,
								true,
								result.getTitle()+"-EMPH-"+dct_parameters.nonlin_y+"_"+dct_parameters.nonlin_c+"_"+dct_parameters.nonlin_corn);
					}

					// temporary convert back to RGB
					idct_data = YPrPbToRBG(idct_data,
							colorProcParameters.kr,        // 0.299;
							colorProcParameters.kb,        // 0.114;
							(tilesX + 1) * dct_parameters.dct_size);

				} else {
					System.out.println("Bypassing nonlinear correction");
				}

				sdfa_instance.showArrays(idct_data,
						(tilesX + 1) * dct_parameters.dct_size,
						(tilesY + 1) * dct_parameters.dct_size,
						true,
						result.getTitle()+"-IDCTDC-RGB");

				// convert to ImageStack of 3 slices
				String [] sliceNames = {"red", "blue", "green"};
				stack = sdfa_instance.makeStack(
						idct_data,
						(tilesX + 1) * dct_parameters.dct_size,
						(tilesY + 1) * dct_parameters.dct_size,
						sliceNames); // or use null to get chn-nn slice names


			} else { // if (this.correctionsParameters.deconvolve) - here use a simple debayer
				System.out.println("Bypassing DCTR-based aberration correction");
				debayer_rbg(stack, 0.25); // simple standard 3x3 kernel debayer
			}

			if (!this.correctionsParameters.colorProc){
				result= new ImagePlus(titleFull, stack);    			  
				eyesisCorrections.saveAndShow(
						result,
						this.correctionsParameters);
				return result;
			}
			System.out.println("before colors.1");
			//Processing colors - changing stack sequence to r-g-b (was r-b-g)
			if (!eyesisCorrections.fixSliceSequence(
					stack,
					debugLevel)){
				if (debugLevel > -1) System.out.println("fixSliceSequence() returned false");
				return null;
			}
			System.out.println("before colors.2");
			if (debugLevel > -1){
				ImagePlus imp_dbg=new ImagePlus(imp_src.getTitle()+"-"+channel+"-preColors",stack);
				eyesisCorrections.saveAndShow(
						imp_dbg,
						this.correctionsParameters);
			}
			System.out.println("before colors.3, scaleExposure="+scaleExposure+" scale = "+(255.0/eyesisCorrections.psfSubpixelShouldBe4/eyesisCorrections.psfSubpixelShouldBe4/scaleExposure));
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
			if (debugLevel > -1) System.out.println("Processed colors to YPbPr, total number of slices="+stack.getSize());
			if (debugLevel > 0){
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
			//			}
			/*
			String [] slice_names_YPrPb={"Y","Pr","Pb"};
			sdfa_instance.showArrays(idct_data,
					stack.getWidth(), // (tilesX + 1) * dct_parameters.dct_size,
					stack.getHeight(), // (tilesY + 1) * dct_parameters.dct_size,
					true,
					result.getTitle()+"-YPrPb",
					slice_names_YPrPb);
			*/		


			if (toRGB) {
				System.out.println("correctionColorProc.YPrPbToRGB");
				stack =  YPrPbToRGB(yPrPb,
						colorProcParameters.kr,        // 0.299;
						colorProcParameters.kb,        // 0.114;
						stack.getWidth());

				/*
				correctionColorProc.YPrPbToRGB(stack,
						colorProcParameters.kr,        // 0.299;
						colorProcParameters.kb,        // 0.114;
						colorProcParameters.useFirstY?9:8,        //  int sliceY,
								6, // int slicePr,
								7// int slicePb
						);

				 */

				title=titleFull; // including "-DECONV" or "-COMBO"
				titleFull=title+"-RGB-float";
				//Trim stack to just first 3 slices
				if (debugLevel > 0){ // 2){
					ImagePlus imp_dbg=new ImagePlus("YPrPbToRGB",stack);
					eyesisCorrections.saveAndShow(
							imp_dbg,
							this.correctionsParameters);
				}
				while (stack.getSize() > 3) stack.deleteLastSlice();
				if (debugLevel > -1) System.out.println("Trimming color stack");
			} else {
				title=titleFull; // including "-DECONV" or "-COMBO"
				titleFull=title+"-YPrPb"; // including "-DECONV" or "-COMBO"
				if (debugLevel > -1) System.out.println("Using full stack, including YPbPr");
			}

			result= new ImagePlus(titleFull, stack);    			  
			// Crop image to match original one (scaled to oversampling)
			if (crop){ // always crop if equirectangular
				System.out.println("cropping");
				stack = eyesisCorrections.cropStack32(stack,splitParameters);
				if (debugLevel > -1) { // 2){
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
				System.out.println("!toRGB && !this.correctionsParameters.jpeg");
				eyesisCorrections.saveAndShow(result, this.correctionsParameters);
				return result;
			} else { // that's not the end result, save if required
				System.out.println("!toRGB && !this.correctionsParameters.jpeg - else");
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
			System.out.println("result.updateAndDraw(), "+titleFull+"-RGB48");

			CompositeImage compositeImage=eyesisCorrections.convertToComposite(result);

			if (!this.correctionsParameters.jpeg && !advanced){ // RGB48 was the end result
				System.out.println("if (!this.correctionsParameters.jpeg && !advanced)");
				eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters);
				return result;
			} else { // that's not the end result, save if required
				System.out.println("if (!this.correctionsParameters.jpeg && !advanced) - else");
				eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters, this.correctionsParameters.save16, false); // save, no show
				//				  eyesisCorrections.saveAndShow(compositeImage, this.correctionsParameters, this.correctionsParameters.save16, true); // save, no show
			}

			imp_RGB=eyesisCorrections.convertRGB48toRGB24(
					stack,
					title+"-RGB24",
					0, 65536, // r range 0->0, 65536->256
					0, 65536, // g range
					0, 65536);// b range
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

		
		public double [][] edge_emphasis(
				final double [][] yPrPb,
				final int         width,
				final int         step,             // just for multi-threading efficiency?
				final double      nonlin_max_y, // =     1.0; // maximal amount of nonlinear line/edge emphasis for Y component
				final double      nonlin_max_c, // =     1.0; // maximal amount of nonlinear line/edge emphasis for C component
				final double      nonlin_y,         //  =        0.01; // amount of nonlinear line/edge emphasis for Y component
				final double      nonlin_c,         //  =        0.01; // amount of nonlinear line/edge emphasis for C component
				final double      nonlin_corn,      // =     0.5;  // relative weight for nonlinear corner elements
				final double      denoise_y,        // =        1.0;  // maximal total smoothing of the Y post-kernel (will compete with edge emphasis)
				final double      denoise_c,        //  =        1.0;  // maximal total smoothing of the color differences post-kernel (will compete with edge emphasis)
				final double      denoise_y_corn,   // =   0.3;  // weight of the 4 corner pixels during denoise y (relative to 4 straight)
				final double      denoise_c_corn,   // =   0.3;  // weight of the 4 corner pixels during denoise y (relative to 4 straight)
				final int         threadsMax,       // maximal number of threads to launch                         
				final int         globalDebugLevel)
		{
			final double [][] yPrPb_new = new double [yPrPb.length][yPrPb[0].length];
			final int         height = yPrPb[0].length/width;
			final int tilesY = (height + step - 1) / step;
			final int tilesX = (width +  step - 1) / step;
			final int nTiles=tilesX*tilesY;
			final int    [][] probes =  {{1,7},{3,5},{2,6},{0,8}}; // indices in [012/345/678] 3x3 square to calculate squared sums of differences
			final int    [][] kerns =   {{ 1,  3,  5,   7},  {1,   3,   5,  7},  {0,   2,   6,  8},  {0,   2,  6,   8}};  // indices in [012/345/678] 3x3 square to convolve data
			final double [][] kernsw =  {{-1.0,1.0,1.0,-1.0},{1.0,-1.0,-1.0,1.0},{1.0,-1.0,-1.0,1.0},{-1.0,1.0,1.0,-1.0}}; // weights of kern elements
			final int    []   neib_indices = {-width-1,-width,-width+1,-1,0,1,width-1,width,width+1};       
			final double [][] kernsw_y = new double [kernsw.length][];
			final double [][] kernsw_c = new double [kernsw.length][];
			final double []   denoise_kern_y = {
					0.125*denoise_y*denoise_y_corn,         0.125*denoise_y*(1.0-denoise_y_corn), 0.125*denoise_y*denoise_y_corn,
					0.125*denoise_y*(1.0-denoise_y_corn),  -0.5*  denoise_y,                      0.125*denoise_y*(1.0-denoise_y_corn),
					0.125*denoise_y*denoise_y_corn,         0.125*denoise_y*(1.0-denoise_y_corn), 0.125*denoise_y*denoise_y_corn};
			final double []   denoise_kern_c = {
					0.125*denoise_c*denoise_c_corn,         0.125*denoise_c*(1.0-denoise_c_corn), 0.125*denoise_c*denoise_c_corn,
					0.125*denoise_c*(1.0-denoise_c_corn),  -0.5*  denoise_c,                      0.125*denoise_c*(1.0-denoise_c_corn),
					0.125*denoise_c*denoise_c_corn,         0.125*denoise_c*(1.0-denoise_c_corn), 0.125*denoise_c*denoise_c_corn};
			
			for (int n = 0; n < kernsw.length; n++){
				kernsw_y[n] = new double [kernsw[n].length];
				kernsw_c[n] = new double [kernsw[n].length];
				for (int i = 0; i < kernsw[n].length; i++){
					double dy = nonlin_y * ((n>=2)? nonlin_corn: 1.0); 
					double dc = nonlin_c * ((n>=2)? nonlin_corn: 1.0); 
					kernsw_y[n][i] = kernsw[n][i] * dy* dy; 
					kernsw_c[n][i] = kernsw[n][i] * dc* dc; 
				}
			}
			if (globalDebugLevel > 0){
		        System.out.println("edge_emphasis(): nonlin_y="+nonlin_y+
		        		", nonlin_c="+nonlin_c+", nonlin_corn="+nonlin_corn);
		        for (int n=0; n<kernsw_y.length; n++){
		        	System.out.print("kernsw_y["+n+"={");
		        	for (int i = 0; i < kernsw_y[n].length; i++){
			        	System.out.print(kernsw_y[n][i]);
			        	if (i == (kernsw_y[n].length - 1)){
			        		System.out.println("}");
			        	} else {
			        		System.out.print(", ");
			        	}
		        		
		        	}
		        }
		        for (int n=0; n<kernsw_c.length; n++){
		        	System.out.print("kernsw_c["+n+"={");
		        	for (int i = 0; i < kernsw_c[n].length; i++){
			        	System.out.print(kernsw_c[n][i]);
			        	if (i == (kernsw_c[n].length - 1)){
			        		System.out.println("}");
			        	} else {
			        		System.out.print(", ");
			        	}
		        	}
		        }
			}
			
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);

			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						int tileY,tileX;
						double [] neibs =   new double[9]; // pixels around current, first Y, then each color diff
						double [] kern_y =  new double[9]; // weights of neighbors to add to the current for Y
						double [] kern_c =  new double[9]; // weights of neighbors to add to the current for colors diffs
						int center_index = 4;
						for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
							tileY = nTile/tilesX;
							tileX = nTile - tileY * tilesX;
							int y1 = (tileY +1) * step;
							if (y1 > height) y1 = height;
							int x1 = (tileX +1) * step;
							if (x1 > width) x1 = width;
							for (int y = tileY * step; y < y1; y++) {
								for (int x = tileX * step; x < x1; x++) {
									int indx = y*width + x;
									for (int n  = 0; n < yPrPb.length; n++) {
										yPrPb_new[n][indx] = yPrPb[n][indx]; // default - just copy old value
									}
									if ((y > 0) && (y < (height - 1)) && (x > 0) && (x < (width - 1))) { // only correct those not on the edge
										// read Y pixels around current
										for (int i = 0; i < neibs.length; i++){
											neibs[i] = yPrPb[0][indx+neib_indices[i]];
										}
										for (int i = 0; i < kern_y.length; i++){
											kern_y[i] = 0.0;
											kern_c[i] = 0.0;
										}
										// calculate amount of each pattern
										for (int n = 0; n < probes.length; n++){
											double squared=0.0;
											for (int i = 0; i < probes[n].length; i++){
												squared += neibs[probes[n][i]];
											}
											squared -= probes[n].length * neibs[center_index];
											squared *= squared; // Now it a square of the sum of differences from the center
											for (int i = 0; i < kernsw[n].length; i++){
												int ni = kerns[n][i];
												kern_y[ni]+= squared*kernsw_y[n][i];
												kern_c[ni]+= squared*kernsw_c[n][i];
											}
										}
										if (denoise_y > 0.0){
											for (int i = 0; i < kern_y.length; i++){
												kern_y[i]+= denoise_kern_y[i];
											}
										}
										if (denoise_c > 0.0){
											for (int i = 0; i < kern_y.length; i++){
												kern_c[i]+= denoise_kern_c[i];
											}
										}

										if (nonlin_max_y != 0){
											double sum = 0;
											for (int i = 0; i < kern_y.length; i++){
												sum += Math.abs(kern_y[i]);
											}
											if (sum > nonlin_max_y){
												sum = nonlin_max_y/sum;
												for (int i = 0; i < kern_y.length; i++){
													kern_y[i] *= sum;
												}
											}
										}
										
										if (nonlin_max_c != 0){
											double sum = 0;
											for (int i = 0; i < kern_c.length; i++){
												sum += Math.abs(kern_c[i]);
											}
											if (sum > nonlin_max_c){
												sum = nonlin_max_c/sum;
												for (int i = 0; i < kern_c.length; i++){
													kern_c[i] *= sum;
												}
											}
										}
										
										for (int i = 0; i<kern_y.length; i++){
											yPrPb_new[0][indx] += neibs[i]*kern_y[i]; 
										}
										for (int n = 1;n < 3; n++){ //color components
											for (int i = 0; i < neibs.length; i++){
												yPrPb_new[n][indx] += yPrPb[0][indx+neib_indices[i]]*kern_c[i]; 
											}
										}
									}
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			return yPrPb_new;
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
		
		
		
}
