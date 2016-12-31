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
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;


public class EyesisDCT {
	public EyesisCorrections eyesisCorrections = null;
	public EyesisCorrectionParameters.CorrectionParameters correctionsParameters=null;
	public EyesisCorrectionParameters.DCTParameters dctParameters = null;
	public DCTKernels [] kernels = null;
	public ImagePlus eyesisKernelImage = null;
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
		public int           numHor =       164; // number of kernel tiles in a row
		public int           dct_size =       8;  // DCT-II size, sym. kernel square side is 2*dct_size-1 
		public int           asym_size =     15; // asymmetrical convolution limits, odd
		public int           asym_nonzdero =  4; // maximal number of non-zero elements in the asymmetrical kernels 
		public double [][]   sym_kernels = null; // per-color channel, DCT kernels in linescan order
		public double [][]   asym_kernels = null; // per-color channel, asymmetrical kernels (all but asym_nonzdero elements are strictly 0)
	}

	public void setKernelImageFile(ImagePlus img_kernels){
		eyesisKernelImage = img_kernels;
	}
	
	public boolean kernelImageSet(){
		return eyesisKernelImage != null;
	}
	
	public boolean createDCTKernels(
			EyesisCorrectionParameters.DCTParameters dct_parameters,
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
				sdfa_instance.showArrays(kernels.sym_kernels,  sym_width, sym_height, true, imp_kernel_sharp.getTitle()+"-sym");

				int asym_width =  kernels.numHor * kernels.asym_size;
				int asym_height = kernels.asym_kernels[0].length /asym_width;
				sdfa_instance.showArrays(kernels.asym_kernels,  asym_width, asym_height, true, imp_kernel_sharp.getTitle()+"-asym");
			}
		}
		
		return true;
		
	}

	  public DCTKernels calculateDCTKernel (
			  final ImageStack kernelStack,  // first stack with 3 colors/slices convolution kernels
			  final int          kernelSize, // 64
			  final EyesisCorrectionParameters.DCTParameters dct_parameters,			  
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
		  // currently each 64x64 kernel corresponds to 16x16 original pixels tile, 2 tiles margin each side
		  final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		  final AtomicInteger ai = new AtomicInteger(0);
		  final int numberOfKernels=     kernelNumHor*kernelNumVert*nChn;
		  final int numberOfKernelsInChn=kernelNumHor*kernelNumVert;
		  final long startTime = System.nanoTime();
		  System.out.println("calculateDCTKernel():numberOfKernels="+numberOfKernels);
		  for (int ithread = 0; ithread < threads.length; ithread++) {
			  threads[ithread] = new Thread() {
				  public void run() {
					  DoubleGaussianBlur gb=null;
					  if (dct_parameters.decimateSigma > 0)	 gb=new DoubleGaussianBlur();
					  float [] kernelPixels= null; // will be initialized at first use NOT yet?
					  double [] kernel=      new double[kernelSize*kernelSize];
					  int targetSize = dct_parameters.asym_size + 2 * dct_parameters.dct_size - 2;
					  double [] target_kernel = new double [targetSize * targetSize];
					  FactorConvKernel factorConvKernel = new FactorConvKernel();
					  factorConvKernel.setDebugLevel       (0); // globalDebugLevel);
					  factorConvKernel.setTargetWindowMode (dct_parameters.dbg_window_mode, dct_parameters.centerWindowToTarget);
					  factorConvKernel.numIterations =     dct_parameters.LMA_steps;
					  factorConvKernel.setAsymCompactness  (dct_parameters.compactness,	dct_parameters.asym_tax_free);
					  
					  int chn,tileY,tileX;
//					  int chn0=-1;
//					  int i;
//					  double sum;
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
						  
						  reformatKernel(
								  kernel,        // will be blured in-place
								  target_kernel, // expand/crop, blur/decimate result
								  kernelSize,
								  targetSize,
								  dct_parameters.decimation,
								  dct_parameters.decimateSigma,
								  gb);
						  // int numAsym = 
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
/*
 * 		  final int kernelNumHor=kernelWidth/kernelSize;
		  final int kernelNumVert=kernelStack.getHeight()/kernelSize;
		  final int nChn=kernelStack.getSize();

		  dct_kernel.sym_kernels =  new double [nChn][kernelNumHor*kernelNumVert*dct_parameters.dct_size * dct_parameters.dct_size];
		  dct_kernel.asym_kernels = new double [nChn][kernelNumHor*kernelNumVert*dct_parameters.asym_size * dct_parameters.asym_size];
 * 
		System.arraycopy(currentfX, 0, convolved, 0, convolved.length);

    			DCT_PARAMETERS.fact_precision,
    			DCT_PARAMETERS.asym_pixels,
    			DCT_PARAMETERS.asym_distance,
    			DCT_PARAMETERS.seed_size);
	
 */
//processChannelImage	
//convolveStackWithKernelStack	
	  // Remove later, copied here as a reference
	  public ImageStack calculateKernelsNoiseGains (
			  final ImageStack kernelStack1, // first stack with 3 colors/slices convolution kernels
			  final ImageStack kernelStack2, // second stack with 3 colors/slices convolution kernels (or null)
			  final int               size, // 128 - fft size, kernel size should be size/2
			  final double       blurSigma,
			  final int          threadsMax,  // maximal number of threads to launch                         
			  final boolean    updateStatus,
			  final int        globalDebugLevel) // update status info
	  {
		  if (kernelStack1==null) return null;
		  final boolean useDiff= (kernelStack2 != null);
		  final int kernelSize=size/2;
		  final int kernelWidth=kernelStack1.getWidth();
		  final int kernelNumHor=kernelWidth/(size/2);
		  final int kernelNumVert=kernelStack1.getHeight()/(size/2);
		  final int length=kernelNumHor*kernelNumVert;
		  final int nChn=kernelStack1.getSize();
		  final float [][] outPixles=new float[nChn][length]; // GLOBAL same as input
		  int i,j;
		  for (i=0;i<nChn;i++) for (j=0;j<length;j++) outPixles[i][j]=0.0f;
		  final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		  final AtomicInteger ai = new AtomicInteger(0);
		  final int numberOfKernels=     kernelNumHor*kernelNumVert*nChn;
		  final int numberOfKernelsInChn=kernelNumHor*kernelNumVert;
		  final long startTime = System.nanoTime();
		  for (int ithread = 0; ithread < threads.length; ithread++) {
			  threads[ithread] = new Thread() {
				  public void run() {
					  DoubleGaussianBlur gb=null;
					  if (blurSigma>0)	 gb=new DoubleGaussianBlur();
					  float [] kernelPixels1= null; // will be initialized at first use
					  float [] kernelPixels2= null; // will be initialized at first use
					  double [] kernel1=      new double[kernelSize*kernelSize];
					  double [] kernel2=      new double[kernelSize*kernelSize];
					  int chn,tileY,tileX;
					  int chn0=-1;
					  int i;
					  double sum;
					  for (int nTile = ai.getAndIncrement(); nTile < numberOfKernels; nTile = ai.getAndIncrement()) {
						  chn=nTile/numberOfKernelsInChn;
						  tileY =(nTile % numberOfKernelsInChn)/kernelNumHor;
						  tileX = nTile % kernelNumHor;
						  if (tileX==0) {
							  if (updateStatus) IJ.showStatus("Processing kernels, channel "+(chn+1)+" of "+nChn+", row "+(tileY+1)+" of "+kernelNumVert);
							  if (globalDebugLevel>2) System.out.println("Processing kernels, channel "+(chn+1)+" of "+nChn+", row "+(tileY+1)+" of "+kernelNumVert+" : "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
						  }
						  
						  if (chn!=chn0) {
							  kernelPixels1=(float[]) kernelStack1.getPixels(chn+1);
							  if (useDiff) kernelPixels2=(float[]) kernelStack2.getPixels(chn+1);
							  chn0=chn;
						  }
						  /* read convolution kernel */
						  extractOneKernel(kernelPixels1, //  array of combined square kernels, each 
								  kernel1, // will be filled, should have correct size before call
								  kernelNumHor, // number of kernels in a row
								  tileX, // horizontal number of kernel to extract
								  tileY); // vertical number of kernel to extract
						  /* optionally read the second convolution kernel */
						  if (useDiff) {extractOneKernel(kernelPixels2, //  array of combined square kernels, each 
								  kernel2, // will be filled, should have correct size before call
								  kernelNumHor, // number of kernels in a row
								  tileX, // horizontal number of kernel to extract
								  tileY); // vertical number of kernel to extract
						     for (i=0; i<kernel1.length;i++) kernel1[i]-=kernel2[i];
						  }
						  if (blurSigma>0) gb.blurDouble(kernel1, kernelSize, kernelSize, blurSigma, blurSigma, 0.01);
						  /* Calculate sum of squared kernel1  elements */
						  sum=0.0;
						  for (i=0; i<kernel1.length;i++) sum+=kernel1[i]*kernel1[i];
						  outPixles[chn][tileY*kernelNumHor+tileX]= (float) (Math.sqrt(sum));
//						  System.out.println("Processing kernels, channel "+(chn+1)+" of "+nChn+", row "+(tileY+1)+" of "+kernelNumVert+" sum="+sum);
					  }
				  }
			  };
		  }		      
		  ImageDtt.startAndJoin(threads);
		  if (globalDebugLevel > 1) System.out.println("Threads done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		  /* prepare result stack to return */
		  ImageStack outStack=new ImageStack(kernelNumHor,kernelNumVert);
		  for (i=0;i<nChn;i++) {
			  outStack.addSlice(kernelStack1.getSliceLabel(i+1), outPixles[i]);
		  }
		  return outStack;
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
			  int       dst_size,  // typical 15
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

}
