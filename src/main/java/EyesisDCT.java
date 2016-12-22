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
	public DCTKernel [] kernels = null;
	public EyesisDCT(
			EyesisCorrections eyesisCorrections,
			EyesisCorrectionParameters.CorrectionParameters correctionsParameters,
			EyesisCorrectionParameters.DCTParameters dctParameters
			){
		this.eyesisCorrections= eyesisCorrections;
		this.correctionsParameters = correctionsParameters;
		this.dctParameters= dctParameters;
	}
	public class DCTKernel{
		public int           size = 32;     // kernel (DCT) size
		public int           img_step = 32; // pixel step in the image for each kernel
		public double [][][] offsets = null; // per color, per kernel,per coord
		public double [][]   kern = null;  // kernel image in linescan order
	}
	public boolean createDCTKernels(
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
			kernels = new DCTKernel[eyesisCorrections.usedChannels.length];
			for (int chn=0;chn<eyesisCorrections.usedChannels.length;chn++){
				kernels[chn] = null;
			}
		}
		for (int chn=0;chn<eyesisCorrections.usedChannels.length;chn++){
			if (eyesisCorrections.usedChannels[chn] && (sharpKernelPaths[chn]!=null) && (kernels[chn]==null)){
				ImagePlus imp_kernel_sharp=new ImagePlus(sharpKernelPaths[chn]);
				if (imp_kernel_sharp.getStackSize()<3) {
					System.out.println("Need a 3-layer stack with kernels");
					sharpKernelPaths[chn]=null;
					continue;
				}
				ImageStack kernel_sharp_stack= imp_kernel_sharp.getStack();
				
			}

		}
			
		
		
		
		return true;
		
	}

	  public DCTKernel calculateDCTKernel (
			  final ImageStack kernelStack,  // first stack with 3 colors/slices convolution kernels
			  final int          kernelSize, // 64
			  final double       blurSigma,
			  final int          scaleDown,  // kernels are saved with higher resolution - scale them down by this value (2)
			  final int          dctSize,
			  final int          threadsMax,  // maximal number of threads to launch                         
			  final boolean      updateStatus,
			  final int          globalDebugLevel) // update status info
	  {
		  if (kernelStack==null) return null;
		  final int kernelWidth=kernelStack.getWidth();
		  final int kernelNumHor=kernelWidth/kernelSize;
		  final int kernelNumVert=kernelStack.getHeight()/kernelSize;
		  final int nChn=kernelStack.getSize();
		  final int length=kernelNumHor*kernelNumVert*dctSize*dctSize;// size of kernel data 
		  final DCTKernel dct_kernel = new DCTKernel();
		  dct_kernel.size = dctSize;
		  dct_kernel.img_step = dctSize; // configure?
		  dct_kernel.offsets = new double[nChn][kernelNumHor*kernelNumVert][2];
		  dct_kernel.kern = new double [nChn][kernelNumHor*kernelNumVert*dctSize*dctSize];
		  // currently each 64x64 kernel corresponds to 16x16 original pixels tile, 2 tiles margin each side
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
					  float [] kernelPixels= null; // will be initialized at first use
					  double [] kernel=      new double[kernelSize*kernelSize];
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
						  /* read convolution kernel */
						  extractOneKernel(kernelPixels, //  array of combined square kernels, each 
								  kernel, // will be filled, should have correct size before call
								  kernelNumHor, // number of kernels in a row
								  tileX, // horizontal number of kernel to extract
								  tileY); // vertical number of kernel to extract
						  if (blurSigma>0) gb.blurDouble(kernel, kernelSize, kernelSize, blurSigma, blurSigma, 0.01);
						  /* Calculate sum of squared kernel1  elements */
						  sum=0.0;
						  for (i=0; i<kernel.length;i++) sum+=kernel[i]*kernel[i];
						  
//						  outPixles[chn][tileY*kernelNumHor+tileX]= (float) (Math.sqrt(sum));
//						  System.out.println("Processing kernels, channel "+(chn+1)+" of "+nChn+", row "+(tileY+1)+" of "+kernelNumVert+" sum="+sum);
					  }
				  }
			  };
		  }		      
		  ImageDtt.startAndJoin(threads);
		  if (globalDebugLevel > 1) System.out.println("Threads done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		  /* prepare result stack to return */
		  ImageStack outStack=new ImageStack(kernelNumHor,kernelNumVert);
		  return dct_kernel;
	  }
	
//processChannelImage	
//convolveStackWithKernelStack	
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
	  
	  void extractOneKernel(float [] pixels, //  array of combined square kernels, each 
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
	
	
	

}
