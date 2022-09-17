package com.elphel.imagej.correction;
/**
** -----------------------------------------------------------------------------**
** CorrectionDenoise.java
**
** De-noise methods for aberration correction for Eyesis4pi
** 
**
** Copyright (C) 2012 Elphel, Inc.
**
** -----------------------------------------------------------------------------**
**  
**  CorrectionDenoise.java is free software: you can redistribute it and/or modify
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

import ij.IJ;
import ij.ImageStack;

import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.cameras.EyesisCorrectionParameters;
import com.elphel.imagej.common.DoubleFHT;
import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.ShowDoubleFloatArrays;


public class CorrectionDenoise {
    public AtomicInteger stopRequested=null; // 1 - stop now, 2 - when convenient
    double [] denoiseMask;
    int       denoiseMaskWidth;
    
    int debugLevel=1;
    
    public CorrectionDenoise(AtomicInteger stopRequested){
    	this.stopRequested=stopRequested;
    }
    double [] getDenoiseMask() {return this.denoiseMask;}
    int       getDenoiseMaskWidth() {return this.denoiseMaskWidth;}
    void setDebug(int debugLevel){this.debugLevel=debugLevel;}
    // uses global OUT_PIXELS to accumulate results

	/* ======================================================================== */
	/* Combine two  3-slice image stacks generated from the same source image - one high-res/high noise, other low-res/low noise 
	 * @param nonlinParameters TODO*/
	  public ImageStack combineLoHiStacks(ImageStack        stack_convolved, // ImageStack with the image, convolved with the reversed PSF (sharp but with high noise)
	                                      ImageStack         stack_gaussian, // ImageStack with the image, convolved with the Gaussian (just lateral compensated)  (blurred, but low noise)
	                                      int                      channel,  // number of channel to apply to the min/max. If <0 - do not apply
	                                      EyesisCorrectionParameters. NonlinParameters nonlinParameters, // show mask generated and used
	                                      final double [][]       noiseMask, // 2-d array of kernelsNoiseGain (divide mask by it)
	                                      final int               noiseStep, // linear pixels per noiseMask pixels (32)
	                                      final int              threadsMax, // maximal step in pixels on the maxRadius for 1 angular step (i.e. 0.5)
	                                      final boolean        updateStatus, // update status info
	                                      final int debugLevel){

	      int i,j,k;
	      int imgWidth= stack_convolved.getWidth();
	      int imgHeight=stack_convolved.getHeight();
	      double [] diffGreens= new double [imgWidth*imgHeight];
	      double [] diffGreens1;
	      double filtMin=nonlinParameters.filtMin;
	      double filtMax=nonlinParameters.filtMax;
	      if (channel>=0){
	    	  filtMin*=nonlinParameters.thresholdCorr[channel];
	    	  filtMax*=nonlinParameters.thresholdCorr[channel];
	      }
	/* find number of the green channel - should be called "green", if none - use last */
	      int greenChn=2;
	      for (i=0;i<3;i++) if (stack_convolved.getSliceLabel(i+1).equals("green")){
	        greenChn=i;
	        break;
	      }
	      double d;
	      double max=0.0f;
//	      double average=0.0f;
		  DoubleGaussianBlur gb=new DoubleGaussianBlur();

	      float [] hipassPixels=(float[]) stack_convolved.getPixels(greenChn+1);
	      float [] lopassPixels=(float[]) stack_gaussian.getPixels(greenChn+1);
	/*
	      for (i=0;i<lopassPixels.length;i++) {
	 
	        d=hipassPixels[i]-lopassPixels[i];
	        diffGreens[i]=d*d;
	        if (max<lopassPixels[i]) max=lopassPixels[i];
	      }
	*/      
	      for (i=0;i<lopassPixels.length;i++) {
	    	  
//	          d=hipassPixels[i]-lopassPixels[i];
//	          diffGreens[i]=d*d;
	    	  diffGreens[i]=hipassPixels[i]-lopassPixels[i];
	      }
	      if ((debugLevel>3) && nonlinParameters.showMask) ShowDoubleFloatArrays.showArrays(diffGreens.clone(), imgWidth, imgHeight,"diffGreens-nofilter");
	      if (nonlinParameters.blurSigma>0)	{
		   	  if (debugLevel>1) System.out.println ( "Applying gaussian blur to difference hi/lo pass, blurSigma="+nonlinParameters.blurSigma);
			  gb.blurDouble(diffGreens, imgWidth, imgHeight, nonlinParameters.blurSigma, nonlinParameters.blurSigma, 0.01);
		  }
	      if ((debugLevel>3) && nonlinParameters.showMask) ShowDoubleFloatArrays.showArrays(diffGreens.clone(), imgWidth, imgHeight,"diffGreens-blurred");
	      for (i=0;i<lopassPixels.length;i++) {
	    	  diffGreens[i]=diffGreens[i]*diffGreens[i];
	      }
	      if ((debugLevel>2) && nonlinParameters.showMask) ShowDoubleFloatArrays.showArrays(lopassPixels.clone(), imgWidth, imgHeight,"lopassPixels");
	      
	      for (i=0;i<lopassPixels.length;i++) {
	    	  if (max<lopassPixels[i]) max=lopassPixels[i];
	      }
	   	  if (debugLevel>1) System.out.println ( "max(lopassPixels)="+max);
//	      max*=((float) NONLIN_PARAMETERS.threshold);
	// Make threshold absolute - when (blured) intensity is below thershold, the divisor is not decreasing
	      max=((float) nonlinParameters.threshold);
	      if ((debugLevel>3) && nonlinParameters.showMask) ShowDoubleFloatArrays.showArrays(diffGreens.clone(), imgWidth, imgHeight,"diffGreens-squared");
	      for (i=0;i<lopassPixels.length;i++) {
	        diffGreens[i]/=(float) Math.max(max,lopassPixels[i]);
	      }
//	      if ((debugLevel>3) && nonlinParameters.showMask) ShowDoubleFloatArrays.showArrays(diffGreens.clone(), imgWidth, imgHeight,"diffG-norm-limited");
	      if ((debugLevel>1) && nonlinParameters.showMask) ShowDoubleFloatArrays.showArrays(diffGreens.clone(), imgWidth, imgHeight,"diffG-norm-limited");
	      if (nonlinParameters.useRejectBlocksFilter) { // use frequency domain filtering
	    	  
	    	  double lowpassSigmaFreq=1.0*nonlinParameters.maskFFTSize/(2*Math.PI*nonlinParameters.lowPassSigma); // low pass sigma in frequency domain
	    	  double [] filterFHT = createFilterForBlockArtifacts(
	    			  nonlinParameters.maskFFTSize, // size of square FHT
	    			  nonlinParameters.blockPeriod, // period (pixels) of the block artifacts to reject (32)
	    			  nonlinParameters.rejectFreqSigma, // sigma of the rejection spots ( 0.0 - just zero a single point)
	    			  lowpassSigmaFreq); // sigma of the low pass filter (frequency units, 0.0 - do not filter)
	    	  if ((debugLevel>3) && nonlinParameters.showMask) ShowDoubleFloatArrays.showArrays(filterFHT,"filterFHT");
	// Extend at least by half sliding window in each direction to reduce border effect 	  
	    	  diffGreens1=extendDoubleArrayForSlidingWindow(
	    			  diffGreens, // input pixel array
	    			  imgWidth,  // width of the image
	    			  nonlinParameters.maskFFTSize/2); // size of sliding step (half of the sliding window size)
	    	  int extendedWidth=  extendDimension(imgWidth, (nonlinParameters.maskFFTSize/2));
	          int extendedHeight= extendDimension(imgHeight,(nonlinParameters.maskFFTSize/2));
	         
	          if ((debugLevel>3) && nonlinParameters.showMask) ShowDoubleFloatArrays.showArrays(diffGreens1.clone(),extendedWidth,  extendedHeight,"diffGreens-extended");
	    	  
	// run block rejection filter
	    	  diffGreens1=filterMaskFromBlockArtifacts(
	    			  diffGreens1, // input pixel array
	    			  extendedWidth, // width of the image
	    			  extendedHeight, // width of the image
	    			  nonlinParameters.maskFFTSize, // size of sliding FHT
	    			  filterFHT, // filter to multiply FHT (created once for the whole filter mask)
	    			  threadsMax, // maximal step in pixels on the maxRadius for 1 angular step (i.e. 0.5)
	    			  updateStatus, // update status info
	    			  debugLevel);
	          if ((debugLevel>3) && nonlinParameters.showMask) ShowDoubleFloatArrays.showArrays(diffGreens1.clone(),extendedWidth,
	        		  extendedHeight,"diffGreens-filtered-extended");/**/
	// cut extra margins, crop to original size
	    	  diffGreens1=reducedDoubleArrayAfterSlidingWindow(
	    			  diffGreens1, // input pixel array
	    			  imgWidth,  // width of the image
	    			  imgHeight,
	    			  nonlinParameters.maskFFTSize/2); // size of sliding step (half of the sliding window size)
	          if ((debugLevel>2) && nonlinParameters.showMask) ShowDoubleFloatArrays.showArrays(diffGreens1.clone(),imgWidth,
	        		  imgHeight,"diffGreens-filtered");
	    	  if (nonlinParameters.combineBothModes) {
//	      		DoubleGaussianBlur gb=new DoubleGaussianBlur();
	    		gb.blurDouble(diffGreens, imgWidth, imgHeight, nonlinParameters.lowPassSigma, nonlinParameters.lowPassSigma, 0.01);
	    		for (i=0;i<diffGreens.length;i++){
	    			d=diffGreens[i]*diffGreens1[i];
	    			diffGreens[i]=(d>0)?Math.sqrt(diffGreens[i]*diffGreens1[i]):0.0;
	    		}
	    	  } else {
	    		  diffGreens=diffGreens1; 
	    	  }
	      } else { // just apply low-pass filter to the mask
//	    		DoubleGaussianBlur gb=new DoubleGaussianBlur();
	    		gb.blurDouble(diffGreens, imgWidth, imgHeight, nonlinParameters.lowPassSigma, nonlinParameters.lowPassSigma, 0.01);
	      }
	      if ((debugLevel>1) && nonlinParameters.showMask) ShowDoubleFloatArrays.showArrays(diffGreens.clone(), imgWidth, imgHeight,"diffGreens-filtered");
	      
//	      final double [][]       noiseMask, // 2-d array of kernelsNoiseGain (divide mask by it)
//	      final int               noiseStep, // linear pixels per noiseMask pixels (32)
	/* divide mask by noiseMask, if defined */
	      if (noiseMask!=null) {
	    	  if (debugLevel>1) System.out.println ( "diffGreens.length="+diffGreens.length+" imgWidth="+imgWidth+" noiseMask.length="+noiseMask.length+" noiseMask[0].length="+noiseMask[0].length);
	    	  
	          for (i=0;i<diffGreens.length;i++) {
	        	  j=(i/imgWidth)/noiseStep;
	        	  k=(i%imgWidth)/noiseStep;
	        	  if (j>=noiseMask.length)    j=noiseMask.length-1;
	        	  if (k>=noiseMask[j].length) k=noiseMask[j].length-1;
	        	  diffGreens[i]/=noiseMask[j][k];
	          }    	  
	      }
	      if ((debugLevel>1) && nonlinParameters.showMask) ShowDoubleFloatArrays.showArrays(diffGreens.clone(), imgWidth, imgHeight,"diffGreens-noise");
	      if (nonlinParameters.useRingFilter) {
	    	  diffGreens=ringFilter(diffGreens,            // mask to be filtered
	    			  imgWidth,                            // mask width
	    			  nonlinParameters.minMaxValue*nonlinParameters.filtMax, // min value for the local maximum to be processed (absolute, not relative)
	    			  nonlinParameters.overRingThreshold,  // Ratio of the local maximum to maximal value in a ring to trigger filter
	    			  nonlinParameters.overRingLimit,      // limit for the pixels in the center ring relative to the maximum in a ring
	    			  nonlinParameters.ringIR,             // ring inner radius
	    			  nonlinParameters.ringOR);            // ring outer radius
	          if ((debugLevel>1) && nonlinParameters.showMask) ShowDoubleFloatArrays.showArrays(diffGreens.clone(), imgWidth, imgHeight,"diffGreens-ring");
	      }
	   	  if (debugLevel>1) System.out.println ( "filtMax="+filtMax+" filtMin="+filtMin);
	      d= (float) ( 1.0/(filtMax-filtMin));
	      if (filtMax>filtMin) {
	        for (i=0;i<diffGreens.length;i++) {
	          if (diffGreens[i]<filtMin) diffGreens[i]=0.0f;
	          else if (diffGreens[i]>filtMax) diffGreens[i]=1.0f;
	          else diffGreens[i]=d*(diffGreens[i]- (float) filtMin);
	        }
	      }
//	      if (nonlinParameters.showMask) {
//	    	  ShowDoubleFloatArrays.showArrays(diffGreens, imgWidth, imgHeight,"mask");
//	      }
	     this.denoiseMask=diffGreens;
	     this.denoiseMaskWidth=imgWidth; 

	/* Combine 2 stacks and a mask */
	      return combineStacksWithMask (stack_gaussian,
	                                   stack_convolved, 
	                                        diffGreens);
	   }
	  
	  public double [] reducedDoubleArrayAfterSlidingWindow(
			  double []     ipixels, // input pixel array
			  int           imgWidth,  // width of the image
			  int           imgHeight,
			  int             step){ // size of sliding step (half of the sliding window size)

		  int width=extendDimension(imgWidth,step);
		  double [] pixels=new double[imgWidth*imgHeight];
		  int i,j,base;
		  int index=0;
		  for (i=0;i< imgHeight;i++) {
			  base=width*(i+step)+step;
			  for (j=0;j<imgWidth;j++) pixels[index++]=ipixels[base++];
		  }
		  return pixels;
	  }


	  
	  public double [] filterMaskFromBlockArtifacts(
	//  public float [] filterMaskFromBlockArtifacts(
			  final double []     pixels, // input pixel array
//			  final float []     pixels, // input pixel array
			  final int         imgWidth, // width of the image
			  final int        imgHeight, // width of the image
			  final int             size, // size of sliding FHT
			  final double []     filter, // filter to multiply FHT (created once for the whole filter mask)
			  final int       threadsMax, // maximal step in pixels on the maxRadius for 1 angular step (i.e. 0.5)
			  final boolean updateStatus, // update status info
			  final int debugLevel)
	  {

		  if (debugLevel>1) System.out.println("filterMaskFromBlockArtifacts, imgWidth="+imgWidth);
		  if (debugLevel>1) System.out.println("filterMaskFromBlockArtifacts, imgHeight="+imgHeight);

		  if (pixels==null) return null;
		  final int length=imgWidth*imgHeight;
		  final int step=size/2;
		  final int tilesX=imgWidth/step-1; // horizontal number of overlapping tiles in the source image (should be expanded from the registered one by "step" in each direction)
		  final int tilesY=imgHeight/step-1; // vertical number of overlapping tiles in the source image (should be expanded from the registered one by "step" in each direction)
		  if (debugLevel>1) System.out.println("filterMaskFromBlockArtifacts, tilesX="+tilesX);
		  if (debugLevel>1) System.out.println("filterMaskFromBlockArtifacts, tilesY="+tilesY);
		  
		  int i; //tileX,tileY;

//		  for (i=0;i<length;i++) MASK_LOHIRES[i]=0.0;
//		  MASK_LOHIRES=new float[length];
		  final double [] maskLoHiRes=new double[length];
		  for (i=0;i<length;i++) maskLoHiRes[i]=0.0;
		  final double [] slidingWindow= getSlidingMask(size); // 256x256?
		  final Thread[] threads = newThreadArray(threadsMax);
		  final AtomicInteger ai = new AtomicInteger(0);
		  final int numberOfKernels=tilesY*tilesX;
		  final long startTime = System.nanoTime();
		  for (int ithread = 0; ithread < threads.length; ithread++) {
			  threads[ithread] = new Thread() {
				  public void run() {
					  double [] tile=        new double[size * size ];
					  int tileY,tileX;
					  DoubleFHT       fht_instance =   new DoubleFHT(); // provide DoubleFHT instance to save on initializations (or null)
					  for (int nTile = ai.getAndIncrement(); nTile < numberOfKernels; nTile = ai.getAndIncrement()) {
						  tileY = nTile /tilesX;
						  tileX = nTile % tilesX;
						  if (tileX==0) {
							  if (updateStatus) IJ.showStatus("Filtering noise rejection mask, row "+(tileY+1)+" of "+tilesY);
							  if (debugLevel>1) System.out.println("Filtering noise rejection mask, row "+(tileY+1)+" of "+tilesY+" : "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
						  }
						  extractSquareTile( pixels, // source pixel array,
									  tile, // will be filled, should have correct size before call
									  slidingWindow, // window (same size as the kernel)
									  imgWidth, // width of pixels array
									  tileX*step, // left corner X
									  tileY*step); // top corner Y
						  fht_instance.swapQuadrants(tile);
						  fht_instance.transform(tile);
						  tile=fht_instance.multiply(tile,filter,false);
						  fht_instance.inverseTransform(tile);
						  fht_instance.swapQuadrants(tile);
						  /* accumulate result */
						  /*This is synchronized method. It is possible to make threads to write to non-overlapping regions of the OUT_PIXELS, but as the accumulation
						   * takes just small fraction of several FHTs, it should be OK - reasonable number of threads will spread and not "stay in line"
						   */
						  accumulateSquareTile(maskLoHiRes, //  float pixels array to accumulate tile
									  tile, // data to accumulate to the pixels array
									  imgWidth, // width of pixels array
									  tileX*step, // left corner X
									  tileY*step); // top corner Y
						  
					  }
				  }
			  };
		  }		      
		  startAndJoin(threads);
		  return maskLoHiRes;
	  }
	 
	  
	  
	  
	/* ======================================================================== */
	    public double [] ringFilter(double [] dmask,       // mask to be filtered
	    		                    int       width,       // mask width
	    		                    double    minMaxValue, // min value for the local maximum to be processed (absolute, not relative)
	    		                    double    overRingThreshold, // Ratio of the local maximum to maximal value in a ring to trigger filter
	    		                    double    overRingLimit, // limit for the pixels in the center ring relative to the maximum in a ring
	    		                    double    ringIR,        // ring inner radius
	    		                    double    ringOR) {      // ring outer radius
	    	if (dmask==null) return null;
	    	int margins=(int)Math.ceil(ringIR);
	    	int height=dmask.length/width;
	    	double [] result= dmask.clone();
	    	int i,j,nr,nc,y,x,n;
	    	double ringIR2=ringIR*ringIR;
	    	double ringOR2=ringOR*ringOR;
	    	double r2;
	    	nc=0;
	    	nr=0;
	    	for (i=-margins+1;i<margins;i++) for (j=-margins+1;j<margins;j++) {
	    		r2=i*i+j*j;
	    		if (r2<ringIR2) nc++;
	    		else if (r2<=ringOR2) nr++;
	    	}
	    	if ((nc==0) || (nr==0)) return result; // do not filter; 
	    	int [] indxc=new int[nc];
	    	int [] indxr=new int[nr];
	    	nc=0;
	    	nr=0;
	    	for (i=-margins+1;i<margins;i++) for (j=-margins+1;j<margins;j++) {
	    		r2=i*i+j*j;
	    		if (r2<ringIR2) {
	    			indxc[nc++]=j+width*i;
	    		}
	    		else if (r2<=ringOR2) {
	    			indxr[nr++]=j+width*i;
	    		}
	    	}
	    	int [] neighb={-width,-width+1,1,width+1,width,width-1,-1,-width-1};
	    	int index;
	    	boolean isMax;
	    	double d,maxInRing;
	    	for (y=margins; y<height-margins;y++) for (x=margins; x<width-margins;x++) {
	    		index=y*width+x;
	    		d=dmask[index];
	    		if (d<minMaxValue) continue; // too small value - don't bother to filter
	    		isMax=true;
	    		for (n=0;n<neighb.length;n++) if (dmask[index+neighb[n]]>d) {
	    			isMax=false;
	    			break;
	    		}
	    		if (!isMax) continue; // only process local maximums
	    		maxInRing=dmask[index+indxr[0]];
	    		for (n=1;n<nr;n++) if (dmask[index+indxr[n]]>maxInRing) maxInRing=dmask[index+indxr[n]];
	    		if (d< (maxInRing*overRingThreshold)) continue; // under threshold, nop
	// limit values in the circle
	    		maxInRing*=overRingLimit;
	    		for (n=0; n<nc;n++) if (dmask[index+indxc[n]]>maxInRing) result[index+indxc[n]]=maxInRing;
	    	}
	    	return result;
	    }
	/* ======================================================================== */
	  
	  

	/* ======================================================================== */

	  /* Combine 2 stacks and a mask */
	  public ImageStack combineStacksWithMask (ImageStack stack_bg,
			  ImageStack stack_fg, 
			  //                                                 float [] mask ) {
			  double [] mask ) {

		  ImageStack stack=new ImageStack(stack_bg.getWidth(),stack_bg.getHeight());
		  int slice,i;
		  float [] fpixels;
		  float [] fpixels_bg;
		  float [] fpixels_fg;
		  for (slice=1; slice <=stack_bg.getSize(); slice++) {
			  fpixels_bg= (float[])stack_bg.getPixels(slice);
			  fpixels_fg= (float[])stack_fg.getPixels(slice);
			  fpixels=new float [fpixels_bg.length];
			  for (i=0;i<fpixels_bg.length;i++) fpixels[i]= (float) (mask[i]*fpixels_fg[i]+(1.0f-mask[i])*fpixels_bg[i]);
			  stack.addSlice(stack_fg.getSliceLabel(slice), fpixels);
		  }
		  return stack;
	  }
	   /* ======================================================================== */
	     public double [] getSlidingMask(int size) {
	       double [] mask = new double [size*size];
	       double [] maskLine=new double [size];
	       double k=2*Math.PI/size;
	       int i,j,index;
	       for (i=0;i<size;i++) maskLine[i]= 0.5*(1.0-Math.cos(i*k));
	       index=0;
	       for (i=0;i<size;i++) for (j=0;j<size;j++) mask[index++]=maskLine[i]*maskLine[j];
	       return mask;
	     }


	  /* ======================================================================== */
	  /* Filters mask that selects between hi-res/high-noise deconvolved image and lo-res/lo-noise image convolved with Gaussian
	   * by rejecting frequencies that correspond to multiples of JPEG blocks (here with the current settings it is 32 pixels - twice 16x16 macroblock
	   */
	  public double [] createFilterForBlockArtifacts(
			  final int             size, // size of square FHT
			  final int     rejectPeriod, // period (pixels) of the block artifacts to reject (32)
			  final double   rejectSigma, // sigma of the rejection spots ( 0.0 - just zero a single point)
			  final double   lopassSigma) // sigma of the low pass filter (frequency units, 0.0 - do not filter)
	  {
		  double [] maskFHT=new double[size*size];
		  int freqPeriod=size/rejectPeriod;
		  int i,j;
		  for (i=0;i<size*size;i++) maskFHT[i]=1.0; // Initialize mask
		  double []rejGauss;
		  double k;
		  int pointX,pointY,x,y;
		  if (rejectSigma>=0.0) {
			  int rejSize= (int) (rejectSigma*4)+1;
			  if (rejSize>(freqPeriod/2)) rejSize=(freqPeriod/2);
			  rejGauss=new double[rejSize];
			  rejGauss[0]=1.0;
			  k=1.0/(2*rejectSigma*rejectSigma);
			  for (i=1;i<rejSize;i++) {
				  rejGauss[i]=Math.exp(-(k*i*i));
			  }
			  for (pointY=0;pointY<size;pointY+=freqPeriod) for (pointX=0;pointX<size;pointX+=freqPeriod) if ((pointY>0) || (pointX>0)) {
				  for (y=-rejSize+1;y<rejSize;y++) for (x=-rejSize+1;x<rejSize;x++) {
					  i=(pointY+y+size)%size;
					  j=(pointX+x+size)%size;
					  maskFHT[i*size+j]=1.0-rejGauss[Math.abs(y)]*rejGauss[Math.abs(x)];

				  }
			  }
		  }

		  if (lopassSigma>0) {
			  double [] maskTmp=maskFHT;
			  maskFHT=new double[size*size];
			  for (i=0;i<size*size;i++) maskFHT[i]=0.0; // Initialize mask

			  int lopassSize= (int) (lopassSigma*4)+1;
			  if (lopassSize>(size/2)) lopassSize=(size/2);
			  rejGauss=new double[lopassSize];
			  rejGauss[0]=1.0;
			  k=1.0/(2*lopassSigma*lopassSigma);
			  for (i=1;i<lopassSize;i++) {
				  rejGauss[i]=Math.exp(-(k*i*i));
			  }
			  for (y=-lopassSize+1;y<lopassSize;y++) for (x=-lopassSize+1;x<lopassSize;x++) {
				  i=(y+size)%size;
				  j=(x+size)%size;
				  maskFHT[i*size+j]=maskTmp[i*size+j]*rejGauss[Math.abs(y)]*rejGauss[Math.abs(x)];
			  }
		  }
		  return maskFHT;
	  }
	  private int extendDimension(
			  int dimension,
			  int step) {
		  return step*((int) Math.ceil(dimension/step) +2);
	  }
	  public double [] extendDoubleArrayForSlidingWindow(
			  double []     ipixels, // input pixel array
			  int           imgWidth,  // width of the image
			  int             step){ // size of sliding step (half of the sliding window size)

		  int imgHeight= ipixels.length/imgWidth;  // width of the image
		  int width=extendDimension(imgWidth,step);
		  int height=extendDimension(imgHeight,step);
		  double [] pixels=new double[width*height];
		  int i,j,k;
		  int index=0;
		  for (i=0;i<height;i++) {
			  if ((i<step) || (i>=(imgHeight+step))) {
				  k=(i+1)*width;
				  for (j=i*width;j<k;j++) pixels[j]=0.0;
			  } else {
				  k=i*width+step;
				  for (j=i*width;j<k;j++) pixels[j]=0.0;
				  k=(i+1)*width;
				  for (j=i*width+step+imgWidth;j<k;j++) pixels[j]=0.0;
				  k=i*width+step+imgWidth;
				  for (j=i*width+step;j<k;j++) pixels[j]=ipixels[index++];

			  }
		  }
		  return pixels;
	  }
	  
	    /* ======================================================================== */
	    /**extract and multiply by window function (same size as kernel itself) */
	     void extractSquareTile(float [] pixels, // source pixel array,
	                             double [] tile, // will be filled, should have correct size before call
	                           double [] window, // window (same size as the kernel)
	                                  int width, // width of pixels array
	                                     int x0, // left corner X
	                                     int y0) { // top corner Y
	       int length=tile.length;
	       int size=(int) Math.sqrt(length);
	       int i,j,x,y;
	       int height=pixels.length/width;
	       int index=0;
	       for (i=0;i<size;i++) {
	         y=y0+i;
	         if ((y>=0) && (y<height)) {
	           index=i*size;
	           for (j=0;j<size;j++) {
	            x=x0+j;
	            if ((x>=0) && (x<width)) tile [index]=pixels[y*width+x]*window[index];
	            index++;
	           }
	         }
	       }
	     }
	   /* ======================================================================== */
	     void extractSquareTile(double [] pixels, // source pixel array,
	   		  double [] tile, // will be filled, should have correct size before call
	   		  double [] window, // window (same size as the kernel)
	   		  int width, // width of pixels array
	   		  int x0, // left corner X
	   		  int y0) { // top corner Y
	   	  int length=tile.length;
	   	  int size=(int) Math.sqrt(length);
	   	  int i,j,x,y;
	   	  int height=pixels.length/width;
	   	  int index=0;
	   	  for (i=0;i<size;i++) {
	   		  y=y0+i;
	   		  if ((y>=0) && (y<height)) {
	   			  index=i*size;
	   			  for (j=0;j<size;j++) {
	   				  x=x0+j;
	   				  if ((x>=0) && (x<width)) tile [index]=pixels[y*width+x]*window[index];
	   				  index++;
	   			  }
	   		  }
	   	  }
	     }

	  
	  
	   /* accumulate square tile to the pixel array (tile may extend beyond the array, will be cropped) */
	     synchronized void  accumulateSquareTile(
	   		  float [] pixels, //  float pixels array to accumulate tile
	   		  double []  tile, // data to accumulate to the pixels array
	   		  int       width, // width of pixels array
	   		  int          x0, // left corner X
	   		  int          y0) { // top corner Y
	   	  int length=tile.length;
	   	  int size=(int) Math.sqrt(length);
	   	  int i,j,x,y;
	   	  int height=pixels.length/width;
	   	  int index=0;
	   	  for (i=0;i<size;i++) {
	   		  y=y0+i;
	   		  if ((y>=0) && (y<height)) {
	   			  index=i*size;
	   			  for (j=0;j<size;j++) {
	   				  x=x0+j;
	   				  if ((x>=0) && (x<width)) pixels[y*width+x]+=tile [index];
	   				  index++;
	   			  }
	   		  }
	   	  }
	     }
	     synchronized void  accumulateSquareTile(
	   		  double [] pixels, //  float pixels array to accumulate tile
	   		  double []  tile, // data to accumulate to the pixels array
	   		  int       width, // width of pixels array
	   		  int          x0, // left corner X
	   		  int          y0) { // top corner Y
	   	  int length=tile.length;
	   	  int size=(int) Math.sqrt(length);
	   	  int i,j,x,y;
	   	  int height=pixels.length/width;
	   	  int index=0;
	   	  for (i=0;i<size;i++) {
	   		  y=y0+i;
	   		  if ((y>=0) && (y<height)) {
	   			  index=i*size;
	   			  for (j=0;j<size;j++) {
	   				  x=x0+j;
	   				  if ((x>=0) && (x<width)) pixels[y*width+x]+=tile [index];
	   				  index++;
	   			  }
	   		  }
	   	  }
	     }

		/* ======================================================================== */
		/* Create a Thread[] array as large as the number of processors available.
			 * From Stephan Preibisch's Multithreading.java class. See:
			 * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD
			 */
			private Thread[] newThreadArray(int maxCPUs) {
				int n_cpus = Runtime.getRuntime().availableProcessors();
				if (n_cpus>maxCPUs)n_cpus=maxCPUs;
				return new Thread[n_cpus];
			}
		/* Start all given threads and wait on each of them until all are done.
			 * From Stephan Preibisch's Multithreading.java class. See:
			 * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD
			 */
			public static void startAndJoin(Thread[] threads)
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
