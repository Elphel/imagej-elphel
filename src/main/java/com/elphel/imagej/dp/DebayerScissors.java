package com.elphel.imagej.dp;
/**
** -----------------------------------------------------------------------------**
** DebayerScissors.java
**
** Frequency-domain debosaic methods for aberration correction for Eyesis4pi
** 
**
** Copyright (C) 2012 Elphel, Inc.
**
** -----------------------------------------------------------------------------**
**  
**  DebayerScissors.java is free software: you can redistribute it and/or modify
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

import javax.swing.SwingUtilities;

import com.elphel.imagej.common.DoubleFHT;
import com.elphel.imagej.common.showDoubleFloatArrays;


public class DebayerScissors {
//	showDoubleFloatArrays SDFA_INSTANCE=   new showDoubleFloatArrays();
    public AtomicInteger stopRequested=null; // 1 - stop now, 2 - when convenient
    double [] debayerEnergy;
    int       debayerEnergyWidth;
    
    int debugLevel=1;
    
    public DebayerScissors(AtomicInteger stopRequested){
    	this.stopRequested=stopRequested;
    }
    double [] getDebayerEnergy() {return this.debayerEnergy;}
    int       getDebayerEnergyWidth() {return this.debayerEnergyWidth;}
    void setDebug(int debugLevel){this.debugLevel=debugLevel;}
    // uses global OUT_PIXELS to accumulate results
    public ImageStack aliasScissorsStack (
  		  final ImageStack                        imageStack,  // stack with 3 colors/slices with the image
  		  final EyesisCorrectionParameters.DebayerParameters          debayerParameters, // 64 - fft size
  		  final boolean                generateDebayerEnergy,
  		  final int                               threadsMax, // maximal step in pixels on the maxRadius for 1 angular step (i.e. 0.5)
  		  final boolean                         updateStatus, // update status info
  		  final int globalDebugLevel)

    {
//  	  final int wasDebugLevel=this.debugLevel;
  	  // this debug of a specific tile will work in a single-threading only, because it uses global DEBUG_LEVEL
  	  final int xTileDebug, yTileDebug;
        if (!debayerParameters.debug || (debayerParameters.xDebug<0) || (debayerParameters.xDebug<0)) {
      	  xTileDebug=-1;
      	  yTileDebug=-1;
        } else {
      	  xTileDebug=(debayerParameters.xDebug>=(debayerParameters.size/4))?((debayerParameters.xDebug-debayerParameters.size/4)/(debayerParameters.size/2)):0;
      	  yTileDebug=(debayerParameters.yDebug>=(debayerParameters.size/4))?((debayerParameters.yDebug-debayerParameters.size/4)/(debayerParameters.size/2)):0;
        }


  	  if (imageStack==null) return null;
  	  final int imgWidth=imageStack.getWidth();
  	  final int imgHeight=imageStack.getHeight();
  	  final int length=imgWidth*imgHeight;
  	  final int step=debayerParameters.size/2;
  	  final int tilesX=imgWidth/step-1; // horizontal number of overlapping tiles in the source image (should be expanded from the registerd one by "step" in each direction)
  	  final int tilesY=imgHeight/step-1; // vertical number of overlapping tiles in the source image (should be expanded from the registerd one by "step" in each direction)
  	  final int nChn=imageStack.getSize();
  	  int i,chn; //tileX,tileY;
  	  /* find number of the green channel - should be called "green", if none - use last */
  	  i=nChn-1;
  	  for (chn=0;chn<nChn;chn++) if (imageStack.getSliceLabel(chn+1).equals("green")){
  		  i=chn;
  		  break;
  	  }
  	  final int greenChn=i;
  	  final float [][] outPixles=new float[nChn][length]; // same as input
  	  debayerEnergy=null;
  	  if (generateDebayerEnergy) {
  		  debayerEnergy=new double[tilesY*tilesX];
  	  }

  	  for (chn=0;chn<nChn;chn++) for (i=0;i<length;i++) outPixles[chn][i]=0.0f;
  	  final double [] slidingWindow= getSlidingMask(debayerParameters.size); // 64x64

//  	  outPixles=new float[nChn][length]; // GLOBAL same as input
  	  final Thread[] threads = newThreadArray(threadsMax);
  	  final AtomicInteger ai = new AtomicInteger(0);
  	  final int numberOfKernels=tilesY*tilesX;
  	  
	  int indx,dx,dy,tx,ty,li;
	  final int [] nonOverlapSeq = new int[numberOfKernels];
	  int [] nextFirstFindex=new int[4];
	  indx = 0;
	  li=0;
	  
	  for (dy=0;dy<2;dy++) for (dx=0;dx<2;dx++) {
		  for (ty=dy; ty < tilesY; ty+=2) for (tx=dx; tx < tilesX; tx+=2){
			  nonOverlapSeq[indx++] = ty*tilesX + tx;
		  }
		  nextFirstFindex[li++] = indx;
	  }
	  final AtomicInteger aStopIndex = new AtomicInteger(0);
  	  final long startTime = System.nanoTime();
  	  final AtomicInteger tilesFinishedAtomic = new AtomicInteger(1); // first finished will be 1
  	  if (updateStatus) IJ.showStatus("Reducing sampling aliases, "+tilesY+" rows, threadsMax="+threadsMax);
  	  if (updateStatus) System.out.println("Reducing sampling aliases, "+tilesY+" rows, threadsMax="+threadsMax);
  	  for (li = 0; li < nextFirstFindex.length; li++){
  		  aStopIndex.set(nextFirstFindex[li]);
  		  if (li>0){
  			  ai.set(nextFirstFindex[li-1]);
  		  }
  		  for (int ithread = 0; ithread < threads.length; ithread++) {
  			  threads[ithread] = new Thread() {
  				  public void run() {
  					  double [][] tile=        new double[nChn][debayerParameters.size * debayerParameters.size ];
  					  double [][] both_masks;
  					  float [][] pixels=       new float[nChn][];
  					  int chn,tileY,tileX,i;
  					  for (chn=0;chn<nChn;chn++) pixels[chn]= (float[]) imageStack.getPixels(chn+1);
  					  DoubleFHT       fht_instance =   new DoubleFHT(); // provide DoubleFHT instance to save on initializations (or null)
  					  showDoubleFloatArrays SDFA_instance=null; // just for debugging?

  					  deBayerScissors debayer_instance=new deBayerScissors( debayerParameters.size, // size of the square array, centar is at size/2, size/2, only top half+line will be used
  							  debayerParameters.polarStep, // maximal step in pixels on the maxRadius for 1 angular step (i.e. 0.5)
  							  debayerParameters.debayerRelativeWidthGreen, // result green mask mpy by scaled default (diamond)
  							  debayerParameters.debayerRelativeWidthRedblue, // result red/blue mask mpy by scaled default (square)
  							  debayerParameters.debayerRelativeWidthRedblueMain, // green mask when applied to red/blue, main (center)
  							  debayerParameters.debayerRelativeWidthRedblueClones);// green mask when applied to red/blue, clones 
  					  //  					  for (int nTile0 = ai.getAndIncrement(); nTile0 < numberOfKernels; nTile0 = ai.getAndIncrement()) {
  					  for (int nTile0 = ai.getAndIncrement(); nTile0 < aStopIndex.get(); nTile0 = ai.getAndIncrement()) {
 						  
  						  int nTile = nonOverlapSeq[nTile0];
  						  tileY = nTile /tilesX;
  						  tileX = nTile % tilesX;
  						  if (tileX < 2) {
  							  int trow=(tileY+((2*(tileY & 1) + (tileX & 1))*tilesY))/4;
  							  if (updateStatus) IJ.showStatus("Reducing sampling aliases, row "+(trow+1)+" of "+tilesY);
  							  //  						  System.out.println("(1)Reducing sampling aliases, row "+(tileY+1)+" of "+tilesY+" ("+nTile+"/"+nTile0+") col="+(tileX+1));
  							  if (globalDebugLevel>2) System.out.println("Reducing sampling aliases, row "+(trow+1)+" of "+tilesY+" : "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
  						  }

  						  for (chn=0;chn<nChn;chn++){
  							  extractSquareTile( pixels[chn], // source pixel array,
  									  tile[chn], // will be filled, should have correct size before call
  									  slidingWindow, // window (same size as the kernel)
  									  imgWidth, // width of pixels array
  									  tileX*step, // left corner X
  									  tileY*step); // top corner Y
  						  }

  						  /* Scale green channel x0.5 as there are twice more pixels there as in red or blue. Or move it somewhere else and multiply to original range ? */
  						  for (i=0;i<tile[greenChn].length;i++) tile[greenChn][i]*=0.5;
  						  if ((tileY==yTileDebug) && (tileX==xTileDebug)) {
  							  if (SDFA_instance==null) SDFA_instance=      new showDoubleFloatArrays();
  							  SDFA_instance.showArrays (tile.clone(),debayerParameters.size,debayerParameters.size, "x"+(tileX*step)+"_y"+(tileY*step));
  						  }
  						  for (chn=0;chn<nChn;chn++){
  							  fht_instance.swapQuadrants(tile[chn]);
  							  fht_instance.transform(tile[chn]);
  						  }
  						  if ((tileY==yTileDebug) && (tileX==xTileDebug) && (SDFA_instance!=null)) SDFA_instance.showArrays (tile.clone(),debayerParameters.size,debayerParameters.size, "tile-fht");
  						  both_masks= debayer_instance.aliasScissors(tile[greenChn], // fht array for green, will be masked in-place
  								  debayerParameters.debayerThreshold, // no high frequencies - use default uniform filter
  								  debayerParameters.debayerGamma, // power function applied to the amplitudes before generating spectral masks
  								  debayerParameters.debayerBonus, // scale far pixels as (1.0+bonus*r/rmax)
  								  debayerParameters.mainToAlias,// relative main/alias amplitudes to enable pixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
  								  debayerParameters.debayerMaskBlur, // for both masks  sigma for Gaussian blur of the binary masks (<0 -do not use "scissors")
  								  debayerParameters.debayerUseScissors, // use "scissors", if false - just apply "diamond" ands "square" with DEBAYER_PARAMETERS.debayerRelativeWidthGreen and DEBAYER_PARAMETERS.debayerRelativeWidthRedblue
  								  ((tileY==yTileDebug) && (tileX==xTileDebug))?4:1);
  						  //                                               1); // internal debug level ((this.debugLevel>2) && (yTile==yTile0) && (xTile==xTile0))?3:1;
  						  if ((tileY==yTileDebug) && (tileX==xTileDebug) && (SDFA_instance!=null)) {
  							  SDFA_instance.showArrays (tile.clone(),debayerParameters.size,debayerParameters.size, "A00");
  							  SDFA_instance.showArrays (both_masks.clone(),debayerParameters.size,debayerParameters.size, "masks");
  						  }
  						  if (debayerEnergy!=null) {
  							  debayerEnergy[tileY*tilesX+tileX]=debayer_instance.getMidEnergy();
  						  }
  						  for (chn=0;chn<nChn;chn++) {
  							  tile[chn]=fht_instance.multiply(tile[chn],both_masks[(chn==greenChn)?0:1],false);
  							  fht_instance.inverseTransform(tile[chn]);
  							  fht_instance.swapQuadrants(tile[chn]);
  							  /* accumulate result */
  							  /*This is (now was) a synchronized method. It is possible to make threads to write to non-overlapping regions of the outPixles, but as the accumulation
  							   * takes just small fraction of several FHTs, it should be OK - reasonable number of threads will spread and not "stay in line"
  							   */

  							  //accumulateSquareTile(
  							  nonSyncAccumulateSquareTile (
  									  outPixles[chn], //  float pixels array to accumulate tile
  									  tile[chn], // data to accumulate to the pixels array
  									  imgWidth, // width of pixels array
  									  tileX*step, // left corner X
  									  tileY*step); // top corner Y
  						  }
  						  if ((tileY==yTileDebug) && (tileX==xTileDebug) && (SDFA_instance!=null)) SDFA_instance.showArrays (tile.clone(),debayerParameters.size,debayerParameters.size, "B00");
  						  
  						  final int numFinished=tilesFinishedAtomic.getAndIncrement();
//  						  final double dprogr= (1.0+numFinished)/numberOfKernels;
  						  if (numFinished % (numberOfKernels/50+1) == 0) {
//  							  System.out.print(numFinished);
  							  SwingUtilities.invokeLater(new Runnable() {
  								  public void run() {
//  									  System.out.println(" --- "+numFinished+"("+numberOfKernels+"): "+dprogr);
  									  IJ.showProgress(numFinished, numberOfKernels);
  								  }
  							  });
  						  }
  					  }
  				  }
  			  };
  		  }		      
  		  startAndJoin(threads);
 	      IJ.showProgress(tilesFinishedAtomic.get(), numberOfKernels);
  	  }
//  	  System.out.println("tilesFinishedAtomic.get()="+tilesFinishedAtomic.get()+", numberOfKernels="+numberOfKernels);
  	  if (updateStatus) IJ.showStatus("Reducing sampling aliases DONE");
  	  IJ.showProgress(1.0);
  	  /* prepare result stack to return */
  	  ImageStack outStack=new ImageStack(imgWidth,imgHeight);
  	  for (chn=0;chn<nChn;chn++) {
  		  outStack.addSlice(imageStack.getSliceLabel(chn+1), outPixles[chn]);
  	  }
  	  debayerEnergyWidth=	 (debayerEnergy!=null)?tilesX:0; // for the image to be displayed externally 
	  if (globalDebugLevel>0) System.out.println("Reducing sampling aliases done in "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
  	  return outStack;
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

     
   /* ======================================================================== */
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

     void  nonSyncAccumulateSquareTile(
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
