package com.elphel.imagej.common;
/**
 **
 ** DoubleFHT - Determine simulation pattern parameters to match
 ** the acquired image
 **
 ** Copyright (C) 2010-2011 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  MatchSimulatedPattern.java is free software: you can redistribute it and/or modify
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
import java.util.ArrayList;
import java.util.List;

import ij.IJ;
import ij.process.FHT;
//import ij.plugin.FFT;
//import ij.plugin.ContrastEnhancer;
//import java.awt.image.ColorModel;
//import ij.plugin.*;
public class DoubleFHT  {

//	private boolean isFrequencyDomain;
	private int maxSize=20;
	private double [] freqMask=null; // frequency mask cache, will be reused if all parameters are the same
//	private int freqMaskN=0;
//	private double freqHighPass=0;
//	private double freqLowPass=0;
	private double [] freqPass=null; //{freqHighPass,freqLowPass};
	private int maxN=-1; // undefined
	private int ln2=-1; // undefined
	private double[] C; // TODO: make arrays for each reasonable FHT size (or use cache)
	private double[] S;
	private int[] bitrev;
	private double[] tempArr;
	public boolean showProgress=false;
	private double [][][] CS_cache=new double[maxSize][][]; // CS cache for different sizes,
	private double [][] freqMask_cache=new double[maxSize][]; // frequency mask cache, will be reused if all parameters are the same
	private double [][] freqPass_cache=new double[maxSize][]; // cache for low/high frequencies
	private int [][] bitrev_cache= new int [maxSize][];
	private double [][] hamming1d=new double[maxSize][]; // will cache 1d Hamming window;
	private double []   hamming1dMin=new double[maxSize]; // will cache 1d Hamming window;
	private double [][] gaussian1d=new double[maxSize][]; // will cache 1d Gaussian window;
	private double [] gaussian1dWidths=new double[maxSize]; // will cache 1d Gaussian window;

	private double [] translateFHT=null; // cached array for sub-pixel translate
	private double [] translateFHTdXY=null;
	private int translateFHTN=0;
	public boolean debug=false;
	private showDoubleFloatArrays SDFA_INSTANCE= new showDoubleFloatArrays(); // just for debugging?
	public DoubleFHT() {
	  this.C = null;
	  this.S = null;
	  for (int i=0;i<CS_cache.length;i++) {
		  this.CS_cache[i]=null;
		  this.freqMask_cache[i]=null;
		  this.freqPass_cache[i]=null;
		  this.bitrev_cache[i]=null;
		  this.hamming1d[i]=null;
		  this.gaussian1d[i]=null;
		  this.gaussian1dWidths[i]=-1.0;
	  }

	}

	public boolean powerOf2Size(double [] data) {
		return powerOf2Size(data.length);
	}
	public boolean powerOf2Size(int len) {
		int i=4;
		while(i<len) i *= 4;
		return i==len;
	}

 	public void swapQuadrants(double [] data) {
	  if (data==null) return;
	  int size= (int) Math.sqrt(data.length);
	  int hsize=size/2;
	  int shift03=(size+1)*hsize;
	  int shift12=(size-1)*hsize;
	  int i,j,index;
	  double d;
	  for (i=0;i<hsize;i++)  for (j=0;j<hsize;j++) {
	    index=i*size+j;
	    d=data[index];
	    data[index]=data[index+shift03];
	    data[index+shift03]=d;
	    index+=hsize;
	    d=data[index];
	    data[index]=data[index+shift12];
	    data[index+shift12]=d;
	  }
	}


	/* Performs a forward transform, converting this image into the frequency domain.
		The image contained in data must be square and its width must be a power of 2. */
	public boolean transform(double [] data) {
		return transform(data, false);
	}

	/* Performs an inverse transform, converting this image into the space domain.
		The image contained in data must be square and its width must be a power of 2. */
	public boolean inverseTransform(double [] data) {
		return transform(data, true);
	}


	/* Returns an inverse transform of this image, which is assumed to be in the frequency domain. */
	public boolean transform(double [] data, boolean inverse) {
		//IJ.log("transform: "+maxN+" "+inverse);
		updateMaxN(data);
//		maxN = (int) Math.sqrt(data.length);
//		if ((S==null) || (S.length!=(maxN/4))) {
//			makeSinCosTables(maxN);
//			makeBitReverseTable(maxN);
//			tempArr = new double[maxN];
//		}
//		float[] fht = (float[])getPixels();
	 	rc2DFHT(data, inverse, this.maxN);
//		isFrequencyDomain = !inverse;
		return true;
	}

    public double []createFrequencyFilter(
   		 double highPass,
   		 double lowPass){
    	return createFrequencyFilter(null,highPass,lowPass);
    }
/*    public double []createFrequencyFilter(
    		double [] data,
      		 double highPass,
      		 double lowPass){
       	return createFrequencyFilter(data.length,highPass,lowPass);
       }
*/
	public double [] createFrequencyFilter(
    		 double [] data, //int n,
    		 double highPass,
    		 double lowPass){
		 if (data !=null) updateMaxN(data);
//		 int n;

		 if ((this.freqMask!=null) && (this.freqPass[0]==highPass)  && (this.freqPass[0]==lowPass)) return this.freqMask;
		 this.freqMask= new double [(this.maxN+1)*this.maxN/2+1];
    	 double [] lo=new double[this.maxN];
    	 double [] hi=new double[this.maxN];
    	 int i,j;
    	 double kHi= (highPass>0)?(0.5/highPass/highPass):0;
    	 double kLo= (lowPass>0)?(1.0/lowPass/lowPass/this.maxN/this.maxN):0;
    	 double mx2;
//    System.out.println("createFrequencyFilter: highPass ="+highPass+ " lowPass= "+lowPass+ " kHi= "+kHi+ " kLo= "+kLo);

    	 for (i=0;i<this.maxN;i++) {
    		 mx2=-i*i;
    		 hi[i]=(kHi>0)?Math.exp(mx2*kHi):0.0;
    		 lo[i]=(kLo>0)?Math.exp(mx2*kLo):1.0;
    	 }
// more precise version should add data from 4 corners, below is the simple version that may have discontinuities in the middle
    	 for (int index=0;index<this.freqMask.length;index++) {
    		 i=index/this.maxN;
    		 j=index%this.maxN;
    		 if (j>this.maxN/2) j=this.maxN-j;
    		 this.freqMask[index]=(1.0-hi[i]*hi[j])*lo[i]*lo[j];
    	 }
//    	 this.freqMaskN=n;
    	 this.freqPass=new double[2];
    	 this.freqPass[0]=highPass;
    	 this.freqPass[1]=lowPass;
//    	 this.freqMask_cache[this.ln2]=this.freqMask;
//    	 this.freqPass_cache[this.ln2]=this.freqPass;
    	 return this.freqMask;
     }
	/**
	 * Shift the data by phase modification in the frequency domain
	 * @param data square 2**N array in line scan order. Data is modified in-place
	 * @param dx amount of the horizontal shift (left) in pixels
	 * @param dy amount of the vertical shift (down) in pixels
	 * @return modified (shifted) data
	 * TODO: add shift+up-sample combination
	 */
	public double [] shift(double [] data, double dx, double dy){
		return shift(data, 1, dx, dy);
		/*
		updateMaxN(data);
		double sX=2*Math.PI*dx/this.maxN;
		double sY=2*Math.PI*dy/this.maxN;
		int halfN=this.maxN/2;
		double [] cosDX = new double[this.maxN];
		double [] sinDX = new double[this.maxN];
		double [] cosDY = new double[halfN+1];
		double [] sinDY = new double[halfN+1];
		for (int i=0;i<=halfN;i++){ // need less?
			cosDX[i]=Math.cos(sX*i);
			sinDX[i]=Math.sin(sX*i);
			cosDY[i]=Math.cos(sY*i);
			sinDY[i]=Math.sin(sY*i);
		}
		for (int i=1;i<halfN;i++){ // need less?
			cosDX[this.maxN-i]= cosDX[i];
			sinDX[this.maxN-i]=-sinDX[i];
		}
		swapQuadrants(data);
		if (!transform(data,false)) return null; // direct FHT
		for (int row =0; row<=halfN; row++) {
			int rowMod = (this.maxN - row) % this.maxN;
			int maxCol=(row<halfN)?(this.maxN-1):halfN;
			for (int col=0; col<=maxCol; col++) {
				int colMod = (this.maxN - col) % this.maxN;
				int index=   row * this.maxN + col;
				int indexMod=rowMod * this.maxN + colMod;
				double re=0.5*(data[index]+data[indexMod]);
				double im=0.5*(data[index]-data[indexMod]);
				if ((col==halfN) || (row==halfN)) im=0;
				double cosDelta= cosDX[col]*cosDY[row] - sinDX[col]*sinDY[row]; // cos(deltaX)*cos(deltaY)-sin(deltaX)*sin(deltaY)
				double sinDelta= sinDX[col]*cosDY[row] + cosDX[col]*sinDY[row]; // sin(deltaX)*cos(deltaY)+cos(deltaX)*sin(deltaY)
				double reMod=re*cosDelta-im*sinDelta;
				double imMod=re*sinDelta+im*cosDelta;
				data[index]=   reMod+imMod;
				data[indexMod]=reMod-imMod;
			}
		}
		if (!transform(data,true)) return null; // inverse FHT
		swapQuadrants(data);
		return data;
		*/
	}
	/**
	 * Upsample input array by padding in the frequency domain
	 * @param first input square power of 2 array
	 * @param scale upsample scale (power of 2)
	 * @return upsampled array, scan
	 */
	public double [] upsample( double [] first, int scale){
		return shift (first, scale, 0.0, 0.0);
/*
		 if (scale <=1) return first.clone();
		 updateMaxN(first);
    	 swapQuadrants(first);
    	 if (!transform(first,false)) return null; // direct FHT
    	 int halfN=this.maxN/2;
    	 int shift=this.maxN*(scale-1);
    	 int scaledN=this.maxN*scale;
    	 double [] result =new double [first.length*scale*scale];
    	 for (int i=0;i<result.length;i++) result [i]=0.0;
    	 double scale2=scale*scale;
    	 for (int i=0;i<first.length;i++){
    		 int iy=i/this.maxN;
    		 int ix=i%this.maxN;
    		 if (ix>halfN) ix+=shift;
    		 if (iy>halfN) iy+=shift;
    		 result[scaledN*iy+ix]=scale2*first[i];
    	 }
		 updateMaxN(result);
    	 if (!transform(result,true)) return null; // inverse FHT
    	 swapQuadrants(result);
    	 return result;
    	 */

	}

	/**
	 * Upsample data and shift it by phase modification in the frequency domain
	 * @param data square 2**N array in line scan order. Data is modified in-place
	 * @param scale upsample scale
	 * @param dx amount of the horizontal shift (left) in pixels (before scaling)
	 * @param dy amount of the vertical shift (down) in pixels (before scaling)
	 * @return modified (shifted) data
	 */

	public double [] shift(double [] data, int scale, double dx, double dy){
		if ((scale==0) && (dx==0.0) && (dy==0.0)) return data;
		updateMaxN(data);
		swapQuadrants(data);
		if (!transform(data,false)) return null; // direct FHT
		double [] result;
		// scale first
		if (scale >1){
			int halfN=this.maxN/2;
			int shift=this.maxN*(scale-1);
			int scaledN=this.maxN*scale;
			result =new double [data.length*scale*scale];
			for (int i=0;i<result.length;i++) result [i]=0.0;
			double scale2=scale*scale;
			for (int i=0;i<data.length;i++){
				int iy=i/this.maxN;
				int ix=i%this.maxN;
				if (ix>halfN) ix+=shift;
				if (iy>halfN) iy+=shift;
				result[scaledN*iy+ix]=scale2*data[i];
			}
			updateMaxN(result);
		} else {
			result =data;
		}
		// now shift
		if ((dx!=0.0) || (dy!=0.0)) {
			dx*=scale;
			dy*=scale;
			double sX=2*Math.PI*dx/this.maxN;
			double sY=2*Math.PI*dy/this.maxN;
			int halfN=this.maxN/2;
			double [] cosDX = new double[this.maxN];
			double [] sinDX = new double[this.maxN];
			double [] cosDY = new double[halfN+1];
			double [] sinDY = new double[halfN+1];
			for (int i=0;i<=halfN;i++){ // need less?
				cosDX[i]=Math.cos(sX*i);
				sinDX[i]=Math.sin(sX*i);
				cosDY[i]=Math.cos(sY*i);
				sinDY[i]=Math.sin(sY*i);
			}
			for (int i=1;i<halfN;i++){ // need less?
				cosDX[this.maxN-i]= cosDX[i];
				sinDX[this.maxN-i]=-sinDX[i];
			}
			for (int row =0; row<=halfN; row++) {
				int rowMod = (this.maxN - row) % this.maxN;
				int maxCol=(row<halfN)?(this.maxN-1):halfN;
				for (int col=0; col<=maxCol; col++) {
					int colMod = (this.maxN - col) % this.maxN;
					int index=   row * this.maxN + col;
					int indexMod=rowMod * this.maxN + colMod;
					double re=0.5*(data[index]+data[indexMod]);
					double im=0.5*(data[index]-data[indexMod]);
					if ((col==halfN) || (row==halfN)) im=0;
					double cosDelta= cosDX[col]*cosDY[row] - sinDX[col]*sinDY[row]; // cos(deltaX)*cos(deltaY)-sin(deltaX)*sin(deltaY)
					double sinDelta= sinDX[col]*cosDY[row] + cosDX[col]*sinDY[row]; // sin(deltaX)*cos(deltaY)+cos(deltaX)*sin(deltaY)
					double reMod=re*cosDelta-im*sinDelta;
					double imMod=re*sinDelta+im*cosDelta;
					data[index]=   reMod+imMod;
					data[indexMod]=reMod-imMod;
				}
			}
		}
		if (!transform(result,true)) return null; // inverse FHT
		swapQuadrants(result);
		return result;
	}


/* calculates correlation, destroys original arrays */

     public double [] correlate (
    			double [] first,
    			double [] second,
    			double highPassSigma,
    			double lowPassSigma,
    			double phaseCoeff){
 		 updateMaxN(first);
    	 createFrequencyFilter(highPassSigma,lowPassSigma); // for repetitive calls will reuse mask
    	 if (this.freqMask.length<first.length/2){
    		 System.out.println("Error: first.length="+first.length+
    				 " second.length="+second.length+
    				 " this.maxN="+this.maxN+
    				 " this.freqMask.length="+this.freqMask.length);
    	 }
    	 return correlate (first, second,phaseCoeff,this.freqMask);
     }
     public double [] correlate (
    			double [] first,
    			double [] second,
    			double phaseCoeff){
    	    	return correlate (first, second, phaseCoeff,null);
    	     }

     public double [] correlate (
		double [] first,
		double [] second,
		double phaseCoeff,
		double [] filter){
    	if (phaseCoeff<=0.0) return correlate (first, second, filter);
    	else return phaseCorrelate (first, second, phaseCoeff,filter);
     }
     //asymmetrical - will divide by squared second amplitude (pattern to match)
     public double [] phaseCorrelate (
 			double [] first,
 			double [] second,
 			double    phaseCoeff,
 			double [] filter){ //  high/low pass filtering
 	    	 if (first.length!=second.length) {
 	   		   IJ.showMessage("Error","Correlation arrays should be the same size");
 		       return null;
 	    	 }
// 	    	 System.out.println("phaseCorrelate, phaseCoeff="+phaseCoeff);

 	    	 swapQuadrants(first);
 	    	 swapQuadrants(second);
 	    	 if (!transform(first,false)) return null; // direct FHT
 	    	 if (!transform(second,false)) return null; // direct FHT
 	    	 first= phaseMultiply(first, second, phaseCoeff); // correlation, not convolution
 	    	 if (filter!=null) multiplyByReal(first, filter);
 	    	 transform(first,true) ; // inverse transform
 	    	 swapQuadrants(first);
 	    	 return first;
 	     }
     public double [] applyFreqFilter(
    		 double [] first,
    		 double [] filter
    		 ){
    	 if (filter==null) filter=this.freqMask;
    	 if (this.freqMask.length<first.length/2){
    		 System.out.println("Error in applyFreqFilter(): first.length="+first.length+
    				 " this.maxN="+this.maxN+
    				 " filter.length="+filter.length);
    	 }
    	 multiplyByReal(first, filter);
    	 return first;
     }

// modify - instead of a sigma - correlationHighPassSigma (0.0 - regular correlation, 1.0 - phase correlation)

     public double [] correlate (
    			double [] first,
    			double [] second){
    	    	return correlate (first, second, null);
    	     }

     public double [] correlate (
    			double [] first,
    			double [] second,
    			double [] filter){ //  high/low pass filtering
//    	 System.out.println("correlate");
    	    	 if (first.length!=second.length) {
    	   		   IJ.showMessage("Error","Correlation arrays should be the same size");
    		       return null;
    	    	 }
    	    	 swapQuadrants(first);
    	    	 swapQuadrants(second);
    	    	 if (!transform(first,false)) return null; // direct FHT
    	    	 if (!transform(second,false)) return null; // direct FHT
    	    	 first= multiply(first, second, true); // correlation, not convolution
    	    	 if (filter!=null) multiplyByReal(first, filter);
    	    	 transform(first,true) ; // inverse transform
    	    	 swapQuadrants(first);
    	    	 return first;
    	     }
     public double [] convolve (
  			double [] first,
  			double [] second) {
    	 return convolve (first, second, null);
     }

     public double [] convolve (
    		 double [] first,
    		 double [] second,
    		 double [] filter){ //  high/low pass filtering
    	 // 	 System.out.println("correlate");
    	 if (first.length!=second.length) {
    		 IJ.showMessage("Error","Correlation arrays should be the same size");
    		 return null;
    	 }
    	 swapQuadrants(first);
    	 swapQuadrants(second);
    	 if (!transform(first,false)) return null; // direct FHT
    	 if (!transform(second,false)) return null; // direct FHT
    	 first= multiply(first, second, false); // convolution, not correlation
    	 if (filter!=null) multiplyByReal(first, filter);
    	 transform(first,true) ; // inverse transform
    	 swapQuadrants(first);
    	 return first;
     }

     /**
      * Translate real array (zero in the center) by a sub-pixel vector (preferably in +/- 0.5 range for each coordinate)
      * by multiplying spectrum by the linear phase gradient. Caches array, so multiple arrays can be shifted by the same vector faster
      * @param data square array to be translated by a sub-pixel value. in the end it holds FHT (before multiplication)!
      * @param dx translation in x-direction
      * @param dy translation in y direction
      * @return modified data array, translated by specified vector (in-place)
      */

     public double [] translateSubPixel (
    		 double [] data,
    		 double dx,
    		 double dy){
    	 if (debug) SDFA_INSTANCE.showArrays(data, "source-"+IJ.d2s(dx,3)+":"+IJ.d2s(dy,3));
    	 swapQuadrants(data);
    	 if (!transform(data,false)) return null; // direct FHT
    	 if ((this.maxN!=translateFHTN) || (dx!=translateFHTdXY[0]) || (dy!=translateFHTdXY[1])){
    		 calcTranslateFHT(dx,dy);
    	 }
    	 if (debug) SDFA_INSTANCE.showArrays(data, "fht-source-"+IJ.d2s(dx,3)+":"+IJ.d2s(dy,3));
    	 data= multiply(data, this.translateFHT, false); // convolution, not correlation
    	 if (debug) SDFA_INSTANCE.showArrays(data, "fht-multiplied-"+IJ.d2s(dx,3)+":"+IJ.d2s(dy,3));
    	 transform(data,true) ; // inverse transform
    	 swapQuadrants(data);
    	 if (debug) SDFA_INSTANCE.showArrays(data, "result-"+IJ.d2s(dx,3)+":"+IJ.d2s(dy,3));
    	 return data;
     }

     /**
      * Calculating (and caching) FHT array for subpixel translation
      * @param dx translation in x-direction
      * @param dy translation in y direction
      */
     private void calcTranslateFHT(double dx, double dy){
    	 double k=2.0*Math.PI/this.maxN;
    	 double kx=k*dx;
    	 double ky=k*dy;
    	 this.translateFHT=new double[this.maxN*this.maxN];
    	 double sqrt2=Math.sqrt(2.0);
    	 double pi4=Math.PI/4;
    	 int hSize=this.maxN/2;
    	 double a,b,p;
    	 int rowMod, colMod;
    	 for (int r=0;r<hSize;r++){
    		 rowMod = (maxN - r) % maxN;
    		 for (int c=0;c<hSize;c++){
    			 colMod = (maxN - c) % maxN;
    			 p=kx*c+ky*r;
    			 a=sqrt2*Math.sin(pi4+p);
    			 b=sqrt2*Math.sin(pi4-p);
    			 this.translateFHT[r*this.maxN+c]=a;
    			 this.translateFHT[rowMod*this.maxN+colMod]=b;
    			 if (c==0) {
        			 this.translateFHT[r*this.maxN+c+hSize]=a;
        			 this.translateFHT[rowMod*this.maxN+colMod+hSize]=b;
    			 } else {
        			 p=-kx*c+ky*r;
        			 a=sqrt2*Math.sin(pi4+p);
        			 b=sqrt2*Math.sin(pi4-p);
        			 this.translateFHT[r*this.maxN+ colMod]=a;
        			 this.translateFHT[rowMod*this.maxN+c]=b;
    			 }
    		 }
    	 }

// update cache
    	 this.translateFHTdXY=new double[2];
    	 this.translateFHTdXY[0]=dx;
    	 this.translateFHTdXY[1]=dy;
    	 this.translateFHTN=this.maxN;
    	 if (debug) SDFA_INSTANCE.showArrays(this.translateFHT, "translateFHT-"+IJ.d2s(dx,3)+":"+IJ.d2s(dy,3));
     }

 	private boolean updateMaxN(double [] data){
 		if (data==null) return false; // do nothing
		if (!powerOf2Size(data)) {
			String msg="Image is not power of 2 size";
			IJ.showMessage("Error",msg);
			throw new IllegalArgumentException (msg);
		}
		int n=(int) Math.sqrt(data.length);
		boolean differentSize=(n!=this.maxN);
		this.maxN =n;
		if (differentSize){
//System.out.println("Changing FFT size: old ln2="+this.ln2+" new maxN="+this.maxN);
			 // save length old mask and mask parameters
			if (this.ln2>=0){
				this.freqMask_cache[this.ln2]=this.freqMask;
				this.freqPass_cache[this.ln2]=this.freqPass;
			}

//			int ln2;
			for (this.ln2=0;maxN>(1<<this.ln2); this.ln2++);
			if ((this.ln2>=CS_cache.length) || (CS_cache[this.ln2]==null)) {
				if (this.ln2>=CS_cache.length) {
	        		String msg="Too large image, increase this.maxSize (it is now "+this.maxSize+", wanted "+this.ln2+")";
	        		IJ.showMessage("Error",msg);
	        		throw new IllegalArgumentException (msg);
				}
				this.freqMask_cache[this.ln2]=null;
				this.freqPass_cache[this.ln2]=new double[2];
				this.freqPass_cache[this.ln2][0]=0.0;
				this.freqPass_cache[this.ln2][1]=0.0;
				makeSinCosTables(this.maxN);
				this.CS_cache[this.ln2]=new double[2][];
				this.CS_cache[this.ln2][0]=this.C;
				this.CS_cache[this.ln2][1]=this.S;
				makeBitReverseTable(this.maxN);
				this.bitrev_cache[this.ln2]=this.bitrev;

			} else {
				this.C=       this.CS_cache[this.ln2][0];
				this.S=       this.CS_cache[this.ln2][1];
				this.bitrev=  this.bitrev_cache[this.ln2];

			}
			this.freqMask=this.freqMask_cache[this.ln2]; // for the new one there will be just null
			this.freqPass=this.freqPass_cache[this.ln2];
			this.tempArr = new double[maxN];


		}
        return differentSize;

	}
 	public double [] getHamming1d(){
 		if (this.maxN>=0) return getHamming1d(this.maxN);
 		return null;
 	}
 	/**
 	 *
 	 * @param n power of 2 FHT size.
 	 * @return one-dimensional array with Hamming window for low-pass frequency filtering (maximum in the center)
 	 */
 	public double [] getHamming1d(int n){
 		return getHamming1d(n,0.08); // normal Hamming window
 	}

 	public double [] getHamming1d(int n, double min){ // min==0.08 - normal Hamming, 0.0 - pure shifted cosine
		int ln2;
		for (ln2=0;n>(1<<ln2); ln2++);
		if (n!=(1<<ln2)){
    		String msg="Not a power of 2 :"+ n;
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
		}
		if (ln2>=this.hamming1d.length) {
    		String msg="Too large array length, increase this.maxSize (it is now "+this.maxSize+", wanted "+ln2+")";
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
		}
		if ((this.hamming1d[ln2]==null) || (this.hamming1dMin[ln2]!=min)){
			this.hamming1d[ln2]=new double[n];
			this.hamming1dMin[ln2]=min;
			double C054=0.5*(1+min);
			double C046=0.5*(1-min);
			for (int i=0; i<=n/2;i++) this.hamming1d[ln2][i]=  (C054-C046*Math.cos((i*2.0*Math.PI)/n));
			for (int i=1; i<=n/2;i++) this.hamming1d[ln2][n-i]= this.hamming1d[ln2][i];
		}
 		return this.hamming1d[ln2];
 	}

/**
 * Creates one-dimensional array for low-pass filtering in the frequency domain, by multiplying by Gaussian
 * Zero is in the center, same as for Hamming above
 * @param lowPass relative frequency. with lowPass==1.0, Gaussian sigma will be sqrt(2)*Nyquist frequency (==2/n)
 * @param n - number of FFT points (2^2)
 * @return one dimensional array, with 1.0 at 0, and minimum in the middle
 */
 	public double [] getGaussian1d(
 			double lowPass){
 		return getGaussian1d(lowPass,this.maxN);
 	}
 	public double [] getGaussian1d(
 			double lowPass,
 			int n){
		int ln2;
		for (ln2=0;n>(1<<ln2); ln2++);
		if (n!=(1<<ln2)){
    		String msg="Not a power of 2 :"+ n;
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
		}
		if (ln2>=this.gaussian1d.length) {
    		String msg="Too large array length, increase this.maxSize (it is now "+this.maxSize+", wanted "+ln2+")";
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
		}
		if (this.gaussian1d[ln2]==null){
			this.gaussian1d[ln2]=new double[n];
			this.gaussian1dWidths[ln2]=-1.0;
		}
		if (this.gaussian1dWidths[ln2]!=lowPass){
			this.gaussian1dWidths[ln2]=lowPass;
	    	 double kLo= (lowPass>0)?(1.0/lowPass/lowPass/n/n):0;
	    	 double mx2;

	    	 for (int i=0;i<=n/2;i++) {
	    		 mx2=-i*i;
	    		 this.gaussian1d[ln2][n/2-i]=(kLo>0.0)?Math.exp(mx2*kLo):1.0;
	    	 }
			for (int i=1; i<=n/2;i++) this.gaussian1d[ln2][n-i]= this.gaussian1d[ln2][i];
		}
 		return this.gaussian1d[ln2];
 	}


 	private void makeSinCosTables(int maxN) {
 		int ln2;
 		for (ln2=0;maxN>(1<<ln2); ln2++);
 		int n = maxN/4;
 		this.C = new double[n];
 		this.S = new double[n];
 		double theta = 0.0;
 		double dTheta = 2.0 * Math.PI/maxN;
 		for (int i=0; i<n; i++) {
 			this.C[i] = Math.cos(theta);
 			this.S[i] = Math.sin(theta);
 			theta += dTheta;
 		}
 	}

 	private void makeBitReverseTable(int maxN) {
		bitrev = new int[maxN];
		int nLog2 = log2(maxN);
		for (int i=0; i<maxN; i++)
			bitrev[i] = bitRevX(i, nLog2);
	}

	/* Performs a 2D FHT (Fast Hartley Transform). */
	public void rc2DFHT(double[] x, boolean inverse, int maxN) {
		//IJ.write("FFT: rc2DFHT (row-column Fast Hartley Transform)");
		for (int row=0; row<maxN; row++)
			dfht3(x, row*maxN, inverse, maxN);
		progress(0.4);
		transposeR(x, maxN);
		progress(0.5);
		for (int row=0; row<maxN; row++)
			dfht3(x, row*maxN, inverse, maxN);
		progress(0.7);
		transposeR(x, maxN);
		progress(0.8);

		int mRow, mCol;
		double A,B,C,D,E;
		for (int row=0; row<=maxN/2; row++) { // Now calculate actual Hartley transform
			for (int col=0; col<=maxN/2; col++) {
				mRow = (maxN - row) % maxN;
				mCol = (maxN - col)  % maxN;
				A = x[row * maxN + col];	//  see Bracewell, 'Fast 2D Hartley Transf.' IEEE Procs. 9/86
				B = x[mRow * maxN + col];
				C = x[row * maxN + mCol];
				D = x[mRow * maxN + mCol];
				E = ((A + D) - (B + C)) / 2;
				x[row * maxN + col] = A - E;
				x[mRow * maxN + col] = B + E;
				x[row * maxN + mCol] = C + E;
				x[mRow * maxN + mCol] = D - E;
			}
		}
		progress(0.95);
	}

	void progress(double percent) {
		if (showProgress)
			IJ.showProgress(percent);
	}

	/* Performs an optimized 1D FHT. */
	public void dfht3 (double[] x, int base, boolean inverse, int maxN) {
		int i, stage, gpNum, gpSize, numGps, Nlog2;
		int bfNum, numBfs;
		int Ad0, Ad1, Ad2, Ad3, Ad4, CSAd;
		double rt1, rt2, rt3, rt4;

		Nlog2 = log2(maxN);
		BitRevRArr(x, base, Nlog2, maxN);	//bitReverse the input array
		gpSize = 2;     //first & second stages - do radix 4 butterflies once thru
		numGps = maxN / 4;
		for (gpNum=0; gpNum<numGps; gpNum++)  {
			Ad1 = gpNum * 4;
			Ad2 = Ad1 + 1;
			Ad3 = Ad1 + gpSize;
			Ad4 = Ad2 + gpSize;
			rt1 = x[base+Ad1] + x[base+Ad2];   // a + b
			rt2 = x[base+Ad1] - x[base+Ad2];   // a - b
			rt3 = x[base+Ad3] + x[base+Ad4];   // c + d
			rt4 = x[base+Ad3] - x[base+Ad4];   // c - d
			x[base+Ad1] = rt1 + rt3;      // a + b + (c + d)
			x[base+Ad2] = rt2 + rt4;      // a - b + (c - d)
			x[base+Ad3] = rt1 - rt3;      // a + b - (c + d)
			x[base+Ad4] = rt2 - rt4;      // a - b - (c - d)
		 }

		if (Nlog2 > 2) {
			 // third + stages computed here
			gpSize = 4;
			numBfs = 2;
			numGps = numGps / 2;
			//IJ.write("FFT: dfht3 "+Nlog2+" "+numGps+" "+numBfs);
			for (stage=2; stage<Nlog2; stage++) {
				for (gpNum=0; gpNum<numGps; gpNum++) {
					Ad0 = gpNum * gpSize * 2;
					Ad1 = Ad0;     // 1st butterfly is different from others - no mults needed
					Ad2 = Ad1 + gpSize;
					Ad3 = Ad1 + gpSize / 2;
					Ad4 = Ad3 + gpSize;
					rt1 = x[base+Ad1];
					x[base+Ad1] = x[base+Ad1] + x[base+Ad2];
					x[base+Ad2] = rt1 - x[base+Ad2];
					rt1 = x[base+Ad3];
					x[base+Ad3] = x[base+Ad3] + x[base+Ad4];
					x[base+Ad4] = rt1 - x[base+Ad4];
					for (bfNum=1; bfNum<numBfs; bfNum++) {
					// subsequent BF's dealt with together
						Ad1 = bfNum + Ad0;
						Ad2 = Ad1 + gpSize;
						Ad3 = gpSize - bfNum + Ad0;
						Ad4 = Ad3 + gpSize;

						CSAd = bfNum * numGps;
						if (((base+Ad2)>=x.length) || ((base+Ad4)>=x.length) ||
						(CSAd>=this.C.length) || (CSAd>=this.S.length)) {

								System.out.println("dfht3 Error: ln2="+this.ln2+" maxN="+this.maxN+
										" x.length="+x.length+" this.C.length="+this.C.length+" this.S.length="+this.S.length+
										" base="+base+" Ad2="+Ad2+" Ad4="+Ad4+" CSAd="+CSAd);
						}

						rt1 = x[base+Ad2] * C[CSAd] + x[base+Ad4] * S[CSAd];
						rt2 = x[base+Ad4] * C[CSAd] - x[base+Ad2] * S[CSAd];

						x[base+Ad2] = x[base+Ad1] - rt1;
						x[base+Ad1] = x[base+Ad1] + rt1;
						x[base+Ad4] = x[base+Ad3] + rt2;
						x[base+Ad3] = x[base+Ad3] - rt2;

					} /* end bfNum loop */
				} /* end gpNum loop */
				gpSize *= 2;
				numBfs *= 2;
				numGps = numGps / 2;
			} /* end for all stages */
		} /* end if Nlog2 > 2 */

		if (inverse)  {
			for (i=0; i<maxN; i++)
			x[base+i] = x[base+i] / maxN;
		}
	}

	void transposeR (double[] x, int maxN) {
		int   r, c;
		double  rTemp;

		for (r=0; r<maxN; r++)  {
			for (c=r; c<maxN; c++) {
				if (r != c)  {
					rTemp = x[r*maxN + c];
					x[r*maxN + c] = x[c*maxN + r];
					x[c*maxN + r] = rTemp;
				}
			}
		}
	}

	int log2 (int x) {
		int count = 15;
		while (!btst(x, count))
			count--;
		return count;
	}

	private boolean btst (int  x, int bit) {
		//int mask = 1;
		return ((x & (1<<bit)) != 0);
	}

	void BitRevRArr (double[] x, int base, int bitlen, int maxN) {
		for (int i=0; i<maxN; i++)
			tempArr[i] = x[base+bitrev[i]];
		for (int i=0; i<maxN; i++)
			x[base+i] = tempArr[i];
	}

	private int bitRevX (int  x, int bitlen) {
		int  temp = 0;
		for (int i=0; i<=bitlen; i++)
			if ((x & (1<<i)) !=0)
				temp  |= (1<<(bitlen-i-1));
		return temp & 0x0000ffff;
	}

	public int [] nonZeroAmpIndices(double [] ampMask){
		int numNonZero=0;
		for (int i=0;i<ampMask.length;i++) if (ampMask[i]>0.0) numNonZero++;
		int [] result= new int [numNonZero];
		numNonZero=0;
		for (int i=0;i<ampMask.length;i++) if (ampMask[i]>0.0) result[numNonZero++]=i;
		return result;
	}
	/**
	 * Calculate feature matching quality
	 * @param h FHT array
	 * @param directionAngle direction of the linear phase gradient
	 * @param distance Distance to the linear feature in the direction directionAngle (determines the phase gradient absolute value)
	 * @param distance Distance to the linear feature in the direction directionAngle (determines the phase gradient absolute value)
	 * @param phaseIntegrationWidth Width of the band corresponding to the linear feature
	 * @param amHPhase array with amplitudes and (modified) phases for each index or null. Used to mask out some of the points
	 * @param phaseStepCorr Zero-based linear phase values plus these values approximate real phases
	 * @param bandMask if not null, multiply amplitudes by this mask
	 * @return pair of correlation and self-RMS (to calculate relative feature strength)
	 */
	public double [] linearFeatureStrength(
			double [] h,
			double directionAngle,
			double distance,
			double phaseIntegrationWidth, // zero - no limit
			int [] nonZeroIndices,
			double [][] amHPhase,
			double [] phaseStepCorr,
			double [] bandMask){
		double halfWidth=(phaseIntegrationWidth>0.0)?phaseIntegrationWidth/2:this.maxN;
		double sin=Math.sin(directionAngle);
		double cos=Math.cos(directionAngle);
		int halfN=this.maxN/2;
		double diffPhase=distance*2*Math.PI/this.maxN;
		double corrSum=0.0;
		double sumSquares=0.0;
		int numSamples=0;
		for (int numPoint=0; numPoint<nonZeroIndices.length;numPoint++) if ((amHPhase== null) || (amHPhase[numPoint][0]>0.0)) {
			int index=nonZeroIndices[numPoint];
		    int x=(index+halfN)%this.maxN-halfN;
		    int y=index/this.maxN;
			double t = cos*x+sin*y;
			if (Math.abs(cos*y-sin*x)<=halfWidth) {
				int indexMod=((this.maxN - (index/this.maxN)) % this.maxN) * this.maxN + ((this.maxN - (index%this.maxN)) % this.maxN);
				double re=0.5*(h[index]+h[indexMod]);
				double im=0.5*(h[index]-h[indexMod]);
			    double amp2=re*re+im*im;
				if (bandMask!=null) amp2*=bandMask[numPoint];
//			    sumSquares+=((y!=0)?2:1)*amp2;
			    sumSquares+=amp2;
				double phase=phaseStepCorr[numPoint]+diffPhase*t;
				double cosPhase=Math.cos(phase);
				double sinPhase=Math.sin(phase);
				corrSum+=re*cosPhase+im*sinPhase;
				// debugprint
				numSamples++;
				if (this.debug) System.out.println("   Used x="+x+" y="+y+" cos*y+sin*x="+(cos*y+sin*x)+ "  cos="+cos+" sin="+sin+" t="+t);
			} else {
				if (this.debug) System.out.println("Dropped x="+x+" y="+y+" cos*y+sin*x="+(cos*y+sin*x)+ "  cos="+cos+" sin="+sin+" t="+t);
			}
		}
		double [] result={corrSum, Math.sqrt(numSamples*sumSquares)}; // second is n*RMS
		return result;
	}
	/**
	 * Reconstruct FHT of the linear feature using the masked amplitudes of the original and linear phase approximation
	 * @param h FHT array
	 * @param directionAngle direction of the linear phase gradient
	 * @param distance Distance to the linear feature in the direction directionAngle (determines the phase gradient absolute value)
	 * @param phaseIntegrationWidth Width of the band corresponding to the linear feature
	 * @param nonZeroIndices array of non-zero indices in the FHT array
	 * @param amHPhase array with amplitudes and (modified) phases for each index or null. Used to mask out some of the points
	 * @param phaseStepCorr Zero-based linear phase values plus these values approximate real phases
	 * @param bandMask if not null, multiply amplitudes by this mask
	 * @return FHT Reconstructed array transformed to space domain with quadrants swapped;
	 */
	public double [] reconstructLinearFeature(
			double [] h,
			double directionAngle,
			double distance,
			double phaseIntegrationWidth, // zero - no limit
			int [] nonZeroIndices,
			double [][] amHPhase,
			double [] phaseStepCorr,
			double [] bandMask){
		double halfWidth=(phaseIntegrationWidth>0.0)?phaseIntegrationWidth/2:this.maxN;
		double sin=Math.sin(directionAngle);
		double cos=Math.cos(directionAngle);
		int halfN=this.maxN/2;
		double diffPhase=distance*2*Math.PI/this.maxN;
		double [] reconstructedFHT=new double [h.length];
		for (int i=0;i<reconstructedFHT.length;i++) reconstructedFHT[i]=0.0;
		for (int numPoint=0; numPoint<nonZeroIndices.length;numPoint++)  if ((amHPhase== null) || (amHPhase[numPoint][0]>0.0)){
			int index=nonZeroIndices[numPoint];
		    int x=(index+halfN)%this.maxN-halfN;
		    int y=index/this.maxN;
		    if (Math.abs(cos*y-sin*x)<=halfWidth) {
				int indexMod=((this.maxN - (index/this.maxN)) % this.maxN) * this.maxN + ((this.maxN - (index%this.maxN)) % this.maxN);
				double re=0.5*(h[index]+h[indexMod]);
				double im=0.5*(h[index]-h[indexMod]);
		    	double amp=Math.sqrt(re*re+im*im);
				if (bandMask!=null) amp*=bandMask[numPoint];
		    	double t = cos*x+sin*y;
		    	double phase=phaseStepCorr[numPoint]+diffPhase*t;
		    	double reconstructedRe=amp*Math.cos(phase);
		    	double reconstructedIm=amp*Math.sin(phase);
		    	reconstructedFHT[index]=   reconstructedRe+reconstructedIm;
		    	reconstructedFHT[indexMod]=reconstructedRe-reconstructedIm;
//			    if (this.debug)	System.out.println("reconstructLinearFeature(): x="+x+" y="+y +" index="+index+ " indexMod="+indexMod+ " re="+re+" im="+im +" amp="+amp+
//			    		" t="+t+ " phase="+ phase+" reconstructedRe="+reconstructedRe+" reconstructedIm="+reconstructedIm);

		    }
		}
		transform(reconstructedFHT,true);
		swapQuadrants(reconstructedFHT);
		return reconstructedFHT; // space domain
	}
	/**
	 * Adds phase shift to stitch symmetrical around 0 areas (there is zero at 0)
	 * @param directionAngle perpendicular to the linear feature (direction of high amplitude in frequency domain)
	 * @param zeroPhase extrapolated phase at 0 for positive direction at angle (half step)
	 * @param distance Estimated distance to line (only used with zeroBinHalfSize)
	 * @param phaseTolerance - if >0.0, will zero amplitude if the phase is too off
	 * @param amHPhase [numLixel]{amplitude, phase} - phase will be modified! May be null - will work faster
	 * @param zeroBinHalfSize - if abs(projection) is less than zeroBinHalfSize, use the best (of 2 options) phase
	 * @param nonZeroIndices
	 * @return applied phase correction, added to the modified phase will result actual phase
	 */
	public double [] compensatePhaseStep(
			double directionAngle,
			double zeroPhase,
			double distance,
			double phaseTolerance,
			double [][] amHPhase,
			int [] nonZeroIndices,
			double zeroBinHalfSize){
		double sin=Math.sin(directionAngle);
		double cos=Math.cos(directionAngle);
		int halfN=this.maxN/2;
		double diffPhase=distance*2*Math.PI/this.maxN;
		double [] phaseCorr=new double [nonZeroIndices.length];
		for (int numPoint=0; numPoint<nonZeroIndices.length;numPoint++) {
			int index=nonZeroIndices[numPoint];
			int x=(index+halfN)%this.maxN-halfN;
			int y=index/this.maxN;
			double t = cos*x+sin*y;
			double dPhase=(t>0)?zeroPhase:(-zeroPhase);
			if (amHPhase!=null) {
				if ( (y!=0) && (Math.abs(t)<zeroBinHalfSize)){
					double flattenedPhase=amHPhase[numPoint][1]-diffPhase*t;
					double diffPlus= flattenedPhase-zeroPhase;
					double diffMinus=flattenedPhase+zeroPhase;
					diffPlus-= (2*Math.PI*Math.round(diffPlus/(2*Math.PI)));
					diffMinus-=(2*Math.PI*Math.round(diffMinus/(2*Math.PI)));
					dPhase= (Math.abs(diffPlus)<Math.abs(diffMinus))?zeroPhase:(-zeroPhase);
					/*
				System.out.println("numPoint="+numPoint+" x="+x+" y="+y+
						" phase(was)="+amHPhase[numPoint][1]+
						" flattenedPhase="+flattenedPhase+
						" diffPlus="+diffPlus+" diffMinus="+diffMinus+
						" dPhase="+dPhase);
					 */
				}

				amHPhase[numPoint][1]-=dPhase;
				amHPhase[numPoint][1]-=(2*Math.PI*Math.round(amHPhase[numPoint][1]/(2*Math.PI)));

				double phaseError=amHPhase[numPoint][1]-diffPhase*t;
				phaseError-=(2*Math.PI*Math.round(phaseError/(2*Math.PI)));

				if ((phaseTolerance>0.0) && (Math.abs(phaseError)>phaseTolerance)){
					amHPhase[numPoint][0]=0.0; // bad phase - make zero;
				}
			}
			phaseCorr[numPoint]=dPhase;
		}
		return phaseCorr;
	}

	public int [] getBandIndices(
			boolean keepDC,
			double ribbonWidth,
			double angle){
		int halfN=this.maxN/2;
		double halfWidth=ribbonWidth/2;
		int length=halfN*this.maxN;
		int [] tmpIndices=new int [length];
		double sin=Math.sin(angle);
		double cos=Math.cos(angle);
		int pointNumber=0;
		for (int index=0;index<length;index++){
		    int x=(index+halfN)%this.maxN-halfN;
		    int y=index/this.maxN;
		    double t= cos*x+sin*y;
		    double d= sin*x-cos*y;
		    int it=(int) Math.round(Math.abs(t));
		    if ((keepDC || (y>0) || (x>0)) &&(Math.abs(d) < halfWidth) && (it < halfN) ){
		    	tmpIndices[pointNumber++]=index;
		    }
		}
		int [] indices=new int [pointNumber];
		for (int i=0;i<pointNumber;i++) indices[i]=tmpIndices[i];
		return indices;
	}
	public double [] getRibbonMask(
			boolean keepDC,
			int [] bandIndices,
			double ribbonWidth,
			double resultHighPass0,
			double angle){
		double resultHighPass=keepDC?0.0:resultHighPass0;
		int halfN=this.maxN/2;
		double halfWidth=ribbonWidth/2;
		double sin=Math.sin(angle);
		double cos=Math.cos(angle);
		if (sin<0){
			sin=-sin;
			cos=-cos;
		}
		double [] ribbonMask=new double [bandIndices.length];
		for (int iIndex=0;iIndex<bandIndices.length;iIndex++){
			int index=bandIndices[iIndex];
		    int x=(index+halfN)%this.maxN-halfN;
		    int y=index/this.maxN;
		    double t= Math.abs(cos*x+sin*y);
		    double d= sin*x-cos*y;
//		    int it=(int) Math.round(Math.abs(t));
		    int it=(int) Math.round(t);
		    if ((keepDC || (y>0) || (x>0)) &&(Math.abs(d) < halfWidth) && (it < halfN) ){
		    	double w=(0.54+0.46*Math.cos(Math.PI*d/halfWidth)); // (masked) amplitude;
		    	if (t<resultHighPass)   w*=0.54+0.46*Math.cos(Math.PI*(resultHighPass-t)/resultHighPass);
		    	else if (t>(halfN-halfWidth)) w*=0.54+0.46*Math.cos(Math.PI*(t-(halfN-halfWidth))/halfWidth);
			    ribbonMask[iIndex]=w;
//			    if (this.debug)	System.out.println("getRibbonMask(): x="+x+" y="+y +" t="+t+" d="+d +" w="+w+ " w0="+ (0.54+0.46*Math.cos(Math.PI*d/halfWidth)));
		    } else {
		    	ribbonMask[iIndex]=0.0;
		    	System.out.println("getRibbonMask(): should not get here");
		    }
		}
	    return ribbonMask;
	}


	/**
	 * Detect linear features in the frequency domain sample
	 * @param ribbonWidth Width of the band for processing phase (data inside will be multiplied by the Hamming filter perpendicular to the direction of the band
	 * @param fht frequency-domain FHT array of the sample
	 * @param minRMSFrac Process only frequencies with amplitude above this ratio of the RMS of all frequencies.  (0 - do not use this threshold)
	 * @param minAbs  Process only frequencies with (absolute) amplitude above this value (0 - do not use this threshold)
	 * @param maxPhaseMismatch Maximal mismatch between the diagonal pair of 2 pixel pairs to consider pixels valid (distance between crossing diagonals)
	 * @param dispertionCost - >0 - use phase dispersion when calculating quality of phase match, 0 - only weights.
	 * @param filter Bitmask enabling different phase half-steps (+1 - 0, +2 - pi/2, +4 - pi, +8 - 3pi/2. Value zero allows arbitrary step
	 *        step 0 corresponds to thin white line, pi - thin black line, +pi/2 - edge black-to-white in the direction of directionAngle,
	 *        3pi/2 - white-to-black in the direction of directionAngle
	 * @return array of the following values:
	 *  {bestAngle, angle of the perpendicular towards the linear feature: 0 in the direction of pixel x (right), pi/2 - in the direction of the pixel Y (down)
	 *  distance, - distance from the sample center to the linear feature in the direction bnestAngle (in pixels)
	 *  phaseAtZero, - extrapolated phase at zero frequency (DC) , see description of the parameter "filter" above
	 *  maxStrength}; strength of the feature (combines multiple factors)
	 */

	public double [] calcPhaseApproximation(
			double ribbonWidth,
			double [] fht,
			double minRMSFrac,
			double minAbs,
			double maxPhaseMismatch,
			double dispertionCost,
			int filter
			){
		 createPolarMasksIndices(ribbonWidth); // all but first time will return immediately.
		int halfN=this.maxN/2;
		double dispersionFatZero=(dispertionCost>0)?(Math.PI/halfN/dispertionCost):0.0;
		double [][] dPHdXdY=new double [5][this.maxN*this.maxN/2]; // amp, phase, phaseWeight, dPhase/dX, dPhase/dY
		for (int n=0;n<dPHdXdY.length;n++) for (int i=0;i<dPHdXdY[0].length;i++) dPHdXdY[n][i]=0.0;
		double sum2=0.0;
		for (int index=0;index<dPHdXdY[0].length;index++){
			int indexMod=((this.maxN - (index/this.maxN)) % this.maxN) * this.maxN + ((this.maxN - (index%this.maxN)) % this.maxN);
			double re=0.5*(fht[index]+fht[indexMod]);
			double im=0.5*(fht[index]-fht[indexMod]);
//		    int x=(index+halfN)%this.maxN-halfN;
//		    int y=index/this.maxN;
		    double amp2=re*re+im*im;
		    sum2+=amp2;
		    dPHdXdY[0][index]=Math.sqrt(amp2);
		    dPHdXdY[1][index]=Math.atan2(im,re);
		}
		double rms=Math.sqrt(sum2/dPHdXdY[0].length);
		double thresholdAmplitude=0.0;
		if (this.debug) System.out.println("calc_dPHdXdY): rms="+rms);
		if (minRMSFrac>0) thresholdAmplitude=minRMSFrac*rms;
		if ((minAbs>0) && ((thresholdAmplitude==0) || (thresholdAmplitude>minAbs))) thresholdAmplitude=minAbs;
		for (int index=0;index<dPHdXdY[0].length;index++) if (dPHdXdY[0][index]<thresholdAmplitude) dPHdXdY[0][index]=0.0;
		double oneThird=1.0/3;
		for (int y=0;y<(halfN-1);y++){
			for (int x=-halfN+1; x<(halfN-1);x++){
//				int index00=y*this.maxN+(x+halfN);
				int index00=y*this.maxN+(x+this.maxN)%this.maxN;
				int indexX0=y*this.maxN+(x+this.maxN+1)%this.maxN;
				int index0Y=index00+this.maxN;
				int indexXY=(y+1)*this.maxN+(x+this.maxN+1)%this.maxN;
				if ((dPHdXdY[0][index00]>0) && (dPHdXdY[0][indexX0]>0) && (dPHdXdY[0][index0Y]>0) && (dPHdXdY[0][indexXY]>0)){
					double dP00X0=dPHdXdY[1][indexX0]-dPHdXdY[1][index00];
					double dP000Y=dPHdXdY[1][index0Y]-dPHdXdY[1][index00];
					double dPX0XY=dPHdXdY[1][indexXY]-dPHdXdY[1][indexX0];
					double dP0YXY=dPHdXdY[1][indexXY]-dPHdXdY[1][index0Y];
					dP00X0-=(2*Math.PI)*Math.round(dP00X0/(2*Math.PI));
					dP000Y-=(2*Math.PI)*Math.round(dP000Y/(2*Math.PI));
					dPX0XY-=(2*Math.PI)*Math.round(dPX0XY/(2*Math.PI));
					dP0YXY-=(2*Math.PI)*Math.round(dP0YXY/(2*Math.PI));
					if (this.debug) {
						System.out.println(" index00="+index00+" indexX0="+indexX0+" index0Y="+index0Y+" indexXY="+indexXY);
						System.out.println(" dPHdXdY[1][index00]="+dPHdXdY[1][index00]+" dPHdXdY[1][indexX0]="+dPHdXdY[1][indexX0]+
								" dPHdXdY[1][index0Y]="+dPHdXdY[1][index0Y]+" dPHdXdY[1][indexXY]="+dPHdXdY[1][indexXY]);
						System.out.println(" dP00X0="+dP00X0+" dP000Y="+dP000Y+" dPX0XY="+dPX0XY+" dP0YXY="+dP0YXY);
					}

					double curv=Math.abs((dPX0XY-dP000Y)/2);
					if ((maxPhaseMismatch==0) || (curv < maxPhaseMismatch)){
						double w=Math.pow(dPHdXdY[0][index00]*dPHdXdY[0][indexX0]*dPHdXdY[0][index0Y], oneThird);
						if (this.debug) System.out.print(">>> x="+x+" y="+y+ " curv="+curv);
						dPHdXdY[2][index00]+=w;
						dPHdXdY[3][index00]+=w*dP00X0;
						dPHdXdY[4][index00]+=w*dP000Y;
						if (this.debug) System.out.print(" index00="+index00+" w="+w);
						w=Math.pow(dPHdXdY[0][index00]*dPHdXdY[0][indexX0]*dPHdXdY[0][indexXY], oneThird);
						dPHdXdY[2][indexX0]+=w;
						dPHdXdY[3][indexX0]+=w*dP00X0;
						dPHdXdY[4][indexX0]+=w*dPX0XY;
						if (this.debug) System.out.print(" indexX0="+indexX0+" w="+w);
						w=Math.pow(dPHdXdY[0][index00]*dPHdXdY[0][index0Y]*dPHdXdY[0][indexXY], oneThird);
						dPHdXdY[2][index0Y]+=w;
						dPHdXdY[3][index0Y]+=w*dP0YXY; //dP000Y;
						dPHdXdY[4][index0Y]+=w*dP000Y; //dP0YXY;
						if (this.debug) System.out.print(" index0Y="+index0Y+" w="+w);
						w=Math.pow(dPHdXdY[0][indexX0]*dPHdXdY[0][index0Y]*dPHdXdY[0][indexXY], oneThird);
						dPHdXdY[2][indexXY]+=w;
						dPHdXdY[3][indexXY]+=w*dP0YXY; //dPX0XY;
						dPHdXdY[4][indexXY]+=w*dPX0XY; //dP0YXY;
						if (this.debug) System.out.print(" indexXY="+indexXY+" w="+w);
					} else {
						if (this.debug) System.out.println(">-- x="+x+" y="+y+ " curv="+curv+" w="+Math.pow(dPHdXdY[0][index00]*dPHdXdY[0][indexX0]*dPHdXdY[0][index0Y], oneThird));

					}

				}
			}
		}
		for (int i=1;i<halfN;i++){
			int j=this.maxN-i;
			if ((dPHdXdY[2][i] +dPHdXdY[2][j])>0.0) {
				dPHdXdY[3][i]=(dPHdXdY[3][i]*dPHdXdY[2][i]+dPHdXdY[3][j]*dPHdXdY[2][j])/(dPHdXdY[2][i]+dPHdXdY[2][j]);
				dPHdXdY[4][i]=(dPHdXdY[4][i]*dPHdXdY[2][i]+dPHdXdY[4][j]*dPHdXdY[2][j])/(dPHdXdY[2][i]+dPHdXdY[2][j]);
				dPHdXdY[2][i]=dPHdXdY[2][i]+dPHdXdY[2][j];
				dPHdXdY[3][j]=dPHdXdY[3][i];
				dPHdXdY[4][j]=dPHdXdY[4][i];
				dPHdXdY[2][j]=dPHdXdY[2][i];
			}
		}
		for (int index=0;index<dPHdXdY[0].length;index++) if (dPHdXdY[2][index]>0){
			dPHdXdY[3][index]/=dPHdXdY[2][index];
			dPHdXdY[4][index]/=dPHdXdY[2][index];
		}
		// generate show image
		if (this.debug){
			double [][]  debugdPHdXdY=new double [5][this.maxN*this.maxN];
			for (int i=0;i<debugdPHdXdY[0].length;i++){
				int x= (i%this.maxN)-halfN;
				int y= (i/this.maxN)-halfN;
				boolean other=(y<0);
				if (y<0) {
					y=-y;
					x=-x;
				}
				if ((y<halfN) && (x>-halfN)  && ( x<halfN)){
					int j= (x+this.maxN)%this.maxN + this.maxN*y;
					for (int n=0;n<debugdPHdXdY.length;n++) {
						debugdPHdXdY[n][i]=(other && ((n==1)))? (-dPHdXdY[n][j]):dPHdXdY[n][j];
					}
				}

			}
			String [] titles={"amp","phase","weight","d/dx","d/dy"};
			SDFA_INSTANCE.showArrays(
					debugdPHdXdY,
					this.maxN,
					this.maxN,
					true,
					"dbgDiffPhases",
					titles);
		}
		int numAngles=this.polarPhaseMasks.length;
		double [] sinAngle=new double[numAngles];
		double [] cosAngle=new double[numAngles];
		double [] S0= new double [numAngles];
		double [] SF= new double [numAngles];
		double [] SF2=new double [numAngles];
		double [] strength=new double [numAngles];
		double [] mean=new double [numAngles];
		for (int iAngle=0;iAngle<numAngles;iAngle++){
			double angle=iAngle*Math.PI/numAngles;
			sinAngle[iAngle]=Math.sin(angle);
			cosAngle[iAngle]=Math.cos(angle);
			S0[iAngle]=0.0;
			SF[iAngle]=0.0;
			SF2[iAngle]=0.0;
			strength[iAngle]=0.0;
			mean[iAngle]=0.0;
		}
		for (int index=0;index<this.polarPhaseIndices.length;index++) if (dPHdXdY[2][index]>0.0){
			if (this.debug){
				int x= ((index+halfN)%this.maxN)-halfN;
				int y= (index/this.maxN);
				System.out.println("\n+++ calc_dPHdXdY()x/y/ |"+x+"|"+y);
			}
			for (int iiAngle=0;iiAngle<this.polarPhaseIndices[index].length;iiAngle++){
				int iAngle=this.polarPhaseIndices[index][iiAngle];
				double w=this.polarPhaseMasks[iAngle][index]*dPHdXdY[2][index];
				S0[iAngle]+=w;
				if (this.debug){
					System.out.println("    calc_dPHdXdY()  iAngle/mask/weight/w: |"+"|"+iAngle+"|"+
							this.polarPhaseMasks[iAngle][index]+"|"+dPHdXdY[2][index]+"|"+w+"|"+S0[iAngle]);
				}
				double F=dPHdXdY[3][index]*cosAngle[iAngle]+dPHdXdY[4][index]*sinAngle[iAngle];
				SF[iAngle]+=w*F;
				SF2[iAngle]+=w*F*F;

			}
		}
		int iAngleMax=-1;
		double maxStrength=0.0;
		for (int iAngle=0;iAngle<numAngles;iAngle++) {
			double disp=0;
			if (S0[iAngle]>0.0){
				mean[iAngle]=SF[iAngle]/S0[iAngle];
				disp=Math.sqrt(S0[iAngle]*SF2[iAngle]-SF[iAngle]*SF[iAngle])/S0[iAngle];
// TODO: reduce dispersion influence if it is small

				strength[iAngle]=(dispersionFatZero>0.0)?(dispersionFatZero* S0[iAngle]/(dispersionFatZero+disp)):S0[iAngle];
			}
			if (this.debug) System.out.println("calc_dPHdXdY(): |"+iAngle+"|"+
					S0[iAngle]+"|"+SF[iAngle]+"|"+SF2[iAngle]+"|"+mean[iAngle]+"|"+disp+"|"+strength[iAngle]);
			if (strength[iAngle]>maxStrength){
				maxStrength=strength[iAngle];
				iAngleMax=iAngle;
			}
		}
		if (iAngleMax<0){
			if (this.debug) System.out.println("calcPhaseApproximation(): Could not find maximal strength direction");
			return null;
		}
		double b=0.5*(strength[(iAngleMax+1)%numAngles]-strength[(iAngleMax+numAngles-1)%numAngles]);
		double a=0.5*(strength[(iAngleMax+1)%numAngles]+strength[(iAngleMax+numAngles-1)%numAngles])-strength[iAngleMax]; // negative
//		double x=-0.5*b/a;
		double dAngle=(-0.5*b/a);
		double bestAngle=iAngleMax+dAngle;
		maxStrength+=a*dAngle*dAngle+b*dAngle;
		if (this.debug) System.out.println("calcPhaseApproximation(): |"+iAngleMax+"|"+dAngle+"|"+bestAngle+"|"+maxStrength);
		int secondAngle= (iAngleMax+((dAngle<0)?-1:1)+numAngles)%numAngles;
		double secondMean=mean[secondAngle];
		if ((dAngle*(secondAngle-iAngleMax))<0) secondMean=-secondMean; // roll over
		double bestMean=mean[iAngleMax]*(1.0-Math.abs(dAngle))+ secondMean*Math.abs(dAngle);
		if (bestMean<0){
			bestAngle+=numAngles;
			bestMean=-bestMean;
		}
		bestAngle*=Math.PI/numAngles;
		bestAngle-=(2*Math.PI)*Math.round(bestAngle/(2*Math.PI));
		double distance=bestMean*this.maxN/(2*Math.PI);

		if (this.debug) System.out.println("calcPhaseApproximation(): bestAngle="+bestAngle+" ("+IJ.d2s(180*bestAngle/Math.PI,2)+"), distance="+distance+" (bestMean="+bestMean+"), maxStrength="+maxStrength);
		// bin weights, Re and Im for the points with non-zero amplitude and within the band in the direction angle along the direction angle
		double [][] wReIm=new double [3][halfN]; // weight/re/im
		double sin=Math.sin(bestAngle);
		double cos=Math.cos(bestAngle);
		double halfWidth=this.polarRibbonWidth/2;
		boolean conj=(sin<0);
		if (conj){
			sin=-sin;
			cos=-cos;
		}

		for (int index=0;index<dPHdXdY[0].length;index++){
			int indexMod=((this.maxN - (index/this.maxN)) % this.maxN) * this.maxN + ((this.maxN - (index%this.maxN)) % this.maxN);
		    int x=(index+halfN)%this.maxN-halfN;
		    int y=index/this.maxN;
		    double t= cos*x+sin*y;
		    double d= sin*x-cos*y;
		    int it=(int) Math.round(Math.abs(t));
		    if (((y>0) || (x>0)) &&(Math.abs(d) < halfWidth) && (it < halfN) ){
				double re=0.5*(fht[index]+fht[indexMod]);
				double im=0.5*(fht[index]-fht[indexMod]);
				if (conj) im=-im;
				if (t<0){
					im=-im;
					t=-t;
				}
		    	double w=dPHdXdY[0][index]*(0.54+0.46*Math.cos(Math.PI*d/halfWidth)); // (masked) amplitude;
		    	if (t<halfWidth)   w*=0.54+0.46*Math.cos(Math.PI*(halfWidth-t)/halfWidth);
		    	else if (t>(halfN-halfWidth)) w*=0.54+0.46*Math.cos(Math.PI*(t-(halfN-halfWidth))/halfWidth);
		    	wReIm[0][it]+=w;
		    	wReIm[1][it]+=w*re;
		    	wReIm[2][it]+=w*im;
		    }
		}
		if (this.debug) {
			System.out.println("calcPhaseApproximation(): binnig w/re/im along selected direction");
			for (int it=0;it<wReIm[0].length;it++){
				System.out.println("calcPhaseApproximation(): |"+it+"|"+wReIm[0][it]+"|"+wReIm[1][it]+"|"+wReIm[2][it]+"|"+Math.atan2(wReIm[2][it],wReIm[1][it]));
			}
		}
		double [] approximatedZeroPhaseDist=approximateLinearPhase(
				wReIm,
				filter);
		if (approximatedZeroPhaseDist==null) return null;
		if (this.debug) {
			System.out.println("calcPhaseApproximation(): initial distance="+distance+" updated distance="+approximatedZeroPhaseDist[1]+
					" phase at zero="+approximatedZeroPhaseDist[0]+" ("+IJ.d2s(approximatedZeroPhaseDist[0]*180.0/Math.PI,1)+")");
			for (int it=0;it<wReIm[0].length;it++){
				System.out.println("calcPhaseApproximation(): |"+it+"|"+wReIm[0][it]+"|"+wReIm[1][it]+"|"+wReIm[2][it]+"|"+Math.atan2(wReIm[2][it],wReIm[1][it]));
			}
		}
		distance=approximatedZeroPhaseDist[1];
		double phaseAtZero=approximatedZeroPhaseDist[0];
		double [] result={bestAngle,distance,phaseAtZero,maxStrength};
		return result;
	}



	private double polarRibbonWidth=-1;
	private double [][] polarPhaseMasks=null;
	private int [][] polarPhaseIndices=null;
	public void  createPolarMasksIndices(
			double ribbonWidth){ // <0 - invalidate
		if (this.polarRibbonWidth==ribbonWidth) return; // already set;
		this.polarRibbonWidth=ribbonWidth;
		if (ribbonWidth<0.0){
			this.polarPhaseMasks=null;
			this.polarPhaseIndices=null;
			return;
		}
		int halfN=this.maxN/2;
		int length=halfN*this.maxN;
		double halfWidth=ribbonWidth/2;

		int numAngles=(int) Math.round(Math.PI*halfN);
		this.polarPhaseMasks=new double [numAngles][length];
		this.polarPhaseIndices=new int [length][];
		for (int i=0;i<this.polarPhaseIndices.length;i++){
			this.polarPhaseIndices[i]=null;
			for (int iAngle=0;iAngle<numAngles;iAngle++) this.polarPhaseMasks[iAngle][i]=0.0;
		}
		for (int iAngle=0;iAngle<numAngles;iAngle++){
			double angle=iAngle*Math.PI/numAngles;
			double sin=Math.sin(angle);
			double cos=Math.cos(angle);
			double sumWeights=0.0;
			for (int index=0;index<length;index++){
			    int x=(index+halfN)%this.maxN-halfN;
			    int y=index/this.maxN;
			    double t= cos*x+sin*y;
			    double d= sin*x-cos*y;
			    if (((y>0) || (x>0)) &&(Math.abs(d) < halfWidth) && (Math.abs(t) < halfN) ){
			    	this.polarPhaseMasks[iAngle][index]=0.54+0.46*Math.cos(d/halfWidth*Math.PI); // normalize
			    	if (Math.abs(t) > (halfN-halfWidth)) this.polarPhaseMasks[iAngle][index]*=0.54+0.46*Math.cos((Math.abs(t)-(halfN-halfWidth))/halfWidth*Math.PI);
			    	sumWeights+=this.polarPhaseMasks[iAngle][index];
			    }
			}
			sumWeights=1.0/sumWeights;
			for (int index=0;index<length;index++) this.polarPhaseMasks[iAngle][index]*=sumWeights;
		}
		for (int index=0;index<length;index++){
			int numUsedAngles=0;
			for (int iAngle=0;iAngle<numAngles;iAngle++) if (this.polarPhaseMasks[iAngle][index]>0.0) numUsedAngles++;
			this.polarPhaseIndices[index]=new int[numUsedAngles]; // may be 0
			int iiAngle=0;
			for (int iAngle=0;iAngle<numAngles;iAngle++) if (this.polarPhaseMasks[iAngle][index]>0.0) this.polarPhaseIndices[index][iiAngle++]=iAngle;
		}

		if (this.debug){
			double [][]  debugPolarPhaseMasks=new double [numAngles+1][this.maxN*this.maxN];
			for (int i=0;i<debugPolarPhaseMasks[0].length;i++){
				int x= (i%this.maxN)-halfN;
				int y= (i/this.maxN)-halfN;
				if (y<0) {
					y=-y;
					x=-x;
				}
				if ((y<halfN) && (x>-halfN)  && ( x<halfN)){
					int j= (x+this.maxN)%this.maxN + this.maxN*y;
					for (int n=0;n<numAngles;n++) {
						debugPolarPhaseMasks[n][i]=this.polarPhaseMasks[n][j];
					}
					debugPolarPhaseMasks[numAngles][i]=this.polarPhaseIndices[j].length;
				}

			}
			String [] titles=new String[numAngles+1];
			for (int n=0;n<numAngles;n++) titles[n]=IJ.d2s(180.0*n/numAngles,1);
			titles[numAngles]="usage number";
			SDFA_INSTANCE.showArrays(
					debugPolarPhaseMasks,
					this.maxN,
					this.maxN,
					true,
					"dbgDiffPhases",
					titles);
		}
	}

	/**
	 * Project and accumulate Re, Im on the direction perpendicular to the linear feature
	 * @param directionAngle direction  of the perpendicular line
	 * @param phaseIntegrationWidth Disregard pixels outside of the band centered along directionAngle
	 * @param h FHT array
	 * @param nonZeroIndices List of the above-threshold indices in the FHT array
	 * @return arrays of {weight, re, Im} along the line with the step equal to 1 pixel
	 */
	public double [][] binReIm(
			double directionAngle,
			double phaseIntegrationWidth, // zero - no limit
			double [] h,
			int [] nonZeroIndices
			){
		int halfN=this.maxN/2;
		double [][] wReIm=new double [3][halfN]; // weight/re/im
		for (int n=0;n<wReIm.length;n++) for (int i=0;i<wReIm[0].length;i++) wReIm[n][i]=0.0;
		double sin=Math.sin(directionAngle);
		double cos=Math.cos(directionAngle);
		double amp,re,im;
		double halfWidth=(phaseIntegrationWidth>0.0)?phaseIntegrationWidth/2:this.maxN;
		for (int numPoint=0; numPoint<nonZeroIndices.length;numPoint++) {
			int index=nonZeroIndices[numPoint];
			int indexMod=((this.maxN - (index/this.maxN)) % this.maxN) * this.maxN + ((this.maxN - (index%this.maxN)) % this.maxN);
			re=0.5*(h[index]+h[indexMod]);
			im=0.5*(h[index]-h[indexMod]);
		    int x=(index+halfN)%this.maxN-halfN;
		    int y=index/this.maxN;
			amp=Math.sqrt(re*re+im*im);
			if (y==0) amp*=0.5;
			double t = cos*x+sin*y;
			if (Math.abs(cos*y-sin*x)<=halfWidth) {
				int itl= (int) Math.floor(t);
				int ith=itl+1;
				double wl=amp*(t-itl);
				double wh=amp*(itl+1-t);
				boolean conj=(itl<0);
				if (conj){
					itl=-itl;
					ith=-ith;
				}
				if (itl<wReIm[0].length){
					wReIm[0][itl]+=wl;
					wReIm[1][itl]+=wl*re;
					if (conj) wReIm[2][itl]-=wl*im;
					else      wReIm[2][itl]+=wl*im;
				}
				if (ith<wReIm[0].length){
					wReIm[0][ith]+=wh;
					wReIm[1][ith]+=wh*re;
					if (conj) wReIm[2][ith]-=wh*im;
					else      wReIm[2][ith]+=wh*im;
				}
			} else {
				System.out.println("binReIm(): dropped point x="+x+" y="+y+" as it is too far from the center line");
			}
		}
		for (int i=0;i<wReIm[0].length;i++){
			if (wReIm[0][i]>0.0){
				wReIm[1][i]/=wReIm[0][i];
				wReIm[2][i]/=wReIm[0][i];
			}
		}
		return wReIm;
	}
	/**
	 * Linear approximate phase, optionally filtering types of features (black line, white line black-to-white edge, white-to-black edge
	 * @param wReIm - array of weights, Re and Im along the selected direction
	 * @param filter Bitmask enabling different phase half-steps (+1 - 0, +2 - pi/2, +4 - pi, +8 - 3pi/2. Value zero allows arbitrary step
	 * step 0 corresponds to thin white line, pi - thin black line, +pi/2 - edge black-to-white in the direction of directionAngle,
	 * 3pi/2 - white-to-black in the direction of directionAngle
	 * @return pair of phase at 0 and slope expressed in distance (in pixels) from the center to the linear feature along the selected
	 *  direction (may be negative)
	 */
	public double [] approximateLinearPhase(
			double [][]wReIm,
			int filter){
		double [] phase=new double [wReIm[0].length];
		for (int i=0;i<phase.length;i++) if (wReIm[0][i]>0.0) phase[i]=Math.atan2(wReIm[2][i],wReIm[1][i]);
		double diffPhase=0.0;
		double sumWeights=0.0;
		for (int i=0;i<(phase.length-1);i++){
			double w=wReIm[0][i]*wReIm[0][i+1];
			if (w>0.0){
				double diff=phase[i+1]-phase[i];
				if (diff>Math.PI) while (diff>Math.PI) diff-=2*Math.PI;
				else if (diff<-Math.PI) while (diff<-Math.PI) diff+=2*Math.PI;
				sumWeights+=w;
				diffPhase+=diff*w;
			}
		}
		if (sumWeights==0.0) return null;
		diffPhase/=sumWeights;
		double [] binPhases = {0.0,Math.PI/2,Math.PI,-Math.PI/2};
		double [] zeroPhaseBins={0.0,0.0,0.0,0.0};
		double zeroPhase=Double.NaN;
		switch (filter){
		case 0:
			double re=0.0,im=0.0;
			for (int n=0;n<phase.length;n++) if (wReIm[0][n]>0.0) {
				double phZero=phase[n]-diffPhase*n;
				re+=wReIm[0][n]*Math.cos(phZero);
				im+=wReIm[0][n]*Math.sin(phZero);
			}
			zeroPhase=Math.atan2(im,re);
			if (Double.isNaN(zeroPhase)) zeroPhase=0.0; //??
			break;
		case 1:zeroPhase=binPhases[0];  break;
		case 2:zeroPhase=binPhases[1];  break;
		case 4:zeroPhase=binPhases[2];  break;
		case 8:zeroPhase=binPhases[3];  break;
		default:
//			double [] zeroPhaseBins={0.0,0.0,0.0,0.0};
			boolean [] binEnabled={(filter & 1)!=0,(filter & 2)!=0,(filter & 4)!=0,(filter & 8)!=0};
			for (int n=0;n<phase.length;n++) if (wReIm[0][n]>0.0) {
				double phZero=phase[n]-diffPhase*n;
				phZero-=2*Math.PI*Math.floor(phZero/(2*Math.PI)); // now in the range 0..2*PI
				double bestBinDif=2*Math.PI;
				int iBestBin=0;

				for (int i=0;i<binEnabled.length;i++) if (binEnabled[i]){
					double d=phZero-binPhases[i];
					d-=2*Math.PI*Math.round(d/(2*Math.PI));
					d=Math.abs(d);
					if (d<bestBinDif){
						bestBinDif=d;
						iBestBin=i;
					}
				}
				zeroPhaseBins[iBestBin]+=wReIm[0][n];
			}
			double maxBinValue=0;
			int iMaxBin=-1;
			for (int i=0;i<zeroPhaseBins.length;i++) if (zeroPhaseBins[i]>maxBinValue){
				maxBinValue=zeroPhaseBins[i];
				iMaxBin=i;
			}
			zeroPhase=iMaxBin*2*Math.PI/zeroPhaseBins.length;
			if (zeroPhase>Math.PI) zeroPhase-=2*Math.PI;
		}
		double [] fullPhase=phase.clone();
		for (int i=0;i<fullPhase.length;i++)if (wReIm[0][i]>0.0){
			double diff=fullPhase[i]-(zeroPhase+diffPhase*i);
			fullPhase[i]-=2*Math.PI*Math.round(diff/(2*Math.PI));
		}
		if (this.debug){
			for (int i=0;i<zeroPhaseBins.length;i++) 	System.out.println("Phase    "+i+" bin: "+zeroPhaseBins[i]);
			System.out.println("Initial estimated distance to line "+(this.maxN*diffPhase/(2*Math.PI))+", zeroPhase="+zeroPhase);

			for (int i=0;i<phase.length;i++) if (wReIm[0][i]>0.0) {
				System.out.println(i+"|"+IJ.d2s(phase[i],3)+"|"+IJ.d2s(fullPhase[i],3)+"|"+IJ.d2s((zeroPhase+diffPhase*i),3));
			}
		}
// re-adjust slope (distance) after zero phase is selected
		double SFX=0.0,SX2=0.0;
		for (int i=0;i<fullPhase.length;i++)if (wReIm[0][i]>0.0){
			SFX+=wReIm[0][i]*i*(fullPhase[i]-zeroPhase);
			SX2+=wReIm[0][i]*i*i;
		}
		diffPhase=SFX/SX2;
		double distance=this.maxN*diffPhase/(2*Math.PI);
		double [] result={zeroPhase,distance};
		if (this.debug){
			System.out.println("Recalculated distance to line "+distance+" (slope="+diffPhase+" rad/pix)");
		}
		return result;
	}


	public double [][] calcIndAmReIm (double [] h, int [] nonZeroIndices){
		double [][] ampReIm=new double[3][nonZeroIndices.length];
		for (int numPoint=0; numPoint<nonZeroIndices.length;numPoint++) {
			int index=nonZeroIndices[numPoint];
			int indexMod=((this.maxN - (index/this.maxN)) % this.maxN) * this.maxN + ((this.maxN - (index%this.maxN)) % this.maxN);
			ampReIm[1][numPoint]=0.5*(h[index]+h[indexMod]);
			ampReIm[2][numPoint]=0.5*(h[index]-h[indexMod]);
			ampReIm[0][numPoint]=Math.sqrt(ampReIm[1][index]*ampReIm[1][index]+ampReIm[2][index]*ampReIm[2][index]);
		}
		return ampReIm;
	}
	public double [][] calcIndAmHPhase (double [] h, int [] nonZeroIndices){
		double re,im;
		double [][] amHPhase=new double[nonZeroIndices.length][2];
		for (int numPoint=0; numPoint<nonZeroIndices.length;numPoint++) {
			int index=nonZeroIndices[numPoint];
			int indexMod=((this.maxN - (index/this.maxN)) % this.maxN) * this.maxN + ((this.maxN - (index%this.maxN)) % this.maxN);
			re=0.5*(h[index]+h[indexMod]);
			im=0.5*(h[index]-h[indexMod]);
			amHPhase[numPoint][0]=Math.sqrt(re*re+im*im);
			if ((index/this.maxN)>0) amHPhase[numPoint][0]*=2; // double weight
			amHPhase[numPoint][1]=Math.atan2(im,re);
			/*
			System.out.println("numPoint="+numPoint+" index="+nonZeroIndices[numPoint]+
					" x="+(nonZeroIndices[numPoint]%this.maxN)+" y="+(nonZeroIndices[numPoint]/this.maxN)+
					" re="+re+" im="+im+" phase="+amHPhase[numPoint][1]);
					*/
		}
		return amHPhase;
	}
	public int updateFullCycles(
			int [] fullCycles, // may be modified
			double px,
			double py,
			int [] nonZeroIndices,
			double [][] amHPhase
			){
		int numChanges=0;
		double e=1E-8;
		double [] linearPhase= tiltedPhase (px, py, nonZeroIndices);
		for (int numPoint=0; numPoint<nonZeroIndices.length;numPoint++) {
//			double delta=amHPhase[numPoint][1]+fullCycles[numPoint]*2*Math.PI -linearPhase[numPoint];
			double delta=linearPhase[numPoint]-(amHPhase[numPoint][1]+fullCycles[numPoint]*2*Math.PI);
			if (Math.abs(delta)>(Math.PI+e)){
				System.out.print("numPoint="+numPoint+" index="+nonZeroIndices[numPoint]+
						" x="+(nonZeroIndices[numPoint]%this.maxN)+" y="+(nonZeroIndices[numPoint]/this.maxN)+
						" fullCycles["+numPoint+"]: "+fullCycles[numPoint]+" -> ");
				fullCycles[numPoint]+=(int)Math.round(delta/(2*Math.PI));
				System.out.println(fullCycles[numPoint]);
				numChanges++;
			}
		}
		return numChanges;
	}
	public double calcWeightedPhaseRMS(
			int [] fullCycles, // may be modified
			double px,
			double py,
			int [] nonZeroIndices,
			double [][] ampPhase
			){
		double [] linearPhase= tiltedPhase (px, py, nonZeroIndices);
		double S0=0.0,S2=0.0;
		for (int numPoint=0; numPoint<ampPhase.length;numPoint++) {

			double delta=ampPhase[numPoint][1]+fullCycles[numPoint]*2*Math.PI -linearPhase[numPoint];
			S0+=ampPhase[numPoint][0];
			S2+=ampPhase[numPoint][0]*delta*delta;
/*
			System.out.println("nonZeroIndices["+numPoint+"]="+nonZeroIndices[numPoint]+
					" x="+(nonZeroIndices[numPoint]%this.maxN)+" y="+(nonZeroIndices[numPoint]/this.maxN)+
					" ampPhase["+numPoint+"][0]="+ampPhase[numPoint][0]+" ampPhase["+numPoint+"][1]="+ampPhase[numPoint][1]+
					" delta="+delta+" S0="+S0+" S2="+S2);
					*/
		}
		double weightedRMS=Math.sqrt(S2/S0);
		return weightedRMS;
	}

	public double [] calcPhaseShift(
			int [] fullCycles, // may be modified
			int [] nonZeroIndices,
			double [][] ampPhase
			){
//		double [] linearPhase= tiltedPhase (px, py, nonZeroIndices);
// minimize weighted RMS with a plane through (0,0)
		int halfN=this.maxN/2;
		double SFX=0.0,SFY=0.0,SX2=0.0,SY2=0.0,SXY=0.0;
		for (int numPoint=0; numPoint<ampPhase.length;numPoint++) {
			int index=nonZeroIndices[numPoint];
		    int x=(index+halfN)%this.maxN-halfN;
		    int y=index/this.maxN;
		    double fXY=ampPhase[numPoint][1]+Math.PI*2*fullCycles[numPoint];
		    double w=ampPhase[numPoint][0];
		    SFX+=w*x*fXY;
		    SFY+=w*y*fXY;
		    SX2+=w*x*x;
		    SY2+=w*y*y;
		    SXY+=w*x*y;
		}
		double denominator=SXY*SXY-SX2*SY2;
		double [] phaseTilt={(SFY*SXY-SFX*SY2)/denominator,(SFX*SXY-SFY*SX2)/denominator};
		return phaseTilt;
	}

	public double  calcPhaseShift(
			double px,
			double py,
			int [] fullCycles, // may be modified
			int [] nonZeroIndices,
			double [][] ampPhase
			){
//		double [] linearPhase= tiltedPhase (px, py, nonZeroIndices);
// minimize weighted RMS with a plane through (0,0)
		int halfN=this.maxN/2;
		double SFX=0.0,SFY=0.0,SX2=0.0,SY2=0.0,SXY=0.0;
		for (int numPoint=0; numPoint<ampPhase.length;numPoint++) {
			int index=nonZeroIndices[numPoint];
		    int x=(index+halfN)%this.maxN-halfN;
		    int y=index/this.maxN;
		    double fXY=ampPhase[numPoint][1]+Math.PI*2*fullCycles[numPoint];
		    double w=ampPhase[numPoint][0];
		    SFX+=w*x*fXY;
		    SFY+=w*y*fXY;
		    SX2+=w*x*x;
		    SY2+=w*y*y;
		    SXY+=w*x*y;
		}
		return (py*SFY+px*SFX-px*py*SXY)/(px*px*SX2+py*py*SY2);
	}


	public double [] tiltedPhase (
			double px,
			double py,
			int [] nonZeroIndices){
		int halfN=this.maxN/2;
		double [] phase=new double [nonZeroIndices.length];
		for (int numPoint=0; numPoint<nonZeroIndices.length;numPoint++) {
			int index=nonZeroIndices[numPoint];
		    int x=(index+halfN)%this.maxN-halfN;
		    int y=index/this.maxN;
		    phase[numPoint]=px*x+py*y;
		}
		return phase;
	}

	/*
	public double [][] calcAmReIm (double [] h, double [] ampMask){
		int rowMod, colMod;
		int halfN=this.maxN/2;
		double [][] ampReIm=new double[3][this.maxN*halfN];
		for (int r =0; r<this.maxN/2; r++) {
			rowMod = (this.maxN - r) % this.maxN;
			for (int c=0; c<this.maxN; c++){
				int index=   r * this.maxN + c;
				if ((ampMask==null) || (ampMask[index]>0.0)){
					colMod = (this.maxN - c) % this.maxN;
					int indexMod=rowMod * this.maxN + colMod;
					ampReIm[1][index]=0.5*(h[index]+h[indexMod]);
					ampReIm[2][index]=0.5*(h[index]-h[indexMod]);
					ampReIm[0][index]=Math.sqrt(ampReIm[1][index]*ampReIm[1][index]+ampReIm[2][index]*ampReIm[2][index]);
				} else {
					ampReIm[1][index]=0.0;
					ampReIm[2][index]=0.0;
					ampReIm[0][index]=0.0;
				}
			}
		}
		return ampReIm;
	}
	*/
	public double [] filterAmplitude(
			double [] h,
			boolean removeIslands,
			double threshold){ // relative to RMS value
		int rowMod, colMod;
		int halfN=this.maxN/2;
		double [] amplitude=new double[this.maxN*halfN];
		double max=0;
		int iMax=-1;
		double sum2=0.0;
		double amp2;
		for (int r =0; r<this.maxN/2; r++) {
			rowMod = (this.maxN - r) % this.maxN;
			for (int c=0; c<this.maxN; c++){
				int index=   r * this.maxN + c;
				if (c==this.maxN/2) {
					amplitude[index]=0.0;
				} else {
					colMod = (this.maxN - c) % this.maxN;
					int indexMod=rowMod * this.maxN + colMod;
					if (r>0) {
						amp2=h[index]*h[index]+h[indexMod]*h[indexMod];  // count twice (for the second symmetrical half)
						amplitude[index]=2.0*amp2; // so square root will be twice
					} else {
						amp2=0.5*(h[index]*h[index]+h[indexMod]*h[indexMod]);
						amplitude[index]=amp2;
					}
					sum2+=amp2;
					if (amp2>max){
						max=amp2;
						iMax=index;
					}
				}
			}
		}
		double minAmp2=threshold*sum2/(this.maxN -1)/(this.maxN -1);
		for (int i=0;i<amplitude.length;i++){
			if (amplitude[i]<minAmp2) amplitude[i]=0.0;
			else amplitude[i]=Math.sqrt(amplitude[i]);
		}
		int numDefined=0;
		if (this.debug) {
			for (int i=0;i<amplitude.length;i++) if (amplitude[i]>0.0) numDefined++;
			System.out.println("quadraticAmplitude(): number of defined cells="+numDefined+" ( of "+amplitude.length+")");
		}

		if (removeIslands) {
			boolean [] confirmed=new boolean [amplitude.length];
			for (int i=0;i<confirmed.length;i++) confirmed[i]=false;
			int [][] dirs={{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}};
		    List <Integer> ampList=new ArrayList<Integer>(confirmed.length);
		    int x=(iMax+halfN)%this.maxN-halfN;
		    int y=iMax/this.maxN;
		    Integer Index=iMax;
		    ampList.add(Index);
		    confirmed[Index]=true;
		    while (ampList.size()>0){
		    	int index=ampList.remove(0);
			    x=(index+halfN)%this.maxN-halfN;
			    y=index/this.maxN;
			    for (int iDir=0;iDir<dirs.length;iDir++){
			    	int x1=x+dirs[iDir][0];
			    	int y1=y+dirs[iDir][1];
			    	if (y1<0){
			    		y1=-y1;
			    		x1=-x1;
			    	}
			    	if ((x1>-halfN) && (x1<halfN) && (y1>=0) && (y1<halfN)){
			    		Index=(x1+this.maxN)%this.maxN + this.maxN*y1;
			    		if (!confirmed[Index] && (amplitude[Index]>0.0)){
			    			confirmed[Index]=true;
			    			ampList.add(Index);
			    		}
			    	}
			    }
		    }
			for (int i=0;i<amplitude.length;i++){
				if (!confirmed[i]) amplitude[i]=0.0;
			}
			if (this.debug) {
				numDefined=0;
				for (int i=0;i<amplitude.length;i++) if (amplitude[i]>0.0) numDefined++;
				System.out.println("quadraticAmplitude(): number of remaining cells="+numDefined+" ( of "+amplitude.length+")");
			}
		}
		/*
		if (this.debug){
			SDFA_INSTANCE.showArrays(
					amplitude,
					this.maxN,
					halfN,
					"fltAmp");
			double [] debugAmplitude=new double [this.maxN*this.maxN];
			for (int i=0;i<debugAmplitude.length;i++){
				debugAmplitude[i]=0.0;
				int x= (i%this.maxN)-halfN;
				int y= (i/this.maxN)-halfN;
				if (y<0) {
					y=-y;
					x=-x;
				}
				if ((y<halfN) && (x>-halfN)  && ( x<halfN)){
					int j= (x+this.maxN)%this.maxN + this.maxN*y;
					debugAmplitude[i]=amplitude[j];
					if (y>0) debugAmplitude[i]*=0.5;
				}

			}
			SDFA_INSTANCE.showArrays(
					debugAmplitude,
					this.maxN,
					this.maxN,
					"fltAmp");
		}
		*/
		return amplitude;
	}
	public double [] quadraticAmplitude(
			double [] amplitude,
			boolean removeIslands,
			double threshold){ // relative to RMS value
		int halfN=this.maxN/2;
		double[] result = new double[3]; // CX2, CY2, CXY
		for (int i=0;i<result.length;i++) result[i]=0.0;
		double s=0.0;
		double amp;
		for (int r =0; r<this.maxN/2; r++) {
			for (int c=0; c<this.maxN; c++) if (c!=this.maxN/2){
				int index=   r * this.maxN + c;
				amp=amplitude[index];
				s+=amp;
				int x=((c+halfN)%this.maxN)-halfN;
				result[0]+=amp*x*x;
				result[1]+=amp*r*r;
				result[2]+=amp*x*r;
			}
		}
		s*=halfN*halfN; // normalize so result will be related to pixels, not to the FFT size
		for (int i=0;i<result.length;i++) result[i]/=s;
		return result;
	}

	/*
	public double [] quadraticAmplitude(
			double [] h){
		int rowMod, colMod;
		double[] result = new double[3]; // CX2, CY2, CXY
		for (int i=0;i<result.length;i++) result[i]=0.0;
		double s=0.0;
		double amp;
		int halfN=this.maxN/2;
		for (int r =0; r<this.maxN/2; r++) {
			rowMod = (this.maxN - r) % this.maxN;
			for (int c=0; c<this.maxN; c++) if (c!=this.maxN/2){
				colMod = (this.maxN - c) % this.maxN;
				int index=   r * this.maxN + c;
				int indexMod=rowMod * this.maxN + colMod;
				amp=Math.sqrt(0.5*(h[index]*h[index]+h[indexMod]*h[indexMod]));
				if (r>0) amp*=2; // count twice (for the second symmetrical half)
				s+=amp;
				int x=((c+halfN)%this.maxN)-halfN;
				result[0]+=amp*x*x;
				result[1]+=amp*r*r;
				result[2]+=amp*x*r;
			}
		}
		s*=halfN*halfN; // normalize so result will be related to pixels, not to the FFT size
		for (int i=0;i<result.length;i++) result[i]/=s;
		return result;
	}
	*/
	public double [] ellipseAmplitude(
			double [] quadratic){
		double[] result = new double[3]; // angle (clockwise from X, large half-diameter, small/large ratio
		double alpha=0.5*Math.atan2(2*quadratic[2], quadratic[0]-quadratic[1]); // 1/2*atan(2*Sxy/(sx2-sy2))
		double cos=Math.cos(alpha), sin=Math.sin(alpha);
		double I1=sin*sin*quadratic[0]+cos*cos*quadratic[1]-2*cos*sin*quadratic[2];
		double I2=cos*cos*quadratic[0]+sin*sin*quadratic[1]+2*cos*sin*quadratic[2];
		boolean rot=(I1>I2);
		result[0]=alpha+(rot?Math.PI:0);
		result[1]=Math.sqrt(rot?I1:I2);
		result[2]=Math.sqrt(rot?(I2/I1):(I1/I2));
		return result;
	}

	public double [] multiply(double [] h1, double [] h2, boolean  conjugate) {
		int rowMod, colMod;
		double h2e, h2o;
		double[] product = new double[maxN*maxN];
		for (int r =0; r<maxN; r++) {
			rowMod = (maxN - r) % maxN;
			for (int c=0; c<maxN; c++) {
				colMod = (maxN - c) % maxN;
				h2e = (h2[r * maxN + c] + h2[rowMod * maxN + colMod]) / 2;
				h2o = (h2[r * maxN + c] - h2[rowMod * maxN + colMod]) / 2;
				if (conjugate)
					product[r * maxN + c] = h1[r * maxN + c] * h2e - h1[rowMod * maxN + colMod] * h2o;
				else
					product[r * maxN + c] = h1[r * maxN + c] * h2e + h1[rowMod * maxN + colMod] * h2o;
			}
		}
		return product;
	}

	public double [] phaseMultiply(double [] h1, double [] h2, double phaseCoeff) {
		int rowMod, colMod;
		double h2e, h2o,d;
		double[] product = new double[maxN*maxN];
		for (int r =0; r<maxN; r++) {
			rowMod = (maxN - r) % maxN;
			for (int c=0; c<maxN; c++) {
				colMod = (maxN - c) % maxN;
				h2e = (h2[r * maxN + c] + h2[rowMod * maxN + colMod]) / 2;
				h2o = (h2[r * maxN + c] - h2[rowMod * maxN + colMod]) / 2;
				d=phaseCoeff*(h2e*h2e+h2o*h2o)+(1.0-phaseCoeff);
				product[r * maxN + c] = (h1[r * maxN + c] * h2e - h1[rowMod * maxN + colMod] * h2o)/d;
			}
		}
		return product;
	}


// Multiply by real array (i.e. filtering in frequency domain). Array m length should be at least maxN*maxN/2+1
	public void multiplyByReal(double [] h, double [] m) {
		int rowMod, colMod;
		int index=0, indexM;
		for (int r =0; r<maxN; r++) {
			rowMod = (maxN - r) % maxN;
			for (int c=0; c<maxN; c++) {
				colMod = (maxN - c) % maxN;
                indexM=rowMod * maxN + colMod;
                if (indexM>index) indexM=index;
                h[index]*=m[indexM];
				index++;
			}
		}
	}


/* nothing to prevent division by small values */
	/* Returns the image resulting from the point by point Hartley division
		of this image by the specified image. Both images are assumed to be in
		the frequency domain. Division in the frequency domain is equivalent
		to deconvolution in the space domain. */
	public double [] divide(double [] h1, double [] h2) {
		int rowMod, colMod;
		double mag, h2e, h2o;
		double[] result = new double[maxN*maxN];
		for (int r=0; r<maxN; r++) {
			rowMod = (maxN - r) % maxN;
			for (int c=0; c<maxN; c++) {
				colMod = (maxN - c) % maxN;
				mag =h2[r*maxN+c] * h2[r*maxN+c] + h2[rowMod*maxN+colMod] * h2[rowMod*maxN+colMod];
				if (mag<1e-20)
					mag = 1e-20;
				h2e = (h2[r*maxN+c] + h2[rowMod*maxN+colMod]);
				h2o = (h2[r*maxN+c] - h2[rowMod*maxN+colMod]);
				double tmp = (h1[r*maxN+c] * h2e - h1[rowMod*maxN+colMod] * h2o);
				result[r*maxN+c] = tmp/mag;
			}
		}
		return result;
	}


  public double [] calculateAmplitude(double [] fht) {
    int size=(int) Math.sqrt(fht.length);
    double[] amp = new double[size*size];
    for (int row=0; row<size; row++) {
      amplitude(row, size, fht, amp);
    }
    swapQuadrants(amp);
    return amp;
  }
  public double [] calculateAmplitude2(double [] fht) {
	    int size=(int) Math.sqrt(fht.length);
	    double[] amp = new double[size*size];
	    for (int row=0; row<size; row++) {
	      amplitude2(row, size, fht, amp);
	    }
	    swapQuadrants(amp);
	    return amp;
	  }
/* Amplitude of one row from 2D Hartley Transform. */
  void amplitude(int row, int size, double[] fht, double[] amplitude) {
    int base = row*size;
    int l;
    for (int c=0; c<size; c++) {
        l = ((size-row)%size) * size + (size-c)%size;
        amplitude[base+c] = Math.sqrt(fht[base+c]*fht[base+c] + fht[l]*fht[l]);
    }
  }

  /* Squared amplitude of one row from 2D Hartley Transform. */
  void amplitude2(int row, int size, double[] fht, double[] amplitude) {
    int base = row*size;
    int l;
    for (int c=0; c<size; c++) {
        l = ((size-row)%size) * size + (size-c)%size;
        amplitude[base+c] = fht[base+c]*fht[base+c] + fht[l]*fht[l];
    }
  }
// Other FHT-related methods moved here
  /* converts FHT results (frequency space) to complex numbers of [fftsize/2+1][fftsize] */
  public double[][][] FHT2FFTHalf (FHT fht, int fftsize) {
     float[] fht_pixels=(float[])fht.getPixels();
     double[][][] fftHalf=new double[(fftsize>>1)+1][fftsize][2];
     int row1,row2,col1,col2;

     for (row1=0;row1<=(fftsize>>1);row1++) {
       row2=(fftsize-row1) %fftsize;
       for (col1=0;col1<fftsize;col1++) {
         col2=(fftsize-col1) %fftsize;
//         fftHalf[row1][col1]=   complex( 0.5*(fht_pixels[row1*fftsize+col1] + fht_pixels[row2*fftsize+col2]),
//                                         0.5*(fht_pixels[row2*fftsize+col2] - fht_pixels[row1*fftsize+col1]));
         fftHalf[row1][col1][0]=   0.5*(fht_pixels[row1*fftsize+col1] + fht_pixels[row2*fftsize+col2]);
         fftHalf[row1][col1][1]=   0.5*(fht_pixels[row2*fftsize+col2] - fht_pixels[row1*fftsize+col1]);
       }
     }
     return fftHalf;
  }

  public double[][][] FHT2FFTHalf (double [] fht_pixels, int fftsize) {
     double[][][] fftHalf=new double[(fftsize>>1)+1][fftsize][2];
     int row1,row2,col1,col2;

     for (row1=0;row1<=(fftsize>>1);row1++) {
       row2=(fftsize-row1) %fftsize;
       for (col1=0;col1<fftsize;col1++) {
         col2=(fftsize-col1) %fftsize;
         fftHalf[row1][col1][0]=   0.5*(fht_pixels[row1*fftsize+col1] + fht_pixels[row2*fftsize+col2]);
         fftHalf[row1][col1][1]=   0.5*(fht_pixels[row2*fftsize+col2] - fht_pixels[row1*fftsize+col1]);
       }
     }
     return fftHalf;
  }

  /* converts FFT arrays of complex numbers of [fftsize/2+1][fftsize] to FHT arrays */
  public float[] floatFFTHalf2FHT (double [][][] fft, int fftsize) {
     float[] fht_pixels=new float [fftsize*fftsize];
     int row1,row2,col1,col2;
     for (row1=0;row1<=(fftsize>>1);row1++) {
       row2=(fftsize-row1) %fftsize;
       for (col1=0;col1 < fftsize;col1++) {
         col2=(fftsize-col1) %fftsize;
         fht_pixels[row1*fftsize+col1]=(float)(fft[row1][col1][0]-fft[row1][col1][1]);
         fht_pixels[row2*fftsize+col2]=(float)(fft[row1][col1][0]+fft[row1][col1][1]);
       }
     }
     return fht_pixels;
  }

  public double[] FFTHalf2FHT (double [][][] fft, int fftsize) {
     double[] fht_pixels=new double [fftsize*fftsize];
     int row1,row2,col1,col2;
     for (row1=0;row1<=(fftsize>>1);row1++) {
       row2=(fftsize-row1) %fftsize;
       for (col1=0;col1 < fftsize;col1++) {
         col2=(fftsize-col1) %fftsize;
         fht_pixels[row1*fftsize+col1]=fft[row1][col1][0]-fft[row1][col1][1];
         fht_pixels[row2*fftsize+col2]=fft[row1][col1][0]+fft[row1][col1][1];
       }
     }
     return fht_pixels;
  }

  // Uses just first half and one line
  public double[] FFTHalf2FHT (double [][]fft, int fftsize) {
	     double[] fht_pixels=new double [fftsize*fftsize];
	     int row1,row2,col1,col2;
	     for (row1=0;row1<=(fftsize>>1);row1++) {
	       row2=(fftsize-row1) %fftsize;
	       for (col1=0;col1 < fftsize;col1++) {
	         col2=(fftsize-col1) %fftsize;
	         fht_pixels[row1*fftsize+col1]=fft[0][row1*fftsize+col1]-fft[1][row1*fftsize+col1];
	         fht_pixels[row2*fftsize+col2]=fft[0][row1*fftsize+col1]+fft[1][row1*fftsize+col1];
	       }
	     }
	     return fht_pixels;
	  }



  /* Amplitude/phase related methods */
  public double [] interpolateFHT (double [] fht0,    // first FHT array
		  double [] fht1,    // second FHT array
		  double   ratio){   // array of interpolation points - 0.0 - fht0, 1.0 - fht1
	  double [] points={ratio};
	  double [][] results= interpolateFHT (fht0,    // first FHT array
			  fht1,    // second FHT array
			  points,  // array of interpolation points - 0.0 - fht0, 1.0 - fht1
			  true);   // do not clone 0.0 and 1.0 (ends)
	  return results[0];
  }

  /* returns array of interpolated FHTs between fht0 and fht1, endpoints if present (0.0, 1.0) are referenced, not cloned */
  public double [][] interpolateFHT (
		  double []   fht0,    // first FHT array
		  double []   fht1,    // second FHT array
		  double [] points){   // array of interpolation points - 0.0 - fht0, 1.0 - fht1
	  return interpolateFHT (fht0,    // first FHT array
			  fht1,    // second FHT array
			  points,   // array of interpolation points - 0.0 - fht0, 1.0 - fht1
			  false);   // do not clone 0.0 and 1.0 (ends)
  }
  public double [][] interpolateFHT (double [] fht0,    // first FHT array
                                     double [] fht1,    // second FHT array
                                   double [] points,   // array of interpolation points - 0.0 - fht0, 1.0 - fht1
                                boolean cloneTrivial   ){
     double [] fht_div=divide(fht1,fht0);
     int size=(int) Math.sqrt(fht0.length);
     int hsize=size/2;
     double [][][] aphase= new double[hsize+1][size][2];
     double [][][] amp01=  new double[hsize+1][size][2]; /* squared amplitudes of fht0 and fht1 */
     double [][]   phase=  new double[hsize+1][size]; /* +/-pi phase of the first array */
     double[][][]fft0=    FHT2FFTHalf (fht0, size);
     double[][][]fft1=    FHT2FFTHalf (fht1, size);
     double[][][]fft_div= FHT2FFTHalf (fht_div, size);
     int i,j,k;
     double a,c,p;
 /* use mul for amplitudes, div - for phases */
     for (i=0;i<=hsize;i++) for (j=0;j<size;j++) {
       amp01[i][j][0]= fft0[i][j][0]*fft0[i][j][0]+fft0[i][j][1]*fft0[i][j][1];
       amp01[i][j][1]= fft1[i][j][0]*fft1[i][j][0]+fft1[i][j][1]*fft1[i][j][1];
       if (amp01[i][j][0]>0.0) phase[i][j]=Math.atan2(fft0[i][j][1],fft0[i][j][0]);
       else phase[i][j]=0.0;
       aphase[i][j][0]=amp01[i][j][0]*amp01[i][j][1]; // product of squared amplitudes is OK for phase restorations, we just need to know where amplitudes are higher
       a=              fft_div[i][j][0]*fft_div[i][j][0]+fft_div[i][j][1]*fft_div[i][j][1];
       if (a>0.0) aphase[i][j][1]=Math.atan2(fft_div[i][j][1],fft_div[i][j][0]);
       else aphase[i][j][1]=0.0;
     }
     aphase[0][0][1]=0.0;
 /* calculate full phases */
     fullPhase(aphase);
     double [][]result=new double[points.length][];
     for (k=0;k<result.length;k++) {
        if      (points[k]==0.0) result[k]=cloneTrivial?fht0.clone():fht0;
        else if (points[k]==1.0) result[k]=cloneTrivial?fht1.clone():fht1;
        else { /* interpolate */
          c=points[k];
          for (i=0;i<=hsize;i++) for (j=0;j<size;j++) {
            if ((amp01[i][j][0]==0.0) || (amp01[i][j][1]==0.0)) a=0.0;
 /* Extrapolation is defined here only in the direction of decreasing of the spectral amplitudes (outside, to the wider PSF), so additional limit to prevent division of small values */
 /* Seems to work, possible improvements: 1-filter spectrum in high-freq areas. 2 - use farther inner points for farther approximation */

            else if ((c<0.0) && (amp01[i][j][0]>amp01[i][j][1])) a=Math.sqrt(amp01[i][j][0]);
            else if ((c>1.0) && (amp01[i][j][0]<amp01[i][j][1])) a=Math.sqrt(amp01[i][j][1]);
            else a= Math.pow(amp01[i][j][0],0.5*(1.0-c))*Math.pow(amp01[i][j][1],0.5*c);
            p= phase[i][j]+c*aphase[i][j][1];
            fft0[i][j][0]=a*Math.cos(p);
            fft0[i][j][1]=a*Math.sin(p);
          }
          result[k]=FFTHalf2FHT(fft0,size);
        }
     }
     return result;
  }
  /**
   *
   * @param fht FHT data
   * @param fullCentered - if false - returns minimal (size*(size/2-1) arrays, if true - fills the full square arrays and swaps quadrants
   * @return 2-d array, first index: 0 - amplitude, 1 - phase. second index range  depends on fullCentered
   */
  public double [][] fht2AmpHase(double [] fht, boolean fullCentered){
	     int size=(int) Math.sqrt(fht.length);
	     int hsize=size/2;
	     double[][][]fft=    FHT2FFTHalf (fht, size);
	     double [][]aphase=new double [2][size*(fullCentered?size:(hsize+1))];
	     int index=0;
	     double a;
	     for (int i=0;i<=hsize;i++) for (int j=0;j<size;j++){
	    	 a=Math.sqrt(fft[i][j][0]*fft[i][j][0]+fft[i][j][1]*fft[i][j][1]);
	    	 aphase[0][index]=a;
	    	 aphase[1][index++]=(a>0.0)?Math.atan2(fft[i][j][1],fft[i][j][0]):0.0;
	     }
	     if (fullCentered) {
		     for (int i=1;i<hsize;i++) for (int j=0;j<size;j++){
		    	 int j1=(size-j)%size;
		    	 aphase[0][(size-i)*size+j]= aphase[0][i*size+j1];
		    	 aphase[1][(size-i)*size+j]=-aphase[1][i*size+j1];
		     }
		     swapQuadrants(aphase[0]);
		     swapQuadrants(aphase[1]);
	     }
	     return aphase;
  }

  public double [][] fht2ReIm(double [] fht, boolean fullCentered){
	     int size=(int) Math.sqrt(fht.length);
	     int hsize=size/2;
	     double[][][]fft=    FHT2FFTHalf (fht, size);
	     double [][]reIm=new double [2][size*(fullCentered?size:(hsize+1))];
	     int index=0;
	     for (int i=0;i<=hsize;i++) for (int j=0;j<size;j++){
	    	 reIm[0][index]=fft[i][j][0];
	    	 reIm[1][index++]=fft[i][j][1];
	     }
	     if (fullCentered) {
		     for (int i=1;i<hsize;i++) for (int j=0;j<size;j++){
		    	 int j1=(size-j)%size;
		    	 reIm[0][(size-i)*size+j]= reIm[0][i*size+j1];
		    	 reIm[1][(size-i)*size+j]=-reIm[1][i*size+j1];
		     }
		     swapQuadrants(reIm[0]);
		     swapQuadrants(reIm[1]);
	     }
	     return reIm;
}

  /* replace +/-pi  phase with the full phase, using amplitude to guide grouth of the covered area, so amplitude and phase does not need to be a pair from the same FFT array */
  public void fullPhase(double [][][] aphase) {
    int size = aphase[0].length;
    int hsize=aphase.length-1;
    boolean [][]  map=    new boolean[hsize+1][size];
    int i,j;
    aphase[0][0][1]=0.0;
    int ix,iy,ix1,iy1,ix1n,iy1n;
    List <Integer> pixelList=new ArrayList<Integer>(100);
    Integer Index;
    int [][] dirs={{-1,0},{-1,-1},{0,-1},{1,-1},{1,0},{1,1},{0,1},{-1,1}};
    ix=0;
    iy=0;
    int clusterSize=0;
    boolean noNew=true;
    int maxX=0;
    int maxY=0;
    int oldX=0;
    int oldY=0;
    boolean oldConj=false;
    int listIndex;
    Index=iy*size + ix;
    pixelList.clear();
    pixelList.add (Index);
    clusterSize++;
    map[iy][ix]=true;
    noNew=true;
    double phase, oldPhase, fullCyclesPhase;
    double maxValue=-1.0;

    while (pixelList.size()>0) {
  /* Find maximal new neighbor */
      maxValue=-1.0;
      listIndex=0;
      while (listIndex<pixelList.size()) {
        Index=pixelList.get(listIndex);
        iy=Index/size;
        ix=Index%size;
        noNew=true;
        for (j=0;j<8;j++) {
          ix1=(ix+dirs[j][0]+size) % size;
          iy1=(iy+dirs[j][1]+size) % size;
          if ((iy1>hsize) || (((iy1==0) || (iy1==hsize)) && (ix1> hsize))) {
            ix1n=(size-ix1)%size;
            iy1n=(size-iy1)%size;
          } else { /* But phase will be opposite sign */
            ix1n=ix1;
            iy1n=iy1;
          }
          if (!map[iy1n][ix1n]) {
            noNew=false;
            if (aphase[iy1n][ix1n][0]>maxValue) {
              maxValue= aphase[iy1n][ix1n][0];
              maxX=ix1n;
              maxY=iy1n;
//      if (DEBUG_LEVEL>4)      System.out.println(" amplPhase(): iy="+iy+ " ix="+ix+" maxY="+maxY+" maxX="+maxX);
            }
          }
        }
        if (noNew) pixelList.remove(listIndex);  //  remove current list element
        else       listIndex++;     // increase list index
      }


      if (pixelList.size()==0) break;
  /* To calculate the phase - find already processed neighbor with the highest amplitude */
      maxValue=-1.0;
      for (j=0;j<8;j++) {
        ix1=(maxX+dirs[j][0]+size) % size;
        iy1=(maxY+dirs[j][1]+size) % size;
        if ((iy1>hsize) || (((iy1==0) || (iy1==hsize)) && (ix1> hsize))) {
          ix1n=(size-ix1)%size;
          iy1n=(size-iy1)%size;
        } else { /* But phase will be opposite sign */
          ix1n=ix1;
          iy1n=iy1;
        }
        if (map[iy1n][ix1n]) {
          if (aphase[iy1n][ix1n][0]>maxValue) {
            maxValue= aphase[iy1n][ix1n][0];
            oldX=ix1n;
            oldY=iy1n;
            oldConj=(iy1!=iy1n) || (ix1!=ix1n); // point on the other half (conjugate)

          }
        }
      }

  /* Calculate the phase from the closest neighbor */
      oldPhase=(oldConj?-1:1)*aphase[oldY][oldX][1];
      fullCyclesPhase=2*Math.PI*Math.floor(oldPhase/(2*Math.PI)+0.5);
      oldPhase-=fullCyclesPhase; // +/- pi
      phase=aphase[maxY][maxX][1];
      if      ((phase - oldPhase) > Math.PI) fullCyclesPhase-=2*Math.PI;
      else if ((oldPhase - phase) > Math.PI) fullCyclesPhase+=2*Math.PI;
      aphase[maxY][maxX][1]+=fullCyclesPhase;
/*      if (DEBUG_LEVEL>3) {
        System.out.println(" amplPhase():Old:["+oldConj+"] ("+pixelList.size()+")  "+oldX+":"+oldY+" "+IJ.d2s(aphase[oldY][oldX][0],2)+":"+ IJ.d2s(aphase[oldY][oldX][1],2)+
                                       " New:"+maxX+":"+maxY+" "+IJ.d2s(aphase[maxY][maxX][0],2)+":"+ IJ.d2s(aphase[maxY][maxX][1],2)+
                                       " Diff="+IJ.d2s((aphase[maxY][maxX][1]-(oldConj?-1:1)*aphase[oldY][oldX][1]),2));
      }
*/
  /* Add this new point to the list */
      Index=maxY*size + maxX;
      pixelList.add (Index);
      clusterSize++;
      map[maxY][maxX]=true;
    } // end of while (pixelList.size()>0)
  /* Fix remaining phases for y=0 and y=hsize */
    for (i=1;i<hsize;i++) {
      aphase[0][size-i][1]=-aphase[0][i][1];
      aphase[hsize][size-i][1]=-aphase[hsize][i][1];
    }
  }





}
