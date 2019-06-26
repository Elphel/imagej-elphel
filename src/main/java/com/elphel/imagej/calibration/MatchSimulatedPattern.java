package com.elphel.imagej.calibration;
/**
 **
 ** MatchSimulatedPattern.java - Determine simulation pattern parameters to match
 ** the acquired image
 **
 ** Copyright (C) 2010-2014 Elphel, Inc.
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

import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.List;
import java.util.Properties;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import javax.swing.SwingUtilities;

import com.elphel.imagej.common.DoubleFHT;
import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.PolynomialApproximation;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.jp4.JP46_Reader_camera;
import com.elphel.imagej.lwir.LwirReaderParameters;

import Jama.LUDecomposition;
import Jama.Matrix;  // Download here: http://math.nist.gov/javanumerics/jama/
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.process.FHT; // get rid, change to double
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;


public class MatchSimulatedPattern {
	private ShowDoubleFloatArrays SDFA_INSTANCE= new ShowDoubleFloatArrays(); // just for debugging?
	public int debugLevel=2;
	public int FFT_SIZE=256;
	public double [][][][] PATTERN_GRID=null; // global to be used with threads? TODO: Same as DIST_ARRAY - merge?
	public int [] minUV={0,0};
	public int [][] reMap=null;               // maps grid coordinates from laser pointers to PATTERN_GRID u,v (2x3 - rotation + translation) - SEEMS NOT USED!
	public int [] UVShiftRot={0,0,0}; // {shift U, shift V, rot (1 of 8 - 4 non-mirrored, 4 - mirrored}
	public int [][][] targetUV=null; // maps PATTERN_GRID cells to the target (absolute) UV
	//    public double [][][] pixelsUV=null; // made of PATTERN_GRID, but does not have any wave vectors. Calculated during laser calibration
	public double [][][] pXYUV=null; // made of PATTERN_GRID, but does not have any wave vectors. Calculated during laser calibration
	public double [][][] gridContrastBrightness=null; //{grid contrast, grid intensity red, grid intensity green, grid intensity blue}[v][u]
	public Rectangle DIST_SELECTION=null;
	public int [] UV_INDEX=null; // array containing index of the pattern UV (scanline order, U first), or -1 for the areas with no pattern
	public int  UV_INDEX_WIDTH=0; // width of UV_INDEX (== full image, not selection, width)
	public int [] debugUV={-1,-1}; // debug the same cell on the second pass
	public int passNumber=0;
	public boolean [] correlationSizesUsed=null;
	public double [] gridFFCorr=null; // array matching greens with the flat field correction for the grid (zero outside of detected grid?)
	public double [] flatFieldForGrid=null; // array matching image pixels, divide the input pixels by these values (if not null)
	public boolean [] focusMask=null; // array matching image pixels, used with focusing (false outside sample areas)
	final private static int [][][] rotations={
			{{ 1, 0},{ 0, 1}}, // not mirrored
			{{ 0, 1},{-1, 0}},
			{{-1, 0},{ 0,-1}},
			{{ 0,-1},{ 1, 0}},

			{{ 1, 0},{ 0,-1}}, // mirrored
			{{ 0, 1},{ 1, 0}},
			{{-1, 0},{ 0, 1}},
			{{ 0,-1},{-1, 0}}};
	// shifts when rotating around unknown center (make it white)
	final private static int [][] dfltShifts={
			{0,0},
			{0,1},
			{0,0},
			{1,0},
			{0,1},
			{0,0},
			{1,0},
			{0,0}};
	final private static int [][] combinedRotations={
			{0,	1,	2,	3,	4,	5,	6,	7},
			{1,	2,	3,	0,	7,	4,	5,	6},
			{2,	3,	0,	1,	6,	7,	4,	5},
			{3,	0,	1,	2,	5,	6,	7,	4},
			{4,	5,	6,	7,	0,	1,	2,	3},
			{5,	6,	7,	4,	3,	0,	1,	2},
			{6,	7,	4,	5,	2,	3,	0,	1},
			{7,	4,	5,	6,	1,	2,	3,	0}};

	public MatchSimulatedPattern (){ }
	public MatchSimulatedPattern (int fft_size){
		this.FFT_SIZE=fft_size;
	}
	// not real clone, just for threads - if there will be FFT - keep individual
	@Override
	public MatchSimulatedPattern clone(){ // used in createPSFMap when creating threads
		MatchSimulatedPattern msp=new MatchSimulatedPattern (this.FFT_SIZE);
		// cloning should be thread safe, when using DoubleFHT - use individual instances
		msp.debugLevel=this.debugLevel;
		msp.PATTERN_GRID=this.PATTERN_GRID; // global to be used with threads? TODO: Same as DIST_ARRAY - merge?
		//    	msp.DIST_ARRAY=this.DIST_ARRAY;
		msp.DIST_SELECTION=this.DIST_SELECTION;
		msp.UV_INDEX=this.UV_INDEX; // array containing index of the pattern UV (scanline order, U first), or -1 for the areas with no pattern
		msp.UV_INDEX_WIDTH=this.UV_INDEX_WIDTH;
		msp.reMap=this.reMap;
		msp.UVShiftRot=this.UVShiftRot.clone();
		//    public int [] UVShiftRot={0,0,0}; // {shift U, shift V, rot (1 of 8 - 4 non-mirrored, 4 - mirrored}

		msp.targetUV=this.targetUV; //
		msp.pXYUV=this.pXYUV;
		msp.flatFieldForGrid=this.flatFieldForGrid;
		msp.focusMask=this.focusMask;
		return msp;
	}

	public MatchSimulatedPattern cloneDeep(
			boolean clonePATTERN_GRID,
			boolean cloneTargetUV,
			boolean clonePixelsUV,
			boolean cloneFlatFieldForGrid,
			boolean cloneFocusMask
			){ // used in createPSFMap when creating threads
		MatchSimulatedPattern msp=new MatchSimulatedPattern (this.FFT_SIZE);
		// cloning should be thread safe, when using DoubleFHT - use individual instances
		msp.debugLevel=this.debugLevel;
		//    	msp.PATTERN_GRID=this.PATTERN_GRID; // global to be used with threads? TODO: Same as DIST_ARRAY - merge?
		if ( clonePATTERN_GRID && (msp.PATTERN_GRID!=null)){
			msp.PATTERN_GRID=new double [this.PATTERN_GRID.length][this.PATTERN_GRID[0].length][][];
			for (int i=0;i<this.PATTERN_GRID.length;i++) for (int j=0;j<this.PATTERN_GRID[i].length;j++) {
				if (this.PATTERN_GRID[i][j]!=null){
					msp.PATTERN_GRID[i][j]=new double [this.PATTERN_GRID[i][j].length][];
					for (int k=0;k<this.PATTERN_GRID[i][j].length;k++){
						if (this.PATTERN_GRID[i][j][k]!=null) msp.PATTERN_GRID[i][j][k]=this.PATTERN_GRID[i][j][k].clone();
						else msp.PATTERN_GRID[i][j][k]=null;
					}
				} else msp.PATTERN_GRID[i][j]=null;
			}
		} else msp.PATTERN_GRID=this.PATTERN_GRID;
		//    	msp.DIST_ARRAY=this.DIST_ARRAY;
		if (this.DIST_SELECTION!=null) msp.DIST_SELECTION=new Rectangle(this.DIST_SELECTION);
		else msp.DIST_SELECTION=null;

		if (this.UV_INDEX!=null) msp.UV_INDEX=this.UV_INDEX.clone(); // array containing index of the pattern UV (scanline order, U first), or -1 for the areas with no pattern
		else msp.UV_INDEX=null;

		msp.UV_INDEX_WIDTH=this.UV_INDEX_WIDTH;

		if (this.reMap!=null) { // probably not used
			msp.reMap=new int [2][];
			msp.reMap[0]=this.reMap[0].clone();
			msp.reMap[1]=this.reMap[1].clone();
		} else msp.reMap=null;

		msp.UVShiftRot=this.UVShiftRot.clone();
		//     public int [] UVShiftRot={0,0,0}; // {shift U, shift V, rot (1 of 8 - 4 non-mirrored, 4 - mirrored}
		//    	msp.targetUV=this.targetUV;
		if (cloneTargetUV && (this.targetUV != null)) {
			msp.targetUV=new int [this.targetUV.length][this.targetUV[0].length][];
			for (int i=0;i<this.targetUV.length;i++) for (int j=0;j<this.targetUV[i].length;j++){
				if (this.targetUV[i][j]!=null) msp.targetUV[i][j]=this.targetUV[i][j].clone();
				else msp.targetUV[i][j]=null;
			}
		} else msp.targetUV=this.targetUV;

		//    	msp.pixelsUV=this.pixelsUV;

		if (clonePixelsUV && (this.pXYUV != null)) {
			msp.pXYUV=new double [this.pXYUV.length][this.pXYUV[0].length][];
			for (int i=0;i<this.pXYUV.length;i++) for (int j=0;j<this.pXYUV[i].length;j++){
				if (this.pXYUV[i][j]!=null) msp.pXYUV[i][j]=this.pXYUV[i][j].clone();
				else msp.pXYUV[i][j]=null;
			}
		} else msp.targetUV=this.targetUV;

		//    	msp.flatFieldForGrid=this.flatFieldForGrid;
		if (cloneFlatFieldForGrid && (this.flatFieldForGrid!=null)) msp.flatFieldForGrid=this.flatFieldForGrid.clone();
		else msp.flatFieldForGrid=this.flatFieldForGrid;

		//		msp.focusMask=this.focusMask;
		if (cloneFocusMask && (this.focusMask!=null)) msp.focusMask=this.focusMask.clone();
		else msp.focusMask=this.focusMask;
		return msp;
	}

	public int [] getUVShiftRot(boolean shift){
		if (!shift) return this.UVShiftRot;
		int [][] reReMap=getRemapMatrix(this.UVShiftRot);
		int [] UVShiftRotCorr=this.UVShiftRot.clone();
		UVShiftRotCorr[0]-=reReMap[0][0]*this.minUV[0]+reReMap[0][1]*this.minUV[1];
		UVShiftRotCorr[1]-=reReMap[1][0]*this.minUV[0]+reReMap[1][1]*this.minUV[1];
		System.out.println("getUVShiftRot(true): minUV[0]="+minUV[0]+" minUV[1]="+minUV[1]);
		return UVShiftRotCorr;
	}

	public static int [][] getRemapMatrix(
			int [] UVShiftRot){
		// Moved shift calculation to calibrateGrid(), here just a regular R,T 2x3 matix
		int [][] reReMap={
				{rotations[UVShiftRot[2]][0][0],rotations[UVShiftRot[2]][0][1], UVShiftRot[0]},
				{rotations[UVShiftRot[2]][1][0],rotations[UVShiftRot[2]][1][1], UVShiftRot[1]}};
		return reReMap;
	}

	public static int [] combineUVShiftRot(
			int [] UVShiftRotA,
			int [] UVShiftRotB){
		int [][] reReMapA=getRemapMatrix(UVShiftRotA);
		int [][] reReMapB=getRemapMatrix(UVShiftRotB);
		int [] UVShiftRot=new int [3];
		UVShiftRot[0] = reReMapB[0][0]*reReMapA[0][2]+ reReMapB[0][1]*reReMapA[1][2] + reReMapB[0][2];
		UVShiftRot[1] = reReMapB[1][0]*reReMapA[0][2]+ reReMapB[1][1]*reReMapA[1][2] + reReMapB[1][2];
		UVShiftRot[2] = combinedRotations[UVShiftRotA[2]][UVShiftRotB[2]];
		return UVShiftRot;
	}

	public double focusQualityR2(
			ImagePlus imp,
			int fftSize,
			int spotSize, //5
			double x0,
			double y0,
			int debugLevel
			){
		int ix0=(int) Math.round(x0);
		int iy0=(int) Math.round(y0);
		Rectangle selection=new Rectangle(ix0-fftSize,iy0-fftSize,fftSize*2,fftSize*2);
		double [][] pixels=splitBayer(imp, selection,true);
		//    	SDFA_INSTANCE.showArrays(pixels, fftSize, fftSize, true,"bayer");
		DoubleFHT fht_instance =new DoubleFHT(); // provide DoubleFHT instance to save on initializations (or null)
		double [] hamming1d=fht_instance.getHamming1d(fftSize);
		for (int c=0;c<pixels.length;c++) if (pixels[c]!=null){
			double sum=0.0;
			for (int i=0;i<pixels[c].length;i++) sum+=pixels[c][i];
			double average=sum/pixels[c].length;
			int index=0;
			for (int y=0;y<fftSize;y++)for (int x=0;x<fftSize;x++){
				pixels[c][index]=(pixels[c][index]-average)*hamming1d[y]*hamming1d[x];
				index++;
			}
		}
		// slightly blur the image to be sure there are no aliases in the low frequencies
		double preBlurSigma=1.5;
		DoubleGaussianBlur gb=new DoubleGaussianBlur();
		for (int c=0;c<pixels.length;c++) if (pixels[c]!=null){
			gb.blurDouble(
					pixels[c],
					fftSize,
					fftSize,
					preBlurSigma,
					preBlurSigma,
					0.01);


		}
		//    	SDFA_INSTANCE.showArrays(pixels, fftSize, fftSize, true,"bayer-winowed");
		for (int c=0;c<pixels.length;c++) if (pixels[c]!=null){
			fht_instance.swapQuadrants(pixels[c]);
			fht_instance.transform(pixels[c]);
			pixels[c]=fht_instance.calculateAmplitude(pixels[c]);
		}
		if (debugLevel>1) SDFA_INSTANCE.showArrays(pixels, fftSize, fftSize, true,"amplitudes");
		int [][][] spectrumMaximums=new int [pixels.length][][];
		double [] spectralContrast=new double [pixels.length];
		int hsize=fftSize/2;
		int radius=spotSize/2; // 2
		int [][]dirs={
				{-1, 0},
				{-1, 1},
				{ 0, 1},
				{ 1, 1},
				{ 1, 0},
				{ 1,-1},
				{ 0,-1},
				{-1,-1}};
		int lowLim= (fftSize*3)/8;
		int highLim=(fftSize*5)/8;
		for (int c=0;c<pixels.length;c++) if (pixels[c]!=null){
			spectrumMaximums[c]=new int [2][2];
			double max=0.0;
			for (int y=lowLim;y<=fftSize/2;y++)for (int x=lowLim;x<highLim;x++){
				if ((y>=(hsize-radius)) && (x>=(hsize-radius)) && (x<=(hsize+radius))) continue; // do not count zero freq
				if (max<pixels[c][y*fftSize+x]) {
					max=pixels[c][y*fftSize+x];
					spectrumMaximums[c][0][0]=x;
					spectrumMaximums[c][0][1]=y;
				}
			}
			max=0.0;
			for (int y=lowLim;y<=fftSize/2;y++)for (int x=lowLim;x<highLim;x++){
				if ((y>=(hsize-radius)) && (x>=(hsize-radius)) && (x<=(hsize+radius))) {
					if ((debugLevel>1) && (c==0))System.out.println("x]="+x+" y="+y+" - too close to 0,0");
					continue; // do not count zero freq
				}
				if ((y>=(spectrumMaximums[c][0][1]-radius)) && (y<=(spectrumMaximums[c][0][1]+radius)) &&
						(x>=(spectrumMaximums[c][0][0]-radius)) && (x<=(spectrumMaximums[c][0][0]+radius))){
					if ((debugLevel>1) && (c==0))System.out.println("x="+x+" y="+y+" - too close to first max "+spectrumMaximums[c][0][0]+":"+spectrumMaximums[c][0][1]);
					continue; // do not count zero freq
				}
				if ((y>=((fftSize-spectrumMaximums[c][0][1])-radius)) && (y<=((fftSize-spectrumMaximums[c][0][1])+radius)) &&
						(x>=((fftSize-spectrumMaximums[c][0][0])-radius)) && (x<=((fftSize-spectrumMaximums[c][0][0])+radius))){
					if ((debugLevel>1) && (c==0))System.out.println("x="+x+" y="+y+" - too close to alias "+(fftSize-spectrumMaximums[c][0][0])+":"+(fftSize-spectrumMaximums[c][0][1]));
					continue; // do not count zero first maximum
				}
				if (max<pixels[c][y*fftSize+x]) {
					max=pixels[c][y*fftSize+x];
					spectrumMaximums[c][1][0]=x;
					spectrumMaximums[c][1][1]=y;
					if ((debugLevel>1) && (c==0))System.out.println("New max at x="+x+" y="+y+": "+max);
				}
			}
			if (debugLevel>1) System.out.println("spectrumMaximums["+c+"]="+(spectrumMaximums[c][0][0]-hsize)+":"+(spectrumMaximums[c][0][1]-hsize)+", "+
					""+(spectrumMaximums[c][1][0]-hsize)+":"+(spectrumMaximums[c][1][1]-hsize));
			//			double s1=0,s2=0;
			int [][] xy3={
					{3*spectrumMaximums[c][0][0]-fftSize,3*spectrumMaximums[c][0][1]-fftSize},
					{3*spectrumMaximums[c][1][0]-fftSize,3*spectrumMaximums[c][1][1]-fftSize}};
			//			if (debugLevel>3){
			if ((c==0) && (debugLevel>1)){
				System.out.println(" xy3[0][0]="+xy3[0][0]+" xy3[0][1]="+xy3[0][1]);
				System.out.println(" xy3[1][0]="+xy3[1][0]+" xy3[1][1]="+xy3[1][1]);
			}
			// make 3-rd  harmonic frequency adjustment
			// May be improved by quadratic maximums
			for (int n=0;n<2;n++) for (int i=0;i<2;i++){
				max=pixels[c][xy3[i][1]*fftSize+xy3[i][0]];
				int d=-1;
				for (int dir=0;dir<dirs.length;dir++) {
					double v=pixels[c][(xy3[i][1]+dirs[dir][1])*fftSize+(xy3[i][0]+dirs[dir][0])];
					if (max<v){
						max=v;
						d=dir;
					}
				}
				if (d>=0) {
					xy3[i][0]+=dirs[d][0];
					xy3[i][1]+=dirs[d][1];
				}
			}
			//			if (debugLevel>2) {
			if ((c==0) && (debugLevel>1)) {
				System.out.println("*xy3[0][0]="+xy3[0][0]+" xy3[0][1]="+xy3[0][1]);
				System.out.println("*xy3[1][0]="+xy3[1][0]+" xy3[1][1]="+xy3[1][1]);
			}
			double [][] xy={
					{(xy3[0][0]-hsize)/3.0,(xy3[0][1]-hsize)/3.0},
					{(xy3[1][0]-hsize)/3.0,(xy3[1][1]-hsize)/3.0}};
			if ((c==0) && (debugLevel>1)) {
				System.out.println(" xy[0][0]="+IJ.d2s(xy[0][0],1)+" xy3[0][1]="+IJ.d2s(xy[0][1],1));
				System.out.println(" xy[1][0]="+IJ.d2s(xy[1][0],1)+" xy3[1][1]="+IJ.d2s(xy[1][1],1));
			}
			// once more - correct 5-th:
			int [][] xy5={
					{hsize+ (int) Math.round(5*xy[0][0]),hsize+ (int) Math.round(5*xy[0][1])},
					{hsize+ (int) Math.round(5*xy[1][0]),hsize+ (int) Math.round(5*xy[1][1])}};
			// just once

			if ((c==0) && (debugLevel>1)) {
				System.out.println(" xy5[0][0]="+xy5[0][0]+" xy5[0][1]="+xy5[0][1]);
				System.out.println(" xy5[1][0]="+xy5[1][0]+" xy5[1][1]="+xy5[1][1]);
			}
			for (int i=0;i<2;i++){
				max=pixels[c][xy5[i][1]*fftSize+xy5[i][0]];
				int d=-1;
				for (int dir=0;dir<dirs.length;dir++) {
					double v=pixels[c][(xy5[i][1]+dirs[dir][1])*fftSize+(xy5[i][0]+dirs[dir][0])];
					if (max<v){
						max=v;
						d=dir;
					}
				}
				if (d>=0) {
					xy5[i][0]+=dirs[d][0];
					xy5[i][1]+=dirs[d][1];
				}
			}
			if ((c==0) && (debugLevel>1)) {
				System.out.println("*xy5[0][0]="+xy5[0][0]+" xy5[0][1]="+xy5[0][1]);
				System.out.println("*xy5[1][0]="+xy5[1][0]+" xy5[1][1]="+xy5[1][1]);
			}
			for (int i=0;i<2;i++) for (int j=0;j<2;j++) xy[i][j]=(xy5[i][j]-hsize)/5.0;
			if ((c==0) && (debugLevel>1)) {
				System.out.println(" xy[0][0]="+IJ.d2s(xy[0][0],1)+" xy3[0][1]="+IJ.d2s(xy[0][1],1));
				System.out.println(" xy[1][0]="+IJ.d2s(xy[1][0],1)+" xy3[1][1]="+IJ.d2s(xy[1][1],1));
			}


			// now generate mask and then blur it
			int maxMode=7; // 1, 3-rd and 5-th harmonics, should be odd (now 7-th included)
			double [][] dirNeg={{-0.5,-0.5},{0.5,-0.5},{-0.5,0.5},{0.5,0.5}}; // direction to negative mask
			double [] mask = new double [fftSize*fftSize];
			for (int i=0;i<mask.length;i++) mask[i]=0.0;
			for (int i=-maxMode;i<=maxMode;i+=2) for (int j=-maxMode;j<=maxMode;j+=2){
				double xpc=-0.5*(xy[0][0]*(i+j)+xy[1][0]*(i-j));
				double ypc=-0.5*(xy[0][1]*(i+j)+xy[1][1]*(i-j));
				double xp=xpc+hsize;
				double yp=ypc+hsize;
				double w=xpc*xpc+ypc*ypc; // increase weight of high frequencies
				int ixp=(int) Math.round(xp);
				int iyp=(int) Math.round(yp);
				if ((ixp<0) ||(ixp>=fftSize) || (iyp<0) ||(iyp>=fftSize)) continue;
				mask[iyp*fftSize+ixp]+=w;
				if ((c==0) && (debugLevel>1)) 	System.out.println(" xp="+IJ.d2s(xp,1)+" xm="+IJ.d2s(yp,1));
				for (int d=0;d<dirNeg.length;d++){
					double xm=xp+dirNeg[d][0]*xy[0][0]+dirNeg[d][1]*xy[1][0];
					double ym=yp+dirNeg[d][0]*xy[0][1]+dirNeg[d][1]*xy[1][1];
					int ixm=(int) Math.round(xm);
					int iym=(int) Math.round(ym);
					ixm=(ixm+fftSize)%fftSize;
					iym=(iym+fftSize)%fftSize;
					mask[iym*fftSize+ixm]-=w/dirNeg.length;
				}
			}
			double sigmaScale=0.2;
			double averageLength=Math.sqrt((xy[0][0]*xy[0][0] + xy[0][1]*xy[0][1] +xy[1][0]*xy[1][0] +xy[1][1]*xy[1][1])/2);
			if ((c==0) && (debugLevel>1)) SDFA_INSTANCE.showArrays(mask, fftSize, fftSize, "mask_color");
			gb.blurDouble(
					mask,
					fftSize,
					fftSize,
					sigmaScale*averageLength,
					sigmaScale*averageLength,
					0.01);
			if ((c==0) && (debugLevel>1)) SDFA_INSTANCE.showArrays(mask, fftSize, fftSize, "mask_color_blured"+IJ.d2s(sigmaScale*averageLength,3));
			double SFR2=0.0,SFP=0.0;
			for (int i=0;i<mask.length;i++){
				int x=(i % fftSize)-hsize;
				int y=(i / fftSize)-hsize;
				SFR2+=pixels[c][i]*mask[i]*(x*x+y*y);
				if (mask[i]>0) SFP+=pixels[c][i]*mask[i]; // sum only positive masks? Or abs value? Does not really matter, it is just a scale
			}
			spectralContrast[c]=Math.sqrt(SFR2/SFP)/averageLength;
			if ((c==0) && (debugLevel>1)) System.out.println("SFR2="+SFR2+" SFP="+SFP+" averageLength="+averageLength);
			if ( (debugLevel>1)) System.out.println("spectrumContrast["+c+"]="+spectralContrast[c]);

		} else{
			spectrumMaximums[c]=null;
			spectralContrast[c]=Double.NaN;
		}
		// 	    return 0.25*(spectralContrast[0]+spectralContrast[1]+spectralContrast[2]+spectralContrast[3]);
		return 0.5*(spectralContrast[0]+spectralContrast[3]); // green only
	}
	public double focusQuality(
			ImagePlus imp,
			int fftSize,
			int spotSize, //5
			double x0,
			double y0,
			int debugLevel
			){
		int ix0=(int) Math.round(x0);
		int iy0=(int) Math.round(y0);
		Rectangle selection=new Rectangle(ix0-fftSize,iy0-fftSize,fftSize*2,fftSize*2);
		double [][] pixels=splitBayer(imp, selection,true);
		//    	SDFA_INSTANCE.showArrays(pixels, fftSize, fftSize, true,"bayer");
		DoubleFHT fht_instance =new DoubleFHT(); // provide DoubleFHT instance to save on initializations (or null)
		double [] hamming1d=fht_instance.getHamming1d(fftSize);
		for (int c=0;c<pixels.length;c++) if (pixels[c]!=null){
			double sum=0.0;
			for (int i=0;i<pixels[c].length;i++) sum+=pixels[c][i];
			double average=sum/pixels[c].length;
			int index=0;
			for (int y=0;y<fftSize;y++)for (int x=0;x<fftSize;x++){
				pixels[c][index]=(pixels[c][index]-average)*hamming1d[y]*hamming1d[x];
				index++;
			}
		}
		// slightly blur the image to be sure there are no aliases in the low frequencies
		double preBlurSigma=1.5;
		DoubleGaussianBlur gb=new DoubleGaussianBlur();
		for (int c=0;c<pixels.length;c++) if (pixels[c]!=null){
			gb.blurDouble(
					pixels[c],
					fftSize,
					fftSize,
					preBlurSigma,
					preBlurSigma,
					0.01);


		}
		//    	SDFA_INSTANCE.showArrays(pixels, fftSize, fftSize, true,"bayer-winowed");
		for (int c=0;c<pixels.length;c++) if (pixels[c]!=null){
			fht_instance.swapQuadrants(pixels[c]);
			fht_instance.transform(pixels[c]);
			pixels[c]=fht_instance.calculateAmplitude(pixels[c]);
		}
		if (debugLevel>1) SDFA_INSTANCE.showArrays(pixels, fftSize, fftSize, true,"amplitudes");
		int [][][] spectrumMaximums=new int [pixels.length][][];
		double [] spectralContrast=new double [pixels.length];
		int hsize=fftSize/2;
		int radius=spotSize/2; // 2
		int [][]dirs={
				{-1, 0},
				{-1, 1},
				{ 0, 1},
				{ 1, 1},
				{ 1, 0},
				{ 1,-1},
				{ 0,-1},
				{-1,-1}};
		int lowLim= (fftSize*3)/8;
		int highLim=(fftSize*5)/8;
		for (int c=0;c<pixels.length;c++) if (pixels[c]!=null){
			spectrumMaximums[c]=new int [2][2];
			double max=0.0;
			for (int y=lowLim;y<=fftSize/2;y++)for (int x=lowLim;x<highLim;x++){
				if ((y>=(hsize-radius)) && (x>=(hsize-radius)) && (x<=(hsize+radius))) continue; // do not count zero freq
				if (max<pixels[c][y*fftSize+x]) {
					max=pixels[c][y*fftSize+x];
					spectrumMaximums[c][0][0]=x;
					spectrumMaximums[c][0][1]=y;
				}
			}
			max=0.0;
			for (int y=lowLim;y<=fftSize/2;y++)for (int x=lowLim;x<highLim;x++){
				if ((y>=(hsize-radius)) && (x>=(hsize-radius)) && (x<=(hsize+radius))) {
					if ((debugLevel>1) && (c==0))System.out.println("x]="+x+" y="+y+" - too close to 0,0");
					continue; // do not count zero freq
				}
				if ((y>=(spectrumMaximums[c][0][1]-radius)) && (y<=(spectrumMaximums[c][0][1]+radius)) &&
						(x>=(spectrumMaximums[c][0][0]-radius)) && (x<=(spectrumMaximums[c][0][0]+radius))){
					if ((debugLevel>1) && (c==0))System.out.println("x="+x+" y="+y+" - too close to first max "+spectrumMaximums[c][0][0]+":"+spectrumMaximums[c][0][1]);
					continue; // do not count zero freq
				}
				if ((y>=((fftSize-spectrumMaximums[c][0][1])-radius)) && (y<=((fftSize-spectrumMaximums[c][0][1])+radius)) &&
						(x>=((fftSize-spectrumMaximums[c][0][0])-radius)) && (x<=((fftSize-spectrumMaximums[c][0][0])+radius))){
					if ((debugLevel>1) && (c==0))System.out.println("x="+x+" y="+y+" - too close to alias "+(fftSize-spectrumMaximums[c][0][0])+":"+(fftSize-spectrumMaximums[c][0][1]));
					continue; // do not count zero first maximum
				}
				if (max<pixels[c][y*fftSize+x]) {
					max=pixels[c][y*fftSize+x];
					spectrumMaximums[c][1][0]=x;
					spectrumMaximums[c][1][1]=y;
					if ((debugLevel>1) && (c==0))System.out.println("New max at x="+x+" y="+y+": "+max);
				}
			}
			if (debugLevel>1) System.out.println("spectrumMaximums["+c+"]="+(spectrumMaximums[c][0][0]-hsize)+":"+(spectrumMaximums[c][0][1]-hsize)+", "+
					""+(spectrumMaximums[c][1][0]-hsize)+":"+(spectrumMaximums[c][1][1]-hsize));
			//			double s1=0,s2=0;
			int [][] xy3={
					{3*spectrumMaximums[c][0][0]-fftSize,3*spectrumMaximums[c][0][1]-fftSize},
					{3*spectrumMaximums[c][1][0]-fftSize,3*spectrumMaximums[c][1][1]-fftSize}};
			//			if (debugLevel>3){
			if ((c==0) && (debugLevel>1)){
				System.out.println(" xy3[0][0]="+xy3[0][0]+" xy3[0][1]="+xy3[0][1]);
				System.out.println(" xy3[1][0]="+xy3[1][0]+" xy3[1][1]="+xy3[1][1]);
			}
			// make 3-rd  harmonic frequency adjustment
			// May be improved by quadratic maximums
			for (int n=0;n<2;n++) for (int i=0;i<2;i++){
				max=pixels[c][xy3[i][1]*fftSize+xy3[i][0]];
				int d=-1;
				for (int dir=0;dir<dirs.length;dir++) {
					double v=pixels[c][(xy3[i][1]+dirs[dir][1])*fftSize+(xy3[i][0]+dirs[dir][0])];
					if (max<v){
						max=v;
						d=dir;
					}
				}
				if (d>=0) {
					xy3[i][0]+=dirs[d][0];
					xy3[i][1]+=dirs[d][1];
				}
			}
			//			if (debugLevel>2) {
			if ((c==0) && (debugLevel>1)) {
				System.out.println("*xy3[0][0]="+xy3[0][0]+" xy3[0][1]="+xy3[0][1]);
				System.out.println("*xy3[1][0]="+xy3[1][0]+" xy3[1][1]="+xy3[1][1]);
			}
			double [][] xy={
					{(xy3[0][0]-hsize)/3.0,(xy3[0][1]-hsize)/3.0},
					{(xy3[1][0]-hsize)/3.0,(xy3[1][1]-hsize)/3.0}};
			if ((c==0) && (debugLevel>1)) {
				System.out.println(" xy[0][0]="+IJ.d2s(xy[0][0],1)+" xy3[0][1]="+IJ.d2s(xy[0][1],1));
				System.out.println(" xy[1][0]="+IJ.d2s(xy[1][0],1)+" xy3[1][1]="+IJ.d2s(xy[1][1],1));
			}
			// once more - correct 5-th:
			int [][] xy5={
					{hsize+ (int) Math.round(5*xy[0][0]),hsize+ (int) Math.round(5*xy[0][1])},
					{hsize+ (int) Math.round(5*xy[1][0]),hsize+ (int) Math.round(5*xy[1][1])}};
			// just once

			if ((c==0) && (debugLevel>1)) {
				System.out.println(" xy5[0][0]="+xy5[0][0]+" xy5[0][1]="+xy5[0][1]);
				System.out.println(" xy5[1][0]="+xy5[1][0]+" xy5[1][1]="+xy5[1][1]);
			}
			for (int i=0;i<2;i++){
				max=pixels[c][xy5[i][1]*fftSize+xy5[i][0]];
				int d=-1;
				for (int dir=0;dir<dirs.length;dir++) {
					double v=pixels[c][(xy5[i][1]+dirs[dir][1])*fftSize+(xy5[i][0]+dirs[dir][0])];
					if (max<v){
						max=v;
						d=dir;
					}
				}
				if (d>=0) {
					xy5[i][0]+=dirs[d][0];
					xy5[i][1]+=dirs[d][1];
				}
			}
			if ((c==0) && (debugLevel>1)) {
				System.out.println("*xy5[0][0]="+xy5[0][0]+" xy5[0][1]="+xy5[0][1]);
				System.out.println("*xy5[1][0]="+xy5[1][0]+" xy5[1][1]="+xy5[1][1]);
			}
			for (int i=0;i<2;i++) for (int j=0;j<2;j++) xy[i][j]=(xy5[i][j]-hsize)/5.0;
			if ((c==0) && (debugLevel>1)) {
				System.out.println(" xy[0][0]="+IJ.d2s(xy[0][0],1)+" xy3[0][1]="+IJ.d2s(xy[0][1],1));
				System.out.println(" xy[1][0]="+IJ.d2s(xy[1][0],1)+" xy3[1][1]="+IJ.d2s(xy[1][1],1));
			}


			// now generate mask and then blur it
			int maxMode=7; // 1, 3-rd and 5-th harmonics, should be odd (now 7-th included)
			double [][] dirNeg={{-0.5,-0.5},{0.5,-0.5},{-0.5,0.5},{0.5,0.5}}; // direction to negative mask
			double [] mask = new double [fftSize*fftSize];
			for (int i=0;i<mask.length;i++) mask[i]=0.0;
			for (int i=-maxMode;i<=maxMode;i+=2) for (int j=-maxMode;j<=maxMode;j+=2){
				double xpc=-0.5*(xy[0][0]*(i+j)+xy[1][0]*(i-j));
				double ypc=-0.5*(xy[0][1]*(i+j)+xy[1][1]*(i-j));
				double xp=xpc+hsize;
				double yp=ypc+hsize;
				double w=xpc*xpc+ypc*ypc; // increase weight of high frequencies
				w*=w; // even sharper
				int ixp=(int) Math.round(xp);
				int iyp=(int) Math.round(yp);
				if ((ixp<0) ||(ixp>=fftSize) || (iyp<0) ||(iyp>=fftSize)) continue;
				mask[iyp*fftSize+ixp]+=w;
				if ((c==0) && (debugLevel>1)) 	System.out.println(" xp="+IJ.d2s(xp,1)+" xm="+IJ.d2s(yp,1));
				for (int d=0;d<dirNeg.length;d++){
					double xm=xp+dirNeg[d][0]*xy[0][0]+dirNeg[d][1]*xy[1][0];
					double ym=yp+dirNeg[d][0]*xy[0][1]+dirNeg[d][1]*xy[1][1];
					int ixm=(int) Math.round(xm);
					int iym=(int) Math.round(ym);
					ixm=(ixm+fftSize)%fftSize;
					iym=(iym+fftSize)%fftSize;
					mask[iym*fftSize+ixm]-=w/dirNeg.length;
				}
			}
			double sigmaScale=0.2;
			double averageLength=Math.sqrt((xy[0][0]*xy[0][0] + xy[0][1]*xy[0][1] +xy[1][0]*xy[1][0] +xy[1][1]*xy[1][1])/2);
			if ((c==0) && (debugLevel>1)) SDFA_INSTANCE.showArrays(mask, fftSize, fftSize, "mask_color");
			gb.blurDouble(
					mask,
					fftSize,
					fftSize,
					sigmaScale*averageLength,
					sigmaScale*averageLength,
					0.01);
			if ((c==0) && (debugLevel>1)){
				SDFA_INSTANCE.showArrays(mask, fftSize, fftSize, "mask_color_blured"+IJ.d2s(sigmaScale*averageLength,3));
				double [] ppixels=new double [fftSize*fftSize];
				for (int i=0;i<ppixels.length;i++) ppixels[i]=pixels[c][i]*mask[i];
				SDFA_INSTANCE.showArrays(ppixels, fftSize, fftSize, "masked-amplitude"+IJ.d2s(sigmaScale*averageLength,3));
			}
			double SFM=0.0,SF=0.0,SM=0;
			for (int i=0;i<mask.length;i++){
				//				int x=(i % fftSize)-hsize;
				//				int y=(i / fftSize)-hsize;
				SFM+=pixels[c][i]*mask[i];
				SF+=pixels[c][i];
				if (mask[i]>0) SM+=mask[i]; // sum only positive masks? Or abs value? Does not really matter, it is just a scale

			}
			int S0=fftSize*fftSize;
			spectralContrast[c]=S0*SFM/SF/SM;
			if ((c==0) && (debugLevel>1)) System.out.println("SFM="+(SFM/S0)+" SF="+(SF/S0)+" SM="+(SM/S0)+" averageLength="+averageLength);
			if ( (debugLevel>1)) System.out.println("spectrumContrast["+c+"]="+spectralContrast[c]);

		} else{
			spectrumMaximums[c]=null;
			spectralContrast[c]=Double.NaN;
		}
		// 	    return 0.25*(spectralContrast[0]+spectralContrast[1]+spectralContrast[2]+spectralContrast[3]);
		return 0.5*(spectralContrast[0]+spectralContrast[3]); // green only
	}

	public void resetCorrelationSizesUsed(){
		this.correlationSizesUsed=null;

	}
	public void setCorrelationSizesUsed(int size){
		int maxSize=20;
		if (this.correlationSizesUsed==null) {
			this.correlationSizesUsed=new boolean [maxSize];
			for (int i=0;i<maxSize;i++) this.correlationSizesUsed[i]=false;
		}
		int ln2;
		for (ln2=0;size>(1<<ln2); ln2++);
		if (size!=(1<<ln2)){
			String msg="Not a power of 2 :"+ size;
			IJ.showMessage("Error",msg);
			throw new IllegalArgumentException (msg);
		}
		if (ln2>=maxSize) {
			String msg="Too large array length, increase maxSize (it is now "+maxSize+", wanted "+ln2+")";
			IJ.showMessage("Error",msg);
			throw new IllegalArgumentException (msg);
		}
		this.correlationSizesUsed[ln2]=true;

	}
	public boolean [] getCorrelationSizesUsed(){
		return this.correlationSizesUsed;
	}

	/* returns array of 3 arrays: first two are 3-element wave vectors (x,y,phase), last - 3-rd order correction coefficients */
	public double[][] findPatternDistorted(
			double [][]             bayer_mono_pixels, // pixel array to process (no windowing!), two greens will be used
			PatternDetectParameters patternDetectParameters,
			double                  min_half_period,
			double                  max_half_period,
			boolean                 greens,  // this is a pattern for combined greens (diagonal), adjust results accordingly
			String                  title){ // title prefix to use for debug  images
		if (bayer_mono_pixels==null) return null;
		if (bayer_mono_pixels.length < 1) return null;
		if (bayer_mono_pixels[0]==null) return null;
		if (greens) {
			if (bayer_mono_pixels.length<4) return null;
			if (bayer_mono_pixels[3]==null) return null;
		} else {
			if (bayer_mono_pixels.length > 1) return null;
		}
		boolean is_mono = bayer_mono_pixels.length == 1;
		int tile_size2=bayer_mono_pixels[0].length;
		int tile_size=(int) Math.sqrt(tile_size2);
		int hsize=tile_size/2;
		int hsize2=hsize*hsize;
		double [] quarterHamming=initWindowFunction(hsize, patternDetectParameters.gaussWidth);
		double [] patternCorr=new double[6]; // second order non-linear pattern correction (perspective+distortion)
		double [] green0= new double[hsize2];
		double [] green3= new double[hsize2];

		double [][]   quarter_pixels =new double[9][];
		double [][][] quarter_patterns =new double[9][][];
		int [] quarterIndex={0, // top left
				tile_size /2, // top right
				tile_size * tile_size/2, // bottom left
				(tile_size + 1) * tile_size/2, // bottom right
				(tile_size + 1) * tile_size/4, // center
				tile_size/4, // top
				tile_size * tile_size/4,// left
				(tile_size + 2) * tile_size/4,   // right
				(2 * tile_size + 1) * tile_size/4};   // bottom

		int i,j,iq;
		int index,qindex;
		if (this.debugLevel>2) SDFA_INSTANCE.showArrays(bayer_mono_pixels, tile_size, tile_size, title+"-bayer");
		for (iq=0; iq<9;iq++) {
			index=quarterIndex[iq];
			qindex=0;
				for (i=0;i<hsize;i++) {
					for (j=0;j<hsize;j++) { //quarter_pixels[iq][qindex++]=input_pixels[index++];
						green0[qindex]=  bayer_mono_pixels[0][index];
						green3[qindex++]=bayer_mono_pixels[3][index++];
					}
				quarter_pixels[iq]=combineDiagonalGreens (green0, green3,  hsize, hsize);
					index+=hsize; // jump to the next line
				}
			quarter_pixels[iq]= normalizeAndWindow (quarter_pixels[iq], quarterHamming);
			if (this.debugLevel>2) SDFA_INSTANCE.showArrays(quarter_pixels[iq],hsize, hsize, title+"-new"+iq);
			quarter_patterns[iq] = findPattern(
					null, // 			DoubleFHT doubleFHT,
					quarter_pixels[iq],
					hsize,
					patternDetectParameters,
					min_half_period,
					max_half_period,
					greens,
					title+"Q_"+iq);

			if (quarter_patterns[iq]==null) return null;
		}
		if (this.debugLevel>2) {
			for (iq=0; iq<9;iq++) {
				System.out.println("Quarter="+iq+
						" W0x="+     IJ.d2s(quarter_patterns[iq][0][0],4)+
						" W0y="+     IJ.d2s(quarter_patterns[iq][0][1],4)+
						" W0_phase="+IJ.d2s(quarter_patterns[iq][0][2],2)+
						" W1x="+     IJ.d2s(quarter_patterns[iq][1][0],4)+
						" W1y="+     IJ.d2s(quarter_patterns[iq][1][1],4)+
						" W1_phase="+IJ.d2s(quarter_patterns[iq][1][2],2));
			}
		}
		/* Filter pattern coefficients to make sure they all match between quadrants (match to the center one)*/
		boolean patternsMatchedInitially=matchPatterns(quarter_patterns,quarter_patterns[4]);    // use center pattern
		if (this.debugLevel>2) {
			System.out.println(patternsMatchedInitially?"All quadrant wave vectors matched initially, no correction needed":"Some quadrant wave vectors were adjusted to match");
		}

		patternCorr=calcPatternNonLinear(quarter_patterns); // divide results by ,(FFT_SIZE/2)^2 - only first 5 patterns are used
		if (this.debugLevel>2) { /* increase LEVEL later */
			System.out.println("Pre- (1000x)   "+
					" Ax="+     IJ.d2s(1000*patternCorr[0]/(FFT_SIZE/2),5)+
					" Bx="+     IJ.d2s(1000*patternCorr[1]/(FFT_SIZE/2),5)+
					" Cx="+     IJ.d2s(1000*patternCorr[2]/(FFT_SIZE/2),5)+
					" Ay="+     IJ.d2s(1000*patternCorr[3]/(FFT_SIZE/2),5)+
					" By="+     IJ.d2s(1000*patternCorr[4]/(FFT_SIZE/2),5)+
					" Cy="+     IJ.d2s(1000*patternCorr[5]/(FFT_SIZE/2),5)+
					" Dx="+     IJ.d2s(1000*patternCorr[6],5)+
					" Ex="+     IJ.d2s(1000*patternCorr[7],5)+
					" Dy="+     IJ.d2s(1000*patternCorr[8],5)+
					" Ey="+     IJ.d2s(1000*patternCorr[9],5));
		}
		patternCorr=refinePatternNonLinear(quarter_patterns, // [tl,tr,bl,br, center][wv0, wv1][x,y,phase]
				patternCorr, //[ax,bx,cx,ay,by,cy]
				hsize ); // distance to quadrats center in sensor pixels ==FFT_SIZE/2



		//    for (i=0;i<patternCorr.length;i++)patternCorr[i]/= hsize;
		for (i=0;i<6;i++)patternCorr[i]/= hsize; /* Not linear Dx,Ex, Dy,Ey! */

		if (this.debugLevel>2) { /* increase LEVEL later */
			System.out.println("Corr (1000x)   "+
					" Ax="+     IJ.d2s(1000*patternCorr[0],5)+
					" Bx="+     IJ.d2s(1000*patternCorr[1],5)+
					" Cx="+     IJ.d2s(1000*patternCorr[2],5)+
					" Ay="+     IJ.d2s(1000*patternCorr[3],5)+
					" By="+     IJ.d2s(1000*patternCorr[4],5)+
					" Cy="+     IJ.d2s(1000*patternCorr[5],5)+
					" Dx="+     IJ.d2s(1000*patternCorr[6],5)+
					" Ex="+     IJ.d2s(1000*patternCorr[7],5)+
					" Dy="+     IJ.d2s(1000*patternCorr[8],5)+
					" Ey="+     IJ.d2s(1000*patternCorr[9],5));
		}

		double [][]result=new double [3][];
		result[0]=quarter_patterns[4][0].clone();
		result[1]=quarter_patterns[4][1].clone();
		result[2]=patternCorr.clone();
		return result;
	}
	/* ======================================================================== */
	public  double [] correlationContrast (
			double [] pixels,       // square pixel array
			double [] widowedGreens, // array to normalize correlation result
			double [][] wVectors,   // wave vectors (same units as the pixels array)
			double  contrastSelectSigmaCenter, // Gaussian sigma to select correlation centers (in PIXELS), 2.0
			double  contrastSelectSigmaOther,  // Gaussian sigma to select correlation off-centers centers (fraction of UV period), 0.1

			double x0,              // center coordinates
			double y0,
			String title){
		// for now - just comparison, later - switch to
		return correlationContrast (
				pixels,       // square pixel array
				wVectors,   // wave vectors (same units as the pixels array)
				contrastSelectSigmaCenter, // Gaussian sigma to select correlation centers (in PIXELS), 2.0
				contrastSelectSigmaOther,  // Gaussian sigma to select correlation off-centers centers (fraction of UV period), 0.1
				x0,              // center coordinates
				y0,
				title, // title base for optional plots names
				this.debugLevel);
	}
	public  double correlationContrastOld ( double [] pixels,       // square pixel array
			double [][] wVectors,   // wave vectors (same units as the pixels array)
			double ringWidth,       // ring (around r=0.5 dist to opposite corr) width
			double x0,              // center coordinates
			double y0,
			String title, // title base for optional plots names
			int debugLevel){
		int size=(int) Math.sqrt(pixels.length);
		double [] xy= new double [2];
		double [] uv;
		double r2,d;
		int i,j;
		/* opposite sign correlation points in uv are at uv=(0,-0.5),(0,0.5), (-0.5,0) and (0.5,0), with radius of (1/2)
    selecting center circle and a ring from 0.25 to 0.75 of the distance to opposite sign correlations */

		double r2WingsOuter= 0.0625*(1.0+ringWidth)*(1.0+ringWidth);
		double r2WingsInner= 0.0625*(1.0-ringWidth)*(1.0-ringWidth);
		double r2Center=0.0625*(ringWidth)*(ringWidth);
		if (debugLevel>2) System.out.println("rWingsOuter="+Math.sqrt(r2WingsOuter)+" rWingsInner="+Math.sqrt(r2WingsInner)+" rCenter="+Math.sqrt(r2Center));

		double valCenter=0.0;
		double valWings=0.0;
		double numCenter=0.0;
		double numWings=0.0;
		for (i=0;i<size;i++) {
			xy[1]=i-size/2-y0;
			for (j=0;j<size;j++) {
				xy[0]=j-size/2-x0;
				uv=matrix2x2_mul(wVectors,xy);
				r2=uv[0]*uv[0]+uv[1]*uv[1];
				if (r2<=r2WingsOuter) {
					d=pixels[i*size+j];
					if (r2<=r2Center){
						valCenter+=d*d;
						numCenter+=1.0;
					} else if (r2>r2WingsInner){
						valWings+=d*d;
						numWings+=1.0;
					}
				}
			}
		}
		if ((numWings==0.0) || (numCenter==0.0)) {
			if (debugLevel>1) System.out.println("Not enough data for correlation contrast: numCenter="+numCenter+" numWings="+numWings+
					" valCenter="+IJ.d2s(valCenter,2)+" valWings="+IJ.d2s(valWings,2));
			return -1.0;
		}
		double contrast=Math.sqrt((valCenter/numCenter)/(valWings/numWings));
		if (debugLevel>2) {
			System.out.println("Correlation contrast is "+contrast);
			double [] maskedPixels=new double[size*size];
			double [] u_value=new double[size*size];
			double [] v_value=new double[size*size];
			int index;
			for (i=0;i<size;i++) {
				xy[1]=i-size/2-y0;
				for (j=0;j<size;j++) {
					xy[0]=j-size/2-x0;
					uv=matrix2x2_mul(wVectors,xy);
					r2=uv[0]*uv[0]+uv[1]*uv[1];
					index=i*size+j;
					u_value[index]=uv[0];
					v_value[index]=uv[1];
					/*          r=Math.sqrt(r2);
          r-=Math.floor(r);
          floatPixels[index]=(float) r;*/
					if (((r2<=r2WingsOuter) && (r2>r2WingsInner)) || (r2<=r2Center)){
						maskedPixels[index]= pixels[index];
					} else {
						maskedPixels[index]=0.0;
					}
				}
			}
			double [][] dbgPixels={pixels,maskedPixels,u_value,v_value};
			String [] titles={"all","masked","u","v"};
			(new ShowDoubleFloatArrays()).showArrays(
					dbgPixels,
					size,
					size,
					true,
					title+"_CORR_MASK",
					titles);
		}
		return contrast;
	}

	public double [] correlationContrast (
			double [] pixels,       // square pixel array
			double [][] wVectors,   // wave vectors (same units as the pixels array)
			double  sigma_center, // GGaussian sigma to select correlation centers (in PIXELS), 1.5
			double  sigma_other,  // Gaussian sigma to select correlation off-centers centers (fraction of UV period), 0.1
			double x0,              // center coordinates
			double y0,
			String title, // title base for optional plots names
			int debugLevel){
		double avg_rad = 2.0/Math.sqrt(wVectors[0][0]*wVectors[0][0]+wVectors[0][1]*wVectors[0][1]+wVectors[1][0]*wVectors[1][0]+wVectors[1][1]*wVectors[1][1]);
		double max_center_sigma = 0.05 * avg_rad; // approximate - make configurable or just relative rather than absolute?
		if (sigma_center > max_center_sigma) {
			sigma_center = max_center_sigma;
		}
		double [] badContrasts={-1.0,-1.0};
		double sigma32_center=9*sigma_center*sigma_center; // in pixels
		double k_center=-0.5/(sigma_center*sigma_center);  // in pixels
		double sigma32_other=9*sigma_other*sigma_other;    // in periods
		double k_other=-0.5/(sigma_other*sigma_other);     // in periods
		double [][] sampleCentersXY={{0.0,0.0},{0.25,0.25},{0.25,-0.25},{-0.25,0.25},{-0.25,-0.25}};
		// System.out.println("avg_rad="+avg_rad);
		int [] sampleTypes = {0,1,1,1,1};
		int size=(int) Math.sqrt(pixels.length);
		double [] xy= new double [2];
		double [] uv;
		double r2;
		int i,j;
		// TODO: limit sigma_center to fit between others;
		double [] dbgMask= new double[size*size];
		for (int n=0;n<dbgMask.length;n++) dbgMask[n]=0.0;
		double [] s={0.0,0.0};
		double [] w={0.0,0.0};
		//		double [] dbg_weights = new double [size*size];
		for (i=0;i<size;i++) {
			xy[1]=i-size/2-y0;
			for (j=0;j<size;j++) {
				int index=i*size+j;
				xy[0]=j-size/2-x0;
				uv=matrix2x2_mul(wVectors,xy);
				for (int np=0;np<sampleCentersXY.length;np++){
					if (sampleTypes[np] == 0) { // center spot, size in pixels
						r2 = xy[0]*xy[0] + xy[1]*xy[1];
						if (r2 < sigma32_center){
							double m=Math.exp(k_center*r2);
							dbgMask[index]+=m;
							w[sampleTypes[np]]+=m;
							double d=m*pixels[index];
							if (sampleTypes[np]>0)  d *= pixels[index]; // squared
							s[sampleTypes[np]]+=d;
						}
					} else {  // between correlation spots, size relative to the  periods
						double dx=uv[0]-sampleCentersXY[np][0];
						double dy=uv[1]-sampleCentersXY[np][1];
						r2=dx*dx+dy*dy;
						if (r2 < sigma32_other){
							double m=Math.exp(k_other*r2);
							dbgMask[index]+=m;
							w[sampleTypes[np]]+=m;
							double d=m*pixels[index];
							if (sampleTypes[np]>0)  d*=pixels[index]; // squared
							s[sampleTypes[np]]+=d;
						}

					}
				}
			}
		}
		if ((w[0]==0.0) || (w[1]==0.0)) {
			if (debugLevel>1) System.out.println("Not enough data for correlation contrast: center - w[0]="+w[0]+" opposite - w[1]="+w[1]);
			return badContrasts;
		}
		double [][] dbg_corr_mask = {pixels, dbgMask};
		double aCenter= s[0]/w[0];
		double aQuiet=Math.sqrt(s[1]/w[1]);
		double rContrast=aCenter/aQuiet;
		double aContrast=aCenter/size/size;

		double [] contrasts={rContrast,aContrast};
		if (debugLevel>2){
			System.out.println("correlationContrast() rContrast="+rContrast+" aContrast="+ aContrast+" aCenter="+aCenter+" aQuiet="+aQuiet+" w[0]="+w[0]+" w[1]="+w[1]+" s[0]="+s[0]+" s[1]="+s[1]);
		}
		if (debugLevel>2) {
			System.out.println("Correlation contrast is: relative="+rContrast+" absolute="+aContrast);
			double [][] dbgPixels={pixels,dbgMask};
			String [] titles={"all","mask"};
			(new ShowDoubleFloatArrays()).showArrays(
					dbgPixels,
					size,
					size,
					true,
					title+"_MASK",
					titles);
		}
		return contrasts;
	}
	public  double correlationContrastOld2 (
			double [] pixels,       // square pixel array
			double [] widowedGreens, // array to normailze correlation result
			double [][] wVectors,   // wave vectors (same units as the pixels array)
			double sigma,
			double sigmaNorm,       // to measure variations for normalization of the contrast
			double x0,              // center coordinates
			double y0,
			String title, // title base for optional plots names
			int debugLevel){
		// TODO: make configurable parameters
		//		double sigma=0.1;
		//		double sigmaNorm=0.5; // to measure variations for normalization of the contrast

		double sigma32=9*sigma*sigma;
		double k=-0.5/(sigma*sigma);

		double sigmaNorm32=9*sigmaNorm*sigmaNorm;
		double kNorm=-0.5/(sigmaNorm*sigmaNorm);

		double [][] sampleCentersXY={{0.0,0.0},{0.0,0.5},{0.5,0.0},{0.0,-0.5},{-0.5,0.0}};
		int [] sampleTypes = {0,1,1,1,1};
		int size=(int) Math.sqrt(pixels.length);
		double [] xy= new double [2];
		double [] uv;
		double r2;
		int i,j;

		/* opposite sign correlation points in uv are at uv=(0,-0.5),(0,0.5), (-0.5,0) and (0.5,0), with radius of (1/2)
    selecting center circle and a ring from 0.25 to 0.75 of the distance to opposite sign correlations */
		double [] dbgMask= new double[size*size];
		for (int n=0;n<dbgMask.length;n++) dbgMask[n]=0.0;
		double [] s={0.0,0.0};
		double [] w={0.0,0.0};
		double S0=0.0,S1=0.0,S2=0.0;
		double SG1=0.0,SG2=0.0;
		// Find measured pixels variations in the window
		for (i=0;i<size;i++) {
			xy[1]=i-size/2-y0;
			for (j=0;j<size;j++) {
				int index=i*size+j;
				xy[0]=j-size/2-x0;
				uv=matrix2x2_mul(wVectors,xy);
				for (int np=0;np<sampleCentersXY.length;np++){
					double dx=uv[0]-sampleCentersXY[np][0];
					double dy=uv[1]-sampleCentersXY[np][1];
					r2=dx*dx+dy*dy;
					if (r2<sigma32){
						double m=Math.exp(k*r2);
						dbgMask[index]+=m;
						w[sampleTypes[np]]+=m;
						s[sampleTypes[np]]+=m*pixels[index];
					}
				}
				r2=uv[0]*uv[0]+uv[1]*uv[1];
				if (r2<sigmaNorm32){
					double m=Math.exp(kNorm*r2);
					S0+=m;
					S1+=m*pixels[index];
					S2+=m*pixels[index]*pixels[index];
					SG1=m*widowedGreens[index];
					SG2=m*widowedGreens[index]*widowedGreens[index];
				}
			}
		}
		if ((w[0]==0.0) || (w[1]==0.0)) {
			if (debugLevel>1) System.out.println("Not enough data for correlation contrast: center - w[0]="+w[0]+" opposite - w[1]="+w[1]);
			return -1.0;
		}
		double ref=Math.sqrt(S2*S0-S1*S1)/S0;
		double refG=Math.sqrt(SG2*S0-SG1*SG1)/S0;
		double contrast=((s[0]/w[0]) -(s[1]/w[1]))/ref;
		double contrastG=((s[0]/w[0]) -(s[1]/w[1]))/refG; ///size;
		if (debugLevel>2){
			System.out.println("correlationContrast() corr_diff="+(((s[0]/w[0]) -(s[1]/w[1])))+" contrast="+contrast+" w[0]="+w[0]+" w[1]="+w[1]+" s[0]="+s[0]+" s[1]="+s[1]);
			System.out.println("correlationContrast() S0="+S0+" S1="+S1+" S2="+S2+" ref="+ref);
			System.out.println("correlationContrast() contrastG="+contrastG+" S0="+S0+" SG1="+SG1+" SG2="+SG2+" refG="+refG);
		}
		//		if (contrast>3.0){
		//			System.out.println("correlationContrast() contrast="+contrast+" w[0]="+w[0]+" w[1]="+w[1]+" s[0]="+s[0]+" s[1]="+s[1]+" S0="+S0+" S1="+S1+" S2="+S2+" ref="+ref);
		//		}
		//		double contrast=Math.sqrt((s[0]/w[0]) /(s[1]/w[1]));
		//		double contrast=((s[0]/w[0]) -(s[1]/w[1]))/(size*size);
		if (debugLevel>2) {
			System.out.println("Correlation contrast is "+contrast);
			double [][] dbgPixels={pixels,dbgMask};
			String [] titles={"all","mask"};
			(new ShowDoubleFloatArrays()).showArrays(
					dbgPixels,
					size,
					size,
					true,
					title+"_CORR_MASK",
					titles);
		}
		return contrast;
	}



	/* ======================================================================== */
	public  double[] correlateWithModel (double [] imagePixels,  // measured pixel array
			double [] modelPixels,  // simulated (model) pixel array)
			double sigma,   // Sigma for high pass filtering TODO: implement!
			String title) { // title base for optional plots names

		if (imagePixels.length!=modelPixels.length) {
			IJ.showMessage("Error","Arrays have different sizes - imagePixels.length="+imagePixels.length+", modelPixels.length="+ modelPixels.length);
			return null;
		}
		int size = (int) Math.sqrt(imagePixels.length);
		ImageProcessor ip,ip_model;
		FHT fht,fht_model;
		double[][][] fft_complex,fft_model;
		int i,j;
		double a;

		float [] floatImagePixels=new float[size*size];
		/* convert to float for image processor; */
		for (i=0;i<(size*size); i++) floatImagePixels[i]=(float) imagePixels[i];
		ip = new FloatProcessor(size,size);
		ip.setPixels(floatImagePixels);
		fht =  new FHT(ip);
		// Swapping quadrants, so the center will be 0,0
		fht.swapQuadrants();
		// get to frequency domain
		fht.transform();
		floatImagePixels=(float []) fht.getPixels();

		if((this.debugLevel>5) && (title!="")) {
			ImageProcessor ip_fht = new FloatProcessor(size,size);
			ip_fht.setPixels(floatImagePixels);
			ip_fht.resetMinAndMax();
			ImagePlus imp_fht= new ImagePlus(title+"_FHT_image", ip_fht);
			imp_fht.show();
		}
		// Convert from FHT to complex FFT
		fft_complex= FHT2FFTHalf (fht,size);
		float [] floatModelPixels=new float[size*size];
		// convert to float for image processor;
		for (i=0;i<(size*size); i++) floatModelPixels[i]=(float) modelPixels[i];
		ip_model = new FloatProcessor(size,size);
		ip_model.setPixels(floatModelPixels);
		fht_model =  new FHT(ip_model);
		// Swapping quadrants, so the center will be 0,0
		fht_model.swapQuadrants();
		// get to frequency domain
		fht_model.transform();
		floatModelPixels=(float []) fht_model.getPixels();
		if ((this.debugLevel>5) && (title!="")) {
			ImageProcessor ip_fht_model = new FloatProcessor(size,size);
			ip_fht_model.setPixels(floatModelPixels);
			ip_fht_model.resetMinAndMax();
			ImagePlus imp_fht_model= new ImagePlus(title+"_FHT_model", ip_fht_model);
			imp_fht_model.show();
		}
		// Convert from FHT to complex FFT
		fft_model= FHT2FFTHalf (fht_model,size);

		// multiply fft_complex by fft_nominator
		for (i=0;i<fft_complex.length; i++) for (j=0;j<fft_complex[0].length;j++) {
			a=                    fft_complex[i][j][0]*fft_model[i][j][0]+fft_complex[i][j][1]*fft_model[i][j][1]; // already changed Im() sign
			fft_complex[i][j][1]=-fft_complex[i][j][0]*fft_model[i][j][1]+fft_complex[i][j][1]*fft_model[i][j][0]; // already changed Im() sign
			fft_complex[i][j][0]=a;
		}
		/* Add sigma high=pass filtering here */
		// Convert fft array back to fht array and
		// set fht_target pixels with new values
		fht.setPixels (floatFFTHalf2FHT (fft_complex,size));
		/// optionally show the result
		if ((this.debugLevel>5) && (title!="")) {
			ImageProcessor ip_fht2 = new FloatProcessor(size,size);
			ip_fht2.setPixels(floatFFTHalf2FHT (fft_complex,size));
			ip_fht2.resetMinAndMax();
			ImagePlus imp_fht2= new ImagePlus(title+"-corr-sigma"+sigma, ip_fht2);
			imp_fht2.show();
		}

		/// transform to space

		fht.inverseTransform();
		fht.swapQuadrants();
		fht.resetMinAndMax();

		//   ImagePlus imp= new ImagePlus(title, ip_fht);
		if ((this.debugLevel>2) && (title!="")) {
			ImagePlus imp_corr= new ImagePlus(title+"_Correlated_filt-"+sigma, fht);
			imp_corr.show();
		}
		//   return direct_target;
		floatImagePixels =(float[])fht.getPixels();
		double [] pixels=new double[floatImagePixels.length];
		for (i=0;i<floatImagePixels.length;i++) pixels[i]=floatImagePixels[i];
		return pixels;
	}


	/**
	Refining non-linear mesh matching by comparing phases in the centers of 4 quadrants and the very center.
	Can only compensate to a fraction of mesh period (TBD - total range, probably +/-half period),
	so non-linear coefficients should be already known to that precision
	9 measurements are used here - top-left, top-right,bottom-left, bottom-right, center, top, left, right,bottom

	 */
	private	double [] refinePatternNonLinear(double [][][] qp, // [tl,tr,bl,br, center][wv0, wv1][x,y,phase]
			double [] nonlin, //[ax,bx,cx,ay,by,cy]
			int size ) { // distance to quadrants center in sensor pixels ==FFT_SIZE/2

		int iq,i,j;
		double [][] xy=new double [qp.length][2];
		double [] uv= new double [2];
		double [] duv=new double [2];

		//  double [][][] wl=new double [5][2][2]; // Wave length vectors - same direction as wavevectors, length=distance between wavefronts
		double [][][] wp=new double [9][2][3]; // pattern vectors (with phase)
		double x1,y1;
		double xq=0;
		double yq=0;
		//probably only wl[4] is needed
		for (iq=0;iq<9;iq++) for (i=0;i<2;i++) for (j=0;j<2;j++) wp[iq]= waveVectorsToPatternVectors(qp[iq][0], qp[iq][1]);

		if (this.debugLevel>2) { /* increase LEVEL later */
			for (iq=0;iq<5;iq++){
				System.out.println(" wp["+iq+"][0][0]="+     IJ.d2s(wp[iq][0][0],4)+
						" wp["+iq+"][0][1]="+     IJ.d2s(wp[iq][0][1],4)+
						" wp["+iq+"][0][2]="+     IJ.d2s(wp[iq][0][2],4)+"("+ IJ.d2s(wp[iq][0][2]/2/Math.PI,4)+")"+
						" wp["+iq+"][1][0]="+     IJ.d2s(wp[iq][1][0],4)+
						" wp["+iq+"][1][1]="+     IJ.d2s(wp[iq][1][1],4)+
						" wp["+iq+"][1][2]="+     IJ.d2s(wp[iq][1][2],4)+"("+ IJ.d2s(wp[iq][1][2]/2/Math.PI,4)+")");
			}
		}
		for (iq=0;iq<9;iq++) {
			if (iq==4) continue; // nothing to calculate in the center
			switch (iq){
			case 0: xq=-1.0;yq=-1.0; break;
			case 1: xq= 1.0;yq=-1.0; break;
			case 2: xq=-1.0;yq= 1.0; break;
			case 3: xq= 1.0;yq= 1.0; break;
			case 4: xq= 0.0;yq= 0.0; break;
			case 5: xq= 0.0;yq=-1.0; break;
			case 6: xq=-1.0;yq= 0.0; break;
			case 7: xq= 1.0;yq= 0.0; break;
			case 8: xq= 0.0;yq= 1.0; break;
			}

			x1=size*(xq + nonlin[0]*xq*xq+ nonlin[1]*yq*yq+ 2* nonlin[2]*xq*yq +nonlin[6]*xq+ nonlin[7]*yq); // in pixels
			y1=size*(yq + nonlin[3]*xq*xq+ nonlin[4]*yq*yq+ 2* nonlin[5]*xq*yq +nonlin[8]*xq+ nonlin[9]*yq); // in pixels

			/* convert x1,y1 into wp vector coordiantes */
			uv[0]=(wp[4][1][1]*x1-wp[4][1][0]*y1)/(wp[4][0][0]*wp[4][1][1]-wp[4][0][1]*wp[4][1][0]); // wl in center vectors, not local !
			uv[1]=(wp[4][0][0]*y1-wp[4][0][1]*x1)/(wp[4][0][0]*wp[4][1][1]-wp[4][0][1]*wp[4][1][0]);
			/* Actually phases seem to be the same? */
			duv[0]=uv[0]-Math.round(uv[0])- (wp[iq][0][2]-wp[4][0][2])/(Math.PI*2);
			duv[1]=uv[1]-Math.round(uv[1])- (wp[iq][1][2]-wp[4][1][2])/(Math.PI*2);

			if (this.debugLevel>2) { /* increase LEVEL later */
				System.out.println("iq=  "+ iq+
						" x1="+     IJ.d2s(x1,4)+
						" y1="+     IJ.d2s(y1,4)+
						" uv[0]"+   IJ.d2s(uv[0],4)+
						" uv[1]"+   IJ.d2s(uv[1],4)+
						" duv[0]"+   IJ.d2s(duv[0],4)+
						" duv[1]"+   IJ.d2s(duv[1],4));
			}
			/* re-normalize phases to be +/- 0.5 (+/-PI) range */
			duv[0]=duv[0]-Math.round(duv[0]);
			duv[1]=duv[1]-Math.round(duv[1]);
			if (this.debugLevel>2) { /* increase LEVEL later */
				System.out.println("iq=  "+ iq+" ------- "+
						" duv[0]"+   IJ.d2s(duv[0],4)+	/* ======================================================================== */

						" duv[1]"+   IJ.d2s(duv[1],4));
			}
			/* Fix half period vertical/half period horizontal shift - is that needed?*/
			if (Math.abs(duv[0])>0.25) {
				duv[0]+=0.5;
				duv[1]+=0.5;
				duv[0]=duv[0]-Math.round(duv[0]);
				duv[1]=duv[1]-Math.round(duv[1]);
				if (this.debugLevel>2) {
					System.out.println("Correct phase shift >0.25 in quadrant "+ iq+ ", now"+
							" duv[0]"+   IJ.d2s(duv[0],4)+
							" duv[1]"+   IJ.d2s(duv[1],4));
				}
			}

			/* Verify here that phase adjustment is within range, fail otherwise */
			if ((Math.abs(duv[0])>0.5) || (Math.abs(duv[1])>0.5)) {
				if (this.debugLevel>0) {
					System.out.println("Error: in quadrant "+ iq+" - attempted to adjust phase too much (>+/- pi/2), keeping initial parematers");
				}
				return nonlin;
			}

			/* convert duv to x,y */
			xy[iq][0]=x1-wp[4][0][0]*duv[0]-wp[4][1][0]*duv[1];
			xy[iq][1]=y1-wp[4][0][1]*duv[0]-wp[4][1][1]*duv[1];


			if (this.debugLevel>2) { /* increase LEVEL later */
				System.out.println(" xy["+iq+"][0]="+   IJ.d2s(xy[iq][0],4)+
						" xy["+iq+"][1]="+   IJ.d2s(xy[iq][1],4));
			}
			/* convert xy to non-linear differences, remove pixels dimensions - quadrats +/- 1.0 */
			xy[iq][0]=(xy[iq][0]-size*xq)/size;
			xy[iq][1]=(xy[iq][1]-size*yq)/size;

			if (this.debugLevel>2) { /* increase LEVEL later */
				System.out.println("Quadrant "+iq+": Original non-linear difference was"+
						" x="+   IJ.d2s(nonlin[0]*xq*xq+ nonlin[1]*yq*yq+ nonlin[2]*xq*yq,4)+
						" y="+   IJ.d2s(nonlin[3]*xq*xq+ nonlin[4]*yq*yq+ nonlin[5]*xq*yq,4));
				System.out.println("Quadrant "+iq+": Refined  non-linear difference is"+
						" x="+   IJ.d2s(xy[iq][0],4)+
						" y="+   IJ.d2s(xy[iq][1],4));
			}

		}

		/* Do the refinement itself - recalculate coefficients minimizing errors in 5 points */
		/**

	    x1=size*(xq + nonlin[0]*xq*xq+ nonlin[1]*yq*yq+ 2* nonlin[2]*xq*yq +nonlin[6]*xq+ nonlin[7]*yq); // in pixels
	    y1=size*(yq + nonlin[3]*xq*xq+ nonlin[4]*yq*yq+ 2* nonlin[5]*xq*yq +nonlin[8]*xq+ nonlin[9]*yq); // in pixels
	  dx=AX*x^2 +BX*y^2 +2*CX*x*y + DX*x +EX *y
	  dy=AY*x^2 +BY*y^2 +2*CY*x*y + DY*x +EY *y

	  dx0= AX +BX +2*CX -DX -EX
	  dy0= AY +BY +2*CY -DY -EY

	  dx1= AX +BX -2*CX +DX -EX
	  dy1= AY +BY -2*CY +DY -EY

	  dx2= AX +BX -2*CX -DX +EX
	  dy2= AY +BY -2*CY -DY +EY

	  dx3= AX +BX +2*CX +DX +EX
	  dy3= AY +BY +2*CY +DY +EY

	  dx5=    +BX           -EX
	  dy5=    +BY           -EY

	  dx6= AX           -DX
	  dy6= AY           -DY

	  dx7= AX           +DX
	  dy7= AY           +DY

	  dx8=    +BX           +EX
	  dy8=    +BY           +EY

	//---------
	-(1)   dx0= AX +BX +2*CX -DX -EX
	-(2)   dx1= AX +BX -2*CX +DX -EX
	+(3)   dx2= AX +BX -2*CX -DX +EX
	+(4)   dx3= AX +BX +2*CX +DX +EX
	-(5)   dx5=    +BX           -EX
	(6)   dx6= AX           -DX
	(7)   dx7= AX           +DX
	+(8)   dx8=    +BX           +EX
	-(1)+(2)-(3)+(4)-(6)+(7) DX=()-dx0+dx1-dx2+dx3-dx6+dx7)/6
	-(1)-(2)+(3)+(4)-(5)+(8) EX=()-dx0-dx1+dx2+dx3-dx5+dx8)/6
	-(1)-(2)+(3)+(4)-(5)+(8) EY=()-dy0-dy1+dy2+dy3-dy5+dy8)/6
	-(1)+(2)-(3)+(4)-(6)+(7) DY=()-dy0+dy1-dy2+dy3-dy6+dy7)/6

	  AX+BX+2*CX = p1x= (dx0+dx3)/2
	  AY+BY+2*CY = p1y= (dy0+dy3)/2
	  AX+BX-2*CX = p2x= (dx1+dx2)/2
	  AY+BY-2*CY = p2y= (dy1+dy2)/2
	     BX      = p3x= (dx5+dx8)/2
	     BY      = p3y= (dy5+dy8)/2
	  AX         = p4x= (dx6+dx7)/2
	  AY         = p4y= (dy6+dy7)/2
	//minimizing sum of squares of errors:
	  CX= (p1x-p2x)/4
	  CY= (p1y-p2y)/4
	  AX= (3*p4x-2*p3x+p1x+p2x)/5
	  AY= (3*p4y-2*p3y+p1y+p2y)/5
	  BX= (3*p3x-2*p4x+p1x+p2x)/5
	  BY= (3*p3y-2*p4y+p1y+p2y)/5
		 */
		double p1x =(xy[0][0]+xy[3][0])/2;
		double p1y =(xy[0][1]+xy[3][1])/2;

		double p2x =(xy[1][0]+xy[2][0])/2;
		double p2y =(xy[1][1]+xy[2][1])/2;

		double p3x =(xy[5][0]+xy[8][0])/2;
		double p3y =(xy[5][1]+xy[8][1])/2;

		double p4x =(xy[6][0]+xy[7][0])/2;
		double p4y =(xy[6][1]+xy[7][1])/2;

		double [] rslt=new double[10];
		rslt[2]=(p1x-p2x)/4; //CX
		rslt[5]=(p1y-p2y)/4; //CY
		rslt[0]= (3*p4x-2*p3x+p1x+p2x)/5 ; //AX
		rslt[3]= (3*p4y-2*p3y+p1y+p2y)/5; // AY
		rslt[1]= (3*p3x-2*p4x+p1x+p2x)/5; // BX
		rslt[4]= (3*p3y-2*p4y+p1y+p2y)/5; // BY
		/*
	-(1)+(2)-(3)+(4)-(6)+(7) DX=(-dx0+dx1-dx2+dx3-dx6+dx7)/6
	-(1)-(2)+(3)+(4)-(5)+(8) EX=(-dx0-dx1+dx2+dx3-dx5+dx8)/6
	-(1)-(2)+(3)+(4)-(5)+(8) EY=(-dy0-dy1+dy2+dy3-dy5+dy8)/6
	-(1)+(2)-(3)+(4)-(6)+(7) DY=(-dy0+dy1-dy2+dy3-dy6+dy7)/6
		 */
		rslt[6]=(-xy[0][0]+xy[1][0]-xy[2][0]+xy[3][0]-xy[6][0]+xy[7][0])/6; // DX=(-dx0+dx1-dx2+dx3-dx6+dx7)/6
		rslt[7]=(-xy[0][0]-xy[1][0]+xy[2][0]+xy[3][0]-xy[5][0]+xy[8][0])/6; // EX=(-dx0-dx1+dx2+dx3-dx5+dx8)/6
		rslt[8]=(-xy[0][1]+xy[1][1]-xy[2][1]+xy[3][1]-xy[6][1]+xy[7][1])/6; // EY=(-dy0-dy1+dy2+dy3-dy5+dy8)/6
		rslt[9]=(-xy[0][1]-xy[1][1]+xy[2][1]+xy[3][1]-xy[5][1]+xy[8][1])/6; // DY=(-dy0+dy1-dy2+dy3-dy6+dy7)/6

		return rslt; // temporary
	}
	/* ======================================================================== */
	public double[][] findPatternUsedNow(
			double [] input_pixels, // pixel array to process
			int size, // FFT size
			PatternDetectParameters patternDetectParameters,
			boolean greens, // this is a pattern for combined greens (diagonal), adjust results accordingly
			String title) { // title prefix to use for debug  images
		double [] pixels=input_pixels.clone();
		double [][]result=new double [2][3];
		//System.out.println("pixels.length="+pixels.length); //4096

		ImageProcessor ip, ip1;
		FHT fht, fht1;
		double[][][] fft_complex,fft_corr;
		double[][] fft_gamma;
		int i,j;
		double DCLevel=0.0;
		double a;
		float []floatPixels=new float[pixels.length];
		for (i=0;i<pixels.length; i++)  DCLevel+=pixels[i];
		DCLevel/=(size*size);
		for (i=0;i<pixels.length; i++)  pixels[i]-=DCLevel;
		// convert to float for image processor;
		for (i=0;i<pixels.length; i++) floatPixels[i]=(float) pixels[i];
		ip = new FloatProcessor(size,size);
		ip.setPixels(floatPixels);
		if (this.debugLevel>8) {
			ip.resetMinAndMax();
			ImagePlus imp_direct=  new ImagePlus(title+"_Direct_"+patternDetectParameters.corrGamma, ip);
			imp_direct.show();
		}
		fht =  new FHT(ip);
		// Swapping quadrants, so the center will be 0,0
		fht.swapQuadrants();
		// get to frequency domain
		fht.transform();
		if (this.debugLevel>5) {
			floatPixels=(float []) fht.getPixels();
			ImageProcessor ip_fht = new FloatProcessor(size,size);
			ip_fht.setPixels(floatPixels);
			ip_fht.resetMinAndMax();
			ImagePlus imp_fht= new ImagePlus(title+"_FHT", ip_fht);
			imp_fht.show();
		}

		// Convert from FHT to complex FFT
		fft_complex= FHT2FFTHalf (fht,size);
		// will need fft_complex  again later for later phase pattern measurements, calculate fft_gamma for correlation (pattern 2 frequencies measurement)
		fft_gamma=new double [size][size];
		floatPixels=new float[pixels.length];
		DCLevel=0.0;
		for (i=0;i<fft_complex.length; i++) for (j=0;j<fft_complex[0].length;j++) {
			fft_gamma[i][j]=Math.pow(fft_complex[i][j][0]*fft_complex[i][j][0]+fft_complex[i][j][1]*fft_complex[i][j][1],patternDetectParameters.corrGamma);
			DCLevel+=fft_gamma[i][j];
			floatPixels[i*size+j]=(float) fft_gamma[i][j];
		}
		DCLevel/=(fft_complex.length*fft_complex[0].length);
		for (i=0;i<fft_complex.length; i++) for (j=0;j<fft_complex[0].length;j++) {
			floatPixels[i*size+j]-=DCLevel;
			if ((i>0)&& (i<(size/2))){
				floatPixels[(size-i)*size+((size-j)%size)]=floatPixels[i*size+j];
			}
		}

		/* TODO:  maybe it is better to find the pattern frequencies just here, without converting back.
    After rejecting low frequencies, there seem to be just 2 nice maximums - easy to extract*/
		// now perform direct FFT of gamma(power spectrum)
		ip1 = new FloatProcessor(size,size);
		ip1.setPixels(floatPixels);
		if (this.debugLevel>7) {
			ip1.resetMinAndMax();
			ImagePlus imp1=  new ImagePlus(title+"_gamma(ps)_"+patternDetectParameters.corrGamma, ip1);
			imp1.show();
		}
		fht1 =  new FHT(ip1);
		// Swapping quadrants, so the center will be 0,0
		fht1.swapQuadrants();
		fht1.transform();
		fft_corr= FHT2FFTHalf (fht1,size);
		double[] highPassFilter=new double[fft_complex[0].length];
		double expK=(patternDetectParameters.corrSigma>0)?(1.0/(2*patternDetectParameters.corrSigma*patternDetectParameters.corrSigma)):0.0;
		for (j=0;j<=fft_complex[0].length/2;j++) {
			highPassFilter[j]=(expK>0.0)?(1.0-Math.exp(-(expK*j*j))):1.0;
			if (j>0) highPassFilter[highPassFilter.length-j]=highPassFilter[j];
		}

		for (i=0;i<fft_complex.length; i++) for (j=0;j<fft_complex[0].length;j++) {
			fft_corr[i][j][0]=highPassFilter[i]*highPassFilter[j]*(fft_corr[i][j][0]*fft_corr[i][j][0]+fft_corr[i][j][1]*fft_corr[i][j][1]);
			fft_corr[i][j][1]=0.0;
		}

		// Convert fft array back to fht array and
		// set fht_target pixels with new values
		fht1.setPixels (floatFFTHalf2FHT (fft_corr,size));   /* FIXME: - done, there is no difference as Im()==0 */
		/// optionally show the result
		if (this.debugLevel>7) {
			ImageProcessor ip_fht2 = new FloatProcessor(size,size);
			ip_fht2.setPixels(floatFFTHalf2FHT (fft_corr,size));
			ip_fht2.resetMinAndMax();
			ImagePlus imp_fht2= new ImagePlus(title+"_fht_corr_"+patternDetectParameters.corrGamma, ip_fht2);
			imp_fht2.show();
		}

		/// transform to space
		fht1.inverseTransform();
		floatPixels=(float []) fht1.getPixels();
		a=1/floatPixels[0];
		for (i=0; i<floatPixels.length; i++){
			floatPixels[i]*=a;
		}
		fht1.setPixels(floatPixels);

		//System.out.println("2:y="+y+" x="+x+" base_b="+base_b+" base="+base);

		fht1.swapQuadrants();

		if (this.debugLevel>2) {
			fht1.resetMinAndMax();
			ImagePlus imp_corr= new ImagePlus(title+"_corr_"+patternDetectParameters.corrGamma, fht1);
			imp_corr.show();
		}
		//   return direct_target;
		floatPixels =(float[])fht1.getPixels();
		for (i=0;i<floatPixels.length;i++) pixels[i]=floatPixels[i];


		int [][] max2OnSpectrum=  findFirst2MaxOnSpectrum (fft_complex, // complex, top half, starting from 0,0
				1,   // skip +- from (0,0) and previous max - add parameter to dialog?
				0.5); // 0.5 - 30deg. orthogonality of 2 vectors - 1.0 - perpendicular, 0.0 - parallel - add parameter to dialog?
		/**TODO:  get out on failure */
		if (max2OnSpectrum==null) {
			if (this.debugLevel>2){
				System.out.println("findPattern() 1: Failed to find a pattern");
				if (this.debugLevel>2){
					SDFA_INSTANCE.showArrays(input_pixels, "failed-findPattern-1-");
				}
			}
			return null;
		}
		/* Trying to filter out unreasonable maximums (if there is no pattern at all) */
		double maxFrequency=0.25*fft_complex.length;
		if ((Math.abs(max2OnSpectrum[0][0])>maxFrequency) ||
				(Math.abs(max2OnSpectrum[0][1])>maxFrequency) ||
				(Math.abs(max2OnSpectrum[1][0])>maxFrequency) ||
				(Math.abs(max2OnSpectrum[1][1])>maxFrequency)) {
			if (this.debugLevel>2) {
				System.out.println("Failed to detect pattern, as frequecy is above limit="+IJ.d2s(maxFrequency,2));
				System.out.println("Maximum 1 on spectrum:  x="+IJ.d2s(max2OnSpectrum[0][0],4)+" y="+IJ.d2s(max2OnSpectrum[0][1],4));
				System.out.println("Maximum 2 on spectrum:  x="+IJ.d2s(max2OnSpectrum[1][0],4)+" y="+IJ.d2s(max2OnSpectrum[1][1],4));
			}
			return null;
		}
		if (this.debugLevel>6) {
			System.out.println("Maximum 1 on spectrum:  x="+IJ.d2s(max2OnSpectrum[0][0],4)+" y="+IJ.d2s(max2OnSpectrum[0][1],4));
			System.out.println("Maximum 2 on spectrum:  x="+IJ.d2s(max2OnSpectrum[1][0],4)+" y="+IJ.d2s(max2OnSpectrum[1][1],4));
		}

		int [][] startPoints={{max2OnSpectrum[0][0]+max2OnSpectrum[1][0], max2OnSpectrum[0][1]+max2OnSpectrum[1][1]},
				{max2OnSpectrum[0][0]-max2OnSpectrum[1][0], max2OnSpectrum[0][1]-max2OnSpectrum[1][1]}};
		if (startPoints[1][1] <0) { /* startPoints[1][1] > 0 anyway */
			startPoints[1][0]= -startPoints[1][0];
			startPoints[1][1]= -startPoints[1][1];
		}
		if (this.debugLevel>2) {
			System.out.println("Predicted correlation maximum 1 from spectrum:  x="+IJ.d2s(startPoints[0][0],4)+" y="+IJ.d2s(startPoints[0][1],4));
			System.out.println("Predicted correlation maximum 2 from spectrum:  x="+IJ.d2s(startPoints[1][0],4)+" y="+IJ.d2s(startPoints[1][1],4));
		}

		double[][] max2=  findFirst2MaxOnCorrelation(
				pixels,
				startPoints,
				patternDetectParameters
				);

		/**TODO:  get out on failure */
		if (max2==null) {
			if (this.debugLevel>2){
				System.out.println("findPattern() 2: Failed to find a pattern");
				if (this.debugLevel>2){
					SDFA_INSTANCE.showArrays(input_pixels, "failed-findPattern-2-");
				}

			}
			return null;
		}
		/* these are combined greens, convert vectors to original pixel space) */
		if (greens) {
			double [][] rotMatrix= {{1.0,-1.0},{1.0,1.0}};
			double [][] max2orig= matrix2x2_mul(max2,rotMatrix);
			for (i=0;i<2;i++) for (j=0;j<2;j++) result[i][j]=max2orig[i][j]; // result is [2][3], max2orig is [2][2]
			if (this.debugLevel>2) {
				System.out.println("Corrected to original pixels[0]  x="+IJ.d2s(result[0][0],4)+" y="+IJ.d2s(result[0][1],4));
				System.out.println("Corrected to original pixels[1]  x="+IJ.d2s(result[1][0],4)+" y="+IJ.d2s(result[1][1],4));
			}
		} else {
			for (i=0;i<2;i++) for (j=0;j<2;j++) result[i][j]=max2[i][j]; // result is [2][3], max2 is [2][2]
		}
		/* Calculate locations of the maximums on FFT (corresponding to the diagonals of the checkerboard pattern) */
		double [][] maxOnFFT = {{2*size*(max2[0][0]-max2[1][0]),2*size*(max2[0][1]-max2[1][1])},
				{2*size*(max2[0][0]+max2[1][0]),2*size*(max2[0][1]+max2[1][1])}};
		/*  We have only one half of the FFT data  so rotate 180-degrees around the center if the point is in the bottom half*/
		double [] maxPhases=new double[2];
		boolean [] invertPhaseSign={false,false};
		int maxIndex; // iterate through the two maximums on FFT for phase measurement
		int [][]      interpolateXY= new int [2][2];
		double [][]   interpolateKxy=new double [2][2];
		double [][][] interpolatePhases=new double [2][2][2];
		double [][][] interpolateAmplitudes=new double [2][2][2];// Maybe use it? if the amplitudes are very different?
		int ix,iy;
		boolean phaseCorr; // phase shift before averaging, to prevent rollover
		for (maxIndex=0;maxIndex<2;maxIndex++) {
			if (maxOnFFT[maxIndex][1]<0) {
				invertPhaseSign[maxIndex]=true;
				maxOnFFT[maxIndex][0]=-maxOnFFT[maxIndex][0];
				maxOnFFT[maxIndex][1]=-maxOnFFT[maxIndex][1];
			}
			interpolateXY [maxIndex][1] = (int) maxOnFFT[maxIndex][1];
			interpolateKxy[maxIndex][1] = maxOnFFT[maxIndex][1] - interpolateXY [maxIndex][1];
			if (maxOnFFT[maxIndex][0]<0) {
				interpolateXY  [maxIndex][0] = (int) (size+ maxOnFFT[maxIndex][0]);
				interpolateKxy [maxIndex][0] = (size+ maxOnFFT[maxIndex][0]) - interpolateXY[maxIndex][0];
			} else {
				interpolateXY  [maxIndex][0] = (int) maxOnFFT[maxIndex][0];
				interpolateKxy [maxIndex][0] = maxOnFFT[maxIndex][0] - interpolateXY[maxIndex][0];
			}
			for (j=0;j<2;j++) {
				ix=(interpolateXY[maxIndex][0]+j) % size;
				for (i=0;i<2;i++) {
					iy=interpolateXY[maxIndex][1]+i;
					interpolateAmplitudes[maxIndex][i][j]= Math.sqrt(fft_complex[iy][ix][0]*fft_complex[iy][ix][0]+fft_complex[iy][ix][1]*fft_complex[iy][ix][1]);
					interpolatePhases[maxIndex][i][j]= Math.atan2((invertPhaseSign[maxIndex]?-1.0:1.0)*fft_complex[iy][ix][1], fft_complex[iy][ix][0]);
					if (this.debugLevel>5) {
						System.out.println("maxIndex="+maxIndex+" ix="+ix+" iy="+iy+" phase="+IJ.d2s(interpolatePhases[maxIndex][i][j],4)+" amplitude="+interpolateAmplitudes[maxIndex][i][j]);
					}
				}
			}
			phaseCorr=false;
			if ((interpolatePhases[maxIndex][0][0]> Math.PI/2) ||
					(interpolatePhases[maxIndex][0][1]> Math.PI/2) ||
					(interpolatePhases[maxIndex][1][0]> Math.PI/2) ||
					(interpolatePhases[maxIndex][1][1]> Math.PI/2) ||
					(interpolatePhases[maxIndex][0][0]<-Math.PI/2) ||
					(interpolatePhases[maxIndex][0][1]<-Math.PI/2) ||
					(interpolatePhases[maxIndex][1][0]<-Math.PI/2) ||
					(interpolatePhases[maxIndex][1][1]<-Math.PI/2)) {
				phaseCorr=true;
				interpolatePhases[maxIndex][0][0]+= (interpolatePhases[maxIndex][0][0]<0)?Math.PI:-Math.PI;
				interpolatePhases[maxIndex][0][1]+= (interpolatePhases[maxIndex][0][1]<0)?Math.PI:-Math.PI;
				interpolatePhases[maxIndex][1][0]+= (interpolatePhases[maxIndex][1][0]<0)?Math.PI:-Math.PI;
				interpolatePhases[maxIndex][1][1]+= (interpolatePhases[maxIndex][1][1]<0)?Math.PI:-Math.PI;
				if (this.debugLevel>5) {
					System.out.println("Shifting phases by PI/2 before averaging (to avoid rollover)");
				}
			}

			maxPhases[maxIndex]=       interpolateKxy[maxIndex][1] *(interpolateKxy[maxIndex][0] * interpolatePhases[maxIndex][1][1] + (1.0-interpolateKxy[maxIndex][0])* interpolatePhases[maxIndex][1][0])+
					(1.0 - interpolateKxy[maxIndex][1])*(interpolateKxy[maxIndex][0] * interpolatePhases[maxIndex][0][1] + (1.0-interpolateKxy[maxIndex][0])* interpolatePhases[maxIndex][0][0]);
			if (phaseCorr) maxPhases[maxIndex]+=(maxPhases[maxIndex]<0)?Math.PI:-Math.PI;
			if (this.debugLevel>5) {
				System.out.println("kx="+IJ.d2s(interpolateKxy [maxIndex][0],4)+ " ky="+IJ.d2s(interpolateKxy [maxIndex][1],4));
			}
			if (this.debugLevel>2) {
				System.out.println("maxIndex="+maxIndex+" phase="+IJ.d2s(maxPhases[maxIndex],4));
			}
		}
		double [] checkerPhases= findCheckerPhases(max2, maxPhases); /* may be different for greens==true . No, the same */
		for (i=0;i<2;i++) result[i][2]=checkerPhases[i];
		if (this.debugLevel>2)  System.out.println();
		return result;
	}

	// Above is currently used for years, below is for experiments

	public double[][] findPattern(
			DoubleFHT doubleFHT,
			double [] input_pixels, // pixel array to process
			int size, // FFT size
			PatternDetectParameters patternDetectParameters,
			double    min_half_period,
			double    max_half_period,
			boolean greens, // this is a pattern for combined greens (diagonal), adjust results accordingly
			String title) { // title prefix to use for debug  images

		if (!patternDetectParameters.use_large_cells) {
			return findPatternUsedNow(
					input_pixels, // pixel array to process
					size, // FFT size
					patternDetectParameters,
					greens, // this is a pattern for combined greens (diagonal), adjust results accordingly
					title); // title prefix to use for debug  images
		}
		int debug_threshold = 4;

//		this.debugLevel = 11;
		if (this.debugLevel > 5) {
			System.out.println("findPattern(): this.debugLevel = "+this.debugLevel);
		}
		if (doubleFHT == null) doubleFHT = new DoubleFHT(); // pass from the caller to re-use
		double [] cpixels=input_pixels.clone();
		if (min_half_period < 2) {
			min_half_period = 2;
		}
		double [] dfht = new double [input_pixels.length]; // will save FHT from phase auto correlation
		double [][]rslt= new double [2][3];
		cpixels = doubleFHT.phaseCorrelate (cpixels, patternDetectParameters.phaseCoeff, 0, patternDetectParameters.lowpass_sigma, dfht);
		if (this.debugLevel > (debug_threshold+1)) {
			double [][] src_corr = {input_pixels, cpixels};
			String [] titles = {"source","corr"};
			(new ShowDoubleFloatArrays()).showArrays(src_corr,size,size,true,"fft_corr-"+patternDetectParameters.phaseCoeff,titles);
		}

		int [][] half_periods = findFirst2MinOnCorrelation( cpixels, min_half_period, max_half_period);
		if (half_periods != null) { // check vectors are not close to colinear
			double sin =  (half_periods[0][0]*half_periods[1][1]-half_periods[0][1]*half_periods[1][0])/
					Math.sqrt((half_periods[0][0]*half_periods[0][0] + half_periods[0][1]*half_periods[0][1])*
							(half_periods[1][0]*half_periods[1][0] + half_periods[1][1]*half_periods[1][1]));
			if (Math.abs(sin) < patternDetectParameters.min_sin) {
				if (this.debugLevel > debug_threshold) System.out.println("Half-period vectors on correlation are too close to colinear - sin="+sin);
				return null;
			}

		}
		if (this.debugLevel>7) {
			if (half_periods == null) {
				if (this.debugLevel > (debug_threshold+0)) System.out.println("Could not find half-periods on correlation");
				return null;
			} else  {
				if (this.debugLevel > (debug_threshold+1)) System.out.println("first half period (x/y) :"+half_periods[0][0]+"/"+half_periods[0][1]+
						", second half period :"+half_periods[1][0]+"/"+half_periods[1][1]);
			}
		}

		double [][] dhp  = getPeriodsFromCorrelation(
				cpixels,
				half_periods, // a pair of half-periods
				patternDetectParameters.min_frac,
				patternDetectParameters.no_crazy) ;
		if (dhp != null) {
			if (this.debugLevel > (debug_threshold+1)) System.out.println(String.format("findPattern(): dhp= [[%.5f,%.5f],[%.5f,%.5f]]",dhp[0][0],dhp[0][1],dhp[1][0],dhp[1][1]));
		} else {
			if (this.debugLevel > (debug_threshold+0)) System.out.println("findPattern():  getPeriodsFromCorrelation() FAILED");
		}

		/* these are combined greens, convert vectors to original pixel space)! */
		if (greens) {
			double [][] rotMatrix= {{1.0,-1.0},{1.0,1.0}};
			double [][] max2orig= matrix2x2_mul(dhp,rotMatrix);
			for (int i=0; i<2; i++) for (int j=0; j<2; j++) rslt[i][j]=max2orig[i][j]; // result is [2][3], max2orig is [2][2]
		} else {
			for (int i=0;i<2;i++) for (int j=0;j<2;j++) rslt[i][j]=2.0 * dhp[i][j]; // result is [2][3], max2 is [2][2]
		}
		if (this.debugLevel > (debug_threshold + 0)) {
			System.out.println("Corrected to original pixels[0]  x="+IJ.d2s(rslt[0][0],4)+" y="+IJ.d2s(rslt[0][1],4));
			System.out.println("Corrected to original pixels[1]  x="+IJ.d2s(rslt[1][0],4)+" y="+IJ.d2s(rslt[1][1],4));
		}



		// Calculate complex FFT - needed to determine pattern phases (autocorrelation does not provide this
		// dfht array has FHT(input_pixels) save from calculation of the phase autocorrelation

		double [][][] fft_cmplx = doubleFHT.FHT2FFTHalf(dfht,size);
		if (this.debugLevel > (debug_threshold+2)) {
			(new ShowDoubleFloatArrays()).showComplex(fft_cmplx,"fft_cmplx");
		}

		/* Calculate locations of the maximums on FFT (corresponding to the diagonals of the checkerboard pattern) */
		double [][] max_on_FFT = {{2*size*(dhp[0][0]-dhp[1][0]),2*size*(dhp[0][1]-dhp[1][1])},
				{2*size*(dhp[0][0]+dhp[1][0]),2*size*(dhp[0][1]+dhp[1][1])}};
		/*  We have only one half of the FFT data  so rotate 180-degrees around the center if the point is in the bottom half*/
		double [] max_phases =    getPatternPhasesFromFFT(fft_cmplx, size, max_on_FFT);
		double [] checker_phases= findCheckerPhases(dhp, max_phases); /* may be different for greens==true . No, the same */
		for (int i=0;i<2;i++) rslt[i][2]=checker_phases[i];
		if (this.debugLevel > (debug_threshold + 0)) System.out.println();
		return rslt;
	}
	/* ======================================================================== */
	// calculate pattern 3D phases by interpolating im/re for the found maximums
	private double [] getPatternPhasesFromFFT(
			double [][][] fft_complex,
			int           size,
			double [][]   maxOnFFT) {
		double [] maxPhases=new double[2];
		boolean [] invertPhaseSign={false,false};
		int maxIndex; // iterate through the two maximums on FFT for phase measurement
		int [][]      interpolateXY= new int [2][2];
		double [][]   interpolateKxy=new double [2][2];
		double [][][] interpolatePhases=new double [2][2][2];
		double [][][] interpolateAmplitudes=new double [2][2][2];// Maybe use it? if the amplitudes are very different?
		int ix,iy;
		boolean phaseCorr; // phase shift before averaging, to prevent rollover
		for (maxIndex=0;maxIndex<2;maxIndex++) {
			if (maxOnFFT[maxIndex][1]<0) {
				invertPhaseSign[maxIndex]=true;
				maxOnFFT[maxIndex][0]=-maxOnFFT[maxIndex][0];
				maxOnFFT[maxIndex][1]=-maxOnFFT[maxIndex][1];
			}
			interpolateXY [maxIndex][1] = (int) maxOnFFT[maxIndex][1];
			interpolateKxy[maxIndex][1] = maxOnFFT[maxIndex][1] - interpolateXY [maxIndex][1];
			if (maxOnFFT[maxIndex][0]<0) {
				interpolateXY  [maxIndex][0] = (int) (size+ maxOnFFT[maxIndex][0]);
				interpolateKxy [maxIndex][0] = (size+ maxOnFFT[maxIndex][0]) - interpolateXY[maxIndex][0];
			} else {
				interpolateXY  [maxIndex][0] = (int) maxOnFFT[maxIndex][0];
				interpolateKxy [maxIndex][0] = maxOnFFT[maxIndex][0] - interpolateXY[maxIndex][0];
			}
			for (int j=0;j<2;j++) {
				ix=(interpolateXY[maxIndex][0]+j) % size;
				for (int i=0;i<2;i++) {
					iy=interpolateXY[maxIndex][1]+i; // next: OOB 2147483647
					interpolateAmplitudes[maxIndex][i][j]= Math.sqrt(fft_complex[iy][ix][0]*fft_complex[iy][ix][0]+fft_complex[iy][ix][1]*fft_complex[iy][ix][1]);
					interpolatePhases[maxIndex][i][j]= Math.atan2((invertPhaseSign[maxIndex]?-1.0:1.0)*fft_complex[iy][ix][1], fft_complex[iy][ix][0]);
					if (this.debugLevel>5) {
						System.out.println("maxIndex="+maxIndex+" ix="+ix+" iy="+iy+" phase="+IJ.d2s(interpolatePhases[maxIndex][i][j],4)+" amplitude="+interpolateAmplitudes[maxIndex][i][j]);
					}
				}
			}
			phaseCorr=false;
			if ((interpolatePhases[maxIndex][0][0]> Math.PI/2) ||
					(interpolatePhases[maxIndex][0][1]> Math.PI/2) ||
					(interpolatePhases[maxIndex][1][0]> Math.PI/2) ||
					(interpolatePhases[maxIndex][1][1]> Math.PI/2) ||
					(interpolatePhases[maxIndex][0][0]<-Math.PI/2) ||
					(interpolatePhases[maxIndex][0][1]<-Math.PI/2) ||
					(interpolatePhases[maxIndex][1][0]<-Math.PI/2) ||
					(interpolatePhases[maxIndex][1][1]<-Math.PI/2)) {
				phaseCorr=true;
				interpolatePhases[maxIndex][0][0]+= (interpolatePhases[maxIndex][0][0]<0)?Math.PI:-Math.PI;
				interpolatePhases[maxIndex][0][1]+= (interpolatePhases[maxIndex][0][1]<0)?Math.PI:-Math.PI;
				interpolatePhases[maxIndex][1][0]+= (interpolatePhases[maxIndex][1][0]<0)?Math.PI:-Math.PI;
				interpolatePhases[maxIndex][1][1]+= (interpolatePhases[maxIndex][1][1]<0)?Math.PI:-Math.PI;
				if (this.debugLevel>5) {
					System.out.println("Shifting phases by PI/2 before averaging (to avoid rollover)");
				}
			}

			maxPhases[maxIndex]=       interpolateKxy[maxIndex][1] *(interpolateKxy[maxIndex][0] * interpolatePhases[maxIndex][1][1] + (1.0-interpolateKxy[maxIndex][0])* interpolatePhases[maxIndex][1][0])+
					(1.0 - interpolateKxy[maxIndex][1])*(interpolateKxy[maxIndex][0] * interpolatePhases[maxIndex][0][1] + (1.0-interpolateKxy[maxIndex][0])* interpolatePhases[maxIndex][0][0]);
			if (phaseCorr) maxPhases[maxIndex]+=(maxPhases[maxIndex]<0)?Math.PI:-Math.PI;
			if (this.debugLevel>5) {
				System.out.println("kx="+IJ.d2s(interpolateKxy [maxIndex][0],4)+ " ky="+IJ.d2s(interpolateKxy [maxIndex][1],4));
			}
			if (this.debugLevel>2) {
				System.out.println("maxIndex="+maxIndex+" phase="+IJ.d2s(maxPhases[maxIndex],4));
			}
		}
		return maxPhases;
	}

	/* ======================================================================== */
	// zero is in the center
	private int [][] findFirst2MinOnCorrelation(
			double [] pixels,
			double    min_half_period, // from white to black
			double    max_half_period) { // from white to black
		double min_fract_period = 0.3; //second minimum should be not closer to the first than this fraction of the first
		// distance from zero
		int size = (int) Math.sqrt(pixels.length);
		int size2 = size >>1;
				if (max_half_period <= 0) max_half_period = size;
				if (min_half_period <  1) min_half_period = 1;
				double mn2 = min_half_period*min_half_period;
				double mx2 = max_half_period*max_half_period;
				int [][] mins={{0,0},{0,0}};
				double mn = 0.0; // minimums are always negative
				// find the first minimum
				for (int row = 0; row <= size2; row++) {
					double dy2 = (row-size2)*(row-size2);
					for (int col = 0; col < size; col++) {
						double d2 = (col-size2)*(col-size2) + dy2;
						if ((d2 >= mn2) && (d2 < mx2) && (pixels[row * size + col] < mn)) {
							mn = pixels[row * size + col];
							mins[0][0] = col-size2;
							mins[0][1] = row-size2; // always <=0
						}
					}
				}
				if (mn ==0.0) {
					return null; // could not find a first minimum
				}
				// find the second minimum (not too close to the first
				int x0 = mins[0][0]+size2;
				int y0 = mins[0][1]+size2;
				int x1 = size - x0;
				int y1 = size - y0;
				mn = 0.0;
				double mn2a = Math.max(mn2, min_fract_period*( mins[0][0]* mins[0][0]+mins[0][1]*mins[0][1]));

				for (int row = 0; row <= size2; row++) {
					double dy2 = (row-size2)*(row-size2);
					for (int col = 0; col < size; col++) {
						double d2 = (col-size2)*(col-size2) + dy2;
						if ((d2 >= mn2) && (d2 < mx2) && (pixels[row * size + col] < mn)) {
							// verify it is far enough from the first minimum and it's mirror
							d2 = (row - y0)*(row -y0)+(col - x0)*(col - x0);
							if (d2 >= 2*mn2a) { // diagonal
								// verify it is far enough from the mirrored first
								d2 = (row - y1)*(row -y1)+(col - x1)*(col - x1);
								if (d2 >= 2 * mn2a) { // diagonal
									mn = pixels[row * size + col];
									mins[1][0] = col-size2;
									mins[1][1] = row-size2; // always <=0
								}
							}
						}
					}
				}
				if (mn ==0.0) {
					return null; // could not find a second minimum
				}
				return mins;
	}
	private int rad_search(int [][]  half_periods) {
		int ahp;
		int amx;
		int dist = -1;
		for (int n = 0; n < half_periods.length; n++) {
			amx = 0;
			for (int i = 0; i < half_periods[0].length; i++) {
				ahp = half_periods[n][i];
				if (ahp < 0) ahp=-ahp;
				if (ahp > amx) amx = ahp;
			}
			if ((dist <0) || (dist > amx)) dist = amx;
		}
		if (dist < 3) return 1;
		return (dist -1) >> 1;
	}

	private boolean isLocalMinMax(double [] pixels, int size, int indx, boolean need_min) {
		int sizesize = size*size;
		if ((indx < 0) || (indx >= sizesize)){
			if (this.debugLevel > 0) {
				System.out.println("isLocalMinMax(): bad indx = "+indx);
			}
			return false;
		}
		for (int dy = -size; dy<=size; dy +=size) {
			for (int dx = -1; dx <=1; dx++) {
				int indx1 = indx + dy +dx;
				if ((indx1 >= 0) && (indx1 < sizesize)){
					if (need_min) {
						if (pixels[indx] > pixels[indx1]) return false;
					} else {
						if (pixels[indx] < pixels[indx1]) return false;
					}
				}
			}
		}
		return true;
	}

	private double [][] getPeriodsFromCorrelation(
			double [] pixels,
			int [][]  half_periods, // a pair of half-periods
			double    min_frac,     // minimal fraction of the center maximum to consider
			boolean   no_crazy ){   // discard min/max with failed interpolation
		// how far to look around expected maximum/minimum
		double [][] rslt = null;
		int dbg_threshold = 3;
		int      search_radius = rad_search (half_periods);
		// below first >= second
		int [][] search_order = {{1,0}, // just half periods, negative
				{1,1}, // diagonals, positive
				{2,0}, // full period, positive
				{2,1}, // 8 (4) directions, negative
				{2,2}}; // positive
		int size = (int) Math.sqrt(pixels.length);
		int size2 = size/2;
		int [][] aper = new int [2][2];
		double [][] dhalf_periods = new double [2][2]; // will refine  half_periods;
		for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) aper[i][j] = Math.abs(half_periods[i][j]);
		for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) dhalf_periods[i][j] =  half_periods[i][j];
		double min_avalue = pixels[size2*(size+1)]*min_frac; // minimal absolute value of min/max
		int max_dist = size2 - 2; //arbitrary ??
		int nord =  0; // to preserve after loop
		layers_iteration:
			for (nord = 0; nord < search_order.length; nord++) {
				int [] xy_ind= search_order[nord];
				// 1 check if too far to fit
				if (    ((xy_ind[0] * aper[0][0] + xy_ind[1] * aper[1][0]) > max_dist ) ||
						((xy_ind[0] * aper[0][1] + xy_ind[1] * aper[1][1]) > max_dist ) ||
						((xy_ind[1] * aper[0][0] + xy_ind[0] * aper[1][0]) > max_dist ) ||
						((xy_ind[1] * aper[0][1] + xy_ind[0] * aper[1][1]) > max_dist )
						) {
					break layers_iteration;
				}
				boolean are_mins = ((xy_ind[0] + xy_ind[1]) % 2) != 0;
				int nspots = ((xy_ind[1]==0) || (xy_ind[0]==xy_ind[1])) ? 2 : 4;
				//			double [][] predicted_centers = new double nsamples
				int [][] lin_comb = new int [nspots][2]; // number of half periods (incl. negative) to combine
				if (xy_ind[1]==0) {
					lin_comb[0][0] = xy_ind[0];  // first spot, first vector
					lin_comb[0][1] = 0;          // first spot, second vector
					lin_comb[1][0] = 0;          // second spot, first vector
					lin_comb[1][1] = xy_ind[0];  // second spot, second vector
				} else if (xy_ind[0]==xy_ind[1]){
					lin_comb[0][0] =  xy_ind[0]; // first spot, first vector
					lin_comb[0][1] =  xy_ind[0]; // first spot, second vector
					lin_comb[1][0] =  xy_ind[0]; // second spot, first vector
					lin_comb[1][1] = -xy_ind[0]; // second spot, second vector
				} else {
					lin_comb[0][0] =  xy_ind[0]; // first spot,  first vector
					lin_comb[0][1] =  xy_ind[1]; // first spot,  second vector
					lin_comb[1][0] =  xy_ind[0]; // second spot, first vector
					lin_comb[1][1] = -xy_ind[1]; // second spot, second vector
					lin_comb[2][0] =  xy_ind[1]; // third spot,  first vector
					lin_comb[2][1] =  xy_ind[0]; // third spot,  second vector
					lin_comb[3][0] = -xy_ind[1]; // fourth spot, first vector
					lin_comb[3][1] =  xy_ind[0]; // fourth spot, second vector
				}
				double [][] spot_centers = new double[nspots][2];
				// iterate through all min/max spots in a layer, break outer loop on any failure (and use previous layer data)
				for (int nspot = 0; nspot < nspots; nspot++) {
					double [] xyc = {
							lin_comb[nspot][0]*dhalf_periods[0][0] + lin_comb[nspot][1]*dhalf_periods[1][0],
							lin_comb[nspot][0]*dhalf_periods[0][1] + lin_comb[nspot][1]*dhalf_periods[1][1]};
					int [] ixyc = {(int) Math.round(xyc[0]),(int) Math.round(xyc[1])};
					if ((ixyc[0] < 0) || (ixyc[1] < 0) || (ixyc[0] >= size) || (ixyc[0] >= size)) {
						break layers_iteration;
					}
					int [] ixym = ixyc.clone();
					for (int y = ixyc[1] - search_radius; y <= ixyc[1] + search_radius; y ++ ) if ((y >= -max_dist) && (y < max_dist)) {
						for (int x = ixyc[0] - search_radius; x <= ixyc[0] + search_radius; x ++ ) if ((x >= -max_dist) && (x < max_dist)) {
							double d = pixels[(y+size2)*size+(x+size2)] - pixels[(ixym[1]+size2)*size+(ixym[0]+size2)];
							if (are_mins) d= -d;
							if (d > 0) {
								ixym[0] = x;
								ixym[1] = y;
							}
						}
					}
					// verify it is local max/min
					if (!isLocalMinMax(pixels, size, (ixym[1]+size2)*size+(ixym[0]+size2), are_mins)) {
						if (this.debugLevel > dbg_threshold) {
							System.out.println("getPeriodsFromCorrelation(): is not a local max/min, nord="+nord);
						}
						break layers_iteration;
					}
					// verify it is strong enough
					if (Math.abs(pixels[(ixym[1]+size2)*size+(ixym[0]+size2)]) < min_avalue) { // BUG: -8030
						if (this.debugLevel > dbg_threshold) {
							System.out.println("getPeriodsFromCorrelation(): min/max is too weak: "+
									pixels[(ixym[1]+size2)*size+(ixym[0]+size2)]+ " < "+ min_avalue+
									" (layer"+nord+")");
						}
						break layers_iteration;
					}

					//use quadratic interpolation to find min/max
					double [][][] data =new double [9][3][];
					int index=0;
					for (int dy = -1; dy <= 1; dy++) {
						for (int dx = -1; dx <= 1; dx++) {
							data[index][0] = new double [2];
							data[index][1]=     new double[1];
							data[index][2]=     new double[1];
							int pix_indx = (ixym[1] + dy + size2) * size + (ixym[0] + dx + size2);
							if ((pix_indx >= 0) && (pix_indx < pixels.length)) {
								data[index][0][0] = dx;
								data[index][0][1] = dy;
								data[index][1][0]=  are_mins ? -pixels[pix_indx] : pixels[pix_indx];
								data[index][2][0]=  1.0;
							} else {
								data[index][2][0]=  0.0;
							}
							index++;
						}
					}
					double [] corrXY=(new PolynomialApproximation()).quadraticMax2d (data);
					if (corrXY == null) {
						if (no_crazy) {
							if (this.debugLevel > dbg_threshold) {
								System.out.println("getPeriodsFromCorrelation(): interpolation failed,"+
										" (layer"+nord+")");
							}
							break layers_iteration;
						}
						corrXY=new double[2];
						corrXY[0] = 0.0;
						corrXY[1] = 0.0;
					}
					// verify that interpolation did not go crazy. If it did - reset to center
					if ((corrXY[0] > 1.5) || (corrXY[1] > 1.5) || (corrXY[0] < -1.5) || (corrXY[1] < -1.5)) {
						if (no_crazy) {
							if (this.debugLevel > dbg_threshold) {
								System.out.println("getPeriodsFromCorrelation(): interpolation failed, "+
										"corrXY = ["+corrXY[0]+", "+corrXY[1]+"]  (layer"+nord+")");
							}
							break layers_iteration;
						}
						corrXY[0] = 0.0;
						corrXY[1] = 0.0;
					}

					spot_centers[nspot][0] = ixym[0] + corrXY[0];
					spot_centers[nspot][1] = ixym[1] + corrXY[1];
				}
				// refine half-period values
				if (xy_ind[1]==0) {
					dhalf_periods[0][0] = spot_centers[0][0]/xy_ind[0];
					dhalf_periods[0][1] = spot_centers[0][1]/xy_ind[0];
					dhalf_periods[1][0] = spot_centers[1][0]/xy_ind[0];
					dhalf_periods[1][1] = spot_centers[1][1]/xy_ind[0];
				} else if (xy_ind[0]==xy_ind[1]){
					dhalf_periods[0][0] = (spot_centers[0][0] + spot_centers[1][0])/2/xy_ind[0];
					dhalf_periods[0][1] = (spot_centers[0][1] + spot_centers[1][1])/2/xy_ind[0];
					dhalf_periods[1][0] = (spot_centers[0][0] - spot_centers[1][0])/2/xy_ind[0];
					dhalf_periods[1][1] = (spot_centers[0][1] - spot_centers[1][1])/2/xy_ind[0];
				} else {
					dhalf_periods[0][0] = (spot_centers[0][0] + spot_centers[1][0] + spot_centers[2][0] - spot_centers[3][0]) / (xy_ind[0] + xy_ind[1])/2;
					dhalf_periods[0][1] = (spot_centers[0][1] + spot_centers[1][1] + spot_centers[2][1] - spot_centers[3][1]) / (xy_ind[0] + xy_ind[1])/2;
					dhalf_periods[1][0] = (spot_centers[0][0] - spot_centers[1][0] + spot_centers[2][0] + spot_centers[3][0]) / (xy_ind[0] + xy_ind[1])/2;
					dhalf_periods[1][1] = (spot_centers[0][1] - spot_centers[1][1] + spot_centers[2][1] + spot_centers[3][1]) / (xy_ind[0] + xy_ind[1])/2;
				}
			} // for (int nord = 0; nord < search_order.length; nord++) {
		//		double [][] pattern_diagonals = {
		//				{dhalf_periods[0][0]+dhalf_periods[1][0], dhalf_periods[0][1]+dhalf_periods[1][1]},
		//				{dhalf_periods[0][0]-dhalf_periods[1][0], dhalf_periods[0][1]-dhalf_periods[1][1]}};
		double [][] scaled_v = new double [2][2];
		double k_scale = 4.0;
		for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) {
			scaled_v[i][j] = k_scale*dhalf_periods[i][j]; // not clear why 4 and not 2
		}
		//		rslt =  vectToWaveVect(pattern_diagonals);
		// make vectors to have positive y:
		for (int i = 0; i < 2; i++) {
			if (scaled_v[i][1] <0) {
				scaled_v[i][0] = - scaled_v[i][0];
				scaled_v[i][1] = - scaled_v[i][1];
			}
		}
		rslt =  vectToWaveVect(scaled_v);
		//		rslt = dhalf_periods;
		// sort so first vector to second vector will be clockwise (positive y is downwards)
		{int j=0, k=1;
		if ((rslt[j][0] * rslt[k][1] - rslt[k][0] * rslt[j][1]) < 0) {
			double [] tmp = rslt[0];
			rslt[0] = rslt[1];
			rslt[1] = tmp;
		}
		}
		for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) {
			if (Double.isNaN(rslt[i][j]) || Double.isInfinite(rslt[i][j])) {
				System.out.println("got NaN/Infinity, nord = "+nord);
				return null;
			}
		}		return rslt;
	}

	//vect[0][] - 1-st vector connecting nodes, vect[1][] - second vector
	// return [0][] - first wave vector, orthogonal to vecgt[0][], return [1][] - to the second vect[1][]
	double [][] vectToWaveVect(double [][] vect){
		double a = 1.0/(vect[0][0]*vect[1][1] -vect[0][1]*vect[1][0]);
		//		double [] rv0 = {-vect[0][1],  vect[0][0]};
		//		double [] rv1 = { vect[1][1], -vect[1][0]};
		//		double [][] wv= {{a*rv0[0], a*rv0[1]},{a*rv1[0], a*rv1[1]}};
		double [][] wv= {{-a*vect[0][1], a*vect[0][0]}, {a*vect[1][1], -a*vect[1][0]}};
		return wv;
	}



	private	int [][] findFirst2MaxOnSpectrum (double [][][] pixels, // complex, top half, starting from 0,0
			/* May need to reduce the skip_around to be able to handle smaller number of pattern periods? Or re-try if failed? Guess somehow?*/
			int skip_around, // skip +- from (0,0) and previous max
			double minOrtho) { // 0.5 - 30deg. orthogonality of 2 vectors - 1.0 - perpendicular, 0.0 - parallel
		int [][] max2={{0,0},{0,0}};
		double thisMax=0.0;
		int x,y,sx1,sx2;
		double p,a;
		/* find first point */
		for (y=0;y<pixels.length;y++) for (x=0;x<pixels[0].length; x++ ) {
			p=pixels[y][x][0]*pixels[y][x][0]+pixels[y][x][1]*pixels[y][x][1];
			if (p>thisMax) {
				if ((y<=skip_around) && ((x<=skip_around) || (x>=pixels[0].length-skip_around))) continue; /* too close to [0,0] */
				max2[0][0]=x;
				max2[0][1]=y;
				thisMax=p;
			}
		}
		thisMax=0.0;
		sx1=(max2[0][0]>(pixels[0].length/2))?(max2[0][0]-pixels[0].length):max2[0][0]; /* y is always positive here */
		/* find second point */
		// Maybe also check if it is a local maximum (not on the border with protected area (0, first point, y=0, ...) ?
		for (y=0;y<pixels.length;y++) for (x=0;x<pixels[0].length; x++ ) {
			p=pixels[y][x][0]*pixels[y][x][0]+pixels[y][x][1]*pixels[y][x][1];
			if (p>thisMax) {
				/* Is this a local maximum? */
				if (y>0) {
					if (                             p< (pixels[y-1][x  ][0]*pixels[y-1][x  ][0]+pixels[y-1][x  ][1]*pixels[y-1][x  ][1]))  continue;
					if ((x>0) &&                    (p< (pixels[y-1][x-1][0]*pixels[y-1][x-1][0]+pixels[y-1][x-1][1]*pixels[y-1][x-1][1]))) continue;
					if ((x<(pixels[0].length-1)) && (p< (pixels[y-1][x+1][0]*pixels[y-1][x+1][0]+pixels[y-1][x+1][1]*pixels[y-1][x+1][1]))) continue;
				}
				if (y< (pixels.length-1)) {
					if (                             p< (pixels[y+1][x  ][0]*pixels[y+1][x  ][0]+pixels[y+1][x  ][1]*pixels[y+1][x  ][1]))  continue;
					if ((x>0) &&                    (p< (pixels[y+1][x-1][0]*pixels[y+1][x-1][0]+pixels[y+1][x-1][1]*pixels[y+1][x-1][1]))) continue;
					if ((x<(pixels[0].length-1)) && (p< (pixels[y+1][x+1][0]*pixels[y+1][x+1][0]+pixels[y+1][x+1][1]*pixels[y+1][x+1][1]))) continue;
				}
				if ((x>0) &&                      (p< (pixels[y  ][x-1][0]*pixels[y  ][x-1][0]+pixels[y  ][x-1][1]*pixels[y  ][x-1][1]))) continue;
				if ((x<(pixels[0].length-1)) &&   (p< (pixels[y  ][x+1][0]*pixels[y  ][x+1][0]+pixels[y  ][x+1][1]*pixels[y  ][x+1][1]))) continue;
				if ((y<=skip_around) && ((x<=skip_around) || (x>=pixels[0].length-skip_around))) {
					if (this.debugLevel>5) {
						System.out.println("rejecting point ["+x+","+y+"] it is too close to [0,0]");
					}
					continue; /* too close to [0,0] */
				}

				if ((y<=skip_around) && ((x<=skip_around) || (x>=pixels[0].length-skip_around))) {
					if (this.debugLevel>5) {
						System.out.println("rejecting point ["+x+","+y+"] it is too close to [0,0]");
					}
					continue; /* too close to [0,0] */
				}
				if (((y<=(max2[0][1]+skip_around)) && (y>=(max2[0][1]-skip_around))) &&
						(((x<=(max2[0][0]+skip_around)) && (x>=(max2[0][0]-skip_around))) ||
								(x>=(max2[0][0]+pixels[0].length-skip_around)) ||
								(x<=(max2[0][0]-pixels[0].length+skip_around)))) {
					if (this.debugLevel>5) {
						System.out.println("rejecting point ["+x+","+y+"] as it is too close to the first one - ["+max2[0][0]+"("+sx1+"),"+max2[0][1]+"]");
					}
					continue; /* too close to first maximum */
				}
				sx2=(x>(pixels[0].length/2))?(x-pixels[0].length):x;
				a=(sx1*y -max2[0][1]*sx2);
				a=a*a/(sx1*sx1+max2[0][1]*max2[0][1])/(sx2*sx2+y*y);
				if (a < (minOrtho*minOrtho)) { /* vectors are too close to parallel */
					if (this.debugLevel>5) {
						System.out.println("rejecting point ["+x+"("+sx2+"),"+y+"] as the vector is too close to the first one - ["+max2[0][0]+"("+sx1+"),"+max2[0][1]+"]");
						System.out.println("pixels.length="+pixels.length+" pixels[0].length="+pixels[0].length);
					}
					continue;
				}
				max2[1][0]=x;
				max2[1][1]=y;
				thisMax=p;
			}
		}
		if (  thisMax==0.0) {
			/*TODO: THat really happens on the real data */
			System.out.println("Failed to find a second maximum");
			return null;
		}
		if (max2[0][0]>(pixels[0].length/2)) max2[0][0]-=pixels[0].length;
		if (max2[1][0]>(pixels[0].length/2)) max2[1][0]-=pixels[0].length;
		return max2;
	}
	/* ======================================================================== */
	/* Can it handle negative y if the refined maximum goes there? (maximal value on positive Y) */
	private double[][] findFirst2MaxOnCorrelation(
			double [] pixels,
			int [][] startPoints,
			PatternDetectParameters patternDetectParameters
			) {
		double reasonbleFrequency=2.0; // reject frequencies below that
		int size =(int) Math.sqrt (pixels.length);
		int [][] imax =startPoints.clone();
		int [][] imax2 =new int [2*patternDetectParameters.multiplesToTry][2];
		boolean  []maxDefined=new boolean [2*patternDetectParameters.multiplesToTry];  //.multiplesToTry = 4;

		double  [] maxValues =new double [startPoints.length];
		double  [] max2Values =new double [2];
		double [][] max2 =new double [2*patternDetectParameters.multiplesToTry][2];
		int lim=size/2-2; // safety measure
		int nmax=0;
		int x=1;
		int y=0;
		int indx;
		int [] dirs = {-1, -size-1, -size, -size+1, 1, size+1, size, size-1};

		boolean isMax;
		int i,j,k, xmn,xmx,ymn,ymx;
		int halfSize=size/2;
		double [] vlengths={Math.sqrt(startPoints[0][0]*startPoints[0][0]+startPoints[0][1]*startPoints[0][1]),
				Math.sqrt(startPoints[1][0]*startPoints[1][0]+startPoints[1][1]*startPoints[1][1])};
		boolean tooFar=false;
		int numVect;

		/* Look for the maximal values around startPoints (+/-diff_spectr_corr )*/
		for (nmax=0; nmax<startPoints.length; nmax++) {
			ymn=imax[nmax][1]-patternDetectParameters.diffSpectrCorr; // .diffSpectrCorr=2
			ymx=imax[nmax][1]+patternDetectParameters.diffSpectrCorr;
			if (ymx>lim) ymx=lim;
			xmn=imax[nmax][0]-patternDetectParameters.diffSpectrCorr;
			if (xmn<-lim) xmx=-lim;
			xmx=imax[nmax][0]+patternDetectParameters.diffSpectrCorr;
			if (xmx>lim) xmx=lim;
			indx=(size+1)*size/2 + imax[nmax][1] *size+imax[nmax][0];
			if ((Math.abs(imax[nmax][0])>lim) || (Math.abs(imax[nmax][1])>lim)) {
				if (this.debugLevel>2) {
					System.out.println("Bad start point  imax[0][0]="+imax[0][0]+" imax[0][1]="+imax[0][1]+" imax[1][0]="+imax[1][0]+" imax[1][1]="+imax[1][1]+ " lim="+lim);
				}
				return null;
			}
			//       if ((indx<0) || (indx>=pixels.length)) {
			//          System.out.println(" imax[0][0]="+imax[0][0]+" imax[0][1]="+imax[0][1]+" imax[1][0]="+imax[1][0]+" imax[1][1]="+imax[1][1]+ " lim="+lim);
			//       }
			maxValues[nmax]=pixels[indx];
			for (y=ymn;y<=ymx;y++) for (x=xmn;x<=xmx; x++) {
				indx=(size+1)*size/2 + y *size+x;
				if (pixels[indx]>maxValues[nmax]) {
					/* Make sure it is closer to this point than to any other or [0,0] */
					tooFar=false;
					for (numVect=0;numVect<startPoints.length;numVect++) {
						if (Math.abs((x-imax[nmax][0])*imax[numVect][0]+(y-imax[nmax][1])*imax[numVect][1])>0.5*vlengths[nmax]*vlengths[numVect]) {
							tooFar=true;
							break;
						}
					}
					if (tooFar) {
						if (this.debugLevel>5) {
							System.out.println("rejecting point ["+x+","+y+"] as the vector is closer to other max than ["+imax[nmax][0]+","+imax[nmax][0]+"]"+
									" in the (+/-) direction: ["+imax[numVect][0]+","+imax[numVect][0]+"]");
						}
						continue;
					}
					maxValues[nmax]=pixels[indx];
					imax[nmax][0]=x;
					imax[nmax][1]=y;
				}
			}
			/* Make sure the maximum in the scanned area is also a local maximum */
			isMax=true;
			indx=(size+1)*size/2 + imax[nmax][1] *size+imax[nmax][0];
			for (j=0;j<7;j++) if (pixels[indx]<pixels[indx+dirs[j]]) {
				isMax=false;
				break;
			}
			if (!isMax) {
				if (this.debugLevel>2) {
					System.out.println("This should not happen:");
					System.out.println("Maximum is not a local maximum - BUG or consider changing patternDetectParameters.diffSpectrCorr="+patternDetectParameters.diffSpectrCorr);
					System.out.println("point #"+(nmax+1)+" (of 2), x0="+startPoints[nmax][0]+" y0="+startPoints[nmax][1]+ " x="+imax[nmax][0]+" y="+imax[nmax][1]);
				}
				/* Maybe return from here with null?*/
				while (!isMax && (indx>0) && (indx>pixels.length)) {
					isMax=true;
					for (j=0;j<7;j++) if (pixels[indx]<pixels[indx+dirs[j]]) {
						isMax=false;
						indx+=dirs[j];
						imax[nmax][1]=(indx / size) - halfSize;
						imax[nmax][0]=(indx % size) - halfSize;
						break;
					}
				}
				if (!isMax) {
					if (this.debugLevel>2) {
						System.out.println("Maximum still not reached, bailing out");
						System.out.println("point #"+(nmax+1)+" (of 2), x0="+startPoints[nmax][0]+" y0="+startPoints[nmax][1]+ " x="+imax[nmax][0]+" y="+imax[nmax][1]);
					}
					return null;
				} else {
					if (this.debugLevel>2) {
						System.out.println("point #"+(nmax+1)+" (of 2), corrected local maximum is  x="+imax[nmax][0]+" y="+imax[nmax][1]);
					}
				}

			}
		}
		/* Sort maximums so first vector to second vector will be clockwise (positive y is downwards) */
		j=0; k=1;
		if ((imax[j][0]*imax[k][1]-imax[k][0]*imax[j][1])<0) {
			j=1;k=0;
		}
		imax2[0][0]=imax[j][0];
		imax2[0][1]=imax[j][1];
		imax2[1][0]=imax[k][0];
		imax2[1][1]=imax[k][1];


		/* Now define maximal radius of cluster (~0.7 of the average distance from 0,0 to the 2 start points */
		int maxX2Y2=0;
		for (i=0;i<2;i++) for (j=0;j<2;j++) maxX2Y2+=imax2[i][j]*imax2[i][j];
		maxX2Y2/=4;
		//System.out.println("maxX2Y2="+maxX2Y2);

		for (i=0;i<2*patternDetectParameters.multiplesToTry;i++)  maxDefined[i]=(i<2);

		nmax=2*patternDetectParameters.multiplesToTry; /* but only the first two are known by now */
		for (i=0;i<2;i++) {
			if (this.debugLevel>5) System.out.println("i="+i+" x="+imax2[i][0]+" y="+imax2[i][1]+" value="+max2Values[i]);
		}

		int clusterNumber;
		List <Integer> pixelList=new ArrayList<Integer>(100);
		Integer Index, NewIndex, HillIndex;
		double cx,cy,cm,minInCluster,f;

		int []clusterMap=new int[pixels.length];
		for (i=0;i<clusterMap.length;i++) clusterMap[i]=0; /// 0 - unused, -1 - "do not use"
		int listIndex;
		boolean noHills;
		int clusterSize;
		boolean isLocalMax;
		int pair;
		for (clusterNumber=0;clusterNumber<nmax;clusterNumber++) {
			pair=clusterNumber/2+1;
			if (!maxDefined[clusterNumber] && maxDefined[clusterNumber-2]) {
				/* We do not know the seed for this maximum, but the previous (of the same direction) may be known */
				x= (int) (max2[clusterNumber-2][0]*pair)/(pair-1);
				y= (int) (max2[clusterNumber-2][1]*pair)/(pair-1);
				if ((x>(-lim+1)) && (x<(lim-1)) && (y>(-lim+1)) && (y<(lim-1))) {
					Index=(y+halfSize)*size + x+halfSize;
					/* there should be local maximum not more than "deviationSteps" steps from the x,y */
					isLocalMax=false;
					i=patternDetectParameters.deviationSteps;
					while ((i>0) && !isLocalMax) {
						isLocalMax=true;
						for (j=0;j<dirs.length;j++) {
							NewIndex=Index+dirs[j];
							if ((NewIndex>=0) && (NewIndex<clusterMap.length) && (pixels[NewIndex]>pixels[Index])) {
								isLocalMax=false;
								Index=NewIndex;
								i--;
								if (this.debugLevel>5) System.out.println("i="+i+" x="+((Index % size) - halfSize)+" y="+((Index / size) - halfSize)+" value="+pixels[Index]);
								break;
							}
						}
					}
					if (isLocalMax && (clusterMap[Index]==0)) { // not yet used
						imax2[clusterNumber][0]= (Index % size) - halfSize;
						imax2[clusterNumber][1]= (Index / size) - halfSize;
						maxDefined[clusterNumber]=true;
					}
				}
			}
			if (maxDefined[clusterNumber]) { // skip if seed for the cluster is not defined
				/* Grow cluster around maximum, find centroid */
				Index=(imax2[clusterNumber][1]+halfSize)*size + imax2[clusterNumber][0]+halfSize;
				pixelList.clear();
				pixelList.add (Index);
				clusterMap[Index]=clusterNumber+1;
				listIndex=0;
				while (listIndex<pixelList.size() ) {
					Index=pixelList.get(listIndex++);
					for (j=0;j<dirs.length;j++) {
						NewIndex=Index+dirs[j];
						if ((NewIndex>=0) && (NewIndex<clusterMap.length) && (clusterMap[NewIndex]==0) && (pixels[NewIndex]<pixels[Index])) {
							/* did we get too far?*/
							y=(NewIndex/size) - halfSize - imax2[clusterNumber][1];
							x=(NewIndex % size) - halfSize - imax2[clusterNumber][0];
							//System.out.println(" dy="+y+" dx="+x+" dx*dx+dy*dy="+(x*x+y*y));
							if ((x*x+y*y) <= maxX2Y2) {
								/* See if there is any neighbor of the new pixel that is higher and not yet marked (prevent rivers flowing between hills) */
								noHills=true;
								for (k=0;k<dirs.length;k++) {
									HillIndex=NewIndex+dirs[k];
									if ((HillIndex>=0) && (HillIndex<clusterMap.length) && (clusterMap[HillIndex]!=(clusterNumber+1)) && (pixels[HillIndex]>pixels[NewIndex])) {
										noHills=false;
										break;
									}
								}
								if (noHills) {
									pixelList.add (NewIndex);
									clusterMap[NewIndex]=clusterNumber+1;
									//System.out.println("NewIndex="+NewIndex+" y="+(NewIndex/size - halfSize)+" x="+((NewIndex % size) - halfSize)+" new pixel="+pixels[NewIndex]+" old pixel="+pixels[Index]);
								}
							}
						}
					}
				}
				/* Shrink clusters to a fraction of initial size */
				// TODO: shring to a value (between min and max) if there is a sharp maximum??
				//, double patternDetectParameters.shrinkClusters
				if (patternDetectParameters.shrinkClusters==0.0) { // use "smart" size
					clusterSize=(int) Math.sqrt(5* pixelList.size()); // use proportional size
				} else if (patternDetectParameters.shrinkClusters<0) {
					clusterSize=(int)(- patternDetectParameters.shrinkClusters ); // use specified size
				} else {
					clusterSize=(int) (pixelList.size()*patternDetectParameters.shrinkClusters); // use proportional size
				}
				if (clusterSize<5) clusterSize=5;
				while (pixelList.size()>clusterSize) {
					i=0;
					f=pixels[pixelList.get(i)];
					for (j=1;j<pixelList.size();j++) if (pixels[pixelList.get(j)]<f){
						i=j;
						f=pixels[pixelList.get(j)];
					}
					clusterMap[pixelList.get(i)]=-1; // Do not use looking for the next cluster
					pixelList.remove(i);
				}



				/* now find centroid of the cluster */
				minInCluster=pixels[pixelList.get(0)];
				for (i=1;i<pixelList.size();i++) if (minInCluster>pixels[pixelList.get(i)]) minInCluster = pixels[pixelList.get(i)];
				cx=0.0; cy=0.0; cm=0.0;
				for (i=0;i<pixelList.size();i++) {
					j=pixelList.get(i);
					y=j / size - halfSize;
					x=j % size - halfSize;
					f=pixels[j]-minInCluster;
					cm+=f;
					cx+=f*x;
					cy+=f*y;
				}
				cx/=cm;
				cy/=cm;
				max2[clusterNumber][0]=cx;
				max2[clusterNumber][1]=cy;
				f=0.0;
				if (pair>1) {
					cx=max2[clusterNumber-2][0]*pair/(pair-1)-max2[clusterNumber][0];
					cy=max2[clusterNumber-2][1]*pair/(pair-1)-max2[clusterNumber][1];
					f=Math.sqrt(cx*cx+cy*cy);
					/* Verify deviation here */
					if (f>patternDetectParameters.deviation) maxDefined[clusterNumber]=false;
				}
				if (this.debugLevel>6) System.out.println("pixelList.size()="+pixelList.size()+" centroid sum="+cm);
				if (this.debugLevel>5) System.out.println("clusterNumber="+clusterNumber+" x="+max2[clusterNumber][0]+" y="+max2[clusterNumber][1] + " x0="+(max2[clusterNumber][0]/pair)+" y0="+(max2[clusterNumber][1]/pair)+" deviat="+f);
				if ((cm==0.0) || (pixelList.size()<3)) maxDefined[clusterNumber]=false;
				/* Filter out unreasonably low frequencies*/
				if ((max2[clusterNumber][0]*max2[clusterNumber][0]+max2[clusterNumber][1]*max2[clusterNumber][1])<(reasonbleFrequency*reasonbleFrequency)) {
					if (this.debugLevel>2) System.out.println("Frequency too low:clusterNumber="+clusterNumber+" x="+max2[clusterNumber][0]+" y="+max2[clusterNumber][1]+ ", minimal allowed frequency is "+reasonbleFrequency);
					maxDefined[clusterNumber]=false;
				}
			}
		}
		/* Average (or just use farthest?) multiple maximums */
		if (this.debugLevel>2){
			float [] dbg_pixels=new float[clusterMap.length];
			for (j=0;j<dbg_pixels.length;j++) dbg_pixels[j]=clusterMap[j];
			dbg_pixels[(size+1)*size/2]=-1; // mark center

			ImageProcessor ip=new FloatProcessor(size,size);
			ip.setPixels(dbg_pixels);
			ip.resetMinAndMax();
			ImagePlus imp=  new ImagePlus("clusters", ip);
			imp.show();
		}

		//    for (i=0;i<2;i++) for (j=0;j<2;j++) max2[i][j]=imax2[i][j];
		double [][]maxFinal=new double[2][2];
		boolean [] definedFinal={false,false};
		for (i=0;i<2;i++) {
			maxFinal[i][0]=0.0;
			maxFinal[i][1]=0.0;
			j=0;
			for (pair=0;pair<nmax/2;pair++) {
				if (maxDefined[i+2*pair]) {
					j++;
					maxFinal[i][0]+=max2[i+2*pair][0]/(pair+1);
					maxFinal[i][1]+=max2[i+2*pair][1]/(pair+1);
					definedFinal[i]=true;
					//          System.out.println("i="+i+" pair="+pair+" j="+j+" maxFinal["+i+"][0]="+maxFinal[i][0]+" maxFinal["+i+"][1]="+maxFinal[i][1]);
				}
			}
			maxFinal[i][0]/=j;
			maxFinal[i][1]/=j;
			if (j==0) definedFinal[i]=false;
			/* These two vectors correspond to checker diagonals and they are calculated for the correlation space.
 Actual frequencies for the checker board will be 1/4 of these, also divide by FFT size so the result will be in cycles per pixel */
			maxFinal[i][0]/=size*4;
			maxFinal[i][1]/=size*4;
		}
		if (this.debugLevel>2) {
			System.out.println("Checkerboard frequency[0]  x="+IJ.d2s(maxFinal[0][0],4)+" y="+IJ.d2s(maxFinal[0][1],4));
			System.out.println("Checkerboard frequency[1]  x="+IJ.d2s(maxFinal[1][0],4)+" y="+IJ.d2s(maxFinal[1][1],4));
			//      System.out.println();
		}
		if (!definedFinal[0] || !definedFinal[1]) {
			if (this.debugLevel>2) {
				System.out.println("Undefined frequency(ies)");
			}
			return null;
		}
		return maxFinal;
	}
	/* ======================================================================== */
	private double [] findCheckerPhases(double [][] WVectors, double [] P) {
		double [][] DWVectors = {{WVectors[0][0]-WVectors[1][0],WVectors[0][1]-WVectors[1][1]},
				{WVectors[0][0]+WVectors[1][0],WVectors[0][1]+WVectors[1][1]}};
		if (this.debugLevel>3)  System.out.println("      DWVectors[0][0]="+IJ.d2s(DWVectors[0][0],4)+"  DWVectors[0][1]="+IJ.d2s(DWVectors[0][1],4));
		if (this.debugLevel>3)  System.out.println("      DWVectors[1][0]="+IJ.d2s(DWVectors[1][0],4)+"  DWVectors[1][1]="+IJ.d2s(DWVectors[1][1],4));
		double [] DWVectorsAbs2={DWVectors[0][0]*DWVectors[0][0]+DWVectors[0][1]*DWVectors[0][1],
				DWVectors[1][0]*DWVectors[1][0]+DWVectors[1][1]*DWVectors[1][1]};
		if (this.debugLevel>3)  System.out.println("      sqrt(DWVectorsAbs2[0])="+IJ.d2s(Math.sqrt(DWVectorsAbs2[0]),4)+"  sqrt(DWVectorsAbs2[1])="+IJ.d2s(Math.sqrt(DWVectorsAbs2[1]),4));
		double [][] DL= {{DWVectors[0][0]/DWVectorsAbs2[0],DWVectors[0][1]/DWVectorsAbs2[0]},
				{DWVectors[1][0]/DWVectorsAbs2[1],DWVectors[1][1]/DWVectorsAbs2[1]}};
		if (this.debugLevel>3)  System.out.println("      DL[0][0]="+IJ.d2s(DL[0][0],4)+"  DL[0][1]="+IJ.d2s(DL[0][1],4));
		if (this.debugLevel>3)  System.out.println("      DL[1][0]="+IJ.d2s(DL[1][0],4)+"  DL[1][1]="+IJ.d2s(DL[1][1],4));

		double v= (DL[0][0]*(DL[1][0]*P[1]-DL[0][0]*P[0]) +  DL[0][1]*(DL[1][1]*P[1]-DL[0][1]*P[0]))/(DL[0][0]*DL[1][1]-DL[0][1]*DL[1][0])/(2*Math.PI);
		if (this.debugLevel>2)  System.out.println("v="+IJ.d2s(v,4));


		double [] WC={-DL[1][0] * P[1] / (2*Math.PI) + DL[1][1] * v,
				-DL[1][1] * P[1] / (2*Math.PI) - DL[1][0] * v};
		if (this.debugLevel>2)  System.out.println("WC[0]="+IJ.d2s(WC[0],4)+"  WC[1]="+IJ.d2s(WC[1],4));

		double [] phases={-2*Math.PI *(WC[0] * WVectors[0][0] + WC[1] * WVectors[0][1]),
				-2*Math.PI *(WC[0] * WVectors[1][0] + WC[1] * WVectors[1][1])};
		if (this.debugLevel>2)  System.out.println("phases[0]="+IJ.d2s(phases[0],4)+"  phases[1]="+IJ.d2s(phases[1],4));
		return phases;
	}
	/* ======================================================================== */
	private	boolean matchPatterns(double [][][] patterns, double [][] targetPattern) {
		int n,i,j;
		double [][] sp=new double[2][2];
		double [] target_lengths= new double[2];
		double [] current_lengths=new double[2];
		double [] swap_wv;
		double    swap;
		boolean noCorrectionWasNeeded=true;
		if (patterns==null) return false;
		if (targetPattern==null) return false;
		for (n=0;n<patterns.length; n++) {
			for (i=0;i<2;i++) {
				target_lengths[i]= Math.sqrt(targetPattern[i][0]*targetPattern[i][0]+targetPattern[i][1]*targetPattern[i][1]);
				current_lengths[i]=Math.sqrt(patterns[n][i][0]*patterns[n][i][0]+patterns[n][i][1]*patterns[n][i][1]);

			}
			for (i=0;i<2;i++) for (j=0;j<2;j++){
				sp[i][j]=    (patterns[n][i][0]*targetPattern[j][0]+patterns[n][i][1]*targetPattern[j][1])/current_lengths[i]/ target_lengths[j];
			}
			/* Swap vectors regardless of sign (parallel/anti-parallel)*/
			if ((Math.abs(sp[0][0])<Math.abs(sp[0][1])) || (Math.abs(sp[0][0])<Math.abs(sp[1][0]))) {
				noCorrectionWasNeeded=false;
				if (this.debugLevel>0) System.out.println("Swapped wave vectors in quadrant "+n);
				swap_wv=patterns[n][0];   patterns[n][0]=patterns[n][1]; patterns[n][1]=swap_wv;
				swap=sp[0][0];            sp[0][0]=sp[1][0];             sp[1][0]=swap;
				swap=sp[0][1];            sp[0][1]=sp[1][1];             sp[1][1]=swap;
			}
			/* Now correct vector signs if needed */
			for (i=0;i<2;i++) {
				if (sp[i][i] <0) {
					noCorrectionWasNeeded=false;
					if (this.debugLevel>0) System.out.println("Changing wave vector "+(i+1)+" direction in quadrant "+n);
					for (j=0;j<patterns[n][i].length;j++) patterns[n][i][j]=-patterns[n][i][j]; /// Will negate phase if available
				}
			}
		}
		return noCorrectionWasNeeded;
	}
	/*	** converts 2 wave vectors (WVx,WVy,phase) into two checker pattern vectors (VPx,VPy, phase)
	phase is in the point x=0,y=0*/
	private	double [][] waveVectorsToPatternVectors(double [] wv0, double [] wv1) {
		double [][] v=new double [2][3];
		double vect_wv0_x_wv1=wv0[0]*wv1[1]-wv0[1]*wv1[0];
		//v[0][0]= wv0[1]/vect_wv0_x_wv1;
		//v[0][1]=-wv0[0]/vect_wv0_x_wv1;
		//v[1][0]=-wv1[1]/vect_wv0_x_wv1;
		//v[1][1]= wv1[0]/vect_wv0_x_wv1;

		v[0][0]= wv1[1]/vect_wv0_x_wv1;
		v[0][1]=-wv1[0]/vect_wv0_x_wv1;
		v[1][0]=-wv0[1]/vect_wv0_x_wv1;
		v[1][1]= wv0[0]/vect_wv0_x_wv1;
		if (wv0.length>2) v[0][2]=wv0[2];
		if (wv1.length>2) v[1][2]=wv1[2];
		/* "white square" center had coordinates
	-wv0[0]
		 */
		return v;
	}
	/**
	  Matching non-linear mesh using wave vectors in the centers of four quadrants, phases are ignored here
	  processes wave vectors for the four quadrants (top-left,tr,bl,br) and calculates second degree polynominal correction
	  X1= Ax*X^2 + Bx*Y^2 + 2Cx*X*Y
	  Y1= Ay*X^2 + By*Y^2 + 2Cy*X*Y
	  where X and Y are normalized so top left corner is (-1,-1), top right - (+1,-1) , bottom left - (-1,+1), bottom right - (+1,+1)
	  returns array of 6 elements {Ax,Bx,Cx,Ay,By,Cy}
	 */
	private	double [] calcPatternNonLinear(double [][][] qp) {
		int iq;
		/* Calculate center WV */

		double [][][] mCorners=    new double [4][][] ; // only two 2x2 are needed, other two - just for debugging
		mCorners[0]= matrixToConvertTwoPairs(qp[0][0], qp[0][1], qp[4][0], qp[4][1]);
		mCorners[1]= matrixToConvertTwoPairs(qp[1][0], qp[1][1], qp[4][0], qp[4][1]);
		mCorners[2]= matrixToConvertTwoPairs(qp[2][0], qp[2][1], qp[4][0], qp[4][1]);
		mCorners[3]= matrixToConvertTwoPairs(qp[3][0], qp[3][1], qp[4][0], qp[4][1]);

		/**
	x,y - image coordinates (distorted), with no center shift (later use effective radius)
	x1,y1 - linear mesh

	x1=x+Ax*x^2+Bx*y^2+2*Cx*x*y+Dx*x+Ex*y
	y1=y+Ay*x^2+By*y^2+2*Cy*x*y+Dy*x+Ey*y

	dx1/dx=1+2*Ax*x+2*Cx*y+Dx
	dx1/dy=0+2*Cx*x+2*Bx*y+Ex

	dy1/dx=0+2*Ay*x+2*Cy*y+Dy
	dy1/dy=1+2*Cy*x+2*By*y+Ey

	| dx1 | = |  1+2*Ax*x+2*Cx*y+Dx     2*Cx*x+2*Bx*y+Ex | * | dx |
	| dy1 |   |    2*Ay*x+2*Cy*y+Dy   1+2*Cy*x+2*By*y+Ey |   | dy |

	for x=-1,y=-1
	M[0]=|  1-2*Ax-2*Cx+Dx    -2*Cx-2*Bx+Ex |
	   |   -2*Ay-2*Cy+Dy   1-2*Cy-2*By+Ey |

	for x=+1,y=-1
	M[1]=|  1+2*Ax-2*Cx+Dx    +2*Cx-2*Bx+Ex |
	   |   +2*Ay-2*Cy+Dy   1+2*Cy-2*By+Ey |

	for x=-1,y=+1
	M[2]=|  1-2*Ax+2*Cx+Dx    -2*Cx+2*Bx+Ex |
	   |   -2*Ay+2*Cy+Dy   1-2*Cy+2*By+Ey |

	for x=+1,y=+1
	M[2]=|  1+2*Ax+2*Cx+Dx    +2*Cx+2*Bx+Ex |
	   |   +2*Ay+2*Cy+Dy   1+2*Cy+2*By+Ey |


		 */
		double [] det=  {mCorners[0][0][0]*mCorners[0][1][1]-mCorners[0][0][1]*mCorners[0][1][0],
				mCorners[1][0][0]*mCorners[1][1][1]-mCorners[1][0][1]*mCorners[1][1][0],
				mCorners[2][0][0]*mCorners[2][1][1]-mCorners[2][0][1]*mCorners[2][1][0],
				mCorners[3][0][0]*mCorners[3][1][1]-mCorners[3][0][1]*mCorners[3][1][0]};

		double [][][] M={{{ mCorners[0][1][1]/det[0], -mCorners[0][1][0]/det[0]},
			{-mCorners[0][0][1]/det[0],  mCorners[0][0][0]/det[0]}},
				{{ mCorners[1][1][1]/det[1], -mCorners[1][1][0]/det[1]},
				{-mCorners[1][0][1]/det[1],  mCorners[1][0][0]/det[1]}},
				{{ mCorners[2][1][1]/det[2], -mCorners[2][1][0]/det[2]},
					{-mCorners[2][0][1]/det[2],  mCorners[2][0][0]/det[2]}},
				{{ mCorners[3][1][1]/det[3], -mCorners[3][1][0]/det[3]},
						{-mCorners[3][0][1]/det[3],  mCorners[3][0][0]/det[3]}}};
		/**
	 Overdefined - 16 equations for 10 unknowns
	1   M[0][0][0]=1-2*Ax-2*Cx+Dx
	2   M[0][0][1]= -2*Cx-2*Bx+Ex
	3   M[0][1][0]= -2*Ay-2*Cy+Dy
	4   M[0][1][1]=1-2*Cy-2*By+Ey
	5   M[1][0][0]=1+2*Ax-2*Cx+Dx
	6   M[1][0][1]=  2*Cx-2*Bx+Ex
	7   M[1][1][0]=  2*Ay-2*Cy+Dy
	8   M[1][1][1]=1+2*Cy-2*By+Ey

	9   M[2][0][0]=1-2*Ax+2*Cx+Dx
	10   M[2][0][1]= -2*Cx+2*Bx+Ex
	11   M[2][1][0]= -2*Ay+2*Cy+Dy
	12   M[2][1][1]=1-2*Cy+2*By+Ey
	13   M[3][0][0]=1+2*Ax+2*Cx+Dx
	14   M[3][0][1]=  2*Cx+2*Bx+Ex
	15   M[3][1][0]=  2*Ay+2*Cy+Dy
	16   M[3][1][1]=1+2*Cy+2*By+Ey

	( 1)   M[0][0][0]=1-2*Ax     -2*Cx                +Dx
	( 2)   M[0][0][1]=      -2*Bx-2*Cx                      +Ex
	( 3)   M[0][1][0]=                -2*Ay     -2*Cy    +Dy
	( 4)   M[0][1][1]=1                    -2*By-2*Cy          +Ey
	( 5)   M[1][0][0]=1+2*Ax     -2*Cx                +Dx
	( 6)   M[1][0][1]=      -2*Bx+2*Cx                      +Ex
	( 7)   M[1][1][0]=                +2*Ay     -2*Cy    +Dy
	( 8)   M[1][1][1]=1                    -2*By+2*Cy          +Ey

	( 9)   M[0][0][0]=1-2*Ax     +2*Cx                +Dx
	(10)   M[0][0][1]=      +2*Bx-2*Cx                      +Ex
	(11)   M[0][1][0]=                -2*Ay     +2*Cy    +Dy
	(12)   M[0][1][1]=1                    +2*By-2*Cy          +Ey
	(13)   M[1][0][0]=1+2*Ax     +2*Cx                +Dx
	(14)   M[1][0][1]=      +2*Bx+2*Cx                      +Ex
	(15)   M[1][1][0]=                +2*Ay     +2*Cy    +Dy
	(16)   M[1][1][1]=1                    +2*By+2*Cy          +Ey

	( 1)+( 5)+( 9)+(13)                     Dx=( M[0][0][0]+M[1][0][0]+M[2][0][0]+M[3][0][0])/4-1
	( 2)+( 6)+(10)+(14)                     Ex=( M[0][0][1]+M[1][0][1]+M[2][0][1]+M[3][0][1])/4
	( 3)+( 7)+(11)+(15)                     Dy=( M[0][1][0]+M[1][1][0]+M[2][1][0]+M[3][1][0])/4
	( 4)+( 8)+(12)+(16)                     Ey=( M[0][1][1]+M[1][1][1]+M[2][1][1]+M[3][1][1])/4-1

	-( 1)+( 5)-( 9)+(13)                     Ax=(-M[0][0][0]+M[1][0][0]-M[2][0][0]+M[3][0][0])/8
	-( 2)-( 6)+(10)+(14)                     Bx=(-M[0][0][1]-M[1][0][1]+M[2][0][1]+M[3][0][1])/8
	-( 2)+( 6)-(10)+(14)-( 1)-( 5)+( 9)+(13) Cx=(-M[0][0][1]+M[1][0][1]-M[2][0][1]+M[3][0][1]-M[0][0][0]-M[1][0][0]+M[2][0][0]+M[3][0][0])/16
	-( 3)+( 7)-(11)+(15)                     Ay=(-M[0][1][0]+M[1][1][0]-M[2][1][0]+M[3][1][0])/8
	-( 4)-( 8)+(12)+(16)                     By=(-M[0][1][1]-M[1][1][1]+M[2][1][1]+M[3][1][1])/8
	-( 4)+( 8)-(12)+(16)-( 3)-( 7)+(11)+(15) Cy=(-M[0][1][1]+M[1][1][1]-M[2][1][1]+M[3][1][1]-M[0][1][0]-M[1][1][0]+M[2][1][0]+M[3][1][0])/16
		 */

		double [] rslt = {(-M[0][0][0]+M[1][0][0]-M[2][0][0]+M[3][0][0])/8, // Ax
				(-M[0][0][1]-M[1][0][1]+M[2][0][1]+M[3][0][1])/8, // Bx
				(-M[0][0][1]+M[1][0][1]-M[2][0][1]+M[3][0][1]-M[0][0][0]-M[1][0][0]+M[2][0][0]+M[3][0][0])/16,  // Cx
				(-M[0][1][0]+M[1][1][0]-M[2][1][0]+M[3][1][0])/8,  // Ay
				(-M[0][1][1]-M[1][1][1]+M[2][1][1]+M[3][1][1])/8,  // By
				(-M[0][1][1]+M[1][1][1]-M[2][1][1]+M[3][1][1]-M[0][1][0]-M[1][1][0]+M[2][1][0]+M[3][1][0])/16,   // Cy
				( M[0][0][0]+M[1][0][0]+M[2][0][0]+M[3][0][0])/4-1, // Dx
				( M[0][0][1]+M[1][0][1]+M[2][0][1]+M[3][0][1])/4,   // Ex
				( M[0][1][0]+M[1][1][0]+M[2][1][0]+M[3][1][0])/4,   // Dy
				( M[0][1][1]+M[1][1][1]+M[2][1][1]+M[3][1][1])/4-1  // Ey
		};

		if (this.debugLevel>2) { /* increase LEVEL later */
			/*          System.out.println("Center   "+
	                           " W0x="+     IJ.d2s(centerWV[0][0],4)+
	                           " W0y="+     IJ.d2s(centerWV[0][1],4)+
	                           " W1x="+     IJ.d2s(centerWV[1][0],4)+
	                           " W1y="+     IJ.d2s(centerWV[1][1],4));*/
			for (iq=0; iq<4;iq++) {
				System.out.println("Matrix= "+iq+
						" M00="+     IJ.d2s(mCorners[iq][0][0],4)+
						" M01="+     IJ.d2s(mCorners[iq][0][1],4)+
						" M10="+     IJ.d2s(mCorners[iq][1][0],4)+
						" M11="+     IJ.d2s(mCorners[iq][1][1],4));
			}
			for (iq=0; iq<4;iq++) {
				System.out.println(" M["+iq+"][0][0]="+  IJ.d2s(M[iq][0][0],4)+
						" M["+iq+"][0][1]="+  IJ.d2s(M[iq][0][1],4)+
						" M["+iq+"][1][0]="+  IJ.d2s(M[iq][1][0],4)+
						" M["+iq+"][1][1]="+  IJ.d2s(M[iq][1][1],4));
			}
		}
		if (this.debugLevel>2) { /* increase LEVEL later */
			System.out.println("Corr   "+
					" Ax="+     IJ.d2s(rslt[0],5)+
					" Bx="+     IJ.d2s(rslt[1],5)+
					" Cx="+     IJ.d2s(rslt[2],5)+
					" Ay="+     IJ.d2s(rslt[3],5)+
					" By="+     IJ.d2s(rslt[4],5)+
					" Cy="+     IJ.d2s(rslt[5],5)+
					" Dx="+     IJ.d2s(rslt[6],5)+
					" Ex="+     IJ.d2s(rslt[7],5)+
					" Dy="+     IJ.d2s(rslt[8],5)+
					" Ey="+     IJ.d2s(rslt[9],5));
		}
		return rslt;
	};
	/* ======================================================================== */
	/* TODO: REPLACE doubleFHT  */
	/* converts FHT results (frequency space) to complex numbers of [fftsize/2+1][fftsize] */
	private double[][][] FHT2FFTHalf (FHT fht, int fftsize) {
		float[] fht_pixels=(float[])fht.getPixels();
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
	/* ======================================================================== */
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
	/* ======================================================================== */
	/* converts FFT arrays of complex numbers of [fftsize/2+1][fftsize] to FHT arrays */
	private float[] floatFFTHalf2FHT (double [][][] fft, int fftsize) {
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
	/* ======================================================================== */
	public	 double[][] normalizeAndWindow (double [][] pixels, double [] windowFunction) {
		return normalizeAndWindow (pixels, windowFunction, true);
	}
	private	 double[] normalizeAndWindow (double [] pixels, double [] windowFunction) {
		return normalizeAndWindow (pixels, windowFunction, true);
	}
	private	 double[][] normalizeAndWindow (double [][] pixels, double [] windowFunction, boolean removeDC) {
		int i;
		for (i=0;i<pixels.length;i++)  if (pixels[i]!=null) pixels[i]=normalizeAndWindow (pixels[i],  windowFunction, removeDC);
		return pixels;
	}
	private	 double[] normalizeAndWindow (double [] pixels, double [] windowFunction, boolean removeDC) {
		int j;
		if (pixels==null) return null;
		double s=0.0,s0=0.0;
		if (removeDC) {
			for (j=0;j<pixels.length;j++){
				s+=pixels[j]*windowFunction[j];
				s0+=windowFunction[j];

			}
			s/=s0;
		}
		for (j=0;j<pixels.length;j++) pixels[j]=(pixels[j]-s)*windowFunction[j];
		return pixels;
	}
	/* ======================================================================== */
	public	double [][] matrix2x2_invert(double [][] m ){
		double det=m[0][0]*m[1][1]-m[0][1]*m[1][0];
		double [][] rslt= {{ m[1][1]/det,  -m[0][1]/det},
				{-m[1][0]/det,   m[0][0]/det}};
		return rslt;
	}
	public	double [][] matrix2x2_mul(double [][] a, double [][] b ){
		double [][] rslt={{a[0][0]*b[0][0]+a[0][1]*b[1][0], a[0][0]*b[0][1]+a[0][1]*b[1][1]},
				{a[1][0]*b[0][0]+a[1][1]*b[1][0], a[1][0]*b[0][1]+a[1][1]*b[1][1]}};
		return rslt;
	}
	public	double []   matrix2x2_mul(double [][] a, double [] b ){
		double [] rslt={a[0][0]*b[0]+a[0][1]*b[1],
				a[1][0]*b[0]+a[1][1]*b[1]};
		return rslt;
	}
	public	double [][] matrix2x2_scale(double [][] a, double  b ){
		double [][] rslt={{a[0][0]*b, a[0][1]*b},
				{a[1][0]*b, a[1][1]*b}};
		return rslt;
	}
	public	double [] matrix2x2_scale(double [] a, double  b ){
		double [] rslt={a[0]*b, a[1]*b};
		return rslt;
	}
	public	double [][] matrix_add(double [][] a, double [][]  b ){
		double [][] rslt= new double [a.length][a[0].length];
		int i,j;
		for (i=0;i<rslt.length;i++) for (j=0;j<rslt[0].length;j++) rslt[i][j]=a[i][j]+b[i][j];
		return rslt;
	}
	public	double []   vector_add(double [] a, double [] b ){
		double [] rslt= new double [a.length];
		int i;
		for (i=0;i<rslt.length;i++) rslt[i]=a[i]+b[i];
		return rslt;
	}
	/* calculates 2x2 matrix that converts two pairs of vectors: u2=M*u1, v2=M*v1*/
	private	double [][] matrixToConvertTwoPairs(double [] u1, double [] v1, double [] u2, double [] v2) {
		double [][] rslt= {{(u2[0]*v1[1]-v2[0]*u1[1])/(u1[0]*v1[1]-v1[0]*u1[1]),
			(v2[0]*u1[0]-u2[0]*v1[0])/(u1[0]*v1[1]-v1[0]*u1[1])},
				{(u2[1]*v1[1]-v2[1]*u1[1])/(u1[0]*v1[1]-v1[0]*u1[1]),
				(v2[1]*u1[0]-u2[1]*v1[0])/(u1[0]*v1[1]-v1[0]*u1[1])}};
		return rslt;
	}

	public double [][] matrix2x2_add(double [][] a, double [][] b ){
		double [][] rslt={{a[0][0]+b[0][0], a[0][1]+b[0][1]},
				{a[1][0]+b[1][0], a[1][1]+b[1][1]}};
		return rslt;
	}

	public double [] matrix2x2_add(double [] a, double [] b ){
		double [] rslt={a[0]+b[0], a[1]+b[1]};
		return rslt;
	}

	public int [][] matrix2x2_add(int [][] a, int [][] b ){
		int [][] rslt={{a[0][0]+b[0][0], a[0][1]+b[0][1]},
				{a[1][0]+b[1][0], a[1][1]+b[1][1]}};
		return rslt;
	}

	public int [] matrix2x2_add(int [] a, int [] b ){
		int [] rslt={a[0]+b[0], a[1]+b[1]};
		return rslt;
	}

	public double [][] matrix2x2_transp(double [][] m ){
		double [][] rslt= {{ m[0][0],  m[1][0]},
				{ m[0][1],  m[1][1]}};
		return rslt;
	}

	/* ======================================================================== */
	private	 double [] combineDiagonalGreens (double [] green0, double []green3, int half_width, int half_height) {
		int y,x,base;
		int base_b=0;
		double [] result= new double [green0.length];
		for (y=0;y<half_height/2; y++){
			base=half_height*half_width/2+ y* (half_width+1);
			for (x=0; x<half_width/2; x++) {
				result[base_b++]=green0[base];
				base-=half_width;
				result[base_b++]=green3[base++];
			}
			base=half_height*half_width/2+ y* (half_width+1);
			for (x=0; x<half_width/2; x++) {
				//System.out.println("2:y="+y+" x="+x+" base_b="+base_b+" base="+base);
				result[base_b++]=green3[base++];
				result[base_b++]=green0[base];
				base-=half_width;
			}
		}
		return result;
	}
	/* ======================================================================== */
	public double[] initWindowFunction(int size, double gaussWidth) {
		return initWindowFunction(size, gaussWidth, 0);
	}
	public double[] initWindowFunction(int size, double gaussWidth, int zeros) {
		double [] windowFunction =new double [size*size];
		double [] windowFunction_line=new double [size];
		double a,k;
		int i,j;
		int size1=size-zeros;
		int i0=(zeros+1)/2;
		if (gaussWidth==0) {
			for (i=0; i<size; i++) windowFunction_line[i]= 1.0;
		} else if (gaussWidth<0) {
			for (i=0; i<size1; i++) windowFunction_line[i+i0]= (0.54-0.46*Math.cos((i*2.0*Math.PI)/size1));
		} else {
			k=2.0/(size*gaussWidth);
			for (i=i0; i<i0+size1; i++) {
				a=(i-size/2)*k;
				windowFunction_line[i]= Math.exp( - a*a);
				if (zeros>0) windowFunction_line[i]*=(0.54-0.46*Math.cos(((i-i0)*2.0*Math.PI)/size1)); // additionally multiply by Hamming
			}
		}
		if (zeros>0){ // make window to be exact zero for certain number of samples (for correlation)
			for (i=0;i<i0;i++) windowFunction_line[i]=0.0;
			for (i=i0+size1;i<size;i++) windowFunction_line[i]=0.0;


		}
		for (i=0; i<size; i++) for (j=0; j<size; j++){
			windowFunction[size*i+j]=windowFunction_line[i]*windowFunction_line[j];
		}
		return windowFunction;
	}
	/* ============================= Distortions ===================================*/
	public void distortionsTest (
			final DistortionParameters distortionParameters, //
			final MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			final	SimulationPattern.SimulParameters  simulParameters,
			final boolean equalizeGreens,
			final ImagePlus imp, // image to process
			final int threadsMax,
			final boolean updateStatus,
			final int debug_level){// debug level used inside loops

		if (imp==null) return;

		Roi roi= imp.getRoi();
		final Rectangle selection;
		if (roi==null){
			selection=new Rectangle(0, 0, imp.getWidth(), imp.getHeight());
		} else {
			selection=roi.getBounds();
		}
		MatchSimulatedPattern matchSimulatedPattern=new MatchSimulatedPattern(distortionParameters.FFTSize);
		matchSimulatedPattern.debugLevel=debugLevel;
		MatchSimulatedPattern matchSimulatedPatternCorr=new MatchSimulatedPattern(distortionParameters.correlationSize);
		matchSimulatedPatternCorr.debugLevel=debugLevel;
		final SimulationPattern.SimulParameters  thisSimulParameters=simulParameters.clone();
		thisSimulParameters.subdiv=     distortionParameters.patternSubdiv;
		thisSimulParameters.bPatternSigma=distortionParameters.bPatternSigma;
		thisSimulParameters.barraySigma=distortionParameters.barraySigma;
		SimulationPattern simulationPattern= new SimulationPattern(thisSimulParameters);
		final double [] bPattern= simulationPattern.patternGenerator(simulParameters); // reuse pattern for next time
		//find center of the selection (to be used to find initial pattern approximation)
		if (debugLevel>2)	System.out.println("bPattern.length="+bPattern.length);
		int xc,yc;
		xc=2*((2*selection.x+selection.width+1)/4);
		yc=2*((2*selection.y+selection.height+1)/4);
		Rectangle initialPatternCell=new Rectangle(xc-distortionParameters.FFTSize,
				yc-distortionParameters.FFTSize,
				2*distortionParameters.FFTSize,2*distortionParameters.FFTSize);
		//create diagonal green selection around xc,yc
		double [][] input_bayer=splitBayer (imp,initialPatternCell,equalizeGreens);
		if (debugLevel>2) SDFA_INSTANCE.showArrays(input_bayer,  true, "selection-bayer-distortionsTest");
		double [] windowFunction=initWindowFunction(distortionParameters.FFTSize, distortionParameters.fftGaussWidth);
		final double [] windowFunctionCorr=initWindowFunction(distortionParameters.correlationSize,distortionParameters.correlationGaussWidth,distortionParameters.zeros);
		double [] greens=normalizeAndWindow (input_bayer[4], windowFunction);

		double [][] pattern=matchSimulatedPattern.findPattern(
				null, // 			DoubleFHT doubleFHT,
				greens,
				distortionParameters.FFTSize,
				patternDetectParameters,
				patternDetectParameters.minGridPeriod/2,
				patternDetectParameters.maxGridPeriod/2,
				true, // this is a pattern for combined greens (diagonal), adjust results accordingly
				"Pattern"); // title - will not be used
		if (pattern==null) {
			System.out.println("Error - pattern not found");
			IJ.showMessage("Error","Failed to find pattern");
			return;
		}
		if (debugLevel>2) System.out.println("FX1="+pattern[0][0]+"  FY1="+pattern[0][1]+"  phase1="+pattern[0][2]);
		if (debugLevel>2) System.out.println("FX2="+pattern[1][0]+"  FY2="+pattern[1][1]+"  phase2="+pattern[1][2]);
		double [] ll2=new double[2];
		double [][] dxy=new double[2][2];
		double [][] phases=new double[2][2];
		int i,j,k;
		for (i=0;i<2;i++) {
			ll2[i]=(pattern[i][0]*pattern[i][0]+pattern[i][1]*pattern[i][1])*2*Math.PI;
		}

		if (debugLevel>2) System.out.println("phase1/2pi="+(pattern[0][2]/2/Math.PI)+"  phase2/2pi="+(pattern[1][2]/2/Math.PI));
		for (k=0;k<2;k++) {
			for (j=0;j<2;j++) {
				phases[k][j]= pattern[j][2]-Math.PI/2+((k>0)? Math.PI:0.0);
				while (phases[k][j]<-Math.PI) phases[k][j]+=2*Math.PI;
				while (phases[k][j]>Math.PI)  phases[k][j]-=2*Math.PI;
			}
			if (debugLevel>2) System.out.println("phase1/2pi="+(phases[k][0]/2/Math.PI)+"  phase2/2pi="+(phases[k][1]/2/Math.PI));
			for (i=0;i<2;i++) {
				dxy[k][i]=0;
				for (j=0;j<2;j++) {
					dxy[k][i]+=(phases[k][j])*pattern[j][i]/ll2[j];
				}
			}
			if (debugLevel>2) System.out.println("dX["+k+"]="+dxy[k][0]+" dY["+k+"]="+dxy[k][1]);
		}
		int phaseSel=((dxy[0][0]*dxy[0][0]+dxy[0][1]*dxy[0][1])<(dxy[1][0]*dxy[1][0]+dxy[1][1]*dxy[1][1]))?0:1;
		if (debugLevel>1) System.out.println("xc="+xc+" yc="+yc);
		if (debugLevel>1) System.out.println("dX="+dxy[phaseSel][0]+" dY="+dxy[phaseSel][1]);
		double [] centerXY0={xc-dxy[phaseSel][0],yc-dxy[phaseSel][1]};
		if (debugLevel>1) System.out.println("+++ Initial center x="+IJ.d2s(centerXY0[0],3)+" y="+ 	IJ.d2s(centerXY0[1],3));

		// debug mode - scan correlation around center point, show result and exit:
		SDFA_INSTANCE.showArrays(simulationPattern.bPattern, "bPattern");
		//			 double [] barray= new double [simulationPattern.barray.length*simulationPattern.barray[0].length];
		//			 for (i=0;i<barray.length;i++) {
		//				 barray[i]=simulationPattern.barray[i/simulationPattern.barray[0].length][i % simulationPattern.barray[0].length];
		//			 }
		//			 SDFA_INSTANCE.showArrays(barray, simulationPattern.barray[0].length, simulationPattern.barray.length,"barray");
		SDFA_INSTANCE.showArrays(simulationPattern.barray, "barray");

		double [][][] scanXY=scanPatternCrossLocation(
				distortionParameters.correlationDx, // range
				(int) Math.round(distortionParameters.correlationDx/distortionParameters.correlationDy)+1,
				centerXY0, // initial coordinates of the pattern cross point
				pattern[0][0],
				pattern[0][1],
				pattern[1][0],
				pattern[1][1],
				imp,                       // image data (Bayer mosaic)
				distortionParameters,      //
				patternDetectParameters,
				matchSimulatedPatternCorr, // correlationSize
				thisSimulParameters,
				equalizeGreens,
				windowFunctionCorr,        // window function
				simulationPattern,
				false, // if true - invert pattern
				null); // will create new instance of DoubleFHT class
		double [][] scanImg= new double [4][scanXY.length*scanXY[0].length];
		//			  System.out.println("scanImg[0].length="+scanImg[0].length);
		for (i=0;i<scanImg[0].length;i++) {
			scanImg[0][i]=scanXY[i/scanXY[0].length][i % scanXY[0].length][0];
			scanImg[1][i]=scanXY[i/scanXY[0].length][i % scanXY[0].length][1];
			scanImg[2][i]=scanXY[i/scanXY[0].length][i % scanXY[0].length][2];
			scanImg[3][i]=scanXY[i/scanXY[0].length][i % scanXY[0].length][3];
		}
		SDFA_INSTANCE.showArrays(scanImg, true, "scan_correlation");
		SDFA_INSTANCE.showArrays(scanImg, false,"scan_correlation");
		return;
	}
	/*
	 *  Try point x,y, test for pattern, return x,y, contrast (or null)
	 *   [0][0] - x
	 *   [0][1] - y
	 *   [0][2] - contrast
	 *   [1][0] - Wave vector 1 x component
	 *   [1][1] - Wave vector 1 y component
	 *   [1][2] - Wave vector 1 phase (not used here)
	 *   [2][0] - Wave vector 2 x component
	 *   [2][1] - Wave vector 2 y component
	 *   [2][2] - Wave vector 2 phase (not used here)
	 */
	public double[][] tryPattern (
			LwirReaderParameters lwirReaderParameters, // null is OK
			DoubleFHT doubleFHT,
			double [] point, // xy to try
			final DistortionParameters distortionParameters, //
			final MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			final double    min_half_period,
			final double    max_half_period,
			final SimulationPattern.SimulParameters  thisSimulParameters,
			final MatchSimulatedPattern matchSimulatedPattern,
			final MatchSimulatedPattern matchSimulatedPatternCorr,
			final SimulationPattern simulationPattern,
			final boolean equalizeGreens,
			final ImagePlus imp, // image to process
			double [] bPattern,
			double [] windowFunction,
			double [] windowFunctionCorr,
			double [] windowFunctionCorr2,
			double [] windowFunctionCorr4,
			double[][] locsNeib, // which neibors to try (here - just the center)
			String dbgStr
			){

//		this.debugLevel = 3;
		if (this.debugLevel == 3) {
			System.out.println("tryPattern(): this.debugLevel = 3");
		}
		int debug_threshold = 2;
		if (imp==null) {
			if (dbgStr!=null) System.out.println(dbgStr+" imp==null");
			return null;
		}
		boolean is_lwir = ((lwirReaderParameters != null) && lwirReaderParameters.is_LWIR(imp));
		int     fft_size = is_lwir ? distortionParameters.FFTSize_lwir : distortionParameters.FFTSize;

		int xc= (int)(2*Math.round(0.5*point[0]));
		int yc= (int)(2*Math.round(0.5*point[1]));
		Roi roi= imp.getRoi();
		final Rectangle selection;
		if (roi==null){
			selection=new Rectangle(0, 0, imp.getWidth(), imp.getHeight());
		} else {
			selection=roi.getBounds();
		}
		Rectangle initialPatternCell=new Rectangle(xc-fft_size,
				yc-fft_size,
				2*fft_size,2*fft_size);
		if (is_lwir) {// move to to all? It requires twice FFT size to be able to select, maybe single is still OK?
			initialPatternCell=new Rectangle(xc-fft_size/2,
					yc-fft_size/2,
					fft_size,fft_size);

		}
		if (!selection.contains(initialPatternCell)) {
			if (dbgStr!=null) System.out.println(dbgStr+" selection ("+
					selection.x+","+selection.y+","+selection.width+","+selection.height+ ") does not contain cell ("+
					initialPatternCell.x+","+initialPatternCell.y+","+initialPatternCell.width+","+initialPatternCell.height+ ")");
			return null; // area for FFT is not inside the initial selection
		}
		double [][] pattern = null;
		if (is_lwir) { // monochrome
			double [] dpixels =  getNoBayer(imp, initialPatternCell);
			double [] dbg_dpixels = dpixels.clone();
			normalizeAndWindow (dpixels, windowFunction);
			if (debugLevel > (debug_threshold + 0)) {
				double [][] dbg_img = {dbg_dpixels, dpixels};
				SDFA_INSTANCE.showArrays(dbg_img,  true, "selection-input"+
				(initialPatternCell.x+initialPatternCell.width/2)+":"+(initialPatternCell.y+initialPatternCell.height/2));
			}
			pattern=matchSimulatedPattern.findPattern(
					doubleFHT,
					dpixels,
					fft_size,
					patternDetectParameters,
					min_half_period,
					max_half_period,
					false, // this is a pattern for combined greens (diagonal), adjust results accordingly
					"Pattern"); // title - will not be used
		} else { // Bayer - extract green
			//create diagonal green selection around xc,yc
			double [][] input_bayer=splitBayer (imp,initialPatternCell,equalizeGreens);
			if (debugLevel > (debug_threshold + 0)) {
				SDFA_INSTANCE.showArrays(input_bayer,  true, "selection--bayer");
			}
			double [] greens=normalizeAndWindow (input_bayer[4], windowFunction);
			pattern=matchSimulatedPattern.findPattern(
					doubleFHT,
					greens,
					fft_size,
					patternDetectParameters,
					min_half_period,
					max_half_period,
					true, // this is a pattern for combined greens (diagonal), adjust results accordingly
					"Pattern"); // title - will not be used
		}
		if (pattern==null) {
			//				System.out.println("Error - pattern not found");
			//			  	IJ.showMessage("Error","Failed to find pattern");
			if (dbgStr!=null) System.out.println(dbgStr+" matchSimulatedPattern.findPattern->null");
			return null;
		}
		if (debugLevel > (debug_threshold + 0)) System.out.println("FX1="+pattern[0][0]+"  FY1="+pattern[0][1]+"  phase1="+pattern[0][2]);
		if (debugLevel > (debug_threshold + 0)) System.out.println("FX2="+pattern[1][0]+"  FY2="+pattern[1][1]+"  phase2="+pattern[1][2]);
		double [] ll2=new double[2];
		double [][] dxy=new double[2][2];
		double [][] phases=new double[2][2];
		int i,j,k;
		for (i=0;i<2;i++) {
			ll2[i]=(pattern[i][0]*pattern[i][0]+pattern[i][1]*pattern[i][1])*2*Math.PI;
		}

		if (debugLevel > (debug_threshold + 0)) System.out.println("phase1/2pi="+(pattern[0][2]/2/Math.PI)+"  phase2/2pi="+(pattern[1][2]/2/Math.PI));
		for (k=0;k<2;k++) {
			for (j=0;j<2;j++) {
				phases[k][j]= pattern[j][2]-Math.PI/2+((k>0)? Math.PI:0.0);
				while (phases[k][j]<-Math.PI) phases[k][j]+=2*Math.PI;
				while (phases[k][j]>Math.PI)  phases[k][j]-=2*Math.PI;
			}
			if (debugLevel > (debug_threshold + 0)) System.out.println("phase1/2pi="+(phases[k][0]/2/Math.PI)+"  phase2/2pi="+(phases[k][1]/2/Math.PI));
			for (i=0;i<2;i++) {
				dxy[k][i]=0;
				for (j=0;j<2;j++) {
					dxy[k][i]+=(phases[k][j])*pattern[j][i]/ll2[j];
				}
			}
			if (debugLevel > (debug_threshold + 0)) System.out.println("dX["+k+"]="+dxy[k][0]+" dY["+k+"]="+dxy[k][1]);
		}
		int phaseSel=((dxy[0][0]*dxy[0][0]+dxy[0][1]*dxy[0][1])<(dxy[1][0]*dxy[1][0]+dxy[1][1]*dxy[1][1]))?0:1;
		if (debugLevel > (debug_threshold + 0)) System.out.println("xc="+xc+" yc="+yc);
		if (debugLevel > (debug_threshold + 0)) System.out.println("dX="+dxy[phaseSel][0]+" dY="+dxy[phaseSel][1]);
		double [] centerXY0={xc-dxy[phaseSel][0],yc-dxy[phaseSel][1]};
		if (debugLevel > (debug_threshold - 1)) System.out.println("+++ Initial center x="+IJ.d2s(centerXY0[0],3)+" y="+ 	IJ.d2s(centerXY0[1],3));
/*
		if (debugLevel > (debug_threshold + 0)){
			double [] test_error_offset = {0.0, 2.0};
			centerXY0[0] += test_error_offset[0];
			centerXY0[1] += test_error_offset[1];
			System.out.println("+++ Initial center with intentional error offset: x="+IJ.d2s(centerXY0[0],3)+" y="+ 	IJ.d2s(centerXY0[1],3));
		}
*/
		double [] centerXY=correctedPatternCrossLocation(
				lwirReaderParameters, // LwirReaderParameters lwirReaderParameters, // null is OK
				centerXY0, // initial coordinates of the pattern cross point
				pattern[0][0],
				pattern[0][1],
				pattern[1][0],
				pattern[1][1],
				null, // correction
				imp,      // image data (Bayer mosaic)
				distortionParameters, //
				patternDetectParameters,
				matchSimulatedPatternCorr, // correlationSize
				thisSimulParameters,
				equalizeGreens,
				windowFunctionCorr,   // window function
				windowFunctionCorr2,   // window function
				windowFunctionCorr4,   // window function
				simulationPattern,
				false, // is_lwir, // false, // if true - invert pattern TODO: Put is_lwir here?
				null, // will create new instance of DoubleFHT class
				distortionParameters.fastCorrelationOnFirstPass,
				locsNeib,
				debugLevel,
				dbgStr);
		if (debugLevel > (debug_threshold - 1)) System.out.println("--- Initial center x="+IJ.d2s(centerXY0[0],3)+" y="+ 	IJ.d2s(centerXY0[1],3)+
				" -> "+((centerXY==null)?" NULL ":(IJ.d2s(centerXY[0],3)+" : "+ 	IJ.d2s(centerXY[1],3))));
		// TODO: check statistics of correction - see if there is any dc bias

/*
		if ((debugLevel > (debug_threshold + 0)) && (centerXY != null)) {
			// trying second time - watch for conversion
			double [] centerXY1=correctedPatternCrossLocation(
					lwirReaderParameters, // LwirReaderParameters lwirReaderParameters, // null is OK
					centerXY, // initial coordinates of the pattern cross point
					pattern[0][0],
					pattern[0][1],
					pattern[1][0],
					pattern[1][1],
					null, // correction
					imp,      // image data (Bayer mosaic)
					distortionParameters, //
					patternDetectParameters,
					matchSimulatedPatternCorr, // correlationSize
					thisSimulParameters,
					equalizeGreens,
					windowFunctionCorr,   // window function
					windowFunctionCorr2,   // window function
					windowFunctionCorr4,   // window function
					simulationPattern,
					false, // is_lwir, // false, // if true - invert pattern TODO: Put is_lwir here?
					null, // will create new instance of DoubleFHT class
					distortionParameters.fastCorrelationOnFirstPass,
					locsNeib,
					debugLevel,
					dbgStr);
			System.out.println("--- Second run: initial center x="+IJ.d2s(centerXY[0],3)+" y="+ 	IJ.d2s(centerXY[1],3)+
					" -> "+((centerXY1==null)?" NULL ":(IJ.d2s(centerXY1[0],3)+" : "+ 	IJ.d2s(centerXY1[1],3))));



		}
*/




		double [][] node = {centerXY,pattern[0],pattern[1]};
		if (dbgStr!=null) {
			if (centerXY==null) System.out.println(dbgStr+" correctedPatternCrossLocation->null");
			else System.out.println(dbgStr+" matchSimulatedPattern.findPattern->{{"+centerXY[0]+","+centerXY[1]+","+centerXY[2]+"}..}");
		}
		return node;
	}

	// ============= end of public double[][] tryPattern() ==================

	/* ================================================================*/
	// Optionally remove the outer (possibly corrupted) layer of the detected pattern nodes, extrapolate new layers of the nodes
	// without pattern matching

	public double[][][][] finalizeDistortionsBorder (
			//			   final double [][][][] patternGrid,
			final DistortionParameters distortionParameters, //
			final boolean updateStatus,
			final int debug_level){// debug level used inside loops
		//		    double[][][][] patternGrid=this.PATTERN_GRID;

		//	        final int [][] directionsUV=  {{1,0},{0,1},{-1,0},{0,-1}}; // should have opposite direction shifted by half
		final int [][] directionsUV8= {{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}}; // first 8 should be the same as in directionsUV

		final List <Integer> waveFrontList=new ArrayList<Integer>(1000);
		// create list of all nodes that have undefined neigbors (up/down/right/left)
		int umax=0, vmax=0, vmin=this.PATTERN_GRID.length, umin=this.PATTERN_GRID[0].length;
		for (int i=0;i<this.PATTERN_GRID.length;i++) for (int j=0;j<this.PATTERN_GRID[i].length;j++) {
			if ((this.PATTERN_GRID[i][j]!=null) && (this.PATTERN_GRID[i][j][0]!=null)) {
				if (vmin > i) vmin = i;
				if (vmax < i) vmax = i;
				if (umin > j) umin = j;
				if (umax < j) umax = j;
			}
		}
		int [] uvNew=new int [2];
		int [] iUV=  new int [2];
		int [] uvdir; // (u,v,direction}
		double [][][] wave;
		for (uvNew[1]=vmin;uvNew[1]<=vmax;uvNew[1]++) for (uvNew[0]=umin;uvNew[0]<=umax;uvNew[0]++) if (isCellDefined(this.PATTERN_GRID,uvNew)){
			for (int dir=0;dir<directionsUV8.length;dir++) {
				iUV[0]=uvNew[0]+directionsUV8[dir][0];
				iUV[1]=uvNew[1]+directionsUV8[dir][1];
				if (!isCellDefined(this.PATTERN_GRID,iUV)){
					putInWaveList (waveFrontList, uvNew, dir); // direction does not matter here
					break;
				}
			}

		}
		final double [][] extrapolationWeights=generateWeights (
				distortionParameters.extrapolationSigma,
				distortionParameters.correlationRadiusScale); //  if 0 - use sigma as radius, inside - 1.0, outside 0.0. If >0 - size of array n*sigma

		if (debugLevel>1) System.out.println("***** finalizeDistortionsBorder, initial wave length="+waveFrontList.size());
		// optionally remove outer (possibly corruopted) layer of nodes
		if (distortionParameters.removeLast) for (int i=0;i<waveFrontList.size();i++){
			if (distortionParameters.numberExtrapolated==0)	invalidatePatternGridCell(this.PATTERN_GRID, getWaveList (waveFrontList,i));
			else 	                                              initPatternGridCell(this.PATTERN_GRID, getWaveList (waveFrontList,i));
		}
		for (int layer=0;(layer<distortionParameters.numberExtrapolated) && (waveFrontList.size()>0);layer++){

			if ((layer>0) || !distortionParameters.removeLast) { // build new layer around the current one
				while (waveFrontList.size()>0) { // will normally break out of the cycle
					uvdir= getWaveList (waveFrontList,0);
					//							if (this.PATTERN_GRID[uvdir[1]][uvdir[0]]==null) break; // finished adding new layer
					if (!isCellDefined(this.PATTERN_GRID,uvdir)) break; // finished adding new layer, hit one of the newely added
					for (int dir=0;dir<directionsUV8.length;dir++) {
						iUV[0]=uvdir[0]+directionsUV8[dir][0];
						iUV[1]=uvdir[1]+directionsUV8[dir][1];
						if ((iUV[0]<0) || (iUV[1]<0) ||
								(iUV[0]>=distortionParameters.gridSize) || (iUV[1]>=distortionParameters.gridSize)) continue; // don't fit into UV grid
						if (!isCellNew(this.PATTERN_GRID,iUV)) continue; // already processed
						// add uv and dir to the list
						putInWaveList (waveFrontList, iUV, dir); //  direction is not used
						initPatternGridCell(this.PATTERN_GRID, iUV);
						//								if (debugLevel>1) System.out.println("-->iUV= "+iUV[0]+",  "+iUV[1]+",  "+((dir+(directionsUV.length/2))%directionsUV.length));
					}
					waveFrontList.remove(0);  // remove first element from the list
					//		    			if (debugLevel>1) System.out.println("xx> remove(0), (waveFrontList.size()="+(waveFrontList.size()));
				}
			}
			// extrapolate x,y for the new layer of pixels (not yet using the new pixels in this layer)
			if (updateStatus) IJ.showStatus("Extrapolating border, layer "+(layer+1)+", length "+waveFrontList.size());
			if (debugLevel>1) System.out.println("Extrapolating border, layer "+(layer+1)+", length "+waveFrontList.size());
			wave = new double [waveFrontList.size()][][];
			for (int i=0;i<wave.length;i++) {
				wave[i]=estimateCell(
						this.PATTERN_GRID,
						getWaveList (waveFrontList,i),
						extrapolationWeights, // quadrant of sample weights
						true, // useContrast
						!distortionParameters.useQuadratic,  // use linear approximation (instead of quadratic)
						1.0E-10,  // thershold ratio of matrix determinant to norm for linear approximation (det too low - fail)
						1.0E-20  // thershold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
						);
				if (wave[i]==null) { // try w/o contrast, just x,y
					wave[i]=estimateCell(
							this.PATTERN_GRID,
							getWaveList (waveFrontList,i),
							extrapolationWeights, // quadrant of sample weights
							false, // do not use Contrast, keep old contrast (even if it is NaN)
							!distortionParameters.useQuadratic,  // use linear approximation (instead of quadratic)
							1.0E-10,  // thershold ratio of matrix determinant to norm for linear approximation (det too low - fail)
							1.0E-20  // thershold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
							);

				}
			}
			// set new values, removed failed cells (normally should not be any)
			for (int i=wave.length-1;i>=0;i--) {
				uvdir=getWaveList (waveFrontList,i);
				this.PATTERN_GRID[uvdir[1]][uvdir[0]]=wave[i]; // null OK
				if (wave[i]==null) {
					if (debugLevel>0) System.out.println("Removing failed node (normally should not happen!), u="+uvdir[0]+", v="+uvdir[1]);
					waveFrontList.remove(i);
				}

			}

		}

		return null;
	}
	/* ================================================================*/
	// it now can start with non-empty Grid
	public int distortions( // returns number of grid cells
			final LwirReaderParameters lwirReaderParameters, // null is OK
			final boolean [] triedIndices,
			final DistortionParameters distortionParameters, //
			final MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			final double    min_half_period,
			final double    max_half_period,
			final SimulationPattern.SimulParameters  simulParameters,
			final boolean equalizeGreens,
			final ImagePlus imp, // image to process
			final int threadsMax,
			final boolean updateStatus,
			final int debug_level, // debug level used inside loops
			final int global_debug_level){
		if (imp==null) return 0;
		final int debugThreshold = 1;
		final boolean is_lwir = ((lwirReaderParameters != null) && lwirReaderParameters.is_LWIR(imp));
		final int      fft_size =             is_lwir ? distortionParameters.FFTSize_lwir :               distortionParameters.FFTSize;
		final int      correlation_size =     is_lwir ? distortionParameters.correlationSizeLwir :        distortionParameters.correlationSize;
		final int      max_correlation_size = is_lwir ? distortionParameters.maximalCorrelationSizeLwir : distortionParameters.maximalCorrelationSize;
		final int      minimal_pattern_cluster = is_lwir ? distortionParameters.minimalPatternClusterLwir : distortionParameters.minimalPatternCluster;
		final int [][] directionsUV=  {{1,0},{0,1},{-1,0},{0,-1}}; // should have opposite direction shifted by half
		final int [][] directionsUV8= {{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}}; // first 8 should be the same as in directionsUV
		final int [] directionsBits8= {1,4,1,4,2,8,2,8}; // should match directionsUV8
		int neibBits;
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger cellNum = new AtomicInteger(0);
		final List <Integer> waveFrontList=new ArrayList<Integer>(1000);
		final int [] centerUV= {distortionParameters.gridSize/2, distortionParameters.gridSize/2};

		final double [][] locsNeib=calcNeibLocsWeights (distortionParameters,false); // no neighbors to average
		Roi roi= imp.getRoi();
		final Rectangle selection;
		if (PATTERN_GRID==null) {
			if (roi==null){
				selection=new Rectangle(0, 0, imp.getWidth(), imp.getHeight());
			} else {
				selection=roi.getBounds();
			}
		} else {
			if ((getImageHeight()!=imp.getHeight()) || (getImageWidth()!=imp.getWidth())){
				String msg="Supplied image does not match in dimensions "+imp.getWidth()+"x"+imp.getHeight()+
						" the one for wich grid was calculated ("+getImageWidth()+"x"+getImageHeight()+")";
				IJ.showMessage("Error",msg);
				throw new IllegalArgumentException (msg);
			}
			selection=new Rectangle(0, 0, getImageWidth(), getImageHeight());
		}
		MatchSimulatedPattern matchSimulatedPattern=new MatchSimulatedPattern(fft_size);
		matchSimulatedPattern.debugLevel=debugLevel;
		MatchSimulatedPattern matchSimulatedPatternCorr=new MatchSimulatedPattern(correlation_size);
		matchSimulatedPatternCorr.debugLevel=debugLevel;
		final SimulationPattern.SimulParameters  thisSimulParameters=simulParameters.clone();
		thisSimulParameters.subdiv=     distortionParameters.patternSubdiv;
		thisSimulParameters.bPatternSigma=distortionParameters.bPatternSigma;
		thisSimulParameters.barraySigma=distortionParameters.barraySigma;
		SimulationPattern simulationPattern= new SimulationPattern(thisSimulParameters);
		final double [] bPattern= simulationPattern.patternGenerator(simulParameters); // reuse pattern for next time
		//find center of the selection (to be used to find initial pattern approximation)
		if (debugLevel>2)	System.out.println("bPattern.length="+bPattern.length);
		double [] windowFunction=initWindowFunction(fft_size, distortionParameters.fftGaussWidth);
		// may need to decrease relative Gaussian width for larger windows
		//			public boolean absoluteCorrelationGaussWidth=false; // do not scale correlationGaussWidth when the FFT size is increased
		///(distortionParameters.absoluteCorrelationGaussWidth?0.5:1.0)
		final double [] windowFunctionCorr= initWindowFunction(
				correlation_size,
				distortionParameters.correlationGaussWidth,
				distortionParameters.zeros);
		final double [] windowFunctionCorr2=initWindowFunction(
				2*correlation_size,
				(distortionParameters.absoluteCorrelationGaussWidth?0.5:1.0)*distortionParameters.correlationGaussWidth,
				distortionParameters.zeros);
		final double [] windowFunctionCorr4=initWindowFunction(
				4*correlation_size,
				(distortionParameters.absoluteCorrelationGaussWidth?0.25:1.0)*distortionParameters.correlationGaussWidth,
				distortionParameters.zeros);

		DistortionParameters thisDistortionParameters=distortionParameters.clone();
		thisDistortionParameters.correlationMaxOffset=0; // no verification of the offset here
		thisDistortionParameters.correlationMinContrast=  distortionParameters.correlationMinInitialContrast; // different contrast minimum here
		thisDistortionParameters.correlationMinAbsoluteContrast=  distortionParameters.correlationMinAbsoluteInitialContrast; // different contrast minimum here
		int was_debug_level=debugLevel;
		int [] iUV=  new int [2];
		final boolean updating=(PATTERN_GRID!=null);
		final boolean useFocusMask=updating && (focusMask!=null); // do not expand wave beyound 1 grid step from needed image regions
		Queue<GridNode> nodeQueue = new ConcurrentLinkedQueue<GridNode>();
		boolean fromVeryBeginning=true;
		if (!updating) {
			for (int i=3;i<triedIndices.length;i++) if (triedIndices[i]){ // do not count first three NPE
				fromVeryBeginning=false;
				break;
			}
			double [] point = new double[2];
			int tryHor=0,tryVert=0;
			//				distortionParameters.searchOverlap=goniometerParameters.searchOverlap;
			// with distortionParameters.searchOverlap==0.5 (default) step will be FFTSize original pixels, so half of the (2xFFTSize) square processed simultaneously
			if (distortionParameters.searchOverlap<0.1) distortionParameters.searchOverlap=0.1;
			int effectiveWidth=(int) (selection.width*0.5/distortionParameters.searchOverlap);
			int effectiveHeight=(int) (selection.height*0.5/distortionParameters.searchOverlap);

			for (int i = fft_size; i < effectiveWidth;i*=2) tryHor++;
			for (int i = fft_size; i < effectiveHeight;i*=2) tryVert++;

			int numTries=triedIndices.length-1; // Should be equal to	   int numTries=1<<(tryHor+tryVert);

			int nbv,nbh,nh,nv,nb;
			if (debugLevel>1) System.out.println("selection.x="+selection.x+" selection.y="+selection.y+" selection.width="+selection.width+" selection.height="+selection.height);
			if (debugLevel>1) System.out.println("numTries="+numTries+" tryHor="+tryHor+" tryVert="+tryVert);
			boolean oldMode=false; //true; // false;

			if (oldMode) { // old (single-threaded) mode
				for (int startScanIndex=3;startScanIndex<=numTries;startScanIndex++) if (!triedIndices[startScanIndex]){
					if (startScanIndex==numTries){
						triedIndices[startScanIndex]=true; // all done
						break;
					}
					nbh=tryHor-1;
					nbv=tryVert-1;
					nh=0;
					nv=0;
					nb=0;
					while (nb<(tryHor+tryVert)) {
						if (nbh>=0) {
							if ((startScanIndex & (1<<nb))!=0) nh |= 1<<nbh;
							nbh--;
							nb++;
						}
						if (nbv>=0) {
							if ((startScanIndex & (1<<nb))!=0) nv |= 1<<nbv;
							nbv--;
							nb++;
						}
					}
					if (debugLevel>2) System.out.println("Searching, n="+startScanIndex+", nv="+nv+", nh="+nh+", nb="+nb );
					if ((nv>0) && (nh>0)) {
						point[0]=(selection.x+nh*selection.width/(1<<tryHor)) & ~1;
						point[1]=(selection.y+nv*selection.height/(1<<tryVert)) & ~1;
						if (debugLevel>2) System.out.println("trying xc="+point[0]+", yc="+point[1]+"(nv="+nv+", nh="+nh+")");
						//						   System.out.println("### trying xc="+point[0]+", yc="+point[1]+"(nv="+nv+", nh="+nh+")");
						if ((debugLevel>2) && (startScanIndex==3)) debugLevel=3; // show debug images for the first point only
						double [][] node=tryPattern (
								lwirReaderParameters,
								null,
								point, // xy to try
								thisDistortionParameters, //no control of the displacement
								patternDetectParameters,
								min_half_period,
								max_half_period,
								thisSimulParameters,
								matchSimulatedPattern,
								matchSimulatedPatternCorr,
								simulationPattern,
								equalizeGreens,
								imp, // image to process
								bPattern,
								windowFunction,
								windowFunctionCorr,
								windowFunctionCorr2,
								windowFunctionCorr4,
								locsNeib, // which neighbors to try (here - just the center)
								null // dbgStr
								);
						debugLevel=was_debug_level;
						if ((node!=null) && (node[0]!=null)) {
							nodeQueue.add(new GridNode(node));
							break;
						}
					}
					triedIndices[startScanIndex]=true;
				}
			} else { // new multithreaded mode
				int startScanIndex=3;
				for (;(startScanIndex<numTries) && triedIndices[startScanIndex];startScanIndex++); // skip already tried indices
				if ((global_debug_level > debugThreshold) && (startScanIndex>3)) System.out.println("distortions(): startScanIndex="+startScanIndex+" > 3 ####");

				if (startScanIndex<numTries) {
					nodeQueue =  findPatternCandidates(
							lwirReaderParameters, // LwirReaderParameters lwirReaderParameters, // null is OK
							triedIndices,
							startScanIndex, // [0] will be updated
							tryHor,
							tryVert,
							//						   numTries,
							selection,
							thisDistortionParameters, //no control of the displacement
							patternDetectParameters,
							min_half_period,
							max_half_period,
							thisSimulParameters,
							matchSimulatedPattern,
							matchSimulatedPatternCorr,
							simulationPattern,
							equalizeGreens,
							imp, // image to process
							bPattern,
							windowFunction,
							windowFunctionCorr,
							windowFunctionCorr2,
							windowFunctionCorr4,
							locsNeib, // which neighbors to try (here - just the center)
							threadsMax,
							updateStatus,
							this.debugLevel
							);
					if (nodeQueue.isEmpty()) { // nodes==null){
						//						   if (debugLevel>1) System.out.println("All start points tried");
						if (global_debug_level > (debugThreshold+1)) {
							System.out.println("All start points tried");
							int numLeft=0;
							for (boolean b:triedIndices) if (!b) numLeft++;
							System.out.println("nodeQueue.isEmpty(), startScanIndex="+startScanIndex+" numTries="+numTries+" numLeft="+numLeft);
						}
						triedIndices[numTries]=true; // all tried
					} else {
						//						   if (debugLevel>1) System.out.println("Found "+nodes.length+" candidates");
						//						   if (debugLevel>1) System.out.println("distortions: Found "+nodeQueue.size()+" candidates");
						if (global_debug_level > debugThreshold){
							System.out.println("distortions: Found "+nodeQueue.size()+" candidates");
						}
					}
				} else {
					System.out.println("All start points tried before - should not get here");
					triedIndices[numTries]=true; // all tried
				}
			}
			//				   if (global_debug_level>0) System.out.println("distortions(): startScanIndex="+startScanIndex);

			//			   if (startScanIndex[0]>=numTries) startScanIndex[0]=-1; // all indices used
			//			   if ((nodes==null) || (nodes[0]==null) || (nodes[0][0]==null)) {
			if ((nodeQueue.isEmpty()) || (nodeQueue.peek().getNode()[0]==null)) {
				if (debugLevel>1) System.out.println("*** Pattern not found");
				this.PATTERN_GRID=null;
				return 0;
			}
			if (global_debug_level> (debugThreshold + 1)) {
				System.out.println("distortions(): found "+nodeQueue.size()+" grid candidates");
				//				   System.out.println("*** distortions: Center x="+IJ.d2s(centerXY[0],3)+" y="+ 	IJ.d2s(centerXY[1],3));
			}

			debugLevel=debug_level; // ????

		} else { // create initial wave from the border nodes of existent grid
			// start with clearing all invalid nodes
			for (iUV[1]=0;iUV[1]<this.PATTERN_GRID.length;iUV[1]++) for (iUV[0]=0;iUV[0]<this.PATTERN_GRID[0].length;iUV[0]++)
				if (!isCellDefined(this.PATTERN_GRID,iUV))clearPatternGridCell(this.PATTERN_GRID,iUV);
			//		   int [] iUV=  new int [2];
			// probably start with clearing all invalid nodes
			//// 	public boolean [] focusMask=null; // array matching image pixels, used with focusing (false outside sample areas)
			int [] iUV1=new int[2];
			for (iUV[1]=0;iUV[1]<this.PATTERN_GRID.length;iUV[1]++) for (iUV[0]=0;iUV[0]<this.PATTERN_GRID[0].length;iUV[0]++)
				if (isCellDefined(this.PATTERN_GRID,iUV)){ // see if it has any new undefined neighbors
					boolean hasNewNeib=false;
					for (int dir=0;dir<directionsUV.length;dir++) {
						iUV1[0]=iUV[0]+directionsUV[dir][0];
						iUV1[1]=iUV[1]+directionsUV[dir][1];
						if (isCellNew(this.PATTERN_GRID,iUV1)){
							hasNewNeib=true;
							break;
						}
					}
					if (hasNewNeib) {
						putInWaveList(waveFrontList, iUV, 0);
					}
				}
			double [][] node={null};
			nodeQueue.add(new GridNode(node)); // will not be used, any element
		}
		int numDefinedCells=0;
		int debug_left=nodeQueue.size();
		for (GridNode gn:nodeQueue){ // trying candidates as grid seeds - until found or nothing left
			debug_left--;
			if (global_debug_level > (debugThreshold + 1)) {
				System.out.println("distortions: nodeQueue has "+(debug_left)+" candidates left (excluding this one)");
			}
			if (!updating){
				double [][] node=gn.getNode();
				double [] centerXY=node[0];
				if (global_debug_level > debugThreshold) {
					//					   System.out.println("distortions: node X/Y are "+centerXY[0]+"/"+centerXY[1]);
					System.out.println("distortions: nodeQueue has "+(debug_left)+" candidates left (excluding this one) :node X/Y are "+centerXY[0]+"/"+centerXY[1]);

				}
				//				   if (debugLevel>1) {
				if (global_debug_level > (debugThreshold + 1)) {
					System.out.println("*** distortions: Center x="+IJ.d2s(centerXY[0],3)+" y="+ 	IJ.d2s(centerXY[1],3));
					System.out.println("*** distortions: setting debugX="+IJ.d2s(centerXY[0],3)+" debugY="+ 	IJ.d2s(centerXY[1],3));
					patternDetectParameters.debugX=centerXY[0]; // Change debug coordinates to the initial node
					patternDetectParameters.debugY=centerXY[1]; //patternDetectParameters.debugRadius);
				}
				debugLevel=debug_level;
				// Reset pattern grid
				this.PATTERN_GRID=setPatternGridArray(distortionParameters.gridSize); // global to be used with threads?
				setPatternGridCell(
						this.PATTERN_GRID,
						centerUV,
						centerXY, // contrast OK?
						node[1],
						node[2]);
				waveFrontList.clear();
				putInWaveList(waveFrontList, centerUV, 0);
				if (global_debug_level > debugThreshold) { //1) {
					System.out.println("putInWaveList(waveFrontList, {"+centerUV[0]+","+centerUV[1]+"}, 0);");
				}
			}

			// Each layer processing may be multi-threaded, they join before going to the next layer
			// When looking for the next cells, the position is estimated knowing the neighbor that has wave vectors defined
			// after the layer pass is over, the wave vectors are calculated from the distances to neighbors (one or both vectors may have
			// to use those from the neighbor?
			if (debugLevel>1) System.out.println("-->centerUV= "+centerUV[0]+",  "+centerUV[1]+",  0");
			int [] uvdir; // (u,v,direction}
			int layer=0;
			int dir;
			double [][][] neibors=new double [8][][]; // uv and xy vectors to 8 neibors (some may be null
			double [][] thisCell;
			double [][] otherCell;
//			final int debugThreshold=2;
			final double [][] extrapolationWeights=generateWeights (
					distortionParameters.correlationWeightSigma,
					distortionParameters.correlationRadiusScale); //  if 0 - use sigma as radius, inside - 1.0, outside 0.0. If >0 - size of array n*sigma

			int umax,vmax,vmin,umin;
			final AtomicInteger addedCells = new AtomicInteger(0); // cells added at cleanup stage
			final AtomicBoolean cleanup=new AtomicBoolean(false); // after the wave dies, it will be restored for all cells with defined neigbors to try again. maybe - try w/o threads?

			final AtomicInteger debugCellSet= new AtomicInteger(0); // cells added at cleanup stage
			// special case (most common, actually) when initial wave has 1  node. Remove it after processing
			//first cell(s) will need large correction and so may fail during "refine", so trying to recalculate it right after the first layer)
			ArrayList<Integer> initialWave=new ArrayList<Integer>();
			for (Integer I:waveFrontList) initialWave.add(I);
			while (waveFrontList.size()>0) {
				// process current list, add new wave layer (moving in one of the 4 directions)
				// proceed until the entry is undefined on the grid (or list is empty
				while (waveFrontList.size()>0) { // will normally break out of the cycle
					uvdir= getWaveList (waveFrontList,0);
					if (this.debugLevel > (debugThreshold + 1)) System.out.println("<--uvdir= "+uvdir[0]+",  "+uvdir[1]+",  "+uvdir[2]);
					if (this.PATTERN_GRID[uvdir[1]][uvdir[0]]==null) break; // finished adding new layer
					if (!isCellDefined(this.PATTERN_GRID,uvdir)) break; // finished adding new layer, hit one of the newely added
					boolean hasNeededNeighbor=true;
					if (useFocusMask) {
						int ix=(int) Math.round(this.PATTERN_GRID[uvdir[1]][uvdir[0]][0][0]);
						int iy=(int) Math.round(this.PATTERN_GRID[uvdir[1]][uvdir[0]][0][1]);
						int indx=iy*getImageWidth()+ix;
						if ((indx<0) || (indx>focusMask.length)){
							System.out.println("distortions(): this.PATTERN_GRID["+uvdir[1]+"]["+uvdir[0]+"][0][0]="+this.PATTERN_GRID[uvdir[1]][uvdir[0]][0][0]);
							System.out.println("distortions(): this.PATTERN_GRID["+uvdir[1]+"]["+uvdir[0]+"][0][1]="+this.PATTERN_GRID[uvdir[1]][uvdir[0]][0][1]);
							System.out.println("distortions(): ix="+ix);
							System.out.println("distortions(): iy="+iy);
							System.out.println("distortions(): focusMask.length="+focusMask.length);
						}
						// TODO: find how it could get negative coordinates
						if ((ix<0) || (iy<0) || (ix>=distortionParameters.gridSize) || (iy>=distortionParameters.gridSize)) hasNeededNeighbor=false; //???
						else hasNeededNeighbor=focusMask[iy*getImageWidth()+ix]; //* OOB -1624 java.lang.ArrayIndexOutOfBoundsException: -1624,  at MatchSimulatedPattern.distortions(MatchSimulatedPattern.java:3063),  at LensAdjustment.updateFocusGrid(LensAdjustment.java:121),  at Aberration_Calibration.measurePSFMetrics(Aberration_Calibration.java:5994)
					}
					for (dir=0;dir<directionsUV.length;dir++) {
						iUV[0]=uvdir[0]+directionsUV[dir][0];
						iUV[1]=uvdir[1]+directionsUV[dir][1];

						if ((iUV[0]<0) || (iUV[1]<0) ||
								(iUV[0]>=distortionParameters.gridSize) || (iUV[1]>=distortionParameters.gridSize)) continue; // don't fit into UV grid
						if (!isCellNew(PATTERN_GRID,iUV)) continue; // already processed (or deleted!)
						// add uv and dir to the list
						// 	public boolean [] focusMask=null; // array matching image pixels, used with focusing (false outside sample areas)
						// New: if it is updating the grid and focusMask is defined - do not go more than 1 step away from the needed image area
						if (hasNeededNeighbor) {
							putInWaveList (waveFrontList, iUV, (dir+(directionsUV.length/2))%directionsUV.length); // opposite direction
							initPatternGridCell(PATTERN_GRID, iUV);
							if (debugLevel>1) System.out.println("-->iUV= "+iUV[0]+",  "+iUV[1]+",  "+((dir+(directionsUV.length/2))%directionsUV.length));
						}
					}
					waveFrontList.remove(0);  // remove first element from the list
					if (debugLevel>1) System.out.println("xx> remove(0), (waveFrontList.size()="+(waveFrontList.size()));
				}
				if (waveFrontList.size()==0) break; // not really needed?
				layer++;
				if (updateStatus) IJ.showStatus("Correlating patterns, layer "+layer+(cleanup.get()?"(cleanup)":"")+", length "+waveFrontList.size());
				//				   if (debugLevel>1) System.out.println("Correlating patterns, layer "+layer+", length "+waveFrontList.size());
				if (global_debug_level > (debugThreshold + 1)) System.out.println("Correlating patterns, layer "+layer+", length "+waveFrontList.size());
				// starting layer
				cellNum.set(0);
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						@Override
						public void run() {
							SimulationPattern simulationPattern= new SimulationPattern(bPattern);
							MatchSimulatedPattern matchSimulatedPatternCorr=new MatchSimulatedPattern(distortionParameters.correlationSize);
							DoubleFHT fht_instance =new DoubleFHT(); // provide DoubleFHT instance to save on initializations (or null)
							String dbgStr="";
							for (int ncell=cellNum.getAndIncrement(); ncell<waveFrontList.size();ncell=cellNum.getAndIncrement()){
								int [] iUVdir=getWaveList (waveFrontList,ncell);
								if (debugLevel>debugThreshold) {
									dbgStr="";
									dbgStr+="<--iUVdir= "+iUVdir[0]+",  "+iUVdir[1]+",  "+iUVdir[2];
								}
								int [] iUVRef=new int[2];
								iUVRef[0]=iUVdir[0]+directionsUV[iUVdir[2]][0];
								iUVRef[1]=iUVdir[1]+directionsUV[iUVdir[2]][1];
								// refCell - is where it came from, but if the initials are disabled, it is null

								double [][] refCell=PATTERN_GRID[iUVRef[1]][iUVRef[0]]; // should never be null as it is an old one
								if (refCell==null){
									System.out.println("**** refCell==null - what does it mean?**** u="+iUVRef[0]+" v="+iUVRef[1]+
											" current="+iUVdir[0]+"/"+iUVdir[1]+" len="+iUVdir.length);
									continue;
								} else if ((refCell[0]!=null) && (refCell[0].length>3)){
									double dbg_contrast=(refCell[0].length>2)?refCell[0][2]:Double.NaN;
									System.out.println("**** refCell was deleted **** u="+iUVRef[0]+" v="+iUVRef[1]+
											" current="+iUVdir[0]+"/"+iUVdir[1]+
											" ncell="+ncell+" waveFrontList.size()="+waveFrontList.size()+
											" ref_x="+IJ.d2s(refCell[0][0],3)+" ref_y="+IJ.d2s(refCell[0][1],3)+
											" contrast="+IJ.d2s(dbg_contrast,3));
								}
								//found reference cell, calculate x/y, make sure it is inside the selection w/o borders
								double [][] wv=new double [2][];
								wv[0]=refCell[1];
								wv[1]=refCell[2];
								double [][]uv2xy=matrix2x2_invert(wv);

								double [] dUV = new double [2];
								dUV[0]=0.5*(iUVdir[0]-iUVRef[0]);
								dUV[1]=0.5*(iUVdir[1]-iUVRef[1]);
								double []dXY=matrix2x2_mul(uv2xy,dUV);
								double [] expectedXY=matrix2x2_add(refCell[0],dXY);

								// Try new extrapolation, debug print both
								//extrapolationWeights
								double [][] estimatedCell=estimateCell(
										PATTERN_GRID,
										iUVdir,
										extrapolationWeights, // quadrant of sample weights
										true, // useContrast
										!distortionParameters.useQuadratic,  // use linear approximation (instead of quadratic)
										1.0E-10,  // thershold ratio of matrix determinant to norm for linear approximation (det too low - fail)
										1.0E-20  // thershold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
										);
								double [][] simulPars=null;
								if (debugLevel>debugThreshold) {
									dbgStr+=" ExpectedXY(old)= "+IJ.d2s(expectedXY[0],3)+" / "+IJ.d2s(expectedXY[1],3)+",  "+
											" vw00="+IJ.d2s(wv[0][0],5)+" vw01="+IJ.d2s(wv[0][1],5)+
											" vw10="+IJ.d2s(wv[1][0],5)+" vw11="+IJ.d2s(wv[1][1],5);
									if (estimatedCell==null) {
										dbgStr+=" -- ExpectedXY(new)= ***** NULL **** ";
									} else {
										dbgStr+=" -- ExpectedXY(new)= "+IJ.d2s(estimatedCell[0][0],3)+" / "+IJ.d2s(estimatedCell[0][1],3)+",  "+
												" vw00="+IJ.d2s(estimatedCell[1][0],5)+" vw01="+IJ.d2s(estimatedCell[1][1],5)+
												" vw10="+IJ.d2s(estimatedCell[2][0],5)+" vw11="+IJ.d2s(estimatedCell[2][1],5);
									}


								}
								if (estimatedCell!=null) {
									expectedXY=estimatedCell[0];
									wv[0]=estimatedCell[1];
									wv[1]=estimatedCell[2];

									simulPars=getSimulationParametersFromGrid(
											PATTERN_GRID,
											iUVdir,          // U,V of the center point (for which the simulation pattern should be built
											expectedXY,          // x,y of the center point (or null to use grid)
											extrapolationWeights, // quadrant of sample weights
											!distortionParameters.useQuadratic,  // use linear approximation (instead of quadratic)
											1.0E-10,  // thershold ratio of matrix determinant to norm for linear approximation (det too low - fail)
											1.0E-20  // thershold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
											);
									if (debugLevel>debugThreshold) {
										dbgStr+=" {"+IJ.d2s(simulPars[0][0],5)+"/"+IJ.d2s(simulPars[0][1],5)+"/"+IJ.d2s(simulPars[0][2],5);
										if (simulPars[0].length>3) dbgStr+="/"+IJ.d2s(simulPars[0][3],7)+"/"+IJ.d2s(simulPars[0][4],7)+"/"+IJ.d2s(simulPars[0][5],7)+"}";
										dbgStr+=" {"+IJ.d2s(simulPars[1][0],5)+"/"+IJ.d2s(simulPars[1][1],5)+"/"+IJ.d2s(simulPars[1][2],5);
										if (simulPars[1].length>3) dbgStr+="/"+IJ.d2s(simulPars[1][3],7)+"/"+IJ.d2s(simulPars[1][4],7)+"/"+IJ.d2s(simulPars[1][5],7)+"}";
									}
								}
								if (!selection.contains((int) Math.round(expectedXY[0]),(int) Math.round(expectedXY[1]))) { // just the center point
									invalidatePatternGridCell(
											PATTERN_GRID,
											iUVdir);
									if (debugLevel>debugThreshold) {
										dbgStr+=" -- not in selection ";
										System.out.println(dbgStr);
									}
									continue; // the correlation selection does not fit into WOI selection
								}
								//Proceed with correlation
								//TODO: add contrast verification ? Maximal distance from expected? (return null if failed)
								double [] centerXY=correctedPatternCrossLocation(
										lwirReaderParameters,
										expectedXY, // initial coordinates of the pattern cross point
										wv[0][0],
										wv[0][1],
										wv[1][0],
										wv[1][1],
										simulPars,
										imp,      // image data (Bayer mosaic)
										distortionParameters, //
										patternDetectParameters,
										matchSimulatedPatternCorr, // correlationSize
										thisSimulParameters,
										equalizeGreens,
										windowFunctionCorr,
										windowFunctionCorr2,
										windowFunctionCorr4,
										simulationPattern,
										((iUVdir[0]^iUVdir[1])&1)!=0, // if true - invert pattern
										fht_instance,
										distortionParameters.fastCorrelationOnFirstPass,
										locsNeib,
										debugLevel,
										null);
								//	    						 System.out.println("*+*debugLevel="+debugLevel);
								if (centerXY==null){
									invalidatePatternGridCell(
											PATTERN_GRID,
											iUVdir);
									if (debugLevel>debugThreshold) {
										dbgStr+=" -- FAILED";
										System.out.println(dbgStr);
									}
									continue; // failed to find pattern in the cell TODO: implement
								}
								if (debugCellSet.getAndIncrement()==0){ // First cell
									if (passNumber==1) {
										debugUV[0]=iUVdir[0];
										debugUV[1]=iUVdir[1];
										if (debugLevel>debugThreshold) System.out.println("debugUV[] set to {"+debugUV[0]+","+debugUV[1]+"}");
										passNumber=2; // global passNumber
									}
								}
								//Found new cell, save info and increment counter
								setPatternGridCell(
										PATTERN_GRID,
										iUVdir,
										centerXY,
										// specify wave vectors from the parent cell, will recalculate (if possible)
										wv[0], //null, //  double [] wv1,
										wv[1]); //null); //  double [] wv2);
								if (cleanup.get()) addedCells.getAndIncrement();
								if (debugLevel>debugThreshold) { //was no "-2"
									dbgStr+="==>added"+iUVdir[0]+"/"+iUVdir[1]+", dir"+iUVdir[2];
									System.out.println(dbgStr);
								}

							}
						}
					};
				}
				startAndJoin(threads);
				// remove invalid cells from the list
				for (int i=waveFrontList.size()-1;i>=0;i--) {
					if (!isCellValid(PATTERN_GRID,getWaveList (waveFrontList,i))) {
						// Make that cell "new", so it will be tried again, until wave will not touch it. So when more neigbors will be defined, previously failed
						// cell will be retried
						clearPatternGridCell(PATTERN_GRID,getWaveList (waveFrontList,i));
						waveFrontList.remove(i);
						if (debugLevel > (debugThreshold + 1)) System.out.println("XXX->clear invalid ("+i+")");
					}
				}
				//If anything was added during the layer - calculate and fill in wave vectors here (they are set to the same as in the parent cell)
				// this code is not needed now, the wave vectors are recalculated from x/y locations, the stored ones are not used
				if (waveFrontList.size()>0) {
					for (int listIndex=0;listIndex<waveFrontList.size();listIndex++) {
						uvdir= getWaveList (waveFrontList,listIndex);
						if (debugLevel > (debugThreshold + 1)) System.out.println("<---= uvdir= "+uvdir[0]+",  "+uvdir[1]+",  "+uvdir[2]);
						thisCell=PATTERN_GRID[uvdir[1]][uvdir[0]];
						neibBits=0;
						for (dir=0;dir<directionsUV8.length;dir++) {
							neibors[dir]=null;
							iUV[0]=uvdir[0]+directionsUV8[dir][0];
							iUV[1]=uvdir[1]+directionsUV8[dir][1];
							if ((iUV[0]<0) || (iUV[1]<0) ||
									(iUV[0]>=distortionParameters.gridSize) || (iUV[1]>=distortionParameters.gridSize)) continue; // don't fit into UV grid
							if (isCellValid(PATTERN_GRID,iUV)) {
								neibors[dir]= new double [2][2];
								otherCell=PATTERN_GRID[iUV[1]][iUV[0]];
								neibors[dir][0][0]=0.5*directionsUV8[dir][0];  // u
								neibors[dir][0][1]=0.5*directionsUV8[dir][1];  // v
								neibors[dir][1][0]=otherCell[0][0]-thisCell[0][0];  // x
								neibors[dir][1][1]=otherCell[0][1]-thisCell[0][1];  // y
								neibBits |= directionsBits8[dir];
							}
						}
						int i=Integer.bitCount(neibBits);
						if (debugLevel > (debugThreshold + 1)) System.out.println("neibBits="+neibBits+", number of bits= "+i);
						if (i>1) {
							double[][] wv= waveVectorsFromNeib(neibors);
							setPatternGridCell(
									PATTERN_GRID,
									uvdir,
									null, // XY already set
									wv[0],
									wv[1]);

							if (debugLevel > (debugThreshold + 1)) System.out.println("==+> number of bits:"+i+
									" vw00="+IJ.d2s(wv[0][0],5)+" vw01="+IJ.d2s(wv[0][1],5)+
									" vw10="+IJ.d2s(wv[1][0],5)+" vw11="+IJ.d2s(wv[1][1],5));
							//				    	wv= WaveVectorsFromNeib(neibors);
							//
							// 		   vectors: [num_vector][0][0] - U
							//                  [num_vector][0][1] - V
							//                  [num_vector][1][0] - X
							//                  [num_vector][1][1] - Y
							//                  [num_vector] == null - skip
							//
						}

					}
				} else if (initialWave!=null){
					if ((global_debug_level > (debugThreshold + 1)) && (initialWave!=null)) {
						System.out.println("No sense to initiate clenaup during first layer"); // problems heer?
					}
				} else if (!cleanup.get() || (addedCells.get()>0)) { // create list of the defined cells on the border (if wave died)
					cleanup.set(true);
					// debug
					//					   if ((global_debug_level>0) && (initialWave!=null)) {
					//						   System.out.println("clenaup during first layer"); // problems heer?
					//						   System.out.println("Added "+addedCells.get()+" during border cleanup on first layer");
					//					   }
					if ((debugLevel > (debugThreshold + 1)) && !cleanup.get())  System.out.println("Added "+addedCells.get()+" during border cleanup"); // can not get here
					addedCells.set(0);
					umax=0;
					vmax=0;
					vmin=PATTERN_GRID.length;
					umin=PATTERN_GRID[0].length;
					for (int i=0;i<PATTERN_GRID.length;i++) for (int j=0;j<PATTERN_GRID[i].length;j++) {
						if ((PATTERN_GRID[i][j]!=null) && (PATTERN_GRID[i][j][0]!=null)) {
							if (vmin > i) vmin = i;
							if (vmax < i) vmax = i;
							if (umin > j) umin = j;
							if (umax < j) umax = j;
						}
					}
					int [] uvNew=new int [2];
					for (uvNew[1]=vmin;uvNew[1]<=vmax;uvNew[1]++) for (uvNew[0]=umin;uvNew[0]<=umax;uvNew[0]++) if (isCellDefined(PATTERN_GRID,uvNew)){
						for (dir=0;dir<directionsUV.length;dir++) {
							iUV[0]=uvNew[0]+directionsUV[dir][0];
							iUV[1]=uvNew[1]+directionsUV[dir][1];
							if (!isCellDefined(PATTERN_GRID,iUV) && !isCellDeleted(PATTERN_GRID,iUV)){
								putInWaveList (waveFrontList, uvNew, dir); // direction does not matter here
								break;
							}
						}
					}
					if (global_debug_level > (debugThreshold + 1)) System.out.println("***** Starting cleanup, wave length="+waveFrontList.size()); //????
				}
				// end of layer - it is a hack below, marking initial wave to recalculate it from neighbors
				if (initialWave!=null){ // just after the first layer (usually one cell) - delete it and add next time - otherwise first one needs large correction
					if (global_debug_level > (debugThreshold + 1)) {
						System.out.println("Removing "+initialWave.size()+" initial wave cells, waveFrontList.size()="+waveFrontList.size());
						for (int listIndex=0;listIndex<waveFrontList.size();listIndex++) {
							int [] dbg_uvdir= getWaveList (waveFrontList,listIndex);
							System.out.println("waveFrontList["+listIndex+"]: "+dbg_uvdir[0]+"/"+dbg_uvdir[1]+" dir="+dbg_uvdir[2]);
						}
					}

					while (initialWave.size()>0){
						uvdir= getWaveList (initialWave,0);
						//						   clearPatternGridCell(PATTERN_GRID, uvdir);
						if (global_debug_level > (debugThreshold + 1))
							System.out.println("Removing x="+uvdir[0]+" y="+uvdir[1]+" dir="+uvdir[2]);
						markDeletedPatternGridCell(PATTERN_GRID, uvdir);
						initialWave.remove(0);
					}
					initialWave=null;
				}
			}//while (waveFrontList.size()>0)
			debugLevel=was_debug_level;
			/*
		   if (updating){
			   return PATTERN_GRID; // no need to crop the array, it should not change
		   }
			 */
			umax=0;
			vmax=0;
			vmin=PATTERN_GRID.length;
			umin=PATTERN_GRID[0].length;
			numDefinedCells=0;
			for (int i=0;i<PATTERN_GRID.length;i++) for (int j=0;j<PATTERN_GRID[i].length;j++) {
				if ((PATTERN_GRID[i][j]!=null) && (PATTERN_GRID[i][j][0]!=null)) {
					if (vmin > i) vmin = i;
					if (vmax < i) vmax = i;
					if (umin > j) umin = j;
					if (umax < j) umax = j;
					numDefinedCells++;
				}
			}

			//			   if (updating){
			//				   return numDefinedCells; // no need to crop the array, it should not change
			//			   }
			if (!updating){
				if (vmin>vmax){
					this.PATTERN_GRID=null;
					continue; // try next in queue if available
					//					   return 0; // null; // nothing found
				}
				// Add extra margins for future extrapolation
				int extra=distortionParameters.numberExtrapolated -((distortionParameters.removeLast)?1:0);
				vmin-=extra;
				if (vmin<0) vmin=0;
				umin-=extra;
				if (umin<0) umin=0;
				vmax+=extra;
				if (vmax>=PATTERN_GRID.length) vmax=PATTERN_GRID.length-1;
				umax+=extra;
				if (umax>=PATTERN_GRID[0].length) umax=PATTERN_GRID[0].length-1;

				// make sure the odd/even uv does not change (and so the cross phases defined by U & V)
				vmin &= ~1;
				umin &= ~1;
				// make width/height even (not needed)
				umax |=1;
				vmax |=1;
				// remove  margins
				double [][][][] result = new double [vmax-vmin+1][umax-umin+1][][];
				for (int i=vmin;i<=vmax;i++) for (int j=umin;j<=umax;j++) {
					if ((PATTERN_GRID[i][j]!=null) && (PATTERN_GRID[i][j][0]!=null)) {
						result[i-vmin][j-umin]=PATTERN_GRID[i][j];
					} else result[i-vmin][j-umin]=null;
				}
				this.debugUV[0]-=umin;
				this.debugUV[1]-=umin;
				if (debugLevel>(debugThreshold+2)) System.out.println("debugUV[] updated to {"+this.debugUV[0]+","+this.debugUV[1]+"}");

				if (debugLevel > (debugThreshold+1)) System.out.println("Total number of defined cells="+numDefinedCells);
				this.PATTERN_GRID=result;
			}
			// more tests here (moved from the caller) that result is good
			double averageGridPeriod=Double.NaN;
			double [] gridPeriods={Double.NaN,Double.NaN};
			if (this.PATTERN_GRID!=null) {
				averageGridPeriod=averageGridPeriod(this.PATTERN_GRID);
				gridPeriods=averageGridPeriods(this.PATTERN_GRID); // {min,max}
			}
			if (debugLevel > debugThreshold){
				System.out.println("Pattern period="+averageGridPeriod+" {"+gridPeriods[0]+","+gridPeriods[1]+"}"+
						" limits are set to :"+patternDetectParameters.minGridPeriod+","+patternDetectParameters.maxGridPeriod);
			}
			if (!Double.isNaN(averageGridPeriod)) {
				if (!Double.isNaN(patternDetectParameters.minGridPeriod) &&
						(patternDetectParameters.minGridPeriod>0.0) &&
						(averageGridPeriod<patternDetectParameters.minGridPeriod)){
					if (debugLevel>0){
						System.out.println("Pattern is too small, period="+averageGridPeriod+
								" minimal="+patternDetectParameters.minGridPeriod);
					}
					continue; // bad grid
				}
				if (!Double.isNaN(patternDetectParameters.maxGridPeriod) &&
						(patternDetectParameters.maxGridPeriod>0.0) &&
						(averageGridPeriod>patternDetectParameters.maxGridPeriod)){
					if (debugLevel>0){
						System.out.println("Pattern is too large, period="+averageGridPeriod+
								" maximal="+patternDetectParameters.maxGridPeriod);
					}
					continue; // bad grid
				}
			}
			if (
					(minimal_pattern_cluster <= 0) || // minimal cluster size is disabled
					(distortionParameters.scaleMinimalInitialContrast<=0) || // minimal cluster size is disabled
					((numDefinedCells  == 0) && fromVeryBeginning)|| // no cells detected at all, starting from the very beginning
					(numDefinedCells   >= minimal_pattern_cluster)  // detected enough cells
					)  {
				return numDefinedCells;
			}
			if ( (numDefinedCells < minimal_pattern_cluster) &&
					(numDefinedCells > 10))   // detected enough cells
			{
				if (global_debug_level > (debugThreshold-1)) {
					System.out.println("***** Initial cluster has "+numDefinedCells+ " cells that is less than "+
							"minimal_pattern_cluster = "+minimal_pattern_cluster+
							" *****");
				}
				//				   return numDefinedCells;
			}


			if (roi!=null){ // don't use this feature with ROI as it can be small
				if (global_debug_level>0) System.out.println("Initial pattern cluster is small ("+numDefinedCells+"), but ROI is set - no retries");
				{
					return numDefinedCells;
				}
			}
		} // next node in queue
		return 0; // none
	}
	// ================= end of public int distortions() ===================

	public double [][] findPatternCandidate_old(
			//			   final int [] startScanIndex, // [0] will be updated
			final boolean [] triedIndices, // which indices are already tried
			final int startScanIndex,
			final int tryHor,
			final int tryVert,
			//			   final int numTries,
			final Rectangle selection,
			final DistortionParameters distortionParameters, //
			final MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			final double    min_half_period,
			final double    max_half_period,
			final SimulationPattern.SimulParameters  thisSimulParameters,
			final MatchSimulatedPattern matchSimulatedPattern,
			final MatchSimulatedPattern matchSimulatedPatternCorr,
			final SimulationPattern simulationPattern,
			final boolean equalizeGreens,
			final ImagePlus imp, // image to process
			final double [] bPattern,
			final double [] windowFunction,
			final double [] windowFunctionCorr,
			final double [] windowFunctionCorr2,
			final double [] windowFunctionCorr4,
			final double[][] locsNeib, // which neibors to try (here - just the center)
			final int threadsMax,
			final boolean updateStatus,
			final int debugLevel
			){
		final Thread[] threads = newThreadArray(threadsMax);
		//		   final AtomicInteger seqNumber = new AtomicInteger(3);
		//		   final AtomicInteger seqNumber = new AtomicInteger(startScanIndex[0]);
		final AtomicInteger seqNumber = new AtomicInteger(startScanIndex);
		//		   startScanIndex
		final AtomicBoolean nodeSet=new AtomicBoolean(false);
		final double [][][] nodeRef= new double[1][][];
		nodeRef[0]=null;
		//		   System.out.println("===== findPatternCandidate(): startScanIndex="+startScanIndex);
		//		   for (int i=0;i<triedIndices.length;i++) System.out.print(triedIndices[i]?"+":"-");
		//		   System.out.println();
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					int nbh, nbv, nh, nv, nb;
					double [] point = new double[2];
					DoubleFHT doubleFHT = new DoubleFHT();
					//					   for (int n=seqNumber.getAndIncrement(); n< numTries;n=seqNumber.getAndIncrement()){
					for (int n=seqNumber.getAndIncrement(); n<(triedIndices.length-1); n=seqNumber.getAndIncrement()) if (!triedIndices[n]){
						if (nodeSet.get()) break; // already set
						nbh=tryHor-1;
						nbv=tryVert-1;
						nh=0;
						nv=0;
						nb=0;
						while (nb<(tryHor+tryVert)) {
							if (nbh>=0) {
								if ((n & (1<<nb))!=0) nh |= 1<<nbh;
								nbh--;
								nb++;
							}
							if (nbv>=0) {
								if ((n & (1<<nb))!=0) nv |= 1<<nbv;
								nbv--;
								nb++;
							}
						}
						if (debugLevel>2) System.out.println("Searching, n="+n+", nv="+nv+", nh="+nh+", nb="+nb);
						if ((nv>0) && (nh>0)) {
							point[0]=(selection.x+nh*selection.width/(1<<tryHor)) & ~1;
							point[1]=(selection.y+nv*selection.height/(1<<tryVert)) & ~1;
							if (debugLevel>2) System.out.println("trying xc="+point[0]+", yc="+point[1]+"(nv="+nv+", nh="+nh+")");
							//							   if ((debugLevel>2) && (n==3)) debugLevel=3; // show debug images for the first point
							double [][] node=tryPattern (
									null, // LwirReaderParameters lwirReaderParameters, // null is OK
									doubleFHT,
									point, // xy to try
									distortionParameters, //no control of the displacement
									patternDetectParameters,
									min_half_period,
									max_half_period,
									thisSimulParameters,
									matchSimulatedPattern,
									matchSimulatedPatternCorr,
									simulationPattern,
									equalizeGreens,
									imp, // image to process
									bPattern,
									windowFunction,
									windowFunctionCorr,
									windowFunctionCorr2,
									windowFunctionCorr4,
									locsNeib, // which neibors to try (here - just the center)
									null // dbgStr
									);
							if ((node!=null) && (node[0]!=null)) {
								if (nodeSet.compareAndSet(false,true)) {
									nodeRef[0]=node;
									//									   startScanIndex[0]=seqNumber.get();
									triedIndices[n]=true; //found and will be processed
									if (debugLevel>1)  System.out.println("probing "+n);
									break;
								} else {
									if (debugLevel>1)  System.out.println("missed "+n);
								}
							} else {
								triedIndices[n]=true; // tried, but nothing found
								//								   System.out.println("empty "+n);
							}
						} else {
							triedIndices[n]=true; // tried, but nothing found
							//							   System.out.println("wrong "+n);
						}
					}
				}
			};
		}
		startAndJoin(threads);
		//		   if (nodeRef[0]==null) startScanIndex[0]=numTries; // all used
		return nodeRef[0];
	}

	class GridNode {
		double [][] node;
		public GridNode(double [][] node){
			this.node=node;
		}
		public double [][] getNode(){
			return this.node;
		}
	}

	private  Queue<GridNode> findPatternCandidates(
			final LwirReaderParameters lwirReaderParameters, // null is OK
			final boolean [] triedIndices, // which indices are already tried
			final int startScanIndex,
			final int tryHor,
			final int tryVert,
			final Rectangle selection,
			final DistortionParameters distortionParameters, //
			final MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			final double    min_half_period,
			final double    max_half_period,
			final SimulationPattern.SimulParameters  thisSimulParameters,
			final MatchSimulatedPattern matchSimulatedPattern,
			final MatchSimulatedPattern matchSimulatedPatternCorr,
			final SimulationPattern simulationPattern,
			final boolean equalizeGreens,
			final ImagePlus imp, // image to process
			final double [] bPattern,
			final double [] windowFunction,
			final double [] windowFunctionCorr,
			final double [] windowFunctionCorr2,
			final double [] windowFunctionCorr4,
			final double[][] locsNeib, // which neibors to try (here - just the center)
			final int threadsMax,
			final boolean updateStatus,
			final int debugLevel
			){
		final int debugThreshold=2; // -1; // 1; ** Restore 1
		if (debugThreshold < 0) {
			System.out.println("findPatternCandidates(): debugThreshold < 0 - restore when done");
		}
		if ((debugLevel>debugThreshold) && ((debugLevel>1) || (startScanIndex>3))) {
			int debugNumLeft=0;
			for (boolean b:triedIndices) if (!b) debugNumLeft++;
			System.out.println("findPatternCandidates(), startScanIndex= "+startScanIndex+",triedIndices.length="+triedIndices.length+" debugNumLeft="+debugNumLeft);
		}

		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger seqNumber = new AtomicInteger(startScanIndex);

		final AtomicInteger debugNumThreadAtomic= new AtomicInteger(0);

		final Queue<GridNode> nodeQueue = new ConcurrentLinkedQueue<GridNode>();
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					int nbh, nbv, nh, nv, nb;
					double [] point = new double[2];
					DoubleFHT doubleFHT = new DoubleFHT();
					int debugNumThread=debugNumThreadAtomic.getAndIncrement();
					for (int n=seqNumber.getAndIncrement(); n<(triedIndices.length-1); n=seqNumber.getAndIncrement()) if (!triedIndices[n]){
						if (!nodeQueue.isEmpty()) break; // already set at least one element - does it work?
						nbh=tryHor-1;
						nbv=tryVert-1;
						nh=0;
						nv=0;
						nb=0;
						while (nb<(tryHor+tryVert)) {
							if (nbh>=0) {
								if ((n & (1<<nb))!=0) nh |= 1<<nbh;
								nbh--;
								nb++;
							}
							if (nbv>=0) {
								if ((n & (1<<nb))!=0) nv |= 1<<nbv;
								nbv--;
								nb++;
							}
						}
						if (debugLevel>2) System.out.println("Searching, n="+n+", nv="+nv+", nh="+nh+", nb="+nb );
						if ((nv>0) && (nh>0)) {
							point[0]=(selection.x+nh*selection.width/(1<<tryHor)) & ~1;
							point[1]=(selection.y+nv*selection.height/(1<<tryVert)) & ~1;
							if (debugLevel>2) System.out.println("trying xc="+point[0]+", yc="+point[1]+"(nv="+nv+", nh="+nh+")");
							if (debugLevel>2) System.out.println(debugNumThread+":"+n+" >> ");
							double [][] node=tryPattern (
									lwirReaderParameters, // LwirReaderParameters lwirReaderParameters, // null is OK
									doubleFHT,
									point, // xy to try
									distortionParameters, //no control of the displacement
									patternDetectParameters,
									min_half_period,
									max_half_period,
									thisSimulParameters,
									matchSimulatedPattern,
									matchSimulatedPatternCorr,
									simulationPattern,
									equalizeGreens,
									imp, // image to process
									bPattern,
									windowFunction,
									windowFunctionCorr,
									windowFunctionCorr2,
									windowFunctionCorr4,
									locsNeib, // which neighbors to try (here - just the center)
									(debugLevel>debugThreshold)?(""+debugNumThread+":"+n+", nv="+nv+", nh="+nh+", nb="+nb+" "+point[0]+"/"+point[1]):null
									);
							if ((node!=null) && (node[0]!=null)) {
								nodeQueue.add(new GridNode(node));
								if (debugLevel>debugThreshold-1)  System.out.println("adding candidate "+n+" x0="+point[0]+" y0="+point[1]+" -> "+ node[0][0]+"/"+node[0][1]+" seqNumber.get()="+seqNumber.get()+" n="+n);
							}
						} else {
							if (debugLevel>debugThreshold) System.out.println("-----"+debugNumThread+":"+n+", nv="+nv+", nh="+nh);
						}
						triedIndices[n]=true; // regardless - good or bad
					}
				}
			};
		}
		startAndJoin(threads);
		if (debugLevel>debugThreshold){
			System.out.println("seqNumber after join is "+seqNumber.get());
		}
		if (seqNumber.get()>=(triedIndices.length-1) ) triedIndices[triedIndices.length-1]=true; // all tried
		return nodeQueue; // never null, may be empty
	}

	/* ================================================================*/
	public void scaleContrast(double scale){
		for (double [][][] patternRow:this.PATTERN_GRID){
			if (patternRow!=null) for (double [][] node:patternRow){
				if ((node!=null) && (node.length>0) && (node[0]!=null) && (node[0].length>2)) {
					node[0][2]*=scale;
				}
			}
		}
	}
	/* ================================================================*/
	public double refineDistortionCorrelation (
			final LwirReaderParameters lwirReaderParameters, // null is OK
			final DistortionParameters distortionParameters, //
			final MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			final SimulationPattern.SimulParameters  simulParameters,
			final boolean equalizeGreens,
			final ImagePlus imp, // image to process
			final double maxCorr, // maximal allowed correction, in pixels (0.0) - any
			final int threadsMax,
			final boolean updateStatus,
			final int debug_level){// debug level used inside loops
		scaleContrast(distortionParameters.scaleFirstPassContrast);
		final double [][][][] patternGrid=this.PATTERN_GRID;
		final int debugThreshold=1;
		final Rectangle selection=new Rectangle(0, 0, imp.getWidth(), imp.getHeight());
		MatchSimulatedPattern matchSimulatedPatternCorr=new MatchSimulatedPattern(distortionParameters.correlationSize);
		matchSimulatedPatternCorr.debugLevel=debugLevel;
		SimulationPattern simulationPattern= new SimulationPattern();
		final SimulationPattern.SimulParameters  thisSimulParameters=simulParameters.clone();
		thisSimulParameters.subdiv=distortionParameters.patternSubdiv;
		final double [] bPattern= simulationPattern.patternGenerator(simulParameters); // reuse pattern for next time
		/*
			final double [] windowFunctionCorr= initWindowFunction(  distortionParameters.correlationSize,distortionParameters.correlationGaussWidth);
			final double [] windowFunctionCorr2=initWindowFunction(2*distortionParameters.correlationSize,
					 (distortionParameters.absoluteCorrelationGaussWidth?0.5:1.0)*distortionParameters.correlationGaussWidth);
			final double [] windowFunctionCorr4=initWindowFunction(4*distortionParameters.correlationSize,
					 (distortionParameters.absoluteCorrelationGaussWidth?0.25:1.0)*distortionParameters.correlationGaussWidth);
		 */
		final double [] windowFunctionCorr= initWindowFunction(
				distortionParameters.correlationSize,
				distortionParameters.correlationGaussWidth,
				distortionParameters.zeros);
		final double [] windowFunctionCorr2=initWindowFunction(
				2*distortionParameters.correlationSize,
				(distortionParameters.absoluteCorrelationGaussWidth?0.5:1.0)*distortionParameters.correlationGaussWidth,
				distortionParameters.zeros);
		final double [] windowFunctionCorr4=initWindowFunction(
				4*distortionParameters.correlationSize,
				(distortionParameters.absoluteCorrelationGaussWidth?0.25:1.0)*distortionParameters.correlationGaussWidth,
				distortionParameters.zeros);
		final int height=patternGrid.length;
		final int width=(height>0)?patternGrid[0].length:0; // oob 0??
		final Thread[] threads = newThreadArray(threadsMax);
		int was_debug_level=debugLevel;
		final int debugOnLevel=(debug_level>0)?3:0;
		final double [][] locsNeib=calcNeibLocsWeights (
				distortionParameters,
				distortionParameters.correlationAverageOnRefine);
		debugLevel=debug_level;
		if (debugLevel>1)  System.out.println("Refining correlations, width= "+width+", height= "+height);
		final double [][] extrapolationWeights=generateWeights (
				distortionParameters.correlationWeightSigma,
				distortionParameters.correlationRadiusScale); //  if 0 - use sigma as radius, inside - 1.0, outside 0.0. If >0 - size of array n*sigma
		int i=-1;
		int [] iUV=  new int [2];

		for (i=0;i<(width*height);i++) {
			iUV[0]=i % width;
			iUV[1]=i / width;
			if (isCellDefined(patternGrid,iUV)) break;
		}
		if (i<0) return Double.NaN; // no defined nodes at all
		final int startCell=i;
		final AtomicInteger cellNum = new AtomicInteger(startCell);
		final AtomicInteger cellNumDoneAtomic = new AtomicInteger(startCell);
		int dc=i;
		//Debug only
		if (debugOnLevel>debugThreshold) {
			final int [][] directionsUV8= {{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}}; // first 8 should be the same as in directionsUV
			for (i=dc;i<(width*height);i++) {
				iUV[0]=i % width;
				iUV[1]=i / width;
				if ((iUV[0]>0) && (iUV[1]>0) && (iUV[0]<(width-1)) && (iUV[1]<(height-1)) &&
						isCellDefined(patternGrid,iUV) &&
						isCellDefined(patternGrid,matrix2x2_add(iUV,directionsUV8[0])) &&
						isCellDefined(patternGrid,matrix2x2_add(iUV,directionsUV8[1])) &&
						isCellDefined(patternGrid,matrix2x2_add(iUV,directionsUV8[2])) &&
						isCellDefined(patternGrid,matrix2x2_add(iUV,directionsUV8[3])) &&
						isCellDefined(patternGrid,matrix2x2_add(iUV,directionsUV8[4])) &&
						isCellDefined(patternGrid,matrix2x2_add(iUV,directionsUV8[5])) &&
						isCellDefined(patternGrid,matrix2x2_add(iUV,directionsUV8[6])) &&
						isCellDefined(patternGrid,matrix2x2_add(iUV,directionsUV8[7]))) {
					dc=i;
					System.out.println("not used: debug U="+iUV[0]+" V="+iUV[1]+" index="+dc);
					break;
				}
			}
		}

		final int debugCell=debugUV[0]+debugUV[1]*width;
		/**
		 * That was wrong to update currently calculated grid, so the newGrid will be calculated instead, then copied altogether
		 */
		final double [][][] newGrid=new double [height][width][];
		final boolean refineInPlace=distortionParameters.refineInPlace;
		for (int v=0;v<height;v++) for (int u=0;u<width;u++) newGrid[v][u]=null;
		IJ.showProgress(0);
		if (updateStatus) IJ.showStatus("Refining correlations");
		//           MinMaxSync  minMaxSync=new MinMaxSync();
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					SimulationPattern simulationPattern= new SimulationPattern(bPattern);
					MatchSimulatedPattern matchSimulatedPatternCorr=new MatchSimulatedPattern(distortionParameters.correlationSize);
					DoubleFHT fht_instance =new DoubleFHT(); // provide DoubleFHT instance to save on initializations (or null)
					int [] iUV=  new int [2];
					boolean nowDebugCell=false;
					for (int ncell=cellNum.getAndIncrement(); ncell<(width*height);ncell=cellNum.getAndIncrement()){
						nowDebugCell=(ncell==debugCell);
						int thisDebug= (nowDebugCell)?debugOnLevel:debugLevel;
						iUV[0]=ncell % width;
						iUV[1]=ncell / width;

						if (nowDebugCell && (thisDebug>1))  System.out.println(">>>>>>>>>>> Debug cell, thisDebug="+thisDebug+", iUV={"+iUV[0]+","+iUV[1]+"}");
						//        				   if ((updateStatus) && (iUV[0]==0)) IJ.showStatus("Refining correlations, row "+(iUV[1]+1)+" of "+height);
						if ((thisDebug>1) && ((iUV[0]==0) || (nowDebugCell))) System.out.println("Refining correlations, row "+(iUV[1]+1)+" of "+height);
						if (!isCellDefined(patternGrid,iUV)) {
							cellNumDoneAtomic.getAndIncrement();
							continue;
						}
						Rectangle centerCross=correlationSelection(
								patternGrid[iUV[1]][iUV[0]][0], // initial coordinates of the pattern cross point
								distortionParameters.correlationSize);
						if (!selection.contains(centerCross)) {
							cellNumDoneAtomic.getAndIncrement();
							continue; // the correlation selection does not fit into WOI selection ??? WOI is now full image
						}
						//Proceed with correlation
						//TODO: add contrast verification ? Maximal distance from expected? (return null if failed)

						double [][] simulPars=getSimulationParametersFromGrid(
								PATTERN_GRID,
								iUV,          // U,V of the center point (for which the simulation pattern should be built
								null,          // x,y of the center point (or null to use grid)
								extrapolationWeights, // quadrant of sample weights
								!distortionParameters.useQuadratic,  // use linear approximation (instead of quadratic)
								1.0E-10,  // thershold ratio of matrix determinant to norm for linear approximation (det too low - fail)
								1.0E-20  // thershold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
								);
						if ((thisDebug>debugThreshold) && (simulPars!=null)) {
							String dbgStr="";
							dbgStr+=" {"+IJ.d2s(simulPars[0][0],5)+"/"+IJ.d2s(simulPars[0][1],5)+"/"+IJ.d2s(simulPars[0][2],5);
							if (simulPars[0].length>3) dbgStr+="/"+IJ.d2s(simulPars[0][3],7)+"/"+IJ.d2s(simulPars[0][4],7)+"/"+IJ.d2s(simulPars[0][5],7)+"}";
							dbgStr+=" {"+IJ.d2s(simulPars[1][0],5)+"/"+IJ.d2s(simulPars[1][1],5)+"/"+IJ.d2s(simulPars[1][2],5);
							if (simulPars[1].length>3) dbgStr+="/"+IJ.d2s(simulPars[1][3],7)+"/"+IJ.d2s(simulPars[1][4],7)+"/"+IJ.d2s(simulPars[1][5],7)+"}";
							System.out.println(dbgStr);
							if (nowDebugCell && (thisDebug>3)) {
								double [] XY={PATTERN_GRID[iUV[1]][iUV[0]][0][0]-32.0,PATTERN_GRID[iUV[1]][iUV[0]][0][1]-32.0};
								System.out.println("iUV[0]="+iUV[0]+"iUV[1]="+iUV[1]);
								System.out.println("CC : "+IJ.d2s(PATTERN_GRID[iUV[1]  ][iUV[0]  ][0][0]-XY[0],3)+"/"+IJ.d2s(PATTERN_GRID[iUV[1]  ][iUV[0]  ][0][1]-XY[1],3));
								System.out.println("TL : "+IJ.d2s(PATTERN_GRID[iUV[1]-1][iUV[0]-1][0][0]-XY[0],3)+"/"+IJ.d2s(PATTERN_GRID[iUV[1]-1][iUV[0]-1][0][1]-XY[1],3)); // sometimes throws
								System.out.println("TC : "+IJ.d2s(PATTERN_GRID[iUV[1]-1][iUV[0]  ][0][0]-XY[0],3)+"/"+IJ.d2s(PATTERN_GRID[iUV[1]-1][iUV[0]  ][0][1]-XY[1],3));
								System.out.println("TR : "+IJ.d2s(PATTERN_GRID[iUV[1]-1][iUV[0]+1][0][0]-XY[0],3)+"/"+IJ.d2s(PATTERN_GRID[iUV[1]-1][iUV[0]+1][0][1]-XY[1],3));
								System.out.println("CR : "+IJ.d2s(PATTERN_GRID[iUV[1]  ][iUV[0]+1][0][0]-XY[0],3)+"/"+IJ.d2s(PATTERN_GRID[iUV[1]  ][iUV[0]+1][0][1]-XY[1],3));
								System.out.println("BR : "+IJ.d2s(PATTERN_GRID[iUV[1]+1][iUV[0]+1][0][0]-XY[0],3)+"/"+IJ.d2s(PATTERN_GRID[iUV[1]+1][iUV[0]+1][0][1]-XY[1],3));
								System.out.println("BC : "+IJ.d2s(PATTERN_GRID[iUV[1]+1][iUV[0]  ][0][0]-XY[0],3)+"/"+IJ.d2s(PATTERN_GRID[iUV[1]+1][iUV[0]  ][0][1]-XY[1],3));
								System.out.println("BL : "+IJ.d2s(PATTERN_GRID[iUV[1]+1][iUV[0]-1][0][0]-XY[0],3)+"/"+IJ.d2s(PATTERN_GRID[iUV[1]+1][iUV[0]-1][0][1]-XY[1],3));
								System.out.println("CL : "+IJ.d2s(PATTERN_GRID[iUV[1]  ][iUV[0]-1][0][0]-XY[0],3)+"/"+IJ.d2s(PATTERN_GRID[iUV[1]  ][iUV[0]-1][0][1]-XY[1],3));

								System.out.println("CC : "+IJ.d2s(PATTERN_GRID[iUV[1]  ][iUV[0]  ][0][0],3)+"/"+IJ.d2s(PATTERN_GRID[iUV[1]  ][iUV[0]  ][0][1],3));
								System.out.println("TL : "+IJ.d2s(PATTERN_GRID[iUV[1]-1][iUV[0]-1][0][0],3)+"/"+IJ.d2s(PATTERN_GRID[iUV[1]-1][iUV[0]-1][0][1],3));
								System.out.println("TC : "+IJ.d2s(PATTERN_GRID[iUV[1]-1][iUV[0]  ][0][0],3)+"/"+IJ.d2s(PATTERN_GRID[iUV[1]-1][iUV[0]  ][0][1],3));
								System.out.println("TR : "+IJ.d2s(PATTERN_GRID[iUV[1]-1][iUV[0]+1][0][0],3)+"/"+IJ.d2s(PATTERN_GRID[iUV[1]-1][iUV[0]+1][0][1],3));
								System.out.println("CR : "+IJ.d2s(PATTERN_GRID[iUV[1]  ][iUV[0]+1][0][0],3)+"/"+IJ.d2s(PATTERN_GRID[iUV[1]  ][iUV[0]+1][0][1],3));
								System.out.println("BR : "+IJ.d2s(PATTERN_GRID[iUV[1]+1][iUV[0]+1][0][0],3)+"/"+IJ.d2s(PATTERN_GRID[iUV[1]+1][iUV[0]+1][0][1],3));
								System.out.println("BC : "+IJ.d2s(PATTERN_GRID[iUV[1]+1][iUV[0]  ][0][0],3)+"/"+IJ.d2s(PATTERN_GRID[iUV[1]+1][iUV[0]  ][0][1],3));
								System.out.println("BL : "+IJ.d2s(PATTERN_GRID[iUV[1]+1][iUV[0]-1][0][0],3)+"/"+IJ.d2s(PATTERN_GRID[iUV[1]+1][iUV[0]-1][0][1],3));
								System.out.println("CL : "+IJ.d2s(PATTERN_GRID[iUV[1]  ][iUV[0]-1][0][0],3)+"/"+IJ.d2s(PATTERN_GRID[iUV[1]  ][iUV[0]-1][0][1],3));
							}
						}

						//							if (nowDebugCell)correctedPatternCrossLocationAverage4(

						double [] centerXY=correctedPatternCrossLocation(
								lwirReaderParameters, // LwirReaderParameters lwirReaderParameters, // null is OK
								patternGrid[iUV[1]][iUV[0]][0], // initial coordinates of the pattern cross point
								patternGrid[iUV[1]][iUV[0]][1][0],
								patternGrid[iUV[1]][iUV[0]][1][1],
								patternGrid[iUV[1]][iUV[0]][2][0],
								patternGrid[iUV[1]][iUV[0]][2][1],
								simulPars,
								imp,      // image data (Bayer mosaic)
								distortionParameters, //
								patternDetectParameters,
								matchSimulatedPatternCorr, // correlationSize
								thisSimulParameters,
								equalizeGreens,
								windowFunctionCorr,
								windowFunctionCorr2,
								windowFunctionCorr4,
								simulationPattern,
								((iUV[0]^iUV[1])&1)!=0, // if true - invert pattern
								fht_instance,
								distortionParameters.fastCorrelationOnFinalPass, //
								locsNeib,
								thisDebug, //thisDebug
								null);


						if (centerXY!=null){
							if (thisDebug>2) System.out.println("==>iUV={"+iUV[0]+",  "+iUV[1]+
									"}. "+patternGrid[iUV[1]][iUV[0]][0][0]+" / "+patternGrid[iUV[1]][iUV[0]][0][1]+" -> "+centerXY[0]+" / "+centerXY[1]);
							if (refineInPlace)	setPatternGridCell(
									patternGrid,
									iUV,
									centerXY,
									null, //  double [] wv1,
									null); //  double [] wv2);
							else newGrid[iUV[1]][iUV[0]]=centerXY.clone();

						} else {
							if (debug_level>0){
								System.out.println("refineDistortionCorrelation(): failed to refine grid for U="+iUV[0]+" V="+iUV[1]+
										" X="+patternGrid[iUV[1]][iUV[0]][0][0]+" Y="+patternGrid[iUV[1]][iUV[0]][0][1]);
							}

						}
						final int numFinished=cellNumDoneAtomic.getAndIncrement();
						SwingUtilities.invokeLater(new Runnable() {
							@Override
							public void run() {
								// Here, we can safely update the GUI
								// because we'll be called from the
								// event dispatch thread
								IJ.showProgress(numFinished,width*height-1);
							}
						});
					}
				}
			};
		}
		startAndJoin(threads);
		IJ.showProgress(1.0); // turn off
		double maxActualCorr=0.0;
		double dx,dy,dist;
		if (!refineInPlace) {
			//maxCorr
			if (maxCorr>0.0) {
				// make sure there are no new undefined cells that were initially defined and no correction more than the limit
				int numUndefined=0,numFar=0;
				//					double maxCorr2=maxCorr*maxCorr;
				for (iUV[1]=0;iUV[1]<height;iUV[1]++) for (iUV[0]=0;iUV[0]<width;iUV[0]++) if (isCellDefined(patternGrid, iUV)){
					if (newGrid[iUV[1]][iUV[0]]==null){
						numUndefined++;
					} else {
						dx=newGrid[iUV[1]][iUV[0]][0] - patternGrid[iUV[1]][iUV[0]][0][0];
						dy=newGrid[iUV[1]][iUV[0]][1] - patternGrid[iUV[1]][iUV[0]][0][1];
						dist=Math.sqrt(dx*dx+dy*dy);
						if (dist>maxActualCorr) maxActualCorr=dist;
						if (dist>maxCorr) {
							numFar++;
							newGrid[iUV[1]][iUV[0]]=null;
						}
					}
				}
				if ((numUndefined>0) || (numFar>0)) {
					if (debug_level>0){
						System.out.println("refineDistortionCorrelation(): failed, number of undefined cells="+numUndefined+", number of too far cells="+numFar+" maxActualCorr="+maxActualCorr );
					}
					if (numUndefined>0) return -numUndefined; // negative - some cells undefined, no info about maximal correction returned
					return maxActualCorr; // no correction performed
				}

			} else {
				// only calculate maximal distance
				for (iUV[1]=0;iUV[1]<height;iUV[1]++) for (iUV[0]=0;iUV[0]<width;iUV[0]++) if (isCellDefined(patternGrid, iUV) && (newGrid[iUV[1]][iUV[0]]!=null)){
					dx=newGrid[iUV[1]][iUV[0]][0] - patternGrid[iUV[1]][iUV[0]][0][0];
					dy=newGrid[iUV[1]][iUV[0]][1] - patternGrid[iUV[1]][iUV[0]][0][1];
					dist=Math.sqrt(dx*dx+dy*dy);
					if (dist>maxActualCorr) maxActualCorr=dist;
				}

			}
			// Copy new values for the grid cells
			for (iUV[1]=0;iUV[1]<height;iUV[1]++) for (iUV[0]=0;iUV[0]<width;iUV[0]++) if (newGrid[iUV[1]][iUV[0]]!=null){
				setPatternGridCell(
						patternGrid,
						iUV,
						newGrid[iUV[1]][iUV[0]],
						null, //  double [] wv1,
						null); //  double [] wv2);

			}
			// correction is only calculated for simultaneous update (not for in-place)
			if (debug_level>1){
				System.out.println("refineDistortionCorrelation(): maximal correction="+ maxActualCorr+" pixels");
			}
		}
		debugLevel=was_debug_level;
		return maxActualCorr;
	}


	public class MinMaxSync {
		private double min;
		private double max;
		private boolean defined;
		public MinMaxSync(){
			defined=false;
			min=Double.NaN;;
			max=Double.NaN;;
		}
		public void reset(){
			defined=false;
		}
		public synchronized void minMax(double d){
			if (!defined){
				min=d;
				max=d;
				defined=true;
			} else {
				if (d>max) max=d;
				else if (d<min) min=d;
			}
		}
		public double getMin() {return min;}
		public double getMax() {return max;}
		public boolean isDefined() {return defined;}
	}
	/* ================================================================*/
	/*
        public boolean flatFieldCorrection=true; // compensate grid uneven intensity (vignetting, illumination)
        public double flatFieldExtarpolate=1.0;  // extrapolate flat field intensity map (relative to the average grid period)
        public double flatFieldBlur=1.0;   // blur the intensity map (relative to the average grid period)

	 */
	public ImagePlus equalizeGridIntensity(
			ImagePlus imp,
			double [][][][] patternGrid,
			DistortionParameters distortionParameters, //
			boolean equalizeGreens,
			int debugLevel,
			boolean updateStatus,
			int threadsMax
			){
		int dbgThreshold=1;
		double [][] gridIntensity=calcGridIntensity(
				4, //bayerComponent
				distortionParameters.correlationSize, // size
				distortionParameters, //
				equalizeGreens,
				imp, // image to process
				patternGrid,
				threadsMax);// debug level used inside loops
		if (debugLevel>(dbgThreshold+2)){
			double [] testGI=new double [gridIntensity.length*gridIntensity[0].length];
			int index=0;
			for (int v=0;v<gridIntensity.length;v++) for (int u=0;u<gridIntensity[0].length;u++)testGI[index++]=gridIntensity[v][u];
			this.SDFA_INSTANCE.showArrays(testGI, gridIntensity[0].length, gridIntensity.length, imp.getTitle()+"-GI");
		}
		double [] fffg=calcFlatFieldForGrid(
				gridIntensity,
				patternGrid,
				imp.getWidth(),
				imp.getHeight());

		double averageGridPeriod=averageGridPeriod( patternGrid);

		int preShrink= (int) (averageGridPeriod * distortionParameters.flatFieldShrink);
		int expand=    (int) (averageGridPeriod * distortionParameters.flatFieldExpand);
		double extrapolateSigma=averageGridPeriod * distortionParameters.flatFieldSigmaRadius;
		double extrapolateKSigma=distortionParameters.flatFieldExtraRadius;
		if (debugLevel>=(dbgThreshold+2)){
			this.SDFA_INSTANCE.showArrays(fffg.clone(), imp.getWidth(), imp.getHeight(), imp.getTitle()+"-fffg");
		}
		extrapolatePatternFlatFieldCorrection(
				fffg, //fieldXY,
				imp.getWidth(),
				preShrink,
				expand,
				extrapolateSigma,
				extrapolateKSigma,
				threadsMax,     //   100; // testing multi-threading, limit maximal number of threads
				updateStatus);

		if (debugLevel>(dbgThreshold+2)){
			this.SDFA_INSTANCE.showArrays(fffg.clone(), imp.getWidth(), imp.getHeight(), imp.getTitle()+"-extrapolated");
		}
		if (distortionParameters.flatFieldBlur>0.0) {
			DoubleGaussianBlur gb=new DoubleGaussianBlur();
			gb.blurDouble(
					fffg,
					imp.getWidth(),
					imp.getHeight(),
					distortionParameters.flatFieldBlur*averageGridPeriod,
					distortionParameters.flatFieldBlur*averageGridPeriod,
					0.01);
		}

		double max=0.0;
		for (int i=0;i<fffg.length;i++) if (max<fffg[i]) max =fffg[i];
		double k=1.0/max;

		for (int i=0;i<fffg.length;i++) {
			fffg[i]*=k;
			if (fffg[i]<distortionParameters.flatFieldMin)fffg[i]=0.0;
		}

		if (debugLevel>1) System.out.println("averageGridPeriod="+averageGridPeriod);

		if (debugLevel>(dbgThreshold+1)){
			this.SDFA_INSTANCE.showArrays(fffg, imp.getWidth(), imp.getHeight(), imp.getTitle()+"-blured");
		}

		this.flatFieldForGrid=fffg;

		ImagePlus imp_eq=  applyFlatField (imp, fffg);
		if (debugLevel>dbgThreshold) imp_eq.show();
		return imp_eq;
	}

	public ImagePlus applyFlatField (ImagePlus imp){
		if (this.PATTERN_GRID==null) return imp;
		if ((getImageHeight()!=imp.getHeight()) || (getImageWidth()!=imp.getWidth())){
			String msg="applyFlatField (): Supplied image does not match in dimensions "+imp.getWidth()+"x"+imp.getHeight()+
					" the one for wich grid was calculated ("+getImageWidth()+"x"+getImageHeight()+")";
			//			   IJ.showMessage("Error",msg);
			throw new IllegalArgumentException (msg);
		}
		return applyFlatField (imp, this.flatFieldForGrid);

	}
	public ImagePlus applyFlatField (ImagePlus imp, double [] ff){
		if (ff==null) return imp; // nothing to apply
		float [] pixels=(float []) imp.getProcessor().getPixels();

		if (pixels.length!=ff.length){
			String msg="Supplied image does not match in dimensions "+(pixels.length)+
					" the one for wich grid was calculated ("+(ff.length)+")";
			//			   IJ.showMessage("Error",msg);
			throw new IllegalArgumentException (msg);
		}



		float [] eqPixels=new float [pixels.length];
		for (int i=0;i<pixels.length;i++) if (ff[i]>0) eqPixels[i]=(float) (pixels[i]/ff[i]);
		else eqPixels[i]=0.0f;

		ImageProcessor ip=new FloatProcessor(imp.getWidth(), imp.getHeight());
		ip.setPixels(eqPixels);
		ip.resetMinAndMax();
		ImagePlus imp_eq=  new ImagePlus(imp.getTitle()+"-flat", ip);
		return imp_eq;
	}

	public double [][][] calcGridIntensities(
			final DistortionParameters distortionParameters, //
			final boolean equalizeGreens,
			final ImagePlus imp, // image to process
			final int threadsMax){
		double dSize=averageGridPeriod(this.PATTERN_GRID)*distortionParameters.averagingAreaScale;
		int size=  ~1 & ((int) Math.round(dSize)); // should be even
		int size4= ~1 & ((int) Math.round(dSize*Math.sqrt(2.0))); // larger when using diagonal greens (component 4)
		int [] bayerIndices={-1,1,4,2}; // -1 - contrast, 1 - R, 4 - G, 2 - B
		this.gridContrastBrightness=new double [bayerIndices.length][][];
		if (this.debugLevel>1) System.out.println("Calculating grid intensities, average period="+IJ.d2s(averageGridPeriod(this.PATTERN_GRID),2)+
				" pixels, using square sample "+size+"x"+size+" for all colors ,but (diagonal) greens, for greens - "+size4+"x"+size4);
		for (int i=0;i<bayerIndices.length;i++){
			this.gridContrastBrightness[i]=calcGridIntensity (
					bayerIndices[i], //final int bayerComponent,
					((bayerIndices[i]==4)?size4:size), //final int size,
					distortionParameters, //final DistortionParameters distortionParameters, //
					equalizeGreens,
					imp, // image to process
					this.PATTERN_GRID,
					threadsMax);
		}
		return this.gridContrastBrightness;
	}


	public double [][] calcGridIntensity (
			final int bayerComponent,
			final int size,
			final DistortionParameters distortionParameters, //
			final boolean equalizeGreens,
			final ImagePlus imp, // image to process
			final double [][][][] patternGrid,
			final int threadsMax){// debug level used inside loops
		final double [][] gridIntensity=new double[patternGrid.length][patternGrid[0].length];
		for (int i=0;i<gridIntensity.length;i++) for (int j=0;j<gridIntensity[0].length;j++)gridIntensity[i][j]=(bayerComponent>=0)?-1.0:0.0; // undefined
		MatchSimulatedPattern matchSimulatedPatternCorr=new MatchSimulatedPattern(distortionParameters.correlationSize);
		matchSimulatedPatternCorr.debugLevel=debugLevel;
		final double [] windowFunctionCorr= initWindowFunction(
				size, // distortionParameters.correlationSize,
				distortionParameters.correlationGaussWidth,
				distortionParameters.zeros);
		final int width=patternGrid[0].length;
		final int height=patternGrid.length;
		final Thread[] threads = newThreadArray(threadsMax);
		if (debugLevel>1)  System.out.println("Calculating average intensity at grid nodes, width= "+width+", height= "+height+ ", bayerComponent="+bayerComponent+", size="+size);
		int [] iUV=  new int [2];
		int i;
		for (i=0;i<(width*height);i++) {
			iUV[0]=i % width;
			iUV[1]=i / width;
			if (isCellDefined(patternGrid,iUV)) break;
		}
		final int startCell=i;
		final AtomicInteger cellNum = new AtomicInteger(startCell);
		final double [][][] newGrid=new double [height][width][];
		for (int v=0;v<height;v++) for (int u=0;u<width;u++) newGrid[v][u]=null;
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					int [] iUV=  new int [2];
					for (int ncell=cellNum.getAndIncrement(); ncell<(width*height);ncell=cellNum.getAndIncrement()){
						iUV[0]=ncell % width;
						iUV[1]=ncell / width;
						if (!isCellDefined(patternGrid,iUV)) continue;
						if (bayerComponent>=0){
							Rectangle centerCross=correlationSelection(
									patternGrid[iUV[1]][iUV[0]][0], // initial coordinates of the pattern cross point
									size); // distortionParameters.correlationSize);
							double[][] input_bayer=splitBayer (imp,centerCross,equalizeGreens);
							double sum=0.0, sumW=0.0;
							for (int i=0;i<input_bayer[bayerComponent].length;i++){
								sum+= input_bayer[bayerComponent][i]*windowFunctionCorr[i];
								sumW+=                  windowFunctionCorr[i];
							}
							gridIntensity[iUV[1]][iUV[0]]=sum/sumW;
						} else {
							// trying alternative
							//							   double [][][][] patternGrid_same=patternGrid;
							gridIntensity[iUV[1]][iUV[0]]=Double.NaN;
							if (isCellDefined(patternGrid,iUV[0],iUV[1])) {
								double [][] patternCell=patternGrid[iUV[1]][iUV[0]];
								if (patternCell[0].length>2) gridIntensity[iUV[1]][iUV[0]]=patternCell[0][2];
							}
							/*
							   gridIntensity[iUV[1]][iUV[0]]=localGridContrast(
									   imp,
									   equalizeGreens,
									   patternGrid,
									   iUV[0],
									   iUV[1]);
							 */
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return gridIntensity;
	}
	/* ======================================================================== */

	public double localGridContrast(
			ImagePlus imp,
			boolean equalizeGreens,
			final double [][][][] patternGrid,
			int u,
			int v){
		if (!isCellDefined(patternGrid,u,v))   return 0.0;
		if (!isCellDefined(patternGrid,u+1,v)) return 0.0;
		if (!isCellDefined(patternGrid,u,v+1)) return 0.0;
		if (!isCellDefined(patternGrid,u-1,v)) return 0.0;
		if (!isCellDefined(patternGrid,u,v-1)) return 0.0;
		double [][] deltas={
				{0.25*(patternGrid[v+1][u][0][0]-patternGrid[v-1][u][0][0]),
					0.25*(patternGrid[v+1][u][0][1]-patternGrid[v-1][u][0][1])},
				{0.25*(patternGrid[v][u+1][0][0]-patternGrid[v][u-1][0][0]),
						0.25*(patternGrid[v][u+1][0][1]-patternGrid[v][u-1][0][1])} };
		double delta=Math.sqrt(0.5*(deltas[0][0]*deltas[0][0]+
				deltas[0][1]*deltas[0][1]+
				deltas[1][0]*deltas[1][0]+
				deltas[1][1]*deltas[1][1]));
		int range=(int) Math.round(0.25*delta); // center of the white/black;
		int [][] iDeltas={
				{(int) Math.round(deltas[0][0]),
					(int) Math.round(deltas[0][1])},
				{(int) Math.round(deltas[1][0]),
						(int) Math.round(deltas[1][1])}};

		int [][] centersOnBayer4=
			{
					{ iDeltas[0][0], iDeltas[0][1]},
					{-iDeltas[0][0],-iDeltas[0][1]},
					{ iDeltas[1][0], iDeltas[1][1]},
					{-iDeltas[1][0],-iDeltas[1][1]} };

		double diff=0.0;
		double sum=0.0;
		int maxDxy=0;
		for (int n=0;n<4;n++) for (int i=0;i<2;i++) if (centersOnBayer4[n][i]>maxDxy) maxDxy=centersOnBayer4[n][i];
		int size=2*(maxDxy+range+1); // this will include all needed pixels
		//		   size+=2;
		int hSize=size/2;
		boolean debug=false; //((u==30) && (v==30));
		Rectangle centerCross=correlationSelection(
				patternGrid[v][u][0], // initial coordinates of the pattern cross point
				size); // distortionParameters.correlationSize);

		//		   Rectangle thisSel=new Rectangle(centerCross.x,centerCross.y,2*size,2*size); // "2" - sensor pixels, befor split to components
		double[][] input_bayer=splitBayer (imp,centerCross,equalizeGreens);
		if (debug) this.SDFA_INSTANCE.showArrays(input_bayer, size,size, true, imp.getTitle()+"-bayer");

		double [] bayer4=input_bayer[4];
		for (int dv=-range;dv<=range;dv++) for (int du=-range;du<=range;du++) {
			int [] indices= new int[4];
			for (int n=0;n<4;n++) indices[n]=size*(centersOnBayer4[n][1]+dv+hSize)+(centersOnBayer4[n][0]+du+hSize);
			if ((indices[0]>bayer4.length) || (indices[1]>bayer4.length) || (indices[2]>bayer4.length) || (indices[3]>bayer4.length)
					|| (indices[0]<0) || (indices[1]<0) || (indices[2]<0) || (indices[3]<0)||
					debug){
				System.out.println("centersOnBayer4[0]={"+centersOnBayer4[0][0]+", "+centersOnBayer4[0][1]+"}");
				System.out.println("centersOnBayer4[1]={"+centersOnBayer4[1][0]+", "+centersOnBayer4[1][1]+"}");
				System.out.println("centersOnBayer4[2]={"+centersOnBayer4[2][0]+", "+centersOnBayer4[2][1]+"}");
				System.out.println("centersOnBayer4[3]={"+centersOnBayer4[3][0]+", "+centersOnBayer4[3][1]+"}");
				System.out.println("range="+range);
				System.out.println("dv="+dv+" du="+du);
				System.out.println("maxDxy="+maxDxy+" size="+size+" hSize="+hSize);
				System.out.println("indices=={"+indices[0]+", "+indices[1]+", "+indices[2]+", "+indices[3]+"}, bayer4.length="+bayer4.length);
			}


			sum+= bayer4[indices[0]]+bayer4[indices[1]]+bayer4[indices[2]]+bayer4[indices[3]];
			diff+=bayer4[indices[0]]+bayer4[indices[1]]-bayer4[indices[2]]-bayer4[indices[3]];

		}
		if (((u ^ v) & 1)!=0) diff=-diff;
		if (sum==0.0) return 0.0;
		return diff/sum;
	}

	/* ======================================================================== */
	public double averageGridPeriod(
			double [][][][] patternGrid){
		int n=0;
		double sum=0.0;
		int [] iUV=new int[2];
		int [] iUV1=new int[2];
		int [][]dirs={{0,1},{1,0}};
		double dx,dy;
		for (iUV[1]=0;iUV[1]<patternGrid.length-1;iUV[1]++) for (iUV[0]=0;iUV[0]<patternGrid[0].length-1;iUV[0]++) if (isCellDefined(patternGrid,iUV)){
			for (int dir=0;dir<dirs.length;dir++){
				iUV1[0]=iUV[0]+dirs[dir][0];
				iUV1[1]=iUV[1]+dirs[dir][1];
				if (isCellDefined(patternGrid,iUV1)){
					//					   dx=patternGrid[iUV1[1]][iUV1[0]][0][0]-patternGrid[iUV1[1]][iUV[0]][0][0]; // old bug, skewed period!
					//					   dy=patternGrid[iUV1[1]][iUV1[0]][0][1]-patternGrid[iUV1[1]][iUV[0]][0][1]; // old bug, skewed period!
					dx=patternGrid[iUV1[1]][iUV1[0]][0][0]-patternGrid[iUV[1]][iUV[0]][0][0];
					dy=patternGrid[iUV1[1]][iUV1[0]][0][1]-patternGrid[iUV[1]][iUV[0]][0][1];
					sum+=dx*dx+dy*dy;
					n++;
				}
			}
		}
		if (n>0) sum/=n;
		return Math.sqrt (sum);
	}
	/* ======================================================================== */
	public double [] averageGridPeriods( // min,max for u,v
			double [][][][] patternGrid){
		double [] result={Double.NaN,Double.NaN};
		//		   int n=0;
		double [] sum={0.0,0.0};
		int [] numSamples={0,0};
		int [] iUV=new int[2];
		int [] iUV1=new int[2];
		int [][]dirs={{0,1},{1,0}};
		double dx,dy;
		for (iUV[1]=0;iUV[1]<patternGrid.length-1;iUV[1]++) for (iUV[0]=0;iUV[0]<patternGrid[0].length-1;iUV[0]++) if (isCellDefined(patternGrid,iUV)){
			for (int dir=0;dir<dirs.length;dir++){
				iUV1[0]=iUV[0]+dirs[dir][0];
				iUV1[1]=iUV[1]+dirs[dir][1];
				if (isCellDefined(patternGrid,iUV1)){
					dx=patternGrid[iUV1[1]][iUV1[0]][0][0]-patternGrid[iUV[1]][iUV[0]][0][0];
					dy=patternGrid[iUV1[1]][iUV1[0]][0][1]-patternGrid[iUV[1]][iUV[0]][0][1];
					sum[dir]+=dx*dx+dy*dy;
					numSamples[dir]++;
				}
			}
		}
		for (int dir=0;dir<dirs.length;dir++){
			if (numSamples[dir]>0) result[dir]=Math.sqrt (sum[dir]/numSamples[dir]);
		}
		if (result[0]>result[1]){
			double tmp=result[0];
			result[0]=result[1];
			result[1]=tmp;
		}
		return result;
	}

	/* ======================================================================== */
	public double [] calcFlatFieldForGrid(
			double [][] gridIntensity,
			double [][][][] patternGrid,
			int sWidth,
			int sHeight){
		int width=patternGrid[0].length;
		int height=patternGrid.length;
		int [][] uvInc={{0,0},{1,0},{0,1},{1,1}}; // four corners as u,v pair
		int [][] cycles={ // counter-clockwise corners bounding the area  (only orthogonal sides?)
				{1,0,2},
				{2,3,1},
				{0,2,3},
				{3,1,0}};
		double [] fffg=   new double [sWidth*sHeight];
		int    [] fffgNum=new int [sWidth*sHeight];
		for (int i=0;i<fffg.length;i++){
			fffg[i]=0.0;
			fffgNum[i]=0;
		}
		int [] iUV=new int[2];
		for (int v=0;v<(height-1); v++) for (int u=0; u<(width-1);u++){
			double [][] cornerXY =new double[4][];
			for (int i=0;i<uvInc.length;i++){
				iUV[0]=u+uvInc[i][0];
				iUV[1]=v+uvInc[i][1];

				if (isCellDefined(patternGrid,iUV)){
					cornerXY[i]=new double[3];
					cornerXY[i][0]=patternGrid[iUV[1]][iUV[0]][0][0];
					cornerXY[i][1]=patternGrid[iUV[1]][iUV[0]][0][1];
					cornerXY[i][2]=gridIntensity[iUV[1]][iUV[0]];
				} else cornerXY[i]=null;
			}
			boolean [] cycleFits=new boolean[cycles.length];
			for (int i=0;i<cycles.length;i++){
				cycleFits[i]=true;
				for (int j=0;j<cycles[i].length;j++) if (cornerXY[cycles[i][j]]==null) {
					cycleFits[i]=false;
					break;
				}
			}
			if (cycleFits[0]&&cycleFits[1]){ // remove overlaps
				cycleFits[2]=false;
				cycleFits[3]=false;
			}
			boolean minMaxUndefined=true;
			double minX=0,maxX=0,minY=0,maxY=0;
			// find bounding rectangle;
			for (int nCycle=0;nCycle<cycles.length;nCycle++) if (cycleFits[nCycle]){
				int [] cycle=cycles[nCycle];
				for (int corner=0; corner<cycle.length;corner++){
					if (minMaxUndefined || (minX>cornerXY[cycle[corner]][0])) minX=cornerXY[cycle[corner]][0];
					if (minMaxUndefined || (maxX<cornerXY[cycle[corner]][0])) maxX=cornerXY[cycle[corner]][0];
					if (minMaxUndefined || (minY>cornerXY[cycle[corner]][1])) minY=cornerXY[cycle[corner]][1];
					if (minMaxUndefined || (maxY<cornerXY[cycle[corner]][1])) maxY=cornerXY[cycle[corner]][1];
					minMaxUndefined=false;
				}
			}
			int iMinX=(int) Math.floor(minX);
			int iMinY=(int) Math.floor(minY);
			int iMaxX=(int) Math.ceil(maxX);
			int iMaxY=(int) Math.ceil(maxY);
			if (iMinX<0) iMinX=0;
			if (iMinY<0) iMinY=0;
			if (iMaxX>=sWidth)  iMaxX=sWidth-1;
			if (iMaxY>=sHeight) iMaxY=sHeight-1;
			double [] originXY=new double [2];
			double [] endXY=new double [2];
			for (int idY=iMinY; idY<=iMaxY;idY++){
				double pY=idY; // in sensor pixels
				for (int idX=iMinX; idX<=iMaxX;idX++){
					double pX=idX; // in sensor pixels
					// scan allowed triangles, usually 2
					for (int nCycle=0;nCycle<cycles.length;nCycle++) if (cycleFits[nCycle]){
						int [] cycle=cycles[nCycle];
						// is this point inside?
						boolean inside=true;
						for (int nEdge=0;nEdge<cycle.length;nEdge++){
							int nextNEdge=(nEdge==(cycle.length-1))?0:(nEdge+1);

							originXY[0]=patternGrid[v+uvInc[cycle[nEdge]][1]][u+uvInc[cycle[nEdge]][0]][0][0];
							originXY[1]=patternGrid[v+uvInc[cycle[nEdge]][1]][u+uvInc[cycle[nEdge]][0]][0][1];
							endXY[0]=   patternGrid[v+uvInc[cycle[nextNEdge]][1]][u+uvInc[cycle[nextNEdge]][0]][0][0];
							endXY[1]=   patternGrid[v+uvInc[cycle[nextNEdge]][1]][u+uvInc[cycle[nextNEdge]][0]][0][1];
							if (((pX-originXY[0])*(endXY[1]-originXY[1]) - (pY-originXY[1])*(endXY[0]-originXY[0]))<0.0){
								inside=false;
								break;
							}
						}
						if (!inside) continue; // point is outside of the interpolation area, try next triangle (if any)
						/* interpolate:
							1. taking cycles[0] as origin and two (non co-linear) edge vectors - V1:from 0 to 1 and V2 from 1 to 2
							    find a1 and a2  so that vector V  (from 0  to pXY) = a1*V1+ a2*V2
							2. if F0 is the value of the interpolated function at cycles[0], F1 and F2 - at cycles[1] and cycles2
							   then F=F0+(F1-F0)*a1 +(F2-F1)*a2
						 */

						double [] XY0={patternGrid[v+uvInc[cycle[0]][1]][u+uvInc[cycle[0]][0]][0][0],patternGrid[v+uvInc[cycle[0]][1]][u+uvInc[cycle[0]][0]][0][1]};
						double [] XY1={patternGrid[v+uvInc[cycle[1]][1]][u+uvInc[cycle[1]][0]][0][0],patternGrid[v+uvInc[cycle[1]][1]][u+uvInc[cycle[1]][0]][0][1]};
						double [] XY2={patternGrid[v+uvInc[cycle[2]][1]][u+uvInc[cycle[2]][0]][0][0],patternGrid[v+uvInc[cycle[2]][1]][u+uvInc[cycle[2]][0]][0][1]};

						double [] V= {pX-XY0[0],pY-XY0[1]};
						double [][] M={
								{XY1[0]-XY0[0],XY2[0]-XY1[0]},
								{XY1[1]-XY0[1],XY2[1]-XY1[1]}};
						double det=M[0][0]*M[1][1]-M[1][0]*M[0][1];
						double [][] MInverse={
								{ M[1][1]/det,-M[0][1]/det},
								{-M[1][0]/det, M[0][0]/det}};
						double [] a12={
								MInverse[0][0]*V[0]+MInverse[0][1]*V[1],
								MInverse[1][0]*V[0]+MInverse[1][1]*V[1]};
						int pCorrIndex=idY*sWidth+idX;
						// some points may be accumulated multiple times - thisPCorr[3] will take care of this
						if (this.debugLevel>3) {
							System.out.println("XY0="+IJ.d2s(XY0[0],3)+":"+IJ.d2s(XY0[1],3));
							System.out.println("XY1="+IJ.d2s(XY1[0],3)+":"+IJ.d2s(XY1[1],3));
							System.out.println("XY2="+IJ.d2s(XY2[0],3)+":"+IJ.d2s(XY2[1],3));
							System.out.println("M00="+IJ.d2s(M[0][0],3)+" M01="+IJ.d2s(M[0][1],3));
							System.out.println("M10="+IJ.d2s(M[1][0],3)+" M11="+IJ.d2s(M[1][1],3));
							System.out.println("MInverse00="+IJ.d2s(MInverse[0][0],5)+" MInverse01="+IJ.d2s(MInverse[0][1],5));
							System.out.println("MInverse10="+IJ.d2s(MInverse[1][0],5)+" MInverse11="+IJ.d2s(MInverse[1][1],5));
							System.out.println("a12="+IJ.d2s(a12[0],3)+":"+IJ.d2s(a12[1],3));
							System.out.println("gridIntensity[v+uvInc[cycle[0]][1]][u+uvInc[cycle[0]][0]]="+
									IJ.d2s(gridIntensity[v+uvInc[cycle[0]][1]][u+uvInc[cycle[0]][0]],3));
							System.out.println("gridIntensity[v+uvInc[cycle[1]][1]][u+uvInc[cycle[1]][0]]="+
									IJ.d2s(gridIntensity[v+uvInc[cycle[1]][1]][u+uvInc[cycle[1]][0]],3));
							System.out.println("gridIntensity[v+uvInc[cycle[2]][1]][u+uvInc[cycle[2]][0]]="+
									IJ.d2s(gridIntensity[v+uvInc[cycle[2]][1]][u+uvInc[cycle[2]][0]],3));
						}
						double val=
								gridIntensity[v+uvInc[cycle[0]][1]][u+uvInc[cycle[0]][0]]+
								(gridIntensity[v+uvInc[cycle[1]][1]][u+uvInc[cycle[1]][0]]-gridIntensity[v+uvInc[cycle[0]][1]][u+uvInc[cycle[0]][0]])*a12[0]+
								(gridIntensity[v+uvInc[cycle[2]][1]][u+uvInc[cycle[2]][0]]-gridIntensity[v+uvInc[cycle[1]][1]][u+uvInc[cycle[1]][0]])*a12[1];
						if (this.debugLevel>3) {
							System.out.println("val="+IJ.d2s(val,3));
						}
						fffg[pCorrIndex]+=val;// error in /data/focus/grid3d/center/1317924548_967543-00.tiff OOB: 5019002
						fffgNum[pCorrIndex]+=1;
					}
				} // idX
				// use same order in calculations, make sure no gaps
			} // idY
		} // finished image
		for (int i=0;i<fffg.length;i++) if (fffgNum[i]>0){
			fffg[i]/=fffgNum[i];
		}
		return fffg;
	}
	/* ======================================================================== */

	/**
	 *  Extrapolates flat-field correction
	 * @param data [nPixels] data to extrapolate
	 * @param sWidth data width
	 * @param preShrink shrink the non-zero data by this number of pixels before extrapolating
	 * @param expand expand the (pre-shrank) data by up to this number of pixels
	 * @param sigma when fitting plane through new point use Gaussian weight function for the neighbors
	 *  (normalized to non-decimated points)
	 * @param ksigma Process pixels in a square with the side 2*sigma*ksigma
	 */
	// TODO: Use threads
	public boolean extrapolatePatternFlatFieldCorrection(
			final double [] data, //fieldXY,
			final int sWidth,
			final int preShrink,
			final int expand,
			final double sigma,
			final double ksigma,
			int    threadsMax,     //   100; // testing multi-threading, limit maximal number of threads
			boolean updateStatus){
		int dbgThreshold=1;
		final int length=data.length;
		final int sHeight=length/sWidth;
		// create mask
		final boolean [] fMask=new boolean[data.length];
		for (int i=0;i<length;i++) fMask[i]= data[i]>0.0;

		final int len= (int) Math.ceil(sigma*ksigma);

		final double [] gaussian=new double[len+1];
		double k=0.5/sigma/sigma;
		for (int i=0;i<=len;i++) gaussian[i]=Math.exp(-i*i*k);
		int [][] dirs={{-1,0},{1,0},{0,-1},{0,1}}; // order matters
		final List <Integer> extList=new ArrayList<Integer>(1000);
		Integer Index, Index2;
		extList.clear();
		// create initial wave
		if (this.debugLevel>2) System.out.println("extrapolatePatternFlatFieldCorrection() sWidth="+sWidth+" sHeight="+sHeight);

		for (int iy=0;iy<sHeight;iy++) for (int ix=0;ix<sWidth;ix++) {
			Index=iy*sWidth+ix;
			if (fMask[Index]) {
				int numNew=0;
				for (int dir=0;dir<dirs.length;dir++){
					int ix1=ix+dirs[dir][0];
					int iy1=iy+dirs[dir][1];
					if ((ix1>=0) && (iy1>=0) && (ix1<sWidth) && (iy1<sHeight)) {
						if (!fMask[iy1*sWidth+ix1]) numNew++;
					}
					if (numNew>0) extList.add(Index); // neighbor will have non-singular matrix
				}
			}
		}
		// now shrink
		// unmask current wave
		for (int i=extList.size()-1; i>=0;i--) fMask[extList.get(i)]=false;
		if (extList.size()==0) return false; // no points
		for (int nShrink=0;nShrink<preShrink;nShrink++){
			int size=extList.size();
			if (size==0) return false; // no points
			// wave step, unmasking
			for (int i=0; i<size;i++) {
				Index=extList.get(0);
				extList.remove(0);
				int iy=Index/sWidth;
				int ix=Index%sWidth;
				for (int dir=0;dir<dirs.length;dir++){
					int ix1=ix+dirs[dir][0];
					int iy1=iy+dirs[dir][1];
					if ((ix1>=0) && (iy1>=0) && (ix1<sWidth) && (iy1<sHeight)){
						Index=iy1*sWidth+ix1;
						if (fMask[Index]){
							extList.add(Index);
							fMask[Index]=false; // restore later?
						}
					}
				}
			}
		}
		// restore mask on the front
		for (int i=extList.size()-1; i>=0;i--) fMask[extList.get(i)]=true;

		if (this.debugLevel>dbgThreshold+1){
			for (int i=0;i<length;i++) if (!fMask[i]) data[i]=0.0;

			this.SDFA_INSTANCE.showArrays(data, sWidth,sHeight, "shrank");
		}


		// repeat with the wave until there is place to move, but not more than "expand" steps

		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger pixInWaveNum = new AtomicInteger();


		int [] dirs2=new int [2];
		for (int n=0; (n<expand) && (extList.size()>0); n++ ){
			if (updateStatus) IJ.showStatus("Expanding, step="+(n+1)+" (of "+expand+"), extList.size()="+extList.size());
			if (this.debugLevel>2) System.out.println("Expanding, step="+n+", extList.size()="+extList.size());
			// move wave front 1 pixel hor/vert
			for (int i=extList.size();i>0;i--){ // repeat current size times
				Index=extList.get(0);
				extList.remove(0);
				int iy=Index/sWidth;
				int ix=Index%sWidth;
				for (int dir=0;dir<dirs.length;dir++){
					int ix1=ix+dirs[dir][0];
					int iy1=iy+dirs[dir][1];
					if ((ix1>=0) && (iy1>=0) && (ix1<sWidth) && (iy1<sHeight)){
						Index=iy1*sWidth+ix1;
						if (!fMask[Index]){
							// verify it has neighbors in the perpendicular direction to dir
							dirs2[0]=(dir+2) & 3;
							dirs2[1]=dirs2[0] ^ 1;
							for (int dir2=0;dir2<dirs2.length;dir2++){
								int ix2=ix+dirs[dirs2[dir2]][0]; // from the old, not the new point!
								int iy2=iy+dirs[dirs2[dir2]][1];
								if ((ix2>=0) && (iy2>=0) && (ix2<sWidth) && (iy2<sHeight)){
									Index2=iy2*sWidth+ix2;
									if (fMask[Index2]){ // has orthogonal neighbor, OK to add
										extList.add(Index);
										fMask[Index]=true; // remove later
										break;
									}
								}
							}
						}
					}
				}
			}
			// now un-mask the pixels in new list new
			for (int i =0;i<extList.size();i++){
				Index=extList.get(i);
				fMask[Index]=false; // now mask is only set for known pixels
			}
			// Calculate values (extrapolate) for the pixels in the list
			/*
	Err = sum (W(x,y)*(f(x,y)-F0-Ax*(x-X0)-Ay*(y-Y0))^2)=
	sum (Wxy*(Fxy^2+F0^2+Ax^2*(x-X0)^2+Ay^2*(y-Y0)^2
	-2*Fxy*F0 -2*Fxy*Ax*(x-X0) - 2*Fxy*Ay*(y-Y0)
	+2*F0*Ax*(x-X0) + 2*F0*Ay*(y-Y0)
	+2*Ax*(x-X0)*Ay*(y-Y0))
	(1)0=dErr/dF0= 2*sum (Wxy*(F0-Fxy+Ax*(x-X0)+Ay(y-Y0)))
	(2)0=dErr/dAx= 2*sum (Wxy*(Ax*(x-X0)^2-Fxy*(x-X0) +F0*(x-X0)+Ay*(x-x0)*(y-Y0)))
	(3)0=dErr/dAy= 2*sum (Wxy*(Ay*(y-y0)^2-Fxy*(y-Y0) +F0*(y-Y0)+Ax*(x-x0)*(y-Y0)))

	S0 = sum(Wxy)
	SF=  sum(Wxy*Fxy)
	SX=  sum(Wxy*(x-X0)
	SY=  sum(Wxy*(y-Y0)
	SFX= sum(Wxy*Fxy*(x-X0)
	SFY= sum(Wxy*Fxy*(y-Y0)
	SX2= sum(Wxy*(x-X0)^2
	SY2= sum(Wxy*(y-Y0)^2
	SXY= sum(Wxy*(x-X0)*(y-Y0)

	(1) F0*S0 - SF + Ax*SX +Ay*Sy = 0
	(2) Ax*SX2-SFX+F0*SX+Ay*SXY = 0
	(3) Ay*Sy2 -SFY + F0*SY +Ax*SXY = 0

	(1) F0*S0  + Ax*SX +Ay*SY = SF
	(2) Ax*SX2+F0*SX+Ay*SXY = SFX
	(3) Ay*Sy2  + F0*SY +Ax*SXY = SFY


	   | F0 |
	V= | Ax |
	   | Ay |

	     | SF  |
	B =  | SFX |
	     | SFY |

	     | S0  SX   SY  |
	M =  | SX  SX2  SXY |
	     | SY  SXY  SY2 |

	M * V = B
			 */
			pixInWaveNum.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {

						//			   for (int i =0;i<extList.size();i++){

						for (int i=pixInWaveNum.getAndIncrement(); i<extList.size();i=pixInWaveNum.getAndIncrement()){

							Integer Indx=extList.get(i);
							int iy=Indx/sWidth;
							int ix=Indx%sWidth;
							double  S0= 0.0;
							double  SF= 0.0;
							double  SX= 0.0;
							double  SY= 0.0;
							double  SFX=0.0;
							double  SFY=0.0;
							double  SX2=0.0;
							double  SY2=0.0;
							double  SXY=0.0;
							int iYmin=iy-len; if (iYmin<0) iYmin=0;
							int iYmax=iy+len; if (iYmax>=sHeight) iYmax=sHeight-1;
							int iXmin=ix-len; if (iXmin<0) iXmin=0;
							int iXmax=ix+len; if (iXmax>=sWidth) iXmax=sWidth-1;
							for (int iy1=iYmin;iy1<=iYmax;iy1++) for (int ix1=iXmin;ix1<=iXmax;ix1++) {
								int ind=ix1+iy1*sWidth;
								if (fMask[ind]){
									double w=gaussian[(iy1>=iy)?(iy1-iy):(iy-iy1)]*gaussian[(ix1>=ix)?(ix1-ix):(ix-ix1)];
									S0+= w;
									SF+= w*data[ind];
									SX+= w*(ix1-ix);
									SY+= w*(iy1-iy);
									SFX+=w*data[ind]*(ix1-ix);
									SFY+=w*data[ind]*(iy1-iy);
									SX2+=w*(ix1-ix)*(ix1-ix);
									SY2+=w*(iy1-iy)*(iy1-iy);
									SXY+=w*(ix1-ix)*(iy1-iy);
								}

							}
							double [][] aB={{SF},{SFX},{SFY}};
							double [][] aM={
									{S0,SX, SY},
									{SX,SX2,SXY},
									{SY,SXY,SY2}
							};
							Matrix B=new Matrix(aB);
							Matrix M=new Matrix(aM);
							if (!(new LUDecomposition(M)).isNonsingular() && (S0!=0.0)){
								data[Indx]=SF/S0;
							} else {
								Matrix V=M.solve(B);  // sometimes singular
								data[Indx]=V.get(0,0);
							}

						}
					}
				};
			}
			startAndJoin(threads);


			// set mask again for the new calculated layer of pixels
			for (int i =0;i<extList.size();i++){
				Index=extList.get(i);
				fMask[Index]=true;
			}
			IJ.showProgress(n+1,expand);
		}
		IJ.showProgress(1.0);

		return true;
	}




	/* ======================================================================== */
	private double [][] calcNeibLocsWeights (
			DistortionParameters distortionParameters,
			boolean useNeib){
		double [][] locsNeib={{0.0,0.0,1.0}};
		if (!useNeib) return locsNeib;
		locsNeib= new double [9][3];
		double [][]dirs={{ 0.0, 0.0},
				{ 1.0, 0.0},
				{ 0.0, 1.0},
				{-1.0, 0.0},
				{ 0.0,-1.0},
				{ 1.0, 1.0},
				{ 1.0,-1.0},
				{-1.0, 1.0},
				{-1.0,-1.0}};
		int i;
		locsNeib[0][2]=1.0-distortionParameters.averageOrthoWeight-distortionParameters.averageOrthoWeight;
		for (i=0;i<4;i++) {
			locsNeib[i+1][2]=0.25*distortionParameters.averageOrthoWeight;
			locsNeib[i+5][2]=0.25*distortionParameters.averageDiagWeight;
		}
		double k=1.0;
		for (i=0;i<9;i++) {
			if (i>0) k=distortionParameters.averageOrthoDist;
			if (i>4) k=distortionParameters.averageDiagDist;
			locsNeib[i][0]=k*dirs[i][0];
			locsNeib[i][1]=k*dirs[i][1];
		}
		return locsNeib;
	}
	/* ======================================================================== */
	public void zeroNaNContrast(){
		for (double [][][]  row:this.PATTERN_GRID){
			for (double [][] node:row){
				if ((node!=null) && (node.length>0) && (node[0]!=null) && (node[0].length>2)){
					if (Double.isNaN(node[0][2])) node[0][2]=0.0;
				}
			}
		}
	}


	public double[][][][] recalculateWaveVectors (
			//			   double[][][][] patternGrid,
			final boolean updateStatus,
			final int debug_level){// debug level used inside loops
		//		   double[][][][] patternGrid=this.PATTERN_GRID;
		int i;
		int [] iuv=new int [2];
		final int [][] directionsUV8= {{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}}; // first 8 should be the same as in directionsUV
		final int [] directionsBits8= {1,4,1,4,2,8,2,8}; // should match directionsUV8
		int neibBits;
		int dir;
		int [] iUV=  new int [2];
		double [][][] neibors=new double [8][][]; // uv and xy vectors to 8 neibors (some may be null
		double [][] thisCell;
		double [][] otherCell;
		int was_debug_level=debugLevel;
		debugLevel=debug_level;
		if (debugLevel>1) System.out.println("Recalculating wave vectors from coordinates...");

		for (iuv[1]=0;iuv[1]<this.PATTERN_GRID.length;iuv[1]++) for (iuv[0]=0;iuv[0]<this.PATTERN_GRID[0].length;iuv[0]++) if (isCellValid(this.PATTERN_GRID,iuv)){
			if (debugLevel>2) System.out.println("<---= iuv= "+iuv[0]+",  "+iuv[1]);
			thisCell=this.PATTERN_GRID[iuv[1]][iuv[0]];
			neibBits=0;
			for (dir=0;dir<directionsUV8.length;dir++) {
				neibors[dir]=null;
				iUV[0]=iuv[0]+directionsUV8[dir][0];
				iUV[1]=iuv[1]+directionsUV8[dir][1];
				if ((iUV[0]<0) || (iUV[1]<0) ||
						(iUV[0]>=this.PATTERN_GRID[0].length) || (iUV[1]>=this.PATTERN_GRID.length)) continue; // don't fit into UV grid
				if (isCellValid(this.PATTERN_GRID,iUV)) {
					neibors[dir]= new double [2][2];
					otherCell=this.PATTERN_GRID[iUV[1]][iUV[0]];
					neibors[dir][0][0]=0.5*directionsUV8[dir][0];  // u
					neibors[dir][0][1]=0.5*directionsUV8[dir][1];  // v
					neibors[dir][1][0]=otherCell[0][0]-thisCell[0][0];  // x
					neibors[dir][1][1]=otherCell[0][1]-thisCell[0][1];  // y
					neibBits |= directionsBits8[dir];
				}
			}
			i=Integer.bitCount(neibBits);
			if (debugLevel>2) System.out.println("neibBits="+neibBits+", number of bits= "+i);
			if (i>1) {
				double[][] wv= waveVectorsFromNeib(neibors);
				setPatternGridCell(
						this.PATTERN_GRID,
						iuv,
						null, // XY already set
						wv[0],
						wv[1]);

				if (debugLevel>2) System.out.println("==+> number of bits:"+i+
						" vw00="+IJ.d2s(wv[0][0],5)+" vw01="+IJ.d2s(wv[0][1],5)+
						" vw10="+IJ.d2s(wv[1][0],5)+" vw11="+IJ.d2s(wv[1][1],5));
				//		    	wv= WaveVectorsFromNeib(neibors);
				//
				//	   vectors: [num_vector][0][0] - U
				//          [num_vector][0][1] - V
				//          [num_vector][1][0] - X
				//          [num_vector][1][1] - Y
				//          [num_vector] == null - skip
				//
			}
		}
		debugLevel=was_debug_level;
		return this.PATTERN_GRID;
	}

	/* ======================================================================== */
	private double [][] waveVectorsFromNeib(double [][][] vectors){
		/*
		 * 		   vectors: [num_vector][0][0] - U
		 *                  [num_vector][0][1] - V
		 *                  [num_vector][1][0] - X
		 *                  [num_vector][1][1] - Y
		 *                  [num_vector] == null - skip
		 *    minimizing sum of squared errors.
		 *
      Ui=Xi*Wv00+Yi*Wv01
      Vi=Xi*Wv10+Yi*Wv11

           sum(Xi^2) *Wv00 +sum(Xi*Yi)*Wv01- sum(Xi*Ui) =0
           sum(Xi*Yi)*Wv00 +sum(Yi^2) *Wv01- sum(Yi*Ui) =0


           sum(Xi^2) *Wv10 +sum(Xi*Yi)*Wv11- sum(Xi*Vi) =0
           sum(Xi*Yi)*Wv10 +sum(Yi^2) *Wv11- sum(Yi*Vi) =0

     S= |  sum(Xi^2)   sum(Xi*Yi) |
        |  sum(Xi*Yi)  sum(Yi^2)  |

     SU= | sum(Xi*Ui) |
         | sum(Yi*Ui) |

     SV= | sum(Xi*Vi) |
         | sum(Yi*Vi) |

     Wv0=| Wv00 |
         | Wv01 |

     Wv1=| Wv10 |
         | Wv11 |

     S * Wv0 = SU
     S * Wv1 = SV

     Wv0 = inv(S) * SU
     Wv1 = inv(S) * SV
		 *
		 */
		int i;
		double [][] S={{0.0,0.0},{0.0,0.0}};
		double []   SU={0.0,0.0};
		double []   SV={0.0,0.0};
		double [][] WV=new double [2][];
		for (i=0;i<vectors.length;i++) if (vectors[i]!=null) {
			//			   if (debugLevel>1) System.out.println("waveVectorsFromNeib: i="+i +": "+
			//					   vectors[i][0][0]+" "+vectors[i][0][1]+" "+vectors[i][1][0]+" "+vectors[i][1][1]+" ");
			S[0][0]+= vectors[i][1][0]*vectors[i][1][0]; // sum(Xi^2)
			S[0][1]+= vectors[i][1][0]*vectors[i][1][1]; // sum(Xi*Yi)
			S[1][1]+= vectors[i][1][1]*vectors[i][1][1]; // sum(Yi^2)
			SU[0]+=   vectors[i][1][0]*vectors[i][0][0]; // sum(Xi*Ui)
			SU[1]+=   vectors[i][1][1]*vectors[i][0][0]; // sum(Yi*Ui)
			SV[0]+=   vectors[i][1][0]*vectors[i][0][1]; // sum(Xi*Vi)
			SV[1]+=   vectors[i][1][1]*vectors[i][0][1]; // sum(Yi*Vi)
		}
		S[1][0]=S[0][1];
		//		   if (debugLevel>1) System.out.println("waveVectorsFromNeib: S00="+S[0][0]+" S01="+S[0][1]);
		//		   if (debugLevel>1) System.out.println("waveVectorsFromNeib: S10="+S[1][0]+" S11="+S[1][1]);
		//		   if (debugLevel>1) System.out.println("waveVectorsFromNeib: SU0="+SU[0]+  " SU1="+SU[1]);
		//		   if (debugLevel>1) System.out.println("waveVectorsFromNeib: SV0="+SV[0]+  " SV1="+SV[1]);

		S=matrix2x2_invert(S);
		//		   if (debugLevel>1) System.out.println("waveVectorsFromNeib: S00="+S[0][0]+" S01="+S[0][1]);
		//		   if (debugLevel>1) System.out.println("waveVectorsFromNeib: S10="+S[1][0]+" S11="+S[1][1]);

		WV[0]=matrix2x2_mul(S,SU);
		WV[1]=matrix2x2_mul(S,SV);
		return WV;
	}
	private void putInWaveList (
			List <Integer> list,
			int [] uv,
			int dir) {
		int l=(Integer.SIZE-2)/2;
		int mask =(1<<l)-1;

		list.add(new Integer((dir & 3) | ((uv[0] & mask) << 2) | ((uv[1] & mask) << (2+l))));
	}
	private int [] getWaveList (
			List <Integer> list,
			int index) {
		int l=(Integer.SIZE-2)/2;
		int mask =(1<<l)-1;

		int d=list.get(index);
		int [] result=new int[3];
		result[2]=(d & 3);
		result[0]=(d >>2 ) & mask;
		result[1]=(d >> (2+l)) & mask;
		return result;
	}


	/* ======================================================================== */
	// set XY coordinates and (optionally) wave vectors of the pattern grid cell
	/*
		   cell==null - new cell, not yet defined
		   cell.length==1 - invalid cell
		   cell.length>1 - initialized:
		   cell[0]==null - undefined
		   cell[0]!=null - defined
	 */
	private void setPatternGridCell(
			double [][][][] grid,
			int [] uv,
			double [] xy, // may be a 3-element, with contrast
			double [] wv1,
			double [] wv2){
		int i;
		initPatternGridCell(grid,uv);
		if (xy!=null) {
			//				  double [] grid_xy= new double[2];
			//				  for (i=0;i<2;i++) grid_xy[i]=xy[i];
			grid[uv[1]][uv[0]][0]= xy.clone(); // grid_xy;
		}
		if (wv1!=null) {
			double [] grid_wv1= new double[2];
			for (i=0;i<2;i++) grid_wv1[i]=wv1[i];
			grid[uv[1]][uv[0]][1]= grid_wv1;
		}
		if (wv2!=null) {
			double [] grid_wv2= new double[2];
			for (i=0;i<2;i++) grid_wv2[i]=wv2[i];
			grid[uv[1]][uv[0]][2]= grid_wv2;
		}
	}
	private void initPatternGridCell(
			double [][][][] grid,
			int [] uv){
		int i;
		if (grid[uv[1]][uv[0]]==null) {
			double [][] grid_cell= new double [3][];
			for (i=0;i<3;i++) grid_cell[i]=null;
			grid[uv[1]][uv[0]]=grid_cell;
		}
	}

	// mark the grid cell as invalid
	private void invalidatePatternGridCell(
			double [][][][] grid,
			int [] uv){
		double [][] cell = new double [1][];
		cell[0]=null;
		grid[uv[1]][uv[0]]=cell;
	}
	private void clearPatternGridCell(
			double [][][][] grid,
			int [] uv){
		grid[uv[1]][uv[0]]=null;
	}

	private void markDeletedPatternGridCell(
			double [][][][] grid,
			int [] uv){
		if ((grid[uv[1]][uv[0]]!=null) && (grid[uv[1]][uv[0]][0]!=null)) {
			double [] newXYC=new double[4];
			for (int i=0;i<newXYC.length;i++){
				if (i<grid[uv[1]][uv[0]][0].length)
					newXYC[i]=grid[uv[1]][uv[0]][0][i];
				else
					newXYC[i]=Double.NaN;
			}
			grid[uv[1]][uv[0]][0]=newXYC;
		}
		//			   grid[uv[1]][uv[0]]=null;

	}

	private boolean isCellDeleted(
			double [][][][] grid,
			int [] uv){
		return ((uv[1]>=0) && (uv[0]>=0) && (uv[1]<grid.length) && (uv[0]<grid[uv[1]].length) &&
				(grid[uv[1]][uv[0]]!=null) && (grid[uv[1]][uv[0]][0]!=null) && (grid[uv[1]][uv[0]][0].length>3));
	}
	private boolean isCellNew( //modified, for invalid uv will return "not new"
			double [][][][] grid,
			int [] uv){
		//	           return (uv[1]>=0) && (uv[0]>=0) && (uv[1]<grid.length) && (uv[0]<grid[uv[1]].length) && (grid[uv[1]][uv[0]]==null);
		// 4-th element is added to mark that the cell is dleted, but keep coordinates
		return (uv[1]>=0) && (uv[0]>=0) && (uv[1]<grid.length) && (uv[0]<grid[uv[1]].length) && ((grid[uv[1]][uv[0]]==null) || (grid[uv[1]][uv[0]].length>3));
	}
	private boolean isCellValid(
			double [][][][] grid,
			int [] uv){
		if ((uv[1]>=0) && (uv[0]>=0) && (uv[1]<grid.length) && (uv[0]<grid[uv[1]].length)) {
			double [][] cell = grid[uv[1]][uv[0]];
			return ((cell!=null) && (cell.length>1));
		}
		return false;
	}
	private boolean isCellDefined(
			double [][][][] grid,
			int [] uv){
		return ((uv[1]>=0) && (uv[0]>=0) && (uv[1]<grid.length) && (uv[0]<grid[uv[1]].length) &&
				(grid[uv[1]][uv[0]]!=null) && (grid[uv[1]][uv[0]][0]!=null));
	}
	private boolean isCellDefined(
			double [][][][] grid,
			int u,
			int v){
		return ((v>=0) && (u>=0) && (v<grid.length) && (u<grid[v].length) &&
				(grid[v][u]!=null) && (grid[v][u][0]!=null));
	}
	private boolean isCellDefined(
			int u,
			int v){
		return isCellDefined(this.PATTERN_GRID,u,v);
	}
	private boolean isCellDefined(
			int [] uv){
		return isCellDefined(this.PATTERN_GRID,uv);
	}

	// with contrast
	private double getCellContrast(double [][][][] grid,
			int [] uv){
		if  ((uv[1]>=0) && (uv[0]>=0) && (uv[1]<grid.length) && (uv[0]<grid[uv[1]].length) &&
				(grid[uv[1]][uv[0]]!=null) && (grid[uv[1]][uv[0]][0]!=null) && (grid[uv[1]][uv[0]][0].length>2)) {
			return grid[uv[1]][uv[0]][0][2];
		} else {
			return Double.NaN;
		}
	}
	private double getCellContrast(double [][][][] grid,
			int u,
			int v){
		if  ((v>=0) && (u>=0) && (v<grid.length) && (u<grid[v].length) &&
				(grid[v][u]!=null) && (grid[v][u][0]!=null) && (grid[v][u][0].length>2)) {
			return grid[v][u][0][2];
		} else {
			return Double.NaN;
		}
	}
	public double getCellContrast(int [] uv){
		return getCellContrast(this.PATTERN_GRID,uv);
	}
	public double getCellContrast(int u, int v){
		return getCellContrast(this.PATTERN_GRID,u,v);
	}
	private boolean isCellDefinedC(
			double [][][][] grid,
			int [] uv){
		return ((uv[1]>=0) && (uv[0]>=0) && (uv[1]<grid.length) && (uv[0]<grid[uv[1]].length) &&
				(grid[uv[1]][uv[0]]!=null) && (grid[uv[1]][uv[0]][0]!=null) && (grid[uv[1]][uv[0]][0].length>2) && !Double.isNaN(grid[uv[1]][uv[0]][0][2]));
	}
	private boolean isCellDefinedC(
			double [][][][] grid,
			int u,
			int v){
		return ((v>=0) && (u>=0) && (v<grid.length) && (u<grid[v].length) &&
				(grid[v][u]!=null) && (grid[v][u][0]!=null) && (grid[v][u][0].length>2) && !Double.isNaN(grid[v][u][0][2]));
	}
	public boolean isCellDefinedC(
			int u,
			int v){
		return isCellDefinedC(this.PATTERN_GRID,u,v);
	}
	public boolean isCellDefinedC(
			int [] uv){
		return isCellDefinedC(this.PATTERN_GRID,uv);
	}



	/*
		   private double [] cellXY(int u, int v){
			   if (!isCellDefined(u,v)) return null;
			   return this.PATTERN_GRID[v][u][0];
		   }
		   private double [] cellXY(int [] uv){
			   if (!isCellDefined(uv)) return null;
			   return this.PATTERN_GRID[uv[1]][uv[0]][0];
		   }
	 */
	private double [] cellXYC(int u, int v){
		if (!isCellDefined(u,v)) return null;
		double [] xyc={
				this.PATTERN_GRID[v][u][0][0],
				this.PATTERN_GRID[v][u][0][1],
				(this.PATTERN_GRID[v][u][0].length>2)?this.PATTERN_GRID[v][u][0][2]:
					((this.gridContrastBrightness==null)?1.0:this.gridContrastBrightness[0][v][u])
		};
		return xyc; // this.PATTERN_GRID[uv[1]][uv[0]][0];
	}
	private double [] cellXYC(int [] uv){
		if (!isCellDefined(uv)) return null;
		double [] xyc={
				this.PATTERN_GRID[uv[1]][uv[0]][0][0],
				this.PATTERN_GRID[uv[1]][uv[0]][0][1],
				(this.PATTERN_GRID[uv[1]][uv[0]][0].length>2)?this.PATTERN_GRID[uv[1]][uv[0]][0][2]:
					((this.gridContrastBrightness==null)?1.0:this.gridContrastBrightness[0][uv[1]][uv[0]])
		};
		return xyc; // this.PATTERN_GRID[uv[1]][uv[0]][0];
	}



	public int numDefinedCells()  {
		return numDefinedCells(this.PATTERN_GRID);
	}

	public int numDefinedCells(double [][][][] grid)  { // calulate/print number of defined nodes in a grid
		int [] iUV=new int [2];
		int numDefinedCells=0;
		for (iUV[1]=0;iUV[1]<grid.length;iUV[1]++) for (iUV[0]=0;iUV[0]<grid[0].length;iUV[0]++)
			if (this.isCellDefined(grid,iUV)) numDefinedCells++;
		return numDefinedCells;
	}
	public int gridUVWidth(){return ((this.PATTERN_GRID==null) || (this.PATTERN_GRID.length==0l) ||(this.PATTERN_GRID[0]==null))?0:this.PATTERN_GRID[0].length;}
	public int gridUVHeight(){return ((this.PATTERN_GRID==null) || (this.PATTERN_GRID.length==0l) ||(this.PATTERN_GRID[0]==null))?0:this.PATTERN_GRID.length;}
	/* ======================================================================== */
	/**
	 * returns number of laser pointers matched (or negative error)
	 * if (this.flatFieldForGrid!=null) it should already be applied !!
	 */
	public int calculateDistortions(
			LwirReaderParameters lwirReaderParameters, // null is OK
			MatchSimulatedPattern.DistortionParameters distortionParameters, //
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			SimulationPattern.SimulParameters  simulParameters,
			boolean equalizeGreens,
			ImagePlus imp, // image to process // has WOI_TOP and possibly - WOI_COMPENSATED
			LaserPointer laserPointer, // LaserPointer object or null
			boolean removeOutOfGridPointers, //
			double [][][] hintGrid, // predicted grid array (or null)
			double        hintGridTolerance, // allowed mismatch (fraction of period) or 0 - orientation only
			int threadsMax,
			boolean updateStatus,
			int global_debug_level, // DEBUG_LEVEL
			int debug_level, // debug level used inside loops
			boolean noMessageBoxes ){
		if (imp==null){
			IJ.showMessage("Error","There are no images open\nProcess canceled");
			return 0;
		}
		final int debugThreshold = 1;
		boolean is_lwir = ((lwirReaderParameters != null) && lwirReaderParameters.is_LWIR(imp));
		double    min_half_period= (is_lwir ? patternDetectParameters.minGridPeriodLwir : patternDetectParameters.minGridPeriod)/2;
		double    max_half_period=(is_lwir ? patternDetectParameters.maxGridPeriodLwir : patternDetectParameters.maxGridPeriod)/2;
		int       minimal_pattern_cluster = is_lwir ? distortionParameters.minimalPatternClusterLwir : distortionParameters.minimalPatternCluster;

		int fft_size = is_lwir ? distortionParameters.FFTSize_lwir : distortionParameters.FFTSize;

		long 	  startTime=System.nanoTime();
		Roi roi= imp.getRoi();
		Rectangle selection;
		if (roi==null){
			setWOI(0, 0, imp.getWidth(), imp.getHeight());
			selection=new Rectangle(0, 0, imp.getWidth(), imp.getHeight());
		} else {
			setWOI(roi.getBounds());
			selection=roi.getBounds();
		}
		this.debugLevel=global_debug_level;
		int patternCells=0;
		// save initial distortionParameters.correlationMinInitialContrast
		double savedCorrelationMinInitialContrast=distortionParameters.correlationMinInitialContrast;
		int reTries= 10; // bail out after these attempts
		boolean foundGoodCluster=false;

		int tryHor=0,tryVert=0;
		//with distortionParameters.searchOverlap==0.5 (default) step will be FFTSize original pixels, so half of the (2xFFTSize) square processed simultaneously
		if (distortionParameters.searchOverlap<0.1) distortionParameters.searchOverlap=0.1;
		int effectiveWidth=(int) (selection.width*0.5/distortionParameters.searchOverlap);
		int effectiveHeight=(int) (selection.height*0.5/distortionParameters.searchOverlap);

		for (int i = fft_size; i < effectiveWidth; i *= 2) tryHor++;
		for (int i = fft_size; i < effectiveHeight; i *= 2) tryVert++;

		int numTries=1<<(tryHor+tryVert);
		boolean [] triedIndices=new boolean[numTries+1]; // last set - all used
		for (int i=0;i<triedIndices.length;i++) triedIndices[i]=(i<3); // mark first 3 as if they are already used


		while (reTries-->0) {
			this.PATTERN_GRID=null;
			invalidateCalibration();

			patternCells=distortions( // calculates matchSimulatedPattern.DIST_ARRAY // invalidates calibration, flatFieldForGrid, resets this.PATTERN_GRID
					lwirReaderParameters, // null is OK
					triedIndices,
					distortionParameters, //
					patternDetectParameters,
					min_half_period,
					max_half_period,
					simulParameters,
					equalizeGreens,
					imp,
					threadsMax,
					updateStatus,
					debug_level,
					global_debug_level); // debug level
					if (global_debug_level>0) System.out.println("Pattern correlation done at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+
							" found "+patternCells+" cells, reTries left: "+reTries);
					if (patternCells>0) {
						foundGoodCluster=true;
						break; // new distortions() code - returns non-zero only if passed other tests
					}

					boolean someLeft=false;
					int startScanIndex=0;
					for (startScanIndex=3;startScanIndex<triedIndices.length;startScanIndex++) if (!triedIndices[startScanIndex]){
						someLeft=true;
						break;
					}

					if (someLeft) {
						if (global_debug_level>0){
							System.out.println("Initial pattern cluster is too small ("+patternCells+
									"), continuing scanning from index "+startScanIndex);
						}
					} else {
						if (global_debug_level>0) System.out.println("--- Tried all - nothing found --- at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
						break;
					}
		}



		// restore initial distortionParameters.correlationMinInitialContrast
		distortionParameters.correlationMinInitialContrast=savedCorrelationMinInitialContrast;
		if (!foundGoodCluster){
			if (global_debug_level > (debugThreshold + 1)) System.out.println("calculateDistortions(): Pattern too small, initial cluster had "+patternCells+" cells");
			if (global_debug_level > (debugThreshold + 2)) IJ.showMessage("Error","Pattern too small: "+patternCells);
			return distortionParameters.errPatternNotFound;
		}
		if (!patternOK()) {
			if (global_debug_level > (debugThreshold + 1)) System.out.println("Pattern not found");
			if (global_debug_level > (debugThreshold+2))   IJ.showMessage("Error","Pattern not found");
			return distortionParameters.errPatternNotFound;
		} else {
			if (global_debug_level > (debugThreshold+1)) System.out.println("Initial pattern cluster has "+patternCells+" cells");
		}
		if (global_debug_level > (debugThreshold+1)) System.out.println("Wave vectors recalculated at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		recalculateWaveVectors (
				updateStatus,
				debug_level);// debug level used inside loops
		ImagePlus imp_eq;
		if (distortionParameters.flatFieldCorrection && (this.flatFieldForGrid==null)) // if it is not null it is already supposed to be applied!
			imp_eq=equalizeGridIntensity(
					imp,
					this.PATTERN_GRID,
					distortionParameters, //
					equalizeGreens,
					global_debug_level,
					updateStatus,
					threadsMax);
		else imp_eq=imp;

		if (distortionParameters.refineCorrelations) {
			refineDistortionCorrelation (
					lwirReaderParameters, // LwirReaderParameters lwirReaderParameters, // null is OK
					distortionParameters, //
					patternDetectParameters,
					simulParameters,
					equalizeGreens,
					imp_eq,
					0.0, //final double maxCorr, // maximal allowed correction, in pixels (0.0) - any
					threadsMax,
					updateStatus,
					debug_level); // debug level

			recalculateWaveVectors (
					updateStatus,
					debug_level);// debug level used inside loops
			if (global_debug_level > (debugThreshold+1)) System.out.println("Second pass over at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		}
		//hack gridSize
		if ((distortionParameters.gridSize & 1)!=0) {
			refineDistortionCorrelation (
					lwirReaderParameters, // LwirReaderParameters lwirReaderParameters, // null is OK
					distortionParameters, //
					patternDetectParameters,
					simulParameters,
					equalizeGreens,
					imp_eq,
					0.0, //final double maxCorr, // maximal allowed correction, in pixels (0.0) - any
					threadsMax,
					updateStatus,
					debug_level); // debug level

			recalculateWaveVectors (
					updateStatus,
					debug_level);// debug level used inside loops
			if (global_debug_level>0) System.out.println("Third pass over at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
			//hack gridSize
		}
		patternCells=numDefinedCells();
		if ((roi!=null) && (patternCells<minimal_pattern_cluster)){
			if (global_debug_level > (debugThreshold+0)) System.out.println("Detected pattern is too small: "+patternCells+
					", minimum is set to "+minimal_pattern_cluster);
			return distortionParameters.errTooFewCells; // -10
		}

		double [] xy0={simulParameters.offsetX,simulParameters.offsetY} ; //debug
		createUV_INDEX(
				imp, // or null - just to determine WOI (when getWOI matches image size)
				xy0, // add to patterGrid xy, null OK
				threadsMax,
				updateStatus,
				global_debug_level, // DEBUG_LEVEL
				debug_level); // debug level used inside loops
		finalizeDistortionsBorder (
				distortionParameters, //
				updateStatus,
				debug_level);// debug level used inside loops
		//Wave vectors are used when calculating PSF
		recalculateWaveVectors (
				updateStatus,
				debug_level);// debug level used inside loops

		int numDifferentFFT=0;
		int maxLn2=0;
		for (int i=0;i<getCorrelationSizesUsed().length;i++) if (getCorrelationSizesUsed()[i]) {
			numDifferentFFT++;
			maxLn2=i;
		}
		if (numDifferentFFT>1){
			String sizesUsed="";
			for (int i=0;i<getCorrelationSizesUsed().length;i++) if (getCorrelationSizesUsed()[i]) sizesUsed+=" "+(1<<i);
			String msg="Different correlation FFT sizes used:"+sizesUsed+". You may consider increasing \"Correlation size\" setting to "+(1<<maxLn2)+" to reduce artifacts";
			if (global_debug_level > (debugThreshold+0)){
				System.out.println(msg);
				if (global_debug_level > (debugThreshold + 1)) IJ.showMessage(msg);
			}

		} else if (numDifferentFFT>0){
			String msg="Single correlation FFT size used: "+(1<<maxLn2);
			if (global_debug_level > (debugThreshold+0)) System.out.println(msg);
		}
		zeroNaNContrast(); // replace grid NaN with 0

		int numPointers=(laserPointer!=null)?laserPointer.laserUVMap.length:0;
		double [][] pointersXY=(numPointers>0)?getPointersXYUV(imp, laserPointer):null;
		if (global_debug_level > (debugThreshold+1)){
			if (pointersXY!=null){
				System.out.println("calculateDistortions() numPointers="+numPointers+" pointersXY.length="+pointersXY.length);
				for (int ii=0;ii<pointersXY.length;ii++) {
					if (pointersXY[ii]!=null){
						System.out.println("calculateDistortions()  pointersXY["+ii+"][0]="+pointersXY[ii][0]);
						System.out.println("                        pointersXY["+ii+"][1]="+pointersXY[ii][1]);
					} else{
						System.out.println("calculateDistortions()  pointersXY["+ii+"]=NULL");
					}
				}
				System.out.println("                                     hintGrid="+((hintGrid==null)?"NULL":"not NULL"));
				System.out.println("                            hintGridTolerance="+hintGridTolerance);
			} else {
				System.out.println("pointersXY == null");
			}
		}
		return combineGridCalibration(
				laserPointer, // LaserPointer object or null
				pointersXY,
				removeOutOfGridPointers, //
				hintGrid, // predicted grid array (or null)
				hintGridTolerance, // alllowed mismatch (fraction of period) or 0 - orientation only
				global_debug_level, // DEBUG_LEVEL
				noMessageBoxes );
	}

	//====== end of calculateDistortions() ==============================================
	/**
	 * Approximate function z(x,y) as a second degree polynomial
	 * f(x,y)=A*x^2+B*y^2+C*x*y+D*x+E*y+F
	 * data array consists of lines of either 2 or 3 vectors:
	 *  2-element vector x,y
	 *  variable length vector z (should be the same for all samples)
	 *  optional 1- element vector w (weight of the sample)
	 *
	 * returns array of vectors or null
	 * each vector (one per each z component) is either 6-element-  (A,B,C,D,E,F) if quadratic is possible and enabled
	 * or 3-element - (D,E,F) if linear is possible and quadratic is not possible or disbled
	 * returns null if not enough data even for the linear approximation
	 */

	public double [][] approximatePSFQuadratic(
			double []       psf,     // PSF function, square array, nominally positive
			double cutoffEnergy,     // fraction of energy in the pixels to be used
			double cutoffLevel,      // minimal level as a fraction of maximal
			int         minArea,      // minimal selected area in pixels
			double      blurSigma,    // optionally blur the selection
			double      maskCutOff,
			int           debugLevel, // debug level
			String        title) {    // prefix used for debug images
		double [] mask=findClusterOnPSF(
				psf,
				cutoffEnergy,
				cutoffLevel,
				minArea,
				blurSigma,
				debugLevel,
				title);
		int numPix=0;
		for (int i=0;i<mask.length;i++)
			if (mask[i]<maskCutOff) mask[i]=0.0;
		for (int i=0;i<mask.length;i++)
			if (mask[i]>0.0) numPix++;
		double [][][]data = new double[numPix][3][];
		numPix=0;
		int size = (int) Math.sqrt(psf.length);
		int hsize=size/2;
		for (int i=0;i<mask.length;i++)  if (mask[i]>0.0) {
			data[numPix][0]=new double[2];
			data[numPix][0][0]=(i % size) - hsize;
			data[numPix][0][1]=(i / size) - hsize;
			data[numPix][1]=new double[1];
			data[numPix][1][0]=psf[i];
			data[numPix][2]=new double[1];
			data[numPix][2][0]=mask[i];
			numPix++;
		}
		return  new PolynomialApproximation(debugLevel).quadraticApproximation(
				data,
				false,  // use linear approximation (instead of quadratic)
				1.0E-10,  // thershold ratio of matrix determinant to norm for linear approximation (det too low - fail)
				1.0E-20);  // thershold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
	}

	//====================================================
	public double [] tangetRadialSizes(
			double ca, // cosine of the center to sample vector
			double sa, // sine of the center to sample vector
			double []       psf,     // PSF function, square array, nominally positive
			double cutoffEnergy,     // fraction of energy in the pixels to be used
			double cutoffLevel,      // minimal level as a fraction of maximal
			int         minArea,      // minimal selected area in pixels
			double      blurSigma,    // optionally blur the selection
			double      maskCutOff,
			int           debugLevel, // debug level
			String        title) {    // prefix used for debug images
		double [] mask=findClusterOnPSF(
				psf,
				cutoffEnergy,
				cutoffLevel,
				minArea,
				blurSigma,
				debugLevel,
				title);
		for (int i=0;i<mask.length;i++)
			if (mask[i]<maskCutOff) mask[i]=0.0;
		int size = (int) Math.sqrt(psf.length);
		int hsize=size/2;
		//			   int nn=0;
		double S0=0.0, SR=0.0, ST=0.0,SR2=0.0,ST2=0.0; //,SRT=0.0;
		for (int i=0;i<mask.length;i++)  if (mask[i]>0.0) {
			double x=(i % size) - hsize;
			double y=(i / size) - hsize;
			double rc= x*ca+ y*sa;
			double tc=-x*sa+ y*ca;
			double d=psf[i]*mask[i];
			S0+=d;
			SR+=d*rc;
			ST+=d*tc;
			SR2+=d*rc*rc;
			ST2+=d*tc*tc;
			//				   nn++;
		}
		if (S0==0.0) return null; // make sure it is OK
		double  [] result={ Math.sqrt(ST2*S0 - ST*ST)/S0, Math.sqrt(SR2*S0 - SR*SR)/S0};
		//			   System.out.println(" mask.length="+mask.length+" nn="+nn+" S0="+S0+" SR="+SR+" ST="+ST+" SR2="+SR2+" ST2="+ST2+
		//					   " result={"+result[0]+","+result[1]+"}");
		return result;
	}

	//====================================================
	public double [] x2y2xySizes(
			double []       psf,     // PSF function, square array, nominally positive
			double cutoffEnergy,     // fraction of energy in the pixels to be used
			double cutoffLevel,      // minimal level as a fraction of maximal
			int         minArea,      // minimal selected area in pixels
			double      blurSigma,    // optionally blur the selection
			double      maskCutOff,
			int           debugLevel, // debug level
			String        title) {    // prefix used for debug images
		double [] mask=findClusterOnPSF(
				psf,
				cutoffEnergy,
				cutoffLevel,
				minArea,
				blurSigma,
				debugLevel,
				title);
		for (int i=0;i<mask.length;i++)
			if (mask[i]<maskCutOff) mask[i]=0.0;
		int size = (int) Math.sqrt(psf.length);
		int hsize=size/2;
		//			   int nn=0;
		double S0=0.0, SX=0.0, SY=0.0,SX2=0.0,SY2=0.0,SXY=0.0;
		for (int i=0;i<mask.length;i++)  if (mask[i]>0.0) {
			double x=(i % size) - hsize;
			double y=(i / size) - hsize;
			double d=psf[i]*mask[i];
			S0+=d;
			SX+=d*x;
			SY+=d*y;
			SX2+=d*x*x;
			SY2+=d*y*y;
			SXY+=d*x*y;
			//				   nn++;
		}
		if (S0==0.0) return null; // make sure it is OK
		double  [] result={
				(SX2*S0 - SX*SX)/S0/S0,
				(SY2*S0 - SY*SY)/S0/S0,
				(SXY*S0 - SX*SY)/S0/S0}; // this may be negative
		//			   System.out.println(" mask.length="+mask.length+" nn="+nn+" S0="+S0+" SX="+SX+" SY="+SY+" SX2="+SXR2+" SY2="+SY2+" SXY="+SXY+
		//					   " result={"+result[0]+","+result[1]+","+result[2]+"}");
		return result;
	}


	//====================================================

	public double [] findClusterOnPSF(
			double []        psf,     // PSF function, square array, nominally positive
			double cutoffEnergy,     // fraction of energy in the pixels to be used
			double cutoffLevel,      // minimal level as a fraction of maximal
			int         minArea,      // minimal selected area in pixels
			double      blurSigma,    // optionally blur the selection
			int           debugLevel, // debug level
			String        title) {    // prefix used for debhug images
		//				int i,j;
		int ix,iy,ix1,iy1;
		List <Integer> pixelList=new ArrayList<Integer>(100);
		Integer Index=0, Index1,IndexMax;
		int size=(int) Math.sqrt(psf.length);
		//				int [][]clusterMap=new int[size][size];
		int [][] dirs={{-1,0},{-1,-1},{0,-1},{1,-1},{1,0},{1,1},{0,1},{-1,1}};
		int len=size*size;
		double [] clusterMap=new double[len];
		double full_energy=0.0;
		double maxValue=0;
		for (int i=0;i<len;i++) {
			clusterMap[i]=0.0;
			if (psf[i]>0.0) full_energy+=psf[i];
			if (maxValue<psf[i]) {
				maxValue=psf[i];
				Index=i;
			}
		}
		if (maxValue<=0.0){
			String msg="psf array does not contain any positive values";
			//	    		  IJ.showMessage("Error",msg);
			System.out.println("Error "+msg);
			throw new IllegalArgumentException (msg);
		}
		ix=Index % size;
		iy=Index / size;
		double theresholdLevel=maxValue*cutoffLevel;
		double theresholdEnergy=full_energy*cutoffEnergy;
		double cluster_energy=0.0;
		int clusterSize=0;

		boolean noNew=true;
		if (debugLevel>1)		System.out.println("findClusterOnPSF(): full_energy="+full_energy+" theresholdEnergy="+theresholdEnergy+
				" maxValue="+maxValue+ " theresholdLevel="+theresholdLevel);
		if (debugLevel>1)		System.out.println("findClusterOnPSF(): ix="+ix+" iy="+iy);
		IndexMax=0;
		int listIndex;
		pixelList.clear();
		pixelList.add (Index);
		clusterSize++;
		clusterMap[Index]=1.0;
		cluster_energy+=psf[Index];
		noNew=true;
		while ((pixelList.size()>0) &&
				((clusterSize<minArea) || (cluster_energy<theresholdEnergy))) { // will break from the loop if  (psf[Index] <theresholdLevel)
			/* Find maximal new neighbor */
			maxValue=0.0;
			listIndex=0;
			while (listIndex<pixelList.size()) {
				Index=pixelList.get(listIndex);
				iy=Index/size;
				ix=Index%size;
				noNew=true;
				for (int j=0;j<8;j++) if (((iy > 0 ) || (dirs[j][1]>=0)) && ((iy < (size-1) ) || (dirs[j][1]<=0))){
					ix1=(ix+dirs[j][0]+size) % size;
					iy1= iy+dirs[j][1];
					Index1=iy1*size+ix1;
					if (clusterMap[Index1]==0.0) {
						noNew=false;
						if (psf[Index1]>maxValue) {
							maxValue= psf[Index1];
							IndexMax=Index1;
						}
					}
				}
				if (noNew) pixelList.remove(listIndex);  //  remove current list element
				else       listIndex++;     // increase list index
			}
			if (maxValue==0.0) 	break; // no positive points left
			if ((clusterSize>=minArea) && (psf[IndexMax]<theresholdLevel)) break; // level is below thershold, minimal size condition met
			/* Add this new point to the list */
			pixelList.add (IndexMax);
			clusterSize++;
			clusterMap[IndexMax]=1.0;
			cluster_energy+=psf[IndexMax];
		} // end of while ((pixelList.size()>0) &&  ...)
		if (debugLevel>3)   System.out.println("findClusterOnPSF: cluster size is "+clusterSize);
		if (debugLevel>3) {
			SDFA_INSTANCE.showArrays(psf, size, size, title+"-psf");
		}
		if (debugLevel>2) {
			SDFA_INSTANCE.showArrays(clusterMap, size, size, title+"-clusterMap");
		}
		if (blurSigma>0.0){
			DoubleGaussianBlur gb=new DoubleGaussianBlur();
			gb.blurDouble(
					clusterMap,
					size,
					size,
					blurSigma,
					blurSigma,
					0.01);
			if (debugLevel>2) {
				SDFA_INSTANCE.showArrays(clusterMap, size, size, title+"-clusterMap-blured");
			}
		}
		return clusterMap;
	}


	/* ======================================================================== */
	/**
	 * Mask (assignes zero) to the flat-field array outside of the sample squares,]
	 * Clears grid nodes that do not have neighbors inside the sample squares.
	 * If (this.flatFieldForGrid==null) - creates mask of 1.0/0.0
	 * @param focusMeasurementParameters - parameters specifying probe points
	 */
	public void maskFocus(
			double x0,   // lens center on the sensor
			double y0,  // lens center on the sensor
			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters){
		if (this.PATTERN_GRID==null) {
			String msg="PATTERN_GRID array does not exist, exiting";
			IJ.showMessage("Error",msg);
			throw new IllegalArgumentException (msg);
		}
		if (this.flatFieldForGrid==null) {
			String msg="Flat field for grid array does not exist, exiting";
			IJ.showMessage("Error",msg);
			throw new IllegalArgumentException (msg);
		}
		int [][] dirs= {{1,0},{0,1},{-1,0},{0,-1},{1,1},{1,-1},{-1,1},{-1,-1}};
		int width=getImageWidth();
		int height=getImageHeight();
		int halfSize=focusMeasurementParameters.sampleSize/2;
		// System.out.println("maskFocus(): width="+width+" height="+height+" _halfSize="+halfSize);

		double [][][] sampleCoord= focusMeasurementParameters.sampleCoordinates(
				x0,   // lens center on the sensor
				y0);  // lens center on the sensor
		this.focusMask =new boolean[this.flatFieldForGrid.length];
		for (int i=0;i<this.focusMask.length;i++) this.focusMask[i]=false;
		for (int i=0;i<focusMeasurementParameters.numSamples[1];i++){
			//System.out.println(i+": y0="+y0);
			for (int j=0;j<focusMeasurementParameters.numSamples[0];j++){
				int xs=(int) sampleCoord[i][j][0];
				int ys=(int) sampleCoord[i][j][1];
				//System.out.println(j+": xs="+xs);
				for (int y=ys-halfSize;y<(ys+halfSize);y++) if ((y>=0) && (y<height)) {
					for (int x=xs-halfSize;x<(xs+halfSize);x++) if ((x>=0) && (x<width)) {
						this.focusMask[y*width+x]=true;
					}
				}
			}

		}
		// Do we really need to zero the this.flatFieldForGrid where this.focusMask is false? We may need
		// the pixels out of WOI to update grid around the needed nodes
		//			   for (int i=0;i<this.focusMask.length;i++) if (!this.focusMask[i]) this.flatFieldForGrid[i]=0.0;
		if (this.PATTERN_GRID.length==0) return;
		boolean [][] maskUV=new boolean[this.PATTERN_GRID.length][this.PATTERN_GRID[0].length];
		int [] iUV={0,0};
		for (iUV[1]=0;iUV[1]<maskUV.length;iUV[1]++) for (iUV[0]=0;iUV[0]<maskUV[0].length;iUV[0]++) {
			maskUV[iUV[1]][iUV[0]]=false;
			if (isCellDefined(this.PATTERN_GRID, iUV)){
				int x= (int) Math.round(this.PATTERN_GRID[iUV[1]][iUV[0]][0][0]);
				int y= (int) Math.round(this.PATTERN_GRID[iUV[1]][iUV[0]][0][1]);
				if ((x>=0) && (x<width) && (y>=0) && (y<height) && this.focusMask[y*width+x]) maskUV[iUV[1]][iUV[0]]=true;
			}
		}
		for (iUV[1]=0;iUV[1]<maskUV.length;iUV[1]++) for (iUV[0]=0;iUV[0]<maskUV[0].length;iUV[0]++) if (!maskUV[iUV[1]][iUV[0]]){
			boolean neibExists=false;
			for (int d=0;d<dirs.length;d++) {
				int [] iUV1={iUV[0]+dirs[d][0],iUV[1]+dirs[d][1]};
				if (isCellDefined(this.PATTERN_GRID, iUV1) && maskUV[iUV1[1]][iUV1[0]]){
					neibExists=true;
					break;
				}
			}
			if (!neibExists) clearPatternGridCell(this.PATTERN_GRID,iUV);

		}
	}
	/* ======================================================================== */
	@Deprecated
	public double[][] getPointersXY(ImagePlus imp, int numPointers){
		// try absolute grid calibratrion
		// read image info to properties (if it was not done yet - should it?
		if ((imp.getProperty("timestamp")==null) || (((String) imp.getProperty("timestamp")).length()==0)) {
			JP46_Reader_camera jp4_instance= new JP46_Reader_camera(false);
			jp4_instance.decodeProperiesFromInfo(imp);
		}
		double [][] pointersXY=new double[numPointers][];
		int numPointerDetected=0;
		for (int i=0;i<pointersXY.length;i++) {
			pointersXY[i]=null;
			if ((imp.getProperty("POINTER_X_"+i)!=null) && (imp.getProperty("POINTER_Y_"+i)!=null)) {
				pointersXY[i]=new double[2];
				pointersXY[i][0]=Double.parseDouble((String) imp.getProperty("POINTER_X_"+i));
				pointersXY[i][1]=Double.parseDouble((String) imp.getProperty("POINTER_Y_"+i));
				numPointerDetected++;
			}
		}
		if (numPointerDetected>0) return pointersXY;
		else return null;
	}

	// laserPointer used as a backup if images do not contain data, and specify if they are needed at all (null - skip)
	public static double[][] getPointersXYUV(ImagePlus imp, LaserPointer laserPointer ){
		int numPointers=(laserPointer!=null)?laserPointer.laserUVMap.length:0;
		// try absolute grid calibratrion
		// read image info to properties (if it was not done yet - should it?
		if ((imp.getProperty("timestamp")==null) || (((String) imp.getProperty("timestamp")).length()==0)) {
			JP46_Reader_camera jp4_instance= new JP46_Reader_camera(false);
			jp4_instance.decodeProperiesFromInfo(imp);
		}
		double [][] pointersXY=new double[numPointers][];
		int numPointerDetected=0;
		// Normally all pointers should either have or not U,V - the smallest one will later be used
		for (int i=0;i<pointersXY.length;i++) {
			pointersXY[i]=null;
			if ((imp.getProperty("POINTER_X_"+i)!=null) && (imp.getProperty("POINTER_Y_"+i)!=null)) {
				if ((imp.getProperty("POINTER_U_"+i)!=null) && (imp.getProperty("POINTER_V_"+i)!=null)) {
					pointersXY[i]=new double[4];
					pointersXY[i][2]=Double.parseDouble((String) imp.getProperty("POINTER_U_"+i));
					pointersXY[i][3]=Double.parseDouble((String) imp.getProperty("POINTER_V_"+i));
				} else {
					pointersXY[i]=new double[2];
				}
				pointersXY[i][0]=Double.parseDouble((String) imp.getProperty("POINTER_X_"+i));
				pointersXY[i][1]=Double.parseDouble((String) imp.getProperty("POINTER_Y_"+i));
				numPointerDetected++;
			}
		}
		if (numPointerDetected>0) return pointersXY;
		else return null;
	}

	public static void setPointersXYUV(ImagePlus imp, double [][] pointersXYUV){

//		JP46_Reader_camera jp4_instance= new JP46_Reader_camera(false);
//		jp4_instance.decodeProperiesFromInfo(imp);
		String prefix = "POINTER_";
		Properties prop = imp.getProperties();
		// Delete all properties starting with "POINTER_"
		Enumeration<Object> enumKey = prop.keys();
		while(enumKey.hasMoreElements()) {
			Object key = enumKey.nextElement();
			if (((String) key).startsWith(prefix)) {
				prop.remove(key);
			}
		}
		// Delete all properties starting with "POINTER_"
		int num_used = 0;
		for (int i=0;i<pointersXYUV.length;i++) if (pointersXYUV[i] != null) {
			prop.setProperty(prefix+"X_"+i, ""+pointersXYUV[i][0]);
			prop.setProperty(prefix+"Y_"+i, ""+pointersXYUV[i][1]);
			if (pointersXYUV[i].length > 2) {
				prop.setProperty(prefix+"U_"+i, ""+pointersXYUV[i][2]);
				prop.setProperty(prefix+"V_"+i, ""+pointersXYUV[i][3]);
			}
			num_used++;
		}
		prop.setProperty("USED_POINTERS",""+num_used);
//		jp4_instance.encodeProperiesToInfo(imp);
//		(new FileSaver(imp)).saveAsTiff(path);

	}




	/* ======================================================================== */

	public String getChannel(ImagePlus imp){
		// read image info to properties (if it was not done yet - should it?
		if ((imp.getProperty("timestamp")==null) || (((String) imp.getProperty("timestamp")).length()==0)) {
			JP46_Reader_camera jp4_instance= new JP46_Reader_camera(false);
			jp4_instance.decodeProperiesFromInfo(imp);
		}
		if (imp.getProperty("channel")==null) return null;
		return (String) imp.getProperty("channel");
	}

	/* ======================================================================== */
	public void showFlatFieldForGrid(){
		if (this.flatFieldForGrid!=null) this.SDFA_INSTANCE.showArrays(this.flatFieldForGrid, getImageWidth(), getImageHeight(), "Flat_field_for_grid");
	}
	public void showFFCorrectedGrid(){
		if (this.gridFFCorr!=null) this.SDFA_INSTANCE.showArrays(this.gridFFCorr, getImageWidth(), getImageHeight(), "Flat_field_corrected_grid");
	}
	public void showFocusMask(){
		if (this.focusMask!=null){
			double [] dfm=new double [this.focusMask.length];
			for (int i=0;i<dfm.length;i++) dfm[i]=this.focusMask[i]?1.0:0.0;
			this.SDFA_INSTANCE.showArrays(dfm, getImageWidth(), getImageHeight(), "Focus_mask");
		}
	}
	public void showUVIndex(){
		if (this.UV_INDEX!=null){
			double [] uv=new double [this.UV_INDEX.length];
			for (int i=0;i<uv.length;i++) uv[i]=this.UV_INDEX[i];
			this.SDFA_INSTANCE.showArrays(uv, getImageWidth(), getImageHeight(), "UV_INDEX");
		}
	}
	public boolean patternOK(){
		return (this.PATTERN_GRID!=null);
	}
	public double [][][][] getDArray(){
		return this.PATTERN_GRID;
	}
	public double [][][] getDArray(int v){
		return this.PATTERN_GRID[v];
	}
	public double [][] getDArray(int v, int u){
		return this.PATTERN_GRID[v][u];
	}
	public double [] getDArray(int v, int u, int n){
		return this.PATTERN_GRID[v][u][n];
	}
	public double getDArray(int v, int u, int n, int k){
		return this.PATTERN_GRID[v][u][n][k];
	}
	public double [] getXY (int v, int u){
		return this.PATTERN_GRID[v][u][0];
	}
	public int getDArrayHeight(){
		if (this.PATTERN_GRID==null) return 0;
		return this.PATTERN_GRID.length;
	}
	public int getDArrayWidth(){
		if (this.PATTERN_GRID==null) return 0;
		return this.PATTERN_GRID[0].length;
	}
	//			matchSimulatedPattern.DIST_SELECTION.width, // image (mask) width
	public Rectangle getWOI(){
		return this.DIST_SELECTION;
	}
	public int getImageWidth(){
		return this.UV_INDEX_WIDTH;
	}
	public int getImageHeight(){
		if (this.UV_INDEX==null) return 0;
		return this.UV_INDEX.length/this.UV_INDEX_WIDTH;
	}
	public void setWOI(Rectangle woi){
		this.DIST_SELECTION=new Rectangle(woi) ;
	}
	public void setWOI(int x, int y, int w, int h){
		this.DIST_SELECTION=new Rectangle(x, y, w, h) ;
	}
	public int [] getUVIndex(){
		return this.UV_INDEX;
	}
	public int getUVIndex(int i){
		if ((i<0) || (i>=this.UV_INDEX.length)) return -1; // let it throw exception
		return this.UV_INDEX[i];
	}
	public int getUVIndex(int x, int y){
		if ((x<0) || (y<0) || (x>=this.PATTERN_GRID[0].length)  || (y>=this.PATTERN_GRID.length)) return -1;
		return this.UV_INDEX[x+UV_INDEX_WIDTH*y];
	}
	/**
	 * Returns a pair of U,V from xy, integer result (does not interpolate
	 * @param x - pixel coordinat X
	 * @param y - pixel coordinat Y
	 * @return {u,v} pair
	 */
	public int [] getUV(int x, int y)    {
		if ((x<0) || (y<0) || (x>=this.UV_INDEX_WIDTH) || (y>=(this.UV_INDEX.length/this.UV_INDEX_WIDTH))) return null; // out of UV_INDEX bounds
		if ((x+this.UV_INDEX_WIDTH*y)>this.UV_INDEX.length) {
			if (this.debugLevel>0){
				System.out.println("getUV("+x+","+y+"): this.UV_INDEX.length="+this.UV_INDEX.length+", this.UV_INDEX_WIDTH="+this.UV_INDEX_WIDTH);
			}
			return null;
		}
		int index=this.UV_INDEX[x+this.UV_INDEX_WIDTH*y];
		if (index<0) return null; // <0 - undefined
		int width=this.PATTERN_GRID[0].length;
		int [] uv={index % width,index / width};
		return uv;
	}
	/**
	 * Map estimated grid hintGrid  to measured
	 * @param hintGrid [v][u][ 0-x, 1-y, 2 - u, 3- v]
	 * @param searchAround  how far (in pixels) to look for the nearest defined one if the specified is undefined
	 * @return array [v][u] [0 - measured u, 1 - measured v, 2 - measured contrast] ([v][u]==null no measured grid there
	 */

	public double [][][] mapEstimatedToMeasured(double [][][] hintGrid, double searchAround){
		double [][][] measuredUV=new double[hintGrid.length][hintGrid[0].length][];
		for (int v=0;v<measuredUV.length;v++) for (int u=0;u<measuredUV[0].length;u++)
			if (hintGrid[v][u]!=null) measuredUV[v][u]=getUVLinear(hintGrid[v][u][0],hintGrid[v][u][1],searchAround);
			else measuredUV[v][u]=null;
		return measuredUV;
	}

	/**
	 * Calculate linear matrix (2x3) parameters measured grid UV from estimated grid UV
	 * @param hintGrid  [v][u][ 0-x, 1-y, 2 - u, 3- v]
	 * @param searchAround how far (in pixels) to look for the nearest defined one if the specified is undefined
	 * @return {{Au, Bu, Cu}.{Av, Bc, Cv}}; where Umeas=Au*Uhint+Bu*Vhint+Cu, Vmeas=Av*Uhint+Bv*Vhint+Cv
	 */
	public double [][] calcGridMatchMatrix (double [][][] hintGrid, double searchAround){
		double [][][] measuredUV=mapEstimatedToMeasured(hintGrid, searchAround); // contrast - minimal of the ones around
		if (this.debugLevel>2){
			double [][] pixels=new double[9][measuredUV.length*measuredUV[0].length];
			int index=0;
			String [] titles={"grid-U","grid-V","contrast","hint-X","hint-Y","hint-U","hint-V","grid-hint-U","grid-hint-V"}; // last 2 only valid for rotation==0
			for (int v=0; v<measuredUV.length;v++) for (int u=0;u<measuredUV[v].length;u++){
				if (measuredUV[v][u]!=null){
					for (int i=0; i<3;i++)	pixels[i][index]=measuredUV[v][u][i];
					for (int i=0; i<4;i++)	pixels[i+3][index]=hintGrid[v][u][i];
					for (int i=0; i<2;i++)	pixels[i+7][index]=measuredUV[v][u][i]-hintGrid[v][u][i+2];
				} else {
					for (int i=0; i<9;i++)	pixels[i][index]=Double.NaN;
				}
				index++;
			}
			(new ShowDoubleFloatArrays()).showArrays(pixels, measuredUV[0].length, measuredUV.length,  true, "measuredUV", titles);

		}

		int numDefined=0;
		for (int v=0;v<measuredUV.length;v++) for (int u=0;u<measuredUV[v].length;u++) if (measuredUV[v][u]!=null) numDefined++;
		double [][][] data =new double [numDefined][3][];//        	   double [][][] data =new double [numDefined][2][2];
		int index=0;
		for (int v=0;v<measuredUV.length;v++) for (int u=0;u<measuredUV[v].length;u++) if (measuredUV[v][u]!=null) {
			data[index][0]=new double[2];
			data[index][1]=new double[2];
			data[index][2]=new double[1];
			data [index][0][0]=hintGrid[v][u][2]; // hinted U
			data [index][0][1]=hintGrid[v][u][3]; // hinted V
			data [index][1][0]=measuredUV[v][u][0]; // measured U
			data [index][1][1]=measuredUV[v][u][1]; // measured V
			data [index][2][0]=measuredUV[v][u][2]; // contrast
			index++;
		}
		if (this.debugLevel>0) {
			System.out.println("calcGridMatchMatrix(), data.length="+data.length+" measuredUV.length="+measuredUV.length+
					((measuredUV.length>0)?(", measuredUV[0].length="+measuredUV[0].length):""));
			if (this.debugLevel>3) for (index=0;index<data.length;index++){
				System.out.println(data[index][0][0]+","+data[index][0][1]+","+data[index][1][0]+","+data[index][1][1]+","+data[index][2][0]);
			}
		}
		//			   int gridRotation=-1; //undefined
		if (data.length<3){
			///        		   if (this.debugLevel>0) System.out.println("calcGridMatchMatrix(), data.length="+data.length+" measuredUV.length="+measuredUV.length+
			///        				   ((measuredUV.length>0)?(", measuredUV[0].length="+measuredUV[0].length):""));
			return null;
		}
		double [][] coeff=new PolynomialApproximation(this.debugLevel).quadraticApproximation(data, true);
		if (coeff!=null) {
			//				   gridRotation=matrixToRot(coeff);
			int rot=matrixToRot(coeff);
			boolean [] flips=rotToFlips(rot);
			double [][] aI={{1,0},{0,1}};
			double [][] aSwap= {{ 0,1},{1,0}};
			double [][] aFlipU={{-1,0},{0,1}};
			double [][] aFlipV={{ 1,0},{0,-1}};
			Matrix M=new Matrix(aI);
			if (flips[0]) M=M.times((new Matrix(aSwap)));
			if (flips[1]) M=M.times((new Matrix(aFlipU)));
			if (flips[2]) M=M.times((new Matrix(aFlipV)));
			// now M reconstructs coeff
			double [][]aM=M.getArray();
			double SMU=0.0,SMV=0.0,SW=0.0;
			for (int i=0;i<data.length;i++){
				SMU+=data[i][2][0]*(data[i][1][0]-(aM[0][0]*data[i][0][0]+aM[0][1]*data[i][0][1]));
				SMV+=data[i][2][0]*(data[i][1][1]-(aM[1][0]*data[i][0][0]+aM[1][1]*data[i][0][1]));
				SW+=data[i][2][0];
			}
			SMU/=SW;
			SMV/=SW;
			if (this.debugLevel>0) {
				System.out.println("coeff[0][0]="+coeff[0][0]+" coeff[0][1]="+coeff[0][1]+" coeff[0][2]="+coeff[0][2]);
				System.out.println("coeff[1][0]="+coeff[1][0]+" coeff[1][1]="+coeff[1][1]+" coeff[1][2]="+coeff[1][2]);
				System.out.println("aM[0][0]="+aM[0][0]+" aM[0][1]="+aM[0][1]+" SMU="+SMU);
				System.out.println("aM[1][0]="+aM[1][0]+" aM[1][1]="+aM[1][1]+" SMV="+SMV);
			}
			coeff[0][0]=aM[0][0];
			coeff[0][1]=aM[0][1];
			coeff[1][0]=aM[1][0];
			coeff[1][1]=aM[1][1];
			coeff[0][2]=SMU;
			coeff[1][2]=SMV;

		}
		return coeff;
	}
	public int matrixToRot(double [][] coeff){
		boolean [] flips = {false,false,false};
		double  [][] aR={{coeff[0][0],coeff[0][1]},{coeff[1][0],coeff[1][1]}};
		double  [][] aSwap={{0,1},{1,0}};
		Matrix R = new Matrix (aR);
		if ((aR[0][0]*aR[0][0]+aR[1][1]*aR[1][1]) < (aR[1][0]*aR[1][0]+aR[0][1]*aR[0][1])){
			flips[0]=true;
			R=(new Matrix(aSwap)).times(R);
		}
		flips[1]=R.getArray()[0][0]<0;
		flips[2]=R.getArray()[1][1]<0;
		return flipsToRot(flips[0],flips[1],flips[2]);
	}

	public int [][] gridMatrixApproximate(double [][] coeff){
		int rot=matrixToRot(coeff);
		boolean [] flips=rotToFlips(rot);
		double [][] aI={{1,0},{0,1}};
		double [][] aSwap= {{ 0,1},{1,0}};
		double [][] aFlipU={{-1,0},{0,1}};
		double [][] aFlipV={{ 1,0},{0,-1}};
		Matrix M=new Matrix(aI);
		if (flips[0]) M=M.times((new Matrix(aSwap)));
		if (flips[1]) M=M.times((new Matrix(aFlipU)));
		if (flips[2]) M=M.times((new Matrix(aFlipV)));
		// now M reconstructs coeff
		double [][]aM=M.getArray();
		// Black/white cells have to be flipped if flipU XOR flipW, regardless of swapUV
		int flipForWhite=(flips[1]^flips[2])?1:0;
		int [][] shifts={
				{2*((int) Math.round(0.5*(coeff[0][2]+0))),   2*((int) Math.round(0.5*(coeff[1][2]+flipForWhite)))    -flipForWhite},
				{2*((int) Math.round(0.5*(coeff[0][2]-1)))+1, 2*((int) Math.round(0.5*(coeff[1][2]+flipForWhite-1)))+1-flipForWhite}
		};
		int shiftSelect=(
				((shifts[0][0]-coeff[0][2])*(shifts[0][0]-coeff[0][2]) + (shifts[0][1]-coeff[1][2])*(shifts[0][1]-coeff[1][2]))>
				((shifts[1][0]-coeff[0][2])*(shifts[1][0]-coeff[0][2]) + (shifts[1][1]-coeff[1][2])*(shifts[1][1]-coeff[1][2])))?1:0;
		if (this.debugLevel>1){
			double d1=Math.sqrt((shifts[0][0]-coeff[0][2])*(shifts[0][0]-coeff[0][2]) + (shifts[0][1]-coeff[1][2])*(shifts[0][1]-coeff[1][2]));
			double d2=Math.sqrt((shifts[1][0]-coeff[0][2])*(shifts[1][0]-coeff[0][2]) + (shifts[1][1]-coeff[1][2])*(shifts[1][1]-coeff[1][2]));
			System.out.println("gridMatrixApproximate(): shifts[0][0]="+shifts[0][0]+" shifts[0][1]="+shifts[0][1]+
					" shifts[1][0]="+shifts[1][0]+" shifts[1][1]="+shifts[1][1]+ " shiftSelect="+shiftSelect+ " d1="+d1+" d2="+d2);
		}
		int [][] iCoeff={
				//        			   {(int) Math.round(aM[0][0]), (int) Math.round(aM[0][1]),(int) Math.round(coeff[0][2])},
				//        			   {(int) Math.round(aM[1][0]), (int) Math.round(aM[1][1]),(int) Math.round(coeff[1][2])}};
				{(int) Math.round(aM[0][0]), (int) Math.round(aM[0][1]),shifts[shiftSelect][0]},
				{(int) Math.round(aM[1][0]), (int) Math.round(aM[1][1]),shifts[shiftSelect][1]}};
		return iCoeff;
	}
	// old version, no distinction between B and W
	public int [][] gridMatrixApproximateNoBW(double [][] coeff){
		int rot=matrixToRot(coeff);
		boolean [] flips=rotToFlips(rot);
		double [][] aI={{1,0},{0,1}};
		double [][] aSwap= {{ 0,1},{1,0}};
		double [][] aFlipU={{-1,0},{0,1}};
		double [][] aFlipV={{ 1,0},{0,-1}};
		Matrix M=new Matrix(aI);
		if (flips[0]) M=M.times((new Matrix(aSwap)));
		if (flips[1]) M=M.times((new Matrix(aFlipU)));
		if (flips[2]) M=M.times((new Matrix(aFlipV)));
		// now M reconstructs coeff
		double [][]aM=M.getArray();
		int [][] iCoeff={
				{(int) Math.round(aM[0][0]), (int) Math.round(aM[0][1]),(int) Math.round(coeff[0][2])},
				{(int) Math.round(aM[1][0]), (int) Math.round(aM[1][1]),(int) Math.round(coeff[1][2])}};
		return iCoeff;
	}

	public double worstGridMatchRotSkew(double [][] coeff){
		int [][] iCoeff=gridMatrixApproximate(coeff);
		double worst=0;
		for (int i=0;i<2;i++) for (int j=0;j<2;j++){
			double d=Math.abs(coeff[i][j]-iCoeff[i][j]);
			if (d>worst) worst=d;
		}
		return worst;
	}
	public double worstGridMatchTranslate(double [][] coeff){ // in grids half-periods, not pixels!
		int [][] iCoeff=gridMatrixApproximate(coeff);
		double worst=0;
		for (int i=0;i<2;i++){
			double d=Math.abs(coeff[i][2]-iCoeff[i][2]);
			if (d>worst) worst=d;
		}
		return worst;
	}
	//searchAround

	/**
	 * Find double u,v from double x,y by linear interpolation from neighbor cells. Requires  this.UV_INDEX to be calculated and
	 * matching this.PATTERN_GRID
	 * @param x pixel coordinate
	 * @param y  pixel coordinate
	 * @param searchAround how far (in pixels) to look for the nearest defined one if the specified is undefined
	 * @return uv pair or null - modified - triplet, last - contrast
	 */
	public double [] getUVLinear(double x, double y, double searchAround){
		int ix0= (int) Math.round(x);
		int iy0= (int) Math.round(y);
		int    ix=ix0,iy=iy0;
		int [] uv0=getUV(ix, iy);
		double best2=searchAround*searchAround+1;
		// if the point is slightly out of grid - find the one near
		if (uv0==null) {
			for (int iy1=iy-((int) Math.round(searchAround));iy1>=iy+searchAround;iy1++)
				for (int ix1=ix-((int) Math.round(searchAround));ix1>=ix+searchAround;ix1++) {
					double d= (ix1-ix)*(ix1-ix)+(iy1-iy)*(iy1-iy);
					if (d<best2)  {
						uv0=getUV(ix1, iy1);
						if (uv0!=null){
							best2=d;
							ix=ix1;
							iy=iy1;
						}
					}
				}

		}
		if (uv0==null) return null; // no grid points near
		int [][] dirDiffs={{1,1},{-1,1},{1,-1},{-1,-1}};
		double [] xy0= cellXYC(uv0);
		double [] xy1=null;
		double [] xy2=null;
		int [] deltaUV=null;
		for (int dir=0;dir<dirDiffs.length;dir++) {
			xy1=cellXYC(uv0[0]+dirDiffs[dir][0],uv0[1]);
			xy2=cellXYC(uv0[0]                 ,uv0[1]+dirDiffs[dir][1]);
			if ((xy1!=null) && (xy2!=null)) {
				deltaUV=new int[2];
				deltaUV[0]=dirDiffs[dir][0];
				deltaUV[1]=dirDiffs[dir][1];
				break;
			}
		}
		double minContrast=Math.min(Math.min(xy1[2], xy2[2]),xy0[2]); // minimal contrast off all 3 points
		if (deltaUV==null) return null; // could not find 2 orthogonal neighbors to interpolate
		/*
x=xy0[0] + dU*deltaUV[0]*(xy1[0]-xy0[0])+dV*deltaUV[1]*(xy2[0]-xy0[0])
y=xy0[1] + dU*deltaUV[0]*(xy1[1]-xy0[1])+dV*deltaUV[1]*(xy2[1]-xy0[1])

		 */
		double [][] aM={
				{deltaUV[0]*(xy1[0]-xy0[0]),deltaUV[1]*(xy2[0]-xy0[0])},
				{deltaUV[0]*(xy1[1]-xy0[1]),deltaUV[1]*(xy2[1]-xy0[1])}};
		Matrix M=new Matrix(aM);
		double [][] aB={{x-xy0[0]},{y-xy0[1]}};
		Matrix B=new Matrix(aB);
		if (!(new LUDecomposition(M)).isNonsingular()){
			System.out.println("getUVLinear("+x+","+y+"): Matix is singular:");
			M.print(10, 6);
			return null;
		}
		double [] dUV=M.solve(B).getRowPackedCopy();
		double [] result={uv0[0]+dUV[0],uv0[1]+dUV[1],minContrast};
		return result;
	}
	public void invalidateAll(){
		invalidateCalibration();
		invalidateFlatFieldForGrid();
		invalidateFocusMask();
	}
	private void invalidateCalibration(){
		this.reMap=null;    // invalidate if any
		this.targetUV=null; // invalidate if any
		this.pXYUV=null; // invalidate if any
		this.passNumber=1;
		resetCorrelationSizesUsed(); // reset which FFT sizes where used in correlation
	}
	public void invalidateFlatFieldForGrid(){
		this.flatFieldForGrid=null; // reset flat field for grid
		this.gridContrastBrightness=null;
	}
	public void invalidateFocusMask(){
		this.focusMask=null;
	}
	/* get height and width of the measured pattern array applies to PATTERN_GRID, targetUV and pixelsUV */

	public int getHeight(){
		if (this.pXYUV==null) return 0;
		return this.pXYUV.length;
	}
	public int getWidth(){
		if (this.pXYUV==null) return 0;
		return this.pXYUV[0].length;
	}
	/* Get physical target UV pair from measured pattern. Requires absolute mapping (by laser spots)
	 *  Pair may be null if no pattern is detected for this node in the image
	 */
	public int [][][] getTargetUV(){
		return this.targetUV;
	}
	/* Get pixel X,Y pair for each node in the measured pattern. Calculated during absolute mapping (by laser spots)
	 *  Pair may be null if no pattern is detected for this node in the image
	 */
	public double [][][] getPXYUV(){
		return this.pXYUV;
	}

	public int restorePatternGridFromGridList(
			double [][][] pixelsXYSet,
			int [][][] pixelsUVSet,
			double [] intensityRange){
		double maxX=0,maxY=0;
		for (int n=0;n<pixelsXYSet.length;n++){
			for (int i=0;i<pixelsXYSet[n].length;i++){
				if (pixelsXYSet[n][i][0]>maxX) maxX=pixelsXYSet[n][i][0];
				if (pixelsXYSet[n][i][1]>maxY) maxY=pixelsXYSet[n][i][1];
			}
		}
		int width=(int) Math.ceil(maxX)+1;
		int height=(int) Math.ceil(maxY)+1;
		return restorePatternGridFromGridList(pixelsXYSet, pixelsUVSet, width,height,intensityRange);
	}

	/**
	 * Restore this.PATTERN_GRID array (no wave vectors) from the lists use in Distortions class
	 * @param pixelsXYSet list of the {x,y} pairs for each grid node (now - a pair of lists used pixels and those that did not fit into physical target)
	 * @param pixelsUVSet list of {u,v} pairs for each grid node (now - a pair of lists)
	 * @param width sensor (image) width
	 * @param height sensor (image) height
	 * @return number of cells in the grid
	 */

	public int restorePatternGridFromGridList(
			double [][][] pixelsXYSet,
			int [][][] pixelsUVSet,
			int width,
			int height,
			double [] intensityRange){
		int numCells=0;
		int minU=0,minV=0,maxU=0,maxV=0;
		this.PATTERN_GRID=null;
		setWOI(0, 0, 0, 0);
		if ((pixelsXYSet!=null) && (pixelsUVSet!=null ) &&  (pixelsXYSet.length>0)){
			setWOI(0, 0, width, height);//  set WOI for the current image
			for (int n=0;n<pixelsXYSet.length;n++)
				for (int i=0;i<pixelsXYSet[n].length;i++) if ((pixelsXYSet[n][i]!=null)&&(pixelsUVSet[n][i]!=null)) {
					if (numCells==0){
						minV=pixelsUVSet[n][i][1];
						maxV=pixelsUVSet[n][i][1];
						minU=pixelsUVSet[n][i][0];
						maxU=pixelsUVSet[n][i][0];
					} else {
						if (minV>pixelsUVSet[n][i][1])  minV=pixelsUVSet[n][i][1];
						else if (maxV<pixelsUVSet[n][i][1])  maxV=pixelsUVSet[n][i][1];
						if (minU>pixelsUVSet[n][i][0])  minU=pixelsUVSet[n][i][0];
						else if (maxU<pixelsUVSet[n][i][0])  maxU=pixelsUVSet[n][i][0];
					}
					numCells++;
				}
			if (numCells>0) {
				// do not break black/white correspondence, always move by even number of cells
				if ((minU & 1)!=0 )minU--;
				if ((minV & 1)!=0 )minV--;
				this.minUV[0]=minU;
				this.minUV[1]=minV; // save shift to restore later
				this.PATTERN_GRID=setPatternGridArray(maxU-minU+1,maxV-minV+1);
				this.gridContrastBrightness=new double[4][this.PATTERN_GRID.length][this.PATTERN_GRID[0].length]; //{grid contrast, grid intensity red, grid intensity green, grid intensity blue}[v][u]
				for (int n=0;n<4;n++)
					for (int v=0;v<this.gridContrastBrightness[0].length;v++)
						for (int u=0;u<this.gridContrastBrightness[0][0].length;u++)
							this.gridContrastBrightness[n][v][u]=0.0;

				for (int n=0;n<pixelsXYSet.length;n++)
					for (int i=0;i<pixelsXYSet[n].length;i++) if ((pixelsXYSet[n][i]!=null)&&(pixelsUVSet[n][i]!=null)) {
						int [] shiftedUV={pixelsUVSet[n][i][0]-minU,pixelsUVSet[n][i][1]-minV};
						setPatternGridCell(
								this.PATTERN_GRID,
								shiftedUV, //pixelsUV[i],
								pixelsXYSet[n][i], // will have extra data
								null,
								null);
						this.gridContrastBrightness[0][shiftedUV[1]][shiftedUV[0]]=pixelsXYSet[n][i][2];
						this.gridContrastBrightness[1][shiftedUV[1]][shiftedUV[0]]=pixelsXYSet[n][i][3]*intensityRange[0];
						this.gridContrastBrightness[2][shiftedUV[1]][shiftedUV[0]]=pixelsXYSet[n][i][4]*intensityRange[1];
						this.gridContrastBrightness[3][shiftedUV[1]][shiftedUV[0]]=pixelsXYSet[n][i][5]*intensityRange[2];

					}
			}
		}
		return numCells;
	}
	/**
	 * Create PATTERN_GRID for calculated grid (for sensor parameters, orientation) for debugging purposes
	 * @param hintGrid  grid array [v][u][0- x,  1 - y, 2 - u, 3 - v] (u,v - not used here)
	 * @param width     image width, in pixels
	 * @param height    image height, in pixels
	 * @return          number of non-empty cells
	 */
	public int restoreSimulatedPatternGridFromHint(double [][][] hintGrid, int width, int height){
		int numCells=0;
		this.PATTERN_GRID=null;
		if (hintGrid==null) return 0;
		setWOI(0, 0, width, height);//  set WOI for the current image
		//        	   this.PATTERN_GRID=setPatternGridArray(hintGrid[0].length+1,hintGrid.length+1);
		this.PATTERN_GRID=setPatternGridArray(hintGrid[0].length,hintGrid.length);
		for (int v=0;v<hintGrid.length;v++) for (int u=0;u<hintGrid[v].length;u++) if (hintGrid[v][u]!=null){
			double [] xy={hintGrid[v][u][0],hintGrid[v][u][1]};
			if ((xy[0]>=0) && (xy[1]>=0) && (xy[0]<width) && (xy[1]<height)) {
				int    [] uv={u,v};
				//            		   if ((xy[0]==0) && (xy[1]==0)) {
				//            			   System.out.println("x==0,y==0 for u="+u+", v="+v+", hintGrid[v][u][2]="+hintGrid[v][u][2]+", hintGrid[v][u][3]="+hintGrid[v][u][3]);
				//            		   }
				setPatternGridCell(
						this.PATTERN_GRID,
						uv,
						xy,
						null,
						null);
				numCells++;
			}
		}
		return numCells;
	}


	/**
	 * restore grid parameters - this.PATTERN_GRID (no wave vectors) only
	 * absolute calibration (if any) is lost, only orientation is preserved
	 * @param imp_grid - grid encoded as image
	 * @return array of laser pointers coordinates - no separate, return number of grid cells
	 */
	public int restorePatternGridFromImage(ImagePlus imp_grid){
		int numCells=0;
		this.PATTERN_GRID=null;
		setWOI(0, 0, 0, 0);
		if (imp_grid!=null){
			setWOI(0, 0, imp_grid.getWidth(), imp_grid.getHeight());//  set WOI for the current image
			ImageStack stack=imp_grid.getStack();
			float [][] pixels=new float[4][];
			if ((stack==null) || (stack.getSize()!=4)) {
				String msg="Expected a 4-slice stack in "+imp_grid.getTitle();
				IJ.showMessage("Error",msg);
				throw new IllegalArgumentException (msg);
			}
			for (int i=0;i<4;i++) pixels[i]= (float[]) stack.getPixels(i+1); // pixel X : negative - no grid here
			int minU=0,minV=0,maxU=0,maxV=0;
			for (int i=0;i<pixels[0].length;i++) if (pixels[0][i]>=0) {
				int u=Math.round(pixels[2][i]);
				int v=Math.round(pixels[3][i]);
				if (numCells==0){
					minV=v;
					maxV=v;
					minU=u;
					maxU=u;
				} else {
					if (minV>v)  minV=v;
					else if (maxV<v)  maxV=v;
					if (minU>u)  minU=u;
					else if (maxU<u)  maxU=u;
				}
				numCells++;
			}
			if (numCells>0) {
				setPatternGridArray(maxU-minU+1,maxV-minV+1);
				double [] xy=new double[2];
				int [] uv=new int[2];
				for (int i=0;i<pixels[0].length;i++) if (pixels[0][i]>=0) {
					uv[0]=(Math.round(pixels[2][i]))-minU;
					uv[1]=(Math.round(pixels[3][i]))-minV;
					xy[0]=pixels[0][i];
					xy[1]=pixels[1][i];
					setPatternGridCell(
							this.PATTERN_GRID,
							uv,
							xy,
							null,
							null);
				}
			}
		}
		return numCells;
	}
	/**
	 * Calculate this.UV_INDEX (and this.UV_INDEX_WIDTH) - map from image pixel to U,V (U*this.UV_INDEX_WIDTH*V). this.DIST_SELECTION should be set
	 * @param imp source image (just to find image size, null - use this.DIST_SELECTION.width, this.DIST_SELECTION.height)
	 * @param shiftXY pattern shift (from debug), null - use {0,0}
	 * @param threadsMax limit on threads to use
	 * @param updateStatus update ImageJ status bar
	 * @param global_debug_level global debug level
	 * @param debug_level loop debug level
	 * @return true - OK, false - failure
	 */
	public boolean createUV_INDEX(
			ImagePlus imp, // or null - just to determine WOI (when getWOI matches image size)
			double [] shiftXY, // add to patterGrid xy, null OK
			int threadsMax,
			boolean updateStatus,
			int global_debug_level, // DEBUG_LEVEL
			int debug_level // debug level used inside loops
			){
		SimulationPattern simulationPattern=new SimulationPattern(); // do not need bitmap array here
		float [] UV_float0= simulationPattern.simulateGrid (
				getDArray(),
				2, // gridFrac, // number of grid steps per pattern full period
				null, //simulParameters,
				getWOI(),
				1, // SIMUL.subdiv/2,
				shiftXY,    // add to patterGrid xy, null OK
				threadsMax,
				updateStatus,
				debug_level); // debug level
		if (UV_float0==null) {
			System.out.println ("BUG: createUV_INDEX(): simulationPattern.simulateGrid() returnerd null");
			System.out.println ("BUG: createUV_INDEX(): getDArray() returnerd "+((getDArray()==null)?"null":"noy null"));
			return false;
		}
		float [] UV_float= simulationPattern.combineWithCanvas(
				-1.0,
				((imp==null)?(getWOI().width):imp.getWidth()),
				((imp==null)?(getWOI().height):imp.getHeight()),
				getWOI(),
				UV_float0 );
		if (global_debug_level>3) SDFA_INSTANCE.showArrays(UV_float0,getWOI().width, getWOI().height, "UV_float0"); // all -1
		if (global_debug_level>2) SDFA_INSTANCE.showArrays(UV_float, ((imp==null)?(getWOI().width):imp.getWidth()), ((imp==null)?(getWOI().height):imp.getHeight()), "UV_float");

		this.UV_INDEX=new int [UV_float.length];
		this.UV_INDEX_WIDTH=((imp==null)?(getWOI().width):imp.getWidth());
		for (int i=0;i<this.UV_INDEX.length;i++) this.UV_INDEX[i]=(int) UV_float[i];
		return true;
	}
	/**
	 * Apply hint grid and laser pointer calibration to the grid. If (hintGrid!= null) && (hintGridTolerance>0) than any result >=0 means match
	 * Create this.pixelsUV, this.targetUV
	 * @param laserPointer Laser pointer parameters DEPRECATE, move to pointersXYUV
	 * @param pointersXY   pairs of detected pointers x,y (or nulls)
	 * @param removeOutOfGridPointers if true - remove pointers if they are outside of the pattern grid
	 * @param hintGrid predicted grid array (or null)
	 * @param hintGridTolerance allowed mismatch (fraction of period) or 0 - orientation only
	 * @param global_debug_level debug level
	 * @param noMessageBoxes do not open (and wait for) dialog boxes
	 * @return >=0 - number of laser pointers used for calibration, <0 - different errors
	 */

	public int combineGridCalibration_deprecated(
			LaserPointer laserPointer, // LaserPointer object or null
			double [][] pointersXY,
			boolean removeOutOfGridPointers, //
			double [][][] hintGrid, // predicted grid array (or null)
			double        hintGridTolerance, // alllowed mismatch (fraction of period) or 0 - orientation only
			int global_debug_level, // DEBUG_LEVEL
			boolean noMessageBoxes
			){
		int acalibrated=0;
		double [][] gridMatchCoeff=null;
		double searchAround=20.0; // how far to look for the grid node
		int gridRotation=-1; //undefined
		int [] iGridTranslateUV=null; // translate UV grid by these integer numbers
		if (hintGrid!=null){
			gridMatchCoeff=calcGridMatchMatrix (hintGrid, searchAround);
			// now already rounds rotate and re-replaces gridMatchCoeff with approximated, refines gridMatchCoeff[0][2] and gridMatchCoeff[1][2]
			if (gridMatchCoeff!=null) {
				gridRotation=matrixToRot(gridMatchCoeff);
				this.debugLevel=global_debug_level;
				int [][] iGridMatchCoeff=gridMatrixApproximate(gridMatchCoeff);
				if (global_debug_level>1){
					System.out.println("gridMatchCoeff[0]={"+IJ.d2s(gridMatchCoeff[0][0],5)+", "+IJ.d2s(gridMatchCoeff[0][1],5)+", "+IJ.d2s(gridMatchCoeff[0][2],5)+"}");
					System.out.println("gridMatchCoeff[1]={"+IJ.d2s(gridMatchCoeff[1][0],5)+", "+IJ.d2s(gridMatchCoeff[1][1],5)+", "+IJ.d2s(gridMatchCoeff[1][2],5)+"}");
					System.out.println("gridRotation="+gridRotation);
					System.out.println("iGridMatchCoeff[0]={"+iGridMatchCoeff[0][0]+", "+iGridMatchCoeff[0][1]+", "+iGridMatchCoeff[0][2]+"}");
					System.out.println("iGridMatchCoeff[1]={"+iGridMatchCoeff[1][0]+", "+iGridMatchCoeff[1][1]+", "+iGridMatchCoeff[1][2]+"}");
					System.out.println("worstGridMatchRotSkew()="+IJ.d2s(worstGridMatchRotSkew(gridMatchCoeff),5));
					System.out.println("worstGridMatchTranslate()="+IJ.d2s(worstGridMatchTranslate(gridMatchCoeff),5));
				}
				// hintGridTolerance==0 - do not try to determine shift from the hint (not reliable yet)
				if (hintGridTolerance>0) {
					if (worstGridMatchTranslate(gridMatchCoeff)<=hintGridTolerance){ // convert to pixels from halfperiods (or just chnage definition of hintGridTolerance)
						if (global_debug_level>1) System.out.println("worstGridMatchTranslate(gridMatchCoeff)= "+worstGridMatchTranslate(gridMatchCoeff)+", hintGridTolerance="+hintGridTolerance);
						iGridTranslateUV=new int[2];
						iGridTranslateUV[0]=iGridMatchCoeff[0][2];
						iGridTranslateUV[1]=iGridMatchCoeff[1][2];
					} else {
						if (global_debug_level>1) System.out.println("*** Warning: combineGridCalibration() failed,  worstGridMatchTranslate(gridMatchCoeff)= "+worstGridMatchTranslate(gridMatchCoeff)+", hintGridTolerance="+hintGridTolerance);
						return -1;
					}
				}
				if (global_debug_level>0){
					System.out.println((((iGridMatchCoeff[0][2]+iGridMatchCoeff[1][2])&1)==0)?"EVEN shift":"ODD shift");
				}
			} else {
				if (global_debug_level>0) System.out.println("*** Warning: combineGridCalibration(): gridMatchCoeff() failed");
				return -1;
			}
		}
		if (((laserPointer!=null) && (laserPointer.laserUVMap.length>0)) ||
				((iGridTranslateUV!=null) && (gridRotation>=0))){ // no laser pointers, but hint grid with specified tolerance
			if ((global_debug_level>1) && (pointersXY==null)) System.out.println("This image does not contain any laser pointer data");
			acalibrated=calibrateGrid( // now should work without laser pointers too
					laserPointer,
					pointersXY,
					removeOutOfGridPointers,
					gridRotation,
					iGridTranslateUV,
					noMessageBoxes,
					global_debug_level);
			if (global_debug_level>1) {
				System.out.println("matchSimulatedPattern.laserCalibrateGrid() returned "+acalibrated+
						((acalibrated>0)?" laser points used.":(((iGridTranslateUV==null) || (acalibrated<0))?" - error code":"none")));
			}
		}
		if (global_debug_level>0) System.out.println("Pattern size is "+getDArrayWidth()+" x "+ getDArrayHeight());
		return acalibrated;

	}

// Modified 06/19 to move laser pointers to files
	public int combineGridCalibration(
			LaserPointer lp, // Only for possible hint on rotations/ flips.LaserPointer object or null
			double [][] pointersXYUV,
			boolean removeOutOfGridPointers, //
			double [][][] hintGrid, // predicted grid array (or null)
			double        hintGridTolerance, // alllowed mismatch (fraction of period) or 0 - orientation only
			int global_debug_level, // DEBUG_LEVEL
			boolean noMessageBoxes
			){
		boolean has_lasers = false;
		if (pointersXYUV != null) {
			for (double[] e:pointersXYUV) if ((e != null) && (e.length>2)) {
				has_lasers = true;
				break;
			}
		}
		int acalibrated=0;
		double [][] gridMatchCoeff=null;
		double searchAround=20.0; // how far to look for the grid node
		int gridRotation=-1; //undefined
		int [] iGridTranslateUV=null; // translate UV grid by these integer numbers
		if (hintGrid!=null){
			gridMatchCoeff=calcGridMatchMatrix (hintGrid, searchAround);
			// now already rounds rotate and re-replaces gridMatchCoeff with approximated, refines gridMatchCoeff[0][2] and gridMatchCoeff[1][2]
			if (gridMatchCoeff!=null) {
				gridRotation=matrixToRot(gridMatchCoeff);
				this.debugLevel=global_debug_level;
				int [][] iGridMatchCoeff=gridMatrixApproximate(gridMatchCoeff);
				if (global_debug_level>1){
					System.out.println("gridMatchCoeff[0]={"+IJ.d2s(gridMatchCoeff[0][0],5)+", "+IJ.d2s(gridMatchCoeff[0][1],5)+", "+IJ.d2s(gridMatchCoeff[0][2],5)+"}");
					System.out.println("gridMatchCoeff[1]={"+IJ.d2s(gridMatchCoeff[1][0],5)+", "+IJ.d2s(gridMatchCoeff[1][1],5)+", "+IJ.d2s(gridMatchCoeff[1][2],5)+"}");
					System.out.println("gridRotation="+gridRotation);
					System.out.println("iGridMatchCoeff[0]={"+iGridMatchCoeff[0][0]+", "+iGridMatchCoeff[0][1]+", "+iGridMatchCoeff[0][2]+"}");
					System.out.println("iGridMatchCoeff[1]={"+iGridMatchCoeff[1][0]+", "+iGridMatchCoeff[1][1]+", "+iGridMatchCoeff[1][2]+"}");
					System.out.println("worstGridMatchRotSkew()="+IJ.d2s(worstGridMatchRotSkew(gridMatchCoeff),5));
					System.out.println("worstGridMatchTranslate()="+IJ.d2s(worstGridMatchTranslate(gridMatchCoeff),5));
				}
				// hintGridTolerance==0 - do not try to determine shift from the hint (not reliable yet)
				if (hintGridTolerance>0) {
					if (worstGridMatchTranslate(gridMatchCoeff)<=hintGridTolerance){ // convert to pixels from halfperiods (or just chnage definition of hintGridTolerance)
						if (global_debug_level>1) System.out.println("worstGridMatchTranslate(gridMatchCoeff)= "+worstGridMatchTranslate(gridMatchCoeff)+", hintGridTolerance="+hintGridTolerance);
						iGridTranslateUV=new int[2];
						iGridTranslateUV[0]=iGridMatchCoeff[0][2];
						iGridTranslateUV[1]=iGridMatchCoeff[1][2];
					} else {
						if (global_debug_level>1) System.out.println("*** Warning: combineGridCalibration() failed,  worstGridMatchTranslate(gridMatchCoeff)= "+worstGridMatchTranslate(gridMatchCoeff)+", hintGridTolerance="+hintGridTolerance);
						return -1;
					}
				}
				if (global_debug_level>0){
					System.out.println((((iGridMatchCoeff[0][2]+iGridMatchCoeff[1][2])&1)==0)?"EVEN shift":"ODD shift");
				}
			} else {
				if (global_debug_level>0) System.out.println("*** Warning: combineGridCalibration(): gridMatchCoeff() failed");
				return -1;
			}
		}
		if (has_lasers ||
				((iGridTranslateUV!=null) && (gridRotation>=0))){ // no laser pointers, but hint grid with specified tolerance
			if ((global_debug_level>1) && (pointersXYUV==null)) System.out.println("This image does not contain any laser pointer data");
			double maxOffsetFromCenter=0.6; // maximal offset of the laser spot from the center, relative to cell radius
			boolean white_only = true;
			if (lp != null) {
				maxOffsetFromCenter = lp.maxOffsetFromCenter;
				white_only = lp.whiteOnly;
			}
			acalibrated=calibrateGrid( // now should work without laser pointers too
					lp, // for flips hint only
					white_only,
					maxOffsetFromCenter,
					pointersXYUV,
					removeOutOfGridPointers,
					gridRotation,
					iGridTranslateUV,
					noMessageBoxes,
					global_debug_level);
			if (global_debug_level>1) {
				System.out.println("matchSimulatedPattern.laserCalibrateGrid() returned "+acalibrated+
						((acalibrated>0)?" laser points used.":(((iGridTranslateUV==null) || (acalibrated<0))?" - error code":"none")));
			}
		}
		if (global_debug_level>0) System.out.println("Pattern size is "+getDArrayWidth()+" x "+ getDArrayHeight());
		return acalibrated;

	}




	public int replaceGridXYWithProjected(double [][][] projectedGrid, String debugTitle){
		int minU=0,minV=0,maxU=0,maxV=0;
		boolean notYetSet=true;
		for (double [][]row:projectedGrid) for (double [] cell:row) if (cell!=null){
			int u = (int) cell[2];
			int v = (int) cell[3];
			if (notYetSet){
				minU=u;
				maxU=u;
				minV=v;
				maxV=v;
				notYetSet=false;
			} else {
				if (minU>u) minU=u;
				if (maxU<u) maxU=u;
				if (minV>v) minV=v;
				if (maxV<v) maxV=v;
			}
		}
		double [][][] grid=new double [maxV-minV+1][maxU-minU+1][];
		//        	   for (double [][]row:grid) for (double [] cell:row) cell=null; // See if this works with "enhanced for loop"
		for (double [][]row:grid) for (int u=0;u<row.length;u++) row[u]=null; // See if this works with "enhanced for loop"

		for (double [][] row:projectedGrid) for (double [] cell:row) if (cell!=null){
			int u = (int) cell[2];
			int v = (int) cell[3];
			double [] xy={cell[0],cell[1]};
			grid[v-minV][u-minU]=xy;
		}
		int numNewDefined=0;
		//        	   System.out.println("this.PATTERN_GRID.length="+this.PATTERN_GRID.length+"this.PATTERN_GRID[0.length="+this.PATTERN_GRID[0].length);
		//        	   System.out.println("this.targetUV.length="+this.targetUV.length+"this.targetUV[0.length="+this.targetUV[0].length);
		if ((debugTitle!=null) && (this.debugLevel>0)){
			double [][] debugReplace=null;
			String [] debugTiltes={"deltaX","deltaY","Contrast", "measX","measY", "targetU","targetV"};
			debugReplace=new double[7][this.PATTERN_GRID.length*this.PATTERN_GRID[0].length];
			for (int i=0;i<debugReplace.length;i++) Arrays.fill(debugReplace[i],Double.NaN);
			for (int v=0;v<this.PATTERN_GRID.length;v++) for (int u=0;u<this.PATTERN_GRID[v].length;u++) {
				double [][] cell=this.PATTERN_GRID[v][u];
				if ((cell !=null) && (cell.length>0) &&(cell[0] !=null) && (cell[0].length>1)){
					int tu=this.targetUV[v][u][0]-minU;
					int tv=this.targetUV[v][u][1]-minV;
					if ((tu>=0) && (tv>=0) && (tv<grid.length) && (tu<grid[tv].length) && (grid[tv][tu]!=null)) {
						int index=v*this.PATTERN_GRID[0].length+u;
						debugReplace[0][index]=grid[tv][tu][0]-this.PATTERN_GRID[v][u][0][0];
						debugReplace[1][index]=grid[tv][tu][1]-this.PATTERN_GRID[v][u][0][1];
						debugReplace[2][index]=this.PATTERN_GRID[v][u][0][2];
						debugReplace[3][index]=this.PATTERN_GRID[v][u][0][0];
						debugReplace[4][index]=this.PATTERN_GRID[v][u][0][1];
						debugReplace[5][index]=tu;
						debugReplace[6][index]=tu;
					}
				}
			}
			SDFA_INSTANCE.showArrays(debugReplace, this.PATTERN_GRID[0].length, this.PATTERN_GRID.length, true, "replaceGridXYWithProjected-"+debugTitle, debugTiltes);
		}

		for (int v=0;v<this.PATTERN_GRID.length;v++) for (int u=0;u<this.PATTERN_GRID[v].length;u++) {
			double [][] cell=this.PATTERN_GRID[v][u];
			if ((cell !=null) && (cell.length>0) &&(cell[0] !=null) && (cell[0].length>1)){
				//                	   System.out.print("v="+v+" u="+u);
				int tu=this.targetUV[v][u][0]-minU;
				int tv=this.targetUV[v][u][1]-minV;
				//                	   System.out.println("  tv="+tv+" tu="+tu);
				if ((tu>=0) && (tv>=0) && (tv<grid.length) && (tu<grid[tv].length) && (grid[tv][tu]!=null)) {
					cell[0][0]=grid[tv][tu][0]; // -81 -.-1
					cell[0][1]=grid[tv][tu][1];
					if (Double.isNaN(cell[0][0]) || Double.isNaN(cell[0][1])){
						this.PATTERN_GRID[v][u]=null; // make it undefined
					} else {
						numNewDefined++;
					}
				} else {
					this.PATTERN_GRID[v][u]=null; // make it undefined
				}
			}
		}
		return numNewDefined;
	}

	/* Get calibrated pattern as a 8-slice image (can be saved as TIFF)
	 * first   slice - pixel X or -1 for undefined
	 * second  slice - pixel Y or -1 for undefined
	 * third   slice - target U (may be negative)
	 * fourth  slice - target V (may be negative)
	 *  other slices - if present
	 * fifth   slice -  local grid contrast (looks for 2  white and 2 blacks around) - can be used to filter
	 * sixth   slice - red intensity of the grid (averaged around the grid node)
	 * seventh slice - green intensity of the grid (averaged around the grid node)
	 * eighth  slice - blue intensity of the grid (averaged around the grid node)
	 */
	public ImagePlus getCalibratedPatternAsImage(String title, int numUsedPointers){
		if ((this.targetUV==null) ||(this.pXYUV==null)) {
//			System.out.println("getCalibratedPatternAsImage(): this.targetUV="+((this.targetUV==null)?"null":"not null")+", this.pixelsUV="+((this.pXYUV==null)?"null":"not null"));
//			System.out.println("Using grid w/o absolute calibration.");
			unCalibrateGrid();
			//        		   return null;
		}
		int numSlices=(this.gridContrastBrightness==null)?4:8;
		float [][] pixels=new float [numSlices][getWidth()*getHeight()];
		ImageStack stack=new ImageStack(getWidth(),getHeight());
		int index=0;
		for (int v=0;v<getHeight();v++) for (int u=0;u<getWidth();u++) {
			if ((this.targetUV[v][u]==null) ||(this.pXYUV[v][u]==null)||
					(this.pXYUV[v][u][0]<0.0) || (this.pXYUV[v][u][1]<0.0)) { // disregard negative sensor pixels
				pixels [0][index]=-1.0f;
				pixels [1][index]=-1.0f;
				pixels [2][index]= 0.0f;
				pixels [3][index]= 0.0f;
				if (numSlices>4){
					pixels [4][index]= 0.0f; // contrast
					pixels [5][index]=-1.0f; // red, undefined
					pixels [6][index]=-1.0f; // green, undefined
					pixels [7][index]=-1.0f; // blue, undefined

				}
			} else {
				pixels [0][index]=(float) this.pXYUV[v][u][0];
				pixels [1][index]=(float) this.pXYUV[v][u][1];
				pixels [2][index]=this.targetUV[v][u][0];
				pixels [3][index]=this.targetUV[v][u][1];
				if (numSlices>4){
					pixels [4][index]=(float) this.gridContrastBrightness[0][v][u]; // grid contrast
					pixels [5][index]=(float) this.gridContrastBrightness[1][v][u]; // red
					pixels [6][index]=(float) this.gridContrastBrightness[2][v][u]; // green
					pixels [7][index]=(float) this.gridContrastBrightness[3][v][u]; // blue
				}
			}
			index++;
		}
		stack.addSlice("pixel-X",  pixels[0]);
		stack.addSlice("pixel-Y",  pixels[1]);
		stack.addSlice("target-U", pixels[2]);
		stack.addSlice("target-V", pixels[3]);
		if (numSlices>4){
			stack.addSlice("contrast",  pixels[4]);
			stack.addSlice("red",       pixels[5]);
			stack.addSlice("green",     pixels[6]);
			stack.addSlice("blue",      pixels[7]);

		}
		ImagePlus imp = new ImagePlus(title, stack);
		imp.setProperty("USED_POINTERS",((numUsedPointers>=0)?numUsedPointers:0)+"");
		//				System.out.println("getCalibratedPatternAsImage(): numUsedPointers="+numUsedPointers+" getProperty(\"USED_POINTERS\")="+imp.getProperty("USED_POINTERS"));
		return imp;
	}
	// searching for single-pixel errors (program bug)
	public ImagePlus getCalibratedPatternCurvatureAsImage(String title){
		if ((this.targetUV==null) ||(this.pXYUV==null)) {
			String msg="this.targetUV="+((this.targetUV==null)?"null":"not null")+", this.pixelsUV="+((this.pXYUV==null)?"null":"not null");
			IJ.showMessage("Error",msg);
			throw new IllegalArgumentException (msg);
		}
		double [][][] curves=     new double [this.targetUV.length][this.targetUV[0].length][2];
		double [][][] diff_curves=new double [this.targetUV.length][this.targetUV[0].length][2];
		boolean [][]   mask_curves=new boolean [this.targetUV.length][this.targetUV[0].length];
		boolean [][]   mask_diff_curves=new boolean [this.targetUV.length][this.targetUV[0].length];
		for (int v=0;v<curves.length;v++) for (int u=0;u<curves[0].length;u++){
			curves[v][u][0]=0.0;
			curves[v][u][1]=0.0;
			diff_curves[v][u][0]=0.0;
			diff_curves[v][u][1]=0.0;
			mask_curves[v][u]=false;
			mask_diff_curves[v][u]=false;
		}
		ImageStack stack=new ImageStack(getWidth(),getHeight());
		int [][] dirs=   {{0,0},{-1,0},{0,-1},{1,0},{0,1},{-1,-1},{-1,1},{1,-1},{1,1}};
		double [] weights={1.0, -0.15,- 0.15, -.15, -.15,   -0.1,  -0.1,  -0.1, -0.1};
		// first pass - calculate "curvature" - difference between pixel dx, dy values and those average (weighted) of 8 neigbors
		for (int v=1;v<getHeight()-1;v++) for (int u=1;u<getWidth()-1;u++) {
			double [] avrg={0.0,0.0};
			boolean valid=true;
			for (int d=0;d<dirs.length;d++){
				int u1=u+dirs[d][0];
				int v1=v+dirs[d][1];
				if ((this.targetUV[v1][u1]==null) ||(this.pXYUV[v1][u1]==null)||
						(this.pXYUV[v1][u1][0]<0.0) || (this.pXYUV[v1][u1][0]<0.0)) { // disregard negative sensor pixels
					valid=false;
					break;
				} else {
					avrg[0]+=weights[d]*this.pXYUV[v1][u1][0];
					avrg[1]+=weights[d]*this.pXYUV[v1][u1][1];
				}
			}
			if (valid) {
				curves[v][u][0]=avrg[0];
				curves[v][u][1]=avrg[1];
				mask_curves[v][u]=true;
			}
		}
		// second pass - calculate difference between curvatives of each node and and those average (weighted) of 8 neigbors
		// This will mostly eliminate the distortion shift, and can be used as a measure of the "noise"
		for (int v=2;v<getHeight()-2;v++) for (int u=2;u<getWidth()-2;u++) {
			double [] avrg={0.0,0.0};
			boolean valid=true;
			for (int d=0;d<dirs.length;d++){
				int u1=u+dirs[d][0];
				int v1=v+dirs[d][1];
				if (!mask_curves[v1][u1]) {
					valid=false;
					break;
				} else {
					avrg[0]+=weights[d]*curves[v1][u1][0];
					avrg[1]+=weights[d]*curves[v1][u1][1];
				}
			}
			if (valid) {
				diff_curves[v][u][0]=avrg[0];
				diff_curves[v][u][1]=avrg[1];
				mask_diff_curves[v][u]=true;
			}
		}

		float [][] pixels=new float [7][getWidth()*getHeight()];
		int numPix=0;
		double sum=0.0;
		int index=0;
		double curvature, diff_curvature;
		for (int v=0;v<getHeight();v++) for (int u=0;u<getWidth();u++) {
			if (mask_curves[v][u]){
				curvature=Math.sqrt(curves[v][u][0]*curves[v][u][0]+curves[v][u][1]*curves[v][u][1]);
				pixels[0][index]=(float) curvature;
				pixels[1][index]=(float) curves[v][u][0];
				pixels[2][index]=(float) curves[v][u][1];
				if (mask_diff_curves[v][u]){
					diff_curvature=Math.sqrt(diff_curves[v][u][0]*diff_curves[v][u][0]+diff_curves[v][u][1]*diff_curves[v][u][1]);
					pixels[3][index]=(float) diff_curvature;
					pixels[4][index]=(float) diff_curves[v][u][0];
					pixels[5][index]=(float) diff_curves[v][u][1];
					pixels[6][index]=1.0f;
					sum+=diff_curvature*diff_curvature;
					numPix++;
				} else {
					pixels[3][index]=0.0f;
					pixels[4][index]=0.0f;
					pixels[5][index]=0.0f;
					pixels[6][index]=0.0f;
				}
			} else {
				pixels[0][index]=0.0f;
				pixels[1][index]=0.0f;
				pixels[2][index]=0.0f;

			}
			index++;
		}
		stack.addSlice("curvature",  pixels[0]);
		stack.addSlice("X-diff",     pixels[1]);
		stack.addSlice("Y-diff",     pixels[2]);
		stack.addSlice("error",      pixels[3]);
		stack.addSlice("X-error",    pixels[4]);
		stack.addSlice("Y-error",    pixels[5]);
		stack.addSlice("mask",       pixels[6]);
		if (numPix>0){
			double rms=Math.sqrt(sum/numPix);
			String msg="Deviation calculated for "+numPix+" grid nodes. RMS="+rms;
			System.out.println(msg);
			IJ.showMessage(msg);
		} else {
			String msg="Zero points to calculate deviation";
			System.out.println(msg);
			IJ.showMessage(msg);

		}
		ImagePlus imp = new ImagePlus(title, stack);
		return imp;
	}
	public ImagePlus getCalibratedPatternAsImage(
			ImagePlus imp_src,
			String prefix, int numUsedPointers){
		//        	   ImagePlus imp_result=getCalibratedPatternAsImage("grid-"+imp_src.getTitle(), numUsedPointers);
		ImagePlus imp_result=getCalibratedPatternAsImage(prefix+imp_src.getTitle(), numUsedPointers);
		if (imp_result == null) {
			System.out.println("getCalibratedPatternAsImage(): Grid is empty !");
			return null;
		}

		// copy all the properties to the new image
		JP46_Reader_camera jp4_instance= new JP46_Reader_camera(false);
		jp4_instance.copyProperties (imp_src,imp_result);
		jp4_instance.encodeProperiesToInfo(imp_result);
		return imp_result;
	}

	public boolean [] rotToFlips(int rot){
		boolean[][] rot2flips={      // swapUV,flipU,flipV for different rotations above
				{false,false,false},
				{true, true, false},
				{false,true, true },
				{true, false,true },
				{false,false,true },
				{true, false,false},
				{false,true, false},
				{true, true, true }};
		return rot2flips[rot];
	}
	public int flipsToRot(boolean swapUV, boolean flipU, boolean flipV) {
		for (int i=0;i<8;i++) if (
				(rotToFlips(i)[0]==swapUV) &&
				(rotToFlips(i)[1]==flipU) &&
				(rotToFlips(i)[2]==flipV)) return i;
		return -1; // never
	}
	//
	// move elsewhere?
	/**
	 * Create this.targetUV and this.pixelsUV for the grid that does not have any laser pointer references
	 */
	public void unCalibrateGrid(){
		// calculate targetUV that maps PATTERN_GRID cells to the target (absolute) UV
		this.targetUV=new int    [this.PATTERN_GRID.length][this.PATTERN_GRID[0].length][];
		this.pXYUV=new double [this.PATTERN_GRID.length][this.PATTERN_GRID[0].length][];
		// Or set it back to original 9do not touch, rotate/shift in the end?
		Arrays.fill(this.UVShiftRot, 0);
		for (int v=0;v<this.PATTERN_GRID.length;v++) for (int u=0;u<this.PATTERN_GRID[v].length;u++){
			if ((this.PATTERN_GRID[v][u]==null) || (this.PATTERN_GRID[v][u][0]==null)) {
				this.targetUV[v][u]=null;
				this.pXYUV[v][u]=null;
			} else {
				this.targetUV[v][u]=new int [2];
				this.targetUV[v][u][0]=u;
				this.targetUV[v][u][1]=v;
				this.pXYUV[v][u]=new double [2];
				this.pXYUV[v][u][0]=PATTERN_GRID[v][u][0][0];
				this.pXYUV[v][u][1]=PATTERN_GRID[v][u][0][1];
			}
		}

	}
	// returns -1 - failure, otherwise - number of points used for calibration array (move default orientation to laserPointer paramerters
	public int calibrateGrid(
			LaserPointer laserPointer,
			double [][] xy, // null and zero length OK
			boolean removeOutOfGridPointers,
			int hintRotation, // rotation (0..7) found from hintGrid, -1 - undefined
			int [] hintTranslateUV, // found from hintGrid: translate UV by this vector or null if undefined
			//				   double [][][] hintGrid, // predicted grid array (or null) - use just direction
			//				   double        hintGridTolerance, // alllowed mismatch (fraction of period) or 0 - orientation only
			boolean noMessageBoxes,
			int debugLevel
			)
	{
		if (xy==null) xy=new double[0][];
		invalidateCalibration();
		double [][]uv=uvFromXY(xy,removeOutOfGridPointers?2.0:-1);
		//        	   if (uv==null) return -1;
		int numPointesLeft=0;
		for (int i=0;i<xy.length;i++) if ((xy[i]!=null) && (uv[i]!=null)) numPointesLeft++;
		if (debugLevel>1){
			int numRemoved=0;
			for (int i=0;i<xy.length;i++) if ((xy[i]!=null) && (uv[i]==null)) numRemoved++;
			System.out.println("Removed "+numRemoved+" out-of-grid pointers, "+numPointesLeft+" pointers remain.");
		}
		// Now remove pointers that are not on white cells

		if ((laserPointer!=null) && laserPointer.whiteOnly){
			int numBad=0;
			for (int i=0;i<uv.length;i++) if (uv[i]!=null) {
				// Verify that laser spots are on the white cells (sum of uv is even)
				if ((((int)(Math.floor(uv[i][0])+Math.floor(uv[i][1]))) & 1)!=0){
					String msg="Laser point "+i+" is not on the white pattern cell, and this check is enforced in the configuration";
					System.out.println("Warning:"+msg);
					if (!noMessageBoxes) IJ.showMessage("Warning",msg);
					uv[i]=null;
					numBad++;
					continue;
				}
			}
			if (numBad>0){
				String msg="Removed "+numBad+" pointers on black cells";
				System.out.println("Warning:"+msg);
			}
		}
		// Later some pointers may be removed even if they are used to determine orientation/shift. But that should not lead
		// to white/black confusion
		/*
        	   int [][][] rotations={
        			   {{ 1, 0},{ 0, 1}}, // not mirrored
        			   {{ 0, 1},{-1, 0}},
        			   {{-1, 0},{ 0,-1}},
        			   {{ 0,-1},{ 1, 0}},

        			   {{ 1, 0},{ 0,-1}}, // mirrored
        			   {{ 0, 1},{ 1, 0}},
        			   {{-1, 0},{ 0, 1}},
        			   {{ 0,-1},{-1, 0}}};
        	   // shifts when rotating around unknown center (make it white)
        	   int [][] dfltShifts={
        			   {0,0},
        			   {0,1},
        			   {0,0},
        			   {1,0},
        			   {0,1},
        			   {0,0},
        			   {1,0},
        			   {0,0}};
		 */
		boolean [] possibleRotations={true,true,true,true,true,true,true,true};
		// If orientation is hinted, remove all other ones from the list of possible ones
		if (hintRotation>=0){ // defind from the hintGrid
			for (int i=0;i<possibleRotations.length;i++) possibleRotations[i]=(i==hintRotation);
		}
		//       	   boolean [] partialPossibleRotations=new boolean [possibleRotations.length];
		boolean pairMatch,allMatch;
		int [] diffUVTable=new int [2]; // difference between points specified in the table
		int [] diffUVMeas= new int [2]; // measured difference (PATTERN_GRID U,V
		int [] rotUVTable=  new int [2]; // rotated 'laser' coordinates difference (should match measured)
		int [] belongsToGoodPair=new int[uv.length];
		int [] belongsToBadPair=new int[uv.length];
		for (int i=0;i<uv.length;i++){
			belongsToGoodPair[i]=0;
			belongsToBadPair[i]=0;
		}
		// pass 0 - process good/bad pairs, do not disable directions if does not match
		// if at least 1 good pair exists - remove all that do not match
		// if no good pairs - remove all bad
		// second pass: if more than 1 good pair - should match all (or error)


		//TODO: When hinted position, remove far pointers before matching pairs

		for (int pass=0;pass<2;pass++) {
			for (int i=0;i<uv.length;i++) if (uv[i]!=null) for (int j=i+1;j<uv.length;j++) if (uv[j]!=null) {
				pairMatch=false;
				allMatch=false;
				diffUVTable[0]=(int) Math.round(laserPointer.laserUVMap[j][0]-laserPointer.laserUVMap[i][0]); // should not get here if uv is {}
				diffUVTable[1]=(int) Math.round(laserPointer.laserUVMap[j][1]-laserPointer.laserUVMap[i][1]);
				diffUVMeas[0]= (int) Math.round(uv[j][0]-uv[i][0]);
				diffUVMeas[1]= (int) Math.round(uv[j][1]-uv[i][1]);
				// see which rotations are possible for this pair of points
				if (debugLevel>2){
					System.out.println("pass="+pass+" i="+i+" j="+j);
					System.out.println("diffUVTable=["+diffUVTable[0]+","+diffUVTable[1]+"]");
					System.out.println("diffUVMeas= ["+diffUVMeas[0]+ ","+diffUVMeas[1]+"]");
				}
				for (int dir=0; dir<rotations.length;dir++) {
					rotUVTable[0]=rotations[dir][0][0]*diffUVTable[0]+rotations[dir][0][1]*diffUVTable[1];
					rotUVTable[1]=rotations[dir][1][0]*diffUVTable[0]+rotations[dir][1][1]*diffUVTable[1];
					if (debugLevel>2){
						System.out.println(dir+": rotUVTable= ["+rotUVTable[0]+ ","+rotUVTable[1]+"]");
					}
					if ((rotUVTable[0]==diffUVMeas[0]) && (rotUVTable[1]==diffUVMeas[1])) {
						pairMatch=true;
						if (possibleRotations[dir]) allMatch=true; // this and hinted direction
					} else {
						if (pass>0){ // do not disable rotation on the first pass
							possibleRotations[dir]=false;
						}
					}
				}
				// TODO: Find maximal number of matching pointers?

				if (pass==0) {
					//        				   if (allMatch) {
					if (pairMatch) {
						belongsToGoodPair[i]++;
						belongsToGoodPair[j]++;
					} else {
						belongsToBadPair[i]++;
						belongsToBadPair[j]++;
					}
				} else { // second pass
					if (!pairMatch) {
						String msg="Laser points\n"+i+" ["+IJ.d2s(laserPointer.laserUVMap[i][0],2)+":"+IJ.d2s(laserPointer.laserUVMap[i][1],2)+
								"] -> ["+IJ.d2s(uv[i][0],2)+":"+IJ.d2s(uv[i][1],2)+"] and \n"+
								j +" ["+IJ.d2s(laserPointer.laserUVMap[j][0],2)+":"+IJ.d2s(laserPointer.laserUVMap[j][1],2)+
								"] -> ["+IJ.d2s(uv[j][0],2)+":"+IJ.d2s(uv[j][1],2)+"] do not match";
						System.out.println(msg);
						if (!noMessageBoxes) IJ.showMessage("Error",msg);
						unCalibrateGrid();
						return -2;
					} else if (!allMatch){
						String msg="Following laser pointers can not be mapped simultaneously:\n";
						for (int k=0;k<=j;k++) if (uv[k]!=null) {
							msg+=k+" ["+IJ.d2s(laserPointer.laserUVMap[k][0],2)+":"+IJ.d2s(laserPointer.laserUVMap[k][1],2)+
									"] -> ["+IJ.d2s(uv[k][0],2)+":"+IJ.d2s(uv[k][1],2)+"]\n";
						}
						System.out.println(msg);
						if (!noMessageBoxes) IJ.showMessage("Error",msg);
						unCalibrateGrid();
						return -1;
					}
				}
			}
			if (pass==0){
				int numInGood=0;
				int numInBad=0;
				for (int i=0;i<uv.length;i++){
					if (belongsToGoodPair[i]>0) numInGood++;
					if (belongsToBadPair[i]>0) numInBad++;
				}
				if (numInBad>0){
					if (numInGood==0){
						String msg="No matching laser points pairs exist, and "+numInBad+" points do not match";
						System.out.println(msg);
						if (!noMessageBoxes) IJ.showMessage("Error",msg);
						/// will report error on the second pass
					} else if (numInBad>0){
						String msg="Matching laser points pair(s) exist(s), but other:";
						for (int i=0;i<uv.length;i++) if ((belongsToBadPair[i]>0) && (belongsToGoodPair[i]==0)){
							msg+=" #"+i+" ("+(i+1)+" of "+uv.length+")";
							uv[i]=null; // remove it from consideration
						}
						msg+=" do not match and will be removed.";
						System.out.println(msg);
						if (!noMessageBoxes) IJ.showMessage("Error",msg);
					}
				}
			}
		}
		//TODO: here at least some rotations match all points. If there ere more than two - try to use closest to the default/previous
		int rotation=(laserPointer!=null)?(flipsToRot(laserPointer.swapUV,laserPointer.flipU,laserPointer.flipV)):0;
		if (!possibleRotations[rotation]) { // current rotation value defined by laserPointer.{swapUV,flipU,flipV} does not match
			// find a new one (first - without mirroring)
			for (int i=0; i<8;i++) if (possibleRotations[(((rotation ^ i) & 4)) | ((rotation+i) & 3)]) {
				rotation=(((rotation ^ i) & 4)) | ((rotation+i) & 3); // first tried in the same half, then - the next one
				break;
			}
			if (!possibleRotations[rotation]) { // Program bug - should not happen
				String msg="Program error - could not find laser point mapping while it should exist\n";
				System.out.println(msg);
				if (!noMessageBoxes) IJ.showMessage("Error",msg);
				unCalibrateGrid();
				return -3;
			}
		}
		// now rotation is the correct one, update laserPointer.{swapUV,flipU,flipV};
		if (laserPointer!=null){
			laserPointer.swapUV=rotToFlips(rotation)[0];
			laserPointer.flipU =rotToFlips(rotation)[1];
			laserPointer.flipV =rotToFlips(rotation)[2];
		}
		//calculate shift
		int [] uvShift=dfltShifts[rotation].clone(); //{0,0};
		for (int i=0;i<uv.length;i++) if (uv[i]!=null) { // laserPointer -> uv=={}
			uvShift[0]=(int) Math.round(uv[i][0]-
					(rotations[rotation][0][0]*laserPointer.laserUVMap[i][0]+
							rotations[rotation][0][1]*laserPointer.laserUVMap[i][1]));
			uvShift[1]=(int) Math.round(uv[i][1]-
					(rotations[rotation][1][0]*laserPointer.laserUVMap[i][0]+
							rotations[rotation][1][1]*laserPointer.laserUVMap[i][1]));
			break;
		}

		// Hinted shift will only be used if no laser pointers are available, otherwise - only verify/warn
		if (hintTranslateUV!=null) {
			//        		   if ((uv.length==0) || (numPointesLeft==0)){
			if (numPointesLeft==0){
				uvShift[0]=hintTranslateUV[0];
				uvShift[1]=hintTranslateUV[1];
				if (debugLevel>1){
					System.out.println("No laser pointers available, using hinted translation");
				}
			} else {
				if ((uvShift[0]==hintTranslateUV[0]) && (uvShift[0]==hintTranslateUV[0])){
					if (debugLevel>1){
						System.out.println("Translation from the laser pointers matches the hinted one");
					}
				} else {
					if (debugLevel>1){
						System.out.println("Translation from the laser pointers does not match the hinted one:");
						System.out.println("Hinted: delta U="+hintTranslateUV[0]+", V="+hintTranslateUV[1]);
						System.out.println("Lasers: delta U="+uvShift[0]+", V="+uvShift[1]);
						System.out.println("Trusting lasers");
					}
				}
			}
		}
		// calculate remap array (rotation+translation) from the target UV to the measured grid UV.
		this.reMap=new int[2][3]; // seems it is never used?
		this.reMap[0][0]= rotations[rotation][0][0];
		this.reMap[0][1]= rotations[rotation][0][1];
		this.reMap[0][2]= uvShift[0];
		this.reMap[1][0]= rotations[rotation][1][0];
		this.reMap[1][1]= rotations[rotation][1][1];
		this.reMap[1][2]= uvShift[1];
		// calculate reverse remap array (rotation+translation) from the the measured grid UV to the target UV
		int reRot=(rotation>=4)?rotation:((4-rotation) & 3); // number of reverse mirror-rotation
		//        	   int [] UVRot={uvShift[0],uvShift[1],reRot};
		int [] UVRot={
				-(rotations[reRot][0][0]*uvShift[0] +rotations[reRot][0][1]*uvShift[1]),
				-(rotations[reRot][1][0]*uvShift[0] +rotations[reRot][1][1]*uvShift[1]),
				reRot};
		return  applyUVShiftRot(
				UVRot, // int [] UVShiftRot,
				uv,   // double [][]uv,
				laserPointer,
				noMessageBoxes);
		/*
        	   int [][] reReMap={
        			   {rotations[reRot][0][0],rotations[reRot][0][1],
        		       -(rotations[reRot][0][0]*uvShift[0] +rotations[reRot][0][1]*uvShift[1])},
        		       {rotations[reRot][1][0],rotations[reRot][1][1],
            		       -(rotations[reRot][1][0]*uvShift[0] +rotations[reRot][1][1]*uvShift[1])}};

        	   if (debugLevel>1){
        		   System.out.println("rotation="+rotation+", reMap= [["+this.reMap[0][0]+","+this.reMap[0][1]+","+this.reMap[0][2]+"]["+
        				   +this.reMap[1][0]+","+this.reMap[1][1]+","+this.reMap[1][2]+"]]");
        		   System.out.println("reRot="+   reRot+",  reReMap= [["+reReMap[0][0]+","+reReMap[0][1]+","+reReMap[0][2]+"]["+
        				   +reReMap[1][0]+","+reReMap[1][1]+","+reReMap[1][2]+"]]");
        	   }
// calculate targetUV that maps PATTERN_GRID cells to the target (absolute) UV
        	   this.targetUV=new int    [this.PATTERN_GRID.length][this.PATTERN_GRID[0].length][];
        	   this.pXYUV=new double [this.PATTERN_GRID.length][this.PATTERN_GRID[0].length][];
        	   for (int v=0;v<this.PATTERN_GRID.length;v++) for (int u=0;u<this.PATTERN_GRID[v].length;u++){
        		   if ((this.PATTERN_GRID[v][u]==null) || (this.PATTERN_GRID[v][u][0]==null)) {
        			   this.targetUV[v][u]=null;
        			   this.pXYUV[v][u]=null;
        		   } else {
        			   this.targetUV[v][u]=new int [2];
        			   this.targetUV[v][u][0]=reReMap[0][0]*u+reReMap[0][1]*v+reReMap[0][2];
        			   this.targetUV[v][u][1]=reReMap[1][0]*u+reReMap[1][1]*v+reReMap[1][2];
//        			   System.out.println("v="+v+", u="+u+", PATTERN_GRID.length="+PATTERN_GRID.length+", PATTERN_GRID[v].length="+PATTERN_GRID[v].length);
//        			   System.out.println("this.pixelsUV.length="+this.pixelsUV.length);
//        			   System.out.println("this.pixelsUV["+v+"].length="+this.pixelsUV[v].length);
        			   this.pXYUV[v][u]=new double [2];
        			   this.pXYUV[v][u][0]=PATTERN_GRID[v][u][0][0];
        			   this.pXYUV[v][u][1]=PATTERN_GRID[v][u][0][1];
        		   }
        	   }
        	   int numGood=0;
        	   int numBad=0;
        	   double [] distUV=new double[2];
        	   double dist;
        	   for (int i=0;i<uv.length;i++) if (uv[i]!=null) { //laserPointer == null > uv={}
        		// Verify that laser spots are inside specified distance from the cell centers
        		   distUV[0]=reReMap[0][0]*uv[i][0]+reReMap[0][1]*uv[i][1]+reReMap[0][2]-laserPointer.laserUVMap[i][0];
        		   distUV[1]=reReMap[1][0]*uv[i][0]+reReMap[1][1]*uv[i][1]+reReMap[1][2]-laserPointer.laserUVMap[i][1];
        		   dist=Math.sqrt(distUV[0]*distUV[0]+distUV[1]*distUV[1]);
            	   if (debugLevel>1){
            		   System.out.println("Laser spot #"+i+", distance from predicted ="+ IJ.d2s(dist,3)+" ("+IJ.d2s(200*dist,3)+
            				   "% of cell radius), du="+IJ.d2s(distUV[0],3)+", dv="+IJ.d2s(distUV[1],3));
            	   }
            	   if ((2*dist)> laserPointer.maxOffsetFromCenter){
        			   String msg="Laser point "+(i+1)+"(of "+uv.length+") is too far from the specified location, and this check is enforced in the configuration\n"+
        			              "measured distance is "+ IJ.d2s(200*dist,1)+"% of the cell radius, specified is "+ IJ.d2s(100*laserPointer.maxOffsetFromCenter,1)+"%";
        			   System.out.println("Warning:"+msg);
        			   if (!noMessageBoxes) IJ.showMessage("Warning",msg);
        			   numBad++;
        			   uv[i]=null;
        			   continue;
            	   }
        		   numGood++;
        	   }
        	   if ((debugLevel>0) && (numBad>0)){
        		   System.out.println("Removed "+numBad+" pointers that are too far from the predicted locations");
        	   }
        	   return numGood;
		 */
	}


	public int calibrateGrid(
			LaserPointer lp, // only as hint for rotations/flips (may update), null is OK
			boolean white_only, // laser pointer only on white
			double maxOffsetFromCenter,
			double [][] xyuv, // null and zero length OK now combines x,y and laserPointer.laserUVMap u,v
			// each non- null xyuv[] should be either 2 or 4 long
			boolean removeOutOfGridPointers,
			int hintRotation, // rotation (0..7) found from hintGrid, -1 - undefined
			int [] hintTranslateUV, // found from hintGrid: translate UV by this vector or null if undefined
			//				   double [][][] hintGrid, // predicted grid array (or null) - use just direction
			//				   double        hintGridTolerance, // alllowed mismatch (fraction of period) or 0 - orientation only
			boolean noMessageBoxes,
			int debugLevel
			)
	{
		if (xyuv==null) xyuv=new double[0][];
		invalidateCalibration();
		boolean has_lasers = false;
		for (double[] e:xyuv) if ((e != null) && (e.length>2)) {
			has_lasers = true;
			break;
		}
		double [][]uv=uvFromXY(xyuv,removeOutOfGridPointers?2.0:-1);
		//        	   if (uv==null) return -1;
		int numPointesLeft=0;
		for (int i=0;i<xyuv.length;i++) if ((xyuv[i]!=null) && (uv[i]!=null)) numPointesLeft++;
		if (debugLevel>1){
			int numRemoved=0;
			for (int i=0;i<xyuv.length;i++) if ((xyuv[i]!=null) && (uv[i]==null)) numRemoved++;
			System.out.println("Removed "+numRemoved+" out-of-grid pointers, "+numPointesLeft+" pointers remain.");
		}
		// Now remove pointers that are not on white cells

		if (has_lasers && white_only){
			int numBad=0;
			for (int i=0;i<uv.length;i++) if (uv[i]!=null) {
				// Verify that laser spots are on the white cells (sum of uv is even)
				if ((((int)(Math.floor(uv[i][0])+Math.floor(uv[i][1]))) & 1)!=0){
					String msg="Laser point "+i+" is not on the white pattern cell, and this check is enforced in the configuration";
					System.out.println("Warning:"+msg);
					if (!noMessageBoxes) IJ.showMessage("Warning",msg);
					uv[i]=null;
					numBad++;
					continue;
				}
			}
			if (numBad>0){
				String msg="Removed "+numBad+" pointers on black cells";
				System.out.println("Warning:"+msg);
			}
		}
		// Later some pointers may be removed even if they are used to determine orientation/shift. But that should not lead
		// to white/black confusion
		/*
        	   int [][][] rotations={
        			   {{ 1, 0},{ 0, 1}}, // not mirrored
        			   {{ 0, 1},{-1, 0}},
        			   {{-1, 0},{ 0,-1}},
        			   {{ 0,-1},{ 1, 0}},

        			   {{ 1, 0},{ 0,-1}}, // mirrored
        			   {{ 0, 1},{ 1, 0}},
        			   {{-1, 0},{ 0, 1}},
        			   {{ 0,-1},{-1, 0}}};
        	   // shifts when rotating around unknown center (make it white)
        	   int [][] dfltShifts={
        			   {0,0},
        			   {0,1},
        			   {0,0},
        			   {1,0},
        			   {0,1},
        			   {0,0},
        			   {1,0},
        			   {0,0}};
		 */
		boolean [] possibleRotations={true,true,true,true,true,true,true,true};
		// If orientation is hinted, remove all other ones from the list of possible ones
		if (hintRotation>=0){ // defind from the hintGrid
			for (int i=0;i<possibleRotations.length;i++) possibleRotations[i]=(i==hintRotation);
		}
		//       	   boolean [] partialPossibleRotations=new boolean [possibleRotations.length];
		boolean pairMatch,allMatch;
		int [] diffUVTable=new int [2]; // difference between points specified in the table
		int [] diffUVMeas= new int [2]; // measured difference (PATTERN_GRID U,V
		int [] rotUVTable=  new int [2]; // rotated 'laser' coordinates difference (should match measured)
		int [] belongsToGoodPair=new int[uv.length];
		int [] belongsToBadPair=new int[uv.length];
		for (int i=0;i<uv.length;i++){
			belongsToGoodPair[i]=0;
			belongsToBadPair[i]=0;
		}
		// pass 0 - process good/bad pairs, do not disable directions if does not match
		// if at least 1 good pair exists - remove all that do not match
		// if no good pairs - remove all bad
		// second pass: if more than 1 good pair - should match all (or error)


		//TODO: When hinted position, remove far pointers before matching pairs

		for (int pass=0;pass<2;pass++) {
			for (int i=0;i<uv.length;i++) if (uv[i]!=null) for (int j=i+1;j<uv.length;j++) if (uv[j]!=null) { // xyuv[i] != null, xyuv[j] != null
				pairMatch=false;
				allMatch=false;
				diffUVTable[0]=(int) Math.round(xyuv[j][2]-xyuv[i][2]); // should not get here if uv is {}
				diffUVTable[1]=(int) Math.round(xyuv[j][3]-xyuv[i][3]);
				diffUVMeas[0]= (int) Math.round(uv[j][0]-uv[i][0]);
				diffUVMeas[1]= (int) Math.round(uv[j][1]-uv[i][1]);
				// see which rotations are possible for this pair of points
				if (debugLevel>2){
					System.out.println("pass="+pass+" i="+i+" j="+j);
					System.out.println("diffUVTable=["+diffUVTable[0]+","+diffUVTable[1]+"]");
					System.out.println("diffUVMeas= ["+diffUVMeas[0]+ ","+diffUVMeas[1]+"]");
				}
				for (int dir=0; dir<rotations.length;dir++) {
					rotUVTable[0]=rotations[dir][0][0]*diffUVTable[0]+rotations[dir][0][1]*diffUVTable[1];
					rotUVTable[1]=rotations[dir][1][0]*diffUVTable[0]+rotations[dir][1][1]*diffUVTable[1];
					if (debugLevel>2){
						System.out.println(dir+": rotUVTable= ["+rotUVTable[0]+ ","+rotUVTable[1]+"]");
					}
					if ((rotUVTable[0]==diffUVMeas[0]) && (rotUVTable[1]==diffUVMeas[1])) {
						pairMatch=true;
						if (possibleRotations[dir]) allMatch=true; // this and hinted direction
					} else {
						if (pass>0){ // do not disable rotation on the first pass
							possibleRotations[dir]=false;
						}
					}
				}
				// TODO: Find maximal number of matching pointers?

				if (pass==0) {
					//        				   if (allMatch) {
					if (pairMatch) {
						belongsToGoodPair[i]++;
						belongsToGoodPair[j]++;
					} else {
						belongsToBadPair[i]++;
						belongsToBadPair[j]++;
					}
				} else { // second pass
					if (!pairMatch) {
						String msg="Laser points\n"+i+" ["+IJ.d2s(xyuv[i][2],2)+":"+IJ.d2s(xyuv[i][3],2)+
								"] -> ["+IJ.d2s(uv[i][0],2)+":"+IJ.d2s(uv[i][1],2)+"] and \n"+
								j +" ["+IJ.d2s(xyuv[j][2],2)+":"+IJ.d2s(xyuv[j][3],2)+
								"] -> ["+IJ.d2s(uv[j][0],2)+":"+IJ.d2s(uv[j][1],2)+"] do not match";
						System.out.println(msg);
						if (!noMessageBoxes) IJ.showMessage("Error",msg);
						unCalibrateGrid();
						return -2;
					} else if (!allMatch){
						String msg="Following laser pointers can not be mapped simultaneously:\n";
						for (int k=0;k<=j;k++) if (uv[k]!=null) {
							msg+=k+" ["+IJ.d2s(xyuv[k][2],2)+":"+IJ.d2s(xyuv[k][3],2)+
									"] -> ["+IJ.d2s(uv[k][0],2)+":"+IJ.d2s(uv[k][1],2)+"]\n";
						}
						System.out.println(msg);
						if (!noMessageBoxes) IJ.showMessage("Error",msg);
						unCalibrateGrid();
						return -1;
					}
				}
			}
			if (pass==0){
				int numInGood=0;
				int numInBad=0;
				for (int i=0;i<uv.length;i++){
					if (belongsToGoodPair[i]>0) numInGood++;
					if (belongsToBadPair[i]>0) numInBad++;
				}
				if (numInBad>0){
					if (numInGood==0){
						String msg="No matching laser points pairs exist, and "+numInBad+" points do not match";
						System.out.println(msg);
						if (!noMessageBoxes) IJ.showMessage("Error",msg);
						/// will report error on the second pass
					} else if (numInBad>0){
						String msg="Matching laser points pair(s) exist(s), but other:";
						for (int i=0;i<uv.length;i++) if ((belongsToBadPair[i]>0) && (belongsToGoodPair[i]==0)){
							msg+=" #"+i+" ("+(i+1)+" of "+uv.length+")";
							uv[i]=null; // remove it from consideration
						}
						msg+=" do not match and will be removed.";
						System.out.println(msg);
						if (!noMessageBoxes) IJ.showMessage("Error",msg);
					}
				}
			}
		}
		//TODO: here at least some rotations match all points. If there ere more than two - try to use closest to the default/previous
		int rotation=(lp != null)?(flipsToRot(lp.swapUV, lp.flipU, lp.flipV)):0;
		if (!possibleRotations[rotation]) { // current rotation value defined by laserPointer.{swapUV,flipU,flipV} does not match
			// find a new one (first - without mirroring)
			for (int i=0; i<8;i++) if (possibleRotations[(((rotation ^ i) & 4)) | ((rotation+i) & 3)]) {
				rotation=(((rotation ^ i) & 4)) | ((rotation+i) & 3); // first tried in the same half, then - the next one
				break;
			}
			if (!possibleRotations[rotation]) { // Program bug - should not happen
				String msg="Program error - could not find laser point mapping while it should exist\n";
				System.out.println(msg);
				if (!noMessageBoxes) IJ.showMessage("Error",msg);
				unCalibrateGrid();
				return -3;
			}
		}
		// now rotation is the correct one, update laserPointer.{swapUV,flipU,flipV};
		if (lp != null){
			lp.swapUV=rotToFlips(rotation)[0];
			lp.flipU =rotToFlips(rotation)[1];
			lp.flipV =rotToFlips(rotation)[2];
		}
		//calculate shift
		int [] uvShift=dfltShifts[rotation].clone(); //{0,0};
		for (int i=0;i<uv.length;i++) if (uv[i]!=null) { // laserPointer -> uv=={}
			uvShift[0]=(int) Math.round(uv[i][0]-
					(rotations[rotation][0][0]*xyuv[i][2]+
							rotations[rotation][0][1]*xyuv[i][3]));
			uvShift[1]=(int) Math.round(uv[i][1]-
					(rotations[rotation][1][0]*xyuv[i][2]+
							rotations[rotation][1][1]*xyuv[i][3]));
			break;
		}

		// Hinted shift will only be used if no laser pointers are available, otherwise - only verify/warn
		if (hintTranslateUV!=null) {
			//        		   if ((uv.length==0) || (numPointesLeft==0)){
			if (numPointesLeft==0){
				uvShift[0]=hintTranslateUV[0];
				uvShift[1]=hintTranslateUV[1];
				if (debugLevel>1){
					System.out.println("No laser pointers available, using hinted translation");
				}
			} else {
				if ((uvShift[0]==hintTranslateUV[0]) && (uvShift[0]==hintTranslateUV[0])){
					if (debugLevel>1){
						System.out.println("Translation from the laser pointers matches the hinted one");
					}
				} else {
					if (debugLevel>1){
						System.out.println("Translation from the laser pointers does not match the hinted one:");
						System.out.println("Hinted: delta U="+hintTranslateUV[0]+", V="+hintTranslateUV[1]);
						System.out.println("Lasers: delta U="+uvShift[0]+", V="+uvShift[1]);
						System.out.println("Trusting lasers");
					}
				}
			}
		}
		// calculate remap array (rotation+translation) from the target UV to the measured grid UV.
		this.reMap=new int[2][3]; // seems it is never used?
		this.reMap[0][0]= rotations[rotation][0][0];
		this.reMap[0][1]= rotations[rotation][0][1];
		this.reMap[0][2]= uvShift[0];
		this.reMap[1][0]= rotations[rotation][1][0];
		this.reMap[1][1]= rotations[rotation][1][1];
		this.reMap[1][2]= uvShift[1];
		// calculate reverse remap array (rotation+translation) from the the measured grid UV to the target UV
		int reRot=(rotation>=4)?rotation:((4-rotation) & 3); // number of reverse mirror-rotation
		//        	   int [] UVRot={uvShift[0],uvShift[1],reRot};
		int [] UVRot={
				-(rotations[reRot][0][0]*uvShift[0] +rotations[reRot][0][1]*uvShift[1]),
				-(rotations[reRot][1][0]*uvShift[0] +rotations[reRot][1][1]*uvShift[1]),
				reRot};
//		return  applyUVShiftRot(
//				UVRot, // int [] UVShiftRot,
//				uv,   // double [][]uv,
//				laserPointer,
//				noMessageBoxes);
		return  applyUVShiftRot(
				UVRot, // int [] UVShiftRot,
				uv, // double [][] uv,
				xyuv, // double [][] xyuv, // [][2], [][3] contain laser pointers u,v
				maxOffsetFromCenter, // double maxOffsetFromCenter,
//				LaserPointer laserPointer,
				noMessageBoxes);
	}

	public int applyUVShiftRot(
			int [] UVShiftRot,
			double [][]uv,
			LaserPointer laserPointer,
			boolean noMessageBoxes
			){
		if (UVShiftRot!=null) this.UVShiftRot=UVShiftRot.clone();
		int [][] reReMap=getRemapMatrix(UVShiftRot);
		if (debugLevel>1){
			//        		   System.out.println("rotation="+rotation+", reMap= [["+this.reMap[0][0]+","+this.reMap[0][1]+","+this.reMap[0][2]+"]["+
			//        				   +this.reMap[1][0]+","+this.reMap[1][1]+","+this.reMap[1][2]+"]]");
			System.out.println("reRot="+   UVShiftRot[2]+",  reReMap= [["+reReMap[0][0]+","+reReMap[0][1]+","+reReMap[0][2]+"]["+
					+reReMap[1][0]+","+reReMap[1][1]+","+reReMap[1][2]+"]]");
		}
		// calculate targetUV that maps PATTERN_GRID cells to the target (absolute) UV
		this.targetUV=new int    [this.PATTERN_GRID.length][this.PATTERN_GRID[0].length][];
		this.pXYUV=new double [this.PATTERN_GRID.length][this.PATTERN_GRID[0].length][];
		for (int v=0;v<this.PATTERN_GRID.length;v++) for (int u=0;u<this.PATTERN_GRID[v].length;u++){
			if ((this.PATTERN_GRID[v][u]==null) || (this.PATTERN_GRID[v][u][0]==null)) {
				this.targetUV[v][u]=null;
				this.pXYUV[v][u]=null;
			} else {
				this.targetUV[v][u]=new int [2];
				this.targetUV[v][u][0]=reReMap[0][0]*u+reReMap[0][1]*v+reReMap[0][2];
				this.targetUV[v][u][1]=reReMap[1][0]*u+reReMap[1][1]*v+reReMap[1][2];
				//        			   System.out.println("v="+v+", u="+u+", PATTERN_GRID.length="+PATTERN_GRID.length+", PATTERN_GRID[v].length="+PATTERN_GRID[v].length);
				//        			   System.out.println("this.pixelsUV.length="+this.pixelsUV.length);
				//        			   System.out.println("this.pixelsUV["+v+"].length="+this.pixelsUV[v].length);
				this.pXYUV[v][u]=new double [2];
				this.pXYUV[v][u][0]=PATTERN_GRID[v][u][0][0];
				this.pXYUV[v][u][1]=PATTERN_GRID[v][u][0][1];
			}
		}
		int numGood=0;
		int numBad=0;
		double [] distUV=new double[2];
		double dist;
		if (laserPointer!=null) {
			for (int i=0;i<uv.length;i++) if (uv[i]!=null) { //laserPointer == null > uv={}
				// Verify that laser spots are inside specified distance from the cell centers
				distUV[0]=reReMap[0][0]*uv[i][0]+reReMap[0][1]*uv[i][1]+reReMap[0][2]-laserPointer.laserUVMap[i][0];
				distUV[1]=reReMap[1][0]*uv[i][0]+reReMap[1][1]*uv[i][1]+reReMap[1][2]-laserPointer.laserUVMap[i][1];
				dist=Math.sqrt(distUV[0]*distUV[0]+distUV[1]*distUV[1]);
				if (debugLevel>1){
					System.out.println("Laser spot #"+i+", distance from predicted ="+ IJ.d2s(dist,3)+" ("+IJ.d2s(200*dist,3)+
							"% of cell radius), du="+IJ.d2s(distUV[0],3)+", dv="+IJ.d2s(distUV[1],3));
				}
				if ((2*dist)> laserPointer.maxOffsetFromCenter){
					String msg="Laser point "+(i+1)+"(of "+uv.length+") is too far from the specified location, and this check is enforced in the configuration\n"+
							"measured distance is "+ IJ.d2s(200*dist,1)+"% of the cell radius, specified is "+ IJ.d2s(100*laserPointer.maxOffsetFromCenter,1)+"%";
					System.out.println("Warning:"+msg);
					if (!noMessageBoxes) IJ.showMessage("Warning",msg);
					numBad++;
					uv[i]=null;
					continue;
				}
				numGood++;
			}
		}
		if ((debugLevel>0) && (numBad>0)){
			System.out.println("Removed "+numBad+" pointers that are too far from the predicted locations");
		}
		return numGood;
	}






	public int applyUVShiftRot(
			int [] UVShiftRot,
			double [][] uv,
			double [][] xyuv, // [][2], [][3] contain laser pointers u,v
			double maxOffsetFromCenter,
//			LaserPointer laserPointer,
			boolean noMessageBoxes
			){
		if (UVShiftRot!=null) this.UVShiftRot=UVShiftRot.clone();
		boolean has_lasers = false;
		for (double[] e:xyuv) if ((e != null) && (e.length>2)) {
			has_lasers = true;
			break;
		}
		int [][] reReMap=getRemapMatrix(UVShiftRot);
		if (debugLevel>1){
			//        		   System.out.println("rotation="+rotation+", reMap= [["+this.reMap[0][0]+","+this.reMap[0][1]+","+this.reMap[0][2]+"]["+
			//        				   +this.reMap[1][0]+","+this.reMap[1][1]+","+this.reMap[1][2]+"]]");
			System.out.println("reRot="+   UVShiftRot[2]+",  reReMap= [["+reReMap[0][0]+","+reReMap[0][1]+","+reReMap[0][2]+"]["+
					+reReMap[1][0]+","+reReMap[1][1]+","+reReMap[1][2]+"]]");
		}
		// calculate targetUV that maps PATTERN_GRID cells to the target (absolute) UV
		this.targetUV=new int    [this.PATTERN_GRID.length][this.PATTERN_GRID[0].length][];
		this.pXYUV=new double [this.PATTERN_GRID.length][this.PATTERN_GRID[0].length][];
		for (int v=0;v<this.PATTERN_GRID.length;v++) for (int u=0;u<this.PATTERN_GRID[v].length;u++){
			if ((this.PATTERN_GRID[v][u]==null) || (this.PATTERN_GRID[v][u][0]==null)) {
				this.targetUV[v][u]=null;
				this.pXYUV[v][u]=null;
			} else {
				this.targetUV[v][u]=new int [2];
				this.targetUV[v][u][0]=reReMap[0][0]*u+reReMap[0][1]*v+reReMap[0][2];
				this.targetUV[v][u][1]=reReMap[1][0]*u+reReMap[1][1]*v+reReMap[1][2];
				//        			   System.out.println("v="+v+", u="+u+", PATTERN_GRID.length="+PATTERN_GRID.length+", PATTERN_GRID[v].length="+PATTERN_GRID[v].length);
				//        			   System.out.println("this.pixelsUV.length="+this.pixelsUV.length);
				//        			   System.out.println("this.pixelsUV["+v+"].length="+this.pixelsUV[v].length);
				this.pXYUV[v][u]=new double [2];
				this.pXYUV[v][u][0]=PATTERN_GRID[v][u][0][0];
				this.pXYUV[v][u][1]=PATTERN_GRID[v][u][0][1];
			}
		}
		int numGood=0;
		int numBad=0;
		double [] distUV=new double[2];
		double dist;
		if (has_lasers) {
			for (int i=0;i<uv.length;i++) if (uv[i]!=null) { //laserPointer == null > uv={}
				// Verify that laser spots are inside specified distance from the cell centers
				distUV[0]=reReMap[0][0]*uv[i][0]+reReMap[0][1]*uv[i][1]+reReMap[0][2]-xyuv[i][2];
				distUV[1]=reReMap[1][0]*uv[i][0]+reReMap[1][1]*uv[i][1]+reReMap[1][2]-xyuv[i][3];
				dist=Math.sqrt(distUV[0]*distUV[0]+distUV[1]*distUV[1]);
				if (debugLevel>1){
					System.out.println("Laser spot #"+i+", distance from predicted ="+ IJ.d2s(dist,3)+" ("+IJ.d2s(200*dist,3)+
							"% of cell radius), du="+IJ.d2s(distUV[0],3)+", dv="+IJ.d2s(distUV[1],3));
				}
				if ((2*dist)> maxOffsetFromCenter){
					String msg="Laser point "+(i+1)+"(of "+uv.length+") is too far from the specified location, and this check is enforced in the configuration\n"+
							"measured distance is "+ IJ.d2s(200*dist,1)+"% of the cell radius, specified is "+ IJ.d2s(100*maxOffsetFromCenter,1)+"%";
					System.out.println("Warning:"+msg);
					if (!noMessageBoxes) IJ.showMessage("Warning",msg);
					numBad++;
					uv[i]=null;
					continue;
				}
				numGood++;
			}
		}
		if ((debugLevel>0) && (numBad>0)){
			System.out.println("Removed "+numBad+" pointers that are too far from the predicted locations");
		}
		return numGood;
	}


	/**
	 * Rotate/flip PATTERN_GRID to match expected
	 * @param hintGrid [v][u][0 - pixel X, 1 - pixel Y, 2 - targetU, 3 - targetV
	 * @return true if possible, false - if not
	 */
	/*
           public boolean applyHintToGrid(double [][][] hintGrid){

           }
	 */
	/**
	 * Calculate grid fractional UV from x,y and grid array
	 * using bi-linear interpolation from the nearest points.
	 * Iterates through all grid points, so it is not optimal
	 * for processing each pixel.
	 */
	public double [][] uvFromXY(
			double [][] xy,
			double maxDist
			){
		if (xy==null) {
			double[][] uv=new double[0][];
			return uv;
		}
		double[][] uv=new double[xy.length][];
		for (int i=0;i<xy.length;i++) {
			uv[i]= uvFromXY(xy[i],maxDist);
			//				   if (uv[i]==null) return null;
		}
		return uv;
	}
	public double [] uvFromXY(
			double [] xy,
			double maxDist
			){
		//			   double [][][][] grid= this.PATTERN_GRID;
		//			   double [] gXY=new double [2];
		//			   int [] iUV=new int [2];
		if (xy==null) return null;
		int [][] iUV=new int[3][2];
		int width=this.PATTERN_GRID[0].length;
		// find closest point to xy
		double dist2,dx,dy,minDist2=-1.0;
		double []dist2Array=new double [this.PATTERN_GRID.length*width];
		for (int v=0;v< this.PATTERN_GRID.length; v++) for (int u=0;u< this.PATTERN_GRID[v].length; u++)
			if ((this.PATTERN_GRID[v][u]!=null) && (this.PATTERN_GRID[v][u][0]!=null)){
				dx=this.PATTERN_GRID[v][u][0][0]-xy[0];
				dy=this.PATTERN_GRID[v][u][0][1]-xy[1];
				dist2=dx*dx+dy*dy;
				dist2Array[v*width+u]=dist2;
				if ((minDist2<0.0) || (minDist2>dist2)) {
					minDist2=dist2;
					iUV[0][0]=u;
					iUV[0][1]=v;
				}
			} else {
				dist2Array[v*width+u]=-1.0;
			}
		// now find two other closest points (not on the same line
		dist2Array[iUV[0][1]*width+iUV[0][0]]=-1.0; // mark used point
		int indx=0;
		minDist2=-1.0;
		for (int i=0; i<dist2Array.length;i++) if ((dist2Array[i]>=0.0) && ((minDist2<0.0) || (minDist2>dist2Array[i]))){
			minDist2=dist2Array[i];
			indx=i;
		}
		iUV[1][0]=indx%width;
		iUV[1][1]=indx/width;
		// mark all points on the same line as iUV[0] and iUV[1]
		// find closest of the remaining points
		indx=0;
		minDist2=-1.0;
		int dU1=iUV[1][0]-iUV[0][0];
		int dV1=iUV[1][1]-iUV[0][1];
		int dU2,dV2;
		for (int i=0; i<dist2Array.length;i++)
			if ((dist2Array[i]>=0.0) && ((minDist2<0.0) || (minDist2>dist2Array[i]))) {
				dU2=i%width-iUV[0][0];
				dV2=i/width-iUV[0][1];
				if (dU1*dV2!=dV1*dU2) {
					minDist2=dist2Array[i];
					indx=i;
				}
			}
		iUV[2][0]=indx%width;
		iUV[2][1]=indx/width;
		// now there are 3 (not co-linear) points to interpolate u,v
		double [][] aMuv={
				{iUV[1][0]-iUV[0][0],iUV[2][0]-iUV[0][0]},
				{iUV[1][1]-iUV[0][1],iUV[2][1]-iUV[0][1]}
		};
		Matrix Muv=new Matrix(aMuv);
		if (
				(this.PATTERN_GRID==null)||
				(this.PATTERN_GRID[iUV[0][1]][iUV[0][0]]==null)||
				(this.PATTERN_GRID[iUV[1][1]][iUV[1][0]]==null)||
				(this.PATTERN_GRID[iUV[2][1]][iUV[2][0]]==null)||
				(this.PATTERN_GRID[iUV[0][1]][iUV[0][0]][0]==null)||
				(this.PATTERN_GRID[iUV[1][1]][iUV[1][0]][0]==null)||
				(this.PATTERN_GRID[iUV[2][1]][iUV[2][0]][0]==null)) return null;
		double [][] aMxy={
				{   this.PATTERN_GRID[iUV[1][1]][iUV[1][0]][0][0]-this.PATTERN_GRID[iUV[0][1]][iUV[0][0]][0][0],
					this.PATTERN_GRID[iUV[2][1]][iUV[2][0]][0][0]-this.PATTERN_GRID[iUV[0][1]][iUV[0][0]][0][0]},
				{   this.PATTERN_GRID[iUV[1][1]][iUV[1][0]][0][1]-this.PATTERN_GRID[iUV[0][1]][iUV[0][0]][0][1],
						this.PATTERN_GRID[iUV[2][1]][iUV[2][0]][0][1]-this.PATTERN_GRID[iUV[0][1]][iUV[0][0]][0][1]}
		};
		Matrix Mxy=new Matrix(aMxy);
		double [][] aVxy={
				{xy[0]-this.PATTERN_GRID[iUV[0][1]][iUV[0][0]][0][0]},
				{xy[1]-this.PATTERN_GRID[iUV[0][1]][iUV[0][0]][0][1]}
		};
		Matrix Vxy=new Matrix(aVxy);
		double [][] aVuv0={
				{iUV[0][0]},
				{iUV[0][1]}
		};
		Matrix Vuv0=new Matrix(aVuv0);
		Matrix Vuv=Vuv0.plus(Muv.times(Mxy.inverse()).times(Vxy));
		double [] result=Vuv.getRowPackedCopy();
		if (this.debugLevel>1) System.out.println("X="+IJ.d2s(xy[0],3)+" Y="+IJ.d2s(xy[1],3));
		if (this.debugLevel>2) System.out.println(" "+
				"Grid["+iUV[0][1]+"]["+iUV[0][0]+"]X="+IJ.d2s(this.PATTERN_GRID[iUV[0][1]][iUV[0][0]][0][0],3)+" "+
				"Grid["+iUV[0][1]+"]["+iUV[0][0]+"]Y="+IJ.d2s(this.PATTERN_GRID[iUV[0][1]][iUV[0][0]][0][1],3)+"\n "+
				"Grid["+iUV[1][1]+"]["+iUV[1][0]+"]X="+IJ.d2s(this.PATTERN_GRID[iUV[1][1]][iUV[1][0]][0][0],3)+" "+
				"Grid["+iUV[1][1]+"]["+iUV[1][0]+"]Y="+IJ.d2s(this.PATTERN_GRID[iUV[1][1]][iUV[1][0]][0][1],3)+"\n "+
				"Grid["+iUV[2][1]+"]["+iUV[2][0]+"]X="+IJ.d2s(this.PATTERN_GRID[iUV[2][1]][iUV[2][0]][0][0],3)+" "+
				"Grid["+iUV[2][1]+"]["+iUV[2][0]+"]Y="+IJ.d2s(this.PATTERN_GRID[iUV[2][1]][iUV[2][0]][0][1],3));
		if (this.debugLevel>1) System.out.println("U="+IJ.d2s(result[0],3)+" V="+IJ.d2s(result[1],3)+"\n");
		minDist2=(result[0]-iUV[0][0])*(result[0]-iUV[0][0])+(result[1]-iUV[0][1])*(result[1]-iUV[0][1]);
		if ((maxDist>0.0) && (minDist2>maxDist*maxDist)) {
			if (this.debugLevel>0) System.out.println("minDist2="+minDist2+" (maxDist="+maxDist+") - pointer too far (x="+xy[0]+" y="+xy[1]+")");
			return null; // pointer too far from the grid (outside of the grid)
		}
		// change test (make sure that all 4 grid points around the result are defined
		int uFloor=(int) Math.floor(result[0]);
		int vFloor=(int) Math.floor(result[1]);
		int extra=(int) Math.round(maxDist)-1;
		if (extra<0) extra=0;
		for (int v=vFloor-extra; v<= vFloor+extra+1;v++)  for (int u=uFloor-extra; u<= uFloor+extra+1;u++)
			if ((v<0) || (u<0) || (v>=this.PATTERN_GRID.length) || (u>=this.PATTERN_GRID[v].length) ||
					(this.PATTERN_GRID[v][u]==null) || (this.PATTERN_GRID[v][u][0]==null)) {
				if (this.debugLevel>1) System.out.println("pointer="+result[0]+":"+result[1]+
						", no grid at "+u+":"+v+" - pointer does not have grid around (x="+xy[0]+" y="+xy[1]+"), extra="+extra+" vFloor="+vFloor+" uFloor="+uFloor);
				for (int iiv=-2;iiv<3;iiv++) {
					for (int iiu=-2;iiu<3;iiu++) {
						boolean iinValid=
								((iiv+vFloor)<0)||
								((iiv+vFloor)>=this.PATTERN_GRID.length) ||
								((iiu+uFloor)<0)||
								((iiu+uFloor)>=this.PATTERN_GRID[0].length) ||
								(this.PATTERN_GRID[iiv+vFloor][iiu+uFloor]==null) ||
								(this.PATTERN_GRID[iiv+vFloor][iiu+uFloor][0]==null);
						if (this.debugLevel>1)  System.out.println((iiu+uFloor)+":"+
								(iiv+vFloor)+ "  "+(iinValid?"---":(IJ.d2s(this.PATTERN_GRID[iiv+vFloor][iiu+uFloor][0][0],1)+":"+IJ.d2s(this.PATTERN_GRID[iiv+vFloor][iiu+uFloor][0][1],1))));
					}
				}
				return null; // pointer too far from the grid (outside of the grid)
			}
		return result;
	}

	/* ======================================================================== */
	/*
    public static MatchSimulatedPattern.LaserPointer LASER_POINTERS= new MatchSimulatedPattern.LaserPointer (

	 */
	//
	/* ======================================================================== */
	private  double [] correctedPatternCrossLocation(
			LwirReaderParameters lwirReaderParameters, // null is OK
			double [] beforeXY, // initial coordinates of the pattern cross point
			double wv0x,
			double wv0y,
			double wv1x,
			double wv1y,
			double [][] correction,
			ImagePlus imp,      // image data (Bayer mosaic)
			DistortionParameters distortionParameters, //
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			MatchSimulatedPattern matchSimulatedPattern, // correlationSize
			SimulationPattern.SimulParameters  thisSimulParameters,
			boolean equalizeGreens,
			double [] window,   // window function
			double [] window2,   // window function - twice FFT size (or null)
			double [] window4,   // window function - 4x FFT size (or null)
			SimulationPattern simulationPattern,
			boolean negative, // invert cross phase
			DoubleFHT fht_instance,
			boolean fast, // use fast measuring of the maximum on the correlation
			double [][] locsNeib, // locations and weights of neighbors to average
			int debug_level,
			String dbgStr){
		if (distortionParameters.legacyMode)
			return correctedPatternCrossLocationOld(
					beforeXY, // initial coordinates of the pattern cross point
					wv0x,
					wv0y,
					wv1x,
					wv1y,
					correction,
					imp,      // image data (Bayer mosaic)
					distortionParameters, //
					patternDetectParameters,
					matchSimulatedPattern, // correlationSize
					thisSimulParameters,
					equalizeGreens,
					window,   // window function
					window2,   // window function - twice FFT size (or null)
					window4,   // window function - 4x FFT size (or null)
					simulationPattern,
					negative, // invert cross phase
					fht_instance,
					fast, // use fast measuring of the maximum on the correlation
					locsNeib, // locations and weights of neighbors to average
					debug_level);
		else
			return correctedPatternCrossLocationAverage4(
					lwirReaderParameters, //  LwirReaderParameters lwirReaderParameters, // null is OK
					beforeXY, // initial coordinates of the pattern cross point
					wv0x,
					wv0y,
					wv1x,
					wv1y,
					correction,
					imp,      // image data (Bayer mosaic)
					distortionParameters, //
					patternDetectParameters,
					matchSimulatedPattern, // correlationSize
					thisSimulParameters,
					equalizeGreens,
					window,   // window function
					window2,   // window function - twice FFT size (or null)
					window4,   // window function - 4x FFT size (or null)
					simulationPattern,
					negative, // invert cross phase
					fht_instance,
					fast, // use fast measuring of the maximum on the correlation
					locsNeib, // locations and weights of neighbors to average
					debug_level,
					dbgStr);
	}
	private  double [] correctedPatternCrossLocationOld(
			double [] beforeXY, // initial coordinates of the pattern cross point
			double wv0x,
			double wv0y,
			double wv1x,
			double wv1y,
			double [][] correction,
			ImagePlus imp,      // image data (Bayer mosaic)
			DistortionParameters distortionParameters, //
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			MatchSimulatedPattern matchSimulatedPattern, // correlationSize
			SimulationPattern.SimulParameters  thisSimulParameters,
			boolean equalizeGreens,
			double [] window,   // window function
			double [] window2,   // window function - twice FFT size (or null)
			double [] window4,   // window function - 4x FFT size (or null)
			SimulationPattern simulationPattern,
			boolean negative, // invert cross phase
			DoubleFHT fht_instance,
			boolean fast, // use fast measuring of the maximum on the correlation
			double [][] locsNeib, // locations and weights of neighbors to average
			int debug_level){

		// Just for testing

		beforeXY[0]+=distortionParameters.correlationDx;  // offset, X (in pixels)
		beforeXY[1]+=distortionParameters.correlationDy; // offset y (in pixels)


		double [][] convMatrix= {{1.0,-1.0},{1.0,1.0}}; // from greens2 to pixel WV
		double [][] invConvMatrix= matrix2x2_scale(matrix2x2_invert(convMatrix),2.0);

		double [] result=new double [3];
		result[0]=beforeXY[0];
		result[1]=beforeXY[1];
		result[2]=0.0; // contrast


		if (fht_instance==null) fht_instance=new DoubleFHT(); // move upstream to reduce number of initializations


		//create diagonal green selection around ixc,iyc

		double [][]wv={{wv0x, wv0y},
				{wv1x, wv1y}};
		double [][] WVgreens=matrix2x2_mul(wv,invConvMatrix);
		if (debug_level>2) System.out.println("WVgreens[0][0]="+IJ.d2s(WVgreens[0][0],3)+
				" WVgreens[0][1]="+IJ.d2s(WVgreens[0][1],3)+
				" WVgreens[1][0]="+IJ.d2s(WVgreens[1][0],3)+
				" WVgreens[1][1]="+IJ.d2s(WVgreens[1][1],3));
		double [] dUV;
		double[][] sim_pix;
		double [] simGreensCentered;
		double [] modelCorr;
		double [] xyCorr={0.0,0.0};
		double [] centerXY;
		double 	contrast;
		int numNeib;
		double []corr=null;
		double [] neibCenter=new double[2];
		if (correction!=null) { // overwrite wave vectors
			wv[0][0]=correction[0][0];
			wv[0][1]=correction[0][1];
			wv[1][0]=correction[1][0];
			wv[1][1]=correction[1][1];
			if (correction[0].length>3) { // enough data for quadratic approximation
				corr=new double[10];
				corr[0]=correction[0][3]/4;
				corr[1]=correction[0][4]/4;
				corr[2]=correction[0][5]/4;
				corr[3]=correction[1][3]/4;
				corr[4]=correction[1][4]/4;
				corr[5]=correction[1][5]/4;
				corr[6]=0.0;
				corr[7]=0.0;
				corr[9]=0.0;
				corr[9]=0.0;
			}
		}
		double u_span=Math.sqrt(wv0x*wv0x+wv0y*wv0y)*distortionParameters.correlationSize;
		double v_span=Math.sqrt(wv1x*wv1x+wv1y*wv1y)*distortionParameters.correlationSize;
		double min_span=Math.min(u_span, v_span);
		int thisCorrelationSize=distortionParameters.correlationSize;
		double [] thisWindow=window;
		double uv_threshold=distortionParameters.minUVSpan*0.25*Math.sqrt(2.0);

		if (
				(min_span<uv_threshold) &&
				(window2!=null) &&
				(thisCorrelationSize<distortionParameters.maximalCorrelationSize)) { // trying to increase only twice
			thisCorrelationSize*=2;
			min_span*=2;
			thisWindow=window2;
			if (
					(min_span<uv_threshold) &&
					(window4!=null) &&
					(thisCorrelationSize<distortionParameters.maximalCorrelationSize)) {
				thisCorrelationSize*=2;
				min_span*=2;
				thisWindow=window4;
			}
		}
		setCorrelationSizesUsed(thisCorrelationSize);
		//		     if (thisCorrelationSize>distortionParameters.correlationSize) System.out.println("**** u/v span too small, increasing FFT size to "+thisCorrelationSize);
		if ((debug_level>0)&&(thisCorrelationSize>distortionParameters.correlationSize)) System.out.println("**** u/v span too small, increasing FFT size to "+thisCorrelationSize);
		Rectangle centerCross=correlationSelection(
				beforeXY, // initial coordinates of the pattern cross point
				//						distortionParameters.correlationSize);
				thisCorrelationSize);

		int ixc=centerCross.x+centerCross.width/2;
		int iyc=centerCross.y+centerCross.height/2;
		double [] diffBeforeXY={beforeXY[0]-ixc, beforeXY[1]-iyc};
		double[][] input_bayer=splitBayer (imp,centerCross,equalizeGreens);

		if (debug_level>3) SDFA_INSTANCE.showArrays(input_bayer,  true, "centered");
		if (debug_level>2) SDFA_INSTANCE.showArrays(input_bayer[4], "greens");
		if (debug_level>2) System.out.println("ixc="+ixc+" iyc="+iyc);
		double [] greens=normalizeAndWindow (input_bayer[4], thisWindow);

		if (debug_level>2) {
			System.out.println(" wv0x="+IJ.d2s(wv0x,5)+" wv0y="+IJ.d2s(wv0y,5));
			System.out.println(" wv1x="+IJ.d2s(wv1x,5)+" wv1y="+IJ.d2s(wv1y,5));
			System.out.println(" u-span="+IJ.d2s(u_span,3)+"  v-span="+IJ.d2s(v_span,3)+" threshold="+IJ.d2s(uv_threshold,3)+" ("+IJ.d2s(distortionParameters.minUVSpan,3)+")");
			if (corr!=null) {
				System.out.println(" Ax="+IJ.d2s(corr[0],8)+" Bx="+IJ.d2s(corr[1],8)+" Cx="+IJ.d2s(corr[2],8)+" Dx="+IJ.d2s(corr[6],8)+" Ex="+IJ.d2s(corr[7],8));
				System.out.println(" Ay="+IJ.d2s(corr[3],8)+" By="+IJ.d2s(corr[4],8)+" Cy="+IJ.d2s(corr[5],8)+" Dy="+IJ.d2s(corr[8],8)+" Ey="+IJ.d2s(corr[9],8));
			}
		}
		for (numNeib=0;numNeib<locsNeib.length;numNeib++) if (locsNeib[numNeib][2]!=0.0) {
			neibCenter[0]=diffBeforeXY[0]+locsNeib[numNeib][0];
			neibCenter[1]=diffBeforeXY[1]+locsNeib[numNeib][1];
			//			 dUV=matrix2x2_scale(matrix2x2_mul(wv,diffBeforeXY),-2*Math.PI);

			dUV=matrix2x2_scale(matrix2x2_mul(wv,neibCenter),-2*Math.PI);
			simulationPattern.simulatePatternFullPattern( // not thread safe
					wv0x,
					wv0y,
					dUV[0]+(negative?(-Math.PI/2):Math.PI/2), // negative?(-Math.PI/2):Math.PI/2,
					wv1x,
					wv1y,
					dUV[1]+Math.PI/2, //Math.PI/2,
					corr, //null, // no mesh distortion here
					thisSimulParameters.subdiv,// SIMUL.subdiv, - do not need high quality here
					thisCorrelationSize,
					true,// center for greens
					false);//boolean mono
			sim_pix= simulationPattern.extractSimulPatterns (
					thisSimulParameters,
					1,       // subdivide output pixels
					thisCorrelationSize,    // number of Bayer cells in width of the square selection (half number of pixels)
					0,
					0);
			if ((debug_level>2) && (numNeib==0)){
				//				 if (debug_level>2){
				System.out.println("==========Showing simul"+ixc+":"+iyc);
				SDFA_INSTANCE.showArrays(sim_pix[4].clone(),  "simul"+ixc+":"+iyc);
			}

			simGreensCentered= normalizeAndWindow (sim_pix[4], thisWindow);
			//				 if ((debug_level>2) && (numNeib==0)){
			if (debug_level>2){
				System.out.println("==========Showing simGreensCentered"+ixc+":"+iyc);

				SDFA_INSTANCE.showArrays(simGreensCentered.clone(),  "simGreensCentered"+ixc+":"+iyc);
				SDFA_INSTANCE.showArrays(greens.clone(),  "greensWidowed"+ixc+":"+iyc);
				//					 System.out.println("debug_level="+debug_level+" *** Remove next line ***");
				//					 sim_pix[14]=null; // make it crash here
			}
			modelCorr=fht_instance.correlate (greens.clone(),  // measured pixel array
					//						 modelCorr=fht_instance.correlate (greens,  // measured pixel array
					simGreensCentered,  // simulated (model) pixel array)
					//	                     distortionParameters.correlationHighPassSigma);
					distortionParameters.correlationHighPassSigma,
					distortionParameters.correlationLowPassSigma,
					distortionParameters.phaseCorrelationFraction);

			//				 if ((debug_level>2) && (numNeib==0)){
			if (debug_level>2){
				System.out.println("==========Showing modelCorr"+ixc+":"+iyc);
				SDFA_INSTANCE.showArrays(modelCorr, "modelCorr"+ixc+":"+iyc);
			}
			//				 xyCorr=new double[2]; //????????????????????
			// Use fast, but less precise method here ?
			//				 if (numNeib==0) System.out.println ("correctedPatternCrossLocation(): debugLevel="+debugLevel+" fast="+fast);
			if (fast) centerXY= correlationMaximum(
					modelCorr,
					distortionParameters.correlationMaxOffset,
					(debug_level>2) && (numNeib==0));
			else      centerXY= correlationMaximum(modelCorr,
					distortionParameters.correlationRadius,
					distortionParameters.correlationThreshold,	//double threshold, // fraction of maximum (slightly less than 1.0) to limit the top part of the maximum for centroid

					distortionParameters.correlationSubdiv,
					distortionParameters.correlationFFTSubdiv,
					fht_instance,
					distortionParameters.correlationMaxOffset,
					0.0, // lowpass filtering already done
					(debug_level>2) && ((numNeib==0) || (passNumber>1)));

			if (centerXY==null) {
				if (debug_level>0) System.out.println("Too far from the center0 ("+beforeXY[0]+"/"+beforeXY[1]+")");
				return null;
			}
			// Verify contrast (if specified) - only for the center sample (numNeib==0)
			if (numNeib==0) {
				double [] contrasts= correlationContrast(
						modelCorr,
						greens,
						WVgreens,    // wave vectors (same units as the pixels array)
						//							distortionParameters.correlationRingWidth,   // ring (around r=0.5 dist to opposite corr) width
						distortionParameters.contrastSelectSigmaCenter, // Gaussian sigma to select correlation centers (pixels, 2.0)
						distortionParameters.contrastSelectSigma, // Gaussian sigma to select correlation centers (fraction of UV period), 0.1
						//TODO: verify that displacement is correct here (sign, direction)
						centerXY[0],    //  x0,              // center coordinates
						centerXY[1],    //y0,
						"test-contrast");   // title base for optional plots names
				contrast=contrasts[0];

				result[2]=contrast;

				if (Double.isNaN(contrasts[0]) || ((distortionParameters.correlationMinContrast>0) && (contrasts[0]<distortionParameters.correlationMinContrast))) {
					if (debug_level>1) System.out.println("Center contrast too low - "+contrasts[0]+"<"+distortionParameters.correlationMinContrast);
					if (debug_level>1) System.out.println("Center contrast "+IJ.d2s(contrasts[0],3)+" ("+distortionParameters.correlationMinContrast+")"+
							" is too low ("+IJ.d2s(beforeXY[0],3)+"/"+IJ.d2s(beforeXY[1],3)+")->"+
							IJ.d2s(centerXY[0],3)+"/"+IJ.d2s(centerXY[1],3));
					return null;
				} else {
					if (debug_level>1) System.out.println("Contrast "+IJ.d2s(contrasts[0],3)+" ("+distortionParameters.correlationMinContrast+")"+
							" is good ("+IJ.d2s(beforeXY[0],3)+"/"+IJ.d2s(beforeXY[1],3)+")->"+
							IJ.d2s(centerXY[0],3)+"/"+IJ.d2s(centerXY[1],3));
				}


				if (Double.isNaN(contrasts[1]) || ((distortionParameters.correlationMinAbsoluteContrast>0) && (contrasts[1]<distortionParameters.correlationMinAbsoluteContrast))) {
					if (debug_level>1) System.out.println("Absolute contrast too low - "+contrasts[1]+"<"+distortionParameters.correlationMinAbsoluteContrast);
					if (debug_level>1) System.out.println("Absolute contrast "+IJ.d2s(contrasts[1],3)+" ("+distortionParameters.correlationMinAbsoluteContrast+")"+
							" is too low ("+IJ.d2s(beforeXY[0],3)+"/"+IJ.d2s(beforeXY[1],3)+")->"+
							IJ.d2s(centerXY[0],3)+"/"+IJ.d2s(centerXY[1],3));
					return null;
				}


				if (debug_level>2) System.out.println("Contarst="+contrast+" (legacy)");
			}
			if (debug_level>2) System.out.println("correctedPatternCrossLocation: Center x="+IJ.d2s(centerXY[0],3)+" y="+ 	IJ.d2s(centerXY[1],3));
			// convert from diagonal greens coordinates to sensor pixel coordinates
			xyCorr[0]+=(-centerXY[0]-centerXY[1])*locsNeib[numNeib][2];
			xyCorr[1]+=( centerXY[0]-centerXY[1])*locsNeib[numNeib][2];/*
				 if (debug_array!=null) {
					 debug_array[9][0]+=(-centerXY[0]-centerXY[1])*locsNeib[numNeib][2];
					 debug_array[9][1]+=( centerXY[0]-centerXY[1])*locsNeib[numNeib][2];
				 }
			 */
			if (debug_level>1) System.out.println("correctedPatternCrossLocation: dist="+IJ.d2s(Math.sqrt(xyCorr[0]*xyCorr[0]+xyCorr[1]*xyCorr[1]),4)+" xyCorr[0]="+IJ.d2s(xyCorr[0],4)+" xyCorr[1]="+ 	IJ.d2s(xyCorr[1],4));
		}
		// average 	xyCorr[]

		//			 result[0]=ixc-xyCorr[0];
		//			 result[1]=iyc-xyCorr[1];

		// disabling correction !!!!!!!!!!!!!!!!!!!!!!!
		result[0]=ixc-xyCorr[0]+diffBeforeXY[0];
		result[1]=iyc-xyCorr[1]+diffBeforeXY[1];
		//			 result[0]=ixc+diffBeforeXY[0];
		//			 result[1]=iyc+diffBeforeXY[1];

		if (debug_level>2) System.out.println("---correctedPatternCrossLocation: before x="+IJ.d2s(beforeXY[0],3)+" y="+IJ.d2s(beforeXY[1],3));
		if (debug_level>2) System.out.println("+++correctedPatternCrossLocation: after  x="+IJ.d2s(result[0],3)+" y="+IJ.d2s(result[1],3));
		//			 if (debug_level>0) System.out.println("---correctedPatternCrossLocation: before x="+IJ.d2s(beforeXY[0],3)+" y="+IJ.d2s(beforeXY[1],3));
		//			 if (debug_level>0) System.out.println("+++correctedPatternCrossLocation: after  x="+IJ.d2s(result[0],3)+" y="+IJ.d2s(result[1],3));
		return result;
	}



	private  double [] correctedPatternCrossLocationAverage4(
			LwirReaderParameters lwirReaderParameters, // null is OK
			double [] beforeXY, // initial coordinates of the pattern cross point
			double wv0x,
			double wv0y,
			double wv1x,
			double wv1y,
			double [][] correction,
			ImagePlus imp,      // image data (Bayer mosaic)
			DistortionParameters distortionParameters, //distortionParameters.refineCorrelations
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			MatchSimulatedPattern matchSimulatedPattern, // correlationSize
			SimulationPattern.SimulParameters  thisSimulParameters,
			boolean equalizeGreens,
			double [] window,   // window function
			double [] window2,   // window function - twice FFT size (or null)
			double [] window4,   // window function - 4x FFT size (or null)
			SimulationPattern simulationPattern,
			boolean negative, // invert cross phase
			DoubleFHT fht_instance,
			boolean fast, // use fast measuring of the maximum on the correlation
			double [][] locsNeib, // locations and weights of neighbors to average
			int debug_level,
			String dbgStr
			){
		if (imp == null) {
			return null;
		}
		boolean is_lwir = ((lwirReaderParameters != null) && lwirReaderParameters.is_LWIR(imp));
		int     correlation_size =     is_lwir ? distortionParameters.correlationSizeLwir :        distortionParameters.correlationSize;
		int     max_correlation_size = is_lwir ? distortionParameters.maximalCorrelationSizeLwir : distortionParameters.maximalCorrelationSize;
		boolean is_mono = false;
		try {
			is_mono = Boolean.parseBoolean((String) imp.getProperty("MONOCHROME"));
		} catch (Exception e) {
		}
		is_mono |= is_lwir;

		int debug_threshold = 3;
		// next print - same for good and bad, correction==null
		if (dbgStr!=null) System.out.println(dbgStr+ ": wv0x="+wv0x+" wv0y="+wv0y+ " wv1x="+wv1x+" wv1y="+wv1y+
				" beforeXY[0]="+beforeXY[0]+", beforeXY[1]="+beforeXY[1]+" correction is "+((correction==null)?"null":"not null"));


		boolean dbgThis=
				(Math.abs(beforeXY[0]-patternDetectParameters.debugX)<patternDetectParameters.debugRadius) &&
				(Math.abs(beforeXY[1]-patternDetectParameters.debugY)<patternDetectParameters.debugRadius);
		// dbgThis=true;
//		dbgThis=true;
		if (dbgThis) {
			System.out.println("correctedPatternCrossLocationAverage4(), beforeXY[0]="+beforeXY[0]+", beforeXY[1]="+beforeXY[1]);
			debug_level+=3;
		}
		// Just for testing
		beforeXY[0]+=distortionParameters.correlationDx;  // offset, X (in pixels)
		beforeXY[1]+=distortionParameters.correlationDy; // offset y (in pixels)

		double [][] invConvMatrix= {{1,0},{0,1}}; // identity
		if (!is_lwir) {
			double [][] convMatrix= {{1.0,-1.0},{1.0,1.0}}; // from greens2 to pixel WV
			invConvMatrix= matrix2x2_scale(matrix2x2_invert(convMatrix),2.0);
		}

		double [] result=new double [3];
		result[0]=beforeXY[0];
		result[1]=beforeXY[1];
		result[2]=0.0; // contrast


		if (fht_instance==null) fht_instance=new DoubleFHT(); // move upstream to reduce number of initializations


		//create diagonal green selection around ixc,iyc

		double [][]wv={{wv0x, wv0y},
				{wv1x, wv1y}};
		double [][] WVgreensMono=matrix2x2_mul(wv,invConvMatrix); // rotated for greens, same dir for mono (lwir)
		if (debug_level > debug_threshold) System.out.println("WVgreensMono[0][0]="+IJ.d2s(WVgreensMono[0][0],3)+
				" WVgreensMono[0][1]="+IJ.d2s(WVgreensMono[0][1],3)+
				" WVgreensMono[1][0]="+IJ.d2s(WVgreensMono[1][0],3)+
				" WVgreensMono[1][1]="+IJ.d2s(WVgreensMono[1][1],3));
		double [] dUV;
//		double [] simGreensCentered;
		double [] simCentered;
		//				 double [] modelCorr;

		double [] centerXY;
		double 	contrast;
		int numNeib;
		double []corr=null;
		double [] neibCenter=new double[2];
		if (correction!=null) { // overwrite wave vectors
			wv[0][0]=correction[0][0];
			wv[0][1]=correction[0][1];
			wv[1][0]=correction[1][0];
			wv[1][1]=correction[1][1];
			if (correction[0].length>3) { // enough data for quadratic approximation
				corr=new double[10];
				corr[0]=correction[0][3]/4;
				corr[1]=correction[0][4]/4;
				corr[2]=correction[0][5]/4;
				corr[3]=correction[1][3]/4;
				corr[4]=correction[1][4]/4;
				corr[5]=correction[1][5]/4;
				corr[6]=0.0;
				corr[7]=0.0;
				corr[9]=0.0;
				corr[9]=0.0;
			}
		}
		double u_span=Math.sqrt(wv0x*wv0x+wv0y*wv0y)*correlation_size;
		double v_span=Math.sqrt(wv1x*wv1x+wv1y*wv1y)*correlation_size;
		double min_span=Math.min(u_span, v_span);
		int thisCorrelationSize=correlation_size;
		double [] thisWindow=window;
		double uv_threshold=distortionParameters.minUVSpan*0.25*Math.sqrt(2.0);

		if (	(min_span < uv_threshold) &&
				(window2 != null) &&
				(thisCorrelationSize < max_correlation_size)) { // trying to increase only twice
			thisCorrelationSize *= 2;
			min_span *= 2;
			thisWindow=window2;
			if (	(min_span < uv_threshold) &&
					(window4 != null) &&
					(thisCorrelationSize < max_correlation_size)) {
				thisCorrelationSize *= 2;
				min_span *= 2;
				thisWindow = window4;
			}
		}

		setCorrelationSizesUsed(thisCorrelationSize);
		if ((debug_level > (debug_threshold - 2))&&(thisCorrelationSize>correlation_size)) System.out.println("**** u/v span too small, increasing FFT size to "+thisCorrelationSize);
		Rectangle centerCross=correlationSelection(
				beforeXY, // initial coordinates of the pattern cross point
				(is_lwir?(thisCorrelationSize/2):(thisCorrelationSize)));

		int ixc = centerCross.x + centerCross.width/2;
		int iyc = centerCross.y + centerCross.height/2;
		double [] diffBeforeXY={beforeXY[0]-ixc, beforeXY[1]-iyc};
		double [] greens_mono; // greens or mono
		if (is_lwir) {
			greens_mono=getNoBayer (imp,centerCross);
			if (debug_level > (debug_threshold +0)) SDFA_INSTANCE.showArrays(greens_mono, "greens_mono");
			if (debug_level > (debug_threshold +0)) System.out.println("ixc="+ixc+" iyc="+iyc);
			normalizeAndWindow (greens_mono, thisWindow);
		} else {
			double[][] input_bayer=splitBayer (imp,centerCross,equalizeGreens);
			if (debug_level > (debug_threshold +1)) SDFA_INSTANCE.showArrays(input_bayer,  true, "centered");
			if (debug_level > (debug_threshold +0)) SDFA_INSTANCE.showArrays(input_bayer[4], "greens");
			if (debug_level > (debug_threshold +0)) System.out.println("ixc="+ixc+" iyc="+iyc);
			greens_mono=normalizeAndWindow (input_bayer[4], thisWindow);
		}
		if (debug_level > (debug_threshold +0)) SDFA_INSTANCE.showArrays(greens_mono, "greens_mono_Windowed");
		// average is not zero - probably
		if (debug_level > (debug_threshold + 0)) {
			System.out.println(" wv0x="+IJ.d2s(wv0x,5)+" wv0y="+IJ.d2s(wv0y,5));
			System.out.println(" wv1x="+IJ.d2s(wv1x,5)+" wv1y="+IJ.d2s(wv1y,5));
			System.out.println(" u-span="+IJ.d2s(u_span,3)+"  v-span="+IJ.d2s(v_span,3)+" threshold="+IJ.d2s(uv_threshold,3)+" ("+IJ.d2s(distortionParameters.minUVSpan,3)+")");
			if (corr!=null) {
				System.out.println(" Ax="+IJ.d2s(corr[0],8)+" Bx="+IJ.d2s(corr[1],8)+" Cx="+IJ.d2s(corr[2],8)+" Dx="+IJ.d2s(corr[6],8)+" Ex="+IJ.d2s(corr[7],8));
				System.out.println(" Ay="+IJ.d2s(corr[3],8)+" By="+IJ.d2s(corr[4],8)+" Cy="+IJ.d2s(corr[5],8)+" Dy="+IJ.d2s(corr[8],8)+" Ey="+IJ.d2s(corr[9],8));
			}
		}
		int [][] gridNeib={{0,0},{0,1},{1,0},{1,1}};
		int numOfNeib=distortionParameters.correlationAverageOnRefine?gridNeib.length:1;
		if (debug_level > (debug_threshold + 0)) {
			System.out.println(" numOfNeib="+numOfNeib+" (distortionParameters.correlationAverageOnRefine="+distortionParameters.correlationAverageOnRefine);
		}
		if (locsNeib.length==1) {
			numOfNeib=1; // on the first pass, from legacy
			if (debug_level > (debug_threshold + 1)) {
				System.out.println("Reduced numOfNeib to "+numOfNeib+" as locsNeib.length="+locsNeib.length);
			}
		}
		if (dbgStr!=null) {
			double dbgSumWindow=0.0;
			for (double dbgD:thisWindow) dbgSumWindow+=dbgD;
			// All he same - good/bad
			System.out.println(dbgStr+ ": thisCorrelationSize="+thisCorrelationSize+" min_span="+min_span+ " dbgSumWindow="+dbgSumWindow+
					"locsNeib.length="+locsNeib.length+" fast="+fast+
					" numOfNeib="+numOfNeib+" (distortionParameters.correlationAverageOnRefine="+distortionParameters.correlationAverageOnRefine);
		}

		double [][] modelCorrs=     new double[numOfNeib][];
		double [][] debugGreens=new double[numOfNeib][0];
		for (numNeib=0;numNeib<numOfNeib;numNeib++) {
			if (is_lwir) { // monochrome, use all pixels
				neibCenter[0]=diffBeforeXY[0]+gridNeib[numNeib][0];
				neibCenter[1]=diffBeforeXY[1]+gridNeib[numNeib][1];
			} else {
				neibCenter[0]=diffBeforeXY[0]+0.5*(gridNeib[numNeib][0]+gridNeib[numNeib][1]);
				neibCenter[1]=diffBeforeXY[1]+0.5*(gridNeib[numNeib][0]-gridNeib[numNeib][1]);
			}
			double [] barray;
			if (is_lwir) {
				//negative=!negative;
				dUV=matrix2x2_scale(matrix2x2_mul(wv,neibCenter),-2*Math.PI);
//				dUV[0] = 0.0; dUV[1] = 0.0;
				if (debug_level > (debug_threshold + 20)){
					double []barray0 = simulationPattern.simulatePatternFullPatternSafe( // Is it the most time-consuming part? should it be done once and then only extraction separate?
							wv0x,
							wv0y,
							(negative?(-Math.PI/2):Math.PI/2), // negative?(-Math.PI/2):Math.PI/2,
							wv1x,
							wv1y,
							Math.PI/2, //Math.PI/2,
							corr, //null, // no mesh distortion here
							thisSimulParameters.subdiv,// SIMUL.subdiv, - do not need high quality here
							thisCorrelationSize,
							false,// false); // center for greens ???
							false);//boolean mono
					double []barray1 = simulationPattern.simulatePatternFullPatternSafe( // Is it the most time-consuming part? should it be done once and then only extraction separate?
							wv0x,
							wv0y,
							(negative?(-Math.PI/2):Math.PI/2), // negative?(-Math.PI/2):Math.PI/2,
							wv1x,
							wv1y,
							Math.PI/2, //Math.PI/2,
							corr, //null, // no mesh distortion here
							thisSimulParameters.subdiv,// SIMUL.subdiv, - do not need high quality here
							thisCorrelationSize,
							true,// false); // center for greens ???
							false);//boolean mono
					double []barray2 = simulationPattern.simulatePatternFullPatternSafe( // Is it the most time-consuming part? should it be done once and then only extraction separate?
							wv0x,
							wv0y,
							(negative?(-Math.PI/2):Math.PI/2), // negative?(-Math.PI/2):Math.PI/2,
							wv1x,
							wv1y,
							Math.PI/2, //Math.PI/2,
							corr, //null, // no mesh distortion here
							thisSimulParameters.subdiv,// SIMUL.subdiv, - do not need high quality here
							thisCorrelationSize,
							true,// false); // center for greens ???
							true);//boolean mono
					double [][] dbg_barray = {barray0, barray1, barray2};
					System.out.println(">=========Showing barray01"+ixc+":"+iyc);
					SDFA_INSTANCE.showArrays(dbg_barray, true, "barray"+ixc+":"+iyc);
					double [] sim_pix0= simulationPattern.extractSimulMono (
							barray0,
							thisSimulParameters,
							1,       // subdivide output pixels
							thisCorrelationSize,    // number of Bayer cells in width of the square selection (half number of pixels)
							0,
							0);
					double [] sim_pix1= simulationPattern.extractSimulMono (
							barray1,
							thisSimulParameters,
							1,       // subdivide output pixels
							thisCorrelationSize,    // number of Bayer cells in width of the square selection (half number of pixels)
							0,
							0);
					double [] sim_pix2= simulationPattern.extractSimulMono (
							barray2,
							thisSimulParameters,
							1,       // subdivide output pixels
							thisCorrelationSize,    // number of Bayer cells in width of the square selection (half number of pixels)
							0,
							0);
					double [][] dbg_sim_pix = {sim_pix0, sim_pix1, sim_pix2};
					System.out.println(">=========Showing barray01"+ixc+":"+iyc);
					SDFA_INSTANCE.showArrays(dbg_sim_pix, true, "sim_pix"+ixc+":"+iyc);

				}

				barray= simulationPattern.simulatePatternFullPatternSafe( // Is it the most time-consuming part? should it be done once and then only extraction separate?
						wv0x,
						wv0y,
						dUV[0]+(negative?(-Math.PI/2):Math.PI/2), // negative?(-Math.PI/2):Math.PI/2,
						wv1x,
						wv1y,
						dUV[1]+Math.PI/2, //Math.PI/2,
						corr, //null, // no mesh distortion here
						thisSimulParameters.subdiv,// SIMUL.subdiv, - do not need high quality here
						thisCorrelationSize,
						true,// false); // center for greens ???
						true);//boolean mono
				if (debug_level > (debug_threshold + 0)){
					System.out.println(">=========Showing barray"+ixc+":"+iyc);
					SDFA_INSTANCE.showArrays(barray, "barray"+ixc+":"+iyc);
				}
				// barray for dUV=={0,0} is symmetrical around center pixel,
				// sim_pix - around {center - 0.5, center - 0.5}
				// for center_for_g2 - sim_pix is symmetrical around [center,center],
				// for !mono && !center_for_g2 - [center+0.5, center+0.5]
				// TODO: reduce size of barray for mono twice in each direction
				double [] sim_pix= simulationPattern.extractSimulMono (
						barray,
						thisSimulParameters,
						1,       // subdivide output pixels
						thisCorrelationSize,    // number of Bayer cells in width of the square selection (half number of pixels)
						0,
						0);
				if (sim_pix==null){
					System.out.println("***** BUG: extractSimulPatterns() FAILED *****");
					return null;
				}
				if (dbgStr!=null) {
					double dbgSumWindow=0.0;
					for (double dbgD:sim_pix) dbgSumWindow+=dbgD;
					System.out.println(dbgStr+ ": SUM of sim_pix="+dbgSumWindow); // First difference good/bad
				}

				simCentered= normalizeAndWindow (sim_pix, thisWindow);

			} else {
				dUV=matrix2x2_scale(matrix2x2_mul(wv,neibCenter),-2*Math.PI);
				barray= simulationPattern.simulatePatternFullPatternSafe( // Is it the most time-consuming part? should it be done once and then only extraction separate?
						wv0x,
						wv0y,
						dUV[0]+(negative?(-Math.PI/2):Math.PI/2), // negative?(-Math.PI/2):Math.PI/2,
						wv1x,
						wv1y,
						dUV[1]+Math.PI/2, //Math.PI/2,
						corr, //null, // no mesh distortion here
						thisSimulParameters.subdiv,// SIMUL.subdiv, - do not need high quality here
						thisCorrelationSize,
						true, // center for greens
						false);//boolean mono
				if (debug_level > (debug_threshold + 0)){
					System.out.println(">=========Showing barray"+ixc+":"+iyc);
					SDFA_INSTANCE.showArrays(barray, "barray"+ixc+":"+iyc);
				}
//		double[][] sim_pix;
				double[][] sim_pix=null;
				if (is_mono) {
					sim_pix= new double[1][];
					sim_pix[0]=  simulationPattern.extractSimulMono ( // TODO: can use twice smaller barray
							barray,
							thisSimulParameters,
							1,  // subdivide output pixels - now 4
							thisCorrelationSize,    // number of Bayer cells in width of the square selection (half number of pixels)
							0,    // selection center, X (in pixels)
							0);
				} else {
					sim_pix= simulationPattern.extractSimulPatterns (
							barray,
							thisSimulParameters,
							1,       // subdivide output pixels
							thisCorrelationSize,    // number of Bayer cells in width of the square selection (half number of pixels)
							0,
							0);
				}
				if (sim_pix==null){
					System.out.println("***** BUG: extractSimulPatterns() FAILED *****");
					return null;
				}
				if (dbgStr!=null) {
					double dbgSumWindow=0.0;
					for (double[] dbgSlice:sim_pix) for (double dbgD:dbgSlice) dbgSumWindow+=dbgD;
					System.out.println(dbgStr+ ": SUM of sim_pix="+dbgSumWindow); // First difference good/bad


				}

				simCentered= normalizeAndWindow (sim_pix[4], thisWindow);
			}



			if (dbgStr!=null) {
				double dbgSumWindow=0.0;
				for (double dbgD:simCentered) dbgSumWindow+=dbgD;
				System.out.println(dbgStr+ ": SUM of simGreensCentered="+dbgSumWindow);
			}

			debugGreens[numNeib]=simCentered.clone();

			modelCorrs[numNeib]=fht_instance.phaseCorrelate (
					greens_mono.clone(),
					simCentered,
					patternDetectParameters.phaseCoeff,
					0,//   distortionParameters.correlationHighPassSigma,
					patternDetectParameters.lowpass_sigma, // (fast?distortionParameters.correlationLowPassSigma:0.0),// moved to decimation via FFT
					null,
					null);

			if (dbgStr!=null) {
				double dbgSumWindow=0.0;
				for (double[] dbgSlice:modelCorrs) for (double dbgD:dbgSlice) dbgSumWindow+=dbgD;
				System.out.println(dbgStr+ ": SUM of modelCorrs="+dbgSumWindow);
			}
		}
		if (debug_level > (debug_threshold + 0)){
			System.out.println(">=========Showing simCentered"+ixc+":"+iyc);
			SDFA_INSTANCE.showArrays(debugGreens, true, "simCentered"+ixc+":"+iyc);
		}
		if (debug_level > (debug_threshold + 0)){
			System.out.println(">=========Showing modelCorrs, passNumber="+passNumber);
			SDFA_INSTANCE.showArrays(modelCorrs, true, "modelCorrs:"+numOfNeib);
		}

		// combine 4 correlations into the double resolution, same output size (so half input size) array
		int halfSize=thisCorrelationSize/2;
		int qSize=thisCorrelationSize/4;
		int thisFFTSubdiv=distortionParameters.correlationFFTSubdiv;
		double thisLowpass=distortionParameters.correlationLowPassSigma;
		double [] modelCorr;
		if (numOfNeib>1) {
			modelCorr=new double [thisCorrelationSize*thisCorrelationSize];
			for (int i=0;i<modelCorr.length;    i++) modelCorr[i]=0.0;

			for (int dy=0;dy<2;dy++) for (int dx=0;dx<2;dx++)  {
				for (int y=0;y<halfSize;y++) for (int x=0;x<halfSize;x++) {
					modelCorr[(2*y+dy)*thisCorrelationSize+(2*x+dx)]+=
							modelCorrs[2*dy+dx][(qSize+y)*thisCorrelationSize+(qSize+x)];
				}
			}
			thisLowpass/=2.0; // the lower the value, the more filtering.  Decimated twice,so low pass filtering - accordingly
			thisFFTSubdiv=(thisFFTSubdiv>1)?(thisFFTSubdiv/2):1;
		} else {
			modelCorr=    modelCorrs[0];     // also - different size
		}

		if (debug_level > (debug_threshold + 0)){
			System.out.println(">==========Showing modelCorr");
			SDFA_INSTANCE.showArrays(modelCorr, thisCorrelationSize,thisCorrelationSize, "modelCorr");
		}

		if (fast) centerXY= correlationMaximum( // maybe twice actual size if
				modelCorr,
				distortionParameters.correlationMaxOffset,
				(debug_level > (debug_threshold + 0)) && (numNeib==0));  // low-pass filtering should already be done
		else      centerXY= correlationMaximum(
				modelCorr,
				distortionParameters.correlationRadius,
				distortionParameters.correlationThreshold,	//double threshold, // fraction of maximum (slightly less than 1.0) to limit the top part of the maximum for centroid

				distortionParameters.correlationSubdiv,
				thisFFTSubdiv,
				fht_instance,
				distortionParameters.correlationMaxOffset,
				thisLowpass, //distortionParameters.correlationLowPassSigma
				//						 (debug_level>2) && (passNumber>1));
				(debug_level > (debug_threshold + 0)));
		if (centerXY==null) {
			if (debug_level > (debug_threshold - 1)) System.out.println("Too far from the center01 ("+beforeXY[0]+"/"+beforeXY[1]+")");
			if (dbgStr!=null) System.out.println(dbgStr+ "- Too far from the center01 ("+beforeXY[0]+"/"+beforeXY[1]+")");
			return null;
		}

		if (numNeib>1){
			centerXY[0]*=0.5;
			centerXY[1]*=0.5;
			for (int i=0;i<2;i++) for (int j=0;j<2;j++) WVgreensMono[i][j]*=0.5;
		}



		double [] contrasts= correlationContrast(
				modelCorr,
				greens_mono,
				WVgreensMono,    // wave vectors (same units as the pixels array)
				distortionParameters.contrastSelectSigmaCenter, // Gaussian sigma to select correlation (pixels, 2.0)
				distortionParameters.contrastSelectSigma, // Gaussian sigma to select correlation centers (fraction of UV period), 0.1
				centerXY[0],    //  x0,              // center coordinates
				centerXY[1],    //y0,
				"test-contrast");   // title base for optional plots names
		if ((debug_level > (debug_threshold - 1))) {
			System.out.println("contrast = "+contrasts[0]);
		}
		contrast=contrasts[0];
		result[2]=contrast;
		if (Double.isNaN(contrasts[0]) || ((distortionParameters.correlationMinContrast>0) && (contrasts[0]<distortionParameters.correlationMinContrast))) {
			if ((debug_level > (debug_threshold - 1))) System.out.println("Contrast too low - "+contrasts[0]+"<"+distortionParameters.correlationMinContrast);
			if (debug_level > (debug_threshold - 1)) System.out.println("Contrast "+IJ.d2s(contrasts[0],3)+" ("+distortionParameters.correlationMinContrast+")"+
					" is TOO LOW ("+IJ.d2s(beforeXY[0],3)+"/"+IJ.d2s(beforeXY[1],3)+")->"+
					IJ.d2s(centerXY[0],3)+"/"+IJ.d2s(centerXY[1],3));
			if (dbgStr!=null) System.out.println(dbgStr+ " - Contrast "+IJ.d2s(contrasts[0],3)+" ("+distortionParameters.correlationMinContrast+")"+
					" is TOO LOW ("+IJ.d2s(beforeXY[0],3)+"/"+IJ.d2s(beforeXY[1],3)+")->"+
					IJ.d2s(centerXY[0],3)+"/"+IJ.d2s(centerXY[1],3));

			return null;
		} else {
			if (debug_level > (debug_threshold - 1)) System.out.println("Contrast "+IJ.d2s(contrasts[0],3)+" ("+distortionParameters.correlationMinContrast+")"+
					" is GOOD ("+IJ.d2s(beforeXY[0],3)+"/"+IJ.d2s(beforeXY[1],3)+")->"+
					IJ.d2s(centerXY[0],3)+"/"+IJ.d2s(centerXY[1],3));
			if (dbgStr!=null) System.out.println(dbgStr+ " - Contrast "+IJ.d2s(contrasts[0],3)+" ("+distortionParameters.correlationMinContrast+")"+
					" is GOOD ("+IJ.d2s(beforeXY[0],3)+"/"+IJ.d2s(beforeXY[1],3)+")->"+
					IJ.d2s(centerXY[0],3)+"/"+IJ.d2s(centerXY[1],3));
		}

		if (Double.isNaN(contrasts[1]) || ((distortionParameters.correlationMinAbsoluteContrast>0) && (contrasts[1]<distortionParameters.correlationMinAbsoluteContrast))) {
			if (debug_level > (debug_threshold - 1)) System.out.println("Absolute contrast too low - "+contrasts[1]+"<"+distortionParameters.correlationMinAbsoluteContrast);
			if (debug_level > (debug_threshold - 1)) System.out.println("Absolute contrast "+IJ.d2s(contrasts[1],3)+" ("+distortionParameters.correlationMinAbsoluteContrast+")"+
					" is too low ("+IJ.d2s(beforeXY[0],3)+"/"+IJ.d2s(beforeXY[1],3)+")->"+
					IJ.d2s(centerXY[0],3)+"/"+IJ.d2s(centerXY[1],3));
			if (dbgStr!=null) System.out.println(dbgStr+ " - Absolute contrast "+IJ.d2s(contrasts[1],3)+" ("+distortionParameters.correlationMinAbsoluteContrast+")"+
					" is too low ("+IJ.d2s(beforeXY[0],3)+"/"+IJ.d2s(beforeXY[1],3)+")->"+
					IJ.d2s(centerXY[0],3)+"/"+IJ.d2s(centerXY[1],3));
			return null;
		} else {
			if (dbgStr!=null) System.out.println(dbgStr+ " - Absolute contrast "+IJ.d2s(contrasts[1],3)+" ("+distortionParameters.correlationMinAbsoluteContrast+")"+
					" is GOOD ("+IJ.d2s(beforeXY[0],3)+"/"+IJ.d2s(beforeXY[1],3)+")->"+
					IJ.d2s(centerXY[0],3)+"/"+IJ.d2s(centerXY[1],3));
		}

		if (debug_level > (debug_threshold - 0))System.out.println(">>>Contrast="+contrasts[0]+"/"+contrasts[1]+" ("+IJ.d2s(beforeXY[0],3)+":"+IJ.d2s(beforeXY[1],3)+")->"+IJ.d2s(result[0],3)+":"+IJ.d2s(result[1],3));

// FIXME: maybe wrong for mono?
		if (is_lwir) {
			result[0]=ixc + diffBeforeXY[0] + centerXY[0];
			result[1]=iyc + diffBeforeXY[1] + centerXY[1];
		} else {
			result[0] = ixc + diffBeforeXY[0] - (-centerXY[0] -centerXY[1]);
			result[1] = iyc + diffBeforeXY[1] - ( centerXY[0] -centerXY[1]);
		}
		if (debug_level > (debug_threshold + 0)) System.out.println(">---correctedPatternCrossLocation: before x="+IJ.d2s(beforeXY[0],3)+" y="+IJ.d2s(beforeXY[1],3));
		if (debug_level > (debug_threshold + 0)) System.out.println(">+++correctedPatternCrossLocation: after  x="+IJ.d2s(result[0],3)+" y="+IJ.d2s(result[1],3));
		return result;
	}

	// ======= end of private  double [] correctedPatternCrossLocationAverage4() ===
	private  double [] correctedPatternCrossLocationAverage4TestOldNew(
			double [] beforeXY, // initial coordinates of the pattern cross point
			double wv0x,
			double wv0y,
			double wv1x,
			double wv1y,
			double [][] correction,
			ImagePlus imp,      // image data (Bayer mosaic)
			DistortionParameters distortionParameters, //distortionParameters.refineCorrelations
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			MatchSimulatedPattern matchSimulatedPattern, // correlationSize
			SimulationPattern.SimulParameters  thisSimulParameters,
			boolean equalizeGreens,
			double [] window,   // window function
			double [] window2,   // window function - twice FFT size (or null)
			double [] window4,   // window function - 4x FFT size (or null)
			SimulationPattern simulationPattern,
			boolean negative, // invert cross phase
			DoubleFHT fht_instance,
			boolean fast, // use fast measuring of the maximum on the correlation
			double [][] locsNeib, // locations and weights of neighbors to average
			int debug_level,
			String dbgStr
			){
		int debug_threshold = 3;
		// next print - same for good and bad, correction==null
		if (dbgStr!=null) System.out.println(dbgStr+ ": wv0x="+wv0x+" wv0y="+wv0y+ " wv1x="+wv1x+" wv1y="+wv1y+
				" beforeXY[0]="+beforeXY[0]+", beforeXY[1]="+beforeXY[1]+" correction is "+((correction==null)?"null":"not null"));


		boolean dbgThis=
				(Math.abs(beforeXY[0]-patternDetectParameters.debugX)<patternDetectParameters.debugRadius) &&
				(Math.abs(beforeXY[1]-patternDetectParameters.debugY)<patternDetectParameters.debugRadius);
		dbgThis=true;
		if (dbgThis) {
			System.out.println("correctedPatternCrossLocationAverage4(), beforeXY[0]="+beforeXY[0]+", beforeXY[1]="+beforeXY[1]);
			debug_level+=3;
		}
		//	System.out.println("correctedPatternCrossLocationAverage4(): beforeXY[0]="+beforeXY[0]+". beforeXY[1]="+beforeXY[1]);
		// Just for testing
		beforeXY[0]+=distortionParameters.correlationDx;  // offset, X (in pixels)
		beforeXY[1]+=distortionParameters.correlationDy; // offset y (in pixels)


		double [][] convMatrix= {{1.0,-1.0},{1.0,1.0}}; // from greens2 to pixel WV
		double [][] invConvMatrix= matrix2x2_scale(matrix2x2_invert(convMatrix),2.0);

		double [] result=new double [3];
		result[0]=beforeXY[0];
		result[1]=beforeXY[1];
		result[2]=0.0; // contrast


		if (fht_instance==null) fht_instance=new DoubleFHT(); // move upstream to reduce number of initializations


		//create diagonal green selection around ixc,iyc

		double [][]wv={{wv0x, wv0y},
				{wv1x, wv1y}};
		double [][] WVgreens=matrix2x2_mul(wv,invConvMatrix);
		if (debug_level > debug_threshold) System.out.println("WVgreens[0][0]="+IJ.d2s(WVgreens[0][0],3)+
				" WVgreens[0][1]="+IJ.d2s(WVgreens[0][1],3)+
				" WVgreens[1][0]="+IJ.d2s(WVgreens[1][0],3)+
				" WVgreens[1][1]="+IJ.d2s(WVgreens[1][1],3));
		double [] dUV;
		double[][] sim_pix;
		double [] simGreensCentered;
		//				 double [] modelCorr;

		double [] centerXY;
		double 	contrast;
		int numNeib;
		double []corr=null;
		double [] neibCenter=new double[2];
		if (correction!=null) { // overwrite wave vectors
			wv[0][0]=correction[0][0];
			wv[0][1]=correction[0][1];
			wv[1][0]=correction[1][0];
			wv[1][1]=correction[1][1];
			if (correction[0].length>3) { // enough data for quadratic approximation
				corr=new double[10];
				corr[0]=correction[0][3]/4;
				corr[1]=correction[0][4]/4;
				corr[2]=correction[0][5]/4;
				corr[3]=correction[1][3]/4;
				corr[4]=correction[1][4]/4;
				corr[5]=correction[1][5]/4;
				corr[6]=0.0;
				corr[7]=0.0;
				corr[9]=0.0;
				corr[9]=0.0;
			}
		}
		double u_span=Math.sqrt(wv0x*wv0x+wv0y*wv0y)*distortionParameters.correlationSize;
		double v_span=Math.sqrt(wv1x*wv1x+wv1y*wv1y)*distortionParameters.correlationSize;
		double min_span=Math.min(u_span, v_span);
		int thisCorrelationSize=distortionParameters.correlationSize;
		double [] thisWindow=window;
		double uv_threshold=distortionParameters.minUVSpan*0.25*Math.sqrt(2.0);

		if (
				(min_span<uv_threshold) &&
				(window2!=null) &&
				(thisCorrelationSize<distortionParameters.maximalCorrelationSize)) { // trying to increase only twice
			thisCorrelationSize*=2;
			min_span*=2;
			thisWindow=window2;
			if (
					(min_span<uv_threshold) &&
					(window4!=null) &&
					(thisCorrelationSize<distortionParameters.maximalCorrelationSize)) {
				thisCorrelationSize*=2;
				min_span*=2;
				thisWindow=window4;
			}
		}

		setCorrelationSizesUsed(thisCorrelationSize);
		if ((debug_level > (debug_threshold - 2))&&(thisCorrelationSize>distortionParameters.correlationSize)) System.out.println("**** u/v span too small, increasing FFT size to "+thisCorrelationSize);
		Rectangle centerCross=correlationSelection(
				beforeXY, // initial coordinates of the pattern cross point
				thisCorrelationSize);

		int ixc=centerCross.x+centerCross.width/2;
		int iyc=centerCross.y+centerCross.height/2;
		double [] diffBeforeXY={beforeXY[0]-ixc, beforeXY[1]-iyc};
		double[][] input_bayer=splitBayer (imp,centerCross,equalizeGreens);

		if (debug_level > (debug_threshold +1)) SDFA_INSTANCE.showArrays(input_bayer,  true, "centered");
		if (debug_level > (debug_threshold +0)) SDFA_INSTANCE.showArrays(input_bayer[4], "greens");
		if (debug_level > (debug_threshold +0)) System.out.println("ixc="+ixc+" iyc="+iyc);
		double [] greens=normalizeAndWindow (input_bayer[4], thisWindow);
		if (debug_level > (debug_threshold +0)) SDFA_INSTANCE.showArrays(greens, "greensWindowed");
		// average is not zero - probably

		if (debug_level > (debug_threshold + 0)) {
			System.out.println(" wv0x="+IJ.d2s(wv0x,5)+" wv0y="+IJ.d2s(wv0y,5));
			System.out.println(" wv1x="+IJ.d2s(wv1x,5)+" wv1y="+IJ.d2s(wv1y,5));
			System.out.println(" u-span="+IJ.d2s(u_span,3)+"  v-span="+IJ.d2s(v_span,3)+" threshold="+IJ.d2s(uv_threshold,3)+" ("+IJ.d2s(distortionParameters.minUVSpan,3)+")");
			if (corr!=null) {
				System.out.println(" Ax="+IJ.d2s(corr[0],8)+" Bx="+IJ.d2s(corr[1],8)+" Cx="+IJ.d2s(corr[2],8)+" Dx="+IJ.d2s(corr[6],8)+" Ex="+IJ.d2s(corr[7],8));
				System.out.println(" Ay="+IJ.d2s(corr[3],8)+" By="+IJ.d2s(corr[4],8)+" Cy="+IJ.d2s(corr[5],8)+" Dy="+IJ.d2s(corr[8],8)+" Ey="+IJ.d2s(corr[9],8));
			}
		}
		int [][] greenNeib={{0,0},{0,1},{1,0},{1,1}};
		int numOfNeib=distortionParameters.correlationAverageOnRefine?greenNeib.length:1;
		if (debug_level > (debug_threshold + 0)) {
			System.out.println(" numOfNeib="+numOfNeib+" (distortionParameters.correlationAverageOnRefine="+distortionParameters.correlationAverageOnRefine);
		}
		if (locsNeib.length==1) {
			numOfNeib=1; // on the first pass, from legacy
			if (debug_level > (debug_threshold + 0)) {
				System.out.println("Reduced numOfNeib to "+numOfNeib+" as locsNeib.length="+locsNeib.length);
			}
		}
		if (dbgStr!=null) {
			double dbgSumWindow=0.0;
			for (double dbgD:thisWindow) dbgSumWindow+=dbgD;
			// All he same - good/bad
			System.out.println(dbgStr+ ": thisCorrelationSize="+thisCorrelationSize+" min_span="+min_span+ " dbgSumWindow="+dbgSumWindow+
					"locsNeib.length="+locsNeib.length+" fast="+fast+
					" numOfNeib="+numOfNeib+" (distortionParameters.correlationAverageOnRefine="+distortionParameters.correlationAverageOnRefine);
		}

		double [][] modelCorrs=     new double[numOfNeib][];
		double [][] modelCorrs_new= new double[numOfNeib][];
		double [][] debugGreens=new double[numOfNeib][0];
		for (numNeib=0;numNeib<numOfNeib;numNeib++) {
			neibCenter[0]=diffBeforeXY[0]+0.5*(greenNeib[numNeib][0]+greenNeib[numNeib][1]);
			neibCenter[1]=diffBeforeXY[1]+0.5*(greenNeib[numNeib][0]-greenNeib[numNeib][1]);
			dUV=matrix2x2_scale(matrix2x2_mul(wv,neibCenter),-2*Math.PI);
			double [] barray= simulationPattern.simulatePatternFullPatternSafe( // Is it the most time-consuming part? should it be done once and then only extraction separate?
					wv0x,
					wv0y,
					dUV[0]+(negative?(-Math.PI/2):Math.PI/2), // negative?(-Math.PI/2):Math.PI/2,
					wv1x,
					wv1y,
					dUV[1]+Math.PI/2, //Math.PI/2,
					corr, //null, // no mesh distortion here
					thisSimulParameters.subdiv,// SIMUL.subdiv, - do not need high quality here
					thisCorrelationSize,
					true, // center for greens
					false);//boolean mono
			sim_pix= simulationPattern.extractSimulPatterns (
					barray,
					thisSimulParameters,
					1,       // subdivide output pixels
					thisCorrelationSize,    // number of Bayer cells in width of the square selection (half number of pixels)
					0,
					0);
			if (sim_pix==null){
				System.out.println("***** BUG: extractSimulPatterns() FAILED *****");
				return null;
			}
			if (dbgStr!=null) {
				double dbgSumWindow=0.0;
				for (double[] dbgSlice:sim_pix) for (double dbgD:dbgSlice) dbgSumWindow+=dbgD;
				System.out.println(dbgStr+ ": SUM of sim_pix="+dbgSumWindow); // First difference good/bad


			}

			simGreensCentered= normalizeAndWindow (sim_pix[4], thisWindow);

			if (dbgStr!=null) {
				double dbgSumWindow=0.0;
				for (double dbgD:simGreensCentered) dbgSumWindow+=dbgD;
				System.out.println(dbgStr+ ": SUM of simGreensCentered="+dbgSumWindow);
			}

			debugGreens[numNeib]=simGreensCentered.clone();

			// testing if phase reversal would exactly inverse result pattern - tested, perfect
			double [] simGreensCenteredClone = simGreensCentered.clone();

			modelCorrs[numNeib]=fht_instance.correlate (greens.clone(),  // measured pixel array
					//						 modelCorr=fht_instance.correlate (greens,  // measured pixel array
					simGreensCentered,  // simulated (model) pixel array)
					//	                     distortionParameters.correlationHighPassSigma);
					distortionParameters.correlationHighPassSigma,
					(fast ? distortionParameters.correlationLowPassSigma : 0.0),// moved to decimation via FFT
					distortionParameters.phaseCorrelationFraction);
			modelCorrs_new[numNeib]=fht_instance.phaseCorrelate (
					greens.clone(),
					simGreensCenteredClone,
					patternDetectParameters.phaseCoeff,
					0,//   distortionParameters.correlationHighPassSigma,
					patternDetectParameters.lowpass_sigma, // (fast?distortionParameters.correlationLowPassSigma:0.0),// moved to decimation via FFT
					null,
					null);

			if (dbgStr!=null) {
				double dbgSumWindow=0.0;
				for (double[] dbgSlice:modelCorrs) for (double dbgD:dbgSlice) dbgSumWindow+=dbgD;
				System.out.println(dbgStr+ ": SUM of modelCorrs="+dbgSumWindow);
			}


		}
		if (debug_level > (debug_threshold + 0)){
			System.out.println(">=========Showing simGreensCentered"+ixc+":"+iyc);
			SDFA_INSTANCE.showArrays(debugGreens, true, "simGreensCentered"+ixc+":"+iyc);
		}

		if (debug_level > (debug_threshold + 0)){
			System.out.println(">=========Showing modelCorrs, passNumber="+passNumber);
			SDFA_INSTANCE.showArrays(modelCorrs, true, "modelCorrs:"+numOfNeib);
			SDFA_INSTANCE.showArrays(modelCorrs_new, true, "modelCorrs_new:"+numOfNeib);
		}

		// combine 4 correlations into the double resolution, same output size (so half input size) array
		int halfSize=thisCorrelationSize/2;
		int qSize=thisCorrelationSize/4;
		int thisFFTSubdiv=distortionParameters.correlationFFTSubdiv;
		double thisLowpass=distortionParameters.correlationLowPassSigma;
		double [] modelCorr;
		double [] modelCorr_new;
		if (numOfNeib>1) {
			modelCorr=new double [thisCorrelationSize*thisCorrelationSize];
			modelCorr_new=new double [thisCorrelationSize*thisCorrelationSize];
			for (int i=0;i<modelCorr.length;    i++) modelCorr[i]=0.0;
			for (int i=0;i<modelCorr_new.length;i++) modelCorr_new[i]=0.0;

			for (int dy=0;dy<2;dy++) for (int dx=0;dx<2;dx++)  {
				for (int y=0;y<halfSize;y++) for (int x=0;x<halfSize;x++) {
					modelCorr[(2*y+dy)*thisCorrelationSize+(2*x+dx)]+=
							modelCorrs[2*dy+dx][(qSize+y)*thisCorrelationSize+(qSize+x)];
					modelCorr_new[(2*y+dy)*thisCorrelationSize+(2*x+dx)]+=
							modelCorrs_new[2*dy+dx][(qSize+y)*thisCorrelationSize+(qSize+x)];
				}
			}
			thisLowpass/=2.0; // the lower the value, the more filtering.  Decimated twice,so low pass filtering - accordingly
			thisFFTSubdiv=(thisFFTSubdiv>1)?(thisFFTSubdiv/2):1;
		} else {
			modelCorr=    modelCorrs[0];     // also - different size
			modelCorr_new=modelCorrs_new[0]; // also - different size
		}



		if (debug_level > (debug_threshold + 0)){
			System.out.println(">==========Showing modelCorr");
			SDFA_INSTANCE.showArrays(modelCorr, thisCorrelationSize,thisCorrelationSize, "modelCorr");
		}

		double [] centerXY_new;

		if (fast) centerXY_new= correlationMaximum( // maybe twice actual size if
				modelCorr_new,
				distortionParameters.correlationMaxOffset,
				(debug_level > (debug_threshold + 0)) && (numNeib==0));  // low-pass filtering should already be done
		else      centerXY_new= correlationMaximum(
				modelCorr_new,
				distortionParameters.correlationRadius,
				distortionParameters.correlationThreshold,	//double threshold, // fraction of maximum (slightly less than 1.0) to limit the top part of the maximum for centroid

				distortionParameters.correlationSubdiv,
				thisFFTSubdiv,
				fht_instance,
				distortionParameters.correlationMaxOffset,
				thisLowpass, //distortionParameters.correlationLowPassSigma
				//						 (debug_level>2) && (passNumber>1));
				(debug_level > (debug_threshold + 0)));

		if (fast) centerXY= correlationMaximum( // maybe twice actual size if
				modelCorr,
				distortionParameters.correlationMaxOffset,
				(debug_level > (debug_threshold + 0)) && (numNeib==0));  // low-pass filtering should already be done
		else      centerXY= correlationMaximum(
				modelCorr,
				distortionParameters.correlationRadius,
				distortionParameters.correlationThreshold,	//double threshold, // fraction of maximum (slightly less than 1.0) to limit the top part of the maximum for centroid

				distortionParameters.correlationSubdiv,
				thisFFTSubdiv,
				fht_instance,
				distortionParameters.correlationMaxOffset,
				thisLowpass, //distortionParameters.correlationLowPassSigma
				//						 (debug_level>2) && (passNumber>1));
				(debug_level > (debug_threshold + 0)));
		if (centerXY==null) {
			if (debug_level > (debug_threshold - 1)) System.out.println("Too far from the center01 ("+beforeXY[0]+"/"+beforeXY[1]+")");
			if (dbgStr!=null) System.out.println(dbgStr+ "- Too far from the center01 ("+beforeXY[0]+"/"+beforeXY[1]+")");
			return null;
		}

		//				 debug_level=3;

		if (numNeib>1){
			centerXY[0]*=0.5;
			centerXY[1]*=0.5;
			centerXY_new[0]*=0.5;
			centerXY_new[1]*=0.5;
			for (int i=0;i<2;i++) for (int j=0;j<2;j++) WVgreens[i][j]*=0.5;
		}

		double [] contrasts_new= correlationContrast(
				modelCorr_new,
				greens,
				WVgreens,    // wave vectors (same units as the pixels array)
				distortionParameters.contrastSelectSigmaCenter, // Gaussian sigma to select correlation centers (pixels, 2.0)
				distortionParameters.contrastSelectSigma, // Gaussian sigma to select correlation centers (fraction of UV period), 0.1
				centerXY[0],    //  x0,              // center coordinates
				centerXY[1],    //y0,
				"test-contrast-new");   // title base for optional plots names

		double [] contrasts= correlationContrast(
				modelCorr,
				greens,
				WVgreens,    // wave vectors (same units as the pixels array)
				distortionParameters.contrastSelectSigmaCenter, // Gaussian sigma to select correlation (pixels, 2.0)
				distortionParameters.contrastSelectSigma, // Gaussian sigma to select correlation centers (fraction of UV period), 0.1
				centerXY[0],    //  x0,              // center coordinates
				centerXY[1],    //y0,
				"test-contrast");   // title base for optional plots names
		if ((debug_level > (debug_threshold - 1))) {
			System.out.println("contrast_new = "+contrasts_new[0]+", contrast = "+contrasts[0]);
		}
		contrast=contrasts[0];
		result[2]=contrast;
		if (Double.isNaN(contrasts[0]) || ((distortionParameters.correlationMinContrast>0) && (contrasts[0]<distortionParameters.correlationMinContrast))) {
			if ((debug_level > (debug_threshold - 1))) System.out.println("Contrast too low - "+contrasts[0]+"<"+distortionParameters.correlationMinContrast);
			if (debug_level > (debug_threshold - 1)) System.out.println("Contrast "+IJ.d2s(contrasts[0],3)+" ("+distortionParameters.correlationMinContrast+")"+
					" is TOO LOW ("+IJ.d2s(beforeXY[0],3)+"/"+IJ.d2s(beforeXY[1],3)+")->"+
					IJ.d2s(centerXY[0],3)+"/"+IJ.d2s(centerXY[1],3));
			if (dbgStr!=null) System.out.println(dbgStr+ " - Contrast "+IJ.d2s(contrasts[0],3)+" ("+distortionParameters.correlationMinContrast+")"+
					" is TOO LOW ("+IJ.d2s(beforeXY[0],3)+"/"+IJ.d2s(beforeXY[1],3)+")->"+
					IJ.d2s(centerXY[0],3)+"/"+IJ.d2s(centerXY[1],3));

			return null;
		} else {
			if (debug_level > (debug_threshold - 1)) System.out.println("Contrast "+IJ.d2s(contrasts[0],3)+" ("+distortionParameters.correlationMinContrast+")"+
					" is GOOD ("+IJ.d2s(beforeXY[0],3)+"/"+IJ.d2s(beforeXY[1],3)+")->"+
					IJ.d2s(centerXY[0],3)+"/"+IJ.d2s(centerXY[1],3));
			if (dbgStr!=null) System.out.println(dbgStr+ " - Contrast "+IJ.d2s(contrasts[0],3)+" ("+distortionParameters.correlationMinContrast+")"+
					" is GOOD ("+IJ.d2s(beforeXY[0],3)+"/"+IJ.d2s(beforeXY[1],3)+")->"+
					IJ.d2s(centerXY[0],3)+"/"+IJ.d2s(centerXY[1],3));
		}


		if (Double.isNaN(contrasts[1]) || ((distortionParameters.correlationMinAbsoluteContrast>0) && (contrasts[1]<distortionParameters.correlationMinAbsoluteContrast))) {
			if (debug_level > (debug_threshold - 1)) System.out.println("Absolute contrast too low - "+contrasts[1]+"<"+distortionParameters.correlationMinAbsoluteContrast);
			if (debug_level > (debug_threshold - 1)) System.out.println("Absolute contrast "+IJ.d2s(contrasts[1],3)+" ("+distortionParameters.correlationMinAbsoluteContrast+")"+
					" is too low ("+IJ.d2s(beforeXY[0],3)+"/"+IJ.d2s(beforeXY[1],3)+")->"+
					IJ.d2s(centerXY[0],3)+"/"+IJ.d2s(centerXY[1],3));
			if (dbgStr!=null) System.out.println(dbgStr+ " - Absolute contrast "+IJ.d2s(contrasts[1],3)+" ("+distortionParameters.correlationMinAbsoluteContrast+")"+
					" is too low ("+IJ.d2s(beforeXY[0],3)+"/"+IJ.d2s(beforeXY[1],3)+")->"+
					IJ.d2s(centerXY[0],3)+"/"+IJ.d2s(centerXY[1],3));
			return null;
		} else {
			if (dbgStr!=null) System.out.println(dbgStr+ " - Absolute contrast "+IJ.d2s(contrasts[1],3)+" ("+distortionParameters.correlationMinAbsoluteContrast+")"+
					" is GOOD ("+IJ.d2s(beforeXY[0],3)+"/"+IJ.d2s(beforeXY[1],3)+")->"+
					IJ.d2s(centerXY[0],3)+"/"+IJ.d2s(centerXY[1],3));
		}

		if (debug_level > (debug_threshold - 0))System.out.println(">>>Contrast="+contrasts[0]+"/"+contrasts[1]+" ("+IJ.d2s(beforeXY[0],3)+":"+IJ.d2s(beforeXY[1],3)+")->"+IJ.d2s(result[0],3)+":"+IJ.d2s(result[1],3));
		result[0]=ixc-(-centerXY[0]-centerXY[1])+diffBeforeXY[0];
		result[1]=iyc-( centerXY[0]-centerXY[1])+diffBeforeXY[1];

		if (debug_level > (debug_threshold + 0)) System.out.println(">---correctedPatternCrossLocation: before x="+IJ.d2s(beforeXY[0],3)+" y="+IJ.d2s(beforeXY[1],3));
		if (debug_level > (debug_threshold + 0)) System.out.println(">+++correctedPatternCrossLocation: after  x="+IJ.d2s(result[0],3)+" y="+IJ.d2s(result[1],3));

		return result;
	}


	/* ======= Debugging only - returns 2-d array of x,y as a function of initial estimation =================== */
	public  double [][][] scanPatternCrossLocation(
			double range, // size of the scanning square
			int    size,  // number of scan points in each direction (total size*size)
			double [] beforeCenterXY, // initial coordinates of the pattern cross point
			double wv0x,
			double wv0y,
			double wv1x,
			double wv1y,
			ImagePlus imp,      // image data (Bayer mosaic)
			DistortionParameters distortionParameters, //
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			MatchSimulatedPattern matchSimulatedPattern, // correlationSize
			SimulationPattern.SimulParameters  thisSimulParameters,
			boolean equalizeGreens,
			double [] window,   // window function
			SimulationPattern simulationPattern,
			boolean negative, // invert cross phase
			DoubleFHT fht_instance
			){
		//		   	double [] result=beforeXY.clone();
		double [][][] result=new double [size][size][4];
		if (fht_instance==null) fht_instance=new DoubleFHT(); // move upstream to reduce number of initializations
		double [] beforeXY=new double[2];
		double [] filter=fht_instance.createFrequencyFilter(
				new double [distortionParameters.correlationSize*distortionParameters.correlationSize], //distortionParameters.correlationSize,
				distortionParameters.correlationHighPassSigma,
				distortionParameters.correlationLowPassSigma);
		if (debugLevel>2){
			double [] maskFull = new double [distortionParameters.correlationSize*distortionParameters.correlationSize];
			for (int i=0;i<maskFull.length;i++) {
				if (i<filter.length) maskFull[i]=filter[i];
				else {
					int rowMod = (distortionParameters.correlationSize - (i/distortionParameters.correlationSize)) % distortionParameters.correlationSize;
					int colMod = (distortionParameters.correlationSize - (i%distortionParameters.correlationSize)) % distortionParameters.correlationSize;
					maskFull[i]=filter[rowMod*distortionParameters.correlationSize+colMod];
				}

			}
			SDFA_INSTANCE.showArrays(maskFull, "filter");
		}


		for (int i=0;i<size;i++) for (int j=0;j<size;j++) {
			beforeXY[1]=beforeCenterXY[1]-range/2+ (range*i)/(size-1);
			beforeXY[0]=beforeCenterXY[0]-range/2+ (range*j)/(size-1);
			Rectangle centerCross=correlationSelection(
					beforeXY, // initial coordinates of the pattern cross point
					distortionParameters.correlationSize);

			int ixc=centerCross.x+centerCross.width/2;
			int iyc=centerCross.y+centerCross.height/2;
			double [] diffBeforeXY={beforeXY[0]-ixc, beforeXY[1]-iyc};

			//create diagonal green selection around ixc,iyc
			double[][] input_bayer=splitBayer (imp,centerCross,equalizeGreens);
			/*
			 double[][] corrWindow=null;
			 if (distortionParameters.correlationRadiusScale>=0.0) {
				 corrWindow=generateWeights (
						 distortionParameters.correlationWeightSigma,
						 distortionParameters.correlationRadiusScale); //  if 0 - use sigma as radius, inside - 1.0, outside 0.0. If >0 - size of array n*sigma

			 }

			 */
			if (debugLevel>3) SDFA_INSTANCE.showArrays(input_bayer,  true, "centered");
			if (debugLevel>1) System.out.println(i+"/"+j+": ixc="+ixc+" iyc="+iyc);
			// alternative way to  generate shifted pattern
			double [][]wv={{wv0x, wv0y},
					{wv1x, wv1y}};
			double [] dUV=matrix2x2_scale(matrix2x2_mul(wv,diffBeforeXY),-2*Math.PI);

			//correlationHighPassSigma
			double [] greens=normalizeAndWindow (input_bayer[4], window);
			simulationPattern.simulatePatternFullPattern(
					wv0x,
					wv0y,
					dUV[0]+(negative?(-Math.PI/2):Math.PI/2),
					wv1x,
					wv1y,
					dUV[1]+Math.PI/2, //0.0,
					null, // no mesh distortion here
					thisSimulParameters.subdiv,// SIMUL.subdiv, - do not need high quality here
					distortionParameters.correlationSize,
					true, // center for greens
					false);//boolean mono

			double[][] sim_pix= simulationPattern.extractSimulPatterns (
					thisSimulParameters,
					1,       // subdivide output pixels
					distortionParameters.correlationSize,    // number of Bayer cells in width of the square selection (half number of pixels)
					0.0, //-diffBeforeXY[0],
					0.0); //-diffBeforeXY[1]);

			double [] simGreensCentered= normalizeAndWindow (sim_pix[4], window);
			if (debugLevel>2) SDFA_INSTANCE.showArrays(greens.clone(), "greens-i"+i+"-j"+j);
			if (debugLevel>2) SDFA_INSTANCE.showArrays(simGreensCentered.clone(), "simGreensCentered-i"+i+"-j"+j);

			double [] modelCorr=fht_instance.correlate (greens,  // measured pixel array
					simGreensCentered,  // simulated (model) pixel array)
					//                     distortionParameters.correlationHighPassSigma);
					filter);
			if (debugLevel>2) SDFA_INSTANCE.showArrays(modelCorr.clone(), "modelCorr-i"+i+"-j"+j);
			double [] xyCorr=new double[2];
			double [] centerXY;
			//			 if (distortionParameters.correlationRadiusScale>=0.0)  centerXY= correlationMaximum(modelCorr,corrWindow);
			if (distortionParameters.correlationRadius>0){
				centerXY= correlationMaximum(modelCorr,
						distortionParameters.correlationRadius,
						distortionParameters.correlationThreshold,
						distortionParameters.correlationSubdiv,
						distortionParameters.correlationFFTSubdiv,
						fht_instance,
						distortionParameters.correlationMaxOffset,
						0.0, // low-pass filtering already done
						(debugLevel>2)
						);
			} 	 else centerXY= correlationMaximum(modelCorr,distortionParameters.correlationMaxOffset,(debugLevel>2));
			if (centerXY==null) {
				centerXY=new double[2];
				centerXY[0]=0.0;
				centerXY[1]=0.0;
			}
			if (debugLevel>2) System.out.println("correctedPatternCrossLocation: Center x="+IJ.d2s(centerXY[0],3)+" y="+ 	IJ.d2s(centerXY[1],3));
			xyCorr[0]=-centerXY[0]-centerXY[1];
			xyCorr[1]= centerXY[0]-centerXY[1];
			if (debugLevel>1) System.out.println("correctedPatternCrossLocation: "+i+"/"+j+": dist="+IJ.d2s(Math.sqrt(xyCorr[0]*xyCorr[0]+xyCorr[1]*xyCorr[1]),4)+" xyCorr[0]="+IJ.d2s(xyCorr[0],4)+" xyCorr[1]="+ 	IJ.d2s(xyCorr[1],4));
			//			 result[0]=ixc-xyCorr[0];
			//			 result[1]=iyc-xyCorr[1];
			result[i][j][0]=ixc-xyCorr[0]+diffBeforeXY[0];
			result[i][j][1]=iyc-xyCorr[1]+diffBeforeXY[1];
			result[i][j][2]=beforeXY[0];
			result[i][j][3]=beforeXY[1];
			if (debugLevel>1) System.out.println("---correctedPatternCrossLocation: "+i+"/"+j+" before x="+IJ.d2s(beforeXY[0],3)+" y="+IJ.d2s(beforeXY[1],3));
			if (debugLevel>1) System.out.println("+++correctedPatternCrossLocation: "+i+"/"+j+" after  x="+IJ.d2s(result[i][j][0],3)+" y="+IJ.d2s(result[i][j][1],3));
		}
		return result;
	}
	/*
			distortionParameters.correlationWeightSigma=  gd.getNextNumber();
			distortionParameters.correlationRadiusScale=  gd.getNextNumber();

	 */

	/* ======================================================================== */
	/**
	 * Interpolate maximum on a square correlation array, return vector from the center
	 */
	// one quater of the weights function to be used to approximate maximum on correlation by a second-degree polynominal
	public double [][] generateWeights (double sigma,
			double n) { //  if 0 - use sigma as radius, inside - 1.0, outside 0.0. If >0 - size of array n*sigma
		double r0=((n>0)?n:1.0)*sigma;
		double r2=r0*r0;
		int size = (int)  Math.ceil(r0);
		double [][] mask=new double [size][size];
		int i,j;
		double [] gaussian=new double [size];
		if (n>0) {
			double k=0.5/sigma/sigma;
			for (i=0;i<size;i++) gaussian[i]=Math.exp(-(k*i*i));
		}
		for (i=0;i<size;i++) for (j=0;j<size;j++){
			if ((i*i+j*j)>r2) mask[i][j]=0.0;
			else if (n>0)     mask[i][j]=gaussian[i]*gaussian[j];
			else              mask[i][j]=1.0;
		}
		return mask;
	}
	public double [] correlationMaximum(
			double [] corr,   // square (correlation) array to find location of the maximum ([0.0,0.0] in the center of the arrray)
			double [][] weights,
			double maxOffset) {
		if ((corr==null) || (corr.length==0)) return null;
		//  	   		double [] corrMax= new double[2];
		int size= (int) Math.sqrt(corr.length);
		int i,j,imax=0,ix,iy, ixc, iyc;
		double max=corr[0];
		for (i=1;i<corr.length;i++) if (max<corr[i]) {
			max=corr[i];
			imax=i;
		}
		iyc=imax/size;
		ixc=imax%size;
		int ixc0=ixc-size/2;
		int iyc0=iyc-size/2;
		if ((maxOffset>0) && (maxOffset*maxOffset<(ixc0*ixc0+iyc0*iyc0))) {
			if (debugLevel>1) System.out.println("Too far from the center1: ixc="+ixc+" iyc="+iyc);
			return null;
		}
		if (debugLevel>1) System.out.println("correlationMaximum: ixc="+ixc+" iyc="+iyc);


		/* ix, iy - the location of the point with maximal value. We'll approximate the vicinity of that maximum using a
		 * second degree polunominal:
   Z(x,y)~=A*x^2+B*y^2+C*x*y+D*x+E*y+F
   by minimizing sum of squared differences between the actual (Z(x,uy)) and approximated values.
   and then find the maximum on the approximated surface. Here is the math:

Z(x,y)~=A*x^2+B*y^2+C*x*y+D*x+E*y+F
minimizing squared error, using W(x,y) as weight function

error=Sum(W(x,y)*((A*x^2+B*y^2+C*x*y+D*x+E*y+F)-Z(x,y))^2)

error=Sum(W(x,y)*(A^2*x^4 + 2*A*x^2*(B*y^2+C*x*y+D*x+E*y+F-Z(x,y)) +(...) )
0=derror/dA=Sum(W(x,y)*(2*A*x^4 + 2*x^2*(B*y^2+C*x*y+D*x+E*y+F-Z(x,y)))
0=Sum(W(x,y)*(A*x^4 + x^2*(B*y^2+C*x*y+D*x+E*y+F-Z(x,y)))

SX4=Sum(W(x,y)*x^4), etc

(1) 0=A*SX4 + B*SX2Y2 + C*SX3Y +D*SX3 +E*SX2Y +F*SX2 - SZX2

derror/dB:

error=Sum(W(x,y)*(B^2*y^4 + 2*B*y^2*(A*x^2+C*x*y+D*x+E*y+F-Z(x,y)) +(...) )
0=derror/dB=Sum(W(x,y)*(2*B*y^4 + 2*y^2*(A*x^2+C*x*y+D*x+E*y+F-Z(x,y)))
0=Sum(W(x,y)*(B*y^4 + y^2*(A*x^2+C*x*y+D*x+E*y+F-Z(x,y)))

(2) 0=B*SY4 + A*SX2Y2 + C*SXY3 +D*SXY2 +E*SY3 +F*SY2 - SZY2
(2) 0=A*SX2Y2 + B*SY4 + C*SXY3 +D*SXY2 +E*SY3 +F*SY2 - SZY2

derror/dC:

error=Sum(W(x,y)*(C^2*x^2*y^2 + 2*C*x*y*(A*x^2+B*y^2+D*x+E*y+F-Z(x,y)) +(...) )
0=derror/dC=Sum(W(x,y)*(2*C*x^2*y^2 + 2*x*y*(A*x^2+B*y^2+D*x+E*y+F-Z(x,y)) )
0=Sum(W(x,y)*(C*x^2*y^2 + x*y*(A*x^2+B*y^2+D*x+E*y+F-Z(x,y)) )

(3) 0= A*SX3Y +  B*SXY3 +  C*SX2Y2 + D*SX2Y + E*SXY2 + F*SXY - SZXY

derror/dD:

error=Sum(W(x,y)*(D^2*x^2 + 2*D*x*(A*x^2+B*y^2+C*x*y+E*y+F-Z(x,y)) +(...) )
0=derror/dD=Sum(W(x,y)*(2*D*x^2 + 2*x*(A*x^2+B*y^2+C*x*y+E*y+F-Z(x,y)) )
0=Sum(W(x,y)*(D*x^2 + x*(A*x^2+B*y^2+C*x*y+E*y+F-Z(x,y)) )

(4) 0= A*SX3 +   B*SXY2 +  C*SX2Y + D*SX2  + E*SXY +  F*SX  - SZX

derror/dE:

error=Sum(W(x,y)*(E^2*y^2 + 2*E*y*(A*x^2+B*y^2+C*x*y+D*x+F-Z(x,y)) +(...) )
0=derror/dE=Sum(W(x,y)*(2*E*y^2 + 2*y*(A*x^2+B*y^2+C*x*y+D*x+F-Z(x,y)) )
0=Sum(W(x,y)*(E*y^2 + y*(A*x^2+B*y^2+C*x*y+D*x+F-Z(x,y)) )
(5) 0= A*SX2Y +  B*SY3 +   C*SXY2 + D*SXY +  E*SY2  + F*SY  - SZY

derror/dF:

error=Sum(W(x,y)*(F^2 +  2*F*(A*x^2+B*y^2+C*x*y+D*x+E*y-Z(x,y)) +(...) )
0=derror/dF=Sum(W(x,y)*(2*F +  2*(A*x^2+B*y^2+C*x*y+D*x+E*y-Z(x,y)) )
0=Sum(W(x,y)*(F +  (A*x^2+B*y^2+C*x*y+D*x+E*y-Z(x,y)) )
(6) 0= A*SX2 +   B*SY2 +   C*SXY +  D*SX +   E*SY   + F*S   - SZ




(1) 0= A*SX4 +   B*SX2Y2 + C*SX3Y +  D*SX3 +  E*SX2Y + F*SX2 - SZX2
(2) 0= A*SX2Y2 + B*SY4 +   C*SXY3 +  D*SXY2 + E*SY3  + F*SY2 - SZY2
(3) 0= A*SX3Y +  B*SXY3 +  C*SX2Y2 + D*SX2Y + E*SXY2 + F*SXY - SZXY
(4) 0= A*SX3 +   B*SXY2 +  C*SX2Y +  D*SX2  + E*SXY  + F*SX  - SZX
(5) 0= A*SX2Y +  B*SY3 +   C*SXY2 +  D*SXY +  E*SY2  + F*SY  - SZY
(6) 0= A*SX2 +   B*SY2 +   C*SXY +   D*SX +   E*SY   + F*S   - SZ


(1) 0= A*S40 + B*S22 + C*S31 + D*S30 + E*S21 + F*S20 - SZ20
(2) 0= A*S22 + B*S04 + C*S13 + D*S12 + E*S03 + F*S02 - SZ02
(3) 0= A*S31 + B*S13 + C*S22 + D*S21 + E*S12 + F*S11 - SZ11
(4) 0= A*S30 + B*S12 + C*S21 + D*S20 + E*S11 + F*S10 - SZ10
(5) 0= A*S21 + B*S03 + C*S12 + D*S11 + E*S02 + F*S01 - SZ01
(6) 0= A*S20 + B*S02 + C*S11 + D*S10 + E*S01 + F*S00 - SZ00


we beed x,y of maximum, so
d(A*x^2+B*y^2+C*x*y+D*x+E*y+F)/dx=0
d(A*x^2+B*y^2+C*x*y+D*x+E*y+F)/dy=0
d()/dx=2*A*x+C*y+D=0
d()/dy=C*x+2*B*y+E=0


  | S40 S22 S31 S30 S21 S20 |   | A |   | SZ20 |
  | S22 S04 S13 S12 S03 S02 |   | B |   | SZ02 |
  | S31 S13 S22 S21 S12 S11 |   | C |   | SZ11 |
  | S30 S12 S21 S20 S11 S10 | * | D | = | SZ10 |
  | S21 S03 S12 S11 S02 S01 |   | E |   | SZ01 |
  | S20 S02 S11 S10 S01 S00 |   | F |   | SZ00 |


  | 2*A    C | * | x | = | -D |
  |   C  2*B |   | Y |   | -E |

		 */

		double S00=0.0,
				S10=0.0,S01=0.0,
				S20=0.0,S11=0.0,S02=0.0,
				S30=0.0,S21=0.0,S12=0.0,S03=0.0,
				S40=0.0,S31=0.0,S22=0.0,S13=0.0,S04=0.0,
				SZ00=0.0,
				SZ10=0.0,SZ01=0.0,
				SZ20=0.0,SZ11=0.0,SZ02=0.0;
		int wsize=weights.length;
		double w,z,x,x2,x3,x4,y,y2,y3,y4,wz;
		for (i=iyc-wsize+1;i<iyc+wsize;i++) if ((i>0) && (i<size)) for (j=ixc-wsize+1;j<ixc+wsize;j++) if ((j>0) && (j<size)) {
			iy=i-iyc;
			ix=j-ixc;
			w=weights[(iy>=0)?iy:-iy][(ix>=0)?ix:-ix];
			if (w>0) {
				z=corr[i*size+j];
				wz=w*z;
				x=ix;
				x2=x*x;
				x3=x2*x;
				x4=x3*x;
				y=iy;
				y2=y*y;
				y3=y2*y;
				y4=y3*y;
				S00+=w;
				S10+=w*x;
				S01+=w*y;
				S20+=w*x2;
				S11+=w*x*y;
				S02+=w*y2;
				S30+=w*x3;
				S21+=w*x2*y;
				S12+=w*x*y2;
				S03+=w*y3;
				S40+=w*x4;
				S31+=w*x3*y;
				S22+=w*x2*y2;
				S13+=w*x*y3;
				S04+=w*y4;
				SZ00+=wz;
				SZ10+=wz*x;
				SZ01+=wz*y;
				SZ20+=wz*x2;
				SZ11+=wz*x*y;
				SZ02+=wz*y2;
			}
		}
		/*
  | S40 S22 S31 S30 S21 S20 |   | A |   | SZ20 |
  | S22 S04 S13 S12 S03 S02 |   | B |   | SZ02 |
  | S31 S13 S22 S21 S12 S11 |   | C |   | SZ11 |
  | S30 S12 S21 S20 S11 S10 | * | D | = | SZ10 |
  | S21 S03 S12 S11 S02 S01 |   | E |   | SZ01 |
  | S20 S02 S11 S10 S01 S00 |   | F |   | SZ00 |

		 */
		double [][] mAarray= {
				{S40,S22,S31,S30,S21,S20},
				{S22,S04,S13,S12,S03,S02},
				{S31,S13,S22,S21,S12,S11},
				{S30,S12,S21,S20,S11,S10},
				{S21,S03,S12,S11,S02,S01},
				{S20,S02,S11,S10,S01,S00}};

		double [] zAarray={SZ20,SZ02,SZ11,SZ10,SZ01,SZ00};
		Matrix M=new Matrix (mAarray);
		Matrix Z=new Matrix (zAarray,6);
		double [] ABCDEF= M.solve(Z).getRowPackedCopy();
		/*
  | 2*A    C | * | x | = | -D |
  |   C  2*B |   | Y |   | -E |
		 */
		double [][] mXYarray= {{2*ABCDEF[0],ABCDEF[2]},{ABCDEF[2],2*ABCDEF[1]}};
		double []   mDEarray={-ABCDEF[3],-ABCDEF[4]};

		Matrix mXY=new Matrix (mXYarray);
		Matrix mDE= new Matrix (mDEarray,2);
		double [] corrMax= mXY.solve(mDE).getRowPackedCopy();
		if (debugLevel>1) System.out.println("correlationMaximum: ixc="+ixc+" iyc="+iyc+" corrMax[0]="+corrMax[0]+" corrMax[1]="+corrMax[1]);
		corrMax[0]+=ixc-size/2;
		corrMax[1]+=iyc-size/2;

		if (debugLevel>2){
			double [] approx=new double [size*size];
			for (i=0;i<approx.length;i++) approx[i]=0.0;
			for (i=iyc-wsize+1;i<iyc+wsize;i++) if ((i>0) && (i<size)) for (j=ixc-wsize+1;j<ixc+wsize;j++) if ((j>0) && (j<size)) {
				iy=i-iyc;
				ix=j-ixc;
				x=ix;
				y=iy;

				//  	   				z=corr[i*size+j];
				approx[i*size+j]=ABCDEF[0]*x*x+
						ABCDEF[1]*y*y+
						ABCDEF[2]*x*y+
						ABCDEF[3]*x+
						ABCDEF[4]*y+
						ABCDEF[5];
			}
			double [][] both= new double[2][];
			both[0]=corr;
			both[1]=approx;
			//corr
			SDFA_INSTANCE.showArrays(both, true, "corr-approx"); // stack
		}
		//  	   		if (debugLevel>2) System.out.println("correlationMaximum: ix="+ix+" iy="+iy);
		//  	   		if (debugLevel>2) System.out.println("correlationMaximum: maxInHor[0] ="+maxInHor[0]+ " maxInHor[1]= "+maxInHor[1]+ " maxInHor[2]= "+maxInHor[2]);
		//  	   		if (debugLevel>2) System.out.println("correlationMaximum: maxInVert[0]="+maxInVert[0]+" maxInVert[1]="+maxInVert[1]+" maxInVert[2]="+maxInVert[2]);
		return corrMax;
	}


	private double [] correlationMaximum(
			double [] corr,
			int dist,       // maximal distance from the maximum to consider
			double threshold, // fraction of maximum (slightly less than 1.0) to limit the top part of the maximum for centroid
			int decimate, // interpolate to finer grid (both FFT and linear)
			int decimateFFT, // should be power of 2
			DoubleFHT fht_instance,
			double maxOffset,
			double lowpass, // relative to original corr size (will be scaled for decimation). Will only be applied if decimateFFT >1!
			boolean showDebug
			){
		if ((corr==null) || (corr.length==0)) return null;
		int size= (int) Math.sqrt(corr.length);
		int i,j,imax=0,ixc,iyc,index;
		//	   	   	if (showDebug) System.out.println("correlationMaximum(), decimateFFT="+decimateFFT);
		/**
		 * Reduces size of the correlation area (using center part) and simultaneously interpolating pixels, so the result is a
		 * scaled version of the center (total FFT suize remains the same)
		 */
		if (showDebug){
			System.out.println("correlationMaximum(): decimate="+decimate+" decimateFFT="+decimateFFT);
		}
		if (decimateFFT>1){
			if (fht_instance==null) fht_instance=new DoubleFHT();
			double scale=decimateFFT*decimateFFT;
			double [] corr1=new double [corr.length];
			/**
			 * As we are interested only in the center part of the image, we'll use flat-top
			 * window based on Hamming.
			 */

			int size1=size/decimateFFT;
			for (i=0;i<corr1.length;i++) corr1[i]=0.0;
			double borderAverage=0;
			index=size*(size+1)*(decimateFFT-1)/decimateFFT/2;
			int i0=index,i1=index+size1,i2=i1+size1*size,i3=i2-size1; // on right and bottom edge goes 1 pixel outside of the used area
			for (i=0;i<size1;i++)	{
				/*	   	   				if (showDebug){
	   	   					System.out.println(":: size="+size+" size1="+size1+ " i="+i+" i0="+i0+" i1="+i1+" i2="+i2+" i3="+i3+" scale="+scale);
	   	   				}
				 */
				borderAverage+=corr[i0++]+corr[i1+=size1]+corr[i2--]+corr[i3-=size1];

			}
			borderAverage/=4*size1;
			double [] preHammingMod=fht_instance.getHamming1d(size1/2);
			double [] hammingMod=new double [size1];
			for (i=0;i<size1/4;i++) hammingMod[i]=preHammingMod[i];
			for (i=1;i<size1/4;i++) hammingMod[size1-i]=preHammingMod[i];
			for (i=size1/4;i<=(size1-size1/4);i++) hammingMod[i]=1.0;

			if (showDebug) System.out.println("scale="+scale+ " borderAverage="+borderAverage);
			for (i=0;i<size1;i++) for (j=0;j<size1;j++) {
				corr1[(i*size+j)*decimateFFT] = scale*(corr[index+i*size+j]-borderAverage)*hammingMod[i]*hammingMod[j]+borderAverage;
			}
			if (showDebug) SDFA_INSTANCE.showArrays(corr1.clone(), "decimatedForFFT");
			fht_instance.swapQuadrants(corr1);
			if (!fht_instance.transform(corr1,false)) return null; // direct FHT
			if (showDebug) SDFA_INSTANCE.showArrays(corr1.clone(), "FFT");
			// zero out aliases
			for (i=0; i<=size1/2;i++) for (j=size1/2+1;j<size -(size1/2);j++) corr1[i*size+j]=0.0;
			for (i=size1/2+1;i<size -(size1/2);i++) for (j=0;j<size;j++) corr1[i*size+j]=0.0;
			for (i=size -(size1/2); i<size;i++) for (j=size1/2+1;j<size -(size1/2);j++) corr1[i*size+j]=0.0;
			// apply window for now - just
			/*
	    	    	if (showDebug) {
	    	    		System.out.println("Getting hamming1d ("+size1+")");
	    	    	}
			 */
			double [] hamming=fht_instance.getHamming1d(size1).clone();
			/*
	    	    	if (showDebug) {
	    	    		for (i=0;i<hamming.length;i++) System.out.println("hamming["+i+"]="+hamming[i]);
	    	    	}
	   	   			if (showDebug) SDFA_INSTANCE.showArrays(corr1.clone(), "NO_ALIAS");
	   	   			// Combine with low-pass Gaussian (if it is >0)
	   	   			if (lowpass>0){
	   	   				double [] gaussian1d=fht_instance.getGaussian1d(lowpass,size1); // no need to divide by /decimateFFT as we use size1, not size
	   	   				for (i=0;i<hammingMod.length;i++) hamming[i]*=gaussian1d[i];
		    	    	if (showDebug) {
		    	    		System.out.println("lowpass="+lowpass);
		    	    		for (i=0;i<gaussian1d.length;i++) System.out.println("gaussian1d["+i+"]="+gaussian1d[i]);
		    	    	}
	   	   			}

	    	    	if (showDebug) {
	    	    		for (i=0;i<hamming.length;i++) System.out.println("hamming["+i+"]="+hamming[i]);
	    	    	}
			 */
			int halfSize1=size1/2, shiftZero=size-halfSize1;
			for (i=0;i<=size1;i++) for (j=0;j<=size1;j++){
				int im=i%size1,jm=j%size1;
				corr1[((i+shiftZero)%size)*size+((j+shiftZero)%size)]*=hamming[im]*hamming[jm];
			}


			if (showDebug) SDFA_INSTANCE.showArrays(corr1.clone(), "FFT-masked");
			if (!fht_instance.transform(corr1,true)) return null; // inverse FHT
			fht_instance.swapQuadrants(corr1);
			if (showDebug) SDFA_INSTANCE.showArrays(corr1.clone(), "decimatedAfterFFT");
			dist*=decimateFFT;
			maxOffset*=decimateFFT;
			decimate/=decimateFFT;
			corr=corr1; // replace
		}


		double max=corr[0];
		for (i=1;i<corr.length;i++) if (max<corr[i]) {
			max=corr[i];
			imax=i;
		}
		iyc=imax/size;
		ixc=imax%size;
		int ixc0=ixc-size/2;
		int iyc0=iyc-size/2;
		if ((maxOffset>0) && (maxOffset*maxOffset<(ixc0*ixc0+iyc0*iyc0))) {
			if (showDebug || (debugLevel>1)) System.out.println("Too far from the center2: ixc="+ixc+" iyc="+iyc+" ixc0="+ixc0+" iyc0="+iyc0+" maxOffset="+maxOffset);
			//					if (showDebug || (debugLevel>0)) System.out.println("Too far from the center2: ixc="+ixc+" iyc="+iyc+" ixc0="+ixc0+" iyc0="+iyc0+" maxOffset="+maxOffset);
			return null;
		}
		if (showDebug || (debugLevel>1)) System.out.println("correlationMaximum: ixc="+ixc+" iyc="+iyc+" ixc0="+ixc0+" iyc0="+iyc0+" maxOffset="+maxOffset);

		// reduce dist if it hits borders
		if (dist>iyc) dist=iyc;
		if (dist>ixc) dist=ixc;
		if (dist>(size-iyc-1)) dist=(size-iyc-1);
		if (dist>(size-ixc-1)) dist=(size-ixc-1);
		int interpSize=2*dist*decimate+1;
		double [][] cell=new double[2][2];
		double [] row = new double [2];
		double [] ki= new double [2];
		double kj;
		double [] interpCorr=new double[interpSize*interpSize];
		int i1,j1;
		int i1Range,j1Range;
		for (i=0;i<2*dist;i++) for (j=0;j<2*dist;j++) {
			index=(iyc-dist+i)*size+(ixc-dist+j);
			cell[0][0]= corr[index];
			cell[0][1]= corr[index+1];
			cell[1][0]= corr[index+size];
			cell[1][1]= corr[index+size+1];
			i1Range=decimate+((i==(2*dist-1))?1:0);
			j1Range=decimate+((j==(2*dist-1))?1:0);
			ki[0]=(cell[1][0]-cell[0][0])/decimate;
			ki[1]=(cell[1][1]-cell[0][1])/decimate;

			for (i1=0;i1<i1Range;i1++){
				row[0]=cell[0][0]+ki[0]*i1;
				row[1]=cell[0][1]+ki[1]*i1;
				kj=(row[1]-row[0])/decimate;
				for (j1=0;j1<j1Range;j1++) {
					interpCorr[(i*decimate+i1)*interpSize+(j*decimate+j1)]=row[0]+kj*j1;
				}
			}
		}
		// Gaussian blur the after linear interpolation, use sigma = 0.75* decimate ?
		// now find the maximal value on the border - it will be a threshold for a wave from the center
		double interpolationBlurSigma=0.75* decimate;
		DoubleGaussianBlur gb=new DoubleGaussianBlur();
		gb.blurDouble(interpCorr, interpSize, interpSize, interpolationBlurSigma, interpolationBlurSigma, 0.01);

		double limit=interpCorr[0];
		for (i=0;i<interpSize;i++) {
			if (limit<interpCorr[i]) limit=interpCorr[i];
			if (limit<interpCorr[interpSize*interpSize-i-1])   limit=interpCorr[interpSize*interpSize-i-1];
			if (limit<interpCorr[interpSize*i])                limit=interpCorr[interpSize*i];
			if (limit<interpCorr[interpSize*i+(interpSize-1)]) limit=interpCorr[interpSize*i+(interpSize-1)];
		}
		// Now modify the limit if it is below threshold*max (sharp maximum)
		if (limit <threshold*max) limit =threshold*max;

		// run wave from the center, border pixels <=limit, so no need to verify array limits

		List <Integer> pixelList=new ArrayList<Integer>(100);
		Integer Index, newIndex;
		int []clusterMap=new int[interpSize*interpSize];
		for (i=0;i<clusterMap.length;i++) clusterMap[i]=0;
		int [] dirs={-1,-interpSize-1,-interpSize,-interpSize+1,1,interpSize+1,interpSize,interpSize-1};
		Index=dist*decimate*(interpSize+1); // center
		pixelList.clear();
		pixelList.add (Index);
		if (showDebug || (debugLevel>1)) System.out.println("correlationMaximum: pixelList.add ("+Index+ "), i="+(Index/interpSize)+ " j= "+(Index%interpSize));
		clusterMap[Index]=1;
		while (pixelList.size()>0) {
			Index=pixelList.remove(0);
			for (i=0;i<dirs.length;i++) {
				newIndex=Index+dirs[i];
				if ((clusterMap[newIndex]==0) && (interpCorr[newIndex]>limit)){
					pixelList.add (newIndex);
					clusterMap[newIndex]=1;
				}
			}
		}
		// Calculate centroid
		double s=0.0,sx=0.0, sy=0.0,x,y,d;
		if (showDebug || (debugLevel>1)) System.out.println("correlationMaximum: dist ="+dist+ " decimate= "+decimate+ " interpSize= "+interpSize);
		for (i=0;i<clusterMap.length;i++) if (clusterMap[i]>0){
			x=(i%interpSize-dist*decimate);
			y=(i / interpSize-dist*decimate);
			d=interpCorr[i]-limit;
			s+=d;
			sx+=x*d;
			sy+=y*d;
		}
		double [] corrXY={sx/s/decimate+ixc-size/2,sy/s/decimate+iyc-size/2};
		if (showDebug || (debugLevel>1)) System.out.println("correlationMaximum: s ="+s+ " sx= "+sx+ " sy= "+sy);
		if (showDebug || (debugLevel>1)) System.out.println("correlationMaximum: sx/s/decimate ="+(sx/s/decimate)+ " sy/s/decimate= "+(sy/s/decimate));
		if (showDebug || (debugLevel>1)) System.out.println("correlationMaximum: dx="+IJ.d2s(corrXY[0],3)+" dy="+IJ.d2s(corrXY[1],3));

		//			    if ((debugLevel>1) && (showDebug)) {
		if (showDebug) {
			double [] decimatedMasked=interpCorr.clone();
			for (i=0;i<decimatedMasked.length;i++) {
				if (clusterMap[i]==0) decimatedMasked[i]=limit;
			}
			double [][] both={interpCorr,decimatedMasked};
			SDFA_INSTANCE.showArrays(both, true, "centerCorr");
		}
		if (decimateFFT>1) {
			corrXY[0]/=decimateFFT;
			corrXY[1]/=decimateFFT;
		}
		return corrXY;
	}

	private double [] correlationMaximum(
			double [] corr,
			double maxOffset,
			boolean showDebug) {
		if ((corr==null) || (corr.length==0)) return null;
		double [] corrMax= new double[2];
		int size= (int) Math.sqrt(corr.length);
		int i,imax=0,ix,iy;
		double max=corr[0];
		for (i=1;i<corr.length;i++) if (max<corr[i]) {
			max=corr[i];
			imax=i;
		}
		iy=imax/size;
		ix=imax%size;

		corrMax[0]=ix-size/2;
		corrMax[1]=iy-size/2;
		if ((maxOffset>0) && (maxOffset*maxOffset<(corrMax[0]*corrMax[0]+corrMax[1]*corrMax[1]))) {
			if (debugLevel>1) System.out.println("Too far from the center3: corrMax[0]="+corrMax[0]+" corrMax[1]="+corrMax[1]);
			return null;
		}

		if ((ix==0) || (iy==0) || (ix==(size-1))  || (iy==(size-1))) return corrMax; // on the border - no interpolation;

		double [] maxInHor= new double[3]; // locations of interpolated maximums for each of 3 rows   (iy-1, iy, iy+1)
		double [] maxInVert=new double[3]; // locations of interpolated maximums for each of 3 columns(ix-1, ix, ix+1)
		if (debugLevel>2) System.out.println("correlationMaximum: ix="+ix+" iy="+iy);

		for (i=0;i<3;i++) {
			maxInHor[i]= -0.5+(corr[imax+size*(i-1)]-corr[imax+size*(i-1)-1])/
					(2*corr[imax+size*(i-1)]-corr[imax+size*(i-1)-1]-corr[imax+size*(i-1)+1]);
			maxInVert[i]=-0.5+(corr[imax+(i-1)]-corr[imax+     (i-1)-size])/
					(2*corr[imax+ (i-1)]-corr[imax+(i-1)-size]-corr[imax+(i-1)+size]);
		}
		if (debugLevel>2) System.out.println("correlationMaximum: maxInHor[0] ="+maxInHor[0]+ " maxInHor[1]= "+maxInHor[1]+ " maxInHor[2]= "+maxInHor[2]);
		if (debugLevel>2) System.out.println("correlationMaximum: maxInVert[0]="+maxInVert[0]+" maxInVert[1]="+maxInVert[1]+" maxInVert[2]="+maxInVert[2]);
		int maxInHorIndex=0;
		int maxInVertIndex=0;
		if ((maxInHor[0] <maxInHor[1] ) && (maxInHor[0] <maxInHor[2]) ) maxInHorIndex=1;
		if ((maxInVert[0]<maxInVert[1]) && (maxInVert[0]<maxInVert[2])) maxInVertIndex=1;
		if (debugLevel>2) System.out.println("correlationMaximum: maxInHorIndex="+maxInHorIndex+" maxInVertIndex="+maxInVertIndex);
		/*
		 * y= (y0+x0(y1-y0))/(1-(x1-x0)(y1-y0))
		 * x= (x0+y0(x1-x0))/(1-(x1-x0)(y1-y0))
		 * d= (1-(x1-x0)(y1-y0))
		 * y= (y0+x0(y1-y0))/d
		 * x= (x0+y0(x1-x0))/d
		 */
		double d=1-(maxInHor[maxInHorIndex+1]-maxInHor[maxInHorIndex])*(maxInVert[maxInVertIndex+1]-maxInVert[maxInVertIndex]);
		corrMax[0]=(maxInHor [maxInHorIndex]+ maxInVert[maxInVertIndex]*(maxInHor [maxInHorIndex +1]-maxInHor [maxInHorIndex ]))/d;
		corrMax[1]=(maxInVert[maxInVertIndex]+maxInHor [maxInHorIndex ]*(maxInVert[maxInVertIndex+1]-maxInVert[maxInVertIndex]))/d;
		if (debugLevel>2) System.out.println("correlationMaximum: corrMax[0]="+corrMax[0]+" corrMax[1]="+corrMax[1]);
		corrMax[0]+=ix-size/2;
		corrMax[1]+=iy-size/2;
		return corrMax;
	}


	/* ======================================================================== */
	public Rectangle correlationSelection(
			double [] beforeXY, // initial coordinates of the pattern cross point
			int size){
		int ixc=2*((int) Math.round(beforeXY[0]/2));
		int iyc=2*((int) Math.round(beforeXY[1]/2));
		Rectangle centerCross=new Rectangle(ixc-size,
				iyc-size,
				2*size,2*size);
		return centerCross;
	}
	/* ======================================================================== */
	// Estimate center xy and wave vectors from the neigbors
	// returns {{x,y},{wv1x,wv1y},{wv2x,wv2y}}


	public double [][] estimateCell(
			double [][][][] grid,
			int [] uv0,
			double [][] weights, // quadrant of sample weights
			boolean useContrast, // do not use cells with undefined contrast
			boolean forceLinear,  // use linear approximation (instead of quadratic)
			double thresholdLin,  // thershold ratio of matrix determinant to norm for linear approximation (det too low - fail)
			double thresholdQuad  // thershold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
			){
		int dist=weights.length-1;
		int size=dist*2+1;
		double [][][] samples0 = new double [size*size][3][];
		int index=0;
		int [] uv=new int[2];
		double w;
		int maxU=-dist-1,minU=dist+1,maxV=-dist-1,minV=dist+1,maxUpV=-2*dist-1,minUpV=2*dist+1,maxUmV=-2*dist-1,minUmV=2*dist+1;
		for (int iDv=-dist;iDv<=dist;iDv++) for (int iDu=-dist;iDu<=dist;iDu++) {
			uv[0]=uv0[0]+iDu;
			uv[1]=uv0[1]+iDv;
			if ((!useContrast && isCellDefined(grid,uv)) || isCellDefinedC(grid,uv)) {
				w=weights[(iDv>=0)?iDv:-iDv][(iDu>=0)?iDu:-iDu];
				if (w!=0.0){
					if (maxU<iDu) maxU=iDu;
					if (minU>iDu) minU=iDu;
					if (maxV<iDv) maxV=iDv;
					if (minV>iDv) minV=iDv;
					if (maxUpV<(iDu+iDv)) maxUpV= iDu+iDv;
					if (minUpV>(iDu+iDv)) minUpV= iDu+iDv;
					if (maxUmV<(iDu-iDv)) maxUmV= iDu-iDv;
					if (minUmV>(iDu-iDv)) minUmV= iDu-iDv;
					samples0[index][0]=new double[2];
					samples0[index][1]=new double[useContrast?3:2];
					samples0[index][2]=new double[1];
					samples0[index][2][0]=w;
					samples0[index][0][0]=iDu;
					samples0[index][0][1]=iDv;
					samples0[index][1][0]=grid[uv[1]][uv[0]][0][0];
					samples0[index][1][1]=grid[uv[1]][uv[0]][0][1];
					if (useContrast){
						samples0[index][1][2]=grid[uv[1]][uv[0]][0][2]; // contrast
					}
					index++;
				}
			}
		}
		if (debugLevel>3) System.out.println(" maxU-minU="+(maxU-minU)+" maxV-minV="+(maxV-minV));
		if (debugLevel>3) System.out.println(" maxUpV-minUpV="+(maxUpV-minUpV)+" maxUmV-minUmV="+(maxUmV-minUmV));

		int diameter=maxU-minU;
		if (diameter>(maxV-minV)) diameter= maxV-minV;
		diameter*=2;
		if (diameter>(maxUpV-minUpV)) diameter= (maxUpV-minUpV);
		if (diameter>(maxUmV-minUmV)) diameter= (maxUmV-minUmV);
		if (debugLevel>3) System.out.println(" diameter="+diameter+" number="+index);
		if (diameter<2) return null;

		double [][][] samples = new double [index][][];
		System.arraycopy(samples0, 0, samples, 0, index);
		double [][] estimatedCell= interpolateQuadraticWithWvAtZero(
				samples,   // see quadraticApproximation()
				forceLinear || (diameter<5),  // use linear approximation diameter <4 should be enough, 5 - just to be safe
				thresholdLin,  // thershold ratio of matrix determinant to norm for linear approximation (det too low - fail)
				thresholdQuad);  // thershold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
		if ((estimatedCell==null) || (estimatedCell[0]==null) || useContrast) return estimatedCell;
		double contrast=Double.NaN;
		if (isCellDefined(grid,uv0)) {
			double [] xycOld=grid[uv0[1]][uv0[0]][0];
			if (xycOld.length>2) contrast=xycOld[2];
		}
		double [] xyc={estimatedCell[0][0],estimatedCell[0][1],contrast};
		estimatedCell[0]=xyc;
		return estimatedCell;
	}

	public double [] interpolateQuadratic(
			double [] xy,         // coordinates for which the interpolation is needed
			double [][][] data,   // see quadraticApproximation()
			boolean forceLinear,  // use linear approximation
			double thresholdLin,  // thershold ratio of matrix determinant to norm for linear approximation (det too low - fail)
			double thresholdQuad){  // thershold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
		double [][] coeff = new PolynomialApproximation(this.debugLevel).quadraticApproximation(
				data,
				forceLinear,  // use linear approximation
				thresholdLin,  // thershold ratio of matrix determinant to norm for linear approximation (det too low - fail)
				thresholdQuad);
		if (coeff==null) return null;
		double [] result = new double [coeff.length];
		int offset=(coeff[0].length>3)?3:0;
		for (int i=0;i<coeff.length;i++) {
			result[i]=coeff[i][offset+0]*xy[0]+coeff[i][offset+1]*xy[1]+coeff[i][offset+2];
			if (offset>0) result[i]+=coeff[i][0]*xy[0]*xy[0]+coeff[i][1]*xy[1]*xy[1]+coeff[i][2]*xy[0]*xy[1];
		}
		return result;
	}

	// returns {{x,y},{wv1x,wv1y},{wv2x,wv2y}}
	public double [][] interpolateQuadraticWithWvAtZero(
			double [][][] data,   // see quadraticApproximation()
			boolean forceLinear,  // use linear approximation
			double thresholdLin,  // thershold ratio of matrix determinant to norm for linear approximation (det too low - fail)
			double thresholdQuad){  // thershold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
		double [][] coeff = new PolynomialApproximation(this.debugLevel).quadraticApproximation(
				data,
				forceLinear,  // use linear approximation
				thresholdLin,  // thershold ratio of matrix determinant to norm for linear approximation (det too low - fail)
				thresholdQuad);
		if (coeff==null) return null;
		if (coeff.length!=2) return null;
		double [][] result = new double [3][2];
		double [][] uv2xy=new double [2][2];
		int offset=(coeff[0].length>3)?3:0;
		for (int i=0;i<2;i++) {
			result[0][i]=coeff[i][offset+2];
			uv2xy[0][i]=coeff[0][offset+i];
			uv2xy[1][i]=coeff[1][offset+i];
		}
		double[][] wv=matrix2x2_invert(matrix2x2_scale(uv2xy,2.0));
		for (int i=0;i<2;i++) {
			result[1][i]=wv[0][i];
			result[2][i]=wv[1][i];
		}
		if ((debugLevel>2) && forceLinear) {
			System.out.println("*************** interpolateQuadraticWithWvAtZero() linear forced, data.length="+data.length);
		}
		if (debugLevel>3) {
			for (int i=0;i<data.length;i++) {
				System.out.println(i+": uv=["+IJ.d2s(data[i][0][0],3)+":"+IJ.d2s(data[i][0][1],3)+"]"+
						" xy=["+IJ.d2s(data[i][1][0],3)+":"+IJ.d2s(data[i][1][1],3)+"]"+
						" weight="+IJ.d2s(data[i][2][0],3));
			}
			String dbgStr="";
			dbgStr+=" [["+IJ.d2s(coeff[0][0],5)+"/"+IJ.d2s(coeff[0][1],5)+"/"+IJ.d2s(coeff[0][2],5);
			if (coeff[0].length>3) dbgStr+="/"+IJ.d2s(coeff[0][3],5)+"/"+IJ.d2s(coeff[0][4],5)+"/"+IJ.d2s(coeff[0][5],5)+"]]";
			dbgStr+=" [["+IJ.d2s(coeff[1][0],5)+"/"+IJ.d2s(coeff[1][1],5)+"/"+IJ.d2s(coeff[1][2],5);
			if (coeff[1].length>3) dbgStr+="/"+IJ.d2s(coeff[1][3],5)+"/"+IJ.d2s(coeff[1][4],5)+"/"+IJ.d2s(coeff[1][5],5)+"]]";
			System.out.println(dbgStr);
			for (int i=0;i<2;i++) {
				System.out.println(i+": uv2xy="+IJ.d2s(uv2xy[i][0],3)+":"+IJ.d2s(uv2xy[i][1],3));
			}
			for (int i=0;i<2;i++) {
				System.out.println(i+": wv="+IJ.d2s(wv[i][0],3)+":"+IJ.d2s(wv[i][1],3));
			}
		}
		return result;
	}
	// calculate simulation parameters for quadratic distortion of the pattern, compatible with SimulationPattern class
	// Returns 2 lines {{wv1x, wv1y, u0, Ax, Bx, Cx},{wv2x, wv2y, v0, Ay, By, Cy}}
	// if quadratic is not possible, only {{wv1x, wv1y, u0},{wv2x, wv2y}} will be returned
	// or just null if even linear is not possible
	// data array consists of lines of either 2 or 3 vectors:
	//  2-element vector x,y
	//  2 element vector u,v
	//  optional 1- element vector w (weight of the sample)
	public double [][] getSimulationParametersFromGrid(
			double [][][][] grid,
			int [] uv0,          // U,V of the center point (for which the simulation pattern should be built
			double [] xy0,          // x,y of the center point (or null to use grid)
			double [][] weights, // quadrant of sample weights
			boolean forceLinear,  // use linear approximation (instead of quadratic)
			double thresholdLin,  // thershold ratio of matrix determinant to norm for linear approximation (det too low - fail)
			double thresholdQuad  // thershold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
			){
		int dist=weights.length-1;
		int size=dist*2+1;
		double [][][] samples0 = new double [size*size][3][];
		int index=0;
		int [] uv=new int[2];
		double w;
		if (xy0==null) {
			if (isCellDefined(grid,uv0)) {
				xy0=new double [2];
				xy0[0]=grid[uv0[1]][uv0[0]][0][0];
				xy0[1]=grid[uv0[1]][uv0[0]][0][1];
			} else {
				return null; //xy of the center is not known
			}
		}
		int maxU=-dist-1,minU=dist+1,maxV=-dist-1,minV=dist+1,maxUpV=-2*dist-1,minUpV=2*dist+1,maxUmV=-2*dist-1,minUmV=2*dist+1;
		for (int iDv=-dist;iDv<=dist;iDv++) for (int iDu=-dist;iDu<=dist;iDu++) {
			uv[0]=uv0[0]+iDu;
			uv[1]=uv0[1]+iDv;
			if ((uv[0]>=0) && (uv[1]>=0) && (uv[1]<grid.length) && (uv[0]<grid[uv[1]].length) && (isCellDefined(grid,uv))) {
				w=weights[(iDv>=0)?iDv:-iDv][(iDu>=0)?iDu:-iDu];
				if (w!=0.0){
					if (maxU<iDu) maxU=iDu;
					if (minU>iDu) minU=iDu;
					if (maxV<iDv) maxV=iDv;
					if (minV>iDv) minV=iDv;
					if (maxUpV<(iDu+iDv)) maxUpV= iDu+iDv;
					if (minUpV>(iDu+iDv)) minUpV= iDu+iDv;
					if (maxUmV<(iDu-iDv)) maxUmV= iDu-iDv;
					if (minUmV>(iDu-iDv)) minUmV= iDu-iDv;
					samples0[index][0]=new double[2];
					samples0[index][1]=new double[2];
					samples0[index][2]=new double[1];
					samples0[index][2][0]=w;
					samples0[index][0][0]=grid[uv[1]][uv[0]][0][0]-xy0[0];
					samples0[index][0][1]=grid[uv[1]][uv[0]][0][1]-xy0[1];
					samples0[index][1][0]=uv[0];
					samples0[index][1][1]=uv[1];

					if (debugLevel>20) {
						System.out.println("iDu="+iDu+" iDv="+iDv+" "+
								" uv[0]="+IJ.d2s(uv[0],3)+" uv[1]="+IJ.d2s(uv[1],3)+" "+
								" samples0["+index+"][0][0]="+IJ.d2s(samples0[index][0][0],3)+" samples0["+index+"][0][1]="+IJ.d2s(samples0[index][0][1],3)+" "+
								" samples0["+index+"][1][0]="+IJ.d2s(samples0[index][1][0],3)+" samples0["+index+"][1][1]="+IJ.d2s(samples0[index][1][1],3));
					}
					index++;
				}
			}
		}
		int diameter=maxU-minU;
		if (diameter>(maxV-minV)) diameter= maxV-minV;
		diameter*=2;
		if (diameter>(maxUpV-minUpV)) diameter= (maxUpV-minUpV);
		if (diameter>(maxUmV-minUmV)) diameter= (maxUmV-minUmV);
		if (debugLevel>2) System.out.println(" diameter="+diameter+" number="+index);
		if (diameter<2) return null;
		double [][][] samples = new double [index][][];
		System.arraycopy(samples0, 0, samples, 0, index);
		double [][] simulParams= getSimulationParametersFromSamples(
				samples,   // see quadraticApproximation()
				forceLinear || (diameter<5),  // use linear approximation diameter <4 should be enough, 5 - just to be safe
				thresholdLin,  // thershold ratio of matrix determinant to norm for linear approximation (det too low - fail)
				thresholdQuad);  // thershold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
		return simulParams;
	}


	public double [][] getSimulationParametersFromSamples(
			double [][][] data,   // see quadraticApproximation()
			boolean forceLinear,  // use linear approximation
			double thresholdLin,  // thershold ratio of matrix determinant to norm for linear approximation (det too low - fail)
			double thresholdQuad){  // thershold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
		double [][] coeff = new PolynomialApproximation(this.debugLevel).quadraticApproximation(
				data,
				forceLinear,  // use linear approximation
				thresholdLin,  // thershold ratio of matrix determinant to norm for linear approximation (det too low - fail)
				thresholdQuad);

		if (debugLevel>2) {
			for (int i=0;i<data.length;i++) {
				System.out.println(i+": xy=["+IJ.d2s(data[i][0][0],3)+":"+IJ.d2s(data[i][0][1],3)+"]"+
						" uv=["+IJ.d2s(data[i][1][0],3)+":"+IJ.d2s(data[i][1][1],3)+"]"+
						" weight="+IJ.d2s(data[i][2][0],3));
			}
			String dbgStr="";
			dbgStr+=" ["+IJ.d2s(coeff[0][0],5)+"/"+IJ.d2s(coeff[0][1],5)+"/"+IJ.d2s(coeff[0][2],5);
			if (coeff[0].length>3) dbgStr+="/"+IJ.d2s(coeff[0][3],5)+"/"+IJ.d2s(coeff[0][4],5)+"/"+IJ.d2s(coeff[0][5],5)+"]";
			dbgStr+=" ["+IJ.d2s(coeff[1][0],5)+"/"+IJ.d2s(coeff[1][1],5)+"/"+IJ.d2s(coeff[1][2],5);
			if (coeff[1].length>3) dbgStr+="/"+IJ.d2s(coeff[1][3],5)+"/"+IJ.d2s(coeff[1][4],5)+"/"+IJ.d2s(coeff[1][5],5)+"]";
			System.out.println(dbgStr);
		}


		if (coeff==null) return null;
		if (coeff.length!=2) return null;
		boolean isQuad=coeff[0].length>3;
		int offset=isQuad?3:0;
		double [][] result = new double [2][isQuad?6:3];
		double [][] xy2uv = new double [2][2];
		for (int i=0;i<2;i++) {
			result[i][2]=coeff[i][offset+2]; // F
			xy2uv[i][0]=coeff[i][offset+0]; // D
			xy2uv[i][1]=coeff[i][offset+1]; // E
			result[i][0]=0.5*xy2uv[i][0]; // 0.5 because uv grid is 0.5 (pos/neg)
			result[i][1]=0.5*xy2uv[i][1]; //
		}
		if (isQuad) {
			double [][] uv2xy=matrix2x2_invert(xy2uv);
			//        	   double [][] ABCuv={{coeff[0][0],coeff[0][1],0.5*coeff[0][2]},
			//        			              {coeff[1][0],coeff[1][1],0.5*coeff[1][2]}};
			double [][] ABCuv={
					{4*coeff[0][0],4*coeff[0][1],2*coeff[0][2]}, // correction that uv grid is 0.5
					{4*coeff[1][0],4*coeff[1][1],2*coeff[1][2]}};
			double [][] ABCxy=new double [2][3];
			for (int i=0;i<2;i++) for (int j=0;j<3;j++) {
				ABCxy[i][j]=0.0;
				for (int k=0;k<2;k++) ABCxy[i][j]+=uv2xy[i][k]*ABCuv[k][j];
			}
			for (int i=0;i<2;i++) for (int j=0;j<3;j++)result[i][j+3]=ABCxy[i][j];
		}
		return result;
	}

	public double [][] findPatternFromGrid(
			int x0, // top-left pixel of the square WOI
			int y0,
			int size, // size of square (pixels)
			double[] halfWindow,
			boolean forceLinear,  // use linear approximation (instead of quadratic)
			double thresholdLin,  // thershold ratio of matrix determinant to norm for linear approximation (det too low - fail)
			double thresholdQuad  // thershold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
			){ // only half-window - half height by half width
		if (this.PATTERN_GRID==null) {
			String msg="PATTERN_GRID is needed, but undefined";
			IJ.showMessage("Error",msg);
			throw new IllegalArgumentException (msg);
		}
		if (this.PATTERN_GRID.length==0) return null;
		double x1=x0,y1=y0;
		double x2=x1+size;
		double y2=y1+size;
		double x,y;
		List <Integer> nodeList=new ArrayList<Integer>(1000);
		Integer Index;
		int len=this.PATTERN_GRID[0].length;
		for (int v=0;v<this.PATTERN_GRID.length;v++) for (int u=0;u<len;u++)
			if ((this.PATTERN_GRID[v][u]!=null) && (this.PATTERN_GRID[v][u][0]!=null)){
				x=this.PATTERN_GRID[v][u][0][0];
				y=this.PATTERN_GRID[v][u][0][1];
				if ((x>=x1) && (x<x2) && (y>=y1) && (y<y2)) {
					Index=v*len+u;
					nodeList.add(Index);
				}

			}
		double [][][] samples =new double [nodeList.size()][3][];
		// pattern parameters are referenced to the center of the square
		double xc=x0+size/2;
		double yc=y0+size/2;
		int halfSize=size/2;
		for (int i=0;i<samples.length;i++){
			int uv=nodeList.get(i);
			int v=uv/len;
			int u=uv%len;
			samples[i][0]=new double [2];
			samples[i][0][0]=this.PATTERN_GRID[v][u][0][0]-xc;
			samples[i][0][1]=this.PATTERN_GRID[v][u][0][1]-yc;
			samples[i][1]=new double [2];
			samples[i][1][0]=u;
			samples[i][1][1]=v;
			samples[i][2]=new double [1];
			int iy=((int) Math.round((this.PATTERN_GRID[v][u][0][1]-y1)/2));
			int ix=((int) Math.round((this.PATTERN_GRID[v][u][0][0]-x1)/2));
			samples[i][2][0]=((iy>=0) && (iy<halfSize) && (ix>=0) && (ix<halfSize))?halfWindow[iy*halfSize+ix]:0.0;
		}
		double [][] simulParams= getSimulationParametersFromSamples(
				samples,   // see quadraticApproximation()
				forceLinear || (halfSize<5),  // use linear approximation diameter <4 should be enough, 5 - just to be safe
				thresholdLin,  // thershold ratio of matrix determinant to norm for linear approximation (det too low - fail)
				thresholdQuad);  // thershold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
		return simulParams;
	}


	/**		Moved to PolyninomialApproximation class */

	/**
	 * Approximate function z(x,y) as a second degree polynomial
	 * f(x,y)=A*x^2+B*y^2+C*x*y+D*x+E*y+F
	 * data array consists of lines of either 2 or 3 vectors:
	 *  2-element vector x,y
	 *  variable length vector z (should be the same for all samples)
	 *  optional 1- element vector w (weight of the sample)
	 *
	 * returns arrrray of vectors or null
	 * each vector (one per each z component) is either 6-element-  (A,B,C,D,E,F) if quadratic is possible and enabled
	 * or 3-element - (D,E,F) if linear is possible and quadratic is not possible or disbled
	 * returns null if not enough data even for the linear approximation

	 */

	/* ======================================================================== */
	/*
	   public double [][] quadraticApproximation(
			   double [][][] data,
			   boolean forceLinear,  // use linear approximation
			   double thresholdLin,  // thershold ratio of matrix determinant to norm for linear approximation (det too low - fail)
			   double thresholdQuad  // thershold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
			   ){
/* ix, iy - the location of the point with maximal value. We'll approximate the vicinity of that maximum using a
	 * second degree polynominal:
   Z(x,y)~=A*x^2+B*y^2+C*x*y+D*x+E*y+F
   by minimizing sum of squared differenceS00between the actual (Z(x,uy)) and approximated values.
   and then find the maximum on the approximated surface. Here iS00the math:

Z(x,y)~=A*x^2+B*y^2+C*x*y+D*x+E*y+F
minimizing squared error, using W(x,y) aS00weight function

error=Sum(W(x,y)*((A*x^2+B*y^2+C*x*y+D*x+E*y+F)-Z(x,y))^2)

error=Sum(W(x,y)*(A^2*x^4 + 2*A*x^2*(B*y^2+C*x*y+D*x+E*y+F-Z(x,y)) +(...) )
0=derror/dA=Sum(W(x,y)*(2*A*x^4 + 2*x^2*(B*y^2+C*x*y+D*x+E*y+F-Z(x,y)))
0=Sum(W(x,y)*(A*x^4 + x^2*(B*y^2+C*x*y+D*x+E*y+F-Z(x,y)))

S40=Sum(W(x,y)*x^4), etc

(1) 0=A*S40 + B*S22 + C*S31 +D*S30 +E*S21 +F*S20 - SZ20

derror/dB:

error=Sum(W(x,y)*(B^2*y^4 + 2*B*y^2*(A*x^2+C*x*y+D*x+E*y+F-Z(x,y)) +(...) )
0=derror/dB=Sum(W(x,y)*(2*B*y^4 + 2*y^2*(A*x^2+C*x*y+D*x+E*y+F-Z(x,y)))
0=Sum(W(x,y)*(B*y^4 + y^2*(A*x^2+C*x*y+D*x+E*y+F-Z(x,y)))

(2) 0=B*S04 + A*S22 + C*S13 +D*S12 +E*S03 +F*SY2 - SZ02
(2) 0=A*S22 + B*S04 + C*S13 +D*S12 +E*S03 +F*SY2 - SZ02

derror/dC:

error=Sum(W(x,y)*(C^2*x^2*y^2 + 2*C*x*y*(A*x^2+B*y^2+D*x+E*y+F-Z(x,y)) +(...) )
0=derror/dC=Sum(W(x,y)*(2*C*x^2*y^2 + 2*x*y*(A*x^2+B*y^2+D*x+E*y+F-Z(x,y)) )
0=Sum(W(x,y)*(C*x^2*y^2 + x*y*(A*x^2+B*y^2+D*x+E*y+F-Z(x,y)) )

(3) 0= A*S31 +  B*S13 +  C*S22 + D*S21 + E*S12 + F*S11 - SZ11

derror/dD:

error=Sum(W(x,y)*(D^2*x^2 + 2*D*x*(A*x^2+B*y^2+C*x*y+E*y+F-Z(x,y)) +(...) )
0=derror/dD=Sum(W(x,y)*(2*D*x^2 + 2*x*(A*x^2+B*y^2+C*x*y+E*y+F-Z(x,y)) )
0=Sum(W(x,y)*(D*x^2 + x*(A*x^2+B*y^2+C*x*y+E*y+F-Z(x,y)) )

(4) 0= A*S30 +   B*S12 +  C*S21 + D*S20  + E*S11 +  F*S10  - SZ10

derror/dE:

error=Sum(W(x,y)*(E^2*y^2 + 2*E*y*(A*x^2+B*y^2+C*x*y+D*x+F-Z(x,y)) +(...) )
0=derror/dE=Sum(W(x,y)*(2*E*y^2 + 2*y*(A*x^2+B*y^2+C*x*y+D*x+F-Z(x,y)) )
0=Sum(W(x,y)*(E*y^2 + y*(A*x^2+B*y^2+C*x*y+D*x+F-Z(x,y)) )
(5) 0= A*S21 +  B*S03 +   C*S12 + D*S11 +  E*SY2  + F*SY  - SZ01

derror/dF:

error=Sum(W(x,y)*(F^2 +  2*F*(A*x^2+B*y^2+C*x*y+D*x+E*y-Z(x,y)) +(...) )
0=derror/dF=Sum(W(x,y)*(2*F +  2*(A*x^2+B*y^2+C*x*y+D*x+E*y-Z(x,y)) )
0=Sum(W(x,y)*(F +  (A*x^2+B*y^2+C*x*y+D*x+E*y-Z(x,y)) )
(6) 0= A*S20 +   B*SY2 +   C*S11 +  D*S10 +   E*SY   + F*S00  - SZ00


(1) 0= A*S40 + B*S22 + C*S31 + D*S30 + E*S21 + F*S20 - SZ20
(2) 0= A*S22 + B*S04 + C*S13 + D*S12 + E*S03 + F*S02 - SZ02
(3) 0= A*S31 + B*S13 + C*S22 + D*S21 + E*S12 + F*S11 - SZ11
(4) 0= A*S30 + B*S12 + C*S21 + D*S20 + E*S11 + F*S10 - SZ10
(5) 0= A*S21 + B*S03 + C*S12 + D*S11 + E*S02 + F*S01 - SZ01
(6) 0= A*S20 + B*S02 + C*S11 + D*S10 + E*S01 + F*S00 - SZ00
	 * /
		   int zDim=data[0][1].length;

		   double w,z,x,x2,x3,x4,y,y2,y3,y4,wz;
		   int i,j,n=0;
		   double S00=0.0,
		   S10=0.0,S01=0.0,
		   S20=0.0,S11=0.0,S02=0.0,
		   S30=0.0,S21=0.0,S12=0.0,S03=0.0,
		   S40=0.0,S31=0.0,S22=0.0,S13=0.0,S04=0.0;
		   double [] SZ00=new double [zDim];
		   double [] SZ01=new double [zDim];
		   double [] SZ10=new double [zDim];
		   double [] SZ11=new double [zDim];
		   double [] SZ02=new double [zDim];
		   double [] SZ20=new double [zDim];
		   for (i=0;i<zDim;i++) {
			   SZ00[i]=0.0;
			   SZ01[i]=0.0;
			   SZ10[i]=0.0;
			   SZ11[i]=0.0;
			   SZ02[i]=0.0;
			   SZ20[i]=0.0;
		   }
		   for (i=0;i<data.length;i++)  {
			   w=(data[i].length>2)? data[i][2][0]:1.0;
			   if (w>0) {
				   n++;
				   x=data[i][0][0];
				   y=data[i][0][1];
				   x2=x*x;
				   y2=y*y;
				   S00+=w;
				   S10+=w*x;
				   S01+=w*y;
				   S11+=w*x*y;
				   S20+=w*x2;
				   S02+=w*y2;
				   if (!forceLinear) {
					   x3=x2*x;
					   x4=x3*x;
					   y3=y2*y;
					   y4=y3*y;
					   S30+=w*x3;
					   S21+=w*x2*y;
					   S12+=w*x*y2;
					   S03+=w*y3;
					   S40+=w*x4;
					   S31+=w*x3*y;
					   S22+=w*x2*y2;
					   S13+=w*x*y3;
					   S04+=w*y4;
				   }
				   for (j=0;j<zDim;j++) {
					   z=data[i][1][j];
					   wz=w*z;
					   SZ00[j]+=wz;
					   SZ10[j]+=wz*x;
					   SZ01[j]+=wz*y;
					   if (!forceLinear) {
						   SZ20[j]+=wz*x2;
						   SZ11[j]+=wz*x*y;
						   SZ02[j]+=wz*y2;
					   }
				   }

			   }
		   }
		   //need to decide if there is enough data for linear and quadratic
		   double [][] mAarrayL= {
				   {S20,S11,S10},
				   {S11,S02,S01},
				   {S10,S01,S00}};
		   Matrix M=new Matrix (mAarrayL);
		   Matrix Z;
 	   	   if (debugLevel>3) System.out.println(">>> n="+n+" det_lin="+M.det()+" norm_lin="+normMatix(mAarrayL));
 	   	   double nmL=normMatix(mAarrayL);
		   if ((nmL==0.0) || (Math.abs(M.det())/nmL<thresholdLin)) return null; // not enough data even for the linear approximation
		   double []zAarrayL=new double [3];
		   double [][] ABCDEF=new double[zDim][];
//		   double [] zAarrayL={SZ10,SZ01,SZ00};
		   for (i=0;i<zDim;i++) {
			   zAarrayL[0]=SZ10[i];
			   zAarrayL[1]=SZ01[i];
			   zAarrayL[2]=SZ00[i];
		       Z=new Matrix (zAarrayL,3);
		       ABCDEF[i]= M.solve(Z).getRowPackedCopy();
		   }
		   if (forceLinear) return ABCDEF;
		   // quote try quadratic approximation
		   double [][] mAarrayQ= {
				   {S40,S22,S31,S30,S21,S20},
				   {S22,S04,S13,S12,S03,S02},
				   {S31,S13,S22,S21,S12,S11},
				   {S30,S12,S21,S20,S11,S10},
				   {S21,S03,S12,S11,S02,S01},
				   {S20,S02,S11,S10,S01,S00}};
		   M=new Matrix (mAarrayQ);
 	   	   if (debugLevel>3) System.out.println("    n="+n+" det_quad="+M.det()+" norm_quad="+normMatix(mAarrayQ)+" data.length="+data.length);
 	   	   double nmQ=normMatix(mAarrayQ);
		   if ((nmQ==0.0) || (Math.abs(M.det())/normMatix(mAarrayQ)<thresholdQuad)) {
			   System.out.println("Using linear approximation, M.det()="+M.det()+" normMatix(mAarrayQ)="+normMatix(mAarrayQ)); //did not happen
			   return ABCDEF; // not enough data for the quadratic approximation, return linear
		   }
//		   double [] zAarrayQ={SZ20,SZ02,SZ11,SZ10,SZ01,SZ00};
		   double [] zAarrayQ=new double [6];
		   for (i=0;i<zDim;i++) {
			   zAarrayQ[0]=SZ20[i];
			   zAarrayQ[1]=SZ02[i];
			   zAarrayQ[2]=SZ11[i];
			   zAarrayQ[3]=SZ10[i];
			   zAarrayQ[4]=SZ01[i];
			   zAarrayQ[5]=SZ00[i];
			   Z=new Matrix (zAarrayQ,6);
			   ABCDEF[i]= M.solve(Z).getRowPackedCopy();
		   }
		   return ABCDEF;
	   }
//	calcualte "volume" made of the matrix row-vectors, placed orthogonally
// to be compared to determinant
	public double normMatix(double [][] a) {
        double d,norm=1.0;
        for (int i=0;i<a.length;i++) {
        	d=0;
        	for (int j=0;j<a[i].length;j++) d+=a[i][j]*a[i][j];
        	norm*=Math.sqrt(d);
        }
		return norm;
	}
	 */
	/* ======================================================================== */
	public double[][][][] setPatternGridArray(int size) {
		return setPatternGridArray(size,size);
	}
	public double[][][][] setPatternGridArray(int width, int height) {
		int i,j;
		double[][][][] result= new double [height][width][][];
		for (i=0;i<height;i++) for (j=0;j<width;j++) result[i][j]=null;
		return result;
	}


	public static class PatternDetectParameters {
		public double gaussWidth; // <=0 - use Hamming window
		public double corrGamma;
		public double corrSigma;
		public int    diffSpectrCorr;
		public double shrinkClusters;
		public int    multiplesToTry;
		public double deviation;
		public int    deviationSteps;
		public double highpass;
		public double corrRingWidth;
		public double minCorrContrast;
		public double minGridPeriod;
		public double maxGridPeriod;
		public double minGridPeriodLwir;
		public double maxGridPeriodLwir;
		public double debugX;
		public double debugY;
		public double debugRadius;
		// added for large cell findPattern (based on phase correlation)
		public boolean use_large_cells = false; // new method based on phase correlation should work with large cells,
		// so only first negative correlation (1/2 period) fits in window
		public double phaseCoeff =      0.5; // "phasiness" of correlation
		public double lowpass_sigma =   .3;  // for phase correlation - frequency fraction of maximal
		public double min_frac =        0.03; // do not use higher order autocorrelation if min/max
		// is weaker than this fraction of the zero maximum
		public double  min_sin =  0.5;  // minimal sine for the angle between two pattern vectors
		public boolean no_crazy = true; // fail if quadratic approximation fails or returns outside of +/- 1.5


		public PatternDetectParameters(
				double gaussWidth,
				double corrGamma,
				double corrSigma,
				int diffSpectrCorr,
				double shrinkClusters,
				int multiplesToTry,
				double deviation,
				int deviationSteps,
				double highpass,
				double corrRingWidth,
				double minCorrContrast,
				double minGridPeriod,
				double maxGridPeriod,
				double minGridPeriodLwir,
				double maxGridPeriodLwir,
				double debugX,
				double debugY,
				double debugRadius,
				boolean use_large_cells,
				double phaseCoeff,
				double lowpass_sigma,
				double min_frac,
				double  min_sin,
				boolean no_crazy
				) {
			this.gaussWidth=gaussWidth;
			this.corrGamma = corrGamma;
			this.corrSigma = corrSigma;
			this.diffSpectrCorr = diffSpectrCorr;
			this.shrinkClusters = shrinkClusters;
			this.multiplesToTry = multiplesToTry;
			this.deviation = deviation;
			this.deviationSteps = deviationSteps;
			this.highpass = highpass;
			this.corrRingWidth = corrRingWidth;
			this.minCorrContrast = minCorrContrast;
			this.minGridPeriod=minGridPeriod;
			this.maxGridPeriod=maxGridPeriod;
			this.minGridPeriodLwir=minGridPeriodLwir;
			this.maxGridPeriodLwir=maxGridPeriodLwir;
			this.debugX=debugX;
			this.debugY=debugY;
			this.debugRadius=debugRadius;
			this.use_large_cells = use_large_cells;
			this.phaseCoeff =      phaseCoeff;
			this.lowpass_sigma =   lowpass_sigma;
			this.min_frac =        min_frac;
			this.min_sin =         min_sin;
			this.no_crazy =        no_crazy;
		}

		public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"gaussWidth",this.gaussWidth+"");
			properties.setProperty(prefix+"corrGamma",this.corrGamma+"");
			properties.setProperty(prefix+"corrSigma",this.corrSigma+"");
			properties.setProperty(prefix+"diffSpectrCorr",this.diffSpectrCorr+"");
			properties.setProperty(prefix+"shrinkClusters",this.shrinkClusters+"");
			properties.setProperty(prefix+"multiplesToTry",this.multiplesToTry+"");
			properties.setProperty(prefix+"deviation",this.deviation+"");
			properties.setProperty(prefix+"deviationSteps",this.deviationSteps+"");
			properties.setProperty(prefix+"highpass",this.highpass+"");
			properties.setProperty(prefix+"corrRingWidth",this.corrRingWidth+"");
			properties.setProperty(prefix+"minCorrContrast",this.minCorrContrast+"");
			properties.setProperty(prefix+"minGridPeriod",this.minGridPeriod+"");
			properties.setProperty(prefix+"maxGridPeriod",this.maxGridPeriod+"");
			properties.setProperty(prefix+"minGridPeriodLwir",this.minGridPeriodLwir+"");
			properties.setProperty(prefix+"maxGridPeriodLwir",this.maxGridPeriodLwir+"");
			properties.setProperty(prefix+"debugX",this.debugX+"");
			properties.setProperty(prefix+"debugY",this.debugY+"");
			properties.setProperty(prefix+"debugRadius",this.debugRadius+"");

			properties.setProperty(prefix+"use_large_cells",this.use_large_cells+"");
			properties.setProperty(prefix+"phaseCoeff",     this.phaseCoeff+"");
			properties.setProperty(prefix+"lowpass_sigma",  this.lowpass_sigma+"");
			properties.setProperty(prefix+"min_frac",       this.min_frac+"");
			properties.setProperty(prefix+"min_sin",        this.min_sin+"");
			properties.setProperty(prefix+"no_crazy",       this.no_crazy+"");
		}

		public void getProperties(String prefix,Properties properties){
			this.gaussWidth=Double.parseDouble(properties.getProperty(prefix+"gaussWidth"));
			this.corrGamma=Double.parseDouble(properties.getProperty(prefix+"corrGamma"));
			this.corrSigma=Double.parseDouble(properties.getProperty(prefix+"corrSigma"));
			this.diffSpectrCorr=Integer.parseInt(properties.getProperty(prefix+"diffSpectrCorr"));
			this.shrinkClusters=Double.parseDouble(properties.getProperty(prefix+"shrinkClusters"));
			this.multiplesToTry=Integer.parseInt(properties.getProperty(prefix+"multiplesToTry"));
			this.deviation=Double.parseDouble(properties.getProperty(prefix+"deviation"));
			this.deviationSteps=Integer.parseInt(properties.getProperty(prefix+"deviationSteps"));
			this.highpass=Double.parseDouble(properties.getProperty(prefix+"highpass"));
			this.corrRingWidth=Double.parseDouble(properties.getProperty(prefix+"corrRingWidth"));
			this.minCorrContrast=Double.parseDouble(properties.getProperty(prefix+"minCorrContrast"));
			if (properties.getProperty(prefix+"minGridPeriod")!=null)
				this.minGridPeriod=Double.parseDouble(properties.getProperty(prefix+"minGridPeriod"));
			else this.minGridPeriod=0.0;
			if (properties.getProperty(prefix+"maxGridPeriod")!=null)
				this.minGridPeriod=Double.parseDouble(properties.getProperty(prefix+"maxGridPeriod"));
			else this.maxGridPeriod=0.0;

			if (properties.getProperty(prefix+"minGridPeriodLwir")!=null)
				this.minGridPeriodLwir=Double.parseDouble(properties.getProperty(prefix+"minGridPeriodLwir"));
			else this.minGridPeriodLwir=0.0;
			if (properties.getProperty(prefix+"maxGridPeriodLwir")!=null)
				this.minGridPeriodLwir=Double.parseDouble(properties.getProperty(prefix+"maxGridPeriodLwir"));
			else this.maxGridPeriodLwir=0.0;

			if (properties.getProperty(prefix+"debugX")!=null)
				this.debugX=Double.parseDouble(properties.getProperty(prefix+"debugX"));
			if (properties.getProperty(prefix+"debugY")!=null)
				this.debugY=Double.parseDouble(properties.getProperty(prefix+"debugY"));
			if (properties.getProperty(prefix+"debugRadius")!=null)
				this.debugRadius=Double.parseDouble(properties.getProperty(prefix+"debugRadius"));

			if (properties.getProperty(prefix+"use_large_cells")!=null) this.use_large_cells=Boolean.parseBoolean(properties.getProperty(prefix+"use_large_cells"));
			if (properties.getProperty(prefix+"phaseCoeff")!=null)      this.phaseCoeff=Double.parseDouble(properties.getProperty(prefix+"phaseCoeff"));
			if (properties.getProperty(prefix+"lowpass_sigma")!=null)   this.lowpass_sigma=Double.parseDouble(properties.getProperty(prefix+"lowpass_sigma"));
			if (properties.getProperty(prefix+"min_frac")!=null)        this.min_frac=Double.parseDouble(properties.getProperty(prefix+"min_frac"));

			if (properties.getProperty(prefix+"min_sin")!=null)         this.min_sin=Double.parseDouble(properties.getProperty(prefix+"min_sin"));
			if (properties.getProperty(prefix+"no_crazy")!=null)        this.no_crazy=Boolean.parseBoolean(properties.getProperty(prefix+"no_crazy"));
		}
	}

	/* ======================================================================== */

	public static class DistortionParameters {
		public int   correlationSize;            // FFTSize/4
		public int   correlationSizeLwir;
		public int   maximalCorrelationSize;     // FFTSize/2
		public int   maximalCorrelationSizeLwir;
		public double correlationGaussWidth; // 0 - no window, <0 - use Hamming
		public boolean absoluteCorrelationGaussWidth=false; // do not scale correlationGaussWidth when the FFT size is increased
		public int zeros; // leave this number of zeros on the margins of the window (toatal from both sides). If correlationGaussWidth>0 will
		// additionally multiply by Hamming
		public int   FFTSize;
		public int   FFTSize_lwir;
		public int   FFTOverlap;      // 32 used for aberration kernels, former FFT_OVERLAP
		public int   FFTOverlap_lwir; // 4

		public double fftGaussWidth;
		public double phaseCorrelationFraction=1.0; // 1.0 - phase correlation, 0.0 - just cross-correlation
		public double correlationHighPassSigma;
		public double correlationLowPassSigma;
		public double correlationRingWidth; // ring (around r=0.5 dist to opposite corr) width , center circle r=0.5*correlationRingWidth
		public double correlationMaxOffset;     // maximal distance between predicted and actual pattern node
		public double correlationMinContrast;   // minimal contrast for the pattern to pass
		public double correlationMinInitialContrast;   // minimal contrast for the pattern of the center (initial point)
		public double correlationMinAbsoluteContrast;   // minimal contrast for the pattern to pass, does not compensate for low ligt
		public double correlationMinAbsoluteInitialContrast;   // minimal contrast for the pattern of the center (initial point)

		public double scaleFirstPassContrast; // Decrease contrast of cells that are too close to the border to be processed in refinement pass
		public double contrastSelectSigmaCenter; // Gaussian sigma to select correlation centers in pixels, 2.0  (center spot)
		public double contrastSelectSigma; // Gaussian sigma to select correlation centers (fraction of UV period), 0.1
		public double contrastAverageSigma; // Gaussian sigma to average correlation variations (as contrast reference) 0.5

		public int    minimalPatternCluster;       //    minimal pattern cluster size (0 - disable retries)
		public int    minimalPatternClusterLwir;   //    minimal pattern cluster size (0 - disable retries)
		public double scaleMinimalInitialContrast; // increase/decrease minimal contrast if initial cluster is >0 but less than minimalPatternCluster

		public double searchOverlap;         // when searching for grid, step this amount of the FFTSize
		public int    patternSubdiv;
		public double correlationDx; // not saved
		public double correlationDy; // not saved
		public int    gridSize;
		public int loop_debug_level;
		public boolean refineCorrelations;
		public boolean fastCorrelationOnFirstPass;
		public boolean fastCorrelationOnFinalPass;
		public double bPatternSigma; // blur bPattern with this sigma
		public double barraySigma; // blur barray with this sigma, multiplied by subdiv
		public double correlationWeightSigma; // sigma (in pixels) for maximum approximation - UNUSED (other maximum methods)
		public double correlationRadiusScale; // maximal radius to consider, in sigmas (if 0 - use sigma as radius) - UNUSED
		public int    correlationRadius;    // radius (green pixel) of the correlation maximum to use for x/y measurement
		public double correlationThreshold; // fraction of the value of the maximum fro the point to be included in centroid calculation
		public int    correlationSubdiv;    // Total subdivision of the correlation maximum (linear and FFT)
		public int    correlationFFTSubdiv; // Increase density of the correlation using FFT
		public boolean correlationAverageOnRefine; // average position between neighbor samples
		public boolean refineInPlace;       // Update coordinates of the grid points as they are recalculated (false - then update all at once)
		public double averageOrthoDist;     // distance to up/down/right left neighbors (0.5)
		public double averageOrthoWeight;   // weight of 4 ortho neighbors (combined) - 0.4), weight of center -s 1.0-averageOrthoWeight-averageDiagWeight
		public double averageDiagDist;     // distance to diagonal neighbors (projection on x/y) (0.5)
		public double averageDiagWeight;   // weight of 4 diagonal neighbors (combined) - 0.4)
		public boolean useQuadratic;       // use quadratic extrapolation to predict position/wave vectors of a new pixel (false - use linear)
		public boolean removeLast;         // remove outer (unreliable) row of nodes
		public int    numberExtrapolated;  // add this number of extrapolated nodes
		public double extrapolationSigma;  // use instead of the correlationWeightSigma during final extrapolation
		public double minUVSpan;           // Minimal u/v span in correlation window that triggers increase of the correlation FFT size
		public boolean flatFieldCorrection=true; // compensate grid uneven intensity (vignetting, illumination)
		public double flatFieldExtarpolate=1.0;  // extrapolate flat field intensity map (relative to the average grid period)
		public double flatFieldBlur=1.0;   // blur the intensity map (relative to the average grid period)
		public double flatFieldMin=0.1;    // do not try to compensate if intensity less than this part of maximal
		public double flatFieldShrink=1.0;     // Shrink before extrapolating intensity map (relative to the average grid period)
		public double flatFieldExpand=3.0;     // Expand during extrapolation (relative to the average grid period)
		public double flatFieldSigmaRadius=1.0;// Extrapolation weight effective radius (relative to the average grid period)
		public double flatFieldExtraRadius=1.5;// Consider pixels in a square with the side twice this (relative to flatFieldSigmaRadius)
		public double averagingAreaScale=  2.0;  // multiply the average grid period to determine the area for averaging the grig brightness


		//        match pointers errors
		public int errTooFewCells=    -10;
		public int errPatternNotFound=-11;
		public boolean legacyMode=false;   // legacy mode

		public DistortionParameters(
				int correlationSize,
				int correlationSizeLwir,
				int maximalCorrelationSize,
				int maximalCorrelationSizeLwir,
				double correlationGaussWidth,
				boolean absoluteCorrelationGaussWidth,
				int zeros,
				int FFTSize,
				int FFTSize_lwir,
				int FFTOverlap,
				int FFTOverlap_lwir,
				double fftGaussWidth,
				double phaseCorrelationFraction,
				double correlationHighPassSigma,
				double correlationLowPassSigma,
				double correlationRingWidth,
				double correlationMaxOffset,     // maximal distance between predicted and actual pattern node
				double correlationMinContrast,   // minimal contrast for the pattern to pass
				double correlationMinInitialContrast,   // minimal contrast for the pattern of the center (initial point)
				double correlationMinAbsoluteContrast,   // minimal contrast for the pattern to pass, does not compensate for low ligt
				double correlationMinAbsoluteInitialContrast,   // minimal contrast for the pattern of the center (initial point)

				double scaleFirstPassContrast, // Decrease contrast of cells that are too close to the border to be processed in refinement pass
				double contrastSelectSigmaCenter, // Gaussian sigma to select correlation centers (fraction of UV period), 0.02 (center spot)
				double contrastSelectSigma, // Gaussian sigma to select correlation centers (fraction of UV period), 0.1
				double contrastAverageSigma, // Gaussian sigma to average correlation variations (as contrast reference) 0.5
				int    minimalPatternCluster,       //    minimal pattern cluster size (0 - disable retries)
				int    minimalPatternClusterLwir,   //    minimal pattern cluster size (0 - disable retries)
				double scaleMinimalInitialContrast, // increase/decrease minimal contrast if initial cluster is >0 but less than minimalPatternCluster
				double searchOverlap,         // when searching for grid, step this amount of the FFTSize
				int patternSubdiv,
				double correlationDx,
				double correlationDy,
				int gridSize,
				int loop_debug_level,
				boolean refineCorrelations,
				boolean fastCorrelationOnFirstPass, // use fast (less precise) correlation on first pass
				boolean fastCorrelationOnFinalPass, // use fast (less precise) correlation on refine pass
				double bPatternSigma,          // blur bPattern with this sigma
				double barraySigma,            // blur barray with this sigma, multiplied by subdiv
				double correlationWeightSigma, // sigma (in pixels) for maximum approximation
				double correlationRadiusScale, // maximal radius to consider, in sigmas (if 0 - use sigma as radius)
				int    correlationRadius,      // radius (green pixel) of the correlation maximum to use for x/y measurement
				double correlationThreshold,   // fraction of the value of the maximum fro the point to be included in centroid calculation
				int    correlationSubdiv,      // Total subdivision of the correlation maximum (linear and FFT)
				int    correlationFFTSubdiv,
				boolean correlationAverageOnRefine, // average position between neighbor samples
				boolean refineInPlace,         // Update coordinates of the grid points as they are recalculated (false - then update all at once)
				double averageOrthoDist,       // distance to up/down/right left neighbors (0.5)
				double averageOrthoWeight,     // weight of 4 ortho neighbors (combined) - 0.4), weight of center -s 1.0-averageOrthoWeight-averageDiagWeight
				double averageDiagDist,        // distance to diagonal neighbors (projection on x/y) (0.5)
				double averageDiagWeight,      // weight of 4 diagonal neighbors (combined) - 0.4)
				boolean useQuadratic,          // use quadratic extrapolation to predict position/wave vectors of a new pixel (false - use linear)
				boolean removeLast,            // remove outer (unreliable) row of nodes
				int    numberExtrapolated,     // add this number of extrapolated nodes
				double extrapolationSigma,     // use instead of the correlationWeightSigma during final extrapolation
				double minUVSpan,              // Minimal u/v span in correlation window that triggers increase of the correlation FFT size
				boolean flatFieldCorrection,   // compensate grid uneven intensity (vignetting, illumination)
				double flatFieldExtarpolate,   // extrapolate flat field intensity map (relative to the average grid period)
				double flatFieldBlur,          // blur the intensity map (relative to the average grid period)
				double flatFieldMin,           // do not try to compensate if intensity less than this part of maximal
				double flatFieldShrink,        // Shrink before extrapolating intensity map (relative to the average grid period)
				double flatFieldExpand,        // Expand during extrapolation (relative to the average grid period)
				double flatFieldSigmaRadius,   // Extrapolation weight effective radius (relative to the average grid period)
				double flatFieldExtraRadius,   // Consider pixels in a square with the side twice this (relative to flatFieldSigmaRadius)
				double averagingAreaScale,     // multiply the average grid period to determine the area for averaging the grig brightness
				boolean legacyMode
				){

			this.correlationSize = correlationSize;
			this.correlationSizeLwir = correlationSizeLwir;
			this.maximalCorrelationSize=maximalCorrelationSize;
			this.maximalCorrelationSizeLwir=maximalCorrelationSizeLwir;
			this.correlationGaussWidth = correlationGaussWidth;
			this.absoluteCorrelationGaussWidth=absoluteCorrelationGaussWidth;
			this.zeros=zeros;
			this.FFTSize = FFTSize;
			this.FFTSize_lwir = FFTSize_lwir;
			this.FFTOverlap = FFTOverlap;
			this.FFTOverlap_lwir = FFTOverlap_lwir;
			this.fftGaussWidth = fftGaussWidth;
			this.phaseCorrelationFraction=phaseCorrelationFraction;
			this.correlationHighPassSigma=correlationHighPassSigma;
			this.correlationLowPassSigma=correlationLowPassSigma;
			this.correlationRingWidth=correlationRingWidth;
			this.correlationMaxOffset=correlationMaxOffset;
			this.correlationMinContrast=correlationMinContrast;
			this.correlationMinInitialContrast=correlationMinInitialContrast;
			this.correlationMinAbsoluteContrast=correlationMinAbsoluteContrast;   // minimal contrast for the pattern to pass, does not compensate for low ligt
			this.correlationMinAbsoluteInitialContrast=correlationMinAbsoluteInitialContrast;   // minimal contrast for the pattern of the center (initial point)
			this.scaleFirstPassContrast=scaleFirstPassContrast; // Decrease contrast of cells that are too close to the border to be processed in refinement pass
			this.contrastSelectSigmaCenter = contrastSelectSigmaCenter; // Gaussian sigma to select correlation centers (pixels, 2.0)
			this.contrastSelectSigma=contrastSelectSigma; // Gaussian sigma to select correlation centers (fraction of UV period), 0.1
			this.contrastAverageSigma=contrastAverageSigma; // Gaussian sigma to average correlation variations (as contrast reference) 0.5
			this.minimalPatternCluster=minimalPatternCluster;        //    minimal pattern cluster size (0 - disable retries)
			this.minimalPatternClusterLwir=minimalPatternClusterLwir;        //    minimal pattern cluster size (0 - disable retries)
			this.scaleMinimalInitialContrast=scaleMinimalInitialContrast; // increase/decrease minimal contrast if initial cluster is >0 but less than minimalPatternCluster
			this.searchOverlap=searchOverlap;         // when searching for grid, step this amount of the FFTSize
			this.patternSubdiv=patternSubdiv;
			this.correlationDx=correlationDx;
			this.correlationDy=correlationDy;
			this.gridSize=gridSize;
			this.loop_debug_level=loop_debug_level;
			this.refineCorrelations=refineCorrelations;
			this.fastCorrelationOnFirstPass=fastCorrelationOnFirstPass;
			this.fastCorrelationOnFinalPass=fastCorrelationOnFinalPass;
			this.bPatternSigma=bPatternSigma; // overwrites SimulationParameters.bPatternSigma
			this.barraySigma=barraySigma;
			this.correlationWeightSigma=correlationWeightSigma;
			this.correlationRadiusScale=correlationRadiusScale;
			this.correlationRadius=correlationRadius;
			this.correlationThreshold=correlationThreshold;
			this.correlationSubdiv=correlationSubdiv;
			this.correlationFFTSubdiv=correlationFFTSubdiv;
			this.correlationAverageOnRefine=correlationAverageOnRefine;
			this.refineInPlace=refineInPlace;
			this.averageOrthoDist=averageOrthoDist;
			this.averageOrthoWeight=averageOrthoWeight;
			this.averageDiagDist=averageDiagDist;
			this.averageDiagWeight=averageDiagWeight;
			this.useQuadratic=useQuadratic;
			this.removeLast=removeLast;
			this.numberExtrapolated=numberExtrapolated;
			this.extrapolationSigma=extrapolationSigma;
			this.minUVSpan=minUVSpan;
			this.flatFieldCorrection=flatFieldCorrection;   // compensate grid uneven intensity (vignetting, illumination)
			this.flatFieldExtarpolate=flatFieldExtarpolate;   // extrapolate flat field intensity map (relative to the average grid period)
			this.flatFieldBlur=flatFieldBlur;          // blur the intensity map (relative to the average grid period)
			this.flatFieldMin=flatFieldMin;
			this.flatFieldShrink=flatFieldShrink;
			this.flatFieldExpand=flatFieldExpand;
			this.flatFieldSigmaRadius=flatFieldSigmaRadius;
			this.flatFieldExtraRadius=flatFieldExtraRadius;
			this.averagingAreaScale=averagingAreaScale;
			this.legacyMode=legacyMode;

		}
		@Override
		public DistortionParameters clone() {
			return new DistortionParameters(
					this.correlationSize,
					this.correlationSizeLwir,
					this.maximalCorrelationSize,
					this.maximalCorrelationSizeLwir,
					this.correlationGaussWidth,
					this.absoluteCorrelationGaussWidth,
					this.zeros,
					this.FFTSize,
					this.FFTSize_lwir,
					this.FFTOverlap,
					this.FFTOverlap_lwir,
					this.fftGaussWidth,
					this.phaseCorrelationFraction,
					this.correlationHighPassSigma,
					this.correlationLowPassSigma,
					this.correlationRingWidth,
					this.correlationMaxOffset,     // maximal distance between predicted and actual pattern node
					this.correlationMinContrast,   // minimal contrast for the pattern to pass
					this.correlationMinInitialContrast,
					this.correlationMinAbsoluteContrast,   // minimal contrast for the pattern to pass, does not compensate for low ligt
					this.correlationMinAbsoluteInitialContrast,   // minimal contrast for the pattern of the center (initial point)
					this.scaleFirstPassContrast, // Decrease contrast of cells that are too close to the border to be processed in refinement pass
					this.contrastSelectSigmaCenter, // Gaussian sigma to select correlation centers (pixels, 2.0)
					this.contrastSelectSigma, // Gaussian sigma to select correlation centers (fraction of UV period), 0.1
					this.contrastAverageSigma, // Gaussian sigma to average correlation variations (as contrast reference) 0.5
					this.minimalPatternCluster,        //    minimal pattern cluster size (0 - disable retries)
					this.minimalPatternClusterLwir,        //    minimal pattern cluster size (0 - disable retries)
					this.scaleMinimalInitialContrast,  // increase/decrease minimal contrast if initial cluster is >0 but less than minimalPatternCluster
					this.searchOverlap,         // when searching for grid, step this amount of the FFTSize
					this.patternSubdiv,
					this.correlationDx,
					this.correlationDy,
					this.gridSize,
					this.loop_debug_level,
					this.refineCorrelations,
					this.fastCorrelationOnFirstPass, // use fast (less precise) correlation on first pass
					this.fastCorrelationOnFinalPass, // use fast (less precise) correlation on refine pass
					this.bPatternSigma, // blur bPattern with this sigma
					this.barraySigma,
					this.correlationWeightSigma, // sigma (in pixels) for maximum approximation
					this.correlationRadiusScale, // maximal radius to consider, in sigmas (if 0 - use sigma as radius)
					this.correlationRadius,    // radius (green pixel) of the correlation maximum to use for x/y measurement
					this.correlationThreshold,
					this.correlationSubdiv,    // Total subdivision of the correlation maximum (linear and FFT)
					this.correlationFFTSubdiv,
					this.correlationAverageOnRefine, // average position between neighbor samples
					this.refineInPlace,       // Update coordinates of the grid points as they are recalculated (false - then update all at once)
					this.averageOrthoDist,    // distance to up/down/right left neighbors (0.5)
					this.averageOrthoWeight,  // weight of 4 ortho neighbors (combined) - 0.4), weight of center -s 1.0-averageOrthoWeight-averageDiagWeight
					this.averageDiagDist,     // distance to diagonal neighbors (projection on x/y) (0.5)
					this.averageDiagWeight,   // weight of 4 diagonal neighbors (combined) - 0.4)
					this.useQuadratic,        // use quadratic extrapolation to predict position/wave vectors of a new pixel (false - use linear)
					this.removeLast,          // remove outer (unreliable) row of nodes
					this.numberExtrapolated,  // add this number of extrapolated nodes
					this.extrapolationSigma,   // use instead of the correlationWeightSigma during final extrapolation
					this.minUVSpan,            // Minimal u/v span in correlation window that triggers increase of the correlation FFT size
					this.flatFieldCorrection,  // compensate grid uneven intensity (vignetting, illumination)
					this.flatFieldExtarpolate, // extrapolate flat field intensity map (relative to the average grid period)
					this.flatFieldBlur,        // blur the intensity map (relative to the average grid period)
					this.flatFieldMin,
					this.flatFieldShrink,
					this.flatFieldExpand,
					this.flatFieldSigmaRadius,
					this.flatFieldExtraRadius,
					this.averagingAreaScale,
					this.legacyMode
					);
		}

		public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"correlationSize",this.correlationSize+"");
			properties.setProperty(prefix+"correlationSizeLwir",this.correlationSizeLwir+"");
			properties.setProperty(prefix+"maximalCorrelationSize",this.maximalCorrelationSize+"");
			properties.setProperty(prefix+"maximalCorrelationSizeLwir",this.maximalCorrelationSizeLwir+"");
			properties.setProperty(prefix+"correlationGaussWidth",this.correlationGaussWidth+"");
			properties.setProperty(prefix+"absoluteCorrelationGaussWidth",this.absoluteCorrelationGaussWidth+"");
			properties.setProperty(prefix+"zeros",this.zeros+"");
			properties.setProperty(prefix+"FFTSize",this.FFTSize+"");
			properties.setProperty(prefix+"FFTSize_lwir",this.FFTSize_lwir+"");
			properties.setProperty(prefix+"FFTOverlap",this.FFTOverlap+"");
			properties.setProperty(prefix+"FFTOverlap_lwir",this.FFTOverlap_lwir+"");
			properties.setProperty(prefix+"fftGaussWidth",this.fftGaussWidth+"");
			properties.setProperty(prefix+"phaseCorrelationFraction",this.phaseCorrelationFraction+"");
			properties.setProperty(prefix+"correlationHighPassSigma",this.correlationHighPassSigma+"");
			properties.setProperty(prefix+"correlationLowPassSigma",this.correlationLowPassSigma+"");
			properties.setProperty(prefix+"correlationRingWidth",this.correlationRingWidth+"");
			properties.setProperty(prefix+"correlationMaxOffset",this.correlationMaxOffset+"");
			properties.setProperty(prefix+"correlationMinContrast",this.correlationMinContrast+"");
			properties.setProperty(prefix+"correlationMinInitialContrast",this.correlationMinInitialContrast+"");
			properties.setProperty(prefix+"correlationMinAbsoluteContrast",this.correlationMinAbsoluteContrast+"");
			properties.setProperty(prefix+"correlationMinAbsoluteInitialContrast",this.correlationMinAbsoluteInitialContrast+"");
			properties.setProperty(prefix+"scaleFirstPassContrast",this.scaleFirstPassContrast+"");
			properties.setProperty(prefix+"contrastSelectSigmaCenter",this.contrastSelectSigmaCenter+"");
			properties.setProperty(prefix+"contrastSelectSigma",this.contrastSelectSigma+"");
			properties.setProperty(prefix+"contrastAverageSigma",this.contrastAverageSigma+"");
			properties.setProperty(prefix+"minimalPatternCluster",this.minimalPatternCluster+"");
			properties.setProperty(prefix+"minimalPatternClusterLwir",this.minimalPatternClusterLwir+"");
			properties.setProperty(prefix+"scaleMinimalInitialContrast",this.scaleMinimalInitialContrast+"");
			properties.setProperty(prefix+"searchOverlap",this.searchOverlap+"");
			properties.setProperty(prefix+"patternSubdiv",this.patternSubdiv+"");
			properties.setProperty(prefix+"correlationDx",this.correlationDx+"");
			properties.setProperty(prefix+"correlationDy",this.correlationDy+"");
			properties.setProperty(prefix+"gridSize",this.gridSize+"");
			properties.setProperty(prefix+"loop_debug_level",this.loop_debug_level+"");
			properties.setProperty(prefix+"refineCorrelations",this.refineCorrelations+"");
			properties.setProperty(prefix+"fastCorrelationOnFirstPass",this.fastCorrelationOnFirstPass+"");
			properties.setProperty(prefix+"fastCorrelationOnFinalPass",this.fastCorrelationOnFinalPass+"");
			properties.setProperty(prefix+"bPatternSigma",this.bPatternSigma+"");
			properties.setProperty(prefix+"barraySigma",this.barraySigma+"");
			properties.setProperty(prefix+"correlationWeightSigma",this.correlationWeightSigma+"");
			properties.setProperty(prefix+"correlationRadiusScale",this.correlationRadiusScale+"");
			properties.setProperty(prefix+"correlationRadius",this.correlationRadius+"");
			properties.setProperty(prefix+"correlationThreshold",this.correlationThreshold+"");
			properties.setProperty(prefix+"correlationSubdiv",this.correlationSubdiv+"");
			properties.setProperty(prefix+"correlationFFTSubdiv",this.correlationFFTSubdiv+"");
			properties.setProperty(prefix+"correlationAverageOnRefine",this.correlationAverageOnRefine+"");
			properties.setProperty(prefix+"refineInPlace",this.refineInPlace+"");
			properties.setProperty(prefix+"averageOrthoDist",this.averageOrthoDist+"");
			properties.setProperty(prefix+"averageOrthoWeight",this.averageOrthoWeight+"");
			properties.setProperty(prefix+"averageDiagDist",this.averageDiagDist+"");
			properties.setProperty(prefix+"averageDiagWeight",this.averageDiagWeight+"");
			properties.setProperty(prefix+"useQuadratic",this.useQuadratic+"");
			properties.setProperty(prefix+"removeLast",this.removeLast+"");
			properties.setProperty(prefix+"numberExtrapolated",this.numberExtrapolated+"");
			properties.setProperty(prefix+"extrapolationSigma",this.extrapolationSigma+"");
			properties.setProperty(prefix+"minUVSpan",this.minUVSpan+"");

			properties.setProperty(prefix+"flatFieldCorrection",this.flatFieldCorrection+"");
			properties.setProperty(prefix+"flatFieldExtarpolate",this.flatFieldExtarpolate+"");
			properties.setProperty(prefix+"flatFieldBlur",this.flatFieldBlur+"");
			properties.setProperty(prefix+"flatFieldMin",this.flatFieldMin+"");

			properties.setProperty(prefix+"flatFieldShrink",this.flatFieldShrink+"");
			properties.setProperty(prefix+"flatFieldExpand",this.flatFieldExpand+"");
			properties.setProperty(prefix+"flatFieldSigmaRadius",this.flatFieldSigmaRadius+"");
			properties.setProperty(prefix+"flatFieldExtraRadius",this.flatFieldExtraRadius+"");
			properties.setProperty(prefix+"averagingAreaScale",this.averagingAreaScale+"");
			properties.setProperty(prefix+"legacyMode",this.minUVSpan+"");
		}
		public void getProperties(String prefix,Properties properties){
//			EProperties properties = (EProperties) pproperties;
			if (properties.getProperty(prefix+"correlationSize")!=null)
				this.correlationSize=Integer.parseInt(properties.getProperty(prefix+"correlationSize"));
			if (properties.getProperty(prefix+"correlationSizeLwir")!=null)
				this.correlationSizeLwir=Integer.parseInt(properties.getProperty(prefix+"correlationSizeLwir"));
			if (properties.getProperty(prefix+"maximalCorrelationSize")!=null)
				this.maximalCorrelationSize=Integer.parseInt(properties.getProperty(prefix+"maximalCorrelationSize"));
			if (properties.getProperty(prefix+"maximalCorrelationSizeLwir")!=null)
				this.maximalCorrelationSizeLwir=Integer.parseInt(properties.getProperty(prefix+"maximalCorrelationSizeLwir"));
			if (properties.getProperty(prefix+"correlationGaussWidth")!=null)
				this.correlationGaussWidth=Double.parseDouble(properties.getProperty(prefix+"correlationGaussWidth"));
			if (properties.getProperty(prefix+"FFTSize")!=null)
				this.FFTSize=Integer.parseInt(properties.getProperty(prefix+"FFTSize"));
			if (properties.getProperty(prefix+"FFTSize_lwir")!=null)
				this.FFTSize_lwir=Integer.parseInt(properties.getProperty(prefix+"FFTSize_lwir"));

			if (properties.getProperty(prefix+"FFTOverlap")!=null)
				this.FFTOverlap=Integer.parseInt(properties.getProperty(prefix+"FFTOverlap"));
			if (properties.getProperty(prefix+"FFTOverlap_lwir")!=null)
				this.FFTOverlap_lwir=Integer.parseInt(properties.getProperty(prefix+"FFTOverlap_lwir"));


// finally shortened :
//			this.FFTOverlap=      properties.getProperty(prefix+"FFTOverlap",     this.FFTOverlap);
//			this.FFTOverlap_lwir= properties.getProperty(prefix+"FFTOverlap_lwir",this.FFTOverlap_lwir);

			if (properties.getProperty(prefix+"absoluteCorrelationGaussWidth")!=null)
				this.absoluteCorrelationGaussWidth=Boolean.parseBoolean(properties.getProperty(prefix+"absoluteCorrelationGaussWidth"));
			if (properties.getProperty(prefix+"zeros")!=null)
				this.zeros=Integer.parseInt(properties.getProperty(prefix+"zeros"));
			if (properties.getProperty(prefix+"fftGaussWidth")!=null)
				this.fftGaussWidth=Double.parseDouble(properties.getProperty(prefix+"fftGaussWidth"));
			if (properties.getProperty(prefix+"phaseCorrelationFraction")!=null)
				this.phaseCorrelationFraction=Double.parseDouble(properties.getProperty(prefix+"phaseCorrelationFraction"));
			if (properties.getProperty(prefix+"correlationHighPassSigma")!=null)
				this.correlationHighPassSigma=Double.parseDouble(properties.getProperty(prefix+"correlationHighPassSigma"));
			if (properties.getProperty(prefix+"correlationLowPassSigma")!=null)
				this.correlationLowPassSigma=Double.parseDouble(properties.getProperty(prefix+"correlationLowPassSigma"));
			if (properties.getProperty(prefix+"correlationRingWidth")!=null)
				this.correlationRingWidth=Double.parseDouble(properties.getProperty(prefix+"correlationRingWidth"));
			if (properties.getProperty(prefix+"correlationMaxOffset")!=null)
				this.correlationMaxOffset=Double.parseDouble(properties.getProperty(prefix+"correlationMaxOffset"));
			if (properties.getProperty(prefix+"correlationMinContrast")!=null)
				this.correlationMinContrast=Double.parseDouble(properties.getProperty(prefix+"correlationMinContrast"));
			if (properties.getProperty(prefix+"correlationMinInitialContrast")!=null)
				this.correlationMinInitialContrast=Double.parseDouble(properties.getProperty(prefix+"correlationMinInitialContrast"));

			if (properties.getProperty(prefix+"correlationMinAbsoluteContrast")!=null)
				this.correlationMinAbsoluteContrast=Double.parseDouble(properties.getProperty(prefix+"correlationMinAbsoluteContrast"));
			if (properties.getProperty(prefix+"correlationMinAbsoluteInitialContrast")!=null)
				this.correlationMinAbsoluteInitialContrast=Double.parseDouble(properties.getProperty(prefix+"correlationMinAbsoluteInitialContrast"));

			if (properties.getProperty(prefix+"scaleFirstPassContrast")!=null)
				this.scaleFirstPassContrast=Double.parseDouble(properties.getProperty(prefix+"scaleFirstPassContrast"));
			if (properties.getProperty(prefix+"contrastSelectSigmaCenter")!=null)
				this.contrastSelectSigmaCenter=Double.parseDouble(properties.getProperty(prefix+"contrastSelectSigmaCenter"));
			if (properties.getProperty(prefix+"contrastSelectSigma")!=null)
				this.contrastSelectSigma=Double.parseDouble(properties.getProperty(prefix+"contrastSelectSigma"));
			if (properties.getProperty(prefix+"contrastAverageSigma")!=null)
				this.contrastAverageSigma=Double.parseDouble(properties.getProperty(prefix+"contrastAverageSigma"));

			if (properties.getProperty(prefix+"minimalPatternCluster")!=null)
				this.minimalPatternCluster=Integer.parseInt(properties.getProperty(prefix+"minimalPatternCluster"));
			if (properties.getProperty(prefix+"minimalPatternClusterLwir")!=null)
				this.minimalPatternClusterLwir=Integer.parseInt(properties.getProperty(prefix+"minimalPatternClusterLwir"));

			if (properties.getProperty(prefix+"scaleMinimalInitialContrast")!=null)
				this.scaleMinimalInitialContrast=Double.parseDouble(properties.getProperty(prefix+"scaleMinimalInitialContrast"));
			if (properties.getProperty(prefix+"searchOverlap")!=null)
				this.searchOverlap=Double.parseDouble(properties.getProperty(prefix+"searchOverlap"));
			if (properties.getProperty(prefix+"patternSubdiv")!=null)
				this.patternSubdiv=Integer.parseInt(properties.getProperty(prefix+"patternSubdiv"));
			if (properties.getProperty(prefix+"correlationDx")!=null)
				this.correlationDx=Double.parseDouble(properties.getProperty(prefix+"correlationDx"));
			if (properties.getProperty(prefix+"correlationDy")!=null)
				this.correlationDy=Double.parseDouble(properties.getProperty(prefix+"correlationDy"));
			if (properties.getProperty(prefix+"gridSize")!=null)
				this.gridSize=Integer.parseInt(properties.getProperty(prefix+"gridSize"));
			if (properties.getProperty(prefix+"loop_debug_level")!=null)
				this.loop_debug_level=Integer.parseInt(properties.getProperty(prefix+"loop_debug_level"));
			if (properties.getProperty(prefix+"refineCorrelations")!=null)
				this.refineCorrelations=Boolean.parseBoolean(properties.getProperty(prefix+"refineCorrelations"));
			if (properties.getProperty(prefix+"fastCorrelationOnFirstPass")!=null)
				this.fastCorrelationOnFirstPass=Boolean.parseBoolean(properties.getProperty(prefix+"fastCorrelationOnFirstPass"));
			if (properties.getProperty(prefix+"fastCorrelationOnFinalPass")!=null)
				this.fastCorrelationOnFinalPass=Boolean.parseBoolean(properties.getProperty(prefix+"fastCorrelationOnFinalPass"));
			if (properties.getProperty(prefix+"bPatternSigma")!=null)
				this.bPatternSigma=Double.parseDouble(properties.getProperty(prefix+"bPatternSigma"));
			if (properties.getProperty(prefix+"barraySigma")!=null)
				this.barraySigma=Double.parseDouble(properties.getProperty(prefix+"barraySigma"));
			if (properties.getProperty(prefix+"correlationWeightSigma")!=null)
				this.correlationWeightSigma=Double.parseDouble(properties.getProperty(prefix+"correlationWeightSigma"));
			if (properties.getProperty(prefix+"correlationRadiusScale")!=null)
				this.correlationRadiusScale=Double.parseDouble(properties.getProperty(prefix+"correlationRadiusScale"));
			if (properties.getProperty(prefix+"correlationRadius")!=null)
				this.correlationRadius=Integer.parseInt(properties.getProperty(prefix+"correlationRadius"));
			if (properties.getProperty(prefix+"correlationThreshold")!=null)
				this.correlationThreshold=Double.parseDouble(properties.getProperty(prefix+"correlationThreshold"));
			if (properties.getProperty(prefix+"correlationSubdiv")!=null)
				this.correlationSubdiv=Integer.parseInt(properties.getProperty(prefix+"correlationSubdiv"));
			if (properties.getProperty(prefix+"correlationFFTSubdiv")!=null)
				this.correlationFFTSubdiv=Integer.parseInt(properties.getProperty(prefix+"correlationFFTSubdiv"));
			if (properties.getProperty(prefix+"correlationAverageOnRefine")!=null)
				this.correlationAverageOnRefine=Boolean.parseBoolean(properties.getProperty(prefix+"correlationAverageOnRefine"));
			if (properties.getProperty(prefix+"refineInPlace")!=null)
				this.refineInPlace=Boolean.parseBoolean(properties.getProperty(prefix+"refineInPlace"));
			if (properties.getProperty(prefix+"averageOrthoDist")!=null)
				this.averageOrthoDist=Double.parseDouble(properties.getProperty(prefix+"averageOrthoDist"));
			if (properties.getProperty(prefix+"averageOrthoWeight")!=null)
				this.averageOrthoWeight=Double.parseDouble(properties.getProperty(prefix+"averageOrthoWeight"));
			if (properties.getProperty(prefix+"averageDiagDist")!=null)
				this.averageDiagDist=Double.parseDouble(properties.getProperty(prefix+"averageDiagDist"));
			if (properties.getProperty(prefix+"correlationRadiusScale")!=null)
				this.averageDiagWeight=Double.parseDouble(properties.getProperty(prefix+"averageDiagWeight"));
			if (properties.getProperty(prefix+"useQuadratic")!=null)
				this.useQuadratic=Boolean.parseBoolean(properties.getProperty(prefix+"useQuadratic"));
			if (properties.getProperty(prefix+"removeLast")!=null)
				this.removeLast=Boolean.parseBoolean(properties.getProperty(prefix+"removeLast"));
			if (properties.getProperty(prefix+"numberExtrapolated")!=null)
				this.numberExtrapolated=Integer.parseInt(properties.getProperty(prefix+"numberExtrapolated"));
			if (properties.getProperty(prefix+"extrapolationSigma")!=null)
				this.extrapolationSigma=Double.parseDouble(properties.getProperty(prefix+"extrapolationSigma"));
			if (properties.getProperty(prefix+"minUVSpan")!=null)
				this.minUVSpan=Double.parseDouble(properties.getProperty(prefix+"minUVSpan"));
			if (properties.getProperty(prefix+"flatFieldCorrection")!=null)
				this.flatFieldCorrection=Boolean.parseBoolean(properties.getProperty(prefix+"flatFieldCorrection"));
			if (properties.getProperty(prefix+"flatFieldExtarpolate")!=null)
				this.flatFieldExtarpolate=Double.parseDouble(properties.getProperty(prefix+"flatFieldExtarpolate"));
			if (properties.getProperty(prefix+"flatFieldBlur")!=null)
				this.flatFieldBlur=Double.parseDouble(properties.getProperty(prefix+"flatFieldBlur"));
			if (properties.getProperty(prefix+"flatFieldMin")!=null)
				this.flatFieldMin=Double.parseDouble(properties.getProperty(prefix+"flatFieldMin"));
			if (properties.getProperty(prefix+"flatFieldShrink")!=null)
				this.flatFieldShrink=Double.parseDouble(properties.getProperty(prefix+"flatFieldShrink"));
			if (properties.getProperty(prefix+"flatFieldExpand")!=null)
				this.flatFieldExpand=Double.parseDouble(properties.getProperty(prefix+"flatFieldExpand"));
			if (properties.getProperty(prefix+"flatFieldSigmaRadius")!=null)
				this.flatFieldSigmaRadius=Double.parseDouble(properties.getProperty(prefix+"flatFieldSigmaRadius"));
			if (properties.getProperty(prefix+"flatFieldExtraRadius")!=null)
				this.flatFieldExtraRadius=Double.parseDouble(properties.getProperty(prefix+"flatFieldExtraRadius"));
			if (properties.getProperty(prefix+"averagingAreaScale")!=null)
				this.averagingAreaScale=Double.parseDouble(properties.getProperty(prefix+"averagingAreaScale"));
			if (properties.getProperty(prefix+"legacyMode")!=null)
				this.legacyMode=Boolean.parseBoolean(properties.getProperty(prefix+"legacyMode"));
		}
	}
	///===== end of public static class DistortionParameters ==============================
	/* Use ROI */
	/* Supply rectangle */

	// Now accepts rectangles not completely contained in the image, pixels will be copied from the image edge
	// getNoBayer() - replacement for splitBayer for monochrome images
	public double[] getNoBayer(ImagePlus imp, Rectangle r) {
		return getNoBayer (imp, 1, r);
	}
	//	private double[][] splitBayer (ImagePlus imp, Rectangle r, boolean equalize_greens) {
	public double[][] splitBayer (ImagePlus imp, Rectangle r, boolean equalize_greens) {
		return splitBayer (imp, 1, r, equalize_greens);
	}
	public double[] getNoBayer (ImagePlus imp,  int sliceNumber, Rectangle r) {
		if (imp==null) return null;
		ImageProcessor ip=null;
		float [] pixels;
		if (imp.getStackSize()>1){
			ip=imp.getStack().getProcessor(sliceNumber);
		} else {
			ip=imp.getProcessor();
		}
		pixels=(float[])ip.getPixels();
		int full_width= imp.getWidth();  // full image width
		int full_height=imp.getHeight(); // full image height
		if (r==null) r=new Rectangle(0,0,full_width,full_height);
		double [] dpixels = new double[r.height * r.width];
		int out_indx = 0;
		for (int y = 0; y < r.height; y++) {
			int y_in = y + r.y;
			if      (y_in < 0)            y_in = 0;
			else if (y_in >= full_height) y_in = full_height - 1;
			int base = y_in * full_width;
			for (int x = 0; x < r.width; x++) {
				int x_in = x + r.x;
				if      (x_in < 0)           x_in = 0;
				else if (x_in >= full_width) x_in = full_width - 1;
				dpixels[out_indx++] = pixels[base+x_in];
			}
		}
		return dpixels;
	}

	public double[][] splitBayer (ImagePlus imp,  int sliceNumber, Rectangle r, boolean equalize_greens) {
		if (imp==null) return null;
		ImageProcessor ip=null;
		float [] pixels;
		if (imp.getStackSize()>1){
			ip=imp.getStack().getProcessor(sliceNumber);
		} else {
			ip=imp.getProcessor();
		}
		pixels=(float[])ip.getPixels();   // null pointer
		int full_width= imp.getWidth();  // full image width
		int full_height=imp.getHeight(); // full image height
		if (r==null) r=new Rectangle(0,0,full_width,full_height);
		if (debugLevel>10) IJ.showMessage("splitBayer","r.width="+r.width+
				"\nr.height="+r.height+
				"\nr.x="+r.x+
				"\nr.y="+r.y+
				"\nlength="+pixels.length);
		if ((debugLevel>2) && ((r.x<0) || (r.y<0) || ((r.x+r.width)>=full_width) || ((r.y+r.height)>=full_height))) System.out.println("r.width="+r.width+
				" r.height="+r.height+
				" r.x="+r.x+
				" r.y="+r.y);
		int x,y,base,base_b,bv,i,j;
		int half_height = (r.height >> 1);
		int half_width =  (r.width >>  1);
		// make them all 0 if not a single pixel falls into the image
		int numColors=(half_height==half_width)?5:4;
		int pixX,pixY;
		double [][] bayer_pixels=new double[numColors][half_height * half_width];
		if ((r.x>=full_width) || (r.y>=full_height) || ((r.x+r.width)<0)  || ((r.y+r.height)<0)) {
			for (i=0;i<bayer_pixels.length;i++) for (j=0;j<bayer_pixels[i].length;j++) bayer_pixels[i][j]=0.0;
			return bayer_pixels;
		}
		//      base=r.width*((y<<1)+bv);
		for (y=0; y<half_height; y++) for (bv=0;bv<2;bv++){
			pixY=(y*2)+bv+r.y;
			base_b=half_width*y;
			//					if ((pixY>=0)

			if (pixY<0) {
				pixY=bv;
			} else if (pixY>=full_height){
				pixY=full_height-2+bv;
			}
			base=full_width*pixY+((r.x>0)?r.x:0);
			//						base=full_width*((y*2)+bv+r.y)+r.x;
			pixX=r.x;
			if (bv==0) for (x=0; x<half_width; x++) {
				if ((pixX<0) || (pixX>=(full_width-2))) {
					bayer_pixels[0][base_b]= pixels[base];
					bayer_pixels[1][base_b]= pixels[base+1];
				} else {
					bayer_pixels[0][base_b]= pixels[base++];
					bayer_pixels[1][base_b]= pixels[base++];
				}
				base_b++;
				pixX+=2;
			} else  for (x=0; x<half_width; x++) {
				if ((pixX<0) || (pixX>=(full_width-2))) {
					bayer_pixels[2][base_b]= pixels[base];
					bayer_pixels[3][base_b]= pixels[base+1];
				} else {
					bayer_pixels[2][base_b]= pixels[base++];
					bayer_pixels[3][base_b]= pixels[base++];
				}
				base_b++;
				pixX+=2;
			}
		}
		if (equalize_greens) {
			double g0=0.0,g3=0.0,g02=0.0,g32=0.0,a0,a3,b0,b3;
			int n=bayer_pixels[0].length;
			for (i=0;i<bayer_pixels[0].length;i++) {
				g0 +=bayer_pixels[0][i];
				g02+=bayer_pixels[0][i]*bayer_pixels[0][i];
				g3 +=bayer_pixels[3][i];
				g32+=bayer_pixels[3][i]*bayer_pixels[3][i];
			}
			g0/=n; // mean value
			g3/=n; // meran value
			g02=g02/n-g0*g0;
			g32=g32/n-g3*g3;
			b0=Math.sqrt(Math.sqrt(g32/g02));
			b3 = 1.0/b0;
			a0= (g0+g3)/2 -b0*g0;
			a3= (g0+g3)/2 -b3*g3;
			if (debugLevel>2) {
				System.out.println("g0= "+g0+ ", g3= "+g3);
				System.out.println("g02="+g02+", g32="+g32);
				System.out.println("a0="+a0+", b0="+b0);
				System.out.println("a3="+a3+", b3="+b3);
			}
			for (i=0;i<bayer_pixels[0].length;i++) {
				bayer_pixels[0][i]=a0+bayer_pixels[0][i]*b0;
				bayer_pixels[3][i]=a3+bayer_pixels[3][i]*b3;
			}

		}

		if (numColors>4) bayer_pixels[4]=combineDiagonalGreens (bayer_pixels[0], bayer_pixels[3],  half_width, half_height);
		return bayer_pixels;
	}

	public double[][] splitBayerOne (ImagePlus imp, Rectangle r, boolean equalize_greens) {
		ImageProcessor ip=imp.getProcessor();
		float [] pixels;
		pixels=(float[])ip.getPixels();
		int full_width= imp.getWidth();  // full image width
		int full_height=imp.getHeight(); // full image height
		if (debugLevel>10) IJ.showMessage("splitBayer","r.width="+r.width+
				"\nr.height="+r.height+
				"\nr.x="+r.x+
				"\nr.y="+r.y+
				"\nlength="+pixels.length);
		if ((debugLevel>2) && ((r.x<0) || (r.y<0) || ((r.x+r.width)>=full_width) || ((r.y+r.height)>=full_height))) System.out.println("r.width="+r.width+
				" r.height="+r.height+
				" r.x="+r.x+
				" r.y="+r.y);
		int x,y,base,base_b,bv,i,j;
		int half_height=r.height>>1;
		int half_width=r.width>>1;
			// make them all 0 if not a single pixel falls into the image
			int numColors=(half_height==half_width)?5:4;
			int pixX,pixY;
			double [][] bayer_pixels=new double[numColors][half_height * half_width];
			if ((r.x>=full_width) || (r.y>=full_height) || ((r.x+r.width)<0)  || ((r.y+r.height)<0)) {
				for (i=0;i<bayer_pixels.length;i++) for (j=0;j<bayer_pixels[i].length;j++) bayer_pixels[i][j]=0.0;
				return bayer_pixels;
			}
			//      base=r.width*((y<<1)+bv);
			for (y=0; y<half_height; y++) for (bv=0;bv<2;bv++){
				pixY=(y*2)+bv+r.y;
				base_b=half_width*y;
				//					if ((pixY>=0)

				if (pixY<0) {
					pixY=bv;
				} else if (pixY>=full_height){
					pixY=full_height-2+bv;
				}
				base=full_width*pixY+((r.x>0)?r.x:0);
				//						base=full_width*((y*2)+bv+r.y)+r.x;
				pixX=r.x;
				if (bv==0) for (x=0; x<half_width; x++) {
					if ((pixX<0) || (pixX>=(full_width-2))) {
						bayer_pixels[0][base_b]= pixels[base];
						bayer_pixels[1][base_b]= pixels[base+1];
					} else {
						bayer_pixels[0][base_b]= pixels[base++];
						bayer_pixels[1][base_b]= pixels[base++];
					}
					base_b++;
					pixX+=2;
				} else  for (x=0; x<half_width; x++) {
					if ((pixX<0) || (pixX>=(full_width-2))) {
						bayer_pixels[2][base_b]= pixels[base];
						bayer_pixels[3][base_b]= pixels[base+1];
					} else {
						bayer_pixels[2][base_b]= pixels[base++];
						bayer_pixels[3][base_b]= pixels[base++];
					}
					base_b++;
					pixX+=2;
				}
			}
			if (equalize_greens) {
				double g0=0.0,g3=0.0,g02=0.0,g32=0.0,a0,a3,b0,b3;
				int n=bayer_pixels[0].length;
				for (i=0;i<bayer_pixels[0].length;i++) {
					g0 +=bayer_pixels[0][i];
					g02+=bayer_pixels[0][i]*bayer_pixels[0][i];
					g3 +=bayer_pixels[3][i];
					g32+=bayer_pixels[3][i]*bayer_pixels[3][i];
				}
				g0/=n; // mean value
				g3/=n; // meran value
				g02=g02/n-g0*g0;
				g32=g32/n-g3*g3;
				b0=Math.sqrt(Math.sqrt(g32/g02));
				b3 = 1.0/b0;
				a0= (g0+g3)/2 -b0*g0;
				a3= (g0+g3)/2 -b3*g3;
				if (debugLevel>2) {
					System.out.println("g0= "+g0+ ", g3= "+g3);
					System.out.println("g02="+g02+", g32="+g32);
					System.out.println("a0="+a0+", b0="+b0);
					System.out.println("a3="+a3+", b3="+b3);
				}
				for (i=0;i<bayer_pixels[0].length;i++) {
					bayer_pixels[0][i]=a0+bayer_pixels[0][i]*b0;
					bayer_pixels[3][i]=a3+bayer_pixels[3][i]*b3;
				}

			}

			if (numColors>4) bayer_pixels[4]=combineDiagonalGreens (bayer_pixels[0], bayer_pixels[3],  half_width, half_height);
			return bayer_pixels;
	}

	public double[][] splitBayerZero (ImagePlus imp, Rectangle r, boolean equalize_greens) {
		ImageProcessor ip=imp.getProcessor();
		float [] pixels;
		pixels=(float[])ip.getPixels();
		if (debugLevel>10) IJ.showMessage("splitBayer","r.width="+r.width+
				"\nr.height="+r.height+
				"\nr.x="+r.x+
				"\nr.y="+r.y+
				"\nlength="+pixels.length);
		int x,y,base,base_b,bv,i;
		int half_height=r.height>>1;
				int half_width=r.width>>1;
				int full_width= imp.getWidth();  // full image width
				int full_height=imp.getHeight(); // full image height
				int numColors=(half_height==half_width)?5:4;
				int pixX,pixY;
				double [][] bayer_pixels=new double[numColors][half_height * half_width];
				//      base=r.width*((y<<1)+bv);
				for (y=0; y<half_height; y++) for (bv=0;bv<2;bv++){
					pixY=(y*2)+bv+r.y;
					base_b=half_width*y;
					if ((pixY<0) || (pixY>=full_height)) {
						if (bv==0) for (x=0; x<half_width; x++) {
							bayer_pixels[0][base_b]= 0.0;
							bayer_pixels[1][base_b]= 0.0;
							base_b++;
						} else  for (x=0; x<half_width; x++) {
							bayer_pixels[2][base_b]= 0.0;
							bayer_pixels[3][base_b]= 0.0;
							base_b++;
						}
					} else {
						base=full_width*((y*2)+bv+r.y)+r.x;
						pixX=r.x;
						if (bv==0) for (x=0; x<half_width; x++) {
							if ((pixX<0) || (pixX>=(full_width-1))) {
								bayer_pixels[0][base_b]= 0.0;
								bayer_pixels[1][base_b]= 0.0;
								base+=2;
							} else {
								bayer_pixels[0][base_b]= pixels[base++];
								bayer_pixels[1][base_b]= pixels[base++];
							}
							base_b++;
							pixX+=2;
						} else  for (x=0; x<half_width; x++) {
							if ((pixX<0) || (pixX>=(full_width-1))) {
								bayer_pixels[2][base_b]= 0.0;
								bayer_pixels[3][base_b]= 0.0;
								base+=2;
							} else {
								bayer_pixels[2][base_b]= pixels[base++];
								bayer_pixels[3][base_b]= pixels[base++];
							}
							base_b++;
							pixX+=2;
						}
					}
				}
				if (equalize_greens) {
					double g0=0.0,g3=0.0,g02=0.0,g32=0.0,a0,a3,b0,b3;
					int n=bayer_pixels[0].length;
					for (i=0;i<bayer_pixels[0].length;i++) {
						g0 +=bayer_pixels[0][i];
						g02+=bayer_pixels[0][i]*bayer_pixels[0][i];
						g3 +=bayer_pixels[3][i];
						g32+=bayer_pixels[3][i]*bayer_pixels[3][i];
					}
					g0/=n; // mean value
					g3/=n; // meran value
					g02=g02/n-g0*g0;
					g32=g32/n-g3*g3;
					b0=Math.sqrt(Math.sqrt(g32/g02));
					b3 = 1.0/b0;
					a0= (g0+g3)/2 -b0*g0;
					a3= (g0+g3)/2 -b3*g3;
					if (debugLevel>2) {
						System.out.println("g0= "+g0+ ", g3= "+g3);
						System.out.println("g02="+g02+", g32="+g32);
						System.out.println("a0="+a0+", b0="+b0);
						System.out.println("a3="+a3+", b3="+b3);
					}
					for (i=0;i<bayer_pixels[0].length;i++) {
						bayer_pixels[0][i]=a0+bayer_pixels[0][i]*b0;
						bayer_pixels[3][i]=a3+bayer_pixels[3][i]*b3;
					}

				}

				if (numColors>4) bayer_pixels[4]=combineDiagonalGreens (bayer_pixels[0], bayer_pixels[3],  half_width, half_height);
				return bayer_pixels;
	}

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
	// Parameters for identifying red laser pointers on the image of the pattern grid

}
