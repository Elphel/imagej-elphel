package com.elphel.imagej.dp;
/**
** -----------------------------------------------------------------------------**
** SFEPhases.java
**
** Detecting SFE phase misalignment
**
** Copyright (C) 2014 Elphel, Inc.
**
** -----------------------------------------------------------------------------**
**  
**  SFEPhases.java is free software: you can redistribute it and/or modify
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
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.io.Opener;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import javax.swing.SwingUtilities;

import com.elphel.imagej.calibration.CalibrationFileManagement;
import com.elphel.imagej.calibration.DistortionProcessConfiguration;
import com.elphel.imagej.calibration.EyesisAberrations;
import com.elphel.imagej.calibration.MatchSimulatedPattern;
import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.jp4.JP46_Reader_camera;

public class SFEPhases {
	public double overexposureThreshold=0.95; // do not consider pixels above this value
	public int tileSize=16; // side of the square tile (each color channel)
	public double binWidth=0.05; // "dark" and "bright" pixels should fit into this wide bins (negative - fraction of max-min)
	public double gapWidth=0.05; // minimal distance between "dark" and "bright" with no pixels in it (negative - fraction of max-min)
	public class Defect{
		public int x;
		public int y;
		public double diff;
		public Defect (int x, int y, double d){
			this.x=x;
			this.y=y;
			this.diff = d;
		}
	}
	public int [] getPhaseStats(
			final ImagePlus imp,
			final int tileHalfSize,
			final int cmask, // bitmask of color channels to process (9 - two greens) 
			final double binWidth, // if negative - relative to max-min
			final double gapWidth,
			final double satLevel, // ignore pixels above this
			final double baseFrac, // fraction of "normal" pixels (0.5)
			final int minPixToProcess, //
			final int numHistBins,
			final boolean forceLinear,  // use linear approximation
			final MatchSimulatedPattern matchSimulatedPattern,
			final EyesisAberrations.ColorComponents colorComponents,
			int threadsMax,
			final int debugLevel){
		final int debugThreshold=0;
		// Get selection from image or full image		
		int imgWidth=imp.getWidth();
		int imgHeight=imp.getHeight();
		Rectangle selection= (imp.getRoi()==null)? (new Rectangle(0, 0, imgWidth, imgHeight)):imp.getRoi().getBounds();
		int numTiles=0;
		for (int i= selection.y/tileHalfSize;i<=((selection.y+selection.height)/tileHalfSize-2);i++){
			int y0=i*tileHalfSize;
			for (int j=selection.x/tileHalfSize; j<=((selection.x+selection.width)/tileHalfSize-2);j++){
				int x0=j*tileHalfSize;
				Rectangle tile=new Rectangle(x0,y0,2*tileHalfSize,2*tileHalfSize);
				if (selection.contains(tile)) numTiles++;
			}
		}
		final int [][] tileList=new int [numTiles][2];
		numTiles=0;
		for (int i= selection.y/tileHalfSize;i<=((selection.y+selection.height)/tileHalfSize-2);i++){
			int y0=i*tileHalfSize;
			for (int j=selection.x/tileHalfSize; j<=((selection.x+selection.width)/tileHalfSize-2);j++){
				int x0=j*tileHalfSize;
				Rectangle tile=new Rectangle(x0,y0,2*tileHalfSize,2*tileHalfSize);
				if (selection.contains(tile)) {
					tileList[numTiles][0]=x0;
					tileList[numTiles][1]=y0;
					numTiles++;
				}
			}
		}
		if (debugLevel>debugThreshold){
			
		}
		final AtomicInteger aNumTile=new AtomicInteger(0);
		final AtomicInteger numOutlayers=new AtomicInteger(0);
		final AtomicInteger numPix=new AtomicInteger(0);
		final Thread[] threads = newThreadArray(threadsMax);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			// Concurrently run in as many threads as CPUs
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = aNumTile.getAndIncrement(); nTile < tileList.length; nTile = aNumTile.getAndIncrement()) {
						int [] tilePhaseStats=getPhaseStatsTile(
								tileList[nTile][0], // int x0, // top left in sensor pixels (even)
								tileList[nTile][1], // int y0,
								imp,
								tileHalfSize,
								cmask, // bitmask of color channels to process (9 - two greens) 
								binWidth, // absolute
								gapWidth, // absolute
								satLevel, // ignore pixels above this
								baseFrac, // fraction of "normal" pixels (0.5)
								minPixToProcess, //
								numHistBins,
								forceLinear,  // use linear approximation
								matchSimulatedPattern,
								colorComponents,
								debugLevel);
						if (tilePhaseStats[0]>0){
							numOutlayers.getAndAdd(tilePhaseStats[0]);
							numPix.getAndAdd(tilePhaseStats[1]);
							if (debugLevel>debugThreshold){
								System.out.println("x0="+tileList[nTile][0]+" y0="+tileList[nTile][1]+" tilePhaseStats[0]="+tilePhaseStats[0]+" tilePhaseStats[1]"+tilePhaseStats[1]);
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
		int [] result = {numOutlayers.get(),numPix.get()};
		return result;
	}
	
	private int [] getPhaseStatsTile(
			int x0, // top left in sensor pixels (even)
			int y0,
			ImagePlus imp,
			int tileHalfSize,
			int cmask, // bitmask of color channels to process (9 - two greens) 
			double binWidth, // absolute
			double gapWidth, // absolute
			double satLevel, // ignore pixels above this
			double baseFrac, // fraction of "normal" pixels (0.5)
			int minPixToProcess, //
			int numHistBins,
			boolean forceLinear,  // use linear approximation
			MatchSimulatedPattern matchSimulatedPattern,
			EyesisAberrations.ColorComponents colorComponents,
			int debugLevel){
		Rectangle tileRect=new Rectangle (x0,y0,2*tileHalfSize,2*tileHalfSize);
		double [][] input_bayer=matchSimulatedPattern.splitBayer(imp,tileRect,colorComponents.equalizeGreens); // does it work trhe same?
		int [] result = {0,0};
		for (int chn=0;chn<input_bayer.length;chn++) if ((input_bayer[chn]!=null) &&((cmask & (1<<chn))!=0)){
			double min=Double.NaN,max=Double.NaN;
			int nPix=0;
			for (int i=0;i<input_bayer[chn].length;i++){
				double d= input_bayer[chn][i];
				if ((satLevel<=0) || (d<satLevel)){
					if (nPix==0){
						min = input_bayer[chn][i];
						max = input_bayer[chn][i];
					} else { 
						if (min>input_bayer[chn][i]) min = input_bayer[chn][i];
						if (max<input_bayer[chn][i]) max = input_bayer[chn][i];
					}
					nPix++;
				}
			}
			if (nPix<minPixToProcess) continue; // not enough non-saturated pixels
			int [] histogram=new int [numHistBins];
			for (int i=0;i<numHistBins;i++)histogram[i]=0;
			int [] histWidth=histogram.clone(); // number of histogram bins, starting from index, having <=baseFrac pixels
			int [] histHeight=histogram.clone(); // number of counts in selected bins
			// build histogram
			for (int i=0;i<input_bayer[chn].length;i++){
				double d= input_bayer[chn][i];
				if ((satLevel<=0) || (d<satLevel)){
					int hi= (int) Math.floor((d-min)*numHistBins/(max-min));
					if (hi<0) hi=0;
					else if (hi>=numHistBins) hi=numHistBins-1;
					histogram[hi]++;
				}
			}
			int iBase= (int) (nPix*baseFrac);
			int startIndex=0, endIndex=0; // last+1)
			int numCounts=0;
			while ((startIndex < numHistBins) && (endIndex<=numHistBins)){
				// move last index to get up to iBase counts total;
				while ((endIndex<numHistBins) && (numCounts<iBase)){
					numCounts+=histogram[endIndex++];
				}
				if ((numCounts>iBase) && (endIndex>(startIndex+1))){ // keep at least 1 bin
					endIndex--;
					numCounts-=histogram[endIndex];
				}
				histWidth[startIndex]=endIndex-startIndex;
				histHeight[startIndex]=numCounts;
				// move start index;
				numCounts-=histogram[startIndex];
				startIndex++;
			}
			int bestHistWidth=0,bestHistHeight=0,bestStartIndex=-1;
			for (int i=0;i<numHistBins;i++) if ((histHeight[i]>0)){
				if ((bestHistWidth==0) || (bestHistWidth>histWidth[i])) {
					bestHistWidth=histWidth[i];
					bestHistHeight=histHeight[i];
					bestStartIndex=i;
				} else if ((bestHistWidth == histWidth[i]) && (bestHistHeight < histHeight[i])){
					bestHistHeight=histHeight[i];
					bestStartIndex=i;
				}
				if ((i+histWidth[i])>=numHistBins) break;
			}
			double baseMin=min+((max-min)*bestStartIndex)/numHistBins;
			double baseMax=min+((max-min)*(bestStartIndex+histWidth[bestStartIndex]))/numHistBins;
			// 2-d linear interpolate "base" (presumably normal) pixels)
			int numBase=0;
			for (int i=1;i<input_bayer[chn].length;i++){
				double d= input_bayer[chn][i];
				if ((d >= baseMin) && (d < baseMax)){
					numBase++;
				}
			}
			if (numBase<((int) (minPixToProcess*baseFrac))){
				continue; // too few points left to interpolate 
			}
			double [][][] data = new double [numBase][2][];
			numBase=0;
			for (int i=1;i<input_bayer[chn].length;i++){
				double d= input_bayer[chn][i];
				if ((d >= baseMin) && (d < baseMax)){
					data[numBase][0]=new double [2];
					data[numBase][0][0]=i % tileHalfSize; // x
					data[numBase][0][1]=i / tileHalfSize; // y
					data[numBase][1]=new double [1];
					data[numBase][1][0]=d; // f()x,y)
					numBase++;
				}
			}
			double[][] polyCoeff=(new PolynomialApproximation()).quadraticApproximation(
					data,
					forceLinear);
			if ((polyCoeff==null) || (polyCoeff.length==0)  || (polyCoeff[0]==null)) continue;

			nPix=0;
			int nOutlayers=0;
			double halfBinWidth=binWidth/2;
			int nGap=0; // should be 0 to count
			for (int i=0;i<input_bayer[chn].length;i++){
				double d= input_bayer[chn][i];
				if ((satLevel<=0) || (d<satLevel)){
					double x=i % tileHalfSize; // x
					double y=i / tileHalfSize; // y
					double f=0;
					if (polyCoeff[0].length>3){
						f=
								polyCoeff[0][0]*x*x+
								polyCoeff[0][1]*x*x+
								polyCoeff[0][2]*x*y+
								polyCoeff[0][3]*x+
								polyCoeff[0][4]*y+
								polyCoeff[0][5];
					} else {
						f=
								polyCoeff[0][0]*x+
								polyCoeff[0][1]*y+
								polyCoeff[0][2];
					}
					double diff=Math.abs(d-f);
					if (diff >= halfBinWidth) {
						if (diff >= (halfBinWidth+gapWidth)){
							nOutlayers++;
						} else {
							nGap++; // do not need to continue, only if for debug
							break;
						}
						
					}
					nPix++;
				}
			}
			if (nGap>0) continue; // gap is not empty, can not tell if the phase is wrong
			result[0]+=nOutlayers;
			result[1]+=nPix;
		}
		return result;
	}

	public double [] getDefectsBayer(
			final ImagePlus imp,
			final int tileClearSize0,
			final int tileMargins0,
			final int cmask, // bitmask of color channels to process (9 - two greens)
			final int numPasses, // number of passes to replace outlayers (will end if none was replaced)
			final int numInBase, // number of neighbors (of 8) to use as a base if they all agree
			final double sigma, // for high-pass filtering
			final boolean processNoReplace, // calculate differences even if no replacements were made
			final double binWidth, // absolute
			final double gapWidth, // absolute
			final int algorithmNumber,
			final int numInBase2, // number of neighbors (of 8) to use as a base if they all agree
			final double binWidth2, // absolute
			final double gapWidth2, // absolute
			final MatchSimulatedPattern matchSimulatedPattern,
			final EyesisAberrations.ColorComponents colorComponents,
			int threadsMax,
			final int debugLevel){
		double [][]defects=getDefects(
				imp,
				tileClearSize0,
				tileMargins0,
				cmask, // bitmask of color channels to process (9 - two greens)
				numPasses, // number of passes to replace outlayers (will end if none was replaced)
				numInBase, // number of neighbors (of 8) to use as a base if they all agree
				sigma, // for high-pass filtering
				processNoReplace, // calculate differences even if no replacements were made
				binWidth, // absolute
				gapWidth, // absolute
				algorithmNumber,
				numInBase2, // number of neighbors (of 8) to use as a base if they all agree
				binWidth2, // absolute
				gapWidth2, // absolute
				matchSimulatedPattern,
				colorComponents,
				threadsMax,
				debugLevel);
		
		int imgWidth=imp.getWidth();
		int imgHeight=imp.getHeight();
		int halfWidth=imgWidth/2;
		double [] pixels=new double [imgWidth*imgHeight];
		for (int i=0;i<pixels.length;i++){
			int y=i/imgWidth;
			int x=i%imgWidth;
			int chn=(y&1)*2+(x&1);
			int index=(y>>1)*halfWidth+(x>>1);
			pixels[i]=(defects[chn]==null)?0.0:defects[chn][index];
		}
		return pixels;
	}
	
	public double [][] getDefects(
			final ImagePlus imp,
			final int tileClearSize0,
			final int tileMargins0,
			final int cmask, // bitmask of color channels to process (9 - two greens)
			final int numPasses, // number of passes to replace outlayers (will end if none was replaced)
			final int numInBase, // number of neighbors (of 8) to use as a base if they all agree
			final double sigma, // for high-pass filtering
			final boolean processNoReplace, // calculate differences even if no replacements were made
			final double binWidth, // absolute
			final double gapWidth, // absolute
			final int algorithmNumber,
			final int numInBase2, // number of neighbors (of 8) to use as a base if they all agree
			final double binWidth2, // absolute
			final double gapWidth2, // absolute
			final MatchSimulatedPattern matchSimulatedPattern,
			final EyesisAberrations.ColorComponents colorComponents,
			int threadsMax,
			final int debugLevel){
		final int debugThreshold=0;
		//both tileMargins and tileClearSize should be even
		final int tileClearSize=2*(tileClearSize0/2);
		final int tileMargins=  2*(tileMargins0/2);
		final int imgWidth=imp.getWidth();
		final int imgHeight=imp.getHeight();
		Rectangle selection= (imp.getRoi()==null)? (new Rectangle(0, 0, imgWidth, imgHeight)):imp.getRoi().getBounds();
		final int x00=2*(selection.x/2);
		final int y00=2*(selection.y/2);
		final int numHor= selection.width/ tileClearSize+((tileClearSize*(selection.width/tileClearSize) < selection.width)?1:0);
		final int numVert=selection.height/tileClearSize+((tileClearSize*(selection.height/tileClearSize) < selection.height)?1:0);
		final double [][] diffs=new double [4][];
		for (int i=0;i<diffs.length;i++){
			if ((cmask& (1<<i))!=0){

				diffs[i]=new double [imgWidth*imgHeight/4];
				for (int j=0;j<diffs[i].length;j++) diffs[i][j]=0.0;
			} else {
				diffs[i]=null;
			}
		}
		if (debugLevel>debugThreshold){
			
		}
		final AtomicInteger aNumTile=new AtomicInteger(0);
		final Thread[] threads = newThreadArray(threadsMax);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			// Concurrently run in as many threads as CPUs
			threads[ithread] = new Thread() {
				public void run() {
					DoubleGaussianBlur doubleGaussianBlur=new DoubleGaussianBlur();
					int tileHalfSize=tileClearSize/2+tileMargins;
					int tileLow=tileMargins/2;
					int tileHigh=(tileMargins+tileClearSize)/2;
					int imgHalfWidth=imgWidth/2;
					int imgHalfHeight=imgHeight/2;
					for (int nTile = aNumTile.getAndIncrement(); nTile < (numHor*numVert); nTile = aNumTile.getAndIncrement()) {
						int nVert=nTile/numHor;
						int nHor=nTile%numHor;
						int x0=x00-tileMargins+(tileClearSize)*nHor;
						int y0=y00-tileMargins+(tileClearSize)*nVert;
						double [][] tileDiffs=getTileDefects( // per processed channel - pixel difference 
								x0, // top left in sensor pixels (even)
								y0,
								imp,
								tileHalfSize,
								cmask, // bitmask of color channels to process (9 - two greens)
								numPasses, // number of passes to replace outlayers (will end if none was replaced)
								numInBase, // number of neighbors (of 8) to use as a base if they all agree
								sigma, // for high-pass filtering
								processNoReplace, // calculate differences even if no replacements were made
								binWidth,
								gapWidth,
								algorithmNumber,
								numInBase2, // number of neighbors (of 8) to use as a base if they all agree (after HPF)
								binWidth2, // (after HPF)
								gapWidth2, // (after HPF)
								matchSimulatedPattern,
								colorComponents,
								doubleGaussianBlur, // or null
								debugLevel);
						// apply result tile
						for (int chn=0;chn<diffs.length;chn++) if ((diffs[chn]!=null) && (tileDiffs[chn]!=null)) {
//							if ((chn==0) &&(debugLevel>debugThreshold)){
//								System.out.println("x0="+x0+" y0="+y0+
//										" tileIndex="+tileLow*tileHalfSize+tileLow);
//							}
							for (int iy=tileLow;iy<tileHigh;iy++){
								int iyImg= y0/2 +iy;
								if ((iyImg>=0.0) && (iyImg<imgHalfHeight)) {
									for (int ix=tileLow;ix<tileHigh;ix++){
										int ixImg= x0/2 +ix;
										if ((ixImg>=0.0) && (ixImg<imgHalfWidth)) {
											int tileIndex=iy*tileHalfSize+ix;
											int destIndex=iyImg*imgHalfWidth+ixImg;
											diffs[chn][destIndex]=tileDiffs[chn][tileIndex];
										}
									}
								}
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return diffs;
	}

	
	
	private double [][] getTileDefects( // per processed channel - pixel difference 
			int x0, // top left in sensor pixels (even)
			int y0,
			ImagePlus imp,
			int tileHalfSize,
			int cmask, // bitmask of color channels to process (9 - two greens)
			int numPasses, // number of passes to replace outlayers (will end if none was replaced)
			int numInBase, // number of neighbors (of 8) to use as a base if they all agree
			double sigma, // for high-pass filtering
			boolean processNoReplace, // calculate differences even if no replacements were made
			double binWidth, // absolute
			double gapWidth, // absolute
			int algorithmNumber,
			int numInBase2, // number of neighbors (of 8) to use as a base if they all agree
			double binWidth2, // absolute
			double gapWidth2, // absolute
			MatchSimulatedPattern matchSimulatedPattern,
			EyesisAberrations.ColorComponents colorComponents,
			DoubleGaussianBlur doubleGaussianBlur, // or null
			int debugLevel){
		int debugThreshold=1;
//		int debugChn=0, debugX=18,debugY=-33;
		int [] dirs={1,tileHalfSize+1,tileHalfSize,tileHalfSize-1,-1,-tileHalfSize-1,-tileHalfSize,-tileHalfSize+1};
		Rectangle tileRect=new Rectangle (x0,y0,2*tileHalfSize,2*tileHalfSize);
		double [][] input_bayer=matchSimulatedPattern.splitBayer(imp,tileRect,colorComponents.equalizeGreens); // does it work the same?
		double [][] pixDiffs=new double [input_bayer.length][];
		for (int i=0;i<pixDiffs.length;i++) pixDiffs[i]=null;
		ArrayList<Integer> replaceList=new ArrayList<Integer>();
		for (int chn=0;chn<input_bayer.length;chn++) if ((input_bayer[chn]!=null) &&((cmask & (1<<chn))!=0)){
			double [] tilePix=input_bayer[chn].clone();
			double [] tilePixTmp=new double [tileHalfSize*tileHalfSize]; //tilePix.clone();
			int numReplacements=0;
			for (int nPass=0;nPass<numPasses;nPass++){ // normally will break from the loop
				replaceList.clear();
				for (int y=1;y<(tileHalfSize-1);y++) for (int x=1;x<(tileHalfSize-1);x++){
//					if ((debugLevel>debugThreshold) && (chn==debugChn) && (y==debugY) && (x==debugX)){
//						System.out.println("getTileDefects("+x0+","+y0+",...) chn="+chn+ " x=" +x+" y="+y);
//					}
					Integer index=y*tileHalfSize+x;
					double [] vals=new double[dirs.length];
					double d=tilePix[index]; // center
					for (int dir=0;dir<dirs.length;dir++) vals[dir]=tilePix[index+dirs[dir]];
					// sort neighbors - just bubble
					boolean outOfOrder;
					do {
						outOfOrder=false;
						for (int i=0;i<dirs.length-1;i++){
							if (vals[i]>vals[i+1]) {
								outOfOrder=true;
								double t=vals[i+1];
								vals[i+1]=vals[i];
								vals[i]=t;
							}
						}
					} while(outOfOrder);
					// find sequence long enough to fit in binWidth
					//int si=0, ei=numInBase;
					int si=0,ei=0;
					label_a:for (si=0;si<(dirs.length-numInBase);si++) for (ei=si+numInBase;ei<=dirs.length;ei++){
						if ((vals[ei-1]-vals[si])>binWidth) continue; // bin too wide
						if ((si>0) && ((vals[si]-vals[si-1])<gapWidth)) continue; // lower too close
						if ((ei<dirs.length) && ((vals[ei]-vals[ei-1])<gapWidth)) continue; // higher too close
						break label_a; // found good
					}
					if (si>=(dirs.length-numInBase)) continue;
					if ((d>=vals[si]) && (d<(vals[ei-1] +gapWidth))) continue; // center pixel not high enough to be replaced
					if ((d<=vals[ei-1]) && (d>(vals[si] -gapWidth))) continue; // center pixel not low enough to be replaced
					double s=0;
					for (int i=si;i<ei;i++) s+=vals[i];
					tilePixTmp[index]=s/(ei-si); // average
					replaceList.add(index);
				}
				if (replaceList.isEmpty()) break;
				numReplacements+=replaceList.size();
				for (Integer index:replaceList) tilePix[index]=tilePixTmp[index]; 
				if (debugLevel>debugThreshold){
					System.out.println("getTileDefects("+x0+"/"+y0+"): pass="+nPass+" replaced: "+ replaceList.size());
				}
			}
			if ((numReplacements==0) && !processNoReplace) continue;
// run low pass filter on the tile with outlayers replaced with the averages 
			if (doubleGaussianBlur==null) doubleGaussianBlur=new DoubleGaussianBlur();
			doubleGaussianBlur.blurDouble(tilePix,tileHalfSize,tileHalfSize,sigma,sigma, 0.01);
			for (int i=0;i<tilePix.length;i++) tilePix[i]=input_bayer[chn][i]-tilePix[i];
			if (numInBase2>0){ // 0 - skip final defects filtering, leave all HPF differences
				replaceList.clear();
				for (int y=1;y<(tileHalfSize-1);y++) for (int x=1;x<(tileHalfSize-1);x++){
					//					if ((debugLevel>debugThreshold) && (chn==debugChn) && (y==debugY) && (x==debugX)){
					//						System.out.println("getTileDefects("+x0+","+y0+",...) chn="+chn+ " x=" +x+" y="+y);
					//					}
					Integer index=y*tileHalfSize+x;
					double [] vals=new double[dirs.length];
					double d=tilePix[index]; // center
					for (int dir=0;dir<dirs.length;dir++) vals[dir]=tilePix[index+dirs[dir]];
					// sort neighbors - just bubble
					switch (algorithmNumber){
					case 1:
						boolean outOfOrder;
						do {
							outOfOrder=false;
							for (int i=0;i<dirs.length-1;i++){
								if (vals[i]>vals[i+1]) {
									outOfOrder=true;
									double t=vals[i+1];
									vals[i+1]=vals[i];
									vals[i]=t;
								}
							}
						} while(outOfOrder);
						// find sequence long enough to fit in binWidth
						//int si=0, ei=numInBase;
						int si=0,ei=0;
						label_a:for (si=0;si<(dirs.length-numInBase2);si++) for (ei=si+numInBase2;ei<=dirs.length;ei++){
							if ((vals[ei-1]-vals[si])>binWidth2) continue; // bin too wide
							if ((si>0) && ((vals[si]-vals[si-1])<gapWidth2)) continue; // lower too close
							if ((ei<dirs.length) && ((vals[ei]-vals[ei-1])<gapWidth2)) continue; // higher too close
							break label_a; // found good
						}
						if (si>=(dirs.length-numInBase2)) continue;
						if ((d>=vals[si]) && (d<(vals[ei-1] +gapWidth2))) continue; // center pixel not high enough to be replaced
						if ((d<=vals[ei-1]) && (d>(vals[si] -gapWidth2))) continue; // center pixel not low enough to be replaced
						replaceList.add(index);
						break;
					case 2:
						double gapDistance=binWidth2/2+gapWidth2;
						if (Math.abs(d)>=gapDistance) {
							int numInBin=0;
							double binDistance=binWidth2/2;
							for (double v:vals){
								double a=Math.abs(v);
								if (a<=binDistance) numInBin++;
								else if (a<gapDistance){
									numInBin=-1;
									break;
								}
							}
							if (numInBin>=numInBase2) {
								replaceList.add(index);
							}
						}
						break;
					}
				}
				if (replaceList.isEmpty()) continue; // no pixel defects in this tile
				for (int i=0;i<tilePixTmp.length;i++) tilePixTmp[i]=0.0;
				for (Integer index:replaceList){
					tilePixTmp[index]=tilePix[index];
				}
				pixDiffs[chn]=tilePixTmp;
				
			} else {
				pixDiffs[chn]=tilePix;
			}
		}
		return pixDiffs;
	}
	
	public class SensorDefects{
		public int [] defectRepetitions;
		public int [] defectPositiveRepetitions;
		public int numberOfImages;
		public int numberOfDefectiveImages;
		public int width;
		public float [] pixels;
		
		public SensorDefects(ImagePlus imp){
			defectRepetitions=null;
			defectPositiveRepetitions=null;
			numberOfImages=0;
			numberOfDefectiveImages=0;
			accummulateImage(imp);
		}
		public SensorDefects(String path){
			defectRepetitions=null;
			defectPositiveRepetitions=null;
			numberOfImages=0;
			numberOfDefectiveImages=0;
			pixels=null;
			accummulateImage(path);
		}

		public void accummulateImage(String path){
			ImagePlus imp=null;
			try {
				imp=new ImagePlus(path);
			} catch (Exception e){
				System.out.println("Failed to open image file: "+path);
				return;
			}
			accummulateImage(imp);
		}

		public void accummulateImage(ImagePlus imp){
			ImageProcessor  ip=imp.getProcessor();
			if (pixels==null){
				pixels=((float[])ip.getPixels()).clone(); 
				width= imp.getWidth();  // full image width
			} else {
				float [] imagePixels=(float[])ip.getPixels();
				for (int i=0;i<pixels.length;i++){
					pixels[i]+=imagePixels[i];
				}
			}
			numberOfImages++;
		}
		
		public ImagePlus getAccummulatedImage(String title){
			if (pixels==null) return null;
			int height=pixels.length/width;
			ImageProcessor ip=new FloatProcessor(width,height);
			ip.setPixels(getPixels());
			ip.resetMinAndMax();
			ImagePlus imp=  new ImagePlus(title, ip);
			return imp;
		}
		
		
		public float [] getPixels(){
			if (pixels==null) return null;
			float [] averagePixels=this.pixels.clone();
			for (int i=0;i<averagePixels.length;i++) averagePixels[i]/=numberOfImages;
			return averagePixels;
		}
		
		public SensorDefects(){
			pixels=null;
			defectRepetitions=null;
			defectPositiveRepetitions=null;
			numberOfImages=0;
			numberOfDefectiveImages=0;
			width=0;
		}
		public SensorDefects(int width,int height){
			pixels=null;
			defectRepetitions=new int[width*height];
			for (int i=0;i<defectRepetitions.length;i++) defectRepetitions[i]=0;
			defectPositiveRepetitions=new int[width*height];
			for (int i=0;i<defectPositiveRepetitions.length;i++) defectPositiveRepetitions[i]=0;
			this.width=width;
		}
		public void incDefect(int index, boolean positive){
			defectRepetitions[index]++;
			if (positive) defectPositiveRepetitions[index]++;
		}
		public void incNumberOfImages(boolean defective){
			numberOfImages++;
			if (defective) numberOfDefectiveImages++;
		}
		public int [] getDefectRepetiotions(){
			return this.defectRepetitions;
		}
		public int [] getDefectPositiveRepetitions(){
			return this.defectPositiveRepetitions;
		}
		public int getNumberOfImages(){
			return numberOfImages;
		}
		public int getNumberOfDefectiveImages(){
			return numberOfDefectiveImages;
		}
		public int [] getWidthHeight(){
			int [] result={this.width, this.defectRepetitions.length/this.width};
			return result;
		}
	}
	public String getChannel(ImagePlus imp){
		// read image info to properties (if it was not done yet - should it?
		if ((imp.getProperty("timestamp")==null) || (((String) imp.getProperty("timestamp")).length()==0)) {
			JP46_Reader_camera jp4_instance= new JP46_Reader_camera(false);
			jp4_instance.decodeProperiesFromInfo(imp);
		}
		if (imp.getProperty("channel")==null) return null;
		return (String) imp.getProperty("channel");
	}
	
	public boolean interactiveDefectivePixelList(
			SensorDefects[] defects,
			boolean enableHot,
			boolean enableCold,
			boolean enableMixed,
			int minConfirmations,
			int debugLevel
			){
		GenericDialog gd = new GenericDialog("Sensor defective pixel list generation");
		gd.addCheckbox    ("Include hot pixels", enableHot);
		gd.addCheckbox    ("Include cold (dark) pixels", enableCold);
		gd.addCheckbox    ("Include mixed (above/below neighbors) pixels", enableMixed);
		gd.addNumericField("Minimal number of images where defect should be present", minConfirmations, 0 ,4, "images");
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		enableHot=      gd.getNextBoolean();
		enableCold=      gd.getNextBoolean();
		enableMixed=      gd.getNextBoolean();
		minConfirmations=    (int) gd.getNextNumber();
		int [][][] defectList=getDefectivePixelList(
				defects,
				enableHot,
				enableCold,
				enableMixed,
				minConfirmations,
				debugLevel);
		System.out.println("Defective pixels lists report: minConfirmations="+minConfirmations+
				" enableHot="+enableHot+" enableCold="+enableCold+" enableMixed="+enableMixed);
		System.out.println("========================================================================");
		int numInLine=8;
		for (int chn=0;chn<defectList.length;chn++) if ((defectList[chn]!=null) && (defectList[chn].length>0)){
			System.out.print("SFE #"+chn+": ");
			for (int i=0;i<defectList[chn].length;i++){
				System.out.print(defectList[chn][i][0]+":"+defectList[chn][i][1]+" ");
				if (((i%numInLine)==(numInLine-1)) || (i == (defectList[chn].length-1))) System.out.println();
			}
			System.out.println();
		}
		return true;
	}
	
	
	public int [][][] getDefectivePixelList(
			SensorDefects[] defects,
			boolean enableHot,
			boolean enableCold,
			boolean enableMixed,
			int minConfirmations,
			int debugLevel
			){
		if (defects==null) return null;
		int [][][] result = new int [defects.length][][];
		for (int chn=0;chn<defects.length;chn++){
			if (defects[chn]!=null) {
				ArrayList<Integer> xyList=new ArrayList<Integer>();
				int [] mapRepeat=       defects[chn].getDefectRepetiotions();
				int [] mapPosRepeat=    defects[chn].getDefectPositiveRepetitions();
				for (int i=0;i<mapRepeat.length;i++)if (mapRepeat[i]>=minConfirmations){
					if (!enableHot   && (mapPosRepeat[i]==mapRepeat[i])) continue;
					if (!enableCold  && (mapPosRepeat[i]==0)) continue;
					if (!enableMixed && (mapPosRepeat[i]!=0) && (mapPosRepeat[i]!=mapRepeat[i])) continue;
					xyList.add(new Integer(i));
				}
				if (!xyList.isEmpty()){
					int width=defects[chn].getWidthHeight()[0];
					result[chn]=new int[xyList.size()][2];
					int index=0;
					for (Integer xy:xyList){
						result[chn][index][0]=  xy%width;
						result[chn][index++][1]=xy/width;
					}
				}
			} else {
				result[chn]=null;
			}
		}
		return result;
	}
	
	// Ask for image selection (multiple directories), accumulate and create array of accumulated images
	public ImagePlus[] getInteractiveAccumulatedImages(
			DistortionProcessConfiguration distortionProcessConfiguration,
			final AtomicInteger stopRequested,
			int threadsMax,
			boolean updateStatus,
			int debugLevel){
		ArrayList <String> sourceList=new ArrayList<String>();
		while (true) {
			String [] files=distortionProcessConfiguration.selectSourceFiles(); // select files - with/without dialog
			if ((files==null) || (files.length==0)) break;
			for (String path:files) sourceList.add(path);
		}
		String [] sourceFilesList=sourceList.toArray(new String[0]);
		SensorDefects[] accumulatedImages= accumulateImages(
				sourceFilesList,
				stopRequested,
				threadsMax,
				updateStatus,
				debugLevel);
		ImagePlus [] imps=new ImagePlus[accumulatedImages.length];
		for (int chn=0;chn<imps.length;chn++){
			if (accumulatedImages[chn]!=null){
				imps[chn]=accumulatedImages[chn].getAccummulatedImage("Accummulated-"+chn);
			} else {
				imps[chn]=null;
			}
		}
		return imps;
	}
	
	

	public Defect [][] interactiveExtractDefectListsFromAccumulatedImages(
			ImagePlus [] imps,
			int tileClearSize,
			int tileMargins,
			int cmask, // bitmask of color channels to process (9 - two greens)
			int numPasses, // number of passes to replace outlayers (will end if none was replaced)
			int numInBase, // number of neighbors (of 8) to use as a base if they all agree
			double sigma, // for high-pass filtering
			boolean processNoReplace, // calculate differences even if no replacements were made
			double binWidth, // absolute
			double gapWidth, // absolute
			int algorithmNumber,
			int numInBase2, // number of neighbors (of 8) to use as a base if they all agree
			double binWidth2, // absolute
			double gapWidth2, // absolute
			boolean processHot,
			boolean processCold,
			boolean updateSensorCalibrationFiles,
			boolean clearDefects, // clear defects if none detected
			MatchSimulatedPattern matchSimulatedPattern,
			EyesisAberrations.ColorComponents colorComponents,
			int threadsMax,
			int debugLevel){
		GenericDialog gd = new GenericDialog("Sensor defects lists generation");
		gd.addNumericField("Tile clearSize", tileClearSize, 0 ,4, "pix");
		gd.addNumericField("Tile margin (extra) width", tileMargins, 0 ,4, "pix");
		gd.addNumericField("Color channel mask", cmask, 0 ,4, "9 - both green colors");
		gd.addNumericField("Number of \"remove outlayers\" passes", numPasses, 0 ,4, "");
		gd.addNumericField("Number of neighbors (of total 8) to base outlayers", numInBase, 0 ,4, "<=8");
		gd.addNumericField("Base low-pass sigma ", sigma, 3 ,7, "double pixels");
		gd.addCheckbox    ("Process tiles with no outlayer replacements", processNoReplace);
		gd.addNumericField("\"Normal\" pixels variation", binWidth, 3 ,7, "");
		gd.addNumericField("Empty level gap between \"normal\" pixels and outlayers", gapWidth, 3 ,7, "");

		gd.addNumericField("Post-HPF algorithm number", algorithmNumber, 0 ,4, "(1 or 2");
		gd.addNumericField("Number of neighbors (of total 8) to base outlayers (after high-pass), 0 - no filter", numInBase2, 0 ,4, "<=8");
		gd.addNumericField("\"Normal\" pixels variation (after high-pass)", binWidth2, 3 ,7, "");
		gd.addNumericField("Empty level gap between \"normal\" pixels and outlayers (after high-pass)", gapWidth2, 3 ,7, "");
		gd.addCheckbox    ("Update sensor calibration files with defects data", updateSensorCalibrationFiles);
		gd.addCheckbox    ("Clear defects in sensor calibration files if none new detected", clearDefects);

		gd.showDialog();
		if (gd.wasCanceled()) return null;
		tileClearSize=    (int) gd.getNextNumber();
		tileMargins=    (int) gd.getNextNumber();
		cmask=           (int) gd.getNextNumber();
		numPasses=       (int) gd.getNextNumber();
		numInBase=       (int) gd.getNextNumber();
		sigma=                 gd.getNextNumber();
		processNoReplace=      gd.getNextBoolean();
		binWidth =             gd.getNextNumber();
		gapWidth =             gd.getNextNumber();

		algorithmNumber=          (int) gd.getNextNumber();
		numInBase2=      (int) gd.getNextNumber();
		binWidth2 =            gd.getNextNumber();
		gapWidth2 =            gd.getNextNumber();
		updateSensorCalibrationFiles=      gd.getNextBoolean();
		clearDefects=          gd.getNextBoolean();
		Defect [][] defects= extractDefectListsFromAccumulatedImages(
				imps,
				tileClearSize,
				tileMargins,
				cmask, // bitmask of color channels to process (9 - two greens)
				numPasses, // number of passes to replace outlayers (will end if none was replaced)
				numInBase, // number of neighbors (of 8) to use as a base if they all agree
				sigma, // for high-pass filtering
				processNoReplace, // calculate differences even if no replacements were made
				binWidth, // absolute
				gapWidth, // absolute
				algorithmNumber,
				numInBase2, // number of neighbors (of 8) to use as a base if they all agree
				binWidth2, // absolute
				gapWidth2, // absolute
				processHot,
				processCold,
				matchSimulatedPattern,
				colorComponents,
				threadsMax,
				debugLevel);
		if (updateSensorCalibrationFiles) {
			String [] extensions={".calib-tiff"};
			CalibrationFileManagement.MultipleExtensionsFileFilter parFilter = new CalibrationFileManagement.MultipleExtensionsFileFilter("",extensions,"distortion calibration .calib-tiff files");
			String anySensorPath=CalibrationFileManagement.selectFile(false,
					"Select one of the sensor calibration files of a set to be updated",
					"Update",
					parFilter,
					""); //String defaultPath
			if ((anySensorPath==null) || (anySensorPath=="")) return defects;
			int numUpdated=applyDefectsMapToSensorFile( // modify to erase sesnor calibration data
					 anySensorPath,
					 clearDefects,
					 defects,
					 debugLevel);
			if (debugLevel>0){
				System.out.println("Updated "+numUpdated+" sensor calibration files");
			}
		}
		return defects;
	}
	
	public Defect [][] extractDefectListsFromAccumulatedImages(
					ImagePlus [] imps,
					int tileClearSize,
					int tileMargins,
					int cmask, // bitmask of color channels to process (9 - two greens)
					int numPasses, // number of passes to replace outlayers (will end if none was replaced)
					int numInBase, // number of neighbors (of 8) to use as a base if they all agree
					double sigma, // for high-pass filtering
					boolean processNoReplace, // calculate differences even if no replacements were made
					double binWidth, // absolute
					double gapWidth, // absolute
					int algorithmNumber,
					int numInBase2, // number of neighbors (of 8) to use as a base if they all agree
					double binWidth2, // absolute
					double gapWidth2, // absolute
					boolean processHot,
					boolean processCold,
					MatchSimulatedPattern matchSimulatedPattern,
					EyesisAberrations.ColorComponents colorComponents,
					int threadsMax,
					int debugLevel){
		Defect [][] resultList = new Defect [imps.length][];
		for (int chn=0;chn<resultList.length;chn++){
			if (imps[chn]==null){
				resultList[chn]=null;
			} else {
				int width=imps[chn].getWidth();
				double [] defects= getDefectsBayer(
						imps[chn],
						tileClearSize,
						tileMargins,
						cmask, // bitmask of color channels to process (9 - two greens)
						numPasses, // number of passes to replace outlayers (will end if none was replaced)
						numInBase, // number of neighbors (of 8) to use as a base if they all agree
						sigma, // for high-pass filtering
						processNoReplace, // calculate differences even if no replacements were made
						binWidth, // absolute
						gapWidth, // absolute
						algorithmNumber,
						numInBase2, // number of neighbors (of 8) to use as a base if they all agree
						binWidth2, // absolute
						gapWidth2, // absolute
						matchSimulatedPattern,
						colorComponents,
						threadsMax,
						debugLevel);
				int numDefects=0;
				for (int i=0;i<defects.length;i++) if (defects[i]!=0) numDefects++;
				if (numDefects==0) {
					resultList[chn]=null;
				} else {
					ArrayList<Defect> dList=new ArrayList<Defect>();
					for (int i=0;i<defects.length;i++) if ((defects[i]!=0) && ((defects[i]>0)?processHot:processCold)) {
						dList.add(new Defect(i%width,i/width,defects[i]));
					}
					// convert index list to array of x/y pairs
					resultList[chn]= dList.toArray(new Defect[0]);
					// sort array of defects
					Arrays.sort(resultList[chn], new Comparator<Defect>() {
						@Override
						public int compare(Defect o1, Defect o2) {
							return ((Double) Math.abs(o2.diff)).compareTo(Math.abs(o1.diff)); 
						}
					});
				}
			}
		}
		return resultList;
	}

	public int applyDefectsMapToSensorFile(
			 String anySensorPath,
			 boolean clearDefects,
			 Defect [][] defects,
			 int debugLevel){
		int debugThreshold=0;
		int numFilesUpdated=0;
    	int indexPeriod=anySensorPath.indexOf('.',anySensorPath.lastIndexOf(Prefs.getFileSeparator()));
    	JP46_Reader_camera jp4=new JP46_Reader_camera(false);

    	for (int chn=0;chn<defects.length;chn++) if ((defects[chn]!=null) || clearDefects){
    		String channelPath=anySensorPath.substring(0,indexPeriod-2)+String.format("%02d",chn)+anySensorPath.substring(indexPeriod);
    		if (debugLevel>debugThreshold) System.out.println("applyDefectsMapToSensorFile(): reading "+channelPath);
    		Opener opener=new Opener();
    		ImagePlus imp=opener.openImage("", channelPath);
        	if (imp==null) {
        		String msg="Failed to read sensor calibration data file "+channelPath;
//        		IJ.showMessage("Error",msg);
        		System.out.println(msg);
        		continue; // next channel
        	}
        	jp4.decodeProperiesFromInfo(imp);
        	// Set/replace Distortions
        	imp.setProperty("comment_defects", "Sensor hot/cold pixels list as x:y:difference");
        	if (defects[chn]!=null) {
        		StringBuffer sb = new StringBuffer();
        		for (Defect d:defects[chn]){
        			if (sb.length()>0) sb.append(" ");
        			sb.append(d.x+":"+d.y+":"+IJ.d2s(d.diff,3));
        		}
        		imp.setProperty("defects", sb.toString());
        	} else {
        		imp.setProperty("defects", null);
        	}
        	jp4.encodeProperiesToInfo(imp);
        	// Save modified file back
        	FileSaver fs=new FileSaver(imp);
        	String msg="Saving sensor calibration with defects list to "+channelPath;
        	if (debugLevel>debugThreshold) System.out.println(msg);
        	try {
            	fs.saveAsTiffStack(channelPath);
        	} catch (Exception e){
        		String msg1="Failed to write sensor calibration file "+channelPath;
//        		IJ.showMessage("Error",msg1);
        		System.out.println(msg1);
        		continue; // next channel
        	}
        	numFilesUpdated++;
    	}
		
		return numFilesUpdated;
	}
	
/*
    	FileSaver fs=new FileSaver(imp);
    	String msg="Saving "+(realData?"":"EMPTY")+" sensor distortions to "+path;
    	if (updateStatus) IJ.showStatus(msg);
    	if (this.debugLevel>0) System.out.println(msg);
    	fs.saveAsTiffStack(path);
    	if (this.pathNames==null){
    		this.pathNames=new String[this.fittingStrategy.distortionCalibrationData.getNumChannels()];
    		for (int i=0;i<this.pathNames.length;i++) this.pathNames[i]=null; 
    	}
    	this.pathNames[numSensor]=path;
    	return imp;
 * 
		Opener opener=new Opener();
		ImagePlus imp=opener.openImage("", path);
    	if (imp==null) {
    		if (!reportProblems) return;
    		String msg="Failed to read sensor calibration data file "+path;
    		IJ.showMessage("Error",msg);
    		System.out.println(msg);
    		throw new IllegalArgumentException (msg);
    	}
		if (this.debugLevel>0) System.out.println("Read "+path+" as a sensor calibration data");
    	(new JP46_Reader_camera(false)).decodeProperiesFromInfo(imp);
    	setDistortionFromImageStack(imp, numSensor, overwriteExtrinsic);
 * 
 D5917: public void setDistortionFromImageStack(String path, boolean overwriteExtrinsic){
    	int indexPeriod=path.indexOf('.',path.lastIndexOf(Prefs.getFileSeparator()));
    	int numSubCameras=fittingStrategy.distortionCalibrationData.eyesisCameraParameters.eyesisSubCameras[0].length;
    	for (int chNum=0;chNum<numSubCameras;chNum++){
    		String channelPath=path.substring(0,indexPeriod-2)+String.format("%02d",chNum)+path.substring(indexPeriod);
    		try { // disable here for now
    			setDistortionFromImageStack(channelPath, chNum, false, overwriteExtrinsic);
    		} catch (Exception e) {
    			System.out.println("setDistortionFromImageStack(): " + e.toString());
    			e.printStackTrace();
    		}
    	}
    }
	
 */
	
	
	
	public SensorDefects[] accumulateImages(
			String [] sourceFilesList,
			final AtomicInteger stopRequested,
			int threadsMax,
			final boolean updateStatus,
			final int debugLevel){
		final long 	  startTime=System.nanoTime();
		final int debugThreshold=1;
		//  sort images into channels
		final Map<Integer,ArrayList<String>> sensorFiles=new HashMap<Integer,ArrayList<String>>();
		final Map<Integer,SensorDefects> sensorDefectsMap=new HashMap<Integer,SensorDefects>();
		for (int numFile=0;numFile<sourceFilesList.length;numFile++){
			Integer nChn=-1;
			try {
				int indexDot=sourceFilesList[numFile].lastIndexOf(".");
				int indexDash=sourceFilesList[numFile].lastIndexOf("-",indexDot);
				// extract channel from path without opening actual image
				nChn=Integer.parseInt(sourceFilesList[numFile].substring(indexDash+1, indexDot));
			} catch (Exception e){
				System.out.println("Failed to extract channel number from path "+sourceFilesList[numFile]);
				continue;
			}
			
			ArrayList<String> fileList=sensorFiles.get(nChn);
			if (fileList==null){
				fileList=new ArrayList<String>();
				sensorFiles.put(nChn,fileList);
				SensorDefects sd=new SensorDefects();
				sensorDefectsMap.put(nChn, sd);
			}
			fileList.add(sourceFilesList[numFile]);
		}
		final Integer [] channels= sensorFiles.keySet().toArray(new Integer[0]);
		final AtomicInteger aIngex=new AtomicInteger(0);
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger processedImages=new AtomicInteger(0);
		final int numImages=sourceFilesList.length;
		final int progressStep=10;
		if (updateStatus) IJ.showStatus("Averaging images ...");
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int index = aIngex.getAndIncrement(); index < channels.length; index = aIngex.getAndIncrement()) {
						ArrayList<String> fileList=sensorFiles.get(channels[index]);
						SensorDefects sd=sensorDefectsMap.get(channels[index]);
						for (String path:fileList) {
							if (stopRequested.get()>0) {
								System.out.println("User requested stop");
								break;
							}
							sd.accummulateImage(path);
							if (updateStatus) {
								final int numFinished=processedImages.getAndIncrement();
								if (numFinished % progressStep == 0) {
									SwingUtilities.invokeLater(new Runnable() {
										public void run() {
											IJ.showStatus("Averaging images "+IJ.d2s(100.0*numFinished/numImages,2)+"%");
//											IJ.showProgress(numFinished,numImages); // Does not work - overwritten by imageJ image load
										}
									});
								}
							}
							if (debugLevel>debugThreshold){
								System.out.println(IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+"s: finished accumulating file "+path);
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
		if (updateStatus) IJ.showStatus("Averaging images ... done");
		if (updateStatus) IJ.showProgress(0);
		int maxChannel=-1;
		for (Integer chn:channels) if (maxChannel<chn) maxChannel=chn;
		SensorDefects[] result=new SensorDefects[maxChannel+1];
		for (int i=0;i<result.length;i++) result[i]=sensorDefectsMap.get(i); // may be null
		System.out.println(IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+"s: finished accumulating files");
		return result;
	}
	
	public SensorDefects[] accummulateSensorDefects(
			DistortionProcessConfiguration distortionProcessConfiguration,
			int tileClearSize,
			int tileMargins,
			int cmask, // bitmask of color channels to process (9 - two greens)
			int numPasses, // number of passes to replace outlayers (will end if none was replaced)
			int numInBase, // number of neighbors (of 8) to use as a base if they all agree
			double sigma, // for high-pass filtering
			boolean processNoReplace, // calculate differences even if no replacements were made
			double binWidth, // absolute
			double gapWidth, // absolute
			int algorithmNumber,
			int numInBase2, // number of neighbors (of 8) to use as a base if they all agree
			double binWidth2, // absolute
			double gapWidth2, // absolute
			MatchSimulatedPattern matchSimulatedPattern,
			EyesisAberrations.ColorComponents colorComponents,
			AtomicInteger stopRequested,
			int threadsMax,
			int debugLevel){

		Map<Integer,SensorDefects> sensorMap=new HashMap<Integer,SensorDefects>();
		String [] sourceFilesList=distortionProcessConfiguration.selectSourceFiles(); // select files - with/without dialog
		long 	  startTime=System.nanoTime();
		ImagePlus imp_sel;
		for (int numFile=0;numFile<sourceFilesList.length;numFile++){
			long 	  startFileTime=System.nanoTime();
			if (debugLevel>0){
				System.out.println(IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+"s: Processing file # "+(numFile+1)+ " (of "+ sourceFilesList.length+"): "+sourceFilesList[numFile]);
			}
			imp_sel=new ImagePlus(sourceFilesList[numFile]); // read source file
			Integer chn =Integer.parseInt(getChannel(imp_sel));
			SensorDefects sensorDefects=sensorMap.get(chn);
			if (sensorDefects==null){
				sensorDefects=new SensorDefects(imp_sel.getWidth(),imp_sel.getHeight());
				sensorMap.put(chn, sensorDefects);
			}
			double [] defects=getDefectsBayer(
					imp_sel,
					tileClearSize,
					tileMargins,
					cmask, // bitmask of color channels to process (9 - two greens)
					numPasses, // number of passes to replace outlayers (will end if none was replaced)
					numInBase, // number of neighbors (of 8) to use as a base if they all agree
					sigma, // for high-pass filtering
					processNoReplace, // calculate differences even if no replacements were made
					binWidth, // absolute
					gapWidth, // absolute
					algorithmNumber,
				    numInBase2,
				    binWidth2,
				    gapWidth2,
					matchSimulatedPattern, 
					colorComponents,
					threadsMax,
					debugLevel);
			int numDefects=0;
			for (int i=0;i<defects.length;i++) if (defects[i]!=0){
				sensorDefects.incDefect(i, defects[i]>0);
				numDefects++;
			}
			sensorDefects.incNumberOfImages(numDefects>0);
			if (debugLevel>0) System.out.println("File "+(numFile+1)+" calculation done at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+" (in "+
					IJ.d2s(0.000000001*(System.nanoTime()-startFileTime),3)+"s ), added "+numDefects+" defects");

			//				
			if (stopRequested.get()>0) {
				System.out.println("User requested stop");
				break;
			}
		}
		int maxChannel=0;
		for (Integer chn:sensorMap.keySet()){
			if (maxChannel<chn)maxChannel = chn;
		}
		SensorDefects[] intSensorMap=new SensorDefects [maxChannel+1];
		for (Integer chn=0;chn<intSensorMap.length;chn++){
			SensorDefects sensorDefects=sensorMap.get(chn);
			if (sensorDefects!=null){
				intSensorMap[chn]=sensorDefects;
			} else {
				intSensorMap[chn]=null;
			}

		}
		if (debugLevel>0) {
			System.out.print("Defects report format:  occurrences_in_images*number ( of them Hot, Cold, Mixed))");

			for (int chn=0;chn<intSensorMap.length;chn++) if (intSensorMap[chn]!=null){
				int maxRepeat=0;
				int [] mapRepeat=       intSensorMap[chn].getDefectRepetiotions();
				int [] mapPosRepeat=    intSensorMap[chn].getDefectPositiveRepetitions();
				int numImages=          intSensorMap[chn].getNumberOfImages();
				int numDefectiveImages= intSensorMap[chn].getNumberOfDefectiveImages();
				
				for (int i=0;i<mapRepeat.length;i++) if (maxRepeat<mapRepeat[i]) maxRepeat = mapRepeat[i];
				System.out.print("Sensor #"+chn+" ");
				if (maxRepeat==0){
					System.out.println("has no defects detected");
				} else {
					int [] histogram=    new int[maxRepeat]; // all - positive, negative, mixed
					int [] histogramHot= new int[maxRepeat]; // positive only
					int [] histogramCold=new int[maxRepeat]; // negative only
					for (int i=0;i<mapRepeat.length;i++) if (mapRepeat[i]>0) {
						histogram[mapRepeat[i]-1]++;
						if (mapPosRepeat[i]==mapRepeat[i]) histogramHot[mapRepeat[i]-1]++;
						if (mapPosRepeat[i]==0)            histogramCold[mapRepeat[i]-1]++;
					}
					double fracDefective=((double) numDefectiveImages)/numImages;
					System.out.print("has "+IJ.d2s(100*fracDefective,1)+"% defective images ("+numDefectiveImages+"/"+numImages+"): ");
					for (int i=0;i<histogram.length;i++) if (histogram[i]!=0){
						System.out.print((i+1)+"*"+histogram[i]+" ( ");
						if (histogramHot[i]>0)  System.out.print("H:"+histogramHot[i]+" ");
						if (histogramCold[i]>0) System.out.print("C:"+histogramCold[i]+" ");
						int m=histogram[i]-histogramHot[i]-histogramCold[i];
						if (m>0) System.out.print("M:"+m+" ");
						System.out.print(") ");
					}
					System.out.println();
				}
			}
		}
		return intSensorMap;
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

/**
 * Approximate function z(x,y) as a second degree polynomial (or just linear)
 * f(x,y)=A*x^2+B*y^2+C*x*y+D*x+E*y+F or f(x,y)=D*x+E*y+F 
 * data array consists of lines of either 2 or 3 vectors:
 *  2-element vector x,y
 *  variable length vector z (should be the same for all samples)
 *  optional 1- element vector w (weight of the sample)
 * 
 * returns array of vectors or null
 * each vector (one per each z component) is either 6-element-  (A,B,C,D,E,F) if quadratic is possible and enabled
 * or 3-element - (D,E,F) if linear is possible and quadratic is not possible or disabled
 * returns null if not enough data even for the linear approximation
 
 */
/*
   public double [][] quadraticApproximation(
		   double [][][] data,
		   boolean forceLinear  // use linear approximation
		   ){
*/


/*
 * 		final Thread[] threads = newThreadArray(threadsMax);

			imp_sel = WindowManager.getCurrentImage();
			if (imp_sel==null){
				IJ.showMessage("Error","There is no image selected");
				return;
			}
			ImagePlus imp_flow=SDFA_INSTANCE.showFlowFromSlices(imp_sel);
			imp_flow.getProcessor().resetMinAndMax();
			imp_flow.show();


*/