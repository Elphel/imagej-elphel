/**
 **
 ** SuperTiles - Process tiles resulted from quad image sets
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  SuperTiles.java is free software: you can redistribute it and/or modify
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
import java.util.concurrent.atomic.AtomicInteger;

//import TilePlanes.PlaneData;

//import TilePlanes.PlaneData;

public class SuperTiles{
	        TileProcessor tileProcessor;
			double      step_far;
			double      step_near;
			double      step_threshold_far;    // relative to min_disparity
			double      step_threshold_near;   // relative to min_disparity
			double      bin_far;               // bin value (before rounding) matching  (min_disparity + step_threshold_far)
			double      bin_near;              // bin value (before rounding) matching  (min_disparity + step_threshold_near)
			double      min_disparity;
			double      max_disparity;
			double      strength_floor;
			double      strength_pow;
			double      stBlurSigma;
			int         numBins;
			double []   bin_centers;
			double [][] disparityHistograms = null;
			double []   stStrength =          null; // per super-tile correlation strength
			double [][][] maxMinMax = null;
			double []   bgDisparity = null;
			double []   bgStrength = null;
//			TileProcessor.CLTPass3d cltPass3d;
			CLTPass3d cltPass3d;
//			ArrayList<TilePlanes.PlaneData>[] planes = null;
			TilePlanes.PlaneData [][] planes = null;
			public SuperTiles(
					CLTPass3d cltPass3d,
					double                  step_near,
					double                  step_far,
					double                  step_threshold,
					double                  min_disparity,
					double                  max_disparity,
					double                  strength_floor,
					double                  strength_pow,
					double                  stBlurSigma)
			{
				this.cltPass3d =           cltPass3d;
				this.tileProcessor =       cltPass3d.getTileProcessor();
				this.step_near =           step_near;
				this.step_far =            step_far;
				this.step_threshold_far = step_threshold;
				this.min_disparity =  min_disparity;
				this.max_disparity =  max_disparity;
				this.strength_floor = strength_floor;
				this.strength_pow =   strength_pow;
				this.stBlurSigma =    stBlurSigma;
//				this.step_threshold_near = this.step_threshold_far * this.step_far / step_near;
				this.step_threshold_near = this.step_threshold_far * step_near / this.step_far ;
				this.bin_far =             this.step_threshold_far / this.step_far; 
				this.bin_near =            this.step_threshold_far / this.step_far * (Math.log(this.step_near/this.step_far) + 1);
				int bin_max = disparityToBin(max_disparity);
				numBins = bin_max + 1;
				bin_centers = new double [numBins];
				int bin = 0;
				for (bin = 0; bin < numBins; bin ++){
					bin_centers[bin] = binToDisparity(bin);  
				}
				
				getDisparityHistograms(null); // calculate and blur supertiles (for all, not just selected?)
				if (tileProcessor.globalDebugLevel > -1){
					System.out.println("SuperTiles(): min_disparity = "+min_disparity+", max_disparity="+max_disparity);
					System.out.println("SuperTiles(): step_far = "+step_far+", step_near="+step_near);
					System.out.println("SuperTiles(): numBins = "+numBins+", bin_far="+bin_far+", bin_near = "+bin_near);
					System.out.println("SuperTiles(): step_threshold_far = "+step_threshold_far+", step_threshold_near="+step_threshold_near);
					for (int i = 0; i<numBins; i++){
						System.out.println(i+": "+bin_centers[i]+", "+disparityToBin(bin_centers[i]));
					}
					//
				}
			}
			
			public int disparityToBin(double disparity){ // not filtered
				double x = disparity - min_disparity;
				int bin;
				if (x < step_threshold_far){
					bin = (int) Math.round(x/step_far);
				} else if (x < step_threshold_near){
					bin = (int) Math.round(step_threshold_far / step_far * (Math.log(x / step_threshold_far) + 1));
				} else  {
					bin = (int) Math.round(bin_near + (x - step_threshold_near) / step_near);
				}
//				if (bin < 0) bin = 0;
//				else if (bin >= numBins) bin = numBins - 1;
				return bin;
			}
			
			public double binToDisparity(double bin){
				double d;
				if (bin < bin_far ) {
					d = bin * step_far;
				} else if (bin < bin_near){
					d =  step_threshold_far * Math.exp(bin * step_far / step_threshold_far - 1.0);
				} else {
					d = step_threshold_near + (bin - bin_near) * step_near; 
				}
				return min_disparity + d;  
			}
			
			public TilePlanes.PlaneData [][] getPlanes()
			{
				return planes;			
			}
			private double [][] getLapWeights(){
				final int superTileSize = tileProcessor.superTileSize;
				final double [][] lapWeight = new double [2 * superTileSize][2 * superTileSize];
				final double [] lapWeight1d = new double [superTileSize];
				final int superTileSize2 = 2 * superTileSize;
				for (int i = 0; i < superTileSize; i++){
					lapWeight1d[i] = 0.5*(1.0 - Math.cos((i + 0.5)* Math.PI/superTileSize));
				}
				for (int i = 0; i < superTileSize; i++){
					for (int j = 0; j < superTileSize; j++){
						lapWeight[i]                     [                     j] = lapWeight1d[i]*lapWeight1d[j]; 
						lapWeight[superTileSize2 - 1 - i][                     j] = lapWeight[i][j];
						lapWeight[i]                     [superTileSize2 - 1 - j] = lapWeight[i][j]; 
						lapWeight[superTileSize2 - 1 - i][superTileSize2 - 1 - j] = lapWeight[i][j]; 
					}
				}
				double s = 0.0;
				for (int i = 0; i < superTileSize2; i++){
					for (int j = 0; j < superTileSize2; j++){
						s+=lapWeight[i][j];
					}
				}
				System.out.println("getLapWeights: sum = "+s);
				return lapWeight;
			}
			// updates disparityHistograms
			public double [][] getDisparityHistograms()
			{
				return getDisparityHistograms(cltPass3d.selected, tileProcessor.globalDebugLevel);
			}
			// updates disparityHistograms
			public double [][] getDisparityHistograms(final boolean [] selected) //  null)
			{
				return getDisparityHistograms(selected, tileProcessor.globalDebugLevel);
			}

			public double [][] getDisparityHistograms(
					final boolean [] selected, // or null
					final int        debugLevel)
			{
				if (this.disparityHistograms != null) return this.disparityHistograms;
				final double step_disparity = step_near; // TODO: implement
				final int tilesX =        tileProcessor.getTilesX();
				final int tilesY =        tileProcessor.getTilesY();
				final int superTileSize = tileProcessor.superTileSize;
				final double  [] disparity = cltPass3d.getDisparity();
				final double  [] strength =  cltPass3d.getStrength();
				final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
				final int stilesY = (tilesY + superTileSize -1)/superTileSize;
				final int nStiles = stilesX * stilesY; 
				final double [][] dispHist =     new double [nStiles][]; // now it will be sparse 
				final double []   strengthHist = new double [nStiles];
				final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
				final AtomicInteger ai = new AtomicInteger(0);
				final int st_start = - superTileSize/2;
				final int superTileSize2 = 2 * superTileSize;
				final double [][] lapWeight = getLapWeights();
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
								int stileY = nsTile / stilesX;  
								int stileX = nsTile % stilesX;  
								double sw = 0.0; // sum weights
								double [] hist = new double [numBins];
								int tY0 = stileY * superTileSize + st_start;
								int tX0 = stileX * superTileSize + st_start;
								for (int tY = 0; tY < superTileSize2; tY++){
									int tileY = tY0 +tY;
									if ((tileY >= 0) && (tileY < tilesY)) {
										for (int tX = 0; tX < superTileSize2; tX++){
											int tileX = tX0 +tX;
											if ((tileX >= 0) && (tileX < tilesX)) {
												int indx = tileY*tilesX + tileX;
												double d = disparity[indx];
												if (!Double.isNaN(d) && ((selected == null) || selected[indx])){
													double w = strength[indx] - strength_floor;
													if (w > 0.0){
														if (strength_pow != 1.0) w = Math.pow(w, strength_pow);
														w *= lapWeight[tY][tX];
														// ignore too near/ too far
														int bin = disparityToBin(d);
														if ((bin >= 0) && (bin < numBins)){ // maybe collect below min and above max somewhere?
															hist[bin] += w; // +1]
															sw +=w;
														}
													}
												}
											}
										}
									}
								}
								strengthHist[nsTile] = sw / superTileSize / superTileSize; // average strength per tile in the super-tile
								if (sw > 0){
									for (int i = 0; i<numBins; i++){
										hist[i] /= sw; 
									}
									dispHist[nsTile] = hist;
								} else {
									dispHist[nsTile] = null;
								}
							}
						}
					};
				}		      
				ImageDtt.startAndJoin(threads);
				this.disparityHistograms = dispHist;
				this.stStrength =          strengthHist;
				if (this.stBlurSigma > 0.0) {
					blurDisparityHistogram(debugLevel);
				}
				return this.disparityHistograms; // dispHist;
			}
			

			
			public void blurDisparityHistogram( // in-place
					final int  debugLevel)
			{
				
				final double [][] dispHist = this.disparityHistograms;
				final double      sigma =    this.stBlurSigma;
				final int sTiles = dispHist.length;
//				final int numBins = dispHist[0].length;
				final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
				final AtomicInteger ai = new AtomicInteger(0);
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							DoubleGaussianBlur gb=new DoubleGaussianBlur();
							for (int nsTile = ai.getAndIncrement(); nsTile < sTiles; nsTile = ai.getAndIncrement()) {
								if (dispHist[nsTile] != null) {
									gb.blur1Direction(dispHist[nsTile], numBins, 1, sigma, 0.01,true);
								}
							}
						}
					};
				}		      
				ImageDtt.startAndJoin(threads);
			}
			// returns odd-length array of max/min (x, strength) pairs
			
			public double [][][] getMaxMinMax(){
				// first find all integer maximums, and if the top is flat - use the middle. If not flat - use 2-nd degree polynomial
				if (disparityHistograms == null) getDisparityHistograms();
				final int globalDebugLevel = tileProcessor.globalDebugLevel;
				maxMinMax = new double [disparityHistograms.length][][];
				final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
				final AtomicInteger ai = new AtomicInteger(0);
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							DoubleGaussianBlur gb=new DoubleGaussianBlur();
							for (int nsTile = ai.getAndIncrement(); nsTile < maxMinMax.length; nsTile = ai.getAndIncrement()) {
								if (disparityHistograms[nsTile] != null) {
									double [] dh = disparityHistograms[nsTile];
									double [][] mmm = new double [numBins][2]; // will definitely be shorter
									int numMax = 0;
									int lo = 0; 
									int hi = 1;
									if ((globalDebugLevel > -1 ) && (nsTile == 359)) {
										System.out.println(nsTile);
									}
									while (hi < numBins) {
										// looking for next max
										while ((hi < numBins) && (dh[hi] >= dh[hi - 1])) hi++; // flat or higher - continue
										if (hi == numBins){ // last
											if (dh[hi - 1] == dh[lo]) break; // no maximums till the very end
											if (dh[hi - 1] > dh[hi-2]) {//  and is higher than previus 
												mmm[numMax * 2][0] = hi - 1;  
											} else { // flat top, but higher than [lo]
												int i = hi - 3;
												while (dh[i] == dh[hi - 1]) i --;
												mmm[numMax * 2][0] = 0.5 * (hi + i); // middle  
											}
											mmm[numMax * 2][1] = dh[hi - 1];
											numMax++;
											break;
										} else { // d[hi] < d[hi-1] // can be 2-nd after higher first, has sharp max at hi-1, have flat max (including from start)
											if ((dh[hi - 1] == dh[lo]) || (dh[hi - 1] == dh[hi-2])){ // just a second with higher 1-st, or a flat initial top
												// flat top
												while (dh[lo] < dh[hi - 1]) lo++;
												mmm[numMax * 2][0] = 0.5 *(hi - 1 + lo);
												mmm[numMax * 2][1] = dh[hi - 1];
											} else { // normal max, use parabola for approximation and max
												double a = 0.5*(dh[hi] -2*dh[hi-1] + dh[hi-2]); 
												double b = 0.5*(dh[hi] - dh[hi-2]);
												double dx =  - b/(2*a);
												// protect agains very low a,b
												if (dx > 1.0){
													dx = 1.0;
													mmm[numMax * 2][1] = dh[hi];
													System.out.println("Max: insufficient precision, dx = "+dx+" > 1.0");
												} else if (dx < -1.0){
													dx = -1.0;
													mmm[numMax * 2][1] = dh[hi - 2];
													System.out.println("Max: insufficient precision, dx = "+dx+" < -1.0");
												} else {
													mmm[numMax * 2][1] = dh[hi-1] - b*b / (4*a);
												}
												mmm[numMax * 2][0] = (hi -1) + dx;
											}
											lo = hi-1;
											numMax++;
										}
										// look for next minimum after maximum
										while ((hi < numBins) && (dh[hi] <= dh[hi - 1])) hi++; // flat or lower - continue
										if (hi == numBins){ // last
											break;// do not need to register minimum after maximum
										} else if (dh[hi - 1] == dh[hi-2]) { // flat top
											while (dh[lo] > dh[hi - 1]) lo++;
											mmm[numMax * 2 - 1][0] = 0.5 *(hi - 1 + lo);
											mmm[numMax * 2 - 1][1] = dh[hi - 1];
										} else { // normal min - use parabola
											double a = 0.5*(dh[hi] -2*dh[hi-1] + dh[hi-2]); 
											double b = 0.5*(dh[hi] - dh[hi-2]); 
											double dx =  - b/(2*a);
											// protect agains very low a,b
											if (dx > 1.0){
												dx = 1.0;
												mmm[numMax * 2 - 1][1] = dh[hi];
												System.out.println("Min: insufficient precision, dx = "+dx+" > 1.0");
											} else if (dx < -1.0){
												dx = -1.0;
												mmm[numMax * 2 - 1][1] = dh[hi - 2];
												System.out.println("Min: insufficient precision, dx = "+dx+" < -1.0");
											} else {
												mmm[numMax * 2 -1][1] = dh[hi-1] - b*b / (4*a);
											}
											mmm[numMax * 2 - 1][0] = (hi -1) + dx;
											lo = hi -1;
										}
									}
									if (numMax > 0) {
									maxMinMax[nsTile] = new double[numMax * 2 - 1][2];
									for (int i = 0; i < maxMinMax[nsTile].length; i++){
										maxMinMax[nsTile][i][0] = binToDisparity(mmm[i][0]); // convert to actual disparity !
										maxMinMax[nsTile][i][1] = mmm[i][1] * stStrength[nsTile];
									}
									if (globalDebugLevel > 0 ) { //************
										System.out.println(nsTile+": "+dh[0]+" "+dh[1]+" "+dh[2]+" "+dh[3]+" "+dh[4]+" "+
												dh[5]+" "+dh[6]+" "+dh[7]+" "+dh[8]+" "+dh[9]+ " "+dh[10]+" "+dh[11]+" "+dh[12]+" "+dh[13]+" "+dh[14]+" "+
												dh[15]+" "+dh[16]+" "+dh[17]+" "+dh[18]+" "+dh[19]+ " "+dh[20]);
										String dbg_str = ""+nsTile+": ";
										for (int i = 0; i < maxMinMax[nsTile].length; i++){
											dbg_str += " "+maxMinMax[nsTile][i][0]+":"+maxMinMax[nsTile][i][1];
										}
										System.out.println(dbg_str);
									}
									} else {
										maxMinMax[nsTile] = null;
									}

								} else {
									maxMinMax[nsTile] = null;
								}
							}
						}
					};
				}		      
				ImageDtt.startAndJoin(threads);
				return maxMinMax;
			}
			
			
			public int showDisparityHistogramWidth()
			{
				final int superTileSize = tileProcessor.superTileSize;
				int sTilesX = (tileProcessor.getTilesX() + superTileSize -1)/superTileSize;  
				return sTilesX * (numBins + 1) + 1;
			}

			public double [] showDisparityHistogram()
			{
				if (disparityHistograms == null){
					getDisparityHistograms(); // calculate and blur with the current settings, specified at instantiation
				}
				return showDisparityHistogram(disparityHistograms);
			}
				
			public double [] showDisparityHistogram(final double [][] dispHist){
				final int superTileSize = tileProcessor.superTileSize;
				int         sTileHeight0 = 0; // vertical pixels for each  histogram (excluding borders). if <= 0 make square cells
				final double [] strengthHist = stStrength;
				final int sTilesX = (tileProcessor.getTilesX() + superTileSize -1)/superTileSize;  
				
				final int sTiles = dispHist.length;
				final int sTilesY = sTiles / sTilesX;
//				final int numBins = dispHist[0].length; // [0] - weight
				final int sTileHeight = (sTileHeight0 > 0)? sTileHeight0 : numBins; 

				final double [] maxHist = new double [sTiles];
				final int width = sTilesX * (numBins + 1) + 1;
				final int height = sTilesY * (sTileHeight + 1) +1;
				double [] rslt = new double [width*height];
				double maxW = 0.0; // use for borders between cells
				for (int i = 0; i < sTiles; i++){
					if (maxW < strengthHist[i]) maxW = strengthHist[i];
					if (dispHist[i] != null) {
						for (int j = 0; j< numBins; j++){
							if (!Double.isNaN( dispHist[i][j]) && (maxHist[i] < dispHist[i][j])) maxHist[i] = dispHist[i][j];
						}
					}
				}
				for (int nsTile = 0; nsTile < sTiles; nsTile++){
					int stileY = nsTile / sTilesX;  
					int stileX = nsTile % sTilesX;  
					int x0 = stileX * (numBins + 1);
					int y0 = stileY * (numBins + 1);
					int indx0 = x0 + y0*width; 
					
					// draw rectangular frame - horisontal dotted lines
					for (int j = 0; j < numBins + 2; j+=2) {
						rslt [indx0 + j] =                            maxW;
						rslt [indx0 + (sTileHeight + 1) * width+ j] = maxW;
					}
					// vertical dotted lines
					for (int j = 0; j < sTileHeight + 2; j+=2) {
						rslt [indx0 + j * width] =                    maxW;
						rslt [indx0 + (numBins + 1)+  j * width] =    maxW;
					}
					// draw normalized histograms, using overall weight as intensity
					if (dispHist[nsTile] != null) {
						for (int bin = 0; bin <numBins; bin ++){
							int h = (int) (sTileHeight * dispHist[nsTile][bin] / maxHist[nsTile]);
							int x = bin + 1;
							for (int j = 0; j <= h;  j++) {
								int y =  sTileHeight + 1 - j;
								rslt [indx0 + y * width + x] =           strengthHist[nsTile];
							}
						}
					}
				}
				return rslt;
			}
			public double [] showMaxMinMax(){
				if (maxMinMax == null){
					getMaxMinMax(); // calculate and blur with the current settings, specified at instantiation
				}
				final int superTileSize = tileProcessor.superTileSize;
				
				int         sTileHeight0 = 0; // vertical pixels for each  histogram (excluding borders). if <= 0 make square cells
				final double [] strengthHist = stStrength;
				final int sTilesX = (tileProcessor.getTilesX() + superTileSize -1)/superTileSize;  
				
//				final int sTiles = disparityHistograms.length;
				final int sTiles = maxMinMax.length;
				final int sTilesY = sTiles / sTilesX;
//				final int numBins = dispHist[0].length; // [0] - weight
				final int sTileHeight = (sTileHeight0 > 0)? sTileHeight0 : numBins; 

				final double [] maxHist = new double [sTiles];
				final int width = sTilesX * (numBins + 1) + 1;
				final int height = sTilesY * (sTileHeight + 1) +1;
				final boolean [] isMax = new boolean [numBins];
				final boolean [] isMin = new boolean [numBins];

				double [] rslt = new double [width*height];
				double maxW = 0.0; // use for borders between cells
				for (int i = 0; i < sTiles; i++){
					if (maxW < strengthHist[i]) maxW = strengthHist[i];
					if (disparityHistograms[i] != null) {
						for (int j = 0; j< numBins; j++){
							if (!Double.isNaN( disparityHistograms[i][j]) && (maxHist[i] < disparityHistograms[i][j])) maxHist[i] = disparityHistograms[i][j];
						}
					}
				}
				for (int nsTile = 0; nsTile < sTiles; nsTile++){
					int stileY = nsTile / sTilesX;  
					int stileX = nsTile % sTilesX;  
					int x0 = stileX * (numBins + 1);
					int y0 = stileY * (numBins + 1);
					int indx0 = x0 + y0*width; 
					
					// draw rectangular frame - horisontal dotted lines
					for (int j = 0; j < numBins + 2; j+=2) {
						rslt [indx0 + j] =                            maxW;
						rslt [indx0 + (sTileHeight + 1) * width+ j] = maxW;
					}
					// vertical dotted lines
					for (int j = 0; j < sTileHeight + 2; j+=2) {
						rslt [indx0 + j * width] =                    maxW;
						rslt [indx0 + (numBins + 1)+  j * width] =    maxW;
					}
					// draw normalized histograms, using overall weight as intensity
					if ((disparityHistograms[nsTile] != null) && (maxMinMax[nsTile] != null)) {
						for (int bin = 0; bin <numBins; bin ++){
							isMax[bin] = false;
							isMin[bin] = false;
						}
						for (int i = 0; i <maxMinMax[nsTile].length; i++){
							int imm = (int) Math.round(maxMinMax[nsTile][i][0]);
							if      (imm <0 )        imm = 0;
							else if (imm >= numBins) imm = numBins - 1;
							if ((i & 1) == 0) isMax[imm] = true;
							else              isMin[imm] = true;
						}
						
						for (int bin = 0; bin <numBins; bin ++){
							int h = (int) (sTileHeight * disparityHistograms[nsTile][bin] / maxHist[nsTile]);
							int x = bin + 1;
							int jMin = isMin[bin] ? h : 0;
							int jMax = isMax[bin] ? h : sTileHeight;
							
							if (isMax[bin] || isMin[bin]) {
								for (int j = jMin; j <= jMax;  j++) {
									int y =  sTileHeight + 1 - j;
									rslt [indx0 + y * width + x] = strengthHist[nsTile];
								}
							}
						}
					}
				}
				return rslt;
			}
			
			// updates bgDisparity, bgStrength
			public double [][] getBgDispStrength(
					final double minBgDisparity, 
					final double minBgFract)
			{
				final double step_disparity = step_near; // TODO: implement
				
				
				if (maxMinMax == null) return null;
				final int superTileSize = tileProcessor.superTileSize;
				final int globalDebugLevel = tileProcessor.globalDebugLevel;
				final int sTiles = maxMinMax.length;
				final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
				final AtomicInteger ai = new AtomicInteger(0);
				bgDisparity = new double[sTiles];
				bgStrength =  new double[sTiles];
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							for (int nsTile = ai.getAndIncrement(); nsTile < sTiles; nsTile = ai.getAndIncrement()) {
								if (nsTile == 49) { // 414){ // 331){
									System.out.println("getBgDispStrength(); nsTIle="+nsTile);
								}
								double [][] mmm = maxMinMax[nsTile];
								bgDisparity[nsTile] = Double.NaN;
								bgStrength[nsTile] =  0.0;
								if (mmm != null){
									int maxNum = 0;
									int numMax = (maxMinMax[nsTile].length + 1) / 2;
									double selStrenth = 0.0;
									int startIndex = 0;
									for (maxNum = 0; maxNum < numMax; maxNum++){
										if (mmm[2 * maxNum][0] >= minBgDisparity){
											if (selStrenth == 0.0) startIndex = maxNum; // will keep first non-zero maximum number
											selStrenth += mmm[2 * maxNum][1]; 
										}
									}
									if (selStrenth > 0.0){
										selStrenth *= minBgFract;
										double accumStrength = 0.0;
										for (maxNum = startIndex; maxNum < numMax; maxNum++){
											accumStrength += mmm[2 * maxNum][1];
											if (accumStrength >= selStrenth){
												break; 
											}
										}
										if (maxNum >= numMax){
											maxNum = numMax - 1; // probably just wrong fraction (>1.0) 
										}
										// if unlikely there are several maximums before minBgFract - use the strongest
										int maxIndex = startIndex;
										if (startIndex < maxNum){
											for (startIndex++; startIndex < numMax ;startIndex++){
												if (mmm[2 * startIndex][1] > mmm[2 * maxIndex][1]) maxIndex = startIndex;
											}
										}
										//maxIndex is what we need. Which strength to use - individual or accumulated? Use individual.
///										bgDisparity[nsTile] = min_disparity + mmm[2 * maxIndex][0] * step_disparity;
										bgDisparity[nsTile] = binToDisparity (mmm[2 * maxIndex][0]);
										bgStrength[nsTile] =  mmm[2 * maxIndex][1];
									}
								}
							}
						}
					};
				}		      
				ImageDtt.startAndJoin(threads);
				final double [][] bgDispStrength = {bgDisparity, bgStrength};
				if (globalDebugLevel > 0) {
					showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
					final int stilesX = (tileProcessor.getTilesX() + superTileSize -1)/superTileSize;  
					final int stilesY = (tileProcessor.getTilesY()  + superTileSize -1)/superTileSize;
					sdfa_instance.showArrays(bgDispStrength, stilesX, stilesY, true, "bgDispStrength");
				}
				return bgDispStrength;
			}
			
			// from per-super-tile disparity/strength interpolate per-tile disparity/strength using same sine-based window
			public double [][] getBgTileDispStrength()
			{
				if ((this.bgDisparity == null) || (this.bgStrength == null)) return null;
				final int tilesX = tileProcessor.getTilesX();
				final int tilesY = tileProcessor.getTilesY();
				final int superTileSize = tileProcessor.superTileSize;
				final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
				final int stilesY = (tilesY + superTileSize -1)/superTileSize;
				final int nStiles = stilesX * stilesY; 
				final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
				final AtomicInteger ai = new AtomicInteger(0);
				final int superTileSize2 = 2 * superTileSize;
				final double [][] lapWeight = getLapWeights();
				final double [] tileDisparity = new double [tilesY * tilesX]; // assuming all 0.0
				final double [] tileStrength =  new double [tilesY * tilesX]; // assuming all 0.0

				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
								int stileY = nsTile / stilesX;  
								int stileX = nsTile % stilesX;  
								int tY0 = stileY * superTileSize;
								int tX0 = stileX * superTileSize;
								double [][] lapWeight_dbg = lapWeight;
								if (nsTile ==  414) { //317) { // 755){
									System.out.println("getBgTileDispStrength(x): stileY="+stileY+", stileX="+stileX+" lapWeight_dbg.length="+lapWeight_dbg.length); 
								}
								for (int tY = 0; tY < superTileSize; tY++){
									int tileY = tY0 +tY;
									if (tileY < tilesY) {
										for (int tX = 0; tX < superTileSize; tX++){
											int tileX = tX0 +tX;
											if (tileX < tilesX) {
												int tIndex =  tileY * tilesX + tileX;
												// iterate for +/- 1 supterile around current, acummulate disparity/strength
												double sd = 0.0, sw = 0.0;
												for (int stDY = -1; stDY <=1; stDY ++){
													int stY =  stileY + stDY;
													int dtY =tY + superTileSize/2 -superTileSize * stDY;
													if ((stY >= 0) && (stY < stilesY) && (dtY >= 0) && (dtY < superTileSize2)){
														for (int stDX = -1; stDX <=1; stDX ++){
															int stX =  stileX + stDX;
															int dtX =tX + superTileSize/2 -superTileSize * stDX;
															if ((stX >= 0) && (stX < stilesX) && (dtX >= 0) && (dtX < superTileSize2)){
																int stIndex = stY * stilesX + stX;
																double w = bgStrength[stIndex] * lapWeight[dtY][dtX];
																if (nsTile ==  415) {
																	System.out.println("tX="+tX+", tY="+tY+", tileX="+tileX+", tileY="+tileY+" stDX="+stDX+" stDY="+stDY+
																			", bgStrength["+stIndex+"] = "+bgStrength[stIndex]+
																			", lapWeight["+dtY+"]["+dtX+"]="+lapWeight[dtY][dtX]+" stY="+stY+" stX="+stX+" w="+w);
																}
																sw += w;
																if (w >0.0) sd += w * bgDisparity[stIndex]; // so NaN will be OK 
															}
														}
													}
												}
												tileStrength[tIndex] = sw;
												if (sw > 0.0) {
													tileDisparity[tIndex] = sd/sw;
												} else {
													tileDisparity[tIndex] = Double.NaN;
												}
												if (nsTile ==  415) {
													System.out.println("sw= "+sw+", sd="+sd);
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
				final double [][] bgTileDispStrength = {tileDisparity, tileStrength};
				return bgTileDispStrength;
			}

			public void processPlanes(

					final boolean [] selected, // or null
					final double     min_disp,
					final boolean    invert_disp, // use 1/disparity
					final double     plDispNorm,
					final int        debugLevel)
			{
				
				final int tilesX =        tileProcessor.getTilesX();
				final int tilesY =        tileProcessor.getTilesY();
				final int superTileSize = tileProcessor.getSuperTileSize();
				final int tileSize =      tileProcessor.getTileSize();
				final double  [] disparity = cltPass3d.getDisparity();
				final double  [] strength =  cltPass3d.getStrength();
				final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
				final int stilesY = (tilesY + superTileSize -1)/superTileSize;
				final int nStiles = stilesX * stilesY; 
				final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
				final AtomicInteger ai = new AtomicInteger(0);
				final int st_start = -superTileSize/2;
				final int superTileSize2 = 2 * superTileSize;
				final double [][] lapWeight = getLapWeights();
				final int len2 = superTileSize2*superTileSize2;
				final double []  double_zero =  new double  [len2];
				final boolean [] boolean_zero = new boolean [len2];

				final int debug_stile = 18 * stilesX + 25; 
				
				
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							TilePlanes tpl = new TilePlanes(tileSize,superTileSize);
							double [] stDisparity = new double [superTileSize2*superTileSize2];
							double [] stStrength =  new double [superTileSize2*superTileSize2];
							boolean [] stSel = new boolean [superTileSize2*superTileSize2];
							for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
								int stileY = nsTile / stilesX;  
								int stileX = nsTile % stilesX;
								double sw = 0.0; // sum weights
								double [] hist = new double [numBins];
								int tY0 = stileY * superTileSize + st_start;
								int tX0 = stileX * superTileSize + st_start;
								System.arraycopy(double_zero,  0, stDisparity, 0, len2);
								System.arraycopy(double_zero,  0, stStrength,  0, len2);
								System.arraycopy(boolean_zero, 0, stSel,       0, len2);
								for (int tY = 0; tY < superTileSize2; tY++){
									int tileY = tY0 +tY;
									if ((tileY >= 0) && (tileY < tilesY)) {
										for (int tX = 0; tX < superTileSize2; tX++){
											int tileX = tX0 +tX;
											if ((tileX >= 0) && (tileX < tilesX)) {
												int indx = tileY*tilesX + tileX;
												double d = disparity[indx];
												if (!Double.isNaN(d) && (d >= min_disp) &&((selected == null) || selected[indx])){
													if (invert_disp){
														d = 1.0/d;
													}
													double w = strength[indx] - strength_floor;
													if (w > 0.0){
														if (strength_pow != 1.0) w = Math.pow(w, strength_pow);
														w *= lapWeight[tY][tX];
														int indx_out =  tY * superTileSize2 + tX;
														stDisparity[indx_out] = d;
														stStrength[indx_out] =  w;
														stSel[indx_out] =       true;
														sw +=w;
													}
												}
											}
										}
									}
								}
								if (sw >0){
									//									int dl = ((nsTile >= debug_stile-1) && (nsTile <= debug_stile+1) ) ? 1 : 0;
									//									int dl = (stileY == 17) ? 1 : 0;
									int dl = (stileY >= 0) ? 1 : 0;
									double [][][] rslt = tpl.getCovar(
											stDisparity,
											stStrength,
											stSel,
											plDispNorm,
											0); // dl); // debugLevel);
									if (dl > 0) {
										int       numPoints =  (int) rslt[2][0][2];
										double    kz =   rslt[2][0][1];
										double    swc =  rslt[2][0][0];
										double [] szxy = rslt[2][1];
										double [][] eig_val =  rslt[0];
										double [][] eig_vect = rslt[1];
										if (swc > 1.0) {
											System.out.println("Processing planes, nsTile="+nsTile+", stileX="+stileX+", stileY="+stileY+", numPoints="+numPoints+
													", kz = "+kz+", sw = "+sw+ ", swc = "+swc+ ", center=["+szxy[0]+","+szxy[1]+","+szxy[2]+"]"+
													", eig_val = {"+eig_val[0][0]+","+eig_val[1][1]+","+eig_val[2][2]+"}"+
													", eig_vect[0] = {"+eig_vect[0][0]+","+eig_vect[1][0]+","+eig_vect[2][0]+"}");
										}
									}

/*
		double [][][] rslt = {
				eig.getD().getArray(),
				eig.getV().getArray(),
				{
					{sw,kz},
					{swz, swx, swy}}};
		return rslt;
									
 */
									
									
									if (dl > 1) {
										System.out.println("Processing planes with average disparity");
										double ss = 0.0, sd = 0.0;
										for (int i = 0; i < stDisparity.length; i++){
											if (!Double.isNaN(stDisparity[i]) && stSel[i]){
												ss += stStrength[i];
												sd += stStrength[i] * stDisparity[i];
											}
										}
										if (ss > 0) {
											sd /= ss; 
										}
										for (int i = 0; i < stDisparity.length; i++){
											stDisparity[i] = sd;
										}
										tpl.getCovar(
												stDisparity,
												stStrength,
												stSel,
												plDispNorm,
												dl); // debugLevel);
									}
								}
							}
						}
					};
				}		      
				ImageDtt.startAndJoin(threads);
			}

			public void processPlanes1(

					final boolean [] selected, // or null
					final double     min_disp,
					final boolean    invert_disp, // use 1/disparity
					final double     plDispNorm,
					final int        debugLevel)
			{
				if (maxMinMax == null) getMaxMinMax();
				final int np_min = 5; // minimal number of points to consider
				final int tilesX =        tileProcessor.getTilesX();
				final int tilesY =        tileProcessor.getTilesY();
				final int superTileSize = tileProcessor.getSuperTileSize();
				final int tileSize =      tileProcessor.getTileSize();
				final double  [] disparity = cltPass3d.getDisparity();
				final double  [] strength =  cltPass3d.getStrength();
				final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
				final int stilesY = (tilesY + superTileSize -1)/superTileSize;
				final int nStiles = stilesX * stilesY; 
				final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
				final AtomicInteger ai = new AtomicInteger(0);
				final int st_start = -superTileSize/2;
				final int superTileSize2 = 2 * superTileSize;
				final double [][] lapWeight = getLapWeights();
				final int len2 = superTileSize2*superTileSize2;
				final double []  double_zero =  new double  [len2];
				final boolean [] boolean_zero = new boolean [len2];

				final int debug_stile = 18 * stilesX + 25; 
				
				
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							TilePlanes tpl = new TilePlanes(tileSize,superTileSize);
							double [] stDisparity = new double [superTileSize2*superTileSize2];
							double [] stStrength =  new double [superTileSize2*superTileSize2];
							boolean [] stSel = new boolean [superTileSize2*superTileSize2];
							for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
								int stileY = nsTile / stilesX;  
								int stileX = nsTile % stilesX;
								double sw = 0.0; // sum weights
								double [] hist = new double [numBins];
								int tY0 = stileY * superTileSize + st_start;
								int tX0 = stileX * superTileSize + st_start;
								System.arraycopy(double_zero,  0, stDisparity, 0, len2);
								System.arraycopy(double_zero,  0, stStrength,  0, len2);
								System.arraycopy(boolean_zero, 0, stSel,       0, len2);
								for (int tY = 0; tY < superTileSize2; tY++){
									int tileY = tY0 +tY;
									if ((tileY >= 0) && (tileY < tilesY)) {
										for (int tX = 0; tX < superTileSize2; tX++){
											int tileX = tX0 +tX;
											if ((tileX >= 0) && (tileX < tilesX)) {
												int indx = tileY*tilesX + tileX;
												double d = disparity[indx];
												if (!Double.isNaN(d) && (d >= min_disp) &&((selected == null) || selected[indx])){
													if (invert_disp){
														d = 1.0/d;
													}
													double w = strength[indx] - strength_floor;
													if (w > 0.0){
														if (strength_pow != 1.0) w = Math.pow(w, strength_pow);
														w *= lapWeight[tY][tX];
														int indx_out =  tY * superTileSize2 + tX;
														stDisparity[indx_out] = d;
														stStrength[indx_out] =  w;
														stSel[indx_out] =       true;
														sw +=w;
													}
												}
											}
										}
									}
								}
								if (sw >0){
									
									//									int dl = ((nsTile >= debug_stile-1) && (nsTile <= debug_stile+1) ) ? 1 : 0;
									//									int dl = (stileY == 17) ? 1 : 0;
									int dl = (stileY >= 0) ? 1 : 0;
									double [][][] rslt = tpl.getCovar(
											stDisparity,
											stStrength,
											stSel,
											plDispNorm,
											0); // dl); // debugLevel);
									double swc_common = 0.0;
									if (dl > 0) {
										int       numPoints =  (int) rslt[2][0][2];
										double    kz =   rslt[2][0][1];
										double    swc =  rslt[2][0][0];
										double [] szxy = rslt[2][1];
										double [][] eig_val =  rslt[0];
										double [][] eig_vect = rslt[1];
										swc_common = swc;
										if (swc > 1.0) {
											System.out.println("Processing planes, nsTile="+nsTile+", stileX="+stileX+", stileY="+stileY+", numPoints="+numPoints+
													", kz = "+kz+", sw = "+sw+ ", swc = "+swc+ ", center=["+szxy[0]+","+szxy[1]+","+szxy[2]+"]"+
													", eig_val = {"+eig_val[0][0]+","+eig_val[1][1]+","+eig_val[2][2]+"}"+
													", eig_vect[0] = {"+eig_vect[0][0]+","+eig_vect[1][0]+","+eig_vect[2][0]+"}");
										}
									}
									double [][] mm = maxMinMax[nsTile];
									
									if (mm.length > 1) { // multiple maximums - separate into multiple selections
										boolean [][] stSels = new boolean[(mm.length +1)/2][stSel.length];
										for (int i = 0; i< stSel.length; i++) if (stSel[i]){
											double d = stDisparity[i];
											loop:{
												for (int m = 0; m < (stSels.length-1); m++){
													if (d < mm[2 * m + 1][0]){
														stSels[m][i] = true;
														break loop;
													}
												}
												stSels[stSels.length-1][i] = true;
											}
										}
										double [][][][] rslts = new double [stSels.length][][][];
										int [] np = new int[stSels.length];
										String dbg_str = "";
										for (int m = 0; m < stSels.length; m++){
											for (int i =0; i < stSels[m].length; i++){
												if (stSels[m][i]) np[m]++; 
											}
											dbg_str +="  "+np[m];
											if (m < stSels.length -1){
												dbg_str +="("+mm[2*m+1][0]+")";
											}
										}
										if (swc_common > 1.0) {
											System.out.println("Processing subplanes:"+dbg_str+", nsTile="+nsTile);
										}
										for (int m = 0; m < stSels.length; m++) if (np[m] > np_min) {

											rslts[m] = tpl.getCovar(
													stDisparity,
													stStrength,
													stSels[m],
													plDispNorm,
													0); // dl); // debugLevel);
											if (dl > 0) {
												if (rslts[m] != null) {
													int       numPoints =  (int) rslts[m][2][0][2];
													double    kz =   rslts[m][2][0][1];
													double    swc =  rslts[m][2][0][0];
													double [] szxy = rslts[m][2][1];
													double [][] eig_val =  rslts[m][0];
													double [][] eig_vect = rslts[m][1];
													if (swc_common > 1.0) { // reduce output
														System.out.println("Processing subplane["+m+"], nsTile="+nsTile+", stileX="+stileX+", stileY="+stileY+", numPoints="+numPoints+
																", kz = "+kz+", sw = "+sw+ ", swc = "+swc+ ", center=["+szxy[0]+","+szxy[1]+","+szxy[2]+"]"+
																", eig_val = {"+eig_val[0][0]+","+eig_val[1][1]+","+eig_val[2][2]+"}"+
																", eig_vect[0] = {"+eig_vect[0][0]+","+eig_vect[1][0]+","+eig_vect[2][0]+"}");
													}
												} else {
													System.out.println("Processing subplane["+m+"], nsTile="+nsTile+", stileX="+stileX+", stileY="+stileY+" RETURNED NULL");
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
			}

			public void processPlanes2(

					final boolean [] selected, // or null
					final double     min_disp,
					final boolean    invert_disp, // use 1/disparity
					final double     plDispNorm,
					final int        plMinPoints, //          =     5;  // Minimal number of points for plane detection
					final double     plTargetEigen, //        =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
					final double     plFractOutliers, //      =   0.3;  // Maximal fraction of outliers to remove
					final int        plMaxOutliers, //        =    20;  // Maximal number of outliers to remove
					final GeometryCorrection geometryCorrection,
					final boolean    correct_distortions, 
					final int        debugLevel,
					final int        dbg_X,
					final int        dbg_Y)
			{
				if (maxMinMax == null) getMaxMinMax();
//				final int np_min = 5; // minimal number of points to consider
				final int tilesX =        tileProcessor.getTilesX();
				final int tilesY =        tileProcessor.getTilesY();
				final int superTileSize = tileProcessor.getSuperTileSize();
				final int tileSize =      tileProcessor.getTileSize();
				
				final double  [] disparity = cltPass3d.getDisparity();
				final double  [] strength =  cltPass3d.getStrength();
				final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
				final int stilesY = (tilesY + superTileSize -1)/superTileSize;
				final int nStiles = stilesX * stilesY; 
				final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
				final AtomicInteger ai = new AtomicInteger(0);
				final int st_start = -superTileSize/2;
				final int superTileSize2 = 2 * superTileSize;
				final double [][] lapWeight = getLapWeights();
				final int len2 = superTileSize2*superTileSize2;
				final double []  double_zero =  new double  [len2];
				final boolean [] boolean_zero = new boolean [len2];
//ArrayList<Individual>[] group = (ArrayList<Individual>[])new ArrayList[4];
//				this.planes = (ArrayList<TilePlanes.PlaneData>[]) new ArrayList[nStiles];
				this.planes = new TilePlanes.PlaneData[nStiles][];
				
//				final int debug_stile = 18 * stilesX + 25; 
//				final int debug_stile = 20 * stilesX + 24; 
//				final int debug_stile = 16 * stilesX + 27; // 10; 
				final int debug_stile = (debugLevel > -1)? (dbg_Y * stilesX + dbg_X):-1; 
				
				
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							TilePlanes tpl = new TilePlanes(tileSize,superTileSize, geometryCorrection);
							double [] stDisparity = new double [superTileSize2*superTileSize2];
							double [] stStrength =  new double [superTileSize2*superTileSize2];
							boolean [] stSel = new boolean [superTileSize2*superTileSize2];
							for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
								int stileY = nsTile / stilesX;  
								int stileX = nsTile % stilesX;
								int [] sTiles = {stileX, stileY};
								double sw = 0.0; // sum weights
								double [] hist = new double [numBins];
								int tY0 = stileY * superTileSize + st_start;
								int tX0 = stileX * superTileSize + st_start;
								System.arraycopy(double_zero,  0, stDisparity, 0, len2);
								System.arraycopy(double_zero,  0, stStrength,  0, len2);
								System.arraycopy(boolean_zero, 0, stSel,       0, len2);
								for (int tY = 0; tY < superTileSize2; tY++){
									int tileY = tY0 +tY;
									if ((tileY >= 0) && (tileY < tilesY)) {
										for (int tX = 0; tX < superTileSize2; tX++){
											int tileX = tX0 +tX;
											if ((tileX >= 0) && (tileX < tilesX)) {
												int indx = tileY*tilesX + tileX;
												double d = disparity[indx];
												if (!Double.isNaN(d) && (d >= min_disp) &&((selected == null) || selected[indx])){
													if (invert_disp){
														d = 1.0/d;
													}
													double w = strength[indx] - strength_floor;
													if (w > 0.0){
														if (strength_pow != 1.0) w = Math.pow(w, strength_pow);
														w *= lapWeight[tY][tX];
														int indx_out =  tY * superTileSize2 + tX;
														stDisparity[indx_out] = d;
														stStrength[indx_out] =  w;
														stSel[indx_out] =       true;
														sw +=w;
													}
												}
											}
										}
									}
								}
								planes[nsTile] = null;
								if (sw >0){
									ArrayList<TilePlanes.PlaneData> st_planes = new ArrayList<TilePlanes.PlaneData>();
//									int dl = ((nsTile >= debug_stile-1) && (nsTile <= debug_stile+1) ) ? 1 : 0;
//									int dl = ((stileY == 17) && (stileX > 4)) ? 1 : 0;
//									int dl = (stileY >= 0) ? 1 : 0;
									int dl1 =  (nsTile == debug_stile) ? 3 : 0;
//									int dl =  (nsTile == debug_stile) ? 3 : 0;
//									int dl = ((stileY >= 15) && (stileY <= 18) && (stileX >= 5) && (stileX <= 31)) ? 1 : 0;
//									int dl = ((stileY == 16) && (stileX == 27) ) ? 3 : 0;
									int dl =  (nsTile == debug_stile) ? 3 : 0;
//		int debugLevel1 = ((sTileXY[0] == 27) && (sTileXY[1] == 16))? 1: 0; // check why v[0][0] <0  
									
									
									boolean [] sel_all = stSel.clone();
									TilePlanes.PlaneData pd = tpl.getPlane(
											sTiles,
											stDisparity,
											stStrength,
											sel_all,
											0); // debugLevel);
									if (pd != null) {
										//correct_distortions
										double swc_common = pd.getWeight();
										if (dl > 0) {
											if (swc_common > 0.3) { // 1.0) {
												System.out.println("Processing planes["+nsTile+"]"+
														", stileX="+stileX+
														", stileY="+stileY+
														", numPoints="+ pd.getNumPoints()+
														", sw = "+swc_common+
														", swc = "+pd.getWeight()+
														", center=["+pd.getZxy()[0]+","+pd.getZxy()[1]+","+pd.getZxy()[2]+"]"+
														", eig_val = {"+pd.getValues()[0]+","+pd.getValues()[1]+","+pd.getValues()[2]+"}"+
														", eig_vect[0] = {"+pd.getVector()[0]+","+pd.getVector()[1]+","+pd.getVector()[2]+"}");
											}
										}
										// now try to remove outliers
										int max_outliers = (int) Math.round(pd.getNumPoints() * plFractOutliers);
										if (max_outliers > plMaxOutliers) max_outliers = plMaxOutliers;
										double targetV = plTargetEigen;
										double z0 = pd.getZxy()[0];
										if ((plDispNorm > 0.0) && (z0 > plDispNorm)) {
											double dd = (plDispNorm + z0)/ plDispNorm; // > 1
											targetV *= dd * dd; // > original
										}
										if (pd.getValues()[0] > targetV) {
											pd = tpl.removePlaneOutliers(
													pd,
													sTiles,
													stDisparity,
													stStrength,
													sel_all,
													targetV, // double     targetEigen, // target eigenvalue for primary axis (is disparity-dependent, so is non-constant)
													max_outliers, // int        maxRemoved,  // maximal number of tiles to remove (not a constant)
													plMinPoints,  // int        minLeft,     // minimal number of tiles to keep
													dl1); // 0); // debugLevel);
											if (dl > 0) {
												if (swc_common > 0.3) { // 1.0) {
													System.out.println("Removed outliers["+nsTile+"]"+
															", stileX="+stileX+
															", stileY="+stileY+
															", numPoints="+ pd.getNumPoints()+
															", sw = "+swc_common+
															", swc = "+pd.getWeight()+
															", center=["+pd.getZxy()[0]+","+pd.getZxy()[1]+","+pd.getZxy()[2]+"]"+
															", eig_val = {"+pd.getValues()[0]+","+pd.getValues()[1]+","+pd.getValues()[2]+"}"+
															", eig_vect[0] = {"+pd.getVector()[0]+","+pd.getVector()[1]+","+pd.getVector()[2]+"}");
												}
											}
										} // nothing to do if already OK
										if (dl > 0) {
											System.out.println("Calculating World normal["+nsTile+"]");
										}
										
										double [] norm_xyz = pd.getWorldXYZ(
												correct_distortions,
												dl);
										if (dl > 0) {
											System.out.println("World normal["+nsTile+"] = {"+
													norm_xyz[0]+", "+norm_xyz[1]+", "+norm_xyz[2]+"}");

										}
										st_planes.add(pd);
										
										// now try for each of the disparity-separated clusters
										
										double [][] mm = maxMinMax[nsTile];

										if (mm.length > 1) { // multiple maximums - separate into multiple selections
											boolean [][] stSels = new boolean[(mm.length +1)/2][stSel.length];
											for (int i = 0; i< stSel.length; i++) if (stSel[i]){
												double d = stDisparity[i];
												loop:{
													for (int m = 0; m < (stSels.length-1); m++){
														if (d < mm[2 * m + 1][0]){
															stSels[m][i] = true;
															break loop;
														}
													}
													stSels[stSels.length-1][i] = true;
												}
											}
											double [][][][] rslts = new double [stSels.length][][][];
											int [] np = new int[stSels.length];
											String dbg_str = "";
											for (int m = 0; m < stSels.length; m++){
												for (int i =0; i < stSels[m].length; i++){
													if (stSels[m][i]) np[m]++; 
												}
												dbg_str +="  "+np[m];
												if (m < stSels.length -1){
													dbg_str +="("+mm[2*m+1][0]+")";
												}
											}
											if ((dl > 0) && (swc_common > 1.0)) {
												System.out.println("Processing subplanes["+nsTile+"]:"+dbg_str);
											}
											for (int m = 0; m < stSels.length; m++) if (np[m] > plMinPoints) {
												
												pd = tpl.getPlane(
														sTiles,
														stDisparity,
														stStrength,
														stSels[m],
														0); // debugLevel);
												if (pd != null) {
													if (dl > 0) {
														if (swc_common > 1.0) {
															System.out.println("Processing subplane["+nsTile+"]["+m+"]"+
																	", stileX="+stileX+
																	", stileY="+stileY+
																	", numPoints="+ pd.getNumPoints()+
																	", sw = "+swc_common+
																	", swc = "+pd.getWeight()+
																	", center=["+pd.getZxy()[0]+","+pd.getZxy()[1]+","+pd.getZxy()[2]+"]"+
																	", eig_val = {"+pd.getValues()[0]+","+pd.getValues()[1]+","+pd.getValues()[2]+"}"+
																	", eig_vect[0] = {"+pd.getVector()[0]+","+pd.getVector()[1]+","+pd.getVector()[2]+"}");
														}
													}
													// now try to remove outliers
													max_outliers = (int) Math.round(pd.getNumPoints() * plFractOutliers);
													if (max_outliers > plMaxOutliers) max_outliers = plMaxOutliers;
													targetV = plTargetEigen;
													z0 = pd.getZxy()[0];
													if ((plDispNorm > 0.0) && (z0 > plDispNorm)) {
														double dd = (plDispNorm + z0)/ plDispNorm; // > 1
														targetV *= dd * dd; // > original
													}
													if (pd.getValues()[0] > targetV) {
														pd = tpl.removePlaneOutliers(
																pd,
																sTiles,
																stDisparity,
																stStrength,
																stSels[m],
																targetV, // double     targetEigen, // target eigenvalue for primary axis (is disparity-dependent, so is non-constant)
																max_outliers, // int        maxRemoved,  // maximal number of tiles to remove (not a constant)
																plMinPoints,  // int        minLeft,     // minimal number of tiles to keep
																dl1); // 0); // debugLevel);
														if (dl > 0) {
															if (swc_common > 1.0) {
																System.out.println("Removed outliers["+nsTile+"]["+m+"]"+
																		", stileX="+stileX+
																		", stileY="+stileY+
																		", numPoints="+ pd.getNumPoints()+
																		", sw = "+swc_common+
																		", swc = "+pd.getWeight()+
																		", center=["+pd.getZxy()[0]+","+pd.getZxy()[1]+","+pd.getZxy()[2]+"]"+
																		", eig_val = {"+pd.getValues()[0]+","+pd.getValues()[1]+","+pd.getValues()[2]+"}"+
																		", eig_vect[0] = {"+pd.getVector()[0]+","+pd.getVector()[1]+","+pd.getVector()[2]+"}");
															}
														}
													}
													norm_xyz = pd.getWorldXYZ(
															correct_distortions);
													st_planes.add(pd);
													if (dl > 0) {
														System.out.println("World normal["+nsTile+"]["+m+"] = {"+
																norm_xyz[0]+", "+norm_xyz[1]+", "+norm_xyz[2]+"}");

													}
												}
											}
										}
									}
									if (st_planes.size() > 0){
										planes[nsTile] = st_planes.toArray(new TilePlanes.PlaneData[0] );
										if (dl >0){
											System.out.println("processPlanes2(): nsTile="+nsTile);
										}
									}
								} // if sw >0
								
							}
						}
					};
				}		      
				ImageDtt.startAndJoin(threads);
			}

			public int [] getShowPlanesWidthHeight()
			{
				final int tilesX =        tileProcessor.getTilesX();
				final int tilesY =        tileProcessor.getTilesY();
				final int superTileSize = tileProcessor.getSuperTileSize();
//				final int tileSize =      tileProcessor.getTileSize();
				final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
				final int stilesY = (tilesY + superTileSize -1)/superTileSize;
				int [] wh = {
						stilesX * superTileSize,// * tileSize,
						stilesY * superTileSize};// * tileSize};
				return wh;
			}

			
			TilePlanes.PlaneData [][] getNeibPlanes(
					final int     dir,     // 0: get from up (N), 1:from NE, ... 7 - from NW
					final boolean dbgMerge, // Combine 'other' plane with current
					final int     dbg_X,
					final int     dbg_Y
					)
			{
				final int tilesX =        tileProcessor.getTilesX();
				final int tilesY =        tileProcessor.getTilesY();
				final int superTileSize = tileProcessor.getSuperTileSize();
				final int tileSize =      tileProcessor.getTileSize();
				final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
				final int stilesY = (tilesY + superTileSize -1)/superTileSize;
				final double [] nan_plane = new double [superTileSize*superTileSize];
				for (int i = 0; i < nan_plane.length; i++) nan_plane[i] = Double.NaN;
				final int [][] dirsYX = {{-1, 0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1}};
				TilePlanes.PlaneData [][] neib_planes = new TilePlanes.PlaneData[stilesY * stilesX][];
				TilePlanes tpl = null; // new TilePlanes(tileSize,superTileSize, geometryCorrection);
//				final int debug_stile = 20 * stilesX + 27;
//				final int debug_stile = 17 * stilesX + 27;
				final int debug_stile = dbg_Y * stilesX + dbg_X;
//				final int debug_stile = -1;

				for (int sty0 = 0; sty0 < stilesY; sty0++){
					int sty = sty0 + dirsYX[dir][0];
					for (int stx0 = 0; stx0 < stilesX; stx0++){
						int nsTile0 = sty0 * stilesX + stx0; // result tile
						int stx = stx0 + dirsYX[dir][1];
						int nsTile = sty * stilesX + stx; // from where to get
						if ((sty >= 0) && (sty < stilesY) && (stx >= 0) && (stx < stilesX) && (planes[nsTile] != null)) {
							neib_planes[nsTile0] = new TilePlanes.PlaneData[planes[nsTile].length];
							for (int np = 0; np < planes[nsTile].length; np++){
								GeometryCorrection geometryCorrection = planes[nsTile][np].getGeometryCorrection();
								int [] sTileXY0 = {stx0,sty0}; 
								if (tpl == null) {
									tpl = new TilePlanes(tileSize,superTileSize, geometryCorrection);
								}
								TilePlanes.PlaneData pd = tpl.new PlaneData(
										sTileXY0,
										tileSize,
										superTileSize,
										geometryCorrection);
								/*
								TilePlanes.PlaneData pd = planes[nsTile0][np].clone();
								pd.setS
								*/
								if (nsTile0 == debug_stile) {
									System.out.println("getNeibPlanes(): nsTile0="+nsTile0+", nsTile="+nsTile+", np = "+np+" dir = "+ dir ); 
								}
								TilePlanes.PlaneData other_pd = pd.getPlaneToThis(
										planes[nsTile][np], // PlaneData otherPd,
										(nsTile0 == debug_stile)? 1:0); // int       debugLevel)
								int num_this = (planes[nsTile0] == null)? 0:(planes[nsTile0].length); 
								int best_index = 0;
								if (dbgMerge && (num_this > 0)) {
									double best_value = Double.NaN;
									if (num_this > 1){
										// find the best candidate for merge
										for (int indx = 1; indx < planes[nsTile0].length; indx++){
											if (nsTile0 == debug_stile) {
												System.out.println("this with this, same weight:");
												TilePlanes.PlaneData this_this_pd = planes[nsTile0][indx].mergePlaneToThis(
														planes[nsTile0][indx], // PlaneData otherPd,
														1.0,      // double    scale_other,
														true,     // boolean   ignore_weights,
														(nsTile0 == debug_stile)? 1:0); // int       debugLevel)
												System.out.println("other with other, same weight:");
												TilePlanes.PlaneData other_other_pd = other_pd.mergePlaneToThis(
														other_pd, // PlaneData otherPd,
														1.0,      // double    scale_other,
														true,     // boolean   ignore_weights,
														(nsTile0 == debug_stile)? 1:0); // int       debugLevel)
												System.out.println("other with this, same weight:");
											}
											TilePlanes.PlaneData merged_pd = planes[nsTile0][indx].mergePlaneToThis(
													other_pd, // PlaneData otherPd,
													1.0,      // double    scale_other,
													true,     // boolean   ignore_weights,
													(nsTile0 == debug_stile)? 1:0); // int       debugLevel)
											if (!(merged_pd.getValue() > best_value)) { // Double.isNaN(best_value) will work too
												best_value = merged_pd.getValue();
												best_index = indx;
											}
										}
									}
									// now merge with weights with the best plane of this supertile
									if (nsTile0 == debug_stile) {
										System.out.println("this with this, proportional weight, np(other):"+np+" best fit(this):"+best_index);
									}									
									TilePlanes.PlaneData merged_pd = planes[nsTile0][best_index].mergePlaneToThis(
											other_pd, // PlaneData otherPd,
											1.0,      // double    scale_other,
											false,    // boolean   ignore_weights,
											(nsTile0 == debug_stile)? 1:0); // int       debugLevel)
									neib_planes[nsTile0][np] = merged_pd;
								} else  {
									neib_planes[nsTile0][np] = other_pd;
								}
							}
						} else {
							neib_planes[nsTile0] = null;
						}
					}
				}
				return neib_planes;
			}
			
			public double mergeRQuality(
					double L1,
					double L2,
					double L,
					double w1,
					double w2)
			{
//				double Lav = Math.sqrt((L1*L1*w1 + L2*L2*w2)/(w1+w2));
				double Lav = (L1*w1 + L2*w2)/(w1+w2);
///				double wors = (L - Lav)*(w1+w2)*(w1+w2) /(Lav*w1*w2);
				double rquality = (L - Lav)*(w1+w2) /(Lav*w1*w2); // ==wors/(w1+w2) to amplify stronger planes
				return rquality;
			}
			
			public void matchPlanes(
					final int dbg_X,
					final int dbg_Y)
			{
				final int tilesX =        tileProcessor.getTilesX();
				final int tilesY =        tileProcessor.getTilesY();
				final int superTileSize = tileProcessor.getSuperTileSize();
//				final int tileSize =      tileProcessor.getTileSize();
				final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
				final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
				final int nStiles =       stilesX * stilesY; 
				final double [] nan_plane = new double [superTileSize*superTileSize];
				for (int i = 0; i < nan_plane.length; i++) nan_plane[i] = Double.NaN;
				final int [][] dirsYX = {{-1, 0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1}};
//				final int debug_stile = 20 * stilesX + 27;
				final int debug_stile = dbg_Y * stilesX + dbg_X;
				
				final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
				final AtomicInteger ai = new AtomicInteger(0);
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							for (int nsTile0 = ai.getAndIncrement(); nsTile0 < nStiles; nsTile0 = ai.getAndIncrement()) {
								int sty0 = nsTile0 / stilesX;  
								int stx0 = nsTile0 % stilesX;
								int dl = (nsTile0 == debug_stile) ? 1:0;
								if ( planes[nsTile0] != null) {
									if (dl > 0){
										System.out.println("matchPlanes(): nsTile0 ="+nsTile0);
									}
									for (int np0 = 0; np0 < planes[nsTile0].length; np0++){ // nu
//										planes[nsTile0][np0].initNeibBest(); // 
										TilePlanes.PlaneData this_plane = planes[nsTile0][np0];
										this_plane.initNeibMatch();
										for (int dir = 0; dir < 4; dir++){ // just half directions - relations are symmetrical
											int stx = stx0 + dirsYX[dir][1];
											int sty = sty0 + dirsYX[dir][0];
//											if ((sty < stilesY) && (sty > 0) && (stx < 0)) {
											if ((stx < stilesX) && (sty < stilesY) && (sty > 0)) {
												int nsTile = sty * stilesX + stx; // from where to get
												if (nsTile >= planes.length){
													System.out.println("BUG!!!!");
												} else {
													TilePlanes.PlaneData [] other_planes = planes[nsTile];
													if (other_planes != null) {
														
														this_plane.initNeibMatch(dir,other_planes.length); // filled with NaN
														for (int np = 0; np < other_planes.length; np ++){
															TilePlanes.PlaneData other_plane = this_plane.getPlaneToThis(
																	other_planes[np],
																	dl); // debugLevel);
															if (other_plane !=null) { // now always, but may add later
																TilePlanes.PlaneData merged_pd = this_plane.mergePlaneToThis(
																		other_plane, // PlaneData otherPd,
																		1.0,         // double    scale_other,
																		false,       // boolean   ignore_weights,
																		dl); // int       debugLevel)
																if (merged_pd !=null) { // now always, but may add later
																	// calculate worsening of the merge (the lower - the better)
																	double rquality = mergeRQuality(
																			this_plane.getValue(), // double L1,
																			other_plane.getValue(), // double L2,
																			merged_pd.getValue(), // double L,
																			this_plane.getWeight(), // double w1,
																			other_plane.getWeight()); // double w2)
																	this_plane.setNeibMatch(dir, np, rquality);
																	
																	// still happens negative but only for very bad planes - investigate?
																	if ((rquality < 0) && (this_plane.getWeight() > 0.0)  && (other_plane.getWeight() > 0.0)  && (merged_pd.getValue() < 1.0) ){ // 
																		System.out.println("nsTile0="+nsTile0+":"+np0+", nsTile="+nsTile+":"+np+", rquality="+rquality+
																				" w1="+this_plane.getWeight()+" w2="+other_plane.getWeight()+
																				" L1="+this_plane.getValue()+" L2="+other_plane.getValue()+" L="+merged_pd.getValue());
																	}
																	if ((merged_pd.getValue() < 0.1) && (this_plane.getWeight() > 1.0)  && (other_plane.getWeight() > 1.0) ){ // fixed
																	System.out.println("nsTile0="+nsTile0+":"+np0+", nsTile="+nsTile+":"+np+", rquality="+rquality+
																			" w1="+this_plane.getWeight()+" w2="+other_plane.getWeight()+
																			" L1="+this_plane.getValue()+" L2="+other_plane.getValue()+" L="+merged_pd.getValue());
																	}
																}
															}
														}
													}
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
				// copy symmetrical relations
				ai.set(0);				
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							TilePlanes.PlaneData [][] dbg_planes = planes;
							for (int nsTile0 = ai.getAndIncrement(); nsTile0 < nStiles; nsTile0 = ai.getAndIncrement()) {
								int sty0 = nsTile0 / stilesX;  
								int stx0 = nsTile0 % stilesX;
								int dl = (nsTile0 == debug_stile) ? 1:0;
								if (dl>0) {
									System.out.println("matchPlanes() nsTile0="+nsTile0);
								}
								if ( planes[nsTile0] != null) {
									for (int np0 = 0; np0 < planes[nsTile0].length; np0++){ // nu
										TilePlanes.PlaneData this_plane = planes[nsTile0][np0];
//										this_plane.initNeibMatch();
										for (int dir = 4; dir < 8; dir++){ // other half - copy from opposite
											int stx = stx0 + dirsYX[dir][1];
											int sty = sty0 + dirsYX[dir][0];
//											if ((stx < stilesX) && (sty < stilesY) && (sty > 0)) {
											if ((sty < stilesY) && (sty > 0) && (stx > 0)) {
												
												int nsTile = sty * stilesX + stx; // from where to get
												TilePlanes.PlaneData [] other_planes = planes[nsTile];
												if (other_planes !=null) {
													this_plane.initNeibMatch(dir,other_planes.length); // filled with NaN
													for (int np = 0; np < other_planes.length; np ++){
														/*
														if ((other_planes[np] != null) && (other_planes[np].getNeibMatch(dir-4) != null)) {
															double [] dbg_nm = other_planes[np].getNeibMatch(dir-4);
															this_plane.setNeibMatch(dir,np, other_planes[np].getNeibMatch(dir-4, np0)); // 
														}
														*/
														if (other_planes[np] != null) { // && (other_planes[np].getNeibMatch(dir-4) != null)) {
															double [] nm = other_planes[np].getNeibMatch(dir-4);
															if (nm != null) {
//																this_plane.setNeibMatch(dir,np, other_planes[np].getNeibMatch(dir-4, np0)); //
																this_plane.setNeibMatch(dir,np, nm[np0]); //
															}
														}
														
													}
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
			}

			
			public void selectNeighborPlanes(
					final double worst_worsening,
					final boolean mutual_only,
					final int debugLevel)
			{
				final int tilesX =        tileProcessor.getTilesX();
				final int tilesY =        tileProcessor.getTilesY();
				final int superTileSize = tileProcessor.getSuperTileSize();
//				final int tileSize =      tileProcessor.getTileSize();
				final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
				final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
				final int nStiles =       stilesX * stilesY; 
				final double [] nan_plane = new double [superTileSize*superTileSize];
				for (int i = 0; i < nan_plane.length; i++) nan_plane[i] = Double.NaN;
				final int [][] dirsYX = {{-1, 0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1}};
//				final int debug_stile = 20 * stilesX + 27;
//				final int debug_stile = 17 * stilesX + 27;
				
				final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
				final AtomicInteger ai = new AtomicInteger(0);
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							for (int nsTile0 = ai.getAndIncrement(); nsTile0 < nStiles; nsTile0 = ai.getAndIncrement()) {
//								int sty0 = nsTile0 / stilesX;  
//								int stx0 = nsTile0 % stilesX;
								// int dl = (nsTile0 == debug_stile) ? 1:0;
								if ( planes[nsTile0] != null) {
									int np0_min = (planes[nsTile0].length > 1) ? 1:0; // Modify if overall plane will be removed
									for (int np0 = np0_min; np0 < planes[nsTile0].length; np0++){ // nu
										TilePlanes.PlaneData this_plane = planes[nsTile0][np0];
										this_plane.initNeibBest();
										if (this_plane.getNeibMatch() != null) {
											for (int dir = 0; dir < 8; dir++){ //
												double [] neib_worse = this_plane.getNeibMatch(dir);
												if (neib_worse != null) {
													int best_index = -1;  
													int np_min = (neib_worse.length > 1) ? 1:0; // Modify if overall plane will be removed
													for (int np = np_min; np < neib_worse.length; np++){
														if (	((worst_worsening == 0.0) || (neib_worse[np] < worst_worsening)) &&
																!Double.isNaN(neib_worse[np]) &&
																((best_index < 0) || (neib_worse[np] < neib_worse[best_index]))){
															best_index = np;
														}
													}													
													if (best_index >= 0){
														this_plane.setNeibBest(dir, best_index);
													}
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
				if (mutual_only) { // resolve love triangles (by combined strength)
					ai.set(0);				
					for (int ithread = 0; ithread < threads.length; ithread++) {
						threads[ithread] = new Thread() {
							public void run() {
								for (int nsTile0 = ai.getAndIncrement(); nsTile0 < nStiles; nsTile0 = ai.getAndIncrement()) {
									int sty0 = nsTile0 / stilesX;  
									int stx0 = nsTile0 % stilesX;
//									int dl = (nsTile0 == debug_stile) ? 1:0;
									if ( planes[nsTile0] != null) {
										for (int np0 = 0; np0 < planes[nsTile0].length; np0++){ // nu
											TilePlanes.PlaneData this_plane = planes[nsTile0][np0];
											if (this_plane.getNeibBest() != null) {
												for (int dir = 0; dir < 4; dir++){ // just half directions - relations are symmetrical
													int other_np = this_plane.getNeibBest(dir);
													if (other_np >= 0){ // may be -1, but the opposite is pointing to this
														int stx = stx0 + dirsYX[dir][1];
														int sty = sty0 + dirsYX[dir][0]; // no need to check limits as other_np >= 0
														int nsTile = sty * stilesX + stx; // from where to get
														TilePlanes.PlaneData [] other_planes = planes[nsTile]; // known != null
														int other_other_np = other_planes[other_np].getNeibBest(dir + 4);
														// see if there is not a mutual love
														if (other_other_np != np0){
															if (debugLevel > -1){
																System.out.println("Planes not mutual: "+nsTile0+":"+np0+":"+dir+" -> "+
																		nsTile+":"+other_np+" -> "+other_other_np);
															}
															if (other_other_np >= 0){
																// Which pair is stronger (not comparing worsening)
																double w_this = this_plane.getWeight()+other_planes[other_np].getWeight();
																double w_other = other_planes[other_np].getWeight() + planes[nsTile0][other_other_np].getWeight();
																if (w_this > w_other){ // kill other's love
//																	other_planes[other_np].setNeibBest(dir + 4, -1);
																	other_planes[other_np].setNeibBest(dir + 4, np0);
																} else { // kill my love
//																	this_plane.setNeibBest(dir, -1);
																	this_plane.setNeibBest(dir, -1);
																}
															} else { // other points nowhere, point it to this
																other_planes[other_np].setNeibBest(dir + 4, np0);
																
															}
														}
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
				}
			}

			
			public void selectNeighborPlanesMutual(
					final double rquality,
//					final boolean mutual_only,
					final double maxEigen, // maximal eigenvalue of planes to consider
					final double minWeight, // minimal pain weight to consider
					final int debugLevel,
					final int dbg_X,
					final int dbg_Y)
			{
				final int tilesX =        tileProcessor.getTilesX();
				final int tilesY =        tileProcessor.getTilesY();
				final int superTileSize = tileProcessor.getSuperTileSize();
//				final int tileSize =      tileProcessor.getTileSize();
				final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
				final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
				final int nStiles =       stilesX * stilesY; 
				final double [] nan_plane = new double [superTileSize*superTileSize];
				for (int i = 0; i < nan_plane.length; i++) nan_plane[i] = Double.NaN;
				final int [][] dirsYX = {{-1, 0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1}};
//				final int debug_stile = 20 * stilesX + 27;
//				final int debug_stile = 17 * stilesX + 27;
//				final int debug_stile = 9 * stilesX + 26;
				final int debug_stile = dbg_Y * stilesX + dbg_X;
				
				final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
				final AtomicInteger ai = new AtomicInteger(0);
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							TilePlanes.PlaneData [][] dbg_planes = planes;
							for (int nsTile0 = ai.getAndIncrement(); nsTile0 < nStiles; nsTile0 = ai.getAndIncrement()) {
								int sty0 = nsTile0 / stilesX;  
								int stx0 = nsTile0 % stilesX;
								int dl = (nsTile0 == debug_stile) ? 1:0;
								if ( planes[nsTile0] != null) {
									if (dl > 0){
										System.out.println("selectNeighborPlanesMutual() nsTile0="+nsTile0);
									}
									int np0_min = (planes[nsTile0].length > 1) ? 1:0; // Modify if overall plane will be removed
									for (int dir = 0; dir < 4; dir++){ //
										int stx = stx0 + dirsYX[dir][1];
										int sty = sty0 + dirsYX[dir][0];
										int nsTile = sty * stilesX + stx; // from where to get
										
										int num_other_planes = 0;
										for (int np0 = np0_min; np0 < planes[nsTile0].length; np0++){
											if (planes[nsTile0][np0].getNeibMatch(dir) != null){
												int l = planes[nsTile0][np0].getNeibMatch(dir).length;
												if (l > num_other_planes)num_other_planes = l;
											}
										}
										if (num_other_planes > 0){ // will eliminate bad margins
											int np_min = (num_other_planes > 1) ? 1:0; // Modify if overall plane will be removed
											boolean [] this_matched = new boolean [planes[nsTile0].length];
											boolean [] other_matched = new boolean [num_other_planes];
											int num_pairs = this_matched.length - np0_min;
											if ((other_matched.length - np_min) < num_pairs) num_pairs = other_matched.length - np_min;
											
											for (int pair = 0; pair < num_pairs; pair ++){
												int [] best_pair = {-1,-1};
												double best_rqual = Double.NaN;
												for (int np0 = np0_min; np0 < this_matched.length; np0++){
													double [] this_rq = planes[nsTile0][np0].getNeibMatch(dir);
													if (!this_matched[np0] &&
															(this_rq != null) &&
															((maxEigen == 0.0) || (planes[nsTile0][np0].getValue() < maxEigen)) &&
															(planes[nsTile0][np0].getWeight() > minWeight)) {
														for (int np = np_min; np < this_rq.length; np++){
															if (!other_matched[np] &&
																	!Double.isNaN(this_rq[np]) &&
																	((maxEigen == 0.0) || (planes[nsTile][np].getValue() < maxEigen)) &&
																	(planes[nsTile][np].getWeight() > minWeight)) {
																if (Double.isNaN(best_rqual) || (this_rq[np] < best_rqual)){ // OK if Double.isNaN(this_rq[np])
																	best_rqual = this_rq[np];
																	best_pair[0]= np0;
																	best_pair[1]= np;
																}
															}															
														}														
													}
												}
												if (Double.isNaN(best_rqual)){
													break; // nothing found
												}
												this_matched[best_pair[0]] = true;
												other_matched[best_pair[1]] = true;
												// neib_best should be initialized as  int [8];
//												planes[nsTile0][0].initNeibBest();
												planes[nsTile0][best_pair[0]].setNeibBest(dir,best_pair[1]);
												planes[nsTile][best_pair[1]].setNeibBest(dir + 4,best_pair[0]);
											}
											// disable remaining neighbors
											for (int np = 0; np < this_matched.length; np++){
												if (!this_matched[np]){
													planes[nsTile0][np].setNeibBest(dir,-1);
												}
											}
											for (int np = 0; np < other_matched.length; np++){
												if (!other_matched[np]){
													planes[nsTile][np].setNeibBest(dir + 4,-1);
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
			}
			
			
			
			public double [][] getShowPlanes(
					TilePlanes.PlaneData [][] planes,
					double minWeight,
					double maxEigen,
					double dispNorm,
					boolean use_NaN,
					double arrow_dark,
					double arrow_white)
			{
				final int tilesX =        tileProcessor.getTilesX();
				final int tilesY =        tileProcessor.getTilesY();
				final int superTileSize = tileProcessor.getSuperTileSize();
//				final int tileSize =      tileProcessor.getTileSize();
				final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
				final int stilesY = (tilesY + superTileSize -1)/superTileSize;
				final int width =  stilesX * superTileSize; // * tileSize;
				final int height = stilesY * superTileSize; // * tileSize;
				final double [] nan_plane = new double [superTileSize*superTileSize];
				for (int i = 0; i < nan_plane.length; i++) nan_plane[i] = Double.NaN;
				final int centerIndex = (superTileSize+1) * superTileSize / 2;
				final int [] dirs = {-superTileSize,
						             -superTileSize + 1,
						             1,
						             superTileSize +  1,
						             superTileSize,
						             superTileSize -  1,
						             -1,
						             -superTileSize - 1};
//				final int debug_stile = 17 * stilesX + 27;
				final int debug_stile = 20 * stilesX + 27;
//				final int debug_stile = 19 * stilesX + 27;
				// scan all get maximal number of satisfying planes
				int num_planes = 0;
				for (int i = 0; i < planes.length; i++) {
					if (i == debug_stile) {
						System.out.println("getShowPlanes():1, tile = "+i);
					}
					if ((planes[i] != null) && (planes[i].length > num_planes)){
						int num_sat = 0;
						for (int j = 0; j < planes[i].length; j++){
							if (planes[i][j].getWeight() >= minWeight){
								double eigVal = planes[i][j].getValue(); // ************** generated does not have value ????????????
								double disp = planes[i][j].getPlaneDisparity(false)[centerIndex];
								if (disp > dispNorm) {
									eigVal *= dispNorm / disp;
								}
								if (eigVal < maxEigen) num_sat ++;
							}
						}
						if (num_sat > num_planes) num_planes = num_sat;
					}
				}
				double [][] data = new double[num_planes][width*height];
				for (int sty = 0; sty < stilesY; sty++){
					for (int stx = 0; stx < stilesX; stx++){
						int nsTile = sty * stilesX + stx;
						if (nsTile == debug_stile) {
							System.out.println("getShowPlanes():2, tile = "+nsTile);
						}
						int ns = 0;
						double [] plane = nan_plane;
						int indx = (sty * superTileSize) * width + (stx * superTileSize);
						boolean debug = (nsTile == debug_stile) ; // (sty == 18) && (stx == 28);
						if (debug) {
							System.out.println("getShowPlanes():"+((planes[nsTile] != null)? planes[nsTile].length : "null"));
						}
						if (planes[nsTile] != null) {
							for (int np = 0; np < planes[nsTile].length; np++){
								if (planes[nsTile][np].getWeight() >= minWeight){
									//&& (planes[nsTile][np].getValue() < maxEigen)){
									double eigVal = planes[nsTile][np].getValue();
//									double disp = planes[nsTile][np].getPlaneDisparity(false)[centerIndex];
									double disp = planes[nsTile][np].getZxy()[0];
									if (disp > dispNorm) {
										eigVal *= dispNorm / disp;
									}
									if (eigVal < maxEigen) {
										plane = planes[nsTile][np].getPlaneDisparity(use_NaN);
										// add connections to neighbors if available
										int [] neib_best = planes[nsTile][np].getNeibBest();
										if (neib_best != null){
											for (int dir = 0; dir < neib_best.length; dir ++) if (neib_best[dir] >= 0){
												// draw a line
												//int center_index = indx + centerIndex;
												for (int i = 0; i < (superTileSize / 2); i++){
													plane[centerIndex + i * dirs[dir]] = ((i & 1) == 0) ? arrow_white:arrow_white; //  arrow_dark : arrow_white; 
												}
											}
										}
										for (int y = 0; y < superTileSize; y++) {
											System.arraycopy(plane, superTileSize * y, data[ns], indx + width * y, superTileSize);
										}
										ns ++;
									}
								}
							}
						}
						// fill same (nearest) plane or NaN if none was available
						for (; ns < num_planes; ns++){
							for (int y = 0; y < superTileSize; y++) {
								System.arraycopy(plane, superTileSize * y, data[ns], indx + width * y, superTileSize);
							}
						}
						
					}
				}
				return data;
			}

		} // end of class SuperTiles
