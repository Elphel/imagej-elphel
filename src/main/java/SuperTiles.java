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
import java.util.concurrent.atomic.AtomicInteger;

public class SuperTiles{
	        TileProcessor tileProcessor;
			double      step_near;
			double      step_far;
			double      step_threshold;
			double      min_disparity;
			double      max_disparity;
			double      strength_floor;
			double      strength_pow;
			double      stBlurSigma;
			int         numBins;
			double [][] disparityHistograms = null;
			double []   stStrength =          null; // per super-tile correlation strength
			double [][][] maxMinMax = null;
			double []   bgDisparity = null;
			double []   bgStrength = null;
//			TileProcessor.CLTPass3d cltPass3d;
			CLTPass3d cltPass3d;
			public SuperTiles(
//					TileProcessor.CLTPass3d cltPass3d,
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
				this.cltPass3d =      cltPass3d;
				this.tileProcessor =  cltPass3d.getTileProcessor();
				this.step_near =      step_near;
				this.step_far =       step_far;
				this.step_threshold = step_threshold;
				this.min_disparity =  min_disparity;
				this.max_disparity =  max_disparity;
				this.strength_floor = strength_floor;
				this.strength_pow =   strength_pow;
				this.stBlurSigma =    stBlurSigma;
//				this.numBins = (int) ((max_disparity - min_disparity)/step_disparity) + 1;
				this.numBins = (int) ((max_disparity - min_disparity)/step_far) + 1;
//				getDisparityHistograms(); // calculate and blur supertileas
				getDisparityHistograms(null); // calculate and blur supertileas (for all, not just selected?)
				
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
//				final int numBins = (int) ((max_disparity - min_disparity)/step_disparity) + 1;
				final double dMin = min_disparity - step_disparity/2; 
				final double dMax = dMin + numBins * step_disparity;
//				final double [][] dispHist =     new double [nStiles][numBins];
				final double [][] dispHist =     new double [nStiles][]; // now it will be sparse 
				final double []   strengthHist = new double [nStiles];
				final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
				final AtomicInteger ai = new AtomicInteger(0);
				final int st_start = - superTileSize/2;
//				final int st_end = st_start + superTileSize;
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
//								for (int tileY = stileY * superTileSize + st_start; tileY < stileY * superTileSize + st_end; tileY++){
								int tY0 = stileY * superTileSize + st_start;
								int tX0 = stileX * superTileSize + st_start;
								for (int tY = 0; tY < superTileSize2; tY++){
									int tileY = tY0 +tY;
									if ((tileY >= 0) && (tileY < tilesY)) {
//										for (int tileX = stileX * superTileSize + st_start; tileX < stileX * superTileSize + st_end; tileX++){
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
														int bin = (int) ((d-dMin)/step_disparity);
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
										maxMinMax[nsTile][i][0] = mmm[i][0];
										maxMinMax[nsTile][i][1] = mmm[i][1] * stStrength[nsTile];
									}
									if (globalDebugLevel > -1 ) {
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
										bgDisparity[nsTile] = min_disparity + mmm[2 * maxIndex][0] * step_disparity;
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
				//				final int numBins = (int) ((max_disparity - min_disparity)/step_disparity) + 1;
				//				final double [][] dispHist =     new double [nStiles][numBins];
				final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
				final AtomicInteger ai = new AtomicInteger(0);
				//				final int st_end = st_start + superTileSize;
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




		} // end of class SuperTiles
