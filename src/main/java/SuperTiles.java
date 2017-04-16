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
import java.util.Collections;
import java.util.Comparator;
import java.util.concurrent.atomic.AtomicInteger;
import java.awt.Point;
import java.util.ArrayList;

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
	TilePlanes.PlaneData [][] planes =  null;
	// 
	TilePlanes.PlaneData [][] planes_mod = null;
	double [][][] fuse_coeff = null; // coefficients to fuse planes in neighbor supertiles

	int [][]    shell_map = null;   // per supertile, per disparity plane - shell index + 1 (0 - none)
	double [][] surfaces;  // per shell, per tile (linescan order) disparity value or NaN in missing supertiles


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
		initFuseCoeff(0.5, false); // true);
	}

	public void initFuseCoeff(
			double scale_diag,
			boolean debug)
	{
		final int superTileSize = tileProcessor.getSuperTileSize();

		fuse_coeff = new double [8][2][superTileSize * superTileSize];
		double [] sin2orthotab = new double [2* superTileSize];
		for (int i = 0; i < 2 * superTileSize; i++){
			sin2orthotab[i] = 0.5 * (1.0 - Math.cos(Math.PI*(i + 0.5)/superTileSize));
		}
		double [] sin2diagtab = new double [4* superTileSize];
		for (int i = 0; i < 4 * superTileSize; i++){
			sin2diagtab[i] = 0.5 * (1.0 - Math.cos(Math.PI*i/(2 * superTileSize)));
		}
		//  ortho neighbors N, E, S, W
		for (int y = 0; y < superTileSize /2 ; y++){
			for (int x = 0; x < superTileSize; x++){
				int [] indx = {
						y * superTileSize +                       x,
						x * superTileSize +                       (superTileSize - y - 1),
						(superTileSize - y - 1) * superTileSize + (superTileSize - x - 1),
						(superTileSize - x - 1) * superTileSize + y,
				};
				for (int dir = 0; dir <4; dir++){
					fuse_coeff[2*dir][0][indx[dir]] = sin2orthotab[x + superTileSize/2] * sin2orthotab[y + superTileSize / 2];  
					fuse_coeff[2*dir][1][indx[dir]] = sin2orthotab[x + superTileSize/2] * sin2orthotab[y + 3 * superTileSize / 2];  
				}
			}
		}
		// Diagonal: NE, SE, SW, NW
		for (int y = 0; y < superTileSize ; y++){
			for (int x = y; x < superTileSize; x++){
				int [] indx = {
						y * superTileSize +                       x,
						x * superTileSize +                       (superTileSize - y - 1),
						(superTileSize - y - 1) * superTileSize + (superTileSize - x - 1),
						(superTileSize - x - 1) * superTileSize + y,
				};

				for (int dir = 0; dir <4; dir++){
					fuse_coeff[2*dir + 1][0][indx[dir]] = scale_diag *
							sin2diagtab[x + y + superTileSize + 1] *
							sin2diagtab[ x - y + 2* superTileSize];  
					fuse_coeff[2*dir + 1][1][indx[dir]] = scale_diag *
							sin2diagtab[x + y + superTileSize + 1] *
							sin2diagtab[ x - y];  
				}
			}
		}
		if (debug) {
			String [] titles = {
					"N_this",  "N-other",
					"NE_this", "NE-other",
					"E_this",  "E-other",
					"SE_this", "SE-other",
					"S_this",  "S-other",
					"SW_this", "SW-other",
					"W_this",  "W-other",
					"NW_this", "NW-other"};
			double [][] dbg_img = {
					fuse_coeff[0][0],fuse_coeff[0][1],
					fuse_coeff[1][0],fuse_coeff[1][1],
					fuse_coeff[2][0],fuse_coeff[2][1],
					fuse_coeff[3][0],fuse_coeff[3][1],
					fuse_coeff[4][0],fuse_coeff[4][1],
					fuse_coeff[5][0],fuse_coeff[5][1],
					fuse_coeff[6][0],fuse_coeff[6][1],
					fuse_coeff[7][0],fuse_coeff[7][1]};
			showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
			sdfa_instance.showArrays(dbg_img, superTileSize, superTileSize, true, "fuse_coeff",titles);
		}
	}

	public double [] getFuseThis(int dir){
		return fuse_coeff[dir][0];
	}
	public double [] getFuseOther(int dir){
		return fuse_coeff[dir][1];
	}

	public void setShellMap(int [][] shells){
		shell_map = shells;
	}

	public int [][] getShellMap(){
		return shell_map;
	}

	public void setSurfaces(double [][] surfaces){
		this.surfaces = surfaces;
	}

	public double [][] getSurfaces(){
		return surfaces;
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

	public TilePlanes.PlaneData [][] getPlanesMod()
	{
		return planes_mod;			
	}

	public void resetPlanesMod()
	{
		planes_mod = null;			
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

	public double corrMaxEigen(
			double maxEigen,
			double dispNorm,
			TilePlanes.PlaneData pd)
	{
		return corrMaxEigen(
				maxEigen,
				dispNorm,
				pd.getZxy()[0]);
	}

	public double corrMaxEigen(
			double maxEigen,
			double dispNorm,
			double disparity)
	{
		double corrV = maxEigen;
		if ((dispNorm > 0.0) && (disparity > dispNorm)) {
			//					double dd = (dispNorm + z0)/ dispNorm; // > 1
			double dd = disparity/ dispNorm; // > 1
			corrV *= dd * dd; // > original
		}
		return corrV;

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
			final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
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
						if (nsTile == debug_stile){
							System.out.println("processPlanes2(): nsTile="+nsTile);
						}
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
						if (sw >0){ // there are some non-zero tiles, process them (all points, not clustered by disparity value)
							ArrayList<TilePlanes.PlaneData> st_planes = new ArrayList<TilePlanes.PlaneData>();
							int dl1 =  (nsTile == debug_stile) ? 3 : 0;
							int dl =  (nsTile == debug_stile) ? 3 : 0;
							boolean [] sel_all = stSel.clone();
							TilePlanes.PlaneData pd = tpl.getPlane(
									sTiles,
									stDisparity,
									stStrength,
									sel_all,
									plPreferDisparity,									
									0); // debugLevel);
							if (pd != null) { // not too few points, probably
								//correct_distortions
								double swc_common = pd.getWeight();
								if (dl > 0) {
									if (swc_common > 0.1) { // 1.0) {
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
								double targetV = corrMaxEigen(
										plTargetEigen,
										plDispNorm,
										pd);
								if (pd.getValue() > targetV) {
									pd = tpl.removePlaneOutliers(
											pd,
											sTiles,
											stDisparity,
											stStrength,
											sel_all,
											targetV, // double     targetEigen, // target eigenvalue for primary axis (is disparity-dependent, so is non-constant)
											max_outliers, // int        maxRemoved,  // maximal number of tiles to remove (not a constant)
											plMinPoints,  // int        minLeft,     // minimal number of tiles to keep
											plPreferDisparity,
											dl1); // 0); // debugLevel);
									if (pd == null) continue;
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
								st_planes.add(pd);  // adding [0] - all supertile tiles, not clustered by disparity value

								// now try for each of the disparity-separated clusters (only for multi-peak histograms)

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
												plPreferDisparity,												
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
														stSels[m],
														targetV, // double     targetEigen, // target eigenvalue for primary axis (is disparity-dependent, so is non-constant)
														max_outliers, // int        maxRemoved,  // maximal number of tiles to remove (not a constant)
														plMinPoints,  // int        minLeft,     // minimal number of tiles to keep
														plPreferDisparity,
														dl1); // 0); // debugLevel);
												if (pd == null) {
													continue;
												}
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
			final boolean preferDisparity,
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
												preferDisparity,
												(nsTile0 == debug_stile)? 1:0); // int       debugLevel)
										System.out.println("other with other, same weight:");
										TilePlanes.PlaneData other_other_pd = other_pd.mergePlaneToThis(
												other_pd, // PlaneData otherPd,
												1.0,      // double    scale_other,
												true,     // boolean   ignore_weights,
												preferDisparity,
												(nsTile0 == debug_stile)? 1:0); // int       debugLevel)
										System.out.println("other with this, same weight:");
									}
									TilePlanes.PlaneData merged_pd = planes[nsTile0][indx].mergePlaneToThis(
											other_pd, // PlaneData otherPd,
											1.0,      // double    scale_other,
											true,     // boolean   ignore_weights,
											preferDisparity,
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
									preferDisparity,
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

	/**
	 * Measure how merging of two plains degrades individual flatness (smaller - better). For comparing divide by (w1+w2) to make strong
	 * planes score better  			
	 * @param L1 smallest eigenvalue of the first plane
	 * @param L2 smallest eigenvalue of the second plane
	 * @param L  smallest eigenvalue of the merged plane
	 * @param w1 weight of the first plane
	 * @param w2 weight of the second plane
	 * @return degrading by merging measure. 0 if both are co-planar, is supposed to be positive. very "bad" planes do produce negative results - 
	 * not yet clear why (related to non-linear coordinate transformation?) 
	 */

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
		///				double rquality = (L - Lav)*(w1+w2) /(Lav*w1*w2); // ==wors/(w1+w2) to amplify stronger planes
		double rquality =  (L - Lav)*(w1+w2)*(w1+w2) /(Lav*w1*w2); 
		return rquality;
	}

	public void matchPlanes(
			final boolean    preferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
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
		// Select best symmetrical match, consider only N, NE, E, SE - later opposite ones will be copied
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
								this_plane.initMergedValue();
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

												this_plane.initMergedValue(dir,other_planes.length); // filled with NaN
												for (int np = 0; np < other_planes.length; np ++){
													TilePlanes.PlaneData other_plane = this_plane.getPlaneToThis(
															other_planes[np],
															dl); // debugLevel);
													if (other_plane !=null) { // now always, but may add later
														TilePlanes.PlaneData merged_pd = this_plane.mergePlaneToThis(
																other_plane, // PlaneData otherPd,
																1.0,         // double    scale_other,
																false,       // boolean   ignore_weights,
																preferDisparity, 
																dl); // int       debugLevel)
														if (merged_pd !=null) { // now always, but may add later
															this_plane.setNeibMatch(dir, np, merged_pd.getValue()); // smallest eigenValue
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
								for (int dir = 4; dir < 8; dir++){ // other half - copy from opposite
									int stx = stx0 + dirsYX[dir][1];
									int sty = sty0 + dirsYX[dir][0];
									if ((sty < stilesY) && (sty > 0) && (stx > 0)) {
										int nsTile = sty * stilesX + stx; // from where to get
										TilePlanes.PlaneData [] other_planes = planes[nsTile];
										if (other_planes !=null) {
											this_plane.initMergedValue(dir,other_planes.length); // filled with NaN
											for (int np = 0; np < other_planes.length; np ++){
												if (other_planes[np] != null) { // && (other_planes[np].getMergedValue(dir-4) != null)) {
													double [] nm = other_planes[np].getMergedValue(dir-4);
													if (nm != null) {
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

	// Oscolete?			
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
								if (this_plane.getMergedValue() != null) {
									for (int dir = 0; dir < 8; dir++){ //
										double [] neib_worse = this_plane.getMergedValue(dir);
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
	public int fillSquares()
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int nStiles =       stilesX * stilesY;
		int num_added = 0;
		for (int stY = 0; stY < (stilesY - 1); stY++ ) {
			for (int stX = 0; stX < (stilesX - 1); stX++ ) {
				int nsTile = stY * stilesX + stX;
				TilePlanes.PlaneData [][] planes4 = {
						planes[nsTile],
						planes[nsTile + 1],
						planes[nsTile + stilesX],
						planes[nsTile + stilesX + 1]};
				if ((planes4[0] != null) && (planes4[1] != null) && (planes4[2] != null) && (planes4[3] != null)){
					for (int np = 0; np < planes4[0].length; np++){
						if (planes4[0][np] == null) continue;
						int [] neibs0 = planes4[0][np].getNeibBest();
						if ((neibs0[2] < 0) || (neibs0[4] < 0)) continue;
						if (planes4[1][neibs0[2]] == null) continue;
						if (planes4[2][neibs0[4]] == null) continue;
						int [] neibs1 = planes4[1][neibs0[2]].getNeibBest();
						int [] neibs2 = planes4[2][neibs0[4]].getNeibBest();
						if ((neibs1[4] < 0) || (neibs2[2] < 0) || (neibs1[4] != neibs2[2])) continue;
						int [] neibs3 = planes4[3][neibs1[4]].getNeibBest();
						// got full square, does it already have diagonals
						if (neibs0[3] < 0) {
							neibs0[3] = neibs1[4];
							neibs3[7] = np;
							num_added++;
						}
						if (neibs1[5] < 0){
							neibs1[5] = neibs0[4];
							neibs2[1] = neibs0[2];
							num_added++;
						}
					}
				}
			}
		}
		return num_added;
	}

	public int cutCorners()
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
//		final int nStiles =       stilesX * stilesY;
		int num_added = 0;
		for (int stY = 0; stY < (stilesY - 1); stY++ ) {
			for (int stX = 0; stX < (stilesX - 1); stX++ ) {
				int nsTile = stY * stilesX + stX;
				TilePlanes.PlaneData [][] planes4 = {
						planes[nsTile],
						planes[nsTile + 1],
						planes[nsTile + stilesX],
						planes[nsTile + stilesX + 1]};
				if (planes4[0] != null) { // && (planes4[1] != null) && (planes4[2] != null)){
					for (int np = 0; np < planes4[0].length; np++){
						if (planes4[0][np] != null) {
							int [] neibs0 = planes4[0][np].getNeibBest();

							if ((neibs0[2] >= 0) && (neibs0[4] < 0)){
								if (planes4[1][neibs0[2]] != null) {
									int [] neibs1 = planes4[1][neibs0[2]].getNeibBest();
									if (neibs1[5] >= 0) {
										neibs0[4] = neibs1[5];
										int [] neibs2 = planes4[2][neibs1[5]].getNeibBest();
										neibs2[0] = np;  //?
										num_added++;
									}
								}
							}

							if ((neibs0[2] < 0) && (neibs0[4] >= 0)){
								if (planes4[2][neibs0[4]] != null) {
									int [] neibs2 = planes4[2][neibs0[4]].getNeibBest();
									if (neibs2[1] >= 0) {
										neibs0[2] = neibs2[1];
										int [] neibs1 = planes4[1][neibs2[1]].getNeibBest();
										neibs1[6] = np;
										num_added++;
									}
								}
							}
							if ((neibs0[2] >= 0) && (neibs0[3] >= 0)){
								int [] neibs1 = planes4[1][neibs0[2]].getNeibBest();
								int [] neibs3 = planes4[3][neibs0[3]].getNeibBest();
								if (neibs1[4] < 0) {
									neibs1[4] = neibs0[3];
									neibs3[0] = neibs0[2];  //?
									num_added++;
								}
							}
							if ((neibs0[3] >= 0) && (neibs0[4] >= 0)){
								int [] neibs3 = planes4[3][neibs0[3]].getNeibBest();
								int [] neibs2 = planes4[2][neibs0[4]].getNeibBest();
								if (neibs2[2] < 0) {
									neibs2[2] = neibs0[3];
									neibs3[6] = neibs0[4];
									num_added++;
								}
							}							
						}
					}
				}
				if (planes4[3] != null) { // && (planes4[1] != null) && (planes4[2] != null)){
					for (int np = 0; np < planes4[3].length; np++){
						if (planes4[3][np] != null){
							int [] neibs3 = planes4[3][np].getNeibBest();
							
							if ((neibs3[0] >= 0) && (neibs3[7] >= 0)){
								int [] neibs0 = planes4[0][neibs3[7]].getNeibBest();
								int [] neibs1 = planes4[1][neibs3[0]].getNeibBest();
								if (neibs0[2] < 0) {
									neibs0[2] = neibs3[0];
									neibs1[6] = neibs3[7];
									num_added++;
								}
							}
							
							
							if ((neibs3[6] >= 0) && (neibs3[7] >= 0)){
								int [] neibs0 = planes4[0][neibs3[7]].getNeibBest();
								int [] neibs2 = planes4[2][neibs3[6]].getNeibBest();
								if (neibs0[4] < 0) {
									neibs0[4] = neibs3[6];
									neibs2[0] = neibs3[7];  //?
									num_added++;
								}
							}
							

							if ((neibs3[7] >= 0) && (neibs3[6] < 0)){
								int [] neibs0 = planes4[0][neibs3[7]].getNeibBest();
								
								if (neibs0[4] >= 0) {
									int [] neibs2 = planes4[2][neibs0[4]].getNeibBest();
									neibs3[6] = neibs0[4];
									neibs2[2] = np;  //?
									num_added++;
								}
							}
							
							if ((neibs3[6] >= 0) && (neibs3[0] < 0)){
								int [] neibs2 = planes4[2][neibs3[6]].getNeibBest();
								
								if (neibs2[1] >= 0) {
									int [] neibs1 = planes4[1][neibs2[1]].getNeibBest();
									neibs3[0] = neibs2[1]; //?
									neibs1[4] = np;
									num_added++;
								}
							}
						}
					}					
				}
			}
		}
		return num_added;

	}
	
	/**
	 * Find mutual links between multi-layer planes for supertiles. requires that for each plane there are calculated smalles eigenvalues
	 * for merging with each plane for each of 8 neighbors 
	 * @param rquality maximal degradation by merging (does not depend on the total weight)
	 * @param maxEigen maximal eigenvalue of each of the merged planes
	 * @param minWeight minimal weight of each of the planes
	 * @param debugLevel debug level
	 * @param dbg_X debug supertile X coordinate
	 * @param dbg_Y debug supertile Y coordinate
	 */
	public void selectNeighborPlanesMutual(
			final double rquality,
			final double weakWorsening,
			final double dispNorm,
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
					for (int nsTile0 = ai.getAndIncrement(); nsTile0 < nStiles; nsTile0 = ai.getAndIncrement()) {
						if ( planes[nsTile0] != null) {
							for (int np0 = 0; np0 < planes[nsTile0].length; np0++){ // nu
								TilePlanes.PlaneData this_plane = planes[nsTile0][np0];
								this_plane.initNeibBest();
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);

		ai.set(0);				
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
									if (planes[nsTile0][np0].getMergedValue(dir) != null){
										int l = planes[nsTile0][np0].getMergedValue(dir).length;
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
											double [] merge_ev = planes[nsTile0][np0].getMergedValue(dir);
											if (!this_matched[np0] &&
													(merge_ev != null) &&
													((maxEigen == 0.0) ||
															(planes[nsTile0][np0].getValue() < corrMaxEigen(
																	maxEigen,
																	dispNorm,
																	planes[nsTile0][np0]))) &&
													(planes[nsTile0][np0].getWeight() > minWeight)) {
												for (int np = np_min; np < merge_ev.length; np++){
													if (!other_matched[np] &&
															!Double.isNaN(merge_ev[np]) &&
															((maxEigen == 0.0) ||
																	(planes[nsTile][np].getValue() < corrMaxEigen(
																			maxEigen,
																			dispNorm,
																			planes[nsTile][np]))) &&
															(planes[nsTile][np].getWeight() > minWeight)) {
														double w1 = planes[nsTile0][np0].getWeight();
														double w2 = planes[nsTile][np].getWeight();
														double this_rq = mergeRQuality(
																planes[nsTile0][np0].getValue(), // double L1,
																planes[nsTile][np].getValue(), // double L2,
																merge_ev[np], // double L,
																w1, // double w1,
																w2); // double w2)
														double this_rq_norm = this_rq;
														if ((w1 + w2) < weakWorsening) this_rq_norm *= (w1 + w2) / weakWorsening; // forgive more for weak planes 
														if (this_rq_norm <= rquality) {   
															this_rq /= (w1 + w2); // for comparision reduce this value for stronger planes 
															if (Double.isNaN(best_rqual) || (this_rq < best_rqual)){ // OK if Double.isNaN(this_rq[np])
																best_rqual = this_rq;
																best_pair[0]= np0;
																best_pair[1]= np;
															}
														}
														if ((merge_ev[np] < 0.4) && (w1 > 1.0)  && (w2 > 1.0) ){
															System.out.println("nsTile0="+nsTile0+":"+np0+", nsTile="+nsTile+":"+np+", this_rq="+this_rq+
																	" w1="+w1+" w2="+w2+
																	" L1="+planes[nsTile0][np0].getValue()+" L2="+planes[nsTile][np].getValue()+" L="+merge_ev[np]);
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
						/*
								if (disp > dispNorm) {
									eigVal *= dispNorm / disp;
								}
								if (eigVal < maxEigen) num_sat ++;
						 */
						if (eigVal < corrMaxEigen(
								maxEigen,
								dispNorm,
								disp)) num_sat ++;
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
							double disp = planes[nsTile][np].getZxy()[0];
							/*
									if (disp > dispNorm) {
										eigVal *= dispNorm / disp;
									}
							 */
							if (eigVal < corrMaxEigen(
									maxEigen,
									dispNorm,
									disp)) {
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

	public double [][] getShowShells(
			TilePlanes.PlaneData [][] planes,
			int [][] shells,
			int max_shells,
			boolean fuse,
			boolean show_connections,
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
		final int [] neib_offsets={
				-stilesX,
				-stilesX + 1,
				1,
				stilesX +  1,
				stilesX,
				stilesX -  1,
				-1,
				-stilesX - 1};
		//				final int debug_stile = 17 * stilesX + 27;
		final int debug_stile = -1; // 18 * stilesX + 30;

		int num_shells = 0;
		for (int nsTile = 0; nsTile < shells.length; nsTile++) if (shells[nsTile] != null){
			for (int np = 0; np < shells[nsTile].length; np++){
				if (shells[nsTile][np] > num_shells) {
					num_shells = shells[nsTile][np]; // starts with 1, not 0 
				}
			}
		}
		if (num_shells > max_shells){
			num_shells = max_shells;
		}
		double [][] data = new double [num_shells][width*height];
		for (int nShell = 0; nShell < num_shells; nShell++) {
			for (int sty = 0; sty < stilesY; sty++){
				for (int stx = 0; stx < stilesX; stx++){
					int nsTile = sty * stilesX + stx;
					boolean debug = (nsTile == debug_stile) ; // (sty == 18) && (stx == 28);
					String [] dbg_titles = {"Center","N","NE","E","SE","S","SW","W","NW"};
					double [][] dbg_img3 = null;
					double [][] dbg_img = null;
					double [][] dbg_this = null;
					double [][] dbg_other = null;
					if (debug) {
						System.out.println("getShowShells():2, tile = "+nsTile);
					}
					double [] plane = nan_plane;
					int indx = (sty * superTileSize) * width + (stx * superTileSize);
					if (debug) {
						System.out.println("getShowShells():"+((planes[nsTile] != null)? planes[nsTile].length : "null"));
					}
					if (planes[nsTile] != null) {
						for (int np = 0; np < planes[nsTile].length; np++){
							if (shells[nsTile][np] == (nShell + 1)) {
								plane = planes[nsTile][np].getPlaneDisparity(use_NaN);
								int [] neib_best = planes[nsTile][np].getNeibBest();
								if (fuse && (neib_best != null)){
									if (debug) {
										dbg_img3 = new double [dbg_titles.length][];
										dbg_img = new double [dbg_titles.length][];
										dbg_this = new double [dbg_titles.length][];
										dbg_other = new double [dbg_titles.length][];
										dbg_img3[0] = planes[nsTile][np].getTriplePlaneDisparity();
										dbg_img[0] = planes[nsTile][np].getPlaneDisparity(false);
									}


									double [] weight_this =  new double [superTileSize*superTileSize];
									double [] weight_other = new double [superTileSize*superTileSize];
									double [] plane_combo =  new double [superTileSize*superTileSize];
									for (int dir = 0; dir < 8; dir++){
										if (neib_best[dir] >=0){
											int nsTileOther = nsTile + neib_offsets[dir];
											int npOther = neib_best[dir];
											int dirOther = (dir + 4) % 8;
											double [] plane_other =	planes[nsTileOther][npOther].getTriplePlaneDisparity(dirOther);
											double [] fuse_this =  getFuseThis(dir);
											double [] fuse_other = getFuseOther(dir);
											if (debug) {
												dbg_img3[dir + 1] =   planes[nsTileOther][npOther].getTriplePlaneDisparity();
												dbg_img[dir + 1] =    planes[nsTileOther][npOther].getTriplePlaneDisparity(dirOther);
												dbg_this[dir + 1] =   fuse_this.clone();
												dbg_other[dir + 1] =  fuse_other.clone();
											}
											for (int i = 0; i < weight_this.length; i++){
												weight_this[i] +=  fuse_this[i];
												weight_other[i] += fuse_other[i];
												plane_combo[i] +=  fuse_other[i] * plane_other[i];
											}
										}
									}
									if (debug) {
										showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
										sdfa_instance.showArrays(dbg_img, superTileSize, superTileSize, true, "plane_"+nsTile,dbg_titles);
										sdfa_instance.showArrays(dbg_img3, 3*superTileSize, 3*superTileSize, true, "plane3_"+nsTile,dbg_titles);
										sdfa_instance.showArrays(dbg_this, superTileSize, superTileSize, true, "weights_this_"+nsTile,dbg_titles);
										sdfa_instance.showArrays(dbg_other, superTileSize, superTileSize, true, "weights_other_"+nsTile,dbg_titles);
									}
									double [] plane_other = planes[nsTile][np].getPlaneDisparity(false);
									for (int i = 0; i < weight_this.length; i++){
										plane_combo[i] += weight_this[i] * plane_other[i];
										double w = weight_this[i] + weight_other[i];
										if (w > 0){
											plane_combo[i] /= w;
										} else {
											plane_combo[i] = plane_other[i];
										}
									}											
									if (use_NaN){
										for (int i = 0; i < plane_combo.length; i++){
											if (Double.isNaN(plane[i])){
												plane_combo[i] = Double.NaN;
											}
										}
									}
									plane = plane_combo;
								}

								if (show_connections) {
									if (neib_best != null){
										for (int dir = 0; dir < neib_best.length; dir ++) if (neib_best[dir] >= 0){
											// draw a line
											//int center_index = indx + centerIndex;
											for (int i = 0; i < (superTileSize / 2); i++){
												plane[centerIndex + i * dirs[dir]] = ((i & 1) == 0) ? arrow_dark : arrow_white; 
											}
										}
									}
								}
								for (int y = 0; y < superTileSize; y++) {
									System.arraycopy(plane, superTileSize * y, data[nShell], indx + width * y, superTileSize);
								}
							}
						}
					}
				}						
			}
		}
		return data;
	}


	public TilePlanes.PlaneData[][] copyPlanes(
			TilePlanes.PlaneData[][] src_planes)
	{
		TilePlanes.PlaneData[][] dst_planes = new TilePlanes.PlaneData[src_planes.length][];
		return copyPlanes(src_planes, dst_planes);
	}

	public TilePlanes.PlaneData[][] copyPlanes(
			final TilePlanes.PlaneData[][] src_planes,
			final TilePlanes.PlaneData[][] dst_planes)
	{
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nsTile = ai.getAndIncrement(); nsTile < src_planes.length; nsTile = ai.getAndIncrement()) {
						if (src_planes[nsTile] != null){
							dst_planes[nsTile] = new TilePlanes.PlaneData[src_planes[nsTile].length];
							for (int np = 0; np < src_planes[nsTile].length; np++){
								if (src_planes[nsTile][np] != null){
									dst_planes[nsTile][np] = src_planes[nsTile][np].clone(); 
								} else {
									dst_planes[nsTile][np] = null;
								}
							}
						} else {
							dst_planes[nsTile] = null;
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return dst_planes;
	}

	public TilePlanes.PlaneData[][] planesSmooth(
			final double      meas_pull,//  relative pull of the original (measured) plane with respect to the average of the neighbors
			final int         num_passes,
			final double      maxDiff, // maximal change in any of the disparity values
			final boolean     preferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
			final int         debugLevel,
			final int         dbg_X,
			final int         dbg_Y)
	{
		if (this.planes_mod == null){
			this.planes_mod =copyPlanes(this.planes); // make always (for now) *********************
		}
		for (int pass = 0; pass < num_passes; pass++){
			double diff = planesSmoothStep(
					meas_pull,       //  relative pull of the original (measured) plane with respect to the average of the neighbors
					this.planes,     // final TilePlanes.PlaneData[][] measured_planes,
					this.planes_mod, // final TilePlanes.PlaneData[][] mod_planes,
					true,            // final boolean calc_diff,
					preferDisparity, 
					(pass < 10)? debugLevel: 0,
							dbg_X,
							dbg_Y);
			if (diff < maxDiff){
				if (debugLevel > -1){
					System.out.println("planesSmooth(): pass:"+pass+" (of "+num_passes+"), rms = "+diff+" < "+maxDiff);
					break;
				}
			}
			if (debugLevel > -1){
				System.out.println("planesSmooth() -  pass:"+pass+" (of "+num_passes+"), rms = "+diff);
			}

		}
		return this.planes_mod;
	}			

	public double planesSmoothStep(
			final double meas_pull,//  relative pull of the original (measured) plane with respect to the average of the neighbors
			final TilePlanes.PlaneData[][] measured_planes,
			final TilePlanes.PlaneData[][] mod_planes,
			final boolean calc_diff,
			final boolean preferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
			final int debugLevel,
			final int dbg_X,
			final int dbg_Y)
	{
		final int [][] dirsYX = {{-1, 0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1}};
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		final TilePlanes.PlaneData[][] new_planes = copyPlanes(mod_planes);
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final int numThreads = threads.length;
		final double [] rslt_diffs = calc_diff ? new double [numThreads] : null; // all 0;
		final AtomicInteger ai_numThread = new AtomicInteger(0);
		final AtomicInteger ai = new AtomicInteger(0);
		final String [] titles= {
				"orig",     // 0
				"measured", // 1
				"n0","n1","n2","n3","n4","n5","n6","n7", // 2 .. 9
				"m0","m1","m2","m3","m4","m5","m6","m7","mm", // 10..18
		"diff"}; // 19

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					double [][] dbg_img=null;
					int numThread = ai_numThread.getAndIncrement(); // unique number of thread to write to rslt_diffs[numThread]
					for (int nsTile0 = ai.getAndIncrement(); nsTile0 < mod_planes.length; nsTile0 = ai.getAndIncrement()) {
						int sty0 = nsTile0 / stilesX;  
						int stx0 = nsTile0 % stilesX;
						int dl = ((debugLevel > -1) && (nsTile0 == debug_stile)) ? 1:0;
						if ( new_planes[nsTile0] != null) {
							if (dl > 0){
								System.out.println("planesSmoothStep nsTile0="+nsTile0);
								dbg_img = new double [titles.length][];
							}
							int np0_min = (new_planes[nsTile0].length > 1) ? 1:0; // Modify if overall plane will be removed
							for (int np0 = np0_min; np0 < new_planes[nsTile0].length; np0 ++){
								TilePlanes.PlaneData this_new_plane = new_planes[nsTile0][np0];
								if (dl > 0) dbg_img[ 0] = this_new_plane.getPlaneDisparity(false);
								if (dl > 0) dbg_img[ 1] = measured_planes[nsTile0][np0].getPlaneDisparity(false);

								int [] neibs = this_new_plane.getNeibBest();
								this_new_plane.setWeight(0.0);
								for (int dir = 0; dir < neibs.length; dir++){
									if (neibs[dir] >= 0) {
										int stx = stx0 + dirsYX[dir][1];
										int sty = sty0 + dirsYX[dir][0];
										int nsTile = sty * stilesX + stx; // from where to get
										TilePlanes.PlaneData other_plane = this_new_plane.getPlaneToThis(
												mod_planes[nsTile][neibs[dir]], // neighbor, previous value
												dl); // debugLevel);
										if (dl > 0) dbg_img[ 2 + dir] = other_plane.getPlaneDisparity(false);
										if (other_plane != null){
											if (this_new_plane.getWeight() > 0.0){
												this_new_plane = this_new_plane.mergePlaneToThis(
														other_plane, // PlaneData otherPd,
														1.0,         // double    scale_other,
														false,       // boolean   ignore_weights,
														preferDisparity,
														dl); // int       debugLevel)
											} else {
												this_new_plane = other_plane; 
											}
											new_planes[nsTile0][np0] = this_new_plane;
											if (dl > 0) dbg_img[10 + dir] = this_new_plane.getPlaneDisparity(false);
										}

									}
								}
								if (    (meas_pull > 0.0) &&
										(measured_planes != null) &&
										(measured_planes[nsTile0] != null) &&
										(measured_planes[nsTile0][np0] != null)){ // merge with "measured"
									if (this_new_plane.getWeight() > 0.0){
										this_new_plane = this_new_plane.mergePlaneToThis(
												measured_planes[nsTile0][np0], // PlaneData otherPd,
												meas_pull,                     // double    scale_other,
												false,                         // boolean   ignore_weights,
												preferDisparity,
												dl);                           // int       debugLevel)
									} else {
										this_new_plane = measured_planes[nsTile0][np0].clone();
									}
									new_planes[nsTile0][np0] = this_new_plane;
									if (dl > 0) dbg_img[18] =  this_new_plane.getPlaneDisparity(false);
								}
								// calculate largest disparity difference between old and new plane
								if (rslt_diffs != null){
									// get plane for both old and new, calc rms of diff
									double [] oldPlane = mod_planes[nsTile0][np0].getPlaneDisparity(
											false); // use_NaN)
									double [] newPlane = new_planes[nsTile0][np0].getPlaneDisparity(
											false); // use_NaN)
									double s = 0.0;
									for (int i = 0; i < oldPlane.length; i++){
										double d = newPlane[i] - oldPlane[i];
										s+= d*d;
									}
									s= Math.sqrt(s/oldPlane.length);
									if (s > rslt_diffs[numThread]){
										rslt_diffs[numThread] = s;
									}
									if (dl > 0) {
										dbg_img[19] =  new double[oldPlane.length];
										for (int i = 0; i < oldPlane.length; i++){
											dbg_img[19][i] = newPlane[i] - oldPlane[i];
										}
									}
								}
								if ((dl > 0) && (debugLevel >1)){
									showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
									sdfa_instance.showArrays(dbg_img, superTileSize, superTileSize, true, "smooth_step_x"+stx0+"_y"+sty0, titles);
								}
							}
						}								
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		copyPlanes (new_planes, mod_planes); // copy back
		if (rslt_diffs == null){
			return Double.NaN;
		}
		double diff = 0.0;
		for (int i = 0; (i < numThreads) ; i++) if (diff < rslt_diffs[i]) diff = rslt_diffs[i];
		return diff; // return maximal difference
	}

	public double [][] planesGetDiff(
			final TilePlanes.PlaneData[][] measured_planes,
			final TilePlanes.PlaneData[][] mod_planes,
			final int debugLevel,
			final int dbg_X,
			final int dbg_Y)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		final TilePlanes.PlaneData[][] new_planes = copyPlanes(mod_planes);
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final double [][] diffs = null;

		return diffs;
	}
	
	/**
	 * Create candidate planes to break a single plane in 2 by splitting consecutive connected
	 * neighbors in 2 groups that make the smallest weighted sum of the eigenvalues
	 * Groups may be optionally supplemented by the center supertile 
	 * @param center_planes [per supertile][per plane] array of plane objects to use as
	 *        the source of connections and optionally other data for the center tiles
	 * @param neib_planes  [per supertile][per plane] array of plane objects to use for
	 *        neighbor data. May be the same as center_planes or different
	 * @param center_pull - merge with center plane when calculating half-planes with this
	 *  relative weight: 0.0 - only neighbors, 1.0 - same weight of the center as each neighbor.
	 * @param min_neibs - minimal number of connected neighbors to work with (>=2)       
	 * @param preferDisparity - the first eigenvalue/vector is the most disparity-like
	 *                          (false - smallest eigenvalue)
	 * @param debugLevel debug level
	 * @param dbg_X supertile horizontal index to show debug information
	 * @param dbg_Y supertile vertical index to show debug information
	 * @return a pair of plane objects for each [supertile][plane][3] and a combination of both
	 *         (to compare eigenvalues)
	 */
	public TilePlanes.PlaneData[][][] breakPlanesToPairs(
			final TilePlanes.PlaneData[][] center_planes, // measured_planes,
			final TilePlanes.PlaneData[][] neib_planes,   //mod_planes,
			final double                   center_pull,
			final int                      min_neibs, // 2
			final boolean                  preferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
			final int                      debugLevel,
			final int                      dbg_X,
			final int                      dbg_Y)
	{
		final int [][] dirsYX = {{-1, 0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1}};
		final int tilesX =        tileProcessor.getTilesX();
//		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
//		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		final TilePlanes.PlaneData[][][] rslt_planes = new TilePlanes.PlaneData[center_planes.length][][];
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final String [] titles= {"center", "first", "second", "all"};

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					double [][] dbg_img=null;
					for (int nsTile0 = ai.getAndIncrement(); nsTile0 < center_planes.length; nsTile0 = ai.getAndIncrement()) {
						int sty0 = nsTile0 / stilesX;  
						int stx0 = nsTile0 % stilesX;
						int dl = ((debugLevel > -1) && (nsTile0 == debug_stile)) ? 1:0;
						if ( center_planes[nsTile0] != null) {
							rslt_planes[nsTile0] = new TilePlanes.PlaneData[center_planes[nsTile0].length][];
							if (dl > 0){
								System.out.println("breakPlanesToPairs nsTile0="+nsTile0);
								dbg_img = new double [titles.length][];
							}
							int np0_min = (center_planes[nsTile0].length > 1) ? 1:0; // Modify if overall plane will be removed
							for (int np0 = np0_min; np0 < center_planes[nsTile0].length; np0 ++){
								TilePlanes.PlaneData center_plane = center_planes[nsTile0][np0];
								if (dl > 0) dbg_img[ 0] = center_plane.getPlaneDisparity(false);
								int [] neibs = center_plane.getNeibBest();
								int num_neibs = 0;
								for (int i = 0; i < neibs.length; i++) if (neibs[i] >= 0) num_neibs++;
								if (num_neibs >= min_neibs) {
									// create all pairs to test
									int num_splits = num_neibs * (num_neibs - 1) / 2;
									double [] split_quality = new double [num_splits]; // weighted sum of eigenvalues of merged
									int [][][] neibs12 = new int [2][num_splits][];
									TilePlanes.PlaneData[][] plane_triads = new TilePlanes.PlaneData[num_splits][3]; 
									int [] neib_index = new int [num_neibs];
									int nn = 0;
									for (int i = 0; i < neibs.length; i++) if (neibs[i] >= 0) neib_index[nn++] = i;
									nn = 0;
									for (int nStart = 1; nStart < num_neibs; nStart++) {
										for (int nEnd = nStart+1; nEnd <= num_neibs; nEnd++){
											neibs12[0][nn] = neibs.clone();
											neibs12[1][nn] = neibs.clone();
											for (int i = 0; i < neib_index.length; i++) {
												if ((i < nStart) || (i >= nEnd)){
													neibs12[0][nn][neib_index[i]] = -1;
												} else {
													neibs12[1][nn][neib_index[i]] = -1;
												}
											}
											nn++;
										}
									}
									
									for (int nSplit = 0; nSplit < num_splits; nSplit++) {
										// Calculate merged plane for each selection (and other), weighted average and store to split_quality
										for (int nis = 0; nis < 2; nis++) {
											plane_triads[nSplit][nis] = center_planes[nsTile0][np0].clone();
											plane_triads[nSplit][nis].setWeight(center_pull * plane_triads[nSplit][nis].getWeight());
											for (int dir = 0; dir <8; dir++){
												int np = neibs12[nis][nSplit][dir];
												if (np >= 0){
													int stx = stx0 + dirsYX[dir][1];
													int sty = sty0 + dirsYX[dir][0];
													int nsTile = sty * stilesX + stx; // from where to get
													TilePlanes.PlaneData other_plane = plane_triads[nSplit][nis].getPlaneToThis(
															neib_planes[nsTile][np],
															dl); // debugLevel);
													// if (dl > 0) dbg_img[ 2 + dir] = other_plane.getPlaneDisparity(false);
													if (other_plane != null){
														if (plane_triads[nSplit][nis].getWeight() > 0.0){
															plane_triads[nSplit][nis] = plane_triads[nSplit][nis].mergePlaneToThis(
																	other_plane, // PlaneData otherPd,
																	1.0,         // double    scale_other,
																	false,       // boolean   ignore_weights,
																	preferDisparity,
																	dl); // int       debugLevel)
														} else {
															plane_triads[nSplit][nis] = other_plane; 
														}
														//if (dl > 0) dbg_img[10 + dir] = this_new_plane.getPlaneDisparity(false);
													} else {
														plane_triads[nSplit][nis] = null;
														break;
													}
												}
											}
										} // for (int nis = 0; nis < 2; nis++) {
										if ((plane_triads[nSplit][0] != null) && (plane_triads[nSplit][1] != null)){
											double w1 = plane_triads[nSplit][0].getWeight();
											double w2 = plane_triads[nSplit][1].getWeight();
											split_quality[nSplit] = (w1 * plane_triads[nSplit][0].getValue() + w2 * plane_triads[nSplit][1].getValue())/
													(w1 + w2);
										} else {
											split_quality[nSplit] = Double.NaN;
										}
									}
									// find minimum in split_quality and generate a pair of plane objects, setting neighbors for each
									// later use these plane pairs to assign tiles to each and generate a new eigenvalues/vectors
									// TODO: How to handle a pair? Any special treatment (like fusing?), 
									// if the plane intersection line is inside supertile - use min/max for snapping
									int best_index = -1;
									for (int nSplit = 0; nSplit < num_splits; nSplit++) {
										if (!Double.isNaN(split_quality[nSplit]) && ((best_index < 0) || (split_quality[nSplit] < split_quality[best_index]))){
											best_index = nSplit;
										}
									}
									if (best_index >= 0) {
										plane_triads[best_index][2] = plane_triads[best_index][0].clone().mergePlaneToThis(
												plane_triads[best_index][1], // PlaneData otherPd,
													1.0,                    // double    scale_other,
													false,                  // boolean   ignore_weights,
													preferDisparity,
													dl);										
										for (int nis = 0; nis < 2; nis++) {
											plane_triads[best_index][nis].setNeibBest(neibs12[nis][best_index]);
										}
										rslt_planes[nsTile0][np0] = plane_triads[best_index];
										if (dl > 0) dbg_img[ 1] = plane_triads[best_index][0].getPlaneDisparity(false);
										if (dl > 0) dbg_img[ 2] = plane_triads[best_index][1].getPlaneDisparity(false);
										if (dl > 0) dbg_img[ 3] = plane_triads[best_index][2].getPlaneDisparity(false);
										
										// merged_eig_val are not set - will have to be recalculated when updating after changing tile selections
										// to make each plane
										// TODO save gain from splitting to planes 
										
									}
									
									if ((dl > 0) && (debugLevel > -1)){
										showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
										sdfa_instance.showArrays(dbg_img, superTileSize, superTileSize, true, "breakPlanesToPairs"+stx0+"_y"+sty0+"_p"+np0, titles);
									}
								}
							}
						}								
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return rslt_planes;
	}




	public boolean [][] selectPlanes(
			final double dispNorm,
			final double maxEigen, // maximal eigenvalue of planes to consider
			final double minWeight, // minimal pain weight to consider
			TilePlanes.PlaneData [][] planes)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		boolean [][] selection = new boolean [planes.length][];
		for (int sty = 0; sty < stilesY; sty++){
			for (int stx = 0; stx < stilesX; stx++){
				int nsTile = sty * stilesX + stx;
				if (planes[nsTile] != null) {
					selection[nsTile] = new boolean [planes[nsTile].length];
					for (int np = 0; np < planes[nsTile].length; np++){
						if (planes[nsTile][np].getWeight() >= minWeight){
							//&& (planes[nsTile][np].getValue() < maxEigen)){
							double eigVal = planes[nsTile][np].getValue();
							double disp = planes[nsTile][np].getZxy()[0];
							selection[nsTile][np] = (eigVal < corrMaxEigen(
									maxEigen,
									dispNorm,
									disp));
						}
					}
				}
			}
		}
		return selection;
	}

	public int [][] createSupertileShells(
			boolean[][]                     selection, // may be null
			boolean                         use_all, // use plane 0 even if there are more than 1
			boolean                         keep_orphans, // single-cell shells
			double                          orphan_strength, // minimal strength no keep orphans (single-tile planes)
			final TilePlanes.PlaneData [][] planes,
			int                             debugLevel)
	{
		class TPW {
			int t;
			int p;
			double w;
			TPW (int t, int p){
				this.t= t;
				this.p= p;
				this.w = planes[t][p].getWeight();
			}
		}
		final int tilesX =        tileProcessor.getTilesX();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
		final int [] dirs= {-stilesX, -stilesX+1, 1, stilesX+1, stilesX, stilesX-1, -1, -stilesX-1};
		boolean [][] field = new boolean [planes.length][];
		for (int nsTile = 0; nsTile < planes.length; nsTile++){
			if (nsTile==724){
				System.out.println("createSupertileShells()1 nsTile="+nsTile);
			}
			if (planes[nsTile] != null){
				field[nsTile] = new boolean [planes[nsTile].length];
				if (selection != null) {
					for (int i = 0; i < selection[nsTile].length; i++){
						if (!selection[nsTile][i]) field[nsTile][i] = true; // used/prohibited
					}
				}
				if (!use_all && (planes[nsTile].length > 1)){
					field[nsTile][0] = true; // used/prohibited
				}
			}
		}
		ArrayList<ArrayList<TPW>> listShells = new ArrayList<ArrayList<TPW>>();
		for (int nsTile = 0; nsTile < field.length; nsTile++){
			if (nsTile==724){
				System.out.println("createSupertileShells()2 nsTile="+nsTile);
			}

			if (field[nsTile] != null){
				for (int np = 0; np < field[nsTile].length; np++){
					if (!field[nsTile][np]){
						// see if 
						ArrayList<TPW> listWave = new ArrayList<TPW>();
						int read_p = 0;
						//								int num_cluster = listShells.size()+1;
						TPW sp = new TPW(nsTile,np);
						field[sp.t][sp.p] = true; // just 1 would be OK
						listWave.add(sp);
						while (read_p < listWave.size()) { // keeping all list, not just front
							TPW tail = listWave.get(read_p++);
							int [] neibs = planes[tail.t][tail.p].getNeibBest();
							for (int dir = 0; dir < neibs.length; dir++){
								if (neibs[dir] >= 0) {
									sp = new TPW(
											tail.t + dirs[dir],
											neibs[dir]);
									if (!field[sp.t][sp.p] && ((selection == null) || (selection[sp.t][sp.p]))){
										field[sp.t][sp.p] = true;
										listWave.add(sp);
									}
								}
							}

						}
						if (listWave.size() > 1) {
							listShells.add(listWave);
						} else if (keep_orphans && (listWave.size() == 1)){ // so far always >= 1
							TPW tail = listWave.get(0);
							if (planes[tail.t][tail.p].getWeight() >= orphan_strength){
								listShells.add(listWave);
							}
						}
						if (debugLevel > 1){
							System.out.println("CreateSupertileShells(): "+(listShells.size())+":"+listWave.size());
						}
					}
				}
			}
		}
		// sort list
		Collections.sort(listShells, new Comparator<ArrayList<TPW>>() {
			@Override
			public int compare(ArrayList<TPW> lhs, ArrayList<TPW> rhs) {
				// -1 - less than, 1 - greater than, 0 - equal, all inverted for descending
				int sgn = lhs.size() > rhs.size() ? -1 : (lhs.size() < rhs.size() ) ? 1 : 0;
				if (sgn != 0) return sgn;
				double wl = 0.0, wr = 0.0;
				for (TPW p:lhs) wl += p.w;
				for (TPW p:rhs) wr += p.w;
				return wl > wr ? -1 : (wl < wr ) ? 1 : 0;
			}
		});


		int [][] shells = new int [planes.length][];
		for (int nsTile = 0; nsTile < planes.length; nsTile++){
			if (planes[nsTile] != null){
				shells[nsTile] = new int [planes[nsTile].length];
			}
		}
		for (int num_shell = 0; num_shell < listShells.size(); num_shell++){
			for (TPW p:listShells.get(num_shell)){
				shells[p.t][p.p] = num_shell+1;
			}
			if (debugLevel > 0){
				double w = 0;
				for (TPW p:listShells.get(num_shell)){
					w += p.w;
				}
				System.out.println("CreateSupertileShells(): "+(num_shell + 1)+":"+listShells.get(num_shell).size()+", total strength = "+w);
			}

		}
		return shells;
	}

	public double normDispADiff(
			double disp1,
			double disp2,
			double disp_norm
			)
	{
		double disp_avg = 0.5* (disp1 + disp2);
		double adiff = Math.abs(disp1 - disp2);
		return (disp_avg <= disp_norm) ? adiff : (adiff * disp_norm / disp_avg);
	}

	public double [] snapDisparity(
			final double []  disparity,
			final double []  strength,
			final boolean [] selection,       // can be null
			final double     dispNorm, // plDispNorm
			final double     snapDispAny,    //    =  .2,   // Maximal (scaled by plDispNorm) disparity difference to snap to plane at any strength
			final double     snapStrengthAny,    //    =  .2,   //  Maximal strength to fit any distance (if does not fit otherwise - treat as zero strength
			final double     snapNegAny,     //    =  .2,   // Maximal negative disparity difference from the best match
			final double     snapDispMax,    //    =  .5,   // Maximal (scaled by plDispNorm) disparity difference to snap to plane at low strength
			final double     snapDispWeight, //    =  .5,   // Maximal disparity diff. by weight product to snap to plane
			final int        snapZeroMode,   // Zero strength snap mode: 0: no special treatment, 1 - strongest, 2 - farthest
			final int        debugLevel)
	{
		final int tilesX =        tileProcessor.getTilesX();
//		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
//		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final double [][] surfaces = getSurfaces();
		final int [][] shell_map =   getShellMap();
		final TilePlanes.PlaneData [][] planes = getPlanesMod(); // to get strength;
		final double [] snap_disp =  disparity.clone();

		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < disparity.length; nTile = ai.getAndIncrement()) {
						if ((selection == null) || selection[nTile]) {
							int tX = nTile % tilesX; 
							int tY = nTile / tilesX;
							int stX = tX / superTileSize;
							int stY = tY / superTileSize;
							int nsTile = stY * stilesX + stX;
							int nTile_st = (tY * stilesX * superTileSize) + tX; // surfaces maybe (is) larger than tilesX *tilesY if 
							// tilesX, tilesY are not multiples of superTileSize !
							if (nsTile == 724){
//								System.out.println("snapDisparity(): tX:tY = "+tX+":"+tY+", stX:stY="+stX+":"+stY+
//										", nTile="+nTile+", nsTile="+nsTile+", nTile_st = "+nTile_st);
							}

							// is there any surface to snap to?
							boolean no_surf = true;
							if (shell_map[nsTile] != null){
								for (int np = 0; no_surf && (np < shell_map[nsTile].length); np++){
									if (shell_map[nsTile][np] > 0) no_surf = false;
								}
							}
							if (!no_surf) {
								double [] surf_disparity = new double[shell_map[nsTile].length];
								for (int np = 0; np < surf_disparity.length; np++){
									if (shell_map[nsTile][np] > 0) {
										surf_disparity[np] = surfaces[shell_map[nsTile][np] -1][nTile_st];
									} else {
										surf_disparity[np] = Double.NaN;
									}
								}
								int best_indx = -1;
								double best_disp_diff = Double.NaN;
								double this_disp = disparity[nTile];
								double this_strength = strength[nTile];
								if (Double.isNaN(this_disp)){
									this_disp = 0.0;
									this_strength = 0.0;
								}
								for (int np = 0; np < surf_disparity.length; np++){
									if (shell_map[nsTile][np] > 0) {
										double disp_diff = normDispADiff (this_disp,surf_disparity[np], dispNorm);
										if ((best_indx < 0) || (disp_diff < best_disp_diff)) {
											best_disp_diff = disp_diff;
											best_indx =      np;
										}
									}											
								}										
								if ((this_strength > 0.0) || (snapZeroMode == 0)){
									if ((best_disp_diff <= snapDispAny) ||
											((best_disp_diff <= snapDispMax) && (best_disp_diff * this_strength < snapDispWeight )) ||  
											((best_disp_diff <= snapNegAny) && (this_disp <  surf_disparity[best_indx]))  // farther than surface - higher tolerance
											){
										snap_disp[nTile] = surf_disparity[best_indx]; // snapped (default - keep source disparity
									} else if ((this_strength < snapStrengthAny) && (snapZeroMode != 0)){ // treat as zero strength?
										this_strength = 0.0;
									}
								}

								if ((this_strength == 0) && (snapZeroMode > 0)) { 
									if (snapZeroMode == 1) {  // zero strength and non-zero snapZeroMode, consider 1 (strongest shell)
										best_indx = -1;
										double best_strength = Double.NaN;
										for (int np = 0; np < surf_disparity.length; np++){
											if (shell_map[nsTile][np] > 0) {
												double w = planes[nsTile][np].getWeight();
												if ((best_indx < 0) || (w >  best_strength)) {
													best_strength = w;
													best_indx =     np;
												}
											}
										}
										snap_disp[nTile] = surf_disparity[best_indx]; // snapped (default - keep source disparity
									} else { // zero strength and snapZeroMode > 1 (now only 2) - farthest surface (smallest disparity)
										best_indx = -1;
										double smallest_disparity = Double.NaN;
										for (int np = 0; np < surf_disparity.length; np++){
											if (shell_map[nsTile][np] > 0) {
												if ((best_indx < 0) || ( surf_disparity[np] <  smallest_disparity)) {
													smallest_disparity = surf_disparity[np];
													best_indx =     np;
												}
											}
										}
										snap_disp[nTile] = surf_disparity[best_indx]; // snapped (default - keep source disparity
									}
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return snap_disp;
	}

	/**
	 * Get number of used surfaces (some may be empty) 
	 * @param snap_surf surface indices (per-tile), 0 - unmatched
	 * @param selected boolean tile selection or null (use all tiles)
	 * @return maximal number of surface used
	 */
	public int getNumSurf(
			int []     snap_surf,
			boolean [] selected) // or null
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
		
		int num_surf = 0;
		for (int nTile = 0; nTile < snap_surf.length; nTile ++) if (snap_surf[nTile] > 0){
			int tX = nTile % tilesX; 
			int tY = nTile / tilesX;
			int stX = tX / superTileSize;
			int stY = tY / superTileSize;
			int nsTile = stY * stilesX + stX;
			int np = snap_surf[nTile] - 1;
			if ((shell_map[nsTile][np] > num_surf) && ((selected == null) || selected[nTile])) {
				num_surf = shell_map[nsTile][np];
			}
		}
		return num_surf;
	}

	/**
	 * Get number of unmatched tiles
	 * @param snap_surf surface indices (per-tile), 0 - unmatched
	 * @param selected boolean tile selection or null (use all tiles)
	 * @return number of unmatched tiles
	 */
	public int getNumNotMatched(
			int []     snap_surf,
			boolean [] selected) // or null
	{
		int num_notmatched = 0;
		for (int i = 0; i < snap_surf.length; i++){
			if ((snap_surf[i] == 0) && ((selected == null) || selected[i])) {
				num_notmatched ++;
			}
		}
		return num_notmatched;
	}
	
	public int getNeibMask(
			int x,
			int y,
			int width,
			int height)
	{
		final int [] corn_side_neib = { // of +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
				0b00011100,    // top left corner
				0b01111100,    // top middle
				0b01110000,    // top right
				0b00011111,    // middle left
				0b11111111,    // middle
				0b11110001,    // middle right
				0b00000111,    // bottom left
				0b11000111,    // bottom middle
				0b11000001};   // bottom right
		int tileType = getNeibType(
				x,
				y,
				width,
				height);
		return corn_side_neib[tileType];
	}
	
	/**
	 * Returns 8-element array - which supertile neighbor enables this tile neighbor (-1 - none)
	 * @param x horizontal position of the tile in supertile 
	 * @param y vertical position of the tile in supertile
	 * @param width supertile width
	 * @param height supertile height
	 * @return 8-element array
	 */
	
	public int[] getNeibBordDir(
			int x,
			int y,
			int width,
			int height)
	{
		final int [][] dir_dir = { // of +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
			    //N  NE   E  SE   S  SW   W  NW	
				{ 0,  0, -1, -1, -1,  6,  6,  7},    // top left corner
				{ 0,  0, -1, -1, -1, -1, -1,  0},    // top middle
				{ 0,  1,  2,  2, -1, -1, -1, -1},    // top right
				{-1, -1, -1, -1, -1,  6,  6,  6},    // middle left
				{-1, -1, -1, -1, -1, -1, -1, -1},    // middle
				{ 0,  2,  2,  2, -1, -1, -1, -1},    // middle right
				{-1, -1, -1,  4,  4,  5,  6,  6},    // bottom left
				{-1, -1, -1,  4,  4,  4, -1, -1},    // bottom middle
				{-1,  2,  2,  3,  4,  4, -1, -1}};   // bottom right
		int tileType = getNeibType(
				x,
				y,
				width,
				height);
		return dir_dir[tileType];
	}

	public int getNeibType(
			int x,
			int y,
			int width,
			int height)
	{
		int tileType = 4;
		int hm1 = height -1;
		int wm1 = width - 1;
		if (y == 0){
			if (x == 0){
				tileType = 0;  
			} else if (x == wm1) {
				tileType = 2;  
			} else {
				tileType = 1;  
			}
		} else if (y == hm1) {
			if (x == 0){
				tileType = 6;  
			} else if (x == wm1) {
				tileType = 8;  
			} else {
				tileType = 7;  
			}
		} else {
			if (x == 0){
				tileType = 3;  
			} else if (x == wm1) {
				tileType = 5;  
			}
		}
		return tileType;
		
	}
	
/**
 * Build tile connections to neighbors form optional selection, per-tile plane indices,
 * and a surface index.
 * @param selected tile selection (or null)
 * @param snap_surf per tile array of supertile plane index plus 1 or 0 for unmatched tiles
 * @param surf_index number of the surface to extract
 * @return per-tile array of -256 for unselected tiles or 8-bit bit mask of the neighbors (0 - N, 1 - NE, ...) 
 */
	public int [] getNeighbors( // creates neighbors mask from bitmask
			final boolean [] selected, // or null
			final int []     snap_surf, // use this size - it matches image, not supertiles
			final int        surf_index)
	{
// TODO: see which are needed		
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
//		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
//		final double [][] surfaces = getSurfaces();
		final int [][] shell_map =   getShellMap();
		final TilePlanes.PlaneData [][] planes = getPlanesMod(); // to get strength and connections
		final int [] neibs = new int [selected.length];
		final int [] dirs8 =   {-tilesX,   -tilesX +  1, 1, tilesX + 1,  tilesX,  tilesX - 1,  -1, -tilesX -  1};
//		final int [] dirs8st = {-stilesX,  -stilesX + 1, 1, stilesX + 1, stilesX, stilesX - 1, -1, -stilesX - 1};
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final int surf_index_plus1 = surf_index + 1;
		final int emptyTile = -256; 
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < snap_surf.length; nTile = ai.getAndIncrement()) {
						neibs[nTile] = emptyTile;
						if ((selected == null) || selected[nTile]) {
							int np = snap_surf[nTile] - 1;
							if (np >= 0) {
								int tX = nTile % tilesX; 
								int tY = nTile / tilesX;
								int stX = tX / superTileSize;
								int stY = tY / superTileSize;
								int fX = tX % superTileSize;
								int fY = tY % superTileSize;
								int nsTile = stY * stilesX + stX;
								if (nsTile == 724){
//									System.out.println("getNeighbors(): tX:tY = "+tX+":"+tY+", stX:stY="+stX+":"+stY+
//											", nTile="+nTile+", nsTile="+nsTile);
								}
								if (shell_map[nsTile][np] == surf_index_plus1) {
									int tileMask = getNeibMask(
											tX,
											tY,
											tilesX,
											tilesY);
									int localMask = getNeibMask(
											fX,
											fY,
											superTileSize,
											superTileSize);
									int [] dir_dir = getNeibBordDir(
											fX,
											fY,
											superTileSize,
											superTileSize);
									neibs[nTile] = tileMask;
									for (int i = 0; i < 8; i++) {
										int b = 1 << i;
										if ((neibs[nTile] & b) != 0){ // that will prevent tile array bounds violation
											int nTile1 = nTile + dirs8[i];
											if (!(((selected == null) || selected[nTile1]) && (snap_surf[nTile1] == snap_surf[nTile]))){
												neibs[nTile] &= ~b; // same supertile, different index
											} else if ((localMask & b) == 0) {// going to the other supertile
												int [] neib_best = planes[nsTile][np].getNeibBest(); // array of connected planes indices/ -1s
												if ((dir_dir[i] < 0) || (neib_best[dir_dir[i]] != (snap_surf[nTile1] -1))) {
													neibs[nTile] &= ~b;
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
		return neibs;

	}
	
	
	
	
	/**
	 * Sort per-tile disparity map values to one of the surfaces	
	 * @param disparity per-tile array of disparity values
	 * @param strength  per-tile array of correlation strength values
	 * @param selection tile selection (or null to process all tiles)
	 * @param dispNorm  disparity difference normalization - differences for average disparities
	 *                  below this value are unchanged, above - proportionally scaled down)
	 * @param snapDispAny tiles with disparity difference below this will snap to the nearest
	 *                  surface for any correlation strength 
	 * @param snapStrengthAny tiles that have lower strength and can not be snapped otherwise are
	 *                  treated as zero-strength and obey snapZeroMode parameter
	 * @param snapNegAny allow snap for tiles that are below (farther) than the closest surface,
	 *                  if the disparity difference is below this value
	 * @param snapDispMax do not snap low-strength tiles if they are farther from the nearest surface
	 *                  than this.
	 * @param snapDispWeight snap tiles in the range snapDispAny and snapDispMax if the product
	 * 					of the disparity difference by the strength is below this
	 * @param snapZeroMode : 0 - no special treatment, 1 snap to strongest surface,
	 *                  2 - snap to the farthest (lowest disparity) surface
	 * @param debugLevel debug level
	 * @return integer array with per-tile elements. 0 - no surface candidate found (keep original
	 *                  disparity), > 0 - number of the plane plus 1 (use shell_map to get
	 *                  the surface number) 
	 */

	public int [] snapSort(
			final double []  disparity,
			final double []  strength,
			final boolean [] selection,       // can be null
			final double     dispNorm, // plDispNorm
			final double     snapDispAny,    //    =  .2,   // Maximal (scaled by plDispNorm) disparity difference to snap to plane at any strength
			final double     snapStrengthAny,    //    =  .2,   //  Maximal strength to fit any distance (if does not fit otherwise - treat as zero strength
			final double     snapNegAny,     //    =  .2,   // Maximal negative disparity difference from the best match
			final double     snapDispMax,    //    =  .5,   // Maximal (scaled by plDispNorm) disparity difference to snap to plane at low strength
			final double     snapDispWeight, //    =  .5,   // Maximal disparity diff. by weight product to snap to plane
			final int        snapZeroMode,   // Zero strength snap mode: 0: no special treatment, 1 - strongest, 2 - farthest
			final int        debugLevel)
	{
		final int tilesX =        tileProcessor.getTilesX();
//		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
//		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final double [][] surfaces = getSurfaces();
		final int [][] shell_map =   getShellMap();
		final TilePlanes.PlaneData [][] planes = getPlanesMod(); // to get strength;
//		final double [] snap_disp =  disparity.clone();
		final int [] snap_sort =  new int [disparity.length];
		
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < disparity.length; nTile = ai.getAndIncrement()) {
						if ((selection == null) || selection[nTile]) {
							int tX = nTile % tilesX; 
							int tY = nTile / tilesX;
							int stX = tX / superTileSize;
							int stY = tY / superTileSize;
							int nsTile = stY * stilesX + stX;
							int nTile_st = (tY * stilesX * superTileSize) + tX; // surfaces maybe (is) larger than tilesX *tilesY if 
							// tilesX, tilesY are not multiples of superTileSize !
							if (nsTile == 724){
//								System.out.println("snapDisparity(): tX:tY = "+tX+":"+tY+", stX:stY="+stX+":"+stY+
//										", nTile="+nTile+", nsTile="+nsTile+", nTile_st = "+nTile_st);
							}
							
							// is there any surface to snap to?
							boolean no_surf = true;
							if (shell_map[nsTile] != null){
								for (int np = 0; no_surf && (np < shell_map[nsTile].length); np++){
									if (shell_map[nsTile][np] > 0) no_surf = false;
								}
							}
							if (!no_surf) {
								double [] surf_disparity = new double[shell_map[nsTile].length];
								for (int np = 0; np < surf_disparity.length; np++){
									if (shell_map[nsTile][np] > 0) {
										surf_disparity[np] = surfaces[shell_map[nsTile][np] -1][nTile_st];
									} else {
										surf_disparity[np] = Double.NaN;
									}
								}
								int best_indx = -1;
								double best_disp_diff = Double.NaN;
								double this_disp = disparity[nTile];
								double this_strength = strength[nTile];
								if (Double.isNaN(this_disp)){
									this_disp = 0.0;
									this_strength = 0.0;
								}
								for (int np = 0; np < surf_disparity.length; np++){
									if (shell_map[nsTile][np] > 0) {
										double disp_diff = normDispADiff (this_disp,surf_disparity[np], dispNorm);
										if ((best_indx < 0) || (disp_diff < best_disp_diff)) {
											best_disp_diff = disp_diff;
											best_indx =      np;
										}
									}											
								}										
								if ((this_strength > 0.0) || (snapZeroMode == 0)){
									if ((best_disp_diff <= snapDispAny) ||
											((best_disp_diff <= snapDispMax) && (best_disp_diff * this_strength < snapDispWeight )) ||  
											((best_disp_diff <= snapNegAny) && (this_disp <  surf_disparity[best_indx]))  // farther than surface - higher tolerance
											){
										snap_sort[nTile] = best_indx+1; // snapped (default - keep source disparity
									} else if ((this_strength < snapStrengthAny) && (snapZeroMode != 0)){ // treat as zero strength?
										this_strength = 0.0;
									}
								}
								
								if ((this_strength == 0) && (snapZeroMode > 0)) { 
									if (snapZeroMode == 1) {  // zero strength and non-zero snapZeroMode, consider 1 (strongest shell)
										best_indx = -1;
										double best_strength = Double.NaN;
										for (int np = 0; np < surf_disparity.length; np++){
											if (shell_map[nsTile][np] > 0) {
												double w = planes[nsTile][np].getWeight();
												if ((best_indx < 0) || (w >  best_strength)) {
													best_strength = w;
													best_indx =     np;
												}
											}
										}
										snap_sort[nTile] = best_indx+1; // snapped (default - keep source disparity
									} else { // zero strength and snapZeroMode > 1 (now only 2) - farthest surface (smallest disparity)
										best_indx = -1;
										double smallest_disparity = Double.NaN;
										for (int np = 0; np < surf_disparity.length; np++){
											if (shell_map[nsTile][np] > 0) {
												if ((best_indx < 0) || ( surf_disparity[np] <  smallest_disparity)) {
													smallest_disparity = surf_disparity[np];
													best_indx =     np;
												}
											}
										}
										snap_sort[nTile] = best_indx+1;
									}
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return snap_sort;
	}


} // end of class SuperTiles
