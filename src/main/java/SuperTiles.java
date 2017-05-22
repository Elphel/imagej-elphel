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
import java.awt.Point;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicInteger;


//import java.awt.Point;
//import java.util.ArrayList;

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
	int         measSel = 1; // bitmask of the selected measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
	
	boolean     smplMode        = true;   // Use sample mode (false - regular tile mode)
	int         smplSide        = 2;      // Sample size (side of a square)
	int         smplNum         = 3;      // Number after removing worst
	double      smplRms         = 0.1;    // Maximal RMS of the remaining tiles in a sample
	
	
	
	MeasuredLayers measuredLayers = null;
	
	CLTPass3d cltPass3d;
	TilePlanes.PlaneData [][] planes =  null;
	TilePlanes.PlaneData [][] planes_mod = null;
	double [][][] fuse_coeff = null; // coefficients to fuse planes in neighbor supertiles

	int [][]    shell_map = null;   // per supertile, per disparity plane - shell index + 1 (0 - none)
	double [][] surfaces;  // per shell, per tile (linescan order) disparity value or NaN in missing supertiles
	
	TileSurface tileSurface = null;	
	

	/**
	 * currently lowest plane for each includes all tiles, so do not use it. May change in the future
	 * @param np total number of planes in a supertile
	 * @return lowest plane to iterate through
	 */
	
	static int LOWEST_PLANE(int np){ // 
		return (np >1) ? 1 : 0;
	}

	public SuperTiles(
			CLTPass3d cltPass3d,
			double                  step_near,
			double                  step_far,
			double                  step_threshold,
			double                  min_disparity,
			double                  max_disparity,
			double                  strength_floor,
			double                  strength_pow,
			double                  stBlurSigma,
			boolean                 smplMode, //        = true;   // Use sample mode (false - regular tile mode)
			int                     smplSide, //        = 2;      // Sample size (side of a square)
			int                     smplNum, //         = 3;      // Number after removing worst
			double                  smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
			int                     measSel)
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
		this.smplMode        = smplMode;   // Use sample mode (false - regular tile mode)
		this.smplSide        = smplSide;   // Sample size (side of a square)
		this.smplNum         = smplNum;    // Number after removing worst
		this.smplRms         = smplRms;    // Maximal RMS of the remaining tiles in a sample
		this.measSel =        measSel;
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
//		initFuseCoeff(0.5, true); // true);
		initFuseCoeff(0.5, false); // true);
		
// Set up MeasuredLayers		
		measuredLayers = new MeasuredLayers(
				4, // combo, quad, hor, vert
				tileProcessor.getTilesX(),
				tileProcessor.getTilesY(),
				tileProcessor.getSuperTileSize());
// MeasuredLayer.setLayer(int, double[], double[], boolean[])
		measuredLayers.setLayer( // combo disparity
				0,                                // int       num_layer,
				cltPass3d.getDisparity(0),        // double [] disparity, 
				cltPass3d.getStrength(),          // double [] strength,
				cltPass3d.getSelected());         // boolean [] selection) // may be null
		measuredLayers.setLayer(                  // quad disparity
				1,                                // int       num_layer,
				cltPass3d.getDisparity(1),        // double [] disparity, 
				cltPass3d.getOriginalStrength(),  // double [] strength,
				cltPass3d.getSelected());         // boolean [] selection) // may be null
		measuredLayers.setLayer(                  // hor disparity
				2,                                // int       num_layer,
				cltPass3d.getDisparity(2),        // double [] disparity, 
				cltPass3d.getHorStrength(),       // double [] strength,
				cltPass3d.getSelected());         // boolean [] selection) // may be null
		measuredLayers.setLayer(                  // vert disparity
				3,                                // int       num_layer,
				cltPass3d.getDisparity(3),        // double [] disparity, 
				cltPass3d.getVertStrength(),      // double [] strength,
				cltPass3d.getSelected());         // boolean [] selection) // may be null

		getDisparityHistograms(
				null,       // double  []   world_plane, // tilt equi-disparity planes to match real world planes (usually horizontal (or null)
				null,       // boolean [][] tile_sel, // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all
				smplMode,   // final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				smplSide,   // final int        smplSide, //        = 2;      // Sample size (side of a square)
				smplNum,    // final int        smplNum,  //         = 3;      // Number after removing worst
				smplRms,    // final double     smplRms,  //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
				measSel);   // calculate and blur supertiles (for all, not just selected?)
		if (tileProcessor.globalDebugLevel > 0){
			System.out.println("SuperTiles(): min_disparity = "+min_disparity+", max_disparity="+max_disparity);
			System.out.println("SuperTiles(): step_far = "+step_far+", step_near="+step_near);
			System.out.println("SuperTiles(): numBins = "+numBins+", bin_far="+bin_far+", bin_near = "+bin_near);
			System.out.println("SuperTiles(): step_threshold_far = "+step_threshold_far+", step_threshold_near="+step_threshold_near);
			for (int i = 0; i<numBins; i++){
				System.out.println(i+": "+bin_centers[i]+", "+disparityToBin(bin_centers[i]));
			}
			//
		}
		if (tileProcessor.globalDebugLevel > 0){
			String [] titles = {"d0","s0","d1","s1","d2","s2","d3","s3"};
			double [][] dbg_img = new double [titles.length][];
			for (int i = 0; i < measuredLayers.getNumLayers(); i++){
				dbg_img[2 * i] =     measuredLayers.getDisparity(i);
				dbg_img[2 * i + 1] = measuredLayers.getStrength(i);
			}
			showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
			sdfa_instance.showArrays(dbg_img,  tileProcessor.getTilesX(), tileProcessor.getTilesY(), true, "measuredLayers",titles);
		}
	}
	public TileProcessor getTileProcessor(){
		return  tileProcessor;
	}
	
	public void setTileSurface( TileSurface tileSurface)
	{
		this.tileSurface = tileSurface;
	}
	public void setTileSurface(
			GeometryCorrection   geometryCorrection)
	{
		this.tileSurface = new TileSurface(
				tileProcessor.getTileSize(),      // int tileSize,
				tileProcessor.getSuperTileSize(), // int superTileSize,
				tileProcessor.getTilesX(),        // int tilesX,
				tileProcessor.getTilesY(),        // int tilesY,
				geometryCorrection,               // GeometryCorrection geometryCorrection,
				tileProcessor.threadsMax);        // int threadsMax); 
	}

	
	public TileSurface getTileSurface()
	{
		return tileSurface;
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

	public void resetDisparityHistograms()
	{
		this.disparityHistograms = null;
	}

	public void setBlurSigma(double blurSigma) // resets histograms, as it needs re-calculating
	{
		this.stBlurSigma =         blurSigma;
		this.disparityHistograms = null;
	}
	
	public double [][] getDisparityHistograms(
			final double [][][][] disparity_strength, // pre-calculated disparity/strength [per super-tile][per-measurement layer][2][tiles] or null
			final boolean [][]    tile_sel, // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all
			final boolean         smplMode, //        = true;   // Use sample mode (false - regular tile mode)
			final int             smplSide, //        = 2;      // Sample size (side of a square)
			final int             smplNum,  //         = 3;      // Number after removing worst
			final double          smplRms,  //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
			final int             measSel)  //
	{
		if (disparity_strength != null) {
			System.out.println("getDisparityHistograms() with non-null disparity_strength"); 
		}
		if ((this.disparityHistograms != null) &&
				(smplMode == this.smplMode) &&
				(smplSide == this.smplSide) &&
				(smplNum  == this.smplNum) &&
				(smplRms  == this.smplRms) &&
				(measSel  == this.measSel)){
				return this.disparityHistograms;
		}
		this.smplMode = smplMode;   // Use sample mode (false - regular tile mode)
		this.smplSide = smplSide;   // Sample size (side of a square)
		this.smplNum  = smplNum;    // Number after removing worst
		this.smplRms  = smplRms;    // Maximal RMS of the remaining tiles in a sample
		this.measSel  = measSel;
		
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int nStiles = stilesX * stilesY; 
		final double [][] dispHist =     new double [nStiles][]; // now it will be sparse 
		final double []   strengthHist = new double [nStiles];
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
//		final TilePlanes tpl = new TilePlanes(tileProcessor.getTileSize(),superTileSize);

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					
					for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
						int stileY = nsTile / stilesX;  
						int stileX = nsTile % stilesX;
						double sw = 0.0; // sum weights
						double [] hist = new double [numBins];
						//						double [][][] disp_strength = new double [measuredLayers.getNumLayers()][][] ;
						for (int nl = 0; nl < measuredLayers.getNumLayers(); nl ++) {
							if (((measSel & (1 << nl)) != 0) && ((tile_sel == null) || ((tile_sel[nl] != null)))) {
								double [][] disp_strength;
								if (disparity_strength == null) {
									if (smplMode) {
										disp_strength =  measuredLayers.getDisparityStrength(
												nl,             // int num_layer,
												stileX,         // int stX,
												stileY,         // int stY,
												(((tile_sel == null) || (tile_sel[nl].length == 0))? null:tile_sel[nl]), // boolean [] sel_in,
												//											null,           // boolean [] sel_in,
												strength_floor, // double strength_floor,
												strength_pow,   // double strength_pow,
												smplSide,       // int        smplSide, // = 2;   // Sample size (side of a square)
												smplNum,        //int        smplNum,   // = 3;   // Number after removing worst (should be >1)
												smplRms,        //double     smplRms,   // = 0.1; // Maximal RMS of the remaining tiles in a sample
												true);          // boolean null_if_none);
									} else {
										disp_strength =  measuredLayers.getDisparityStrength(
												nl,             // int num_layer,
												stileX,         // int stX,
												stileY,         // int stY,
												(((tile_sel == null) || (tile_sel[nl].length == 0))? null:tile_sel[nl]), // boolean [] sel_in,
												//											null,           // boolean [] sel_in,
												strength_floor, // double strength_floor,
												strength_pow,   // double strength_pow,
												true);          // boolean null_if_none);
									}
								} else {
									disp_strength =  (disparity_strength[nsTile] != null) ? disparity_strength[nsTile][nl] : null;
								}
								if (disp_strength != null) {
									for (int indx = 0; indx < disp_strength[1].length; indx++) {
										double w = disp_strength[1][indx]; 
										if ( w > 0.0){
											double d = disp_strength[0][indx]; 
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
			blurDisparityHistogram(0); // debugLevel);
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

	public double [][][] getMaxMinMax(
			final double [][][][] disparity_strength, // pre-calculated disparity/strength [per super-tile][per-measurement layer][2][tiles] or null
			final boolean [][] tile_sel // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all
			){
		// first find all integer maximums, and if the top is flat - use the middle. If not flat - use 2-nd degree polynomial
		if (disparityHistograms == null) getDisparityHistograms(
//				world_plane, // double  []   world_plane, // tilt equi-disparity planes to match real world planes (usually horizontal (or null)
				disparity_strength, // pre-calculated disparity/strength [per super-tile][per-measurement layer][2][tiles] or null				
				tile_sel,   // boolean [][] tile_sel, // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all
				this.smplMode,   // final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				this.smplSide,   // final int        smplSide, //        = 2;      // Sample size (side of a square)
				this.smplNum,    // final int        smplNum,  //         = 3;      // Number after removing worst
				this.smplRms,    // final double     smplRms,  //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
				this.measSel);
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
			getDisparityHistograms(
					null, // double  []   world_plane, // tilt equi-disparity planes to match real world planes (usually horizontal (or null)
					null,   // boolean [][] tile_sel, // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all
					this.smplMode,   // final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
					this.smplSide,   // final int        smplSide, //        = 2;      // Sample size (side of a square)
					this.smplNum,    // final int        smplNum,  //         = 3;      // Number after removing worst
					this.smplRms,    // final double     smplRms,  //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
					this.measSel); // calculate and blur with the current settings, specified at instantiation
		}
		return showDisparityHistogram(disparityHistograms);
	}

	
	public double [] showDisparityHistogram(
//			double  []   world_plane, // tilt equi-disparity planes to match real world planes (usually horizontal (or null)
			final double [][][][] disparity_strength, // pre-calculated disparity/strength [per super-tile][per-measurement layer][2][tiles] or null
			boolean [][] tile_sel, // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all
			boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
			int        smplSide, //        = 2;      // Sample size (side of a square)
			int        smplNum,  //         = 3;      // Number after removing worst
			double     smplRms,  //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
			int        measSel) // bitmask of the selected measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
	{
		getDisparityHistograms( // will recalculate if does not exist or some parameters changed
//				world_plane, // double  []   world_plane, // tilt equi-disparity planes to match real world planes (usually horizo				
				disparity_strength, // pre-calculated disparity/strength [per super-tile][per-measurement layer][2][tiles] or nullntal (or null)
				tile_sel,   // boolean [][] tile_sel, // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all
				smplMode,   // final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				smplSide,   // final int        smplSide, //        = 2;      // Sample size (side of a square)
				smplNum,    // final int        smplNum,  //         = 3;      // Number after removing worst
				smplRms,    // final double     smplRms,  //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
				measSel); // calculate and blur with the current settings, specified at instantiation
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
			return showMaxMinMax( // calculate and blur with the current settings, specified at instantiation
					null,
					null);
	}
	public double [] showMaxMinMax(
//			double  []   world_plane, // tilt equi-disparity planes to match real world planes (usually horizontal (or null)
			double [][][][] disparity_strength, // pre-calculated disparity/strength [per super-tile][per-measurement layer][2][tiles] or null
			
			boolean [][] tile_sel // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all
		){
		if (maxMinMax == null){ 
			getMaxMinMax( // calculate and blur with the current settings, specified at instantiation
//					world_plane, // tilt equi-disparity planes to match real world planes (usually horizontal (or null)
					disparity_strength, // pre-calculated disparity/strength [per super-tile][per-measurement layer][2][tiles] or null
					tile_sel); // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all
					
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

			// draw rectangular frame - horizontal dotted lines
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


	public void processPlanes3(
			final boolean [] selected, // or null
			final double     min_disp,
//			final boolean    invert_disp, // use 1/disparity
			final int        stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			final double     plDispNorm,
			final int        plMinPoints, //          =     5;  // Minimal number of points for plane detection
			final double     plTargetEigen, //        =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
			final double     plFractOutliers, //      =   0.3;  // Maximal fraction of outliers to remove
			final int        plMaxOutliers, //        =    20;  // Maximal number of outliers to remove
			final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
			final GeometryCorrection geometryCorrection,
			final boolean    correct_distortions,
			final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
			final int        smplSide, //        = 2;      // Sample size (side of a square)
			final int        smplNum, //         = 3;      // Number after removing worst
			final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
			
			final int        debugLevel,
			final int        dbg_X,
			final int        dbg_Y)
	{
		if (maxMinMax == null) getMaxMinMax(null, null);
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int tileSize =      tileProcessor.getTileSize();

		final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int nStiles = stilesX * stilesY; 
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		this.planes = new TilePlanes.PlaneData[nStiles][];
		final int debug_stile = (debugLevel > -1)? (dbg_Y * stilesX + dbg_X):-1;
		
//		final boolean [][] dflt_select = {{}, null, null, null, null}; // use layer 0 (combo) only
		final boolean [][] dflt_select = new boolean [measuredLayers.getNumLayers()][];
		for (int i = 0; i < dflt_select.length; i++){
			if ((stMeasSel & (1 << i)) !=0){
				dflt_select[i] = new boolean[0];
			} else {
				dflt_select[i] = null;
			}
		}
		
		// TODO: Remove when promoting PlaneData
		final TilePlanes tpl = new TilePlanes(tileSize,superTileSize, geometryCorrection);
		
//		final double  [] disparity = cltPass3d.getDisparity();
//		final double  [] strength =  cltPass3d.getStrength();
		
		measuredLayers.setLayer (
				0, // int       num_layer,
				cltPass3d.getDisparity(), // double [] disparity,
				cltPass3d.getStrength(), // double [] strength,
				null); // boolean [] selection) // may be null
		if (debugLevel > -1) {
			String [] titles = {"d0","s0","d1","s1","d2","s2","d3","s3","s","d"};
			double [][] dbg_img = new double [titles.length][];
			for (int i = 0; i < measuredLayers.getNumLayers(); i++){
				dbg_img[2 * i] =     measuredLayers.getDisparity(i);
				dbg_img[2 * i + 1] = measuredLayers.getStrength(i);
			}
			dbg_img[8] = cltPass3d.getDisparity();
			dbg_img[9] = cltPass3d.getStrength();
			showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
			sdfa_instance.showArrays(dbg_img,  tileProcessor.getTilesX(), tileProcessor.getTilesY(), true, "measuredLayers",titles);
			
			
		}
//		if (maxMinMax == null) 
			getMaxMinMax(null, null);

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
						if (nsTile == debug_stile){
							System.out.println("processPlanes3(): nsTile="+nsTile);
						}
						int stileY = nsTile / stilesX;  
						int stileX = nsTile % stilesX;
						int [] sTiles = {stileX, stileY};
						planes[nsTile] = null;
						// first make a plane from all tiles
						TilePlanes.PlaneData pd0 = tpl.new  PlaneData (
								sTiles, // int [] sTileXY, 
								tileSize, // int tileSize,
								geometryCorrection, // GeometryCorrection   geometryCorrection,
								correct_distortions,
								measuredLayers,     // MeasuredLayers measuredLayers,
								plPreferDisparity);   // boolean preferDisparity)
						boolean [][] tile_sel = dflt_select.clone();
						for (int i = 0; i < dflt_select.length; i++){
							if (dflt_select[i] != null) tile_sel[i] = dflt_select[i].clone();
						}
//						int dl1 =  (nsTile == debug_stile) ? 3 : 0;
						int dl =  (nsTile == debug_stile) ? 4 : 0;

						double [][][] disp_strength = pd0.getPlaneFromMeas(
								tile_sel,     // boolean [][] tile_sel, // null - do not use, {} use all (will be modified)
								null,
								min_disp,     // double       disp_far, // minimal disparity to select (or NaN)
								Double.NaN,   // double       disp_near, // maximal disparity to select (or NaN)
								0.0,          // plDispNorm,  // double       dispNorm,   //  Normalize disparities to the average if above
								0.0,          // double       min_weight,
								plMinPoints,  // int          min_tiles,
								strength_floor, // 
								strength_pow, // double       strength_pow,
								smplMode,
								smplSide,
								smplNum,
								smplRms,
								dl);          // int          debugLevel)


						if (disp_strength != null){ // there are some non-zero tiles, process them (all points, not clustered by disparity value)
							boolean OK;
							TilePlanes.PlaneData pd0_full = pd0.clone();  //
							ArrayList<TilePlanes.PlaneData> st_planes = new ArrayList<TilePlanes.PlaneData>();
							// now try to remove outliers
							int max_outliers = (int) Math.round(pd0.getNumPoints() * plFractOutliers);
							if (max_outliers > plMaxOutliers) max_outliers = plMaxOutliers;
							double targetV = corrMaxEigen(
									plTargetEigen,
									plDispNorm,
									pd0);
							if (pd0.getValue() > targetV) {
								OK = pd0.removeOutliers( // getPlaneFromMeas should already have run
										disp_strength,
										targetV,      // double     targetEigen, // target eigenvalue for primary axis (is disparity-dependent, so is non-constant)
										max_outliers, // int        maxRemoved,  // maximal number of tiles to remove (not a constant)
										dl); // int        debugLevel)
								if (!OK) continue;
								if (dl > 0) {
									if (pd0.getWeight() > 0.3) { // 1.0) {
										System.out.println("Removed outliers["+nsTile+"]"+
												", stileX="+stileX+
												", stileY="+stileY+
												", numPoints="+ pd0.getNumPoints()+
												", swc = "+pd0.getWeight()+
												", center=["+pd0.getZxy()[0]+","+pd0.getZxy()[1]+","+pd0.getZxy()[2]+"]"+
												", eig_val = {"+pd0.getValues()[0]+","+pd0.getValues()[1]+","+pd0.getValues()[2]+"}"+
												", eig_vect[0] = {"+pd0.getVector()[0]+","+pd0.getVector()[1]+","+pd0.getVector()[2]+"}");
									}
								}
							} // nothing to do if already OK
							if (dl > 0) {
								System.out.println("Calculating World normal["+nsTile+"]");
							}
							double swc_common = pd0.getWeight();
							double [] norm_xyz = pd0.getWorldXYZ(
									correct_distortions,
									dl);
							if (dl > 0) {
								System.out.println("World normal["+nsTile+"] = {"+
										norm_xyz[0]+", "+norm_xyz[1]+", "+norm_xyz[2]+"}");

							}
							st_planes.add(pd0);  // adding [0] - all supertile tiles, not clustered by disparity value

							// now try for each of the disparity-separated clusters (only for multi-peak histograms)


							double [][] mm = maxMinMax[nsTile];
							if (mm == null){
//								double [][][] dbg_min_max = maxMinMax;
//								System.out.println("maxMinMax["+nsTile+"] == null");
							}

							if ((mm!= null) &&  (mm.length > 1)) { // multiple maximums - separate into multiple selections // null pointer
								for (int m = 0; m < (mm.length +1)/2; m++){
									double [] far_near = new double [2];
									if (m == 0) {
										far_near[0] = Double.NaN; // pd0 already filtered by min_disp
									} else {
										far_near[0] = mm[2 * m - 1][0];
									}
									if (m == (mm.length -1)/2) {
										far_near[1] = Double.NaN; // pd0 already filtered by min_disp
									} else {
										far_near[1] = mm[2 * m + 1][0];
									}


									TilePlanes.PlaneData pd = pd0_full.clone(); 
									OK = (pd.getPlaneFromMeas(
											null, // tile_sel,       // boolean [][] tile_sel, // null - do not use, {} use all (will be modified)
											disp_strength, 
											far_near[0],    // double       disp_far, // minimal disparity to select (or NaN)
											far_near[1],    // double       disp_near, // maximal disparity to select (or NaN)
											0.0,            // plDispNorm,  // double       dispNorm,   //  Normalize disparities to the average if above
											0.0,            // double       min_weight,
											plMinPoints,    // int          min_tiles,
											strength_floor, // 
											strength_pow,   // double       strength_pow,
											smplMode,
											smplSide,
											smplNum,
											smplRms,
											dl) != null);            // int          debugLevel)
									if (OK) {
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
											OK = pd.removeOutliers( // getPlaneFromMeas should already have run
													disp_strength, 
													targetV,      // double     targetEigen, // target eigenvalue for primary axis (is disparity-dependent, so is non-constant)
													max_outliers, // int        maxRemoved,  // maximal number of tiles to remove (not a constant)
													dl); // int        debugLevel)
											if (!OK) {
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
							if (st_planes.size() > 0){
								planes[nsTile] = st_planes.toArray(new TilePlanes.PlaneData[0] );
								if (dl >0){
									System.out.println("processPlanes3(): nsTile="+nsTile);
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}

	public int [][] getTransMatrix(
			boolean [][][] selections ) // for each plane should have the same non-null ml
	{
		int num_ml = measuredLayers.getNumLayers();
		int num_p = selections.length;
		int superTileSize = tileProcessor.getSuperTileSize();
		int st2 = 2 * superTileSize; 
		int [][] trans_mat = new int [num_p][num_p];
		for (int ml = 0; ml < num_ml; ml++) if (selections[0][ml] != null){
			for (int y = 0; y < st2; y++){
				for (int x = 0; x < (st2 - 1); x++){
					int indx1 = y * st2 + x;
					int indx2 = y * st2 + x + 1;
					for (int np1 = 0; np1 < num_p; np1++){
						for (int np2 = 0; np2 < num_p; np2++){
							if (selections[np1][ml][indx1] && selections[np2][ml][indx2]){
								trans_mat[np1][np2]++;
							}
						}
					}
				}
			}
			for (int x = 0; x < st2; x++){
				for (int y = 0; y < (st2 - 1); y++){
					int indx1 = y * st2 + x;
					int indx2 = y * st2 + x + st2;
					for (int np1 = 0; np1 < num_p; np1++){
						for (int np2 = 0; np2 < num_p; np2++){
							if (selections[np1][ml][indx1] && selections[np2][ml][indx2]){
								trans_mat[np1][np2]++;
							}
						}
					}
				}
			}
		}
		return trans_mat;
	}
	
	public double [][] getTransRel(
			int [][] trans_matrix)
	{
		double [][] trans_rel = new double [trans_matrix.length][trans_matrix.length];
		for (int i = 0; i < trans_matrix.length; i++){
			for (int j = 0; j < trans_matrix[i].length; j++){
				if ((trans_matrix[i][i] + trans_matrix[j][j]) > 0){
					trans_rel[i][j] = (trans_matrix[i][j] + trans_matrix[j][i]);
					trans_rel[i][j] /= (trans_matrix[i][i] + trans_matrix[j][j]);
				}
			}
		}
		return trans_rel;
	}
	
	public String [] showSupertileSeparationTitles(
			double [][][]  disp_strength,
			boolean [][][] selections)
	{
		return 	showSupertileSeparationTitles(disp_strength, selections, null);
	}
	
	public double [][] showSupertileSeparation(
			boolean        useWorld,
			double [][][]  disp_strength,
			boolean [][][] selections)
	{
		return 	showSupertileSeparation(useWorld, disp_strength, selections, null);
	}

	public String [] showSupertileSeparationTitles(
			double [][][]  disp_strength,
			boolean [][][] selections,
			TilePlanes.PlaneData [] planes)
	{
		int num_ml = disp_strength.length;
		int num_p  = (selections == null) ? 0: selections.length;
		int num_pm = num_ml * num_p;
		int num_pd = (planes != null) ? (planes.length - LOWEST_PLANE(planes.length)) : 0; 
		String [] titles = new String [num_pm + 3 * num_ml + 4 * num_pd];
		for (int np = 0; np < num_p; np++){
			for (int ml = 0; ml < num_ml; ml++){
				titles [np * num_ml + ml] = "p"+np+"_l"+ml;
			}
		}
		for (int ml = 0; ml < num_ml; ml++){
			titles [num_pm + 0 * num_ml + ml] = "disp_"+ml;
			titles [num_pm + 1 * num_ml + ml] = "strn_"+ml;
			titles [num_pm + 2 * num_ml + ml] = "sel_"+ml;
		}
		for (int npd = 0; npd < num_pd; npd++){
			titles [num_pm + 3 * num_ml + 0 * num_pd + npd] = "pd_"+npd;
			titles [num_pm + 3 * num_ml + 1 * num_pd + npd] = "pdm_"+npd;
			titles [num_pm + 3 * num_ml + 2 * num_pd + npd] = "pdd_"+npd;
			titles [num_pm + 3 * num_ml + 3 * num_pd + npd] = "pds_"+npd;
		}

	// TilePlanes.PlaneData	
		return titles;
	}

	public double [][] showSupertileSeparation(
			boolean useWorld,
			double [][][]  disp_strength,
			boolean [][][] selections,
			TilePlanes.PlaneData [] planes
			)
	{
		int superTileSize = tileProcessor.getSuperTileSize();
		int num_ml = disp_strength.length;
		int num_p  = (selections == null) ? 0: selections.length;
		int num_pm = num_ml * num_p;
		int num_pd = (planes != null) ? (planes.length - LOWEST_PLANE(planes.length)) : 0; 
		double [][] data = new double [num_pm + 3 * num_ml + 4 * num_pd ][]; // 4* superTileSize*superTileSize];
		for (int np = 0; np < num_p; np++) if (selections [np] != null){
			for (int ml = 0; ml < num_ml; ml++) if ((disp_strength[ml]!=null) && (selections[np][ml] != null)){
				int nd = np * num_ml + ml; 
				data[nd] = new double[4 * superTileSize * superTileSize];
				for (int i = 0; i < data[nd].length; i++){
					if (selections[np][ml][i] && (disp_strength[ml][1][i] > 0.0)){
						data[nd][i] = disp_strength[ml][0][i];		
					} else {
						data[nd][i] = Double.NaN;		
					}
				}
			}
		}
		for (int ml = 0; ml < num_ml; ml++) if (disp_strength[ml]!=null){
			int nd = num_pm + 2 * num_ml + ml;	
		
			data [nd] = new double [4* superTileSize*superTileSize];
			for (int i = 0; i < data[nd].length; i++){
				data [nd][i] = Double.NaN; 		
				for (int np = 0; np < num_p; np++) if ((selections [np] != null) && (selections [np][ml] != null) && selections [np][ml][i]){
					data [nd][i] = np + 1; 		
					break;
				}
			}
			data [num_pm + 0 * num_ml + ml] = disp_strength[ml][0];
			data [num_pm + 1 * num_ml + ml] = disp_strength[ml][1];
		}
		for (int npd = 0; npd < num_pd; npd++) if (planes[npd +LOWEST_PLANE(planes.length)] != null){
			double [][] ellipsoids = planes[npd +LOWEST_PLANE(planes.length)].getDoublePlaneDisparityStrength(
					useWorld,
					null, // double [] window,
					true, // boolean   use_sel,
					true, // boolean   divide_by_area,
					1.5, // double   scale_projection,
					1); // int       debugLevel)
			
			data [num_pm + 3 * num_ml + 0 * num_pd + npd] = planes[npd +LOWEST_PLANE(planes.length)].getDoublePlaneDisparity(
					useWorld,
					false);
			
			boolean [][] msel = planes[npd +LOWEST_PLANE(planes.length)].getMeasSelection();
			int ntiles = data [num_pm + 3 * num_ml + 0 * num_pd + npd].length; 
			data [num_pm + 3 * num_ml + 1 * num_pd + npd] = new double [ntiles];
			for (int i = 0; i < ntiles; i++) {
				data [num_pm + 3 * num_ml + 1 * num_pd + npd][i] = Double.NaN;
				for (int ml = 0; ml < num_ml; ml++) if ((disp_strength[ml]!=null) && (msel[ml] != null)){
					if (msel[ml][i] && (disp_strength[ml][1][i] > 0.0)){
						data [num_pm + 3 * num_ml + 1 * num_pd + npd][i] = data [num_pm + 3 * num_ml + 0 * num_pd + npd][i];
					}					
				}				
			}	
			data [num_pm + 3 * num_ml + 2 * num_pd + npd] = ellipsoids[0];
			data [num_pm + 3 * num_ml + 3 * num_pd + npd] = ellipsoids[1];
		}
//		public boolean [] getMeasSelection(int nl){
		
		return data;
	}
	
	// calculate "tilted" disparity, so planes parallel to the same world plane would have the same disparity
	// also produces non-tilted, if world_plane_norm == null
	// Protecting from behind the horizon - set strength of all tiles in the negative disparity area to 0
	public double [][][][] getPlaneDispStrengths(
			final double []  world_plane_norm, // real world normal vector to a suggested plane family (0,1,0) for horizontal planes
			final int        stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert

			final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
			final GeometryCorrection geometryCorrection,
			final boolean    correct_distortions,

			final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
			final int        smplSide, //        = 2;      // Sample size (side of a square)
			final int        smplNum, //         = 3;      // Number after removing worst
			final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
			
			final int        debugLevel,
			final int        dbg_X,
			final int        dbg_Y)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int tileSize =      tileProcessor.getTileSize();

		final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int nStiles = stilesX * stilesY; 
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final int debug_stile = (debugLevel > -1)? (dbg_Y * stilesX + dbg_X):-1;
		// TODO: Remove when promoting PlaneData
		final TilePlanes tpl = new TilePlanes(tileSize,superTileSize, geometryCorrection);
		final double [][][][] plane_disp_strength = new double [nStiles][][][];
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
						if (nsTile == debug_stile){
							System.out.println("getPlaneDispStrengths(): nsTile="+nsTile);
						}
						int stileY = nsTile / stilesX;  
						int stileX = nsTile % stilesX;
						int [] sTiles = {stileX, stileY};
						// first make a plane from all tiles
						TilePlanes.PlaneData pd0 = tpl.new  PlaneData (
								sTiles, // int [] sTileXY, 
								tileSize, // int tileSize,
								geometryCorrection, // GeometryCorrection   geometryCorrection,
								correct_distortions,
								measuredLayers,     // MeasuredLayers measuredLayers,
								plPreferDisparity);   // boolean preferDisparity)

						int dl =  (nsTile == debug_stile) ? 3 : 0;

						plane_disp_strength[nsTile] = new double[measuredLayers.getNumLayers()][][];
						
						for (int ml = 0; ml < plane_disp_strength[nsTile].length; ml++) if ((stMeasSel & ( 1 << ml)) != 0){
							
							if (smplMode) {
								plane_disp_strength[nsTile][ml] =  measuredLayers.getDisparityStrength(
										ml,             // int num_layer,
										stileX,         // int stX,
										stileY,         // int stY,
										null,           // boolean [] sel_in,
										strength_floor, // double strength_floor,
										strength_pow,   // double strength_pow,
										smplSide,       // int        smplSide, // = 2;   // Sample size (side of a square)
										smplNum,        //int        smplNum,   // = 3;   // Number after removing worst (should be >1)
										smplRms,        //double     smplRms,   // = 0.1; // Maximal RMS of the remaining tiles in a sample
										true);          // boolean null_if_none);
							} else {
								plane_disp_strength[nsTile][ml] =  measuredLayers.getDisparityStrength(
										ml,             // int num_layer,
										stileX,         // int stX,
										stileY,         // int stY,
										null,           // boolean [] sel_in,
										strength_floor, // double strength_floor,
										strength_pow,   // double strength_pow,
										true);          // boolean null_if_none);
							}
						}
						if (world_plane_norm != null) {
							// find average disparity for the supertile (improve?)
							double sd = 0.0, sw = 0.0;
							for (int ml = 0; ml < plane_disp_strength[nsTile].length; ml++) if ((stMeasSel & ( 1 << ml)) != 0){
								if ((plane_disp_strength[nsTile] !=null) && (plane_disp_strength[nsTile][ml] !=null)  && (plane_disp_strength[nsTile][ml][1] !=null)) {
									for (int i = 0; i < plane_disp_strength[nsTile][ml][1].length; i++){
										double w = plane_disp_strength[nsTile][ml][1][i];
										double d = plane_disp_strength[nsTile][ml][0][i];
										sd += w * d;
										sw += w;
									}
								}
							}
							if (sw > 0) {
								if (dl > 0) {
									System.out.println("Plane tilted disparity for stileX = "+stileX+" stileY="+stileY+", average disparity "+(sd/sw));
								}
								if (dl>1) {
									String [] dbg_titles = showSupertileSeparationTitles( plane_disp_strength[nsTile], null);
									double [][] dbg_img = showSupertileSeparation(false, plane_disp_strength[nsTile], null); // plane_sels);
									showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
									sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "plane_separation_original_disp"+nsTile,dbg_titles);
									dbg_img = showSupertileSeparation(true, plane_disp_strength[nsTile], null); // plane_sels);
									sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "plane_separation_original_world"+nsTile,dbg_titles);
								}
								plane_disp_strength[nsTile] = pd0.getDisparityToPlane(
										world_plane_norm, // double []     world_normal_xyz,
										sd / sw, // average disparity // double        disp_center,
										null, // boolean [][]  tile_sel, // null - do not use, {} use all (will be modified)
										plane_disp_strength[nsTile], // double [][][] disp_str, // calculate just once if null 
										dl); // 1); // int           debugLevel);
								if (dl>1) {
									String [] dbg_titles = showSupertileSeparationTitles( plane_disp_strength[nsTile], null);
									double [][] dbg_img = showSupertileSeparation(false, plane_disp_strength[nsTile], null); // plane_sels);
									showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
									sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "plane_separation_disp"+nsTile,dbg_titles);
									dbg_img = showSupertileSeparation(true, plane_disp_strength[nsTile], null); // plane_sels);
									sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "plane_separation_world"+nsTile,dbg_titles);
								}
							} else {
								plane_disp_strength[nsTile] = null;
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return plane_disp_strength;
	}

	// use histogram data (min/max/strength) assign tiles to clusters
	// returns [supertile] [ plane number] [measurement layer] [tile index]
	public boolean [][][][] dispClusterize(
			final double [][][][] disparity_strengths, // either normal or tilted disparity/strengths
			final double [][][] hist_max_min_max, // histogram data: per tile array of odd number of disparity/strengths pairs, starting with first maximum  
			final boolean [][][] selected,   // tiles OK to be assigned [supertile][measurement layer] [tile index] or null (or null or per-measurement layer)
			final boolean [][][] prohibited, // already assigned tiles [supertile][measurement layer] [tile index] or null
			final int        stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			final double     max_diff,  // maximal disparity difference (to assign to a cluster (of Double.NaN)
			final int        plMinPoints, //          =     5;  // Minimal number of points for plane detection
			final double     smallDiff,  //       = 0.4;   // Consider merging initial planes if disparity difference below
			final double     highMix,    //stHighMix         = 0.4;   // Consider merging initial planes if jumps between ratio above
			final int        debugLevel,
			final int        dbg_X,
			final int        dbg_Y)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();

		final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int nStiles = stilesX * stilesY; 
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		this.planes = new TilePlanes.PlaneData[nStiles][];
		final int debug_stile = (debugLevel > -1)? (dbg_Y * stilesX + dbg_X):-1;

		final boolean [][][][] plane_selections = new boolean [nStiles][][][]; // [supertile] [ plane number] [measurement layer] [tile index]
		final double max_diff2 = Double.isNaN(max_diff)? Double.NaN: (max_diff*max_diff);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
						if (disparity_strengths[nsTile] != null){
							if (nsTile == debug_stile){
								System.out.println("dispClusterize(): nsTile="+nsTile);
							}
							int stileY = nsTile / stilesX;  
							int stileX = nsTile % stilesX;
							int dl =  (nsTile == debug_stile) ? 3 : 0;

							double[][][] disp_strength = new double[measuredLayers.getNumLayers()][][];

							for (int ml = 0; ml < disp_strength.length; ml++) if ((stMeasSel & ( 1 << ml)) != 0){
								disp_strength[ml] = disparity_strengths[nsTile][ml]; // will we change it  - no, no need to clone 
							}

							if (hist_max_min_max[nsTile] == null){
								continue;
							}
							double [][] max_only = new double [(hist_max_min_max[nsTile].length + 1)/2][2];
							for (int i = 0; i < max_only.length; i++){
								max_only[i] = hist_max_min_max[nsTile][2 * i];
							}
							boolean [][][] plane_sels = null;
							int num_ml = disp_strength.length;
							int num_tiles = 4 * superTileSize * superTileSize;
							int [] num_sel;
							boolean [][] tile_en = new boolean [num_ml][];
							for (int ml = 0; ml < num_ml; ml++) if (
									(disp_strength[ml] != null) &&
									((stMeasSel & ( 1 << ml)) != 0) &&
									((selected == null) || ((selected[nsTile] != null) && (selected[nsTile][ml] != null)))
									) {
								tile_en[ml] = new boolean [num_tiles];
								// apply all restrictions here
								for (int i = 0; i < num_tiles; i++){
									if (    (disp_strength[ml][1][i] > 0.0) &&
											((selected == null) || (selected[nsTile][ml].length == 0) || selected[nsTile][ml][i]) &&
											((prohibited == null) || (prohibited[nsTile] == null) || (prohibited[nsTile][ml] == null) || !prohibited[nsTile][ml][i])
											) {
										tile_en[ml][i] = true;
									}
								}
							}

							for (int iter = 0; iter < 2; iter ++){
								int num_p  = max_only.length;
								plane_sels = new boolean[num_p][num_ml][];
								num_sel = new int [num_p];
								for (int np = 0; np < num_p; np++) {
									for (int ml = 0; ml < num_ml; ml++) if (
											(disp_strength[ml] != null) &&
											((stMeasSel & ( 1 << ml)) != 0) &&
											((selected == null) || ((selected[nsTile] != null) && (selected[nsTile][ml] != null)))
											) {
										plane_sels[np][ml] = new boolean[num_tiles];
									}
								}
								// compare closest to be able to use tilted planes later
								for (int ml = 0; ml < num_ml; ml++) if (tile_en[ml] != null) {
									if ( dl >1) {
										System.out.println("dispClusterize(), nsTile="+nsTile+", tile_en["+ml+"] + iter = "+iter);
										for (int i = 0; i < 2 * superTileSize; i ++){
											for (int j = 0; j < 2 * superTileSize; j ++){
												int indx = 2 * superTileSize * i + j;
												System.out.print(tile_en[ml][indx]?" +":" .");
											}
											System.out.println();
										}
									}
									for (int indx = 0; indx < num_tiles; indx++) if (tile_en[ml][indx]){
										int best_plane = -1;
										double best_d2 = Double.NaN;
										for (int np = 0; np < num_p; np++) {
											double d2 = max_only[np][0] - disp_strength[ml][0][indx];
											// add disp_norm correction here?
											d2 *= d2;
											if (!(d2 >= best_d2)){
												best_d2 = d2;
												best_plane = np;
											}
										}
										if (best_plane >= 0){ // compare to max diff here too
											if  (!(best_d2 > max_diff2)) { // works if max_diff2 is Double.NaN
												plane_sels[best_plane][ml][indx] = true; // so far exclusive
											}
										}
									}
									if ( dl >1) {
										System.out.println("dispClusterize(), nsTile="+nsTile+", plane_sels[?]["+ml+"] ,iter = "+iter);
										for (int i = 0; i < 2 * superTileSize; i ++){
											for (int j = 0; j < 2 * superTileSize; j ++){
												int indx = 2 * superTileSize * i + j;
												label1: {
													for (int np = 0; np < num_p; np++) {
														if (plane_sels[np][ml][indx]) {
															System.out.print(" " +np);
															break label1;
														}
													}
													System.out.print(" .");
												}
											}
											System.out.println();
										}
									}
									
									
									
								}
								// recalculate average disparities for each plane and show number of tiles in each in debug mode
								for (int np = 0; np < num_p; np++) {
									double sd = 0.0, sw = 0.0;
									for (int ml = 0; ml < num_ml; ml++) if (disp_strength[ml] != null) {
										for (int indx = 0; indx < num_tiles; indx++) if (plane_sels[np][ml][indx]){
											double w = disp_strength[ml][1][indx];
											sd += w * disp_strength[ml][0][indx];
											sw += w;
											num_sel[np]++;
										}
									}
									if (sw > 0) {
										sd /= sw;
									}
									if (dl > 0) {
										System.out.println("plane num_sel["+np+"] = "+num_sel[np]+" disp "+max_only[np][0]+"->"+sd+
												", weight "+max_only[np][1]+"->"+sw);
									}
									max_only[np][0] = sd;
									max_only[np][1] = sw;
								}
								// calculate transitions matrix (to find candidates for merge
								int [][]    trans_mat =  getTransMatrix(plane_sels);
								double [][] rel_trans =  getTransRel(trans_mat);


								if (dl > 0) {
									System.out.println("trans_mat = ");
									for (int i = 0; i < trans_mat.length; i++){
										System.out.print(i+": ");
										for (int j = 0; j < trans_mat[i].length; j++){
											System.out.print(trans_mat[i][j]+" ");
										}
										System.out.println();
									}
									System.out.println("rel_trans = ");
									for (int i = 0; i < rel_trans.length; i++){
										System.out.print(i+": ");
										for (int j = 0; j < rel_trans[i].length; j++){
											System.out.print(rel_trans[i][j]+" ");
										}
										System.out.println();
									}
								}
								if ((iter > 0 ) && (num_p > 1)){ // remove /join bad 
									int windx = 0;
									int remove_indx = -1;
									for (int i = 1; i < num_p; i++)	if (num_sel[i] < num_sel[windx]) windx = i;
									if (num_sel[windx] < plMinPoints) {
										if ((debugLevel > 2) || (dl > 0)){
											System.out.println ("dispClusterize(): stileX = "+stileX+" stileY="+stileY+
													": removing plane "+windx+" with "+num_sel[windx]+" tiles ( <"+plMinPoints+")");
										}
										remove_indx = windx;
									}
									if (remove_indx < 0) {
										// find candidates for merge
										windx = -1; 
										for (int i = 0; i < (num_p - 1); i++)	{
											if (((max_only[i+1][0] - max_only[i][0]) < smallDiff) &&  // close enough to consider merging
													(rel_trans[i][i+1] > highMix)) {
												if ((windx < 0) || (rel_trans[i][i+1] > rel_trans[windx][windx+1])) windx = i;
											}
										}
										if (windx >=0 ) {
											if ((debugLevel > 2) || (dl > 0)) {
												System.out.println ("dispClusterize(): stileX = "+stileX+" stileY="+stileY+
														": merging plane "+windx+" with " + (windx + 1)+": "+
														num_sel[windx] + " and "+num_sel[windx+1]+" tiles, "+
														" rel_trans="+rel_trans[windx][windx + 1]+ " ( > " + highMix+"),"+
														" diff="+ (max_only[windx + 1][0]- max_only[windx][0]) + " ( < " + smallDiff+" ),"+
														" disp1 = "+max_only[windx][0]+" disp2 = "+max_only[windx + 1][0]);
											}
											double sum_w = max_only[windx][1] + max_only[windx + 1][1];
											max_only[windx+1][0] = (max_only[windx][0]*max_only[windx][1] +  max_only[windx+1][0]*max_only[windx+1][1]) / sum_w;
											max_only[windx+1][1] = sum_w;
											remove_indx = windx;
										}
									}
									if (remove_indx >= 0){
										double [][] max_only_copy = max_only.clone();
										for (int i = 0; i < max_only.length; i++) max_only_copy[i] = max_only[i];
										max_only = new double [max_only.length - 1][];
										int indx = 0;
										for (int i = 0; i < max_only_copy.length; i++) if (i != remove_indx) max_only[indx++] =max_only_copy[i]; 
										iter = 0;
										continue; // restart from 0

									}

									// Show other candidates for merge
									if ((debugLevel > 2) || (dl > 0)){
										double max_sep = 0.2;
										if (iter  > 0) {
											for (int i = 0; i < (num_p-1); i++){
												if (rel_trans[i][i+1] > max_sep) {
													System.out.println("dispClusterize() stileX = "+stileX+" stileY="+stileY+" lowplane = "+i+
															" num_sel1 = "+num_sel[i] + " num_sel2 = "+num_sel[i+1] +
															" rel_trans="+rel_trans[i][i+1]+
															" diff="+ (max_only[i+1][0]- max_only[i][0]) +
															" disp1 = "+max_only[i][0]+" disp2 = "+max_only[i+1][0]);
												}
											}
										}
									}								
								}
							}

							if (dl > 3) {
								String [] dbg_titles = showSupertileSeparationTitles( disp_strength, plane_sels);
								double [][] dbg_img = showSupertileSeparation(false,disp_strength, plane_sels);
								showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
								sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "disp_clusterize_disp"+nsTile+"-"+debugLevel,dbg_titles);
								dbg_img = showSupertileSeparation(true, disp_strength, plane_sels);
								sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "disp_clusterize_world"+nsTile+"-"+debugLevel,dbg_titles);
							}
							plane_selections[nsTile] = plane_sels;
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return plane_selections;
	}

	// use both horizontal and const disparity tiles to create tile clusters
	// Add max_diff (maximal disparity difference while extracting initial tile selection) and max_tries (2..3) parameters
	
	// Add separate method to create + remove outliers from all planes (2 different ones)?
	// TODO later re-assign pixels according to existing plane parameters
	// Sort plane data by center (plane or supertile) disparity
	
	public boolean [][][][]  initialDiscriminateTiles(
			final int        growSelection,                     // grow initial selection before processing 
			final int        stMeasSel,       //      = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			final double     plDispNorm,
			final int        plMinPoints, //          =     5;  // Minimal number of points for plane detection
			final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
			final GeometryCorrection geometryCorrection,
			final boolean    correct_distortions,

			final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
			final int        smplSide, //        = 2;      // Sample size (side of a square)
			final int        smplNum, //         = 3;      // Number after removing worst
			final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample

			final double     bin_blur_hor,   // Blur disparity histograms for horizontal clusters by this sigma (in bins)
			final double     bin_blur_vert,  // Blur disparity histograms for constant disparity clusters by this sigma (in bins)
			final double     max_diff_hor,   // maximal disparity difference (to assign to a cluster (of Double.NaN) at first run for horizontal planes
			final double     max_diff_vert,  // maximal disparity difference (to assign to a cluster (of Double.NaN) at first run for vertical plane

			final int        max_tries,       // on last run - assign all rfemaining pixels to some cluster (disregard max_diff)
			final double     smallDiff,  //       = 0.4;   // Consider merging initial planes if disparity difference below
			final double     highMix,    //stHighMix         = 0.4;   // Consider merging initial planes if jumps between ratio above
			final double []  world_hor, // horizontal plane normal (default [0.0, 1.0, 0.0])
			final int        debugLevel,
			final int        dbg_X,
			final int        dbg_Y)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
//		final int tileSize =      tileProcessor.getTileSize();

		final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int nStiles = stilesX * stilesY; 
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		this.planes = new TilePlanes.PlaneData[nStiles][];
		final int debug_stile = (debugLevel > -1)? (dbg_Y * stilesX + dbg_X):-1;
		// TODO: Remove when promoting PlaneData
//		final TilePlanes tpl = new TilePlanes(tileSize,superTileSize, geometryCorrection);
		final int num_layers = measuredLayers.getNumLayers();
		final int num_tiles = 4 * superTileSize * superTileSize;
		boolean [] grown_selection = null;
		if (growSelection >= 0) {
			grown_selection = cltPass3d.getSelected();
			if (growSelection > 0) {
				grown_selection = grown_selection.clone();
				measuredLayers.growSelection(
						growSelection, // grow,
						grown_selection, // tiles,
						null); // prohibit);
			}
		}
		measuredLayers.setLayer ( // overwrite?
				0, // int       num_layer,
				cltPass3d.getDisparity(), // double [] disparity,
				cltPass3d.getStrength(), // double [] strength,
				grown_selection); // null); // boolean [] selection) // may be null
		if (debugLevel > -1) {
			String [] titles = {"d0","s0","d1","s1","d2","s2","d3","s3","s","d","selection"};
			boolean [] dbg_sel= grown_selection; // cltPass3d.getSelected(); 
			double [][] dbg_img = new double [titles.length][];
			for (int i = 0; i < measuredLayers.getNumLayers(); i++){
				dbg_img[2 * i] =     measuredLayers.getDisparity(i);
				dbg_img[2 * i + 1] = measuredLayers.getStrength(i);
			}
			dbg_img[8] = cltPass3d.getDisparity();
			dbg_img[9] = cltPass3d.getStrength();
			dbg_img[10]= new double [dbg_sel.length];
			for (int i = 0; i < dbg_sel.length; i++){
				dbg_img[10][i] = ((dbg_sel == null) || dbg_sel[i])? 1.0:0.0;
			}
			showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
			sdfa_instance.showArrays(dbg_img,  tileProcessor.getTilesX(), tileProcessor.getTilesY(), true, "measuredLayers",titles);
		}
		
//		double [] world_hor = {0.0, 1.0, 0.0};

		final double [][][][] hor_disp_strength = getPlaneDispStrengths(
				world_hor,           // final double []  world_plane_norm, // real world normal vector to a suggested plane family (0,1,0) for horizontal planes
				stMeasSel,           //final int        stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
				plPreferDisparity,   // final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
				geometryCorrection,  // final GeometryCorrection geometryCorrection,
				correct_distortions, // final boolean    correct_distortions,
				smplMode,            // final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				smplSide,            //final int        smplSide, //        = 2;      // Sample size (side of a square)
				smplNum,             //final int        smplNum, //         = 3;      // Number after removing worst
				smplRms,             //final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample

				debugLevel,
				dbg_X,
				dbg_Y);

		final double [][][][] vert_disp_strength = getPlaneDispStrengths(
				null,           // final double []  world_plane_norm, // real world normal vector to a suggested plane family (0,1,0) for horizontal planes
				stMeasSel,           //final int        stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
				plPreferDisparity,   // final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
				geometryCorrection,  // final GeometryCorrection geometryCorrection,
				correct_distortions, // final boolean    correct_distortions,
				smplMode,            // final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				smplSide,            //final int        smplSide, //        = 2;      // Sample size (side of a square)
				smplNum,             //final int        smplNum, //         = 3;      // Number after removing worst
				smplRms,             //final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
				debugLevel,
				dbg_X,
				dbg_Y);


		String [] dbg_hist_titles = {"all","hor","mm_vert","mm_hor"};
		double [][] dbg_hist = new double [dbg_hist_titles.length][];
		final boolean [][][] used = new boolean [nStiles][num_layers][]; // num_tiles
		final boolean [][][][] planes_selections = new boolean [nStiles][][][]; // num_tiles
		for (int pass = 0; pass < max_tries; pass ++) {
//			final int fpass = pass;
			// get both horizontal planes and constant disparity planes histograms, 
//			resetDisparityHistograms();
			setBlurSigma(bin_blur_hor);
			final double [][][] mmm_hor = getMaxMinMax(
					hor_disp_strength,  // final double [][][][] disparity_strength, // pre-calculated disparity/strength [per super-tile][per-measurement layer][2][tiles] or null
					null); // final boolean [][] tile_sel // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all
			if (debugLevel > -1) {
				dbg_hist[1] = showDisparityHistogram();
				dbg_hist[3] = showMaxMinMax();
			}
//			resetDisparityHistograms();
			setBlurSigma(bin_blur_vert);
			final double [][][] mmm_vert = getMaxMinMax(
					vert_disp_strength,  // final double [][][][] disparity_strength, // pre-calculated disparity/strength [per super-tile][per-measurement layer][2][tiles] or null
					null); // final boolean [][] tile_sel // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all

			if (debugLevel > 0) {
				dbg_hist[0] = showDisparityHistogram();
				dbg_hist[2] = showMaxMinMax();
			}
			if (debugLevel > 0) {
				int hist_width0 =  showDisparityHistogramWidth();
				int hist_height0 = dbg_hist[0].length/hist_width0;
				showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
				sdfa_instance.showArrays(dbg_hist, hist_width0, hist_height0, true, "vert_hor_histograms_"+pass,dbg_hist_titles);
			}
			// try to independently (same selections) clusterize both ways
			if (debugLevel > -1){
				System.out.println("initialDiscriminateTiles(): before new_planes_hor, pass =" + (pass + 1) + " ( of "+max_tries+" )");
			}
			
			final boolean [][][][] new_planes_hor = dispClusterize(
					hor_disp_strength, // final double [][][][] disparity_strengths, // either normal or tilted disparity/strengths
					mmm_hor,           // final double [][][] hist_max_min_max, // histogram data: per tile array of odd number of disparity/strengths pairs, starting with first maximum  
					null,              // final boolean [][][] selected,   // tiles OK to be assigned [supertile][measurement layer] [tile index] or null (or null or per-measurement layer)
					used,              // final boolean [][][] prohibited, // already assigned tiles [supertile][measurement layer] [tile index] or null
					stMeasSel,         // final int        stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
					((pass < (max_tries - 1)) ? max_diff_hor : Double.NaN), // final double     max_diff,  // maximal disparity difference (to assign to a cluster (of Double.NaN)
					plMinPoints,       // final int        plMinPoints, //          =     5;  // Minimal number of points for plane detection
					smallDiff,         // final double     smallDiff,  //       = 0.4;   // Consider merging initial planes if disparity difference below
					highMix,          // final double     highMix,    //stHighMix         = 0.4;   // Consider merging initial planes if jumps between ratio above
					1, // debugLevel,
					dbg_X,
					dbg_Y);

			if (debugLevel > -1){
				System.out.println("initialDiscriminateTiles(): before new_planes_vert, pass =" + (pass + 1) + " ( of "+max_tries+" )");
			}
			final boolean [][][][] new_planes_vert = dispClusterize(
					vert_disp_strength, // final double [][][][] disparity_strengths, // either normal or tilted disparity/strengths
					mmm_vert,           // final double [][][] hist_max_min_max, // histogram data: per tile array of odd number of disparity/strengths pairs, starting with first maximum  
					null,              // final boolean [][][] selected,   // tiles OK to be assigned [supertile][measurement layer] [tile index] or null (or null or per-measurement layer)
					used,              // final boolean [][][] prohibited, // already assigned tiles [supertile][measurement layer] [tile index] or null
					stMeasSel,         // final int        stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
					((pass < (max_tries - 1)) ? max_diff_vert : Double.NaN), // final double     max_diff,  // maximal disparity difference (to assign to a cluster (of Double.NaN)
					plMinPoints,       // final int        plMinPoints, //          =     5;  // Minimal number of points for plane detection
					smallDiff,         // final double     smallDiff,  //       = 0.4;   // Consider merging initial planes if disparity difference below
					highMix,          // final double     highMix,    //stHighMix         = 0.4;   // Consider merging initial planes if jumps between ratio above
					2, // debugLevel,
					dbg_X,
					dbg_Y);
			
			// compare which vert or hor provide stronger clusters, add them to selections, recalculate "used"
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
							int dl =  ((debugLevel > -1) && (nsTile == debug_stile)) ? 3 : 0;
							if (dl > 0){
								System.out.println("initialDiscriminateTiles() selecting: nsTile="+nsTile);
							}
							int stileY = nsTile / stilesX;  
							int stileX = nsTile % stilesX;
//							int [] sTiles = {stileX, stileY};
							double [][][][] ds = {vert_disp_strength[nsTile],hor_disp_strength[nsTile]};
							boolean [][][][] sels_all = {new_planes_vert[nsTile],new_planes_hor[nsTile]}; // make possible to iterate
							class SelStrength{
								int type;
								int indx;
								double strength;
								SelStrength(int type, int indx, double strength){
									this.type = type;
									this.indx = indx;
									this.strength = strength;
								}
							}
							ArrayList<SelStrength> selStrengthList = new ArrayList <SelStrength>();
							for (int pType = 0; pType < sels_all.length; pType++){
								if (sels_all[pType] == null) sels_all[pType] = new boolean[0][][];
								// calculate strength of each selection;
								for (int np = 0; np < sels_all[pType].length; np ++) {
									double sw = 0.0;
									for (int ml = 0; ml < num_layers; ml++) if ((sels_all[pType][np] != null) && (sels_all[pType][np][ml] != null)){
										for (int i = 0; i < num_tiles; i++) {
											if (sels_all[pType][np][ml][i]) { 
												sw +=  ds[pType][ml][1][i];
											}
										}
									}
									selStrengthList.add(new SelStrength(pType, np, sw));
								}
							}
							if (dl > 0){
								System.out.println("initialDiscriminateTiles() got list of clusters for "+nsTile);
								for (int i = 0; i < selStrengthList.size(); i++){
									System.out.println(i+": type = "+selStrengthList.get(i).type+
											" index = "+selStrengthList.get(i).indx +
											" strength = "+selStrengthList.get(i).strength);
								}
							}
							if (selStrengthList.size() > 0) {
								Collections.sort(selStrengthList, new Comparator<SelStrength>() {
									@Override
									public int compare(SelStrength lhs, SelStrength rhs) {
										// -1 - less than, 1 - greater than, 0 - equal, all inverted for descending
										return (lhs.strength > rhs.strength) ? -1 : (lhs.strength < rhs.strength ) ? 1 : 0;
									}
								});
								if (dl > 0){
									System.out.println("initialDiscriminateTiles() sorted list of clusters for "+nsTile);
									for (int i = 0; i < selStrengthList.size(); i++){
										System.out.println(i+": type = "+selStrengthList.get(i).type+
												" index = "+selStrengthList.get(i).indx +
												" strength = "+selStrengthList.get(i).strength);
									}
								}
								
								int best_type = selStrengthList.get(0).type;
								int num_planes = 0;
								for (; num_planes < selStrengthList.size(); num_planes++) {
									if (selStrengthList.get(num_planes).type != best_type) {
										break;
									}
								}
								int num_old_planes = (planes_selections[nsTile] == null) ? 0: planes_selections[nsTile].length;
								if (dl > 0){
									System.out.println("initialDiscriminateTiles() num_old_planes= "+num_old_planes);
								}
								boolean [][][] new_planes_selections = new boolean [num_old_planes + num_planes][][];
								int np = 0;
								for (; np < num_old_planes; np++){
									new_planes_selections[np] = planes_selections[nsTile][np];
								}
								for (int nnp = 0; nnp < num_planes; nnp++) { // SelStrength ss:selStrengthList){
									SelStrength ss = selStrengthList.get(nnp);
									new_planes_selections[np++] = sels_all[ss.type][ss.indx];
									for (int ml = 0; ml < num_layers; ml++) if (sels_all[ss.type][ss.indx][ml] != null){
										if (used[nsTile][ml] == null){
											used[nsTile][ml] = new boolean[num_tiles];
										}
										for (int i = 0; i < num_tiles; i++){
											used[nsTile][ml][i] |= sels_all[ss.type][ss.indx][ml][i];
										}
									}
								}
								planes_selections[nsTile] = new_planes_selections;
								if (dl > 0){
									System.out.println("initialDiscriminateTiles() new_planes_selections.length= "+new_planes_selections.length);
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}
		return planes_selections;
	}

	public TilePlanes.PlaneData [][] createPlanesFromSelections(
			final boolean [][][][] plane_selections, //  = new boolean [nStiles][][][]; // num_tiles
			final double  [][][][] disp_strength,			
			final double     plDispNorm,
			final int        plMinPoints, //          =     5;  // Minimal number of points for plane detection
			final double     plTargetEigen, //        =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
			final double     plFractOutliers, //      =   0.3;  // Maximal fraction of outliers to remove
			final int        plMaxOutliers, //        =    20;  // Maximal number of outliers to remove
//			final double     plVertWors,    //        =    1.5  // if rotating plane vertical does not increase 'eigenvalue' more, use vertical  
			final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
			final GeometryCorrection geometryCorrection,
			final boolean    correct_distortions,

			final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
			final int        smplSide, //        = 2;      // Sample size (side of a square)
			final int        smplNum, //         = 3;      // Number after removing worst
			final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample

			final int        debugLevel,
			final int        dbg_X,
			final int        dbg_Y)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int tileSize =      tileProcessor.getTileSize();

		final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int nStiles = stilesX * stilesY; 
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final TilePlanes.PlaneData [][] result_planes = new TilePlanes.PlaneData[nStiles][];
		//		this.planes = new TilePlanes.PlaneData[nStiles][];
		final int debug_stile = (debugLevel > -1)? (dbg_Y * stilesX + dbg_X):-1;
		// TODO: Remove when promoting PlaneData
		final TilePlanes tpl = new TilePlanes(tileSize,superTileSize, geometryCorrection);

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
						if (nsTile == debug_stile){
							System.out.println("createPlanesFromSelections(): nsTile="+nsTile);
						}
						if (plane_selections[nsTile] != null) {
							int stileY = nsTile / stilesX;  
							int stileX = nsTile % stilesX;
							int [] sTiles = {stileX, stileY};
							int dl =  (nsTile == debug_stile) ? 3 : 0;
							result_planes[nsTile] = null;
							// first make a plane from all tiles
							ArrayList<TilePlanes.PlaneData> st_planes = new ArrayList<TilePlanes.PlaneData>();
							TilePlanes.PlaneData pd0 = tpl.new  PlaneData (
									sTiles, // int [] sTileXY, 
									tileSize, // int tileSize,
									geometryCorrection, // GeometryCorrection   geometryCorrection,
									correct_distortions,
									measuredLayers,     // MeasuredLayers measuredLayers,
									plPreferDisparity);   // boolean preferDisparity)
							// iterate through all plane selections
							for (int ps = 0; ps < plane_selections[nsTile].length; ps++) {
								TilePlanes.PlaneData pd = pd0.clone(); 
								boolean OK = (pd.getPlaneFromMeas(
										plane_selections[nsTile][ps], // tile_sel,       // boolean [][] tile_sel, // null - do not use, {} use all (will be modified)
										disp_strength[nsTile],
										Double.NaN,    // double       disp_far, // minimal disparity to select (or NaN)
										Double.NaN,    // double       disp_near, // maximal disparity to select (or NaN)
										plDispNorm,    // 0.0,            // plDispNorm,  // double       dispNorm,   //  Normalize disparities to the average if above
										0.0,            // double       min_weight,
										plMinPoints,    // int          min_tiles,
										strength_floor, // 
										strength_pow,   // double       strength_pow,
										// update !
										smplMode,
										smplSide,
										smplNum,
										smplRms,
										dl) != null);            // int          debugLevel)
								if (OK) {
									if (dl > 0) {
										if (pd.getWeight() > 1.0) {
											System.out.println("Processing subplane["+nsTile+"]["+ps+"]"+
													", stileX="+stileX+
													", stileY="+stileY+
													", numPoints="+ pd.getNumPoints()+
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
										OK = pd.removeOutliers( // getPlaneFromMeas should already have run
												disp_strength[nsTile], 
												targetV,      // double     targetEigen, // target eigenvalue for primary axis (is disparity-dependent, so is non-constant)
												max_outliers, // int        maxRemoved,  // maximal number of tiles to remove (not a constant)
												dl); // int        debugLevel)
										if (!OK) {
											continue;
										}
										if (dl > 0) {
											if (pd.getWeight() > 1.0) {
												System.out.println("Removed outliers["+nsTile+"]["+ps+"]"+
														", stileX="+stileX+
														", stileY="+stileY+
														", numPoints="+ pd.getNumPoints()+
														", swc = "+pd.getWeight()+
														", center=["+pd.getZxy()[0]+","+pd.getZxy()[1]+","+pd.getZxy()[2]+"]"+
														", eig_val = {"+pd.getValues()[0]+","+pd.getValues()[1]+","+pd.getValues()[2]+"}"+
														", eig_vect[0] = {"+pd.getVector()[0]+","+pd.getVector()[1]+","+pd.getVector()[2]+"}");
											}
										}
									}
									double [] norm_xyz = pd.getWorldXYZ(
											correct_distortions);
									st_planes.add(pd);
									if (dl > 0) {
										System.out.println("World normal["+nsTile+"]["+ps+"] = {"+
												norm_xyz[0]+", "+norm_xyz[1]+", "+norm_xyz[2]+"}");

									}
								}

							}
							if (st_planes.size() > 0){
								// sort planes by increasing disparity (tile center or plane center ? ) Using plane center
								Collections.sort(st_planes, new Comparator<TilePlanes.PlaneData>() {
									@Override
									public int compare(TilePlanes.PlaneData lhs, TilePlanes.PlaneData rhs) {
										// -1 - less than, 1 - greater than, 0 - equal
										return (rhs.getZxy()[0] > lhs.getZxy()[0]) ? -1 : (rhs.getZxy()[0] < lhs.getZxy()[0] ) ? 1 : 0;
									}
								});
								if (LOWEST_PLANE(2) > 0) st_planes.add(0, st_planes.get(0)); // insert dummy at pos 0;
								result_planes[nsTile] = st_planes.toArray(new TilePlanes.PlaneData[0] );
								if (LOWEST_PLANE(2) > 0) result_planes[nsTile][0] = null; // remove dummy
								if (dl >0){
									System.out.println("createPlanesFromSelections(): nsTile="+nsTile);
								}
							}
							if (dl > 2) {
								String [] dbg_titles = showSupertileSeparationTitles( disp_strength[nsTile], plane_selections[nsTile], result_planes[nsTile]);
								double [][] dbg_img =  showSupertileSeparation(false, disp_strength[nsTile], plane_selections[nsTile], result_planes[nsTile]);
								showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
								sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "create_planes_disp-"+nsTile+"-"+debugLevel,dbg_titles);
								dbg_img =  showSupertileSeparation(true, disp_strength[nsTile], plane_selections[nsTile], result_planes[nsTile]);
								sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "create_planes_world-"+nsTile+"-"+debugLevel,dbg_titles);
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return result_planes;		
	}

	
	public void processPlanes5(
			final int        growSelection,                     // grow initial selection before processing 
			final int        stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			final double     plDispNorm,
			final int        plMinPoints, //          =     5;  // Minimal number of points for plane detection
			final double     plTargetEigen, //        =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
			final double     plFractOutliers, //      =   0.3;  // Maximal fraction of outliers to remove
			final int        plMaxOutliers, //        =    20;  // Maximal number of outliers to remove
			final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
			final GeometryCorrection geometryCorrection,
			final boolean    correct_distortions,

			final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
			final int        smplSide, //        = 2;      // Sample size (side of a square)
			final int        smplNum, //         = 3;      // Number after removing worst
			final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample

			final double     bin_blur_hor,   // Blur disparity histograms for horizontal clusters by this sigma (in bins)
			final double     bin_blur_vert,  // Blur disparity histograms for constant disparity clusters by this sigma (in bins)
			final double     max_diff_hor,   // maximal disparity difference (to assign to a cluster (of Double.NaN) at first run for horizontal planes
			final double     max_diff_vert,  // maximal disparity difference (to assign to a cluster (of Double.NaN) at first run for vertical plane
			final int        max_tries,       // on last run - assign all remaining pixels to some cluster (disregard max_diff)

			final boolean    msUseSel,        // final boolean                   use_sel,
			final boolean    msDivideByArea,  // final boolean                   divide_by_area,
			final double     msScaleProj,     // final double                    scale_projection,
			
			final double     smallDiff,  //       = 0.4;   // Consider merging initial planes if disparity difference below
			final double     highMix,    //stHighMix         = 0.4;   // Consider merging initial planes if jumps between ratio above
			final double []  world_hor, // horizontal plane normal (default [0.0, 1.0, 0.0])
			final int        debugLevel,
			final int        dbg_X,
			final int        dbg_Y)
	{
		// use both horizontal and const disparity tiles to create tile clusters
		// Add max_diff (maximal disparity difference while extracting initial tile selection) and max_tries (2..3) parameters
		
		// Add separate method to create + remove outliers from all planes (2 different ones)?
		// TODO later re-assign pixels according to existing plane parameters
		// Sort plane data by center (plane or supertile) disparity
		
		boolean [][][][]  plane_selections = initialDiscriminateTiles(
				growSelection,        // final int        growSelection,                     // grow initial selection before processing 
				stMeasSel,            // final int        stMeasSel,       //      = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
				plDispNorm,           // final double     plDispNorm,
				plMinPoints,          // final int        plMinPoints, //          =     5;  // Minimal number of points for plane detection
				plPreferDisparity,    // final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
				geometryCorrection,  // final GeometryCorrection geometryCorrection,
				correct_distortions, // final boolean    correct_distortions,

				smplMode,            // final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				smplSide,            //  final int        smplSide, //        = 2;      // Sample size (side of a square)
				smplNum, //  final int        smplNum, //         = 3;      // Number after removing worst
				smplRms, //  final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample

				bin_blur_hor,   // final double     bin_blur_hor,   // Blur disparity histograms for horizontal clusters by this sigma (in bins)
				bin_blur_vert,  // final double     bin_blur_vert,  // Blur disparity histograms for constant disparity clusters by this sigma (in bins)
				max_diff_hor,   // final double     max_diff_hor,   // maximal disparity difference (to assign to a cluster (of Double.NaN) at first run for horizontal planes
				max_diff_vert,  // final double     max_diff_vert,  // maximal disparity difference (to assign to a cluster (of Double.NaN) at first run for vertical plane
				max_tries,       //final int        max_tries,       // on last run - assign all rfemaining pixels to some cluster (disregard max_diff)
				smallDiff,      // final double     smallDiff,  //       = 0.4;   // Consider merging initial planes if disparity difference below
				highMix,    //final double     highMix,    // stHighMix         = 0.4;   // Consider merging initial planes if jumps between ratio above
				world_hor, // final double []  world_hor, // horizontal plane normal (default [0.0, 1.0, 0.0])
				debugLevel, // final int        debugLevel,
				dbg_X, // final int        dbg_X,
				dbg_Y); // final int        dbg_Y)

		// get per-tile disparity strength again (may consider using non-filtered data here)
		
		double [][][][] disp_strength = getPlaneDispStrengths( // here use actual disparity, not tilted
				null,                // final double []  world_plane_norm, // real world normal vector to a suggested plane family (0,1,0) for horizontal planes
				stMeasSel,           // final int        stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert

				plPreferDisparity,   // final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
				geometryCorrection,  // final GeometryCorrection geometryCorrection,
				correct_distortions, // final boolean    correct_distortions,

				smplMode,            // final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				smplSide,            // final int        smplSide, //        = 2;      // Sample size (side of a square)
				smplNum,             // final int        smplNum, //         = 3;      // Number after removing worst
				smplRms,             // final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
				
				debugLevel,          // final int        debugLevel,
				dbg_X,               // final int        dbg_X,
				dbg_Y);              // final int        dbg_Y)
        // Create plane data (ellipsoids) from the disparity/strength data and tile selections, 
		// remove outliers, order each supertile result in ascending disparity (from far to near) 		
		// TODO: consider 1) assigning of the non-assigned tiles to clusters and 2) re-assigning all clusters one by one to the "best" plane
		TilePlanes.PlaneData [][] new_planes = createPlanesFromSelections(
				plane_selections,   // final boolean [][][][] plane_selections, //  = new boolean [nStiles][][][]; // num_tiles
				disp_strength,       // final double  [][][][] disp_strength,			
				plDispNorm,          // final double     plDispNorm,
				plMinPoints,         // final int        plMinPoints, //          =     5;  // Minimal number of points for plane detection
				plTargetEigen,       // final double     plTargetEigen, //        =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
				plFractOutliers,     // final double     plFractOutliers, //      =   0.3;  // Maximal fraction of outliers to remove
				plMaxOutliers,       // final int        plMaxOutliers, //        =    20;  // Maximal number of outliers to remove
				plPreferDisparity,   // final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
				geometryCorrection,  // final GeometryCorrection geometryCorrection,
				correct_distortions, // final boolean    correct_distortions,

				smplMode,            // final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				smplSide,            // final int        smplSide, //        = 2;      // Sample size (side of a square)
				smplNum,             // final int        smplNum, //         = 3;      // Number after removing worst
				smplRms,             // final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample

				debugLevel,          // final int        debugLevel,
				dbg_X,               // final int        dbg_X,
				dbg_Y);              // final int        dbg_Y)
		this.planes = new_planes; // save as "measured" (as opposed to "smoothed" by neighbors) planes
		
	}
	
	
	public void processPlanes4(
//			final boolean [] selected, // or null
//			final double     min_disp,
			final int        stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			final double     plDispNorm,
			final int        plMinPoints, //          =     5;  // Minimal number of points for plane detection
			final double     plTargetEigen, //        =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
			final double     plFractOutliers, //      =   0.3;  // Maximal fraction of outliers to remove
			final int        plMaxOutliers, //        =    20;  // Maximal number of outliers to remove
			final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
			final GeometryCorrection geometryCorrection,
			final boolean    correct_distortions,

			final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
			final int        smplSide, //        = 2;      // Sample size (side of a square)
			final int        smplNum, //         = 3;      // Number after removing worst
			final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample

			final double     bin_blur_hor,   // Blur disparity histograms for horizontal clusters by this sigma (in bins)
			final double     bin_blur_vert,  // Blur disparity histograms for constant disparity clusters by this sigma (in bins)

			final double     smallDiff,  //       = 0.4;   // Consider merging initial planes if disparity difference below
			final double     highMix,    //stHighMix         = 0.4;   // Consider merging initial planes if jumps between ratio above
//			final double []  vertical_xyz, // real world up unit vector in camera CS (x - right, y - up, z - to camera};
			final int        debugLevel,
			final int        dbg_X,
			final int        dbg_Y)
	{
		if (maxMinMax == null) getMaxMinMax(null, null); // so far - no planes, no selection
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int tileSize =      tileProcessor.getTileSize();

		final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int nStiles = stilesX * stilesY; 
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		this.planes = new TilePlanes.PlaneData[nStiles][];
		final int debug_stile = (debugLevel > -1)? (dbg_Y * stilesX + dbg_X):-1;
		
//		final boolean [][] dflt_select = {{}, null, null, null, null}; // use layer 0 (combo) only
/*		
		final boolean [][] dflt_select = new boolean [measuredLayers.getNumLayers()][];
		for (int i = 0; i < dflt_select.length; i++){
			if ((stMeasSel & (1 << i)) !=0){
				dflt_select[i] = new boolean[0];
			} else {
				dflt_select[i] = null;
			}
		}
*/
		// TODO: Remove when promoting PlaneData
		final TilePlanes tpl = new TilePlanes(tileSize,superTileSize, geometryCorrection);
		
//		final double  [] disparity = cltPass3d.getDisparity();
//		final double  [] strength =  cltPass3d.getStrength();
		
		measuredLayers.setLayer (
				0, // int       num_layer,
				cltPass3d.getDisparity(), // double [] disparity,
				cltPass3d.getStrength(), // double [] strength,
				null); // boolean [] selection) // may be null
		if (debugLevel > -1) {
			String [] titles = {"d0","s0","d1","s1","d2","s2","d3","s3","s","d"};
			double [][] dbg_img = new double [titles.length][];
			for (int i = 0; i < measuredLayers.getNumLayers(); i++){
				dbg_img[2 * i] =     measuredLayers.getDisparity(i);
				dbg_img[2 * i + 1] = measuredLayers.getStrength(i);
			}
			dbg_img[8] = cltPass3d.getDisparity();
			dbg_img[9] = cltPass3d.getStrength();
			showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
			sdfa_instance.showArrays(dbg_img,  tileProcessor.getTilesX(), tileProcessor.getTilesY(), true, "measuredLayers",titles);
		}
//		getMaxMinMax(
//		null,  // final double [][][][] disparity_strength, // pre-calculated disparity/strength [per super-tile][per-measurement layer][2][tiles] or null
//		null); // final boolean [][] tile_sel // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all
		
		double [] world_hor = {0.0, 1.0, 0.0};

		final double [][][][] plane_disp_strength = getPlaneDispStrengths(
				world_hor,           // final double []  world_plane_norm, // real world normal vector to a suggested plane family (0,1,0) for horizontal planes
				stMeasSel,           //final int        stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
				plPreferDisparity,   // final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
				geometryCorrection,  // final GeometryCorrection geometryCorrection,
				correct_distortions, // final boolean    correct_distortions,
				smplMode,            // final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				smplSide,            //final int        smplSide, //        = 2;      // Sample size (side of a square)
				smplNum,             //final int        smplNum, //         = 3;      // Number after removing worst
				smplRms,             //final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample

				debugLevel,
				dbg_X,
				dbg_Y);

		String [] dbg_hist_titles = {"all","hor","mm_all","mm_hor"};
		double [][] dbg_hist = new double [dbg_hist_titles.length][];
		
		System.out.println("Calculating histograms for hoirizontal planes");
//		resetDisparityHistograms();
		setBlurSigma(bin_blur_hor);

		final double [][][] mmm_hor = getMaxMinMax(
				plane_disp_strength,  // final double [][][][] disparity_strength, // pre-calculated disparity/strength [per super-tile][per-measurement layer][2][tiles] or null
				null); // final boolean [][] tile_sel // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all
		
		if (debugLevel > -1) {
			dbg_hist[1] = showDisparityHistogram().clone();
			dbg_hist[3] = showMaxMinMax().clone();
		}
//		resetDisparityHistograms();
		setBlurSigma(bin_blur_vert);
		final double [][][] mmm_all = getMaxMinMax(
				null,  // final double [][][][] disparity_strength, // pre-calculated disparity/strength [per super-tile][per-measurement layer][2][tiles] or null
				null); // final boolean [][] tile_sel // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all

		if (debugLevel > -1) {
			dbg_hist[0] = showDisparityHistogram().clone();
			dbg_hist[2] = showMaxMinMax().clone();
		}
		if (debugLevel > -1) {
			int hist_width0 =  showDisparityHistogramWidth();
			int hist_height0 = dbg_hist[0].length/hist_width0;
			showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
			sdfa_instance.showArrays(dbg_hist, hist_width0, hist_height0, true, "all_hor_histograms",dbg_hist_titles);
		}

		
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
						if (nsTile == debug_stile){
							System.out.println("processPlanes4(): nsTile="+nsTile);
						}
						int stileY = nsTile / stilesX;  
						int stileX = nsTile % stilesX;
						int [] sTiles = {stileX, stileY};
						planes[nsTile] = null;
						// first make a plane from all tiles
						TilePlanes.PlaneData pd0 = tpl.new  PlaneData (
								sTiles, // int [] sTileXY, 
								tileSize, // int tileSize,
								geometryCorrection, // GeometryCorrection   geometryCorrection,
								correct_distortions,
								measuredLayers,     // MeasuredLayers measuredLayers,
								plPreferDisparity);   // boolean preferDisparity)
/*
						boolean [][] tile_sel = dflt_select.clone();
						for (int i = 0; i < dflt_select.length; i++){
							if (dflt_select[i] != null) tile_sel[i] = dflt_select[i].clone();
						}
*/
//						int dl1 =  (nsTile == debug_stile) ? 3 : 0;
						int dl =  (nsTile == debug_stile) ? 3 : 0;
						// plane_disp_strength
						if (dl > 2) {
							String [] dbg_titles = showSupertileSeparationTitles(plane_disp_strength[nsTile], null);
							double [][] dbg_img = showSupertileSeparation(false,plane_disp_strength[nsTile], null);
							showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
							sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "HOR_SEP_DISP"+nsTile,dbg_titles);
							dbg_img = showSupertileSeparation(true, plane_disp_strength[nsTile], null);
							sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "HOR_SEP_WORLD"+nsTile,dbg_titles);
						}
						
						ArrayList<TilePlanes.PlaneData> st_planes = new ArrayList<TilePlanes.PlaneData>();

						double[][][] disp_strength = new double[measuredLayers.getNumLayers()][][];
						
						for (int ml = 0; ml < disp_strength.length; ml++) if ((stMeasSel & ( 1 << ml)) != 0){
							
							if (smplMode) {
								disp_strength[ml] =  measuredLayers.getDisparityStrength(
										ml,             // int num_layer,
										stileX,         // int stX,
										stileY,         // int stY,
										null,           // boolean [] sel_in,
										strength_floor, // double strength_floor,
										strength_pow,   // double strength_pow,
										smplSide,       // int        smplSide, // = 2;   // Sample size (side of a square)
										smplNum,        //int        smplNum,   // = 3;   // Number after removing worst (should be >1)
										smplRms,        //double     smplRms,   // = 0.1; // Maximal RMS of the remaining tiles in a sample
										true);          // boolean null_if_none);
							} else {
								disp_strength[ml] =  measuredLayers.getDisparityStrength(
										ml,             // int num_layer,
										stileX,         // int stX,
										stileY,         // int stY,
										null,           // boolean [] sel_in,
										strength_floor, // double strength_floor,
										strength_pow,   // double strength_pow,
										true);          // boolean null_if_none);
							}
							
						}
						
						
						boolean OK;
						double [][] mm = maxMinMax[nsTile];
						if (mm == null){
//							double [][][] dbg_min_max = maxMinMax;
//							System.out.println("maxMinMax["+nsTile+"] == null");
							continue;
						}
						double [][] max_only = new double [(mm.length + 1)/2][2];
						for (int i = 0; i < max_only.length; i++){
							max_only[i] = mm[2 * i];
						}
						boolean [][][] plane_sels = null;
						
						int num_ml = disp_strength.length;
//						int num_p  = max_only.length;
						int num_tiles = 4 * superTileSize * superTileSize;
						int [] num_sel; 
						for (int iter = 0; iter < 2; iter ++){
							int num_p  = max_only.length;
							plane_sels = new boolean[num_p][num_ml][];
							num_sel = new int [num_p];
							for (int np = 0; np < num_p; np++) {
								for (int ml = 0; ml < num_ml; ml++) if (disp_strength[ml] != null) {
									plane_sels[np][ml] = new boolean[num_tiles];
								}
							}
							// compare closest to be able to use tilted planes later
							for (int ml = 0; ml < num_ml; ml++) if (disp_strength[ml] != null) {
								for (int indx = 0; indx < num_tiles; indx++) if (disp_strength[ml][1][indx] > 0.0){
									int best_plane = -1;
									double best_d2 = Double.NaN;
									for (int np = 0; np < num_p; np++) {
										double d2 = max_only[np][0] - disp_strength[ml][0][indx];
										// add disp_norm correction here?
										d2 *= d2;
										if (!(d2 >= best_d2)){
											best_d2 = d2;
											best_plane = np;
										}
									}
									if (best_plane >= 0){ // compare to max diff here too
										plane_sels[best_plane][ml][indx] = true; // so far exclusive
									}
								}
							}
							// recalculate average disparities for each plane and show number of tiles in each in debug mode
							for (int np = 0; np < num_p; np++) {
								double sd = 0.0, sw = 0.0;
								int nt = 0;
								for (int ml = 0; ml < num_ml; ml++) if (disp_strength[ml] != null) {
									for (int indx = 0; indx < num_tiles; indx++) if (plane_sels[np][ml][indx]){
										double w = disp_strength[ml][1][indx];
										sd += w * disp_strength[ml][0][indx];
										sw += w;
										num_sel[np]++;
									}
								}
								if (sw > 0) {
									sd /= sw;
								}
								if (dl > 0) {
									System.out.println("plane num_sel["+np+"] = "+num_sel[np]+" disp "+max_only[np][0]+"->"+sd+
											", weight "+max_only[np][1]+"->"+sw);
								}
								max_only[np][0] = sd;
								max_only[np][1] = sw;
							}
							// calculate transitions matrix (to find candidates for merge
							int [][]    trans_mat =  getTransMatrix(plane_sels);
							double [][] rel_trans =  getTransRel(trans_mat);
							
							
							if (dl > 0) {
								System.out.println("trans_mat = ");
								for (int i = 0; i < trans_mat.length; i++){
									System.out.print(i+": ");
									for (int j = 0; j < trans_mat[i].length; j++){
										System.out.print(trans_mat[i][j]+" ");
									}
									System.out.println();
								}
								System.out.println("rel_trans = ");
								for (int i = 0; i < rel_trans.length; i++){
									System.out.print(i+": ");
									for (int j = 0; j < rel_trans[i].length; j++){
										System.out.print(rel_trans[i][j]+" ");
									}
									System.out.println();
								}
							}
							if ((iter > 0 ) && (num_p > 1)){ // remove /join bad 
								int windx = 0;
								int remove_indx = -1;
								for (int i = 1; i < num_p; i++)	if (num_sel[i] < num_sel[windx]) windx = i;
								if (num_sel[windx] < plMinPoints) {
									if (debugLevel > 0){
										System.out.println ("processPlanes(): stileX = "+stileX+" stileY="+stileY+
												": removing plane "+windx+" with "+num_sel[windx]+" tiles ( <"+plMinPoints+")");
									}
									remove_indx = windx;
								}
								if (remove_indx < 0) {
									// find candidates for merge
									windx = -1; 
									for (int i = 0; i < (num_p - 1); i++)	{
										if (((max_only[i+1][0] - max_only[i][0]) < smallDiff) &&  // close enough to consider merging
												(rel_trans[i][i+1] > highMix)) {
											if ((windx < 0) || (rel_trans[i][i+1] > rel_trans[windx][windx+1])) windx = i;
										}
									}
									if (windx >=0 ) {
										if (debugLevel > 0){
											System.out.println ("processPlanes(): stileX = "+stileX+" stileY="+stileY+
													": merging plane "+windx+" with " + (windx + 1)+": "+
													num_sel[windx] + " and "+num_sel[windx+1]+" tiles, "+
													" rel_trans="+rel_trans[windx][windx + 1]+ " ( > " + highMix+"),"+
													" diff="+ (max_only[windx + 1][0]- max_only[windx][0]) + " ( < " + smallDiff+" ),"+
													" disp1 = "+max_only[windx][0]+" disp2 = "+max_only[windx + 1][0]);
										}
										double sum_w = max_only[windx][1] + max_only[windx + 1][1];
										max_only[windx+1][0] = (max_only[windx][0]*max_only[windx][1] +  max_only[windx+1][0]*max_only[windx+1][1]) / sum_w;
										max_only[windx+1][1] = sum_w;
										remove_indx = windx;
									}
								}
								if (remove_indx >= 0){
									double [][] max_only_copy = max_only.clone();
									for (int i = 0; i < max_only.length; i++) max_only_copy[i] = max_only[i];
									max_only = new double [max_only.length - 1][];
									int indx = 0;
									for (int i = 0; i < max_only_copy.length; i++) if (i != remove_indx) max_only[indx++] =max_only_copy[i]; 
									iter = 0;
									continue; // restart from 0

								}
								
								// Show other candidates for merge
								if (debugLevel > 0){
									double max_sep = 0.2;
									if (iter  > 0) {
										for (int i = 0; i < (num_p-1); i++){
											if (rel_trans[i][i+1] > max_sep) {
												System.out.println("processPlanes4() stileX = "+stileX+" stileY="+stileY+" lowplane = "+i+
														" num_sel1 = "+num_sel[i] + " num_sel2 = "+num_sel[i+1] +
														" rel_trans="+rel_trans[i][i+1]+
														" diff="+ (max_only[i+1][0]- max_only[i][0]) +
														" disp1 = "+max_only[i][0]+" disp2 = "+max_only[i+1][0]);
											}
										}
									}
								}								
							}
						}
						
						
						if (dl > 2) {
							String [] dbg_titles = showSupertileSeparationTitles( disp_strength, plane_sels);
							double [][] dbg_img = showSupertileSeparation(false,disp_strength, plane_sels);
							showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
							sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "initial_separation_disp"+nsTile,dbg_titles);
							dbg_img = showSupertileSeparation(true, disp_strength, plane_sels);
							sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "initial_separation_world"+nsTile,dbg_titles);
						}

						if (dl > 2) {
							double [] world_hor = {0.0, 1.0, 0.0};
							double sd = 0.0, sw = 0.0;
							for (int i = 0; i < max_only.length; i++){
								sd += max_only[i][0] * max_only[i][1];
								sw += max_only[i][1];
							}
							if (sw > 0) {
								System.out.println("Horizontally tilted disparity for stileX = "+stileX+" stileY="+stileY+", average disparity "+(sd/sw));
								double [][][] hor_disp_strength = pd0.getDisparityToPlane(
										world_hor, // double []     world_normal_xyz,
										sd / sw, // average disparity // double        disp_center,
										null, // boolean [][]  tile_sel, // null - do not use, {} use all (will be modified)
										disp_strength, // double [][][] disp_str, // calculate just once if null 
										1); // int           debugLevel);
								String [] dbg_titles = showSupertileSeparationTitles( hor_disp_strength, plane_sels);
								double [][] dbg_img = showSupertileSeparation(false, hor_disp_strength, plane_sels);
								showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
								sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "hor_separation_disp"+nsTile,dbg_titles);
								dbg_img = showSupertileSeparation(true, hor_disp_strength, plane_sels);
								sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "hor_separation_world"+nsTile,dbg_titles);
							}
						}
						
						for (int m = 0; m < max_only.length; m++) {
							//									TilePlanes.PlaneData pd = pd0_full.clone(); 
							TilePlanes.PlaneData pd = pd0.clone(); 
							OK = (pd.getPlaneFromMeas(
									plane_sels[m], // tile_sel,       // boolean [][] tile_sel, // null - do not use, {} use all (will be modified)
									disp_strength,
									Double.NaN,    // double       disp_far, // minimal disparity to select (or NaN)
									Double.NaN,    // double       disp_near, // maximal disparity to select (or NaN)
									0.0,            // plDispNorm,  // double       dispNorm,   //  Normalize disparities to the average if above
									0.0,            // double       min_weight,
									plMinPoints,    // int          min_tiles,
									strength_floor, // 
									strength_pow,   // double       strength_pow,
									// update !
									smplMode,
									smplSide,
									smplNum,
									smplRms,
									dl) != null);            // int          debugLevel)
							if (OK) {
								if (dl > 0) {
									if (pd.getWeight() > 1.0) {
										System.out.println("Processing subplane["+nsTile+"]["+m+"]"+
												", stileX="+stileX+
												", stileY="+stileY+
												", numPoints="+ pd.getNumPoints()+
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
									OK = pd.removeOutliers( // getPlaneFromMeas should already have run
											disp_strength, 
											targetV,      // double     targetEigen, // target eigenvalue for primary axis (is disparity-dependent, so is non-constant)
											max_outliers, // int        maxRemoved,  // maximal number of tiles to remove (not a constant)
											dl); // int        debugLevel)
									if (!OK) {
										continue;
									}
									if (dl > 0) {
										if (pd.getWeight() > 1.0) {
											System.out.println("Removed outliers["+nsTile+"]["+m+"]"+
													", stileX="+stileX+
													", stileY="+stileY+
													", numPoints="+ pd.getNumPoints()+
													", swc = "+pd.getWeight()+
													", center=["+pd.getZxy()[0]+","+pd.getZxy()[1]+","+pd.getZxy()[2]+"]"+
													", eig_val = {"+pd.getValues()[0]+","+pd.getValues()[1]+","+pd.getValues()[2]+"}"+
													", eig_vect[0] = {"+pd.getVector()[0]+","+pd.getVector()[1]+","+pd.getVector()[2]+"}");
										}
									}
								}
								double [] norm_xyz = pd.getWorldXYZ(
										correct_distortions);
								st_planes.add(pd);
								if (dl > 0) {
									System.out.println("World normal["+nsTile+"]["+m+"] = {"+
											norm_xyz[0]+", "+norm_xyz[1]+", "+norm_xyz[2]+"}");

								}
							}

						}
						if (st_planes.size() > 0){
							st_planes.add(0, st_planes.get(0)); // insert dummy at pos 0;
							planes[nsTile] = st_planes.toArray(new TilePlanes.PlaneData[0] );
							planes[nsTile][0] = null; // remove dummy
							if (dl >0){
								System.out.println("processPlanes4(): nsTile="+nsTile);
							}
						}
						//						}
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

	// TODO: Obsolete (never used)?
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
						boolean correct_distortions = planes[nsTile][np].getCorrectDistortions();
						int [] sTileXY0 = {stx0,sty0}; 
						if (tpl == null) {
							tpl = new TilePlanes(tileSize,superTileSize, geometryCorrection);
						}
						TilePlanes.PlaneData pd = tpl.new PlaneData(
								sTileXY0,
								tileSize,
								superTileSize,
								geometryCorrection,
								correct_distortions);
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
												true, // boolean   sum_weights,
												preferDisparity,
												(nsTile0 == debug_stile)? 1:0); // int       debugLevel)
										System.out.println("other with other, same weight:");
										TilePlanes.PlaneData other_other_pd = other_pd.mergePlaneToThis(
												other_pd, // PlaneData otherPd,
												1.0,      // double    scale_other,
												true,     // boolean   ignore_weights,
												true, // boolean   sum_weights,
												preferDisparity,
												(nsTile0 == debug_stile)? 1:0); // int       debugLevel)
										System.out.println("other with this, same weight:");
									}
									TilePlanes.PlaneData merged_pd = planes[nsTile0][indx].mergePlaneToThis(
											other_pd, // PlaneData otherPd,
											1.0,      // double    scale_other,
											true,     // boolean   ignore_weights,
											true, // boolean   sum_weights,
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
									true, // boolean   sum_weights,
									preferDisparity,
									(nsTile0 == debug_stile)? 1:0); // int       debugLevel)
							if (merged_pd != null) {
								merged_pd.scaleWeight(0.5);
							}
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
		final int debug_stile = dbg_Y * stilesX + dbg_X;

		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		// Select best symmetrical match, consider only N, NE, E, SE - later opposite ones will be copied
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					TilePlanes.PlaneData [] dbg_planes = null; 
					for (int nsTile0 = ai.getAndIncrement(); nsTile0 < nStiles; nsTile0 = ai.getAndIncrement()) {
						int sty0 = nsTile0 / stilesX;  
						int stx0 = nsTile0 % stilesX;
						int dl = ((debugLevel > -1) && (nsTile0 == debug_stile)) ? 1:0;
						if ( planes[nsTile0] != null) {
							if (dl > 0){
								System.out.println("matchPlanes(): nsTile0 ="+nsTile0);
								dbg_planes = planes[nsTile0];
							}
							for (int np0 = 0; np0 < planes[nsTile0].length; np0++){ // nu
								//										planes[nsTile0][np0].initNeibBest(); // 
								TilePlanes.PlaneData this_plane = planes[nsTile0][np0];
								if (this_plane != null) {
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
														if (other_planes[np] != null) {
															TilePlanes.PlaneData other_plane = this_plane.getPlaneToThis(
																	other_planes[np],
																	dl-1); // debugLevel);
															if (other_plane !=null) { // now always, but may add later
																TilePlanes.PlaneData merged_pd = this_plane.mergePlaneToThis(
																		other_plane, // PlaneData otherPd,
																		1.0,         // double    scale_other,
																		false,       // boolean   ignore_weights,
																		true, // boolean   sum_weights,
																		preferDisparity, 
																		dl-1); // int       debugLevel)

																if (merged_pd !=null) { // now always, but may add later
																	///															merged_pd.scaleWeight(0.5);
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
							if (dl > 0){
								System.out.println("matchPlanes(): nsTile0 ="+nsTile0+ " Done.");
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
						int dl = ((debugLevel > -1) && (nsTile0 == debug_stile)) ? 1:0;
						if (dl>0) {
							System.out.println("matchPlanes() nsTile0="+nsTile0);
						}
						if ( planes[nsTile0] != null) {
							for (int np0 = 0; np0 < planes[nsTile0].length; np0++){ // nu
								TilePlanes.PlaneData this_plane = planes[nsTile0][np0];
								if (this_plane != null) {
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
						if (dl>0) {
							System.out.println("matchPlanes() nsTile0="+nsTile0);
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
											int np_min = LOWEST_PLANE(neib_worse.length);
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

	public void testResoveTriangle(
			double     rquality,
			double     weakWorsening,
			double     okMergeEigen,
			double     maxWorldSin2,
			double     dispNorm,
			double     maxEigen, // maximal eigenvalue of planes to consider
			boolean    preferDisparity,
			int [][][] conflicts,
			int        debugLevel,
			int        dbg_X,
			int        dbg_Y)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int nsTile = dbg_Y * stilesX + dbg_X;
		final TileSurface.TileNeibs tnSurface = tileSurface.new TileNeibs(stilesX, stilesY);
		if (conflicts[nsTile] == null){
			System.out.println("There are no conflicts for nsTile = "+nsTile);
			return;
		}
		System.out.println("There are "+ conflicts[nsTile].length+" conflicts for nsTile = "+nsTile);
		for (int nConfl = 0; nConfl < conflicts[nsTile].length; nConfl++){
			String sdir_mask = "";
			for (int i = 7; i >=0; i--) sdir_mask +=((conflicts[nsTile][nConfl][2] & (1 << i)) != 0) ? i:"."; 
			System.out.println("Resolving conflict "+ nConfl+":"+
					" nl1 = "+conflicts[nsTile][nConfl][0]+
					" nl2 = "+conflicts[nsTile][nConfl][1]+
					" dir_mask = "+ sdir_mask+" ("+conflicts[nsTile][nConfl][2]+")");
			int win_layer = resolveTriangularConflict(
					nsTile,
					conflicts[nsTile][nConfl][0], // int nl1,
					conflicts[nsTile][nConfl][1], // int nl2,
					conflicts[nsTile][nConfl][2], // int dir_mask,
					rquality,
					weakWorsening,
					okMergeEigen,
					maxWorldSin2,
					dispNorm,
					maxEigen, // maximal eigenvalue of planes to consider
					tnSurface,
					preferDisparity,
					debugLevel);
			if (win_layer <0) {
				System.out.println("Conflict "+ nConfl+":"+
						" nl1 = "+conflicts[nsTile][nConfl][0]+
						" nl2 = "+conflicts[nsTile][nConfl][1]+
						" dir_mask = "+ sdir_mask+" ("+conflicts[nsTile][nConfl][2]+") is NOT RESOLVED");
			} else {
				System.out.println("Conflict "+ nConfl+":"+
						" nl1 = "+conflicts[nsTile][nConfl][0]+
						" nl2 = "+conflicts[nsTile][nConfl][1]+
						" dir_mask = "+ sdir_mask+" ("+conflicts[nsTile][nConfl][2]+") is RESOLVED - use layer "+win_layer);
				
			}
		}
	}

	public int [] resolveDiagonalTriangularConflicts(
			int [][][] conflicts,
			Conflicts conflict_stats, // to be updated after applying resolution
			int        starSteps, // How far to look around when calculationg connection cost
			double     orthoWeight,
			double     diagonalWeight,
			double     starPwr, // Divide cost by number of connections to this power
			double     dblTriLoss, //  When resolving double triangles allow minor degradation (0.0 - strict)
			boolean    preferDisparity,
			int        debugLevel,
			int        dbg_X,
			int        dbg_Y)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int dbgTile = dbg_Y * stilesX + dbg_X;
		final TileSurface.TileNeibs tnSurface = tileSurface.new TileNeibs(stilesX, stilesY);
		int [] rslt = {0,0};
		for (int nsTile = 0; nsTile < conflicts.length; nsTile++) if (conflicts[nsTile] != null) {
			// conflicts may disappear after being fixed, recheck for null
			for (int nConfl = 0; (conflicts[nsTile] != null) && (nConfl < conflicts[nsTile].length); nConfl++){
				int dl = ((debugLevel > -1) && (nsTile == dbgTile)) ? 3 : 0;
				Conflict confl = new Conflict(nsTile, conflicts[nsTile][nConfl]);
				for (int dir_orient = 0; dir_orient < 16; dir_orient++){
					boolean right = (dir_orient & 8) != 0;
					boolean reverse = (dir_orient & 4) != 0;
					int start_dir4 = dir_orient & 3;
					if (confl.conflictExists(
							start_dir4,
							right,
							reverse)){
						boolean OK = resolveDiagonalTriangularConflict(
								confl,           // Conflict conflict,
								start_dir4,      // int start_dir4,
								right,           // boolean right,
								reverse,         // boolean reverse,
								tnSurface,       // TileSurface.TileNeibs tnSurface,
								conflicts,       // int [][][] conflicts,
								conflict_stats,  // Conflicts conflict_stats, // to be updated after applying resolution
								starSteps, // How far to look around when calculationg connection cost
								orthoWeight,     // double     orthoWeight,
								diagonalWeight,  // double     diagonalWeight,
								starPwr,         // double     starPwr, // Divide cost by number of connections to this power
								dblTriLoss,      //  double     dblTriLoss, //  When resolving double triangles allow minor degradation (0.0 - strict)
								preferDisparity, // boolean    preferDisparity,
								dl);             // int        debugLevel)
						if (OK) rslt[0]++; // will not count if some conflicts disappeared as a result
						else    rslt[1]++;
					}
				}
			}
		}
		return rslt;
	}
	
	
	
	public boolean resolveDiagonalTriangularConflict(

			Conflict conflict,
			int start_dir4,
			boolean right,
			boolean reverse,
			TileSurface.TileNeibs tnSurface,
			int [][][] conflicts,
			Conflicts conflict_stats, // to be updated after applying resolution
			int        starSteps, // How far to look around when calculationg connection cost
			double     orthoWeight,
			double     diagonalWeight,
			double     starPwr, // Divide cost by number of connections to this power
			double     dblTriLoss, //  When resolving double triangles allow minor degradation (0.0 - strict)
			boolean    preferDisparity,
			int        debugLevel)
	{
		if (debugLevel > 1) {
			System.out.println("resolveDiagonalTriangularConflict(): "+conflict.toString()+": start_dir4 = "+ start_dir4 + ", right = "+right+", reverse = "+reverse);
		}
		
		if (!conflict.conflictExists(
				start_dir4,
				right,
				reverse)) {
			System.out.println("conflict does not exist");
			return false; // should be tested before to distinguish between non-existent and failure to resolve
		}
		Conflicts iconflicts = new Conflicts(this); 
		int dir01 = 2 * start_dir4;
		int dir12 = 2 * ((start_dir4 + (right? 1 : 3)) % 4);
		int dir20 = (dir12 + (right? 3 : 5)) % 8;
		int dir02 = (dir20 + 4) % 8;


		int nsTile = conflict.getSTile();
		int nl = conflict.getStartLayer(reverse);
		int nl_end = conflict.getEndLayer(reverse);
		int nsTile1 = tnSurface.getNeibIndex(nsTile, dir01);
		int nl1 = planes[nsTile][nl].getNeibBest(dir01);
		int nsTile2 = tnSurface.getNeibIndex(nsTile1, dir12);
		if (nl1 <0){
			System.out.println(conflict.toString()+" : start_dir4 = "+ start_dir4 + ", right = "+right+", reverse = "+reverse+" nl1 < 0 !");
			return false;
		}
		
		if (debugLevel > 1) {
			if (debugLevel > 1) {
				System.out.println("resolveDiagonalTriangularConflict(): nsTile = "+nsTile+" nsTile1 = "+nsTile1+" nsTile2="+nsTile2);
			}
		}
		
		int nl2 = planes[nsTile1][nl1].getNeibBest(dir12);
		// now nsTile2:nl2 is supposed to be connected to nsTile:nl_end != nl

		int [] mod_supertiles = {nsTile,nsTile2};
		int [] nsTiles = getInvolvedSupertiles( // first mod_supertiles.length entries will be mod_supertiles[]
				mod_supertiles,
				tnSurface);
		if (debugLevel > 1) {
			System.out.println("Involved supertiles:");
			for (int i = 0; i < nsTiles.length; i++){
				System.out.println(i+":"+nsTiles[i]);
			}
		}
		HashMap<Integer,Integer> replacement_tiles = new HashMap<Integer,Integer>();
		for (int i = 0; i < mod_supertiles.length; i++){
			replacement_tiles.put(mod_supertiles[i], new Integer(i));
		}
		ConnectionCosts connectionCosts = new ConnectionCosts(
				orthoWeight,
				diagonalWeight,
				starPwr,    // Divide cost by number of connections to this power
				starSteps,
				this.planes,
				tnSurface,
				preferDisparity);
		
		int [][][] neibs_prev = connectionCosts.initConnectionCosts(mod_supertiles);
		
		int [][][] neibs_prev_old = new int [mod_supertiles.length][][];
		double [][][]  val_weights = new double [mod_supertiles.length][][];
		
		// Calculate original costs and neighhbors
		updateConnectionsCost (
				mod_supertiles,      // int []         nsTiles,
				null,         // int [][][]     neibs_prev,
				neibs_prev_old,   // int [][][]     neibs,		// should be initialized at top dimension if neibs_prev==null	
				val_weights,  // double [][][]  val_weights, // should be initialized at top dimension if neibs_prev==null
				orthoWeight,
				diagonalWeight,
				starPwr,        // double     starPwr, // Divide cost by number of connections to this power
				tnSurface,
				preferDisparity,
				debugLevel);

		if ((neibs_prev == null) || (neibs_prev[0] == null) || (neibs_prev[1] == null)){
			System.out.println("resolveDiagonalTriangularConflict(): BUG");
			System.out.println("neibs_prev="+neibs_prev);
			System.out.println("neibs_prev[0]="+neibs_prev[0]);
			System.out.println("neibs_prev[1]="+neibs_prev[1]);
			return false;
		}
		int [][][] neibs = {neibs_prev[0].clone(),neibs_prev[1].clone()}; // both non-null
		for (int nt = 0; nt < neibs.length; nt++){
			for (int l = 0; l < neibs_prev[nt].length; l++) if (neibs_prev[nt][l]!= null){
				neibs[nt][l] = neibs_prev[nt][l].clone();
			}
		}
		int [][][] neibs_old = {neibs_prev[0].clone(),neibs_prev[1].clone()}; // both non-null
		for (int nt = 0; nt < neibs.length; nt++){
			for (int l = 0; l < neibs_prev[nt].length; l++) if (neibs_prev[nt][l]!= null){
				neibs_old[nt][l] = neibs_prev_old[nt][l].clone();
			}
		}
		int [][][] conflicts_old = new int [nsTiles.length][][];
		for (int isTile = 0; isTile < nsTiles.length; isTile++){
			conflicts_old[isTile] = iconflicts.detectTriangularTileConflicts(
					nsTiles[isTile],   // int nsTile0,
					replacement_tiles, //HashMap<Integer,Integer> replacement_tiles, //
					neibs_prev, // int [][][] replacement_neibs,
					tnSurface); // TileSurface.TileNeibs tnSurface)
		}
		if (debugLevel > 1) {
			System.out.println("Calculating original conflicts");
		}
		iconflicts.addConflicts(conflicts_old,
				debugLevel - 1); // debugLevel);

		if (nl2 <0){
			System.out.println(conflict.toString()+" : start_dir4 = "+ start_dir4 + ", right = "+right+", reverse = "+reverse+" nl2 < 0 !");
			return false;
		}
		
		// now swap diagonal connection
		neibs[0][nl][dir02] = nl2;
		neibs[1][nl2][dir20] = nl;
		neibs[0][nl_end][dir02] = neibs_prev[0][nl][dir02];
		if (neibs_prev[0][nl][dir02] >= 0) {
			neibs[1][neibs_prev[0][nl][dir02]][dir20] = nl_end;
		}

		neibs_old[0][nl][dir02] = nl2;
		neibs_old[1][nl2][dir20] = nl;
		neibs_old[0][nl_end][dir02] = neibs_prev_old[0][nl][dir02];
		if (neibs_prev_old[0][nl][dir02] >= 0) {
			neibs_old[1][neibs_prev_old[0][nl][dir02]][dir20] = nl_end;
		}
		
		
		
		int [][][] new_conflicts = new int [nsTiles.length][][];
		
		double new_costs_diff = 	connectionCosts.getConnectionsCostDiff(
				neibs,
				debugLevel);
		
		if (debugLevel > -1) {
			System.out.println("resolveDiagonalTriangularConflict(): resolving conflict "+conflict.toString()+
					", start dir4 = "+start_dir4+
					", orientation = "+(right?"right":"left") +
					", improvement (negative diff) = "+new_costs_diff);
			
		}
		double new_costs_diff_old = 	updateConnectionsCost (
				mod_supertiles,      // int []         nsTiles,
				neibs_prev_old,          // int [][][]     neibs_prev,
				neibs_old,               // int [][][]     neibs,		// should be initialized at top dimension if neibs_prev==null	
				val_weights,         // double [][][]  val_weights, // should be initialized at top dimension if neibs_prev==null
				orthoWeight,
				diagonalWeight,
				starPwr,        // double     starPwr, // Divide cost by number of connections to this power
				tnSurface,
				preferDisparity,
				debugLevel);
				
		for (int isTile = 0; isTile < nsTiles.length; isTile++){
			new_conflicts[isTile] = iconflicts.detectTriangularTileConflicts(
					nsTiles[isTile],   // int nsTile0,
					replacement_tiles, //HashMap<Integer,Integer> replacement_tiles, //
					neibs, // neibs, // int [][][] replacement_neibs,
					tnSurface); // TileSurface.TileNeibs tnSurface)
		}
		Conflicts new_conflicts_stats = new Conflicts(
				new_conflicts,
				this,
				debugLevel - 1); // debugLevel);
		new_conflicts_stats.subConflicts(iconflicts); // subtract old number of different types of conflicts
	
		if (debugLevel > -1) {
			System.out.println("resolveDiagonalTriangularConflict() OLD: resolving conflict "+conflict.toString()+
					", start dir4 = "+start_dir4+
					", orientation = "+(right?"right":"left") +
					", improvement (negative diff) = "+new_costs_diff_old);
			
		}
		if (debugLevel > 0) { // -1) {
			new_conflicts_stats.printConflictSummary(
					"Conflicts difference after resolution:", true, true, true);
		}
		// See if it is good
		if (debugLevel > 0) { // -1) {
			if (new_conflicts_stats.getNumConflicts() > 0) {
				if (debugLevel > -1) System.out.println("FAILURE: number of conflicts increased");
			} else if (new_costs_diff > dblTriLoss){
				if (debugLevel > -1) System.out.println("FAILURE: cost increased too much");
			} else if ((new_costs_diff >= 0.0) && (new_conflicts_stats.getNumConflicts() == 0)){
				if (debugLevel > -1) System.out.println("FAILURE: cost increased, but conflicts number is the same");
			} else if (new_costs_diff >= 0.0){
				if (debugLevel > -1) System.out.println("SUCCESS: cost did not decrease, but conflicts did");
			} else if (new_conflicts_stats.getNumConflicts() == 0){
				if (debugLevel > -1) System.out.println("SUCCESS: conflicts number is the same, but connection cost decreased");
			} else {
				if (debugLevel > -1) System.out.println("SUCCESS: both conflicts and connection cost decreased");
			}
			if (debugLevel > 1) {
				System.out.println("resolveDiagonalTriangularConflict(): nsTile = "+conflict.getSTile());
			}
		}

		if (new_conflicts_stats.getNumConflicts() > 0) {
			return false;
		}
		if (new_costs_diff > dblTriLoss){
			return false;
		}
		if ((new_costs_diff >= 0.0) && (new_conflicts_stats.getNumConflicts() == 0)){
			return false;
		}
		boolean apply = true;
		// apply changes
		if (apply) {
			// update statistics
			conflict_stats.addConflicts(new_conflicts_stats);
			// update conflict
			for (int i = 0; i < nsTiles.length; i++){
				conflicts[nsTiles[i]]=  new_conflicts[i];
			}

			// apply resolution
			for (int i = 0; i < mod_supertiles.length; i++){
				for (int l = 0; l < neibs[i].length; l ++) if (neibs[i][l] != null){
					if (planes[mod_supertiles[i]].length <= l) {
						System.out.println("BUG: (planes["+mod_supertiles[i]+"].length="+planes[mod_supertiles[i]].length+" <= "+l);
						return false;
					}
					planes[mod_supertiles[i]][l].setNeibBest(neibs[i][l]);
				}
			}
		}
		return true;		
	}
	

	public boolean resolveMultiTriangularConflict(
			int nsTile,
			int nl1,
			int nl2,
			int dir_mask,
			TileSurface.TileNeibs tnSurface,
			int [][][] conflicts,
			Conflicts conflict_stats, // to be updated after applying resolution
			int        starSteps, // How far to look around when calculationg connection cost
			double     orthoWeight,
			double     diagonalWeight,
			double     starPwr, // Divide cost by number of connections to this power
			double     dblTriLoss, //  When resolving double triangles allow minor degradation (0.0 - strict)
			boolean    preferDisparity,
			int        debugLevel)
	{
		Conflicts iconflicts = new Conflicts(this); 
		// neibs_raw may have unused directions
		int [][][][] neibs_raw = getTriangularResolutionVariants(
				nsTile,
				nl1,
				nl2,
				dir_mask,
				tnSurface,
				debugLevel);
		if (neibs_raw.length == 0){
			System.out.println("resolveMultiTriangularConflict() - zero variants. A bug?");
			return false;
		}
		
		int num_tiles = 0;
		for (int i = 0; i < neibs_raw[0].length; i++){
			if (neibs_raw[0][i] != null) num_tiles++;
 		}
		int [] mod_supertiles = new int [num_tiles];
		int [] indices = new int[num_tiles];
		int indx = 0;
		indices[indx++] = 4;
		for (int i = 0; (i < neibs_raw[0].length) && (indx < indices.length); i++){
			if (neibs_raw[0][i] != null) indices[indx++] = i;
 		}
		for (int i = 0; i < indices.length; i++){
			int dir4 = indices[i];
			if (dir4 > 3) {
				mod_supertiles[i] = nsTile;
			} else {
				mod_supertiles[i] = tnSurface.getNeibIndex(nsTile, 2 * dir4);
			}
		}
		int [][][][] neibs_vars = new int [neibs_raw.length][mod_supertiles.length][][];
		for (int ivar = 0; ivar < neibs_raw.length; ivar++){
			for (int i = 0; i < indices.length; i++){
				neibs_vars[ivar][i] = neibs_raw[ivar][indices[i]];
			}
		}

		
		
		// See how this application will influence number of conflicts
		// All supertiles that may have different conflicts 
		int [] nsTiles = getInvolvedSupertiles( // first mod_supertiles.length entries will be mod_supertiles[]
				mod_supertiles,
				tnSurface);
		HashMap<Integer,Integer> replacement_tiles = new HashMap<Integer,Integer>();
		for (int i = 0; i < mod_supertiles.length; i++){
			replacement_tiles.put(mod_supertiles[i], new Integer(i));
		}
		
		ConnectionCosts connectionCosts = new ConnectionCosts(
				orthoWeight,
				diagonalWeight,
				starPwr,    // Divide cost by number of connections to this power
				starSteps,
				this.planes,
				tnSurface,
				preferDisparity);
		
		int [][][] neibs_prev = connectionCosts.initConnectionCosts(mod_supertiles);
		
/** */		
		int [][][] neibs_prev_old = new int [mod_supertiles.length][][];
		double [][][]  val_weights_original = new double [mod_supertiles.length][][];
		updateConnectionsCost (
				mod_supertiles,      // int []         nsTiles,
				null,         // int [][][]     neibs_prev,
				neibs_prev_old,   // int [][][]     neibs,		// should be initialized at top dimension if neibs_prev==null	
				val_weights_original,  // double [][][]  val_weights, // should be initialized at top dimension if neibs_prev==null
				orthoWeight,
				diagonalWeight,
				starPwr,        // double     starPwr, // Divide cost by number of connections to this power
				tnSurface,
				preferDisparity,
				debugLevel);
/** */
		int [][][] conflicts_old = new int [nsTiles.length][][];
		for (int isTile = 0; isTile < nsTiles.length; isTile++){
			conflicts_old[isTile] = iconflicts.detectTriangularTileConflicts(
					nsTiles[isTile],   // int nsTile0,
					replacement_tiles, //HashMap<Integer,Integer> replacement_tiles, //
					neibs_prev, // int [][][] replacement_neibs,
					tnSurface); // TileSurface.TileNeibs tnSurface)
		}
		if (debugLevel > 1) {
			System.out.println("Involved supertiles:");
			for (int i = 0; i < nsTiles.length; i++){
				System.out.println(i+":"+nsTiles[i]);
			}
		}

		if (debugLevel > 1) {
			System.out.println("Calculating original conflicts");
		}

		iconflicts.addConflicts(conflicts_old,
				debugLevel - 1); // debugLevel);
		// After getting old (referfence data) iterate through all variants, for each calculate cost and number of conflicts.
		// First - just collect data - cost and full statistics of the conflicts, print them
		// Then select the best one (not incrementing number of conflicts? and reducing cost)
		int [][][][] variant_conflicts = new int [neibs_vars.length][nsTiles.length][][];
		
		Conflicts [] variant_conflicts_stats = new Conflicts [neibs_vars.length];

		
		double [] variant_costs_diff = new double [neibs_vars.length];
		for (int variant = 0; variant < neibs_vars.length; variant ++){
			if (debugLevel > 0) {
				System.out.println("resolveMultiTriangularConflict(): resolving conflict for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +
						", dir_mask = "+dir_mask+" variant = "+variant);
				
			}
			variant_costs_diff[variant] = 	connectionCosts.getConnectionsCostDiff(
					neibs_vars[variant],
					debugLevel);
			
			if (debugLevel > -1) {
				System.out.println("resolveMultiTriangularConflict(): resolving conflict for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +
						", dir_mask = "+dir_mask+" variant = "+variant+" improvement (negative diff) = "+variant_costs_diff[variant]);
				
			}
			
			double [] variant_costs_diff_old = new double [neibs_vars.length];
			double [][][]  val_weights = val_weights_original.clone();
			for (int nt = 0; nt < val_weights_original.length; nt++) if (val_weights_original[nt] != null){
				val_weights[nt] = val_weights_original[nt].clone();
				for (int nl = 0; nl < val_weights_original[nt].length; nl++) if (val_weights_original[nt][nl] != null){
					val_weights[nt][nl] = val_weights_original[nt][nl].clone();
				}
			}
			variant_costs_diff_old[variant] = 	updateConnectionsCost (
					mod_supertiles,      // int []         nsTiles,
					neibs_prev_old,          // int [][][]     neibs_prev,
					neibs_vars[variant], // int [][][]     neibs,		// should be initialized at top dimension if neibs_prev==null	
					val_weights,         // double [][][]  val_weights, // should be initialized at top dimension if neibs_prev==null
					orthoWeight,
					diagonalWeight,
					starPwr,        // double     starPwr, // Divide cost by number of connections to this power
					tnSurface,
					preferDisparity,
					debugLevel);
		
			
			
			for (int isTile = 0; isTile < nsTiles.length; isTile++){
				variant_conflicts[variant][isTile] = iconflicts.detectTriangularTileConflicts(
						nsTiles[isTile],   // int nsTile0,
						replacement_tiles, //HashMap<Integer,Integer> replacement_tiles, //
						neibs_vars[variant], // neibs, // int [][][] replacement_neibs,
						tnSurface); // TileSurface.TileNeibs tnSurface)
			}
			variant_conflicts_stats[variant] = new Conflicts(
					variant_conflicts[variant],
					this,
					debugLevel - 1); // debugLevel);
			variant_conflicts_stats[variant].subConflicts(iconflicts); // subtract old number of different types of conflicts
			
			if (debugLevel > -1) {
				System.out.println("resolveMultiTriangularConflict() OLD: resolving conflict for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +
						", dir_mask = "+dir_mask+" variant = "+variant+" improvement (negative diff) = "+variant_costs_diff_old[variant]);
				
			}
			if (debugLevel > -1) {
				variant_conflicts_stats[variant].printConflictSummary(
						"Conflicts difference after resolution:", true, true, false);
			}
		}
		// How to compare? 1 attempt: none of the conflicts get worse, some get better or cost goes down
		int best_variant = -1;
		int [][] num_better_worse = new int [neibs_vars.length][2];
		for (int variant = 0; variant < neibs_vars.length; variant ++){
			int num_worse = variant_conflicts_stats[variant].numBetterWorse(
					false, // boolean better,
					false, // boolean use_all,
					true, // boolean use_odo,
					false); // ); // boolean use_ood)

			int num_better = variant_conflicts_stats[variant].numBetterWorse(
					true, // boolean better,
					false, // boolean use_all,
					true, // boolean use_odo,
					false); // ); // boolean use_ood)
			num_better_worse[variant][0] = num_better;
			num_better_worse[variant][1] = num_worse;
			if ((num_worse == 0) &&
					(variant_costs_diff[variant] <= dblTriLoss) && // not too worse
					((variant_costs_diff[variant] < 0) || (num_better > 0)) && // either 
					((best_variant < 0) ||
							( num_better_worse[variant][0] > num_better_worse[best_variant][0]) ||
							(( num_better_worse[variant][0] == num_better_worse[best_variant][0]) && 
									(variant_costs_diff[variant] < variant_costs_diff[best_variant])))){
				best_variant = variant;
			}
		}

		if (debugLevel > 1){
			System.out.println("resolveMultiTriangularConflict(): for tile "+nsTile);
		}
		
		if ((best_variant < 0) || (variant_costs_diff[best_variant] > dblTriLoss)){
			if (debugLevel > -1) {
				System.out.println("resolveMultiTriangularConflict(): FAILED find a sutable solution for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +
						", dir_mask = "+dir_mask+" of "+ neibs_vars.length+" variants");
				return false;
			}
		} else {
			if (debugLevel > -1) {
				System.out.println("resolveMultiTriangularConflict(): SUCCESS to find a sutable solution for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +
						", dir_mask = "+dir_mask+". Of "+ neibs_vars.length+" variants - use variant # " + best_variant+
						" cost difference (negative) = "+variant_costs_diff[best_variant]+" num conflict reductions = "+num_better_worse[best_variant][0]);
				variant_conflicts_stats[best_variant].printConflictSummary(
						"Conflicts number change per type: ",
						true, // use_all,
						true, //use_odo,
						true); // use_ood);
				iconflicts.printConflictSummary(
						"Conflicts before resolution: ",
						true, // use_all,
						true, //use_odo,
						true); // use_ood);
				// update statistics
				conflict_stats.addConflicts(variant_conflicts_stats[best_variant]);

				
				// update conflict
				for (int i = 0; i < nsTiles.length; i++){
					conflicts[nsTiles[i]]=  variant_conflicts[best_variant][i];
				}
				
				
				// apply resolution
				for (int i = 0; i < mod_supertiles.length; i++){
					for (int nl = 0; nl < neibs_vars[best_variant][i].length; nl ++) if (neibs_vars[best_variant][i][nl] != null){
						planes[mod_supertiles[i]][nl].setNeibBest(neibs_vars[best_variant][i][nl]);
					}
				}
			}
		}
		return true;
	}
	
	
	public int [] resolveMultiTriangularConflicts(
			int [][][] conflicts,
			Conflicts conflict_stats, // to be updated after applying resolution
			int        starSteps, // How far to look around when calculationg connection cost
			double     orthoWeight,
			double     diagonalWeight,
			double     starPwr, // Divide cost by number of connections to this power
			double     dblTriLoss, //  When resolving double triangles allow minor degradation (0.0 - strict)
			boolean    preferDisparity,
			int        debugLevel,
			int        dbg_X,
			int        dbg_Y)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int dbgTile = dbg_Y * stilesX + dbg_X;
		final TileSurface.TileNeibs tnSurface = tileSurface.new TileNeibs(stilesX, stilesY);
		int [] rslt = {0,0};
		for (int nsTile = 0; nsTile < conflicts.length; nsTile++) if (conflicts[nsTile] != null) {
			// conflicts may disappear after being fixed, recheck for null
			for (int nConfl = 0; (conflicts[nsTile] != null) && (nConfl < conflicts[nsTile].length); nConfl++){
				int dl = ((debugLevel > 0) && (nsTile == dbgTile)) ? 3 : 0;
				boolean OK = resolveMultiTriangularConflict(
						nsTile,
						conflicts[nsTile][nConfl][0], // int nl1,
						conflicts[nsTile][nConfl][1], // int nl2,
						conflicts[nsTile][nConfl][2], // int dir_mask,
						tnSurface,
						conflicts,
						conflict_stats, // to be updated after applying resolution
//						maxEigen, // maximal eigenvalue of planes to consider
						starSteps, // How far to look around when calculationg connection cost
						orthoWeight,
						diagonalWeight,
						starPwr,        // double     starPwr, // Divide cost by number of connections to this power
						dblTriLoss, //  When resolving double triangles allow minor degradation (0.0 - strict)
						preferDisparity,
						dl); // debugLevel,
				if (OK) rslt[0]++;
				else    rslt[1]++;
			}
		}
		return rslt;
	}

	
	public int [] resolveStarConflicts(
			int [][][] conflicts,
			Conflicts conflict_stats, // to be updated after applying resolution
			int        starSteps, // How far to look around when calculationg connection cost
			double     orthoWeight,
			double     diagonalWeight,
			double     starPwr, // Divide cost by number of connections to this power
			double     dblTriLoss, //  When resolving double triangles allow minor degradation (0.0 - strict)
			boolean    newConfl, // Allow more conflicts if overall cost is reduced
			boolean    preferDisparity,
			int        debugLevel,
			int        dbg_X,
			int        dbg_Y)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int dbgTile = dbg_Y * stilesX + dbg_X;
		final TileSurface.TileNeibs tnSurface = tileSurface.new TileNeibs(stilesX, stilesY);
		int [] rslt = {0,0};
		for (int nsTile = 0; nsTile < conflicts.length; nsTile++) if (conflicts[nsTile] != null) {
			// conflicts may disappear after being fixed, recheck for null
			for (int nConfl = 0; (conflicts[nsTile] != null) && (nConfl < conflicts[nsTile].length); nConfl++){
				int dl = ((debugLevel > 0) && (nsTile == dbgTile)) ? 3 : 0;
				boolean OK = resolveStarConflict(
						nsTile,
						conflicts[nsTile][nConfl][0], // int nl1,
						conflicts[nsTile][nConfl][1], // int nl2,
						tnSurface,
						conflicts,
						conflict_stats, // to be updated after applying resolution
						starSteps, // How far to look around when calculating connection cost
						orthoWeight,
						diagonalWeight,
						starPwr,        // double     starPwr, // Divide cost by number of connections to this power
						dblTriLoss, //  When resolving double triangles allow minor degradation (0.0 - strict)
						newConfl, // boolean    newConfl, // Allow more conflicts if overall cost is reduced
						preferDisparity,
						dl); // debugLevel,
				if (OK) rslt[0]++;
				else    rslt[1]++;
			}
		}
		return rslt;
	}

	
	public boolean resolveStarConflict0(
			int nsTile,
			int nl1,
			int nl2,
			TileSurface.TileNeibs tnSurface,
			int [][][] conflicts,
			Conflicts conflict_stats, // to be updated after applying resolution
			int        starSteps, // How far to look around when calculationg connection cost
			double     orthoWeight,
			double     diagonalWeight,
			double     starPwr, // Divide cost by number of connections to this power
			double     dblTriLoss, //  When resolving double triangles allow minor degradation (0.0 - strict)
			boolean    newConfl, // Allow more conflicts if overall cost is reduced
			boolean    preferDisparity,
			int        debugLevel)
	{
		
		if (newConfl && (dblTriLoss > 0.0)){
			dblTriLoss = 0.0; // require strict reducing of the cost if conflict increase is allowed
		}
		Conflicts iconflicts = new Conflicts(this); 

		TwoLayerNeighbors twoLayerNeighbors = new TwoLayerNeighbors();
		
		for (int dir = -1; dir < 8; dir++){
			int nt = tnSurface.getNeibIndex(nsTile, dir);
			if ((nt >= 0) && (planes[nt] != null)) {
				int [][] neibs =  new int [planes[nt].length][];
				for (int nl = 0; nl < planes[nt].length; nl++) if ( planes[nt][nl] != null){
					neibs[nl] =  planes[nt][nl].getNeibBest();
				}
				twoLayerNeighbors.setNeighbors(neibs,dir);
			}
		}
		twoLayerNeighbors.setLayers(nl1, nl2);
		if (debugLevel > 1) {
			System.out.println("resolveStarConflict(): nsTile ="+nsTile+" nl1="+nl1+" nl2="+nl2);
		}

		int [][][][] neibs_vars_dir = twoLayerNeighbors.getNeighborVariants(debugLevel);		

		int [] mod_supertiles = {nsTile};
		mod_supertiles = getInvolvedSupertiles( // first mod_supertiles.length entries will be mod_supertiles[]
				mod_supertiles,
				tnSurface);

		HashMap<Integer,Integer> replacement_tiles = new HashMap<Integer,Integer>();
		for (int i = 0; i < mod_supertiles.length; i++){
			replacement_tiles.put(mod_supertiles[i], new Integer(i));
		}
		
		// up to 9 tiles
		int [] indexToDir = new int[mod_supertiles.length];
		for (int i = 0; i < indexToDir.length; i++) indexToDir[i] = -1;
		for (int dir = -1; dir < 8; dir++){
			int nindx = (dir < 0) ? 8 : dir; 
			int nt = tnSurface.getNeibIndex(nsTile, dir);
			int indx = -1;
			Integer Indx = replacement_tiles.get(nt);
			if (Indx == null){
				System.out.println("resolveStarConflict(): nsTile = "+nsTile+" nindx="+nindx+" Indx == null ");
			} else {
				indx = Indx;
			}
			if (indx >= 0) {
				indexToDir[indx] = nindx;
			}
		}
		int [][][][] neibs_vars = new int [neibs_vars_dir.length][][][];
		for (int variant = 0; variant < neibs_vars_dir.length; variant ++){
			neibs_vars[variant] = new int [indexToDir.length][][];
			for (int i = 0; i <indexToDir.length; i++){
				if (indexToDir[i] >=0){
					neibs_vars[variant][i] = neibs_vars_dir[variant][indexToDir[i]];
				} else {
					System.out.println("resolveStarConflict(): a BUG:  indexToDir["+i+"] = "+indexToDir[i]);
				}
			}
		}
		
		// See how this application will influence number of conflicts
		// All supertiles that may have different conflicts 
		int [] nsTiles = getInvolvedSupertiles( // first mod_supertiles.length entries will be mod_supertiles[]
				mod_supertiles,
				tnSurface);
		
		ConnectionCosts connectionCosts = new ConnectionCosts(
				orthoWeight,
				diagonalWeight,
				starPwr,    // Divide cost by number of connections to this power
				starSteps,
				this.planes,
				tnSurface,
				preferDisparity);
		
		int [][][] neibs_prev = connectionCosts.initConnectionCosts(mod_supertiles);
		
		int [][][] conflicts_old = new int [nsTiles.length][][];
		for (int isTile = 0; isTile < nsTiles.length; isTile++){
			conflicts_old[isTile] = iconflicts.detectTriangularTileConflicts(
					nsTiles[isTile],   // int nsTile0,
					replacement_tiles, //HashMap<Integer,Integer> replacement_tiles, //
					neibs_prev, // int [][][] replacement_neibs,
					tnSurface); // TileSurface.TileNeibs tnSurface)
		}
		if (debugLevel > 1) {
			System.out.println("Involved supertiles:");
			for (int i = 0; i < nsTiles.length; i++){
				System.out.println(i+":"+nsTiles[i]);
			}
		}

		if (debugLevel > 1) {
			System.out.println("Calculating original conflicts");
		}

		iconflicts.addConflicts(conflicts_old,
				debugLevel - 1); // debugLevel);
		// After getting old (referfence data) iterate through all variants, for each calculate cost and number of conflicts.
		// First - just collect data - cost and full statistics of the conflicts, print them
		// Then select the best one (not incrementing number of conflicts? and reducing cost)
		int [][][][] variant_conflicts = new int [neibs_vars.length][nsTiles.length][][];
		
		Conflicts [] variant_conflicts_stats = new Conflicts [neibs_vars.length];

		
		double [] variant_costs_diff = new double [neibs_vars.length];
		for (int variant = 0; variant < neibs_vars.length; variant ++){
			for (int isTile = 0; isTile < nsTiles.length; isTile++){
				variant_conflicts[variant][isTile] = iconflicts.detectTriangularTileConflicts(
						nsTiles[isTile],   // int nsTile0,
						replacement_tiles, //HashMap<Integer,Integer> replacement_tiles, //
						neibs_vars[variant], // neibs, // int [][][] replacement_neibs,
						tnSurface); // TileSurface.TileNeibs tnSurface)
			}
			variant_conflicts_stats[variant] = new Conflicts(
					variant_conflicts[variant],
					this,
					debugLevel - 1); // debugLevel);
			variant_conflicts_stats[variant].subConflicts(iconflicts); // subtract old number of different types of conflicts
			if (debugLevel > 1) {
				System.out.println("resolveStarConflict(): resolving conflict for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +
						", variant = "+variant+
						", num_conflicts diff = "+variant_conflicts_stats[variant].getNumConflicts()
						);
				
			}
			if (!newConfl && (variant_conflicts_stats[variant].getNumConflicts() > 0)){
				continue;
			}
			variant_costs_diff[variant] = 	connectionCosts.getConnectionsCostDiff(
					neibs_vars[variant],
					debugLevel);
			
			if (debugLevel > 0) {
				System.out.println("resolveStarConflict(): resolving conflict for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +
						", variant = "+variant+" improvement (negative diff) = "+variant_costs_diff[variant]+
						", num_conflicts diff = "+variant_conflicts_stats[variant].getNumConflicts());
			}
			
			
			if (debugLevel > 0) { // -1) {
				variant_conflicts_stats[variant].printConflictSummary(
						"Conflicts difference after resolution:", true, true, false);
			}
		}
		// How to compare? 1 attempt: none of the conflicts get worse, some get better or cost goes down
		int best_variant = -1;
		int best_ignore_conflicts = -1;
		int [][] num_better_worse = new int [neibs_vars.length][2];
		/*
		for (int variant = 0; variant < neibs_vars.length; variant ++){
			int num_worse = variant_conflicts_stats[variant].numBetterWorse(
					false, // boolean better,
					false, // boolean use_all,
					true, // boolean use_odo,
					false); // ); // boolean use_ood)

			int num_better = variant_conflicts_stats[variant].numBetterWorse(
					true, // boolean better,
					false, // boolean use_all,
					true, // boolean use_odo,
					false); // ); // boolean use_ood)
			num_better_worse[variant][0] = num_better;
			num_better_worse[variant][1] = num_worse;
			if ((num_worse == 0) &&
					(variant_costs_diff[variant] <= dblTriLoss) && // not too worse
					((variant_costs_diff[variant] < 0) || (num_better > 0)) && // either 
					((best_variant < 0) ||
							( num_better_worse[variant][0] > num_better_worse[best_variant][0]) ||
							(( num_better_worse[variant][0] == num_better_worse[best_variant][0]) && 
									(variant_costs_diff[variant] < variant_costs_diff[best_variant])))){
				best_variant = variant;
			}
			if ((best_ignore_conflicts <0) || (variant_costs_diff[variant] < variant_costs_diff[best_ignore_conflicts])){
				best_ignore_conflicts = variant;
			}
		}
*/
		int [] num_var_conflicts = new int [neibs_vars.length];
		for (int variant = 0; variant < neibs_vars.length; variant ++){
//			int num_conflicts = variant_conflicts_stats[variant].getNumConflicts();
			num_var_conflicts[variant] = variant_conflicts_stats[variant].getNumConflicts();
			if ((num_var_conflicts[variant] <= 0) &&
					(variant_costs_diff[variant] <= dblTriLoss) && // not too worse
					((variant_costs_diff[variant] < 0) || (num_var_conflicts[variant] < 0)) && // either 
					((best_variant < 0) ||
							(num_var_conflicts[variant] < num_var_conflicts[best_variant]) ||
							((num_var_conflicts[variant] == num_var_conflicts[best_variant]) && 
									(variant_costs_diff[variant] < variant_costs_diff[best_variant])))){
				best_variant = variant;
			}
			if ((best_ignore_conflicts <0) || (variant_costs_diff[variant] < variant_costs_diff[best_ignore_conflicts])){
				best_ignore_conflicts = variant;
			}
		}

		if (debugLevel > 1){
			System.out.println("resolveStarConflict(): for tile "+nsTile);
		}
		
		if ((best_variant < 0) && newConfl && (variant_costs_diff[best_ignore_conflicts] < 0)){ // should be cost improvement
			best_variant = best_ignore_conflicts;
			if (debugLevel > -1) {
				System.out.println("resolveMultiTriangularConflict(): conflicts increase but cost decreases "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +" of "+ neibs_vars.length+" variants");
			}
		}
		if ((best_variant < 0) || (variant_costs_diff[best_variant] > dblTriLoss)){
			if (debugLevel > -1) {
				System.out.println("resolveMultiTriangularConflict(): FAILED find a sutable solution for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +" of "+ neibs_vars.length+" variants");
				return false;
			}
		} else {
			if (debugLevel > -1) {
				System.out.println("resolveStarConflict(): SUCCESS to find a sutable solution for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +". Of "+ neibs_vars.length+" variants - use variant # " + best_variant+
						" cost difference (negative) = "+variant_costs_diff[best_variant]+
						" num conflict reductions = "+(-num_var_conflicts[best_variant])+" (old "+num_better_worse[best_variant][0]+")");
				variant_conflicts_stats[best_variant].printConflictSummary(
						"Conflicts number change per type: ",
						true, // use_all,
						true, //use_odo,
						true); // use_ood);
				iconflicts.printConflictSummary(
						"Conflicts before resolution: ",
						true, // use_all,
						true, //use_odo,
						true); // use_ood);
				// update statistics
				conflict_stats.addConflicts(variant_conflicts_stats[best_variant]);

				
				// update conflict
				for (int i = 0; i < nsTiles.length; i++){
					conflicts[nsTiles[i]]=  variant_conflicts[best_variant][i];
				}
				
				
				// apply resolution
				for (int i = 0; i < mod_supertiles.length; i++){
					if (debugLevel > 1){
						System.out.println("resolveStarConflict(): nsTile = "+nsTile+ "mod_supertiles["+i+"]="+mod_supertiles[i]);
					}
					if (neibs_vars[best_variant][i] != null) {
						for (int nl = 0; nl < neibs_vars[best_variant][i].length; nl ++){
							if (debugLevel > 1){
								System.out.println("resolveStarConflict(): nl= = "+nl);
							}
							if (neibs_vars[best_variant][i][nl] != null){
								planes[mod_supertiles[i]][nl].setNeibBest(neibs_vars[best_variant][i][nl]);
							}
						}
					}
				}
			}
		}
		return true;
	}

	public boolean resolveStarConflict(
			int nsTile,
			int nl1,
			int nl2,
			TileSurface.TileNeibs tnSurface,
			int [][][] conflicts,
			Conflicts conflict_stats, // to be updated after applying resolution
			int        starSteps, // How far to look around when calculationg connection cost
			double     orthoWeight,
			double     diagonalWeight,
			double     starPwr, // Divide cost by number of connections to this power
			double     dblTriLoss, //  When resolving double triangles allow minor degradation (0.0 - strict)
			boolean    newConfl, // Allow more conflicts if overall cost is reduced
			boolean    preferDisparity,
			int        debugLevel)
	{
		
		if (newConfl && (dblTriLoss > 0.0)){
			dblTriLoss = 0.0; // require strict reducing of the cost if conflict increase is allowed
		}
		Conflicts iconflicts = new Conflicts(this); 

		TwoLayerNeighbors twoLayerNeighbors = new TwoLayerNeighbors();
		
		for (int dir = -1; dir < 8; dir++){
			int nt = tnSurface.getNeibIndex(nsTile, dir);
			if ((nt >= 0) && (planes[nt] != null)) {
				int [][] neibs =  new int [planes[nt].length][];
				for (int nl = 0; nl < planes[nt].length; nl++) if ( planes[nt][nl] != null){
					neibs[nl] =  planes[nt][nl].getNeibBest();
				}
				twoLayerNeighbors.setNeighbors(neibs,dir);
			}
		}
		twoLayerNeighbors.setLayers(nl1, nl2);
		if (debugLevel > 1) {
			System.out.println("resolveStarConflict(): nsTile ="+nsTile+" nl1="+nl1+" nl2="+nl2);
		}

		int [][][][] neibs_vars_dir = twoLayerNeighbors.getNeighborVariants(debugLevel);		

		int [] mod_supertiles = {nsTile};
		mod_supertiles = getInvolvedSupertiles( // first mod_supertiles.length entries will be mod_supertiles[]
				mod_supertiles,
				tnSurface);

		HashMap<Integer,Integer> replacement_tiles = new HashMap<Integer,Integer>();
		for (int i = 0; i < mod_supertiles.length; i++){
			replacement_tiles.put(mod_supertiles[i], new Integer(i));
		}
		
		// up to 9 tiles
		int [] indexToDir = new int[mod_supertiles.length];
		for (int i = 0; i < indexToDir.length; i++) indexToDir[i] = -1;
		for (int dir = -1; dir < 8; dir++){
			int nindx = (dir < 0) ? 8 : dir; 
			int nt = tnSurface.getNeibIndex(nsTile, dir);
			int indx = -1;
			Integer Indx = replacement_tiles.get(nt);
			if (Indx == null){
				System.out.println("resolveStarConflict(): nsTile = "+nsTile+" nindx="+nindx+" Indx == null ");
			} else {
				indx = Indx;
			}
			if (indx >= 0) {
				indexToDir[indx] = nindx;
			}
		}
		int [][][][] neibs_vars = new int [neibs_vars_dir.length][][][];
		for (int variant = 0; variant < neibs_vars_dir.length; variant ++){
			neibs_vars[variant] = new int [indexToDir.length][][];
			for (int i = 0; i <indexToDir.length; i++){
				if (indexToDir[i] >=0){
					neibs_vars[variant][i] = neibs_vars_dir[variant][indexToDir[i]];
				} else {
					System.out.println("resolveStarConflict(): a BUG:  indexToDir["+i+"] = "+indexToDir[i]);
				}
			}
		}
		
		// See how this application will influence number of conflicts
		// All supertiles that may have different conflicts 
		int [] nsTiles = getInvolvedSupertiles( // first mod_supertiles.length entries will be mod_supertiles[]
				mod_supertiles,
				tnSurface);

		ConnectionCosts connectionCosts = new ConnectionCosts(
				orthoWeight,
				diagonalWeight,
				starPwr,    // Divide cost by number of connections to this power
				starSteps,
				this.planes,
				tnSurface,
				preferDisparity);
		
		int [][][] neibs_prev = connectionCosts.initConnectionCosts(mod_supertiles);
// connectionCosts now contains last calculated	val/weight pairs for broader array of tile data	
		
		int [][][] conflicts_old = new int [nsTiles.length][][];
		double conflicts_old_cost = 0.0;
		for (int isTile = 0; isTile < nsTiles.length; isTile++){
			conflicts_old[isTile] = iconflicts.detectTriangularTileConflicts(
					nsTiles[isTile],   // int nsTile0,
					replacement_tiles, // HashMap<Integer,Integer> replacement_tiles, //
					neibs_prev, // int [][][] replacement_neibs,
					tnSurface); // TileSurface.TileNeibs tnSurface)
// calculate cost of conflicts			
			
			conflicts_old_cost += iconflicts.getConflictsCost(
					nsTiles[isTile],       // int           nsTile0,
					1.0,                   // double        scaleStartEnd, // include start and and layer tiles in the center in overall cost for each triangle (1.0)
					conflicts_old[isTile], // int [][]      conflicts, //
					replacement_tiles,     // HashMap<Integer,Integer> replacement_tiles, // null is OK
					neibs_prev,            // int [][][]    replacement_neibs,               // null OK if  replacement_tiles == null
					connectionCosts,       // ConnectionCosts connectionCosts,
					tnSurface); // TileSurface.TileNeibs tnSurface)

			
		}
		if (debugLevel > 1) {
			System.out.println("Sum of old conflict costs around " + nsTile + " is "+conflicts_old_cost);
			System.out.println("Involved supertiles:");
			for (int i = 0; i < nsTiles.length; i++){
				System.out.println(i+":"+nsTiles[i]);
			}
		}

		if (debugLevel > 1) {
			System.out.println("Calculating original conflicts");
		}

		iconflicts.addConflicts(conflicts_old,
				debugLevel - 1); // debugLevel);
		// After getting old (referfence data) iterate through all variants, for each calculate cost and number of conflicts.
		// First - just collect data - cost and full statistics of the conflicts, print them
		// Then select the best one (not incrementing number of conflicts? and reducing cost)
		int [][][][] variant_conflicts = new int [neibs_vars.length][nsTiles.length][][];
		
		Conflicts [] variant_conflicts_stats = new Conflicts [neibs_vars.length];

		
		double [] variant_costs_diff = new double [neibs_vars.length];
		double [] conflicts_var_cost = new double [neibs_vars.length];
		for (int variant = 0; variant < neibs_vars.length; variant ++){
			// connections costs are needed before conflicts 
			variant_costs_diff[variant] = 	connectionCosts.getConnectionsCostDiff(
					neibs_vars[variant],
					debugLevel);
			conflicts_var_cost[variant] = 0.0; 
			for (int isTile = 0; isTile < nsTiles.length; isTile++){
				variant_conflicts[variant][isTile] = iconflicts.detectTriangularTileConflicts(
						nsTiles[isTile],   // int nsTile0,
						replacement_tiles, //HashMap<Integer,Integer> replacement_tiles, //
						neibs_vars[variant], // neibs, // int [][][] replacement_neibs,
						tnSurface); // TileSurface.TileNeibs tnSurface)
				// connectionCosts now contains last calculated	val/weight pairs for broader array of tile data
				conflicts_var_cost[variant] += iconflicts.getConflictsCost(
						nsTiles[isTile],                    // int           nsTile0,
						1.0,                                // double        scaleStartEnd, // include start and and layer tiles in the center in overall cost for each triangle (1.0)
						variant_conflicts[variant][isTile], // int [][]      conflicts, //
						replacement_tiles,                  // HashMap<Integer,Integer> replacement_tiles, // null is OK
						neibs_vars[variant],                // int [][][]    replacement_neibs,               // null OK if  replacement_tiles == null
						connectionCosts,                    // ConnectionCosts connectionCosts,
						tnSurface);                         // TileSurface.TileNeibs tnSurface)
			}
			variant_conflicts_stats[variant] = new Conflicts(
					variant_conflicts[variant],
					this,
					debugLevel - 1); // debugLevel);
			variant_conflicts_stats[variant].subConflicts(iconflicts); // subtract old number of different types of conflicts
			if (debugLevel > 1) {
				System.out.println("resolveStarConflict(): resolving conflict for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +
						", variant = "+variant+
						"  conflit cost change = "+ (conflicts_var_cost[variant] - conflicts_old_cost)+
						", num_conflicts diff = "+variant_conflicts_stats[variant].getNumConflicts()
						);
				
// Initially do not process conflicts cost difference, just view it				
			}
/*			
			if (!newConfl && (variant_conflicts_stats[variant].getNumConflicts() > 0)){
				continue;
			}
*/			
			
			if (debugLevel > 0) {
				System.out.println("resolveStarConflict(): resolving conflict for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +
						", variant = "+variant+" improvement (negative diff) = "+variant_costs_diff[variant]+
						"  conflit cost change = "+ (conflicts_var_cost[variant] - conflicts_old_cost)+
						", num_conflicts diff = "+variant_conflicts_stats[variant].getNumConflicts());
			}
			
			
			if (debugLevel > 0) { // -1) {
				variant_conflicts_stats[variant].printConflictSummary(
						"Conflicts difference after resolution:", true, true, false);
			}
		}
		// How to compare? 1 attempt: none of the conflicts get worse, some get better or cost goes down
		int best_variant = -1;
		int best_ignore_conflicts = -1;
//		int [][] num_better_worse = new int [neibs_vars.length][2];
		int [] num_var_conflicts = new int [neibs_vars.length];
		
		
		for (int variant = 0; variant < neibs_vars.length; variant ++){
			if (variant_costs_diff[variant] <= dblTriLoss) {
				if  (conflicts_var_cost[variant] < conflicts_old_cost) {
					if ((best_variant < 0) || (conflicts_var_cost[variant] < conflicts_var_cost[best_variant])) {
						best_variant = variant;
					}
				}
				if ((best_ignore_conflicts < 0) || (variant_costs_diff[variant] < variant_costs_diff[best_ignore_conflicts])){
					best_ignore_conflicts = variant;
				}
			}
		}

		if (debugLevel > 1){
			System.out.println("resolveStarConflict(): for tile "+nsTile);
		}
		
		if ((best_variant < 0) && newConfl && (variant_costs_diff[best_ignore_conflicts] < 0)){ // should be cost improvement
			best_variant = best_ignore_conflicts;
			if (debugLevel > -1) {
				System.out.println("resolveMultiTriangularConflict(): conflicts increase but cost decreases "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +" of "+ neibs_vars.length+" variants");
			}
		}
		if ((best_variant < 0) || (variant_costs_diff[best_variant] > dblTriLoss)){
			if (debugLevel > -1) {
				System.out.println("resolveMultiTriangularConflict(): FAILED find a sutable solution for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +" of "+ neibs_vars.length+" variants");
				return false;
			}
		} else {
			if (debugLevel > -1) {
				System.out.println("resolveStarConflict(): SUCCESS to find a sutable solution for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +". Of "+ neibs_vars.length+" variants - use variant # " + best_variant+
						"  conflit cost change = "+ (conflicts_var_cost[best_variant] - conflicts_old_cost)+
						" cost difference (negative) = "+variant_costs_diff[best_variant]+
						" num conflict reductions = "+(-num_var_conflicts[best_variant]));
				variant_conflicts_stats[best_variant].printConflictSummary(
						"Conflicts number change per type: ",
						true, // use_all,
						true, //use_odo,
						true); // use_ood);
				iconflicts.printConflictSummary(
						"Conflicts before resolution: ",
						true, // use_all,
						true, //use_odo,
						true); // use_ood);
				// update statistics
				conflict_stats.addConflicts(variant_conflicts_stats[best_variant]);

				
				// update conflict
				for (int i = 0; i < nsTiles.length; i++){
					conflicts[nsTiles[i]]=  variant_conflicts[best_variant][i];
				}
				
				
				// apply resolution
				for (int i = 0; i < mod_supertiles.length; i++){
					if (debugLevel > 1){
						System.out.println("resolveStarConflict(): nsTile = "+nsTile+ "mod_supertiles["+i+"]="+mod_supertiles[i]);
					}
					if (neibs_vars[best_variant][i] != null) {
						for (int nl = 0; nl < neibs_vars[best_variant][i].length; nl ++){
							if (debugLevel > 1){
								System.out.println("resolveStarConflict(): nl= = "+nl);
							}
							if (neibs_vars[best_variant][i][nl] != null){
								planes[mod_supertiles[i]][nl].setNeibBest(neibs_vars[best_variant][i][nl]);
							}
						}
					}
				}
				// recalculate starValueWeights for and around tiles with modified  neighbors (no outside connections changed )nsTiles
				updateStarValueStrength(
						nsTiles,          // final int []         mod_supertiles,
						orthoWeight,      // final double         orthoWeight,
						diagonalWeight,   // final double         diagonalWeight,
						starPwr,          // final double         starPwr,    // Divide cost by number of connections to this power
						starSteps,        // final int            steps,
						planes,           // final TilePlanes.PlaneData [][] planes,
						preferDisparity); // final boolean        preferDisparity)
			}
		}
		return true;
	}
	
	public void calcStarValueStrength(
			final double         orthoWeight,
			final double         diagonalWeight,
			final double         starPwr,    // Divide cost by number of connections to this power
			final int            steps,
			final TilePlanes.PlaneData [][] planes,
			final boolean        preferDisparity)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int nStiles =       stilesX * stilesY; 
		final TileSurface.TileNeibs tnSurface = tileSurface.new TileNeibs(stilesX, stilesY);
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					ConnectionCosts connectionCosts = new ConnectionCosts(
							orthoWeight,    // double         orthoWeight,
							diagonalWeight, // double         diagonalWeight,
							starPwr,        // double         starPwr,    // Divide cost by number of connections to this power
							steps,          // int            steps,
							planes,         // TilePlanes.PlaneData [][] planes,
							tnSurface,      // TileSurface.TileNeibs tnSurface,
							preferDisparity); // boolean preferDisparity)
					int [] mod_supertiles = new int[1];
					for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
						if ( planes[nsTile] != null) {
							mod_supertiles[0] = nsTile;
							connectionCosts.initConnectionCosts(mod_supertiles);
							double [][][] val_weights = connectionCosts.getValWeights();
							for (int np = 0; np < planes[nsTile].length; np++){ // nu
								if (planes[nsTile][np] != null) {
									planes[nsTile][np].setStarValueWeight(val_weights[0][np]);
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}

	public void updateStarValueStrength(
			final int []         mod_supertiles,
			final double         orthoWeight,
			final double         diagonalWeight,
			final double         starPwr,    // Divide cost by number of connections to this power
			final int            steps,
			final TilePlanes.PlaneData [][] planes,
			final boolean        preferDisparity)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final TileSurface.TileNeibs tnSurface = tileSurface.new TileNeibs(stilesX, stilesY);
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					ConnectionCosts connectionCosts = new ConnectionCosts(
							orthoWeight,    // double         orthoWeight,
							diagonalWeight, // double         diagonalWeight,
							starPwr,        // double         starPwr,    // Divide cost by number of connections to this power
							steps,          // int            steps,
							planes,         // TilePlanes.PlaneData [][] planes,
							tnSurface,      // TileSurface.TileNeibs tnSurface,
							preferDisparity); // boolean preferDisparity)
					int [] mod_supertile = new int[1];
					for (int isTile = ai.getAndIncrement(); isTile < mod_supertiles.length; isTile = ai.getAndIncrement()) {
						int nsTile = mod_supertiles[isTile];
						if ((nsTile >= 0) && ( planes[nsTile] != null)) {
							mod_supertile[0] = nsTile;
							connectionCosts.initConnectionCosts(mod_supertile);
							double [][][] val_weights = connectionCosts.getValWeights();
							for (int np = 0; np < planes[nsTile].length; np++){ // nu
								if (planes[nsTile][np] != null) {
									planes[nsTile][np].setStarValueWeight(val_weights[0][np]);
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}
	
	
	//	
	/**
	 * Generate variants for changing connections while preserving number of connections between each pair of tiles
	 * each variant is free of own conflicts, but may lead to conflicts on other layers or tiles
	 * These variants should be checked on number of conflicts and overall cost to make sure conflict resolution leads
	 * to improvement - either in number of conflicts or global connection cost
	 * @param nsTile linescan tile index of the supertile
	 * @param nl1 first layer at nsTile to have triangle (2 ortho and a diagonal between them) conflict
	 * @param nl2 second layer at nsTile to have a triangle conflict, nl2 > nl1
	 * @param dir_mask start directions for the triangle conflicts 1 - N, 2 - E, 4 - S, 8 - W (starting from nl1),
	 * next 4 bits - similar fro starting from nl2 
	 * @param tnSurface ileSurface.TileNeibs instance to generate indices of the tiles at certain direction from the current
	 * @param debugLevel debug level
	 * @return [variant][5][nl][8] [5] = 4 directions and center ([4]), some may be null - should be removed by the caller
	 */

	public int [][][][] getTriangularResolutionVariants(
			int nsTile,
			int nl1,
			int nl2,
			int dir_mask,
			TileSurface.TileNeibs tnSurface,
			int        debugLevel)
	{
		int [][] layers0 = new int[5][2]; // 4 - center, 0 - N, 1 - E, 2 - S, 3 - w. 4 has 2, each other 0,1 or 2 layers
		for (int i = 0; i < (layers0.length - 1); i++){
			layers0[i][0] = -1;
			layers0[i][1] = -1;
		}
		layers0[4][0] = nl1;
		layers0[4][1] = nl2;
		int [] neibs1 = planes[nsTile][nl1].getNeibBest();
		int [] neibs2 = planes[nsTile][nl2].getNeibBest();
		for (int dir4 = 0; dir4 < 4; dir4++){
			if ((dir_mask & (1 << dir4)) != 0) {
				layers0[dir4][0] = neibs1[2 * dir4];
				layers0[(dir4 + 1) % 4][1] = neibs2[(2 * (dir4 + 1)) % 8];
			}
			if ((dir_mask & (16 << dir4)) != 0) {
				layers0[(dir4 + 1) % 4][0] = neibs1[(2 * (dir4 + 1)) % 8];
				layers0[dir4][1] = neibs2[2 * dir4];
			}
		}
		int [][] layers = new int[5][];
		int [] nsTiles = new int [5];
		nsTiles[4] = nsTile;
		for (int i = 0; i < 5; i++){
			if ((layers0[i][0] >=0 ) || (layers0[i][1] >=0 )) {
				if ((layers0[i][0] >=0 ) && (layers0[i][1] >=0 )){
					layers[i] = layers0[i].clone();
				} else {
					layers[i] = new int[1];
					layers[i][0] = (layers0[i][0] >=0 ) ? layers0[i][0] : layers0[i][1]; 
				}
				
			} else {
				layers[i] = new int[0];
				nsTiles[i] = -1;
			}
		}
		for (int dir4 = 0; dir4 < 4; dir4++){
			nsTiles[dir4] = tnSurface.getNeibIndex(nsTile, 2 * dir4);
		}
		int [] diag_connections = new int [4];// number of variants for diagonal connections (1 or 2 where applicable)
		for (int dir4 = 0; dir4 < 4; dir4 ++){
			int dir4_next = (dir4 + 1) % 4;
			int dir_diag = (2 * dir4 + 3) % 8;
			diag_connections[dir4] = 0;
			for (int il1 = 0; il1 < layers[dir4].length; il1++){
				int [] neibs_s = planes[nsTiles[dir4]][layers[dir4][il1]].getNeibBest();
				for (int il2 = 0; il2 < layers[dir4_next].length; il2++){
					if (neibs_s[dir_diag] == layers[dir4_next][il2]) diag_connections[dir4] ++;
				}
			}
		}
		int num_vars = 1;
		for (int dir4 = 0; dir4 < 4; dir4 ++){
			if (layers[dir4].length > 0) num_vars *=2;
			// number of ways to connect diagonals, preserving same number of links
			if ((layers[dir4].length == 2) && (layers[(dir4 + 1) % 4].length == 2) && (diag_connections[dir4] == 1)) num_vars *= 2;
		}
		if (debugLevel > -1){
			System.out.println("getTriangularResolutionVariants() nsTile = "+nsTile+" variants: "+num_vars+" (may be actually less)");
		}
		int [][][][] neib_vars = new int [num_vars][][][]; // [5][][];
		int [][][] neib_var_zero = new int [nsTiles.length][][];
		for (int i = 0; i < nsTiles.length; i++){
			if (layers[i].length > 0){
				neib_var_zero[i] = new int [planes[nsTiles[i]].length][];
				for (int nl = 0; nl < planes[nsTiles[i]].length; nl++){
					if (planes[nsTiles[i]][nl] != null){
						neib_var_zero[i][nl] = planes[nsTiles[i]][nl].getNeibBest(); 
					}
				}
			}
		}
		 // [i][0] - which of the neighbors layers to connect to nl1, [i][1] which of the diagonal connections (if choice) to keep
		boolean [][] selection = new boolean [4][2];
		int var_num = 0;
		while(true){
			int [][][] neib_var = neib_var_zero.clone();
			for (int i = 0; i < neib_var_zero.length; i++) if (neib_var_zero[i] != null){
				neib_var[i] = neib_var_zero[i].clone();
				for (int nl = 0; nl < neib_var_zero[i].length; nl++) if (neib_var_zero[i][nl] != null){
					neib_var[i][nl] = neib_var_zero[i][nl].clone();
				}			
			}			
			// use current selection to build neib_vars[var_num]. Increment [var_num] if OK, not to skip invalid variants
			// start with connecting no nl1, then - to nl2 and finally - diagonals
			good_variant:
			{
				for (int dir4 = 0; dir4 < 4; dir4 ++){
					int dir4_next = (dir4 + 1) % 4;
					int dir_f = 2 * dir4;
					int dir_r = (dir_f +4 ) % 8;
					int dir12_f = (2 * dir4 + 3 ) % 8;
					int dir12_r = (dir12_f + 4 ) % 8;

					if (layers[dir4].length > 0) {
						if (layers[dir4].length > 1){
							neib_var[4][nl1][dir_f] = selection[dir4][0] ? layers[dir4][1] : layers[dir4][0]; 
							neib_var[4][nl2][dir_f] = selection[dir4][0] ? layers[dir4][0] : layers[dir4][1]; 
							neib_var[dir4][layers[dir4][0]][dir_r] = selection[dir4][0] ? nl2 : nl1; 
							neib_var[dir4][layers[dir4][1]][dir_r] = selection[dir4][0] ? nl1 : nl2; 

						} else { // only 1 connection
							int other_layer = (neib_var_zero[4][nl1][dir_f] == layers[dir4][0]) ? neib_var_zero[4][nl2][dir_f] : neib_var_zero[4][nl1][dir_f];
							neib_var[4][nl1][dir_f] = selection[dir4][0] ? other_layer : layers[dir4][0]; 
							neib_var[4][nl2][dir_f] = selection[dir4][0] ? layers[dir4][0] : other_layer; 
							neib_var[dir4][layers[dir4][0]][dir_r] = selection[dir4][0] ? nl2 : nl1;
							if (other_layer >= 0) { // only if it is connected
								neib_var[dir4][other_layer][dir_r] = selection[dir4][0] ? nl1 : nl2; // null
							}
						}
						// diagonal connections - some make this variant invalid
						if (layers[dir4_next].length > 0) { // only add and verify connections if there are both radial connections
							int last_index_1 = (layers[dir4].length > 1) ? 1 :0;
							int last_index_2 = (layers[dir4_next].length > 1) ? 1 :0;
							int [][] neibs1_x = {neib_var_zero[dir4][layers[dir4][0]], // layers connected the first/single tile at dir4 from the center
									neib_var_zero[dir4][layers[dir4][last_index_1]]}; // layers connected the second/single tile at dir4 from the center
							int [][] neibs2_x = {
									neib_var_zero[dir4_next][layers[dir4_next][0]], // layers connected the first/single tile at dir4_next from the center
									neib_var_zero[dir4_next][layers[dir4_next][last_index_2]]}; // layers connected the last/single tile at dir4_next from the center
							// if there is a single variant. ( already connected)
							//if ((layers[dir4].length => 1)))
							if ((layers[dir4].length == 1) && (layers[dir4_next].length == 1)){ // 1 + 1
								if (selection[dir4][0] == selection[dir4_next][0]) continue; // no more checks - connected or not
								if (neibs1_x[0][dir12_f] != layers[dir4_next][0]) continue;     // not connected - keep whatever connections they had
								// connected, but radials are to the to different layers - invalid 
								break good_variant;
							} else if ((layers[dir4].length == 2) && (layers[dir4_next].length == 2)){ // 2 + 2
								if      (diag_connections[dir4] == 0) continue;
								// want to connect layers[dir4][0] (in the direction dir12_f) to :
								int want_10 = (selection[dir4][0] ^ selection[dir4_next][0]) ? layers[dir4_next][1] : layers[dir4_next][0];
								// want to connect layers[dir4][1] (in the direction dir12_f) to :
								int want_11 = (selection[dir4][0] ^ selection[dir4_next][0]) ? layers[dir4_next][0] : layers[dir4_next][1];
								// want to connect layers[dir4_next][0] (in the direction dir12_r) to :
								int want_20 = (selection[dir4][0] ^ selection[dir4_next][0]) ? layers[dir4][1] : layers[dir4][0];
								// want to connect layers[dir4_next][1] (in the direction dir12_rf) to :
								int want_21 = (selection[dir4][0] ^ selection[dir4_next][0]) ? layers[dir4][0] : layers[dir4][1];

								if (diag_connections[dir4] == 2) {
									neib_var[dir4][layers[dir4][0]][dir12_f] = want_10;
									neib_var[dir4][layers[dir4][1]][dir12_f] = want_11;
									neib_var[dir4_next][layers[dir4_next][0]][dir12_r] = want_20; 
									neib_var[dir4_next][layers[dir4_next][1]][dir12_r] = want_21; 
								} else { // 1 connection, (possibly) 2 variants
									// if correct connection, skip second variant. If not - selection[dir4][1]==false - swap first pair, true - the second
									if ((want_10 == neibs1_x[0][dir12_f]) || (want_11 == neibs1_x[1][dir12_f])){
										// Already connected correctly, so just a single variant - skip the second one
										if (!selection[dir4_next][1]) continue;
										break good_variant; // prevent duplicate
									}
									if (!selection[dir4_next][1]) { 
										neib_var[dir4][layers[dir4][0]][dir12_f] = neibs1_x[1][dir12_f];
										neib_var[dir4][layers[dir4][1]][dir12_f] = neibs1_x[0][dir12_f];
										// back links only if they existed
										if (neibs1_x[0][dir12_f] >= 0) neib_var[dir4_next][neibs1_x[0][dir12_f]][dir12_r] = layers[dir4][1]; 
										if (neibs1_x[1][dir12_f] >= 0) neib_var[dir4_next][neibs1_x[1][dir12_f]][dir12_r] = layers[dir4][0]; 
									} else {
										neib_var[dir4_next][layers[dir4_next][0]][dir12_r] = neibs2_x[1][dir12_r];
										neib_var[dir4_next][layers[dir4_next][1]][dir12_r] = neibs2_x[0][dir12_r];
										// back links only if they existed
										if (neibs2_x[0][dir12_r] >= 0) neib_var[dir4_next][neibs2_x[0][dir12_r]][dir12_f] = layers[dir4_next][1]; 
										if (neibs2_x[1][dir12_r] >= 0) neib_var[dir4_next][neibs2_x[1][dir12_r]][dir12_f] = layers[dir4_next][0]; 
									}
								}
							}  else {// remaining 2 + 1 and 1 + 2
								if (diag_connections[dir4] == 0) continue; // no links existed - nothing to do
								// here there was one link
								int conn_indx = (selection[dir4][0] ^ selection[dir4_next][0]) ? 1 : 0;
								int other_indx = 1 - conn_indx;
								if (layers[dir4].length == 1){
									if (neibs1_x[0][dir12_f] == layers[dir4_next][conn_indx]){
										continue; // already matched
									} else { // swap at dir4_next
										neib_var[dir4]     [layers[dir4][0]][dir12_f] = layers[dir4_next][conn_indx];
										neib_var[dir4_next][layers[dir4_next][conn_indx]][dir12_r] = layers[dir4][0];
										neib_var[dir4_next][layers[dir4_next][other_indx]][dir12_r] = neibs2_x[conn_indx][dir12_r];
										if (neibs2_x[conn_indx][dir12_r] >= 0) neib_var[dir4][neibs2_x[conn_indx][dir12_r]][dir12_f] = layers[dir4_next][other_indx];  
									}
								} else {
									if (neibs2_x[0][dir12_r] == layers[dir4][conn_indx]){
										continue; // already matched
									} else { // swap at dir4
										neib_var[dir4_next][layers[dir4_next][0]]         [dir12_r] = layers[dir4]     [conn_indx];
										neib_var[dir4]     [layers[dir4]     [conn_indx]] [dir12_f] = layers[dir4_next][0];
										neib_var[dir4]     [layers[dir4]     [other_indx]][dir12_f] = neibs1_x[conn_indx][dir12_f];
										if (neibs1_x[conn_indx][dir12_f] >= 0) neib_var[dir4_next][neibs1_x[conn_indx][dir12_f]][dir12_r] = layers[dir4][other_indx];  
									}
								}
							}

						}
					}
				}
				neib_vars[var_num++] = neib_var;
			}
			// calculate next selection
			next_sel:{
				for (int dir4 = 3; dir4 >= 0; dir4--){
					if (!selection[dir4][1] && (layers[dir4].length == 2) && (layers[(dir4 + 1) % 4].length == 2) && (diag_connections[dir4] == 1)){
						selection[dir4][1] = true;
						break next_sel;
//					} else if (!selection[dir4][0] && (layers[dir4].length > 1)){
					} else if (!selection[dir4][0] && (layers[dir4].length > 0)){
						selection[dir4][0] = true;
						selection[dir4][1] = false;
						break next_sel;
					} else {
						selection[dir4][0] = false;
						selection[dir4][1] = false;
					}
				}
				break; // while(true) - all done
			}
		}
// trim neib_vars
		if (debugLevel > -1){
			System.out.println("getTriangularResolutionVariants() nsTile = "+nsTile+" ACTUAL variants: "+var_num);
		}
		if (var_num == neib_vars.length) return neib_vars;
		int [][][][] neib_vars_trimmed = new int [var_num][][][];
		for (int i = 0; i < var_num; i++ ) neib_vars_trimmed[i] = neib_vars[i];
		
		return neib_vars_trimmed; 
	}

	public void resolveConflicts(
			double     maxEigen,      // maximal eigenvalue of planes to consider
	  		boolean    conflDualTri,  // Resolve dual triangles conflict (odoodo)
	  		boolean    conflMulti,    // Resolve multiple odo triangles conflicts
	  		boolean    conflDiag,     // Resolve diagonal (ood) conflicts
	  		boolean    conflStar,     // Resolve all conflicts around a supertile 
			int        starSteps, // How far to look around when calculationg connection cost
			double     orthoWeight,
			double     diagonalWeight,
			double     starPwr, // Divide cost by number of connections to this power
			double     dblTriLoss, //  When resolving double triangles allow minor degradation (0.0 - strict)
			boolean    newConfl, // Allow more conflicts if overall cost is reduced
			boolean    preferDisparity, 
			int        debugLevel,
			int        dbg_X,
			int        dbg_Y)
			
	{
		calcStarValueStrength(
				orthoWeight,      // final double         orthoWeight,
				diagonalWeight,   // final double         diagonalWeight,
				starPwr,          // final double         starPwr,    // Divide cost by number of connections to this power
				starSteps,        // final int            steps,
				this.planes,      // final TilePlanes.PlaneData [][] planes,
				preferDisparity); // final boolean        preferDisparity)
		
		Conflicts iconflicts0 = new Conflicts(this); 
		int [][][] conflicts0 = iconflicts0.detectTriangularConflicts(
				1); // final int debugLevel)
		Conflicts conflicts0_stats =  new Conflicts(
				conflicts0,
				this,
				-1); // debugLevel); 	
		
		for (int pass = 0; pass < 10; pass ++) {
			int [] dual_tri_results =            {0,0};
			int [] multi_resoultion_results =    {0,0};
			int [] diagonal_resoultion_results = {0,0};	
			int [] conflict_star_results = 	     {0,0};
			if (conflDualTri) {
				dual_tri_results = resolveDualTriangularConflicts(
						conflicts0, // int [][][] conflicts,
						conflicts0_stats,
						maxEigen,
						starSteps, // How far to look around when calculationg connection cost
						orthoWeight, // double     orthoWeight,
						diagonalWeight, // double     diagonalWeight,
						starPwr,        // double     starPwr, // Divide cost by number of connections to this power
						dblTriLoss, // double     diagonalWeight,
						preferDisparity,
						debugLevel, // 1, //  final int debugLevel)
						dbg_X,
						dbg_Y);
				System.out.println("Pass "+(pass+1)+": dual_tri_results (success/failures) = "+dual_tri_results[0]+" / "+dual_tri_results[1]);
			}
			if (conflMulti) {
				multi_resoultion_results = resolveMultiTriangularConflicts(
						conflicts0, // int [][][] conflicts,
						conflicts0_stats,
						starSteps, // How far to look around when calculationg connection cost
						orthoWeight, // double     orthoWeight,
						diagonalWeight, // double     diagonalWeight,
						starPwr,        // double     starPwr, // Divide cost by number of connections to this power
						dblTriLoss, // double     diagonalWeight,
						preferDisparity,
						debugLevel, // 1, //  final int debugLevel)
						dbg_X,
						dbg_Y);
				System.out.println("Pass "+(pass+1)+": multi_tri_results (success/failures) = "+multi_resoultion_results[0]+" / "+multi_resoultion_results[1]);
			}
			if (conflDiag) {
				diagonal_resoultion_results = resolveDiagonalTriangularConflicts(
						conflicts0, // int [][][] conflicts,
						conflicts0_stats,
						starSteps, // How far to look around when calculationg connection cost
						orthoWeight, // double     orthoWeight,
						diagonalWeight, // double     diagonalWeight,
						starPwr,        // double     starPwr, // Divide cost by number of connections to this power
						dblTriLoss, // double     diagonalWeight,
						preferDisparity,
						((pass == 3)? 1: 0), // debugLevel, // 1, //  final int debugLevel)
						dbg_X,
						dbg_Y);
				System.out.println("Pass "+(pass+1)+": resolveDiagonalTriangularConflicts (success/failures) = "+diagonal_resoultion_results[0]+" / "+diagonal_resoultion_results[1]);
			}			
			if (conflStar) {
				conflict_star_results = 	resolveStarConflicts(
						conflicts0, // int [][][] conflicts,
						conflicts0_stats,
						starSteps, // How far to look around when calculationg connection cost
						orthoWeight, // double     orthoWeight,
						diagonalWeight, // double     diagonalWeight,
						starPwr,        // double     starPwr, // Divide cost by number of connections to this power
						dblTriLoss, // double     diagonalWeight,
						newConfl, // Allow more conflicts if overall cost is reduced
						preferDisparity,
						debugLevel, // 1, //  final int debugLevel)
						dbg_X,
						dbg_Y);
				System.out.println("Pass "+(pass+1)+": resolveStarConflicts (success/failures) = "+conflict_star_results[0]+" / "+conflict_star_results[1]);
			}
			
			if (    (dual_tri_results[0] == 0) && 
					(multi_resoultion_results[0] == 0) &&
					(diagonal_resoultion_results[0] == 0) &&
					(conflict_star_results[0] ==       0)) {
				System.out.println("No more improvements");
				break;
			}
		}
		Conflicts conflicts1_stats =  new Conflicts(
				conflicts0,
				this,
				1); // -1); // debugLevel); 	

		conflicts1_stats.printConflictSummary("Detected conflicts (all):", true,false,false);
		conflicts1_stats.printConflictSummary("Detected conflicts (ortho-diag-ortho):", false, true,false);
		conflicts1_stats.printConflictSummary("Detected conflicts(ortho-ortho-diag):", false, false, true);
		
		iconflicts0.resetConflicts();
		conflicts0 = iconflicts0.detectTriangularConflicts(
				1); // final int debugLevel)
		conflicts1_stats.printConflictSummary("Recounted conflicts (all):", true,false,false);
		conflicts1_stats.printConflictSummary("Recounted conflicts (ortho-diag-ortho):", false, true,false);
		conflicts1_stats.printConflictSummary("Recounted conflicts(ortho-ortho-diag):", false, false, true);
		
		
		/*
		testResoveTriangle(
				clt_parameters.plWorstWorsening, // final double worst_worsening,
				clt_parameters.plWeakWorsening,  // final double worst_worsening,
				maxEigen,   // final double okMergeEigen,
				clt_parameters.plMaxWorldSin2,   // final double maxWorldSin2,
				clt_parameters.plDispNorm,
				clt_parameters.plMaxEigen,
				clt_parameters.plPreferDisparity,
				conflicts0, // int [][][] conflicts,
				1, // final int debugLevel)
				clt_parameters.tileX,
				clt_parameters.tileY);
      */				
		
	}
	
	
	
	
	
	
	
	
	public int [] resolveDualTriangularConflicts(
			int [][][] conflicts,
			Conflicts conflict_stats, // to be updated after applying resolution
			double     maxEigen, // maximal eigenvalue of planes to consider
			int        starSteps, // How far to look around when calculationg connection cost
			double     orthoWeight,
			double     diagonalWeight,
			double     starPwr, // Divide cost by number of connections to this power
			double     dblTriLoss, //  When resolving double triangles allow minor degradation (0.0 - strict)
			boolean    preferDisparity,
			int        debugLevel,
			int        dbg_X,
			int        dbg_Y)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int dbgTile = dbg_Y * stilesX + dbg_X;
		final TileSurface.TileNeibs tnSurface = tileSurface.new TileNeibs(stilesX, stilesY);
		final Conflicts iconflicts = new Conflicts(this); 

		int [] rslt = {0,0};
		for (int nsTile = 0; nsTile < conflicts.length; nsTile++) if (conflicts[nsTile] != null) {
			// conflicts may disappear after being fixed, recheck for null
			for (int nConfl = 0; (conflicts[nsTile] != null) && (nConfl < conflicts[nsTile].length); nConfl++){
				Conflict confl = new Conflict(nsTile, conflicts[nsTile][nConfl]); // getDualTriOrthoDiagOrthoConflicts
//				if (numDualTri(conflicts[nsTile][nConfl][2]) > 0) {
				if (confl.getDualTriOrthoDiagOrthoConflicts() > 0) {
					for (int dir4 = 0; dir4 < 4; dir4 ++ ){
						// conflicts may disappear after being fixed, recheck for null
						if ((conflicts[nsTile] != null) && (conflicts[nsTile][nConfl] != null) && ((conflicts[nsTile][nConfl][2] & (1 << dir4)) > 0) && ((conflicts[nsTile][nConfl][2] & (16 << dir4)) > 0)) {
							if (debugLevel > 0) { //1){
								System.out.println("resolveDualTriangularConflicts(): resolving dual triangular conflict for tile "+nsTile+
										", nl1 = "+conflicts[nsTile][nConfl][0]+
										", nl2 = "+conflicts[nsTile][nConfl][1]+
										", dir4 = "+dir4);
								if (nsTile == dbgTile){
									System.out.println("resolveDualTriangularConflicts(): nsTile == dbgTile");
								}
							}
							boolean [] dual_resol = resolveDualTriangularConflict(
									nsTile,
									conflicts[nsTile][nConfl][0],
									conflicts[nsTile][nConfl][1],
									dir4,
									maxEigen, // maximal eigenvalue of planes to consider
									tnSurface,
									preferDisparity,
									debugLevel + ((nsTile == dbgTile)? 1 : 0));
							if (dual_resol == null){
								if (debugLevel > 1) {
									System.out.println("resolveDualTriangularConflicts(): resolving dual triangular conflict for tile "+nsTile+
											", nl1 = "+conflicts[nsTile][nConfl][0]+
											", nl2 = "+conflicts[nsTile][nConfl][1]+
											", dir4 = "+dir4+"FAILED");
									rslt[1]++;
								}
							} else {
								if (debugLevel > 1) {
									System.out.println("resolveDualTriangularConflicts(): resolving dual triangular conflict for tile "+nsTile+
											", nl1 = "+conflicts[nsTile][nConfl][0]+
											", nl2 = "+conflicts[nsTile][nConfl][1]+
											", dir4 = "+dir4+" RESOLVED:");
									if (dual_resol[0])	System.out.println("swap links in dir4 = "+dir4);
									if (dual_resol[1])	System.out.println("swap links in dir4 = "+((dir4 + 1) % 4));
									if (dual_resol[2])	System.out.println("swap links in a digonal between dir4 = "+dir4+" and "+((dir4 + 1) % 4));
								}
								// calculate cost of swap
								int [] mod_supertiles = {
										nsTile,
										tnSurface.getNeibIndex(nsTile, 2 * dir4),
										tnSurface.getNeibIndex(nsTile, 2 * ((dir4 + 1) % 4))};
								int nl1 = conflicts[nsTile][nConfl][0];
								int nl2 = conflicts[nsTile][nConfl][1];

								ConnectionCosts connectionCosts = new ConnectionCosts(
										orthoWeight,
										diagonalWeight,
										starPwr,    // Divide cost by number of connections to this power
										starSteps,
										this.planes,
										tnSurface,
										preferDisparity);
								int [][][] neibs_prev = connectionCosts.initConnectionCosts(mod_supertiles);

								int [][][] neibs_prev_old = new int [mod_supertiles.length][][];
								double [][][]  val_weights = new double [mod_supertiles.length][][];
								updateConnectionsCost (
										mod_supertiles,      // int []         nsTiles,
										null,         // int [][][]     neibs_prev,
										neibs_prev_old,   // int [][][]     neibs,		// should be initialized at top dimension if neibs_prev==null	
										val_weights,  // double [][][]  val_weights, // should be initialized at top dimension if neibs_prev==null
										orthoWeight,
										diagonalWeight,
										starPwr,        // double     starPwr, // Divide cost by number of connections to this power
										tnSurface,
										preferDisparity,
										debugLevel);
										
								int [][][] neibs = neibs_prev.clone();
								for (int i = 0; i < neibs.length; i++) if (neibs_prev[i] != null){
									neibs[i] = neibs_prev[i].clone();
									for (int j = 0; j < neibs[i].length; j++) if (neibs_prev[i][j] != null){
										neibs[i][j] = neibs_prev[i][j].clone();
									}
								}
								//apply swaps:
//								int [] indices = {tile_index.get
								int dir1f = 2 * dir4; 
								int dir1r = (2 * dir4 + 4) % 8;
								int dir2f = (2 * dir4 + 3) % 8; 
								int dir2r = (2 * dir4 + 7) % 8;
								int dir3f = (2 * dir4 + 6) % 8; 
								int dir3r = (2 * dir4 + 2) % 8;
								int nl1_1 = neibs_prev[0][nl1][dir1f]; // layer on the first that was connected from nl1->N
								int nl2_1 = neibs_prev[0][nl2][dir1f]; // layer on the first that was connected from nl2->N
								int nl1_2 = neibs_prev[1][nl1_1][dir2f]; // layer on second that was connected from nl1->N->SE
								int nl2_2 = neibs_prev[1][nl2_1][dir2f]; // layer on second that was connected from nl2->N->SE
								//  2->0 was already swapped

								if (dual_resol[0])	{
									if (debugLevel > 0) System.out.println("swapping links in dir4 = "+dir4);
									neibs[0][nl1][dir1f] = neibs_prev[0][nl2][dir1f]; 
									neibs[0][nl2][dir1f] = neibs_prev[0][nl1][dir1f];
									neibs[1][nl1_1][dir1r] = neibs_prev[1][nl2_1][dir1r]; 
									neibs[1][nl2_1][dir1r] = neibs_prev[1][nl1_1][dir1r]; 
									
								}
								if (dual_resol[1])	{
									if (debugLevel > 0) System.out.println("swapping links in dir4 = "+((dir4 + 1) % 4));
									neibs[0][nl1][dir3r] = neibs_prev[0][nl2][dir3r]; 
									neibs[0][nl2][dir3r] = neibs_prev[0][nl1][dir3r];
									neibs[2][nl1_2][dir3f] = neibs_prev[2][nl2_2][dir3f]; 
									neibs[2][nl2_2][dir3f] = neibs_prev[2][nl1_2][dir3f]; 
								}
								if (dual_resol[2]){
									if (debugLevel > 0) System.out.println("swapping links in a digonal between dir4 = "+dir4+" and "+((dir4 + 1) % 4));
									neibs[1][nl1_1][dir2f] = neibs_prev[1][nl2_1][dir2f]; 
									neibs[1][nl2_1][dir2f] = neibs_prev[1][nl1_1][dir2f];
									
									neibs[2][nl1_2][dir2r] = neibs_prev[2][nl2_2][dir2r]; 
									neibs[2][nl2_2][dir2r] = neibs_prev[2][nl1_2][dir2r]; 
								}
								
								double cost_diff = 	connectionCosts.getConnectionsCostDiff(
										neibs,
										debugLevel);
								
								if (debugLevel > 0) {
									System.out.println("resolveDualTriangularConflicts(): resolving dual triangular conflict for tile "+nsTile+
											", nl1 = "+conflicts[nsTile][nConfl][0]+
											", nl2 = "+conflicts[nsTile][nConfl][1]+
											", dir4 = "+dir4+" improvement (negative diff) = "+cost_diff);
								}
								
								double cost_diff_old = 	updateConnectionsCost (
										mod_supertiles, // int []         nsTiles,
										neibs_prev,     // int [][][]     neibs_prev,
										neibs,          // int [][][]     neibs,		// should be initialized at top dimension if neibs_prev==null	
										val_weights,    // double [][][]  val_weights, // should be initialized at top dimension if neibs_prev==null
										orthoWeight,
										diagonalWeight,
										starPwr,        // double     starPwr, // Divide cost by number of connections to this power
										tnSurface,
										preferDisparity,
										debugLevel);
										
								if (debugLevel > 0) {
									System.out.println("resolveDualTriangularConflicts()_OLD: resolving dual triangular conflict for tile "+nsTile+
											", nl1 = "+conflicts[nsTile][nConfl][0]+
											", nl2 = "+conflicts[nsTile][nConfl][1]+
											", dir4 = "+dir4+" improvement (negative diff) = "+cost_diff_old);
								}
								if (cost_diff < dblTriLoss) {
									if (debugLevel > 1) {
										if (cost_diff > 0 ) {
											System.out.println("resolveDualTriangularConflicts(): Score is WORSE, but degrading is less than "+
													dblTriLoss+", considered OK.");
										} else {
											System.out.println("resolveDualTriangularConflicts(): SUCCESS");
										}
									}
									// See how this application will influence number of conflicts
									// All supertiles that may have different conflicts 
									int [] nsTiles = getInvolvedSupertiles( // first 3 entries will be mod_supertiles[]
											mod_supertiles,
											tnSurface);
									HashMap<Integer,Integer> replacement_tiles = new HashMap<Integer,Integer>();
									for (int i = 0; i < mod_supertiles.length; i++){
										replacement_tiles.put(mod_supertiles[i], new Integer(i));
									}
									int [][][] conflicts_old = new int [nsTiles.length][][];
									int [][][] conflicts_new = new int [nsTiles.length][][];
									for (int isTile = 0; isTile < nsTiles.length; isTile++){
										conflicts_old[isTile] = iconflicts.detectTriangularTileConflicts(
												nsTiles[isTile],   // int nsTile0,
												replacement_tiles, //HashMap<Integer,Integer> replacement_tiles, //
												neibs_prev, // int [][][] replacement_neibs,
												tnSurface); // TileSurface.TileNeibs tnSurface)
										conflicts_new[isTile] = iconflicts.detectTriangularTileConflicts(
												nsTiles[isTile],   // int nsTile0,
												replacement_tiles, //HashMap<Integer,Integer> replacement_tiles, //
												neibs, // int [][][] replacement_neibs,
												tnSurface); // TileSurface.TileNeibs tnSurface)
									}
									if (debugLevel > 1) {
										System.out.println("Involved supertiles:");
										for (int i = 0; i < nsTiles.length; i++){
											System.out.println(i+":"+nsTiles[i]);
										}
									}

									if (debugLevel > 1) {
										System.out.println("Calculating original conflicts");
									}

									iconflicts.addConflicts(conflicts_old,
											debugLevel - 1); // debugLevel);
									if (debugLevel > 1) {
										System.out.println("Calculating resolved conflicts");
									}
									Conflicts conflict_stats_new = new Conflicts(
											conflicts_new,
											this,
											debugLevel - 1); // debugLevel);
									conflict_stats_new.subConflicts(iconflicts);
									
									if (debugLevel > 0) {
										/*
										System.out.print("Conflicts difference after resolution:");
										printConflictSummary(conflict_stats_diff);
										*/
										conflict_stats_new.printConflictSummary(
												"Conflicts difference after resolution:", true, true, false);
										
									}
									if ((conflict_stats_new.getNumOrthoIncompat() >= 0) && (cost_diff > -dblTriLoss)){
										//getNumOrthoIncompat
										if (debugLevel > 1) {
											System.out.println("Number of incompatible triangles is not reduced and no significant cost reduction, resolution WILL NOT BE APPLIED");
										}
										rslt[1]++;
									} else {
										if (debugLevel > 0) {
											if (conflict_stats_new.getNumOrthoIncompat() < 0) {
												System.out.println("Number of incompatible triangles is reduced, resolution WILL BE APPLIED");
											} else {
												System.out.println("Number of incompatible triangles is not reduced but a good cost reduction, resolution WILL BE APPLIED");
											}
										}
										// update statistics
										conflict_stats.addConflicts(conflict_stats_new);

										// update conflict
										for (int i = 0; i < nsTiles.length; i++){
											conflicts[nsTiles[i]]=  conflicts_new[i];
										}
										// apply resolution
										for (int i = 0; i < mod_supertiles.length; i++){
											for (int nl = 0; nl < neibs[i].length; nl ++) if (neibs[i][nl] != null){
												planes[mod_supertiles[i]][nl].setNeibBest(neibs[i][nl]);
											}
										}
										rslt[0]++;
									}
								} else {
									if (debugLevel > 0) {
										System.out.println("resolveDualTriangularConflicts(): FAILURE");
									}
									rslt[1]++;
								}
							}
						}
					}
				}
			}
		}
		return rslt;
	}
	
	
	/**
	 *  Resolve dual triangle, where NL1 -> N -> SE -> W -> NL2 -> N -> SE -> W -> NL1 (or similar)
	 * @param nsTile supertile index
	 * @param nl1 first conflicting layer (lower of the too
	 * @param nl2 second conflicting layer
	 * @param dir4 start direction to (turn right)  0 - N, 1 - E, 2 - S, 3 - W
	 * @param maxEigen maximal sum of 2 eigenvalues of the tri-plane combinations
	 * @param tnSurface TileNeibs instance to navigate tile index and control array borders
	 * @param preferDisparity - the first eigenvalue/vector is the most disparity-like
	 *                          (false - smallest eigenvalue)
	 * @param debugLevel
	 * @return links to swap : [0] first in dir, [1] - first in dir + 1, [2] - other diagonal. Or null if failed
	 * total number of bits set is odd  
	 */
	public boolean [] resolveDualTriangularConflict(
			int nsTile,
			int nl1,
			int nl2,
			int dir4,
			double maxEigen, // maximal eigenvalue of planes to consider
			TileSurface.TileNeibs tnSurface,
			boolean preferDisparity,
			int debugLevel)
	{
		int [] neibs1 = planes[nsTile][nl1].getNeibBest();
		int [] neibs2 = planes[nsTile][nl2].getNeibBest();
		int dir1 = 2 * dir4;
		int dir2 = (2 * dir4 + 2) % 8;
		int nsTile1 = tnSurface.getNeibIndex(nsTile, dir1);
		int nsTile2 = tnSurface.getNeibIndex(nsTile, dir2);
		
		TilePlanes.PlaneData [][] tri_pd = {
				{planes[nsTile][nl1], planes[nsTile][nl2]},
				{planes[nsTile1][neibs1[dir1]], planes[nsTile1][neibs2[dir1]]},
				{planes[nsTile2][neibs1[dir2]], planes[nsTile2][neibs2[dir2]]}};
		TilePlanes.PlaneData [][] tri_merged = new TilePlanes.PlaneData [4][2]; // keep to the very end for debugging
		// initial connections: [0][0] -> [1][0] -> [2][1] -> [0][1] -> [1][1]-> [2][0] -> [0][0]
		// 4 variants for 3-supertile sets. get with smallest sum of the 2 main eigenvalues, as long as each
		// is below maxEigen
		if (debugLevel > 1) {
			System.out.println ("resolveDualTriangularConflict() + nsTile = "+nsTile+", nl1 = "+nl1+", nl2 = "+nl2+" dir4="+dir4);
		}
		double [][] ev_variants = new double [4][2];
		int best_pair = -1;
		double best_value = Double.NaN; 
		for (int variant = 0; variant < 4; variant++){
			TilePlanes.PlaneData [][] tri_sets = {
					{tri_pd[0][0], tri_pd[1][    (variant & 1)], tri_pd[2][    ((variant >> 1) & 1)]},
					{tri_pd[0][1], tri_pd[1][1 - (variant & 1)], tri_pd[2][1 - ((variant >> 1) & 1)]}};
			for (int nt = 0; nt < tri_sets.length; nt++){
				TilePlanes.PlaneData tri_plane = tri_sets[nt][0];
				for (int indx = 1; indx <  tri_sets[nt].length; indx++){
					TilePlanes.PlaneData other_plane = tri_sets[nt][0].getPlaneToThis(  // layer here does not matter
							tri_sets[nt][indx],
							debugLevel - 2); // debugLevel);
					tri_plane = tri_plane.mergePlaneToThis(
							other_plane,     // PlaneData otherPd,
							1.0,             // double    scale_other,
							false,           // boolean   ignore_weights,
							true,            // boolean   sum_weights,
							preferDisparity, 
							debugLevel - 2); // int       debugLevel)
				}
				ev_variants[variant][nt] = tri_plane.getValue();
				tri_merged[variant][nt] = tri_plane;
			}
			if ((best_pair < 0) || ((ev_variants[variant][0] + ev_variants[variant][1]) < best_value)){
				best_value = ev_variants[variant][0] + ev_variants[variant][1];
				best_pair = variant;
			}
			if (debugLevel > 1) {
				System.out.println ("variant="+variant+" value1="+ev_variants[variant][0]+
						" value2="+ev_variants[variant][1]+" sum value = "+(ev_variants[variant][0]+ev_variants[variant][1]));
			}
		}
		if (debugLevel > 1) {
			System.out.println ("resolveDualTriangularConflict() + nsTile = "+nsTile+", nl1 = "+nl1+", nl2 = "+nl2+" dir4="+dir4+" best variant ="+ best_pair);
		}
		if (best_value < maxEigen) {
			boolean [] rslt = {
					((best_pair & 1) != 0),
					((best_pair & 2) != 0),
					false};
			rslt[2] = !rslt[0] ^ rslt[1];
			return rslt;
		}
		return null;
	}	
	
	/**
	 * Calculate main eigenvalue of the current plane and all connected ones - used to estimate advantage of connection swap
	 * @param nsTile supertile index
	 * @param nl surface layer
	 * @param neibs array of 8 neighbors layers (N,NE,...NW), -1 - not connected
	 * @param orthoWeight multiply contribution of ortho neighbors
	 * @param diagonalWeight  multiply contribution of diagonal neighbors
	 * @param diagonalWeight  divide value by number of connections to this power (if !=0)
	 * @param tnSurface TileNeibs instance to navigate tile index and control array borders
	 * @param preferDisparity - the first eigenvalue/vector is the most disparity-like
	 *                          (false - smallest eigenvalue)
	 * @param debugLevel
	 * @return a pair of eigenvalue of the combine plane and its weight
	 */
	
	public double [] getStarValueWeight(
			int    nsTile,
			int    nl,
			int [] neibs,
			double orthoWeight,
			double diagonalWeight,
			double starPwr,    // Divide cost by number of connections to this power
			TileSurface.TileNeibs tnSurface,
			boolean preferDisparity,
			int    debugLevel)
	{
		TilePlanes.PlaneData merged_plane = planes[nsTile][nl];
		for (int dir = 0; dir < 8; dir++){
			if (neibs[dir] >= 0){
				double other_weight = ((dir & 1) != 0) ? diagonalWeight : orthoWeight;
				TilePlanes.PlaneData other_plane = merged_plane.getPlaneToThis(  // layer here does not matter
						planes[tnSurface.getNeibIndex(nsTile, dir)][neibs[dir]],
						debugLevel - 1); // debugLevel);
				merged_plane = merged_plane.mergePlaneToThis(
						other_plane,     // PlaneData otherPd,
						other_weight,    // double    scale_other,
						false,           // boolean   ignore_weights,
						true,            // boolean   sum_weights,
						preferDisparity, 
						debugLevel - 1); // int       debugLevel)
			}
		}
		double [] value_weight = {merged_plane.getValue(),merged_plane.getWeight()};
		if (starPwr != 0){
			value_weight[0] /= (Math.pow((planes[nsTile][nl].getNumNeibBest() + 1.0), starPwr));
		}
		return value_weight;
	}

	/**
	 * Calculate array of supertile indices that need to have connection cost recalculated when they are updated
	 * first entries of the result will be the same in input array
	 * @param mod_supertiles array of supertile indices that will be modified
	 * @param tnSurface TileNeibs instance to navigate tile index and control array borders
	 * @return array of supertile indices to watch connection cost
	 */
	
	
	public int [] getInvolvedSupertiles(
			int [] mod_supertiles,
			TileSurface.TileNeibs tnSurface)
	{
		HashSet<Integer> stiles_set = new HashSet<Integer>();
		for (int i = 0; i < mod_supertiles.length; i++){
			stiles_set.add(new Integer(mod_supertiles[i]));
			for (int dir = 0; dir < 8; dir++){
				Integer nsTile = tnSurface.getNeibIndex(mod_supertiles[i], dir);
				if (nsTile >= 0) stiles_set.add (nsTile);
			}
		}
		int [] stiles = new int [stiles_set.size()];
		int indx = 0;
		for (; indx < mod_supertiles.length; indx++){
			stiles_set.remove(new Integer(mod_supertiles[indx]));
			stiles[indx] = mod_supertiles[indx];
		}
		for (Integer nsTile: stiles_set){
			stiles[indx++] = nsTile;
		}
		return stiles;
	}
	
	/**
	 * Update cost of involved connections after swapping some of them. Total number of links should probably stay
	 *  the same.
	 * @param nsTiles array of indices of the supertiles to watch
	 * @param neibs_prev 3-d array [per supertile index][per layer][per direction] of previous connections for which
	 * the cost was already calculated. Or null - in that case established connections will be used and returned in
	 * neibs array
	 * @param neibs new connections to calculate cost for. If neibs_prev == null, this array should be initialized
	 * by the caller as new int [nsTiles.length][][]
	 * @param val_weights array of values*weights and weights per supertile index, per layer of {val*weight, weight}
	 * @param orthoWeight multiply contribution of ortho neighbors
	 * @param diagonalWeight  multiply contribution of diagonal neighbors
	 * @param tnSurface TileNeibs instance to navigate tile index and control array borders
	 * @param preferDisparity - the first eigenvalue/vector is the most disparity-like
	 *                          (false - smallest eigenvalue)
	 * @param debugLevel
	 * @return difference between new cost and old cost, negative means improvement
	 */
	
	
	public double updateConnectionsCost (
			int []         nsTiles,
			int [][][]     neibs_prev,
			int [][][]     neibs,		// should be initialized at top dimension if neibs_prev==null	
			double [][][]  val_weights, // should be initialized at top dimension if neibs_prev==null
			double         orthoWeight,
			double         diagonalWeight,
			double         starPwr,    // Divide cost by number of connections to this power
			TileSurface.TileNeibs tnSurface,
			boolean preferDisparity,
			int    debugLevel)
	{
		if (neibs_prev == null){
//			neibs_prev = new int [nsTiles.length][][];
			for (int isTile = 0; isTile < nsTiles.length; isTile++){
				int nsTile = nsTiles[isTile];
				if (planes[nsTile] != null){
					neibs[isTile] = new int [planes[nsTile].length][];
					val_weights[isTile] = new double [planes[nsTile].length][];
					for (int nl = 0; nl < planes[nsTile].length; nl++) if ( planes[nsTile][nl] != null){
						neibs[isTile][nl] =  planes[nsTile][nl].getNeibBest();
						//getNumNeibBest()
						val_weights[isTile][nl] = new double[2];
					}
				}
			}
		}
		// calculate old cost
		double old_value = 0.0;
		double old_weight = 0.0; // should not change during update
		for (int isTile = 0; isTile < nsTiles.length; isTile++){
			if (val_weights[isTile] != null){
				for (int nl = 0; nl < val_weights[isTile].length; nl++) if ( val_weights[isTile][nl] != null){
					old_value += val_weights[isTile][nl][0] * val_weights[isTile][nl][1];
					old_weight += val_weights[isTile][nl][1];
				}
			}
		}
		// now re-calculate val_weights where neibs are different from neibs_prev
		for (int isTile = 0; isTile < nsTiles.length; isTile++){
			int nsTile = nsTiles[isTile];
			if (planes[nsTile] != null) {
				if (neibs[isTile] != null){
					for (int nl = 0; nl < planes[nsTile].length; nl++) if (planes[nsTile][nl] != null){
						if  (neibs[isTile][nl] != null) {
							boolean neibs_changed = false;
							if ((neibs_prev == null) || (neibs_prev[isTile] == null)  || (neibs_prev[isTile][nl] == null)){
								neibs_changed = true;
							} else {
								for (int dir = 0; dir < 8; dir++) if (neibs[isTile][nl][dir] != neibs_prev[isTile][nl][dir]){
									neibs_changed = true;
									break;
								}
							}
							if (neibs_changed){
								val_weights[isTile][nl] = getStarValueWeight(
										nsTile,
										nl,
										neibs[isTile][nl],
										orthoWeight,
										diagonalWeight,
										starPwr, // double         starPwr,    // Divide cost by number of connections to this power
										tnSurface,
										preferDisparity,
										-1); // debugLevel);
							}
						} else {
							val_weights[isTile][nl] = null;
						}
					}				
				} else {
					val_weights[isTile] = null;
				}
			}
		}
		// calculate new cost
		double new_value = 0.0;
		double new_weight = 0.0; // should not change during update
		for (int isTile = 0; isTile < nsTiles.length; isTile++){
			if (val_weights[isTile] != null){
				for (int nl = 0; nl < val_weights[isTile].length; nl++) if ( val_weights[isTile][nl] != null){
					new_value += val_weights[isTile][nl][0] * val_weights[isTile][nl][1];
					new_weight += val_weights[isTile][nl][1];
				}
			}
		}
		// Likely weight will never change except first run, but we will still normalize by weight
		if (old_weight != 0.0) old_value /= old_weight;
		if (new_weight != 0.0) new_value /= new_weight;
		return new_value - old_value; // negative - improvement
	}
	
	
	
	/**
	 * Find out which of the layers to use to resolve the conflict. Multiple conflict may need to
	 * be pre-filtered or split to prevent using simultaneously different neighbors for the same
	 * direction
	 * @param nsTile supertile index
	 * @param nl1 first conflicting layer (lower of the too
	 * @param nl2 second conflicting layer
	 * @param dir_mask direction bit mask: bit 0: nl1->N->SE->W->nl2, bit 1: nl1->E->SW->N->nl2,...,
	 * bit 4 : nl1->E->NW->S->nl2 (backwards of bit1), bit 5: nl1->S->NE->W->nl2 (backwards of bit1)
	 *  	 * @param rquality maximal degradation by merging (does not depend on the total weight)
	 * @param okMergeEigen if result eigenvalue of the merged planes is below, OK to bypass worst worsening
	 * @param maxWorldSin2 maximal square of the sine of the angle between the planes to allow merge (>= 1.0 - disable)
	 * @param maxEigen maximal eigenvalue of each of the merged planes
	 * @param tnSurface TileNeibs instance to navigate tile index and control array borders
	 * @param preferDisparity - the first eigenvalue/vector is the most disparity-like
	 *                          (false - smallest eigenvalue)
	 * @param debugLevel
	 * @return winning layer or -1 if not possible
	 * TODO: add more parameters - some thresholds 
	 */
	public int resolveTriangularConflict(
			int nsTile,
			int nl1,
			int nl2,
			int dir_mask,
			double rquality,
			double weakWorsening,
			double okMergeEigen,
			double maxWorldSin2,
			double dispNorm,
			double maxEigen, // maximal eigenvalue of planes to consider
			TileSurface.TileNeibs tnSurface,
			boolean preferDisparity,
			int debugLevel)
	{
		ArrayList<Point> surf_list = new ArrayList<Point>();
		int [] neibs1 = planes[nsTile][nl1].getNeibBest();
		int [] neibs2 = planes[nsTile][nl2].getNeibBest();
		for (int dir = 0; dir < 8; dir +=2){
			if ((dir_mask & (1 << (dir / 2))) != 0){
				int nsTile1 = tnSurface.getNeibIndex(nsTile, dir);
				int nsTile2 = tnSurface.getNeibIndex(nsTile, (dir + 2) % 8);
				
				surf_list.add(new Point(nsTile1,neibs1[dir]));
				surf_list.add(new Point(nsTile2,neibs2[(dir + 2) % 8]));
			}
			if ((dir_mask & (1 << ((dir / 2) + 4))) != 0){
				int nsTile1 = tnSurface.getNeibIndex(nsTile, dir);
				int nsTile2 = tnSurface.getNeibIndex(nsTile, (dir + 2) % 8);
				surf_list.add(new Point(nsTile1,neibs2[dir]));
				surf_list.add(new Point(nsTile2,neibs1[(dir + 2) % 8]));
			}
		} // TODO: do not use duplicates? Or keep to increase weight?
		TilePlanes.PlaneData merged_plane = null;
		TilePlanes.PlaneData this_plane1 = planes[nsTile][nl1];
		TilePlanes.PlaneData this_plane2 = planes[nsTile][nl2];
		for (Point p : surf_list) {
			TilePlanes.PlaneData other_plane = this_plane1.getPlaneToThis(  // layer here does not matter
					planes[p.x][p.y],
					debugLevel - 1); // debugLevel);
			if (merged_plane == null){
				merged_plane = other_plane; 
			} else {
				merged_plane = merged_plane.mergePlaneToThis(
						other_plane, // PlaneData otherPd,
						1.0,         // double    scale_other,
						false,       // boolean   ignore_weights,
						true, // boolean   sum_weights,
						preferDisparity, 
						debugLevel - 1); // int       debugLevel)
			}
		}
		TilePlanes.PlaneData [] these = {this_plane1, this_plane2};
		TilePlanes.PlaneData [] merged = new TilePlanes.PlaneData [2];

		// now find which is better and if any of them fits at all. Print other criteria also - maybe use them later
		boolean [] fit = {true,true}; // initially both fit
		double [] this_rq = new double[2];
		for (int np = 0; np < 2; np++){
			merged[np] = these[np].mergePlaneToThis(
					merged_plane, // PlaneData otherPd,
					1.0,          // double    scale_other,
					false,        // boolean   ignore_weights,
					true,         // boolean   sum_weights,
					preferDisparity, 
					debugLevel - 1); // int       debugLevel)

			// verify result plane smallest eigenvalue is OK  
			if ((maxEigen != 0.0) &&
					(merged[np].getValue() > corrMaxEigen(
							maxEigen,
							dispNorm,
							merged[np]))) {
				fit[np] = false;
				continue;
			}
			double w1 = merged_plane.getWeight();
			double w2 = merged[np].getWeight();
			this_rq[np] = mergeRQuality(
					merged_plane.getValue(), // double L1,
					merged[np].getValue(), // double L2,
					merged[np].getValue(), // double L,
					w1, // double w1,
					w2); // double w2)
			double this_rq_norm = this_rq[np];
			if ((w1 + w2) < weakWorsening) this_rq_norm *= (w1 + w2) / weakWorsening; // forgive more for weak planes

			// verify that eigenvalue worsening by the merge is OK
			if (    (this_rq_norm > rquality) &&
					(merged[np].getValue() > okMergeEigen)){
				fit[np] = false;
				continue;
			}
			// verify
			if ((maxWorldSin2 < 1.0) && (merged_plane.getWorldSin2(these[np]) > maxWorldSin2)) { // use these[np], no merged[np], invert comparison
				fit[np] = false;
				continue;
			}
			this_rq[np] /= (w1 + w2); // for comparison reduce this value for stronger planes 
			if (debugLevel >0) {
				System.out.println("nsTile="+nsTile+":"+np+", nsTile="+nsTile+":"+np+", this_rq="+this_rq[np]+
						" w1="+w1+" w2="+w2+
						" L1="+merged_plane.getValue()+" L2="+these[np].getValue()+" L="+merged[np].getValue());

				System.out.println("nsTile="+nsTile+":"+np+", world sin2 ="+
						merged_plane.getWorldSin2(these[np]));

				System.out.println("nsTile="+nsTile+":"+np+
						", world dist this="+ Math.sqrt(merged_plane.getWorldPlaneDist2(these[np]))+
						", world dist other="+Math.sqrt(these[np].getWorldPlaneDist2(merged_plane))+
						", world dist sum="+Math.sqrt(merged_plane.getWorldPlaneDist2(these[np])+
								these[np].getWorldPlaneDist2(merged_plane)));
			}
		}
		if (fit[0] || fit[1]){
			if (!fit[0]) return nl2;
			if (!fit[1]) return nl1;
			
			if (this_rq[0] < this_rq[1]) return nl1;
			else return nl2;
		}

		return -1;
	}
	
	
	
	/**
	 * Find mutual links between multi-layer planes for supertiles. requires that for each plane there are calculated smalles eigenvalues
	 * for merging with each plane for each of 8 neighbors 
	 * @param rquality maximal degradation by merging (does not depend on the total weight)
	 * @param okMergeEigen if result eigenvalue of the merged planes is below, OK to bypass worst worsening
	 * @param maxWorldSin2 maximal square of the sine of the angle between the planes to allow merge (>= 1.0 - disable)
	 * @param maxEigen maximal eigenvalue of each of the merged planes
	 * @param minWeight minimal weight of each of the planes
	 * @param debugLevel debug level
	 * @param dbg_X debug supertile X coordinate
	 * @param dbg_Y debug supertile Y coordinate
	 */
	public void selectNeighborPlanesMutual(
			final double rquality,
			final double weakWorsening,
			final double okMergeEigen,
			final double maxWorldSin2,
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
								if (this_plane != null) {
									this_plane.initNeibBest();
								}
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
									if ((planes[nsTile0][np0] != null) && (planes[nsTile0][np0].getMergedValue(dir) != null)){
										int l = planes[nsTile0][np0].getMergedValue(dir).length;
										if (l > num_other_planes)num_other_planes = l;
									}
								}
								if (num_other_planes > 0){ // will eliminate bad margins
									int np_min = LOWEST_PLANE(num_other_planes);
									
									boolean [] this_matched = new boolean [planes[nsTile0].length];
									boolean [] other_matched = new boolean [num_other_planes];
									int num_pairs = this_matched.length - np0_min;
									if ((other_matched.length - np_min) < num_pairs) num_pairs = other_matched.length - np_min;

									for (int pair = 0; pair < num_pairs; pair ++){
										int [] best_pair = {-1,-1};
										double best_rqual = Double.NaN;
										for (int np0 = np0_min; np0 < this_matched.length; np0++) if (planes[nsTile0][np0] != null){
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
															(planes[nsTile][np] != null) &&
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
														if ((this_rq_norm <= rquality) ||(merge_ev[np] <= okMergeEigen)) {
															if ((maxWorldSin2 >= 1.0) || (planes[nsTile0][np0].getWorldSin2(planes[nsTile][np]) <=maxWorldSin2)) {
																this_rq /= (w1 + w2); // for comparison reduce this value for stronger planes 
																if (Double.isNaN(best_rqual) || (this_rq < best_rqual)){ // OK if Double.isNaN(this_rq[np])
																	best_rqual = this_rq;
																	best_pair[0]= np0;
																	best_pair[1]= np;
																}
															}
														}
														if (debugLevel > 0){
															if ((merge_ev[np] < 0.4) && (w1 > 1.0)  && (w2 > 1.0) ){
																System.out.println("nsTile0="+nsTile0+":"+np0+", nsTile="+nsTile+":"+np+", this_rq="+this_rq+
																		", this_rq*(w1+w2)="+(this_rq * (w1 + w2))+
																		" w1="+w1+" w2="+w2+
																		" L1="+planes[nsTile0][np0].getValue()+" L2="+planes[nsTile][np].getValue()+" L="+merge_ev[np]);
															}
														}
														if (dl >0) {
															System.out.println("nsTile0="+nsTile0+":"+np0+", nsTile="+nsTile+":"+np+", this_rq="+this_rq+
																	" w1="+w1+" w2="+w2+
																	" L1="+planes[nsTile0][np0].getValue()+" L2="+planes[nsTile][np].getValue()+" L="+merge_ev[np]);
															System.out.println("nsTile0="+nsTile0+":"+np0+", nsTile="+nsTile+":"+np+", world sin2 ="+
																	planes[nsTile0][np0].getWorldSin2(planes[nsTile][np]));
															System.out.println("nsTile0="+nsTile0+":"+np0+", nsTile="+nsTile+":"+np+
																	", world dist this="+ Math.sqrt(planes[nsTile0][np0].getWorldPlaneDist2(planes[nsTile][np]))+
																	", world dist other="+Math.sqrt(planes[nsTile][np].getWorldPlaneDist2(planes[nsTile0][np0]))+
																	", world dist sum="+Math.sqrt(planes[nsTile0][np0].getWorldPlaneDist2(planes[nsTile][np])+
																			planes[nsTile][np].getWorldPlaneDist2(planes[nsTile0][np0])));
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
									for (int np = 0; np < this_matched.length; np++) if (planes[nsTile0][np] != null){
										if (!this_matched[np]){
											planes[nsTile0][np].setNeibBest(dir,-1);
										}
									}
//									for (int np = 0; np < other_matched.length; np++) if (planes[nsTile][np] != null){
									for (int np = 0; np < other_matched.length; np++) if ((planes[nsTile][np] != null) && (planes[nsTile][np].getWeight() > 0.0)){ // disregard 0-weight planes
										if (!other_matched[np]){
											planes[nsTile][np].setNeibBest(dir + 4,-1);
										}
									}
									if (dl >0) {
										for (int np = 0; np < this_matched.length; np++) if (planes[nsTile0][np] != null){
											int [] bn = planes[nsTile0][np].getNeibBest();
											System.out.println("nsTile0="+nsTile0+":"+np+" : ["+
											bn[0]+","+bn[1]+","+bn[2]+","+bn[3]+","+bn[4]+","+bn[5]+","+bn[6]+","+bn[7]+"]");
										}
//										for (int np = 0; np < other_matched.length; np++) if (planes[nsTile][np] != null){
										for (int np = 0; np < other_matched.length; np++) if ((planes[nsTile][np] != null) && (planes[nsTile][np].getWeight() > 0.0)){ // disregard 0-weight planes
											int [] bn = planes[nsTile][np].getNeibBest();
											System.out.println("nsTile="+nsTile+":"+np+" best neighbors : ["+
											bn[0]+","+bn[1]+","+bn[2]+","+bn[3]+","+bn[4]+","+bn[5]+","+bn[6]+","+bn[7]+"]");
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
				for (int j = 0; j < planes[i].length; j++) if (planes[i][j] != null){
//					if (planes[i][j].getWeight() >= minWeight){ /// Later does not compare too minWeight
						double eigVal = planes[i][j].getValue(); // ************** generated does not have value ????????????
						double disp = planes[i][j].getSinglePlaneDisparity(false)[centerIndex];
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
//					}
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
//					for (int np = 0; np < planes[nsTile].length; np++) if (planes[nsTile][np] != null){
					for (int np = 0; np < planes[nsTile].length; np++) if ((planes[nsTile][np] != null) && (planes[nsTile][np].getWeight() > 0.0)){ // disregard 0-weight planes
					
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
								plane = planes[nsTile][np].getSinglePlaneDisparity(use_NaN);
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
			int  nlayer, // over multi-layer - do not render more than nlayer on top of each other
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
						int num_rendered = 0;
						for (int np = 0; np < planes[nsTile].length; np++){
							if ((shells[nsTile][np] == (nShell + 1)) && ((nlayer == 0) || (num_rendered < nlayer)) ) {
								num_rendered ++;
								plane = planes[nsTile][np].getSinglePlaneDisparity(use_NaN);
								int [] neib_best = planes[nsTile][np].getNeibBest();
								if (fuse && (neib_best != null)){
									if (debug) {
										dbg_img3 = new double [dbg_titles.length][];
										dbg_img = new double [dbg_titles.length][];
										dbg_this = new double [dbg_titles.length][];
										dbg_other = new double [dbg_titles.length][];
										dbg_img3[0] = planes[nsTile][np].getTriplePlaneDisparity();
										dbg_img[0] = planes[nsTile][np].getSinglePlaneDisparity(false);
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
									double [] plane_other = planes[nsTile][np].getSinglePlaneDisparity(false);
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

	
	public double [][][] getOverlappingShells(
			int  nlayer, // over multi-layer - do not render more than nlayer on top of each other
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
						int num_rendered = 0;
						for (int np = 0; np < planes[nsTile].length; np++){
							if ((shells[nsTile][np] == (nShell + 1)) && ((nlayer == 0) || (num_rendered < nlayer)) ) {
								num_rendered ++;
								plane = planes[nsTile][np].getSinglePlaneDisparity(use_NaN);
								int [] neib_best = planes[nsTile][np].getNeibBest();
								if (fuse && (neib_best != null)){
									if (debug) {
										dbg_img3 = new double [dbg_titles.length][];
										dbg_img = new double [dbg_titles.length][];
										dbg_this = new double [dbg_titles.length][];
										dbg_other = new double [dbg_titles.length][];
										dbg_img3[0] = planes[nsTile][np].getTriplePlaneDisparity();
										dbg_img[0] = planes[nsTile][np].getSinglePlaneDisparity(false);
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
									double [] plane_other = planes[nsTile][np].getSinglePlaneDisparity(false);
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
		return  null; // data;
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
			final double      maxValue, // do not combine with too bad planes with primary eigenvalue above this value ( 0 any OK)
			final int         num_passes,
			final boolean     stopBadNeib, // do not update supertile if any of connected neighbors is not good (false: just skip that neighbor)
			final double      normPow,
			
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
					maxValue,        // final double maxValue, // do not combine with too bad planes
					stopBadNeib,     // do not update supertile if any of connected neighbors is not good (false: just skip that neighbor)
					normPow,         // 0.0 - 8 neighbors pull 8 times more than 1. 1.0 - same as one
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
			final double maxValue, // do not combine with too bad planes
			final boolean stopBadNeib, // do not update supertile if any of connected neighbors is not good (false: just skip that neighbor)
			final double normPow,
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
//		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
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
								if (this_new_plane == null){
									System.out.println("Bug?  new_planes["+nsTile0+"]["+np0+"] == null");
									continue;
								}
								if (dl > 0) dbg_img[ 0] = this_new_plane.getSinglePlaneDisparity(false);
								if (dl > 0) dbg_img[ 1] = measured_planes[nsTile0][np0].getSinglePlaneDisparity(false);

								int [] neibs = this_new_plane.getNeibBest();
								this_new_plane.setWeight(0.0); //
								double num_merged = 0.0; // double to add fractional pull weight of the center
								for (int dir = 0; dir < neibs.length; dir++){
									if (neibs[dir] >= 0) {
										int stx = stx0 + dirsYX[dir][1];
										int sty = sty0 + dirsYX[dir][0];
										int nsTile = sty * stilesX + stx; // from where to get
										TilePlanes.PlaneData other_plane = this_new_plane.getPlaneToThis(
												mod_planes[nsTile][neibs[dir]], // neighbor, previous value
												dl - 1); // debugLevel);
										if (dl > 0) dbg_img[ 2 + dir] = other_plane.getSinglePlaneDisparity(false);
										if ((other_plane != null) && ((other_plane.getValue() <= maxValue) || (maxValue == 0))) {
											if (this_new_plane.getWeight() > 0.0){
												this_new_plane = this_new_plane.mergePlaneToThis(
														other_plane, // PlaneData otherPd,
														1.0,         // double    scale_other,
														false,       // boolean   ignore_weights,
														true,        // boolean   sum_weights,
														preferDisparity,
														dl - 1); // int       debugLevel)
												if (this_new_plane != null){
													num_merged += 1.0;
												}
												
											} else {
												this_new_plane = other_plane; 
											}
											new_planes[nsTile0][np0] = this_new_plane;
											if (dl > 0) dbg_img[10 + dir] = this_new_plane.getSinglePlaneDisparity(false);
										} else if (stopBadNeib){
											this_new_plane = null;
											break;
										}
									}
								}
								if (this_new_plane != null) {
									// average weight over participating directions, so the relative pull
									// does not depend on number of neighbors
									if ((num_merged > 0.0) && (this_new_plane != null) && (normPow > 1.0)){
										double scale = Math.pow(num_merged, normPow);
										this_new_plane.scaleWeight(1.0/scale);
										num_merged /=scale;
									}
									
									if (    (meas_pull > 0.0) &&
											(measured_planes != null) &&
											(measured_planes[nsTile0] != null) &&
											(measured_planes[nsTile0][np0] != null) &&
//											((measured_planes[nsTile0][np0].getValue() < maxValue) || (maxValue == 0))
											((measured_planes[nsTile0][np0].getValue() < maxValue) || (maxValue == 0) || 
													(this_new_plane == null) || (this_new_plane.getWeight() == 0.0)) // keep measured if impossible to merge
											){ // merge with "measured"

										if (this_new_plane.getWeight() > 0.0){
											this_new_plane = this_new_plane.mergePlaneToThis(
													measured_planes[nsTile0][np0], // PlaneData otherPd,
													meas_pull,                     // double    scale_other,
													false,                         // boolean   ignore_weights,
													true, // boolean   sum_weights,
													preferDisparity,
													dl - 1);                           // int       debugLevel)
											if (this_new_plane != null){
												num_merged += meas_pull;	// num_merged was 1.0 and weight is averaged over all neighbors 
											}
										} else {
											this_new_plane = measured_planes[nsTile0][np0].clone();
											num_merged = 1.0;	
										}
										new_planes[nsTile0][np0] = this_new_plane;
										if (dl > 0) dbg_img[18] =  this_new_plane.getSinglePlaneDisparity(false);
									}
									if ((num_merged > 0.0) && (this_new_plane != null)){
										this_new_plane.scaleWeight(1.0/num_merged);
									}
									// Revert if the result value is higher than imposed maximum
									if ((this_new_plane.getValue() > maxValue)  && (maxValue != 0)){
										this_new_plane = mod_planes[nsTile0][np0].clone();
										new_planes[nsTile0][np0] = this_new_plane;
									}

									// calculate largest disparity difference between old and new plane
									if (rslt_diffs != null){ // filter out outliers here?
										// get plane for both old and new, calc rms of diff
										double [] oldPlane = mod_planes[nsTile0][np0].getSinglePlaneDisparity(
												false); // use_NaN)
										double [] newPlane = new_planes[nsTile0][np0].getSinglePlaneDisparity(
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
											if (debugLevel > 0) {
												System.out.println("planesSmoothStep() nsTile0="+nsTile0+" rms = "+s);
											}
										}
										if (debugLevel > -1) {
											//										if ((s > 5.0) || (dl > 0)){
											if ((s > 5.0) || (dl > 0)){
												System.out.println("planesSmoothStep() nsTile0="+nsTile0+
														" num_merged="+num_merged + " rms = "+s+
														" new_weight = "+new_planes[nsTile0][np0].getWeight()+
														" old_weight = "+mod_planes[nsTile0][np0].getWeight()+
														" new_value = "+new_planes[nsTile0][np0].getValue()+
														" old_value = "+mod_planes[nsTile0][np0].getValue()+
														" measured_weight = "+measured_planes[nsTile0][np0].getWeight()+
														" measured_value = "+measured_planes[nsTile0][np0].getValue());
											}
										}

									}
									if ((dl > 0) && (debugLevel > 0)){
										showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
										sdfa_instance.showArrays(dbg_img, superTileSize, superTileSize, true, "smooth_step_x"+stx0+"_y"+sty0, titles);
									}
								} else { // if (this_new_plane != null)
									this_new_plane = mod_planes[nsTile0][np0].clone();
									new_planes[nsTile0][np0] = this_new_plane;
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
	 * Prepare visualization of the plane separation lines 
	 * @param split_planes per supertile, per plane sets of 3 plane data instances
	 *  (first, second, merged)
	 * @param debugLevel debug level
	 * @param dbg_X supertile horizontal index to show debug information
	 * @param dbg_Y supertile vertical index to show debug information
	 * @return array to be visualized as width=superTileSize*stilesX,
	 *  height = superTileSize*stilesY data. Lines value match plane index,
	 *  length corresponds to triple-sized supertiles
	 */
	
	public double [] showSplitLines(
			final TilePlanes.PlaneData[][][] split_planes,
			final int                        debugLevel,
			final int                        dbg_X,
			final int                        dbg_Y)
	{
		final int [][] dirsYX = {{-1, 0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1}};
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int width = superTileSize * stilesX;
		final int height = superTileSize * stilesY;
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		double [] split_lines = new double[width * height];
		final int  ss3 = 3 * superTileSize;
		for (int nsTile = 0; nsTile < split_planes.length; nsTile++) {
			if (split_planes[nsTile] != null) {
				int sty = nsTile / stilesX;  
				int stx = nsTile % stilesX;
				int dl = ((debugLevel > -1) && (nsTile == debug_stile)) ? 1:0;
				int min_np = (split_planes[nsTile].length > 1)? 1 : 0;
				for (int np = min_np; np < split_planes[nsTile].length; np++){
					if (split_planes[nsTile][np] != null){
						TilePlanes.PlaneData split_plane1 = split_planes[nsTile][np][0];
						TilePlanes.PlaneData split_plane2 = split_planes[nsTile][np][1];
						if ((split_plane1 != null) && ( split_plane2 != null)){
							double [] tpd =  split_plane1.getTriplePlaneDisparity();
							double [] tpd1 = split_plane2.getTriplePlaneDisparity();
							// difference of the 2 split planes disparities, show where
							// they have lowest value in any of the 2 orthogonal directions 
							for (int i = 0; i < tpd.length; i++) {
								tpd[i] = Math.abs(tpd[i] - tpd1[i]);
							}

							int [] vby = new int [ss3];
							for (int x = 0; x < ss3; x ++){
								vby[x] = 0;
								for (int y = 1; y < ss3; y ++){
									if (tpd[y * ss3 + x] < tpd[vby[x] * ss3 + x]){
										vby[x] = y;
									}
								}
							}
							int [] vbx = new int [ss3];
							for (int y = 0; y < ss3; y ++){
								vbx[y] = 0;
								for (int x = 1; x < ss3; x ++){
									if (tpd[y * ss3 + x] < tpd[y * ss3 + vbx[y]]){
										vbx[y] = x;
									}
								}
							}

							for (int dx = 0; dx < ss3; dx ++){
								if ((vby[dx] > 0) && (vby[dx] < (ss3-1))) {
									int y = vby[dx] + superTileSize * (sty - 1);
									int x = dx + superTileSize * (stx - 1);
									if ((y >= 0) && (y < height) && (x >= 0) && (x < width)){
										split_lines[ y* width + x]	= np - min_np + 1;		
									}
								}
							}

							for (int dy = 0; dy < ss3; dy ++){
								if ((vbx[dy] > 0) && (vbx[dy] < (ss3-1))) {
									int x = vbx[dy] + superTileSize * (stx - 1);
									int y = dy + superTileSize * (sty - 1);
									if ((y >= 0) && (y < height) && (x >= 0) && (x < width)){
										split_lines[ y* width + x]	= np - min_np + 1;		
									}
								}
							}
						}
					}
				}
			}
		}		      
		return split_lines;
	}
	
	/**
	 * Re-generate "planes" (ellipsoids) from the measured data according to provided
	 * plane pairs for some supertiles.
	 * Current version destroys connections between neighbor planes, they need to be
	 * re-calculated
	 * TODO: consider amplifying difference (angle) between planes before/during/after
	 * running this method. Like 3d "sharpening" to compensate for being washed out.
	 * @param planes per-supertile (in linescan order), per-plane data to be replaced
	 * @param brokenPd per-supertile , per-plane 3 PD instances (two first make the pair).
	 *        generated by breakPlanesToPairs method. This data will not be modified.
	 * @param max_diff maximal normalized disparity difference from the plane to consider
	 * @param other_diff maximal difference of the added tile ratio to the average
	 *  disparity difference of the exclusively selected tiles
	 * @param non_exclusive allow tiles to belong to both planes of the pair
	 * @param use_other_planes allow other tiles from the same supertile
	 * @param measSel (with use_other_planes) select measurements for supertiles:
	 *  +1 - combo, +2 - quad +4 - hor +8 - vert
	 * @param allow_parallel allow parallel shift of the specified planes before adding
	 *  more tiles (from the same pair and/or other tiles of the supertile (see
	 *  non_exclusive, use_other_planes
	 * @param splitXY generate XY separation masks for split planes, so tiles will not
	 *        be allowed to snap to the other plane in pair
	 * @param splitXYTolerance XY separation disparity tolerance (there will be some
	 *  overlapping if > 0.0
	 * @param disp_far disparity low limit 
	 * @param disp_near disparity high limit (Double.NaN - no limit)
	 * @param dispNorm normalize disparity difference (keep relative) for disparities above
	 * @param min_weight minimal total weight of the plane
	 * @param min_tiles minimal number of tiles in a plane
	 * @param targetEigen target primary eigenvalue when removing outliers
	 * @param fractOutliers maximal fraction of outliers from all tiles to remove
	 * @param maxOutliers maximal number of outliers to remove
	 * @param strength_floor subtract from correlation strength (and limit by 0) to use
	 *  correlation strength as weight
	 * @param strength_pow raise correlation strength (after subtracting strength_floor) to
	 *  this power before using as weight
	 * @param smplMode use square sample mode, false - single-tile samples 
	 * @param smplSide size of the square sample side 
	 * @param smplNum number of averaged samples (should be <= smplSide * smplSide and > 1) 
	 * @param smplRms maximal square root of variance (in disparity pixels) to accept the result
	 * @param debugLevel debug level
	 * @param dbg_X supertile horizontal index to show debug information
	 * @param dbg_Y supertile vertical index to show debug information
	 * @return number of planes replaced
	 */
	// TODO: Mark replaced planes and reduce eigenvalue requirements for them (in case they are small)
	// TODO:  planesSmooth() add power function to modify pull/neighbors relation (8 pull stronger than 1, but not 8 times stronger
	public int replaceBrokenPlanes(
			final TilePlanes.PlaneData[][]   planes,
			final TilePlanes.PlaneData[][][] brokenPd,
			final double                     max_diff, // maximal disparity difference (0 - any), will be normalized by dispNorm
			final double                     other_diff, 
			final boolean                    non_exclusive,
			final boolean                    use_other_planes, // TODO:
			final int                        measSel, // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			final boolean                    allow_parallel,
			
			final boolean                    splitXY,
			final double                     splitXYTolerance,
			
			// parameters to generate ellipsoids			
			final double                     disp_far, // minimal disparity to select (or NaN)
			final double                     disp_near, // maximal disparity to select (or NaN)
			final double                     dispNorm,   //  Normalize disparities to the average if above
			final double                     min_weight,
			final int                        min_tiles,
			// parameters to reduce outliers			
			final double                     targetEigen,   //     =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
			final double                     fractOutliers, //     =   0.3;  // Maximal fraction of outliers to remove
			final int                        maxOutliers,   //     =   20;  // Maximal number of outliers to remove
			final double                     strength_floor,
			final double                     strength_pow,
			
			final boolean                    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
			final int                        smplSide, //        = 2;      // Sample size (side of a square)
			final int                        smplNum, //         = 3;      // Number after removing worst
			final double                     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
			
			final int                        debugLevel,
			final int                        dbg_X,
			final int                        dbg_Y)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger ai_numThread = new AtomicInteger(0);
		final int [] replaced = new int[threads.length];
		final int [][] dbg_dirsYX = {{-1, 0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1}};
		
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					double [][] dbg_img = null;
					String [] dbg_titles = null;
					TilePlanes.PlaneData[][]   dbg_planes = planes;
					int numThread = ai_numThread.getAndIncrement(); // unique number of thread to write to rslt_diffs[numThread]
					for (int nsTile = ai.getAndIncrement(); nsTile < planes.length; nsTile = ai.getAndIncrement()) {
						int sty = nsTile / stilesX;  
						int stx = nsTile % stilesX;
						int dl = ((debugLevel > -1) && (nsTile == debug_stile)) ? 4:0;
						int num_bplanes = 2; // brokenPd[nsTile].length;
						int num_meas_layers = measuredLayers.getNumLayers();
						int [] dbg_indices = {
								0,
								num_bplanes,
								num_bplanes + 3 * num_meas_layers,
								num_bplanes + 3 * num_meas_layers + 2 * (num_bplanes * num_meas_layers),
								num_bplanes + 3 * num_meas_layers + 2 * (num_bplanes * num_meas_layers) + num_bplanes,
						};

						if ( brokenPd[nsTile] != null) {
							if (dl > 0){
								System.out.println("replaceBrokenPlanes nsTile="+nsTile);
								dbg_titles = new String[dbg_indices[4]];
								for (int i = 0; i < num_bplanes; i++){
									dbg_titles[i + dbg_indices[0]]="src_plane_"+i;
								}
								for (int j = 0; j < num_meas_layers; j++){
									dbg_titles[3* j + dbg_indices[1] + 0] =     "disp_"+j;
									dbg_titles[3* j + dbg_indices[1] + 1] =     "mask_disp_"+j;
									dbg_titles[3* j + dbg_indices[1] + 2] =     "str_"+j;
								}
								for (int i = 0; i < num_bplanes; i++){
									for (int j = 0; j < num_meas_layers; j++){
										dbg_titles[2 * (num_meas_layers * i +j) + dbg_indices[2] + 0] =  "t"+i+"_m"+j;
										dbg_titles[2 * (num_meas_layers * i +j) + dbg_indices[2] + 1] =  "ft"+i+"_m"+j;
									}
								}
								for (int i = 0; i < num_bplanes; i++){
									dbg_titles[i + dbg_indices[3]]="plane_"+i;
								}
							}
							
							TilePlanes.PlaneData[][] bpd = new TilePlanes.PlaneData[brokenPd[nsTile].length][];
							for (int np = 0; np < bpd.length; np++){
								if (brokenPd[nsTile][np] != null){
									bpd[np] = new TilePlanes.PlaneData[2];
									bpd[np][0] = brokenPd[nsTile][np][0].clone();
									bpd[np][1] = brokenPd[nsTile][np][1].clone(); // no need to clone [2]
								}
							}
							if (splitXY) {
								for (int np = 0; np < bpd.length; np++){
									if (brokenPd[nsTile][np] != null){
										// if it will fail (return false) no selection masks, still possible to separate
										planes[nsTile][np].calcSelMasks(
												bpd[np][0],
												bpd[np][1],
												splitXYTolerance, // double    tolerance,
												min_tiles, 
												dl); // int       debugLevel
									}
								}
							}
							
							// now works on the pairs in bpd, null them out in case of failure
							for (int np = 0; np < bpd.length; np++){
								if ((dl > 1) && (bpd[np] !=null)) {
									int ss2 = 2 * superTileSize;
									dbg_img = new double [dbg_titles.length][];
									for (int ni = 0; ni < num_bplanes ; ni++) {
										dbg_img[ni + dbg_indices[0]] =  bpd[np][ni].getDoublePlaneDisparity(false);
										boolean [] sm = bpd[np][ni].getSelMask();
										// checkerboard pattern over disabled part of the plane
										if (sm != null) {
											for (int i = 0; i < sm.length; i++){
												if (!sm[i] && ((((i % ss2)+ (i / ss2)) & 1) != 0)){
													dbg_img[ni + dbg_indices[0]][i] = Double.NaN;
												}
											}
										}
									}
									// add connections to neighbors
									for (int dir = 0; dir <8; dir++){
										int indx = (ss2 / 2) * (ss2 + 1);
										int dindx = dbg_dirsYX[dir][0] * ss2 + dbg_dirsYX[dir][1];
										for (int l = 0; l < (ss2 / 2 - 2); l++){
											if ( l > 2){ // keep center not modified
												for (int ni = 0; ni < num_bplanes ; ni++) {
													if (bpd[np][ni].getNeibBest(dir) >= 0){
														dbg_img[ni + dbg_indices[0]][indx] = Double.NaN;
													}
												}
											}
											indx += dindx; 
										}
									}
									// show disparity (all and masked for this plane) and strength
									for (int ml = 0; ml < num_meas_layers; ml++) {
										double [][] dbg_disp_str = planes[nsTile][np].getMeasuredLayers().getDisparityStrength(
												ml,                     // int num_layer,
												stx,        // int stX,
												sty,        // int stY,
												null,       // boolean [] sel_in,
												strength_floor,         //  double strength_floor,
												strength_pow,  // double strength_pow,
												true);                  // boolean null_if_none);
										if (dbg_disp_str != null){
											boolean [] dbg_msel = planes[nsTile][np].getMeasSelection(ml);
											dbg_img[3 * ml + dbg_indices[1] + 0] = dbg_disp_str[0];
											dbg_img[3 * ml + dbg_indices[1] + 2] = dbg_disp_str[1];
											dbg_img[3 * ml + dbg_indices[1] + 1] = dbg_disp_str[0].clone();
											if (dbg_msel != null) {
												for (int i = 0; i < dbg_msel.length; i++){
													if (!dbg_msel[i]) dbg_img[3 * ml + dbg_indices[1] + 1][i] = Double.NaN;
												}
											}
										}
									}
									
								} // (dl > 1)
								
								if (brokenPd[nsTile][np] != null){
									// Is it a single set to replace current plane? If yes, use all tiles, not just original selection
									int np_min = LOWEST_PLANE(planes[nsTile].length);
									boolean single_plane = (np_min == (planes[nsTile].length - 1));
									
									Boolean OK = planes[nsTile][np].splitPlaneTiles ( // now always OK TODO: add layer mask or make use_other_planes int
											bpd[np],          // PlaneData [] pd_set,
											single_plane,     // boolean      single_plane,
											max_diff,         // double       max_diff, // maximal disparity difference (0 - any), will be normalized by dispNorm
											other_diff,
											non_exclusive,    // boolean      non_exclusive,
											use_other_planes, // boolean      use_other_planes,
											smplMode,
											smplSide,
											smplNum,
											smplRms,
											measSel,          // int          measSel, // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
											allow_parallel,   //boolean      allow_parallel,
											dl); // int          debugLevel)
									if (!OK) {
										brokenPd[nsTile][np] = null;
										continue;
									}
									if (dl > 1) {// show initial selection (made of disparity clone - dbg_img[3 * ml + dbg_indices[1] + 0] = dbg_disp_str[0];
										for (int ni = 0; ni < num_bplanes ; ni++) {
											boolean [][] msel = new boolean [num_meas_layers][];
											for (int ml = 0; ml < num_meas_layers; ml++) {
												msel[ml] = bpd[np][ni].getMeasSelection(ml);
												if (msel[ml] != null) {
													dbg_img[2 * (num_meas_layers * ni +ml) + dbg_indices[2] + 0] =  dbg_img[3 * ml + dbg_indices[1] + 0].clone();
													for (int i = 0; i < msel[ml].length; i++){
														if (!msel[ml][i]) dbg_img[2 * (num_meas_layers * ni +ml) + dbg_indices[2] + 0][i] = Double.NaN;
													}
												}
											}
										}
									}
									
									for (int npip = 0; npip < bpd[np].length; npip++) {
										double [][][] disp_strength = bpd[np][npip].getPlaneFromMeas(
												null,           // boolean [][] tile_sel, // null - do not use, {} use all (will be modified)
												null,
												disp_far,       // double       disp_far, // minimal disparity to select (or NaN)
												disp_near,      // double       disp_near, // maximal disparity to select (or NaN)
												dispNorm,       // double       dispNorm,   //  Normalize disparities to the average if above
												min_weight,     // double       min_weight,
												min_tiles,      // int          min_tiles,
												strength_floor, // double       strength_floor,
												strength_pow,   // double       strength_pow,
// OK?												
												smplMode,
												smplSide,
												smplNum,
												smplRms,
												dl);            // int          debugLevel)
										if (disp_strength == null) break;
										// remove outliers //removeOutliers
										// now try to remove outliers
										int max_outliers = (int) Math.round(bpd[np][npip].getNumPoints() * fractOutliers);
										if (max_outliers > maxOutliers) max_outliers = maxOutliers;
										double targetV = corrMaxEigen(
												targetEigen,
												dispNorm,
												bpd[np][npip]);
										if (bpd[np][npip].getValue() > targetV) {
											OK = bpd[np][npip].removeOutliers( // getPlaneFromMeas should already have run
													disp_strength,
													targetV,      // double     targetEigen, // target eigenvalue for primary axis (is disparity-dependent, so is non-constant)
													max_outliers, // int        maxRemoved,  // maximal number of tiles to remove (not a constant)
													dl); // int        debugLevel)
											if (!OK) break;
										}										
									}
									if (!OK) {
										brokenPd[nsTile][np] = null;
										continue;
									}
									
									if (dl > 1) {// filtered selection (made of disparity clone - dbg_img[3 * ml + dbg_indices[1] + 0] = dbg_disp_str[0];
										for (int ni = 0; ni < num_bplanes ; ni++) {
											boolean [][] msel = new boolean [num_meas_layers][];
											for (int ml = 0; ml < num_meas_layers; ml++) {
												msel[ml] = bpd[np][ni].getMeasSelection(ml);
												if (msel[ml] != null) {
													dbg_img[2 * (num_meas_layers * ni +ml) + dbg_indices[2] +1] =  dbg_img[3 * ml + dbg_indices[1] + 0].clone();
													for (int i = 0; i < msel[ml].length; i++){
														if (!msel[ml][i]) dbg_img[2 * (num_meas_layers * ni +ml) + dbg_indices[2] + 1][i] = Double.NaN;
													}
												}
											}
										}
										// show new planes
										for (int ni = 0; ni < num_bplanes ; ni++) {
											dbg_img[ni + dbg_indices[3]] =  bpd[np][ni].getDoublePlaneDisparity(false);
										}
										// actually show the image
										
										if (debugLevel > -1){
											showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
											sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2 * superTileSize, true, "replaceBrokenPlanes"+stx+"_y"+sty+"_p"+np, dbg_titles);
										}
									}
								}
							}
							// now see if anything is left
							int npairs = 0;
							for (int np = 0; np < bpd.length; np++){
								if (bpd[np] != null) npairs++;
							}
							
							// Fixing that first plane in multi=-plane is not processed
							int np_start = LOWEST_PLANE(2); // planes[nsTile].length);

							if (npairs > 0){
								TilePlanes.PlaneData[] old_planes = planes[nsTile];
								int old_len = old_planes.length - LOWEST_PLANE(old_planes.length); 
								planes[nsTile] = new TilePlanes.PlaneData[old_len + npairs + np_start];
								int npr = 0;
								for (int np = LOWEST_PLANE(old_planes.length); np < np_start;  np++) {
									planes[nsTile][npr++] = null; // see if there will be any conflicts, fix or replace with old_planes[0];  
								}
								for (int np = 0; np < bpd.length; np++){
									if (bpd[np] != null) {
										planes[nsTile][npr++] = bpd[np][0]; 
										planes[nsTile][npr++] = bpd[np][1]; 
									} else {
										planes[nsTile][npr++] = old_planes[np]; 
									}
								}
							}
							replaced[numThread] += npairs;
						}								
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		int replaced_planes = 0;
		for (int i = 0; i < replaced.length; i++) replaced_planes += replaced[i];
		return replaced_planes;
	}
	
	/**
	 * Create candidate planes to break a single plane in 2 by splitting consecutive connected
	 * neighbors in 2 groups that make the smallest weighted sum of the eigenvalues
	 * Groups may be optionally supplemented by the center supertile
	 * Result weights are sums of participating supertile weights, for normalization  (if needed)
	 * divide by number of neighbors plus center_pull 
	 * @param center_planes [per supertile][per plane] array of plane objects to use as
	 *        the source of connections and optionally other data for the center tiles
	 * @param neib_planes  [per supertile][per plane] array of plane objects to use for
	 *        neighbor data. May be the same as center_planes or different
	 * @param center_pull - merge with center plane when calculating half-planes with this
	 *  relative weight: 0.0 - only neighbors, 1.0 - same weight of the center as each neighbor.
	 * @param min_neibs - minimal number of connected neighbors to work with (>=2)       
	 * @param splitMinWeight -  minimal weight of the tile to split       
	 * @param splitMinQuality - minimal quality of split (EV of the single plane over weighted average
	 *        of each half-plane        
	 * @param preferDisparity - the first eigenvalue/vector is the most disparity-like
	 *                          (false - smallest eigenvalue)
	 * @param debugLevel debug level
	 * @param dbg_X supertile horizontal index to show debug information
	 * @param dbg_Y supertile vertical index to show debug information
	 * @return a pair of plane objects for each [supertile][plane][3] and a combination of both
	 *         (to compare eigenvalues)
	 */
	public TilePlanes.PlaneData[][][] breakPlanesToPairs(
			final TilePlanes.PlaneData[][] center_planes,   // measured_planes,
			final TilePlanes.PlaneData[][] neib_planes,     //mod_planes,
			final double                   center_pull,
			final int                      min_neibs,       // 2
			final double                   splitMinWeight,  //     =  2.0;  // Minimal weight of split plains to show
			final double                   splitMinQuality, //    =  1.1;  // Minimal split quality to show
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
		final String [] titles= {"center", "first", "second", "all", "f-s"};

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
								if (center_plane.getWeight() >= splitMinWeight) {
									if (dl > 0) dbg_img[ 0] = center_plane. getTriplePlaneDisparity();
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
												plane_triads[nSplit][nis].scaleWeight(center_pull);
												for (int dir = 0; dir <8; dir++){
													int np = neibs12[nis][nSplit][dir];
													if (np >= 0){
														int stx = stx0 + dirsYX[dir][1];
														int sty = sty0 + dirsYX[dir][0];
														int nsTile = sty * stilesX + stx; // from where to get
														TilePlanes.PlaneData other_plane = plane_triads[nSplit][nis].getPlaneToThis(
																neib_planes[nsTile][np],
																dl-1); // debugLevel);
														if (other_plane != null){
															if (plane_triads[nSplit][nis].getWeight() > 0.0){
																plane_triads[nSplit][nis] = plane_triads[nSplit][nis].mergePlaneToThis(
																		other_plane, // PlaneData otherPd,
																		1.0,         // double    scale_other,
																		false,       // boolean   ignore_weights,
																		true, // boolean   sum_weights,
																		preferDisparity,
																		dl-1); // int       debugLevel)
															} else {
																plane_triads[nSplit][nis] = other_plane; 
															}
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
												if (dl >0){
													plane_triads[nSplit][2] = plane_triads[nSplit][0].clone().mergePlaneToThis(
															plane_triads[nSplit][1], // PlaneData otherPd,
															1.0,                    // double    scale_other,
															false,                  // boolean   ignore_weights,
															true, // boolean   sum_weights,
															preferDisparity,
															dl-1);										
												}
											} else {
												split_quality[nSplit] = Double.NaN;
											}
											if (dl >0){
												System.out.println("breakPlanesToPairs(): split_quality["+nSplit+"] = "+split_quality[nSplit]);
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
													1.0,                     // double    scale_other,
													false,                   // boolean   ignore_weights,
													true,                    // boolean   sum_weights,
													preferDisparity,
													dl-1);
											//  check quality
											if ((plane_triads[best_index][2] != null) &&
													(plane_triads[best_index][2].getValue()/split_quality[best_index] > splitMinQuality)) {
												for (int nis = 0; nis < 2; nis++) {
													plane_triads[best_index][nis].setNeibBest(neibs12[nis][best_index]);
												}
												rslt_planes[nsTile0][np0] = plane_triads[best_index];
												if (dl > 0) {
													int ss3 = 3 * superTileSize;
													dbg_img[ 1] = plane_triads[best_index][0].getTriplePlaneDisparity();
													dbg_img[ 2] = plane_triads[best_index][1].getTriplePlaneDisparity();
													dbg_img[ 3] = plane_triads[best_index][2].getTriplePlaneDisparity();
													dbg_img[ 4] = new double [dbg_img[ 1].length];
													for (int i = 0; i< dbg_img[ 1].length; i++){
														dbg_img[ 4][i] = dbg_img[ 1][i] - dbg_img[ 2][i];
													}
													int [] vby = new int [ss3];
													double [] dbg_split = dbg_img[ 4];
													for (int dx = 0; dx < ss3; dx ++){
														vby[dx] = 0;
														for (int dy = 1; dy < ss3; dy ++){
															int indx = dy * ss3 + dx;
															if (Math.abs(dbg_split[indx]) < Math.abs(dbg_split[vby[dx] * ss3 + dx])){
																vby[dx] = dy; // indx;
															}
														}
													}
													for (int dy = 0; dy < ss3; dy ++){
														int bx = 0;
														for (int dx = 1; dx < ss3; dx ++){
															if (Math.abs(dbg_split[dy * ss3 + dx]) < Math.abs(dbg_split[dy * ss3 + bx])){
																bx = dx;
															}
														}
														if ((bx > 0) && (bx < (ss3-1))) {
															dbg_split[dy * ss3 + bx] = Double.NaN;
															dbg_img[ 3][dy * ss3 + bx] = Double.NaN;
														}
													}

													for (int dx = 0; dx < ss3; dx ++){
														int by = vby[dx];
														if ((by > 0) && (by < (ss3-1))) {
															dbg_split[by * ss3 + dx] = Double.NaN; 
															dbg_img[ 3][by * ss3 + dx] = Double.NaN;
														}
													}											

													for (int dir = 0; dir <8; dir++){
														int indx = (ss3 / 2) * (ss3 + 1);
														int dindx = dirsYX[dir][0] * ss3 + dirsYX[dir][1];

														for (int l = 0; l < (ss3 / 2 - 2); l++){
															if ( l > 2){ // keep center not modified
																if (plane_triads[best_index][0].getNeibBest(dir) >= 0){
																	dbg_img[ 1][indx] = Double.NaN;
																}
																if (plane_triads[best_index][1].getNeibBest(dir) >= 0){
																	dbg_img[ 2][indx] = Double.NaN;
																}
															}
															indx += dindx; 
														}
													}
													System.out.println("breakPlanesToPairs(): num_splits="+num_splits+", best_index="+best_index);
													for (int nSplit = 0; nSplit < num_splits; nSplit++) {
														if ((plane_triads[nSplit][0] != null) && (plane_triads[nSplit][1] != null) && (plane_triads[nSplit][2] != null)){
															System.out.print("breakPlanesToPairs():"+nSplit+" ");
															for (int dir = 0; dir < 8; dir++){
																System.out.print((neibs12[0][nSplit][dir] >=0)?"+":"-");
															}
															System.out.print(" ");
															for (int dir = 0; dir < 8; dir++){
																System.out.print((neibs12[1][nSplit][dir] >=0)?"+":"-");
															}
															System.out.print(" ");
															for (int dir = 0; dir < 8; dir++){
																System.out.print((plane_triads[nSplit][2].getNeibBest(dir) >=0)?"+":"-");
															}
															System.out.println(" : "+split_quality[nSplit]+
																	" "+plane_triads[nSplit][0].getValue()+
																	" "+plane_triads[nSplit][1].getValue()+
																	" "+plane_triads[nSplit][2].getValue()+
																	" w"+plane_triads[nSplit][0].getWeight()+
																	" w"+plane_triads[nSplit][1].getWeight()+
																	" w"+plane_triads[nSplit][2].getWeight()
																	);
														}
													}

													// merged_eig_val are not set - will have to be recalculated when updating after changing tile selections
													// to make each plane
													// TODO save gain from splitting to planes 

													if (debugLevel > 3){
														showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
														sdfa_instance.showArrays(dbg_img, 3*superTileSize, 3*superTileSize, true, "breakPlanesToPairs"+stx0+"_y"+sty0+"_p"+np0, titles);
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
						if ((planes[nsTile][np] != null) && (planes[nsTile][np].getWeight() >= minWeight)){
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
			final int         debugLevel,
			final int         dbg_X,
			final int         dbg_Y)
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

	public double [][][] getDisparityStrengths(
			int        stMeasSel) //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert)
	{
		int numMeasLayers = measuredLayers.getNumLayers();
		double [][][] ds = new double[numMeasLayers][][];
		for (int ml = 0; ml <  numMeasLayers; ml++) if ((stMeasSel & ( 1 << ml)) != 0) {
			ds[ml] = 	measuredLayers.getDisparityStrength (
					ml, // int num_layer,
					this.strength_floor, // double strength_floor,
					this.strength_pow); //  double strength_pow)
		}
		return ds;
	}

	public boolean [][] getMeasurementSelections(
			int        stMeasSel) //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert)
	{
		int numMeasLayers = measuredLayers.getNumLayers();
		boolean [][] sels = new boolean[numMeasLayers][];
		for (int ml = 0; ml <  numMeasLayers; ml++) if ((stMeasSel & ( 1 << ml)) != 0) {
			sels[ml] = 	measuredLayers.getSelection (ml);
		}
		return sels;
	}
	
} // end of class SuperTiles
