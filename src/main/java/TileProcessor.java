/**
 **
 ** TileProcessor - Process tiles resulted from quad image sets
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  TileProcessor.java is free software: you can redistribute it and/or modify
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
import java.util.concurrent.atomic.AtomicInteger;

public class TileProcessor {
	public ArrayList <CLTPass3d> clt_3d_passes = null;
	public int       clt_3d_passes_size = 0; //clt_3d_passes size after initial processing
	private int      tilesX;
	private int      tilesY;
	private double   corr_magic_scale =    0.85;  // reported correlation offset vs. actual one (not yet understood)
	private double   trustedCorrelation =  4.0;   // trusted measured disparity difference (before scaling)    
	private int      tileSize =           8; // number of linear pixels in a tile (tile is  square tileSize*tileSize)   
	int               superTileSize =      8; // number of linear tiles in a super-tile (supertile is  square superTileSize*superTileSize tiles
	                                        // or (superTileSize*tileSize) * (superTileSize*tileSize) pixels, currently 64x64 pixels)
	public int       threadsMax =       100; // maximal number of frames to run  
	public int       globalDebugLevel = 0;
	
	// All parameters are set only once, during instantiation
	
	public TileProcessor(
			int tilesX,
			int tilesY,
			int tileSize,   
			int superTileSize,
			double scale,
			double trustedCorrelation,
			int threadsMax)
	{
		this.tilesX = tilesX;
		this.tilesY = tilesY;
		this.tileSize = tileSize;   
		this.superTileSize = superTileSize;
		this.corr_magic_scale = scale;
		this.trustedCorrelation = trustedCorrelation;
		this.threadsMax = threadsMax; 
	}
	public int getTilesX() {return tilesX;};
	public int getTilesY() {return tilesY;};
	public int getTileSize() {return tileSize;};
	public int getSuperTileSize() {return superTileSize;};
	public void setTrustedCorrelation(double trustedCorrelation)
	{
		this.trustedCorrelation = trustedCorrelation;
	}
	public double getTrustedCorrelation()
	{
		return this.trustedCorrelation;
	}
	public double getMagicScale()
	{
		return this.corr_magic_scale;
	}
	
//	public void setMagicScale (double scale)
//	{
//		this.corr_magic_scale = scale;
//	}

	public void resetCLTPasses(){
		clt_3d_passes = new ArrayList<CLTPass3d>();
		clt_3d_passes_size = 0;
	}
	
	public void saveCLTPasses(){
		clt_3d_passes_size = clt_3d_passes.size();
	}

	public void trimCLTPasses(){
		while (clt_3d_passes.size() > clt_3d_passes_size){
			clt_3d_passes.remove(clt_3d_passes_size);
		}
	}

	
	
	
	/**
	 * Basic combining: find smallest residual disparity and use the tile data from it
	 * Copy link to texture tile from the same pass, "forced" bit in tile_op is copied too
	 * Even when this method compares calculated values, it still only copies raw ones, all derivatives should
	 * be re-calculated for the new combined pass 
	 * 
	 * Calculates max_tried_disparity that shows maximal tried dieparity for each tile, regardless of the reulsts/strength
	 * @param passes list of passes to merge
	 * @param firstPass first index in the list to use
	 * @param lastPass last index in the list to use
	 * @param debugLevel debug level
	 * @return combined pass, contains same data as after the measuremnt of the actual one
	 */
	public CLTPass3d combinePasses(
			 final ArrayList <CLTPass3d> passes,
			 final int                   firstPass,
			 final int                   lastPassPlus1,
			 final boolean               skip_combo, // do not process other combo scans
			 final boolean               use_last,   // use last scan data if nothing better 
			 // TODO: when useCombo - pay attention to borders (disregard)
			 final boolean               usePoly,  // use polynomial method to find max), valid if useCombo == false
			 final double                minStrength, // ignore too weak tiles
			 final boolean               show_combined)
	{
		CLTPass3d combo_pass =new CLTPass3d(this);
		final int tlen = tilesX * tilesY;
		combo_pass.disparity =            new double [tilesY][tilesX];
		combo_pass.tile_op =              new int [tilesY][tilesX];
		combo_pass.disparity_map =        new double [ImageDtt.DISPARITY_TITLES.length][tlen];
		combo_pass.texture_tiles =        new double [tilesY][tilesX][][];
		combo_pass.max_tried_disparity =  new double [tilesY][tilesX];
		combo_pass.is_combo =             true;
		showDoubleFloatArrays sdfa_instance = null;
		String[] titles = null;
		int dbg_tile = -1; // 27669;
		double [][] dbg_data = null;
		if (show_combined) {
			sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
			int numScans = lastPassPlus1 - firstPass;
			titles = new String [3 * (numScans + 1) + 1];
			dbg_data = new double [titles.length][tlen];
			for (int i = 0; i < numScans; i++) {
				CLTPass3d dbg_pass = passes.get(firstPass + i);
				if (dbg_pass.disparity_map != null) {
					for (int ty = 0; ty < tilesY; ty++) for (int tx = 0; tx < tilesX; tx++){
						int nt = ty * tilesX +tx;
						dbg_data[i                     ][nt] = dbg_pass.disparity[ty][tx];
						dbg_data[i + 1 * (numScans + 1)][nt] = dbg_pass.disparity_map[ImageDtt.DISPARITY_INDEX_CM][nt];
						dbg_data[i + 2 * (numScans + 1)][nt] = dbg_pass.disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX][nt];
					}
					titles[i                     ] = "disparity_"+i;
					titles[i + 1 * (numScans + 1)] = "cm_disparity_"+i;
					titles[i + 2 * (numScans + 1)] = "strength_"+i;
				}
			}
		}

		
		for (int ty = 0; ty < tilesY; ty ++) for (int tx = 0; tx < tilesX; tx ++) combo_pass.texture_tiles[ty][tx] = null;
		for (int ty = 0; ty < tilesY; ty ++) {
			for (int tx = 0; tx < tilesX; tx ++){
				int nt = ty * tilesX + tx;
				int best_index = -1;
				int best_weak_index = -1;
				// for "strong" result (above minStrength) the best fit is with smallest residual disparity
				// for weak ones - the strongest.
				// TODO: handle ortho (if requested so)
				// after refine - use last, after extend - use strongest of all
				
				double adiff_best =      Double.NaN;
				double strongest_weak = 0.0;
				combo_pass.max_tried_disparity[ty][tx] = Double.NaN;
				for (int ipass = firstPass;  ipass <lastPassPlus1; ipass++ ){
					CLTPass3d pass = passes.get(ipass);
					if (nt == dbg_tile) {
						System.out.println("combinePasses(): ipass = "+ipass+" nt = "+nt+" pass.tile_op["+ty+"]["+tx+"]="+pass.tile_op[ty][tx]+
								" pass.isCombo()="+(pass.isCombo())+" pass.isProcessed()="+(pass.isProcessed()));
					}
					if (    pass.isMeasured() &&
							(pass.tile_op[ty][tx] != 0 )) {
						if (	(Double.isNaN(combo_pass.max_tried_disparity[ty][tx]) ||
								(pass.disparity[ty][tx] > combo_pass.max_tried_disparity[ty][tx]))){
							combo_pass.max_tried_disparity[ty][tx] = pass.disparity[ty][tx];
						}

						double adiff, strength;
						if (nt == dbg_tile){
//							System.out.println("combinePasses(): pass.calc_disparity["+nt+"]="+pass.calc_disparity[nt]+
//									" pass.disparity["+ty+"]["+tx+"] = "+pass.disparity[ty][tx]);
							System.out.println("combinePasses(): pass.disparity["+ty+"]["+tx+"] = "+pass.disparity[ty][tx]);
						}
						if (usePoly) { // just an amplitude of the polynomial maximum calculated disparity
							adiff =    Math.abs(pass.disparity_map[ImageDtt.DISPARITY_INDEX_POLY][nt]); // polynomial method
						} else { // just an amplitude of center of mass calculated disparity
							adiff = Math.abs(pass.disparity_map[ImageDtt.DISPARITY_INDEX_CM][nt]); // center mass method
						}
						strength = Math.abs(pass.disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX][nt]);
						if ((strength > 0.0) &&
								!Double.isNaN(adiff) &&
								((best_weak_index < 0) || (strongest_weak < strength))) {
							best_weak_index = ipass;
							strongest_weak = strength;
						}
						if ((strength > minStrength) &&
								!Double.isNaN(adiff) && 
								((best_index < 0) || (adiff < adiff_best))) {
							best_index = ipass;
							adiff_best = adiff;
						}
						if (nt == dbg_tile){
							System.out.println("combinePasses(): strength="+strength+" best_weak_index="+best_weak_index+
									" best_index="+best_index+" adiff_best="+adiff_best+" ipass="+ipass+"adiff="+adiff);
						}
					}
				}
				if (best_index < 0) {
					if (use_last && (passes.get(lastPassPlus1 - 1).tile_op[ty][tx] != 0)) { // last is prevferred and has some data
						best_index = lastPassPlus1 - 1;
					}	else { // just use result with strongest correlation
						best_index = best_weak_index;
					}
				}
				if (best_index >= 0){
					CLTPass3d pass = passes.get(best_index);
					combo_pass.tile_op[ty][tx] =          pass.tile_op[ty][tx];
					combo_pass.disparity[ty][tx] =        pass.disparity[ty][tx];
					if ((pass.texture_tiles == null) ||(combo_pass.texture_tiles == null)) {
						if ((ty==0) && (tx==0)) {
							System.out.println("BUG: best_index="+best_index);
						}
					} else {
						combo_pass.texture_tiles[ty][tx] =    pass.texture_tiles[ty][tx];
						for (int i = 0; i < ImageDtt.DISPARITY_TITLES.length; i++){
							combo_pass.disparity_map[i][nt] = pass.disparity_map[i][nt]; 
						}
					}

					// do not copy any of the calculated values - they should be re-calculated
				}
			}
		}
		if (show_combined) {
			int numScans = lastPassPlus1 - firstPass;
			for (int ty = 0; ty < tilesY; ty++) for (int tx = 0; tx < tilesX; tx++){
				int nt = ty * tilesX +tx;
				dbg_data[numScans                     ][nt] = combo_pass.disparity[ty][tx];
				dbg_data[numScans + 1 * (numScans + 1)][nt] = combo_pass.disparity_map[ImageDtt.DISPARITY_INDEX_CM][nt];
				dbg_data[numScans + 2 * (numScans + 1)][nt] = combo_pass.disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX][nt];
				dbg_data[3 * (numScans + 1)][nt] = combo_pass.max_tried_disparity[ty][tx];
			}
			titles[numScans                     ] = "disparity_combo";
			titles[numScans + 1 * (numScans + 1)] = "cm_disparity_combo";
			titles[numScans + 2 * (numScans + 1)] = "strength_combo";
			titles[3 * (numScans + 1)] =            "max_tried_disparity";
			sdfa_instance.showArrays(dbg_data,  tilesX, tilesY, true, "combo_scan_"+lastPassPlus1,titles);

		}
		return combo_pass;
	}

	/**
	 * Calculates calc_disparity, calc_disparity_hor,  calc_disparity_vert, strength, strength_hor, strength_vert,
	 * max_tried_disparity from the subset of a list of measurement passes (skipping non-measured)
	 * disparity, strength, *_hor and vert may come from the different scans.
	 * tile_op is set to 0/non-0 to show which tiles contain valid data
	 * disparity (target disparity) array is set to null as it can not be made valid (different values for disparity,
	 * disparity_hor and disparity_vert. 
	 * @param passes          list of scan passes to take data from
	 * @param firstPass       index of the first pass to use
	 * @param lastPassPlus1   index plus 1 of the last pass to use
	 * @param trustedCorrelation maximal absolute value of measured correlation (no scaling) to trust (may use global trustedCorrelation) 
	 * @param disp_far        lowest disparity value to consider (does not apply to max_tried_disparity)
	 * @param disp_near       highest disparity value to consider (does not apply to max_tried_disparity)
	 * @param minStrength     full correlation strength to consider data to be reliable 
	 * @param minStrengthHor  horizontal (for vertical features) correlation strength to consider data to be reliable
	 * @param minStrengthVert vertical (for horizontal features) correlation strength to consider data to be reliable
	 * @param use_last        use last scan data if nothing strong enough (false - use the strongest)
	 * @param usePoly         use polynomial method to find max for full correlation, false - use center of mass
	 * @param copyDebug       copy data tyhat is only needed for debug purposes
	 * @return new composite scan pass (not added to the list
	 */
	public CLTPass3d compositeScan(
			 final ArrayList <CLTPass3d> passes,
			 final int                   firstPass,
			 final int                   lastPassPlus1,
//			 final boolean               skip_combo, // do not process other combo scans
			 final double                trustedCorrelation,
			 final double                disp_far,   // limit results to the disparity range
			 final double                disp_near,
			 final double                minStrength, 
			 final double                minStrengthHor,
			 final double                minStrengthVert,
			 final boolean               use_last,   //  
			 // TODO: when useCombo - pay attention to borders (disregard)
			 final boolean               usePoly,  // use polynomial method to find max), valid if useCombo == false
			 final boolean               copyDebug,
			 final int                   debugLevel)
	{
		final int dbg_tile = (debugLevel > 0)? 35114: -1; // x = 122, y= 108; -1; // 27669;
		CLTPass3d combo_pass =new CLTPass3d(this);
		final int tlen = tilesX * tilesY;
		final int disparity_index = usePoly ? ImageDtt.DISPARITY_INDEX_POLY : ImageDtt.DISPARITY_INDEX_CM;
		combo_pass.tile_op =              new int [tilesY][tilesX]; // for just non-zero
		combo_pass.disparity_map =        new double [ImageDtt.DISPARITY_TITLES.length][];
		for (int i = 0; i< ImageDtt.QUAD; i++) combo_pass.disparity_map[ImageDtt.IMG_DIFF0_INDEX + i] = new double[tlen];
		if (copyDebug){
			combo_pass.disparity_map[ImageDtt.DISPARITY_INDEX_CM] =            new double[tlen];
			combo_pass.disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX] =      new double[tlen];
			combo_pass.disparity_map[ImageDtt.DISPARITY_INDEX_HOR] =           new double[tlen];
			combo_pass.disparity_map[ImageDtt.DISPARITY_INDEX_HOR_STRENGTH] =  new double[tlen];
			combo_pass.disparity_map[ImageDtt.DISPARITY_INDEX_VERT] =          new double[tlen];
			combo_pass.disparity_map[ImageDtt.DISPARITY_INDEX_VERT_STRENGTH] = new double[tlen];
		}
		
		
		// for now - will copy from the best full correlation measurement
		combo_pass.texture_tiles =        new double [tilesY][tilesX][][];
		combo_pass.max_tried_disparity =  new double [tilesY][tilesX];
		combo_pass.is_combo =             true;
		combo_pass.calc_disparity =       new double [tlen];
		combo_pass.calc_disparity_hor =   new double [tlen];
		combo_pass.calc_disparity_vert =  new double [tlen];
		combo_pass.strength =             new double [tlen];
		combo_pass.strength_hor =         new double [tlen];
		combo_pass.strength_vert =        new double [tlen];
		
//		 * Calculates calc_disparity, calc_disparity_hor,  calc_disparity_vert, strength, strength_hor, strength_vert,
		
		
		for (int ty = 0; ty < tilesY; ty ++) for (int tx = 0; tx < tilesX; tx ++) combo_pass.texture_tiles[ty][tx] = null;
		for (int ty = 0; ty < tilesY; ty ++) {
			for (int tx = 0; tx < tilesX; tx ++){
				int nt = ty * tilesX + tx;
				int best_index = -1;
				int best_index_hor = -1;
				int best_index_vert = -1;
				int best_weak_index = -1;
				int best_weak_index_hor = -1;
				int best_weak_index_vert = -1;
				// for "strong" result (above minStrength) the best fit is with smallest residual disparity
				// for weak ones - the strongest.
				// TODO: handle ortho (if requested so)
				// after refine - use last, after extend - use strongest of all
				
				double adiff_best =          Double.NaN;
				double adiff_best_hor =      Double.NaN;
				double adiff_best_vert =     Double.NaN;
				double strongest_weak =      0.0;
				double strongest_weak_hor =  0.0;
				double strongest_weak_vert = 0.0;
				
				combo_pass.max_tried_disparity[ty][tx] = Double.NaN;
				if (nt == dbg_tile) {
					System.out.println("compositeScan(): nt = "+nt);
				}
				for (int ipass = firstPass;  ipass <lastPassPlus1; ipass++ ){
					CLTPass3d pass = passes.get(ipass);
					if (nt == dbg_tile) {
						System.out.println("compositeScan(): ipass = "+ipass+" nt = "+nt+" pass.tile_op["+ty+"]["+tx+"]="+pass.tile_op[ty][tx]+
								" pass.isCombo()="+(pass.isCombo())+" pass.isProcessed()="+(pass.isProcessed()));
					}
					if ( pass.isMeasured() && (pass.tile_op[ty][tx] != 0 )) { // current tile has valid data
						if (	(Double.isNaN(combo_pass.max_tried_disparity[ty][tx]) ||
								(pass.disparity[ty][tx] > combo_pass.max_tried_disparity[ty][tx]))){
							combo_pass.max_tried_disparity[ty][tx] = pass.disparity[ty][tx];
						}
						boolean last = (ipass == (lastPassPlus1-1)) && use_last;
						double mdisp =         pass.disparity_map[disparity_index][nt];
						double mdisp_hor =     pass.disparity_map[ImageDtt.DISPARITY_INDEX_HOR][nt];
						double mdisp_vert =    pass.disparity_map[ImageDtt.DISPARITY_INDEX_VERT][nt];
						
						double strength =      pass.disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX][nt];
						double strength_hor =  pass.disparity_map[ImageDtt.DISPARITY_INDEX_HOR_STRENGTH][nt];
						double strength_vert = pass.disparity_map[ImageDtt.DISPARITY_INDEX_VERT_STRENGTH][nt];

						double adiff =      Math.abs(mdisp);
						double adiff_hor =  Math.abs(mdisp_hor);
						double adiff_vert = Math.abs(mdisp_vert);
						
						if (adiff <= trustedCorrelation){
							double disp = mdisp/corr_magic_scale +  pass.disparity[ty][tx];
							if ((disp >= disp_far) && (disp <= disp_near) && !Double.isNaN(adiff)){
								if (strength >= minStrength) { 
									if (!(adiff >= adiff_best)){ // adiff_best == Double.NaN works too 
										adiff_best = adiff;
										best_index = ipass;
									}
								} else {
									if ((last && (strength > 0.0)) || (strength > strongest_weak)){ 
										strongest_weak = strength;
										best_weak_index = ipass;
									}
								}
							}
						}

						if (adiff_hor <= trustedCorrelation){
							double disp_hor = mdisp_hor/corr_magic_scale +  pass.disparity[ty][tx];
							if ((disp_hor >= disp_far) && (disp_hor <= disp_near) && !Double.isNaN(adiff_hor)){
								if (strength_hor >= minStrength) { 
									if (!(adiff_hor >= adiff_best_hor)){ // adiff_best == Double.NaN works too 
										adiff_best_hor = adiff_hor;
										best_index_hor = ipass;
									}
								} else {
									if ((last && (strength_hor > 0.0))  || (strength_hor > strongest_weak_hor)){ 
										strongest_weak_hor = strength_hor;
										best_weak_index_hor = ipass;
									}
								}
							}
						}
						
						if (adiff_vert <= trustedCorrelation){
							double disp_vert = mdisp_vert/corr_magic_scale +  pass.disparity[ty][tx];
							if ((disp_vert >= disp_far) && (disp_vert <= disp_near) && !Double.isNaN(adiff_vert)){
								if (strength_vert >= minStrength) { 
									if (!(adiff_vert >= adiff_best_vert)){ // adiff_best == Double.NaN works too 
										adiff_best_vert = adiff_vert;
										best_index_vert = ipass;
									}
								} else {
									if ((last && (strength_vert > 0.0))  || (strength_vert > strongest_weak_vert)){ 
										strongest_weak_vert = strength_vert;
										best_weak_index_vert = ipass;
									}
								}
							}
						}
					}
				}
				if (best_index < 0)      best_index =      best_weak_index;						
				if (best_index_hor < 0)  best_index_hor =  best_weak_index_hor;						
				if (best_index_vert < 0) best_index_vert = best_weak_index_vert;

				if (best_index >= 0){
					CLTPass3d pass = passes.get(best_index);
					combo_pass.tile_op[ty][tx] =            pass.tile_op[ty][tx];
					if (pass.texture_tiles != null) {
						combo_pass.texture_tiles[ty][tx] =  pass.texture_tiles[ty][tx];
					}
					combo_pass.calc_disparity[nt] =         pass.disparity_map[disparity_index][nt]/corr_magic_scale +  pass.disparity[ty][tx];
					combo_pass.strength[nt] =               pass.disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX][nt];
					for (int i = 0; i< ImageDtt.QUAD; i++)  combo_pass.disparity_map[ImageDtt.IMG_DIFF0_INDEX + i][nt] = pass.disparity_map[ImageDtt.IMG_DIFF0_INDEX + i][nt];
					if (copyDebug){
						combo_pass.disparity_map[ImageDtt.DISPARITY_INDEX_CM][nt] =            pass.disparity_map[ImageDtt.DISPARITY_INDEX_CM][nt];
						combo_pass.disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX][nt] =      pass.disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX][nt];
					}
				} else {
					combo_pass.calc_disparity[nt] = Double.NaN;
					combo_pass.strength[nt] =       0.0; 
					if (copyDebug){
						combo_pass.disparity_map[ImageDtt.DISPARITY_INDEX_CM][nt] =            Double.NaN;
						combo_pass.disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX][nt] =      0.0;
					}
				}


				if (best_index_hor >= 0){
					CLTPass3d pass = passes.get(best_index_hor);
					combo_pass.tile_op[ty][tx] =            pass.tile_op[ty][tx]; // just non-zero
					combo_pass.calc_disparity_hor[nt] =     pass.disparity_map[ImageDtt.DISPARITY_INDEX_HOR][nt]/corr_magic_scale +  pass.disparity[ty][tx];
					combo_pass.strength_hor[nt] =           pass.disparity_map[ImageDtt.DISPARITY_INDEX_HOR_STRENGTH][nt];
					if (copyDebug){
						combo_pass.disparity_map[ImageDtt.DISPARITY_INDEX_HOR][nt] =           pass.disparity_map[ImageDtt.DISPARITY_INDEX_HOR][nt];
						combo_pass.disparity_map[ImageDtt.DISPARITY_INDEX_HOR_STRENGTH][nt] =  pass.disparity_map[ImageDtt.DISPARITY_INDEX_HOR_STRENGTH][nt];
					}
				} else {
					combo_pass.calc_disparity_hor[nt] = Double.NaN;
					combo_pass.strength_hor[nt] =       0.0; 
					if (copyDebug){
						combo_pass.disparity_map[ImageDtt.DISPARITY_INDEX_HOR][nt] =           Double.NaN;
						combo_pass.disparity_map[ImageDtt.DISPARITY_INDEX_HOR_STRENGTH][nt] =  0.0;
					}
					
				}

				if (best_index_vert >= 0){
					CLTPass3d pass = passes.get(best_index_vert);
					combo_pass.tile_op[ty][tx] =            pass.tile_op[ty][tx]; // just non-zero
					combo_pass.calc_disparity_vert[nt] =     pass.disparity_map[ImageDtt.DISPARITY_INDEX_VERT][nt]/corr_magic_scale +  pass.disparity[ty][tx];
					combo_pass.strength_vert[nt] =           pass.disparity_map[ImageDtt.DISPARITY_INDEX_VERT_STRENGTH][nt];
					if (copyDebug){
						combo_pass.disparity_map[ImageDtt.DISPARITY_INDEX_VERT][nt] =          pass.disparity_map[ImageDtt.DISPARITY_INDEX_VERT][nt];
						combo_pass.disparity_map[ImageDtt.DISPARITY_INDEX_VERT_STRENGTH][nt] = pass.disparity_map[ImageDtt.DISPARITY_INDEX_VERT_STRENGTH][nt];
					}
				} else {
					combo_pass.calc_disparity_vert[nt] = Double.NaN;
					combo_pass.strength_vert[nt] =       0.0; 
					if (copyDebug){
						combo_pass.disparity_map[ImageDtt.DISPARITY_INDEX_VERT][nt] =          Double.NaN;
						combo_pass.disparity_map[ImageDtt.DISPARITY_INDEX_VERT_STRENGTH][nt] = 0.0;
					}
				}
			}
		}
		combo_pass.fixNaNDisparity(); // mostly for debug, measured disparity should be already fixed from NaN
		return combo_pass;
	}
//trustedCorrelation	
	
/**
 * Calculate max_tried_disparity from the provided scans (use only actually measured ones)
 * @param passes list of passes
 * @param firstPass index of the first pass to look at
 * @param lastPassPlus1 last index plus 1 of the pass to look at (may include the current one too if it is measured - will be used)
 * @param new_scan current scan to attach result to
 */
	public void calcMaxTried(
			final ArrayList <CLTPass3d> passes,
			final int                   firstPass,
			final int                   lastPassPlus1, 
			final CLTPass3d             new_scan)
	{
		new_scan.max_tried_disparity =  new double [tilesY][tilesX];
		//		new_scan.is_combo =             true;
		int dbg_tile = -1; // 27669;

		for (int ty = 0; ty < tilesY; ty ++) {
			for (int tx = 0; tx < tilesX; tx ++){
				int nt = ty * tilesX + tx;
				if (nt == dbg_tile){
					System.out.println("calcMaxTried(): nt = "+nt+", tx = "+tx+", ty = "+ty);
				}
				new_scan.max_tried_disparity[ty][tx] = Double.NaN;
				for (int ipass = firstPass;  ipass <lastPassPlus1; ipass++ ){
					CLTPass3d pass = passes.get(ipass);
					if (    pass.isMeasured() &&
							(pass.tile_op[ty][tx] != 0 ) &&
							(Double.isNaN(new_scan.max_tried_disparity[ty][tx]) ||
									(pass.disparity[ty][tx] > new_scan.max_tried_disparity[ty][tx]))){
						new_scan.max_tried_disparity[ty][tx] = pass.disparity[ty][tx];
					}
				}
			}
		}
	}
	
	
	
	public int [] makeUnique(
			 final ArrayList <CLTPass3d> passes,
			 final int                   firstPass,
			 final int                   lastPassPlus1,
			 final CLTPass3d             new_scan,
			 final double                unique_tolerance,
			 final boolean               show_unique)
	{
		int [][] dbg_tile_op = null;
		if (show_unique){
			dbg_tile_op = new int [tilesY][];
			for (int ty = 0; ty < tilesY; ty ++){
				dbg_tile_op[ty] = new_scan.tile_op[ty].clone();
			}
		}
		int removed = 0, total = 0;
		for (int ty = 0; ty < tilesY; ty ++) {
			for (int tx = 0; tx < tilesX; tx ++){
				if (new_scan.tile_op[ty][tx] != 0){
					total ++;
					for (int ipass = firstPass;  ipass <lastPassPlus1; ipass++ ){
						CLTPass3d pass = passes.get(ipass);
						if ((pass.tile_op[ty][tx] != 0) && pass.isMeasured()){
							if (Math.abs(new_scan.disparity[ty][tx] - pass.disparity[ty][tx]) < unique_tolerance){
								new_scan.tile_op[ty][tx] = 0;
								removed++;
								break;
							}
						}
					}
				}
			}
		}
		if (show_unique){
			showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
			double [] dbg_data = new double [tilesY*tilesX];
			for (int ty = 0; ty < tilesY; ty ++) for (int tx = 0; tx < tilesX; tx ++){
				dbg_data[ty * tilesX + tx] = (dbg_tile_op[ty][tx] == 0) ? 0: ((new_scan.tile_op[ty][tx] == 0.0) ? 1.0 : 2.0 );
			}
			sdfa_instance.showArrays(dbg_data, tilesX, tilesY, "unique_scan_"+lastPassPlus1);
		}
		int [] rslt = {total-removed, removed};
		return rslt;
	}
	
	public int setupExtendDisparity( // returns number of new tiles to try
			final CLTPass3d   scan,            // combined scan with max_tried_disparity, will be modified to re-scan
			final CLTPass3d   last_scan,       // last prepared tile - can use last_scan.disparity, .border_tiles and .selected
			final CLTPass3d   bg_scan,         // background scan data
			final int         grow_sweep,      // 8; // Try these number of tiles around known ones 
			final double      grow_disp_max,   //  =   50.0; // Maximal disparity to try
			final double      tried_margin,    //  =  4.0; // consider alrdeady tried if within this margin from already tried 
			final double      grow_disp_step,  //  =   6.0; // Increase disparity (from maximal tried) if nothing found in that tile // TODO: handle enclosed dips?  
			final double      grow_min_diff,   // =    0.5; // Grow more only if at least one channel has higher variance from others for the tile
			EyesisCorrectionParameters.CLTParameters           clt_parameters,
			GeometryCorrection geometryCorrection,
			final boolean     show_debug,
			final int         threadsMax,  // maximal number of threads to launch                         
			final boolean     updateStatus,
			final int         debugLevel)
	{
		final int dbg_tile = 54627; // ty=159, tx = 249
		final int tlen = tilesY * tilesX;
		double [][] dbg_img = null;
		String [] dbg_titles = null; 
		showDoubleFloatArrays sdfa_instance = null;
		DisparityProcessor dp = new DisparityProcessor(this, clt_parameters.transform_size * geometryCorrection.getScaleDzDx());
		double [] disparity = scan.getDisparity(); // to modify in-place
		boolean [] these_no_border = new boolean [tlen];
		for (int i = 0; i < these_no_border.length; i++) {
			these_no_border[i] = last_scan.selected[i] && !last_scan.border_tiles[i]; 
		}
		
		boolean [] dbg_last_selected = last_scan.selected.clone();
		boolean [] dbg_last_border =   last_scan.border_tiles.clone();
		boolean [] dbg_no_border =     these_no_border.clone();
		
		boolean [] known_tiles = these_no_border.clone();
		// known are background or these tiles
		for (int i = 0; i < known_tiles.length; i++) {
			known_tiles[i] |= bg_scan.selected[i];
		}
		// set combo disparity from  last prepared
		for (int nt = 0; nt < known_tiles.length; nt++){
			int ty = nt / tilesX;
			int tx = nt % tilesX;
			disparity[nt] = 0.0;
///			if (last_scan.selected[nt]) disparity[nt] = last_scan.disparity[ty][tx];  
			if (these_no_border[nt]) disparity[nt] = last_scan.disparity[ty][tx];  
		}
		
		
		boolean [] grown = known_tiles.clone();
		growTiles(
				2*grow_sweep,          // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				grown,      // boolean [] tiles,
				null);     // boolean [] prohibit)
		boolean [] border = grown.clone();
		for (int i = 0; i < border.length; i++) border[i] &= !known_tiles[i];
//		double [] dbg_before = scan.getDisparity().clone();
		scan.fixNaNDisparity(
				border, // boolean [] select,   // which tiles to correct (null - all)
				scan.getDisparity(), // double [] disparity,
				scan.getStrength()); // double [] strength)
//		double [] dbg_after = scan.getDisparity().clone();
		
		int [] neighbors = dp.getNeighbors( // creates neighbors mask from bitmask
				grown, // these_tiles, // grown, // these_tiles, // boolean [] selected,
				tilesX);
/*		
		if (clt_parameters.show_neighbors) {
			double [] dbg_neib = dp.dbgShowNeighbors(
					grown, // these_tiles, // grown, // these_tiles,
					neighbors, // _orig, // int [] neighbors,
					clt_parameters.transform_size, // int    tile_size,
					-1.0, // double bgnd,
					1.0); // double fgnd)
			sdfa_instance.showArrays(dbg_neib,tilesX*clt_parameters.transform_size, tilesY*clt_parameters.transform_size,"XXneighbors");
		}
*/
		dp.smoothDisparity(
				clt_parameters.tiDispPull,   // final double     dispPull, // clt_parameters.tiDispPull or 0.0
				2, // 2, // 3,                           // final int        mask,     // 1 - work on internal elements, 2 - on border elements, 3 - both (internal first);
				clt_parameters.tiIterations, //  final int        num_passes,
				Math.pow(10.0,  -clt_parameters.tiPrecision), // final double     maxDiff, // maximal change in any of the disparity values
				neighbors,                                    // final int     [] neighbors, // +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
				scan.getDisparity(),                          // final double  [] disparity,          // current disparity value
				scan.getDisparity().clone(),                  // final double  [] measured_disparity, // measured disparity
				scan.getStrength(),                           // final double  [] strength,
				null, // this_hor_disparity,                  // final double     hor_disparity, // not yet used
				null, // hor_strength_conv,                   // final double     hor_strength, // not yet used
				known_tiles, // these_tiles,                                  // grown, // these_tiles,                 // final boolean [] selected,
				border,                                       // final boolean [] border,
				clt_parameters,
				threadsMax,                                   // maximal number of threads to launch                         
				debugLevel);
		
		// increase disparity if it
		boolean [] tried_before = new boolean [tlen];
		int num_tried_before = 0;
		int num_extended = 0;
		for (int num_scan = 0; num_scan < clt_3d_passes.size(); num_scan++){
			CLTPass3d pass = clt_3d_passes.get(num_scan);
//			if (pass != scan) { // not the same reference
			if (pass.isMeasured()) { // measured with pass.disparity[ty][tx]
				
				for (int ty = 0; ty < tilesY; ty++) for (int tx = 0; tx < tilesX; tx++){
					int nt = ty * tilesX + tx;
					if ((debugLevel >- 1) && (nt == dbg_tile)){
						System.out.println("setupExtendDisparity(): nt = "+nt);
					}
					if (border[nt]) {
						num_extended++;
						if (	(pass.tile_op[ty][tx] != 0) &&
								(Math.abs(pass.disparity[ty][tx] - disparity[nt]) < tried_margin)) {
							if (!tried_before[nt]) num_tried_before ++;
							tried_before[nt] = true; 
						}
					}
				}
			}
		}
		if (debugLevel > -1) {
			System.out.println("Number of extended tiles: "+num_extended+", tiles already tried (within +/-"+tried_margin+") :"+num_tried_before);
		}
		if (num_extended == 0){
			if (debugLevel > -1) {
				System.out.println("No new tiles wanted");
			}
			return 0;
		}
		int num_replaced = 0;
		if (num_tried_before > 0){
			for (int ty = 0; ty < tilesY; ty++) for (int tx = 0; tx < tilesX; tx++){
				int nt = ty * tilesX + tx;
				if (tried_before[nt]){
					double new_disp = scan.max_tried_disparity[ty][tx] + grow_disp_step;
					if (new_disp > grow_disp_max){
						// border[nt] = false; // let it repeat, will be filtered out if too close anyway. And removing border needs to remove selected too
					} else {
						num_replaced ++;
						disparity[nt] = new_disp;
// TODO: Make it smarter to give up looking at differences between images						
					}
				}
			}
		}
		if (debugLevel > -1) {
			System.out.println("Replaced "+num_replaced+" already tried tiles (of "+num_tried_before+"), abandoned "+(num_tried_before -num_replaced));
		}
		num_extended += num_replaced - num_tried_before;
		if (num_extended == 0){
			if (debugLevel > -1) {
				System.out.println("No new tiles wanted");
			}
			return 0;
		}
		scan.selected =     grown;
		scan.border_tiles = border;
		scan.disparity = new double [tilesY][tilesX];
		scan.tile_op =       new int [tilesY][tilesX];
		int op = ImageDtt.setImgMask(0, 0xf);
		op =     ImageDtt.setPairMask(op,0xf);
		op =     ImageDtt.setForcedDisparity(op,true);

		for (int ty = 0; ty < tilesY; ty++) for (int tx = 0; tx <tilesX; tx++){
			int indx =  tilesX * ty + tx;
			if (scan.selected[indx]) {
				scan.disparity[ty][tx] = disparity[indx];
				scan.tile_op[ty][tx] = op;
			} else {
				scan.disparity[ty][tx] = 0.0;
				scan.tile_op[ty][tx] = 0;
			}
		}
		
		
		if (show_debug){
			String [] dbg_titles0 = {"tried","disparity","bgnd","these","known_in", "known","last_sel_border", "no_border","tried_before"};
			dbg_titles = dbg_titles0; 
			dbg_img = new double[dbg_titles.length][];
			dbg_img[0] = new double[tilesY * tilesX];
//			dbg_img[1] = scan.getDisparity();
			dbg_img[1] = new double[tilesY * tilesX];
			
			dbg_img[2] = new double[tilesY * tilesX];
			dbg_img[3] = new double[tilesY * tilesX];
			dbg_img[4] = new double[tilesY * tilesX];
			dbg_img[5] = new double[tilesY * tilesX];

			dbg_img[6] = new double[tilesY * tilesX];
			dbg_img[7] = new double[tilesY * tilesX];
			dbg_img[8] = new double[tilesY * tilesX];

			for (int i = 0; i <tlen; i++){
				int ty = i / tilesX;
				int tx = i % tilesX;
				dbg_img[0][i] = scan.max_tried_disparity[ty][tx];
				dbg_img[1][i] = scan.disparity[ty][tx];
//				dbg_img[1][i] = last_scan.disparity[ty][tx];
				dbg_img[2][i] = ((bg_scan.selected != null) &&    (bg_scan.selected[i]))? 1.0:0.0 ;
				dbg_img[3][i] = ((last_scan.selected != null) &&  (last_scan.selected[i]))? 1.0:0.0 ;
				dbg_img[4][i] = dbg_img[2][i] + dbg_img[3][i];
				dbg_img[5][i] = (scan.selected[i]?1:0)+ (scan.border_tiles[i]?2:0);
				
				dbg_img[6][i] = (dbg_last_selected[i]?1:0)+ (dbg_last_border[i]?2:0);
				dbg_img[7][i] = dbg_no_border[i]?1:0;
				dbg_img[8][i] = tried_before[i]?1:0;
			}

			sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
			sdfa_instance.showArrays(dbg_img,  tilesX, tilesY, true, "extend_disparity",dbg_titles);
		}
/*
		boolean [] dbg_last_selected = last_scan.selected.clone();
		boolean [] dbg_last_border =   last_scan.border_tiles.clone();
		boolean [] dbg_no_border =     these_no_border.clone();
		
 */
		return num_extended;
	}
	
	public void showScan(
			CLTPass3d   scan,
			String title)
	{
		showDoubleFloatArrays sdfa_instance = null;
		sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		String [] titles = {
				"tile_op",       //  0
				"final",         //  1 - calculated, filetered, combined disparity
				"disparity",     //  2
				"disp_cm",       //  3
				"disp_hor",      //  4
				"disp_vert",     //  5
				"final_strength",//  6
				"strength",      //  7
				"strength_hor",  //  8
				"strength_vert", //  9
				"selection",     // 10
				"border_tiles",  // 11
				"max_tried"};    // 12
		int tlen = tilesX*tilesY;
//		double [][] dbg_img = new double[titles.length][tlen];
		double [][] dbg_img = new double[titles.length][];
		if (scan.tile_op != null)              dbg_img[ 0] = new double[tlen];
		if (scan.disparity !=  null)           dbg_img[ 2] = new double[tlen];
		if (scan.selected != null)             dbg_img[10] = new double[tlen];
		if (scan.border_tiles != null)         dbg_img[11] = new double[tlen];
		if (scan.max_tried_disparity != null)  dbg_img[12] = new double[tlen];
		
		for (int ty = 0; ty < tilesY; ty++) for (int tx = 0; tx < tilesX; tx++){
			int nt = ty*tilesX + tx;
			if (scan.tile_op != null)              dbg_img[ 0][nt] = scan.tile_op[ty][tx];
			if (scan.disparity !=  null)           dbg_img[ 2][nt] = scan.disparity[ty][tx];
			if (scan.selected != null)             dbg_img[10][nt] = scan.selected[nt]? 1.0:0.0;
			if (scan.border_tiles != null)         dbg_img[11][nt] = scan.border_tiles[nt]? 1.0:0.0;
			if (scan.max_tried_disparity != null)  dbg_img[12][nt] = scan.max_tried_disparity[ty][tx];
		}
		if (scan.disparity_map != null){
			dbg_img[3] = scan.disparity_map[ImageDtt.DISPARITY_INDEX_CM];
			dbg_img[4] = scan.disparity_map[ImageDtt.DISPARITY_INDEX_HOR];
			dbg_img[5] = scan.disparity_map[ImageDtt.DISPARITY_INDEX_VERT];
			dbg_img[7] = scan.disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX];
			dbg_img[8] = scan.disparity_map[ImageDtt.DISPARITY_INDEX_HOR_STRENGTH];
			dbg_img[9] = scan.disparity_map[ImageDtt.DISPARITY_INDEX_VERT_STRENGTH];
		}
		dbg_img[1] = scan.calc_disparity_combo;
		dbg_img[6] = scan.strength;
		
		title += "-";
		if (scan.isMeasured())  title += "M";
		if (scan.isProcessed()) title += "P";
		if (scan.isCombo())     title += "C";
		sdfa_instance.showArrays(dbg_img,  tilesX, tilesY, true, title,titles);
		System.out.println("showScan("+title+"): isMeasured()="+scan.isMeasured()+", isProcessed()="+scan.isProcessed()+", isCombo()="+scan.isCombo());
		
	}
	
	//sure_smth	  
	public boolean [] getBackgroundMask( // which tiles do belong to the background
			double     bgnd_range,       // disparity range to be considered background
			double     bgnd_sure,        // minimal strength to be considered definitely background
			double     bgnd_maybe,       // maximal strength to ignore as non-background
			double     sure_smth,        // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
			int        min_clstr_seed,   // number of tiles in a cluster to seed (just background?)
			int        min_clstr_block,  // number of tiles in a cluster to block (just non-background?)
			int        disparity_index,  // index of disparity value in disparity_map == 2 (0,2 or 4)
			boolean    show_bgnd_nonbgnd,
			int        debugLevel
			){
		boolean [] bgnd_tiles =       new boolean [tilesY * tilesX];
		boolean [] nonbgnd_tiles =    new boolean [tilesY * tilesX];
		boolean [] block_propagate =  new boolean [tilesY * tilesX];
		int quad = 4;
		showDoubleFloatArrays sdfa_instance = null;
		if (debugLevel > -1) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?

		//		  boolean [] new_tiles =        new boolean [tilesY * tilesX]; // grow selection by 1 tile over non-background?
		CLTPass3d bgnd_data = clt_3d_passes.get(0); 
		//		  double [][][][] texture_tiles = bgnd_data.texture_tiles;
		double [][]     disparity_map=  bgnd_data.disparity_map;
		for (int tileY = 0; tileY < tilesY; tileY++){
			for (int tileX = 0; tileX < tilesX; tileX++){
				int tindx = tileY * tilesX + tileX;

				if (Math.abs(disparity_map[disparity_index][tindx]) < bgnd_range){
					if (disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX][tindx] >= bgnd_sure){
						bgnd_tiles[tindx] = true;  
					}
				} else {
					if (disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX][tindx] > bgnd_maybe){ // maybe non-bkgnd
						if (disparity_map[disparity_index][tindx] > 0.0) { // disregard negative disparity points
							nonbgnd_tiles[tindx] = true;  
						}
					}
				}
				// see if the second worst variation exceeds sure_smth (like a window), really close object
				int imax1 = 0;
				for (int i = 1; i< quad; i++){
					if (disparity_map[ImageDtt.IMG_DIFF0_INDEX+i][tindx] > disparity_map[ImageDtt.IMG_DIFF0_INDEX + imax1][tindx]) imax1 = i;
				}
				int imax2 = (imax1 == 0)? 1 : 0;
				for (int i = 0; i< quad; i++) if (i != imax1) {
					if (disparity_map[ImageDtt.IMG_DIFF0_INDEX+i][tindx] > disparity_map[ImageDtt.IMG_DIFF0_INDEX + imax2][tindx]) imax2 = i;
				}
				block_propagate[tindx] = (disparity_map[ImageDtt.IMG_DIFF0_INDEX + imax2][tindx] > sure_smth);
			}
		}
		// TODO: check if minimal cluster strengh should be limited here
		if (min_clstr_seed > 1){
			removeSmallClusters(
					true,            // boolean    diag_en,   // enable diagonal directions, false only up, dowm, right,left
					bgnd_tiles,      // boolean [] tiles_src,    // selected tiles, will modified
					null,            // double []   weights_src,   // or null
					min_clstr_seed,  // int        min_area,  // minimal number of pixels
					0.0, //clt_parameters.min_clstr_weight,  // double     min_weight // minimal total weight of the cluster
					0.0); // clt_parameters.min_clstr_max);    // double     min_max_weight // minimal value of the maximal strengh in the cluster
		}		  

		if (min_clstr_block > 1){
			removeSmallClusters(
					true,            // boolean    diag_en,   // enable diagonal directions, false only up, dowm, right,left
					nonbgnd_tiles,   // boolean [] tiles_src,    // selected tiles, will modified
					null,            // double []   weights_src,   // or null
					min_clstr_block, // int        min_area,  // minimal number of pixels
					0.0, //clt_parameters.min_clstr_weight,  // double     min_weight // minimal total weight of the cluster
					0.0); // clt_parameters.min_clstr_max);    // double     min_max_weight // minimal value of the maximal strengh in the cluster
		}

		if ((sdfa_instance != null) && show_bgnd_nonbgnd) {
			String [] titles = {"bgnd","nonbgnd","block","strength","disparity"};
			double [][] dbg_img = new double[titles.length][tilesY * tilesX];
			for (int i = 0; i<dbg_img[0].length;i++){
				dbg_img[0][i] =    bgnd_tiles  [i] ? 1 : 0;
				dbg_img[1][i] = nonbgnd_tiles  [i] ? 1 : 0;
				dbg_img[2][i] = block_propagate[i] ? 1 : 0;
				dbg_img[3][i] = disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX][i];
				dbg_img[4][i] = disparity_map[disparity_index][i];
			}
			sdfa_instance.showArrays(dbg_img,  tilesX, tilesY, true, "bgnd_nonbgnd",titles);
		}

		for (int gain = 1; gain > 0;){
			gain = 0;
			// extend to the right
			for (int tileY = 0; tileY < tilesY; tileY++){
				for (int tileX = 0; tileX < (tilesX - 1); tileX++){
					int tindx = tileY * tilesX + tileX;
					if (bgnd_tiles[tindx] && !bgnd_tiles[tindx + 1] && !nonbgnd_tiles[tindx + 1] && !block_propagate[tindx + 1]){
						bgnd_tiles[tindx+1] = true;
						gain++;
					}
				}
			}
			// extend to the left
			for (int tileY = 0; tileY < tilesY; tileY++){
				for (int tileX = tilesX - 1; tileX > 1; tileX--){
					int tindx = tileY * tilesX + tileX;
					if (bgnd_tiles[tindx] && !bgnd_tiles[tindx - 1] && !nonbgnd_tiles[tindx -1]  && !block_propagate[tindx - 1]){
						bgnd_tiles[tindx - 1] = true;
						gain++;
					}
				}
			}
			// extend down
			for (int tileX = 0; tileX < tilesX; tileX++){
				for (int tileY = 0; tileY < (tilesY -1); tileY++){
					int tindx = tileY * tilesX + tileX;
					if (bgnd_tiles[tindx] && !bgnd_tiles[tindx + tilesX] && !nonbgnd_tiles[tindx + tilesX] && !block_propagate[tindx + tilesX]){
						bgnd_tiles[tindx + tilesX] = true;
						gain++;
					}
				}
			}
			// extend up
			for (int tileX = 0; tileX < tilesX; tileX++){
				for (int tileY = tilesY - 1; tileY > 1; tileY--){
					int tindx = tileY * tilesX + tileX;
					if (bgnd_tiles[tindx] && !bgnd_tiles[tindx - tilesX] && !nonbgnd_tiles[tindx - tilesX] && !block_propagate[tindx - tilesX]){
						bgnd_tiles[tindx - tilesX] = true;
						gain++;
					}
				}
			}
			if (debugLevel > -1) {
				System.out.println("getBackgroundMask(), gain="+gain);
			}
			if (sdfa_instance!=null){
				double [] dbg_img = new double[tilesY * tilesX];
				for (int i = 0; i<dbg_img.length;i++){
					dbg_img[i] =    bgnd_tiles[i]?1:0;
				}
				sdfa_instance.showArrays(dbg_img,  tilesX, tilesY, "tiles");
			}

		}
		return bgnd_tiles;
	}
	
	public int extendIndex(int indx){
		int tileY = indx / tilesX;
		int tileX = indx % tilesX;
		return (tilesX + 2) * (tileY + 1) + (tileX + 1);
		
	}
	
	
/**
 * Create tasks for individual clusters (with flaps). Border_fixed will have alpha = 0 and provided disparity, border float - adjusted disparity from neighbors
 * No ordering here, just filtering
 * @param maxClusters
 * @param minClusterArea
 * @param disparity
 * @param strength
 * @param minStrength
 * @param overlap_clusters
 * @param minDisparity
 * @param maxDisparity
 * @param debugLevel
 * @return
 */
	public int createTileOverlapTasks(
			int         maxClusters,
			int         minClusterArea,
			double []   disparity_in, 
			double []   strength_in,
			double      maxChange,   // adjust border disparity until change is below this.
			double      minStrength,
			int [][][]  overlap_clusters, // [cluster] {cd.internal, cd.border_fixed, cd.border_float}
			double    minDisparity,
			double    maxDisparity,
			boolean   show_shells,
			int       debugLevel)
	{
		final int maxrep=1000;
		showDoubleFloatArrays sdfa_instance = null;
		if (debugLevel > -1) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		int tilesX2 = tilesX + 2;
		int tilesY2 = tilesY + 2;
		int tlen2 = tilesX2*tilesY2;
		double [] disparity0 = new double [tlen2];
		double [] strength0 = new double [tlen2];

		for (int i = 0; i < tilesY; i++) {
			for (int j = 0; j < tilesX; j++) {
				strength0[(i + 1) * tilesX2 + j + 1] =  strength_in[i * tilesX + j];
				disparity0[(i + 1) * tilesX2 + j + 1] = disparity_in[i * tilesX + j];
			}
		}

		int [] dirs8 = {-tilesX2,  -tilesX2 + 1, 1, tilesX2 +1, tilesX2, tilesX2 - 1, -1, -tilesX2 - 1};
		int [] dirs = dirs8;

		int op = ImageDtt.setImgMask(0, 0xf);
		op =     ImageDtt.setPairMask(op,0xf);
		op =     ImageDtt.setForcedDisparity(op,true);

//		int tlen = clusters.length;
		int numClusters = 0;
//		int leng = disparity_in.length;
		double  [] disparity = new double  [tlen2];
		boolean [] selected =  new boolean [tlen2]; 
		boolean [] border =    new boolean [tlen2]; 
		boolean [] fixed =     new boolean [tlen2];
		boolean [] none4 =     new boolean [tlen2];
		int startPass = clt_3d_passes.size();
		for (int ncl = 0; (ncl < overlap_clusters.length) && ((maxClusters == 0) || (numClusters < maxClusters)); ncl++){ // need to break;
			if (overlap_clusters[ncl][0].length < minClusterArea) {
				if (debugLevel > -1) {
					System.out.println("createTileOverlapTasks(): skipping small cluster "+ncl+", internal area "+(overlap_clusters[ncl][0].length)+" < "+minClusterArea);
				}
				continue; // skip this cluster
			}

			System.arraycopy(disparity0, 0, disparity, 0, tlen2);
			System.arraycopy(none4,      0, selected,  0, tlen2);
			System.arraycopy(none4,      0, border,    0, tlen2);
			System.arraycopy(none4,      0, fixed,     0, tlen2);
			double clusterMaxStrength = 0.0;
			for (int i = 0; i < overlap_clusters[ncl][0].length; i++) {
				int indx = extendIndex(overlap_clusters[ncl][0][i]);
				selected [indx] = true;
				fixed [ indx] = true;
				if (strength0[indx] > clusterMaxStrength) clusterMaxStrength = strength0[indx];
			}
			if (clusterMaxStrength < minStrength) {
				if (debugLevel > -1) {
					System.out.println("createTileOverlapTasks(): skipping weak cluster "+ncl+", strength "+clusterMaxStrength+" < "+minStrength);
				}
				continue; // skip this cluster
			}
			
			for (int i = 0; i < overlap_clusters[ncl][1].length; i++) {
				int indx = extendIndex(overlap_clusters[ncl][1][i]);
				border [indx] = true;
				fixed [ indx] = true;
			}
			
			if (overlap_clusters[ncl][2].length > 0) {
				int [] indx2 = new int [overlap_clusters[ncl][2].length];
				for (int i = 0; i< indx2.length; i++) indx2[i] = extendIndex(overlap_clusters[ncl][2][i]);
				
			
				for (int i = 0; i < indx2.length; i++) {
					int indx = indx2[i];
					border [indx] = true;
					double sw = 0.0, sd = 0.0;
					for (int dir = 0; dir < dirs.length; dir++){
						if (fixed[indx + dirs[dir]]){
							sw += 1.0;
							sd += disparity[indx + dirs[dir]];
						}
					}
					if (sw >0) sd/= sw;
//////					disparity[indx] = sw;
					disparity[indx] = sd;
				}
				double [] new_disparity = new double [indx2.length];
				for (int nrep = 0; nrep < maxrep; nrep++){
					for (int i = 0; i < indx2.length; i++) {
						int indx =indx2[i];
						double sw = 0.0, sd = 0.0;
						for (int dir = 0; dir < dirs.length; dir++){
							if ((fixed[indx + dirs[dir]]) || (border[indx + dirs[dir]])) {
								sw += 1.0;
								sd += disparity[indx + dirs[dir]];
							}
						}
						if (sw > 0) sd/= sw;
						new_disparity[i] = sd;
					}
					double adiff = 0.0;
					for (int i = 0; i < indx2.length; i++) {
						int indx =indx2[i];
						double ad = Math.abs(new_disparity[i] - disparity[indx]);
						if (ad > adiff) adiff = ad;
						disparity[indx] = new_disparity[i];
					}
					if (adiff < maxChange) break;
				}
			}
			
			// now fill the task
			// return to original dimensions
			double [][] disparityTask = new double [tilesY][tilesX];
			int [][]    tile_op =       new int [tilesY][tilesX];
			boolean [] borderTiles =    new boolean[tilesY*tilesX]; // to zero alpha in the images
			for (int ty = 0; ty < tilesY; ty++) for (int tx = 0; tx <tilesX; tx++){
				int indx =  tilesX * ty + tx;
				int indx4 = tilesX2 * (ty + 1) + (tx + 1);
				if (selected[indx4] || border[indx4]){
					if      (disparity[indx4] < minDisparity) disparity[indx4] = minDisparity;
					else if (disparity[indx4] > maxDisparity) disparity[indx4] = maxDisparity;
					
					disparityTask[ty][tx] = disparity[indx4];
					tile_op[ty][tx] = op;
					borderTiles[indx] = border[indx4];
				} else {
					disparityTask[ty][tx] = 0.0;
					tile_op[ty][tx] =       0;
					borderTiles[indx] =     false;
				}
			}
			
			// Create FPGA task for this cluster
//		public int             dbg_index;

			CLTPass3d scan_next =new CLTPass3d(this);
			scan_next.disparity =     disparityTask;
			scan_next.tile_op =       tile_op;
			scan_next.border_tiles =  borderTiles;
			scan_next.dbg_index =     ncl; // matching overlapClusters[] 
			clt_3d_passes.add(scan_next);

			if (debugLevel > -1) {
				System.out.println("createTileOverlapTasks(): Added cluster "+ncl+", internal area = "+(overlap_clusters[ncl][0].length)+ ", max strength=clusterMaxStrength"+
			", fixed border: "+(overlap_clusters[ncl][1].length)+ ", float border: "+(overlap_clusters[ncl][2].length));
			}
			
			numClusters++;
			
		}
		if (debugLevel > -1) {
			System.out.println("createTileOverlapTasks(): prepared "+numClusters+" clusters");
		}
		if ((debugLevel > -1) && show_shells) {
			double [][] dbg_shells = new double [clt_3d_passes.size() - startPass][tilesY*tilesX+1];
			double ampl = dbg_shells.length;
			for (int ns = 1; ns < dbg_shells.length; ns++) {
				CLTPass3d dbg_scan = clt_3d_passes.get(ns -1 + startPass);
				for (int ty = 0; ty < tilesY; ty++){
					for (int tx = 0; tx < tilesX; tx++){
						int indx = ty*tilesX + tx;
						if (dbg_scan.tile_op[ty][tx] !=0){
							if (dbg_scan.border_tiles[indx]) {
								dbg_shells[ns][indx] = 0.5*ampl;
							} else {
								dbg_shells[ns][indx] = 1.0*ampl;
								dbg_shells[ 0][indx] = ns;
							}
						}
					}
				}
			}
			String [] titles = new String[dbg_shells.length];
			titles[0] = "index";
			for (int i = 1; i <titles.length; i++) titles[i] = "shell "+i;
			sdfa_instance.showArrays(dbg_shells, tilesX, tilesY, true, "shells",titles);
		}

		return numClusters;  
	}
/**
 * Return fade alpha for the border. First index of the array (0..15) is a bitmask (1 - N, 2 - E, 4 - S, 8 - W) of a solid
 * (non-border) tiles around this border ones
 * Solid tiles on 2 opposite sides return alpha = 1.0;
 * Solid on 2 corner ones - a "roof" - use maximum of 2
 * Solid on 1 side - 0.5 * (cosine + 1.0)
 * @param transform_size
 * @return
 */
	public double [][] getAlphaFade(
			int transform_size)
	{
		int ts2 = 2 * transform_size;
		int ts2m1 = ts2-1;
		double [][] alphaFade = new double[16][ts2*ts2];
		double [] fade1d = new double [ts2];
		for (int i = 0; i <  ts2; i++){
			fade1d[i] = 0.5 * (1.0 - Math.cos(Math.PI * (i +0.5) /ts2));
		}
		for (int i = 0; i < ts2; i++){
			for (int j = 0; j < ts2; j++){
				int indx = i * ts2 + j;
				for (int m = 0; m < 16; m++){
					switch (m){
					case  1:
						alphaFade[m][indx] = fade1d[ts2m1 - i];
						break;
					case  2:
						alphaFade[m][indx] = fade1d[j];
						break;
					case  4:
						alphaFade[m][indx] = fade1d[i];
						break;
					case  8:
						alphaFade[m][indx] = fade1d[ts2m1 - j];
						break;
					case  3:
						alphaFade[m][indx] = (j > ts2m1 - i) ? fade1d[j]: fade1d[ts2m1 - i];
						break;
					case  6:
						alphaFade[m][indx] = (j > i) ?         fade1d[j]: fade1d[i];
						break;
					case  9:
						alphaFade[m][indx] = (j > i) ?         fade1d[ts2m1 - i]: fade1d[ts2m1 - j];
						break;
					case  12:
						alphaFade[m][indx] = (i > ts2m1 - j) ? fade1d[i]: fade1d[ts2m1 - j];
						break;
					case  5:
					case  7:
					case 10:
					case 11:
					case 13:
					case 14:
					case 15:
						alphaFade[m][indx] = 1.0;
						break;
// 0 is impossible, let alpha b0 0 then						
					}
				}
				
			}
		}
		return alphaFade;
	}
	
	

	public int createTileTasks(
			int       maxClusters,
			int       minClusterArea,
			double [] disparity_in,
			double [] strength_in,
			double    minStrength,
			int []    clusters_in,
			double    minDisparity,
			double    maxDisparity,
			int       debugLevel)
	{
		showDoubleFloatArrays sdfa_instance = null;
		if (debugLevel > -1) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?

		// adding 1-tile frame around to avoid checking for the borders
		int tilesX4 = tilesX + 4;
		int tilesY4 = tilesY + 4;
		int tlen4 = tilesX4*tilesY4;
		double [] disparity0 = new double [tlen4];
		double [] strength0 = new double [tlen4];
		int [] clusters = new int [tlen4];
		for (int i = 0; i < tlen4; i++) {
			clusters[i] = 0;
			strength0[i] = 0.0;
		}

		for (int i = 0; i < tilesY; i++) {
			for (int j = 0; j < tilesX; j++) {
				clusters[(i + 2) * tilesX4 + j + 2] =  clusters_in[i * tilesX + j];
				strength0[(i + 2) * tilesX4 + j + 2] =  strength_in[i * tilesX + j];
				disparity0[(i + 2) * tilesX4 + j + 2] = disparity_in[i * tilesX + j];
			}
		}


		int [] dirs8 = {-tilesX4,  -tilesX4 + 1, 1, tilesX4 +1, tilesX4, tilesX4 - 1, -1, -tilesX4 - 1};
		int [] dirs = dirs8;
		int op = ImageDtt.setImgMask(0, 0xf);
		op =     ImageDtt.setPairMask(op,0xf);
		op =     ImageDtt.setForcedDisparity(op,true);

		int tlen = clusters.length;
		int numClusters = 0;
		for (int ncl = 0; (maxClusters == 0) || (ncl < maxClusters); ncl++){ // need to break;

			double [] strength = strength0.clone();
			double [] disparity = disparity0.clone();
			boolean [] bcluster = new boolean [tlen];
			int clusterSize = 0;
			int cindx = ncl + 1;
			for (int i = 0; i < tlen; i++) if (clusters[i] == cindx) {
				clusterSize++;
				bcluster[i] = true;
			}
			if (debugLevel > -1) {
				System.out.println("createTileTasks(), cluster #"+ncl+", cluster size = "+clusterSize);
			}
			if ((clusterSize < minClusterArea) || (clusterSize == 0)) break;
			boolean [] grown_clusters = bcluster.clone();
			growTiles(
					2,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
					grown_clusters,
					null,
					tilesX4,
					tilesY4);
			// fill missing disparities by a simple method - actually should be done before
			int numFailures = 0;
			int numSuccesses = 0;
			if ( (ncl ==0) && (debugLevel > 0)) {
				String [] titles = {"strength","disparity","clusters","grown"};
				double [][] dbg_img = new double [titles.length][tilesX4*tilesY4];
				for (int i = 0; i < clusters.length; i++){
					dbg_img[2][i] = clusters[i];
					dbg_img[3][i] = grown_clusters[i]?1:0;
				}
				dbg_img[0] = strength0;
				dbg_img[1] = disparity0;
				sdfa_instance.showArrays(dbg_img,  tilesX4, tilesY4,true, "createTileTasks",titles);
			}

			while (true){
				numFailures = 0;
				numSuccesses = 0;
				for (int i = 0; i < tlen ; i++){
					if (bcluster[i] && (Double.isNaN(disparity[i]) || (strength[i] < minStrength) || (disparity[i] < minDisparity)  || (disparity[i] > maxDisparity))){
						double sw =0.0, sd = 0.0;
						int n = 0;
						for (int d = 0; d < dirs.length; d++){
							int indx = i + dirs[d];
							if (bcluster[indx] && !Double.isNaN(disparity[indx]) && (strength[indx] > minStrength)){
								sw += strength[indx];
								sd += strength[indx] * disparity[indx];
								n++;
							}
						}
						if (sw >0){
							numSuccesses++;
							disparity[i] = sd/sw;
							strength[i] = sw/n;
						} else {
							numFailures++;
						}
					}
				}
				if (debugLevel > 0) {
					System.out.println("createTileTasks(), numFailures="+numFailures+", numSuccesses = "+numSuccesses);
				}

				if ((numFailures == 0) || (numSuccesses == 0))  break;
			}



			if (numFailures > 0){
				System.out.println("**** Could not fill some missing tiles for cluster # " + cindx + ", numFailures = "+numFailures+" ****");
				continue; // with next cluster
			}

			// now fill border pixels disparity disregarding known for that tile
			for (int i = 0; i < tlen ; i++){
				if (!bcluster[i] && grown_clusters[i]){ // only border
					double sw =0.0, sd = 0.0;
					int n = 0;
					for (int d = 0; d < dirs.length; d++){
						int indx = i + dirs[d];
						if (bcluster[indx]){
							sw += strength[indx];
							sd += strength[indx] * disparity[indx];
							n++;
						}
					}
					if (sw >0){
						disparity[i] = sd/sw;
						strength[i] = sw/n; // will not actually be used
					} else {
						System.out.println("**** Program BUG, should not happen ****");
						disparity[i] = 0.0;
						strength[i] = 0.0;
					}
				}
			}
			// return to original dimensions
			double [][] disparityTask = new double [tilesY][tilesX];
			int [][]    tile_op =       new int [tilesY][tilesX];
			boolean [] borderTiles =    new boolean[tilesY*tilesX]; // to zero alpha in the images
			for (int ty = 0; ty < tilesY; ty++) for (int tx = 0; tx <tilesX; tx++){
				int indx =  tilesX * ty + tx;
				int indx4 = tilesX4 * (ty+2) + (tx + 2);
				if (grown_clusters[indx4]){
					disparityTask[ty][tx] = disparity[indx4];
					tile_op[ty][tx] = op;
					borderTiles[indx] = !bcluster[indx4];
				} else {
					disparityTask[ty][tx] = 0.0;
					tile_op[ty][tx] =       0;
					borderTiles[indx] =     false;
				}
			}
			// Create FPGA task for this cluster

			CLTPass3d scan_next =new CLTPass3d(this);
			scan_next.disparity =     disparityTask;
			scan_next.tile_op =       tile_op;
			scan_next.border_tiles =  borderTiles;
			clt_3d_passes.add(scan_next);
			numClusters++;
		}
		return numClusters;  
	}

	
	

	public int createTileTasks(
			int       maxClusters,
			int       minClusterArea,
			double [] disparity_in,
			int []    clusters_in,
			int       debugLevel)
	{
		showDoubleFloatArrays sdfa_instance = null;
		if (debugLevel > -1) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?

		// adding 1-tile frame around to avoid checking for the borders
		int tilesX4 = tilesX + 4;
		int tilesY4 = tilesY + 4;
		int tlen4 = tilesX4*tilesY4;
		double [] disparity0 = new double [tlen4];
		int [] clusters = new int [tlen4];
		for (int i = 0; i < tlen4; i++) {
			clusters[i] = 0;
		}

		for (int i = 0; i < tilesY; i++) {
			for (int j = 0; j < tilesX; j++) {
				clusters[(i + 2) * tilesX4 + j + 2] =  clusters_in[i * tilesX + j];
				disparity0[(i + 2) * tilesX4 + j + 2] = disparity_in[i * tilesX + j];
			}
		}


		int [] dirs8 = {-tilesX4,  -tilesX4 + 1, 1, tilesX4 +1, tilesX4, tilesX4 - 1, -1, -tilesX4 - 1};
		int [] dirs = dirs8;
		int op = ImageDtt.setImgMask(0, 0xf);
		op =     ImageDtt.setPairMask(op,0xf);
		op =     ImageDtt.setForcedDisparity(op,true);

		int tlen = clusters.length;
		int numClusters = 0;
		for (int ncl = 0; (maxClusters == 0) || (ncl < maxClusters); ncl++){ // need to break;
			double [] disparity = disparity0.clone();
			boolean [] bcluster = new boolean [tlen];
			int clusterSize = 0;
			int cindx = ncl + 1;
			for (int i = 0; i < tlen; i++) if (clusters[i] == cindx) {
				clusterSize++;
				bcluster[i] = true;
			}
			if (debugLevel > -1) {
				System.out.println("createTileTasks(), cluster #"+ncl+", cluster size = "+clusterSize);
			}
			if ((clusterSize < minClusterArea) || (clusterSize == 0)) break;
			boolean [] grown_clusters = bcluster.clone();
			growTiles(
					2,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
					grown_clusters,
					null,
					tilesX4,
					tilesY4);

			// now fill border pixels disparity disregarding known for that tile
			for (int i = 0; i < tlen ; i++){
				if (!bcluster[i] && grown_clusters[i]){ // only border
					double sd = 0.0;
					int n = 0;
					for (int d = 0; d < dirs.length; d++){
						int indx = i + dirs[d];
						if (bcluster[indx]){
							sd += disparity[indx];
							n++;
						}
					}
					if (n > 0){
						disparity[i] = sd/n;
					} else {
						System.out.println("**** Program BUG, should not happen ****");
						disparity[i] = 0.0;
					}
				}
			}
			// return to original dimensions
			double [][] disparityTask = new double [tilesY][tilesX];
			int [][]    tile_op =       new int [tilesY][tilesX];
			boolean [] borderTiles =    new boolean[tilesY*tilesX]; // to zero alpha in the images
			for (int ty = 0; ty < tilesY; ty++) for (int tx = 0; tx <tilesX; tx++){
				int indx =  tilesX * ty + tx;
				int indx4 = tilesX4 * (ty+2) + (tx + 2);
				if (grown_clusters[indx4]){
					disparityTask[ty][tx] = disparity[indx4];
					tile_op[ty][tx] = op;
					borderTiles[indx] = !bcluster[indx4];
				} else {
					disparityTask[ty][tx] = 0.0;
					tile_op[ty][tx] =       0;
					borderTiles[indx] =     false;
				}
			}
			// Create FPGA task for this cluster

			CLTPass3d scan_next =new CLTPass3d(this);
			scan_next.disparity =     disparityTask;
			scan_next.tile_op =       tile_op;
			scan_next.border_tiles =  borderTiles;
			clt_3d_passes.add(scan_next);
			numClusters++;
		}
		return numClusters;  
	}

	
	
	
	
	
	// temporary  - simple separation until continuous tile areas are split along disparity steps
	public int [] enumerateClusters(
			boolean diag_en,
			boolean [] tiles_src)
	{
		// adding 1-tile frame around to avoid checking for the borders
		int tilesX2 = tilesX + 2;
		int tilesY2 = tilesY + 2;
		boolean [] tiles = new boolean [tilesX2*tilesY2];
		int [] enum_clust2 = new int[tiles.length];
		for (int i = 0; i < enum_clust2.length; i++) enum_clust2[i] = 0;
		int numClust = 0;
		for (int i = 0; i < tilesY2; i++) {
			tiles[tilesX2 * 0 + i] =         false;
			tiles[tilesX2*(tilesY2-1) + i] = false;
			tiles[tilesX2* i + 0] =          false;
			tiles[tilesX2* i + tilesX2 -1] = false;
		}
		for (int i = 0; i < tilesY; i++) {
			for (int j = 0; j < tilesX; j++) {
				tiles[(i+1)*tilesX2 + j + 1]= tiles_src[i * tilesX + j];
			}
		}

		int [] dirs4 = {-tilesX2, 1, tilesX2,-1};
		int [] dirs8 = {-tilesX2,  -tilesX2 + 1, 1, tilesX2 +1, tilesX2, tilesX2 - 1, -1, -tilesX2 - 1};
		int [] dirs = diag_en? dirs8 : dirs4;
		int [] waves = new int [tiles.length];
		for (int i = 0; i < tiles.length; i++) waves[i] = tiles[i] ? 0: -1;
		ArrayList<Integer> front = new ArrayList<Integer>();
		for (int start_indx = 0; start_indx < tiles.length; start_indx++) if (waves[start_indx] == 0){ // found first pixel of a new cluster
			numClust ++;
			Integer ipx = start_indx;
			Integer ipx1;
			front.clear();
			int area = 1;
			waves[ipx] = area;
			enum_clust2[ipx] = numClust;
			front.add(ipx);
			while (!front.isEmpty()) {
				ipx = front.remove(0);// get oldest element
				for (int d = 0; d < dirs.length; d++){
					ipx1 = ipx + dirs[d];
					if (waves[ipx1] == 0) {
						area++;
						waves[ipx1] = area;
						enum_clust2[ipx1] = numClust;
						front.add(ipx1);
					}
				}
			}
		}
		// return to original dimensions
		int [] enum_clust = new int[tiles_src.length];
		for (int i = 0; i < tilesY; i++) {
			for (int j = 0; j < tilesX; j++) {
				enum_clust[i * tilesX + j] =enum_clust2[(i+1)*tilesX2 + j + 1];
			}
		}
		// count cluster
		int []clustSizes = new int [numClust];
		for (int i = 0; i < clustSizes.length; i++) clustSizes[i] = 0;
		for (int i = 0; i < enum_clust.length; i++) if (enum_clust[i] > 0) clustSizes[enum_clust[i]-1]++;

		class Pair implements Comparable<Pair> {
			public final int index;
			public final int value;

			public Pair(int index, int value) {
				this.index = index;
				this.value = value;
			}

			@Override
			public int compareTo(Pair other) {
				return Integer.valueOf(this.value).compareTo(other.value);
			}
		}

		Pair[] pairs = new Pair[numClust];
		for (int i = 0; i < clustSizes.length; i++) pairs [i] = new Pair (i, clustSizes[i]);
		Arrays.sort(pairs);
		int [] revIndex = new int [numClust];
		for (int i = 0; i < revIndex.length; i++) revIndex [pairs[i].index] = (numClust - i); // array was in accending order
		int [] enum_clust_ordered = new int[tiles_src.length];
		for (int i=0; i < enum_clust_ordered.length; i++){
			enum_clust_ordered[i] = (enum_clust[i] > 0) ? revIndex[enum_clust[i] - 1] : 0; 
		}
		return enum_clust_ordered;  
	}

	public void fillGaps( // grows, then shrinks
			int depth, // same as grow - odd - 4 directions, even - 8
			boolean    poison, // do not fill gaps that even touch prohibited
			boolean [] tiles,
			boolean [] prohibit){
		boolean [] orig_tiles = tiles.clone();
		if (poison && (prohibit != null)) {
			
			growTiles(
					depth,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
					tiles,
					null);
			boolean [] poisoned = prohibit.clone();
			for (int i =0; i< tiles.length; i++) poisoned[i] &= tiles[i] && !orig_tiles[i];
			growTiles(
					depth,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
					poisoned,
					orig_tiles);
			
			// invert selection, so "grow" will be "shrink"
			for (int i =0; i< tiles.length; i++) tiles[i] = !tiles[i];
			growTiles(
					depth,
					tiles,
					null);
			// invert selection again (restore)
			for (int i =0; i< tiles.length; i++) tiles[i] = !tiles[i] && !poisoned[i];

		} else {
			growTiles(
					depth,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
					tiles,
					prohibit);
			// invert selection, so "grow" will be "shrink"
			for (int i =0; i< tiles.length; i++) tiles[i] = !tiles[i];
			growTiles(
					depth,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
					tiles,
					prohibit);
			// invert selection again (restore)
			for (int i =0; i< tiles.length; i++) tiles[i] = !tiles[i];
			//		  for (int i =0; i< tiles.length; i++) tiles[i] = !tiles[i] || orig_tiles[i]; // and or with original (why result does not have some of original - only with prohibit?
			if (prohibit != null) { // second pass w/o prohibit, then and results
				growTiles(
						depth,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
						orig_tiles,
						null);
				// invert selection, so "grow" will be "shrink"
				for (int i =0; i< tiles.length; i++) orig_tiles[i] = !orig_tiles[i];
				growTiles(
						depth,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
						orig_tiles,
						null);
				// invert selection again (restore)
				for (int i =0; i< tiles.length; i++) tiles[i] &= !orig_tiles[i];
			}
		}
	}
	public void growTiles(
			int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
			boolean [] tiles,
			boolean [] prohibit)
	{
		growTiles(
				grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				tiles,
				prohibit,
				this.tilesX,
				this.tilesY);
	}


	public void growTiles(
			int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
			boolean [] tiles,
			boolean [] prohibit,
			int tilesX,
			int tilesY)
	{
		boolean [] src_tiles = tiles.clone(); // just in case
		// grow
		boolean hor = true;
		for (; grow > 0; grow--){
			boolean single = (grow ==1) && hor;
			src_tiles = tiles.clone();
			if (hor){
				for (int tileY = 0; tileY < tilesY; tileY++){
					for (int tileX = 0; tileX < (tilesX - 1); tileX++){
						int tindx = tileY * tilesX + tileX;
						if ((prohibit == null) || (!prohibit[tindx] && !prohibit[tindx + 1])) {
							tiles[tindx + 1] |= src_tiles[tindx];  
						}

					}
					for (int tileX = 1; tileX < tilesX; tileX++){
						int tindx = tileY * tilesX + tileX;
						if ((prohibit == null) || (!prohibit[tindx] && !prohibit[tindx - 1])) {
							tiles[tindx - 1] |= src_tiles[tindx];
						}
					}
				}
			}
			if (!hor || single){ // do vertically, but from previous state
				for (int tileX = 0; tileX < tilesX; tileX++){
					for (int tileY = 0; tileY < (tilesY - 1); tileY++){
						int tindx = tileY * tilesX + tileX;
						if ((prohibit == null) || (!prohibit[tindx] && !prohibit[tindx + tilesX])) {
							tiles[tindx + tilesX] |= src_tiles[tindx];
						}

					}
					for (int tileY = 1; tileY < tilesY; tileY++){
						int tindx = tileY * tilesX + tileX;
						if ((prohibit == null) || (!prohibit[tindx] && !prohibit[tindx - tilesX])) {
							tiles[tindx - tilesX] |= src_tiles[tindx];
						}
					}
				}
			}
			hor = !hor;
		}
	}

	public void removeLoneClusters(
			boolean    diag_en,   // enable diagonal directions, false only up, dowm, right,left
			boolean [] tiles_src,    // selected tiles, will modified
			double []  weights_src,   // or null
			int        min_area,  // minimal number of pixels
			double     min_weight, // minimal total weight of the cluster
			double     min_max_weight // minimal value of the maximal strengh in the cluster
			){
		removeLoneClusters(
				diag_en,   // enable diagonal directions, false only up, dowm, right,left
				tiles_src,    // selected tiles, will modified
				weights_src,   // or null
				min_area,  // minimal number of pixels
				min_weight, // minimal total weight of the cluster
				min_max_weight,
				0// minimal value of the maximal strengh in the cluster
				);
	}
	
	public void removeLoneClusters(
			boolean    diag_en,   // enable diagonal directions, false only up, dowm, right,left
			boolean [] tiles,     // selected tiles, will modified
			double []  weights_src,   // or null
			int        min_area,  // minimal number of pixels
			double     min_weight, // minimal total weight of the cluster (expanded!
			double     min_max_weight, // minimal value of the maximal strengh in tghe cluster
			int        debugLevel
			){
		boolean [] grown_by_1 = tiles.clone();
		growTiles(2, // 1, // 2,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				grown_by_1,
				null);
		int grown_min = (int) (min_area + 4 * Math.sqrt(min_area) + 4); // would be for compact cluster (floor - forgiving)
		removeSmallClusters(
				diag_en,   // enable diagonal directions, false only up, dowm, right,left
				grown_by_1,    // selected tiles, will modified
				weights_src,   // or null
				grown_min,  // minimal number of pixels
				min_weight, // minimal total weight of the cluster
				min_max_weight, // minimal value of the maximal strengh in tghe cluster
				debugLevel);
		for (int i = 0; i< tiles.length; i++) tiles[i] &= grown_by_1[i]; 
	}

	public void removeSmallClusters(
			boolean    diag_en,   // enable diagonal directions, false only up, dowm, right,left
			boolean [] tiles_src,    // selected tiles, will modified
			double []  weights_src,   // or null
			int        min_area,  // minimal number of pixels
			double     min_weight, // minimal total weight of the cluster
			double     min_max_weight // minimal value of the maximal strengh in the cluster
			){
		removeSmallClusters(
				diag_en,   // enable diagonal directions, false only up, dowm, right,left
				tiles_src,    // selected tiles, will modified
				weights_src,   // or null
				min_area,  // minimal number of pixels
				min_weight, // minimal total weight of the cluster
				min_max_weight,
				0// minimal value of the maximal strengh in the cluster
				);
	}
	public void removeSmallClusters(
			boolean    diag_en,   // enable diagonal directions, false only up, dowm, right,left
			boolean [] tiles_src,    // selected tiles, will modified
			double []  weights_src,   // or null
			int        min_area,  // minimal number of pixels
			double     min_weight, // minimal total weight of the cluster
			double     min_max_weight, // minimal value of the maximal strengh in the cluster
			int        debugLevel
			){
		// adding 1-tile frame around to avoid checking for the borders
		int dbg_tile2 = 50095; // 217+326*153
				

		int tilesX2 = tilesX+2;
		int tilesY2 = tilesY+2;
		boolean [] tiles = new boolean [tilesX2*tilesY2];
		double [] weights = null;
		for (int i = 0; i < tilesY2; i++) {
			tiles[tilesX2 * 0 + i] =         false;
			tiles[tilesX2*(tilesY2-1) + i] = false;
			tiles[tilesX2* i + 0] =          false;
			tiles[tilesX2* i + tilesX2 -1] = false;
		}
		for (int i = 0; i < tilesY; i++) {
			for (int j = 0; j < tilesX; j++) {
				tiles[(i+1)*tilesX2 + j + 1]= tiles_src[i * tilesX + j];
			}
		}
		if (weights_src != null) {
			weights = new double [tiles.length];
			for (int i = 0; i < tilesY; i++) {
				for (int j = 0; j < tilesX; j++) {
					weights[(i+1)*tilesX2 + j + 1]= weights_src[i * tilesX + j];
				}
			}
		}


		int [] dirs4 = {-tilesX2, 1, tilesX2,-1};
		int [] dirs8 = {-tilesX2,  -tilesX2 + 1, 1, tilesX2 +1, tilesX2, tilesX2 - 1, -1, -tilesX2 - 1};
		int [] dirs = diag_en? dirs8 : dirs4;
		int [] waves = new int [tiles.length];
		for (int i = 0; i < tiles.length; i++) waves[i] = tiles[i] ? 0: -1;
		ArrayList<Integer> front = new ArrayList<Integer>();
		for (int start_indx = 0; start_indx < tiles.length; start_indx++) if (waves[start_indx] == 0){ // found first pixel of a new cluster
			if ((debugLevel > 1) && (start_indx == dbg_tile2)) {
				System.out.println("removeSmallClusters(): start_indx="+start_indx);
			}
			Integer ipx = start_indx;
			Integer ipx1;
			front.clear();
			int area = 1;
			double weight = 0.0;
			double max_weight = 0.0;
			if (weights != null) {
				weight += weights[ipx];
				if (weights[ipx] > max_weight) max_weight = weights[ipx];
			}
			waves[ipx] = area;
			front.add(ipx);
			while (!front.isEmpty()) {
				ipx = front.remove(0);// get oldest element
				for (int d = 0; d < dirs.length; d++){
					ipx1 = ipx + dirs[d];
					if (waves[ipx1] == 0) {
						area++;
						if (weights != null) {
							weight += weights[ipx1];
							if (weights[ipx1] > max_weight) max_weight = weights[ipx1];
						}
						waves[ipx1] = area;
						front.add(ipx1);
					}
				}
			}
			if ((area < min_area) || ((weights != null) && ((weight < min_weight) || (max_weight < min_max_weight)))){
				waves[ipx] = -1;
				tiles[ipx] = false;
				front.add(ipx);
				while (!front.isEmpty()) {
					ipx = front.remove(0);// get oldest element
					for (int d = 0; d < dirs.length; d++){
						ipx1 = ipx + dirs[d];
						if (waves[ipx1] > 0) {
							waves[ipx1] = -1;
							tiles[ipx1] = false;
							front.add(ipx1);
						}
					}
				}
			}
		}
		// return to original dimensions
		for (int i = 0; i < tilesY; i++) {
			for (int j = 0; j < tilesX; j++) {
				tiles_src[i * tilesX + j] =tiles[(i+1)*tilesX2 + j + 1];
			}
		}

	}


	public int [][] setSameTileOp(
			EyesisCorrectionParameters.CLTParameters           clt_parameters,
			int        op,
			int        debugLevel
			)
	{
		int [][] tile_op = new int [tilesY][tilesX]; // all zero
		int txl =  clt_parameters.tile_task_wl;
		int txr =  txl + clt_parameters.tile_task_ww;

		int tyt =  clt_parameters.tile_task_wt;
		int tyb =  tyt + clt_parameters.tile_task_wh;
		if      (txl < 0)       txl = 0;
		else if (txl >= tilesX) txl = tilesX - 1;

		if      (txr <= txl)    txr = txl + 1;
		else if (txr >  tilesX) txr = tilesX;

		if      (tyt < 0)       tyt = 0;
		else if (tyt >= tilesY) tyt = tilesY - 1;

		if      (tyb <= tyt)    tyb = tyt + 1;
		else if (tyb >  tilesY) tyb = tilesY;

		for (int i = tyt; i < tyb; i++) {
			for (int j = txl; j < txr; j++) {
				tile_op[i][j] = op;
			}
		}
		if (debugLevel > -1){
			System.out.println("clt_parameters.tile_task_wl="+clt_parameters.tile_task_wl );
			System.out.println("clt_parameters.tile_task_wt="+clt_parameters.tile_task_wt );
			System.out.println("clt_parameters.tile_task_ww="+clt_parameters.tile_task_ww );
			System.out.println("clt_parameters.tile_task_wh="+clt_parameters.tile_task_wh );
			System.out.println("getImgMask("+op+")="+ImageDtt.getImgMask(op) );
			System.out.println("getPairMask("+op+")="+ImageDtt.getPairMask(op) );
			System.out.println("getForcedDisparity("+op+")="+ImageDtt.getForcedDisparity(op) );
			System.out.println("getOrthoLines("+op+")="+ImageDtt.getOrthoLines(op) );
			System.out.println("tyt="+tyt+" tyb="+tyb+" txl="+txl+" txr="+txr );
		}
		return tile_op;
	}

	public double [][] setSameDisparity(
			double disparity)
	{
		double [][] disp_array = new double [tilesY][tilesX];
		for (int i = 0; i< tilesY; i++){
			for (int j = 0; j< tilesX; j++){
				disp_array[i][j] = disparity;
			}
		}
		return disp_array;
	}
	
	public boolean [] FilterScan(
			final CLTPass3d   scan,
			final boolean  [] bg_tiles,          // get from selected in clt_3d_passes.get(0);
			final int      [] horVertMod, // +1 - modified by hor correlation, +2 - modified by vert correlation (or null)
			final double      disparity_far,    //
			final double      disparity_near,   // 
			final double      this_sure,        // minimal strength to be considered definitely background
			final double      this_maybe,       // maximal strength to ignore as non-background
			final double      sure_smth,        // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
//			final EyesisCorrectionParameters.CLTParameters clt_parameters,
//			final int         threadsMax,  // maximal number of threads to launch
			
			final boolean     show_filter_scan,
			final int         min_clstr_seed, // clt_parameters.min_clstr_seed
			final double      min_clstr_weight,  // double     min_weight // minimal total weight of the cluster
			final double      min_clstr_max,     // double     min_max_weight // minimal value of the maximal strengh in the cluster
			final int         min_clstr_lone,    // int        min_area,  // minimal number of pixels
			final int         fill_gaps,         // int depth, // same as grow - odd - 4 directions, even - 8
			final int         poison_gaps,       // second (usually larger than fill_gaps) when any touch of prohibited cell invalidates patching
			final boolean     zero_gap_strength, // set strength to zero when covering gaps
			final int         debugLevel
	)
	{
/*
 					clt_parameters.min_clstr_weight,  // double     min_weight // minimal total weight of the cluster
					clt_parameters.min_clstr_max,    // double     min_max_weight // minimal value of the maximal strengh in the cluster
					clt_parameters.min_clstr_lone,  // int        min_area,  // minimal number of pixels
						clt_parameters.fill_gaps, // int depth, // same as grow - odd - 4 directions, even - 8
		
 */
		
		
		final int dbg_tile = 49468; // x = 220, y = 152 (pavement line)
		final boolean show_scan = show_filter_scan || (debugLevel > 1);
		showDoubleFloatArrays sdfa_instance = null;
		if ((debugLevel > -1) || show_scan) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		if (debugLevel > 0){
			System.out.println("FilterScan(,,"+disparity_far+", " +disparity_near+", "+ sure_smth);
		}
		final int tlen = tilesY * tilesX;
		double [] this_disparity     = scan.getDisparity(); // currently calculated, including ortho 
		double [] this_strength =      scan.getStrength(); // cloned, can be modified/ read back 
		double  [][] these_diffs =     scan.getDiffs();

//		double [] orig_strength =      scan.getOriginalStrength(); // to compare clusters 

		boolean [] these_tiles =       new boolean [tlen];
		boolean [] near_tiles =        new boolean [tlen];
		boolean [] far_tiles =         new boolean [tlen];
		boolean [] block_propagate =   new boolean [tlen];

		for (int i = 0; i <tlen; i++) if (!Double.isNaN(this_disparity[i])){
			if ((debugLevel > 0) && (i==dbg_tile)){
				System.out.println("FilterScan(): this_disparity["+i+"] = "+this_disparity[i]+"this_strength["+i+"] = "+this_strength[i]);
			}
			
			if (this_disparity[i] < disparity_far) {
				if (this_strength[i] > this_maybe){
					if (bg_tiles[i]) { // far can only be among previously selected for bgnd?
						far_tiles[i] = true;
					}
				}
			} else if (this_disparity[i] > disparity_near){
				//				  if ((this_strength[i] > this_maybe) && ! bg_tiles[i]){ // can not be farther if selected for near?
				//				  if ((this_strength[i] > this_maybe) || used_hor[i] || used_vert[i] ){ // can not be farther if selected for near?
				if ((this_strength[i] > this_maybe)){ // can not be farther if selected for near?
					near_tiles[i] = true;
				}
			} else { // in range
				if ((this_strength[i] > this_sure) || ((horVertMod != null) && (horVertMod[i] != 0))){
					these_tiles[i] = true;
				}
			}
			// see if the second worst variation exceeds sure_smth (like a window), really close object
			int imax1 = 0;
			for (int ip = 1; ip < these_diffs.length; ip++){
				if (these_diffs[ip][i] > these_diffs[imax1][i]) imax1 = ip;
			}
			int imax2 = (imax1 == 0)? 1 : 0;
			for (int ip = 0; ip< these_diffs.length; ip++) if (ip != imax1) {
				if (these_diffs[ip][i] > these_diffs[imax2][i]) imax2 = ip;
			}
			block_propagate[i] = (these_diffs[imax2][i] > sure_smth);
		}
		boolean[] prohibit = null; // TBD
		boolean[] dbg_before_small = null;
		boolean[] dbg_before_lone =  null;
		boolean[] dbg_before_gaps =  null;
		if (min_clstr_seed > 1){
			
			// TODO: check - now no limit on the strength of the offending selections, only on these onses
			
			removeSmallClusters(
					false, //true,            // boolean    diag_en,   // enable diagonal directions, false only up, dowm, right,left
					far_tiles,      // boolean [] tiles_src,    // selected tiles, will modified
					null,            // double []   weights_src,   // or null
					min_clstr_seed,  // int        min_area,  // minimal number of pixels
					0.0,  // clt_parameters.min_clstr_weight,  // double     min_weight // minimal total weight of the cluster
					0.0); // clt_parameters.min_clstr_max);    // double     min_max_weight // minimal value of the maximal strengh in the cluster
			removeSmallClusters(
					false, // true,            // boolean    diag_en,   // enable diagonal directions, false only up, dowm, right,left
					near_tiles,      // boolean [] tiles_src,    // selected tiles, will modified
					null,            // double []   weights_src,   // or null
					min_clstr_seed,  // int        min_area,  // minimal number of pixels
					0.0, // clt_parameters.min_clstr_weight,  // double     min_weight // minimal total weight of the cluster
					0.0); // clt_parameters.min_clstr_max);    // double     min_max_weight // minimal value of the maximal strengh in the cluster
			if ((sdfa_instance!=null) && show_scan) dbg_before_small = these_tiles.clone();
			// only remove far outstanding clusters
			removeSmallClusters( // remove single-tile clusters - anywhere
					false, // true,            // boolean    diag_en,   // enable diagonal directions, false only up, dowm, right,left
					these_tiles,     // boolean [] tiles_src,    // selected tiles, will modified
					this_strength, // orig_strength, // null,            // double []   weights_src,   // or null
					min_clstr_seed, // 2,               // int        min_area,  // minimal number of pixels
					min_clstr_weight,  // double     min_weight // minimal total weight of the cluster
					min_clstr_max,    // double     min_max_weight // minimal value of the maximal strengh in the cluster
					debugLevel);
			if ((sdfa_instance!=null) && show_scan) dbg_before_lone = these_tiles.clone();
			removeLoneClusters(
					false, // true,            // boolean    diag_en,   // enable diagonal directions, false only up, dowm, right,left
					these_tiles,      // boolean [] tiles_src,    // selected tiles, will modified
					this_strength, // orig_strength, // null,            // double []   weights_src,   // or null
					min_clstr_lone,  // int        min_area,  // minimal number of pixels
					min_clstr_weight,  // double     min_weight // minimal total weight of the cluster
					min_clstr_max,    // double     min_max_weight // minimal value of the maximal strengh in the cluster
					debugLevel);

			if ((sdfa_instance!=null) && show_scan) dbg_before_gaps = these_tiles.clone();
			boolean [] pre_gap = null;
			if ((fill_gaps > 0) || (poison_gaps > 0)) {
				if (zero_gap_strength) pre_gap = these_tiles.clone();
				prohibit = far_tiles.clone(); // do not fill gaps over known background/far tiles
				if (fill_gaps > 0) {
					fillGaps( // grows, then shrinks
							fill_gaps,    // int depth, // same as grow - odd - 4 directions, even - 8
							false,         // boolean    poison, // do not fill gaps that even touch prohibited
							these_tiles,               // boolean [] tiles,
							prohibit);
				}

				if (poison_gaps > 0) { // only block for when poison_gaps is enabled (last filtering with high fill_gap value)
					for (int i = 0; i< prohibit.length; i++) prohibit[i] |=  block_propagate[i];
					fillGaps( // grows, then shrinks
							poison_gaps, // int depth, // same as grow - odd - 4 directions, even - 8
							true, // boolean    poison, // do not fill gaps that even touch prohibited
							these_tiles,               // boolean [] tiles,
							prohibit);
				}
				if (zero_gap_strength) {
					for (int i = 0; i < these_tiles.length; i++) if (these_tiles[i] && !pre_gap[i]){
						this_strength[i] = 0.0;
					}
				}
			}
			//zero_gap_strength
		}		  
		
		double [] this_disparity_masked = this_disparity.clone();
		for (int i = 0; i < this_disparity.length; i++){
			if (!these_tiles[i])this_disparity_masked[i] = Double.NaN;
		}
		
		
		if ((sdfa_instance!=null) && show_scan){
			
			int [] enum_clusters = enumerateClusters(
					true, // boolean diag_en,
					these_tiles); // boolean [] tiles_src)
			
			String [] titles = {"masked","map","orig_map","hor_map","vert_map","bg_sel","far",
					"before_small","before_lone","before_gaps","these","near","block",
					"strength","hor-strength","vert-strength",
					"diff0","diff1","diff2","diff3", "enum_clusters", "disp_cm", "disp_poly", "disp_hor", "disp_vert", "poison"};
			double [][] dbg_img = new double[titles.length][tilesY * tilesX];
			for (int i = 0; i<dbg_img[0].length;i++){
				dbg_img[ 0][i] =   this_disparity_masked[i];
				dbg_img[ 1][i] =   this_disparity[i];
//				dbg_img[ 2][i] =   dbg_orig_disparity[i];
//				dbg_img[ 3][i] =   this_hor_disparity[i];
//				dbg_img[ 4][i] =   this_vert_disparity[i];
				dbg_img[ 5][i] =    bg_tiles    [i] ? 1 : -1;
				dbg_img[ 6][i] =     far_tiles  [i] ? 1 : -1;
				dbg_img[ 7][i] =   dbg_before_small  [i] ? 1 : -1;
				dbg_img[ 8][i] =   dbg_before_lone  [i] ? 1 : -1;
				
				
				dbg_img[ 9][i] =   dbg_before_gaps  [i] ? 1 : -1;
				dbg_img[10][i] =   these_tiles  [i] ? 1 : -1;
				dbg_img[11][i] =    near_tiles  [i] ? 1 : -1;
				dbg_img[12][i] = block_propagate[i] ? 1 : -1;
				dbg_img[13][i] =   this_strength[i];
//				dbg_img[14][i] =   hor_strength[i];
//				dbg_img[15][i] =   vert_strength[i];
				
				dbg_img[16][i] =   these_diffs[0][i];
				dbg_img[17][i] =   these_diffs[1][i];
				dbg_img[18][i] =   these_diffs[2][i];
				dbg_img[19][i] =   these_diffs[3][i];
				dbg_img[20][i] =   enum_clusters[i];
				dbg_img[25][i] =   prohibit[i] ? 1 : -1;
			}
			dbg_img[ 2] =   scan.getDisparity(1);
			dbg_img[ 3] =   scan.getDisparity(2);
			dbg_img[ 4] =   scan.getDisparity(3);
			dbg_img[14] =   scan.getHorStrength();
			dbg_img[15] =   scan.getVertStrength();
			
			dbg_img[21] =  scan. disparity_map[ImageDtt.DISPARITY_INDEX_CM];
			dbg_img[22] =  scan. disparity_map[ImageDtt.DISPARITY_INDEX_POLY];
			dbg_img[23] =  scan. disparity_map[ImageDtt.DISPARITY_INDEX_HOR];
			dbg_img[24] =  scan. disparity_map[ImageDtt.DISPARITY_INDEX_VERT];
			
			sdfa_instance.showArrays(dbg_img,  tilesX, tilesY, true, "FilterScan"+clt_3d_passes.size(),titles);
			System.out.println("FilterScan"+clt_3d_passes.size());
		}
		
		return these_tiles;
	}
	
	
	public int [] combineOrthoDisparity(
			final CLTPass3d   scan,        // scan data
			final boolean     or_hor,       // true;  // Apply ortho correction to horizontal correlation (vertical features)
			final boolean     or_vert,      // true;  // Apply ortho correction to vertical correlation (horizontal features)
			final double      or_sigma,     // 2.0;   // Blur sigma: verically for horizontal correlation, horizontally - for vertically
			final double      or_sharp,     // 0.5;   // 3-point sharpening (-k, +2k+1, -k)
			final double      or_scale,     // 2.0;   // Scale ortho correletion strength relative to 4-directional one
			final double      or_offset,    // 0.1;   // Subtract from scaled correlation strength, limit by 0
			final double      or_asym ,     // 1.5;   // Minimal ratio of orthogonal strengths required for dis[parity replacement
			final double      or_threshold, // 1.5;   // Minimal scaled offsetg ortho strength to normal strength needed for replacement
	  		final double      or_absHor,    // 0.15;  // Minimal horizontal absolute scaled offset ortho strength needed for replacement
	  		final double      or_absVert,   // 0.19;  // Minimal vertical absolute scaled offset ortho strength needed for replacement
			final int         debugLevel)
	{
		double [] disparity =      scan.getDisparity(0);  // calculated, to be modified
		double [] strength =       scan.getStrength();    // combo, to be modified
		double [] disparity_hor =  scan.getDisparity(2); // .clone(); // calculated
		double [] disparity_vert = scan.getDisparity(3); // .clone(); // calculated
		double [] strength_hor =   scan.getHorStrength(); // .clone();
		double [] strength_vert =  scan.getVertStrength(); // .clone();
		int [] replaced = new int[strength.length];
		if (or_hor){
			for (int i = 0; i < strength_hor.length; i++){
				strength_hor[i] *= or_scale;
			}
			SharpBlurPair(
					disparity_hor,  // double [] data,     // data array for in-place modification
					strength_hor,   // double [] strength, // data weights array for in-place modification
					or_sigma,       // double sigma,       // blur sigma
					or_sharp,       // double k,           // sharpen in orthogonal direction with (-k,2*k-1,-k). 0 - no sharpening
					or_offset,      // double offset,      // subtract from strength, limit by 0.0
					false);         // boolean vert)       // true - sharpen vertically,  blur horizontally. False - sharpen horizontally, blur vertically
		}
		if (or_vert){
			for (int i = 0; i < strength_vert.length; i++){
				strength_vert[i] *= or_scale;
			}
			SharpBlurPair(
					disparity_vert, // double [] data,     // data array for in-place modification
					strength_vert,  // double [] strength, // data weights array for in-place modification
					or_sigma,       // double sigma,       // blur sigma
					or_sharp,       // double k,           // sharpen in orthogonal direction with (-k,2*k-1,-k). 0 - no sharpening
					or_offset,      // double offset,      // subtract from strength, limit by 0.0
					true);          // boolean vert)       // true - sharpen vertically,  blur horizontally. False - sharpen horizontally, blur vertically
		}
		double ko = (or_threshold - 1) * or_offset;
		double ao = (or_asym - 1) * or_offset;
		if (or_hor){
			for (int i = 0; i < strength_hor.length; i++){
				if (
						(strength_hor[i] > or_absHor) && 
						(strength_hor[i] > or_threshold * strength[i] - ko) &&
						(!or_vert || (strength_hor[i] > or_asym * strength_vert[i] - ao))){
					strength[i] =  strength_hor[i];
					disparity[i] = disparity_hor[i];
					replaced[i] |= 1;
				}
			}
		}
		if (or_vert){
			for (int i = 0; i < strength_vert.length; i++){
				if (
						(strength_vert[i] > or_absVert) &&
						(strength_vert[i] > or_threshold * strength[i] - ko) &&
						(!or_hor || (strength_vert[i] > or_asym * strength_hor[i] - ao))){
					strength[i] =  strength_vert[i];
					disparity[i] = disparity_vert[i];
					replaced[i] |= 2;
				}
			}
		}
		return replaced;
	}
	
	public boolean [] combineHorVertDisparity(
			final CLTPass3d   scan,
			final boolean [] bg_tiles,          // get from selected in clt_3d_passes.get(0); 
			final double      disparity_far,    //
			final double      disparity_near,   // 
			final double      this_sure,        // minimal strength to be considered definitely background
			final double      this_maybe,       // maximal strength to ignore as non-background
			final double      sure_smth,        // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
			final EyesisCorrectionParameters.CLTParameters clt_parameters,
//			final int         threadsMax,  // maximal number of threads to launch                         
//			final boolean     updateStatus,
			final int         debugLevel
	)
	{
		showDoubleFloatArrays sdfa_instance = null;
		if (debugLevel > -1) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		
		// scale and shift all 3 disparities - combo, vert and hor
		final int tlen = tilesY * tilesX;
		double [] this_disparity     = scan.getDisparity(); // returns a copy of the FPGA-generated disparity combined with the target one
		double [] this_hor_disparity=  scan.getDisparity(2);
		double [] this_vert_disparity= scan.getDisparity(3);
		double [] this_strength =      scan.getStrength(); // cloned, can be modified/ read back 
		double [] strength =           scan.getOriginalStrength(); // reference, not cloned 
		double [] hor_strength =       scan.getHorStrength();
		double [] vert_strength =      scan.getVertStrength();
		double [] dbg_orig_disparity = scan.getDisparity(1); // unmodified
		double  [][] these_diffs =     scan.getDiffs();

		boolean [] these_tiles =       new boolean [tlen];
		boolean [] near_tiles =        new boolean [tlen];
		boolean [] far_tiles =         new boolean [tlen];
		
		boolean [] block_propagate =   new boolean [tlen];
		boolean [] used_hor =          new boolean [tlen];
		boolean [] used_vert =         new boolean [tlen];

		// convolve hor/vert strengths to emphasize multi-tile features
		double []  hor_strength_conv = detectOrtho(
				hor_strength, // double [] tiles,
				false, // boolean vert,    // true for vertical correlation and _horizontal features
				clt_parameters.ortho_half_length, // int radius,      // one dimension - [-radius, +radius], other dimension -1,0, 1
				debugLevel);
		double []  vert_strength_conv = detectOrtho(
				vert_strength, // double [] tiles,
				true, // boolean vert,    // true for vertical correlation and _horizontal features
				clt_parameters.ortho_half_length, // int radius,      // one dimension - [-radius, +radius], other dimension -1,0, 1
				debugLevel);

		
		
		
		
		
		
		// now mix raw with convolved, but keep raw - it will be needed for weighted average
		for (int i = 0; i< tlen; i++){
			hor_strength_conv[i] =  (clt_parameters.ortho_mix * hor_strength_conv[i] + (1.0 - clt_parameters.ortho_mix ) * hor_strength[i]);
			vert_strength_conv[i] = (clt_parameters.ortho_mix * vert_strength_conv[i] + (1.0 - clt_parameters.ortho_mix ) * vert_strength[i]);
		}

		//getOriginalStrength()

		for (int i = 0; i < tilesY; i++){
			for (int j = 0; j < tilesX; j++){
				int indx = i * tilesX + j;
				if (!Double.isNaN(this_disparity[i])){
					double disp = this_disparity[indx];
					// Enhance foreground detection: compare horizontal-only (for vertical features) and vertical (for horizontal) and replace 4-pair disparity
					// with that value
					if ((hor_strength_conv[indx] >= clt_parameters.ortho_min_hor) &&
							(hor_strength_conv[indx] / vert_strength_conv[indx]>= clt_parameters.ortho_asym) &&
							(this_hor_disparity[indx] > disp)) {
						disp = this_hor_disparity[indx];
						used_hor[indx] = true;
					}
					if ((vert_strength_conv[indx] >= clt_parameters.ortho_min_vert) &&
							(vert_strength_conv[indx] / hor_strength_conv[indx]>= clt_parameters.ortho_asym) &&
							(this_vert_disparity[indx] > disp)) {
						disp = this_vert_disparity[indx];
						used_vert[indx] = true;
					}
				}
			}
		}
		
		boolean [] dbg_used_hor =  used_hor.clone();
		boolean [] dbg_used_vert = used_vert.clone();
		
		if (clt_parameters.min_clstr_seed > 1){
			// TODO: do we need to limit the strenghs of the clusters?
			removeSmallClusters(
					false, //true,                 // boolean    diag_en,   // enable diagonal directions, false only up, dowm, right,left
					used_hor,                      // boolean [] tiles_src,    // selected tiles, will modified
					null,                          // double []   weights_src,   // or null
					clt_parameters.min_clstr_seed, // int        min_area,  // minimal number of pixels
					0.0, //clt_parameters.min_clstr_weight,  // double     min_weight // minimal total weight of the cluster
					0.0); // clt_parameters.min_clstr_max);    // double     min_max_weight // minimal value of the maximal strengh in the cluster
			
			removeSmallClusters(
					false,  // true,               // boolean    diag_en,   // enable diagonal directions, false only up, dowm, right,left
					used_vert,                     // boolean [] tiles_src,    // selected tiles, will modified
					null,                          // double []   weights_src,   // or null
					clt_parameters.min_clstr_seed, // int        min_area,  // minimal number of pixels
					0.0, //clt_parameters.min_clstr_weight,  // double     min_weight // minimal total weight of the cluster
					0.0); // clt_parameters.min_clstr_max);    // double     min_max_weight // minimal value of the maximal strengh in the cluster
		}		  

		// bridge over small gaps in horizontal/vertical features
		int numHorBridged = bridgeFgndOrthoGap(
				clt_parameters,
				false, // vert, // verical pairs, horizontal features
				true,  // boolean disp_interpolate, // fill interpolated disparity for bridged over tiles, false - copy from  ortho_disparity (if not null)
				true,  // boolean closer_only, // only update disparity if laregr than was 
				used_hor, // boolean [] used_ortho,
				strength,
				hor_strength_conv, //  for comparing
				hor_strength, //  for weighted average
				this_disparity, // will be modified
				this_hor_disparity // may be null
				);
		if (debugLevel > -1) System.out.println("Bridged over "+numHorBridged+" tiles along vertical features");
		// bridge over small gaps in horizontal/vertical features
		
		int numVertBridged = bridgeFgndOrthoGap(
				clt_parameters,
				true, // vert, // verical pairs, horizontal features
				false, // true,  // boolean disp_interpolate, // fill interpolated disparity for bridged over tiles, false - copy from  ortho_disparity (if not null)
				true,  // boolean closer_only, // only update disparity if laregr than was 
				used_vert, // boolean [] used_ortho,
				strength,
				vert_strength_conv, //  double [] ortho_strength,
				vert_strength, //  double [] ortho_strength,
				this_disparity, // will be modified
				this_vert_disparity // may be null
				);
		if (debugLevel > -1) System.out.println("Bridged over "+numVertBridged+" tiles along horizontal features");

		
		for (int i = 0; i <tlen; i++) if (!Double.isNaN(this_disparity[i])){
			if (this_disparity[i] < disparity_far) {
				if (this_strength[i] > this_maybe){
					if (bg_tiles[i]) { // far can only be among previously selected for bgnd?
						far_tiles[i] = true;
					}
				}
			} else if (this_disparity[i] > disparity_near){
				//				  if ((this_strength[i] > this_maybe) && ! bg_tiles[i]){ // can not be farther if selected for near?
				//				  if ((this_strength[i] > this_maybe) || used_hor[i] || used_vert[i] ){ // can not be farther if selected for near?
				if ((this_strength[i] > this_maybe)){ // can not be farther if selected for near?
					near_tiles[i] = true;
				}
			} else { // in range
				if ((this_strength[i] > this_sure) || used_hor[i] || used_vert[i])
					these_tiles[i] = true;
			}
			// see if the second worst variation exceeds sure_smth (like a window), really close object
			int imax1 = 0;
			for (int ip = 1; ip < these_diffs.length; ip++){
				if (these_diffs[ip][i] > these_diffs[imax1][i]) imax1 = ip;
			}
			int imax2 = (imax1 == 0)? 1 : 0;
			for (int ip = 0; ip< these_diffs.length; ip++) if (ip != imax1) {
				if (these_diffs[ip][i] > these_diffs[imax2][i]) imax2 = ip;
			}
			block_propagate[i] = (these_diffs[imax2][i] > sure_smth);
		}
		boolean[] prohibit = null; // TBD
		boolean[] dbg_before_gaps = null;
		if (clt_parameters.min_clstr_seed > 1){

			// TODO: check - now no limit on the strength of the offending selections, only on these onses
			
			removeSmallClusters(
					false, //true,            // boolean    diag_en,   // enable diagonal directions, false only up, dowm, right,left
					far_tiles,      // boolean [] tiles_src,    // selected tiles, will modified
					null,            // double []   weights_src,   // or null
					clt_parameters.min_clstr_seed,  // int        min_area,  // minimal number of pixels
					0.0, //clt_parameters.min_clstr_weight,  // double     min_weight // minimal total weight of the cluster
					0.0); // clt_parameters.min_clstr_max);    // double     min_max_weight // minimal value of the maximal strengh in the cluster
			removeSmallClusters(
					false, // true,            // boolean    diag_en,   // enable diagonal directions, false only up, dowm, right,left
					near_tiles,      // boolean [] tiles_src,    // selected tiles, will modified
					null,            // double []   weights_src,   // or null
					clt_parameters.min_clstr_seed,  // int        min_area,  // minimal number of pixels
					0.0, //clt_parameters.min_clstr_weight,  // double     min_weight // minimal total weight of the cluster
					0.0); // clt_parameters.min_clstr_max);    // double     min_max_weight // minimal value of the maximal strengh in the cluster
			// only remove far outstanding clusters
			removeSmallClusters( // remove single-tile clusters - anywhere
					false, // true,            // boolean    diag_en,   // enable diagonal directions, false only up, dowm, right,left
					these_tiles,     // boolean [] tiles_src,    // selected tiles, will modified
					null,            // double []   weights_src,   // or null
					clt_parameters.min_clstr_seed, // 2,               // int        min_area,  // minimal number of pixels
					clt_parameters.min_clstr_weight,  // double     min_weight // minimal total weight of the cluster
					clt_parameters.min_clstr_max);    // double     min_max_weight // minimal value of the maximal strengh in the cluster

			removeLoneClusters(
					false, // true,            // boolean    diag_en,   // enable diagonal directions, false only up, dowm, right,left
					these_tiles,      // boolean [] tiles_src,    // selected tiles, will modified
					null,            // double []   weights_src,   // or null
					clt_parameters.min_clstr_lone,  // int        min_area,  // minimal number of pixels
					clt_parameters.min_clstr_weight,  // double     min_weight // minimal total weight of the cluster
					clt_parameters.min_clstr_max);    // double     min_max_weight // minimal value of the maximal strengh in the cluster
			dbg_before_gaps = these_tiles.clone();
			prohibit = far_tiles.clone(); // do not fill gaps over known background/far tiles
			for (int i = 0; i< prohibit.length; i++) prohibit[i] |=  block_propagate[i];
			
			if (clt_parameters.fill_gaps > 0) {
				fillGaps( // grows, then shrinks
						clt_parameters.fill_gaps, // int depth, // same as grow - odd - 4 directions, even - 8
						false, // boolean    poison, // do not fill gaps that even touch prohibited
						these_tiles,               // boolean [] tiles,
						prohibit);
			}
		}		  
		
		
		double [] this_disparity_masked = this_disparity.clone();
		for (int i = 0; i < this_disparity.length; i++){
			if (!these_tiles[i])this_disparity_masked[i] = Double.NaN;
		}

		if ((sdfa_instance!=null) && clt_parameters.show_bgnd_nonbgnd) {
			
			int [] enum_clusters = enumerateClusters(
					true, // boolean diag_en,
					these_tiles); // boolean [] tiles_src)
			
			String [] titles = {"masked","map","orig_map","hor_map","vert_map","bg_sel","far","these_gaps","these","near","block",
					"strength","hor-strength","hor-conv-strength","vert-strength","vert-conv-strength",
					"hor","hor-bridged","vert","vert-bridged","diff0","diff1","diff2","diff3", "enum_clusters", "disp_cm", "disp_poly", "disp_hor", "disp_vert"};
			double [][] dbg_img = new double[titles.length][tilesY * tilesX];
			for (int i = 0; i<dbg_img[0].length;i++){
				dbg_img[ 0][i] =   this_disparity_masked[i];
				dbg_img[ 1][i] =   this_disparity[i];
				dbg_img[ 2][i] =   dbg_orig_disparity[i];
				dbg_img[ 3][i] =   this_hor_disparity[i];
				dbg_img[ 4][i] =   this_vert_disparity[i];
				dbg_img[ 5][i] =    bg_tiles    [i] ? 1 : -1;
				dbg_img[ 6][i] =     far_tiles  [i] ? 1 : -1;
				dbg_img[ 7][i] =   dbg_before_gaps  [i] ? 1 : -1;
				dbg_img[ 8][i] =   these_tiles  [i] ? 1 : -1;
				dbg_img[ 9][i] =    near_tiles  [i] ? 1 : -1;
				dbg_img[10][i] = block_propagate[i] ? 1 : -1;
				dbg_img[11][i] =   this_strength[i];
				dbg_img[12][i] =   hor_strength[i];
				dbg_img[13][i] =   hor_strength_conv[i];
				dbg_img[14][i] =   vert_strength[i];
				dbg_img[15][i] =   vert_strength_conv[i];
				dbg_img[16][i] =   dbg_used_hor[i]?  1 : -1;
				dbg_img[17][i] =   used_hor[i]?      1 : -1;
				dbg_img[18][i] =   dbg_used_vert[i]? 1 : -1;
				dbg_img[19][i] =   used_vert[i]?     1 : -1;
				
				dbg_img[20][i] =   these_diffs[0][i];
				dbg_img[21][i] =   these_diffs[1][i];
				dbg_img[22][i] =   these_diffs[2][i];
				dbg_img[23][i] =   these_diffs[3][i];
				dbg_img[24][i] =   enum_clusters[i];
			}
			dbg_img[25] =  scan. disparity_map[ImageDtt.DISPARITY_INDEX_CM];
			dbg_img[26] =  scan. disparity_map[ImageDtt.DISPARITY_INDEX_POLY];
			dbg_img[27] =  scan. disparity_map[ImageDtt.DISPARITY_INDEX_HOR];
			dbg_img[28] =  scan. disparity_map[ImageDtt.DISPARITY_INDEX_VERT];
			
			sdfa_instance.showArrays(dbg_img,  tilesX, tilesY, true, "bgnd_nonbgnd"+clt_3d_passes.size(),titles);
		}
		
		return these_tiles;
		
	}
	

	public CLTPass3d refinePassSetup( // prepare tile tasks for the second pass based on the previous one(s)
			//			  final double [][][]       image_data, // first index - number of image in a quad
			EyesisCorrectionParameters.CLTParameters           clt_parameters,
			// disparity range - differences from
			boolean           use_supertiles,
			int               bg_scan_index,
			double            disparity_far,    //
			double            disparity_near,   // 
			double            this_sure,        // minimal strength to be considered definitely background
			double            this_maybe,       // maximal strength to ignore as non-background
			double            sure_smth,        // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
			int               disparity_index,  // index of disparity value in disparity_map == 2 (0,2 or 4)
			GeometryCorrection geometryCorrection,
			final int         threadsMax,  // maximal number of threads to launch                         
			final boolean     updateStatus,
			final int         debugLevel)
	{
		CLTPass3d scan_prev = clt_3d_passes.get(clt_3d_passes.size() -1);
		CLTPass3d scan_bg =   clt_3d_passes.get(bg_scan_index); // 
		showDoubleFloatArrays sdfa_instance = null;
		if (debugLevel > -1) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		//TODO: for next passes - combine all selected for previous passes (all passes with smaller disparity)
		int [] replaced =  null; // +1 - hor, +2 - vert
		int [] replaced0 = null; // +1 - hor, +2 - vert
//		if (clt_parameters.or_hor || clt_parameters.or_vert) {
			// TODO:			add filtering before/after
		double [][] dbg_img = null;
		String [] dbg_titles = null;
		boolean show_ortho = clt_parameters.show_ortho_combine || (debugLevel > 1);
		boolean show_super = clt_parameters.show_refine_supertiles || (debugLevel > 1);
		boolean show_st =    clt_parameters.stShow || (debugLevel > 1);
		if (show_ortho){
			String [] dbg_titles0 = {
					"combo_disparity",    //  0
					"orig_disparity",     //  1
					"hor_disparity",      //  2
					"hor_orig_disparity", //  3
					"vert_disparity",     //  4
					"vert_orig_disparity",//  5
					"combo_strength",     //  6
					"orig_strength",      //  7
					"hor_strength",       //  8
					"hor_orig_strength",  //  9
					"vert_strength",      // 10
					"vert_orig_strength", // 11
					"replaced0",          // 12
					"replaced",           // 13
					"selection",          // 14  
			        "tilesHor"};          // 15
			dbg_titles = dbg_titles0;
			dbg_img = new double [dbg_titles.length][];
			dbg_img[ 1] = scan_prev.getDisparity(1).clone();
			dbg_img[ 3] = scan_prev.getDisparity(2).clone();
			dbg_img[ 5] = scan_prev.getDisparity(3).clone();
			dbg_img[ 7] = scan_prev.getStrength().clone();
			dbg_img[ 9] = scan_prev.getHorStrength().clone();
			dbg_img[11] = scan_prev.getVertStrength().clone();
			dbg_img[14] = new double [scan_prev.getDisparity().length];
			dbg_img[15] = new double [scan_prev.getDisparity().length];
		}			
		replaced = combineOrthoDisparity(
				scan_prev, // final CLTPass3d   scan,        // scan data
				clt_parameters.or_hor,       // true;  // Apply ortho correction to horizontal correlation (vertical features)
				clt_parameters.or_vert,      // true;  // Apply ortho correction to vertical correlation (horizontal features)
				clt_parameters.or_sigma,     // 2.0;   // Blur sigma: verically for horizontal correlation, horizontally - for vertically
				clt_parameters.or_sharp,     // 0.5;   // 3-point sharpening (-k, +2k+1, -k)
				clt_parameters.or_scale,     // 2.0;   // Scale ortho correletion strength relative to 4-directional one
				clt_parameters.or_offset,    // 0.1;   // Subtract from scaled correlation strength, limit by 0
				clt_parameters.or_asym ,     // 1.5;   // Minimal ratio of orthogonal strengths required for dis[parity replacement
				clt_parameters.or_threshold, // 1.5;   // Minimal scaled offsetg ortho strength to normal strength needed for replacement
				clt_parameters.or_absHor,    // 0.15;  // Minimal horizontal absolute scaled offset ortho strength needed for replacement
				clt_parameters.or_absVert,   // 0.19;  // Minimal vertical absolute scaled offset ortho strength needed for replacement
				debugLevel);

		if (clt_parameters.poles_fix) {
			boolean [] selection = new boolean [replaced.length];
			boolean [] tilesHor =   new boolean [replaced.length];
			double  [] disparity = scan_prev.getDisparity();
			for (int i = 0; i < tilesHor.length; i++){
				tilesHor[i] = (replaced[i] & 1) != 0;
				selection[i] = !Double.isNaN(disparity[i]) && (disparity[i] >= disparity_far) && (disparity[i] <= disparity_near);
			}
			if (show_ortho){
				for (int i = 0; i < tilesHor.length; i++){
					dbg_img[14][i] = selection[i]?1.0:0.0;
					dbg_img[15][i] = tilesHor[i]?1.0:0.0;
				}
			}
			int numFixed = fixVerticalPoles( // return number of replaced cells
					scan_prev, // CLTPass3d   scan,             // scan data to use
					selection,                         // start with only from selections (if not null, continue regardless)
					tilesHor,                          // horizontal correlation tiles used for composite disparity/strength;
					clt_parameters.poles_len ,         // int         max_len,          // maximal length to cover
					clt_parameters.poles_ratio,        // Maximal ratio of invisible to visible pole length
					clt_parameters.poles_min_strength, // double      min_new_strength, // set strength to hor_strength, but not less than this
					clt_parameters.poles_force_disp,   // boolean     force_disparity   // copy disparity down (false - use horDisparity
					true); 
			if (debugLevel > -1){
				System.out.println("fixVerticalPoles() replaced "+ numFixed+ " tiles.");
			}
			replaced0 = replaced.clone();
			for (int i = 0; i < replaced.length; i++){
				if (tilesHor[i]) replaced[i] |= 1;
			}
		}

		if (show_ortho){

			dbg_img[ 0] = scan_prev.getDisparity(0);
			dbg_img[ 2] = scan_prev.getDisparity(2);
			dbg_img[ 4] = scan_prev.getDisparity(3);
			dbg_img[ 6] = scan_prev.getStrength();
			dbg_img[ 8] = scan_prev.getHorStrength();
			dbg_img[10] = scan_prev.getVertStrength();

			double [] dreplaced0 = new double [replaced.length];
			double [] dreplaced =  new double [replaced.length];
			for (int i = 0; i < dreplaced.length; i++){
				dreplaced0[i] = replaced0[i];
				dreplaced[i] =  replaced[i];
			}
			dbg_img[12] = dreplaced0;
			dbg_img[13] = dreplaced;
			sdfa_instance.showArrays(dbg_img, tilesX, tilesY, true, "ortho_combine",dbg_titles);
		}

		boolean [] these_tiles = FilterScan(
				scan_prev,        // final CLTPass3d   scan,
				scan_bg.selected, // get from selected in clt_3d_passes.get(0);
				replaced,         // final int      [] horVertMod, // +1 - modified by hor correlation, +2 - modified by vert correlation (or null)
				disparity_far,    // final double      disparity_far,    //
				disparity_near,   // final double      disparity_near,   // 
				this_sure,        // final double      this_sure,        // minimal strength to be considered definitely background
				this_maybe,       // final double      this_maybe,       // maximal strength to ignore as non-background
				sure_smth,        // final double      sure_smth,        // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
				clt_parameters.show_filter_scan,
				clt_parameters.min_clstr_seed, // clt_parameters.min_clstr_seed
				clt_parameters.min_clstr_weight,  // double     min_weight // minimal total weight of the cluster
				clt_parameters.min_clstr_max,    // double     min_max_weight // minimal value of the maximal strengh in the cluster
				clt_parameters.min_clstr_lone,  // int        min_area,  // minimal number of pixels
				clt_parameters.fill_gaps, // int depth, // same as grow - odd - 4 directions, even - 8
				0, // final int     poison_gaps,      // Do not fill gaps that have even single "poisoned" tile 
				false,             // final boolean     zero_gap_strength, // set strength to zero when covering gaps
				debugLevel);


		//		}

		/*		

		boolean [] these_tiles = combineHorVertDisparity(
				scan_prev,                     // final CLTPass3d   scan,
				scan_prev.selected, // clt_3d_passes.get(0).selected, // final boolean [] bg_tiles,          // get from selected in clt_3d_passes.get(0); 
				disparity_far,    //
				disparity_near,   // 
				this_sure,        // minimal strength to be considered definitely background
				this_maybe,       // maximal strength to ignore as non-background
				sure_smth,        // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
				clt_parameters,
				debugLevel);

		scan_prev.combineHorVertStrength(true, false); // strength now max of original and horizontal. Use scale  instead of boolean?
		 */



		//		double [] this_disparity     = scan_prev.getDisparity(); // returns a copy of the FPGA-generated disparity combined with the target one
		//		double [] this_strength =      scan_prev.getStrength(); // cloned, can be modified/ read back 
		//************************************************
		// Show supertiles histograms

		//		if (clt_parameters.stShow){

		// try renovated supertiles. Do twice to show both original and blured histograms
		double [] dbg_orig_disparity =  null;
		double [] dbg_with_super_disp = null;
		double [] dbg_outlayers =       null;

		if (use_supertiles) {	
			String [] dbg_st_titles = {"raw", "blurred"+clt_parameters.stSigma,"max-min-max"}; 
			double [][] dbg_hist = new double[dbg_st_titles.length][];

			scan_prev.setSuperTiles(
					clt_parameters.stStepNear,       // double     step_disparity,
					clt_parameters.stStepFar,        // double     step_near,
					clt_parameters.stStepThreshold,  // double     step_threshold,
					clt_parameters.stMinDisparity,   // double     min_disparity,
					clt_parameters.stMaxDisparity,   // double     max_disparity,
					clt_parameters.stFloor,          // double     strength_floor,
					clt_parameters.stPow,            // double     strength_pow,
					0.0,// NO BLUR double     stBlurSigma)
					false, //clt_parameters.stSmplMode,  // Use sample mode (false - regular tile mode)
					clt_parameters.stSmplSide,  // Sample size (side of a square)
					clt_parameters.stSmplNum,   // Number after removing worst
					clt_parameters.stSmplRms,   // Maximal RMS of the remaining tiles in a sample
					clt_parameters.stMeasSel); // bitmask of the selected measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			dbg_hist[0] = scan_prev.showDisparityHistogram();

			scan_prev.setSuperTiles(
					clt_parameters.stStepNear,       // double     step_disparity,
					clt_parameters.stStepFar,        // double     step_near,
					clt_parameters.stStepThreshold,  // double     step_threshold,
					clt_parameters.stMinDisparity,   // double     min_disparity,
					clt_parameters.stMaxDisparity,   // double     max_disparity,
					clt_parameters.stFloor,          // double     strength_floor,
					clt_parameters.stPow,            // double     strength_pow,
					clt_parameters.stSigma, // with blur double     stBlurSigma)
					false, //clt_parameters.stSmplMode,  // Use sample mode (false - regular tile mode)
					clt_parameters.stSmplSide,  // Sample size (side of a square)
					clt_parameters.stSmplNum,   // Number after removing worst
					clt_parameters.stSmplRms,   // Maximal RMS of the remaining tiles in a sample
					clt_parameters.stMeasSel); // bitmask of the selected measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			dbg_hist[1] = scan_prev.showDisparityHistogram();

			dbg_hist[2] = scan_prev.showMaxMinMax();

			int hist_width0 =  scan_prev.showDisparityHistogramWidth();
			int hist_height0 = dbg_hist[0].length/hist_width0;
			if (show_st){
				sdfa_instance.showArrays(dbg_hist, hist_width0, hist_height0, true, "disparity_supertiles_histograms",dbg_st_titles);
			}

			scan_prev.getBgDispStrength( // calculate (check non-null)?
					clt_parameters.stMinBgDisparity, // final double minBgDisparity, 
					clt_parameters.stMinBgFract); // final double minBgFract);


			dbg_orig_disparity = scan_prev.getDisparity().clone();
			
			// combine weak with supertiles
			dbg_with_super_disp = scan_prev.combineSuper(
					true, //boolean updateStrength, // use ST strength if true, keep original (update disparity only) if false
					clt_parameters.stStrengthScale, // Multiply st strength if used instead of regular strength (only if updateStrength)
					clt_parameters.stUseDisp); //.15;  // Use background disparity from supertiles if tile strength is less 
			
			
			if (dbg_with_super_disp != null) dbg_with_super_disp = dbg_with_super_disp.clone(); // else no super disparity available
		}	
		// replace weak outlaye tiles with weighted averages (modifies disparity)
		boolean[] outlayers = scan_prev.replaceWeakOutlayers(
				null, // final boolean [] selection,
				clt_parameters.outlayerStrength , //final double weakStrength,    // strength to be considered weak, subject to this replacement
				clt_parameters.outlayerDiff, // final double maxDiff)
				clt_parameters.outlayerDiffPos, // final double maxDiff)
				clt_parameters.outlayerDiffNeg, // final double maxDiff)
				0.5 * disparity_far,
				2.0 * disparity_near,
				debugLevel);
		dbg_outlayers = new double[outlayers.length];

		for (int i = 0; i < outlayers.length; i++){
			dbg_outlayers[i] = outlayers[i]? 1.0:0.0;
		}			
		
		// set disparity for border pixels (may be overkill)

		DisparityProcessor dp = new DisparityProcessor(this, clt_parameters.transform_size * geometryCorrection.getScaleDzDx());
		boolean [] grown = these_tiles.clone();
		growTiles(
				2,          // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				grown,      // boolean [] tiles,
				null);     // boolean [] prohibit)
		boolean [] border = grown.clone();
		for (int i = 0; i < border.length; i++) border[i] &= !these_tiles[i];
//		double [] dbg_before = scan_prev.getDisparity().clone();
		scan_prev.fixNaNDisparity(
				grown, // border, // boolean [] select,   // which tiles to correct (null - all)
				scan_prev.getDisparity(), // double [] disparity,
				scan_prev.getStrength()); // double [] strength)
//		double [] dbg_after = scan_prev.getDisparity().clone();

		int [] neighbors = dp.getNeighbors( // creates neighbors mask from bitmask
				grown, // these_tiles, // grown, // these_tiles, // boolean [] selected,
				tilesX);
//		double [] dbg_after1 = scan_prev.getDisparity().clone();
		if (clt_parameters.show_neighbors && (debugLevel > 1)) {
			double [] dbg_neib = dp.dbgShowNeighbors(
					grown, // these_tiles, // grown, // these_tiles,
					neighbors, // _orig, // int [] neighbors,
					clt_parameters.transform_size, // int    tile_size,
					-1.0, // double bgnd,
					1.0); // double fgnd)
			sdfa_instance.showArrays(dbg_neib,tilesX*clt_parameters.transform_size, tilesY*clt_parameters.transform_size,"XXneighbors");
		}
//		double [] dbg_after2 = scan_prev.getDisparity().clone();

		
		
		dp.smoothDisparity(
				clt_parameters.tiDispPull,   // final double     dispPull, // clt_parameters.tiDispPull or 0.0
				2, // 2, // 3,                           // final int        mask,     // 1 - work on internal elements, 2 - on border elements, 3 - both (internal first);
				clt_parameters.tiIterations, //  final int        num_passes,
				Math.pow(10.0,  -clt_parameters.tiPrecision), // final double     maxDiff, // maximal change in any of the disparity values
				neighbors,                   // final int     [] neighbors, // +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
				scan_prev.getDisparity(),                // final double  [] disparity,          // current disparity value
				scan_prev.getDisparity().clone(),        // final double  [] measured_disparity, // measured disparity
				scan_prev.getStrength(),               // final double  [] strength,
				null, // this_hor_disparity,          // final double     hor_disparity, // not yet used
				null, // hor_strength_conv,           // final double     hor_strength, // not yet used
				these_tiles,                 // grown, // these_tiles,                 // final boolean [] selected,
				border, // final boolean [] border,
				clt_parameters,
				threadsMax,                  // maximal number of threads to launch                         
				debugLevel); // (clt_3d_passes.size() == 3)? 2: 0); // debugLevel);
//		double [] dbg_after3 = scan_prev.getDisparity().clone();
		
		/*		
		double [] measured_disparity = 	dp.dbgRescaleToPixels(
				this_disparity,
				clt_parameters.transform_size); // int    tile_size)
		 */		
		// prepare new task and run
		double [][] disparityTask = new double [tilesY][tilesX];
		int [][]    tile_op =       new int [tilesY][tilesX];
		boolean [] borderTiles =    new boolean[tilesY*tilesX]; // to zero alpha in the images
		int op = ImageDtt.setImgMask(0, 0xf);
		op =     ImageDtt.setPairMask(op,0xf);
		op =     ImageDtt.setForcedDisparity(op,true);
		double [] prev_disparity = scan_prev.getDisparity();
		for (int ty = 0; ty < tilesY; ty++) for (int tx = 0; tx <tilesX; tx++){
			int indx =  tilesX * ty + tx;
			if ( grown[indx]) {
				borderTiles[indx] =  !these_tiles[indx];
				disparityTask[ty][tx] = prev_disparity[indx];
				tile_op[ty][tx] = op;
			} else {
				disparityTask[ty][tx] = 0.0;
				tile_op[ty][tx] = 0;
				borderTiles[indx] = false;
			}
		}


		double [] masked_filtered = scan_prev.getDisparity().clone();
		for (int i = 0; i < masked_filtered.length; i++){
			if (!grown[i])  masked_filtered[i] = Double.NaN;
		}
		if (show_super){
			String [] dbg_disp_tiltes={"masked", "filtered", "disp_combo", "disparity","st_disparity", "strength",
					"st_strength","outlayers","these","border","border_tiles"}; // ,"before","after","after1","after2","after3","neib"};
			double [][] dbg_disp = new double [dbg_disp_tiltes.length][];
			dbg_disp[0] = masked_filtered;            // +
			dbg_disp[1] = scan_prev.getDisparity();   // +
			dbg_disp[2] = dbg_with_super_disp;        // -
			dbg_disp[3] = dbg_orig_disparity;         // -
			dbg_disp[4] = scan_prev.getBgDisparity(); // -
			dbg_disp[5] = scan_prev.getStrength();    // +
			dbg_disp[6] = scan_prev.getBgStrength();  // -
			dbg_disp[7] = dbg_outlayers;              // +
			dbg_disp[8] = new double [masked_filtered.length];
			dbg_disp[9] = new double [masked_filtered.length];
			dbg_disp[10] = new double [masked_filtered.length];
//			dbg_disp[16] = new double [masked_filtered.length];
			for (int i = 0; i < dbg_disp[8].length; i++){
				dbg_disp[8][i] = these_tiles[i]? 1.0: 0.0;
				dbg_disp[9][i] = border[i]? 1.0: 0.0;
				dbg_disp[10][i] = borderTiles[i]? 1.0: 0.0;
//				dbg_disp[16][i] = neighbors[i];
			}
//			dbg_disp[11] = dbg_before;
//			dbg_disp[12] = dbg_after;
//			dbg_disp[13] = dbg_after1;
//			dbg_disp[14] = dbg_after2;
//			dbg_disp[15] = dbg_after3;

			sdfa_instance.showArrays(dbg_disp, tilesX, tilesY, true, "refine_disparity_supertiles-"+clt_3d_passes.size(),dbg_disp_tiltes);
		}

		CLTPass3d scan_next =new CLTPass3d(this);
		scan_next.disparity =     disparityTask;
		scan_next.tile_op =       tile_op;
		scan_next.border_tiles =  borderTiles;
		scan_next.selected =      grown; // includes border_tiles
		clt_3d_passes.add(scan_next);
		//		}
		return scan_next;
	}


//======================
	public void showPlanes(
			EyesisCorrectionParameters.CLTParameters           clt_parameters,
			GeometryCorrection geometryCorrection,
			final int         threadsMax,  // maximal number of threads to launch                         
			final boolean     updateStatus,
			final int         debugLevel)
	{
		trimCLTPasses(); // make possible to run this method multiple time - remove extra passes added by it last time		
		CLTPass3d scan_prev = clt_3d_passes.get(clt_3d_passes.size() -1); // get last one

		boolean show_st =    clt_parameters.stShow || (debugLevel > 1);
		// recalculate supertiles (may be removed later)
		//		if (use_supertiles || show_st) {	
		String [] dbg_st_titles = {"raw", "sampled", "blurred"+clt_parameters.stSigma,"max-min-max"}; 
		double [][] dbg_hist = new double[dbg_st_titles.length][];
		if (show_st) { // otherwise only blured version is needed 
			scan_prev.setSuperTiles(
					clt_parameters.stStepNear,       // double     step_disparity,
					clt_parameters.stStepFar,        // double     step_near,
					clt_parameters.stStepThreshold,  // double     step_threshold,
					clt_parameters.stMinDisparity,   // double     min_disparity,
					clt_parameters.stMaxDisparity,   // double     max_disparity,
					clt_parameters.stFloor,          // double     strength_floor,
					clt_parameters.stPow,            // double     strength_pow,
					0.0, // NO BLUR double     stBlurSigma)
					false, //clt_parameters.stSmplMode,  // Use sample mode (false - regular tile mode)
					clt_parameters.stSmplSide,  // Sample size (side of a square)
					clt_parameters.stSmplNum,   // Number after removing worst
					clt_parameters.stSmplRms,   // Maximal RMS of the remaining tiles in a sample
					clt_parameters.stMeasSel); // bitmask of the selected measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			dbg_hist[0] = scan_prev.showDisparityHistogram();
			scan_prev.setSuperTiles(
					clt_parameters.stStepNear,       // double     step_disparity,
					clt_parameters.stStepFar,        // double     step_near,
					clt_parameters.stStepThreshold,  // double     step_threshold,
					clt_parameters.stMinDisparity,   // double     min_disparity,
					clt_parameters.stMaxDisparity,   // double     max_disparity,
					clt_parameters.stFloor,          // double     strength_floor,
					clt_parameters.stPow,            // double     strength_pow,
					0.0, // NO BLUR double     stBlurSigma)
					clt_parameters.stSmplMode,  // Use sample mode (false - regular tile mode)
					clt_parameters.stSmplSide,  // Sample size (side of a square)
					clt_parameters.stSmplNum,   // Number after removing worst
					clt_parameters.stSmplRms,   // Maximal RMS of the remaining tiles in a sample
					clt_parameters.stMeasSel); // bitmask of the selected measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			dbg_hist[1] = scan_prev.showDisparityHistogram();
		}

//		SuperTiles st = 
				scan_prev.setSuperTiles(
				clt_parameters.stStepNear,       // double     step_disparity,
				clt_parameters.stStepFar,        // double     step_near,
				clt_parameters.stStepThreshold,  // double     step_threshold,
				clt_parameters.stMinDisparity,   // double     min_disparity,
				clt_parameters.stMaxDisparity,   // double     max_disparity,
				clt_parameters.stFloor,          // double     strength_floor,
				clt_parameters.stPow,            // double     strength_pow,
				clt_parameters.stSigma,          // with blur double     stBlurSigma)
				false, //clt_parameters.stSmplMode,  // Use sample mode (false - regular tile mode)
				clt_parameters.stSmplSide,  // Sample size (side of a square)
				clt_parameters.stSmplNum,   // Number after removing worst
				clt_parameters.stSmplRms,   // Maximal RMS of the remaining tiles in a sample
				clt_parameters.stMeasSel); // bitmask of the selected measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
		if (show_st) { // otherwise only blured version is needed 
			dbg_hist[2] = scan_prev.showDisparityHistogram();
			dbg_hist[3] = scan_prev.showMaxMinMax();
		}


		if (show_st){
			int hist_width0 =  scan_prev.showDisparityHistogramWidth();
			int hist_height0 = dbg_hist[0].length/hist_width0;
			showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
			sdfa_instance.showArrays(dbg_hist, hist_width0, hist_height0, true, "disparity_supertiles_histograms",dbg_st_titles);
		}
		
		if (clt_parameters.stSmplMode){ // just temporary making it last
			scan_prev.setSuperTiles(
					clt_parameters.stStepNear,       // double     step_disparity,
					clt_parameters.stStepFar,        // double     step_near,
					clt_parameters.stStepThreshold,  // double     step_threshold,
					clt_parameters.stMinDisparity,   // double     min_disparity,
					clt_parameters.stMaxDisparity,   // double     max_disparity,
					clt_parameters.stFloor,          // double     strength_floor,
					clt_parameters.stPow,            // double     strength_pow,
					0.0, // NO BLUR double     stBlurSigma)
					clt_parameters.stSmplMode,  // Use sample mode (false - regular tile mode)
					clt_parameters.stSmplSide,  // Sample size (side of a square)
					clt_parameters.stSmplNum,   // Number after removing worst
					clt_parameters.stSmplRms,   // Maximal RMS of the remaining tiles in a sample
					clt_parameters.stMeasSel); // bitmask of the selected measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
		}
		
		
		
		SuperTiles st = scan_prev.getSuperTiles();
		
// moved here
		if (clt_parameters.dbg_migrate) {
		st.processPlanes4(
				null, // final boolean [] selected, // or null
				0.3, // final double     min_disp,
				clt_parameters.stMeasSel, //            =     1   //Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
				clt_parameters.plDispNorm, //           =   2.0;  // Normalize disparities to the average if above
				clt_parameters.plMinPoints, //          =     5;  // Minimal number of points for plane detection
				clt_parameters.plTargetEigen, //        =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
				clt_parameters.plFractOutliers, //      =   0.3;  // Maximal fraction of outliers to remove
				clt_parameters.plMaxOutliers, //        =    20;  // Maximal number of outliers to remove\
				clt_parameters.plPreferDisparity,
				geometryCorrection,
				clt_parameters.correct_distortions,
				clt_parameters.stSmplMode , // final boolean                    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				clt_parameters.stSmplSide , // final int                        smplSide, //        = 2;      // Sample size (side of a square)
				clt_parameters.stSmplNum ,  // final int                        smplNum, //         = 3;      // Number after removing worst
				clt_parameters.stSmplRms ,  // final double                     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
				clt_parameters.stSmallDiff, //       = 0.4;   // Consider merging initial planes if disparity difference below
				clt_parameters.stHighMix, //stHighMix         = 0.4;   // Consider merging initial planes if jumps between ratio above
				clt_parameters.vertical_xyz,
				0, // -1, // debugLevel,                  // final int        debugLevel)
				clt_parameters.tileX,
				clt_parameters.tileY);
		} else {
			st.processPlanes3(
					null, // final boolean [] selected, // or null
					0.3, // final double     min_disp,
					clt_parameters.stMeasSel, //            =     1   //Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
					clt_parameters.plDispNorm, //           =   2.0;  // Normalize disparities to the average if above
					clt_parameters.plMinPoints, //          =     5;  // Minimal number of points for plane detection
					clt_parameters.plTargetEigen, //        =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
					clt_parameters.plFractOutliers, //      =   0.3;  // Maximal fraction of outliers to remove
					clt_parameters.plMaxOutliers, //        =    20;  // Maximal number of outliers to remove\
					clt_parameters.plPreferDisparity,
					geometryCorrection,
					clt_parameters.correct_distortions,
//					false, // clt_parameters.dbg_migrate && clt_parameters.stSmplMode , // final boolean                    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
					clt_parameters.stSmplMode , // final boolean                    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
					clt_parameters.stSmplSide , // final int                        smplSide, //        = 2;      // Sample size (side of a square)
					clt_parameters.stSmplNum , // final int                        smplNum, //         = 3;      // Number after removing worst
					clt_parameters.stSmplRms , // final double                     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
					
					0, // -1, // debugLevel,                  // final int        debugLevel)
					clt_parameters.tileX,
					clt_parameters.tileY);
			/*
			st.processPlanes2(
					null, // final boolean [] selected, // or null
					0.3, // final double     min_disp,
					false, // final boolean    invert_disp, // use 1/disparity
					clt_parameters.plDispNorm, //            =   2.0;  // Normalize disparities to the average if above
					clt_parameters.plMinPoints, //           =     5;  // Minimal number of points for plane detection
					clt_parameters.plTargetEigen, //         =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
					clt_parameters.plFractOutliers, //      =   0.3;  // Maximal fraction of outliers to remove
					clt_parameters.plMaxOutliers, //        =    20;  // Maximal number of outliers to remove\
					clt_parameters.plPreferDisparity,
					geometryCorrection,
					clt_parameters.correct_distortions,
					0, // -1, // debugLevel,                  // final int        debugLevel)
					clt_parameters.tileX,
					clt_parameters.tileY);
					*/
		}
		
		showDoubleFloatArrays sdfa_instance = null;
		if (debugLevel > -1) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?

		st.matchPlanes(
				clt_parameters.plPreferDisparity,
				clt_parameters.tileX,
				clt_parameters.tileY); 

		st.selectNeighborPlanesMutual(
				clt_parameters.plWorstWorsening, // final double worst_worsening,
				clt_parameters.plWeakWorsening, // final double worst_worsening,
				clt_parameters.plDispNorm,
				clt_parameters.plMaxEigen,
				clt_parameters.plMinStrength,
				0, // final int debugLevel)
				clt_parameters.tileX,
				clt_parameters.tileY);
		
		if (clt_parameters.plSplitApply) {
			while (true) {
				int num_added = 0;
				num_added += st.fillSquares();
				if (debugLevel > -1) {
					System.out.println("after fillSquares() added "+num_added);
				}
				num_added += st.cutCorners();
				if (debugLevel > -1) {
					System.out.println("after plCutCorners() added (cumulative) "+num_added);
				}
				if (num_added == 0) break;
			}
		}		
		
		TilePlanes.PlaneData[][][]       split_planes =   // use original (measured planes. See if smoothed are needed here)
				st.breakPlanesToPairs(
				st.getPlanes(), // Mod(),             // final TilePlanes.PlaneData[][] center_planes, // measured_planes,
				st.getPlanes(), // Mod(),             // final TilePlanes.PlaneData[][] neib_planes,   //mod_planes,
				clt_parameters.plSplitPull ,          // final double                   center_pull,
				clt_parameters.plSplitMinNeib ,       // min_neibs, // 2
				clt_parameters.plSplitMinWeight,      // final double splitMinWeight,  //     =  2.0;  // Minimal weight of split plains to show
				clt_parameters.plSplitMinQuality,     // final double splitMinQuality, //    =  1.1;  // Minimal split quality to show

				clt_parameters.plPreferDisparity,
				1,                                   // final int debugLevel)
				clt_parameters.tileX,
				clt_parameters.tileY); 
		
		if (clt_parameters.plSplitApply) {
			int numSplitPlanes = st.replaceBrokenPlanes(
					st.getPlanes(),                 // final TilePlanes.PlaneData[][]   planes,
					split_planes,                    // final TilePlanes.PlaneData[][][] brokenPd,
					clt_parameters.plMaxDiff,        //  final double                     max_diff, // maximal disparity difference (0 - any), will be normalized by dispNorm
					clt_parameters.plOtherDiff,      //  final double                     other_diff, // other_diff maximal difference of the added tile ratio to the average disparity difference
					clt_parameters.plNonExclusive,   // final boolean                    non_exclusive,
					clt_parameters.plUseOtherPlanes, // final boolean                    use_other_planes, // TODO:
					clt_parameters.stMeasSel,        // final int                        measSel, // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
					clt_parameters.plAllowParallel,  // final boolean                    allow_parallel,
					clt_parameters.plSplitXY,        // final boolean                    splitXY,
					clt_parameters.plSplitXYTolerance, // final double                   splitXYTolerance,
					// parameters to generate ellipsoids			
					0.0, // 3,                             // final double                     disp_far, // minimal disparity to select (or NaN)
					Double.NaN,                      // final double                     disp_near, // maximal disparity to select (or NaN)
					clt_parameters.plDispNorm,       // final double                     dispNorm,   //  Normalize disparities to the average if above
					0.0,                             // final double                     min_weight,
					clt_parameters.plMinPoints,      // final int                        min_tiles,
					// parameters to reduce outliers			
					clt_parameters.plTargetEigen,    // final double                     targetEigen,   //     =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
					clt_parameters.plFractOutliers,  // final double                     fractOutliers, //     =   0.3;  // Maximal fraction of outliers to remove
					clt_parameters.plMaxOutliers,    // final int                        maxOutliers,   //     =   20;  // Maximal number of outliers to remove
					clt_parameters.stFloor,          // final double                     strength_floor,
					clt_parameters.stPow,            // final double                     strength_pow,
					clt_parameters.dbg_migrate && clt_parameters.stSmplMode , // final boolean                    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
					clt_parameters.stSmplSide , // final int                        smplSide, //        = 2;      // Sample size (side of a square)
					clt_parameters.stSmplNum , // final int                        smplNum, //         = 3;      // Number after removing worst
					clt_parameters.stSmplRms , // final double                     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
					1,                               // final int debugLevel)
					clt_parameters.tileX,
					clt_parameters.tileY);
			if (debugLevel > -1){
				System.out.println("replaceBrokenPlanes(): Replaced " + numSplitPlanes + " planes by pairs.");
			}
			// now re-create connections and best neighbors
			st.matchPlanes(
					clt_parameters.plPreferDisparity,
					clt_parameters.tileX,
					clt_parameters.tileY); 

			st.selectNeighborPlanesMutual(
					clt_parameters.plWorstWorsening, // final double worst_worsening,
					clt_parameters.plWeakWorsening, // final double worst_worsening,
					clt_parameters.plDispNorm,
					clt_parameters.plMaxEigen,
					clt_parameters.plMinStrength,
					0, // final int debugLevel)
					clt_parameters.tileX,
					clt_parameters.tileY);
		}

		while (true) {
			int num_added = 0;
			if (clt_parameters.plFillSquares){
				num_added += st.fillSquares();
			}
			if (debugLevel > -1) {
				System.out.println("after fillSquares() added "+num_added);
			}
			if (clt_parameters.plCutCorners){
				num_added += st.cutCorners();
			}
			if (debugLevel > -1) {
				System.out.println("after plCutCorners() added (cumulative) "+num_added);
			}
			if (num_added == 0) break;
		}
		
		
		TilePlanes.PlaneData [][] planes_mod = null;

		// smooth planes (by averaging with neighbors and the "measured" one with variable "pull")  
		if (clt_parameters.plIterations > 0) {
			st.resetPlanesMod(); // clean start
			planes_mod = st.planesSmooth(
					clt_parameters.plPull,                        // final double      meas_pull,//  relative pull of the original (measured) plane with respect to the average of the neighbors
					clt_parameters.plMaxEigen,                    // final double      maxValue, // do not combine with too bad planes
					clt_parameters.plIterations,                  // final int         num_passes,
					clt_parameters.plStopBad,                     // Do not update supertile if any of connected neighbors is not good (false: just skip that neighbor)
					clt_parameters.plNormPow,                     // 0.0: 8 neighbors pull 8 times as 1, 1.0 - same as 1
					Math.pow(10.0,  -clt_parameters.plPrecision), // final double      maxDiff, // maximal change in any of the disparity values
					clt_parameters.plPreferDisparity,
					0, // final int debugLevel)
					clt_parameters.tileX,
					clt_parameters.tileY); 
		}
		// filter out weak planes, create boolean array [per-supertile][per disparity plane]
		boolean [][] selected_planes =	st.selectPlanes(
				clt_parameters.plDispNorm,
				clt_parameters.plMaxEigen, //  maxEigen,
				clt_parameters.plMinStrength, //  minWeight,
				st.getPlanesMod());
        // detect connected "shells" by running wave algorithm over multi-plane supertiles, create integer array (same structure as planes and selected_planes
		// Each element [per supertile][per disparity plane] is either 0 (not used) or unique shell index (started from 1)
		int [][] shells = st.createSupertileShells(
				selected_planes, // boolean[][]               selection, // may be null
				false,             // boolean                   use_all, // use plane 0 even if there are more than 1
				clt_parameters.plKeepOrphans,             // true, // boolean                   keep_orphans, // single-cell shells
				clt_parameters.plMinOrphan, // orphan_strength // add separate parameter
				st.getPlanesMod(), // TilePlanes.PlaneData [][] planes,
				1); // int                       debugLevel)
		// save shell indices with SuperTiles instance
		st.setShellMap(shells); // persistent
		// Create array [per shell][per tile] of shell disparities. tiles of unused supertiles are Double.NaN
		double [][] surfaces = st.getShowShells(
				st.getPlanesMod(),   // TilePlanes.PlaneData [][] planes,
				st.getShellMap(), // shells,              // int [][] shells,
				1000,                 // int max_shells,
				clt_parameters.plFuse,// boolean fuse,
				false,               // boolean show_connections,
				false,               // boolean use_NaN,
				10.0,                 // double arrow_dark,
				10.0);               // double arrow_white)
		// save surfaces with SuperTiles instance. They can be used to snap to for the per-tile disparity maps.
		st.setSurfaces(surfaces);

/*
  		TilePlanes.PlaneData[][][]       split_planes =  
				st.breakPlanesToPairs(
				st.getPlanes(), // Mod(),             // final TilePlanes.PlaneData[][] center_planes, // measured_planes,
				st.getPlanes(), // Mod(),             // final TilePlanes.PlaneData[][] neib_planes,   //mod_planes,
				clt_parameters.plSplitPull ,          // final double                   center_pull,
				clt_parameters.plSplitMinNeib ,       // min_neibs, // 2
				clt_parameters.plSplitMinWeight,      // final double splitMinWeight,  //     =  2.0;  // Minimal weight of split plains to show
				clt_parameters.plSplitMinQuality,     // final double splitMinQuality, //    =  1.1;  // Minimal split quality to show

				clt_parameters.plPreferDisparity,
				1,                                   // final int debugLevel)
				clt_parameters.tileX,
				clt_parameters.tileY); 
*/		
		
		if (clt_parameters.show_planes){
			double [] split_lines = st.showSplitLines(
					split_planes,
					1,                              // final int debugLevel)
					clt_parameters.tileX,
					clt_parameters.tileY); 

			int [] wh = st.getShowPlanesWidthHeight();
			double [][] plane_data_nonan = st.getShowPlanes(
					(planes_mod != null) ? st.getPlanesMod():st.getPlanes(),
					clt_parameters.plMinStrength, //  minWeight,
					clt_parameters.plMaxEigen, //  maxEigen,
					clt_parameters.plDispNorm,
					false, //boolean use_NaN)
					0.0,
					10.0);
			double [][] plane_data_nan = st.getShowPlanes(
					(planes_mod != null) ? st.getPlanesMod():st.getPlanes(),
					clt_parameters.plMinStrength, //  minWeight,
					clt_parameters.plMaxEigen, //  maxEigen,
					clt_parameters.plDispNorm,
					true, //boolean use_NaN)
					0.0,
					10.0);
			double [][] plane_data = new double [plane_data_nonan.length + plane_data_nan.length + 2][];
			int indx = 0;
			for (int i = 0; i < plane_data_nonan.length; i++){
				plane_data[indx++] = plane_data_nonan[i];
			}
			for (int i = 0; i < plane_data_nan.length; i++){
				plane_data[indx++] = plane_data_nan[i];
			}
			plane_data[indx++] = split_lines;
			plane_data[indx] = plane_data[indx-2].clone(); // java.lang.ArrayIndexOutOfBoundsException: -1
			for (int i = 0; i < plane_data[indx].length;i++){
				if (Double.isNaN(plane_data[indx][i])) plane_data[indx][i] = 0.0; 
				if (plane_data[indx-1][i] > 0) plane_data[indx][i] = Double.NaN; 
			}
			
			sdfa_instance.showArrays(plane_data, wh[0], wh[1], true, "plane_data");
			
			plane_data = st.getShowShells(
					st.getPlanesMod(),   // TilePlanes.PlaneData [][] planes,
					st.getShellMap(), // shells,              // int [][] shells,
					100,                 // int max_shells,
					clt_parameters.plFuse,// boolean fuse,
					false,               // boolean show_connections,
					false,               // boolean use_NaN,
					10.0,                 // double arrow_dark,
					10.0);               // double arrow_white)
			
			sdfa_instance.showArrays(plane_data, wh[0], wh[1], true, "shells_all");
			
			plane_data = st.getShowShells(
					st.getPlanesMod(),   // TilePlanes.PlaneData [][] planes,
					st.getShellMap(), // shells,              // int [][] shells,
					100,                 // int max_shells,
					clt_parameters.plFuse,// boolean fuse,
					true,                // boolean show_connections,
					false,               // boolean use_NaN,
					10.0,                 // double arrow_dark,
					10.0);               // double arrow_white)
			sdfa_instance.showArrays(plane_data, wh[0], wh[1], true, "shells_all_arrows");

			
			plane_data = st.getShowShells(
					st.getPlanesMod(),   // TilePlanes.PlaneData [][] planes,
					st.getShellMap(), // shells,              // int [][] shells,
					100,                 // int max_shells,
					clt_parameters.plFuse,// boolean fuse,
					false,               // boolean show_connections,
					true,               // boolean use_NaN,
					10.0,                 // double arrow_dark,
					10.0);               // double arrow_white)
			sdfa_instance.showArrays(plane_data, wh[0], wh[1], true, "shells");
// also test snap here			
			int [] snap_surf =  st.snapSort(
					st.cltPass3d.getDisparity(),                  // final double []  disparity,
					st.cltPass3d.getStrength(),                   // final double []  strength,
					st.cltPass3d.getSelected(),      // null, // this_selection,         // final boolean [] selection,       // can be null
					clt_parameters.plDispNorm,       //  final double     dispNorm, // plDispNorm
					clt_parameters.plSnapDispAny,    //final double     snapDispAny,    //    =  .2,   // Maximal (scaled by plDispNorm) disparity difference to snap to plane at any strength
					clt_parameters.plSnapStrengthAny,//final double     snapStrengthAny,    //    =  .2,   //  Maximal strength to fit any distance (if does not fit otherwise - treat as zero strength
					clt_parameters.plSnapNegAny,     //final double     snapNegAny,     //    =  .2,   // Maximal negative disparity difference from the best match
					clt_parameters.plSnapDispMax,    //  final double     snapDispMax,    //    =  .5,   // Maximal (scaled by plDispNorm) disparity difference to snap to plane at low strength
					clt_parameters.plSnapDispWeight, //  final double     snapDispWeight, //    =  .5,   // Maximal disparity diff. by weight product to snap to plane
					clt_parameters.plSnapZeroMode,   // final int        snapZeroMode,   // Zero strength snap mode: 0: no special treatment, 1 - strongest, 2 - farthest
					1); // final int        debugLevel)
			
			double [] snap_disparity =  st.snapDisparity(
					st.cltPass3d.getDisparity(),                  // final double []  disparity,
					st.cltPass3d.getStrength(),                   // final double []  strength,
					null, // this_selection,         // final boolean [] selection,       // can be null
					clt_parameters.plDispNorm,       //  final double     dispNorm, // plDispNorm
					clt_parameters.plSnapDispAny,    //final double     snapDispAny,    //    =  .2,   // Maximal (scaled by plDispNorm) disparity difference to snap to plane at any strength
					clt_parameters.plSnapStrengthAny,//final double     snapStrengthAny,    //    =  .2,   //  Maximal strength to fit any distance (if does not fit otherwise - treat as zero strength
					clt_parameters.plSnapNegAny,     //final double     snapNegAny,     //    =  .2,   // Maximal negative disparity difference from the best match
					clt_parameters.plSnapDispMax,    //  final double     snapDispMax,    //    =  .5,   // Maximal (scaled by plDispNorm) disparity difference to snap to plane at low strength
					clt_parameters.plSnapDispWeight, //  final double     snapDispWeight, //    =  .5,   // Maximal disparity diff. by weight product to snap to plane
					clt_parameters.plSnapZeroMode,   // final int        snapZeroMode,   // Zero strength snap mode: 0: no special treatment, 1 - strongest, 2 - farthest
					1); // final int        debugLevel)
			double [] snap_disparity_masked = snap_disparity.clone();
			double [] disparity_masked =    st.cltPass3d.getDisparity().clone();
			double [] strength_masked =     st.cltPass3d.getStrength().clone();
			boolean [] selected =           st.cltPass3d.getSelected();
			for (int i = 0; i < snap_disparity_masked.length; i++){
				if ( !selected[i]){
					disparity_masked[i] = Double.NaN;
					snap_disparity_masked[i] = Double.NaN;
					strength_masked[i] = Double.NaN;
				}
			}
			String [] snap_titles = {"disp","snap_disp","strength","disp_masked","snap_disp_masked","unmatched","strength_masked"};
			double [][] snap_img = new double [snap_titles.length][];
			snap_img[5] = new double [snap_surf.length];
			for (int i = 0; i < snap_surf.length; i++){
				snap_img[5][i] = (selected[i] && (snap_surf[i] == 0)) ? 10.0: 0.0;
			}
			snap_img[0] = st.cltPass3d.getDisparity();
			snap_img[1] = snap_disparity;
			snap_img[2] = st.cltPass3d.getStrength();
			snap_img[3] = disparity_masked;
			snap_img[4] = snap_disparity_masked;
			snap_img[6] = strength_masked;
			sdfa_instance.showArrays(snap_img, tilesX, tilesY, true, "snap",snap_titles);
			

			boolean [] these_tiles = scan_prev.getSelected();



			DisparityProcessor dp = new DisparityProcessor(this, clt_parameters.transform_size * geometryCorrection.getScaleDzDx());
			boolean [] grown = these_tiles.clone();
			growTiles(
					2,          // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
					grown,      // boolean [] tiles,
					null);     // boolean [] prohibit)
			boolean [] border = grown.clone();
			for (int i = 0; i < border.length; i++) border[i] &= !these_tiles[i];


			int num_unmatched = st.getNumNotMatched(
					snap_surf,
					these_tiles);
			int num_surf = st.getNumSurf(
					snap_surf,
					null); // these_tiles);
			if (debugLevel > -1){
				System.out.println("thirdPassSetup(): surfaces = "+num_surf+" left unmatched="+num_unmatched);
			}
			int [][] neib_surf =   new int [num_surf][];
			double [][] dbg_neibs = new double [num_surf][];
			for (int ns = 0 ; ns < num_surf; ns ++){

				neib_surf[ns] = st.getNeighbors( // creates neighbors mask from bitmask
						these_tiles, // final boolean [] selected, // or null
						snap_surf,   // final int []     snap_surf, // use this size - it matches image, not supertiles
						ns);        // final int        surf_index)

				dbg_neibs[ns]= dp.dbgShowNeighbors(
						these_tiles, // grown, // these_tiles, 
						neib_surf[ns], // _orig, // int [] neighbors,
						clt_parameters.transform_size, // int    tile_size,
						-1.0, // Double.NaN, // double bgnd,
						10.0); // double fgnd)
			}

			sdfa_instance.showArrays(dbg_neibs,
					tilesX*clt_parameters.transform_size,
					tilesY*clt_parameters.transform_size,
					true,
					"surf-neighbors");		
			
		}
	}
	
	public void secondPassSetup( // prepare tile tasks for the second pass based on the previous one(s)
			//			  final double [][][]       image_data, // first index - number of image in a quad
			EyesisCorrectionParameters.CLTParameters           clt_parameters,
			boolean           use_supertiles,
			int               bg_scan_index,
			// disparity range - differences from 
			double            disparity_far,    //
			double            disparity_near,   // 
			double            this_sure,        // minimal strength to be considered definitely background
			double            this_maybe,       // maximal strength to ignore as non-background
			double            sure_smth,        // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
			int               disparity_index,  // index of disparity value in disparity_map == 2 (0,2 or 4)
			GeometryCorrection geometryCorrection,
			final int         threadsMax,  // maximal number of threads to launch                         
			final boolean     updateStatus,
			final int         debugLevel)
	{
		CLTPass3d scan_bg =   clt_3d_passes.get(bg_scan_index); // 
		CLTPass3d scan_prev = clt_3d_passes.get(clt_3d_passes.size() -1); // get last one
//		CLTPass3d scan_next =new CLTPass3d(this);
		showDoubleFloatArrays sdfa_instance = null;
		if (debugLevel > -1) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		//TODO: for next passes - combine all selected for previous passes (all passes with smaller disparity)
		boolean show_st =    clt_parameters.stShow || (debugLevel > 1);
		
		boolean [] these_tiles;
		int [] replaced =  null; // +1 - hor, +2 - vert
		int [] replaced0 = null; // +1 - hor, +2 - vert
		
		if (!clt_parameters.ortho_old) {
			replaced = combineOrthoDisparity(
					scan_prev, // final CLTPass3d   scan,        // scan data
					clt_parameters.or_hor,       // true;  // Apply ortho correction to horizontal correlation (vertical features)
					clt_parameters.or_vert,      // true;  // Apply ortho correction to vertical correlation (horizontal features)
					clt_parameters.or_sigma,     // 2.0;   // Blur sigma: verically for horizontal correlation, horizontally - for vertically
					clt_parameters.or_sharp,     // 0.5;   // 3-point sharpening (-k, +2k+1, -k)
					clt_parameters.or_scale,     // 2.0;   // Scale ortho correletion strength relative to 4-directional one
					clt_parameters.or_offset,    // 0.1;   // Subtract from scaled correlation strength, limit by 0
					clt_parameters.or_asym ,     // 1.5;   // Minimal ratio of orthogonal strengths required for dis[parity replacement
					clt_parameters.or_threshold, // 1.5;   // Minimal scaled offsetg ortho strength to normal strength needed for replacement
					clt_parameters.or_absHor,    // 0.15;  // Minimal horizontal absolute scaled offset ortho strength needed for replacement
					clt_parameters.or_absVert,   // 0.19;  // Minimal vertical absolute scaled offset ortho strength needed for replacement
					debugLevel);

			if (clt_parameters.poles_fix) {
				boolean [] selection = new boolean [replaced.length];
				boolean [] tilesHor =   new boolean [replaced.length];
				double  [] disparity = scan_prev.getDisparity();
				for (int i = 0; i < tilesHor.length; i++){
					tilesHor[i] = (replaced[i] & 1) != 0;
					selection[i] = !Double.isNaN(disparity[i]) && (disparity[i] >= disparity_far) && (disparity[i] <= disparity_near);
				}
				int numFixed = fixVerticalPoles( // return number of replaced cells
						scan_prev, // CLTPass3d   scan,             // scan data to use
						selection,                         // start with only from selections (if not null, continue regardless)
						tilesHor,                          // horizontal correlation tiles used for composite disparity/strength;
						clt_parameters.poles_len ,         // int         max_len,          // maximal length to cover
						clt_parameters.poles_ratio,        // Maximal ratio of invisible to visible pole length
						clt_parameters.poles_min_strength, // double      min_new_strength, // set strength to hor_strength, but not less than this
						clt_parameters.poles_force_disp,   // boolean     force_disparity   // copy disparity down (false - use horDisparity
						true); 
				if (debugLevel > -1){
					System.out.println("fixVerticalPoles() replaced "+ numFixed+ " tiles.");
				}
				replaced0 = replaced.clone();
				for (int i = 0; i < replaced.length; i++){
					if (tilesHor[i]) replaced[i] |= 1;
				}
			}
			these_tiles = FilterScan(
					scan_prev,        // final CLTPass3d   scan,
					scan_bg.selected, // get from selected in clt_3d_passes.get(0);
					replaced,         // final int      [] horVertMod, // +1 - modified by hor correlation, +2 - modified by vert correlation (or null)
					disparity_far,    // final double      disparity_far,    //
					disparity_near,   // final double      disparity_near,   // 
					this_sure,        // final double      this_sure,        // minimal strength to be considered definitely background
					this_maybe,       // final double      this_maybe,       // maximal strength to ignore as non-background
					sure_smth,        // final double      sure_smth,        // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
					true, // clt_parameters.show_filter_scan,
					clt_parameters.min_clstr_seed, // clt_parameters.min_clstr_seed
					clt_parameters.min_clstr_weight,  // double     min_weight // minimal total weight of the cluster
					clt_parameters.min_clstr_max,    // double     min_max_weight // minimal value of the maximal strengh in the cluster
					clt_parameters.min_clstr_lone,  // int        min_area,  // minimal number of pixels
					clt_parameters.fill_gaps, // fill_gaps, // int depth, // same as grow - odd - 4 directions, even - 8
					clt_parameters.fill_final,             // final boolean     poison_gaps,      // Do not fill gaps that have even single "poisoned" tile
					true,             // final boolean     zero_gap_strength, // set strength to zero when covering gaps
					//					final int         threadsMax,  // maximal number of threads to launch                         
					//					final boolean     updateStatus,
					2); //debugLevel);
		} else {
			these_tiles= combineHorVertDisparity(
				scan_prev,                     // final CLTPass3d   scan,
				scan_bg.selected, // clt_3d_passes.get(0).selected, // final boolean [] bg_tiles,          // get from selected in clt_3d_passes.get(0); 
				disparity_far,    //
				disparity_near,   // 
				this_sure,        // minimal strength to be considered definitely background
				this_maybe,       // maximal strength to ignore as non-background
				sure_smth,        // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
				clt_parameters,
				debugLevel);
		
			scan_prev.combineHorVertStrength(true, false); // strength now max of original and horizontal. Use scale  instead of boolean?
		}
		
		
		
		double [] this_disparity     = scan_prev.getDisparity(); // returns a copy of the FPGA-generated disparity combined with the target one
		double [] this_strength =      scan_prev.getStrength(); // cloned, can be modified/ read back
		
		String []  dbg_test_titles = {"was_selected","these_tiles", "disparity","strength"};
		double [][] dbg_test_img = new double [dbg_test_titles.length][];
		dbg_test_img[0] = new double [tilesX * tilesY];
		dbg_test_img[1] = new double [tilesX * tilesY];
		boolean [] dbg_test_selected = scan_prev.getSelected();
		for (int i = 0; i < dbg_test_img[0].length; i++){
			if (dbg_test_selected != null) dbg_test_img[0][i] =  dbg_test_selected[i]?1.0:0.0; 
			dbg_test_img[1][i] =  these_tiles[i]?1.0:0.0; 
		}
		dbg_test_img[2] = this_disparity;
		dbg_test_img[3] = this_strength;
		if (sdfa_instance != null) sdfa_instance.showArrays(dbg_test_img, tilesX, tilesY, true, "dbg_test_img",dbg_test_titles);
		
		scan_prev.setSelected(these_tiles); // New
		
		
		//************************************************
		// Show supertiles histograms

		//		if (clt_parameters.stShow){

		// try renovated supertiles. Do twice to show both original and blured histograms
		double [] dbg_orig_disparity =  null;
		double [] dbg_with_super_disp = null;
		double [] dbg_outlayers =       null;
//		boolean [] grown = these_tiles.clone();

		if (use_supertiles || show_st) {	
			String [] dbg_st_titles = {"raw", "blurred"+clt_parameters.stSigma,"max-min-max"}; 
			double [][] dbg_hist = new double[dbg_st_titles.length][];

			scan_prev.setSuperTiles(
					clt_parameters.stStepNear,       // double     step_disparity,
					clt_parameters.stStepFar,        // double     step_near,
					clt_parameters.stStepThreshold,  // double     step_threshold,
					clt_parameters.stMinDisparity,   // double     min_disparity,
					clt_parameters.stMaxDisparity,   // double     max_disparity,
					clt_parameters.stFloor,          // double     strength_floor,
					clt_parameters.stPow,            // double     strength_pow,
					0.0, // NO BLUR double     stBlurSigma)
					false, //clt_parameters.stSmplMode,  // Use sample mode (false - regular tile mode)
					clt_parameters.stSmplSide,  // Sample size (side of a square)
					clt_parameters.stSmplNum,   // Number after removing worst
					clt_parameters.stSmplRms,   // Maximal RMS of the remaining tiles in a sample
					clt_parameters.stMeasSel); // bitmask of the selected measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			dbg_hist[0] = scan_prev.showDisparityHistogram();

			SuperTiles st = scan_prev.setSuperTiles(
					clt_parameters.stStepNear,       // double     step_disparity,
					clt_parameters.stStepFar,        // double     step_near,
					clt_parameters.stStepThreshold,  // double     step_threshold,
					clt_parameters.stMinDisparity,   // double     min_disparity,
					clt_parameters.stMaxDisparity,   // double     max_disparity,
					clt_parameters.stFloor,          // double     strength_floor,
					clt_parameters.stPow,            // double     strength_pow,
					clt_parameters.stSigma, // with blur double     stBlurSigma)
					false, //clt_parameters.stSmplMode,  // Use sample mode (false - regular tile mode)
					clt_parameters.stSmplSide,  // Sample size (side of a square)
					clt_parameters.stSmplNum,   // Number after removing worst
					clt_parameters.stSmplRms,   // Maximal RMS of the remaining tiles in a sample
					clt_parameters.stMeasSel); // bitmask of the selected measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			dbg_hist[1] = scan_prev.showDisparityHistogram();

			dbg_hist[2] = scan_prev.showMaxMinMax();

			int hist_width0 =  scan_prev.showDisparityHistogramWidth();
			int hist_height0 = dbg_hist[0].length/hist_width0;

			if (show_st){
				sdfa_instance.showArrays(dbg_hist, hist_width0, hist_height0, true, "disparity_supertiles_histograms",dbg_st_titles);
			}
/**			
			st.processPlanes2(
					null, // final boolean [] selected, // or null
					0.3, // final double     min_disp,
					false, // final boolean    invert_disp, // use 1/disparity
					clt_parameters.plDispNorm, //            =   2.0;  // Normalize disparities to the average if above
					clt_parameters.plMinPoints, //           =     5;  // Minimal number of points for plane detection
					clt_parameters.plTargetEigen, //         =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
					clt_parameters.plFractOutliers, //      =   0.3;  // Maximal fraction of outliers to remove
					clt_parameters.plMaxOutliers, //        =    20;  // Maximal number of outliers to remove
					geometryCorrection,
					clt_parameters.correct_distortions,
					0, // -1, // debugLevel,                  // final int        debugLevel)
					clt_parameters.tileX,
					clt_parameters.tileY);
			
			
			if (use_supertiles) {
				scan_prev.getBgDispStrength( // calculate (check non-null)?
						clt_parameters.stMinBgDisparity, // final double minBgDisparity, 
						clt_parameters.stMinBgFract); // final double minBgFract);


				dbg_orig_disparity = scan_prev.getDisparity().clone();

				// combine weak with supertiles
				dbg_with_super_disp = scan_prev.combineSuper(
						true, //boolean updateStrength, // use ST strength if true, keep original (update disparity only) if false
						clt_parameters.stStrengthScale, // Multiply st strength if used instead of regular strength (only if updateStrength)
						clt_parameters.stUseDisp); //.15;  // Use background disparity from supertiles if tile strength is less


				if (dbg_with_super_disp != null) dbg_with_super_disp = dbg_with_super_disp.clone(); // else no super disparity available
			}
*/			
		}
		
//		return scan_next;
	}
	
//==================
	public void thirdPassSetup( // prepare tile tasks for the second pass based on the previous one(s)
			EyesisCorrectionParameters.CLTParameters           clt_parameters,
			double            disparity_far,    //
			double            disparity_near,   // 
			GeometryCorrection geometryCorrection,
			final int         threadsMax,  // maximal number of threads to launch                         
			final boolean     updateStatus,
			final int         debugLevel)
	{
		CLTPass3d scan_prev = clt_3d_passes.get(clt_3d_passes.size() -1); // get last one
///		double [] dbg_orig_disparity =  null;
///		double [] dbg_with_super_disp = null;
///		double [] dbg_outlayers =       null;
		showDoubleFloatArrays sdfa_instance = null;
		if (debugLevel > -1) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		boolean [] these_tiles = scan_prev.getSelected();
//		double [] this_disparity     = scan_prev.getDisparity(); // returns a copy of the FPGA-generated disparity combined with the target one
//		double [] this_strength =      scan_prev.getStrength(); // cloned, can be modified/ read back
		double [] this_disparity     = scan_prev.getDisparity().clone(); // returns a copy of the FPGA-generated disparity combined with the target one
		double [] this_strength =      scan_prev.getStrength().clone(); // cloned, can be modified/ read back
		
		scan_prev.fixNaNDisparity(
				null, // border, // boolean [] select,   // which tiles to correct (null - all)
				this_disparity, // scan_prev.getDisparity(), // double [] disparity,
				this_strength); // scan_prev.getStrength()); // double [] strength)
		
///		dbg_orig_disparity = this_disparity.clone();
		
		SuperTiles st = scan_prev.getSuperTiles();

//		if (clt_parameters.plSnapDispAny >= 0.0) {
			double [] snap_disparity =  st.snapDisparity(
					this_disparity,                  // final double []  disparity,
					this_strength,                   // final double []  strength,
					null, // this_selection,         // final boolean [] selection,       // can be null
					clt_parameters.plDispNorm,       //  final double     dispNorm, // plDispNorm
					clt_parameters.plSnapDispAny,    //final double     snapDispAny,    //    =  .2,   // Maximal (scaled by plDispNorm) disparity difference to snap to plane at any strength
					clt_parameters.plSnapStrengthAny,//final double     snapStrengthAny,    //    =  .2,   //  Maximal strength to fit any distance (if does not fit otherwise - treat as zero strength
					clt_parameters.plSnapNegAny,     //final double     snapNegAny,     //    =  .2,   // Maximal negative disparity difference from the best match
					clt_parameters.plSnapDispMax,    //  final double     snapDispMax,    //    =  .5,   // Maximal (scaled by plDispNorm) disparity difference to snap to plane at low strength
					clt_parameters.plSnapDispWeight, //  final double     snapDispWeight, //    =  .5,   // Maximal disparity diff. by weight product to snap to plane
					clt_parameters.plSnapZeroMode,   // final int        snapZeroMode,   // Zero strength snap mode: 0: no special treatment, 1 - strongest, 2 - farthest
					1); // final int        debugLevel)
			System.arraycopy(snap_disparity,  0, this_disparity, 0, snap_disparity.length);
//		}
/**		
// Replace weak outlayers
		if (clt_parameters.replaceWeakOutlayers) {
//			boolean[] outlayers = 
					scan_prev.replaceWeakOutlayers(
					null, // final boolean [] selection,
					clt_parameters.outlayerStrength , //final double weakStrength,    // strength to be considered weak, subject to this replacement
					clt_parameters.outlayerDiff, // final double maxDiff)
					clt_parameters.outlayerDiffPos, // final double maxDiff)
					clt_parameters.outlayerDiffNeg, // final double maxDiff)
					0.5 * disparity_far,
					2.0 * disparity_near,
					debugLevel);

		}
*/		
		int [] snap_surf =  st.snapSort(
				this_disparity,                  // final double []  disparity,
				this_strength,                   // final double []  strength,
				null, // this_selection,         // final boolean [] selection,       // can be null
				clt_parameters.plDispNorm,       //  final double     dispNorm, // plDispNorm
				clt_parameters.plSnapDispAny,    //final double     snapDispAny,    //    =  .2,   // Maximal (scaled by plDispNorm) disparity difference to snap to plane at any strength
				clt_parameters.plSnapStrengthAny,//final double     snapStrengthAny,    //    =  .2,   //  Maximal strength to fit any distance (if does not fit otherwise - treat as zero strength
				clt_parameters.plSnapNegAny,     //final double     snapNegAny,     //    =  .2,   // Maximal negative disparity difference from the best match
				clt_parameters.plSnapDispMax,    //  final double     snapDispMax,    //    =  .5,   // Maximal (scaled by plDispNorm) disparity difference to snap to plane at low strength
				clt_parameters.plSnapDispWeight, //  final double     snapDispWeight, //    =  .5,   // Maximal disparity diff. by weight product to snap to plane
				clt_parameters.plSnapZeroMode,   // final int        snapZeroMode,   // Zero strength snap mode: 0: no special treatment, 1 - strongest, 2 - farthest
				1); // final int        debugLevel)
		
		//

		DisparityProcessor dp = new DisparityProcessor(this, clt_parameters.transform_size * geometryCorrection.getScaleDzDx());
		boolean [] grown = these_tiles.clone();
		growTiles(
				2,          // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				grown,      // boolean [] tiles,
				null);     // boolean [] prohibit)
		boolean [] border = grown.clone();
		for (int i = 0; i < border.length; i++) border[i] &= !these_tiles[i];
		

		int num_unmatched = st.getNumNotMatched(
				snap_surf,
				these_tiles);
		int num_surf = st.getNumSurf(
				snap_surf,
				null); // these_tiles);
		if (debugLevel > -1){
			System.out.println("thirdPassSetup(): surfaces = "+num_surf+" left unmatched="+num_unmatched);
		}
		int [][] neib_surf =   new int [num_surf][];
		double [][] dbg_neibs = new double [num_surf][];
		for (int ns = 0 ; ns < num_surf; ns ++){
			
			neib_surf[ns] = st.getNeighbors( // creates neighbors mask from bitmask
					these_tiles, // final boolean [] selected, // or null
					snap_surf,   // final int []     snap_surf, // use this size - it matches image, not supertiles
					ns);        // final int        surf_index)
			
			dbg_neibs[ns]= dp.dbgShowNeighbors(
					these_tiles, // grown, // these_tiles, 
					neib_surf[ns], // _orig, // int [] neighbors,
					clt_parameters.transform_size, // int    tile_size,
					-1.0, // Double.NaN, // double bgnd,
					10.0); // double fgnd)
		}
		
		sdfa_instance.showArrays(dbg_neibs,
				                 tilesX*clt_parameters.transform_size,
				                 tilesY*clt_parameters.transform_size,
				                 true,
				                 "surf-neighbors");		
		int numScans0 = 0;
		for (int ns = 0 ; ns < num_surf; ns ++){
//		for (int ns = 1 ; (ns < num_surf) && (ns < 2); ns ++){
			boolean [] surface_selection = new boolean [these_tiles.length];
			for (int i = 0; i < surface_selection.length; i++) {
				surface_selection[i] = these_tiles[i] && (neib_surf[ns][i] >= 0);
			}
			
			
			
			grown = surface_selection.clone();
			growTiles(
					2,          // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
					grown,      // boolean [] tiles,
					null);     // boolean [] prohibit)
			border = grown.clone();
			for (int i = 0; i < border.length; i++) border[i] &= !these_tiles[i];
			
			int [] enum_clusters = enumerateClusters(
					true, // boolean diag_en,
					surface_selection); // grown); // these_tiles); // boolean [] tiles_src)
			double [] new_disparity0 = this_disparity.clone();

			numScans0 +=  createTileTasks(
					50, // int       maxClusters,
					0,  // int       minClusterArea,
					new_disparity0, // this_disparity, //  [] disparity_in, masked ok too
					enum_clusters, // int []    clusters_in,
					(ns == 1)?2:0); // debugLevel);
			
			if (debugLevel > -1){
				System.out.println("thirdPassSetup(): surface=" + ns + " created "+ numScans0+ " FPGA passes.");
			}
			
		}		
		
	}
	
	
	
	
	public void thirdPassSetupOld( // prepare tile tasks for the second pass based on the previous one(s)
			EyesisCorrectionParameters.CLTParameters           clt_parameters,
			double            disparity_far,    //
			double            disparity_near,   // 
			GeometryCorrection geometryCorrection,
			final int         threadsMax,  // maximal number of threads to launch                         
			final boolean     updateStatus,
			final int         debugLevel)
	{
		CLTPass3d scan_prev = clt_3d_passes.get(clt_3d_passes.size() -1); // get last one
		double [] dbg_orig_disparity =  null;
		double [] dbg_with_super_disp = null;
		double [] dbg_outlayers =       null;
		showDoubleFloatArrays sdfa_instance = null;
		if (debugLevel > -1) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		boolean [] these_tiles = scan_prev.getSelected();
		double [] this_disparity     = scan_prev.getDisparity(); // returns a copy of the FPGA-generated disparity combined with the target one
		double [] this_strength =      scan_prev.getStrength(); // cloned, can be modified/ read back
		
		
		dbg_orig_disparity = this_disparity.clone();
		
		SuperTiles st = scan_prev.getSuperTiles();
		if (clt_parameters.plSnapDispAny >= 0.0) {
			double [] snap_disparity =  st.snapDisparity(
					this_disparity,                  // final double []  disparity,
					this_strength,                   // final double []  strength,
					null, // this_selection,         // final boolean [] selection,       // can be null
					clt_parameters.plDispNorm,       //  final double     dispNorm, // plDispNorm
					clt_parameters.plSnapDispAny,    //final double     snapDispAny,    //    =  .2,   // Maximal (scaled by plDispNorm) disparity difference to snap to plane at any strength
					clt_parameters.plSnapStrengthAny,//final double     snapStrengthAny,    //    =  .2,   //  Maximal strength to fit any distance (if does not fit otherwise - treat as zero strength
					clt_parameters.plSnapNegAny,     //final double     snapNegAny,     //    =  .2,   // Maximal negative disparity difference from the best match
					clt_parameters.plSnapDispMax,    //  final double     snapDispMax,    //    =  .5,   // Maximal (scaled by plDispNorm) disparity difference to snap to plane at low strength
					clt_parameters.plSnapDispWeight, //  final double     snapDispWeight, //    =  .5,   // Maximal disparity diff. by weight product to snap to plane
					clt_parameters.plSnapZeroMode,   // final int        snapZeroMode,   // Zero strength snap mode: 0: no special treatment, 1 - strongest, 2 - farthest
					1); // final int        debugLevel)
			System.arraycopy(snap_disparity,  0, this_disparity, 0, snap_disparity.length);
		}
		if (clt_parameters.replaceWeakOutlayers) {
			boolean[] outlayers = scan_prev.replaceWeakOutlayers(
					null, // final boolean [] selection,
					clt_parameters.outlayerStrength , //final double weakStrength,    // strength to be considered weak, subject to this replacement
					clt_parameters.outlayerDiff, // final double maxDiff)
					clt_parameters.outlayerDiffPos, // final double maxDiff)
					clt_parameters.outlayerDiffNeg, // final double maxDiff)
					0.5 * disparity_far,
					2.0 * disparity_near,
					debugLevel);
			dbg_outlayers = new double[outlayers.length];

			for (int i = 0; i < outlayers.length; i++){
				dbg_outlayers[i] = outlayers[i]? 1.0:0.0;
			}			

			double [] masked_filtered = scan_prev.getDisparity().clone();
			if (clt_parameters.stShow){
				String [] dbg_disp_tiltes={"masked", "filtered", "disp_combo", "disparity","st_disparity", "strength", "st_strength","outlayers"};
				double [][] dbg_disp = new double [dbg_disp_tiltes.length][];
				dbg_disp[0] = masked_filtered;            //+
				dbg_disp[1] = scan_prev.getDisparity();   //+
				dbg_disp[2] = dbg_with_super_disp;        // -
				dbg_disp[3] = dbg_orig_disparity;         //+ 
				dbg_disp[4] = scan_prev.getBgDisparity(); 
				dbg_disp[5] = scan_prev.getStrength();    //+
				dbg_disp[6] = scan_prev.getBgStrength();
				dbg_disp[7] = dbg_outlayers;              // + (all 0)
				sdfa_instance.showArrays(dbg_disp, tilesX, tilesY, true, "disparity_supertiles",dbg_disp_tiltes);
			}
		}

		//clt_parameters.transform_size;		
		DisparityProcessor dp = new DisparityProcessor(this, clt_parameters.transform_size * geometryCorrection.getScaleDzDx());

		boolean [] grown = these_tiles.clone();
		growTiles(
				2,          // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				grown,      // boolean [] tiles,
				null);     // boolean [] prohibit)
		boolean [] border = grown.clone();
		for (int i = 0; i < border.length; i++) border[i] &= !these_tiles[i];
		//	double [] dbg_before = scan_prev.getDisparity().clone();
		scan_prev.fixNaNDisparity(
				grown, // border, // boolean [] select,   // which tiles to correct (null - all)
				scan_prev.getDisparity(), // double [] disparity,
				scan_prev.getStrength()); // double [] strength)
		
		//	double [] dbg_after = scan_prev.getDisparity().clone();
///		for (int i = 0; i < masked_filtered.length; i++){
///			if (!grown[i])   masked_filtered[i] = Double.NaN;
///		}

		int [] neighbors = dp.getNeighbors( // creates neighbors mask from bitmask
				grown, // these_tiles, // grown, // these_tiles, // boolean [] selected,
				tilesX);
		//	int [] neighbors_orig = neighbors.clone();
		double [] dbg_neib = dp.dbgShowNeighbors(
				grown, // these_tiles, // grown, // these_tiles,
				neighbors, // _orig, // int [] neighbors,
				clt_parameters.transform_size, // int    tile_size,
				-1.0, // double bgnd,
				1.0); // double fgnd)

		double [] new_disparity = this_disparity.clone();
		double [][]dbgDeriv = new double [2][]; // [these_tiles.length];
		//	sdfa_instance.showArrays(dbg_neib,tilesX*clt_parameters.transform_size, tilesY*clt_parameters.transform_size,"neighbors");

		dp.smoothDisparity(
				clt_parameters.tiDispPull,   // final double     dispPull, // clt_parameters.tiDispPull or 0.0
				3,                           // final int        mask,     // 1 - work on internal elements, 2 - on border elements, 3 - both (internal first);
				clt_parameters.tiIterations, //  final int        num_passes,
				Math.pow(10.0,  -clt_parameters.tiPrecision), // final double     maxDiff, // maximal change in any of the disparity values
				neighbors,                   // final int     [] neighbors, // +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
				new_disparity,               // final double  [] disparity,          // current disparity value
				this_disparity,              // final double  [] measured_disparity, // measured disparity
				this_strength,               // final double  [] strength,
				null, // this_hor_disparity,          // final double     hor_disparity, // not yet used
				null, // hor_strength_conv,           // final double     hor_strength, // not yet used
				these_tiles,                 // grown, // these_tiles,                 // final boolean [] selected,
				border, // final boolean [] border,
				clt_parameters,
				//			dbgDeriv,                    // final double [][] dbgDeriv, //double [2][len] or null;
				threadsMax,                  // maximal number of threads to launch                         
				debugLevel);

		double [] measured_disparity = 	dp.dbgRescaleToPixels(
				this_disparity,
				clt_parameters.transform_size); // int    tile_size)
		//break once
		double [][] stresses =           new double [clt_parameters.tiNumCycles + 1][];
		double [][] neibs_broken =       new double [clt_parameters.tiNumCycles][];
		double [][] smooth_disparities = new double [clt_parameters.tiNumCycles + 1][];
		smooth_disparities[0] = dp.dbgRescaleToPixels(
				new_disparity,
				clt_parameters.transform_size);
		boolean [] too_far =  new boolean [this_strength.length];
		boolean [] too_near = new boolean [this_strength.length];
		double  [] true_strength = this_strength.clone(); // this strength will be modified to remove too near tiles (too far - TBD)

/** */		
//		for (int numCycle = 0; numCycle <  clt_parameters.tiNumCycles; numCycle++) { //        = 5;     // Number of cycles break-smooth (after the first smooth)
		for (int numCycle = 0; numCycle <  1; numCycle++) { //        = 5;     // Number of cycles break-smooth (after the first smooth)

			// more smoothing
			double disp_pull = clt_parameters.tiDispPull;
			if (numCycle >= (clt_parameters.tiNumCycles -2)) {
				disp_pull = clt_parameters.tiDispPullPreFinal; // 0.1 * clt_parameters.tiDispPull;
				if (numCycle >= (clt_parameters.tiNumCycles -1)) {
					disp_pull = clt_parameters.tiDispPullPreFinal; // 0.01 * clt_parameters.tiDispPull;
				}
			}

			double breakScale = disp_pull / clt_parameters.tiDispPull; // 1.0;
/** */				
			
			if ((clt_parameters.tiBreakMode & 1) != 0) {
				dp.breakDisparity( // break using derivatives
						clt_parameters.tiBreak3 * breakScale, // clt_parameters.tiBreak,     // final double      break4, // clt_parameters.tiBreak/0 allow disconnecting from neighbors (fg/bg)
						0,      // mode, // 0: 3-rd derivative, 1: - third * first (compare positive threshold only), 2: second by abs(first) (compare positive threshold only)
						neighbors,                  // final int     []  neighbors, // UPDATED +1 - up (N), +2 - up-right - NE, ... +0x80 - NW 
						new_disparity,              // final double  []  disparity,          // current disparity value
						grown, // these_tiles,                // final boolean []  selected,
						true,                       // final boolean     extend_flat, // if the tile is on the hor/vert edge, assume same disparity on the other side
						clt_parameters.tiBreakSame, //  final double      k_same,
						clt_parameters.tiBreakTurn, // final double      k_turn,
						clt_parameters,
						dbgDeriv,                   // final double [][] dbgDeriv, //double [2][len] or null;
						threadsMax,      // maximal number of threads to launch                         
						debugLevel);
			}

			if ((clt_parameters.tiBreakMode & 2) != 0) {
				dp.breakDisparity( // break using derivatives
						clt_parameters.tiBreak31 * breakScale, // clt_parameters.tiBreak,     // final double      break4, // clt_parameters.tiBreak/0 allow disconnecting from neighbors (fg/bg)
						1,      // mode, // 0: 3-rd derivative, 1: - third * first (compare positive threshold only), 2: second by abs(first) (compare positive threshold only)
						neighbors,                  // final int     []  neighbors, // UPDATED +1 - up (N), +2 - up-right - NE, ... +0x80 - NW 
						new_disparity,              // final double  []  disparity,          // current disparity value
						grown, // these_tiles,                // final boolean []  selected,
						true,                       // final boolean     extend_flat, // if the tile is on the hor/vert edge, assume same disparity on the other side
						clt_parameters.tiBreakSame, //  final double      k_same,
						clt_parameters.tiBreakTurn, // final double      k_turn,
						clt_parameters,
						dbgDeriv,                   // final double [][] dbgDeriv, //double [2][len] or null;
						threadsMax,      // maximal number of threads to launch                         
						debugLevel);
			}

			if ((clt_parameters.tiBreakMode & 4) != 0) {
				dp.breakDisparity( // break using derivatives
						clt_parameters.tiBreak21 * breakScale, // clt_parameters.tiBreak,     // final double      break4, // clt_parameters.tiBreak/0 allow disconnecting from neighbors (fg/bg)
						2,      // mode, // 0: 3-rd derivative, 1: - third * first (compare positive threshold only), 2: second by abs(first) (compare positive threshold only)
						neighbors,                  // final int     []  neighbors, // UPDATED +1 - up (N), +2 - up-right - NE, ... +0x80 - NW 
						new_disparity,              // final double  []  disparity,          // current disparity value
						grown, // these_tiles,                // final boolean []  selected,
						true,                       // final boolean     extend_flat, // if the tile is on the hor/vert edge, assume same disparity on the other side
						clt_parameters.tiBreakSame, //  final double      k_same,
						clt_parameters.tiBreakTurn, // final double      k_turn,
						clt_parameters,
						dbgDeriv,                   // final double [][] dbgDeriv, //double [2][len] or null;
						threadsMax,      // maximal number of threads to launch                         
						debugLevel);
			}

			if ((clt_parameters.tiBreakMode & 8) != 0) { // remove too far
				boolean [] found = dp.findNearFar(
						true, // boolean findFar,
						clt_parameters.tiBreakFar, //  double  threshold, // select tiles with non-zero strength that are far/near
						new_disparity,          // current (approximated) disparity value
						this_disparity, // measured disparity
						this_strength, // final double  []  strength,           // masked by previously removed tiles as far/near
						grown);//these_tiles); // final boolean []  selected);
				if (found != null){
					for (int i=0;i < too_far.length; i++){
						too_far[i] |=found[i];
						if (too_far[i]) this_strength[i] =0.0;
					}
				}
			}


			if ((clt_parameters.tiBreakMode & 16) != 0) { // remove too near
				boolean [] found = dp.findNearFar(
						false, // boolean findFar,
						clt_parameters.tiBreakNear, //  double  threshold, // select tiles with non-zero strength that are far/near
						new_disparity,          // current (approximated) disparity value
						this_disparity, // measured disparity
						this_strength, // final double  []  strength,           // masked by previously removed tiles as far/near
						grown);// these_tiles); // final boolean []  selected);
				if (found != null){
					for (int i=0;i < too_near.length; i++){
						too_near[i] |=found[i];
						if (too_near[i]) this_strength[i] =0.0;
					}
				}
			}
 /* */
			if ((dbgDeriv != null) && (dbgDeriv[0] != null)) {
				stresses[numCycle] = dp.dbgShowStress(
						dbgDeriv, // double [][] stress,
						clt_parameters.transform_size); // int    tile_size)
			} else {
				stresses[numCycle] = null;
			}

			// Now heal some broken gaps back

			double healOrtho = 0.0;
			if (numCycle >= (clt_parameters.tiNumCycles -2)) {
				healOrtho = clt_parameters.tiHealPreLast;
				if (numCycle >= (clt_parameters.tiNumCycles -1)) {
					healOrtho =clt_parameters.tiHealLast;
				}
			}
			if (healOrtho > 0.0){
				dp.reconnectDisparity( // connect disconnected tiles if they have close approximated disparity
						healOrtho, // final double      maxDiffOrto,
						0, // healOrtho*1.5, // final double      maxDiffDiagonal,
						neighbors, // final int     []  neighbors, // +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
						new_disparity, // final double  []  disparity,          // current disparity value
						grown, // these_tiles, // final boolean []  selected,
						threadsMax,      // maximal number of threads to launch                         
						debugLevel);

			}
			if (numCycle >= (clt_parameters.tiNumCycles -2)) { // only heal during last 2 passes?
				if (clt_parameters.tiHealSame > 0){
					int numHealed = dp.healSame( // returns number of new ortho connections // ty = 325
							neighbors,
							clt_parameters.tiHealSame , // int maxlen,
							// just to fill in diagonals
							grown, // these_tiles, // final boolean []  selected,        
							threadsMax,                         
							debugLevel);
					if (debugLevel > -1){
						System.out.println("Healed "+numHealed+" connections");
					}
				}
			}			

			neibs_broken[numCycle] = dp.dbgShowNeighbors(
					grown, // these_tiles, // grown, // these_tiles,
					neighbors, // _orig, // int [] neighbors,
					clt_parameters.transform_size, // int    tile_size,
					-1.0, // double bgnd,
					1.0); // double fgnd)

			dp.smoothDisparity(
					disp_pull, // clt_parameters.tiDispPull,   // final double     dispPull, // clt_parameters.tiDispPull or 0.0
					3,                           // final int        mask,     // 1 - work on internal elements, 2 - on border elements, 3 - both (internal first);
					clt_parameters.tiIterations, //  final int        num_passes,
					Math.pow(10.0,  -clt_parameters.tiPrecision), // final double     maxDiff, // maximal change in any of the disparity values
					neighbors,                   // final int     [] neighbors, // +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
					new_disparity,               // final double  [] disparity,          // current disparity value
					this_disparity,              // final double  [] measured_disparity, // measured disparity
					this_strength,               // final double  [] strength,
					null, // this_hor_disparity,          // final double     hor_disparity, // not yet used
					null, // hor_strength_conv,           // final double     hor_strength, // not yet used
					grown, //  these_tiles,    //????             // grown, // these_tiles,                 // final boolean [] selected,
					border, // final boolean [] border,
					clt_parameters,
					threadsMax,                  // maximal number of threads to launch                         
					debugLevel);
			smooth_disparities[numCycle+1] = dp.dbgRescaleToPixels(
					new_disparity,
					clt_parameters.transform_size);
		}
		
/** */
		

		// Just calculate stress, do not actually break (after last smoothing)
		dp.breakDisparity(
				0.0,                        // final double      break4, // clt_parameters.tiBreak/0 allow disconnecting from neighbors (fg/bg)
				0, // clt_parameters.tiBreakMode, // mode, // 0: 3-rd derivative, 1: - third * first (compare positive threshold only), 2: second by abs(first) (compare positive threshold only)
				neighbors,                  // final int     []  neighbors, // UPDATED +1 - up (N), +2 - up-right - NE, ... +0x80 - NW 
				new_disparity,              // final double  []  disparity,          // current disparity value
				grown, // these_tiles,                // final boolean []  selected,
				true,                       // final boolean     extend_flat, // if the tile is on the hor/vert edge, assume same disparity on the other side
				clt_parameters.tiBreakSame, //  final double      k_same,
				clt_parameters.tiBreakTurn, // final double      k_turn,
				clt_parameters,
				dbgDeriv,                   // final double [][] dbgDeriv, //double [2][len] or null;
				threadsMax,      // maximal number of threads to launch                         
				debugLevel);
		stresses[clt_parameters.tiNumCycles] = dp.dbgShowStress(
				dbgDeriv, // double [][] stress,
				clt_parameters.transform_size); // int    tile_size)

		//	String [] titles0 = {"neib", "neib_broken", "stress", "stress1", "disp-measured", "disp-result","disp-result1"};
		double [] dbg_far_near=new double[too_far.length];
		for (int i = 0; i < dbg_far_near.length; i++){
			if (!these_tiles[i])  dbg_far_near[i]= Double.NaN;
			else if (too_far[i])  dbg_far_near[i]= -1.0;
			else if (too_near[i]) dbg_far_near[i]= 1.0;
			else                  dbg_far_near[i]= 0.0;
		}
		double [][] disp_diff = new double [3][this_disparity.length];
		for (int i = 0; i< this_disparity.length; i++) {
			if (this_strength[i] > 0.0) disp_diff[0][i] = this_disparity[i]-new_disparity[i];
			else  disp_diff[0][i] = Double.NaN;
			if (too_far[i])             disp_diff[1][i] = this_disparity[i]-new_disparity[i];
			else  disp_diff[1][i] = Double.NaN;
			if (too_near[i])            disp_diff[2][i] = this_disparity[i]-new_disparity[i];
			else  disp_diff[2][i] = Double.NaN;
		}

		if (clt_parameters.show_neighbors) {
			int numImages = 2 + 3*clt_parameters.tiNumCycles+2+3 + 3;
			String [] titles_all = new String[numImages];
			double [][] dbg_img = new double[numImages][];
			int indx = 0;
			titles_all[indx] = "neib";
			dbg_img[indx++]= dbg_neib;
			for (int i = 0; i < neibs_broken.length; i++){
				titles_all[indx] = "neib_"+i;
				dbg_img[indx++]= neibs_broken[i];
			}
			for (int i = 0; i < stresses.length; i++){
				titles_all[indx] = "stress_"+i;
				dbg_img[indx++]= stresses[i];
			}
			titles_all[indx] = "disp_meas";
			dbg_img[indx++]= measured_disparity;

			for (int i = 0; i < smooth_disparities.length; i++){
				titles_all[indx] = "disp_"+i;
				dbg_img[indx++]= smooth_disparities[i];
			}
			titles_all[indx] = "far-/near+";
			dbg_img[indx++] =  dp.dbgRescaleToPixels(
					dbg_far_near,
					clt_parameters.transform_size);
			//		double [][] dbg_img = {dbg_neib, dbg_neib_broken, stress, stress1, measured_disparity, smooth_disparity,smooth_disparity1};

			//		sdfa_instance.showArrays(dbg_neib,tilesX*clt_parameters.transform_size, tilesY*clt_parameters.transform_size,"neighbors");
			titles_all[indx] = "strength_orig";
			dbg_img[indx++] =  dp.dbgRescaleToPixels(
					true_strength,
					clt_parameters.transform_size);
			titles_all[indx] = "strength_mod";
			dbg_img[indx++] =   dp.dbgRescaleToPixels(
					this_strength,
					clt_parameters.transform_size);
			titles_all[indx] = "diff_this";
			dbg_img[indx++] =   dp.dbgRescaleToPixels(
					disp_diff[0],
					clt_parameters.transform_size);
			titles_all[indx] = "diff_far";
			dbg_img[indx++] =   dp.dbgRescaleToPixels(
					disp_diff[1],
					clt_parameters.transform_size);
			titles_all[indx] = "diff_near";
			dbg_img[indx++] =   dp.dbgRescaleToPixels(
					disp_diff[2],
					clt_parameters.transform_size);
			sdfa_instance.showArrays(dbg_img, tilesX*clt_parameters.transform_size, tilesY*clt_parameters.transform_size,
					true, "neighbors", titles_all);
		}
		//disp_diff
		//************************************************
		int [][] flaps = dp.createOverlapGeometry(
				neighbors,       // +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
				grown, // these_tiles,     // final boolean []  selected, // only inner?
				border, // final boolean []  border,
				threadsMax,      // maximal number of threads to launch                         
				debugLevel);
		String [] titleFlaps = {"neib","N","NE","E","SE","S","SW","W","NW"};
		if (clt_parameters.show_flaps_dirs){
			double [][] dbg_flaps =  dp.dbgShowOverlaps(
					//				boolean [] selected,
					flaps, // int [][] flaps,
					clt_parameters.transform_size, // int    tile_size,
					-1.0, // double bgnd,
					1.0); // double fgnd)
			double [] dbg_neibs = dp.dbgShowNeighbors(
					grown, // these_tiles, // grown, // these_tiles,
					neighbors, // _orig, // int [] neighbors,
					clt_parameters.transform_size, // int    tile_size,
					-1.0, // double bgnd,
					1.0); // double fgnd)
			double [][] dbg_flaps_all = {dbg_neibs,dbg_flaps[0],dbg_flaps[1],dbg_flaps[2],dbg_flaps[3],dbg_flaps[4],dbg_flaps[5],dbg_flaps[6],dbg_flaps[7]};
			sdfa_instance.showArrays(dbg_flaps_all, tilesX*clt_parameters.transform_size, tilesY*clt_parameters.transform_size,
					true, "flaps-dirs", titleFlaps);
		}
		int [][][] clustersNO= dp.extractNonOlerlap(
				true, // diag_en,
				neighbors, // +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
				grown, // these_tiles, // final boolean []  selected, // only inner?
				border,          // border should be diagonal!
				threadsMax,      // maximal number of threads to launch                         
				debugLevel);
		if (clt_parameters.show_first_clusters){

			int dbg_max_cluster_show = 50;
			if (dbg_max_cluster_show > clustersNO.length) dbg_max_cluster_show = clustersNO.length;

			double [][] dbg_clusters_show = new double [titleFlaps.length+dbg_max_cluster_show][neighbors.length];
			for (int i = 0; i < neighbors.length; i++){
				dbg_clusters_show[0][i] = neighbors[i];
				for (int j = 0; j < 8; j++){
					if (flaps[i] != null) dbg_clusters_show[1+j][i] = flaps[i][j];
				}
			}

			//		String [] titleFlaps = {"neib","N","NE","E","SE","S","SW","W","NW"};
			String [] titleClusters = new String [titleFlaps.length+dbg_max_cluster_show];
			//		int indxClust = 0;
			for (int i = 0; i < titleFlaps.length; i++){
				titleClusters[i] = titleFlaps[i];
			}

			for (int i = 0; i < dbg_max_cluster_show; i++){
				titleClusters[titleFlaps.length+i] = "C"+i;
			}
			for (int nClust = 0; nClust <dbg_max_cluster_show; nClust++){
				for (int i = 0; i < clustersNO[nClust][0].length; i++){
					dbg_clusters_show[nClust + 9][clustersNO[nClust][0][i]] += 32.0; // 1.0;
				}
				for (int i = 0; i < clustersNO[nClust][1].length; i++){
					dbg_clusters_show[nClust + 9][clustersNO[nClust][1][i]] += 64.0; // 2.0;
				}
				for (int i = 0; i < clustersNO[nClust][2].length; i++){
					dbg_clusters_show[nClust + 9][clustersNO[nClust][2][i]] += 128.0;// 4.0;
				}
			}
			sdfa_instance.showArrays(dbg_clusters_show, tilesX, tilesY, true, "first "+dbg_max_cluster_show+" clusters",titleClusters);
		}
		int numScans= 0;

		if (clt_parameters.shUseFlaps) {
			numScans =  createTileOverlapTasks(
					clt_parameters.max_clusters, // 50,                                           // int       maxClusters,
					clt_parameters.shMinArea,                     // int       minClusterArea,
					new_disparity,                                // [] disparity_in, masked ok too
					this_strength,  // double [] strength_in,
					Math.pow(10.0,  -clt_parameters.tiPrecision), // double      maxChange,   // adjust border disparity until change is below this.
					clt_parameters.shMinStrength,                 // double    minStrength,
					clustersNO,                                   // int []    clusters_in,
					disparity_far,
					disparity_near,
					clt_parameters.show_shells,
					debugLevel);

		} else {
			int [] enum_clusters = enumerateClusters(
					true, // boolean diag_en,
					grown); // these_tiles); // boolean [] tiles_src)

			numScans =  createTileTasks(
					50, // int       maxClusters,
					0,  // int       minClusterArea,
					new_disparity, // this_disparity, //  [] disparity_in, masked ok too
					this_strength,  // double [] strength_in,
					0.0,            // double    minStrength,
					enum_clusters, // int []    clusters_in,
					disparity_far,
					disparity_near,
					debugLevel);
		}

		if (debugLevel > -1){
			System.out.println("secondPassSetup(): created "+ numScans+ " FPGA passes.");
		}
		if (debugLevel > 10){ // null pointer
			String [] titles = new String [clt_3d_passes.size()];
			double [][] disparities = new double [titles.length][tilesX*tilesY];
			for (int i = 0; i < titles.length; i++) {
				titles[i] = i+"_scan";
				double [][] disparityTiles =  clt_3d_passes.get(i).disparity;
				for (int ty = 0; ty < tilesY; ty++) for (int tx = 0; tx < tilesX; tx++){
					disparities[i][ty*tilesX+tx] = disparityTiles[ty][tx]; // null pointer
				}
			}
			sdfa_instance.showArrays(disparities,  tilesX, tilesY, true, "disparities_scans",titles);
		}
		//	return scan_next;
	}
	
	
	
	
	
	
//==================	
	public void secondPassSetupOld( // prepare tile tasks for the second pass based on the previous one(s)
			//			  final double [][][]       image_data, // first index - number of image in a quad
			EyesisCorrectionParameters.CLTParameters           clt_parameters,
			boolean           use_supertiles,
			int               bg_scan_index,
			// disparity range - differences from 
			double            disparity_far,    //
			double            disparity_near,   // 
			double            this_sure,        // minimal strength to be considered definitely background
			double            this_maybe,       // maximal strength to ignore as non-background
			double            sure_smth,        // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
			int               disparity_index,  // index of disparity value in disparity_map == 2 (0,2 or 4)
			GeometryCorrection geometryCorrection,
			final int         threadsMax,  // maximal number of threads to launch                         
			final boolean     updateStatus,
			final int         debugLevel)
	{
		CLTPass3d scan_bg =   clt_3d_passes.get(bg_scan_index); // 
		CLTPass3d scan_prev = clt_3d_passes.get(clt_3d_passes.size() -1); // get last one
//		CLTPass3d scan_next =new CLTPass3d(this);
		showDoubleFloatArrays sdfa_instance = null;
		if (debugLevel > -1) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		//TODO: for next passes - combine all selected for previous passes (all passes with smaller disparity)
		
		
		boolean [] these_tiles = combineHorVertDisparity(
				scan_prev,                     // final CLTPass3d   scan,
				scan_bg.selected, // clt_3d_passes.get(0).selected, // final boolean [] bg_tiles,          // get from selected in clt_3d_passes.get(0); 
				disparity_far,    //
				disparity_near,   // 
				this_sure,        // minimal strength to be considered definitely background
				this_maybe,       // maximal strength to ignore as non-background
				sure_smth,        // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
				clt_parameters,
				debugLevel);
		
		scan_prev.combineHorVertStrength(true, false); // strength now max of original and horizontal. Use scale  instead of boolean?
		
		double [] this_disparity     = scan_prev.getDisparity(); // returns a copy of the FPGA-generated disparity combined with the target one
		double [] this_strength =      scan_prev.getStrength(); // cloned, can be modified/ read back 
		//************************************************
		// Show supertiles histograms

		//		if (clt_parameters.stShow){

		// try renovated supertiles. Do twice to show both original and blured histograms
		double [] dbg_orig_disparity =  null;
		double [] dbg_with_super_disp = null;
		double [] dbg_outlayers =       null;
		boolean [] grown = these_tiles.clone();

		if (use_supertiles) {	
			String [] dbg_st_titles = {"raw", "blurred"+clt_parameters.stSigma,"max-min-max"}; 
			double [][] dbg_hist = new double[dbg_st_titles.length][];

			scan_prev.setSuperTiles(
					clt_parameters.stStepNear,       // double     step_disparity,
					clt_parameters.stStepFar,        // double     step_near,
					clt_parameters.stStepThreshold,  // double     step_threshold,
					clt_parameters.stMinDisparity,   // double     min_disparity,
					clt_parameters.stMaxDisparity,   // double     max_disparity,
					clt_parameters.stFloor,          // double     strength_floor,
					clt_parameters.stPow,            // double     strength_pow,
					0.0,// NO BLUR double     stBlurSigma)
					false, //clt_parameters.stSmplMode,  // Use sample mode (false - regular tile mode)
					clt_parameters.stSmplSide,  // Sample size (side of a square)
					clt_parameters.stSmplNum,   // Number after removing worst
					clt_parameters.stSmplRms,   // Maximal RMS of the remaining tiles in a sample
					clt_parameters.stMeasSel); // bitmask of the selected measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			dbg_hist[0] = scan_prev.showDisparityHistogram();

			scan_prev.setSuperTiles(
					clt_parameters.stStepNear,       // double     step_disparity,
					clt_parameters.stStepFar,        // double     step_near,
					clt_parameters.stStepThreshold,  // double     step_threshold,
					clt_parameters.stMinDisparity,   // double     min_disparity,
					clt_parameters.stMaxDisparity,   // double     max_disparity,
					clt_parameters.stFloor,          // double     strength_floor,
					clt_parameters.stPow,            // double     strength_pow,
					clt_parameters.stSigma, // with blur double     stBlurSigma)
					false, //clt_parameters.stSmplMode,  // Use sample mode (false - regular tile mode)
					clt_parameters.stSmplSide,  // Sample size (side of a square)
					clt_parameters.stSmplNum,   // Number after removing worst
					clt_parameters.stSmplRms,   // Maximal RMS of the remaining tiles in a sample
					clt_parameters.stMeasSel); // bitmask of the selected measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			dbg_hist[1] = scan_prev.showDisparityHistogram();

			dbg_hist[2] = scan_prev.showMaxMinMax();

			if (clt_parameters.stShow){		
				int hist_width0 =  scan_prev.showDisparityHistogramWidth();
				int hist_height0 = dbg_hist[0].length/hist_width0;
				sdfa_instance.showArrays(dbg_hist, hist_width0, hist_height0, true, "disparity_supertiles_histograms",dbg_st_titles);
			}
			scan_prev.getBgDispStrength( // calculate (check non-null)?
					clt_parameters.stMinBgDisparity, // final double minBgDisparity, 
					clt_parameters.stMinBgFract); // final double minBgFract);

//			st_grown = these_tiles.clone();
			growTiles(
					2,          // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
					grown,      // boolean [] tiles,
					null);     // boolean [] prohibit)
			dbg_orig_disparity = scan_prev.getDisparity().clone();

			// combine weak with supertiles
			dbg_with_super_disp = scan_prev.combineSuper(
					true, //boolean updateStrength, // use ST strength if true, keep original (update disparity only) if false
					clt_parameters.stStrengthScale, // Multiply st strength if used instead of regular strength (only if updateStrength)
					clt_parameters.stUseDisp);
			
			if (dbg_with_super_disp != null) dbg_with_super_disp = dbg_with_super_disp.clone(); // else no super disparity available
		} else {
			growTiles(
					2,          // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
					grown,      // boolean [] tiles,
					null);     // boolean [] prohibit)

		}
		// replace weak outlaye tiles with weighted averages (modifies disparity)
		boolean[] outlayers = scan_prev.replaceWeakOutlayers(
				null, // final boolean [] selection,
				clt_parameters.outlayerStrength , //final double weakStrength,    // strength to be considered weak, subject to this replacement
				clt_parameters.outlayerDiff, // final double maxDiff)
				clt_parameters.outlayerDiffPos, // final double maxDiff)
				clt_parameters.outlayerDiffNeg, // final double maxDiff)
				0.5 * disparity_far,
				2.0 * disparity_near,
				debugLevel);

		dbg_outlayers = new double[outlayers.length];

		for (int i = 0; i < outlayers.length; i++){
			dbg_outlayers[i] = outlayers[i]? 1.0:0.0;
		}			
		
		double [] masked_filtered = scan_prev.getDisparity().clone();
		for (int i = 0; i < masked_filtered.length; i++){
			if (!grown[i])   masked_filtered[i] = Double.NaN;
		}
		if (clt_parameters.stShow){
			String [] dbg_disp_tiltes={"masked", "filtered", "disp_combo", "disparity","st_disparity", "strength", "st_strength","outlayers"};
			double [][] dbg_disp = new double [dbg_disp_tiltes.length][];
			dbg_disp[0] = masked_filtered;
			dbg_disp[1] = scan_prev.getDisparity();
			dbg_disp[2] = dbg_with_super_disp;
			dbg_disp[3] = dbg_orig_disparity;
			dbg_disp[4] = scan_prev.getBgDisparity();
			dbg_disp[5] = scan_prev.getStrength();
			dbg_disp[6] = scan_prev.getBgStrength();
			dbg_disp[7] = dbg_outlayers;
			sdfa_instance.showArrays(dbg_disp, tilesX, tilesY, true, "disparity_supertiles",dbg_disp_tiltes);
		}	
		// prepare new task and run
		double [][] disparityTask = new double [tilesY][tilesX];
		int [][]    tile_op =       new int [tilesY][tilesX];
		boolean [] borderTiles =    new boolean[tilesY*tilesX]; // to zero alpha in the images
		int op = ImageDtt.setImgMask(0, 0xf);
		op =     ImageDtt.setPairMask(op,0xf);
		op =     ImageDtt.setForcedDisparity(op,true);
		double [] prev_disparity = scan_prev.getDisparity();
		for (int ty = 0; ty < tilesY; ty++) for (int tx = 0; tx <tilesX; tx++){
			int indx =  tilesX * ty + tx;
			if (grown[indx]) {
				borderTiles[indx] =  !these_tiles[indx];
				disparityTask[ty][tx] = prev_disparity[indx];
				tile_op[ty][tx] = op;
			} else {
				disparityTask[ty][tx] = 0.0;
				tile_op[ty][tx] = 0;
				borderTiles[indx] = false;
			}
		}
		CLTPass3d scan_next =new CLTPass3d(this);
		scan_next.disparity =     disparityTask;
		scan_next.tile_op =       tile_op;
		scan_next.border_tiles =  borderTiles;
		clt_3d_passes.add(scan_next);


		//		}

		
		
//clt_parameters.transform_size;		
		DisparityProcessor dp = new DisparityProcessor(this, clt_parameters.transform_size * geometryCorrection.getScaleDzDx());
/*		
		boolean [] grown = these_tiles.clone();
		
		growTiles(
//				1,          // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				2,          // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				grown,      // boolean [] tiles,
				null);     // boolean [] prohibit)
*/
		boolean [] border = grown.clone();
		for (int i = 0; i < border.length; i++) border[i] &= !these_tiles[i];
		
		scan_prev.fixNaNDisparity(
				grown, // border, // boolean [] select,   // which tiles to correct (null - all)
				scan_prev.getDisparity(), // double [] disparity,
				scan_prev.getStrength()); // double [] strength)
		
		int [] neighbors = dp.getNeighbors( // creates neighbors mask from bitmask
				these_tiles, // grown, // these_tiles, // boolean [] selected,
				tilesX);
//		int [] neighbors_orig = neighbors.clone();
		double [] dbg_neib = dp.dbgShowNeighbors(
				these_tiles, // grown, // these_tiles,
				neighbors, // _orig, // int [] neighbors,
				clt_parameters.transform_size, // int    tile_size,
				-1.0, // double bgnd,
				1.0); // double fgnd)

		double [] new_disparity = this_disparity.clone();
		double [][]dbgDeriv = new double [2][]; // [these_tiles.length];
//		sdfa_instance.showArrays(dbg_neib,tilesX*clt_parameters.transform_size, tilesY*clt_parameters.transform_size,"neighbors");
		
		dp.smoothDisparity(
				clt_parameters.tiDispPull,   // final double     dispPull, // clt_parameters.tiDispPull or 0.0
				3,                           // final int        mask,     // 1 - work on internal elements, 2 - on border elements, 3 - both (internal first);
				clt_parameters.tiIterations, //  final int        num_passes,
				Math.pow(10.0,  -clt_parameters.tiPrecision), // final double     maxDiff, // maximal change in any of the disparity values
				neighbors,                   // final int     [] neighbors, // +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
				new_disparity,               // final double  [] disparity,          // current disparity value
				this_disparity,              // final double  [] measured_disparity, // measured disparity
				this_strength,               // final double  [] strength,
				null, // this_hor_disparity,          // final double     hor_disparity, // not yet used
				null, // hor_strength_conv,           // final double     hor_strength, // not yet used
				these_tiles,                 // grown, // these_tiles,                 // final boolean [] selected,
				border, // final boolean [] border,
				clt_parameters,
//				dbgDeriv,                    // final double [][] dbgDeriv, //double [2][len] or null;
				threadsMax,                  // maximal number of threads to launch                         
				debugLevel);
		
		double [] measured_disparity = 	dp.dbgRescaleToPixels(
				this_disparity,
				clt_parameters.transform_size); // int    tile_size)
// break once
		double [][] stresses =           new double [clt_parameters.tiNumCycles + 1][];
		double [][] neibs_broken =       new double [clt_parameters.tiNumCycles][];
		double [][] smooth_disparities = new double [clt_parameters.tiNumCycles + 1][];
		smooth_disparities[0] = dp.dbgRescaleToPixels(
				new_disparity,
				clt_parameters.transform_size);
		boolean [] too_far =  new boolean [this_strength.length];
		boolean [] too_near = new boolean [this_strength.length];
		double [] true_strength = this_strength.clone(); // this strength will be modified to remove too near tiles (too far - TBD)

		for (int numCycle = 0; numCycle <  clt_parameters.tiNumCycles; numCycle++) { //        = 5;     // Number of cycles break-smooth (after the first smooth)
			
			
			double breakScale = 1.0;
			if ((clt_parameters.tiBreakMode & 1) != 0) {
				dp.breakDisparity( // break using derivatives
						clt_parameters.tiBreak3 * breakScale, // clt_parameters.tiBreak,     // final double      break4, // clt_parameters.tiBreak/0 allow disconnecting from neighbors (fg/bg)
						0,      // mode, // 0: 3-rd derivative, 1: - third * first (compare positive threshold only), 2: second by abs(first) (compare positive threshold only)
						neighbors,                  // final int     []  neighbors, // UPDATED +1 - up (N), +2 - up-right - NE, ... +0x80 - NW 
						new_disparity,              // final double  []  disparity,          // current disparity value
						these_tiles,                // final boolean []  selected,
						true,                       // final boolean     extend_flat, // if the tile is on the hor/vert edge, assume same disparity on the other side
						clt_parameters.tiBreakSame, //  final double      k_same,
						clt_parameters.tiBreakTurn, // final double      k_turn,
						clt_parameters,
						dbgDeriv,                   // final double [][] dbgDeriv, //double [2][len] or null;
						threadsMax,      // maximal number of threads to launch                         
						debugLevel);
			}

			if ((clt_parameters.tiBreakMode & 2) != 0) {
				dp.breakDisparity( // break using derivatives
						clt_parameters.tiBreak31 * breakScale, // clt_parameters.tiBreak,     // final double      break4, // clt_parameters.tiBreak/0 allow disconnecting from neighbors (fg/bg)
						1,      // mode, // 0: 3-rd derivative, 1: - third * first (compare positive threshold only), 2: second by abs(first) (compare positive threshold only)
						neighbors,                  // final int     []  neighbors, // UPDATED +1 - up (N), +2 - up-right - NE, ... +0x80 - NW 
						new_disparity,              // final double  []  disparity,          // current disparity value
						these_tiles,                // final boolean []  selected,
						true,                       // final boolean     extend_flat, // if the tile is on the hor/vert edge, assume same disparity on the other side
						clt_parameters.tiBreakSame, //  final double      k_same,
						clt_parameters.tiBreakTurn, // final double      k_turn,
						clt_parameters,
						dbgDeriv,                   // final double [][] dbgDeriv, //double [2][len] or null;
						threadsMax,      // maximal number of threads to launch                         
						debugLevel);
			}

			if ((clt_parameters.tiBreakMode & 4) != 0) {
				dp.breakDisparity( // break using derivatives
						clt_parameters.tiBreak21 * breakScale, // clt_parameters.tiBreak,     // final double      break4, // clt_parameters.tiBreak/0 allow disconnecting from neighbors (fg/bg)
						2,      // mode, // 0: 3-rd derivative, 1: - third * first (compare positive threshold only), 2: second by abs(first) (compare positive threshold only)
						neighbors,                  // final int     []  neighbors, // UPDATED +1 - up (N), +2 - up-right - NE, ... +0x80 - NW 
						new_disparity,              // final double  []  disparity,          // current disparity value
						these_tiles,                // final boolean []  selected,
						true,                       // final boolean     extend_flat, // if the tile is on the hor/vert edge, assume same disparity on the other side
						clt_parameters.tiBreakSame, //  final double      k_same,
						clt_parameters.tiBreakTurn, // final double      k_turn,
						clt_parameters,
						dbgDeriv,                   // final double [][] dbgDeriv, //double [2][len] or null;
						threadsMax,      // maximal number of threads to launch                         
						debugLevel);
			}

			if ((clt_parameters.tiBreakMode & 8) != 0) { // remove too far
				boolean [] found = dp.findNearFar(
						true, // boolean findFar,
						clt_parameters.tiBreakFar, //  double  threshold, // select tiles with non-zero strength that are far/near
						new_disparity,          // current (approximated) disparity value
						this_disparity, // measured disparity
						this_strength, // final double  []  strength,           // masked by previously removed tiles as far/near
						these_tiles); // final boolean []  selected);
				if (found != null){
					for (int i=0;i < too_far.length; i++){
						too_far[i] |=found[i];
						if (too_far[i]) this_strength[i] =0.0;
					}
				}
			}
			
			
			if ((clt_parameters.tiBreakMode & 16) != 0) { // remove too near
				boolean [] found = dp.findNearFar(
						false, // boolean findFar,
						clt_parameters.tiBreakNear, //  double  threshold, // select tiles with non-zero strength that are far/near
						new_disparity,          // current (approximated) disparity value
						this_disparity, // measured disparity
						this_strength, // final double  []  strength,           // masked by previously removed tiles as far/near
						these_tiles); // final boolean []  selected);
				if (found != null){
					for (int i=0;i < too_near.length; i++){
						too_near[i] |=found[i];
						if (too_near[i]) this_strength[i] =0.0;
					}
				}
			}
			
			if ((dbgDeriv != null) && (dbgDeriv[0] != null)) {
				stresses[numCycle] = dp.dbgShowStress(
						dbgDeriv, // double [][] stress,
						clt_parameters.transform_size); // int    tile_size)
			} else {
				stresses[numCycle] = null;
			}

			// Now heal some broken gaps back
			
			double healOrtho = 0.0;
			if (numCycle >= (clt_parameters.tiNumCycles -2)) {
				healOrtho = clt_parameters.tiHealPreLast;
					if (numCycle >= (clt_parameters.tiNumCycles -1)) {
						healOrtho =clt_parameters.tiHealLast;
					}
			}
			if (healOrtho > 0.0){
				dp.reconnectDisparity( // connect disconnected tiles if they have close approximated disparity
						healOrtho, // final double      maxDiffOrto,
						0, // healOrtho*1.5, // final double      maxDiffDiagonal,
						neighbors, // final int     []  neighbors, // +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
						new_disparity, // final double  []  disparity,          // current disparity value
						these_tiles, // final boolean []  selected,
						threadsMax,      // maximal number of threads to launch                         
						debugLevel);

			}
			if (numCycle >= (clt_parameters.tiNumCycles -2)) { // only heal during last 2 passes?
				if (clt_parameters.tiHealSame > 0){
					int numHealed = dp. healSame( // returns number of new ortho connections
							neighbors,
							clt_parameters.tiHealSame , // int maxlen,
							// just to fill in diagonals
							these_tiles, // final boolean []  selected,        
							threadsMax,                         
							debugLevel);
					if (debugLevel > -1){
						System.out.println("Healed "+numHealed+" connections");
					}
				}
			}			

			neibs_broken[numCycle] = dp.dbgShowNeighbors(
					these_tiles, // grown, // these_tiles,
					neighbors, // _orig, // int [] neighbors,
					clt_parameters.transform_size, // int    tile_size,
					-1.0, // double bgnd,
					1.0); // double fgnd)

			// more smoothing
			double disp_pull = clt_parameters.tiDispPull;
			if (numCycle >= (clt_parameters.tiNumCycles -2)) {
				 disp_pull = 0.1 * clt_parameters.tiDispPull;
					if (numCycle >= (clt_parameters.tiNumCycles -1)) {
						 disp_pull = 0.01 * clt_parameters.tiDispPull;
					}
			}
			dp.smoothDisparity(
					disp_pull, // clt_parameters.tiDispPull,   // final double     dispPull, // clt_parameters.tiDispPull or 0.0
					3,                           // final int        mask,     // 1 - work on internal elements, 2 - on border elements, 3 - both (internal first);
					clt_parameters.tiIterations, //  final int        num_passes,
					Math.pow(10.0,  -clt_parameters.tiPrecision), // final double     maxDiff, // maximal change in any of the disparity values
					neighbors,                   // final int     [] neighbors, // +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
					new_disparity,               // final double  [] disparity,          // current disparity value
					this_disparity,              // final double  [] measured_disparity, // measured disparity
					this_strength,               // final double  [] strength,
					null, // this_hor_disparity,          // final double     hor_disparity, // not yet used
					null, // hor_strength_conv,           // final double     hor_strength, // not yet used
					these_tiles,                 // grown, // these_tiles,                 // final boolean [] selected,
					border, // final boolean [] border,
					clt_parameters,
					threadsMax,                  // maximal number of threads to launch                         
					debugLevel);
			smooth_disparities[numCycle+1] = dp.dbgRescaleToPixels(
					new_disparity,
					clt_parameters.transform_size);
			

		}
		
		
		
		
		// Just calculate stress, do not actually break (after last smoothing)
		dp.breakDisparity(
				0.0,                        // final double      break4, // clt_parameters.tiBreak/0 allow disconnecting from neighbors (fg/bg)
				0, // clt_parameters.tiBreakMode, // mode, // 0: 3-rd derivative, 1: - third * first (compare positive threshold only), 2: second by abs(first) (compare positive threshold only)
				neighbors,                  // final int     []  neighbors, // UPDATED +1 - up (N), +2 - up-right - NE, ... +0x80 - NW 
				new_disparity,              // final double  []  disparity,          // current disparity value
				these_tiles,                // final boolean []  selected,
				true,                       // final boolean     extend_flat, // if the tile is on the hor/vert edge, assume same disparity on the other side
				clt_parameters.tiBreakSame, //  final double      k_same,
				clt_parameters.tiBreakTurn, // final double      k_turn,
				clt_parameters,
				dbgDeriv,                   // final double [][] dbgDeriv, //double [2][len] or null;
				threadsMax,      // maximal number of threads to launch                         
				debugLevel);
		stresses[clt_parameters.tiNumCycles] = dp.dbgShowStress(
				dbgDeriv, // double [][] stress,
				clt_parameters.transform_size); // int    tile_size)
		
//		String [] titles0 = {"neib", "neib_broken", "stress", "stress1", "disp-measured", "disp-result","disp-result1"};
		double [] dbg_far_near=new double[too_far.length];
		for (int i = 0; i < dbg_far_near.length; i++){
			if (!these_tiles[i])  dbg_far_near[i]= Double.NaN;
			else if (too_far[i])  dbg_far_near[i]= -1.0;
			else if (too_near[i]) dbg_far_near[i]= 1.0;
			else                  dbg_far_near[i]= 0.0;
		}
		double [][] disp_diff = new double [3][this_disparity.length];
		for (int i = 0; i< this_disparity.length; i++) {
			if (this_strength[i] > 0.0) disp_diff[0][i] = this_disparity[i]-new_disparity[i];
			else  disp_diff[0][i] = Double.NaN;
			if (too_far[i])             disp_diff[1][i] = this_disparity[i]-new_disparity[i];
			else  disp_diff[1][i] = Double.NaN;
			if (too_near[i])            disp_diff[2][i] = this_disparity[i]-new_disparity[i];
			else  disp_diff[2][i] = Double.NaN;
		}
		
		if (clt_parameters.show_neighbors) {
			int numImages = 2 + 3*clt_parameters.tiNumCycles+2+3 + 3;
			String [] titles_all = new String[numImages];
			double [][] dbg_img = new double[numImages][];
			int indx = 0;
			titles_all[indx] = "neib";
			dbg_img[indx++]= dbg_neib;
			for (int i = 0; i < neibs_broken.length; i++){
				titles_all[indx] = "neib_"+i;
				dbg_img[indx++]= neibs_broken[i];
			}
			for (int i = 0; i < stresses.length; i++){
				titles_all[indx] = "stress_"+i;
				dbg_img[indx++]= stresses[i];
			}
			titles_all[indx] = "disp_meas";
			dbg_img[indx++]= measured_disparity;

			for (int i = 0; i < smooth_disparities.length; i++){
				titles_all[indx] = "disp_"+i;
				dbg_img[indx++]= smooth_disparities[i];
			}
			titles_all[indx] = "far-/near+";
			dbg_img[indx++] =  dp.dbgRescaleToPixels(
					dbg_far_near,
					clt_parameters.transform_size);
			//		double [][] dbg_img = {dbg_neib, dbg_neib_broken, stress, stress1, measured_disparity, smooth_disparity,smooth_disparity1};

			//		sdfa_instance.showArrays(dbg_neib,tilesX*clt_parameters.transform_size, tilesY*clt_parameters.transform_size,"neighbors");
			titles_all[indx] = "strength_orig";
			dbg_img[indx++] =  dp.dbgRescaleToPixels(
					true_strength,
					clt_parameters.transform_size);
			titles_all[indx] = "strength_mod";
			dbg_img[indx++] =   dp.dbgRescaleToPixels(
					this_strength,
					clt_parameters.transform_size);
			titles_all[indx] = "diff_this";
			dbg_img[indx++] =   dp.dbgRescaleToPixels(
					disp_diff[0],
					clt_parameters.transform_size);
			titles_all[indx] = "diff_far";
			dbg_img[indx++] =   dp.dbgRescaleToPixels(
					disp_diff[1],
					clt_parameters.transform_size);
			titles_all[indx] = "diff_near";
			dbg_img[indx++] =   dp.dbgRescaleToPixels(
					disp_diff[2],
					clt_parameters.transform_size);
			sdfa_instance.showArrays(dbg_img, tilesX*clt_parameters.transform_size, tilesY*clt_parameters.transform_size,
					true, "neighbors", titles_all);
		}
//disp_diff
//************************************************
		int [][] flaps = dp.createOverlapGeometry(
				neighbors,       // +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
				these_tiles,     // final boolean []  selected, // only inner?
				border, // final boolean []  border,
				threadsMax,      // maximal number of threads to launch                         
				debugLevel);
		String [] titleFlaps = {"neib","N","NE","E","SE","S","SW","W","NW"};
		if (clt_parameters.show_flaps_dirs){
			double [][] dbg_flaps =  dp.dbgShowOverlaps(
					//				boolean [] selected,
					flaps, // int [][] flaps,
					clt_parameters.transform_size, // int    tile_size,
					-1.0, // double bgnd,
					1.0); // double fgnd)
			double [] dbg_neibs = dp.dbgShowNeighbors(
					these_tiles, // grown, // these_tiles,
					neighbors, // _orig, // int [] neighbors,
					clt_parameters.transform_size, // int    tile_size,
					-1.0, // double bgnd,
					1.0); // double fgnd)
			double [][] dbg_flaps_all = {dbg_neibs,dbg_flaps[0],dbg_flaps[1],dbg_flaps[2],dbg_flaps[3],dbg_flaps[4],dbg_flaps[5],dbg_flaps[6],dbg_flaps[7]};
			sdfa_instance.showArrays(dbg_flaps_all, tilesX*clt_parameters.transform_size, tilesY*clt_parameters.transform_size,
					true, "flaps-dirs", titleFlaps);
		}
		int [][][] clustersNO= dp.extractNonOlerlap(
				true, // diag_en,
				neighbors, // +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
				these_tiles, // final boolean []  selected, // only inner?
				border,          // border should be diagonal!
				threadsMax,      // maximal number of threads to launch                         
				debugLevel);
		if (clt_parameters.show_first_clusters){

			int dbg_max_cluster_show = 50;
			if (dbg_max_cluster_show > clustersNO.length) dbg_max_cluster_show = clustersNO.length;

			double [][] dbg_clusters_show = new double [titleFlaps.length+dbg_max_cluster_show][neighbors.length];
			for (int i = 0; i < neighbors.length; i++){
				dbg_clusters_show[0][i] = neighbors[i];
				for (int j = 0; j < 8; j++){
					if (flaps[i] != null) dbg_clusters_show[1+j][i] = flaps[i][j];
				}
			}

			//		String [] titleFlaps = {"neib","N","NE","E","SE","S","SW","W","NW"};
			String [] titleClusters = new String [titleFlaps.length+dbg_max_cluster_show];
			//		int indxClust = 0;
			for (int i = 0; i < titleFlaps.length; i++){
				titleClusters[i] = titleFlaps[i];
			}

			for (int i = 0; i < dbg_max_cluster_show; i++){
				titleClusters[titleFlaps.length+i] = "C"+i;
			}
			for (int nClust = 0; nClust <dbg_max_cluster_show; nClust++){
				for (int i = 0; i < clustersNO[nClust][0].length; i++){
					dbg_clusters_show[nClust + 9][clustersNO[nClust][0][i]] += 32.0; // 1.0;
				}
				for (int i = 0; i < clustersNO[nClust][1].length; i++){
					dbg_clusters_show[nClust + 9][clustersNO[nClust][1][i]] += 64.0; // 2.0;
				}
				for (int i = 0; i < clustersNO[nClust][2].length; i++){
					dbg_clusters_show[nClust + 9][clustersNO[nClust][2][i]] += 128.0;// 4.0;
				}
			}
			sdfa_instance.showArrays(dbg_clusters_show, tilesX, tilesY, true, "first "+dbg_max_cluster_show+" clusters",titleClusters);
		}
		int numScans= 0;
		
		if (clt_parameters.shUseFlaps) {
			numScans =  createTileOverlapTasks(
					clt_parameters.max_clusters, // 50,                                           // int       maxClusters,
					clt_parameters.shMinArea,                     // int       minClusterArea,
					new_disparity,                                // [] disparity_in, masked ok too
					this_strength,  // double [] strength_in,
					Math.pow(10.0,  -clt_parameters.tiPrecision), // double      maxChange,   // adjust border disparity until change is below this.
					clt_parameters.shMinStrength,                 // double    minStrength,
					clustersNO,                                   // int []    clusters_in,
					disparity_far,
					disparity_near,
					clt_parameters.show_shells,
					debugLevel);
			
		} else {
			int [] enum_clusters = enumerateClusters(
					true, // boolean diag_en,
					these_tiles); // boolean [] tiles_src)

			numScans =  createTileTasks(
					50, // int       maxClusters,
					0,  // int       minClusterArea,
					new_disparity, // this_disparity, //  [] disparity_in, masked ok too
					this_strength,  // double [] strength_in,
					0.0,            // double    minStrength,
					enum_clusters, // int []    clusters_in,
					disparity_far,
					disparity_near,
					debugLevel);
		}
		
		if (debugLevel > -1){
			System.out.println("secondPassSetup(): created "+ numScans+ " FPGA passes.");
		}
		if (debugLevel > 0){
			String [] titles = new String [clt_3d_passes.size()];
			double [][] disparities = new double [titles.length][tilesX*tilesY];
			for (int i = 0; i < titles.length; i++) {
				titles[i] = i+"_scan";
				double [][] disparityTiles =  clt_3d_passes.get(i).disparity;
				for (int ty = 0; ty < tilesY; ty++) for (int tx = 0; tx < tilesX; tx++){
					disparities[i][ty*tilesX+tx] = disparityTiles[ty][tx];
				}
			}
			sdfa_instance.showArrays(disparities,  tilesX, tilesY, true, "disparities_scans",titles);
		}
//		return scan_next;
	}

	
	
	
	public int fixVerticalPoles( // return number of replaced cells
			CLTPass3d   scan,             // scan data to use
			boolean []  selection,        // start with only from selections (if not null, continue regardless)
			boolean []  tilesHor,         // horizontal correlation tiles used for composite disparity/strength;
			int         max_len,          // maximal length to cover
			double      poles_ratio,      // Maximal ratio of invisible to visible pole length
			double      min_new_strength, // set strength to hor_strength, but not less than this
			boolean     force_disparity,  // copy disparity down (false - use horDisparity
			boolean     keepStrength      // do not reduce composite strength from what it was before replacement 
			){
		double [] disparity =     scan.getDisparity(); 
		double [] hor_disparity = scan.getDisparity(2);
		double [] strength =      scan.getStrength();
		double [] hor_strength =  scan.getHorStrength();
		int tlen = tilesX * tilesY;
		int num_replaced = 0;
		for (int nTile = 0; nTile < tlen - tilesX; nTile++) {
			if (	tilesHor[nTile] && 
					!tilesHor[nTile + tilesX] &&
					(disparity[nTile] > disparity[nTile + tilesX]) && // NaN will fail - OK 
					((selection == null) || selection[nTile])) {
				// see how far to go
				int nTileEnd;
				for (nTileEnd = nTile + tilesX; nTileEnd < tlen; nTileEnd += tilesX){
					if (    tilesHor[nTileEnd] ||
							(disparity[nTileEnd] > disparity[nTile])){
						if (((nTileEnd - nTile) <= (max_len * tilesX)) || (max_len == 0)){
							// Calculate length of visible pole (above break)
							if (poles_ratio > 0.0){
								int pole_length = 1;
								for (int nt = nTile - tilesX; nt >=0 ; nt -= tilesX){
									if (!tilesHor[nt]) break;
									pole_length ++;
								}
								if ((nTileEnd - nTile) > (poles_ratio * pole_length * tilesX)){
									break; // too long invisible part
								}
							}

							for (int nt = nTile + tilesX; nt < nTileEnd; nt += tilesX){
								disparity[nt] = force_disparity?disparity[nTile]: hor_disparity[nt];
								if (!keepStrength || (strength[nt] < hor_strength[nt])) {
									strength[nt] = hor_strength[nt];
								}
								if (strength[nt] < min_new_strength) strength[nt] = min_new_strength;
								tilesHor[nt] = true;
								num_replaced ++;
							}
//							break;
						}
						break;
					}
				}
			}
		}
		
		return num_replaced;
	}
	
	
	public int bridgeFgndOrthoGap(
			EyesisCorrectionParameters.CLTParameters           clt_parameters,
			boolean vert, // verical pairs, horizontal features
			boolean disp_interpolate, // fill interpolated disparity for bridged over tiles, false - copy from  ortho_disparity (if not null)
			boolean closer_only, // only update disparity if laregr than was 
			boolean [] used_ortho,
			double [] this_strength,
			double [] ortho_strength_conv,
			double [] ortho_strength_raw,
			double [] this_disparity, // may be null - will not be updated
			double [] ortho_disparity // may be null
			){
		int numNew = 0;
		int nxt = vert? 1: tilesX; 
		int rm = tilesX - (vert? 1:0);
		int bm = tilesY - (vert? 0:1);
		int gap_2 = clt_parameters.ortho_bridge + 2;
		//		  boolean [] used_ortho_original = used_ortho.clone();
		if (clt_parameters.ortho_bridge > 0) {
			for (int ty = 0; ty < bm; ty++){
				for (int tx = 0; tx < rm; tx++){
					int indx = ty*tilesX + tx;
					if (used_ortho[indx] && !used_ortho[indx+nxt]){ // last before gap
						// look for end of gap
						int bridge_other = 0;
						for (int i = 1; i < gap_2; i++){
							if ( vert && (i+tx >= tilesX)) break;
							if (!vert && (i+ty >= tilesY)) break;
							int indx1 = indx + i*nxt;
							if (used_ortho[indx1]) {
								bridge_other = i;
								break;
							}
							// bridge over sufficient ortho strength or too uniform areas 
							if ((this_strength[indx1] > 0.0) && (ortho_strength_conv[indx1] < clt_parameters.ortho_sustain)) break;
						}
						if (bridge_other > 0){ // bridge the gap
							// verify bridge has sufficient support on at least one side
							boolean good_support = false;
							label_start_run:{
								for (int i = 1; i < clt_parameters.ortho_run; i++){
									if ( vert && (tx - i < 0))     break label_start_run;
									if (!vert && (ty - i < 0))     break label_start_run;
									if (!used_ortho[indx - i*nxt]) break label_start_run;
								}
								good_support = true;
							}
							if (!good_support){
								label_end_run:{
								for (int i = 1; i < clt_parameters.ortho_run; i++){
									if ( vert && (i+tx >= tilesX)) break label_end_run;
									if (!vert && (i+ty >= tilesY)) break label_end_run;
									if (!used_ortho[indx + i*nxt]) break label_end_run;
								}
								good_support = true;
							}
							}
							if (!good_support) continue; 

							double disp_start = 0, disp_end = 0;
							if ((ortho_disparity !=null) && disp_interpolate) { // interpolate ortho disparity, not yet the result one
								disp_start = ortho_disparity[indx];
								disp_end =   ortho_disparity[indx + bridge_other * nxt];
							}
							for (int i = 1; i < bridge_other; i++){
								int indx1 = indx + i*nxt;
								used_ortho[indx1] = true;
								numNew++;
								if ((ortho_disparity !=null) && disp_interpolate) { // interpolate ortho disparity, not yet the result one
									ortho_disparity[indx1]= disp_start + i* (disp_end - disp_start)/ bridge_other;
								}
							}
						}
					}
				}
			}
		}
		// now remove remaining short runs
		if (vert){
			for (int ty = 0; ty < tilesY; ty++){
				for (int tx = 0; tx < tilesX; tx++){
					for (; (tx < tilesX) && !used_ortho[ty*tilesX + tx]; tx++); // skip to first selected
					if (tx < tilesX) {
						double mx = ortho_strength_conv[ty*tilesX + tx];
						int tx1 = tx+1;
						for (; (tx1 < tilesX) && used_ortho[ty*tilesX + tx1];tx1++){
							if (ortho_strength_conv[ty*tilesX + tx1] > mx) mx = ortho_strength_conv[ty*tilesX + tx1];
						}
						if ((tx1 <= tx + clt_parameters.ortho_run) || (mx < clt_parameters.ortho_minmax)) {
							for (; tx <tx1; tx++){
								used_ortho[ty*tilesX + tx] = false; // remove short run
							}
						} else {
							tx = tx1;
						}
					}
				}
			}
		} else {
			for (int tx = 0; tx < tilesX; tx++){
				for (int ty = 0; ty < tilesY; ty++){
					for (; (ty < tilesY) && !used_ortho[ty*tilesX + tx]; ty++); // skip to first selected
					if (ty <tilesY) {
						double mx = ortho_strength_conv[ty*tilesX + tx];
						int ty1 = ty+1;
						for (; (ty1 < tilesY) && used_ortho[ty1*tilesX + tx]; ty1++){
							if (ortho_strength_conv[ty1*tilesX + tx] > mx) mx = ortho_strength_conv[ty1*tilesX + tx];
						}
						if ((ty1 <= ty + clt_parameters.ortho_run) || (mx < clt_parameters.ortho_minmax)) {
							for (; ty < ty1; ty++){
								used_ortho[ty*tilesX + tx] = false; // remove short run
							}
						} else {
							ty = ty1;
						}
					}
				}
			}

		}
		// update disparity from ortho when appropriate
		if ((this_disparity != null) && (ortho_disparity != null)) {
			for (int i = 0; i < ortho_disparity.length; i++) if (used_ortho[i] && (!closer_only || (ortho_disparity[i] > this_disparity[i]))){
				this_disparity[i] = ortho_disparity[i];
			}

			// average disparity along the runs, if they are close enough (remove outlayers?)
			if (clt_parameters.ortho_rms > 0){
				double ortho_rms2 = clt_parameters.ortho_rms * clt_parameters.ortho_rms; 

				if (vert){
					for (int ty = 0; ty < tilesY; ty++){
						for (int tx = 0; tx < tilesX; tx++){
							for (; (tx < tilesX) && !used_ortho[ty*tilesX + tx]; tx++); // skip to first selected
							if (tx < tilesX) {
								int tx1; // +1;
								double s0=0.0, s1=0.0, s2 = 0.0; 
								for (tx1 = tx; (tx1 < tilesX) && used_ortho[ty*tilesX + tx1];tx1++){
									int indx = ty*tilesX + tx1;
									double w = ortho_strength_raw[indx];
									if (w > 0) { // Double.NaN is not > 0.0
										s0 += w; 
										s1 += w * this_disparity[indx]; 
										s2 += w * this_disparity[indx] * this_disparity[indx];
									}
								}
								System.out.println ("vert: tx="+tx+" ty="+ty+" disp avg="+(s1/s0)+" rms = "+Math.sqrt((s2*s0 - s1*s1)/(s0*s0)));
								if ((s0 > 0) && ((s2*s0 - s1*s1)/(s0*s0) < ortho_rms2)) {
									s1 /= s0;
									for (tx1 = tx; (tx1 < tilesX) && used_ortho[ty*tilesX + tx1];tx1++){
										this_disparity[ty*tilesX + tx1] = s1;
									}									  
								}
								tx = tx1;
							}
						}
					}
				} else {
					for (int tx = 0; tx < tilesX; tx++){
						for (int ty = 0; ty < tilesY; ty++){
							for (; (ty < tilesY) && !used_ortho[ty*tilesX + tx]; ty++); // skip to first selected
							if (ty <tilesY) {
								int ty1;
								double s0=0.0, s1=0.0, s2 = 0.0; 
								for (ty1 = ty; (ty1 < tilesY) && used_ortho[ty1*tilesX + tx]; ty1++){
									int indx = ty1*tilesX + tx;
									double w = ortho_strength_raw[indx];
									if (w > 0) { // Double.NaN is not > 0.0
										s0 += w; 
										s1 += w * this_disparity[indx]; 
										s2 += w * this_disparity[indx] * this_disparity[indx];
									}
								}
								System.out.println ("hor: tx="+tx+" ty="+ty+" disp avg="+(s1/s0)+" rms = "+Math.sqrt((s2*s0 - s1*s1)/(s0*s0)));
								if ((s0 > 0) && ((s2*s0 - s1*s1)/(s0*s0) < ortho_rms2)) {
									s1 /= s0;
									for (ty1 = ty; (ty1 < tilesY) && used_ortho[ty1*tilesX + tx]; ty1++){
										this_disparity[ty1*tilesX + tx] = s1;
									}									  
								}
								ty = ty1;
							}
						}
					}
				}
			}
		}		  
		return numNew;  
	}

	// detect multi-tile vertical/horizontal features by convolving with averaging in one direction, sharpening in the other

	public double [] detectOrtho(
			double [] tiles,
			boolean vert,    // true for vertical correlation and _horizontal features
			int radius,      // one dimension - [-radius, +radius], other dimension -1,0, 1
			int debug)
	{
		double [][] kernel = null;
		if (vert) kernel = new double[3][2*radius+1];
		else      kernel = new double[2*radius+1][3];
		for (int i = 0; i < 2 * radius+1; i++){
			double w = Math.sin(Math.PI*(i+1)/2/(radius+1));
			for (int j = 0; j < 3; j++){
				double d = ((j == 1) ? 2.0: -0.5)* w;
				if (vert) kernel[j][i] = d;
				else      kernel[i][j] = d;
			}
		}
		if (debug > -1){
			System.out.println("detectOrtho("+vert+","+radius+")");
			for (int i = 0; i< kernel.length; i++){
				System.out.print(i+": ");
				for (int j = 0; j < kernel[i].length; j++){
					System.out.print(kernel[i][j]+" ");
				}
				System.out.println();
			}
		}
		return convolveTiles( tiles, kernel);
	}
	
	public void SharpBlurPair(
			double [] data,     // data array for in-place modification
			double [] strength, // data weights array for in-place modification
			double sigma,       // blur sigma
			double k,           // sharpen in orthogonal direction with (-k,2*k-1,-k). 0 - no sharpening
			double offset,      // subtract from strength, limit by 0.0
			boolean vert)       // true - sharpen vertically,  blur horizontally. False - sharpen horizontally, blur vertically
	{
		if (offset != 0.0) {
			for (int i = 0; i < data.length; i++) {
				strength[i] -= offset;
				if (strength[i] < 0.0) strength[i] = 0.0;
			}
		}
		
		for (int i = 0; i < data.length; i++) data[i] *= strength[i];
		SharpBlurTiles(strength, sigma, k, vert);

		// Maybe not needed? Can NaN be blured?
		for (int i = 0; i < data.length; i++) {
			if (Double.isNaN (data[i])) data[i] = 0.0;
		}
		
		
		SharpBlurTiles(data,     sigma, k, vert);
		for (int i = 0; i < data.length; i++) {
			if (strength[i] != 0.0) data[i] /= strength[i];
			else data[i] = Double.NaN;
		}
		if (offset != 0.0) {
			for (int i = 0; i < data.length; i++) {
				strength[i] += offset;
				if (strength[i] < 0.0) strength[i] = 0.0; // may be after sharpening
			}
		}
	}
	
	
	public void SharpBlurTiles(
			double [] tiles, // data array for in-place modification
			double sigma,    // blur sigma
			double k,        // sharpen in orthogonal direction with (-k,2*k-1,-k). 0 - no sharpening
			boolean vert)    // true - sharpen vertically,  blur horizontally. False - sharpen horizontally, blur vertically
	{
		if (k != 0){
			sharpTiles1d(
					tiles,
					k,     // 0.5 : -0.5/+2.0/-0.5
					vert);
		}
		blurTiles1d(
				tiles,
				sigma,
				!vert);
	}
	
	
	public void blurTiles1d(
			double [] tiles,
			double sigma,
			boolean vert)
	{
		(new DoubleGaussianBlur()).blur1Direction(
				tiles, // double [] pixels,
                tilesX, // int        width,
                tilesY, // int       height,
                sigma, // double     sigma,
                0.01,  // double   accuracy,
                !vert); // boolean xDirection
	}

	public void sharpTiles1d(
			double [] tiles,
			double k,     // 0.5 : -0.5/+2.0/-0.5
			boolean vert)
	{
		double [] src = tiles.clone();
		double k2 = 1.0 + 2 * k;
		if (vert){
			for (int tx = 0; tx < tilesX; tx++){
				for (int ty = 1; ty < (tilesY - 1); ty ++){
					int indx = ty*tilesX + tx;
					tiles[indx] = k2 * src[indx] - k * (src[indx - tilesX] + src[indx + tilesX]);
				}
			}
			
		} else {
			for (int ty = 0; ty < tilesY; ty++){
				for (int tx = 1; tx < (tilesX - 1); tx ++){
					int indx = ty*tilesX + tx;
					tiles[indx] = k2 * src[indx] - k * (src[indx - 1] + src[indx + 1]);
				}
			}
		}
	}

	
	
	
	
	
	public double [] convolveTiles(
			double [] tiles,
			double [][] kernel_in) // should be odd * odd, with non-zero sum (result will be normalized 
	{
		double [] rslt = tiles.clone(); // will just skip margins
		int half_height = (kernel_in.length - 1) / 2; 
		int half_width = (kernel_in[0].length - 1) / 2;

		System.out.println("convolveTiles(), half_height = "+half_height+", half_width = "+half_width);

		double [][] kernel = new double [2* half_height + 1] [2* half_width + 1];
		double s = 0;
		for (int i = 0; i < kernel.length; i++) for (int j = 0; j < kernel[0].length; j++) s+=kernel_in[i][j];
		s = 1.0/s;
		for (int i = 0; i < kernel.length; i++) for (int j = 0; j < kernel[0].length; j++) kernel[i][j] = s* kernel_in[i][j];

		for (int i = half_height; i < (tilesY - half_height); i++) for (int j = half_width; j < (tilesX - half_width); j++){
			double d = 0;
			for (int y = -half_height; y <= half_height; y++) for (int x = -half_width; x <= half_width; x++){
				d+= kernel[half_height-y][half_width-x] * tiles[tilesX *(i+y) + (j+x)];
			}
			rslt[tilesX *i + j] = d;
		}
		return rslt;
	}

	public int [][] getCoordIndices( // starting with 0, -1 - not selected
			Rectangle bounds,
			boolean [] selected
			)
	{
		int [][] indices = new int [bounds.height][bounds.width];
		int indx = 0;
		for (int y = 0; y < bounds.height; y++) {
			for (int x = 0; x < bounds.width; x++){
				if (selected[this.tilesX * (bounds.y + y) + (bounds.x + x)]){
					indices[y][x] = indx++;
				} else {
					indices[y][x] = -1;
				}
			}
		}
		return indices;
	}

	public double [][] getTexCoords( // get texture coordinates for indices
			int [][] indices)
	{
		int maxIndex = -1;
		int height = indices.length;
		int width = indices[0].length;
		outer_label:{
			for (int y = height - 1 ; y >= 0; y--) {
				for (int x = width - 1; x >= 0; x--){
					if (indices[y][x] >=0){
						maxIndex = indices[y][x]; 
						break outer_label;
					}
				}
			}
		}
		double [][] textureCoordinate = new double [maxIndex+1][2];
		int indx = 0;
		for (int y = 0;  indx <= maxIndex; y++) {
			for (int x = 0; (x < width) && (indx <= maxIndex); x++){
				if (indices[y][x] >=0){
					textureCoordinate[indx][0] = (x + 0.5)/width;
					textureCoordinate[indx][1] = (height - y - 0.5) / height; // y is up
					indx ++;
				}
			}
		}
		return textureCoordinate;
	}
	public double [] getIndexedDisparities( // get disparity for each index
			double [] disparity,
			double min_disparity,
			double max_disparity,
			Rectangle bounds,
			int [][]  indices,
			int       tile_size)
	{
		//		  int height = indices.length;
		int width = indices[0].length;
		int maxIndex = getMaxIndex(indices);
		double [] indexedDisparity = new double [maxIndex+1];
		int indx = 0;
		for (int y = 0;  indx <= maxIndex; y++) {
			for (int x = 0; (x < width) && (indx <= maxIndex); x++){
				if (indices[y][x] >=0){
					// center coordinates for 8*8 tile is [3.5,3.5]
					double disp = disparity[(bounds.y + y) * tilesX + (bounds.x + x)];
					if      (disp < min_disparity) disp = min_disparity;
					else if (disp > max_disparity) disp = max_disparity;
					indexedDisparity[indx] =disp;
					indx ++;
				}
			}
		}
		return indexedDisparity;
	}

	
	public double [][] getCoords( // get world XYZ in meters for indices
			double [] disparity,
			double min_disparity,
			double max_disparity,
			Rectangle bounds,
			int [][]  indices,
			int       tile_size,
			boolean   correctDistortions, // requires backdrop image to be corrected also
			GeometryCorrection geometryCorrection)
	{
		//		  int height = indices.length;
		int width = indices[0].length;
		int maxIndex = getMaxIndex(indices);
		double [][] coordinate = new double [maxIndex+1][];
		int indx = 0;
		for (int y = 0;  indx <= maxIndex; y++) {
			for (int x = 0; (x < width) && (indx <= maxIndex); x++){
				if (indices[y][x] >=0){
					// center coordinates for 8*8 tile is [3.5,3.5]
					double px = (bounds.x + x + 0.5) * tile_size - 0.5;
					double py = (bounds.y + y + 0.5) * tile_size - 0.5;
					double disp = disparity[(bounds.y + y) * tilesX + (bounds.x + x)];
					if      (disp < min_disparity) disp = min_disparity;
					else if (disp > max_disparity) disp = max_disparity;
					coordinate[indx] = geometryCorrection.getWorldCoordinates(
							px,
							py,
							disp,
							correctDistortions);
					indx ++;
				}
			}
		}
		return coordinate;
	}

	int getMaxIndex(int [][] indices)
	{
		int height = indices.length;
		int width = indices[0].length;
		for (int y = height - 1 ; y >= 0; y--) {
			for (int x = width - 1; x >= 0; x--){
				if (indices[y][x] >= 0){
					return indices[y][x]; 
				}
			}
		}
		return -1;
	}

	public int [][] filterTriangles(
			int  [][] triangles,
			double [] disparity, // disparities per vertex index
			double    maxDispDiff) // maximal disparity difference in a triangle
{
			final double min_avg = 0.5; // minimal average disparity to normalize triangle
		class Triangle {
			int [] points = new int [3];
			Triangle (int i1, int i2, int i3){
				points[0] = i1;
				points[1] = i2;
				points[2] = i3;
			}
		}
		ArrayList<Triangle> triList = new ArrayList<Triangle>();
		for (int i = 0; i < triangles.length; i++){
			double disp_avg = (triangles[i][0] + triangles[i][1]+ triangles[i][2])/3.0;
			if (disp_avg < min_avg) disp_avg = min_avg;
			loop:{
			for (int j = 0; j < 3; j++){
				int j1 = (j + 1) % 3;
				if (Math.abs(disparity[triangles[i][j]] - disparity[triangles[i][j1]]) > (disp_avg* maxDispDiff)) break loop;
			}
			triList.add(new Triangle(
					triangles[i][0],
					triangles[i][1],
					triangles[i][2]));
		}
		}
		int [][] filteredTriangles = new int [triList.size()][3];
		for (int i = 0; i < filteredTriangles.length; i++){
			filteredTriangles[i] = triList.get(i).points;
		}
		return filteredTriangles;
}
	

	public int [][] triangulateCluster(
			int [][]  indices)
	{
		int height = indices.length;
		int width = indices[0].length;
		//		  int [][] ind = new int [height][];
		//		  for (int i = 0; i < width; i++) ind[i] = indices[i].clone();
		class Triangle {
			int [] points = new int [3];
			Triangle (int i1, int i2, int i3){
				points[0] = i1;
				points[1] = i2;
				points[2] = i3;
			}
		}
		ArrayList<Triangle> triList = new ArrayList<Triangle>();
		for (int y = 0;  y < (height - 1); y++){
			for (int x = 0; x < width; x++){
				if (indices[y][x] >= 0){
					if ((x > 0) && (indices[y + 1][x - 1] >= 0) && (indices[y + 1][x] >= 0)){
						triList.add(new Triangle(
								indices[y][x],
								indices[y + 1][x],
								indices[y + 1][x - 1]));
					}
					if (x < (width - 1)) {
						if (indices[y + 1][x] >= 0){
							if (indices[y][x + 1] >= 0){
								triList.add(new Triangle(
										indices[y][x],
										indices[y][x + 1],
										indices[y + 1][x]));

							} else if (indices[y + 1][x + 1] >= 0){
								triList.add(new Triangle(
										indices[y][x],
										indices[y + 1][x + 1],
										indices[y + 1][x]));

							}
						} else if ((indices[y][x + 1] >= 0) && (indices[y + 1][x + 1] >= 0)) {
							triList.add(new Triangle(
									indices[y][x],
									indices[y][x + 1],
									indices[y + 1][x + 1]));
						}
					}
				}
			}
		}
		int [][] triangles = new int [triList.size()][3];
		for (int i = 0; i < triangles.length; i++){
			triangles[i] = triList.get(i).points;
		}
		return triangles;
	}

	int iSign (int a) {return (a > 0) ? 1 : ((a < 0)? -1 : 0);}
	public void testTriangles(
			String     texturePath,
			Rectangle  bounds,
			boolean [] selected,
			double []  disparity,
			int        tile_size,
			int [][]  indices,
			int [][]  triangles)
	{
		showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
		String [] titles = {"disparity","triangles"};
		double [][] dbg_img = new double [titles.length][tilesX*tilesY*tile_size*tile_size];
		for (int i = 0; i < selected.length; i++ ){
			double d = selected[i]?disparity[i]:Double.NaN;
			int y = i / tilesX;
			int x = i % tilesX;
			for (int dy = 0; dy <tile_size; dy ++){
				for (int dx = 0; dx <tile_size; dx ++){
					dbg_img[0][(y * tile_size + dy)*(tile_size*tilesX) + (x * tile_size + dx)] = d;
				}
			}
		}
		int maxIndex = getMaxIndex(indices);
		int [][] pxy = new int [maxIndex+1][2];
		int height = indices.length;
		int width = indices[0].length;

		for (int y = 0; y < height; y++){
			for (int x = 0; x < width; x++){
				if (indices[y][x] >= 0){
					pxy[indices[y][x]][0] =  (bounds.x + x)*tile_size + (tile_size/2);
					pxy[indices[y][x]][1] =  (bounds.y + y)*tile_size + (tile_size/2);
				}
			}
		}

		for (int i = 0; i < triangles.length; i++ ){
			for (int side = 0; side < triangles[i].length; side++){
				int [] pntIndx = {
						triangles[i][side],
						(side == (triangles[i].length -1)? triangles[i][0]:triangles[i][side+1])};
				int dx = iSign(pxy[pntIndx[1]][0] - pxy[pntIndx[0]][0]); 
				int dy = iSign(pxy[pntIndx[1]][1] - pxy[pntIndx[0]][1]);
				for (int j = 0; j < tile_size; j++){
					int x = pxy[pntIndx[0]][0] + dx*j;
					int y = pxy[pntIndx[0]][1] + dy*j;
					//					  int indx = y * tile_size * tilesX + x;
					//					  if (indx < dbg_img[1].length) {
					dbg_img[1][y * tile_size * tilesX + x] = 10.0; //1711748
					//					  } else {
					//						  indx += 0;
					//					  }
				}
			}
		}
		sdfa_instance.showArrays(dbg_img,
				tilesX * tile_size,
				tilesY * tile_size,
				true,
				"triangles-"+texturePath,
				titles);
	}
	/* Create a Thread[] array as large as the number of processors available.
	 * From Stephan Preibisch's Multithreading.java class. See:
	 * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD
	 */
	static Thread[] newThreadArray(int maxCPUs) {
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